import argparse
import sys
from pathlib import Path
from typing import Dict, Tuple
import SimpleITK as sitk

import numpy as np

# ---------------------------------------------------------------------
# MetaImage (.mhd/.raw) reading
# ---------------------------------------------------------------------
def load_and_flatten_metaimage(mhd_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    # Read the MHD file
    image = sitk.ReadImage(mhd_path)
    # print("Image size: ", image.GetSize())
    # print("Image spacing: ", image.GetSpacing())
    # print("Image origin: ", image.GetOrigin()) 
    # print("Image direction: ", image.GetDirection())
    # print("Image pixel type: ", image.GetPixelIDTypeAsString())
    # print("Image number of components per pixel: ", image.GetNumberOfComponentsPerPixel())
    # print("Image dimension: ", image.GetDimension())

    # Get the voxel values as a NumPy array
    image_array = sitk.GetArrayFromImage(image)
    # print("Original Image shape: ", image_array.shape)

    # print(image_array)
    # Flatten all but the last dimension
    # this leaves a two dim array with shape (Nvox, channels)
    image_array = image_array.reshape(-1, image_array.shape[-1])
    # print("Final Image shape with running numbers: ", image_array.shape)
    # print(image_array)

    return image_array

# ----------------- Feature computation ----------------- #

def fwhm_1d(x: np.ndarray, y: np.ndarray) -> float:
    """
    Compute approximate FWHM (full width at half maximum) of main peak.

    Returns 0.0 if it cannot be reasonably computed (e.g. flat or nonpositive).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    if x.size < 2 or y.size < 2:
        return 0.0

    max_idx = int(np.argmax(y))
    y_max = float(y[max_idx])
    if y_max <= 0.0:
        return 0.0

    half = 0.5 * y_max
    left = None
    right = None

    # left crossing
    for i in range(max_idx - 1, -1, -1):
        if (y[i] - half) * (y[i + 1] - half) <= 0:
            x1, x2 = x[i], x[i + 1]
            y1, y2 = y[i], y[i + 1]
            if y2 != y1:
                t = (half - y1) / (y2 - y1)
                left = x1 + t * (x2 - x1)
            else:
                left = x1
            break

    # right crossing
    for i in range(max_idx, y.size - 1):
        if (y[i] - half) * (y[i + 1] - half) <= 0:
            x1, x2 = x[i], x[i + 1]
            y1, y2 = y[i], y[i + 1]
            if y2 != y1:
                t = (half - y1) / (y2 - y1)
                right = x1 + t * (x2 - x1)
            else:
                right = x2
            break

    if left is None or right is None:
        return 0.0

    return float(max(right - left, 0.0))

# ---------------------------------------------------------------------
# Per-voxel macroscopic feature computation
# ---------------------------------------------------------------------

def compute_features_per_voxel(
    img: np.ndarray,
) -> Dict[str, np.ndarray]:
    """
    Compute spectral macroscopic features per voxel.

    Args:
        img       : N-D array, with spectral dimension at `spec_axis`.
        x_spec    : 1D array of spectral coordinates along that axis.
        spec_axis : axis index in `img` corresponding to x_spec.

    Returns:
        dict of feature_name -> N-D array (same shape as img without spectral axis).
    """

    Nvox = img.shape[0]
    base_shape = img.shape[:-1]
    y = img  # shape: (Nvox, num_channels)
    x_spec = np.arange(y.shape[1])  # spectral coordinates

    feats_flat: Dict[str, np.ndarray] = {}

    # Basic stats
    feats_flat["mean"] = np.mean(y, axis=1)
    feats_flat["median"] = np.median(y, axis=1)
    feats_flat["q10"] = np.percentile(y, 10.0, axis=1)
    feats_flat["q90"] = np.percentile(y, 90.0, axis=1)

    # Integral (area) using trapz
    feats_flat["integral"] = np.trapezoid(y, x_spec, axis=1)

    # Peak x & y
    peak_idx = np.argmax(y, axis=1)                 # shape: (Nvox,)
    feats_flat["peak_y"] = y[np.arange(Nvox), peak_idx]
    feats_flat["peak_x"] = x_spec[peak_idx]

    # Center of mass in x (spectral COM)
    y_sum = np.sum(y, axis=1)
    # avoid division by zero
    with np.errstate(divide="ignore", invalid="ignore"):
        com_x = np.sum(y * x_spec, axis=1) / np.where(y_sum != 0, y_sum, 1.0)
    feats_flat["com_x"] = com_x

    # FWHM per voxel (loop is fine here)
    fwhm = np.zeros(Nvox, dtype=np.float64)
    for i in range(Nvox):
        fwhm[i] = fwhm_1d(x_spec, y[i, :])
    feats_flat["fwhm"] = fwhm

    # Reshape all feature volumes back to base_shape
    feats: Dict[str, np.ndarray] = {}
    for name, arr_flat in feats_flat.items():
        feats[name] = arr_flat.reshape(base_shape)

    # print("Computed features: ", list(feats.keys()))
    # print(feats)
    return feats


# ---------------------------------------------------------------------
# Feature comparison (per voxel, per feature)
# ---------------------------------------------------------------------

def compare_feature_volumes(
    ref_feats: Dict[str, np.ndarray],
    test_feats: Dict[str, np.ndarray],
    rel_tols: Dict[str, float],
    abs_tols: Dict[str, float],
    min_pass_fraction: float,
) -> bool:
    """
    Compare per-voxel feature volumes between reference and test.

    Condition per voxel, per feature:
        |test - ref| <= abs_tol + rel_tol * |ref|

    For each feature:
        - compute pass_fraction = (#voxels passing) / (#voxels total)
        - print summary metrics

    Overall PASS if all features have pass_fraction >= min_pass_fraction.
    """
    feature_names = sorted(rel_tols.keys())
    overall_pass = True

    print("=== Per-voxel spectral comparison ===")
    # print(f"Using per-voxel condition: |test - ref| <= {abs_tols} + {rel_tols} * |ref|")
    # print(f"Required min_pass_fraction per feature: {min_pass_fraction}")
    # print("")

    for name in feature_names:

        abs_tol = abs_tols[name]
        rel_tol = rel_tols[name]
        
        ref = ref_feats[name].astype(np.float64)
        test = test_feats[name].astype(np.float64)

        if ref.shape != test.shape:
            overall_pass = False
            raise ValueError(
            f"Shape mismatch in feature '{name}': ref {ref.shape}, test {test.shape}"
            )

        diff = test - ref
        abs_diff = np.abs(diff)
        allowed = abs_tol + rel_tol * np.abs(ref)

        passed = abs_diff <= allowed
        total_voxels = ref.size
        passed_voxels = int(np.count_nonzero(passed))
        pass_fraction = passed_voxels / float(total_voxels)

        overall_pass = overall_pass and (pass_fraction >= min_pass_fraction)

        print(f"Feature: {name}: \t\t pass_fraction:\t{pass_fraction:.2f} "
              f"({passed_voxels}/{total_voxels})")

        # Display entries that failed
        if not passed.all():
            failed_voxels = np.where(~passed)[0]
            print(f"  Failed voxels: {failed_voxels.tolist()}")

    return overall_pass



# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------



def testMultiDimContent(ref_path,test_path,rel_tols,calcVoxelPercent=0.1):
    print("Comparing histogram MetaImage files:")
    print(f" Reference: {ref_path} \t Test: {test_path}")
    if not ref_path.exists():
        print(f"Reference file not found: {ref_path}")
        return False
    if not test_path.exists():
        print(f"Test file not found: {test_path}")
        return False

    ref_img=load_and_flatten_metaimage(ref_path)
    test_img=load_and_flatten_metaimage(test_path)
    if ref_img.shape != test_img.shape:
        raise ValueError(
            f"Image shape mismatch: ref {ref_img.shape}, test {test_img.shape}"
        )
    
    # Keep only first 15% of voxels
    num_voxels = ref_img.shape[0]
    num_keep = max(1, int(num_voxels * calcVoxelPercent))
    print(f"Comparing only first {num_keep} voxels out of {num_voxels} ({calcVoxelPercent*100}%)")
    ref_img = ref_img[:num_keep]
    test_img = test_img[:num_keep]
    
    # Compute per-voxel features for reference and test
    ref_feats = compute_features_per_voxel(ref_img)
    test_feats = compute_features_per_voxel(test_img)

    # Here is a set
    # rel_tols = {
    #     "integral": 0.20,   # 20% in area
    #     "com_x":    0.20,   # 20% in center of mass
    #     "fwhm":     0.20,   # 25% in FWHM
    #     "peak_x":   0.20,   # ~20% in peak position (in x-units)
    #     "q10":      0.20,   # 20% in lower tail level
    #     "q90":      0.20,   # 20% in upper tail level
    #     # "median":   0.20,
    #     "mean":     0.20,
    #     # "peak_y":   1.15,
    # }

    # Small absolute tolerances to catch tiny-reference cases
    abs_tols = {name: 1e-6 for name in ref_feats.keys()}

    ok = compare_feature_volumes(
        ref_feats=ref_feats,
        test_feats=test_feats,
        rel_tols=rel_tols,
        abs_tols=abs_tols,
        min_pass_fraction=0.99,
    )
    # print("Overall result:", "PASS" if ok else "FAIL")
    return ok


#####################################################################################
#!/usr/bin/env python3
"""
Compare two 1D MetaImage (.mhd/.raw) curves using macroscopic features.

Features compared by default:
  - integral      (area under curve)
  - com_x         (center of mass in x)
  - fwhm          (full width at half maximum of main peak)
  - peak_x        (x position of max)
  - q10, q90      (10th & 90th percentile of y)
  - median, mean  (optional; leave looser tolerance)

Result:
  - Prints a table of features, differences, and per-feature pass/fail
  - Exits with code 0 if ALL selected features are within tolerance, else 1
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


# ----------------- I/O: read MetaImage as 1D ----------------- #

def parse_mhd(path: Path) -> Dict[str, str]:
    info: Dict[str, str] = {}
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = [s.strip() for s in line.split("=", 1)]
                info[k] = v
    return info


def load_metaimage_1d(mhd_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load .mhd/.raw as a 1D curve y(x).

    Assumes 1D data stored as a 3D image with DimSize = (N, 1, 1)
    and uses first dimension as x.
    """
    hdr = parse_mhd(mhd_path)
    dimsize = list(map(int, hdr["DimSize"].split()))
    spacing = list(map(float, hdr["ElementSpacing"].split()))
    offset = list(map(float, hdr["Offset"].split()))
    data_file = hdr["ElementDataFile"]
    element_type = hdr["ElementType"]

    # Only a few types for simplicity; extend as needed
    type_map = {
        "MET_DOUBLE": "<f8",
        "MET_FLOAT": "<f4",
        "MET_SHORT": "<i2",
        "MET_USHORT": "<u2",
    }
    if element_type not in type_map:
        raise ValueError(f"Unsupported ElementType: {element_type}")
    dtype = type_map[element_type]

    n_voxels = dimsize[0] * dimsize[1] * dimsize[2]
    raw_path = mhd_path.with_name(data_file)

    data = np.fromfile(raw_path, dtype=dtype, count=n_voxels)
    if data.size != n_voxels:
        raise ValueError(
            f"Expected {n_voxels} values in {raw_path}, got {data.size}"
        )

    # Flatten and treat as 1D along first dimension
    y = data.reshape(dimsize[2], dimsize[1], dimsize[0]).ravel()

    # x from offset and spacing in first dimension
    x0 = offset[0]
    dx = spacing[0]
    x = x0 + dx * np.arange(dimsize[0], dtype=np.float64)

    return x, y


# ----------------- Feature computation ----------------- #

def safe_fwhm(x: np.ndarray, y: np.ndarray) -> float:
    """
    Approximate FWHM (full width at half maximum) of the main peak.
    Returns 0.0 if it can’t be computed reasonably.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if len(x) < 2:
        return 0.0

    max_idx = int(np.argmax(y))
    y_max = float(y[max_idx])
    if y_max <= 0:
        return 0.0

    half = 0.5 * y_max
    left = None
    right = None

    # Left crossing
    for i in range(max_idx - 1, -1, -1):
        if (y[i] - half) * (y[i + 1] - half) <= 0:
            x1, x2 = x[i], x[i + 1]
            y1, y2 = y[i], y[i + 1]
            if y2 != y1:
                t = (half - y1) / (y2 - y1)
                left = x1 + t * (x2 - x1)
            else:
                left = x1
            break

    # Right crossing
    for i in range(max_idx, len(x) - 1):
        if (y[i] - half) * (y[i + 1] - half) <= 0:
            x1, x2 = x[i], x[i + 1]
            y1, y2 = y[i], y[i + 1]
            if y2 != y1:
                t = (half - y1) / (y2 - y1)
                right = x1 + t * (x2 - x1)
            else:
                right = x2
            break

    if left is None or right is None:
        return 0.0

    return float(max(right - left, 0.0))


def compute_features(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    feats: Dict[str, float] = {}
    feats["mean"] = float(np.mean(y))
    feats["median"] = float(np.median(y))
    feats["q10"] = float(np.percentile(y, 10))
    feats["q90"] = float(np.percentile(y, 90))

    # Integral under curve
    feats["integral"] = float(np.trapezoid(y, x)) if len(x) > 1 else 0.0

    # Peak position and height
    peak_idx = int(np.argmax(y))
    feats["peak_x"] = float(x[peak_idx])
    feats["peak_y"] = float(y[peak_idx])

    # Center of mass in x
    y_sum = float(np.sum(y))
    feats["com_x"] = float(np.sum(x * y) / y_sum) if y_sum != 0.0 else float(np.mean(x))

    # FWHM
    feats["fwhm"] = safe_fwhm(x, y)

    return feats


# ----------------- Comparison logic ----------------- #

def compare_features(
    ref_feats: Dict[str, float],
    test_feats: Dict[str, float],
    rel_tols: Dict[str, float],
    abs_tols: Dict[str, float],
) -> bool:
    """
    Compare features with per-feature relative/absolute tolerances.

    Condition per feature:
        |test - ref| <= abs_tol + rel_tol * |ref|
    """
    print("=== Macroscopic feature comparison ===")
    print(f"{'Feature':<10} {'Ref':>12} {'Test':>12} {'Diff':>12} "
          f"{'RelErr':>10} {'Pass':>6}")

    overall_pass = True
    eps = 1e-12

    for name in sorted(rel_tols.keys()):
        r = ref_feats[name]
        t = test_feats.get(name, np.nan)
        diff = t - r
        denom = abs(r) if abs(r) > eps else 1.0
        rel_err = abs(diff) / denom

        rel_tol = rel_tols.get(name, 0.1)       # default 10%
        abs_tol = abs_tols.get(name, 1e-6)

        allowed = abs_tol + rel_tol * abs(r)
        passed = abs(diff) <= allowed

        print(
            f"{name:<10} {r:>12.6g} {t:>12.6g} "
            f"{diff:>12.6g} {rel_err:>10.3g} {str(passed):>6}"
        )

        overall_pass = overall_pass and passed

    print("\nOVERALL RESULT:", "PASS" if overall_pass else "FAIL")
    return overall_pass


# ----------------- CLI ----------------- #



def testOneDimensionalContent(ref_path,test_path,rel_tols):

    print("Comparing 1D MetaImage files:")
    print(f" Reference: {ref_path} \t Test: {test_path}")

    if not ref_path.exists():
        print(f"Reference file not found: {ref_path}")
        return False
    if not test_path.exists():
        print(f"Test file not found: {test_path}")
        return False
    x_ref, y_ref = load_metaimage_1d(ref_path)
    x_test, y_test = load_metaimage_1d(test_path)

    ref_feats = compute_features(x_ref, y_ref)
    test_feats = compute_features(x_test, y_test)

    # Here is a set 
    # rel_tols = {
    #     "integral": 0.20,   # 20% in area
    #     "com_x":    0.20,   # 20% in center of mass
    #     "fwhm":     0.20,   # 25% in FWHM
    #     "peak_x":   0.20,   # ~20% in peak position (in x-units)
    #     "q10":      0.20,   # 20% in lower tail level
    #     "q90":      0.20,   # 20% in upper tail level
    #     "median":   0.20,
    #     "mean":     0.20,
    #     # "peak_y":   1.15,
    # }

    # Small absolute tolerances to catch tiny-reference cases
    abs_tols = {name: 1e-6 for name in ref_feats.keys()}

    ok = compare_features(ref_feats, test_feats, rel_tols, abs_tols)
    return ok

###############################################################################################
