import opengate_core as g4
import itk
import numpy as np
from .base import ActorBase
from ..utility import g4_units, ensure_filename_is_str
from ..exception import fatal, warning
from ..image import (
    create_3d_image,
    update_image_py_to_cpp,
    get_cpp_image,
    align_image_with_physical_volume,
    get_info_from_image,
    get_origin_wrt_images_g4_position,
)


class AMDMActor(g4.GateAMDMActor, ActorBase):
    """
    AMDMActor: compute a 3D edep/dose map for deposited
    energy/absorbed dose in the attached volume

    The dose map is parameterized with:
        - size (number of voxels)
        - spacing (voxel size)
        - translation (according to the coordinate system of the "attachedTo" volume)
        - no rotation

    Position:
    - by default: centered according to the "attachedTo" volume center
    - if the attachedTo volume is an Image AND the option "img_coord_system" is True:
        the origin of the attachedTo image is used for the output dose.
        Hence, the dose can be superimposed with the attachedTo volume

    Options
        - edep only for the moment
        - later: add dose, uncertainty, squared etc

    """

    type_name = "AMDMActor"

    def set_default_user_info(user_info):
        ActorBase.set_default_user_info(user_info)
        # required user info, default values
        mm = g4_units.mm
        user_info.size = [10, 10, 10]
        user_info.spacing = [1 * mm, 1 * mm, 1 * mm]
        user_info.output = "AMDM.mhd"  # FIXME change to 'output' ?
        user_info.translation = [0, 0, 0]
        user_info.img_coord_system = None
        user_info.output_origin = None
        # user_info.uncertainty = True
        # user_info.gray = False
        user_info.LUTfilename = "AMDM_LUT.txt"
        user_info.physical_volume_index = None
        user_info.hit_type = "random"

    def __init__(self, user_info):
        ActorBase.__init__(self, user_info)
        g4.GateAMDMActor.__init__(self, user_info.__dict__)
        # attached physical volume (at init)
        self.g4_phys_vol = None
        # default images (py side)
        self.py_restricted_edep_image = None
        # self.py_amdm_delta_image = None
        # self.py_amdm_omega_image = None
        # default uncertainty
        # self.uncertainty_image = None
        # internal states
        self.img_origin_during_run = None
        self.first_run = None
        self.output_origin = None
        self.LUTfilename = None
        self.AMDM_Bins = 6

    def __str__(self):
        u = self.user_info
        s = f'AMDMActor "{u.name}": dim={u.size} spacing={u.spacing} {u.output} tr={u.translation}'
        return s

    def __getstate__(self):
        # superclass getstate
        ActorBase.__getstate__(self)
        # do not pickle itk images
        self.py_restricted_edep_image = None
        # self.py_amdm_delta_image = None
        # self.py_amdm_omega_image = None
        return self.__dict__

    def initialize(self, volume_engine=None):
        """
        At the start of the run, the image is centered according to the coordinate system of
        the mother volume. This function computes the correct origin = center + translation.
        Note that there is a half-pixel shift to align according to the center of the pixel,
        like in ITK.
        """
        super().initialize(volume_engine)
        # create itk image (py side)
        size = np.array(self.user_info.size)
        spacing = np.array(self.user_info.spacing)

        self.py_restricted_edep_image = create_3d_image(size, spacing, "double")

        # self.py_amdm_delta_image = create_4d_image(
        #     size, spacing, FourthDimension=10
        # )
        # self.py_amdm_omega_image = create_4d_image(
        #     size, spacing, FourthDimension=10
        # )

        # compute the center, using translation and half pixel spacing
        self.img_origin_during_run = (
            -size * spacing / 2.0 + spacing / 2.0 + self.user_info.translation
        )
        # for initialization during the first run
        self.first_run = True

    def StartSimulationAction(self):
        # init the origin and direction according to the physical volume
        # (will be updated in the BeginOfRun)
        attached_to_volume = self.volume_engine.get_volume(self.user_info.mother)
        if self.user_info.physical_volume_index is None:
            physical_volume_index = 0
        else:
            physical_volume_index = self.user_info.physical_volume_index
        try:
            self.g4_phys_vol = attached_to_volume.g4_physical_volumes[
                physical_volume_index
            ]
        except IndexError:
            fatal(
                f"Error in the AMDMActor {self.user_info.name}. "
                f"Could not find the physical volume with index {physical_volume_index} "
                f"in volume '{self.user_info.mother}' to which this actor is attached. "
            )

        align_image_with_physical_volume(
            attached_to_volume,
            self.py_restricted_edep_image,
            initial_translation=self.user_info.translation,
        )
        # align_image_with_physical_volume(
        #     self.g4_phys_vol.GetName(),
        #     self.py_amdm_delta_image,
        #     self.user_info.translation,
        # )
        # align_image_with_physical_volume(
        #     self.g4_phys_vol.GetName(),
        #     self.py_amdm_omega_image,
        #     self.user_info.translation,
        # )

        # Set the real physical volume name
        self.fPhysicalVolumeName = str(self.g4_phys_vol.GetName())

        # print("AMDMActor: StartSimulationAction")
        # print("translation = ", self.user_info.translation)
        # print("origin = ", self.img_origin_during_run)

        # FIXME for multiple run and motion
        if not self.first_run:
            warning(f"Not implemented yet: AMDMActor with several runs")
        # send itk image to cpp side, copy data only the first run.
        update_image_py_to_cpp(
            self.py_restricted_edep_image,
            self.cpp_amdm_restricted_edep_image,
            self.first_run,
        )
        # print("AMDMActor: StartSimulationAction - creating image to cpp")
        # print(self.py_amdm_delta_image)
        # print(self.cpp_amdm_delta_image)
        # update_image_py_to_cpp(
        #     self.py_amdm_delta_image, self.cpp_amdm_delta_image, self.first_run
        # )
        # update_image_py_to_cpp(
        #     self.py_amdm_omega_image, self.cpp_amdm_gamma_image, self.first_run
        # )
        # print("AMDMActor: StartSimulationAction - after image to cpp")

        # now, indicate the next run will not be the first
        self.first_run = False

        # If attached to a voxelized volume, we may want to use its coord system.
        # So, we compute in advance what will be the final origin of the dose map
        attached_to_volume = self.simulation.volume_manager.volumes[
            self.user_info.mother
        ]
        vol_type = attached_to_volume.volume_type
        self.output_origin = self.img_origin_during_run

        # vol_name = self.user_info.mother
        # vol_type = self.simulation.get_volume_user_info(vol_name).type_name
        # vol_name = self.simulation.volume_manager.volumes[self.user_info.mother]
        # vol_type = attached_to_volume.volume_type

        # print(f"AMDMActor: StartSimulationAction - vol_type = {vol_type}")

        if vol_type == "ImageVolume":
            if self.user_info.img_coord_system:
                vol = self.volume_engine.g4_volumes[vol_name]
                # Translate the output dose map so that its center correspond to the image center.
                # The origin is thus the center of the first voxel.
                img_info = get_info_from_image(attached_to_volume.itk_image)
                dose_info = get_info_from_image(self.py_restricted_edep_image)
                self.output_origin = get_origin_wrt_images_g4_position(
                    img_info, dose_info, self.user_info.translation
                )
        else:
            if self.user_info.img_coord_system:
                warning(
                    f'AMDMActor "{self.user_info.name}" has '
                    f"the flag img_coord_system set to True, "
                    f"but it is not attached to an Image "
                    f'volume ("{attached_to_volume.name}", of type "{vol_type}"). '
                    f"So the flag is ignored."
                )
        # user can set the output origin
        if self.user_info.output_origin is not None:
            if self.user_info.img_coord_system:
                warning(
                    f'AMDMActor "{self.user_info.name}" has '
                    f"the flag img_coord_system set to True, "
                    f"but output_origin is set, so img_coord_system ignored."
                )
            self.output_origin = self.user_info.output_origin
        # print(
        #     f"AMDMActor: StartSimulationAction - output_origin = {self.output_origin}"
        # )

    def EndSimulationAction(self):
        # print("AMDMActor: EndSimulationAction")
        g4.GateAMDMActor.EndSimulationAction(self)

        # Get the itk image from the cpp side
        # Currently a copy. Maybe latter as_pyarray ?
        self.py_restricted_edep_image = get_cpp_image(
            self.cpp_amdm_restricted_edep_image
        )
        # self.py_amdm_delta_image = get_cpp_image(self.cpp_amdm_delta_image)
        # self.py_amdm_omega_image = get_cpp_image(self.cpp_amdm_gamma_image)

        # set the property of the output image:
        # in the coordinate system of the attached volume
        # FIXME no direction for the moment ?
        self.py_restricted_edep_image.SetOrigin(self.output_origin)
        # self.py_amdm_delta_image.SetOrigin(self.output_origin)
        # self.py_amdm_omega_image.SetOrigin(self.output_origin)

        # # omega image
        # n = check_filename_type(self.user_info.output).replace(
        #     ".mhd", "_omega.mhd"
        # )
        # itk.imwrite(self.py_amdm_omega_image, n)

        # # delta image
        # n = check_filename_type(self.user_info.output).replace(
        #     ".mhd", "_delta.mhd"
        # )
        # itk.imwrite(self.py_amdm_delta_image, n)

        # edep image
        n = ensure_filename_is_str(self.user_info.output).replace(
            ".mhd", "-restrictedEdep.mhd"
        )
        itk.imwrite(self.py_restricted_edep_image, n)
