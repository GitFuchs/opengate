/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   ------------------------------------ -------------- */

#include "GateAMDMActor.h"
#include "G4Navigator.hh"
#include "G4RandomTools.hh"
#include "G4RunManager.hh"
#include "GateAMDMActor_lookUpTable.hh"
#include "GateHelpers.h"
#include "GateHelpersDict.h"
#include "GateHelpersImage.h"
#include "itkDivideImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkPasteImageFilter.h"
#include <string>

// Mutex that will be used by thread to write in the edep/dose image
G4Mutex SetPixelMutexAMDM = G4MUTEX_INITIALIZER;

GateAMDMActor::GateAMDMActor(py::dict &user_info)
    : GateVActor(user_info, true) {
  // Create the image pointer
  // (the size and allocation will be performed on the py side)
  // cpp_amdm_delta_image = VectorImageType::New();
  // cpp_amdm_gamma_image = VectorImageType::New();

  // define the number of bins to be used for the AMDM model
  fAMDM_total_bin_number = DictGetInt(user_info, "AMDM_Bins");
  // std::cout << "fAMDM_total_bin_number: " << fAMDM_total_bin_number
  //           << std::endl;

  fstoreMergingData = DictGetBool(user_info, "AMDM_Bins");
  // std::cout << "fstoreMergingData: " << fstoreMergingData << std::endl;

  cpp_amdm_delta_image = ImageType4D::New();
  cpp_amdm_gamma_image = ImageType4D::New();
  cpp_amdm_restricted_edep_image = ImageType::New();
  // initialize the image
  ImageType4D::IndexType start;
  start[0] = 0; // first index on X
  start[1] = 0; // first index on Y
  start[2] = 0; // first index on Z
  start[3] = 0; // first index on T

  fImageSpacing = DictGetG4ThreeVector(user_info, "spacing");
  fImageSize = DictGetG4ThreeVector(user_info, "size");

  ImageType4D::SizeType size;
  size[0] = fImageSize[0];          // size along X
  size[1] = fImageSize[1];          // size along Y
  size[2] = fImageSize[2];          // size along Z
  size[3] = fAMDM_total_bin_number; // size along T

  ImageType4D::SpacingType spacing;
  spacing[0] = fImageSpacing[0]; // spacing along X
  spacing[1] = fImageSpacing[1]; // spacing along Y
  spacing[2] = fImageSpacing[2]; // spacing along Z
  spacing[3] = 1;                // spacing along T
  // std::cout
  //     << "AMDM image size: " << size[0] << " " << size[1] << " " << size[2]
  //     << " " << size[3] << std::endl;

  ImageType4D::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  cpp_amdm_delta_image->SetRegions(region);
  cpp_amdm_delta_image->SetSpacing(spacing);
  cpp_amdm_delta_image->Allocate();
  cpp_amdm_delta_image->FillBuffer(0.);

  cpp_amdm_gamma_image->SetRegions(region);
  cpp_amdm_gamma_image->SetSpacing(spacing);
  cpp_amdm_gamma_image->Allocate();
  cpp_amdm_gamma_image->FillBuffer(0.);

  // // Action for this actor: during stepping
  fActions.insert("SteppingAction");
  fActions.insert("BeginOfRunAction");
  fActions.insert("EndOfRunAction");
  fActions.insert("EndSimulationAction");
  // // Option: compute uncertainty
  // fUncertaintyFlag = DictGetBool(user_info, "uncertainty");
  // Option: compute dose in Gray
  // fGrayFlag = DictGetBool(user_info, "gray");
  // translation
  fInitialTranslation = DictGetG4ThreeVector(user_info, "translation");
  // Hit type (random, pre, post etc)
  fHitType = DictGetStr(user_info, "hit_type");
  sOutputFileNameBase = DictGetStr(user_info, "output");
  sOutputFileNameBase = GateAMDMActor::remove_extension(sOutputFileNameBase);

  // // intialize and read in look up tables
  fLUTfilename = DictGetStr(user_info, "LUTfilename");
  if (fLUTfilename.empty()) {
    std::cout << "The file name for the AMDM LUT is empty. Aborting."
              << std::endl;
    exit(-1);
  }

  famdm_tables = new lookUpTable();

  famdm_tables->readLookUpTable(fLUTfilename);
  // famdm_tables->printLookUpTable();

  std::cout << "AMDM actor created" << std::endl;
}

void GateAMDMActor::ActorInitialize(py::dict &user_info) {
  // // intialize and read in look up tables
  fLUTfilename = DictGetStr(user_info, "LUTfilename");
  if (fLUTfilename.empty()) {
    std::cout << "The file name for the AMDM LUT is empty. Aborting."
              << std::endl;
    exit(-1);
  }

  lookUpTable *famdm_tables = new lookUpTable();
  famdm_tables->readLookUpTable(fLUTfilename);
  famdm_tables->printLookUpTable();
  // std::cout << "AMDM actor initialized" << std::endl;
}

void GateAMDMActor::BeginOfRunAction(const G4Run *) {
  // Important ! The volume may have moved, so we re-attach each run
  AttachImageToVolume<ImageType4D>(cpp_amdm_delta_image, fPhysicalVolumeName,
                                   fInitialTranslation);
  AttachImageToVolume<ImageType4D>(cpp_amdm_gamma_image, fPhysicalVolumeName,
                                   fInitialTranslation);
  //   AttachImageToVolume<
  //   AttachImageToVolume<ImageType4D>(cpp_amdm_restricted_edep_image,
  //   fPhysicalVolumeName,
  // >(cpp_amdm_delta_image, fPhysicalVolumeName,
  //                                        fInitialTranslation);
  //   AttachImageToVolume<
  //   AttachImageToVolume<ImageType4D>(cpp_amdm_restricted_edep_image,
  //   fPhysicalVolumeName,
  // >(cpp_amdm_gamma_image, fPhysicalVolumeName,
  //                                        fInitialTranslation);
  AttachImageToVolume<ImageType>(cpp_amdm_restricted_edep_image,
                                 fPhysicalVolumeName, fInitialTranslation);

  // // compute volume of a dose voxel
  // auto sp = cpp_amdm_image->GetSpacing();
  // fVoxelVolume = sp[0] * sp[1] * sp[2];
  // std::cout << "AMDM actor starting run"
  //           << std::endl;
}

void GateAMDMActor::SteppingAction(G4Step *step) {
  // std::cout << "a step"
  //           << std::endl;
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();

  // FIXME If the volume has multiple copy, touchable->GetCopyNumber(0) ?

  // consider random position between pre and post
  auto position = postGlobal;
  if (fHitType == "pre") {
    position = preGlobal;
  } else if (fHitType == "random") {
    auto x = G4UniformRand();
    auto direction = postGlobal - preGlobal;
    position = preGlobal + x * direction;
  } else if (fHitType == "middle") {
    auto direction = postGlobal - preGlobal;
    position = preGlobal + 0.5 * direction;
  } else if (fHitType == "post") {
    position = postGlobal;
  } else {
    std::cout << "Unknown hit type " << fHitType << " Aborting." << std::endl;
    exit(-1);
  }

  // get local position
  auto localPosition =
      touchable->GetHistory()->GetTransform(0).TransformPoint(position);

  // convert G4ThreeVector to itk PointType
  ImageType::PointType point;
  point[0] = localPosition[0];
  point[1] = localPosition[1];
  point[2] = localPosition[2];

  // get edep in MeV (take weight into account)
  auto w = step->GetTrack()->GetWeight();
  auto edep = step->GetTotalEnergyDeposit() / CLHEP::MeV * w;

  // get pixel index
  ImageType::IndexType index;
  // ImageType4D::IndexType index4D;
  bool isInside = cpp_amdm_restricted_edep_image->TransformPhysicalPointToIndex(
      point, index);

  // set value
  if (isInside) {
    // With mutex (thread)
    G4AutoLock mutex(&SetPixelMutexAMDM);

    G4int charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
    G4int atomic_mass = step->GetTrack()->GetDefinition()->GetBaryonNumber();

    // skip if charge is smaller or equal to zero, as outside of applicability
    // of MCF MKM model
    if (charge <= 0)
      return;

    // double cut = DBL_MAX;
    // important add units !
    // G4double density = step->GetTrack()->GetMaterial()->GetDensity() /
    // (CLHEP::g / CLHEP::cm3);
    G4double kinEnergy = step->GetTrack()->GetKineticEnergy() / (CLHEP::MeV);
    // if nan, that is if atomic mass is zero
    if (isnan(kinEnergy / atomic_mass)) {
      // std::cout << "kinetic Energy: " << kinEnergy << " atomic mass: " <<
      // atomic_mass << std::endl;
      return;
    }
    // std::cout << "edep: " << edep << "kinetic Energy: " << kinEnergy << "
    // atomic mass: " << atomic_mass << std::endl;

    bool success = false;
    std::vector<double> entry;
    // famdm_tables->printLookUpTable();
    // std::cout << "starting lookup" << std::endl;
    success = famdm_tables->findEntriesByEnergy(charge, kinEnergy / atomic_mass,
                                                entry);
    // first fAMDM_total_bin_number entries for gamma, second
    // fAMDM_total_bin_number entries for delta std::cout << "entry size: " <<
    // entry.size() << std::endl; std::cout
    // << "entry 0: " << entry[0] << std::endl; std::cout << "success is not yet
    // there" << std::endl;
    if (success) {
      // std::cout << "success is true" << std::endl;
      // store restricted edep
      ImageAddValue<ImageType>(cpp_amdm_restricted_edep_image, index, edep);
      // std::cout << "added edep to edep image with index: " << index <<
      // std::endl;

      G4double gamma_d = 0;
      G4double delta_e = 0;
      // loop over all bins
      // first half of lookUpTable contains gamma, second half contains delta
      int secondHalf = fAMDM_total_bin_number;

      // using PixelType = itk::Vector<float, 10>;

      // ImageType::IndexType index;
      // VectorImageType::PixelType BinnedPixel_delta_e;
      // VectorImageType::PixelType BinnedPixel_gamma_d;

      // BinnedPixel_delta_e.SetSize(fAMDM_total_bin_number);
      // BinnedPixel_gamma_d.SetSize(fAMDM_total_bin_number);

      // // convert G4ThreeVector to itk PointType
      ImageType4D::PointType point4D;
      ImageType4D::IndexType index4D;

      // std::cout << "Dimension of cpp_amdm_delta_image " <<
      // cpp_amdm_delta_image->GetLargestPossibleRegion().GetSize() <<
      // std::endl; std::cout << "Dimension of cpp_amdm_gamma_image " <<
      // cpp_amdm_gamma_image->GetLargestPossibleRegion().GetSize() <<
      // std::endl; std::cout << "Origin of cpp_amdm_delta_image" <<
      // cpp_amdm_delta_image->GetOrigin() << std::endl; std::cout << "Origin of
      // cpp_amdm_gamma_image" << cpp_amdm_gamma_image->GetOrigin() <<
      // std::endl; std::cout << "Spacing of cpp_amdm_delta_image" <<
      // cpp_amdm_delta_image->GetSpacing() << std::endl; std::cout << "Spacing
      // of cpp_amdm_gamma_image" << cpp_amdm_gamma_image->GetSpacing() <<
      // std::endl;

      point4D[0] = localPosition[0];
      point4D[1] = localPosition[1];
      point4D[2] = localPosition[2];
      point4D[3] = 0; // bin index of AMDM

      index4D = cpp_amdm_delta_image->TransformPhysicalPointToIndex(point4D);
      // std::cout << "point4D: " << point4D << std::endl;
      // std::cout << "index4D: " << index4D << std::endl;

      // for every bin, we have two LUT entries going from 0 to 9=kMaxMCFBins-1
      // LUT entries going from 0 to 9=kMaxMCFBins-1
      // we loop through all bins and create the pixel value for each bin
      for (int i = 0; i < fAMDM_total_bin_number; i++) {
        index4D[3] = i;
        // std::cout << "point4D: " << point4D << std::endl;
        // std::cout << "index4D: " << index4D << std::endl;

        // index of edep image reuseable for delta and gamma image

        delta_e = entry[i + secondHalf] * edep;
        gamma_d = entry[i] * delta_e;
        // std::cout << "storing gamma_d and delta_e in image" << std::endl;

        // store gamma_d and delta_e
        ImageAddValue<ImageType4D>(cpp_amdm_delta_image, index4D, delta_e);
        ImageAddValue<ImageType4D>(cpp_amdm_gamma_image, index4D, gamma_d);
      }
    }
    // else
    // {
    //   std::cout << "INFO: no data found in LUT for particle: charge: " <<
    //   charge
    //             << " kin energy: " << kinEnergy << " Skipping." << std::endl;
    // }
  }
}

// Called every time a Run ends
void GateAMDMActor::EndOfRunAction(const G4Run * /*unused*/) {

  // std::cout << "end of run action"
  //           << std::endl;
  // moved to end of simulation to avoid multiple threads
}

GateAMDMActor::ImageType4D::Pointer
GateAMDMActor::divide4dby3dImage(ImageType4D::Pointer imageToProcess,
                                 ImageType::Pointer denominatorImage) {
  ImageType4D::IndexType pixelIndex;
  pixelIndex[0] = 0; // x position of the pixel
  pixelIndex[1] = 0; // y position of the pixel
  pixelIndex[2] = 0; // z position of the pixel
  pixelIndex[3] = 0; // bin position of the pixel
  // std::cout << "in function before division: " <<
  // imageToProcess->GetPixel(pixelIndex) << std::endl;

  // define the extract filter
  using ExtractFilterType = itk::ExtractImageFilter<ImageType4D, ImageType>;
  auto extractFilter = ExtractFilterType::New();
  extractFilter->SetDirectionCollapseToSubmatrix();

  // set up the extraction region [one slice]
  ImageType4D::RegionType inputRegion = imageToProcess->GetBufferedRegion();

  ImageType4D::SizeType size = inputRegion.GetSize();
  ImageType4D::IndexType indexToProcess = inputRegion.GetIndex();
  // std::cout << "indexToProcess: " << indexToProcess << std::endl;
  // std::cout << "imageToProcess size: " << size << std::endl;
  int bins = size[3];

  for (int i = 0; i < bins; i++) {
    // ImageType::IndexType pixelIndex3D;
    // pixelIndex3D[0] = 0; // x position of the pixel
    // pixelIndex3D[1] = 0; // y position of the pixel
    // pixelIndex3D[2] = 0; // z position of the pixel

    // std::cout << "in function loop before division: " <<
    // imageToProcess->GetPixel(pixelIndex) << std::endl; std::cout << "in
    // function loop to divide with: " <<
    // denominatorImage->GetPixel(pixelIndex3D) << std::endl; pixelIndex[3] = i;
    // // bin position of the pixel for debugging output

    // std::cout << "i: " << i << std::endl;
    // process bin number i
    indexToProcess[3] = i;
    // std::cout << "indexToProcess: " << indexToProcess << std::endl;

    ImageType4D::RegionType desiredRegion;
    size[3] = 0; // we extract along 4th dimension

    desiredRegion.SetSize(size);
    desiredRegion.SetIndex(indexToProcess);

    // initialize the extraction region
    extractFilter->SetExtractionRegion(desiredRegion);
    extractFilter->SetDirectionCollapseToSubmatrix();
    extractFilter->SetInput(imageToProcess);
    // std::cout << "in function loop 3D before division: " <<
    // extractFilter->GetOutput()->GetPixel(pixelIndex3D) << std::endl;

    // the paste filter puts the extracted/processed slice back into its
    // original place
    using PasteFilterType = itk::PasteImageFilter<ImageType4D, ImageType>;
    auto pasteFilter = PasteFilterType::New();
    // pasteFilter->SetInPlace(false);
    // std::cout << "can in place " << pasteFilter->GetInPlace() << std::endl;

    // one 3d image is extracted and processed at a time
    using DivideImageFilterType =
        itk::DivideImageFilter<ImageType, ImageType, ImageType>;
    auto divideImageFilter = DivideImageFilterType::New();
    divideImageFilter->SetInput1(extractFilter->GetOutput());
    divideImageFilter->SetInput2(denominatorImage);
    divideImageFilter->Update();
    ImageType::Pointer divisionImage = divideImageFilter->GetOutput();
    // std::cout << "after division " << divisionImage->GetPixel(pixelIndex3D)
    // << std::endl;

    pasteFilter->SetSourceImage(divisionImage);
    pasteFilter->SetSourceRegion(divisionImage->GetLargestPossibleRegion());
    pasteFilter->SetDestinationImage(imageToProcess);
    pasteFilter->SetDestinationIndex(indexToProcess);
    pasteFilter->UpdateLargestPossibleRegion();
    pasteFilter->UpdateOutputInformation();
    pasteFilter->Update();

    imageToProcess = pasteFilter->GetOutput();
  }
  // std::cout << "at end of division function data at voxel: " <<
  // imageToProcess->GetPixel(pixelIndex) << std::endl;

  return imageToProcess;
}

GateAMDMActor::ImageType4D::Pointer
GateAMDMActor::divide4dImage(ImageType4D::Pointer imageToProcess,
                             ImageType4D::Pointer denominatorImage) {
  // ImageType4D::IndexType pixelIndex;
  // pixelIndex[0] = 0; // x position of the pixel
  // pixelIndex[1] = 0; // y position of the pixel
  // pixelIndex[2] = 0; // z position of the pixel
  // pixelIndex[3] = 0; // bin position of the pixel
  // std::cout << "in function before division: " <<
  // imageToProcess->GetPixel(pixelIndex) << std::endl;

  // one 4d image is divided by the other
  using DivideImageFilterType =
      itk::DivideImageFilter<ImageType4D, ImageType4D, ImageType4D>;
  auto divideImageFilter = DivideImageFilterType::New();
  divideImageFilter->SetInput1(imageToProcess);
  divideImageFilter->SetInput2(denominatorImage);
  divideImageFilter->Update();

  imageToProcess = divideImageFilter->GetOutput();
  // std::cout << "in function after division: " <<
  // imageToProcess->GetPixel(pixelIndex) << std::endl;
  return imageToProcess;
}

void GateAMDMActor::EndSimulationAction() {
  // std::cout << "end of simulation action"
  //           << std::endl;

  // post processing
  // divide gamma by delta

  // // create index for pixel for debugging output only
  // ImageType4D::IndexType pixelIndex;
  // pixelIndex[0] = 0; // x position of the pixel
  // pixelIndex[1] = 0; // y position of the pixel
  // pixelIndex[2] = 0; // z position of the pixel
  // pixelIndex[3] = 0; // bin position of the pixel

  // std::cout << "gamma before division: " <<
  // cpp_amdm_gamma_image->GetPixel(pixelIndex) << std::endl;

  std::string outputFileName;

  if (fstoreMergingData) {
    outputFileName =
        sOutputFileNameBase + "-unprocessedForMergingOnly-delta.mhd";
    // std::cout << "outputFileName: " << outputFileName << std::endl;
    GateAMDMActor::write4DImage(cpp_amdm_delta_image, outputFileName);
    outputFileName =
        sOutputFileNameBase + "-unprocessedForMergingOnly-gamma.mhd";
    // std::cout << "outputFileName: " << outputFileName << std::endl;
    GateAMDMActor::write4DImage(cpp_amdm_gamma_image, outputFileName);
  }

  // divide gamma by delta
  cpp_amdm_gamma_image =
      GateAMDMActor::divide4dImage(cpp_amdm_gamma_image, cpp_amdm_delta_image);
  // std::cout << "gamma after division: " <<
  // cpp_amdm_gamma_image->GetPixel(pixelIndex) << std::endl;

  // divide delta by restricted edep, every bin by the same edep

  // std::cout << "delta before division: " <<
  // cpp_amdm_delta_image->GetPixel(pixelIndex) << std::endl;
  cpp_amdm_delta_image = GateAMDMActor::divide4dby3dImage(
      cpp_amdm_delta_image, cpp_amdm_restricted_edep_image);
  // std::cout << "delta after division: " <<
  // cpp_amdm_delta_image->GetPixel(pixelIndex) << std::endl; write images to
  // file

  outputFileName = sOutputFileNameBase + "-delta.mhd";
  // std::cout << "outputFileName: " << outputFileName << std::endl;
  GateAMDMActor::write4DImage(cpp_amdm_delta_image, outputFileName);

  outputFileName = sOutputFileNameBase + "-gamma.mhd";
  GateAMDMActor::write4DImage(cpp_amdm_gamma_image, outputFileName);
}

std::string GateAMDMActor::remove_extension(const std::string &filename) {
  size_t lastdot = filename.find_last_of(".");
  if (lastdot == std::string::npos)
    return filename;
  return filename.substr(0, lastdot);
}

void GateAMDMActor::write4DImage(const ImageType4D::Pointer image,
                                 std::string filename) {
  using WriterType = itk::ImageFileWriter<ImageType4D>;
  auto writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  try {
    writer->Update();
  } catch (const itk::ExceptionObject &error) {
    std::cerr << "Error: " << error << std::endl;
    exit(EXIT_FAILURE);
  }
}
