/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateAMDMActor_h
#define GateAMDMActor_h

#include "G4VPrimitiveScorer.hh"
#include "GateAMDMActor_lookUpTable.hh"
#include "GateVActor.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

class GateAMDMActor : public GateVActor {

public:
  // Image type is 4D float by default
  typedef itk::Image<float, 3> ImageType;
  typedef itk::Image<float, 4> ImageType4D;

  // Constructor
  GateAMDMActor(py::dict &user_info);

  virtual void ActorInitialize(py::dict &user_info);

  // Main function called every step in attached volume
  virtual void SteppingAction(G4Step *);

  // Called every time a Run starts (all threads)
  virtual void BeginOfRunAction(const G4Run *run);

  virtual void EndOfRunAction(const G4Run *);

  virtual void EndSimulationAction();

  ImageType4D::Pointer
  divide4dImage(const ImageType4D::Pointer imageToProcess,
                const ImageType4D::Pointer denominatorImage);
  ImageType4D::Pointer
  divide4dby3dImage(const ImageType4D::Pointer imageToProcess,
                    const ImageType::Pointer denominatorImage);
  std::string remove_extension(const std::string &filename);

  // typedef itk::Vector<float, 10> PixelType;

  // typedef itk::VariableLengthVector<float> PixelType;
  // using PixelType = itk::VariableLengthVector<double>;
  // typedef itk::VectorImage<PixelType, 3> VectorImageType;
  // using PixelType = itk::Vector<float, 10>;

  // using VectorImageType = itk::Image<PixelType, 3>;

  // The image is accessible on py side (shared by all threads)
  ImageType::Pointer cpp_amdm_restricted_edep_image;
  ImageType4D::Pointer cpp_amdm_delta_image;
  ImageType4D::Pointer cpp_amdm_gamma_image;

  // VectorImageType::Pointer cpp_amdm_delta_image;
  // VectorImageType::Pointer cpp_amdm_gamma_image;

  const int kMaxAMDMBins = 10;

  // // Option: indicate if we must compute uncertainty
  // bool fUncertaintyFlag;

  // Option: indicate if we must compute dose in Gray also
  // bool fGrayFlag;

  // For normalization, we need temporary images
  // ImageType::Pointer cpp_restricted_edep_image;
  // ImageType::Pointer cpp_temp_image;
  // ImageType::Pointer cpp_last_id_image;
  // ImageType::Pointer cpp_dose_image;
  // double fVoxelVolume;

  std::string fPhysicalVolumeName;

  G4ThreeVector fInitialTranslation;
  G4ThreeVector fImageSize;
  G4ThreeVector fImageSpacing;
  std::string fHitType;
  std::string fLUTfilename;
  std::string sOutputFileNameBase;

  lookUpTable *famdm_tables;
};

#endif // GateAMDMActor_h
