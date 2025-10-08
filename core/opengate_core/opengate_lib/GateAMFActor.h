/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateAMFActor_h
#define GateAMFActor_h

#include "G4Cache.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4DataVector.hh"

#include "GateHelpersImage.h"
#include "GateVActor.h"

#include "itkImage.h"
#include "itkImageFileWriter.h"

#include <pybind11/stl.h>
#include <cmath>


namespace py = pybind11;

class GateAMFActor : public GateVActor {

public:
  // Image type is 4D float by default
  typedef itk::Image<double, 3> Image3DType;
  typedef itk::Image<double, 5> Image5DType;

  // Constructor
  GateAMFActor(py::dict &user_info);
  ~GateAMFActor();
  void BeginOfRunAction(const G4Run *);
  void SteppingAction(G4Step *step);
  // void EndOfEventAction(const G4Event *event);

  void InitializeUserInfo(py::dict &user_info) override;

  std::vector<std::pair<double, double>> calculateMicrodosimetricFunction(double izz, double iAA, double ene, double dEdx);
  G4double getDose(G4Step *step);
  G4double GetStoppingPower(G4Step *step);
  void getAparaion(const double& CelDiam, const double& ene, const int& iAA, const int& izz, double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1);
  double sedmean(double x, double depev, int ic1, int ie1, int ip1, double ratioc, double ratioe, double ratiop, double Apara[]);
  double sedfunc(double x, double depev, const double Apara[], size_t size);
  void loadIonData();

  std::string GetPhysicalVolumeName() const { return fPhysicalVolumeName; }

  void InitializeCpp();
  void SetPhysicalVolumeName(std::string s) { fPhysicalVolumeName = s; }
  void GetVoxelPosition(G4Step *step, G4ThreeVector &position, bool &isInside,
                        Image3DType::IndexType &index) const;


  // void EndSimulationAction();

  // Image3DType::SizeType size_edep{};
  double CelDiam = 0.6; 
  static constexpr int nybin = 400;
  std::vector<std::vector<double>> IonData;
  std::map<G4int, std::vector<G4double>> totalSpectra; // Key: bin index, Value: spectra (vector of 180 values)
  std::map<G4int, G4double> cumulativeDose; 

  //  The image is accessible on py side (shared by all threads)
  Image3DType::Pointer cpp_amf_dose_image;
  Image5DType::Pointer cpp_amf_microdosimetric_spectra;
  Image3DType::Pointer cpp_amf_mean_lineal_energy;
  Image3DType::Pointer cpp_amf_dose_averaged_lineal_energy;


  double fVoxelVolume{};
  std::string fPhysicalVolumeName;
  std::string fHitType;
  int NbOfThreads = 0;
  G4ThreeVector fImageSize;
  G4ThreeVector fImageSpacing;



  G4ThreeVector fTranslation;

private:
  std::vector<double> yhig, yfy, ydy;
  static constexpr int mparased = 9;
  static constexpr int iunit = 2;
  static constexpr int ROWS = 576;
  static constexpr int COLS = 9;
    // Define constants outside the class
	const std::vector<double> eincion = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 100.0, 300.0, 999.0};
	const std::vector<double> cdiamion = {0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0};
	const std::vector<int> izion = {1, 2, 6, 10, 14, 26};
};

#endif // GateBeamQualityActor_h
