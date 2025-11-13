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
// #include <iterator> 
// #include <map> 
// #include <vector>
// #include <iostream>

namespace py = pybind11;

class MicrodosimetricCalculator;

class GateAMFActor : public GateVActor {

public:
  // Image type is 4D float by default
  typedef itk::Image<double, 3> Image3DType;

  // typedef itk::VariableLengthVector<double> VectorPixelType;
  // typedef itk::Image<VectorPixelType, 3> ImageVectorType;

  typedef itk::VectorImage<double, 3> ImageVectorType;
  typedef ImageVectorType::PixelType VectorPixelType;


  // Constructor
  GateAMFActor(py::dict &user_info);
  ~GateAMFActor();
  void BeginOfRunAction(const G4Run *);
  void SteppingAction(G4Step *step);
  int EndOfRunActionMasterThread(int run_id) override;
    // Called every time a Run ends (all threads)
  void EndOfRunAction(const G4Run *run) override;
  // void EndOfEventAction(const G4Event *event);

  void InitializeUserInfo(py::dict &user_info) override;

  // void calculateDoseWeightedMicrodosimetricFunction(VectorPixelType& microDosSpectra,double izz, double iAA, double ene, double dEdx, double dose,double& LinealEnergy_Dose, double& LinealEnergyS);
  // void getBinValueAndContent(const VectorPixelType& vec, 
  //                                         size_t index, 
  //                                         double& binValue, 
  //                                         double& binContent)const ;

  // void setBinValueAndContent(VectorPixelType& vec, 
  //                                         size_t index, 
  //                                         double binValue, 
  //                                         double binContent);
  G4double getDose(G4Step *step);
  G4double GetStoppingPower(G4Step *step);
  // void getAparaion(const double& CelDiam, const double& ene, const int& iAA, const int& izz, double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1);
  // double sedmean(double x, double depev, int ic1, int ie1, int ip1, double ratioc, double ratioe, double ratiop, double Apara[]);
  // double sedfunc(double x, double depev, const double Apara[], size_t size);
  // void loadIonData();
  G4bool IsMaster() const;

  std::string GetPhysicalVolumeName() const { return fPhysicalVolumeName; }
  void writeVectorImage(const ImageVectorType::Pointer image,
                                 std::string filename);
  void writeVectorToTextFile(const std::vector<double> &vec,
                                         const std::string &filename,
                                         const std::string &header);                             

  void InitializeCpp();
  void SetPhysicalVolumeName(std::string s) { fPhysicalVolumeName = s; }
  void GetVoxelPosition(G4Step *step, G4ThreeVector &position, bool &isInside,
                        Image3DType::IndexType &index) const;
                    
  bool GetLinealEnergySpectraFlag() const { return flinealEnergySpectra; }
  void SetLinealEnergySpectraFlag(const bool b) { flinealEnergySpectra = b; }

  bool GetMeanLinealEnergyFlag() const { return fmeanLinealEnergy; }
  void SetMeanLinealEnergyFlag(const bool b) { fmeanLinealEnergy = b; }

  bool GetDoseAveragedLinealEnergyFlag() const { return fdoseAveragedLinealEnergy; }
  void SetDoseAveragedLinealEnergyFlag(const bool b) { fdoseAveragedLinealEnergy = b; }



  // void EndSimulationAction();

  // Image3DType::SizeType size_edep{};
  std::vector<std::vector<double>> IonData;
  std::map<G4int, std::vector<G4double>> totalSpectra; // Key: bin index, Value: spectra (vector of 180 values)
  std::map<G4int, G4double> cumulativeDose; 

  //  The image is accessible on py side (shared by all threads)
  Image3DType::Pointer cpp_amf_dose_image;
  ImageVectorType::Pointer cpp_amf_microdosimetric_spectra;
  Image3DType::Pointer cpp_amf_mean_lineal_energy;
  Image3DType::Pointer cpp_amf_dose_averaged_lineal_energy;


  double fVoxelVolume{};
  std::string fPhysicalVolumeName;
  std::string fHitType;
  std::string fTSEDfilename;
  std::string fSpectraOutputFileName;
  int NbOfThreads = 0;
  G4ThreeVector fImageSize;
  G4ThreeVector fImageSpacing;
  double fDomainRadius;
  int NbOfEvent = 0;
  double fNucleusRadius;
  double fBetaRef;
  bool fRanOnce = false;
  std::vector<double> histo_x_labels;
  VectorPixelType microdosimetricSpectra;

  MicrodosimetricCalculator* calculator;

  // Option: indicate we must calculate lineal energy spectra
  bool flinealEnergySpectra{};
  // Option: indicate we must calculate mean lineal energy
  bool fmeanLinealEnergy{};
  // Option: indicate we must calculate dose-averaged lineal energy
  bool fdoseAveragedLinealEnergy{};



  G4ThreeVector fTranslation;

private:
  std::vector<double> yhig, yfy, ydy;
  static constexpr int nybin = 400;
  static constexpr int mparased = 9;
  static constexpr int iunit = 2;
  // static constexpr int ROWS = 576;
  // static constexpr int COLS = 9;
    // Define constants outside the class
	// const std::vector<double> eincion = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 100.0, 300.0, 999.0};
	// const std::vector<double> cdiamion = {0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0};
	// const std::vector<int> izion = {1, 2, 6, 10, 14, 26};
};


class MicrodosimetricCalculator {
  private:
    using ImageVectorType = itk::VectorImage<double, 3>;
    using VectorPixelType = ImageVectorType::PixelType;

    // Member variables
    size_t nybin{};
    double CelDiam{};
    double fDomainRadius{};
    double fNucleusRadius{};
    double fBetaRef{};
    int iunit{};
    int mparased{};
    double factor{};
    std::string fTSEDfilename;


    std::vector<double> yhig, yfy, ydy;
    std::vector<double> ymid, ywid, eventmid, Z;
    std::vector<std::vector<double>> IonData;

    double Apara[9];  // Fixed size based on mparased constant

    static constexpr int ROWS = 576;
    static constexpr int COLS = 9;

    // Constants (class-level)
    // static const std::array<double, 12> kEinCion;
    // static const std::array<double, 8>  kCdiamion;
    // static const std::array<int, 6>     kIzion;
    const std::vector<double> eincion = {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 100.0, 300.0, 999.0};
	  const std::vector<double> cdiamion = {0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0};
	  const std::vector<int> izion = {1, 2, 6, 10, 14, 26};

    std::vector<double> histo_x_labels;

    double ypower{};            // initialized in initialize()
    static constexpr double ystep = 0.02;
    double unitconv{};
    double binsperDecade{};
    double y0{};

    // Option: indicate we must calculate lineal energy spectra
    bool flinealEnergySpectra=true;
    // Option: indicate we must calculate mean lineal energy
    bool fmeanLinealEnergy=true;
    // Option: indicate we must calculate dose-averaged lineal energy
    bool fdoseAveragedLinealEnergy=true;

    // Private initialization function
    void initialize();

public:
    // Pixel type alias exposed (if callers need it)
    using PixelVectorType = VectorPixelType;

    // Constructor
    MicrodosimetricCalculator(size_t nybin_val, double celDiam, double domainRadius, 
                             double nucleusRadius, double betaRef, int iunit_val, int mparased_val);

    // Public reinitialization function if parameters need to be updated
    void reinitialize(size_t nybin_val,
                      double celDiam,
                      double domainRadius,
                      double nucleusRadius,
                      double betaRef,
                      int iunit_val,
                      int mparased_val);

    void get_Histo_X_Labels(std::vector<double>& labels) const;

    void calculateDoseWeightedMicrodosimetricFunction(
        VectorPixelType& microDosSpectra,
        double izz, double iAA, double ene, double dEdx, double dose,
        double& LinealEnergy_Dose, double& LinealEnergyS);

    // These methods must be implemented/linked by the user
    void getAparaion(const double& CelDiam, const double& ene, const int& iAA, const int& izz,
                     double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1);

    double sedmean(double x, double depev, int ic1, int ie1, int ip1,
                   double ratioc, double ratioe, double ratiop, double Apara[]);
    double sedfunc(double x, double depev, const double Apara[], size_t size);
    void setTSEDfilename(const std::string& filename);
    void loadIonData();
    void setCalculationFlags(bool linealEnergySpectra, bool meanLinealEnergy, bool doseAveragedLinealEnergy) {
        flinealEnergySpectra = linealEnergySpectra;
        fmeanLinealEnergy = meanLinealEnergy;
        fdoseAveragedLinealEnergy = doseAveragedLinealEnergy;
    } 


};



#endif // GateBeamQualityActor_h
