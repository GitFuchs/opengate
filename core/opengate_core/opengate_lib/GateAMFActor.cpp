/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   ------------------------------------ -------------- */

#include "GateAMFActor.h"
#include "G4LinInterpolation.hh"
#include "G4Navigator.hh"
#include "G4RandomTools.hh"
#include "G4RunManager.hh"
#include "GateHelpers.h"
#include "GateHelpersDict.h"
#include "GateHelpersImage.h"
#include <itkImageRegionIterator.h>

#include "G4EmCalculator.hh"


#include "G4Deuteron.hh"
#include "G4Electron.hh"
#include "G4EmCalculator.hh"
#include "G4Gamma.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"


#include <cmath>
#include <filesystem>
// #include <iterator> 
// #include <map> 
// #include <vector>
// #include <iostream>

G4Mutex AMFMutex = G4MUTEX_INITIALIZER;


GateAMFActor::~GateAMFActor() {
}


GateAMFActor::GateAMFActor(py::dict &user_info) : GateVActor(user_info, true) {
  fActions.insert("SteppingAction");
  fActions.insert("BeginOfRunAction");
  fActions.insert("EndOfRunAction");
    }

void GateAMFActor::InitializeUserInfo(py::dict &user_info) {
  // IMPORTANT: call the base class method
  GateVActor::InitializeUserInfo(user_info);
    // Hit type (random, pre, post etc)
  fHitType = DictGetStr(user_info, "hit_type");

  fSpectraOutputFileName = DictGetStr(user_info, "microdosimetric_spectra_file_name");
  std::cout << "fOutputFileName: " << fSpectraOutputFileName << std::endl;

  // translation
  fTranslation = DictGetG4ThreeVector(user_info, "translation");

  fImageSpacing = DictGetG4ThreeVector(user_info, "spacing");
  fImageSize = DictGetG4ThreeVector(user_info, "size");
  std::cout << "fImageSize: " << fImageSize << std::endl;
  std::cout << "fImageSpacing: " << fImageSpacing << std::endl;  
  fTSEDfilename = DictGetStr(user_info, "tsed_file_name");
  std::cout << "fTSEDfilename: " << fTSEDfilename << std::endl;


    // DomainRadius = fPm->GetDoubleParameter(GetFullParmName("DomainRadius"), "Length");
    // DomainRadius
    // G4cout << "DomainRadius: " << DomainRadius / um << " um" << G4endl;
    // DomainRadius = DomainRadius / um;   // Convert DomainRadius from mm to um   
    fDomainRadius = 0.3; // in um
    double CelDiam = 2.0 * fDomainRadius;
    fNucleusRadius = 0.8 * fDomainRadius; // in um
    fBetaRef = 0.5 * fDomainRadius; // in um
        // Perform unit conversions 
    // DomainRadius = DomainRadius / um;   // Convert DomainRadius from mm to um 
    // NucleusRadius = NucleusRadius / um; // Convert NucleusRadius from mm to um 
    // BetaRef = BetaRef / (1. / (gray * gray)) ;             // Convert BetaRef to /Gy2
    calculator = new MicrodosimetricCalculator(nybin,
                              CelDiam,
                              fDomainRadius,
                              fNucleusRadius,
                              fBetaRef, iunit, mparased);
    calculator->setTSEDfilename(fTSEDfilename);

    calculator->setCalculationFlags(flinealEnergySpectra, fmeanLinealEnergy, fdoseAveragedLinealEnergy);
}

void GateAMFActor::InitializeCpp() {
  GateVActor::InitializeCpp();
  NbOfThreads = G4Threading::GetNumberOfRunningWorkerThreads();

  if (flinealEnergySpectra)
  {
          // Create the image pointers
        cpp_amf_microdosimetric_spectra = ImageVectorType::New();

        ImageVectorType::IndexType start;
        start[0] = 0; // first index on X
        start[1] = 0; // first index on Y
        start[2] = 0; // first index on Z

        ImageVectorType::SizeType size;
        size[0] = fImageSize[0];          // size along X
        size[1] = fImageSize[1];          // size along Y
        size[2] = fImageSize[2];          // size along Z

        ImageVectorType::SpacingType spacing;
        spacing[0] = fImageSpacing[0]; // spacing along X
        spacing[1] = fImageSpacing[1]; // spacing along Y
        spacing[2] = fImageSpacing[2]; // spacing along Z

        ImageVectorType::RegionType region;
        region.SetSize(size);
        region.SetIndex(start);

        cpp_amf_microdosimetric_spectra->SetRegions(region);
        cpp_amf_microdosimetric_spectra->SetSpacing(spacing);
        //image only contains histogram bins, not the labels, these are in histo_x_labels
        cpp_amf_microdosimetric_spectra->SetNumberOfComponentsPerPixel(nybin);

        cpp_amf_microdosimetric_spectra->Allocate();

        // Optionally initialize pixels
            ImageVectorType::PixelType pixel;
            pixel.SetSize(nybin);
            pixel.Fill(0.0);
            cpp_amf_microdosimetric_spectra->FillBuffer(pixel);

            //initialize histo_x_labels
            histo_x_labels.assign(nybin, 0.0);
  }

  if (fmeanLinealEnergy)
    {
            // Create the image pointers
            cpp_amf_mean_lineal_energy = Image3DType::New();
    }
  if (fdoseAveragedLinealEnergy)
        {
            // Create the image pointers
            cpp_amf_dose_averaged_lineal_energy = Image3DType::New();
        }
  // Create the image pointers
  // (the size and allocation will be performed on the py side)

    cpp_amf_dose_image = Image3DType::New();

//  std::cout << "AMF image size: " << size[0] << " " << size[1] << " " << size[2] << std::endl;
//  std::cout << "AMF image spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;
 std::cout << "End of InitializeCpp" << std::endl;
}

// void GateAMFActor::calculateDoseWeightedMicrodosimetricFunction(VectorPixelType& microDosSpectra, double izz, double iAA, double ene, double dEdx, double dose, double& LinealEnergy_Dose, double& LinealEnergyS) {
//     // This is now a wrapper that uses the optimized calculator
//     // You should create a member variable MicrodosimetricCalculator* calculator in the class
//     // and initialize it once in the constructor or InitializeCpp()
    
//     // For now, keeping original implementation as fallback
//     // TODO: Replace with calculator->calculateDoseWeightedMicrodosimetricFunction(...)
// }

// void GateAMFActor::getBinValueAndContent(const VectorPixelType& vec, 
//                                           size_t index, 
//                                           double& binValue, 
//                                           double& binContent) const {
//     if (index >= vec.Size() / 2) {
//         std::cerr << "Error: index out of range" << std::endl;
//         binValue = 0.0;
//         binContent = 0.0;
//         return;
//     }
//     binValue = vec[2 * index];
//     binContent = vec[2 * index + 1];
// }

// void GateAMFActor::setBinValueAndContent(VectorPixelType& vec, 
//                                           size_t index, 
//                                           double binValue, 
//                                           double binContent) {
//     if (index >= vec.Size() / 2) {
//         std::cerr << "Error: index out of range" << std::endl;
//         return;
//     }
//     vec[2 * index] = binValue;
//     vec[2 * index + 1] = binContent;
// }





G4double GateAMFActor::getDose(G4Step *step) {
  // get edep in MeV (take weight into account)
  auto w = step->GetTrack()->GetWeight();
//   auto edep = step->GetTotalEnergyDeposit() / CLHEP::MeV * w;
  auto edep = step->GetTotalEnergyDeposit() / joule * w;

  double dose;
  double density;

  auto *current_material = step->GetPreStepPoint()->GetMaterial();
  density = current_material->GetDensity()/(kg/mm3); // ensure density is in kg/mm3, by default in G4 it is in internal units
//   std::cout<< "Material name: " << current_material->GetName() << std::endl;
//   std::cout<< "Material formula: " << current_material->GetChemicalFormula() << std::endl;      
  dose = edep / (density*fVoxelVolume); // in Gy (J/kg)
  return dose;
}

G4double GateAMFActor::GetStoppingPower(G4Step *step) {
  // get the kinetic energy in MeV
  // G4double kinEnergy = step->GetTrack()->GetKineticEnergy() / (CLHEP::MeV);
  G4double kinEnergy = step->GetTrack()->GetKineticEnergy();

  auto* track = step->GetTrack();                      // the track associated with this step
  auto* particle_definition = track->GetParticleDefinition();  // G4ParticleDefinition*
  auto* mat    = step->GetPreStepPoint()->GetMaterial();

  G4EmCalculator emcalc;

  auto total_dEdx = emcalc.ComputeTotalDEDX(kinEnergy, particle_definition, mat);              // MeV*cm2/g
  return total_dEdx;
}

void GateAMFActor::GetVoxelPosition(G4Step *step, G4ThreeVector &position,
                                     bool &isInside,
                                     Image3DType::IndexType &index) const {
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();

  // consider random position between pre and post
  if (fHitType == "pre") {
    position = preGlobal;
  }
  if (fHitType == "random") {
    auto x = G4UniformRand();
    auto direction = postGlobal - preGlobal;
    position = preGlobal + x * direction;
  }
  if (fHitType == "middle") {
    auto direction = postGlobal - preGlobal;
    position = preGlobal + 0.5 * direction;
  }

  auto localPosition =
      touchable->GetHistory()->GetTransform(0).TransformPoint(position);

  // convert G4ThreeVector to itk PointType
  Image3DType::PointType point;
  point[0] = localPosition[0];
  point[1] = localPosition[1];
  point[2] = localPosition[2];

  isInside = cpp_amf_dose_image->TransformPhysicalPointToIndex(point, index);
}



void GateAMFActor::BeginOfRunAction(const G4Run *) {

  std::cout << "AMF actor starting run BeginOfRunActionMasterThread"
  << std::endl;
    //   std::cout << "fPhysicalVolumeName: " << fPhysicalVolumeName << std::endl;
    //   std::cout << "fInitialTranslation: " << fTranslation << std::endl;  

    if (flinealEnergySpectra){
        AttachImageToVolume<ImageVectorType>(cpp_amf_microdosimetric_spectra, fPhysicalVolumeName,
                                    fTranslation);
    }
    if (fmeanLinealEnergy){
    AttachImageToVolume<Image3DType>(cpp_amf_mean_lineal_energy, fPhysicalVolumeName,
                                    fTranslation);
    }
    if (fdoseAveragedLinealEnergy){
        AttachImageToVolume<Image3DType>(cpp_amf_dose_averaged_lineal_energy, fPhysicalVolumeName,
                                    fTranslation);
    }

        // Create the image pointers}
      // Important ! The volume may have moved, so we re-attach each run
  AttachImageToVolume<Image3DType>(cpp_amf_dose_image, fPhysicalVolumeName,
                                   fTranslation);


  auto sp = cpp_amf_dose_image->GetSpacing();
  fVoxelVolume = sp[0] * sp[1] * sp[2];
//   std::cout << "Voxel spacing: " << sp << " mm" << std::endl;
//   std::cout << "Voxel volume: " << fVoxelVolume << " mm3" << std::endl;
  std::cout << "end of BeginOfRunActionMasterThread"
  << std::endl;
}

void GateAMFActor::SteppingAction(G4Step *step) {
    //  std::cout << "Begin of SteppingAction" << std::endl;

  auto event_id =
      G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();

  double dEdx, LinealEnergy_Dose, LinealEnergyS;

  G4double dose = GateAMFActor::getDose(step);

  // Get the voxel index
  G4ThreeVector position;
  bool isInside;
  Image3DType::IndexType index;
  GetVoxelPosition(step, position, isInside, index);
//   std::cout << "Voxel index: " << index << std::endl;
//   std::cout << "position: " << position << std::endl;
//   std::cout << "dose: " << dose << std::endl;
//   std::cout << "isInside: " << isInside << std::endl;

    // If the position is not inside the image, return
  if (!isInside)
    return;

    if (dose > 0.) {
        G4double density = step->GetPreStepPoint()->GetMaterial()->GetDensity();
        G4double iAA = step->GetTrack()->GetDefinition()->GetAtomicMass();
        G4double kenergy = step->GetPreStepPoint()->GetKineticEnergy();
        G4double energyPerNucleon = kenergy / iAA; // energy per nucleon
        G4double izz = step->GetTrack()->GetDefinition()->GetAtomicNumber();

        if (izz >= 1 && izz <= 18 && energyPerNucleon >= 0.025) {
            // G4int G4binIndex = GetIndex(step);

            dEdx = GetStoppingPower(step);
            LinealEnergy_Dose = 0.0;
            LinealEnergyS = 0.0;
            
            if (flinealEnergySpectra || fdoseAveragedLinealEnergy || fmeanLinealEnergy){ 
                calculator->calculateDoseWeightedMicrodosimetricFunction(microdosimetricSpectra, izz, iAA, energyPerNucleon, dEdx, dose, LinealEnergy_Dose, LinealEnergyS);
            }

            // calculateDoseWeightedMicrodosimetricFunction(microdosimetricSpectra, izz, iAA, energyPerNucleon, dEdx, dose, LinealEnergy_Dose, LinealEnergyS);
            // VectorPixelType microdosimetricSpectra;
            //     microdosimetricSpectra.SetSize(nybin);
            //     microdosimetricSpectra.Fill(0.0);

            
            // LinealEnergyS=12.2;
            // std::cout << "dEdx: " << dEdx << std::endl;
            // std::cout << "LinealEnergy_Dose: " << LinealEnergy_Dose << std::endl;
            // std::cout << "LinealEnergyS: " << LinealEnergyS << std::endl;
            // std::cout << "dose: " << dose << std::endl;
            // std::cout << "Voxel index: " << index << std::endl;
            // if (LinealEnergyS>0){
            //     std::cout << "LinealEnergyS: " << LinealEnergyS << std::endl;
            //     std::cout << "Voxel index: " << index << std::endl;

            // }
  
            if (flinealEnergySpectra){
                ImageAddValue<ImageVectorType>(cpp_amf_microdosimetric_spectra, index, microdosimetricSpectra);
            }
            if (fdoseAveragedLinealEnergy){
                ImageAddValue<Image3DType>(cpp_amf_dose_averaged_lineal_energy, index, LinealEnergyS);

            }
            if (fmeanLinealEnergy){
                ImageAddValue<Image3DType>(cpp_amf_mean_lineal_energy, index, LinealEnergy_Dose);
            }
            ImageAddValue<Image3DType>(cpp_amf_dose_image, index, dose);


            // std::cout<<"Voxel index: " << index << std::endl;

            // auto pixelValue = cpp_amf_microdosimetric_spectra->GetPixel(index);
            // std::cout << "Pixel value: " << pixelValue << std::endl;
            // std::cout << "End of SteppingAction" << std::endl;

            return ;
        }
    }
        //  std::cout << "End of SteppingAction" << std::endl;
    
    return ;
}


int GateAMFActor::EndOfRunActionMasterThread(int run_id)
{
      std::cout << "begin of EndOfRunActionMasterThread"
            << std::endl;
    return 0;}

//   std::cout << "AMF actor ending run BeginOfRunActionMasterThread"
//   << std::endl;

//     writeVectorImage(cpp_amf_microdosimetric_spectra, fSpectraOutputFileName);
//     return 0;
// }

  // Called every time a Run ends (all threads)
void GateAMFActor::EndOfRunAction(const G4Run *run)
{
    // // Run the final, once-per-run action only on master
    // if (IsMaster()) {
    //     std::cout << "begin of EndOfRunAction MasterThread detected"
    //         << std::endl;
    //   writeVectorImage(cpp_amf_microdosimetric_spectra, fSpectraOutputFileName);

    //   // If using G4AnalysisManager:
    //   // auto am = G4AnalysisManager::Instance();
    //   // am->Write();
    //   // am->CloseFile();
    // } else {
    //   std::cout << "begin of EndOfRunAction WorkerThread"
    //         << std::endl;
    // }

    std::cout << "begin of EndOfRunAction" << std::endl;
    {
    G4AutoLock mutex(&AMFMutex);
    {
        if (flinealEnergySpectra){
            if (!fRanOnce) {

                writeVectorImage(cpp_amf_microdosimetric_spectra, fSpectraOutputFileName);

                // Remove file extension from fSpectraOutputFileName
                namespace fs = std::filesystem;
                fs::path p(fSpectraOutputFileName);
                std::string baseFilename = (p.parent_path() / p.stem()).string();
                std::string histo_x_labels_filename = baseFilename + "_histo_x_labels.txt";
                writeVectorToTextFile(histo_x_labels, histo_x_labels_filename, "#Histogram x-axis labels (lineal energy in keV/um)");
                fRanOnce = true;
            }
        }
    }
    // Image3DType::IndexType index;
    // index[0] = 0;
    // index[1] = 0;
    // index[2] = 0;
    // auto pixelValue = cpp_amf_microdosimetric_spectra->GetPixel(index);
    // std::cout << "Pixel value: " << pixelValue << std::endl;
    // std::cout << "cpp_amf_microdosimetric_spectra: " <<cpp_amf_microdosimetric_spectra->GetNumberOfComponentsPerPixel() << std::endl;


    }
    std::cout << "end of EndOfRunAction" << std::endl;
}

void GateAMFActor::writeVectorImage(const ImageVectorType::Pointer image,
                                 std::string filename)
{
  using WriterType = itk::ImageFileWriter<ImageVectorType>;
  auto writer = WriterType::New();
  std::cout << "Writing Vector image to " << filename << std::endl;
  writer->SetFileName(filename);
  writer->SetInput(image);
  try
  {
    writer->Update();
  }
  catch (const itk::ExceptionObject &error)
  {
    std::cerr << "Error: " << error << std::endl;
    exit(EXIT_FAILURE);
  }
}

void GateAMFActor::writeVectorToTextFile(const std::vector<double> &vec,
                                         const std::string &filename,
                                         const std::string &header) {
  std::ofstream outFile(filename);
  if (!outFile.is_open()) {
    std::cerr << "Failed to open " << filename << " for writing." << std::endl;
    return;
  }
  
  if (!header.empty()) {
    outFile << header << std::endl;
  }
  
  for (const auto &value : vec) {
    outFile << value << std::endl;
  }
  
  outFile.close();
}

  // Prefer this helper (true for master in MT, true in sequential mode as well)
G4bool GateAMFActor::IsMaster() const {
    return G4Threading::IsMasterThread();
  }




  

// New optimized class for microdosimetric calculations
// Public reinitialization function if parameters need to be updated
void MicrodosimetricCalculator::reinitialize(size_t nybin_val, double celDiam, double domainRadius, 
                    double nucleusRadius, double betaRef, int iunit_val, int mparased_val) {
    nybin = nybin_val;
    CelDiam = celDiam;
    fDomainRadius = domainRadius;
    fNucleusRadius = nucleusRadius;
    fBetaRef = betaRef;
    iunit = iunit_val;
    mparased = mparased_val;
    initialize();
}

void MicrodosimetricCalculator::initialize() {
    ypower = -3.0;  // Initialize here instead
    
    std::fill(Apara, Apara + 9, 0.0);

    yhig.resize(nybin + 1);
    yfy.resize(nybin);
    ydy.resize(nybin);

    ymid.resize(nybin);
    ywid.resize(nybin);
    eventmid.resize(nybin);
    Z.resize(nybin);
    histo_x_labels.resize(nybin);

    for (size_t i = 0; i < yhig.size(); ++i) {
        yhig[i] = std::pow(10.0, ypower);
        ypower += ystep;
    }

    if (iunit <= 1) {
        unitconv = 1.0;
    } else if (iunit == 2) {
        unitconv = 1.0e-3 * (2.0 / 3.0 * CelDiam);
    } else if (iunit == 3) {
        unitconv = 4.0 / 3.0 * M_PI * std::pow(CelDiam / 2.0, 3) * 1.0e-15 / 1.602e-13;
    }

    // Calculate y0 first as it's needed in the loop
    y0 = (M_PI * fDomainRadius * std::pow(fNucleusRadius, 2)) / (std::sqrt(fBetaRef * (std::pow(fDomainRadius, 2) + std::pow(fNucleusRadius, 2))) * 0.16022);

    for (size_t i = 0; i < nybin; ++i) {
        double ymid_val = (yhig[i] + yhig[i + 1]) / 2.0;
        ymid[i] = ymid_val;
        ywid[i] = yhig[i + 1] - yhig[i];
        eventmid[i] = ymid_val * factor * unitconv;
        Z[i] = 1 - std::exp(-std::pow(ymid_val, 2) / std::pow(y0, 2));
        histo_x_labels[i] = ymid_val;
    }

    // Calculate the bins per decade
    binsperDecade = nybin / (std::log10(yhig.back() / yhig[0]));
}


    // Constructor with all parameters
    MicrodosimetricCalculator::MicrodosimetricCalculator(size_t nybin_val, double celDiam, double domainRadius, 
                             double nucleusRadius, double betaRef, int iunit_val, int mparased_val)
        : nybin(nybin_val), CelDiam(celDiam), fDomainRadius(domainRadius),
          fNucleusRadius(nucleusRadius), fBetaRef(betaRef), iunit(iunit_val), mparased(mparased_val),
          factor(0.0), unitconv(0.0), binsperDecade(0.0), y0(0.0), ypower(-3.0) {
        initialize();
    }


    
void MicrodosimetricCalculator::get_Histo_X_Labels(std::vector<double>& labels) const {
    labels = histo_x_labels;
}

void MicrodosimetricCalculator::calculateDoseWeightedMicrodosimetricFunction(VectorPixelType& microDosSpectra, double izz, double iAA, double ene, double dEdx, double dose, double& LinealEnergy_Dose, double& LinealEnergyS) {
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;

    if (microDosSpectra.Size() != nybin) {
        microDosSpectra.SetSize(nybin);
        microDosSpectra.Fill(0.0);
    }
    else {
        microDosSpectra.Fill(0.0);
    }
    int ic1, ie1, ip1;
    double ratioc, ratioe, ratiop;

    double erg = ene * iAA;
    double depev = std::min(dEdx * CelDiam * 1.0e3, erg * 1.0e6);

    getAparaion(CelDiam, ene, iAA, izz, ratioc, ratioe, ratiop, ic1, ie1, ip1);    
    sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
    
    factor = (iunit == 0) ? 1.0 : 1.0e6 / Apara[8];
    for (size_t i = 0; i < nybin; ++i) {
        yfy[i] = ymid[i] * sedmean(eventmid[i], depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
        ydy[i] = yfy[i] * ymid[i];
        sum0 += yfy[i] * ywid[i] / ymid[i];
        sum1 += yfy[i] * ywid[i];
        sum2 += yfy[i] * ywid[i] * ymid[i];
    }


    if (flinealEnergySpectra){
        double normalization_factor = (binsperDecade / std::log(10)) / std::accumulate(ydy.begin(), ydy.end(), 0.0);

        for (size_t i = 0; i < ydy.size(); ++i) {
            ydy[i] *= normalization_factor;
            microDosSpectra[i] = ydy[i]*dose;
        }
    }

    if (fmeanLinealEnergy){
        if (sum1 == 0) {
            std::cout << "Warning: sum1 is zero, returning zero vector." << std::endl;
            LinealEnergy_Dose = 0.0;
            return;
        }
        LinealEnergy_Dose = sum2/sum1;
    }

    if (fdoseAveragedLinealEnergy){
        double LinealEnergy_Freq = sum1 / sum0;
        
        double sumNumerator = 0.0;
        double sumDenominator = 0.0;

        for (size_t i = 0; i < nybin; ++i) {
            sumNumerator += yfy[i] * Z[i];
            sumDenominator += yfy[i];
        }
        LinealEnergyS = 0.0;
        if (sumDenominator > 0.0) {
            LinealEnergyS = ((sumNumerator / sumDenominator) / LinealEnergy_Freq ) * std::pow(y0, 2);
        }
    }

    return;
}


void MicrodosimetricCalculator::getAparaion(const double& CelDiam, const double& ene, const int& iAA, const int& izz, double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1) {
    double erg = ene * iAA;
    int modifiedIzz = (izz > 26) ? 26 : izz;
    double erg_AA = erg / iAA;
    double CD = std::abs(CelDiam);

    int ic = 0;
    for (ic = 1; ic <= 9; ++ic) {
        if (cdiamion[ic - 1] >= CD) {
            break;
        }
    }

    if (ic == 1) {
        ic1 = 1;
        ratioc = 0.0;
    } else {
        ic1 = ic - 1;
        ratioc = std::min(1.0, (std::log10(CD) - std::log10(cdiamion[ic1-1])) / (std::log10(cdiamion[ic1]) - std::log10(cdiamion[ic1-1])));
    }

    int ie = 0;
    for (ie = 1; ie <= 12; ++ie) {
        if (eincion[ie-1] >= erg_AA) {
            break;
        }
    }

    if (ie == 1) {
        ie1 = 1;
        ratioe = 0.0;
    } else {
        ie1 = ie - 1;
        ratioe = std::min(1.0, (std::log10(erg_AA) - std::log10(eincion[ie1-1])) / (std::log10(eincion[ie1]) - std::log10(eincion[ie1-1])));
    }

    int ip = 1;
    for (ip = 2; ip <= 6; ++ip) {
        if (izion[ip-1] >= izz) {
            break;
        }
    }

    ip1 = ip - 1;
    ratiop = std::min(1.0, static_cast<double>(izz - izion[ip1-1]) / (izion[ip1] - izion[ip1-1]));
}


double MicrodosimetricCalculator::sedfunc(double x, double depev, const double Apara[], size_t size) {
    double getfirst = 0.0, getsecond = 0.0, getthird = 0.0;

    if (Apara[0] > 0.0) {
        double tmp;
        if (depev == 0.0) {
            tmp = std::min(50.0, std::pow(std::abs(x - Apara[1]), Apara[2]) / (2 * Apara[1]));
            getfirst = Apara[0] * std::exp(-tmp);
        } else {
            double cst1 = depev / Apara[8];
            tmp = std::min(50.0, Apara[1] * (x - cst1 * Apara[2]));
            getfirst = Apara[0] * x / (std::exp(tmp) + 1) * (2.0 / std::pow(cst1 * Apara[2], 2));
        }
    }

    if (Apara[3] > 0.0) {
        double tmp = std::min(50.0, std::pow(std::abs(x - Apara[4]), Apara[5]) / (2 * Apara[4]));
        getsecond = Apara[3] * std::exp(-tmp);
    }

    if (Apara[6] > 0.0) {
        getthird = Apara[6] / (Apara[7] - 1.0) * std::pow((Apara[7] - 1.0) / Apara[7], x);
    }

    double sedfunc = getfirst + getsecond + getthird;
    if (sedfunc < 1.0e-10) {
        sedfunc = 0.0;
    }
    return sedfunc;
}

double MicrodosimetricCalculator::sedmean(double x, double depev, int ic1, int ie1, int ip1, double ratioc, double ratioe, double ratiop, double Apara[]) {
    double sedmean = 0.0;
    double A9 = 0.0;

    for (int ip = ip1; ip <= ip1 + 1; ip++) {
        double Rp = (ip == ip1) ? (1.0 - ratiop) : ratiop;
        for (int ie = ie1; ie <= ie1 + 1; ie++) {
            double Re = (ie == ie1) ? (1.0 - ratioe) : ratioe;
            for (int ic = ic1; ic <= ic1 + 1; ic++) {
                double Rc = (ic == ic1) ? (1.0 - ratioc) : ratioc;
                int index = ((ip-1) * 96) + ((ie-1) * 8) + (ic-1);
                for (int i = 0; i < mparased; i++) {
                    Apara[i] = IonData[index][i];
                }
                double wei = Rp * Re * Rc;
                double sedfuncResult = sedfunc(x, depev, Apara, mparased);
                sedmean += sedfuncResult * wei;
                A9 += Apara[8] * wei;
            }
        }
    }

    Apara[8] = A9;
    return sedmean;
}

void MicrodosimetricCalculator::setTSEDfilename(const std::string& filename) {
    fTSEDfilename = filename;
    loadIonData();
}

void MicrodosimetricCalculator::loadIonData() {
    IonData.resize(ROWS, std::vector<double>(COLS));
    std::ifstream file(fTSEDfilename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << fTSEDfilename << " for reading. Aborting." << std::endl;
        exit(1);
        // return;
    }

    std::string line;
    int row = 0;

    while (std::getline(file, line) && row < ROWS) {
        std::istringstream iss(line);
        double value;
        int col = 0;

        while (iss >> value && col < COLS) {
            IonData[row][col] = value;
            col++;
        }
        row++;
    }

    file.close();
}






