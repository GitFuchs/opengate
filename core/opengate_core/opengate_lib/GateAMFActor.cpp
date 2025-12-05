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

  std::cout << "fdomainRadiusInUm: "<< fdomainRadiusInUm <<std::endl;

    double CelDiam = 2.0 * fdomainRadiusInUm; // in um, fDomainRadiusInUm is in um
    double nucleusRadius = 0.8 * fdomainRadiusInUm; // in um
    double betaRef = 0.5 * fdomainRadiusInUm; // in um

    calculator = new MicrodosimetricCalculator(nybin,
                              CelDiam,
                              fdomainRadiusInUm,
                              nucleusRadius,
                              betaRef, iunit, mparased);
    calculator->setTSEDfilename(fTSEDfilename);

    
    // std::cout << "InitializeUserInfo - Reinitializing MicrodosimetricCalculator with parameters:" << std::endl;
    // std::cout << "  nybin: " << nybin << std::endl;
    // std::cout << "  CelDiam: " << CelDiam << std::endl;
    // std::cout << "  fdomainRadiusInUm: " << fdomainRadiusInUm << std::endl;
    // std::cout << "  fNucleusRadius: " << fNucleusRadius << std::endl;
    // std::cout << "  fBetaRef: " << fBetaRef << std::endl;
    // std::cout << "  iunit: " << iunit << std::endl;
    // std::cout << "  mparased: " << mparased << std::endl;
    // calculator->reinitialize(nybin, CelDiam, fdomainRadiusInUm, 
    //                 fNucleusRadius, fBetaRef, iunit, mparased);

    // calculator->setCalculationFlags(flinealEnergySpectra, fmeanLinealEnergy, fdoseAveragedLinealEnergy);
    // std::cout << "InitializeUserInfo Calculation flags set - linealEnergySpectra: " << flinealEnergySpectra
    //               << ", meanLinealEnergy: " << fmeanLinealEnergy
    //               << ", doseAveragedLinealEnergy: " << fdoseAveragedLinealEnergy << std::endl;
}

void GateAMFActor::InitializeCpp() {
  GateVActor::InitializeCpp();
  NbOfThreads = G4Threading::GetNumberOfRunningWorkerThreads();

  calculator->setCalculationFlags(flinealEnergySpectra, fmeanLinealEnergy, fdoseAveragedLinealEnergy);


    // std::cout << "BeginOfRunActionMasterThread Calculation flags set - linealEnergySpectra: " << flinealEnergySpectra
    //               << ", meanLinealEnergy: " << fmeanLinealEnergy
    //               << ", doseAveragedLinealEnergy: " << fdoseAveragedLinealEnergy << std::endl;

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
//  std::cout << "End of InitializeCpp" << std::endl;
}

G4double GateAMFActor::getDose(G4Step *step) {
//in Joule/kg = Gy
//  std::cout << "Calculating dose..." << std::endl;
    // joule and kg/mm3 are defined in G4SystemOfUnits.hh
  // get edep in MeV (take weight into account)
  auto w = step->GetTrack()->GetWeight();
//   auto edep = step->GetTotalEnergyDeposit() / CLHEP::MeV * w;
  auto edep = step->GetTotalEnergyDeposit() / joule * w;
//   std::cout << "edep default unit: " << step->GetTotalEnergyDeposit() << " default" << std::endl;  
//   std::cout << "edep joule unit: " << step->GetTotalEnergyDeposit()/joule << " J" << std::endl;  
//   std::cout << "edep joule MeV: " << step->GetTotalEnergyDeposit()/CLHEP::MeV << " J" << std::endl;  
  double dose;
  double density;

  auto *current_material = step->GetPreStepPoint()->GetMaterial();
  density = current_material->GetDensity()/(kg/mm3); // ensure density is in kg/mm3, by default in G4 it is in internal units

//   std::cout << "density default unit: " << current_material->GetDensity() << " default" << std::endl;  
//   std::cout << "density kg/mm3 unit: " << current_material->GetDensity()/(kg/mm3) << " (kg/mm3)" << std::endl;  
  //   std::cout<< "Material name: " << current_material->GetName() << std::endl;
//   std::cout<< "Material formula: " << current_material->GetChemicalFormula() << std::endl;      
  dose = edep / (density*fVoxelVolume); // in Gy (J/kg)
//   std::cout << "dose Gy unit: " << dose << " Gy" << std::endl;
//   std::cout << "fVoxelVolume: " << fVoxelVolume << " (mm3)" << std::endl;
//   std::cout << "dose int units" << step->GetTotalEnergyDeposit()/(current_material->GetDensity()*fVoxelVolume)<<std::endl;
//   dose=step->GetTotalEnergyDeposit()/(current_material->GetDensity()*fVoxelVolume);
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

  auto total_dEdx = emcalc.ComputeTotalDEDX(kinEnergy, particle_definition, mat) / (keV/um);              // keV/um
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

//   std::cout << "AMF actor starting run BeginOfRunActionMasterThread"
//   << std::endl;
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

        // Create the image pointers
      // Important ! The volume may have moved, so we re-attach each run
        AttachImageToVolume<Image3DType>(cpp_amf_dose_image, fPhysicalVolumeName,
                                        fTranslation);


  auto sp = cpp_amf_dose_image->GetSpacing();
  fVoxelVolume = sp[0] * sp[1] * sp[2];
//   std::cout << "Voxel spacing: " << sp << " mm" << std::endl;
//   std::cout << "Voxel volume: " << fVoxelVolume << " mm3" << std::endl;
//   std::cout << "end of BeginOfRunActionMasterThread"
//   << std::endl;
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
        // G4double density = step->GetPreStepPoint()->GetMaterial()->GetDensity();
        G4double iAA = step->GetTrack()->GetDefinition()->GetAtomicMass();
        // G4double kenergy = step->GetPreStepPoint()->GetKineticEnergy();

        // Calculate the step's mean kinetic energy
        G4double eKinPre = step->GetPreStepPoint()->GetKineticEnergy() / (MeV); 

        G4double eKinPost = step->GetPostStepPoint()->GetKineticEnergy() / (MeV);
        G4double eKinMean = (eKinPre + eKinPost) * 0.5; // in MeV

        G4double energyPerNucleon = eKinMean / iAA; // energy per nucleon
        G4double izz = step->GetTrack()->GetDefinition()->GetAtomicNumber();



        // izz=6;
        // iAA=12;
        // energyPerNucleon=100.0;
        // dEdx=10.0;
        // dose=1.0;

        if (izz >= 1 && izz <= 18 && energyPerNucleon >= 0.025) {

            dEdx = GetStoppingPower(step) ; // in keV/um
            LinealEnergy_Dose = 0.0;
            LinealEnergyS = 0.0;
            // std::cout <<"Energy per nucleon: " << energyPerNucleon << " MeV/u" << ", Z: " << izz << ", A: " << iAA <<" dEdx: " << dEdx <<" dose: " << dose << std::endl;

            
            if (flinealEnergySpectra || fdoseAveragedLinealEnergy || fmeanLinealEnergy){ 
                // std::cout << "Calculating microdosimetric spectra for Z=" << izz << ", A=" << iAA << ", E/A=" << energyPerNucleon << " MeV/u, dE/dx=" << dEdx << " keV/um" << std::endl;
                // std::cout << "Dose: " << dose << " Gy" << std::endl;
                calculator->calculateDoseWeightedMicrodosimetricFunction(microdosimetricSpectra, izz, iAA, energyPerNucleon, dEdx, dose, LinealEnergy_Dose, LinealEnergyS);
                                    // Debug: print microdosimetric spectra content
                // for (size_t i = 0; i < microdosimetricSpectra.Size(); ++i) {
                //     std::cout << "microdosimetricSpectra[" << i << "] = " << microdosimetricSpectra[i] << std::endl;
                // }
                // Write microdosimetricSpectra to text file
                // calculator->get_Histo_X_Labels(histo_x_labels);
                // {
                //     std::ofstream spectraFile("microdosimetric_spectra_dump.txt", std::ios::app);
                //     if (spectraFile.is_open()) {
                //         for (size_t i = 0; i < microdosimetricSpectra.Size(); ++i) {
                //             spectraFile << histo_x_labels[i] << "\t" << microdosimetricSpectra[i] << std::endl;
                //         }
                //         spectraFile.close();
                //     }
                // }
            }




            
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
            // std::cout<<"LinealEnergyS: " << LinealEnergyS << "Dose averaged lineal energy: " << LinealEnergy_Dose << std::endl;
            
            {
                G4AutoLock mutex(&AMFMutex);
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

            } // end of G4AutoLock


            // std::cout<<"Voxel index: " << index << std::endl;

            // auto pixelValue = cpp_amf_microdosimetric_spectra->GetPixel(index);
            // std::cout << "Pixel value: " << pixelValue << std::endl;
            // std::cout << "End of SteppingAction" << std::endl;


            // writeVectorImage(cpp_amf_microdosimetric_spectra, fSpectraOutputFileName);
            // exit(0);

            return ;
        }
    }
        //  std::cout << "End of SteppingAction" << std::endl;
    
    return ;
}


int GateAMFActor::EndOfRunActionMasterThread(int run_id)
{
    //   std::cout << "begin of EndOfRunActionMasterThread"
    //         << std::endl;
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

    // std::cout << "begin of EndOfRunAction" << std::endl;
    {
    G4AutoLock mutex(&AMFMutex);
    {
        if (flinealEnergySpectra){
            if (!fRanOnce) {

                //divide cpp_amf_microdosimetric_spectra by dose
                divideVectorImageByScalarImage(cpp_amf_microdosimetric_spectra, cpp_amf_dose_image);

                writeVectorImage(cpp_amf_microdosimetric_spectra, fSpectraOutputFileName);

                // Remove file extension from fSpectraOutputFileName
                namespace fs = std::filesystem;
                fs::path p(fSpectraOutputFileName);
                calculator->get_Histo_X_Labels(histo_x_labels);
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
    // std::cout << "end of EndOfRunAction" << std::endl;
}

void GateAMFActor::divideVectorImageByScalarImage(const ImageVectorType::Pointer vectorImage,
                                        const Image3DType::Pointer scalarImage)
{
    itk::ImageRegionIterator<ImageVectorType> vecIt(vectorImage, vectorImage->GetRequestedRegion());
    itk::ImageRegionIterator<Image3DType> scalIt(scalarImage, scalarImage->GetRequestedRegion());

    for (vecIt.GoToBegin(), scalIt.GoToBegin(); !vecIt.IsAtEnd() && !scalIt.IsAtEnd(); ++vecIt, ++scalIt)
    {
        ImageVectorType::PixelType vecPixel = vecIt.Get();
        double scalPixel = scalIt.Get();

        if (scalPixel != 0)
        {
            for (unsigned int i = 0; i < vecPixel.GetSize(); ++i)
            {
                vecPixel[i] /= scalPixel;
            }
            vecIt.Set(vecPixel);
        }
        else
        {
            // Handle division by zero if necessary
            // For example, set the vector pixel to zero
            vecPixel.Fill(0.0);
            vecIt.Set(vecPixel);
        }
    }
}

void GateAMFActor::writeVectorImage(const ImageVectorType::Pointer image,
                                 std::string filename)
{
  using WriterType = itk::ImageFileWriter<ImageVectorType>;
  auto writer = WriterType::New();
//   std::cout << "Writing Vector image to " << filename << std::endl;
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
// // Public reinitialization function if parameters need to be updated
// void MicrodosimetricCalculator::reinitialize(size_t nybin_val, double celDiam, double domainRadius, 
//                     double nucleusRadius, double betaRef, int iunit_val, int mparased_val) {
//     nybin = nybin_val;
//     CelDiam = celDiam;
//     fDomainRadius = domainRadius;
//     fNucleusRadius = nucleusRadius;
//     fBetaRef = betaRef;
//     iunit = iunit_val;
//     mparased = mparased_val;
//     initialize();
// }

void MicrodosimetricCalculator::initialize() {
    ypower = -3.0;  // Initialize here instead

    if (yhig.size() != nybin + 1) yhig.resize(nybin + 1);
    if (yfy.size() != nybin) yfy.resize(nybin);
    if (ydy.size() != nybin) ydy.resize(nybin); // Resize ydy

    if (ymid.size() != nybin) ymid.resize(nybin);
    if (ywid.size() != nybin) ywid.resize(nybin);
    if (eventmid.size() != nybin) eventmid.resize(nybin);
    if (Z.size() != nybin) Z.resize(nybin);
    if (histo_x_labels.size() != nybin) histo_x_labels.resize(nybin);


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

    // // Calculate y0
    y0 = (M_PI * fDomainRadius * std::pow(fNucleusRadius, 2)) / (std::sqrt(fBetaRef * (std::pow(fDomainRadius, 2) + std::pow(fNucleusRadius, 2))) * 0.16022);

    for (size_t i = 0; i < nybin; ++i) {
        double ymid_val = (yhig[i] + yhig[i + 1]) / 2.0;
        ymid[i] = ymid_val;
        ywid[i] = yhig[i + 1] - yhig[i];
        // eventmid[i] = ymid_val * factor * unitconv;
        Z[i] = 1 - std::exp(-std::pow(ymid_val, 2) / std::pow(y0, 2));
        histo_x_labels[i] = ymid_val;
    // std::cout << "in initialize eventmid[" << i << "]: " << eventmid[i] << std::endl;
	// std::cout << "unit conv: " << unitconv << std::endl;
	// std::cout << "factor: " << factor << std::endl;
	// std::cout << "ymid[" << i << "]: " << ymid[i] << std::endl;
    }

    // // Calculate the bins per decade
    binsperDecade = nybin / (std::log10(yhig.back() / yhig[0]));
}


    // Constructor with all parameters
    MicrodosimetricCalculator::MicrodosimetricCalculator(size_t nybin_val, double celDiam, double domainRadius, 
                             double nucleusRadius, double betaRef, int iunit_val, int mparased_val)
        : nybin(nybin_val), CelDiam(celDiam), fDomainRadius(domainRadius),
          fNucleusRadius(nucleusRadius), fBetaRef(betaRef), iunit(iunit_val), mparased(mparased_val),
          factor(0.0), unitconv(0.0), binsperDecade(0.0), y0(0.0), ypower(-3.0) {
        initialize();
        std::cout << "MicrodosimetricCalculator Constructor - Parameters:" << std::endl;
        std::cout << "  nybin: " << nybin << std::endl;
        std::cout << "  CelDiam: " << CelDiam << std::endl;
        std::cout << "  fDomainRadius: " << fDomainRadius << std::endl;
        std::cout << "  fNucleusRadius: " << fNucleusRadius << std::endl;
        std::cout << "  fBetaRef: " << fBetaRef << std::endl;
        std::cout << "  iunit: " << iunit << std::endl;
        std::cout << "  mparased: " << mparased << std::endl;
    }


    
void MicrodosimetricCalculator::get_Histo_X_Labels(std::vector<double>& labels) const {
    labels = histo_x_labels;
}


void MicrodosimetricCalculator::calculateDoseWeightedMicrodosimetricFunction(VectorPixelType& microDosSpectra, double izz, double iAA, double energyPerNucleon, double dEdx, double dose, double& LinealEnergy_Dose, double& LinealEnergyS)
{
    if (microDosSpectra.Size() != nybin) {
        microDosSpectra.SetSize(nybin);
        microDosSpectra.Fill(0.0);
    }

    // std::cout << "calculateDoseWeightedMicrodosimetricFunction inputs:" << std::endl;
    // std::cout << "  izz: " << izz << std::endl;
    // std::cout << "  iAA: " << iAA << std::endl;
    // std::cout << "  energyPerNucleon: " << energyPerNucleon << std::endl;
    // std::cout << "  dEdx: " << dEdx << std::endl;
    // std::cout << "  dose: " << dose << std::endl;
    // std::cout << "  microDosSpectra.Size(): " << microDosSpectra.Size() << std::endl;
    // std::cout << "  nybin: " << nybin << std::endl;
    // std::cout << "  CelDiam: " << CelDiam << std::endl;
    // std::cout << "  fDomainRadius: " << DomainRadius << std::endl;
    // std::cout << "  iunit: " << iunit << std::endl;
    // std::cout << "  mparased: " << mparased << std::endl;

    // microDosSpectra.SetSize(nybin);
    // microDosSpectra.Fill(0.0);

    // double unitconv;
    double factor;
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    double Apara[mparased] = {0.0};

    // if (yhig.size() != nybin + 1) yhig.resize(nybin + 1);
    // if (yfy.size() != nybin) yfy.resize(nybin);
    // if (ydy.size() != nybin) ydy.resize(nybin); // Resize ydy

    // in initialize
    // double ypower = -3.0;
    // const double ystep = 0.02;	
    // for (size_t i = 0; i < yhig.size(); ++i) {
    // 	yhig[i] = std::pow(10.0, ypower);
    // 	ypower += ystep;
    // }
  
    // if (iunit <= 1) {
    //     unitconv = 1.0;
    // } else if (iunit == 2) {
    //     unitconv = 1.0e-3 * (2.0 / 3.0 * CelDiam);
    // } else if (iunit == 3) {
    //     unitconv = 4.0 / 3.0 * M_PI * std::pow(CelDiam / 2.0, 3) * 1.0e-15 / 1.602e-13;
    // }

    int ic1, ie1, ip1;
    double ratioc, ratioe, ratiop;
    double erg = energyPerNucleon * iAA;
    // std::cout << "dEdx: " << dEdx << " CelDiam: " << CelDiam<<" erg: " << erg << std::endl;
    double depev = std::min(dEdx * CelDiam * 1.0e3, erg * 1.0e6);
    // std::cout << "depev: " << depev << std::endl;

    getAparaion(CelDiam, energyPerNucleon, iAA, izz, ratioc, ratioe, ratiop, ic1, ie1, ip1);
    // std::cout << "getAparaion inputs - CelDiam: " << CelDiam << ", energyPerNucleon: " << energyPerNucleon << ", iAA: " << iAA << ", izz: " << izz << std::endl;
    
    sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
    // std::cout << "sedmean inputs - x: 1.0, depev: " << depev << ", ic1: " << ic1 << ", ie1: " << ie1 << ", ip1: " << ip1 << ", ratioc: " << ratioc << ", ratioe: " << ratioe << ", ratiop: " << ratiop << std::endl;
    factor = (iunit == 0) ? 1.0 : 1.0e6 / Apara[8];
//    double tmp = sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop);
//    std::cout << "factor: " << factor << std::endl;
//    std::cout << "Apara[8]: " << Apara[8] << std::endl;

    double sumYdy = 0.0;
    double zNumerator = 0.0;
    double zDenominator = 0.0;

    for (size_t i = 0; i < nybin; ++i) {
//         double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
// //        std::cout << "ymid_bin: " << ymid << std::endl;
//         double ywid = yhig[i + 1] - yhig[i];
        double eventmid = ymid[i] * factor * unitconv;
	// std::cout << "in calc eventmid[" << i << "]: "<< eventmid[i] << std::endl;
	// std::cout << "unit conv: " << unitconv << std::endl;
	// std::cout << "factor: " << factor << std::endl;
	// std::cout << "ymid[" << i << "]: " << ymid[i] << std::endl;
    // std::cout << "yfy calculation inputs - depev: " << depev << ", ic1: " << ic1 << ", ie1: " << ie1 << ", ip1: " << ip1 << ", ratioc: " << ratioc << ", ratioe: " << ratioe << ", ratiop: " << ratiop << std::endl;
        yfy[i] = ymid[i] * sedmean(eventmid, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
//        std::cout << "yfy_bin: " << ydy[i] << std::endl;
        ydy[i] = yfy[i] * ymid[i]; // Calculate ydy
        sum0 += yfy[i] * ywid[i] / ymid[i];
        sum1 += yfy[i] * ywid[i];
        sum2 += yfy[i] * ywid[i] * ymid[i];
        sumYdy += ydy[i];   // for normalization
        if (fdoseAveragedLinealEnergy) {
            zNumerator   += yfy[i] * Z[i];
            zDenominator += yfy[i];
        }

        // std::cout << "ymid: "<<ymid<<" ydy[" << i << "] before normalization: " << ydy[i] << std::endl;
    }

    // in initialize
    // Calculate the bins per decade and normalization factor
    // double binsperDecade = nybin / (std::log10(yhig.back() / yhig[0]));
    // std::cout << "Bins per decade: " << binsperDecade << std::endl;

    double normalization_factor = (binsperDecade / std::log(10)) / sumYdy;

    // double normalization_factor = (binsperDecade / std::log(10)) / std::accumulate(ydy.begin(), ydy.end(), 0.0);
    // std::cout << "Normalization factor: " << normalization_factor << std::endl;
    //    std::cout << "yF: " << sum1/sum0 << std::endl;
//    std::cout << "yD: " << sum2/sum1 << std::endl;


    for (size_t i = 0; i < ydy.size(); ++i) {
        ydy[i] *= normalization_factor;
        microDosSpectra[i] = ydy[i]*dose;
    }

    // yD calculation
    if (fmeanLinealEnergy){
        if (sum1 == 0) {
            std::cout << "Warning: sum1 is zero, returning zero vector." << std::endl;
            LinealEnergy_Dose = 0.0;
            return;
        }
        LinealEnergy_Dose = sum2/sum1;
        LinealEnergy_Dose *= dose;
    }
    // yS calculation
    if (fdoseAveragedLinealEnergy && zDenominator > 0.0) {
        double LinealEnergy_Freq = sum1 / sum0;
        LinealEnergyS = ((zNumerator / zDenominator) / LinealEnergy_Freq) * (y0 * y0);
        LinealEnergyS *= dose;
    }

    return;
}



void MicrodosimetricCalculator::getAparaion(const double& CelDiam, const double& energyPerNucleon, const int& iAA, const int& izz, double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1) {
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
        if (eincion[ie-1] >= energyPerNucleon) {
            break;
        }
    }

    if (ie == 1) {
        ie1 = 1;
        ratioe = 0.0;
    } else {
        ie1 = ie - 1;
        ratioe = std::min(1.0, (std::log10(energyPerNucleon) - std::log10(eincion[ie1-1])) / (std::log10(eincion[ie1]) - std::log10(eincion[ie1-1])));
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


inline double MicrodosimetricCalculator::sedfunc(double x, double depev, const double Apara[], size_t size) {
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

inline double MicrodosimetricCalculator::sedmean(double x, double depev, int ic1, int ie1, int ip1, double ratioc, double ratioe, double ratiop, double Apara[]) {
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
 //                   std::cout << "ip1: "<< ip1 <<", ic1: " << ic1 <<", ie1: " << ie1 << std::endl;
                }
                double wei = Rp * Re * Rc;
                double sedfuncResult = sedfunc(x, depev, Apara, mparased);
                sedmean += sedfuncResult * wei;
                A9 += Apara[8] * wei;
            }
        } 
    }
    Apara[8] = A9;
//    std::cout << "Apara0: " <<Apara[0] << ", Apara8: " << Apara[8] << std::endl;
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
