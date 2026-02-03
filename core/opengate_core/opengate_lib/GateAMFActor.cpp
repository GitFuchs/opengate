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
  fTranslation = DictGetG4ThreeVector(user_info, "translation");
  fImageSpacing = DictGetG4ThreeVector(user_info, "spacing");
  fImageSize = DictGetG4ThreeVector(user_info, "size");
  fTSEDfilename = DictGetStr(user_info, "tsed_file_name");

    double CelDiam = 2.0 * fdomainRadiusInUm; // in um, fDomainRadiusInUm is in um
    // double nucleusRadius = 0.8 * fdomainRadiusInUm; // in um
    // double betaRef = 0.5 * fdomainRadiusInUm; // in um


    calculator = new MicrodosimetricCalculator(nybin,
                              CelDiam,
                              fdomainRadiusInUm,
                              fNucleusRadiusInUm,
                              fBetaRefinGyminus2, iunit, mparased);
    calculator->setTSEDfilename(fTSEDfilename);
}

void GateAMFActor::InitializeCpp() {
  GateVActor::InitializeCpp();
  NbOfThreads = G4Threading::GetNumberOfRunningWorkerThreads();

  calculator->setCalculationFlags(fMicrodosimetricSpectra, fdoseAveragedLinealEnergySaturationCorrected, fdoseAveragedLinealEnergy);

  if (fMicrodosimetricSpectra)
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

  if (fdoseAveragedLinealEnergySaturationCorrected)
    {
            // Create the image pointers
            cpp_amf_dose_averaged_lineal_energy_saturation_corrected = Image3DType::New();
    }
  if (fdoseAveragedLinealEnergy)
        {
            // Create the image pointers
            cpp_amf_dose_averaged_lineal_energy = Image3DType::New();
        }
  // Create the image pointers
  // (the size and allocation will be performed on the py side)
    cpp_amf_dose_image = Image3DType::New();
}

G4double GateAMFActor::getDose(G4Step *step) {
    // joule and kg/mm3 are defined in G4SystemOfUnits.hh
  // get edep in MeV (take weight into account)
  auto w = step->GetTrack()->GetWeight();
//   auto edep = step->GetTotalEnergyDeposit() / CLHEP::MeV * w;
  auto edep = step->GetTotalEnergyDeposit() / joule * w;
  double dose;
  double density;

  auto *current_material = step->GetPreStepPoint()->GetMaterial();
  density = current_material->GetDensity()/(kg/mm3); // ensure density is in kg/mm3, by default in G4 it is in internal units   
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

    if (fMicrodosimetricSpectra){
        AttachImageToVolume<ImageVectorType>(cpp_amf_microdosimetric_spectra, fPhysicalVolumeName,
                                    fTranslation);
    }
    if (fdoseAveragedLinealEnergySaturationCorrected){
    AttachImageToVolume<Image3DType>(cpp_amf_dose_averaged_lineal_energy_saturation_corrected, fPhysicalVolumeName,
                                    fTranslation);
    }
    if (fdoseAveragedLinealEnergy){
        AttachImageToVolume<Image3DType>(cpp_amf_dose_averaged_lineal_energy, fPhysicalVolumeName,
                                    fTranslation);
    }

    // Important ! The volume may have moved, so we re-attach each run
    AttachImageToVolume<Image3DType>(cpp_amf_dose_image, fPhysicalVolumeName,
                                        fTranslation);


  auto sp = cpp_amf_dose_image->GetSpacing();
  fVoxelVolume = sp[0] * sp[1] * sp[2];
}

void GateAMFActor::SteppingAction(G4Step *step) {
    //  std::cout << "Begin of SteppingAction" << std::endl;

  auto event_id =
      G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();

  double dEdx, LinealEnergy_Dose, LinealEnergy_Dose_saturation_correctedS;

  G4double dose = GateAMFActor::getDose(step);

  // Get the voxel index
  G4ThreeVector position;
  bool isInside;
  Image3DType::IndexType index;
  GetVoxelPosition(step, position, isInside, index);

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

        G4double energyPerNucleon = eKinMean / iAA; // energy per nucleon in MeV/u
        G4double izz = step->GetTrack()->GetDefinition()->GetAtomicNumber();

        if (izz >= 1 && izz <= 18 && energyPerNucleon >= 0.025) {

            dEdx = GetStoppingPower(step) ; // in keV/um
            LinealEnergy_Dose = 0.0;
            LinealEnergy_Dose_saturation_correctedS = 0.0;

            if (fMicrodosimetricSpectra || fdoseAveragedLinealEnergy || fdoseAveragedLinealEnergySaturationCorrected){ 
                // calculator->calculateDoseWeightedMicrodosimetricFunction(microdosimetricSpectra, izz, iAA, energyPerNucleon, dEdx, dose, LinealEnergy_Dose, LinealEnergy_Dose_saturation_correctedS);
            
                calculator->calculateDoseWeightedMicrodosimetricFunctionFast(microdosimetricSpectra, izz, iAA, energyPerNucleon, dEdx, dose, LinealEnergy_Dose, LinealEnergy_Dose_saturation_correctedS);

            }

            {
                G4AutoLock mutex(&AMFMutex);
                if (fMicrodosimetricSpectra){
                    ImageAddValue<ImageVectorType>(cpp_amf_microdosimetric_spectra, index, microdosimetricSpectra);
                }
                if (fdoseAveragedLinealEnergy){
                    ImageAddValue<Image3DType>(cpp_amf_dose_averaged_lineal_energy, index,LinealEnergy_Dose );

                }
                if (fdoseAveragedLinealEnergySaturationCorrected){
                    ImageAddValue<Image3DType>(cpp_amf_dose_averaged_lineal_energy_saturation_corrected, index, LinealEnergy_Dose_saturation_correctedS);
                }
                ImageAddValue<Image3DType>(cpp_amf_dose_image, index, dose);

            } // end of G4AutoLock

            return ;
        }
    }
        //  std::cout << "End of SteppingAction" << std::endl;
    
    return ;
}


int GateAMFActor::EndOfRunActionMasterThread(int run_id)
{
    return 0;}


  // Called every time a Run ends (all threads)
void GateAMFActor::EndOfRunAction(const G4Run *run)
{
    // std::cout << "begin of EndOfRunAction" << std::endl;
    {
        G4AutoLock mutex(&AMFMutex);
    
        if (fMicrodosimetricSpectra){
            if (!fRanOnce) {

                //divide cpp_amf_microdosimetric_spectra by dose
                divideVectorImageByScalarImage(cpp_amf_microdosimetric_spectra, cpp_amf_dose_image);

                //divide cpp_amf_dose_averaged_lineal_energy_saturation_corrected by cpp_amf_dose_image
                divideImage3DByImage3D(cpp_amf_dose_averaged_lineal_energy_saturation_corrected, cpp_amf_dose_image);

                //divide cpp_amf_dose_averaged_lineal_energy by cpp_amf_dose_image
                divideImage3DByImage3D(cpp_amf_dose_averaged_lineal_energy, cpp_amf_dose_image);

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
    // std::cout << "end of EndOfRunAction" << std::endl;
}

void GateAMFActor::divideImage3DByImage3D(Image3DType::Pointer numeratorImage,
                                          const Image3DType::Pointer denominatorImage)
{
    itk::ImageRegionIterator<Image3DType> numIt(numeratorImage, numeratorImage->GetRequestedRegion());
    itk::ImageRegionIterator<Image3DType> denomIt(denominatorImage, denominatorImage->GetRequestedRegion());

    for (numIt.GoToBegin(), denomIt.GoToBegin(); !numIt.IsAtEnd() && !denomIt.IsAtEnd(); ++numIt, ++denomIt)
    {
        double numPixel = numIt.Get();
        double denomPixel = denomIt.Get();

        if (denomPixel != 0.0)
        {
            numIt.Set(numPixel / denomPixel);
        }
        else
        {
            numIt.Set(0.0);
        }
    }
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




  
// #################################################################################
// #################################################################################
// #################################################################################
// Dedicated, optimized class for microdosimetric calculations
// #################################################################################
// #################################################################################
// #################################################################################

// Constructor with all parameters
MicrodosimetricCalculator::MicrodosimetricCalculator(size_t nybin_val, double celDiam, double domainRadius, 
                            double nucleusRadius, double betaRef, int iunit_val, int mparased_val)
    : nybin(nybin_val), CelDiam(celDiam), fDomainRadius(domainRadius),
        fNucleusRadius(nucleusRadius), fBetaRef(betaRef), iunit(iunit_val), mparased(mparased_val),
        factor(0.0), unitconv(0.0), binsperDecade(0.0), y0(0.0), ypower(-3.0) {
    initialize();
}


inline int MicrodosimetricCalculator::buildIonParamCombos(
    double depev,
    int ic1, int ie1, int ip1,
    double ratioc, double ratioe, double ratiop,
    IonParamCombo combos[8],          // output
    double& weightedA8                // output: Σ weight * A8
) const {
    int idx = 0;
    weightedA8 = 0.0;
    const bool depevZero = (depev == 0.0);

    for (int ip = ip1; ip <= ip1 + 1; ++ip) {
        double Rp = (ip == ip1) ? (1.0 - ratiop) : ratiop;
        for (int ie = ie1; ie <= ie1 + 1; ++ie) {
            double Re = (ie == ie1) ? (1.0 - ratioe) : ratioe;
            for (int ic = ic1; ic <= ic1 + 1; ++ic) {
                double Rc = (ic == ic1) ? (1.0 - ratioc) : ratioc;
                double w  = Rp * Re * Rc;

                int index = ((ip - 1) * 96) + ((ie - 1) * 8) + (ic - 1);
                const auto& row = IonData[index];

                IonParamCombo& c = combos[idx++];
                c.weight = w;
                c.A0 = row[0];
                c.A1 = row[1];
                c.A2 = row[2];
                c.A3 = row[3];
                c.A4 = row[4];
                c.A5 = row[5];
                c.A6 = row[6];
                c.A7 = row[7];
                c.A8 = row[8];

                weightedA8 += w * c.A8;

                if (!depevZero && c.A8 > 0.0 && c.A2 != 0.0) {
                    c.cst1    = depev / c.A8;
                    double denom = c.cst1 * c.A2;
                    c.firstPref = 2.0 / (denom * denom);
                } else {
                    c.cst1     = 0.0;
                    c.firstPref = 0.0;
                }

                if (c.A6 > 0.0 && c.A7 != 0.0) {
                    c.logBase7 = std::log((c.A7 - 1.0) / c.A7);
                } else {
                    c.logBase7 = 0.0;
                }
            }
        }
    }
    return idx; // should always be 8
}

inline double MicrodosimetricCalculator::sedmeanFast(
    double x,
    bool depevZero,
    const IonParamCombo* combos,
    int nCombos
) const {
    double sed = 0.0;

    for (int k = 0; k < nCombos; ++k) {
        const IonParamCombo& c = combos[k];
        const double w = c.weight;

        double getfirst = 0.0;
        double getsecond = 0.0;
        double getthird = 0.0;

        // First component
        if (c.A0 > 0.0) {
            if (depevZero) {
                double dx  = std::abs(x - c.A1);
                double tmp = (dx > 0.0)
                           ? std::pow(dx, c.A2) / (2.0 * c.A1)
                           : 0.0;
                if (tmp > 50.0) tmp = 50.0;
                getfirst = c.A0 * std::exp(-tmp);
            } else {
                double tmp = c.A1 * (x - c.cst1 * c.A2);
                if (tmp > 50.0) tmp = 50.0;
                double denom = std::exp(tmp) + 1.0;
                getfirst = c.A0 * x / denom * c.firstPref;
            }
        }

        // Second component
        if (c.A3 > 0.0) {
            double dx  = std::abs(x - c.A4);
            double tmp = (dx > 0.0)
                       ? std::pow(dx, c.A5) / (2.0 * c.A4)
                       : 0.0;
            if (tmp > 50.0) tmp = 50.0;
            getsecond = c.A3 * std::exp(-tmp);
        }

        // Third component
        if (c.A6 > 0.0) {
            // pow(base, x) -> exp(log(base) * x), log(base) precomputed
            double tmp = std::exp(c.logBase7 * x);
            getthird = c.A6 / (c.A7 - 1.0) * tmp;
        }

        double total = getfirst + getsecond + getthird;
        if (total > 1.0e-10) {
            sed += w * total;
        }
        // (else: original code clamps very small sedfunc to 0)
    }

    return sed;
}

void MicrodosimetricCalculator::calculateDoseWeightedMicrodosimetricFunctionFast(
    VectorPixelType& microDosSpectra,
    double izz,
    double iAA,
    double energyPerNucleon,
    double dEdx,
    double dose,
    double& LinealEnergy_Dose,
    double& LinealEnergy_Dose_saturation_correctedS)
{
    // Resize only if needed – no Fill(), we overwrite every element.
    if (microDosSpectra.Size() != nybin) {
        microDosSpectra.SetSize(nybin);
    }

    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    double sumYdy = 0.0;
    double zNumerator = 0.0;
    double zDenominator = 0.0;

    // Energy in eV
    const double erg = energyPerNucleon * iAA;
    const double depev = std::min(dEdx * CelDiam * 1.0e3, erg * 1.0e6);
    const bool depevZero = (depev == 0.0);

    int ic1, ie1, ip1;
    double ratioc, ratioe, ratiop;

    getAparaion(CelDiam, energyPerNucleon, static_cast<int>(iAA),
                static_cast<int>(izz),
                ratioc, ratioe, ratiop, ic1, ie1, ip1);

    // Precompute parameter combinations once
    IonParamCombo combos[8];
    double weightedA8 = 0.0;
    int nCombos = buildIonParamCombos(
        depev, ic1, ie1, ip1,
        ratioc, ratioe, ratiop,
        combos, weightedA8);

    // Factor using weighted A8 (matches original sedmean behaviour)
    double factorLocal;
    if (iunit == 0 || weightedA8 == 0.0) {
        factorLocal = 1.0;
    } else {
        factorLocal = 1.0e6 / weightedA8;
    }

    // main loop over y-bins
    for (size_t i = 0; i < nybin; ++i) {
        const double ymid_val = ymid[i];
        const double ywid_val = ywid[i];

        const double eventmid = ymid_val * factorLocal * unitconv;

        // y * sedmean(y)
        const double y_sed = ymid_val *
                             sedmeanFast(eventmid, depevZero, combos, nCombos);

        const double ydy_val = y_sed * ymid_val;   // y^2 * sedmean(y)

        // Accumulate integrals
        sum0   += y_sed * ywid_val / ymid_val;    // ∫ f(y) dy
        sum1   += y_sed * ywid_val;               // ∫ y f(y) dy
        sum2   += y_sed * ywid_val * ymid_val;    // ∫ y^2 f(y) dy
        sumYdy += ydy_val;                        // Σ ydy (for normalization)

        if (fdoseAveragedLinealEnergySaturationCorrected) {
            zNumerator   += y_sed * Z[i];
            zDenominator += y_sed;
        }

        // Store un-normalized ydy; we normalize in a second pass
        ydy[i] = ydy_val;
    }

    // Normalization factor (already precomputed binsperDecade in initialize())
    const double normalization_factor =
        (binsperDecade / std::log(10.0)) / sumYdy;

    for (size_t i = 0; i < nybin; ++i) {
        ydy[i] *= normalization_factor;
        microDosSpectra[i] = ydy[i] * dose;
    }
    // yD calculation
    if (fdoseAveragedLinealEnergy) {
        if (sum1 == 0.0) {
            std::cout << "Warning: sum1 is zero, returning zero vector." << std::endl;
            LinealEnergy_Dose = 0.0;
        } else {
            LinealEnergy_Dose = (sum2 / sum1) * dose;
        }
    }

    // yS calculation (saturation-corrected)
    if (fdoseAveragedLinealEnergySaturationCorrected && zDenominator > 0.0) {
        const double LinealEnergy_Freq = sum1 / sum0;
        LinealEnergy_Dose_saturation_correctedS =
            ((zNumerator / zDenominator) / LinealEnergy_Freq) * (y0 * y0) * dose;
    }
}


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
        Z[i] = 1 - std::exp(-std::pow(ymid_val, 2) / std::pow(y0, 2));
        histo_x_labels[i] = ymid_val;
    }

    // // Calculate the bins per decade
    binsperDecade = nybin / (std::log10(yhig.back() / yhig[0]));
}

  
void MicrodosimetricCalculator::get_Histo_X_Labels(std::vector<double>& labels) const {
    labels = histo_x_labels;
}


void MicrodosimetricCalculator::calculateDoseWeightedMicrodosimetricFunction(VectorPixelType& microDosSpectra, double izz, double iAA, double energyPerNucleon, double dEdx, double dose, double& LinealEnergy_Dose, double& LinealEnergy_Dose_saturation_correctedS)
{
    if (microDosSpectra.Size() != nybin) {
        microDosSpectra.SetSize(nybin);
        // microDosSpectra.Fill(0.0);
    }

    double factor;
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    double Apara[mparased] = {0.0};
    int ic1, ie1, ip1;
    double ratioc, ratioe, ratiop;
    double erg = energyPerNucleon * iAA;
    double depev = std::min(dEdx * CelDiam * 1.0e3, erg * 1.0e6);

    getAparaion(CelDiam, energyPerNucleon, iAA, izz, ratioc, ratioe, ratiop, ic1, ie1, ip1);    
    sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
    factor = (iunit == 0) ? 1.0 : 1.0e6 / Apara[8];

    double sumYdy = 0.0;
    double zNumerator = 0.0;
    double zDenominator = 0.0;

    for (size_t i = 0; i < nybin; ++i) {
        double eventmid = ymid[i] * factor * unitconv;
        yfy[i] = ymid[i] * sedmean(eventmid, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
        ydy[i] = yfy[i] * ymid[i]; // Calculate ydy
        sum0 += yfy[i] * ywid[i] / ymid[i];
        sum1 += yfy[i] * ywid[i];
        sum2 += yfy[i] * ywid[i] * ymid[i];
        sumYdy += ydy[i];   // for normalization
        if (fdoseAveragedLinealEnergySaturationCorrected) {
            zNumerator   += yfy[i] * Z[i];
            zDenominator += yfy[i];
        }
    }

    double normalization_factor = (binsperDecade / std::log(10)) / sumYdy;

    for (size_t i = 0; i < ydy.size(); ++i) {
        ydy[i] *= normalization_factor;
        microDosSpectra[i] = ydy[i]*dose;
    }

    // yD calculation
    if (fdoseAveragedLinealEnergy){
        if (sum1 == 0) {
            std::cout << "Warning: sum1 is zero, returning zero vector." << std::endl;
            LinealEnergy_Dose = 0.0;
            return;
        }
        LinealEnergy_Dose = sum2/sum1;
        LinealEnergy_Dose *= dose;
    }
    // yS calculation
    if (fdoseAveragedLinealEnergySaturationCorrected && zDenominator > 0.0) {
        double LinealEnergy_Freq = sum1 / sum0;
        LinealEnergy_Dose_saturation_correctedS = ((zNumerator / zDenominator) / LinealEnergy_Freq) * (y0 * y0);
        LinealEnergy_Dose_saturation_correctedS *= dose;
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
