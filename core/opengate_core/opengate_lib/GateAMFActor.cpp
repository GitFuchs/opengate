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
  loadIonData();

    // DomainRadius = fPm->GetDoubleParameter(GetFullParmName("DomainRadius"), "Length");
    // DomainRadius
    // G4cout << "DomainRadius: " << DomainRadius / um << " um" << G4endl;
    // DomainRadius = DomainRadius / um;   // Convert DomainRadius from mm to um   
    fDomainRadius = 0.3; // in um
    CelDiam = 2.0 * fDomainRadius;
    fNucleusRadius = 0.8 * fDomainRadius; // in um
    fBetaRef = 0.5 * fDomainRadius; // in um
        // Perform unit conversions 
    // DomainRadius = DomainRadius / um;   // Convert DomainRadius from mm to um 
    // NucleusRadius = NucleusRadius / um; // Convert NucleusRadius from mm to um 
    // BetaRef = BetaRef / (1. / (gray * gray)) ;             // Convert BetaRef to /Gy2
}

void GateAMFActor::InitializeCpp() {
  GateVActor::InitializeCpp();
  NbOfThreads = G4Threading::GetNumberOfRunningWorkerThreads();

  // Create the image pointers
  // (the size and allocation will be performed on the py side)

    cpp_amf_dose_image = Image3DType::New();
    cpp_amf_mean_lineal_energy = Image3DType::New();
    cpp_amf_dose_averaged_lineal_energy = Image3DType::New();
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
  cpp_amf_microdosimetric_spectra->SetNumberOfComponentsPerPixel(nybin*2);

  cpp_amf_microdosimetric_spectra->Allocate();

  // Optionally initialize pixels
    ImageVectorType::PixelType pixel;
    pixel.SetSize(nybin*2);
    pixel.Fill(0.0);
    cpp_amf_microdosimetric_spectra->FillBuffer(pixel);

//   itk::VariableLengthVector<double> v;
//     v.SetSize(nybin*2);
//     v.Fill(0.0);
//    cpp_amf_microdosimetric_spectra->FillBuffer(v);

//  std::cout << "AMF image size: " << size[0] << " " << size[1] << " " << size[2] << std::endl;
//  std::cout << "AMF image spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;
 std::cout << "End of InitializeCpp" << std::endl;
}


GateAMFActor::VectorPixelType GateAMFActor::calculateDoseWeightedMicrodosimetricFunction(double izz, double iAA, double ene, double dEdx, double dose, double& LinealEnergy_Dose, double& LinealEnergyS) {
    double unitconv, factor;
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    double Apara[mparased] = {0.0};

    // itk::VariableLengthVector<double> v;
    // ImageVectorType::PixelType v;
    VectorPixelType v;
    v.SetSize(nybin*2);
    v.Fill(0.0);


    yhig.resize(nybin + 1);
    yfy.resize(nybin);
    ydy.resize(nybin); // Resize ydy

    double ypower = -3.0;
    const double ystep = 0.02;
	
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

    int ic1, ie1, ip1;
    double ratioc, ratioe, ratiop;
    // std::vector<std::pair<double, double>> microdosimetricSpectra(nybin);

    double erg = ene * iAA;
    double depev = std::min(dEdx * CelDiam * 1.0e3, erg * 1.0e6);

    getAparaion(CelDiam, ene, iAA, izz, ratioc, ratioe, ratiop, ic1, ie1, ip1);
    
    // std::cout << "CelDiam: " << CelDiam << std::endl;
    // std::cout << "ene: " << ene << std::endl;
    // std::cout << "iAA: " << iAA << std::endl;
    // std::cout << "izz: " << izz << std::endl;
    // std::cout << "ratioc: " << ratioc << std::endl;
    // std::cout << "ratioe: " << ratioe << std::endl;
    // std::cout << "ratiop: " << ratiop << std::endl;
    // std::cout << "ic1: " << ic1 << std::endl;
    // std::cout << "ie1: " << ie1 << std::endl;
    // std::cout << "ip1: " << ip1 << std::endl;

    // std::cout << "Apara values before sedmean: ";
    // for (int i = 0; i < mparased; ++i) {
    //     std::cout << Apara[i] << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "IonData matrix (" << IonData.size() << " x " << (IonData.empty() ? 0 : IonData[0].size()) << "):" << std::endl;
    // for (size_t i = 0; i < IonData.size(); ++i) {
    //     std::cout << "Row " << i << ": ";
    //     for (size_t j = 0; j < IonData[i].size(); ++j) {
    //         std::cout << IonData[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    
    double temp=sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
    // std::cout << "Return value of sedmean: " << temp << std::endl;
    // std::cout << "Apara values after sedmean: ";
    // for (int i = 0; i < mparased; ++i) {
    //     std::cout << Apara[i] << " ";
    // }
    
    // std::cout << std::endl;
    // exit(1);
    factor = (iunit == 0) ? 1.0 : 1.0e6 / Apara[8];
    for (size_t i = 0; i < nybin; ++i) {
        double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
        double ywid = yhig[i + 1] - yhig[i];
        double eventmid = ymid * factor * unitconv;
        yfy[i] = ymid * sedmean(eventmid, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
        ydy[i] = yfy[i] * ymid; // Calculate ydy
        // std::cout << "ymid: " << ymid << std::endl;
        // std::cout << "yfy[i]: " << yfy[i] << std::endl;
        // std::cout << "ydy[i]: " << ydy[i] << std::endl;
        sum0 += yfy[i] * ywid / ymid;
        sum1 += yfy[i] * ywid;
        sum2 += yfy[i] * ywid * ymid;
    }


    // Calculate the bins per decade and normalization factor
    double binsperDecade = nybin / (std::log10(yhig.back() / yhig[0]));
    double normalization_factor = (binsperDecade / std::log(10)) / std::accumulate(ydy.begin(), ydy.end(), 0.0);

    // Apply normalization to ydy
    for (size_t i = 0; i < ydy.size(); ++i) {
        ydy[i] *= normalization_factor;
    }

    // for (size_t i = 0; i + 1 < yhig.size(); ++i) {
    //     double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
    //     microdosimetricSpectra[i] = std::make_pair(ymid, ydy[i]); // Use ydy instead of yfy
    // }
    for (size_t k = 0; k < nybin; ++k) {
        double ymid = (yhig[k] + yhig[k + 1]) / 2.0;
        // v[2*k]   = ymid;
        // v[2*k+1] = ydy[k]*dose; // Use ydy instead of yfy
        setBinValueAndContent(v, k, ymid, ydy[k]*dose);

    // v[2*k]   = x[k];
    // v[2*k+1] = y[k];
    }
    if (sum1 == 0)
    {
        std::cout << "Warning: sum1 is zero, returning zero vector." << std::endl;
        LinealEnergy_Dose = 0.0;
        LinealEnergyS = 0.0;
        return v;
    }
    LinealEnergy_Dose = sum2/sum1;
    // std::cout << "LinealEnergy_Dose: " << LinealEnergy_Dose << std::endl;

    double LinealEnergy_Freq = sum1 / sum0;
    if (true) {
        
        // std::cout << "LinealEnergy_Dose: " << LinealEnergy_Dose << std::endl;
        // yS calculation
        double y0 = (M_PI * fDomainRadius * std::pow(fNucleusRadius, 2)) / (std::sqrt(fBetaRef * (std::pow(fDomainRadius, 2) + std::pow(fNucleusRadius, 2))) * 0.16022);
        std::vector<double> Z(nybin);
        
        for (size_t i = 0; i < nybin; ++i) {
            double F2 = (yhig[i] + yhig[i + 1]) / 2.0; // Bin center
            Z[i] = 1 - std::exp(-std::pow(F2, 2) / std::pow(y0, 2));
        }

        double sumNumerator = 0.0;
        double sumDenominator = 0.0;

        for (size_t i = 0; i < nybin; ++i) {
            sumNumerator += yfy[i] * Z[i];
            sumDenominator += yfy[i];
        }
        double LinealEnergyS = 0.0;
        if (sumDenominator > 0.0) {
            LinealEnergyS = ((sumNumerator / sumDenominator) / LinealEnergy_Freq ) * std::pow(y0, 2);
        }
        // std::cout << "LinealEnergyS: " << LinealEnergyS << std::endl;
        }

    return v;
}

void GateAMFActor::getBinValueAndContent(const VectorPixelType& vec, 
                                          size_t index, 
                                          double& binValue, 
                                          double& binContent) const {
    if (index >= vec.Size() / 2) {
        std::cerr << "Error: index out of range" << std::endl;
        binValue = 0.0;
        binContent = 0.0;
        return;
    }
    binValue = vec[2 * index];
    binContent = vec[2 * index + 1];
}

void GateAMFActor::setBinValueAndContent(VectorPixelType& vec, 
                                          size_t index, 
                                          double binValue, 
                                          double binContent) {
    if (index >= vec.Size() / 2) {
        std::cerr << "Error: index out of range" << std::endl;
        return;
    }
    vec[2 * index] = binValue;
    vec[2 * index + 1] = binContent;
}


double GateAMFActor::sedfunc(double x, double depev, const double Apara[], size_t size) {
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

double GateAMFActor::sedmean(double x, double depev, int ic1, int ie1, int ip1, double ratioc, double ratioe, double ratiop, double Apara[]) {
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

void GateAMFActor::loadIonData() {
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



void GateAMFActor::getAparaion(const double& CelDiam, const double& ene, const int& iAA, const int& izz, double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1) {
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
  std::cout << "fPhysicalVolumeName: " << fPhysicalVolumeName << std::endl;
  std::cout << "fInitialTranslation: " << fTranslation << std::endl;  

      // Important ! The volume may have moved, so we re-attach each run
  AttachImageToVolume<Image3DType>(cpp_amf_dose_image, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<Image3DType>(cpp_amf_mean_lineal_energy, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<Image3DType>(cpp_amf_dose_averaged_lineal_energy, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<ImageVectorType>(cpp_amf_microdosimetric_spectra, fPhysicalVolumeName,
                                   fTranslation);

//  Image3DType::PointType currentOrigin = cpp_amf_dose_image->GetOrigin();
//      std::cout << "cpp_amf_dose_image Origin: "
//               << currentOrigin[0] << ", "
//               << currentOrigin[1] << ", "
//               << currentOrigin[2] << std::endl;

//  auto currentOriginTwo = cpp_amf_microdosimetric_spectra->GetOrigin();
//      std::cout << "cpp_amf_dose_image cpp_amf_microdosimetric_spectra: "
//               << currentOriginTwo[0] << ", "
//               << currentOriginTwo[1] << ", "
//               << currentOriginTwo[2] << std::endl;

  auto sp = cpp_amf_dose_image->GetSpacing();
  fVoxelVolume = sp[0] * sp[1] * sp[2];
//   std::cout << "Voxel spacing: " << sp << " mm" << std::endl;
//   std::cout << "Voxel volume: " << fVoxelVolume << " mm3" << std::endl;
  std::cout << "end of BeginOfRunActionMasterThread"
  << std::endl;
}

void GateAMFActor::SteppingAction(G4Step *step) {
     std::cout << "Begin of SteppingAction" << std::endl;

  auto event_id =
      G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();


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

            double dEdx = GetStoppingPower(step);
            double LinealEnergy_Dose;
            double LinealEnergyS;
            auto microdosimetricSpectra = calculateDoseWeightedMicrodosimetricFunction(izz, iAA, energyPerNucleon, dEdx, dose, LinealEnergy_Dose, LinealEnergyS);

            // Set every 2nd value (i.e., index 1, 3, 5, ...) to zero, e.g. the indexes corresponding to the contents
            if (fRanOnce) {
                // for (size_t k = 0; k < nybin; ++k) {
                // double binValue, binContent;
                // getBinValueAndContent(microdosimetricSpectra, k, binValue, binContent);
                // std::cout << "Bin " << k << ": Value = " << binValue << ", Content = " << binContent << std::endl;
                // }

                // std::cout << "Setting every 2nd value of microdosimetricSpectra to zero." << std::endl;
                for (unsigned int i = 0; i < microdosimetricSpectra.GetSize(); i += 2) {
                microdosimetricSpectra[i] = 0.0;
                }

            }
            fRanOnce = true;
            // for (size_t k = 0; k < nybin; ++k) {
            //     double binValue, binContent;
            //     getBinValueAndContent(microdosimetricSpectra, k, binValue, binContent);
            //     std::cout << "Bin " << k << ": Value = " << binValue << ", Content = " << binContent << std::endl;
            // }

            
            // LinealEnergyS=12.2;
            std::cout << "dEdx: " << dEdx << std::endl;
            std::cout << "LinealEnergy_Dose: " << LinealEnergy_Dose << std::endl;
            std::cout << "LinealEnergyS: " << LinealEnergyS << std::endl;
            std::cout << "dose: " << dose << std::endl;
            // std::cout << "Voxel index: " << index << std::endl;
  
            ImageAddValue<Image3DType>(cpp_amf_dose_image, index, dose);
            ImageAddValue<Image3DType>(cpp_amf_dose_averaged_lineal_energy, index, LinealEnergyS);
            ImageAddValue<Image3DType>(cpp_amf_mean_lineal_energy, index, LinealEnergy_Dose);

            ImageAddValue<ImageVectorType>(cpp_amf_microdosimetric_spectra, index, microdosimetricSpectra);

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
    writeVectorImage(cpp_amf_microdosimetric_spectra, fSpectraOutputFileName);
    Image3DType::IndexType index;
    index[0] = 0;
    index[1] = 0;
    index[2] = 0;
    auto pixelValue = cpp_amf_microdosimetric_spectra->GetPixel(index);
    std::cout << "Pixel value: " << pixelValue << std::endl;
    std::cout << "cpp_amf_microdosimetric_spectra: " <<cpp_amf_microdosimetric_spectra->GetNumberOfComponentsPerPixel() << std::endl;


    }
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

  // Prefer this helper (true for master in MT, true in sequential mode as well)
G4bool GateAMFActor::IsMaster() const {
    return G4Threading::IsMasterThread();
  }