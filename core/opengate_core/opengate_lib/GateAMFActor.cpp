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

#include <cmath>

GateAMFActor::~GateAMFActor() {
}


GateAMFActor::GateAMFActor(py::dict &user_info)
    : GateVActor(user_info) {}

void GateAMFActor::InitializeUserInfo(py::dict &user_info) {
  // IMPORTANT: call the base class method
  GateVActor::InitializeUserInfo(user_info);
  loadIonData();

}

std::vector<std::pair<double, double>> GateAMFActor::calculateMicrodosimetricFunction(double izz, double iAA, double ene, double dEdx) {
    double unitconv, factor;
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    double Apara[mparased] = {0.0};

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
    std::vector<std::pair<double, double>> microdosimetricSpectra(nybin);

    double erg = ene * iAA;
    double depev = std::min(dEdx * CelDiam * 1.0e3, erg * 1.0e6);

    getAparaion(CelDiam, ene, iAA, izz, ratioc, ratioe, ratiop, ic1, ie1, ip1);
    sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
    factor = (iunit == 0) ? 1.0 : 1.0e6 / Apara[8];

    for (size_t i = 0; i < nybin; ++i) {
        double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
        double ywid = yhig[i + 1] - yhig[i];
        double eventmid = ymid * factor * unitconv;
        yfy[i] = ymid * sedmean(eventmid, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
        ydy[i] = yfy[i] * ymid; // Calculate ydy
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

    for (size_t i = 0; i + 1 < yhig.size(); ++i) {
        double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
        microdosimetricSpectra[i] = std::make_pair(ymid, ydy[i]); // Use ydy instead of yfy
    }

    return microdosimetricSpectra;
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
    std::ifstream file("tsed.dat");
    if (!file.is_open()) {
        std::cerr << "Failed to open tsed.dat for reading." << std::endl;
        return;
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
  auto edep = step->GetTotalEnergyDeposit() / CLHEP::MeV * w;
  double dose;
  double density;

  auto *current_material = step->GetPreStepPoint()->GetMaterial();
  density = current_material->GetDensity();
      
  dose = edep / density;
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



void GateAMFActor::BeginOfRunAction(const G4Run *) {

}

void GateAMFActor::SteppingAction(G4Step *step) {

  G4double dose = GateAMFActor::getDose(step);

    if (dose > 0.) {
        G4double density = step->GetPreStepPoint()->GetMaterial()->GetDensity();
        G4double iAA = step->GetTrack()->GetDefinition()->GetAtomicMass();
        G4double kenergy = step->GetPreStepPoint()->GetKineticEnergy();
        G4double ene = kenergy / iAA; // energy per nucleon
        G4double energy = ene;
        G4double izz = step->GetTrack()->GetDefinition()->GetAtomicNumber();

        if (izz >= 1 && izz <= 18 && ene >= 0.025) {
            G4int binIndex = GetIndex(step);

            double dEdx = GetStoppingPower(step);

            auto microdosimetricSpectra = calculateMicrodosimetricFunction(izz, iAA, ene, dEdx);

//             if (std::isnan(microdosimetricSpectra[0].second)) {
//  /*               G4cout << "Debug Info - NAN Detected: "
//                        << "Atomic Number: " << izz
//                        << ", Energy: " << energy
//                        << ", dEdx: " << dEdx
//                        << G4endl; */
//             }

            if (totalSpectra.find(binIndex) == totalSpectra.end()) {
                totalSpectra[binIndex] = std::vector<G4double>(microdosimetricSpectra.size(), 0.0);
            } else if (totalSpectra[binIndex].size() < microdosimetricSpectra.size()) {
                totalSpectra[binIndex].resize(microdosimetricSpectra.size(), 0.0);
            }

           for (size_t i = 0; i < microdosimetricSpectra.size(); ++i) {
                totalSpectra[binIndex][i] += microdosimetricSpectra[i].second * dose;
           }

            cumulativeDose[binIndex] += dose;
            return ;
        }
    }
    return ;
}

void GateAMFActor::EndOfEventAction(const G4Event *event) {
  for (auto actor : fEndOfEventAction_actors) {
    actor->EndOfEventAction(event);
  }
}
