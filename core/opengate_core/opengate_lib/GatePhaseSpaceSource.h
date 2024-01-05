/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GatePhaseSpaceSource_h
#define GatePhaseSpaceSource_h

#include "GateSPSVoxelsPosDistribution.h"
#include "GateSingleParticleSource.h"
#include "GateVSource.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class GatePhaseSpaceSource : public GateVSource {

public:
  // signature of the callback function in Python
  // that will generate particles (from phsp file)
  using ParticleGeneratorType = std::function<int(GatePhaseSpaceSource *)>;

  GatePhaseSpaceSource();

  ~GatePhaseSpaceSource() override;

  void InitializeUserInfo(py::dict &user_info) override;

  void PrepareNextRun() override;

  double PrepareNextTime(double current_simulation_time) override;

  void GeneratePrimaries(G4Event *event,
                         double current_simulation_time) override;

  void GenerateOnePrimary(G4Event *event, double current_simulation_time);

  void AddOnePrimaryVertex(G4Event *event, const G4ThreeVector &position,
                           const G4ThreeVector &momentum_direction,
                           double energy, double time, double w);

  void SetGeneratorFunction(ParticleGeneratorType &f) const;

  bool ParticleIsPrimary();

  // virtual void SetGeneratorInfo(py::dict &user_info);

  void GenerateBatchOfParticles();

  G4ParticleDefinition *fParticleDefinition;
  G4ParticleTable *fParticleTable;

  float fCharge;
  float fMass;
  bool fGlobalFag;
  bool fUseParticleTypeFromFile;
  bool fVerbose;

  // unsigned long fMaxN;
  long fNumberOfGeneratedEvents;
  size_t fCurrentBatchSize;

  void SetPDGCodeBatch(const py::array_t<std::int32_t> &fPDGCode) const;

  void SetEnergyBatch(const py::array_t<std::float_t> &fEnergy) const;

  void SetWeightBatch(const py::array_t<std::float_t> &fWeight) const;

  void SetPositionXBatch(const py::array_t<std::float_t> &fPositionX) const;

  void SetPositionYBatch(const py::array_t<std::float_t> &fPositionY) const;

  void SetPositionZBatch(const py::array_t<std::float_t> &fPositionZ) const;

  void SetDirectionXBatch(const py::array_t<std::float_t> &fDirectionX) const;

  void SetDirectionYBatch(const py::array_t<std::float_t> &fDirectionY) const;

  void SetDirectionZBatch(const py::array_t<std::float_t> &fDirectionZ) const;

  // For MT, all threads local variables are gathered here
  struct threadLocalTPhsp {

    bool fgenerate_until_next_primary;
    int fprimary_PDGCode;
    float fprimary_lower_energy_threshold;

    ParticleGeneratorType fGenerator;
    unsigned long fNumberOfGeneratedEvents;
    size_t fCurrentIndex;
    size_t fCurrentBatchSize;

    int *fPDGCode;

    float *fPositionX;
    float *fPositionY;
    float *fPositionZ;

    float *fDirectionX;
    float *fDirectionY;
    float *fDirectionZ;

    float *fEnergy;
    float *fWeight;
    // double * fTime;
  };
  G4Cache<threadLocalTPhsp> fThreadLocalDataPhsp;
};

#endif // GatePhaseSpaceSource_h
