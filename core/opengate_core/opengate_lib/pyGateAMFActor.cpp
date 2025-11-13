/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "GateAMFActor.h"

class PyGateAMFActor : public GateAMFActor {
public:
  // Inherit the constructors
  using GateAMFActor::GateAMFActor;

  void BeginOfRunActionMasterThread(int run_id) override {
    PYBIND11_OVERLOAD(void, GateAMFActor, BeginOfRunActionMasterThread, run_id);
  }

  int EndOfRunActionMasterThread(int run_id) override {
    PYBIND11_OVERLOAD(int, GateAMFActor, EndOfRunActionMasterThread, run_id);
  }
};

void init_GateAMFActor(py::module &m) {
  py::class_<GateAMFActor, PyGateAMFActor,
             std::unique_ptr<GateAMFActor, py::nodelete>,
             GateVActor>(m, "GateAMFActor")
      .def(py::init<py::dict &>())
      .def("BeginOfRunActionMasterThread",
           &GateAMFActor::BeginOfRunActionMasterThread)
      .def("EndOfRunActionMasterThread",
           &GateAMFActor::EndOfRunActionMasterThread)
      .def_readwrite("cpp_amf_dose_image", &GateAMFActor::cpp_amf_dose_image)
      .def_readwrite("cpp_amf_mean_lineal_energy",
                     &GateAMFActor::cpp_amf_mean_lineal_energy)
      .def_readwrite("cpp_amf_dose_averaged_lineal_energy",
                     &GateAMFActor::cpp_amf_dose_averaged_lineal_energy)
      // // .def_readwrite("NbOfEvent", &GateAMFActor::NbOfEvent)
      .def("GetPhysicalVolumeName", &GateAMFActor::GetPhysicalVolumeName)
      .def("SetPhysicalVolumeName", &GateAMFActor::SetPhysicalVolumeName)
      .def_readwrite("NbOfEvent", &GateAMFActor::NbOfEvent)
      .def("GetLinealEnergySpectraFlag", &GateAMFActor::GetLinealEnergySpectraFlag)
      .def("SetLinealEnergySpectraFlag", &GateAMFActor::SetLinealEnergySpectraFlag)
      .def("GetMeanLinealEnergyFlag", &GateAMFActor::GetMeanLinealEnergyFlag)
      .def("SetMeanLinealEnergyFlag", &GateAMFActor::SetMeanLinealEnergyFlag)
      .def("GetDoseAveragedLinealEnergyFlag", &GateAMFActor::GetDoseAveragedLinealEnergyFlag)
      .def("SetDoseAveragedLinealEnergyFlag", &GateAMFActor::SetDoseAveragedLinealEnergyFlag)
       .def_readwrite("fPhysicalVolumeName",
       &GateAMFActor::fPhysicalVolumeName);
}
