/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "GateAMDMActor.h"

void init_GateAMDMActor(py::module &m) {
  py::class_<GateAMDMActor, std::unique_ptr<GateAMDMActor, py::nodelete>,
             GateVActor>(m, "GateAMDMActor")
      .def(py::init<py::dict &>())
      //  .def_readwrite("cpp_amdm_delta_image",
      //  &GateAMDMActor::cpp_amdm_delta_image)
      //  .def_readwrite("cpp_amdm_gamma_image",
      //  &GateAMDMActor::cpp_amdm_gamma_image)
      .def_readwrite("cpp_amdm_restricted_edep_image",
                     &GateAMDMActor::cpp_amdm_restricted_edep_image)
      .def_readwrite("fPhysicalVolumeName",
                     &GateAMDMActor::fPhysicalVolumeName);
}
