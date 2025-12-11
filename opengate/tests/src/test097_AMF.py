#!/usr/bin/env python3

# #################################################################
# This software is distributed under the terms
# of the GNU Lesser General Public Licence (LGPL)
# Hermann Fuchs, Medical University of Vienna
# 10.08.2024
# #################################################################

import opengate as gate

import pathlib
import itk
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.optimize import curve_fit
import logging
# import argparse
import os

from opengate.tests import utility

import test097_AMF_helpers as helpers



def simulation(
    data_path,
    output_path,
    number_of_particles=1,
    init_particle_energy=120.0,
    number_of_threads=1,
    file_prefix="sim",
):
    """Run Gate simulation"""
    # units
    km = gate.g4_units.km
    m = gate.g4_units.m
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
    eV = gate.g4_units.eV
    MeV = gate.g4_units.MeV
    um = gate.g4_units.um
    nm = gate.g4_units.nm
    km = gate.g4_units.km
    gcm3 = gate.g4_units.g / gate.g4_units.cm3
    deg = gate.g4_units.deg
    mrad = gate.g4_units.mrad

    number_of_particles = number_of_particles / number_of_threads

    # simulation object
    sim = gate.Simulation()
    sim.running_verbose_level = 0
    sim.g4_verbose = False
    sim.visu = False
    sim.visu_verbose = False
    sim.random_engine = "MersenneTwister"
    sim.random_seed = 123456789
    sim.number_of_threads = number_of_threads
    # Materials
    sim.volume_manager.add_material_database(data_path / "../GateMaterials.db")

    # ######################################################################
    # # Defining geometry world plus phantom
    # ######################################################################
    world = sim.world
    world.size = [2 * m, 2 * m, 2 * m]
    world.material = "G4_AIR"
    # phantom cube
    phantom = sim.add_volume("Box", "phantom")
    phantom.material = "G4_WATER"
    phantom.mother = "world"  # by default
    phantom.size = [3.8 * cm, 3 * cm, 3 * cm]
    phantom.translation = [-phantom.size[0] / 2, 0 * mm, 0 * mm]
    rotation_matrix = (Rotation.from_euler("y", 180, degrees=True).as_matrix(),)
    phantom.rotation = rotation_matrix
    
    # ######################################################################
    # # Defining actors
    # ######################################################################

    # base_name for output files
    base_name = file_prefix

    # calculate actor size and spacing based on phantom size and resolution
    target_voxel_size = [1 * mm, phantom.size[1] * mm, phantom.size[2] * mm]
    is_dimensions = [
        int(phantom.size[0] / target_voxel_size[0]),
        int(phantom.size[1] / target_voxel_size[1]),
        int(phantom.size[2] / target_voxel_size[2]),
    ]

    amf_actor = sim.add_actor("AMFActor", "amf_actor")
    amf_actor.attached_to = "phantom"  
    amf_actor.spacing = target_voxel_size
    amf_actor.size = is_dimensions
    amf_actor.tsed_file_name = str(data_path / "tsed.dat")
    amf_actor.microdosimetric_spectra_file_name = os.path.join(output_path, base_name + "AMF_microdosimetricSpectra" + ".mhd")
    amf_actor.hit_type = "middle"
    amf_actor.DoseAveragedLinealEnergySaturationCorrected.active = True
    amf_actor.DoseAveragedLinealEnergy.active = True
    amf_actor.MicrodosimetricSpectra = True
    amf_actor.DomainRadius = 0.3 * um
    amf_actor.output_filename = os.path.join(output_path, base_name + "AMF" + ".mhd")


    ###############################
    # source
    ###############################
    source = sim.add_source("GenericSource", "Default")
    source.particle = 'ion 6 12'  # Carbon
    source.energy.mono = init_particle_energy*12* MeV
    source.position.radius = 1 * um
    source.direction.type = "momentum"
    source.direction.momentum = [-1, 0, 0]
    source.n = number_of_particles
    # ######################################################################
    # # Defining physics
    # ######################################################################
    sim.physics_manager.physics_list_name = "ShieldingLIQMD_HP_EMZ"
    global_cut = 1000000 * km
    sim.physics_manager.global_production_cuts.gamma = global_cut
    sim.physics_manager.global_production_cuts.electron = global_cut
    sim.physics_manager.global_production_cuts.positron = global_cut
    sim.physics_manager.global_production_cuts.proton = global_cut

    # adding cuts for waterphantom
    reg = sim.physics_manager.add_region("reg")
    reg.max_step_size = 50*um
    reg.production_cuts.gamma = 50*um
    reg.production_cuts.electron = 1000 * m
    reg.production_cuts.positron = 50*um
    reg.production_cuts.proton = 50*um
    reg.associate_volume("phantom")
    sim.physics_manager.set_user_limits_particles("all")

    #start simulation
    output_sim = sim.run(start_new_process=True)

    return output_sim



def main():

    paths = utility.get_default_test_paths(__file__, gate_folder="test_0097_AMF", output_folder="test0097")


    # logging.basicConfig(level=logging.DEBUG)
    logging.basicConfig(level=logging.INFO)
    print("Running Gate Simulation")
    data_path = paths.data / "test097"
    output_path = paths.output

    simulation(
        data_path,
        output_path,
        number_of_particles=50,
        init_particle_energy=120.0,
        number_of_threads=1,
        file_prefix="sim_",
    )
    print("Simulation finished here")

    all_results =[]

# Now compare the generated results
    ref_path = paths.output_ref / "ref_120_AMF_microdosimetricSpectra.mhd"
    test_path = paths.output / "sim_AMF_microdosimetricSpectra.mhd"

    rel_tols = {
        "integral": 0.20,   # 20% in area
        "com_x":    0.20,   # 20% in center of mass
        "fwhm":     0.20,   # 25% in FWHM
        "peak_x":   0.20,   # ~20% in peak position (in x-units)
        "q10":      0.20,   # 20% in lower tail level
        "q90":      0.20,   # 20% in upper tail level
        # "median":   0.20,
        "mean":     0.20,
        # "peak_y":   1.15,
    }

    mDimIsOK=helpers.testMultiDimContent(ref_path,test_path,rel_tols,calcVoxelPercent=0.1)
    all_results.append(mDimIsOK)

    rel_tols = {
        "integral": 0.20,   # 20% in area
        "com_x":    0.20,   # 20% in center of mass
        "fwhm":     0.20,   # 25% in FWHM
        "peak_x":   0.20,   # ~20% in peak position (in x-units)
        "q10":      0.20,   # 20% in lower tail level
        "q90":      0.20,   # 20% in upper tail level
        "median":   0.20,
        "mean":     0.20,
        # "peak_y":   1.15,
    }

    ref_files = [paths.output_ref / "ref_120_AMF_doseaveragedlinealenergy.mhd",
                 paths.output_ref / "ref_120_AMF_doseaveragedlinealenergysaturationcorrected.mhd",]
    test_files = [paths.output / "sim_AMF_doseaveragedlinealenergy.mhd",
                 paths.output / "sim_AMF_doseaveragedlinealenergysaturationcorrected.mhd",]

    for ref_file, test_file in zip(ref_files, test_files):
        oDimIsOK = helpers.testOneDimensionalContent(ref_file, test_file, rel_tols)
        all_results.append(oDimIsOK)

    print("===================================")
    print("all results: ", all_results )
    overall_results = all(all_results)
    if overall_results:
        print("Overall result: PASS")
    else:
        print("Overall result: FAIL")

    utility.test_ok(overall_results)
   

if __name__ == "__main__":
    main()
