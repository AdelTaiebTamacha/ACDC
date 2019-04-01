#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import math
import numpy

from MARILib.earth import environment as earth

#===========================================================================================================
def eval_turboprop_pylon_mass(aircraft):
    """
    Turbofan pylon mass & CG estimation
    """

    nacelle = aircraft.turboprop_nacelle
    engine = aircraft.turboprop_engine

    pylon = aircraft.turboprop_pylon

    pylon.mass = 0.

    pylon.c_g = 0.

    return


#===========================================================================================================
def eval_turboprop_engine_design(aircraft):
    """
    Turboprop architecture design
    """

    engine = aircraft.turboprop_engine

    engine.rating_factor = (0.800,0.688,0.624,0.560,0.100)

    disa = 0.
    altp = 0.
    mach = 0.25

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    Vsnd = earth.sound_speed(tamb)
    Vair = Vsnd*mach

    engine.reference_power = engine.reference_thrust*(Vair/engine.propeller_efficiency)

    # assuming a fan disk load of 3000 N/m2
    engine.propeller_diameter = math.sqrt((4./math.pi)*(engine.reference_thrust/3000))

    return


#===========================================================================================================
def eval_turboprop_nacelle_design(aircraft):
    """
    Turboprop architecture design
    """

    fuselage = aircraft.fuselage
    wing = aircraft.wing
    engine = aircraft.turboprop_engine

    nacelle = aircraft.turboprop_nacelle

    nacelle.width = 0.25*(engine.reference_power/1.e3)**0.2        # statistical regression

    nacelle.length = 0.84*(engine.reference_power/1.e3)**0.2       # statistical regression

    nacelle.net_wetted_area = (2.3*(engine.reference_power/1.e3)**0.2)*engine.n_engine     # statistical regression

    nacelle.y_ext = 0.5 * fuselage.width + 0.8 * engine.propeller_diameter      # statistical regression

    nacelle.z_ext =  - 0.5 * fuselage.height \
                     + wing.z_root + (nacelle.y_ext-wing.y_root)*math.tan(wing.dihedral)

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + math.tan(wing.sweep)

    nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

    return


#===========================================================================================================
def eval_turboprop_nacelle_mass(aircraft):
    """
    Turboprop architecture design
    """

    engine = aircraft.turboprop_engine
    nacelle = aircraft.turboprop_nacelle

    nacelle.mass = (1.266*(engine.reference_power/1.e3)**0.9)*engine.n_engine       # statistical regression

    nacelle.c_g = nacelle.x_ext + 0.7 * nacelle.length      # statistical regression

    return


