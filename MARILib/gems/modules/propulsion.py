#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import math
import numpy


#===========================================================================================================
def eval_turbofan_pylon_mass(aircraft):
    """
    Turbofan pylon mass & CG estimation
    """

    nacelle = aircraft.turbofan_nacelle
    engine = aircraft.turbofan_engine

    pylon = aircraft.turbofan_pylon

    pylon.mass = 0.0031*engine.reference_thrust*engine.n_engine

    pylon.c_g = nacelle.x_ext + nacelle.length

    return


#===========================================================================================================
def eval_turbofan_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    engine = aircraft.turbofan_engine

    engine.rating_factor = (0.800,0.688,0.624,0.560,0.100)

    return


#===========================================================================================================
def eval_turbofan_nacelle_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    fuselage = aircraft.fuselage
    wing = aircraft.wing
    engine = aircraft.turbofan_engine

    nacelle = aircraft.turbofan_nacelle

    nacelle.width = 0.49 * engine.bpr ** 0.67 + 4.8E-6 * engine.reference_thrust        # statistical regression

    nacelle.length = 0.86 * nacelle.width + engine.bpr ** 0.37      # statistical regression

    Knac = math.pi * nacelle.width * nacelle.length

    nacelle.net_wetted_area = Knac *(1.48 - 0.0076*Knac)*engine.n_engine        # statistical regression

    if nacelle.attachment == 1 :
        nacelle.y_ext = 0.7 * fuselage.width + 1.4 * nacelle.width      # statistical regression

        nacelle.z_ext = - 0.5 * fuselage.height \
                    + (nacelle.y_ext - 0.5 * fuselage.width) * math.tan(wing.dihedral) \
                    - 0.5*nacelle.width
    else:
        nacelle.y_ext = 0.5 * fuselage.width + 0.6 * nacelle.width      # statistical regression

        nacelle.z_ext = 0.5 * fuselage.height

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

    return


#===========================================================================================================
def eval_turbofan_nacelle_mass(aircraft):
    """
    Thermal propulsive nacelle mass estimation
    """

    engine = aircraft.turbofan_engine

    nacelle = aircraft.turbofan_nacelle

    nacelle.mass = (1250. + 0.021*engine.reference_thrust)*engine.n_engine       # statistical regression

    nacelle.c_g = nacelle.x_ext + 0.7 * nacelle.length      # statistical regression

    return


