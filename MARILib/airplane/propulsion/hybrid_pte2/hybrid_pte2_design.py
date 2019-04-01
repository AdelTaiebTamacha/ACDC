#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import math
import numpy

from scipy.optimize import fsolve
from MARILib.tools.math import lin_interp_1d

from MARILib.earth import environment as earth

from MARILib.aircraft_model.airplane import aerodynamics as craft_aero

from MARILib.airplane.propulsion import propulsion_models as propu

from MARILib.airplane.propulsion.hybrid_pte1.hybrid_pte1_design import eval_bli_nacelle_design


#===========================================================================================================
def eval_hybrid_turbofan_engine_design(aircraft):
    """
    Thermal hybrid propulsive architecture design
    """

    engine = aircraft.turbofan_engine

    engine.rating_factor = (1.162,1.015,0.660,0.582,0.100)      # WARNING to be investigated

    return


#===========================================================================================================
def eval_hybrid_body_nacelle_design(aircraft):
    """
    Hybrid propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing

    propulsion = aircraft.propulsion

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle
    body = aircraft.body_nacelle

    e_engine = aircraft.electric_engine
    e_nacelle = aircraft.electric_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    #-----------------------------------------------------------------------------------------------------------
    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    body.y_ext = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

    body.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.5*body.length

    body.z_ext = - 0.5 * fuselage.height \
                    + (nacelle.y_ext - 0.5 * fuselage.width) * math.tan(wing.dihedral) \
                    - 0.5*nacelle.width

    body.net_wetted_area = 2.7*body.length*body.width*engine.n_engine   # All bodies wetted area

    # Turbofan nacelles geometry is designed according to cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    dISA = 0.
    Altp = design_driver.ref_cruise_altp
    Mach = design_driver.cruise_mach
    nei = 0

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(Altp,dISA)

    fn,data = propu.turbofan_thrust(aircraft,Pamb,Tamb,Mach,MCR,nei)
    (fn_core,fn_fan0,fn0,shaft_power0) = data

    shaft_power1 = (1-e_engine.mcr_e_power_ratio)*shaft_power0     # Shaft power dedicated to the fan

    body_hub_width = body.hub_width     # Diameter of the fan hub

    body_length = body.length
    body_width = body.width

    eval_bli_nacelle_design(nacelle,Pamb,Tamb,Mach,shaft_power1,body_hub_width,body_length,body_width)

    nacelle.y_ext = body.y_ext
    nacelle.x_ext = body.x_ext + body.length
    nacelle.z_ext = body.z_ext

    # Electric nacelle is design for cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    dISA = 0.
    Altp = design_driver.ref_cruise_altp
    Mach = design_driver.cruise_mach

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(Altp,dISA)

    shaft_power = e_engine.mcr_e_shaft_power

    efan_hub_width = 0.5     # Diameter of the e fan hub
    body_length = fuselage.length
    body_width = fuselage.width

    eval_bli_nacelle_design(e_nacelle,Pamb,Tamb,Mach,shaft_power,efan_hub_width,body_length,body_width)

    e_nacelle.x_ext = fuselage.length + 0.2*e_nacelle.width
    e_nacelle.y_ext = 0
    e_nacelle.z_ext = 0.91*fuselage.height - 0.51*fuselage.height

    # Engine performance update
    #-----------------------------------------------------------------------------------------------------------
    fd = e_engine.flight_data

    e_fan_thrust = numpy.zeros(5)

    for rating in propulsion.rating_code:

        altp = fd.get("altp")[rating]
        disa = fd.get("disa")[rating]
        mach = fd.get("mach")[rating]
        nei = fd.get("nei")[rating]

        (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(altp,disa)
        (fn,sec,data) = propu.hybrid_thrust(aircraft,Pamb,Tamb,mach,rating,nei)
        (fn_core,fn_fan1,fn_fan2,dVbli_o_V,shaft_power2,fn0,shaft_power0) = data

        e_fan_thrust[rating] = fn_fan2

    e_engine.mto_e_fan_thrust = e_fan_thrust[MTO]
    e_engine.mcn_e_fan_thrust = e_fan_thrust[MCN]
    e_engine.mcl_e_fan_thrust = e_fan_thrust[MCL]
    e_engine.mcr_e_fan_thrust = e_fan_thrust[MCR]
    e_engine.fid_e_fan_thrust = e_fan_thrust[FID]

    Vair = Mach*earth.sound_speed(Tamb)

    (eFanFnBli,q1,dVbli) = propu.fan_thrust_with_bli(e_nacelle,Pamb,Tamb,Mach,Vair,shaft_power)
    (eFanFn,q0) = propu.fan_thrust(e_nacelle,Pamb,Tamb,Mach,Vair,shaft_power)
    propulsion.bli_e_thrust_factor = eFanFnBli / eFanFn     # Thrust increase due to BLI at iso shaft power

    (FanFnBli,q1,dVbli) = propu.fan_thrust_with_bli(nacelle,Pamb,Tamb,Mach,Vair,shaft_power)
    (FanFn,q0) = propu.fan_thrust(nacelle,Pamb,Tamb,Mach,Vair,shaft_power)
    propulsion.bli_thrust_factor = FanFnBli / FanFn     # Thrust increase due to BLI at iso shaft power

    return


#===========================================================================================================
def eval_hybrid_body_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimation
    """

    fuselage = aircraft.fuselage

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle
    body = aircraft.body_nacelle

    e_engine = aircraft.electric_engine
    e_nacelle = aircraft.electric_nacelle

    power_elec = aircraft.power_elec_chain

    # Body mass
    # -----------------------------------------------------------------------
    kbody = numpy.pi*body.length*body.width

    structure_mass = 5.0*kbody**1.2      # WARNING: ONE BODY MASS

    tank_mass = 0   # WARNING: ONE BODY MASS, TO BE UPDATED

    power_elec_mass = 0   # WARNING: ONE BODY MASS, TO BE UPDATED

    body.mass = (structure_mass + tank_mass + power_elec_mass)*engine.n_engine        # Betteries (if any) are accounted separetly (see aircraft.battery)

    # -----------------------------------------------------------------------
    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    e_shaft_power = numpy.array([e_engine.mto_e_shaft_power,
                                 e_engine.mcn_e_shaft_power,
                                 e_engine.mcl_e_shaft_power,
                                 e_engine.mcr_e_shaft_power,
                                 e_engine.fid_e_shaft_power])

    shaftPowerMax = max(e_shaft_power)

    turboFanMass0 = 1250 + 0.021*engine.reference_thrust # Statistical regression

    turboFanMass1 = 1250 + 0.021*engine.reference_thrust*engine.kfn_off_take

    kTurboFanMass = turboFanMass1 / turboFanMass0

    kMass = kTurboFanMass + engine.core_weight_ratio*(1-kTurboFanMass)     # Assuming core mass remains unchanged

    nacelle.mass = body.mass + engine.n_engine * turboFanMass0 * kMass     # Total engine mass

    power_elec.mass = (  1/power_elec.generator_pw_density + 1/power_elec.rectifier_pw_density \
                       + 1/power_elec.wiring_pw_density + 1/power_elec.cooling_pw_density \
                       ) * shaftPowerMax

    e_nacelle.mass = (  1/e_nacelle.controler_pw_density + 1/e_nacelle.motor_pw_density \
                      + 1/e_nacelle.nacelle_pw_density \
                      ) * shaftPowerMax

    # Propulsion system CG
    # ------------------------------------------------------------------------
    body.c_g = body.x_ext + 0.5*body.length

    nacelle.c_g = ( (nacelle.x_ext + 0.70*nacelle.length)*(nacelle.mass-body.mass) \
                   + body.c_g*body.mass \
                  )/nacelle.mass

    power_elec.c_g = 0.70*nacelle.c_g + 0.30*fuselage.length

    e_nacelle.c_g = fuselage.length + 0.5*e_nacelle.length

    return


#===========================================================================================================
def eval_body_tank_data(aircraft):
    """
    Body tank predesign
    """

    propulsion = aircraft.propulsion
    body = aircraft.body_nacelle

    fuselage = aircraft.fuselage
    wing = aircraft.wing

    tanks = aircraft.tanks

    tanks.body_volume = 0.      # TO BE UPDATED

    tanks.fuel_body_cg = 0.      # TO BE UPDATED

    tanks.fuel_density = earth.fuel_density(propulsion.fuel_type)

    tanks.mfw_volume_limited = tanks.body_volume*tanks.fuel_density

    # TO BE UPDATED (TO BE DELETED WHEN ABOVE CODE IS UPDATED)
    tanks.cantilever_volume = 0.2 * (wing.area*wing.mac*(0.5*wing.t_o_c_r + 0.3*wing.t_o_c_k + 0.2*wing.t_o_c_t))
    tanks.central_volume = 1.3 * fuselage.width * wing.t_o_c_r * wing.mac**2
    tanks.mfw_volume_limited = (tanks.central_volume + tanks.cantilever_volume)*tanks.fuel_density
    tanks.fuel_cantilever_cg =  0.25*(wing.x_root + 0.40*wing.c_root) \
                              + 0.65*(wing.x_kink + 0.40*wing.c_kink) \
                              + 0.10*(wing.x_tip + 0.40*wing.c_tip)
    tanks.fuel_central_cg = wing.x_root + 0.30*wing.c_root

    # TO BE UPDATED
    tanks.fuel_max_fwd_cg = tanks.fuel_central_cg    # Fuel max forward CG, central tank is forward only within backward swept wing
    tanks.fuel_max_fwd_mass = tanks.central_volume*tanks.fuel_density

    # TO BE UPDATED
    tanks.fuel_max_bwd_cg = tanks.fuel_cantilever_cg    # Fuel max Backward CG
    tanks.fuel_max_bwd_mass = tanks.cantilever_volume*tanks.fuel_density

    return


#===========================================================================================================
def eval_body_battery_cg(aircraft):
    """
    Body battery predesign
    """

    body = aircraft.body_nacelle

    battery = aircraft.battery

    battery.c_g = body.c_g      # TO BE UPDATED

    return









