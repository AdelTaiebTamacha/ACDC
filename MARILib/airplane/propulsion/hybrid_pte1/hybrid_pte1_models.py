#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import numpy

from MARILib.earth import environment as earth

from MARILib.airplane.propulsion import propulsion_models as propu

#===========================================================================================================
def hybrid_sfc(aircraft,pamb,tamb,mach,rating,nei):
    """
    Bucket SFC for a turbofan
    """

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    power_elec = aircraft.power_elec_chain
    e_engine = aircraft.electric_engine
    e_nacelle = aircraft.electric_nacelle

    power_ratio = numpy.array([e_engine.mto_e_power_ratio,
                               e_engine.mcn_e_power_ratio,
                               e_engine.mcl_e_power_ratio,
                               e_engine.mcr_e_power_ratio,
                               e_engine.fid_e_power_ratio])

    sfc0 = ( 0.4 + 1/engine.bpr**0.895 )/36000

    kC = engine.core_thrust_ratio
    kW = power_ratio[rating]

    if (propulsion.bli_effect>0):
        kBLIe = propulsion.bli_e_thrust_factor
    else:
        kBLIe = 1.

    eff_prop = nacelle.efficiency_prop
    eff_e_prop = e_nacelle.efficiency_prop
    eff_chain = power_elec.overall_efficiency

    eff_h = kC + (1-kC)*( kW*kBLIe*(eff_e_prop/eff_prop)*eff_chain + (1-kW) )

    sfc = sfc0 / eff_h

    return sfc


#===========================================================================================================
def hybrid_thrust(aircraft,Pamb,Tamb,Mach,rating,nei):

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    battery = aircraft.battery
    power_elec = aircraft.power_elec_chain
    e_engine = aircraft.electric_engine
    e_nacelle = aircraft.electric_nacelle

    power_ratio = numpy.array([e_engine.mto_e_power_ratio,
                               e_engine.mcn_e_power_ratio,
                               e_engine.mcl_e_power_ratio,
                               e_engine.mcr_e_power_ratio,
                               e_engine.fid_e_power_ratio])

    # Battery power feed is used in temporary phases only
    battery_power_feed = numpy.array([1,0,1,0,0])*battery.power_feed \
                                                 *e_nacelle.controler_efficiency \
                                                 *e_nacelle.motor_efficiency

    fn,data = propu.turbofan_thrust(aircraft,Pamb,Tamb,Mach,rating,nei)
    (fn_core,fn_fan0,fn0,shaft_power0) = data

    Vsnd = earth.sound_speed(Tamb)

    Vair = Vsnd*Mach

    shaft_power1 = (1-power_ratio[rating])*shaft_power0     # Shaft power dedicated to the fan

    fn_fan1 = nacelle.efficiency_prop*shaft_power1/Vair     # Effective fan thrust

    shaft_power2 = power_ratio[rating]*shaft_power0*(engine.n_engine - nei)     # Shaft power dedicated to electric generator

    # Effective eFan shaft power
    pw_shaft2 =   power_elec.overall_efficiency*shaft_power2 \
                + e_nacelle.motor_efficiency*e_nacelle.controler_efficiency*battery_power_feed[rating]

    if (pw_shaft2 > 0.):

        if (propulsion.bli_effect>0):
            (fn_fan2,q1,dVbli) = propu.fan_thrust_with_bli(e_nacelle,Pamb,Tamb,Mach,Vair,pw_shaft2)
            dVbli_o_V = dVbli/Vair
        else:
            (fn_fan2,q0) = propu.fan_thrust(e_nacelle,Pamb,Tamb,Mach,Vair,pw_shaft2)
            dVbli_o_V = 0.

        sec = (pw_shaft2/e_nacelle.motor_efficiency)/fn_fan2

    else:

        dVbli_o_V  = 0.
        fn_fan2 = 0.
        sec = 0

    fn = (fn_core + fn_fan1)*(engine.n_engine - nei) + fn_fan2

    data = (fn_core,fn_fan1,fn_fan2,dVbli_o_V,shaft_power2,fn0,shaft_power0)

    return (fn,sec,data)


#===========================================================================================================
def electric_nacelle_drag(aircraft,nacelle,Re,Mach):
    """
    Turbofan nacelle drag
    """

    wing = aircraft.wing

    fac = (1 + 0.126*Mach**2)

    # All nacelle drag
    nac_nwa = nacelle.net_wetted_area
    nac_cxf =   1.15*((0.455/fac)*(numpy.log(10)/numpy.log(Re*nacelle.length))**2.58)*nac_nwa/wing.area

    return nac_cxf,nac_nwa


