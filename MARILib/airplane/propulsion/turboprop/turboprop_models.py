#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import numpy

from MARILib.earth import environment as earth


#===========================================================================================================
def turboprop_sfc(aircraft,pamb,tamb,mach,rating,nei):
    """
    Bucket SFC for a turbofan
    """

    engine = aircraft.turboprop_engine

    Vsnd = earth.sound_speed(tamb)
    Vair = Vsnd*mach

    psfc = 0.4*1.68969e-07   # 0.4 lb/shp/h

    tsfc = psfc*Vair/engine.propeller_efficiency

    return tsfc


#===========================================================================================================
def turboprop_thrust(aircraft,Pamb,Tamb,Mach,rating,nei):
    """
    Calculation of thrust for pure turboprop airplane
    Warning : ALL engine thrust returned
    """

    engine = aircraft.turboprop_engine

    factor = engine.rating_factor       # [MTO,MCN,MCL,MCR,FID]

    (rho,sig) = earth.air_density(Pamb,Tamb)

    shaft_power = factor[rating]*engine.reference_power*sig

    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*Mach

    fn0 = shaft_power*engine.propeller_efficiency/Vair

    fn = fn0*(engine.n_engine - nei)        # All turboprop thrust

    data = (fn0,shaft_power)   # Data for ONE turboprop engine

    return fn,data


#===========================================================================================================
def turboprop_nacelle_drag(aircraft,nacelle,Re,Mach):
    """
    Turboprop nacelle drag
    """

    wing = aircraft.wing

    fac = (1 + 0.126*Mach**2)

    # All nacelle drag
    nac_nwa = nacelle.net_wetted_area
    nac_cxf =   1.15*((0.455/fac)*(numpy.log(10)/numpy.log(Re*nacelle.length))**2.58)*nac_nwa/wing.area

    return nac_cxf,nac_nwa


#===========================================================================================================
def turboprop_oei_drag(aircraft,nacelle,pamb,tamb):
    """
    Inoperative engine drag coefficient (to be updated)
    """

    wing = aircraft.wing

    dCx = 0.06*nacelle.width**2 / wing.area

    return dCx



