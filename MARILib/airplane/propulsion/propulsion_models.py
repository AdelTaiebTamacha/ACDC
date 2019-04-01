#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import math
from scipy.optimize import fsolve
from MARILib.tools.math import lin_interp_1d

from MARILib.earth import environment as earth

from MARILib.aircraft_model.airplane import aerodynamics as craft_aero, regulation as regul

from MARILib.airplane.propulsion.turbofan.turbofan_models \
    import turbofan_sfc, turbofan_thrust, turbofan_nacelle_drag, \
        turbofan_oei_drag

from MARILib.airplane.propulsion.turboprop.turboprop_models \
    import turboprop_sfc, turboprop_thrust, turboprop_nacelle_drag, \
        turboprop_oei_drag

from MARILib.airplane.propulsion.hybrid_pte1.hybrid_pte1_models \
    import hybrid_sfc, hybrid_thrust, electric_nacelle_drag

from MARILib.airplane.propulsion.hybrid_pte2.hybrid_pte2_models \
    import hybrid_body_sfc, hybrid_body_thrust, body_nacelle_drag


#===========================================================================================================
def sfc(aircraft,pamb,tamb,mach,rating,nei):
    """
    Bucket SFC for a turbofan
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture==1):

        sfc = turbofan_sfc(aircraft,pamb,tamb,mach,rating,nei)

    elif (propulsion.architecture==2):

        sfc = hybrid_sfc(aircraft,pamb,tamb,mach,rating,nei)

    elif (propulsion.architecture==3):

        sfc = hybrid_body_sfc(aircraft,pamb,tamb,mach,rating,nei)

    elif (propulsion.architecture==4):

        sfc = turboprop_sfc(aircraft,pamb,tamb,mach,rating,nei)

    else:
        raise Exception("propulsion.architecture index is out of range")

    return sfc


#===========================================================================================================
def thrust(aircraft,Pamb,Tamb,Mach,rating,nei):
    """
    Calculation of thrust for pure turbofan airplane
    Warning : ALL engine thrust returned
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture==1):

        fn,data = turbofan_thrust(aircraft,Pamb,Tamb,Mach,rating,nei)

    elif (propulsion.architecture==2):

        fn,sec,data = hybrid_thrust(aircraft,Pamb,Tamb,Mach,rating,nei)

    elif (propulsion.architecture==3):

        fn,sec,data = hybrid_body_thrust(aircraft,Pamb,Tamb,Mach,rating,nei)

    elif (propulsion.architecture==4):

        fn,data = turboprop_thrust(aircraft,Pamb,Tamb,Mach,rating,nei)

    else:
        raise Exception("propulsion.architecture index is out of range")

    return fn,data


#===========================================================================================================
def nacelle_drag(aircraft,Re,Mach):
    """
    All nacelle drag coefficient
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture==1):

        nacelle = aircraft.turbofan_nacelle
        nacelle_cxf,nacelle_nwa = turbofan_nacelle_drag(aircraft,nacelle,Re,Mach)

    elif (propulsion.architecture==2):

        nacelle = aircraft.turbofan_nacelle
        t_nacelle_cxf,t_nacelle_nwa = turbofan_nacelle_drag(aircraft,nacelle,Re,Mach)

        nacelle = aircraft.electric_nacelle
        e_nacelle_cxf,e_nacelle_nwa = electric_nacelle_drag(aircraft,nacelle,Re,Mach)

        nacelle_cxf = t_nacelle_cxf + e_nacelle_cxf
        nacelle_nwa = t_nacelle_nwa + e_nacelle_nwa

    elif (propulsion.architecture==3):

        nacelle = aircraft.turbofan_nacelle
        t_nacelle_cxf,t_nacelle_nwa = turbofan_nacelle_drag(aircraft,nacelle,Re,Mach)

        body = aircraft.body_nacelle
        body_cxf,body_nwa = body_nacelle_drag(aircraft,body,Re,Mach)

        nacelle = aircraft.electric_nacelle
        e_nacelle_cxf,e_nacelle_nwa = electric_nacelle_drag(aircraft,nacelle,Re,Mach)

        nacelle_cxf = body_cxf + t_nacelle_cxf + e_nacelle_cxf
        nacelle_nwa = body_nwa + t_nacelle_nwa + e_nacelle_nwa

    elif (propulsion.architecture==4):

        nacelle = aircraft.turboprop_nacelle
        nacelle_cxf,nacelle_nwa = turboprop_nacelle_drag(aircraft,nacelle,Re,Mach)

    else:
        raise Exception("propulsion.architecture index is out of range")

    return nacelle_cxf,nacelle_nwa


#===========================================================================================================
def oei_drag(aircraft,pamb,tamb):
    """
    Inoperative engine drag coefficient
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture==1):

        nacelle = aircraft.turbofan_nacelle
        dcx = turbofan_oei_drag(aircraft,nacelle,pamb,tamb)

    elif (propulsion.architecture==2):

        nacelle = aircraft.turbofan_nacelle
        dcx = turbofan_oei_drag(aircraft,nacelle,pamb,tamb)

    elif (propulsion.architecture==3):

        nacelle = aircraft.turbofan_nacelle
        dcx = turbofan_oei_drag(aircraft,nacelle,pamb,tamb)

    elif (propulsion.architecture==4):

        nacelle = aircraft.turboprop_nacelle
        dcx = turbofan_oei_drag(aircraft,nacelle,pamb,tamb)

    else:
        raise Exception("propulsion.architecture index is out of range")

    return dcx


#===========================================================================================================
def thrust_pitch_moment(aircraft,fn,pamb,mach,dcx_oei):

    propulsion = aircraft.propulsion
    wing = aircraft.wing

    gam = earth.heat_ratio()

    if (propulsion.architecture==1):

        nacelle = aircraft.turbofan_nacelle

    elif (propulsion.architecture==2):

        nacelle = aircraft.turbofan_nacelle

    elif (propulsion.architecture==3):

        nacelle = aircraft.turbofan_nacelle

    elif (propulsion.architecture==4):

        nacelle = aircraft.turboprop_nacelle

    else:
        raise Exception("propulsion.architecture index is out of range")

    cm_prop = nacelle.z_ext*(dcx_oei - fn/(0.5*gam*pamb*mach**2*wing.area))

    return cm_prop


#===========================================================================================================
def thrust_yaw_moment(aircraft,fn,pamb,mach,dcx_oei):
    """
    Assumed right engine inoperative
    """

    propulsion = aircraft.propulsion
    wing = aircraft.wing

    gam = earth.heat_ratio()

    if (propulsion.architecture==1):

        nacelle = aircraft.turbofan_nacelle

    elif (propulsion.architecture==2):

        nacelle = aircraft.turbofan_nacelle

    elif (propulsion.architecture==3):

        nacelle = aircraft.turbofan_nacelle

    elif (propulsion.architecture==4):

        nacelle = aircraft.turboprop_nacelle

    else:
        raise Exception("propulsion.architecture index is out of range")

    cn_prop = (nacelle.y_ext/wing.mac)*(fn/(0.5*gam*pamb*mach**2*wing.area) + dcx_oei)

    return cn_prop


#===========================================================================================================
def fan_thrust_with_bli(nacelle,Pamb,Tamb,Mach,Vair,PwShaft):
    """
    Compute the thrust of a fan of a given geometry swallowing
    the boundary layer (BL) of a body of a given geometry
    The amount of swallowed BL depends on the given shaft power and flying
    conditions.
    """

    bnd_layer = nacelle.bnd_layer

    gam = earth.heat_ratio()
    r = earth.gaz_constant()
    Cp = earth.heat_constant(gam,r)

    Re = craft_aero.reynolds_number(Pamb,Tamb,Mach)
    (rho,sig) = earth.air_density(Pamb,Tamb)
    Vsnd = earth.sound_speed(Tamb)

    d0 = craft_aero.boundary_layer(Re,nacelle.body_length)      # theorical thickness of the boundary layer without taking account of fuselage tapering
    r1 = 0.5*nacelle.hub_width      # Radius of the hub of the eFan nacelle
    d1 = lin_interp_1d(d0,bnd_layer[:,0],bnd_layer[:,1])     # Using the precomputed relation

    #===========================================================================================================
    def fct_power_bli(y,PwShaft,Tamb,Pamb,rho,Mach,Vair,Vsnd,r1,d1,nozzle_area):

        Ttot = earth.total_temperature(Tamb,Mach)      # Stagnation temperature at inlet position
        (q0,q1,q2,Vinlet,dVbli) = craft_aero.air_flows(rho,Vair,r1,d1,y)
        Tstat = Ttot - 0.5*Vinlet**2/Cp     # Static temperature at inlet position
        Vsnd_inlet = earth.sound_speed(Tstat)       # Sound speed at inlet position
        MachInlet = Vinlet/Vsnd_inlet        # Mean Mach number at inlet position
        PwInput = nacelle.efficiency_fan * PwShaft
        Vjet = math.sqrt(2.*PwInput/q1 + Vinlet**2)
        TtotJet = Ttot + PwShaft/(q1*Cp)        # Stagnation temperature increases due to introduced work
        Tstat = TtotJet - 0.5*Vjet**2/Cp        # Static temperature
        VsndJet = earth.sound_speed(Tstat)     # Sound speed at nozzle exhaust
        MachJet = Vjet/VsndJet                  # Mach number at nozzle output
        PtotJet = earth.total_pressure(Pamb,MachJet)       # total pressure at nozzle exhaust (P = Pamb)
        CQoA1 = craft_aero.corrected_air_flow(PtotJet,TtotJet,MachJet)    # Corrected air flow per area at fan position
        q = CQoA1*nozzle_area

        y = q1 - q

        return y
    #-----------------------------------------------------------------------------------------------------------

    nozzle_area = nacelle.nozzle_area

    fct_arg = (PwShaft,Tamb,Pamb,rho,Mach,Vair,Vsnd,r1,d1,nozzle_area)

    # Computation of y1 : thikness of the vein swallowed by the inlet
    output_dict = fsolve(fct_power_bli, x0=0.50, args=fct_arg, full_output=True)

    y = output_dict[0][0]

    Ttot = earth.total_temperature(Tamb,Mach)      # Stagnation temperature at inlet position
    (q0,q1,q2,Vinlet,dVbli) = craft_aero.air_flows(rho,Vair,r1,d1,y)
    Tstat = Ttot - 0.5*Vinlet**2/Cp     # Static temperature at inlet position
    Vsnd_inlet = earth.sound_speed(Tstat)       # Sound speed at inlet position
    MachInlet = Vinlet/Vsnd_inlet        # Mean Mach number at inlet position
    PwInput = nacelle.efficiency_fan * PwShaft
    Vjet = math.sqrt(2.*PwInput/q1 + Vinlet**2)

    eFn = q1*(Vjet - Vinlet)

    return (eFn,q1,dVbli)


#===========================================================================================================
def fan_thrust(nacelle,Pamb,Tamb,Mach,Vair,PwShaft):
    """
    Compute the thrust of a fan of given geometry swallowing
    the boundary layer (BL) of a body of given geometry
    The amount of swallowed BL depends on the given shaft power
    and flying conditions
    """

    bnd_layer = nacelle.bnd_layer

    Re = craft_aero.reynolds_number(Pamb,Tamb,Mach)
    (rho,sig) = earth.air_density(Pamb,Tamb)
    Vsnd = earth.sound_speed(Tamb)

    d0 = craft_aero.boundary_layer(Re,nacelle.body_length)      #Theoritical thickness of the boundary layer without taking account of fuselage tapering
    d1 = lin_interp_1d(d0,bnd_layer[:,0],bnd_layer[:,1])    # Using the precomputed relation
    r1 = 0.5*nacelle.hub_width        # Radius of the hub of the eFan nacelle

    PwInput = nacelle.efficiency_fan*PwShaft

    #===========================================================================================================
    def fct_power(q,PwInput,Vair,Vsnd,eNozzleArea):

        Vinlet = Vair
        MachInlet = Vinlet/Vsnd     # Mean Mach number at inlet position
        Ptot = earth.total_pressure(Pamb,MachInlet)        #total pressure at inlet position
        Ttot = earth.total_temperature(Tamb,MachInlet)     # Total temperature at inlet position
        Vjet = math.sqrt(2.*PwInput/q + Vinlet**2)
        MachJet = Vjet/Vsnd
        CQoA1 = craft_aero.corrected_air_flow(Ptot,Ttot,MachJet)       # Corrected air flow per area at fan position
        q0 = CQoA1*eNozzleArea

        y = q - q0

        return y
    #-----------------------------------------------------------------------------------------------------------

    eNozzleArea = nacelle.nozzle_area

    fct_arg = (PwInput,Vair,Vsnd,eNozzleArea)

    (q0init,q1,q2,V1,dV) = craft_aero.air_flows(rho,Vair,r1,d1,1.00)

    # Computation of y1 : thikness of the vein swallowed by the inlet
    output_dict = fsolve(fct_power, x0=q0init, args=fct_arg, full_output=True)

    q0 = output_dict[0][0]

    Vinlet = Vair
    Vjet = math.sqrt(2.*PwInput/q0 + Vinlet**2)

    eFn = q0*(Vjet - Vinlet)

    return (eFn,q0)

