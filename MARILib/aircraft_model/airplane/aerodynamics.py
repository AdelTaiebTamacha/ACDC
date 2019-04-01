#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

import numpy
import math

from MARILib.tools.math import maximize_1d, lin_interp_1d

from MARILib.earth import environment as earth
from MARILib.airplane.propulsion import propulsion_models as propu


#===========================================================================================================
def reynolds_number(pamb,tamb,mach):
    """
    Reynolds number
    """
    fac = ( 1 + 0.126*mach**2 )
    re = 47899*pamb*mach*(fac*tamb + 110.4) / (tamb**2 * fac**2.5)
    return re


#===========================================================================================================
def high_lift(wing,hld_conf):
    """
# 0 =< HLDtype =< 10
# 0 =< HLDconf =< 1
# Typically : HLDconf = 1 ==> CzmaxLD
#           : HLDconf = 0.1 to 0.5 ==> CzmaxTO
    """

    cz_max_2d = {
            0 : 1.45 ,  # Clean
            1 : 2.25 ,  # Flap only, Rotation without slot
            2 : 2.60 ,  # Flap only, Rotation with slot      (ATR)
            3 : 2.80 ,  # Flap only, Rotation double slot
            4 : 2.80 ,  # Flap only, Fowler
            5 : 2.00 ,  # Slap only
            6 : 2.45 ,  # Slat + Flap rotation double slot
            7 : 2.60 ,  # Slat + Flap rotation with slot
            8 : 2.90 ,  # Slat + Flap rotation double slot
            9 : 3.00 ,  # Slat + Fowler                      (A320)
            10 : 3.20,  # Slat + Fowler + Fowler double slot (A321)
            }.get(wing.hld_type, "Erreur - high_lift_, HLDtype out of range")    # 9 is default if x not found

    cz_max_ld = cz_max_2d * numpy.cos(wing.sweep)

    if(wing.hld_type<5):
        cz_max_base = 1.45  # Flap only
    else:
        if(hld_conf==0):
            cz_max_base = 1.45  # Clean
        else:
            cz_max_base = 2.00  # Slat + Flap

    cz_max = (1-hld_conf)*cz_max_base + hld_conf*cz_max_ld

    cz_0 = cz_max - cz_max_base  # Assumed the Lift vs AoA is just translated upward and Cz0 clean equal to zero

    return cz_max, cz_0


#===========================================================================================================
def drag(aircraft, pamb, tamb, mach, cz):
    """
    Total aircraft drag with the assumption that the wing takes all the lift
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing
    htp = aircraft.horizontal_tail
    vtp = aircraft.vertical_tail

    # Form and friction drag
    #-----------------------------------------------------------------------------------------------------------
    re = reynolds_number(pamb,tamb,mach)

    nac_cxf,nac_nwa = propu.nacelle_drag(aircraft,re,mach)

    fac = ( 1 + 0.126*mach**2 )

    fuse_cxf = 1.05*((0.455/fac)*(numpy.log(10)/numpy.log(re*fuselage.length))**2.58 ) * fuselage.net_wetted_area / wing.area
    wing_cxf = 1.4*((0.455/fac)*(numpy.log(10)/numpy.log(re*wing.mac))**2.58) * wing.net_wetted_area / wing.area
    htp_cxf = 1.4*((0.455/fac)*(numpy.log(10)/numpy.log(re*htp.mac))**2.58) * htp.net_wetted_area / wing.area
    vtp_cxf = 1.4*((0.455/fac)*(numpy.log(10)/numpy.log(re*vtp.mac))**2.58) * vtp.net_wetted_area / wing.area

    aircraft_cxf = fuse_cxf + wing_cxf + htp_cxf + vtp_cxf + nac_cxf

    # Parasitic drag (seals, antennas, sensors, ...)
    #-----------------------------------------------------------------------------------------------------------
    aircraft_net_wetted_area =   fuselage.net_wetted_area + wing.net_wetted_area + htp.net_wetted_area + vtp.net_wetted_area \
                               + nac_nwa

    aircraft_knwa = aircraft_net_wetted_area/1000

    aircraft_kp = (0.0247*aircraft_knwa - 0.11)*aircraft_knwa + 0.166       # Parasitic drag factor

    aircraft_cx_parasitic = aircraft_cxf*aircraft_kp

    # Additional drag
    #-----------------------------------------------------------------------------------------------------------
    X = numpy.array([1.0, 1.5, 2.4, 3.3, 4.0, 5.0])
    Y = numpy.array([0.036, 0.020, 0.0075, 0.0025, 0, 0])

    param = fuselage.tail_cone_length/fuselage.width

    fuse_cx_tapered = lin_interp_1d(param,X,Y)     # Tapered fuselage drag (tail cone)

    # Total zero lift drag
    #-----------------------------------------------------------------------------------------------------------
    aircraft_cx0 = aircraft_cxf + aircraft_cx_parasitic + fuse_cx_tapered

    # Induced drag
    #-----------------------------------------------------------------------------------------------------------
    ki = ((fuselage.width / wing.span)**2 + 1.05 )  / (numpy.pi * wing.aspect_ratio)
    aircraft_cxi = ki*cz**2  # Induced drag

    # Compressibility drag
    #-----------------------------------------------------------------------------------------------------------
    # Freely inspired from Korn equation
    cz_design = 0.5
    aircraft_mach_div = 0.03 + design_driver.cruise_mach + 0.1*(cz_design-cz)

    aircraft_cxc = 0.0025 * numpy.exp(40*(mach - aircraft_mach_div) )

    # Sum up
    #-----------------------------------------------------------------------------------------------------------
    aircraft_cx = aircraft_cx0 + aircraft_cxi + aircraft_cxc

    aircraft_lod  = cz/aircraft_cx

    return aircraft_cx, aircraft_lod


#===========================================================================================================
def lod_max(aircraft,pamb,tamb,mach):
    """
    Maximum lift to drag ratio
    """

    #=======================================================================================
    def fct_lod_max(cz,aircraft,pamb,tamb,mach):
        [cx,lod] = drag(aircraft,pamb,tamb,mach,cz)
        return lod
    #---------------------------------------------------------------------------------------
    cz_ini = 0.5
    dcz = 0.05

    fct = [fct_lod_max, 1,aircraft,pamb,tamb,mach]
    [lod_max_cz,lod_max,rc] = maximize_1d(cz_ini,dcz,fct)

    return lod_max,lod_max_cz


#===========================================================================================================
def corrected_air_flow(Ptot,Ttot,Mach):
    """
    Computes the corrected air flow per square meter
    """

    R = earth.gaz_constant()
    gam = earth.heat_ratio()

    f_M = Mach*(1. + 0.5*(gam-1)*Mach**2)**(-(gam+1.)/(2.*(gam-1.)))

    CQoA = (math.sqrt(gam/R)*Ptot/math.sqrt(Ttot))*f_M

    return CQoA


#===========================================================================================================
def boundary_layer(re,x_length):
    """
    Thickness of a turbulent boundary layer which developped turbulently from its starting point
    """

    d = (0.385*x_length)/(re*x_length)**(1/5)

    return d


#===========================================================================================================
def air_flows(rho,v_air,r,d,y):
    """
    Air flows and speeds at rear end of a cylinder of radius rear_radius mouving at v_air in the direction of its axes
    y is the elevation upon the surface of the cylinder : 0 < y < inf
    """

    # exponent in the formula of the speed profile inside a turbulent BL of thickness bly : Vy/Vair = (y/d)**(1/7)
    n = 1/7

    # Cumulated air flow at y_elev, without BL
    q0 = (2.*numpy.pi)*(rho*v_air)*(r*y + (1/2)*y**2)

    ym = min(y,d)

    # Cumulated air flow at y_elev, with BL
    q1 = (2.*numpy.pi)*(rho*v_air)*d*( (r/(n+1))*(ym/d)**(n+1) + (d/(n+2))*(ym/d)**(n+2) )

    if (y>d):
        # Add to Q1 the air flow outside the BL
        q1 = q1 + q0 - (2.*numpy.pi)*(rho*v_air)*( r*d + (1/2)*d**2 )

    q2 = q1 - q0        # Cumulated air flow at y_elev, inside the BL (going speed wise)

    v1 = v_air*(q1/q0)     # Mean speed of q1 air flow at y_elev

    dv = v_air - v1       # Mean air flow speed variation at y_elev

    return q0,q1,q2,v1,dv


#===========================================================================================================
def specific_air_flows(r,d,y):
    """
    Specific air flows and speeds at rear end of a cylinder of radius R
    mouving at Vair in the direction of Qs = Q/(rho*Vair)     Vs = V/Vair
    its axes, y is the elevation upon the surface of the cylinder :
                              0 < y < inf
    WARNING : even if all mass flows are positive,
    Q0 and Q1 are going backward in fuselage frame, Q2 is going forward
    in ground frame
    """

    n = 1/7 # exponent in the formula of the speed profile inside a turbulent
            # BL of thickness d : Vy/Vair = (y/d)^(1/7)

    q0s = (2.*math.pi)*( r*y + (1/2)*y**2 )
                            # Cumulated specific air flow at y, without BL
    ym = min(y,d)

    q1s = (2.*math.pi)*d*( (r/(n+1))*(ym/d)**(n+1) + (d/(n+2))*(ym/d)**(n+2) )
                            # Cumulated specific air flow at y, without BL
    if y>d:

        q1s = q1s + q0s - (2.*math.pi)*( r*d + (1/2)*d**2 )
                            # Add to Q1 the specific air flow outside the BL
    q2s = q0s - q1s
        # Cumulated specific air flow at y, inside the BL (going speed wise)
    v1s = (q1s/q0s) # Mean specific speed of Q1 air flow at y

    dVs = (1 - v1s) # Mean specific air flow spped variation at y

    return (q0s,q1s,q2s,v1s,dVs)

