#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

import numpy
import math

from MARILib.tools import units as unit

from MARILib.earth import environment as earth

from MARILib.aircraft_model.airplane import aerodynamics as craft_aero
from MARILib.airplane.airframe import airframe_models as frame_aero


#===========================================================================================================
def eval_cabin_design(aircraft):
    """
    Cabin pre design
    """

    cabin = aircraft.cabin

    # Distance between nose and cabin forward limit (Design rule)
    cabin.fwd_limit = 4.0

    # Statistical regression
    cabin.width = 0.38*cabin.n_pax_front + 1.05*cabin.n_aisle + 0.15

    # Statistical regression
    cabin.length = 6.3*(cabin.width + 0.4) \
                   + 0.005*(cabin.n_pax_ref/cabin.n_pax_front)**2.25 \
                   - cabin.fwd_limit

    # Factor 0.95 accounts for tappered parts
    cabin.floor_area = 0.95*cabin.length*cabin.width

    return


#===========================================================================================================
def eval_cabin_mass(aircraft):
    """
    Cabin mass & CG estimation
    """

    design_driver = aircraft.design_driver

    cabin = aircraft.cabin

    cabin.m_furnishing = 0.063*cabin.n_pax_ref**2 + 9.76*cabin.n_pax_ref        # Furnishings

    cabin.m_op_item = 5.2*(cabin.n_pax_ref*design_driver.design_range*1e-6)     # Operator items

    cabin.cg_furnishing = cabin.fwd_limit + 0.55*cabin.length                   # Rear cabin is heavier because of higher density

    cabin.cg_op_item = cabin.cg_furnishing

    return


#===========================================================================================================
def eval_fuselage_design(aircraft):
    """
    Fuselage pre design
    """

    cabin = aircraft.cabin

    fuselage = aircraft.fuselage

    # Fuselage walls are supposed 0.2m thick
    fuselage.width = cabin.width + 0.4

    fuselage.height = 1.25*(cabin.width - 0.15)

    fuselage.length = cabin.fwd_limit + cabin.length + 1.5*fuselage.width

    fuselage.tail_cone_length = 3.45*fuselage.width

    fuselage.net_wetted_area = 2.7*fuselage.length*numpy.sqrt(fuselage.width*fuselage.height)

    return


#===========================================================================================================
def eval_fuselage_mass(aircraft):
    """
    Fuselage mass & CG estimation
    """

    fuselage = aircraft.fuselage

    kfus = numpy.pi*fuselage.length*numpy.sqrt(fuselage.width*fuselage.height)  ;

    fuselage.mass = 5.47*kfus**1.2 ;  # Statistical regression versus fuselage built surface

    fuselage.c_g = 0.50*fuselage.length ;  # Middle of the fuselage

    return


#===========================================================================================================
def eval_pre_design_vtp(aircraft):
    """
    VTP predesign
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    nacelle = aircraft.turbofan_nacelle
    engine = aircraft.turbofan_engine

    vtp = aircraft.vertical_tail

    vtp.taper_ratio = 0.40      # Design rule
    vtp.aspect_ratio = 1.7      # Design rule
    vtp.t_o_c = 0.10            # Design rule

    wing_sweep = 1.6*max(0,(design_driver.cruise_mach-0.5))     # Empirical law
    vtp.sweep = wing_sweep + unit.rad_deg(10)        # Design rule

#    if(HQsizing==0):  # Statistical sizing

    vtp.volume = 0.4   # Design rule
    vtp.lever_arm = 3.4e-4*fuselage.length**2 + 0.42*fuselage.length
    vtp.area = vtp.volume*(1e-3*engine.reference_thrust*nacelle.y_ext)/vtp.lever_arm

#    else:  # HQ sizing
#        vtp.lever_arm = VtpLarm0 - dXwing ;
#        vtp.area = VtpArea0 + dVtpArea ;
#        vtp.volume = (vtp.area*vtp.lever_arm)/(1e-3*aircraft.reference_thrust*nac_y_ext_i) ;

    vtp.net_wetted_area = 2.01*vtp.area     # Factor 2.01 accounts for carmans

    vtp.height = numpy.sqrt(vtp.aspect_ratio*vtp.area)
    vtp.c_root = 2*vtp.area/(vtp.height*(1+vtp.taper_ratio))
    vtp.c_tip = vtp.taper_ratio*vtp.c_root

    vtp_x_wise_anchor = 0.825       # Locate VTP versus end fuselage length

    vtp.x_root = fuselage.length*(1-fuselage.tail_cone_length/fuselage.length*(1-vtp_x_wise_anchor)) - vtp.c_root
    vtp.x_tip = vtp.x_root + 0.25*(vtp.c_root-vtp.c_tip) + vtp.height*numpy.tan(vtp.sweep)

    vtp.z_root = fuselage.height
    vtp.z_tip = vtp.z_root+vtp.height

    vtp.mac = vtp.height*(vtp.c_root**2+vtp.c_tip**2+vtp.c_root*vtp.c_tip)/(3*vtp.area)
    vtp.x_mac = vtp.x_root+(vtp.x_tip-vtp.x_root)*vtp.height*(2*vtp.c_tip+vtp.c_root)/(6*vtp.area)

    return


#===========================================================================================================
def eval_vtp_mass(aircraft):
    """
    VTP mass & CG estimation
    """

    htp = aircraft.horizontal_tail

    vtp = aircraft.vertical_tail

    if(htp.attachment==1):        # Classical
        vtp.mass = 25 * vtp.area
    elif(htp.attachment==2):      # T-tail
        vtp.mass = 28 * vtp.area  # Heavier because of HTP on top

    vtp.c_g = vtp.x_mac + 0.20*vtp.mac

    return


#===========================================================================================================
def eval_pre_design_htp(aircraft):
    """
    HTP predesign
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing
    vtp = aircraft.vertical_tail

    htp = aircraft.horizontal_tail

    htp.taper_ratio = 0.35      # Design rule
    htp.aspect_ratio = 5        # Design rule
    htp.t_o_c = 0.10            # Design rule

    wing_sweep = 1.6*max(0,(design_driver.cruise_mach-0.5))     # Empirical law
    htp.sweep = wing_sweep + unit.rad_deg(5)     # Design rule

    htp.dihedral = unit.rad_deg(5)       # HTP dihedral

#    if(HQsizing==0):  # Statistical sizing

    htp.volume = 0.94       # Design rule
    if(htp.attachment==1):    # Classical
        htp.lever_arm = vtp.lever_arm - 0.20*numpy.tan(vtp.sweep)*vtp.height
    elif(htp.attachment==2):  # T-tail
        htp.lever_arm = vtp.lever_arm + 0.80*numpy.tan(vtp.sweep)*vtp.height
    else:
        print("airframe_predesign_, HtpType value is out of range")

    htp.area = htp.volume*(wing.area*wing.mac/htp.lever_arm)

#    else:  # HQ sizing
#        htp.area = HtpArea0 + dHtpArea ;
#        htp.lever_arm = HtpLarm0 - dXwing ;
#        htp.volume = (htp.area*htp.lever_arm)/(wing.area*wing.mac)

    if(htp.attachment==1):      # Classical
        htp.net_wetted_area = 1.63*htp.area
    elif(htp.attachment==2):        # T-tail
        htp.net_wetted_area = 2.01*htp.area
    else:
        print("airframe_predesign_, htp.attachment value is out of range")

    htp.span = numpy.sqrt(htp.aspect_ratio*htp.area)

    htp.c_axe = 2*htp.area/(htp.span*(1+htp.taper_ratio))
    htp.c_tip = htp.taper_ratio*htp.c_axe
    htp.y_tip = 0.5*htp.span

    htp.mac = htp.span*(htp.c_axe**2+htp.c_tip**2+htp.c_axe*htp.c_tip)/(3*htp.area)
    htp.y_mac = htp.y_tip**2*(2*htp.c_tip+htp.c_axe)/(3*htp.area)
    x_tip_local = 0.25*(htp.c_axe-htp.c_tip) + htp.y_tip*numpy.tan(htp.sweep)
    x_mac_local = htp.y_tip*x_tip_local*(htp.c_tip*2.+htp.c_axe)/(3*htp.area)

    htp.x_mac = vtp.x_mac + 0.25*vtp.mac - vtp.lever_arm + htp.lever_arm - 0.25*htp.mac
    htp.x_axe = htp.x_mac - x_mac_local
    htp.x_tip = htp.x_axe + x_tip_local

    htp_z_wise_anchor = 0.80       # Locate HTP versus end fuselage height

    if(htp.attachment==1):      # Classical
        htp.z_axe = htp_z_wise_anchor*fuselage.height
    elif(htp.attachment==2):    # T-tail
        htp.z_axe = fuselage.height + vtp.height
    else:
        print("airframe_predesign_, htp.attachment value is out of range")

    htp.z_tip = htp.z_axe + htp.y_tip*math.tan(htp.dihedral)

    return


#===========================================================================================================
def eval_htp_mass(aircraft):
    """
    HTP mass & CG estimation
    """

    htp = aircraft.horizontal_tail

    htp.mass = 22 * htp.area      # Statistical regression

    htp.c_g = htp.x_mac + 0.20*htp.mac

    return


#===========================================================================================================
def eval_pre_design_wing(aircraft):
    """
    Wing predesign
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    vtp = aircraft.vertical_tail
    nacelle = aircraft.turbofan_nacelle
    weights = aircraft.weights

    wing = aircraft.wing

    wing.t_o_c_t = 0.10
    wing.t_o_c_k = wing.t_o_c_t + 0.01
    wing.t_o_c_r = wing.t_o_c_k + 0.03

    wing.sweep = 1.6*max(0,(design_driver.cruise_mach-0.5))     # Empirical law

    wing.dihedral = unit.rad_deg(5)

    if(wing.morphing==1):   # Aspect ratio is driving parameter
        wing.span = numpy.sqrt(wing.aspect_ratio*wing.area)
    elif(wing.morphing==2): # Span is driving parameter
        wing.aspect_ratio = wing.span**2/wing.area
    else:
        print("geometry_predesign_, wing_morphing index is unkown")

    # Correlation between span loading and tapper ratio
    wing.taper_ratio = 0.3 - 0.025*(1e-3*weights.mtow/wing.span)

    # Factor 1.64 accounts for the part of HTP ref area hidden in the fuselage
    wing.net_wetted_area = 1.64*wing.area

    wing.y_kink = nacelle.y_ext
    wing.y_root = 0.5*fuselage.width
    wing.y_tip = 0.5*wing.span

    if(15<unit.deg_rad(wing.sweep)):  # With kink
      Phi100intTE = max( 0 , 2*(wing.sweep-unit.rad_deg(32)) )
      tan_phi100 = numpy.tan(Phi100intTE)
      A = ((1-0.25*wing.taper_ratio)*wing.y_kink+0.25*wing.taper_ratio*wing.y_root-wing.y_tip) / (0.75*wing.y_kink+0.25*wing.y_root-wing.y_tip)
      B = (numpy.tan(wing.sweep)-tan_phi100) * ((wing.y_tip-wing.y_kink)*(wing.y_kink-wing.y_root)) / (0.25*wing.y_root+0.75*wing.y_kink-wing.y_tip)
      wing.c_root = (wing.area-B*(wing.y_tip-wing.y_root)) / (wing.y_root+wing.y_kink+A*(wing.y_tip-wing.y_root)+wing.taper_ratio*(wing.y_tip-wing.y_kink))
      wing.c_kink = A*wing.c_root + B
      wing.c_tip = wing.taper_ratio*wing.c_root

    else:		# Without kink
      wing.c_root = 2*wing.area / (2*wing.y_root*(1-wing.taper_ratio) + (1+wing.taper_ratio)*numpy.sqrt(wing.aspect_ratio*wing.area))
      wing.c_tip = wing.taper_ratio*wing.c_root
      wing.c_kink = ((wing.y_tip-wing.y_kink)*wing.c_root + (wing.y_kink-wing.y_root)*wing.c_tip) / (wing.y_tip-wing.y_root)


    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    wing.mac = 2.*(3*wing.y_root*wing.c_root**2 \
            +(wing.y_kink-wing.y_root)*(wing.c_root**2+wing.c_kink**2+wing.c_root*wing.c_kink) \
            +(wing.y_tip-wing.y_kink)*(wing.c_kink**2+wing.c_tip**2+wing.c_kink*wing.c_tip) \
            )/(3*wing.area)

    wing.y_mac = ( 3*wing.c_root*wing.y_root**2 \
              +(wing.y_kink-wing.y_root)*(wing.c_kink*(wing.y_root+wing.y_kink*2.)+wing.c_root*(wing.y_kink+wing.y_root*2.)) \
              +(wing.y_tip-wing.y_kink)*(wing.c_tip*(wing.y_kink+wing.y_tip*2.)+wing.c_kink*(wing.y_tip+wing.y_kink*2.)) \
              )/(3*wing.area)

    x_mac_local = ( (wing.y_kink-wing.y_root)*tan_phi0*((wing.y_kink-wing.y_root)*(wing.c_kink*2.+wing.c_root) \
                   +(wing.y_tip-wing.y_kink)*(wing.c_kink*2.+wing.c_tip))+(wing.y_tip-wing.y_root)*tan_phi0*(wing.y_tip-wing.y_kink)*(wing.c_tip*2.+wing.c_kink) \
                  )/(3*wing.area)

    wing.x_root = vtp.x_mac + 0.25*vtp.mac - vtp.lever_arm - 0.25*wing.mac - x_mac_local

    wing.x_kink = wing.x_root + (wing.y_kink-wing.y_root)*tan_phi0
    wing.x_tip = wing.x_root + (wing.y_tip-wing.y_root)*tan_phi0

    wing.x_mac = wing.x_root+( (wing.x_kink-wing.x_root)*((wing.y_kink-wing.y_root)*(wing.c_kink*2.+wing.c_root) \
                        +(wing.y_tip-wing.y_kink)*(wing.c_kink*2.+wing.c_tip))+(wing.x_tip-wing.x_root)*(wing.y_tip-wing.y_kink)*(wing.c_tip*2.+wing.c_kink) \
                       )/(wing.area*3.)
    if (wing.attachment==1):
        wing.z_root = 0
    else:
        wing.z_root = fuselage.height - 0.5*wing.t_o_c_r*wing.c_root

    wing.z_kink = wing.z_root+(wing.y_kink-wing.y_root)*numpy.tan(wing.dihedral)
    wing.z_tip = wing.z_root+(wing.y_tip-wing.y_root)*numpy.tan(wing.dihedral)

    # Wing setting
    #-----------------------------------------------------------------------------------------------------------
    g = earth.gravity()
    gam = earth.heat_ratio()

    disa = 0
    rca = design_driver.ref_cruise_altp
    mach = design_driver.cruise_mach
    mass = 0.95*weights.mtow

    pamb,tamb,tstd,dtodz = earth.atmosphere(rca,disa)

    cza_wo_htp = frame_aero.cza_wo_htp(mach, fuselage.width, wing.aspect_ratio, wing.span, wing.sweep)

    # AoA = 2.5° at cruise start
    wing.setting = (0.97*mass*g)/(0.5*gam*pamb*mach**2*wing.area*cza_wo_htp) - unit.rad_deg(2.5)

    return


#===========================================================================================================
def eval_wing_mass(aircraft):
    """
    Wing mass & CG estimation
    """

    weights = aircraft.weights
    aerodynamics = aircraft.aerodynamics

    wing = aircraft.wing

    (cz_max_ld,cz0) = craft_aero.high_lift(wing, aerodynamics.hld_conf_ld)

    A = 32*wing.area**1.1
    B = 4*wing.span**2 * numpy.sqrt(weights.mtow*weights.mzfw)
    C = 1.1e-6*(1+2*wing.aspect_ratio)/(1+wing.aspect_ratio)
    D = (0.6*wing.t_o_c_r+0.3*wing.t_o_c_k+0.1*wing.t_o_c_t) * (wing.area/wing.span)
    E = numpy.cos(wing.sweep)**2
    F = 1200*(cz_max_ld - 1.8)**1.5

    wing.mass = A + (B*C)/(D*E) + F   # Shevell formula + high lift device regression

    wing.c_g =  0.25*(wing.x_root + 0.40*wing.c_root) \
              + 0.55*(wing.x_kink + 0.40*wing.c_kink) \
              + 0.20*(wing.x_tip + 0.40*wing.c_tip)

    return


#===========================================================================================================
def eval_landing_gear_mass(aircraft):
    """
    Landing Gears mass & CG estimation
    """
    wing = aircraft.wing
    weights = aircraft.weights

    ldg = aircraft.landing_gears

    ldg.mass = 0.02*weights.mtow**1.03 + 0.012*weights.mlw    # Landing gears

    ldg.c_g = wing.x_root + 0.70*wing.c_root

    return


#===========================================================================================================
def eval_tank_data(aircraft):
    """
    Tank predesign
    """

    propulsion = aircraft.propulsion
    fuselage = aircraft.fuselage
    wing = aircraft.wing

    tanks = aircraft.tanks

    tanks.cantilever_volume = 0.2 * (wing.area*wing.mac*(0.5*wing.t_o_c_r + 0.3*wing.t_o_c_k + 0.2*wing.t_o_c_t))

    tanks.central_volume = 1.3 * fuselage.width * wing.t_o_c_r * wing.mac**2

    tanks.fuel_density = earth.fuel_density(propulsion.fuel_type)

    tanks.mfw_volume_limited = (tanks.central_volume + tanks.cantilever_volume)*tanks.fuel_density

    tanks.fuel_cantilever_cg =  0.25*(wing.x_root + 0.40*wing.c_root) \
                              + 0.65*(wing.x_kink + 0.40*wing.c_kink) \
                              + 0.10*(wing.x_tip + 0.40*wing.c_tip)

    tanks.fuel_central_cg = wing.x_root + 0.30*wing.c_root

    tanks.fuel_max_fwd_cg = tanks.fuel_central_cg    # Fuel max forward CG, central tank is forward only within backward swept wing
    tanks.fuel_max_fwd_mass = tanks.central_volume*tanks.fuel_density

    tanks.fuel_max_bwd_cg = tanks.fuel_cantilever_cg    # Fuel max Backward CG
    tanks.fuel_max_bwd_mass = tanks.cantilever_volume*tanks.fuel_density

    return


