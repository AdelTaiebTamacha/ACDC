#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

import numpy as np
from scipy.optimize import fsolve,minimize
from MARILib.tools import units as unit

from MARILib.earth import environment as earth
from MARILib.aircraft_model.airplane import airplane_design as airplane, aerodynamics as craft_aero, \
    regulation as regul
from MARILib.aircraft_model import initialization as init

from MARILib.airplane.airframe import airframe_design as airframe

from MARILib.airplane.propulsion import propulsion_design as propulsion
from MARILib.airplane.propulsion import propulsion_models as propu

from MARILib.aircraft_model.operations import handling_qualities as h_q, \
                                              mission as perfo, \
                                              environmental_impact as environ, \
                                              pricing_and_costing as costing

from MARILib.processes import solvers as sub_proc


#===========================================================================================================
def aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine):
    """
    Initialize a generic aircraft
    """

    aircraft.propulsion.architecture = propu_config

    aircraft.propulsion.fuel_type = init.fuel_type()

    aircraft.name = "my_test_airplane"
    aircraft.design_driver.design_range = design_range        # TLR
    aircraft.design_driver.cruise_mach = cruise_mach          # TLR
    aircraft.cabin.n_pax_ref = n_pax_ref                      # TLR

    aircraft.design_driver.ref_cruise_altp = init.ref_cruise_altp(propu_config)        # TLR
    aircraft.design_driver.top_of_climb_altp = init.top_of_climb_altp(propu_config)    # TLR

    aircraft.aerodynamics.hld_conf_clean = init.hld_conf_clean()
    aircraft.aerodynamics.hld_conf_ld = init.hld_conf_ld()

    aircraft.low_speed.altp_tofl = init.altp_tofl()
    aircraft.low_speed.disa_tofl = init.disa_tofl()
    aircraft.low_speed.kvs1g_tofl = regul.kvs1g_min_take_off()       # Regulation
    aircraft.low_speed.req_tofl = init.req_tofl(design_range)        # TLR

    aircraft.low_speed.altp_app_speed = init.altp_app_speed()
    aircraft.low_speed.disa_app_speed = init.disa_app_speed()
    aircraft.low_speed.kvs1g_app_speed = regul.kvs1g_min_landing()   # Regulation
    aircraft.low_speed.req_app_speed = init.req_app_speed()          # TLR

    aircraft.low_speed.disa_oei = init.disa_oei()
    aircraft.low_speed.req_oei_path = regul.ceil_oei_min_path(n_engine)     # Regulation
    aircraft.low_speed.req_oei_altp = init.req_oei_altp(propu_config)       # TLR

    aircraft.high_speed.disa_climb = init.disa_climb()
    aircraft.high_speed.req_vz_climb = init.req_vz_climb()           # TLR
    aircraft.high_speed.req_vz_cruise = init.req_vz_cruise()         # TLR
    aircraft.high_speed.req_toc_altp = init.top_of_climb_altp(propu_config)
    aircraft.high_speed.cas1_ttc = init.cas1_ttc(propu_config)
    aircraft.high_speed.cas2_ttc = init.cas2_ttc(propu_config)
    aircraft.high_speed.req_ttc = init.req_ttc()                     # TLR

    aircraft.cost_mission.disa = init.cost_mission_disa()
    aircraft.cost_mission.range = init.cost_mission_range(design_range)

    aircraft.economics.fuel_price = init.fuel_price()
    aircraft.economics.elec_price = init.elec_price()
    aircraft.economics.battery_price = init.battery_price()
    aircraft.economics.labor_cost = init.labor_cost()
    aircraft.economics.irp = init.irp()
    aircraft.economics.period = init.period()
    aircraft.economics.interest_rate = init.interest_rate()
    aircraft.economics.utilisation = init.utilisation(design_range)

    aircraft.environmental_impact.CO2_index = earth.emission_index("CO2")
    aircraft.environmental_impact.H2O_index = earth.emission_index("H2O")
    aircraft.environmental_impact.SO2_index = earth.emission_index("SO2")
    aircraft.environmental_impact.NOx_index = earth.emission_index("NOx")
    aircraft.environmental_impact.CO_index = earth.emission_index("CO")
    aircraft.environmental_impact.HC_index = earth.emission_index("HC")
    aircraft.environmental_impact.sulfuric_acid_index = earth.emission_index("sulfuric_acid")
    aircraft.environmental_impact.nitrous_acid_index = earth.emission_index("nitrous_acid")
    aircraft.environmental_impact.nitric_acid_index = earth.emission_index("nitric_acid")
    aircraft.environmental_impact.soot_index = earth.emission_index("soot")

    aircraft.weights.mzfw = init.mzfw(n_pax_ref,design_range)
    aircraft.weights.mtow = init.mtow(n_pax_ref,design_range)
    aircraft.weights.mlw = init.mlw(aircraft.weights.mtow,aircraft.weights.mzfw)

    aircraft.cabin.n_pax_front = init.n_pax_front(n_pax_ref)
    aircraft.cabin.n_aisle = init.n_aisle(aircraft.cabin.n_pax_front)

    aircraft.payload.m_pax_nominal = init.m_pax_nominal(design_range)      # TLR
    aircraft.payload.m_pax_max = init.m_pax_max(design_range)              # TLR

    aircraft.wing.attachment = init.wing_attachment()
    aircraft.wing.morphing = init.wing_morphing()
    aircraft.wing.hld_type = init.hld_type()

    aircraft.wing.area = init.wing_area(n_pax_ref,design_range)                                              # Main design variable
    aircraft.wing.aspect_ratio = init.wing_aspect_ratio()
    aircraft.wing.span = init.wing_span(aircraft.wing.area,aircraft.wing.aspect_ratio)

    aircraft.horizontal_tail.attachment = init.htp_attachment()

    if (aircraft.propulsion.architecture<4):
        aircraft.turbofan_nacelle.attachment = init.nacelle_attachment()
        aircraft.turbofan_nacelle.efficiency_fan = init.efficiency_fan()
        aircraft.turbofan_nacelle.efficiency_prop = init.efficiency_prop()

        aircraft.turbofan_engine.n_engine = n_engine
        aircraft.turbofan_engine.bpr = init.bpr()
        aircraft.turbofan_engine.reference_thrust = init.reference_thrust(n_pax_ref,design_range,n_engine)                                            # Main design variable

        aircraft.turbofan_engine.core_thrust_ratio = init.core_thrust_ratio()
        aircraft.turbofan_engine.core_width_ratio = init.core_width_ratio()
        aircraft.turbofan_engine.core_weight_ratio = init.core_weight_ratio()

    if (aircraft.propulsion.architecture==3):
        aircraft.body_nacelle.length = init.nacelle_body_length()
        aircraft.body_nacelle.width = init.nacelle_body_width()
        aircraft.body_nacelle.hub_width = init.nacelle_body_hub_width()

    if (aircraft.propulsion.architecture==4):
        aircraft.turboprop_engine.n_engine = n_engine
        aircraft.turboprop_engine.reference_thrust = init.reference_thrust(n_pax_ref,design_range,n_engine)                                            # Main design variable
        aircraft.turboprop_engine.propeller_efficiency = init.propeller_efficiency()

    aircraft.power_elec_chain.mto = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.mcn = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.mcl = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.mcr = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.fid = 0.01

    aircraft.battery.strategy = init.battery_strategy()
    aircraft.battery.power_feed = init.battery_power_feed()
    aircraft.battery.time_feed = init.battery_time_feed()
    aircraft.battery.energy_cruise = init.battery_energy_cruise()
    aircraft.battery.energy_density = init.battery_energy_density()
    aircraft.battery.power_density = init.battery_power_density()

    aircraft.power_elec_chain.overall_efficiency = init.e_chain_efficiency()
    aircraft.power_elec_chain.generator_pw_density = init.generator_power_density()
    aircraft.power_elec_chain.rectifier_pw_density = init.rectifier_pw_density()
    aircraft.power_elec_chain.wiring_pw_density = init.wiring_pw_density()
    aircraft.power_elec_chain.cooling_pw_density = init.cooling_pw_density()

    aircraft.electric_nacelle.efficiency_fan = init.efficiency_fan()
    aircraft.electric_nacelle.efficiency_prop = init.efficiency_prop()
    aircraft.electric_nacelle.motor_efficiency = init.e_motor_efficiency()
    aircraft.electric_nacelle.controler_efficiency = init.controler_efficiency()
    aircraft.electric_nacelle.controler_pw_density = init.controler_pw_density()
    aircraft.electric_nacelle.motor_pw_density = init.e_motor_pw_density()
    aircraft.electric_nacelle.nacelle_pw_density = init.e_nacelle_pw_density()

    aircraft.propulsion.bli_effect = init.boundary_layer_effect()

    return


#===========================================================================================================
def aircraft_pre_design(aircraft):
    """
    Perform geometrical pre design
    Solves the coupling carried by nacelle geometry
    """

    airframe.eval_cabin_design(aircraft)
    airframe.eval_fuselage_design(aircraft)

    #===========================================================================================================
    def fct_aircraft_pre_design(x_in,aircraft):

        ac = aircraft

        ac.turbofan_nacelle.width = x_in[0]                         # Coupling variable
        ac.turbofan_nacelle.y_ext = x_in[1]                         # Coupling variable

        airframe.eval_pre_design_vtp(ac)
        airframe.eval_pre_design_wing(ac)
        airframe.eval_pre_design_htp(ac)

        propulsion.eval_propulsion_design(ac)

        y_out = np.array([x_in[0] - ac.turbofan_nacelle.width, \
                          x_in[1] - ac.turbofan_nacelle.y_ext])

        return y_out
    #-----------------------------------------------------------------------------------------------------------

    bpr = aircraft.turbofan_engine.bpr
    reference_thrust = aircraft.turbofan_engine.reference_thrust
    nacelle_attachment = aircraft.turbofan_nacelle.attachment
    fuselage_width = aircraft.fuselage.width

    nacelle_width_i = init.turbofan_nacelle_width(bpr,reference_thrust)
    nacelle_y_ext_i = init.turbofan_nacelle_y_ext(nacelle_attachment,fuselage_width,nacelle_width_i)

    x_ini = np.array([nacelle_width_i,nacelle_y_ext_i])

    fct_arg = aircraft

    output_dict = fsolve(fct_aircraft_pre_design, x0=x_ini, args=fct_arg, full_output=True)

    aircraft.turbofan_nacelle.width = output_dict[0][0]                         # Coupling variable
    aircraft.turbofan_nacelle.y_ext = output_dict[0][1]                         # Coupling variable

    airframe.eval_pre_design_vtp(aircraft)
    airframe.eval_pre_design_wing(aircraft)
    airframe.eval_pre_design_htp(aircraft)

    propulsion.eval_propulsion_design(aircraft)

    airplane.eval_aerodynamics_design(aircraft)

    return


#===========================================================================================================
def mass_estimation(aircraft):
    """
    Estimate mass and CGs of the airplane
    Takes MTOW as input but solves the coupling carried by MZFW and MLW
    """

    #===========================================================================================================
    def fct_mass(x_in,aircraft):

        ac = aircraft

        ac.weights.mlw = x_in[0]                         # Coupling variable
        ac.weights.mzfw = x_in[1]                         # Coupling variable

        # Mass
        #------------------------------------------------------------------------------------------------------
        airframe.eval_cabin_mass(ac)
        airframe.eval_fuselage_mass(ac)
        airframe.eval_vtp_mass(ac)
        airframe.eval_wing_mass(ac)
        airframe.eval_htp_mass(ac)
        airframe.eval_landing_gear_mass(ac)

        propulsion.eval_propulsion_mass(ac)
        propulsion.eval_battery_mass(ac)
        propulsion.eval_tank_data(ac)

        airplane.eval_system_mass(ac)
        airplane.eval_payload_mass(ac)
        airplane.eval_aircraft_weights(ac)

        y_out = np.array([x_in[0] - ac.weights.mlw, \
                          x_in[1] - ac.weights.mzfw])

        return y_out
    #-----------------------------------------------------------------------------------------------------------

    n_pax_ref = aircraft.cabin.n_pax_ref
    design_range = aircraft.design_driver.design_range

    mtow_i = aircraft.weights.mtow
    mzfw_i = aircraft.weights.mzfw
    mlw_i = aircraft.weights.mlw

    x_ini = np.array([mlw_i,mzfw_i])

    fct_arg = aircraft

    output_dict = fsolve(fct_mass, x0=x_ini, args=fct_arg, full_output=True)

    aircraft.weights.mlw = output_dict[0][0]                          # Coupling variable
    aircraft.weights.mzfw = output_dict[0][1]                         # Coupling variable

    # Update mass
    #------------------------------------------------------------------------------------------------------
    airframe.eval_cabin_mass(aircraft)
    airframe.eval_fuselage_mass(aircraft)
    airframe.eval_vtp_mass(aircraft)
    airframe.eval_wing_mass(aircraft)
    airframe.eval_htp_mass(aircraft)
    airframe.eval_landing_gear_mass(aircraft)

    propulsion.eval_propulsion_mass(aircraft)
    propulsion.eval_battery_mass(aircraft)
    propulsion.eval_tank_data(aircraft)

    airplane.eval_system_mass(aircraft)
    airplane.eval_payload_mass(aircraft)
    airplane.eval_aircraft_weights(aircraft)
    airplane.eval_aircraft_cg(aircraft)

    return


#===========================================================================================================
def nominal_mission(aircraft):
    """
    Compute nominal mission with range as input
    """

    disa = 0
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach
    nei = 0

    aircraft.nominal_mission.payload = aircraft.payload.nominal
    aircraft.nominal_mission.range = aircraft.design_driver.design_range
    aircraft.nominal_mission.tow = aircraft.weights.mtow

    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code
    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    lod_max,_ = craft_aero.lod_max(aircraft, pamb, tamb, mach)
    sfc = propu.sfc(aircraft,pamb,tamb,mach,MCR,nei)

    aircraft.high_speed.cruise_lod = lod_max
    aircraft.high_speed.cruise_sfc = sfc

    payload = aircraft.nominal_mission.payload
    range = aircraft.nominal_mission.range
    tow = aircraft.nominal_mission.tow

    block_fuel,block_time,total_fuel = perfo.mission(aircraft,range,tow,altp,mach,disa)

    aircraft.nominal_mission.block_fuel = block_fuel
    aircraft.nominal_mission.block_time = block_time
    aircraft.nominal_mission.total_fuel = total_fuel

    return


#===========================================================================================================
def mass_mission_adaptation(aircraft):
    """
    Perform mass - mission adaptation and update mass and CGs
    """

    #===========================================================================================================
    def fct_mass_mission(x_in,aircraft):

        ac = aircraft

        ac.weights.mtow = x_in[0]                         # Coupling variable
        ac.weights.mlw = x_in[1]                         # Coupling variable
        ac.weights.mzfw = x_in[2]                         # Coupling variable

        # Mass
        #------------------------------------------------------------------------------------------------------
        airframe.eval_cabin_mass(ac)
        airframe.eval_fuselage_mass(ac)
        airframe.eval_vtp_mass(ac)
        airframe.eval_wing_mass(ac)
        airframe.eval_htp_mass(ac)
        airframe.eval_landing_gear_mass(ac)

        propulsion.eval_propulsion_mass(ac)
        propulsion.eval_battery_mass(ac)
        propulsion.eval_tank_data(ac)

        airplane.eval_system_mass(ac)
        airplane.eval_payload_mass(ac)
        airplane.eval_aircraft_weights(ac)

        # Mission
        #------------------------------------------------------------------------------------------------------
        disa = 0
        altp = ac.design_driver.ref_cruise_altp
        mach = ac.design_driver.cruise_mach

        ac.nominal_mission.payload = ac.payload.nominal
        ac.nominal_mission.range = ac.design_driver.design_range
        ac.nominal_mission.tow = ac.weights.mtow

        payload = ac.nominal_mission.payload
        range = ac.nominal_mission.range
        tow = ac.nominal_mission.tow

        block_fuel,block_time,total_fuel = perfo.mission(ac,range,tow,altp,mach,disa)

        ac.nominal_mission.block_fuel = block_fuel
        ac.nominal_mission.block_time = block_time
        ac.nominal_mission.total_fuel = total_fuel

        required_mtow = ac.weights.owe + payload + total_fuel

        y_out = np.array([x_in[0] - required_mtow, \
                          x_in[1] - ac.weights.mlw, \
                          x_in[2] - ac.weights.mzfw])

        return y_out
    #-----------------------------------------------------------------------------------------------------------

    n_pax_ref = aircraft.cabin.n_pax_ref
    design_range = aircraft.design_driver.design_range

    mtow_i = init.mtow(n_pax_ref,design_range)
    mzfw_i = init.mzfw(n_pax_ref,design_range)
    mlw_i = init.mlw(mtow_i,mzfw_i)

    x_ini = np.array([mtow_i,mlw_i,mzfw_i])

    fct_arg = aircraft

    output_dict = fsolve(fct_mass_mission, x0=x_ini, args=fct_arg, full_output=True)

    aircraft.weights.mtow = output_dict[0][0]                         # Coupling variable
    aircraft.weights.mlw = output_dict[0][1]                          # Coupling variable
    aircraft.weights.mzfw = output_dict[0][2]                         # Coupling variable

    # Update mass
    #------------------------------------------------------------------------------------------------------
    airframe.eval_cabin_mass(aircraft)
    airframe.eval_fuselage_mass(aircraft)
    airframe.eval_vtp_mass(aircraft)
    airframe.eval_wing_mass(aircraft)
    airframe.eval_htp_mass(aircraft)
    airframe.eval_landing_gear_mass(aircraft)

    propulsion.eval_propulsion_mass(aircraft)
    propulsion.eval_battery_mass(aircraft)
    propulsion.eval_tank_data(aircraft)

    airplane.eval_system_mass(aircraft)
    airplane.eval_payload_mass(aircraft)
    airplane.eval_aircraft_weights(aircraft)
    airplane.eval_aircraft_cg(aircraft)

    return


#===========================================================================================================
def payload_range_analysis(aircraft):
    """
    Compute Payload - Range diagram corner points
    """
    disa = 0
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach

    # Max payload mission
    #------------------------------------------------------------------------------------------------------
    tow = aircraft.weights.mtow
    payload = aircraft.payload.nominal

    aircraft.max_payload_mission.tow = tow
    aircraft.max_payload_mission.payload = payload

    range,block_fuel,block_time,total_fuel = sub_proc.mission_range(aircraft,tow,payload,altp,mach,disa)

    aircraft.max_payload_mission.range = range
    aircraft.max_payload_mission.block_fuel = block_fuel
    aircraft.max_payload_mission.block_time = block_time
    aircraft.max_payload_mission.total_fuel = total_fuel

    # Max fuel mission
    #------------------------------------------------------------------------------------------------------
    tow = aircraft.weights.mtow
    total_fuel = aircraft.weights.mfw

    aircraft.max_fuel_mission.tow = tow
    aircraft.max_fuel_mission.total_fuel = total_fuel

    range,payload,block_fuel,block_time = sub_proc.mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa)

    aircraft.max_fuel_mission.payload = payload
    aircraft.max_fuel_mission.range = range
    aircraft.max_fuel_mission.block_fuel = block_fuel
    aircraft.max_fuel_mission.block_time = block_time

    # zero fuel mission
    #------------------------------------------------------------------------------------------------------
    total_fuel = aircraft.weights.mfw
    tow = aircraft.weights.owe + total_fuel

    aircraft.zero_payload_mission.tow = tow
    aircraft.zero_payload_mission.total_fuel = total_fuel

    range,payload,block_fuel,block_time = sub_proc.mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa)

    aircraft.zero_payload_mission.range = range
    aircraft.zero_payload_mission.block_fuel = block_fuel
    aircraft.zero_payload_mission.block_time = block_time

    return


#===========================================================================================================
def performance_analysis(aircraft):
    """
    Compute operational performances
    """

    # Nominal mission (here : range is an output)
    #------------------------------------------------------------------------------------------------------
    disa = 0.
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach
    nei = 0.

    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code
    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    lod_max,_ = craft_aero.lod_max(aircraft, pamb, tamb, mach)
    sfc = propu.sfc(aircraft,pamb,tamb,mach,MCR,nei)

    aircraft.high_speed.cruise_lod = lod_max
    aircraft.high_speed.cruise_sfc = sfc

    aircraft.nominal_mission.payload = aircraft.payload.nominal
    aircraft.nominal_mission.range = aircraft.design_driver.design_range
    aircraft.nominal_mission.tow = aircraft.weights.mtow

    payload = aircraft.nominal_mission.payload
    tow = aircraft.nominal_mission.tow

    range,block_fuel,block_time,total_fuel = sub_proc.mission_range(aircraft,tow,payload,altp,mach,disa)

    aircraft.nominal_mission.range = range

    aircraft.nominal_mission.block_fuel = block_fuel
    aircraft.nominal_mission.block_time = block_time
    aircraft.nominal_mission.total_fuel = total_fuel

    # Ceilings
    #------------------------------------------------------------------------------------------------------
    toc = aircraft.design_driver.top_of_climb_altp
    oei_ceil_req = aircraft.low_speed.req_oei_altp

    vz_clb,vz_crz,oei_path,oei_mach = perfo.ceilings(aircraft,toc,oei_ceil_req)

    aircraft.high_speed.eff_vz_climb = vz_clb
    aircraft.high_speed.eff_vz_cruise = vz_crz
    aircraft.low_speed.eff_oei_path = oei_path

    # Time to climb to requested altitude
    #------------------------------------------------------------------------------------------------------
    toc = aircraft.high_speed.req_toc_altp
    disa = 0
    mass = aircraft.weights.mtow
    vcas1 = aircraft.high_speed.cas1_ttc
    vcas2 = aircraft.high_speed.cas2_ttc
    mach = aircraft.design_driver.cruise_mach

    ttc = perfo.time_to_climb(aircraft,toc,disa,mass,vcas1,vcas2,mach)

    aircraft.high_speed.eff_ttc = ttc

    # Take off field length
    #------------------------------------------------------------------------------------------------------
    altp = aircraft.low_speed.altp_tofl
    disa = aircraft.low_speed.disa_tofl
    mass = aircraft.weights.mtow
    hld_conf_to = aircraft.aerodynamics.hld_conf_to

    tofl,seg2_path,eff_kvs1g,limitation = perfo.take_off_field_length(aircraft,altp,disa,mass,hld_conf_to)

    aircraft.low_speed.eff_tofl = tofl
    aircraft.low_speed.eff_kvs1g = eff_kvs1g
    aircraft.low_speed.seg2_path = seg2_path
    aircraft.low_speed.limitation = limitation

    # Approach speed
    #------------------------------------------------------------------------------------------------------
    altp = aircraft.low_speed.altp_app_speed
    disa = aircraft.low_speed.disa_app_speed
    mass = aircraft.weights.mlw
    hld_conf_ld = aircraft.aerodynamics.hld_conf_to

    app_speed = perfo.approach_speed(aircraft,altp,disa,mass,hld_conf_ld)

    aircraft.low_speed.eff_app_speed = app_speed

    # Environment
    #------------------------------------------------------------------------------------------------------
    CO2_metric,rgf = environ.fuel_efficiency_metric(aircraft)

    aircraft.environmental_impact.rgf = rgf
    aircraft.environmental_impact.CO2_metric = CO2_metric

    # Cost mission
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach

    disa = aircraft.cost_mission.disa
    range = aircraft.cost_mission.range

    payload = aircraft.payload.nominal

    aircraft.cost_mission.payload = payload

    tow,block_fuel,block_time,total_fuel = sub_proc.mission_tow(aircraft,payload,range,altp,mach,disa)

    aircraft.cost_mission.block_fuel = block_fuel
    aircraft.cost_mission.block_time = block_time
    aircraft.cost_mission.total_fuel = total_fuel

    aircraft.cost_mission.block_CO2 = block_fuel * aircraft.environmental_impact.CO2_index

    # Economics
    #------------------------------------------------------------------------------------------------------
    direct_op_cost,cash_op_cost,block_fuel,engine_price,airplane_price = costing.operating_costs(aircraft,block_fuel,block_time)

    aircraft.economics.engine_price = engine_price
    aircraft.economics.airplane_price = airplane_price

    aircraft.economics.direct_operating_cost = direct_op_cost
    aircraft.economics.cash_operating_cost = cash_op_cost

    return


#===========================================================================================================
def handling_qualities_analysis(aircraft):
    """
    Compute CG limits from handling qualities point of view
    """
    # Forward limit : trim landing
    #------------------------------------------------------------------------------------------------------
    altp = unit.m_ft(0)
    disa = 0
    nei = 0
    speed_mode = 1
    hld_conf = aircraft.aerodynamics.hld_conf_ld
    mass = aircraft.center_of_gravity.max_fwd_mass

    cg_max_fwd_stall,speed,fn,aoa,ih,c_z,cx_trimmed = h_q.forward_cg_stall(aircraft,altp,disa,nei,hld_conf,speed_mode,mass)

    aircraft.center_of_gravity.max_fwd_trim_cg = cg_max_fwd_stall         # Forward cg limit

    # Backward limit : static stability
    #------------------------------------------------------------------------------------------------------
    stability_margin = 0.05

    cg_max_bwd_stab = h_q.backward_cg_stab(aircraft,stability_margin)

    aircraft.center_of_gravity.max_bwd_stab_cg = cg_max_bwd_stab          # Backward cg limit

    return


#===========================================================================================================
def eval_optim_data(x_in,aircraft,crit_index):
    """
    Compute criterion and constraints
    """

    ac = aircraft

    if (ac.propulsion.architecture==1):

        ac.turbofan_engine.reference_thrust = x_in[0]

    elif (ac.propulsion.architecture==2):

        ac.turbofan_engine.reference_thrust = x_in[0]

    elif (ac.propulsion.architecture==3):

        ac.turbofan_engine.reference_thrust = x_in[0]

    elif (ac.propulsion.architecture==4):

        ac.turboprop_engine.reference_thrust = x_in[0]

    else:
        raise Exception("propulsion.architecture index is out of range")

    ac.wing.area = x_in[1]

    aircraft_pre_design(ac)
    mass_mission_adaptation(ac)
    performance_analysis(ac)

    cst = np.zeros(6)

    # Constraints are violated is negative
    #------------------------------------------------------------------------------------------------------
    cst[0] = (ac.high_speed.eff_vz_climb - ac.high_speed.req_vz_climb) * 1.
    cst[1] = (ac.high_speed.eff_vz_cruise - ac.high_speed.req_vz_cruise) * 1.
    cst[2] = (ac.low_speed.eff_oei_path - ac.low_speed.req_oei_path) * 1.

    cst[3] = (ac.high_speed.req_ttc - ac.high_speed.eff_ttc) * 1.
    cst[4] = (ac.low_speed.req_tofl - ac.low_speed.eff_tofl) * 1.
    cst[5] = (ac.low_speed.req_app_speed - ac.low_speed.eff_app_speed) * 1.

    crt = np.zeros(5)

    # All criteria have to be minimized
    #------------------------------------------------------------------------------------------------------
    crt[0] = ac.weights.mtow
    crt[1] = ac.cost_mission.block_fuel
    crt[2] = ac.environmental_impact.CO2_metric
    crt[3] = ac.economics.cash_operating_cost
    crt[4] = ac.economics.direct_operating_cost

    crit = crt[crit_index]

    return crit,cst


#===========================================================================================================
def eval_optim_cst(x_in,aircraft,crit_index):
    """
    Retrieve constraints
    """

    crit,cst = eval_optim_data(x_in,aircraft,crit_index)

    print("cst :",cst)

    return cst


#===========================================================================================================
def eval_optim_crt(x_in,aircraft,crit_index):
    """
    Retreve criteria
    """

    crit,cst = eval_optim_data(x_in,aircraft,crit_index)

    print("Design :",x_in)
    print("Crit :",crit)

    return crit


#===========================================================================================================
def optimization(aircraft,search_domain,criterion):
    """
    Compute criterion and constraints
    """

    if (aircraft.propulsion.architecture==1):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==2):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==3):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==4):
        start_value = (aircraft.turboprop_engine.reference_thrust,aircraft.wing.area)
    else:
        raise Exception("propulsion.architecture index is out of range")


    if (criterion=="MTOW"):
        crit_index = 0
    elif (criterion=="cost_fuel"):
        crit_index = 1
    elif (criterion=="CO2_metric"):
        crit_index = 2
    elif (criterion=="COC"):
        crit_index = 3
    elif (criterion=="DOC"):
        crit_index = 4
    else:
        raise Exception("Criterion name is unknown")


    res = minimize(eval_optim_crt, start_value, args=(aircraft,crit_index,), method="SLSQP", bounds=search_domain,
                   constraints={"type":"ineq","fun":eval_optim_cst,"args":(aircraft,crit_index,)},
                   options={"maxiter":20,"ftol":0.025})

    #res = minimize(eval_optim_crt, x_in, args=(aircraft,), method="COBYLA", bounds=((110000,140000),(120,160)),
    #               constraints={"type":"ineq","fun":eval_optim_cst,"args":(aircraft,)},
    #               options={"maxiter":100,"tol":0.1,"catol":0.0002,'rhobeg': 1.0})
    print(res)

    return res


