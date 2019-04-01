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
from MARILib.aircraft_model.airplane import airplane_design as airplane, \
                                            aerodynamics as craft_aero, \
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
def aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config):
    """
    Initialize a generic aircraft
    """

    aircraft.propulsion.architecture = propu_config

    n_engine = init.n_engine()

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

    aircraft.high_speed.req_toc_altp = init.top_of_climb_altp(propu_config)
    aircraft.high_speed.disa_climb = init.disa_climb()
    aircraft.high_speed.req_vz_climb = init.req_vz_climb()           # TLR
    aircraft.high_speed.req_vz_cruise = init.req_vz_cruise()         # TLR
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

    aircraft.turbofan_nacelle.attachment = init.nacelle_attachment()
    aircraft.turbofan_nacelle.efficiency_fan = init.efficiency_fan()
    aircraft.turbofan_nacelle.efficiency_prop = init.efficiency_prop()

    aircraft.turbofan_engine.n_engine = n_engine
    aircraft.turbofan_engine.bpr = init.bpr()
    aircraft.turbofan_engine.reference_thrust = init.reference_thrust(n_pax_ref,design_range,n_engine)                                            # Main design variable

    aircraft.turbofan_engine.core_thrust_ratio = init.core_thrust_ratio()
    aircraft.turbofan_engine.core_width_ratio = init.core_width_ratio()
    aircraft.turbofan_engine.core_weight_ratio = init.core_weight_ratio()

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

    aircraft.electric_nacelle.hub_width = init.e_hub_width()
    aircraft.electric_nacelle.efficiency_fan = init.efficiency_fan()
    aircraft.electric_nacelle.efficiency_prop = init.efficiency_prop()
    aircraft.electric_nacelle.motor_efficiency = init.e_motor_efficiency()
    aircraft.electric_nacelle.controler_efficiency = init.controler_efficiency()
    aircraft.electric_nacelle.controler_pw_density = init.controler_pw_density()
    aircraft.electric_nacelle.motor_pw_density = init.e_motor_pw_density()
    aircraft.electric_nacelle.nacelle_pw_density = init.e_nacelle_pw_density()

    aircraft.propulsion.bli_effect = init.boundary_layer_effect()

    return


