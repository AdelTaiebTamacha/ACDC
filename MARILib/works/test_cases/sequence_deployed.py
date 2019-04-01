#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""


import pickle

from MARILib.tools import units as unit

from MARILib.aircraft_data.aircraft_description import Aircraft

from MARILib.earth import environment as earth

from MARILib.airplane.airframe import airframe_design as airframe

from MARILib.airplane.propulsion import propulsion_design as propulsion
from MARILib.airplane.propulsion import propulsion_models as propu

from MARILib.aircraft_model.operations import handling_qualities as h_q,\
                                              mission as perfo, \
                                              environmental_impact as environ, \
                                              pricing_and_costing as costing

from MARILib.aircraft_model.airplane import regulation as regul, \
    airplane_design as airplane, \
                                            aerodynamics as craft_aero, viewer as show
from MARILib.aircraft_model import initialization as init

from MARILib.processes import solvers as sub_proc


# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
ac = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propulsive_architecture = init.prop_architecture()
ac.propulsion.architecture = propulsive_architecture

ac.propulsion.fuel_type = init.fuel_type()

n_engine = init.n_engine()

ac.name = "my_test_airplane"
ac.design_driver.design_range = design_range        # TLR
ac.design_driver.cruise_mach = cruise_mach          # TLR
ac.cabin.n_pax_ref = n_pax_ref                      # TLR


ac.design_driver.ref_cruise_altp = init.ref_cruise_altp(propulsive_architecture)        # TLR
ac.design_driver.top_of_climb_altp = init.top_of_climb_altp(propulsive_architecture)    # TLR

ac.aerodynamics.hld_conf_clean = init.hld_conf_clean()
ac.aerodynamics.hld_conf_ld = init.hld_conf_ld()

ac.low_speed.altp_tofl = init.altp_tofl()
ac.low_speed.disa_tofl = init.disa_tofl()
ac.low_speed.kvs1g_tofl = regul.kvs1g_min_take_off()       # Regulation
ac.low_speed.req_tofl = init.req_tofl(design_range)        # TLR

ac.low_speed.altp_app_speed = init.altp_app_speed()
ac.low_speed.disa_app_speed = init.disa_app_speed()
ac.low_speed.kvs1g_app_speed = regul.kvs1g_min_landing()   # Regulation
ac.low_speed.req_app_speed = init.req_app_speed()          # TLR

ac.low_speed.disa_oei = init.disa_oei()
ac.low_speed.req_oei_path = regul.ceil_oei_min_path(n_engine)  # Regulation
ac.low_speed.req_oei_altp = init.req_oei_altp(propulsive_architecture)                # TLR

ac.high_speed.disa_climb = init.disa_climb()
ac.high_speed.req_vz_climb = init.req_vz_climb()           # TLR
ac.high_speed.req_vz_cruise = init.req_vz_cruise()         # TLR
ac.high_speed.req_toc_altp = init.top_of_climb_altp(propulsive_architecture)
ac.high_speed.cas1_ttc = init.cas1_ttc(propulsive_architecture)
ac.high_speed.cas2_ttc = init.cas2_ttc(propulsive_architecture)
ac.high_speed.req_ttc = init.req_ttc()                     # TLR

ac.cost_mission.disa = init.cost_mission_disa()
ac.cost_mission.range = init.cost_mission_range(design_range)

ac.economics.fuel_price = init.fuel_price()
ac.economics.elec_price = init.elec_price()
ac.economics.battery_price = init.battery_price()
ac.economics.labor_cost = init.labor_cost()
ac.economics.irp = init.irp()
ac.economics.period = init.period()
ac.economics.interest_rate = init.interest_rate()
ac.economics.utilisation = init.utilisation(design_range)

ac.environmental_impact.CO2_index = earth.emission_index("CO2")
ac.environmental_impact.H2O_index = earth.emission_index("H2O")
ac.environmental_impact.SO2_index = earth.emission_index("SO2")
ac.environmental_impact.NOx_index = earth.emission_index("NOx")
ac.environmental_impact.CO_index = earth.emission_index("CO")
ac.environmental_impact.HC_index = earth.emission_index("HC")
ac.environmental_impact.sulfuric_acid_index = earth.emission_index("sulfuric_acid")
ac.environmental_impact.nitrous_acid_index = earth.emission_index("nitrous_acid")
ac.environmental_impact.nitric_acid_index = earth.emission_index("nitric_acid")
ac.environmental_impact.soot_index = earth.emission_index("soot")

mzfw_i = init.mzfw(n_pax_ref,design_range)
mtow_i = init.mtow(n_pax_ref,design_range)
mlw_i = init.mlw(mtow_i,mzfw_i)

ac.weights.mzfw = mzfw_i                            # Coupling variable of mass perfo adaptation
ac.weights.mlw = mlw_i                              # Coupling variable of mass perfo adaptation
ac.weights.mtow = mtow_i                            # Coupling variable of mass perfo adaptation

n_pax_front_i = init.n_pax_front(n_pax_ref)
ac.cabin.n_aisle = init.n_aisle(n_pax_front_i)         # TLR
ac.cabin.n_pax_front = n_pax_front_i                    # TLR

ac.payload.m_pax_nominal = init.m_pax_nominal(design_range)
ac.payload.m_pax_max = init.m_pax_max(design_range)

ac.wing.attachment = init.wing_attachment()
ac.wing.morphing = init.wing_morphing()
ac.wing.hld_type = init.hld_type()

wing_area_i = init.wing_area(n_pax_ref,design_range)
wing_aspect_ratio_i = init.wing_aspect_ratio()

ac.wing.area = wing_area_i                                              # Main design variable
ac.wing.aspect_ratio = wing_aspect_ratio_i
ac.wing.span = init.wing_span(wing_area_i,wing_aspect_ratio_i)

ac.horizontal_tail.attachment = init.htp_attachment()

nacelle_attachment_i = init.nacelle_attachment()
ac.turbofan_nacelle.attachment = nacelle_attachment_i
ac.turbofan_nacelle.efficiency_fan = init.efficiency_fan()
ac.turbofan_nacelle.efficiency_prop = init.efficiency_prop()

bpr_i = init.bpr()
reference_thrust_i = init.reference_thrust(n_pax_ref,design_range,n_engine)

ac.turbofan_engine.n_engine = n_engine
ac.turbofan_engine.bpr = bpr_i
ac.turbofan_engine.reference_thrust = reference_thrust_i                                            # Main design variable
ac.turbofan_engine.core_thrust_ratio = init.core_thrust_ratio()
ac.turbofan_engine.core_width_ratio = init.core_width_ratio()
ac.turbofan_engine.core_weight_ratio = init.core_weight_ratio()

ac.body_nacelle.length = init.nacelle_body_length()
ac.body_nacelle.width = init.nacelle_body_width()
ac.body_nacelle.hub_width = init.nacelle_body_hub_width()

ac.turboprop_engine.n_engine = n_engine
ac.turboprop_engine.reference_thrust = reference_thrust_i                                            # Main design variable
ac.turboprop_engine.propeller_efficiency = init.propeller_efficiency()

e_power = 1e6       # Watts, electric motor power

ac.power_elec_chain.mto = e_power
ac.power_elec_chain.mcn = e_power
ac.power_elec_chain.mcl = e_power
ac.power_elec_chain.mcr = e_power
ac.power_elec_chain.fid = 0.01

ac.battery.strategy = init.battery_strategy()
ac.battery.power_feed = init.battery_power_feed()
ac.battery.time_feed = init.battery_time_feed()
ac.battery.energy_cruise = init.battery_energy_cruise()
ac.battery.energy_density = init.battery_energy_density()
ac.battery.power_density = init.battery_power_density()

ac.power_elec_chain.overall_efficiency = init.e_chain_efficiency()
ac.power_elec_chain.generator_pw_density = init.generator_power_density()
ac.power_elec_chain.rectifier_pw_density = init.rectifier_pw_density()
ac.power_elec_chain.wiring_pw_density = init.wiring_pw_density()
ac.power_elec_chain.cooling_pw_density = init.cooling_pw_density()

ac.electric_nacelle.efficiency_fan = init.efficiency_fan()
ac.electric_nacelle.efficiency_prop = init.efficiency_prop()
ac.electric_nacelle.motor_efficiency = init.e_motor_efficiency()
ac.electric_nacelle.controler_efficiency = init.controler_efficiency()
ac.electric_nacelle.controler_pw_density = init.controler_pw_density()
ac.electric_nacelle.motor_pw_density = init.e_motor_pw_density()
ac.electric_nacelle.nacelle_pw_density = init.e_nacelle_pw_density()

ac.propulsion.bli_effect = init.boundary_layer_effect()

print("-------------------------------------------")
print("Initialization : done")

# airplane geometry pre design
#------------------------------------------------------------------------------------------------------
airframe.eval_cabin_design(ac)
airframe.eval_fuselage_design(ac)

nacelle_width_i = init.turbofan_nacelle_width(bpr_i,reference_thrust_i)
ac.turbofan_nacelle.width = nacelle_width_i                         # Coupling variable

fuselage_width_i = ac.fuselage.width
nacelle_y_ext_i = init.turbofan_nacelle_y_ext(nacelle_attachment_i,fuselage_width_i,nacelle_width_i)
ac.turbofan_nacelle.y_ext = nacelle_y_ext_i                         # Coupling variable

airframe.eval_pre_design_vtp(ac)
airframe.eval_pre_design_wing(ac)
airframe.eval_pre_design_htp(ac)

propulsion.eval_propulsion_design(ac)

delta_nacelle_width = ac.turbofan_nacelle.width - nacelle_width_i   # Should be driven to zero
delta_nacelle_y_ext = ac.turbofan_nacelle.y_ext - nacelle_y_ext_i   # Should be driven to zero

print('Coupling : delta_nacelle_width = ',"%.3f"%delta_nacelle_width,' (= 0 ?)')
print('Coupling : delta_nacelle_y_ext = ',"%.3f"%delta_nacelle_y_ext,' (= 0 ?)')

airplane.eval_aerodynamics_design(ac)

print("-------------------------------------------")
print("Pre design : done")

# airplane mass & CG estimation
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
airplane.eval_aircraft_cg(ac)

delta_mzfw = mzfw_i - ac.weights.mzfw                               # Should be positive or zero
delta_mlw = mlw_i - ac.weights.mlw                                  # Should be positive or zero

print('Coupling : delta_mzfw = ',"%.2f"%delta_mzfw,' kg, (= 0 ?)')
print('Coupling : delta_mlw = ',"%.2f"%delta_mlw,'kg, (= 0 ?)')

print("-------------------------------------------")
print("Mass and CG estimation : done")

# Nominal mission
#------------------------------------------------------------------------------------------------------
disa = 0
altp = ac.design_driver.ref_cruise_altp
mach = ac.design_driver.cruise_mach
nei = 0

(MTO,MCN,MCL,MCR,FID) = ac.propulsion.rating_code
pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
lod_max,_ = craft_aero.lod_max(ac, pamb, tamb, mach)
sfc = propu.sfc(ac,pamb,tamb,mach,MCR,nei)

ac.high_speed.cruise_lod = lod_max
ac.high_speed.cruise_sfc = sfc

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

mtow_eff = ac.weights.owe + payload + total_fuel

delta_mtow = mtow_i - mtow_eff                                  # Should be positive or zero

print("-------------------------------------------")
print("Nominal mission : done")
print("")
print('Coupling : delta_mtow = ',"%.2f"%delta_mtow,'kg, (= 0 ?)')

# Ceilings
#------------------------------------------------------------------------------------------------------
toc = ac.design_driver.top_of_climb_altp
oei_ceil_req = ac.low_speed.req_oei_altp

vz_clb,vz_crz,oei_path,oei_mach = perfo.ceilings(ac,toc,oei_ceil_req)

ac.high_speed.eff_vz_climb = vz_clb
ac.high_speed.eff_vz_cruise = vz_crz
ac.low_speed.eff_oei_path = oei_path
ac.low_speed.eff_oei_mach = oei_mach

delta_vz_clb = ac.high_speed.req_vz_climb - vz_clb
delta_vz_crz = ac.high_speed.req_vz_cruise - vz_crz
delta_oei_path = ac.low_speed.req_oei_path - oei_path

print('')
print('Constraint : delta_vz_clb = ',"%.2f"%unit.ftpmin_mps(delta_vz_clb),'ft/min, (=< 0 ?)')
print('Constraint : delta_vz_crz = ',"%.2f"%unit.ftpmin_mps(delta_vz_crz),'ft/min, (=< 0 ?)')
print('Constraint : delta_oei_path = ',"%.2f"%(delta_oei_path*100),'%, (=< 0 ?)')

# Time to climb to requested altitude
#------------------------------------------------------------------------------------------------------
toc = ac.high_speed.req_toc_altp
disa = 0
mass = ac.weights.mtow
vcas1 = ac.high_speed.cas1_ttc
vcas2 = ac.high_speed.cas2_ttc
mach = ac.design_driver.cruise_mach

ttc = perfo.time_to_climb(ac,toc,disa,mass,vcas1,vcas2,mach)

ac.high_speed.eff_ttc = ttc

delta_ttc = ttc - ac.high_speed.req_ttc

print('Constraint : delta_ttc = ',"%.2f"%unit.min_s(delta_ttc),'min (=< 0 ?)')

# Take off field length
#------------------------------------------------------------------------------------------------------
altp = ac.low_speed.altp_tofl
disa = ac.low_speed.disa_tofl
mass = ac.weights.mtow
hld_conf_to = ac.aerodynamics.hld_conf_to

tofl,seg2_path,eff_kvs1g,limitation = perfo.take_off_field_length(ac,altp,disa,mass,hld_conf_to)

ac.low_speed.eff_tofl = tofl
ac.low_speed.eff_kvs1g = eff_kvs1g
ac.low_speed.seg2_path = seg2_path
ac.low_speed.limitation = limitation

delta_tofl = tofl - ac.low_speed.req_tofl

print('Constraint : delta_tofl = ',"%.2f"%delta_tofl,'m, (=< 0 ?)')

# Approach speed
#------------------------------------------------------------------------------------------------------
altp = ac.low_speed.altp_app_speed
disa = ac.low_speed.disa_app_speed
mass = ac.weights.mlw
hld_conf_ld = ac.aerodynamics.hld_conf_to

app_speed = perfo.approach_speed(ac,altp,disa,mass,hld_conf_ld)

ac.low_speed.eff_app_speed = app_speed

delta_app_speed = app_speed - ac.low_speed.req_app_speed

print('Constraint : delta_app_speed = ',"%.2f"%unit.kt_mps(delta_app_speed),'kt, (=< 0 ?)')

# Carbon dioxide metric
#------------------------------------------------------------------------------------------------------
CO2_metric,rgf = environ.fuel_efficiency_metric(ac)

ac.environmental_impact.rgf = rgf
ac.environmental_impact.CO2_metric = CO2_metric

# Cost mission
#-----------------------------------------------------------------------------------------------------------------------------------------------
altp = ac.design_driver.ref_cruise_altp
mach = ac.design_driver.cruise_mach

disa = ac.cost_mission.disa
range = ac.cost_mission.range

payload = ac.payload.nominal

ac.cost_mission.payload = payload

tow,block_fuel,block_time,total_fuel = sub_proc.mission_tow(ac,payload,range,altp,mach,disa)

ac.cost_mission.block_fuel = block_fuel
ac.cost_mission.block_time = block_time
ac.cost_mission.total_fuel = total_fuel

block_CO2 = block_fuel * ac.environmental_impact.CO2_index

ac.cost_mission.block_CO2 = block_CO2

# Economics
#------------------------------------------------------------------------------------------------------
direct_op_cost,cash_op_cost,block_fuel,engine_price,aircraft_price = costing.operating_costs(ac,block_fuel,block_time)

print('')
print('Criterion : CO2_metric = ',"%.4f"%(CO2_metric*1000),'kg/km/m0.48 (minimize)')
print('Criterion : block_CO2 = ',"%.0f"%block_CO2,'kg, (minimize)')
print('Criterion : block_fuel = ',"%.2f"%block_fuel,'kg, (minimize)')
print('Criterion : mtow_eff = ',"%.1f"%mtow_eff,'kg, (minimize)')
print('Criterion : cash_op_cost = ',"%.2f"%cash_op_cost,'$/trip, (minimize)')
print('Criterion : direct_op_cost = ',"%.2f"%direct_op_cost,'$/trip, (minimize)')

print("-------------------------------------------")
print("Performance analysis : done")

# Handling_Qualities
#------------------------------------------------------------------------------------------------------
altp = unit.m_ft(0)
disa = 0
nei = 0
speed_mode = 1
hld_conf = ac.aerodynamics.hld_conf_ld
mass = ac.center_of_gravity.max_fwd_mass

cg_max_fwd_stall,speed,fn,aoa,ih,c_z,cx_trimmed = h_q.forward_cg_stall(ac,altp,disa,nei,hld_conf,speed_mode,mass)

ac.center_of_gravity.max_fwd_trim_cg = cg_max_fwd_stall         # Forward cg limit


stability_margin = 0.05

cg_max_bwd_stab = h_q.backward_cg_stab(ac,stability_margin)

ac.center_of_gravity.max_bwd_stab_cg = cg_max_bwd_stab          # Backward cg limit


delta_fwd_cg = ac.center_of_gravity.max_fwd_trim_cg - ac.center_of_gravity.max_fwd_req_cg
delta_bwd_cg = ac.center_of_gravity.max_bwd_req_cg - ac.center_of_gravity.max_bwd_stab_cg

print('')
print('Constraint : delta_fwd_cg = ',"%.2f"%delta_fwd_cg,'m, (=<0 ?)')
print('Constraint : delta_bwd_cg = ',"%.2f"%delta_bwd_cg,'m, (=<0 ?)')


d_vtp_area = h_q.vertical_tail_sizing(ac)

print('')
print('VTP design : vtp_area = ',"%.2f"%ac.vertical_tail.area,'m2')
print('VTP design : d_vtp_area = ',"%.2f"%d_vtp_area,'m2, (=<0 ?)')

print("-------------------------------------------")
print("Handling quality analysis : done")

# Payload - Range diagram
#------------------------------------------------------------------------------------------------------
disa = 0
altp = ac.design_driver.ref_cruise_altp
mach = ac.design_driver.cruise_mach

# Max payload mission
tow = ac.weights.mtow
payload = ac.payload.nominal

ac.max_payload_mission.tow = tow
ac.max_payload_mission.payload = payload

range,block_fuel,block_time,total_fuel = sub_proc.mission_range(ac,tow,payload,altp,mach,disa)

ac.max_payload_mission.range = range
ac.max_payload_mission.block_fuel = block_fuel
ac.max_payload_mission.block_time = block_time
ac.max_payload_mission.total_fuel = total_fuel

# Max fuel mission
tow = ac.weights.mtow
total_fuel = ac.weights.mfw

ac.max_fuel_mission.tow = tow
ac.max_fuel_mission.total_fuel = total_fuel

range,payload,block_fuel,block_time = sub_proc.mission_fuel_limited(ac,tow,total_fuel,altp,mach,disa)

ac.max_fuel_mission.payload = payload
ac.max_fuel_mission.range = range
ac.max_fuel_mission.block_fuel = block_fuel
ac.max_fuel_mission.block_time = block_time

# zero fuel mission
total_fuel = ac.weights.mfw
tow = ac.weights.owe + total_fuel

ac.zero_payload_mission.tow = tow
ac.zero_payload_mission.total_fuel = total_fuel

range,payload,block_fuel,block_time = sub_proc.mission_fuel_limited(ac,tow,total_fuel,altp,mach,disa)

ac.zero_payload_mission.range = range
ac.zero_payload_mission.block_fuel = block_fuel
ac.zero_payload_mission.block_time = block_time

print("-------------------------------------------")
print("Payload-Range analysis : done")

print("-------------------------------------------")
print("Propeller =",ac.turboprop_engine.propeller_diameter)

# airplane 3D view
#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(ac,"sequence deployed","This plane")


file_name = "Aircraft.pkl"
output = open(file_name,"wb")
pickle.dump(ac,output)
output.close()
