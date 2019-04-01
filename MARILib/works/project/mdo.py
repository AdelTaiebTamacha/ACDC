#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry

This study illustrates the integration of a partial turbo-electric airplane
"""

from MARILib.tools import units as unit
from MARILib.earth import environment as earth
from MARILib.airplane.propulsion import propulsion_models as propu

from MARILib.aircraft_data.aircraft_description import Aircraft

from MARILib.processes import design as run
from MARILib.aircraft_model.airplane import viewer as show

# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
aircraft = Aircraft()

# Design drivers
#------------------------------------------------------------------------------------------------------
design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propulsion_architecture = 3 # 1:turbofan, 2:partial turboelectric first model (TD2), 3:partial turboelectric second model (project)
n_engine = 2

# Initialize all input data
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsion_architecture, n_engine)


# Study input
#======================================================================================================
study_name = "This airplane"

# Input to be adjusted according to requirements
#------------------------------------------------------------------------------------------------------
aircraft.propulsion.fuel_type = 1       # 1:kerosene, 2:hydrogene

aircraft.cabin.n_pax_front = 6
aircraft.cabin.n_aisle = 1

# Input to be optimised according to requirements and selected criterion
#------------------------------------------------------------------------------------------------------
aircraft.turbofan_engine.reference_thrust = 120000  # Newtons
aircraft.wing.area = 155                            # m2


aircraft.body_nacelle.width = 1.5       # m
aircraft.body_nacelle.length = 4.       # m


e_power = 1.0e6       # Watts, electric motor power

aircraft.power_elec_chain.mto = e_power     # required power during take off
aircraft.power_elec_chain.mcn = e_power     # required power in case of engine failure
aircraft.power_elec_chain.mcl = e_power     # required power during climb
aircraft.power_elec_chain.mcr = e_power     # required power during cruise


aircraft.battery.strategy = 1           # 1:battery mass driven by power_feed, time_feed & energy_cruise, 2: battery mass is input
aircraft.battery.power_feed = 0.        # J, Power delivered to e-fan(s) at take off and(or) climb during a total of time_feed
aircraft.battery.time_feed = 0.         # s, Maximum duration of the power_feed delivered to e-fan(s)
aircraft.battery.energy_cruise = 0.     #unit.J_kWh(140)     # J, energy stored in the battery dedicated to the cruise
aircraft.battery.mass = 0.              # kg, # Battery mass, input when battery_strategy==2

# Input to be used for sensitivity studies
#------------------------------------------------------------------------------------------------------
aircraft.propulsion.bli_effect = 1                      # boundary layer effect, 0:inactive, 1:active
aircraft.battery.energy_density = unit.J_kWh(0.2)       # J/kg, # Battery energy density
aircraft.power_elec_chain.overall_efficiency = 0.90     # 0.90 from init.e_chain_efficiency()


aircraft.economics.fuel_price = 2/unit.liter_usgal(1)   # 2 $/USgal
aircraft.economics.elec_price = 0.15/unit.J_kWh(1)      # 0.05 $/kWh


# Design process
#======================================================================================================

# Automatic optimization of the airplane
#------------------------------------------------------------------------------------------------------
thrust_bnd = (110000,150000)        # TO BE UPDATED (so that bound constraints are not active)
area_bnd = (100,200)                # TO BE UPDATED (so that bound constraints are not active)

search_domain = (thrust_bnd,area_bnd)

# Criterion to be chosen among  "MTOW", "cost_fuel", "CO2_metric", "COC", "DOC"
criterion = "MTOW"

run.optimization(aircraft,search_domain,criterion)


# Result printing
#======================================================================================================

# Print relevant output  TO BE UPDATED  (according to study needs)
#------------------------------------------------------------------------------------------------------

print("--------------------------------------------------------------")
print("Engine effective reference thrust = ","%.1f"%(aircraft.propulsion.reference_thrust_effective/10)," daN")
print("Wing area = ","%.1f"%aircraft.wing.area," m2")
print("MTOW = ","%.0f"%aircraft.weights.mtow," kg")
print("OWE = ","%.0f"%aircraft.weights.owe," kg")

print("")
print("LoD cruise = ","%.2f"%(aircraft.aerodynamics.cruise_lod_max)," no_dim")
print("SFC cruise = ","%.3f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
print("SEC cruise = ","%.3f"%(aircraft.propulsion.sec_cruise_ref/100)," kW/daN")

print("")
print("MTO refernce thrust = ","%.1f"%(aircraft.propulsion.mto_thrust_ref/10)," daN")
print("MCN refernce thrust = ","%.1f"%(aircraft.propulsion.mcn_thrust_ref/10)," daN")
print("MCL refernce thrust = ","%.1f"%(aircraft.propulsion.mcl_thrust_ref/10)," daN")
print("MCR refernce thrust = ","%.1f"%(aircraft.propulsion.mcr_thrust_ref/10)," daN")

print("")
if (aircraft.propulsion.bli_thrust_factor!=None): print("Fan BLI thrust factor = ","%.4f"%aircraft.propulsion.bli_thrust_factor," no_dim")
if (aircraft.turbofan_nacelle.nozzle_width!=None): print("Turbofan nozzle width = ","%.2f"%aircraft.turbofan_nacelle.nozzle_width," m")
if (aircraft.turbofan_nacelle.fan_width!=None): print("Turbofan fan width = ","%.2f"%aircraft.turbofan_nacelle.fan_width," m")
print("Turbofan nacelle width = ","%.2f"%aircraft.turbofan_nacelle.width," m")
print("Turbofan nacelle length = ","%.2f"%aircraft.turbofan_nacelle.length," m")
print("Turbofan nacelle mass = ","%.1f"%aircraft.turbofan_nacelle.mass," kg")

print("")
if (aircraft.propulsion.bli_e_thrust_factor!=None): print("eFan BLI thrust factor = ","%.4f"%aircraft.propulsion.bli_e_thrust_factor," no_dim")
if (aircraft.electric_nacelle.width!=None): print("Electric nacelle width = ","%.1f"%aircraft.electric_nacelle.width," m")
if (aircraft.electric_nacelle.length!=None): print("Electric nacelle length = ","%.1f"%aircraft.electric_nacelle.length," m")
if (aircraft.electric_nacelle.mass!=None): print("Electric nacelle mass = ","%.1f"%aircraft.electric_nacelle.mass," kg")
if (aircraft.power_elec_chain.mass!=None): print("Power electric mass = ","%.1f"%aircraft.power_elec_chain.mass," kg")
if (aircraft.battery.mass!=None): print("Battery mass = ","%.1f"%aircraft.battery.mass," kg")

print("")
print("--------------------------------------------------------------")
print("Nominal mission fuel = ","%.1f"%(aircraft.nominal_mission.block_fuel)," kg")
print("Maximum fuel weight = ","%.1f"%(aircraft.weights.mfw)," kg")
print("")
print("Take off field length required = ","%.1f"%aircraft.low_speed.req_tofl," m")
print("Take off field length effective = ","%.1f"%aircraft.low_speed.eff_tofl," m")
print("")
print("Approach speed required = ","%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)," kt")
print("Approach speed effective = ","%.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)," kt")
print("")
print("Flight path required OEI = ","%.2f"%(aircraft.low_speed.req_oei_path*100)," %")
print("light path effective OEI = ","%.2f"%(aircraft.low_speed.eff_oei_path*100)," %")
print("")
print("Vertical speed required with MCL rating = ","%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)," ft/min")
print("Vertical speed effective with MCL rating = ","%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)," ft/min")
print("")
print("Vertical speed required with MCR rating = ","%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_cruise)," ft/min")
print("Vertical speed effective with MCR rating = ","%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_cruise)," ft/min")
print("")
print("Time to climb required = ","%.1f"%unit.min_s(aircraft.high_speed.req_ttc)," min")
print("Time to climb effective = ","%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)," min")

print("")
print("--------------------------------------------------------------")
print("MTOW = ","%.0f"%aircraft.weights.mtow," kg")
print("Cost mission fuel = ","%.1f"%aircraft.cost_mission.block_fuel," kg")
print("Cash Operating Cost = ","%.1f"%aircraft.economics.cash_operating_cost," $/trip")
print("Carbon dioxid emission = ","%.1f"%(aircraft.cost_mission.block_CO2)," kg/trip")
print("Fuel efficiency metric = ","%.4f"%(aircraft.environmental_impact.CO2_metric*1e7)," 10-7kg/km/m0.48")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
#show.draw_3d_view(aircraft,"study_n5",study_name)

aircraft.export_to_ini_file("Output.txt")
