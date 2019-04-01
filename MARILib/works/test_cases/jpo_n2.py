#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry

This study allows to illustrate the integration of a partial turbo-electric airplane
"""

from MARILib.tools import units as unit

from MARILib.aircraft_data.aircraft_description import Aircraft

from MARILib.processes import design as run

# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
aircraft = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propu_config = 2
n_engine = 2

# Initialize all input data
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach,propu_config, n_engine)


# Possibility to modify initial values
#------------------------------------------------------------------------------------------------------
study_name = "This airplane"

aircraft.battery.strategy = 2

aircraft.propulsion.bli_effect = 1                      # 1: with, 0: without
aircraft.power_elec_chain.overall_efficiency = 0.90     # 0.90 from init.e_chain_efficiency()

aircraft.economics.battery_price = 20                   # $/kg

aircraft.economics.elec_price = 0.15/unit.J_kWh(1)      # 0.05 $/kWh
aircraft.economics.fuel_price = 2/unit.liter_usgal(1)   # 2 $/USgal

aircraft.battery.energy_density = unit.J_kWh(0.2)*4     # J/kg, # Battery energy density

aircraft.battery.mass = 0   # kg, # Battery mass  [0 to 2000 kg]

e_power = 1.0e6     # Watts, electric motor power [0 to 1 MW]

aircraft.power_elec_chain.mto = e_power
aircraft.power_elec_chain.mcn = e_power
aircraft.power_elec_chain.mcl = e_power
aircraft.power_elec_chain.mcr = e_power


#------------------------------------------------------------------------------------------------------
thrust_bnd = (110000,150000)
area_bnd = (100,200)
search_domain = (thrust_bnd,area_bnd)

run.optimization(aircraft,search_domain)


# Print relevant output
#------------------------------------------------------------------------------------------------------
print("")
print("Engine thrust = ","%.1f"%(aircraft.propulsion.reference_thrust_effective/10)," daN")
print("Wing area = ","%.1f"%aircraft.wing.area," m2")
print("MTOW = ","%.0f"%aircraft.weights.mtow," kg")
print("OWE = ","%.0f"%aircraft.weights.owe," kg")

print("")
print("MCR power ratio = ","%.4f"%(aircraft.electric_engine.mcr_e_power_ratio)," no_dim")
print("Turbofan nacelle mass = ","%.1f"%aircraft.turbofan_nacelle.mass," kg")
print("Electric nacelle mass = ","%.1f"%aircraft.electric_nacelle.mass," kg")
print("Power electric mass = ","%.1f"%aircraft.power_elec_chain.mass," kg")
print("Battery mass = ","%.1f"%aircraft.battery.mass," kg")

print("")
print("LoD cruise = ","%.2f"%(aircraft.aerodynamics.cruise_lod_max)," no_dim")
print("SFC cruise = ","%.3f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
print("SEC cruise = ","%.3f"%(aircraft.propulsion.sec_cruise_ref/100)," kW/daN")

print("")
print("Take off field length required = ","%.1f"%aircraft.low_speed.req_tofl," m")
print("Take off field length effective = ","%.1f"%aircraft.low_speed.eff_tofl," m")
print("")
print("Approach speed required = ","%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)," kt")
print("Approach speed effective = ","%.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)," kt")
print("")
print("Vertical speed required = ","%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)," ft/min")
print("Vertical speed effective = ","%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)," ft/min")
print("")
print("Time to climb required = ","%.1f"%unit.min_s(aircraft.high_speed.req_ttc)," min")
print("Time to climb effective = ","%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)," min")

print("")
print("Cash Operating Cost = ","%.1f"%aircraft.economics.cash_operating_cost," $/trip")
print("Carbon dioxid emission = ","%.1f"%(aircraft.cost_mission.block_CO2)," kg/trip")
print("")
print("Fuel efficiency metric = ","%.1f"%(aircraft.environmental_impact.CO2_metric*1e7)," 10e-7kg/km/m0.48")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
#show.draw_3d_view(aircraft,"study_n5",study_name)

