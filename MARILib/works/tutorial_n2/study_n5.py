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
from MARILib.aircraft_model.airplane import aerodynamics as craft_aero

from MARILib.processes import design as run
from MARILib.aircraft_model.airplane import viewer as show

# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
aircraft = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propulsion_architecture = 2
n_engine = 2

# Initialize all input data
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsion_architecture, n_engine)


# Possibility to modify initial values
#------------------------------------------------------------------------------------------------------
study_name = "This airplane"

e_power = 1.0e6       # Watts, electric motor power

aircraft.power_elec_chain.mto = e_power
aircraft.power_elec_chain.mcn = e_power
aircraft.power_elec_chain.mcl = e_power
aircraft.power_elec_chain.mcr = e_power

aircraft.battery.energy_cruise = unit.J_kWh(140)     # J, energy stored in the battery dedicated to the cruise
aircraft.battery.energy_density = unit.J_kWh(0.2)    # J/kg, # Battery energy density

aircraft.propulsion.bli_effect = 1                      #init.boundary_layer_effect()
aircraft.power_elec_chain.overall_efficiency = 0.90     # 0.90 from init.e_chain_efficiency()

# Solve the geometric coupling between airframe and engines
#------------------------------------------------------------------------------------------------------
run.aircraft_pre_design(aircraft)

# Estimate all mass and CGs
#------------------------------------------------------------------------------------------------------
run.mass_mission_adaptation(aircraft)

# Calculate all airplane performances
#------------------------------------------------------------------------------------------------------
run.performance_analysis(aircraft)

# Print relevant output
#------------------------------------------------------------------------------------------------------
altp = unit.m_ft(35000)
disa = 0
pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
(MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code
nei = 0
Fn,Data = propu.thrust(aircraft,pamb,tamb,cruise_mach,MTO,nei)
vsnd = earth.sound_speed(tamb)
tas = vsnd*cruise_mach
print("")
print("True air speed in cruise","%.1f"%tas," m/s")
print("Totalthrust in cruise","%.0f"%Fn," N")

print("")
print("Engine thrust = ","%.1f"%(aircraft.propulsion.reference_thrust_effective/10)," daN")
print("Wing area = ","%.1f"%aircraft.wing.area," m2")
print("MTOW = ","%.0f"%aircraft.weights.mtow," kg")
print("OWE = ","%.0f"%aircraft.weights.owe," kg")

print("")
print("Turbofan nacelle mass = ","%.1f"%aircraft.turbofan_nacelle.mass," kg")
if (aircraft.propulsion.architecture==2):
    print("Electric nacelle mass = ","%.1f"%aircraft.electric_nacelle.mass," kg")
    print("Power electric mass = ","%.1f"%aircraft.power_elec_chain.mass," kg")
    print("Battery mass = ","%.1f"%aircraft.battery.mass," kg")

print("")
print("LoD cruise = ","%.2f"%(aircraft.aerodynamics.cruise_lod_max)," no_dim")
print("SFC cruise = ","%.3f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
print("SEC cruise = ","%.3f"%(aircraft.propulsion.sec_cruise_ref/100)," kW/daN")
print("Nominal mission fuel = ","%.1f"%(aircraft.nominal_mission.block_fuel)," kg")

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
print("Fuel efficiency metric = ","%.4f"%(aircraft.environmental_impact.CO2_metric*1e7)," 10-7kg/km/m0.48")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
#show.draw_3d_view(aircraft,"study_n5",study_name)

