#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry

This study illustrates the performance optimisation process
"""

from MARILib.tools import units as unit

from MARILib.aircraft_data.aircraft_description import Aircraft

from MARILib.processes import design as run

from MARILib.aircraft_model.airplane import viewer as show

# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
aircraft = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propu_config = 1 # 1: turbofan, 2: partial turboelectric, 3: turboprop
n_engine = 2

# Initialize all input data
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)


# Possibility to modify initial values
#------------------------------------------------------------------------------------------------------
reference_thrust_i = 120000 #aircraft.turbofan_engine.reference_thrust
area_i = 140   #aircraft.wing.area
print("")
print("Engine thrust = ","%.1f"%(reference_thrust_i/10)," daN")
print("Wing area = ","%.1f"%area_i," m2")


# Reloading updated values
#------------------------------------------------------------------------------------------------------
aircraft.turbofan_engine.reference_thrust = reference_thrust_i
aircraft.wing.area = area_i


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
print("")
print("MTOW = ","%.1f"%aircraft.weights.mtow," kg")

print("")
print("Take off field length required = "+"%.1f"%aircraft.low_speed.req_tofl+" m")
print("Take off field length effective = "+"%.1f"%aircraft.low_speed.eff_tofl+" m")
print("")
print("Approach speed required = "+"%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)+" kt")
print("Approach speed effective = "+"%.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)+" kt")
print("")
print("Vertical speed required = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)+" ft/min")
print("Vertical speed effective = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)+" ft/min")
print("")
print("Time to climb required = "+"%.1f"%unit.min_s(aircraft.high_speed.req_ttc)+" min")
print("Time to climb effective = "+"%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)+" min")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
#show.draw_3d_view(aircraft,"study_n3","This plane")

