#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""


from MARILib.tools import units as unit

from MARILib.aircraft_data.aircraft_description import Aircraft

from MARILib.processes import design as perform

from MARILib.aircraft_model.airplane import viewer as show


# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
aircraft = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propu_config = 1
n_engine = 2

#------------------------------------------------------------------------------------------------------
perform.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

#aircraft.turboprop_engine.reference_thrust = 75000
#aircraft.wing.area = 90

print("-------------------------------------------")
print("Initialization : done")

#------------------------------------------------------------------------------------------------------
perform.aircraft_pre_design(aircraft)

print("-------------------------------------------")
print("Pre design : done")

#------------------------------------------------------------------------------------------------------
#perform.mass_estimation(aircraft)
perform.mass_mission_adaptation(aircraft)

print("-------------------------------------------")
print("Mass & CG estimation : done")

#------------------------------------------------------------------------------------------------------
perform.performance_analysis(aircraft)

print("-------------------------------------------")
print("Nominal mission & Performance analysis : done")

#------------------------------------------------------------------------------------------------------
perform.payload_range_analysis(aircraft)

print("-------------------------------------------")
print("Payload-Range analysis : done")

#------------------------------------------------------------------------------------------------------
perform.handling_qualities_analysis(aircraft)

print("-------------------------------------------")
print("Handling qualities analysis : done")

print("")
print("Rayon d'action demandé = ","%.0f"%unit.NM_m(aircraft.design_driver.design_range)," NM")
print(" . . . . . . .  effectif = ","%.0f"%unit.NM_m(aircraft.nominal_mission.range)," NM")
print("")
print("Longueur de piste au décollage demandée = "+"%.0f"%aircraft.low_speed.req_tofl+" m")
print(" . . . . . . . . . . . . . . . effective = "+"%.0f"%aircraft.low_speed.eff_tofl+" m")
print("")
print("Vitesse d'approche demandée = "+"%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)+" kt")
print(" . . . . . . . . .  effective = "+"%.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)+" kt")
print("")
print("Vitesse de monté demandé = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)+" ft/min")
print(" . . . . . . . . effective = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)+" ft/min")
print("")
print("Temps de monté demandé = "+"%.1f"%unit.min_s(aircraft.high_speed.req_ttc)+" min")
print(" . . . . . . . effectif = "+"%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)+" min")
print("")
print("Coût d'un voyage = "+"%.0f"%aircraft.economics.direct_operating_cost+" $")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(aircraft,"sequence packaged","This plane")

