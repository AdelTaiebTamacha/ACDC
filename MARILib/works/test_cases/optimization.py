#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from MARILib.tools import units as unit

from MARILib.aircraft_model.airplane import viewer as show
from MARILib.aircraft_model import initialization as init

from MARILib.aircraft_data.aircraft_description import Aircraft

from MARILib.processes import design as perform

# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
architecture = init.prop_architecture()

aircraft = Aircraft(architecture)

design_range = 14800000 #unit.m_NM(8000)
cruise_mach = 0.85
n_pax_ref = 240

propu_config = 1
n_engine = 4

#------------------------------------------------------------------------------------------------------
perform.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

aircraft.wing.attachment = 1

print("-------------------------------------------")
print("Initialization : done")

# Modify TLRs here

#------------------------------------------------------------------------------------------------------
thrust_bnd = (110000,800000)
area_bnd = (100,800)
search_domain = (thrust_bnd,area_bnd)

perform.optimization(aircraft,search_domain)

print("-------------------------------------------")
print("Optimization : done")

print("")
print("Wing area = ","%.0f"%aircraft.wing.area," m2")
print("Reference thrust = ","%.0f"%aircraft.propulsion.reference_thrust_effective," N")
if (propu_config==3):
    print("Reference power = ","%.0f"%aircraft.turboprop_engine.reference_power," W")
print("")
print("OWE structure = ","%.2f"%aircraft.weights.owe," kg")
print("MTOW structure = ","%.2f"%aircraft.weights.mtow," kg")

print("")
print("Rayon d'action demandé = ","%.0f"%unit.NM_m(aircraft.design_driver.design_range)," NM")
print(" . . . . . .  effectif = ","%.0f"%unit.NM_m(aircraft.nominal_mission.range)," NM")
print("")
print("Longueur de piste au décollage demandée = "+"%.0f"%aircraft.low_speed.req_tofl+" m")
print(" . . . . . . . . . . . . . .  effective = "+"%.0f"%aircraft.low_speed.eff_tofl+" m")
print("")
print("Vitesse d'approche demandée = "+"%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)+" kt")
print(" . . . . . . . .  effective = "+"%.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)+" kt")
print("")
print("Pente un moteur en panne demandée = "+"%.1f"%(aircraft.low_speed.req_oei_path*100)+" %")
print(" . . . . . . . . . . .  effective = "+"%.1f"%(aircraft.low_speed.eff_oei_path*100)+" %")
print("")
print("Vitesse de monté demandé (MCL) = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)+" ft/min")
print(" . . . . . . . effective (MCL) = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)+" ft/min")
print("")
print("Vitesse de monté demandé (MCR) = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_cruise)+" ft/min")
print(" . . . . . . . effective (MCR) = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_cruise)+" ft/min")
print("")
print("Temps de monté demandé = "+"%.1f"%unit.min_s(aircraft.high_speed.req_ttc)+" min")
print(" . . . . . .  effectif = "+"%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)+" min")
print("")
print("Coût d'un voyage = "+"%.0f"%aircraft.economics.direct_operating_cost+" $")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(aircraft,"sequence packaged","This plane")
