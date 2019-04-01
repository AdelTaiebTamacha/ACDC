#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""

import numpy

from MARILib.tools import units as unit
from MARILib.aircraft_data.aircraft_description import Aircraft
from MARILib.processes import design as perform
from MARILib.aircraft_model.airplane import viewer as show

# Configuration data
#---------------------------------------------------------------------------
input = {"title":("Poussée moteur de référence","Surface de l'aile","Masse maximale"), \
         "unit":("daN","m2","kg"), \
         "default":(12000,140,74000)}

output = {"x_axis":("Poussée de référence","daN",(10000,20000)), \
          "y_axis":("Surface de l'aile","m2",(100,200)), \
          "title":("Rayon d'action demandé",
                   "Rayon d'action effectif",
                   "Longueur de piste au décollage demandée",
                   "Longueur de piste au décollage effective",
                   "Vitesse d'approche demandée",
                   "Vitesse d'approche effective",
                   "Vitesse de monté demandé",
                   "Vitesse de monté effective",
                   "Temps de monté demandé",
                   "Temps de monté effectif",
                   "Coût d'un vol"), \
         "unit":("km","km","m","m","km/h","km/h","m/s","m/s","minutes","minutes","$"), \
         "default":(5400,5400,2000,2000,240,240,1.8,1.8,25,25,12000)}

info_dict = {"input":input, "output":output}


# Initialize aircraft data structure
#---------------------------------------------------------------------------
aircraft = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propu_config = 1
n_engine = 2

perform.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

aircraft.turbofan_engine.reference_thrust = 125000  # Newtons
aircraft.wing.area = 140                # m2
aircraft.weights.mtow = 74000           # kg

perform.aircraft_pre_design(aircraft)
perform.mass_estimation(aircraft)
perform.performance_analysis(aircraft)

# Output data
#---------------------------------------------------------------------------
coordinate = (aircraft.weights.mtow, aircraft.turbofan_engine.reference_thrust, aircraft.wing.area)

performance = (aircraft.design_driver.design_range,
               aircraft.nominal_mission.range, \
               aircraft.low_speed.req_tofl, \
               aircraft.low_speed.eff_tofl, \
               aircraft.low_speed.req_app_speed, \
               aircraft.low_speed.eff_app_speed, \
               aircraft.high_speed.req_vz_climb, \
               aircraft.high_speed.eff_vz_climb, \
               aircraft.high_speed.req_ttc, \
               aircraft.high_speed.eff_ttc, \
               aircraft.economics.direct_operating_cost)

condition = ((aircraft.nominal_mission.range - aircraft.design_driver.design_range)>0, \
             (aircraft.low_speed.req_tofl - aircraft.low_speed.eff_tofl)>0, \
             (aircraft.low_speed.req_app_speed - aircraft.low_speed.eff_app_speed)>0, \
             (aircraft.high_speed.eff_vz_climb - aircraft.high_speed.req_vz_climb)>0, \
             (aircraft.high_speed.req_ttc - aircraft.high_speed.eff_ttc)>0)

if (all(condition)):
    color = "green"
else:
    color = "orange"

out_dict = {"coordinate":coordinate, "color":color, "performance":performance}

# Print data
#---------------------------------------------------------------------------
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
print("Vitesse de monté demandé = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)+" ft/min")
print(" . . . . . . . effective = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)+" ft/min")
print("")
print("Temps de monté demandé = "+"%.1f"%unit.min_s(aircraft.high_speed.req_ttc)+" min")
print(" . . . . . .  effectif = "+"%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)+" min")
print("")
print("Coût d'un voyage = "+"%.0f"%aircraft.economics.direct_operating_cost+" $")

#show.draw_3d_view(aircraft,"Test","Mon avion")
