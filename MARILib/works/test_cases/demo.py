#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""

from MARILib.tools import units as unit
from MARILib.aircraft_data.aircraft_description import Aircraft
from MARILib.processes import design as perform
from MARILib.aircraft_model.operations.handling_qualities import vertical_tail_sizing
from MARILib.aircraft_model.airplane import viewer as show

# Initialize aircraft data structure
#---------------------------------------------------------------------------
aircraft = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propu_config = 1
n_engine = 2

perform.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

aircraft.turbofan_engine.reference_thrust = 130000  # Newtons
aircraft.wing.area = 155                # m2
aircraft.wing.aspect_ratio = 9                # m2
aircraft.weights.mtow = 77000   # kg

aircraft.cost_mission.range = unit.m_NM(500)

perform.aircraft_pre_design(aircraft)
#perform.mass_estimation(aircraft)
perform.mass_mission_adaptation(aircraft)
perform.performance_analysis(aircraft)
perform.payload_range_analysis(aircraft)
perform.handling_qualities_analysis(aircraft)
d_vtp_area = vertical_tail_sizing(aircraft)

print("d_vtp_area = ",d_vtp_area)

print("")
print("LoD cruise (LoD max) = ","%.2f"%aircraft.high_speed.cruise_lod," no_dim")
print("VTP area = ","%.2f"%aircraft.vertical_tail.area," m2")
print("HTP area = ","%.2f"%aircraft.horizontal_tail.area," m2")
print("HTP lever arm = ","%.2f"%aircraft.horizontal_tail.lever_arm," m")
print("Y nacelle = ","%.2f"%aircraft.turbofan_nacelle.y_ext," m")
print("Wing span = ","%.2f"%aircraft.wing.span," m")
print("Wing cg = ","%.2f"%aircraft.wing.c_g," m")
print("Wing mass = ","%.2f"%aircraft.wing.mass," kg")
print("Wing xnp = ","%.2f"%(aircraft.wing.x_mac+0.25*aircraft.wing.mac)," kg")
print("Fuselage cg = ","%.2f"%aircraft.fuselage.c_g," m")
print("Fuselage mass = ","%.2f"%aircraft.fuselage.mass," kg")
print("Fuselage length = ","%.2f"%aircraft.fuselage.length," m")
print("Propulsion cg = ","%.2f"%aircraft.propulsion.c_g," m")
print("Propulsion mass = ","%.2f"%aircraft.propulsion.mass," kg")
print("System cg = ","%.2f"%aircraft.systems.c_g," m")
print("System mass = ","%.2f"%aircraft.systems.mass," kg")
print("LDG cg = ","%.2f"%aircraft.landing_gears.c_g," m")
print("LDG mass = ","%.2f"%aircraft.landing_gears.mass," kg")
print("HTP cg = ","%.2f"%aircraft.horizontal_tail.c_g," m")
print("HTP mass = ","%.2f"%aircraft.horizontal_tail.mass," kg")
print("VTP cg = ","%.2f"%aircraft.vertical_tail.c_g," m")
print("VTP mass = ","%.2f"%aircraft.vertical_tail.mass," kg")
print("VTP AR = ","%.2f"%aircraft.vertical_tail.aspect_ratio," no_dim")
print("")
print("MTOW = ","%.2f"%aircraft.weights.mtow," kg")
print("OWE = ","%.2f"%aircraft.weights.owe," kg")
print("Total mission fuel = ","%.2f"%aircraft.nominal_mission.total_fuel," kg")
print("Payload = ","%.3f"%aircraft.payload.nominal," kg")

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
print("COC = "+"%.0f"%aircraft.economics.cash_operating_cost+" $")

show.draw_3d_view(aircraft,"Test","My plane")

