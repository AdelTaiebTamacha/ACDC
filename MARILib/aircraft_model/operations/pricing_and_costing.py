#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""


import numpy as numpy

from MARILib.earth import environment as earth


#===========================================================================================================
def gear_price(aircraft):
    """
    Typical value
    """

    ldg = aircraft.landing_gears

    gear_price = 720 * ldg.mass

    return gear_price


#===========================================================================================================
def engine_price(aircraft):
    """
    Regression on catalog prices
    """

    enigne = aircraft.turbofan_engine

    engine_price = ((2.115e-4*enigne.reference_thrust + 78.85)*enigne.reference_thrust)

    return engine_price


#===========================================================================================================
def airframe_price(aircraft):
    """
    Regression on catalog prices corrected with engine prices
    """

    weights = aircraft.weights
    
    airframe_price = 0.7e3*(9e4 + 1.15*weights.mwe - 1.8e9/(2e4 + weights.mwe**0.94))

    return airframe_price


#===========================================================================================================
def operating_costs(aircraft, block_fuel,block_time):
    """
    Computes Cash and Direct Operating Costs per flight (based on AAE 451 Spring 2004)
    """

    cabin = aircraft.cabin
    propulsion = aircraft.propulsion
    battery = aircraft.battery
    engine = aircraft.turbofan_engine
    weights = aircraft.weights
    cost_mission = aircraft.cost_mission
    economics = aircraft.economics

    labor_cost = economics.labor_cost
    irp = economics.irp
    period = economics.period
    interest_rate = economics.interest_rate
    utilisation = economics.utilisation

    fuel_density = earth.fuel_density(propulsion.fuel_type)

    # Cash Operating Cost
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    fuel_cost =   (block_fuel*(economics.fuel_price*1e3)/fuel_density) \
                + battery.energy_cruise*economics.elec_price

    b_h = block_time/3600
    t_t = b_h + 0.25
    w_f = (10000 + weights.mwe-propulsion.mass)*1e-5

    labor_frame = ((1.26+1.774*w_f-0.1071*w_f**2)*t_t + (1.614+0.7227*w_f+0.1204*w_f**2))*labor_cost
    matrl_frame = (12.39+29.8*w_f+0.1806*w_f**2)*t_t + (15.20+97.330*w_f-2.8620*w_f**2)
    frame_mc = labor_frame + matrl_frame

    t_h = 0.05*(propulsion.reference_thrust_effective/4.4482198)*1e-4

    labor_engine = engine.n_engine*(0.645*t_t+t_h*(0.566*t_t+0.434))*labor_cost
    matrl_engine = engine.n_engine*(25*t_t+t_h*(0.62*t_t+0.38))
    engine_mc = labor_engine + matrl_engine

    w_g = weights.mtow*1e-3

    cockpit_crew = b_h*2*(440-0.532*w_g)

    cabin_crew = b_h*numpy.ceil(cabin.n_pax_ref/50)*labor_cost

    ldg_fees = 8.66*(weights.mtow*1e-3)

    nav_fees = 57*(cost_mission.range/185200)*numpy.sqrt((weights.mtow/1000)/50)

    catering = 3.07 * cabin.n_pax_ref

    pax_handling = 2 * cabin.n_pax_ref

    ramp_handling = 8.70 * cabin.n_pax_ref

    std_op_cost = fuel_cost + frame_mc + engine_mc + cockpit_crew + ldg_fees + nav_fees 

    cash_op_cost = std_op_cost + cabin_crew + catering + pax_handling + ramp_handling

    # DirectOperating Cost
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    eng_price = engine_price(aircraft)

    ldg_price = gear_price(aircraft)       # added by NP

    frame_price = airframe_price(aircraft)

    battery_price = economics.battery_price*battery.mass

    aircraft_price = frame_price + eng_price * engine.n_engine + ldg_price + battery_price

    total_investment = frame_price * 1.06 + engine.n_engine * eng_price * 1.025

    interest = (total_investment/(utilisation*period)) * (irp * 0.04 * (((1 + interest_rate)**irp)/((1 + interest_rate)**irp - 1)) - 1)

    insurance = 0.0035 * aircraft_price/utilisation

    depreciation = 0.99 * (total_investment / (utilisation * period))     # Depreciation

    direct_op_cost = cash_op_cost + interest + depreciation + insurance

    return direct_op_cost,cash_op_cost,block_fuel,engine_price,aircraft_price



