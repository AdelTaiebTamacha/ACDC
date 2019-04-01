#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""


from scipy.optimize import fsolve

from MARILib.aircraft_model.operations import mission as perfo


#===========================================================================================================
def mission_payload(aircraft,tow,range,altp,mach,disa):
    """
    Mission simulation (payload weight is output)
    """

    weights = aircraft.weights

    [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)

    payload = tow - weights.owe - total_fuel

    return payload,block_fuel,block_time,total_fuel


#===========================================================================================================
def mission_tow(aircraft,payload,range,altp,mach,disa):
    """
    Mission simulation (take off weight is output)
    """

    weights = aircraft.weights

    def fct_mission(tow,aircraft,payload,range,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)
        Y = tow - weights.owe - payload - total_fuel
        return Y
    #---------------------------------------------------------------------------------------

    tow_ini = weights.owe + payload + 2000

    fct_arg = (aircraft,payload,range,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = tow_ini, args=fct_arg, full_output=True)

    tow = output_dict[0][0]

    [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)

    return tow,block_fuel,block_time,total_fuel


#===========================================================================================================
def mission_range(aircraft,tow,payload,altp,mach,disa):
    """
    Mission simulation (range is output)
    """

    design_driver = aircraft.design_driver

    def fct_mission(range,aircraft,tow,payload,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)
        Y = tow - weights.owe - payload - total_fuel
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,payload,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)

    return range,block_fuel,block_time,total_fuel


#===========================================================================================================
def mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa):
    """
    Mission fuel limited (range & payload are output)
    """

    design_driver = aircraft.design_driver
    weights = aircraft.weights

    def fct_mission(range,aircraft,tow,total_fuel,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)
        Y = total_fuel - fuel
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,total_fuel,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)

    payload = tow - weights.owe - total_fuel

    return range,payload,block_fuel,block_time





