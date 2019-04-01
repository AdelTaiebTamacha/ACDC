#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class TurbopropPylon():
    """
    Turboprop pylon mass characteristics
    """
    def __init__(self, mass = None,
                        c_g = None):
        """
        Constructor :
            :param mass: Uu kg - OoM 10^3 - Equiped mass of the pylons
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the pylons
        """
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class TurbopropNacelle():
    """
    Turbofan nacelle characteristics
    """
    def __init__(self, width = None,
                       length = None,
                       x_ext = None,
                       y_ext = None,
                       z_ext = None,
                       net_wetted_area = None,
                       mass = None,
                       c_g = None):
        """
        Constructor

            :param width: Uu m - OoM 10^0 - Maximum width of the nacelles
            :param length: Uu m - OoM 10^0 - Length of the fan cowl
            :param x_ext: Uu m - OoM 10^1 - Longitudinal position of the center of the air inlet
            :param y_ext: Uu m - OoM 10^1 - Span wise position of the center of the air inlet
            :param z_ext: Uu m - OoM 10^0 - Vertical position of the center of the air inlet
            :param net_wetted_area: Uu m2 - OoM 10^1 - Total net wetted area of the nacelles (fan cowls)
            :param mass: Uu kg - OoM 10^3 - Equiped mass of the nacelles (including engine masses)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the nacelles
        """
        self.width = width
        self.length = length
        self.x_ext = x_ext
        self.y_ext = y_ext
        self.z_ext = z_ext
        self.net_wetted_area = net_wetted_area
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class TurbopropEngine():
    """
    Turbofan engine characteristics
    """
    def __init__(self, n_engine = None,
                       reference_thrust = None,
                       reference_power = None,
                       rating_code = None,
                       rating_factor = None,
                       propeller_diameter = None,
                       propeller_efficiency = None):
        """
        Constructor :
            :param n_engine: Uu int - OoM 10^0 - Number of turboprops
            :param reference_thrust: Uu daN - OoM 10^4 - Design Reference Thrust of the engines
            :param power: Uu kW - OoM 10^4 - Reference power of the turboprop
            :param rating_code: Uu int - OoM 10^0 - Array of rating codes [0:MTO, 1:MCN, 2:MCL, 3:MCR, 4:FID]
            :param rating_factor: Uu int - OoM 10^0 - Array of rating factors versus maximum power
            :param propeller_diameter: Uu m - OoM 10^0 - Diameter of the propeller
            :param propeller_efficiency: Uu no_dim - OoM 10^0 - Propeller efficiency
        """
        self.n_engine = n_engine
        self.reference_thrust = reference_thrust
        self.reference_power = reference_power
        self.rating_code = rating_code
        self.rating_factor = rating_factor
        self.propeller_diameter = propeller_diameter
        self.propeller_efficiency = propeller_efficiency

