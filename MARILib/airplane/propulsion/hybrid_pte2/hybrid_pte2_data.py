#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class BodyNacelle():
    """
    Body nacelle characteristics
    """
    def __init__(self, width = None,
                       length = None,
                       x_axe = None,
                       y_axe = None,
                       z_axe = None,
                       net_wetted_area = None,
                       mass = None,
                       c_g = None):
        """
        Constructor :
            :param width: Uu m - OoM 10^0 - Maximum width of the electric fan cowl
            :param length: Uu m - OoM 10^0 - Length of the electric fan cowl
            :param x_axe: Uu m - OoM 10^1 - Longitudinal position of the center of the electric nacelle air inlet
            :param y_axe: Uu m - OoM 10^1 - Span wise position of the center of the electric nacelle air inlet
            :param z_axe: Uu m - OoM 10^0 - Vertical position of the center of the electric nacelle air inlet
            :param net_wetted_area: Uu m2 - OoM 10^1 - Total net wetted area of the electric fan nacelle (fan cowl)
            :param mass: Uu kg - OoM 10^2 - Equiped mass of the nacelle of the electric fan (including the controler, motor and nacelle)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the electric nacelle
        """
        self.width = width
        self.length = length
        self.x_axe = x_axe
        self.y_axe = y_axe
        self.z_axe = z_axe
        self.net_wetted_area = net_wetted_area
        self.mass = mass
        self.c_g = c_g

