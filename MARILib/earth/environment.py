#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

import numpy


#===========================================================================================================
def gravity():
    """
    Reference gravity acceleration
    """
    g = 9.80665     # Gravity acceleration at sea level
    return g

#===========================================================================================================
def sea_level_density():
    """
    Reference air density at sea level
    """
    rho0 = 1.225    # (kg/m3) Air density at sea level
    return rho0

#===========================================================================================================
def gaz_constant():
    """
    Ideal gaz constant
    """
    R = 287.05      # (J/kg/K) Ideal gaz constant for the air
    return R

#===========================================================================================================
def heat_ratio():
    """
    Reference air heat capacity ratio Cp/Cv
    """
    gam = 1.4       # Heat capacity ratio for the air
    return gam

#===========================================================================================================
def heat_constant(gam,R):
    """
    Reference air heat at constant pressure
    """
    Cp = gam*R/(gam-1) # Heat constant at constant pressure
    return Cp

#===========================================================================================================
def atmosphere(altp,disa):
    """
    Atmosphere data from ground to 50 km
    """
    if (altp<=11000):
        dtodz = -0.0065
        tstd = 288.15 + dtodz*altp
        tamb = tstd + disa
        pamb = (8.961962852-altp*2.021612304E-4)**5.255880
    elif (altp<=20000):
        dtodz = 0
        tstd = 216.65
        tamb = tstd + disa
        pamb = 128244.6928*numpy.exp(altp*(-1.576885E-4))
    elif (altp<=32000):
        dtodz = 0.0010
        tstd = 216.65 + dtodz*(altp-20000)
        tamb = tstd + disa
        pamb = (0.7055184555+altp*3.587686018E-6)**(-34.16322)
    elif (altp<=47000):
        dtodz = 0.0028
        tstd = 228.65 + dtodz*(altp-32000)
        tamb = tstd + disa
        pamb = (0.3492686141+altp*7.033096869E-6)**(-12.20115)
    elif (altp<=50000):
        dtodz = 0
        tstd = 270.65
        tamb = tstd + disa
        pamb = 41828.42421*numpy.exp(altp*(-1.2622656E-4))
    else:
        raise Exception("atmosphere_, altitude cannot exceed 50km")

    return pamb,tamb,tstd,dtodz

#===========================================================================================================
def pressure_altitude(pamb):
    """
    Pressure altitude from ground to 50 km
    """
    if (pamb>=22632.):  altp = 44330.76923 - 4946.546863 * pamb**0.1902630958
    elif (pamb>=5474.87): altp = 74588.16198 - 6341.616541 * numpy.log(pamb)
    elif (pamb>=868.014): altp = -196650 + 278731.1919 / pamb**0.02927124551
    elif (pamb>=110.906): altp = -49660.71227 + 142184.8751 / pamb**0.08195948743
    elif (pamb>=75.9443): altp = 84303.42544 - 7922.262953 * numpy.log(pamb)
    return altp

#===========================================================================================================
def air_density(pamb,tamb):
    """
    Ideal gaz density
    """
    R = gaz_constant()
    rho0 = sea_level_density()
    rho = pamb / ( R * tamb )
    sig = rho / rho0
    return rho, sig

#===========================================================================================================
def sound_speed(tamb):
    """
    Sound speed for ideal gaz
    """
    R = gaz_constant()
    gam = heat_ratio()
    vsnd = numpy.sqrt( gam * R * tamb )
    return vsnd

#===========================================================================================================
def total_temperature(tamb,mach):
    gam = heat_ratio()
    ttot = tamb*(1+(gam-1)/2*mach**2)
    return ttot

#===========================================================================================================
def total_pressure(pamb,mach):
    gam = heat_ratio()
    ptot = pamb*(1+(gam-1)/2*mach**2)**(gam/(gam-1))
    return ptot

#===========================================================================================================
def vtas_from_mach(altp,disa,mach):
    """
    Subsonic only
    """
    pamb,tamb,tstd,dtodz = atmosphere(altp,disa)
    vsnd = sound_speed(tamb)
    vtas = vsnd*mach
    return vtas

#===========================================================================================================
def mach_from_vcas(pamb,Vcas):
    """
    Subsonic only
    """
    mach = numpy.sqrt(((((0.2*(Vcas/340.29)**2+1)**3.5-1)*101325/pamb+1)**(1/3.5)-1)/0.2) ;
    return mach

#===========================================================================================================
def vcas_from_mach(pamb,mach):
    """
    Subsonic only
    """
    vcas = 340.29*numpy.sqrt(5*(((pamb*((1+0.2*mach**2)**3.5-1)/101325)+1)**(0.4/1.4)-1)) ;
    return vcas

#===========================================================================================================
def cross_over_altp(Vcas,mach):
    """
    Altitude where constant Vcas meets constant mach
    Subsonic only
    """
    pamb = ((1+0.2*(Vcas/340.29)**2)**3.5-1)*101325/((1+0.2*mach**2)**3.5-1)
    altp = pressure_altitude(pamb)
    return altp

#===========================================================================================================
def climb_mode(speed_mode,dtodz,tstd,disa,mach):
    """
    Acceleration factor depending on speed driver
    WARNING : input is mach number whatever SpeedMode
    """
    acc_factor = {
            1 : 1 + (((1+0.2*mach**2)**3.5-1)/(1+0.2*mach**2)**2.5) + 20.49029021*mach**2*tstd/(tstd+disa)*dtodz ,
            2 : 1 + 20.49029021*mach**2*tstd/(tstd+disa)*dtodz
            }.get(speed_mode, "Erreur")    # 9 is default if x not found
    return acc_factor

#===========================================================================================================
def fuel_density(fuel_type):
    """
    Reference kerosene density
    """
    if (fuel_type==1):
        fuel_density = 803. # Kerosene : between 775-840 kg/m3
    elif (fuel_type==2):
        fuel_density = 70.8 # Hydrogene : TO BE UPDATED
    else:
        raise Exception("fuel_type index is out of range")
    return fuel_density

#===========================================================================================================
def fuel_heat(fuel_type):
    """
    Reference fuel lower heating value
    """
    if (fuel_type==1):
        fuel_heat = 43.1e6 # J/kg, kerosene
    elif (fuel_type==2):
        fuel_heat = 142.1e6 # J/kg, hydrogene : TO BE UPDATED
    else:
        raise Exception("fuel_type index is out of range")
    return fuel_heat

#===========================================================================================================
def emission_index(compound):
    index = {"CO2" : 3140./1000.,
             "H2O" : 1290./1000.,
             "SO2" : 0.8/1000.,
             "NOx" : 14./1000.,
             "CO" : 3./1000.,
             "HC" : 0.4/1000.,
             "sulfuric_acid" : 0.04/1000.,
             "nitrous_acid" : 0.4/1000.,
             "nitric_acid" : 0.2/1000.,
             "soot" : 2.5e12}
    return index[compound]




#===========================================================================================================
# Acceleration factor depending on speed driver
# WARNING : input is mach number whatever SpeedMode
# La Solution ci-dessous s'appuie sur un dictionnaire qui renplace une fonction
# TO BE CALLED BY : CLIMB_MODE_[SpeedMode](dtodz,tstd,disa,mach)
#-----------------------------------------------------------------------------------------------------------
#CLIMB_MODE_ = {1 : "climb_mode_constant_CAS_",
#               2 : "climb_mode_constant_Mach_"}
#
#def climb_mode_constant_CAS(dtodz,tstd,disa,mach):
#    accfactor = 1 + (((1+0.2*mach**2)**3.5-1)/(1+0.2*mach**2)**2.5) + 20.49029021*mach**2*tstd/(tstd+disa)*dtodz ;
#    return accfactor
#
#def climb_mode_constant_mach(dtodz,tstd,disa,mach):
#    accfactor = 1 + 20.49029021*mach**2*tstd/(tstd+disa)*dtodz ;
#    return accfactor
#-----------------------------------------------------------------------------------------------------------
