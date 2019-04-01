#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class Aerodynamics():
    """
    Main aerodynamic characteristics
    """
    def __init__(self, cruise_lod_max = None,
                       cz_cruise_lod_max = None,
                       hld_conf_clean = None,
                       cz_max_clean = None,
                       hld_conf_to = None,
                       cz_max_to = None,
                       hld_conf_ld = None,
                       cz_max_ld = None):
        """
        Constructor :
            :param cruise_lod_max: Uu no_dim - OoM 10^1 - Maximum lift over drag ratio in cruise
            :param cz_cruise_lod_max: Uu no_dim - OoM 10^0 - Lift coefficient corresponding to maximum lift over drag
            :param hld_conf_clean: Uu no_dim - OoM 10^0 - High lift device setting in clean configuration (0 by definition)
            :param cz_max_clean: Uu no_dim - OoM 10^0 - Maximum lift coefficient in clean wing configuration
            :param hld_conf_to: Uu no_dim - OoM 10^0 - High lift device setting in take off configuration (0 < hld_conf < 0,6)
            :param cz_max_to: Uu no_dim - OoM 10^0 - Maximum lift coefficient in take off configuration
            :param hld_conf_ld: Uu no_dim - OoM 10^0 - High lift device setting in landing configuration (nominal value is 1)
            :param cz_max_ld: Uu no_dim - OoM 10^0 - Maximum lift coefficient in landing configuration
        """
        self.cruise_lod_max = cruise_lod_max
        self.cz_cruise_lod_max = cz_cruise_lod_max
        self.hld_conf_clean = hld_conf_clean
        self.cz_max_clean = cz_max_clean
        self.hld_conf_to = hld_conf_to
        self.cz_max_to = cz_max_to
        self.hld_conf_ld = hld_conf_ld
        self.cz_max_ld = cz_max_ld

#--------------------------------------------------------------------------------------------------------------------------------
class Propulsion():
    """
    Turbomachine driven thrust in given conditions
    """
    def __init__(self, architecture = None,
                       fuel_type = None,
                       reference_thrust_effective = None,
                       sfc_cruise_ref = None,
                       sec_cruise_ref = None,
                       bli_effect = None,
                       bli_e_thrust_factor = None,
                       bli_thrust_factor = None,
                       rating_code = None,
                       mto_thrust_ref = None,
                       mcn_thrust_ref = None,
                       mcl_thrust_ref = None,
                       mcr_thrust_ref = None,
                       fid_thrust_ref = None,
                       mass = None,
                       c_g = None):
        """
        Constructor :
            :param architecture: Uu int - OoM 10^0 - Propulsion architecture, 1:turbofan, 2:partial turbo electric n°1, 3:partial turbo electric n°2, 4:turboprop
            :param fuel_type: Uu int - OoM 10^0 - Type of fuel, 1:kerosene, 2:hydrogene
            :param reference_thrust_effective: Uu daN - OoM 10^5 - Effective reference_thrust computed as max thrust(Mach = 0.25, ISA+15, Sea Level) / 0.8
            :param sfc_cruise_ref: Uu kg/daN/h - OoM 10^0 - Specific Fuel Consumption in cruise condition, isa, ref_cruise_altp, cruise_mach
            :param sec_cruise_ref: Uu kg/daN/h - OoM 10^0 - Specific Energy Consumption of the electric chain (if any) in cruise condition, isa, ref_cruise_altp, cruise_mach
            :param bli_effect: Uu int - OoM 10^0 - BLI effect switch, 0: without, 1: with
            :param bli_e_thrust_factor: Uu no_dim - OoM 10^0 - Thrust factor at constant power due to boundary layer ingestion of the e fan
            :param bli_thrust_factor: Uu no_dim - OoM 10^0 - Thrust factor at constant power due to boundary layer ingestion of other fans
            :param rating_code: Uu int - OoM 10^0 - Array of rating codes [0:MTO, 1:MCN, 2:MCL, 3:MCR, 4:FID]
            :param mto_thrust_ref: Uu daN - OoM 10^4 - Turbofan thrust in take off rating (one engine), Sea Level, ISA+15, Mach 0,25
            :param mcn_thrust_ref: Uu daN - OoM 10^4 - Turbofan thrust in maxi continuous rating (one engine), Required ceiling altitude, ISA, cruise Mach
            :param mcl_thrust_ref: Uu daN - OoM 10^4 - Turbofan thrust in max climb rating (one engine), Required Top of Climb altitude, ISA, cruise Mach
            :param mcr_thrust_ref: Uu daN - OoM 10^4 - Turbofan thrust in max cruise rating (one engine), Reference cruise altitude, ISA, cruise Mach
            :param fid_thrust_ref: Uu daN - OoM 10^4 - Turbofan thrust in flight idle rating (one engine), Reference cruise altitude, ISA, cruise Mach
            :param mass: Uu kg - OoM 10^3 - Total mass of the propulsion system (pylons, nacelles, engines, ...)
            :param c_g: Uu m - OoM 10^1 - Global CG position for thr whole propulsion system (pylons, nacelles, engines, ...)
        """
        self.architecture = architecture
        self.fuel_type = fuel_type
        self.reference_thrust_effective = reference_thrust_effective
        self.sfc_cruise_ref = sfc_cruise_ref
        self.sec_cruise_ref = sec_cruise_ref
        self.bli_effect = bli_effect
        self.bli_e_thrust_factor = bli_e_thrust_factor
        self.bli_thrust_factor = bli_thrust_factor
        self.rating_code = rating_code
        self.mto_thrust_ref = mto_thrust_ref
        self.mcn_thrust_ref = mcn_thrust_ref
        self.mcl_thrust_ref = mcl_thrust_ref
        self.mcr_thrust_ref = mcr_thrust_ref
        self.fid_thrust_ref = fid_thrust_ref
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class CharacteristicWeight():
    """
    Aircraft characteristic weights
    """
    def __init__(self, mwe = None,
                       owe = None,
                       mzfw = None,
                       mlw = None,
                       mtow = None,
                       mfw = None):
        """
        Constructor :
            :param mwe: Uu kg - OoM 10^4 - Manufacturer Weight Empty
            :param owe: Uu kg - OoM 10^4 - Operating Weight Empty (= mwe + m_op_item + m_cont_pallet)
            :param mzfw: Uu kg - OoM 10^4 - Maximum Zero Fuel Weight (= owe + n_pax_ref.m_pax_max)
            :param mlw: Uu kg - OoM 10^4 - Maximum Landing Weight (close or equal to 1,07mzfw except for small aircraft where mlw = mtow)
            :param mtow: Uu kg - OoM 10^4 - Maximum Take Off Weight
            :param mfw: Uu kg - OoM 10^4 - Maximum Fuel Weight
        """
        self.mwe = mwe
        self.owe = owe
        self.mzfw = mzfw
        self.mlw = mlw
        self.mtow = mtow
        self.mfw = mfw

#--------------------------------------------------------------------------------------------------------------------------------
class CenterOfGravity():
    """
    Required operational positions of the center of gravity
    """
    def __init__(self, mwe = None,
                       owe = None,
                       max_fwd_mass = None,
                       max_fwd_req_cg = None,
					   max_fwd_trim_cg = None,
                       max_bwd_mass = None,
                       max_bwd_req_cg = None,
                       max_bwd_stab_cg = None):
        """
        Constructor :
            :param mwe: Uu m - OoM 10^1 - Longitudinal position of MWE CG
            :param owe: Uu m - OoM 10^1 - Longitudinal position of OWE CG
            :param max_fwd_mass: Uu kg - OoM 10^2 - Aircraft mass at maximum forward CG
            :param max_fwd_req_cg: Uu m - OoM 10^1 - Required maximum forward aircraft CG
            :param max_fwd_trim_cg: Uu m - OoM 10^1 - Maximum trim-able forward CG
            :param max_bwd_mass: Uu kg - OoM 10^2 - Aircraft mass at maximum backward payload CG
            :param max_bwd_req_cg: Uu m - OoM 10^1 - Required maximum backward aircraft CG
            :param max_bwd_stab_cg: Uu m - OoM 10^1 - Maximum backward CG
        """
        self.mwe = mwe
        self.owe = owe
        self.max_fwd_mass = max_fwd_mass
        self.max_fwd_req_cg = max_fwd_req_cg
        self.max_fwd_trim_cg = max_fwd_trim_cg
        self.max_bwd_mass = max_bwd_mass
        self.max_bwd_req_cg = max_bwd_req_cg
        self.max_bwd_stab_cg = max_bwd_stab_cg

