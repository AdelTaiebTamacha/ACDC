#===========================================================================================================
def eval_hybrid_body_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimation
    """

    fuselage = aircraft.fuselage

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle
    body = aircraft.body_nacelle

    e_engine = aircraft.electric_engine
    e_nacelle = aircraft.electric_nacelle

    power_elec = aircraft.power_elec_chain

    # Body mass
    # -----------------------------------------------------------------------
    kbody = numpy.pi*body.length*body.width

    structure_mass = 5.0*kbody**1.2      # WARNING: ONE BODY MASS

    tank_mass = 0   # WARNING: ONE BODY MASS, TO BE UPDATED

    body.mass = (structure_mass + tank_mass)*engine.n_engine        # Betteries (if any) are accounted separetly (see aircraft.battery)

    # -----------------------------------------------------------------------
    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    e_shaft_power = numpy.array([e_engine.mto_e_shaft_power,
                                 e_engine.mcn_e_shaft_power,
                                 e_engine.mcl_e_shaft_power,
                                 e_engine.mcr_e_shaft_power,
                                 e_engine.fid_e_shaft_power])

    shaftPowerMax = max(e_shaft_power)

    turboFanMass0 = 1250 + 0.021*engine.reference_thrust # Statistical regression

    turboFanMass1 = 1250 + 0.021*engine.reference_thrust*engine.kfn_off_take

    kTurboFanMass = turboFanMass1 / turboFanMass0

    kMass = kTurboFanMass + engine.core_weight_ratio*(1-kTurboFanMass)     # Assuming core mass remains unchanged

    nacelle.mass = engine.n_engine * (body.mass + turboFanMass0 * kMass)     # Total engine mass

    power_elec.mass = (  1/power_elec.generator_pw_density + 1/power_elec.rectifier_pw_density \
                       + 1/power_elec.wiring_pw_density + 1/power_elec.cooling_pw_density \
                       ) * shaftPowerMax

    e_nacelle.mass = (  1/e_nacelle.controler_pw_density + 1/e_nacelle.motor_pw_density \
                      + 1/e_nacelle.nacelle_pw_density \
                      ) * shaftPowerMax

    # Propulsion system CG
    # ------------------------------------------------------------------------
    body.c_g = body.x_ext + 0.5*body.length

    nacelle.c_g = ( (nacelle.x_ext + 0.70*nacelle.length)*(nacelle.mass-body.mass) \
                   + body.c_g*body.mass \
                  )/nacelle.mass

    power_elec.c_g = nacelle.x_ext

    e_nacelle.c_g = fuselage.length + 0.5*e_nacelle.length

    return


