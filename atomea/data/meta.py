from enum import Enum


class Cadence(Enum):
    """Frequency at which a field's data changes."""

    ENSEMBLE = 0
    """
    An ensemble describes the desired conditions of an atomistic system we are studying.
    When learning statistical mechanics, you will cover various ensembles, including microcanonical, canonical, isothermal-isobaric, and grand canonical.
    Terms like these are just shorthand for specifying a thermodynamic state involving:

    - Temperature ($T$): Average kinetic energy of all particles in the system;
    - Pressure ($P$): Average force per unit area that particles exert on barriers of the system;
    - Volume ($V$): Physical space the system occupies;
    - Composition: Identity and amount of all chemical species in the system.

    An atomistic system in the canonical ensembles has constant (1) number and identity of particles, (2) volume, and (3) temperature.
    We also specify the ensemble by stating which thermodynamic variables are held constant; in this case, $NVT$.
    This could be pure oxygen gas sealed in a 351 mL container submerged in a bath maintained at a constant temperature of 63 &deg;F (290.3722 K).
    Observing this system can provide insights (e.g., pressure, distribution functions) that are only valid under those specific conditions.
    Slowly raising the temperature of the bath to 75 &deg;F (297.0389 K) would affect our observations because its physical state has changed.
    Thus, studies across computational and chemistry have to specify the exact ensemble they are working with.

    Data that would distinctly change the ensemble are given this `ENSEMBLE` cadence.
    """

    RUN = 1
    """
    The goal of all computational studies is to reproduce an experimental, real-world ensemble.
    Achieving this is virtually impossible.
    Our computational representations of our system always have some approximations for tractability purposes.

    Furthermore, many computational approaches have parameters that would impact our observations or impact reproducibility.
    For example, we often run multiple independent trajectories of molecular simulations by generating new initial velocities to improve sampling.
    Each trajectory samples from the same ensemble, but their dynamics are not concateable.
    At other times, we are testing the effect of a particular parameter or setting.
    Changing the exhaustiveness of many docking programs could result in entirely different sets of results.

    All runs share the same desired [ensemble][data.meta.Cadence.ENSEMBLE]; however, there may be some computational differences that should keep these data separate.
    Data that _could_ change between runs are given this `RUN` cadence.
    """

    MICROSTATE = 2
    """
    A microstate is a single, distinguishable configuration of particles (i.e., atoms) sampled from a single [run][data.meta.Cadence.RUN].
    The ordering of these microstates is defined by the underlying computational technique used.
    For example, a microstate could be a

    - frame of a molecular dynamics simulation trajectory;
    - protein-ligand pose in a docking calculation;
    - transition state of a chemical reaction.

    Data that _could_ change between microstates, such as atomistic coordinates, energy, docking score, instantaneous temperature or pressure, dipole moment, electronic state, etc., are given a cadence of `MICROSTATE`.
    """
