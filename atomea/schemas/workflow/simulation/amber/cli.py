from pydantic import BaseModel

from ....io import YamlIO


class AmberCLIBase(BaseModel, YamlIO):
    mdin: str | None = None
    """
    Path to input file for controlling AMBER calculations and operations. Options
    specified in [`AmberInputsBase`]
    [schemas.workflow.simulation.amber.inputs.AmberInputsBase] should be in this file.
    """

    mdout: str | None = None
    """Output user readable state info and diagnostics. `-o stdout` will send output to
    stdout (i.e., the terminal) instead of to a file. This stream will contain the
    following information.

    **Characteristics of the Amber release.**

    ```text

            -------------------------------------------------------
            Amber 22 PMEMD                              2022
            -------------------------------------------------------

    | PMEMD implementation of SANDER, Release 22

    |  Compiled date/time: Thu Apr 14 14:06:37 2022
    ```

    **File paths.**

    ```text
    File Assignments:
    |   MDIN: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    |  MDOUT: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    | INPCRD: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    |   PARM: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    | RESTRT: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    |   REFC: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    |  MDVEL: mdvel
    |   MDEN: mden
    |  MDCRD: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    | MDINFO: /bgfs/jdurrant/amm503/oasci/metalflare/study/data/001-rogfp-md/simulat
    |LOGFILE: logfile
    |  MDFRC: mdfrc
    ```

    Then it will provide any notes and information about the chosen methods and system.

    Afterwords, system information will be printed every
    [`ntpr`][schemas.workflow.simulation.amber.inputs.AmberInputsBase.ntpr] steps.

    ```text
    NSTEP =      500   TIME(PS) =    1021.000  TEMP(K) =   300.64  PRESS =     0.0
    Etot   =   -103950.9351  EKtot   =     19726.5940  EPtot      =   -123677.5291
    BOND   =       650.1569  ANGLE   =      1869.0536  DIHED      =      1260.6252
    UB     =         0.0000  IMP     =         0.0000  CMAP       =       155.1250
    1-4 NB =       787.8972  1-4 EEL =     10500.5192  VDWAALS    =     18651.1452
    EELEC  =   -157552.0512  EHBOND  =         0.0000  RESTRAINT  =         0.0000
    EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    317418.4248
                                                        Density    =         1.0361
    Ewald error estimate:   0.1408E-04
    ------------------------------------------------------------------------------
    ```

    The meaning of each term is below for convenience.
    You can find more information in the Amber manual.

    -   `NSTEP`: Current step number.
    -   `TIME(PS)`: Cumulative simulation time in picoseconds.
        The previous simulation time is included here if the simulation is a restart.
    -   `TEMP(K)`: Instantaneous system temperature in Kelvin.
    -   `PRESS`: System pressure.
    -   `Etot`: Total energy as the sum of kinetic and potential energies.
    -   `EKtot`: Total kinetic energy.
    -   `EPtot`: Total potential energy.
    -   `BOND`: Sum of bond energies.
    -   `ANGLE`: Sum of angle energies.
    -   `DIHED`: Sum of dihedral energies.
    -   `UB`: Urey-Bradley energy term in the CHARMM force field.
    -   `IMP`: Improper energy term as in the CHARMM force field.
    -   `CMAP`: Correction map for dihedral energies
        [first introduced in CHARMM](https://doi.org/10.1002/jcc.20082) in the
        [extended in Amber's ff19SB](https://doi.org/10.1021/acs.jctc.9b00591)
        (replaces cosine-based $\phi$/$\psi$ dihedral terms in ff14SB).
    -   `1-4 NB`: 1-4 non-bonded (i.e., van der Waal) energy between atoms separated
        by three consecutive bonds. These interactions scale unlike the standard terms,
        so they are treated separately.
    -   `1-4 EEL`: 1-4 electrostatic energy between atoms separated by three consecutive bonds.
        These interactions scale unlike the standard terms, so they are treated separately.
    -   `VDWAALS`: van der Waals energy.
    -   `EELEC`: Electrostatic energy.
    -   `EHBOND`: Depreciated hydrogen bonding contributions.
    -   `RESTRAINT`: Energy of any `ntr` constraints.
    -   `EKCMT`: ?
    -   `VIRIAL`: ?
    -   `VOLUME`: System volume in Å<sup>3</sup>.
    -   `Density`: System density in g/cm<sup>3</sup>.
    """

    mdinfo: str | None = None
    """output latest mdout-format energy info
    """

    prmtop: str | None = None
    """input molecular topology, force field, periodic box type, atom and residue names
    """

    inpcrd: str | None = None
    """input initial coordinates and (optionally) velocities and periodic box size
    """
    refc: str | None = None
    """input (optional) reference coords for position restraints; also used for targeted MD
    """
    mtmd: str | None = None
    """input (optional) containing list of files and parameters for targeted MD to
    multiple targets
    """
    mdcrd: str | None = None
    """output coordinate sets saved over trajectory
    """
    inptraj: str | None = None
    """input coordinate sets in trajectory format, when imin=5 or 6
    """
    mdvel: str | None = None
    """output velocity sets saved over trajectory
    """
    mdfrc: str | None = None
    """output force sets saved over trajectory
    """
    mden: str | None = None
    """output extensive energy data over trajectory (not synchronized with mdcrd or mdvel)
    """
    restrt: str | None = None
    """output final coordinates, velocity, and box dimensions if any - for restarting run
    """
    inpdip: str | None = None
    """input polarizable dipole file, when indmeth=3
    """
    rstdip: str | None = None
    """output polarizable dipole file, when indmeth=3
    """
    cpin: str | None = None
    """input protonation state definitions
    """
    cprestrt: str | None = None
    """protonation state definitions, final protonation states for restart
    (same format as cpin)
    """
    cpout: str | None = None
    """output protonation state data saved over trajectory
    """
    cein: str | None = None
    """input redox state definitions
    """
    cerestrt: str | None = None
    """redox state definitions, final redox states for restart (same format as cein)
    """
    ceout: str | None = None
    """output redox state data saved over trajectory
    """
    evbin: str | None = None
    """input input for EVB potentials
    """
    suffix: str | None = None
    """output this string will be added to all unspecified output files that are
    printed (for multisander runs, it will append this suffix to all output files)
    """
