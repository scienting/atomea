from pydantic import BaseModel

from ....io import YamlIO


class AmberCLIBase(BaseModel, YamlIO):
    mdin: str | None = None
    """input control data for the min/md run"""
    mdout: str | None = None
    """output user readable state info and diagnostics -o stdout will send output to
    stdout (to the terminal) instead of to a file.
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
