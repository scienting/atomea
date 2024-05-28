"""Simulation contexts for Amber"""
from typing import Literal

from pydantic import BaseModel, Field

from ...io import IOBase


class AmberSchemaBase(BaseModel, IOBase):
    r"""Validate Amber contexts."""

    imin: Literal[0, 1] = Field(default=0)
    """Flag for running the energy minimization procedure.

    -   `0`: Perform molecular dynamics (MD) simulation. This mode generates
        configurations by integrating Newtonian equations of motion, allowing the
        system to sample more configurational space and cross small potential energy
        barriers.
    -   `1`: Perform energy minimization. This mode iteratively moves the atoms down
        the energy gradient to relax the structure until a sufficiently low average
        gradient is obtained. Minimization is useful for preparing a system before
        MD simulations to remove bad contacts and high-energy configurations.
    """

    irest: Literal[0, 1] = Field(default=0)
    """Flag to restart a simulation from a previously saved restart file.

    -   `0`: Do not restart the simulation. The simulation will start from initial
        coordinates and velocities as specified in the input files.
    -   `1`: Restart the simulation from a previous run using data from the restart
        file (e.g., coordinates, velocities, box dimensions). This is useful for
        continuing long simulations in segments or for resuming simulations after
        an interruption.

    Restarting a simulation can help in managing large simulations by breaking them
    into smaller, manageable segments and allows for continued simulation from the
    point of interruption without losing progress.
    """

    ntx: Literal[1, 2] = Field(default=1)
    """Option to read the coordinates from the “inpcrd” file. Only options 1 and 2 are
    supported in this releases. Other options will cause pbsa to issue a warning though
    it does not affect the energy calculation.

    -   `1`: File is read with no initial velocity information. Suitable for
        starting simulations with new systems where velocities are generated based
        on `tempi`.
    -   `2`: File is read unformatted with no initial velocity information. This is
        less common and mainly used for specific needs when dealing with unformatted
        coordinate files.
    """

    ntmin: Literal[0, 1, 2, 3, 4] = Field(default=1)
    """
    Flag for selecting minimization type. Determines the algorithm used for energy
    minimization.

    -   `0`: Full conjugate gradient minimization. The first four cycles are steepest
        descent at the start of the run and after every nonbonded pair list update.
        Conjugate gradient is slower than steepest descent when far from a minimum
        but becomes more efficient when close to the minimum.
    -   `1`: For `ncyc` cycles, the steepest descent method is used, then the conjugate
        gradient is switched on. This option combines the robustness of steepest
        descent with the efficiency of conjugate gradient, making it a recommended
        choice for many scenarios.
    -   `2`: Only the steepest descent method is used. This algorithm is popular
        because it is robust and easy to implement. It is generally effective for
        systems far from equilibrium.
    -   `3`: The XMIN method is used. This method leverages advanced optimization
        algorithms for more efficient minimization, especially useful for large or
        complex systems.
    -   `4`: The LMOD method is used. This approach uses Low-Mode Conformational
        Search combined with XMIN for energy relaxation and minimization, particularly
        effective for exploring conformational space in flexible molecules.
    """

    maxcyc: int = Field(default=9999)
    """
    Maximum number of minimization cycles allowed. This parameter sets the upper limit
    on the number of cycles the minimization algorithm will perform.
    Values typically range from `1000` to `50000`. Lower values may be sufficient for
    small systems or those close to their minimum energy state, while larger values
    can help ensure convergence in more complex or strained systems.
    The default value is set to `9999`, which is generally appropriate for a wide
    range of systems, balancing computational efficiency with the need for thorough
    energy minimization.

    tip:
        -   If the system is large or highly strained, consider increasing the value to
            ensure sufficient minimization.
        -   For smaller or well-equilibrated systems, a lower value might be adequate,
            saving computational resources.
        -   Monitor the energy convergence; if the energy stabilizes before reaching
            `maxcyc`, the minimization can be considered complete.
    """

    ncyc: int = Field(default=10)
    """
    If `ntmin = 1`, then the minimization method will be switched from steepest descent
    to conjugate gradient after NCYC cycles. Values typically range from `5` to `100`.
    Lower values will switch to conjugate gradient sooner, which can be more efficient
    for systems nearing their minimum energy state. Higher values will keep the
    minimization in steepest descent mode longer, which can be beneficial for systems
    far from equilibrium. The default value is set to `10`, which provides a balance
    between the robustness of steepest descent and the efficiency of conjugate gradient
    methods.

    tip:
        -   Adjust the value based on the size and complexity of your system. For larger
            or more complex systems, a higher value might be necessary to ensure
            sufficient initial minimization.
        -   For systems that are relatively close to their equilibrium state, a lower
            value may expedite the minimization process by switching to the more
            efficient conjugate gradient method sooner.
        -   Monitor the minimization process; if convergence is slow, consider
            increasing `ncyc` to allow more initial steepest descent cycles.
    """

    ig: int = Field(default=-1)
    """
    The seed for the pseudo-random number generator. This affects the starting velocities
    for MD simulations if `tempi` is nonzero. If this is `-1`, a random seed will be
    computed based on the current date and time. This should almost always be `-1` unless
    you are reproducing a run.

    tip:
        -   Set to `-1` to ensure that each simulation run starts with a different seed,
            providing varied initial conditions for better sampling.
        -   Use a specific integer value if you need to reproduce an exact simulation run
            for debugging or verification purposes.
        -   Ensure that `tempi` is set appropriately when using this option to affect
            starting velocities.
    """

    dt: float = Field(default=0.002)
    """
    The time step in picoseconds. This parameter defines the interval of time between
    each step in the simulation.

    -   **With SHAKE (`ntc = 2`):** The maximum value is `0.002` ps (2 fs). SHAKE constrains
        bond lengths involving hydrogen atoms, allowing for a larger time step.
    -   **Without SHAKE:** The maximum value is `0.001` ps (1 fs). Without bond constraints,
        a smaller time step is needed to maintain numerical stability.

    Since longer simulations are usually desired, the maximum value is typically used.
    However, values lower than the maximum can be used if necessary for the phenomena of
    interest, such as high-frequency motions or other specific needs.

    tip:
        -   Use `0.002` ps when using SHAKE (`ntc = 2`) for most standard simulations, as this
            allows for longer simulation times.
        -   Use `0.001` ps if SHAKE is not enabled or if the system requires a smaller time step
            for stability.
        -   Consider the nature of your simulation; if capturing high-frequency motions
            is critical, a smaller time step might be required.
    """

    nstlim: int = Field(default=1, gt=0)
    """
    Number of MD steps to perform. Multiply `nstlim` and `dt` to get the duration of the
    MD simulation in picoseconds.

    - **Calculation:** The total simulation time in picoseconds is obtained by multiplying
      `nstlim` (number of steps) by `dt` (time step). For example, with `dt = 0.002` and
      `nstlim = 10000000`, the total simulation time is 20000 ps (20 ns).

    tip:
        -   For short equilibration runs or quick tests, use lower values of `nstlim`, such as
            `5000` to `50000`.
        -   For production runs aimed at sampling equilibrium states or studying slower processes,
            higher values of `nstlim` are recommended, typically in the range of `1000000` to
            `50000000` (corresponding to 2 ns to 100 ns for `dt = 0.002`).
        -   Adjust `nstlim` based on the specific requirements of your study, considering both
            the computational resources available and the timescale of the phenomena you are
            investigating.
    """

    nscm: int = Field(default=1000)
    """
    Flag for the removal of translational and rotational center-of-mass motion every `nscm` steps.
    For periodic simulations (`ntb` is `1` or `2`), only the translational center-of-mass motion
    is removed.
    Reasonable values are between `100` and `2000`. Lower values ensure
    more frequent removal of center-of-mass motion, which can be beneficial for stability in
    certain simulations.

    tip:
        -   Use lower values (closer to `100`) for systems that may experience significant drift
            or when higher precision is required.
        -   Higher values (up to `2000`) can be used for more stable systems or to reduce computational
            overhead in very long simulations.
        -   Monitor the simulation to ensure that the chosen value effectively manages the center-of-mass
            motion without introducing artifacts or instability.
    """

    ntr: Literal[0, 1] = Field(default=0)
    """
    Flag for restraining positions of specified atoms using a harmonic potential.

    -   `0`: No constraints. The positions of all atoms are free to move according to the simulation dynamics.
    -   `1`: Constrain atoms specified in `restraintmask`. This applies a harmonic potential to the atoms
        defined in `restraintmask`, effectively fixing their positions relative to the rest of the system.

    tip:
        -   Use `0` for fully flexible simulations where no positional restraints are needed.
        -   Use `1` when specific atoms need to be restrained, such as in cases where you want to focus on a
            particular region of the system while keeping another region fixed or minimally perturbed.
        -   Ensure that `restraintmask` is properly defined to specify the atoms that require constraints.
    """

    restraint_wt: float = Field(default=4.0, gt=0.0)
    """
    The weight (in kcal mol<sup>-1</sup> Å<sup>-2</sup>) when `ntr = 1`. The form of the restraint is
    \( k (\Delta x)^2 \) where \( \Delta x \) is the deviation of the atom's coordinate from the reference
    position.

    -   Reasonable values are between `0.1` and `10.0`.
        - Lower values (`0.1` to `1.0`) allow for more movement, providing a looser restraint.
        - Higher values (`5.0` to `10.0`) significantly restrict movement, enforcing a tighter restraint.
    -   The default value is set to `4.0`, which provides a moderate restraint, balancing
        stability and flexibility.

    tip:
        -   Use lower values for systems where slight positional deviations are acceptable or desired.
        -   Use higher values for systems requiring strict positional maintenance of specified atoms.
        -   Adjust `restraint_wt` based on the specific requirements of your simulation, considering the
            nature and behavior of the atoms under restraint.
    """

    restraintmask: str = Field(default="!(@H=)")
    """
    Strings that specify the restricted atoms when `ntr = 1`. To see what atoms will be restrained, you can use
    `ambmask -p mol.prmtop -c mol.inpcrd -out pdb -find "RESTRAINT_STRING"` in `ambertools`. Here are some examples
    of `restraintmask`s and their descriptions.

    -   `"!(@H=)"`: Restrain all atoms except hydrogens.
    -   `"!(:WAT) & !(@H=)"`: Restrain all atoms except hydrogens and water molecules.
    -   `"!(:WAT) & !(@H=) & !(:Na+,Cl-)"`: Same as above, but does not restrain Na<sup>+</sup> and Cl<sup>-</sup> ions.
    -   `"!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')"`: Restrains protein, DNA, and RNA backbone atoms.

    We must include the `!(:WAT)` to avoid restraining oxygen atoms in the water molecules.

    tip:
        -   Use `"!(@H=)"` to exclude hydrogen atoms from the restraint, which is common to allow their flexibility.
        -   Use `"!(:WAT) & !(@H=)"` to exclude both hydrogens and water molecules from the restraint.
        -   Use `"!(:WAT) & !(@H=) & !(:Na+,Cl-)"` if you also want to exclude ions like Na<sup>+</sup> and Cl<sup>-</sup>.
        -   Use `"!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')"` to focus restraints on backbone atoms in proteins, DNA, and RNA.
        -   Modify the `restraintmask` string according to the specific requirements of your simulation and the atoms you want to restrict.
    """

    ntb: Literal[0, 1, 2] = Field(default=0)
    """
    Flag for periodic boundary conditions when computing non-bonded interactions.
    Bonds that cross the boundary are not supported.

    -   `0`: No periodic boundary conditions. This is suitable for simulations where boundary effects are not a concern, such as in isolated systems or gas-phase simulations.
    -   `1`: Constant volume. This maintains a fixed simulation box size, appropriate for systems where volume changes are not expected or desired.
    -   `2`: Constant pressure. This allows the simulation box to fluctuate in size to maintain constant pressure, suitable for more realistic simulations of condensed-phase systems where pressure control is needed.
    """

    ntf: Literal[1, 2, 3] = Field(default=1)
    """
    Force evaluation type. This parameter determines which interactions are considered during the force calculations.

    -   `1`: All contributions. This option includes all types of interactions in the force evaluation and is required for minimization.
    -   `2`: Ignore bond interactions involving hydrogens. This option is typically used when `ntc = 2`, meaning constraints are applied to bonds involving hydrogens (e.g., SHAKE algorithm).
    -   `3`: All bond interactions are omitted. This option is used when `ntc = 3`, which implies constraints are applied to all bonds.
    """

    ntc: Literal[1, 2, 3] = Field(default=1)
    """
    Flag for SHAKE to perform bond length constraints.

    -   `1`: No SHAKE. This option does not apply any bond length constraints and is required for minimization.
    -   `2`: Bonds involving hydrogen are constrained. This is the recommended setting for molecular dynamics (MD) simulations as it allows for larger time steps while maintaining stability.
    -   `3`: All bonds are constrained. This setting applies constraints to all bonds, providing the most rigid structure.

    (Not available for parallel or qmmm runs in sander).
    """

    cut: float = Field(default=8.0)
    """
    Specifies the nonbonded cutoff in Angstroms. This parameter sets the distance beyond which nonbonded interactions
    (such as van der Waals and electrostatic interactions) are not explicitly calculated.

    For Particle Mesh Ewald (PME), this value limits the direct space sum.
    PME is a method used to compute long-range electrostatic interactions
    efficiently by splitting them into short-range (direct space) and long-range
    (reciprocal space) components. The cutoff distance specifies the range within
    which the direct space interactions are computed explicitly, while interactions
    beyond this distance are handled in reciprocal space.

    -   `8.0` Å is commonly used in many simulations as a balance between accuracy and
        computational efficiency. It ensures that the majority of significant
        interactions are captured while keeping the computational cost manageable.
    -   `10.0` Å can provide improved accuracy, especially in systems where
        long-range interactions  are important. However, it increases the computational
        load.

    When using an implicit solvent model, the nonbonded cutoff should be set to
    `9999.0` Å, effectively making it infinite. This is because implicit solvent models
    treat solvation effects as a continuous medium, and all interactions need to be
    considered without a cutoff.
    """

    ntt: Literal[0, 1, 2, 3, 9, 10, 11] = Field(default=3)
    """Switch for temperature scaling.

    -   `0`: Constant total energy (NVE).
    -   `1`: Constant temperature using the weak-coupling algorithm.
        A single scaling factor is used for all atoms.
        Generally not recommended.
    -   `2`: Andersen-like temperature coupling scheme, in which imaginary "collisions"
        are performed with heat bath of temperature `temp0` every `vrand` steps.
    -   `3`: Use Langevin dynamics with the collision frequency `gamma_ln`.
        Since Langevin simulations are highly susceptible to "synchronization"
        artifacts, you should explicitly set `ig` to a different value every
        restart (e.g., `-1`). **[ recommended ]**
    -   `9`: Optimized Isokinetic Nose-Hoover chain ensemble (OIN).
        Implemented mainly for 3D-RISM and RESPA simulations.
    -   `10`: Stochastic Isokinetic Nose-Hoover RESPA integrator.
        Mainly used for RESPA simulations.
    -   `11`: Stochastic version of Berendsen thermostat, also known as the
        Bussi thermostat.
    """

    tempi: float = Field(default=100.0)
    """
    Initialization temperature in Kelvin. This parameter sets the initial temperature
    for the system at the start of the simulation.

    -   If `ntx = 1`, the initial velocities of the atoms are assigned from a
        Maxwellian distribution corresponding to the temperature `tempi`. This is
        typically used to start a new simulation where the initial conditions need to
        be defined.
    -   If `ntx = 5`, this parameter has no effect because the velocities are read from
        the input coordinates file, meaning the system continues from a previously
        equilibrated state.

    tip:
        -   Use a lower initial temperature (e.g., `100.0` K) to gradually equilibrate
            the system and avoid potential instabilities or large initial forces.
        -   For systems where you want to rapidly equilibrate to the target temperature,
            you can set `tempi` closer to the target simulation temperature.
        -   Ensure that `tempi` is appropriate for the system and the specific phase of
            your simulation (e.g., initial equilibration vs. production run).
    """

    temp0: float = Field(default=300.0)
    """
    Reference temperature at which the system will be kept in Kelvin. This parameter
    defines the target temperature for the simulation, around which the thermostat will
    regulate the system. It is crucial for simulations where temperature control is
    needed to mimic real-world conditions, such as simulations of biological systems
    or materials at specific temperatures.

    The default value is set to `300.0` K, which corresponds to approximately room
    temperature. This is a common target temperature for many biological simulations.
    Other typical values might include `310.0` K for physiological conditions
    (human body temperature) or other specific temperatures relevant to your study.
    """

    gamma_ln: float = Field(default=2.0)
    """
    The collision frequency, $\gamma$, when `ntt = 3`. This parameter is used in
    Langevin dynamics to control the rate of coupling between the system and a
    heat bath, thereby regulating the temperature.

    The default value is set to `2.0` ps<sup>-1</sup>, which is commonly used in many
    simulations to provide effective temperature control without overly damping the
    system's dynamics. Values typically range from `2.0` to `5.0` ps<sup>-1</sup>.
    While the physical collision frequency is about `50` ps<sup>-1</sup>, using such
    high values can overly dampen the system. Smaller values are preferred to maintain
    more realistic dynamics.
    """

    ntp: Literal[0, 1, 2] = Field(default=1)
    """
    Flag for constant pressure dynamics. This parameter controls how pressure is
    managed during the simulation.

    -   `0`: No pressure scaling. The system is run at constant volume, and no
        adjustments are made to maintain a specific pressure.
    -   `1`: Isotropic position scaling. This is the recommended setting for most
        simulations as it scales the simulation box uniformly in all directions to
        maintain constant pressure.
    -   `2`: Anisotropic pressure scaling. This option allows different scaling factors
        for each dimension and can only be used for orthogonal boxes. It is typically
        used in membrane simulations where different surface tensions exist in the
        x, y, and z directions. Solutes dissolved in water should not use this
        setting as it can introduce artifacts.
    """

    barostat: Literal[1, 2] = Field(default=2)
    """
    Flag used to control which barostat is used to regulate the pressure during the
    simulation. The barostat setting determines how the simulation box's volume is
    adjusted to maintain the desired pressure. Choosing the right barostat is important
    for the accuracy and stability of the simulation.

    -   `1`: Berendsen barostat. This method scales the box dimensions and atomic
        coordinates to achieve the desired pressure. It is simpler but can lead to
        less accurate pressure control and may not generate a true NPT ensemble.
    -   `2`: Monte Carlo barostat. This method uses Monte Carlo moves to adjust the
        box volume and is generally more accurate for maintaining constant pressure.
        It provides better sampling of the NPT ensemble and is recommended for most
        simulations.
    """

    pres0: float = Field(default=1.01325)
    """
    Reference pressure, in bar, at which the system is maintained. This is the target
    pressure for the barostat. This value is almost always used in simulations to
    mimic standard atmospheric conditions. Default value is`1.01325` bar, which is the
    standard atmospheric pressure.
    """

    mcbarint: int = Field(default=100)
    """
    Number of steps between volume change attempts performed as part of the Monte Carlo
    barostat. This determines how frequently the volume of the simulation box is
    adjusted during pressure regulation.
    """

    comp: float = Field(default=44.6)
    """
    Compressibility of the system when `npt > 0` in units of
    10<sup>-6</sup> bar<sup>-1</sup>. 44.6` x 10<sup>-6</sup> bar<sup>-1</sup>,
    appropriate for water. This value is used in pressure regulation to account for the
    compressibility of the solvent or system being simulated.
    """

    taup: float = Field(default=1.0)
    """
    Pressure relaxation time in picoseconds when `ntp > 0`. Recommended values are
    between `1.0` and `5.0` ps. This parameter controls how quickly the pressure
    adjusts to the target value. Start with `1.0` ps. If your simulations are unstable,
    consider increasing this value.
    """

    ntxo: Literal[1, 2] = Field(default=2)
    """
    Format of the final coordinates, velocities, and box size.

    - `1`: ASCII.
    - `2`: Binary NetCDF file. **[ recommended ]**

    tip:
        Use `2` for binary NetCDF files for efficient storage and faster read/write
        operations. Choose `1` for ASCII format if compatibility with older software or
        easier human readability is needed.
    """

    ntwr: int = Field(default=1000)
    """
    Write the restart file every `ntwr` steps.`1000` steps. Restart files allow
    resuming simulations from intermediate states and can help in recovering from
    interruptions. Use values between `100` and `2000` to balance file size and the
    need for frequent checkpoints. Adjust based on the length and criticality of your
    simulation; shorter intervals provide more frequent recovery points but larger file sizes.
    """

    ntpr: int = Field(default=1000)
    """
    Print energy information every `ntpr` steps in a human-readable form.
    This parameter controls how often energy and other summary information is written
    to the output file.
    """

    ntwx: int = Field(default=1000)
    """
    Coordinates are written every `ntwx` steps to the `mdcrd` file. This parameter
    controls how often the coordinates are saved, which can be used for trajectory
    analysis.
    """

    ioutfm: Literal[0, 1] = Field(default=1)
    """
    Format of coordinate and velocity trajectory files.

    - **`0`**: ASCII.
    - **`1`**: Binary NetCDF files. **[ recommended ]**

    tip:
        Use `2` for binary NetCDF files for efficient storage and faster read/write
        operations. Choose `1` for ASCII format if compatibility with older software or
        easier human readability is needed.
    """

    iwrap: Literal[0, 1] = Field(default=1)
    """
    Flag for wrapping coordinates around the periodic boundary.

    -   **`0`**: No wrapping.
    -   **`1`**: If atoms move across the periodic boundary, Amber will wrap them around
        to the other side. **[ recommended ]**
    """
