# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import pytest


@pytest.fixture
def results_step_data():
    return b""" NSTEP =      500   TIME(PS) =    1021.000  TEMP(K) =   302.32  PRESS =     0.0
 Etot   =   -109458.5000  EKtot   =     20735.8122  EPtot      =   -130194.3121
 BOND   =       678.7182  ANGLE   =      1926.4968  DIHED      =      1263.9976
 UB     =         0.0000  IMP     =         0.0000  CMAP       =       159.8573
 1-4 NB =       800.8682  1-4 EEL =     10477.9827  VDWAALS    =     19893.6001
 EELEC  =   -165395.8331  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    332147.5960
                                                    Density    =         1.0361
 Ewald error estimate:   0.2092E-04
 ------------------------------------------------------------------------------"""


@pytest.fixture
def complete_amber_file():
    """A minimal but complete Amber output file"""
    return b"""          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

File Assignments:
|   MDIN: input.in
|  MDOUT: output.out

 Here is the input file:

Test simulation
&cntrl
    nstlim=1000,
&end

--------------------------------------------------------------------------------
   1.  RESOURCE   USE: 
--------------------------------------------------------------------------------
 NATOM  =   33582 NTYPES =      19 NBONH =   31714 MBONA  =    1852
 NRES   =   10270

--------------------------------------------------------------------------------
   2.  CONTROL  DATA  FOR  THE  RUN
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
   3.  ATOMIC COORDINATES AND VELOCITIES
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------

 NSTEP =      500   TIME(PS) =    1021.000  TEMP(K) =   302.32  PRESS =     0.0
 Etot   =   -109458.5000  EKtot   =     20735.8122  EPtot      =   -130194.3121

 NSTEP =     1000   TIME(PS) =    1022.000  TEMP(K) =   300.63  PRESS =     0.0
 Etot   =   -109561.8383  EKtot   =     20619.8891  EPtot      =   -130181.7274

      A V E R A G E S   O V E R    1000 S T E P S

--------------------------------------------------------------------------------
   5.  TIMINGS
--------------------------------------------------------------------------------

|  Master Total wall time:        1230    seconds     0.34 hours
"""
