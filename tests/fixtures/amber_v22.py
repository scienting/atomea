# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/atomea
#
# See the LICENSE.md file for full license terms.

import pytest


@pytest.fixture
def amber_header():
    return b"""          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Thu Apr 14 14:06:37 2022
"""


@pytest.fixture
def results_step_data():
    return b"""--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------

 ---------------------------------------------------
 APPROXIMATING switch and d/dx switch using CUBIC SPLINE INTERPOLATION
 using   5000.0 points per unit in tabled values
 TESTING RELATIVE ERROR over r ranging from 0.0 to cutoff
| CHECK switch(x): max rel err =   0.2738E-14   at   2.422500
| CHECK d/dx switch(x): max rel err =   0.8314E-11   at   2.736960
 ---------------------------------------------------
|---------------------------------------------------
| APPROXIMATING direct energy using CUBIC SPLINE INTERPOLATION
|  with   50.0 points per unit in tabled values
| Relative Error Limit not exceeded for r .gt.   2.33
| APPROXIMATING direct force using CUBIC SPLINE INTERPOLATION
|  with   50.0 points per unit in tabled values
| Relative Error Limit not exceeded for r .gt.   2.80
|---------------------------------------------------

 NSTEP =      500   TIME(PS) =    1021.000  TEMP(K) =   302.32  PRESS =     0.0
 Etot   =   -109458.5000  EKtot   =     20735.8122  EPtot      =   -130194.3121
 BOND   =       678.7182  ANGLE   =      1926.4968  DIHED      =      1263.9976
 UB     =         0.0000  IMP     =         0.0000  CMAP       =       159.8573
 1-4 NB =       800.8682  1-4 EEL =     10477.9827  VDWAALS    =     19893.6001
 EELEC  =   -165395.8331  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    332147.5960
                                                    Density    =         1.0361
 Ewald error estimate:   0.2092E-04
 ------------------------------------------------------------------------------

| MC Barostat: Decreasing size of volume moves

 NSTEP =     1000   TIME(PS) =    1022.000  TEMP(K) =   300.63  PRESS =     0.0
 Etot   =   -109561.8383  EKtot   =     20619.8891  EPtot      =   -130181.7274
 BOND   =       662.5735  ANGLE   =      1907.5771  DIHED      =      1284.3333
 UB     =         0.0000  IMP     =         0.0000  CMAP       =       149.8874
 1-4 NB =       829.1809  1-4 EEL =     10493.5073  VDWAALS    =     19275.6156
 EELEC  =   -164784.4025  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    332147.5960
                                                    Density    =         1.0361
 Ewald error estimate:   0.3496E-04
 ------------------------------------------------------------------------------


 NSTEP =     1500   TIME(PS) =    1023.000  TEMP(K) =   299.65  PRESS =     0.0
 Etot   =   -109729.4502  EKtot   =     20553.3061  EPtot      =   -130282.7563
 BOND   =       633.6869  ANGLE   =      1953.8346  DIHED      =      1286.2385
 UB     =         0.0000  IMP     =         0.0000  CMAP       =       150.5271
 1-4 NB =       806.4296  1-4 EEL =     10479.0387  VDWAALS    =     19484.2413
 EELEC  =   -165076.7531  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    332147.5960
                                                    Density    =         1.0361
 Ewald error estimate:   0.3773E-04
 ------------------------------------------------------------------------------

| MC Barostat: Decreasing size of volume moves

 NSTEP =     2000   TIME(PS) =    1024.000  TEMP(K) =   300.74  PRESS =     0.0
 Etot   =   -109711.6013  EKtot   =     20627.5611  EPtot      =   -130339.1625
 BOND   =       669.1943  ANGLE   =      1912.8495  DIHED      =      1239.1669
 UB     =         0.0000  IMP     =         0.0000  CMAP       =       140.7227
 1-4 NB =       788.7926  1-4 EEL =     10467.1885  VDWAALS    =     19323.3816
 EELEC  =   -164880.4585  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    332147.5960
                                                    Density    =         1.0361
 Ewald error estimate:   0.4189E-04
 ------------------------------------------------------------------------------
"""


@pytest.fixture
def system_info_data():
    return b"""--------------------------------------------------------------------------------
   1.  RESOURCE   USE: 
--------------------------------------------------------------------------------
 getting box info from netcdf restart file
 NATOM  =   33582 NTYPES =      19 NBONH =   31714 MBONA  =    1852
 NTHETH =    4022 MTHETA =    2505 NPHIH =    8139 MPHIA  =    7890
 NHPARM =       0 NPARM  =       0 NNB   =   59724 NRES   =   10270
 NBONA  =    1852 NTHETA =    2505 NPHIA =    7890 NUMBND =      89
 NUMANG =     205 NPTRA  =     205 NATYP =      47 NPHB   =       0
 IFBOX  =       1 NMXRS  =      36 IFCAP =       0 NEXTRA =       0
 NCOPY  =       0

| Coordinate Index Table dimensions:    10   11   12
| Direct force subcell size =     6.1902    6.3323    6.4194

     BOX TYPE: RECTILINEAR
"""


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
 BOND   =       678.7182  ANGLE   =      1926.4968  DIHED      =      1263.9976
 UB     =         0.0000  IMP     =         0.0000  CMAP       =       159.8573
 1-4 NB =       800.8682  1-4 EEL =     10477.9827  VDWAALS    =     19893.6001
 EELEC  =   -165395.8331  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    332147.5960
                                                    Density    =         1.0361
 Ewald error estimate:   0.2092E-04
 ------------------------------------------------------------------------------

| MC Barostat: Decreasing size of volume moves

 NSTEP =     1000   TIME(PS) =    1022.000  TEMP(K) =   300.63  PRESS =     0.0
 Etot   =   -109561.8383  EKtot   =     20619.8891  EPtot      =   -130181.7274
 BOND   =       662.5735  ANGLE   =      1907.5771  DIHED      =      1284.3333
 UB     =         0.0000  IMP     =         0.0000  CMAP       =       149.8874
 1-4 NB =       829.1809  1-4 EEL =     10493.5073  VDWAALS    =     19275.6156
 EELEC  =   -164784.4025  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    332147.5960
                                                    Density    =         1.0361
 Ewald error estimate:   0.3496E-04
 ------------------------------------------------------------------------------

      A V E R A G E S   O V E R    1000 S T E P S

--------------------------------------------------------------------------------
   5.  TIMINGS
--------------------------------------------------------------------------------

|  Master Total wall time:        1230    seconds     0.34 hours
"""
