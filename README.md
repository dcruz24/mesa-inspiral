# Mesa Inspiral
Compiles under these versions..
1. pyMesa Version   2.0.1
2. gfort2py Version 2.4.2
3. Mesa Version     r23.05.1
4. Python Version   3.11

## Description
Evolve common envleope (CE) binary star systems using MESA hydrodynamics. 

1. Strict Circular Orbit: Euler Integration
2. Two body integrator:   Euler Integration

ce-circular.py:
    - Contains, main inspiral code, and where you set runtimeParameters
    
MESAGen.py: 
Contains core functions to interface between python and MESA.    
1.  Class:    pyMesa_context()             pyMesa object calling relevant libaries..
2.  Function: getEosResult_IE()            
3.  Function: getEosResult_Temperature()   
4.  Function: addInternalEnergy()          Core Function: Adds energy for all zones
5.  Function: RunMesa()                    Runs Mesa via './rn'
6.  Function: readModfile()
7.  Function: writeModfile()
8.  Function: getIsotopes()               Reads istopes from mod file and mass fractions for eos calls
9.  Function: writeInspiralParamters()
10. Function: includeCompactPotential()
11. Function: includeMassLoss()
12. Class:    fileInfo()
13. Class:    RunTimeParameters()

EnergyInjections.py:
Two core functions to run MESA.
1. Function: callMesaPart1() Takes in energy and zone arrays to add energy to the modfile.
2. Function: callMesaPart2() Increase inlist max_age and runs Mesa
        
EnergyHelper.py
Helper Functions in calculating energy distritbutions and zones to depost energy.
1. Function: getCompanionIndex()
2. Function: getAccretionZones()
3. Function: getGaussEnergy()            *No longer used...
4. Function: getUniformEnergy()          *No longer used...
5. Function: getEnergyDistribution()     *Added uniform and gaussian distributions with others here..

CEphysics.py:
Contains helpful physics and different drag descriptions..
1. getAcc2019Drag(): De, Macleod et al 2019 APJ 897 130
2. getNS2015Drag():  Macelod 2015 APJ 803:41,   * 3D Wind Tunnels, Macelod 2015 APJL 798:L19  * NS Accretion Growth
3. getJetfeedback(): Grichener 021, APJ 922 61
4. getBronnerDrag(): Bronner2023, Kim 2010 APJ 725: 1069-1081, Kim/Kim 2007, APJ 665:432,
5. getOst1999Drag(): Ostriker 1999, APJ 513 252
6. getHybirdDrag():  Fragos 2019, APJ 883:L45
        
InlistHelper.py: 
Helper Functions focused on reading and writting text, lines, and saving data. 
    
StarModels: 
Contains MESA mod, profile, and history files for stars. * May not have any due to git data limits.

pyMesaInitData.py:
Script to cache all relevant mesa libaries to use with pyMesa

pyMesaStarTest.py
Sample so far to run Mesa within Python...

Things to Consider in case things don't work...
1. Make sure Mesa is built using USE_SHARED=YES in $MESA_DIR/utils/makefile_header to work with python
2. MesaWorkDir is compiled, './mk'
3. Did you change the pyMesaUtils.py file to enable caching? Look at PyMesaStarTest.py

## Authors and acknowledgment
David Cruz Lopez  \ 
Miguel A. Holgado \
Paul M. Ricker    \

## License
For open source projects, say how it is licensed.