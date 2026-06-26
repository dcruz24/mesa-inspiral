import numpy as np
import subprocess
import csv
import pandas as pd
import os
from src.io import inlist as ih
from src.config import constants as const

# -------------------------------------------------------------------------------------------------------------
# path: path to your MESA workspace
def RunMesa(path, numthreads):
    work_dir = os.path.abspath(path)
    rn_path = os.path.join(work_dir, 'rn')
    star_path = os.path.join(work_dir, 'star')

    if not os.path.exists(rn_path):
        raise FileNotFoundError(f"MESA run script not found: {rn_path}")
    if not os.access(rn_path, os.X_OK):
        raise PermissionError(f"MESA run script is not executable: {rn_path}")
    if not os.path.exists(star_path):
        raise FileNotFoundError(
            f"MESA executable not found: {star_path}. Build the work directory first with `cd {work_dir} && ./mk`."
        )
    if not os.access(star_path, os.X_OK):
        raise PermissionError(f"MESA executable is not executable: {star_path}")

    env = dict(os.environ)
    env['OMP_NUM_THREADS'] = str(numthreads)
    with open('MESA_OUTPUT.txt', 'w') as f:
        proc = subprocess.run(
            ['./rn'],
            cwd=work_dir,
            stdout=f,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            check=False,
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"MESA run failed with exit code {proc.returncode}. See MESA_OUTPUT.txt for details."
        )
    return

from src.mesa.eos import (
    PyMesaEOS,
    getEosResult_IE,
    getEosResult_Temperature,
    pyMesa_context,
)


## New mesa-inspiral internal energy modfier for pyMesa v2.0.0 NEEDS TO BE FIXED
def addInternalEnergy(pyMesaObject, modData, zones, Ezones, isotopeData, totalStarMass, args):
    dqLst   = args[0]
    tempLst = args[1]
    rhoLst  = args[2]
    isotopes = list(isotopeData.keys())
    xaLst   = list(isotopeData.values())
    species = len(isotopes)
    chem_ids = np.zeros(species)

    #print('Isotopes List: ', isotopes)
    #print('Zones: ', zones)
    #print('xalst: ', xalst)
    #print('ENergy to deposit: ', Ezones)
    
    for i, name in enumerate(isotopes):
        chem_ids[i] = pyMesaObject.chem_lib.chem_get_iso_id(name).result

    net_iso = chem_ids
    
    for i, zone in enumerate(zones):
        dq     = dqLst[zone]    # Mass Fraction
        T      = tempLst[zone]  # Tempertaure
        Told   = T
        rho    = rhoLst[zone]   # Density
        Einput = Ezones[i]            # Energy to deposit
        #xa     = xaLst[i]           # Isotope Mass Fraction
        #xa     = np.array(xa)
        xa = np.array([xaLst[j][i] for j in range(species)])
        xa = ih.changeDEModfile(xa, 'D', 'E')
        #print("xa:       ", xa)
        #print("sum: ", sum(xa))
        #print("net_iso:  ", net_iso)
        #print("chem_ids: ", chem_ids)
        
        # Calls Mesa EOS to get Internal Energy
        eos_result = getEosResult_IE(pyMesaObject, rho, T, species, chem_ids, net_iso, xa)
        res        = eos_result.args["res"]
        i_lnE      = pyMesaObject.eos_def.i_lne - 1  # index for Internal Energy
        IE         = np.exp(res[i_lnE])
        #print('Internal Energy:   {} erg/g'.format(IE))
        
        # Adding internal energy
        newIE = IE + (Einput/(dq*totalStarMass))
        #print('TESTING: ', IE)
        #print('RHS:     ', (Einput/(dq*totalStarMass)))
        #print('Einput:  ', Einput)
        #print('dq:      ', dq)
        #newIE = IE # Test that we get the same internal energy
        #print('New Internal Energy: {} erg/g'.format(newIE))
    
        # Calls Mesa EOS to get new Temperature.... Need to write as nautural log for mod file
        T = getEosResult_Temperature(pyMesaObject, rho, T, species, chem_ids, net_iso, xa, IE, newIE, i_lnE)
        Tnew = T
        #print('New addded Ln T: ', np.log(T))
        
        # Debugging
        deltaT = Tnew - Told
        rel_change = deltaT / Told
        #f abs(rel_change) > 0.2:
        #    print(f"[WARN] Zone {zone}: ΔT/T = {rel_change:.2f}, T_old = {Told:.3e}, T_new = {Tnew:.3e}")

        #modData['lnT'][zones[a]] = str(np.log(T))
        modData['lnT'][zone] = f"{np.log(T):.16E}"
    return modData
            
# -------------------------------------------------------------------------------------------------------------
# Grealty improved from previous version... 
def readModfile(filename):
    
    params          = {}
    lines           = ih.getTxtArrays(filename)           # List of lines of text...
    startheaderInd  = ih.getTxtIndex(lines, 'lnd')        # Index where header data begins...
    headerVar       = lines[startheaderInd].split()    # Array containing header names...
    numVar          = len(headerVar)                   # Total number of header variables...
    startcomments   = lines[0:startheaderInd-1]        # Gets the inital summary comments...
    endheaderInd    = ih.getTxtIndex(lines, 'previous') - 3  # this will get line index of "previous model"? 
    lastComments    = lines[endheaderInd:len(lines)]      # Gets last lines of summary data        
    
    # Saves useful summary variables... currently misses the "error" and "num_tries" lines...
    # Not important right now..
    for a in range(4, startheaderInd-3):
        var, val = lines[a].split()
        params[var] = val
    
    # Main Focus: Gets lines containing the data we want to save..
    mainData = lines[startheaderInd+1: endheaderInd-1]
    
    # Create new array structure of data to easily save...
    newMainData = []
    for a in range(0, len(mainData)):
        varVal = mainData[a].split()[1:]
        newMainData.append(np.array(varVal))
    newMainData = np.array(newMainData)
    
    # Saves data into dict() params... 
    for a in range(0, numVar):
        datavar = headerVar[a]
        params[datavar] = newMainData[:,a]
            
    # Saving data into dict()
    params['headerVariables'] = lines[startheaderInd]
    params['Index']           = np.arange(1,len(newMainData)+1) # create array for indices
    params['startComments']   = startcomments
    params['lastComments']    = lastComments
    
    return params

def writeModfile(data, newfile):
    t = open(newfile, "w")
    
    # Writing the begining summary comments...
    for a in range(0, len(data['startComments'])):
        t.write(data['startComments'][a] + "\n")
    t.write("\n")
    t.write(data['headerVariables'] + "\n")
    
    # Writing main array columns..
    writer    = csv.writer(t, delimiter = ' ')
    headerVar = data['headerVariables'].split()
    outData   = []
    for a in range(0, len(headerVar)):
        outData.append(data[headerVar[a]])
    
    mainDat = pd.DataFrame(outData, columns = data['Index'])
    mainDat = mainDat.transpose().to_string(header=False)
    t.write(mainDat)
    
    # Writing ending summary comments
    t.write("\n\n")
    for a in range(0, len(data['lastComments'])):
        t.write(data['lastComments'][a] + "\n")    
    t.close()

def getIsotopes(modData, zones):
    presentIsotopes = [iso for iso in const.commonIsotopes if iso in modData]
    isotopeData = {iso: modData[iso][zones] for iso in presentIsotopes}

    # Use for debugging, make sure we are depositing in the correct radial zones...
    #isotopeData['lnR']   = modData['lnR'][zones]
    #isotopeData['Index'] = zones
    
    return isotopeData

    