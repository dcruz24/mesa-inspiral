import sys
import numpy as np
import pyMesa as pym        # pyMesa v2.0.0
from subprocess import run
import subprocess
import csv
import pandas as pd
import os
import InlistHelper as ih
import CEconstants as const

# -------------------------------------------------------------------------------------------------------------
# path: path to your MESA workspace
def RunMesa(path, numthreads):
    env = dict(os.environ)
    env['OMP_NUM_THREADS'] = str(numthreads)
    with open('MESA_OUTPUT.txt', 'w') as f:
        run(["./rn"], cwd = path, stdout = f, text = True, env = env)
    return

# Object Container for pyMesa v2.0.0
class pyMesa_context():
    def __init__(self, cacheDir):
        
        self.eos_lib, self.eos_def     = pym.loadMod("eos", cacheDir)
        self.const_lib, self.const_def = pym.loadMod("const", cacheDir)
        self.math_lib, _               = pym.loadMod("math", cacheDir)
        self.chem_lib, self.chem_def   = pym.loadMod("chem", cacheDir)
        
        self.math_lib.math_init()
        ierr=0
        self.const_lib.const_init(pym.MESA_DIR, ierr)
        ierr=0
        self.chem_lib.chem_init('isotopes.data', ierr)
        ierr=0
        self.eos_lib.eos_init(pym.EOSDT_CACHE, True, ierr)
        
        self.eos_inlist = 'eos.inlist'
        with open(self.eos_inlist, 'r') as f:
            print(self.eos_inlist)
            print("eos.inlist contents:")
            print(f.read())#with open(self.eos_inlist,'w') as f:
        #    print('&eos',file=f)
        #    # Set options here
        #    print('/',file=f)
            
        ierr = 0
        self.res = self.eos_lib.alloc_eos_handle_using_inlist(self.eos_inlist, ierr)
        
        if self.res.args['ierr'] != 0:
            raise ValueError("Ierr not zero from alloc eos_handle")
        self.eos_handle = self.res.result
        print('Loaded Modules for eos, const, math, chem...')
        return
        
# =======================================================================================================
# Using pyMesa, calls MESA EOS to get Internal Energy
def getEosResult_IE(context, rho, T, species, chem_ids, net_iso, xa):
    #print('--------- Calling EOS to get IE ---------')
    logRho  = np.log10(rho)
    logT    = np.log10(T)
    #print('Current Log10 T:    {}'.format(logT))
    #print('Current Log10 Rho:  {}'.format(logRho))
    nres   = context.eos_def.num_eos_basic_results
    res    = np.zeros(nres)
    d_dlnRho = np.zeros(nres)
    d_dlnT = np.zeros(nres)
    d_dxa  = np.zeros((nres,species))
    ierr   = 0
    eos_result = context.eos_lib.eosdt_get(context.eos_handle, species, chem_ids, net_iso, xa,
                                     10**logRho, logRho, 10**logT, logT, res,
                                     d_dlnRho, d_dlnT, d_dxa, ierr)
    if(eos_result.args['ierr'] == -1):
        print('Error with getting interntal energy... context.eos_lib.eosdt_get(...)')
        return
    return eos_result
    
# =======================================================================================================
# Using pyMesa, calls Mesa EOS to get temperature...
# **Note** 
#  For the test case using T and rho value and calling getEosResult_IE() to get the internal energy...
#  Then using the same rho and IE value to calculate the same temperature...
#  Increasing the max_iter seems to worsen the result

def getEosResult_Temperature(context, rho, T, species, chem_ids, net_iso, xa, oldIE, IE, i_lnE):
   # print('--------- Calling EOS to get Temperature ---------')
    logRho, logT = np.log10(rho), np.log10(T)
    #print('Current Log10 T:   {}'.format(logT))
    #print('Current Log10 Rho: {}'.format(logRho))
    #print('species:             ', species)
    #print('xa:             ', sum(xa))
    
    arg_not_provided = -9e99
    which_other      = 2     # refers to index defined in eos_def. IE=2     #i_lnE + 1
    other_value      = np.log(IE)
    logT_tol         = 1e-8
    other_tol        = 1e-8
    max_iter         = 20
    ierr = eos_calls = 0
    logT_bnd1 = logT_bnd2 = other_at_bnd1 = other_at_bnd2 = arg_not_provided
    logT_guess = logT_result = logT

    deltaE = abs(IE-oldIE)/IE
    #print('Zone IE:     ', oldIE)
    #print('New Zone IE: ', IE)
    #print('DELTA IE:    ', IE-oldIE)
    #logT_guess = logT + logT_shift
    #logT_guess = logT + np.log10(1.1)#* deltaE
    logT_guess = logT
    #print('Guess of log10T: ', logT_guess)
    
    nres = context.eos_def.num_eos_basic_results
    res = d_dlnRho_const_T = d_dlnT_const_Rho = np.zeros(nres)
    d_dxa_const_TRho = np.zeros((nres, species))
    
    eosT_result = context.eos_lib.eosdt_get_t(
        context.eos_handle, species, chem_ids, net_iso, xa,
        logRho, which_other, other_value,
        logT_tol, other_tol, max_iter, logT_guess,
        logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2,    # INS
        logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho,  # INOUT
        d_dxa_const_TRho, eos_calls, ierr)
    
    newLogT = eosT_result.args['logt_result']
    #print('Final Resulted  Log10T: ', newLogT)
    
    # For debugging to original temperature...
    absError = abs(newLogT- logT)
    relError = (newLogT- logT) / logT
    #print('Resulting Log10T:      {}'.format(newLogT))
    #print('EOS Calls:             {}'.format(eosT_result.args['eos_calls']))
    #print('WHAT: ', eosT_result)
    #print('Abs error:             {}'.format(absError ))
    #print('Rel error:             {}'.format(relError))
    #print('% error:               {}'.format( relError * 100))

    if(eosT_result.args['ierr'] == -1):
        print('Error getting Temperature...context.eos_lib.eosdt_get_t()')
        return
        
    return 10**newLogT
# =======================================================================================================

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

    