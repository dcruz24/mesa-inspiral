import sys
import numpy as np
import pyMesa as pym        # pyMesa v2.0.0
from subprocess import run
import subprocess
import os
from src.io import inlist as ih
import mesa_reader as mr
from src.config import constants as const
import gfort2py
import pandas as pd

"""
    Description: Contains general I/O functions to help with writing out data files and 
                 simulation parameters to reproduce.
"""


def createDir(storDir, subDir):
    """
        Helper function to create directory to save for a run...
    """
    fullPath = os.path.join(storDir, subDir)        
    if os.path.exists(fullPath) and os.path.isdir(fullPath):
        print(f"'{subDir}' already exists within '{storDir}'.")
    else:
        os.makedirs(storDir + '/' + subDir)
        print('Created ' + subDir)
    return storDir + '/' + subDir


def saveInspiralFile(fileName, filePath, count, string):
    """
        Saves the inspiral data file?
    """
    subprocess.run('cp {} '.format(fileName) + filePath + '{}'.format('{}'.format(count) + string + fileName)
              ,shell = True)
    return


def saveCEData(simParameters, runParameters):
    if (simParameters.count % runParameters.saveFreq ==0):
        row = {
            'Time': simParameters.t2 / const.yr_cgs,
            'a':  simParameters.a / const.Rsun_cgs,
            'M1': simParameters.m1 / const.Msun_cgs,
            'M2': simParameters.m2 / const.Msun_cgs,
            'v': simParameters.v,
            'Eorb': simParameters.Eorb,
            'Edrag': simParameters.Edrag,
            'Ejet': simParameters.Ejet,
            }
        if hasattr(simParameters, 'x'):
            row.update({
                'x': simParameters.x / const.Rsun_cgs,
                'y': simParameters.y / const.Rsun_cgs,
                'vx': simParameters.vx,
                'vy': simParameters.vy,
            })
        save1dRowData(row, 'simpleInspiralData.csv', runParameters)
    return
    
def save1dRowData(row_dict, fileName, runParameters):
    filePath = os.path.join(runParameters.pathAccretor, fileName)
    row_df = pd.DataFrame([row_dict])
    header_flag = not os.path.exists(filePath)
    row_df.to_csv(filePath, mode='a', header=header_flag, index=False)
    return
    

def save1dData(data, fileName, headerName, runParameters):
    np.savetxt(runParameters.pathAccretor + fileName, data, fmt='%f', header = headerName)
    return


def printSimulationOutput(runParameters, simParameters, donorStar):
    if(simParameters.count % runParameters.saveFreq == 0):
        table_data = [
            ['Iteration:               ', simParameters.count], 
            ['SaveNum:                 ', simParameters.savenum], 
            ['Time                 [yr]', simParameters.t1/const.yr_cgs],
            ['Donor Mass         [Msun]', simParameters.m2/const.Msun_cgs],
            ['Accretor Mas       [Msun]', simParameters.m1/const.Msun_cgs],
            ['Donors Max Radius  [Rsun]', donorStar.R[-1]/const.Rsun_cgs],
            ['Accretor Position  [Rsun]', simParameters.a/const.Rsun_cgs],
            ['Wall Clock               ', simParameters.wallclock]
            ]
        
        print("==========================================================")
        for row in table_data:
            print("{: >10} | {: >10}".format(*row))
    return

# Need to create function that saves appending data and gathers it to 


def checkPointSave(runParameters, simParameters):
    if (simParameters.chkpoint < runParameters.ckptMax):
        saveInspiralFile(runParameters.profile, 
                         runParameters.pathDonor, simParameters.chkpoint, '_chkpt_')
        saveInspiralFile(runParameters.modFile, 
                         runParameters.pathDonor, simParameters.chkpoint, '_chkpt_')
        simParameters.chkpoint += 1
    elif(simParameters.chkpoint > runParameters.ckptMax):
        simParameters.chkpoint = 0
    return

def restartPointSave(runParameters):
    print('Saving Last Mod & Profile')
    saveInspiralFile(runParameters.profile, runParameters.pathDonor, runParameters.restart, '_restartPt_')   # Restart File
    saveInspiralFile(runParameters.modFile, runParameters.pathDonor, runParameters.restart, '_restartPt_')   # Restart File
    return

def writeInspiralParamters(runParameters, simParameters, donorStar):
    """
    Writes out a log containing CE parameters and relevant directory information
    using the InspiralParameters class.
    """
    m10    = simParameters.m1
    a0     = simParameters.a
    t0     = simParameters.t1
    Om0    = simParameters.Omega
    m20    = simParameters.m2
    Edrag0 = simParameters.Edrag
    Ejet0  = simParameters.Ejet
    outR   = donorStar.R[-1]

    t = open('InspiralRuntimeParameters', "w")

    t.write("***********************************************************************\n")
    t.write('-----------System & Library Information--------------\n')
    t.write('Python version:                 {}\n'.format(sys.version))
    t.write('pyMesa Directory:               {}\n'.format(os.path.dirname(pym.__file__)))
    t.write('gFort2py Directory:             {}\n'.format(os.path.dirname(gfort2py.__file__)))
    t.write('MESA_DIR:                       {}\n'.format(pym.MESA_DIR))
    t.write('MESASDK_ROOT:                   {}\n'.format(pym.MESASDK_ROOT))
    t.write('OMP_NUM_THREADS:                {}\n'.format(runParameters.numThreads))
    t.write("\n\n")

    t.write('----------- Mesa Inlist File/Path Names --------------\n')
    t.write('Base Profile:                   {}\n'.format(runParameters.profile))
    t.write('Input Mod File:                 {}\n'.format(runParameters.modFile))
    t.write('OutPut Mod File:                {}\n'.format(runParameters.finalMod))
    t.write('Mesa Inlist:                    {}\n'.format(runParameters.inlistFile))
    t.write('Saving Run Data to:             {}\n'.format(runParameters.pathAccretor))
    t.write('Saving Mesa Profiles:           {}\n'.format(runParameters.pathDonor))
    t.write('Mesa Terminal File:             {}\n'.format(runParameters.terminalOutput))
    t.write('Inspiral Working Dir:           {}\n'.format(runParameters.inspiralPath))
    t.write('Mesa Working Dir:               {}\n'.format(runParameters.pathMesa))
    t.write("\n\n")

    t.write('----------- Inspiral Initial Parameters --------------\n')
    t.write('Restart Run:                    {}\n'.format(runParameters.restart))
    t.write('Star Type:                      {}\n'.format(runParameters.starType))
    t.write('Donor Mass m1:       [Msun]     {}\n'.format(round(m20 / const.Msun_cgs, 3)))
    t.write('Donor Radius:        [Rsun]     {}\n\n'.format(round(outR / const.Rsun_cgs, 3)))

    t.write('Accretor Mass m2:    [Msun]     {}\n'.format(round(m10 / const.Msun_cgs, 2)))
    t.write('Accretor Radius:     [km]       {}\n\n'.format(runParameters.accretorRadi / 1e5))

    t.write('dt:                  [yrs]      {}\n'.format(runParameters.dtType / const.yr_cgs))
    t.write('Orbit Integrator:               {}\n'.format(runParameters.orbitIntegrator))
    t.write('Initial Time:        [yrs]      {}\n'.format(t0 / const.yr_cgs))
    t.write('Initial Separation   [Rsun]     {}\n'.format(round(a0 / const.Rsun_cgs, 3)))
    t.write('Initial Drag Energy  [ergs]     {}\n'.format(Edrag0))
    t.write('Initial Jet Energy   [ergs]     {}\n\n'.format(Ejet0))

    t.write('Include MESA:                   {}\n'.format(runParameters.mesaRun))
    t.write('Deposit Drag Energy:            {}\n'.format(runParameters.includeDrag))
    t.write('Deposit Jet Energy:             {}\n'.format(runParameters.includeJet))
    t.write('Include Compact Potential:      {}\n'.format(runParameters.gravFlag))
    t.write('Include Mass Loss:              {}\n'.format(runParameters.windFlag))
    t.write('Include Donor Rotation:         {}\n\n'.format(runParameters.rotateFlag))

    t.write('Drag Description:               {}\n'.format(runParameters.dragMethod))
    t.write('Drag Energy Kernel:             {}\n'.format(runParameters.dragKernel))
    t.write('Cd Drag Parameter:              {}\n'.format(runParameters.freeCd))
    t.write('Ch Drag Heating Parameter:      {}\n\n'.format(runParameters.freeCh))

    t.write('Jet Energy:                     {}\n'.format(runParameters.jetMethod))
    t.write('Jet Energy Kernel:              {}\n'.format(runParameters.jetKernel))
    t.write('jetEff:                         {}\n\n'.format(runParameters.jetEff))

    t.write('Outer Zone Energy Limiter:      {}\n'.format(runParameters.zoneLimit))
    t.write("\n\n")
    t.write('----------- OutputFiles --------------\n')
    t.write('Mesa Profile Save Frequency:    {}\n'.format(runParameters.saveFreq))
    t.write('Max # CheckPoint Mesa Profiles: {}\n'.format(runParameters.ckptMax))
    t.write("\n")
    t.close()
