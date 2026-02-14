import sys
sys.path.append('../ExtraFiles/')
sys.path.append('../../RunData/') # Gets access to ReadInspiralData: Needed for restarts
sys.path.append('../')
import pandas as pd
import numpy as np
import time
import os
import mesa_reader as mr

# Needed for MESA integration
import MESAGen as          mg
import EnergyInjections as ei
import EnergyHelper as     mh
import InlistHelper as     ih
import utility as          rd    # Needed for Restart Runs [Need to Work]
import CEIO as ceIO
import CEphysics as cePhys
import CEParameters as ceParam
import CEconstants  as const

class InspiralParameters:
    def __init__(self, separation, donor_profile, accretor_mass):
        Rs, rhos, _ms, css, gams, Ebs = donor_profile
        self.m1 = accretor_mass * const.Msun_cgs
        self.Rs, self.rhos, self._ms, self.css, self.gams, self.Ebs = Rs, rhos, _ms, css, gams, Ebs
        self.a = separation*const.Rsun_cgs * self.Rs[-1]

        self.m2 = np.interp(separation, Rs[1:], _ms)
        self.rho = np.interp(separation, Rs, rhos)
        self.cs  = np.interp(separation, Rs, css)
        self.gamma = np.interp(separation, Rs, gams)
        self.Eb = np.interp(separation, Rs, Ebs)
        self.v = np.sqrt(const.G_cgs * (self.m1 + self.m2) / separation)
        self.P = np.sqrt(4 * np.pi**2 / (const.G_cgs * (self.m1 + self.m2)) * separation**3)
        self.Omega = 2 * np.pi / self.P
        self.Eorb = -0.5 * const.G_cgs * self.m1 * self.m2 / separation
        self.Ra = cePhys.getRaccHL(self.m1, self.v)
        self.Mach = self.v / self.cs
        return




# ---------------- SETUP------------------------------
# Could probably refactor to make it nicer...
runParameters = ceParam.RunTimeParameters()
pathAccretor = ceIO.createDir('../../RunData', 'NSAccretorTEST')
pathDonor    = ceIO.createDir('../../RunData', 'RSGDonorTEST')
starModelPath    = '/home/david101/Desktop/Research/CommonEnvelope/mesa-inspiral/StarModels'
starProfileModel = '/Fragos/12MsunV1Profile.data'
starModModel     = '/Fragos/12MsunV1.mod'

# Sets up general profile and mod files to match star
starDat = mr.MesaData(starModelPath + starProfileModel)
starMass = round(starDat.header_data['star_mass'], 2)
runParameters.profile, runParameters.modFile, runParameters.finalMod = [str(starMass) + 'Msun' + string for string in ['Profile.data', 'Input.mod', 'Output.mod']]

# Setting Up name prameters
runParameters.pathDonor    = pathDonor +'/'
runParameters.pathAccretor = pathAccretor + '/'

# Copyinig starting mod and mesa profiles
mg.writeModfile(mg.readModfile(starModelPath + starModModel), runParameters.modFile)
ih.copyFile(starModelPath + starProfileModel, runParameters.inspiralPath, runParameters.profile)
# ---------------------------------------------------
# Sets up MESA stuff.
ih.setupMesaInlistPaths(runParameters)     # MESA Working Dir reads current Dir inlist...
ih.setupInspiralInlistPaths(runParameters) # Modify inlist_project setting up file/path names...
ih.setInitialInlistAge(runParameters)      # Modify inlist to match Mesa profile's current age...

donorDat = mg.getMESAData(runParameters.profile)

print(runParameters.accretorMass)
inspiralParm = InspiralParameters(runParameters.a0, donorDat, runParameters.accretorMass)
ceIO.writeInspiralParamters(runParameters, inspiralParm)














