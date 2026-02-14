import numpy as np
import os
import CEconstants as const

"""
    Description: Set up runtime parameters you can choose at the start of run. You can choose
                 to run predfined MESA models from MESA, type of drag or jet feedback 
"""
class RunTimeParameters:
    def __init__(self):
        """
            Choose donor and accretor defining radi and mass
        """
        self.starType         = 0         # 0=User. Consider preset
        self.starModelPath    = '/home/david101/Desktop/Research/CommonEnvelope/mesa-inspiral/StarModels'
        self.starProfileModel = '/Fragos/12MsunV1Profile.data'
        self.starModModel     = '/Fragos/12MsunV1.mod'
        
        self.accretorType  = 'NS'      # NS / BH / WD etc.
        self.accretorRadi  = 1.2e6     # Bronner=3.1*Rsun,Griencher = 1.2e6[12km]
        self.accretorMass  = 1.4       # [Msun]

        """
            Choose simulation paramters...
        """
        self.mesaRun      = False        # Run With MESA?
        self.a0           = 0.7          # x donor max radius
        self.Rends        = 0.1*const.Rsun_cgs # Sepeartion Termination
        self.dtType       = 1.e-3*const.yr_cgs     # Timestep, 0 = Uses "secular" timestep...
        self.includeDrag  = True         # Include Drag in MESA
        self.includeJet   = False        # Includes Jet in MESA

        """
            Choose drag and/or jet feedback and revelant parameters
        """
        self.dragMethod   = 'Macleod1' # Macleod2019/2015,Ostriker1999,Bronner2023,Fragos2019
        self.dragKernel   = 'Gauss'   # Uniform, Gauss, Parabolic, Cosine, Quatric' etc.
        self.freeCd       = 0.1       # Drag Coefficient Free Parameter for Bronner
        self.freeCh       = 0.5         # freeCh*Ra for Upper and lower heating limits
        self.jetMethod    = 1         # 1 = Grichener 
        self.jetKernel    = 'Uniform' # Uniform, Gauss, Parabolic, Cosine, Quatric' etc.
        self.jetEff       = 0.4       # Jet Efficiency

        """
            Choose to call MESA Hooks
        """
        self.gravFlag     = False     # Include Graviational Potential hook routine
        self.rotateFlag   = False     # Include Rotation NEED TODO
        self.windFlag     = False     # Includes wind mass hook routine

        """
            Threads, I/O output frequecy, restarts, and debug helpers
        """
        self.numThreads   = 16        # Number of Threads
        self.saveFreq     = 10        # Save Freq. of data, (~10 or 50 works for longer runs)
        self.ckptMax      = 1         # Max checkpoint files
        self.zoneLimit    = 1         # Number of outer zones we do not deposit energy
        self.nmax          = int(10e9)
        self.restart      = 0         # TODO

        """
            Can change these things as needed
        """
        self.inlistFile       = 'inlist_project'  # Name of MESA inlist
        self.pathMesa         = './MesaWorkDir'   # Path to Mesa working directory
        self.saveDataPath     = '../..'
        self.saveDonorName    = 'RSGDonor1'
        self.saveAccretorName = 'NSAccretor1'
        
        """
            File and directory names. Don't modify. These are set during set up 
        """
        self.profile        = ''
        self.modFile        = ''
        self.finalMod       = ''
        self.pathDonor      = ''
        self.pathAccretor   = ''
        self.terminalOutput   = 'MesaTerminalOutPut' # MESA ongoing summary output
        self.inspiralPath   = os.getcwd()

        