import numpy as np
import CEconstants  as const
import CEphysics as cephys
class InspiralParameters:
    """
        Contains the local variables where the accretor is located at
        the envelope of the donor
    """
    def __init__(self, parameters, donorModel, accretorModel):

        # These are setting up the initial conditions depending on where the accretor is 
        # at the location of envelope of the donor
        self.a     = parameters.a0  * donorModel.R[-1]
        self.m1    = accretorModel.m
        self.m2    = np.interp(self.a, donorModel.R, donorModel.Menc)
        self.rho   = np.interp(self.a, donorModel.R, donorModel.rho)
        self.cs    = np.interp(self.a, donorModel.R, donorModel.cs)
        self.gamma = np.interp(self.a, donorModel.R, donorModel.gamma)
        self.Eb    = np.interp(self.a, donorModel.R, donorModel.Ebind)
        self.v     = np.sqrt(const.G_cgs * (self.m1 + self.m2) / self.a)
        self.P     = np.sqrt(4.0 * np.pi**2.0 / (const.G_cgs * (self.m1 + self.m2)) * self.a**3)
        self.Omega = 2.0 * np.pi / self.P
        self.Eorb  = -0.5 * const.G_cgs * self.m1 * self.m2 / self.a
        self.Ra    = cephys.getRaccHL(self.m1, self.v)
        self.Mach  = self.v / self.cs
        self.Edrag = 0.0  # [ergs]
        self.Ejet  = 0.0  # [ergs]
        self.Mdot  = 0.0
        self.th    = 0.0
        self.t1    = 0.0
        self.t2    = 0.0
        self.th    = 0.0
        # Setting initial parameters
        self.f      = 50.              # Frequency
        self.Om     = 2.*np.pi/self.P  # Angular Velocity
        
        # Other values that would be usefult to be defined throughout
        self.dEdrag  = 0.0
        self.dEjet   = 0.0
        self.da_dt   = 0.0
        self.dth_dt  = 0.0
        self.epGrad  = 0.0


        self.runFlag       = True
        self.mesaFlag      = True
        self.MESA_maxAge  = 0.0
        self.MESA_currAge = 0.0
        self.wallclock     = 0.0

        self.count         = 0
        self.savenum       = 0
        self.chkpoint      = 0

        return

    # Helper function whhen we get the local values interpolating as the accretor inspirals
    # within the envelope
    def getLocalVariable(self, donorModel):
        self.m2    = np.interp(self.a, donorModel.R, donorModel.Menc)
        self.rho   = np.interp(self.a, donorModel.R, donorModel.rho)
        self.v     = np.sqrt(const.G_cgs * (self.m1 + self.m2) / self.a)
        self.P     = np.sqrt(4.0 * np.pi**2.0 / (const.G_cgs * (self.m1 + self.m2)) * self.a**3)
        self.Ra    = cephys.getRaccHL(self.m1, self.v)
        self.cs    = np.interp(self.a, donorModel.R, donorModel.cs)
        self.gamma = np.interp(self.a, donorModel.R, donorModel.gamma)
        self.Mach  = self.v / self.cs
        self.Eorb  = -0.5 * const.G_cgs * self.m1 * self.m2 / self.a
        self.Eb    = np.interp(self.a, donorModel.R, donorModel.Ebind)
        return

    def calculateDragRadius(self, runParameters):
        self.RupperDrag = self.a + runParameters.freeCh*self.Ra
        self.RlowerDrag = self.a - runParameters.freeCh*self.Ra
        return

    def calculateJetRadius(self, runParameters, donorStar):
        self.RupperJet = 1.0*donorStar.R[-1]
        self.RlowerJet = self.a
        return


