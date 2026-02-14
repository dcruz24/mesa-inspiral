import numpy as np
import CEconstants  as const


class accretorStar:
    def __init__(self, runParameters):
        self.m    = runParameters.accretorMass * const.Msun_cgs
        self.r    = runParameters.accretorRadi # [km]
        self.Bind = True   # I think this is when calcualting the binding energy?
        return

    # Do this at the start setting arrays of data.
    def loadEOS(self):
        eos = '../eos_tables/SLY4.dat'
        M,Rn,dM,dR,enu,bI,h0 = np.loadtxt(eos,unpack=True)
        self.M_EOS   = M*const.Msun_cgs
        self.Rn_EOS  = Rn/const.R_km_geo*1.e5
        self.dM_EOS  = dM*const.Msun_cgs
        self.dR_EOS  = dR/const.R_km_geo*1.e5
        self.I_EOS   = bI*(M*const.Msun_cgs)**3.*const.G_cgs*const.G_cgs/const.c_cgs**4.0
        self.enu_EOS = enu
        self.h0_EOS  = h0
        return

    # need to feed in Om, and determine if this only the mass of the accretor?
    def getmTOV(self, simParameters):
        Om = simParameters.Om
        Omks  = np.sqrt(const.G_cgs * self.M_EOS / self.Rn_EOS**3)
        mtov1 = self.m +(Om/Omks)**2.0*self.dM_EOS
        self.mTOV  = np.interp(self.m, mtov1, self.M_EOS)
        return

    def getstaticRadius(self):
        self.R = np.interp(self.mTOV, self.M_EOS, self.Rn_EOS)
        return

    def getInertial(self):
        self.I = np.interp(self.mTOV, self.M_EOS, self.I_EOS)
        return
        
    def getDiffRadi(self):
        self.dR = np.interp(self.mTOV,self.M_EOS,self.dR_EOS)
        return
        
    def getdIdM(self):
        self.dI_dM = np.interp(self.mTOV,self.M_EOS,np.gradient(self.I_EOS, self.M_EOS))
        return

    def getEnu(self):
        self.enu = np.interp(self.mTOV, self.M_EOS, self.enu_EOS)
        return
        
    def getH0(self):
        self.h0 = np.interp(self.mTOV, self.M_EOS, self.h0_EOS)
        return
        
    def getFk(self):
        self.fk = 1./(2.*np.pi)*np.sqrt(const.G_cgs*self.mTOV/self.R**3.)
        return
        
    def getAccretorValues(self, simParameters):
        self.getmTOV(simParameters)
        self.getstaticRadius()
        self.getInertial()
        self.getDiffRadi()
        self.getEnu()
        self.getH0()
        self.getFk()
        self.getdIdM()
        self.f = simParameters.Om / (2.0 *np.pi)
        self.eps = self.f / self.fk
        self.eps2 = self.eps * self.eps
        self.Rep  = self.R + self.eps2 * self.dR
        self.Rs_Ra = self.Rep/simParameters.Ra        
        return
        
    def checkBind(self, dMdot):
        if (self.Bind):
            self.hep   = self.eps2*self.h0
            self.phi   = np.sqrt(self.enu* (1.0 + 2.0*self.hep))
            self.dm_dt = self.phi * dMdot
        return

    def getAngMom(self, simParameters):
        Om = simParameters.Om
        denom = (self.I * (np.sqrt(const.G_cgs*self.m*self.Rep) - Om * self.dI_dM))
        self.dOm_dt = self.dm_dt / denom
        return
        