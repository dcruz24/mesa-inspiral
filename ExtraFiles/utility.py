import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
import pandas as pd
import time
import mesa_reader as mr

G    = 6.6743e-8
c    = 2.99792458e10     # [cm/sec]
Msun = 1.9891e33      # [g]
Rsun = 6.957e10       # [cm]
pc   = 3.08567758128e18 # [cm]
kpc  = 1.e3*pc
Mpc  = 1.e6*pc
au   = 1.49598073e13   # [cm]
yr   = 3.1556e7       # [secs]
Kyr  = 1.e3*yr
Gyr  = 1.e9*yr
# ================================= Contains the following ================================================ 
# Class: readInspiralData
# Class: readTwoBodyInspiralData: Needs to be worked on...
# Class: EnvelopeData
# ---------------------------- For Grientcher Paper ------------------------------------------
#  get_CompareChi()        : Calculates ratio of envelope density at the time the companion is
#                            to the same locaiton in the unperturbed case 
#  getZeta()               : Calculates ration of accretion onto the companion at the time
#                            it is in the envelope to the unperturbed case
# -------------------------------------------------------------------------------------------

class readInspiralData:
    def __init__(self, starPath, donorPath):
        data = np.loadtxt(starPath + 'ceMain', dtype = float, unpack=True)
        self.tl        = data[0] # Time             [yr]
        self.a         = data[1] # Seperation       [Rsun]
        self.m1        = data[2] # Donor Star       [Msun]
        self.m2        = data[3] # Companion Star   [Msun]
        self.dmbs      = data[4]
        self.omn       = data[5]
        self.Edrag     = data[6]
        self.Ejet      = data[7]
        self.v_        = data[8]
        self.P         = data[9]
        self.Rho       = data[10]
        self.Ra        = data[11]
        self.Csound    = data[12]
        self.Gamma     = data[13]
        self.Mach      = data[14]
        self.Eorb      = data[15]
        self.Ebind     = data[16]
        self.Fd        = data[17]
        self.egrad     = data[18]
        self.th        = data[19]
        self.RaccUpper = self.a + .5*self.Ra
        self.RaccLower = self.a - .5*self.Ra
        self.pathAccretor = starPath
        self.pathDonor    = donorPath
        
# Needs to be developed....
class read2BodyData:
    def __init__(self, starPath, donorPath):
        data = np.loadtxt(starPath + 'ce2BODYMain', dtype = float, unpack=True)
        self.tl    = data[0]
        self.a     = data[1]
        self.m1    = data[2]
        self.m2    = data[3]
        self.dmbs  = data[4]
        self.Omn   = data[5]
        self.Edrag = data[6]
        self.Ejet  = data[7]
        self.x1, self.y1, self.z1    = data[8],data[9],data[10]
        self.x2, self.y2, self.z2    = data[11],data[12],data[13]
        self.vx1, self.vy1, self.vz1 = data[14],data[15],data[16]
        self.vx2, self.vy2, self.vz2 = data[17],data[18],data[19]
        self.xs,self.ys,self.zs      = self.x1-self.x2,self.y1-self.y2,self.z1-self.z2
        self.vxs,self.vys,self.vzs   = self.vx1-self.vx2,self.vy1-self.vy2,self.vz1-self.vz2
        self.ra1 = np.sqrt(self.x1**2+self.y1**2+self.z1)
        self.ra2 = np.sqrt(self.x2**2+self.y2**2+self.z2)        
        self.rs   = np.sqrt(self.xs**2+self.ys**2+self.zs**2)
        self.vs   = np.sqrt(self.vxs**2+self.vys**2+self.vzs**2)
        self.v_   = np.sqrt(G*(self.m2+self.m1)/self.rs)
        self.Racc = 2.0*G*self.m1/(self.v_*self.v_)
        self.pathAccretor = starPath
        self.pathDonor    = donorPath
        return

# Calculates the density ratio of the primary's envelope 
# as the companion enters to the first primary mesa profile 
# (Analysis for Grientcher Paper...)
def getCompareChi(ns1, initMesa, filename):
    chi, tl = [], []
    starInit_rho = 10**initMesa.logRho[::-1]
    starInit_Rs  = 10**initMesa.logR[::-1]*Rsun
    
    print('Needed Iteration: ', len(ns1.ts))    
    for a in range(0, len(ns1.ts)):
        if(a % 1 ==0):
            file = ns1.pathStar2 + '{}_'.format(a) + filename
            df       = pd.read_csv(file,encoding="latin-1",delim_whitespace=True,skip_blank_lines=True,skiprows=5)
            Rs       = 10.**df['logR'].to_numpy()[::-1]*Rsun
            rhos     = 10.**df['logRho'].to_numpy()[::-1]
            
            rho_init = np.interp(ns1.rs[a]*Rsun, starInit_Rs, starInit_rho)  # Density at the intial mesa profile
            rho_cur  = np.interp(ns1.rs[a]*Rsun, Rs, rhos)                   # Density Currently
            tl.append(ns1.ts[a])
            chi.append(rhoj/rho0)
    return np.array(tl), np.array(chi)

# Compares the bondi rate in the perturbed envelope (our run) to the unperturbed envelope(inital mesa profile)
# (Analysis Grienthcer Paper..., will be different since we evolve hydrodnamically opposed to hydrodstatic
def getZeta(ns1,initMesa):
    
    Racc = lambda M, v: 2.*G*M/(v*v)                             # Accretion Radius
    Mdot_HL = lambda M, rho, v: np.pi*Racc(M,v)*Racc(M,v)*rho*v  # Bondi Hoyle Lyttleon 

    radi  = 10**initMesa.logR[::-1]*Rsun
    rho   = 10**initMesa.logRho[::-1]
    Macc  = ns1.dmbn # Mass Accretion
    eta   = 0.1       # Jet Efficeny
    Mblst = []
    m11   = 1.4*Msun
    
    for a in range(0, len(ns1.dmbn)):
        rhoi =np.interp(ns1.rs[a]*Rsun, radi,rho)
        _v = np.sqrt(G*(m11 + ns1._ms[a]*Msun)/(ns1.rs[a]*Rsun))        
        Mbondi = Mdot_HL(m11,rhoi,_v)
        m11 = m11 + Mbondi
        Mblst.append(Mbondi)
    Mbon = np.array(Mblst)
    zeta = eta*(Macc/Mbon)
    return zeta, Macc, Mbon