import numpy as np
from src.physics import basic as ce_basic
from src.config import constants as const

"""
    Contains values from fit functions for drag models for stellar mass and compact objects.
    Papers are referenced...
"""

def chooseDragMethod(runParameters, simParameters):
    m1   = simParameters.m1
    m2   = simParameters.m2
    rho  = simParameters.rho
    v    = simParameters.v
    # Need to check if this is valid for the fit that uses mass ratios
    # Is that using donor total mass or enclosed mass?
    q     = m1/m2
    Ra    = simParameters.Ra
    Mach  = simParameters.Mach
    gamma = simParameters.gamma
    eGrad = simParameters.epGrad
    fitType = runParameters.dragMethod

    # Need to see how different his is from using a static defined accretor radius?
    Rs_Ra = 1.0
    
    if (fitType == 'Macleod1'):
        Mdot, Fd = getAcc2019Drag(m1, rho,v, q, Mach, gamma, Rs_Ra)
    #elif (fitType == 'Macleod2'):
    #    Mdot, Fd = cePhys.getAcc2019Drag(m1, rho, v, Ra, eGrad)
    #elif (fitType == 'Ostriker'):
    #    Fd    = cePhys.getOst1999Drag(m1, Mach, rho, v, Rs_Ra, dt)
    #    Mdot = cePhys.Mdot_HL(m1, rho, v)
    #elif(fitType == 'Bronner'):
    #    I      = cePhys.getCoulombLog(Macg, R, Rs_Ra)
    #    _Fd    = cePhys.getBronnerDrag(m1, _cs, _Ma,_rho, freeCd, _R, _v, I)
    #    dmb_dt = cePhys.Mdot_HL(m1, _rho, _v)
    #elif (fitType == 'Hybrid'):
    #    Mdot, Fd   = cePhys.getHybirdDrag(m1, _Ma, _rho, _v, a, _Ra, _cs, rMin, dt, epGrad)
    else:
        print('No Drag Model Picked')
    dEdrag = -Fd*v
    
    return Mdot, Fd, dEdrag


# Cad43 and Cad53: drag coefficients
def Cad43(q, Ma):
    """
    Gamma = 4/3
    """
    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 = 0.8169,-9.9784,0.1382,-0.4803,46.9755,-0.3330,0.6713,-4.1620,-58.9379,0.0379
    d1,d2,d3,d4,d5,d6,d7,d8,d9,d10 = 0.5510,0.4502,-0.6741,0.5946,-14.95,0.3159,-0.0203,-1.7043,30.4494,-0.0309
    log10Ca = a1 + a2*q + a3*Ma + a4*q*Ma + a5*q*q + a6*Ma*Ma + a7*q*Ma*Ma + a8*q*q*Ma + a9*q**3 + a10*Ma**3
    log10Cd = d1 + d2*q + d3*Ma + d4*q*Ma + d5*q*q + d6*Ma*Ma + d7*q*Ma*Ma + d8*q*q*Ma + d9*q**3 + d10*Ma**3
    return 10**log10Ca, 10**log10Cd

def Cad53(q,Ma):
    """
    Gamma = 5/3
    """
    a1,a2,a3,a4,a5,a6 = 0.9184,-0.9619,-1.2057,1.2247,-2.480,0.1150233
    d1,d2,d3,d4,d5,d6 = -0.1552,-3.0323,0.2756,0.1976,1.4186,-0.0092
    log10Ca = a1 + a2*q + a3*Ma + a4*q*Ma + a5*q*q + a6*Ma*Ma
    log10Cd = d1 + d2*q + d3*Ma + d4*q*Ma + d5*q*q + d6*Ma*Ma
    Ca, Cd  = 10.**log10Ca, 10.**log10Cd
    return Ca,Cd

# ----------------------------------------------------------------------------------------------------
def getAcc2019Drag(m,rho,v,q,Ma,gam, Rs_Ra, fit=False):
    """
    De, Macleod et. 2019 APJ 897 130, Appendix A
    Note: Cad43' coefficents best applicable for systems where donor Mass > 10Msun and Accretor is a black hole
    q: Mass Ratio
    Ma: Mach Number
    """
    dMHL = ce_basic.getMdotHL(m, rho, v)
    FdHL = ce_basic.getFdHL(m, rho, v)
    Ca43,Cd43  = Cad43(q,Ma)
    Ca53,Cd53  = Cad53(q,Ma)
    
    # Radiative Pressure Dominated State?
    if(gam < 4./3.): 
        w43, w53 = 1., 0.        
    # A mixture bewteen both radiative and gas presusure state?
    else: 
        w43, w53 = 1.-3.*(gam-4./3.), 1.-3.*(5./3.-gam)
        
    Ca,Cd = w43*Ca43 + w53*Ca53, w43*Cd43 + w53*Cd53
    
    # See Figure 7th
    if(Ca > 1.):
        Ca = 1.
    if(Cd > 10.):
        Cd = 10.
        
    if(Ma >= 1.7):
        alphaM = 0.33
    elif(Ma < 1.7):
        alphaM = (0.62*Ma) - 0.72
        
    #Ca = Ca*(Rs_Ra/0.05)**alphaM # Section 4.0
    Ca = Ca*(Rs_Ra/0.05)**0.33 # Section 4.0
    dM, Fd = Ca*dMHL, Cd*FdHL
    return dM,Fd

# ----------------------------------------------------------------------------------------------------
def getNS2015Drag(m,rho,v,Rs,ep,fit=False):
    """
    Macleod 2015 APJ 803:41    (3D studies wind tunnels)
    Macelod 2015 APJL 798:L19 (Accretion during CE for NS)
    Appendix B

    Note: Best applicable fo MS embedded inside CE
          Caution when applied for compact objects
    ep: 0.3-3.0
    m:   Mass of Accretor
    rho: local envelope density
    Ra:  Accretion Radius
    ep:  density gradient range
    """

    Racc = ce_basic.getRaccHL(m,v)
    
    # need a better way to seperate between MS and compact object
    # Rs/Ra = 0.05
    if (m >  3*const.Msun_cgs):
        c1,c2,c3    = 1.98255197, -1.33691133, 0.62963326
        b1, b2, b3  = 0.06549409, -6.87212261, 27.02371844
        a1,a2,a3,a4 = -1.65171739, 1.49979486, 0.10226072, 3.93190671
    else:
        # Rs/Ra = 0.01
        c1,c2,c3    = 1.9179196, -1.5281498, 0.75992092
        b1, b2, b3 = 1.86818916e-2, -6.42396570, 34.0135578
        a1,a2,a3,a4 = -2.14034214, 1.94694764, 1.19007536, 1.05762477

    FdragFit    = c1 + (c2*ep) + (c3*ep**2)
    mdotDragFit = 10**(a1 + a2/(1+(a3*ep) + (a4*ep**2) ))
    LdotDragFit = (b1**ep)  / (1+b2*ep + b3*ep**2)

    if fit == True:
        return mdotDragFit, FdragFit, LdotDragFit
            
    FdHL  = ce_basic.getFdHL(m,rho,v)
    Fdrag = FdragFit * FdHL
    
    mdotHL      = ce_basic.getMdotHL(m,rho,v)
    mdotDrag    = mdotDragFit * mdotHL

    LdotDrag    = mdotHL*Racc*v*LdotDragFit
    
    return mdotDrag, Fdrag, LdotDrag
    
# ----------------------------------------------------------------------------------------------------
def getCoulombLog(Ma,r, rmin):
    """
    Kim 2010 APJ 725:1069-1081
    Kim/Kim 2007, APJ 665:432
    Bronner 2023, Ostriker 1999
    Ma:   Mach number 
    r:    Radius of Orbit.
    rmin: Size of the Accretor
    """
    if (Ma < 1.0):
        I = 0.7706*np.log((1+Ma)/(1.0004-0.9185*Ma)) - 1.473 * Ma
    elif (Ma >= 1.0 and Ma < 4.4):
        I = np.log(330.0*(r/rmin) * ( ((Ma-0.71)**5.72)/ (Ma**9.58) ) )
    elif  (Ma >= 4.4):
        I = np.log( (r/rmin)/ (0.11*Ma + 1.65) )
    return I

# ----------------------------------------------------------------------------------------------------
def getBronnerDrag(m2,cs, Ma,rho, Cd, r, v, I):
    """
    Bronner 2023
    """
    Beta  = (const.G_cgs*m2)/ ((cs**2)*r)
    nBeta = Beta/((Ma**2)-1)
    rhoNL = rho*(1+ ( (0.46*(Beta**1.1))/ (((Ma**2)-1)**0.11)) )

    if (nBeta > 0.1 and Ma > 1.01):
        Fd = Cd*(4*np.pi*rhoNL*(const.G_cgs*m2)**2)/(v**2) * (0.7/(nBeta**0.5))
    else:
        Fd = Cd*(4*np.pi*rho*(const.G_cgs*m2)**2)/(v**2) * (I)
    return Fd

# ----------------------------------------------------------------------------------------------------
def getOst1999Drag(m2, Ma, rho, v, rmin, dt):
    """
    Ostriker 1999, APJ 513 252
    """
    if (Ma > 1): # Supersonic Regime
        I = 0.5 * np.log(1-(1/Ma**2)) + np.log(v*dt/rmin)
        # Note... If Mach >> 1, I = np.log(v*dt/rmin)
    else: # Sub-sonic...
        I = 0.5 * np.log((1+Ma)/(1-Ma)) - Ma
        # Note... If Mach <<< 1, I=M**3/3
    Fd = (4.0*np.pi*(const.G_cgs*m2)**2 * rho) / (v**2)
    Fd = Fd * I
    return Fd

# ----------------------------------------------------------------------------------------------------
def getHybirdDrag(m2, Ma,rho, v, a, Ra, cs, rmin, dt, epGrad):
    """
    Create a Hybird drag Formualism for different Mach regions
    Ma > 1.1 Apply Macleod 2015
    Ma < 0.9 Apply Ostriker 1999
    """
    if(Ma > 1.1):
        dM, Fd = getNS2015Drag(m2, rho, v, Ra, epGrad)
    elif (Ma < 1.1):
        Fd = getOst1999Drag(m2, Ma, rho, v, rmin, dt)
        dM = ce_basic.getMdotHL(m2, rho, v)
    return dM, Fd

