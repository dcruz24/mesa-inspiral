import numpy as np
import CEconstants as const

"""
    Contains common compuatations physics estimations..
"""
def getRaccHL(M, v):
    """
    Hoyle–Lyttleton accretion radius [cm]
    M: accretor mass [g]
    v: relative velocity [cm/s]
    """
    return 2. * const.G_cgs * M / v**2

def getMdotHL(M, rho, v):
    """
    Hoyle–Lyttleton mass accretion rate [g/s]
    rho: ambient density [g/cm^3]
    v: relative velocity [cm/s]
    """
    Racc = getRaccHL(M, v)
    return np.pi * Racc**2 * rho * v

def getFdHL(M, rho, v):
    """
    Hoyle–Lyttleton drag force [dyne]
    """
    return getMdotHL(M, rho, v) * v

def getMdotBHL(M, rho, v, cs):
    """
    Bondi–Hoyle–Lyttleton accretion rate [g/s]
    cs: sound speed [cm/s]
    """
    return 4 * np.pi * const.G_cgs**2 * M**2 * rho / ((cs**2 + v**2)**1.5)

def getRaccBHL(M, v, cs):
    """
    BHL accretion radius [cm]
    """
    return 2 * const.G_cgs * M / (v**2 + cs**2)

def getMu(m1, m2):
    """
    Reduced mass [g]
    """
    return m1 * m2 / (m1 + m2)

def getEdotHL(M, rho, v):
    """
    Accretion energy rate (kinetic energy) [erg/s]
    """
    return getMdotHL(M, rho, v) * v**2

def getDedA(M, m, rho, a):
    """
    Approximate orbital energy gradient [erg/cm]
    Assumes background density `rho` and secondary mass `m`.
    """
    return -0.5 * const.G_cgs * M / a * (4 * np.pi * rho * a**2 - m / a)

def getAccLedd(M):
    """
    Eddington luminosity [erg/s] for mass M [g]
    """
    return 4 * np.pi * const.G_cgs * M * const.m_p_cgs * const.c/const.sigmaT

def getMdotEdd(M, eff=0.1):
    """
    Eddington accretion rate [g/s] assuming radiative efficiency `eff`
    """
    Lum = getAccLedd(M)
    return Lum / (eff * const.c_cgs**2)

def getNSLum(Mdot, eff=0.1):
    """
    Luminosity from accretion onto a neutron star [erg/s]
    """
    return eff * Mdot * const.c_cgs**2

def getStellarL(R, T):
    """
    Stellar luminosity using Stefan–Boltzmann law [erg/s]
    """
    return 4 * np.pi * R**2 * const.sigma_SB_cgs * T**4

def getCFL(cs, dt, dx):
    """
    Return Courant–Friedrichs–Lewy number (dimensionless)
    """
    return cs * dt / dx

def getRsch(M):
    """
    Schwarzchild Radius [cm]
    """
    return 2.0*const.G_cgs*M/c**2

# Adding some prior expressions, need to see papers where they were expressed
#def getMdotEdd(r,k):
#    [Msun yr**-1]
#    return 2e-8*(r/12)*(k/0.34)**(-1)

#def getMdotHyper(k):
#    return 1.9e-4*(k/0.34)**(-0.73)

#def getTrapingRadius(mdot, k):
#    return 5.8e13 * (mdot)*(k/0.34)

#def getShockRadius(mdot):
#    return 1.6e8*(mdot)**(-0,37)

def getReduceMass(m1, m2):
    return m1*m2/(m1+m2)

def getKeplerianPeriod(M, a):
    return 2*np.pi*np.sqrt(a**3/ (const.G*M))

def HamiltonianGrad(M,q,r):
    return const.G*M*q/r**3

def getVelocity(M,a):
    return np.sqrt(const.G*M/a)


# Find what these were
# What are the two following equations?
#r0 = (1.-e0)/(1.+q0)*a0
#v0 = 1./(1.+q0)*np.sqrt((1.+e0)/(1.-e0))*np.sqrt(G*M0/a0)
# Work on Bondi HOyle Tests... 
#rbhl = lambda M, v, c: G*M/(v**2 + c **2) # wind, sound speed
#vr = lambda M, v, r, l: -np.sqrt(v**2 + (2*G*M/r) - l**2 * (v**2/r**2) )
#vtheta = lambda l,v,r : l*v/r
#r = lambda M,l,v,theta: l**2*v**2 /(G*M*(1+np.cos(theta)) + l*v**2*np.sin(theta))
#rho = lambda rhor, l, r, theta : rhor*l**2/ (r*np.sin(theta)*(2*l-r*np.sin(theta)))
#Mdotbhl = lambda M, rhor, c, v: 2*np.pi*G**2*M**2*rhor/ (c**2+v**2)**(3.0/2.0)

