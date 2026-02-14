import numpy as np
from bisect import bisect_right
from scipy.stats import triang, norm
import CEconstants as const
import MESAGen as mg
import InlistHelper as ih

# ----------------------------------------------------------------------------------------------
def getPointIndex(R, x):
    R = R[::-1]
    realIndex = bisect_right(R, x)
    modIndex = len(R)-realIndex
    return modIndex

# ----------------------------------------------------------------------------------------------
def getAccretionZones(R, a, xhigh, xlow, limit=1):
    highIndex = getPointIndex(R, xhigh)
    lowIndex  = getPointIndex(R, xlow)
    if (highIndex < limit):
         return np.arange(limit, lowIndex+1)
    return np.arange(highIndex, lowIndex+1)

# ----------------------------------------------------------------------------------------------------
def getEnergyDistribution(R, zones, dq, starM, Einput, kernel):
    firstZone = zones[0]
    lastZone  = zones[-1]
    accRadius = R[firstZone-1:lastZone]
    dq        = dq[firstZone - 1: lastZone] * starM
    totalM    = sum(dq)
    fraction  = np.linspace(0, 1, len(dq))  # Create an array of fractions from 0 to 1

    if (kernel == 'Triangle'):
        # Triangular distribution
        c = 0.5  # Center of the triangular distribution (adjust as needed)
        weights = triang.pdf(fraction, c)
    elif(kernel == 'Parabolic'):
        # Parabolic distribution
        weights = 1 - (2 * fraction - 1)**2
    elif(kernel == 'Cosine'):
        mu = 0.5  # Center of the cosine distribution (adjust as needed)
        weights = 0.5 * (1 + np.cos((fraction - mu) * np.pi))
    elif(kernel == 'Quartic'):
        mu = 0.5  # Center of the quartic distribution (adjust as needed)
        weights = (1 - (fraction - mu)**2)**2
    elif(kernel == 'Triweight'):
        weights = (35 / 32) * (1 - (fraction**2))**3
    elif(kernel == 'Gauss'):
        mu = 0.5  # Mean of the Gaussian distribution (adjust as needed)
        sigma = 0.1  # Standard deviation of the Gaussian distribution (adjust as needed)
        weights = norm.pdf(fraction, loc=mu, scale=sigma)
    elif(kernel == 'Uniform'):
        weights = np.ones_like(fraction)
    elif(kernel == 'Os'):
        sigma = 0.1
        gaussian_center = np.exp(-(fraction - 0.5)**2 / (2 * sigma**2))
        weights = gaussian_center * (0.5 + 0.5 * np.sin(2 * np.pi * 50.0 * fraction))
        
    Ezones = Einput * (weights / sum(weights))
    return Ezones, accRadius
    
# ----------------------------------------------------------------------------------------------------
def getEPGrad(rhos, rs, a, Ra, rhoLocal, simParameters):
    
    RaccUpper = simParameters.RupperDrag
    RaccLower = simParameters.RlowerDrag
    
    Upperind = bisect_right(rs,RaccUpper)
    Lowerind = bisect_right(rs,RaccLower)
    
    deltaR   = rs[Upperind-1]- rs[Lowerind]
    deltaRho = rhos[Upperind-1]-rhos[Lowerind]
    Hp       = -1.0*rhoLocal*(deltaR/deltaRho)
    ep       = np.abs(Ra/Hp)
    return ep
