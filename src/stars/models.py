import mesa_reader as mr
from src.config import constants as const
import numpy as np
"""
    Currently setting up to read MESA files using mesa_reader. Would there be
    other ways to determine stellar structures?

    Holds the 1D arrays from the MESA profile and helper function to read and 
    updated the variables we need

"""
class donorStar:
    def __init__(self, filename):
        starData = mr.MesaData(filename)
        self.R     = 10**starData.data('logR')[::-1] * const.Rsun_cgs
        self.rho   = 10**starData.data('logRho')[::-1]
        self.P     = 10**starData.data('logP')[::-1]
        self.dm    = starData.data('dm')[::-1]
        self.TotalE= starData.data('total_energy')[::-1]
        self.Hrho  = starData.data('scale_height')[::-1] * const.Rsun_cgs
        self.Hp    = starData.data('pressure_scale_height')[::-1] * const.Rsun_cgs
        self.gamma = starData.data('gamma1')[::-1]
        self.Ebind = self.TotalE * self.dm
        self.Menc  = starData.data('m_grav')[::-1] * const.Msun_cgs
        #self.Menc  = spi.cumulative_trapezoid(
        #    y = 4. * np.pi * self.R**2 * self.rho,
        #    x = self.R,
        #    initial = 0.0
        #)
        self.Mass = round(starData.header_data['star_mass'], 4)
        self.starAge = float(starData.header_data['star_age'])
        # Speed of sound squared approximation from pressure and density gradients
        self.cs2 = self.P / self.rho * self.Hrho / self.Hp
        self.cs  = np.sqrt(self.cs2)
        return
        
    def readMesaProfile(self, filename):
        starData = mr.MesaData(filename)
        self.R     = 10**starData.data('logR')[::-1] * const.Rsun_cgs
        self.rho   = 10**starData.data('logRho')[::-1]
        self.Menc  = starData.data('m_grav')[::-1] * const.Msun_cgs
        self.Mass  = self.Menc[-1]
        self.starAge = float(starData.header_data['star_age'])
        self.P     = 10**starData.data('logP')[::-1]
        self.dm    = starData.data('dm')[::-1]
        self.Hrho  = starData.data('scale_height')[::-1] * const.Rsun_cgs
        self.Hp    = starData.data('pressure_scale_height')[::-1] * const.Rsun_cgs
        self.TotalE= starData.data('total_energy')[::-1]
        self.gamma = starData.data('gamma1')[::-1]

        self.cs2 = self.P / self.rho * self.Hrho / self.Hp
        self.cs  = np.sqrt(self.cs2)
        self.Ebind = self.TotalE * self.dm
        return
        