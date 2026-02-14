import numpy as np
import CEconstants as const


def chooseJetMethod(runParameters, simParameters, accretorStar):
    m1     = simParameters.m1
    Mdot   = simParameters.Mdot
    jetEff = runParameters.jetEff
    dEjet  = getJetfeedback(m1, Mdot, accretorStar.R, jetEff)
    return dEjet
    
# ----------------------------------------------------------------------------------------------------
def getJetfeedback(M, dMdot, r, jetEff):
    """"
    Grichener 2021, APJ 922, 61
    M: Mass of Accretor
    r: Radius of Accretor
    jetEff: Efficiency paramter taken from Schroder et al, 2020
    """
    return (jetEff * const.G_cgs * M * dMdot) / r
