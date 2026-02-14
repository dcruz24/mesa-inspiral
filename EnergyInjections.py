import time
import MESAGen as ms
import InlistHelper as ih
import CEconstants as const
import numpy as np
import EnergyHelper as eh

def updateRunMesa(pyMesaObject, simParameters, runParameters, donorStar, modData):
    
    RupperDrag    = simParameters.RupperDrag
    RlowerDrag    = simParameters.RlowerDrag
    a             = simParameters.a
    Edrag         = abs(simParameters.Edrag)
    Ejet          = abs(simParameters.Ejet)
    dragKernel    = runParameters.dragKernel
    jetKernel     = runParameters.jetKernel
    dt            = runParameters.dtType
    totalStarMass = donorStar.Mass

    #print("Drag Energy: ", Edrag)
    #print("Star Mass: ", totalStarMass/const.Msun_cgs)
    
    # 1.Get log arrays from zone and change format to edit
    lnR      = ih.changeDEModfile(modData['lnR'], 'D', 'E')
    dq       = ih.changeDEModfile(modData['dq'], 'D', 'E')
    lnT      = ih.changeDEModfile(modData['lnT'], 'D', 'E')
    lnd      = ih.changeDEModfile(modData['lnd'], 'D', 'E')
    starAge  = float(ih.replaceChar(modData['star_age'], 'D', 'E'))

    #print('Before...')
    #print('Star Age: ', starAge)
    #print('Inlist Max Age: ', ih.readInlistAge(runParameters.inlistFile))
    
    R        = np.exp(lnR)
    dq       = dq
    T        = np.exp(lnT)
    rho      = np.exp(lnd)
    args = [dq, T, rho]
    
    # 2. Get the zones to modify form the modfile based on the ranges
    #    a) This step is importnat to onyl call modfile and create data once since
    #       this takes up a lot time if have fucntions calling modfile to get data.
    zones = eh.getAccretionZones(R, a, RupperDrag, RlowerDrag)
    #print('zones: ', zones[0], zones[-1])
    # Get energy distribution.
    Ezones, AccRegion = eh.getEnergyDistribution(R, zones, dq, totalStarMass,
                                             Edrag, dragKernel)

    # fig, ax = plt.subplots()
    #ax.plot(zones, Ezones)

    # 3. Get the isotopes header defined, need to be done each iteartion since istopes can change
    isotopeData = ms.getIsotopes(modData, zones)

    # 4. Deposit the energy needed modifying the MESA Mod Data
    finalModData = ms.addInternalEnergy(pyMesaObject, modData, zones, Ezones, isotopeData, 
                     totalStarMass, args)

    lnT = np.array([float(val.replace('D', 'E')) for val in finalModData['lnT']])
    T = np.exp(lnT)
    for i in range(1, len(T)):
        rel_jump = abs(T[i] - T[i-1]) / T[i]
        if rel_jump > 0.3:
            print(f"Zone {i-1} → {i}: T = {T[i-1]:.3e} → {T[i]:.3e}, ΔT/T = {rel_jump:.2f}")

    
    # 5. Now Write the Modfied Mod data back into MOD file
    ms.writeModfile(finalModData, runParameters.modFile)  # Writing modfile with energy deposited...

    # 6. Run Mesa.
    ih.addAgeInlist(runParameters.inlistFile, dt)  # Increase inlist max_age...

    #print('after add age...')
    #print('Star Age: ', starAge)
    #print('Inlist Max Age: ', ih.readInlistAge(runParameters.inlistFile))

    
    mesaStartTime = time.time()
    ms.RunMesa(runParameters.pathMesa, runParameters.numThreads)
    mesaFinalTime = time.time() - mesaStartTime
    print("Mesa Called...            {} [secs]".format(mesaFinalTime), flush = True)
    
    # At this point MESA has a new modfile present.
    # Maybe just return final mod data after reading once, and 
    #   set it up outside for the next iteration
    # since int setp 5, we always write it out prep for the next iteration
    # overwriting previous.
    finalModData = ms.readModfile(runParameters.finalMod)
    nextStarAge  = float(ih.replaceChar(finalModData['star_age'], 'D', 'E'))
    #ms.writeModfile(newModFile, runParameters.modFile) #

    if (nextStarAge == starAge):
        print('Old Age: ', starAge)
        print('New Age: ', nextStarAge)
        print("MESA: Star did not evolve, check MesaOutPut.txt")
        simParameters.mesaFlag = False
        return finalModData

    return finalModData