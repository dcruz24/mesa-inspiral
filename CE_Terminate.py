import CEconstants as const

def checkTermination(runParameters, simParameters):
    
    if(simParameters.a  < runParameters.Rends):
        simParameters.runFlag = False
        print("Reach Final sepeartion")
        
    elif (simParameters.count > runParameters.nmax):
        simParameters.runFlag = False
        print("Max Iteartions")
        
    elif(simParameters.m1 > 1.9*const.Msun_cgs):
        simParameters.runFlag = False
        print('NS Mass Cut OFF M >1.9 Msun', simParameters.m1/const.Msun_cgs)
        
    elif (simParameters.mesaFlag == False):
        simParameters.runFlag = simParameters.mesaFlag
        print('Mesa Stopped Running')
    return

def changeDt(runParameters, sysParm, donorStar):
    m1 = sysParm.m1
    a  = sysParm.a
    m2 = np.interp(a, donorStar.R, donorStar.Menc)
    P = np.sqrt(4.0*np.pi**2 / const.G_cgs*(m1+m2)*a**3.0)
    dt = max([0.005*P, runParameters.dtType])
    runParameters.dtType = dt
    return
