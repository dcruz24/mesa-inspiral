import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
import pandas as pd
import time
import utility as ut
import os

G    = 6.6743e-8
c    = 2.99792458e10
Msun = 1.9891e33
Rsun = 6.957e10
pc   = 3.08567758128e18
kpc, Mpc = 1.e3*pc, 1.e6*pc
au   = 1.49598073e13
yr   = 31556952.

# Note: Iterates per mesa profile based on user defined interval, in case there are too many mesa profiles outputs..(default=1)
# Choose interval number of profiles to read. Dont usually need to change
class KippenHahn1D:
    def __init__(self, starPath, donorPath, profileNum = 1):
        self.inspiralData = ut.readInspiralData(starPath, donorPath) # Hold all the data from the inspiral run
        self.donorPath = donorPath
        self.newts = []
        self.ts    = self.inspiralData.tl
        print('Number of Total Profiles: ', len(self.ts))
        self.profileNum = profileNum
        for a in range(0, len(self.ts)):
            if(a % self.profileNum == 0):
                self.newts.append(self.ts[a])
        print('New number of profiles to iterate: ', len(self.newts))
        self.ts = self.newts
        print('Initialized....')
        
    def setGrid(self, rgridL, rgridR, sgrid = 2000):
        self.ns   = np.arange(0 ,len(self.ts), step=1)
        self.Rz   = np.logspace(np.log10(rgridL*Rsun),
                                np.log10(rgridR*Rsun),
                                sgrid)
        self.varg = np.zeros((sgrid, len(self.ns)))
        self.tgtest , self.rgtest = np.meshgrid(self.ts, self.Rz)
        print('Grid Set Up')

    # name: Profile variable of interest
    def varInterpolate(self, name, grad, vesc, basefile):
        t1 = time.time()
        temp = 0
        for n in range(0, len(self.ts)):
            #if(n % 1 == 0):
            file = self.donorPath + '{}_'.format(n) + basefile
            df   = pd.read_csv(file,encoding="latin-1", 
                                   delim_whitespace=True,skip_blank_lines=True,skiprows=5)
                
            Rs   = 10.**df['logR'].to_numpy()[::-1]*Rsun
                
            if(name[0:3] == 'log'):
                var = 10**df[name].to_numpy()[::-1]
            else:
                var = df[name].to_numpy()[::-1]

            # Calculates dr/drho
            if(grad == True and name  == 'logRho'): 
                var =  -1.*Rs/var*np.gradient(var,Rs)
            
            # Calculates escape velocity
            if(vesc == True and name == 'velocity'):
                Mgrav   = df['m_grav'].to_numpy()[::-1]*Msun*1e-3
                escapeV = np.sqrt((2*G*Mgrav)/Rs)
                var     = var/escapeV
                
            vari = np.interp(self.Rz,Rs,var)
            inds = np.argwhere(self.Rz > Rs[-1])
            vari[inds] = 0.
            self.varg[:,temp] = vari
            temp += 1
        print('Total interpolation time: {}'.format(time.time()-t1))        
        return self.varg
    
# Generates kippenhahn data for list of interested variables
def SaveKippenhanData(self, variables, filename, path, run, grad = False, vesc = False):
    storDir = path + filename + '_'
    np.savetxt(storDir + 'tgtest_' + str(run), self.tgtest); print('Saved Time 2D Arr.')
    np.savetxt(storDir + 'rgtest_' + str(run), self.rgtest); print('Saved Rado 2D Arr.')
    
    for a in range(0, len(variables)):
        var     = self.varInterpolate(variables[a], False, False, filename)
        np.savetxt(storDir + variables[a] + '_' + str(run), var)        ; print('Saved '+ variables[a]+' 2D Arr.')
        if(variables[a] == 'logRho' and grad == True ):
            var = self.varInterpolate(variables[a], grad, vesc, filename)
            np.savetxt(storDir + 'drho_dr' + '_' + str(run), var)       ; print('Saved drho/dr 2D Arr.')            
        if(variables[a] == 'velocity' and vesc == True ):
            var = self.varInterpolate(variables[a], grad, vesc, filename)
            np.savetxt(storDir + 'vesc' + '_' + str(run), var)          ; print('Saved Escape Velocity 2D Arr.')
    
class readKippenHahnData:
    def __init__(self, directoryPath, run):
        fileList = os.listdir(directoryPath)
        parsedList = [item.split('_')[1] for item in fileList]
        print('Found Data for the following...')
        print('\n'.join(parsedList))
        for a in range(0, len(parsedList)):
            fileName = fileList[a]
            setattr(self, parsedList[a], np.loadtxt(directoryPath + fileName))


# ---------------------------------------- NEEDS DEVELOPMENT --------------------------------
class KippenHahnTwoBody:
    def __init__(self, pathns,pathrsg):
        self.ts, self.m1s, self._ms = np.loadtxt(pathns + '1',dtype = float, unpack=True)
        self.xa, self.ya, self.za, self.xb, self.yb, self.zb, self.vxa, self.vya,self.vza, self.vxb, self.vyb,self.vzb = np.loadtxt(pathns +'2', dtype = float, unpack=True)
        x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2 = self.xa, self.ya, self.za, self.xb, self.yb, self.zb, self.vxa, self.vya, self.vza, self.vxb, self.vyb, self.vzb
        xs,ys,zs = x1-x2,y1-y2,z1-z2
        vxs,vys,vzs = vx1-vx2,vy1-vy2,vz1-vz2
        self.rs,self.vs = np.sqrt(xs*xs+ys*ys+zs*zs),np.sqrt(vxs*vxs+vys*vys+vzs*vzs)        
        self.pathStar1 = pathns
        self.pathStar2 = pathrsg
        self.v_ = np.sqrt(G*(self.m1s+self._ms)/self.rs)
        self.rc = 2.*G*self.m1s/(self.v_*self.v_)
        self.r_ac = [self.rs + .5*self.rc, self.rs - .5*self.rc]
                
    def set_grid(self,rgridL,rgridR,sgrid):
        na, nb = 0, len(self.ts)
        self.ns = np.arange(na,nb,step=1)
        nx = int(sgrid)
        Ra,Rb = rgridL*Rsun, rgridR*Rsun
        self.Rz = np.logspace(np.log10(Ra),np.log10(Rb),nx)
        nt = len(self.ns) ## Size of our time array
        self.varg = np.zeros((nx,nt))
        self.tgtest , self.rgtest = np.meshgrid((self.ts[0:len(self.ns-1)]), self.Rz)
        print('Grid Set Up')

    def var_interpolate(self,name,grad,vesc,basefile):
        t1 = time.time()
        for n in self.ns:
            file = self.pathStar2 + '{}_'.format(n) + basefile
            df = pd.read_csv(file,encoding="latin-1",delim_whitespace=True,skip_blank_lines=True,skiprows=5)
            Rs = 10.**df['logR'].to_numpy()[::-1]*Rsun
            
            if(name[0:3] == 'log'):var = 10**df[name].to_numpy()[::-1]
            else:var = df[name].to_numpy()[::-1]
            
            if(grad == True and name  == 'logRho'): var =  -1.*Rs/var*np.gradient(var,Rs)            
            if(vesc == True and name == 'velocity'):
                Mgrav = df['m_grav'].to_numpy()[::-1]*Msun
                vesc = np.sqrt((2*G*Mgrav)/Rs)
                var = var/vesc
                
            vari = np.interp(self.Rz,Rs,var)
            inds = np.argwhere(self.Rz>Rs[-1])
            vari[inds] = 0.
            self.varg[:,n] = vari
            
        print('Finished Interpolate')
        print('Total Time: ', time.time() - t1)
        return self.varg