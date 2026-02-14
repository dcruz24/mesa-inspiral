import numpy as np
import plothelpers
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import CEconstants as const


eos = '../eos_tables/SLY4.dat'
M,Rn,dM,dR,enu,bI,h0 = np.loadtxt(eos,unpack=True)
M,Rn,dM,dR = M*const.Msun_cgs, Rn/const.R_km_geo*1.e5, dM*const.Msun_cgs, dR/const.R_km_geo*1.e5
I  = bI*(M*const.Msun_cgs)**3.*const.G_cgs*const.G_cgs/const.c_cgs**4.0

fig, ax = plt.subplots(facecolor='w')
ax.set(xlabel = 'Radius', ylabel = 'M')
ax.plot(Rn, M/const.Msun_cgs)
ax.grid(True)

fig, ax = plt.subplots(facecolor='w')
ax.set(xlabel = 'Radius', ylabel = 'Energy Density', )
ax.plot(Rn, enu)
ax.grid(True)

fig, ax = plt.subplots(facecolor='w')
ax.set(xlabel = 'Radius', ylabel = '')
ax.plot(Rn, bI)
ax.grid(True)

fig, ax = plt.subplots(facecolor='w')
ax.set(xlabel = 'Radius', ylabel = '')
ax.plot(Rn, h0)
ax.grid(True)

plt.show()