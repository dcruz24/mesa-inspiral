import numpy as np
import matplotlib.pyplot as plt
import plothelpers
import sys
sys.path.append('../')

import CEdragModels as cedrag
import CEphysics    as cephy
import CEconstants   as const


m1  = 1.4 * const.Msun
rho = 1e-7
v   = 1e7
cs  = 1e6
Ra  = 0.5
epArr = np.linspace(0.3,3.0, 100)
mdotFit, fdFit, lfit = cedrag.getNS2015Drag(m1, rho, v, Ra, epArr, True)
mdotFit2, fdFit2, lfit2= cedrag.getNS2015Drag(5*m1, rho, v, Ra, epArr, True)


fig, ax = plt.subplots(facecolor='w')
ax.set(yscale = 'log',
       ylabel = 'Accretion Drag',
       xlabel = r'$\epsilon_p$',
       title = 'APJ 803: 41, Fits from CE Wind Tunnels',
#       ylim = (1e-3, 10e0)
      )
ax.plot(epArr, fdFit, label = r'$F / F_{d, HL}, R_s/R_a = 0.01$')
ax.plot(epArr, mdotFit, label = r'$ \dot M / \dot M_{HL}, R_s/R_a = 0.01$')

ax.plot(epArr, fdFit2, label = r'$F / F_{d, HL}, R_s/R_a = 0.05$')
ax.plot(epArr, mdotFit2, label = r'$ \dot M / \dot M_{HL}, R_s/R_a = 0.05$')

ax.axvline(0.3,color = 'grey')
ax.axvline(3.0,color = 'grey')
ax.legend(loc='best')
ax.grid(True)
plt.tight_layout()
plt.savefig('WindTunnelFit.png', dpi = 300)
#plt.show()


qarr = [0.1, 0.143, 0.2, 0.25, 0.333]  # mass ratio array
Mach = np.linspace(1, 5.0, 100)  # Mach number array

fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
colors = plt.cm.magma(np.linspace(0.2, 0.8, len(qarr)))  # Nice color gradient

for i, q in enumerate(qarr):
    # Get coefficients from fit functions
    Ca43, Cd43 = cedrag.Cad43(q, Mach)
    Ca53, Cd53 = cedrag.Cad53(q, Mach)

    # Plotting
    axes[0, 0].plot(Mach, Ca43, color=colors[i], label=fr'$q_r = {q}$')
    axes[0, 1].plot(Mach, Cd43, color=colors[i], label=fr'$q_r = {q}$')
    axes[1, 0].plot(Mach, Ca53, color=colors[i])
    axes[1, 1].plot(Mach, Cd53, color=colors[i])

# Axis labels
axes[1, 0].set_xlabel(r'$\mathcal{M}_\infty$', fontsize=14)
axes[1, 1].set_xlabel(r'$\mathcal{M}_\infty$', fontsize=14)
axes[0, 0].set_ylabel(r'$C_a$', fontsize=14)
axes[1, 0].set_ylabel(r'$C_a$', fontsize=14)
axes[0, 1].set_ylabel(r'$C_d$', fontsize=14)
axes[1, 1].set_ylabel(r'$C_d$', fontsize=14)

# Titles
axes[0, 0].set_title(r'$\gamma = 4/3$', fontsize=14)
axes[0, 1].set_title(r'$\gamma = 4/3$', fontsize=14)
axes[1, 0].set_title(r'$\gamma = 5/3$', fontsize=14)
axes[1, 1].set_title(r'$\gamma = 5/3$', fontsize=14)

# Legends
axes[0, 0].legend(loc='best', fontsize=10)
axes[0, 1].legend(loc='best', fontsize=10)

# Grid and layout
for ax in axes.flatten():
    ax.grid(True)
    ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig('WindTunnel2.png', dpi=300)
plt.show()


