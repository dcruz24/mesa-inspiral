from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


import matplotlib.pyplot as plt
import numpy as np

from src.config import constants as const

OUT_DIR = Path(__file__).resolve().parent / 'plots'
OUT_DIR.mkdir(parents=True, exist_ok=True)

EOS_TABLE = Path(__file__).resolve().parents[2] / 'eos_tables' / 'SLY4.dat'
M, Rn, dM, dR, enu, bI, h0 = np.loadtxt(EOS_TABLE, unpack=True)
M = M * const.Msun_cgs
Rn = Rn / const.R_km_geo * 1.0e5
dM = dM * const.Msun_cgs
dR = dR / const.R_km_geo * 1.0e5
I = bI * M**3.0 * const.G_cgs**2.0 / const.c_cgs**4.0

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
axes = axes.ravel()
axes[0].plot(Rn, M / const.Msun_cgs)
axes[0].set(xlabel='Radius [cm]', ylabel='M [Msun]')
axes[1].plot(Rn, enu)
axes[1].set(xlabel='Radius [cm]', ylabel='Energy density')
axes[2].plot(Rn, I)
axes[2].set(xlabel='Radius [cm]', ylabel='Moment of inertia [cgs]')
axes[3].plot(Rn, h0)
axes[3].set(xlabel='Radius [cm]', ylabel='h0')
for ax in axes:
    ax.grid(True)
fig.tight_layout()
fig.savefig(OUT_DIR / 'SLY4_eos_curves.png', dpi=300)
plt.close(fig)
