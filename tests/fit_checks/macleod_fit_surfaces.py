from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

OUT_DIR = Path(__file__).resolve().parent / 'plots'
OUT_DIR.mkdir(parents=True, exist_ok=True)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.physics import drag as cedrag  # Ensure this module contains Cad43 and Cad53

from matplotlib import cm

# Meshgrid for Mach and q
q_vals = np.linspace(0.01, 1.0, 100)
Mach_vals = np.linspace(1.0, 5.0, 100)
Q, M = np.meshgrid(q_vals, Mach_vals)

# Evaluate your original Cad43/Cad53 (non-log values)
Ca43, Cd43 = cedrag.Cad43(Q, M)
Ca53, Cd53 = cedrag.Cad53(Q, M)

# Take log10 to match the paper's A1–A4 presentation
logCa43 = np.log10(Ca43)
logCd43 = np.log10(Cd43)
logCa53 = np.log10(Ca53)
logCd53 = np.log10(Cd53)

def plot_surface(X, Y, Z, title, zlabel, filename, zmin, zmax):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis, edgecolor='none', vmin=zmin, vmax=zmax)
    ax.set_xlabel(r'Mach Number $\mathcal{M}$', fontsize=12)
    ax.set_ylabel(r'Mass Ratio $q$', fontsize=12)
    ax.set_zlabel(zlabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_zlim(zmin, zmax)
    fig.colorbar(surf, shrink=0.5, aspect=10, label=zlabel)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# Generate and save plots with fixed color limits
plot_surface(M, Q, logCa43, r'$\log_{10}(C_a)$ for $\gamma=4/3$', r'$\log_{10}(C_a)$',
             OUT_DIR / 'A1_logCa43.png', zmin=-2, zmax=0.5)
plot_surface(M, Q, logCd43, r'$\log_{10}(C_d)$ for $\gamma=4/3$', r'$\log_{10}(C_d)$',
             OUT_DIR / 'A2_logCd43.png', zmin=-0.2, zmax=1.2)
plot_surface(M, Q, logCa53, r'$\log_{10}(C_a)$ for $\gamma=5/3$', r'$\log_{10}(C_a)$',
             OUT_DIR / 'A3_logCa53.png', zmin=-2, zmax=0.0)
plot_surface(M, Q, logCd53, r'$\log_{10}(C_d)$ for $\gamma=5/3$', r'$\log_{10}(C_d)$',
             OUT_DIR / 'A4_logCd53.png', zmin=0.0, zmax=1.0)