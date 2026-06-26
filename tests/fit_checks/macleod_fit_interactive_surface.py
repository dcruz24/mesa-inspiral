from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

OUT_DIR = Path(__file__).resolve().parent / 'plots'
OUT_DIR.mkdir(parents=True, exist_ok=True)

import numpy as np
import matplotlib.pyplot as plt

from src.physics import drag as cedrag
from src.physics import basic as cephy
from src.config import constants as const

from mpl_toolkits.mplot3d import Axes3D
try:
    import plotly.graph_objects as go
except ImportError:
    go = None

if go is None:
    print('Plotly is not installed; run macleod_fit_surfaces.py for static PNG surfaces.')
    raise SystemExit(0)

# Grid for q and Mach
q_vals = np.linspace(0.01, 0.5, 100)
Mach_vals = np.linspace(1.0, 5.0, 100)
Q, M = np.meshgrid(q_vals, Mach_vals)

# Evaluate log values of Ca and Cd from your fit functions
Ca43, Cd43 = cedrag.Cad43(Q, M)
logCa43 = np.log10(Ca43)
logCd43 = np.log10(Cd43)

def plot_interactive_surface(
    X, Y, Z,
    title,
    zlabel,
    cmin, cmax,
    zmin, zmax,
    xrange=None,
    yrange=None,
    filename=None
):
    fig = go.Figure(data=[
        go.Surface(
            z=Z,
            x=X,
            y=Y,
            colorscale='Earth',
            cmin=cmin,
            cmax=cmax,
            colorbar=dict(title=zlabel)
        )
    ])

    scene_dict = {
        'xaxis_title': 'Mach Number (𝓜)',
        'yaxis_title': 'Mass Ratio (q)',
        'zaxis_title': zlabel,
        'zaxis': dict(range=[zmin, zmax])
    }

    if xrange:
        scene_dict['xaxis'] = dict(range=list(xrange))
    if yrange:
        scene_dict['yaxis'] = dict(range=list(yrange))

    fig.update_layout(
        scene=scene_dict,
        scene_camera=dict(eye=dict(x=2.0, y=-2.0, z=1.0)),
        margin=dict(l=0, r=0, b=0, t=40),
        title=title
    )

    if filename:
        fig.write_html(filename)
        print(f"Plot saved to {filename}")

# Call the plotting function and save as PNG
plot_interactive_surface(
    X=M, Y=Q, Z=logCa43,
    title='log₁₀(Cₐ) for γ=4/3 (Direct Polynomial Fit)',
    zlabel='log₁₀(Cₐ)',
    cmin=-2, cmax=0.5,
    zmin=-2, zmax=1,
    xrange=(1.0, 5.0),
    yrange=(0.0, 0.5),
    filename=OUT_DIR / 'A1_logCa43_interactive.html'
)
