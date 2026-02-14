import numpy as np
import matplotlib.pyplot as plt
import plothelpers
import sys
sys.path.append('../')

import CEdragModels as cedrag
import CEphysics    as cephy
import CEconstants  as const

from mpl_toolkits.mplot3d import Axes3D
import plotly.io as pio
pio.renderers.default = 'notebook_connected'
import plotly.graph_objects as go

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

    fig.show()

    if filename:
        fig.write_image(filename, width=1000, height=800)
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
    filename="A1_logCa43.png"
)
