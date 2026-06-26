import numpy as np
from src.config import constants as const


def initialize_orbit_state(simParameters):
    """Create relative two-body state from the scalar circular initial state."""
    if not hasattr(simParameters, 'x'):
        simParameters.x = simParameters.a
        simParameters.y = 0.0
        simParameters.vx = 0.0
        simParameters.vy = simParameters.v
    return


def step_orbit(runParameters, simParameters, accretorStar, dt):
    orbit_method = getattr(runParameters, 'orbitIntegrator', 'circular')
    if orbit_method == 'circular':
        step_circular_euler(simParameters, accretorStar, dt)
    elif orbit_method == 'two_body':
        step_two_body_euler(simParameters, accretorStar, dt)
    else:
        raise ValueError(f"Unknown orbitIntegrator '{orbit_method}'")
    return


def step_circular_euler(simParameters, accretorStar, dt):
    """Original circular-orbit explicit Euler update."""
    simParameters.m1 = simParameters.m1 + dt * accretorStar.dm_dt
    simParameters.a = simParameters.a + dt * simParameters.da_dt
    simParameters.th = simParameters.th + dt * simParameters.th_dt
    simParameters.Om = simParameters.Om + dt * accretorStar.dOm_dt
    return


def step_two_body_euler(simParameters, accretorStar, dt):
    """Explicit Euler step for the relative two-body orbit.

    The state is the accretor-donor relative vector. Gravity uses the enclosed
    donor mass at the current separation plus the accretor mass. Drag acts
    opposite the relative velocity with magnitude supplied by the selected drag
    model.
    """
    initialize_orbit_state(simParameters)

    x = simParameters.x
    y = simParameters.y
    vx = simParameters.vx
    vy = simParameters.vy

    r = np.hypot(x, y)
    speed = np.hypot(vx, vy)
    if r <= 0.0:
        raise ValueError('Two-body orbit radius became non-positive')

    total_mass = simParameters.m1 + simParameters.m2
    ax = -const.G_cgs * total_mass * x / r**3
    ay = -const.G_cgs * total_mass * y / r**3

    if speed > 0.0:
        mu = simParameters.m1 * simParameters.m2 / total_mass
        drag_acc = getattr(simParameters, 'Fd', 0.0) / mu
        ax -= drag_acc * vx / speed
        ay -= drag_acc * vy / speed

    simParameters.x = x + dt * vx
    simParameters.y = y + dt * vy
    simParameters.vx = vx + dt * ax
    simParameters.vy = vy + dt * ay

    simParameters.m1 = simParameters.m1 + dt * accretorStar.dm_dt
    simParameters.a = np.hypot(simParameters.x, simParameters.y)
    simParameters.v = np.hypot(simParameters.vx, simParameters.vy)
    simParameters.th = np.arctan2(simParameters.y, simParameters.x)
    simParameters.Om = simParameters.Om + dt * accretorStar.dOm_dt
    return
