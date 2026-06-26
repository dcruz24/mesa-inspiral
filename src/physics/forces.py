from src.energy import distribution as eh
from src.physics import drag as cedm
from src.physics import jets as cejm


def update_force_rates(runParameters, simParameters, donorStar, accretorStar):
    """Update local envelope values and drag/jet rates for the current orbit state.

    The circular path preserves the original velocity prescription. For the
    two-body path, the orbit integrator owns the relative velocity and this
    routine only refreshes local stellar quantities from the current radius.
    """
    orbit_method = getattr(runParameters, 'orbitIntegrator', 'circular')
    use_circular_velocity = orbit_method == 'circular'

    simParameters.getLocalVariable(donorStar, circular_velocity=use_circular_velocity)
    simParameters.calculateDragRadius(runParameters)
    simParameters.calculateJetRadius(runParameters, donorStar)
    accretorStar.getAccretorValues(simParameters)

    simParameters.epGrad = eh.getEPGrad(
        donorStar.rho,
        donorStar.R,
        simParameters.a,
        simParameters.Ra,
        simParameters.rho,
        simParameters,
    )

    dMdot, Fd, dEdrag = cedm.chooseDragMethod(runParameters, simParameters)
    dEjet = cejm.chooseJetMethod(runParameters, simParameters, accretorStar)

    simParameters.Mdot = dMdot
    simParameters.Fd = Fd
    simParameters.dEdrag = dEdrag
    simParameters.dEjet = dEjet
    simParameters.Edrag = dEdrag * runParameters.dtType

    accretorStar.checkBind(dMdot)
    accretorStar.getAngMom(simParameters)

    simParameters.da_dt = simParameters.a * (dEdrag / abs(simParameters.Eorb))
    simParameters.th_dt = 2.0 * 3.141592653589793 / simParameters.P
    simParameters.dth_dt = simParameters.th_dt
    return
