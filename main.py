import time

from src.config import runtime as cp
from src.stars import models as sm
from src.accretors import models as accretorModel
from src.evolution import state as ip
from src.energy import injection as energy_injection
from src.mesa import interface as mesa
from src.io import inlist as inlist_io
from src.io import output as output_io
from src.evolution import termination
from src.physics import forces
from src.physics import orbits


def build_initial_state():
    run_parameters = cp.RunTimeParameters()
    py_mesa = mesa.pyMesa_context('pyMesaStarCache/')

    donor_star = sm.donorStar(run_parameters.starModelPath + run_parameters.starProfileModel)
    accretor_star = accretorModel.accretorStar(run_parameters)
    accretor_star.loadEOS()
    sim_parameters = ip.InspiralParameters(run_parameters, donor_star, accretor_star)

    if run_parameters.orbitIntegrator == 'two_body':
        orbits.initialize_orbit_state(sim_parameters)

    inlist_io.initInspiralFiles(run_parameters, donor_star)
    inlist_io.initializeMesaInlist(run_parameters, sim_parameters, donor_star)
    output_io.writeInspiralParamters(run_parameters, sim_parameters, donor_star)

    return run_parameters, sim_parameters, donor_star, accretor_star, py_mesa


def evolve(run_parameters, sim_parameters, donor_star, accretor_star, py_mesa):
    output_io.saveCEData(sim_parameters, run_parameters)

    while sim_parameters.runFlag:
        itertime = time.time()
        donor_star.readMesaProfile(run_parameters.profile)
        dt = run_parameters.dtType
        sim_parameters.t2 = sim_parameters.t1 + dt

        forces.update_force_rates(run_parameters, sim_parameters, donor_star, accretor_star)
        orbits.step_orbit(run_parameters, sim_parameters, accretor_star, dt)

        energy_injection.updateRunMesa(
            py_mesa,
            sim_parameters,
            run_parameters,
            donor_star,
        )

        termination.checkTermination(run_parameters, sim_parameters)

        output_io.printSimulationOutput(run_parameters, sim_parameters, donor_star)
        output_io.saveCEData(sim_parameters, run_parameters)
        output_io.saveInspiralFile(
            run_parameters.profile,
            run_parameters.pathDonor,
            sim_parameters.savenum,
            '_',
        )
        output_io.checkPointSave(run_parameters, sim_parameters)

        sim_parameters.savenum += 1
        sim_parameters.count += 1
        sim_parameters.t1 = sim_parameters.t2
        sim_parameters.wallclock += time.time() - itertime

    output_io.restartPointSave(run_parameters)


def copy_run_inputs(run_parameters):
    inlist_io.copyFile(
        run_parameters.inspiralPath + '/InspiralRuntimeParameters',
        run_parameters.pathAccretor,
        'InspiralRuntimeParameters.txt',
    )
    inlist_io.copyFile(
        run_parameters.inspiralPath + '/inlist_project',
        run_parameters.pathAccretor,
        'inlist_project.txt',
    )


def main():
    run_parameters, sim_parameters, donor_star, accretor_star, py_mesa = build_initial_state()
    try:
        evolve(run_parameters, sim_parameters, donor_star, accretor_star, py_mesa)
        copy_run_inputs(run_parameters)
    finally:
        energy_injection.shutdown_mesa_runtime()


if __name__ == '__main__':
    main()
