import os
import shutil
import sys
import time
from src.mesa import interface as ms
from src.io import inlist as ih
from src.mesa import inlist_controls as ms_inlist


def _import_py_run_star_support(runParameters):
    candidate_paths = []
    configured_path = getattr(runParameters, 'pyMesaMainSrcPath', None)
    if configured_path:
        candidate_paths.append(os.path.abspath(configured_path))
    candidate_paths.append(os.path.abspath(os.path.join(runParameters.inspiralPath, '..', 'pyMesaMain', 'src')))

    for path in candidate_paths:
        if os.path.isdir(path):
            if path not in sys.path:
                sys.path.insert(0, path)
            try:
                from pyRunLib import pyRunStarSupport
                return pyRunStarSupport
            except ImportError:
                continue

    raise ImportError(
        'Could not import pyRunStarSupport from pyMesaMain. Set runParameters.pyMesaMainSrcPath '
        'to the pyMesaMain/src directory or set mesaMode=0.'
    )


def _run_python_interval_backend(runParameters):
    pyRunStarSupport = _import_py_run_star_support(runParameters)
    os.environ['PYMESA_LOAD_MODEL_FILENAME'] = os.path.abspath(runParameters.modFile)
    runner = pyRunStarSupport(
        cacheDir=os.path.join(runParameters.inspiralPath, 'pyMesaStarCache'),
        workDir=runParameters.pathMesa,
        omp_threads=runParameters.numThreads,
    )
    return runner.run1_star(
        inlist='inlist',
        allocate=True,
        free=True,
        allow_restart=False,
        restart=False,
    )


def shutdown_mesa_runtime():
    return


def updateRunMesa(pyMesaObject, simParameters, runParameters, donorStar):
    """Advance MESA with drag feedback applied by the Fortran other_energy hook."""

    RupperDrag = simParameters.RupperDrag
    RlowerDrag = simParameters.RlowerDrag
    Edrag = abs(simParameters.Edrag)
    dragKernel = runParameters.dragKernel
    dt = runParameters.dtType

    starAge = float(donorStar.starAge)

    drag_enabled = bool(runParameters.includeDrag) and Edrag > 0.0
    ms_inlist.setDragEnergyInjection(
        runParameters.inlistFile,
        drag_enabled,
        Edrag,
        RlowerDrag,
        RupperDrag,
        dt,
        kernel=dragKernel,
    )

    if not runParameters.mesaRun:
        if simParameters.count == 0:
            print('MESA run skipped because runParameters.mesaRun is False.', flush=True)
        return None

    ih.addAgeInlist(runParameters.inlistFile, dt)

    mesa_mode = int(getattr(runParameters, 'mesaMode', 0))

    mesaStartTime = time.time()
    if mesa_mode == 0:
        ms.RunMesa(runParameters.pathMesa, runParameters.numThreads)
        backend_name = './rn'
    elif mesa_mode == 1:
        _run_python_interval_backend(runParameters)
        backend_name = 'py-run1-star'
    else:
        raise ValueError(f'Unsupported runParameters.mesaMode={mesa_mode}. Expected 0 or 1.')

    mesaFinalTime = time.time() - mesaStartTime
    print("Mesa Called...            {} [secs] ({})".format(mesaFinalTime, backend_name), flush=True)

    if not os.path.exists(runParameters.finalMod):
        raise FileNotFoundError(
            f"MESA completed but did not write expected model file: {runParameters.finalMod}"
        )
    if not os.path.exists(runParameters.profile):
        raise FileNotFoundError(
            f"MESA completed but did not write expected profile file: {runParameters.profile}"
        )

    # Preserve the MESA output model as the input model for the next coupling interval
    # without parsing or rewriting the .mod file in Python.
    if os.path.abspath(runParameters.finalMod) != os.path.abspath(runParameters.modFile):
        shutil.copy2(runParameters.finalMod, runParameters.modFile)

    donorStar.readMesaProfile(runParameters.profile)
    nextStarAge = float(donorStar.starAge)

    if nextStarAge == starAge:
        print('Old Age: ', starAge)
        print('New Age: ', nextStarAge)
        print("MESA: Star did not evolve, check MesaOutPut.txt")
        simParameters.mesaFlag = False
        return None

    return None
