from src.io import inlist as ih

# Input (string, float, float): inlistFile, m, a
# Output (None)               : Writes out inlist file
# Adds compact gravitationl potential... calling other_momentum hook...
def includeCompactPotential(inlistFile, flag, m, a):
    if (flag == True):
        lines   = ih.getTxtArrays(inlistFile)
        lines   = ih.replaceLines(lines, 'use_other_momentum', '.true.')
        lines   = ih.replaceLines(lines, 'x_ctrl(1)', str(m))  # [grams]
        lines   = ih.replaceLines(lines, 'x_ctrl(2)', str(a))  # [cm]
        ih.writeFile(lines, inlistFile)
    return
    
# Input (string) : inlistFile
# Output(None)   : Writes out inlist file
# Turns on wind mass loss... calling other_wind hook...
# Removes outer layer mass based on escape velocity
def includeMassLoss(inlistFile, flag):
    if (flag == True):
        lines   = ih.getTxtArrays(inlistFile)
        lines   = ih.replaceLines(lines, 'use_other_wind', '.true.')  # [cm]
        ih.writeFile(lines, inlistFile)
    return
    
def addRotation(inlistFile, flag, value):
    if (flag == True):
        lines   = ih.getTxtArrays(inlistFile)
        lines   = ih.replaceLines(lines, 'set_omega', '.true.')  # [cm]
        lines   = ih.replaceLines(lines, 'new_omega', str(value))  # [grams]
        lines   = ih.replaceLines(lines, 'new_rotation_flag', '.true.')  # [cm]
        lines   = ih.replaceLines(lines, 'change_rotation_flag', '.true.')  # [cm]
        lines   = ih.replaceLines(lines, 'relax_omega', '.true.')
        ih.writeFile(lines, inlistFile)
    return
    
def turnOffOmega(inlistFile):
    lines   = ih.getTxtArrays(inlistFile)
    lines   = ih.replaceLines(lines, 'set_omega', '.false.')
    lines   = ih.replaceLines(lines, 'relax_omega', '.false.')
    ih.writeFile(lines, inlistFile)
    return

def _drag_kernel_id(kernel):
    normalized = str(kernel).strip().lower()
    if normalized in {"uniform", "flat"}:
        return 0
    if normalized in {"gauss", "gaussian"}:
        return 1
    raise ValueError("drag heating hook supports only Uniform or Gauss/Gaussian kernels")


def setDragEnergyInjection(
    inlistFile,
    flag,
    energy_erg,
    r_inner_cm,
    r_outer_cm,
    duration_sec,
    kernel="Uniform",
    gaussian_sigma_fraction=0.2,
):
    """Configure run_star_extras.other_energy radial drag heating controls.

    x_ctrl mapping used by MesaWorkDir/src/run_star_extras.f90:
      x_ctrl(3): total drag energy for this coupling step [erg]
      x_ctrl(4): inner radius of heated region [cm]
      x_ctrl(5): outer radius of heated region [cm]
      x_ctrl(6): duration over which x_ctrl(3) is deposited [s]
      x_ctrl(7): Gaussian sigma as a fraction of radial interval width
      x_integer_ctrl(1): 0 uniform, 1 Gaussian
    """
    lines = ih.getTxtArrays(inlistFile)
    if flag:
        r_low = min(float(r_inner_cm), float(r_outer_cm))
        r_high = max(float(r_inner_cm), float(r_outer_cm))
        lines = ih.replaceLines(lines, 'use_other_energy', '.true.')
        lines = ih.replaceLines(lines, 'x_ctrl(3)', str(abs(float(energy_erg))))
        lines = ih.replaceLines(lines, 'x_ctrl(4)', str(r_low))
        lines = ih.replaceLines(lines, 'x_ctrl(5)', str(r_high))
        lines = ih.replaceLines(lines, 'x_ctrl(6)', str(float(duration_sec)))
        lines = ih.replaceLines(lines, 'x_ctrl(7)', str(float(gaussian_sigma_fraction)))
        lines = ih.replaceLines(lines, 'x_integer_ctrl(1)', str(_drag_kernel_id(kernel)))
    else:
        lines = ih.replaceLines(lines, 'use_other_energy', '.false.')
        lines = ih.replaceLines(lines, 'x_ctrl(3)', '0d0')
    ih.writeFile(lines, inlistFile)
    return

