import numpy as np
import InlistHelper as ih

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
        lines   = ih.getTxtArr ays(inlistFile)
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