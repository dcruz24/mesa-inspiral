import subprocess
from src.config import constants as const

# Function: writeFile()
# Function: replaceInlistVar()
# Function: replaceLines()
# Function: copyFile()
# Function: replaceChar()
# Function: readInlistAge()
# Function: addAgeInlist()
# Function: changeDEModfile()

def copyFile(filePath, targetPath, outputFile):
    subprocess.run("cp " + filePath + ' ' + targetPath + '/' + outputFile,
              shell = True)
    return
    
def getTxtArrays(filename):
    with open(filename) as f:
        # returns list of strings
        lines = [line.rstrip() for line in f]
    return lines

# Input (List[strings], string): lines, targetName
# Output (Int)                 : Index where the string is located
def getTxtIndex(lines, targetName):
    # searching for index where each inlist initalized is located...
    for a in range(0, len(lines)):
        if(len(lines[a].split()) > 2):
            _name = lines[a].split()[0]
            if(_name == targetName): 
                index = a                
    return index

# Input (list[strings], string) : lines, filename
# Output: None [Writes outfile with filename]
def writeFile(lines, filename):
    newfile = open(filename, "w")    
    for a in range(0,len(lines)):
        newfile.write(lines[a] + "\n")
    newfile.close()

# Input (string, string, string) : inlistfile, variable, target
# Output: (List[strings])        : Modify values in inlist
def replaceInlistVar(inlistfile, variable, target):
    lines = getTxtArrays(inlistfile)
    Ind = getTxtIndex(lines, variable)
    temp = lines[Ind].split()
    temp[-1] = str(target)
    temp = '      ' + ' '.join(temp)
    lines[Ind] = temp
    return lines

# SIMILAR to above, but takes in the lines if you want edit multiple values....
# Input (List[strings], string, string)
# Output: (List[strings])        : Modify values in inlist
def replaceLines(lines, variable, target):
    Ind = getTxtIndex(lines, variable)
    temp = lines[Ind].split()
    temp[-1] = str(target)
    temp = '      ' + ' '.join(temp)
    lines[Ind] = temp
    return lines
    
# Replaces a single character in a string 
# current: The character to change within the string
# target: The desired character to replace with
def replaceChar(string, current, target):
    parsedStr = list(string)
    correctStr = ''
    # Go through each char in string
    for a in range(0, len(parsedStr)):
        if(parsedStr[a] == current):
            parsedStr[a] = target
        correctStr += parsedStr[a]
    return correctStr
# =========================================================================================================
# Reads star age from mesa inlist
def readInlistAge(inlistfile):
    lines = getTxtArrays(inlistfile)
    ind = getTxtIndex(lines, 'max_age')
    return lines[ind].split()[-1]

# Gets starAge from mod file
def getstarAge(modData):
    starAge = float(replaceChar(modData['star_age'], 'D', 'E'))
    return starAge

# Increments Mesa Inlist star age by dt
def addAgeInlist(inlistfile, dt):
    inlistAge = readInlistAge(inlistfile)
    newAge = float(inlistAge) + float(dt)/const.yr_cgs
    newAge = str(newAge)
    lines = replaceInlistVar(inlistfile, 'max_age', newAge)    
    writeFile(lines, inlistfile)
    return

# Returns an mod file array as float types to modify... Needed for internal Mesa uses...
def changeDEModfile(modDataVar, currentChar, targetChar):
    try:
        import importlib

        np = importlib.import_module('numpy')
    except ImportError as exc:
        raise ImportError(
            "numpy is required for changeDEModfile but is not installed"
        ) from exc

    newvalLst = []
    for a in range(0, len(modDataVar)):
        if(targetChar == 'D'):
            val = replaceChar(str(modDataVar[a]), currentChar, targetChar)
            newvalLst.append(val)
        elif(targetChar == 'E' or targetChar == 'e'):
            val = replaceChar(modDataVar[a], currentChar, targetChar)
            newvalLst.append(val)
    
    if(targetChar == 'D'):
        newvalLst = np.array(str(newvalLst))
    elif(targetChar == 'E' or targetChar == 'e'):
        newvalLst = np.array(newvalLst).astype(float)
    return newvalLst


# Resets inlist star age to match the current MESA profile age [used for fresh runs]
def setInitialInlistAge(runParameters, simParameters, starAge):
    starAge = float(starAge)
    simParameters.MESA_currAge = starAge
    simParameters.MESA_maxAge  = starAge
    currAge = str(starAge)
    lines   = replaceInlistVar(runParameters.inlistFile, 'max_age', currAge)
    writeFile(lines, runParameters.inlistFile) # write new file
    return
    
# Edits inlist_header to point to your working dir
def setupMesaInlistPaths(runParameters):    
    mesaInlist = runParameters.pathMesa +'/inlist'    
    inspiralInlist = "'" + runParameters.inspiralPath + '/' + runParameters.inlistFile + "'"
    lines = getTxtArrays(mesaInlist)
    
    #Gets indicies of the inlist file header files...
    lines   = getTxtArrays(mesaInlist)
    lines   = replaceLines(lines, 'extra_star_job_inlist_name(1)', inspiralInlist)
    lines   = replaceLines(lines, 'extra_eos_inlist_name(1)', inspiralInlist)
    lines   = replaceLines(lines, 'extra_kap_inlist_name(1)', inspiralInlist)
    lines   = replaceLines(lines, 'extra_controls_inlist_name(1)', inspiralInlist)
    lines   = replaceLines(lines, 'extra_controls_inlist_name(1)', inspiralInlist)
    writeFile(lines, mesaInlist)
    return
# Edits inlist_project to set up file names and point to inspiral dir....
def setupInspiralInlistPaths(runParameters):
    
    lines = getTxtArrays(runParameters.inlistFile)
    lines = replaceLines(lines,'load_model_filename',
                         "'" + runParameters.inspiralPath + '/' + runParameters.modFile + "'")
    lines = replaceLines(lines,'save_model_filename',
                         "'" + runParameters.inspiralPath + '/' + runParameters.finalMod + "'")
    lines = replaceLines(lines,'filename_for_profile_when_terminate',
                         "'" + runParameters.inspiralPath + '/' + runParameters.profile + "'")
    lines = replaceLines(lines,'extra_terminal_output_file',
                         "'" + runParameters.inspiralPath + '/' + runParameters.terminalOutput + "'")    
    lines = replaceLines(lines,'star_history_name',
                         "'" +  runParameters.profile+'.history'  + "'")
    writeFile(lines, runParameters.inlistFile)
    return


def initInspiralFiles(runParameters, donorStar):
    from src.io import output as ceio

    # Create directories to save data
    pathAccretor = ceio.createDir(runParameters.saveDataPath, 
                                  runParameters.saveAccretorName)
    pathDonor    = ceio.createDir(runParameters.saveDataPath, 
                                  runParameters.saveDonorName)
    print("Finished Creating Accretor and Donor star save directories...")

    # Setting Up name prameters
    runParameters.pathDonor    = pathDonor +'/'
    runParameters.pathAccretor = pathAccretor + '/'

    # Sets up general profile and mod files to match star
    runParameters.profile, runParameters.modFile, runParameters.finalMod = [str(round(donorStar.Mass, 2)) + 'Msun' + string for string in ['Profile.data', 'Input.mod', 'Output.mod']]

    # Copyinig starting mod and mesa profiles
    copyFile(runParameters.starModelPath + runParameters.starModModel,
             runParameters.inspiralPath, runParameters.modFile)
    copyFile(runParameters.starModelPath + runParameters.starProfileModel, 
             runParameters.inspiralPath, runParameters.profile)
    print("Finished copying inital MESA mod and profile files...")
    return

def initializeMesaInlist(runParameters, simParameters, donorStar):
    # MESA Working dir. reads inlist in this directory
    setupMesaInlistPaths(runParameters)
    print("MESA reading {} inlist file..".format(runParameters.inlistFile))

    # Inlist in thie directory has file/path names set up
    setupInspiralInlistPaths(runParameters)
    print("MESA inlist file and path names set...")

    # Setup Inlist current age match star profile current age
    setInitialInlistAge(runParameters, simParameters, donorStar.starAge)
    print("MESA inlist current age set to match MESA star profile...")
    
    return
    
# =========================================================================================================