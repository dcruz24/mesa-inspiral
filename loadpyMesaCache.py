import pyMesa as pym
import os
"""
    Description: Contains function that uses pyMesa load MESA modules and creates a local cache
                 making it easier to faster access. Creating the cache for the first time may 
                 take a while depending on module size.
        
    Note: Can't remeber if there is a line of code that needs to be adjusted in the pyMesa source
          code. Maybe not and just requires pym.loadMod() to be given desired cache directory.    
"""
def loadpyMesaCache():
    cacheDir1            = os.getcwd() + '/pyMesaStarCache'
    print("pyMesa Cache Directory: {}".format(cacheDir1))
    
    # These are needed for energy deposit functions
    eos_lib, eos_def     = pym.loadMod("eos", cacheDir1)
    print("Loaded MESA EOS modules...")

    const_lib, const_def = pym.loadMod("const", cacheDir1)
    print("Loaded MESA const modules...")

    math_lib, _          = pym.loadMod("math", cacheDir1)
    print("Loaded MESA math modules...")
    
    chem_lib,chem_def    = pym.loadMod("chem", cacheDir1)
    print("Loaded MESA chem modules...")
    
    
    #_, star_data_def     = pym.loadMod("star_data", cacheDir1)
    #star_lib, star_def   = pym.loadMod("star", cacheDir1)
    #rse, _               = pym.loadMod('run_star_extras', cacheDir1)
    #star, _              = pym.loadMod("run_star_support", cacheDir1)
    # atm_lib, atm_def    = pym.loadMod("atm",cacheDir1)
    # col_lib,col_def     = pym.loadMod("colors",cacheDir1)
    #kap_lib,kap_def      = pym.loadMod("kap",cacheDir1)
    #net_lib, net_def     = pym.loadMod("net",cacheDir1)
    #rates_lib, rates_def = pym.loadMod("rates",cacheDir1)
    #ion_lib, ion_def     = pym.loadMod("ionization",cacheDir1)
    #neu_lib,neu_def      = pym.loadMod("neu",cacheDir1)
    #utils_lib, utils_def = pym.loadMod("utils",cacheDir1)
    return
    
# Uncomment below..
loadpyMesaCache()
