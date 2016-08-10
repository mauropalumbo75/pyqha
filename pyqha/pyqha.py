#!/usr/bin/env python3
#encoding: UTF-8


import time

from fitEtot import fitEtot, fitEtotV
from fitFvib import fitFvib, fitFvibV

def get_input_parameters():
    """
    Get postprocessing input parameters using argparse.
    """
    import argparse 
    
    parser = argparse.ArgumentParser(description='QE post processing')
    
    parser.add_argument('-what', type=int, nargs='?', default=0, choices=range(0, 4),
    help="""selects what to do:\n
    0  = fit Etot(V) (isotropic), energy data must be as a function of volume\n
    1  = fit Etot anisotropic, energy data must be as a function of lattice parameters\n
    2  = fit Etot(V)+Fvib(V) (isotropic), energy data must be as a function of volume\n
    3  = fit Etot+Fvib anisotropic, energy data must be as a function of lattice parameters\n
    """
        )
    parser.add_argument('-ibrav', type=int, nargs='?', default=4, choices=range(0, 2),
    help="""Bravais lattice type, as in QE:\n
    0  = not implemented \n
    1  = 
    """
        )
    parser.add_argument('-fittypeEtot', type=int, nargs='?', default=0, choices=range(0, 2),
    help="""Bravais lattice type, as in QE:\n
    0  = Murnaghan EOS\n
    1  = quadratic polynomial\n
    2  = quartic polynomial 
    """
        )
    parser.add_argument('-fittypeFvib', type=int, nargs='?', default=0, choices=range(0, 2),
    help="""Bravais lattice type, as in QE:\n
    0  = Murnaghan EOS\n
    1  = quadratic polynomial\n
    2  = quartic polynomial 
    """
        )
        
    default_inputdir = "../tests/Os_iso"
    #default_inputdir = "../tests/Os_aniso"
    parser.add_argument('-inputdir', type=str, nargs='?', default=default_inputdir,
                    help='prefix of files saved by program pw.x')
    parser.add_argument('-outputdir', type=str, nargs='?', default="./",
                    help='directory containing the input data, i.e. the same as in pw.x')
    args = parser.parse_args()
    
    return args


if __name__ == "__main__":
    
    start_time = time.time()
    
    # get the input parameters
    pars = get_input_parameters()
    
    
    if (pars.what==0):   
        fitEtotV(pars.inputdir,pars.outputdir)
        
    elif (pars.what==1):   
        if (pars.fittypeEtot==0):
            print ("Murnaghan fit not possible in the anisotropic case. Set properly fittypeEtot.")
        else:
            fitEtot(pars.inputdir,pars.outputdir,pars.ibrav,pars.fittypeEtot)
 
    elif (pars.what==2):    
        fitFvibV(pars.inputdir,pars.outputdir,pars.fittypeEtot,pars.fittypeFvib)
            
    elif (pars.what==3):   
        if (pars.fittypeEtot==0 or pars.fittypeFvib==0):
            print ("Murnaghan fit not possible in the anisotropic case")
        else:
            fitFvib(pars.inputdir,pars.outputdir,pars.ibrav,pars.fittypeEtot,pars.fittypeFvib)
            
    else:
        print ("Not implemented yet")

  
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")

