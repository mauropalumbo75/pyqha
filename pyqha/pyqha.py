#!/usr/bin/env python3
#encoding: UTF-8


import time

from fitEtot import fitEtot, fitEtotV
from fitFvib import fitFvib, fitFvibV
from fitC import fitCT
from alphagruneisenp import compute_alpha_gruneisein


def get_input_parameters():
    """
    Get postprocessing input parameters using argparse.
    """
    import argparse 
    
    parser = argparse.ArgumentParser(description='Python quasi-harmonic postprocessing')
    
    parser.add_argument('-what', type=int, nargs='?', default=0, choices=range(0, 7),
    help="""selects what to do:\n
    0  = fit Etot(V) (isotropic), energy data must be as a function of volume\n
    1  = fit Etot anisotropic, energy data must be as a function of lattice parameters\n
    2  = fit Etot(V)+Fvib(V) (isotropic), energy data must be as a function of volume\n
    3  = fit Etot+Fvib anisotropic, energy data must be as a function of lattice parameters\n
    4  = fit C (elastic constant tensor) for an anisotropic system as a function of T. 
        Not implemented yet.\n
    5  = fit C (elastic constant tensor) for an anisotropic system as a function of T\n
    6  = derive the anisotropic thermal expansions from the Gruneisen parameters
        (additional parameters are -alpha_grun and -nproc)\n
    """
        )
        
    parser.add_argument('-alpha_grun', type=int, nargs='?', default=0, choices=range(0, 3),
    help="""Compute the anisotropic thermal expansions from the Gruneisen parameters
    using:\n
    0  = V, gamma and S at 0 K\n
    1  = V(T), gamma(T) but S at 0 K\n
    2  = V(T), gamma(T) and S(T) (the latter calculated as for what=5
    """
        )
        
    parser.add_argument('-nproc', type=int, nargs='?', default=1, choices=range(1, 9),
    help="""Number of processor (1-8) to be used for parallel processing. 
    Currently implemented for what=6\n
    """
        )
        
    parser.add_argument('-ibrav', type=int, nargs='?', default=4, choices=range(1, 5),
    help="""Bravais lattice type, as in QE:\n
    0  = not implemented \n
    1  = simple cubic\n
    2  = bcc\n
    3  = fcc\n
    4  = hcp\n
    """
        )
        
    parser.add_argument('-fittypeEtot', type=int, nargs='?', default=1, choices=range(0, 3),
    help="""Type of fitting function to use for Etot:\n
    0  = Murnaghan EOS\n
    1  = quadratic polynomial\n
    2  = quartic polynomial 
    """
        )
        
    parser.add_argument('-fittypeFvib', type=int, nargs='?', default=1, choices=range(0, 3),
    help="""Type of fitting function to use for Fvib:\n
    0  = Murnaghan EOS\n
    1  = quadratic polynomial\n
    2  = quartic polynomial 
    """
        )
    
    parser.add_argument('-fittypefreq', type=int, nargs='?', default=1, choices=range(1, 3),
    help="""Type of fitting function to use for the frequencies:\n
    1  = quadratic polynomial\n
    2  = quartic polynomial 
    """
        )
        
    parser.add_argument('-fittypeC', type=int, nargs='?', default=1, choices=range(1, 3),
    help="""Type of fitting function to use for the elastic constants (or compliances):\n
    1  = quadratic polynomial\n
    2  = quartic polynomial 
    """
        )
        
    parser.add_argument('-Tmin', type=float, nargs='?', default=1, 
    help="""Min temperature
    """
        )
        
    parser.add_argument('-Tmax', type=float, nargs='?', default=100, 
    help="""Max temperature
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
            
    elif (pars.what==4):   
        print ("Not implemented yet")  
                         
    elif (pars.what==5):   
        if (pars.fittypeEtot==0 or pars.fittypeFvib==0):
            print ("Murnaghan fit not possible in the anisotropic case")
        else:
            fitCT(pars.inputdir,pars.outputdir,pars.ibrav,pars.fittypeEtot,pars.fittypeFvib,pars.fittypeC)
            
    elif (pars.what==6):   
        if (pars.fittypeEtot==0 or pars.fittypeFvib==0):
            print ("Murnaghan fit not possible in the anisotropic case")
        else:
            compute_alpha_gruneisein(pars.inputdir,pars.outputdir,
            pars.fittypeEtot,pars.fittypeFvib,pars.fittypeC,pars.fittypefreq,
            pars.ibrav,pars.alpha_grun,range(pars.Tmin,pars.Tmax),pars.nproc)
    else:
        print ("Not implemented yet")
  
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")

