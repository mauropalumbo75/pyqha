#!/usr/bin/env python
#encoding: UTF-8

"""
This is an example of a quasi-harmonic calculation using the Murnaghan EOS.
"""
    
if __name__ == "__main__":

    from pyqha import RY_KBAR
    from pyqha import gen_TT, read_Etot, read_dos_geo, compute_thermo_geo, read_thermo, rearrange_thermo, fitFvib, write_celldmsT, write_alphaT
    from pyqha import simple_plot_xy, multiple_plot_xy
    from pyqha import read_elastic_constants_geo, write_C_geo, write_CT, rearrange_Cx, fitCxx, fitCT

    fEtot = "./Etot.dat"
    celldmsx, Ex = read_Etot(fEtot)  # since the fitFvib does not return Etot data, you must read them from the original file
        
    # this part is for calculating the thermodynamic properties from the dos
    fdos="dos_files/output_dos.dat.g"	# base name for the dos files (numbers will be added as postfix)
    fthermo = "thermo"		# base name for the output files (numbers will be added as postfix)

    ngeo = 25	# this is the number of volumes for which a dos has been calculated  

    #TT =gen_TT(1,1000)	# generate the numpy array of temperatures for which the properties will be calculated
    #T, Evib, Fvib, Svib, Cvib, ZPE, modes = compute_thermo_geo(fdos,fthermo,ngeo,TT)
    #nT = len(T)
    
    # Alternatively, read the thermodynamic data from files if you have already
    # done the calculations
    T1, Evib1, Fvib1, Svib1, Cvib1 = read_thermo( fthermo, ngeo )
    nT, T, Evib, Fvib, Svib, Cvib = rearrange_thermo( T1, Evib1, Fvib1, Svib1, Cvib1, ngeo )   
    
    fEtot = "./Etot.dat"
    thermodata = nT, T, Evib, Fvib, Svib, Cvib 
    TT, Fmin, celldmsminT, alphaT, a0, chi, aT, chi = fitFvib(fEtot,thermodata,typeEtot="quartic",typeFvib="quartic",defaultguess=[5.12374914,0.0,8.19314311,0.0,0.0,0.0])


    # Now start the quasi-static calculation    
    fC = "./elastic_constants/output_el_cons.g"
    
    # Read the elastic constants and compliances from files 
    Cx, Sx = read_elastic_constants_geo(fC, ngeo)    
    Cxx = rearrange_Cx(Cx,ngeo)     # rearrange them in the proper order for fitting
    
    # Optionally save them
    write_C_geo(celldmsx, Cxx, ibrav=4, fCout="./elastic_constants/")
    
    # Fit the elastic constants as a function of celldmsx
    aC, chiC = fitCxx(celldmsx, Cxx, ibrav=4,typeC="quadratic")
  
    T, CT = fitCT(aC, chiC, TT, celldmsminT, ibrav=4, typeC="quadratic")
           
    write_CT(TT,CT,fCout="./elastic_constants/")
    
    simple_plot_xy(TT,CT[:,0,0],xlabel="T (K)",ylabel="C11 (kbar)")
    simple_plot_xy(TT,CT[:,0,1],xlabel="T (K)",ylabel="C12 (kbar)")
    simple_plot_xy(TT,CT[:,0,2],xlabel="T (K)",ylabel="C13 (kbar)")
    simple_plot_xy(TT,CT[:,2,2],xlabel="T (K)",ylabel="C33 (kbar)")
    
    # plot now 4 elastic constants in the same plot
    import numpy as np
    pCxx = np.zeros((len(T),4))
    pCxx[:,0] = CT[:,0,0]
    pCxx[:,1] = CT[:,0,1]
    pCxx[:,2] = CT[:,0,2]
    pCxx[:,3] = CT[:,2,2]
    Clabels = ["C11","C12","C13","C33"]
    multiple_plot_xy(T,pCxx,xlabel="T (K)",ylabel="Cxx (kbar)",labels=Clabels)

