#!/usr/bin/env python
#encoding: UTF-8

"""
This is an example of a quasi-harmonic calculation using the Murnaghan EOS.
"""
    
if __name__ == "__main__":

    from pyqha import RY_KBAR
    from pyqha import gen_TT, read_dos_geo, compute_thermo_geo, read_thermo, rearrange_thermo, fitFvibV, write_xy
    from pyqha import simple_plot_xy, multiple_plot_xy, print_eos_data

    # this part is for calculating the thermodynamic properties from the dos
    fdos="dos_files/output_dos.dat.g"	# base name for the dos files (numbers will be added as postfix)
    fthermo = "thermo"		# base name for the output files (numbers will be added as postfix)

    ngeo = 9	# this is the number of volumes for which a dos has been calculated			

    TT =gen_TT(1,1000)	# generate the numpy array of temperatures for which the properties will be calculated
    T, Evib, Fvib, Svib, Cvib, ZPE, modes = compute_thermo_geo(fdos,fthermo,ngeo,TT)
    nT = len(T)

    
    # Alternatively, read the thermodynamic data from files, if you have already
    # done the calculations. Uncomment the following 2 lines and delete the previous 3 lines
    #T1, Evib1, Fvib1, Svib1, Cvib1 = read_thermo( fthermo, ngeo )
    #T, T, Evib, Fvib, Svib, Cvib = rearrange_thermo( T1, Evib1, Fvib1, Svib1, Cvib1, ngeo )    

    
    fEtot = "./Etot.dat"
    thermodata = nT, T, Evib, Fvib, Svib, Cvib 
    TT, Fmin, Vmin, B0, betaT, Cv, Cp, aT, chi = fitFvibV(fEtot,thermodata)
    
    fig1 = simple_plot_xy(TT,Fmin,xlabel="T (K)",ylabel="Fmin (Ry/cell)")
    fig2 = simple_plot_xy(TT,Vmin,xlabel="T (K)",ylabel="Vmin (a.u.^3)")
    fig3 = simple_plot_xy(TT,B0,xlabel="T (K)",ylabel="B0 (kbar)")
    fig4 = simple_plot_xy(TT,betaT,xlabel="T (K)",ylabel="beta")
    fig5 = simple_plot_xy(TT,Cp,xlabel="T (K)",ylabel="Cp (Ry/cell/K")
    fig1.savefig("figure_1.png")
    fig2.savefig("figure_2.png")
    fig3.savefig("figure_3.png")
    fig4.savefig("figure_4.png")
    fig5.savefig("figure_5.png")    

    # save the results in a file if you want...
    write_xy("Fmin.dat",T,Fmin,"T (K)","Fmin (Ryd/cell)")
    write_xy("Vmin.dat",T,Vmin,"T (K)","Vmin (a.u.^3)")
    write_xy("B0.dat",T,B0*RY_KBAR,"T (K)","B0 (kbar)")
    write_xy("beta.dat",T, betaT,"T (K)","Beta=1/V dV/dT (1/K)")
    
    import numpy as np
    CvCp = np.zeros((len(T),2))
    CvCp[:,0] = Cv
    CvCp[:,1] = Cp
    fig6 = multiple_plot_xy(TT,CvCp,xlabel="T (K)",ylabel="Cv/Cp (Ry/cell/K")
    fig6.savefig("figure_6.png")

    #print_eos_data(V,E+Fvib[i],a,chi,"E")  # print full detail at each T

