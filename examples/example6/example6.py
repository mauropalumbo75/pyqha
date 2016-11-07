#!/usr/bin/env python
#encoding: UTF-8

"""
This is an example of a quasi-harmonic calculation for an hexagonal system using a quadratic polynomial.
"""
    
if __name__ == "__main__":

    from pyqha import RY_KBAR
    from pyqha import gen_TT, read_Etot, read_dos_geo, compute_thermo_geo, read_thermo, rearrange_thermo, fitFvib, write_celldmsT, write_alphaT
    from pyqha import simple_plot_xy, plot_Etot, plot_Etot_contour

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
    TT, Fmin, celldmsminT, alphaT, a0, chi, aT, chi = fitFvib(fEtot,thermodata,minoptions={'gtol': 1e-7})

    fig1 = simple_plot_xy(TT,Fmin,xlabel="T (K)",ylabel="Fmin (Ry/cell)")
    fig2 = simple_plot_xy(TT,celldmsminT[:,0],xlabel="T (K)",ylabel="a_min (a.u.)")
    fig3 = simple_plot_xy(TT,celldmsminT[:,2],xlabel="T (K)",ylabel="c_min (a.u.)")
    fig4 = simple_plot_xy(TT,celldmsminT[:,2]/celldmsminT[:,0],xlabel="T (K)",ylabel="c/a ")    
    fig5 = simple_plot_xy(TT,alphaT[:,0],xlabel="T (K)",ylabel="alpha_xx (1/K)")
    fig6 = simple_plot_xy(TT,alphaT[:,2],xlabel="T (K)",ylabel="alpha_zz (1/K)")

    # write a(T) and c(T) on a file
    write_celldmsT("celldmminT",T,celldmsminT,ibrav=4)
    # write alpha_xx(T) and alpha_zz(T) on a file
    write_alphaT("alphaT",T,alphaT,ibrav=4)
    
    # Plot several quantities at T=998+1 K as an example
    celldmsx, Ex = read_Etot(fEtot)  # since the fitFvib does not return Etot data, you must read them from the original file
    iT=998                  # this is the index of the temperatures array, not the temperature itself
    print("T= ",TT[iT]," (K)")
    # 3D plot only with fitted energy (Etot+Fvib)
    fig7 = plot_Etot(celldmsx,Ex=None,n=(5,0,5),nmesh=(50,0,50),fittype="quadratic",ibrav=4,a=a0+aT[iT])
    # 3D plot fitted energy and points
    fig8 = plot_Etot(celldmsx,Ex+Fvib[iT],n=(5,0,5),nmesh=(50,0,50),fittype="quadratic",ibrav=4,a=a0+aT[iT])
    # 3D plot with fitted energy Fvib only
    fig9 = plot_Etot(celldmsx,Ex=None,n=(5,0,5),nmesh=(50,0,50),fittype="quadratic",ibrav=4,a=aT[iT])
    # 2D contour plot with fitted energy (Etot+Fvib)
    fig10 = plot_Etot_contour(celldmsx,nmesh=(50,0,50),fittype="quadratic",ibrav=4,a=a0+aT[iT])
    # 2D contour plot with fitted energy Fvib only
    fig11 = plot_Etot_contour(celldmsx,nmesh=(50,0,50),fittype="quadratic",ibrav=4,a=aT[iT])   

    # Save all plots
    fig1.savefig("figure_1.png")
    fig2.savefig("figure_2.png")
    fig3.savefig("figure_3.png")
    fig4.savefig("figure_4.png")
    fig5.savefig("figure_5.png")    
    fig6.savefig("figure_6.png")      
    fig7.savefig("figure_7.png") 
    fig8.savefig("figure_8.png")
    fig9.savefig("figure_9.png")
    fig10.savefig("figure_10.png")
    fig11.savefig("figure_11.png")

