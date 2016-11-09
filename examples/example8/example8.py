#!/usr/bin/env python
#encoding: UTF-8

"""
This is an example showing several possible numerical issues when fitting/minimizing quasi-harmonic results.
"""
    
if __name__ == "__main__":

    from pyqha import RY_KBAR
    from pyqha import gen_TT, read_Etot, read_dos_geo, compute_thermo_geo, read_thermo, rearrange_thermo, fitFvib, write_celldmsT, write_alphaT
    from pyqha import simple_plot_xy, plot_Etot, plot_Etot_contour, multiple_plot_xy
    import numpy as np

    # this part is for calculating the thermodynamic properties from the dos
    fdos="dos_files/output_dos.dat.g"	# base name for the dos files (numbers will be added as postfix)
    fthermo = "thermo"		# base name for the output files (numbers will be added as postfix)

    ngeo = 25	# this is the number of volumes for which a dos has been calculated			

    #TT =gen_TT(1,2000)	# generate the numpy array of temperatures for which the properties will be calculated
    #T, Evib, Fvib, Svib, Cvib, ZPE, modes = compute_thermo_geo(fdos,fthermo,ngeo,TT)
    #nT = len(T)
    
    # Alternatively, read the thermodynamic data from files if you have already
    # done the calculations
    T1, Evib1, Fvib1, Svib1, Cvib1 = read_thermo( fthermo, ngeo )
    nT, T, Evib, Fvib, Svib, Cvib = rearrange_thermo( T1, Evib1, Fvib1, Svib1, Cvib1, ngeo )    

    
    fEtot = "./Etot.dat"
    thermodata = nT, T, Evib, Fvib, Svib, Cvib 

    # Fit and minimize with default options, quadratic polynomials for both Etot and Fvib, minimization method="BFGS", minoptions={'gtol': 1e-5}
    #TT, Fmin, celldmsminT, alphaT, a0, chi, aT, chiT = fitFvib(fEtot,thermodata)
    res1 = fitFvib(fEtot,thermodata)
    fig1 = simple_plot_xy(res1[0],res1[2][:,0],xlabel="T (K)",ylabel="a_min (a.u.)")
    fig2 = simple_plot_xy(res1[0],res1[2][:,2],xlabel="T (K)",ylabel="c_min (a.u.)")
    fig3 = simple_plot_xy(res1[0],res1[3][:,0],xlabel="T (K)",ylabel="alpha_xx (1/K)")
    fig4 = simple_plot_xy(res1[0],res1[3][:,2],xlabel="T (K)",ylabel="alpha_zz (1/K)")

    # Fit and minimize with quadratic polynomials for both Etot and Fvib, minimization method="BFGS", increasing gtol
    res2 = fitFvib(fEtot,thermodata,method="BFGS",minoptions={'gtol': 1e-6})
    res3 = fitFvib(fEtot,thermodata,method="BFGS",minoptions={'gtol': 1e-7})
    res4 = fitFvib(fEtot,thermodata,method="BFGS",minoptions={'gtol': 1e-8})

    # plot together the 3 resulting thermal expansions
    y = np.zeros((len(res1[0]),3))
    y[:,0] = res2[2][:,0]	
    y[:,1] = res3[2][:,0]
    y[:,2] = res4[2][:,0]
    fig5 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="a_min (a.u.)",labels=["1e-6","1e-7","1e-8"])
    y[:,0] = res2[2][:,2]	
    y[:,1] = res3[2][:,2]
    y[:,2] = res4[2][:,2]
    fig6 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="c_min (a.u.)",labels=["1e-6","1e-7","1e-8"])
    y[:,0] = res2[3][:,0]	
    y[:,1] = res3[3][:,0]
    y[:,2] = res4[3][:,0]
    fig7 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="alpha_xx (1/K)",labels=["1e-6","1e-7","1e-8"])
    y[:,0] = res2[3][:,2]	
    y[:,1] = res3[3][:,2]
    y[:,2] = res4[3][:,2]
    fig8 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="alpha_zz (1/K)",labels=["1e-6","1e-7","1e-8"])



    # Refit the quadratic polynomial for Etot and quadratic for Fvib, default minimization method="Newton-CG", higher minoptions={'gtol': 1e-7}
    res5 = fitFvib(fEtot,thermodata,method="Newton-CG",minoptions={'gtol': 1e-7})

    # Fit the quartic polynomial for Etot and quadratic for Fvib, default minimization method="Newton-CG", higher minoptions={'gtol': 1e-7}
    res6 = fitFvib(fEtot,thermodata,method="Newton-CG",typeEtot="quartic",minoptions={'gtol': 1e-7})

    # Fit the quartic polynomial for Etot and quartic for Fvib, default minimization method="Newton-CG", higher minoptions={'gtol': 1e-7}
    res7 = fitFvib(fEtot,thermodata,method="Newton-CG",typeEtot="quartic",typeFvib="quartic",minoptions={'gtol': 1e-7})

    # plot together the 3 resulting lattice parameters and thermal expansions
    y[:,0] = res5[2][:,0]	
    y[:,1] = res6[2][:,0]
    y[:,2] = res7[2][:,0]
    fig9 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="a_min (a.u.)",labels=["quad+quad fit","quart+quad fit","quart+quart fit"])
    y[:,0] = res5[2][:,2]	
    y[:,1] = res6[2][:,2]
    y[:,2] = res7[2][:,2]
    fig10 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="c_min (a.u.)",labels=["quad+quad fit","quart+quad fit","quart+quart fit"])
    y[:,0] = res5[3][:,0]	
    y[:,1] = res6[3][:,0]
    y[:,2] = res7[3][:,0]
    fig11 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="alpha_xx (1/K)",labels=["quad+quad fit","quart+quad fit","quart+quart fit"])
    y[:,0] = res5[3][:,2]	
    y[:,1] = res6[3][:,2]
    y[:,2] = res7[3][:,2]
    fig12 = multiple_plot_xy(res1[0],y,xlabel="T (K)",ylabel="alpha_zz (1/K)",labels=["quad+quad fit","quart+quad fit","quart+quart fit"])


    # The minimum shifts with temperature, so does the quality of the fit (for example the chi^2)
    celldmsx, Ex = read_Etot(fEtot)  # since the fitFvib does not return Etot data, you must read them from the original file
    iT=1                  # this is the index of the temperatures array, not the temperature itself
    print("T= ",res7[0][iT]," (K)")
    # 2D contour plot with fitted energy (Etot+Fvib)
    fig13 = plot_Etot_contour(celldmsx,nmesh=(50,0,50),fittype="quartic",ibrav=4,a=res7[4]+res7[6][iT]) 
    iT=1998                  # this is the index of the temperatures array, not the temperature itself
    print("T= ",res7[0][iT]," (K)")
    # 2D contour plot with fitted energy (Etot+Fvib)
    fig14 = plot_Etot_contour(celldmsx,nmesh=(50,0,50),fittype="quartic",ibrav=4,a=res7[4]+res7[6][iT]) 
    # plot the chi^2 as a function of temperature
    fig15 = simple_plot_xy(res7[0],res7[5]+res7[7],xlabel="T (K)",ylabel="chi^2")   


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
    fig12.savefig("figure_12.png")
    fig13.savefig("figure_13.png")
    fig14.savefig("figure_14.png")
    fig15.savefig("figure_15.png")

