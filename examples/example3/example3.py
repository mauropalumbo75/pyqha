#!/usr/bin/env python
#encoding: UTF-8

"""
This is an example of computing the vibrational energy, free energy, entropy and enthalpy from the phonon DOS.
"""
    
if __name__ == "__main__":

    from pyqha import gen_TT, read_dos, compute_thermo, write_thermo, RY_TO_CMM1

    fin = "./dos.dat"
    fout = "./thermo"  
    
    TT = gen_TT(1,1000,0.5)	# create a numpy array of temperatures from 1 to 1000 step 0.5
    
    E, dos = read_dos(fin)	# read the dos file. It returns the energies and dos values.

    
    T, Evib, Svib, Cvib, Fvib, ZPE, modes = compute_thermo(E/RY_TO_CMM1,dos*RY_TO_CMM1,TT)
    write_thermo(fout,T, Evib, Fvib, Svib, Cvib, ZPE, modes)  
    
    from pyqha import simple_plot_xy, multiple_plot_xy
    # plot the original phonon DOS
    fig1 = simple_plot_xy(E,dos,xlabel="E (Ryd/cell)",ylabel="phonon DOS (Ry/cell)^{-1}")
    fig1.savefig("figure_1.png")
    # create several plots for the thermodynamic quantities computed
    fig2 = simple_plot_xy(T,Evib,xlabel="T (K)",ylabel="Evib (Ry/cell)")
    fig2.savefig("figure_2.png")
    fig3 = simple_plot_xy(T,Fvib,xlabel="T (K)",ylabel="Fvib (Ry/cell)")
    fig3.savefig("figure_3.png")
    fig4 = simple_plot_xy(T,Svib,xlabel="T (K)",ylabel="Svib (Ry/cell/K)")
    fig4.savefig("figure_4.png")
    fig5 = simple_plot_xy(T,Cvib,xlabel="T (K)",ylabel="Cvib (Ry/cell/K)")
    fig5.savefig("figure_5.png")

