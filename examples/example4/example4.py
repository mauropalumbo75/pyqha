#!/usr/bin/env python
#encoding: UTF-8

"""
This is an example of computing the vibrational energy, free energy, entropy and enthalpy from several phonon DOS stored in different files. Each
"""
    
if __name__ == "__main__":

    from pyqha import gen_TT, read_dos_geo, compute_thermo_geo
    from pyqha import simple_plot_xy, multiple_plot_xy

    fin = "dos_files/output_dos.dat.g"	# base name for the dos files (numbers will be added as postfix)
    fout = "thermo"		# base name for the output files (numbers will be added as postfix)

    ngeo = 9				

    gE, gdos = read_dos_geo(fin,ngeo) 	# read ngeo=9 dos files


    # plot the first 5 phonon dos    
    multiple_plot_xy(gE[:,0:5],gdos[:,0:5],xlabel="E (Ryd/cell)",ylabel="phonon DOS (cell/Ryd)")

    TT =gen_TT(1,1000)	# generate the numpy array of temperatures for which the properties will be calculated

    # compute the thermodynamic properties for all ngeo dos files and write them in fout+"i" files, where is an int from 1 to ngeo
    T, ggEvib, ggFvib, ggSvib, ggCvib, ggZPE, ggmodes = compute_thermo_geo(fin,fout,ngeo,TT)
    
    # plot the vibrational Helmholtz energy for the first 5 phonon dos    
    multiple_plot_xy(T,ggFvib[:,0:5],xlabel="T (K)",ylabel="Cvib (Ry/cell/K")

    # plot the vibrational entropy for the first 5 phonon dos    
    multiple_plot_xy(T,ggSvib[:,0:5],xlabel="T (K)",ylabel="Cvib (Ry/cell/K")

    # plot the vibrational heat capacity for the first 5 phonon dos    
    multiple_plot_xy(T,ggCvib[:,0:5],xlabel="T (K)",ylabel="Cvib (Ry/cell/K")

