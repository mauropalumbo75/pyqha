#!/usr/bin/env python3
#encoding: UTF-8

import sys
import time
from read import read_EtotV, read_Etot
from eos import fit_Murn, print_data, calculate_fitted_points, write_Etotfitted
from fitutils import fit_anis
from minutils import find_min


def fitEtotV(fin,fout):
    """
    This function reads E(V) data from a file, fits them with a Murnaghan EOS,
    prints the results on the screen and in the file "Etotfitted.dat"
    """
    inputfileEtot=fin+"/energy_files/Etot.dat"
    outputfileEtot=fout+"Etotfitted.dat"
    
    V, E = read_EtotV(inputfileEtot)   
    a, cov, chi = fit_Murn(V,E)
    print_data(V,E,a,chi,"Etot")
    write_Etotfitted(outputfileEtot,V,E,a,chi,"Etot")
    
    # Plotting using matplotlib
    Vdense, Edensefitted = calculate_fitted_points(V,a)
    import matplotlib.pyplot as plt
    plt.plot(V, E, 'o', label='Etot data', markersize=10)
    plt.plot(Vdense, Edensefitted, 'r', label='Fitted EOS')
    plt.legend()
    plt.xlabel('V (a.u.^3)')
    plt.ylabel('E (a.u.) ')
    plt.show()


def fitEtot(fin,fout,ibrav,fittype):
    """
    This function reads E(celldms) data from a file, fits them with a quartic or 
    quadratic polynomial, prints the results on the screen and in the file "Etotfitted.dat"
    """
    guess=[5.12374914,0.0,8.19314311,0.0,0.0,0.0]
    
    inputfileEtot=fin+"/energy_files/Etot.dat"
    outputfileEtot=fout+"Etotfitted.dat"
    
    # Read the energies 
    celldmsx, Ex = read_Etot(inputfileEtot)

    fittypeEtot="quadratic"
    if (fittype==2):
        fittypeEtot="quartic"
    # Fit and find the minimun at 0 K
    a0, chia0 = fit_anis(celldmsx, Ex, ibrav, out=True, type=fittypeEtot)
    if chia0!=None:
        mina0, fmina0 = find_min(a0, ibrav, type=fittypeEtot, guess=guess)
    
