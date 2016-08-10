#!/usr/bin/env python3
#encoding: UTF-8

import sys
import numpy as np
from read import read_Etot, read_qha_elastic_constants
from write import write_qha_C, write_qha_CT
from fitutils import fit_anis, print_polynomial
from minutils import find_min, fquadratic, fquartic
from fitFvib import fitFvib

################################################################################

def rearrange_Cx(Cx,ngeo):
    """
    This function rearrange the input numpy matrix Cx into an equivalent matrix Cxx
    for the subsequent fitting.
    Cx is a ngeo*6*6 matrix, each Cx[i] is the 6*6 C matrix for a given geometry (i)
    Cxx is a 6*6*ngeo matrix, each Cxx[i][j] is a vector with all values for different
    geometries of the Cij elastic constant matrix element. For example, Cxx[0][0]
    is the vector with ngeo values of the C11 elastic constant and so on.
    """
    Cxx = []
    for i in range(0,6):
        for j in range(0,6):
            temp = []
            for k in range(0,ngeo):
                temp.append(Cx[k][i][j])
            Cxx.append(temp)    
  
    temp = np.array(Cxx)
    temp.shape = (6,6,ngeo)
    return temp


################################################################################

def fitCxx(celldmsx,Cxx,ibrav,out,typeCx):
    """
    This function fits the elastic constant elements of Cxx as a function of the
    grid of lattice parameters.

    It returns a 6*6 matrix, whose element [i,j] is the set of coefficients of the 
    polynomial fit and another 6*6 matrix, whose element [i,j] is the corresponding
    chi squared. If the chi squared is zero, the fitting procedure was NOT succesful
    """
    Ca = []
    Cchi = []
    for i in range(0,6):
        for j in range(0,6):
            Clabel="C"+str(i+1)+str(j+1)
            a, chi = fit_anis(celldmsx, Cxx[i][j], ibrav, out, typeCx, Clabel)
            Ca.append(a)
            Cchi.append(chi)
    
    tempCa = np.array(Ca)
    tempCchi = np.array(Cchi)

    # Reshape the matrixes so that tempCa[i][j] contains the polinomial coefficients
    # for the Cij elastic constant and tempCchi[i][j] the corresponding chi squared
    tempCa.shape = (6,6,len(tempCa[0]))
    tempCchi.shape = (6,6,1)
    
    return tempCa, tempCchi


################################################################################

def fitC(inputfileEtot, inputpathCx, ibrav, typeCx="quadratic"):
    """
    An auxiliary function for fitting the elastic constant elements of C over a 
    grid of values of the lattice parameters.
    """
    # Read the energies (this is necessary to read the celldmsx)
    celldmsx, Ex = read_Etot(inputfileEtot)

    ngeo = len(Ex)
    Cx, Sx = read_qha_elastic_constants(ngeo, inputpathCx)
    
    Cxx = rearrange_Cx(Cx,ngeo)
    write_qha_C(celldmsx, Cxx, ibrav, inputpathCx)
    
    aC, chiC = fitCxx(celldmsx, Cxx, ibrav, True, typeCx)
    
    return aC, chiC


################################################################################

def fitCT(indir, outdir, ibrav, fittypeEtot, fittypeFvib, fittypeC):
    """
    This function calculates the elastic constants tensor as a function of
    temperatature in the quasi-static approximation.
    First the free energy (Etot+Fvib) is minimized to find the equilibrium values
    of the lattice parameters. Then the elastic constants are fitted over a grid
    of lattice parameters values using a quadratic or quartic polynomial. 
    Finally the elastic constants at each T are obtained from the the fitted 
    polynomial at the values of the lattice parameters which minimize the free
    energy at each temperature.
    """
    
    inputfileEtot=indir+"/energy_files/Etot.dat"
    inputfileFvib=indir+"/therm_files/output_therm.dat"
    inputpathCx=indir+"/elastic_constants/"
    
    # Change the fit type from numeric to string 
    typeC="quadratic"
    if (fittypeC==2):
        typeC="quartic" 
    
    # Call the processing function
    TT, minT, fminT = fitFvib(indir,outdir,ibrav,fittypeEtot,fittypeFvib)

    aC, chiC = fitC(inputfileEtot, inputpathCx, ibrav, typeC)

    nT = 0
    CTtemp = []
    # Find the elastic constants
    for T in range(0,len(TT)):
        C = []
        for i in range(0,6):
            Ccol = []
            for j in range(0,6):
                if typeC=="quadratic":
                    Ctemp = fquadratic(minT[T],aC[i][j],ibrav=4)
                elif typeC=="quartic":
                    Ctemp = fquartic(minT[T],aC[i][j],ibrav=4)  
                Ccol.append(Ctemp)
            C.append(Ccol)
        CTtemp.append(C)
        nT += 1
    
    CT = np.array(CTtemp)
    CT.shape = (nT,6,6)     
    
    write_qha_CT(TT,CT,outdir)

