#!/usr/bin/env python3
#encoding: UTF-8

import sys
import numpy as np
from read import read_qha_elastic_constants
from read import read_Etot
from write import write_qha_C
from fitutils import fit_anis
from fitutils import print_polynomial
from minutils import find_min

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
    An auxiliary function for fitting the elastic constant elements of Cxx 
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
#   MAIN, for testing the module
################################################################################
#
# sys.argv[1] -> name of the input file (default="output_energy1")
# sts.argv[2] -> ibrav, lattice type as in Quantum Espresso (default=4 hexagonal for now)

if __name__ == "__main__":
  
    # Default command line parameters
#    inputfileEtot = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
#    inputpathCx = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/"  
    
    inputfileEtot = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    inputpathCx = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/"  
    
    ibrav = 4  # default value, to be later changed to 1
    typeCx = "quartic"
        
    # Read and process command line parameters, if any
    print ("Possible options are: -fileEtot -ibrav -pathCx -fittypeCx -guess\n")
    print ("Assuming default values for unspecified options: \n")
    if len(sys.argv)>1:
        for i in range(1,len(sys.argv)):
            if sys.argv[i] == "-fileEtot": inputfileEtot = sys.argv[i+1]
            elif sys.argv[i] == "-ibrav": ibrav = int(sys.argv[i+1])
            elif sys.argv[i] == "-pathCx": inputpathCx = int(sys.argv[i+1])
            elif sys.argv[i] == "-fittypeCx": typeEtot = sys.argv[i+1]
            elif sys.argv[i] == "-guess": guess = sys.argv[i+1]

    aC, chiC = fitC(inputfileEtot, inputpathCx, ibrav, typeCx)

