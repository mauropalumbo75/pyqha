#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

import sys
import time
from read import read_EtotV, read_Etot
from eos import fit_Murn, print_eos_data, write_Etotfitted
from fitutils import fit_anis
from minutils import find_min

   
def fitEtotV(fin,fout=None):
    """
    This function reads :math:`E(V)` data from the input file *fin*, fits them with a Murnaghan EOS,
    prints the results on the *stdout* and write them in the file "fout".
    It returns the volumes and energies read from the input file, the fitted coefficients 
    of the EOS and the corresponding :math:`\chi^2`.
    """
    
    V, E = read_EtotV(fin)   
    a, cov, chi = fit_Murn(V,E)
    print_eos_data(V,E,a,chi,"Etot")
    if (fout!=None):
        write_Etotfitted(fout,V,E,a,chi,"Etot")
    
    return V, E, a, chi
    

def fitEtot(fin, out=True, ibrav=4,fittype="quadratic",guess=None):
    """
    This function reads the file *fin* containing the energies as a function
    of the lattice parameters :math:`E(a,b,c)` and fits them with a quartic (*fittype="quartic"*) or 
    quadratic (*fittype="quadratic"*) polynomial. Then it finds the minimun energy
    and the corresponding lattice parameters. 
    ibrav is the Bravais lattice, guess is an initial guess for the minimization.
    Depending on ibrav, a different number of lattice parameters is considered.
    It prints fitting results on the screen (which can be redericted to *stdout*)
    if *out=True*.
    It returns the lattice parameters and energies as in the input file *fin*,
    the fitted coefficients of the polynomial, the corresponding :math:`\chi^2`,
    the lattice parameters at the minimum and the minimun energy.
    
    Note: for cubic systems use fitEtotV instead.
    """
    
    # Read the energies 
    celldmsx, Ex = read_Etot(fin,ibrav)

    # Fit and find the minimun at 0 K
    a0, chia0 = fit_anis(celldmsx, Ex, ibrav, out, type=fittype)
    if chia0!=None:
        mincelldms, fmin = find_min(a0, ibrav, type=fittype, guess=guess)
        
    return celldmsx, Ex, a0, chia0, mincelldms, fmin 

