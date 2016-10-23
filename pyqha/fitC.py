#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

import sys
import numpy as np
from read import read_Etot, read_elastic_constants_geo, read_thermo
from write import write_C_geo, write_CT
from fitutils import fit_anis, print_polynomial
from minutils import find_min, fquadratic, fquartic
from fitFvib import fitFvib
from thermo import gen_TT, compute_thermo_geo, rearrange_thermo
from plotutils import simple_plot_xy, multiple_plot_xy

################################################################################

def rearrange_Cx(Cx,ngeo):
    """
    This function rearrange the input numpy matrix *Cx* into an equivalent matrix *Cxx*
    for fitting it.
    *Cx* is a :math:`ngeo*6*6` matrix, each *Cx[i]* is the 6*6 *C* matrix for a given geometry ( *i* )
    *Cxx* is a Lmath:`6*6*ngeo` matrix, each *Cxx[i][j]* is a vector with all values for different
    geometries of the *Cij* elastic constant matrix element. For example, *Cxx[0,0]*
    is the vector with ngeo values of the *C11* elastic constant and so on.
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

def fitCxx(celldmsx,Cxx,ibrav=4,typeC="quadratic"):
    """
    This function fits the elastic constant elements of *Cxx* as a function of the
    grid of lattice parameters :math:`(a,b,c)`. 
    The real number of lattice parameters depends on *ibrav*, for example for 
    hexagonal systems (*ibrav=4*) you have only (a,c) values. *ibrav* identifies
    the Bravais lattice, as in Quantum Espresso.

    It returns a 6*6 matrix, each element *[i,j]* being the set of coefficients of the 
    polynomial fit and another 6*6 matrix, each element *[i,j]* being the corresponding
    :math:`\chi^2`. If the chi squared is zero, the fitting procedure was NOT succesful
    """
    
    # the most general way (in view of possible extensions) to determine the  
    # number of fitting coefficients (for any possible ibrav and typeCx) is to
    # fit one and then use len()
    atemp, chitemp = fit_anis(celldmsx, Cxx[0][0], ibrav, False, typeC, "")      
    Ca = np.zeros((6,6,len(atemp)))
    Cchi = np.zeros((6,6,1))
    
    for i in range(0,6):
        for j in range(0,6):
            Clabel="C"+str(i+1)+str(j+1)
            Ca[i,j], Cchi[i,j] = fit_anis(celldmsx, Cxx[i][j], ibrav, False, typeC, Clabel)
    
    return Ca, Cchi



################################################################################

def fitCT(aC, chiC, T, minT, ibrav=4, typeC="quadratic"):
    """
    This function calculates the elastic constants tensor *CT* as a function of
    temperatature in the quasi-static approximation.
    It takes in input *aC* and *chiC*, the fitted coefficients of the elastic 
    constants as a function of :math:`(a,b,c)` and the corresponding :math:`\chi^2`.
    It also takes in input an array of temperatures *T* and the corresponding
    lattice parameters *minT*, i.e. :math:`(a_{min},b_{min},c_{min})` from a 
    previous quasi-harmonic calculations (as in example6). 
    It also needs in input the Bravais lattice ( *ibrav* ) and the type of polynomial
    ( *typeC* ) used for fitting the input *aC*.
    
    The function uses the coefficients *aC* to compute the elastic tensor at
    each temperature in the array *T* from the corresponding lattice parameters
    :math:`(a_{min},b_{min},c_{min})` in *minT*.
    
    It returns the temperature array and the a matrix *CT* with all the elastic
    tensors at each T ( *CT[i]* is the elastic constants matrix for the 
    temperature *T[i]*)
    
    .. Warning::
       The coefficients *aC* must be the result of fitting the elastic constants
       over the same :math:`(a,b,c)` grid used in the quasi-harmonic calculations
       corresponding to *minT* values! (See example7) 
    """
    
    CT = np.zeros((len(T),6,6))
    # Find the elastic constants
    for iT in range(0,len(T)):
        C = []
        for i in range(0,6):
            Ccol = []
            for j in range(0,6):
                if typeC=="quadratic":
                    Ctemp = fquadratic(minT[iT],aC[i,j],ibrav)
                elif typeC=="quartic":
                    Ctemp = fquartic(minT[iT],aC[i,j],ibrav)  
                Ccol.append(Ctemp)
            C.append(Ccol)
        CT[iT] = C
    
    return T, CT


