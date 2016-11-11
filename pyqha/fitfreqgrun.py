#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
This submodule groups all functions relevant for fitting phonon frequencies and
computing Gruneisen parameters. 
"""

import sys
import time
import numpy as np
from read import read_Etot, read_freq_geo
from fitutils import fit_anis, print_polynomial
from minutils import fquadratic, fquartic, fquadratic_der, fquartic_der, find_min
from write import write_freq

################################################################################

def rearrange_freqx(freqx):
    """
    This function rearranges the input numpy matrix *freqx* into an equivalent 
    matrix *freqxx* for the subsequent fitting.
    *freqx* is a :math:`ngeo*nq*modes` matrix, each *freqx[i]* is the :math:`nq*modes`
    frequency matrix for a given geometry *i*.
    freqxx is a :math:`nq*modes*ngeo` matrix, each *freqxx[i,j]* is a vector with
    all values for different geometries of the frequencies at point *q=i* and *mode=j*.
    For example, *freqxx[0,0]* is the vector with ngeo values of the frequencies
    at the first q-point and first mode so on.
    """
    ngeo=freqx.shape[0]
    nq=freqx.shape[1]
    modes=freqx.shape[2]
    freqxx = []
    for i in range(0,nq):
        for j in range(0,modes):
            temp = []
            for k in range(0,ngeo):
                temp.append(freqx[k][i][j])
            freqxx.append(temp)    
  
    temp = np.array(freqxx)
    temp.shape = (nq,modes,ngeo)
    return temp


################################################################################

def fitfreqxx(celldmsx,freqxx,ibrav,out,typefreq):
    """
    This function fits the frequencies in freqxx as a function of the
    grid of lattice parameters :math:`(a,b,c)`.

    It returns a :math:`nq*modes` matrix, whose element [i,j] is the set of coefficients of the 
    polynomial fit and another :math:`nq*modes` matrix, whose element [i,j] is the corresponding
    :math:`\chi^2`. If the :math:`\chi^2=0`, the fitting procedure was NOT succesful
    """    
    freqa = []
    freqchi = []
    nq = freqxx.shape[0] # total number of q points read
    modes = freqxx.shape[1]  # number of frequency modes
    #nq=10   # for testing quickly
    count = 0
    for i in range(0,nq):
        for j in range(0,modes):
            freqlabel="f(nq="+str(i+1)+",mode="+str(j+1)+")"
            a, chi = fit_anis(celldmsx, freqxx[i][j], ibrav, False, typefreq, freqlabel)          
            freqa.append(a)
            freqchi.append(chi)
            count += 1
            if (count%1000==0):
                print (str(count)+" frequencies fitted...")
    
    tempfreqa = np.array(freqa)
    tempfreqchi = np.array(freqchi)

    # Reshape the matrixes so that tempCa[i][j] contains the polinomial coefficients
    # for the Cij elastic constant and tempCchi[i][j] the corresponding chi squared
    tempfreqa.shape = (nq,modes,len(tempfreqa[0]))
    tempfreqchi.shape = (nq,modes,1)
    
    return tempfreqa, tempfreqchi


################################################################################

#
def freqmin(afreq, min0, nq, modes, ibrav, typefreq):
    """
    This function calculates the frequencies from the fitted polynomials coeffients (one
    for each q point and mode) at the minimun point *min0* given in input. 
    *afreq* is a :math:`nq*modes` numpy matrix containing the fitted polynomial coefficients.
    It can be obtained from :py:func:`fitfreqxx`.

    It returns a :math:`nq*modes` matrix, each element [i,j] being the fitted frequency 
    """
    f = np.zeros((nq,modes))
    for i in range(0,nq):
        for j in range(0,modes):
            if typefreq=="quadratic":
                f[i,j] = fquadratic(min0,afreq[i][j],ibrav=4)
            elif typefreq=="quartic":
                f[i,j] = fquartic(min0,afreq[i][j],ibrav=4) 
                
    return f

################################################################################

def freqmingrun(afreq, min0, nq, modes, ibrav, typefreq):
    """
    This function calculates the frequencies and the Gruneisen parameters
    from the fitted polynomials coeffients (one
    for each q point and mode) at the minimun point *min0* given in input. 
    *afreq* is a :math:`nq*modes` numpy matrix containing the fitted polynomial coefficients.
    It can be obtained from :py:func:`fitfreqxx`.

    It returns a :math:`nq*modes` matrix, each element [i,j] being the fitted frequency 
    In addition, it returns a :math:`nq*modes*6` with the Gruneisein parameters.
    Each element [i,j,k] is the the Gruneisein parameter at *nq=i*, *mode=j* and direction
    *k* (for example, in hex systems *k=0* is *a* direction, *k=2* is *c* direction, others are zero)

    Note that the Gruneisein parameters are not multiplied for the lattice parameters. 
    """
    #nq = 10 # for testing quickly
    f = np.zeros((nq,modes))
    grun = np.zeros((6,nq,modes))   # initialize the Gruneisein parameters matrix
    for i in range(0,nq):
        for j in range(0,modes):
            if typefreq=="quadratic":
                f[i,j] = fquadratic(min0,afreq[i][j],ibrav=4)
                if f[i,j]<1E-5:
                    grun[:,i,j] = 0.0
                else:
                    grun[:,i,j] = fquadratic_der(min0,afreq[i][j],ibrav)/f[i,j]*min0
            elif typefreq=="quartic":
                f[i,j] = fquartic(min0,afreq[i][j],ibrav=4) 
                if f[i,j]<1E-5:
                    grun[:,i,j] = 0.0
                else:
                    grun[:,i,j] = fquartic_der(min0,afreq[i][j],ibrav)/f[i,j]*min0
    
    return f, grun 


################################################################################
 
def fitfreq(celldmsx, min0, filefreq, ibrav=4, typefreq="quadratic", compute_grun=False):
    """
    An auxiliary function for fitting the frequencies. 
    
    *celldmsx* is the matrix of lattice parameters :math:`(a,b,c)` where the total
    energies where computed. *min0* is the a set of :math:`(a,b,c)`. *filefreq*
    defines the input files (*filefreq1*, *filefreq2*, etc.) containing the 
    frequencies for different geometries. The number of geometries is determined
    from the size of *celldmsx*. *ibrav* is the usual Bravais lattice. 
    *typefreq* can be "quadratic" (default) or "quartic", i.e. the kind of 
    polynomial to be used for fitting. *compute_grun* defines if the Gruneisen
    parameters must be calculated (True) or not (False, default). 
    
    It returns a matrix of :math:`nq*modes` frequencies obtained for the fitted 
    polynomial coefficients (quadratic or quartic) at the 
    minimun point *min0*. It also returns the weigths of each q-point where the 
    frequencies are available.
    """
    # get the weigths and the frequencies from files 
    weightsx, freqx = read_freq_geo(filefreq,celldmsx.shape[0])
        
    print ("Rearranging frequencies...")
    freqxx = rearrange_freqx(freqx)
    print ("Done!")
    del freqx
    
    print ("Fitting frequencies...")
    afreq, chifreq = fitfreqxx(celldmsx, freqxx, ibrav, True, typefreq)
    print ("Done!")  
    
    f, grun = freqmingrun(afreq, min0, freqxx.shape[0],freqxx.shape[1], ibrav, typefreq)
    
    return weightsx[0,:], f, grun	# weights for different geometries are supposed to be the same, return only one

