#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

import sys
import time
import numpy as np
from read import read_Etot, read_freq_ext_geo
from fitutils import fit_anis, print_polynomial
from minutils import fquadratic, fquartic, fquadratic_der, fquartic_der, find_min
from write import write_freq_ext



################################################################################

def rearrange_freqx(freqx):
    """
    This function rearrange the input numpy matrix freqx into an equivalent matrix freqxx
    for the subsequent fitting.
    freqx is a ngeo*nq*modes matrix, each freqx[i] is the nq*modes freq matrix for a given geometry (i)
    freqxx is a nq*modes*ngeo matrix, each freqxx[i][j] is a vector with all values for different
    geometries of the frequencies at point q=i and mode=j. For example, freqxx[0][0]
    is the vector with ngeo values of the frequencies at the first q-point and first mode so on.
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
    grid of lattice parameters.

    It returns a nq*modes matrix, whose element [i,j] is the set of coefficients of the 
    polynomial fit and another nq*modes matrix, whose element [i,j] is the corresponding
    chi squared. If the chi squared is zero, the fitting procedure was NOT succesful
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
    This function calculate the frequencies from the fitted polynomials at the 
    minimun point min0. afreq contains the fitted polynomial coefficients. 

    It returns a nq*modes matrix, whose element [i,j] is the fitted frequency 
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
    This function calculate the frequencies and the gruneisen parameters 
    from the fitted polynomials at the minimun point min0. afreq contains the 
    fitted polynomial coefficients. 

    It returns a nq*modes matrix, whose element [i,j] is the fitted frequency 
    In addition, it returns a nq*modes*6 with the Gruneisein parameters.
    Each element [i,j,k] is the the Gruneisein parameter at nq=i, mode=j and direction
    k (for example, in hex systems k=0 is a direction, k=2 is c direction, other are zero)

    Note that the Gruneisein parameters are not multiplied for the lattice parameters 
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
 
def fitfreq(celldmsx, min0, inputfilefreq, ibrav=4, typefreq="quadratic", compute_grun=False):
    """
    An auxiliary function for fitting the frequencies. It returns a matrix of nq*modes
    frequencies obtained for the fitted polynomial (quadratic or quartic) at the 
    minimun point min0. It also returns the weigths of each q point where the 
    frequencies are available.
    """
    # get the weigths and the frequencies from files 
    weightsx, freqx = read_freq_ext_geo(inputfilefreq,range(1,celldmsx.shape[0]+1))
        
    print ("Rearranging frequencies...")
    freqxx = rearrange_freqx(freqx)
    print ("Done!")
    del freqx
    
    print ("Fitting frequencies...")
    afreq, chifreq = fitfreqxx(celldmsx, freqxx, ibrav, True, typefreq)
    print ("Done!")  
    
    f, grun = freqmingrun(afreq, min0, freqxx.shape[0],freqxx.shape[1], ibrav, typefreq)
    
    return weightsx[0,:], f, grun	# weights for different geometries are supposed to be the same, return only one


################################################################################
#   MAIN, to be used for testing the module
################################################################################
#

if __name__ == "__main__":
  
    # Default command line parameters
    inputfileEtot = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    inputfilefreq = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/save_frequencies.dat.g"  
    
    #inputfileEtot = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    #inputfilefreq = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/save_frequencies.dat.g"  
    
    ibrav = 4  # default value, to be later changed to 1
    typeEtot = "quartic"
    typefreq = "quadratic"
    guess = [5.13327423, 0.0, 1.57852513, 0.0, 0.0, 0.0]
        
    start_time = time.time()
    
    # Read the energies (this is necessary to read the celldmsx)
    celldmsx, Ex = read_Etot(inputfileEtot)
    # Fit and find the minimun at 0 K
    a0, chia0 = fit_anis(celldmsx, Ex, ibrav, out=True, type=typeEtot)
    if chia0!=None:
        min0, fmin0 = find_min(a0, ibrav, type=typeEtot, guess=guess)
    
    # fit the frequencies over the grid of lattice parameters and get the fitted
    # ones at the minimum point from Etot (min0)
    weights, f, grun = fitfreq(celldmsx, min0, inputfilefreq, ibrav, typefreq="quadratic", compute_grun=True)
    write_freq_ext(weights,f,"average_freqPython")
    write_freq_ext(weights,grun[0],"output_grun_along_a_ext3Dfit")
    write_freq_ext(weights,grun[2],"output_grun_along_c_ext3Dfit")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")
