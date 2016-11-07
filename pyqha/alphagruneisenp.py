#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.


import numpy as np
import math
from constants import RY_KBAR, K_BOLTZMANN_RY, KB1
from fitfreqgrun import freqmingrun
from properties_anis import compute_volume
from fitC import fS

from multiprocessing import Pool

################################################################################
# 
# Try to import a c version (faster) of the function c_qv, if available
# You should install it with "sudo python3 setupmodule.py install" in your 
# system before starting this program to make this module available
#
is_grunc = False
try:
    from grunc import c_qvc  # This is the same routine c_qv implemented in C to speed it up
    is_grunc = True
except:
    print ("Warning: C module c_qvc not imported, using python (slower) version.")
   

################################################################################
# 
 
def c_qv_python(T,omega):
    """
    This function calculates the mode contribution to the heat capacity at a given T
    and omega. A similar (faster) function should be available as C extension.
    """
    #print ("Python c_qv")
    if (T<1E-9 or omega<1E-9):
        return 0.0
    x = omega * KB1 / T 
    expx = math.exp(-x)   # exponential term
    x2 = math.pow(x,2)
    if expx>1E-3:           # compute normally
        return x2*K_BOLTZMANN_RY*expx/math.pow(expx-1.0,2)
    else:                   # Taylor series
        return K_BOLTZMANN_RY*expx* (x/math.pow(x-0.5*math.pow(x,2)+ 
        0.16666666666666667*math.pow(x,3)+0.04166666666666666667*math.pow(x,4),2))

################################################################################
# 
# If available use a c version of the function c_qv, else use the (slower)
# Python version
# 
if is_grunc:
    c_qv=c_qvc
else: 
    c_qv=c_qv_python


# Same as c_qv but no if. Slightly more efficient, roughly a 30% faster
def c_qv2(T,omega):
    x = omega * KB1 / T 
    expx = math.exp(-x)   # exponential term
    x2 = math.pow(x,2)

    return x2*K_BOLTZMANN_RY*expx/math.pow(expx-1.0,2)




################################################################################
# 

def compute_alpha_grun(T,V,S,weights,freq,grun,ibrav=4):    
    """
    This function computes the thermal expansions alpha using the Gruneisein 
    parameters at a given temperature *T*. *V* is the unit cell volume, *S* is
    the elastic compliances matrix in Voigt notation, *freq* and *weights*  are
    the phonon frequencies in a grid of q-point and their corresponding weights.
    *grun* are the Gruneisen parameters. *ibrav* identifies the Bravais lattice.
    
    It implements the following equation:
    
    :math:``
    
    more comments to be added

    First with min0, freq and grun T-independent

    More ibrav types to be implemented
    """
    nq = freq.shape[0]      # total number of q points
    modes = freq.shape[1]   # number of frequency modes
    alpha = np.zeros(6)     # inizializations
    alphaaux = np.zeros(6)
        
    # compute the Cqv*gruneisen terms, weights for each q-point, and sum
    # for each ibrav (crystalline system) proceed in the proper way
    if ibrav ==1:
        for iq in range(0,nq):     
            for mode in range(0,modes):
                alphaaux[0] += c_qv(T,freq[iq,mode]) * weights[iq] * grun[0,iq,mode]
            
        alphaaux[0] = alphaaux[0] / 3.0 
        alphaaux[1] = alphaaux[0]
        alphaaux[2] = alphaaux[0]
            
    if ibrav ==4:
        for iq in range(0,nq):     
            for mode in range(0,modes):
                temp = c_qv(T,freq[iq,mode]) * weights[iq]   # should be quicker with this additional variable
                alphaaux[0] += temp * grun[0,iq,mode]
                alphaaux[2] += temp * grun[2,iq,mode]
            
        alphaaux[0] = alphaaux[0] / 2.0 
        alphaaux[1] = alphaaux[0]

    else:
        print ("Not implemented yet")
            
    # multiply for the elastic compliances
    for i in range(0,6):
        for j in range(0,6):
            alpha[i] += alphaaux[j]*S[i,j]

    alpha = -alpha/V 
            
    return alpha


################################################################################
# 

def compute_alpha_gruneisen_loopparallel(it):
    """
    This function implements the parallel loop where the alpha are computed. It  
    essentially calls the compute_alpha_grun function with the proper parameters
    according to the selected option. "it" is a list with all function parameters
    for the call in pool.map
    
    it[0]== 0 -> use V, S, average frequencies and gruneisen parameters at 0 K 
    it[0]== 1 -> use S at 0 K, calculate V, average frequencies and gruneisen parameters at each T 
    it[0]== 2 ->  calculate S, V, average frequencies and gruneisen parameters at each T 
 
    """
    
    alphaT = []
    if (it[0]==0):  # use V, S, average frequencies and gruneisen parameters at 0 K 
        for i in range(it[1],it[2]+1):
            alpha = compute_alpha_grun(it[3][i],it[4],it[5],it[6],it[7],it[8],it[9])
            print ("T= "+str(it[3][i])+"\nalpha = ", alpha)
            alphaT.append(alpha)   
    elif (it[0]==1): # use S at 0 K, calculate V, average frequencies and gruneisen parameters at each T             
        for i in range(it[1],it[2]+1):                
            V=compute_volume(it[11][i],it[9])
            freq, grun = freqmingrun(it[10],it[11][i],it[12],it[13],it[9],it[14])
            alpha = compute_alpha_grun(it[3][i],V,it[5],it[6],freq,grun,it[9])
            print ("V = ",str(V),"\nT= "+str(it[3][i])+"\nalpha = ", alpha)
            alphaT.append(alpha)   
    elif (it[0]==2): # calculate V, average frequencies, gruneisen parameters ans S at each T             
        for i in range(it[1],it[2]+1):                
            V=compute_volume(it[11][i],it[9])
            S = fS(it[15],it[11][i],it[16])
            S = S * RY_KBAR    # convert elastic compliances in (Ryd/au)^-1
            freq, grun = freqmingrun(it[10],it[11][i],it[12],it[13],it[9],it[14])
            alpha = compute_alpha_grun(it[3][i],V,S,it[6],freq,grun,it[9])
            print ("V = ",str(V),"\nElastic compliaces matrix: \n",S,"\nT= "+str(it[3][i])+"\nalpha = ", alpha)
            alphaT.append(alpha)   
    else:
        exit()
    
      
    return alphaT

################################################################################

def join(partials):
    """
    This function simply puts together the different temperature ranges where alpha
    was calculated in parallel into a single numpy array
    """
    res = []
    for i in range(0,len(partials)):
        for j in range(0,len(partials[i])):
            res.append(partials[i][j])
            
    return np.array(res)


################################################################################
# 

def compute_alpha_gruneisein(TT,ibrav,celldmsx,min0=None,S=None,weights=None,freq=None,\
     grun=None,minT=None,afreq=None,fittypefreq=None,aS=None,fittypeS=None,nproc=1):

    """
    This function computes the thermal expansions at different temperatures from
    the Gruneisen parameters and other parameters. 
    
    *TT* is a numpy array of temperatures for which the thermal expansions are
    computed. *ibrav* identifies the Bravais lattice. *celldmsx* contains the lattice
    parameters of the computed grid :math:`(a,b,c)` as read with :py:func:`read_Etot`.
    Different calculation options are possible, according to the optional 
    input parameters given in input:

    +---------------------------------------+----------------------------------------------------------------------+
    | Input parameters                      | Meaning                                                              | 
    +=======================================+======================================================================+
    | *min0,S,weights,freq,grun*            | Use temperarature-independent (0 K) volume, elastic compliances,     |
    |                                       | average frequencies and Gruneisen parameters. The volume is obtained |
    |                                       | from *min0*, the lattice parameters at the minimun 0 K energy. *S* is|   
    |                                       | the elastic compliances tensor at 0 K calculated at *min0*. *weights*|
    |                                       | ,*freq*,*grun* are the weights, phonon frequencies and Gruneisen     |
    |                                       | parameters at 0 K calculated at *min0*.                              |
    +---------------------------------------+----------------------------------------------------------------------+    
    | *S,weights,minT,afreq,fittypefreq*    | Use temperarature-independent (0 K) elastic compliances but          |    
    |                                       | temperature-dependent volume, average frequencies and Gruneisen      |
    |                                       | parameters. The volume is obtained from *minT*, a numpy array with   |   
    |                                       | the lattice parameters at the minimun free energy at each temperature|  
    |                                       | . *S* is the elastic compliances tensor at 0 K calculated at *min0*. |
    |                                       | *weights* are the weights of the phonon frequencies.  *afreq*,       |
    |                                       | *fittypefreq* are the coefficients of the fitted polynomials for the |
    |                                       | phonon frequencies and the type of polynomial. They are used here to |
    |                                       | compute the average frequencies and Gruneisen parameters at each     |
    |                                       | temperature.                                                         |
    +---------------------------------------+----------------------------------------------------------------------+
    |*weights,afreq,fittypefreq,aS,fittypeS*| Use temperarature-dependent (0 K) volume, elastic compliances,       |
    |                                       | average frequencies and Gruneisen parameters. The volume is obtained |
    |                                       | from *minT*, a numpy array with the lattice parameters at the minimun|   
    |                                       | free energy at each temperature. *weights* are the weights of the    |  
    |                                       | phonon frequencies. *afreq*,                                         |
    |                                       | *fittypefreq* are the coefficients of the fitted polynomials for the |
    |                                       | phonon frequencies and the type of polynomial. They are used here to |
    |                                       | the average frequencies and Gruneisen parameters at each temperature.|
    |                                       | *aS,fittypeS* are the coefficients of the fitted polynomials for the |
    |                                       | elastic compliances and the polynomial type. They are used here to   |
    |                                       | compute the elastic compliances at each temperature as in the quasi-s|
    |                                       |tatic approximation.                                                  |
    +---------------------------------------+----------------------------------------------------------------------+

    By default, all the input parameters in the above table are ==None and if nothing
    is given in input the function will return None without computing anything.
    
    *nproc* is the number of processes to be run in parallel.
    """
    
    if (min0!=None):        # if the lattice parameters at the minimun 0 K energy are given
        V=compute_volume(min0,ibrav)  # compute eq. volume at 0 K
        print ("Using the volume at 0 K = ",str(V))
    
    if (S!=None):           # S at 0 K given in input
        S = S * RY_KBAR     # convert elastic compliances in (Ryd/au)^-1
     
    ############################################################################
    # Beginning of parallel part

    n_T_in_total = len(TT)
    pool = Pool(processes=nproc)
    n_T_per_process = int(n_T_in_total / nproc)
    print ("Computing {} temperatures per process".format(n_T_per_process))
    
    it = []
    startT = [0]
    endT = [startT[0]+n_T_per_process - 1]
    for i in range(1,nproc):
        startT.append(n_T_per_process + startT[i-1])
        endT.append(n_T_per_process + endT[i-1])
    endT[nproc-1] = n_T_in_total-1
    
    print ("Temperature ranges: ")
    print (startT)
    print (endT)
    
    # pack all parameter in "it" for the parallel call to pool.map,
    # according to different options
    if (min0!=None and S!=None and weights!=None and freq!=None and grun!=None):
        for i in range(0,nproc):
            it.append([0,startT[i],endT[i],TT,V,S,weights,freq,grun,ibrav])
    elif (S!=None and weights!=None and minT!=None and afreq!=None and fittypefreq!=None):
        for i in range(0,nproc):
            it.append([1,startT[i],endT[i],TT,None,S,weights,None,None,ibrav,afreq, minT, afreq.shape[0],afreq.shape[1], fittypefreq])
#            it.append([1,startT[i],endT[i],TT,None,S,weights,None,None,ibrav,afreq, minT, freqxx.shape[0],freqxx.shape[1], fittypefreq])
    elif (weights!=None and afreq!=None and fittypefreq!=None and aS!=None and fittypeS!=None):
        for i in range(0,nproc):
            it.append([2,startT[i],endT[i],TT,None,None,weights,None,None,ibrav,afreq, minT, afreq.shape[0],afreq.shape[1], fittypefreq, aS, fittypeS])   
#            it.append([2,startT[i],endT[i],TT,None,None,weights,None,None,ibrav,afreq, minT, freqxx.shape[0],freqxx.shape[1], fittypefreq, aS, fittypeS])       
    else:
        print ("Option not implemented. Exiting...")
        return None
        
    # Call pool.map to distribute the work in parallel. Results are collected
    # in alphaTpartials which is a list with the results in different 
    # temperature ranges
    alphaTpartials = pool.map(compute_alpha_gruneisen_loopparallel,it)

    # Join all partial results in the numpy array alphaT
    alphaT = join(alphaTpartials)
    
    # End of parallel part
    ############################################################################

    return TT, alphaT
         
