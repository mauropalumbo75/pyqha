# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
from constants import RY_KBAR

################################################################################ 

def compute_volume(celldms,ibrav=4):
    """
    Compute the volume given the celldms, only for ibrav=4 for now
    """
    if ibrav==4:
        return 0.866025404*celldms[0]*celldms[0]*celldms[2]


################################################################################
# 

def compute_alpha(minT,ibrav):
    """
    This function calculate the thermal expansion alphaT at different temperatures
    from the input minT matrix by computing the numerical derivatives with numpy.
    The input matrix minT has shape nT*6, where the first index is the temperature 
    and the second the lattice parameter. For example, minT[i,0] and minT[i,2] are
    the lattice parameters a and c at the temperature i
 
    More ibrav types must be implemented
    """
    
    grad=np.gradient(np.array(minT))  # numerical derivatives with numpy
    alphaT = np.array(grad[0])  # grad[0] contains the derivatives with respect to T, which is the first axis in minT
                                # also convert to np.array format
    
    # now normalize the alpha properly. It must be different for different ibrav
    # to avoid a divide by 0 error (minT is zero for lattice parameters not defined
    # in the system)
    if ibrav==4:
        alphaT[:,0] = alphaT[:,0]/minT[:,0]
        alphaT[:,2] = alphaT[:,2]/minT[:,2]
        
    return alphaT


################################################################################

def compute_alpha_splines(TT,minT,ibrav):
    """
    This function calculate the thermal expansion alphaT at different temperatures
    as the previous function but with splines

    """
    alphaT = np.zeros(len(TT)*6)
    alphaT.shape = (len(TT),6)
    
    x = np.array(TT)
    y = np.array(minT[:,0])
    tck = interpolate.splrep(x, y, k=3, s=5E-11)
    ynew = interpolate.splev(x, tck, der=0)
    alphaT[:,0] = interpolate.splev(x, tck, der=1)
    
    y = np.array(minT[:,2])
    tck = interpolate.splrep(x, y, k=3, s=2E-10)
    ynew = interpolate.splev(x, tck, der=0)
    alphaT[:,2] = interpolate.splev(x, tck, der=1)
    
    # now normalize the alpha properly. It must be different for different ibrav
    # to avoid a divide by 0 error (minT is zero for lattice parameters not defined
    # in the system)
    if ibrav==4:
        alphaT[:,0] = alphaT[:,0]/minT[:,0]
        alphaT[:,2] = alphaT[:,2]/minT[:,2]
        
    return alphaT

################################################################################

def compute_heat_capacity(TT,minT,alphaT,C,ibrav=4):
    """
    This function calculate the difference between the constant stress heat capacity
    C_sigma and the constant strain heat capacity C_epsilon from the V, the thermal
    expansions and the elastic constant tensor C
    """
    C = C / RY_KBAR
    CT = np.zeros(len(TT))
    for i in range(1,len(TT)):
        V = compute_volume(minT[i],ibrav)
        for l in range(0,6):
            for m in range(0,6):
                temp = alphaT[i,l] * C[l,m] * alphaT[i,m]
        CT[i] = V * TT[i] * temp    # this is C_sigma-C_epsilon at a given T
  
    return CT
