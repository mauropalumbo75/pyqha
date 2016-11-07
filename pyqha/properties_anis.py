#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

import numpy as np
from scipy import interpolate 
from constants import RY_KBAR
from fitutils import fit_anis
from minutils import fquadratic, fquartic

################################################################################ 

def compute_volume(celldms,ibrav=4):
    """
    Compute the volume given the *celldms*. Only for ibrav=4 for now, else 
    returns 0.
    """
    if ibrav==4:
        return 0.866025404*celldms[0]*celldms[0]*celldms[2]

    return 0

################################################################################
# 

def compute_alpha(minT,ibrav):
    """
    This function calculates the thermal expansion alphaT at different temperatures
    from the input minT matrix by computing the numerical derivatives with numpy.
    The input matrix minT has shape nT*6, where the first index is the temperature 
    and the second the lattice parameter. For example, minT[i,0] and minT[i,2] are
    the lattice parameters a and c at the temperature i.
 
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

def compute_alpha_splines(TT,minT,ibrav,splinesoptions):
    """
    This function calculates the thermal expansions alphaT at different temperatures
    as the previous function but using spline interpolation as implemented in
    scipy.interpolate.

    """
    alphaT = np.zeros(len(TT)*6)
    alphaT.shape = (len(TT),6)
    
    x = np.array(TT)
    y0 = np.array(minT[:,0])
    y1 = np.array(minT[:,1])
    y2 = np.array(minT[:,2])
        
    if (splinesoptions=={}):
        tck0 = interpolate.splrep(x, y0)
        tck1 = interpolate.splrep(x, y1)
        tck2 = interpolate.splrep(x, y2)  
    else:
        tck0 = interpolate.splrep(x, y0, k=splinesoptions['k0'], s=splinesoptions['s0'])
        tck1 = interpolate.splrep(x, y1, k=splinesoptions['k1'], s=splinesoptions['s1'])
        tck2 = interpolate.splrep(x, y2, k=splinesoptions['k2'], s=splinesoptions['s2'])        
        
    ynew0 = interpolate.splev(x, tck0, der=0)
    alphaT[:,0] = interpolate.splev(x, tck0, der=1)
    ynew1 = interpolate.splev(x, tck1, der=0)
    alphaT[:,1] = interpolate.splev(x, tck1, der=1)
    ynew2 = interpolate.splev(x, tck2, der=0)
    alphaT[:,2] = interpolate.splev(x, tck2, der=1)
    
    # now normalize the alphaTs properly. It must be different for different ibrav
    # to avoid a divide by 0 error (minT is zero for lattice parameters not defined
    # in the system)
    if ibrav==4:
        alphaT[:,0] = alphaT[:,0]/minT[:,0]
        alphaT[:,2] = alphaT[:,2]/minT[:,2]
        
    return alphaT

################################################################################

def compute_S(min0,celldmsx,T,Svib,ibrav=4,typeSvib="quadratic"):
    """
    This function calculates the entropy as a function of temperature. By definition
    :math:`S = -(dF/dT)_{\epsilon}`. To avoid the numerical derivation, within the quasi-harmonic
    approximation it is better to derive it from fitting the harmonic entropy
    results on the grid :math:`(a,b,c)` at the equilibrium lattic parameters given
    in *min0*. *celldms* is the grid :math:`(a,b,c)`, *Svib* are the harmonic 
    entropies on the grid.
    The procedure is the same as the for the :math:`E_{tot}+F_{vib}` in the quasi-harmonic
    calculation but without the minimization step.
    
    Note: a better way would be to do a full harmonic calculation at exactly *min0*.
    The difference with the above way is usually negligible.
    
    **Important**: the above procedure relies on the quasi-harmonic approximation,
    i.e. on the fact that anharmonic contribution are only due to the change of
    phonon frequencies with the lattice parameters. In reality, this is not the 
    case and the entropy so obtained can only be taken as an approximation of the
    real one.
    """
    
    S = np.zeros(len(T))
    if (ibrav==4):
        for iT in range(0,len(T)): 
            # Fit Svib with a quadratic or quartic polynomial, as for Fvib
            aTtemp, chiTtemp = fit_anis(celldmsx,Svib[iT],ibrav,type=typeSvib, ylabel="Svib")
            
            if typeSvib=="quadratic":
                S[iT] = fquadratic(min0,aTtemp)
            elif typeSvib=="quartic":
                S[iT] = fquartic(min0,aTtemp)
        
        return S
    else:
        return None



def compute_Ceps(min0,celldmsx,T,Cvib,ibrav=4,typeCvib="quadratic"):
    """
    This function calculates the constant strain heat capacity :math:`C_{\epsilon}`
    as a function of temperature. 
    By definition :math:`C_{\epsilon} = -T(dS/dT)_{\epsilon}=-T(d^2F/dT^2)_{\epsilon}`. 
    To avoid the numerical derivation, within the quasi-harmonic
    approximation it is better to derive it from fitting the harmonic heat capacities
    results on the grid :math:`(a,b,c)` at the equilibrium lattic parameters given
    in *min0*. *celldms* is the grid :math:`(a,b,c)`, *Cvib* are the harmonic 
    heat capacity on the grid.
    The procedure is the same as the for the :math:`E_{tot}+F_{vib}` in the quasi-harmonic
    calculation but without the minimization step.
    
    Note: a better way would be to do a full harmonic calculation at exactly *min0*.
    The difference with the above way is usually negligible.
    
    **Important**: the above procedure relies on the quasi-harmonic approximation,
    i.e. on the fact that anharmonic contribution are only due to the change of
    phonon frequencies with the lattice parameters. In reality, this is not the 
    case and the entropy so obtained can only be taken as an approximation of the
    real one.
    """
    
    Ceps = np.zeros(len(T))
    if (ibrav==4):
        for iT in range(0,len(T)): 
            # Fit Svib with a quadratic or quartic polynomial, as for Fvib
            aTtemp, chiTtemp = fit_anis(celldmsx,Cvib[iT],ibrav,type=typeCvib, ylabel="Ceps")
            
            if typeCvib=="quadratic":
                Ceps[iT] = fquadratic(min0,aTtemp)
            elif typeCvib=="quartic":
                Ceps[iT] = fquartic(min0,aTtemp)
        
        return Ceps
    else:
        return None


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


def compute_Csigma(TT,Ceps,minT,alphaT,C,ibrav=4):
    """
    This function calculates the constant strain heat capacity :math:`C_{\sigma}`
    as a function of temperature. 
    By definition :math:`C_{\sigma} = -T(dS/dT)_{\sigma}=-T(d^2F/dT^2)_{\epsilon,\sigma}`. 
    To avoid the numerical derivation, within the quasi-harmonic
    approximation it is better to derive it from fitting the harmonic heat capacities
    results on the grid :math:`(a,b,c)` at the equilibrium lattic parameters given
    in *min0*. *celldms* is the grid :math:`(a,b,c)`, *Cvib* are the harmonic 
    heat capacity on the grid.
    The procedure is the same as the for the :math:`E_{tot}+F_{vib}` in the quasi-harmonic
    calculation but without the minimization step.
    
    Note: a better way would be to do a full harmonic calculation at exactly *min0*.
    The difference with the above way is usually negligible.
    
    **Important**: the above procedure relies on the quasi-harmonic approximation,
    i.e. on the fact that anharmonic contribution are only due to the change of
    phonon frequencies with the lattice parameters. In reality, this is not the 
    case and the entropy so obtained can only be taken as an approximation of the
    real one.
    """
    
    Csigma = np.zeros(len(T))
    C = C / RY_KBAR
    Ctemp = np.zeros(len(TT))
    for i in range(1,len(TT)):
        V = compute_volume(minT[i],ibrav)
        for l in range(0,6):
            for m in range(0,6):
                temp = alphaT[i,l] * C[l,m] * alphaT[i,m]
        Ctemp[i] = V * TT[i] * temp    # this is C_sigma-C_epsilon at a given T
  
    Csigma = Ceps + Ctemp
    return Csigma