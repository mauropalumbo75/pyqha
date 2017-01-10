#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
This submodule groups all functions relevant for the total and vibrational energies. 
"""

import sys
import numpy as np
from .constants import RY_KBAR, RYDBERG_SI, NA
from .read import read_Etot, read_EtotV, read_thermo, read_dos_geo
from .fitutils import print_polynomial, fit_anis, expand_quadratic_to_quartic
from .minutils import find_min
from .eos import fit_Murn, print_eos_data, calculate_fitted_points, compute_beta, compute_beta_splines, compute_Cv, compute_Cp
from .properties_anis import compute_alpha, compute_alpha_splines
from .write import write_celldmsT, write_alphaT, write_xy
from .thermo import gen_TT, compute_thermo_geo, rearrange_thermo
from .plotutils import simple_plot_xy, multiple_plot_xy


################################################################################

def fitFvib(fEtot,thermodata,ibrav=4,typeEtot="quadratic",typeFvib="quadratic",defaultguess=[0.0,0.0,0.0,0.0,0.0,0.0],\
    method="BFGS",minoptions={},splinesoptions=None):
    """
    
    This function computes quasi-harmonic quantities from the 
    :math:`E_{tot}(a,b,c)+F_{vib}(a,b,c,T)` as a function of temperature with Murnaghan's
    EOS. :math:`E_{tot}(a,b,c)` is read from the *fin* file. :math:`F_{vib}(a,b,c,T)`
    are given in *thermodata* which is a list containing the number of temperatures
    ( *nT* ) for which the calculations are done and the numpy matrices for 
    temperatures, vibrational energy, Helmholtz energy, entropy and
    heat capacity. All these quantities are for each (a,b,c) as in *fin* file. The 
    real number of lattice parameters depends on *ibrav*, for example for 
    hexagonal systems (*ibrav=4*) you have only (a,c) values. *ibrav* identifies
    the Bravais lattice, as in Quantum Espresso.
        
    The function fits :math:`E_{tot}(a,b,c)+F_{vib}(a,b,c,T)` with a quadratic
    or quartic polynomial (as defined by *typeEtot* and *typeFvib* ) at each
    temperature in *thermodata* and then stores the fitted coefficients.    
    Note that you can chose a different polynomial type for fitting :math:`E_{tot}(a,b,c)`
    and :math:`F_{vib}(a,b,c)`. Then it computes the minimun energy :math:`E_{tot}+F_{vib}`
    and the corresponding lattice parameters :math:`(a_{min},b_{min},c_{min})` 
    at each temperature by miniimizing the energy.
    
    It also computes the linear thermal expansion tensor (as a numerical derivative of
    the minimum lattice parameters as a function of temperature (:py:func:`compute_alpha`).
    
    It returns the numpy arrays and matrices containing the temperatures (as in input), the
    minimun energy, minimun lattice parameters, linear thermal expansions. It also
    returns the fitted coefficients and the :math:`\chi^2` for :math:`E_{tot}(a,b,c)` 
    only (at T=0 K) and the fitted coefficients and the :math:`\chi^2` for 
    :math:`E_{tot}(a,b,c)+F_{vib}(a,b,c,T)` at each temperature. 
    
    .. Warning::
       The quantities in *thermodata* are usually obtained from :py:func:`compute_thermo_geo`
       or from :py:func:`read_thermo` and :py:func:`rearrange_thermo`. It is
       important that the order in the total energy file *fin* and the order of
       the thermodynamic data in *thermodata* is the same!  See also *example6* and 
       the tutorial.
       
        
    Advanced input parameters:
    
    *guess*, an initial guess for the minimization. It is a 6 elements list 
    [a,b,c,0,0,0]. 
    
    *method*, the method to be used in the minimization procedure, as in the
    scipy.optimize.minimize. See its documentation for details. Note that the 
    methods which usually gives better results for quasi-harmonic calculations
    are the "BFGS" or Newton-CG". Default is "BFGS".
    
    *minoptions*, a dictionary with additional options for the minimization 
    procedure, as in the scipy.optimize.minimize. See its documentation for details.
    Note the the options are different for different methods.
    
    *splinesoptions*, determines whether to use or not splines to reduce the noise
    on numerical derivatives (thermal expansions). If *splinesoptions*==None, use 
    finete differences for derivatives, else use splines as implemented in
    scipy.interpolate (see documentation). In the latter case, *splinesoptions*
    must be a dictionary. This dictionary must contains the keywords *k0*, *s0*,
    *k1*, *s1*, *k2*, *s2* which are passed to :py:func:`scipy.interpolate.splrep`, 
    one couple for each set of thermal expansions (alpha_xx, alpha_yy, alpha_zz).
    *k* is the order of the spline (default=3), s a smoothing condition (default=None).
    If *splinesoptions=={}* use the default options of :py:func:`scipy.interpolate.splrep`
    Note: use this option with care
       
    """
        
    print ("Minimization method: ",method,"\t options: ",minoptions)
    
    # Read the Etot at 0 K
    celldmsx, Ex = read_Etot(fEtot)
    
    ngeo = len(Ex)   # total number of geometries

    nT, T, Evib, Fvib, Svib, Cvib = thermodata
    
    # Fit and find the minimun at 0 K
    a0temp, chi0 = fit_anis(celldmsx, Ex, ibrav, out=True, type=typeEtot, ylabel="Etot")
    if chi0!=None:
        min0, fmin0 = find_min(a0temp, ibrav, type=typeEtot, guess=defaultguess, method=method,minoptions=minoptions)
     
    # if a mixed type fitting has been used, expand the vector of polynomial coefficients
    if (typeEtot=="quadratic" and typeFvib=="quartic"):
        a0=expand_quadratic_to_quartic(a0temp)
    else:
        a0=a0temp
        
    # Fit and find the minimun at T>0 K for each T 
    minT = np.zeros((nT,6))
    fminT= np.zeros(nT)
    aT = np.zeros((nT,len(a0)))
    chiT= np.zeros(nT)
    minguess = min0
    for i in range(0,nT): 
        print ("####################################################################")
        print ("T=",T[i])
  
        # Fit Fvib with a quadratic or quartic polynomial, then add the fitted 
        # quadratic or quartic Etot at 0 K and find the minimun
        aTtemp, chiT[i] = fit_anis(celldmsx,Fvib[i],ibrav,type=typeFvib, ylabel="Fvib")
        if (chiT[i]>0):
            print (typeEtot+" Etot + "+typeFvib+" Fvib polynomials")
            if (typeEtot=="quartic" and typeFvib=="quadratic"):
                aT[i,:]=expand_quadratic_to_quartic(aTtemp)
                typemix="quartic"       # find the minimum as a quadratic
            else:
                aT[i,:]=aTtemp       
                typemix=typeFvib         # typeFvib and typeFtot are the same
            print_polynomial(a0+aT[i,:])

            if (T[i]>50):   # Tests have shown that for temperatures > 50~K  you get better results 
                            # using as initial guess the minimum at the previous temperature
                minguess=minT[i-1,:]
            minT[i,:], fminT[i] = find_min(a0+aT[i,:],ibrav,type=typemix,guess=minguess,method=method,minoptions=minoptions)            

    # Calculate the thermal expansion(s)
    if (splinesoptions==None):
        alphaT = compute_alpha(minT,ibrav)
    else:
        alphaT = compute_alpha_splines(T,minT,ibrav,splinesoptions)
    
    return T, fminT, minT, alphaT, a0, chi0, aT, chiT


################################################################################

def fitFvibV(fin,thermodata,splinesoptions={},lm_pars={},verbosity="low"):
    """
    This function computes quasi-harmonic quantities from the 
    :math:`E_{tot}(V)+F_{vib}(V,T)` as a function of temperature with Murnaghan's
    EOS. :math:`E_{tot}(V)` is read from the *fin* file. :math:`F_{vib}(V,T)`
    are given in *thermodata* which is a list containing the number of temperatures
    ( *nT* ) for which the calculations are done and the numpy matrices for 
    temperatures, vibrational energy, Helmholtz energy, entropy and
    heat capacity. All these quantities are for each volume as in *fin* file.
        
    The function fits :math:`E_{tot}(V)+F_{vib}(V,T)` with a Murnaghan's EOS
    at each temperature in *thermodata* and then stores the fitted coefficients.
    It also computes the volume thermal expansion as a numerical derivative of
    the minimum volume as a function of temperature (:py:func:`compute_beta`), the
    constant volume heat capacity at the minimum volume at each T
    (:py:func:`compute_Cv`) and the constant pression heat capacity (:py:func:`compute_Cp`).
    
    It returns the numpy 1D arrays containing the temperatures (as in input), the
    minimun energy,lm_pars minimun volume, bulk modulus, volume thermal expansion, constant
    volume and constant pressure heat capacities, one matrix with all fitted 
    coefficients at each T and finally an array with the :math:`\chi^2` at each T.
    
    .. Warning::
       The quantities in *thermodata* are usually obtained from :py:func:`compute_thermo_geo`
       or from :py:func:`read_thermo` and :py:func:`rearrange_thermo`. It is
       important that the order in the total energy file *fin* and the order of
       the thermodynamic data in *thermodata* is the same!  See also *example5* and 
       the tutorial. 
       
    """
    
    V, E = read_EtotV(fin)
    
    a0, cov, chi = fit_Murn(V,E)
   
    print_eos_data(V,E,a0,chi,"E")

    ngeo = len(E)   # total number of geometries
    
    nT, T, Evib, Fvib, Svib, Cvib = thermodata  # import list of thermodynamic data
   
    # Fit and find the minimun at T>0 K for each T                    
    aT = np.zeros((nT,4))    # Murnaghan EOS always has 4 parameters, store one for eact T
    chiT = np.zeros((nT))
    for i in range(0,nT): 
        print ("####################################################################")
        print ("T=",T[i])
  
        a, cov, chi = fit_Murn(V,E+Fvib[i],a0,lm_pars)
        if (verbosity=="high"):
            print_eos_data(V,E+Fvib[i],a,chi,"E+Fvib")  # print full detail at each T
             
        print ("Fmin = {:.10e}".format(a[0])+"\tVmin = {:.10e}".format(a[1])
        +"\tB0 = {:.10e}".format(a[2]*RY_KBAR)+"\t chi^2= {:.10e}".format(chi))
        
        aT [i,:] = a
        chiT[i] = chi
        
      
    # As in Murnaghan EOS, aT[:,0] is Fmin(T), aT[:,1] is Vmin(T), aT[:,2] is B0(T)
    # They are arrays, all function of T
    
    # Compute the volume thermal expansion as numerical derivative of Vmin(T)
    if (splinesoptions=={}):
        betaT = compute_beta(T, aT[:,1])
        print ("Computing beta with finite differences...")
    else:
        betaT = compute_beta_splines(T, aT[:,1],splinesoptions)
        print ("Computing beta with splines...")
    
    # compute Cv(Vmin(T))
    Cv = compute_Cv(T,aT[:,1],V,Cvib)
    
    # compute Cp as Cp = Cv + TV beta^2 B0
    Cp = compute_Cp(T,Cvib[:,4],aT[:,1],aT[:,2],betaT)
    
    return T, aT[:,0], aT[:,1], aT[:,2]*RY_KBAR, betaT, Cv*RYDBERG_SI*NA, Cp*RYDBERG_SI*NA, aT, chiT
    
    
