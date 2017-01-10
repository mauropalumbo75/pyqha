#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
This submodule groups several functions for calculating isotropic quasi-harmonic
quantities using the Murnaghan EOS.
"""

from .constants import RY_KBAR
from math import pow
import numpy as np
from scipy.optimize import curve_fit
from scipy import interpolate

################################################################################
# Murnaghan EOS functions 
#
# This one is in the format ideal for fitting, not the canonical one in textbooks 
def E_MurnV(V,a0,a1,a2,a3):
    """
    This function implements the Murnaghan EOS (in a form which is best for fitting).
    Returns the energy at the volume *V* using the coefficients *a0,a1,a2,a3* 
    from the equation:
    
    .. math::
       E = a_0 - (a_2*a_1)/(a_3-1.0) V a_2/a_3 ( a_1/V^{a_3})/(a_3-1.0) +1.0 )
    
    """
    res=np.zeros(len(V))
    for i in range(0,len(V)):
        res[i]=a0 - a2*a1/(a3-1.0) + V[i]*a2/a3*( pow(a1/V[i],a3)/(a3-1.0)+1.0 )
    return res

# Other functions
def E_Murn(V,a):
    """
    As :py:func:`E_MurnV` but input parameters are given as a single list 
    *a=[a0,a1,a2,a3]*.
    """
    return a[0] - a[2]*a[1]/(a[3]-1.0) + V*a[2]/a[3]*( pow(a[1]/V,a[3])/(a[3]-1.0)+1.0 )

def P_Murn(V,a):
    """
    As :py:func:`E_MurnV` but input parameters are given as a single list 
    *a=[a0,a1,a2,a3]* and it returns the pressure not the energy from the EOS.
    """
    return a[2]/a[3]*(pow(a[1]/V,a[3])-1.0)

def H_Murn(V,a):
    """ 
    As :py:func:`E_MurnV` but input parameters are given as a single list 
    *a=[a0,a1,a2,a3]* and it returns the enthalpy not the energy from the EOS.
    """
    return E_Murn(V,a)+P_Murn(V,a)*V


################################################################################

def print_eos_data(x,y,a,chi,ylabel="Etot"):   
    """
    Print the data and the fitted results using the Murnaghan EOS. It can be used for
    different fitted quantities using the proper ylabel. ylabel can be "Etot", 
    "Fvib", etc.
    """
    print ("# Murnaghan EOS \t\t chi squared= {:.10e}".format(chi))
    print ("# "+ylabel+"min= {:.10e} Ry".format(a[0])+"\t Vmin= {:.10e} a.u.^3".format(a[1])+"\t B0= {:.10e} kbar".format(a[2]*RY_KBAR)
    +"\t dB0/dV= {:.10e}".format(a[3]))
    print (80*"#")
    print ("# V *a.u.^3)","\t\t",ylabel," (Ry)\t\t",ylabel+"fit"," (Ry)\t\t",ylabel+"-"+ylabel+"fit (Ry)\tP (kbar)")
    for i in range(0,len(y)):
        print ("{:.10e}".format(x[i]),"\t", "{:.10e}".format(y[i])+
        "\t {:.10e}".format(E_Murn(x[i],a))+
        "\t {:.10e}".format(y[i]-E_Murn(x[i],a))+
        "\t {:.10e}".format(P_Murn(x[i],a)*RY_KBAR))

################################################################################

def write_Etotfitted(filename,x,y,a,chi,ylabel="E"): 
    """
    Write in filename the data and the fitted results using the Murnaghan EOS. It can be used for
    different fitted quantities using the proper ylabel. ylabel can be "Etot", 
    "Fvib", etc.
    """
    fout=open(filename, "w")
    fout.write("# Murnaghan EOS \t\t chi squared= {:.10e}".format(chi))
    fout.write("# E0= {:.10e} Ry".format(a[1])+"\t V0= {:.10e} a.u.^3".format(a[1])+"\t B0= {:.10e} kbar".format(a[2]*RY_KBAR)
    +"\t dB0/dV= {:.10e}".format(a[3]))
    fout.write(80*"#")
    print ("# V *a.u.^3)","\t\t",ylabel," (Ry)\t\t",ylabel+"fit"," (Ry)\t\t",ylabel+"-"+ylabel+"fit (Ry)\tP (kbar)")
    for i in range(0,len(y)):
        fout.write("{:.10e}".format(x[i])+"\t"+"{:.10e}".format(y[i])+
        "\t {:.10e}".format(E_Murn(x[i],a))+
        "\t {:.10e}".format(y[i]-E_Murn(x[i],a))+
        "\t {:.10e}".format(P_Murn(x[i],a)*RY_KBAR))

################################################################################
#
def calculate_fitted_points(V,a):
    """
    Calculates a denser mesh of E(V) points (1000) for plotting.
    """
    Vstep = (V[len(V)-1]-V[0])/1000
    Vdense = np.zeros(1000)
    Edensefitted = np.zeros(1000)
    for i in range(0,1000):
        Vdense[i] = V[0] + Vstep*i 
        Edensefitted[i] = E_Murn(Vdense[i],a)
        
    return Vdense, Edensefitted


################################################################################

def fit_Murn(V,E,guess=[0.0,0.0,900/RY_KBAR,1.15],lm_pars={}):
    """
    This is the function for fitting with the Murnaghan EOS as a function of volume only.

    The input variable *V* is an 1D array of volumes, *E* are the corresponding 
    energies (or other analogous quantity to be fitted with the Murnaghan EOS.
    *a*
    
    Note: volumes must be in a.u.^3 and energies in Rydberg.
    
    """
    # reasonable initial guesses for EOS parameters
    if guess[0]==0.0:
        guess[0] = E[len(E) / 2]
    if guess[1]==0.0:
        guess[1] = V[len(V) / 2]

    a, pcov = curve_fit(E_MurnV, V, E, p0=guess, **lm_pars)
    
    chi = 0
    for i in range(0,len(V)):
        chi += (E[i]-E_Murn(V[i],a))**2
    
    return a, pcov, chi



def compute_beta(TT, minT):
    """
    This function computes the volumetric thermal expansion as a numerical
    derivative of the volume as a function of temperature V(T) given in the
    input array *minT*. This array can obtained
    from the free energy minimization which should be done before.
    """
    grad=np.gradient(np.array(minT))  # numerical derivatives with numpy
    betaT = np.array(grad)  # grad contains the derivatives with respect to T
                                # also convert to np.array format    
    return betaT/minT


def compute_beta_splines(TT, minT, splinesoptions={}):
    """
    This function computes the volumetric thermal expansion as a numerical
    derivative of the volume as a function of temperature V(T) given in the
    input array *minT*. This array can obtained
    from the free energy minimization which should be done before.
    This function uses splines for smoothing the derivative.
    """

    betaT = np.zeros(len(TT))

    x = np.array(TT)
    y0 = np.array(minT)

    if (splinesoptions == {}):
        tck0 = interpolate.splrep(x, y0)
    else:
        tck0 = interpolate.splrep(x, y0, k=splinesoptions['k0'], s=splinesoptions['s0'])

    ynew0 = interpolate.splev(x, tck0, der=0)
    betaT = interpolate.splev(x, tck0, der=1)

    return betaT/minT


def compute_Cv(T,Vmin,V,Cvib):
    """
    This function computes the isocoric heat capacity as a function of temperature.
    From *Cvib*, which is a matrix with *Cvib(T,V)* as from the harmonic calculations
    determines the *Cv* at each temperature by linear interpolation between the values
    at the two volumes closest to Vmin(T). Vmin(T) is from the minimization of F(V,T)
    and *V* is the array of volumes used for it.
    Returns *Cv(T)*.
    
    Work in progress... for now it uses all volumes in the interpolation.
    
    """
    Cv = np.zeros(len(T))
    for iT in range(0,len(T)):
        Cv_interpolated = np.interp(Vmin[iT], V, Cvib[iT,:])
        Cv[iT] = Cv_interpolated
        
    return Cv

def compute_Cp(T,Cv,V,B0,beta):
    """
    This function computes the isobaric heat capacity from the equation:
    
    :math:`Cp - Cv = T V beta^2 B0`

    where *Cp,Cv* are the isobaric and isocoric heat capacities respectively,
    *T* is the temperature, *V* the unit cell volume, *beta* the volumetric
    thermal expansion and *B0* the isothermal bulk modulus.
    
    """
    Cp = Cv + T * V * beta * beta * B0
    return Cp
