#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
A collection of functions for computing harmonic quantities from phonon DOS.
"""

from .constants import RY_TO_CMM1, K_BOLTZMANN_RY
from math import tanh, sinh, log, exp
import numpy as np
from .read import read_dos, read_dos_geo
from .write import write_thermo


def dos_integral(E,dos,m=0):
    """
    A function to compute the integral of an input phonon DOS (*dos*) with the 3/8 Simpson method.
    *m* is the moment of the integral, if :math:`m>0` different moments can be calculated.
    For example, with :math:`m=0` (default) it returns the number of modes from the dos, 
    with :math:`m=1` it returns the ZPE. The input energy (*E*) and phonon DOS (*dos*) are expected to be in
    :math:`cm^{-1}`.
    """
    somma = 0.0
    h = 0.5*(E[2]-E[0])
    for j in range(0,len(dos)-3,3):
        somma += 3.0*pow(E[j],m)*dos[j]+3.0*pow(E[j+1],m)*dos[j+1]+2.0*pow(E[j+2],m)*dos[j+2]
        
    return h*somma*3.0/8.0;


def compute_thermo(E,dos,TT):
    """
    This function computes the vibrational energy, Helmholtz energy, entropy and
    heat capacity in the harmonic approximation from the input numpy arrays *E* 
    and *dos* containing the phonon DOS(E). The calculation is done over a set of
    temperatures given in input as a numpy array *TT*.
    It also computes the number of phonon modes obtained from the input DOS (which
    must be approximately equal to :math:`3*N`, with *N* the number of atoms per cell)
    and the ZPE. The input energy and dos are expected to be in 1/cm-1. 
    It returns numpy arrays for the following quantities (in this order):
    temperatures, vibrational energy, Helmholtz energy, entropy, heat capacity.
    Plus it returns the ZPE and number of phonon modes obtained from the input DOS.
    """
    if (len(dos)<3):
        print ("Not enough points in the phonon DOS!")
        return None
    
    ZPE = 0.5*dos_integral(E,dos,1)
    modes = dos_integral(E,dos)
    
    EvibT = np.zeros(len(TT))
    SvibT = np.zeros(len(TT))
    CvibT = np.zeros(len(TT))
    FvibT = np.zeros(len(TT))
    for i in range(0,len(TT)):
        h = 0.5*(E[2]-E[0])
        arg = K_BOLTZMANN_RY*TT[i]
        arg2 = 2.0 * arg
        Evib = 0.0
        Svib = 0.0
        Cvib = 0.0
        for j in range(0,len(dos)-3,3):

            Evib += 3.0*E[j]/tanh(E[j]/(arg2))*dos[j]+\
            3.0*E[j+1]/tanh(E[j+1]/(arg2))*dos[j+1]+\
            2.0*E[j+2]/tanh(E[j+2]/(arg2))*dos[j+2]
        
            Svib += 3.0*(E[j]/arg2/tanh(E[j]/arg2)-log(2.0*sinh(E[j]/arg2)))*dos[j]+\
            3.0*(E[j+1]/arg2/tanh(E[j+1]/arg2)-log(2.0*sinh(E[j+1]/arg2)))*dos[j+1]+\
            2.0*(E[j+2]/arg2/tanh(E[j+2]/arg2)-log(2.0*sinh(E[j+2]/arg2)))*dos[j+2]

            try:  # avoid overflow error for arg very small
                Cvib += 3.0*pow(E[j]/arg,2)/( 4.0*pow(sinh(E[j]/(arg2)),2) )*dos[j]+\
                    3.0*pow(E[j+1]/arg,2)/( 4.0*pow(sinh(E[j+1]/(arg2)),2) )*dos[j+1]+\
                    2.0*pow(E[j+2]/arg,2)/( 4.0*pow(sinh(E[j+2]/(arg2)),2) )*dos[j+2]
            except:
                Cvib += 0.0

        EvibT[i] = h*0.5*Evib*3.0/8.0  #  h is the integration step, 0.5 comes from the equation for E,
                                                    # the factor 3.0/8.0 comes from the Simpson 3/8 rule
        SvibT[i] = h*K_BOLTZMANN_RY*Svib*3.0/8.0
        CvibT[i] = h*K_BOLTZMANN_RY*Cvib*3.0/8.0
    FvibT = EvibT - SvibT * TT

    print ()
    return TT, EvibT, SvibT, CvibT, FvibT, ZPE, modes


def compute_thermo_geo(fin,fout=None,ngeo=1,TT=np.array([1])):
    """
    This function reads the input dos file(s) from *fin+i*, with *i* a number from
    1 to *ngeo* + 1 and computes vibrational energy, Helmholtz energy, entropy and
    heat capacity in the harmonic approximation. Then writes the output on file(s)
    if fout!=None.
    Output file(s) have the following format:

    +------+-----------------+-----------------+-----------------+-----------------+
    | T    | :math:`E_{vib}` | :math:`F_{vib}` | :math:`S_{vib}` | :math:`C_{vib}` | 
    +======+=================+=================+=================+=================+
    | 1    | ...             | ...             | ...             | ...             |
    +------+-----------------+-----------------+-----------------+-----------------+   

    and are names *fout* +1, *fout* +2,... for each geometry.

    Returning values are (len(TT),ngeo) numpy matrices (T,gEvib,gFvib,gSvib,gCvib,gZPE,gmodes) 
    containing the 
    temperatures and the above mentioned thermodynamic functions as for example:
    Fvib[T,geo] -> Fvib at the temperature "T" for the geometry "geo"
    """

    gEvib=np.zeros((len(TT),ngeo))
    gFvib=np.zeros((len(TT),ngeo))
    gSvib=np.zeros((len(TT),ngeo))
    gCvib=np.zeros((len(TT),ngeo))
    gZPE=np.zeros((ngeo))
    gmodes=np.zeros((ngeo))
    for i in range(0,ngeo):
        E, dos = read_dos(fin+str(i+1))        
        T, Evib, Svib, Cvib, Fvib, ZPE, modes = compute_thermo(E/RY_TO_CMM1,dos*RY_TO_CMM1,TT)
        if (fout!=None):
            write_thermo(fout+str(i+1),T, Evib, Fvib, Svib, Cvib, ZPE, modes)
        
        gEvib[:,i]=Evib
        gFvib[:,i]=Fvib
        gSvib[:,i]=Svib
        gCvib[:,i]=Cvib
        gZPE[i]=ZPE
        gmodes[i]=modes
    
    return TT, gEvib, gFvib, gSvib, gCvib, gZPE, gmodes


def gen_TT(Tstart=1,Tend=1000,Tstep=1):
    """
    A simple function to generate a numpy array of temperatures, starting from
    *Tstart* and ending to *Tend* (or the closest *T<Tend* accorinding to the *Tstep* )
    with step *Tstep* .
    """
    Tsteps = int((Tend-Tstart)//Tstep)
    TT = np.zeros(Tsteps)
    for i in range(0,Tsteps):
        TT[i] = Tstart + i*Tstep
        
    return TT


################################################################################
def rearrange_thermo(T,Evib,Fvib,Svib,Cvib,ngeo=1):
    """
    This function just rearranges the order of the elements in the input matrices
    The first index of the returning matrices *X* now gives all geometries at a given
    *T*, i.e. *X[0]* is the vector of the property *X* a *T=T[0,0]* . *X[0,0]* for the first 
    geometry, *X[0,1]* the second geometry and so on.
    """
    Evib2 = np.zeros((len(T[0]),ngeo))
    Fvib2 = np.zeros((len(T[0]),ngeo))
    Svib2 = np.zeros((len(T[0]),ngeo))
    Cvib2 =np.zeros((len(T[0]),ngeo)) 

    for i in range(0,len(T[0])):
        for j in range(0,ngeo):
            Evib2 [i,j] = Evib[j][i]
            Fvib2 [i,j] = Fvib[j][i]
            Svib2 [i,j] = Svib[j][i]
            Cvib2 [i,j] = Cvib[j][i]
    
    return len(T[0]), T[0], Evib2, Fvib2, Svib2, Cvib2

