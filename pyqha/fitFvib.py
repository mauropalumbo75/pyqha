#!/usr/bin/env python3
#encoding: UTF-8

import sys
import numpy as np
from constants import RY_KBAR
from read import read_Etot, read_EtotV, read_thermo
from fitutils import print_polynomial, fit_anis, expand_quadratic_to_quartic
from minutils import find_min
from eos import fit_Murn, print_data,calculate_fitted_points, compute_beta
from properties_anis import compute_alpha, compute_alpha_splines,compute_heat_capacity
from write import write_celldmsT, write_alphaT, write_y_T


################################################################################

def rearrange_thermo(T,Evib,Fvib,Svib,Cvib,ngeo=1):
    """
    This function just rearranges the order of the elements in the input matrices
    The first index of the returning matrices X now gives all geometries at a given
    T, i.e. X[0] is the vector of the property X a T=T[0][0]. X[0][0] for the first 
    geometry, X[0][1] the second geometry and so on.
    """
    T2=[]
    Evib2=[]
    Fvib2=[]
    Svib2=[]
    Cvib2=[]
    #for i in range(0,len(T[0])):
    # stupid test
    for i in range(0,len(T[0])):
        TT=[]
        TEvib=[]
        TFvib=[]
        TSvib=[]
        TCvib=[]
        for j in range(0,ngeo):
            TT.append(T[j][i])
            TEvib.append(Evib[j][i])
            TFvib.append(Fvib[j][i])
            TSvib.append(Svib[j][i])
            TCvib.append(Cvib[j][i])
      
        T2.append(TT)
        Evib2.append(TEvib)
        Fvib2.append(TFvib)
        Svib2.append(TSvib)
        Cvib2.append(TCvib)
    
    return len(T2), np.array(T2), Evib2, Fvib2, Svib2, Cvib2


################################################################################

def fitFvib(fin,fout,ibrav,fittypeEtot,fittypeFvib,defaultguess=[0.0,0.0,0.0,0.0,0.0,0.0]):
    """
    An auxiliary function for handling the qha calculations with Etot+Fvib    
    """
    inputfileEtot=fin+"/energy_files/Etot.dat"
    inputfiledos=fin+"/dos_files/output_therm.dat"
    inputfileFvib=fin+"/therm_files/output_therm.dat"
    
    # Change the fit type from numeric to string
    typeEtot="quadratic"
    if (fittypeEtot==2):
        typeEtot="quartic"
    typeFvib="quadratic"
    if (fittypeFvib==2):
        typeFvib="quartic"   
        
    # Read the Etot at 0 K
    celldmsx, Ex = read_Etot(inputfileEtot)
    
    ngeo = len(Ex)   # total number of geometries
    
    #compute_thermo_files()    # to be implemented
    
    # Read the vibrational energies from the directory "therm_files" if present
    T1, Evib1, Fvib1, Svib1, Cvib1 = read_thermo( inputfileFvib, ngeo )
    nT, T, Evib, Fvib, Svib, Cvib = rearrange_thermo( T1, Evib1, Fvib1, Svib1, Cvib1, ngeo )
    
    # Fit and find the minimun at 0 K
    a0temp, chi0 = fit_anis(celldmsx, Ex, ibrav, out=True, type=typeEtot, ylabel="Etot")
    if chi0!=None:
        min0, fmin0 = find_min(a0temp, ibrav, type=typeEtot, guess=defaultguess)
     
    # if a mixed type fitting has been used, expand the vector of polynomial coefficients
    if (typeEtot=="quadratic" and typeFvib=="quartic"):
        a0=expand_quadratic_to_quartic(a0temp)
    else:
        a0=a0temp
        
    # Fit and find the minimun at T>0 K for each T                    
    minT=[]
    fminT=[]
    for i in range(0,nT): 
        print ("####################################################################")
        print ("T=",T[i][0])
  
        # Fit Fvib with a quadratic or quartic polynomial, then add the fitted 
        # quadratic or quartic Etot at 0 K and find the minimun
        aTtemp, chiT = fit_anis(celldmsx,Fvib[i],ibrav,type=typeFvib, ylabel="Fvib")
        if (chiT>0):
            print (typeEtot+" Etot + "+typeFvib+" Fvib polynomials")
            if (typeEtot=="quartic" and typeFvib=="quadratic"):
                aT=expand_quadratic_to_quartic(aTtemp)
                typemix="quartic"       # find the minimum as a quadratic
            else:
                aT=aTtemp       
                typemix=typeFvib         # typeFvib and typeFtot are the same
            print_polynomial(a0+aT)

            minTtemp, fminTtemp = find_min(a0+aT,ibrav,type=typemix,guess=min0)
                
            minT.append(minTtemp)
            fminT.append(fminTtemp)  


    TT = T[:,0]
    minT = np.array(minT)
    fminT = np.array(fminT)
    write_celldmsT(fout+"/celldmsT.dat",TT,minT,ibrav)
   
    # Calculate the thermal expansion(s)
    alphaT = compute_alpha(minT,ibrav)
    #alphaT = compute_alpha_splines(TT,minT,ibrav)
    write_alphaT(fout+"/alphaT.dat",TT,alphaT,ibrav)


################################################################################

def fitFvibV(fin,fout,fittypeEtot,fittypeFvib):
    """
    An auxiliary function for handling the qha calculations with Etot+Fvib   
    """

    
    inputfileEtot=fin+"/energy_files/Etot.dat"
    inputfiledos=fin+"/dos_files/output_therm.dat"
    inputfileFvib=fin+"/therm_files/output_therm.dat"
    
    V, E = read_EtotV(inputfileEtot)
    
    a0, cov, chi = fit_Murn(V,E)
   
    print_data(V,E,a0,chi,"E")

    ngeo = len(E)   # total number of geometries
    
    # Read the vibrational energies
    T1, Evib1, Fvib1, Svib1, Cvib1 = read_thermo( inputfileFvib, ngeo )
    nT, T, Evib, Fvib, Svib, Cvib = rearrange_thermo( T1, Evib1, Fvib1, Svib1, Cvib1, ngeo )
   
    # Fit and find the minimun at T>0 K for each T                    
    aTtemp=[]
    chiTtemp=[]
    for i in range(0,nT): 
        print ("####################################################################")
        print ("T=",T[i,0])
  
        a, cov, chi = fit_Murn(V,E+Fvib[i])   
        #print_data(V,E+Fvib[i],a,chi,"E")
             
        print ("Fmin = {:.10e}".format(a[0])+"\tVmin = {:.10e}".format(a[1])
        +"\tB0 = {:.10e}".format(a[2]*RY_KBAR)+"\t chi^2= {:.10e}".format(chi))
        
        aTtemp.append(a)
        chiTtemp.append(chi)
        
    aT = np.array(aTtemp)
    aT.shape = (nT,len(a0))
      
    write_y_T(fout+"/Fmin.dat",T[:,0],aT[:,0],"Fmin")
    write_y_T(fout+"/Vmin.dat",T[:,0],aT[:,1],"Vmin")
    write_y_T(fout+"/B0.dat",T[:,0],aT[:,2]*RY_KBAR,"B0")
    
    # Derived quantities
    betaT = compute_beta(aT[:,1])
    write_y_T(fout+"/beta.dat",T[:,0], betaT,"Beta=1/V dV/dT")
    
    return 
    
    import matplotlib.pyplot as plt
    
    plt.figure(1)
    plt.subplot(311)
    plt.xlabel('T (K)')
    plt.ylabel('Fmin (Ryd) ')
    plt.plot(T[:,0], aT[:,0], 'r')
    plt.subplot(312)
    plt.xlabel('T (K)')
    plt.ylabel('V (a.u.^3) ')
    plt.plot(T[:,0], aT[:,1], 'r')
    plt.subplot(313)
    plt.xlabel('T (K)')
    plt.ylabel('B0 ')
    plt.plot(T[:,0], aT[:,2], 'r')
    plt.show()
    
    plt.figure(2)
    plt.xlabel('T (K)')
    plt.ylabel(r'$\beta$ ($10^6 K^{-1}$)')
    plt.plot(T[:,0], betaT*1E6, 'r')
    plt.show()