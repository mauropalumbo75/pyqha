#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

# This part is still very messy...

import numpy as np
import time
import math
import sys
from read import read_Etot, read_freq, read_freq_ext, read_elastic_constants, \
                 read_elastic_constants_geo, read_freq_ext_geo
from write import write_freq, write_freq_ext, write_alphaT, write_qha_C, write_qha_CT
from constants import RY_KBAR, K_BOLTZMANN_RY, kb1
from fitutils import fit_anis
from minutils import find_min, fquadratic, fquartic
from fitfreqgrun import fitfreq, fitfreqxx, freqmingrun, rearrange_freqx
from fitFvib import fitFvib
from fitC import rearrange_Cx, fitCxx

from grunc import c_qvc  # This is the same routine c_qv implemented in C to speed it up


################################################################################
# 
# Compute the volume given the celldms, only for ibrav=4 for now
def compute_volume(celldms,ibrav=4):
    if ibrav==4:
        return 0.866025404*celldms[0]*celldms[0]*celldms[2]
        #return 0.866025404*celldms[0]*celldms[0]*celldms[0]*celldms[2]
    

################################################################################
# 
# Function to calculate the mode contribution to the heat capacity at a given T
# and omega
# This is a possible bottleneck as it is implemented in Python. It would be 
# better to write it in C and link it to CPython or similar
#
# 
def c_qv(T,omega):
    if (T<1E-9 or omega<1E-9):
        return 0.0
    x = omega * kb1 / T 
    expx = math.exp(-x)   # exponential term
    x2 = math.pow(x,2)
    if expx>1E-3:           # compute normally
        return x2*K_BOLTZMANN_RY*expx/math.pow(expx-1.0,2)
    else:                   # Taylor series
        return K_BOLTZMANN_RY*expx* (x/math.pow(x-0.5*math.pow(x,2)+ 
        0.16666666666666667*math.pow(x,3)+0.04166666666666666667*math.pow(x,4),2))

# Same as c_qv but no if. Slightly more efficient, roughly a 30% faster
def c_qv2(T,omega):
    x = omega * kb1 / T 
    expx = math.exp(-x)   # exponential term
    x2 = math.pow(x,2)

    return x2*K_BOLTZMANN_RY*expx/math.pow(expx-1.0,2)


################################################################################
# 
# This function computes the thermal expansions alpha using the Gruneisein 
# parameters
# more comments to be added
# First with min0, freq and grun T-independent
#
# More ibrav types to be implemented
def compute_alpha_grun(T,V,S,weights,freq,grun,ibrav=4):    
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
                temp = c_qvc(T,freq[iq,mode]) * weights[iq]   # should be quicker with this additional variable
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


def compute_alpha_gruneisein(inputfileEtot,inputfileC,inputfilefreq,rangeT,typeEtot,typefreq,ibrav):
    # Read the energies 
    celldmsx, Ex = read_Etot(inputfileEtot)
    # Fit and find the minimun at 0 K
    a0, chia0 = fit_anis(celldmsx, Ex, ibrav, out=True, type=typeEtot)
    if chia0!=None:
        min0, fmin0 = find_min(a0, ibrav, type=typeEtot, guess=guess)
    
    # First read the elastic compliances which are need for the thermal expansions
    print ("Reading elastic constants and compliances from file "+inputfileC+"...")
    C, S = read_elastic_constants(inputfileC)
    #print (S)
    
    # Compute the Gruneisen parameters 
    weights, freq, grun = fitfreq(celldmsx, min0, inputfilefreq, ibrav, typefreq="quadratic", compute_grun=True)
    
    # Alternatively, we can read the gruneisen parameters from files (already written before)
    #weights, freq = read_freq_ext("average_freq0K")
    #weights, gruntemp1 = read_freq_ext("output_grun_along_a_ext3Dfit1.0")
    #weights, gruntemp2 = read_freq_ext("output_grun_along_c_ext3Dfit1.0")
    #nq = gruntemp1.shape[0]
    #modes = gruntemp1.shape[1]
    #grun = np.zeros((6,nq,modes))
    #grun[0] = gruntemp1
    #grun[1] = gruntemp1
    #grun[2] = gruntemp2
    
    V=compute_volume(min0,ibrav)  # eq. volume at 0 K
    print ("V = ",str(V))
    S = S * RY_KBAR    # convert elastic compliances in (Ryd/au)^-1

    alphaT= np.zeros((len(rangeT),6))
    counterT=0
    for T in rangeT:
        alpha = compute_alpha_grun(T,V,S,weights,freq,grun)
        alphaT[counterT]=alpha
        counterT += 1
        
        print ("T= "+str(T)+"\t"+str(alpha[0])+"\t"+str(alpha[2]))

    write_alphaT("alpha_gruneisen",rangeT,alphaT,4)
    



def compute_alpha_gruneiseinT(inputfileEtot,inputfileFvib,inputfileC,inputfilefreq,typeEtot,typeFvib,typefreq,ibrav,guess):
    # Read the energies 
    celldmsx, Ex = read_Etot(inputfileEtot)
    
    T, minT, fminT = fitFvib(inputfileEtot,inputfileFvib,ibrav,typeEtot,typeFvib,guess)    

    # First read the elastic compliances which are need for the thermal expansions
    print ("Reading elastic constants and compliances from file "+inputfileC+"...")
    C, S = read_elastic_constants(inputfileC)
    print (S)
    S = S * RY_KBAR    # convert elastic compliances in (Ryd/au)^-1

    # get the weigths and the frequencies from files 
    weightsx, freqx = read_freq_ext_geo(inputfilefreq,range(1,celldmsx.shape[0]+1))
    weights = weightsx[0,:]

    print ("Rearranging frequencies...")
    freqxx = rearrange_freqx(freqx)
    print ("Done!")
    del freqx
    
    print ("Fitting frequencies...")
    afreq, chifreq = fitfreqxx(celldmsx, freqxx, ibrav, True, typefreq)
    print ("Done!")  
    
    alphaT= np.zeros((len(T),6))
    for i in range(0,len(T)):
        # Compute the Gruneisen parameters, the average frequencies and alpha at each T
        V=compute_volume(minT[i],ibrav)
        print ("V = ",str(V))
        freq, grun = freqmingrun(afreq, minT[i], freqxx.shape[0],freqxx.shape[1], ibrav, typefreq)
    
        #write_freq_ext(weights,freq,"average_freqPython"+str(T[i]))
        #write_freq_ext(weights,grun[0],"output_grun_along_a_ext3Dfit"+str(T[i]))
        #write_freq_ext(weights,grun[2],"output_grun_along_c_ext3Dfit"+str(T[i]))

        alpha = compute_alpha_grun(T[i],V,S,weights,freq,grun)
        print ("T= "+str(T[i]))
        print (alpha)
        alphaT[i,:] = alpha
        
    write_alphaT("alpha_gruneisenT",T,alphaT,4)

################################################################################
# 
# This function is only meant to test the Cqv modes. It has to be removed later...
#
def testCqv(inputfilefreq, rangeT, out="Cqvtest"):
    weights, freq = read_freq_ext(inputfilefreq)
    nq = freq.shape[0] # total number of q points read
    modes = freq.shape[1]  # number of frequency modes
    
    for T in rangeT:
        Cqv = []
        for iq in range(0,nq):
            Cqvq=[]
            for ifreq in range(0,modes):
                temp = c_qv2(T,freq[iq,ifreq])
                Cqvq.append(temp)
            Cqv.append(Cqvq)
    
        Cqv = np.array(Cqv)
        outT = out+str(T)
        write_freq_ext(weights,Cqv,outT)


################################################################################
# An auxiliary function for fitting the elastic constant elements of Sxx 
#
#
def fitS(inputfileEtot, inputpathCx, ibrav, typeSx="quadratic"):
    # Read the energies (this is necessary to read the celldmsx)
    celldmsx, Ex = read_Etot(inputfileEtot)

    ngeo = len(Ex)
    Cx, Sx = read_elastic_constants_geo(ngeo, inputpathCx)
    
    # This function works for both C and S, here I use it for S 
    Sxx = rearrange_Cx(Sx,ngeo)
    write_qha_C(celldmsx, Sxx, ibrav, inputpathCx)	# Write the S as a function of T for reference
    
    aS, chiS = fitCxx(celldmsx, Sxx, ibrav, True, typeSx)
    
    return aS, chiS

def fitST(aS,mintemp,typeCx):
    S = np.zeros((6,6))
    for i in range(0,6):
       for j in range(0,6):
           if typeCx=="quadratic":
               S[i,j] = fquadratic(mintemp,aS[i,j],ibrav=4)
           elif typeCx=="quartic":
               S[i,j] = fquartic(mintemp,aS[i,j],ibrav=4)  
    return S
   


def compute_alpha_gruneiseinCT(inputfileEtot,inputfileFvib,inputpathCx,inputfilefreq,typeEtot,typeFvib,typeSx,typefreq,ibrav,guess):
    # Read the energies 
    celldmsx, Ex = read_Etot(inputfileEtot)
    
    T, minT, fminT = fitFvib(inputfileEtot,inputfileFvib,ibrav,typeEtot,typeFvib,guess)    

    # Get the polynomial coefficients aS from fitting the elastic compliances (to be used later to get S(T))
    aS, chiS = fitS(inputfileEtot, inputpathCx, ibrav, typeSx)

    # Now get the polynomial coeffients afreq from fitting the frequencies (to be used later to get average frequencies and 
    # gruneisen parameters as a function of T)
    weightsx, freqx = read_freq_ext_geo(inputfilefreq,range(1,celldmsx.shape[0]+1))
    weights = weightsx[0,:]

    print ("Rearranging frequencies...")
    freqxx = rearrange_freqx(freqx)
    print ("Done!")
    del freqx
    
    print ("Fitting frequencies...")
    afreq, chifreq = fitfreqxx(celldmsx, freqxx, ibrav, True, typefreq)
    print ("Done!")  
    
    alphaT= np.zeros((len(T),6))
    for i in range(0,len(T)):
        # Compute the Gruneisen parameters, the average frequencies and alpha at each T
        V=compute_volume(minT[i],ibrav)
        print ("V = ",str(V))

        S = fitST(aS,minT[i],typeSx)
        print (S)
        S = S * RY_KBAR    # convert elastic compliances in (Ryd/au)^-1

        freq, grun = freqmingrun(afreq, minT[i], freqxx.shape[0],freqxx.shape[1], ibrav, typefreq)
    
        #write_freq_ext(weights,freq,"average_freqPython"+str(T[i]))
        #write_freq_ext(weights,grun[0],"output_grun_along_a_ext3Dfit"+str(T[i]))
        #write_freq_ext(weights,grun[2],"output_grun_along_c_ext3Dfit"+str(T[i]))

        alpha = compute_alpha_grun(T[i],V,S,weights,freq,grun)
        print ("T= "+str(T[i]))
        print (alpha)
        alphaT[i,:] = alpha
        
    write_alphaT("alpha_gruneisenT",T,alphaT,4)
        
        
################################################################################
#   MAIN
################################################################################
#

if __name__ == "__main__":
    
    start_time = time.time()
  
    # Default command line parameters
    #inputfileEtot = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    #inputfileFvib = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/therm_files/output_therm.dat"  
    #inputfilefreq = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/save_frequencies.dat.g"
    #inputfileC = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/output_el_cons.dat"
    #inputfileC = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/" 
    
    #inputfileEtot = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    #inputfileFvib = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/therm_files/output_therm.dat"
    #inputfilefreq = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/save_frequencies.dat.g"
    #inputfileC = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/"
    
    inputfileEtot = "Tc.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    inputfileFvib = "Tc.pz-spn-kjpaw_psl.1.0.0.UPF/therm_files/output_therm.dat"
    inputfilefreq = "Tc.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/save_frequencies.dat.g"
    inputfileC = "Tc.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/"
    
    ibrav = 4  # default value, to be later changed to 1
    typeEtot = "quartic"
    typeFvib = "quartic"
    typefreq = "quadratic"
    typeSx = "quartic"
    # default guess for Os pz
    #guess=[5.12376396,0.0,8.19355157,0.0,0.0,0.0]
    # default guess for Re pz
    #guess=[5.16916792,0.0,8.33684861,0.0,0.0,0.0]
    # default guess for Tc pz
    guess=[5.12374914,0.0,8.19314311,0.0,0.0,0.0]
    
    option=0    # possible choices: 0, 1, 2   
    
    # Read and process command line parameters, if any
    print ("Possible options are: -fileEtot -fileFvib -filefreq -fileC "+ 
    "-ibrav -fittypeEtot -fittypeFvib -fittypefreq -fittypeSx -option -guess\n")
    print ("Assuming default values for unspecified options: \n")
    if len(sys.argv)>1:
        for i in range(1,len(sys.argv)):
            if sys.argv[i] == "-fileEtot": inputfileEtot = sys.argv[i+1]
            elif sys.argv[i] == "-fileFvib": inputfileFvib = sys.argv[i+1]
            elif sys.argv[i] == "-filefreq": inputfileFvib = sys.argv[i+1]
            elif sys.argv[i] == "-fileC": inputfileFvib = sys.argv[i+1]
            elif sys.argv[i] == "-ibrav": ibrav = int(sys.argv[i+1])
            elif sys.argv[i] == "-fittypeEtot": typeEtot = sys.argv[i+1]
            elif sys.argv[i] == "-fittypeFvib": typeFvib = sys.argv[i+1]
            elif sys.argv[i] == "-fittypefreq": typefreq = sys.argv[i+1]
            elif sys.argv[i] == "-fittypeSx": typeSx = sys.argv[i+1]
            elif sys.argv[i] == "-option": option = int(sys.argv[i+1])
            elif sys.argv[i] == "-guess": guess = sys.argv[i+1]

    # Add ibrav parameter in all 3
    compute_alpha_gruneisein(inputfileEtot,inputfileC,inputfilefreq,range(1,2000,1),typeEtot,typefreq,ibrav) 
    #compute_alpha_gruneiseinT(inputfileEtot,inputfileFvib,inputfileC,inputfilefreq,typeEtot,typeFvib,typefreq,ibrav,guess)
    #compute_alpha_gruneiseinCT(inputfileEtot,inputfileFvib,inputpathCx,inputfilefreq,typeEtot,typeFvib,typeSx,typefreq,ibrav,guess)

    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")

