#!/usr/bin/env python3
#encoding: UTF-8
#
# This part is still very messy...

import numpy as np
import time
import math
import sys
from read import read_Etot, read_freq, read_freq_ext, read_elastic_constants, \
                 read_qha_elastic_constants, read_freq_ext_geo
from write import write_freq, write_freq_ext, write_alphaT, write_qha_C, write_qha_CT
from constants import RY_KBAR, K_BOLTZMANN_RY, kb1
from fitutils import fit_anis
from minutils import find_min, fquadratic, fquartic
from fitfreqgrun import fitfreq, fitfreqxx, freqmingrun, rearrange_freqx
from fitFvib import fitFvib
from fitC import rearrange_Cx, fitCxx

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
# Compute the volume given the celldms, only for ibrav=4 for now
#
def compute_volume(celldms,ibrav=4):
    if ibrav==4:
        return 0.866025404*celldms[0]*celldms[0]*celldms[2]
        #return 0.866025404*celldms[0]*celldms[0]*celldms[0]*celldms[2]
    

################################################################################
# 
 
def c_qv_python(T,omega):
    """
    This function calculates the mode contribution to the heat capacity at a given T
    and omega. A similar (faster) function should be available as C extension.
    """
    print ("Python c_qv")
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
    x = omega * kb1 / T 
    expx = math.exp(-x)   # exponential term
    x2 = math.pow(x,2)

    return x2*K_BOLTZMANN_RY*expx/math.pow(expx-1.0,2)


################################################################################

def fitS(inputfileEtot, inputpathCx, ibrav, typeSx="quadratic"):
    """
    An auxiliary function for fitting the elastic compliances elements over a
    grid of lattice parameters, i.e. over different geometries.
    """
    # Read the energies (this is necessary to read the celldmsx)
    celldmsx, Ex = read_Etot(inputfileEtot)

    ngeo = len(Ex)
    Cx, Sx = read_qha_elastic_constants(ngeo, inputpathCx)
    
    # This function works for both C and S, here I use it for S 
    Sxx = rearrange_Cx(Sx,ngeo)
    write_qha_C(celldmsx, Sxx, ibrav, inputpathCx)	# Write the S as a function of T for reference
    
    aS, chiS = fitCxx(celldmsx, Sxx, ibrav, True, typeSx)
    
    return aS, chiS


def fitST(aS,mintemp,typeCx):
    """
    An auxiliary function returning the elastic compliances 6x6 tensor at the
    set of lattice parameters given in input as mintemp. These should be the
    lattice parameters at a given temperature obtained from the free energy
    minimization, so that S(T) can be obtained.
    Before calling this function, the polynomial coefficients resulting from 
    fitting the elastic compliances over a grid of lattice parameters, i.e. over
    different geometries, must be obtained and passed as input in aS. 
    """
    S = np.zeros((6,6))
    for i in range(0,6):
       for j in range(0,6):
           if typeCx=="quadratic":
               S[i,j] = fquadratic(mintemp,aS[i,j],ibrav=4)
           elif typeCx=="quartic":
               S[i,j] = fquartic(mintemp,aS[i,j],ibrav=4)  
    return S


################################################################################
# 

def compute_alpha_grun(T,V,S,weights,freq,grun,ibrav=4):    
    """
    This function computes the thermal expansions alpha using the Gruneisein 
    parameters at a given temperature T
    
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
    elif (it[0]==2): # use S at 0 K, calculate V, average frequencies and gruneisen parameters at each T             
        for i in range(it[1],it[2]+1):                
            V=compute_volume(it[11][i],it[9])
            S = fitST(it[15],it[11][i],it[16])
            #S = S * RY_KBAR    # convert elastic compliances in (Ryd/au)^-1
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
    This function simply put together the different temperature ranges where alpha
    was calculated in parallel into a single numpy array
    """
    res = []
    for i in range(0,len(partials)):
        for j in range(0,len(partials[i])):
            res.append(partials[i][j])
            
    return np.array(res)


################################################################################
# 

def compute_alpha_gruneisein(indir,outdir,typeEtot,typeFvib,typeSx,typefreq,\
    ibrav,option=0,rangeT=range(1,10),nproc=1,guess=None):
    """
    This function computes the thermal expansions at different temperatures from
    the gruneisen parameters. Different options are possible as listed below:

    option == 0 -> use V, S, average frequencies and gruneisen parameters at 0 K 
    option == 1 -> use S at 0 K, calculate V, average frequencies and gruneisen parameters at each T 
    option == 2 ->  calculate S, V, average frequencies and gruneisen parameters at each T 

    It also determines how to run the temperature loop in parallel to speed the calculations up
    nproc is the number of processes to be run in parallel
    """
    
    inputfileEtot=indir+"/energy_files/Etot.dat"
    inputfilefreq=indir+"/frequencies"
    
    fittypeEtot="quadratic"
    if (typeEtot==2):
        fittypeEtot="quartic"

    # Read the energies at 0 K
    celldmsx, Ex = read_Etot(inputfileEtot)
    # Fit and find the minimun at 0 K
    a0, chia0 = fit_anis(celldmsx, Ex, ibrav, out=True, type=fittypeEtot)
    if chia0!=None:
        min0, fmin0 = find_min(a0, ibrav, type=fittypeEtot, guess=guess)
    
    # if no guess is given, use the 0 K minimum
    if guess==None:
        guess=min0
    
    if (option!=0): # get the lattice parameters at each T from minimization of F
        T, minT, fminT = fitFvib(indir,outdir,ibrav,typeEtot,typeFvib,guess) 
        rangeT=T
        rangeT=range(1,10)
    else:
        V=compute_volume(min0,ibrav)  # eq. volume at 0 K
        print ("Using the volume at 0 K = ",str(V))
    
    
    if (option==2):   # fit S and (later) get the average values at each T
        # Get the polynomial coefficients aS from fitting the elastic compliances (to be used later to get S(T))
        aS, chiS = fitS(inputfileEtot, inputfileC, ibrav, typeSx)
        aS = aS * RY_KBAR    # convert elastic compliances in (Ryd/au)^-1
    else:           # read S at 0 K from file
        print ("Reading elastic constants and compliances from file "+inputfileC+"...")
        C, S = read_elastic_constants(inputfileC+"output_el_cons.dat")
        print (S)
        S = S * RY_KBAR    # convert elastic compliances in (Ryd/au)^-1
    
    if (option==0):
        # Compute the Gruneisen parameters at 0 K once
        weights, freq, grun = fitfreq(celldmsx, min0, inputfilefreq, ibrav, typefreq="quadratic", compute_grun=True)
        
        # Alternatively, for testing option 0, read the gruneisen parameters from files (already written before)
        #weights, freq = read_freq_ext("average_freq0K")
        #weights, gruntemp1 = read_freq_ext("output_grun_along_a_ext3Dfit")
        #weights, gruntemp2 = read_freq_ext("output_grun_along_c_ext3Dfit")
        #nq = gruntemp1.shape[0]
        #modes = gruntemp1.shape[1]
        #grun = np.zeros((6,nq,modes))
        #grun[0] = gruntemp1
        #grun[1] = gruntemp1
        #grun[2] = gruntemp2
    else:          
        # get the weigths and the frequencies from files 
        weightsx, freqx = read_freq_ext_geo(inputfilefreq,range(1,celldmsx.shape[0]+1))
        weights = weightsx[0,:]

        print ("Rearranging frequencies...")
        freqxx = rearrange_freqx(freqx)
        print ("Done!")
        del freqx       # try to free some memory, but not ideal with the garbage collector
    
        print ("Fitting frequencies...")
        afreq, chifreq = fitfreqxx(celldmsx, freqxx, ibrav, True, typefreq)
        print ("Done!")  
     
    ############################################################################
    # Beginning of parallel part

    n_T_in_total = len(rangeT)
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
    if (option==0):
        for i in range(0,nproc):
            it.append([option,startT[i],endT[i],rangeT,V,S,weights,freq,grun,ibrav])
    elif (option==1):
        for i in range(0,nproc):
            it.append([option,startT[i],endT[i],rangeT,None,S,weights,None,None,ibrav,afreq, minT, freqxx.shape[0],freqxx.shape[1], typefreq])
    elif (option==2):
        for i in range(0,nproc):
            it.append([option,startT[i],endT[i],rangeT,None,None,weights,None,None,ibrav,afreq, minT, freqxx.shape[0],freqxx.shape[1], typefreq, aS, typeSx])     
    else:
        print ("Option not implemented. Exiting...")
        exit()
        
    # Call pool.map to distribute the work in parallel. Results are collected
    # in alphaTpartials which is a list with the results in different 
    # temperature ranges
    alphaTpartials = pool.map(compute_alpha_gruneisen_loopparallel,it)

    # Join all partial results in the numpy array alphaT
    alphaT = join(alphaTpartials)
    
    # End of parallel part
    ############################################################################
 
    fout = outdir+"/alpha_gruneisen"
    write_alphaT(fout,rangeT,alphaT,ibrav)
         
        
################################################################################
#   MAIN, for testing
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
    nproc=2     # number of defauls processes
    
    # Read and process command line parameters, if any
    print ("Possible options are: -fileEtot -fileFvib -filefreq -fileC "+ 
    "-ibrav -fittypeEtot -fittypeFvib -fittypefreq -fittypeSx -option -nproc -guess\n")
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
            elif sys.argv[i] == "-nproc": nproc = int(sys.argv[i+1])
            elif sys.argv[i] == "-guess": guess = sys.argv[i+1]

    print (nproc)
    exit()
    compute_alpha_gruneisein(inputfileEtot,inputfileFvib,inputfileC,inputfilefreq,
    typeEtot,typeFvib,typeSx,typefreq,ibrav,guess,option,range(1,10),nproc)
      
    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")

