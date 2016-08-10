# -*- coding: utf-8 -*-
__author__ = "mauro"
__date__ = "$30-Jan-2016 21:25:36$"

import numpy as np

################################################################################

def read_elastic_constants(fname):
    """
    This function reads and returns the elastic constants and compliances.
    Elastic constants (and elastic compliances) are stored in Voigt notation 
    They are then 6x6 matrices, stored as numpy matrices of shape [6,6]
    So, the elastic constant C11 is in C[0][0], C12 in C[0][1] and so on.
    Same for the elastic compliances
    """
    C=[]          # C, S matrices for elastic constants and elastic compliances
    S=[]
    c=[]          # c, s vectors of partial elastic constants and elastic compliances
    s=[]
    i=0           # initialize counters
    countline=1
    with open(fname, "r") as lines:
        for line in lines:
            linesplit=line.split()
            if countline<13:                  # Read the first 12 lines with the elastic constants
                if (i%2)==0:
                    for j in range(0,4):
                        c.append(float(linesplit[j]))
                    i += 1
                else:
                    for j in range(0,2):
                        c.append(float(linesplit[j])) 
                    i=0
                    C.append(c)
                    c=[]
            elif (countline>13):             # Read the last 12 lines with the elastic compliances (1 line in between)
                if (i%2)==0:
                    for j in range(0,4):
                        s.append(float(linesplit[j]))
                    i += 1
                else:
                    for j in range(0,2):
                        s.append(float(linesplit[j])) 
                    i=0
                    S.append(s)
                    s=[]
            countline += 1
   
    return np.array(C),np.array(S)


################################################################################

def read_qha_elastic_constants(ngeo, path=""):
    """
    Read elastic constants calculated on a multidimensional grid of lattice parameters
    ngeo defines the total number of geometries evaluated
    Note: the order must be the same as for the total energies!
    """
    Cx = []
    Sx = []
    for i in range(1,ngeo+1):
        fname = path + "output_el_cons.g" + str(i) + ".dat"
        C, S = read_elastic_constants(fname)
        Cx.append(C)
        Sx.append(S)
        
    return np.array(Cx), np.array(Sx)

################################################################################

def read_EtotV(fname):
    """
    Read volumes and energies for a isotropic calculation. Values are taken from the file "fname"
    ibrav is as in Quantum Espresso and is needed in input (default is cubic)
    """
    Vx=[]
    Ex=[]
    
    with open(fname, "r") as lines:
        for line in lines:
            linesplit=line.split()
            V=float(linesplit[0])
            E=float(linesplit[1])
            Vx.append(V)
            Ex.append(E)
    
    return np.array(Vx), np.array(Ex)


################################################################################

def read_Etot(fname, ibrav=4):
    """
    Read cell parameters and energies from input file fname. 
    Each celldms is a vector of lenght 6 containing a,b,c,alpha,beta,gamma respectively
    celldmsx and Ex contains the grid of values of celldms and E so that:
    celldmsx[0] = celldms0      Ex[0] = E0
    celldmsx[1] = celldms1      Ex[1] = E1
    celldmsx[2] = celldms2      Ex[2] = E2
    ........
    values are taken from the file "fname"
    ibrav is as in Quantum Espresso and is needed in input (default is cubic)
 
    Values of lattice parameters from QE are converted in standard lattice parameters,
    i.e. b and c instead of factors b/a and c/a
    """
    celldmsx=[]
    Ex=[]
    
    if ibrav==4:    # Hexagonal systems
        with open(fname, "r") as lines:
            for line in lines:
                celldms=np.zeros([6])
                linesplit=line.split()
                celldms[0]=float(linesplit[0])
                #celldms[2]=float(linesplit[1]) # as it is
                celldms[2]=float(linesplit[1])*celldms[0]   # convert from c/a to c
                E=float(linesplit[2])
                celldmsx.append(celldms)
                Ex.append(E)
    
    return np.array(celldmsx), np.array(Ex)

################################################################################
# 
 
def read_thermo(fname, ngeo=1):
    """
    Read cell vibrational thermodynamic functions (Evib, Fvib, Svib, Cvib) as a 
    function of temperature (nT temperatures as in input files) and for different
    geometries (ngeo) 
    Input file(s) must have the following format:
    T   Evib    Fvib    Svib    Cvib    
    1   ...     ...     ....    ...
    ........
    values are taken from the file(s) fname.g1, fname.g2, etc. for each geometry.

    Returning values are nT*ngeo numpy matrices (T,Evib,Fvib,Svib,Cvib) containing the 
    temperatures and the above mentioned thermodynamic functions as for example:
    Fvib[T][geo] -> Fvib at the temperature "T" for the geometry "geo"
    """
    T=[]
    Evib=[]
    Fvib=[]
    Svib=[]
    Cvib=[]
    for i in range(0,ngeo):
        gT=[]
        gEvib=[]
        gFvib=[]
        gSvib=[]
        gCvib=[]
        fname2=fname+".g"+str(i+1)
        with open(fname2, "r") as lines:
            for line in lines:
                linesplit=line.split()
                if line!="" and linesplit[0]!="#":
                    gT.append( float(linesplit[0]) )
                    gEvib.append( float(linesplit[1]) )
                    gFvib.append( float(linesplit[2]) )
                    gSvib.append( float(linesplit[3]) )
                    gCvib.append( float(linesplit[4]) )
        T.append(gT)
        Evib.append(gEvib)
        Fvib.append(gFvib)
        Svib.append(gSvib)
        Cvib.append(gCvib)
    
    return np.array(T), np.array(Evib), np.array(Fvib), np.array(Svib), np.array(Cvib)


################################################################################

def read_freq(filename):
    """
    Read the phonon frequencies at each q point from a frequency file.  
    Input file (for example from thermo_pw) has the following format:
    &plot nbnd=     6, nks=   342 /
            0.000000  0.000000  0.000000
    -0.0000   -0.0000    0.0000  200.1106  200.1106  299.8824
            0.022222  0.000000  0.000000
    9.8401   11.1499   17.7609  200.2470  200.3982  299.7714
            0.044444  0.000000  0.000000
    ....

    Returning values are a nq*3 matrix q, each q[i] being a q point (vector of 3 elements)
    and a nq*modes matrix freq, each element freq[i] being the phonon frequencies
    (vector of modes elements)
    """
    q=[]
    freq=[]
    count=0
    print ("Reading file of frequencies "+filename+"...")
    with open(filename) as lines:
        for line in lines:
            linesplit=line.split()
            if linesplit[0]=="&plot" and count==0:   # read first line with modes and total number of q points
                modes=int(linesplit[2].strip(","))
                nq=int(linesplit[4])
                count = count + 1
            elif count==1:  # it's a q point line
                q.append([float(linesplit[0]),float(linesplit[1]),float(linesplit[2])])
                count = count + 1
            elif count==2:  # it's a frequency line
                frequencies=[]
                for i in range(0,modes):
                    frequencies.append(float(linesplit[i]))   # list of frequencies in this line (for a given q point), 3*modes
                freq.append(frequencies)
                count = 1
    return np.array(q), np.array(freq)


################################################################################
 

def read_freq_ext(filename):
    """
    Read the phonon frequencies at each q point from a frequency file. The format 
    of this file is different from the one read by the function read_freq and
    contains usually more frequencies, each with a weight, but no qpoint coordinates.
    Input file (for example from thermo_pw) has the following format:
         2       192       192       192    307393
         F
         0.141285083912966E-06
    -0.7142368571499E-05
    -0.5893225492475E-05
    0.1387907734881E-05
    0.2001106263587E+03
    0.2001106263588E+03
    0.2998824026349E+03
         0.282570167825932E-06
    0.1584937719099E+01
    0.1584937719102E+01
    0.3082066539423E+01
    0.2001045986544E+03
    0.2001045986545E+03
    0.2998676623586E+03
    ....
    First line contains n. atoms, nqx, nqy, nqz, nq total.
    Second line not read.
    Third line: weight of the first qpoint
    Following lines: phonon frequencies (their number is modes=3*n. atoms), one per line
    then again: weight of the next qpoint, phonon frequencies (3*modes), one per line, etc.

    Weights are diffent because of simmetry

    Returning values are a nq vector weights, each weights[i] being the weight of a q point 
    and a nq*modes matrix freq, each element freq[i] being the phonon frequencies
    (vector of modes elements) at the qpoint i
    """
    countline=0 
    weights=[]
    freq=[]
    frequencies=[]
    print ("Reading file of frequencies "+filename+"...")
    with open(filename, "r") as lines:
        for line in lines:
            linesplit=line.split()
            if countline==0:                            # this is the very first line with natoms and nq
                modes=3*int(linesplit[0].strip())
                nq=int(linesplit[4])
                countline = countline + 1
            elif countline==1:  # it's the second line, just skip it
                countline = countline + 1
            elif countline==2:  # it's a weight line
                weights.append(float(linesplit[0]))
                countline = countline + 1
            elif countline>2:  # it's a frequency line
                frequencies.append( float(linesplit[0]) )   # collect the frequencies from each line
                countline = countline + 1
                if (countline==(modes+3)): 		     # this is the last line with a frequency, next line is again a weight
                    freq.append(frequencies)
                    frequencies=[]
                    countline = 2

    return np.array(weights), np.array(freq)


################################################################################ 

def read_freq_geo(inputfilefreq,rangegeo):
    """
    Read the frequencies for all geometries where the gruneisen parameters must be
    calculated. Start, stop, step must be given accordingly. It can be used to read
    the frequencies only at some geometries from a larger set, if necessary, 
    providing the proper start, stop and step values.

    Notes:
    nq = qgeo.shape[1] -> total number of q points read
    modes = freqgeo.shape[2] -> number of frequency modes
    """
    qgeo=[]
    freqgeo=[]
    for i in rangegeo:
        q, freq = read_freq(inputfilefreq+str(i))
        qgeo.append(q)
        freqgeo.append(freq)  
        
    return np.array(qgeo), np.array(freqgeo)


################################################################################


def read_freq_ext_geo(inputfilefreq,rangegeo):
    """
    Read the frequencies for all geometries where the gruneisen parameters must be
    calculated. Start, stop, step must be given accordingly. It can be used to read
    the frequencies only at some geometries from a larger set, if necessary, 
    providing the proper start, stop and step values.

    Notes:
    nq = qgeo.shape[1] -> total number of q points read
    modes = freqgeo.shape[2] -> number of frequency modes
    """
    weightsgeo=[]
    freqgeo=[]
    for i in rangegeo:
        weights, freq = read_freq_ext(inputfilefreq+str(i))
        weightsgeo.append(weights)
        freqgeo.append(freq)  
        
    return np.array(weightsgeo), np.array(freqgeo)


################################################################################
################################################################################
# Work in progress from here... These routines are not ready yet to be used.


def read_celldmt_hex(filename):
  T=[]
  celldmst1=[]
  celldmst3=[]
  with open(filename, "r") as lines:
    for line in lines:
        linesplit=line.split()
        if (linesplit[0]!="#"):
          T.append(float(linesplit[0]))
          celldmst1.append(float(linesplit[1]))
          celldmst3.append(float(linesplit[2]))
  return T, celldmst1, celldmst3


def read_alpha(fname):
  T=[]
  a=[]
  ca=[]
  alphaxx=[]
  alphazz=[]
  with open(fname, "r") as lines:
    for line in lines:
      linesplit=line.split()
      if line!="" and linesplit[0]!="#":
        T.append( float(linesplit[0]) )
        a.append( float(linesplit[1]) )
        ca.append( float(linesplit[2]) )
        alphaxx.append( float(linesplit[4]) )
        alphazz.append( float(linesplit[5]) )
  x = [a,ca]
  alpha = [alphaxx, alphazz]
  return T,x,alpha
