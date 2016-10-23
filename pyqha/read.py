#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

import numpy as np

################################################################################

def read_elastic_constants(fname):
    """
    This function reads and returns the elastic constants and compliances from
    the file *fname* .
    Elastic constants (and elastic compliances) are stored in Voigt notation 
    They are then 6x6 matrices, stored as numpy matrices of shape [6,6]
    So, the elastic constant C11 is in C[0,0], C12 in C[0,1] and so on.
    Same for the elastic compliances.
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

def read_elastic_constants_geo(fC,ngeo):
    """
    Read elastic constants calculated on a multidimensional grid of lattice parameters
    *ngeo* defines the total number of geometries evaluated
    Note: the order must be the same as for the total energies!
    """
    Cx = []
    Sx = []
    for i in range(1,ngeo+1):
        fname = fC + str(i)
        C, S = read_elastic_constants(fname)
        Cx.append(C)
        Sx.append(S)
        
    return np.array(Cx), np.array(Sx)

################################################################################

def read_EtotV(fname):
    """
    Read cell volumes and the corresponding energies from input file *fname*
    (1st col, volumes, 2nd col energies). Units must be :math:`a.u.^3` and 
    :math:`Ryd/cell`
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

def read_Etot(fname, ibrav=4, bc_as_a_ratio = True):
    """
    Read cell parameters *(a,b,c)* and the corresponding energies from input file *fname*. 
    Each set of cell parameters is stored in a numpy array of lenght 6 
    for *(a,b,c,alpha,beta,gamma)* respectively. This is done for a future possible 
    extension but for now only the first 3 elements are used (the others are always 0).
    All sets are stored in *celldmsx* and *Ex*, the former is a nE*6 matrix, 
    the latter is a nE array.
    
    *ibrav* identifies the Bravais lattice as in Quantum Espresso and is needed 
    in input (default is 4, i.e. hexagonal cell). The input file format depends
    on *ibrav*, for example in the hex case, the first two columns are for *a* and
    *c* and the third is for the energies.
    
    If *bc_as_a_ratio=True*, the input data are assumed to be given as 
    :math:`(a,b/a,c/a)` in the input file and hence converted into :math:`(a,b,c)`
    which is how they are always stored internally in :py:mod:`pyqha`.
    
    Units must be :math:`a.u.` and :math:`Ryd/cell`
 
    """
    celldmsx=[]
    Ex=[]
    
    if ibrav==4:    # Hexagonal systems
        with open(fname, "r") as lines:
            for line in lines:
                celldms=np.zeros([6])
                linesplit=line.split()
                celldms[0]=float(linesplit[0])
                if (bc_as_a_ratio==True):
                    celldms[2]=float(linesplit[1])*celldms[0]   # convert from c/a to c
                else:
                    celldms[2]=float(linesplit[1])              # as it is
                E=float(linesplit[2])
                celldmsx.append(celldms)
                Ex.append(E)
    else:
        return
    
    return np.array(celldmsx), np.array(Ex)

################################################################################
# 
 
def read_thermo(fname, ngeo=1):
    """
    Read vibrational thermodynamic functions (Evib, Fvib, Svib, Cvib) as a 
    function of temperature from the input file *fname*. *ngeo* is the number
    of input files to read, corresponding for example to different geometries
    in a quasi-harmonic calculation.
    If *ngeo>1* reads from the files *fname1*, *fname2*, etc. up to *ngeo*  
    Input file(s) have the following format:

    +------+-----------------+-----------------+-----------------+-----------------+
    | T    | :math:`E_{vib}` | :math:`F_{vib}` | :math:`S_{vib}` | :math:`C_{vib}` | 
    +======+=================+=================+=================+=================+
    | 1    | ...             | ...             | ...             | ...             |
    +------+-----------------+-----------------+-----------------+-----------------+  
    
    Lines starting with "#" are not read (comments).    

    Returning values are :math:`nT*ngeo` numpy matrices (T,Evib,Fvib,Svib,Cvib) containing the 
    temperatures and the above mentioned thermodynamic functions as for example:
    Fvib[T,geo] -> Fvib at the temperature *T* for the geometry *geo*
    
    Units must be *K* for temperature, *Ryd/cell* for energies, *Ryd/cell/K* for
    entropy and heat capacity. 
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
        if (ngeo>1):
            fname2=fname+str(i+1)
        else:
            fname2=fname
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
    This funcstion reads the phonon frequencies at each *q* point from a frequency file.  
    Input file has the following format (to be done).

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
    Input file has the following format:

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
    calculated. 

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


def read_dos(filename):
    """
    Read the phonon density of states (y axis) and the corresponding energies (x axis)
    from the input file *filename* (1st col energies, 2nd col DOS) and store it
    in two numpy arrays which are returned. 

    """
    E = []
    dens = []
    print ("Reading phonon dos file "+filename+"...")
    with open(filename, "r") as lines:
        for line in lines:
            linesplit=line.split()
            if (float(linesplit[0])<1.0E-10):
                E.append(1.0E-10)
            else:
                E.append(float(linesplit[0]))
            dens.append(float(linesplit[1]))

    return np.array(E), np.array(dens)


def read_dos_geo(fin,ngeo):
    """
    Read the phonon density of states and energies as in :py:func:`read_dos` from *ngeo* input files
    *fin1*, *fin2*, etc. and store it in two numpy matrices which are returned. 
    """
    
    # read the first file to determine the number of DOS points
    E, dos = read_dos(fin+"1") 
     
    gE = np.zeros((len(E),ngeo))
    gdos = np.zeros((len(E),ngeo))   
    gE[:,0] = E
    gdos[:,0] = dos
    
    # then read the rest of files 
    for i in range(1,ngeo):
        EE, ddos = read_dos(fin+str(i+1)) 
        gE[:,i] = EE
        gdos[:,i] = ddos
    
    return gE, gdos

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


################################################################################
#   MAIN, for testing
################################################################################
#

if __name__ == "__main__":
    pass
