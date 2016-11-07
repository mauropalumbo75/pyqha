#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

import numpy as np

def write_elastic_constants(C,S,fname):
    """
    Elastic constants (and elastic compliances) are stored in Voigt notation 
    They are then 6x6 matrices, stored as numpy matrices of shape [6,6]
    So, the elastic constant C11 is in C[0][0], C12 in C[0][1] and so on.
    Same for the elastic compliances   
    """
    
    fout = open(fname, "w")
    # first write the elastic constants
    for i in range(0,6):
        for j in range(0,6):
            fout.write("{:.10e}".format(C[i,j])+"\t")
        fout.write("\n")
        
    fout.write("\n")
    
    # now write the elastic compliances
    for i in range(0,6):
        for j in range(0,6):
            fout.write("{:.10e}".format(S[i,j])+"\t")
        fout.write("\n")  
    fout.close()

################################################################################

# 
def write_C_geo(celldmsx,C,ibrav=4,fCout=""):
    """
    Write elastic constants calculated on a multidimensional grid of lattice parameters
    ngeo defines the total number of geometries evaluated
    Note: the order must be the same as for the total energies in the quasi-harmonic calculations!  
    """
    for i in range(0,6):
        for j in range(0,6):
            C_label = "C" + str(i+1) + str(j+1)
            fname = fCout + C_label + ".dat"
            fout = open(fname, "w")
            if ibrav==4:    # Hexagonal systems
                fout.write("a\t\t\tc/a\t\t\t"+C_label+"\n")
                for k in range(0,len(celldmsx)):
                    fout.write("{:.10e}".format(celldmsx[k,0])+"\t"+"{:.10e}".format(celldmsx[k,2])
                    +"\t"+"{:.10e}".format(C[i][j][k])+"\n")
            fout.close()
                    
################################################################################

# 
def write_CT(Ts,CT,fCout=""):
    """
    Write elastic constants calculated on a multidimensional grid of lattice parameters
    ngeo defines the total number of geometries evaluated
    Note: the order must be the same as for the total energies!
    """
    for i in range(0,6):
        for j in range(0,6):
            C_label = "C" + str(i+1) + str(j+1)
            fname = fCout + C_label + "T.dat"
            fout = open(fname, "w")
            fout.write("T (K)\t\t\t"+C_label+"\n")
            for T in range(0,CT.shape[0]):
                fout.write("{:.10e}".format(Ts[T])+"\t"+"{:.10e}".format(CT[T,i,j])+"\n")
            fout.close()

################################################################################
# 


def write_Etot(celldmsx,Ex,fname,ibrav=4):
    """
    Read cell parameters (a,b,c,alpha,beta,gamma) and energies for a grid of cell
    parameters values from file output_energy1. 
    Each celldms is a vector of lenght 6 containing a,b,c,alpha,beta,gamma respectively
    celldmsx and Ex contains the grid of values of celldms and E so that:
    celldmsx[0] = celldms0      Ex[0] = E0
    celldmsx[1] = celldms1      Ex[1] = E1
    celldmsx[2] = celldms2      Ex[2] = E2
    ........
    values are taken from the file "fname"
    ibrav is the Bravais lattice as in Quantum Espresso and is needed in input (default is cubic)
    """
    fout = open(fname, "w")
    if ibrav==4:    # Hexagonal systems
        fout.write("a\t\t\tc\t\t\tE_tot\n")
        for i in range(0,len(Ex)):
            fout.write("{:.10e}".format(celldmsx[i][0])+"\t"+"{:.10e}".format(celldmsx[i][2])+"\t"+"{:.10e}".format(Ex[i])+"\n")
        fout.close()
        

def write_celldmsT(fname,T,x,ibrav=4):
    fout=open(fname,"w")
    if ibrav==4:    # Hexagonal systems    
        fout.write("# T (K)\t\ta\t\tc\n") 
        for i in range(0,len(T)): 
            fout.write("{:.10e}".format(T[i])+"\t"+"{:.10e}".format(x[i][0])+"\t"+
            "{:.10e}".format(x[i][2])+"\t"+"\n") 
    fout.close()


def write_alphaT(fname,T,alphaT,ibrav=4):
    fout=open(fname,"w")
    if ibrav==4:    # Hexagonal systems    
        fout.write("# T (K) \t\talpha_xx\t\talpha_zz\n") 
        for i in range(0,len(T)): 
            fout.write("{:.10e}".format(T[i])+"\t"
            +"\t"+"{:.10e}".format(alphaT[i,0])+"\t"+"{:.10e}".format(alphaT[i,2])+"\n") 
    fout.close()
  

################################################################################

def write_xy(fname,x,y,labelx,labely):
    """
    This function writes a quantity y versus quantity x into the file fname.
    y and x are arrays and should have the same lenght. labelx and labely
    are the axis labels (possibly with units), written in the header of the file
    (first line).
    """
    fout=open(fname,"w")
    header = "# "+labelx+"\t"+labely+"\n"
    fout.write(header)     
    for i in range(0,len(x)): 
        line = "{:.10e}".format(x[i])+"\t\t"+"{:.10e}".format(y[i])+"\n"
        fout.write(line)
    fout.close()


################################################################################

def write_freq(weights,freq,filename):
    """
    Write frequencies (or Gruneisen parameters) on an extended mesh in a file.
    In this format, q points coordinates are NOT written but the weight of each point yes.
    It can be used to write the Gruneisen mode parameters, giving them in input as freq
    Write the gruneisen parameters 
    """
    nq = freq.shape[0] # total number of q points read
    modes = freq.shape[1]  # number of frequency modes
    fout=open(filename, "w")
    fout.write("         "+str(modes//3)+" 192     192    192   "+str(nq)+"\n")
    fout.write("         F\n")
    for iq in range(0,nq):
        fout.write("              "+"{:.15e}".format(weights[iq])+"\n")
        for ifreq in range(0,modes):
            fout.write("{:.15e}".format(freq[iq,ifreq])+"\n")
    fout.close()
    print ("\nFrequencies (or Gruneisen parameters) written in file "+filename+"\n")
    
    
def write_thermo(fname,T, Evib, Fvib, Svib, Cvib,ZPE,modes):
    fout=open(fname,"w")
    fout.write('# total modes from dos = {:.10e}\n'.format(modes))
    fout.write('# ZPE = {:.10e} Ry/cell\n'.format(ZPE))
    fout.write("# Multiply by 13.6058 to have energies in eV/cell etc..\n\
# Multiply by 13.6058 x 23060.35 = 313 754.5 to have energies in cal/(N mol).\n\
# Multiply by 13.6058 x 96526.0 = 1 313 313 to have energies in J/(N mol).\n\
# N is the number of formula units per cell.\n#\n")
    
    fout.write("# T (K) \tEvib (Ry/cell)\tFvib (Ry/cell)\tSvib (Ry/cell/K)\tCvib (Ry/cell/K)\n") 
    for i in range(0,len(Evib)): 
        fout.write('{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\t{:.10e}\n'.format(T[i],Evib[i],Fvib[i],Svib[i],Cvib[i]))
    fout.close()

################################################################################
# This function may not be necessary, format not used

def write_freq_old(qgeo,freq,filename):
    """
    Write frequencies (or Gruneisen parameters) in a file. In this format also q points
    coordinates are written but not the weight of each point.
    It can be used to write the Gruneisen mode parameters, giving them in input as freq
    """
    nq = qgeo.shape[1] # total number of q points read
    modes = freq.shape[1]  # number of frequency modes
    fout=open(filename, "w")
    fout.write(" &plot nbnd=   "+str(modes)+", nks= "+str(nq)+" /")
    fout.write("\n")
    for iq in range(0,nq):
        fout.write("              ")
        for iqxyz in range(0,3):
            fout.write("{:.6e}".format(qgeo[0][iq][iqxyz])+"    ")    # qpoints coordinates
        fout.write("\n    ")
        for ifreq in range(0,modes):
            fout.write("{:.15e}".format(freq[iq][ifreq])+"\t")  # Gruneisein parameters (one for each mode)
        fout.write("\n")
    fout.close()
    print ("\nFrequencies (or Gruneisen parameters) written in file "+filename+"\n")


################################################################################