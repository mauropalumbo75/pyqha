#!/usr/bin/env python3
#encoding: UTF-8

################################################################################
# 
# This file contains some functions to compute the Gruneisein parameters
# using 1-dimensional polynomials along each necessary direction according to the
# unit cell of the system. For example, for hexagonal systems it is along a and c
# directions where a and c are the vectors of the hex cell.
# 
# The results obtained in this way are less accurate than performing the fit
# using a multi-dimensional polynomial over the whole grid of lattice parameters.
# The latter is thus preferable (implemented in fitfreqgrun.py).
# Note: the functions in this file would need to be tested and possibly improved

import numpy as np
import time
from read import read_Etot
from read import read_freq
from read import read_freq_ext
from write import write_freq
from write import write_freq_ext
from read import read_freq_ext_geo


################################################################################
# 
# Find the center geometries. Remember indexex in lists starts from 0...
def find_geocenters(ngeo):
    cgeo = np.zeros(6,dtype=np.int)
    for i in range(0,len(ngeo)):
        if (ngeo[i] % 2)==0:
            cgeo[i] = ngeo[i]//2-1
        else:
            cgeo[i] = ngeo[i]//2-1
            
    return cgeo+1


################################################################################
# 
# Compute the Gruneisen parameters along one direction.
# This function uses a 1-dimensional polynomial of fourth degree to fit the 
# frequencies along a certain direction (along a and c axis in hexagonal systems
# for example). 
#
def compute_grun_along_one_direction(nq,modes,ngeo,cgeo,celldmsx,freqgeo,rangegeo,xindex=0):
        # set a numpy array of volumes for the fit (n=5)
        xtemp=[]
        for igeo in rangegeo:
            xtemp.append(celldmsx[igeo,xindex])
        x=np.array(xtemp)

        grun=[]
        for iq in range(0,nq):
            grunq=[]
            for ifreq in range(0,modes):
                ytemp=[]
                for igeo in rangegeo:
                    ytemp.append(freqgeo[igeo,iq,ifreq])
                y=np.array(ytemp)    
                z=np.polyfit(x, y, 4)
                p=np.poly1d(z)
                pderiv=np.polyder(p)
                if freqgeo[cgeo[xindex],iq,ifreq]<1E-3:
                    grunq.append(0.0)
                else:
                    grunq.append(pderiv(celldmsx[cgeo[xindex],xindex])/freqgeo[cgeo[xindex],iq,ifreq]) #*celldmsx[cgeo[xindex],xindex])
                #print (x,y,z,p(x),pderiv)
            grun.append(grunq)
            
        return np.array(grun)    
    

################################################################################
# 
# Read the frequencies for all geometries where the gruneisen parameters must be
# calculated. This depends on the direction (along a, along c, etc.)
# According to the direction chosen, start,stop,step must be given to loop over 
# all geometries as listed in the file containing the energies 
#
#
# More work to do: entend to other ibrav types, etc.
def compute_grun(ngeo,celldmsx,inputfilefreq,ibrav=4,ext=False):
    if ibrav==4:    # Hexagonal systems
        cgeo = find_geocenters(ngeo)
        
        # read frequencies and q points for every geometry
        startgeo = 1
        endgeo = ngeo[0]*ngeo[2]+1
        stepgeo = 1
        if ext:
            weightsgeo, freqgeo = read_freq_ext_geo(inputfilefreq,range(startgeo,endgeo,stepgeo))
            nq = weightsgeo.shape[1] # total number of q points read
        else:
            qgeo, freqgeo = read_freq_geo(inputfilefreq,range(startgeo,endgeo,stepgeo))
            nq = qgeo.shape[1] # total number of q points read
        
        modes = freqgeo.shape[2]  # number of frequency modes
                
        ################################################################################
        # First along the a axis
        startgeo = cgeo[2]*ngeo[0]
        endgeo = cgeo[2]*ngeo[0]+ngeo[0]
        stepgeo = 1
        
        print ("\n"+80*"#"+"\n"+80*"#")
        print ("Calculating Gruneisen parameter along a axis... \n")
        grun_along_a = compute_grun_along_one_direction(nq,modes,ngeo,cgeo,celldmsx,freqgeo,range(startgeo,endgeo,stepgeo),0)
        print ("Done! \n")
        if ext:
            write_freq_ext(weightsgeo[0,:],grun_along_a,"output_grun_along_a_ext1Dfit")
            pass
        else:
            write_freq(qgeo,grun_along_a,"output_grun_along_a1Dfit")
            
 
        ################################################################################
        # now along c
        startgeo = cgeo[0]
        endgeo = ngeo[0]*ngeo[2]
        stepgeo = ngeo[0]
            
        print ("\n"+80*"#"+"\n"+80*"#")
        print ("Calculating Gruneisen parameter along c axis... \n")    
        grun_along_c = compute_grun_along_one_direction(nq,modes,ngeo,cgeo,celldmsx,freqgeo,range(startgeo,endgeo,stepgeo),2)
        print ("Done! \n")
        if ext:
            write_freq_ext(weightsgeo[0,:],grun_along_c,"output_grun_along_c_ext1Dfit")
            pass
        else:
            write_freq(qgeo,grun_along_c,"output_grun_along_c1Dfit")


    
################################################################################
#   MAIN
################################################################################
#

if __name__ == "__main__":
    
    start_time = time.time()
  
    # Default command line parameters
    inputfileEtot = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    inputfilefreq = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/save_frequencies.dat.g"
    
    #inputfileEtot = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    #inputfilefreq = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/save_frequencies.dat.g"
    
    #inputfilefreq = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/frequencies/output_frq.dat.g"
    
    # Read the energies 
    celldmsx, Ex = read_Etot(inputfileEtot)

    ngeo = np.zeros(6,dtype=np.int)
    ngeo[0]=5
    ngeo[2]=5
  
    # Compute the Gruneisen parameters if you haven't done it before
    compute_grun(ngeo,celldmsx,inputfilefreq,ibrav=4,ext=True)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print ("Finished. Elapsed time: "+str(elapsed_time)+" s.")

