#!/usr/bin/env python3
#encoding: UTF-8

import sys
import numpy as np
from fitFvib import fitFvib
from fitC import fitC
from minutils import fquadratic, fquartic
from write import write_qha_CT

################################################################################
#   MAIN
################################################################################
#

if __name__ == "__main__":
  
    # Default command line parameters
#    inputfileEtot = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
#    inputfileFvib = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/therm_files/output_therm.dat" 
#    inputpathCx = "Os.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/"  
    
#    inputfileEtot = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
#    inputfileFvib = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/therm_files/output_therm.dat" 
#    inputpathCx = "Re.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/"  
    
    inputfileEtot = "Tc.pz-spn-kjpaw_psl.1.0.0.UPF/energy_files/output_energy1"
    inputfileFvib = "Tc.pz-spn-kjpaw_psl.1.0.0.UPF/therm_files/output_therm.dat" 
    inputpathCx = "Tc.pz-spn-kjpaw_psl.1.0.0.UPF/elastic_constants/"      
    
    ibrav = 4  # default value, to be later changed to 1
#    typeEtot = "quadratic"
#    typeFvib = "quadratic"
#    typeCx = "quadratic"
    typeEtot = "quartic"
    typeFvib = "quartic"
    typeCx = "quartic"
    # Initial guess for Re
#    defaultguess=[5.16928989,  1.61267114]
    # default guess for Tc    
    defaultguess=[ 5.12376396,  1.59903566]
        
    # Read and process command line parameters, if any
    print ("Possible options are: -fileEtot -fileFvib -pathCx -ibrav -fittypeEtot -fittypeFvib -fittypeCx -guess\n")
    print ("Assuming default values for unspecified options: \n")
    if len(sys.argv)>1:
        for i in range(1,len(sys.argv)):
            if sys.argv[i] == "-fileEtot": inputfileEtot = sys.argv[i+1]
            elif sys.argv[i] == "-fileFvib": inputfileFvib = sys.argv[i+1]
            elif sys.argv[i] == "-pathCx": inputpathCx = sys.argv[i+1]
            elif sys.argv[i] == "-ibrav": ibrav = int(sys.argv[i+1])
            elif sys.argv[i] == "-fittypeEtot": typeEtot = sys.argv[i+1]
            elif sys.argv[i] == "-fittypeFvib": typeFvib = sys.argv[i+1]
            elif sys.argv[i] == "-fittypeCx": typeCx = sys.argv[i+1]
            elif sys.argv[i] == "-guess": guess = sys.argv[i+1]


    # Call the processing function
    TT, minT, fminT = fitFvib(inputfileEtot,inputfileFvib,ibrav,typeEtot,typeFvib,defaultguess)

    aC, chiC = fitC(inputfileEtot, inputpathCx, ibrav, typeCx)

    nT = 0
    CTtemp = []
    # Find the elastic constants
    for T in range(0,len(TT)):
        C = []
        for i in range(0,6):
            Ccol = []
            for j in range(0,6):
                if typeCx=="quadratic":
                    Ctemp = fquadratic(minT[T],aC[i][j],ibrav=4)
                elif typeCx=="quartic":
                    Ctemp = fquartic(minT[T],aC[i][j],ibrav=4)  
                Ccol.append(Ctemp)
            C.append(Ccol)
        CTtemp.append(C)
        nT += 1
    
    CT = np.array(CTtemp)
    CT.shape = (nT,6,6)     
    
    write_qha_CT(TT,CT,inputpathCx)
