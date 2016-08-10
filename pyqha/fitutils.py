# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import minimize


################################################################################

def print_polynomial(a,ibrav=4):
    """
    This function prints the fitted polynomial, either quartic or quadratic
    """
    if ibrav==4:
        print ("\nFitted polynomial is: \n")
        if (len(a)==6):
            print ("p(x1,x2) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2 + "
            +str(a[3])+" * x2 + "+str(a[4])+" * x2^2 + "+str(a[5])+" *x1*x2\n")
        elif (len(a)==15):
            print ("p(x1,x2) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2 + "
            +str(a[3])+" * x1^3 + "+str(a[4])+" * x1^4 +\n" 
            +str(a[5])+" *x2 + "+str(a[6])+" *x2^2 + "+str(a[7])+" *x2^3 + "+str(a[8])+" *x2^4 +\n" 
            +str(a[9])+" *x1*x2 + "+str(a[10])+" *x1*x2^2 + "+str(a[11])+" *x1*x2^3 + "+str(a[12])
            +" *x1^2*x2 + "+str(a[13])+" *x1^2*x2^2 +\n"+str(a[14])+" *x1^3*x2\n")
        else:
            print ("Polynomial order not implemented")


################################################################################

def print_data(x,y,results,A,ibrav,ylabel="E"):   
    """
    This function prints the data and the fitted results 
    ylabel can be "E", "Fvib", "Cxx", etc. so that can be used for different
    fitted quantities
    """
    if ibrav==4:
        print ("a","\t\t\t","c","\t\t\t",ylabel,"\t\t\t",ylabel+"fit","\t\t\t",ylabel+"-"+ylabel+"fit")
        for i in range(0,len(y)):
            s=sum(results[0]*A[i])
            print ("{:.10e}".format(x[i,0]),"\t", "{:.10e}".format(x[i,2]),"\t", 
            "{:.10e}".format(y[i]),"\t", "{:.10e}".format(s),"\t", "{:.10e}".format(y[i]-s))
    else:
        print ("ibrav not implememnted yet")

################################################################################

def fit_anis(celldmsx, Ex, ibrav=1, out=False, type="quadratic", ylabel="Etot"):  
    """
    An auxiliary function for handling fitting in the anisotropic case
    """
    if out:
        print (type+" fit")
        if type=="quadratic":
            a, chi = fit_quadratic(celldmsx, Ex, ibrav, out, ylabel)
        elif type=="quartic":
            a, chi = fit_quartic(celldmsx, Ex, ibrav, out, ylabel)
        else:
            print ("Fitting type not implemented")
            return None, None

        if chi!=None:
            print_polynomial(a)
            print ("Chi squared: ",chi,"\n")
            return a, chi 

        return a, None
    
    else:
        if type=="quadratic":
            a, chi = fit_quadratic(celldmsx, Ex, ibrav, False, ylabel)
        elif type=="quartic":
            a, chi = fit_quartic(celldmsx, Ex, ibrav, False, ylabel)
        else:
            return None, None

        if chi!=None:
            return a, chi 

        return a, None
    
################################################################################

def fit_quadratic(x,y,ibrav=4,out=False,ylabel="E"):
    """
    This is the function for fitting with a quadratic polynomial

    The most general fitting multidimensional quadratic polynomial for a triclinic
    system is:
    a1 + a2 x1 + a3 x1^2 + a4  x2 + a5  x2^2 + a6 x1*x2 +    
    + a7  x3 + a8  x3^2 + a9  x1*x3 + a10 x2*x3 +         
    + a11 x4 + a12 x4^2 + a13 x1*x4 + a14 x2*x4 + a15 x3*x4 +        
    + a16 x5 + a17 x5^2 + a18 x1*x5 + a19 x2*x5 + a20 x3*x5 + a21 x4*x5 
    + a22 x6 + a23 x6^2 + a24 x1*x6 + a25 x2*x6 + a26 x3*x6 + a27 x4*x6 + a28 x5*x6

    ONLY THE HEXAGONAL AND GENERAL CASE ARE IMPLEMENTED, more to be done

    The input variable x is a matrix ngeo*6, where
    x[:,0] is the set of a values  
    x[:,1] is the set of b values  
    x[:,2] is the set of c values  
    x[:,3] is the set of alpha values  
    x[:,4] is the set of beta values   
    x[:,5] is the set of gamma values  
 
    """
    
    # Create the auxiliary A numpy matrix for the fitting coefficients (quadratic polynomial)
    if ibrav==1 or ibrav==2 or ibrav==3:  # cubic
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0]])   
    elif ibrav==4:   # hexagonal systems
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,2], x[:,2]*x[:,2], x[:,0]*x[:,2]])       
    elif ibrav==14:          # triclinic systems (possibly)
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,1], x[:,1]*x[:,1], x[:,0]*x[:,1],
                    x[:,2], x[:,2]*x[:,2], x[:,0]*x[:,2], x[:,1]*x[:,2],
                    x[:,3], x[:,3]*x[:,3], x[:,0]*x[:,3], x[:,1]*x[:,3], x[:,2]*x[:,3],
                    x[:,4], x[:,4]*x[:,4], x[:,0]*x[:,4], x[:,1]*x[:,4], x[:,2]*x[:,4], x[:,3]*x[:,4],
                    x[:,5], x[:,5]*x[:,5], x[:,0]*x[:,5], x[:,1]*x[:,5], x[:,2]*x[:,5], x[:,3]*x[:,5], x[:,4]*x[:,5]])
    else:          # ibrav not implemented
        return None, None
    
    A = A.T

    # Fit with the numpy linear solver... returns the polynomial coefficients in a
    results = np.linalg.lstsq(A, y)
  
    # Print some extra output if requested
    if out:
        print_data(x,y,results,A,ibrav,ylabel)
     
    if len(results[1])==0:   # if not residual is return, np.linalg.lstsq has not converged
        print ("WARNING! Problems in quadratic fit... convergence not achieved\n")
        return results[0], None

    return results[0], results[1][0]   # first parameter is the vector of coefficients,
                                       # second is the chi squared


################################################################################

def fit_quartic(x,y,ibrav=4,out=False,ylabel="E"):
    """
    This is the function for fitting with a quartic polynomial

    The most general fitting multidimensional quadratic polynomial for a triclinic
    system is:

    a1 + a2  x1 + a3  x1^2 + a4  x1^3 + a5 x1^4        
      + a6  x2 + a7  x2^2 + a8  x2^3 + a9 x2^4       
      + a10 x1*x2 + a11 x1*x2^2 + a12  x1*x2^3
                           + a13 x1^2*x2 + a14  x1^2*x2^2 
                           + a15 x1^3*x2 
      + a16 x3 + a17 x3^2 + a18 x3^3 + a19 x3^4       
      + a20 x1*x3 + a21 x1*x3^2 + a22  x1*x3^3
                           + a23 x1^2*x3 + a24  x1^2*x3^2 
                           + a25 x1^3*x3 
      + a26 x2*x3 + a27 x2*x3^2 + a28  x2*x3^3
                           + a29 x2^2*x3 + a30  x2^2*x3^2 
                           + a31 x2^3*x3 
      + a32 x1 * x2 * x3 + a33 x1 * x2^2 * x3
      + a34 x1 * x2 * x3^2 + a35 x1^2 * x2 * x3
      + a36 x4 + a37 x4^2 + a38 x4^3 + a39 x4^4       
      + a40 x1*x4 + a41 x1*x4^2 + a42  x1*x4^3
                           + a43 x1^2*x4 + a44  x1^2*x4^2 
                           + a45 x1^3*x4 
      + a46 x2*x4 + a47 x2*x4^2 + a48  x2*x4^3
                           + a49 x2^2*x4 + a50  x2^2*x4^2 
                           + a51 x2^3*x4 
      + a52 x3*x4 + a53 x3*x4^2 + a54  x3*x4^3
                           + a55 x3^2*x4 + a56  x3^2*x4^2 
                           + a57 x3^3*x4 
      + a58 x1 * x2 * x4 + a59 x1 * x2^2 * x4
      + a60 x1 * x2 * x4^2 + a61 x1^2 * x2 * x4
      + a62 x1 * x3 * x4 + a63 x1 * x3^2 * x4
      + a64 x1 * x3 * x4^2 + a65 x1^2 * x3 * x4
      + a66 x2 * x3 * x4 + a67 x2 * x3^2 * x4
      + a68 x2 * x3 * x4^2 + a69 x2^2 * x3 * x4
      + a70 x1 * x2 * x3 * x4
      + a71 x5 + a72 x5^2 + a73 x5^3 + a74 x5^4       
      + a75 x1*x5 + a76 x1*x5^2 + a77 x1*x5^3
                           + a78 x1^2*x5 + a79 x1^2*x5^2 
                           + a80 x1^3*x5 
      + a81 x2*x5 + a82 x2*x5^2 + a83 x2*x5^3
                           + a84 x2^2*x5 + a85 x2^2*x5^2 
                           + a86 x2^3*x5 
      + a87 x3*x5 + a88 x3*x5^2 + a89 x3*x5^3
                           + a90 x3^2*x5 + a91 x3^2*x5^2 
                           + a92 x3^3*x5 
      + a93 x4*x5 + a94 x4*x5^2 + a95 x4*x5^3
                           + a96 x4^2*x5 + a97 x4^2*x5^2 
                           + a98 x4^3*x5 
      + a99  x1 * x2 * x5 + a100 x1 * x2^2 * x5
      + a101 x1 * x2 * x5^2 + a102 x1^2 * x2 * x5
      + a103 x1 * x3 * x5 + a104 x1 * x3^2 * x5
      + a105 x1 * x3 * x5^2 + a106 x1^2 * x3 * x5
      + a107 x1 * x4 * x5 + a108 x1 * x4^2 * x5
      + a109 x1 * x4 * x5^2 + a110 x1^2 * x4 * x5
      + a111 x2 * x3 * x5 + a112 x2 * x3^2 * x5
      + a113 x2 * x3 * x5^2 + a114 x2^2 * x3 * x5
      + a115 x2 * x4 * x5 + a116 x2 * x4^2 * x5
      + a117 x2 * x4 * x5^2 + a118 x2^2 * x4 * x5
      + a119 x3 * x4 * x5 + a120 x3 * x4^2 * x5
      + a121 x3 * x4 * x5^2 + a122 x3^2 * x4 * x5
      + a123 x1 * x2 * x3 * x5
      + a124 x1 * x2 * x4 * x5
      + a125 x1 * x3 * x4 * x5
      + a126 x2 * x3 * x4 * x5
      + a127 x6 + a128 x6^2 + a129 x6^3 + a130 x6^4   
      + a131 x1*x6 + a132 x1*x6^2 + a133 x1*x6^3
                         + a134 x1^2*x6 + a135 x1^2*x6^2 
                         + a136 x1^3*x6 
      + a137 x2*x6 + a138 x2*x6^2 + a139 x2*x6^3
                         + a140 x2^2*x6 + a141 x2^2*x6^2 
                         + a142 x2^3*x6 
      + a143 x3*x6 + a144 x3*x6^2 + a145 x3*x6^3
                         + a146 x3^2*x6 + a147 x3^2*x6^2 
                         + a148 x3^3*x6 
      + a149 x4*x6 + a150 x4*x6^2 + a151 x4*x6^3
                         + a152 x4^2*x6 + a153 x4^2*x6^2 
                         + a154 x4^3*x6 
      + a155 x5*x6 + a156 x5*x6^2 + a157 x5*x6^3
                         + a158 x5^2*x6 + a159 x5^2*x6^2 
                         + a160 x5^3*x6 
      + a161 x1 * x2 * x6 + a162 x1 * x2^2 * x6
      + a163 x1 * x2 * x6^2 + a164 x1^2 * x2 * x6
      + a165 x1 * x3 * x6 + a166 x1 * x3^2 * x6
      + a167 x1 * x3 * x6^2 + a168 x1^2 * x3 * x6
      + a169 x1 * x4 * x6 + a170 x1 * x4^2 * x6
      + a171 x1 * x4 * x6^2 + a172 x1^2 * x4 * x6
      + a173 x1 * x5 * x6 + a174 x1 * x5^2 * x6
      + a175 x1 * x5 * x6^2 + a176 x1^2 * x5 * x6
      + a177 x2 * x3 * x6 + a178 x2 * x3^2 * x6
      + a179 x2 * x3 * x6^2 + a180 x2^2 * x3 * x6
      + a181 x2 * x4 * x6 + a182 x2 * x4^2 * x6
      + a183 x2 * x4 * x6^2 + a184 x2^2 * x4 * x6
      + a185 x2 * x5 * x6 + a186 x2 * x5^2 * x6
      + a187 x2 * x5 * x6^2 + a188 x2^2 * x5 * x6
      + a189 x3 * x4 * x6 + a190 x3 * x4^2 * x6
      + a191 x3 * x4 * x6^2 + a192 x3^2 * x4 * x6
      + a193 x3 * x5 * x6 + a194 x3 * x5^2 * x6
      + a195 x3 * x5 * x6^2 + a196 x3^2 * x5 * x6
      + a197 x4 * x5 * x6 + a198 x4 * x5^2 * x6
      + a199 x4 * x5 * x6^2 + a200 x4^2 * x5 * x6
      + a201 x1 * x2 * x3 * x6
      + a202 x1 * x2 * x4 * x6
      + a203 x1 * x2 * x5 * x6
      + a204 x1 * x3 * x4 * x6
      + a205 x1 * x3 * x5 * x6
      + a206 x1 * x4 * x5 * x6
      + a207 x2 * x3 * x4 * x6
      + a208 x2 * x3 * x5 * x6
      + a209 x2 * x4 * x5 * x6
      + a210 x3 * x4 * x5 * x6

    ONLY THE HEXAGONAL CASE IS IMPLEMENTED, more to be done


    The input variable x is a matrix ngeo*6, where
    x[:,0] is the set of a values  
    x[:,1] is the set of b values  
    x[:,2] is the set of c or c/a values  
    x[:,3] is the set of alpha values  
    x[:,4] is the set of beta values   
    x[:,5] is the set of gamma values  
 
    """
    
    # Create the auxiliary A numpy matrix for the fitting coefficients (quadratic polynomial)
    if ibrav==4:   # hexagonal systems
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0]*x[:,0],
               x[:,2], x[:,2]*x[:,2], x[:,2]*x[:,2]*x[:,2], x[:,2]*x[:,2]*x[:,2]*x[:,2], 
               x[:,0]*x[:,2], x[:,0]*x[:,2]*x[:,2], x[:,0]*x[:,2]*x[:,2]*x[:,2], x[:,0]*x[:,0]*x[:,2], x[:,0]*x[:,0]*x[:,2]*x[:,2],
               x[:,0]*x[:,0]*x[:,0]*x[:,2]])
    else:          # ibrav not implemented
        return None, None
    
      # Create the auxiliary A numpy matrix for the fitting coefficients (quadratic polynomial)

    A = A.T

    # Fit with the numpy linear solver... returns the polynomial coefficients in a
    results = np.linalg.lstsq(A, y)
  
    # Print some extra output if requested
    if out:
        print_data(x,y,results,A,ibrav,ylabel)
     
    if len(results[1])==0:   # if not residual is return, np.linalg.lstsq has not converged
        print ("WARNING! Problems in quadratic fit... convergence not achieved\n")
        return results[0], None

    return results[0], results[1][0]   # first parameter is the vector of coefficients,
                                       # second is the chi squared

################################################################################

def expand_quadratic_to_quartic(a):
    """
    This function gets a vector of coefficients from a quadratic fit and turns it
    into a vector of coeffients as from a quartic fit (extra coeffients are set to zero)
    """
    b=[a[0], a[1], a[2], 0.0, 0.0, 
    a[3], a[4], 0.0, 0.0,
    a[5], 0.0, 0.0, 0.0, 0.0,
    0.0]
    return b