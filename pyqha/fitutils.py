#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.
import numpy as np
from scipy.optimize import minimize


################################################################################

def print_polynomial(a,ibrav=4):
    """
    This function prints the fitted polynomial, either quartic or quadratic
    """
    if ibrav in (1,2,3):
        print ("\nFitted polynomial is: \n")
        if (len(a)==3):
            print ("p(x1) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2\n")
        elif (len(a)==5):
            print ("p(x1) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2 + "
            +str(a[3])+" * x1^3 + "+str(a[4])+" * x1^4 +\n")
        else:
            print ("Polynomial order not implemented")
    elif ibrav==4:
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
    if ibrav in (1,2,3):
        print ("a or V","\t\t\t",ylabel,"\t\t\t",ylabel+"fit","\t\t\t",ylabel+"-"+ylabel+"fit")
        for i in range(0,len(y)):
            s=sum(results[0]*A[i])
            print ("{:.10e}".format(x[i,0]),"\t",  
            "{:.10e}".format(y[i]),"\t", "{:.10e}".format(s),"\t", "{:.10e}".format(y[i]-s))
    elif ibrav==4:
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
    into a vector of coeffients as from a quartic fit (extra coeffients are set to zero).
    """
    b=[a[0], a[1], a[2], 0.0, 0.0, 
    a[3], a[4], 0.0, 0.0,
    a[5], 0.0, 0.0, 0.0, 0.0,
    0.0]
    return b