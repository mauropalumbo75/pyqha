#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
This submodule contains all functions relevant for fitting anistotropic quantities
(energies, frequencies, etc.). 
"""

import numpy as np
from scipy.optimize import minimize


################################################################################

def print_polynomial(a,ibrav=4):
    """
    This function prints the fitted polynomial, either quartic or quadratic
    """
    if ibrav in (1,2,3):                # cubic systems (a,a,a)
        print ("\nFitted polynomial is: \n")
        if (len(a)==3):
            print ("p(x1) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2\n")
        elif (len(a)==5):
            print ("p(x1) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2 + "
            +str(a[3])+" * x1^3 + "+str(a[4])+" * x1^4 +\n")
        else:
            print ("Polynomial order not implemented")
    elif ibrav in (4,6,7):              # hexagonal or tetragonal systems (a,a,c)
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
    elif ibrav in (8,9,10,11):          # orthorombic systems (a,b,c)
        if (len(a)==10):
            print ("p(x1,x2,x3) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2 + "
            +str(a[3])+" * x2 + "+str(a[4])+" * x2^2 + "+str(a[5])+" *x1*x2\n")
        elif (len(a)==35):
            print ("p(x1,x2,x3) = "+str(a[0])+" + "+str(a[1])+" * x1 + "+str(a[2])+" * x1^2 + "
            +str(a[3])+" * x1^3 + "+str(a[4])+" * x1^4 +\n" 
            +str(a[5])+" *x2 + "+str(a[6])+" *x2^2 + "+str(a[7])+" *x2^3 + "+str(a[8])+" *x2^4 +\n" 
            +str(a[9])+" *x1*x2 + "+str(a[10])+" *x1*x2^2 + "+str(a[11])+" *x1*x2^3 + "+str(a[12])
            +" *x1^2*x2 + "+str(a[13])+" *x1^2*x2^2 +\n"+str(a[14])+" *x1^3*x2\n"
            +str(a[15])+" *x3"+str(a[16])+" *x3^2"+str(a[17])+" *x3^3"+str(a[18])+" *x3^4\n"
            +str(a[19])+" *x1*x3"+str(a[20])+" *x1*x3^2"+str(a[21])+" *x1*x3^3"+str(a[22])+" *x1^2*x3"
            +str(a[23])+" *x1^2*x3^2"+str(a[24])+" *x1^3*x3"+str(a[25])+" *x2*x3"+str(a[26])+" *x2*x3^2\n"
            +str(a[27])+" *x2*x3^3"+str(a[28])+" *x2^2*x3"+str(a[29])+" *x2^2*x3^2"+str(a[30])+" *x2^3*x3\n"
            +str(a[31])+" *x1*x2*x3"+str(a[32])+" *x1^2*x2*x3"+str(a[33])+" *x1*x2^2x3"+str(a[34])+" *x1*x2*x3^2\n"
            )
        else:
            print ("Polynomial order not implemented")
    else:
        print ("ibrav not implememnted yet")

################################################################################

def print_data(x,y,results,A,ibrav,ylabel="E"):   
    """
    This function prints the data and the fitted results 
    ylabel can be "E", "Fvib", "Cxx", etc. so that can be used for different
    fitted quantities
    """
    if ibrav in (1,2,3):                # cubic systems (a,a,a)
        print ("a or V","\t\t\t",ylabel,"\t\t\t",ylabel+"fit","\t\t\t",ylabel+"-"+ylabel+"fit")
        for i in range(0,len(y)):
            s=sum(results[0]*A[i])
            print ("{:.10e}".format(x[i,0]),"\t",  
            "{:.10e}".format(y[i]),"\t", "{:.10e}".format(s),"\t", "{:.10e}".format(y[i]-s))
    elif ibrav in (4,6,7):              # hexagonal or tetragonal systems (a,a,c)
        print ("a","\t\t\t","c","\t\t\t",ylabel,"\t\t\t",ylabel+"fit","\t\t\t",ylabel+"-"+ylabel+"fit")
        for i in range(0,len(y)):
            s=sum(results[0]*A[i])
            print ("{:.10e}".format(x[i,0]),"\t", "{:.10e}".format(x[i,2]),"\t", 
            "{:.10e}".format(y[i]),"\t", "{:.10e}".format(s),"\t", "{:.10e}".format(y[i]-s))
    elif ibrav in (8,9,10,11):          # orthorombic systems (a,b,c)
        print ("a","\t\t\t","b","\t\t\t","c","\t\t\t",ylabel,"\t\t\t",ylabel+"fit","\t\t\t",ylabel+"-"+ylabel+"fit")
        for i in range(0,len(y)):
            s=sum(results[0]*A[i])
            print ("{:.10e}".format(x[i,0]),"\t", "{:.10e}".format(x[i,1]),"\t", "{:.10e}".format(x[i,2]),"\t", 
            "{:.10e}".format(y[i]),"\t", "{:.10e}".format(s),"\t", "{:.10e}".format(y[i]-s))
    else:
        print ("ibrav not implememnted yet")

################################################################################

def fit_anis(celldmsx, Ex, ibrav=4, out=False, type="quadratic", ylabel="Etot"):  
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
    This function fits the *y* values (energies) with a quadratic polynomial of
    up to 3 variables :math:`x1,x2,x3` corresponding to the lattice parameters
    :math:`(a,b,c)`. In the most general form the polynomial is:

    .. math::
      a_1 + a_2 * x_1 + a_3 * x_1^2 + a_4 * x_2 + a_5 * x_2^2 + a_6 * x_1*x_2 +    
      + a_7 * x_3 + a_8 * x_3^2 + a_9 * x_1*x_3 + a_10 * x_2*x_3

    The input variable *x* is a matrix :math:`ngeo*6`, where:
    
    *x[:,0]* is the set of a values  
    
    *x[:,1]* is the set of b values  
    
    *x[:,2]* is the set of c values  
    
    and x[:,3], x[:,4], x[:,5] are all zeros. 
    *ibrav* defines the Bravais lattice, *out* set the output verbosity (*out=True* 
    verbose output), *ylabel* set a label for the quantity in *y* (can be :math:`E_{tot}`
    or :math:`E_{tot}+F_{vib}` for example).
    
    Note 1: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)  
    
    Note 2: the polynomial fit is done using :py:func:`numpy.linalg.lstsq`. 
    Please refer to numpy documentation for further details.
    """
    
    # Create the auxiliary A numpy matrix for the fitting coefficients (quadratic polynomial)
    # A dimension is different according to ibrav
    if ibrav in (1,2,3):                # cubic systems (a,a,a)
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0]])   
    elif ibrav in (4,6,7):              # hexagonal or tetragonal systems (a,a,c)
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,2], x[:,2]*x[:,2], x[:,0]*x[:,2]])       
    elif ibrav in (8,9,10,11):          # orthorombic systems (a,b,c)
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,1], x[:,1]*x[:,1], x[:,0]*x[:,1],
                    x[:,2], x[:,2]*x[:,2], x[:,0]*x[:,2], x[:,1]*x[:,2]])
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
    This function fits the *y* values (energies) with a quartic polynomial of
    up to 3 variables :math:`x1,x2,x3` corresponding to the lattice parameters
    :math:`(a,b,c)`. In the most general form the polynomial is:

    .. math::
      a_1 + a_2 * x_1 + a_3 * x_1^2 + a_4 * x_1^3 + a_5 * x_1^4 + a_6 * x_2 + a_7 * x_2^2 + 
      a_8 * x_2^3 + a_9 * x_2^4 + a_10 * x_1*x_2 + a_11 * x_1*x_2^2 + a_12 * x_1*x_2^3 +
      a_13 * x_1^2*x_2 + a_14 * x_1^2*x_2^2 + a_15 * x_1^3*x_2 + a_16 * x_3 + a_17 * x_3^2 + a_18 * x_3^3 +
    .. math::
      a_19 * x_3^4 + a_20 * x_1*x_3 + a_21 * x_1*x_3^2 + a_22 * x_1*x_3^3 + a_23 * x_1^2*x_3 + a_24 * x_1^2*x_3^2 +
      a_25 * x_1^3*x_3 + a_26 * x_2*x_3 + a_27 * x_2*x_3^2 + a_28 * x_2*x_3^3 + a_29 * x_2^2*x_3 +
      a_30 * x_2^2*x_3^2 + a_31 * x_2^3*x_3 + a_32 * x_1*x_2*x_3 + a_33 * x_1^2*x_2*x_3 +
      a_34 * x_1*x_2^2*x_3 + a_35 * x_1*x_2*x_3^2 

    The input variable *x* is a matrix :math:`ngeo*6`, where:
    
    *x[:,0]* is the set of a values  
    
    *x[:,1]* is the set of b values  
    
    *x[:,2]* is the set of c values  
    
    and x[:,3], x[:,4], x[:,5] are all zeros. 
    *ibrav* defines the Bravais lattice, *out* set the output verbosity (*out=True* 
    verbose output), *ylabel* set a label for the quantity in *y* (can be :math:`E_{tot}`
    or :math:`E_{tot}+F_{vib}` for example).
    
    Note 1: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)  
    
    Note 2: the polynomial fit is done using :py:func:`numpy.linalg.lstsq`. 
    Please refer to numpy documentation for further details.
    """
    
    # Create the auxiliary A numpy matrix for the fitting coefficients (quadratic polynomial)
    # A dimension is different according to ibrav
    if ibrav in (1,2,3):                # cubic systems (a,a,a)
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0]*x[:,0]])        
    elif ibrav in (4,6,7):              # hexagonal or tetragonal systems (a,a,c)
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0]*x[:,0],
               x[:,2], x[:,2]*x[:,2], x[:,2]*x[:,2]*x[:,2], x[:,2]*x[:,2]*x[:,2]*x[:,2], 
               x[:,0]*x[:,2], x[:,0]*x[:,2]*x[:,2], x[:,0]*x[:,2]*x[:,2]*x[:,2], x[:,0]*x[:,0]*x[:,2], x[:,0]*x[:,0]*x[:,2]*x[:,2],
               x[:,0]*x[:,0]*x[:,0]*x[:,2]])
    elif ibrav in (8,9,10,11):          # orthorombic systems (a,b,c)
        A = np.vstack([np.ones(len(x[:,0])), x[:,0], x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0], x[:,0]*x[:,0]*x[:,0]*x[:,0],
               x[:,1], x[:,1]*x[:,1], x[:,1]*x[:,1]*x[:,1], x[:,1]*x[:,1]*x[:,1]*x[:,1], 
               x[:,0]*x[:,1], x[:,0]*x[:,1]*x[:,1], x[:,0]*x[:,1]*x[:,1]*x[:,1], x[:,0]*x[:,0]*x[:,1], x[:,0]*x[:,0]*x[:,1]*x[:,1],
               x[:,0]*x[:,0]*x[:,0]*x[:,1], x[:,2], x[:,2]*x[:,2], x[:,2]*x[:,2]*x[:,2], x[:,2]*x[:,2]*x[:,2]*x[:,2],
               x[:,0]*x[:,2], x[:,0]*x[:,2]*x[:,2], x[:,0]*x[:,2]*x[:,2]*x[:,2], x[:,0]*x[:,0]*x[:,2], x[:,0]*x[:,0]*x[:,2]*x[:,2],
               x[:,0]*x[:,0]*x[:,0]*x[:,2],
               x[:,1]*x[:,2], x[:,1]*x[:,2]*x[:,2], x[:,1]*x[:,2]*x[:,2]*x[:,2], x[:,1]*x[:,1]*x[:,2], x[:,1]*x[:,1]*x[:,2]*x[:,2],
               x[:,1]*x[:,1]*x[:,1]*x[:,2],
               x[:,0]*x[:,1]*x[:,2], x[:,0]*x[:,0]*x[:,1]*x[:,2], x[:,0]*x[:,1]*x[:,1]*x[:,2], x[:,0]*x[:,1]*x[:,2]*x[:,2]])    
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