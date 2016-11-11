#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
This submodule groups all functions relevant for minimizing the energy in the 
anisotropic case. It also contains the functions implementing the fitting
polynomials, which can be "quadratic" or "quartic" and depend on up to 3 variables.
"""


import numpy as np
from scipy.optimize import minimize


################################################################################

def expand_vector(x, ibrav=4):
    """
    Utility function: expands a vector *x*, len(x)<6, into a 6-dim vector according
    to the Bravais lattice type as in *ibrav*
    
    Note: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)
    """
    if (len(x)>6 or x==None):
        return None
    elif (len(x)==6):
        return x
    else:
        if (ibrav==1 or ibrav==2 or ibrav==3):      # cubic systems (a,a,a)
            xnew = np.array([x[0],0.0,0.0,0.0,0.0,0.0])
        elif (ibrav==4 or ibrav==6 or ibrav==7):      # hexagonal or tetragonal systems (a,a,c)
            xnew = np.array([x[0],0.0,x[1],0.0,0.0,0.0])
        elif (ibrav>=8 and ibrav<=11):              # orthorombic systems (a,b,c)
            xnew = np.array([x[0],x[1],x[2],0.0,0.0,0.0])
        else:
            pass
        return xnew
            
    
################################################################################

def contract_vector(x, ibrav=4):
    """
    Utility function: contract a vector *x*, len(x)=6, into a x-dim vector (x<6) 
    according to the Bravais lattice type as in *ibrav*
    
    Note: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)
    """
    if (len(x)!=6):
        return None
    else:
        if (ibrav==1 or ibrav==2 or ibrav==3):      # cubic systems (a,a,a)
            xnew = np.array([x[0]])
        elif (ibrav==4 or ibrav==6 or ibrav==7):    # hexagonal or tetragonal systems (a,a,c)
            xnew = np.array([x[0],x[2]])
        elif (ibrav>=8 and ibrav<=11):              # orthorombic systems (a,b,c)
            xnew = np.array([x[0],x[1],x[2]])        
        else:
            pass
        return xnew

################################################################################
           
def fquadratic(x,a,ibrav=4):
    """
    This function implements the quadratic polynomials for fitting and miminizing
    according to the Bravais lattice type as in *ibrav*.
    *x* is the vector with input coordinates :math:`(a,b,c,alpha,beta,gamma)`,
    *a* is the vector
    with the polynomial coeffients. The dimension of *a* depends on
    *ibrav*. *x* has always 6 elements with zeros for those not used according to
    *ibrav*.
    
    Note: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)  
    """
    if (ibrav==1 or ibrav==2 or ibrav==3):          # cubic systems (a,a,a)
        return a[0]+a[1]*x[0]+a[2]*x[0]**2
    elif (ibrav==4 or ibrav==6 or ibrav==7):        # hexagonal or tetragonal systems (a,a,c)
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[2]+a[4]*x[2]**2+a[5]*x[0]*x[2]
    elif (ibrav>=8 and ibrav<=11):                  # orthorombic systems (a,b,c)
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
        +a[6]*x[2]+a[7]*x[2]**2+a[8]*x[0]*x[2]+a[9]*x[1]*x[2]
    else:
        return None
    
        # Not implemented, just for future reference
        # This is the most general quadratic polynomial for the triclinic case with (a,b,c,alpha,beta,gamma) variables
        #a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
        #+a[6]*x[2]+a[7]*x[2]**2+a[8]*x[0]*x[2]+a[9]*x[1]*x[2]        
        #+a[10]*x[3]+a[11]*x[3]**2+a[12]*x[0]*x[3]+a[13]*x[1]*x[3]+a[14]*x[2]*x[3]
        #+a[15]*x[4]+a[16]*x[4]**2+a[17]*x[0]*x[4]+a[18]*x[1]*x[4]+a[19]*x[2]*x[4]+a[20]*x[3]*x[4]
        #+a[21]*x[5]+a[22]*x[5]**2+a[23]*x[0]*x[5]+a[24]*x[1]*x[5]+a[25]*x[2]*x[5]+a[26]*x[3]*x[5]+a[26]*x[4]*x[5]  


def fquadratic_der(x,a,ibrav=4):
    """
    This function implements the first derivatives of the quadratic polynomials 
    for fitting and miminizing according to the Bravais lattice type as in *ibrav*. 
    *x* is the vector with input coordinates :math:`(a,b,c,alpha,beta,gamma)`,
    *a* is the vector
    with the polynomial coeffients. The dimension of *a* depends on
    *ibrav*. *x* has always 6 elements with zeros for those not used according to
    *ibrav*.
    The derivatives are returned as a numpy vector of 6 elements, each element 
    being a derivative with respect to :math:`(a,b,c,alpha,beta,gamma)`.
    
    Note: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)
    """
    if (ibrav==1 or ibrav==2 or ibrav==3):          # cubic systems (a,a,a)
        return np.array([
        a[1]+2.0*a[2]*x[0],
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
        ])
    elif (ibrav==4 or ibrav==6 or ibrav==7):        # hexagonal or tetragonal systems (a,a,c)
        return np.array([
        a[1]+2.0*a[2]*x[0]+a[5]*x[2],
        0.0,
        a[3]+2.0*a[4]*x[2]+a[5]*x[0],
        0.0,
        0.0,
        0.0
        ])
    elif (ibrav>=8 and ibrav<=11):                  # orthorombic systems (a,b,c)
        return np.array([ 
        a[1]+2.0*a[2]*x[0]+a[5]*x[1]+a[8]*x[2],        
        a[3]+2.0*a[4]*x[1]+a[5]*x[0]+a[9]*x[2], 
        a[6]+2.0*a[7]*x[2]+a[8]*x[0]+a[9]*x[1], 
        0.0,
        0.0, 
        0.0
        ]) 
    else:
        return None
    
          
def fquartic(x,a,ibrav=4):
    """
    This function implements the quartic polynomials for fitting and miminizing
    according to the Bravais lattice type as in *ibrav*.
    *x* is the vector with input coordinates :math:`(a,b,c,alpha,beta,gamma)`,
    *a* is the vector
    with the polynomial coeffients. The dimension of *a* depends on
    *ibrav*. *x* has always 6 elements with zeros for those not used according to
    *ibrav*.
    
    Note: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)  
    """
    if (ibrav==1 or ibrav==2 or ibrav==3):          # cubic systems (a,a,a)
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4)    
    elif (ibrav==4 or ibrav==6 or ibrav==7):        # hexagonal or tetragonal systems (a,a,c)
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4+ 
           a[5]*x[2]+a[6]*x[2]**2+a[7]*x[2]**3+a[8]*x[2]**4+ 
           a[9]*x[0]*x[2]+a[10]*x[0]*x[2]**2+a[11]*x[0]*x[2]**3+a[12]*x[0]**2*x[2]+a[13]*x[0]**2*x[2]**2+ 
           a[14]*x[0]**3*x[2])
    elif (ibrav>=8 and ibrav<=11):                  # orthorombic systems (a,b,c)
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4+ 
           a[5]*x[1]+a[6]*x[1]**2+a[7]*x[1]**3+a[8]*x[1]**4+ 
           a[9]*x[0]*x[1]+a[10]*x[0]*x[1]**2+a[11]*x[0]*x[1]**3+a[12]*x[0]**2*x[1]+a[13]*x[0]**2*x[1]**2+ 
           a[14]*x[0]**3*x[1]+
           a[15]*x[2]+a[16]*x[2]**2+a[17]*x[2]**3+a[18]*x[2]**4+
           a[19]*x[0]*x[2]+a[20]*x[0]*x[2]**2+a[21]*x[0]*x[2]**3+a[22]*x[0]**2*x[2]+a[23]*x[0]**2*x[2]**2+
           a[24]*x[0]**3*x[2]+
           a[25]*x[1]*x[2]+a[26]*x[1]*x[2]**2+a[27]*x[1]*x[2]**3+a[28]*x[1]**2*x[2]+a[29]*x[1]**2*x[2]**2+
           a[30]*x[1]**3*x[2]+
           a[31]*x[0]*x[1]*x[2]+a[32]*x[0]**2*x[1]*x[2]+a[33]*x[0]*x[1]**2*x[2]+a[34]*x[0]*x[1]*x[2]**2
           )       
    else:
        return None
    
    
def fquartic_der(x,a,ibrav=4):
    """
    This function implements the first derivatives of the quadratic polynomials 
    for fitting and miminizing according to the Bravais lattice type as in *ibrav*. 
    The derivatives are returned as a numpy vector of 6 elements, each element 
    being a derivative with respect to :math:`(a,b,c,alpha,beta,gamma)`.
    *x* is the vector with input coordinates :math:`(a,b,c,alpha,beta,gamma)`,
    *a* is the vector
    with the polynomial coeffients. The dimension of *a* depends on
    *ibrav*. *x* has always 6 elements with zeros for those not used according to
    *ibrav*.
    
    Note: implemented for cubic (*ibrav=1,2,3*), hexagonal (*ibrav=4*), 
    tetragonal (*ibrav=6,7*), orthorombic (*ibrav=8,9,10,11*)
    """
    if (ibrav==1 or ibrav==2 or ibrav==3):          # cubic systems (a,a,a)
        return np.array([
        a[1]+2.0*a[2]*x[0]+3.0*a[3]*x[0]**2+4.0*a[4]*x[0]**3,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
        ])    
    elif (ibrav==4 or ibrav==6 or ibrav==7):        # hexagonal or tetragonal systems (a,a,c)
        return np.array([
        a[1]+2.0*a[2]*x[0]+3.0*a[3]*x[0]**2+4.0*a[4]*x[0]**3+ 
        a[9]*x[2]+a[10]*x[2]**2+a[11]*x[2]**3+2.0*a[12]*x[0]*x[2]+2.0*a[13]*x[0]*x[2]**2+ 
        3.0*a[14]*x[0]**2*x[2],
        0.0,
        a[5]+2.0*a[6]*x[2]+3.0*a[7]*x[2]**2+4.0*a[8]*x[2]**3+ 
        a[9]*x[0]+2.0*a[10]*x[0]*x[2]+3.0*a[11]*x[0]*x[2]**2+a[12]*x[0]**2+2.0*a[13]*x[0]**2*x[2]+ 
        a[14]*x[0]**3,
        0.0,
        0.0,
        0.0
        ])
    elif (ibrav>=8 and ibrav<=11):                  # orthorombic systems (a,b,c)
        return np.array([
        a[1]+2.0*a[2]*x[0]+3.0*a[3]*x[0]**2+4.0*a[4]*x[0]**3+ 
        a[9]*x[1]+a[10]*x[1]**2+a[11]*x[1]**3+2.0*a[12]*x[0]*x[1]+2.0*a[13]*x[0]*x[1]**2+ 
        3.0*a[14]*x[0]**2*x[1]+
        a[19]*x[2]+a[20]*x[2]**2+a[21]*x[2]**3+2.0*a[22]*x[0]*x[2]+2.0*a[23]*x[0]*x[2]**2+
        3.0*a[24]*x[0]**2*x[2]+
        a[31]*x[1]*x[2]+2.0*a[32]*x[0]*x[1]*x[2]+a[33]*x[1]**2*x[2]+a[34]*x[1]*x[2]**2,
        
        a[5]+2.0*a[6]*x[1]+3.0*a[7]*x[1]**2+4.0*a[8]*x[1]**3+ 
        a[9]*x[0]+2.0*a[10]*x[0]*x[1]+3.0*a[11]*x[0]*x[1]**2+a[12]*x[0]**2+2.0*a[13]*x[0]**2*x[1]+ 
        a[14]*x[0]**3+
        a[25]*x[2]+a[26]*x[2]**2+a[27]*x[2]**3+2.0*a[28]*x[1]*x[2]+2.0*a[29]*x[1]*x[2]**2+
        3.0*a[30]*x[1]**2*x[2]+
        a[31]*x[0]*x[2]+a[32]*x[0]**2*x[2]+2.0*a[33]*x[0]*x[1]*x[2]+a[34]*x[0]*x[2]**2, 
        
        a[15]+2.0*a[16]*x[2]+3.0*a[17]*x[2]**2+4.0*a[18]*x[2]**3+
        a[19]*x[0]+2.0*a[20]*x[0]*x[2]+3.0*a[21]*x[0]*x[2]**2+a[22]*x[0]**2+2.0*a[23]*x[0]**2*x[2]+
        a[24]*x[0]**3+
        a[25]*x[1]+2.0*a[26]*x[1]*x[2]+3.0*a[27]*x[1]*x[2]**2+a[28]*x[1]**2+2.0*a[29]*x[1]**2*x[2]+
        a[30]*x[1]**3+
        a[31]*x[0]*x[1]+a[32]*x[0]**2*x[1]+a[33]*x[0]*x[1]**2+2.0*a[34]*x[0]*x[1]*x[2],
        
        0.0,
        0.0,
        0.0
        ])
    else:
        return None
  
     
################################################################################

def find_min(a,ibrav,type,guess=None,method="BFGS",minoptions={}):
    """
    An auxiliary function for handling the minimum search.    
    """
    if type=="quadratic":
        xmin=find_min_quadratic(a,ibrav,guess,method,minoptions)
        fmin=fquadratic(xmin,a)
    elif type=="quartic":
        xmin=find_min_quartic(a,ibrav,guess,method,minoptions)
        fmin=fquartic(xmin,a)
    else:
        return None
    print ("Minimun "+type+": ",xmin,"\tEnergy at the minimum: "+"{:.20e}".format(fmin)+"\n") 

    return xmin, fmin

################################################################################

def find_min_quadratic(a,ibrav,guess,method,minoptions):
    """
    This is the function for finding the minimum of the quadratic polynomial
    """
    
    # Redefine here some functions for the quadratic polynomials and gradients.
    # Here the input vector x is reduce in dimension according to ibrav
    # Not ideal, but it works for now
    def fquadratic1(x):             # cubic systems (a,a,a)
        return a[0]+a[1]*x[0]+a[2]*x[0]**2
    
    def fquadratic4(x):             # hexagonal or tetragonal systems (a,a,c)
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
    
    def fquadratic8(x):             # orthorombic systems (a,b,c)
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
        +a[6]*x[2]+a[7]*x[2]**2+a[8]*x[0]*x[2]+a[9]*x[1]*x[2] 
        
    # Gradients:
    def fquadratic_der_1(x):         # cubic systems (a,a,a)
        return np.array([a[1]+2.0*a[2]*x[0]])

    def fquadratic_der_4(x):         # hexagonal or tetragonal systems (a,a,c)
        return np.array([a[1]+2.0*a[2]*x[0]+a[5]*x[1], a[3]+2.0*a[4]*x[1]+a[5]*x[0]])
    
    def fquadratic_der_8(x):         # orthorombic systems (a,b,c)
        return np.array([a[1]+2.0*a[2]*x[0]+a[5]*x[1]+a[8]*x[2],        
        a[3]+2.0*a[4]*x[1]+a[5]*x[0]+a[9]*x[2], 
        a[6]+2.0*a[7]*x[2]+a[8]*x[0]+a[9]*x[1], ])
        
    # Hessians:  
    def fquadratic_der2_1(x):         # cubic systems (a,a,a)
        return np.array([2.0*a[2]])

    def fquadratic_der2_4(x):         # hexagonal or tetragonal systems (a,a,c)
        return np.array([
        [ 2.0*a[2] , a[5]       ],
        [ a[5]     , 2.0*a[4]   ]
        ])
    
    def fquadratic_der2_8(x):         # orthorombic systems (a,b,c)
        return np.array([
        [ 2.0*a[2]   , a[5]        , a[8]      ],        
        [ a[5]       , 2.0*a[4]    , a[9]      ],  
        [ a[8]       , a[9]        , 2.0*a[7]  ]
        ])
        
        
    # First set the initial set, if given
    if guess!=None: 
        x_in = contract_vector(guess,ibrav)
    else:
        x_in = contract_vector([0.,0.,0.,0.,0.,0.],ibrav)

    #  Find the minimun using minimize from scipy.optimize with the gradient
    #  whatever algorithm with the gradient is more stable and I recommend using it 
    if (ibrav==1 or ibrav==2 or ibrav==3):          # cubic systems (a,a,a)
        res_quadratic = minimize(fquadratic1,x_in, jac=fquadratic_der_1, hess=fquadratic_der2_1, method=method, options=minoptions)
    elif (ibrav==4 or ibrav==6 or ibrav==7):        # hexagonal or tetragonal systems (a,a,c)
        res_quadratic = minimize(fquadratic4,x_in, jac=fquadratic_der_4, hess=fquadratic_der2_4, method=method, options=minoptions)
    elif (ibrav>=8 and ibrav<=11):                  # orthorombic systems (a,b,c)
        res_quadratic = minimize(fquadratic8,x_in, jac=fquadratic_der_8, hess=fquadratic_der2_8, method=method, options=minoptions)
  
    if (not res_quadratic.success):
        print ("WARNING! Problems in extremum finding: ",res_quadratic.message)
  
    return expand_vector(res_quadratic.x,ibrav)


################################################################################

def find_min_quartic(a,ibrav,guess,method,minoptions):
    """
    This is the function for finding the minimum of the quartic polynomial
    """
    
    # Redefine here some functions for the quadratic polynomials and gradients.
    # Not ideal, but it works for now
    def fquartic1(x):
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4)
    
    def fquartic4(x):
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4+ 
        a[5]*x[1]+a[6]*x[1]**2+a[7]*x[1]**3+a[8]*x[1]**4+ 
        a[9]*x[0]*x[1]+a[10]*x[0]*x[1]**2+a[11]*x[0]*x[1]**3+a[12]*x[0]**2*x[1]+a[13]*x[0]**2*x[1]**2+ 
        a[14]*x[0]**3*x[1])
    
    def fquartic8(x):
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4+ 
           a[5]*x[1]+a[6]*x[1]**2+a[7]*x[1]**3+a[8]*x[1]**4+ 
           a[9]*x[0]*x[1]+a[10]*x[0]*x[1]**2+a[11]*x[0]*x[1]**3+a[12]*x[0]**2*x[1]+a[13]*x[0]**2*x[1]**2+ 
           a[14]*x[0]**3*x[1]+
           a[15]*x[2]+a[16]*x[2]**2+a[17]*x[2]**3+a[18]*x[2]**4+
           a[19]*x[0]*x[2]+a[20]*x[0]*x[2]**2+a[21]*x[0]*x[2]**3+a[22]*x[0]**2*x[2]+a[23]*x[0]**2*x[2]**2+
           a[24]*x[0]**3*x[2]+
           a[25]*x[1]*x[2]+a[26]*x[1]*x[2]**2+a[27]*x[1]*x[2]**3+a[28]*x[1]**2*x[2]+a[29]*x[1]**2*x[2]**2+
           a[30]*x[1]**3*x[2]+
           a[31]*x[0]*x[1]*x[2]+a[32]*x[0]**2*x[1]*x[2]+a[33]*x[0]*x[1]**2*x[2]+a[34]*x[0]*x[1]*x[2]**2
           )

    # Gradients:
    def fquartic_der_1(x):
        return np.array([a[1]+2.0*a[2]*x[0]+3.0*a[3]*x[0]**2+4.0*a[4]*x[0]**3])
        
    def fquartic_der_4(x):
        return np.array([a[1]+2.0*a[2]*x[0]+3.0*a[3]*x[0]**2+4.0*a[4]*x[0]**3+ 
        a[9]*x[1]+a[10]*x[1]**2+a[11]*x[1]**3+2.0*a[12]*x[0]*x[1]+2.0*a[13]*x[0]*x[1]**2+ 
        3.0*a[14]*x[0]**2*x[1],           
        a[5]+2.0*a[6]*x[1]+3.0*a[7]*x[1]**2+4.0*a[8]*x[1]**3+ 
        a[9]*x[0]+2.0*a[10]*x[0]*x[1]+3.0*a[11]*x[0]*x[1]**2+a[12]*x[0]**2+2.0*a[13]*x[0]**2*x[1]+ 
        a[14]*x[0]**3])
        
    def fquartic_der_8(x):
        return np.array([
        a[1]+2.0*a[2]*x[0]+3.0*a[3]*x[0]**2+4.0*a[4]*x[0]**3+ 
        a[9]*x[1]+a[10]*x[1]**2+a[11]*x[1]**3+2.0*a[12]*x[0]*x[1]+2.0*a[13]*x[0]*x[1]**2+ 
        3.0*a[14]*x[0]**2*x[1]+
        a[19]*x[2]+a[20]*x[2]**2+a[21]*x[2]**3+2.0*a[22]*x[0]*x[2]+2.0*a[23]*x[0]*x[2]**2+
        3.0*a[24]*x[0]**2*x[2]+
        a[31]*x[1]*x[2]+2.0*a[32]*x[0]*x[1]*x[2]+a[33]*x[1]**2*x[2]+a[34]*x[1]*x[2]**2,
        
        a[5]+2.0*a[6]*x[1]+3.0*a[7]*x[1]**2+4.0*a[8]*x[1]**3+ 
        a[9]*x[0]+2.0*a[10]*x[0]*x[1]+3.0*a[11]*x[0]*x[1]**2+a[12]*x[0]**2+2.0*a[13]*x[0]**2*x[1]+ 
        a[14]*x[0]**3+
        a[25]*x[2]+a[26]*x[2]**2+a[27]*x[2]**3+2.0*a[28]*x[1]*x[2]+2.0*a[29]*x[1]*x[2]**2+
        3.0*a[30]*x[1]**2*x[2]+
        a[31]*x[0]*x[2]+a[32]*x[0]**2*x[2]+2.0*a[33]*x[0]*x[1]*x[2]+a[34]*x[0]*x[2]**2, 
        
        a[15]+2.0*a[16]*x[2]+3.0*a[17]*x[2]**2+4.0*a[18]*x[2]**3+
        a[19]*x[0]+2.0*a[20]*x[0]*x[2]+3.0*a[21]*x[0]*x[2]**2+a[22]*x[0]**2+2.0*a[23]*x[0]**2*x[2]+
        a[24]*x[0]**3+
        a[25]*x[1]+2.0*a[26]*x[1]*x[2]+3.0*a[27]*x[1]*x[2]**2+a[28]*x[1]**2+2.0*a[29]*x[1]**2*x[2]+
        a[30]*x[1]**3+
        a[31]*x[0]*x[1]+a[32]*x[0]**2*x[1]+a[33]*x[0]*x[1]**2+2.0*a[34]*x[0]*x[1]*x[2]
        ])

    # Hessians:
    def fquartic_der2_1(x):
        return np.array([2.0*a[2]+6.0*a[3]*x[0]+12.0*a[4]*x[0]**2])
        
    def fquartic_der2_4(x):
        return np.array([
    
        [  2.0*a[2]+6.0*a[3]*x[0]+12.0*a[4]*x[0]**2+2.0*a[12]*x[1]+2.0*a[13]*x[1]**2+6.0*a[14]*x[0]*x[1], 
        
           a[9]+2.0*a[10]*x[1]+3.0*a[11]*x[1]**2+2.0*a[12]*x[0]+4.0*a[13]*x[0]*x[1]+3.0*a[14]*x[0]**2       ],     
           
        [  a[9]+2.0*a[10]*x[1]+3.0*a[11]*x[1]**2+2.0*a[12]*x[0]+4.0*a[13]*x[0]*x[1]+3.0*a[14]*x[0]**2, 
        
           2.0*a[6]+6.0*a[7]*x[1]+12.0*a[8]*x[1]**2+2.0*a[10]*x[0]+6.0*a[11]*x[0]*x[1]**1+2.0*a[13]*x[0]**2 ]
           
        ])
        
    def fquartic_der2_8(x):
        return np.array([
        
        [  2.0*a[2]+6.0*a[3]*x[0]+12.0*a[4]*x[0]**2+2.0*a[12]*x[1]+2.0*a[13]*x[1]**2+6.0*a[14]*x[0]*x[1]+
           +2.0*a[22]*x[2]+2.0*a[23]*x[2]**2+6.0*a[24]*x[0]*x[2]+2.0*a[32]*x[1]*x[2],
           
           a[9]+2.0*a[10]*x[1]+3.0*a[11]*x[1]**2+2.0*a[12]*x[0]+4.0*a[13]*x[0]*x[1]+3.0*a[14]*x[0]**2+
           a[31]*x[2]+2.0*a[32]*x[0]*x[2]+2.0*a[33]*x[1]*x[2]+a[34]*x[2]**2,       
        
           a[19]+2.0*a[20]*x[2]+3.0*a[21]*x[2]**2+2.0*a[22]*x[0]+4.0*a[23]*x[0]*x[2]+3.0*a[24]*x[0]**2+
           a[31]*x[1]+2.0*a[32]*x[0]*x[1]+a[33]*x[1]**2+a[34]*x[1]*x[2]                                     ], 
 
        [  a[9]+2.0*a[10]*x[1]+3.0*a[11]*x[1]**2+2.0*a[12]*x[0]+4.0*a[13]*x[0]*x[1]+3.0*a[14]*x[0]**2+
           a[31]*x[2]+2.0*a[32]*x[0]*x[2]+2.0*a[33]*x[1]*x[2]+a[34]*x[2]**2, 
           
           2.0*a[6]+6.0*a[7]*x[1]+12.0*a[8]*x[1]**2+2.0*a[10]*x[0]+6.0*a[11]*x[0]*x[1]+2.0*a[13]*x[0]**2+ 
           2.0*a[28]*x[2]+2.0*a[29]*x[2]**2+6.0*a[30]*x[1]*x[2]+2.0*a[33]*x[0]*x[2],       
        
           a[25]+2.0*a[26]*x[2]+3.0*a[27]*x[2]**2+2.0*a[28]*x[1]+4.0*a[29]*x[1]*x[2]+3.0*a[30]*x[1]**2+
           a[31]*x[0]+a[32]*x[0]**2+2.0*a[33]*x[0]*x[1]+2.0*a[34]*x[0]*x[2]                                  ],
          
        
        [  a[19]+2.0*a[20]*x[2]+3.0*a[21]*x[2]**2+2.0*a[22]*x[0]+4.0*a[23]*x[0]*x[2]+3.0*a[24]*x[0]**2+
           a[31]*x[1]+2.0*a[32]*x[0]*x[1]+a[33]*x[1]**2+2.0*a[34]*x[1]*x[2], 
           
           a[25]+2.0*a[26]*x[2]+3.0*a[27]*x[2]**2+2.0*a[28]*x[1]+4.0*a[29]*x[1]*x[2]+3.0*a[30]*x[1]**2+
           a[31]*x[0]+a[32]*x[0]**2+2.0*a[33]*x[0]*x[1]+2.0*a[34]*x[0]*x[2],       
        
           +2.0*a[16]+6.0*a[17]*x[2]+12.0*a[18]*x[2]**2+2.0*a[20]*x[0]+6.0*a[21]*x[0]*x[2]+2.0*a[23]*x[0]**2+
           2.0*a[26]*x[1]+6.0*a[27]*x[1]*x[2]+2.0*a[29]*x[1]**2+2.0*a[34]*x[0]*x[1]                          ],
        
        ])
    
    # First set the initial set, if given
    if guess!=None: 
        x_in = contract_vector(guess,ibrav)  
    else:
        x_in = contract_vector([0.,0.,0.,0.,0.,0.],ibrav)

    #  Find the minimun using minimize from scipy.optimize with the gradient
    #  whatever algorithm with the gradient is more stable and I recommend using it 
    if (ibrav==1 or ibrav==2 or ibrav==3):      # cubic systems (a,a,a)
        res_quartic = minimize(fquartic1,x_in,method=method, jac=fquartic_der_1, hess=fquartic_der2_1,  options=minoptions)
    elif (ibrav==4 or ibrav==6 or ibrav==7):    # hexagonal or tetragonal systems (a,a,c)
        res_quartic = minimize(fquartic4,x_in,method=method, jac=fquartic_der_4, hess=fquartic_der2_4,  options=minoptions)
    elif (ibrav>=8 and ibrav<=11):              # orthorombic systems (a,b,c)
        res_quartic = minimize(fquartic8,x_in,method=method, jac=fquartic_der_8, hess=fquartic_der2_8,  options=minoptions)  
    
    if (not res_quartic.success):
        print ("WARNING! Problems in extremum finding: ",res_quartic.message)
  
    return expand_vector(res_quartic.x,ibrav)


################################################################################
#
def calculate_fitted_points_anis(celldmsx,nmesh, fittype="quadratic",ibrav=4,a=None):
    """
    Calculates a denser mesh of Efitted(celldmsx) points for plotting. nmesh = (nx,ny,nz)
    gives the dimensions of the mesh. 
    """
    na = nmesh[0]
    nb = nmesh[1]
    nc = nmesh[2]
    
    if (fittype=="quadratic"):
        fit_fun=fquadratic
    else:
        fit_fun=fquartic
    
    if (ibrav==4):      # hex
        astep = (celldmsx[:,0].max()-celldmsx[:,0].min())/na
        cstep = (celldmsx[:,2].max()-celldmsx[:,2].min())/nc
        celldmsxdense = np.zeros((na*nc,6))
        Edense = np.zeros(na*nc)
        for i in range(0,na):
            for j in range(0,nc):
                index = i*na+j
                celldmsxdense[index,0] = celldmsx[:,0].min() + astep*i 
                celldmsxdense[index,2] = celldmsx[:,2].min() + cstep*j 
                Edense[index] = fit_fun(celldmsxdense[index,:],a)
        
    return celldmsxdense, Edense