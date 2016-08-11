# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import minimize


################################################################################

def expand_vector(x, ibrav=4):
    """
    Utility function: expands a vector x, len(x)<6, into a 6-dim vector according
    to ibrav
    Note: not all ibrav are implemented yet
    """
    if (len(x)>6 or x==None):
        return None
    elif (len(x)==6):
        return x
    else:
        if ibrav==1 or ibrav==2 or ibrav==3:
            xnew = np.array([x[0],0.0,0.0,0.0,0.0,0.0])
        elif ibrav==4:
            xnew = np.array([x[0],0.0,x[1],0.0,0.0,0.0])
        else:
            pass
        return xnew
            
    
################################################################################

def contract_vector(x, ibrav=4):
    """
    Utility function: contract a vector x, len(x)=6, into a x-dim vector (x<6) 
    according to ibrav
    Note: not all ibrav are implemented yet
    """
    if (len(x)!=6):
        return None
    else:
        if ibrav==1 or ibrav==2 or ibrav==3:
            xnew = np.array([x[0]])
        elif ibrav==4:
            xnew = np.array([x[0],x[2]])
        else:
            pass
        return xnew

################################################################################
           
def fquadratic(x,a,ibrav=4):
    """
    Implemented polynomials for fitting and miminizing

    only ibrav=4 and the most general case are implemented for now     
    """
    if ibrav==4:
        #return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[2]+a[4]*x[2]**2+a[5]*x[0]*x[2]
    else:
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
        +a[6]*x[2]+a[7]*x[2]**2+a[8]*x[0]*x[2]+a[9]*x[1]*x[2]
        +a[10]*x[3]+a[11]*x[3]**2+a[12]*x[0]*x[3]+a[13]*x[1]*x[3]+a[14]*x[2]*x[3]
        +a[15]*x[4]+a[16]*x[4]**2+a[17]*x[0]*x[4]+a[18]*x[1]*x[4]+a[19]*x[2]*x[4]+a[20]*x[3]*x[4]
        +a[21]*x[5]+a[22]*x[5]**2+a[23]*x[0]*x[5]+a[24]*x[1]*x[5]+a[25]*x[2]*x[5]+a[26]*x[3]*x[5]+a[26]*x[4]*x[5]

def fquadratic_der(x,a,ibrav=4):
    if ibrav==4:
        return np.array([
        a[1]+2.0*a[2]*x[0]+a[5]*x[2],
        0.0,
        a[3]+2.0*a[4]*x[2]+a[5]*x[0],
        0.0,
        0.0,
        0.0
        ])
    else:
        return np.array([ 
        a[1]+2.0*a[2]*x[0]+a[5]*x[1]+a[8]*x[2]+a[12]*x[3]+a[17]*x[4]+a[23]*x[5],        
        a[3]+2.0*a[4]*x[1]+a[5]*x[1]+a[9]*x[2]+a[13]*x[3]+a[18]*x[4]+a[24]*x[5], 
        a[6]+2.0*a[7]*x[2]+a[8]*x[0]+a[9]*x[1]+a[14]*x[3]+a[19]*x[4]+a[25]*x[5], 
        a[10]+2.0*a[11]*x[3]+a[12]*x[0]+a[13]*x[1]+a[14]*x[2]+a[20]*x[4]+a[26]*x[5],
        a[15]+2.0*a[16]*x[4]+a[17]*x[0]+a[18]*x[1]+a[19]*x[2]+a[20]*x[3]+a[26]*x[5], 
        a[21]+2.0*a[22]*x[5]+a[23]*x[0]+a[24]*x[1]+a[25]*x[2]+a[26]*x[3]+a[26]*x[4]
        ])  

# only ibrav=4 is implemented for now  
def fquartic(x,a,ibrav=4):
    if ibrav==4:
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4+ 
           a[5]*x[2]+a[6]*x[2]**2+a[7]*x[2]**3+a[8]*x[2]**4+ 
           a[9]*x[0]*x[2]+a[10]*x[0]*x[2]**2+a[11]*x[0]*x[2]**3+a[12]*x[0]**2*x[2]+a[13]*x[0]**2*x[2]**2+ 
           a[14]*x[0]**3*x[2])
        #return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4+ 
        #   a[5]*x[1]+a[6]*x[1]**2+a[7]*x[1]**3+a[8]*x[1]**4+ 
        #   a[9]*x[0]*x[1]+a[10]*x[0]*x[1]**2+a[11]*x[0]*x[1]**3+a[12]*x[0]**2*x[1]+a[13]*x[0]**2*x[1]**2+ 
        #   a[14]*x[0]**3*x[1])
    else:
        return None
    
def fquartic_der(x,a,ibrav=4):
    if ibrav==4:
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
    else:
        return None
    
    
    
################################################################################

def find_min(a,ibrav,type,guess=None):
    """
    An auxiliary function for handling the minimum search    
    """
    if type=="quadratic":
        xmin=find_min_quadratic(a,ibrav,guess)
        fmin=fquadratic(xmin,a)
    elif type=="quartic":
        xmin=find_min_quartic(a,ibrav,guess)
        fmin=fquartic(xmin,a)
    else:
        return None
    print ("Minimun "+type+": ",xmin,"\tEnergy at the minimum: "+"{:.20e}".format(fmin)+"\n") 

    return xmin, fmin

################################################################################

def find_min_quadratic(a,ibrav=4,guess=None):
    """
    This is the function for finding the minimum of the quadratic polynomial
    """
    
    # Define here the polynomial and gradient functions... not ideal, better with global variables?
    def fquadratic4(x):
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
    
    def fquadratic(x):
        return a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[1]+a[4]*x[1]**2+a[5]*x[0]*x[1]
        +a[6]*x[2]+a[7]*x[2]**2+a[8]*x[0]*x[2]+a[9]*x[1]*x[2]
        +a[10]*x[3]+a[11]*x[3]**2+a[12]*x[0]*x[3]+a[13]*x[1]*x[3]+a[14]*x[2]*x[3]
        +a[15]*x[4]+a[16]*x[4]**2+a[17]*x[0]*x[4]+a[18]*x[1]*x[4]+a[19]*x[2]*x[4]+a[20]*x[3]*x[4]
        +a[21]*x[5]+a[22]*x[5]**2+a[23]*x[0]*x[5]+a[24]*x[1]*x[5]+a[25]*x[2]*x[5]+a[26]*x[3]*x[5]+a[26]*x[4]*x[5]

    def fquadratic_der4(x):
        return np.array([a[1]+2.0*a[2]*x[0]+a[5]*x[1], a[3]+2.0*a[4]*x[1]+a[5]*x[0]])
    
    def fquadratic_der(x):
        return np.array( 
        a[1]+2.0*a[2]*x[0]+a[5]*x[1]+a[8]*x[2]+a[12]*x[3]+a[17]*x[4]+a[23]*x[5],        
        a[3]+2.0*a[4]*x[1]+a[5]*x[1]+a[9]*x[2]+a[13]*x[3]+a[18]*x[4]+a[24]*x[5], 
        a[6]+2.0*a[7]*x[2]+a[8]*x[0]+a[9]*x[1]+a[14]*x[3]+a[19]*x[4]+a[25]*x[5], 
        a[10]+2.0*a[11]*x[3]+a[12]*x[0]+a[13]*x[1]+a[14]*x[2]+a[20]*x[4]+a[26]*x[5],
        a[15]+2.0*a[16]*x[4]+a[17]*x[0]+a[18]*x[1]+a[19]*x[2]+a[20]*x[3]+a[26]*x[5], 
        a[21]+2.0*a[22]*x[5]+a[23]*x[0]+a[24]*x[1]+a[25]*x[2]+a[26]*x[3]+a[26]*x[4])            

    # First set the initial set, if given
    if guess!=None: 
        x_in = contract_vector(guess,ibrav)
    else:
        x_in = contract_vector([0.,0.,0.,0.,0.,0.],ibrav)

    #  Find the minimun using minimize from scipy.optimize with the gradient
    #  whatever algorithm with the gradient is more stable and I recommend using it 
    if ibrav==4:
        #res_quadratic = minimize(fquadratic4,x_in,jac=fquadratic_der4)
        res_quadratic = minimize(fquadratic4,x_in,method="BFGS", jac=fquadratic_der4,  options={'gtol': 1e-5})
    else:
        res_quadratic = minimize(fquadratic,x_in, jac=fquadratic_der)
        #res_quadratic = minimize(fquadratic,x_in, method="BFGS", jac=fquadratic_der,  options={'gtol': 1e-12}) 
  
    if (not res_quadratic.success):
        print ("WARNING! Problems in extremum finding: ",res_quadratic.message)
  
    return expand_vector(res_quadratic.x,ibrav)


################################################################################

def find_min_quartic(a,ibrav=4,guess=None):
    """
    This is the function for finding the minimum of the quartic polynomial
    """
    # Define here the polynomial and gradient functions... not ideal, better with global variables?
    def fquartic4(x):
        return (a[0]+a[1]*x[0]+a[2]*x[0]**2+a[3]*x[0]**3+a[4]*x[0]**4+ 
        a[5]*x[1]+a[6]*x[1]**2+a[7]*x[1]**3+a[8]*x[1]**4+ 
        a[9]*x[0]*x[1]+a[10]*x[0]*x[1]**2+a[11]*x[0]*x[1]**3+a[12]*x[0]**2*x[1]+a[13]*x[0]**2*x[1]**2+ 
        a[14]*x[0]**3*x[1])

    def fquartic_der4(x):
        return np.array([a[1]+2.0*a[2]*x[0]+3.0*a[3]*x[0]**2+4.0*a[4]*x[0]**3+ 
        a[9]*x[1]+a[10]*x[1]**2+a[11]*x[1]**3+2.0*a[12]*x[0]*x[1]+2.0*a[13]*x[0]*x[1]**2+ 
        3.0*a[14]*x[0]**2*x[1],           
        a[5]+2.0*a[6]*x[1]+3.0*a[7]*x[1]**2+4.0*a[8]*x[1]**3+ 
        a[9]*x[0]+2.0*a[10]*x[0]*x[1]+3.0*a[11]*x[0]*x[1]**2+a[12]*x[0]**2+2.0*a[13]*x[0]**2*x[1]+ 
        a[14]*x[0]**3])
        
    def fquartic_der2_4(x):
        return np.array([[2.0*a[2]+6.0*a[3]*x[0]+12.0*a[4]*x[0]**2+2.0*a[12]*x[1]
        +2.0*a[13]*x[1]**2+6.0*a[14]*x[0]*x[1],     
        a[9]+2.0*a[10]*x[1]+3.0*a[11]*x[1]**2+2.0*a[12]*x[0]+4.0*a[13]*x[0]*x[1]+ 
        3.0*a[14]*x[0]**2],        
        [a[9]+2.0*a[10]*x[1]+3.0*a[11]*x[1]**2+2.0*a[12]*x[0]+4.0*a[13]*x[0]*x[1]+ 
        3.0*a[14]*x[0]**2,      
        2.0*a[6]+6.0*a[7]*x[1]+12.0*a[8]*x[1]**2+ 
        +2.0*a[10]*x[0]+6.0*a[11]*x[0]*x[1]**1+2.0*a[13]*x[0]**2]])
    
    # First set the initial set, if given
    if guess!=None: 
        x_in = contract_vector(guess,ibrav)  
    else:
        x_in = contract_vector([0.,0.,0.,0.,0.,0.],ibrav)

    #  Find the minimun using minimize from scipy.optimize with the gradient
    #  whatever algorithm with the gradient is more stable and I recommend using it 
    if ibrav==4:
        #res_quartic = minimize(fquartic4,x_in,jac=fquartic_der4)
        #res_quartic = minimize(fquartic4,x_in,method="BFGS", jac=fquartic_der4,  options={'gtol': 1e-7})
        #print (res_quartic.hess_inv)
        #print (fquartic_der2_4(res_quartic.x))
        #res_quartic = minimize(fquartic4,x_in,method="CG", jac=fquartic_der4,  options={'gtol': 1e-8})
        res_quartic = minimize(fquartic4,x_in,method="Newton-CG", jac=fquartic_der4, hess=fquartic_der2_4,  options={'xtol': 1e-8})
        #res_quartic = minimize(fquartic4,x_in,method="dogleg", jac=fquartic_der4, hess=fquartic_der2_4,  options={'gtol': 1e-8})
#    else:
#        res_quartic = minimize(fquartic,x_in,jac=fquartic_der)        
  
    if (not res_quartic.success):
        print ("WARNING! Problems in extremum finding: ",res_quartic.message)
  
    return expand_vector(res_quartic.x,ibrav)