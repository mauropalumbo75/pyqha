#encoding: UTF-8

from constants import RY_KBAR
from math import pow
import numpy as np
from scipy.optimize import curve_fit

################################################################################
# Murnaghan EOS functions 
#
# This one is in the format for the fit 
def E_MurnV(V,a0,a1,a2,a3):
    res=np.zeros(len(V))
    for i in range(0,len(V)):
        res[i]=a0 - a2*a1/(a3-1.0) + V[i]*a2/a3*( pow(a1/V[i],a3)/(a3-1.0)+1.0 )
    return res

# Other functions
def E_Murn(V,a):
    return a[0] - a[2]*a[1]/(a[3]-1.0) + V*a[2]/a[3]*( pow(a[1]/V,a[3])/(a[3]-1.0)+1.0 )

def P_Murn(V,a):
    return a[2]/a[3]*(pow(a[1]/V,a[3])-1.0)

def H_Murn(V,a):
    return E_Murn(V,a)+P_Murn(V,a)*V


################################################################################
# Print the data and the fitted results 
# ylabel can be "E", "Fvib", etc.
def print_data(x,y,a,chi,ylabel="E"):    
    print ("# Murnaghan EOS \t\t chi squared= {:.10e}".format(chi))
    print ("# E0= {:.10e}".format(a[1])+"\t V0= {:.10e}".format(a[1])+"\t B0= {:.10e}".format(a[2]*RY_KBAR)
    +"\t dB0/dV= {:.10e}".format(a[3]))
    print (80*"#")
    print ("# V","\t\t\t",ylabel,"\t\t\t",ylabel+"fit","\t\t\t",ylabel+"-"+ylabel+"fit")
    for i in range(0,len(y)):
        print ("{:.10e}".format(x[i]),"\t", "{:.10e}".format(y[i])+
        "\t {:.10e}".format(E_Murn(x[i],a))+
        "\t {:.10e}".format(y[i]-E_Murn(x[i],a))+
        "\t {:.10e}".format(P_Murn(x[i],a)*RY_KBAR))

################################################################################
# Write the data and the fitted results as in the Murnaghan EOS 
# ylabel can be "E", "Fvib", etc.
def write_Etotfitted(filename,x,y,a,chi,ylabel="E"): 
    fout=open(filename, "w")
    fout.write("# Murnaghan EOS \t\t chi squared= {:.10e}".format(chi))
    fout.write("# E0= {:.10e}".format(a[1])+"\t V0= {:.10e}".format(a[1])+"\t B0= {:.10e}".format(a[2]*RY_KBAR)
    +"\t dB0/dV= {:.10e}".format(a[3]))
    fout.write(80*"#")
    fout.write("# V"+"\t\t\t"+ylabel+"\t\t\t"+ylabel+"fit"+"\t\t\t"+ylabel+"-"+ylabel+"fit")
    for i in range(0,len(y)):
        fout.write("{:.10e}".format(x[i])+"\t"+"{:.10e}".format(y[i])+
        "\t {:.10e}".format(E_Murn(x[i],a))+
        "\t {:.10e}".format(y[i]-E_Murn(x[i],a))+
        "\t {:.10e}".format(P_Murn(x[i],a)*RY_KBAR))

################################################################################
# Calculated a denser mesh of E(V) points for plotting...
#
def calculate_fitted_points(V,a):
    Vstep = (V[len(V)-1]-V[0])/1000
    Vdense = np.zeros(1000)
    Edensefitted = np.zeros(1000)
    for i in range(0,1000):
        Vdense[i] = V[0] + Vstep*i 
        Edensefitted[i] = E_Murn(Vdense[i],a)
        
    return Vdense, Edensefitted


################################################################################
# This is the function for fitting with a EOS as a function of volume only
#
# The input variable x is an 1D array of volumes, y are the corresponding 
# energies
# 
def fit_Murn(V,E):
    
    # reasonable initial guesses for EOS parameters
    a0=E[len(E)/2]
    a1=V[len(V)/2]
    a2=500/RY_KBAR
    a3=5.0
    
    # Create the auxiliary A numpy matrix for the fitting coefficients (quadratic polynomial)
    a, pcov = curve_fit(E_MurnV, V, E, p0=[a0,a1,a2,a3])
    
    chi = 0
    for i in range(0,len(V)):
        chi += (E[i]-E_Murn(V[i],a))**2
    
    return a, pcov, chi



def compute_beta(minT):
    grad=np.gradient(np.array(minT))  # numerical derivatives with numpy
    betaT = np.array(grad)  # grad contains the derivatives with respect to T
                                # also convert to np.array format    
    return betaT/minT