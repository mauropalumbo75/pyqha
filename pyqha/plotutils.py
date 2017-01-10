#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
A collection of wrappers for the *matplotlib* functions. 

.. Note::
  All functions return a *matplotlib* which can be modified by the user.
"""


try:
    wx
    from matplotlib import use
    use('WXAgg')
except:
    pass
import matplotlib.pyplot as plt

    
import numpy as np


def simple_plot_xy(x,y,xlabel="",ylabel=""):
    """
    This function generates a simple xy plot with matplotlib.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    ax.plot(x, y, 'r')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    
    return fig

   
def multiple_plot_xy(x,y,xlabel="",ylabel="",labels=""):
    """
    This function generates a simple xy plot with matplotlib overlapping several
    lines as in the matrix y. y second index refers to a line in the plot, the first 
    index is for the array to be plotted.
    """
    
    if (len(y[0,:])>7):
        print ("Too many data on y axis!")
        return
    
    colors = ['k','r','b','g','c','m','y']
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    if (labels==""):     
        try:    # try if there are multiple data on x axis
            for i in range(0,len(y[0,:])):
                ax.plot(x[:,i], y[:,i], colors[i])
        except: # if not use a single x axis
            for i in range(0,len(y[0,:])):
                ax.plot(x, y[:,i], colors[i])
    else:
        try:    # try if there are multiple data on x axis
            for i in range(0,len(y[0,:])):
                ax.plot(x[:,i], y[:,i], colors[i],label=labels[i])
        except: # if not use a single x axis
            for i in range(0,len(y[0,:])):
                ax.plot(x, y[:,i], colors[i],label=labels[i])       
        ax.legend() 
        
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()
    
    return fig
    

def plot_EV(V,E,a=None,labely="Etot"):
    """
    This function plots with matplotlib E(V) data and if a is given it also plot
    the fitted results
    """
    
    from .eos import calculate_fitted_points
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    ax.plot(V, E, 'o', label=labely+" data", markersize=10)
    if (a!=None):
        Vdense, Edensefitted = calculate_fitted_points(V,a)
        ax.plot(Vdense, Edensefitted, 'r', label='Fitted EOS')
    ax.legend()
    ax.set_xlabel('V (a.u.^3)')
    ax.set_ylabel('E (Ry)')
    plt.show()
 
    return fig
    
    
def plot_Etot(celldmsx,Ex,n,nmesh=(50,50,50),fittype="quadratic",ibrav=4,a=None):
    """
    This function makes a 3D plot with matplotlib Ex(celldmsx) data and if a is given it also plot
    the fitted results. The plot type depends on ibrav.
    """
    
    from mpl_toolkits.mplot3d import axes3d
    from matplotlib import cm
    from .minutils import calculate_fitted_points_anis
    
    if (Ex==None) and (a==None):
        return
    
    fig = plt.figure() 
    if (ibrav==4):  # hex case

        ax = fig.gca(projection='3d')       
        
        if (Ex!=None):
            na=n[0]
            nc=n[2]
            X = np.zeros((na,nc))
            Y = np.zeros((na,nc))
            Z = np.zeros((na,nc))

            for j in range(0,nc):
                for i in range(0,na): 
                    index = j*na+i
                    X[i,j] = celldmsx[index,0]
                    Y[i,j] = celldmsx[index,2]
                    Z[i,j] = Ex[index]
                    #print (index,X[i,j],Y[i,j],Z[i,j])
                    
            ax.set_xlim(X.min(),X.max())
            ax.set_ylim(Y.min(),Y.max())
            ax.set_zlim(Z.min(),Z.max())
            ax.scatter(X,Y,Z,c='r',marker='o')
    
  
        if (a!=None):
            celldmsxdense, Edensefitted = calculate_fitted_points_anis(celldmsx,nmesh,fittype,ibrav,a)
            Xd = np.zeros((nmesh[0],nmesh[2]))
            Yd = np.zeros((nmesh[0],nmesh[2]))
            Zd = np.zeros((nmesh[0],nmesh[2]))
            for i in range(0,nmesh[0]):
                for j in range(0,nmesh[2]): 
                    index = i*nmesh[0]+j
                    Xd[i,j] = celldmsxdense[index,0]
                    Yd[i,j] = celldmsxdense[index,2]
                    Zd[i,j] = Edensefitted[index]
    
            ax.set_xlim(Xd.min(),Xd.max())
            ax.set_ylim(Yd.min(),Yd.max())
            ax.set_zlim(Zd.min(),Zd.max())
            ax.plot_surface(Xd, Yd, Zd, rstride=1, cstride=1, alpha=0.3)
            cset = ax.contour(Xd, Yd, Zd, zdir='z', offset=Zd.min(), cmap=cm.coolwarm)
       
        ax.set_xlabel("a (a.u.)")
        ax.set_ylabel("c (a.u.)")
        ax.set_zlabel("Etot (Ry)")
        plt.show()
    
    return fig


def plot_Etot_contour(celldmsx,nmesh=(50,50,50),fittype="quadratic",ibrav=4,a=None):
    """
    This function makes a countour plot with matplotlib of Ex(celldmsx) fitted results. 
    The plot type depends on ibrav.
    """

    from .minutils import calculate_fitted_points_anis
    
    if a==None:
        return
    
    fig = plt.figure() 
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    if (ibrav==4):
        celldmsxdense, Edensefitted = calculate_fitted_points_anis(celldmsx,nmesh,fittype,ibrav,a)
        Xd = np.zeros((nmesh[0],nmesh[2]))
        Yd = np.zeros((nmesh[0],nmesh[2]))
        Zd = np.zeros((nmesh[0],nmesh[2]))
        for i in range(0,nmesh[0]):
            for j in range(0,nmesh[2]): 
                index = i*nmesh[0]+j
                Xd[i,j] = celldmsxdense[index,0]
                Yd[i,j] = celldmsxdense[index,2]
                Zd[i,j] = Edensefitted[index]

        CS = ax.contour(Xd, Yd, Zd)
        plt.clabel(CS, inline=1, fontsize=10)
        
        CS.ax.set_xlabel("a (a.u.)")
        CS.ax.set_ylabel("c (a.u.)")
        plt.show()
    
    return fig
