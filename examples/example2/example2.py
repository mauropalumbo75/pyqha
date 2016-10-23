#!/usr/bin/env python
#encoding: UTF-8

"""
This is a simple example of fitting the total energy E as a function of the cell
parameters (a,c) of an hexagonal cell with a quartic polynomial.
"""
    
if __name__ == "__main__":

    from pyqha import fitEtot, plot_Etot, plot_Etot_contour

    fin = "./Etot.dat"  	# contains the input energies
    
    # fits the energies and returns the coeffients a0 and the chi squared chia0
    # the fit is done with a quartic polynomial
    celldmsx, Ex, a0, chia0, mincelldms, fmin = fitEtot(fin,fittype="quartic",guess=[5.12374914,0.0,8.19314311,0.0,0.0,0.0])
    
    # 3D plot only with fitted energy
    plot_Etot(celldmsx,Ex=None,n=(5,0,5),nmesh=(50,0,50),fittype="quartic",ibrav=4,a=a0)
    # 3D plot fitted energy and points
    plot_Etot(celldmsx,Ex,n=(5,0,5),nmesh=(50,0,50),fittype="quartic",ibrav=4,a=a0)
    # 2D contour plot with fitted energy 
    plot_Etot_contour(celldmsx,nmesh=(50,0,50),fittype="quartic",ibrav=4,a=a0)
    
