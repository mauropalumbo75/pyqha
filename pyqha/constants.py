#encoding: UTF-8
# Copyright (C) 2016 Mauro Palumbo
# This file is distributed under the terms of the # MIT License. 
# See the file `License' in the root directory of the present distribution.

"""
Some useful standard constants for conversions and calculations.
""" 

PI     = 3.14159265358979323846 
TPI    = 2.0 * PI
FPI    = 4.0 * PI

C_SI             = 2.99792458E+8    # m sec^-1
H_PLANCK_SI      = 6.62606896E-34   # J s
K_BOLTZMANN_SI   = 1.3806504E-23    # J K^-1 
HARTREE_SI       = 4.35974394E-18   # J
BOHR_RADIUS_SI   = 0.52917720859E-10 # m
EV               = 1.602176565e-19 # electronvolt, J
NA               = 6.02214129e23  # Avogadro constant, mol-1
RYDBERG_SI       = HARTREE_SI/2.0   # J
AU_SEC           = H_PLANCK_SI/TPI/HARTREE_SI
AU_PS            = AU_SEC * 1.0E+12
AU_TERAHERTZ  = AU_PS
AU_GPA           = HARTREE_SI / BOHR_RADIUS_SI ** 3 / 1.0E+9 
K_BOLTZMANN_RY   = K_BOLTZMANN_SI / RYDBERG_SI
RY_TO_THZ = 1.0 / AU_TERAHERTZ / FPI
RY_TO_GHZ = RY_TO_THZ*1000.0
RY_TO_CMM1 = 1.0E+10 * RY_TO_THZ / C_SI
RY_KBAR          = 10.0 * AU_GPA / 2.0

KB1=1.0/K_BOLTZMANN_RY/RY_TO_CMM1 # inverse Boltzmann constant in cm^{-1}/K

KB_TO_EV         = K_BOLTZMANN_SI/EV  # Boltzmann constant in eV/K
EV_TO_J          = EV 
INVERSEM_TO_EV   = C_SI*H_PLANCK_SI             # m-1 to J