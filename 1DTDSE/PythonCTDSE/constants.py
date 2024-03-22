"""
Constants
=========

This module contains physical constants for PythonCTDSE.
"""

import numpy as np

### Physical constants
Ip_HeV = 27.21138602 
"""Ionisation potential"""
hbar = 1.054571800e-34
"""Planck constant"""
alpha_fine = 1/137.035999139
"""Fine structure constant"""
c_light = 299792458.
"""Speed of light in vacuum"""
elcharge = 1.602176565e-19
"""Electron charge"""
elmass = 9.10938356e-31
"""Electron mass"""
mu0 = 4.0*np.pi*1e-7
"""Magnetic permeability"""
eps0 = 1.0/(mu0*c_light*c_light)
"""Vacuum permitivity"""
r_Bohr = 4.0*np.pi*eps0*hbar*hbar/(elmass*elcharge*elcharge)
"""Bohr radius"""
TIMEau = (elmass*r_Bohr*r_Bohr)/hbar
"""Atomic unit of time"""
EFIELDau = hbar*hbar/(elmass*r_Bohr*r_Bohr*r_Bohr*elcharge)
"""Atomic unit of field intensity"""
k_Boltz = 1.38064852e-23
"""Boltzmann constant"""
absolute_zero = -273.15
"""Absolute zero temperature"""
torr2SI = 101325./760.
"""Torr to SI conversion factor"""