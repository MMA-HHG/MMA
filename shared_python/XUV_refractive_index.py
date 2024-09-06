"""
This module creates the interpolation functions from tabulated scattering factors
(in the XUV range) for noble gases stored in the external 'XUV_refractive_index_tables.h5'
h5-file (the available tables are Henke: https://henke.lbl.gov/optical_constants/asf.html
and NIST: https://physics.nist.gov/PhysRefData/FFast/html/form.html). These functions are:  
    
    getf, getf1, getf2
    
Next, there are other functions to access directly polarisabilities, susceptibilities,
absorption lengths, ... (see their descriptions).  
The functions are:
    
    dispersion_function, beta_factor_ref, L_abs, susc_ref, polarisability

The module also provides the reference particle density 'N_ref_default'  for
p = 1 bar & T = 20 °C

Note: The reference particle density is given by the ideal-gas law. Assuming this law,
it is directly scalable to any pressure p temperature T by p[bar]*((273.15+20)/T[Kelvin])

-------
Jan Vabek   
ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2021 - 2022)   
ELI ERIC (2023)
"""

import numpy as np
from scipy import interpolate
import h5py
import os
import units
import mynumerics as mn

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# Load tabulated scattering factors stored in 'XUV_refractive_index_tables.h5'
# and create interpolating functions from them.
source_archive = os.path.join(THIS_DIR, 'XUV_refractive_index_tables.h5')
index_funct = {}
with h5py.File(source_archive, 'r') as SourceFile: # access option http://docs.h5py.org/en/stable/high/file.html#file
    gases = list(SourceFile.keys())
    index_table = {}
    print(gases)
    for gas in gases:
        local_table = {
            'Energy_f1': SourceFile[gas]['Energy_f1'][:],
            'Energy_f2': SourceFile[gas]['Energy_f2'][:],
            'f1': SourceFile[gas]['f1'][:],
            'f2': SourceFile[gas]['f2'][:]
        }
        index_table.update({gas: local_table})
        local_table = {
            'f1': interpolate.interp1d(SourceFile[gas]['Energy_f1'][:], SourceFile[gas]['f1'][:]),
            'f2': interpolate.interp1d(SourceFile[gas]['Energy_f2'][:], SourceFile[gas]['f2'][:])
        }
        index_funct.update({gas: local_table})


## FUNCTIONS PROVIDING THE SCATTERING FACTORS

def getf(g,E):
    """
    Returns tabulated scattering factors for a given XUV photon.

    Parameters
    ----------
    g : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.
    E : scalar
        the energy of the incident photon [eV]

    Returns
    -------
    (f1, f2): the values of the scattering factors

    """
    return index_funct[g]['f1'](E)[()], index_funct[g]['f2'](E)[()]

def getf1(g,E):
    """
    Returns tabulated scattering factor 'f1' for a given XUV photon.

    Parameters
    ----------
    g : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.
    E : scalar
        the energy of the incident photon [eV]

    Returns
    -------
    f1: The scattering factor

    """
    return index_funct[g]['f1'](E)[()]

def getf2(g,E):
    """
    Returns tabulated scattering factor 'f2' for a given XUV photon.

    Parameters
    ----------
    g : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.
    E : scalar
        the energy of the incident photon [eV]

    Returns
    -------
    f1: The scattering factor

    """
    return index_funct[g]['f2'](E)[()]


## VARIOUS FUNCTIONS TO PROVIDE DIRECTLY POLARISABILITIES, SUSCEPTIBILITIES, ...

# N_ref_default = 1e5/(units.Boltzmann_constant*(273.15+20.)) # reference gas number density (p = 1 bar & T = 20 °C)
N_ref_default = 2.7e25

def dispersion_function(omega, pressure, gas, n_IR=1., N_ref=N_ref_default):
    """
    Returns the part of the dephasing caused by the different phase velocities
    of the IR and XUV fields. The output quantity is  
              (1/phase_velocity_IR) - (1/phase_velocity_XUV)  
    phase_velocity_XUV is computed from tabulated scattering factors  
    phase_velocity_IR = c_light/n_IR.
    
    Possbile usage in linear medium (after the distance 'z') is  
              phase = omega*z*dispersion_function  

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]  
    pressure : pressure [bar]   
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    n_IR : scalar, optional
        The refractive index in the IR range. The default is 1 (i.e. vacuum progation).  
    N_ref : scalar, optional
        gas number density for  (see the module description).
        The default is N_ref_default.

    Returns
    -------
    (1/phase_velocity_IR) - (1/phase_velocity_XUV)
    """
    f1_value = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    # print(f1_value)
    if (-9998. > np.asarray(f1_value)).any(): raise ValueError('f1-scattering coefficient not listed for this frequency.')
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    nXUV     = 1.0 - pressure*N_ref*units.r_electron_classical * ((lambdaSI**2)*f1_value/(2.0*np.pi))           
    phase_velocity_XUV  = units.c_light / nXUV    
    phase_velocity_IR = units.c_light / n_IR
    return ((1./phase_velocity_IR) - (1./phase_velocity_XUV))


def beta_factor_ref(omega, gas, N_ref=N_ref_default):
    """
    It returns the imaginary poart 'beta' of the refractive index  
                 n = n0 + 1j*beta   
    using tabulated scattering factors

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]  
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    N_ref : scalar, optional
        gas number particle density
        The default is N_ref_default (p = 1 bar & T = 20 °C)

    Returns
    -------
    beta_factor : scalar
        
    Notes
    -------
    See Chapter 3.1, Eqs. (3.12) and (3.13) of 'D. Attwood; SOFT X-RAYS AND
    EXTREME ULTRAVIOLET RADIATION, Cambridge University Press, 1st Edition (2000)'
    """
    f2_value    = getf2(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaXUV    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    beta_factor = N_ref*units.r_electron_classical * \
                  ((lambdaXUV**2)*f2_value/(2.0*np.pi))
    return beta_factor


def L_abs(omega, pressure, gas, N_ref=N_ref_default):
    """
    Returns absorption length in XUV range.

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]   
    pressure : pressure [bar]   
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    N_ref : scalar, optional
        gas number particle density
        The default is N_ref_default (p = 1 bar & T = 20 °C)

    Returns
    -------
    L_abs [m]
    """
    f2_value    = getf2(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaXUV   = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 1.0 / (2.0 * pressure * N_ref * units.r_electron_classical * lambdaXUV * f2_value) 


def susc_ref(omega, gas, N_ref=N_ref_default):
    """
    Returns susceptibility for p = 1 bar & T = 20 °C.

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]   
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    N_ref : scalar, optional
        gas number particle density
        The default is N_ref_default (p = 1 bar & T = 20 °C)

    Returns
    -------
    susceptibility
    """
    f1 = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    nXUV_ref = 1.0 - N_ref*units.r_electron_classical*(mn.ConvertPhoton(omega,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    return nXUV_ref**2 - 1


def nXUV(omega,gas,pressure,complex=True, N_ref=N_ref_default):
    f1 = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    if complex: f2 = getf2(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    else:       f2 = 0.

    return 1.0 - pressure*N_ref*units.r_electron_classical*(mn.ConvertPhoton(omega,'omegaSI','lambdaSI')**2)*(f1+1j*f2)/(2.0*np.pi)


def polarisability(omega, gas, N_ref=N_ref_default):
    """
    Returns polarisability.

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]   
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    N_ref : scalar, optional
        gas number density for atmospheric pressure (see the module description).
        The default is N_ref_default.

    Returns
    -------
    polarisability
    """
    f1 = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    nXUV_ref = 1.0 - N_ref*units.r_electron_classical*(mn.ConvertPhoton(omega,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    susc_XUV_ref = nXUV_ref**2 - 1
    pol_XUV = susc_XUV_ref/N_ref
    return pol_XUV