"""
This module creates the interpolation functions from tabulated scattering factors
(in the XUV range) for noble gases stored in the external 'XUV_refractive_index_tables.h5'
h5-file (see detaild in the h5-file constructor). These functions are:  
    
    getf, getf1, getf2
    
Next, there are other functions to access directly polarisabilities, susceptibilities,
absorption lengths, ... (see their descriptions).  
The functions are:
    
    dispersion_function, beta_factor_atm, L_abs, susc_atm, polarisability
    
Note: The reference is provided in 'atmospheric pressure' characterised by the
Loschmidt constant (https://en.wikipedia.org/wiki/Loschmidt_constant),
the default value is N_atm = 2.7e25 m-3.

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

N_atm_default = 2.7e19*1e6

def dispersion_function(omega, pressure, gas, n_IR=1., N_atm=N_atm_default):
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
    N_atm : scalar, optional
        gas number density for atmospheric pressure (see the module description).
        The default is N_atm_default.

    Returns
    -------
    (1/phase_velocity_IR) - (1/phase_velocity_XUV)
    """
    f1_value = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    nXUV     = 1.0 - pressure*N_atm*units.r_electron_classical * ((lambdaSI**2)*f1_value/(2.0*np.pi))           
    phase_velocity_XUV  = units.c_light / nXUV    
    phase_velocity_IR = units.c_light / n_IR
    return ((1./phase_velocity_IR) - (1./phase_velocity_XUV))


def beta_factor_atm(omega, gas, N_atm=N_atm_default):
    """
    It returns the imaginary part 'beta' of the refractive index  
                 n = n0 + 1j*beta   
    at atmospúheric pressure using tabulated scattering factors. Usage for
    arbitrary pressure:
                 n = n0(p) + 1j*(p/p_atm)*beta_atm . 

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]  
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    N_atm : scalar, optional
        gas number density for atmospheric pressure (see the module description).
        The default is N_atm_default.

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
    beta_factor = N_atm*units.r_electron_classical * \
                  ((lambdaXUV**2)*f2_value/(2.0*np.pi))
    return beta_factor


def L_abs(omega, pressure, gas, N_atm=N_atm_default):
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
    N_atm : scalar, optional
        gas number density for atmospheric pressure (see the module description).
        The default is N_atm_default.

    Returns
    -------
    L_abs [m]
    """
    f2_value    = getf2(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaXUV   = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 1.0 / (2.0 * pressure * N_atm * units.r_electron_classical * lambdaXUV * f2_value) 


def susc_atm(omega, gas, N_atm=N_atm_default):
    """
    Returns susceptibility for atmospheric pressure.

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]   
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    N_atm : scalar, optional
        gas number density for atmospheric pressure (see the module description).
        The default is N_atm_default.

    Returns
    -------
    susceptibility
    """
    f1 = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    return nXUV_atm**2 - 1


def polarisability(omega, gas, N_atm=N_atm_default):
    """
    Returns polarisability.

    Parameters
    ----------
    omega : scalar
        The frequency of the incident field [rad/s]   
    gas : string
        The specifier of gas and used tables, it has the form {'He', 'Ne', 'Ar',
        'Kr', 'Xe'}+'_'+{'NIST','Henke'}. For example gas='Ar_NIST'.  
    N_atm : scalar, optional
        gas number density for atmospheric pressure (see the module description).
        The default is N_atm_default.

    Returns
    -------
    polarisability
    """
    f1 = getf1(gas,mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    susc_XUV_atm = nXUV_atm**2 - 1
    pol_XUV = susc_XUV_atm/N_atm
    return pol_XUV