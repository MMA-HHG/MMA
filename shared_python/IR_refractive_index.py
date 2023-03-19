"""
This module provides suscpetibilities and polarisabilities for noble gases
(He, Ne, Ar, Kr, Xe) given by the analytic formulae from
https://royalsocietypublishing.org/doi/epdf/10.1098/rspa.1960.0237

The functions are:
    
    getf, getf1, getf2
    
Next, there are other functions to access directly polarisabilities, susceptibilities,
absorption lengths, ... (see their descriptions) The functions are:
    
    dispersion_function, beta_factor_atm, L_abs, susc_atm, polarisability
    
Note: The reference is provided in 'atmospheric pressure' characterised by the
Loschmidt constant (https://en.wikipedia.org/wiki/Loschmidt_constant),
the default value is N_atm = 2.7e25 m-3.

-------
Jan Vabek   
ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2021 - 2022)   
ELI ERIC (2023)  
"""
susc_Ar = lambda x: (5.547e-4)*(1.0 + (5.15e5)/(x**2) + (4.19e11)/(x**4) + (4.09e17)/(x**6) + (4.32e23)/(x**8))
susc_Kr = lambda x: (8.377e-4)*(1.0 + (6.7e5)/(x**2) + (8.84e11)/(x**4) + (1.49e18)/(x**6) + (2.74e24)/(x**8) + (5.1e30)/(x**10))
susc_He = lambda x: (6.927e-5)*(1.0 + (2.24e5)/(x**2) + (5.94e10)/(x**4) + (1.72e16)/(x**6))
susc_Ne = lambda x: (1.335e-4)*(1.0 + (2.24e5)/(x**2) + (8.09e10)/(x**4) + (3.56e16)/(x**6))
susc_Xe = lambda x: (1.366e-3)*(1.0 + (9.02e5)/(x**2) + (1.81e12)/(x**4) + (4.89e18)/(x**6) + (1.45e25)/(x**8) + (4.34e31)/(x**10))

susc_funct = {
    'Ar': susc_Ar,
    'Kr': susc_Kr,
    'He': susc_Kr,
    'Ne': susc_Kr,
    'Xe': susc_Kr
}

def getsusc(g,lambd):
    """
    Returns susceptibility for atmospheric pressure.

    Parameters
    ----------
    g : string
        available gases: 'He', 'Ne', 'Ar', 'Kr', 'Xe'  
    lambd : laser wavelength [m]

    Returns
    -------
    susceptibility
    """
    return susc_funct[g](1e10*lambd)

susc_atm = getsusc

def getpol(g,lambd,N_atm=2.7e19*1e6):
    """
    Returns polarisability.
    
    Parameters
    ----------
    g : string
        available gases: 'He', 'Ne', 'Ar', 'Kr', 'Xe'  
    lambd : laser wavelength [m]   
    N_atm : float, optional
        Number density [m-3] at atmospheric pressure. The default is 2.7e15.
        (https://en.wikipedia.org/wiki/Loschmidt_constant)

    Returns
    -------
    polarisability
    
    Notes
    -------
    The original paper tabulates susceptibilities at 'standard pressure'.

    """
    return getsusc(g,lambd)/N_atm

polarisability = getpol