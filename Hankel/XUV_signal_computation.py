"""
This module xxx.



parameters = {'XUV_table_type_dispersion': 'Henke' or 'NIST',
              'XUV_table_type_absorption': 'Henke' or 'NIST',
              'gas_type': available gases 'He', 'Ne', 'Ar', 'Kr', 'Xe'  
              'omegaSI': fundamental laser frequency [rad/s],
              'Aq' : Amplitude of the harmonic response}

-------
Jan Vabek - ELI ERIC (2023)
"""
import numpy as np
import units
import mynumerics as mn
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index



N_atm = XUV_index.N_atm_default # gas number density (atmospheric pressure)

def Phi_2pi_decider(Phi, tol = 8.*np.finfo(float).eps):
    """
    Auxilliary function to decide phase giving singular absorption-free &
    perfectly phase-matched generation. (Some analytic expressions are singular
    and need special treatment).

    Parameters
    ----------
    Phi : complex phase [rad]
    tol : float, optional
        The tolerance for considering phase giving absorption-free & perfectly
        phase matched. The default is 8.*np.finfo(float).eps (doubled for the
        real part mod 2*pi).

    Returns
    -------
    boolean
    """
    return ((np.imag(Phi) <= tol ) and
                (
                ((np.real(Phi) % (2.0*np.pi)) <= 2.*tol)
              or
                ( abs( (np.real(Phi) % (2.0*np.pi)) - 2.0*np.pi) <= 2.*tol)
            ))


def single_period_S1(pressure, zeta, ionisation_ratio, l1, Horder, parameters, include_absorption = True):
    """
    Compute the signal from a gas cell. The legth of the cell is l1, it includes
    dispersion, absorption, geometrical (linear) phase and plasma due to ionisation.

    Parameters
    ----------
    pressure : pressure [bar]   
    zeta : geometrical phase factor [-]   
    ionisation_ratio : ionisation degree [-]   
    l1 : length of one gas-medium [m]   
    Horder : harmonic order [-]   
    parameters : dict (see documentation of the module)     
    include_absorption : boolean, optional.  The default is True.  

    Returns
    -------
    S1 : output signal (complexified field)  
    delta_k1 : wavenumer mismatch (possibly complex)  
    L_coh : coherence length [m]  
    L_abs : absorption length [m]
        np.inf for absorption-free medium
    """
    
    gas_type = parameters['gas_type']
    XUV_table_type_absorption = parameters['XUV_table_type_absorption']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']
    
    k0 = omegaSI /units.c_light
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)    
    
    polarisability_IR = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    polarisability_XUV = XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion)    
    
    delta_k1 = Horder * k0 * (0.5*pressure*N_atm*( (polarisability_IR - polarisability_XUV) - ionisation_ratio*plasma_constant) - zeta )
    L_coh = np.abs(np.pi/delta_k1) 
       
    # add absorption
    if include_absorption:
        beta_factor_atm = XUV_index.beta_factor_atm(Horder*omegaSI, gas_type + '_' + XUV_table_type_absorption)
        delta_k1 = delta_k1 + 1j*Horder*k0 * pressure * beta_factor_atm
    
    ## Here we compute the generated field after the distance l1, it is given by
    #         Aq* 1j * (np.exp(1j * delta_k1 * l1)-1.0) / delta_k1)
    # all the tedious construction is here to deal with singular absorption-free 
    # perfecly-phase-matched cases for vectorised inputs.
    
    if hasattr(delta_k1 * l1, "__len__"):
        phase = delta_k1 * l1; S1 = []
        l1_list = hasattr(l1, "__len__")
        pressure_list = hasattr(pressure, "__len__")
        for k1 in range(len(phase)):
            if Phi_2pi_decider(phase[k1]): # singular case of perfect phase-matching
              if pressure_list:
                S1.append(pressure[k1] * parameters['Aq'] * 1j * (1j*l1))
              elif l1_list:
                S1.append(pressure * parameters['Aq'] * 1j * (1j*l1[k1]))
              else:
                S1.append(pressure * parameters['Aq'] * 1j * (1j*l1))
            else:
              if pressure_list:
                S1.append(pressure[k1] * parameters['Aq'] * 1j * (np.exp(1j * delta_k1 * l1)-1.0) / delta_k1)
              elif l1_list:
                S1.append(pressure * parameters['Aq'] * 1j * (np.exp(1j * delta_k1 * l1[k1])-1.0) / delta_k1)
              else:
                S1.append(pressure * parameters['Aq'] * 1j * (np.exp(1j * delta_k1[k1] * l1)-1.0) / delta_k1)
        S1 = np.asarray(S1)
    else:            
        if Phi_2pi_decider(delta_k1 * l1): # singular case of perfect phase-matching
            S1 = pressure * parameters['Aq'] * 1j * (1j*l1)
        else:
            S1 = pressure * parameters['Aq'] * 1j * (np.exp(1j * delta_k1 * l1)-1.0) / delta_k1
    
    # compute L_abs
    if include_absorption: L_abs = XUV_index.L_abs(Horder*omegaSI, pressure, gas_type + '_' + XUV_table_type_absorption)
    else: L_abs = np.inf
        
    return S1, delta_k1, L_coh, L_abs


def compute_Phi(pressure, zeta, l1, xi, ionisation_ratio, Horder, parameters, include_absorption = True):
    """
    The complex phase characterising the geometry (see Section III).

    Parameters
    ----------
    pressure : pressure [bar]   
    zeta : geometrical phase factor [-]   
    l1 : length of one gas-medium [m]   
    xi : ratio l2/l1 (vacuum/gas) [-]    
    ionisation_ratio : ionisation degree [-]   
    Horder : harmonic order [-]   
    parameters : dict (see documentation of the module)  
    include_absorption : boolean, optional. The default is True.  
    
    Returns
    -------
    Phi : The (complex) phase characterising the chain of the media.
        See Eq. (12c).

    """
    
    gas_type = parameters['gas_type']
    XUV_table_type_absorption = parameters['XUV_table_type_absorption']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']
    
    k0 = omegaSI /units.c_light
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
    
    polarisability_IR = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    polarisability_XUV = XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion) 
 
    
    Phi = Horder*l1*k0*(
                    0.5*pressure*N_atm*( (polarisability_IR - polarisability_XUV) - ionisation_ratio*plasma_constant) -
                    zeta * (1.0 + xi)
                    )
    
    # add absorption
    if include_absorption:
        beta_factor_atm = XUV_index.beta_factor_atm(Horder*omegaSI, gas_type + '_' + XUV_table_type_absorption)
        Phi = Phi + 1j*Horder*l1*k0 * pressure * beta_factor_atm
    
    
    return Phi


def periodic_medium_sum(pressure, zeta, l1, xi, ionisation_ratio, Horder, m_max, parameters, include_absorption = True):
    """
    Compute the sum (not including the single emitter) that modulates the signal
    from the chain of peridically repating gas-media.
    
    Parameters
    ----------
    pressure : pressure [bar]   
    zeta : geometrical phase factor [-]   
    l1 : length of one gas-medium [m]   
    xi : ratio l2/l1 (vacuum/gas) [-]    
    ionisation_ratio : ionisation degree [-]   
    Horder : harmonic order [-]   
    m_max : the number of periods gas-vacuum (integer)  
    parameters : dict (see documentation of the module)  
    include_absorption : boolean, optional. The default is True.  

    Returns
    -------
    output signal (complexified)  
        The signal from a chain of media.
    Phi : phase characterising the periodic medium
        See function 'compute_Phi'.
        
    """
    
    
    Phi = compute_Phi(pressure, zeta, l1, xi, ionisation_ratio, Horder, parameters, include_absorption=include_absorption)
    
   
    # Deal with singular vectorised cases, see the comment inside 'single_period_S1'
    if hasattr(Phi, "__len__"):
        signal = []
        for k1 in range(len(Phi)):
            if Phi_2pi_decider(Phi[k1]):
                signal.append(m_max+1)
            else:
                signal.append(
                    (np.exp(1j*Phi[k1]*(m_max+1)) - 1.0)/ (np.exp(1j*Phi[k1]) - 1.0)
                    )
        signal = np.asarray(signal)
        return signal, Phi
    else:
        if Phi_2pi_decider(Phi): 
            return (m_max+1), Phi
        else:
            return (np.exp(1j*Phi*(m_max+1)) - 1.0)/ (np.exp(1j*Phi) - 1.0), Phi


def periodic_medium_signal(pressure, zeta, l1, xi, ionisation_ratio, Horder, m_max, parameters, include_absorption = True):
    """
    The main computation routine provide the total signal from the pariodic medium.
    The total signal is given as the product of the signal from a single period
    modulated by the chain of media. In other words, it's the coherent sum of all
    the periods. See Eqs. (8) and (14).

    Parameters
    ----------
    pressure : pressure [bar]   
    
    zeta : geometrical phase factor [-]   
    l1 : length of one gas-medium [m]   
    xi : ratio l2/l1 (vacuum/gas) [-]    
    ionisation_ratio : ionisation degree [-]   
    Horder : harmonic order [-]   
    m_max : the number of periods gas-vacuum (integer)  
    parameters : dict (see documentation of the module)  
    include_absorption : boolean, optional. The default is True.  

    Returns
    -------
    signal : The total (complexified) signal from the periodic medium.
    signal2 : |signal|^2
        
    Note
    -------
    The complexified signal and its |·|^2 are computed independetly by Eqs. (7) and (14).
    """
      
    S1 = single_period_S1(pressure, zeta, ionisation_ratio, l1, Horder, parameters, include_absorption=include_absorption)
    chain = periodic_medium_sum(pressure, zeta, l1, xi, ionisation_ratio, Horder, m_max, parameters, include_absorption=include_absorption)
    
    signal = S1[0]*chain[0] # This is already the required signal    
    
    
    ## In the second part, we analytically compute the |·|^2 of the signal using
    # analytic expression.
    
    # Deal with singular vectorised cases for a single segmant applied for |·|^2.
    # See the comment inside 'single_period_S1'
    phase = l1*S1[1]
    if hasattr(phase, "__len__"):
        l1_list = hasattr(l1, "__len__")
        abs_S1_2 = []
        k1r = np.real(S1[1])
        k1i = np.imag(S1[1])
        for k1 in range(len(phase)):
            if Phi_2pi_decider(phase[k1]):
              if l1_list:
                abs_S1_2.append(l1[k1]**2)
              else:
                abs_S1_2.append(l1**2)
            else:
              if l1_list:
                abs_S1_2.append(np.exp(-k1i*l1[k1]) * ( (np.sinh(k1i*l1[k1]))**2 + (np.sin(k1r*l1[k1]))**2) / (k1r**2 + k1i**2))
              else:
                abs_S1_2.append(np.exp(-k1i[k1]*l1) * ( (np.sinh(k1i[k1]*l1))**2 + (np.sin(k1r[k1]*l1))**2) / (k1r[k1]**2 + k1i[k1]**2))
                
        abs_S1_2 = np.asarray(abs_S1_2)
    else:
        if Phi_2pi_decider(l1*S1[1]):
            abs_S1_2 = l1**2
        else:
            k1r = np.real(S1[1])
            k1i = np.imag(S1[1])
            abs_S1_2 = np.exp(-k1i*l1) * ( (np.sinh(k1i*l1))**2 + (np.sin(k1r*l1))**2) / (k1r**2 + k1i**2)
            
    
    # Deal with singular vectorised cases for the chain applied for |·|^2.
    # See the comment inside 'single_period_S1'
    phase = chain[1]
    if hasattr(phase, "__len__"):
        abs_chain_2 = []
        for k1 in range(len(phase)):
            if Phi_2pi_decider(phase[k1]): 
                abs_chain_2.append(m_max**2)
            else:
                Phir = np.real(phase[k1])
                Phii = np.imag(phase[k1])
                abs_chain_2.append(
                                   np.exp(-(m_max-1) * Phii) *(((np.sinh(0.5*m_max*Phii))**2 + (np.sin(0.5*m_max*Phir))**2)/
                                                          ((np.sinh(0.5*Phii))**2 + (np.sin(0.5*Phir))**2))
                                   )
        abs_chain_2 = np.asarray(abs_chain_2)                            
    else:
        if Phi_2pi_decider(chain[1]): 
            abs_chain_2 = m_max**2
        else:
            Phir = np.real(chain[1])
            Phii = np.imag(chain[1])
            abs_chain_2 = np.exp(-(m_max-1) * Phii) *(((np.sinh(0.5*m_max*Phii))**2 + (np.sin(0.5*m_max*Phir))**2)/
                                                  ((np.sinh(0.5*Phii))**2 + (np.sin(0.5*Phir))**2))
    
    # Computte |·|^2 of the signal.
    signal2 = (pressure * parameters['Aq'])**2 * abs_S1_2 * abs_chain_2 
    
    return signal, signal2



## Functions to find optimising parameters of the scheme.

def eta_opt(Horder, parameters):
    """
    Compute the optimal ionisation degree for phase matching.
    (It compensates the dispersion of both IR and XUV.)

    Parameters
    ----------
    Horder : harmonic order [-]  
    parameters : dict (see documentation of the module)  

    Returns
    -------
    optimal ionisation degree [-]
    """
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']
        
    delta_polarisability = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI')) - \
                           XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion)

    return (omegaSI**2) *units.eps0*units.elmass*delta_polarisability/(units.elcharge**2)


def zeta_single_segment_pm(pressure, Horder, ionisation_ratio, parameters):
    """
    Compute the geometrical phase 'zeta' for perfect phase matching within
    the gas.

    Parameters
    ----------
    pressure : pressure [bar]      
    Horder : harmonic order [-]  
    ionisation_ratio : ionisation degree [-]   
    parameters : dict (see documentation of the module)

    Returns
    -------
    zeta [-]
    """
    
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)

    delta_polarisability = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI')) - \
                           XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion)   

    zeta =  0.5 * pressure * N_atm * ( delta_polarisability - ionisation_ratio*plasma_constant) 
    
    return zeta


def xi_chain_pm(delta_phi, pressure, l1, zeta, ionisation_ratio, Horder, parameters):
    """
    Compute the stride characterised by 'xi' to ensure the phase jump by
    'delta_phi' within one elementary segment. The optimal value to select
    'Horder' is delta_phi = pi/(Horder*n); n = -1, -2, -3, ...

    Parameters
    ----------
    delta_phi : requred phase jump [rad]  
    pressure : pressure [bar]  
    l1 : length of one gas-medium [m]  
    zeta : geometrical phase factor [-]   
    ionisation_ratio : ionisation degree [-]   
    Horder : harmonic order [-]   
    parameters : dict (see documentation of the module)  

    Returns
    -------
    xi [-]
    """
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']    
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
        

    delta_polarisability = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI')) - \
                           XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion)        
    
    # other parameters
    k0 = omegaSI /units.c_light
    
    
    xi = (1/zeta) * (0.5 * pressure * N_atm * (delta_polarisability - ionisation_ratio*plasma_constant) -
                     2.0 * delta_phi / (k0*l1)
                     ) - 1.0
    
    return xi


def zeta_chain_pm(delta_phi, pressure, l1, xi, ionisation_ratio, Horder, parameters):
    """
    Compute the geometrical phase 'zeta' to ensure the phase jump by
    'delta_phi' within one elementary segment. The optimal value to select
    'Horder' is delta_phi = pi/(Horder*n); n = -1, -2, -3, ...

    Parameters
    ----------
    delta_phi : requred phase jump [rad]    
    pressure : pressure [bar]    
    l1 : length of one gas-medium [m]    
    xi : ratio l2/l1 (vacuum/gas) [-]    
    ionisation_ratio : ionisation degree [-]    
    Horder : harmonic order [-]    
    parameters : dict (see documentation of the module)

    Returns
    -------
    zeta [-]
    """
    
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']  
    k0 = omegaSI /units.c_light
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
    
    polarisability_IR = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    polarisability_XUV = XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion) 

    zeta = (1.0/(1.0+xi)) * (
            0.5 * pressure * N_atm * ( (polarisability_IR - polarisability_XUV) - ionisation_ratio*plasma_constant) -
            2.0 * delta_phi/ (k0*l1))
    
    return zeta


## Transformations of the stride characterised by 'xi' or 'r'
xi2r = lambda xi : 1.0/(1.0 + xi)
r2xi = lambda r : (1.0-r)/r


## back compatibility, kept for instant
zeta_calc = zeta_chain_pm 
xi_calc_pm = xi_chain_pm
zeta_single_segment = zeta_single_segment_pm

compute_chain_abs = periodic_medium_signal
compute_sum_abs = periodic_medium_sum
compute_S1_abs = single_period_S1