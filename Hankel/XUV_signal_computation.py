"""
This module xxx.
"""
import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn
# import mynumerics as mn
import matplotlib.pyplot as plt

import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index





## Parameters

# # gas
# gas_type = 'He'
# XUV_table_type_absorption = 'Henke' # {Henke, NIST}    
# XUV_table_type_dispersion = 'NIST'

# # laser  
# omegaSI = mn.ConvertPhoton(800e-9, 'lambdaSI', 'omegaSI') 
# Horder_init = 95


# control
# l1 = 1.5e-3
# xi = 1.0
# pressure = 42e-3
# ionisation_ratio = 0.03
# zeta = 1.8e-6

## Calculation

N_atm = XUV_index.N_atm_default #2.7e19 * 1e6
# N_atm = 2.6867774e25


def Phi_2pi_decider(Phi, tol = 8.*np.finfo(float).eps):
    return ((np.imag(Phi) <= tol ) and
            (
            ((np.real(Phi) % (2.0*np.pi)) <= 2.*tol)
            or
            ( abs( (np.real(Phi) % (2.0*np.pi)) - 2.0*np.pi) <= 2.*tol)
            ))


def compute_S1_abs(pressure, zeta, ionisation_ratio, l1, Horder, parameters, include_absorption = True):
    
    gas_type = parameters['gas_type']
    XUV_table_type_absorption = parameters['XUV_table_type_absorption']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']
    
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
    
    # # polarisability XUV (including absorption)
    # def polarisability_XUV(Horder_foo):
    #     f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(Horder_foo*omegaSI, 'omegaSI', 'eV'))
    #     nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(Horder_foo*omegaSI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    #     susc_XUV_atm = nXUV_atm**2 - 1
    #     pol_XUV = susc_XUV_atm/N_atm
    #     return pol_XUV
    
    # def f2_funct(E):
    #     return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)

    # def beta_factor_atm(Horder_foo):
    #     omegaXUV    = Horder_foo*omegaSI
    #     f2_value    = f2_funct(mn.ConvertPhoton(omegaXUV, 'omegaSI', 'eV'))
    #     lambdaXUV    = mn.ConvertPhoton(omegaXUV, 'omegaSI', 'lambdaSI')
    #     beta_factor = N_atm*units.r_electron_classical * \
    #                   ((lambdaXUV**2)*f2_value/(2.0*np.pi))
    #     return beta_factor

    # def L_abs(Horder_foo):
    #     omegaXUV    = Horder_foo*omegaSI
    #     f2_value    = f2_funct(mn.ConvertPhoton(omegaXUV, 'omegaSI', 'eV'))
    #     lambdaXUV    = mn.ConvertPhoton(omegaXUV, 'omegaSI', 'lambdaSI')
    #     return 1.0 / (2.0 * pressure * N_atm * units.r_electron_classical * lambdaXUV * f2_value) 
    
    # # polarisability IR
    # susc_IR_atm = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    # # print('susIR anal:', susc_IR_atm)
    # polarisability_IR = susc_IR_atm/N_atm
    
    # # print('nXUV anal:', np.sqrt(1.+ pressure*N_atm*polarisability_XUV(Horder)))
    # # print('nIR anal:', np.sqrt(1.+ pressure*N_atm*polarisability_IR)) 
    
    # susc_IR_atm = IR_index.susc_atm(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    # susc_XUV_atm = XUV_index.susc_atm(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion) 
    
    polarisability_IR = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    polarisability_XUV = XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion) 
    
    
    # other parameters
    k0 = omegaSI /units.c_light
    
    ## delta k1
    delta_k1 = Horder * k0 * (0.5*pressure*N_atm*( (polarisability_IR - polarisability_XUV) - ionisation_ratio*plasma_constant) - zeta )
    
    # not working for vectorised inputs
    # if (abs(delta_k1)<=np.finfo(float).eps): L_coh = float('inf')
    # else: L_coh = np.abs(np.pi/delta_k1)
    
    L_coh = np.abs(np.pi/delta_k1) 
    
    # add absorption
    if include_absorption:
        beta_factor_atm = XUV_index.beta_factor_atm(Horder*omegaSI, gas_type + '_' + XUV_table_type_absorption)
        delta_k1 = delta_k1 + 1j*Horder*k0 * pressure * beta_factor_atm
    
    ## field from 1 segment
    # deal with listed inputs, quite tedius, just need to decide which variable makes it a list
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
    
    if include_absorption:
        return S1, delta_k1, L_coh, XUV_index.L_abs(Horder*omegaSI, pressure, gas_type + '_' + XUV_table_type_absorption)#L_abs(Horder)
    else:
        return S1, delta_k1, L_coh







def compute_Phi_abs(pressure, zeta, l1, xi, ionisation_ratio, Horder, parameters, include_absorption = True):
    
    gas_type = parameters['gas_type']
    XUV_table_type_absorption = parameters['XUV_table_type_absorption']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']
    
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
    
    # # polarisability XUV
    # def polarisability_XUV(Horder_foo):
    #     f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(Horder_foo*omegaSI, 'omegaSI', 'eV'))
    #     nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(Horder_foo*omegaSI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    #     susc_XUV_atm = nXUV_atm**2 - 1
    #     pol_XUV = susc_XUV_atm/N_atm
    #     return pol_XUV
 
    # def f2_funct(E):
    #     return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)

    # def beta_factor_atm(Horder_foo):
    #     omegaXUV    = Horder_foo*omegaSI
    #     f2_value    = f2_funct(mn.ConvertPhoton(omegaXUV, 'omegaSI', 'eV'))
    #     lambdaXUV    = mn.ConvertPhoton(omegaXUV, 'omegaSI', 'lambdaSI')
    #     beta_factor = N_atm*units.r_electron_classical * \
    #                   ((lambdaXUV**2)*f2_value/(2.0*np.pi))
    #     return beta_factor


    # # polarisability IR
    # susc_IR_atm = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    # polarisability_IR = susc_IR_atm/N_atm


    polarisability_IR = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    polarisability_XUV = XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion) 
    
    # other parameters
    k0 = omegaSI /units.c_light
    # plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
    # Aq = 1.0 # normalised XUV gain   
    
    Phi = Horder*l1*k0*(
                    0.5*pressure*N_atm*( (polarisability_IR - polarisability_XUV) - ionisation_ratio*plasma_constant) -
                    zeta * (1.0 + xi)
                    )
    
    # add absorption
    if include_absorption:
        beta_factor_atm = XUV_index.beta_factor_atm(Horder*omegaSI, gas_type + '_' + XUV_table_type_absorption)
        Phi = Phi + 1j*Horder*l1*k0 * pressure * beta_factor_atm
    
    
    return Phi


def compute_sum_abs(pressure, zeta, l1, xi, ionisation_ratio, Horder, m_max, parameters, include_absorption = True):
    
    Phi = compute_Phi_abs(pressure, zeta, l1, xi, ionisation_ratio, Horder, parameters, include_absorption=include_absorption)
    
   
    # singular case when we have almost perfect phase matching
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



def compute_chain_abs(pressure, zeta, l1, xi, ionisation_ratio, Horder, m_max, parameters, include_absorption = True):
      
    S1 = compute_S1_abs(pressure, zeta, ionisation_ratio, l1, Horder, parameters, include_absorption=include_absorption)
    chain = compute_sum_abs(pressure, zeta, l1, xi, ionisation_ratio, Horder, m_max, parameters, include_absorption=include_absorption)
    signal = S1[0]*chain[0]
    
    
    
    # singular case of perfect phase-matching in one segment
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
 
        
 
    
    # singular case of perfect phase-matching
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
    
    signal2 = (pressure * parameters['Aq'])**2 * abs_S1_2 * abs_chain_2 
    
    return signal, signal2



def zeta_calc(delta_phi, pressure, l1, xi, ionisation_ratio, Horder, parameters):
    
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']    
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
    
    # # polarisability XUV
    # def polarisability_XUV(Horder_foo):
    #     f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(Horder_foo*omegaSI, 'omegaSI', 'eV'))
    #     nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(Horder_foo*omegaSI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    #     susc_XUV_atm = nXUV_atm**2 - 1
    #     pol_XUV = susc_XUV_atm/N_atm
    #     return pol_XUV

    # # polarisability IR
    # susc_IR_atm = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    # polarisability_IR = susc_IR_atm/N_atm    


    polarisability_IR = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    polarisability_XUV = XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion) 
    
    # other parameters
    k0 = omegaSI /units.c_light
    
    
    zeta = (1.0/(1.0+xi)) * (
            0.5 * pressure * N_atm * ( (polarisability_IR - polarisability_XUV) - ionisation_ratio*plasma_constant) -
            2.0 * delta_phi/ (k0*l1))
    
    return zeta



def eta_opt(Horder, parameters):
    
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']   

    # susc_IR = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))



    # # N_atm = 2.6867774e25
    # def susc_XUV(omega_in_SI):
    #     f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(omega_in_SI, 'omegaSI', 'eV'))
    #     nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega_in_SI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    #     return nXUV_atm**2 - 1
    #     # susc_XUV = nXUV_atm**2 - 1
        
    # def delta_susc(omega_in_SI):
    #     return (susc_IR - susc_XUV(omega_in_SI)) /N_atm
        
    delta_polarisability = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI')) - \
                           XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion)
    

    return (omegaSI**2) *units.eps0*units.elmass*delta_polarisability/(units.elcharge**2)



def zeta_single_segment(pressure, Horder, ionisation_ratio, parameters):
    
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']   

    # susc_IR = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)



    # # N_atm = 2.6867774e25
    # def susc_XUV(omega_in_SI):
    #     f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(omega_in_SI, 'omegaSI', 'eV'))
    #     nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega_in_SI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    #     return nXUV_atm**2 - 1
    #     # susc_XUV = nXUV_atm**2 - 1
        
    # def delta_susc(omega_in_SI):
    #     return (susc_IR - susc_XUV(omega_in_SI)) /N_atm

    delta_polarisability = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI')) - \
                           XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion)        
    

    zeta =  0.5 * pressure * N_atm * ( delta_polarisability - ionisation_ratio*plasma_constant) 
    
    return zeta

def xi_calc_pm(delta_phi, pressure, l1, zeta, ionisation_ratio, Horder, parameters):
    
    gas_type = parameters['gas_type']
    XUV_table_type_dispersion = parameters['XUV_table_type_dispersion']
    omegaSI = parameters['omegaSI']    
    plasma_constant = units.elcharge**2 / (units.eps0 * units.elmass * omegaSI**2)
    
    # # polarisability IR
    # # susc_IR_atm = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    # # polarisability_IR = susc_IR_atm/N_atm    
    # susc_IR = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
    
    # # N_atm = 2.6867774e25
    # def susc_XUV(omega_in_SI):
    #     f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(omega_in_SI, 'omegaSI', 'eV'))
    #     nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega_in_SI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    #     return nXUV_atm**2 - 1
    #     # susc_XUV = nXUV_atm**2 - 1
        
    # def delta_susc(omega_in_SI):
    #     return (susc_IR - susc_XUV(omega_in_SI)) /N_atm
    

    delta_polarisability = IR_index.polarisability(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI')) - \
                           XUV_index.polarisability(Horder*omegaSI, gas_type+'_'+XUV_table_type_dispersion)        
    
    # other parameters
    k0 = omegaSI /units.c_light
    
    
    xi = (1/zeta) * (0.5 * pressure * N_atm * (delta_polarisability - ionisation_ratio*plasma_constant) -
                     2.0 * delta_phi / (k0*l1)
                     ) - 1.0
    
    
            # (1.0/(1.0+xi)) * (
            # 0.5 * pressure * N_atm * ( (polarisability_IR - polarisability_XUV(Horder)) - ionisation_ratio*plasma_constant) -
            # 2.0 * delta_phi/ (k0*l1))
    
    return xi

xi2r = lambda xi : 1.0/(1.0 + xi)
r2xi = lambda r : (1.0-r)/r