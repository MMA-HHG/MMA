from scipy import special
from scipy import integrate
from scipy import interpolate
import numpy as np
import struct
import array
import os
import time
import warnings
# import ray
# import matlab.engine
# import string
import multiprocessing as mp
import math
# import joblib
# from mpi4py import MPI
# import oct2py
import shutil
import h5py
import sys
import units
import mynumerics as mn

import XUV_refractive_index as XUV_index

from Hankel_tools import get_propagation_pre_factor_function

import copy

def HankelTransform(ogrid, rgrid, FField, distance, rgrid_FF,
                    integrator = integrate.trapz,
                    near_field_factor = True,
                    pre_factor = 1.):
    """
    It computes Hankel transform with an optional near-field factor.
    

    Parameters
    ----------
    ogrid : array_like
        grid of FField in frequencies [SI]
    rgrid : array_like
        grid of FField in the radial coordinate [SI]
    FField : 2D array
        The source terms on ogrid and rgrid. (FField[omega,r])       
    distance : scalar
       The distance of the generating plane from the observational screen
    rgrid_FF : array_like
        The grid used to investigate the transformed field
    integrator : function handle, optional
        The function used for the integration. It's called by
        'integrator(integrand,rgrid)'. It shall be extended by a list of
        [args, kwargs].
        The default is integrate.trapz (from scipy).
    near_field_factor : logical, optional
        Include near field factor. The default is True.

    Returns
    -------
    FField_FF : 2D array
         The far-field spectra on ogrid and rgrid_FF

    """
    
    # print('pre-factor shape ', np.shape(pre_factor))
    # print('rgrid_int len ', len(rgrid))
    # print('ogrid     len ', len(ogrid))
    
    
    t_start = time.perf_counter()
    
    apply_radial_factor = (len(np.shape(pre_factor))==2)
    
    
    No = len(ogrid); Nr = len(rgrid); Nr_FF = len(rgrid_FF)
    FField_FF = np.empty((No,Nr_FF), dtype=np.cdouble)
    
    integrands = np.empty((No,Nr_FF,Nr), dtype=np.cdouble)
    rgrids = np.empty((No,Nr_FF,Nr), dtype=np.cdouble)
    
    integrand = np.empty((Nr), dtype=np.cdouble)
    for k1 in range(No):
        k_omega = ogrid[k1] / units.c_light # ogrid[k3] / units.c_light; ogrid[k1] * units.alpha_fine  # ??? units
        for k2 in range(Nr_FF):
            for k3 in range(Nr):
                if near_field_factor:
                    if apply_radial_factor:  radial_factor_local = pre_factor[k3,k1]
                    else:                    radial_factor_local = 1.
                    integrand[k3] = radial_factor_local *\
                                    np.exp(-1j * k_omega * (rgrid[k3] ** 2) / (2.0 * distance)) * rgrid[k3] *\
                                    FField[k1,k3] * special.jn(0, k_omega * rgrid[k3] * rgrid_FF[k2] / distance)
                else:
                    if apply_radial_factor:  radial_factor_local = pre_factor[k3,k1]
                    else:                    radial_factor_local = 1.
                    integrand[k3] = radial_factor_local *\
                                    rgrid[k3] * FField[k1,k3] * special.jn(0, k_omega * rgrid[k3] * rgrid_FF[k2] / distance)
                                    
            FField_FF[k1,k2] = integrator(integrand,rgrid)
            integrands[k1,k2,:] = integrand
            rgrids[k1,k2,:] = rgrid
            

    # if   (len(np.shape(pre_factor))==1):
    #     FField_FF *= np.outer(pre_factor,np.ones(FField_FF.shape[1]))
    # elif (len(np.shape(pre_factor))==0):
    #     FField_FF *= pre_factor
    
    # if (len(np.shape(pre_factor))==0):
    #     FField_FF *= pre_factor


    print('time spent only in the integrator ', time.perf_counter()-t_start)
    
    return FField_FF , [integrands, rgrids, FField_FF, copy.copy(FField_FF)]


def HankelTransform_long(target, # FSourceTerm(r,z,omega)
                         distance, rgrid_FF,
                         preset_gas = 'vacuum',
                         pressure = 1.,
                         absorption_tables = 'Henke',
                         include_absorption = True,
                         dispersion_tables = 'Henke',
                         include_dispersion = True,
                         effective_IR_refrective_index = 1.,
                         integrator_Hankel = integrate.trapz,
                         integrator_longitudinal = 'trapezoidal',
                         near_field_factor = True,
                         store_cummulative_result = False,
                         frequencies_to_trace_maxima = None,
                         store_entry_plane_transform = False
                         ):
    # """
    # It computes XUV propagation using a sum of Hankel transforms along the medium.

    # Parameters
    # ----------
    # target : an instance of 'FSource_provider' class
    #     Provides all the field in the long medium together with its respective grids [SI]
        
    # pressure : scalar or dict [SI]
    #     scalar is used to scale the result and compute the dispersion & absorption, 
    #     dict: 'zgrid' - pressure modulation zgrid; 'value' - the modulation \n
    # NOTE: grid might need to be finer to properly retrieve the XUV phase
                    

    #     order: (r,z,omega); there is no inter check of dimensions, it has to match
    # distance : scalar
    #     The distance from the first point of the medium [SI]
    # rgrid_FF : array_like
    #     the required radial grid in far-field [SI]
    # dispersion_function : function, optional
    #     The dispersion is factored as exp(i*z*omega*dispersion_function(omega)). The default is None.
    # absorption_function : function, optional
    #     Analogical to dispersion. The default is None.
    # integrator_Hankel : function, optional
    #     used method to integrate in the radial direction. The default is integrate.trapz.
    # integrator_longitudinal : function, optional
    #     used method to integrate in the longitudinal direction. The default is 'trapezoidal'.
    # near_field_factor : logical, optional
    #     False for far field without the Fresnel term. The default is True.    
    # frequencies_to_trace_maxima : list of 2D-array_like, optional
    #     If present, these windows given by the 2-D-array-likes are used to trace maxima of respective planes of integration.
    #     The default is None.

    # Raises
    # ------
    # NotImplementedError
    #     In the case a non-implemented integration rule is inputed.

    # Returns
    # -------
    # result : the field in the radial grid of the observation plane
    
    # result , planes_maxima: if frequencies_to_trace_maxima are present
    #     .

    # """

    

    if not(
            ((preset_gas+'_'+absorption_tables in XUV_index.gases)
             and
             (preset_gas+'_'+dispersion_tables in XUV_index.gases)
             )
             or (preset_gas == 'vacuum')
        ): raise ValueError('Wrongly specified preset gas (or tables).')

    No = len(target.ogrid);
    Nz = len(target.zgrid);
    Nr_FF = len(rgrid_FF)
    # include_dispersion = not(dispersion_function is None)
    # include_absorption = not(absorption_function is None)
    # trace_maxima_log = not(frequencies_to_trace_maxima is None)
    
    
    diagnostics = [[],[],[],[],[],[],[]] #integrands, intagrals, prefactor, cummulative, H-args, integrands inside, ad hoc integration #integrands
    
    
    # init pre_factor
    FF_integrated = 0.
    FF_integrated2 = 0.
    
    
    pre_factor, abs_factor = get_propagation_pre_factor_function(
                                    target.zgrid,
                                    target.rgrid,
                                    target.ogrid,
                                    preset_gas = preset_gas,
                                    pressure = pressure,
                                    absorption_tables = absorption_tables,
                                    include_absorption = include_absorption,
                                    dispersion_tables = dispersion_tables,
                                    include_dispersion = include_dispersion,
                                    effective_IR_refrective_index = effective_IR_refrective_index)

                          
                          

            
    # we keep the data for now, consider on-the-fly change
    print('Computing Hankel from planes')
    t_start  = time.perf_counter()
    t_check1 = t_start
    
    integrands_plane = next(target.Fsource_plane)
    
    diagnostics[0].append(integrands_plane)
    
    Fsource_plane1, integrands = HankelTransform(target.ogrid,
                                     target.rgrid,
                                     integrands_plane,
                                     distance-target.zgrid[0],
                                     rgrid_FF,
                                     integrator = integrator_Hankel,
                                     near_field_factor = near_field_factor,
                                     pre_factor = 1.)
    
    diagnostics[4].append([target.ogrid,target.rgrid,integrands_plane,
                           distance-target.zgrid[0],rgrid_FF,integrator_Hankel,
                           near_field_factor])
    diagnostics[5].append(integrands)
    
    
    diagnostics[1].append(copy.copy(Fsource_plane1))
    
    diagnostics[2].append(np.outer(pre_factor(0)[0,:],np.ones(Fsource_plane1.shape[1])))
    
    Fsource_plane1 *= np.outer(pre_factor(0)[0,:],np.ones(Fsource_plane1.shape[1]))
    
    
    
    if store_cummulative_result:
         cummulative_field = np.empty((Nz-1,) + Fsource_plane1.shape, dtype=np.cdouble)
    
    
    for k1 in range(Nz-1):
        t_check2 = time.perf_counter()
        print('plane', k1, 'time:', t_check2-t_start, 'this iteration: ', t_check2-t_check1)
        t_check1 = t_check2
        
        
        # ad-hoc integration
        pl1 = HankelTransform(target.ogrid,
                                                 target.rgrid,
                                                 integrands_plane,
                                                 distance-target.zgrid[k1],
                                                 rgrid_FF,
                                                 integrator = integrator_Hankel,
                                                 near_field_factor = near_field_factor,
                                                 pre_factor = 1.)[0]
        
        integrands_plane = next(target.Fsource_plane)
        
        pl2 = HankelTransform(target.ogrid,
                                                 target.rgrid,
                                                 integrands_plane,
                                                 distance-target.zgrid[k1+1],
                                                 rgrid_FF,
                                                 integrator = integrator_Hankel,
                                                 near_field_factor = near_field_factor,
                                                 pre_factor = 1.)[0]
        
        
        diagnostics[0].append(integrands_plane)
        
        Fsource_plane2, integrands = HankelTransform(target.ogrid,
                                         target.rgrid,
                                         integrands_plane,
                                         distance-target.zgrid[k1+1],
                                         rgrid_FF,
                                         integrator = integrator_Hankel,
                                         near_field_factor = near_field_factor,
                                         pre_factor = 1.)  
        
        diagnostics[4].append([target.ogrid,target.rgrid,integrands_plane,
                               distance-target.zgrid[k1+1],rgrid_FF,integrator_Hankel,
                               near_field_factor])
        diagnostics[5].append(integrands)
        
        
        diagnostics[1].append(copy.copy(Fsource_plane2))        
        
        
        diagnostics[2].append(np.outer(pre_factor(k1+1)[0,:],np.ones(Fsource_plane2.shape[1])))
        Fsource_plane2 *= np.outer(pre_factor(k1+1)[0,:],np.ones(Fsource_plane2.shape[1]))        
         
        
        FF_integrated += 0.5*(target.zgrid[k1+1]-target.zgrid[k1])*(Fsource_plane1 + Fsource_plane2)
        
        
        pl1 = np.outer(pre_factor(k1)[0,:],np.ones(pl1.shape[1]))*pl1
        
        # pl1 = np.outer(pre_factor(k1)[0,:],np.ones(pl1.shape[1]))*pl1/pressure
        
        pl2 = np.outer(pre_factor(k1+1)[0,:],np.ones(pl2.shape[1]))*pl2
        
        # pl2 = np.outer(pre_factor(k1+1)[0,:],np.ones(pl2.shape[1]))*pl2/pressure
        
        FF_integrated2 += 0.5*(target.zgrid[k1+1]-target.zgrid[k1])*(
                          pl1
                          +
                          pl2)
        
        
        diagnostics[6].append(copy.copy(FF_integrated2))
        
        
        
        diagnostics[3].append(copy.copy(FF_integrated))
        

    
        
        
        if store_cummulative_result:
            cummulative_field[k1,:,:] = copy.copy(FF_integrated)


        
        
        Fsource_plane1 = copy.copy(Fsource_plane2)
        

        

    
    
    
    if store_cummulative_result:
        return FF_integrated, cummulative_field, diagnostics
    else:
        return FF_integrated
        
        
        
        
        
        
        
        
        
    #     FSourceTerm_select = np.squeeze(FSourceTerm[:,k1,:]).T
    #     FField_FF = HankelTransform(ogrid,
    #                               rgrid,
    #                               FSourceTerm_select,
    #                               distance-zgrid[k1],
    #                               rgrid_FF,
    #                               integrator = integrator_Hankel,
    #                               near_field_factor = near_field_factor)
        
    #     if (k1 == 0): # allocate space
    #         FField_FF_z = np.zeros( (Nz,) + FField_FF.shape,dtype=np.cdouble) 
        
    #     if (include_dispersion or include_absorption):  
    #         FField_FF_z[k1,:,:] = np.outer(factor_e[k1,:],np.ones(FField_FF.shape[1]))*FField_FF
    #     else:
    #         FField_FF_z[k1,:,:] = FField_FF # (z,omega,r)
    
    # if store_cummulative_result:
    #     cummulative_field = np.empty((Nz-1,) + FField_FF.shape, dtype=np.cdouble)
        
    # if (integrator_longitudinal == 'trapezoidal'):        
    #     for k1 in range(Nz-1):    
    #         k_step = 1
    #         if (k1 == 0):
    #             dum = 0.5*(zgrid[(k1+1)*k_step]-zgrid[k1*k_step]) * \
    #                   (FField_FF_z[k1*k_step,:,:] + FField_FF_z[(k1+1)*k_step,:,:])
    #         else:
    #             dum = dum + \
    #                   0.5*(zgrid[(k1+1)*k_step]-zgrid[k1*k_step]) * \
    #                   (FField_FF_z[k1*k_step,:,:] + FField_FF_z[(k1+1)*k_step,:,:])
                      
    #         if store_cummulative_result:
    #             if include_absorption:
    #                 # we need renormalise the end of the medium
    #                 exp_renorm = np.exp( (zgrid[-1]-zgrid[k1]) * absorption_factor)
    #                 for k2 in range(No):
    #                     for k3 in range(Nr_FF):
    #                         cummulative_field[k1,k2,k3] = exp_renorm[k2]*dum[k2,k3]
    #             else:
    #                 cummulative_field[k1,:,:] = dum
                
    # else:
    #     raise NotImplementedError('Only trapezoidal rule implemented now')
        
    # if trace_maxima_log:
        
    #     frequency_indices = []
    #     planes_maxima = []
    #     for frequency_list in frequencies_to_trace_maxima:
    #         try:
    #             frequency_indices.append(mn.FindInterval(ogrid, frequency_list))
    #             planes_maxima.append([])
    #         except:
    #             warnings.warn("A frequency from frequencies_to_trace_maxima doesn't match ogrid.")
        
    #     if (len(frequency_indices)>0):
    #         for k1 in range(Nz):
    #             for k2 in range(len(frequency_indices)):
    #                 planes_maxima[k2].append(np.max(abs(
    #                                   FField_FF_z[k1,frequency_indices[k2][0]:frequency_indices[k2][1],:]
    #                                         )))
                    
    #         for k1 in range(len(frequency_indices)):
    #             planes_maxima[k1] = np.asarray(planes_maxima[k1])
        
    #     if store_cummulative_result:
    #         return dum , planes_maxima, cummulative_field
    #     else:
    #         return dum , planes_maxima
            
    # else: 
    #     if store_cummulative_result:
    #         return dum, cummulative_field
    #     else:
    #         return dum
                
                

def Signal_cum_integrator(ogrid, zgrid, FSourceTerm,
                         integrator = integrate.cumulative_trapezoid):
    
    No = len(ogrid); Nz = len(zgrid)
    signal = np.zeros((No,Nz), dtype=np.cdouble)
    integrand = FSourceTerm
    for k1 in range(No):
        signal[k1,1:] = integrator(integrand[k1,:],zgrid)
    return signal    
    












    # pressure_modulation = interpolate.interp1d(
    #                                 pressure['zgrid'],
    #                                 pressure['value'],
    #                                 bounds_error = False,
    #                                 fill_value = (pressure['value'][0],
    #                                               pressure['value'][-1]),
    #                                 copy = False
    #                                 )(zgrid)
        
    # phase_factor = np.empty(Nz,No); absorption_factor = np.empty(Nz,No)
    # for k1 in range(No):
    #     integral_for_phase_factor = integrate.cumtrapz(
    #                                             pressure['zgrid'],
    #                                             XUV_index.dispersion_function(
    #                                                 omega, 
    #                                                 pressure['value'],
    #                                                 preset_gas+'_'+dispersion_tables,
    #                                                 n_IR = effective_IR_refrective_index),
    #                                             initial=0.
    #                                             )
        
    #     phase_factor[:,k1] =  ogrid[k1] * \
    #                           interpolate.interp1d(
    #                                 pressure['zgrid'],
    #                                 integral_for_phase_factor,
    #                                 bounds_error = False,
    #                                 fill_value = (integral_for_phase_factor[0],
    #                                               integral_for_phase_factor[-1]),
    #                                 copy = False
    #                                 )(zgrid)
    
        
    #     integral_beta_factor = XUV_index.beta_factor_ref(
    #                                 omega,
    #                                 preset_gas+'_'+dispersion_tables) * \
    #                            integrate.cumtrapz(
    #                                 pressure['zgrid'],
    #                                 pressure['value'],
    #                                 initial=0.
    #                                 )
                               
    #     absorption_factor[:,k1] = ogrid[k1] * \ 
    #                               interpolate.interp1d(
    #                                 pressure['zgrid'],
    #                                 integral_beta_factor,
    #                                 bounds_error = False,
    #                                 fill_value = (integral_beta_factor[0],
    #                                               integral_beta_factor[-1]),
    #                                 copy = False
    #                               )(zgrid)

    
    
    # if include_dispersion:
    #     dispersion_factor = np.empty(No)
    #     for k1 in range(No):
    #         dispersion_factor[k1] = ogrid[k1]*dispersion_function(ogrid[k1])
        
    # if include_absorption:
    #     absorption_factor = np.empty(No)
    #     for k1 in range(No):
    #         absorption_factor[k1] = ogrid[k1]*absorption_function(ogrid[k1])
      
    # # compute z-evolution of the factors        
    # if (include_dispersion and include_absorption):
    #     factor_e = np.exp(
    #                       1j*np.outer(zgrid,dispersion_factor) +
    #                       np.outer(zgrid-zgrid[-1] ,absorption_factor)
    #                       )
    # elif include_dispersion:
    #     factor_e = np.exp(1j*np.outer(zgrid,dispersion_factor))

    # elif include_absorption:
    #     factor_e = np.exp(np.outer(zgrid-zgrid[-1] ,absorption_factor)) 