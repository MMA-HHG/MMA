"""
Created on Sun May 19 00:38:22 2024

@author: vabek
"""
import numpy as np
import time
import units
import XUV_refractive_index as XUV_index
from scipy import integrate
from scipy import special

from Hankel_tools import get_propagation_pre_factor_function



# for testing now
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

        
        
class Hankel_long:
    def __init__(self,
                 target, # FSourceTerm(r,z,omega)
                 distance,
                 rgrid_FF,
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
                 store_non_normalised_cummulative_result = False,
                 store_entry_and_exit_plane_transform = True
                 ):
        
        
        # keep some inputs for diagnostics
        self.include_dispersion = include_dispersion
        if include_dispersion: self.dispersion_tables = dispersion_tables
        self.include_absorption = include_absorption
        if include_absorption: self.absorption_tables = absorption_tables
        self.effective_IR_refrective_index = effective_IR_refrective_index
        self.integrator_Hankel = integrator_Hankel
        self.integrator_longitudinal = integrator_longitudinal
        self.near_field_factor = near_field_factor
        
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
        
        
        pre_factor, abs_factor_omega = get_propagation_pre_factor_function(
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
        
        if store_entry_and_exit_plane_transform:
            self.entry_plane_transform = np.copy(Fsource_plane1)
        
        diagnostics[4].append([target.ogrid,target.rgrid,integrands_plane,
                               distance-target.zgrid[0],rgrid_FF,integrator_Hankel,
                               near_field_factor])
        diagnostics[5].append(integrands)
        
        
        diagnostics[1].append(copy.copy(Fsource_plane1))
        
        diagnostics[2].append(np.outer(pre_factor(0)[0,:],np.ones(Fsource_plane1.shape[1])))
        
        Fsource_plane1 *= np.outer(pre_factor(0)[0,:],np.ones(Fsource_plane1.shape[1]))
        
        
        
        if store_cummulative_result:
             cummulative_field = np.empty((Nz-1,) + Fsource_plane1.shape, dtype=np.cdouble)
        if store_non_normalised_cummulative_result:
            cummulative_field_no_norm = np.empty((Nz-1,) + Fsource_plane1.shape, dtype=np.cdouble)
        
        
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
                if isinstance(pressure,dict): 
                  if ('rgrid' in pressure.keys()):
                    raise NotImplementedError('Renormalisation of the signal is not implemented for radially modulated density.')
                
                exp_renorm = np.exp( (target.zgrid[-1]-target.zgrid[k1]) * abs_factor_omega)
                for k2 in range(No):
                    for k3 in range(Nr_FF):
                        cummulative_field[k1,k2,k3] = exp_renorm[k2]*FF_integrated[k2,k3]
                
                # cummulative_field[k1,:,:] = np.copy(FF_integrated)
            if store_non_normalised_cummulative_result:
                cummulative_field_no_norm[k1,:,:]  =  np.copy(FF_integrated)
    
    
            
            
            Fsource_plane1 = copy.copy(Fsource_plane2)
            
    
            
    
        
        self.FF_integrated = FF_integrated
        
        if store_entry_and_exit_plane_transform:
            self.exit_plane_transform = np.copy(Fsource_plane2)
        
        if store_cummulative_result:
            self.cummulative_field = cummulative_field
            
        if store_non_normalised_cummulative_result:
            self.cummulative_field_no_norm = cummulative_field_no_norm
            
        self.diagnostics = diagnostics
         
        # if store_cummulative_result:
        #     return FF_integrated, cummulative_field, diagnostics
        # else:
        #     return FF_integrated
