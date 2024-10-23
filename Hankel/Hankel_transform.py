"""
Created on Sun May 19 00:38:22 2024

@author: vabek
"""
import numpy as np
import mynumerics as mn
import time
import units
import XUV_refractive_index as XUV_index
from scipy import interpolate
from scipy import integrate
from scipy import special




class FSources_provider:
    """
    This class provides all the necessary inputs related to the
    field in a long medium:
        ogrid: frequency grid of the field
        rgrid: radial grid
        zgrid: longitudinal grid
        FSource: the source term in Fourier domain in the form of a generator.
                 This generator yields the field planes along the z-grid.
                 
        This structure is chosen because it allows flexible use for large inputs:
        The FField can be a whole static array in the memory, or it can be read
        plane-by-plne from the input hdf5.
        ! NOTE: if 'dynamic' is used, the source data-stream must be available
                (e.g. inside a 'with' block)
    """
    
    def __init__(self, # static=None,dynamic=None,
                 zgrid, rgrid, ogrid,
                 data_source = 'static',
                 FSource = None,
                 h5_handle = None,
                 h5_path = None,
                 ko_min  =  0,
                 ko_step =  1,
                 ko_max  = 'end',
                 kr_step =  1,
                 kr_max  = 'end',
                 kz_step =  1):

        # if (Nproc == 1)
        if (ko_max  == 'end'): ko_max = len(ogrid)
        if (kr_max  == 'end'): kr_max = len(rgrid)
        self.zgrid = zgrid[0:-1:kz_step]
        self.rgrid = rgrid[0:kr_max:kr_step]
        self.ogrid = ogrid[ko_min:ko_max:ko_step]
        
        if (data_source == 'static'):
            def FSource_plane_():
                for k1 in range(len(self.zgrid)):
                    yield FSource[k1*kz_step,ko_min:ko_max:ko_step,0:kr_max:kr_step]
            self.Fsource_plane = FSource_plane_()
        elif (data_source == 'dynamic'):
            def FSource_plane_():
                for k1 in range(len(self.zgrid)):
                    yield np.squeeze(
                            h5_handle[h5_path][k1*kz_step,0:kr_max:kr_step,ko_min:ko_max:ko_step,0]
                            +
                            1j*h5_handle[h5_path][k1*kz_step,0:kr_max:kr_step,ko_min:ko_max:ko_step,1]).T
            self.Fsource_plane = FSource_plane_()
        else:
            raise ValueError('Wrongly specified input of the class.')
            
            
            

    
    
    

def get_propagation_pre_factor_function(zgrid,
                                rgrid,
                                ogrid,
                                preset_gas = 'vacuum',
                                pressure = 1.,
                                absorption_tables = 'Henke',
                                include_absorption = True,
                                dispersion_tables = 'Henke',
                                include_dispersion = True,
                                effective_IR_refrective_index = 1.):
    """

    Returns
    -------
    int
        DESCRIPTION.

    """
    
    Nz = len(zgrid); Nr = len(rgrid); No = len(ogrid)
    
    
    # switch over geometries    
    if isinstance(pressure,dict):
        if (('zgrid' in pressure.keys()) and not('rgrid' in pressure.keys())):
            # usual case: omega-and-z-integrals, r on-the-fly
            print('z modulation')            
            
            pre_factor_value = np.empty((Nz,No),dtype=np.cdouble)
            absorption_factor_omega = np.empty((Nz,No),dtype=np.double)
            
            Nz_table = len(pressure['zgrid'])            
            pressure_modulation_local = interpolate.interp1d(
                                            pressure['zgrid'],
                                            pressure['value'],
                                            bounds_error = False,
                                            fill_value = (pressure['value'][0],
                                                          pressure['value'][-1]),
                                            copy = False
                                            )(zgrid)
            
            for k1 in range(No):
                if include_dispersion: 
                    
                    integrand = XUV_index.dispersion_function(
                        ogrid[k1], 
                        pressure['value'],                
                        preset_gas+'_'+dispersion_tables,
                        n_IR = effective_IR_refrective_index)
                    
                    integral_for_phase_factor = integrate.cumulative_trapezoid(
                                                        integrand,
                                                        x=pressure['zgrid'],
                                                        initial=0.
                                                        )
                    
                    dispersion_factor = 1j *\
                                        interpolate.interp1d(                               # dispersion factor
                                          pressure['zgrid'],
                                          integral_for_phase_factor,
                                          bounds_error = False,
                                          fill_value = (integral_for_phase_factor[0],
                                                        integral_for_phase_factor[-1]),
                                          copy = False
                                          )(zgrid)
                else:
                    # dispersion_factor = 0.    
                    # correction with respect to the co-moving frame according to the formula:
                    # dispersion_factor = ((1./v_co_moving) - (1./units.c_light)), v_co_moving = n*c
                    dispersion_factor = 1j * (effective_IR_refrective_index - 1.)/units.c_light 
                  
                    
                if include_absorption:
                    integral_beta_factor = XUV_index.beta_factor_ref(
                                            ogrid[k1],
                                            preset_gas+'_'+absorption_tables) * \
                                       integrate.cumulative_trapezoid(
                                            pressure['value'],
                                            x=pressure['zgrid'],
                                            initial=0.
                                            )
                                       
                    absorption_factor = (1./units.c_light) *\
                                        interpolate.interp1d(
                                          pressure['zgrid'],
                                          integral_beta_factor,
                                          bounds_error = False,
                                          fill_value = (integral_beta_factor[0],
                                                        integral_beta_factor[-1]),
                                          copy = False
                                        )(zgrid)
                     
                    absorption_factor_omega[:,k1] = ogrid[k1]*absorption_factor
                    absorption_factor = absorption_factor - absorption_factor[-1]
                                        
                else:   
                    absorption_factor = 0.
                    absorption_factor_omega[:,k1] = 0.
                    
                
                                       
                pre_factor_value[:,k1] = pressure_modulation_local *\
                                            np.exp(
                                                ogrid[k1] *
                                                    (dispersion_factor
                                                    + 
                                                    absorption_factor
                                                    )
                                            )
            
            def exp_renorm(kz):
                # return np.exp((zgrid[-1]-zgrid[kz]) * abs_factor_omega)
                return np.exp(absorption_factor_omega[-1,:] - absorption_factor_omega[kz,:])
            
            def pre_factor(kz):
                return np.outer(np.ones(len(rgrid)),np.squeeze(pre_factor_value[kz,:]))
            
            return pre_factor, exp_renorm
                
            
            
        elif (not('zgrid' in pressure.keys()) and ('rgrid' in pressure.keys())):
            print('r modulation')
            
            
            # First, we obtain the pressure modulation on the "(pressure z-grid, computational r-grid)"
            # pressure_my_rgrid = np.empty((Nr,))
            pressure_my_rgrid = interpolate.interp1d(
                                          pressure['rgrid'],
                                          pressure['value'][k1,:],
                                          bounds_error = False,
                                          fill_value = (pressure['value'][k1,0],
                                                        pressure['value'][k1,-1]),
                                          copy = False
                                          )(rgrid)
               
            dispersion_factor_omega = np.empty((Nr,No),dtype=np.cdouble)
            
            if include_dispersion:   
                print('dispersion applied')
                for k1 in range(Nr):
                    dispersion_factor_omega[k1,:] = 1j * ogrid * \
                                                    XUV_index.dispersion_function(
                                                            ogrid, 
                                                            pressure['value'][k1],             
                                                            preset_gas+'_'+dispersion_tables,
                                                            n_IR = effective_IR_refrective_index)
            else:
                print('no dispersion')
                # dispersion_factor = 0.    
                # correction with respect to the co-moving frame according to the formula:
                # dispersion_factor = ((1./v_co_moving) - (1./units.c_light)), v_co_moving = n*c
                for k1 in range(Nr):
                    dispersion_factor_omega[k1,:] = 1j * ogrid * \
                                                    (effective_IR_refrective_index - 1.)/units.c_light 
            
            absorption_factor_omega = np.zeros((Nr,No),dtype=np.double)   
            if include_absorption: 
                print('absorption applied')                
                for k1 in range(Nr):
                    absorption_factor_omega[k1,:] =  ogrid * \
                                                (pressure['value'][k1]/units.c_light) *\
                                                XUV_index.beta_factor_ref(
                                                    ogrid,
                                                    preset_gas+'_'+absorption_tables)
            else:
                print('no absorption')
                pass # done in the allocation

                
            def pre_factor(kz):
                return pressure * np.exp(
                                          (zgrid[kz]*dispersion_factor
                                           +
                                          (zgrid[kz]-zgrid[-1])*absorption_factor)
                                         )
                                                  

            
            return pre_factor, None # renormalisation not implemented

            
            
        elif (('zgrid' in pressure.keys()) and ('rgrid' in pressure.keys())):
            print('zr modulation')
            # pass
            # full case megamatrix: we store pre-factor in memory
            pre_factor_value = np.empty((Nz,No,Nr),dtype=np.cdouble)
            
            # First, we obtain the pressure modulation on the "(pressure z-grid, computational r-grid)"
            Nz_table = len(pressure['zgrid'])
            pressure_my_rgrid = np.empty((Nz_table,Nr))
            for k1 in range(Nz_table):
                pressure_my_rgrid[k1,:] = interpolate.interp1d(
                                              pressure['rgrid'],
                                              pressure['value'][k1,:],
                                              bounds_error = False,
                                              fill_value = (pressure['value'][k1,0],
                                                            pressure['value'][k1,-1]),
                                              copy = False
                                              )(rgrid)
            
                
            
            
            # We compute partial integrals for each r
            for k1 in range(Nr):
                
                # get pressure modulation for this particular r on the computational zgrid
                pressure_modulation_local = interpolate.interp1d(
                                                pressure['zgrid'],
                                                pressure_my_rgrid[:,k1],
                                                bounds_error = False,
                                                fill_value = (pressure_my_rgrid[0,k1],
                                                             pressure_my_rgrid[-1,k1]),
                                                copy = False
                                                )(zgrid)                        
                
                for k2 in range(No):                    
                    if include_dispersion: 
                        integral_for_phase_factor = integrate.cumulative_trapezoid(
                                                            XUV_index.dispersion_function(
                                                                ogrid[k2], 
                                                                pressure_my_rgrid[:,k1],
                                                                preset_gas+'_'+dispersion_tables,
                                                                n_IR = effective_IR_refrective_index),
                                                            x = pressure['zgrid'],
                                                            initial=0.
                                                            )
                        
                        dispersion_factor = 1j *\
                                            interpolate.interp1d(                                                
                                              pressure['zgrid'],
                                              integral_for_phase_factor,                                              
                                              bounds_error = False,
                                              fill_value = (integral_for_phase_factor[0],
                                                            integral_for_phase_factor[-1]),
                                              copy = False
                                              )(zgrid)
                    else:
                        # dispersion_factor = 0.    
                        # correction with respect to the co-moving frame according to the formula:
                        # dispersion_factor = ((1./v_co_moving) - (1./units.c_light)), v_co_moving = n*c
                        dispersion_factor = 1j * (effective_IR_refrective_index - 1.)/units.c_light 
                      
                        
                    if include_absorption:
                        integral_beta_factor = XUV_index.beta_factor_ref(
                                                ogrid[k2],
                                                preset_gas+'_'+dispersion_tables) * \
                                           integrate.cumulative_trapezoid(
                                                pressure_my_rgrid[:,k1],
                                                x = pressure['zgrid'],
                                                initial=0.
                                                )
                                           
                        absorption_factor = (1./units.c_light) *\
                                            interpolate.interp1d(
                                              pressure['zgrid'],
                                              integral_beta_factor,
                                              bounds_error = False,
                                              fill_value = (integral_beta_factor[0],
                                                            integral_beta_factor[-1]),
                                              copy = False
                                            )(zgrid)                   
                    else:   
                        absorption_factor = 0.
                                           
                    pre_factor_value[:,k2,k1] = pressure_modulation_local *\
                                                np.exp(
                                                    ogrid[k2] *
                                                        (dispersion_factor
                                                        + 
                                                        absorption_factor
                                                        )
                                                )
            def pre_factor(kz):
                return np.squeeze(pre_factor_value[kz,:,:])
            
            return pre_factor, None # renormalisation not implemented
            
        else:
            raise ValueError('Pressure modulation wrongly specified.')
    else:
        
        print('no modulation')
        if include_dispersion:   
            print('dispersion applied')
            dispersion_factor = 1j * XUV_index.dispersion_function(
                                            ogrid, 
                                            pressure,                                   # scalar               
                                            preset_gas+'_'+dispersion_tables,
                                            n_IR = effective_IR_refrective_index)
        else:
            print('no dispersion')
            # dispersion_factor = 0.    
            # correction with respect to the co-moving frame according to the formula:
            # dispersion_factor = ((1./v_co_moving) - (1./units.c_light)), v_co_moving = n*c
            dispersion_factor = 1j * (effective_IR_refrective_index - 1.)/units.c_light 
           
        if include_absorption: 
            print('absorption applied')
            absorption_factor = (pressure/units.c_light) *\
                                XUV_index.beta_factor_ref(
                                    ogrid,
                                    preset_gas+'_'+absorption_tables)
        else:
            print('no absorption')
            absorption_factor = 0.
        
        abs_factor_omega = ogrid * absorption_factor
        
        def exp_renorm(kz):
            return np.exp((zgrid[-1]-zgrid[kz]) * abs_factor_omega)
            
        def pre_factor(kz):
            return pressure * np.outer(np.ones(len(rgrid)),
                                        np.exp(
                                              ogrid * (zgrid[kz]*dispersion_factor
                                                        +
                                                        (zgrid[kz]-zgrid[-1])*absorption_factor)
                                              )
                                        )
        
        return pre_factor, exp_renorm 
        
        
        
def HankelTransform(ogrid, rgrid, FField, distance, rgrid_FF,
                    integrator = lambda y, x: integrate.trapezoid(y,x=x),
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
       
    
    t_start = time.perf_counter()
    
    apply_radial_factor = (len(np.shape(pre_factor))==2)
    
    
    No = len(ogrid); Nr = len(rgrid); Nr_FF = len(rgrid_FF)
    FField_FF = np.empty((No,Nr_FF), dtype=np.cdouble)
    
    
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


    print('time spent only in the integrator ', time.perf_counter()-t_start)
    
    return FField_FF

        
        
class Hankel_long:
    def __init__(self,
                 target, 
                 distance,
                 rgrid_FF,
                 preset_gas = 'vacuum',
                 pressure = 1.,
                 absorption_tables = 'Henke',
                 include_absorption = True,
                 dispersion_tables = 'Henke',
                 include_dispersion = True,
                 effective_IR_refrective_index = 1.,
                 integrator_Hankel = lambda y, x: integrate.trapezoid(y,x=x), # integrate.trapz,
                 integrator_longitudinal = 'trapezoidal',
                 near_field_factor = True,
                 store_cummulative_result = False,
                 store_non_normalised_cummulative_result = False,
                 store_entry_and_exit_plane_transform = True
                 ):
        
        
        # keep some inputs to keep data packed together
        self.include_dispersion = include_dispersion
        if include_dispersion: self.dispersion_tables = dispersion_tables
        self.include_absorption = include_absorption
        if include_absorption: self.absorption_tables = absorption_tables
        self.effective_IR_refrective_index = effective_IR_refrective_index
        self.integrator_Hankel = integrator_Hankel
        self.integrator_longitudinal = integrator_longitudinal
        self.near_field_factor = near_field_factor
        
        self.rgrid = rgrid_FF
        self.ogrid = np.copy(target.ogrid)
        
        if (store_cummulative_result or store_non_normalised_cummulative_result):
            self.zgrid = np.copy(target.zgrid)
        
        if not(
                ((preset_gas+'_'+absorption_tables in XUV_index.gases)
                 and
                 (preset_gas+'_'+dispersion_tables in XUV_index.gases)
                 )
                 or (preset_gas == 'vacuum')
            ): raise ValueError('Wrongly specified preset gas (or tables).')
    
        Nz = len(target.zgrid)
        Nr_FF = len(rgrid_FF)    

        
        # init pre_factor
        pre_factor, renorm_factor = get_propagation_pre_factor_function(
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
        Fsource_plane1 = HankelTransform(target.ogrid,
                                         target.rgrid,
                                         integrands_plane,
                                         distance-target.zgrid[0],
                                         rgrid_FF,
                                         integrator = integrator_Hankel,
                                         near_field_factor = near_field_factor,
                                         pre_factor = pre_factor(0)).T

        if store_cummulative_result:
             cummulative_field = np.empty((Nz-1,) + Fsource_plane1.shape, dtype=np.cdouble)
        if store_non_normalised_cummulative_result:
            cummulative_field_no_norm = np.empty((Nz-1,) + Fsource_plane1.shape, dtype=np.cdouble)        
        if store_entry_and_exit_plane_transform:
            self.entry_plane_transform = np.copy(Fsource_plane1)        
        
        FF_integrated = 0.
        for k1 in range(Nz-1):
            t_check2 = time.perf_counter()
            print('plane', k1, 'time:', t_check2-t_start, 'this iteration: ', t_check2-t_check1)
            t_check1 = t_check2
            
          
            integrands_plane = next(target.Fsource_plane) 
            Fsource_plane2 = HankelTransform(target.ogrid,
                                             target.rgrid,
                                             integrands_plane,
                                             distance-target.zgrid[k1+1],
                                             rgrid_FF,
                                             integrator = integrator_Hankel,
                                             near_field_factor = near_field_factor,
                                             pre_factor = pre_factor(k1+1)).T

            FF_integrated += 0.5*(target.zgrid[k1+1]-target.zgrid[k1])*(Fsource_plane1 + Fsource_plane2)
            
            if store_cummulative_result:
                
                if isinstance(pressure,dict): 
                  if ('rgrid' in pressure.keys()):
                    raise NotImplementedError('Renormalisation of the signal is not implemented for radially modulated density.')
                cummulative_field[k1,:,:] = np.outer(np.ones(Nr_FF),renorm_factor(k1))*FF_integrated
                
                
            if store_non_normalised_cummulative_result:
                cummulative_field_no_norm[k1,:,:]  =  np.copy(FF_integrated)
    
            Fsource_plane1 = np.copy(Fsource_plane2)
            

        self.FF_integrated = FF_integrated
        
        if store_entry_and_exit_plane_transform:
            self.exit_plane_transform = np.copy(Fsource_plane2)
        
        if store_cummulative_result:
            self.cummulative_field = cummulative_field
            
        if store_non_normalised_cummulative_result:
            self.cummulative_field_no_norm = cummulative_field_no_norm
           
            
           
        
def Signal_cum_integrator(ogrid, zgrid, FSourceTerm,
                         integrator = integrate.cumulative_trapezoid):
    
    No = len(ogrid); Nz = len(zgrid)
    signal = np.zeros((No,Nz), dtype=np.cdouble)
    integrand = FSourceTerm
    for k1 in range(No):
        signal[k1,1:] = integrator(integrand[k1,:],zgrid)
    return signal     