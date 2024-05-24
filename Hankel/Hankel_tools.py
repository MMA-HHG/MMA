import numpy as np
from scipy import integrate
from scipy import interpolate
import units
# import h5py
import XUV_refractive_index as XUV_index
import MMA_administration as MMA



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
                    # !!!! THIS IS SUPER SLOW, try:
                    # https://docs.h5py.org/en/stable/high/dataset.html#:~:text=)%0Adataset%20inaccessible-,read_direct,-(array%2C
                    # https://stackoverflow.com/a/72397433
            self.Fsource_plane = FSource_plane_()
        else:
            raise ValueError('Wrongly specified input of the class.')



# this function will provide a fuction that gives the pre-factor for the integration in long media:
# it loads preset gases - easy use for physical situations
#  

    
    
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
            print('z modulation')
            # pass
            # usual case: omega-and-z-integrals, r on-the-fly
            
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
                    
                    integral_for_phase_factor = integrate.cumtrapz(
                                                        integrand,
                                                        pressure['zgrid'],
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
                    dispersion_factor = 0.
                  
                    
                if include_absorption:
                    integral_beta_factor = XUV_index.beta_factor_ref(
                                            ogrid[k1],
                                            preset_gas+'_'+absorption_tables) * \
                                       integrate.cumtrapz(
                                            pressure['value'],
                                            pressure['zgrid'],
                                            initial=0.
                                            )
                                       
                    absorption_factor = (1./units.c_light) *\
                                        interpolate.interp1d(                               # absorption factor
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
            def pre_factor(kz):
                return 1.
            
            return pre_factor
            # pass
            # no integrals, only r-scaling, z on-the-fly
            
            
        elif (('zgrid' in pressure.keys()) and ('rgrid' in pressure.keys())):
            print('zr modulation')
            # pass
            # full case megamatrix
            # NOW: ASSUMIG (Z,R)-order
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
                
            pre_factor_value = np.empty((Nz,No,Nr),dtype=np.cdouble)
            for k1 in range(Nr):
                
                # get pressure modulation for this particular r on our zgrid
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
                        integral_for_phase_factor = integrate.cumtrapz(
                                                            pressure['zgrid'],
                                                            XUV_index.dispersion_function(
                                                                ogrid[k2], 
                                                                pressure_my_rgrid[:,k1],                # value - matrix
                                                                preset_gas+'_'+dispersion_tables,
                                                                n_IR = effective_IR_refrective_index),
                                                            initial=0.
                                                            )
                        
                        dispersion_factor = 1j *\
                                            interpolate.interp1d(                               # dispersion factor
                                              integral_for_phase_factor,
                                              pressure['zgrid'],
                                              bounds_error = False,
                                              fill_value = (integral_for_phase_factor[0],
                                                            integral_for_phase_factor[-1]),
                                              copy = False
                                              )(zgrid)
                    else:
                        dispersion_factor = 0.
                      
                        
                    if include_absorption:
                        integral_beta_factor = XUV_index.beta_factor_ref(
                                                ogrid[k2],
                                                preset_gas+'_'+dispersion_tables) * \
                                           integrate.cumtrapz(
                                                pressure_my_rgrid[:,k1],
                                                pressure['zgrid'],
                                                initial=0.
                                                )
                                           
                        absorption_factor = (1./units.c_light) *\
                                            interpolate.interp1d(                               # absorption factor
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
            
            return pre_factor
            
            
                    # phase_factor[:,k1] =  ogrid[k1] * \
                    #                       interpolate.interp1d(
                    #                             pressure['zgrid'],
                    #                             integral_for_phase_factor,
                    #                             bounds_error = False,
                    #                             fill_value = (integral_for_phase_factor[0],
                    #                                           integral_for_phase_factor[-1]),
                    #                             copy = False
                    #                             )(zgrid)
                                          
                    
                    # pre_factor_value[:,k2,k1] = pressure['value'][:,k1] *\
                    #                             np.exp(
                    #                                 ogrid[k2] *
                    #                                     (
                    #                                     1j *
                    #                                     interpolate.interp1d(                               # dispersion factor
                    #                                       pressure['zgrid'],
                    #                                       integral_for_phase_factor,
                    #                                       bounds_error = False,
                    #                                       fill_value = (integral_for_phase_factor[0],
                    #                                                     integral_for_phase_factor[-1]),
                    #                                       copy = False
                    #                                       )(zgrid)
                    #                                     + 
                    #                                     interpolate.interp1d(                               # absorption factor
                    #                                       pressure['zgrid'],
                    #                                       integral_beta_factor,
                    #                                       bounds_error = False,
                    #                                       fill_value = (integral_beta_factor[0],
                    #                                                     integral_beta_factor[-1]),
                    #                                       copy = False
                    #                                     )(zgrid)
                    #                                     )
                    #                             )
        else:
            raise ValueError('Pressure modulation wrongly specified.')
    else:
        # from previous version: no integrals, np.outer-stuff
        # pass
    
        print('no modulation')
        if not(preset_gas == 'vacuum'):
            if include_dispersion:   
                print('dispersion applied')
                dispersion_factor = 1j * XUV_index.dispersion_function(
                                                ogrid, 
                                                pressure,                                   # scalar               
                                                preset_gas+'_'+dispersion_tables,
                                                n_IR = effective_IR_refrective_index)
            else:
                print('no dispersion')
                dispersion_factor = 0.    
     
        
        
            if include_absorption: 
                print('absorption applied')
                absorption_factor = (pressure/units.c_light) *\
                                    XUV_index.beta_factor_ref(
                                        ogrid,
                                        preset_gas+'_'+absorption_tables)
            else:
                print('no absorption')
                absorption_factor = 0.
                
            
            lin_prop_factor = ogrid * (dispersion_factor + absorption_factor)
            
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
        
        else:
            def pre_factor(kz):
                return 1.
            
            return pre_factor

    # def pre_factor_empty(kz):
    #     return 1.
    
    # return pre_factor_empty



def pressure_constructor(h5_handle):
    
    # obtain the base unmodulated pressure
    pressure_base = h5_handle[MMA.paths['CUPRAD_inputs']+'/medium_pressure_in_bar'][()]
    # check whether the density modulation is applied
    if ('density_mod' in h5_handle[MMA.paths['global_inputs']].keys()):
        dens_mod_grp = h5_handle[MMA.paths['global_inputs']+'/'+'density_mod']
        pressure = {'value' : pressure_base*dens_mod_grp['table'][:]}
        if ('zgrid' in dens_mod_grp.keys()):
            pressure['zgrid'] = dens_mod_grp['zgrid'][:]
        if ('rgrid' in dens_mod_grp.keys()):
            pressure['rgrid'] = dens_mod_grp['rgrid'][:]
    else:
        pressure = pressure_base
        
    return pressure



# class linear_propagation_e_factor():
#     """

#     Returns
#     -------
#     int
#         DESCRIPTION.

#     """
#     def __init__(self,static=None,dynamic=None):
#         pass
    
#     def e_factor(kz):
#         pass
    
    
    


     
# np.transpose(
#         np.squeeze(
#             dynamic['FSource'][:,k1,:,0] + 1j*dynamic['FSource'][:,k1,:,1]),
#     axes=())