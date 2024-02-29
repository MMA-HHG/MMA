import numpy as np
from scipy import integrate
from scipy import interpolate
import units
# import h5py
import XUV_refractive_index as XUV_index



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
    
    def __init__(self,static=None,dynamic=None,
                 ko_step =  1,
                 No_max  = 'end',
                 kr_step =  1,
                 Nr_max  = 'end',
                 kz_step =  1):
        
        
        # ko_step = 1; No_max = -1
        # kr_step = 1; Nr_max = -1
        # kz_step = 1
        # if not(reduce_resolution is None):
        #     if ('ko_step' in reduce_resolution.keys()): ko_step = reduce_resolution['ko_step']
        #     if ('No_max'  in reduce_resolution.keys()): No_max  = reduce_resolution['No_max']
        #     if ('kr_step' in reduce_resolution.keys()): kr_step = reduce_resolution['kr_step']
        #     if ('Nr_max'  in reduce_resolution.keys()): Nr_max  = reduce_resolution['Nr_max']
        #     if ('kz_step' in reduce_resolution.keys()): kz_step = reduce_resolution['kz_step']
        
        
        if (isinstance(static,dict) and (dynamic is None)):
            
            if (No_max  == 'end'): No_max = len(static['ogrid'])
            if (Nr_max  == 'end'): Nr_max = len(static['rgrid'])
            
            self.ogrid = static['ogrid'][0:No_max:ko_step]
            self.rgrid = static['rgrid'][0:Nr_max:kr_step]
            self.zgrid = static['zgrid'][0:-1:kz_step]
            def FSource_plane_():
                for k1 in range(len(self.zgrid)):
                    yield static['FSource'][k1*kz_step,0:No_max:ko_step,0:Nr_max:kr_step]
            self.Fsource_plane = FSource_plane_()
            
        elif ((static is None) and isinstance(dynamic,dict)):
            
            if (No_max  == 'end'): No_max = len(dynamic['ogrid'])
            if (Nr_max  == 'end'): Nr_max = len(dynamic['rgrid'])
            
            self.ogrid = dynamic['ogrid']
            self.rgrid = dynamic['rgrid']
            self.zgrid = dynamic['zgrid']
            def FSource_plane_():
                for k1 in range(len(self.zgrid)):
                    yield np.squeeze(
                            dynamic['h5_file'][dynamic['Fsource_path']][0:Nr_max:kr_step,k1*kz_step,0:No_max:ko_step,0]
                            +
                            1j*dynamic['h5_file'][dynamic['Fsource_path']][0:Nr_max:kr_step,k1*kz_step,0:No_max:ko_step,1]).T
                    
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
            # pass
            # usual case: omega-and-z-integrals, r on-the-fly
            
            pre_factor_value = np.empty((Nz,No),dtype=np.cdouble)
            
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
                    integral_for_phase_factor = integrate.cumtrapz(
                                                        pressure['zgrid'],
                                                        XUV_index.dispersion_function(
                                                            ogrid[k1], 
                                                            pressure['value'],                
                                                            preset_gas+'_'+dispersion_tables,
                                                            n_IR = effective_IR_refrective_index),
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
                                            preset_gas+'_'+dispersion_tables) * \
                                       integrate.cumtrapz(
                                            pressure['zgrid'],
                                            pressure['value'],
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
                                       
                pre_factor_value[:,k1] = pressure_modulation_local *\
                                            np.exp(
                                                ogrid[k1] *
                                                    (dispersion_factor
                                                    + 
                                                    absorption_factor
                                                    )
                                            )
                
            def pre_factor(kz):
                return np.outer(np.ones(len(rgrid)),np.squeeze(pre_factor_value[kz,:]))
                
            
            
        elif (not('zgrid' in pressure.keys()) and ('rgrid' in pressure.keys())):
            pass
            # no integrals, only r-scaling, z on-the-fly
            
            
        elif (('zgrid' in pressure.keys()) and ('rgrid' in pressure.keys())):
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
                                                ogrid[k2],
                                                preset_gas+'_'+dispersion_tables) * \
                                           integrate.cumtrapz(
                                                pressure['zgrid'],
                                                pressure_my_rgrid[:,k1],
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
        if include_dispersion:                
                dispersion_factor = 1j * XUV_index.dispersion_function(
                                            ogrid, 
                                            pressure,                                   # scalar               
                                            preset_gas+'_'+dispersion_tables,
                                            n_IR = effective_IR_refrective_index)
        else:
            dispersion_factor = 0.    
 
    
    
        if include_absorption:                
            absorption_factor = (pressure/units.c_light) *\
                                XUV_index.beta_factor_ref(
                                    ogrid,
                                    preset_gas+'_'+dispersion_tables)
        else:
            absorption_factor = 0.
            
        
        lin_prop_factor = ogrid * (dispersion_factor + absorption_factor)
        
        
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
         
            
        def pre_factor(kz):
            return pressure * np.exp(
                                zgrid[kz]*np.outer(
                                    np.ones(len(rgrid)),lin_prop_factor
                                    )
                                )            
        
        # pre_factor_value = np.empty(Nz,No,Nr)
        # for k1 in range(Nr):
        #     for k2 in range(No):
        #         for k3 in range(Nz):
        #             pre_factor_value[k3,k2,k1] = ogrid[k2]*dispersion_function(ogrid[k2])
        
            
            
    
    # def pre_factor(kz):
    #     pressure_modulation = interpolate.interp1d(
    #                                     pressure['zgrid'],
    #                                     pressure['value'],
    #                                     bounds_error = False,
    #                                     fill_value = (pressure['value'][0],
    #                                                   pressure['value'][-1]),
    #                                     copy = False
    #                                     )(zgrid)
            
    #     phase_factor = np.empty(Nz,No); absorption_factor = np.empty(Nz,No)
    #     for k1 in range(No):
    #         integral_for_phase_factor = integrate.cumtrapz(
    #                                                 pressure['zgrid'],
    #                                                 XUV_index.dispersion_function(
    #                                                     omega, 
    #                                                     pressure['value'],
    #                                                     preset_gas+'_'+dispersion_tables,
    #                                                     n_IR = effective_IR_refrective_index),
    #                                                 initial=0.
    #                                                 )
            
    #         phase_factor[:,k1] =  ogrid[k1] * \
    #                               interpolate.interp1d(
    #                                     pressure['zgrid'],
    #                                     integral_for_phase_factor,
    #                                     bounds_error = False,
    #                                     fill_value = (integral_for_phase_factor[0],
    #                                                   integral_for_phase_factor[-1]),
    #                                     copy = False
    #                                     )(zgrid)
        
            
    #         integral_beta_factor = XUV_index.beta_factor_ref(
    #                                     omega,
    #                                     preset_gas+'_'+dispersion_tables) * \
    #                                integrate.cumtrapz(
    #                                     pressure['zgrid'],
    #                                     pressure['value'],
    #                                     initial=0.
    #                                     )
                                   
    #         absorption_factor[:,k1] = ogrid[k1] * \
    #                                   interpolate.interp1d(
    #                                     pressure['zgrid'],
    #                                     integral_beta_factor,
    #                                     bounds_error = False,
    #                                     fill_value = (integral_beta_factor[0],
    #                                                   integral_beta_factor[-1]),
    #                                     copy = False
    #                                   )(zgrid)
                                      
    #     return pre_factor_
    
    return pre_factor

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