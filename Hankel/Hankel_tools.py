import numpy as np
from scipy import interpolate
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
    
    def __init__(self,static=None,dynamic=None):
        if (isinstance(static,dict) and (dynamic is None)):
            self.ogrid = static['ogrid']
            self.rgrid = static['rgrid']
            self.zgrid = static['zgrid']
            def FSource_plane_():
                for k1 in range(len(self.zgrid)):
                    yield static['FSource'][k1,:,:]
            self.Fsource_plane = FSource_plane_()
            
        elif ((static is None) and isinstance(dynamic,dict)):
            self.ogrid = dynamic['ogrid']
            self.rgrid = dynamic['rgrid']
            self.zgrid = dynamic['zgrid']
            def FSource_plane_():
                for k1 in range(len(self.zgrid)):
                    yield np.squeeze(
                            dynamic['FSource'][:,k1,:,0] + 1j*dynamic['FSource'][:,k1,:,1])
                    
            self.Fsource_plane = FSource_plane_()
            
        else:
            raise ValueError('Wrongly specified input of the class.')



# this function will provide a fuction that gives the pre-factor for the integration in long media:
# it loads preset gases - easy use for physical situations
#  

def get_propagation_pre_factor():
    """

    Returns
    -------
    int
        DESCRIPTION.

    """
    
    # switch over geometries
    
    if isinstance(pressure,dict):
        if (('zgrid' in pressure.keys()) and not('rgrid' in pressure.keys())):
            pass
            # usual case: omega-and-z-integrals, r on-the-fly
        elif (not('zgrid' in pressure.keys()) and ('rgrid' in pressure.keys())):
            pass
            # no integrals, only r-scaling, z on-the-fly
        elif (('zgrid' in pressure.keys()) and ('rgrid' in pressure.keys())):
            pass
            # full case megamatrix
        else:
            raise ValueError('Pressure modulation wrongly specified.')
    else:
        pass
        # from previous version: no integrals, np.outer-stuff
            
            
    
    def pre_factor(kz):
        pressure_modulation = interpolate.interp1d(
                                        pressure['zgrid'],
                                        pressure['value'],
                                        bounds_error = False,
                                        fill_value = (pressure['value'][0],
                                                      pressure['value'][-1]),
                                        copy = False
                                        )(zgrid)
            
        phase_factor = np.empty(Nz,No); absorption_factor = np.empty(Nz,No)
        for k1 in range(No):
            integral_for_phase_factor = integrate.cumtrapz(
                                                    pressure['zgrid'],
                                                    XUV_index.dispersion_function(
                                                        omega, 
                                                        pressure['value'],
                                                        preset_gas+'_'+dispersion_tables,
                                                        n_IR = effective_IR_refrective_index),
                                                    initial=0.
                                                    )
            
            phase_factor[:,k1] =  ogrid[k1] * \
                                  interpolate.interp1d(
                                        pressure['zgrid'],
                                        integral_for_phase_factor,
                                        bounds_error = False,
                                        fill_value = (integral_for_phase_factor[0],
                                                      integral_for_phase_factor[-1]),
                                        copy = False
                                        )(zgrid)
        
            
            integral_beta_factor = XUV_index.beta_factor_ref(
                                        omega,
                                        preset_gas+'_'+dispersion_tables) * \
                                   integrate.cumtrapz(
                                        pressure['zgrid'],
                                        pressure['value'],
                                        initial=0.
                                        )
                                   
            absorption_factor[:,k1] = ogrid[k1] * \
                                      interpolate.interp1d(
                                        pressure['zgrid'],
                                        integral_beta_factor,
                                        bounds_error = False,
                                        fill_value = (integral_beta_factor[0],
                                                      integral_beta_factor[-1]),
                                        copy = False
                                      )(zgrid)
                                      
        return pre_factor_
    
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