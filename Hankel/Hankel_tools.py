import numpy as np
import h5py



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
     
            
     
# np.transpose(
#         np.squeeze(
#             dynamic['FSource'][:,k1,:,0] + 1j*dynamic['FSource'][:,k1,:,1]),
#     axes=())