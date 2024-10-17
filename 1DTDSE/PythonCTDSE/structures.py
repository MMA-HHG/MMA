"""
Structures
==========

This module contains the C-compatible structures definitions for the Python-C binding. 
"""


from PythonCTDSE.ctypes_helper import *
from typing import Any
import h5py
from PythonCTDSE.constants import *
import MMA_administration as MMA

### Define structures

### Field structure
class Efield_var(Structure):
    _fields_ = [
        ("tgrid", POINTER(c_double)),
        ("Field", POINTER(c_double)),
        ("dt", c_double),
        ("omega", c_double),
        ("E0", c_double),
        ("phi", c_double),
        ("ton", c_double),
        ("toff", c_double),
        ("Nt", c_int),
        ("nc", c_int)
    ]

class trg_def(Structure):
    _fields_ = [
        ("a", c_double)
    ]

class analy_def(Structure):
    _fields_ = [
        ("tprint", c_double),
        ("writewft", c_int)
    ]

class output_print_def(Structure):
    _fields_ = [
        ("Efield", c_int),
        ("FEfield", c_int),
        ("sourceterm", c_int),
        ("Fsourceterm", c_int),
        ("FEfieldM2", c_int),
        ("FsourceTermM2", c_int),
        ("PopTot", c_int),
        ("tgrid", c_int),
        ("omegagrid", c_int),
        ("PopInt", c_int),
        ("expval_x", c_int)
    ]

class inputs_def(Structure):
    """
    Input structure with input data for the 1DTDSE computation.

    Attributes:
    -----------
    trg:
        Target definition structure.
    Efield:
        Electric field structure.
    Eguess:
        GS energy.
    Einit:
        Initial guess for the GS computation.
    tmin:
        Minimum time.
    Nt:
        TDSE temporal resolution.
    num_t:
        Number of points per 800nm field cycle.
    dt:
        Time step.
    num_r:
        TDSE spatial resolution.
    dx:
        Spatial step.
    psi0:
        Ground state (GS) wavefunction.
    x:
        Spatial grid.
    gauge:
        TDSE gauge: 0 == length, 1 == velocity (not implemented yet!)
    x_int:
        Electron density in range (x-x_int, x+x_int)
    analy:
        Analytical values print.
    InterpByDTorNT:
        Interpolate by timestep (dt == 1) or refine init dt by ```Ntinterp``` points (== 0).
    Ntinterp:
        Number of points for initial dt refinement.
    Print:
        Output prints structure.
    CV:
        Convergence of the GS.
    Precision:
        Floating point precision.
    """

    def __init__(self, *args: Any, **kw: Any):
        super().__init__(*args, **kw)
        self.ptr = byref(self)

    _fields_ = [
        ("trg", trg_def),
        ("Efield", Efield_var),
        ("Eguess", c_double),
        ("Einit", c_double),
        ("tmin", c_double),
        ("Nt", c_int),
        ("num_t", c_int),
        ("dt", c_double),
        ("num_r", c_int),
        ("dx", c_double),
        ("psi0", POINTER(c_double)),
        ("x", POINTER(c_double)),
        ("gauge", c_int),
        ("x_int", c_double),
        ("analy", analy_def),
        ("InterpByDTorNT", c_int),
        ("Ntinterp", c_int),
        ("Print", output_print_def),
        ("CV", c_double),
        ("precision", c_char * 2)
    ]

    def load_from_hdf5(self, filename):
        """
        Initializes input structure from an HDF5 archive.

        Parameters:
        -----------
        filename: str
            Path to hdf5 archive.
        """
        with h5py.File(filename, "r") as f:
            self.Eguess = c_double(f["TDSE_inputs/Eguess"][()])
            self.num_r = c_int(f["TDSE_inputs/N_r_grid"][()])
            self.dx = c_double(f["TDSE_inputs/dx"][()])
            self.InterpByDTorNT = c_int(f["TDSE_inputs/InterpByDTorNT"][()])
            self.dt = c_double(f["TDSE_inputs/dt"][()])
            self.Ntinterp = c_int(f["TDSE_inputs/Ntinterp"][()])
            self.analy.writewft = c_int(f["TDSE_inputs/analy_writewft"][()])
            self.analy.tprint = c_double(f["TDSE_inputs/analy_tprint"][()])
            self.x_int = c_double(f["TDSE_inputs/x_int"][()])
            self.trg.a = c_double(f["TDSE_inputs/trg_a"][()])
            self.CV = c_double(f["TDSE_inputs/CV_criterion_of_GS"][()])
            self.gauge = c_int(f["TDSE_inputs/gauge_type"][()])
            precision = 'd'
            self.precision = precision.encode('utf-8')

            try: 
                x_grid = np.array(f["TDSE_inputs/x_grid"][()])
                self.x = ctypes_arr_ptr(c_double, self.num_r+1, x_grid)
                psi0 = np.array(f["TDSE_inputs/psi0"][()]).flatten()
                self.psi0 = ctypes_arr_ptr(c_double, 2*(self.num_r+1), psi0)
                self.Einit = c_double(f["TDSE_inputs/Einit"][()])
            except KeyError:
                pass

            try: 
                self.Efield.Nt = c_int(f["TDSE_inputs/Nt"][()])
                Efield = np.array(f["TDSE_inputs/Efield"][()])
                self.Efield.Field = ctypes_arr_ptr(c_double, self.Efield.Nt, Efield)
                tgrid = np.array(f["TDSE_inputs/tgrid"][()])
                self.Efield.tgrid = ctypes_arr_ptr(c_double, self.Efield.Nt, tgrid)
            except KeyError:
                pass

    def init_default_inputs(self,
                            Eguess = -1.,
                            num_r = 16000,
                            dx = 0.4,
                            InterpByDTorNT = 0,
                            dt = 0.25,
                            trg_a = 1.3677,
                            CV = 1e-25,
                            gauge = 0,
                            Ntinterp = 1,
                            writewft = 0,
                            tprint = 10.,
                            x_int = 2.,
                            precision = np.bytes_('d')
                            ):
        """
        Initializes default inputs for running 1D-TDSE with custom parameters 
        within Python API.

        Parameters:
        -----------
        Eguess: float, optional, default {-1.}
            Initial guess for the GS computation.
        num_r: int, optional, default {16000}
            Spatial grid resolution. ! Corresponds to Nx in the next version. !
        dx: float, optional, default {0.4}
            Spatial grid stepsize.
        InterpByDTorNT: int, optional, default {0}
            Interpolate by 'dt' (0) or number of points 'Ntinterp' (1).
        dt: float, optional, default {0.25}
            Temporal step size.
        trg_a: float, optional, default {1.3677}
            Rare gas parameter: H {sqrt(2)}, He {0.6950}, Ne {0.8161}, Ar {1.1893}, Kr {1.3676}, Xe {1.6171} 
            [Dissertation thesis Jan Vabek, tab. 7.1]
        CV: float, optional, default {1e-25}
            Convergence value for the GS computation using resolvent
        gauge: int, optional, default {0}
            Selection of gauge, length (0), velocity (1) <--- NOT IMPLEMENTED YET
        Ntinterp: int, optional, default {1}
            Number of points for the interpolation.
        writewft: int, optional, default {0}
            Store the wavefunction during the propagation (0 == No), (1 == Yes).
        tprint: float, optional, default {10.}
            Store the wavefunction every 'tprint' units of time (a.u.). 
            If 'tprint' is larger than half of the temporal grid, only the last wavefunction is returned.
        x_int: float, optional, default {2.}
            Integration limit for the ionization computation. 
        """
        
        self.Eguess = c_double(Eguess)
        self.num_r = c_int(num_r)
        self.dx = c_double(dx)
        self.InterpByDTorNT = c_int(InterpByDTorNT)
        self.dt = c_double(dt)
        self.Ntinterp = c_int(Ntinterp)
        self.analy.writewft = c_int(writewft)
        self.analy.tprint = c_double(tprint)
        self.x_int = c_double(x_int)
        self.trg.a = c_double(trg_a)
        self.CV = c_double(CV)
        self.gauge = c_int(gauge)
        self.precision = precision
        
    def init_prints(self, path_to_DLL):
        """
        Sets all prints to HDF5 to 1.
        """
        DLL = CDLL(path_to_DLL)
        set_prints = DLL.Set_all_prints
        set_prints.restype = output_print_def
        self.Print = set_prints()

    def init_time_and_field(self, DLL, filename = "", z_i = 0, r_i = 0, E = None, t = None):
        """
        Initializes field and temporal grid from custom arrays or from an hdf5 archive.

        Parameters:
        -----------
        DLL: TDSE_DLL class
            TDSE_DLL instance.
        filename: str, optional, default {""}
            HDF5 filename.
        z_i: int, optional, default {0}
            Index along z-axis in CUPRAD field.
        r_i: int, optional, default {0}
            Index along r-axis in CUPRAD field.
        E: optional, default {None}
            Electric field array.
        t: optional, default {None}
            Time array.
        """
        if (filename != "") and (E is None or t is None):
            f = h5py.File(filename, "r")
            field_shape = f[MMA.paths["CUPRAD_outputs"]+"/output_field"].shape
            if (z_i < 0) or (z_i >= field_shape[0]):
                print("Incorrect z-grid dimension selection. Select z in range (0, {})".format(field_shape[0]-1))
                f.close()
                return
            if (r_i < 0) or (r_i >= field_shape[2]):
                print("Incorrect r-grid dimension selection. Select r in range (0, {})".format(field_shape[2]-1))
                f.close()
                return        
            ### Load tgrid
            tgrid = f[MMA.paths["CUPRAD_outputs"]+"/tgrid"][()]/TIMEau
            ### Load field and convert to a.u.
            field = f[MMA.paths["CUPRAD_outputs"]+"/output_field"][z_i, r_i, :][()]/EFIELDau
            
            Nt = len(tgrid)
            self.Efield.Nt = Nt
            ### Init temporal grid
            DLL.set_time_and_field(self.ptr, tgrid, field, Nt)
            #self.Efield.tgrid = ctypes_arr_ptr(c_double, Nt, tgrid)
            #self.Efield.Field = ctypes_arr_ptr(c_double, Nt, field)
            f.close()

        else:
            Nt = len(t)
            assert(Nt == len(E))
            self.Efield.Nt = Nt
            DLL.set_time_and_field(self.ptr, t, E, Nt)
            #DLL.set_arr_vals(self.Efield.tgrid, ctypes_arr_ptr(c_double, Nt, t), Nt)
            #DLL.set_arr_vals(self.Efield.Field, ctypes_arr_ptr(c_double, Nt, E), Nt)
            #self.Efield.tgrid = ctypes_arr_ptr(c_double, Nt, t)
            #self.Efield.Field = ctypes_arr_ptr(c_double, Nt, E)

    def save_to_hdf5(self, filename):
        """
        Saves all available inputs from the `inputs_def` class into an HDF5 file. 
        Unavailable or unallocated inputs are neglected (e.g. `Field`, `tgrid`, `psi0`, ..).

        Parameters:
        -----------
        filename: str
            Name of the HDF5 archive for writing. 
        """
        f = h5py.File(filename, "a")
        try:
            f.create_group('TDSE_inputs')
        except ValueError:
            pass
        
        ### Write default inputs
        try:
            f.create_dataset("TDSE_inputs/trg_a", dtype="f", data=self.trg.a)
            f.create_dataset("TDSE_inputs/dx", dtype="f", data=self.dx)
            f.create_dataset("TDSE_inputs/dt", dtype="f", data=self.dt)
            f.create_dataset("TDSE_inputs/Eguess", dtype="f", data=self.Eguess)
            f.create_dataset("TDSE_inputs/N_r_grid", dtype="i", data=self.num_r)
            f.create_dataset("TDSE_inputs/gauge_type", dtype="i", data=self.gauge)
            f.create_dataset("TDSE_inputs/Ntinterp", dtype="i", data=self.Ntinterp)
            f.create_dataset("TDSE_inputs/InterpByDTorNT", dtype="i", data=self.InterpByDTorNT)
            f.create_dataset("TDSE_inputs/analy_writewft", dtype="i", data=self.analy.writewft)
            f.create_dataset("TDSE_inputs/analy_tprint", dtype="f", data=self.analy.tprint)
            f.create_dataset("TDSE_inputs/CV_criterion_of_GS", dtype="f", data=self.CV)
            f.create_dataset("TDSE_inputs/x_int", dtype="f", data=self.x_int)
            f.create_dataset("TDSE_inputs/num_t", dtype="i", data=self.num_t)
        except ValueError:
            pass

        ### Write field and time grid
        if self.Efield.Nt != 0:
            try:
                f.create_dataset("TDSE_inputs/Nt", dtype="i", data=self.Efield.Nt)
                f.create_dataset("TDSE_inputs/Efield", dtype="f", data=ctype_arr_to_numpy(self.Efield.Field, self.Efield.Nt))
                f.create_dataset("TDSE_inputs/tgrid", dtype="f", data=ctype_arr_to_numpy(self.Efield.tgrid, self.Efield.Nt))
            except ValueError:
                pass


        ### Write ground state, GS energy and x grid
        if self.Einit != 0.:
            try:
                f.create_dataset("TDSE_inputs/psi0", dtype="f", data=np.array([self.get_GS().real, self.get_GS().imag]).transpose())
                f.create_dataset("TDSE_inputs/x_grid", dtype="f", data=self.get_xgrid())
                f.create_dataset("TDSE_inputs/Einit", dtype="f", data=self.Einit)
            except ValueError:
                pass

        f.close()

    def get_xgrid(self):
        """
        Returns spatial grid.
        """
        return ctype_arr_to_numpy(self.x, self.num_r+1)
    
    def get_GS(self):
        """
        Returns ground state.
        """
        return ctype_cmplx_arr_to_numpy(self.psi0, self.num_r+1)
    
    def get_tgrid(self):
        """
        Returns temporal grid.
        """
        return ctype_arr_to_numpy(self.Efield.tgrid, self.Efield.Nt)
    
    def get_Efield(self):
        """
        Returns electric field.
        """
        return ctype_arr_to_numpy(self.Efield.Field, self.Efield.Nt)
        
    def delete(self, DLL):
        """
        Frees structure memory.

        Warning: if called twice on the same input, the kernel crashes.

        Parameters:
        -----------
        DLL: TDSE_DLL
            TDSE dynamic library.
        """
        DLL.free_inputs(self.ptr)
            
class outputs_def(Structure):
    """
    Output structure

    Attributes:
    -----------
    tgrid:
        Temporal grid.
    Efield:
        Electric field.
    sourceterm:
        Source term for Maxwell eqs.: <-grad V> - E term
    omegagrid:
        Frequency grid.
    FEfield:
        Field spectrum.
    FEfield_data:
        Field frequencies (positive frequencies only).
    FSourceterm:
        Source term spectrum (positive frequencies only).
    FSourceterm_data:
        Source term frequencies (positive frequencies only).
    FEfieldM2:
        Modulus squared of field spectrum (positive frequencies only)
    FsourcetermM2:
        Modulus squared of source term spectrum (positive frequencies only)
    PopTot:
        Total population of the ground state.
    PopInt:
        Ionization probability.
    expval:
        Expectation value of electron position.
    Nt:
        Temporal resolution.
    Nomega:
        Frequency grid resolution.
    psi:
        Wavefunction.
    
    """
    def __init__(self, *args: Any, **kw: Any):
        super().__init__(*args, **kw)
        self.ptr = byref(self)

    _fields_ = [
        ("tgrid", POINTER(c_double)),
        ("Efield", POINTER(c_double)),
        ("sourceterm", POINTER(c_double)),
        ("omegagrid", POINTER(c_double)),
        ("FEfield", POINTER(c_double)),
        ("Fsourceterm", POINTER(c_double)),
        ("FEfieldM2", POINTER(c_double)),
        ("FsourcetermM2", POINTER(c_double)),
        ("PopTot", POINTER(c_double)),
        ("PopInt", POINTER(c_double)),
        ("expval", POINTER(c_double)),
        ("Nt", c_int),
        ("Nomega", c_int),
        ("psi", POINTER(POINTER(c_double)))
    ]

    def save_to_hdf5(self, filename, inputs = None):
        """
        Saves the outputs into an HDF5 file. 

        To enable saving the wavefunction, user must supply the `inputs_def` 
        class as `inputs` argument, `inputs_def.analy.write_wft = c_int(1)` must be set.

        Warning: the wavefunction (if available) is stored as real and imaginary 
        part separately within the HDF5 file - psi_re, psi_im.


        Parameters:
        -----------
        filename: str
            Name of the HDF5 archive for writing. 
        inputs: inputs_def, optional, default {None}
            If included, the wavefunction can be stored into the HDF5 archive. 
        """
        f = h5py.File(filename, "a")
        try:
            f.create_group('TDSE_outputs')
        except ValueError:
            pass
        
        ### Write outputs
        try:
            f.create_dataset("TDSE_outputs/Nomega", dtype="i", data=self.Nomega)
            f.create_dataset("TDSE_outputs/Nt", dtype="i", data=self.Nt)
            f.create_dataset("TDSE_outputs/tgrid", dtype="f", data=self.get_tgrid())
            f.create_dataset("TDSE_outputs/Efield", dtype="f", data=self.get_Efield())
            f.create_dataset("TDSE_outputs/sourceterm", dtype="f", data=self.get_sourceterm())
            f.create_dataset("TDSE_outputs/omegagrid", dtype="f", data=self.get_omegagrid())
            f.create_dataset("TDSE_outputs/FEfield", dtype="f", data=np.array([self.get_FEfield().real, self.get_FEfield().imag]).transpose())
            f.create_dataset("TDSE_outputs/Fsourceterm", dtype="f", data=np.array([self.get_Fsourceterm().real, self.get_Fsourceterm().imag]).transpose())
            f.create_dataset("TDSE_outputs/PopTot", dtype="f", data=self.get_PopTot())
            f.create_dataset("TDSE_outputs/PopInt", dtype="f", data=self.get_PopInt())
            f.create_dataset("TDSE_outputs/expval", dtype="f", data=self.get_expval())
        except ValueError:
            pass
        
        ### Write wavefunction
        if inputs != None:
            try:
                wf = self.get_wavefunction(inputs, grids=False)
                wf_re = wf.real
                wf_im = wf.imag
                f.create_dataset("TDSE_outputs/psi_re", dtype="f", data = wf_re)
                f.create_dataset("TDSE_outputs/psi_im", dtype="f", data = wf_im)
            except:
                pass

        f.close()

    def load_from_hdf5(self, filename):
        """
        Loads the contents of the HDF5 archive into the `outputs_def` class.

        Parameters:
        -----------
        filename: str
            Name of the HDF5 archive for loading.         
        """
        with h5py.File(filename, "r") as f:
            try:
                self.Nomega = c_int(f["TDSE_outputs/Nomega"][()])
                self.Nt = c_int(f["TDSE_outputs/Nt"][()])
                self.tgrid = ctypes_arr_ptr(c_double, self.Nt, f["TDSE_outputs/tgrid"][()])
                self.Efield = ctypes_arr_ptr(c_double, self.Nt, f["TDSE_outputs/Efield"][()])
                self.sourceterm = ctypes_arr_ptr(c_double, self.Nt, f["TDSE_outputs/sourceterm"][()])
                self.expval = ctypes_arr_ptr(c_double, self.Nt, f["TDSE_outputs/expval"][()])
                self.omegagrid = ctypes_arr_ptr(c_double, self.Nomega, f["TDSE_outputs/omegagrid"][()])
                self.FEfield = ctypes_arr_ptr(c_double, 2*self.Nomega, np.array(f["TDSE_outputs/FEfield"][()]).flatten())
                self.Fsourceterm = ctypes_arr_ptr(c_double, 2*self.Nomega, np.array(f["TDSE_outputs/Fsourceterm"][()]).flatten())
                self.PopTot = ctypes_arr_ptr(c_double, self.Nt, f["TDSE_outputs/PopTot"][()])
                self.PopInt = ctypes_arr_ptr(c_double, self.Nt, f["TDSE_outputs/PopInt"][()])
            except KeyError:
                print("No output data stored in the HDF5 file.")
                return
            
            try:
                psi_re = f["TDSE_outputs/psi_re"][()]
                psi_im = f["TDSE_outputs/psi_im"][()]

                psi = np.array([np.array([[psi_r, psi_i] for psi_r, psi_i in 
                                          zip(psi_re[i], psi_im[i])]).flatten() 
                                          for i in range(len(psi_re))])
                
                self.psi = ctypes_mtrx_ptr(c_double, psi.shape, psi)
            except KeyError:
                pass
                


    def get_wavefunction(self, inputs, grids = True):
        """
        Returns complex ND-array storing the wavefunction in time. 

        Parameters:
        -----------
        inputs: inputs_def
            Input structure with inputs corresponding to outputs.
        grids: bool, optional, default {True}
            If true, returns the time and space grids
        
        Returns:
        --------
        tuple: tgrid, x, wavefunction (if grids == True)
        """
        if inputs.analy.writewft == 0:
            print("No wavefunction is stored!")
            return
        
        ### Number of steps per dt for printing in the temporal grid
        steps_per_dt = np.floor(inputs.analy.tprint/(self.tgrid[1]-self.tgrid[0]))
        ### Number of wavefunctions in the final grid
        size = int(self.Nt/steps_per_dt)

        wavefunction = get_wavefunction(self.psi, size, inputs.num_r+1)
        
        if grids:
            if size == 1:
                tgrid = np.array([self.tgrid[self.Nt-1]])
            else:
                t = self.get_tgrid()
                tgrid = np.linspace(t[0], t[-1], size)
            #x = np.linspace(-inputs.dx*inputs.num_r/2, inputs.dx*inputs.num_r/2, inputs.num_r+1)
            x = inputs.get_xgrid()
            return tgrid, x, wavefunction
        
        return wavefunction
    
    def get_tgrid(self):
        """
        Returns temporal grid.
        """
        return ctype_arr_to_numpy(self.tgrid, self.Nt)
    
    def get_Efield(self):
        """
        Returns electric field.
        """
        return ctype_arr_to_numpy(self.Efield, self.Nt)
    
    def get_sourceterm(self):
        """
        Returns source term <grad V> + E.
        """
        return ctype_arr_to_numpy(self.sourceterm, self.Nt)
    
    def get_PopTot(self):
        """
        Returns population of the ground state.
        """
        return ctype_arr_to_numpy(self.PopTot, self.Nt)
    
    def get_PopInt(self):
        """
        Returns integrated population.
        """
        return ctype_arr_to_numpy(self.PopInt, self.Nt)
    
    def get_Fsourceterm(self):
        """
        Returns source term <grad V> + E spectrum.
        """
        return ctype_cmplx_arr_to_numpy(self.Fsourceterm, self.Nomega)
    
    def get_expval(self):
        """
        Returns expectation value of x
        """
        return ctype_arr_to_numpy(self.expval, self.Nt)
    
    def get_omegagrid(self):
        """
        Returns omega grid.
        """
        return ctype_arr_to_numpy(self.omegagrid, self.Nomega)
    
    def get_FEfield(self):
        """
        Returns electric field positive spectrum.
        """
        return ctype_cmplx_arr_to_numpy(self.FEfield, self.Nomega)

    def delete(self, DLL):
        """
        Frees structure memory.

        Warning: if called twice on the same input, the kernel crashes.

        Parameters:
        -----------
        DLL: TDSE_DLL
            TDSE dynamic library.
        """
        DLL.free_outputs(self.ptr)