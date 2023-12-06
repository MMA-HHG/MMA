"""
Python TDSE
===========
"""

from ctypes import *
from typing import Any
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

### Return ctypes array to pointer
def ctypes_arr_ptr(ctype_, size_, arr_):
    """
    Creates ctypes array from a numpy array

    Parameters:
    -----------
    ctype_:
        Data type of the array.
    size_: int
        Size of the ctypes array.
    arr_:
        1D list or numpy array containing the data.

    Returns:
    --------
    ctypes pointer to an array
    """
    arr = ctype_ * size_
    return arr(*arr_)

def ctypes_mtrx_ptr(ctype_, size_, mtrx_):
    """
    Creates ctypes array from a numpy array

    Parameters:
    -----------
    ctype_:
        Data type of the matrix.
    size_: tuple
        Size tuple of the ctypes matrix.
    mtrx_:
        2D list or numpy array containing the data.

    Returns:
    --------
    ctypes pointer to an array
    """
    row_ptr = (POINTER(c_double)*size_[0])()
    arr_ptr = c_double*size_[1]
    
    for i in range(size_[0]):
        row_ptr[i] = arr_ptr()
        for j in range(size_[1]):
            row_ptr[i][j] = mtrx_[i][j]

    return row_ptr

### return numpy array from c_types array
def ctype_arr_to_numpy(c_arr, size):
    """
    Returns numpy array from a Ctypes array.

    Parameters:
    -----------
    c_arr: 
        Pointer to C array.
    size: int
        Array size.

    Returns:
    --------
    numpy array with the C array data
    """
    return np.array([c_arr[i] for i in range(size)])

### return numpy array from c_types array of complex numbers
def ctype_cmplx_arr_to_numpy(c_arr, size):
    """
    Returns complex numpy array from a Ctypes array. Data in complex ctype array
    is stored as: z = c_arr[2*j] + i * c_arr[2*j+1], j = 0, .., size-1.

    Parameters:
    -----------
    c_arr: 
        Pointer to C array.
    size: int
        Array size.

    Returns:
    --------
    numpy array with the C array data
    """
    return np.array([c_arr[2*i]+1j*c_arr[2*i+1] for i in range(size)])

### Return numpy array from c_types matrix of wavefunction values
def get_wavefunction(c_arr, timesteps, wf_size):
    """
    Returns complex numpy ND-array wavefunction from a Ctypes matrix. Data in complex ctype array
    is stored as: z = c_arr[2*j] + i * c_arr[2*j+1], j = 0, .., size-1.

    Parameters:
    -----------
    c_arr: 
        Pointer to C array.
    timesteps: int
        Number of stored wavefunctions in time.
    wf_size: int
        Size of a single wavefunction.

    Returns:
    --------
    numpy ND-array with the C array data
    """
    return np.array([[c_arr[i][2*j] + 1j*c_arr[i][2*j+1] for j in range(wf_size)] for i in range(timesteps)])

def ctype_mtrx_to_numpy(c_arr, N_rows, N_cols):
    """
    Returns numpy ND-array from a Ctypes matrix. 

    Parameters:
    -----------
    c_arr: 
        Pointer to C array.
    N_rows: int
        Number of rows.
    N_cols: int
        Number of columns.

    Returns:
    --------
    numpy ND-array with the C array data
    """
    return np.array([[c_arr[i][j] for j in range(N_cols)] for i in range(N_rows)])


### Physical constants
Ip_HeV = 27.21138602 
"""Ionisation potential"""
hbar = 1.054571800e-34
"""Planck constant"""
alpha_fine = 1/137.035999139
"""Fine structure constant"""
c_light = 299792458.
"""Speed of light in vacuum"""
elcharge = 1.602176565e-19
"""Electron charge"""
elmass = 9.10938356e-31
"""Electron mass"""
mu0 = 4.0*np.pi*1e-7
"""Magnetic permeability"""
eps0 = 1.0/(mu0*c_light*c_light)
"""Vacuum permitivity"""
r_Bohr = 4.0*np.pi*eps0*hbar*hbar/(elmass*elcharge*elcharge)
"""Bohr radius"""
TIMEau = (elmass*r_Bohr*r_Bohr)/hbar
"""Atomic unit of time"""
EFIELDau = hbar*hbar/(elmass*r_Bohr*r_Bohr*r_Bohr*elcharge)
"""Atomic unit of field intensity"""
k_Boltz = 1.38064852e-23
"""Boltzmann constant"""
absolute_zero = -273.15
"""Absolute zero temperature"""
torr2SI = 101325./760.
"""Torr to SI conversion factor"""
    
def plot_colormap(
        y_axis, 
        x_axis, 
        z, 
        plot_scale = "log", 
        z_min = 1e-14,
        z_max = 1,
        figsize = (12, 3),
        y_label = "",
        x_label = "",
        cmap = "jet"
    ):
    """
    Plots the colormap of arbitrary matrix for arbitrary distribution.
    
    Parameters:
    -----------
    y_axis (numpy.ndarray):
        Values on the y-axis on the colormap.
    x_axis (numpy.ndarray):
        Values on the y-axis on the colormap.
    z (numpy.ndarray):
        Matrix for plotting.
    plot_scale: {'log', 'linear'}, optional, default: 'log'
        Plot with logarithmic or linear scale.
    z_min (float):
        Minimum plotting value on the colormesh.
    z_max (float):
        Maximum threshold for the colorbar
    figsize (tuple), optional, default: (12, 3)
        Set figure size in inches - (width, height)
    y_label (str), optional, default: ""
        Label of the y-axis
    x_label (str), optional, default: ""
        Label of the x-axis
    cmap (str), optional, default: "jet"
        Colormap name from the Matplotlib colormap library

    Returns:
    --------
    None
    """
    fig, ax = plt.subplots()
    
    if plot_scale == 'log':
        col = LogNorm(vmin = z_min, vmax = z_max)
    elif plot_scale == 'linear':
        col = Normalize(vmin = z_min, vmax = z_max)
    else:
        raise ValueError("Unknown scale '" + plot_scale +"'. "
                         "Available options are 'log' and 'linear'")
        
    c = ax.pcolormesh(x_axis, y_axis, z,
        cmap = cmap,
        shading = 'gouraud',
        norm = col
    )

    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)

    fig.colorbar(c, ax=ax)
    fig.dpi = 300
    fig.set_size_inches(figsize)

    plt.show()

def plot(
        *args,
        plot_scale = "linear", 
        figsize = (5, 3),
        y_label = "",
        x_label = "",
        **kwargs
        ):
    """
    Simple plot method.

    Parameters:
    -----------
    *args
        Arguments for the plot command, i.e. x, y or x1, y1, x2, y2, .. .
    plot_scale: {'log', 'linear'}, optional, default: 'linear'
        Plot with logarithmic or linear scale.
    figsize (tuple), optional, default: (5, 3)
        Set figure size in inches - (width, height)
    y_label (str), optional, default: ""
        Label of the y-axis
    x_label (str), optional, default: ""
        Label of the x-axis
    **kwargs
        Keyword arguments for the plot command.

    Returns:
    --------
    None
    """
    
    fig, ax = plt.subplots()

    if plot_scale == 'log':
        plt.semilogy(*args, **kwargs)
    elif plot_scale == 'linear':
        plt.plot(*args, **kwargs)
    else:
        raise ValueError("Unknown scale '" + plot_scale +"'. "
                         "Available options are 'log' and 'linear'")
    fig.dpi = 300
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    fig.set_size_inches(figsize)
    
    plt.show()

### Define structures
class sin2_definition(Structure):
    _fields_ = [
        ("oc", c_double),
        ("o", c_double),
        ("A0", c_double),
        ("nc", c_double),
        ("phi", c_double),
        ("phi0", c_double),
        ("ti", c_double),
        ("E0", c_double)
    ]

class trap_definition(Structure):
    _fields_ = [
        ("omega", c_double),
        ("E0", c_double),
        ("phi", c_double),
        ("ton", c_double),
        ("toff", c_double),
        ("nc", c_int)
    ]


class Esin2_definition(Structure):
    _fields_ = [
        ("oc", c_double),
        ("o", c_double),
        ("A0", c_double),
        ("nc", c_double),
        ("phi", c_double),
        ("phi0", c_double),
        ("ti", c_double),
        ("E0", c_double)
    ]

class flattop1_definition(Structure):
    _fields_ = [
        ("ton", c_double),
        ("toff", c_double),
        ("T", c_double),
        ("o", c_double),
        ("phi", c_double),
        ("A", c_double),
        ("ti", c_double)
    ]

class flattop1chirp_definition(Structure):
    _fields_ = [
        ("ton", c_double),
        ("toff", c_double),
        ("T", c_double),
        ("o", c_double),
        ("phi", c_double),
        ("A", c_double),
        ("ti", c_double),
        ("b", c_double),
        ("c", c_double)
    ]
    
### Field structure
class Efield_var(Structure):
    _fields_ = [
        ("fieldtype", c_int),
        ("fieldtype2", c_int),
        ("tgrid", POINTER(c_double)),
        ("Field", POINTER(c_double)),
        ("dt", c_double),
        ("omega", c_double),
        ("E0", c_double),
        ("phi", c_double),
        ("ton", c_double),
        ("toff", c_double),
        ("Nt", c_int),
        ("Nftlt1", c_int),
        ("Nftlt1ch", c_int),
        ("Nsin2", c_int),
        ("NEsin2", c_int),
        ("nc", c_int),
        ("trap", trap_definition),
        ("sin2", POINTER(sin2_definition)),
        ("Esin2", POINTER(Esin2_definition)),
        ("flt1", POINTER(flattop1_definition)),
        ("flt1ch", POINTER(flattop1chirp_definition)),
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

    def init_inputs(self, filename):
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
                            tprint = 10,
                            x_int = 2.,
                            precision = np.string_('d')
                            ):
        """
        Initializes default inputs for running 1D-TDSE with custom parameters 
        within Python API.

        Parameters:
        -----------
        Eguess: float, optional, default {-1.}
            Initial guess for the GS computation.
        num_r: int, optional, default {16000}
            Spatial grid resolution.
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
        tprint: float, optional, default {10}
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
        #precision = precision
        self.precision = precision
        
    def init_prints(self, path_to_DLL):
        """
        Sets all prints to HDF5 to 1.
        """
        DLL = CDLL(path_to_DLL)
        set_prints = DLL.Set_all_prints
        set_prints.restype = output_print_def
        self.Print = set_prints()

    def init_time_and_field(self, filename = "", z_i = 0, r_i = 0, E = None, t = None):
        """
        Initializes field and temporal grid from custom arrays or from an hdf5 archive.

        Parameters:
        -----------
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
            field_shape = f["outputs/output_field"].shape
            if (z_i < 0) or (z_i >= field_shape[0]):
                print("Incorrect z-grid dimension selection. Select z in range (0, {})".format(field_shape[0]-1))
                f.close()
                return
            if (r_i < 0) or (r_i >= field_shape[2]):
                print("Incorrect r-grid dimension selection. Select r in range (0, {})".format(field_shape[2]-1))
                f.close()
                return        
            ### Load tgrid
            tgrid = f["outputs/tgrid"][()]/TIMEau
            ### Load field and convert to a.u.
            field = f["outputs/output_field"][z_i, :, r_i][()]/EFIELDau
            
            Nt = len(tgrid)
            self.Efield.Nt = Nt
            ### Init temporal grid
            self.Efield.tgrid = ctypes_arr_ptr(c_double, Nt, tgrid)
            self.Efield.Field = ctypes_arr_ptr(c_double, Nt, field)
            f.close()

        else:
            Nt = len(t)
            assert(Nt == len(E))
            self.Efield.Nt = Nt
            self.Efield.tgrid = ctypes_arr_ptr(c_double, Nt, t)
            self.Efield.Field = ctypes_arr_ptr(c_double, Nt, E)
            ### Do not interpolate
            #self.InterpByDTorNT = c_int(1)

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

class TDSE_DLL:
    """
    Contains functions and routines from the CTDSE.

    Atrributes:
    -----------
    DLL: a CDLL instance

    Methods:
    --------
    * __init__ - class initialization
    * init_GS - initialization of the ground state
    * call1DTDSE - propagates the ground state according the field
    * compute_PES - computes photoelectron spectrum of wavefunction
    * gabor_transform - computes Gabor transform of a signal
    * free_mtrx - frees double matrix in C
    * free_outputs - frees output structure in C
    * free_inputs - frees input structure in C
    """
    def __init__(self, path_to_DLL):
        """
        Class initialization.

        Parameters:
        -----------
        path_to_DLL: str
            Path to the dynamic library.
        """
        self.DLL = CDLL(path_to_DLL)

    def init_GS(self, inputs):
        """
        Initialization of the ground state.

        Parameters:
        -----------
        inputs: inputs_def
            Input structure
        """
        ### Find ground state and init grids
        init_grid = self.DLL.Initialise_grid_and_ground_state
        init_grid.restype = None
        init_grid.argtypes = [POINTER(inputs_def)]
        init_grid(inputs.ptr)

    def call1DTDSE(self, inputs, outputs):
        """
        Propagates the ground state according the field.

        Parameters:
        -----------
        inputs: 
            Ctypes inputs structure
        outputs: 
            Ctypes outputs structure
        """
        ### Do the propagation
        TDSE = self.DLL.call1DTDSE
        TDSE.restype = None
        TDSE.argtypes = [POINTER(inputs_def), POINTER(outputs_def)]
        TDSE(inputs.ptr, outputs.ptr)

    def compute_PES(self, inputs, psi, E_start = -0.6, num_E = 10000, dE = 5e-4, Estep = 5e-4):
        """
        Computes photoelectron spectrum (PES) from a wavefunction.

        Parameters:
        -----------
        inputs:
            Ctypes input structure instance
        psi:
            Ctypes pointer to C array with wavefunction
        E_start: float, optional, default E_start = -0.6
            Initial energy for the PES
        num_E: int, optional, default num_E = 10000
            Energy grid resolution
        dE: float, optional, default dE = 5e-4
            Integration energy step
        E_step: float, optional, default E_step = 5e-4
            Energy step for the PES
        """
        PES = self.DLL.window_analysis
        PES.restype = POINTER(c_double)
        PES.argtypes = [inputs_def, POINTER(c_double), c_int, c_double, c_double, c_double]
        res = PES(inputs, psi, c_int(num_E), c_double(dE), c_double(Estep), c_double(E_start))
        E_grid = np.linspace(E_start, E_start+(num_E-1)*Estep, num_E)
        res_np = ctype_arr_to_numpy(res, num_E)
        self.free_arr(res)
        return E_grid, res_np

    def gabor_transform(self, signal, dt, N, omega_max, t_min, t_max, N_G, a = 8.):
        """
        Computes fast Gabor transform of a signal.

        Parameters:
        -----------
        signal: ctypes double array pointer, list or numpy array
            Signal for the Gabor transform.
        dt: float
            Time step of the signal.
        N: int
            Size of the signal.
        omega_max: float
            Maximum frequency for the transform in a.u..
        t_min: float
            Minimum time for Gabor.
        t_max: float
            Maximum time for Gabor.
        N_G: int
            Resolution for Gabor.
        a: float, optional, default a = 8.
            Width of the Gabor window, in a.u..
        """
        gabor = self.DLL.GaborTransform
        gabor.restype = POINTER(POINTER(c_double))
        gabor.argtypes = [POINTER(c_double), c_double, c_int, c_int, c_int, c_double, c_double, c_double]
        omegas = 2*np.pi*np.fft.fftfreq(N, dt)[0:N//2]
        omega_range = (omegas <= omega_max)
        N_freq = len(omegas[omega_range])
        if ((type(signal) == list) or (type(signal) == np.ndarray)):
            signal = ctypes_arr_ptr(c_double, len(signal), signal)
        res = gabor(signal, c_double(dt), c_int(N), c_int(N_freq), c_int(N_G), 
                    c_double(t_min), c_double(t_max), c_double(a))
        gabor_res = ctype_mtrx_to_numpy(res, N_G, N_freq)
        self.free_mtrx(res, N_G)
        return np.linspace(t_min, t_max, N_G), omegas[omega_range], np.transpose(gabor_res)

    def free_mtrx(self, buffer_ptr, N_rows):
        """
        Frees 2-D C array.

        Parameters:
        -----------
        buffer_ptr: ctypes pointer instance
            Pointer to the buffer to be freed
        N_rows: int
            Number of rows in the 2-D array
        """
        self.DLL.free_mtrx.restype = None
        self.DLL.free_mtrx.argtypes = [POINTER(POINTER(c_double)), c_int]
        self.DLL.free_mtrx(buffer_ptr, c_int(N_rows))

    def free_arr(self, buffer_ptr):
        """
        Frees 1-D C array.

        Parameters:
        -----------
        buffer_ptr: ctypes pointer instance
            Pointer to the buffer to be freed
        """
        self.DLL.free_arr.restype = None
        self.DLL.free_arr.argtypes = [POINTER(c_double)]
        self.DLL.free_arr(buffer_ptr)

    def free_outputs(self, out_ptr):
        """
        Frees outputs structure.

        Warning: if called twice on the same input, the kernel crashes.

        Parameters:
        -----------
        out_ptr: ctypes pointer instance
            Pointer to the structure to be freed
        """
        self.DLL.outputs_destructor.restype = None
        self.DLL.outputs_destructor.argtypes = [POINTER(outputs_def)]
        self.DLL.outputs_destructor(out_ptr)

    def free_inputs(self, in_ptr):
        """
        Frees inputs structure.

        Warning: if called twice on the same input, the kernel crashes.

        Parameters:
        -----------
        out_ptr: ctypes pointer instance
            Pointer to the structure to be freed
        """
        self.DLL.inputs_destructor.restype = None
        self.DLL.inputs_destructor.argtypes = [POINTER(inputs_def)]
        self.DLL.inputs_destructor(in_ptr)

