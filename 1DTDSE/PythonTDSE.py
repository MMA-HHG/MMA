"""
Python TDSE
===========
"""

from ctypes import *
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

### Return ctypes array to pointer
def ctypes_arr_ptr(ctype_, size_, arr_):
    arr = ctype_ * size_
    return arr(*arr_)

### return numpy array from c_types array
def ctype_arr_to_numpy(c_arr, size):
    return np.array([c_arr[i] for i in range(size)])

### return numpy array from c_types array of complex numbers
def ctype_cmplx_arr_to_numpy(c_arr, size):
    return np.array([c_arr[2*i]+1j*c_arr[2*i+1] for i in range(size)])

### Return numpy array from c_types matrix of wavefunction values
def get_wavefunction(c_arr, timesteps, wf_size):
    return np.array([[c_arr[i][2*j] + 1j*c_arr[i][2*j+1] for j in range(wf_size)] for i in range(timesteps)])

def ctype_mtrx_to_numpy(c_arr, N_rows, N_cols):
    return np.array([[c_arr[i][j] for j in range(N_cols)] for i in range(N_rows)])


### Physical constants
Ip_HeV = 27.21138602
hbar = 1.054571800e-34
alpha_fine = 1/137.035999139
c_light = 299792458.
elcharge = 1.602176565e-19
elmass = 9.10938356e-31
mu0 = 4.0*np.pi*1e-7
eps0 = 1.0/(mu0*c_light*c_light)
r_Bohr = 4.0*np.pi*eps0*hbar*hbar/(elmass*elcharge*elcharge)
TIMEau = (elmass*r_Bohr*r_Bohr)/hbar
EFIELDau = hbar*hbar/(elmass*r_Bohr*r_Bohr*r_Bohr*elcharge)
k_Boltz = 1.38064852e-23
absolute_zero = -273.15
torr2SI = 101325./760.
    
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
        x, 
        y, 
        plot_scale = "linear", 
        figsize = (5, 3),
        y_label = "",
        x_label = "",
        ):
    """
    Simple plot method.

    Parameters:
    -----------
    x (numpy.ndarray):
        x axis.
    y (numpy.ndarray):
        y axis.
    plot_scale: {'log', 'linear'}, optional, default: 'linear'
        Plot with logarithmic or linear scale.
    figsize (tuple), optional, default: (5, 3)
        Set figure size in inches - (width, height)
    y_label (str), optional, default: ""
        Label of the y-axis
    x_label (str), optional, default: ""
        Label of the x-axis

    Returns:
    --------
    None
    """
    
    fig, ax = plt.subplots()
    fig.dpi = 300
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    fig.set_size_inches(figsize)

    if plot_scale == 'log':
        plt.semilogy(x, y)
    elif plot_scale == 'linear':
        plt.plot(x, y)
    else:
        raise ValueError("Unknown scale '" + plot_scale +"'. "
                         "Available options are 'log' and 'linear'")
    
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
    """
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
        ("num_exp", c_int),
        ("dx", c_double),
        ("psi0", POINTER(c_double)),
        ("psi", POINTER(c_double)),
        ("x", POINTER(c_double)),
        ("ton", c_double),
        ("toff", c_double),
        ("timet", POINTER(c_double)),
        ("dipole", POINTER(c_double)),
        ("gauge", c_int),
        ("transformgauge", c_int),
        ("x_int", c_double),
        ("analy", analy_def),
        ("InterpByDTorNT", c_int),
        ("Ntinterp", c_int),
        ("PrintGaborAndSpectrum", c_int),
        ("PrintOutputMethod", c_int),
        ("textend", c_double),
        ("dtGabor", c_double),
        ("tmin1window", c_double),
        ("tmin2window", c_double),
        ("tmax1window", c_double),
        ("tmax2window", c_double),
        ("a_Gabor", c_double),
        ("omegaMaxGabor", c_double),
        ("Print", output_print_def),
        ("CV", c_double),
        ("precision", c_char * 2)
    ]

    def init_inputs(self, filename):
        with h5py.File(filename, "r") as f:
            self.Eguess = c_double(f["TDSE_inputs/Eguess"][()])
            self.num_r = c_int(f["TDSE_inputs/N_r_grid"][()])
            self.num_exp = c_int(f["TDSE_inputs/N_r_grid_exp"][()])
            self.dx = c_double(f["TDSE_inputs/dx"][()])
            self.InterpByDTorNT = c_int(f["TDSE_inputs/InterpByDTorNT"][()])
            self.dt = c_double(f["TDSE_inputs/dt"][()])
            self.Ntinterp = c_int(f["TDSE_inputs/Ntinterp"][()])
            self.textend = c_double(f["TDSE_inputs/textend"][()])
            self.analy.writewft = c_int(f["TDSE_inputs/analy_writewft"][()])
            self.analy.tprint = c_double(f["TDSE_inputs/analy_tprint"][()])
            self.x_int = c_double(f["TDSE_inputs/x_int"][()])
            self.PrintGaborAndSpectrum = c_int(f["TDSE_inputs/PrintGaborAndSpectrum"][()])
            self.a_Gabor = c_double(f["TDSE_inputs/a_Gabor"][()])
            self.omegaMaxGabor = c_double(f["TDSE_inputs/omegaMaxGabor"][()])
            self.dtGabor = c_double(f["TDSE_inputs/dtGabor"][()])
            self.tmin1window = c_double(f["TDSE_inputs/tmin1window"][()])
            self.tmax1window = c_double(f["TDSE_inputs/tmax1window"][()])
            self.tmin2window = c_double(f["TDSE_inputs/tmin2window"][()])
            self.tmax2window = c_double(f["TDSE_inputs/tmax2window"][()])
            self.PrintOutputMethod = c_int(f["TDSE_inputs/PrintOutputMethod"][()])
            self.trg.a = c_double(f["TDSE_inputs/trg_a"][()])
            self.CV = c_double(f["TDSE_inputs/CV_criterion_of_GS"][()])
            self.gauge = c_int(f["TDSE_inputs/gauge_type"][()])
            precision = 'd'
            self.precision = precision.encode('utf-8')

    def init_default_inputs(self,
                            Eguess = -1.,
                            num_r = 16000,
                            dx = 0.4,
                            InterpByDTorNT = 0,
                            dt = 0.25,
                            trg_a = 1.3677,
                            CV = 1e-25,
                            gauge = 0,
                            num_exp = 0,
                            Ntinterp = 1,
                            textend = 200.0,
                            writewft = 0,
                            tprint = 10,
                            x_int = 2.,
                            PrintOutputMethod = 1,
                            precision = np.string_('d')
                            ):
        
        self.Eguess = c_double(Eguess)
        self.num_r = c_int(num_r)
        self.num_exp = c_int(num_exp)
        self.dx = c_double(dx)
        self.InterpByDTorNT = c_int(InterpByDTorNT)
        self.dt = c_double(dt)
        self.Ntinterp = c_int(Ntinterp)
        self.textend = c_double(textend)
        self.analy.writewft = c_int(writewft)
        self.analy.tprint = c_double(tprint)
        self.x_int = c_double(x_int)
        self.PrintOutputMethod = c_int(PrintOutputMethod)
        self.trg.a = c_double(trg_a)
        self.CV = c_double(CV)
        self.gauge = c_int(gauge)
        #precision = precision
        self.precision = precision
        
    def init_prints(self, path_to_DLL):
        DLL = CDLL(path_to_DLL)
        set_prints = DLL.Set_all_prints
        set_prints.restype = output_print_def
        self.Print = set_prints()

    def init_time_and_field(self, filename = "", z_i = 0, r_i = 0, E = None, t = None):

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
            #t = c_double * Nt
            #t = t(*tgrid)
            #self.Efield.tgrid = t
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
            

class outputs_def(Structure):
    _fields_ = [
        ("tgrid", POINTER(c_double)),
        ("Efield", POINTER(c_double)),
        ("sourceterm", POINTER(c_double)),
        ("omegagrid", POINTER(c_double)),
        ("FEfield", POINTER(POINTER(c_double))),
        ("FEfield_data", POINTER(POINTER(c_double))),
        ("Fsourceterm", POINTER(POINTER(c_double))),
        ("Fsourceterm_data", POINTER(c_double)),
        ("FEfieldM2", POINTER(c_double)),
        ("FsourcetermM2", POINTER(c_double)),
        ("PopTot", POINTER(c_double)),
        ("sourcetermfiltered", POINTER(c_double)),
        ("PopInt", POINTER(c_double)),
        ("expval", POINTER(c_double)),
        ("Nt", c_int),
        ("Nomega", c_int),
        ("psi", POINTER(POINTER(c_double)))
    ]

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
        init_grid(byref(inputs))

    def call1DTDSE(self, in_ptr, out_ptr):
        """
        Propagates the ground state according the field.

        Parameters:
        -----------
        in_ptr: 
            Ctypes pointer instance of inputs structure
        out_ptr: 
            Ctypes pointer instance of outputs structure
        """
        ### Do the propagation
        TDSE = self.DLL.call1DTDSE
        TDSE.restype = None
        TDSE.argtypes = [POINTER(inputs_def), POINTER(outputs_def)]
        TDSE(in_ptr, out_ptr)

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
        E_grid = np.linspace(E_start, E_start+(num_E-1)*dE, num_E)
        return E_grid, ctype_arr_to_numpy(res, num_E)

    def gabor_transform(self, signal, dt, N, omega_max, t_min, t_max, N_t, a = 8.):
        """
        Computes fast Gabor transform of a signal.

        Parameters:
        -----------
        signal: ctypes double array pointer
            Signal for the Gabor transform
        dt: float
            Time step of the signal
        N: int
            Size of the signal
        omega_max: float
            Maximum frequency for the transform in a.u.
        t_min: float
            Minimum time for Gabor
        t_max: float
            Maximum time for Gabor
        N_t: int
            Resolution for Gabor
        a: float, optional, default a = 8.
            Width of the Gabor window, in a.u.
        """
        gabor = self.DLL.GaborTransform
        gabor.restype = POINTER(POINTER(c_double))
        gabor.argtypes = [POINTER(c_double), c_double, c_int, c_int, c_int, c_double, c_double, c_double]
        omegas = np.fft.fftfreq(N, dt)[0:N//2]
        omega_range = (omegas <= omega_max)
        N_freq = len(omegas[omega_range])
        res = gabor(signal, c_double(dt), c_int(N), c_int(N_freq), c_int(N_t), 
                    c_double(t_min), c_double(t_max), c_double(a))
        gabor_res = ctype_mtrx_to_numpy(res, N_t, N_freq)
        self.free_mtrx(res, N_t)
        return np.linspace(t_min, t_max, N_t), omegas[omega_range], np.transpose(gabor_res)

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

    def free_outputs(self, out_ptr):
        """
        Frees outputs structure.

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

        Parameters:
        -----------
        out_ptr: ctypes pointer instance
            Pointer to the structure to be freed
        """
        self.DLL.inputs_destructor.restype = None
        self.DLL.inputs_destructor.argtypes = [POINTER(inputs_def)]
        self.DLL.inputs_destructor(in_ptr)