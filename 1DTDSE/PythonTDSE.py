"""
Python CTDSE
============

Python CTDSE provides the Python interface for executing the 1D-TDSE C code (CTDSE).
This module contains the functions for executing the CTDSE and the C-accelerated 
analysis functions.
"""

from PythonCTDSE.constants import *
from PythonCTDSE.ctypes_helper import *
from PythonCTDSE.plotting import *
from PythonCTDSE.structures import *
from typing import Any
import numpy as np


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

    def compute_PES(self, inputs, psi, E_start = -0.6, num_E = 10000, epsilon = 5e-4, Estep = 5e-4):
        """
        Computes photoelectron spectrum (PES) from a wavefunction.

        Parameters:
        -----------
        inputs:
            Ctypes input structure instance
        psi:
            Numpy array with wavefunction or ctypes pointer to C array with wavefunction
        E_start: float, optional, default E_start = -0.6
            Initial energy for the PES
        num_E: int, optional, default num_E = 10000
            Energy grid resolution
        epsilon: float, optional, default epsilon = 5e-4
            Window size for resolvent
        E_step: float, optional, default E_step = 5e-4
            Energy step for the PES
        """
        PES = self.DLL.window_analysis
        PES.restype = POINTER(c_double)
        PES.argtypes = [inputs_def, POINTER(c_double), c_int, c_double, c_double, c_double]
        if type(psi) is np.ndarray:
            psi = np.array([[x.real, x.imag] for x in psi]).flatten()
            res = PES(inputs, ctypes_arr_ptr(c_double, len(psi), psi), c_int(num_E), c_double(epsilon), c_double(Estep), c_double(E_start))
        else:
            res = PES(inputs, psi, c_int(num_E), c_double(epsilon), c_double(Estep), c_double(E_start))
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

    def set_time_and_field(self, in_ptr, time, field, N):
        """
        Sets time and fields from Numpy arrays using C allocation.

        Parameters:
        -----------
        in_ptr: ctypes pointer instance
            Pointer to the input structure
        time: ctypes pointer instance
            Pointer to the time array
        field: ctypes pointer instance
            Pointer to the field array
        """
        self.DLL.set_time_and_field.restype = None
        self.DLL.set_time_and_field.argtypes = [POINTER(inputs_def), POINTER(c_double), POINTER(c_double), c_int]
        time_ptr = ctypes_arr_ptr(c_double, N, time)
        field_ptr = ctypes_arr_ptr(c_double, N, field)
        self.DLL.set_time_and_field(in_ptr, time_ptr, field_ptr, N)
