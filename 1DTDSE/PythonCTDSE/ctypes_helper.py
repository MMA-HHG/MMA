"""
ctypes helper
=============

This module contains helper functions for arrays and matrices manipulation using
Python and C.
"""

from ctypes import *
import numpy as np


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
