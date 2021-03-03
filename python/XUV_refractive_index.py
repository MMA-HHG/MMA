####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2021)
#
# This module creates the interpolation functions from tables loaded from the archive
# The functions are called by f1, f2 = XUV_refractive_index.getf(gas,Energy), Energy is in eV
# alternatively, the full handle is XUV_refractive_index.index_funct[gas][f1/f2](Energy), Energy is in eV


import numpy as np
from scipy import interpolate
import h5py
import sys
import os
import shutil

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


source_archive = os.path.join(THIS_DIR, 'XUV_refractive_index_tables.h5')
index_funct = {}

with h5py.File(source_archive, 'r') as SourceFile: # access option http://docs.h5py.org/en/stable/high/file.html#file
    gases = list(SourceFile.keys())
    index_table = {}
    print(gases)
    for gas in gases:
        local_table ={
            'Energy': SourceFile[gas]['Energy'][:],
            'f1': SourceFile[gas]['f1'][:],
            'f2': SourceFile[gas]['f2'][:]
        }
        index_table.update({gas: local_table})
        local_table ={
            'f1': interpolate.interp1d(SourceFile[gas]['Energy'][:],SourceFile[gas]['f1'][:]),
            'f2': interpolate.interp1d(SourceFile[gas]['Energy'][:], SourceFile[gas]['f2'][:])
        }
        index_funct.update({gas: local_table})

def getf(g,E):
    return index_funct[g]['f1'](E)[()], index_funct[g]['f2'](E)[()]