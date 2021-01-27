####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2021)
#
# convert cxro tables (https://henke.lbl.gov/optical_constants/asf.html) to HDF5
import numpy as np
import h5py
import sys
import os
import shutil

gases = ['Ar', 'Kr']
filenames = {
    'Ar' : 'ar.nff',
    'Kr' : 'kr.nff'
}
tables_path = 'cxro_tables'
target_archive = 'XUV_refractive_index_tables.h5'


with h5py.File(target_archive, 'w') as GeneratedFile: # access option http://docs.h5py.org/en/stable/high/file.html#file
    for gas in gases:
        inputfilename = os.path.join(tables_path,filenames[gas])
        grp = GeneratedFile.create_group(gas)
        with open(inputfilename, "r") as InputFile:
            omegagrid = []
            f1 = []
            f2 = []
            lines = InputFile.readlines()
            First = True
            for line in lines:
                sep_line = line.split()
                if (First or (len(sep_line) == 0)):
                    First = False
                else:
                    omegagrid.append(sep_line[0])
                    f1.append(sep_line[1])
                    f2.append(sep_line[2])

        omegagrid = np.asarray(omegagrid, dtype='d')
        f1 = np.asarray(f1, dtype='d')
        f2 = np.asarray(f2, dtype='d')
        dset_id = grp.create_dataset('Energy', data=omegagrid)
        dset_id.attrs['units'] = np.string_('[eV]')
        dset_id = grp.create_dataset('f1', data=f1)
        dset_id.attrs['units'] = np.string_('[-]')
        dset_id = grp.create_dataset('f2', data=f2)
        dset_id.attrs['units'] = np.string_('[-]')

# print(omegagrid)


print('done')
