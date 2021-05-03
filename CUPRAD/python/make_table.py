####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2021)
#
# convert cxro tables (https://henke.lbl.gov/optical_constants/asf.html) to HDF5
import numpy as np
import h5py
import sys
import os
import shutil

gases = ['Ar_Henke', 'Kr_Henke', 'Kr_NIST']
filenames = {
    'Ar_Henke' : ['ar.nff'],
    'Kr_Henke' : ['kr.nff'],
    'Kr_NIST' : ['kr_f1_NIST.txt', 'kr_f2_NIST.txt']
}
tables_path = 'XUV_tables'
target_archive = 'XUV_refractive_index_tables.h5'


with h5py.File(target_archive, 'w') as GeneratedFile: # access option http://docs.h5py.org/en/stable/high/file.html#file
    for gas in gases:
        if len(filenames[gas]) == 1:
            inputfilename = os.path.join(tables_path,filenames[gas][0])
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
            dset_id = grp.create_dataset('Energy_f1', data=omegagrid)
            dset_id.attrs['units'] = np.string_('[eV]')
            dset_id = grp.create_dataset('Energy_f2', data=omegagrid)
            dset_id.attrs['units'] = np.string_('[eV]')
            dset_id = grp.create_dataset('f1', data=f1)
            dset_id.attrs['units'] = np.string_('[-]')
            dset_id = grp.create_dataset('f2', data=f2)
            dset_id.attrs['units'] = np.string_('[-]')
            
        elif len(filenames[gas]) == 2:
            inputfilename1 = os.path.join(tables_path,filenames[gas][0])
            inputfilename2 = os.path.join(tables_path,filenames[gas][1])
            grp = GeneratedFile.create_group(gas)
            with open(inputfilename1, "r") as InputFile1, open(inputfilename2, "r") as InputFile2:
                omegagrid1 = []
                omegagrid2 = []
                f1 = []
                f2 = []
                
                lines = InputFile1.readlines()
                k1 = 1
                for line in lines:
                    sep_line = line.split()
                    if (((k1 <= 3)) or (len(sep_line) == 0)):
                        k1 = k1 + 1
                    else:
                        omegagrid1.append(sep_line[0])
                        f1.append(sep_line[1])

                lines = InputFile2.readlines()
                k1 = 1
                for line in lines:
                    sep_line = line.split()
                    if (((k1 <= 3)) or (len(sep_line) == 0)):
                        k1 = k1 + 1
                    else:
                        omegagrid2.append(sep_line[0])
                        f2.append(sep_line[1])
    
            omegagrid1 = np.asarray(omegagrid1, dtype='d')
            omegagrid2 = np.asarray(omegagrid2, dtype='d')
            f1 = np.asarray(f1, dtype='d')
            f2 = np.asarray(f2, dtype='d')
            dset_id = grp.create_dataset('Energy_f1', data=1e3*omegagrid1)
            dset_id.attrs['units'] = np.string_('[eV]')
            dset_id = grp.create_dataset('Energy_f2', data=1e3*omegagrid2)
            dset_id.attrs['units'] = np.string_('[eV]')
            dset_id = grp.create_dataset('f1', data=f1)
            dset_id.attrs['units'] = np.string_('[-]')
            dset_id = grp.create_dataset('f2', data=f2)
            dset_id.attrs['units'] = np.string_('[-]')
# print(omegagrid)


print('done')
