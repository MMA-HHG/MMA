import h5py
import os
import glob
import subprocess
import argparse

import MMA_administration as MMA

### Argument parser
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--paramfile", required=True, help="Input parameter file.")
ap.add_argument("-o", "--outhdf5", required=True, help="Output HDF5 file. ")
args = vars(ap.parse_args())

param_file = args['paramfile']
print("Parameter file: ", param_file)
results_file = glob.glob(args['outhdf5'])
print("Results file: ", results_file[0])

try:
    with h5py.File(results_file[0], 'r') as InputArchive:
        Nz = len(InputArchive[MMA.paths['CUPRAD_output']+'/zgrid'][:])
except FileNotFoundError:
    print("File '{}' not found, check the repository!".format(results_file[0]))


content = 'Nz_max\t'+str(Nz)+'\tI\t-\n' # Nz_max  1250 I   -

with open(param_file,'r') as inp_part, open('TDSE.inp.tmp','w') as inp_tmp:
    inp_tmp.write(inp_part.read()+'\n'+content)


run_args = ['python3', os.path.join(os.environ['UNIV_INPUT_PATH'],'create_universal_HDF5.py'),
            '-i', 'TDSE.inp.tmp',
            '-ihdf5', results_file[0],
            '-ohdf5', results_file[0],
            '-g',     MMA.paths['CTDSE_inputs']]

subprocess.run(run_args)

print("HDF5 file ready for MPI-TDSE.")
