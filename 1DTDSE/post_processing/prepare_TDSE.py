import h5py
import os
import glob
import subprocess
import argparse

### Argument parser
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--paramfile", required=True, help="Input parameter file.")
ap.add_argument("-o", "--outhdf5", required=False, help="Output HDF5 file. " \
                "Default path is '../CUPRAD/build/results.h5'")
args = vars(ap.parse_args())

#param_file = None
param_file = args['paramfile']
print("Parameter file", param_file)
results_file = None

#try:
#    h5_file = args['outhdf5']
#    print(h5_file)
#    results_file = glob.glob(h5_file)
#except KeyError:
if args['outhdf5'] == None:
    print("Looking up for an h5 file '../CUPRAD/build/results.h5'...")
    results_file = glob.glob('../CUPRAD/build/results.h5')
else:
    results_file = glob.glob(args['outhdf5'])

#try: 
#    param_file = sys.argv[1]
#except IndexError:
#    print("Append the name of the parameter file!")

#try:
#    results_file = glob.glob(sys.argv[2])
#except IndexError:
#    print("Looking up for an h5 file '../CUPRAD/build/results.h5'...")
#    results_file = glob.glob('../CUPRAD/build/results.h5')
#    pass


#results_file = glob.glob('../results_*.h5')
try:
    print("Found file", results_file[0])
except IndexError:
    print("HDF5 file not found!")

try:
    with h5py.File(results_file[0], 'r') as InputArchive:
        Nz = len(InputArchive['/outputs/zgrid'][:])
except FileNotFoundError:
    print("File '{}' not found, check the repository!".format(results_file[0]))


inp_part_path = os.path.join(os.environ['TDSE_1D_HOME'], param_file)
print("Reading parameter file: ", inp_part_path)

content = 'Nz_max\t'+str(Nz)+'\tI\t-\n' # Nz_max  1250 I   -

with open(inp_part_path,'r') as inp_part, open('TDSE.inp.tmp','w') as inp_tmp:
    inp_tmp.write(inp_part.read()+'\n'+content)


run_args = ['python3', os.path.join(os.environ['UNIV_INPUT_PATH'],'create_universal_HDF5.py'),
            '-i', 'TDSE.inp.tmp',
            '-ihdf5', results_file[0],
            '-ohdf5', 'results.h5',
            '-g', 'TDSE_inputs']

# print(run_args)
subprocess.run(run_args)

print('done')
