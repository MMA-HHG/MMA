import numpy as np
import h5py
import sys
import os
import shutil
import re
import glob
import subprocess

import matplotlib.pyplot as plt



# Load argument for input file - the name of the file is argument 0, so we want argument 1


## dealing with the input arguments
# output_file = ""
arguments = sys.argv
arg_index = arguments.index("-inpath")
inpath = arguments[arg_index+1]

files = glob.glob(inpath+'*')



cwd = os.getcwd()
for file in files:
    # get name and create directory
    print(file)
    print('t1')
    ind1 = file.rfind('mbar_') + len('mbar_')
    ind2 = file.rfind('.inp')
    specification = file[ind1:ind2]
    ind1 = file.rfind('TDSE_') + len('TDSE_')
    ind2 = file.rfind('mbar')
    pressure = file[ind1:ind2]
    os.mkdir(pressure+'_'+specification)
    
    os.chdir(pressure+'_'+specification)
    results_file = 'results_'+pressure+'_'+specification+'.h5'
    run_args = ['python3', os.environ['UNIV_INPUT_PATH']+'/create_universal_HDF5.py',
                '-i', file,
                '-ohdf5', results_file,
                '-g', 'inputs']
    
    subprocess.run(run_args)
    
    subprocess.run(os.environ['MULTISCALE_SCRIPTS']+'/run_multiscale.sh')
    
    # subprocess.run(program_path)
    # submit all jobs
    
    os.chdir(cwd)
    
    


# # arg_index = arguments.index("-ohdf5")
# # target_archive = arguments[arg_index+1]
    


# # inputfilename = 'HHGRun83.1.Kr.40A.35mbar.Z=1.dat'


# target_archive = inputfilename.replace('.dat','.h5')

# # ion_mult = 100.
# # r_point = 0.
# # separator = ','

# pressure_regexp = re.search('(\d+)(?=\s*mbar)', inputfilename)
# pressure_string = pressure_regexp.group(0)
# pressure_value = float(pressure_regexp.group(0))

# data_of_interest={
#     'tgrid': [],
#     'rgrid': [],
#     'preion': []
#     }

# First_Line = True
# with open(inputfilename, "r") as InputFile, h5py.File(target_archive, 'w') as GeneratedFile: # access option http://docs.h5py.org/en/stable/high/file.html#file
#     lines = InputFile.readlines()

#     # read data
#     for line in lines:
#         if (separator == 'default'):
#             sep_line = line.split()
#         else:
#             sep_line = line.split(separator) # re.split("\s|(?<!\d)[,.](?!\d)", line) #line.split()  # separate the line   
        
#         if First_Line: # first line is difficult to regex
#             First_Line = False
#             continue
            
        
#         data_of_interest['tgrid'].append(float(sep_line[0]))
#         data_of_interest['rgrid'].append(np.single(sep_line[1]))
#         data_of_interest['preion'].append(float(sep_line[-1])/ion_mult)
        
#     # reshape data
#     Ndata = len(data_of_interest['tgrid'])
#     tgrid = np.unique(data_of_interest['tgrid']); Nt = len(tgrid) # only t is uniform
#     data_of_interest['rgrid'] = np.asarray(data_of_interest['rgrid'])
#     data_of_interest['preion'] = np.asarray(data_of_interest['preion'])

   
#     preion = [] 
#     for k1 in range(Nt):
#         t_indices = np.where(data_of_interest['tgrid'] == tgrid[k1])[0]
#         rgrid_loc = data_of_interest['rgrid'][t_indices]
#         preion_loc = data_of_interest['preion'][t_indices]
#         preion.append(np.interp(r_point, rgrid_loc, preion_loc))
#     preion = np.asarray(preion)
 
#     GeneratedFile.create_dataset('pressure', data = pressure_value)
#     GeneratedFile.create_dataset('tgrid', data = tgrid)
#     GeneratedFile.create_dataset('preion', data = preion)
    
    
# # plt.plot(tgrid,preion)
# # plt.show()      


print('done')
