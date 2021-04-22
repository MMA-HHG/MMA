import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
sys.path.append('D:\git\python_modules')
import units
import mynumerics as mn
import HHG
# import mynumerics as mn
import matplotlib.pyplot as plt
import re
import glob
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index
from contextlib import ExitStack
import dataformat_CUPRAD as dfC

# There may be a large memory overload in this implementation, it's for testing putposes

# =============================================================================
# Inputs of the script
arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
else:
    # results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
    # results_path = os.path.join("D:\data", "Discharges", "f_scan_Kerr")
    results_path = os.path.join("D:\data", "Discharges", "preion2")


cwd = os.getcwd()
    
try:
    with h5py.File('inputs_tmp.h5', 'r') as Parameters:
        gas_type = mn.readscalardataset(Parameters, 'inputs/gas_type', 'S')
        XUV_table_type = mn.readscalardataset(Parameters, 'inputs/XUV_table_type', 'S') # 'NIST' # {Henke, NIST}
        Horders = Parameters['inputs/Horders'][:].tolist() # mn.readscalardataset(Parameters, 'inputs/gas_type', 'S') # [15, 17, 19, 21, 23] # [19, 21, 23, 25, 27]
        Lcoh_saturation = mn.readscalardataset(Parameters, 'inputs/Lcoh_saturation', 'N') #0.02
        Lcoh_zero = mn.readscalardataset(Parameters, 'inputs/Lcoh_zero', 'N') # 0.0
        H_shift_mask = mn.readscalardataset(Parameters, 'inputs/H_shift_mask', 'N') # 4.0
        tlim = Parameters['inputs/tlim'][:].tolist() # [-60.0,60.0]
        t_probe = Parameters['inputs/t_probe'][:]
        rmax = mn.readscalardataset(Parameters, 'inputs/rmax', 'N')
        dr = mn.readscalardataset(Parameters, 'inputs/dr', 'N')
        
        print_parameters = Parameters['inputs/print_parameters'][:].tolist()
        print_parameters = [k1.decode() for k1 in print_parameters]
        
except:
    print('error reading hdf5 file, using defaults') 
    
    gas_type = 'Kr'
    XUV_table_type = 'NIST' # {Henke, NIST}
    Horders = [15, 17, 19, 21, 23]# [19, 21, 23, 25, 27]
    Lcoh_saturation = 0.02
    Lcoh_zero = 0.0
    H_shift_mask = 4.0
    tlim = [-60.0,60.0]
    t_probe = 1e-15*np.asanyarray([-5.0, 0.0, 5.0])
    rmax = 130e-6 # only for analyses
    dr = rmax/40.0

vacuum_frame = True


files = ['results_1.h5', 'results_2.h5', 'results_3.h5', 'results_4.h5']
files = ['results_5.h5', 'results_6.h5', 'results_7.h5', 'results_8.h5']
# files = ['results_1.h5', 'results_2.h5']
         

full_resolution = False

OutPath = 'outputs'

# Get FSPA 
with h5py.File(os.path.join(cwd,'FSPA_tables_Krypton_test.h5'),'r') as h5_FSPA_tables:
    interp_FSPA_short = HHG.FSPA.get_dphase(h5_FSPA_tables,'Igrid','Hgrid','short/dphi')


# =============================================================================
# The body of the script
if os.path.exists(OutPath) and os.path.isdir(OutPath):
  shutil.rmtree(OutPath)
  print('deleted previous results')
os.mkdir(OutPath)

os.chdir(OutPath)



Nfiles = len(files)
res = []

with ExitStack() as stack:    
    InArch = [stack.enter_context(h5py.File(os.path.join(results_path,fname), 'r')) for fname in files]
    for k1 in range(Nfiles):
        res.append(dfC.get_data(InArch[k1], r_resolution = [full_resolution, dr, rmax])) # load data 
        ### Shift to the vacuum frame
        if vacuum_frame:
            res[k1].vacuum_shift()
            

# tlim = [-20.0,0.0]

fig1, ax1 = plt.subplots()
for k1 in range(4):
    ax1.plot(1e15*res[k1].tgrid, res[k1].E_trz[:,0,-1],label=res[k1].preionisation_string)

ax1.set_xlim(tlim)
ax1.legend(loc='upper right') # 'best'
ax1.set_xlabel('t [fs]'); ax1.set_ylabel('E [V/m]'); # ax1.set_title('onaxis intensity'+title_string)   
fig1.savefig('Efield.png', dpi = 600)
# ax1.plot(1e3*res.zgrid,q*(grad_z_phase[k1,0,:] + res.k0_wave*(nXUV[k2]-1)),label='H'+str(q))
         
    # color=colors_plt[k2], linestyle=linestyles_plt[k1])

if showplots: plt.show()
# plt.close()            

fig2, ax2 = plt.subplots()
ax2.plot(res[0].E_trz[:,0,0])
ax2.plot(res[0].E_trz[:,0,res[0].Nz//2])
ax2.plot(res[0].E_trz[:,0,-1])
# ax1.plot(res[1].E_trz[:,0,0])
# ax1.plot(res[2].E_trz[:,0,0])
# ax1.plot(1e3*res.zgrid,q*(grad_z_phase[k1,0,:] + res.k0_wave*(nXUV[k2]-1)),label='H'+str(q))
         
    # color=colors_plt[k2], linestyle=linestyles_plt[k1])

if showplots: plt.show()
# plt.close()        

os.chdir(cwd)
print('done')