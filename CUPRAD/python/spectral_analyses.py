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
import matplotlib.pyplot as plt
import re
# import copy
import glob
import pandas as pd
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index
import subprocess
import dataformat_CUPRAD as dfC

import gc

try:
    univ_input_path = os.environ['UNIV_INPUT_PATH']
except:
    univ_input_path = os.path.join("D:\git", "universal_input")
        

print('path:', univ_input_path)
print('test1')
# sys.exit(0)
# There may be a large memory overload in this implementation, it's for testing putposes

# =============================================================================
# Inputs of the script

arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
else:
    # results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
    # results_path = os.path.join("D:\data", "Discharges")
    # results_path = os.path.join("D:\data", "Discharges", "f_scan_Kerr")
    # results_path = os.path.join("D:\data", "Discharges","TDSE", "t6")
    # results_path = os.path.join("D:\data", "Discharges", "preion2")
    results_path = os.path.join("D:\data", "Discharges","CUPRAD_tests")

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'

# results_path = os.path.join("D:\TEMP", "OCCIGEN_CUPRAD", "tests", "sim")

cwd = os.getcwd()
os.chdir(results_path)
files = glob.glob('results_*.h5')
os.chdir(cwd)

vacuum_frame = True

# files = ['results_1.h5','results_25.h5','results_2.h5']
# files = ['results_1.h5','results_2.h5', 'results_3.h5']
# files = ['results.h5']
# files = ['results_19.h5']

# files = ['results_10.h5']

# files = ['results_1.h5','results_10.h5', 'results_15.h5']

# files = ['results_10.h5']

out_h5name = 'analyses.h5'

# create inputs from the free-form inp
# run_args = ['python3', univ_input_path +'/create_universal_HDF5.py',
#             '-i', 'map3.inp', '-ohdf5', 'inputs_tmp.h5',
#             '-g', 'inputs']

# subprocess.run(run_args) # problematic in Spyder
# # runfile(run_args)

print('archive created')
print('my dir is:',os.getcwd())




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
    Horders = [15, 17, 19, 21, 23] # [19, 21, 23, 25, 27]
    Lcoh_saturation = 0.02
    Lcoh_zero = 0.0
    H_shift_mask = 4.0
    tlim = [-60.0,60.0]
    t_probe = 1e-15*np.asanyarray([-5.0, 0.0, 5.0])
    rmax = 130e-6 # only for analyses
    dr = rmax/40.0

fix_saturation = True

full_resolution = False

invoke_garbage_collector = True

OutPath = 'outputs_spectral'

# t_fix = 0.0e-15 # the time of our interest to inspect e.g. phase
# q = 23 # harmonic of our interest
# Coherence_length = True
# Beam_analysis = False
# Efield_analysis = False
# Gaussian_curvature = True # print Gaussian curvature, it is applied only in the first run
# fluence_source = 'computed' # options: 'file', 'computed'
# sys.exit(0)

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

prop_cycle = plt.rcParams['axes.prop_cycle']
colors_plt = prop_cycle.by_key()['color']
linestyles_plt = 10*['--','-','-.'] # horrible hot-fix


NH = len(Horders)
Nt_probe = len(t_probe)

# df_delta_t_columns = pd.MultiIndex.from_arrays([['pressure', 'ionisation', 'Intensity'        ,'IR phase','IR group','XUV phase'],
#                                                 ['{[bar]}',     '[\%]',      '$10^{14}$ [W/cm2]','[fs]',    '[fs]',    '[fs]']])


#  ['pressure \\\ [mbar]', 'ionisation \\\ [\%]', 'Intensity \\\ $10^{14}$ [W/cm2]', 'IR phase \\\ [fs]', 'IR group \\\ [fs]', 'XUV phase \\\ [fs]']
# df_delta_t = pd.DataFrame(columns = df_delta_t_columns)
# df_delta_t.loc[len(df_delta_t)+1] = [10, 15, 20]

with h5py.File(out_h5name,'w') as OutFile: # this file contains numerical analyses
# This is the archive, where we store numerical results at the instant. It may
# be further extended by storing also figure sources etc.
    for fname in files:
    # Here we loop over all result files in the destiantion folder.
        filename = fname
        
        # Using regexp to obtain the rusult number as a literal.
        k_sim = re.findall(r'\d+', filename)
        k_sim = k_sim[0] 
        
        # Prepare group in OutFile and open the repsective result file
        # grp = OutFile.create_group('results_' + str(k_sim))   
        file_path = os.path.join(results_path,filename)
        print('processing:', file_path)             
        with h5py.File(file_path, 'r') as InputArchive:
            
            # load data
            res = dfC.get_data(InputArchive, r_resolution = [full_resolution, dr, rmax])
            # print(res.Gaussian_zR)
            
            ### Shift to the vacuum frame
            if vacuum_frame:
                res.vacuum_shift()
         

            # ===============================================
            # Get the spectra
            
            title_string = dfC.create_param_string(print_parameters, res)
            
            res.compute_spectrum(compute_dE_domega = True)

            # fig, ax = plt.subplots()     
            # plt.plot(res.ogrid/res.omega0,np.abs(res.FE_trz[:,0,0]))
            # plt.show()
            
            fig, ax = plt.subplots()  
            ax.semilogy(res.ogrid/res.omega0,np.abs(res.FE_trz[:,0,0])**2, label='z='+"{:.1f}".format(1e3*res.zgrid[0]) )
            ax.semilogy(res.ogrid/res.omega0,np.abs(res.FE_trz[:,0,(res.Nz-1)//2])**2, label='z='+"{:.1f}".format(1e3*res.zgrid[(res.Nz-1)//2]) )
            ax.semilogy(res.ogrid/res.omega0,np.abs(res.FE_trz[:,0,-1])**2, label='z='+"{:.1f}".format(1e3*res.zgrid[-1]))
            ax.legend(loc='best')
            ax.set_xlabel('omega [-]'); ax.set_ylabel('|E(omega)|^2 [arb. u.]'); ax.set_title('FEfield'+title_string)  
            fig.savefig('FEfield_onax_sim'+str(k_sim)+'.png', dpi = 600)
            plt.show()  
 
            
  
            # ax10.plot(1e6*res.rgrid,Cutoff[k_t,:,(res.Nz-1)//2], color=colors_plt[1],) 
            # ax10.plot(1e6*res.rgrid,Cutoff[k_t,:,-1], color=colors_plt[2],)
            
            
            fig, ax = plt.subplots()     
            ax.plot(res.energy_zgrid,1e3*res.energy)
            ax.set_xlabel('z [-]'); ax.set_ylabel('E [mJ]'); ax.set_title('FEfield'+title_string) 
            fig.savefig('Energy_z_sim'+str(k_sim)+'.png', dpi = 600)
            plt.show()   
            
            
            spectrum_logscale = np.log10(abs(res.FE_trz[:,0,:])**2);
            vmin = np.max(spectrum_logscale) - 8.0# FF_orders_plot
            
            fig7, ax7 = plt.subplots()
            map1 = ax7.pcolor(1e3*res.zgrid, res.ogrid/res.omega0, spectrum_logscale, shading='auto', vmin=vmin, cmap='plasma')
            # map1 = ax7.pcolor(1e3*res.zgrid, res.ogrid/res.omega0, np.abs(res.FE_trz[:,0,:]), shading='auto')
            # map2 = ax7.contour(1e3*res.zgrid, 1e6*res.rgrid, Cutoff[t_probe_ind[k1],:,:], Horders, colors = "black")
            ax7.set_xlabel('z [mm]'); ax7.set_ylabel('omega [-]'); ax7.set_title('onax spectrum'+title_string)
            ax7.set_ylim([0.5,1.5])
            fig7.colorbar(map1) 
            
            fig7.savefig('FEfield_onax_map_sim'+str(k_sim)+'.png', dpi = 600)
            
            # dE_domega = res.dE_domega
            
            fig, ax = plt.subplots()  
            plt.semilogy(res.ogrid/res.omega0,res.dE_domega[:,0], label='z='+"{:.1f}".format(1e3*res.zgrid[0]) )
            plt.semilogy(res.ogrid/res.omega0,res.dE_domega[:,(res.Nz-1)//2], label='z='+"{:.1f}".format(1e3*res.zgrid[(res.Nz-1)//2]) )
            plt.semilogy(res.ogrid/res.omega0,res.dE_domega[:,-1], label='z='+"{:.1f}".format(1e3*res.zgrid[-1]))
            ax.legend(loc='best')
            ax.set_xlabel('omega [-]'); ax.set_ylabel('|E(omega)|^2 [arb. u.]'); ax.set_title('dE/domega'+title_string)  
            fig.savefig('dE_domega_sim'+str(k_sim)+'.png', dpi = 600)
            plt.show()
            
            
            spectrum_logscale = np.log10(abs(res.dE_domega));
            vmin = np.max(spectrum_logscale) - 8.0# FF_orders_plot
            
            fig7, ax7 = plt.subplots()
            map1 = ax7.pcolor(1e3*res.zgrid, res.ogrid/res.omega0, spectrum_logscale, shading='auto', vmin=vmin, cmap='plasma')
            # map1 = ax7.pcolor(1e3*res.zgrid, res.ogrid/res.omega0, np.abs(res.FE_trz[:,0,:]), shading='auto')
            # map2 = ax7.contour(1e3*res.zgrid, 1e6*res.rgrid, Cutoff[t_probe_ind[k1],:,:], Horders, colors = "black")
            ax7.set_xlabel('z [mm]'); ax7.set_ylabel('omega [-]'); ax7.set_title('dE/domega'+title_string)
            ax7.set_ylim([0.5,1.5])
            fig7.colorbar(map1) 
            
            fig7.savefig('dE_domega_map_sim'+str(k_sim)+'.png', dpi = 600)
            
            if showplots: plt.show()
            plt.close()          
            # ogrid, dum, Nt = mn.fft_t(res.tgrid, res.E_trz[:,0,0])
            # fig, ax = plt.subplots()     
            # plt.plot(ogrid/res.omega0,np.abs(dum))
            # plt.show()
            
            
            # FE_trz = np.zeros((len(ogrid),len(res.rgrid),len(res.zgrid)),dtype=complex)
            
            # for k1 in range(res.Nz):
            #     for k2 in range(res.Nr):
            #         FE_trz[:,k2,k1] = mn.fft_t(res.tgrid, res.E_trz[:,k2,k1])[1]

            # fig, ax = plt.subplots()     
            # plt.plot(ogrid/res.omega0,np.abs(FE_trz[:,0,-1]))
            # plt.show()            
            
            
            # sys.exit()

            
        # if invoke_garbage_collector:
        #     del res
        #     gc.collect()
        #     plt.close('all')
        

os.chdir(cwd)
print('done')