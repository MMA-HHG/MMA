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
import plot_presets as pp
import string

import gc

try:
    univ_input_path = os.environ['UNIV_INPUT_PATH']
except:
    univ_input_path = os.path.join("D:\git", "universal_input")
        

alphabet_list = list(string.ascii_lowercase)

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
    results_path = os.path.join("D:\data", "Discharges","scan_pI0_fine")

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

tlim2 = [-20,5]

map_scale = 'nipy_spectral'

fix_saturation = True

full_resolution = False

invoke_garbage_collector = True

OutPath = 'Cutoff_maps_PhD'

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
            
            
            
            ### Shift to the vacuum frame
            if vacuum_frame:
                res.vacuum_shift()
         

            # ===============================================
            # Complexify the fields: E(r,z,t) = Re(E_cmplx(r,z,t))
            E_trz_cmplx_envel = res.complexify_envel()


            # ===============================================
            # Retrieve the phase at the time of our interest
            t_probe_ind = mn.FindInterval(res.tgrid, t_probe) # find time indices
            
            phase = np.angle(E_trz_cmplx_envel[t_probe_ind,:,:]) # ordering (r,z)
            Intens = mn.FieldToIntensitySI(abs(E_trz_cmplx_envel))
            

            

                    
            nXUV = [];
            for k1 in range(NH):
                q = Horders[k1]
                f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(q*res.omega0, 'omegaSI', 'eV'))
                nXUV.append(1.0 - res.rho0_init*units.r_electron_classical*(mn.ConvertPhoton(q*res.omega0,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi))
            
                
            # Cut-off maps
            Cutoff = HHG.ComputeCutoff(Intens/units.INTENSITYau,
                                       mn.ConvertPhoton(res.omega0,'omegaSI','omegaau'),
                                       mn.ConvertPhoton(res.Ip_eV,'eV','omegaau')
                                       )[1]
            
            Cutoff_max = np.empty(Cutoff.shape[1:3])
            for k1 in range(res.Nr):
                for k2 in range(res.Nz):
                    Cutoff_max[k1,k2] = np.max(Cutoff[:,k1,k2])
                    
            
            
            image = pp.figure_driver()
            image.set_fontsizes = 'doublet+'     
            image.sf = [pp.plotter() for k1 in range(2)]
            
            image.sf[0].method = plt.pcolor    
            image.sf[0].args = [1e3*res.zgrid,1e6*res.rgrid,Cutoff_max]    
            image.sf[0].kwargs = {'shading' : 'auto', 'cmap' : 'plasma'}             
            image.sf[0].colorbar.show = True
            image.sf[0].colorbar.kwargs = {'label': r'Cutoff [-]'}
            
            image.sf[1].method = plt.contour
            image.sf[1].args = [1e3*res.zgrid, 1e6*res.rgrid, Cutoff_max, Horders]    
            image.sf[1].kwargs = {'colors' : 'black'} 
            image.sf[1].colorbar.show_contours = True
             
              
            image.xlabel = r'$z$ [mm]'; image.ylabel = r'$\rho$ [$\mu$m]'
            
            image.title = 'Maximum'
            
            # outfpath = os.path.join(OutPath,'cutoffmax') #+'{:.0f}'.format(100*preions[choice1[0]]))
            outfpath = 'cutoffmax' 
            image.savefigs_args = [[outfpath + '.pdf'], [outfpath + '.png']]
            image.savefigs_kwargs = [{'bbox_inches' : 'tight'} for k1 in range(2)]
            
            image.annotation = [['(d)'],
                {'xy' : (0.05, .9),
                  'xycoords' : 'axes fraction',
                  'color' : 'w'}]  
            # image.sf[0].colorbar.kwargs['location'] = 'top'
            
            pp.plot_preset(image)
            
            for k1 in range(Nt_probe):
                # t_string = ', '+"{:.1f}".format(1e15*res.tgrid[t_probe_ind[k1]])+' fs'
                t_string = "{:.1f}".format(1e15*res.tgrid[t_probe_ind[k1]])+'$ fs'
                
                image = pp.figure_driver()
                image.set_fontsizes = 'doublet+'     
                image.sf = [pp.plotter() for k1 in range(2)]
                
                image.sf[0].method = plt.pcolor    
                image.sf[0].args = [1e3*res.zgrid,1e6*res.rgrid,Cutoff[t_probe_ind[k1],:,:]]    
                image.sf[0].kwargs = {'shading' : 'auto', 'cmap' : 'plasma'}                
                image.sf[0].colorbar.show = True
                image.sf[0].colorbar.kwargs = {'label': r'Cutoff [-]'}
                
                image.sf[1].method = plt.contour
                image.sf[1].args = [1e3*res.zgrid, 1e6*res.rgrid, Cutoff[t_probe_ind[k1],:,:], Horders]    
                image.sf[1].kwargs = {'colors' : 'black'} 
                image.sf[1].colorbar.show_contours = True
                  
                image.xlabel = r'$z$ [mm]'; image.ylabel = r'$\rho$ [$\mu$m]'
                
                image.title = r'$t='+t_string
                
                # outfpath = os.path.join(OutPath,'cutoffmax') #+'{:.0f}'.format(100*preions[choice1[0]]))
                outfpath = 'cutofftfix' + str(k1)
                image.savefigs_args = [[outfpath + '.pdf'], [outfpath + '.png']]
                image.savefigs_kwargs = [{'bbox_inches' : 'tight'} for k1 in range(2)]
                
                image.annotation = [['('+alphabet_list[k1]+')'],
                {'xy' : (0.05, .9),
                  'xycoords' : 'axes fraction',
                  'color' : 'w'}]
                
                pp.plot_preset(image)
            
            
            

            
            
            
        if invoke_garbage_collector:
            del res
            # del res2
            gc.collect()
            # plt.close('all')
        

os.chdir(cwd)
print('done')