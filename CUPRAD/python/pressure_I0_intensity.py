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
    # results_path = os.path.join("D:\data", "Discharges","CUPRAD_tests")
    results_path = os.path.join("D:\data", "Discharges","scan_pressures1")

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'

# results_path = os.path.join("D:\TEMP", "OCCIGEN_CUPRAD", "tests", "sim")

cwd = os.getcwd()
os.chdir(results_path)
files = glob.glob('results_*.h5')
driving_h5file = 'results.h5'
os.chdir(cwd)

# ensure 'results_1.h5' being the first, to be uniquely defined as glob doesn't impose ordering
res1_index = files.index('results_1.h5')
files[0], files[res1_index] = files[res1_index], files[0]

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

fix_saturation = True

full_resolution = False

invoke_garbage_collector = True

OutPath = 'outputs_pressureI0_scan'

# t_fix = 0.0e-15 # the time of our interest to inspect e.g. phase
# q = 23 # harmonic of our interest
# Coherence_length = True
# Beam_analysis = False
# Efield_analysis = False
# Gaussian_curvature = True # print Gaussian curvature, it is applied only in the first run
# fluence_source = 'computed' # options: 'file', 'computed'
# sys.exit(0)

# # Get FSPA 
# with h5py.File(os.path.join(cwd,'FSPA_tables_Krypton_test.h5'),'r') as h5_FSPA_tables:
#     interp_FSPA_short = HHG.FSPA.get_dphase(h5_FSPA_tables,'Igrid','Hgrid','short/dphi')

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

# load driving params
file_path = os.path.join(results_path,driving_h5file)
with h5py.File(file_path, 'r') as InputArchive:
    pressure_list_mbar = np.unique(1e3*InputArchive['/inputs/medium_pressure_in_bar'][:])
    I0_list = np.unique(InputArchive['/inputs/laser_focus_intensity_Gaussian'][:])
    


if not('-here' in arguments):
    pressure_list_mbar = pressure_list_mbar[0:3]  

N_press = len(pressure_list_mbar) 
N_I0 = len(I0_list) 

pressure_5mbar_steps = mn.get_divisible_interior_points(
                       [pressure_list_mbar[0],pressure_list_mbar[-1]],
                       5) 
    
# sys.exit()
    
    

Firstrun = True    
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
            res.get_plasma(InputArchive, r_resolution = [full_resolution, dr, rmax])
            
            if vacuum_frame: res.vacuum_shift() # Shift to the vacuum frame    
            E_trz_cmplx_envel = res.complexify_envel() # Complexify the fields: E(r,z,t) = Re(E_cmplx(r,z,t))
            Intens = mn.FieldToIntensitySI(abs(E_trz_cmplx_envel))
            
            k_press = np.where(pressure_list_mbar == res.pressure_mbar)[0][0] #pressure_list.index(res.pressure_mbar)
            k_I0 = np.where(I0_list == res.Intensity_Gaussian_focus)[0][0]
            
            
            
            print(k_press)
            
            if Firstrun: # create outfiles etc.
                Intens_tmax_pI0map = np.empty((N_press,N_I0,res.Nr,res.Nz))
                zgrid_ref = 1.*res.zgrid; Nz_ref = len(zgrid_ref)
                rgrid_ref = 1.*res.rgrid; Nr_ref = len(rgrid_ref)
                k_w0 = mn.FindInterval(rgrid_ref, res.w0_entry)
                plasma_end_pI0map = np.empty((N_press,N_I0,res.Nr,res.Nz))
                plasma_tmax_pI0map = np.empty((N_press,N_I0,res.Nr,res.Nz))
                z_half = 0.5*zgrid_ref[-1]
                kz_half = mn.FindInterval(zgrid_ref, z_half)
                Firstrun = False
            
            
            Intens_tmax = np.empty(Intens.shape[1:3])  
            plasma_end = np.empty(res.plasma.value_trz.shape[1:3]) 
            plasma_tmax = np.empty(res.plasma.value_trz.shape[1:3])  
            for k1 in range(res.Nz):
              for k2 in range(res.Nr):
                index_of_max = np.argmax(Intens[:,k2,k1])
                Intens_tmax[k2,k1] = Intens[index_of_max,k2,k1]
                plasma_end[k2,k1] = 100*res.plasma.value_trz[-1,k2,k1]/res.rho0_init
                plasma_tmax[k2,k1] = 100*res.plasma.value_trz[index_of_max,k2,k1]/res.rho0_init
                
                
            # interp on the reference grid            
            for k1 in range(res.Nr):
                Intens_tmax_pI0map[k_press,k_I0,k1,:] = np.interp(zgrid_ref,res.zgrid,Intens_tmax[k1,:])
                plasma_end_pI0map[k_press,k_I0,k1,:] = np.interp(zgrid_ref,res.zgrid, plasma_end[k1,:])
                plasma_tmax_pI0map[k_press,k_I0,k1,:] = np.interp(zgrid_ref,res.zgrid,plasma_tmax[k1,:])
                
                
    
    
    
    Cutoff_tmax_pmap = HHG.ComputeCutoff(
                         Intens_tmax_pI0map/units.INTENSITYau,
                         mn.ConvertPhoton(res.omega0,'omegaSI','omegaau'),
                         mn.ConvertPhoton(res.Ip_eV,'eV','omegaau')
                         )[1]




    ## Print intensities
    ######################
    
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,0,-1].T, shading='auto', cmap='plasma')
    map2 = ax1.contour(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,0,-1].T, Horders, colors = "black")
    
    # ax1.plot(1e3*res.zgrid, 1e6*radius_RMS, '--', linewidth=1, color = 'k')
    # ax1.plot(1e3*res.zgrid, 1e6*radius_inv_e2, ':', linewidth=1, color = 'k')
    
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Max cutoff, exit')
    fig1.colorbar(map1) 
    fig1.savefig('Cutoff_exit.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()                    

    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,0,0].T, shading='auto', cmap='plasma')  
    map2 = ax1.contour(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,0,0].T, Horders, colors = "black")
              
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Max cutoff, entry')
    fig1.colorbar(map1) 
    
    fig1.savefig('Cutoff_entry.png', dpi = 600)
    if showplots: plt.show()
    # plt.close() 
    

    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,0,kz_half+1].T, shading='auto', cmap='plasma') 
    map2 = ax1.contour(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,0,kz_half+1].T, Horders, colors = "black")
              
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Max cutoff, middle')
    fig1.colorbar(map1) 
    
    fig1.savefig('Cutoff_middle.png', dpi = 600)
    if showplots: plt.show()
    # plt.close() 
    

    # fig1, ax1 = plt.subplots()
    # map1 = ax1.pcolor(pressure_list_mbar, 1e3*zgrid_ref, Cutoff_tmax_pmap[:,0,:].T, shading='auto', cmap='plasma')  
    # map2 = ax1.contour(pressure_list_mbar, 1e3*zgrid_ref, Cutoff_tmax_pmap[:,0,:].T, Horders, colors = "black")
              
    # # ax1.set_ylim([0,1e6*rmax])
    # ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('z [mm]'); ax1.set_title('Max cutoff, on-axis')
    # fig1.colorbar(map1) 
    
    # fig1.savefig('Cutoff_onax.png', dpi = 600)
    # if showplots: plt.show()
    # # plt.close()   
    
    
    
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,k_w0,-1].T, shading='auto', cmap='plasma')  
    map2 = ax1.contour(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,k_w0,-1].T, Horders, colors = "black")
              
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Max cutoff, w0')
    fig1.colorbar(map1) 
    
    fig1.savefig('Cutoff_w0.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()  
 
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,k_w0//2,-1].T, shading='auto', cmap='plasma')  
    map2 = ax1.contour(pressure_list_mbar, I0_list, Cutoff_tmax_pmap[:,:,k_w0//2,-1].T, Horders, colors = "black")
              
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Max cutoff, w0/2')
    fig1.colorbar(map1) 
    
    fig1.savefig('Cutoff_w0_half.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()  
   
    
   
    
    ## Print plasmas - end
    ########################
 
    
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_end_pI0map[:,:,0,-1].T, shading='auto', cmap='plasma')
    
    # ax1.plot(1e3*res.zgrid, 1e6*radius_RMS, '--', linewidth=1, color = 'k')
    # ax1.plot(1e3*res.zgrid, 1e6*radius_inv_e2, ':', linewidth=1, color = 'k')
    
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma end, exit')
    fig1.colorbar(map1) 
    fig1.savefig('Plasma_end_exit.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()                    

    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_end_pI0map[:,:,0,0].T, shading='auto', cmap='plasma')
              

    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma end, entry')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_end_entry.png', dpi = 600)
    if showplots: plt.show()
    # plt.close() 
    

    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list ,plasma_end_pI0map[:,:,0,kz_half+1].T, shading='auto', cmap='plasma')
              

    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma end, middle')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_end_middle.png', dpi = 600)
    if showplots: plt.show()
    # plt.close() 
    

    # fig1, ax1 = plt.subplots()
    # map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_end_pI0map[:,0,:].T, shading='auto', cmap='plasma')
              
    # # ax1.set_ylim([0,1e6*rmax])
    # ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma end, on-axis')
    # fig1.colorbar(map1) 
    
    # fig1.savefig('Plasma_end_onax.png', dpi = 600)
    # if showplots: plt.show()
    # # plt.close()   
    
    
    
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_end_pI0map[:,:,k_w0,-1].T, shading='auto', cmap='plasma')
              
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma end, w0')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_end_w0.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()  
 
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_end_pI0map[:,:,k_w0//2,-1].T, shading='auto', cmap='plasma')
              
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma end, w0/2')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_end_w0_half.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()   
    
    
    # ## Print plasmas - tmax
    # ########################
    
    
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_tmax_pI0map[:,:,0,-1].T, shading='auto', cmap='plasma')
    
    # ax1.plot(1e3*res.zgrid, 1e6*radius_RMS, '--', linewidth=1, color = 'k')
    # ax1.plot(1e3*res.zgrid, 1e6*radius_inv_e2, ':', linewidth=1, color = 'k')
    
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma tmax, cutoff, exit')
    fig1.colorbar(map1) 
    fig1.savefig('Plasma_tmax_exit.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()                    

    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_tmax_pI0map[:,:,0,0].T, shading='auto', cmap='plasma')
              
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma tmax, entry')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_tmax_entry.png', dpi = 600)
    if showplots: plt.show()
    # plt.close() 
    

    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list ,plasma_tmax_pI0map[:,:,0,kz_half+1].T, shading='auto', cmap='plasma')
              
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma tmax, middle')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_tmax_middle.png', dpi = 600)
    if showplots: plt.show()
    # plt.close() 
    

    # fig1, ax1 = plt.subplots()
    # map1 = ax1.pcolor(pressure_list_mbar, 1e3*zgrid_ref, plasma_tmax_pmap[:,0,:].T, shading='auto', cmap='plasma')
              
    # # ax1.set_ylim([0,1e6*rmax])
    # ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('z [mm]'); ax1.set_title('Plasma tmax, on-axis')
    # fig1.colorbar(map1) 
    
    # fig1.savefig('Plasma_tmax_onax.png', dpi = 600)
    # if showplots: plt.show()
    # # plt.close()   
    
    
    
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_tmax_pI0map[:,:,k_w0,-1].T, shading='auto', cmap='plasma')
              
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma tmax, w0')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_tmax_w0.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()  
 
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(pressure_list_mbar, I0_list, plasma_tmax_pI0map[:,:,k_w0//2,-1].T, shading='auto', cmap='plasma')
              
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]'); ax1.set_title('Plasma tmax, w0/2')
    fig1.colorbar(map1) 
    
    fig1.savefig('Plasma_tmax_w0_half.png', dpi = 600)
    if showplots: plt.show()
    # plt.close()           
    
    
    
    
    # # plot cuts
    # #############    
    
    # k_presses = mn.FindInterval(pressure_list_mbar, pressure_5mbar_steps)
    # if not(hasattr(k_presses, "__len__")): k_presses =[k_presses]
    
    # # intens
    # fig1, ax1 = plt.subplots()
    # for k1 in k_presses:
    #     ax1.plot(1e3*zgrid_ref,Cutoff_tmax_pmap[k1,0,:]) 
    
    # # ax1.legend(loc='best')
    # ax1.set_xlabel('z [mm]'); ax1.set_ylabel('I [cutoff]'); ax1.set_title('Intensity, onax, tmax')   
    # fig1.savefig('Intens_cuts.png', dpi = 600)  
    
    # # plasma - end
    # fig1, ax1 = plt.subplots()
    # for k1 in k_presses:
    #     ax1.plot(1e3*zgrid_ref,plasma_end_pmap[k1,0,:]) 
    
    # # ax1.legend(loc='best')
    # ax1.set_xlabel('z [mm]'); ax1.set_ylabel('P_ion [%]'); ax1.set_title('Plasma end, onax')   
    # fig1.savefig('plasma_end_cuts.png', dpi = 600)  


    # # plasma - tmax
    # fig1, ax1 = plt.subplots()
    # for k1 in k_presses:
    #     ax1.plot(1e3*zgrid_ref,plasma_tmax_pmap[k1,0,:]) 
    
    # # ax1.legend(loc='best')
    # ax1.set_xlabel('z [mm]'); ax1.set_ylabel('P_ion [%]'); ax1.set_title('Plasma tmax, onax')   
    # fig1.savefig('plasma_tmax_cuts.png', dpi = 600)  


os.chdir(cwd)
print('done')