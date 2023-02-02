import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn
import Hfn
import Hfn2

import re
import glob

import warnings


import matplotlib.pyplot as plt

  


arguments = sys.argv
if ('-here' in arguments):
    results_path = os.getcwd()
else:
    # results_path = os.path.join("D:\data", "Discharges", "TDSE","scan_mix2","exit")
    results_path = os.path.join("D:\data", "Discharges", "TDSE","convergence")

# cwd = os.getcwd()
# os.chdir(results_path)
# files = glob.glob('*.h5')
# os.chdir(cwd)

files = ['Hankel_all_cummulative.h5', 'Hankel_all_cummulative_dr2.h5']

N_files = len(files)

FF_orders_plot = 4

# p_grid = []; preion_grid = []
# for fname in files:
#     numbers = re.findall(r'\d+', fname)
#     p_grid.append(float(numbers[0])); preion_grid.append(float(numbers[1]))

# p_grid = np.unique(p_grid); preion_grid = np.unique(preion_grid); 
# N_press = len(p_grid); N_preion = len(preion_grid)
 
# sys.exit()

Firstrun = True
k1 = -1
for fname in files: # Here we loop over all result files in the destiantion folder.    
    numbers = re.findall(r'\d+', fname)
    k1 = k1+1
    filename = fname 
    
    
    # sys.exit()
    file_path = os.path.join(results_path,filename)
    
    with h5py.File(file_path, 'r') as InputArchive:
       numbers = re.findall(r'\d+', filename)
        # load data
       data_group = InputArchive['XUV']
       available_data = list(data_group.keys())
        
       # Maxima = InputArchive['XUV/Maxima_of_planes'][:] #/np.pi
       
       # Phases_onax = InputArchive['XUV/Phase_on_axis'][:] #/np.pi
       # Phases_first = InputArchive['XUV/Phase_first_plane'][:]
       
       FField_FF = (InputArchive['XUV/Spectrum_on_screen'][:,:,0] + \
                   1j*InputArchive['XUV/Spectrum_on_screen'][:,:,1])
       
       
       zgrid_integration = InputArchive['XUV/zgrid_integration'][:]; Nz = len(zgrid_integration)
       
                 
       L_abs_analyse = ('L_abs_Hgrid' in available_data)
       if L_abs_analyse:
           L_abs_Hgrid = InputArchive['XUV/L_abs_Hgrid'][:]
           
       if Firstrun:
           rgrid_FF = InputArchive['XUV/rgrid_FF'][:]
           Hgrid = InputArchive['XUV/Hgrid_select'][:]; NH = len(Hgrid)
           Hgrid_study = InputArchive['XUV/Maxima_Hgrid'][:]; NH_study = len(Hgrid_study)
           FField_FF_pp = np.empty((N_files,) + FField_FF.shape,dtype=np.cdouble)
           Firstrun = False

       FField_FF_pp[k1,:,:] = FField_FF

error_FF_abs = abs(FField_FF_pp[1,:,:]-FField_FF_pp[0,:,:])/np.max(abs(FField_FF_pp[0,:,:]))
    
OutPath = 'outputs_XUV_gain'

# sys.exit()               
 
# zgrid_integration_midpoints = 0.5*(zgrid_integration[1:]+zgrid_integration[:-1])


# Compute integrated spectrum
dE_dH = np.zeros((N_files, NH))
XUV_energy_pp = np.empty((N_files, NH_study))
for k1 in range(N_files):
    for k3 in range(NH):
      dE_dH[k1,k3] = np.trapz(np.abs(FField_FF_pp[k1,k3,:])**2)

# sys.exit()


# for k1 in range(N_press):
#   for k2 in range(N_preion):
    for k3 in range(NH_study):
      XUV_energy_pp[k1,k3] = mn.integrate_subinterval(
                                 dE_dH[k1,:],
                                 Hgrid,
                                 [Hgrid_study[k3]-0.5 , Hgrid_study[k3]+0.5]
                                 )


fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(error_FF_abs.T);
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
map1 = ax.pcolor(Hgrid,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field error log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
plt.show()


fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF_pp[0,:,:].T));
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
map1 = ax.pcolor(Hgrid,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field error log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
plt.show()

fig, ax = plt.subplots()   
FF_spectrum_logscale = np.log10(abs(FField_FF_pp[1,:,:].T));
vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
map1 = ax.pcolor(Hgrid,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
fig.colorbar(map1)
plt.title('Far-field error log')
plt.xlabel('H [-]')
plt.ylabel('r [m]')
plt.show()
# if showplots
    
sys.exit()

# store results
if os.path.exists(OutPath) and os.path.isdir(OutPath):
  shutil.rmtree(OutPath)
  print('deleted previous results')
os.mkdir(OutPath)

os.chdir(OutPath)



for k1 in range(NH_study):
    
    fig, ax = plt.subplots()     
    plt.plot(p_grid, XUV_energy_pp[:,0,k1])
    plt.plot(p_grid, XUV_energy_pp[:,1,k1])
    ax.set_xlabel('p [mbar]'); ax.set_ylabel('E [arb. u.]');
    ax.set_title('Energy , H'+str(Hgrid_study[k1]))
    fig.savefig('Energy_H'+str(Hgrid_study[k1])+'.png', dpi = 600)
    plt.show()
    
    
    fig, ax = plt.subplots()     
    plt.plot(p_grid, XUV_energy_pp[:,1,k1]/XUV_energy_pp[:,0,k1])
    ax.set_xlabel('p [mbar]'); ax.set_ylabel('E8/E0 [-]');
    ax.set_title('Amplification, H'+str(Hgrid_study[k1]))
    fig.savefig('Amplification_H'+str(Hgrid_study[k1])+'.png', dpi = 600)
    plt.show()

# fig, ax = plt.subplots()     
# plt.plot(p_grid, XUV_energy_pp[:,0,1])
# plt.plot(p_grid, XUV_energy_pp[:,1,1])
# plt.show()

# fig, ax = plt.subplots()     
# plt.plot(p_grid, XUV_energy_pp[:,0,2])
# plt.plot(p_grid, XUV_energy_pp[:,1,2])
# plt.show()

# kp = 4
for kp in range(N_press):
    fig, ax = plt.subplots()     
    plt.plot(Hgrid, dE_dH[kp,0,:])
    plt.plot(Hgrid, dE_dH[kp,1,:])
    ax.set_xlabel('H [-]'); ax.set_ylabel('dE/dH [arb. u.]');
    ax.set_title('dE/dH, p = '+str(p_grid[kp])+' mbar')
    fig.savefig('dE_dH_p_'+str(p_grid[kp])+'.png', dpi = 600)
    plt.show()


# sys.exit()

fig, ax = plt.subplots()     
plt.plot(p_grid, XUV_energy_pp[:,0,1]/np.max(abs(XUV_energy_pp[:,0,1])))
# plt.plot(Hgrid, dE_dH[kp,1,:])
# ax.set_xlabel('H [-]'); ax.set_ylabel('dE/dH [arb. u.]');
# ax.set_title('dE/dH, p = '+str(p_grid[kp])+' mbar')
# fig.savefig('dE_dH_p_'+str(p_grid[kp])+'.png', dpi = 600)
plt.show()


   
   


kp = 4
kpre = 0

kp_kpre = [[0,0],[0,1],[2,0],[2,1],[4,0],[4,1]]

for k1 in range(len(kp_kpre)):
    fig, ax = plt.subplots()   
    FF_spectrum_logscale = np.log10(abs(FField_FF_pp[kp_kpre[k1][0],kp_kpre[k1][1],:,:].T)**2);
    vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
    # map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
    map1 = ax.pcolor(Hgrid,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
    # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
    fig.colorbar(map1)
    plt.title('Far-field spectrum, log')
    plt.xlabel('H [-]')
    plt.ylabel('r [m]')
    plt.show()
    # if showplots: plt.show()
    # plt.close(fig)
    # # sys.exit()

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# # FF_spectrum_logscale = np.log10(abs(FField_FF.T)**2);
# # vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# # map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# map1 = ax.pcolor(Hgrid,rgrid_FF,abs(FField_FF.T)**2, shading='auto')
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), integrated')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # if showplots: plt.show()
# # plt.close(fig)
# # sys.exit()

# out_h5name = 'XUV_gains.h5'
# with h5py.File(out_h5name,'w') as OutFile: # this file contains numerical analyses
#     OutFile.create_dataset('XUV_energy_pp',
#                                            data = np.stack((FField_FF_pp.real, FField_FF_pp.imag),axis=-1)
#                            )
#     OutFile.create_dataset('Hgrid', data = Hgrid)
#     OutFile.create_dataset('Hgrid_study', data = Hgrid_study)
#     OutFile.create_dataset('rgrid_FF', data = rgrid_FF)
#     OutFile.create_dataset('pressure_grid', data = p_grid)
#     OutFile.create_dataset('preionisation_grid', data = preion_grid)
#     OutFile.create_dataset('XUV_energy', data = XUV_energy_pp)
#     OutFile.create_dataset('dE_dH', data = dE_dH)
    
os.chdir(cwd)