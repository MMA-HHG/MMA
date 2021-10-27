import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn


import warnings

import matplotlib
# matplotlib.rcParams['text.usetex'] = True
# import mynumerics as mn
import matplotlib.pyplot as plt

arguments = sys.argv

showplots = not('-nodisplay' in arguments)


# results_path = os.path.join("D:\data", "Discharges", "I0_p","preion_8")
results_path = os.path.join("D:\data", "Discharges", "I0_p","preion_8")

filename = 'analyses.h5'

FF_orders_plot = 4

filename_path = os.path.join(results_path,filename)
 
        
with h5py.File(filename_path, 'r') as InputArchive:
    # load data
   available_data = list(InputArchive.keys())
   
   Intens_map = InputArchive['Intensity_tmax_SI_p_I0_r_z'][:]
   Lcoh_map = InputArchive['Lcoh'][:]
   Lcoh_no_FSPA_map = InputArchive['Lcoh_no_FSPA'][:]
   zgrid = InputArchive['z_grid'][:]
   rgrid = InputArchive['r_grid'][:]
   I0_grid = InputArchive['I0_grid'][:]
   p_grid = InputArchive['p_grid'][:]
   Hgrid = InputArchive['Lcoh_Hgrid'][:]


contours = 1e3*np.asarray([0.0075, 0.015, 0.03, 0.06])
def plot_p_I0_map(data_map, title, fname, cmap='plasma',vmax=None):
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(p_grid, I0_grid, data_map.T, shading='auto', cmap=cmap, vmax=vmax)  
    
    map2 = ax1.contour(p_grid, I0_grid, data_map.T, contours, colors = "black", linestyles = ['dashdot','dotted','dashed','solid'])
    
    # ax1.set_xlim(left = 10)
    # ['solid', 'dashed', 'dashdot', 'dotted' ]
          
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]');
    ax1.set_title(title)
    cbar = fig1.colorbar(map1) 
    cbar.set_label(r'$L_{coh}$ [mm]')
    
    # fig1.savefig(fname, dpi = 600)
    if showplots: plt.show()
               
plot_p_I0_map(1e3*Lcoh_map[1,:,:,0,-1], 'H17', 'test.png',vmax=60, cmap = 'plasma')





############## plot intensity curves ############################

I0_indices = [4,13,19]
p_indices = [4,13,19]
colors = ["tab:orange","tab:blue","tab:green"]
linestyles = ['-','--',':']

pressures_round = np.round(p_grid[p_indices])
I0s_round = 1e18*np.round(1e-18*np.asarray(I0_grid[p_indices]),decimals=1)

pressures_leg = [str(pressure_round)+' mbar' for pressure_round in pressures_round]
I0s_leg = [str(I0_round)+' W/m2' for I0_round in I0s_round]

fig, ax = plt.subplots()    
for k1 in range(len(I0_indices)):
    for k2 in range(len(p_indices)):
        ax.plot(1e3*zgrid, Intens_map[p_indices[k2],I0_indices[k1],0,:],
                color=colors[k1],
                linestyle=linestyles[k2],
                linewidth=3)    
        
        

 
# ax.plot(1e3*zgrid, Intens_map[0,0,0,:], color="tab:orange", linewidth=3)
# ax.plot(1e3*zgrid, Intens_map[0,9,0,:], color="tab:blue", linewidth=3)
# ax.plot(1e3*zgrid, Intens_map[0,4,0,:], color="tab:green", linewidth=3)

# ax.plot(1e3*zgrid, Intens_map[4,0,0,:], color="tab:orange", linewidth=3, linestyle="--")
# ax.plot(1e3*zgrid, Intens_map[4,9,0,:], color="tab:blue", linewidth=3, linestyle="--")
# ax.plot(1e3*zgrid, Intens_map[4,4,0,:], color="tab:green", linewidth=3, linestyle="--")

# ax.plot(1e3*zgrid, Intens_map[9,0,0,:], color="tab:orange", linewidth=3, linestyle=":")
# ax.plot(1e3*zgrid, Intens_map[9,9,0,:], color="tab:blue", linewidth=3, linestyle=":")
# ax.plot(1e3*zgrid, Intens_map[9,4,0,:], color="tab:green", linewidth=3, linestyle=":")

# ax.set_ylabel("Intensity [W/m2]")
# ax.legend(loc=1, ncol=3)
from matplotlib.lines import Line2D
custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
                Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
                Line2D([0], [0], color="tab:blue", lw=3),
                Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
                Line2D([0], [0], color="tab:green", lw=3),                
                Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

ax.legend(custom_lines, [I0s_leg[0],
                         pressures_leg[0],
                         I0s_leg[1],
                         pressures_leg[1],
                         I0s_leg[2],
                         pressures_leg[2]],
          loc=1, ncol=3)

ax.set_title("On-axis defocusing")
ax.set_xlabel('z [mm]')
ax.tick_params(axis="both")
ax.set_ylabel("Intensity [W/m2]")

plt.show()


############## plot \Delta k curves ############################

I0_indices = [4,13,19]
p_indices = [4,13,19]
colors = ["tab:orange","tab:blue","tab:green"]
linestyles = ['-','--',':']

pressures_round = np.round(p_grid[p_indices])
I0s_round = 1e18*np.round(1e-18*np.asarray(I0_grid[p_indices]),decimals=1)

pressures_leg = [str(pressure_round)+' mbar' for pressure_round in pressures_round]
I0s_leg = [str(I0_round)+' W/m2' for I0_round in I0s_round]

fig, ax = plt.subplots()    
for k1 in range(len(I0_indices)):
    for k2 in range(len(p_indices)):
        ax.plot(1e3*zgrid, np.pi/Lcoh_map[1,p_indices[k2],I0_indices[k1],0,:],
                color=colors[k1],
                linestyle=linestyles[k2],
                linewidth=3)    



# ax.legend(loc=1, ncol=3)
from matplotlib.lines import Line2D
custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
                Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
                Line2D([0], [0], color="tab:blue", lw=3),
                Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
                Line2D([0], [0], color="tab:green", lw=3),                
                Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

ax.legend(custom_lines, [I0s_leg[0],
                         pressures_leg[0],
                         I0s_leg[1],
                         pressures_leg[1],
                         I0s_leg[2],
                         pressures_leg[2]],
          loc=1, ncol=3)

ax.set_title("H17")
ax.set_xlabel('z [mm]')
ax.tick_params(axis="both")
ax.set_ylabel("|$\Delta$ k| [1/m]")
ax.set_ylim([0,500])

plt.show()


############## plot \Delta k curves ############################

I0_indices = [4,13,19]
p_indices = [4,13,19]
colors = ["tab:orange","tab:blue","tab:green"]
linestyles = ['-','--',':']

pressures_round = np.round(p_grid[p_indices])
I0s_round = 1e18*np.round(1e-18*np.asarray(I0_grid[p_indices]),decimals=1)

pressures_leg = [str(pressure_round)+' mbar' for pressure_round in pressures_round]
I0s_leg = [str(I0_round)+' W/m2' for I0_round in I0s_round]

fig, ax = plt.subplots()    
for k1 in range(len(I0_indices)):
    for k2 in range(len(p_indices)):
        ax.plot(1e3*zgrid, np.pi/Lcoh_no_FSPA_map[1,p_indices[k2],I0_indices[k1],0,:],
                color=colors[k1],
                linestyle=linestyles[k2],
                linewidth=3)    



# ax.legend(loc=1, ncol=3)
from matplotlib.lines import Line2D
custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
                Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
                Line2D([0], [0], color="tab:blue", lw=3),
                Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
                Line2D([0], [0], color="tab:green", lw=3),                
                Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

ax.legend(custom_lines, [I0s_leg[0],
                         pressures_leg[0],
                         I0s_leg[1],
                         pressures_leg[1],
                         I0s_leg[2],
                         pressures_leg[2]],
          loc=1, ncol=3)

ax.set_title("H17, no FSPA")
ax.set_xlabel('z [mm]')
ax.tick_params(axis="both")
ax.set_ylabel("|$\Delta$ k| [1/m]")
ax.set_ylim([0,500])

plt.show()

# intensity map


contours = 1e3*np.asarray([15, 17, 19])
def plot_p_I0_map(data_map, title, fname, cmap='plasma',vmax=None):
    fig1, ax1 = plt.subplots()
    map1 = ax1.pcolor(p_grid, I0_grid, data_map.T, shading='auto', cmap=cmap, vmax=vmax)  
    
    map2 = ax1.contour(p_grid, I0_grid, data_map.T, contours, colors = "black")
    
    # ax1.set_xlim(left = 10)
    # ['solid', 'dashed', 'dashdot', 'dotted' ]
          
    # ax1.set_ylim([0,1e6*rmax])
    ax1.set_xlabel('p [mbar]'); ax1.set_ylabel('I0 [SI]');
    # ax1.set_title(title)
    cbar = fig1.colorbar(map1) 
    cbar.set_label(r'$I$ [cutoff]')
    
    # fig1.savefig(fname, dpi = 600)
    if showplots: plt.show()
               
plot_p_I0_map(1e3*Intens_map[:,:,0,-1], 'H17', 'test.png',vmax=None, cmap = 'plasma')

# fig, ax = plt.subplots()     
# ax.plot(zgrid, Intens_map[0,0,0,:], color="tab:orange", linewidth=3, label="5 mbar")
# ax.plot(zgrid, Intens_map[0,9,0,:], color="tab:blue", linewidth=3, label="5 mbar")
# ax.plot(zgrid, Intens_map[0,4,0,:], color="tab:green", linewidth=3, label="5 mbar")

# ax.plot(zgrid, Intens_map[4,0,0,:], color="tab:orange", linewidth=3, linestyle="--", label="25 mbar")
# ax.plot(zgrid, Intens_map[4,9,0,:], color="tab:blue", linewidth=3, linestyle="--", label="25 mbar")
# ax.plot(zgrid, Intens_map[4,4,0,:], color="tab:green", linewidth=3, linestyle="--", label="25 mbar")

# ax.plot(zgrid, Intens_map[9,0,0,:], color="tab:orange", linewidth=3, linestyle=":", label="50 mbar")
# ax.plot(zgrid, Intens_map[9,9,0,:], color="tab:blue", linewidth=3, linestyle=":", label="50 mbar")
# ax.plot(zgrid, Intens_map[9,4,0,:], color="tab:green", linewidth=3, linestyle=":", label="50 mbar")
    
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

# sys.exit()
    
# zgrid_integration_midpoints = 0.5*(zgrid_integration[1:]+zgrid_integration[:-1])
   
# Maxima = Maxima/np.max(Maxima)  
   
# fig, ax = plt.subplots()     
# plt.plot(Maxima[0,:])
# plt.plot(Maxima[1,:])
# plt.plot(Maxima[2,:])
# plt.plot(Maxima[3,:])
# plt.show()

# fig, ax = plt.subplots()     
# plt.semilogy(Maxima[0,:])
# plt.semilogy(Maxima[1,:])
# plt.semilogy(Maxima[2,:])
# plt.semilogy(Maxima[3,:])
# plt.show()

# fig, ax = plt.subplots()     
# plt.plot(Phases_onax[0,:])
# plt.plot(Phases_onax[1,:])
# plt.plot(Phases_onax[2,:])
# plt.plot(Phases_onax[3,:])
# plt.show()

# fig, ax = plt.subplots()     
# plt.plot(Phases_first[0,:])
# plt.plot(Phases_first[1,:])
# plt.plot(Phases_first[2,:])
# plt.plot(Phases_first[3,:])
# plt.show()

# # vmin = np.max(np.log(Gaborr))-6.
# fig, ax = plt.subplots()   
# FF_spectrum_logscale = np.log10(abs(FField_FF.T)**2);
# vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# # map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# map1 = ax.pcolor(Hgrid,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# fig.colorbar(map1)
# plt.title('Far-field spectrum (30 cm), integrated, log')
# plt.xlabel('H [-]')
# plt.ylabel('r [m]')
# plt.show()
# # if showplots: plt.show()
# # plt.close(fig)
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


# # Compute integrated spectrum
# dE_dH = np.zeros((NH,))
# for k1 in range(NH):
#     dE_dH[k1] = np.trapz(np.abs(FField_FF[k1,:])**2)


# fig, ax = plt.subplots()     
# plt.plot(Hgrid,dE_dH)
# plt.xlabel('H [-]')
# plt.ylabel('dE/dH [arb. u.]')
# plt.show()



# print(np.max(abs(FField_FF.T)**2))

# Hs_to_trace_maxima = []
# for k1 in range(len(Hgrid_study)):
#     Hs_to_trace_maxima.append([Hgrid_study[k1]-0.5 , Hgrid_study[k1]+0.5])

# H_indices = []
# planes_maxima = []
# XUV_beams_energy = []
# for H_list in Hs_to_trace_maxima:
#     try:
#         H_indices.append(mn.FindInterval(Hgrid, H_list))
#         planes_maxima.append([]); XUV_beams_energy.append([])
#     except:
#         warnings.warn("A frequency from frequencies_to_trace_maxima doesn't match ogrid.")

# dE_dH_z = np.zeros((NH,))
    
# if (len(H_indices)>0):
#     for k1 in range(Nz-1):
#         for k2 in range(NH): dE_dH_z[k2] = np.trapz(np.abs(FField_FF_cummulative[k1,k2,:])**2) 
#         for k2 in range(len(H_indices)):
#             planes_maxima[k2].append(np.max(abs(
#                               FField_FF_cummulative[k1,H_indices[k2][0]:H_indices[k2][1],:]
#                                     )))
#             XUV_beams_energy[k2].append(
#                               mn.integrate_subinterval(dE_dH_z,Hgrid,[Hgrid_study[k2]-0.5 , Hgrid_study[k2]+0.5])
#                                     )
            
#     for k1 in range(len(H_indices)):
#         planes_maxima[k1] = np.asarray(planes_maxima[k1])
#         XUV_beams_energy[k1] = np.asarray(XUV_beams_energy[k1])

# planes_maxima = np.asarray(planes_maxima)
# XUV_beams_energy = np.asarray(XUV_beams_energy)

# # find maxima in the case we need them
# fig, ax = plt.subplots()     
# plt.plot(planes_maxima[0,:])
# plt.plot(planes_maxima[1,:])
# plt.plot(planes_maxima[2,:])
# plt.plot(planes_maxima[3,:])
# plt.show()

# fig, ax = plt.subplots()     
# plt.plot(zgrid_integration_midpoints, planes_maxima[0,:])
# plt.plot(zgrid_integration_midpoints, planes_maxima[1,:])
# plt.plot(zgrid_integration_midpoints, planes_maxima[2,:])
# plt.plot(zgrid_integration_midpoints, planes_maxima[3,:])
# plt.show()

# fig, ax = plt.subplots()     
# plt.plot(zgrid_integration_midpoints, XUV_beams_energy[0,:])
# plt.plot(zgrid_integration_midpoints, XUV_beams_energy[1,:])
# plt.plot(zgrid_integration_midpoints, XUV_beams_energy[2,:])
# plt.plot(zgrid_integration_midpoints, XUV_beams_energy[3,:])
# plt.show()
   
   
   
# if L_abs_analyse:  
#     fig, ax = plt.subplots()     
#     plt.plot(Hgrid,1e3*L_abs_Hgrid)
#     plt.xlabel('H [-]')
#     plt.ylabel('Labs [mm]')
#     plt.show()
   
   
