import numpy as np
import os
import time
import copy
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

import plot_presets as pp
  

import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index

from matplotlib.lines import Line2D
  


omegaSI = mn.ConvertPhoton(792e-9, 'lambdaSI', 'omegaSI') 
# Horder = 17


gas_type = 'Kr'
XUV_table_type_absorption = 'Henke' # {Henke, NIST}    
XUV_table_type_dispersion = 'NIST'
def f2_funct(E):
    return XUV_index.getf2(gas_type + '_' + XUV_table_type_absorption, E)

def sigma(omega):
    f2_value    = f2_funct(mn.ConvertPhoton(omega, 'omegaSI', 'eV'))
    lambdaSI    = mn.ConvertPhoton(omega, 'omegaSI', 'lambdaSI')
    return 2.0 * units.r_electron_classical * lambdaSI * f2_value


susc_IR = IR_index.getsusc(gas_type, mn.ConvertPhoton(omegaSI,'omegaSI','lambdaSI'))
# susc_IR = 0.0008395761480897157 # 800 nm from https://refractiveindex.info/?shelf=main&book=Kr&page=Borzsonyi

N_atm = 2.7e19 * 1e6
# N_atm = 2.6867774e25
def susc_XUV(omega_in_SI):
    f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type_dispersion, mn.ConvertPhoton(omega_in_SI, 'omegaSI', 'eV'))
    nXUV_atm = 1.0 - N_atm*units.r_electron_classical*(mn.ConvertPhoton(omega_in_SI,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
    return nXUV_atm**2 - 1
    # susc_XUV = nXUV_atm**2 - 1
    
def delta_susc(omega_in_SI):
    return (susc_IR - susc_XUV(omega_in_SI)) /N_atm
    
# delta_susc = (susc_IR - susc_XUV) /N_atm

def ionisation_ratio(A,I,q,omega):
    x1 = (q*omega/units.c_light)**2
    x2 = units.elcharge**2/((omega**2)*units.eps0*units.elmass)
    
    return ( (delta_susc(q*omega) - np.sqrt( ( (4.*A/I) - (sigma(q*omega))**2)/x1 ) )/x2 ,
             (delta_susc(q*omega) + np.sqrt( ( (4.*A/I) - (sigma(q*omega))**2)/x1 ) )/x2  )

def ionisation_ratio_optimal(q):
    return (omegaSI**2) *units.eps0*units.elmass*delta_susc(q*omegaSI)/(units.elcharge**2)
# ionisation_ratio_optimal = (omegaSI**2) *units.eps0*units.elmass*delta_susc/(units.elcharge**2)

def A_norm(I0,q,omega):
    return I0*((sigma(q*omega))**2 + (q*omega*delta_susc(q*omega)/units.c_light)**2)/4.      

def IntensXUV(eta,q,omega,A):
    return 4.*A/ ((sigma(q*omega))**2 + (q*omega*(delta_susc(q*omega) - eta * units.elcharge**2 /((omega**2)*units.eps0*units.elmass) )/units.c_light)**2)



arguments = sys.argv


# compare more results at once, we assume all parmaters to match (# of pressures, ...), so there are only one master p_grid, etc.
# results_paths = [os.path.join("D:\data", "Discharges", "TDSE","FF1","scan3"),
#                  os.path.join("D:\data", "Discharges", "TDSE","FF1","scan4")]
results_paths = [os.path.join("D:\data", "Discharges", "TDSE","scan3"),
                 os.path.join("D:\data", "Discharges", "TDSE","scan4")]

cwd = os.getcwd()

files = []
for results_path in results_paths:
    os.chdir(results_path)
    files.append(glob.glob('Hankel_all_cummulative*.h5'))
    os.chdir(cwd)

# filtered = filter(lambda file: not(contains_specifiers(file,['noabs','nodispersion'])), files)
# list(filtered)
for k1 in range(len(files)):
    files[k1] = list(filter(lambda file: not(mn.contains_substrings(file,['noabs','nodispersion'])), files[k1]))
# files = list(filter(lambda file: mn.contains_substrings(file,['nodispersion']), files))

# sys.exit()
filter_value = 0.1
FF_orders_plot = 4

filter_field_type = None #'low_abs2_cut'
filter_field_value = 0.05

filter_integral_type = None
filter_integral_value = 0.05


p_grid = [];
preion_extensions = ['no','half','end']
# preion_extensions = ['no.','half','end']
for fname in files[0]:
    numbers = re.findall(r'\d+', fname)
    p_grid.append(float(numbers[0]))

p_grid = np.unique(p_grid)
N_press = len(p_grid); N_preion = len(preion_extensions)
 
# sys.exit()

# load ionisations
ionisations = []
for results_path in results_paths:
    with h5py.File(os.path.join(results_path,'ionisations.h5'), 'r') as InputArchive:
        ionisations.append({
                    'half_init': InputArchive['ionisation_half_init'][:],
                    'end_init': InputArchive['ionisation_end_init'][:],
                    'half': InputArchive['ionisation_half_tmax'][:],  
                    'end': InputArchive['ionisation_end_tmax'][:]
                    })


FField_FF_pp = []
for k_1 in range(len(files)):
    Firstrun = True
    for fname in files[k_1]: # Here we loop over all result files in the destiantion folder.    
        numbers = re.findall(r'\d+', fname)
        p_value = float(numbers[0]); preion_value = float(numbers[1])
        filename = fname 
        
        k_press = np.where(p_grid == p_value)[0][0] #pressure_list.index(res.pressure_mbar)
        for k1 in range(N_preion):
            if ('_'+preion_extensions[k1]) in filename: k_preion=k1; break
        
        # sys.exit()
        # k_preion = np.where(preion_grid == preion_value)[0][0]
        
        # sys.exit()
        file_path = os.path.join(results_paths[k_1],filename)
        
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
               FField_FF_pp.append(np.empty((N_press, N_preion) + FField_FF.shape,dtype=np.cdouble))
               Firstrun = False
    
           FField_FF_pp[k_1][k_press,k_preion,:,:] = FField_FF
    
OutPath = 'outputs_XUV_gain'             
 
# zgrid_integration_midpoints = 0.5*(zgrid_integration[1:]+zgrid_integration[:-1])

def choice_to_label(choice_foo):
    if (choice_foo[2] == 0): eta0 = 0
    elif (choice_foo[2] == 1): eta0 = ionisations[choice_foo[0]]['half_init'][choice_foo[1]]
    elif (choice_foo[2] == 2): eta0 = ionisations[choice_foo[0]]['end_init'][choice_foo[1]]
    
    local_title = r'$p$='+ '{:.1f}'.format(p_grid[choice_foo[1]]) + ' mbar, ' + \
                  r'$\eta_0$='+ '{:.1f}'.format(eta0) + ' %'
    return local_title

# Compute integrated spectrum
FField_FF_pp_filtered = copy.deepcopy(FField_FF_pp)
kHs = [0,*mn.FindInterval(Hgrid,
                     mn.get_divisible_interior_points([Hgrid[0],Hgrid[-1]],2)),
       len(Hgrid)-1]

for k_1 in range(len(files)):
  for k1 in range(N_press):
    for k2 in range(N_preion):
      for k3 in range(len(kHs)-1):
        if not(kHs[k3] == kHs[k3+1]):
          FField_FF_pp_filtered[k_1][k1,k2,kHs[k3]:kHs[k3+1],:] = mn.clamp_array(FField_FF_pp_filtered[k_1][k1,k2,kHs[k3]:kHs[k3+1],:],
                                                                              filter_field_type,filter_field_value)
          
dE_dH = []; dE_dH_filtered = []
XUV_energy_pp = []; XUV_energy_pp_filtered = []
for k_1 in range(len(files)):
    dE_dH.append(np.zeros((N_press, N_preion, NH))); dE_dH_filtered.append(np.zeros((N_press, N_preion, NH)));
    XUV_energy_pp.append(np.empty((N_press, N_preion, NH_study))); XUV_energy_pp_filtered.append(np.empty((N_press, N_preion, NH_study)))
    
    for k1 in range(N_press):
      for k2 in range(N_preion):
        for k3 in range(NH):
          dE_dH[k_1][k1,k2,k3] = np.trapz(rgrid_FF*np.abs(FField_FF_pp[k_1][k1,k2,k3,:])**2)
          dE_dH_filtered[k_1][k1,k2,k3] = np.trapz(mn.clamp_array(
                                          rgrid_FF*np.abs(FField_FF_pp_filtered[k_1][k1,k2,k3,:])**2,
                                          filter_integral_type,filter_integral_value))
              

                                                                              
                                                   
                                                   


    # sys.exit()
    
    
    # for k1 in range(N_press):
    #   for k2 in range(N_preion):
        for k3 in range(NH_study):
          XUV_energy_pp[k_1][k1,k2,k3] = mn.integrate_subinterval(
                                     dE_dH[k_1][k1,k2,:],
                                     Hgrid,
                                     [Hgrid_study[k3]-0.5 , Hgrid_study[k3]+0.5]
                                     # [Hgrid_study[k3]-0.98 , Hgrid_study[k3]+0.98]
                                     )
          XUV_energy_pp_filtered[k_1][k1,k2,k3] = mn.integrate_subinterval(
                                     1.0*dE_dH_filtered[k_1][k1,k2,:],
                                     Hgrid,
                                     [Hgrid_study[k3]-0.5 , Hgrid_study[k3]+0.5]
                                     # [Hgrid_study[k3]-0.98 , Hgrid_study[k3]+0.98]
                                     )

      
# sys.exit()  

# store results
if os.path.exists(OutPath) and os.path.isdir(OutPath):
  shutil.rmtree(OutPath)
  print('deleted previous results')
os.mkdir(OutPath)

os.chdir(OutPath)


## plot H17 + ionisations
k1 = 1 # index to access given harmonic plot


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(11)]



image.sf[0].args = [p_grid, XUV_energy_pp[0][:,0,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
image.sf[1].args = [p_grid, XUV_energy_pp[0][:,1,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}

image.sf[2].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge/2'}
# image.sf[3].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.max(XUV_energy_pp[1][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}



# A0 = A_norm(np.mean(XUV_energy_pp[0][:,0,k1]),Hgrid_study[k1],omegaSI)

# image.sf[2].args = [p_grid, XUV_energy_pp[:,2,k1]/np.max(XUV_energy_pp[:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}



# image.sf[9].args = [p_grid, IntensXUV(1e-2*ionisations['half_init'],17,omegaSI,A0)/np.max(XUV_energy_pp[:,0,k1]),'g--'];
# image.sf[9].kwargs = {'label' : 'T_discharge/2 from analytical estimate'}    

for k2 in range(3,8): image.sf[k2].axis_location = 'right'
image.sf[3].args = [p_grid, ionisations[0]['half_init'], 'b:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[4].args = [p_grid, ionisations[0]['half'], 'b--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
image.sf[5].args = [p_grid, ionisations[1]['half_init'], 'r:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[6].args = [p_grid, ionisations[1]['half'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 

image.sf[7].args = [p_grid, 100*ionisation_ratio_optimal(Hgrid_study[k1])*np.ones(len(p_grid)), 'g--']
  
# # image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
# # image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  




# preions = ionisation_ratio(A0,XUV_energy_pp[0][:,1,k1],Hgrid_study[k1],omegaSI)
# image.sf[7].args = [p_grid, 100*preions[0][:], 'c--'];
# image.sf[8].args = [p_grid, 100*preions[1][:], 'c--'];
# image.sf[7].method = None
# image.sf[8].method = None

## custom legend
# custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
#                 Line2D([0], [0], color="tab:blue", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
#                 Line2D([0], [0], color="tab:green", lw=3),                
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

custom_lines = [Line2D([1], [0], color="k"),
                Line2D([0], [0], color="tab:grey", linestyle=":"),
                Line2D([0], [0], color="b"),
                Line2D([0], [0], color="tab:grey", linestyle="--"),
                Line2D([0], [0], color="r"),                
                Line2D([0], [0], color="g", linestyle="--")]

# ax.legend(custom_lines, [I0s_leg[0],
#                          pressures_leg[0],
#                          I0s_leg[1],
#                          pressures_leg[1],
#                          I0s_leg[2],
#                          pressures_leg[2]],
#           loc=1, ncol=3)

image.legend_args = [custom_lines,['no preion', r'$\eta_0$', '40 A', r'$\eta_{las.}$','50 A', '$\eta_{opt.}$']]
image.legend_kwargs = {'loc': 1, 'ncol': 3}

# image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
image.xlabel = 'p [mbar]'; image.ylabel = 'I_XUV [arb. u.]'; image.right_ylabel = 'ionisation [%]'

image.title = 'H'+str(Hgrid_study[k1]) + ', T_discharge/2'

image.savefig_args = ['compare1.pdf']
image.savefig_kwargs = {'dpi' : 600,'bbox_inches' : 'tight'}

pp.plot_preset(image)

# os.chdir(cwd)
# sys.exit()


###################################################################################################################
###################################################################################################################
###################################################################################################################


# ## plot H19
k1 = 2 # index to access given harmonic plot


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(11)]



image.sf[0].args = [p_grid, XUV_energy_pp[0][:,0,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
image.sf[1].args = [p_grid, XUV_energy_pp[0][:,1,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}

image.sf[2].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge/2'}
# image.sf[3].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.max(XUV_energy_pp[1][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}



# A0 = A_norm(np.mean(XUV_energy_pp[0][:,0,k1]),Hgrid_study[k1],omegaSI)

# image.sf[2].args = [p_grid, XUV_energy_pp[:,2,k1]/np.max(XUV_energy_pp[:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}



# image.sf[9].args = [p_grid, IntensXUV(1e-2*ionisations['half_init'],17,omegaSI,A0)/np.max(XUV_energy_pp[:,0,k1]),'g--'];
# image.sf[9].kwargs = {'label' : 'T_discharge/2 from analytical estimate'}    

for k2 in range(3,8): image.sf[k2].axis_location = 'right'
image.sf[3].args = [p_grid, ionisations[0]['half_init'], 'b:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[4].args = [p_grid, ionisations[0]['half'], 'b--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
image.sf[5].args = [p_grid, ionisations[1]['half_init'], 'r:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[6].args = [p_grid, ionisations[1]['half'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 

image.sf[7].args = [p_grid, 100*ionisation_ratio_optimal(Hgrid_study[k1])*np.ones(len(p_grid)), 'g--']
  
# # image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
# # image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  




# preions = ionisation_ratio(A0,XUV_energy_pp[0][:,1,k1],Hgrid_study[k1],omegaSI)
# image.sf[7].args = [p_grid, 100*preions[0][:], 'c--'];
# image.sf[8].args = [p_grid, 100*preions[1][:], 'c--'];
# image.sf[7].method = None
# image.sf[8].method = None

## custom legend
# custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
#                 Line2D([0], [0], color="tab:blue", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
#                 Line2D([0], [0], color="tab:green", lw=3),                
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

custom_lines = [Line2D([1], [0], color="k"),
                Line2D([0], [0], color="tab:grey", linestyle=":"),
                Line2D([0], [0], color="b"),
                Line2D([0], [0], color="tab:grey", linestyle="--"),
                Line2D([0], [0], color="r"),                
                Line2D([0], [0], color="g", linestyle="--")]

# ax.legend(custom_lines, [I0s_leg[0],
#                          pressures_leg[0],
#                          I0s_leg[1],
#                          pressures_leg[1],
#                          I0s_leg[2],
#                          pressures_leg[2]],
#           loc=1, ncol=3)

image.legend_args = [custom_lines,['no preion', r'$\eta_0$', '40 A', r'$\eta_{las.}$','50 A', '$\eta_{opt.}$']] 
image.legend_kwargs = {'loc': 1, 'ncol': 3}

# image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
image.xlabel = 'p [mbar]'; image.ylabel = 'I_XUV [arb. u.]'; image.right_ylabel = 'ionisation [%]'


image.title = 'H'+str(Hgrid_study[k1]) + ', T_discharge/2'

image.savefig_args = ['compare2.pdf']
image.savefig_kwargs = {'dpi' : 600,'bbox_inches' : 'tight'}

pp.plot_preset(image)



###################################################################################################################
###################################################################################################################
###################################################################################################################


# ## plot H17 ionisation analyses 40 A
k1 = 1 # index to access given harmonic plot


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(11)]



image.sf[0].args = [p_grid, XUV_energy_pp[0][:,0,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
image.sf[1].args = [p_grid, XUV_energy_pp[0][:,1,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}

# image.sf[2].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge/2'}
# image.sf[3].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.max(XUV_energy_pp[1][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}



A0 = A_norm(np.mean(XUV_energy_pp[0][:,0,k1]),Hgrid_study[k1],omegaSI)
image.sf[2].args = [p_grid, IntensXUV(1e-2*ionisations[0]['half_init'],Hgrid_study[k1],omegaSI,A0)/np.mean(XUV_energy_pp[0][:,0,k1]),'r:'];

# image.sf[2].args = [p_grid, XUV_energy_pp[:,2,k1]/np.max(XUV_energy_pp[:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}



# image.sf[9].args = [p_grid, IntensXUV(1e-2*ionisations['half_init'],17,omegaSI,A0)/np.max(XUV_energy_pp[:,0,k1]),'g--'];
# image.sf[9].kwargs = {'label' : 'T_discharge/2 from analytical estimate'}    

for k2 in range(3,8): image.sf[k2].axis_location = 'right'
# image.sf[3].args = [p_grid, ionisations[0]['half_init'], 'b:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[3].args = [p_grid, ionisations[0]['half'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
# image.sf[5].args = [p_grid, ionisations[1]['half_init'], 'r:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
# image.sf[6].args = [p_grid, ionisations[1]['half'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 

# image.sf[7].args = [p_grid, 100*ionisation_ratio_optimal*np.ones(len(p_grid)), 'g--']
  
# # image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
# # image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  

preions = ionisation_ratio(A0,XUV_energy_pp[0][:,1,k1],Hgrid_study[k1],omegaSI)
image.sf[4].args = [p_grid, 100*preions[0][:], 'g--'];
image.sf[5].args = [p_grid, 100*preions[1][:], 'g--'];



custom_lines = [Line2D([1], [0], color="k"),
                Line2D([0], [0], color="b"),
                Line2D([0], [0], color="r", linestyle=":"),
                Line2D([0], [0], color="r", linestyle="--"),
                Line2D([0], [0], color="g", linestyle="--")]


image.legend_args = [custom_lines,['reference', 'numerical', 'estimate', 'num. ionisation','est. ionisation']]
image.legend_kwargs = {'loc': 1, 'ncol': 2}

# image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
image.xlabel = 'p [mbar]'; image.ylabel = 'I_XUV [arb. u.]'; image.right_ylabel = 'ionisation [%]'


image.title = 'H'+str(Hgrid_study[k1]) + ', 40 A, analytic benchmark'

image.savefig_args = ['compare3.pdf']
image.savefig_kwargs = {'dpi' : 600,'bbox_inches' : 'tight'}

pp.plot_preset(image)


###################################################################################################################
###################################################################################################################
###################################################################################################################

# ## plot H17 ionisation analyses 50 A
k1 = 1 # index to access given harmonic plot


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(11)]



image.sf[0].args = [p_grid, XUV_energy_pp[1][:,0,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
image.sf[1].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}

# image.sf[2].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge/2'}
# image.sf[3].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.max(XUV_energy_pp[1][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}



A0 = A_norm(np.mean(XUV_energy_pp[1][:,0,k1]),Hgrid_study[k1],omegaSI)
image.sf[2].args = [p_grid, IntensXUV(1e-2*ionisations[1]['half_init'],Hgrid_study[k1],omegaSI,A0)/np.mean(XUV_energy_pp[1][:,0,k1]),'r:'];

# image.sf[2].args = [p_grid, XUV_energy_pp[:,2,k1]/np.max(XUV_energy_pp[:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}



# image.sf[9].args = [p_grid, IntensXUV(1e-2*ionisations['half_init'],17,omegaSI,A0)/np.max(XUV_energy_pp[:,0,k1]),'g--'];
# image.sf[9].kwargs = {'label' : 'T_discharge/2 from analytical estimate'}    

for k2 in range(3,8): image.sf[k2].axis_location = 'right'
# image.sf[3].args = [p_grid, ionisations[0]['half_init'], 'b:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[3].args = [p_grid, ionisations[1]['half'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
# image.sf[5].args = [p_grid, ionisations[1]['half_init'], 'r:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
# image.sf[6].args = [p_grid, ionisations[1]['half'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 

# image.sf[7].args = [p_grid, 100*ionisation_ratio_optimal*np.ones(len(p_grid)), 'g--']
  
# # image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
# # image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  

preions = ionisation_ratio(A0,XUV_energy_pp[1][:,1,k1],Hgrid_study[k1],omegaSI)
image.sf[4].args = [p_grid, 100*preions[0][:], 'g--'];
image.sf[5].args = [p_grid, 100*preions[1][:], 'g--'];



custom_lines = [Line2D([1], [0], color="k"),
                Line2D([0], [0], color="b"),
                Line2D([0], [0], color="r", linestyle=":"),
                Line2D([0], [0], color="r", linestyle="--"),
                Line2D([0], [0], color="g", linestyle="--")]


image.legend_args = [custom_lines,['reference', 'numerical', 'estimate', 'num. ionisation','est. ionisation']]
image.legend_kwargs = {'loc': 1, 'ncol': 2}

# image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
image.xlabel = 'p [mbar]'; image.ylabel = 'I_XUV [arb. u.]'; image.right_ylabel = 'ionisation [%]'


image.title = 'H'+str(Hgrid_study[k1]) + ', 50 A, analytic benchmark'

image.savefig_args = ['compare4.pdf']
image.savefig_kwargs = {'dpi' : 600,'bbox_inches' : 'tight'}

pp.plot_preset(image)

###################################################################################################################
###################################################################################################################
###################################################################################################################


## plot H17 + high ionisations
k1 = 1 # index to access given harmonic plot


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(11)]



image.sf[0].args = [p_grid, XUV_energy_pp[0][:,0,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
image.sf[1].args = [p_grid, XUV_energy_pp[0][:,2,k1]/np.mean(XUV_energy_pp[0][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}

image.sf[2].args = [p_grid, XUV_energy_pp[1][:,2,k1]/np.mean(XUV_energy_pp[1][:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge/2'}
# image.sf[3].args = [p_grid, XUV_energy_pp[1][:,1,k1]/np.max(XUV_energy_pp[1][:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}



# A0 = A_norm(np.mean(XUV_energy_pp[0][:,0,k1]),Hgrid_study[k1],omegaSI)

# image.sf[2].args = [p_grid, XUV_energy_pp[:,2,k1]/np.max(XUV_energy_pp[:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}



# image.sf[9].args = [p_grid, IntensXUV(1e-2*ionisations['half_init'],17,omegaSI,A0)/np.max(XUV_energy_pp[:,0,k1]),'g--'];
# image.sf[9].kwargs = {'label' : 'T_discharge/2 from analytical estimate'}    

for k2 in range(3,8): image.sf[k2].axis_location = 'right'
image.sf[3].args = [p_grid, ionisations[0]['end_init'], 'b:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[4].args = [p_grid, ionisations[0]['end'], 'b--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
image.sf[5].args = [p_grid, ionisations[1]['end_init'], 'r:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[6].args = [p_grid, ionisations[1]['end'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 

image.sf[7].args = [p_grid, 100*ionisation_ratio_optimal(Hgrid_study[k1])*np.ones(len(p_grid)), 'g--']
  
# # image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
# # image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  




# preions = ionisation_ratio(A0,XUV_energy_pp[0][:,1,k1],Hgrid_study[k1],omegaSI)
# image.sf[7].args = [p_grid, 100*preions[0][:], 'c--'];
# image.sf[8].args = [p_grid, 100*preions[1][:], 'c--'];
# image.sf[7].method = None
# image.sf[8].method = None

## custom legend
# custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
#                 Line2D([0], [0], color="tab:blue", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
#                 Line2D([0], [0], color="tab:green", lw=3),                
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

custom_lines = [Line2D([1], [0], color="k"),
                Line2D([0], [0], color="tab:grey", linestyle="--"),
                Line2D([0], [0], color="b"),
                Line2D([0], [0], color="tab:grey", linestyle=":"),
                Line2D([0], [0], color="r"),                
                Line2D([0], [0], color="g", linestyle="--")]

# ax.legend(custom_lines, [I0s_leg[0],
#                          pressures_leg[0],
#                          I0s_leg[1],
#                          pressures_leg[1],
#                          I0s_leg[2],
#                          pressures_leg[2]],
#           loc=1, ncol=3)

image.legend_args = [custom_lines,['no preion', r'$\eta_0$', '40 A', r'$\eta_{las.}$','50 A', '$\eta_{opt.}$']]
image.legend_kwargs = {'loc': 1, 'ncol': 3}

# image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
image.xlabel = 'p [mbar]'; image.ylabel = 'I_XUV [arb. u.]'; image.right_ylabel = 'ionisation [%]'


image.title = 'H'+str(Hgrid_study[k1]) + ', T_discharge'

image.savefig_args = ['compare5.pdf']
image.savefig_kwargs = {'dpi' : 600,'bbox_inches' : 'tight'}

pp.plot_preset(image)


###################################################################################################################
###################################################################################################################
###################################################################################################################


## Ionisations
k1 = 1 # index to access given harmonic plot


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(11)]




image.sf[0].args = [p_grid, ionisations[0]['end_init'], 'b:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[1].args = [p_grid, ionisations[0]['end'], 'b--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
image.sf[2].args = [p_grid, ionisations[1]['end_init'], 'r:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[3].args = [p_grid, ionisations[1]['end'], 'r--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'}

image.sf[4].args = [p_grid, ionisations[0]['half_init'], 'g:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[5].args = [p_grid, ionisations[0]['half'], 'g--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'} 
image.sf[6].args = [p_grid, ionisations[1]['half_init'], 'c:']; # image.sf[3].kwargs = {'label' : 'by discharge'}   
image.sf[7].args = [p_grid, ionisations[1]['half'], 'c--']; # image.sf[4].kwargs = {'label' : 'by discharge + transient'}  

# image.sf[7].args = [p_grid, 100*ionisation_ratio_optimal(Hgrid_study[k1])*np.ones(len(p_grid)), 'g--']
  
# # image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
# # image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  




# preions = ionisation_ratio(A0,XUV_energy_pp[0][:,1,k1],Hgrid_study[k1],omegaSI)
# image.sf[7].args = [p_grid, 100*preions[0][:], 'c--'];
# image.sf[8].args = [p_grid, 100*preions[1][:], 'c--'];
# image.sf[7].method = None
# image.sf[8].method = None

## custom legend
# custom_lines = [Line2D([1], [0], color="tab:orange", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="-"),
#                 Line2D([0], [0], color="tab:blue", lw=3),
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle="--"),
#                 Line2D([0], [0], color="tab:green", lw=3),                
#                 Line2D([0], [0], color="tab:grey", lw=3, linestyle=":")]

# custom_lines = [Line2D([1], [0], color="k"),
#                 Line2D([0], [0], color="tab:grey", linestyle="--"),
#                 Line2D([0], [0], color="b"),
#                 Line2D([0], [0], color="tab:grey", linestyle=":"),
#                 Line2D([0], [0], color="r"),                
#                 Line2D([0], [0], color="g", linestyle="--")]

# ax.legend(custom_lines, [I0s_leg[0],
#                          pressures_leg[0],
#                          I0s_leg[1],
#                          pressures_leg[1],
#                          I0s_leg[2],
#                          pressures_leg[2]],
#           loc=1, ncol=3)

image.legend_args = [custom_lines,['no preion', r'$\eta_0$', '40 A', r'$\eta_{las.}$','50 A', '$\eta_{opt.}$']]
image.legend_kwargs = {'loc': 1, 'ncol': 3}

# image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
image.xlabel = 'p [mbar]';
image.ylabel = 'ionisation [%]'
# image.right_ylabel = 'ionisation [%]'


# image.title = 'H'+str(Hgrid_study[k1]) + ', T_discharge'

# image.savefig_args = ['compare5.pdf']
# image.savefig_kwargs = {'dpi' : 600}

pp.plot_preset(image)




############# PLOT SPECTRA ETC.

## plot dE/dH
image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(2)]
image.sf[0].args = [Hgrid, dE_dH[0][1,0,:]/np.max(dE_dH[0][1,0,:])]; image.sf[0].kwargs = {'label' : 'no_preion'}    
image.sf[1].args = [Hgrid, dE_dH[0][1,1,:]/np.max(dE_dH[0][1,0,:])]; image.sf[1].kwargs = {'label' : 'T_discharge/2'}

pp.plot_preset(image)



# spectra in logscale

# choices = [[0,0,0],[0,0,1],[0,2,0],[0,2,1],[0,4,0],[0,4,1]]   #  [[0,0]] # # organisation [scan, pressure, preionT]

### !!! The run showing evolution with ionisation
# choices = [[0,0,0],[0,0,1],[0,1,1],[0,2,1],[0,3,1],[0,4,1],[0,5,1]]   #  [[0,0]] # # organisation [scan, pressure, preionT]

# choices = [[0,0,0],[0,1,0],[0,2,0],[0,3,0],[0,4,0],[0,5,0]]   #  [[0,0]] # # organisation [scan, pressure, preionT]

# choices = [[1,1,2],[1,2,2],[1,3,2],[1,4,2],[1,5,2]]

### !!! The run showing evolution with pressure without pre-ionisation
# choices = [[0,0,0],[0,1,0],[0,2,0],[0,3,0],[0,4,0],[0,5,0]]

# choices = [[0,0,0],[0,2,0],[0,4,0]]

# choices = [[0,0,1],[0,0,1],[0,0,2],[1,0,1],[1,0,1],[1,0,2]]

### !!! show two cases for 20 mbar
# choices = [[0,1,0],[0,1,1]]

### !!! show two cases for 20 mbar
choices = [[0,1,0],[0,1,1],[0,1,2],[1,1,1],[1,1,2]]



apply_global_norm = True
global_norm_log = np.max(np.log10(abs(FField_FF_pp[choices[0][0]][choices[0][1],choices[0][2],:,:].T)**2))
global_norm_lin = np.max(abs(FField_FF_pp[choices[0][0]][choices[0][1],choices[0][2],:,:].T)**2)

plot_scale = 'log'
include_average_dE_dH = True
include_first_in_dE_dH_plot = True

plot_zoom = True
zoom_xlim = [15,17]; zoom_ylim = [0,0.004];



# plot all linear averaged spectra


# dE_dH_avrg = 1.*dE_dH[0][0,0,:]
# for k1 in range(1,len(p_grid)):
#     dE_dH_avrg += dE_dH[0][k1,0,:]
# dE_dH_avrg = dE_dH_avrg/len(p_grid)
dE_dH_avrg = np.mean(dE_dH[0][:,0,:], axis = 0)
dE_dH_filtered_avrg = np.mean(dE_dH_filtered[0][:,0,:], axis = 0)

if include_average_dE_dH: k_start_k_end = [1,len(choices)+1]
else: k_start_k_end = [0,len(choices)]


image = pp.figure_driver()    
image.sf = [pp.plotter() for k2 in range(k_start_k_end[1])]
if include_average_dE_dH:
    norm = np.max(dE_dH_filtered_avrg)
    image.sf[0].args = [Hgrid, dE_dH_filtered_avrg/norm]
    image.sf[0].kwargs = {'label' : 'average'}
else: norm = np.max(dE_dH[choices[0][0]][choices[0][1],choices[0][2],:])    

k2 = 0    
for k1 in range(*k_start_k_end):
    # norm = np.max(dE_dH[choices[0][0]][choices[0][1],choices[0][2],:])
    image.sf[k1].args = [Hgrid, dE_dH_filtered[choices[k2][0]][choices[k2][1],choices[k2][2],:]/norm]
    image.sf[k1].kwargs = {'label' : choice_to_label(choices[k2])}
    if ((k2==0) and not(include_first_in_dE_dH_plot)): image.sf[k1].method = None
    k2 += 1

image.xlabel = 'H [-]'
image.ylabel = '$\mathrm{d}E/\mathrm{d}H$ [arb. u.]'

image.savefigs_args = [['dEdH.pdf'],['dEdH.png']]
image.savefigs_kwargs = [{'bbox_inches' : 'tight'} for k1 in range(2)]
    
pp.plot_preset(image)
    

for k1 in range(len(choices)):
    FF_spectrum_logscale = np.log10(abs(FField_FF_pp[choices[k1][0]][choices[k1][1],choices[k1][2],:,:].T)**2)
    FF_spectrum_linscale = abs(FField_FF_pp[choices[k1][0]][choices[k1][1],choices[k1][2],:,:].T)**2
    
    if apply_global_norm: FF_spectrum_logscale = FF_spectrum_logscale -  global_norm_log # normalise
    else: FF_spectrum_logscale = FF_spectrum_logscale -  np.max(FF_spectrum_logscale[:])# normalise
    vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
        
    
    # if (choices[k1][2] == 0): eta0 = 0
    # elif (choices[k1][2] == 1): eta0 = ionisations[choices[k1][0]]['half_init'][choices[k1][1]]
    # elif (choices[k1][2] == 2): eta0 = ionisations[choices[k1][0]]['end_init'][choices[k1][1]]
    
    # local_title = r'$p$='+ '{:.1f}'.format(p_grid[choices[k1][1]]) + ' mbar, ' + \
    #               r'$\eta_0$='+ '{:.1f}'.format(eta0) + ' %'
     
    # local_title = choice_to_label(choices[k1])              
    
    
    image = pp.figure_driver()    
    # image.sf.append(pp.plotter())
    image.sf.append(pp.plotter())
    image.sf[0].method = plt.pcolor
    
    if (plot_scale == 'log'):
        image.sf[0].args = [Hgrid,rgrid_FF,FF_spectrum_logscale];
        image.sf[0].kwargs = {'shading': 'auto', 'vmin': vmin, 'cmap' : 'plasma', 'rasterized' : True}          
    elif (plot_scale == 'lin'):
      if apply_global_norm:
        image.sf[0].args = [Hgrid,rgrid_FF,FF_spectrum_linscale/global_norm_lin]
      else:
        image.sf[0].args = [Hgrid,rgrid_FF,FF_spectrum_linscale/np.max(FF_spectrum_linscale)]  
      image.sf[0].kwargs = {'shading': 'auto', 'cmap' : 'plasma', 'rasterized' : True}           
        
    image.sf[0].colorbar.show = True
    
    image.xlabel = 'H [-]'; image.ylabel = 'r [m]'
    # image.title = 'Far-field spectrum, log'
    image.title = choice_to_label(choices[k1])
    
    if (plot_scale == 'log'):
        image.sf[0].colorbar.kwargs = {'label': r'$\mathrm{log}_{10}|\mathcal{E}(\omega,\rho)|^2$ [arb. u.]'}
    elif (plot_scale == 'lin'):
        image.sf[0].colorbar.kwargs = {'label': r'$|\mathcal{E}(\omega,\rho)|^2$ [arb. u.]'}
        
    fname = 'Spectrum_sim' + str(choices[k1][0]) + \
            '_press_' + str(choices[k1][1]) + \
            '_preion_' + str(choices[k1][2])
                
    image.savefigs_args = [[fname+'.pdf'],[fname+'.png']]
    image.savefigs_kwargs = [{'bbox_inches' : 'tight'} for k2 in range(2)]
    
    pp.plot_preset(image)
    
    if plot_zoom:
        image.xlim_args = [zoom_xlim]
        image.ylim_args = [zoom_ylim]
        fname2 = fname + '_zoom'
        image.savefigs_args = [[fname2+'.pdf'],[fname2+'.png']]  
        pp.plot_preset(image)        
    
    
    


    # if showplots: plt.show()
    # plt.close(fig)
    # # sys.exit()

# os.chdir(cwd)
# sys.exit()


# fig, ax = plt.subplots()     
# plt.plot(p_grid, XUV_energy_pp[:,0,1],'k',label = 'no preion')
# plt.plot(p_grid, XUV_energy_pp[:,1,1],'b',label = '40A')


# ax2 = ax.twinx()
# ax2.plot(p_grid, ionisations['half_init'], 'b:', label = 'by discharge') # 'b-'
# ax2.plot(p_grid, ionisations['half'], 'b--', label = 'by discharge + transient') # 'b-'

# ax.legend(loc='upper right')
# ax2.legend(loc='upper left')

# ax2.set_ylabel('ionisation [%]')
# ax.set_xlabel('p [mbar]'); ax.set_ylabel('E_XUV [arb. u.]');
# ax.set_title('Energy , H'+str(Hgrid_study[1]))
# fig.savefig('Energy_H'+str(Hgrid_study[1])+'_partial.png', dpi = 600)
# plt.show()


# for k1 in range(NH_study):
    
#     fig, ax = plt.subplots()     
#     plt.plot(p_grid, XUV_energy_pp[:,0,k1],'k',label = 'no preion')
#     plt.plot(p_grid, XUV_energy_pp[:,1,k1],'b',label = 'T_discharge/2')
#     plt.plot(p_grid, XUV_energy_pp[:,2,k1],'r',label = 'T_discharge')
    
    
    
#     ax2 = ax.twinx()
#     ax2.plot(p_grid, ionisations['half_init'], 'b:', label = 'by discharge') # 'b-'
#     ax2.plot(p_grid, ionisations['half'], 'b--', label = 'by discharge + transient') # 'b-'
#     ax2.plot(p_grid, ionisations['end_init'], 'r:') # 'b-'
#     ax2.plot(p_grid, ionisations['end'], 'r--') # 'b-'
    
#     ax.legend(loc='upper right')
#     ax2.legend(loc='upper left')
    
#     ax2.set_ylabel('ionisation [%]')
#     ax.set_xlabel('p [mbar]'); ax.set_ylabel('E_XUV [arb. u.]');
#     ax.set_title('Energy , H'+str(Hgrid_study[k1]))
#     fig.savefig('Energy_H'+str(Hgrid_study[k1])+'.png', dpi = 600)
#     plt.show()
    
#     # preset plot
#     image = pp.figure_driver()    
#     image.sf = [pp.plotter() for k2 in range(11)]
    
    
    
#     image.sf[0].args = [p_grid, XUV_energy_pp[:,0,k1]/np.max(XUV_energy_pp[:,0,k1]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
#     image.sf[1].args = [p_grid, XUV_energy_pp[:,1,k1]/np.max(XUV_energy_pp[:,0,k1]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}
#     image.sf[2].args = [p_grid, XUV_energy_pp[:,2,k1]/np.max(XUV_energy_pp[:,0,k1]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}
    
#     A0 = A_norm(np.mean(XUV_energy_pp[:,0,k1]),Hgrid_study[k1],omegaSI)
    
#     image.sf[9].args = [p_grid, IntensXUV(1e-2*ionisations['half_init'],17,omegaSI,A0)/np.max(XUV_energy_pp[:,0,k1]),'g--'];
#     image.sf[9].kwargs = {'label' : 'T_discharge/2 from analytical estimate'}    
    
#     for k2 in range(3,9): image.sf[k2].axis_location = 'right'
#     image.sf[3].args = [p_grid, ionisations['half_init'], 'b:']; image.sf[3].kwargs = {'label' : 'by discharge'}   
#     image.sf[4].args = [p_grid, ionisations['half'], 'b--']; image.sf[4].kwargs = {'label' : 'by discharge + transient'}   
#     image.sf[5].args = [p_grid, ionisations['end_init'], 'r:']; # image.sf[5].kwargs = {'label' : 'by discharge'}   
#     image.sf[6].args = [p_grid, ionisations['end'], 'r--']; # image.sf[6].kwargs = {'label' : 'by discharge'}  
    

    

#     preions = ionisation_ratio(A0,XUV_energy_pp[:,1,k1],Hgrid_study[k1],omegaSI)
#     image.sf[7].args = [p_grid, 100*preions[0][:], 'c--'];
#     image.sf[8].args = [p_grid, 100*preions[1][:], 'c--'];
#     # image.sf[7].method = None
#     # image.sf[8].method = None

    
    
#     image.legend_kwargs = {'loc':'upper right'}; image.right_axis_legend_kwargs = {'loc':'upper left'} 
#     image.xlabel = 'p [mbar]'; image.ylabel = 'E_XUV [arb. u.]'; image.right_ylabel = 'ionisation [%]'
#     image.title = 'Energy , H'+str(Hgrid_study[k1])
    
#     pp.plot_preset(image)
    
    
#     fig, ax = plt.subplots()     
#     plt.plot(p_grid, XUV_energy_pp[:,1,k1]/XUV_energy_pp[:,0,k1],label = 'T_discharge/2')
#     plt.plot(p_grid, XUV_energy_pp[:,2,k1]/XUV_energy_pp[:,0,k1],label = 'T_discharge')    
#     ax.legend()
#     ax.set_xlabel('p [mbar]'); ax.set_ylabel('E (end) /E0 [-]');
#     ax.set_title('Amplification, H'+str(Hgrid_study[k1]))
#     fig.savefig('Amplification_H'+str(Hgrid_study[k1])+'.png', dpi = 600)
#     plt.show()

# # fig, ax = plt.subplots()     
# # plt.plot(p_grid, XUV_energy_pp[:,0,1])
# # plt.plot(p_grid, XUV_energy_pp[:,1,1])
# # plt.show()

# # fig, ax = plt.subplots()     
# # plt.plot(p_grid, XUV_energy_pp[:,0,2])
# # plt.plot(p_grid, XUV_energy_pp[:,1,2])
# # plt.show()

# # kp = 4
# for kp in range(N_press):
#     fig, ax = plt.subplots()     
#     plt.plot(Hgrid, dE_dH[kp,0,:],'k',label = 'no preion')
#     plt.plot(Hgrid, dE_dH[kp,1,:],'b',label = 'T_discharge/2')
#     plt.plot(Hgrid, dE_dH[kp,2,:],'r',label = 'T_discharge/2')
#     ax.legend()
#     ax.set_xlabel('H [-]'); ax.set_ylabel('dE/dH [arb. u.]');
#     ax.set_title('dE/dH, p = '+str(p_grid[kp])+' mbar')
#     fig.savefig('dE_dH_p_'+str(p_grid[kp])+'.png', dpi = 600)
#     plt.show()

# for kp in range(N_press):
#     image = pp.figure_driver()    
#     # image.sf.append(pp.plotter())
#     image.sf = [pp.plotter() for k1 in range(3)]
#     image.sf[0].args = [Hgrid, dE_dH[kp,0,:]/np.max(dE_dH[kp,0,:]),'k']; image.sf[0].kwargs = {'label' : 'no_preion'}    
#     image.sf[1].args = [Hgrid, dE_dH[kp,1,:]/np.max(dE_dH[kp,0,:]),'b']; image.sf[1].kwargs = {'label' : 'T_discharge/2'}
#     image.sf[2].args = [Hgrid, dE_dH[kp,2,:]/np.max(dE_dH[kp,0,:]),'r']; image.sf[2].kwargs = {'label' : 'T_discharge'}
    
#     image.xlabel = 'H [-]'; image.ylabel = 'dE/dH [arb. u.]'
#     image.title = 'dE/dH, p = '+str(p_grid[kp])+' mbar v2'
    
#     pp.plot_preset(image)
    
# # sys.exit()

# fig, ax = plt.subplots()     
# plt.plot(p_grid, XUV_energy_pp[:,0,1]/np.max(abs(XUV_energy_pp[:,0,1])))
# # plt.plot(Hgrid, dE_dH[kp,1,:])
# # ax.set_xlabel('H [-]'); ax.set_ylabel('dE/dH [arb. u.]');
# # ax.set_title('dE/dH, p = '+str(p_grid[kp])+' mbar')
# # fig.savefig('dE_dH_p_'+str(p_grid[kp])+'.png', dpi = 600)
# plt.show()


   
   


# kp = 4
# kpre = 0

# kp_kpre = [[0,0]] # [[0,0],[0,1],[2,0],[2,1],[4,0],[4,1]]

# for k1 in range(len(kp_kpre)):
#     fig, ax = plt.subplots()   
#     FF_spectrum_logscale = np.log10(abs(FField_FF_pp[kp_kpre[k1][0],kp_kpre[k1][1],:,:].T)**2);
#     vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
#     # map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
#     map1 = ax.pcolor(Hgrid,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
#     # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
#     fig.colorbar(map1)
#     plt.title('Far-field spectrum, log')
#     plt.xlabel('H [-]')
#     plt.ylabel('r [m]')
#     plt.show()
#     # if showplots: plt.show()
#     # plt.close(fig)
#     # # sys.exit()

# for k1 in range(len(kp_kpre)):
#     FF_spectrum_logscale = np.log10(abs(FField_FF_pp[kp_kpre[k1][0],kp_kpre[k1][1],:,:].T)**2);
#     vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
        
    
#     image = pp.figure_driver()    
#     # image.sf.append(pp.plotter())
#     image.sf.append(pp.plotter())
#     image.sf[0].method = plt.pcolor
#     image.sf[0].args = [Hgrid,rgrid_FF,FF_spectrum_logscale];
#     image.sf[0].kwargs = {'shading': 'auto', 'vmin': vmin}
#     image.sf[0].colorbar = True
    
#     image.xlabel = 'H [-]'; image.ylabel = 'r [m]'
#     image.title = 'Far-field spectrum, log'
    
#     pp.plot_preset(image)
    
#     # if showplots: plt.show()
#     # plt.close(fig)
#     # # sys.exit()
    
    
# # # vmin = np.max(np.log(Gaborr))-6.
# # fig, ax = plt.subplots()   
# # # FF_spectrum_logscale = np.log10(abs(FField_FF.T)**2);
# # # vmin = np.max(FF_spectrum_logscale)-FF_orders_plot
# # # map1 = ax.pcolor(Hgrid_select,rgrid_FF,FF_spectrum_logscale, shading='auto',vmin=vmin)
# # map1 = ax.pcolor(Hgrid,rgrid_FF,abs(FField_FF.T)**2, shading='auto')
# # # plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# # fig.colorbar(map1)
# # plt.title('Far-field spectrum (30 cm), integrated')
# # plt.xlabel('H [-]')
# # plt.ylabel('r [m]')
# # plt.show()
# # # if showplots: plt.show()
# # # plt.close(fig)
# # # sys.exit()

# out_h5name = 'XUV_gains.h5'
# with h5py.File(out_h5name,'w') as OutFile: # this file contains numerical analyses
#     OutFile.create_dataset('XUV_FField',
#                                            data = np.stack((FField_FF_pp.real, FField_FF_pp.imag),axis=-1)
#                            )
#     OutFile.create_dataset('Hgrid', data = Hgrid)
#     OutFile.create_dataset('Hgrid_study', data = Hgrid_study)
#     OutFile.create_dataset('rgrid_FF', data = rgrid_FF)
#     OutFile.create_dataset('pressure_grid', data = p_grid)
#     # OutFile.create_dataset('preionisation_grid', data = preion_grid)
#     OutFile.create_dataset('XUV_energy_integrated', data = XUV_energy_pp)
#     OutFile.create_dataset('dE_dH', data = dE_dH)
    
os.chdir(cwd)