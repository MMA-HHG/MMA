import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
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

import plot_presets as pp


   
cwd = os.getcwd()

vacuum_frame = True

# base_path = os.path.join("C:\data", "JZ","density_mod")
base_path = os.path.join("D:\data", "JZ","density_mod")

sims_to_analyse = []

# sims_to_analyse.append( os.path.join("halving_test1","reference") )
# sims_to_analyse.append( os.path.join("halving_test1","half_simple") )
# sims_to_analyse.append( os.path.join("halving_test1","half_table") )
# sims_to_analyse.append( os.path.join("halving_test2","reference") )
# sims_to_analyse.append( os.path.join("halving_test2","half_simple") )
# sims_to_analyse.append( os.path.join("halving_test2","half_table") )

# sims_to_analyse.append( os.path.join("100","100test1","simple") )
# sims_to_analyse.append( os.path.join("100","100test1","table") )


# sims_to_analyse.append( os.path.join("100","100test2","simple_short") )
# sims_to_analyse.append( os.path.join("100","100test2","table_short") )

# sims_to_analyse.append( os.path.join("100","100test2","simple") )
# sims_to_analyse.append( os.path.join("100","100test2","table") )

# sims_to_analyse.append( os.path.join("series6","test1_modT2") )
# sims_to_analyse.append( os.path.join("series3","test1_mod_incT2") )

# sims_to_analyse.append( os.path.join("series6","test1_GaussT2") )
# sims_to_analyse.append( os.path.join("series6","test2_GaussT2") )

# sims_to_analyse.append( os.path.join("series6","test1_modT2") )
# sims_to_analyse.append( os.path.join("series6","test1_mod_incT2") )
# sims_to_analyse.append( os.path.join("series6","test1_mod_decT2") )

# sims_to_analyse.append( os.path.join("series7","test1_modT2") )
# sims_to_analyse.append( os.path.join("series7","test1_mod_incT2") )
# sims_to_analyse.append( os.path.join("dec1","test1_mod_decT2") )

# sims_to_analyse.append( os.path.join("series8","test1_modT1") )
# sims_to_analyse.append( os.path.join("series8","test1_modT2") )
# sims_to_analyse.append( os.path.join("series8","test1_modT3") )

## !!!!
# sims_to_analyse.append( os.path.join("series8","test3_modT2") )
# sims_to_analyse.append( os.path.join("series8","test3_mod_incT2") )
# sims_to_analyse.append( os.path.join("series8","test3_mod_decT2") )

# sims_to_analyse.append( os.path.join("series7","test1_modT2") )
# sims_to_analyse.append( os.path.join("series8","test3_modT2") )
# sims_to_analyse.append( os.path.join("series8","test4_modT2") )


# sims_to_analyse.append( os.path.join("series7","test1_mod_decT2") )
# sims_to_analyse.append( os.path.join("series8","test3_mod_decT2") )
# sims_to_analyse.append( os.path.join("series8","test4_mod_decT2") )

sims_to_analyse.append( os.path.join("series9","test7_modT2") )
sims_to_analyse.append( os.path.join("series9","test7_mod_incT2") )
sims_to_analyse.append( os.path.join("series9","test7_mod_decT2") )

# sims_to_analyse.append( os.path.join("series9","test1_modT2") )
# sims_to_analyse.append( os.path.join("series9","test7_modT2") )



# sims_to_analyse.append( os.path.join("series9","test1_modT2") )
# sims_to_analyse.append( os.path.join("series9","oldcodeT2") )


# sims_to_analyse=[os.path.join("C:\data", "JZ","density_mod","100","100test1","simple"),
#                  os.path.join("E:\data", "JZ","density_mod","series3","test1_modT2")]


results_filename = "results.h5"

plot_vacuum = False
plot_onax = False
plot_density = False
fluence_analysis = False


# plot_vacuum = True
# plot_onax = True
# plot_density = True
# fluence_analysis = True



file_paths = [os.path.join(base_path,inter_path,results_filename) for inter_path in sims_to_analyse]


full_resolution = False
rmax = 130e-6 # only for analyses
dr = rmax/40.0
    
linestyles = ['','--',':','-.']



with ExitStack() as stack:   
    
    InArch = [stack.enter_context(h5py.File(fpath, 'r')) for fpath in file_paths]
    
    res = [dfC.get_data(arch, r_resolution = [full_resolution, dr, rmax]) for arch in InArch]
    Nsim = len(res) 
    
    try: dens_profile = [arch['test_density'][:] for arch in InArch]
    except: dens_profile = np.asarray([[1]])
    
    print('')
    # characteristics
    for k1 in range(Nsim):
        if hasattr(res[k1], 'density_mod_profile_mbar'):
            pressure_str = "{:.1f}".format(res[k1].density_mod_profile_mbar[0])+' mbar (table applied)'
        else:
            pressure_str = res[k1].pressure_string
            
        print('simulation '+str(k1)+'\n'+
              'intensity '+ res[k1].Intensity_entry_string+'\n'+
              'pressure '+ pressure_str+'\n'+
              'ref density', res[k1].rho0_init, '\n'+
              'group vel ', res[k1].VG_IR, '(diff to c {:.2e}'.format(abs(1e2*(res[k1].VG_IR-units.c_light)/units.c_light)) + ' %)\n'
            )


    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'endfield'
    for k1 in range(Nsim): image.sf[k1].args = [1e15*res[k1].tgrid,res[k1].E_trz[:,0,-1],linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  
    

    z_fix = 2.3e-3 # SI
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'fix-field'
    for k1 in range(Nsim):
        kz_loc = mn.FindInterval(res[k1].zgrid, z_fix)
        image.sf[k1].args = [1e15*res[k1].tgrid,res[k1].E_trz[:,0,kz_loc],linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  


    for k1 in range(Nsim): res[k1].get_plasma(InArch[k1], r_resolution = [full_resolution, dr, rmax])
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'fix-plasma'
    for k1 in range(Nsim):
        kz_loc = mn.FindInterval(res[k1].zgrid, z_fix)
        image.sf[k1].args = [1e15*res[k1].plasma.tgrid,res[k1].plasma.value_trz[:,0,kz_loc]/(res[k1].rho0_init*dens_profile[k1][0,0]),linestyles[k1%len(linestyles)]]
            
                        
    pp.plot_preset(image)  
    
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'startfield'
    for k1 in range(Nsim): image.sf[k1].args = [1e15*res[k1].tgrid,res[k1].E_trz[:,0,0],linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  
   
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'startfield + 1'
    for k1 in range(Nsim): image.sf[k1].args = [1e15*res[k1].tgrid,res[k1].E_trz[:,0,1],linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  
    
    kz = mn.FindInterval(res[k1].zgrid, 2e-3)
    print('kz', kz)
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(res[0].Nz)]
    image.title = 'field'
    for k1 in range(20,60): image.sf[k1].args = [1e15*res[0].tgrid,res[0].E_trz[:,0,k1]]                
    pp.plot_preset(image)  
    


    for k1 in range(Nsim): res[k1].compute_spectrum()
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'endspect'
    for k1 in range(Nsim): image.sf[k1].args = [res[k1].ogrid,np.abs(res[k1].FE_trz[:,0,-1])**2,linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  
    

    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'endspect log'
    for k1 in range(Nsim):
        image.sf[k1].args = [res[k1].ogrid,np.abs(res[k1].FE_trz[:,0,-1])**2,linestyles[k1%len(linestyles)]]        
        image.sf[k1].method = plt.semilogy        
    pp.plot_preset(image)  
    
 
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'endspect real'
    for k1 in range(Nsim): image.sf[k1].args = [res[k1].ogrid,res[k1].FE_trz[:,0,-1].real,linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  
    
    
    
    for k1 in range(Nsim): res[k1].get_plasma(InArch[k1], r_resolution = [full_resolution, dr, rmax])
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'start ion'
    for k1 in range(Nsim): image.sf[k1].args = [1e15*res[k1].plasma.tgrid,res[k1].plasma.value_trz[:,0,0]/(res[k1].rho0_init*dens_profile[k1][0,0]),linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  

    
    for k1 in range(Nsim): res[k1].get_plasma(InArch[k1], r_resolution = [full_resolution, dr, rmax])
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    image.title = 'end ion'
    for k1 in range(Nsim): image.sf[k1].args = [1e15*res[k1].plasma.tgrid,res[k1].plasma.value_trz[:,0,-1]/(res[k1].rho0_init*dens_profile[k1][0,0]),linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  
    
  
    for k1 in range(Nsim):
        
        image = pp.figure_driver()  
        image.title = 'sim'+str(k1)+' plasma end'
        image.xlabel = 'z [mm]'
        image.ylabel = r'r [$\mu$m]'
        
        image.sf = [pp.plotter() for k1 in range(2)]
        
        image.sf[0].method = plt.pcolormesh    
        image.sf[0].args = [1e3*res[k1].plasma.zgrid,1e6*res[k1].plasma.rgrid, res[k1].plasma.value_trz[-1,:,:]]    
        image.sf[0].kwargs = {'shading' : 'auto', 'cmap' : 'plasma'}  
       
        
        image.sf[0].colorbar.show = True  
        # image.sf[0].colorbar.kwargs = {'label': r'$L_{coh}$ [mm]'}       
        
        
        # image.annotation = [['(b)'],
        #                 {'xy' : (0.025, .9),
        #                  'xycoords' : 'axes fraction',
        #                  'color' : 'w'}]
        pp.plot_preset(image)
  
    
  
    if fluence_analysis:
        for k1 in range(Nsim):
            res[k1].get_Fluence(InArch[k1], fluence_source='computed')
            
            image = pp.figure_driver()  
            image.title = 'sim'+str(k1)+' fluence'
            image.xlabel = 'z [mm]'
            image.ylabel = r'r [$\mu$m]'
            
            image.sf = [pp.plotter() for k1 in range(2)]
            
            image.sf[0].method = plt.pcolormesh    
            image.sf[0].args = [1e3*res[k1].Fluence.zgrid,1e6*res[k1].Fluence.rgrid, res[k1].Fluence.value]    
            image.sf[0].kwargs = {'shading' : 'auto', 'cmap' : 'plasma'}  
           
            
            image.sf[0].colorbar.show = True  
            # image.sf[0].colorbar.kwargs = {'label': r'$L_{coh}$ [mm]'}       
            
            
            # image.annotation = [['(b)'],
            #                 {'xy' : (0.025, .9),
            #                  'xycoords' : 'axes fraction',
            #                  'color' : 'w'}]
            pp.plot_preset(image)
  
    if plot_onax:
        for k1 in range(Nsim):
            res[k1].complexify_envel(output='add')

            image = pp.figure_driver()  
            image.title = 'sim'+str(k1)+' onax |envel|'
            image.ylabel = 't [fs]'
            image.xlabel = 'z [mm]'
            
            image.sf = [pp.plotter() for k1 in range(2)]
            
            image.sf[0].method = plt.pcolormesh    
            image.sf[0].args = [1e3*res[k1].zgrid, 1e15*res[k1].tgrid, np.abs(np.squeeze(res[k1].E_trz_cmplx_envel[:,0,:]))]    
            image.sf[0].kwargs = {'shading' : 'auto', 'cmap' : 'plasma'}           
            
            image.sf[0].colorbar.show = True  
            pp.plot_preset(image)


            image = pp.figure_driver()  
            image.title = 'sim'+str(k1)+' onax field'
            image.ylabel = 't [fs]'
            image.xlabel = 'z [mm]'
            
            image.sf = [pp.plotter() for k1 in range(2)]
            
            image.sf[0].method = plt.pcolormesh    
            image.sf[0].args = [1e3*res[k1].zgrid,1e15*res[k1].tgrid, np.abs(np.squeeze(res[k1].E_trz[:,0,:]))]    
            image.sf[0].kwargs = {'shading' : 'auto', 'cmap' : 'plasma'}           
            
            image.sf[0].colorbar.show = True  
            pp.plot_preset(image)
            
            
    if plot_vacuum:
        # for foo in res: foo.vacuum_shift()
        for k1 in range(Nsim): res[k1].vacuum_shift(output='add')
        
        image = pp.figure_driver()
        
        image.title = 'shift to vacuum frame'
        image.xlabel = 't [fs]'
            
        image.sf = [pp.plotter() for k1 in range(2*Nsim)]
        
        for k1 in range(Nsim): image.sf[k1].args = [1e15*res[k1].tgrid,res[k1].E_trz[:,0,-1],linestyles[k1%len(linestyles)]]  
        for k1 in range(Nsim): image.sf[k1].kwargs = {'label' : 'orig '+str(k1)}
        
        for k1 in range(Nsim): image.sf[k1+Nsim].args = [1e15*res[k1].tgrid,res[k1].E_trz_vac[:,0,-1],linestyles[(k1+Nsim)%len(linestyles)]]               
        for k1 in range(Nsim): image.sf[k1+Nsim].kwargs = {'label' : 'vacuum '+str(k1)}
        
        pp.plot_preset(image)  
        
        
        # zgrid = [arch['/outputs/zgrid'][:] for arch in InArch]; Nz = [len(foo) for foo in zgrid]
        # tgrid = [arch['/outputs/tgrid'][:] for arch in InArch]
        # E_slice = [arch['/outputs/output_field'][N-1,:,0] for N, arch in zip(Nz, InArch)]



    if plot_density:
        image = pp.figure_driver()
        image.sf = [pp.plotter() for k1 in range(Nsim)]
        image.title = 'density profile (onax)'
        for k1 in range(Nsim): image.sf[k1].args = [1e3*res[k1].zgrid,dens_profile[k1][0,:],linestyles[k1%len(linestyles)]]                
        pp.plot_preset(image)  



# for k1 in range(len(zgrid)): print('sim'+str(k1)+ ' zend', zgrid[k1][-1])

# Nsim = len(zgrid)   

# ogrid, FE_slice = zip(*[mn.fft_t(tgrid[k1], E_slice[k1])[0:2] for k1 in range(Nsim)])


# linestyles = ['','--',':','-.']


# image = pp.figure_driver()
# image.sf = [pp.plotter() for k1 in range(Nsim)]

# for k1 in range(Nsim): image.sf[k1].args = [tgrid[k1],E_slice[k1],linestyles[k1%len(linestyles)]]
            
# pp.plot_preset(image)      



# image = pp.figure_driver()
# image.sf = [pp.plotter() for k1 in range(Nsim)]
# for k1 in range(Nsim):
#     image.sf[k1].args = [ogrid[k1],abs(FE_slice[k1])**2,linestyles[k1%len(linestyles)]]
# pp.plot_preset(image)

# image = pp.figure_driver()
# image.sf = [pp.plotter() for k1 in range(Nsim)]
# for k1 in range(Nsim):
#     image.sf[k1].args = [ogrid[k1],abs(FE_slice[k1])**2,linestyles[k1%len(linestyles)]]
#     image.sf[k1].method = plt.semilogy
# pp.plot_preset(image)

# image = pp.figure_driver()
# image.sf = [pp.plotter() for k1 in range(Nsim)]
# for k1 in range(Nsim): image.sf[k1].args = [ogrid[k1],FE_slice[k1].real,linestyles[k1%len(linestyles)]]          
# pp.plot_preset(image)

