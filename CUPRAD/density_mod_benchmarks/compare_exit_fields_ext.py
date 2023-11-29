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

base_path = os.path.join("C:\data", "JZ","density_mod")

sims_to_analyse = []

# sims_to_analyse.append( os.path.join("halving_test1","reference") )
# sims_to_analyse.append( os.path.join("halving_test1","half_simple") )
# sims_to_analyse.append( os.path.join("halving_test1","half_table") )
# sims_to_analyse.append( os.path.join("halving_test2","reference") )
# sims_to_analyse.append( os.path.join("halving_test2","half_simple") )
# sims_to_analyse.append( os.path.join("halving_test2","half_table") )

sims_to_analyse.append( os.path.join("100","100test1","simple") )
sims_to_analyse.append( os.path.join("100","100test1","table") )

# sims_to_analyse.append( os.path.join("100","100test2","simple") )
# sims_to_analyse.append( os.path.join("100","100test2","table") )

results_filename = "results.h5"


file_paths = [os.path.join(base_path,inter_path,results_filename) for inter_path in sims_to_analyse]


full_resolution = False
rmax = 130e-6 # only for analyses
dr = rmax/40.0
    
linestyles = ['','--',':','-.']

with ExitStack() as stack:   
    
    InArch = [stack.enter_context(h5py.File(fpath, 'r')) for fpath in file_paths]
    
    res = [dfC.get_data(arch, r_resolution = [full_resolution, dr, rmax]) for arch in InArch]
    Nsim = len(res) 

    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(Nsim)]
    for k1 in range(Nsim): image.sf[k1].args = [res[k1].tgrid,res[k1].E_trz[:,0,-1],linestyles[k1%len(linestyles)]]                
    pp.plot_preset(image)  
    
    # for foo in res: foo.vacuum_shift()
    for k1 in range(Nsim): res[k1].vacuum_shift(output='add')
    
    image = pp.figure_driver()
    image.sf = [pp.plotter() for k1 in range(2*Nsim)]
    
    for k1 in range(Nsim): image.sf[k1].args = [res[k1].tgrid,res[k1].E_trz[:,0,-1],linestyles[k1%len(linestyles)]]  
    for k1 in range(Nsim): image.sf[k1].kwargs = {'label' : 'orig '+str(k1)}
    
    for k1 in range(Nsim): image.sf[k1+Nsim].args = [res[k1].tgrid,res[k1].E_trz_vac[:,0,-1],linestyles[(k1+Nsim)%len(linestyles)]]               
    for k1 in range(Nsim): image.sf[k1+Nsim].kwargs = {'label' : 'vacuum '+str(k1)}
    
    pp.plot_preset(image)  
    
    
    # zgrid = [arch['/outputs/zgrid'][:] for arch in InArch]; Nz = [len(foo) for foo in zgrid]
    # tgrid = [arch['/outputs/tgrid'][:] for arch in InArch]
    # E_slice = [arch['/outputs/output_field'][N-1,:,0] for N, arch in zip(Nz, InArch)]







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

