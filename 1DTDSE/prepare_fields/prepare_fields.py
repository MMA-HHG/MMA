import numpy as np
import os
import time
import shutil
import h5py
import sys
import units
import mynumerics as mn
import re
# import copy
import glob
import plot_presets as pp

arguments = sys.argv
showplots = not('-nodisplay' in arguments)




mypulse = mn.pulse_types('sin2')        
    

## testplot
if showplots:
    pass
    # image = pp.figure_driver()    
    # image.sf = [pp.plotter() for k1 in range(1)]
    # image.sf[0].args = [tgrid, Efield]
    # pp.plot_preset(image)   


# read parameters
inputfilename = 'TDSE_create_fields.inp'
with open(inputfilename, 'r') as InputMP:
    myparams3 = mn.parameters_selector(*mn.multiparameters_lines2dict( InputMP.readlines() ))

# for k1 in range(myparams3.N_combinations):
#     print(myparams3.ret(k1))
## store

tgrid = mypulse.construct_tgrid(myparams3.param_grids['dt'],
                                myparams3.param_grids['T_FWHM'])

## list of fields
Efields = []
for k1 in range(myparams3.N_combinations):
    Efields.append(
        mypulse.pulse(tgrid,
                      *mypulse.inputs_converter(
                                                *myparams3.ret(k1), 
                                                given_inps=myparams3.assumed_output_order
                                                )                
                      )
                    )


# ## testplot
# if showplots:
#     # pass
#     image = pp.figure_driver()    
#     image.sf = [pp.plotter() for k1 in range(myparams3.N_combinations)]
#     for k1 in range(myparams3.N_combinations):
#         image.sf[k1].args = [tgrid, Efields[k1]]
#     pp.plot_preset(image)  
    
    
out_h5name = 'field_input.h5'
# omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')

with h5py.File(out_h5name,'w') as OutFile:
    mn.adddataset(OutFile, 'IRField/tgrid', tgrid , '[a.u.]' )    
    dset=OutFile.create_dataset('IRField/Eields_table',(myparams3.N_combinations,len(tgrid)), 'd')
    dset.attrs['units']=np.string_('[a.u.]')
    for k1 in range(myparams3.N_combinations):
        dset[k1,:] = mypulse.pulse(tgrid,
                                   *mypulse.inputs_converter(
                                             *myparams3.ret(k1), 
                                             given_inps=myparams3.assumed_output_order
                                             )                
                               )
    
    grp = OutFile.create_group('params')
    
    myparams3.store_to_h5(grp)
    
    # image = pp.figure_driver()    
    # image.sf = [pp.plotter() for k1 in range(5)]

    # image.sf[0].args = [tgrid, OutFile['IRField/Eields_table'][5,:]]
    # image.sf[1].args = [tgrid, OutFile['IRField/Eields_table'][6,:]]
    # image.sf[3].args = [tgrid, OutFile['IRField/Eields_table'][7,:]]
    # image.sf[4].args = [tgrid, OutFile['IRField/Eields_table'][9500,:]]
    # pp.plot_preset(image)
    # Store param grid


# ## testplot
# omegac = 2.0*np.arccos(2**(-0.25))/(tFWHMSI/units.TIMEau)
# tgrid = np.linspace(0,np.pi/omegac,Nt)

# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(1)]
# image.sf[0].args = [tgrid, sin2pulse(tgrid,mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau'),omegac,E0_max,0.0)]
# pp.plot_preset(image)   



# out_h5name = 'fields_table.h5'
# omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')

# with h5py.File(out_h5name,'w') as OutFile:
#     shape = [N,Nt]
#     dset = OutFile.create_dataset('fields_list', shape, precision)
#     dset.attrs['units'] = np.string_('[a.u.]')
    
#     for k1 in range(N):
#         MultInd = np.unravel_index(k1,param_dims)
#         dset[k1,:] = sin2pulse(tgrid,omega0,omegac,param_grids[1][MultInd[1]],param_grids[0][MultInd[0]])

#     OutFile.create_dataset('param_list',data=np.string_(varying_list)) 
    
#     mn.adddataset(OutFile, 'tgrid', tgrid , '[a.u.]' )
    
#     # store grids
#     for k1 in range(len(varying_list)):
#         # OutFile.create_dataset('param_'+str(k1), data=param_grids[k1] )
#         mn.adddataset(OutFile, 'param_'+str(k1), param_grids[k1] , parameters[varying_list[k1]][3] )
#         # .create_dataset('param_'+str(k1), data=param_grids[k1] )

print('done')


### Initial point for simulations

# Conditions are:
# No envelop temporally
# Integration on 1 optical cycle
# Trajectory selection (Lewenstein model)
# I0=1e14 W/cm2 // future 7e13
# w0=25 microns // future 10
# Lambda0=800 nm // future 1030 nm
# Spatial mode= TEM00 gaussian
# Pressure=0.1 bar
# Atom: argon// future xenon
# Absorption=1cm response
# Length of medium= 50 micron
# Position of target= @ focus and before z=-100 micron and after z=+100
# Calculation plan: @ farfield this 50 cm from medium output

# + times

# 30 fs - 800 nm
# 200 fs - 1030 nm