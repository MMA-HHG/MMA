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



groupname = arguments[arguments.index("-g")+1] if ("-g" in arguments) else 'fields_list'
out_h5name = arguments[arguments.index("-ohdf5")+1] if ("-ohdf5" in arguments) else 'field_input.h5'
inputfilename = arguments[arguments.index("-i")+1] if ("-i" in arguments) else 'TDSE_create_fields.inp'
groupname_params = arguments[arguments.index("-g_params")+1] if ("-g_params" in arguments) else 'grids_for_scan'

    


# read parameters

with open(inputfilename, 'r') as InputMP:
    myparams3 = mn.parameters_selector(*mn.multiparameters_lines2dict( InputMP.readlines() ))

mypulse = mn.pulse_types(myparams3.param_grids['pulse_type'])    

tgrid = mypulse.construct_tgrid(myparams3.param_grids['dt'],
                                myparams3.param_grids['T_FWHM'])

# ## list of fields
# Efields = []
# for k1 in range(myparams3.N_combinations):
#     Efields.append(
#         mypulse.pulse(tgrid,
#                       *mypulse.inputs_converter(
#                                                 *myparams3.ret(k1), 
#                                                 given_inps=myparams3.assumed_output_order
#                                                 )                
#                       )
#                     )
  
    
    


with h5py.File(out_h5name,'w') as OutFile:
    mn.adddataset(OutFile, groupname+'/tgrid', tgrid , '[a.u.]' )    
    dset=OutFile.create_dataset(groupname+'/Efields_table',(myparams3.N_combinations,len(tgrid)), 'd')
    dset.attrs['units']=np.string_('[a.u.]')
    for k1 in range(myparams3.N_combinations):
        dset[k1,:] = mypulse.pulse(tgrid,
                                   *mypulse.inputs_converter(
                                             *myparams3.ret(k1), 
                                             given_inps=myparams3.assumed_output_order
                                             )                
                               )
    
    grp = OutFile.create_group(groupname_params)    
    myparams3.store_to_h5(grp)
    


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