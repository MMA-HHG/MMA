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

I0_max = 1e18 # SI
lambdaSI = 800e-9
tFWHMSI = 30e-15

E0_max = np.sqrt(I0_max/units.INTENSITYau)
NE = 100

Nt = 1000;
precision = 'd'

# # read parameters
# inputfilename = 'TDSE_create_fields.inp'
# with open(inputfilename, 'r') as InputMP:
#     lines = InputMP.readlines()
    

def multiparameters_lines2dict(lines):
    processing_fixed_inputs = False;  processing_varying_inputs = False
    varying_params = []; fixed_params = []; params_dict = {}
    params_dict['_units']={}
    for line in lines:
        sep_line = line.split()  # separate the line
        
        if ((len(sep_line) == 0) or (sep_line[0] == '#') or (sep_line[0] == '##')):
            pass  # print('empty or commented line')
        else:
          if (sep_line[0] == '$fixed'):
              processing_fixed_inputs = True; processing_varying_inputs = False
          elif (sep_line[0] == '$varying'):
              processing_fixed_inputs = False; processing_varying_inputs = True
          elif processing_fixed_inputs: 
              fixed_params.append(sep_line[0])
              if (sep_line[2] == 'R'): params_dict[sep_line[0]] = float(sep_line[1])
              elif (sep_line[2] == 'I'): params_dict[sep_line[0]] = int(sep_line[1])
              elif (sep_line[2] == 'S'): params_dict[sep_line[0]] = sep_line[1]
              else: raise TypeError('line:' + line + '\n Specify datatype by R or I.')
              params_dict['_units'][sep_line[0]] = sep_line[3]
          elif processing_varying_inputs: 
              varying_params.append(sep_line[0]) 
              params_dict[sep_line[0]] = [float(sep_line[2]), float(sep_line[3]),
                                          int(sep_line[4])]
              params_dict['_units'][sep_line[0]] = sep_line[1]
          else: raise NotImplementedError('The input file must follow the $fixed - $varying structure now')
          
    return varying_params, fixed_params, params_dict


class parameters_selector:
    def __init__(self, 
                 varying_params, fixed_params, params_constructor,
                 dtypes={}, assumed_output_order=None):
        self.varying_params = varying_params
        self.fixed_params = fixed_params
        param_grids = {}
        self.N_combinations = 1; self.N_fixed = len(varying_params); self.N_varying = len(varying_params)
        self.varying_params_lengths = []
        
        # param = varying_list[k1]
        # N_points = parameters[param][2] + 2
        # if (parameters[param][2] == -1): param_grids.append(
        #                                             np.asarray([parameters[param][1]])
        #                                             )
        # else: param_grids.append( np.linspace( *parameters[param][:2], N_points ) )
        # N *= N_points; param_dims.append(N_points)  
        
        for param in varying_params:   
            dtype = None if not(param in dtypes.keys()) else dtypes[param]
            if (params_constructor[param][2] == -1):                
                param_grids[param] = np.asarray([params_constructor[param][1]],dtype=dtype)
            else:
                N_points = params_constructor[param][2] + 2
                param_grids[param] = np.linspace( *params_constructor[param][:2], N_points, dtype=dtype)
                self.N_combinations *= N_points
            self.varying_params_lengths.append( len(param_grids[param]) )
        for param in fixed_params:
            param_grids[param] = params_constructor[param]
        self.param_grids = param_grids
                
                    
                

        self.assumed_output_order = (list(varying_params+fixed_params)
                                     if (assumed_output_order is None) 
                                     else list(assumed_output_order))
        
        if ('_units' in params_constructor.keys()): self.units = params_constructor['_units']
        
        
    def ret(self,N,variables=None):
        MultInd = np.unravel_index(N, self.varying_params_lengths)
        output_ordered = {}; ouput_required = []
        for k1 in range(self.N_varying):
            output_ordered[self.varying_params[k1]] = self.param_grids[self.varying_params[k1]][MultInd[k1]]
        if (variables == None): variables=self.assumed_output_order
        all_possible_ouputs = output_ordered|{key:self.param_grids[key] for key in self.fixed_params}
        for var in variables:
            ouput_required.append(all_possible_ouputs[var])
        return ouput_required
    
    def store_to_h5(self,h_path):
        pass
        
    

class pulse_types:
    def __init__(self, pulse_type):
        if (pulse_type == 'sin2'): # sin^2 - envelope
            def sin2pulse(t,omega0,omegac,E0,phi0):
                phi_central = - (np.pi * omega0) / (2.0*omegac) # use the peak of the pulse as the reference for the cosine pulse
                return (t>=0) * (t<=np.pi/omegac) * E0*((np.sin(omegac*t))**2) *  np.cos(omega0*t + phi_central + phi0)
            self.pulse = sin2pulse
            # self.inputs_list_direct = ['omega0', 'omegac', 'E0', 'phi0']
            # self.inputs_list = ['lambda', 'T_FWHM', 'E0', 'phi0']
            
            def inputs_converter1(lambdaSI,tFWHMSI,E0,phi0):
                omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')
                omegac = 2.0*np.arccos(2**(-0.25))/(tFWHMSI/units.TIMEau)
                return omega0, omegac, E0, phi0
            self.inputs_converter1 = inputs_converter1
            
            # def construct_tgrid(Nt, **kwargs):
            #     omegac = 2.0*np.arccos(2**(-0.25))/(kwargs['tFWHMSI']/units.TIMEau)
            #     tgrid = np.linspace(0,np.pi/omegac,Nt)
            #     return tgrid
            # self.construct_tgrid = construct_tgrid
            
            def make_field(Nt,**kwargs):
                tgrid = construct_tgrid(Nt, **kwargs)
                Efield = sin2pulse(tgrid,*inputs_converter(**kwargs))
                return tgrid, Efield
            self.make_field = make_field
            
            def inputs_converter(*args,given_inps=['omega0','omegac','E0','phi0']):
                
                # omega0
                if ('omega0' in given_inps): omega0 = args[given_inps.index('omega0')]
                elif ('lambda0' in given_inps): omega0 = mn.ConvertPhoton(args[given_inps.index('lambda0')], 'lambdaau', 'omegaau') 
                else: raise NotImplementedError("omega0 doesn't have input")
                
                # omegac
                if ('omegac' in given_inps): omegac = args[given_inps.index('omegac')]
                elif ('T_FWHM' in given_inps):
                    omegac = 2.0*np.arccos(2.0**(-0.25))/args[given_inps.index('T_FWHM')]
                elif ('Ncyc' in given_inps):
                    omegac = omega0/(2.0*args[given_inps.index('Ncyc')])
                else: raise NotImplementedError("omegac doesn't have input")
                
                # E0
                if ('E0' in given_inps): E0 = args[given_inps.index('E0')]
                else: raise NotImplementedError("E0 doesn't have input")
                
                # phi0
                if ('phi0' in given_inps): phi0 = args[given_inps.index('phi0')]
                else: raise NotImplementedError("phi0 doesn't have input")
                
                return omega0, omegac, E0, phi0
            self.inputs_converter = inputs_converter
            
            def construct_tgrid(*args, N_points_control='dt', duration_definition='T_FWHM'):
                
                if (duration_definition=='T_FWHM'):
                    omegac = 2.0*np.arccos(2.0**(-0.25))/args[1]
                elif (duration_definition=='omegac'):
                    omegac = args[1]
                elif (duration_definition=='Ncyc'):
                    omegac = args[1]/(2.0*args[2])
                else: raise NotImplementedError("srongly sepcified 'duration_definition'")
                
                if (N_points_control=='dt'):
                    tgrid = np.arange(0, (np.pi/omegac)+args[0], args[0])
                elif (N_points_control=='Nt'):
                    tgrid = np.linspace(0,np.pi/omegac,args[0])
                else: raise NotImplementedError("srongly sepcified 'N_points_control'")
                return tgrid
            self.construct_tgrid = construct_tgrid
            

        # elif (pulse_type == 'Gaussian'): # sin^2 - envelope
        #     def Gaussian_pulse(t,omega0,tFWHM,E0,phi0):
        #         return E0* np.exp(-(2.0*np.log(2.0)*t/tFWHM)**2)*  np.cos(omega0*t + phi0)
        #     self.pulse = Gaussian_pulse
        #     self.inputs_list_direct = ['omega0', 'T_FWHM', 'E0', 'phi0']
        #     self.inputs_list = ['lambda', 'T_FWHM', 'E0', 'phi0']
            
        #     def inputs_converter(lambdaSI,tFWHMSI,E0,phi0):
        #         omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')
        #         return omega0, tFWHMSI/units.TIMEau, E0, phi0
        #     self.inputs_converter = inputs_converter
            
        #     def construct_tgrid(Nt, **kwargs):
        #         tmax = kwargs['t_expand'] * kwargs['tFWHMSI']/units.TIMEau
        #         tgrid = np.linspace(-0.5*tmax, 0.5*tmax,Nt)
        #         return tgrid
        #     self.construct_tgrid = construct_tgrid
            
        #     def make_field(Nt,**kwargs):
        #         tgrid = construct_tgrid(Nt, **kwargs)
        #         Efield = sin2pulse(tgrid,*inputs_converter(**kwargs))
        #         return tgrid, Efield
        #     self.make_field = make_field
                
        else: raise NotImplementedError('The input fiel must follow the $fixed - $varying structure now')


# dp = pulse_types(fixed['pulse_type'])

# def inputs_wrapper(k, inputs_list):    
#     MultInd = np.unravel_index(k, param_dims)
#     inputs = []
#     for inp in inputs_list:
#         if (inp in fixed_params): inputs.append(fixed[inp])
#         elif (inp in varying_list):
#             k1 = varying_list.index(inp)
#             inputs.append(param_grids[k1][MultInd[k1]])
#     return inputs

mypulse = pulse_types('sin2')        
    
tgrid = mypulse.construct_tgrid(0.5, tFWHMSI/units.TIMEau)
# Efield = mypulse.pulse(tgrid,*mypulse.inputs_converter(lambdaSI, tFWHMSI, E0_max, 0.0))


Efield = mypulse.pulse(tgrid,
                       # *mypulse.inputs_converter(lambdaSI, tFWHMSI, E0_max, 0.0)
                       *mypulse.inputs_converter(lambdaSI/units.LENGTHau, tFWHMSI/units.TIMEau, E0_max, 0.0, 
                                         given_inps=['lambda0','T_FWHM','E0','phi0'])
                       )

## testplot
if showplots:
    # pass
    image = pp.figure_driver()    
    image.sf = [pp.plotter() for k1 in range(1)]
    image.sf[0].args = [tgrid, Efield]
    pp.plot_preset(image)   


# read parameters
inputfilename = 'TDSE_create_fields.inp'
with open(inputfilename, 'r') as InputMP:
    myparams3 = parameters_selector(*multiparameters_lines2dict( InputMP.readlines() ))

for k1 in range(myparams3.N_combinations):
    print(myparams3.ret(k1))
## store


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


## testplot
if showplots:
    # pass
    image = pp.figure_driver()    
    image.sf = [pp.plotter() for k1 in range(myparams3.N_combinations)]
    for k1 in range(myparams3.N_combinations):
        image.sf[k1].args = [tgrid, Efields[k1]]
    pp.plot_preset(image)  
    
    
out_h5name = 'field_input.h5'
omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')

with h5py.File(out_h5name,'w') as OutFile:
    mn.adddataset(OutFile, 'IRField/tgrid', tgrid , '[a.u.]' )
    mn.adddataset(OutFile, 'IRField/Field', Efield , '[a.u.]' )
    mn.adddataset(OutFile, 'IRField/Eields', np.asarray(Efields) , '[a.u.]' )
    
    
    ## store params
    grp = OutFile.create_group('params')
    for k1 in range(myparams3.N_varying):
        # grp.create_dataset('param_'+str(k1),data=myparams3.param_grids[myparams3.varying_params[k1]]) 
        try:
            mn.adddataset(grp, 'param_'+str(k1) , myparams3.param_grids[myparams3.varying_params[k1]] , '['+myparams3.units[myparams3.varying_params[k1]] +']' )
        except:
            mn.adddataset(grp, 'param_'+str(k1) , myparams3.param_grids[myparams3.varying_params[k1]] , '[?]' )
    grp.create_dataset('varying_params',data=np.string_(myparams3.varying_params))
    
    for k1 in range(myparams3.N_fixed):
        # grp.create_dataset(myparams3.fixed_params[k1], data=myparams3.param_grids[myparams3.fixed_params[k1]])
        try:
            mn.adddataset(grp, myparams3.fixed_params[k1] , myparams3.param_grids[myparams3.fixed_params[k1]] , '['+myparams3.units[myparams3.fixed_params[k1]] +']' )
        except:
            mn.adddataset(grp, myparams3.fixed_params[k1] , myparams3.param_grids[myparams3.fixed_params[k1]] , '[?]' )
    grp.create_dataset('fixed_params',data=np.string_(myparams3.fixed_params))
    
    mn.adddataset(grp, 'test', ['a','aaa'] , '[a.u.]' )
    mn.adddataset(grp, 'test2', np.string_(['a','aaa']) , '[a.u.]' )
    mn.adddataset(grp, 'test3', 3.0, '[a.u.]' )
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