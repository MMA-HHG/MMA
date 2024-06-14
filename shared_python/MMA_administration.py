# this is an administration module for the MMA model

CUPRAD_group = 'CUPRAD'
CTDSE_group = 'CTDSE'
Hankel_group = 'Hankel'
global_inputs_group = 'global_inputs'
paths={'CUPRAD'             : CUPRAD_group,
       'CUPRAD_inputs'      : CUPRAD_group + '/inputs',
       'CUPRAD_outputs'     : CUPRAD_group + '/outputs',
       'CUPRAD_logs'        : CUPRAD_group + '/logs',
       
       'CTDSE '             : CTDSE_group,
       'CTDSE_inputs'       : CTDSE_group + '/inputs',
       'CTDSE_outputs'      : CTDSE_group + '/outputs',

       'global_inputs'      : global_inputs_group}





def pressure_constructor(h5_handle):
    
    # obtain the base unmodulated pressure
    pressure_base = h5_handle[paths['global_inputs']+'/medium_pressure_in_bar'][()]
    # check whether the density modulation is applied
    if ('density_mod' in h5_handle[paths['global_inputs']].keys()):
        dens_mod_grp = h5_handle[paths['global_inputs']+'/'+'density_mod']
        pressure = {'value' : pressure_base*dens_mod_grp['table'][:]}
        if ('zgrid' in dens_mod_grp.keys()):
            pressure['zgrid'] = dens_mod_grp['zgrid'][:]
        if ('rgrid' in dens_mod_grp.keys()):
            pressure['rgrid'] = dens_mod_grp['rgrid'][:]
    else:
        pressure = pressure_base
        
    return pressure
