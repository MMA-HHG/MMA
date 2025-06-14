# this is an administration module for the MMA model

CUPRAD_group = 'CUPRAD'
CTDSE_group = 'CTDSE'
Hankel_group = 'Hankel'
global_inputs_group = 'global_inputs'
global_inputs_pre_ionised_subgroup = 'pre_ionised'
paths={'CUPRAD'                     : CUPRAD_group,
       'CUPRAD_inputs'              : CUPRAD_group + '/inputs',
       'CUPRAD_outputs'             : CUPRAD_group + '/outputs',
       'CUPRAD_logs'                : CUPRAD_group + '/logs',
       'CUPRAD_pre-processed'       : CUPRAD_group + '/pre-processed',
       'CUPRAD_ionisation_model'    : CUPRAD_group + '/ionisation_model',
       
       'CTDSE'                      : CTDSE_group,
       'CTDSE_inputs'               : CTDSE_group + '/inputs',
       'CTDSE_outputs'              : CTDSE_group + '/outputs',
       
       'Hankel'                     : Hankel_group,
       'Hankel_inputs'              : Hankel_group + '/inputs',
       'Hankel_outputs'             : Hankel_group + '/outputs',

       'global_inputs'              : global_inputs_group,
       'global_inputs_pre_ionised'  : global_inputs_group +'/'+ global_inputs_pre_ionised_subgroup}


filenames={'ionisation_tables'      : 'ionisation_tables.h5'}



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


global_variable_type_lists ={
    'S' : ['gas_preset'],
    'R' : ['medium_pressure_in_bar'],
    'I' : [], 'R-array' : []}

CUPRAD_variable_type_lists ={
    'I' : ['numerics_number_of_points_in_r', 'numerics_number_of_points_in_t',
           'numerics_operators_t_t-1'],
    'R' : ['laser_wavelength', 'laser_pulse_duration_in_1_e_Efield',
           'laser_focus_intensity_Gaussian', 'laser_focus_beamwaist_Gaussian',
           'laser_focus_position_Gaussian', 'medium_physical_distance_of_propagation',
           'numerics_physical_first_stepwidth', 'numerics_phase_threshold_for_decreasing_delta_z',
           'numerics_length_of_window_for_r_normalized_to_beamwaist',
           'numerics_length_of_window_for_t_normalized_to_pulse_duration',
           'numerics_number_of_absorber_points_in_time', 'numerics_physical_output_distance_for_plasma_and_Efield',
           'numerics_output_distance_in_z-steps_for_fluence_and_power',
           'numerics_radius_for_diagnostics', 'numerics_run_time_in_hours'],
    'S' : ['ionization_model'],
    'R-array': []}

CTDSE_variable_type_lists ={
    'I' : ['Nr_max', 'kr_step', 'kz_step', 'Nx_max']+
          ['print_'+foo for foo in ('GS', 'Efield', 'F_Efield', 'Source_Term', 
                                    'F_Source_Term', 'GS_population',
                                    'integrated_population', 'x_expectation_value')],
    'R' : ['x_int', 'dx', 'dt', 'CV_criterion_of_GS'],
    'S' : [], 'R-array' : []}

Hankel_variable_type_lists ={
    'I' : ['Nr_FF', 'kr_step', 'ko_step', 'Nr_max', 'Nthreads', 'store_cumulative_result'],
    'R' : ['distance_FF', 'rmax_FF'],
    'S' : ['XUV_table_type_dispersion', 'XUV_table_type_absorption'],
    'R-array': ['Harmonic_range']}

