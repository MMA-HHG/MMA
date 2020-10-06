PROGRAM make_start
  USE HDF5
  USE write_listing
  USE HDF5_helper

  IMPLICIT NONE
  PRINT *,"Pre-processor started"
  PRINT*, 'Specify name of parameterfile (HDF5 format)' 
  READ(5,*) filename

  ! Open FORTRAN HDF5 interface
  CALL h5open_f(error)
  
  ! Open the file
  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
  !_______________________________!

  CALL read_dset(file_id, 'inputs/number_of_processors', num_proc)
  !num_proc = 2
  CALL read_dset(file_id, 'inputs/run_time_in_hours', time_limit)
  time_limit = 0.48d0
  CALL read_dset(file_id, 'inputs/length_of_window_for_t_normalized_to_pulse_duration', lt)
  CALL read_dset(file_id, 'inputs/number_of_points_in_t', dim_t)
  CALL read_dset(file_id, 'inputs/length_of_window_for_r_normalized_to_beamwaist', lr)
  CALL read_dset(file_id, 'inputs/number_of_points_in_r', dim_r)
  CALL read_dset(file_id, 'inputs/number_of_absorber_points_in_time', absorb)
  CALL read_dset(file_id, 'inputs/phase_threshold_for_decreasing_delta_z', decrease)
  CALL read_dset(file_id, 'inputs/physical_distance_of_propagation', proplength_m_phys)
  proplength_m_phys = 5.d-3
  CALL read_dset(file_id, 'inputs/physical_output_distance_for_matlab_files', outlength_m_phys)
  outlength_m_phys = 1.d-3
  CALL read_dset(file_id, 'inputs/output_distance_in_z-steps_for_fluence_and_power', rhodist)
  rhodist = 2
  CALL read_dset(file_id, 'inputs/radius_for_diagnostics', rfil_mm_phys)
  CALL read_dset(file_id, 'inputs/physical_first_stepwidth', delta_z_mm_phys)
  CALL read_dset(file_id, 'inputs/operators_t_t-1', switch_T)
  
  if(switch_T.GT.4) then
      write(6,*) 'You have selected a bad value for the type of equation'
      write(6,*) ' You have to choose between 1, 2 or 3'
      write(6,*) ' The code will be stopped'
      STOP
  ENDIF                                     
  
  CALL read_dset(file_id, 'inputs/laser_wavelength', lambda0_cm_phys)
  CALL read_dset(file_id, 'inputs/beamwaist', w0_cm_phys)
  CALL read_dset(file_id, 'inputs/degree_of_supergaussian', super_N)
  CALL read_dset(file_id, 'inputs/pulse_duration_in_1_e', tp_fs_phys)
  CALL read_dset(file_id, 'inputs/degree_of_supergaussian_in_time', super_t)
  CALL read_dset(file_id, 'inputs/ratio_pin_pcr', numcrit)
  numcrit = 1.d0
  CALL read_dset(file_id, 'inputs/input', switch_start)

  if(switch_start.GT.4) then
    write(6,*) 'You have selected a bad value for the type of input beamshape'
    write(6,*) ' You have to choose between 1 or 4'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  CALL read_dset(file_id, 'inputs/filename_for_method_2', inputfilename_t) 
  CALL read_dset(file_id, 'inputs/filename_for_method_3', inputfilename_c)
  CALL read_dset(file_id, 'inputs/amplituderatio_for_method_3', restartamp)
  CALL read_dset(file_id, 'inputs/spatial_noise_on_the_input_shape', noise_s)
  CALL read_dset(file_id, 'inputs/temporal_noise_on_the_input_shape', noise_t)
  CALL read_dset(file_id, 'inputs/noise_on_the_input_shape', noise)
  CALL read_dset(file_id, 'inputs/focal_length_in_the_medium_cm', f_cm_phys) ! this dataset in HDF5 has wrong unit [0 for no lense]
  CALL read_dset(file_id, 'inputs/initial_chirp_phase', chirp_factor)
  CALL read_dset(file_id, 'inputs/pressure_in_bar', pressure)
  CALL read_dset(file_id, 'inputs/type_of_dispersion_law', switch_dispersion)
  
  if(switch_dispersion.GT.8) then
    write(6,*) 'You have selected a bad value for the dispersion law'
    write(6,*) ' You have to choose in integer between 1 or 7'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  CALL read_dset(file_id, 'inputs/filename_for_method_4', dispfilename)
  CALL read_dset(file_id, 'inputs/linear_refractive_index', n0)
  CALL read_dset(file_id, 'inputs/inverse_gv_coefficient-n0_c', delta_k_p_fs_per_cm_phys)
  CALL read_dset(file_id, 'inputs/gvd_coefficient', k_pp_fs2_per_cm_phys)
  CALL read_dset(file_id, 'inputs/third_order_dispersion_coefficient', k_ppp_fs3_per_cm_phys)
  CALL read_dset(file_id, 'inputs/fourth_order_dispersion_coefficient', k_pppp_fs4_per_cm_phys)
  CALL read_dset(file_id, 'inputs/fifth_order_dispersion_coefficient', k_ppppp_fs5_per_cm_phys)
  CALL read_dset(file_id, 'inputs/nonlinear_refractive_index_kerr_coefficient', n2_phys)
  CALL read_dset(file_id, 'inputs/type_of_delayed_kerr_response', switch_dKerr)

  if(switch_dKerr.GT.3) then
    write(6,*) 'You have selected a bad value for the type of Delayed Kerr response'
    write(6,*) ' You have to choose between 1, 2 or 3'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  CALL read_dset(file_id, 'inputs/ratio_of_delayed_kerr_xdk', xdk)
  CALL read_dset(file_id, 'inputs/time_of_delayed_kerr_tdk', tdk_fs_phys)
  CALL read_dset(file_id, 'inputs/frequency_in_delayed_kerr_wr', raman_phys)
  CALL read_dset(file_id, 'inputs/chi5_coefficient', n4_phys)
  CALL read_dset(file_id, 'inputs/effective_density_of_neutral_molecules', rhont_cm3_phys)
  CALL read_dset(file_id, 'inputs/ionization_poential_of_neutral_molecules', Ui_eV_phys)
  CALL read_dset(file_id, 'inputs/initial_electron_density', rho0_phys)
  CALL read_dset(file_id, 'inputs/type_of_ionization_method', switch_rho)
  if(switch_rho.GT.8) then 
    write(6,*) 'You have selected a bad value for the type ionization method'
    write(6,*) ' You have to choose in integer between 1 and 8'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF
  switch_rho = 8
  CALL read_dset(file_id, 'inputs/mpi_cross_section_for_method_1-2', sigmak_phys)
  CALL read_dset(file_id, 'inputs/angular_momentum_for_method_3_7', angular_momentum)
  CALL read_dset(file_id, 'inputs/effective_residue_charge_for_method_3-4_7', residue_charge)
  CALL read_dset(file_id, 'inputs/reduced_mass_of_hole-electron_for_method_5', reduced_mass)
  CALL read_dset(file_id, 'inputs/number_of_photons_to_ionize_from_slg1', KKp)
  CALL read_dset(file_id, 'inputs/cross_section_to_ionize_from_slg1', sigmakp_phys)
  CALL read_dset(file_id, 'inputs/density_of_defects_slg1', rhoslg1_phys)
  CALL read_dset(file_id, 'inputs/cross_section_to_ionize_from_slg2', sigma_phys)
  CALL read_dset(file_id, 'inputs/sigmacvref_for_i_ref', sigmacv_ref_phys)
  CALL read_dset(file_id, 'inputs/reference_intensity_i_ref', I_ref_phys)
  CALL read_dset(file_id, 'inputs/reference_exponent_exp_ref', exp_ref)
  CALL read_dset(file_id, 'inputs/number_of_photons_to_populate_slg2', KKpp)
  CALL read_dset(file_id, 'inputs/cross_section_to_populate_slg2', sigmakpp_phys)
  CALL read_dset(file_id, 'inputs/saturation_density_for_slg2', rhosat_phys)
  CALL read_dset(file_id, 'inputs/effective_density_of_neutral_n2_molecules', rhont_N2_cm3_phys)
  CALL read_dset(file_id, 'inputs/ionization_poential_of_neutral_n2_molecules', Ui_N2_eV_phys)
  CALL read_dset(file_id, 'inputs/angular_momentum_n2_for_method_7', angular_momentum_N2)
  CALL read_dset(file_id, 'inputs/effective_residue_charge_n2_for_method_7', residue_charge_N2)
  CALL read_dset(file_id, 'inputs/initial_free_electron_temperature', T_init_eV_phys)
  CALL read_dset(file_id, 'inputs/electron_colision_time', tauc_fs_phys)
  CALL read_dset(file_id, 'inputs/linear_recombination_coefficient', alpha_fs_phys)
  CALL read_dset(file_id, 'inputs/linear_recombination_coefficient_(6:_slg1_elec.)', alpha1_fs_phys)
  CALL read_dset(file_id, 'inputs/linear_recombination_coefficient_(6:_holes)', alphah_fs_phys)
  CALL read_dset(file_id, 'inputs/quadratic_recombination_(gasses)', alphaquad_fscm3_phys)
  CALL read_dset(file_id, 'inputs/number_of_photons_involved_in_the_n-absorption', NN)
  CALL read_dset(file_id, 'inputs/the_n-photon_absoption_cross_section', sigman_phys)
  CALL read_dset(file_id, 'inputs/density_of_absorbing_molecules', rhoabs_cm3_phys)

  !------------------------------------!
  ! CALL read_dset(file_id, 'inputs/', )!
  ! CALL read_dset(file_id, 'inputs/', ) !
  !------------------------------------!
  ! The original program closed the file here, so I'm closing it too
  !------------------------------------!
  ! Close the file.
  ! CALL h5fclose_f(file_id, error)
  
  PRINT*, 'Specify name of indexfile, or 1 to load from given HDF5 file, or 0 to ignore'
  READ(5,*) indexfile

  IF (indexfile.EQ.'0') THEN
    i_x_max=1
    i_z_max=1
    ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
    xx_mum(i_x_max)=0.D0
    zz_mum(i_z_max)=0.D0
    Indice(i_x_max,i_z_max)=1.D0
  ELSE IF (indexfile.EQ.'1') THEN
    CALL h5open_f(error)
    CALL h5fopen_f ("calculated_tables.h5", H5F_ACC_RDWR_F, file_id, error)
    CALL ask_for_size_1D(file_id, "/indexes/r_vector", i_x_max)
    CALL ask_for_size_1D(file_id, "/indexes/z_vector", i_z_max)
    ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
    CALL read_dset(file_id, "/indexes/indexes", Indice, i_x_max, i_z_max)
    CALL read_dset(file_id, "/indexes/r_vector", xx_mum, i_x_max)
    CALL read_dset(file_id, "/indexes/z_vector", zz_mum, i_z_max)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
  ELSE
    OPEN(12,FILE=indexfile,STATUS='UNKNOWN',FORM='UNFORMATTED')
    READ(12) i_x_max, i_z_max
    ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
    READ(12) (xx_mum(i_x),i_x=1,i_x_max)
    DO i_z = 1, i_z_max
      READ(12) zz_mum(i_z)
      READ(12) (Indice(i_x,i_z),i_x=1,i_x_max)
    ENDDO
    CLOSE(12)
  ENDIF

  PRINT*, 'acdispersion'

  CALL compute_dispersion(switch_dispersion)
  CALL compute_parameters
  CALL write_listingfile
  CALL h5gcreate_f(file_id, output_groupname, group_id, error)  
  ALLOCATE(e_full(dim_t,dim_r))
  PRINT*, 'numproc', num_proc
  DO p=0,num_proc-1
    PRINT*, 'calc1', p
    CALL calc_startingfield(switch_start,p) 
    PRINT*, 'write1'
    CALL write_startingfile(p)
  ENDDO
  CALL h5gclose_f(group_id, error)

  CALL h5fclose_f(file_id, error)
  ! Close FORTRAN HDF5 interface.
  CALL h5close_f(error)  
  PRINT *,"Pre-processor done"

END PROGRAM make_start
