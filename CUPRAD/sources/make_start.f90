PROGRAM make_start
  USE HDF5
  USE write_listing
  USE default_inputs
  USE HDF5_helper_serial

  IMPLICIT NONE
  integer :: st
  character(100) :: dumstring, dumstring2, dumstring3
  logical :: testingmode=.FALSE., exists_h5_refractive_index, dumlog
  integer :: test_number = 1
  real(8) :: Energy
  real(8) :: dumreal

  integer, parameter :: available_dispersions(6) = (/1, 3, 4, 6, 7, 9/), available_ionisations(4) = (/1, 2, 3, 8/), &
                        available_beams(2) = (/1, 2/)


  PRINT *,"Pre-processor started"
  PRINT*, 'Specify name of parameterfile (HDF5 format)' 
  READ(5,*) filename


  ! Open FORTRAN HDF5 interface
  CALL h5open_f(error)

  ! IF ( (filename == "test") .or. (filename == "test2")) THEN ! this is an option for developers to run the code with pre-defaults inputs; if driving file exists, it's superior
  IF ( ANY( filename == available_tests ) ) THEN ! this is an option for developers to run the code with pre-defaults inputs; if driving file exists, it's superior
    test_number = get_test_number(filename)
    filename = "results.h5"
    INQUIRE(FILE=filename, EXIST=dumlog)
    
    IF (.NOT.dumlog) THEN
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)   
    ELSE
      CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
    ENDIF

    CALL h5lexists_f(file_id, CUPRAD_grp, dumlog, error)
    IF (.NOT.dumlog) THEN
      CALL h5gcreate_f(file_id, CUPRAD_grp, group_id, error)
      CALL h5gclose_f(group_id, error)   
    ENDIF

    CALL h5lexists_f(file_id, global_inps_grp, dumlog, error)
    IF (.NOT.dumlog) THEN
      CALL h5gcreate_f(file_id, global_inps_grp, group_id, error)
      CALL h5gclose_f(group_id, error)   
    ENDIF

    CALL h5lexists_f(file_id, in_grpname, dumlog, error)    
    IF (.NOT.dumlog) THEN
      CALL h5gcreate_f(file_id, in_grpname, group_id, error)
      CALL h5gclose_f(group_id, error)   
    ENDIF

    CALL h5fclose_f(file_id, error)
    testingmode = .TRUE.
  ENDIF


  OPEN(11,FILE='msg.tmp')
  WRITE(11,'(a)') TRIM(filename)
  CLOSE(11)
  
  ! Open the file
  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
  CALL h5gcreate_f(file_id, in_grpname//'/calculated', group_id, error)


  !-------------!
  ! Testingmode !
  !-------------!

  ! Testing mode presets all parameters and the code may run with an empty input.
  ! If some inputs are preent in the input file, they are superior.
  IF (testingmode) THEN
    CALL create_dset(file_id, in_grpname//'/calculated/test_test_number', test_number)
    CALL testing_values(test_number)
  ENDIF

  !---------------------------!
  ! Preset default parameters !
  !---------------------------!

  ! All the parameters of the medium (dispersion, ionization, Kerr, absrption, ...) are defined if this directive is present
  CALL h5lexists_f(file_id, global_inps_grp//'/gas_preset', dumlog, error)
  IF (dumlog) THEN
    CALL read_dset(file_id, global_inps_grp//'/gas_preset', gas_preset)
    CALL h5lexists_f(file_id, in_grpname//'/ionization_model', dumlog, error)    
    IF (dumlog) THEN
      CALL read_dset(file_id, in_grpname//'/ionization_model', ionization_model)
    ELSE
      print *, 'WARNING: gas is specified, but not the ionization model, PPT used'
      ionization_model = 'PPT'
    ENDIF
    gas_preset = TRIM(gas_preset)//'_'//TRIM(ionization_model)  
    print *, 'using preset gas model: ',  gas_preset
    CALL preset_parameters_gas
  ENDIF

  ! some parameters for a standard operation
  CALL preset_numerics
  CALL preset_laser


  !----------------------!
  ! Numerical properties !
  !----------------------!

  ! CALL save_or_replace(file_id, in_grpname//'/numerics_number_of_processors', num_proc, error, units_in = '[-]')
  num_proc = 4 ! kept from obsolete design, where the I/O handling depended on num_proc, no real use now
  CALL save_or_replace(file_id, in_grpname//'/numerics_run_time_in_hours', time_limit, error, units_in = '[h]')
  CALL save_or_replace(file_id, in_grpname//'/numerics_length_of_window_for_t_normalized_to_pulse_duration', lt, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_number_of_points_in_t', dim_t, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_length_of_window_for_r_normalized_to_beamwaist', lr, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_number_of_points_in_r', dim_r, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_number_of_absorber_points_in_time', absorb, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_phase_threshold_for_decreasing_delta_z', decrease, error)
  
  CALL save_or_replace(file_id, in_grpname//'/numerics_output_distance_in_z-steps_for_fluence_and_power', rhodist, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_radius_for_diagnostics', rfil_mm_phys, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_physical_first_stepwidth', delta_z_mm_phys, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_operators_t_t-1', switch_T, error)
  
  IF (switch_T.GT.4) THEN
      write(6,*) 'You have selected a bad value for the type of equation'
      write(6,*) ' You have to choose between 1, 2 or 3'
      write(6,*) ' The code will be stopped'
      STOP
  ENDIF

  CALL save_or_replace(file_id, in_grpname//'/numerics_type_of_input_beam', switch_start, error)

  !IF (NOT(ANY(available_beams == switch_start))) THEN
  IF (.NOT.ANY(available_beams == switch_start)) THEN
    write(6,*) 'You have selected a bad value for the type of input beamshape'
    write(6,*) ' You have to choose between 1 and 2'
    write(6,*) ' Continuation not implemented yet'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF   

  IF (2 == switch_start) THEN
    CALL save_or_replace(file_id, in_grpname//'/numerics_filename_for_method_2', inputfilename_t, error) 
  ELSEIF (3 == switch_start) THEN
    ! CALL save_or_replace(file_id, 'inputs/filename_for_method_3', inputfilename_c, error)
    CALL save_or_replace(file_id, in_grpname//'/numerics_amplituderatio_for_method_3', restartamp, error)
  ENDIF  

  CALL save_or_replace(file_id, in_grpname//'/numerics_spatial_noise_on_the_input_shape', noise_s, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_temporal_noise_on_the_input_shape', noise_t, error)
  CALL save_or_replace(file_id, in_grpname//'/numerics_noise_on_the_input_shape', noise, error)


  CALL save_or_replace(file_id, in_grpname//'/numerics_physical_output_distance_for_plasma_and_Efield', outlength_m_phys, error)

  CALL h5lexists_f(file_id, in_grpname//'/numerics_physical_output_distance_for_Efield_only', out_Efield, error)
  IF (out_Efield)  CALL save_or_replace(file_id, in_grpname//'/numerics_physical_output_distance_for_Efield_only', &
                                        outlength_Efield_m_phys, error)


  !--------!
  ! Medium !
  !--------!

  CALL save_or_replace(file_id, in_grpname//'/medium_physical_distance_of_propagation', proplength_m_phys, error)
  CALL save_or_replace(file_id, global_inps_grp//'/medium_pressure_in_bar', pressure, error)
  CALL save_or_replace(file_id, in_grpname//'/medium_pressure_in_bar', pressure, error)      ! Thisd ensures priority of the local quantity
  CALL save_or_replace(file_id, in_grpname//'/medium_effective_atmospheric_density_of_neutral_molecules', rhont_cm3_phys, &
                       error, units_in = '[1/cm3]')
  !CALL read_dset(file_id, 'inputs/effective_density_of_neutral_molecules', rhont_cm3_phys)
  !rhont_cm3_phys = 0.5

  CALL save_or_replace(group_id, 'medium_effective_density_of_neutral_molecules', pressure*rhont_cm3_phys, error, &
                       units_in = '[1/cm3]')


  !--------------------------------!
  ! Dispersion + Kerr + ionisation !
  !--------------------------------!

  CALL save_or_replace(file_id, in_grpname//'/dispersion_type_of_dispersion_law', switch_dispersion, error)
  
  IF (.NOT.ANY(available_dispersions == switch_dispersion)) THEN
    write(6,*) 'You have selected a bad value for the dispersion law'
    !write(6,*) ' You have to choose in integer between 1 or 7'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  SELECT CASE(switch_dispersion)
  CASE(1) ! Taylor
    CALL save_or_replace(file_id, in_grpname//'/dispersion_linear_refractive_index', n0, error)
    CALL save_or_replace(file_id, in_grpname//'/dispersion_inverse_gv_coefficient-n0_c', delta_k_p_fs_per_cm_phys, error)
    CALL save_or_replace(file_id, in_grpname//'/dispersion_gvd_coefficient', k_pp_fs2_per_cm_phys, error)
    CALL save_or_replace(file_id, in_grpname//'/dispersion_third_order_dispersion_coefficient', k_ppp_fs3_per_cm_phys, error)
    CALL save_or_replace(file_id, in_grpname//'/dispersion_fourth_order_dispersion_coefficient', k_pppp_fs4_per_cm_phys, error)
    CALL save_or_replace(file_id, in_grpname//'/dispersion_fifth_order_dispersion_coefficient', k_ppppp_fs5_per_cm_phys, error)
  CASE(4) ! from file
    CALL save_or_replace(file_id, in_grpname//'/dispersion_filename_for_method_4', dispfilename, error)
  END SELECT

  ! Kerr
  CALL save_or_replace(file_id, in_grpname//'/Kerr_ionised_atoms_relative_Kerr_response', ions_Kerr_ratio, error) 

  CALL save_or_replace(file_id, in_grpname//'/Kerr_nonlinear_refractive_index_kerr_coefficient', n2_phys, error)
  CALL save_or_replace(file_id, in_grpname//'/Kerr_chi5_coefficient', n4_phys, error)

  ! delayed Kerr
  CALL save_or_replace(file_id, in_grpname//'/Kerr_type_of_delayed_kerr_response', switch_dKerr, error)

  IF (switch_dKerr.GT.3) THEN
    write(6,*) 'You have selected a bad value for the type of Delayed Kerr response'
    write(6,*) ' You have to choose between 1, 2 or 3'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  IF (ANY( (/2, 3/) ==  switch_dKerr )) THEN
    CALL save_or_replace(file_id, in_grpname//'/Kerr_ratio_of_delayed_kerr_xdk', xdk, error)
    CALL save_or_replace(file_id, in_grpname//'/Kerr_time_of_delayed_kerr_tdk', tdk_fs_phys, error)
    IF (3==switch_dKerr) THEN
      CALL save_or_replace(file_id, in_grpname//'/Kerr_frequency_in_delayed_kerr_wr', raman_phys, error)
    ENDIF
  ENDIF


 CALL save_or_replace(file_id, in_grpname//'/ionization_type_of_ionization_method', switch_rho, error, units_in = '[-]')

  IF (.NOT.ANY(available_ionisations == switch_rho)) THEN
    write(6,*) 'You have selected a bad value for the type ionization method'
    write(6,*) ' You have to choose in integer between 1 and 8'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  IF (ANY( (/1, 2/) ==  switch_rho)) THEN
    CALL save_or_replace(file_id, in_grpname//'/ionization_mpi_cross_section_for_method_1-2', sigmak_phys, error, &
                         units_in = '[s-1cm2K/WK]')
  ELSEIF (3 == switch_rho) THEN
    CALL save_or_replace(file_id, in_grpname//'/ionization_ionization_potential_of_neutral_molecules', Ui_eV_phys, &
                         error, units_in = '[eV]')
    CALL save_or_replace(file_id, in_grpname//'/ionization_angular_momentum_for_method_3_7', angular_momentum, error, &
                         units_in = '[-]')
    CALL save_or_replace(file_id, in_grpname//'/ionization_effective_residue_charge_for_method_3-4_7', residue_charge, &
                         error, units_in = '[-]')
  ENDIF

  
  CALL save_or_replace(file_id, in_grpname//'/plasma_initial_electron_density', rho0_phys, error, units_in = '[1/cm3]')
  CALL save_or_replace(file_id, in_grpname//'/plasma_electron_colision_time', tauc_fs_phys, error, units_in = '[fs]')
  CALL save_or_replace(file_id, in_grpname//'/plasma_linear_recombination_coefficient', alpha_fs_phys, error, units_in = '[fs-1]')
  CALL save_or_replace(file_id, in_grpname//'/plasma_quadratic_recombination_(gasses)', alphaquad_fscm3_phys, error, &
                       units_in = '[fs-1cm3]')
  CALL save_or_replace(file_id, in_grpname//'/plasma_number_of_photons_involved_in_the_n-absorption', NN, error, units_in = '[-]')
  CALL save_or_replace(file_id, in_grpname//'/plasma_the_n-photon_absoption_cross_section', sigman_phys, error, &
                       units_in = '[s-1cm2N/Wn]')
  CALL save_or_replace(file_id, in_grpname//'/plasma_density_of_absorbing_molecules', rhoabs_cm3_phys, error, units_in = '[1/cm3]?')


  !-------!
  ! Laser !
  !-------!  
  
  CALL save_or_replace(file_id, in_grpname//'/laser_wavelength', lambda0_cm_phys, error)

  CALL save_or_replace(file_id, in_grpname//'/laser_degree_of_supergaussian', super_N, error)
  IF (super_N .NE. 1) THEN
    write(6,*) 'Warning, chosen beam shape differs from the Gaussian.'
    write(6,*) 'The Gaussian reference may not be representative.'
  ENDIF

  CALL save_or_replace(file_id, in_grpname//'/laser_degree_of_supergaussian_in_time', super_t, error)
  IF (super_N .NE. 1) THEN
    write(6,*) 'Warning, chosen temporal profile differs from the Gaussian.'
    write(6,*) 'The Gaussian reference may not be representative.'
  ENDIF

  CALL h5lexists_f(file_id, in_grpname//'/laser_pulse_duration_in_FWHM_Efield', dumlog, error)
  IF (dumlog) THEN
    dumstring = 'FWHM'; dumstring2 = 'Efield'; dumstring3 = in_grpname//'/laser_pulse_duration_in_FWHM_Efield'
  ENDIF
  CALL h5lexists_f(file_id, in_grpname//'/laser_pulse_duration_in_FWHM_Intensity', dumlog, error)
  IF (dumlog) THEN
    dumstring = 'FWHM'; dumstring2 = 'Intensity'; dumstring3 = in_grpname//'/laser_pulse_duration_in_FWHM_Intensity'
  ENDIF
  CALL h5lexists_f(file_id, in_grpname//'/laser_pulse_duration_in_rms_Efield', dumlog, error)
  IF (dumlog) THEN
    dumstring = 'rms'; dumstring2 = 'Efield'; dumstring3 = in_grpname//'/laser_pulse_duration_in_rms_Efield'
  ENDIF
  CALL h5lexists_f(file_id, in_grpname//'/laser_pulse_duration_in_rms_Intensity', dumlog, error)
  IF (dumlog) THEN
    dumstring = 'rms'; dumstring2 = 'Intensity'; dumstring3 = in_grpname//'/laser_pulse_duration_in_rms_Intensity'
  ENDIF
  CALL h5lexists_f(file_id, in_grpname//'/laser_pulse_duration_in_1_e_Intensity', dumlog, error)
  IF (dumlog) THEN
    dumstring = '1/e'; dumstring2 = 'Intensity'; dumstring3 = in_grpname//'/laser_pulse_duration_in_1_e_Intensity'
  ENDIF
  CALL h5lexists_f(file_id, in_grpname//'/laser_pulse_duration_in_1_e_Efield', dumlog, error)
  IF (dumlog) THEN
    dumstring = '1/e'; dumstring2 = 'Efield'; dumstring3 = in_grpname//'/laser_pulse_duration_in_1_e_Efield'
  ENDIF

  ! CALL h5lexists_f(file_id, 'inputs/laser_pulse_duration_in_rms', dumlog, error)
  ! IF (dumlog) dumstring = 'rms'
  ! CALL h5lexists_f(file_id, 'inputs/laser_pulse_duration_in_1_e', dumlog, error)
  ! IF (dumlog) dumstring = '1/e'

  CALL read_dset(file_id, dumstring3, tp_fs_phys)
  SELECT CASE(dumstring)
  CASE('FWHM')    
    !tp_fs_phys = FWHM2e_inv(tp_fs_phys)
    tp_fs_phys = Convert_pulse_duration(tp_fs_phys, dumstring, '1/e', type2_in = dumstring2, type2_out = 'Efield')
    CALL save_or_replace(group_id, 'laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')
  CASE('rms')
    ! CALL read_dset(file_id, 'inputs/laser_pulse_duration_in_rms', tp_fs_phys)
    !tp_fs_phys = FWHM2e_inv(tp_fs_phys)
    tp_fs_phys = Convert_pulse_duration(tp_fs_phys, dumstring, '1/e', type2_in = dumstring2, type2_out = 'Efield')
    CALL save_or_replace(group_id, 'laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')
  CASE('1/e')
    IF (dumstring2 == 'Intensity') THEN
      tp_fs_phys = Convert_pulse_duration(tp_fs_phys, dumstring, '1/e', type2_in = dumstring2, type2_out = 'Efield')
      CALL save_or_replace(group_id, 'laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')
    ENDIF
  CASE DEFAULT
    STOP "wrong pulse duration specification"
  END SELECT

  IF (debug_print) THEN
    print *, 'Efield'
    dumreal = Convert_pulse_duration(tp_fs_phys, '1/e', '1/e', type2_in = 'Efield', type2_out = 'Efield')
    print *, dumreal
    dumreal = Convert_pulse_duration(tp_fs_phys, '1/e', 'FWHM', type2_in = 'Efield', type2_out = 'Efield')
    print *, dumreal
    dumreal = Convert_pulse_duration(tp_fs_phys, '1/e', 'rms', type2_in = 'Efield', type2_out = 'Efield')
    print *, dumreal

    print *, 'Intensity'
    dumreal = Convert_pulse_duration(tp_fs_phys, '1/e', '1/e', type2_in = 'Efield', type2_out = 'Intensity')
    print *, dumreal
    dumreal = Convert_pulse_duration(tp_fs_phys, '1/e', 'FWHM', type2_in = 'Efield', type2_out = 'Intensity')
    print *, dumreal
    dumreal = Convert_pulse_duration(tp_fs_phys, '1/e', 'rms', type2_in = 'Efield', type2_out = 'Intensity')
    print *, dumreal
    !stop
  ENDIF


  CALL save_or_replace(group_id, 'laser_pulse_duration_in_1_e_Efield', Convert_pulse_duration(tp_fs_phys, '1/e', '1/e', &
                      type2_in = 'Efield', type2_out = 'Efield'), error, units_in = '[fs]')
  CALL save_or_replace(group_id, 'laser_pulse_duration_in_rms_Efield', Convert_pulse_duration(tp_fs_phys, '1/e', 'rms', &
                      type2_in = 'Efield', type2_out = 'Efield'), error, units_in = '[fs]')
  CALL save_or_replace(group_id, 'laser_pulse_duration_in_FWHM_Efield', Convert_pulse_duration(tp_fs_phys, '1/e', 'FWHM', &
                      type2_in = 'Efield', type2_out = 'Efield'), error, units_in = '[fs]')

  CALL save_or_replace(group_id, 'laser_pulse_duration_in_1_e_Intensity', Convert_pulse_duration(tp_fs_phys, '1/e', '1/e', &
                      type2_in = 'Efield', type2_out = 'Intensity'), error, units_in = '[fs]')

  CALL save_or_replace(group_id, 'laser_pulse_duration_in_rms_Intensity', Convert_pulse_duration(tp_fs_phys, '1/e', 'rms', &
                       type2_in = 'Efield', type2_out = 'Intensity'), error, units_in = '[fs]')
  CALL save_or_replace(group_id, 'laser_pulse_duration_in_FWHM_Intensity', Convert_pulse_duration(tp_fs_phys, '1/e', 'FWHM', &
                       type2_in = 'Efield', type2_out = 'Intensity'), error, units_in = '[fs]')


  
  CALL save_or_replace(file_id, in_grpname//'/laser_initial_chirp_phase', chirp_factor, error)
  
  CALL h5lexists_f(file_id, in_grpname//'/laser_beamwaist_entry', dumlog, error)
  IF (dumlog) THEN
    CALL save_or_replace(file_id, in_grpname//'/laser_focal_length_in_the_medium_cm', f_cm_phys, error) !! input given at the entrance plane, it's equal to the radius of the curvature
                                                                                                  ! this dataset in HDF5 has wrong unit [0 for no lense]
    
    Curvature_radius_entry = f_cm_phys * 1.D-2                              !!!!!!!!!!!!!!! CHECK SIGN
    
    CALL save_or_replace(file_id, in_grpname//'/laser_beamwaist_entry', w0_m_phys, error)
    CALL h5lexists_f(file_id, in_grpname//'/laser_ratio_pin_pcr', dumlog, error)
    IF (dumlog) THEN ! input given in the P_in/P_cr
      CALL save_or_replace(file_id, in_grpname//'/laser_ratio_pin_pcr', numcrit, error)
      ! convert
      Intensity_entry = ratio_Pin_Pcr_entry2I_entry(numcrit,w0_m_phys,pressure*n2_phys*1.D-4,lambda0_cm_phys*1.D-2)
      CALL save_or_replace(group_id, 'laser_intensity_entry', Intensity_entry, error, units_in = '[SI]')
    ELSE ! input given in the entrance intensity
      CALL save_or_replace(file_id, in_grpname//'/laser_intensity_entry', Intensity_entry, error)
      ! convert

	    print *, 'in make_start'
	    print *, 'n2_phys', n2_phys
	    print *, 'pressure', pressure
	    print *, 'all', pressure*n2_phys*1.D-4

      numcrit = I_entry2ratio_Pin_Pcr_entry(Intensity_entry,w0_m_phys,pressure*n2_phys*1.D-4,lambda0_cm_phys*1.D-2)
      CALL save_or_replace(group_id, 'laser_ratio_pin_pcr', numcrit, error, units_in = '[-]')
    ENDIF

    ! convert to focus values 
    IF (Curvature_radius_entry == 0.D0) THEN ! cumbersome as the inverse of the radius is not used... For compatibility.
      invCurvature_radius_entry = 0.D0
    ELSE
      invCurvature_radius_entry = 1.D0 / Curvature_radius_entry
    ENDIF   
    CALL Gaussian_entry2Gaussian_focus(Intensity_entry,w0_m_phys,invCurvature_radius_entry,Intensity_focus, waist_focus, &
                                       focus_position, lambda0_cm_phys*1.D-2)
  
    ! Store the reference Gaussian beam
    CALL save_or_replace(group_id, 'laser_focus_beamwaist_Gaussian', waist_focus, error, units_in = '[SI]')
    CALL save_or_replace(group_id, 'laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')
    CALL save_or_replace(group_id, 'laser_focus_position_Gaussian', focus_position, error, units_in = '[SI]')
  
  ELSE

    ! Input given by the reference Gaussian beam
    CALL save_or_replace(file_id, in_grpname//'/laser_focus_beamwaist_Gaussian', waist_focus, error)
    CALL save_or_replace(file_id, in_grpname//'/laser_focus_intensity_Gaussian', Intensity_focus, error)
    CALL save_or_replace(file_id, in_grpname//'/laser_focus_position_Gaussian', focus_position, error) !!!!!!!!!!!!!!! CHECK SIGN

    ! Convert
    CALL Gaussian_focus2Gaussian_entry(Intensity_focus,waist_focus,focus_position,Intensity_entry,w0_m_phys,&
                                       Curvature_radius_entry,lambda0_cm_phys*1.D-2)
    IF (Curvature_radius_entry == 0.D0) THEN ! cumbersome as the inverse of the radius is not used... For compatibility.
      f_cm_phys = 0.D0
    ELSE
      Curvature_radius_entry = 1.D0 / Curvature_radius_entry
      f_cm_phys = Curvature_radius_entry*1.D2
    ENDIF
    ! w0_cm_phys = w0_cm_phys*1.D2
    

    numcrit = I_entry2ratio_Pin_Pcr_entry(Intensity_entry,w0_m_phys,pressure*n2_phys*1.D-4,lambda0_cm_phys*1.D-2) ! intensity -> ratio

    CALL save_or_replace(group_id, 'laser_intensity_entry', Intensity_entry, error, units_in = '[SI]')
    CALL save_or_replace(group_id, 'laser_ratio_pin_pcr', numcrit, error, units_in = '[-]')
    CALL save_or_replace(group_id, 'laser_beamwaist_entry', w0_m_phys, error, units_in = '[m]')
    CALL save_or_replace(group_id, 'laser_focal_length_in_the_medium_cm', f_cm_phys, error, units_in = '[cm]') !! input given at the entrance plane, it's equal to the radius of the curvature    

  ENDIF

  w0_cm_phys = 1.D2*w0_m_phys ! It is used in the normalisation module, kept to minimise testing therein

  ! test Energy
  Energy = ratio_Pin_Pcr_entry2Energy(numcrit,pressure*n2_phys*1.D-4,lambda0_cm_phys*1.D-2,tp_fs_phys*1.D-15)
  IF (debug_print) print *, 'Energy=', Energy
  CALL save_or_replace(group_id, 'laser_Energy', Energy, error, units_in = '[J]')

  IF (debug_print) THEN
    print *, '========= Testing waists ========='
    print *, 'waist_entry =', w0_m_phys, 'm'
    print *, 'waist_Gauss =', waist_focus, 'm'
    print *, 'Pin/Pcr =', numcrit
    print *, '=================================='
  ENDIF


  CALL h5gclose_f(group_id, error)  

  !------------------------------------!
  ! CALL read_dset(file_id, 'inputs/', )!
  ! CALL read_dset(file_id, 'inputs/', ) !
  !------------------------------------!
  ! The original program closed the file here, so I'm closing it too
  !------------------------------------!
  ! Close the file.
  ! CALL h5fclose_f(file_id, error)

  CALL h5lexists_f(file_id,'refractive_index_table',exists_h5_refractive_index,error)

  IF (exists_h5_refractive_index) THEN
    !CALL h5open_f(error)
    !CALL h5fopen_f ("calculated_tables.h5", H5F_ACC_RDWR_F, file_id, error)
    CALL ask_for_size_1D(file_id, "/refractive_index_table/rgrid", i_x_max)
    CALL ask_for_size_1D(file_id, "/refractive_index_table/zgrid", i_z_max)
    ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
    CALL read_dset(file_id, "/refractive_index_table/refractive_index_shift", Indice, i_x_max, i_z_max)
    CALL read_dset(file_id, "/refractive_index_table/rgrid", xx_mum, i_x_max)
    CALL read_dset(file_id, "/refractive_index_table/zgrid", zz_mum, i_z_max)
    ! convert to microns
    xx_mum = 1.D6*xx_mum; zz_mum = 1.D6*zz_mum
    !CALL h5fclose_f(file_id, error)
    !CALL h5close_f(error)
    PRINT*, 'refractive index loaded from the input hdf5 archive'
  ELSE
    ! PRINT*, 'Specify name of indexfile, or 0 to ignore'
    ! READ(5,*) indexfile
    indexfile='0'
    IF (indexfile.EQ.'0') THEN
      i_x_max=1
      i_z_max=1
      ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
      xx_mum(i_x_max)=0.D0
      zz_mum(i_z_max)=0.D0
      Indice(i_x_max,i_z_max)=1.D0
    ELSE
      OPEN(12,FILE=indexfile,STATUS='UNKNOWN',FORM='FORMATTED')
      PRINT*, indexfile
      READ(unit=12,fmt=*,iostat=st) i_x_max, i_z_max
      PRINT*, i_x_max, i_z_max
      ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
      READ(unit=12,fmt=*,iostat=st) (xx_mum(i_x),i_x=1,i_x_max)
      DO i_z = 1, i_z_max
        READ(unit=12,fmt=*,iostat=st) zz_mum(i_z)
        READ(unit=12,fmt=*,iostat=st) (Indice(i_x,i_z),i_x=1,i_x_max)
      ENDDO
      CLOSE(12)
    ENDIF
  ENDIF
	 
  IF (debug_print) THEN
    PRINT*, 'index_rgrid'
    PRINT*, xx_mum
    PRINT*, 'index_zgrid'
    PRINT*, zz_mum

    PRINT*, 'indexes:'
    PRINT*, Indice(:,1)
    !PRINT*, Indice(:,2)

    PRINT*, 'acdispersion'
  ENDIF

  CALL compute_dispersion(switch_dispersion)
  CALL compute_parameters
  CALL write_listingfile

  ! save reference values
  ! CALL h5gcreate_f(file_id, 'logs', group_id, error)
  ! CALL save_or_replace(group_id, 'pulse_duration_1_e_Efield', tp_fs_phys*1.d-15, error, units_in = '[s]')
  ! CALL save_or_replace(group_id, 'Critical_power',Pcr_phys, error, units_in = '[W]')
  ! CALL h5gclose_f(group_id, error)

  CALL h5gcreate_f(file_id, output_groupname, group_id, error)  
  ALLOCATE(e_full(dim_t,dim_r))
  IF (debug_print) PRINT*, 'numproc', num_proc
  DO p=0,num_proc-1
    IF (debug_print) PRINT*, 'calc1', p
    CALL calc_startingfield(switch_start,p) 
    IF (debug_print) PRINT*, 'write1'
    CALL write_startingfile(p)
  ENDDO
  CALL h5gclose_f(group_id, error)

  CALL h5fclose_f(file_id, error)
  ! Close FORTRAN HDF5 interface.
  CALL h5close_f(error)  
  PRINT *,"Pre-processor done"

END PROGRAM make_start
