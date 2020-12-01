PROGRAM make_start
  USE HDF5
  USE write_listing
  USE default_inputs
  USE HDF5_helper

  IMPLICIT NONE
  integer :: st
  logical :: testingmode=.FALSE., exists_h5_refractive_index, dumlog
  real(8) :: Intensity_entry, Intensity_focus, waist_focus, Curvature_radius_entry, focus_position


  PRINT *,"Pre-processor started"
  PRINT*, 'Specify name of parameterfile (HDF5 format)' 
  READ(5,*) filename


  ! Open FORTRAN HDF5 interface
  CALL h5open_f(error)

  IF (filename == "test") THEN ! this is an option for developers to run the code with pre-defaults inputs; if driving file exists, it's superior
    filename = "results.h5"
    INQUIRE(FILE=filename, EXIST=dumlog)
    IF (NOT(dumlog)) THEN
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)   
    ELSE
      CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
    ENDIF

    CALL h5lexists_f(file_id, 'inputs', dumlog, error)
    IF (NOT(dumlog)) THEN
      CALL h5gcreate_f(file_id, 'inputs', group_id, error)
      CALL h5gclose_f(group_id, error)   
    ENDIF

    CALL h5fclose_f(file_id, error)
    testingmode = .TRUE.
  ENDIF


  
  ! Open the file
  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

  !_______________________________!
  ! CALL preset_numerics
  IF (testingmode) THEN
    CALL testing_values
  ENDIF

  CALL preset_laser


  CALL  h5lexists_f(file_id, 'inputs/gas_preset', dumlog, error)
  IF (dumlog) THEN
    CALL read_dset(file_id, 'inputs/gas_preset', gas_preset)
    CALL preset_parameters_gas
  ENDIF


  CALL save_or_replace(file_id, 'inputs/number_of_processors', num_proc, error, units_in = '[-]')
  !num_proc = 2
  CALL save_or_replace(file_id, 'inputs/run_time_in_hours', time_limit, error, units_in = '[h]')
  !time_limit = 0.48d0
  CALL save_or_replace(file_id, 'inputs/length_of_window_for_t_normalized_to_pulse_duration', lt, error)
  CALL save_or_replace(file_id, 'inputs/number_of_points_in_t', dim_t, error)
  CALL save_or_replace(file_id, 'inputs/length_of_window_for_r_normalized_to_beamwaist', lr, error)
  CALL save_or_replace(file_id, 'inputs/number_of_points_in_r', dim_r, error)
  CALL save_or_replace(file_id, 'inputs/number_of_absorber_points_in_time', absorb, error)
  CALL save_or_replace(file_id, 'inputs/phase_threshold_for_decreasing_delta_z', decrease, error)
  CALL save_or_replace(file_id, 'inputs/physical_distance_of_propagation', proplength_m_phys, error)
  !proplength_m_phys = 5.d-3
  CALL save_or_replace(file_id, 'inputs/physical_output_distance_for_matlab_files', outlength_m_phys, error)
  !outlength_m_phys = 1.d-3
  CALL save_or_replace(file_id, 'inputs/output_distance_in_z-steps_for_fluence_and_power', rhodist, error)
  !rhodist = 2
  CALL save_or_replace(file_id, 'inputs/radius_for_diagnostics', rfil_mm_phys, error)
  CALL save_or_replace(file_id, 'inputs/physical_first_stepwidth', delta_z_mm_phys, error)
  CALL save_or_replace(file_id, 'inputs/operators_t_t-1', switch_T, error)
  
  if(switch_T.GT.4) then
      write(6,*) 'You have selected a bad value for the type of equation'
      write(6,*) ' You have to choose between 1, 2 or 3'
      write(6,*) ' The code will be stopped'
      STOP
  ENDIF                                     
  
  CALL save_or_replace(file_id, 'inputs/laser_wavelength', lambda0_cm_phys, error)
  ! CALL save_or_replace(file_id, 'inputs/beamwaist', w0_cm_phys, error)
  ! This will be recomputed once n2 is known
  CALL save_or_replace(file_id, 'inputs/degree_of_supergaussian', super_N, error)
  CALL h5lexists_f(file_id, 'inputs/pulse_duration_in_1_e', dumlog, error)
  IF (dumlog) THEN
    CALL read_dset(file_id, 'inputs/pulse_duration_in_1_e', tp_fs_phys)
    CALL save_or_replace(file_id, 'inputs/pulse_duration_in_FWHM', e_inv2FWHM(tp_fs_phys), error, units_in = '[fs]')
  ELSE
    CALL read_dset(file_id, 'inputs/pulse_duration_in_FWHM', tp_fs_phys)
    tp_fs_phys = FWHM2e_inv(tp_fs_phys)
    CALL save_or_replace(file_id, 'inputs/pulse_duration_in_1_e', tp_fs_phys, error, units_in = '[fs]')
  ENDIF
  CALL save_or_replace(file_id, 'inputs/degree_of_supergaussian_in_time', super_t, error)
  ! CALL save_or_replace(file_id, 'inputs/ratio_pin_pcr', numcrit, error)
  ! This will be recomputed once n2 is known
  !numcrit = 1.d0
  CALL save_or_replace(file_id, 'inputs/input', switch_start, error)

  if(switch_start.GT.4) then
    write(6,*) 'You have selected a bad value for the type of input beamshape'
    write(6,*) ' You have to choose between 1 or 4'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  CALL save_or_replace(file_id, 'inputs/filename_for_method_2', inputfilename_t, error) 
  CALL save_or_replace(file_id, 'inputs/filename_for_method_3', inputfilename_c, error)
  CALL save_or_replace(file_id, 'inputs/amplituderatio_for_method_3', restartamp, error)
  CALL save_or_replace(file_id, 'inputs/spatial_noise_on_the_input_shape', noise_s, error)
  CALL save_or_replace(file_id, 'inputs/temporal_noise_on_the_input_shape', noise_t, error)
  CALL save_or_replace(file_id, 'inputs/noise_on_the_input_shape', noise, error)
  !CALL save_or_replace(file_id, 'inputs/focal_length_in_the_medium_cm', f_cm_phys, error) ! this dataset in HDF5 has wrong unit [0 for no lense]
  CALL save_or_replace(file_id, 'inputs/initial_chirp_phase', chirp_factor, error)
  
  
  CALL save_or_replace(file_id, 'inputs/pressure_in_bar', pressure, error)
  !pressure = 0.5d0
  CALL save_or_replace(file_id, 'inputs/type_of_dispersion_law', switch_dispersion, error)
  
  if(switch_dispersion.GT.9) then
    write(6,*) 'You have selected a bad value for the dispersion law'
    write(6,*) ' You have to choose in integer between 1 or 7'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  CALL save_or_replace(file_id, 'inputs/filename_for_method_4', dispfilename, error)
  CALL save_or_replace(file_id, 'inputs/linear_refractive_index', n0, error)
  CALL save_or_replace(file_id, 'inputs/inverse_gv_coefficient-n0_c', delta_k_p_fs_per_cm_phys, error)
  CALL save_or_replace(file_id, 'inputs/gvd_coefficient', k_pp_fs2_per_cm_phys, error)
  CALL save_or_replace(file_id, 'inputs/third_order_dispersion_coefficient', k_ppp_fs3_per_cm_phys, error)
  CALL save_or_replace(file_id, 'inputs/fourth_order_dispersion_coefficient', k_pppp_fs4_per_cm_phys, error)
  CALL save_or_replace(file_id, 'inputs/fifth_order_dispersion_coefficient', k_ppppp_fs5_per_cm_phys, error)


  CALL save_or_replace(file_id, 'inputs/nonlinear_refractive_index_kerr_coefficient', n2_phys, error)
  CALL save_or_replace(file_id, 'inputs/type_of_delayed_kerr_response', switch_dKerr, error)

  if(switch_dKerr.GT.3) then
    write(6,*) 'You have selected a bad value for the type of Delayed Kerr response'
    write(6,*) ' You have to choose between 1, 2 or 3'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF

  CALL save_or_replace(file_id, 'inputs/ratio_of_delayed_kerr_xdk', xdk, error)
  CALL save_or_replace(file_id, 'inputs/time_of_delayed_kerr_tdk', tdk_fs_phys, error)
  CALL save_or_replace(file_id, 'inputs/frequency_in_delayed_kerr_wr', raman_phys, error)
  CALL save_or_replace(file_id, 'inputs/chi5_coefficient', n4_phys, error)

  ! here we know n2 needed for the Pin/Pcr
  CALL h5lexists_f(file_id, 'inputs/beamwaist', dumlog, error)
  IF (dumlog) THEN
    CALL save_or_replace(file_id, 'inputs/focal_length_in_the_medium_cm', f_cm_phys, error) !! input given at the entrance plane, it's equal to the radius of the curvature
    
    Curvature_radius_entry = f_cm_phys * 1.D-2                              !!!!!!!!!!!!!!! CHECK SIGN
    
    CALL save_or_replace(file_id, 'inputs/beamwaist', w0_cm_phys, error)
    CALL h5lexists_f(file_id, 'inputs/ratio_pin_pcr', dumlog, error)
    IF (dumlog) THEN ! input given in the P_in/P_cr
      CALL save_or_replace(file_id, 'inputs/ratio_pin_pcr', numcrit, error)
      Intensity_entry = ratio_Pin_Pcr_entry2I_entry(numcrit,w0_cm_phys*1.D-2,pressure*n2_phys,lambda0_cm_phys*1.D-2)
      CALL save_or_replace(file_id, 'inputs/intensity_entry', Intensity_entry, error, units_in = '[SI]')
    ELSE ! input given in the entrance intensity
      CALL save_or_replace(file_id, 'inputs/intensity_entry', Intensity_entry, error)
      ! convert
      numcrit = I_entry2ratio_Pin_Pcr_entry(Intensity_entry,w0_cm_phys*1.D-2,pressure*n2_phys,lambda0_cm_phys*1.D-2)
      CALL save_or_replace(file_id, 'inputs/ratio_pin_pcr', numcrit, error, units_in = '[-]')
    ENDIF

    ! convert to focus values    
    CALL Gaussian_entry2Gaussian_focus(Intensity_entry,w0_cm_phys*1.D-2,Curvature_radius_entry,Intensity_focus, waist_focus, focus_position, lambda0_cm_phys*1.D-2)
  ELSE

    ! Input given by the reference Gaussian beam
    CALL save_or_replace(file_id, 'inputs/focus_beamwaist_Gaussian', waist_focus, error)
    CALL save_or_replace(file_id, 'inputs/focus_intensity_Gaussian', Intensity_focus, error)
    CALL save_or_replace(file_id, 'inputs/focus_position_Gaussian', focus_position, error) !!!!!!!!!!!!!!! CHECK SIGN

    ! Convert
    CALL Gaussian_focus2Gaussian_entry(Intensity_focus,waist_focus,focus_position,Intensity_entry,w0_cm_phys,Curvature_radius_entry,lambda0_cm_phys*1.D-2)
    w0_cm_phys = w0_cm_phys*1.D2
    f_cm_phys = Curvature_radius_entry*1.D2

    numcrit = I_entry2ratio_Pin_Pcr_entry(Intensity_entry,w0_cm_phys*1.D-2,pressure*n2_phys,lambda0_cm_phys*1.D-2) ! intensity -> ratio

    CALL save_or_replace(file_id, 'inputs/intensity_entry', Intensity_entry, error, units_in = '[SI]')
    CALL save_or_replace(file_id, 'inputs/ratio_pin_pcr', numcrit, error, units_in = '[-]')
    CALL save_or_replace(file_id, 'inputs/beamwaist', w0_cm_phys, error, units_in = '[cm]')
    CALL save_or_replace(file_id, 'inputs/focal_length_in_the_medium_cm', f_cm_phys, error, units_in = '[cm]') !! input given at the entrance plane, it's equal to the radius of the curvature    
    
    !CALL save_or_replace(file_id, 'inputs/ratio_pin_pcr', numcrit, error)
  ENDIF


  CALL save_or_replace(file_id, 'inputs/effective_density_of_neutral_molecules', rhont_cm3_phys, error, units_in = '[1/cm3]')
  !CALL read_dset(file_id, 'inputs/effective_density_of_neutral_molecules', rhont_cm3_phys)
  !rhont_cm3_phys = 0.5
  CALL save_or_replace(file_id, 'inputs/ionization_poential_of_neutral_molecules', Ui_eV_phys, error, units_in = '[eV]')
  CALL save_or_replace(file_id, 'inputs/initial_electron_density', rho0_phys, error, units_in = '[1/cm3]')
  CALL save_or_replace(file_id, 'inputs/type_of_ionization_method', switch_rho, error, units_in = '[-]')
  if(switch_rho.GT.8) then 
    write(6,*) 'You have selected a bad value for the type ionization method'
    write(6,*) ' You have to choose in integer between 1 and 8'
    write(6,*) ' The code will be stopped'
    STOP
  ENDIF
  switch_rho = 3 !!! testing PPT

  CALL save_or_replace(file_id, 'inputs/mpi_cross_section_for_method_1-2', sigmak_phys, error, units_in = '[s-1cm2K/WK]')
  CALL save_or_replace(file_id, 'inputs/angular_momentum_for_method_3_7', angular_momentum, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/effective_residue_charge_for_method_3-4_7', residue_charge, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/reduced_mass_of_hole-electron_for_method_5', reduced_mass, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/number_of_photons_to_ionize_from_slg1', KKp, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/cross_section_to_ionize_from_slg1', sigmakp_phys, error, units_in = '[s-1cm2Kp/WKp]')
  CALL save_or_replace(file_id, 'inputs/density_of_defects_slg1', rhoslg1_phys, error, units_in = '[cm-3]')
  CALL save_or_replace(file_id, 'inputs/cross_section_to_ionize_from_slg2', sigma_phys, error, units_in = '[s-1cm2/W]')
  CALL save_or_replace(file_id, 'inputs/sigmacvref_for_i_ref', sigmacv_ref_phys, error, units_in = '[cm2/s]')
  CALL save_or_replace(file_id, 'inputs/reference_intensity_i_ref', I_ref_phys, error, units_in = '[W/cm^2]')
  CALL save_or_replace(file_id, 'inputs/reference_exponent_exp_ref', exp_ref, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/number_of_photons_to_populate_slg2', KKpp, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/cross_section_to_populate_slg2', sigmakpp_phys, error, units_in = '[s-1cm2Kpp/WKpp]')
  CALL save_or_replace(file_id, 'inputs/saturation_density_for_slg2', rhosat_phys, error, units_in = '[cm-3]')
  CALL save_or_replace(file_id, 'inputs/effective_density_of_neutral_n2_molecules', rhont_N2_cm3_phys, error, units_in = '[cm-3]')
  CALL save_or_replace(file_id, 'inputs/ionization_poential_of_neutral_n2_molecules', Ui_N2_eV_phys, error, units_in = '[eV]')
  CALL save_or_replace(file_id, 'inputs/angular_momentum_n2_for_method_7', angular_momentum_N2, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/effective_residue_charge_n2_for_method_7', residue_charge_N2, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/initial_free_electron_temperature', T_init_eV_phys, error, units_in = '[eV]')

  CALL save_or_replace(file_id, 'inputs/electron_colision_time', tauc_fs_phys, error, units_in = '[fs]')
  CALL save_or_replace(file_id, 'inputs/linear_recombination_coefficient', alpha_fs_phys, error, units_in = '[fs-1]')
  CALL save_or_replace(file_id, 'inputs/linear_recombination_coefficient_(6:_slg1_elec.)', alpha1_fs_phys, error, units_in = '[fs-1]')
  CALL save_or_replace(file_id, 'inputs/linear_recombination_coefficient_(6:_holes)', alphah_fs_phys, error, units_in = '[fs-1]')
  CALL save_or_replace(file_id, 'inputs/quadratic_recombination_(gasses)', alphaquad_fscm3_phys, error, units_in = '[fs-1cm3]')

  CALL save_or_replace(file_id, 'inputs/number_of_photons_involved_in_the_n-absorption', NN, error, units_in = '[-]')
  CALL save_or_replace(file_id, 'inputs/the_n-photon_absoption_cross_section', sigman_phys, error, units_in = '[s-1cm2N/Wn]')
  CALL save_or_replace(file_id, 'inputs/density_of_absorbing_molecules', rhoabs_cm3_phys, error, units_in = '[1/cm3]?')

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
    PRINT*, 'Specify name of indexfile, or 0 to ignore'
    READ(5,*) indexfile
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
	 
  PRINT*, 'index_rgrid'
  PRINT*, xx_mum
  PRINT*, 'index_zgrid'
  PRINT*, zz_mum

  PRINT*, 'indexes:'
  PRINT*, Indice(:,1)
  PRINT*, Indice(:,2)

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
