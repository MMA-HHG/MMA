PROGRAM make_start
  USE write_listing
  IMPLICIT NONE

  CHARACTER(LEN=15)      :: null_char,filename

  PRINT*, 'Specify name of parameterfile'
  READ(5,*) filename

  OPEN(UNIT=2, FILE=filename, FORM='FORMATTED')
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, i15)') num_proc
  READ(2, '(t57, e15.8)') time_limit
  READ(2, '(t57, e15.8)')lt
  READ(2, '(t57, i15)') dim_t
  READ(2, '(t57, e15.8)') lr
  READ(2, '(t57, i15)') dim_r
  READ(2, '(t57, i15)') absorb
  READ(2, '(t57, e15.8)') decrease
  READ(2, '(t57, e15.8)') proplength_m_phys
  READ(2, '(t57, e15.8)') outlength_m_phys
  READ(2, '(t57, i15)') rhodist
  READ(2, '(t57, e15.8)') rfil_mm_phys
  READ(2, '(t57, e15.8)') delta_z_mm_phys
  READ(2, '(t57, i15)') switch_T
  if(switch_T.GT.4) then
     write(6,*) 'You have selected a bad value for the type of equation'
     write(6,*) ' You have to choose between 1, 2 or 3'
     write(6,*) ' The code will be stopped'
     STOP
  ENDIF
  READ(2, '(t57, a15)') null_char 
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') lambda0_cm_phys
  READ(2, '(t57, e15.8)') w0_cm_phys
  READ(2, '(t57, i15)') super_N
  READ(2, '(t57, e15.8)') tp_fs_phys
  READ(2, '(t57, i15)') super_t
  READ(2, '(t57, e15.8)') numcrit
  READ(2, '(t57, i15)') switch_start
  if(switch_start.GT.4) then
     write(6,*) 'You have selected a bad value for the type of input beamshape'
     write(6,*) ' You have to choose between 1 or 4'
     write(6,*) ' The code will be stopped'
     STOP
  ENDIF
  READ(2, '(t57, a15)') null_char 
  READ(2, '(t57, a15)') inputfilename_t
  READ(2, '(t57, a15)') inputfilename_c
  READ(2, '(t57, e15.8)') restartamp
  READ(2, '(t57, e15.8)') noise_s
  READ(2, '(t57, e15.8)') noise_t
  READ(2, '(t57, e15.8)') noise
  READ(2, '(t57, e15.8)') f_cm_phys
  READ(2, '(t57, e15.8)') chirp_factor
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') pressure
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, i15)') switch_dispersion
  if(switch_dispersion.GT.8) then
     write(6,*) 'You have selected a bad value for the dispersion law'
     write(6,*) ' You have to choose in integer between 1 or 7'
     write(6,*) ' The code will be stopped'
     STOP
  ENDIF
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, a15)') dispfilename
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') n0
  READ(2, '(t57, e15.8)') delta_k_p_fs_per_cm_phys
  READ(2, '(t57, e15.8)') k_pp_fs2_per_cm_phys
  READ(2, '(t57, e15.8)') k_ppp_fs3_per_cm_phys
  READ(2, '(t57, e15.8)') k_pppp_fs4_per_cm_phys
  READ(2, '(t57, e15.8)') k_ppppp_fs5_per_cm_phys
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') n2_phys
  READ(2, '(t57, i15)') switch_dKerr
  if(switch_dKerr.GT.3) then
     write(6,*) 'You have selected a bad value for the type of Delayed Kerr response'
     write(6,*) ' You have to choose between 1, 2 or 3'
     write(6,*) ' The code will be stopped'
     STOP
  ENDIF
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') xdk
  READ(2, '(t57, e15.8)') tdk_fs_phys
  READ(2, '(t57, e15.8)') raman_phys
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') n4_phys
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') rhont_cm3_phys
  READ(2, '(t57, e15.8)') Ui_eV_phys
  READ(2, '(t57, e15.8)') rho0_phys
  READ(2, '(t57, i15)') switch_rho
  if(switch_rho.GT.8) then
     write(6,*) 'You have selected a bad value for the type ionization method'
     write(6,*) ' You have to choose in integer between 1 and 8'
     write(6,*) ' The code will be stopped'
     STOP
  ENDIF
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') sigmak_phys
  READ(2, '(t57, i15)') angular_momentum
  READ(2, '(t57, e15.8)') residue_charge
  READ(2, '(t57, e15.8)') reduced_mass
  READ(2, '(t57, i15)') KKp
  READ(2, '(t57, e15.8)') sigmakp_phys
  READ(2, '(t57, e15.8)') rhoslg1_phys
  READ(2, '(t57, e15.8)') sigma_phys
  READ(2, '(t57, e15.8)') sigmacv_ref_phys
  READ(2, '(t57, e15.8)') I_ref_phys
  READ(2, '(t57, e15.8)') exp_ref
  READ(2, '(t57, i15)') KKpp
  READ(2, '(t57, e15.8)') sigmakpp_phys
  READ(2, '(t57, e15.8)') rhosat_phys
  READ(2, '(t57, e15.8)') rhont_N2_cm3_phys
  READ(2, '(t57, e15.8)') Ui_N2_eV_phys
  READ(2, '(t57, i15)') angular_momentum_N2
  READ(2, '(t57, e15.8)') residue_charge_N2
  READ(2, '(t57, e15.8)') T_init_eV_phys
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, e15.8)') tauc_fs_phys
  READ(2, '(t57, e15.8)') alpha_fs_phys
  READ(2, '(t57, e15.8)') alpha1_fs_phys
  READ(2, '(t57, e15.8)') alphah_fs_phys
  READ(2, '(t57, e15.8)') alphaquad_fscm3_phys
  READ(2, '(t57, a15)') null_char
  READ(2, '(t57, i15)') NN
  READ(2, '(t57, e15.8)') sigman_phys
  READ(2, '(t57, e15.8)') rhoabs_cm3_phys

  CLOSE(2)

  PRINT*, 'Specify name of indexfile, or 0 to ignore'
  READ(5,*) indexfile

  IF (indexfile.NE.'0') THEN
     OPEN(12,FILE=indexfile,STATUS='UNKNOWN',FORM='UNFORMATTED')
     READ(12) i_x_max, i_z_max
     ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
     READ(12) (xx_mum(i_x),i_x=1,i_x_max)
     DO i_z = 1, i_z_max
        READ(12) zz_mum(i_z)
        READ(12) (Indice(i_x,i_z),i_x=1,i_x_max)
     ENDDO
     CLOSE(12)
  ELSE
     i_x_max=1
     i_z_max=1
     ALLOCATE(xx_mum(i_x_max),zz_mum(i_z_max),Indice(i_x_max,i_z_max))
     xx_mum(i_x_max)=0.D0
     zz_mum(i_z_max)=0.D0
     Indice(i_x_max,i_z_max)=1.D0
  ENDIF

  CALL compute_dispersion(switch_dispersion)
  CALL compute_parameters
  CALL write_listingfile
  DO p=0,num_proc-1
     CALL calc_startingfield(switch_start,p)
     CALL write_startingfile(p)
  ENDDO
  
END PROGRAM make_start
