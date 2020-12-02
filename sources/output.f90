! There are modules used for saving the ouputs
!
! - "write_output": it saves fields and plasma-density during the propagation 
!
! - "field_out": it stores the end-plane solution in the form of the code input for possible continuation
!
! - "linked_list_out": it saves the data buffered in the linked list at the end of the code
!
! The original procedures were developed by Stefan Skupin
! The change to HDF5 was designed by Jan Vabek and 
! co-implemented by Jakub Jelinek

MODULE output
  USE fields
  USE parameters
  USE mpi_stuff
  USE long_step
  USE run_status
  USE normalization
  USE HDF5
  USE HDF5_helper
  USE pre_ionised
CONTAINS
  
  SUBROUTINE write_output
    USE fft
    IMPLICIT NONE
    
    ! General purpose variables: looping, dummy variables
    INTEGER(4) k1,k2
    
    INTEGER :: field_dimensions ! Dataset rank & # of points in z
    ! the kind of this variable has to correspond with the precision stored in HDF5-file
    REAL(4), ALLOCATABLE :: fields_array(:,:,:), plasma_array(:,:,:), rgrid(:), tgrid(:)  
    REAL(4), ALLOCATABLE :: spect_array_1(:),  spect_array_2(:,:), spect_array_3(:,:), spect_array_4(:,:), spect_array_5(:,:) 
    INTEGER(HSIZE_T)               :: r_offset

    INTEGER(4) j,l
    REAL(8) rhotemp,r,mpa
    COMPLEX(8) help
    LOGICAL :: first = .FALSE.
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: group_id      ! Group identifier 
    INTEGER(HID_T) :: h5parameters  ! Property list identifier 
    INTEGER(HSIZE_T), DIMENSION(3) :: dims, offset, ccount
    INTEGER(HSIZE_T), DIMENSION(2) :: dims_2d, offset_2d, ccount_2d
    INTEGER                        :: error
    LOGICAL                        :: group_status
    CHARACTER(LEN=15) :: h5_filename="results.h5"
    CHARACTER(LEN=15) :: groupname="outputs"
    CHARACTER(LEN=25) :: field_dset_name="outputs/output_field"
    CHARACTER(LEN=25) :: plasma_dset_name="outputs/output_plasma"
    CHARACTER(LEN=25) :: spect_1d_dset_name_1="outputs/omegagrid"
    CHARACTER(LEN=60) :: spect_1d_dset_name_2="outputs/spectral_intensity_integrated_over_the_numerical_box"
    CHARACTER(LEN=69) :: spect_1d_dset_name_3="outputs/spectral_intensity_integrated_over_cylinder_with_radius_rfill"
    CHARACTER(LEN=26) :: spect_1d_dset_name_4="outputs/spectral_intensity"
    CHARACTER(LEN=22) :: spect_1d_dset_name_5="outputs/spectral_phase"
    CHARACTER(LEN=13) :: zgrid_dset_name = "outputs/zgrid"
    CHARACTER(LEN=13) :: tgrid_dset_name = "outputs/tgrid"
    CHARACTER(LEN=13) :: rgrid_dset_name = "outputs/rgrid"

      
    field_dimensions = 3
    allocate(fields_array(1,dim_r_local,dim_t))

    r_offset = dim_r_start(num_proc)-1
    DO k1=1, dim_r_local
    DO k2=1, dim_t
      fields_array(1,k1,k2) = REAL( efield_factor*REAL( (efield_osc(k2)*e(k2,r_offset+k1)) ) , 4 ) ! SINGLE PRECISION, corresponding H5T_NATIVE_REAL (REAL(.,8) corresponds to H5T_NATIVE_DOUBLE)
      ! e(t,r)
    ENDDO
    ENDDO

    allocate(plasma_array(1,dim_r_local,dim_t))
    k1 = 1
    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
      e_2=ABS(e(1:dim_t,l))**2
      e_2KK=e_2**KK
      IF (apply_pre_ionisation) THEN
         rhotemp = initial_electron_density_tip(r,z,l,dim_r_start(num_proc))
      ELSE
         rhotemp = 0.D0
      ENDIF
      rhompi=0.D0
      rho1=0.D0
      rho2=0.D0
      rhoth=0.D0
      rhotr=0.D0
      rhofh=0.D0
      rhoslg2=0.D0
      rhoav=0.D0
      rhoO2=rho0
      rhoN2=0.D0
      Tev=T_init_eV_phys
      DO j=1,dim_t
         e_2KKm2(j)=rhotemp
         IF (j.NE.dim_t) THEN
            CALL calc_rho(rhotemp,mpa,e_2(j),e_2(j+1))
         ENDIF
      ENDDO
      plasma_array(1,k1,:) = REAL(plasma_normalisation_factor_m3*e_2KKm2,4) ! SI units
      !plasma_array(1,k1,:) = REAL(e_2KKm2,4) ! computational units
      k1 = k1 + 1
    ENDDO

    IF ( output_write_count .EQ. 1) THEN
      first = .TRUE.
    ELSE
      first = .FALSE.
    ENDIF

    ! Calculate dimensions for the field to be preallocated, the offset and the hyperslab size
    dims = (/int(Nz_points,HSIZE_T),int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/)
    offset = (/int(output_write_count-1,HSIZE_T),int(dim_r_start(num_proc)-1,HSIZE_T),int(0,HSIZE_T)/)
    ccount = (/int(1,HSIZE_T), int(dim_r_local,HSIZE_T) , int(dim_t,HSIZE_T)/)
    CALL h5open_f(error) 
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) ! create HDF5 access parameters
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! set parameters for MPI access
    CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) ! Open collectivelly the file
    CALL h5pclose_f(h5parameters,error) ! close the parameters


    IF ( first ) THEN
      !Create group for the output if it does not already exist
      CALL h5lexists_f(file_id, groupname, group_status, error)
      IF ( group_status .EQV. .FALSE. ) THEN
        CALL h5gcreate_f(file_id, groupname, group_id, error) 
        CALL h5gclose_f(group_id, error)
      ENDIF
          
      ! Call writing routine
      CALL create_3D_array_real_dset_p(file_id, field_dset_name, fields_array, dims, offset, ccount)
      CALL create_3D_array_real_dset_p(file_id, plasma_dset_name, plasma_array, dims, offset, ccount)

      ! Terminate
      CALL h5fclose_f(file_id,error)

      IF (my_rank.EQ.0) THEN ! single-write start
        CALL h5fopen_f (h5_filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file.
        CALL h5_add_units_1D(file_id, field_dset_name, '[V/m]') 
	CALL h5_add_units_1D(file_id, plasma_dset_name, '[m^(-3)]') 

        ! r and t grids saved in the first run
	allocate(tgrid(dim_t),rgrid(dim_r)) ! space for grids: first itration, proc # 0
        DO k1=1, dim_t
          tgrid(k1) = REAL( tps*(tlo+REAL(k1,8)*delta_t) , 4)
        ENDDO
        DO k1=1, dim_r
          rgrid(k1) = REAL( w0m*(REAL(k1-1,8)*delta_r) , 4)
        ENDDO

        CALL create_dset(file_id, rgrid_dset_name, rgrid, dim_r)
        CALL h5_add_units_1D(file_id, rgrid_dset_name, '[m]')
        CALL create_dset(file_id, tgrid_dset_name, tgrid, dim_t)
        CALL h5_add_units_1D(file_id, tgrid_dset_name, '[m]')
        deallocate(tgrid,rgrid)

        CALL create_1D_dset_unlimited(file_id, zgrid_dset_name, (/REAL(four_z_Rayleigh*z,4)/), 1) ! the actual z-coordinate in SI units 
        CALL h5_add_units_1D(file_id, zgrid_dset_name, '[m]')

        CALL h5fclose_f(file_id, error) ! close the file
      ENDIF ! single-write end


    ELSE !!!! APPENDING THE DATA IN NEXT ITERATIONS
      CALL write_hyperslab_to_dset_p(file_id, field_dset_name, fields_array, offset, ccount)
      CALL write_hyperslab_to_dset_p(file_id, plasma_dset_name, plasma_array, offset, ccount)
      CALL h5fclose_f(file_id,error)

      IF (my_rank.EQ.0) THEN ! only one worker is extending the zgrid
        CALL h5open_f(error)  !Initialize HDF5
        CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error)
        ! only z-grid in 1D
        !dumh51D = (/int(HDF5write_count-1,HSIZE_T)/) ! offset
        !dumh51D2 = (/int(1,HSIZE_T)/) ! count
        !dumr4(1) = REAL(four_z_Rayleigh*z,4) ! the actual z-coordinate in SI units 
        CALL extend_1D_dset_unlimited(file_id, zgrid_dset_name, (/REAL(four_z_Rayleigh*z,4)/), new_dims=(/int(output_write_count,HSIZE_T)/), & 
          memspace_dims=(/int(1,HSIZE_T)/), offset=(/int(output_write_count-1,HSIZE_T)/), hyperslab_size=(/int(1,HSIZE_T)/))
        CALL h5fclose_f(file_id,error)
      ENDIF ! single-write end
    ENDIF

    ! Deallocate no longer needed arrays
    deallocate(fields_array)
    deallocate(plasma_array)
    CALL h5close_f(error) ! Close fortran H5 interface

    etemp=CSHIFT(e,dim_t/2-1,1)
    CALL dfftw_execute(plan_spec)
    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
     DO  j=1,dim_th
        help=etemp(j+dim_t/2,l)
        etemp(j+dim_t/2,l)=etemp(j,l)
        etemp(j,l)=help
     ENDDO
    ENDDO

    e_2(1:dim_t)=0.D0
    e_2KKm2(1:dim_t)=0.D0
    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
      r=REAL(l-1)*delta_r
      e_2(1:dim_t)=e_2(1:dim_t)+ABS(etemp(1:dim_t,l))**2*REAL(l-1,8)
      IF (rfil.GT.r) e_2KKm2(1:dim_t)=e_2KKm2(1:dim_t)+ABS(etemp(1:dim_t,l))**2*REAL(l-1,8)
    ENDDO


    CALL MPI_REDUCE(e_2(1:dim_t),e_2KK(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(e_2KKm2(1:dim_t),e_2(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    IF (my_rank.EQ.0) THEN
      dims_2d = (/int(Nz_points,HSIZE_T), int(dim_t,HSIZE_T)/)
      offset_2d = (/int(output_write_count-1,HSIZE_T),int(0,HSIZE_T)/)
      ccount_2d = (/int(1,HSIZE_T),int(dim_t,HSIZE_T)/)
      IF ( first ) THEN
        allocate(spect_array_1(1:dim_t),spect_array_2(1,1:dim_t),spect_array_3(1,1:dim_t), &
          spect_array_4(1,1:dim_t),spect_array_5(1,1:dim_t))
      ELSE
        allocate(spect_array_2(1,1:dim_t),spect_array_3(1,1:dim_t), spect_array_4(1,1:dim_t),spect_array_5(1,1:dim_t))
      ENDIF
      CALL h5open_f(error)
      CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error)
      IF ( first ) THEN
        DO j=1,dim_t
          spect_array_1(j) = REAL(k_t*REAL(j-1-dim_th,8)+omega_uppe,4) 
          spect_array_2(1,j) = REAL(e_2KK(j),4)
          spect_array_3(1,j) = REAL(e_2(j),4)
          spect_array_4(1,j) = REAL(ABS(etemp(j,1))**2,4)
          spect_array_5(1,j) = REAL(ATAN2(AIMAG(etemp(j,1)),REAL(etemp(j,1))+1.D-20),4)
        ENDDO
        CALL create_1D_array_real_dset(file_id, spect_1d_dset_name_1, spect_array_1, dim_t)
	CALL h5_add_units_1D(file_id, spect_1d_dset_name_1, '[?]')
        CALL create_and_preallocate_2D_array_real_dset(file_id, spect_1d_dset_name_2, spect_array_2, dims_2d, offset_2d, ccount_2d)
	CALL h5_add_units_1D(file_id, spect_1d_dset_name_2, '[arb.u.]')
        CALL create_and_preallocate_2D_array_real_dset(file_id, spect_1d_dset_name_3, spect_array_3, dims_2d, offset_2d, ccount_2d)
        CALL h5_add_units_1D(file_id, spect_1d_dset_name_3, '[arb.u.]')
        CALL create_and_preallocate_2D_array_real_dset(file_id, spect_1d_dset_name_4, spect_array_4, dims_2d, offset_2d, ccount_2d)
        CALL h5_add_units_1D(file_id, spect_1d_dset_name_4, '[arb.u.]')
        CALL create_and_preallocate_2D_array_real_dset(file_id, spect_1d_dset_name_5, spect_array_5, dims_2d, offset_2d, ccount_2d)
        CALL h5_add_units_1D(file_id, spect_1d_dset_name_5, '[arb.u.]')
        deallocate(spect_array_1, spect_array_2, spect_array_3, spect_array_4, spect_array_5)
      ELSE
        DO j=1,dim_t
          spect_array_2(1,j) = REAL(e_2KK(j),4)
          spect_array_3(1,j) = REAL(e_2(j),4)
          spect_array_4(1,j) = REAL(ABS(etemp(j,1))**2,4)
          spect_array_5(1,j) = REAL(ATAN2(AIMAG(etemp(j,1)),REAL(etemp(j,1))+1.D-20),4)
        ENDDO
        CALL write_hyperslab_to_2D_dset(file_id, spect_1d_dset_name_2, spect_array_2, offset_2d, ccount_2d)
        CALL write_hyperslab_to_2D_dset(file_id, spect_1d_dset_name_3, spect_array_3, offset_2d, ccount_2d)
        CALL write_hyperslab_to_2D_dset(file_id, spect_1d_dset_name_4, spect_array_4, offset_2d, ccount_2d)
        CALL write_hyperslab_to_2D_dset(file_id, spect_1d_dset_name_5, spect_array_5, offset_2d, ccount_2d)
        deallocate(spect_array_2, spect_array_3, spect_array_4, spect_array_5)
      ENDIF
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error) ! close the HDF5 workspace
    ENDIF
    output_write_count = output_write_count + 1 !increase counter in all cases
    RETURN
  END SUBROUTINE  write_output

  SUBROUTINE  field_out
    USE ppt
    USE normalization
    USE HDF5
    USE HDF5_helper
    IMPLICIT  NONE

    INTEGER(4) k,k1,k2
    CHARACTER*10  iz,filename
    INTEGER(HSIZE_T)               :: r_offset
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: group_id      ! Group identifier 
    INTEGER(HID_T) :: field_group_id! Field out group identifier
    INTEGER(HID_T) :: h5parameters  ! Property list identifier 
    INTEGER(HSIZE_T), DIMENSION(2) :: dims, offset, ccount
    INTEGER :: error
    LOGICAL :: group_status
    CHARACTER(LEN=15) :: h5_filename="results.h5"
    CHARACTER(LEN=15) :: groupname="/outputs"
    CHARACTER(LEN=25) :: field_out_groupname="/outputs/field_out_group"
    INTEGER(HID_T) :: indexes_group_id
    CHARACTER(LEN=50) :: indexes_groupname="/outputs/field_out_group/indexes_group"    
    REAL(4), ALLOCATABLE :: real_e(:,:),imag_e(:,:)
    
    WRITE(iz,920) z
    DO  k=1,10
       IF (iz(k:k).EQ.' ') iz(k:k)='0'
       IF (iz(k:k).EQ.'.') iz(k:k)='_'
    ENDDO
    IF (my_rank.EQ.0) THEN
       filename='non'
       OPEN(unit_logfile,FILE='PROP_RAD.LOG',STATUS='UNKNOWN')
       DO
          READ(unit_logfile,*,END=999) filename
       ENDDO
999    CONTINUE
       CLOSE(unit_logfile)
    ENDIF
    CALL MPI_BCAST(filename,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    IF (filename.NE.iz) THEN
       CALL h5open_f(error)
       CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) ! create HDF5 access parameters
       CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! set parameters for MPI access
       CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) ! Open collectivelly the file
       CALL h5pclose_f(h5parameters,error) ! close the parameters

       !Create group for the output if it does not already exist
       CALL h5lexists_f(file_id, groupname, group_status, error)
       IF ( group_status .EQV. .FALSE. ) THEN
         CALL h5gcreate_f(file_id, groupname, group_id, error) 
       ELSE
         CALL h5gopen_f(file_id, groupname, group_id, error)
       ENDIF
       CALL h5lexists_f(file_id, field_out_groupname, group_status, error)
       IF ( group_status.EQV..FALSE.) THEN
         CALL h5gcreate_f(file_id, field_out_groupname, field_group_id, error) 
       ELSE
         CALL h5gopen_f(file_id, field_out_groupname, field_group_id, error)
       ENDIF
       IF(my_rank.EQ.0) THEN
          OPEN(unit_logfile,FILE='PROP_RAD.LOG',STATUS='UNKNOWN',POSITION='APPEND')
          WRITE(unit_logfile,*) iz    
          CLOSE(unit_logfile)
       ENDIF
       CALL create_dset(field_group_id,'num_proc',num_proc)
       CALL create_dset(field_group_id,'dim_t',dim_t)
       CALL create_dset(field_group_id,'dim_r',dim_r)
       CALL create_dset(field_group_id,'rek0',rek0)
       CALL create_dset(field_group_id,'rekp',rekp)
       CALL create_dset(field_group_id,'c3',c3)
       CALL create_dset(field_group_id,'c5',c5)
       CALL create_dset(field_group_id,'gamma1',gamma1)
       CALL create_dset(field_group_id,'gamma2',gamma2)
       CALL create_dset(field_group_id,'muk',muk)
       CALL create_dset(field_group_id,'betainv2KK',beta_inv_2KK)
       CALL create_dset(field_group_id,'KK',KK)
       CALL create_dset(field_group_id,'rho0',rho0)
       CALL create_dset(field_group_id,'nu',nu)
       CALL create_dset(field_group_id,'alpha',alpha)
       CALL create_dset(field_group_id,'alphaquad',alphaquad)
       CALL create_dset(field_group_id,'rhoat_inv',rhoat_inv)
       CALL create_dset(field_group_id,'xdk',xdk)
       CALL create_dset(field_group_id,'tdk',tdk)
       CALL create_dset(field_group_id,'raman',raman)
       CALL create_dset(field_group_id,'omega',omega)
       CALL create_dset(field_group_id,'komega',komega(1:dim_t),dim_t)
       CALL create_dset(field_group_id,'NN',NN)
       CALL create_dset(field_group_id,'eta1',eta1)
       CALL create_dset(field_group_id,'eta2',eta2)
       CALL create_dset(field_group_id,'lt',lt)
       CALL create_dset(field_group_id,'lr',lr)
       CALL create_dset(field_group_id,'proplength',proplength)
       CALL create_dset(field_group_id,'outlength',outlength)
       CALL create_dset(field_group_id,'delta_z',delta_z)
       CALL create_dset(field_group_id,'z',z)
       CALL create_dset(field_group_id,'z_out',z_out)
       CALL create_dset(field_group_id,'rfil',rfil)
       CALL create_dset(field_group_id,'switch_rho',switch_rho)
       CALL create_dset(field_group_id,'switchKerr',switch_dKerr)
       CALL create_dset(field_group_id,'switch_T',switch_T)
       CALL create_dset(field_group_id,'absorb',absorb)
       CALL create_dset(field_group_id,'increase',increase)
       CALL create_dset(field_group_id,'decrease',decrease)
       CALL create_dset(field_group_id,'rhodist',rhodist)
       CALL create_dset(field_group_id,'timelimit',timelimit)
       CALL create_dset(field_group_id,'photenergy',photon_energy)
       CALL create_dset(field_group_id,'pulsedurat',pulse_duration)
       CALL create_dset(field_group_id,'critpower',critical_power)
       CALL create_dset(field_group_id,'beam_waist',beam_waist)
       CALL create_dset(field_group_id,'ionpot',ionisation_potential)
       CALL create_dset(field_group_id,'rescharge',residue_charge)
       CALL create_dset(field_group_id,'n0_indice',n0_indice)
       CALL create_dset(field_group_id,'critdens',critical_density)
       CALL create_dset(field_group_id,'atomdens',atomic_density)
       CALL create_dset(field_group_id,'reducmass',reduced_mass)
       CALL create_dset(field_group_id,'angmom',angular_momentum)
       CALL create_dset(field_group_id,'KKp',KKp)
       CALL create_dset(field_group_id,'beta_inv_2KKp',beta_inv_2KKp)
       CALL create_dset(field_group_id,'mukp',mukp)
       CALL create_dset(field_group_id,'beta_inv_2',beta_inv_2)
       CALL create_dset(field_group_id,'mu',mu)
       CALL create_dset(field_group_id,'KKpp',KKpp)
       CALL create_dset(field_group_id,'beta_inv_2KKpp',beta_inv_2KKpp)
       CALL create_dset(field_group_id,'mukpp',mukpp)
       CALL create_dset(field_group_id,'eti_ref',eti_ref)
       CALL create_dset(field_group_id,'exp_ref',exp_ref)
       CALL create_dset(field_group_id,'alpha1',alpha1)
       CALL create_dset(field_group_id,'alpha2',alpha2)
       CALL create_dset(field_group_id,'alphah',alphah)
       CALL create_dset(field_group_id,'rhosat',rhosat)
       CALL create_dset(field_group_id,'finished',finished)
       CALL create_dset(field_group_id,'omega_uppe',omega_uppe)
       CALL create_dset(field_group_id,'gamma1e',gamma1e)
       CALL create_dset(field_group_id,'nuO2',nuO2)
       CALL create_dset(field_group_id,'nuN2',nuN2)
       CALL create_dset(field_group_id,'T_init_eV_phys',T_init_eV_phys)
       CALL create_dset(field_group_id,'nukB',nukB)
       CALL create_dset(field_group_id,'nucp',nucp)
       CALL create_dset(field_group_id,'nucO2',nucO2)
       CALL create_dset(field_group_id,'nucN2',nucN2)
       CALL create_dset(field_group_id,'rhoat_N2_inv',rhoat_N2_inv)
       CALL create_dset(field_group_id,'ionpotN2',ionisation_potential_N2)
       CALL create_dset(field_group_id,'rescharge_N2',residue_charge_N2)
       CALL create_dset(field_group_id,'atomdens_N2',atomic_density_N2)
       CALL create_dset(field_group_id,'angmom_N2',angular_momentum_N2)
       efield_factor = SQRT(critical_power*1.D9*3.D8*4.D0*3.1415D-7/(4.D0*3.1415D0*beam_waist**2*1.D-4*2.D0*n0_indice))*2.D0 ! normalization factor electric field V/m
       ALLOCATE(real_e(dim_t,dim_r/num_proc),imag_e(dim_t,dim_r/num_proc))
       r_offset = dim_r/num_proc*my_rank
       DO k1=1, dim_t
         DO k2=1, dim_r/num_proc
            real_e(k1,k2) = REAL(REAL((efield_factor*efield_osc(k2)*e(k1,k2+r_offset))),4)
            imag_e(k1,k2) = REAL(AIMAG((efield_factor*efield_osc(k2)*e(k1,k2+r_offset))),4)
         ENDDO
       ENDDO
       dims = (/int(dim_t,HSIZE_T), int(dim_r,HSIZE_T)/)
       offset = (/int(0,HSIZE_T),int(dim_r/num_proc*my_rank,HSIZE_T)/)
       ccount = (/int(dim_t,HSIZE_T),int(dim_r/num_proc,HSIZE_T)/)
       CALL create_2D_array_real_dset_p(field_group_id, "startfield_r", real_e, dims, offset, ccount)
       CALL create_2D_array_real_dset_p(field_group_id, "startfield_i", imag_e, dims, offset, ccount)
       DEALLOCATE(real_e,imag_e)
       CALL h5gcreate_f(file_id, indexes_groupname, indexes_group_id, error)
       CALL create_dset(indexes_group_id, "r_vector", REAL(xx(1:i_x_max),4), i_x_max)
       CALL create_dset(indexes_group_id, "z_vector", REAL(zz(1:i_z_max),4), i_z_max)
       CALL create_2D_array_real_dset(indexes_group_id, "indexes", REAL(Indice_norm(1:i_x_max, 1:i_z_max),8), i_x_max, i_z_max)
       CALL h5gclose_f(indexes_group_id, error)
       CALL h5gclose_f(field_group_id, error)
       CALL h5gclose_f(group_id, error)
       CALL h5fclose_f(file_id, error)
       CALL h5close_f(error)
    ENDIF

920 FORMAT (F10.6)

    RETURN
  END SUBROUTINE field_out
  
  SUBROUTINE linked_list_out
    USE linked_list
    USE longstep_vars
    USE HDF5_helper

    IMPLICIT NONE
    
    INTEGER        :: i
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: h5parameters  ! Property list identifier 
    INTEGER(HSIZE_T), DIMENSION(2) :: dims, offset, ccount
    INTEGER                        :: error
    CHARACTER(LEN=15) :: h5_filename="results.h5"
    CHARACTER(LEN=25) :: fluence_dset_name="longstep/fluence"
    CHARACTER(LEN=23) :: plasma_channel_dset_name="longstep/plasma_channel"
    CHARACTER(LEN=22) :: losses_plasma_dset_name="longstep/losses_plasma"
    CHARACTER(LEN=26) :: losses_ionization_dset_name="longstep/losses_ionization"
    REAL(4), ALLOCATABLE :: fluence_part(:,:)
    REAL(4), ALLOCATABLE :: plasma_channel_part(:,:)
    REAL(4), ALLOCATABLE :: losses_plasma_part(:,:)
    REAL(4), ALLOCATABLE :: losses_ionization_part(:,:)
    TYPE(list_t), POINTER      :: next_fluence_ll
    TYPE(list_t), POINTER      :: next_plasma_channel_ll
    TYPE(list_t), POINTER      :: next_losses_plasma_ll
    TYPE(list_t), POINTER      :: next_losses_ionization_ll
    
    ! Open HDF5 archive collectivelly
    CALL h5open_f(error)
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) ! create HDF5 access parameters
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! set parameters for MPI access
    CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) ! Open collectivelly the file
    CALL h5pclose_f(h5parameters,error) ! close the parameters
    
    ! Calculate the size of the dataset
    dims = (/int(length_of_linked_list,HSIZE_T),int(dim_r,HSIZE_T)/)
    ALLOCATE(fluence_part(1,dim_r_local))
    ALLOCATE(plasma_channel_part(1,dim_r_local))
    ALLOCATE(losses_plasma_part(1,dim_r_local))
    ALLOCATE(losses_ionization_part(1,dim_r_local))

    DO i=1,length_of_linked_list
      offset = (/int(i-1,HSIZE_T),int(dim_r_start(num_proc)-1,HSIZE_T)/)
      ccount = (/int(1,HSIZE_T), int(dim_r_local,HSIZE_T)/)
      IF (i .EQ. 1) THEN
        fluence_part(1,:) = transfer(list_get(fluence_ll), fluence_part(1,:))   
        next_fluence_ll => list_next(fluence_ll)
        CALL create_2D_array_real_dset_p(file_id, fluence_dset_name, fluence_part, dims, offset, ccount)

        plasma_channel_part(1,:) = transfer(list_get(plasma_channel_ll), plasma_channel_part(1,:))   
        next_plasma_channel_ll => list_next(plasma_channel_ll)
        CALL create_2D_array_real_dset_p(file_id, plasma_channel_dset_name, plasma_channel_part, dims, offset, ccount)
        
        losses_plasma_part(1,:) = transfer(list_get(losses_plasma_ll), losses_plasma_part(1,:))   
        next_losses_plasma_ll => list_next(losses_plasma_ll)
        CALL create_2D_array_real_dset_p(file_id, losses_plasma_dset_name, losses_plasma_part, dims, offset, ccount)
        
        losses_ionization_part(1,:) = transfer(list_get(losses_ionization_ll), losses_ionization_part(1,:))   
        next_losses_ionization_ll => list_next(losses_ionization_ll)
        CALL create_2D_array_real_dset_p(file_id, losses_ionization_dset_name, losses_ionization_part, dims, offset, ccount)
      ELSE
        fluence_part(1,:) = transfer(list_get(next_fluence_ll), fluence_part(1,:))
        CALL write_hyperslab_to_2D_dset(file_id, fluence_dset_name, fluence_part, offset, ccount)
        next_fluence_ll => list_next(next_fluence_ll)
        
        plasma_channel_part(1,:) = transfer(list_get(next_plasma_channel_ll), plasma_channel_part(1,:))   
        CALL write_hyperslab_to_2D_dset(file_id, plasma_channel_dset_name, plasma_channel_part, offset, ccount)
        next_plasma_channel_ll => list_next(next_plasma_channel_ll)
        
        losses_plasma_part(1,:) = transfer(list_get(next_losses_plasma_ll), losses_plasma_part(1,:))   
        CALL write_hyperslab_to_2D_dset(file_id, losses_plasma_dset_name, losses_plasma_part, offset, ccount)
        next_losses_plasma_ll => list_next(next_losses_plasma_ll)
        
        losses_ionization_part(1,:) = transfer(list_get(next_losses_ionization_ll), losses_ionization_part(1,:))  
        CALL write_hyperslab_to_2D_dset(file_id, losses_ionization_dset_name, losses_ionization_part, offset, ccount)
        next_losses_ionization_ll => list_next(next_losses_ionization_ll)
      ENDIF
    END DO
    DEALLOCATE(fluence_part, plasma_channel_part, losses_plasma_part, losses_ionization_part)
    CALL list_free(fluence_ll)
    CALL list_free(plasma_channel_ll)
    CALL list_free(losses_plasma_ll)
    CALL list_free(losses_ionization_ll)
    
    ! Terminate collective access
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

    ! Add units (not implemented for collective access)
      IF (my_rank.EQ.0) THEN ! single-write start
        CALL h5open_f(error)
        CALL h5fopen_f (h5_filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file.
        CALL h5_add_units_1D(file_id, fluence_dset_name, '[C.U.]') 
	CALL h5_add_units_1D(file_id, plasma_channel_dset_name, '[C.U.]') 
        CALL h5_add_units_1D(file_id, losses_plasma_dset_name, '[C.U.]') 
	CALL h5_add_units_1D(file_id, losses_ionization_dset_name, '[C.U.]')
        CALL h5fclose_f(file_id, error) ! close the file
        CALL h5close_f(error)
      ENDIF 

  END SUBROUTINE linked_list_out

END MODULE output



!  SUBROUTINE matlab_out
!     USE fft
!     IMPLICIT NONE

!     INTEGER(4) j,k,l
!     REAL(8) rhotemp,r,mpa
!     COMPLEX(8) help
!     CHARACTER*10 iz,filename
        

!     WRITE(iz,920) z
!     DO k=1,10
!        IF (iz(k:k).EQ.' ') iz(k:k)='0'
!        IF (iz(k:k).EQ.'.') iz(k:k)='_'
!     ENDDO
!     IF (my_rank.EQ.0) THEN
!        filename='non'
!        OPEN(unit_logfile,FILE='MERGE_RAD.LOG',STATUS='UNKNOWN')
!        DO
!           READ(unit_logfile,*,END=999) filename
!        ENDDO
! 999    CONTINUE
!        CLOSE(unit_logfile)
!     ENDIF

!     CALL MPI_BCAST(filename,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!     IF (filename.NE.iz) THEN

!        IF(my_rank.EQ.0) THEN
!           OPEN(unit_logfile,FILE='MERGE_RAD.LOG',STATUS='UNKNOWN',POSITION='APPEND')
!           WRITE(unit_logfile,*) iz
!           CLOSE(unit_logfile)
!        ENDIF
!        OPEN(unit_field,FILE=iz//'_FIELD_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
!        WRITE(unit_field) dim_t,dim_r,num_proc
!        WRITE(unit_field) REAL(delta_t,4),REAL(delta_r,4),REAL(tlo,4)
!        DO j=dim_r_start(num_proc),dim_r_end(num_proc)
!           WRITE(unit_field) CMPLX(e(1:dim_t,j),KIND=4)
!        ENDDO
!        CLOSE(unit_field)

!        OPEN(unit_field,FILE=iz//'_PLASMA_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
!        WRITE(unit_field) dim_t,dim_r,num_proc
!        WRITE(unit_field) REAL(delta_t,4),REAL(delta_r,4),REAL(tlo,4)
!        DO l=dim_r_start(num_proc),dim_r_end(num_proc)
!           e_2=ABS(e(1:dim_t,l))**2
!           e_2KK=e_2**KK
!           rhotemp=rho0
!           rhompi=0.D0
!           rho1=0.D0
!           rho2=0.D0
!           rhoth=0.D0
!           rhotr=0.D0
!           rhofh=0.D0
!           rhoslg2=0.D0
!           rhoav=0.D0
!           rhoO2=rho0
!           rhoN2=0.D0
!           Tev=T_init_eV_phys
!           DO j=1,dim_t
!              e_2KKm2(j)=rhotemp
!              IF (j.NE.dim_t) THEN
!                 CALL calc_rho(rhotemp,mpa,e_2(j),e_2(j+1))
!              ENDIF
!           ENDDO
!           WRITE(unit_field) REAL(e_2KKm2,4)
!        ENDDO
!        CLOSE(unit_field)

!        etemp=CSHIFT(e,dim_t/2-1,1)
!        CALL dfftw_execute(plan_spec)
!        DO l=dim_r_start(num_proc),dim_r_end(num_proc)
!           DO  j=1,dim_th
!              help=etemp(j+dim_t/2,l)
!              etemp(j+dim_t/2,l)=etemp(j,l)
!              etemp(j,l)=help
!           ENDDO
!        ENDDO
       
!        e_2(1:dim_t)=0.D0
!        e_2KKm2(1:dim_t)=0.D0
!        DO l=dim_r_start(num_proc),dim_r_end(num_proc)
!           r=REAL(l-1)*delta_r
!           e_2(1:dim_t)=e_2(1:dim_t)+ABS(etemp(1:dim_t,l))**2*REAL(l-1,8)
!           IF (rfil.GT.r) e_2KKm2(1:dim_t)=e_2KKm2(1:dim_t)+ABS(etemp(1:dim_t,l))**2*REAL(l-1,8)
!        ENDDO


!        CALL MPI_REDUCE(e_2(1:dim_t),e_2KK(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!        CALL MPI_REDUCE(e_2KKm2(1:dim_t),e_2(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!        IF (my_rank.EQ.0) THEN
!           OPEN(unit_field,FILE=iz//'_spect_1d.dat',STATUS='UNKNOWN',RECL=2**10)
!           DO j=1,dim_t
!              WRITE(unit_field,*) REAL(k_t*REAL(j-1-dim_th,8)+omega_uppe,4),REAL(e_2KK(j),4),REAL(e_2(j),4), &
!                   REAL(ABS(etemp(j,1))**2,4),REAL(ATAN2(AIMAG(etemp(j,1)),REAL(etemp(j,1))+1.D-20),4)
!           ENDDO
!           CLOSE(unit_field)
!        ENDIF
!     ENDIF

! 920 FORMAT (F10.6)
!     RETURN
!   END SUBROUTINE  matlab_out
