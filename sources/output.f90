MODULE output
  USE fields
  USE parameters
  USE mpi_stuff
  USE long_step
  USE run_status
  USE normalization
  USE HDF5
  USE HDF5_helper
CONTAINS
  
  SUBROUTINE write_output
    USE fft
    IMPLICIT NONE
    
    ! General purpose variables: looping, dummy variables
    INTEGER(4) k1,k2
    REAL(4) dumr4
    INTEGER(HSIZE_T), DIMENSION(1):: dumh51D, dumh51D2
    
    INTEGER :: field_dimensions ! Dataset rank & # of points in z
    ! the kind of this variable has to correspond with the precision stored in HDF5-file
    REAL(4), ALLOCATABLE :: fields_array(:,:,:), plasma_array(:,:,:), spect_array(:,:,:) 
    
    INTEGER(HSIZE_T)               :: r_offset

    INTEGER(4) j,k,l
    REAL(8) rhotemp,r,mpa
    COMPLEX(8) help
    CHARACTER*10 iz
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: group_id      ! Group identifier 
    INTEGER(HID_T) :: h5parameters  ! Property list identifier 
    INTEGER(HSIZE_T), DIMENSION(3) :: dims
    INTEGER(HSIZE_T), DIMENSION(3) :: ccount  
    INTEGER(HSIZE_T), DIMENSION(3) :: offset 
    INTEGER :: error
    CHARACTER(LEN=15) :: h5_filename="results.h5"
    CHARACTER(LEN=15) :: groupname="outputs"
    CHARACTER(LEN=25) :: field_dset_name="outputs/output_field"
    CHARACTER(LEN=25) :: plasma_dset_name="outputs/output_plasma"
    CHARACTER(LEN=25) :: spect_1d_dset_name="outputs/spect_1d"
  

      
    field_dimensions = 3
    allocate(fields_array(1,dim_r_local,dim_t))

    r_offset = dim_r_start(num_proc)-1
    DO k1=1, dim_r_local
    DO k2=1, dim_t
      fields_array(1,k1,k2) = REAL( REAL( (efield_factor*efield_osc(k2)*e(k2,r_offset+k1)) ) , 4 ) ! SINGLE PRECISION, corresponding H5T_NATIVE_REAL (REAL(.,8) corresponds to H5T_NATIVE_DOUBLE)
      ! e(t,r)
    ENDDO
    ENDDO

    allocate(plasma_array(1,dim_r_local,dim_t))
    k1 = 1
    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
      e_2=ABS(e(1:dim_t,l))**2
      e_2KK=e_2**KK
      rhotemp=rho0
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
      plasma_array(1,k1,:) = REAL(e_2KKm2,4)
      k1 = k1 + 1
    ENDDO

    dims = (/int(Nz_points,HSIZE_T),int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/)
    offset = (/int(output_write_count-1,HSIZE_T),int(dim_r_start(num_proc)-1,HSIZE_T),int(0,HSIZE_T)/)
    ccount = (/int(1,HSIZE_T), int(dim_r_local,HSIZE_T) , int(dim_t,HSIZE_T)/)
    CALL h5open_f(error) 
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) ! create HDF5 access parameters
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! set parameters for MPI access
    CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) ! Open collectivelly the file
    CALL h5pclose_f(h5parameters,error) ! close the parameters
    IF ( output_write_count == 1) THEN 
      !Create group for the output
      CALL h5gcreate_f(file_id, groupname, group_id, error) 
      CALL h5gclose_f(group_id, error)
          
      ! Call writing routine
      CALL create_3D_array_real_dset_p(file_id, field_dset_name, fields_array, dims, offset, ccount)
      CALL create_3D_array_real_dset_p(file_id, plasma_dset_name, plasma_array, dims, offset, ccount)

      ! Terminate
      CALL h5fclose_f(file_id,error)

      IF (my_rank.EQ.0) THEN
        CALL h5fopen_f (h5_filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file.
        CALL h5_add_units_1D(file_id, field_dset_name, '[V/m]') ! add units
        CALL h5fclose_f(file_id, error)
      ENDIF ! single-write end
    ELSE !!!! APPENDING THE DATA IN NEXT ITERATIONS
      CALL write_hyperslab_to_dset_p(file_id, field_dset_name, fields_array, offset, ccount)
      CALL write_hyperslab_to_dset_p(file_id, plasma_dset_name, plasma_array, offset, ccount)
      CALL h5fclose_f(file_id,error)
    ENDIF


    deallocate(fields_array)
    deallocate(plasma_array)
    CALL h5close_f(error)
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
      dims = (/int(Nz_points,HSIZE_T),int(5,HSIZE_T), int(dim_t,HSIZE_T)/)
      offset = (/int(output_write_count-1,HSIZE_T),int(0,HSIZE_T),int(0,HSIZE_T)/)
      ccount = (/int(1,HSIZE_T),int(5,HSIZE_T),int(dim_t,HSIZE_T)/)
      allocate(spect_array(1,5,1:dim_t))
      CALL h5open_f(error)
      CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error)
      DO j=1,dim_t
         spect_array(1,1,j) = REAL(k_t*REAL(j-1-dim_th,8)+omega_uppe,4) 
         spect_array(1,2,j) = REAL(e_2KK(j),4)
         spect_array(1,3,j) = REAL(e_2(j),4)
         spect_array(1,4,j) = REAL(ABS(etemp(j,1))**2,4)
         spect_array(1,5,j) = REAL(ATAN2(AIMAG(etemp(j,1)),REAL(etemp(j,1))+1.D-20),4)
      ENDDO
      IF ( output_write_count .EQ. 1) THEN
        CALL create_3D_array_real_dset(file_id, spect_1d_dset_name, spect_array, dims, offset, ccount)
      ELSE
        CALL write_hyperslab_to_dset(file_id, spect_1d_dset_name, spect_array, offset, ccount)
      ENDIF
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error) ! close the HDF5 workspace
      deallocate(spect_array)
    ENDIF
    output_write_count = output_write_count + 1 !increase counter in all cases
    RETURN
  END SUBROUTINE  write_output

  SUBROUTINE matlab_out
    USE fft
    IMPLICIT NONE

    INTEGER(4) j,k,l
    REAL(8) rhotemp,r,mpa
    COMPLEX(8) help
    CHARACTER*10 iz,filename
        

    WRITE(iz,920) z
    DO k=1,10
       IF (iz(k:k).EQ.' ') iz(k:k)='0'
       IF (iz(k:k).EQ.'.') iz(k:k)='_'
    ENDDO
    IF (my_rank.EQ.0) THEN
       filename='non'
       OPEN(unit_logfile,FILE='MERGE_RAD.LOG',STATUS='UNKNOWN')
       DO
          READ(unit_logfile,*,END=999) filename
       ENDDO
999    CONTINUE
       CLOSE(unit_logfile)
    ENDIF

    CALL MPI_BCAST(filename,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    IF (filename.NE.iz) THEN

       IF(my_rank.EQ.0) THEN
          OPEN(unit_logfile,FILE='MERGE_RAD.LOG',STATUS='UNKNOWN',POSITION='APPEND')
          WRITE(unit_logfile,*) iz
          CLOSE(unit_logfile)
       ENDIF
       OPEN(unit_field,FILE=iz//'_FIELD_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
       WRITE(unit_field) dim_t,dim_r,num_proc
       WRITE(unit_field) REAL(delta_t,4),REAL(delta_r,4),REAL(tlo,4)
       DO j=dim_r_start(num_proc),dim_r_end(num_proc)
          WRITE(unit_field) CMPLX(e(1:dim_t,j),KIND=4)
       ENDDO
       CLOSE(unit_field)

       OPEN(unit_field,FILE=iz//'_PLASMA_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
       WRITE(unit_field) dim_t,dim_r,num_proc
       WRITE(unit_field) REAL(delta_t,4),REAL(delta_r,4),REAL(tlo,4)
       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
          e_2=ABS(e(1:dim_t,l))**2
          e_2KK=e_2**KK
          rhotemp=rho0
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
          WRITE(unit_field) REAL(e_2KKm2,4)
       ENDDO
       CLOSE(unit_field)

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
          OPEN(unit_field,FILE=iz//'_spect_1d.dat',STATUS='UNKNOWN',RECL=2**10)
          DO j=1,dim_t
             WRITE(unit_field,*) REAL(k_t*REAL(j-1-dim_th,8)+omega_uppe,4),REAL(e_2KK(j),4),REAL(e_2(j),4), &
                  REAL(ABS(etemp(j,1))**2,4),REAL(ATAN2(AIMAG(etemp(j,1)),REAL(etemp(j,1))+1.D-20),4)
          ENDDO
          CLOSE(unit_field)
       ENDIF
    ENDIF

920 FORMAT (F10.6)
    RETURN
  END SUBROUTINE  matlab_out

  SUBROUTINE  field_out
    USE ppt
    IMPLICIT  NONE

    INTEGER(4) j,k,i_x,i_z
    CHARACTER*10  iz,filename
    CHARACTER*15  id

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
       IF(my_rank.EQ.0) THEN
          OPEN(unit_logfile,FILE='PROP_RAD.LOG',STATUS='UNKNOWN',POSITION='APPEND')
          WRITE(unit_logfile,*) iz    
          CLOSE(unit_logfile)
       ENDIF
       OPEN(unit_field,FILE=iz//'_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
       id='num_proc'
       WRITE(unit_field) id,num_proc
       id='dim_t'
       WRITE(unit_field) id,dim_t
       id='dim_r'
       WRITE(unit_field) id,dim_r
       id='rek0'
       WRITE(unit_field) id,rek0
       id='rekp'
       WRITE(unit_field) id,rekp
       id='c3'
       WRITE(unit_field) id,c3
       id='c5'
       WRITE(unit_field) id,c5
       id='gamma1'
       WRITE(unit_field) id,gamma1
       id='gamma2'
       WRITE(unit_field) id,gamma2
       id='muk'
       WRITE(unit_field) id,muk
       id='betainv2KK'
       WRITE(unit_field) id,beta_inv_2KK
       id='KK'
       WRITE(unit_field) id,KK
       id='rho0'
       WRITE(unit_field) id,rho0
       id='nu'
       WRITE(unit_field) id,nu
       id='alpha'
       WRITE(unit_field) id,alpha
       id='alphaquad'
       WRITE(unit_field) id,alphaquad
       id='rhoat_inv'
       WRITE(unit_field) id,rhoat_inv
       id='xdk'
       WRITE(unit_field) id,xdk
       id='tdk'
       WRITE(unit_field) id,tdk
       id='raman'
       WRITE(unit_field) id,raman
       id='omega'
       WRITE(unit_field) id,omega
       id='komega'
       WRITE(unit_field) id,komega(1:dim_t)
       id='NN'
       WRITE(unit_field) id,NN
       id='eta1'
       WRITE(unit_field) id,eta1
       id='eta2'
       WRITE(unit_field) id,eta2
       id='lt'
       WRITE(unit_field) id,lt
       id='lr'
       WRITE(unit_field) id,lr
       id='proplength'
       WRITE(unit_field) id,proplength
       id='outlength'
       WRITE(unit_field) id,outlength
       id='delta_z'
       WRITE(unit_field) id,delta_z
       id='z'
       WRITE(unit_field) id,z
       id='z_out'
       WRITE(unit_field) id,z_out
       id='rfil'
       WRITE(unit_field) id,rfil
       id='switch_rho'
       WRITE(unit_field) id,switch_rho
       id='switchKerr'
       WRITE(unit_field) id,switch_dKerr
       id='switch_T'
       WRITE(unit_field) id,switch_T
       id='absorb'
       WRITE(unit_field) id,absorb
       id='increase'
       WRITE(unit_field) id,increase
       id='decrease'
       WRITE(unit_field) id,decrease
       id='rhodist'
       WRITE(unit_field) id,rhodist
       id='timelimit'
       WRITE(unit_field) id,timelimit
       id='photenergy'
       WRITE(unit_field) id,photon_energy
       id='pulsedurat'
       WRITE(unit_field) id,pulse_duration
       id='critpower'
       WRITE(unit_field) id,critical_power
       id='beam_waist'
       WRITE(unit_field) id,beam_waist
       id='ionpot'
       WRITE(unit_field) id,ionisation_potential
       id='rescharge'
       WRITE(unit_field) id,residue_charge
       id='n0_indice'
       WRITE(unit_field) id,n0_indice
       id='critdens'
       WRITE(unit_field) id,critical_density
       id='atomdens'
       WRITE(unit_field) id,atomic_density
       id='reducmass'
       WRITE(unit_field) id,reduced_mass
       id='angmom'
       WRITE(unit_field) id,angular_momentum
       id='KKp'
       WRITE(unit_field) id,KKp
       id='beta_inv_2KKp'
       WRITE(unit_field) id,beta_inv_2KKp
       id='mukp'
       WRITE(unit_field) id,mukp
       id='beta_inv_2'
       WRITE(unit_field) id,beta_inv_2
       id='mu'
       WRITE(unit_field) id,mu
       id='KKpp'
       WRITE(unit_field) id,KKpp
       id='beta_inv_2KKpp'
       WRITE(unit_field) id,beta_inv_2KKpp
       id='mukpp'
       WRITE(unit_field) id,mukpp
       id='eti_ref'
       WRITE(unit_field) id,eti_ref
       id='exp_ref'
       WRITE(unit_field) id,exp_ref
       id='alpha1'
       WRITE(unit_field) id,alpha1
       id='alpha2'
       WRITE(unit_field) id,alpha2
       id='alphah'
       WRITE(unit_field) id,alphah
       id='rhosat'
       WRITE(unit_field) id,rhosat
       id='finished'
       WRITE(unit_field) id,finished
       id='omega_uppe'
       WRITE(unit_field) id,omega_uppe
       id='gamma1e'
       WRITE(10) id,gamma1e
       id='nuO2'
       WRITE(10) id,nuO2
       id='nuN2'
       WRITE(10) id,nuN2
       id='T_init_eV_phys'
       WRITE(10) id,T_init_eV_phys
       id='nukB'
       WRITE(10) id,nukB
       id='nucp'
       WRITE(10) id,nucp
       id='nucO2'
       WRITE(10) id,nucO2
       id='nucN2'
       WRITE(10) id,nucN2
       id='rhoat_N2_inv'
       WRITE(10) id,rhoat_N2_inv
       id='ionpotN2'
       WRITE(10) id,ionisation_potential_N2
       id='rescharge_N2'
       WRITE(10) id,residue_charge_N2
       id='atomdens_N2'
       WRITE(10) id,atomic_density_N2
       id='angmom_N2'
       WRITE(10) id,angular_momentum_N2
       id='startfield'
       WRITE(unit_field) id
       DO j=dim_r_start(num_proc),dim_r_end(num_proc)
          WRITE(unit_field) e(1:dim_t,j)
       ENDDO
       id='index'
       WRITE(unit_field) id
       WRITE(unit_field) i_x_max, i_z_max
       WRITE(unit_field) (xx(i_x),i_x=1,i_x_max)
       DO i_z = 1, i_z_max
          WRITE(unit_field) zz(i_z)
          WRITE(unit_field) (Indice_norm(i_x,i_z),i_x=1,i_x_max)
       ENDDO
       CLOSE(unit_field)
       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    ENDIF

920 FORMAT (F10.6)

    RETURN
  END SUBROUTINE field_out

END MODULE output
