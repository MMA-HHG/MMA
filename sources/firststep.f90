! There are preparatory subroutines called before the main propagation
!
! "calc_propagator" computes the matrix representation of propagation 
! operators according to various models (and is accessed if step-size changed)
! 
! "initialise" does all preparation calculation, mainly:
! Reads and calculate code parameters.
! Calls ionisation models.
!
! It provides an interface with other modules, so it was touched by
! the authors of most of them.

MODULE first_step
  USE fields
  USE parameters
  USE mpi_stuff
  USE run_status
  USE normalization
CONTAINS

  SUBROUTINE calc_propagator
    IMPLICIT NONE

    INTEGER(4)  :: j,k
    REAL(8) t

    delta_zh=0.5D0*delta_z

    SELECT CASE (switch_T)
    CASE(1)
       DO j=1,dim_th
          p_t(j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(dim_th+j))
          p_t(dim_th+j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(j))
       ENDDO
       hfac=1.D0
    CASE(2)
       DO j=1,dim_th
          p_t(j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(dim_th+j))
          p_t(dim_th+j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(j))
       ENDDO
       hfac=1.D0
    CASE(3)
       DO j=1,dim_th
          p_t(j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(dim_th+j))
          p_t(dim_th+j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(j))
       ENDDO
       DO j=1,dim_t
          t=tlo+REAL(j,8)*delta_t
          hfac(j,0)=exp(CMPLX(0.D0,(omega-omega_uppe)*t,8))
          hfac(j,1)=CMPLX(0.D0,c3i*delta_zh/3.D0)*exp(CMPLX(0.D0,(omega_uppe-3.D0*omega)*t,8))
          hfac(j,2)=CMPLX(0.D0,-c5*delta_zh/2.D0)*exp(CMPLX(0.D0,(omega_uppe-3.D0*omega)*t,8))
          hfac(j,3)=CMPLX(0.D0,-c5*delta_zh/10.D0)*exp(CMPLX(0.D0,(omega_uppe-5.D0*omega)*t,8))
          hfac(j,4)=exp(CMPLX(0.D0,(omega_uppe+omega)*t,8))
       ENDDO
    CASE(4)
       DO j=1,dim_th
          p_t(j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(dim_th+j))
          p_t(dim_th+j)=exp(CMPLX(0.D0,delta_zh,8)*komega_red(j))
       ENDDO
       hfac=1.D0
    END SELECT

    delta_rel=op_t_inv*delta_zh/delta_r**2
    DO k=dim_t_start(num_proc),dim_t_end(num_proc)  
       DU(1,k)=CMPLX(0.D0,-2.D0,8)*delta_rel(k)
       DO j=1,dim_r-2
          DU(j+1,k)=delta_rel(k)*CMPLX(0.D0,-1.D0*(0.25D0/REAL(j,8)+0.5D0),8)
          DL(j,k)=delta_rel(k)*CMPLX(0.D0,(0.25D0/REAL(j,8)-0.5D0),8)
       ENDDO
       D(1,k)=CMPLX(1.D0,0.D0,8)+CMPLX(0.D0,2.D0,8)*delta_rel(k)
       DO j=1,dim_r-2
          D(j+1,k)=CMPLX(1.D0,0.D0,8)+CMPLX(0.D0,1.D0,8)*delta_rel(k)
       ENDDO
       D(dim_r,k)=CMPLX(1.D0,0.D0,8)
    ENDDO

    RETURN
  END SUBROUTINE calc_propagator




  SUBROUTINE initialize
    USE ppt
    USE External_ionisation_table 
    USE fft
    USE HDF5
    USE HDF5_helper
    USE pre_ionised
    IMPLICIT NONE

    INTEGER(4)  j,k,help,i_x,i_z,k1
    REAL(8) absorb_factor
    LOGICAL ext
    CHARACTER*10 filename,id
    CHARACTER(LEN=10), PARAMETER :: hdf5_input = "results.h5"  ! File name for the HDF5 input file
    CHARACTER(LEN = *), PARAMETER :: output_groupname = "pre-processed" 
    CHARACTER(LEN = *), PARAMETER :: input_groupname = "inputs" 
    INTEGER(HID_T) :: file_id, group_id     ! File identifier 
    INTEGER        :: error
    REAL(8), ALLOCATABLE :: real_e(:,:),imag_e(:,:)
    REAL(8) :: PI

    PI = 4.d0 * Atan(1.d0) ! still local

    HDF5write_count = 1
    output_write_count = 1
    CALL MPI_Init(ierr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
    CALL MPI_Comm_size(MPI_COMM_WORLD,num_proc,ierr)
    IF (my_rank.EQ.0) THEN
       OPEN(unit_rho,FILE='STOP',STATUS='UNKNOWN')
       CLOSE(unit_rho)
       OPEN(unit_logfile,FILE='PROP_RAD.LOG',STATUS='OLD')
       DO
          READ(unit_logfile,*,END=999) filename
       ENDDO
999    CONTINUE
       CLOSE(unit_logfile)
    ENDIF

  IF (my_rank.EQ.0) THEN
    print *, "init started "
  ENDIF

    CALL MPI_BCAST(filename,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    WRITE(ip,930) my_rank
    DO k=1,3
       IF (ip(k:k).EQ.' ') ip(k:k)='0'
    ENDDO
   
! OPEN HDF5 interface
    CALL h5open_f(error) 
    CALL h5fopen_f (hdf5_input, H5F_ACC_RDONLY_F, file_id, error)

    ! direct code inputs
    CALL h5gopen_f(file_id, input_groupname, group_id, error) 
   ! CALL read_dset(group_id, 'ionised_atoms_relative_Kerr_response', ions_Kerr_ratio) 
    CALL h5gclose_f(group_id, error) ! all pre-processed inputs read
    ions_Kerr_ratio = 1.D0/3.D0

    ! inputs from the pre-processor
    CALL h5gopen_f(file_id, output_groupname, group_id, error) 
    CALL read_dset(group_id, 'num_proc', num_proc)
    CALL read_dset(group_id, 'dim_t', dim_t)
    CALL read_dset(group_id, 'dim_r', dim_r)
    CALL read_dset(group_id, 'rek0', rek0)
    CALL read_dset(group_id, 'rekp',rekp)
    CALL read_dset(group_id, 'c3',c3)
    CALL read_dset(group_id, 'c5',c5)
    CALL read_dset(group_id, 'gamma1',gamma1)
    CALL read_dset(group_id, 'gamma2',gamma2)
    CALL read_dset(group_id, 'muk',muk)
    CALL read_dset(group_id, 'beta_inv_2KK',beta_inv_2KK)
    CALL read_dset(group_id,'KK',KK)
    CALL read_dset(group_id, 'rho0',rho0)
    CALL read_dset(group_id, 'nu',nu)
    CALL read_dset(group_id, 'alpha',alpha)
    CALL read_dset(group_id, 'alphaquad',alphaquad)
    CALL read_dset(group_id, 'rhoat_inv',rhoat_inv)
    CALL read_dset(group_id, 'xdk',xdk)
    CALL read_dset(group_id, 'tdk',tdk)
    CALL read_dset(group_id, 'raman',raman)
    CALL read_dset(group_id, 'omega',omega)
    ALLOCATE(komega(dim_t),komega_red(dim_t))
    CALL read_dset(group_id,'komega',komega,dim_t)
    CALL read_dset(group_id, 'NN',NN)
    CALL read_dset(group_id, 'eta1',eta1)
    CALL read_dset(group_id, 'eta2',eta2)
    CALL read_dset(group_id, 'lt',lt)
    CALL read_dset(group_id, 'lr',lr)
    CALL read_dset(group_id, 'proplength',proplength)
    CALL read_dset(group_id, 'outlength',outlength)
    CALL read_dset(group_id, 'delta_z',delta_z)
    CALL read_dset(group_id, 'z',z)
    CALL read_dset(group_id, 'z_out',z_out)
    CALL read_dset(group_id, 'rfil',rfil)
    CALL read_dset(group_id,'switch_rho', switch_rho)
    switch_rho = 8
    CALL read_dset(group_id,'switchKerr', switch_dKerr)
    CALL read_dset(group_id,'switch_T',switch_T)
    CALL read_dset(group_id,'absorb',absorb)
    CALL read_dset(group_id, 'increase',increase)
    CALL read_dset(group_id, 'decrease',decrease)
    CALL read_dset(group_id,'rhodist',rhodist)
    CALL read_dset(group_id, 'timelimit',timelimit)
    CALL read_dset(group_id, 'photenergy',photon_energy)
    CALL read_dset(group_id, 'pulsedurat',pulse_duration)
    CALL read_dset(group_id, 'critpower',critical_power)
    CALL read_dset(group_id, 'beam_waist',beam_waist)
    CALL read_dset(group_id, 'ionpot',ionisation_potential)
    CALL read_dset(group_id, 'rescharge',residue_charge)
    CALL read_dset(group_id, 'n0_indice',n0_indice)
    CALL read_dset(group_id, 'critdens',critical_density)
    CALL read_dset(group_id, 'atomdens',atomic_density)
    CALL read_dset(group_id, 'reducmass',reduced_mass)
    CALL read_dset(group_id,'angmom',angular_momentum)
    CALL read_dset(group_id,'KKp',KKp)
    CALL read_dset(group_id, 'beta_inv_2KKp',beta_inv_2KKp)
    CALL read_dset(group_id, 'mukp',mukp)
    CALL read_dset(group_id, 'beta_inv_2',beta_inv_2)
    CALL read_dset(group_id, 'mu',mu)
    CALL read_dset(group_id,'KKpp',KKpp)
    CALL read_dset(group_id, 'beta_inv_2KKpp',beta_inv_2KKpp)
    CALL read_dset(group_id, 'mukpp',mukpp)
    CALL read_dset(group_id, 'eti_ref',eti_ref)
    CALL read_dset(group_id, 'exp_ref',exp_ref)
    CALL read_dset(group_id, 'alpha1',alpha1)
    CALL read_dset(group_id, 'alpha2',alpha2)
    CALL read_dset(group_id, 'alphah',alphah)
    CALL read_dset(group_id, 'rhosat',rhosat)
    CALL read_dset(group_id, 'finished',finished)
    CALL read_dset(group_id, 'omega_uppe',omega_uppe)
    CALL read_dset(group_id, 'gamma1e',gamma1e)
    CALL read_dset(group_id, 'nuO2',nuO2)
    CALL read_dset(group_id, 'nuN2',nuN2)
    CALL read_dset(group_id, 'T_init_eV_phys',T_init_eV_phys)
    CALL read_dset(group_id, 'nukB',nukB)
    CALL read_dset(group_id, 'nucp',nucp)
    CALL read_dset(group_id, 'nucO2',nucO2)
    CALL read_dset(group_id, 'nucN2',nucN2)
    CALL read_dset(group_id, 'rhoat_N2_inv',rhoat_N2_inv)
    CALL read_dset(group_id, 'ionpotN2',ionisation_potential_N2)
    CALL read_dset(group_id, 'rescharge_N2',residue_charge_N2)
    CALL read_dset(group_id, 'atomdens_N2',atomic_density_N2)
    CALL read_dset(group_id,'angmom_N2',angular_momentum_N2)

    ! Prepare the fourier transforms
    CALL fft_init

    delta_t=lt/REAL(dim_t,8)
    delta_r=lr/REAL(dim_r,8)
    delta_t_inv=1.D0/delta_t
    tlo=-0.5D0*lt
    count=0
    k_t=8.D0*DATAN(1.D0)/lt
    dim_th=dim_t/2
    DO j=1,dim_t
       komega_red(j)=komega(j)-CMPLX(rek0+rekp*(k_t*(REAL(j-dim_th-1,8))+omega_uppe-omega),0.D0,8)
    ENDDO
    SELECT CASE (switch_T)
    CASE(1)
       help=2*NINT(omega/k_t)
    CASE(2)
       help=2*NINT(omega/k_t)
    CASE(3)
       omega_offset(1)=NINT(omega_uppe/k_t)-dim_th
       omega_offset(2)=NINT(omega/k_t)-(NINT(omega_uppe/k_t)-dim_th)
       IF(c5.EQ.0.D0) THEN 
          help=3*NINT(omega/k_t)-(NINT(omega_uppe/k_t)-dim_th)+1
       ELSE
          help=5*NINT(omega/k_t)-(NINT(omega_uppe/k_t)-dim_th)+1
       ENDIF
    CASE(4)
       help=2*NINT(omega/k_t)
    END SELECT
    help=MIN(help,dim_t)
    delta_z_max=4.D0*delta_r**2
    IF (MAXVAL(ABS(REAL(komega_red(1:help)))).GT.0.D0) THEN
       delta_z_max=MIN(delta_z_max,4.D0/MAXVAL(ABS(REAL(komega_red(1:help)))))
    ENDIF
    IF (switch_T.EQ.3) THEN
       IF (ABS(REAL(komega_red(help))).GT.0.D0) THEN
          delta_z_max=MIN(delta_z_max,8.D-2*DATAN(1.D0)/ABS(REAL(komega_red(help))))
       ENDIF
       IF (ABS(rek0-rekp*omega).GT.0.D0) THEN
          delta_z_max=MIN(delta_z_max,2.D-2*DATAN(1.D0)/ABS(rek0-rekp*omega))
       ENDIF
    ENDIF
    delta_z=MIN(delta_z,2.D0*delta_z_max)
    IF(my_rank.EQ.0) THEN
       OPEN(unit_rho,FILE='ZSTEP.DAT',STATUS='UNKNOWN',POSITION='APPEND')
       WRITE(unit_rho,*) 'z=',REAL(z,4),' delta_z=',REAL(delta_z,4)
       CLOSE(unit_rho)
    ENDIF
    ALLOCATE(bound_t(dim_t))
    bound_t=1.D0
    IF (absorb.GT.0) THEN
       absorb_factor=5.3D0/REAL(absorb,8)
       DO j=1,MIN(absorb*10,dim_t-1)
          bound_t(j)=bound_t(j)-1.D0/cosh(absorb_factor*REAL(j,8))
          bound_t(dim_t+1-j)=bound_t (dim_t+1-j)-1.D0/cosh(absorb_factor*REAL(j-1,8))
       ENDDO
    ENDIF
    ALLOCATE(p_t(dim_t),delta_rel(dim_t),op_t(dim_t),op_t_inv(dim_t),hfac(dim_t,0:4))
    ALLOCATE(DL(dim_r-2,dim_t_start(num_proc):dim_t_end(num_proc)),D(dim_r,dim_t_start(num_proc):dim_t_end(num_proc)), & 
    DU(dim_r-1,dim_t_start(num_proc):dim_t_end(num_proc)))
    SELECT CASE (switch_T)
    CASE(1)
       op_t=1.D0
       op_t_inv=1.D0
    CASE(2)
       DO j=1,dim_th
          op_t(j)=CMPLX(omega_uppe+k_t*REAL(j-1,8),0.D0,8)/omega
          op_t_inv(j)=CMPLX(omega/(omega_uppe+k_t*REAL(j-1,8)),0.D0,8)
          op_t(dim_th+j)=CMPLX(omega_uppe+k_t*REAL(j-dim_th-1,8),0.D0,8)/omega
          op_t_inv(dim_th+j)=CMPLX(omega/(omega_uppe+k_t*REAL(j-dim_th-1,8)),0.D0,8)
       ENDDO
    CASE(3)
       DO j=1,dim_th
          op_t(j)=rek0/omega**2*(omega_uppe+k_t*REAL(j-1,8))**2/komega(dim_th+j)
          op_t_inv(j)=rek0/komega(dim_th+j)
          op_t(dim_th+j)=rek0/omega**2*(omega_uppe+k_t*REAL(j-dim_th-1,8))**2/komega(j)
          op_t_inv(dim_th+j)=rek0/komega(j)
       ENDDO
    CASE(4)
       DO j=1,dim_th
          op_t(j)=rek0/omega**2*(omega_uppe+k_t*REAL(j-1,8))**2/komega(dim_th+j)
          op_t_inv(j)=rek0/komega(dim_th+j)
          op_t(dim_th+j)=rek0/omega**2*(omega_uppe+k_t*REAL(j-dim_th-1,8))**2/komega(j)
          op_t_inv(dim_th+j)=rek0/komega(j)
       ENDDO
    END SELECT
    IF (my_rank.EQ.0)  THEN
       OPEN(unit_logfile,FILE='OP_T.DAT',STATUS='UNKNOWN',RECL=2**8)
       DO j=1,dim_t
          WRITE(unit_logfile,*) j,REAL(op_t(j)),AIMAG(op_t(j)),REAL(op_t_inv(j)),AIMAG(op_t_inv(j))
       ENDDO
       CLOSE(unit_logfile)
    ENDIF
    SELECT CASE (switch_dKerr)
    CASE(1)
       c3i = c3
       c3d = 0.d0
    CASE(2)
       c3i=c3*(1.D0-xdk)
       c3d=c3*xdk     
       expt1=exp(-delta_t/tdk)
       expt2=-expt1+tdk/delta_t-tdk/delta_t*expt1
       expt3=1.D0-tdk/delta_t+tdk/delta_t*expt1
    CASE(3)
       c3i=c3*(1.D0-xdk)
       c3d=c3*xdk
       expt1=-(-exp(-delta_t/tdk)*cos(raman*delta_t)*tdk**3*raman**3*delta_t &
            -exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t*raman**2*tdk**2 &
            -exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t &
            -exp(-delta_t/tdk)*cos(raman*delta_t)*tdk*raman*delta_t)/(tdk*(raman**2*tdk**2+1)*delta_t*raman)
       expt2=-(-exp(-delta_t/tdk)*sin(raman*delta_t)*tdk**3*delta_t*raman**2 &
            -exp(-delta_t/tdk)*sin(raman*delta_t)*tdk*delta_t)/(tdk*(raman**2*tdk**2+1)*delta_t*raman)
       expt3=-(-exp(-delta_t/tdk)*sin(raman*delta_t)*raman**2*tdk**3 &
            +exp(-delta_t/tdk)*sin(raman*delta_t)*tdk-2*tdk**2*raman &
            +exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t*raman**2*tdk**2 &
            +exp(-delta_t/tdk)*cos(raman*delta_t)*tdk**3*raman**3*delta_t &
            +exp(-delta_t/tdk)*cos(raman*delta_t)*tdk*raman*delta_t &
            +2*exp(-delta_t/tdk)*cos(raman*delta_t)*tdk**2*raman &
            +exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t)/(tdk*(raman**2*tdk**2+1)*delta_t*raman)
       expt4=-(-tdk**3*delta_t*raman**3-2*exp(-delta_t/tdk)*cos(raman*delta_t)*tdk**2*raman &
            +exp(-delta_t/tdk)*sin(raman*delta_t)*raman**2*tdk**3-tdk*delta_t*raman &
            -exp(-delta_t/tdk)*sin(raman*delta_t)*tdk+2*tdk**2*raman)/(tdk*(raman**2*tdk**2+1)*delta_t*raman)
       expt1p=-(+exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t*raman**2*tdk**2+exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t)/(raman*delta_t*tdk**2)
       expt2p=-(-exp(-delta_t/tdk)*cos(raman*delta_t)*tdk**2*raman*delta_t+exp(-delta_t/tdk)*sin(raman*delta_t)*tdk*delta_t)/(raman*delta_t*tdk**2)
       expt3p=-(-exp(-delta_t/tdk)*cos(raman*delta_t)*tdk**2*raman &
            -exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t*raman**2*tdk**2+tdk**2*raman-exp(-delta_t/tdk)*sin(raman*delta_t)*tdk-exp(-delta_t/tdk)*sin(raman*delta_t)*delta_t)/(raman*delta_t*tdk**2)
       expt4p=-(+exp(-delta_t/tdk)*sin(raman*delta_t)*tdk-tdk**2*raman+exp(-delta_t/tdk)*cos(raman*delta_t)*tdk**2*raman)/(raman*delta_t*tdk**2)
    END SELECT
    CALL calc_propagator
    ALLOCATE(real_e(dim_t,dim_r/num_proc),imag_e(dim_t,dim_r/num_proc))
    CALL read_dset(group_id,'startfield_r',real_e,dim_t,dim_r,dim_t,dim_r/num_proc,0,(dim_r/num_proc)*my_rank)
    CALL read_dset(group_id,'startfield_i',imag_e,dim_t,dim_r,dim_t,dim_r/num_proc,0,(dim_r/num_proc)*my_rank)
    efield_factor = SQRT(critical_power*1.D9*3.D8*4.D0*3.1415D-7/(4.D0*3.1415D0*beam_waist**2*1.D-4*2.D0*n0_indice))*2.D0 ! normalization factor electric field V/m
    ALLOCATE(efield_osc(dim_t))
    DO j=1,dim_t
       efield_osc(j) = exp(CMPLX(0.D0,-omega_uppe*(tlo+REAL(j,8)*delta_t),8)) ! fast oscillating term exp(-i*omegauppe*t)
    ENDDO
    e = CMPLX(real_e,imag_e)
    DO k1=1, dim_t
      !DO k2=1, dim_r/num_proc
      !  real_part = 1/efield_factor*(cos(efield_osc(k1))*real_e(k1,k2) + sin(efield_osc(k1))*imag_e(k1,k2))
      !  imag_part = 1/efield_factor*(cos(efield_osc(k1))*imag_e(k1,k2) - sin(efield_osc(k1))*real_e(k1,k2))
      !  real_e(k1,k2) = real_part
      !  imag_e(k1,k2) = imag_part
      !ENDDO
      e(k1,:) = (1/efield_factor)*CONJG(efield_osc(k1))*e(k1,:)
    ENDDO
    !e = CMPLX(real_e,imag_e)

    !CALL read_dset(group_id,'startfield',e,dim_r,dim_t,dim_r/num_proc,dim_t,(dim_r/num_proc)*my_rank,0)
    CALL ask_for_size_1D(group_id, "indexes_group/r_vector", i_x_max)
    CALL ask_for_size_1D(group_id, "indexes_group/z_vector", i_z_max)
    ALLOCATE(xx(i_x_max),zz(i_z_max),Indice_norm(i_x_max,i_z_max))
    CALL read_dset(group_id, "indexes_group/indexes", Indice_norm, i_x_max, i_z_max)
    CALL read_dset(group_id, "indexes_group/r_vector", xx, i_x_max)
    CALL read_dset(group_id, "indexes_group/z_vector", zz, i_z_max)
    
    CALL h5gclose_f(group_id, error) ! all pre-processed inputs read
   
   
    ! DIRECT INPUTS (it bypasses the pre-processor) 
    ! normalisation factors
    CALL h5gopen_f(file_id, output_groupname, group_id, error) 
    CALL read_dset(group_id,'critdens',rhoc_cm3_phys)     
    CALL h5gclose_f(group_id, error)


    ! PREPARATION OF INPUTS
    ! compute normalization factors
 
    lambdanm = 6.634D-34*3.D17/photon_energy/4.359d-18 ! center wavelength in nm
    plasma_normalisation_factor_m3 = rhoc_cm3_phys/(4.0d0*PI**2 * (beam_waist**2) / ((lambdanm*1.0D-7)**2 )) ! in cm^(-3)
    plasma_normalisation_factor_m3 = 1.0D6 * plasma_normalisation_factor_m3 ! in m^(-3)
    tps = pulse_duration*1.D-15 ! pulse duration in s (tpfs*1.e-15 in octace files)
    w0m = beam_waist*1.D-2 ! beam width in m (w0cm*1e-2 in octave files)
    ! n0_indice : refractive index at center frequency (n0 in octave files)
    ! electric field: REAL(efield_factor*e(:,k)*efield_osc,4) : one temporal profile in GV/m    
    four_z_Rayleigh = 4.D0*3.1415D0*n0_indice/(lambdanm*1.D-9)*(beam_waist/100.D0)**2 ! 4 times the rayleigh length in m (normalization factor for z)
    Nz_points = CEILING(proplength/outlength)+1 ! expected number of hdf5 output along z (with safety)


    ! pre-ionisation
    CALL h5lexists_f(file_id,'pre_ionised',apply_pre_ionisation,error) ! it finds only if it's applied, the rest is fully encapsulated in the module        
    IF (apply_pre_ionisation) CALL init_pre_ionisation(file_id)

    
    ! CLOSE HDF5 interface (ionisation models will re-open again)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

    i_x_old=2
    i_z_old=2
    ALLOCATE(peakmax(rhodist),rhomax(rhodist),rhoabs_max(rhodist),energy(rhodist),z_buff(rhodist),energy_fil(rhodist),rhoO2max(rhodist),rhoN2max(rhodist),Tevmax(rhodist))
    ALLOCATE(rho(dim_r_start(num_proc):dim_r_end(num_proc)),fluence(dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(losses_ionization(dim_r_start(num_proc):dim_r_end(num_proc)),losses_plasma(dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(rhoabs(dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(e_2(dim_t),e_2KK(dim_t),e_2KKm2(dim_t))

   ! select the ionisation model used in the code
    SELECT CASE (switch_rho)
    CASE(1,2,6)
       CONTINUE
    CASE(8)
       CALL RESCALE_TABLE_EXT
    CASE(3)
       CALL INITIALISE_PPT('PPT')
       CALL FIND_INTENSITY_AREA('PPT')
       CALL FILL_TABLE('PPT')
    CASE(4)
       CALL INITIALISE_PPT('ADK')
       CALL FIND_INTENSITY_AREA('ADK')
       CALL FILL_TABLE('ADK')
    CASE(5)
       CALL INITIALISE_PPT('IRC')
       CALL FIND_INTENSITY_AREA('IRC')
       CALL FILL_TABLE('IRC')
    CASE(7)
       CALL INITIALISE_PPT('PPT')
       CALL FIND_INTENSITY_AREA('PPT')
       CALL FILL_TABLE('PPT')
       CALL INITIALIZE_PPT_N2(ionisation_potential_N2,residue_charge_N2,n0_indice,critical_density,atomic_density_N2,0.D0,angular_momentum_N2)
    END SELECT

    930 FORMAT (I3)

    print *, "subroutine initialize for proc ", ip, " done"
    RETURN
  END SUBROUTINE initialize

END MODULE first_step
