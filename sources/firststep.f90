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
    USE Complex_rotation 
    USE fft
    USE HDF5
    IMPLICIT NONE

    INTEGER(4)  j,k,help,i_x,i_z
    REAL(8) absorb_factor
    LOGICAL ext
    CHARACTER*10 filename,id
    CHARACTER(LEN=10), PARAMETER :: hdf5_input = "test.h5"  ! File name for the HDF5 input file

    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER        :: error

    HDF5write_count = 1
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
    CALL MPI_BCAST(filename,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    WRITE(ip,930) my_rank
    DO k=1,3
       IF (ip(k:k).EQ.' ') ip(k:k)='0'
    ENDDO

    CALL h5open_f(error) 
    CALL h5fopen_f (hdf5_input, H5F_ACC_RDONLY_F, file_id, error)
    CALL h5fclose_f(file_id, error)
    ! Close FORTRAN HDF5 interface.
    print *,"I'm here"
    CALL h5close_f(error)

    OPEN(unit_field,FILE=filename//'_'//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
    print *, filename, ip
    READ(unit_field) id,num_proc
    READ(unit_field) id,dim_t
    READ(unit_field) id,dim_r
    READ(unit_field) id,rek0
    READ(unit_field) id,rekp
    READ(unit_field) id,c3
    READ(unit_field) id,c5
    READ(unit_field) id,gamma1
    READ(unit_field) id,gamma2
    READ(unit_field) id,muk
    READ(unit_field) id,beta_inv_2KK
    READ(unit_field) id,KK
    READ(unit_field) id,rho0
    READ(unit_field) id,nu
    READ(unit_field) id,alpha
    READ(unit_field) id,alphaquad
    READ(unit_field) id,rhoat_inv
    READ(unit_field) id,xdk
    READ(unit_field) id,tdk
    READ(unit_field) id,raman
    READ(unit_field) id,omega
    ALLOCATE(komega(dim_t),komega_red(dim_t))
    READ(unit_field) id,komega(1:dim_t)
    READ(unit_field) id,NN
    READ(unit_field) id,eta1
    READ(unit_field) id,eta2
    READ(unit_field) id,lt
    READ(unit_field) id,lr
    READ(unit_field) id,proplength
    READ(unit_field) id,outlength
    READ(unit_field) id,delta_z
    READ(unit_field) id,z
    READ(unit_field) id,z_out
    READ(unit_field) id,rfil
    READ(unit_field) id,switch_rho
    READ(unit_field) id,switch_dKerr
    READ(unit_field) id,switch_T
    READ(unit_field) id,absorb
    READ(unit_field) id,increase
    READ(unit_field) id,decrease
    READ(unit_field) id,rhodist
    READ(unit_field) id,timelimit
    READ(unit_field) id,photon_energy
    READ(unit_field) id,pulse_duration
    READ(unit_field) id,critical_power
    READ(unit_field) id,beam_waist
    READ(unit_field) id,ionisation_potential
    READ(unit_field) id,residue_charge
    READ(unit_field) id,n0_indice
    READ(unit_field) id,critical_density
    READ(unit_field) id,atomic_density
    READ(unit_field) id,reduced_mass
    READ(unit_field) id,angular_momentum
    READ(unit_field) id,KKp
    READ(unit_field) id,beta_inv_2KKp
    READ(unit_field) id,mukp
    READ(unit_field) id,beta_inv_2
    READ(unit_field) id,mu
    READ(unit_field) id,KKpp
    READ(unit_field) id,beta_inv_2KKpp
    READ(unit_field) id,mukpp
    READ(unit_field) id,eti_ref
    READ(unit_field) id,exp_ref
    READ(unit_field) id,alpha1
    READ(unit_field) id,alpha2
    READ(unit_field) id,alphah
    READ(unit_field) id,rhosat
    READ(unit_field) id,finished
    READ(unit_field) id,omega_uppe
    READ(unit_field) id,gamma1e
    READ(unit_field) id,nuO2
    READ(unit_field) id,nuN2
    READ(unit_field) id,T_init_eV_phys
    READ(unit_field) id,nukB
    READ(unit_field) id,nucp
    READ(unit_field) id,nucO2
    READ(unit_field) id,nucN2
    READ(unit_field) id,rhoat_N2_inv
    READ(unit_field) id,ionisation_potential_N2
    READ(unit_field) id,residue_charge_N2
    READ(unit_field) id,atomic_density_N2
    READ(unit_field) id,angular_momentum_N2
    print *,"I loaded most of the stuff"
    READ(unit_field) id
    IF (id.NE.'startfield') THEN
       PRINT*, 'Error reading parameters'
       READ(5,*)
    ENDIF
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
    print *,"242"
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
    print *,"292"
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
    DO j=dim_r_start(num_proc),dim_r_end(num_proc)
       READ(unit_field) e(1:dim_t,j)
    ENDDO
    READ(unit_field) id
    IF (id.NE.'index') THEN
       PRINT*, 'Error reading index variations'
       READ(5,*)
    ENDIF
    READ(unit_field) i_x_max, i_z_max
    ALLOCATE(xx(i_x_max),zz(i_z_max),Indice_norm(i_x_max, i_z_max))
    READ(unit_field) (xx(i_x),i_x=1,i_x_max)
    DO i_z = 1, i_z_max
       READ(unit_field) zz(i_z)
       READ(unit_field) (Indice_norm(i_x,i_z),i_x=1,i_x_max)
    ENDDO  
    CLOSE(unit_field)
    i_x_old=2
    i_z_old=2
    ALLOCATE(peakmax(rhodist),rhomax(rhodist),rhoabs_max(rhodist),energy(rhodist),z_buff(rhodist),energy_fil(rhodist),rhoO2max(rhodist),rhoN2max(rhodist),Tevmax(rhodist))
    INQUIRE(FILE='FLUENCE_'//ip//'.DAT',EXIST=ext)
    IF (.NOT.ext) THEN
       OPEN(unit_rho,FILE='FLUENCE_'//ip//'.DAT',STATUS='NEW',FORM='UNFORMATTED')
       WRITE(unit_rho) num_proc,dim_r
       WRITE(unit_rho) (REAL(REAL(j-1,8)*delta_r,4),j=dim_r_start(num_proc),dim_r_end(num_proc))
       CLOSE(unit_rho)
    ENDIF
    print *,"355"
    INQUIRE(FILE='PLASMACHANNEL_'//ip//'.DAT',EXIST=ext)
    IF (.NOT.ext) THEN
       OPEN(unit_rho,FILE='PLASMACHANNEL_'//ip//'.DAT',STATUS='NEW',FORM='UNFORMATTED')
       WRITE(unit_rho) num_proc,dim_r
       WRITE(unit_rho) (REAL(REAL(j-1,8)*delta_r,4),j=dim_r_start(num_proc),dim_r_end(num_proc))
       CLOSE(unit_rho)
    ENDIF
    INQUIRE(FILE='LOSSES_IONIZATION_'//ip//'.DAT',EXIST=ext)
    IF (.NOT.ext) THEN
       OPEN(unit_rho,FILE='LOSSES_IONIZATION_'//ip//'.DAT',STATUS='NEW',FORM='UNFORMATTED')
       WRITE(unit_rho) num_proc,dim_r
       WRITE(unit_rho) (REAL(REAL(j-1,8)*delta_r,4),j=dim_r_start(num_proc),dim_r_end(num_proc))
       CLOSE(unit_rho)
    ENDIF
    INQUIRE(FILE='LOSSES_PLASMA_'//ip//'.DAT',EXIST=ext)
    IF (.NOT.ext) THEN
       OPEN(unit_rho,FILE='LOSSES_PLASMA_'//ip//'.DAT',STATUS='NEW',FORM='UNFORMATTED')
       WRITE(unit_rho) num_proc,dim_r
       WRITE(unit_rho) (REAL(REAL(j-1,8)*delta_r,4),j=dim_r_start(num_proc),dim_r_end(num_proc))
       CLOSE(unit_rho)
    ENDIF
    IF (my_rank.EQ.0)  THEN
       INQUIRE(FILE='ONAX_T.DAT',EXIST=ext)
       IF (.NOT.ext) THEN
          OPEN(unit_rho,FILE='ONAX_T.DAT',STATUS='NEW',FORM='UNFORMATTED')
          WRITE(unit_rho) dim_t
          WRITE(unit_rho) (REAL(tlo+REAL(j-1,8)*delta_t,4),j=1,dim_t)
          CLOSE(unit_rho)
       ENDIF
    ENDIF
    ALLOCATE(rho(dim_r_start(num_proc):dim_r_end(num_proc)),fluence(dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(losses_ionization(dim_r_start(num_proc):dim_r_end(num_proc)),losses_plasma(dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(rhoabs(dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(e_2(dim_t),e_2KK(dim_t),e_2KKm2(dim_t))
    print *,"391"
    SELECT CASE (switch_rho)
    CASE(1,2,6)
       CONTINUE
    CASE(8)
!       CALL INITIALISE_CPR
!       CALL FIND_INTENSITY_AREA_CPR
       CALL RESCALE_TABLE_CPR
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

! compute normalization factors
    tps = pulse_duration*1.D-15 ! pulse duration in s (tpfs*1.e-15 in octace files)
    w0m = beam_waist*1.D-2 ! beam width in m (w0cm*1e-2 in octave files)
    ! n0_indice : refractive index at center frequency (n0 in octave files)
    efield_factor = SQRT(critical_power*1.D9*3.D8*4.D0*3.1415D-7/(4.D0*3.1415D0*beam_waist**2*1.D-4*2.D0*n0_indice))*2.D0 ! normalization factor electric field V/m
    ALLOCATE(efield_osc(dim_t))
    DO j=1,dim_t
       efield_osc(j) = exp(CMPLX(0.D0,-omega_uppe*(tlo+REAL(j,8)*delta_t),8)) ! fast oscillating term exp(-i*omegauppe*t)
    ENDDO
    ! electric field: REAL(efield_factor*e(:,k)*efield_osc,4) : one temporal profile in GV/m
    lambdanm = 6.634D-34*3.D17/photon_energy/4.359d-18 ! center wavelength in nm
    four_z_Rayleigh = 4.D0*3.1415D0*n0_indice/(lambdanm*1.D-9)*(beam_waist/100.D0)**2 ! 4 times the rayleigh length in m (normalization factor for z)
    Nz_points = CEILING(proplength/outlength)+1 ! expected number of hdf5 output along z (with safety)
    print *, "end"
    RETURN
  END SUBROUTINE initialize

END MODULE first_step
