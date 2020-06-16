MODULE output
  USE fields
  USE parameters
  USE mpi_stuff
  USE long_step
  USE run_status
CONTAINS

  SUBROUTINE matlab_out
    USE fft
    IMPLICIT NONE

    INTEGER(4) j,k,l, k1, k2
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

    print *, "mtl, bbcast", my_rank

    CALL MPI_BCAST(filename,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    print *, "mtl, abcast", my_rank

    IF (filename.NE.iz) THEN

		print *, "mtl, if accessed", my_rank

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
	print *, "mtl, field written", my_rank
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
       print *, "mtl, plasma written", my_rank
       etemp=CSHIFT(e,dim_t/2-1,1)

	print *, "mtl, bfft", my_rank
        print *, my_rank, 'SIZE(etemp)', SIZE(etemp)

!       k1 = 1;
!       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
!	  k2 = 1;
!          DO  j=1,dim_th
!             etemp_test(k2,k1)=etemp(j,l)
!             k2 = k2 + 1;
!          ENDDO
!          k1 = k1 + 1;
!       ENDDO

       etemp_test = etemp;	

       CALL dfftw_execute(plan_spec)
	print *, "mtl, afft", my_rank


       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
          DO  j=1,dim_th
             help=etemp(j+dim_t/2,l)
             etemp(j+dim_t/2,l)=etemp(j,l)
             etemp(j,l)=help
          ENDDO
       ENDDO

	print *, "mtl, bfield", my_rank
       
       e_2(1:dim_t)=0.D0
       e_2KKm2(1:dim_t)=0.D0
       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
          r=REAL(l-1)*delta_r
          e_2(1:dim_t)=e_2(1:dim_t)+ABS(etemp(1:dim_t,l))**2*REAL(l-1,8)
          IF (rfil.GT.r) e_2KKm2(1:dim_t)=e_2KKm2(1:dim_t)+ABS(etemp(1:dim_t,l))**2*REAL(l-1,8)
       ENDDO

	print *, "mtl, breduce", my_rank

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

    print *, "mtl, finished", my_rank
    RETURN
  END SUBROUTINE  matlab_out

  SUBROUTINE  field_out
    USE ppt
    IMPLICIT  NONE

    INTEGER(4) j,k,i_x,i_z
    CHARACTER*10  iz,id,filename

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
