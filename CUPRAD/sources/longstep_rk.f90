!> @brief This module contains the main routine for one propagation step \ref propagation. Spectral operator splitting is used. The linear step in subroutine \ref mult_propagator is executed as two half-steps, in order to to second order accuracy of the splitting scheme. The nonlinear step in subroutine \ref mult_phase is executed as multiplication by exponential or second order Runge Kutta, depending whether higher-order SVEA corrections (T operators, switch_T=2) are taken into account or not (no T operator, switch_T=1)
!!
!! @author Jan Vábek
!! @author Stefan Skupin

! This is the main computational module stating explicite operations within every step
! 
! TABLE OF CONTENTS:
! "index_interpolation": inrepolate the refractive index if it varies in the medium
!
! "absorbation": absorb exiting radion to prevent reflections
!
! supplementary procedures of the evolution
!
! "mult_phase": the core of the propagation
!
! "propagation" drives the single-step calculation
!
!
! It provides an interface with other modules, so it was touched by
! the authors of most of them.


MODULE long_step
  USE constants
  USE parameters
  USE density_module
  REAL(8) rhompi,rho1,rho2,rhoth,rhotr,rhofh,rhoslg2,rhoav
CONTAINS

  SUBROUTINE index_interpolation(phase_index,r) ! refractive index lookup
    USE fields
    IMPLICIT NONE

    INTEGER(4) i_x,i_z,help
    REAL(8) phase_index,r

    IF ((z.GE.zz(1)).AND.(z.LE.zz(i_z_max)).AND.(r.LE.xx(i_x_max))) THEN
       IF (z.LT.zz(i_z_old)) THEN
          i_z_old=2
       ENDIF
       IF (r.LT.xx(i_x_old)) THEN
          i_x_old=2
       ENDIF
       DO i_z=i_z_old,i_z_max
          IF (z.LE.zz(i_z)) THEN
             help=i_z
             EXIT
          ENDIF
       ENDDO
       i_z_old=help   
       DO i_x=i_x_old,i_x_max
          IF (r.LE.xx(i_x)) THEN
             help=i_x
             EXIT
          ENDIF
       ENDDO
       i_x_old=help
       IF ((z.LT.zz(i_z_old-1)).OR.(z.GT.zz(i_z_old))) THEN
          print*, 'Error in z index interpolation', r, z
          READ(5)
       ENDIF
       IF ((r.LT.xx(i_x_old-1)).OR.(r.GT.xx(i_x_old))) THEN
          print*, 'Error in r index interpolation', r, z
          READ(5)
       ENDIF              
       phase_index=(zz(i_z_old)-z)*(xx(i_x_old)-r)*Indice_norm(i_x_old-1,i_z_old-1)
       phase_index=phase_index+(z-zz(i_z_old-1))*(xx(i_x_old)-r)*Indice_norm(i_x_old-1,i_z_old)
       phase_index=phase_index+(zz(i_z_old)-z)*(r-xx(i_x_old-1))*Indice_norm(i_x_old,i_z_old-1)
       phase_index=phase_index+(z-zz(i_z_old-1))*(r-xx(i_x_old-1))*Indice_norm(i_x_old,i_z_old)
       phase_index=phase_index/((zz(i_z_old)-zz(i_z_old-1))*(xx(i_x_old)-xx(i_x_old-1)))
    ELSE
       phase_index=0.D0
    ENDIF

    RETURN
  END SUBROUTINE index_interpolation

  
  SUBROUTINE absorbation ! apply absorption using cosh in time on the boundary (see firstep.f90) 
    USE fields
    USE mpi_stuff
    IMPLICIT NONE

    INTEGER(4) j,l

    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
       DO j=1,dim_t
          e(j,l)=e(j,l)*bound_t(j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE absorbation


  SUBROUTINE tridag(dim,DLn,D,DL,DU,B)
    IMPLICIT NONE

    INTEGER(4) dim
    COMPLEX(8) DLn
    COMPLEX(8) D(dim),DL(2:dim-1),DU(dim-1),B(dim)
    INTEGER(4) j
    COMPLEX(8) bet,gam(2:dim)

    bet=D(1)
    B(1)=B(1)/bet
    DO j=2,dim-1
       gam(j)=DU(j-1)/bet
       bet=D(j)-DL(j)*gam(j)
       B(j)=(B(j)-DL(j)*B(j-1))/bet
    ENDDO
    gam(dim)=DU(dim-1)/bet
    bet=D(dim)-DLn*gam(dim)
    B(dim)=(B(dim)-DLn*B(dim-1))/bet
    DO j=dim-1,1,-1
       B(j)=B(j)-gam(j+1)*B(j+1)
    ENDDO

    RETURN
  END SUBROUTINE tridag

  SUBROUTINE cn(B,k)
    USE fields
    USE parameters
    IMPLICIT NONE

    INTEGER(4) k
    COMPLEX(8) B(dim_r)
    INTEGER(4) j
    COMPLEX(8) DLn,help_1,help_2

    help_1=B(1)
    B(1)=(CMPLX(1.D0,0.D0,8)+CMPLX(0.D0,-2.D0,8)*delta_rel(1,k))*B(1)+CMPLX(0.D0,2.D0,8)*delta_rel(1,k)*B(2)
    DO j=2,dim_r-1
       help_2=help_1
       help_1=B(j)
       B(j)=-DL(j-1,k)*help_2+(CMPLX(1.D0,0.D0,8)+CMPLX(0.D0,-1.D0,8)*delta_rel(j,k))*B(j)-DU(j,k)*B(j+1)
    ENDDO
    B(dim_r)=CMPLX(0.D0,0.D0,8)*B(dim_r)

    IF (help_2.NE.(0.D0,0.D0)) THEN
       DLn=-help_1/help_2
       IF (AIMAG(DLn).GT.0.D0) DLn=CMPLX(-1.D0,0.D0,8)
    ELSE
       DLn=(-1.D0,0.D0)
    ENDIF

    CALL tridag(dim_r,DLn,D(:,k),DL(:,k),DU(:,k),B)

    RETURN
  END SUBROUTINE cn
  
  !> @brief Linear propagation step using Cranck-Nicholson for diffraction (\ref cn) and mutltiplication by exponential for linear dispersion.

  SUBROUTINE mult_propagator ! it does the propagation in z (linear dispersion - linear part of the propagator)
    USE fields
    USE mpi_stuff
    IMPLICIT NONE

    INTEGER(4) k,l

    DO k=dim_t_start(num_proc),dim_t_end(num_proc)
       CALL cn(efft(1:dim_r,k),k)
       DO l=1,dim_r
         efft(l,k)=efft(l,k)*p_t(l,k) ! p_t - the dispersion "p_t =  exp(i*(k(omega)*z-k' * omega*z)" (co-moving frame taken)
                                              ! k(omega) = (omega/c)*sqrt(eps(omega))
                                              ! k' evaluated at omega0 defines our co-moving frame, in fact exp(i*(k(omega)*z-(omega*z/vg) )
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mult_propagator

  SUBROUTINE calc_absorption(rhoabs, mediumabs, eti, etip1)    ! N-photon absorption (to be removed)
    REAL(8) :: rhoabs, mediumabs, eti, etip1
    REAL(8) :: intF

    IF ( eta1.NE.0.D0 ) THEN
       intF = -eta2 * 0.5d0 *( eti**Nn + etip1**NN) * delta_t  
       rhoabs = 1.d0 - (1.D0 - rhoabs) * exp(intF) ! rate eq. for excited molecules
       IF (rhoabs.LT.0.D0) rhoabs = 0.D0
       mediumabs = eta1 * etip1**(NN-1)*(1.D0-rhoabs)
    ENDIF

    RETURN
  END SUBROUTINE calc_absorption

  SUBROUTINE calc_rho(rho,mpa,eti,etip1,l)
    USE PPT
    USE External_ionisation_table
    IMPLICIT NONE

    INTEGER(4) :: l
    REAL(8) :: rho,mpa
    REAL(8) :: eti,etip1
    REAL(8) intF,var1,rhosave

    rhosave=rho
    SELECT CASE (switch_rho)
    CASE(8)
       CALL interpolate_ext(var1,mpa,etip1)
       intF=(nu*0.5d0*(etip1+eti)*density_mod(l)-rhoat_inv*var1-alpha)*delta_t
       rho=rho*exp(intF)+var1*density_mod(l)*delta_t
       mpa=mpa*(density_mod(l)-rhosave*rhoat_inv)
    CASE(1)
       rho=rho+(nu*rho*eti*density_mod(l)+beta_inv_2KK*eti**KK*(density_mod(l)-rho*rhoat_inv)-alpha*rho)*delta_t
       mpa=muk*etip1**(KK-1)*(density_mod(l)-rhosave*rhoat_inv)
    CASE(2)
       var1=0.5d0*beta_inv_2KK*(etip1**KK+eti**KK)
       intF=(nu*0.5d0*(etip1+eti)*density_mod(l)-rhoat_inv*var1-alpha)*delta_t
       rho=rho*exp(intF)+var1*density_mod(l)*delta_t
       mpa=muk*etip1**(KK-1)*(density_mod(l)-rhosave*rhoat_inv)
   !  CASE(3,4)
    CASE(3)
       CALL interpolate_ppt(var1,mpa,etip1)
       intF=(nu*0.5d0*(etip1+eti)*density_mod(l)-rhoat_inv*var1-alpha)*delta_t
       rho=rho*exp(intF)+var1*density_mod(l)*delta_t
       mpa=mpa*(density_mod(l)-rhosave*rhoat_inv)
    END SELECT
    !print *, 'alpha4', alphaquad
    rho=rho-alphaquad*rhosave**2*delta_t
    IF (rho.LT.0.D0) rho=0.D0

    IF (density_mod(l) .GT. 0.D0) THEN
      IF ( (rho*rhoat_inv/density_mod(l)) .GT.1.D0) rho= density_mod(l)/rhoat_inv
    ENDIF

    IF (rhompi.LT.0.D0) rhompi=0.D0
    IF (rho1.LT.0.D0) rho1=0.D0
    IF (rho2.LT.0.D0) rho2=0.D0
    IF (rhoth.LT.0.D0) rhoth=0.D0
    IF (rhotr.LT.0.D0) rhotr=0.D0
    IF (rhofh.LT.0.D0) rhofh=0.D0
    IF (rhoslg2.LT.0.D0) rhoslg2=0.D0
    IF (rhoav.LT.0.D0) rhoav=0.D0

    RETURN
  END SUBROUTINE calc_rho

  SUBROUTINE calc_delkerr(delkerr,delkerrp, eti,etip1)
    IMPLICIT NONE

    REAL(8)  :: delkerr, delkerrp
    REAL(8)   :: eti
    REAL(8)   :: etip1

    SELECT CASE (switch_dKerr)
    CASE(1)
       CONTINUE
    CASE(2)
       delkerr=expt1*delkerr+expt2*eti+expt3*etip1
    CASE(3)
       delkerr=expt1*delkerr+expt2*delkerrp+expt3*eti+expt4*etip1
       delkerrp=expt1p*delkerr+expt2p*delkerrp+expt3p*eti+expt4p*etip1
    END SELECT

    RETURN
  END SUBROUTINE calc_delkerr

  !================================== 
  ! The main computational subroutine
  !==================================

 !> @brief Nonlinear propagation step. For standard SVEA (no T operators, switch_T=1) all nonlinear terms can be treated by exponential factor. The exponential factor is evaluated for an intermediate electric field at half-step. For higher-order SVEA (with T operators, switch_T=2), a second order Runge Kutta scheme is used, taking into account intermediate field at half-step. Note that for switch_T=2 the final propagation nonlinear propagation step is performed in \ref fft::fft_forward_inplace.
  
  SUBROUTINE mult_phase
    USE fields
    USE fft
    USE mpi_stuff
    USE HDF5
    USE HDF5_helper
    USE longstep_vars
    USE linked_list
    USE ll_data
    USE ppt
    USE normalization
    USE pre_ionised
    USE h5namelist
    IMPLICIT NONE

    INTEGER(4) j,l,k,k1
    REAL(8)  phase,maxphase_part,peakmax_part,energy_part,energy_fil_part,rhomax_part,delkerr,delkerrp,rhotemp,mpa,r,phase_p,phase_j,losses_j,phase_index
    REAL(8) mediumabs , rhoabs_max_part, rhoabstemp

    ! For storing in HDF5
    INTEGER(HID_T)    :: file_id       ! File identifier 
    INTEGER(HID_T)    :: group_id      ! Group identifier 
    INTEGER           :: error         ! hdferr
    LOGICAL           :: group_status  ! boolean - does the group exists?
   !  CHARACTER(LEN=15) :: h5_filename="results.h5" ! hdf5 file name
    
    ! Names of dsets and groups
    CHARACTER(*), PARAMETER :: groupname = longstep_grpname     
    CHARACTER(*), PARAMETER :: rhoabs_max_dset_name=longstep_grpname//"/rhoexcmax_max"
    CHARACTER(*), PARAMETER :: peakmax_dset_name=longstep_grpname//"/peakmax"
    CHARACTER(*), PARAMETER :: energy_dset_name=longstep_grpname//"/energy"
    CHARACTER(*), PARAMETER :: energy_fil_dset_name=longstep_grpname//"/energy_fil"
    CHARACTER(*), PARAMETER :: rhomax_dset_name=longstep_grpname//"/rhomax"
    CHARACTER(*), PARAMETER :: powmax_dset_name=longstep_grpname//"/powmax"
    CHARACTER(*), PARAMETER :: z_buff_dset_name=longstep_grpname//"/z_buff"
    CHARACTER(*), PARAMETER :: every_rhodist_z_dset_name=longstep_grpname//"/zgrid_analyses2"
    CHARACTER(*), PARAMETER :: onax_t_dset_name=longstep_grpname//"/Efied_onaxis"
    
    ! Arrays needed for temporary storing of dataset content and their writing
    INTEGER(HSIZE_T), DIMENSION(1) :: new_dims, memspace_dims, offset, hyperslab_size
    REAL,DIMENSION(1,2) :: powmax_data
    REAL,ALLOCATABLE  :: onax_t_data(:,:)

    ! Variables needed for linked list buffering
    REAL, ALLOCATABLE, TARGET  :: fluence_data(:)
    REAL, ALLOCATABLE, TARGET  :: plasma_channel_data(:)
    REAL, ALLOCATABLE, TARGET  :: losses_plasma_data(:)
    REAL, ALLOCATABLE, TARGET  :: losses_ionization_data(:)
   
    !=====================
    ! INITIALISATION PHASE
    !=====================
    
    ! Prepare the physical effects (ionisation, diffraction, ...)

    rho=0.D0
    rhoabs = 0.D0 
    rhoabs_max_part = 0.D0
    losses_plasma=0.D0
    losses_ionization=0.D0
    peakmax_part=0.D0
    rhomax_part=0.D0
    energy_part=0.D0
    energy_fil_part=0.D0
    maxphase_part=0.D0

    DO l=dim_r_start(num_proc),dim_r_end(num_proc) !
       r=REAL(l-1)*delta_r
       e_2=ABS(e(1:dim_t,l))
       peakmax_part=MAX(peakmax_part,MAXVAL(e_2))
       e_2=e_2**2
       fluence(l)=SUM(e_2)                                                    ! get fluence
       energy_part=energy_part+fluence(l)*REAL(l-1,8)                         ! get energy
       IF (rfil.GT.r) energy_fil_part=energy_fil_part+fluence(l)*REAL(l-1,8)  ! fluence(l)*REAL(l-1,8) for the radial Jacobian
       delkerr=0.D0
       delkerrp=0.d0
      IF (apply_pre_ionisation) THEN                                          ! initial electron density
         rhotemp = initial_electron_density(r,z,l,dim_r_start(num_proc))
      ELSE
         rhotemp = 0.D0
      ENDIF
      !print *, 'rhotmp', rhotemp, 'myrank', my_rank
      !print *, 'rhoatm', 1.d0/rhoat_inv, 'myrank', my_rank
      !stop
       rhompi=0.D0
       rho1=0.D0
       rho2=0.D0
       rhoth=0.D0
       rhotr=0.D0
       rhofh=0.D0
       rhoslg2=0.D0
       rhoav=0.D0


       rhoabstemp=0.D0
       mediumabs=0.D0

       ! Physical effects
       CALL index_interpolation(phase_index,r) ! refractive index
       phase_index=phase_index*delta_zh
       CALL calc_rho(rhotemp,mpa,0.D0,e_2(1),l) ! compute ionisation
       CALL calc_absorption(rhoabstemp, mediumabs, 0.D0,e_2(1)) ! N-photon absorption (to be removed)

       DO j=1,dim_t
          phase_p=(c3i*e_2(j)+c3d*delkerr-c5*e_2(j)**2)*((density_mod(l)-rhotemp*rhoat_inv)+ions_Kerr_ratio*rhotemp*rhoat_inv)*delta_zh ! phase_p is polarisation
          ! polarisation = i*phase_p
          ! c3i: instanataneous Kerr
          ! c3d: delayed Kerr (Raman)
          phase_j=-gamma2*rhotemp*delta_zh
          losses_j=-gamma1*density_mod(l)*rhotemp*delta_zh ! ~ losses_plasma(l)
          rho(l)=MAX(rho(l),rhotemp)
          rhoabs(l) = MAX(rhoabs(l),rhoabstemp)
          losses_plasma(l)=losses_plasma(l)+2.D0*e_2(j)*gamma1*density_mod(l)*rhotemp  ! losses for diagnostic
          losses_ionization(l)=losses_ionization(l)+2.D0*e_2(j)*mpa     ! losses for diagnostic (from ionization rate eq.)
          phase=phase_p+phase_j+phase_index
          maxphase_part=MAX(maxphase_part,ABS(phase))
          
! till here

          SELECT CASE (switch_T) ! decision over various propagators
          CASE(1)
             etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase,8))
          CASE(2)
             ptemp(j,l)=e(j,l)*CMPLX(0.D0,phase_p)
             jtemp(j,l)=e(j,l)*CMPLX(0.D0,phase_j)
             etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase_index,8)) ! applying the losses, T-operator does not affect
!           CASE(3)
!              ptemp(j,l)=e(j,l)*CMPLX(0.D0,phase_p) + (hfac(j,4)*CONJG(hfac(j,0)*e(j,l))*CMPLX(0.D0,phase_p) &
!                   + (hfac(j,1)*(hfac(j,0)*e(j,l))**3 + hfac(j,2)*(hfac(j,0)*e(j,l))**3*e_2(j)) &
!                   + (hfac(j,3)*(hfac(j,0)*e(j,l))**5))*((density_mod(l)-rhotemp*rhoat_inv)+ions_Kerr_ratio*rhotemp*rhoat_inv)
!             ! terms comes from the expansion (E+E*)^3
!             ! e(j,l)*CMPLX(0.D0,phase_p) - the same as above (usual Kerr)
!             ! hfac(j,4)*CONJG(hfac(j,0)*e(j,l))*CMPLX(0.D0,phase_p) 3E*^2E (-omega process)
!             ! (hfac(j,1)*(hfac(j,0)*e(j,l))**3 - E^3 (3omega)
!             ! hfac(j,2)*(hfac(j,0)*e(j,l))**3*e_2(j)) - 5th-order process
!             ! (hfac(j,3)*(hfac(j,0)*e(j,l))**5) - 5th-order process
!              jtemp(j,l)=e(j,l)*CMPLX(0.D0,phase_j)                                                  ! same treatment of the losses
!              etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase_index,8)) 
!           CASE(4)                                                                                   ! to be removed
!              ptemp(j,l)=e(j,l)*CMPLX(0.D0,phase_p)
!              jtemp(j,l)=e(j,l)*CMPLX(0.D0,phase_j)
!              etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase_index,8))
          END SELECT
          IF (j.NE.dim_t) THEN
             CALL calc_rho(rhotemp,mpa,e_2(j),e_2(j+1),l) ! update ionization
             CALL calc_absorption(rhoabstemp, mediumabs, e_2(j), e_2(j+1))
             CALL calc_delkerr(delkerr,delkerrp,e_2(j),e_2(j+1)) ! 2nd ored eq.; it requires derivative
          ENDIF
       ENDDO
    ENDDO


    !==================================
    ! THE APPLICATION OF THE PROPAGATOR (physics in omega-domain)
    !==================================

    ! This is the core of the propagation. As the code works in the omega-domain:
    !  1) it performs FFT initialised in fft.f90:fft_init
    !  2) it applies the operator in the omega-domain
    !  3) it goes back to time domain

    SELECT CASE (switch_T) ! decision over various propagators
    CASE(1) 
       continue ! non-lin step is done becasue the full exponential was applied in the previous switch
    CASE(2) 
       CALL dfftw_execute(plan_forward_erk)
       CALL dfftw_execute(plan_p)
       CALL dfftw_execute(plan_j)
       DO k=dim_r_start(num_proc),dim_r_end(num_proc)
          DO l=1,dim_t
             etemp(l,k)=etemp(l,k)+op_t(l,k)*ptemp(l,k)+op_t_inv(l,k)*jtemp(l,k) ! Euler step in omega-space (intermediate step of Runge-Kutta)
          ENDDO
       ENDDO
       CALL dfftw_execute(plan_backward_erk) ! field to time domain
       etemp=diminv*etemp ! fft - normalization
!     CASE(3) ! possibly merge with the previous case (CASE(2,3))
!        CALL dfftw_execute(plan_forward_erk) 
!        CALL dfftw_execute(plan_p)
!        CALL dfftw_execute(plan_j)
!        DO k=dim_r_start(num_proc),dim_r_end(num_proc)
!           DO l=dim_th+1,dim_t ! different range - this is the cutting of high freq. (low-pass filter)
!              etemp(l,k)=etemp(l,k)+op_t(l,k)*ptemp(l,k)+op_t_inv(l,k)*jtemp(l,k)
!           ENDDO
!        ENDDO
!        CALL dfftw_execute(plan_backward_erk)
!        etemp=diminv*etemp
!     CASE(4)
!        CALL dfftw_execute(plan_forward_erk)
!        CALL dfftw_execute(plan_p)
!        CALL dfftw_execute(plan_j)
!        DO k=dim_r_start(num_proc),dim_r_end(num_proc)
!           DO l=1,dim_t
!              etemp(l,k)=etemp(l,k)+op_t(l,k)*ptemp(l,k)+op_t_inv(l,k)*jtemp(l,k)
!           ENDDO
!        ENDDO
!        CALL dfftw_execute(plan_backward_erk)
!        etemp=diminv*etemp
    END SELECT


    !=======================
    ! PHYSICS IN TIME DOMAIN
    !=======================

    ! The second step of the Runge-Kutta (fft (correxponding to the previous case) is done in the main cuprad routine)

    ! Physical effects in time domain and on-the-fly analyses are computed here:

    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
       r=REAL(l-1)*delta_r
       e_2=ABS(etemp(1:dim_t,l))
       e_2=e_2**2
       delkerr=0.D0
       delkerrp=0.d0
      IF (apply_pre_ionisation) THEN
         rhotemp = initial_electron_density(r,z,l,dim_r_start(num_proc))
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

       rhoabstemp=0.D0
       mediumabs=0.D0
       CALL index_interpolation(phase_index,r)
       phase_index=phase_index*delta_z
       CALL calc_rho(rhotemp,mpa,0.D0,e_2(1),l)
       CALL calc_absorption(rhoabstemp, mediumabs, 0.D0,e_2(1))
       DO j=1,dim_t
          phase_p=(c3i*e_2(j)+c3d*delkerr-c5*e_2(j)**2)*((density_mod(l)-rhotemp*rhoat_inv)+ions_Kerr_ratio*rhotemp*rhoat_inv)*delta_z
          !phase_p=(c3i*e_2(j)+c3d*delkerr-c5*e_2(j)**2)*delta_z
          phase_j=-gamma2*rhotemp*delta_z
          losses_j=-gamma1*density_mod(l)*rhotemp*delta_z
          phase=phase_p+phase_j+phase_index
          SELECT CASE (switch_T)
          CASE(1)
             e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase,8))
          CASE(2)
             ptemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_p)
             jtemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_j)
             e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase_index,8))
!           CASE(3)
!              ptemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_p) + (hfac(j,4)*CONJG(hfac(j,0)*etemp(j,l))*CMPLX(0.D0,phase_p) &
!                   + 2.D0*(hfac(j,1)*(hfac(j,0)*etemp(j,l))**3 + hfac(j,2)*(hfac(j,0)*etemp(j,l))**3*e_2(j)) &
!                   + 2.D0*(hfac(j,3)*(hfac(j,0)*etemp(j,l))**5))*((density_mod(l)-rhotemp*rhoat_inv)+ions_Kerr_ratio*rhotemp*rhoat_inv)
!              jtemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_j)
!              e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase_index,8))
!           CASE(4)
!              ptemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_p)
!              jtemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_j)
!              e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase_index,8))
          END SELECT
          IF (j.NE.dim_t) THEN
             CALL calc_rho(rhotemp,mpa,e_2(j),e_2(j+1),l)
             CALL calc_absorption(rhoabstemp, mediumabs, e_2(j), e_2(j+1))
             CALL calc_delkerr(delkerr,delkerrp,e_2(j),e_2(j+1))
          ENDIF
       ENDDO
    ENDDO



    !=============================================!
    !                                             !
    !                SAVE OUTPUTS                 !
    !                                             !
    !=============================================!

    ! The following procedures treats only the storing of some characteristics of the pulse.
    ! There are no computationally relevant procedures till the end of the 'mult_phase' procedure.



    ! The aggregated data are distributed to the writer (proc 0).
    ! Some data are dtored directly in HDF5-archive and some are 
    ! buffered in the linked list.
    
    CALL MPI_REDUCE (maxphase_part,maxphase,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(maxphase,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    count=count+1
    rhomax_part=MAXVAL(rho(dim_r_start(num_proc):dim_r_end(num_proc)))
    rhoabs_max_part = MAXVAL(rhoabs(dim_r_start(num_proc):dim_r_end(num_proc)))
    CALL MPI_REDUCE(peakmax_part,peakmax(count),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(rhomax_part,rhomax(count),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(rhoabs_max_part, rhoabs_max(count),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr) 
    CALL MPI_REDUCE(energy_part,energy(count),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(energy_fil_part,energy_fil(count),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    IF(my_rank.EQ.0) z_buff(count)=four_z_Rayleigh*z
    ! If any writting is going to happen, allocate arrays with flag target for addition to linked list buffers
    IF (count.GE.rhodist .OR. z.LE.delta_z) THEN
      ALLOCATE(onax_t_data(1,1:dim_t))
      ALLOCATE(fluence_data(1:dim_r/num_proc))
      ALLOCATE(plasma_channel_data(1:dim_r/num_proc))
      ALLOCATE(losses_plasma_data(1:dim_r/num_proc))
      ALLOCATE(losses_ionization_data(1:dim_r/num_proc))
    ENDIF

    IF(count.GE.rhodist) THEN
       e_2=0.D0
       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
          e_2=e_2+ABS(e(1:dim_t,l))**2*REAL(l-1,8)
       ENDDO
       CALL MPI_REDUCE(e_2(1:dim_t),e_2KK(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF(my_rank.EQ.0) THEN
          CALL h5open_f(error) ! Prepare the HDF5 writing
          CALL h5fopen_f(main_h5_fname, H5F_ACC_RDWR_F, file_id, error ) ! Open the HDF5 file
          ! Create long_step group if it does not exist yet
          CALL h5lexists_f(file_id, longstep_grpname, group_status, error)
          IF ( group_status .EQV. .FALSE. ) THEN
            CALL h5gcreate_f(file_id, longstep_grpname, group_id, error) 
            CALL h5gclose_f(group_id, error)
          ENDIF

          ! if z was not lower or equal than delta_z and rhodist was reached for the first time
          IF (longstep_write_count.EQ.0) THEN
            ! z_buff dataset creation and writting, originally the z_buff was paired with each of the values from following
            ! datasets, here it has it's own dataset, the dimensions are equivalent
            CALL create_1D_dset_unlimited(file_id, z_buff_dset_name, REAL(z_buff(1:rhodist),4), rhodist)

            ! rhoabs_max dataset creation and writting
            CALL create_1D_dset_unlimited(file_id, rhoabs_max_dset_name, REAL(rhoabs_max(1:rhodist),4), rhodist) 
            ! the following normalisation code is commented out, because of unsolved variable linking
            ! rhoabs_max(1:rhodist)=rhoc_cm3_phys/(4.D0*PI**2*beam_waist**2/(wavelength*1.D-7)**2)*rhoabs_max(1:rhodist)
            CALL create_1D_dset_unlimited(file_id, TRIM(rhoabs_max_dset_name)//"_normalised", REAL(rhoabs_max(1:rhodist),4), rhodist)
           
            ! peakmax creation of the dataset, writting to it and doing the same after the normalisation
            CALL create_1D_dset_unlimited(file_id, peakmax_dset_name, REAL(peakmax(1:rhodist),4), rhodist)
            peakmax(1:rhodist) = critical_power/(1000.d0*4.d0*PI*beam_waist**2)*peakmax(1:rhodist)*peakmax(1:rhodist)
            CALL create_1D_dset_unlimited(file_id, TRIM(peakmax_dset_name)//"_normalised", REAL(peakmax(1:rhodist),4), rhodist)
            
            ! rhomax dataset creation and writting to it
            CALL create_1D_dset_unlimited(file_id, rhomax_dset_name, REAL(rhomax(1:rhodist),4), rhodist)
            
            ! energy dataset creation and writting to it, following the same process with normalised data
            ! CALL create_1D_dset_unlimited(file_id, energy_dset_name, &
            !   REAL(2.D0*PI*energy(1:rhodist)*delta_t*delta_r**2,4), rhodist) 

            energy(1:rhodist)=1.D-6*critical_power*1.D9/(4.D0*PI)*pulse_duration*1.D-15*1.D6*energy(1:rhodist)
            CALL create_1D_dset_unlimited(file_id, energy_dset_name, &
              REAL(2.D0*PI*energy(1:rhodist)*delta_t*delta_r**2,4), rhodist)

            CALL h5_add_units_1D(file_id, energy_dset_name, '[J]')
            
            ! energy_fil dataset creation and writting to it, following the same process with normalised data
            ! CALL create_1D_dset_unlimited(file_id, energy_fil_dset_name, &
            !   REAL(2.D0*PI*energy_fil(1:rhodist)*delta_t*delta_r**2,4), rhodist)

            energy_fil(1:rhodist)=1.D-6*critical_power*1.D9/(4.D0*PI)*pulse_duration*1.D-15*1.D6*energy_fil(1:rhodist)
            CALL create_1D_dset_unlimited(file_id, energy_fil_dset_name, &
              REAL(2.D0*PI*energy_fil(1:rhodist)*delta_t*delta_r**2,4), rhodist)
         
            CALL h5_add_units_1D(file_id, energy_fil_dset_name, '[J]')
           

            ! write max power with corresponding z to a variable and prepare for writting to a dataset
            !powmax_data(1,1) = REAL(z,4)
            !powmax_data(1,2) = REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)
            ! write on axis in time tensor to a variable and prepare for writting to a dataset
            !z_data(1) = REAL(z,4)
            !onax_t_data(1,:) = REAL(ABS(e(1:dim_t,1)),4)
            DO k1=1, dim_t
              onax_t_data(1,k1) = REAL( efield_factor * REAL( efield_osc(k1)*e(k1,1) ) , 4 )
              !fields_array(1,k1,k2) = REAL( REAL( (efield_factor*efield_osc(k2)*e(k2,r_offset+k1)) ) , 4 ) ! SINGLE PRECISION, corresponding H5T_NATIVE_REAL (REAL(.,8) corresponds to H5T_NATIVE_DOUBLE)
              ! e(t,r)
            ENDDO
            ! create the datasets if they do not exist yet 
            IF ( dset_write_count .EQ. 0 ) THEN
              CALL create_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)/), 1)
              CALL h5_add_units_1D(file_id, powmax_dset_name, '[C.U.]')
              CALL create_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, dim_t)
              CALL h5_add_units_1D(file_id, onax_t_dset_name, '[V/m]?')
              CALL create_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh*z,4)/), 1)
              CALL h5_add_units_1D(file_id, every_rhodist_z_dset_name, '[m]')
            ! extend datasets if they do exist
            ELSE
              CALL extend_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)/), &
              new_dims=(/int(dset_write_count + 1,HSIZE_T)/), memspace_dims=(/int(1,HSIZE_T)/), &
              offset = (/int(dset_write_count, HSIZE_T)/), hyperslab_size = (/int(1,HSIZE_T)/))
              CALL extend_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, & 
                new_dims = (/int(dset_write_count + 1, HSIZE_T), int(dim_t, HSIZE_T)/), & 
                memspace_dims = (/int(1,HSIZE_T), int(dim_t, HSIZE_T)/), & 
                offset = (/int(dset_write_count,HSIZE_T),int(0,HSIZE_T)/), & 
                hyperslab_size = (/int(1,HSIZE_T), int(dim_t, HSIZE_T)/))
              CALL extend_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh*z,4)/), &
                new_dims=(/int(dset_write_count + 1, HSIZE_T)/), memspace_dims=(/int(1,HSIZE_T)/), &
                offset = (/int(dset_write_count, HSIZE_T)/), hyperslab_size = (/int(1,HSIZE_T)/))
            ENDIF

            ! increase write counters
            dset_write_count = dset_write_count + 1
            longstep_write_count = longstep_write_count + 1
            ! store the original rhodist size to prepare for the last iteration which forces the data storage by `rhodist = count`
            original_rhodist = rhodist
          ! if the rhodist was already reached or z was lower of equal to delta_z
          ELSE
            ! rhodist is used for the extensions of the unlimited datasets
            IF (rhodist.NE.original_rhodist) THEN
              new_dims = (/int(longstep_write_count*original_rhodist+rhodist, HSIZE_T)/)
            ELSE
              new_dims = (/int((longstep_write_count+1)*rhodist, HSIZE_T)/)
            ENDIF
            ! dataspace size
            memspace_dims = (/int(rhodist,HSIZE_T)/)
            ! offset (the size of data present in the datasets)
            offset = (/int(longstep_write_count*original_rhodist,HSIZE_T)/)
            ! size of the current extension to the datasets
            hyperslab_size = (/int(rhodist,HSIZE_T)/)
            ! prevent writting of datasets of size 0
            IF (count .EQ. 0) THEN
              print *,"Last write skipped (longstep_rk.f90)"
            ELSE
              ! extend z_buff dataset
              CALL extend_1D_dset_unlimited(file_id, z_buff_dset_name, & 
                REAL(z_buff(1:count),4), new_dims, &
                memspace_dims, offset, hyperslab_size)
              
              ! extend rhoabs_max dataset
              CALL extend_1D_dset_unlimited(file_id, rhoabs_max_dset_name, & 
                REAL(rhoabs_max(1:count),4), new_dims, &
                memspace_dims, offset, hyperslab_size)
              
              ! extend peakmax dataset and the normalised peakmax dataset
              CALL extend_1D_dset_unlimited(file_id, peakmax_dset_name, & 
                REAL(peakmax(1:count),4), new_dims, &
                memspace_dims, offset, hyperslab_size)
              peakmax(1:rhodist) = critical_power/(1000.D0*4.D0*PI*beam_waist**2)*peakmax(1:rhodist)*peakmax(1:rhodist)
              CALL extend_1D_dset_unlimited(file_id, TRIM(peakmax_dset_name)//"_normalised", & 
                REAL(peakmax(1:count),4), new_dims, &
                memspace_dims, offset, hyperslab_size)
              
              ! extend rhomax dataset
              CALL extend_1D_dset_unlimited(file_id, rhomax_dset_name, & 
                REAL(rhomax(1:count),4), new_dims, &
                memspace_dims, offset, hyperslab_size)

              ! extend energy dataset and the normalised energy dataset
            !   CALL extend_1D_dset_unlimited(file_id, energy_dset_name, & 
            !     REAL(2.D0*PI*energy(1:count)*delta_t*delta_r**2,4), new_dims, &
            !     memspace_dims, offset, hyperslab_size)

              energy(1:rhodist)=1.D-6*critical_power*1.D9/(4.D0*PI)*pulse_duration*1.D-15*1.D6*energy(1:rhodist)
              CALL extend_1D_dset_unlimited(file_id, energy_dset_name, & 
                REAL(2.D0*PI*energy(1:count)*delta_t*delta_r**2,4), new_dims, &
                memspace_dims, offset, hyperslab_size)
              
              ! extend energy_fil dataset and the normalised energy_fil dataset
            !   CALL extend_1D_dset_unlimited(file_id, energy_fil_dset_name, &
            !     REAL(2.D0*PI*energy_fil(1:count)*delta_t*delta_r**2,4), new_dims, &
            !     memspace_dims, offset, hyperslab_size)

              energy_fil(1:rhodist)=1.D-6*critical_power*1.D9/(4.D0*PI)*pulse_duration*1.D-15*1.D6*energy_fil(1:rhodist)
              CALL extend_1D_dset_unlimited(file_id, energy_fil_dset_name, & 
                REAL(2.D0*PI*energy_fil(1:count)*delta_t*delta_r**2,4), new_dims, &
                memspace_dims, offset, hyperslab_size)

              
            ENDIF
            ! extend powmax and on axis data dataset
            !powmax_data(1,1) = REAL(z,4)
            !powmax_data(1,2) = REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)
            !z_data(1) = REAL(z,4)
            !onax_t_data(1,:) = REAL(ABS(e(1:dim_t,1)),4)
            DO k1=1, dim_t
              onax_t_data(1,k1) = REAL( efield_factor * REAL( efield_osc(k1)*e(k1,1) ) , 4 )
              !fields_array(1,k1,k2) = REAL( REAL( (efield_factor*efield_osc(k2)*e(k2,r_offset+k1)) ) , 4 ) ! SINGLE PRECISION, corresponding H5T_NATIVE_REAL (REAL(.,8) corresponds to H5T_NATIVE_DOUBLE)
              ! e(t,r)
            ENDDO 
            CALL extend_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)/), &
              new_dims=(/int(dset_write_count + 1,HSIZE_T)/), memspace_dims=(/int(1,HSIZE_T)/), &
              offset = (/int(dset_write_count, HSIZE_T)/), hyperslab_size = (/int(1,HSIZE_T)/))
            CALL extend_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, & 
              new_dims = (/int(dset_write_count + 1, HSIZE_T), int(dim_t, HSIZE_T)/), & 
              memspace_dims = (/int(1,HSIZE_T), int(dim_t, HSIZE_T)/), & 
              offset = (/int(dset_write_count,HSIZE_T),int(0,HSIZE_T)/), & 
              hyperslab_size = (/int(1,HSIZE_T), int(dim_t, HSIZE_T)/))
            CALL extend_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh*z,4)/), &
              new_dims=(/int(dset_write_count + 1,HSIZE_T)/), memspace_dims=(/int(1,HSIZE_T)/), &
              offset = (/int(dset_write_count, HSIZE_T)/), hyperslab_size = (/int(1,HSIZE_T)/))
            dset_write_count = dset_write_count + 1
            longstep_write_count = longstep_write_count + 1
          ENDIF
          
          ! Terminate HDF5 file access
          CALL h5fclose_f(file_id, error)
          CALL h5close_f(error)
       ENDIF

       ! Copy arrays to arrays with target flag
       fluence_data(:) = REAL(fluence*delta_t,4)
       plasma_channel_data(:) = REAL(rho,4)
       losses_plasma_data(:) = REAL(losses_plasma*delta_t,4)  
       losses_ionization_data(:) = REAL(losses_ionization*delta_t,4) 

       ! Create pointers to the arrays with target flag
       ptr_f => fluence_data
       ptr_p => plasma_channel_data
       ptr_lp => losses_plasma_data
       ptr_li => losses_ionization_data
       ! fluence_ll, plasma_channel_ll , etc store the pointers to the first link of the linked list of each variable
       IF (length_of_linked_list .EQ. 0) THEN
         ! initialize the linked lists
         call list_init(fluence_ll, DATA=transfer(ptr_f, list_data))
         call list_init(plasma_channel_ll, DATA=transfer(ptr_p, list_data))
         call list_init(losses_plasma_ll, DATA=transfer(ptr_lp, list_data))
         call list_init(losses_ionization_ll, DATA=transfer(ptr_li, list_data))
       ELSE
         ! append next link to the end of the linked lists
         call list_append(fluence_ll, DATA=transfer(ptr_f, list_data))
         call list_append(plasma_channel_ll, DATA=transfer(ptr_p, list_data))
         call list_append(losses_plasma_ll, DATA=transfer(ptr_lp, list_data))
         call list_append(losses_ionization_ll, DATA=transfer(ptr_li, list_data))
       ENDIF

       ! deallocate pointers to solve outlive warning
       !PRINT*, 'nullify'
       NULLIFY(ptr_f,ptr_p,ptr_lp,ptr_li)

       ! Increase the length counter
       length_of_linked_list = length_of_linked_list + 1
       ! Reset the counter
       count=0
    ELSE IF (z.LE.delta_z) THEN
       ! Copy arrays to arrays with target flag
       fluence_data(:) = REAL(fluence*delta_t,4)
       plasma_channel_data(:) = REAL(rho,4)
       losses_plasma_data(:) = REAL(losses_plasma*delta_t,4)  
       losses_ionization_data(:) = REAL(losses_ionization*delta_t,4) 

       ! Create pointers to the arrays with target flag
       ptr_f => fluence_data
       ptr_p => plasma_channel_data
       ptr_lp => losses_plasma_data
       ptr_li => losses_ionization_data
       ! fluence_ll, plasma_channel_ll , etc store the pointers to the first link of the linked list of each variable
       IF (length_of_linked_list .EQ. 0) THEN
         ! initialize the linked lists
         call list_init(fluence_ll, DATA=transfer(ptr_f, list_data))
         call list_init(plasma_channel_ll, DATA=transfer(ptr_p, list_data))
         call list_init(losses_plasma_ll, DATA=transfer(ptr_lp, list_data))
         call list_init(losses_ionization_ll, DATA=transfer(ptr_li, list_data))
       ELSE
         ! append next link to the end of the linked lists
         call list_append(fluence_ll, DATA=transfer(ptr_f, list_data))
         call list_append(plasma_channel_ll, DATA=transfer(ptr_p, list_data))
         call list_append(losses_plasma_ll, DATA=transfer(ptr_lp, list_data))
         call list_append(losses_ionization_ll, DATA=transfer(ptr_li, list_data))
       ENDIF

       ! deallocate pointers to solve outlive warning
       !PRINT*, 'nullify'
       NULLIFY(ptr_f,ptr_p,ptr_lp,ptr_li)

       ! Increase the length counter
       length_of_linked_list = length_of_linked_list + 1
       e_2=0.D0
       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
          e_2=e_2+ABS(e(1:dim_t,l))**2*REAL(l-1,8)
       ENDDO
       CALL MPI_REDUCE(e_2(1:dim_t),e_2KK(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF(my_rank.EQ.0) THEN
          CALL h5open_f(error) ! Prepare the HDF5 writing
          CALL h5fopen_f(main_h5_fname, H5F_ACC_RDWR_F, file_id, error ) ! Open the HDF5 file
          ! Create long_step group if it does not exist yet
          CALL h5lexists_f(file_id, longstep_grpname, group_status, error)
          IF ( group_status .EQV. .FALSE. ) THEN
            CALL h5gcreate_f(file_id, longstep_grpname, group_id, error) 
            CALL h5gclose_f(group_id, error)
          ENDIF
          ! store data of maximal power to a variable
          powmax_data(1,1) = REAL(z,4)
          powmax_data(1,2) = REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)
          !z_data(1) = REAL(four_z_Rayleigh*z,4)
          ! calculate on axis data
          !onax_t_data(1,:) = REAL(ABS(e(1:dim_t,1)),4)
          DO k1=1, dim_t
            onax_t_data(1,k1) = REAL( efield_factor * REAL( efield_osc(k1)*e(k1,1) ) , 4 )
            !fields_array(1,k1,k2) = REAL( REAL( (efield_factor*efield_osc(k2)*e(k2,r_offset+k1)) ) , 4 ) ! SINGLE PRECISION, corresponding H5T_NATIVE_REAL (REAL(.,8) corresponds to H5T_NATIVE_DOUBLE)
            ! e(t,r)
          ENDDO 
          ! write to datasets, either create
          IF ( dset_write_count .EQ. 0 ) THEN
            CALL create_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)/), 1)
            CALL create_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, dim_t)
            CALL create_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh*z,4)/), 1)
            CALL h5_add_units_1D(file_id, powmax_dset_name, '[C.U.]2')
            CALL h5_add_units_1D(file_id, onax_t_dset_name, '[V/m]?2')
            CALL h5_add_units_1D(file_id, every_rhodist_z_dset_name, '[m]2')
          ! or extend the existing ones
          ELSE
            CALL extend_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0*PI*MAXVAL(e_2KK)*delta_r**2,4)/), &
              new_dims=(/int(dset_write_count + 1,HSIZE_T)/), memspace_dims=(/int(1,HSIZE_T)/), &
              offset = (/int(dset_write_count, HSIZE_T)/), hyperslab_size = (/int(1,HSIZE_T)/))
            CALL extend_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, & 
              new_dims = (/int(dset_write_count + 1, HSIZE_T), int(dim_t, HSIZE_T)/), & 
              memspace_dims = (/int(1,HSIZE_T), int(dim_t, HSIZE_T)/), & 
              offset = (/int(dset_write_count,HSIZE_T),int(0,HSIZE_T)/), & 
              hyperslab_size = (/int(1,HSIZE_T), int(dim_t, HSIZE_T)/))
            CALL extend_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh*z,4)/), &
              new_dims=(/int(dset_write_count + 1,HSIZE_T)/), memspace_dims=(/int(1,HSIZE_T)/), &
              offset = (/int(dset_write_count, HSIZE_T)/), hyperslab_size = (/int(1,HSIZE_T)/))
          ENDIF
          ! increase the counter
          dset_write_count = dset_write_count + 1
          
          ! Terminate HDF5 file access
          CALL h5fclose_f(file_id, error)
          CALL h5close_f(error)
       ENDIF

       ! deallocate arrays used for writting to linked list
       DEALLOCATE(fluence_data, plasma_channel_data, losses_plasma_data, losses_ionization_data)
    ENDIF
    RETURN
  END SUBROUTINE mult_phase
  
  
!> @brief One propagation step using second order spectral operator splitting. In between two linear half-steps (\ref mult_propagator), the nonlinear step is executed (\ref mult_phase). The FFTs can be found in the module \ref fft. The subroutine \ref absorption implements absorbing boundary conditions in time.

  SUBROUTINE propagation
    USE fields
    USE fft
    USE mpi_stuff

    IMPLICIT NONE
    
    CALL fft_forward_inplace(.FALSE.)
    CALL mult_propagator
    CALL fft_backward2_inplace
    z=z+delta_zh
    CALL mult_phase
    CALL fft_forward_inplace(.TRUE.)
    CALL mult_propagator
    CALL fft_backward2_inplace
    CALL absorbation
    z=z+delta_zh

    RETURN
  END SUBROUTINE propagation

END MODULE long_step
