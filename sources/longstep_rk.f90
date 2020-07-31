MODULE long_step
  USE parameters
  REAL(8) rhompi,rho1,rho2,rhoth,rhotr,rhofh,rhoslg2,rhoav,rhoO2,rhoN2,Tev
CONTAINS

  SUBROUTINE index_interpolation(phase_index,r)
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
  
  SUBROUTINE absorbation
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
    B(1)=(CMPLX(1.D0,0.D0,8)+CMPLX(0.D0,-2.D0,8)*delta_rel(k))*B(1)+CMPLX(0.D0,2.D0,8)*delta_rel(k)*B(2)
    DO j=2,dim_r-1
       help_2=help_1
       help_1=B(j)
       B(j)=-DL(j-1,k)*help_2+(CMPLX(1.D0,0.D0,8)+CMPLX(0.D0,-1.D0,8)*delta_rel(k))*B(j)-DU(j,k)*B(j+1)
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

  SUBROUTINE mult_propagator
    USE fields
    USE mpi_stuff
    IMPLICIT NONE

    INTEGER(4) k

    DO k=dim_t_start(num_proc),dim_t_end(num_proc)
       CALL cn(efft(1:dim_r,k),k)
       efft(1:dim_r,k)=efft(1:dim_r,k)*p_t(k)
    ENDDO

    RETURN
  END SUBROUTINE mult_propagator

  SUBROUTINE calc_absorption(rhoabs, mediumabs, eti, etip1)
    REAL(8) :: rhoabs, mediumabs, eti, etip1
    REAL(8) :: intF

    IF ( eta1.NE.0.D0 ) THEN
       intF = -eta2 * 0.5d0 *( eti**Nn + etip1**NN) * delta_t  
       rhoabs = 1.d0 - (1.D0 - rhoabs) * exp(intF)
       IF (rhoabs.LT.0.D0) rhoabs = 0.D0
       mediumabs = eta1 * etip1**(NN-1)*(1.D0-rhoabs)
    ENDIF

    RETURN
  END SUBROUTINE calc_absorption

  SUBROUTINE calc_rho(rho,mpa,eti,etip1)
    USE PPT
    USE Complex_rotation
    IMPLICIT NONE

    REAL(8) :: rho,mpa,mpa_N2
    REAL(8) :: eti,etip1
    REAL(8) intF,var1,var1_N2,rhosave,trapped,colfreqO2p,colfreqN2p,colfreqO2,colfreqN2

    rhosave=rho
    SELECT CASE (switch_rho)
    CASE(8)
       CALL interpolate_cpr(var1,mpa,etip1)
       intF=(nu*0.5d0*(etip1+eti)-rhoat_inv*var1-alpha)*delta_t
       rho=rho*exp(intF)+var1*delta_t
       mpa=mpa*(1.D0-rhosave*rhoat_inv)
    CASE(1)
       rho=rho+(nu*rho*eti+beta_inv_2KK*eti**KK*(1.D0-rho*rhoat_inv)-alpha*rho)*delta_t
       mpa=muk*etip1**(KK-1)*(1.D0-rhosave*rhoat_inv)
    CASE(2)
       var1=0.5d0*beta_inv_2KK*(etip1**KK+eti**KK)
       intF=(nu*0.5d0*(etip1+eti)-rhoat_inv*var1-alpha)*delta_t
       rho=rho*exp(intF)+var1*delta_t
       mpa=muk*etip1**(KK-1)*(1.D0-rhosave*rhoat_inv)
    CASE(3,4)
       CALL interpolate_ppt(var1,mpa,etip1)
       intF=(nu*0.5d0*(etip1+eti)-rhoat_inv*var1-alpha)*delta_t
       rho=rho*exp(intF)+var1*delta_t
       mpa=mpa*(1.D0-rhosave*rhoat_inv)
    CASE(5)
       CALL interpolate_ppt(var1,mpa,etip1)
       intF=(nu*0.5d0*(etip1+eti)-rhoat_inv*var1-alpha)*delta_t
       rho=rho*exp(intF)+var1*delta_t
       mpa=mpa*(1.D0-rhosave*rhoat_inv)
    CASE(6)
       rhompi=rhompi+(beta_inv_2KK*eti**KK*(1.D0-rho*rhoat_inv)-alpha*rhompi)*delta_t
       rho1=rho1+(beta_inv_2KKp*eti**KKp-alpha1*rho1)*delta_t
       trapped=MIN(rhoth-rhotr,alpha2*(eti/eti_ref+1.D-6)**exp_ref*rho2*(rhoth-rhotr)*delta_t)
       rho2=rho2+beta_inv_2*eti*rhoslg2*delta_t-trapped
       rhotr=rhotr+trapped
       rhoth=rhoth+alphah*rhofh*delta_t
       rhofh=rhofh+(beta_inv_2*eti*rhoslg2-alphah*rhofh)*delta_t
       mpa=muk*etip1**(KK-1)*(1.D0-rho*rhoat_inv)+mukp*etip1**(KKp-1)+mu*rhoslg2+mukpp*etip1**(KKpp-1)
       rhoslg2=MIN(rhoslg2+(beta_inv_2KKpp*eti**KKpp-beta_inv_2*eti*rhoslg2)*delta_t,rhosat)
       rhoav=rhoav+(nu*rho*eti-alpha*rhoav)*delta_t
       rho=rhompi+rhoav+rho1+rho2
    CASE(7)
       CALL interpolate_ppt(var1,mpa,etip1)
       CALL interpolate_ppt_N2(var1_N2,mpa_N2,etip1)
       colfreqO2p=nucp*rhoO2*Tev**(-1.5D0)
       colfreqN2p=nucp*rhoN2*Tev**(-1.5D0)
       colfreqO2=nucO2*(1.D0-rhoO2*rhoat_inv)*SQRT(Tev)
       colfreqN2=nucN2*(1.D0-rhoN2*rhoat_N2_inv)*SQRT(Tev)
       mpa=mpa*(1.D0-rhoO2*rhoat_inv)+mpa_N2*(1.D0-rhoN2*rhoat_N2_inv)
       intF=(-rhoat_inv*var1-alpha)*delta_t
       rhoO2=rhoO2*exp(intF)+var1*delta_t+nuO2*colfreqO2*rho*0.5d0*(etip1+eti)*delta_t
       intF=(-rhoat_N2_inv*var1_N2-alpha)*delta_t
       rhoN2=rhoN2*exp(intF)+var1_N2*delta_t+nuN2*colfreqN2*rho*0.5d0*(etip1+eti)*delta_t
       gamma1=gamma1e*(colfreqO2p+colfreqN2p+colfreqO2+colfreqN2)
       intF=(-(nuO2*colfreqO2+nuN2*colfreqN2)*0.5d0*(etip1+eti))*delta_t
       Tev=Tev*exp(intF)+nukB*(colfreqO2p+colfreqN2p)*0.5d0*(etip1+eti)*delta_t
       IF (rhoO2*rhoat_inv.GT.1.D0) rhoO2=1.D0/rhoat_inv
       IF (rhoN2*rhoat_N2_inv.GT.1.D0) rhoN2=1.D0/rhoat_N2_inv
       rho=rhoO2+rhoN2
    END SELECT
    rho=rho-alphaquad*rhosave**2*delta_t
    IF (rho.LT.0.D0) rho=0.D0
    IF (rhoO2.LT.0.D0) rhoO2=0.D0
    IF (rhoN2.LT.0.D0) rhoN2=0.D0
    IF (Tev.LT.0.D0) Tev=0.D0
    IF (switch_rho.NE.7) THEN
       IF (rho*rhoat_inv.GT.1.D0) rho=1.D0/rhoat_inv
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

  SUBROUTINE mult_phase
    USE fields
    USE fft
    USE mpi_stuff
    USE HDF5
    USE HDF5_helper
    USE longstep_vars

    IMPLICIT NONE

    INTEGER(4) j,l,k
    REAL(8)  phase,maxphase_part,peakmax_part,energy_part,energy_fil_part,rhomax_part,delkerr,delkerrp,rhotemp,mpa,r,phase_p,phase_j,losses_j,phase_index
    REAL(8) mediumabs , rhoabs_max_part, rhoabstemp, rhoO2max_part, rhoN2max_part, Tevmax_part

    ! For storing in HDF5
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: group_id      ! Group identifier 
    INTEGER(HID_T) :: h5parameters  ! Property list identifier 
    INTEGER        :: error
    LOGICAL        :: group_status
    CHARACTER(LEN=15) :: h5_filename="results.h5"
    CHARACTER(LEN=15) :: groupname="longstep"
    CHARACTER(LEN=25) :: rhoabs_max_dset_name="longstep/rhoexcmax_max"
    CHARACTER(LEN=25) :: peakmax_dset_name="longstep/peakmax"
    CHARACTER(LEN=25) :: energy_dset_name="longstep/energy"
    CHARACTER(LEN=25) :: energy_fil_dset_name="longstep/energy_fil"
    CHARACTER(LEN=25) :: rhomax_dset_name="longstep/rhomax"
    CHARACTER(LEN=25) :: powmax_dset_name="longstep/powmax"
    INTEGER(HSIZE_T), DIMENSION(2) :: new_dims, memspace_dims, offset, hyperslab_size

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
    phaserelation(1)=exp(CMPLX(0.D0,2.D0*(rek0-rekp*omega)*z,8))
    phaserelation(2)=exp(CMPLX(0.D0,4.D0*(rek0-rekp*omega)*z,8))
    rhoO2max_part=0.D0
    rhoN2max_part=0.D0
    Tevmax_part=0.D0

    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
       r=REAL(l-1)*delta_r
       e_2=ABS(e(1:dim_t,l))
       peakmax_part=MAX(peakmax_part,MAXVAL(e_2))
       e_2=e_2**2
       fluence(l)=SUM(e_2)
       energy_part=energy_part+fluence(l)*REAL(l-1,8)
       IF (rfil.GT.r) energy_fil_part=energy_fil_part+fluence(l)*REAL(l-1,8)
       delkerr=0.D0
       delkerrp=0.d0
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
       rhoabstemp=0.D0
       mediumabs=0.D0
       CALL index_interpolation(phase_index,r)
       phase_index=phase_index*delta_zh
       CALL calc_rho(rhotemp,mpa,0.D0,e_2(1))
       CALL calc_absorption(rhoabstemp, mediumabs, 0.D0,e_2(1))
       DO j=1,dim_t
          IF (switch_rho.EQ.7) THEN
             phase_p=(c3i*e_2(j)+c3d*delkerr)*((1.D0-rhoO2*rhoat_inv)/3.D0+2.D0*(1.D0-rhoN2*rhoat_N2_inv)/3.D0)*delta_zh
          ELSE
             phase_p=(c3i*e_2(j)+c3d*delkerr-c5*e_2(j)**2)*((1.D0-rhotemp*rhoat_inv)+rhotemp*rhoat_inv/3.0d0)*delta_zh
          ENDIF
          phase_j=-gamma2*rhotemp*delta_zh
          losses_j=-gamma1*rhotemp*delta_zh
          rho(l)=MAX(rho(l),rhotemp)
          rhoabs(l) = MAX(rhoabs(l),rhoabstemp)
          losses_plasma(l)=losses_plasma(l)+2.D0*e_2(j)*gamma1*rhotemp
          losses_ionization(l)=losses_ionization(l)+2.D0*e_2(j)*mpa
          phase=phase_p+phase_j+phase_index
          maxphase_part=MAX(maxphase_part,ABS(phase))
          rhoO2max_part=MAX(rhoO2max_part,rhoO2)
          rhoN2max_part=MAX(rhoN2max_part,rhoN2)
          Tevmax_part=MAX(Tevmax_part,Tev)
          SELECT CASE (switch_T)
          CASE(1)
             etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase,8))
          CASE(2)
             ptemp(j,l)=e(j,l)*CMPLX(0.D0,phase_p)
             jtemp(j,l)=e(j,l)*CMPLX(0.D0,phase_j)
             etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase_index,8))
          CASE(3)
             ptemp(j,l)=e(j,l)*CMPLX(0.D0,phase_p) + hfac(j,4)*CONJG(hfac(j,0)*e(j,l)*phaserelation(1))*CMPLX(0.D0,phase_p) &
                  + (hfac(j,1)*(hfac(j,0)*e(j,l))**3 + hfac(j,2)*(hfac(j,0)*e(j,l))**3*e_2(j))*phaserelation(1) &
                  + (hfac(j,3)*(hfac(j,0)*e(j,l))**5)*phaserelation(2)
             jtemp(j,l)=e(j,l)*CMPLX(0.D0,phase_j)
             etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase_index,8))
          CASE(4)
             ptemp(j,l)=e(j,l)*CMPLX(0.D0,phase_p)
             jtemp(j,l)=e(j,l)*CMPLX(0.D0,phase_j)
             etemp(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_zh*(mpa+mediumabs),phase_index,8))
          END SELECT
          IF (j.NE.dim_t) THEN
             CALL calc_rho(rhotemp,mpa,e_2(j),e_2(j+1))
             CALL calc_absorption(rhoabstemp, mediumabs, e_2(j), e_2(j+1))
             CALL calc_delkerr(delkerr,delkerrp,e_2(j),e_2(j+1))
          ENDIF
       ENDDO
    ENDDO

    SELECT CASE (switch_T)
    CASE(1)
       continue
    CASE(2)
       CALL dfftw_execute(plan_forward_erk)
       CALL dfftw_execute(plan_p)
       CALL dfftw_execute(plan_j)
       DO k=dim_r_start(num_proc),dim_r_end(num_proc)
          DO l=1,dim_t
             etemp(l,k)=etemp(l,k)+op_t(l)*ptemp(l,k)+op_t_inv(l)*jtemp(l,k)
          ENDDO
       ENDDO
       CALL dfftw_execute(plan_backward_erk)
       etemp=diminv*etemp
    CASE(3)
       CALL dfftw_execute(plan_forward_erk) 
       CALL dfftw_execute(plan_p)
       CALL dfftw_execute(plan_j)
       DO k=dim_r_start(num_proc),dim_r_end(num_proc)
          DO l=dim_th+1,dim_t
             etemp(l,k)=etemp(l,k)+op_t(l)*ptemp(l,k)+op_t_inv(l)*jtemp(l,k)
          ENDDO
       ENDDO
       CALL dfftw_execute(plan_backward_erk)
       etemp=diminv*etemp
    CASE(4)
       CALL dfftw_execute(plan_forward_erk)
       CALL dfftw_execute(plan_p)
       CALL dfftw_execute(plan_j)
       DO k=dim_r_start(num_proc),dim_r_end(num_proc)
          DO l=1,dim_t
             etemp(l,k)=etemp(l,k)+op_t(l)*ptemp(l,k)+op_t_inv(l)*jtemp(l,k)
          ENDDO
       ENDDO
       CALL dfftw_execute(plan_backward_erk)
       etemp=diminv*etemp
    END SELECT

    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
       r=REAL(l-1)*delta_r
       e_2=ABS(etemp(1:dim_t,l))
       e_2=e_2**2
       delkerr=0.D0
       delkerrp=0.d0
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
       rhoabstemp=0.D0
       mediumabs=0.D0
       CALL index_interpolation(phase_index,r)
       phase_index=phase_index*delta_z
       CALL calc_rho(rhotemp,mpa,0.D0,e_2(1))
       CALL calc_absorption(rhoabstemp, mediumabs, 0.D0,e_2(1))
       DO j=1,dim_t
          IF (switch_rho.EQ.7) THEN
             phase_p=(c3i*e_2(j)+c3d*delkerr)*((1.D0-rhoO2*rhoat_inv)/3.D0+2.D0*(1.D0-rhoN2*rhoat_N2_inv)/3.D0)*delta_z
          ELSE
             phase_p=(c3i*e_2(j)+c3d*delkerr-c5*e_2(j)**2)*((1.D0-rhotemp*rhoat_inv)+rhotemp*rhoat_inv/3.0d0)*delta_z
          ENDIF
!          phase_p=(c3i*e_2(j)+c3d*delkerr-c5*e_2(j)**2)*delta_z
          phase_j=-gamma2*rhotemp*delta_z
          losses_j=-gamma1*rhotemp*delta_z
          phase=phase_p+phase_j+phase_index
          SELECT CASE (switch_T)
          CASE(1)
             e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase,8))
          CASE(2)
             ptemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_p)
             jtemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_j)
             e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase_index,8))
          CASE(3)
             ptemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_p) + hfac(j,4)*CONJG(hfac(j,0)*etemp(j,l)*phaserelation(1))*CMPLX(0.D0,phase_p) &
                  + 2.D0*(hfac(j,1)*(hfac(j,0)*etemp(j,l))**3 + hfac(j,2)*(hfac(j,0)*etemp(j,l))**3*e_2(j))*phaserelation(1) &
                  + 2.D0*(hfac(j,3)*(hfac(j,0)*etemp(j,l))**5)*phaserelation(2)
             jtemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_j)
             e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase_index,8))
          CASE(4)
             ptemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_p)
             jtemp(j,l)=etemp(j,l)*CMPLX(0.D0,phase_j)
             e(j,l)=e(j,l)*exp(CMPLX(losses_j-delta_z*(mpa+mediumabs),phase_index,8))
          END SELECT
          IF (j.NE.dim_t) THEN
             CALL calc_rho(rhotemp,mpa,e_2(j),e_2(j+1))
             CALL calc_absorption(rhoabstemp, mediumabs, e_2(j), e_2(j+1))
             CALL calc_delkerr(delkerr,delkerrp,e_2(j),e_2(j+1))
          ENDIF
       ENDDO
    ENDDO

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
    CALL MPI_REDUCE(rhoO2max_part,rhoO2max(count),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(rhoN2max_part,rhoN2max(count),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(Tevmax_part,Tevmax(count),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    IF(my_rank.EQ.0) z_buff(count)=z
    IF(count.GE.rhodist) THEN
       e_2=0.D0
       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
          e_2=e_2+ABS(e(1:dim_t,l))**2*REAL(l-1,8)
       ENDDO
       CALL MPI_REDUCE(e_2(1:dim_t),e_2KK(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF(my_rank.EQ.0) THEN
          CALL h5open_f(error) ! Prepare the HDF5 writing
          CALL h5fopen_f(h5_filename, H5F_ACC_RDWR_F, file_id, error ) ! Open the HDF5 file
          ! Create long_step group if it does not exist yet
          CALL h5lexists_f(file_id, groupname, group_status, error)
          IF ( group_status .EQV. .FALSE. ) THEN
            CALL h5gcreate_f(file_id, groupname, group_id, error) 
            CALL h5gclose_f(group_id, error)
          ENDIF
          OPEN(unit_rhoabs_max, FILE='rhoexcmax.dat',STATUS='UNKNOWN',POSITION='APPEND') 
          OPEN(unit_peakmax,FILE='peakmax.dat',STATUS='UNKNOWN',POSITION='APPEND')
          OPEN(unit_energy,FILE='energy.dat',STATUS='UNKNOWN',POSITION='APPEND')
          OPEN(unit_rhomax,FILE='rhomax.dat',STATUS='UNKNOWN',POSITION='APPEND')
          OPEN(unit_rho,FILE='powmax.dat',STATUS='UNKNOWN',POSITION='APPEND')
          DO j=1,count
             WRITE(unit_rhoabs_max,*) REAL(z_buff(j),4) ,REAL(rhoabs_max(j),4)
             WRITE(unit_peakmax,*) REAL(z_buff(j),4),REAL(peakmax(j),4)
             WRITE(unit_rhomax,*) REAL(z_buff(j),4) ,REAL(rhomax(j),4)
             WRITE(unit_energy,*) REAL(z_buff(j),4),REAL(6.2831853D0*energy(j)*delta_t*delta_r**2,4), &
                  REAL(6.2831853D0*energy_fil(j)*delta_t*delta_r**2,4)
             print *,z_buff(j)
          ENDDO
          IF (number_of_steps.EQ.0) THEN
            ALLOCATE(data_to_write(1,1:rhodist))
            data_to_write(1,1:rhodist) = REAL(rhoabs_max(1:rhodist),4)
            CALL create_2D_dset_unlimited(file_id, rhoabs_max_dset_name, data_to_write, (count))
            data_to_write(1,1:rhodist) = REAL(peakmax(1:rhodist),4)
            CALL create_2D_dset_unlimited(file_id, peakmax_dset_name, data_to_write, (count))
            data_to_write(1,1:rhodist) = REAL(rhomax(1:rhodist),4)
            CALL create_2D_dset_unlimited(file_id, rhomax_dset_name, data_to_write, (count))
            data_to_write(1,1:rhodist) = REAL(6.2831853D0*energy(1:rhodist)*delta_t*delta_r**2,4)
            CALL create_2D_dset_unlimited(file_id, energy_dset_name, data_to_write, (count))
            data_to_write(1,1:rhodist) = REAL(6.2831853D0*energy_fil(1:rhodist)*delta_t*delta_r**2,4)
            CALL create_2D_dset_unlimited(file_id, energy_fil_dset_name, data_to_write, (count))
            number_of_steps = number_of_steps + 1
            original_rhodist = rhodist
          ELSE
            number_of_steps = number_of_steps + 1
            IF (rhodist.NE.original_rhodist) THEN
              new_dims = (/int(number_of_steps, 8), int(original_rhodist, 8)/)
            ELSE
              new_dims = (/int(number_of_steps, 8), int(rhodist, 8)/)
            ENDIF
            memspace_dims = (/int(1,HSIZE_T), int(original_rhodist,HSIZE_T)/)
            offset = (/int(number_of_steps-1,HSIZE_T),int(0,HSIZE_T)/)
            hyperslab_size = (/int(1,HSIZE_T),int(original_rhodist,HSIZE_T)/)
            IF (count .EQ. 0) THEN
              print *,"last write skipped"
            ELSE
              data_to_write(1,1:count) = REAL(rhoabs_max(1:count),4)
              IF (count .NE. original_rhodist) THEN
                DO j=count+1,original_rhodist
                  data_to_write(1,j) = 0d0
                ENDDO
              ENDIF
              CALL extend_2D_dset_unlimited(file_id, rhoabs_max_dset_name, data_to_write, new_dims, &
                memspace_dims, offset, hyperslab_size)
              
              data_to_write(1,1:count) = REAL(peakmax(1:count),4)
              CALL extend_2D_dset_unlimited(file_id, peakmax_dset_name, data_to_write, new_dims, &
                memspace_dims, offset, hyperslab_size)
              
              data_to_write(1,1:count) = REAL(rhomax(1:count),4)
              CALL extend_2D_dset_unlimited(file_id, rhomax_dset_name, data_to_write, new_dims, &
                memspace_dims, offset, hyperslab_size)

              data_to_write(1,1:count) = REAL(6.2831853D0*energy(1:count)*delta_t*delta_r**2,4)
              CALL extend_2D_dset_unlimited(file_id, energy_dset_name, data_to_write, new_dims, &
                memspace_dims, offset, hyperslab_size)
              
              data_to_write(1,1:count) = REAL(6.2831853D0*energy_fil(1:count)*delta_t*delta_r**2,4)
              CALL extend_2D_dset_unlimited(file_id, energy_fil_dset_name, data_to_write, new_dims, &
                memspace_dims, offset, hyperslab_size)
            ENDIF
          ENDIF
          WRITE(unit_rho,*) REAL(z,4),REAL(6.2831853D0*MAXVAL(e_2KK)*delta_r**2,4)
          CLOSE(unit_rhoabs_max)
          CLOSE(unit_peakmax)
          CLOSE(unit_rhomax)
          CLOSE(unit_energy)
          CLOSE(unit_rho)
          IF (switch_rho.EQ.7) THEN
             OPEN(unit_rho,FILE='rho_pavel.dat',STATUS='UNKNOWN',POSITION='APPEND')
             DO j=1,count
                WRITE(unit_rho,*) REAL(z_buff(j),4) ,REAL(rhoO2max(j),4) ,REAL(rhoN2max(j),4) ,REAL(Tevmax(j),4)
             ENDDO
             CLOSE(unit_rho)
          ENDIF
          OPEN(unit_rho,FILE='ONAX_T.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
          WRITE(unit_rho) REAL(z,4),REAL(ABS(e(1:dim_t,1)),4)
          CLOSE(unit_rho)
          ! Terminate HDF5 file access
          CALL h5fclose_f(file_id, error)
          CALL h5close_f(error)
       ENDIF
       OPEN(unit_rho,FILE='FLUENCE_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(fluence*delta_t,4)
       CLOSE(unit_rho)
       OPEN(unit_rho,FILE='PLASMACHANNEL_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(rho,4)
       CLOSE(unit_rho)
       OPEN(unit_rho,FILE='LOSSES_PLASMA_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(losses_plasma*delta_t,4)
       CLOSE(unit_rho)
       OPEN(unit_rho,FILE='LOSSES_IONIZATION_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(losses_ionization*delta_t,4)
       CLOSE(unit_rho)
       count=0
    ENDIF
    IF (z.LE.delta_z) THEN
       OPEN(unit_rho,FILE='FLUENCE_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(fluence*delta_t,4)
       CLOSE(unit_rho)
       OPEN(unit_rho,FILE='PLASMACHANNEL_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(rho,4)
       CLOSE(unit_rho)
       OPEN(unit_rho,FILE='LOSSES_PLASMA_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(losses_plasma*delta_t,4)
       CLOSE(unit_rho)
       OPEN(unit_rho,FILE='LOSSES_IONIZATION_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
       WRITE(unit_rho) REAL(z,4),REAL(losses_ionization*delta_t,4)
       CLOSE(unit_rho)
       e_2=0.D0
       DO l=dim_r_start(num_proc),dim_r_end(num_proc)
          e_2=e_2+ABS(e(1:dim_t,l))**2*REAL(l-1,8)
       ENDDO
       CALL MPI_REDUCE(e_2(1:dim_t),e_2KK(1:dim_t),dim_t,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       IF(my_rank.EQ.0) THEN
          OPEN(unit_rho,FILE='powmax.dat',STATUS='UNKNOWN',POSITION='APPEND')
          WRITE(unit_rho,*) REAL(z,4),REAL(6.2831853D0*MAXVAL(e_2KK)*delta_r**2,4)
          CLOSE(unit_rho)
          OPEN(unit_rho,FILE='ONAX_T.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
          WRITE(unit_rho) REAL(z,4),REAL(ABS(e(1:dim_t,1)),4)
          CLOSE(unit_rho)
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE mult_phase

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
