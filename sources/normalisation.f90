MODULE normalisation
  USE calc_start

  REAL(8) n2_phys,n4_phys,lambda0_cm_phys,Ui_eV_phys,proplength_m_phys,outlength_m_phys,delta_z_mm_phys,sigmak_phys,sigman_phys
  REAL(8) sigma_cm2_phys,alpha_fs_phys,alphaquad_fscm3_phys,tdk_fs_phys,raman_phys,rfil_mm_phys,pressure,tauc_fs_phys
  REAL(8) alpha1_fs_phys,alphah_fs_phys,rho0_phys
  REAL(8) delta_k_p_fs_per_cm_phys,k_p_fs_per_cm_phys,k_pp_fs2_per_cm_phys,k_ppp_fs3_per_cm_phys,k_pppp_fs4_per_cm_phys
  REAL(8) k_ppppp_fs5_per_cm_phys,rhoabs_cm3_phys,f_cm_phys,betak_phys
  REAL(8) z_rayleigh_cm_phys,k0_phys
  REAL(8) deltak3omega(3),deltak5omega(3)
  REAL(8), ALLOCATABLE :: xx_mum(:),zz_mum(:),Indice(:,:)
  CHARACTER(15) dispfilename

  INTEGER(4) dim_chi
  REAL(8), ALLOCATABLE :: omegachi(:)
  COMPLEX(8), ALLOCATABLE :: chi(:)

CONTAINS
  
  SUBROUTINE compute_dispersion(switch_dispersion)
    IMPLICIT NONE

    INTEGER(4) switch_dispersion,j,i,thg,fhg,startcut,endcut
    REAL(8) c,delta,tod,fod,vod,k_t,help1,help2
    LOGICAL(4) mess1,mess2,cut

    c = 3.d10 !speed of light in the vacuum      cm/s
    k0_phys = 2.D0*3.1415D0/lambda0_cm_phys  !central wave number in vacuum     cm-1
    omega=c*k0_phys*tp_fs_phys*1.d-15 ! adimensioned frequency
    i=CEILING(omega*lt/(8.D0*DATAN(1.D0)))
    lt=8.D0*DATAN(1.D0)*REAL(i,8)/omega
    i=2*INT(1.D-1*omega*lt/(8.D0*DATAN(1.D0)))
    omega_uppe=MAX(omega,4.D0*DATAN(1.D0)*REAL(dim_t+i)/lt)
    ALLOCATE(komega(dim_t)) ! omega0 at dim_t/2+1

    SELECT CASE(switch_dispersion)

    CASE(1)
       IF (pressure.NE.1.D0) THEN
          PRINT*, 'WARNING: pressure dependency of taylor expanded',' dispersion law valid for n~1 only'
       ENDIF
       help1=(n0**2-1.D0)*pressure
       n0=sqrt(1.D0+help1)
       k_p_fs_per_cm_phys=(delta_k_p_fs_per_cm_phys)*pressure+1.D0/(c*1.D-15)
       k_pp_fs2_per_cm_phys=k_pp_fs2_per_cm_phys*pressure
       k_ppp_fs3_per_cm_phys=k_ppp_fs3_per_cm_phys*pressure
       k_pppp_fs4_per_cm_phys=k_pppp_fs4_per_cm_phys*pressure
       k_ppppp_fs5_per_cm_phys=k_ppppp_fs5_per_cm_phys*pressure
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       delta = 2.D0*z_rayleigh_cm_phys*(k_pp_fs2_per_cm_phys)/((tp_fs_phys)**2)  !adimensionned coefficient for GVD
       tod = (2.D0/3.D0)*z_rayleigh_cm_phys*(k_ppp_fs3_per_cm_phys)/((tp_fs_phys)**3)
       fod = (2.D0/12.D0)*z_rayleigh_cm_phys*(k_pppp_fs4_per_cm_phys)/((tp_fs_phys)**4)
       vod = (2.D0/60.D0)*z_rayleigh_cm_phys*(k_ppppp_fs5_per_cm_phys)/((tp_fs_phys)**5)
       rek0 = 4.D0*z_rayleigh_cm_phys*k0_phys*n0 !adimensioned central wavenumber in the medium
       rekp = 4.D0*z_rayleigh_cm_phys*k_p_fs_per_cm_phys/tp_fs_phys !adimensioned inverse group velocity
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = CMPLX(rek0 + rekp * k_t  + delta * k_t**2 + tod * k_t**3 + fod * k_t**4 + vod * k_t**5, 0.D0, 8)
          IF ((k_t.GT.5.D0*omega).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       k_t=2.D0*omega
       deltak3omega(1)=3.D0*rek0 - (rek0 + rekp * k_t  + delta * k_t**2 + tod * k_t**3 + fod * k_t**4 + vod * k_t**5) 
       k_t=4.D0*omega
       deltak5omega(1)=5.D0*rek0 - (rek0 + rekp * k_t  + delta * k_t**2 + tod * k_t**3 + fod * k_t**4 + vod * k_t**5)

    CASE(2)
       IF (pressure.NE.1.D0) THEN
          PRINT*, 'WARNING: pressure dependency for silica does not really make sense. Sure what you are doing?'
       ENDIF
       dim_chi=dim_t
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       DO j=dim_chi/2,dim_chi/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = 0.69661663D0 * (1.d0/0.0684043D0**2) / ( (1.d0/0.0684043D0**2) - omegachi(j)**2) + &
               0.4079426D0 * (1.D0/0.1162414D0**2) / ( (1.D0/0.1162414D0**2) - omegachi(j)**2) + &
               0.8974794D0 * (1.d0/9.896161D0**2) / ( (1.d0/9.896161D0**2) -  omegachi(j)**2 )      ! Sellmeier formula silica, Agrawal
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
       ENDDO
       n0 = REAL(chi(dim_t/2+1)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0))!adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = 0.69661663D0 * (1.d0/0.0684043D0**2) / ( (1.d0/0.0684043D0**2) - omegachi(j)**2) + &
               0.4079426D0 * (1.D0/0.1162414D0**2) / ( (1.D0/0.1162414D0**2) - omegachi(j)**2) + &
               0.8974794D0 * (1.d0/9.896161D0**2) / ( (1.d0/9.896161D0**2) -  omegachi(j)**2 )      ! Sellmeier formula silica, Agrawal
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_chi
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = 0.69661663D0 * (1.d0/0.0684043D0**2) / ( (1.d0/0.0684043D0**2) - omegachi(j)**2) + &
               0.4079426D0 * (1.D0/0.1162414D0**2) / ( (1.D0/0.1162414D0**2) - omegachi(j)**2) + &
               0.8974794D0 * (1.d0/9.896161D0**2) / ( (1.d0/9.896161D0**2) -  omegachi(j)**2 )      ! Sellmeier formula silica, Agrawal
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          IF ((omegachi(j).GT.1.D0/0.185D0).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       DEALLOCATE(chi,omegachi)

    CASE(3)
       dim_chi=dim_t
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       DO j=dim_chi/2,dim_chi/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 5.547d-4*(1.D0+5.15d5*omegachi(j)**2+4.19d11*omegachi(j)**4+4.09d17*omegachi(j)**6+4.32d23*omegachi(j)**8) !argon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       n0 = REAL(chi(dim_t/2+1)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0))!adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 5.547d-4*(1.D0+5.15d5*omegachi(j)**2+4.19d11*omegachi(j)**4+4.09d17*omegachi(j)**6+4.32d23*omegachi(j)**8) !argon, dalgarno and kingston 
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_chi
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 5.547d-4*(1.D0+5.15d5*omegachi(j)**2+4.19d11*omegachi(j)**4+4.09d17*omegachi(j)**6+4.32d23*omegachi(j)**8) !argon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          IF ((omegachi(j).GT.7.D-4).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       DEALLOCATE(chi,omegachi)

    CASE(4)
       dim_chi=2048
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       OPEN(10,FILE=dispfilename,STATUS='OLD')
       DO j=1,dim_chi
          READ(10,*) omegachi(j),help1,help2
          chi(j)=CMPLX(help1,help2,8)
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       CLOSE(10)
       omegachi=omegachi*tp_fs_phys*1.d-15-c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15
       n0 = REAL(compute_n(0.D0)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))* &
            compute_n(k_t)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0)) !adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))* &
            compute_n(k_t)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))* &
            compute_n(k_t)
       ENDDO
       mess1=.TRUE.
       mess2=.TRUE.
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          IF (k_t.GE. omegachi(dim_chi)) THEN
             komega(j)=CMPLX(REAL(komega(j)),AIMAG(komega(j-1)),8)
             IF (mess1) THEN
                PRINT*, 'Warning: frequency-range of ',dispfilename,' exceeded (to the right), values extrapolated'
                mess1=.FALSE.
             ENDIF
          ENDIF
       ENDDO
       DO j=dim_t,1,-1
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          IF (k_t.LE. omegachi(1)) THEN
             komega(j)=CMPLX(REAL(komega(j)),AIMAG(komega(j+1)),8)
             IF (mess2) THEN
                PRINT*, 'Warning: frequency-range of ',dispfilename,' exceeded (to the left), values extrapolated'
                mess2=.FALSE.
             ENDIF
          ENDIF
       ENDDO
       DEALLOCATE(chi,omegachi)

    CASE(5)
       dim_chi=dim_t
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       DO j=dim_chi/2,dim_chi/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = (1.D0+1.D-8*(8060.51+2480990/(132.274-omegachi(j)**2)+17455.7/(39.32957-omegachi(j)**2)))**2-1.D0 !air, peck and reeder
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       n0 = REAL(chi(dim_t/2+1)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0))!adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = (1.D0+1.D-8*(8060.51+2480990/(132.274-omegachi(j)**2)+17455.7/(39.32957-omegachi(j)**2)))**2-1.D0 !air, peck and reeder
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_chi
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = (1.D0+1.D-8*(8060.51+2480990/(132.274-omegachi(j)**2)+17455.7/(39.32957-omegachi(j)**2)))**2-1.D0 !air, peck and reeder
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          IF ((omegachi(j).GT.1.D0/0.185D0).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       DEALLOCATE(chi,omegachi)

    CASE(6)
       dim_chi=dim_t
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       DO j=dim_chi/2,dim_chi/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 1.366d-3*(1.D0+9.02d5*omegachi(j)**2+1.81d12*omegachi(j)**4+4.89d18*omegachi(j)**6+1.45d25*omegachi(j)**8) !xenon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       n0 = REAL(chi(dim_t/2+1)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0))!adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 1.366d-3*(1.D0+9.02d5*omegachi(j)**2+1.81d12*omegachi(j)**4+4.89d18*omegachi(j)**6+1.45d25*omegachi(j)**8) !xenon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_chi
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 1.366d-3*(1.D0+9.02d5*omegachi(j)**2+1.81d12*omegachi(j)**4+4.89d18*omegachi(j)**6+1.45d25*omegachi(j)**8) !xenon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          IF ((omegachi(j).GT.7.D-4).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       DEALLOCATE(chi,omegachi)

    CASE(7)
       dim_chi=dim_t
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       DO j=dim_chi/2,dim_chi/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 1.335d-4*(1.D0+2.24d5*omegachi(j)**2+8.09d10*omegachi(j)**4+3.56d16*omegachi(j)**6) !neon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       n0 = REAL(chi(dim_t/2+1)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0))!adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 1.335d-4*(1.D0+2.24d5*omegachi(j)**2+8.09d10*omegachi(j)**4+3.56d16*omegachi(j)**6) !neon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_chi
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 1.335d-4*(1.D0+2.24d5*omegachi(j)**2+8.09d10*omegachi(j)**4+3.56d16*omegachi(j)**6) !neon, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          IF ((omegachi(j).GT.7.D-4).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       DEALLOCATE(chi,omegachi)

    CASE(8)
       IF (pressure.NE.1.D0) THEN
          PRINT*, 'WARNING: pressure dependency for KDP does not really make sense. Sure what you are doing?'
       ENDIF
       dim_chi=dim_t
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       DO j=dim_chi/2,dim_chi/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = 2.259276D0 + 0.01008956D0/(omegachi(j)**(-2)-0.0129426D0) + &
            13.00522D0*omegachi(j)**(-2)/(omegachi(j)**(-2)-400.D0) - 1.D0 ! KDPo, Handbook of Optics
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
       ENDDO
       n0 = REAL(chi(dim_t/2+1)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0))!adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = 2.259276D0 + 0.01008956D0/(omegachi(j)**(-2)-0.0129426D0) + & 
            13.00522D0*omegachi(j)**(-2)/(omegachi(j)**(-2)-400.D0) - 1.D0 ! KDPo, Handbook of Optics
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_chi
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          chi(j) = 2.259276D0 + 0.01008956D0/(omegachi(j)**(-2)-0.0129426D0) + &
            13.00522D0*omegachi(j)**(-2)/(omegachi(j)**(-2)-400.D0) - 1.D0 ! KDPo, Handbook of Optics
          chi(j) = 1.D0+chi(j)*pressure
          IF ((k_t.LT.0.D0).AND.(REAL(chi(j)).LT.0.1D0)) THEN
             PRINT*, 'remove singularity below omega at',k_t
             chi(j)=CMPLX(0.1D0,0.D0,8)
          ENDIF
          chi(j) = sqrt(chi(j))
       ENDDO
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)   ! omegachi = w
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d4)
          IF ((omegachi(j).GT.1.D0/0.185D0).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       DEALLOCATE(chi,omegachi)
    CASE(9) ! Krypton
       dim_chi=dim_t
       ALLOCATE(chi(dim_chi),omegachi(dim_chi))
       DO j=dim_chi/2,dim_chi/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 8.377d-4*(1.D0+6.70d5*omegachi(j)**2+8.84d11*omegachi(j)**4+1.49d18*omegachi(j)**6+2.74d24*omegachi(j)**8+5.10d30*omegachi(j)**10) !krypton, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       n0 = REAL(chi(dim_t/2+1)) !refractive index at lambda0
       z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
       DO j=dim_t/2,dim_t/2+2
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       rek0 = REAL(komega(dim_t/2+1)) !adimensioned central wavenumber in the medium
       rekp = REAL(komega(dim_t/2+2)-komega(dim_t/2))*lt/(16.D0*DATAN(1.D0))!adimensioned group velocity

       DO j=1,2
          k_t=REAL(j,8)*2.D0*omega   ! k_t = w_adim - w0_adim
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 8.377d-4*(1.D0+6.70d5*omegachi(j)**2+8.84d11*omegachi(j)**4+1.49d18*omegachi(j)**6+2.74d24*omegachi(j)**8+5.10d30*omegachi(j)**10) !krypton, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))                ! chi = E(w) =  sqrt(1+X(w))
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       deltak3omega(1)=3.D0*rek0 - REAL(komega(1)) 
       deltak5omega(1)=5.D0*rek0 - REAL(komega(2))

       DO j=1,dim_chi
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          chi(j) = 8.377d-4*(1.D0+6.70d5*omegachi(j)**2+8.84d11*omegachi(j)**4+1.49d18*omegachi(j)**6+2.74d24*omegachi(j)**8+5.10d30*omegachi(j)**10) !krypton, dalgarno and kingston
          chi(j) = 1.D0+chi(j)*pressure
          chi(j) = sqrt(chi(j))
       ENDDO
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
       ENDDO
       startcut=dim_t
       endcut=dim_t
       DO j=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt+omega_uppe-omega
          komega(j) = (c*2.D0*3.1415D0/lambda0_cm_phys*tp_fs_phys*1.d-15+k_t)/(c*tp_fs_phys*1.d-15/(4*z_rayleigh_cm_phys))*chi(j)
          omegachi(j) = c*2.D0*3.1415D0/lambda0_cm_phys+k_t/(tp_fs_phys*1.d-15)
          omegachi(j) = omegachi(j)/(c*2.D0*3.1415D0*1.d8)
          IF ((omegachi(j).GT.7.D-4).AND.(startcut.EQ.dim_t)) startcut=j  
       ENDDO
       IF (startcut.LT.endcut) THEN
          CALL artifdisp(startcut,endcut,REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)* &
            REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega)))
          CALL artifabs(startcut,endcut,10.D0*REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(dim_t/2)/lt+omega_uppe-omega)))
       ENDIF
       DEALLOCATE(chi,omegachi)

    END SELECT

    mess1=.TRUE.
    mess2=.TRUE.
    DO j=1,dim_t
       if (AIMAG(komega(j)).LT.0.D0) then
          IF (mess1) THEN
             PRINT*, 'Warning: negative imaginary part of k,',' changed to 0, check your dispersion law'
             PRINT*, 'Check komega_original.dat and komega_original_reduced.dat'
             OPEN(10,FILE='komega_original.dat',STATUS='UNKNOWN')
             DO i=1,dim_t
                k_t=8.D0*DATAN(1.D0)*REAL(i-dim_t/2-1)/lt
                write(10,*) k_t+omega_uppe,REAL(komega(i)),AIMAG(komega(i))
             ENDDO
             CLOSE(10)
             OPEN(10,FILE='komega_original_reduced.dat',STATUS='UNKNOWN')
             DO i=1,dim_t
                k_t=8.D0*DATAN(1.D0)*REAL(i-dim_t/2-1)/lt
                write(10,*) k_t+omega_uppe,REAL(komega(i)-rek0-rekp*(k_t+omega_uppe-omega)),AIMAG(komega(i))
             ENDDO
             CLOSE(10)
             mess1=.FALSE.
          ENDIF
          komega(j) = CMPLX(REAL(komega(j)), 0.D0, 8)
       endif
    ENDDO
    help1=0.D0
    DO j=1,dim_t
       k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
       help1=MAX(help1,ABS(REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega))))
    ENDDO
    PRINT*, 'Warning: maximal absolute reduced real part of k is ',help1
    PRINT*, 'Check komega_original.dat and komega_original_reduced.dat'
    IF (mess1) THEN
       OPEN(10,FILE='komega_original.dat',STATUS='UNKNOWN')
       DO i=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(i-dim_t/2-1)/lt
          write(10,*) k_t+omega_uppe,REAL(komega(i)),AIMAG(komega(i))
       ENDDO
       CLOSE(10)
       OPEN(10,FILE='komega_original_reduced.dat',STATUS='UNKNOWN')
       DO i=1,dim_t
          k_t=8.D0*DATAN(1.D0)*REAL(i-dim_t/2-1)/lt
          write(10,*) k_t+omega_uppe,REAL(komega(i)-rek0-rekp*(k_t+omega_uppe-omega)),AIMAG(komega(i))
       ENDDO
       CLOSE(10)
    ENDIF
    PRINT*, 'For artifical absorption type level of k, or 0 to ignore'
    READ(5,*) help1
    IF (help1.GT.0.D0) THEN
       PRINT*, 'absorption at ',help1,', see komega.dat and komega_reduced.dat'
       cut=.FALSE.
       DO j=1,dim_t 
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          if (( ABS(REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega))).GT.help1).AND.(cut.EQV..FALSE.)) then
             startcut=j
             cut=.TRUE.
          ENDIF
          if (( ABS(REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega))).LE.help1).AND.(cut.EQV..TRUE.)) then
             endcut=j-1
             cut=.FALSE.
             CALL artifabs(startcut,endcut,help1)
          endif
          if (( ABS(REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega))).GT.help1).AND.(cut.EQV..TRUE.).AND.(j.EQ.dim_t)) then
             endcut=dim_t
             CALL artifabs(startcut,endcut,help1)
          endif
       ENDDO
    ENDIF
    PRINT*, 'For cutting the reduced real part of k type level, or 0 to ignore'
    READ(5,*) help1
    IF (help1.GT.0.D0) THEN
       PRINT*, 'cutting at ',help1,', see komega.dat and komega_reduced.dat'
       cut=.FALSE.
       DO j=1,dim_t 
          k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
          if (( ABS(REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega))).GT.help1).AND.(cut.EQV..FALSE.)) then
             startcut=j
             cut=.TRUE.
          ENDIF
          if (( ABS(REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega))).LE.help1).AND.(cut.EQV..TRUE.)) then
             endcut=j-1
             cut=.FALSE.
             CALL artifdisp(startcut,endcut,help1)
          endif
          if (( ABS(REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega))).GT.help1).AND.(cut.EQV..TRUE.).AND.(j.EQ.dim_t)) then
             endcut=dim_t
             CALL artifdisp(startcut,endcut,help1)
          endif
       ENDDO
    ENDIF
    OPEN(10,FILE='komega.dat',STATUS='UNKNOWN')
    DO j=1,dim_t
       k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
       write(10,*) k_t+omega_uppe,REAL(komega(j)),AIMAG(komega(j))
    ENDDO
    CLOSE(10)
    OPEN(10,FILE='komega_reduced.dat',STATUS='UNKNOWN')
    DO j=1,dim_t
       k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
       write(10,*) k_t+omega_uppe,REAL(komega(j)-rek0-rekp*(k_t+omega_uppe-omega)),AIMAG(komega(j))
    ENDDO
    CLOSE(10)
    OPEN(10,FILE='nomega.dat',STATUS='UNKNOWN')
    DO j=1,dim_t
       k_t=8.D0*DATAN(1.D0)*REAL(j-dim_t/2-1)/lt
       write(10,*) k_t+omega_uppe,c*tp_fs_phys*1.d-15/(4.D0*z_rayleigh_cm_phys)*REAL(komega(j)/(k_t+omega_uppe)), &
            c*tp_fs_phys*1.d-15/(4.D0*z_rayleigh_cm_phys)*AIMAG(komega(j)/(k_t+omega_uppe))
    ENDDO
    CLOSE(10)

    IF (switch_T.EQ.3) THEN
       k_t=8.D0*DATAN(1.D0)/lt
       thg=3*NINT(omega/k_t)-(NINT(omega_uppe/k_t)-dim_t/2)+1
       fhg=5*NINT(omega/k_t)-(NINT(omega_uppe/k_t)-dim_t/2)+1
       IF (thg.GT.dim_t/2) THEN
          PRINT*, 'trusted spectral range does not include THG'
          deltak3omega(2)=0.D0
          deltak3omega(3)=0.D0
       ELSE      
          k_t=8.D0*DATAN(1.D0)*REAL(thg-dim_t/2-1)/lt+omega_uppe-omega
          deltak3omega(2)=REAL(komega(thg)) - rek0 - rekp * k_t
          deltak3omega(3)=2.D0*(rekp * omega - rek0)
       ENDIF
       IF (fhg.GT.dim_t/2) THEN
          PRINT*, 'trusted spectral range does not include FHG'
          deltak5omega(2)=0.D0
          deltak5omega(3)=0.D0
       ELSE
          k_t=8.D0*DATAN(1.D0)*REAL(fhg-dim_t/2-1)/lt+omega_uppe-omega
          deltak5omega(2)=REAL(komega(fhg)) - rek0 - rekp * k_t
          deltak5omega(3)=4.D0*(rekp * omega - rek0)
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE compute_dispersion

  COMPLEX(8) FUNCTION compute_n(omegab)
     IMPLICIT NONE

     REAL(8) omegab
     INTEGER(4) i

     IF(omegab.LT.omegachi(1)) THEN
        compute_n=chi(1)
     ELSE IF(omegab.GT.omegachi(dim_chi)) THEN
        compute_n=chi(dim_chi)
     ELSE
        DO i=2,dim_chi
           IF(omegachi(i).GE.omegab) THEN
              compute_n=((omegab-omegachi(i-1))*chi(i)+(omegachi(i)-omegab)*chi(i-1))/(omegachi(i)-omegachi(i-1))
              EXIT
           ENDIF
        ENDDO
     ENDIF
  END FUNCTION compute_n

  SUBROUTINE artifabs(startcut,endcut,help1)
    IMPLICIT NONE

    INTEGER(4) startcut,endcut,l
    REAL(8) help1

    IF(endcut-startcut.GT.500) THEN
       IF((startcut.EQ.1).AND.(ABS(REAL(komega(dim_t)-rek0-rekp*(8.D0*DATAN(1.D0)* &
         REAL(dim_t/2-1)/lt+omega_uppe-omega))).GT.help1)) THEN
          DO l=startcut,startcut+249
             komega(l) = komega(l)+CMPLX(0.D0,1.D6,8)
          ENDDO
       ELSE
          DO l=startcut,startcut+249
             komega(l) = komega(l)+CMPLX(0.D0,1.D6*exp(-(4.5D0*(REAL(l-startcut,8)-250.D0)/500.D0)**4),8)
          ENDDO
       ENDIF
       DO l=startcut+250,endcut-250
          komega(l) = komega(l)+CMPLX(0.D0,1.D6,8)
       ENDDO
       IF((endcut.EQ.dim_t).AND.(ABS(REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)* &
         REAL(-dim_t/2)/lt+omega_uppe-omega))).GT.help1)) THEN
          DO l=endcut-249,endcut
             komega(l) = komega(l)+CMPLX(0.D0,1.D6,8)
          ENDDO
       ELSE
          DO l=endcut-249,endcut
             komega(l) = komega(l)+CMPLX(0.D0,1.D6*exp(-(4.5D0*(REAL(l-endcut,8)+250.D0)/500.D0)**4),8)
          ENDDO
       ENDIF
    ELSE
       IF((startcut.EQ.1).AND.(ABS(REAL(komega(dim_t)-rek0-rekp*(8.D0*DATAN(1.D0)* &
         REAL(dim_t/2-1)/lt+omega_uppe-omega))).GT.help1)) THEN
          DO l=startcut,endcut
             komega(l) = komega(l)+CMPLX(0.D0,1.D6*exp(-(2.25D0*REAL(l-startcut,8)/REAL(endcut-startcut,8))**4),8)
          ENDDO
       ELSEIF ((endcut.EQ.dim_t).AND.(ABS(REAL(komega(1)-rek0-rekp*(8.D0*DATAN(1.D0)* &
           REAL(-dim_t/2)/lt+omega_uppe-omega))).GT.help1)) THEN
          DO l=startcut,endcut
             komega(l) = komega(l)+CMPLX(0.D0,1.D6*exp(-(2.25D0*(REAL(l-startcut,8)-REAL(endcut-startcut,8))/ &
               REAL(endcut-startcut,8))**4),8)
          ENDDO
       ELSE
          DO l=startcut,endcut
             komega(l) = komega(l)+CMPLX(0.D0,1.D6*exp(-(4.5D0*(REAL(l-startcut,8)-0.5D0*REAL(endcut-startcut,8))/ &
               REAL(endcut-startcut,8))**4),8)
          ENDDO
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE artifabs  
  
  SUBROUTINE artifdisp(startcut,endcut,help1)
    IMPLICIT NONE

    INTEGER(4) startcut,endcut,l
    REAL(8) help1,c1,c2,c3,c4,f1,fp1,f2,fp2

    if (startcut.EQ.1) then
       f1=help1
       f2=REAL(komega(endcut)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(endcut-dim_t/2-1)/lt+omega_uppe-omega))
       fp2=REAL(komega(endcut+1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(endcut+1-dim_t/2-1)/lt+omega_uppe-omega))-f2
       c3=(f1-f2+fp2*REAL(endcut-1,8))/(1.D0-REAL(endcut,8)**2+2.D0*REAL(endcut,8)*REAL(endcut-1,8))
       c2=fp2-2.D0*c3*REAL(endcut,8)
       c1=f1-c2-c3
       DO l=startcut,endcut
          komega(l) = CMPLX( c1+c2*REAL(l,8)+c3*REAL(l,8)**2 + rek0 + rekp*(8.D0*DATAN(1.D0)* &
            REAL(l-dim_t/2-1)/lt+omega_uppe-omega), &
               AIMAG(komega(l)) , 8)
       ENDDO
    elseif (endcut.EQ.dim_t) then
       f1=REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega))
       fp1=f1-REAL(komega(startcut-1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(startcut-1-dim_t/2-1)/lt+omega_uppe-omega))
       f2=help1
       c3=(f2-f1-fp1*REAL(endcut-startcut,8))/(REAL(endcut,8)**2-REAL(startcut,8)**2-2.D0*REAL(startcut,8)*REAL(endcut-startcut,8))
       c2=fp1-2.D0*c3*REAL(startcut,8)
       c1=f1-c2*REAL(startcut,8)-c3*REAL(startcut,8)**2
       DO l=startcut,endcut
          komega(l) = CMPLX( c1+c2*REAL(l,8)+c3*REAL(l,8)**2 + rek0 + rekp*(8.D0*DATAN(1.D0)* &
            REAL(l-dim_t/2-1)/lt+omega_uppe-omega), &
               AIMAG(komega(l)), 8)
       ENDDO
    else
       f1=REAL(komega(startcut)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(startcut-dim_t/2-1)/lt+omega_uppe-omega))
       fp1=f1-REAL(komega(startcut-1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(startcut-1-dim_t/2-1)/lt+omega_uppe-omega))
       f2=REAL(komega(endcut)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(endcut-dim_t/2-1)/lt+omega_uppe-omega))
       fp2=REAL(komega(endcut+1)-rek0-rekp*(8.D0*DATAN(1.D0)*REAL(endcut+1-dim_t/2-1)/lt+omega_uppe-omega))-f2
       c4=(-2.D0*f1+2.D0*f2-REAL(endcut,8)*fp2+REAL(startcut,8)*fp2+REAL(startcut,8)*fp1-REAL(endcut,8)*fp1) &
            /(-REAL(endcut,8)**3+3.D0*REAL(startcut,8)*REAL(endcut,8)**2+ &
            REAL(startcut,8)**3-3.D0*REAL(endcut,8)*REAL(startcut,8)**2)
       c3=(fp1-fp2-3.D0*c4*(REAL(startcut,8)**2-REAL(endcut,8)**2))/(2.D0*(REAL(startcut,8)-REAL(endcut,8)))
       c2=fp1-2.D0*c3*REAL(startcut,8)-3.D0*c4*REAL(startcut,8)**2
       c1=f1-c2*REAL(startcut,8)-c3*REAL(startcut,8)**2-c4*REAL(startcut,8)**3
       DO l=startcut,endcut
          komega(l) = CMPLX( c1+c2*REAL(l,8)+c3*REAL(l,8)**2 +c4*REAL(l,8)**3 &
               + rek0 + rekp*(8.D0*DATAN(1.D0)*REAL(l-dim_t/2-1)/lt+omega_uppe-omega), AIMAG(komega(l)) , 8)
       ENDDO
    endif

    RETURN
  END SUBROUTINE artifdisp

  SUBROUTINE compute_parameters
    IMPLICIT NONE

    INTEGER(4) KKnew
    REAL(8) c,h

    c = 3.d10 !speed of light in the vacuum      cm/s
    h = 6.634d-34/(2*3.1415) !Planck constant   J.s
    n2_phys = n2_phys*pressure
    n4_phys = n4_phys*pressure
    rhont_cm3_phys = rhont_cm3_phys*pressure

    tauc_fs_phys = tauc_fs_phys/pressure
    rhoc_cm3_phys = 1.11d13/lambda0_cm_phys**2 ! critical plasma censity  cm-3
    sigma_cm2_phys =  3.535d-12*(omega*tauc_fs_phys/tp_fs_phys)/(k0_phys*n0*(1.D0+(omega*tauc_fs_phys/tp_fs_phys)**2)) !cross section for inverse bremsstrahlung  cm2
    photon_energy_au_phys = h*c*2*3.1415D0/lambda0_cm_phys/4.359d-18 ! photon energy  au
    Ui_au_phys = Ui_eV_phys/27.2116d0 ! gap potential for ionization of oxygen molecules  au

    z_rayleigh_cm_phys  = 3.1415D0*w0_cm_phys**2*n0/lambda0_cm_phys  !Rayleigh lenght      cm
    KK = int(Ui_au_phys/photon_energy_au_phys)+1
    betak_phys =(KK*h*2*3.1415D0*c/(lambda0_cm_phys))*rhont_cm3_phys*sigmak_phys !coefficient of multiphoton absorption   cm2K-3/WK-1
    
    IF (switch_rho.LE.2) THEN
       PRINT*, 'From given ionization potential we get KK = ', KK, 'type alternative value or 0 to keep it'
       READ(5,*) KKnew
       IF (KKnew.GE.0) THEN
          KK=KKnew
          betak_phys =(Ui_au_phys/photon_energy_au_phys*h*2*3.1415D0*c/(lambda0_cm_phys))*rhont_cm3_phys*sigmak_phys
       ENDIF
    ENDIF
    
    Pcr_phys = ((lambda0_cm_phys)**2)/(2*3.1415*n0*n2_phys)  !critical power          W
    print *, 'aPcr ', 'n0', n0, 'n2', n2_phys
    k0_phys = 2.D0*3.1415D0/lambda0_cm_phys  !central wave number in vacuum     cm-1
    proplength = proplength_m_phys*100.D0/(4.D0*z_rayleigh_cm_phys)  !adimensionned distance of propagation
    outlength = outlength_m_phys*100.D0/(4.D0*z_rayleigh_cm_phys)  !adimmensionned output distance for whole field
    delta_z = delta_z_mm_phys/(40.D0*z_rayleigh_cm_phys)  !adimmensionned (first) stepwidth for whole field
    z=0.D0
    z_out=z
    rfil=0.1D0*rfil_mm_phys/w0_cm_phys ! adimensioned radius for diagnostics
    c3=1.d0  !adimensionned coefficient for kerr effect
    c5=(n4_phys)/(n2_phys)*(Pcr_phys)/(4*3.1415*w0_cm_phys**2)  !adimensionned coefficient for chi5 effect
    gamma1=sigma_cm2_phys*n0*rhoc_cm3_phys/k0_phys  !adimensionned coefficient for losses due to normalized conductivity
    gamma2=1.d0  !adimensionned coefficient for plasma defocusing
    muk=2.D0*z_rayleigh_cm_phys*betak_phys*(Pcr_phys/(4.d0*3.1415))**(KK-1)*w0_cm_phys**(2*(1-KK))  !adimensionned coefficient for MPA effect
    beta_inv_2KK=k0_phys**2*(sigmak_phys)*1.d-15*tp_fs_phys*(rhont_cm3_phys/rhoc_cm3_phys)* &
      (Pcr_phys/(4.d0*3.1415))**KK*w0_cm_phys**(2*(1-KK))  !adimensionned MPI coefficient
    eta1 = 2.d0 * z_rayleigh_cm_phys * sigman_phys * NN * ( h * 2.d0 *3.1415 *c/ lambda0_cm_phys)  * rhoabs_cm3_phys * & 
          (Pcr_phys/(4.d0*3.1415))**(NN-1)* w0_cm_phys**(2*(1-NN))                               ! adimensionned coefficient for absorption  
    eta2 = sigman_phys * 1.d-15*tp_fs_phys * ( Pcr_phys /(4.d0*3.1415 *w0_cm_phys**2))**NN       ! adimensionned coefficient for excited molecules 
    rho0=(rho0_phys)/(rhoc_cm3_phys*n0/(2.D0*z_rayleigh_cm_phys*k0_phys))  !adimensionned initial electron density
    nu=(sigma_cm2_phys*tp_fs_phys*1.d-15/(w0_cm_phys**2*Ui_eV_phys*1.6d-19))*(Pcr_phys/(4.d0*3.1415))  !adimensionned avalanche coefficient    
    alpha=alpha_fs_phys*tp_fs_phys !adimensionned linear recombination coefficient
    alphaquad=alphaquad_fscm3_phys*tp_fs_phys*(rhoc_cm3_phys*n0/(2*z_rayleigh_cm_phys*k0_phys)) !adimensionned quadratic recombination coefficient
    rhoat_inv=rhoc_cm3_phys*n0/(rhont_cm3_phys*2.d0*z_rayleigh_cm_phys*k0_phys) !adimensioned inverse density of neutral molekules
    tdk=tdk_fs_phys/tp_fs_phys  !adimensionned coefficient for delay time in kerr effect 
    raman=raman_phys*tp_fs_phys   ! adimensioned frequency for Raman response  
    if (f_cm_phys.EQ.0.D0 ) then
       lense_factor = 0.D0
    else
       lense_factor=3.1415D0*w0_cm_phys**2*n0/(f_cm_phys*lambda0_cm_phys)
    endif
    increase=decrease/2.5D0               !phase threshold for increasing delta_z, should be <= increase/2.5


    


    


    alpha1=alpha1_fs_phys*tp_fs_phys !adimensionned linear recombination for SLG1 electrons coefficient
    
    alphah=alphah_fs_phys*tp_fs_phys !adimensionned linear recombination for holes coefficient

    gamma1e=3.535d-12/(k0_phys*n0*omega)*n0*rhoc_cm3_phys/k0_phys  !factor for adimensionned coefficient for losses due to normalized conductivity




    density_normalisation_factor = (2.D0*z_rayleigh_cm_phys*k0_phys)/(rhoc_cm3_phys*n0)

    ALLOCATE(xx(i_x_max),zz(i_z_max),Indice_norm(i_x_max,i_z_max))
    xx = xx_mum*1.D-4/w0_cm_phys
    zz = zz_mum*1.D-4/(4.D0*z_rayleigh_cm_phys)
    Indice_norm = 4.D0*z_rayleigh_cm_phys*k0_phys*(Indice-n0)

    IF ((i_x_max.GT.1).OR.(i_z_max.GT.1)) THEN
       print*, 'The maximum deviation from background index on index array boundaries is ',MAX(MAXVAL(ABS(Indice(i_x_max,:)-n0)), &
         MAXVAL(ABS(Indice(:,1)-n0)),MAXVAL(ABS(Indice(:,i_z_max)-n0))),' hit enter to continue'
       READ(5,*)
    ENDIF
    
    RETURN
  END SUBROUTINE compute_parameters

END MODULE normalisation
