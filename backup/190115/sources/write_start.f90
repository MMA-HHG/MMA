MODULE write_start

  INTEGER(4) num_proc,dim_t,dim_r,KK,NN,switch_rho,switch_dKerr,absorb,rhodist,angular_momentum,switch_T,KKp,KKpp,angular_momentum_N2
  INTEGER(4) i_x_max, i_z_max, i_x, i_z
  REAL(8) rek0,rekp,c3,c5,gamma1,gamma2,muk,beta_inv_2KK,rho0,nu,alpha,alphaquad,rhoat_inv,xdk,tdk,raman,omega,eta1,eta2
  REAL(8) beta_inv_2KKp,eti_ref,exp_ref,beta_inv_2,mukp,mu,mukpp,beta_inv_2KKpp,alpha1,alpha2,alphah,rhosat,gamma1e,nuO2,nuN2,T_init_eV_phys,nukB
  REAL(8) omega_uppe,nucp,nucO2,nucN2,rhoat_N2_inv
  REAL(8) lt,lr,proplength,outlength,delta_z,z,z_out,rfil,increase,decrease,time_limit
  REAL(8) photon_energy_au_phys,tp_fs_phys,Pcr_phys,w0_cm_phys
  REAL(8) Ui_au_phys,residue_charge,n0,rhoc_cm3_phys,rhont_cm3_phys,reduced_mass
  REAL(8) Ui_au_phys_N2,residue_charge_N2,rhont_N2_cm3_phys
  REAL(8), ALLOCATABLE :: xx(:),zz(:),Indice_norm(:,:)
  COMPLEX(8), ALLOCATABLE :: e(:,:),komega(:)

CONTAINS

  SUBROUTINE write_startingfile(p)
    IMPLICIT NONE

    INTEGER(4) p,j,l,k
    CHARACTER(LEN = 3) :: ip
    CHARACTER(LEN = 10):: iz,id

    WRITE(iz,920) z
    DO k=1,10
       IF (iz(k:k).EQ.' ') iz(k:k)='0'
       IF (iz(k:k).EQ.'.') iz(k:k)='_'
    ENDDO

    IF (p.EQ.0) THEN    
       OPEN(10,FILE='PROP_RAD.LOG',STATUS='NEW')
       WRITE(10,*) iz
       CLOSE(10)
       OPEN(10,FILE='STOP',STATUS='NEW')
       CLOSE(10)
       OPEN(10,FILE='MERGE_RAD.LOG',STATUS='NEW')
       CLOSE(10)
    ENDIF

    WRITE(ip,930) p
    DO l=1,3
       IF (ip(l:l).EQ.' ') ip(l:l)='0'
    ENDDO
    OPEN(10,FILE=iz//'_'//ip//'.DAT',STATUS='NEW',FORM='UNFORMATTED')
    PRINT*, 'write starting file for proc ',p
    id='num_proc'
    WRITE(10) id,num_proc
    id='dim_t'
    WRITE(10) id,dim_t
    id='dim_r'
    WRITE(10) id,dim_r
    id='rek0'
    WRITE(10) id,rek0
    id='rekp'
    WRITE(10) id,rekp
    id='c3'
    WRITE(10) id,c3
    id='c5'
    WRITE(10) id,c5
    id='gamma1'
    WRITE(10) id,gamma1
    id='gamma2'
    WRITE(10) id,gamma2
    id='muk'
    WRITE(10) id,muk
    id='betainv2KK'
    WRITE(10) id,beta_inv_2KK
    id='KK'
    WRITE(10) id,KK
    id='rho0'
    WRITE(10) id,rho0
    id='nu'
    WRITE(10) id,nu
    id='alpha'
    WRITE(10) id,alpha
    id='alphaquad'
    WRITE(10) id,alphaquad
    id='rhoat_inv'
    WRITE(10) id,rhoat_inv
    id='xdk'
    WRITE(10) id,xdk
    id='tdk'
    WRITE(10) id,tdk
    id='raman'
    WRITE(10) id,raman
    id='omega'
    WRITE(10) id,omega
    id='komega'
    WRITE(10) id,komega(1:dim_t)
    id='NN'
    WRITE(10) id,NN
    id='eta1'
    WRITE(10) id,eta1
    id='eta2'
    WRITE(10) id,eta2
    id='lt'
    WRITE(10) id,lt
    id='lr'
    WRITE(10) id,lr
    id='proplength'
    WRITE(10) id,proplength
    id='outlength'
    WRITE(10) id,outlength
    id='delta_z'
    WRITE(10) id,delta_z
    id='z'
    WRITE(10) id,z
    id='z_out'
    WRITE(10) id,z_out
    id='rfil'
    WRITE(10) id,rfil
    id='switch_rho'
    WRITE(10) id,switch_rho
    id='switchKerr'
    WRITE(10) id,switch_dKerr
    id='switch_T'
    WRITE(10) id,switch_T
    id='absorb'
    WRITE(10) id,absorb
    id='increase'
    WRITE(10) id,increase
    id='decrease'
    WRITE(10) id,decrease
    id='rhodist'
    WRITE(10) id,rhodist
    id='timelimit'
    WRITE(10) id,time_limit
    id='photenergy'
    WRITE(10) id,photon_energy_au_phys
    id='pulsedurat'
    WRITE(10) id,tp_fs_phys
    id='critpower'
    WRITE(10) id,Pcr_phys*1.D-9
    id='beam_waist'
    WRITE(10) id,w0_cm_phys
    id='ionpot'
    WRITE(10) id,Ui_au_phys
    id='rescharge'
    WRITE(10) id,residue_charge
    id='n0_indice'
    WRITE(10) id,n0
    id='critdens'
    WRITE(10) id,rhoc_cm3_phys
    id='atomdens'
    WRITE(10) id,rhont_cm3_phys
    id='reducmass'
    WRITE(10) id,reduced_mass
    id='angmom'
    WRITE(10) id,angular_momentum
    id='KKp'
    WRITE(10) id,KKp
    id='beta_inv_2KKp'
    WRITE(10) id,beta_inv_2KKp
    id='mukp'
    WRITE(10) id,mukp
    id='beta_inv_2'
    WRITE(10) id,beta_inv_2
    id='mu'
    WRITE(10) id,mu
    id='KKpp'
    WRITE(10) id,KKpp
    id='beta_inv_2KKpp'
    WRITE(10) id,beta_inv_2KKpp
    id='mukpp'
    WRITE(10) id,mukpp
    id='eti_ref'
    WRITE(10) id,eti_ref
    id='exp_ref'
    WRITE(10) id,exp_ref
    id='alpha1'
    WRITE(10) id,alpha1
    id='alpha2'
    WRITE(10) id,alpha2
    id='alphah'
    WRITE(10) id,alphah
    id='rhosat'
    WRITE(10) id,rhosat
    id='finished'
    WRITE(10) id,.FALSE.
    id='omega_uppe'
    WRITE(10) id,omega_uppe
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
    WRITE(10) id,Ui_au_phys_N2
    id='rescharge_N2'
    WRITE(10) id,residue_charge_N2
    id='atomdens_N2'
    WRITE(10) id,rhont_N2_cm3_phys
    id='angmom_N2'
    WRITE(10) id,angular_momentum_N2
    id='startfield'
    WRITE(10) id
    DO j=p*(dim_r/num_proc)+1,(p+1)*(dim_r/num_proc)
       WRITE(10) e(1:dim_t,j)
    ENDDO
    id='index'
    WRITE(10) id
    WRITE(10) i_x_max, i_z_max
    WRITE(10) (xx(i_x),i_x=1,i_x_max)
    DO i_z = 1, i_z_max
       WRITE(10) zz(i_z)
       WRITE(10) (Indice_norm(i_x,i_z),i_x=1,i_x_max)
    ENDDO  
    CLOSE(10)
    DEALLOCATE(e)

920 FORMAT (F10.6)
930 FORMAT (I3)

    RETURN
  END SUBROUTINE write_startingfile

END MODULE write_start
