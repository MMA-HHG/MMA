MODULE write_start
  USE hdf5_helper_serial
  USE HDF5
  USE constants
  USE h5namelist

  ! REAL(8), PARAMETER  :: Pi = acos(-1.0d0)

  INTEGER(4) num_proc,dim_t,dim_r,KK,NN,switch_rho,switch_dKerr,absorb,rhodist,angular_momentum,switch_T
  INTEGER(4) i_x_max, i_z_max, i_x, i_z
  REAL(8) rek0,rekp,c3,c5,gamma1,gamma2,muk,beta_inv_2KK,rho0,nu,alpha,alphaquad,rhoat_inv,xdk,tdk,raman,omega,eta1,eta2

  REAL(8) omega_uppe
  REAL(8) lt,lr,proplength,outlength,outlength_Efield,delta_z,z,z_out,z_out_Efield,rfil,increase,decrease,time_limit
  LOGICAL out_Efield
  REAL(8) photon_energy_au_phys,tp_fs_phys,Pcr_phys,w0_m_phys,w0_cm_phys
  REAL(8) Ui_au_phys,residue_charge,n0,rhoc_cm3_phys,rhont_cm3_phys

  REAL(8) z_rayleigh_cm_phys

  REAL(8) ions_Kerr_ratio

  REAL(8) density_normalisation_factor
  REAL(8) delta_t, tlo 
  REAL(8), ALLOCATABLE :: xx(:),zz(:),Indice_norm(:,:),real_e(:,:),imag_e(:,:)
  COMPLEX(8), ALLOCATABLE :: e(:,:),e_full(:,:),komega(:)
  INTEGER(HID_T) :: file_id, group_id, group_id2
  INTEGER :: error
  CHARACTER(LEN = *), PARAMETER :: output_groupname = pre_proc_grpname
  CHARACTER(100) :: filename  ! File name variable, which is assigned a value by user eg.: results.h5
  INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims       
  
CONTAINS

  SUBROUTINE write_startingfile(p)
    IMPLICIT NONE
    INTEGER(HSIZE_T) :: r_offset
    INTEGER(4) p,j,l,k,k1,k2
    CHARACTER(LEN = 3) :: ip
    CHARACTER(LEN = 10):: iz,id
    REAL(8), DIMENSION(1:2):: test
    REAL(8) :: efield_factor ! normalization factor electric field V/m
    COMPLEX(8), ALLOCATABLE  :: efield_osc(:) ! fast oscillating term exp(-i*omegauppe*t)
    INTEGER(HID_T) :: indexes_group_id
    CHARACTER(LEN=35) :: indexes_groupname="/pre-processed/indexes_group"    
    INTEGER(HID_T) :: tables_file_id
    CHARACTER(LEN=25) :: tables_filename="calculated_tables.h5"
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
    IF(p.EQ.(num_proc-1))THEN
      PRINT*, 'final writing'
      CALL create_dset(group_id, 'num_proc', num_proc)
      CALL create_dset(group_id, 'dim_t', dim_t)
      CALL create_dset(group_id, 'dim_r', dim_r) 
      CALL create_dset(group_id, 'rek0', rek0)
      CALL create_dset(group_id, 'rekp', rekp)
      CALL create_dset(group_id, 'c3', c3)
      CALL create_dset(group_id, 'c5', c5)
      CALL create_dset(group_id, 'gamma1', gamma1)
      CALL create_dset(group_id, 'gamma2', gamma2)
      CALL create_dset(group_id, 'muk', muk)
      CALL create_dset(group_id, 'beta_inv_2KK', beta_inv_2KK)
      CALL create_dset(group_id,'KK',KK)
      CALL create_dset(group_id,'rho0',rho0)
      CALL create_dset(group_id,'nu',nu)
      CALL create_dset(group_id,'alpha',alpha)
      CALL create_dset(group_id,'alphaquad', alphaquad)
      CALL create_dset(group_id,'rhoat_inv', rhoat_inv)
      CALL create_dset(group_id,'xdk', xdk)
      CALL create_dset(group_id,'tdk', tdk)
      CALL create_dset(group_id,'raman', raman)
      CALL create_dset(group_id,'omega', omega)
      CALL create_dset(group_id,'komega', komega(1:dim_t), dim_t)
      CALL create_dset(group_id,'NN', NN)
      CALL create_dset(group_id,'eta1', eta1)
      CALL create_dset(group_id,'eta2', eta2)
      CALL create_dset(group_id,'lt', lt)
      CALL create_dset(group_id,'lr', lr)
      CALL create_dset(group_id,'proplength', proplength)
      CALL create_dset(group_id,'outlength', outlength)
      CALL create_dset(group_id,'delta_z', delta_z)
      CALL create_dset(group_id,'z', z)
      CALL create_dset(group_id,'z_out', z_out)
      CALL create_dset(group_id,'rfil', rfil)
      CALL create_dset(group_id,'switch_rho', switch_rho)
      CALL create_dset(group_id,'switchKerr', switch_dKerr)
      CALL create_dset(group_id,'switch_T', switch_T)
      CALL create_dset(group_id,'absorb', absorb)
      CALL create_dset(group_id,'increase', increase)
      CALL create_dset(group_id,'decrease', decrease)
      CALL create_dset(group_id,'rhodist', rhodist)
      CALL create_dset(group_id,'timelimit', time_limit)
      CALL create_dset(group_id,'photenergy',photon_energy_au_phys)
      CALL create_dset(group_id,'pulsedurat',tp_fs_phys)
      CALL create_dset(group_id,'critpower',Pcr_phys*1.D-9)
      CALL create_dset(group_id,'beam_waist',w0_cm_phys)
      CALL create_dset(group_id,'ionpot',Ui_au_phys)
      CALL create_dset(group_id,'rescharge',residue_charge)
      CALL create_dset(group_id,'n0_indice',n0)
      CALL create_dset(group_id,'critdens',rhoc_cm3_phys)
      CALL create_dset(group_id,'atomdens',rhont_cm3_phys)

      CALL create_dset(group_id,'angmom',angular_momentum)

      CALL create_dset(group_id,'finished',.FALSE.)
      CALL create_dset(group_id,'omega_uppe', omega_uppe)


      CALL create_dset(group_id,'out_Efield', out_Efield)
      IF (out_Efield) THEN
        CALL create_dset(group_id,'z_out_Efield', z_out_Efield)
        CALL create_dset(group_id,'outlength_Efield', outlength_Efield)
      ENDIF

      

      CALL create_dset(group_id,'density_normalisation_factor',density_normalisation_factor)

      CALL create_dset(group_id,'four_z_rayleigh_cm_phys', 4.D0*z_rayleigh_cm_phys)

      !r_offset = dim_r_start(num_proc)-1
      efield_factor = SQRT(Pcr_phys*1.D-9*1.D9*c_light*4.D0*PI*1.D-7/(4.D0*PI*w0_cm_phys**2*1.D-4*2.D0*n0))*2.D0 ! normalization factor electric field V/m
      print *, 'efield_factor', efield_factor, 'w0', w0_cm_phys, 'Pcr', Pcr_phys

      print *, 'old efield fact', SQRT(Pcr_phys*1.D-9*1.D9*3.D8*4.D0*PI*1.D-7/(4.D0*PI*w0_cm_phys**2*1.D-4*2.D0*n0))*2.D0
      print *, 'new efield fact', SQRT(Pcr_phys*1.D-9*1.D9*c_light*4.D0*PI*1.D-7/(4.D0*PI*w0_cm_phys**2*1.D-4*2.D0*n0))*2.D0


      ALLOCATE(efield_osc(dim_t))
      PRINT*, 'beosc'
      DO j=1,dim_t
        efield_osc(j) = exp(CMPLX(0.D0,-omega_uppe*(tlo+REAL(j,8)*delta_t),8)) ! fast oscillating term exp(-i*omegauppe*t)
      ENDDO
      PRINT*, 'aeosc'
      
      ALLOCATE(real_e(dim_t,dim_r),imag_e(dim_t,dim_r))
      PRINT*, 'brealimag'
      DO k1=1, dim_t
        DO k2=1, dim_r
           real_e(k1,k2) = REAL( REAL( (efield_factor*efield_osc(k1)*e_full(k1,k2)) ) , 8 )
           imag_e(k1,k2) = REAL( AIMAG( (efield_factor*efield_osc(k1)*e_full(k1,k2)) ) , 8 )
           !real_e(k2,k1) = REAL( REAL( (e_full(k1,k2)) ) , 8 )
           !imag_e(k2,k1) = REAL( IMAG( (e_full(k1,k2)) ) , 8 )
        ENDDO
      ENDDO
      CALL create_dset(group_id, "startfield_r", real_e, dim_t, dim_r)
      PRINT*, 'efieldwritten'
      CALL create_dset(group_id, "startfield_i", imag_e, dim_t, dim_r)
      ! CALL create_dset(group_id,"startfield",e_full,dim_t,dim_r)
      PRINT*, 'bindicesgroup'
      CALL h5gcreate_f(file_id, indexes_groupname, indexes_group_id, error)
      CALL create_dset(indexes_group_id, "r_vector", REAL(xx(1:i_x_max),4), i_x_max)
      CALL create_dset(indexes_group_id, "z_vector", REAL(zz(1:i_z_max),4), i_z_max)
      CALL create_2D_array_real_dset(indexes_group_id, "indexes", REAL(Indice_norm(1:i_x_max, 1:i_z_max),8), i_x_max, i_z_max)
      CALL h5gclose_f(indexes_group_id, error)
      ! Open file calculated_tables.h5
      !CALL h5fopen_f(tables_filename, H5F_ACC_RDWR_F, tables_file_id, error)
      !CALL h5gcreate_f(tables_file_id, "indexes", indexes_group_id, error)
      !CALL create_dset(indexes_group_id, "r_vector", REAL(xx(1:i_x_max),4), i_x_max)
      !CALL create_dset(indexes_group_id, "z_vector", REAL(zz(1:i_z_max),4), i_z_max)
      !CALL create_2D_array_real_dset(indexes_group_id, "indexes", REAL(Indice_norm(1:i_x_max, 1:i_z_max),8), i_x_max, i_z_max)
      !CALL h5gclose_f(indexes_group_id, error)
      !CALL h5fclose_f(tables_file_id, error)
    ENDIF
    DEALLOCATE(e)

920 FORMAT (F10.6)
930 FORMAT (I3)

    RETURN
  END SUBROUTINE write_startingfile

END MODULE write_start
