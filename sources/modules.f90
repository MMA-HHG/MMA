MODULE fields
  INTEGER(4)               :: dim_t,dim_r,dim_t_local,dim_r_local,omega_offset(2), i_x_max, i_z_max
  COMPLEX(8) phaserelation(3)
  INTEGER(4), ALLOCATABLE  :: dim_t_start(:),dim_t_end(:),dim_r_start(:),dim_r_end(:),num_ex(:)
  INTEGER(4), ALLOCATABLE  :: send_e(:),send_etemp(:),send_efft(:),recv_e(:),recv_etemp(:),recv_efft(:)
  REAL(8) , ALLOCATABLE    :: bound_t(:),e_2(:),e_2KK(:),e_2KKm2(:),rho(:),fluence(:),rhoabs(:),losses_ionization(:),losses_plasma(:)
  REAL(8), ALLOCATABLE     :: peakmax(:),rhomax(:),energy(:),z_buff(:),energy_fil(:), rhoabs_max(:), xx(:), zz(:), Indice_norm(:,:),rhoO2max(:),rhoN2max(:),Tevmax(:)
  COMPLEX(8), ALLOCATABLE  :: e(:,:),etemp(:,:),efft(:,:),ptemp(:,:),jtemp(:,:),D(:,:),DL(:,:),DU(:,:),p_t(:),op_t(:),op_t_inv(:),pharm(:,:),hfac(:,:)
  LOGICAL, ALLOCATABLE     :: send_first(:)
END MODULE fields

MODULE parameters
  REAL(8) :: rek0,rekp,c3,c5,gamma1,gamma2,muk,beta_inv_2KK,omega, eta1, eta2, omega_uppe
  REAL(8) :: beta_inv_2KKp,eti_ref,exp_ref,beta_inv_2,mukp,mu,mukpp,beta_inv_2KKpp
  REAL(8) :: gamma1e,nuO2,nuN2,T_init_eV_phys,nukB,nucp,nucO2,nucN2,rhoat_N2_inv
  COMPLEX(8), ALLOCATABLE :: komega(:),komega_red(:)
  INTEGER(4) :: KK,KKp,KKpp
  INTEGER(4) :: NN
  REAL(8) :: rho0
  REAL(8) :: nu,alpha,alphaquad,rhoat_inv
  REAL(8) :: alpha1,alpha2,alphah,rhosat
  REAL(8)    :: xdk,tdk,raman, expt1,expt2,expt3,expt4,expt1p,expt2p,expt3p,expt4p,c3i,c3d
  INTEGER(4) :: switch_dKerr
  REAL(8) timelimit
  INTEGER(4) :: dim_th,count,i_x_old,i_z_old
  REAL(8)    :: tlo,lt,lr
  REAL(8)    :: proplength,outlength,z,delta_z,z_out,increase,decrease,delta_zh,delta_t,delta_r,delta_t_inv,delta_z_max,rfil
  REAL(8)    :: k_t
  COMPLEX(8), ALLOCATABLE :: delta_rel(:)
  INTEGER(4) :: switch_rho,absorb,rhodist,switch_T
  REAL(8)    :: maxphase
  INTEGER, parameter :: unit_peakmax=7,unit_rho=8,unit_logfile=9,unit_field=10,unit_energy=11,unit_rhomax=12, unit_rhoabs_max = 13
  REAL(8)    :: ionisation_potential_N2,residue_charge_N2,atomic_density_N2
  INTEGER(4) :: angular_momentum_N2
  INTEGER(4) :: HDF5write_count
END MODULE parameters

MODULE mpi_stuff
  include 'mpif.h'
  INTEGER(4) my_rank,num_proc,ierr,MPI_SUBARRAY,MPI_SUBARRAY_TRANSPOSED
  INTEGER(4) status(MPI_STATUS_SIZE)
  CHARACTER(3) ip
END MODULE mpi_Stuff

MODULE run_status
  LOGICAL finished
END MODULE run_status

