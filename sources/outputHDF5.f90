MODULE output
  USE fields
  USE parameters
  USE mpi_stuff
  USE long_step
  USE run_status
CONTAINS

  SUBROUTINE HDF5_out
    USE fft
    USE HDF5
    IMPLICIT NONE

    INTEGER(4) j,k,l
    REAL(8) rhotemp,r,mpa
    COMPLEX(8) help
    CHARACTER*10 iz,filename




	!
	CHARACTER(LEN=10), PARAMETER :: filename = "results.h5"

	!
	!dataset rank is 2 and name is "ExtendibleArray"
	!
	CHARACTER(LEN=15), PARAMETER :: dsetname = "micro/FieldsForTDSEreal"
	CHARACTER(LEN=15), PARAMETER :: dsetname = "micro/FieldsForTDSEimag"
	INTEGER :: RANK = 2

	INTEGER(HID_T) :: file_id       ! File identifier 
	INTEGER(HID_T) :: dset_id       ! Dataset identifier 
	INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
	INTEGER(HID_T) :: memspace      ! Memory dataspace identifier 
	INTEGER(HID_T) :: h5parameters      ! Dataset creation property identifier 

    !!! in the first run, create extendible dataset and fill data
	field_dimensions = 3;
	IF (!FIRST RUN)


	!Initialize HDF5
	CALL h5open_f(error) 

	!define parameters of HDF5 workflow for MPI-access
	CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5parameters, error) ! create access parameters
        CALL h5pset_dxpl_mpio_f(h5parameters, H5FD_MPIO_COLLECTIVE_F, error) ! allow MPI access (should it be here?)

	!Open collectivelly the file
	CALL h5fopen_f(filename, H5F_ACC_RDWR, file_id, error, access_prp = h5parameters ) ! open file collectivelly
	CALL h5pclose(h5parameters) ! parameters were used for MPI open, close them
	

	!Create the dataspace with unlimited dimension in z. ! again, what should I use for parallel access?
	maxdims = (/H5S_UNLIMITED_F, dim_r, dim_t/)
	dims = (/1,dim_r,dim_t/)
	CALL h5screate_simple_f(field_dimensions, dims, dataspace, error, maxdims) ! create dataspace, the dataspace is its identifier !!!(global dataspace for all the data)
	dims = (1, dim_r_end(num_proc)-dim_r_start(num_proc), dim_t) ! dimension of my field
	CALL h5screate_simple_f(field_dimensions, dims, dataspacelocal, error, maxdims) ! create dataspace, the dataspace is its identifier !!!(local dataspace for local data)


	! we create the dataset collectivelly
	CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_FLOAT, filespace, dset_id, error) ! according to smilei, propabably optional args 

	!we use hyperslab to assign part of the global dataset
	!chunk data for each worker
	offset = (/0,dim_r_start(num_proc),0/)
	ccount = (/1, dim_r_end(num_proc) - dim_r_start(num_proc) , dim_t/)
	CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error) ! we should have access to its part of the dataset for each worker

	!Finally, write data
	CALL h5pcreate_f(H5P_DATASET_XFER_F, h5parameters, error)
	CALL h5pset_dxpl_mpio_f(h5parameters, H5FD_MPIO_COLLECTIVE_F, error) ! collective writting
!	CALL h5pset_chunk_f(crp_list, field_dimensions, dimsc, error) ???????????? Do we need chunk it?
	CALL h5dwrite_vl_f(dset_id , H5T_NATIVE_FLOAT, Fields, dimsfi, error, memspace,file_space_id=filespace,xfer_prp = h5parameters)! data are written !!!( probably variable length)
	CALL h5pclose(h5parameters) ! parameters were used for MPI open, close them
	
	DO k1=1,dim_t	
	DO k2=1,dim_t
		Fields(1,k1,k2) = REAL(e(k1,k2));
	ENDDO
	ENDDO



	
	CALL h5close_f(error) ! close the HDF5 workspace
   
    !!! now, just append data 

	ELSE (!ANOTHER RUN)

	!Open file for reading and writing ! ? WHAT SHOULD I USE FOR PARALLEL
	CALL h5fopen_f(filename, H5F_ACC_RDWR file_id, error)


	! open the dataset
  	CALL h5dopen_f(file_id, dsetname, dset_id, error) ! Open an existing dataset.

	! recognise the dataset size now
        CALL h5sget_simple_extent_dims(dataspace,dims,maxdims,error)
	zdim = dims(1)
	zdim = zdim + 1

	newsize = (/zdim, dim_r, dim_t/)

	!extend the dataset
	CALL h5dset_extent_f(dset_id, newsize, error)


	! Now I use hyperslab to select proper plane
!	offset = (/zdim,0,0/)
!	ccount = (/1,dim_r,dim_t/)
	!chunk data for each worker
	offset = (/zdim,dim_r_start(num_proc),0/)
	ccount = (/1, dim_r_end(num_proc) - dim_r_start(num_proc) , dim_t/)

	data_dims = (/1,dim_r,dim_t/)
	CALL h5dget_space_f(dset_id, dataspace, error)
	CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

	do k1 = dim_r_start(num_proc),dim_r_end(num_proc)
		do k2 = 1, dim_t
			data_append(k1,k2) = REAL(e(k2,k1))
		enddo
	enddo

	!Write data to the dataset - hyperslab next plane
	CALL H5dwrite_f(dset_id, H5T_NATIVE_FLOAT, data_append, data_dims, error, memspace, dataspace)





	ENDIF

!       OPEN(unit_field,FILE=iz//'_FIELD_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
!       WRITE(unit_field) dim_t,dim_r,num_proc
!       WRITE(unit_field) REAL(delta_t,4),REAL(delta_r,4),REAL(tlo,4)
!       DO j=dim_r_start(num_proc),dim_r_end(num_proc)
!          WRITE(unit_field) CMPLX(e(1:dim_t,j),KIND=4)
!       ENDDO
!       CLOSE(unit_field)


    RETURN
  END SUBROUTINE  HDF5_out

!  SUBROUTINE  field_out
!    USE ppt
!    IMPLICIT  NONE

!    INTEGER(4) j,k,i_x,i_z
!    CHARACTER*10  iz,id,filename

!    WRITE(iz,920) z
!    DO  k=1,10
!       IF (iz(k:k).EQ.' ') iz(k:k)='0'
!       IF (iz(k:k).EQ.'.') iz(k:k)='_'
!    ENDDO
!    IF (my_rank.EQ.0) THEN
!       filename='non'
!       OPEN(unit_logfile,FILE='PROP_RAD.LOG',STATUS='UNKNOWN')
!       DO
!          READ(unit_logfile,*,END=999) filename
!       ENDDO
!999    CONTINUE
!       CLOSE(unit_logfile)
!    ENDIF
!    CALL MPI_BCAST(filename,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!    IF (filename.NE.iz) THEN
!       IF(my_rank.EQ.0) THEN
!          OPEN(unit_logfile,FILE='PROP_RAD.LOG',STATUS='UNKNOWN',POSITION='APPEND')
!          WRITE(unit_logfile,*) iz    
!          CLOSE(unit_logfile)
!       ENDIF
!       OPEN(unit_field,FILE=iz//'_'//ip//'.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
!       id='num_proc'
!       WRITE(unit_field) id,num_proc
!       id='dim_t'
!       WRITE(unit_field) id,dim_t
!       id='dim_r'
!       WRITE(unit_field) id,dim_r
!       id='rek0'
!       WRITE(unit_field) id,rek0
!       id='rekp'
!       WRITE(unit_field) id,rekp
!       id='c3'
!       WRITE(unit_field) id,c3
!       id='c5'
!       WRITE(unit_field) id,c5
!       id='gamma1'
!       WRITE(unit_field) id,gamma1
!       id='gamma2'
!       WRITE(unit_field) id,gamma2
!       id='muk'
!       WRITE(unit_field) id,muk
!       id='betainv2KK'
!       WRITE(unit_field) id,beta_inv_2KK
!       id='KK'
!       WRITE(unit_field) id,KK
!       id='rho0'
!       WRITE(unit_field) id,rho0
!       id='nu'
!       WRITE(unit_field) id,nu
!       id='alpha'
!       WRITE(unit_field) id,alpha
!       id='alphaquad'
!       WRITE(unit_field) id,alphaquad
!       id='rhoat_inv'
!       WRITE(unit_field) id,rhoat_inv
!       id='xdk'
!       WRITE(unit_field) id,xdk
!       id='tdk'
!       WRITE(unit_field) id,tdk
!       id='raman'
!       WRITE(unit_field) id,raman
!       id='omega'
!       WRITE(unit_field) id,omega
!       id='komega'
!       WRITE(unit_field) id,komega(1:dim_t)
!       id='NN'
!       WRITE(unit_field) id,NN
!       id='eta1'
!       WRITE(unit_field) id,eta1
!       id='eta2'
!       WRITE(unit_field) id,eta2
!       id='lt'
!       WRITE(unit_field) id,lt
!       id='lr'
!       WRITE(unit_field) id,lr
!       id='proplength'
!       WRITE(unit_field) id,proplength
!       id='outlength'
!       WRITE(unit_field) id,outlength
!       id='delta_z'
!       WRITE(unit_field) id,delta_z
!       id='z'
!       WRITE(unit_field) id,z
!       id='z_out'
!       WRITE(unit_field) id,z_out
!       id='rfil'
!       WRITE(unit_field) id,rfil
!       id='switch_rho'
!       WRITE(unit_field) id,switch_rho
!       id='switchKerr'
!       WRITE(unit_field) id,switch_dKerr
!       id='switch_T'
!       WRITE(unit_field) id,switch_T
!       id='absorb'
!       WRITE(unit_field) id,absorb
!       id='increase'
!       WRITE(unit_field) id,increase
!       id='decrease'
!       WRITE(unit_field) id,decrease
!       id='rhodist'
!       WRITE(unit_field) id,rhodist
!       id='timelimit'
!       WRITE(unit_field) id,timelimit
!       id='photenergy'
!       WRITE(unit_field) id,photon_energy
!       id='pulsedurat'
!       WRITE(unit_field) id,pulse_duration
!       id='critpower'
!       WRITE(unit_field) id,critical_power
!       id='beam_waist'
!       WRITE(unit_field) id,beam_waist
!       id='ionpot'
!       WRITE(unit_field) id,ionisation_potential
!       id='rescharge'
!       WRITE(unit_field) id,residue_charge
!       id='n0_indice'
!       WRITE(unit_field) id,n0_indice
!       id='critdens'
!       WRITE(unit_field) id,critical_density
!       id='atomdens'
!       WRITE(unit_field) id,atomic_density
!       id='reducmass'
!       WRITE(unit_field) id,reduced_mass
!       id='angmom'
!       WRITE(unit_field) id,angular_momentum
!       id='KKp'
!       WRITE(unit_field) id,KKp
!       id='beta_inv_2KKp'
!       WRITE(unit_field) id,beta_inv_2KKp
!       id='mukp'
!       WRITE(unit_field) id,mukp
!       id='beta_inv_2'
!       WRITE(unit_field) id,beta_inv_2
!       id='mu'
!       WRITE(unit_field) id,mu
!       id='KKpp'
!       WRITE(unit_field) id,KKpp
!       id='beta_inv_2KKpp'
!       WRITE(unit_field) id,beta_inv_2KKpp
!       id='mukpp'
!       WRITE(unit_field) id,mukpp
!       id='eti_ref'
!       WRITE(unit_field) id,eti_ref
!       id='exp_ref'
!       WRITE(unit_field) id,exp_ref
!       id='alpha1'
!       WRITE(unit_field) id,alpha1
!       id='alpha2'
!       WRITE(unit_field) id,alpha2
!       id='alphah'
!       WRITE(unit_field) id,alphah
!       id='rhosat'
!       WRITE(unit_field) id,rhosat
!       id='finished'
!       WRITE(unit_field) id,finished
!       id='omega_uppe'
!       WRITE(unit_field) id,omega_uppe
!       id='gamma1e'
!       WRITE(10) id,gamma1e
!       id='nuO2'
!       WRITE(10) id,nuO2
!       id='nuN2'
!       WRITE(10) id,nuN2
!       id='T_init_eV_phys'
!       WRITE(10) id,T_init_eV_phys
!       id='nukB'
!       WRITE(10) id,nukB
!       id='nucp'
!       WRITE(10) id,nucp
!       id='nucO2'
!       WRITE(10) id,nucO2
!       id='nucN2'
!       WRITE(10) id,nucN2
!       id='rhoat_N2_inv'
!       WRITE(10) id,rhoat_N2_inv
!       id='ionpotN2'
!       WRITE(10) id,ionisation_potential_N2
!       id='rescharge_N2'
!       WRITE(10) id,residue_charge_N2
!       id='atomdens_N2'
!       WRITE(10) id,atomic_density_N2
!       id='angmom_N2'
!       WRITE(10) id,angular_momentum_N2
!       id='startfield'
!       WRITE(unit_field) id
!       DO j=dim_r_start(num_proc),dim_r_end(num_proc)
!          WRITE(unit_field) e(1:dim_t,j)
!       ENDDO
!       id='index'
!       WRITE(unit_field) id
!       WRITE(unit_field) i_x_max, i_z_max
!       WRITE(unit_field) (xx(i_x),i_x=1,i_x_max)
!       DO i_z = 1, i_z_max
!          WRITE(unit_field) zz(i_z)
!          WRITE(unit_field) (Indice_norm(i_x,i_z),i_x=1,i_x_max)
!       ENDDO
!       CLOSE(unit_field)
!       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
!    ENDIF

!920 FORMAT (F10.6)

!    RETURN
!  END SUBROUTINE field_out

END MODULE output
