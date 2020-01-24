module hdf5_stuff
USE HDF5
IMPLICIT NONE
CONTAINS

SUBROUTINE h5_add_units_1D(dset_id, units_value,error)
 
    IMPLICIT NONE
	INTEGER(HID_T) :: dset_id       ! Dataset identifier 
	character(*) :: units_value
	INTEGER :: error
	INTEGER(HSIZE_T), DIMENSION(1):: dumh51D
	! ATTRIBUTES OF HDF5 DATASETS
	 CHARACTER(LEN=5), PARAMETER :: aname = "units"   ! Attribute name
     INTEGER(HID_T) :: attr_id       ! Attribute identifier
     INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
     INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
     CHARACTER(LEN=10), DIMENSION(1) ::  attr_data  ! Attribute data

	
	! add attributes ( https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/fortran/examples/h5_crtatt.f90 )

    dumh51D = (/int(1,HSIZE_T)/) ! attribute dimension
    CALL h5screate_simple_f(1, dumh51D, aspace_id, error) ! Create scalar data space for the attribute. 1 stands for the rank
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) ! Create datatype for the attribute.
    CALL h5tset_size_f(atype_id, int(10,HSIZE_T), error) ! 10 is attribute length	
    CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error) ! Create dataset attribute.
    dumh51D = (/int(1,HSIZE_T)/) ! dimension of attributes
	attr_data(1) = units_value 
    CALL h5awrite_f(attr_id, atype_id, attr_data, dumh51D, error)
    CALL h5aclose_f(attr_id, error)  ! Close the attribute.
    CALL h5tclose_f(atype_id, error)  ! Close the attribute datatype.
    CALL h5sclose_f(aspace_id, error) ! Terminate access to the attributes data space.

	! attributes added

END SUBROUTINE  h5_add_units_1D

end module hdf5_stuff


!!!!!!!!!!!!!!!!!! THE MODULE ITSELF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! version created by Jan Vabek, 22/01/2020
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE outputHDF5
  USE fields
  USE parameters
  USE mpi_stuff
  USE long_step
  USE run_status
  USE normalization
  USE hdf5_stuff
CONTAINS

  SUBROUTINE HDF5_out
    USE fft
    USE HDF5
    IMPLICIT NONE

    ! General purpose variables: looping, dummy variables
    INTEGER(4) k1,k2
	REAL(4) dumr4
	INTEGER(HSIZE_T), DIMENSION(1):: dumh51D, dumh51D2
	
    ! HDF5 general purpose variables
     INTEGER(HID_T) :: file_id       ! File identifier 
     INTEGER(HID_T) :: dset_id       ! Dataset identifier 
     INTEGER(HID_T) :: dataspace     ! Dataspace identifier in file 
     INTEGER(HID_T) :: filespace     ! Filespace identifier
     INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
     INTEGER(HID_T) :: h5parameters  ! Property list identifier 
     INTEGER(HSIZE_T), DIMENSION(3) :: dims,dimsfi
     INTEGER(HSIZE_T), DIMENSION(3) :: ccount  
     INTEGER(HSIZE_T), DIMENSION(3) :: offset 
	 INTEGER :: error ! Error flags

	! HDF specific variables
     INTEGER(HSIZE_T)               :: r_offset

    ! code variables: the physical quantities etc.
     INTEGER :: field_dimensions ! Dataset rank & # of points in z
	 REAL(4), ALLOCATABLE :: Fields(:,:,:), rgrid(:), tgrid(:) ! the kind of this variable has to correspond with the precision stored in HDF5-file


    ! file & dataset names
     CHARACTER(LEN=10), PARAMETER :: filename = "results.h5"  ! File name
	 CHARACTER(LEN=17), PARAMETER :: Fields_dset_name = "IRprop/Fields_rzt" ! Dataset name
	 CHARACTER(LEN=12), PARAMETER :: zgrid_dset_name = "IRprop/zgrid" ! Dataset name
	 CHARACTER(LEN=12), PARAMETER :: tgrid_dset_name = "IRprop/tgrid" ! Dataset name
	 CHARACTER(LEN=12), PARAMETER :: rgrid_dset_name = "IRprop/rgrid" ! Dataset name

    ! THE CODE

    
	 IF (my_rank.EQ.0) THEN ! still in development mode, keep for the instant
	   print *, "HDF5 writting interation: ", HDF5write_count
	 ENDIF

    
	!!! ALLOCATING SPACE FOR GRIDS AND FIELDS
	field_dimensions = 3; ! total # of dims
	allocate(Fields(1,dim_r_local,dim_t)) ! space for fields in every itaration
	IF ( ( my_rank .EQ. 0 ) .AND. ( HDF5write_count .EQ. 1) ) THEN
	   allocate(tgrid(dim_t),rgrid(dim_r)) ! space for grids: first itration, proc # 0
	ENDIF
	

	! AFTER discussion with Jiri Vyskocil, he pointed out the column-majorness or row-majorness could a serious performance issue (note that c and FORTAN differ, h5 is a c-lib)
	r_offset = dim_r_start(num_proc)-1
	DO k1=1, dim_r_local	
	DO k2=1, dim_t
		Fields(1,k1,k2) = REAL( REAL( (efield_factor*efield_osc(k2)*e(k2,r_offset+k1)) ) , 4 ) ! SINGLE PRECISION, corresponding H5T_NATIVE_REAL (REAL(.,8) corresponds to H5T_NATIVE_DOUBLE)
	    ! e(t,r)
	ENDDO
	ENDDO

	! IF ( ( my_rank .EQ. 0 ) .AND. ( HDF5write_count .EQ. 1) ) THEN ! fill tables during the first call, proc # 0
	!    	DO k1=1, dim_t
	! 		tgrid(k1) = REAL( tps*(tlo+REAL(k1,8)*delta_t) , 4)
	!     ENDDO
	! 	DO k1=1, dim_r
	! 		rgrid(k1) = REAL( w0m*(REAL(k1-1,8)*delta_r) , 4)
	!     ENDDO
	! ENDIF



! This implementation almost straightforwadly follows tutorials from HDF5 portal. The only extension is our 3-dimensionality of the code.
! The code is in a "raw" form, without encapsulating in functions (see comment below).
! This is a stable version. However, it's probably not the final version. THere rest some serious things to test (row-, column-majorness etc.).
! We use pre-allocation of the file. The reason is that appending would require chunking, that is a rather advanced HDF5 topic, especially for performance issues.

!!!!!!!!!!
! The idea to extend this to writing it in multiple files is not to use ctrl-c--ctrl-v for new quantities. Almost all operations may be done in loops on various files, except hereogeneous writing. 
! When changing the code, please encapsulate the following in functions
!!!!!!!!!!

IF ( HDF5write_count == 1) THEN 

    !! THE WRITING OF THE FIELDS IS HERE

	CALL h5open_f(error) 
    
	CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) ! create HDF5 access parameters
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! set parameters for MPI access
	CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) ! Open collectivelly the file
	CALL h5pclose_f(h5parameters,error) ! close the parameters

  IF (my_rank.EQ.0) THEN ! only proc # 0 writes zgrid and other stuff in the first iteration

    !!!! HERE WE WRITE z-grid, appended in each iteration
	!!! extendible dataset for single-writter (following the tuto https://portal.hdfgroup.org/display/HDF5/Examples+from+Learning+the+Basics#ExamplesfromLearningtheBasics-changingex https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5/browse/fortran/examples/h5_extend.f90?at=89fbe00dec8187305b518d91c3ddb7d910665f79&raw )

	! CALL h5open_f(error) ! HDF5 initialise
	! CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file. 

    dumh51D = (/int(1,HSIZE_T)/) ! dim
	dumh51D2 = (/H5S_UNLIMITED_F/) ! maxdim
	CALL h5screate_simple_f(1, dumh51D, dataspace, error, dumh51D2 ) ! Create the data space with unlimited dimensions.
	CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5parameters, error) ! Modify dataset creation properties, i.e. enable chunking
	CALL h5pset_chunk_f(h5parameters, 1, dumh51D, error) ! enable chunking (1 is the dimension of the dataset)
	CALL h5dcreate_f(file_id, zgrid_dset_name, H5T_NATIVE_REAL, dataspace, dset_id, error, h5parameters) ! create dataset
	dumh51D = (/int(1,HSIZE_T)/) ! dimension of data
	dumr4 = REAL(four_z_Rayleigh*z,4) ! the actual z-coordinate in SI units 
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, dumr4, dumh51D, error)
    CALL h5_add_units_1D(dset_id, '[SI]',error)
	CALL h5sclose_f(dataspace, error)
    CALL h5pclose_f(h5parameters, error)
    CALL h5dclose_f(dset_id, error)
	! we will still be working with the file

	!!! we create rgrid and tgrid only in the first iteration
	dumh51D = (/int(dim_r,HSIZE_T)/) ! lengtsh of the grid
	CALL h5screate_simple_f(1, dumh51D, dataspace, error) ! Create the dataspace.
	CALL h5dcreate_f(file_id, rgrid_dset_name, H5T_NATIVE_REAL, dataspace, dset_id, error) ! create dataset
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, rgrid, dumh51D, error) ! write the data
	CALL h5_add_units_1D(dset_id, '[SI]',error)
	CALL h5sclose_f(dataspace, error)
	CALL h5dclose_f(dset_id, error)

    dumh51D = (/int(dim_t,HSIZE_T)/) ! lengtsh of the grid
	CALL h5screate_simple_f(1, dumh51D, dataspace, error) ! Create the dataspace.
	CALL h5dcreate_f(file_id, tgrid_dset_name, H5T_NATIVE_REAL, dataspace, dset_id, error) ! create dataset
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, tgrid, dumh51D, error) ! write the data
	CALL h5_add_units_1D(dset_id, '[SI]',error)
	CALL h5sclose_f(dataspace, error)
	CALL h5dclose_f(dset_id, error)

	! CALL h5fclose_f(file_id, error)
	! CALL h5close_f(error) ! Close FORTRAN interface.
	deallocate(rgrid,tgrid)

  ENDIF ! single-write end


    !!!!!
	!!!!! An extendible dataset seems to be a serious issue. We stick to pre-computing dataset size for now (see the piece ofthe code at the end of this file for details)
	!!!!! 
	dims = (/int(Nz_points,HSIZE_T),int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/) ! full prealocated dimension
	CALL h5screate_simple_f(field_dimensions, dims, filespace, error) ! Create the dataspace for the  dataset	
	CALL h5dcreate_f(file_id, Fields_dset_name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly

	offset = (/int(HDF5write_count-1,HSIZE_T),int(dim_r_start(num_proc)-1,HSIZE_T),int(0,HSIZE_T)/) ! (c-indexing from 0)
	ccount = (/int(1,HSIZE_T), int(dim_r_local,HSIZE_T) , int(dim_t,HSIZE_T)/) ! size of the chunk used by this MPI-worker

	CALL h5screate_simple_f(field_dimensions, ccount, memspace, error) ! dataset dimensions in memory (this worker)
	
	CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, ccount, error) ! hyperslab = part of the array acessed by this MPI-worker
	CALL h5pcreate_f(H5P_DATASET_XFER_F, h5parameters, error) ! Create access parametwers for writing
	CALL h5pset_dxpl_mpio_f(h5parameters, H5FD_MPIO_COLLECTIVE_F, error) ! specify the collective writting
	dimsfi = (/int(Nz_points,HSIZE_T),int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/) ! according to the tuto, it seems that whole dataset dimension is required
	CALL h5dwrite_f(dset_id , H5T_NATIVE_REAL, Fields, dimsfi, error,file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5parameters) ! Write the data collectivelly (we may try also to do it independently.... I think it could avoid some broadcast?)

    ! IF (my_rank.EQ.0) THEN !! testing
    !    OPEN(15,FILE='test.DAT',STATUS='UNKNOWN')
	!    DO k1 = 1,dim_r_local
	!       WRITE(15,*) FIELDS(1,k1,1024)
	!    ENDDO
    !    CLOSE(15)
    ! ENDIF

	!close the files etc.
	CALL h5sclose_f(filespace,error)
	CALL h5sclose_f(memspace,error)
	CALL h5dclose_f(dset_id,error)
	CALL h5pclose_f(h5parameters,error)
	CALL h5fclose_f(file_id,error)	
	CALL h5close_f(error) ! close the HDF5 workspace

    !!! attributes are probably not well handled by MPI... ( https://forum.hdfgroup.org/t/write-attributes-collectively-in-mpi-run/4902/2 ), just add them once by one worker
    IF (my_rank.EQ.0) THEN
    CALL h5open_f(error) 
	CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file.
	CALL h5dopen_f(file_id, Fields_dset_name, dset_id, error)   !Open the  dataset
	CALL h5_add_units_1D(dset_id, '[V/m]',error) ! add units
	CALL h5dclose_f(dset_id,error)
	CALL h5fclose_f(file_id,error)
	CALL h5close_f(error) ! close the HDF5 workspace
	ENDIF

ELSE !!!! APPENDING THE DATA IN NEXT ITERATIONS

!    IF ( (my_rank.EQ.0) .AND. (HDF5write_count .EQ. 2)) THEN !! testing
!        OPEN(15,FILE='test2.DAT',STATUS='UNKNOWN')
! 	   DO k1 = 1,dim_r_local
! 	      WRITE(15,*) FIELDS(1,k1,1024)
! 	   ENDDO
!        CLOSE(15)
!     ENDIF



    
	CALL h5open_f(error)  !Initialize HDF5
	CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) !define parameters of HDF5 workflow for MPI-access
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! allow MPI access
	CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) !Open collectivelly the file
	CALL h5pclose_f(h5parameters,error) ! parameters were used for MPI open, close them

	IF (my_rank.EQ.0) THEN ! only one worker is extending the zgrid

      ! only z-grid in 1D
	  CALL h5open_f(error)
	  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file.
	  CALL h5dopen_f(file_id, zgrid_dset_name, dset_id, error)   !Open the  dataset

	  dumh51D = (/int(HDF5write_count,HSIZE_T)/) ! new dimension of the dataset
	  CALL h5dset_extent_f(dset_id, dumh51D, error) ! extend the dataset
	  CALL h5dget_space_f(dset_id, dataspace, error) ! get the dataspace of the dataset

	  dumh51D = (/int(1,HSIZE_T)/) ! dimension of the memspace (it's the chunk that is appended to the dataset)
	  CALL h5screate_simple_f (1, dumh51D, memspace, error) ! create memory space

      dumh51D = (/int(HDF5write_count-1,HSIZE_T)/) ! offset
	  dumh51D2 = (/int(1,HSIZE_T)/) ! count
	  CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, dumh51D, dumh51D2, error) ! choose the hyperslab in the file
	  dumh51D = (/int(1,HSIZE_T)/) ! the dimension of data written
	  dumr4 = REAL(four_z_Rayleigh*z,4) ! the actual z-coordinate in SI units 
	  CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, dumr4, dumh51D, error, memspace, dataspace) ! wrtiting the data

	  CALL h5sclose_f(memspace, error)
	  CALL h5sclose_f(dataspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5fclose_f(file_id, error)
	  CALL h5close_f(error) ! Close FORTRAN interface.

    ENDIF ! single-write end
	
	CALL h5dopen_f(file_id, Fields_dset_name, dset_id, error) ! open the dataset (already created)
    CALL h5dget_space_f(dset_id,filespace,error) ! filespace from the dataset (get instead of create)

	offset = (/int(HDF5write_count-1,HSIZE_T),int(dim_r_start(num_proc)-1,HSIZE_T),int(0,HSIZE_T)/) ! (c-indexing from 0)
	ccount = (/int(1,HSIZE_T), int(dim_r_local,HSIZE_T) , int(dim_t,HSIZE_T)/) ! size of the chunk used by this MPI-worker
    CALL h5screate_simple_f(field_dimensions, ccount, memspace, error) ! dataset dimensions in memory (this worker)

	CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, ccount, error) ! hyperslab = part of the array acessed by this MPI-worker	
	CALL h5pcreate_f(H5P_DATASET_XFER_F, h5parameters, error) ! Create access parametwers
	CALL h5pset_dxpl_mpio_f(h5parameters, H5FD_MPIO_COLLECTIVE_F, error) ! collective writting
	dimsfi = (/int(Nz_points,HSIZE_T),int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/) ! according to the tuto, it seems that whole dataset dimension is required
	CALL h5dwrite_f(dset_id , H5T_NATIVE_REAL, Fields, dimsfi, error,file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5parameters) ! write data

    ! closing (all hdf5 objects)
	CALL h5dclose_f(dset_id,error)
	CALL h5sclose_f(filespace,error)
	CALL h5sclose_f(memspace,error)
	CALL h5pclose_f(h5parameters,error)
	CALL h5fclose_f(file_id,error)	
	CALL h5close_f(error) ! close the HDF5 workspace

ENDIF


	HDF5write_count = HDF5write_count + 1 !increase counter in all cases
	deallocate(Fields)

	IF (my_rank.EQ.0) THEN
	  print *, "finished", my_rank, "iteration+1", HDF5write_count
	ENDIF

    RETURN
  END SUBROUTINE  HDF5_out

END MODULE outputHDF5




!!!!! 2d attributes version
    ! dumh51D = (/int(2,HSIZE_T)/) ! attribute dimension
    ! CALL h5screate_simple_f(1, dumh51D, aspace_id, error) ! Create scalar data space for the attribute. 1 stands for the rank
    ! CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) ! Create datatype for the attribute.
    ! CALL h5tset_size_f(atype_id, int(80,HSIZE_T), error) ! 80 is attribute length	
    ! CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error) ! Create dataset attribute.
    ! dumh51D = (/int(2,HSIZE_T)/) ! for specification and its values
	! attr_data(1) = "units"
    ! attr_data(2) = "[SI]"    
    ! CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims, error)
    ! CALL h5aclose_f(attr_id, error)  ! Close the attribute.
    ! CALL h5tclose_f(atype_id, error)  ! Close the attribute datatype.
    ! CALL h5sclose_f(aspace_id, error) ! Terminate access to the attributes data space.








	    ! This is the research I did about the implentation of extendible datesets:
	! maxdims = (/int(3,HSIZE_T), int(128,HSIZE_T), int(2,HSIZE_T)/) ! maxdims = (/H5S_UNLIMITED_F, int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/) 
	! CALL h5screate_simple_f(field_dimensions, dims, filespace, error, maxdims) ! Create the data space for the  dataset. 
	! ! the caused error are discussed in the following forums
	! https://github.com/h5py/h5py/issues/115
	! https://www.mathworks.com/matlabcentral/answers/161779-why-do-i-get-an-error-from-the-hdf5-library-when-i-try-to-create-an-hdf5-dataset-with-one-of-the-dim
	! http://hdf-forum.184993.n3.nabble.com/Errors-when-creating-dataset-when-dataspace-has-a-max-dimension-of-H5S-UNLIMITED-td2900828.html
	! According to a discussion with Jiri Vyskocil, the solution can be also in the chunking
	! (see also Jiri's recommends': 
	! https://forum.hdfgroup.org/t/writing-to-an-extendable-dataset-in-a-loop-c-c/5793?fbclid=IwAR3M_kX4L6ZqC1teRWGQuslaToFCIYi7qPUCb8uA5Gerh4rrlfgJvco5dFk
	! https://forum.hdfgroup.org/t/parallel-dataset-resizing-strategies/3826?fbclid=IwAR1DlG6xXlYLp7XmAET3IrP4gbcdDIzaTTKI31ZsDW9h_AHScEDmLrdOFek
	! )

	!!! EXTENDIBLE DATASETS ARE PROBABLY A PROBLEM THIS is an error message from this version
! 	HDF5-DIAG: Error detected in HDF5 (1.10.5) MPI-process 21:
!   #000: H5D.c line 145 in H5Dcreate2(): unable to create dataset
!     major: Dataset
!     minor: Unable to initialize object
!   #001: H5Dint.c line 329 in H5D__create_named(): unable to create and link to dataset
!     major: Dataset
!     minor: Unable to initialize object
!   #002: H5L.c line 1557 in H5L_link_object(): unable to create new link to object
!     major: Links
!     minor: Unable to initialize object
!   #003: H5L.c line 1798 in H5L__create_real(): can't insert link
!     major: Links
!     minor: Unable to insert object
!   #004: H5Gtraverse.c line 851 in H5G_traverse(): internal path traversal failed
!     major: Symbol table
!     minor: Object not found
!   #005: H5Gtraverse.c line 627 in H5G__traverse_real(): traversal operator failed
!     major: Symbol table
!     minor: Callback failed
!   #006: H5L.c line 1604 in H5L__link_cb(): unable to create object
!     major: Links
!     minor: Unable to initialize object
!   #007: H5Oint.c line 2453 in H5O_obj_create(): unable to open object
!     major: Object header
!     minor: Can't open object
!   #008: H5Doh.c line 300 in H5O__dset_create(): unable to create dataset
!     major: Dataset
!     minor: Unable to initialize object
!   #009: H5Dint.c line 1274 in H5D__create(): unable to construct layout information
!     major: Dataset
!     minor: Unable to initialize object
!   #010: H5Dcontig.c line 402 in H5D__contig_construct(): extendible contiguous non-external dataset not allowed
!     major: Dataset
!     minor: Feature is unsupported
