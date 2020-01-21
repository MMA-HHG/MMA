MODULE outputHDF5
  USE fields
  USE parameters
  USE mpi_stuff
  USE long_step
  USE run_status
  USE normalization
CONTAINS

  SUBROUTINE HDF5_out
    USE fft
    USE HDF5
    IMPLICIT NONE

    INTEGER(4) j,k,l
    INTEGER(4) k1,k2,k3,k4,k5
    REAL(8) rhotemp,r,mpa
    COMPLEX(8) help
    CHARACTER*10 iz !,filename



     CHARACTER(LEN=10), PARAMETER :: filename = "results.h5"  ! File name
     CHARACTER(LEN=23), PARAMETER :: dsetname = "micro/FieldsForTDSEreal" ! Dataset name

     INTEGER(HID_T) :: file_id       ! File identifier 
     INTEGER(HID_T) :: dset_id       ! Dataset identifier 
     INTEGER(HID_T) :: dataspace     ! Dataspace identifier in file 
     INTEGER(HID_T) :: filespace     ! Filespace identifier
     INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
     INTEGER(HID_T) :: h5parameters  ! Property list identifier 

!     INTEGER(HSIZE_T), DIMENSION(3) :: dimsf  ! Dataset dimensions.
!!     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
!     INTEGER(HSIZE_T), DIMENSION(3) :: dimsfi 
     INTEGER(HSIZE_T), DIMENSION(3) :: dims,dimsfi,maxdims 
	 INTEGER(HSIZE_T), DIMENSION(1):: dumh51D, dumh51D2

     INTEGER(HSIZE_T), DIMENSION(3) :: ccount  
     INTEGER(HSIZE_T), DIMENSION(3) :: offset 
     INTEGER(HSIZE_T), DIMENSION(3) :: stride
     INTEGER(HSIZE_T), DIMENSION(3) :: cblock
     INTEGER(HSIZE_T)               :: r_offset
     INTEGER :: field_dimensions ! Dataset rank & # of points in z

     INTEGER :: error, error_n  ! Error flags
     !
     ! MPI definitions and calls.
     !

    !  INTEGER :: comm, info


     ! code variables
     REAL(4), ALLOCATABLE :: Fields(:,:,:) ! the kind of this variable has to correspond with the precision stored in HDF5-file
     INTEGER :: Nz_dim_old

	 ! testing variables
	 INTEGER, DIMENSION(4,6) :: dset_data, data_out ! Data buffers
     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
	 CHARACTER(LEN=7), PARAMETER :: filename2 = "test.h5" ! Dataset name
	 CHARACTER(LEN=16), PARAMETER :: dsetname2 = "TestCUPRADSingle" ! Dataset name
	 CHARACTER(LEN=18), PARAMETER :: dsetname3 = "TestCUPRADParallel" ! Dataset name
	 CHARACTER(LEN=5), PARAMETER :: dsetname4 = "zgrid" ! Dataset name


     !  comm = MPI_COMM_WORLD
     !  info = MPI_INFO_NULL


	!  print *, "HDF5 output accessed, proc", my_rank

	 IF (my_rank.EQ.0) THEN
	   print *, "writting interation: ", HDF5write_count
	 ENDIF


    !!! in the first run, create dataset and fill random data
	field_dimensions = 3;
	allocate(Fields(1,1,2)) !allocate(Fields(1,dim_r_end(num_proc)-dim_r_start(num_proc),dim_t))

	! AFTER discussion with Jiri Vyskocil, he pointed out the column-majorness or row-majorness could a serious performance issue (note that c and FORTAN differs, h5 is a c-lib)
	r_offset = dim_r_start(num_proc)-1
	DO k1=1, 1 !( dim_r_end(num_proc)-dim_r_start(num_proc) )	
	DO k2=1, 2 !dim_t
		! Fields(1,k1,k2) = REAL(my_rank+HDF5write_count+k1+k2,8) !REAL(e(k2,r_offset+k1));
		Fields(1,k1,k2) = REAL(my_rank,4) !REAL(e(k2,r_offset+k1));  ! SINGEL PRECISION, corresponding H5T_NATIVE_REAL
		! Fields(1,k1,k2) = REAL(k2,8) !REAL(e(k2,r_offset+k1)); ! SINGEL PRECISION, corresponding H5T_NATIVE_DOUBLE
	ENDDO
	ENDDO

	! print *, "fields allocated, proc", my_rank





! At the end, this implementation almost straightforwadly follows tutorial from HDF5 page. The only extension is our 3-dimensionality of the code
! The idea to extend this to writing it in multiple files is not to use ctrl-c--ctrl-v for new quantities. Almost all operations may be done in loops on various files, except hereogeneous writing
	IF ( HDF5write_count == 1) THEN


  !!!!!!!!!!!! HDF5 testing 
  ! I test HDF dunctionality just in single-writer opearation
  ! in one relase we wrote data here, we keep it now as a dummy write
  IF (my_rank.EQ.0) THEN ! only one worker

    print *, "HDF5 testfile IF accessed"
	! Initialize the dset_data array.
	DO k1 = 1, 4
		DO k2 = 1, 6
			dset_data(k1,k2) = (k1-1)*6 + k2
		END DO
	END DO

	data_dims(1) = 4
	data_dims(2) = 6

	!
	! Initialize FORTRAN interface.
	CALL h5open_f(error)
	! CALL h5fcreate_f(filename2, H5F_ACC_TRUNC_F, file_id, error) ! create test file 
	CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file. 
    CALL h5screate_simple_f(2, data_dims, dataspace, error) ! Create the dataspace.
	CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_INTEGER, dataspace, dset_id, error) ! create the dataset
	! CALL h5dopen_f(file_id, dsetname2, dset_id, error)  ! Open an existing dataset.

	CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error) ! Write the dataset.
	CALL h5dclose_f(dset_id, error) ! Close the dataset.
	CALL h5fclose_f(file_id, error) ! Close the file.
	CALL h5close_f(error) ! Close FORTRAN interface.
	!!!!!!!!!!!! HDF5 testing 

	!!! extendible dataset for single-writter (following the tuto https://portal.hdfgroup.org/display/HDF5/Examples+from+Learning+the+Basics#ExamplesfromLearningtheBasics-changingex https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5/browse/fortran/examples/h5_extend.f90?at=89fbe00dec8187305b518d91c3ddb7d910665f79&raw )
	CALL h5open_f(error)
	CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file. 

    dumh51D = (/int(1,HSIZE_T)/)
	dumh51D2 = (/H5S_UNLIMITED_F/) !Create the data space with unlimited dimensions.
	CALL h5screate_simple_f(1, dumh51D, dataspace, error, dumh51D2 )
	CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5parameters, error)   !Modify dataset creation properties, i.e. enable chunking
	CALL h5pset_chunk_f(h5parameters, 1, dumh51D, error)
	CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_REAL, dataspace, dset_id, error, h5parameters)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, REAL(HDF5write_count,4), data_dims, error)
	CALL h5sclose_f(dataspace, error)
    CALL h5pclose_f(h5parameters, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5fclose_f(file_id, error)
	CALL h5close_f(error) ! Close FORTRAN interface.


	! !! test appending at the instant
	! CALL h5open_f(error)
	! CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file. 

	! ! maxdims = (/H5S_UNLIMITED_F/) !Create the data space with unlimited dimensions.
	! ! CALL h5screate_simple_f(RANK, dims, dataspace, error, maxdims)
	! ! CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5parameters, error)   !Modify dataset creation properties, i.e. enable chunking
	! ! dumh51D = (/int(1,HSIZE_T)/)
	! ! CALL h5pset_chunk_f(h5parameters, 1, dumh51D, error)
	! ! CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_REAL, dataspace, dset_id, error, h5parameters)
	! CALL h5dopen_f(file_id, dsetname4, dset_id, error)   !Open the  dataset
	! CALL h5dset_extent_f(dset_id, (/int(HDF5write_count+1,HSIZE_T)/), error)
	! dumh51D = (/int(1,HSIZE_T)/)
	! CALL h5screate_simple_f (1, dumh51D, memspace, error)

	! CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, (/int(HDF5write_count,HSIZE_T)/), (/int(1,HSIZE_T)/), error)
	! CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, REAL(HDF5write_count,4), data_dims, error, memspace, dataspace)
	! CALL h5sclose_f(dataspace, error)
    ! CALL h5dclose_f(dset_id, error)
    ! CALL h5fclose_f(file_id, error)
	! CALL h5close_f(error) ! Close FORTRAN interface.


  ENDIF ! single-write end



    ! CALL MPI_Barrier(MPI_COMM_WORLD,ierr) !! try barrier here
	!Initialize HDF5

	! print *, "before h5 init, proc", my_rank
	CALL h5open_f(error) 
    
	!define parameters of HDF5 workflow for MPI-access
	! print *, "before h5 param create, proc", my_rank  
	CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) ! create access parameters
	! print *, "before h5 param set, proc", my_rank
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! allow MPI access (should it be here?)

	!Open collectivelly the file
	! print *, "before h5 filecreation, proc", my_rank
	! CALL h5fcreate_f(filename2, H5F_ACC_TRUNC_F, file_id, error, access_prp = h5parameters) ! we again first test creating the file collectivelly
	! CALL h5fcreate_f(filename2, H5F_ACC_TRUNC_F, file_id, error)  ! single - version
	CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) ! open file collectivelly

!CINES correction:	CALL h5pclose(h5parameters) ! parameters were used for MPI open, close them
!CINES h5pclose_f is the correct interface for fortran90
	CALL h5pclose_f(h5parameters,error) ! parameters were used for MPI open, close them
	



	! Extendible dataset seems to be a serious issue. We stick to pre-computing dataset size for now	
	dims = (/int(Nz_points,HSIZE_T),int(128,HSIZE_T), int(2,HSIZE_T)/) !dims = (/int(1,HSIZE_T),int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/) ! only line per proc. now, code runned on 128
	! print *, "before h5 dataspace creation, proc", my_rank	
	CALL h5screate_simple_f(field_dimensions, dims, filespace, error) ! Create the data space for the  dataset. ! maybe problem with the exension??? !!!!!!
	

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


    !! Code continuation


	! Maybe we don't need do this
	! dims = (/1, 1, dim_t/) !dims = (/1, dim_r_end(num_proc)-dim_r_start(num_proc), dim_t/) ! dimension of my field
	! CALL h5screate_simple_f(field_dimensions, dims, dataspace, error, maxdims) ! dataset dimensions in memory (this worker)

    ! print *, "before h5 dataset creation, proc", my_rank
	! we create the dataset collectivelly
	CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_REAL, filespace, dset_id, error)
!CINES correction	CALL h5sclose(filespace,error)
	! CALL h5sclose_f(filespace,error)
	! CALL h5dclose_f(dset_id,error) !! remove this now

  

	!we use hyperslab to assign part of the global dataset
	!chunk data for each worker

	! the following two lines should be inspected for performance
	! stride = (/1,1,1/) ! we write a block of data directly in file, no skipped lines, rows, ...
	! cblock = (/1,1,2/) ! cblock = (/1,dim_r_end(num_proc) - dim_r_start(num_proc),dim_t/) ! allows flush data at once 

	offset = (/int(HDF5write_count-1,HSIZE_T),int(my_rank,HSIZE_T),int(0,HSIZE_T)/) ! offset = (/0,dim_r_start(num_proc),0/) ! c-indexing from 0
	ccount = (/int(1,HSIZE_T), int(1,HSIZE_T) , int(2,HSIZE_T)/) ! ccount = (/1, dim_r_end(num_proc) - dim_r_start(num_proc) , dim_t/)

	! memory space allocated for each worker is here
	CALL h5screate_simple_f(field_dimensions, ccount, memspace, error) ! dataset dimensions in memory (this worker)
	
	! CALL h5dget_space_f(dset_id,filespace,error)
	! print *, "before h5 hyperslab, proc", my_rank
	CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, ccount, error) ! we should have access to its part of the dataset for each worker
	!CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, ccount, error, stride, cblock) ! these extra parameters should be possible in a genralised case


	!!!Finally, write data
	! Create access parametwers
	CALL h5pcreate_f(H5P_DATASET_XFER_F, h5parameters, error)
	! print *, "before h5 dxpl mpio, proc", my_rank
	CALL h5pset_dxpl_mpio_f(h5parameters, H5FD_MPIO_COLLECTIVE_F, error) ! collective writting

!	CALL h5pset_chunk_f(crp_list, field_dimensions, dimsc, error) ???????????? Do we need chunk it?

	! Write the data collectivelly (we may try also to do it independently.... I think it could avoid some broadcast?)
	dimsfi = (/Nz_points,128,2/) ! dimsfi = (/1,dim_r,dim_t/) ! according to the tuto, it seems that whole dataset dimension is required
	! print *, "before h5 parallel write, proc", my_rank
	CALL h5dwrite_f(dset_id , H5T_NATIVE_REAL, Fields, dimsfi, error,file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5parameters)
	! CALL h5dwrite_f(dset_id , H5T_NATIVE_REAL, Fields, dimsfi, error,file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5parameters)! general params, data are written !!!( probably variable length)


	!close the files etc.
	CALL h5sclose_f(filespace,error)
	CALL h5sclose_f(memspace,error)
	CALL h5dclose_f(dset_id,error)
	CALL h5pclose_f(h5parameters,error)
	CALL h5fclose_f(file_id,error)	
	CALL h5close_f(error) ! close the HDF5 workspace

	! print *, "after h5 fileclose, proc", my_rank
   
!     !!! now, just append data 


! 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	ELSE ! We're just appending the data now

  IF (my_rank.EQ.0) THEN ! only one worker

    print *, "HDF5 z-extension accessed"

	!!! extendible dataset for single-writter 
	CALL h5open_f(error)
	CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file. 

	! maxdims = (/H5S_UNLIMITED_F/) !Create the data space with unlimited dimensions.
	! CALL h5screate_simple_f(RANK, dims, dataspace, error, maxdims)
	! CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5parameters, error)   !Modify dataset creation properties, i.e. enable chunking
	! dumh51D = (/int(1,HSIZE_T)/)
	! CALL h5pset_chunk_f(h5parameters, 1, dumh51D, error)
	! CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_REAL, dataspace, dset_id, error, h5parameters)
	CALL h5dopen_f(file_id, dsetname4, dset_id, error)   !Open the  dataset
	CALL h5dset_extent_f(dset_id, (/int(HDF5write_count,HSIZE_T)/), error)
	dumh51D = (/int(1,HSIZE_T)/)
	CALL h5screate_simple_f (1, dumh51D, memspace, error)
	CALL h5dget_space_f(dset_id, dataspace, error)
	CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, (/int(HDF5write_count-1,HSIZE_T)/), (/int(1,HSIZE_T)/), error)
	CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, REAL(HDF5write_count,4), data_dims, error, memspace, dataspace)
	CALL h5sclose_f(memspace, error)
	CALL h5sclose_f(dataspace, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5fclose_f(file_id, error)
	CALL h5close_f(error) ! Close FORTRAN interface.


	! CALL h5dget_space_f(dset_id, dataspace, error)  !Get dataset's dataspace handle.

	! !CALL h5sget_simple_extent_ndims_f(dataspace, rankr, error)   !Get dataspace's rank. ! that's in the tuto, we already know it's 1
    
    ! CALL h5sget_simple_extent_dims_f(dataspace, dumh51D, dumh51D2, error)   !Get dataspace's dimensions.


  ENDIF ! single-write end

    
	CALL h5open_f(error)  !Initialize HDF5
	CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) !define parameters of HDF5 workflow for MPI-access
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! allow MPI access
	CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) !Open collectivelly the file
	CALL h5pclose_f(h5parameters,error) ! parameters were used for MPI open, close them
	

    ! The code should differ now, since dataset is already created, we should be able to get filespace from it without creating it
	! first, we open the dataset
    ! print *, "before h5 dataset open, proc", my_rank, "iteration", HDF5write_count
	CALL h5dopen_f(file_id, dsetname3, dset_id, error) ! this should be enough !h5dcreate_f(file_id, dsetname3, H5T_NATIVE_REAL, filespace, dset_id, error)

    CALL h5dget_space_f(dset_id,filespace,error) ! filespace shoulb be obtained now from the file ! CALL h5screate_simple_f(field_dimensions, dims, filespace, error) ! Create the data space for the  dataset. ! maybe problem with the exension??? !!!!!!

    ! the only change is the offset, but linked with the cummulative HDF5write_count
	offset = (/int(HDF5write_count-1,HSIZE_T),int(my_rank,HSIZE_T),int(0,HSIZE_T)/) ! offset = (/0,dim_r_start(num_proc),0/) ! c-indexing from 0
	ccount = (/int(1,HSIZE_T), int(1,HSIZE_T) , int(2,HSIZE_T)/) ! ccount = (/1, dim_r_end(num_proc) - dim_r_start(num_proc) , dim_t/)

	! we select the hyperslab
	CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, ccount, error) ! we should have access to its part of the dataset for each worker

    ! print *, "before h5 writing, proc", my_rank, "iteration", HDF5write_count

	!!!Finally, write data	
	CALL h5pcreate_f(H5P_DATASET_XFER_F, h5parameters, error) ! Create access parametwers
	CALL h5pset_dxpl_mpio_f(h5parameters, H5FD_MPIO_COLLECTIVE_F, error) ! collective writting
	dimsfi = (/Nz_points,128,2/) ! dimsfi = (/1,dim_r,dim_t/) ! according to the tuto, it seems that whole dataset dimension is required
	CALL h5dwrite_f(dset_id , H5T_NATIVE_REAL, Fields, dimsfi, error,file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5parameters)

    ! closing
	CALL h5dclose_f(dset_id,error)
	CALL h5sclose_f(filespace,error)
	CALL h5sclose_f(memspace,error)
	CALL h5dclose_f(dset_id,error)
	CALL h5pclose_f(h5parameters,error)
	CALL h5fclose_f(file_id,error)	
	CALL h5close_f(error) ! close the HDF5 workspace

	ENDIF
	HDF5write_count = HDF5write_count + 1 !increase counter in all cases
	deallocate(Fields)

	IF (my_rank.EQ.0) THEN
	  print *, "finished", my_rank, "iteration+1", HDF5write_count
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

END MODULE outputHDF5
