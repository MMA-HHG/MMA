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
  USE hdf5_helper
CONTAINS

  SUBROUTINE HDF5_out
  ! This implementation almost straightforwadly follows tutorials from HDF5 portal. The only extension is our 3-dimensionality of the code.
  ! The code is in a "raw" form, without encapsulating in functions (see comment below).
  ! This is a stable version. However, it's probably not the final version. THere rest some serious things to test (row-, column-majorness etc.).
  ! We use pre-allocation of the file. The reason is that appending would require chunking, that is a rather advanced HDF5 topic, especially for performance issues.
  ! There is need to call some operations collectivelly https://portal.hdfgroup.org/display/HDF5/Collective+Calling+Requirements+in+Parallel+HDF5+Applications
  !   - It seems that barriers are induced by these commands.
  !   - We encountered a problem when we tried to do single-worker operation before the collective open. Working with one file has to be done with cre if the file is accessed both sequentially and parallelly. We encountered some problems with that, doing collective first and sequential after seems stable.

    USE fft
    USE HDF5
    IMPLICIT NONE

    ! General purpose variables: looping, dummy variables
    INTEGER(4) k1,k2
    REAL(4) dumr4
    INTEGER(HSIZE_T), DIMENSION(1):: dumh51D, dumh51D2

    ! HDF5 general purpose variables
    INTEGER(HID_T) :: file_id       ! File identifier 
    INTEGER(HID_T) :: group_id      ! Group identifier 
    INTEGER(HID_T) :: h5parameters  ! Property list identifier 
    INTEGER(HSIZE_T), DIMENSION(3) :: dims
    INTEGER(HSIZE_T), DIMENSION(3) :: ccount  
    INTEGER(HSIZE_T), DIMENSION(3) :: offset 
    INTEGER :: error ! Error flags

    ! HDF specific variables
    INTEGER(HSIZE_T)               :: r_offset

    ! code variables: the physical quantities etc.
    INTEGER :: field_dimensions ! Dataset rank & # of points in z
    ! the kind of this variable has to correspond with the precision stored in HDF5-file
    REAL(4), ALLOCATABLE :: fields_array(:,:,:), rgrid(:), tgrid(:) 

    ! file & dataset names
     CHARACTER(LEN=10), PARAMETER :: filename = "results.h5"  ! File name
     CHARACTER(LEN=17), PARAMETER :: Fields_dset_name = "IRprop/Fields_rzt" ! Dataset name
     CHARACTER(LEN=12), PARAMETER :: zgrid_dset_name = "IRprop/zgrid" ! Dataset name
     CHARACTER(LEN=12), PARAMETER :: tgrid_dset_name = "IRprop/tgrid" ! Dataset name
     CHARACTER(LEN=12), PARAMETER :: rgrid_dset_name = "IRprop/rgrid" ! Dataset name
     CHARACTER(LEN=12), PARAMETER :: groupname = "IRprop"
    ! THE CODE

    
    IF (my_rank.EQ.0) THEN ! still in development mode, keep for the instant
      print *, "HDF5 writting iteration: ", HDF5write_count
    ENDIF

    !!! ALLOCATING SPACE FOR GRIDS AND FIELDS
    field_dimensions = 3; ! total # of dims
    allocate(fields_array(1,dim_r_local,dim_t)) ! space for fields in every itaration
    IF ( ( my_rank .EQ. 0 ) .AND. ( HDF5write_count .EQ. 1) ) THEN
      allocate(tgrid(dim_t),rgrid(dim_r)) ! space for grids: first itration, proc # 0
    ENDIF


    ! AFTER discussion with Jiri Vyskocil, he pointed out the column-majorness or row-majorness could a serious performance issue (note that c and FORTAN differ, h5 is a c-lib)
    r_offset = dim_r_start(num_proc)-1
    DO k1=1, dim_r_local
      DO k2=1, dim_t
        fields_array(1,k1,k2) = REAL( REAL( (efield_factor*efield_osc(k2)*e(k2,r_offset+k1)) ) , 4 ) ! SINGLE PRECISION, corresponding H5T_NATIVE_REAL (REAL(.,8) corresponds to H5T_NATIVE_DOUBLE)
        ! e(t,r)
      ENDDO
    ENDDO

    IF ( ( my_rank .EQ. 0 ) .AND. ( HDF5write_count .EQ. 1) ) THEN ! fill tables during the first call, proc # 0
      DO k1=1, dim_t
        tgrid(k1) = REAL( tps*(tlo+REAL(k1,8)*delta_t) , 4)
      ENDDO
      DO k1=1, dim_r
        rgrid(k1) = REAL( w0m*(REAL(k1-1,8)*delta_r) , 4)
      ENDDO
    ENDIF



    !!!!!!!!!!
    ! The idea to extend this to writing it in multiple files is not to use ctrl-c--ctrl-v for new quantities. Almost all operations may be done in loops on various files, except hereogeneous writing. 
    ! When changing the code, please encapsulate the following in functions
    !!!!!!!!!!
 
   ! Initialize dims arrays
    dims = (/int(Nz_points,HSIZE_T),int(dim_r,HSIZE_T), int(dim_t,HSIZE_T)/)
    offset = (/int(HDF5write_count-1,HSIZE_T),int(dim_r_start(num_proc)-1,HSIZE_T),int(0,HSIZE_T)/)
    ccount = (/int(1,HSIZE_T), int(dim_r_local,HSIZE_T) , int(dim_t,HSIZE_T)/)
      
    CALL h5open_f(error) 
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, h5parameters, error) ! create HDF5 access parameters
    CALL h5pset_fapl_mpio_f(h5parameters, MPI_COMM_WORLD, MPI_INFO_NULL, error) ! set parameters for MPI access
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = h5parameters ) ! Open collectivelly the file
    CALL h5pclose_f(h5parameters,error) ! close the parameters

    IF ( HDF5write_count == 1) THEN 
      !Create group for the output
      CALL h5gcreate_f(file_id, groupname, group_id, error) 
      CALL h5gclose_f(group_id, error)
          
      ! Call writing routine
      CALL create_3D_array_real_dset_p(file_id, Fields_dset_name, fields_array, dims, offset, ccount)

      ! Terminate
      CALL h5fclose_f(file_id,error)

      !!! attributes are probably not well handled by MPI... ( https://forum.hdfgroup.org/t/write-attributes-collectively-in-mpi-run/4902/2 ), just add them once by one worker
      IF (my_rank.EQ.0) THEN
        CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Open an existing file.
        CALL h5_add_units_1D(file_id, Fields_dset_name, '[V/m]') ! add units

        !!!! HERE WE WRITE z-grid, appended in each iteration
        ! we will still be working with the file
        dumr4 = REAL(four_z_Rayleigh*z,4) ! the actual z-coordinate in SI units 
        CALL create_1D_dset_unlimited(file_id, zgrid_dset_name, dumr4)
        CALL h5_add_units_1D(file_id, zgrid_dset_name, '[SI]')

        CALL create_dset(file_id, rgrid_dset_name, rgrid, dim_r)
        CALL h5_add_units_1D(file_id, rgrid_dset_name, '[SI]')
        CALL create_dset(file_id, tgrid_dset_name, tgrid, dim_t)
        CALL h5_add_units_1D(file_id, tgrid_dset_name, '[SI]')

        CALL h5fclose_f(file_id, error)
        deallocate(rgrid,tgrid)

      ENDIF ! single-write end
    ELSE !!!! APPENDING THE DATA IN NEXT ITERATIONS
      CALL write_hyperslab_to_dset(file_id, Fields_dset_name, fields_array, offset, ccount)
      CALL h5fclose_f(file_id,error)
      IF (my_rank.EQ.0) THEN ! only one worker is extending the zgrid
        CALL h5open_f(error)  !Initialize HDF5
        CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
        ! only z-grid in 1D
        dumh51D = (/int(HDF5write_count-1,HSIZE_T)/) ! offset
        dumh51D2 = (/int(1,HSIZE_T)/) ! count
        dumr4 = REAL(four_z_Rayleigh*z,4) ! the actual z-coordinate in SI units 
        CALL extend_1D_dset_unlimited(file_id, zgrid_dset_name, dumr4, new_dims=(/int(HDF5write_count,HSIZE_T)/), & 
          memspace_dims=(/int(1,HSIZE_T)/), offset=dumh51D, hyperslab_size=dumh51D2)
        CALL h5fclose_f(file_id,error)
      ENDIF ! single-write end
    ENDIF
    CALL h5close_f(error) ! close the HDF5 workspace


    HDF5write_count = HDF5write_count + 1 !increase counter in all cases
    deallocate(fields_array)

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
