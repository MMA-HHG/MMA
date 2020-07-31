MODULE hdf5_helper
  USE HDF5
  INTERFACE read_dset
    PROCEDURE readint, readreal, readbool, readstring, read_array_complex_dset, read_2D_array_complex_dset, &
      read_2D_array_complex_dset_slice, read_2D_array_real_dset_slice
  END INTERFACE
  INTERFACE create_dset
    PROCEDURE create_scalar_boolean_dset, create_scalar_int_dset, create_scalar_real_dset, create_1D_array_real_dset, &
      create_array_complex_dset, create_2D_array_complex_dset, create_2D_array_real_dset
  END INTERFACE
  CONTAINS

    !******!
    ! READ !
    !******!

    SUBROUTINE readint(file_id, name, var)
      INTEGER(4)     :: file_id
      CHARACTER(*)   :: name
      INTEGER(4)     :: var
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      !print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    SUBROUTINE readreal(file_id, name, var)
      INTEGER(4)     :: file_id
      CHARACTER(*)   :: name
      REAL(8)        :: var
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      !print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    SUBROUTINE readbool(file_id, name, var)
      INTEGER(4)     :: file_id
      CHARACTER(*)   :: name
      LOGICAL        :: var
      INTEGER        :: value
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, value, data_dims, error)
      IF(value.EQ.1)THEN
        var = .TRUE.
      ELSE
        var = .FALSE.
      ENDIF
      !print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE
  
   SUBROUTINE readstring(file_id, name, var)
      INTEGER(4)                       :: file_id
      CHARACTER(*)                     :: name
      INTEGER(SIZE_T), PARAMETER       :: length = 15
      CHARACTER(LEN=length)            :: var
      INTEGER(HID_T)                   :: dset_id, memtype, filetype, space
      INTEGER                          :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, maxdims
      INTEGER(SIZE_T)                  :: size
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dget_type_f(dset_id, filetype, error)
      CALL h5tget_size_f(filetype, size, error)

      IF(size.GT.length)THEN
        PRINT*,'ERROR: Character LEN is to small'
        STOP
      ENDIF

      CALL h5dget_space_f(dset_id, space, error)
      CALL h5sget_simple_extent_dims_f(space, dims, maxdims, error)
      CALL h5tcopy_f(H5T_FORTRAN_S1, memtype, error)
      CALL h5tset_size_f(memtype, length, error)
      CALL h5dread_f(dset_id, memtype, var, dims, error, space)
      !print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE
    
    SUBROUTINE read_array_complex_dset(file_id, name, var, dims_y)
      COMPLEX(8), DIMENSION(:) :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: dims_y, error
      INTEGER(HID_T) :: dset_id
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      REAL(8), DIMENSION(2,dims_y) :: res
      CALL h5dopen_f(file_id, name, dset_id, error)
      data_dims(1) = 2
      data_dims(2) = dims_y
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      DO i = 1, dims_y
        var(i) = CMPLX(res(1,i),res(2,i))
      END DO
      !print *,name," loaded"
    END SUBROUTINE read_array_complex_dset
    
    SUBROUTINE read_2D_array_complex_dset(file_id, name, var, dims_x, dims_y)
      COMPLEX(8), DIMENSION(:,:) :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: dims_x, dims_y, error
      INTEGER(HID_T) :: dset_id
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      REAL(8), DIMENSION(dims_y,dims_x,2) :: res
      CALL h5dopen_f(file_id, name, dset_id, error)
      data_dims(1) = dims_y
      data_dims(2) = dims_x
      data_dims(3) = 2
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      !print *,"res loaded ok", dims_x, dims_y
      DO i = 1, dims_y
        DO j = 1, dims_x
          var(i,j) = CMPLX(res(j,i,1),res(j,i,2))
        END DO
      END DO
    END SUBROUTINE read_2D_array_complex_dset
    
    SUBROUTINE read_2D_array_complex_dset_slice(file_id, name, var, dims_x, dims_y, slice_x, slice_y, offset_x, offset_y)
      COMPLEX(8), DIMENSION(:,:) :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: dims_x, dims_y, offset_x, offset_y, slice_x, slice_y, rank, error
      INTEGER(HID_T)           :: dset_id,dataspace,memspace
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, slice_dims
      REAL(8), DIMENSION(dims_y,dims_x,2) :: res
      INTEGER(HSIZE_T), DIMENSION(1:3) :: count  ! Size of hyperslab
      INTEGER(HSIZE_T), DIMENSION(1:3) :: offset ! Hyperslab offset
      INTEGER(HSIZE_T), DIMENSION(1:3) :: stride = (/1,1,1/) ! Hyperslab stride
      INTEGER(HSIZE_T), DIMENSION(1:3) :: block = (/1,1,1/)  ! Hyperslab block size
      count = (/2,slice_x,slice_y/)
      offset = (/0,offset_x,offset_y/)
      rank = 3
      slice_dims = (/2,slice_x,slice_y/)
      data_dims = (/2,dims_x,dims_y/)
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dget_space_f(dset_id, dataspace, error)
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, BLOCK)
      CALL h5screate_simple_f(rank, slice_dims, memspace, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, slice_dims, error, memspace, dataspace)
      CALL h5sclose_f(dataspace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      !print *,"res loaded ok", dims_x, dims_y
      DO i = 1, slice_y
        DO j = 1, slice_x
           var(i,j) = CMPLX(res(j,i,1),res(j,i,2))
        END DO
      END DO
      !print *,"e:", var(1:5,1:5)
    END SUBROUTINE read_2D_array_complex_dset_slice


    SUBROUTINE read_2D_array_real_dset_slice(file_id, name, var, dims_x, dims_y, slice_x, slice_y, offset_x, offset_y)
      REAL(8), DIMENSION(:,:) :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: dims_x, dims_y, offset_x, offset_y, slice_x, slice_y, rank, error
      INTEGER(HID_T)           :: dset_id,dataspace,memspace
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, slice_dims, dims, maxdims
      INTEGER(HSIZE_T), DIMENSION(1:2) :: count  ! Size of hyperslab
      INTEGER(HSIZE_T), DIMENSION(1:2) :: offset ! Hyperslab offset
      INTEGER(HSIZE_T), DIMENSION(1:2) :: stride = (/1,1/) ! Hyperslab stride
      INTEGER(HSIZE_T), DIMENSION(1:2) :: block = (/1,1/)  ! Hyperslab block size
      !print *, shape(var)
      count = (/slice_x,slice_y/)
      offset = (/offset_x,offset_y/)
      rank = 2
      slice_dims = (/slice_x,slice_y/)
      data_dims = (/dims_x,dims_y/)
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dget_space_f(dset_id, dataspace, error)
      CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, error)
      !print *, name, dims
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, BLOCK)
      CALL h5screate_simple_f(rank, slice_dims, memspace, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, slice_dims, error, memspace, dataspace)
      CALL h5sclose_f(dataspace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      !print *,"res loaded ok", dims_x, dims_y
    END SUBROUTINE read_2D_array_real_dset_slice

    SUBROUTINE ask_for_size(file_id, name, var)
      INTEGER(4)               :: file_id, dset_id, dataspace
      CHARACTER(*)             :: name
      INTEGER(HSIZE_T), DIMENSION(2) :: var, maxdims
      INTEGER                  :: error
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dget_space_f(dset_id, dataspace, error)
      CALL h5sget_simple_extent_dims_f(dataspace, var, maxdims, error)
      CALL h5sclose_f(dataspace, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE ask_for_size

    !*******!
    ! WRITE !
    !*******!
    
    SUBROUTINE create_scalar_boolean_dset(file_id, name, var)
      INTEGER(4)    :: file_id
      CHARACTER(*)  :: name
      LOGICAL       :: var
      INTEGER       :: value = 0
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      IF (var)THEN
        value = 1
      ENDIF
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, value, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_boolean_dset

    SUBROUTINE create_scalar_int_dset(file_id, name, var)
      INTEGER(4)     :: file_id
      CHARACTER(*)   :: name
      INTEGER(4)     :: var
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_int_dset

    SUBROUTINE create_scalar_real_dset(file_id, name, var)
      INTEGER(4)     :: file_id
      CHARACTER(*)   :: name
      REAL(8)     :: var
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_real_dset
    
    SUBROUTINE create_1D_array_real_dset(file_id, name, var, dims)
      REAL, DIMENSION(:)       :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: error, dims
      INTEGER                  :: rank = 1
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      data_dims = (/dims/)
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_1D_array_real_dset
  
    ! This subroutine supports only arrays of rank 1
    SUBROUTINE create_array_complex_dset(file_id, name, var, dims_y)
      COMPLEX(8), DIMENSION(:) :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: dims_y, error
      INTEGER                  :: rank = 2
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER(HSIZE_T), DIMENSION(2) :: dims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      REAL(8), DIMENSION(2,dims_y) :: res
      dims = (/2, dims_y/)
      DO i = 1, dims_y
        res(1,i) = real(var(i))
        res(2,i) = imag(var(i))
      END DO
      CALL h5screate_simple_f(rank, dims, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      data_dims(1) = 2
      data_dims(2) = dims_y
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_array_complex_dset
    
    SUBROUTINE create_2D_array_complex_dset(file_id, name, var, dims_x, dims_y)
      COMPLEX(8), DIMENSION(:,:) :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: dims_x, dims_y, error
      INTEGER                  :: rank_of_space
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER(HSIZE_T), DIMENSION(3) :: dims
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      REAL(8), DIMENSION(dims_y,dims_x,2) :: res
      rank_of_space = rank(var) + 1
      dims = (/2, dims_x, dims_y/)
      DO i = 1, dims_y
        DO j = 1, dims_x
          res(i,j,1) = real(var(i,j))
          res(i,j,2) = imag(var(i,j))
        END DO
      END DO
      CALL h5screate_simple_f(rank_of_space, dims, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      data_dims(1) = dims_y
      data_dims(2) = dims_x
      data_dims(3) = 2
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_2D_array_complex_dset

    SUBROUTINE create_2D_array_real_dset(file_id, name, var, dims_x, dims_y)
      REAL(8), DIMENSION(:,:)  :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER                  :: dims_x, dims_y, error
      INTEGER                  :: rank = 2
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      data_dims = (/dims_x, dims_y/)
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_2D_array_real_dset
    
    SUBROUTINE create_and_preallocate_2D_array_real_dset(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL, DIMENSION(:,:)   :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, filespace, memspace
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, offset, hyperslab_size
      INTEGER                  :: rank, error
      rank = 2
      CALL h5screate_simple_f(rank, data_dims, filespace, error) ! Create the dataspace for the  dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory (this worker)  
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
    END SUBROUTINE create_and_preallocate_2D_array_real_dset
    
    SUBROUTINE create_3D_array_real_dset(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)   :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, filespace, memspace
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size
      INTEGER                  :: rank, error
      rank = 3
      CALL h5screate_simple_f(rank, data_dims, filespace, error) ! Create the dataspace for the  dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory (this worker)  
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
    END SUBROUTINE create_3D_array_real_dset

    SUBROUTINE h5_add_units_1D(file_id, name, units_value)
      INTEGER(HID_T) :: file_id, dset_id       ! File and Dataset identifiers
      CHARACTER(*)             :: name
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
      CALL h5dopen_f(file_id, name, dset_id, error)
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
      CALL h5dclose_f(dset_id, error)
      ! attributes added
    END SUBROUTINE  h5_add_units_1D
    
    SUBROUTINE create_2D_array_real_dset_p(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL(4), DIMENSION(:,:)   :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, filespace, memspace, h5_parameters
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, offset, hyperslab_size
      INTEGER                  :: rank, error
      rank = 2
      CALL h5screate_simple_f(rank, data_dims, filespace, error) ! Create the dataspace for the  dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory (this worker)  
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, h5_parameters, error) ! Create access parameters for writing
      CALL h5pset_dxpl_mpio_f(h5_parameters, H5FD_MPIO_COLLECTIVE_F, error) ! specify the collective writting
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5_parameters)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
      CALL h5pclose_f(h5_parameters,error)
    END SUBROUTINE create_2D_array_real_dset_p

    ! This subroutine creates 3D array of real numbers parallelly
    ! File has to be opened collectively
    SUBROUTINE create_3D_array_real_dset_p(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)   :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, filespace, memspace, h5_parameters
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size
      INTEGER                  :: rank, error
      rank = 3
      CALL h5screate_simple_f(rank, data_dims, filespace, error) ! Create the dataspace for the  dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory (this worker)  
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, h5_parameters, error) ! Create access parameters for writing
      CALL h5pset_dxpl_mpio_f(h5_parameters, H5FD_MPIO_COLLECTIVE_F, error) ! specify the collective writting
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5_parameters)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
      CALL h5pclose_f(h5_parameters,error)
    END SUBROUTINE create_3D_array_real_dset_p

    SUBROUTINE create_1D_dset_unlimited(file_id, name, var, dim)
      REAL                     :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, dataspace, h5_parameters
      INTEGER                  :: error, dim
      INTEGER(HSIZE_T), DIMENSION(1):: dumh51D, dumh51D2
      dumh51D = (/int(dim,HSIZE_T)/) ! dim
      dumh51D2 = (/H5S_UNLIMITED_F/) ! maxdim
      CALL h5screate_simple_f(1, dumh51D, dataspace, error, dumh51D2 ) ! Create the data space with unlimited dimensions.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5_parameters, error) ! Modify dataset creation properties, i.e. enable chunking
      CALL h5pset_chunk_f(h5_parameters, 1, dumh51D, error) ! enable chunking (1 is the dimension of the dataset)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, dataspace, dset_id, error, h5_parameters) ! create dataset
      dumh51D = (/int(1,HSIZE_T)/) ! dimension of data
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, dumh51D, error)
      CALL h5sclose_f(dataspace, error)
      CALL h5pclose_f(h5_parameters, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE create_1D_dset_unlimited 
    
    SUBROUTINE create_2D_dset_unlimited(file_id, name, var, size)
      REAL,DIMENSION(:,:)        :: var
      REAL,ALLOCATABLE         :: data(:,:)
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, dataspace, h5_parameters
      INTEGER                  :: rank, error, size
      INTEGER, DIMENSION(1)    :: dim
      INTEGER(HSIZE_T), DIMENSION(2):: data_dims, max_dims
      INTEGER(HSIZE_T), DIMENSION(1:2) :: dimsc = (/2,5/)
      ALLOCATE(data(1,size))
      rank = 2
      data_dims = (/1,size/) ! dim
      !Create the data space with unlimited dimensions.
      max_dims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
      CALL h5screate_simple_f(RANK, data_dims, dataspace, error, max_dims)
      !Modify dataset creation properties, i.e. enable chunking
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5_parameters, error)
      CALL h5pset_chunk_f(h5_parameters, RANK, dimsc, error)
      !Create a dataset with 3X3 dimensions using cparms creation propertie .
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, dataspace, &
        dset_id, error, h5_parameters )
      CALL h5sclose_f(dataspace, error)
      !Write data array to dataset
      data_dims(1:2) = (/1,100/)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error)
      CALL h5pclose_f(h5_parameters,error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE create_2D_dset_unlimited 

!!! extendible dataset for single-writter (following the tuto https://portal.hdfgroup.org/display/HDF5/Examples+from+Learning+the+Basics#ExamplesfromLearningtheBasics-changingex https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5/browse/fortran/examples/h5_extend.f90?at=89fbe00dec8187305b518d91c3ddb7d910665f79&raw )

    SUBROUTINE extend_1D_dset_unlimited(file_id, name, var, new_dims, memspace_dims, offset, hyperslab_size)
      REAL                     :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, dataspace
      INTEGER                  :: error
      INTEGER(HSIZE_T), DIMENSION(1):: new_dims, memspace_dims, offset, hyperslab_size
      CALL h5dopen_f(file_id, name, dset_id, error)   !Open the  dataset
      CALL h5dset_extent_f(dset_id, new_dims, error) ! extend the dataset
      CALL h5dget_space_f(dset_id, dataspace, error) ! get the dataspace of the dataset
      CALL h5screate_simple_f (1, memspace_dims, memspace, error) ! create memory space
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, hyperslab_size, error) ! choose the hyperslab in the file
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, hyperslab_size, error, memspace, dataspace) ! wrtiting the data
      CALL h5sclose_f(memspace, error)
      CALL h5sclose_f(dataspace, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE extend_1D_dset_unlimited
    
    SUBROUTINE extend_2D_dset_unlimited(file_id, name, var, new_dims, memspace_dims, offset, hyperslab_size)
      REAL,DIMENSION(:,:)        :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, dataspace
      INTEGER                  :: error, rank
      INTEGER(HSIZE_T), DIMENSION(2):: new_dims, memspace_dims, offset, hyperslab_size
      rank = 2
      CALL h5dopen_f(file_id, name, dset_id, error)   !Open the  dataset
      CALL h5dset_extent_f(dset_id, new_dims, error) ! extend the dataset
      CALL h5dget_space_f(dset_id, dataspace, error) ! get the dataspace of the dataset
      CALL h5screate_simple_f (rank, memspace_dims, memspace, error) ! create memory space
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, hyperslab_size, error) ! choose the hyperslab in the file
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, hyperslab_size, error, memspace, dataspace) ! wrtiting the data
      CALL h5sclose_f(memspace, error)
      CALL h5sclose_f(dataspace, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE extend_2D_dset_unlimited
    
    SUBROUTINE write_hyperslab_to_dset(file_id, name, var, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)       :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, filespace, memspace
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size
      INTEGER                  :: rank, error
      rank = 3
      CALL h5dopen_f(file_id, name, dset_id, error) ! open the dataset (already created)
      CALL h5dget_space_f(dset_id,filespace, error) ! filespace from the dataset (get instead of create)
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory 
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, & 
        file_space_id=filespace,mem_space_id=memspace)
      CALL h5dclose_f(dset_id,error)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
    END SUBROUTINE write_hyperslab_to_dset

    SUBROUTINE write_hyperslab_to_2D_dset(file_id, name, var, offset, hyperslab_size)
      REAL, DIMENSION(:,:)     :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, filespace, memspace
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, offset, hyperslab_size
      INTEGER                  :: rank, error
      rank = 2
      CALL h5dopen_f(file_id, name, dset_id, error) ! open the dataset (already created)
      CALL h5dget_space_f(dset_id,filespace, error) ! filespace from the dataset (get instead of create)
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory 
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, & 
        file_space_id=filespace,mem_space_id=memspace)
      CALL h5dclose_f(dset_id,error)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
    END SUBROUTINE write_hyperslab_to_2D_dset

    SUBROUTINE write_hyperslab_to_dset_p(file_id, name, var, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)       :: var
      INTEGER(4)               :: file_id
      CHARACTER(*)             :: name
      INTEGER(HID_T)           :: dset_id, filespace, memspace, h5_parameters
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size
      INTEGER                  :: rank, error
      rank = 3
      CALL h5dopen_f(file_id, name, dset_id, error) ! open the dataset (already created)
      CALL h5dget_space_f(dset_id,filespace, error) ! filespace from the dataset (get instead of create)
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory 
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, h5_parameters, error) ! Create access parametwers
      CALL h5pset_dxpl_mpio_f(h5_parameters, H5FD_MPIO_COLLECTIVE_F, error) ! collective writting
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, & 
        file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5_parameters)
      CALL h5dclose_f(dset_id,error)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5pclose_f(h5_parameters, error)
    END SUBROUTINE write_hyperslab_to_dset_p
END MODULE hdf5_helper
