MODULE hdf5_helper
  CONTAINS

    !******!
    ! READ !
    !******!

    SUBROUTINE readint(file_id, name, var)
      USE HDF5
      INTEGER(8)     :: file_id
      CHARACTER(*)   :: name
      INTEGER(4)     :: var
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    SUBROUTINE readreal(file_id, name, var)
      USE HDF5
      INTEGER(8)     :: file_id
      CHARACTER(*)   :: name
      REAL(8)        :: var
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    SUBROUTINE readbool(file_id, name, var)
      USE HDF5
      INTEGER(8)     :: file_id
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
      print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE
  
   SUBROUTINE readstring(file_id, name, var)
      USE HDF5
      INTEGER(8)                       :: file_id
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
      print *, name, var
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    !*******!
    ! WRITE !
    !*******!
    
    SUBROUTINE create_scalar_boolean_dset(file_id, name, var)
      USE HDF5
      INTEGER(8)    :: file_id
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
      USE HDF5
      INTEGER(8)     :: file_id
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
      USE HDF5
      INTEGER(8)     :: file_id
      CHARACTER(*)   :: name
      REAL(8)     :: var
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_real_dset
  
    ! This subroutine supports only arrays of rank 1
    SUBROUTINE create_array_complex_dset(file_id, name, var, dims_y)
      USE HDF5
      COMPLEX(8), DIMENSION(:) :: var
      INTEGER(8)               :: file_id
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
      USE HDF5
      COMPLEX(8), DIMENSION(:,:) :: var
      INTEGER(8)               :: file_id
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

    ! DO NOT USE
    SUBROUTINE create_array_real_dset(file_id, name, var, dims_y)
      USE HDF5
      COMPLEX(8), DIMENSION(:) :: var
      INTEGER(8)               :: file_id
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
    END SUBROUTINE create_array_real_dset

END MODULE hdf5_helper
