MODULE hdf5_helper
  CONTAINS
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
    
    ! The group in which you want to create the dataset has to be passed with group_id
    SUBROUTINE create_scalar_int_dset(group_id, name, var)
      USE HDF5
      INTEGER(8)     :: group_id
      CHARACTER(*)   :: name
      INTEGER(4)     :: var
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      CALL h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      CALL h5sclose_f(dataspace_id, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE create_scalar_int_dset

    SUBROUTINE create_scalar_real_dset(group_id, name, var)
      USE HDF5
      INTEGER(8)     :: group_id
      CHARACTER(*)   :: name
      REAL(8)     :: var
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      CALL h5dcreate_f(group_id, name, H5T_NATIVE_REAL, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error)
      CALL h5sclose_f(dataspace_id, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE create_scalar_real_dset

    SUBROUTINE create_array_int_dset(group_id, name, var, dims_x, dims_y)
      USE HDF5
      INTEGER(8)     :: group_id
      CHARACTER(*)   :: name
      INTEGER(4)     :: var
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error, rank = 2
      INTEGER        :: dims_x, dims_y
      INTEGER(HSIZE_T), DIMENSION(dims_x:dims_y) :: data_dims
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      CALL h5dcreate_f(group_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      CALL h5sclose_f(dataspace_id, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE create_array_int_dset

    SUBROUTINE create_array_real_dset(group_id, name, var, dims_x, dims_y)
      USE HDF5
      INTEGER(8)     :: group_id
      CHARACTER(*)   :: name
      REAL(8)        :: var
      INTEGER(HID_T) :: dset_id, dataspace_id
      INTEGER        :: error, rank = 2
      INTEGER        :: dims_x, dims_y
      INTEGER(HSIZE_T), DIMENSION(dims_x:dims_y) :: data_dims
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      CALL h5dcreate_f(group_id, name, H5T_NATIVE_REAL, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error)
      CALL h5sclose_f(dataspace_id, error)
      CALL h5dclose_f(dset_id, error)
  END SUBROUTINE create_array_real_dset  
   
END MODULE hdf5_helper
