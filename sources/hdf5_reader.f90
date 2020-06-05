MODULE hdf5_reader
  CONTAINS
    SUBROUTINE readint(file_id, name, var)
      USE HDF5
      INTEGER(8)     :: file_id
      CHARACTER(*)   :: name
      INTEGER(4)     :: var
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    SUBROUTINE readreal(file_id, name, var)
      USE HDF5
      INTEGER(8)     :: file_id
      CHARACTER(*)   :: name
      REAL(8)        :: var
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
   END SUBROUTINE

   SUBROUTINE readstring(file_id, name, var)
      USE HDF5
      INTEGER(8)     :: file_id
      CHARACTER(*)   :: name
      CHARACTER(15)  :: var
      INTEGER(HID_T) :: dset_id
      INTEGER        :: error
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      CALL h5dopen_f(file_id, name, dset_id, error)
      CALL h5dread_f(dset_id, H5T_FORTRAN_S1, var, data_dims, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

END MODULE hdf5_reader
