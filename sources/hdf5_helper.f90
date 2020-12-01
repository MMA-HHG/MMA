! This is the collection of HDF5 operations on the files used in the code
! This module was designed by Jakub Jelinek & Jan Vabek

MODULE hdf5_helper
  USE HDF5
  ! USE H5LT ! lite
  ! Create an interface for reading a dset. It includes most of the reading subroutines, the ones that are uniquely identifiable by their parameters.
  INTERFACE read_dset
    PROCEDURE readint, readreal, readbool, readstring, read_array_complex_dset, read_2D_array_complex_dset, &
     read_array_real_dset, read_2D_array_real_dset, read_2D_array_complex_dset_slice, read_2D_array_real_dset_slice
  END INTERFACE
  ! Create an interface for creating a dset. It includes most of the writing subroutines, the ones that are uniquely identifiable by their parameters.
  INTERFACE create_dset
    PROCEDURE create_scalar_boolean_dset, create_scalar_int_dset, create_scalar_real_dset, create_1D_array_real_dset, &
      create_1D_array_real8_dset, create_1D_array_complex_dset, create_2D_array_complex_dset, create_2D_array_real_dset, &
      create_scalar_string_dset
  END INTERFACE
  ! Check if a given dataset exists. It reads the value in such a case, it cretes this dataset and savs the input value otherwise.
  INTERFACE save_or_replace
    PROCEDURE save_or_replace_real8, save_or_replace_int, save_or_replace_bool, save_or_replace_string
  END INTERFACE
  CONTAINS

    !******!
    ! READ !
    !******!

    ! This subroutine reads an integer from scalar dataset.
    SUBROUTINE readint(file_id, name, var)
      INTEGER(4)     :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)   :: name ! name of the dataset to be read
      INTEGER(4)     :: var ! variable for that stores the value read from the dataset
      INTEGER(HID_T) :: dset_id ! id of the dataset to be opened
      INTEGER        :: error ! error variable
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims ! dimension of the dataset
      ! Open the dataset and store dataset id in dset_id
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Read the value from dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    ! This subroutine reads an integer from scalar dataset
    SUBROUTINE readreal(file_id, name, var)
      INTEGER(4)     :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)   :: name ! name of the dataset to be read
      REAL(8)        :: var ! variable for that stores the value read from the dataset
      INTEGER(HID_T) :: dset_id ! id of the dataset to be opened
      INTEGER        :: error ! error variable
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims ! dimension of the dataset
      ! Open the dataset and store dataset id in dset_id
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Read the value from dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE

    ! This subroutine reads an integer (either 0 or 1) and evaluates it as a boolean
    SUBROUTINE readbool(file_id, name, var)
      INTEGER(4)     :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)   :: name ! name of the dataset to be read
      LOGICAL        :: var ! variable for that stores the value read from the dataset
      INTEGER(HID_T) :: dset_id ! id of the dataset to be opened
      INTEGER        :: error ! error variable
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims ! dimension of the dataset
      INTEGER        :: value ! Variable for the value to be analysed
      ! Open the dataset and store dataset id in dset_id
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Read the value from dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, value, data_dims, error)
      ! If value is 1, store .TRUE. in var, store .FALSE. otherwise
      IF(value.EQ.1)THEN
        var = .TRUE.
      ELSE
        var = .FALSE.
      ENDIF
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE
  
    ! This subroutine reads a string from scalar dataset
    SUBROUTINE readstring(file_id, name, var)
      INTEGER(4)                       :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)                     :: name ! name of the dataset to be read
      INTEGER(SIZE_T), PARAMETER       :: length = 15 ! length of the string, it can be enlarged if needed
      CHARACTER(LEN=length)            :: var ! variable for that stores the value read from the dataset
      INTEGER(HID_T)                   :: dset_id, memtype, filetype, space ! define necessary identifiers
      INTEGER                          :: error ! error variable
      INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, maxdims ! dimension of the dataset and the 
      INTEGER(SIZE_T)                  :: size ! real size of the string in the dataset
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Get the dataspace of the dataset
      CALL h5dget_type_f(dset_id, filetype, error)
      ! Get the real size of the string
      CALL h5tget_size_f(filetype, size, error) 
      ! Check if the length of the variable is large enough
      IF(size.GT.length)THEN
        PRINT*,'ERROR: Character LEN is to small'
        STOP
      ENDIF
      ! Get the dataspace of the dataset
      CALL h5dget_space_f(dset_id, space, error)
      ! Get the size and maximum sizes of each dimension of a dataspace
      CALL h5sget_simple_extent_dims_f(space, dims, maxdims, error)
      ! The following returns a modifiable transient datatype which is a copy of H5T_FORTRAN_S1
      CALL h5tcopy_f(H5T_FORTRAN_S1, memtype, error)
      ! Set size of the dataspaceto the length
      CALL h5tset_size_f(memtype, length, error)
      ! Read the value from the dataset
      CALL h5dread_f(dset_id, memtype, var, dims, error, space)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE
    
    ! This subroutine reads a 1D real dataset
    SUBROUTINE read_array_real_dset(file_id, name, var, dims_y)
      REAL(8), DIMENSION(:) :: var ! variable for that stores the value read from the dataset
      INTEGER(4)               :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)             :: name ! name of the dataset to be read
      INTEGER                  :: dims_y, error ! dims_y stores the size of the dataset to be read, error stores error messages
      INTEGER(HID_T)           :: dset_id ! dataset id
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims ! data dimensions array
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Initialize data dimensions variable 
      data_dims(1) = dims_y
      ! Read the dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE read_array_real_dset
    
    ! This subroutine reads a 2D real dataset
    SUBROUTINE read_2D_array_real_dset(file_id, name, var, dims_x, dims_y)
      REAL(8), DIMENSION(:,:) :: var ! variable for that stores the value read from the dataset
      INTEGER(4)               :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)             :: name ! name of the dataset to be read
      INTEGER                  :: dims_x, dims_y, error ! dims_x and dims_y store the size of the dataset to be read, error stores error messages
      INTEGER(HID_T)           :: dset_id ! dataset id
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims ! data dimensions array
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Initialize data dimensions variable 
      data_dims(1) = dims_x
      data_dims(2) = dims_y
      ! Read the dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE read_2D_array_real_dset
    
    ! This subroutine reads a 1D complex dataset
    SUBROUTINE read_array_complex_dset(file_id, name, var, dims_y)
      COMPLEX(8), DIMENSION(:) :: var ! variable for that stores the value read from the dataset
      INTEGER(4)               :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)             :: name ! name of the dataset to be read
      INTEGER                  :: dims_y, error ! dims_y stores the size of the dataset to be read, error stores error messages
      INTEGER(HID_T)           :: dset_id ! dataset id
      ! As HDF5 does not support complex numbers, the dataset is a 2 by dims_y dataset of reals
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims ! data dimensions array has two values because of that
      REAL(8), DIMENSION(2,dims_y) :: res ! this variable stores the actual data from the dataset
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Initialize data dimensions variable 
      data_dims(1) = 2
      data_dims(2) = dims_y
      ! Read the dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
      ! Use the temporary variable to fill in the output variable
      DO i = 1, dims_y
        var(i) = CMPLX(res(1,i),res(2,i))
      END DO
    END SUBROUTINE read_array_complex_dset
    
    ! This subroutine reads a 2D complex dataset
    SUBROUTINE read_2D_array_complex_dset(file_id, name, var, dims_x, dims_y)
      COMPLEX(8), DIMENSION(:,:) :: var ! variable for that stores the value read from the dataset
      INTEGER(4)               :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)             :: name ! name of the dataset to be read
      INTEGER                  :: dims_x, dims_y, error ! dims_x and dims_y store the size of the dataset to be read, error stores error messages
      INTEGER(HID_T)           :: dset_id ! dataset id
      ! As HDF5 does not support complex numbers, the dataset is a 2 by dims_x by dims_y dataset of reals
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims ! data dimensions array has three values because of that
      REAL(8), DIMENSION(dims_y,dims_x,2) :: res ! this variable stores the actual data from the dataset
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Initialize data dimensions variable 
      data_dims(1) = dims_y
      data_dims(2) = dims_x
      data_dims(3) = 2
      ! Read the dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
      ! Use the temporary variable to fill in the output variable      
      DO i = 1, dims_y
        DO j = 1, dims_x
          var(i,j) = CMPLX(res(j,i,1),res(j,i,2))
        END DO
      END DO
    END SUBROUTINE read_2D_array_complex_dset
    
    ! This subroutine reads a part of the 2D complex dataset
    SUBROUTINE read_2D_array_complex_dset_slice(file_id, name, var, dims_x, dims_y, slice_x, slice_y, offset_x, offset_y)
      COMPLEX(8), DIMENSION(:,:) :: var ! variable for that stores the value read from the dataset
      INTEGER(4)                 :: file_id ! the id of the already opened h5 file or a group in a file
      CHARACTER(*)               :: name ! name of the dataset to be read
      ! dims_x and dims_y store the size of the dataset to be read
      ! offset_x and offset_y store the offset from the left top corner of the dataset
      ! rank of the dataset, error stores error messages
      INTEGER                    :: dims_x, dims_y, offset_x, offset_y, slice_x, slice_y, rank, error
      INTEGER(HID_T)             :: dset_id,dataspace,memspace ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, slice_dims ! dimensions variables
      REAL(8), DIMENSION(dims_y,dims_x,2) :: res ! temporary variable for storing read data
      INTEGER(HSIZE_T), DIMENSION(1:3) :: count  ! Size of hyperslab
      INTEGER(HSIZE_T), DIMENSION(1:3) :: offset ! Hyperslab offset
      INTEGER(HSIZE_T), DIMENSION(1:3) :: stride = (/1,1,1/) ! Hyperslab stride
      INTEGER(HSIZE_T), DIMENSION(1:3) :: block = (/1,1,1/)  ! Hyperslab block size
      ! Initialize arrays
      count = (/2,slice_x,slice_y/)
      offset = (/0,offset_x,offset_y/)
      slice_dims = (/2,slice_x,slice_y/)
      data_dims = (/2,dims_x,dims_y/)
      ! Initialize rank of the dataset
      rank = 3
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Get dataspace of the dataset
      CALL h5dget_space_f(dset_id, dataspace, error)
      ! Select the wanted hyperslab
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, BLOCK)
      ! Create a dataspace with a size of the slice
      CALL h5screate_simple_f(rank, slice_dims, memspace, error)
      ! Read the dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, slice_dims, error, memspace, dataspace)
      ! Close dataspaces and dataset
      CALL h5sclose_f(dataspace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      ! Use the temporary variable to fill in the output variable
      DO i = 1, slice_y
        DO j = 1, slice_x
           var(i,j) = CMPLX(res(j,i,1),res(j,i,2))
        END DO
      END DO
    END SUBROUTINE read_2D_array_complex_dset_slice

    ! This subroutine reads a part of the 2D complex dataset
    SUBROUTINE read_2D_array_real_dset_slice(file_id, name, var, dims_x, dims_y, slice_x, slice_y, offset_x, offset_y)
      REAL(8), DIMENSION(:,:) :: var ! variable for that stores the value read from the dataset
      INTEGER(4)               :: file_id ! identifier of the file
      CHARACTER(*)             :: name ! dataset name
      ! dims_x and dims_y store the size of the dataset to be read
      ! offset_x and offset_y store the offset from the left top corner of the dataset
      ! rank of the dataset, error stores error messages
      INTEGER                  :: dims_x, dims_y, offset_x, offset_y, slice_x, slice_y, rank, error
      ! Create necessary identifiers
      INTEGER(HID_T)           :: dset_id,dataspace,memspace
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, slice_dims, dims, maxdims ! Create necessary arrays
      INTEGER(HSIZE_T), DIMENSION(1:2) :: count  ! Size of hyperslab
      INTEGER(HSIZE_T), DIMENSION(1:2) :: offset ! Hyperslab offset
      INTEGER(HSIZE_T), DIMENSION(1:2) :: stride = (/1,1/) ! Hyperslab stride
      INTEGER(HSIZE_T), DIMENSION(1:2) :: block = (/1,1/)  ! Hyperslab block size
      ! Initialize arrays and rank
      count = (/slice_x,slice_y/)
      offset = (/offset_x,offset_y/)
      rank = 2
      slice_dims = (/slice_x,slice_y/)
      data_dims = (/dims_x,dims_y/)
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Get the dataspace
      CALL h5dget_space_f(dset_id, dataspace, error)
      ! Get the size and maximum sizes of each dimension of a dataspace
      CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, error)
      ! Select hyperslab from the dataset
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, BLOCK)
      ! Create a dataspace with needed dimensions
      CALL h5screate_simple_f(rank, slice_dims, memspace, error)
      ! Read the dataset
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, slice_dims, error, memspace, dataspace)
      ! Close dataspaces and dataset
      CALL h5sclose_f(dataspace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE read_2D_array_real_dset_slice
    
    ! This subroutine returns a size of a 1D dataset with given name
    SUBROUTINE ask_for_size_1D(file_id, name, var)
      INTEGER(4)               :: file_id, dset_id, dataspace ! necessary identifiers
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HSIZE_T), DIMENSION(1) :: tmp_var, maxdims ! necessary arrays
      INTEGER(4)               :: var ! variable to be written to
      INTEGER                  :: error ! error message
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Get the dataspace
      CALL h5dget_space_f(dset_id, dataspace, error)
      ! Get the dimensions of the dataspace
      CALL h5sget_simple_extent_dims_f(dataspace, tmp_var, maxdims, error)
      var = INT(tmp_var(1), 4)
      ! Close the dataspace
      CALL h5sclose_f(dataspace, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE ask_for_size_1D

    ! This subroutine returns a size of a 2D dataset with given name
    SUBROUTINE ask_for_size_2D(file_id, name, var)
      INTEGER(4)               :: file_id, dset_id, dataspace ! necessary identifiers
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HSIZE_T), DIMENSION(2) :: var, maxdims ! necessary arrays
      INTEGER                  :: error ! error message
      ! Open the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)
      ! Get the dataspace
      CALL h5dget_space_f(dset_id, dataspace, error)
      ! Get the dimensions of the dataspace
      CALL h5sget_simple_extent_dims_f(dataspace, var, maxdims, error)
      ! Close the dataspace
      CALL h5sclose_f(dataspace, error)
      ! Close the dataset
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE ask_for_size_2D

    !*******!
    ! WRITE !
    !*******!
    
    ! This subroutine creates a scalar dataset of type boolean (as HDF5 does not support booleans, the real part is an integer)
    SUBROUTINE create_scalar_boolean_dset(file_id, name, var)
      INTEGER(4)    :: file_id ! identifier of the file or group in which the dataset is supposed to be created
      CHARACTER(*)  :: name ! name of the dataset to be created
      LOGICAL       :: var ! variable to be written into the dataset
      INTEGER       :: value = 0 ! temporary value, which stores the integer equivalent of the boolean variable
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER        :: error ! error stores error messages of HDF5 interface
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims ! dimensions of the dataset
      ! evaluate var and find integer equivalent - .FALSE. -> 0, .TRUE. -> 1
      IF (var) THEN
        value = 1
      ENDIF
      ! create scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      ! create integer dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
      ! write value to dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, value, data_dims, error)
      ! close dataset and dataspace, terminate the subroutine
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_boolean_dset

    ! This subroutine creates a scalar dataset of type integer
    SUBROUTINE create_scalar_int_dset(file_id, name, var)
      INTEGER(4)     :: file_id ! identifier of the file or group in which the dataset is supposed to be created 
      CHARACTER(*)   :: name ! name of the dataset to be created
      INTEGER(4)     :: var ! variable to be written into the dataset
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER        :: error ! error stores error messages of HDF5 interface
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims ! dimensions of the dataset
      ! create scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      ! create integer dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
      ! write var to dset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      ! close dataset and dataspace, terminate subroutine
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_int_dset

    ! This subroutine creates a scalar dataset of type real (double precision)
    SUBROUTINE create_scalar_real_dset(file_id, name, var)
      INTEGER(4)     :: file_id ! identifier of the file or group in which the dataset is supposed to be created 
      CHARACTER(*)   :: name ! name of the dataset to be created
      REAL(8)        :: var ! variable to be written into the dataset
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER        :: error ! error stores error messages of HDF5 interface
      INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims ! dimensions of the dataset
      ! create scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      ! create double precision dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      ! write var to dset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      ! close dataset and dataspace, terminate subroutine
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_real_dset
   
    ! This subroutine creates one dimensional dataset of type real (double precision)
    SUBROUTINE create_1D_array_real_dset(file_id, name, var, dims)
      REAL, DIMENSION(:)       :: var ! variable to be written into the dataset
      INTEGER(4)               :: file_id ! file or group identifier to write the dataset to
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER                  :: error, dims ! error stores error messages of the HDF5 interface, dims stores the size of the dataset
      INTEGER                  :: rank = 1 ! rank of the dataset
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims ! dimensions array
      data_dims = (/dims/) ! assign value to dimensions array
      ! create dataspace of rank 1 and size data_dims
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      ! create dataset of type double precision
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      ! write to dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error)
      ! close dataset and dataspace, terminate subroutine
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_1D_array_real_dset

    ! This subroutine creates one dimensional dataset of type real (double precision)
    SUBROUTINE create_1D_array_real8_dset(file_id, name, var, dims)
      REAL(8), DIMENSION(:)       :: var ! variable to be written into the dataset
      INTEGER(4)               :: file_id ! file or group identifier to write the dataset to
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER                  :: error, dims ! error stores error messages of the HDF5 interface, dims stores the size of the dataset
      INTEGER                  :: rank = 1 ! rank of the dataset
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims ! dimensions array
      data_dims = (/dims/) ! assign value to dimensions array
      ! create dataspace of rank 1 and size data_dims
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      ! create dataset of type double precision
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      ! write to dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      ! close dataset and dataspace, terminate subroutine
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_1D_array_real8_dset
  
    ! This subroutine creates one dimensional dataset of type complex (as HDF5 does not support complex numbers, the real value is
    ! double precision and the dataset is of size 2 by dims)
    SUBROUTINE create_1D_array_complex_dset(file_id, name, var, dims_y)
      COMPLEX(8), DIMENSION(:) :: var ! complex variable to be written to the dataset
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER                  :: dims_y, error ! dims_y stores the size of the dataset, error stores the error messages from HDF5 interface
      INTEGER                  :: rank = 2 ! rank of the dataset - it is rank of the complex array + 1 as the real and imaginary part have to be seprated
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims ! dimensions of the dataset
      REAL(8), DIMENSION(2,dims_y) :: res ! temporary variable which stores real and imag part separately
      ! assign dimensions a value 
      data_dims(1) = 2
      data_dims(2) = dims_y
      ! separate the real and imaginary part of each complex number and store them into temporary variable
      DO i = 1, dims_y
        res(1,i) = real(var(i))
        res(2,i) = imag(var(i))
      END DO
      ! create dataspace of rank 2 and size dims
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      ! create dataset of type double precision
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      ! write data to dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      ! close the dataset and dataspace
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_1D_array_complex_dset
   
    ! This subroutine creates two dimensional dataset of type complex (as HDF5 does not support complex numbers, the real value is
    ! double precision and the dataset is of size 2 by dims_x by dims_y)
    SUBROUTINE create_2D_array_complex_dset(file_id, name, var, dims_x, dims_y)
      COMPLEX(8), DIMENSION(:,:) :: var ! complex variable to be written to the dataset
      INTEGER(4)               :: file_id ! file or group id to create the dataset in
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER                  :: dims_x, dims_y, error ! dims_x and dims_y are passed to initialize the dimensions, error stores error messages of the HDF5 interface
      INTEGER                  :: rank ! rank of the dataspace
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims ! data dimensions array
      REAL(8), DIMENSION(dims_y,dims_x,2) :: res ! temporary variable storing the real and imag part separately
      rank = 3 ! initialize rank
      ! fill in the data_dims array
      data_dims = (/2, dims_x, dims_y/)
      ! separate the real and imaginary part of each complex number and store them into temporary variable
      DO i = 1, dims_y
        DO j = 1, dims_x
          res(i,j,1) = real(var(i,j))
          res(i,j,2) = imag(var(i,j))
        END DO
      END DO
      ! create dataspace of rank 3 and size data_dims
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      ! create dataset of type double precision
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      ! write to dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
      ! close the dataset and dataspace, terminate subroutine
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_2D_array_complex_dset

    ! this subroutine creates two dimensional real dataset
    SUBROUTINE create_2D_array_real_dset(file_id, name, var, dims_x, dims_y)
      REAL(8), DIMENSION(:,:)  :: var ! variable to be written into the dataset
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER                  :: dims_x, dims_y, error ! dims_x and dims_y store the size of the dataset, error stores error messages of the HDF5 interface
      INTEGER                  :: rank = 2 ! rank of the dataset
      INTEGER(HID_T) :: dset_id, dataspace_id ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims ! data dimensions array
      data_dims = (/dims_x, dims_y/) ! initialize data dimensions array with size of the dataset
      ! create dataspace of rank 2 and size of the size of the dataset
      CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
      ! create a dataset of type double precision
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
      ! write var to the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
      ! close the dataset and dataspace, terminate the subroutine
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_2D_array_real_dset
   
    ! this subroutine creates a two dimensional dataset and writes to part of it
    SUBROUTINE create_and_preallocate_2D_array_real_dset(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL, DIMENSION(:,:)     :: var ! variable to be written to the first part of the dataset  
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HID_T)           :: dset_id, filespace, memspace ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, offset, hyperslab_size ! dimensions of data, data offset and the hyperslab size
      INTEGER                  :: rank, error ! rank stores the rank of the dataset, error stores error messages of HDF5 interface
      rank = 2 ! assign a value to rank
      ! Create the dataspace for the dataset
      CALL h5screate_simple_f(rank, data_dims, filespace, error)
      ! Create the dataset collectivelly
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)
      ! Dataset dimensions in memory (this worker)
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error)
      ! Select hyperslab to be written to
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      ! Write to the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace)
      ! Close both dataspaces and the dataset, terminate subroutine
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
    END SUBROUTINE create_and_preallocate_2D_array_real_dset
   
    ! This subroutine creates a three dimensional dataset and writes to part of it
    SUBROUTINE create_and_preallocate_3D_array_real_dset(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)   :: var ! variable to be written to the first part of the dataset 
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HID_T)           :: dset_id, filespace, memspace ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size ! dimensions of data, data offset and the hyperslab size
      INTEGER                  :: rank, error ! rank is the rank of the dataset, error stores the error messages of HDF5 interface
      rank = 3 ! initialize rank to be 3
      ! create dataspace with the rank of 3 and size of data_dims
      CALL h5screate_simple_f(rank, data_dims, filespace, error) ! Create the dataspace for the  dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory (this worker)  
      ! select the hyperslab wanted
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      ! write to the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace)
      ! close the dataspaces and the dataset
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
    END SUBROUTINE create_and_preallocate_3D_array_real_dset

    ! This subroutine adds units to a dataset
    SUBROUTINE h5_add_units_1D(file_id, name, units_value)
      INTEGER(HID_T) :: file_id, dset_id       ! File and Dataset identifiers
      CHARACTER(*)   :: name ! name of the dataset
      character(*)   :: units_value ! units to be written to the attribute
      INTEGER        :: error ! error stores error messages of the HDF5 interface
      INTEGER(HSIZE_T), DIMENSION(1):: dumh51D ! temporary variable
      ! ATTRIBUTES OF HDF5 DATASETS
      CHARACTER(LEN=5), PARAMETER :: aname = "units"   ! Attribute name
      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      CHARACTER(LEN=20), DIMENSION(1) ::  attr_data  ! Attribute data (solve the length by a correct allocation)
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
    
    ! This subroutine creates two dimensional real dataset for parallel writting
    SUBROUTINE create_2D_array_real_dset_p(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL(4), DIMENSION(:,:)   :: var ! variable to be written to the dataset
      INTEGER(4)                :: file_id ! file or group identifier 
      CHARACTER(*)              :: name ! name of the dataset
      INTEGER(HID_T)            :: dset_id, filespace, memspace, h5_parameters ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, offset, hyperslab_size ! dimensions of data, offset of the data and the hyperslab size
      INTEGER                   :: rank, error ! rank is the rank of the dataset, error stores error messages of the HDF5 interface
      rank = 2 ! assign rank a value of 2
      CALL h5screate_simple_f(rank, data_dims, filespace, error) ! Create the dataspace for the  dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory (this worker)  
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error) ! select hyperslab
      CALL h5pcreate_f(H5P_DATASET_XFER_F, h5_parameters, error) ! Create access parameters for writing
      CALL h5pset_dxpl_mpio_f(h5_parameters, H5FD_MPIO_COLLECTIVE_F, error) ! specify the collective writting
      ! write to the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5_parameters)
      ! close dataspaces and the dataset
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
      CALL h5pclose_f(h5_parameters,error)
    END SUBROUTINE create_2D_array_real_dset_p

    ! This subroutine creates 3D array of real numbers parallelly
    ! File has to be opened collectively
    SUBROUTINE create_3D_array_real_dset_p(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)   :: var ! variable to be written to the dataset
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HID_T)           :: dset_id, filespace, memspace, h5_parameters ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size ! dimensions of data, offset of the data and the hyperslab size
      INTEGER                  :: rank, error ! rank is the rank of the dataset, error stores error messages of the HDF5 interface
      rank = 3 ! rank of the dataset
      CALL h5screate_simple_f(rank, data_dims, filespace, error) ! Create the dataspace for the  dataset
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, filespace, dset_id, error)  ! create the dataset collectivelly
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory (this worker) 
      ! select the correct hyperslab
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, h5_parameters, error) ! Create access parameters for writing
      CALL h5pset_dxpl_mpio_f(h5_parameters, H5FD_MPIO_COLLECTIVE_F, error) ! specify the collective writting
      ! write to the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, &
        file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5_parameters)
      ! Close the dataspaces, the dataset and the h5 parameter
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5dclose_f(dset_id,error)
      CALL h5pclose_f(h5_parameters,error)
    END SUBROUTINE create_3D_array_real_dset_p

    ! This subroutine creates a one dimensional unlimited dataset
    SUBROUTINE create_1D_dset_unlimited(file_id, name, var, dim)
      REAL,DIMENSION(:)        :: var ! variable which is supposed to be initialized with the dataset
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset to be created
      INTEGER(HID_T)           :: dset_id, dataspace, h5_parameters ! necessary identifiers
      INTEGER                  :: error, dim ! error stores error messages of the HDF5 interface, dim is the size of the var array
      INTEGER(HSIZE_T), DIMENSION(1):: dumh51D, dumh51D2 ! these are temporary variables for storing the dimensions of the dataset
      dumh51D = (/int(dim,HSIZE_T)/) ! dim
      dumh51D2 = (/H5S_UNLIMITED_F/) ! maxdim
      CALL h5screate_simple_f(1, dumh51D, dataspace, error, dumh51D2 ) ! Create the data space with unlimited dimensions.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, h5_parameters, error) ! Modify dataset creation properties, i.e. enable chunking
      CALL h5pset_chunk_f(h5_parameters, 1, dumh51D, error) ! enable chunking (1 is the dimension of the dataset)
      CALL h5dcreate_f(file_id, name, H5T_NATIVE_REAL, dataspace, dset_id, error, h5_parameters) ! create dataset
      ! write to dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, dumh51D, error)
      ! close dataspace, h5 parameter and the dataset
      CALL h5sclose_f(dataspace, error)
      CALL h5pclose_f(h5_parameters, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE create_1D_dset_unlimited 
    
    ! this subroutine creates a two dimensional unlimited dataset of type real
    SUBROUTINE create_2D_dset_unlimited(file_id, name, var, size)
      REAL,DIMENSION(:,:)      :: var ! variable to be written to the dataset
      REAL,ALLOCATABLE         :: data(:,:) ! 
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
        dset_id, error, h5_parameters)
      ! close the dataspace
      CALL h5sclose_f(dataspace, error)
      ! Write data array to dataset
      data_dims = (/1,size/)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error)
      ! close the h5 parameters and the dataset
      CALL h5pclose_f(h5_parameters,error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE create_2D_dset_unlimited 

!!! extendible dataset for single-writter (following the tuto https://portal.hdfgroup.org/display/HDF5/Examples+from+Learning+the+Basics#ExamplesfromLearningtheBasics-changingex https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5/browse/fortran/examples/h5_extend.f90?at=89fbe00dec8187305b518d91c3ddb7d910665f79&raw )

    ! this subroutine extends the unlimited dataset with given name
    SUBROUTINE extend_1D_dset_unlimited(file_id, name, var, new_dims, memspace_dims, offset, hyperslab_size)
      REAL,DIMENSION(:)        :: var ! value which is supposed to be appended to the dataset
      INTEGER(4)               :: file_id ! file or group identifier containing the dataset
      CHARACTER(*)             :: name ! name of the dataset to be extended
      INTEGER(HID_T)           :: dset_id, dataspace ! necessary identifiers
      INTEGER                  :: error ! error stores error messages of the HDF5 interface
      INTEGER(HSIZE_T), DIMENSION(1):: new_dims, memspace_dims, offset, hyperslab_size ! new dimensions of the dataset, of the dataspace, the offset and the hyperslab size
      CALL h5dopen_f(file_id, name, dset_id, error)   !Open the  dataset
      CALL h5dset_extent_f(dset_id, new_dims, error) ! extend the dataset
      CALL h5dget_space_f(dset_id, dataspace, error) ! get the dataspace of the dataset
      CALL h5screate_simple_f (1, memspace_dims, memspace, error) ! create memory space
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, hyperslab_size, error) ! choose the hyperslab in the file
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, hyperslab_size, error, memspace, dataspace) ! wrtiting the data
      ! close the dataspaces and the dataset
      CALL h5sclose_f(memspace, error)
      CALL h5sclose_f(dataspace, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE extend_1D_dset_unlimited
    
    ! this subroutine extends a two dimensional unlimited dataset
    SUBROUTINE extend_2D_dset_unlimited(file_id, name, var, new_dims, memspace_dims, offset, hyperslab_size)
      REAL,DIMENSION(:,:)      :: var ! variable which is supposed to be appended to the dataset
      INTEGER(4)               :: file_id ! file or group identifier containing the dataset
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HID_T)           :: dset_id, dataspace ! necessary identifiers
      INTEGER                  :: error, rank ! error stores error messages of the HDF5 interface, rank stores dataset rank
      INTEGER(HSIZE_T), DIMENSION(2):: new_dims, memspace_dims, offset, hyperslab_size ! new dimensions of the dataset, of the dataspace, the offset and the hyperslab size

      rank = 2 ! rank of the dataset
      CALL h5dopen_f(file_id, name, dset_id, error)   !Open the  dataset
      CALL h5dset_extent_f(dset_id, new_dims, error) ! extend the dataset
      CALL h5dget_space_f(dset_id, dataspace, error) ! get the dataspace of the dataset
      CALL h5screate_simple_f (rank, memspace_dims, memspace, error) ! create memory space
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, hyperslab_size, error) ! choose the hyperslab in the file
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, hyperslab_size, error, memspace, dataspace) ! wrtiting the data
      ! close the dataset and the dataspaces
      CALL h5sclose_f(memspace, error)
      CALL h5sclose_f(dataspace, error)
      CALL h5dclose_f(dset_id, error)
    END SUBROUTINE extend_2D_dset_unlimited
    
    ! this subroutine writes a hyperslab to a preallocated one dimensional dataset
    SUBROUTINE write_hyperslab_to_dset(file_id, name, var, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)   :: var ! variable to written
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HID_T)           :: dset_id, filespace, memspace ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size ! dimensions of data, offset and hyperslab to be written
      INTEGER                  :: rank, error ! rank is the rank of the dataset, error stores error messages of the HDF5 interface
      rank = 3 ! give rank a value
      CALL h5dopen_f(file_id, name, dset_id, error) ! open the dataset (already created)
      CALL h5dget_space_f(dset_id,filespace, error) ! filespace from the dataset (get instead of create)
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory 
      ! select a hyperslab
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      ! write to the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, & 
        file_space_id=filespace,mem_space_id=memspace)
      ! close the dataset and both dataspaces
      CALL h5dclose_f(dset_id,error)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
    END SUBROUTINE write_hyperslab_to_dset

    ! this subroutine writes a hyperslab to a preallocated two dimensional dataset
    SUBROUTINE write_hyperslab_to_2D_dset(file_id, name, var, offset, hyperslab_size)
      REAL, DIMENSION(:,:)     :: var ! variable to be written to the hyperslab
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HID_T)           :: dset_id, filespace, memspace ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, offset, hyperslab_size ! dimensions of the data, the data offset and the size of the hyperslab
      INTEGER                  :: rank, error ! rank is the rank of the dataset (2), error stores error messages of the HDF5 interface
      rank = 2 ! assign a value to rank
      CALL h5dopen_f(file_id, name, dset_id, error) ! open the dataset (already created)
      CALL h5dget_space_f(dset_id,filespace, error) ! filespace from the dataset (get instead of create)
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory 
      ! select the correct hyperslab
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      ! write to the chosen hyperslab of the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, & 
        file_space_id=filespace,mem_space_id=memspace)
      ! close the dataspaces and the dataset
      CALL h5dclose_f(dset_id,error)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
    END SUBROUTINE write_hyperslab_to_2D_dset

    ! this subroutine writes a hyperslab to a three dimensional dataset parallelly
    SUBROUTINE write_hyperslab_to_dset_p(file_id, name, var, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)   :: var ! variable to be written to the hyperslab
      INTEGER(4)               :: file_id ! file or group identifier
      CHARACTER(*)             :: name ! name of the dataset
      INTEGER(HID_T)           :: dset_id, filespace, memspace, h5_parameters ! necessary identifiers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, offset, hyperslab_size ! dimensions of the data, the data offset and the size of the hyperslab
      INTEGER                  :: rank, error ! rank is the rank of the dataset, error stores error messages of the HDF5 interface
      rank = 3 ! assign a value to rank
      CALL h5dopen_f(file_id, name, dset_id, error) ! open the dataset (already created)
      CALL h5dget_space_f(dset_id,filespace, error) ! filespace from the dataset (get instead of create)
      CALL h5screate_simple_f(rank, hyperslab_size, memspace, error) ! dataset dimensions in memory 
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, hyperslab_size, error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, h5_parameters, error) ! Create access parametwers
      CALL h5pset_dxpl_mpio_f(h5_parameters, H5FD_MPIO_COLLECTIVE_F, error) ! collective writting
      ! write to the dataset
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, var, data_dims, error, & 
        file_space_id=filespace,mem_space_id=memspace,xfer_prp = h5_parameters)
      ! close the dataset, the dataspaces and the h5 parameters
      CALL h5dclose_f(dset_id,error)
      CALL h5sclose_f(filespace,error)
      CALL h5sclose_f(memspace,error)
      CALL h5pclose_f(h5_parameters, error)
    END SUBROUTINE write_hyperslab_to_dset_p


    !*************!
    ! READ & WRITE!
    !*************!
    SUBROUTINE save_or_replace_real8(file_id, name, var, error, units_in) ! allow to add units, reading currently unsupported (h5aexists_f)
      REAL(8)                 :: var
      INTEGER(4)              :: file_id ! file or group identifier
      CHARACTER(*)            :: name
      INTEGER                 :: error
      CHARACTER(*), OPTIONAL  :: units_in

      LOGICAL         :: exists_dataset

      CALL  h5lexists_f(file_id,name,exists_dataset,error)

      IF (exists_dataset) THEN
        CALL read_dset(file_id, name, var)
      ELSE
        CALL create_dset(file_id, name, var)
        IF (PRESENT(units_in)) CALL h5_add_units_1D(file_id, name, units_in)
      ENDIF
    END SUBROUTINE save_or_replace_real8

    SUBROUTINE save_or_replace_int(file_id, name, var, error,units_in)
      INTEGER                 :: var
      INTEGER(4)              :: file_id ! file or group identifier
      CHARACTER(*)            :: name
      INTEGER                 :: error
      CHARACTER(*), OPTIONAL  :: units_in

      LOGICAL         :: exists_dataset

      CALL  h5lexists_f(file_id,name,exists_dataset,error)

      IF (exists_dataset) THEN
        CALL read_dset(file_id, name, var)
      ELSE
        CALL create_dset(file_id, name, var)
        IF (PRESENT(units_in)) CALL h5_add_units_1D(file_id, name, units_in)
      ENDIF
    END SUBROUTINE save_or_replace_int

    SUBROUTINE save_or_replace_bool(file_id, name, var, error, units_in)
      LOGICAL                 :: var
      INTEGER(4)              :: file_id ! file or group identifier
      CHARACTER(*)            :: name
      INTEGER                 :: error
      CHARACTER(*), OPTIONAL  :: units_in

      LOGICAL         :: exists_dataset

      CALL  h5lexists_f(file_id,name,exists_dataset,error)

      IF (exists_dataset) THEN
        CALL read_dset(file_id, name, var)
      ELSE
        CALL create_dset(file_id, name, var)
        IF (PRESENT(units_in)) CALL h5_add_units_1D(file_id, name, units_in)
      ENDIF
    END SUBROUTINE save_or_replace_bool


    ! This subroutine creates a scalar dataset of type character
    SUBROUTINE create_scalar_string_dset(file_id, name, var)
      USE ISO_C_BINDING ! saving string using c-pointer
      INTEGER(4)     :: file_id ! identifier of the file or group in which the dataset is supposed to be created 
      CHARACTER(*)   :: name ! name of the dataset to be created
      CHARACTER(*)   :: var ! variable to be written into the dataset
      INTEGER(HID_T) :: dset_id, dataspace_id, type_id ! necessary identifiers
      INTEGER        :: error ! error stores error messages of HDF5 interface

      TYPE(C_PTR) :: f_ptr
      TYPE(C_PTR), DIMENSION(1), TARGET :: wdata 
      CHARACTER(len = 50, kind=c_char), DIMENSION(1), TARGET  :: c_var ! adjust length...


      ! INTEGER(HSIZE_T), DIMENSION(1):: dumh51D ! temporary variable
      ! ! ATTRIBUTES OF HDF5 DATASETS
      ! CHARACTER(LEN=5), PARAMETER :: aname = "units"   ! Attribute name
      ! INTEGER(HID_T) :: attr_id       ! Attribute identifier
      ! INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      ! INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
      ! CHARACTER(LEN=20), DIMENSION(1) ::  attr_data  ! Attribute data (solve the length by a correct allocation)
      ! ! add attributes ( https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/fortran/examples/h5_crtatt.f90 )
      ! CALL h5dopen_f(file_id, name, dset_id, error)
      ! dumh51D = (/int(1,HSIZE_T)/) ! attribute dimension
      ! CALL h5screate_simple_f(1, dumh51D, aspace_id, error) ! Create scalar data space for the attribute. 1 stands for the rank
      ! CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) ! Create datatype for the attribute.
      ! CALL h5tset_size_f(atype_id, int(10,HSIZE_T), error) ! 10 is attribute length	
      ! CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error) ! Create dataset attribute.
      ! dumh51D = (/int(1,HSIZE_T)/) ! dimension of attributes
      ! attr_data(1) = units_value 
      ! CALL h5awrite_f(attr_id, atype_id, attr_data, dumh51D, error)
      ! CALL h5aclose_f(attr_id, error)  ! Close the attribute.
      ! CALL h5tclose_f(atype_id, error)  ! Close the attribute datatype.
      ! CALL h5sclose_f(aspace_id, error) ! Terminate access to the attributes data space.


      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, error)
      CALL h5tset_size_f(type_id, int(len(TRIM(var)),HSIZE_T), error)
      CALL h5screate_simple_f(1, (/ INT(1,HSIZE_T) /), dataspace_id, error)

      CALL h5dcreate_f(file_id, name, type_id, dataspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, type_id, (/ TRIM(var) /), (/ INT(1,HSIZE_T) /), error)

      CALL h5dclose_f(dset_id , error)
      CALL h5sclose_f(dataspace_id, error)
      CALL H5Tclose_f(type_id, error)
      

      ! c_var(1) = TRIM(var)//C_NULL_CHAR

      !   ! Create file and memory datatypes.  For this example we will save
      ! ! the strings as C variable length strings, H5T_STRING is defined
      ! ! as a variable length string.
      ! !
      ! CALL H5Tcopy_f(H5T_STRING, type_id, error)
      ! !
      ! ! Create dataspace.
      ! !
      ! ! CALL h5screate_simple_f(1, dims, dataspace_id, error)
      ! CALL h5screate_simple_f(1, (/ INT(1,HSIZE_T) /), dataspace_id, error)
      ! !
      ! ! Create the dataset and write the variable-length string data to
      ! ! it.
      ! !
      ! CALL h5dcreate_f(file_id, name, type_id, dataspace_id, dset_id, error)

      ! f_ptr = C_LOC(c_var(1))
      ! CALL h5dwrite_f(dset_id, type_id, f_ptr, error)
      ! !
      ! ! Close and release resources.
      ! !
      ! CALL h5dclose_f(dset_id , error)
      ! CALL h5sclose_f(dataspace_id, error)
      ! CALL H5Tclose_f(type_id, error)
      ! ! CALL h5fclose_f(file , hdferr)

      ! ! create scalar dataspace
      ! CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
      ! ! create integer dataset
      ! CALL h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
      ! ! write var to dset
      ! CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
      ! ! close dataset and dataspace, terminate subroutine
      ! CALL h5dclose_f(dset_id, error)
      ! CALL h5sclose_f(dataspace_id, error)
    END SUBROUTINE create_scalar_string_dset



    SUBROUTINE save_or_replace_string(file_id, name, var, error, units_in)
      CHARACTER(*)            :: var
      INTEGER(4)              :: file_id ! file or group identifier
      CHARACTER(*)            :: name
      INTEGER                 :: error
      CHARACTER(*), OPTIONAL  :: units_in

      LOGICAL         :: exists_dataset

      CALL  h5lexists_f(file_id,name,exists_dataset,error)

      IF (exists_dataset) THEN
        CALL read_dset(file_id, name, var)
      ELSE
        ! https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5-examples/browse/1_10/FORTRAN/H5T/h5ex_t_stringC_F03.f90
        ! print *, 'not writing strings not implemented yet'


        CALL create_scalar_string_dset(file_id, name, var)
        IF (PRESENT(units_in)) CALL h5_add_units_1D(file_id, name, units_in)
      ENDIF
    END SUBROUTINE save_or_replace_string
END MODULE hdf5_helper
