MODULE hdf5_helper
  USE HDF5
  ! Create an interface for reading a dset. It includes most of the reading subroutines, the ones that are uniquely identifiable by their parameters.
  INTERFACE read_dset
    PROCEDURE readint, readreal, readbool, readstring, read_array_complex_dset, read_2D_array_complex_dset, &
     read_array_real_dset, read_2D_array_real_dset, read_2D_array_complex_dset_slice, read_2D_array_real_dset_slice
  END INTERFACE
  ! Create an interface for creating a dset. It includes most of the writing subroutines, the ones that are uniquely identifiable by their parameters.
  INTERFACE create_dset
    PROCEDURE create_scalar_boolean_dset, create_scalar_int_dset, create_scalar_real_dset, create_1D_array_real_dset, &
      create_array_complex_dset, create_2D_array_complex_dset, create_2D_array_real_dset
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
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, slice_dims ! dimesnions variables
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
      REAL,DIMENSION(:)        :: var
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
      REAL,DIMENSION(:)        :: var
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
