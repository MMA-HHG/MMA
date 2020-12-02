! This is the collection of HDF5 operations on the files used in the code
! This module was designed by Jakub Jelinek & Jan Vabek

MODULE hdf5_helper_parallel
  USE HDF5
  ! USE H5LT ! lite
  ! Create an interface for reading a dset. It includes most of the reading subroutines, the ones that are uniquely identifiable by their parameters.

  IMPLICIT NONE

      ! This subroutine creates two dimensional real dataset for parallel writting
    SUBROUTINE create_2D_array_real_dset_p(file_id, name, var, data_dims, offset, hyperslab_size)
      REAL(4), DIMENSION(:,:)   :: var ! variable to be written to the dataset
      INTEGER(HID_T)            :: file_id ! file or group identifier 
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
      INTEGER(HID_T)           :: file_id ! file or group identifier
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
  
    ! this subroutine writes a hyperslab to a three dimensional dataset parallelly
    SUBROUTINE write_hyperslab_to_dset_p(file_id, name, var, offset, hyperslab_size)
      REAL, DIMENSION(:,:,:)   :: var ! variable to be written to the hyperslab
      INTEGER(HID_T)           :: file_id ! file or group identifier
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

END MODULE hdf5_helper_parallel

MODULE hdf5_helper
  USE hdf5_helper_serial
  USE hdf5_helper_parallel
  IMPLICIT NONE
END MODULE hdf5_helper
