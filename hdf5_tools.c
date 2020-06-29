#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "mpi.h"


void readreal(hid_t file_id, char *dset_name, herr_t *h5error, double *value)
{
  hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); // open dataset
  hid_t datatype  = H5Dget_type(dset_id);
  *h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
  *h5error = H5Dclose(dset_id);
}

void readint(hid_t file_id, char *dset_name, herr_t *h5error, int *value)
{
  hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); // open dataset
  hid_t datatype  = H5Dget_type(dset_id);
  *h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
  *h5error = H5Dclose(dset_id);
}


double * readreal1Darray_fort(hid_t file_id, char *dset_name, herr_t *h5error, int *N_points) // fort is for the extra diemnsion due to fortran
{
	hid_t dset_id = H5Dopen2 (file_id, dset_name, H5P_DEFAULT); // open dataset	     
	hid_t dspace_id = H5Dget_space (dset_id); // Get the dataspace ID     
	const int ndims = H5Sget_simple_extent_ndims(dspace_id); // number of dimensions in the tgrid
	hsize_t dims[ndims]; // we need the size to allocate tgrid for us
	H5Sget_simple_extent_dims(dspace_id, dims, NULL); // get dimensions

	hid_t datatype  = H5Dget_type(dset_id);     // we get the type of data (SINGLE, DOUBLE, etc. from HDF5)


	// double tgrid[dims[0]]; // allocate the grid
  double *array = malloc((int)dims[0],sizeof(double));
  
	/*see https://stackoverflow.com/questions/10575544/difference-between-array-type-and-array-allocated-with-malloc
	      https://stackoverflow.com/questions/216259/is-there-a-max-array-length-limit-in-c/216731#216731  */
	*h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array); // read the grid
	// if ( ( comment_operation == 1 ) && ( myrank == 0 ) ){printf("(t_init,t_end) = (%e,%e) \n",tgrid[0],tgrid[dims[0]-1]);}
	*h5error = H5Dclose(dset_id);

  *N_points = (int)dims[0];
  return array;
}

//int linkexists(hid_t file_id, char *link_name, herr_t *h5error)
//{
//  hid_t dset_id = H5Lexists (file_id, link_name, H5P_DEFAULT); // open dataset
//  hid_t datatype  = H5Dget_type(dset_id);
//  *h5error = H5Dread(dset_id,  datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
//  *h5error = H5Dclose(dset_id);
//}

void addone(int *val ){*val=*val+1;} // to test pointers



  //   SUBROUTINE readint(file_id, name, var)
  //     USE HDF5
  //     INTEGER(4)     :: file_id
  //     CHARACTER(*)   :: name
  //     INTEGER(4)     :: var
  //     INTEGER(HID_T) :: dset_id
  //     INTEGER        :: error
  //     INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
  //     !print *, name, var
  //     CALL h5dclose_f(dset_id, error)
  //   END SUBROUTINE

  //   SUBROUTINE readreal(file_id, name, var)
  //     USE HDF5
  //     INTEGER(4)     :: file_id
  //     CHARACTER(*)   :: name
  //     REAL(8)        :: var
  //     INTEGER(HID_T) :: dset_id
  //     INTEGER        :: error
  //     INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
  //     !print *, name, var
  //     CALL h5dclose_f(dset_id, error)
  //   END SUBROUTINE

  //   SUBROUTINE readbool(file_id, name, var)
  //     USE HDF5
  //     INTEGER(4)     :: file_id
  //     CHARACTER(*)   :: name
  //     LOGICAL        :: var
  //     INTEGER        :: value
  //     INTEGER(HID_T) :: dset_id
  //     INTEGER        :: error
  //     INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, value, data_dims, error)
  //     IF(value.EQ.1)THEN
  //       var = .TRUE.
  //     ELSE
  //       var = .FALSE.
  //     ENDIF
  //     !print *, name, var
  //     CALL h5dclose_f(dset_id, error)
  //   END SUBROUTINE
  
  //  SUBROUTINE readstring(file_id, name, var)
  //     USE HDF5
  //     INTEGER(4)                       :: file_id
  //     CHARACTER(*)                     :: name
  //     INTEGER(SIZE_T), PARAMETER       :: length = 15
  //     CHARACTER(LEN=length)            :: var
  //     INTEGER(HID_T)                   :: dset_id, memtype, filetype, space
  //     INTEGER                          :: error
  //     INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, maxdims
  //     INTEGER(SIZE_T)                  :: size
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     CALL h5dget_type_f(dset_id, filetype, error)
  //     CALL h5tget_size_f(filetype, size, error)

  //     IF(size.GT.length)THEN
  //       PRINT*,'ERROR: Character LEN is to small'
  //       STOP
  //     ENDIF

  //     CALL h5dget_space_f(dset_id, space, error)
  //     CALL h5sget_simple_extent_dims_f(space, dims, maxdims, error)
  //     CALL h5tcopy_f(H5T_FORTRAN_S1, memtype, error)
  //     CALL h5tset_size_f(memtype, length, error)
  //     CALL h5dread_f(dset_id, memtype, var, dims, error, space)
  //     !print *, name, var
  //     CALL h5dclose_f(dset_id, error)
  //   END SUBROUTINE
    
  //   SUBROUTINE read_array_complex_dset(file_id, name, var, dims_y)
  //     USE HDF5
  //     COMPLEX(8), DIMENSION(:) :: var
  //     INTEGER(4)               :: file_id
  //     CHARACTER(*)             :: name
  //     INTEGER                  :: dims_y, error
  //     INTEGER(HID_T) :: dset_id
  //     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
  //     REAL(8), DIMENSION(2,dims_y) :: res
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     data_dims(1) = 2
  //     data_dims(2) = dims_y
  //     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     DO i = 1, dims_y
  //       var(i) = CMPLX(res(1,i),res(2,i))
  //     END DO
  //     !print *,name," loaded"
  //   END SUBROUTINE read_array_complex_dset
    
  //   SUBROUTINE read_2D_array_complex_dset(file_id, name, var, dims_x, dims_y)
  //     USE HDF5
  //     COMPLEX(8), DIMENSION(:,:) :: var
  //     INTEGER(4)               :: file_id
  //     CHARACTER(*)             :: name
  //     INTEGER                  :: dims_x, dims_y, error
  //     INTEGER(HID_T) :: dset_id
  //     INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
  //     REAL(8), DIMENSION(dims_y,dims_x,2) :: res
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     data_dims(1) = dims_y
  //     data_dims(2) = dims_x
  //     data_dims(3) = 2
  //     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     !print *,"res loaded ok", dims_x, dims_y
  //     DO i = 1, dims_y
  //       DO j = 1, dims_x
  //         var(i,j) = CMPLX(res(j,i,1),res(j,i,2))
  //       END DO
  //     END DO
  //   END SUBROUTINE read_2D_array_complex_dset
    
  //   SUBROUTINE read_2D_array_complex_dset_slice(file_id, name, var, dims_x, dims_y, slice_x, slice_y, offset_x, offset_y)
  //     USE HDF5
  //     COMPLEX(8), DIMENSION(:,:) :: var
  //     INTEGER(4)               :: file_id
  //     CHARACTER(*)             :: name
  //     INTEGER                  :: dims_x, dims_y, offset_x, offset_y, slice_x, slice_y, rank, error
  //     INTEGER(HID_T)           :: dset_id,dataspace,memspace
  //     INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, slice_dims
  //     REAL(8), DIMENSION(dims_y,dims_x,2) :: res
  //     INTEGER(HSIZE_T), DIMENSION(1:3) :: count  ! Size of hyperslab
  //     INTEGER(HSIZE_T), DIMENSION(1:3) :: offset ! Hyperslab offset
  //     INTEGER(HSIZE_T), DIMENSION(1:3) :: stride = (/1,1,1/) ! Hyperslab stride
  //     INTEGER(HSIZE_T), DIMENSION(1:3) :: block = (/1,1,1/)  ! Hyperslab block size
  //     count = (/2,slice_x,slice_y/)
  //     offset = (/0,offset_x,offset_y/)
  //     rank = 3
  //     slice_dims = (/2,slice_x,slice_y/)
  //     data_dims = (/2,dims_x,dims_y/)
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     CALL h5dget_space_f(dset_id, dataspace, error)
  //     CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, BLOCK)
  //     CALL h5screate_simple_f(rank, slice_dims, memspace, error)
  //     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, res, slice_dims, error, memspace, dataspace)
  //     CALL h5sclose_f(dataspace, error)
  //     CALL h5sclose_f(memspace, error)
  //     CALL h5dclose_f(dset_id, error)
  //     !print *,"res loaded ok", dims_x, dims_y
  //     DO i = 1, slice_y
  //       DO j = 1, slice_x
  //          var(i,j) = CMPLX(res(j,i,1),res(j,i,2))
  //       END DO
  //     END DO
  //     !print *,"e:", var(1:5,1:5)
  //   END SUBROUTINE read_2D_array_complex_dset_slice


  //   SUBROUTINE read_2D_array_real_dset_slice(file_id, name, var, dims_x, dims_y, slice_x, slice_y, offset_x, offset_y)
  //     USE HDF5
  //     REAL(8), DIMENSION(:,:) :: var
  //     INTEGER(4)               :: file_id
  //     CHARACTER(*)             :: name
  //     INTEGER                  :: dims_x, dims_y, offset_x, offset_y, slice_x, slice_y, rank, error
  //     INTEGER(HID_T)           :: dset_id,dataspace,memspace
  //     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims, slice_dims
  //     REAL(8), DIMENSION(dims_y,dims_x) :: res
  //     INTEGER(HSIZE_T), DIMENSION(1:2) :: count  ! Size of hyperslab
  //     INTEGER(HSIZE_T), DIMENSION(1:2) :: offset ! Hyperslab offset
  //     INTEGER(HSIZE_T), DIMENSION(1:2) :: stride = (/1,1/) ! Hyperslab stride
  //     INTEGER(HSIZE_T), DIMENSION(1:2) :: block = (/1,1/)  ! Hyperslab block size
  //     count = (/slice_x,slice_y/)
  //     offset = (/offset_x,offset_y/)
  //     rank = 2
  //     slice_dims = (/slice_x,slice_y/)
  //     data_dims = (/dims_x,dims_y/)
  //     CALL h5dopen_f(file_id, name, dset_id, error)
  //     CALL h5dget_space_f(dset_id, dataspace, error)
  //     CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, BLOCK)
  //     CALL h5screate_simple_f(rank, slice_dims, memspace, error)
  //     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var, slice_dims, error, memspace, dataspace)
  //     CALL h5sclose_f(dataspace, error)
  //     CALL h5sclose_f(memspace, error)
  //     CALL h5dclose_f(dset_id, error)
  //     !print *,"res loaded ok", dims_x, dims_y
  //   END SUBROUTINE read_2D_array_real_dset_slice
  //   !*******!
  //   ! WRITE !
  //   !*******!
    
  //   SUBROUTINE create_scalar_boolean_dset(file_id, name, var)
  //     USE HDF5
  //     INTEGER(4)    :: file_id
  //     CHARACTER(*)  :: name
  //     LOGICAL       :: var
  //     INTEGER       :: value = 0
  //     INTEGER(HID_T) :: dset_id, dataspace_id
  //     INTEGER        :: error
  //     INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
  //     IF (var)THEN
  //       value = 1
  //     ENDIF
  //     CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
  //     CALL h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
  //     CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, value, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     CALL h5sclose_f(dataspace_id, error)
  //   END SUBROUTINE create_scalar_boolean_dset

  //   SUBROUTINE create_scalar_int_dset(file_id, name, var)
  //     USE HDF5
  //     INTEGER(4)     :: file_id
  //     CHARACTER(*)   :: name
  //     INTEGER(4)     :: var
  //     INTEGER(HID_T) :: dset_id, dataspace_id
  //     INTEGER        :: error
  //     INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
  //     CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
  //     CALL h5dcreate_f(file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dset_id, error)
  //     CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, var, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     CALL h5sclose_f(dataspace_id, error)
  //   END SUBROUTINE create_scalar_int_dset

  //   SUBROUTINE create_scalar_real_dset(file_id, name, var)
  //     USE HDF5
  //     INTEGER(4)     :: file_id
  //     CHARACTER(*)   :: name
  //     REAL(8)     :: var
  //     INTEGER(HID_T) :: dset_id, dataspace_id
  //     INTEGER        :: error
  //     INTEGER(HSIZE_T), DIMENSION(1:1) :: data_dims
  //     CALL h5screate_f(H5S_SCALAR_F, dataspace_id, error)
  //     CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
  //     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     CALL h5sclose_f(dataspace_id, error)
  //   END SUBROUTINE create_scalar_real_dset
  
  //   ! This subroutine supports only arrays of rank 1
  //   SUBROUTINE create_array_complex_dset(file_id, name, var, dims_y)
  //     USE HDF5
  //     COMPLEX(8), DIMENSION(:) :: var
  //     INTEGER(4)               :: file_id
  //     CHARACTER(*)             :: name
  //     INTEGER                  :: dims_y, error
  //     INTEGER                  :: rank = 2
  //     INTEGER(HID_T) :: dset_id, dataspace_id
  //     INTEGER(HSIZE_T), DIMENSION(2) :: dims
  //     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
  //     REAL(8), DIMENSION(2,dims_y) :: res
  //     dims = (/2, dims_y/)
  //     DO i = 1, dims_y
  //       res(1,i) = real(var(i))
  //       res(2,i) = imag(var(i))
  //     END DO
  //     CALL h5screate_simple_f(rank, dims, dataspace_id, error)
  //     CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
  //     data_dims(1) = 2
  //     data_dims(2) = dims_y
  //     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     CALL h5sclose_f(dataspace_id, error)
  //   END SUBROUTINE create_array_complex_dset
    
  //   SUBROUTINE create_2D_array_complex_dset(file_id, name, var, dims_x, dims_y)
  //     USE HDF5
  //     COMPLEX(8), DIMENSION(:,:) :: var
  //     INTEGER(4)               :: file_id
  //     CHARACTER(*)             :: name
  //     INTEGER                  :: dims_x, dims_y, error
  //     INTEGER                  :: rank_of_space
  //     INTEGER(HID_T) :: dset_id, dataspace_id
  //     INTEGER(HSIZE_T), DIMENSION(3) :: dims
  //     INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
  //     REAL(8), DIMENSION(dims_y,dims_x,2) :: res
  //     rank_of_space = rank(var) + 1
  //     dims = (/2, dims_x, dims_y/)
  //     DO i = 1, dims_y
  //       DO j = 1, dims_x
  //         res(i,j,1) = real(var(i,j))
  //         res(i,j,2) = imag(var(i,j))
  //       END DO
  //     END DO
  //     CALL h5screate_simple_f(rank_of_space, dims, dataspace_id, error)
  //     CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
  //     data_dims(1) = dims_y
  //     data_dims(2) = dims_x
  //     data_dims(3) = 2
  //     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, res, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     CALL h5sclose_f(dataspace_id, error)
  //   END SUBROUTINE create_2D_array_complex_dset

  //   SUBROUTINE create_2D_array_real_dset(file_id, name, var, dims_x, dims_y)
  //     USE HDF5
  //     REAL(8), DIMENSION(:,:)  :: var
  //     INTEGER(4)               :: file_id
  //     CHARACTER(*)             :: name
  //     INTEGER                  :: dims_x, dims_y, error
  //     INTEGER                  :: rank = 2
  //     INTEGER(HID_T) :: dset_id, dataspace_id
  //     INTEGER(HSIZE_T), DIMENSION(2) :: dims
  //     INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
  //     data_dims = (/dims_x, dims_y/)
  //     CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
  //     CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
  //     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, var, data_dims, error)
  //     CALL h5dclose_f(dset_id, error)
  //     CALL h5sclose_f(dataspace_id, error)
  //   END SUBROUTINE create_2D_array_real_dset


  //   ! Remove from final version
  //   SUBROUTINE test(process_rank)
  //       USE HDF5
  //       INTEGER       :: process_rank
  //       CHARACTER(LEN=10), PARAMETER :: filename = "sds.h5"  ! File name
  //       CHARACTER(LEN=10), PARAMETER :: name = "testfield"
  //       INTEGER(HID_T) :: file_id       ! File identifier 
  //       INTEGER(HID_T) :: plist_id      ! Property list identifier 
  //       INTEGER(HID_T) :: dset_id, dataspace_id
  //       INTEGER        :: error
  //       REAL(8), DIMENSION(2) :: testing_field
  //       INTEGER :: comm, info, rank
  //       INTEGER(HSIZE_T), DIMENSION(2) :: data_dims 
  //       print *,"my_rank = ", process_rank
  //       testing_field(1:process_rank+1) = process_rank
  //       rank = 2
  //       comm = MPI_COMM_WORLD
  //       info = MPI_INFO_NULL
  //       data_dims(1) = 1
  //       data_dims(2) = 2
  //       ! Initialize FORTRAN predefined datatypes
  //       CALL h5open_f(error) 

  //       ! Setup file access property list with parallel I/O access.
  //       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  //       CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

  //       ! H5F_ACC_EXCL_F - fail if file already exists, ..._TRUNC_F - overwrites existing file
  //       print *, "proc ", process_rank
          
  //       CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
  //       print *,"file created"
  //       CALL h5screate_simple_f(rank, data_dims, dataspace_id, error)
  //       print *,"dataspce created"
  //       CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dset_id, error)
  //       print *,"written ",testing_field
  //       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, testing_field, data_dims, error)
  //       CALL h5dclose_f(dset_id, error)
  //       CALL h5sclose_f(dataspace_id, error)
  //       CALL h5fclose_f(file_id, error)
          
  //       ! Close property list
  //       CALL h5pclose_f(plist_id, error)
  //       ! Close FORTRAN interface
  //       CALL h5close_f(error)
  //   END SUBROUTINE test
