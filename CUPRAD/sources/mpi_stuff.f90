MODULE mpi_stuff
  !include 'mpif.h'
  use mpi
  INTEGER(4) my_rank,num_proc,ierr,MPI_SUBARRAY,MPI_SUBARRAY_TRANSPOSED
  INTEGER(4) status(MPI_STATUS_SIZE)
  CHARACTER(3) ip

  REAL(8) start_time_MPI

END MODULE mpi_Stuff