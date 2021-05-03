void MPE_Counter_create(MPI_Comm, int, MPI_Win *);
int MPE_Counter_nxtval(MPI_Win, int, int *);
int MPE_Mutex_acquire(MPI_Win, int, int);
int MPE_Mutex_release(MPI_Win, int, int);

double t_mpi[10]; 


