#include "mpi.h"


// !!!!!!!!!!! see: https://www.mcs.anl.gov/research/projects/mpi/usingmpi/examples-advmpi/rma2/

// extern int MPE_MUTEX_KEYVAL;

extern int MPE_COUNTER_KEYVAL;

void MPE_Counter_create( MPI_Comm old_comm, MPI_Win *counter_win )
{
int size, rank, *counter_mem, i, *myval_p;
MPI_Comm_rank(old_comm, &rank);
MPI_Comm_size(old_comm, &size);
if (rank == 0) {
MPI_Alloc_mem( size * sizeof(int), MPI_INFO_NULL,
&counter_mem );
for (i=0; i<size; i++) counter_mem[i] = 0;
MPI_Win_create( counter_mem, size * sizeof(int), sizeof(int),
MPI_INFO_NULL, old_comm, counter_win );
}
else {
MPI_Win_create( NULL, 0, 1, MPI_INFO_NULL, old_comm,
counter_win );
}
/* Create my local counter */
if (MPE_COUNTER_KEYVAL == MPI_KEYVAL_INVALID) {
MPI_Win_create_keyval( MPI_WIN_NULL_COPY_FN,
MPI_WIN_NULL_DELETE_FN,
&MPE_COUNTER_KEYVAL, NULL );
}
myval_p = (int *)malloc( sizeof(int) );
MPI_Win_set_attr( *counter_win, MPE_COUNTER_KEYVAL, myval_p );
}


int MPE_Counter_nxtval( MPI_Win counter_win, int *value )
{
MPI_Group group;
int rank, size, myval, flag, i, *val, one = 1;
MPI_Aint *myval_p;
MPI_Win_get_group( counter_win, &group );
MPI_Group_rank( group, &rank );
MPI_Group_size( group, &size );
MPI_Group_free( &group );
MPI_Win_get_attr( counter_win, MPE_COUNTER_KEYVAL, &myval_p,
&flag );
myval = *myval_p;
val = (int *)malloc ( size * sizeof(int) );
MPI_Win_lock( MPI_LOCK_EXCLUSIVE, 0, 0, counter_win );
for (i=0; i<size; i++) {
if (i == rank)
MPI_Accumulate( &one, 1, MPI_INT, 0, i, 1, MPI_INT,
MPI_SUM, counter_win );
else
MPI_Get( &val[i], 1, MPI_INT, 0, i, 1, MPI_INT,
counter_win );
}
MPI_Win_unlock( 0, counter_win );
/* Add to our contribution */
*myval_p = *myval_p + 1;
/* Compute the overall value.
Storing *myval_p into val[rank] and starting *value at zero
would eliminate the if test */
*value = myval;
for (i=0; i<size; i++) {
if (i != rank) *value = *value + val[i];
}
free ( val );
return 0;
}


// int MPE_Mutex_acquire(MPI_Win mutex_win, int num) {
// int mone = -1, one=1, oldval;
// int lrank, flag, size, *attrval;
// MPI_Aint lidx;
// /* Compute the location of the counter */
// MPI_Win_get_attr(mutex_win, MPE_MUTEX_KEYVAL, &attrval, &flag);
// if (!flag) return -1; /* Error: counterWin not setup */
// size = (int)(MPI_Aint)attrval; /* We stored the integer as a
// pointer */
// lrank = num % size; lidx = num / size;
// MPI_Win_lock(MPI_LOCK_SHARED, lrank, 0, mutex_win);
// do {
// MPI_Fetch_and_op(&one, &oldval, MPI_INT,
// lrank, lidx, MPI_SUM, mutex_win);
// MPI_Win_flush(lrank, mutex_win);
// if (oldval == 0) break;
// MPI_Accumulate(&mone, 1, MPI_INT, lrank, lidx, 1, MPI_INT,
// MPI_SUM, mutex_win);
// MPI_Win_flush(lrank, mutex_win);
// /* We could wait a little bit, depending on oldval */
// } while (1);
// MPI_Win_unlock(lrank, mutex_win);
// return 0;
// }

// int MPE_Mutex_release(MPI_Win mutex_win, int num)
// {
// int mone = -1;
// int lrank, flag, size, *attrval;
// MPI_Aint lidx;
// /* Compute the location of the counter */
// MPI_Win_get_attr(mutex_win, MPE_MUTEX_KEYVAL, &attrval, &flag);
// if (!flag) return -1; /* Error: counterWin setup */
// size = (int)(MPI_Aint)attrval; /* We stored the integer as a
// pointer */
// lrank = num % size; lidx = num / size;
// MPI_Win_lock(MPI_LOCK_SHARED, lrank, 0, mutex_win);
// MPI_Accumulate(&mone, 1, MPI_INT, lrank, lidx, 1, MPI_INT,
// MPI_SUM, mutex_win);
// MPI_Win_unlock(lrank, mutex_win);
// return 0;
// }