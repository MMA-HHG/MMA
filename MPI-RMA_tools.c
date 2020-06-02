#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "mpi.h"


extern int MPEi_CounterFree(MPI_Win counter_win, int keyval, void *attr_val, void *extra_state);


int MPE_Counter_create(MPI_Comm comm, int num, MPI_Win *counter_win) // MPI-3 version
{
static int MPE_COUNTER_KEYVAL = MPI_KEYVAL_INVALID;
int size, rank, lnum, lleft, i, *counterMem=0;
MPI_Aint counterSize;
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size);
lnum = num / size;
lleft = num % size;
if (rank < lleft) lnum++;
counterSize = lnum * sizeof(int);
if (counterSize > 0) {
MPI_Alloc_mem(counterSize, MPI_INFO_NULL, &counterMem);
for (i=0; i<lnum; i++) {counterMem[i] = 0;}
}
/* By using MPI_Alloc_mem first, we ensure that the initial value of the counters are zero. See text */
MPI_Win_create(counterMem, counterSize, sizeof(int),MPI_INFO_NULL, comm, counter_win);
/* Create key if necessary and store the number of counters */
if (MPE_COUNTER_KEYVAL == MPI_KEYVAL_INVALID) {
MPI_Win_create_keyval(MPI_WIN_NULL_COPY_FN, MPEi_CounterFree, &MPE_COUNTER_KEYVAL, NULL);
}
MPI_Win_set_attr(*counter_win, MPE_COUNTER_KEYVAL, (void*)(MPI_Aint)num);
return MPE_COUNTER_KEYVAL;
}


int MPEi_CounterFree(MPI_Win counter_win, int keyval, void *attr_val, void *extra_state)
{
int counter_flag, *counterMem;
MPI_Win_get_attr(counter_win, MPI_WIN_BASE,&counterMem, &counter_flag);
/* Free the memory used by the counter */
if (counter_flag && counterMem)
MPI_Free_mem(counterMem);
return MPI_SUCCESS;
}


int MPE_Counter_nxtval(MPI_Win counterWin, int counterNum, int *value, int MPE_COUNTER_KEYVAL) // MPI-3 version 
{
const int one = 1;
int lrank, flag, size, *attrval;
MPI_Aint lidx;
/* Compute the location of the counter */
MPI_Win_get_attr(counterWin, MPE_COUNTER_KEYVAL, &attrval,&flag);
if (!flag) return -1; /* Error: counterWin not correctly setup */
size = (MPI_Aint)attrval; /* We stored the integer as a pointer */
lrank = counterNum % size;
lidx = counterNum / size;
/* Update and return the counter */
MPI_Win_lock(MPI_LOCK_SHARED, 0, lrank, counterWin);
MPI_Fetch_and_op(&one, value, MPI_INT,lrank, lidx, MPI_SUM, counterWin);
MPI_Win_unlock(lrank, counterWin);
return 0;
}


int MPE_Mutex_acquire(MPI_Win mutex_win, int num, int MPE_MUTEX_KEYVAL) {
int mone = -1, one=1, oldval;
int lrank, flag, size, *attrval;
MPI_Aint lidx;
/* Compute the location of the counter */
MPI_Win_get_attr(mutex_win, MPE_MUTEX_KEYVAL, &attrval, &flag);
if (!flag) return -1; /* Error: counterWin not setup */
size = (int)(MPI_Aint)attrval; /* We stored the integer as a pointer */
lrank = num % size; lidx = num / size;
MPI_Win_lock(MPI_LOCK_SHARED, lrank, 0, mutex_win);
do {
MPI_Fetch_and_op(&one, &oldval, MPI_INT,lrank, lidx, MPI_SUM, mutex_win);
MPI_Win_flush(lrank, mutex_win);
if (oldval == 0) break;
MPI_Accumulate(&mone, 1, MPI_INT, lrank, lidx, 1, MPI_INT,MPI_SUM, mutex_win);
MPI_Win_flush(lrank, mutex_win);
/* We could wait a little bit, depending on oldval */
} while (1);
MPI_Win_unlock(lrank, mutex_win);
return 0;
}


int MPE_Mutex_release(MPI_Win mutex_win, int num, int MPE_MUTEX_KEYVAL)
{
int mone = -1;
int lrank, flag, size, *attrval;
MPI_Aint lidx;
/* Compute the location of the counter */
MPI_Win_get_attr(mutex_win, MPE_MUTEX_KEYVAL, &attrval, &flag);
if (!flag) return -1; /* Error: counterWin setup */
size = (int)(MPI_Aint)attrval; /* We stored the integer as a pointer */
lrank = num % size; lidx = num / size;
MPI_Win_lock(MPI_LOCK_SHARED, lrank, 0, mutex_win);
MPI_Accumulate(&mone, 1, MPI_INT, lrank, lidx, 1, MPI_INT,
MPI_SUM, mutex_win);
MPI_Win_unlock(lrank, mutex_win);
return 0;
}



// collective Mutex: Barrier imposes synchronisation, first worker then acquire mutex, must be called collectivelly
int MPE_Mutex_acquire_collective(MPI_Win mutex_win, int num, int MPE_MUTEX_KEYVAL) {

int numproc, rank;
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &numproc);


if (rank == 0) { // only first worker acquires mutices
    int oldval;
    int lrank, flag, size, *attrval;
    int mnumproc = -numproc;

    MPI_Aint lidx;
    /* Compute the location of the counter */
    MPI_Win_get_attr(mutex_win, MPE_MUTEX_KEYVAL, &attrval, &flag);
    if (!flag) return -1; /* Error: counterWin not setup */
    size = (int)(MPI_Aint)attrval; /* We stored the integer as a pointer */
    lrank = num % size; lidx = num / size;
    MPI_Win_lock(MPI_LOCK_SHARED, lrank, 0, mutex_win);
    do {
    MPI_Fetch_and_op(&numproc, &oldval, MPI_INT,lrank, lidx, MPI_SUM, mutex_win);
    MPI_Win_flush(lrank, mutex_win);
    if (oldval == 0) break;
    MPI_Accumulate(&mnumproc, 1, MPI_INT, lrank, lidx, 1, MPI_INT,MPI_SUM, mutex_win);
    MPI_Win_flush(lrank, mutex_win);
    /* We could wait a little bit, depending on oldval */
    } while (1);
    MPI_Win_unlock(lrank, mutex_win);
}

MPI_Barrier(MPI_COMM_WORLD); // Parallel epoque starts, should be here

return 0;
}