/**
 * @file tools_MPI-RMA.c
 * @brief Contains MPI routines for MPI tasks management.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "tools_MPI-RMA.h"

int MPE_COUNTER_KEYVAL = MPI_KEYVAL_INVALID;
extern int MPEi_CounterFree(MPI_Win counter_win, int keyval, void *attr_val, void *extra_state);


void MPE_Counter_create(MPI_Comm comm, int num, MPI_Win *counter_win) // MPI-3 version
{
    //static int MPE_COUNTER_KEYVAL = MPI_KEYVAL_INVALID;
    //int MPE_COUNTER_KEYVAL = MPI_KEYVAL_INVALID;
    int size, rank, lnum, lleft, i, *counterMem=0;
    MPI_Aint counterSize;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    lnum = num / size;
    lleft = num % size;

    if (rank < lleft) 
        lnum++;
    counterSize = lnum * sizeof(int);

    if (counterSize > 0) {
        MPI_Alloc_mem(counterSize, MPI_INFO_NULL, &counterMem);
        for (i=0; i<lnum; i++) 
            counterMem[i] = 0;
    }

    printf("bc keyval %i\n",MPE_COUNTER_KEYVAL);
    /* By using MPI_Alloc_mem first, we ensure that the initial value of the counters are zero. See text */
    MPI_Win_create(counterMem, counterSize, sizeof(int),MPI_INFO_NULL, comm, counter_win);
    /* Create key if necessary and store the number of counters */
    if (MPE_COUNTER_KEYVAL == MPI_KEYVAL_INVALID) {
        MPI_Win_create_keyval(MPI_WIN_NULL_COPY_FN, MPEi_CounterFree, &MPE_COUNTER_KEYVAL, NULL);
    }

    MPI_Win_set_attr(*counter_win, MPE_COUNTER_KEYVAL, (void*)(MPI_Aint)num);
    printf("ac keyval %i\n",MPE_COUNTER_KEYVAL);
}


int MPEi_CounterFree(MPI_Win counter_win, int keyval, void *attr_val, void *extra_state) {
    int counter_flag, *counterMem;
    MPI_Win_get_attr(counter_win, MPI_WIN_BASE,&counterMem, &counter_flag);
    /* Free the memory used by the counter */
    if (counter_flag && counterMem)
        MPI_Free_mem(counterMem);
    return MPI_SUCCESS;
}


int MPE_Counter_nxtval(MPI_Win counterWin, int counterNum, int *value) // MPI-3 version 
{
    const int one = 1;
    int lrank, flag, size, *attrval;
    MPI_Aint lidx;
    /* Compute the location of the counter */
    MPI_Win_get_attr(counterWin, MPE_COUNTER_KEYVAL, &attrval,&flag);
    if (!flag) 
        return -1; /* Error: counterWin not correctly setup */
    size = (MPI_Aint)attrval; /* We stored the integer as a pointer */
    lrank = counterNum % size;
    lidx = counterNum / size;
    printf("lrank %i\n",lrank);
    /* Update and return the counter */
    printf("cc keyval %i\n",MPE_COUNTER_KEYVAL);
    t_mpi[4] = MPI_Wtime();
    MPI_Win_lock(MPI_LOCK_SHARED, 0, lrank, counterWin); //shared
    t_mpi[5] = MPI_Wtime();
    MPI_Fetch_and_op(&one, value, MPI_INT,lrank, lidx, MPI_SUM, counterWin);
    t_mpi[2] = MPI_Wtime();
    MPI_Win_unlock(lrank, counterWin);
    t_mpi[3] = MPI_Wtime();
    return 0;
}
