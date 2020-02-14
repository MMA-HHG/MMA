/*
We test MPI-windows here
*/
#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
// #include "util.h"


int MPE_COUNTER_KEYVAL;

/* see desrption pg 198 MPI-2 */

// github.com/shawfdong/ams250/blob/master/examples/mpi/Using_MPI/advanced/nxtval-test.c

int main(int argc, char *argv[])
{

int myid, numprocs, i;

int counter_value = 0;

// int MPE_COUNTER_KEYVAL; // how to fix this?

MPI_Win counter_win; // this is memory window for the counter
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);

// MPE_setKeyval(0); // attempt

// create counter
MPE_COUNTER_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD,1, &counter_win );

printf("I am node %d of %d and my counter value before addition is %d \n", myid, numprocs, counter_value);

MPI_Win_fence(0, counter_win);

printf("I am node %d of %d and my key %d \n", myid, numprocs, MPE_COUNTER_KEYVAL);

MPE_Counter_nxtval( counter_win,1, &counter_value, MPE_COUNTER_KEYVAL );

MPI_Win_fence(0, counter_win);
printf("I am node %d of %d and my counter value is %d \n", myid, numprocs, counter_value);
printf("I am node %d of %d and my key %d \n", myid, numprocs, MPE_COUNTER_KEYVAL);

printf("fence \n");
fflush(stdout);

MPE_Counter_nxtval( counter_win,1, &counter_value, MPE_COUNTER_KEYVAL );

MPI_Win_fence(0, counter_win);
printf("I am node %d of %d and my counter value is %d \n", myid, numprocs, counter_value);
printf("fence \n");
fflush(stdout);

MPE_Counter_nxtval( counter_win,1, &counter_value, MPE_COUNTER_KEYVAL );

printf("I am node %d of %d and my counter value is %d \n", myid, numprocs, counter_value);



// if (myid == 0){
// printf("proc 0: Printing mine counter array \n");
// for (i = 0; i < numprocs; i++)
// {
//     printf("the value %d \n",counter[i]);
// }
// }

MPI_Finalize();
return 0;
}

// done in the code
// int counter[numprocs];
// for (i = 0; i < numprocs; i++) { counter[i] = 0; } // set to 0;

// if (myid == 0) {
// // the counter memory sits in the first process's memory
// // tuto says its size dimension should contain all the ranks
// MPI_Win_create(&counter, numprocs*sizeof(int), 1, MPI_INFO_NULL, MPI_COMM_WORLD, &counter_win);
// }
// else {
// MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &counter_win);
// }
// MPI_Win_fence(0, counter_win); // we have the counter window set