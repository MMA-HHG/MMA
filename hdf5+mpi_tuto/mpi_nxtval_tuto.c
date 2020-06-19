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
MPE_Counter_nxtval( counter_win,0, &counter_value, MPE_COUNTER_KEYVAL );
printf("I am node %d of %d and my counter value is %d \n", myid, numprocs, counter_value);
MPE_Counter_nxtval( counter_win,0, &counter_value, MPE_COUNTER_KEYVAL );
printf("I am node %d of %d and my counter value is %d \n", myid, numprocs, counter_value);
MPE_Counter_nxtval( counter_win,0, &counter_value, MPE_COUNTER_KEYVAL );
printf("I am node %d of %d and my counter value is %d \n", myid, numprocs, counter_value);

// window should be deleted

MPI_Finalize();
return 0;
}
