/*
We test MPI3-mutex here
*/
#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
// #include "util.h"


int MPE_MUTEX_KEYVAL;

/* see desrption pg 198 MPI-2 */

// github.com/shawfdong/ams250/blob/master/examples/mpi/Using_MPI/advanced/nxtval-test.c

int main(int argc, char *argv[])
{

int myid, numprocs, i;


// int MPE_COUNTER_KEYVAL; // how to fix this?

MPI_Win mutex_win; // this is memory window for the counter
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);

if (myid == 0) { printf("Non-mutexed outputs are here\n"); fflush(stdout);}
MPI_Barrier( MPI_COMM_WORLD );

printf("1 I am node %d of %d\n", myid, numprocs); fflush(stdout);
sleep(2);
printf("2 I am node %d of %d\n", myid, numprocs); fflush(stdout);
sleep(2);
printf("3 I am node %d of %d\n", myid, numprocs); fflush(stdout);


MPI_Barrier( MPI_COMM_WORLD );
if (myid == 0) {printf("Mutexed outputs are here\n"); fflush(stdout);}
MPI_Barrier( MPI_COMM_WORLD );

MPE_MUTEX_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD,1, &mutex_win );

MPE_Mutex_acquire(mutex_win, 0, MPE_MUTEX_KEYVAL);

printf("1 I am node %d of %d\n", myid, numprocs); fflush(stdout);
sleep(2);
printf("2 I am node %d of %d\n", myid, numprocs); fflush(stdout);
sleep(2);
printf("3 I am node %d of %d\n", myid, numprocs); fflush(stdout);

MPE_Mutex_release(mutex_win, 0, MPE_MUTEX_KEYVAL);


// window should be deleted

MPI_Finalize();
return 0;
}
