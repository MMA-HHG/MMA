/*
We test MPI3-mutex together with a counter, both placed in a shared window, we also test there is no interference

this is the final tutorial that will be used in the TDSE procedure:

mutexes will be there for writing, the counter will be used for assigning simulations to workers at the moment they finish their work.

*/
#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
// #include "util.h"


int MPE_MC_KEYVAL;


int main(int argc, char *argv[])
{

int myid, numprocs, i;
int counter_value;


MPI_Win mc_win; // this is the shared window

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);

MPE_MC_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD,2, &mc_win );

if (myid == 0) { printf("Non-mutexed outputs are here\n"); fflush(stdout);}
MPI_Barrier( MPI_COMM_WORLD );

    MPE_Counter_nxtval( mc_win,0, &counter_value, MPE_MC_KEYVAL );
    printf("1 I am node %d of %d, counter value is %d \n", myid, numprocs, counter_value); fflush(stdout);
    sleep(2);

    MPE_Counter_nxtval( mc_win,0, &counter_value, MPE_MC_KEYVAL );
    printf("2 I am node %d of %d, counter value is %d \n", myid, numprocs, counter_value); fflush(stdout);
    sleep(2);

    MPE_Counter_nxtval( mc_win,0, &counter_value, MPE_MC_KEYVAL );
    printf("3 I am node %d of %d, counter value is %d \n", myid, numprocs, counter_value); fflush(stdout);


MPI_Barrier( MPI_COMM_WORLD ); 
if (myid == 0) {printf("Mutexed outputs are here\n"); fflush(stdout);}
MPI_Barrier( MPI_COMM_WORLD );



MPE_Mutex_acquire(mc_win, 1, MPE_MC_KEYVAL); // mutex sits on the second element

    MPE_Counter_nxtval( mc_win,0, &counter_value, MPE_MC_KEYVAL );
    printf("1 I am node %d of %d, counter value is %d \n", myid, numprocs, counter_value); fflush(stdout);
    sleep(2);

    MPE_Counter_nxtval( mc_win,0, &counter_value, MPE_MC_KEYVAL );
    printf("2 I am node %d of %d, counter value is %d \n", myid, numprocs, counter_value); fflush(stdout);
    sleep(2);

    MPE_Counter_nxtval( mc_win,0, &counter_value, MPE_MC_KEYVAL );
    printf("3 I am node %d of %d, counter value is %d \n", myid, numprocs, counter_value); fflush(stdout);
    sleep(2);

MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);


// window should be deleted

MPI_Finalize();
return 0;
}
