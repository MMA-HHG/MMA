/*
this program writes some outputs with waiting
*/
#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
/*#include "hdf5.h"*/
#include "util.h"



/* Compute pi by numerical integration, RMA version */
int main(int argc, char *argv[])
{

int n, myid, numprocs, i;
MPI_Win nwin, piwin;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);


		printf("1 I am node %d of %d\n", myid, numprocs);
        sleep(2);
		printf("2 I am node %d of %d\n", myid, numprocs);
        sleep(2);
		printf("3 I am node %d of %d\n", myid, numprocs);


MPI_Finalize();
return 0;
}

