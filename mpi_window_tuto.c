/*
We test MPI-windows here

It implements the tuto for numerical computing pi. origin "Using MPI-2, chapter 2.3.2"
*/
#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "util.h"



/* Compute pi by numerical integration, RMA version */
int main(int argc, char *argv[])
{
int n, myid, numprocs, i;
double PI25DT = 3.141592653589793238462643;
double mypi, pi, h, sum, x;
MPI_Win nwin, piwin;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
MPI_Comm_rank(MPI_COMM_WORLD,&myid);
if (myid == 0) {
MPI_Win_create(&n, sizeof(int), 1, MPI_INFO_NULL,
MPI_COMM_WORLD, &nwin);
MPI_Win_create(&pi, sizeof(double), 1, MPI_INFO_NULL,
MPI_COMM_WORLD, &piwin);
}
else {
MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL,
MPI_COMM_WORLD, &nwin);
MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL,
MPI_COMM_WORLD, &piwin);
}
