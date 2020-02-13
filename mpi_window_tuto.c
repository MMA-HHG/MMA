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
    //  simulate the user inputs
    int[10] user_inputs; // = {15, 10, 5, 25, 10, 15, 50, 60, 0, 15};
    user_inputs[0] = 15; user_inputs[1] = 10; user_inputs[2] = 5; user_inputs[3] = 25; user_inputs[4] = 10; user_inputs[5] = 15;
    user_inputs[6] = 50; user_inputs[7] = 60; user_inputs[8] = 0; user_inputs[9] = 15;
    int iterator = 0;

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

MPI_Win_fence(0, nwin);
while (1) {
if (myid == 0) {
printf("Enter the number of intervals: (0 quits) ");
// fflush(stdout);
// scanf("%d",&n);
n = user_inputs[iterator];
pi = 0.0;
}
MPI_Win_fence(0, nwin);
if (myid != 0)
MPI_Get(&n, 1, MPI_INT, 0, 0, 1, MPI_INT, nwin);
MPI_Win_fence(0, nwin);
if (n == 0)
break;
else {
h = 1.0 / (double) n;
sum = 0.0;
for (i = myid + 1; i <= n; i += numprocs) {
x = h * ((double)i - 0.5);
sum += (4.0 / (1.0 + x*x));
}
mypi = h * sum;
MPI_Win_fence( 0, piwin);
MPI_Accumulate(&mypi, 1, MPI_DOUBLE, 0, 0, 1, MPI_DOUBLE,
MPI_SUM, piwin);
MPI_Win_fence(0, piwin);
if (myid == 0)
printf("pi is approximately %. 16f, Error is %. 16f/n",
pi, fabs(pi - PI25DT));
}
++iterator;
}
MPI_Win_free(&nwin);
MPI_Win_free(&piwin);
MPI_Finalize();
return 0;
}

