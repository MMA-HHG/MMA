#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<unistd.h>
#include"util_mpi.h"


int k1, k2, k3;
double t_mpi[10]; 


clock_t start_main, finish2_main, finish1_main, finish3_main, finish4_main;


int main(int argc, char *argv[]) 
{
	// vars:
	// dummy
	int dum3int[3];

	// Processing the queue


	int comment_operation = 1, disp_tasks = 50;
	

	// Initialise MPI
	int myrank, nprocs;
	MPI_Win mc_win, c_win, m_win; // this is the shared window, it is used both  for mutices and counter
	MPI_Init(&argc,&argv);
	t_mpi[0] = MPI_Wtime(); // the clock
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	if (comment_operation == 1 ){printf("Proc %i started the program\n",myrank);}

	////////////////////////
	// PREPARATION PAHASE //
	////////////////////////
		

    int Ntot = 6;
    int t_job = 10; // [s]
    int Nsim; 
    int local_counter = 0;
    printf("invalid keyval %i\n",MPI_KEYVAL_INVALID);

	

	//////////////////////////
	// COMPUTATIONAL PAHASE //
	//////////////////////////

    MPE_Counter_create(MPI_COMM_WORLD, 2, &c_win);




    MPE_Counter_nxtval(c_win, 0, &Nsim); 
    sleep(t_job);
    t_mpi[1] = MPI_Wtime();
    printf("Proc %i, point 1, value %i\narrived:           %f sec\nattempted to lock: %f sec\nhaving lock:       %f sec\nincreased value:   %f sec\nunlocked window:   %f sec\n\n",
	myrank,Nsim,t_mpi[1]-t_mpi[0],t_mpi[4]-t_mpi[0],t_mpi[5]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);
    MPE_Counter_nxtval(c_win, 0, &Nsim); 
    sleep(t_job);
    t_mpi[1] = MPI_Wtime();
    //printf("Proc %i, 2  : %f sec\n",myrank,t_mpi[1]-t_mpi[0]);
    //printf("Proc %i, point 2\n  arrived: %f sec\nincreased value %f sec\nunlocked window: %f sec\n\n",myrank,t_mpi[1]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);
    printf("Proc %i, point 2, value %i\narrived:           %f sec\nattempted to lock: %f sec\nhaving lock:       %f sec\nincreased value:   %f sec\nunlocked window:   %f sec\n\n",
	myrank,Nsim,t_mpi[1]-t_mpi[0],t_mpi[4]-t_mpi[0],t_mpi[5]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);

    MPE_Counter_nxtval(c_win, 0, &Nsim); 
    sleep(t_job);
    t_mpi[1] = MPI_Wtime();
    //printf("Proc %i, 3  : %f sec\n",myrank,t_mpi[1]-t_mpi[0]);
    //printf("Proc %i, point 3\n  arrived: %f sec\nincreased value %f sec\nunlocked window: %f sec\n\n",myrank,t_mpi[1]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);
    printf("Proc %i, point 3, value %i\narrived:           %f sec\nattempted to lock: %f sec\nhaving lock:       %f sec\nincreased value:   %f sec\nunlocked window:   %f sec\n\n",
	myrank,Nsim,t_mpi[1]-t_mpi[0],t_mpi[4]-t_mpi[0],t_mpi[5]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);

	
 


    MPI_Barrier(MPI_COMM_WORLD); // Barrier
    MPI_Finalize();
    return 0;	
}

