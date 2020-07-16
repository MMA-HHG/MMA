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
	int MPE_MC_KEYVAL, MPE_C_KEYVAL, MPE_M_KEYVAL; // this is used to address the mutex and counter
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

	

	//////////////////////////
	// COMPUTATIONAL PAHASE //
	//////////////////////////

	// create counter and mutex in one pointer
	//MPE_MC_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 2, &mc_win); // first is counter, second mutex
	//MPE_M_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 1, &m_win); // first is counter, second mutex
	MPE_C_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 1, &c_win); // first is counter, second mutex


	// first process is preparing the file and the rest may do their own work (the file is locked in the case they want to write); it creates the resulting dataset
	// an empty dataset is prepared to be filled with the data
	
	// if ( myrank == 0 ){
	// 	//MPE_Mutex_acquire(m_win, 0, MPE_M_KEYVAL); // first process get mutex and hold it to ensure to prepare the file before the others
	// 	MPE_Counter_nxtval(c_win, 0, &Nsim, MPE_C_KEYVAL); // get my first task
	// }
	// MPI_Barrier(MPI_COMM_WORLD); // Barrier
	// printf("Proc %i abarier \n",myrank);
	// if ( myrank == 0 ){
	// 	sleep(t_job);
	// 	local_counter++;

	// 	//MPE_Mutex_release(m_win, 0, MPE_M_KEYVAL);
	// }

    MPE_Counter_nxtval(c_win, 0, &Nsim, MPE_C_KEYVAL); 
    sleep(t_job);
    t_mpi[1] = MPI_Wtime();
    printf("Proc %i, point 1, value %i\narrived:           %f sec\nattempted to lock: %f sec\nhaving lock:       %f sec\nincreased value:   %f sec\nunlocked window:   %f sec\n\n",
	myrank,Nsim,t_mpi[1]-t_mpi[0],t_mpi[4]-t_mpi[0],t_mpi[5]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);
    MPE_Counter_nxtval(c_win, 0, &Nsim, MPE_C_KEYVAL); 
    sleep(t_job);
    t_mpi[1] = MPI_Wtime();
    //printf("Proc %i, 2  : %f sec\n",myrank,t_mpi[1]-t_mpi[0]);
    //printf("Proc %i, point 2\n  arrived: %f sec\nincreased value %f sec\nunlocked window: %f sec\n\n",myrank,t_mpi[1]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);
    printf("Proc %i, point 2, value %i\narrived:           %f sec\nattempted to lock: %f sec\nhaving lock:       %f sec\nincreased value:   %f sec\nunlocked window:   %f sec\n\n",
	myrank,Nsim,t_mpi[1]-t_mpi[0],t_mpi[4]-t_mpi[0],t_mpi[5]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);

    MPE_Counter_nxtval(c_win, 0, &Nsim, MPE_C_KEYVAL); 
    sleep(t_job);
    t_mpi[1] = MPI_Wtime();
    //printf("Proc %i, 3  : %f sec\n",myrank,t_mpi[1]-t_mpi[0]);
    //printf("Proc %i, point 3\n  arrived: %f sec\nincreased value %f sec\nunlocked window: %f sec\n\n",myrank,t_mpi[1]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);
    printf("Proc %i, point 3, value %i\narrived:           %f sec\nattempted to lock: %f sec\nhaving lock:       %f sec\nincreased value:   %f sec\nunlocked window:   %f sec\n\n",
	myrank,Nsim,t_mpi[1]-t_mpi[0],t_mpi[4]-t_mpi[0],t_mpi[5]-t_mpi[0],t_mpi[2]-t_mpi[0],t_mpi[3]-t_mpi[0]);

	// t_mpi[6] = MPI_Wtime();
    // sleep(t_job);
	// printf("Proc %i, reached the point 1  : %f sec\n",myrank,t_mpi[6]-t_mpi[0]);
	
 


	MPI_Barrier(MPI_COMM_WORLD); // Barrier
	MPI_Finalize();
	return 0;	
}

