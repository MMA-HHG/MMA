#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>


int k1, k2, k3;


clock_t start_main, finish2_main, finish1_main, finish3_main, finish4_main;


int main(int argc, char *argv[]) 
{
	// vars:
	// dummy
	int dum3int[3];

	// Processing the queue


	int comment_operation = 1;
	double t_mpi[10]; 

	// Initialise MPI
	int myrank, nprocs;
	int MPE_MC_KEYVAL, MPE_C_KEYVAL, MPE_M_KEYVAL; // this is used to address the mutex and counter
	MPI_Win mc_win, c_win, m_win; // this is the shared window, it is used both  for mutices and counter
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	if (comment_operation == 1 ){printf("Proc %i started the program\n",myrank);}

	////////////////////////
	// PREPARATION PAHASE //
	////////////////////////
	t_mpi[0] = MPI_Wtime();	start_main = clock(); // the clock	

    int Ntot = 6;
    int t_job = 15; // [s]
    int Nsim; 

	

	//////////////////////////
	// COMPUTATIONAL PAHASE //
	//////////////////////////

	// create counter and mutex in one pointer
	//MPE_MC_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 2, &mc_win); // first is counter, second mutex
	MPE_M_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 1, &m_win); // first is counter, second mutex
	MPE_C_KEYVAL = MPE_Counter_create(MPI_COMM_WORLD, 1, &c_win); // first is counter, second mutex


	// first process is preparing the file and the rest may do their own work (the file is locked in the case they want to write); it creates the resulting dataset
	// an empty dataset is prepared to be filled with the data
	
	if ( myrank == 0 ){
		MPE_Mutex_acquire(m_win, 0, MPE_M_KEYVAL); // first process get mutex and hold it to ensure to prepare the file before the others
		MPE_Counter_nxtval(c_win, 0, &Nsim, MPE_C_KEYVAL); // get my first task
	}
	MPI_Barrier(MPI_COMM_WORLD); // Barrier
	printf("Proc %i abarier \n",myrank);
	if ( myrank == 0 ){
		sleep(t_job);

		MPE_Mutex_release(m_win, 0, MPE_M_KEYVAL);
	}

	t_mpi[6] = MPI_Wtime();
	printf("Proc %i, reached the point 1  : %f sec\n",myrank,t_mpi[6]-t_mpi[0]);
	
 
	t_mpi[4] = MPI_Wtime();	
		

	// we now process the MPI queue
	MPE_Counter_nxtval(c_win, 0, &Nsim, MPE_C_KEYVAL); // get my first task (every worker calls, second call for proc 0 due to the preparation phase)

	t_mpi[7] = MPI_Wtime();
	printf("Proc %i, reached the point 2  : %f sec\n",myrank,t_mpi[7]-t_mpi[0]);

	t_mpi[5] = MPI_Wtime(); 

    
	while (Nsim < Ntot){ // run till queue is not treated

		// MPE_Mutex_release(mc_win, 1, MPE_MC_KEYVAL);

    	t_mpi[3] = MPI_Wtime(); finish3_main = clock();
		sleep(t_job); 
		t_mpi[1] = MPI_Wtime(); finish1_main = clock();

		// write results
		MPE_Mutex_acquire(m_win, 0, MPE_M_KEYVAL); // mutex is acquired

		t_mpi[2] = MPI_Wtime(); finish2_main = clock();
		if ( ( comment_operation == 1 ) && ( Nsim < 20 ) ){
			printf("Proc %i job %i\n",myrank,Nsim);
			printf("Proc %i, returned mutex last time   : %f sec\n",myrank,t_mpi[4]-t_mpi[0]);
			printf("Proc %i, obtained counter last time : %f sec\n",myrank,t_mpi[5]-t_mpi[0]);
			printf("Proc %i, before job started         : %f sec\n",myrank,t_mpi[3]-t_mpi[0]);
			printf("Proc %i, clock the umnutexed value  : %f sec\n",myrank,t_mpi[1]-t_mpi[0]);
			printf("Proc %i, clock in the mutex block   : %f sec\n",myrank,t_mpi[2]-t_mpi[0]);
			fflush(NULL); // force write
    	}
    


		MPE_Mutex_release(m_win, 0, MPE_M_KEYVAL);
    	t_mpi[4] = MPI_Wtime();
		
		// outputs_destructor(outputs); // free memory
		MPE_Counter_nxtval(c_win, 0, &Nsim, MPE_C_KEYVAL); // get my next task
		t_mpi[5] = MPI_Wtime();
	}


	
 
//  if ( myrank == 0 )  // time
//	{
//  finish2 = clock();
//  printf("\nFirst processor measured time: %f sec\n\n",(double)(finish2 - start2) / CLOCKS_PER_SEC);
//  }
  // MPI_Win_free(&mc_win);
	MPI_Finalize();
	return 0;	
}

