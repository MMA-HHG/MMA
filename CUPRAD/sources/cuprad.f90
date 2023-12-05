!=============================================================!
! Here is the main cuprad program, using the propagation loop.
! It prepares the calculation and then loop until finshed or 
! estimated time limit is reached.                              !!! NOT REIMPLEMENTED NOW
!
! 
! The main developer of the code is Stefan Skupin. The code uses
! many modules co-developed by others contributors, some of
! them are mentioned in the respective modules.
!=============================================================!

PROGRAM cuprad
  USE fields
  USE parameters
  USE first_step
  USE output
  USE long_step
  USE mpi_stuff
  USE run_status
  USE longstep_vars
  USE density_module

  USE normalization

  IMPLICIT NONE 
  

  INTEGER(4) deltatime, limit_s
!  external time
  INTEGER(4) :: tcount, count_rate, count_max
  !  REAL(8)  TREM

  REAL(8) :: local_time_MPI

!=====================
! INITIALISATION PHASE
!=====================

!  starttime = time()
  call system_clock(tcount, count_rate, count_max)
  write(*,*) tcount, count_rate, count_max
  CALL initialize
  limit_s=timelimit*3600

  start_time_MPI = MPI_Wtime()

!==========================
! MAIN COMPUTATIONAL PHASE
!==========================

  IF (.NOT.finished) THEN ! it's here in the case one runs cuprads consquitively (e.g.) changing medium. TO BE REINTRODUCED
 
     ! it loops until the end of the medium is reached or timelimit is reached (exit statement in the body)
     DO WHILE (z.LT.proplength) 

         ! there are two independent writes (see manual) 
         IF( out_Efield .AND. (z_out_Efield .LE. z)) THEN ! only-field print ! ADD an extra logic to disable this option at all
           CALL Efield_out
           z_out_Efield = z_out_Efield + outlength_Efield
         ENDIF

         IF(z_out.LE.z) THEN
           local_time_MPI  = MPI_Wtime()
           IF (my_rank.EQ.0) THEN
            print *, '-------------------------------------------------------------------------------'
            print *, "printing number:", output_write_count, ":"
            print *, "z[m]=", z*four_z_Rayleigh
            print *, "before printing:", local_time_MPI - start_time_MPI
           ENDIF
           CALL write_output
           local_time_MPI  = MPI_Wtime()
           IF (my_rank.EQ.0) THEN
            print *, "after printing:", local_time_MPI - start_time_MPI
           ENDIF           
           z_out = z_out + outlength
         ENDIF

        ! Standard operation of the code: call propagation until end of medium reached
        CALL propagation
        
       IF (apply_density_mod) THEN
         CALL calc_density_mod(z)
         IF (is_density_changed) THEN
            call calc_time_propagator
            call calc_cn_propagator
         ENDIF
       ENDIF

        ! Adaptive step-size controlling + the duration of the calculation
        ! Step size is changed, propagation operator recalculated and step-size stored in outputs.

        IF (maxphase.GT.decrease) THEN ! derease step size
           delta_z=0.5D0*(decrease+increase)/maxphase*delta_z
           IF(my_rank.EQ.0) THEN
              CALL write_extended_dz    ! save the length of the steps
           ENDIF
           call calc_time_propagator
           call calc_cn_propagator
        ENDIF
        IF ((maxphase.LT.increase).AND.(delta_z.LT.delta_z_max)) THEN ! increase step size (not above maximally allowed)
           delta_z=MIN(delta_z_max,0.5D0*(decrease+increase)/maxphase*delta_z)
           IF(my_rank.EQ.0) THEN
              CALL write_extended_dz
           ENDIF
           call calc_time_propagator
           call calc_cn_propagator
        ENDIF

        ! Check if time limit for the calculation reached
        !deltatime = time() - starttime
        CALL MPI_BCAST(deltatime,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        IF (deltatime.GE.limit_s) EXIT
        !     IF (my_rank.EQ.0) CALL TREMAIN(TREM)
        !     CALL MPI_BCAST(TREM,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !     IF (TREM.LE.3600) EXIT     
     ENDDO


    !====================
    ! FINALISATION PHASE
    !====================

     IF (z.GE.proplength) THEN ! The end of the medium is reached
        IF (out_Efield) CALL Efield_out         
        finished = .TRUE. ! prevents possible consequitive slurm job to be executed
        IF (my_rank.EQ.0) THEN
           OPEN(unit_rho,FILE='STOP',STATUS='OLD')
           CLOSE(unit_rho,STATUS='DELETE')
        ENDIF
     ENDIF
     
     ! end-plane outputs are written in all cases at the end of the code
     CALL write_output 
     !CALL matlab_out 
     rhodist=count
     CALL propagation ! ??? WHY DO WE PROPAGATE EVEN ONCE MORE?
     CALL field_out
     CALL linked_list_out
     IF (my_rank.EQ.0) THEN ! write output to check if the code is not completely wrong
         print *, 'test:'
         print *, '(4.05444499579252,2.107331137558244E-002)'
         print *, e(dim_t/2,2)
     ENDIF
     CALL finalize
     PRINT*, "program finished"
  ENDIF

END PROGRAM cuprad
