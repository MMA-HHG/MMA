PROGRAM cuprad
  USE fields
  USE parameters
  USE first_step
  USE output
  USE outputHDF5
  USE long_step
  USE mpi_stuff
  USE run_status
  USE longstep_vars

  IMPLICIT NONE

!  INTEGER(4) time, starttime, deltatime, limit_s
  INTEGER(4) time, starttime, deltatime, limit_s
!  external time
  INTEGER(4) :: tcount, count_rate, count_max
  !  REAL(8)  TREM

!  starttime = time()
  call system_clock(tcount, count_rate, count_max)
  write(*,*) tcount, count_rate, count_max
  CALL initialize
  limit_s=timelimit*3600
  IF (.NOT.finished) THEN
     DO WHILE (z.LT.proplength)

        IF(z_out.LE.z) THEN ! HDF5 printing
           CALL HDF5_out! that's the printing
!           z_outHD5=z_outHD5+outlengthHD5 !!!! WILL BE USED WHEN GRIDS DISATTACHED
        ENDIF

        IF(z_out.LE.z) THEN
           !CALL matlab_out ! that's the printing
           CALL write_output
           z_out=z_out+outlength
        ENDIF


        CALL propagation

        IF (maxphase.GT.decrease) THEN
           delta_z=0.5D0*(decrease+increase)/maxphase*delta_z
           IF(my_rank.EQ.0) THEN
              OPEN(unit_rho,FILE='ZSTEP.DAT',STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(unit_rho,*) 'z=',REAL(z,4),' delta_z=',REAL(delta_z,4),REAL(maxphase,4)
              CLOSE(unit_rho)
           ENDIF
           call calc_propagator
        ENDIF
        IF ((maxphase.LT.increase).AND.(delta_z.LT.delta_z_max)) THEN
           delta_z=MIN(delta_z_max,0.5D0*(decrease+increase)/maxphase*delta_z)
           IF(my_rank.EQ.0) THEN
              OPEN(unit_rho,FILE='ZSTEP.DAT',STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(unit_rho,*) 'z=',REAL(z,4),' delta_z=',REAL(delta_z,4),REAL(maxphase,4)
              CLOSE(unit_rho)
           ENDIF
           call calc_propagator
        ENDIF
!        deltatime = time() - starttime
        CALL MPI_BCAST(deltatime,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        IF (deltatime.GE.limit_s) EXIT
        !     IF (my_rank.EQ.0) CALL TREMAIN(TREM)
        !     CALL MPI_BCAST(TREM,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !     IF (TREM.LE.3600) EXIT     
     ENDDO
     IF (z.GE.proplength) THEN
        CALL HDF5_out
        finished = .TRUE.
     ELSE
        IF (my_rank.EQ.0) THEN
           OPEN(unit_rho,FILE='STOP',STATUS='OLD')
           CLOSE(unit_rho,STATUS='DELETE')
        ENDIF
     ENDIF
     CALL write_output 
     !CALL matlab_out 
     rhodist=count
     CALL propagation
     CALL field_out
     CALL linked_list_out
     CALL finalize
  ENDIF

END PROGRAM cuprad
