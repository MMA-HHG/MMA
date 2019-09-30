SUBROUTINE finalize
  USE fields
  USE parameters
  USE fft
  IMPLICIT NONE

  CALL fft_exit
  CLOSE(unit_peakmax)
  CLOSE(unit_energy)
  CLOSE(unit_rhomax)
  CLOSE(unit_rho)
  CLOSE(unit_logfile)
  call MPI_Finalize(ierr)

  RETURN
END SUBROUTINE finalize
