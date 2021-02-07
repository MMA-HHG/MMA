! This is the essential module to compute FFT for a matrix distributed among the workers.
! The send-receive scheme used here was designed by Stefan Skupin

MODULE fft
  USE fields
  USE parameters
  USE mpi_stuff

  INTEGER(8) plan_forward1,plan_forward2,plan_backward1,plan_backward2,plan_spec,plan_p,plan_j,plan_pharm
  REAL(8)  diminv

CONTAINS

  SUBROUTINE fft_init
    IMPLICIT NONE

    INTEGER(4) schema(num_proc),i,j,power,help,m,n,s
    INTEGER(4), ALLOCATABLE  ::  array_of_blocklengths(:),array_of_displacements(:)

    DO j=1,num_proc
       schema(j)=j
    ENDDO
    power=0
    DO j=1,16
       IF (2**j.EQ.num_proc) power=j
    ENDDO
    IF (power.EQ.0) WRITE(6,*) 'ERROR in number of procs'
    IF ((dim_t/num_proc)*num_proc.NE.dim_t) WRITE(6,*) 'ERROR in dim_t'
    IF ((dim_r/num_proc)*num_proc.NE.dim_r) WRITE(6,*) 'ERROR in dim_r'
    help=my_rank
    DO j=power-1,0,-1
       IF (help.GE.2**j) THEN
          help=help-2**j
          CALL vertausch(schema,num_proc,2**j )
       ENDIF
    ENDDO
    ALLOCATE (dim_t_start(num_proc),dim_t_end(num_proc),dim_r_start(num_proc),dim_r_end(num_proc))
    ALLOCATE (num_ex(num_proc-1),send_first(num_proc-1))
    DO j=1,num_proc
       IF  (schema(j).GT.1) THEN
          dim_t_start(schema(j)-1)=dim_t/num_proc*(j-1)+1
          dim_t_end(schema(j)-1)=dim_t/num_proc*j
          dim_r_start(schema(j)-1)=dim_r/num_proc*(j-1)+1
          dim_r_end(schema(j)-1)=dim_r/num_proc*j
          num_ex(schema(j)-1)=j-1
          send_first(schema(j)-1)=(num_ex(schema(j)-1).GT.my_rank)
       ENDIF
    ENDDO
    dim_t_start(num_proc)=dim_t/num_proc*my_rank+1
    dim_t_end(num_proc)=dim_t/num_proc*(my_rank+1)
    dim_r_start(num_proc)=dim_r/num_proc*my_rank+1
    dim_r_end(num_proc)=dim_r/num_proc*(my_rank+1)
    dim_t_local=dim_t/num_proc
    dim_r_local=dim_r/num_proc
    diminv=1.D0/REAL(dim_t,8)

    ALLOCATE(e(dim_t,dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(efft(dim_r,dim_t_start(num_proc):dim_t_end(num_proc)),etemp(dim_t,dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(ptemp(dim_t,dim_r_start(num_proc):dim_r_end(num_proc)),jtemp(dim_t,dim_r_start(num_proc):dim_r_end(num_proc)))
    ALLOCATE(pharm(dim_t,dim_r_start(num_proc):dim_r_end(num_proc)))

    CALL MPI_TYPE_vector(dim_r/num_proc,2*dim_t/num_proc,2*dim_t,MPI_DOUBLE_PRECISION,MPI_SUBARRAY,ierr)
    CALL MPI_TYPE_COMMIT(MPI_SUBARRAY,ierr)

    ALLOCATE(array_of_blocklengths(dim_t_local*dim_r_local),array_of_displacements(dim_t_local*dim_r_local))
    array_of_blocklengths=2
    DO j=0,dim_r_local-1
       DO i=0,dim_t_local-1
          array_of_displacements(dim_t_local*j+i+1)=2*j+2*dim_r*i
       ENDDO
    ENDDO
    CALL MPI_TYPE_indexed(dim_t_local*dim_r_local,array_of_blocklengths,array_of_displacements, &
         MPI_DOUBLE_PRECISION,MPI_SUBARRAY_TRANSPOSED,ierr)
    DEALLOCATE(array_of_blocklengths,array_of_displacements)
    CALL MPI_TYPE_COMMIT(MPI_SUBARRAY_TRANSPOSED,ierr)

    ALLOCATE (send_e(num_proc-1),send_etemp(num_proc-1),send_efft(num_proc-1))
    ALLOCATE (recv_e(num_proc-1),recv_etemp(num_proc-1),recv_efft(num_proc-1))
    DO j=1,num_proc-1
       CALL MPI_SSEND_INIT(e(dim_t_start(j),dim_r_start(num_proc)),1,MPI_SUBARRAY,num_ex(j),1,MPI_COMM_WORLD,send_e(j),ierr)
       CALL MPI_SSEND_INIT(efft(dim_r_start(j),dim_t_start(num_proc)),1,MPI_SUBARRAY_TRANSPOSED,num_ex(j),1,MPI_COMM_WORLD,send_efft(j),ierr)
       CALL MPI_RECV_INIT(e(dim_t_start(j),dim_r_start(num_proc)),1,MPI_SUBARRAY,num_ex(j),1,MPI_COMM_WORLD,recv_e(j),ierr)
       CALL MPI_RECV_INIT(efft(dim_r_start(j),dim_t_start(num_proc)),1,MPI_SUBARRAY_TRANSPOSED,num_ex(j),1,MPI_COMM_WORLD,recv_efft(j),ierr)
    ENDDO

    m=dim_r_local
    n=dim_t
    s=1
    CALL dfftw_plan_many_dft(plan_forward1,1,dim_t,m,e,dim_t,s,n,e,dim_t,s,n,1,0)
    CALL dfftw_plan_many_dft(plan_backward2,1,dim_t,m,e,dim_t,s,n,e,dim_t,s,n,-1,0)
    CALL dfftw_plan_many_dft(plan_p,1,dim_t,m,ptemp,dim_t,s,n,ptemp,dim_t,s,n,1,0)
    CALL dfftw_plan_many_dft(plan_j,1,dim_t,m,jtemp,dim_t,s,n,jtemp,dim_t,s,n,1,0)
    CALL dfftw_plan_many_dft(plan_pharm,1,dim_t,m,pharm,dim_t,s,n,pharm,dim_t,s,n,1,0)
    m=dim_r_local
    n=dim_t
    s=1
    CALL dfftw_plan_many_dft(plan_spec,1,dim_t,m,etemp,dim_t,s,n,etemp,dim_t,s,n,1,0)
    RETURN
  END SUBROUTINE fft_init

  SUBROUTINE vertausch(vector,n,m)
    IMPLICIT NONE

    INTEGER(4)  n,m,vector(n),i,j,help

    DO i=0,n/(2*m)-1
       DO j=1,m
          help=vector(2*m*i+j)
          vector(2*m*i+j)=vector(2*m*i+j+m)
          vector(2*m*i+j+m)=help
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE vertausch

  SUBROUTINE fft_forward_inplace
    IMPLICIT NONE

    INTEGER(4) k,l,i

    CALL dfftw_execute(plan_forward1)
    SELECT CASE (switch_T)
    CASE(1)
       continue
    CASE(2)
       CALL dfftw_execute(plan_p)
       CALL dfftw_execute(plan_j)
       DO k=dim_r_start(num_proc),dim_r_end(num_proc)
          DO l=1,dim_t
             e(l,k)=e(l,k)+(op_t(l)-1.D0)*ptemp(l,k)+(op_t_inv(l)-1.D0)*jtemp(l,k)
          ENDDO
       ENDDO
    CASE(3)
       CALL dfftw_execute(plan_p)
       CALL dfftw_execute(plan_j)
       CALL dfftw_execute(plan_pharm)
       DO k=dim_r_start(num_proc),dim_r_end(num_proc)
          DO l=dim_th+1,dim_t
             e(l,k)=e(l,k)+(op_t(l)-1.D0)*ptemp(l,k)+(op_t_inv(l)-1.D0)*jtemp(l,k)+op_t(l)*pharm(l,k)
          ENDDO
          DO l=dim_th+1,dim_t-2*omega_offset(1)+1
             e(l,k)=e(l,k)-op_t(l)*CONJG(ptemp(dim_t-2*(omega_offset(1)-1)-l,k))
          ENDDO
          DO l=dim_t-2*(omega_offset(1)-1),dim_t
             e(l,k)=e(l,k)-op_t(l)*CONJG(ptemp(2*dim_t-2*(omega_offset(1)-1)-l,k))
          ENDDO
!          DO l=1,dim_th
!             e(l,k)=e(l,k)-phaserelation(3)*op_t(l)*CONJG(ptemp(dim_t-2*(omega_offset(1)-1)-l,k)*phaserelation(1))
!          ENDDO
       ENDDO
    CASE(4)
       CALL dfftw_execute(plan_p)
       CALL dfftw_execute(plan_j)
       DO k=dim_r_start(num_proc),dim_r_end(num_proc)
          DO l=1,dim_t
             e(l,k)=e(l,k)+(op_t(l)-1.D0)*ptemp(l,k)+(op_t_inv(l)-1.D0)*jtemp(l,k)
          ENDDO
       ENDDO
    END SELECT
    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
       DO k=dim_t_start(num_proc),dim_t_end(num_proc)
          efft(l,k)=e(k,l)
       ENDDO
    ENDDO
    DO i=1,num_proc-1
       IF (send_first(i)) THEN
          CALL MPI_START(send_e(i),ierr)
          CALL MPI_WAIT(send_e(i),status,ierr)
          CALL MPI_START(recv_efft(i),ierr)
          CALL MPI_WAIT(recv_efft(i),status,ierr)
       ELSE
          CALL MPI_START(recv_efft(i) ,ierr)
          CALL MPI_WAIT(recv_efft(i),status,ierr)
          CALL MPI_START(send_e(i),ierr)
          CALL MPI_WAIT(send_e(i),status,ierr)
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE fft_forward_inplace

  SUBROUTINE fft_backward2_inplace
    IMPLICIT NONE

    INTEGER(4)  k,l,i

    DO l=dim_r_start(num_proc),dim_r_end(num_proc)
       DO k=dim_t_start(num_proc),dim_t_end(num_proc)
          e(k,l)=efft(l,k)
       ENDDO
    ENDDO
    DO i=1,num_proc-1
       IF (send_first(i)) THEN
          CALL MPI_START(send_efft(i),ierr)
          CALL MPI_WAIT(send_efft(i),status,ierr)
          CALL MPI_START(recv_e(i),ierr)
          CALL MPI_WAIT(recv_e(i),status,ierr)
       ELSE
          CALL MPI_START(recv_e(i),ierr)
          CALL MPI_WAIT(recv_e(i),status,ierr)
          CALL MPI_START(send_efft(i),ierr)
          CALL MPI_WAIT(send_efft(i),status,ierr)
       ENDIF
    ENDDO
    CALL dfftw_execute(plan_backward2)
    e=diminv*e

    RETURN
  END SUBROUTINE fft_backward2_inplace

  SUBROUTINE fft_exit
    IMPLICIT NONE

    CALL dfftw_destroy_plan(plan_forward1)
    CALL dfftw_destroy_plan(plan_p)
    CALL dfftw_destroy_plan(plan_j)
    CALL dfftw_destroy_plan(plan_backward2)
    CALL dfftw_destroy_plan(plan_spec)

    RETURN
  END SUBROUTINE fft_exit

END MODULE fft
