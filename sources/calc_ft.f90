PROGRAM CALC_FT
  USE constants
  USE fft

  IMPLICIT NONE

  INTEGER(8) plan_ex,plan_ey
  REAL(4) k_r
  REAL(8) delta_r_inv
  REAL(8), ALLOCATABLE :: r(:)
  COMPLEX(4), ALLOCATABLE :: cbuffer(:,:)
  COMPLEX(8), ALLOCATABLE :: ex(:),ey(:)
  CHARACTER(10) filename
  LOGICAL existence

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,num_proc,ierr)
  CALL READPARA
  CALL fft_init
  existence=.FALSE.
  DO WHILE (.NOT.existence)
     CALL READFILE
     IF (.NOT.existence) THEN
        CALL COMPUTE_FT
        CALL WRITEFILE
     ENDIF
  ENDDO

  CALL fft_exit
  CALL MPI_Finalize(ierr)

CONTAINS

  SUBROUTINE READPARA
    IMPLICIT NONE

    INTEGER(4) j,i
    REAL(8) lambda0_cm_phys,tp_fs_phys,c,k0_phys
    CHARACTER(50) dummy

    IF (my_rank.EQ.0) THEN
       OPEN(unit=100,FILE='listing',form='formatted',STATUS='OLD')
       READ(100,*)
       READ(100,*)
       READ(100,'(a,t50,i5)') dummy,dim_t
       READ(100,'(a,t50,i5)') dummy,dim_r
       READ(100,'(a,t50,es12.4)') dummy,lt
       READ(100,'(a,t50,es12.4)') dummy,lr
       DO j=1,13
          READ(100,*)
       ENDDO
       READ(100,'(a,t50,es12.4)') dummy,lambda0_cm_phys
       READ(100,*)
       READ(100,*)
       READ(100,*)
       READ(100,'(a,t50,es12.4)') dummy,tp_fs_phys
       CLOSE(100)
       c = 3.d10 !speed of light in the vacuum      cm/s
       k0_phys = 2.D0*PI/lambda0_cm_phys  !central wave number in vacuum     cm-1
       omega=c*k0_phys*tp_fs_phys*1.d-15 ! adimensioned frequency
       i=1.D-1*omega*lt/(4.D0*DATAN(1.D0))
       omega_uppe=MAX(omega,4.D0*DATAN(1.D0)*REAL(dim_t+i)/lt)
    ENDIF
    CALL MPI_BCAST(dim_t,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(dim_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(omega_uppe,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    ALLOCATE (cbuffer(dim_t,dim_r),r(dim_r),ex(2*dim_r),ey(2*dim_r))
    k_t=8.D0*DATAN(1.D0)/lt
    k_r=8.D0*DATAN(1.D0)/(2.D0*lr)
    delta_r=lr/REAL(dim_r,8)
    delta_r_inv=1.D0/delta_r
    DO j=1,dim_r
       r(j)=delta_r*REAL(j-1)
    ENDDO
    CALL dfftw_plan_dft_1d(plan_ex,2*dim_r,ex,ex,-1,0)
    CALL dfftw_plan_dft_1d(plan_ey,2*dim_r,ey,ey,-1,0)

    RETURN
  END SUBROUTINE READPARA

  SUBROUTINE READFILE
    IMPLICIT NONE

    INTEGER(4) j,k

    IF (my_rank.EQ.0) THEN
       OPEN(10,FILE='FT_RAD.LOG',STATUS='OLD')
       DO
          existence=.TRUE.
          READ(10,*,END=999) filename
          INQUIRE(FILE=filename//'_field_k_omega.dat',EXIST=existence)
          IF (.NOT.existence) EXIT
       ENDDO
999    CLOSE(10)
       IF (.NOT.existence) THEN
          OPEN(11,FILE=filename//'_field_r_t.dat',STATUS='OLD',FORM='UNFORMATTED')
          READ(11) 
          READ(11)
          DO k=1,dim_r
             READ(11)
             READ(11) (cbuffer(j,k),j=1,dim_t)
          ENDDO
          CLOSE(11)
       ENDIF
    ENDIF
    CALL MPI_BCAST(existence,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    IF (.NOT.existence) THEN
       CALL MPI_BCAST(cbuffer,dim_t*dim_r,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
       e(1:dim_t,dim_r_start(num_proc):dim_r_end(num_proc))=CSHIFT(cbuffer(1:dim_t,dim_r_start(num_proc):dim_r_end(num_proc)),dim_t/2-1,1)
    ENDIF

    RETURN
  END SUBROUTINE READFILE

  SUBROUTINE COMPUTE_FT
    IMPLICIT NONE

    INTEGER(4) k,l

    switch_T=1
    CALL fft_forward_inplace
    DO k=dim_t_start(num_proc),dim_t_end(num_proc)
       DO l=1,dim_r
          CALL INTERPOL(efft(:,k),ex,r(l))
          CALL dfftw_execute(plan_ex)
          ey(dim_r+l-1)=ex(1)
          ey(dim_r-l+1)=ey(dim_r+l-1)
       ENDDO
       ey(2*dim_r)=ey(2*dim_r-1)
       CALL dfftw_execute(plan_ey)
       efft(:,k)=ey(1:dim_r)
    ENDDO

    RETURN
  END SUBROUTINE COMPUTE_FT

  SUBROUTINE INTERPOL(vecin,vecout,y)
    IMPLICIT NONE

    REAL(8) y
    COMPLEX(8) vecin(dim_r),vecout(2*dim_r)

    INTEGER(4) i,j,k
    REAL(8) rint

    j=2
    DO k=1,dim_r
       rint=sqrt(r(k)**2+y**2)
       IF(rint.GT.r(dim_r)) THEN
          vecout(dim_r+k-1)=CMPLX(0.D0,0.D0,8) !vecin(dim_r)
       ELSE IF(rint.LE.0.D0) THEN
          vecout(dim_r)=vecin(1)
       ELSE
          DO i=j,dim_r
             IF(r(i).GE.rint) THEN
                vecout(dim_r+k-1)=((rint-r(i-1))*vecin(i)+(r(i)-rint)*vecin(i-1))*delta_r_inv
                j=i
                EXIT
             ENDIF
          ENDDO
       ENDIF
       vecout(dim_r-k+1)=vecout(dim_r+k-1)
    ENDDO
    vecout(2*dim_r)=vecout(2*dim_r-1)

    RETURN
  END SUBROUTINE INTERPOL

  SUBROUTINE WRITEFILE
    IMPLICIT NONE

    INTEGER(4) l,k,j
    COMPLEX(4) help

    DO l=1,dim_r
       DO k=dim_t_start(num_proc),dim_t_end(num_proc)
          cbuffer(k,l)=efft(l,k)
       ENDDO
       CALL MPI_GATHER(cbuffer(dim_t_start(num_proc),l),dim_t_local,MPI_COMPLEX,cbuffer(1,l),dim_t_local,MPI_COMPLEX,0,MPI_COMM_WORLD,ierr)
    ENDDO
    IF (my_rank.EQ.0) THEN
       OPEN(11,FILE=filename//'_field_k_omega.dat',STATUS='NEW',FORM='UNFORMATTED')
       WRITE(11) dim_t,dim_r
       WRITE(11) (REAL(omega_uppe,4)+REAL(j-dim_t/2-1,4)*REAL(k_t,4),j=1,dim_t)
       DO k=1,dim_r 
          DO j=1,dim_t/2
             help=cbuffer(j+dim_t/2,k)
             cbuffer(j+dim_t/2,k)=cbuffer(j,k)
             cbuffer(j,k)=help
          ENDDO
          WRITE(11) k_r*REAL(k-1,4)
          WRITE(11) (cbuffer(j,k),j=1,dim_t)
       ENDDO
       CLOSE(11)
    ENDIF

    RETURN
  END SUBROUTINE WRITEFILE
END PROGRAM CALC_FT
