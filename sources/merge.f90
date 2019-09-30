PROGRAM MERGE
  IMPLICIT NONE

  INTEGER(4) j,k,l,p,num_proc,dim_t,dim_r,dim_z, zcount
  REAL(4) delta_t,delta_r,tlo, zmax
  REAL(4), ALLOCATABLE :: buffer(:,:)
  COMPLEX(4), ALLOCATABLE :: cbuffer(:,:)
  CHARACTER(3) ip
  CHARACTER(10) filename
  LOGICAL existence


  OPEN(10,FILE='MERGE_RAD.LOG',STATUS='OLD')
  DO
     READ(10,*,END=999) filename
     CALL CONVERT('_PLASMA_','_plasma',.FALSE.)
     CALL CONVERT('_FIELD_','_field_r_t',.TRUE.)
  ENDDO
999 CLOSE(10)

  OPEN(11,FILE='FLUENCE_000.DAT',STATUS='OLD',FORM='UNFORMATTED')
  OPEN(12,FILE='fluence.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
  READ(11) num_proc,dim_r
  READ(11)
  dim_z=0
  DO
     READ(11,END=990)
     dim_z=dim_z+1
  ENDDO
990 CONTINUE
  CLOSE(11)
  ALLOCATE (buffer(dim_r+1,dim_z+1))
  DO p=0,num_proc-1
     WRITE(ip,'(I3)') p
     DO l=1,3
        IF (ip(l:l).EQ.' ') ip(l:l)='0'
     ENDDO
     OPEN(11,FILE='FLUENCE_'//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
     READ(11)
     READ(11) buffer(p*dim_r/num_proc+1:(p+1)*dim_r/num_proc,1)
     DO j=2,dim_z+1
        READ(11) buffer(1,j),buffer(p*dim_r/num_proc+2:(p+1)*dim_r/num_proc+1,j)
     ENDDO
     CLOSE(11)
  ENDDO
  zmax = MAX( buffer(1,2), 0.D0)
  zcount = 1
  DO j = 3, dim_z+1
     IF( buffer(1,j).GT.zmax) THEN
        zcount = zcount + 1
        zmax = MAX(buffer(1,j), zmax)
     ENDIF
  ENDDO
  WRITE(12) dim_r,zcount
  WRITE(12) buffer(1:dim_r,1)                ! r values
  WRITE(12) buffer(1,2)                      ! z = 0
  zmax = MAX( buffer(1,2), 0.)
  WRITE(12) buffer(2:dim_r+1,2)              ! fluence(z=0)
  DO k=3,dim_z+1
     IF ( buffer(1,k).GT.zmax )  THEN
        WRITE(12) buffer(1,k)                ! z
        WRITE(12) buffer(2:dim_r+1,k)        ! fluence (z)
        zmax = MAX( buffer(1,k), zmax)
     ENDIF
  ENDDO
  CLOSE(12)
  DEALLOCATE(buffer)

  OPEN(11,FILE='PLASMACHANNEL_000.DAT',STATUS='OLD',FORM='UNFORMATTED')
  OPEN(12,FILE='plasmachannel.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
  READ(11) num_proc,dim_r
  READ(11)
  dim_z=0
  DO
     READ(11,END=992)
     dim_z=dim_z+1
  ENDDO
992 CONTINUE
  CLOSE(11)
  ALLOCATE (buffer(dim_r+1,dim_z+1))
  DO p=0,num_proc-1
     WRITE(ip,'(I3)') p
     DO l=1,3
        IF (ip(l:l).EQ.' ') ip(l:l)='0'
     ENDDO
     OPEN(11,FILE='PLASMACHANNEL_'//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
     READ(11)
     READ(11) buffer(p*dim_r/num_proc+1:(p+1)*dim_r/num_proc,1)
     DO j=2,dim_z+1
        READ(11) buffer(1,j),buffer(p*dim_r/num_proc+2:(p+1)*dim_r/num_proc+1,j)
     ENDDO
     CLOSE(11)
  ENDDO
  zmax = MAX( buffer(1,2), 0.D0)
  zcount = 1
  DO j = 3, dim_z+1
     IF( buffer(1,j).GT.zmax) THEN
        zcount = zcount + 1
        zmax = MAX(buffer(1,j), zmax)
     ENDIF
  ENDDO
  WRITE(12) dim_r,zcount
  WRITE(12) buffer(1:dim_r,1)                ! r values
  WRITE(12) buffer(1,2)                      ! z = 0
  zmax = MAX( buffer(1,2), 0.)
  WRITE(12) buffer(2:dim_r+1,2)              ! fluence(z=0)
  DO k=3,dim_z+1
     IF ( buffer(1,k).GT.zmax )  THEN
        WRITE(12) buffer(1,k)                ! z
        WRITE(12) buffer(2:dim_r+1,k)        ! fluence (z)
        zmax = MAX( buffer(1,k), zmax)
     ENDIF
  ENDDO
  CLOSE(12)
  DEALLOCATE(buffer)

  OPEN(11,FILE='LOSSES_PLASMA_000.DAT',STATUS='OLD',FORM='UNFORMATTED')
  OPEN(12,FILE='losses_plasma.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
  READ(11) num_proc,dim_r
  READ(11)
  dim_z=0
  DO
     READ(11,END=993)
     dim_z=dim_z+1
  ENDDO
993 CONTINUE
  CLOSE(11)
  ALLOCATE (buffer(dim_r+1,dim_z+1))
  DO p=0,num_proc-1
     WRITE(ip,'(I3)') p
     DO l=1,3
        IF (ip(l:l).EQ.' ') ip(l:l)='0'
     ENDDO
     OPEN(11,FILE='LOSSES_PLASMA_'//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
     READ(11)
     READ(11) buffer(p*dim_r/num_proc+1:(p+1)*dim_r/num_proc,1)
     DO j=2,dim_z+1
        READ(11) buffer(1,j),buffer(p*dim_r/num_proc+2:(p+1)*dim_r/num_proc+1,j)
     ENDDO
     CLOSE(11)
  ENDDO
  zmax = MAX( buffer(1,2), 0.D0)
  zcount = 1
  DO j = 3, dim_z+1
     IF( buffer(1,j).GT.zmax) THEN
        zcount = zcount + 1
        zmax = MAX(buffer(1,j), zmax)
     ENDIF
  ENDDO
  WRITE(12) dim_r,zcount
  WRITE(12) buffer(1:dim_r,1)                ! r values
  WRITE(12) buffer(1,2)                      ! z = 0
  zmax = MAX( buffer(1,2), 0.)
  WRITE(12) buffer(2:dim_r+1,2)              ! fluence(z=0)
  DO k=3,dim_z+1
     IF ( buffer(1,k).GT.zmax )  THEN
        WRITE(12) buffer(1,k)                ! z
        WRITE(12) buffer(2:dim_r+1,k)        ! fluence (z)
        zmax = MAX( buffer(1,k), zmax)
     ENDIF
  ENDDO
  CLOSE(12)
  DEALLOCATE(buffer)

  OPEN(11,FILE='LOSSES_IONIZATION_000.DAT',STATUS='OLD',FORM='UNFORMATTED')
  OPEN(12,FILE='losses_ionization.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
  READ(11) num_proc,dim_r
  READ(11)
  dim_z=0
  DO
     READ(11,END=994)
     dim_z=dim_z+1
  ENDDO
994 CONTINUE
  CLOSE(11)
  ALLOCATE (buffer(dim_r+1,dim_z+1))
  DO p=0,num_proc-1
     WRITE(ip,'(I3)') p
     DO l=1,3
        IF (ip(l:l).EQ.' ') ip(l:l)='0'
     ENDDO
     OPEN(11,FILE='LOSSES_IONIZATION_'//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
     READ(11)
     READ(11) buffer(p*dim_r/num_proc+1:(p+1)*dim_r/num_proc,1)
     DO j=2,dim_z+1
        READ(11) buffer(1,j),buffer(p*dim_r/num_proc+2:(p+1)*dim_r/num_proc+1,j)
     ENDDO
     CLOSE(11)
  ENDDO
  zmax = MAX( buffer(1,2), 0.D0)
  zcount = 1
  DO j = 3, dim_z+1
     IF( buffer(1,j).GT.zmax) THEN
        zcount = zcount + 1
        zmax = MAX(buffer(1,j), zmax)
     ENDIF
  ENDDO
  WRITE(12) dim_r,zcount
  WRITE(12) buffer(1:dim_r,1)                ! r values
  WRITE(12) buffer(1,2)                      ! z = 0
  zmax = MAX( buffer(1,2), 0.)
  WRITE(12) buffer(2:dim_r+1,2)              ! fluence(z=0)
  DO k=3,dim_z+1
     IF ( buffer(1,k).GT.zmax )  THEN
        WRITE(12) buffer(1,k)                ! z
        WRITE(12) buffer(2:dim_r+1,k)        ! fluence (z)
        zmax = MAX( buffer(1,k), zmax)
     ENDIF
  ENDDO
  CLOSE(12)
  DEALLOCATE(buffer)

  OPEN(11,FILE='ONAX_T.DAT',STATUS='OLD',FORM='UNFORMATTED')
  OPEN(13,FILE='intensity_onax_t.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
  READ(11) dim_t
  READ(11)
  dim_z=0
  DO
     READ(11,END=991)
     dim_z=dim_z+1
  ENDDO
991 CONTINUE
  CLOSE(11)
  ALLOCATE (buffer(dim_t+1,dim_z+1))
  OPEN(11,FILE='ONAX_T.DAT',STATUS='OLD',FORM='UNFORMATTED')
  READ(11) dim_t
  READ(11) buffer(2:dim_t+1,1)
  zmax = 0.d0
  zcount = 0
  DO j=1,dim_z
     READ(11) buffer(1:dim_t+1,j+1)
     IF ( buffer(1,j+1).GT.zmax) THEN
        zcount = zcount + 1
        zmax = MAX( buffer(1,j+1), zmax)
     ENDIF
  ENDDO
  CLOSE(11)
  WRITE(13) dim_t,zcount
  WRITE(13) buffer(2:dim_t+1,1)
  zmax = 0.D0
  DO j=1,dim_z
     IF(  buffer(1,j+1).GT.zmax) THEN
        WRITE(13) buffer(1,j+1)
        WRITE(13) buffer(2:dim_t+1,j+1)
        zmax = MAX( buffer(1,j+1), zmax)
     ENDIF
  ENDDO
  CLOSE(13)
  DEALLOCATE(buffer)

CONTAINS

  SUBROUTINE CONVERT(inname,outname,complex)
    IMPLICIT NONE

    CHARACTER(*) inname,outname
    LOGICAL complex

    INQUIRE(FILE=filename//inname//'000.DAT',EXIST=existence)
    IF (.NOT.existence) GOTO 100
    OPEN(11,FILE=filename//inname//'000.DAT',STATUS='OLD',FORM='UNFORMATTED')
    READ(11) dim_t,dim_r,num_proc
    READ(11) delta_t,delta_r,tlo
    IF (complex) THEN
       ALLOCATE (cbuffer(dim_t,dim_r))
       DO j=1,dim_r/num_proc
          READ(11) cbuffer(1:dim_t,j)
       ENDDO
    ELSE
       ALLOCATE (buffer(dim_t,dim_r))
       DO j=1,dim_r/num_proc
          READ(11) buffer(1:dim_t,j)
       ENDDO
    ENDIF
    CLOSE(11)
    DO p=1,num_proc-1
       WRITE(ip,'(I3)') p
       DO l=1,3
          IF (ip(l:l).EQ.' ') ip(l:l)='0'
       ENDDO
       OPEN(11,FILE=filename//inname//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
       READ(11)
       READ(11)
       IF (complex) THEN
          DO j=dim_r/num_proc*p+1,dim_r/num_proc*(p+1)
             READ(11) cbuffer (1:dim_t,j)
          ENDDO
       ELSE
          DO j=dim_r/num_proc*p+1,dim_r/num_proc*(p+1)
             READ(11) buffer (1:dim_t,j)
          ENDDO
       ENDIF
       CLOSE(11)
    ENDDO
    OPEN(11,FILE=filename//outname//'.dat',STATUS='NEW',FORM='UNFORMATTED')
    WRITE(11) dim_t,dim_r
    WRITE(11) (tlo+REAL(j,4)*delta_t,j=1,dim_t)
    IF (complex) THEN
       DO k=1,dim_r
          WRITE(11) REAL(k-1,4)*delta_r
          WRITE(11) (cbuffer(j,k),j=1,dim_t)
       ENDDO
    ELSE
       DO k=1,dim_r
          WRITE(11) REAL(k-1,4)*delta_r
          WRITE(11) (buffer(j,k),j=1,dim_t)
       ENDDO
    ENDIF
    CLOSE(11)
    DO p=0,num_proc-1
       WRITE(ip,'(I3)') p
       DO l=1,3
          IF (ip(l:l).EQ.' ') ip(l:l)='0'
       ENDDO
       OPEN(11,FILE=filename//inname//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
       CLOSE(11,STATUS='DELETE')
    ENDDO
    IF (complex) THEN
       DEALLOCATE(cbuffer)
    ELSE
       DEALLOCATE(buffer)
    ENDIF
100 CONTINUE

    RETURN
  END SUBROUTINE CONVERT

END PROGRAM MERGE
