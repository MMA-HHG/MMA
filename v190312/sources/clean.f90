PROGRAM CLEAN
  IMPLICIT NONE

  INTEGER(4) l,p,num_proc
  CHARACTER(3) ip
  CHARACTER*10 filename,save1,save2
  LOGICAL existence

  OPEN(10,FILE='FLUENCE_000.DAT',STATUS='OLD',FORM='UNFORMATTED')
  READ(10) num_proc
  CLOSE(10)
  OPEN(10,FILE='PROP_RAD.LOG',STATUS='OLD')
  DO
     READ(10,*,END=888) filename
     save2=save1
     save1=filename
  ENDDO
888 CLOSE(10)
  OPEN(10,FILE='PROP_RAD.LOG',STATUS='OLD')
  DO
     READ(10,*,END=999) filename
     INQUIRE(FILE=filename//'_000.DAT',EXIST=existence)
     IF (existence.AND.(filename.NE.save1).AND.(filename.NE.save2)) THEN
        PRINT*, filename
        DO p=0,num_proc-1
           WRITE(ip,930) p
           DO l=1,3
              IF (ip(l:l).EQ.' ') ip(l:l)='0'
           ENDDO
           OPEN(11,FILE=filename//'_'//ip//'.DAT',STATUS='OLD',FORM='UNFORMATTED')
           CLOSE(11,STATUS='DELETE')
        ENDDO
     ENDIF
  ENDDO
999 CLOSE(10)

930 FORMAT (I3) 

END PROGRAM CLEAN
