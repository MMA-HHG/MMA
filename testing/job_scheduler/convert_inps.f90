PROGRAM test


IMPLICIT NONE


Complex*16,PARAMETER :: xi=(0.0d0,1.0d0);
real*8,PARAMETER :: pi=acos(-1.0d0);
Complex*16 :: c; 
!real*4 :: lambdaf, w0zf, focal_lengthf;
real*16 :: lambda16, w0z16, focal_length16;
real*8 :: lambda, w0z, focal_length, tau;

integer :: k,N;

!read(*,*)  lambda16, w0z16, focal_length16;
!lambda=real(lambda16,8); w0z=real(w0z16,8); focal_length=real(focal_length16,8);

read(*,*)  lambda, w0z, focal_length, tau;
 

!write(*,*)  lambda, w0z, focal_length;

open(UNIT=1,FILE="converted.dat",FORM="FORMATTED",action='write');

write(1,'(ES8.2)') 100.0d0*lambda; write(*,'(ES8.2)') 100.0d0*lambda;
write(1,'(ES8.2)') 100.d0*w0z; write(*,'(ES8.2)') 100.d0*w0z;
write(1,'(ES8.2)') 100.0d0*focal_length; write(*,'(ES8.2)') 100.0d0*focal_length;
write(1,'(ES8.2)') tau; write(*,'(ES8.2)') tau;

 close(1);



!write(*,*) "the work is done";

END PROGRAM test
