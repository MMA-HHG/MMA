PROGRAM test


IMPLICIT NONE


Complex*16,PARAMETER :: xi=(0.0d0,1.0d0);
real*8,PARAMETER :: pi=acos(-1.0d0);
Complex*16 :: c; 
real*8 :: I0, press, w0, focus, tau, n2atm, lambda;
Real*8 :: zR, Iz, w0z, Pcr, Pin, focal_length, n2p;
integer :: k,N;

read(*,*) I0, press, w0, focus, tau, n2atm, lambda;

zR = pi*w0*w0/(lambda);
w0z = w0*sqrt(1.0d0 + (focus/zR)**2 );
n2p = n2atm*press;
!Iz = I0*((w0/w0z)**2)*;
Pin = I0*w0*w0*pi/2.0d0;!in focus, independent on z
Pcr = 0.5d0*(lambda**2)/(pi*n2p);

focal_length = focus + zR*zR/focus;





open(UNIT=1,FILE="recalcul.dat",FORM="FORMATTED",action='write');
write(1,'(ES8.2)') Pin/Pcr; write(*,'(ES8.2)') Pin/Pcr;
write(1,*) w0z; write(*,*) w0z;
write(1,*) focal_length; write(*,*) focal_length;


 close(1);



!write(*,*) "the work is done";

END PROGRAM test
