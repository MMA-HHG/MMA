module util
Complex*16,PARAMETER :: xi=(0.0d0,1.0d0);
real*8,PARAMETER :: pi = acos(-1.0d0);
end module util


module tools

use util;

IMPLICIT NONE;

CONTAINS

!%%%%%%%%% Definition of the fields
!Complex*16 function VecPot(A0,phi,t)

!	real*8,intent(in) :: A0, phi;
!	complex*16,intent(in) :: t;

!end function Vecpot


!%%%%%%%%%   Calculation of the Action
Complex*16  FUNCTION actionS(A0,nc,w0,phi,ti,tr,wq,Ip,ps)


  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr,ps;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;
  
  actionS = (Ip+0.5d0*ps*ps)*(tr-ti)+ps*int_A_vect(A0,nc,w0,phi,ti,tr)+0.5d0*int_A2_vect(A0,nc,w0,phi,ti,tr);

end function actionS


!%%%%%%%%%   Calculation of the vector potential
Complex*16  FUNCTION dipz(p)

  Complex*16,intent(in) :: p;

  real*8 :: alpha;


  alpha = 0.8d0*0.5d0; 
  
  !dipz = 4.0d0*sqrt(6.0d0)*pi**(-1.5d0)*p*(1.0d0+p*p)**(-3.0d0);
  dipz = xi*(alpha*pi)**(-3.0d0/4.0d0)*p*exp(-p*p*(2.0d0*alpha)**(-1.0d0))/alpha;

end function dipz

!%%%%%%%%%   Calculation of the determinant of the action
Complex*16  FUNCTION det_S(A0,nc,w0,phi,ti,tr,wq,Ip,ps)

  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr,ps;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;

  !det_S = (ps+A_vect(A0,nc,w0,phi,ti))*(ps+A_vect(A0,nc,w0,phi,tr))/(tr-ti);
  !det_S = det_S*det_S-(-2.0d0*(wq-Ip)/(tr-ti)-E_vect(A0,nc,w0,phi,tr)*(ps+A_vect(A0,nc,w0,phi,tr)))*(2.0d0*Ip/(tr-ti)+E_vect(A0,nc,w0,phi,ti)*(ps+A_vect(A0,nc,w0,phi,ti)));

  det_S = E_vect(A0,nc,w0,phi,tr)*(ps+A_vect(A0,nc,w0,phi,tr))-E_vect(A0,nc,w0,phi,ti)*(ps+A_vect(A0,nc,w0,phi,ti));  

end function det_S


!%%%%%%%%%   Calculation of the vector potential
Complex*16  FUNCTION A_vect(A0,nc,w0,phi,t)

  !IMPLICIT NONE;

  integer, intent(in) :: nc;
  Complex*16,intent(in) :: t;
  real*8,intent(in) :: A0,w0,phi;
  
  real*8 :: wc;

  wc = w0/(2.0d0*dfloat(nc));

  !A_vect = A0*cos(w0*t+phi)*(sin(wc*t))**(2.0d0);

  A_vect = -A0*sin(w0*t);

end function A_vect


!%%%%%%%%%   Calculation of the elctric field potential
Complex*16  FUNCTION E_vect(A0,nc,w0,phi,t)

  !IMPLICIT NONE;

  integer, intent(in) :: nc;
  Complex*16,intent(in) :: t;
  real*8,intent(in) :: A0,w0,phi;
  
  real*8 :: wc;

  wc = w0/(2.0d0*dfloat(nc));

  !E_vect = A0*w0*sin(wc*t)*(sin(wc*t)*sin(w0*t+phi)-2.0d0*(wc/w0)*cos(wc*t)*cos(w0*t+phi));

  E_vect = A0*w0*cos(w0*t);

end function E_vect


!%%%%%%%%%   Calculation of the integral of the vector potential
Complex*16  FUNCTION int_A_vect(A0,nc,w0,phi,t1,t2)

  !IMPLICIT NONE;

  integer, intent(in) :: nc;
  Complex*16, intent(in) :: t1,t2;
  real*8, intent(in) :: A0,w0,phi;
  
  real*8 :: wc;

  wc = w0/(2.0d0*dfloat(nc));

  !write(*,*) "w0 :",w0;

  !int_A_vect = A0*0.5d0*(sin(w0*t2+phi)-sin(w0*t1+phi))/w0;
  !int_A_vect = int_A_vect-A0*0.25d0*(sin((2.0d0*wc+w0)*t2+phi)-sin((2.0d0*wc+w0)*t1+phi))/(2.0d0*wc+w0);
  !int_A_vect = int_A_vect-A0*0.25d0*(sin((2.0d0*wc-w0)*t2-phi)-sin((2.0d0*wc-w0)*t1-phi))/(2.0d0*wc-w0);

  int_A_vect = (A0/w0)*(cos(w0*t2)-cos(w0*t1));

end function int_A_vect


!%%%%%%%%%   Calculation of the integral of the vector potential squared
Complex*16  FUNCTION int_A2_vect(A0,nc,w0,phi,t1,t2)

  !IMPLICIT NONE;

  integer, intent(in) :: nc;
  Complex*16, intent(in) :: t1,t2;
  real*8, intent(in) :: A0,w0,phi;
  
  real*8 :: wc;

  wc = w0/(2.0d0*dfloat(nc));

  !write(*,*) "w0 :",w0;


  !int_A2_vect = 3.0d0*(t2-t1)/16.0d0;
  !int_A2_vect = int_A2_vect - 2.0d0*(sin(2.0d0*wc*t2)-sin(2.0d0*wc*t1))/(16.0d0*wc);
  !int_A2_vect = int_A2_vect + (sin(4.0d0*wc*t2)-sin(4.0d0*wc*t1))/(16.0d0*4.0d0*wc);
  !int_A2_vect = int_A2_vect + 3.0d0*(sin(2.0d0*w0*t2+2.0d0*phi)-sin(2.0d0*w0*t1+2.0d0*phi))/(16.0d0*2.0d0*w0);
  !int_A2_vect = int_A2_vect - (sin(2.0d0*(w0+wc)*t2+2.0d0*phi)-sin(2.0d0*(w0+wc)*t1+2.0d0*phi))/(16.0d0*(wc+w0));
  !int_A2_vect = int_A2_vect - (sin(2.0d0*(-w0+wc)*t2-2.0d0*phi)-sin(2.0d0*(-w0+wc)*t1-2.0d0*phi))/(16.0d0*(wc-w0));
  !int_A2_vect = int_A2_vect + (sin(2.0d0*(w0+2.0d0*wc)*t2+2.0d0*phi)-sin(2.0d0*(w0+2.0d0*wc)*t1+2.0d0*phi))/(64.0d0*(2.0d0*wc+w0));
  !int_A2_vect = int_A2_vect + (sin(2.0d0*(-w0+2.0d0*wc)*t2-2.0d0*phi)-sin(2.0d0*(-w0+2.0d0*wc)*t1-2.0d0*phi))/(64.0d0*(2.0d0*wc-w0));
  !int_A2_vect = A0*A0*int_A2_vect;


  int_A2_vect = 0.5d0*A0*A0*(t2-t1)-A0*A0*(sin(2.0d0*w0*t2)-sin(2.0d0*w0*t1))/(4.0d0*w0);

end function int_A2_vect


!%%%%%%%%%   Calculation of 
Complex*16  FUNCTION f(A0,nc,w0,phi,ti,tr,wq,Ip)

  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;
  

  f = (int_A_vect(A0,nc,w0,phi,ti,tr)/(tr-ti)-A_vect(A0,nc,w0,phi,tr))**(2.0d0)-2.0d0*(wq-Ip);

end function f


!%%%%%%%%%   Calculation of 
Complex*16  FUNCTION g(A0,nc,w0,phi,ti,tr,wq,Ip)


  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;
  

  g = (int_A_vect(A0,nc,w0,phi,ti,tr)/(tr-ti)-A_vect(A0,nc,w0,phi,ti))**(2.0d0)+2.0d0*Ip;

end function g


!%%%%%%%%%   Calculation of 
Complex*16  FUNCTION der_f_x(A0,nc,w0,phi,ti,tr,wq,Ip)

  

  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;


  !write(*,*) A0,nc,w0,phi,ti,tr
  !write(*,*) int_A_vect(A0,nc,w0,phi,ti,tr)

  !der_f_x = 2.0d0*(int_A_vect(A0,nc,w0,phi,ti,tr)-(tr-ti)*A_vect(A0,nc,w0,phi,tr))*(A_vect(A0,nc,w0,phi,tr)-A_vect(A0,nc,w0,phi,ti));
  !der_f_x = der_f_x + 4.0d0*(tr-ti)*(wq-Ip);

  der_f_x = 2.0d0*(-A_vect(A0,nc,w0,phi,ti)*(tr-ti)+int_A_vect(A0,nc,w0,phi,ti,tr))
  der_f_x = der_f_x*(int_A_vect(A0,nc,w0,phi,ti,tr)-A_vect(A0,nc,w0,phi,tr)*(tr-ti));
  

  der_f_x = der_f_x/(tr-ti)**(3.0d0);

end function der_f_x


!%%%%%%%%%   Calculation of 
Complex*16  FUNCTION der_f_y(A0,nc,w0,phi,ti,tr,wq,Ip)


  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;
  
  !der_f_y = 2.0d0*(int_A_vect(A0,nc,w0,phi,ti,tr)-(tr-ti)*A_vect(A0,nc,w0,phi,tr))*(tr-ti)*E_vect(A0,nc,w0,phi,tr);
  !der_f_y = der_f_y - 4.0d0*(tr-ti)*(wq-Ip);

  der_f_y = 2.0d0*(A_vect(A0,nc,w0,phi,tr)*(tr-ti)-int_A_vect(A0,nc,w0,phi,ti,tr)+E_vect(A0,nc,w0,phi,tr)*(tr-ti)**(2.0d0));
  der_f_y = der_f_y*(int_A_vect(A0,nc,w0,phi,ti,tr)-A_vect(A0,nc,w0,phi,tr)*(tr-ti));

  der_f_y = der_f_y/(tr-ti)**(3.0d0);

end function der_f_y


!%%%%%%%%%   Calculation of 
Complex*16  FUNCTION der_g_x(A0,nc,w0,phi,ti,tr,wq,Ip)


  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;

  !der_g_x = 2.0d0*(int_A_vect(A0,nc,w0,phi,ti,tr)-(tr-ti)*A_vect(A0,nc,w0,phi,ti))*(tr-ti)*(E_vect(A0,nc,w0,phi,ti));
  !der_g_x = der_g_x - 4.0d0*(tr-ti)*Ip; 

  der_g_x = 2.0d0*(-A_vect(A0,nc,w0,phi,ti)*(tr-ti)+int_A_vect(A0,nc,w0,phi,ti,tr)+E_vect(A0,nc,w0,phi,ti)*(tr-ti)**(2.0d0));
  der_g_x = der_g_x*(int_A_vect(A0,nc,w0,phi,ti,tr)-A_vect(A0,nc,w0,phi,ti)*(tr-ti));

  der_g_x = der_g_x/(tr-ti)**(3.0d0);


end function der_g_x


!%%%%%%%%%   Calculation of 
Complex*16  FUNCTION der_g_y(A0,nc,w0,phi,ti,tr,wq,Ip)


  integer, intent(in) :: nc;
  Complex*16,intent(in) :: ti,tr;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;

  !der_g_y = 2.0d0*(int_A_vect(A0,nc,w0,phi,ti,tr)-(tr-ti)*A_vect(A0,nc,w0,phi,ti))*(A_vect(A0,nc,w0,phi,tr)-A_vect(A0,nc,w0,phi,ti));
  !der_g_y = der_g_y + 4.0d0*(tr-ti)*Ip;

  der_g_y = 2.0d0*(A_vect(A0,nc,w0,phi,tr)*(tr-ti)-int_A_vect(A0,nc,w0,phi,ti,tr));
  der_g_y = der_g_y*(int_A_vect(A0,nc,w0,phi,ti,tr)-A_vect(A0,nc,w0,phi,ti)*(tr-ti));

  der_g_y = der_g_y/(tr-ti)**(3.0d0);

end function der_g_y



subroutine SDP_procedure(A0,nc,w0,phi,ti,tr,wq,Ip)

  integer, intent(in) :: nc;
  Complex*16 :: ti,tr;
  real*8,intent(in) :: A0,w0,phi,wq,Ip;

  Complex*16 :: det,xxi,yi; 
  integer :: i;
  real*8 :: eps,norm;

  complex*16 :: gdum,fdum,fxdum,fydum,gxdum,gydum; 

  xxi = ti; yi = tr;

  !write(*,*) "init : ", xxi,yi;

  eps = 1.0d-10; i=0; norm = 1.0d0;

  do while (norm.gt.eps)

  !do i = 1 , 80
	
	gdum=g(A0,nc,w0,phi,ti,tr,wq,Ip);
	fdum=f(A0,nc,w0,phi,ti,tr,wq,Ip);
	fxdum=der_f_x(A0,nc,w0,phi,ti,tr,wq,Ip);
	fydum=der_f_y(A0,nc,w0,phi,ti,tr,wq,Ip);
	gxdum=der_g_x(A0,nc,w0,phi,ti,tr,wq,Ip);
	gydum=der_g_y(A0,nc,w0,phi,ti,tr,wq,Ip);

     det = fxdum*gydum;
     det = det - gxdum*fydum;

     det = det + 1.0d-20;

     xxi = xxi - (gydum*fdum-fydum*gdum)/det;
     yi = yi - (-gxdum*fdum+fxdum*gdum)/det;


    !write(*,*) "new ti:",xxi;
    !write(*,*) "new tr:",yi;
    !write(*,*) "det:",det;
    !write(*,*) "";
    
    norm = real((xxi-ti)*conjg(xxi-ti)+(yi-tr)*conjg(yi-tr));
    norm = dsqrt(norm);

    ti = xxi; tr = yi;  i = i+1;
  
    !write(*,*) "sol:",xxi,det,norm;

    if ( i == 10**4 ) then !!! output to screen if it doesn't converge
    write(*,*) "not converged A0="
    write(*,*) A0, i;
    exit
    end if


  enddo

  write(7,*) "# of steps for convergence", i;

  !write(*,*) "number of iterations: ",i;

end subroutine SDP_procedure





end module tools
