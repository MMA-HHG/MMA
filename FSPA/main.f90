
!%%%%%%%%%%%%%%%%%  Program  %%%%%%%%%%%%%%
PROGRAM SDP

  use util; use tools; 

  IMPLICIT NONE


  !%%%%%%%%%  Calculation of the saddle points

  Complex*16 :: ti,tr,ps,dip_s,dip_l, phasePhi, dipole;
  Complex*16 :: ti_long, tr_long,t_exc;
  Complex*16 :: a11,a22,a12,detS;

  real*8 :: w0,phi,wq,Ip,period,order,eps,norm;
  
  real*8 :: A0,A01,A02,dA;
  real*8 :: dipole_mod, dipole_ph;

  real*8 :: ph_long_dum,ph_long0,ph_short_dum,ph_short0;

  integer :: nc,n_init,n_rec,j,nb_pts,i;

  integer :: st;
  
  
  open(UNIT=1,FILE="param.inp",FORM="FORMATTED", ACCESS="SEQUENTIAL", status='old', action='read');
  open(UNIT=2,FILE="times_long.dat",FORM="FORMATTED",action='write');
  open(UNIT=3,FILE="times_short.dat",FORM="FORMATTED",action='write');
  open(UNIT=4,FILE="phase_long.dat",FORM="FORMATTED",action='write');
  open(UNIT=5,FILE="phase_short.dat",FORM="FORMATTED",action='write');
  !open(UNIT=6,FILE="errors.dat",FORM="FORMATTED",action='write');
  open(UNIT=7,FILE="log.txt",FORM="FORMATTED",action='write');
  
  


  !open(UNIT=1,FILE="psi_end.dat",FORM="FORMATTED",action='write');
write(*,*) "test";
!!!INPUT PARAMETERS (all in the atomic  units)
  read(UNIT=1,fmt=*, IOSTAT=st) Ip; !Ionization potential of the ground state
  read(UNIT=1,fmt=*, IOSTAT=st) w0; !photon energy
  read(UNIT=1,fmt=*, IOSTAT=st) order; !harmonic order
  read(UNIT=1,fmt=*, IOSTAT=st) A01; !minimal vector potential
  read(UNIT=1,fmt=*, IOSTAT=st) dA; !minimal vector potential ! spacing in the vector potential
  read(UNIT=1,fmt=*, IOSTAT=st) nb_pts; ! # of points


 ! read(UNIT=1,fmt=*, IOSTAT=st) A02; !maximal vector potential
  
 ! read(UNIT=1,fmt=*, IOSTAT=st) nb_pts; !# of points
 
!  A01=0.1d0; ! 0.5d0
 ! order=55.0d0;
!  nb_pts=1000;

  phi=0.0d0;
  nc=1; ! # of cycles in the enevelope - not used now

  close(1);


  
   !Ip = 0.792482; !Neon 0.792482
   !A0 = 1.0d0; nc = 10; w0 = 0.0568d0; phi = 0.0d0; !n_init = 9; n_rec = 10; 
   
 
!!! COMPUTED PARAMETERS
   period = 2.0d0*pi/w0;
   wq=order*w0;

!!! POSSIBLE IMPROVEMENT
! the long trajectories looks as stable ones according to results. Next, the recollision times doesn't coincide. Thus, if a short trajectory converges to the "long-result", initial condition is modified and rerunned. If this procedure fails, an error is written into the output error file. 


!!! THIS LOOP CALCULATE PHASE IN GIVEN RANGE OF FREQUENCIES

! dA=0.01d0;!(A02-A01)/dfloat(nb_pts);
n_init=0; n_rec=1; !inherited from the previous code, used in initial conditions

write(*,*) Ip, w0, order;

write(7,*) "I_p, photon energy, harmonic order";
write(7,*) Ip, w0, order;
write(7,*)

ph_long_dum = 0.0d0; ph_long0 = 0.0d0; ph_short_dum = 0.0d0; ph_short0 = 0.0d0;

do j = 1, (nb_pts+1); 



	
	A0=A01+dfloat(j-1)*dA;
write(7,*) " ";
write(7,*) "j =",j, "I =", (A0*w0)**2.0d0, "A0=", A0;
write(7,*) "expected cutoff:", (3.17*A0**2.0d0/4.0d0 + Ip)/w0;


    !!! LONG TRAJECTORIES PART
!   write(*,*)
!   write(*,*)


!   write(*,*) "phi_i_init = ",ti*w0;
!   write(*,*) "phi_r_init = ",tr*w0;

	!!!!! Init for long trajectories
	ti = 0.15d0*period+2.0d0*period*xi + phi/w0; !! 0.15d0*period+2.0d0*period*xi+ dfloat(n_init)*period + phi/w0;
	tr = ti + 0.6d0*period - 0.0d0*period*xi; !! ti + 0.6d0*period + dfloat(n_rec-1)*period - 0.0d0*period*xi;

	call SDP_procedure(A0,nc,w0,phi,ti,tr,wq,Ip);
	ti_long=ti; ! ionization time for long traj.
	tr_long=tr; ! recombination time for long traj.
	ps = -int_A_vect(A0,nc,w0,phi,ti,tr)/(tr-ti); 
	phasePhi=wq*tr-actionS(A0,nc,w0,phi,ti,tr,wq,Ip,ps);

	!the determinant from the Saddle point part
	!detS=a12^2-a11*a22 (see the referenced article)
	t_exc=tr-ti; !excursion time
	a12 = (ps+A_vect(A0,nc,w0,phi,tr))*(ps+A_vect(A0,nc,w0,phi,ti))/t_exc ;
	a11 = 2.0d0*Ip/t_exc + E_vect(A0,nc,w0,phi,ti)*(ps+A_vect(A0,nc,w0,phi,ti)) ;
	a22 = -2.0d0*(wq-Ip)/t_exc - E_vect(A0,nc,w0,phi,tr)*(ps+A_vect(A0,nc,w0,phi,ti)) ;
	detS=a12**2.d0-a11*a22;
	

!!! REMARK: The transition matrix dipole elements has no significant influence for the amplitude, but what about the accumulated phase? I've found that SFA model is very sensitive for this issues (only change i -> -i in square root completely destroyed the result, etc.)

	dipole = E_vect(A0,nc,w0,phi,ti)*(t_exc**(-1.5d0))*exp(xi*PhasePhi)/sqrt(detS) ; !  E_vect(A0,nc,w0,phi,ti)*(t_exc**(-1.5d0))*exp(xi*PhasePhi)/sqrt(detS) ;!including significant terms from the first function e^(i*PhasePhi)
	dipole_mod = sqrt(real(dipole*conjg(dipole)));
	dipole_ph = atan2(dimag(dipole),real(dipole));

	!!!PHASE UNWRAP
	if ( ( abs(dipole_ph-ph_long_dum) > pi ) .and. j > 1 ) then
		ph_long0 = ph_long0 - SIGN(2.0d0*pi,(dipole_ph-ph_long_dum)); !ph_long0 = ph_long0 - 2.0d0*pi;		
	endif
	ph_long_dum=dipole_ph;
	dipole_ph=dipole_ph+ph_long0;
	 
	write(UNIT=2,fmt='(XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6)') (w0*A0)**2.0d0 ,real(w0*ti),& !recombination times for a given intensity (harmonic field)
	&dimag(w0*ti),real(w0*tr),dimag(w0*tr),real(ps),dimag(ps);
 
	
	write(UNIT=4,fmt='(XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6)') (w0*A0)**2.0d0 ,real(phasePhi), dimag(phasePhi),&
	&real(dipole), dimag(dipole), exp(-dimag(PhasePhi)) ,dipole_mod, dipole_ph; !! dipole phase for given times as a function of intensity (A0*w0)^2 and dipole, last two ough to be the same


!   write(*,*)

   write(7,*) "phi_i_sol = ",ti*w0;
   write(7,*) "phi_r_sol = ",tr*w0;

 !  write(*,*) "f = ",f(A0,nc,w0,phi,ti,tr,wq,Ip);
!   write(*,*) "g = ",g(A0,nc,w0,phi,ti,tr,wq,Ip);

!   write(*,*) "l1 = ",(ps+A_vect(A0,nc,w0,phi,tr))**2.0d0;
!   write(*,*) "l2 = ",2.0d0*(wq-Ip);

!   write(*,*) "l1 = ",(ps+A_vect(A0,nc,w0,phi,ti))**2.0d0;
!   write(*,*) "l2 = ",2.0d0*(-Ip);

!   write(*,*) "p à la recollision: ",int_A_vect(A0,nc,w0,phi,ti,tr)/(tr-ti); 




    !!! SHORT TRAJECTORIES PART
 !  write(*,*)
!   write(*,*)

!   write(*,*) "phi_i_init = ",ti*w0;
!   write(*,*) "phi_r_init = ",tr*w0;


	!!!!! Init for short trajectories	
	ti = 0.3d0*period+2.0d0*period*xi+ dfloat(n_init)*period + phi/w0; 
	tr = ti + 0.14d0*period + dfloat(n_rec-1)*period - 0.0d0*period*xi;

	call SDP_procedure(A0,nc,w0,phi,ti,tr,wq,Ip);

	!control if this solution is properly converged (see REMARK)
		!eps=1.0d-4; !ought to somehow* correspond with eps in the SDP_procedure (* we don't know characer of convergence 1/n, 1/n^2, etc.; 10-4 looks sufficiently from simulations
  
		!norm = real((ti_long-ti)*conjg(ti_long-ti)+(tr_long-tr)*conjg(tr_long-tr));
		!norm = w0*dsqrt(norm);

		!i=0;
		!do while (norm < eps)
			
			
		!	ti = 0.3d0*period+0.5d0*period*xi+ dfloat(n_init)*period + phi/w0; 
		!	tr = ti + (0.13d0+dfloat(i)*1.0d-2)*period + dfloat(n_rec-1)*period - 0.0d0*period*xi;
		!	call SDP_procedure(A0,nc,w0,phi,ti,tr,wq,Ip);
		!	norm = real((ti_long-ti)*conjg(ti_long-ti)+(tr_long-tr)*conjg(tr_long-tr));
		!	norm = w0*dsqrt(norm);
		!	if (i==20) then
		!		write(* , * ) "short traj. problem for j=",j,"vecpot A0=", A0;
		!		exit;
		!	end if;
		!	i=i+1;
		!enddo
	 !end of the control

	
 
	ps = -int_A_vect(A0,nc,w0,phi,ti,tr)/(tr-ti); 
	phasePhi=wq*tr-actionS(A0,nc,w0,phi,ti,tr,wq,Ip,ps);

	!the determinant from the Saddle point part
	!detS=a12^2-a11*a22 (see the referenced article)
	t_exc=tr-ti; !excursion time
	a12 = (ps+A_vect(A0,nc,w0,phi,tr))*(ps+A_vect(A0,nc,w0,phi,ti))/t_exc ;
	a11 = 2.0d0*Ip/t_exc + E_vect(A0,nc,w0,phi,ti)*(ps+A_vect(A0,nc,w0,phi,ti)) ;
	a22 = -2.0d0*(wq-Ip)/t_exc - E_vect(A0,nc,w0,phi,tr)*(ps+A_vect(A0,nc,w0,phi,ti)) ;
	detS=a12**2.d0-a11*a22;
	



	dipole = E_vect(A0,nc,w0,phi,ti)*(t_exc**(-1.5d0))*exp(xi*PhasePhi)/sqrt(detS);!including significant terms from the first function e^(i*PhasePhi) E_vect(A0,nc,w0,phi,ti)*(t_exc**(-1.5d0))*exp(xi*PhasePhi)/sqrt(detS)
	dipole_mod = sqrt(real(dipole*conjg(dipole)));
	dipole_ph = atan2(dimag(dipole),real(dipole));	


	!!!PHASE UNWRAP
	if ( ( abs(dipole_ph-ph_short_dum) > pi ) .and. j > 1 ) then
		ph_short0 = ph_short0 - SIGN(2.0d0*pi,(dipole_ph-ph_short_dum)); !ph_short0 = ph_short0 - 2.0d0*pi;		
	endif
	ph_short_dum=dipole_ph;
	dipole_ph=dipole_ph+ph_short0;

	
	write(UNIT=3,fmt='(XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6)') (w0*A0)**2.0d0 ,real(w0*ti),& !recombination times for a given intensity (harmonic field)
	&dimag(w0*ti),real(w0*tr),dimag(w0*tr),real(ps),dimag(ps);

	

	write(UNIT=5,fmt='(XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6,XE12.6)') (w0*A0)**2.0d0, real(phasePhi),  dimag(phasePhi),&
	&real(dipole), dimag(dipole), exp(-dimag(PhasePhi)), dipole_mod, dipole_ph; !! dipole phase for given times as a function of intensity (A0*w0)^2


!   write(*,*)

   write(7,*) "phi_i_sol = ",ti*w0;
   write(7,*) "phi_r_sol = ",tr*w0;

!   write(*,*) "f = ",f(A0,nc,w0,phi,ti,tr,wq,Ip);
!   write(*,*) "g = ",g(A0,nc,w0,phi,ti,tr,wq,Ip);

!   write(*,*) "l1 = ",(ps+A_vect(A0,nc,w0,phi,tr))**2.0d0;
!   write(*,*) "l2 = ",2.0d0*(wq-Ip);

 !  write(*,*) "l1 = ",(ps+A_vect(A0,nc,w0,phi,ti))**2.0d0;
!   write(*,*) "l2 = ",2.0d0*(-Ip);

!   write(*,*) "p à la recollision: ",int_A_vect(A0,nc,w0,phi,ti,tr)/(tr-ti); 


enddo


  close(1);close(2);close(3);close(4);close(5);

 ! write(*,*) "expected cutoff:", (3.17*A0**2.0d0/4 + Ip)/w0;
 
write(*,*) nc;
write(*,*) "the work is done"
END PROGRAM SDP














