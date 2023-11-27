!> @brief This module contains the lookup subroutine \ref findinterval and
!! the linear interpolation subroutine \ref interpolate_lin, the procedures are
!! available up to 2 dimensions. Optionally, bookkeeping is possible (see
!! their descriptions for details).
!! 
!! Precisely, the subroutines are available via the aforementioned generic 
!! interfaces. For \ref interpolate_lin, the specific subroutines \ref
!! interpolate1D_lin and \ref interpolate2D_lin are available as well.
!! For the others, suppplementary procedures are encapsulatd within this
!! module and not publicly visible.
!!
!! @author Jan Vábek
module array_helper

implicit none
private
public  :: findinterval, interpolate_lin
public  :: interpolate2D_lin, interpolate1D_lin
public  :: interpolate1D_decomposed_eq
public  :: interpolate2D_decomposed_eq


!> @brief returns interval where x0-value is placed, ordering <..)<..)..<..>, extrapolation 0, n
!!
!! This suboroutine finds the (left) index of the interval where 'x0' is placed within the array 'x'.
!! The arguments are simply doubled for the 2D case (kx, ky, x0, y0, ...), (details in \ref findinterval_1D
!! and \ref findinterval_2D). \n 
!!  It uses the bisection lookup and initial guess(es) 'k_guess' can be used.
!!
!! @param[out]      k1          the seeked index
!! @param[in]       x0          
!! @param[in]       x           array specifying the partial intervals
!! @param[in]       n           the length of x
!! @param[in]       k_guess     initial guess (optional)
INTERFACE findinterval
    procedure findinterval_1D, findinterval_2D
END INTERFACE


INTERFACE interpolate_lin
    procedure interpolate1D_lin, interpolate2D_lin
END INTERFACE


CONTAINS

!> @brief (Available via \ref findinterval, not publicly visible)
!! Returns (the left index of) the interval where is placed x0 value,
!! ordering <..)<..)..<..>, extrapolation 0, n
!!
!! This suboroutine finds the index of the interval where 'x0' is placed within the array 'x'
!! It uses the bisection lookup and allows for an initial guess.
!!
!! @param[out]      k1          the seeked index
!! @param[in]       x0          
!! @param[in]       x           array specifying the partial intervals
!! @param[in]       n           length of 'x'
!! @param[in]       k_guess     initial guess
subroutine findinterval_1D(k1,x0,x,k_guess)
    integer, intent(out)                :: k1
    integer                             :: n
	real(8), intent(in)                 :: x0, x(:)
 	integer, optional       :: k_guess

    integer :: k2, length;
   
    n = size(x)

    k2 = n    

    ! check boundary cases
    if (x0 >= x(n)) then
	    if (x0 == x(n)) then
            k1 = n-1
        else
            k1 = n
        endif
        return
    elseif (x0 < x(1)) then
        k1 = 0
        return
    endif

    ! procedure with doubling the interval starting from the initial guess
    ! + eliminating subintervals (Nelder-Mead-ish style)
    if (present(k_guess)) then
        k1 = k_guess
        length = 1
        if ( x0 >= x(k_guess) ) then
            k2 = k_guess + length
            do while ( x(k2) <= x0 )
                k1 = k2                     
                length = 2*length
                k2 = min(k_guess+length,n)
            enddo
        else
            k2 = k1
            k1 = k_guess - length
            do while ( x(k1) > x0 )
                k2 = k1
                length = 2*length
                k1 = max(k_guess-length,1)
            enddo
        endif
        length = k2-k1
    else
        k1 = 1
        length = n
    endif

    ! bisection
    do while (length > 1)
        if ( x0 < x(k1 + (length/2)) ) then
            k2 = k1 + (length/2)
        else
            k1 = k1 + (length/2)
        endif
        length = k2 - k1
    enddo

end subroutine findinterval_1D


!> @brief (Available via \ref findinterval, not publicly visible)
!! Returns (the left indices of) the intervals where ('x0', 'y0')
!! are placed. Ordering <..)<..)..<..>, extrapolation 0, n.
!!
!! \ref findinterval_1D is used in each dimension.
!!
!! @param      kx           [out]          the seeked 'x'-index
!! @param      ky           [out]          the seeked 'y'-index
!! @param      x0           [in]       
!! @param      y0           [in]        
!! @param      x            [in]           array specifying the partial intervals
!! @param      y            [in]           array specifying the partial intervals
!! @param      Nx           [in]           len(x)
!! @param      Ny           [in]           len(y)
!! @param      k_guess      [in]           initial guess (integer)
subroutine findinterval_2D(kx,ky,x0,y0,x,y,kx_guess,ky_guess) ! returns interval where is placed x value, if it is out of the range, 0 is used
!intervals are ordered: <..)<..)<..)...<..>
    integer, intent(out)                :: kx,ky
    integer                             :: Nx,Ny
	real(8), intent(in)                 :: x0, y0, x(:), y(:)
 	integer, intent(in), optional       :: kx_guess, ky_guess
	
    Nx = size(x); Ny = size(y)

	if (present(kx_guess)) then
        call findinterval_1D(kx,x0,x, k_guess=kx_guess)
    else
        call findinterval_1D(kx,x0,x)
    endif

    if (present(ky_guess)) then
        call findinterval_1D(ky,y0,y, k_guess=ky_guess)
    else
        call findinterval_1D(ky,y0,y)
    endif

end subroutine findinterval_2D


!subroutine interpolate(n,x,y,x_grid,y_grid) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
subroutine interpolate1D_decomposed_eq(k,x,fx,xgrid,fxgrid,n,tol) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
	real(8), intent(out)    :: fx
    integer, intent(in)     :: k,n
	real(8), intent(in)     :: x, xgrid(n), fxgrid(n)	
    real(8), optional       :: tol
    real(8), parameter      :: eps_def = EPSILON(1.D0)
    real(8)                 :: eps
    

    if (present(tol)) then
        eps = tol
    else
        eps = eps_def
    endif


    if ( (k>1) .and. ( k<n ) ) then
        if ( abs(xgrid(k)-xgrid(k-1)) < eps  ) then
            fx = fxgrid(k)
            return
        elseif ( abs(xgrid(k)-xgrid(k+1)) < eps  ) then
            fx = fxgrid(k+1)
            return
        else
            fx = fxgrid(k)+(x-xgrid(k))*(fxgrid(k+1)-fxgrid(k))/(xgrid(k+1)-xgrid(k))
            return
        endif
    elseif (k==n) then
        fx = fxgrid(n)
        return
    elseif (k==0) then
        fx = fxgrid(1)
        return
    elseif (k==1) then
        if ( abs(xgrid(2)-xgrid(1)) < eps ) then
            fx = fxgrid(2)
            return
        else
            fx = fxgrid(1)+(x-xgrid(1))*(fxgrid(2)-fxgrid(1))/(xgrid(2)-xgrid(1))
            return
        endif
    endif

	
	! call findinterval(n,x,x_grid,k1,k2);
	! !write(*,*) k1;
	! if ( k1==0 ) then
	! 	y=y_grid(1);
	! elseif ( k2 == 0 ) then
	! 	y=y_grid(n);
	! else
	! 	y=y_grid(k1)+(x-x_grid(k1))*(y_grid(k2)-y_grid(k1))/(x_grid(k2)-x_grid(k1));
	! endif

end subroutine interpolate1D_decomposed_eq

subroutine interpolate1D_lin(x,fx,xgrid,fxgrid,n,k_known,tol) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
	real(8), intent(out)    :: fx
    integer, intent(in)     :: n
	real(8), intent(in)     :: x, xgrid(n), fxgrid(n)	
    real(8), optional       :: tol
    integer, optional       :: k_known

    integer                 :: k1
    
    if (present(k_known)) then
        k1 = k_known
    else
        call findinterval_1D(k1,x,xgrid)
    endif

    if (present(tol)) then
        call interpolate1D_decomposed_eq(k1,x,fx,xgrid,fxgrid,n,tol=tol)
    else
        call interpolate1D_decomposed_eq(k1,x,fx,xgrid,fxgrid,n)
    endif

end subroutine interpolate1D_lin

subroutine interpolate2D_decomposed_eq(kx,ky,x,y,fxy,xgrid,ygrid,fxygrid,Nx,Ny) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
	real(8), intent(out)    :: fxy
    integer, intent(in)     :: kx,ky,Nx,Ny
	real(8), intent(in)     :: x,y,xgrid(Nx),ygrid(Ny),fxygrid(Nx,Ny)

    real(8), parameter          :: eps = EPSILON(1.D0)
    real(8)                     :: fx1,fx2


    ! first interpolate in y, then interpolate in x

    ! the same decision logic in x, encapsulates y
    if ( (kx>1) .and. ( kx< Nx ) ) then
        if ( abs(xgrid(kx)-xgrid(kx-1)) < eps  ) then
            call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(kx,:),  Ny)
            return
        elseif ( abs(xgrid(kx)-xgrid(kx+1)) < eps  ) then
            call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(kx+1,:),  Ny)
            return
        else
            call interpolate1D_decomposed_eq(ky,y,fx1,ygrid,fxygrid(kx,:),  Ny)
            call interpolate1D_decomposed_eq(ky,y,fx2,ygrid,fxygrid(kx+1,:),Ny)
            fxy = fx1+(x-xgrid(kx))*(fx2-fx1)/(xgrid(kx+1)-xgrid(kx))
            return
        endif
    elseif (kx==Nx) then
        call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(Nx,:),  Ny)
        return
    elseif (kx==0) then
        call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(1,:),  Ny)
        return
    elseif (kx==1) then
        if ( abs(xgrid(2)-xgrid(1)) < eps ) then
            call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(2,:),  Ny)
            return
        else
            call interpolate1D_decomposed_eq(ky,y,fx1,ygrid,fxygrid(1,:),Ny)
            call interpolate1D_decomposed_eq(ky,y,fx2,ygrid,fxygrid(2,:),Ny)
            fxy = fx1+(x-xgrid(kx))*(fx2-fx1)/(xgrid(kx+1)-xgrid(kx))
            return
        endif
    endif

end subroutine interpolate2D_decomposed_eq

subroutine interpolate2D_lin(x,y,fxy,xgrid,ygrid,fxygrid,Nx,Ny,kx_known,ky_known) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
	real(8), intent(out)                :: fxy
    integer, intent(in)                 :: Nx,Ny
	real(8), intent(in)                 :: x,y,xgrid(Nx),ygrid(Ny),fxygrid(Nx,Ny)
    integer, intent(in), optional       :: kx_known, ky_known

    integer                             :: kx, ky
    
    if (present(kx_known)) then
        kx = kx_known
    else
        call findinterval_1D(kx,x,xgrid)
    endif

    if (present(ky_known)) then
        ky = ky_known
    else
        call findinterval_1D(ky,y,ygrid)
    endif

    call interpolate2D_decomposed_eq(kx,ky,x,y,fxy,xgrid,ygrid,fxygrid,Nx,Ny)

end subroutine interpolate2D_lin

end module array_helper
