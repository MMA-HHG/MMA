!> @brief This module contains the lookup subroutine \ref findinterval and
!! the linear interpolation subroutine \ref interpolate_lin, the procedures are
!! available up to 2 dimensions. Optionally, bookkeeping is possible (see
!! their descriptions for details).
!! 
!! Precisely, the subroutines are available via the aforementioned generic 
!! interfaces. For \ref interpolate_lin, the specific subroutines \ref
!! interpolate1D_lin and \ref interpolate2D_lin are available as well.
!! Others suppplementary procedures are encapsulatd within this module
!! and not publicly visible.
!!
!! @author Jan Vábek
module array_helper

implicit none
private
public  :: findinterval, interpolate_lin
! public  :: interpolate2D_lin, interpolate1D_lin

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


!> @brief Interpolates a functions based on tabulated values. (Available in 1D & 2D)
!!
!! This subroutines returns the interpolated value 'fx=f(x)' ('fxy=f(x,y)' in 2D, respectively). It uses linear
!! interpolation (bilinear for 2D). Extrapolation is by the boundary values. It supports discontinuous 'xgrid',
!! the interpolation is then by the right value.* \n
!! Note: The interval where 'x' is placed is find within the procedure. This interval can be specified by 'kx_known'.
!! There's no (!) further check whether it's correct. \n
!! Note 2: See the source to reach the private functions of the module encapsulated by this interval.
!!
!! *Example: interpolation using xgrid \f$(0,1,1,2,3,3,4)\f$ and tabulated values \f$(0,0,1,2,3,0,0)\f$ provides \f$f(1)=1\f$, \f$f(3)=0\f$.
!!
!! @param[in]       x           ('x','y' in 2D)           
!! @param[out]      fx          ('fxy' in 2D)
!! @param[in]       xgrid       array for tabulated values ( ... , 'ygrid' in 2D)
!! @param[in]       fxgrid      tabulated values ('fxygrid' in 2D)
!! @param[in]       Nx          the length of xgrid ( ... , 'Ny' in 2D)
!! @param[in]       kx_known    left-boundary of the interval, where we interpolate (optional)
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


!> @cond INCLUDE_ARRAY_HELPER_INTERNALS
! 1D interpolation, the table lookup performed before, see 'interpolate1D_lin'
subroutine interpolate1D_decomposed_eq(k,x,fx,xgrid,fxgrid,n,tol)
	real(8), intent(out)    :: fx
    integer, intent(in)     :: k,n
	real(8), intent(in)     :: x, xgrid(n), fxgrid(n)	
    real(8), optional       :: tol
    real(8), parameter      :: eps_def = EPSILON(1.D0)
    real(8)                 :: eps
    
    ! tolerence used for neigbouring grid points to be equal
    if (present(tol)) then
        eps = tol
    else
        eps = eps_def
    endif


    if ( (k>1) .and. (k<n) ) then
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

end subroutine interpolate1D_decomposed_eq

! 1D-linear interpolation
subroutine interpolate1D_lin(x,fx,xgrid,fxgrid,n,k_known,tol)
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

! 2D-bilinear interpolation, the table lookup performed before, see 'interpolate2D_lin'
subroutine interpolate2D_decomposed_eq(kx,ky,x,y,fxy,xgrid,ygrid,fxygrid,Nx,Ny)
	real(8), intent(out)    :: fxy
    integer, intent(in)     :: kx,ky,Nx,Ny
	real(8), intent(in)     :: x,y,xgrid(Nx),ygrid(Ny),fxygrid(Nx,Ny)

    real(8), parameter          :: eps = EPSILON(1.D0)
    real(8)                     :: fx1,fx2


    ! first interpolate in y, then interpolate in x
    if ( (kx>1) .and. (kx< Nx) ) then
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

! 2D-bilinear interpolation
subroutine interpolate2D_lin(x,y,fxy,xgrid,ygrid,fxygrid,Nx,Ny,kx_known,ky_known)
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
!> @endcond 


end module array_helper
