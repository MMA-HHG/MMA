! It was developed by Jan Vabek
module array_helper

implicit none
private
public  :: findinterval, interpolate_lin
public  :: interpolate1D_decomposed_eq
public  :: interpolate2D_decomposed_eq

INTERFACE findinterval
    procedure findinterval_1D, findinterval_2D
END INTERFACE

INTERFACE interpolate_lin
    procedure interpolate1D_lin, interpolate2D_lin
END INTERFACE


CONTAINS

subroutine findinterval_1D(k1,x0,x,n,k_tip) ! returns interval where is placed x value, if it is out of the range, 0 is used
!intervals are ordered: <..)<..)<..)...<..>
    integer, intent(out)                :: k1
    integer, intent(in)                 :: n
	real(8), intent(in)                 :: x0, x(n)
 	integer,optional       :: k_tip

    integer :: k2, length;
   
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

    ! procedure with doubling the interval using initial tip
    if (present(k_tip)) then
        k1 = k_tip
        length = 1
        if ( x0 >= x(k_tip) ) then
            k2 = k_tip + length
            do while ( x(k2) <= x0 )
                k1 = k2
                length = 2*length
                k2 = min(k_tip+length,n)
            enddo
        else
            k2 = k1
            k1 = k_tip - length
            do while ( x(k1) > x0 )
                k2 = k1
                length = 2*length
                k1 = max(k_tip-length,1)
            enddo
        endif
        length = k2-k1
    else
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


subroutine findinterval_2D(kx,ky,x0,y0,x,y,Nx,Ny,kx_tip,ky_tip) ! returns interval where is placed x value, if it is out of the range, 0 is used
!intervals are ordered: <..)<..)<..)...<..>
    integer, intent(out)                :: kx,ky
    integer, intent(in)                 :: Nx,Ny
	real(8), intent(in)                 :: x0, y0, x(Nx), y(Ny)
 	integer, intent(in), optional       :: kx_tip, ky_tip
	
	if (present(kx_tip)) then
        call findinterval_1D(kx,x0,x,Nx, k_tip=kx_tip)
    else
        call findinterval_1D(kx,x0,x,Nx)
    endif

    if (present(ky_tip)) then
        call findinterval_1D(ky,y0,y,Ny, k_tip=ky_tip)
    else
        call findinterval_1D(ky,y0,y,Ny)
    endif

end subroutine findinterval_2D


!subroutine interpolate(n,x,y,x_grid,y_grid) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
subroutine interpolate1D_decomposed_eq(k,x,fx,xgrid,fxgrid,n,tol) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
	real(8), intent(out)    :: fx
    integer, intent(in)     :: k,n
	real(8), intent(in)     :: x, xgrid(n), fxgrid(n)	
    real(8), optional       :: tol
    real(8), parameter      :: eps_def = EPSILON(1.D0) ! definition here invokes save attribute
    real(8)                 :: eps
    

    if (present(tol)) then
        eps = tol
    else
        eps = eps_def
    endif


    if ( (k>1) .and. ( k<n ) ) then
        if ( abs(xgrid(k)-xgrid(k-1)) < eps  ) then
            fx = fxgrid(k-1)
            return
        elseif ( abs(xgrid(k)-xgrid(k+1)) < eps  ) then
            fx = fxgrid(k)
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
            fx = fxgrid(1)
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

subroutine interpolate1D_lin(x,fx,xgrid,fxgrid,n,tol,k_tip) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
	real(8), intent(out)    :: fx
    integer, intent(in)     :: n
	real(8), intent(in)     :: x, xgrid(n), fxgrid(n)	
    real(8), optional       :: tol
    integer, optional       :: k_tip

    integer                 :: k1
    
    if (present(k_tip)) then
        call findinterval_1D(k1,x,xgrid,n,k_tip=k_tip)
    else
        call findinterval_1D(k1,x,xgrid,n)
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
            call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(kx-1,:),  Ny)
            return
        elseif ( abs(xgrid(kx)-xgrid(kx+1)) < eps  ) then
            call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(kx,:),  Ny)
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
            call interpolate1D_decomposed_eq(ky,y,fxy,ygrid,fxygrid(1,:),  Ny)
            return
        else
            call interpolate1D_decomposed_eq(ky,y,fx1,ygrid,fxygrid(1,:),Ny)
            call interpolate1D_decomposed_eq(ky,y,fx2,ygrid,fxygrid(2,:),Ny)
            fxy = fx1+(x-xgrid(kx))*(fx2-fx1)/(xgrid(kx+1)-xgrid(kx))
            return
        endif
    endif

end subroutine interpolate2D_decomposed_eq

subroutine interpolate2D_lin(x,y,fxy,xgrid,ygrid,fxygrid,Nx,Ny,kx_tip,ky_tip) !inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
	real(8), intent(out)                :: fxy
    integer, intent(in)                 :: Nx,Ny
	real(8), intent(in)                 :: x,y,xgrid(Nx),ygrid(Ny),fxygrid(Nx,Ny)
    integer, intent(in), optional       :: kx_tip, ky_tip

    integer                 :: k1, k2
    
    if (present(kx_tip) .and. present(ky_tip)) then
        call findinterval_2D(k1,k2,x,y,xgrid,ygrid,Nx,Ny,kx_tip=kx_tip,ky_tip=ky_tip)
    elseif (present(kx_tip) .and. .not.(present(ky_tip))) then
        call findinterval_2D(k1,k2,x,y,xgrid,ygrid,Nx,Ny,kx_tip=kx_tip)
    elseif (present(ky_tip)) then
        call findinterval_2D(k1,k2,x,y,xgrid,ygrid,Nx,Ny,ky_tip=ky_tip)
    else
        call findinterval_2D(k1,k2,x,y,xgrid,ygrid,Nx,Ny)
    endif

    call interpolate2D_decomposed_eq(k1,k2,x,y,fxy,xgrid,ygrid,fxygrid,Nx,Ny)

end subroutine interpolate2D_lin

end module array_helper
