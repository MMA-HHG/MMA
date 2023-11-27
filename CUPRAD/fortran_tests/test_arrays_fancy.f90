module procedures
implicit none
contains
    subroutine determine_size(vec)
        integer     :: vec(:)

        print *, size(vec)

    end subroutine
end module procedures


program test_arrays_fancy
use procedures

implicit none



integer                                       :: k1
integer, parameter                            :: Nx = 10, Ny = 10 
integer, dimension(Nx,Ny), parameter          :: arr = reshape(  (/ ( k1 , k1=0, Nx*Ny-1) /) ,(/Nx,Ny/)) 



    print *, arr(1,:)
    print *, 'size', size(arr(1,:))
    print *, 'shape', shape(arr(1,:))
    print *, '---'


    print *, arr(:,1) 
    print *, 'size', size(arr(:,1))
    print *, 'shape', shape(arr(:,1))
    print *, '---'


    print *, arr(1:2,:) 
    print *, 'size', size(arr(1:2,:))
    print *, 'shape', shape(arr(1:2,:))
    print *, '---'


    print *, arr(1,3:6)
    print *, 'size', size(arr(1,3:6))
    print *, 'shape', shape(arr(1,3:6))
    print *, '---'

    print *, 'subrout'
    call determine_size(arr(1,3:6))
   

end program test_arrays_fancy