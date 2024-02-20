program test_bisection_search

use array_helper

implicit none

integer, parameter  :: N_arrs = 3, N_elem=6, N1_test = 2*N_elem + 1 ! 11
integer             :: k1, k2, k3, k4
real(8), parameter  :: array5(N_elem) = (/ (real(k1,8), k1=0, N_elem-1) /)

real(8)              :: testvals(N1_test) = (/ (-.5d0 + .5*k1, k1=0, N1_test-1) /)

integer              :: indexes_noguess(N_arrs,N1_test), indexes_guess(N_arrs,N1_test)
logical              :: works_for_all


! print *, 'array32:', array32
! print *, 'array10:', array10
! print *,

k2=0; k4 = 0
do k1 = 1, N1_test

    print *, '---------------------------------------'
    print *, 'arr5', 'value', testvals(k1)
    call findinterval_1D(k2,testvals(k1),array5,N_elem)

    works_for_all = .true.
    do k3 = 1, N_elem-1
        print *, '---- next guess ----'
        call findinterval_1D(k4,testvals(k1),array5,N_elem,k_guess=k3)
        if (k2 /= k4)  then
            works_for_all = .false.
            print *, 'problem for the guess:', k3
        endif
    enddo
    print *, 'interval', k2
    if (works_for_all) print *, 'Works for all guesses.'
    print *,

enddo




end program test_bisection_search