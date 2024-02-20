program test_bisection_search

use array_helper

implicit none

integer, parameter  :: N1_test = 12, N_arrs = 3
integer             :: k1, k2
real(8), parameter  :: array32(32) = (/ (real(k1,8), k1=0,31) /), &
                       array10(10) = (/ (real(k1,8), k1=0,9) /), &
                       array13(13) = (/ (real(k1,8), k1=0,12) /)

real(8)              :: testvals(N1_test) = (/ 0.5d0, 2.0d0, 5.0d0, 7.5d0, 0.0d0, 1.0d0, 2.0d0, 1.99999d0, 8.5d0, 6.6d0, 7.1d0, real(9,8) /)

integer              :: indexes_noguess(N_arrs,N1_test), indexes_guess(N_arrs,N1_test)


! print *, 'array32:', array32
! print *, 'array10:', array10
! print *,

k2=0
do k1 = 1, N1_test

    print *, 'arr32', 'value', testvals(k1)
    call findinterval_1D(k2,testvals(k1),array32,32)
    print *, 'interval', k2
    indexes_noguess(1,k1) = k2

    print *,

    print *, 'arr10', 'value', testvals(k1)
    call findinterval_1D(k2,testvals(k1),array10,10)
    print *, 'interval', k2 
    indexes_noguess(2,k1) = k2

    print *,

    print *, 'arr13', 'value', testvals(k1)
    call findinterval_1D(k2,testvals(k1),array13,13)
    print *, 'interval', k2 
    indexes_noguess(3,k1) = k2

    print *,

enddo



print *, 'with guess'
print *, '---------------------------------------------------'

k2=0
do k1 = 1, N1_test

    print *, 'arr32', 'value', testvals(k1)
    call findinterval_1D(k2,testvals(k1),array32,32,k_guess=5)
    print *, 'interval', k2 
    indexes_guess(1,k1) = k2
    print *, indexes_noguess(1,k1) == indexes_guess(1,k1)

    print *,

    print *, 'arr10', 'value', testvals(k1)
    call findinterval_1D(k2,testvals(k1),array10,10,k_guess=5)
    print *, 'interval', k2 
    indexes_guess(2,k1) = k2
    print *, indexes_noguess(2,k1) == indexes_guess(2,k1)

    print *,

    print *, 'arr13', 'value', testvals(k1)
    call findinterval_1D(k2,testvals(k1),array13,13,k_guess=5)
    print *, 'interval', k2 
    indexes_guess(3,k1) = k2
    print *, indexes_noguess(3,k1) == indexes_guess(3,k1)

    print *,

enddo

print *, indexes_noguess == indexes_guess

end program test_bisection_search