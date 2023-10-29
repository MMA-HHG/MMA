program test_bisection_search

use array_helper

implicit none

integer, parameter          :: N_arrs = 3, N_elem=6, N1_test = 2*N_elem + 1, N_elem_max = 4
integer                     :: k1, k2, k3, k4, k_found, k_found_guess
real(8), allocatable        :: testvals(:), array(:) ! (N_elem) = (/ (real(k1,8), k1=0, N_elem-1) /)

! real(8), allocatable(:)     :: testvals = (/ (-.5d0 + .5*k1, k1=0, N1_test-1) /)

integer                     :: indexes_noguess(N_arrs,N1_test), indexes_guess(N_arrs,N1_test)
logical                     :: works_for_all_guesses, found_correctly


! print *, 'array32:', array32
! print *, 'array10:', array10
! print *,

k2=0; k4 = 0
k_found = 0; k_found_guess = 0
do k1 = 2, N_elem_max

    array = (/ (real(k2,8), k2=0, k1-1) /)
    testvals = (/ (-.5d0 + .5*k2, k2=0, 2*k1 + 1) /)

    print *, array
    print *, '---'
    print *, testvals

    print *, '---------------------------------------'

    found_correctly = .true.
    do k2 = 1, 2*k1 + 1
        call findinterval_1D(k_found,testvals(k2),array,k1)

        print *, 'value', testvals(k2) ,'interval', k_found, (k_found == k2/2), k2, k1

        if (k2 /= 2*k1) then
            if (k_found == k2/2) then
                print *, 'correct'
            else
                print *, 'fail'
            endif
        else
            if (k_found == (k1-1)) then
                print *, 'correct'
            else
                print *, 'fail'
            endif
        endif

        works_for_all_guesses = .true.
        do k3 = 1, k1
            call findinterval_1D(k_found_guess,testvals(k2),array,k1,k_guess=k3)
            if (k_found /= k_found_guess)  then
                works_for_all_guesses = .false.
                found_correctly = .false.
                print *, 'problem for the guess:', k3
            endif
        enddo
        if (works_for_all_guesses) print *, 'Works for all guesses.'

    enddo

    ! print *, 'arr5', 'value', testvals(k1)
    ! call findinterval_1D(k2,testvals(k1),array5,N_elem)

    ! works_for_all = .true.
    ! do k3 = 1, N_elem-1
    !     print *, '---- next guess ----'
    !     call findinterval_1D(k4,testvals(k1),array5,N_elem,k_guess=k3)
    !     if (k2 /= k4)  then
    !         works_for_all = .false.
    !         print *, 'problem for the guess:', k3
    !     endif
    ! enddo
    ! print *, 'interval', k2
    ! if (works_for_all) print *, 'Works for all guesses.'
    ! print *,

enddo




end program test_bisection_search