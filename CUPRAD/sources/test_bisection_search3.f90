program test_bisection_search

use array_helper

implicit none

integer, parameter          :: N_elem_max = 2000
integer                     :: k1, k2, k3, k4, k_found, k_found_guess
real(8), allocatable        :: testvals(:), array(:)


logical                     :: works_for_all_guesses, found_correctly, test_passed


! print *, 'array32:', array32
! print *, 'array10:', array10
! print *,

integer clock_start, clock_end, clock_rate

call system_clock(clock_start, clock_rate)

k2=0; k4 = 0
k_found = 0; k_found_guess = 0
test_passed = .true.
do k1 = 2, N_elem_max

    array = (/ (real(k2,8), k2=0, k1-1) /)
    testvals = (/ (-.5d0 + .5*k2, k2=0, 2*k1 + 1) /)

    print *, '---------------------------------------'
    print *, 'testing N_interval=', k1

    found_correctly = .true.
    do k2 = 1, 2*k1 + 1
        call findinterval(k_found,testvals(k2),array)

        ! print *, 'value', testvals(k2) ,'interval', k_found, (k_found == k2/2), k2, k1

        if (k2 /= 2*k1) then
            if (k_found /= k2/2) then
                print *, '!!! FAIL while finding the value', testvals(k2)
                found_correctly = .false.
            endif
        else
            if (k_found /= (k1-1)) then
                print *, '!!! FAIL while finding the value', testvals(k2)
                found_correctly = .false.
            endif
        endif

        works_for_all_guesses = .true.
        do k3 = 1, k1
            call findinterval(k_found_guess,testvals(k2),array,k_guess=k3)
            if (k_found /= k_found_guess)  then
                works_for_all_guesses = .false.
                found_correctly = .false.
                print *, '!!! FAIL for the guess:', k3, 'while finding the value', testvals(k2)
            endif
        enddo
    enddo
    if (found_correctly) then
        print *, 'Works fine'
    else
        print *, '!!! TABLE LOOKUP FAILED, look above for details'
        test_passed = .false.
    endif

enddo

print *, '---------------------------------------'
if (test_passed) then
    print *, 'SUCCESS! Test works for all values and guesses'
else
    print *, '!!! TEST FAILED !!!'
endif
print *, '---------------------------------------'
call system_clock(clock_end)

print *, "total elapsed time: ", real(clock_end - clock_start) / real(clock_rate)

end program test_bisection_search