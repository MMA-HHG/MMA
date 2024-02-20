program test_arrays_constructors

implicit none



integer, parameter            :: Ndim = 10
integer                       :: k1
integer, parameter            :: arr_int(Ndim) = [1,2,(k1,k1=3,7),8,9,10]
character(len=255), parameter :: arr_char(Ndim) = [character(len=255) :: "a", &
                                                   'bb', ("ccc", k1 = 3, 7), &
                                                   "dddd", "eeeee", "ffffff"] 

#if INTEL
character(*), parameter          :: arr_char2(Ndim) = (/ "a", &                   ! THIS IS FINE FOR INTEL
                                                      'bb', ("ccc", k1 = 3, 7), &
                                                      "dddd", "eeeee", "ffffff" /)
#endif

! character(*), parameter          :: arr_char2(Ndim) = [character(len=*) :: "a", &                  
!                                                       'bb', ("ccc", k1 = 3, 7), &
!                                                       "dddd", "eeeee", "ffffff" ]


    print *, arr_int
    print *, arr_char
#if INTEL
    print *, arr_char2
#endif


end program test_arrays_constructors