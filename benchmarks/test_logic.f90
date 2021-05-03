PROGRAM test_modules

  IMPLICIT NONE
  
  integer, parameter :: Np = 100;
  integer :: k1,k2
  integer :: iarray4(4) = (/1, 2, 3, 4/)
  logical :: dumlog

  character(*), parameter :: sarray6(6) = (/"a", "ab", "abc", "abcd", "abcde", "abcdf"/)


  print *,"Test any mask on integers"

  
  print *, .true.

  dumlog = ANY(iarray4 == 1)

  print *,"arr ni 1", dumlog

  dumlog = ANY(1 == iarray4 )

  print *,"1 in arr", dumlog

  dumlog = ANY(iarray4 == 5)

  print *,"arr ni 5", dumlog

  dumlog = ANY(5 == iarray4 )

  print *,"5 in arr", dumlog  

  dumlog = ANY( "abc" == (/"a", "ab", "abc", "abcd", "abcde", "abcdf"/) )

  print *,"test array ", dumlog  

  dumlog = ANY( "aaa" == (/"a", "ab", "abcd", "abcde", "abcdef"/) )

  print *,"test array ", dumlog  

  dumlog = ANY( "abc" == sarray6 )

  print *,"test array ", dumlog  

  dumlog = ANY( "aaa" == sarray6 )

  print *,"test array ", dumlog  

END PROGRAM test_modules
