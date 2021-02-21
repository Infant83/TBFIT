#include "alias.inc"
subroutine test()
  use iso_fortran_env
  use do_math
  implicit none
  integer*4, allocatable::   a(:)
  integer*4   b(3,2 )
  integer*4   na
  b = 1 
  na = 3
  b(1,1) = -1 
  b(1,2) = -2

  write(6,*)"ZZ ", b

  write(6,*)"VV ", abs(b)

  write(6,*)"CXX", sum(abs(b))

 stop
endsubroutine


subroutine aa(a, b, na)
    integer*4   a(na)
    integer*4   b(na)

    b = 99

    return
endsubroutine
