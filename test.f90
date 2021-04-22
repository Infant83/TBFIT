#include "alias.inc"
subroutine test()
! use iso_fortran_env
  use do_math
  implicit none
  integer*4, allocatable::   a(:)
  integer*4   b(3,2 )
  integer*4   na
  integer*8   nb
  integer*4   nn(80000,80000,1)

    write(6,*)"BBB ",  2**31-1  - size(nn,kind=8)

 stop
endsubroutine
