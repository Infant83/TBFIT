#include "alias.inc"
subroutine test()
  implicit none
  real*8   A(2,2)

  A=9.0d0

  write(6,*)"ZZZ ", A
  stop 
endsubroutine

