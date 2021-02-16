#include "alias.inc"
subroutine test()
  use iso_fortran_env
  use do_math
  implicit none
  real*8 aaa(3,3)

    aaa(1,:) = (/1d0,2d0, 2d0/)
    aaa(2,:) = (/2d0,7d0, 6d0/)
    aaa(3,:) = (/2d0,6d0,10d0/)
    write(6,*)"CCCCC  ", inv(aaa)
stop


endsubroutine

