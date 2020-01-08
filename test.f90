#include "alias.inc"
subroutine test()
  implicit none
  real*8       cc(4), dd(4), ee(4)
  integer*4    irange(2)
  irange(1:2) = (/1,3/)

  dd = (/1d0,2d0,3d0,4d0/)
  write(6,*)"EE1", dd(:)
  write(6,*)"EE2", dd(irange)
  stop

endsubroutine
