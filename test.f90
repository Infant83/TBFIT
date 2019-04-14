#include "alias.inc"
subroutine test()
  implicit none
  integer*4    A(10), B(5)
  integer*4    i

  A = 0
  do i = 1,5
    A(i) = i
  enddo
  B=(/1,3,5,7,9/)

  write(6,*)"XXX ", A( (/1,3,5,7,9/) )
  write(6,*)"BBB ", A( B )
  
  stop

endsubroutine
