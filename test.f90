#include "alias.inc"
subroutine test()
  implicit none
  integer*4    A
  integer*8    B
  integer*8    C 

  A=200000000
  B=300000000

  write(6,*)"XX ", A, B, A * B

  C = int8(A) * int8(A)

  call ww(A*B, A )
 
  stop

endsubroutine

subroutine ww(AA, BB)
  integer*4 AA
  integer*8 BB

  write(6,*)"XXXX ", AA, BB * BB

endsubroutine
