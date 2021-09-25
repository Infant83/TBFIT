#include "alias.inc"
subroutine test()
  use random_mod
  implicit none
    integer i
    real*8 r
    real*8 a(5)
    real*8 b(5)

    a=(/1.d0,2.d0,3.d0,4.d0,5.d0/)
    b=0d0
    
    b(4:) = a(4:)
    write(6,*)"ZXX ", b
 stop
endsubroutine
