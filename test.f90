#include "alias.inc"
subroutine test()
  use random_mod
  implicit none
    integer i
    real*8 r
    real*8 pi , a, b
    character*20 temp

    open(11,file='foo.txt',status='unknown')
    read(11,*) a
    read(11,*) b

    close(11)

    write(6,'(A,2F12.6)')"EEE ", a, b

 stop
end subroutine
