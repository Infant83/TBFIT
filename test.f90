#include "alias.inc"
subroutine test()
  use random_mod
  implicit none
    integer i
    real*8 r
    real*8 pi 
    character*20 temp
    call execute_command_line('grep VERSION= makefile > _foo')

    open(11,file='_foo',status='unknown')
    read(11,*) temp

    close(11)

    write(6,*)"EEE ", trim(temp)

 stop
end subroutine
