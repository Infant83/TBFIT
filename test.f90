#include "alias.inc"
subroutine test()
  implicit none
  integer*4    A
  integer*8    B
  integer*8    C
  integer*4    i, irecl, irecl_temp
  real*8       CC(3), DD(3)

  irecl = 3 + 1
  CC = (/1d0,2d0,3d0/)


  open(10,file='test.dat',form='unformatted',access='direct', status='unknown', recl=irecl)

  write(6,*)"xxx ", CC

  write(10, rec=1) irecl
  do i = 2, 4
    write(10, rec=i) CC(i)
  enddo

  close(10)

  open(100,file='test.dat',form='unformatted',access='direct',status='unknown',recl=1)
  read(100,rec=1) irecl_temp
  write(6,*)"XXX ", irecl_temp
  close(100)

  open(100,file='test.dat',form='unformatted',access='direct',status='unknown',recl=irecl_temp)
  
  read(100,rec=2) irecl_temp
  

  stop

endsubroutine
