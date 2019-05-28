#include "alias.inc"
subroutine test()
  implicit none
  integer*4    A
  integer*8    B
  integer*8    C
  integer*4    i, irecl, irecl_temp
  real*8       CC(3), DD(3)
  real*4       EE(3), FF(3)
  complex*16   GG(3,2)
  complex*16   KK(3,2)
  complex*8    KK_(3,2)

  irecl = 3 + 1
  CC = (/.0001d0,2d0,3d0/) / 5.3234d0
  EE = real(CC,4)
  GG(1:3,1) = (/(1d0,3d0), (5d0,3d0), (.0008d0,9d0)/)
  GG(1:3,2) = GG(1:3,1) * 2.5d0

  open(10,file='test.dat',form='unformatted', status='unknown')

  write(6,*)"xxx ", GG(:,2)
  write(10) (cmplx((/GG(1,i),GG(3,i)/), kind=4),i=1,2)
  close(10)

  open(10, file='test.dat', form='unformatted', status='old')
  read(10) (KK_(1:2,i),i=1,2)
  KK = 0d0
  KK = KK_
  write(6,*)"GGG", KK(1:2,2), KK_(1:2,2), sizeof(KK), sizeof(KK_)

  close(10)

! write(6,*)"XXXX", cmplx(GG,kind=4)
! write(10, rec=1) irecl
! do i = 2, 4
!   write(10, rec=i) CC(i)
! enddo

! close(10)

! open(100,file='test.dat',form='unformatted',access='direct',status='unknown',recl=1)
! read(100,rec=1) irecl_temp
! write(6,*)"XXX ", irecl_temp
! close(100)

! open(100,file='test.dat',form='unformatted',access='direct',status='unknown',recl=irecl_temp)
! 
! read(100,rec=2) irecl_temp
  

  stop

endsubroutine
