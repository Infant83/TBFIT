#include "alias.inc"
subroutine test()
  use print_io
  implicit none
  real*8       cc(4), dd(4), ee(4)
  integer*4    irange(2), ik,is, i, j
  real*8       AA(10)
  character*10 fm
  character*20, external :: int2str
  character*80 fname_header, header

  irange(1) = 1
  irange(2) = 3
  cc(1) = 11.1d0
  cc(2) = 22.1d0
  cc(3) = 33.1d0
  cc(4) = 44.1d0

  open(99,file='test.dat',form='unformatted', status='unknown')
  do i = 1, 1 
  write(99,'(i3, *(F10.5))', advance='no') irange(i), (cc(j),j=1,irange(i))
  write(6,'(i3, *(F10.5))', advance='no') irange(i), (cc(j),j=1,irange(i))
  enddo
  write(99,'(i3, *(F10.5))'              ) irange(2), (cc(j),j=1,irange(2))
  write(6,'(i3, *(F10.5))'               ) irange(2), (cc(j),j=1,irange(2))

  close(99)
  read(99,'(i3,*(F10.5))',advance='no') irange(i), (cc(i),i=1,irange(i)) 

  header = ''
  fname_header = trim('band_structure_TBA'//trim(header))

  stop

endsubroutine
