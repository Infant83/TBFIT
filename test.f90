#include "alias.inc"
subroutine test()
  use print_io
  use sorting
  use do_math, only : idx2Di, idx2Dj
  implicit none
  character*2  A
  character*2  B

  B='xx'

  A='ab'

  B=A(2:2)//A(1:1)
   
  write(6,*)"ZZZ ", B

  stop

endsubroutine
