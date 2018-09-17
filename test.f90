subroutine test ()
   implicit none
   integer*4  W(3,3)
   real*8     T(3)
   real*8     R(3)
   real*8     R_(3)

   W(1,:) = (/0, -1, 0/)
   W(2,:) = (/1,  1, 0/)
   W(3,:) = (/0,  0, 1/)

   T(1)   = 0d0
   T(2)   = 1d0/2d0
   T(3)   = 0d0
  
   R(1)   = 1d0/6d0
   R(2)   = 1d0/6d0
   R(3)   = 0d0

   R_ = matmul(W,R) + T
 
   write(6,*)"XXXX", mod(R_ + 1d0,1d0)

stop
return
endsubroutine
