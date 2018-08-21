module kronecker_prod
 implicit none
 public :: kproduct
 public :: kproduct_r
contains

function kproduct(A,B,row_a,col_a,row_b,col_b) result (AB)
 implicit none
 integer*4, intent(in) ::  row_a,col_a, row_b,col_b
 complex*16 A(row_a,col_a)
 complex*16 B(row_b,col_b)
 complex*16 AB(row_a*row_b,col_a*col_b)
 integer*4  i, j , k , l
 integer*4  m, n , p , q

 do i = 1,row_a
  do j = 1,col_a
   n=(i-1)*row_b + 1
   m=n+row_b -1
   p=(j-1)*col_b + 1
   q=p+col_b -1

   AB(n:m,p:q) = A(i,j)*B
  enddo
 enddo

endfunction

function kproduct_r(A,B,row_a,col_a,row_b,col_b) result (AB)
 implicit none
 integer*4, intent(in) ::  row_a,col_a, row_b,col_b
 real*8     A(row_a,col_a)
 real*8     B(row_b,col_b)
 real*8     AB(row_a*row_b,col_a*col_b)
 integer*4  i, j , k , l
 integer*4  m, n , p , q

 do i = 1,row_a
  do j = 1,col_a
   n=(i-1)*row_b + 1
   m=n+row_b -1
   p=(j-1)*col_b + 1
   q=p+col_b -1
   if( A(i,j) .eq. 0) then
     AB(n:m,p:q) = 0d0
   else
     AB(n:m,p:q) = A(i,j)*B
   endif
  enddo
 enddo

endfunction

end module
