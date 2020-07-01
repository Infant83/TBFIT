module kronecker_prod
 use parameters, only : zi, spmat
 implicit none
 public :: kproduct
 public :: kproduct_r
contains

subroutine kproduct_pauli_y_CSR(CSR)
 implicit none
 type(spmat)  ::  CSR
 integer*4        nnz, m
 integer*4        I(CSR%msize*2 + 1)
 integer*4        J(CSR%nnz*2)
 complex*16       H(CSR%nnz*2)

 nnz = CSR%nnz
 m   = CSR%msize

 H(1:nnz) = -zi * CSR%H
 H(nnz+1:nnz*2) = zi * CSR%H 

 I(1:m)     = CSR%I(1:m)
 I(m+1:m*2) = CSR%I(1:m) + nnz
 I(m*2 + 1) = nnz * 2 + 1

 J(1:nnz)       = CSR%J + m
 J(nnz+1:nnz*2) = CSR%J 

 deallocate(CSR%H)
 deallocate(CSR%I)
 deallocate(CSR%J)
   allocate(CSR%H(nnz * 2))
   allocate(CSR%I(m * 2+1))
   allocate(CSR%J(nnz * 2))

 CSR%H = H
 CSR%I = I
 CSR%J = J
 CSR%msize=m * 2
 CSR%nnz  =nnz * 2

return
endsubroutine

subroutine kproduct_pauli_y_sparse(A,IA,JA,na_scr,na,B,IB,JB)
 implicit none
 integer*4, intent(in) :: na_scr, na
 complex*16 A(na_scr)   ! input  sparse matrix array
 complex*16 B(na_scr*2) ! output sparse matrix array
 integer*4  i, j
 integer*4  IA(na+1)
 integer*4  JA(na_scr)
 integer*4  IB(na*2+1)
 integer*4  JB(na_scr*2)
 ! this subroutine performs Kronecker product between 
 ! pauli matrix sigma_z and sparse matrix A and return B as a result.
 ! Array index IA, and JA is also transformed into corresponding IB and JB

 B(1:na_scr) = -zi * A
 B(na_scr+1:na_scr*2) = zi * A 

 IB(1:na) = IA(1:na)
 IB(na+1:na*2) = IA(1:na) + na_scr
 IB(na*2+1) = na_scr * 2 + 1

 JB(1:na_scr) = JA + na
 JB(na_scr+1:na_scr*2) = JA 

endsubroutine

subroutine kproduct_pauli_x_CSR(CSR)
 implicit none
 type(spmat)  ::  CSR
 integer*4        nnz, m
 integer*4        I(CSR%msize*2 + 1)
 integer*4        J(CSR%nnz*2)
 complex*16       H(CSR%nnz*2)



 nnz = CSR%nnz
 m   = CSR%msize

 H(1:nnz) = CSR%H
 H(nnz+1:nnz*2) = CSR%H

 I(1:m)     = CSR%I(1:m)
 I(m+1:m*2) = CSR%I(1:m) + nnz
 I(m*2 + 1) = nnz * 2 + 1

 J(1:nnz)       = CSR%J + m
 J(nnz+1:nnz*2) = CSR%J

 deallocate(CSR%H)
 deallocate(CSR%I)
 deallocate(CSR%J)
   allocate(CSR%H(nnz * 2))
   allocate(CSR%I(m * 2+1))
   allocate(CSR%J(nnz * 2))

 CSR%H = H
 CSR%I = I
 CSR%J = J
 CSR%msize=m * 2
 CSR%nnz  =nnz * 2

return
endsubroutine

subroutine kproduct_pauli_x_sparse(A,IA,JA,na_scr,na,B,IB,JB)
 implicit none
 integer*4, intent(in) :: na_scr, na
 complex*16 A(na_scr)   ! input  sparse matrix array
 complex*16 B(na_scr*2) ! output sparse matrix array
 integer*4  i, j
 integer*4  IA(na+1)
 integer*4  JA(na_scr)
 integer*4  IB(na*2+1)
 integer*4  JB(na_scr*2)
 ! this subroutine performs Kronecker product between 
 ! pauli matrix sigma_z and sparse matrix A and return B as a result.
 ! Array index IA, and JA is also transformed into corresponding IB and JB

 B(1:na_scr) = A
 B(na_scr+1:na_scr*2) = A

 IB(1:na) = IA(1:na)
 IB(na+1:na*2) = IA(1:na) + na_scr
 IB(na*2+1) = na_scr * 2 + 1

 JB(1:na_scr) = JA + na
 JB(na_scr+1:na_scr*2) = JA 

endsubroutine

subroutine kproduct_pauli_z_sparse(A,IA,JA,na_scr,na,B,IB,JB)
 implicit none
 integer*4, intent(in) :: na_scr, na
 complex*16 A(na_scr)   ! input  sparse matrix array
 complex*16 B(na_scr*2) ! output sparse matrix array
 integer*4  i, j
 integer*4  IA(na+1)
 integer*4  JA(na_scr)
 integer*4  IB(na*2+1)
 integer*4  JB(na_scr*2)
 ! this subroutine performs Kronecker product between 
 ! pauli matrix sigma_z and sparse matrix A and return B as a result.
 ! Array index IA, and JA is also transformed into corresponding IB and JB

 B(1:na_scr) = A
 B(na_scr+1:na_scr*2) = A * (-1d0)

 IB(1:na) = IA(1:na)
 IB(na+1:na*2) = IA(1:na) + na_scr
 IB(na*2+1) = na_scr * 2 + 1

 JB(1:na_scr) = JA
 JB(na_scr+1:na_scr*2) = JA + na

endsubroutine

subroutine kproduct_r_pauli_0_sparse(A,IA,JA,na_scr,na,B,IB,JB)
 implicit none
 integer*4, intent(in) :: na_scr, na
 real*8     A(na_scr)   ! input  sparse matrix array
 real*8     B(na_scr*2) ! output sparse matrix array
 integer*4  i, j
 integer*4  IA(na+1)
 integer*4  JA(na_scr)
 integer*4  IB(na*2+1)
 integer*4  JB(na_scr*2)
 ! this subroutine performs Kronecker product between 
 ! identity matrix I and sparse matrix A and return B as a result.
 ! Array index IA, and JA is also transformed into corresponding IB and JB

 B(1:na_scr) = A
 B(na_scr+1:na_scr*2) = A

 IB(1:na) = IA(1:na)
 IB(na+1:na*2) = IA(1:na) + na_scr
 IB(na*2+1) = na_scr * 2 + 1

 JB(1:na_scr) = JA
 JB(na_scr+1:na_scr*2) = JA + na

endsubroutine

subroutine kproduct_pauli_z_CSR(CSR)
 implicit none
 type(spmat)  ::  CSR
 integer*4        nnz, m
 integer*4        I(CSR%msize*2 + 1)
 integer*4        J(CSR%nnz*2)
 complex*16       H(CSR%nnz*2)

 nnz = CSR%nnz
 m   = CSR%msize

 H(1:nnz) = CSR%H 
 H(nnz+1:nnz*2) = -1 * CSR%H

 I(1:m)     = CSR%I(1:m)
 I(m+1:m*2) = CSR%I(1:m) + nnz
 I(m*2 + 1) = nnz * 2 + 1

 J(1:nnz)       = CSR%J
 J(nnz+1:nnz*2) = CSR%J + m

 deallocate(CSR%H)
 deallocate(CSR%I)
 deallocate(CSR%J)
   allocate(CSR%H(nnz * 2))
   allocate(CSR%I(m * 2+1))
   allocate(CSR%J(nnz * 2))

 CSR%H = H
 CSR%I = I
 CSR%J = J
 CSR%msize=m * 2
 CSR%nnz  =nnz * 2

return
endsubroutine

subroutine kproduct_pauli_0_CSR(CSR)
 implicit none
 type(spmat)  ::  CSR
 integer*4        nnz, m
 integer*4        I(CSR%msize*2 + 1)
 integer*4        J(CSR%nnz*2)
 complex*16       H(CSR%nnz*2)

 nnz = CSR%nnz
 m   = CSR%msize

 H(1:nnz) = CSR%H
 H(nnz+1:nnz*2) = CSR%H

 I(1:m)     = CSR%I(1:m)
 I(m+1:m*2) = CSR%I(1:m) + nnz
 I(m*2 + 1) = nnz * 2 + 1

 J(1:nnz)       = CSR%J 
 J(nnz+1:nnz*2) = CSR%J + m

 deallocate(CSR%H)
 deallocate(CSR%I)
 deallocate(CSR%J)
   allocate(CSR%H(nnz * 2))
   allocate(CSR%I(m * 2+1))
   allocate(CSR%J(nnz * 2))

 CSR%H = H
 CSR%I = I
 CSR%J = J
 CSR%msize=m * 2
 CSR%nnz  =nnz * 2

return
endsubroutine
subroutine kproduct_pauli_0_sparse(A,IA,JA,na_scr,na,B,IB,JB)
 implicit none
 integer*4, intent(in) :: na_scr, na
 complex*16 A(na_scr)   ! input  sparse matrix array
 complex*16 B(na_scr*2) ! output sparse matrix array
 integer*4  i, j
 integer*4  IA(na+1)
 integer*4  JA(na_scr)
 integer*4  IB(na*2+1)
 integer*4  JB(na_scr*2)
 ! this subroutine performs Kronecker product between 
 ! identity matrix I and sparse matrix A and return B as a result.
 ! Array index IA, and JA is also transformed into corresponding IB and JB

 B(1:na_scr) = A
 B(na_scr+1:na_scr*2) = A

 IB(1:na) = IA(1:na)
 IB(na+1:na*2) = IA(1:na) + na_scr
 IB(na*2+1) = na_scr * 2 + 1

 JB(1:na_scr) = JA
 JB(na_scr+1:na_scr*2) = JA + na

endsubroutine

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
