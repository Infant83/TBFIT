#include "alias.inc"
module sparse_tool
   use do_math
   use mpi_setup
   use parameters, only: spmat
#ifdef MKL_SPARSE
   use ISO_C_BINDING
   use MKL_SPBLAS

   interface
   function MKL_SPARSE_Z_EXPORT_CSR_CF(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) &
              bind(C, name='MKL_SPARSE_Z_EXPORT_CSR')
              use, intrinsic :: iso_c_binding , only : c_int, c_double_complex, c_ptr
              import SPARSE_MATRIX_T
              type(SPARSE_MATRIX_T) , intent(in) :: source
              integer(c_int), intent(inout) :: indexing
              integer       , intent(inout) :: rows
              integer       , intent(inout) :: cols
              type(c_ptr)   , intent(inout) :: rows_start
              type(c_ptr)   , intent(inout) :: rows_end
              type(c_ptr)   , intent(inout) :: col_indx
              type(c_ptr)   , intent(inout) :: values
              integer(c_int) MKL_SPARSE_Z_EXPORT_CSR_CF
   endfunction 
   endinterface 
#endif

contains

#ifdef MKL_SPARSE
  subroutine sparse_create_coo_handle(S, COO) ! create sparse handle for COO format
    implicit none
    type(SPARSE_MATRIX_T) :: S
    type(spmat)           :: COO
    integer*4                indx
    integer*4                m, nnz
    integer*4                istat

    indx = SPARSE_INDEX_BASE_ONE
    m    = COO%msize
    nnz  = COO%nnz

    ! create COO handle
    istat = MKL_SPARSE_Z_CREATE_COO(S, indx, m, m, nnz, COO%I, COO%J, COO%H)
    call sparse_error_report('MKL_SPARSE_Z_CREATE_COO: S sparse_create_coo_handle ', istat)

    return
  endsubroutine
  subroutine sparse_create_csr_handle(S, CSR) ! create sparse handle for CSR format
    implicit none 
    type(SPARSE_MATRIX_T) :: S
    type(spmat)           :: CSR
    integer*4                indx
    integer*4                m, nnz
    integer*4                ierr, istat

    ! Only accpet for square matrix
    indx = SPARSE_INDEX_BASE_ONE 
    m    = CSR%msize
    nnz  = CSR%nnz

    ! create 3-array CSR handle
    istat = MKL_SPARSE_Z_CREATE_CSR(S, indx, m, m, CSR%I, CSR%I(2), CSR%J, CSR%H)
    call sparse_error_report('MKL_SPARSE_Z_CREATE_CSR: S in sparse_create_csr_handle ', istat)

    return
  endsubroutine

  subroutine sparse_export_csr(S, CSR) ! create CSR from sparse handle
    use parameters, only : t1, t0
    use time
    implicit none
    type(SPARSE_MATRIX_T) :: S, SS
    type(spmat)           :: CSR
    type(c_ptr)           :: cptrB, cptrE, cptrJ, cptrH
    integer*4, pointer    :: ptrB(:), ptrE(:), ptrJ(:)
    complex*16,pointer    :: ptrH(:)
    integer*4                indx
    integer*4                m, nnz
    integer*4                ierr, istat

    ! Only accpet for square matrix
    indx = SPARSE_INDEX_BASE_ONE

    if(allocated(CSR%H)) deallocate(CSR%H)
    if(allocated(CSR%I)) deallocate(CSR%I)
    if(allocated(CSR%J)) deallocate(CSR%J)

    ! generate 4-variant CSR pointer
    istat = MKL_SPARSE_CONVERT_CSR(S,SPARSE_OPERATION_NON_TRANSPOSE,SS)
    call sparse_error_report('MKL_SPARSE_CONVERT_CSR: S->SS in sparse_export_csr ', istat)

    istat = MKL_SPARSE_Z_EXPORT_CSR_CF(SS,indx,m,m,    cptrB,    cptrE,    cptrJ,    cptrH)
    call sparse_error_report('MKL_SPARSE_Z_EXPORT_CSR: SS in sparse_export_csr ', istat)

    call c_f_pointer(    cptrB,     ptrB, [m])
    call c_f_pointer(    cptrE,     ptrE, [m])
    nnz =     ptrE(m) - 1
    CSR%nnz = nnz
    CSR%msize = m
    call c_f_pointer(    cptrJ,     ptrJ, [nnz])
    call c_f_pointer(    cptrH,     ptrH, [nnz])
 
    allocate(CSR%H(nnz))
    allocate(CSR%I(m+1))
    allocate(CSR%J(nnz))
    CSR%H(1:nnz) =     ptrH(1:nnz)
    CSR%J(1:nnz) =     ptrJ
    CSR%I(1:m)   =     ptrB(1:m)
    CSR%I(m+1)   = nnz + 1

    istat = MKL_SPARSE_DESTROY(SS)
    call sparse_error_report('MKL_SPARSE_DESTROY: SS in sparse_export_csr ', istat)
    istat = MKL_SPARSE_DESTROY(S)
    call sparse_error_report('MKL_SPARSE_DESTROY: S in sparse_export_csr ', istat)

    return
  endsubroutine

  subroutine sparse_convert_coo_csr(COO, CSR)
    implicit none
    type(SPARSE_MATRIX_T) :: S_COO
    type(spmat)           ::   COO,   CSR
    integer*4                indx, m, nnz
    integer*4                istat

    call sparse_create_coo_handle(S_COO, COO)
    call sparse_export_csr(S_COO, CSR)

    return
  endsubroutine

  subroutine sparse_error_report(func_name, istat)
    implicit none
    integer*4    istat
    character(*) func_name 

    if(istat .ne. SPARSE_STATUS_SUCCESS) then
      write(6,'(2A)')  '    !WARN! ERROR OCCURED IN MKL_SPARSE ROUTINES: ', trim(func_name)
      write(6,'(A,I0)')'           ISTAT = ', istat 
!     if(istat .gt. 1) stop
      stop
    endif

    return
  endsubroutine

#endif

subroutine csr_dns(CSR, H)
  use parameters, only : spmat
  implicit none
  type(spmat) :: CSR
  integer*4  i, j, k, msize, ierr
  complex*16 H(CSR%msize, CSR%msize)

  msize = CSR%msize

  ierr = 0
  H = 0.0D+00
  do i = 1, msize
    do k = CSR%I(i), CSR%I(i+1)-1
      j = CSR%J(k)
      if ( msize < j ) then
        ierr = i
        return
      end if
      H(i,j) = CSR%H(k)
    end do
  end do

  return
endsubroutine

subroutine print_sparse(CSR, sname)
   implicit none
   type(spmat)  :: CSR
   character(*), optional, intent(in) :: sname

   if(present(sname)) then
     write(6,'(2A)')' Sparse Matrix: ', trim(sname)
   endif

   write(6,'(A,*(F4.0,"+i",F4.0))')'VALUES=',CSR%H
   write(6,'(A,*(I4))')'I=',CSR%I
   write(6,'(A,*(I4))')'J=',CSR%J
   write(6,'(A,I4)')'NNZ=',CSR%I(size(CSR%I))-1

return
endsubroutine
subroutine print_sparse_3(A, I, J, nnz, msize, sname)
   implicit none
   integer*4    nnz !number of non-zero elements
   integer*4    msize ! matrix size
   complex*16   A(nnz) ! non-zero element arrays
   integer*4    I(msize+1) ! column index of first non-zero element appears
                           ! msize+1 th element represents nnz
   integer*4    J(nnz) ! column index for each non-zero elements
   character(*), optional, intent(in) :: sname

   if(present(sname)) then
     write(6,'(2A)')' Sparse Matrix: ', trim(sname)
   endif

   write(6,'(A,*(F3.0,"+i",F3.0))')'VALUES=',A
   write(6,'(A,*(I3))')'I=',I
   write(6,'(A,*(I3))')'J=',J
   write(6,'(A,I3)')'NNZ=',I(size(I))-1

return
endsubroutine

subroutine print_sparse_4(A, B, E, J, nnz, msize, sname)
   implicit none
   integer*4    nnz
   integer*4    msize
   complex*16   A(nnz)
   integer*4    B(msize)
   integer*4    E(msize)
   integer*4    J(nnz)   
   character(*), optional, intent(in) :: sname

   if(present(sname)) then
     write(6,'(2A)')' Sparse Matrix: ', trim(sname)
   endif

   write(6,'(A,*(F3.0,"+i",F3.0))')'VALUES=',A
   write(6,'(A,*(I3))')'B=',B
   write(6,'(A,*(I3))')'E=',E
   write(6,'(A,*(I3))')'J=',J
   write(6,'(A,I3)')'NNZ=',E(size(E))-1

return
endsubroutine

! "dnscsr, csrdns, and aplb"
! has been adopted from the SPARSEKIT package 
! written by Youcef Saad (University of Minnesota, saad@cs.umn.edu).
! https://www-users.cs.umn.edu/~saad/software/SPARSKIT/
subroutine dnscsr_r ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )
!*****************************************************************************80
!! DNSCSR converts Dense to Compressed Row Sparse format.
!  Discussion:
!    This routine converts a densely stored matrix into a row orientied
!    compactly sparse matrix.  It is the reverse of CSRDNS.
!    This routine does not check whether an element is small.  It considers 
!    that A(I,J) is zero only if it is exactly equal to zero.
!  Modified: 07 January 2004
!  Author: Youcef Saad
!  Parameters:
!    Input, integer*4 NROW, the row dimension of the matrix.
!    Input, integer*4 NCOL, the column dimension of the matrix.
!    Input, integer*4 NZMAX, the maximum number of nonzero elements allowed.  
!                            This should be set to be the lengths of the arrays A and JA.
!    Input, real*8    DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!    Input, integer*4 NDNS, the first dimension of DNS, which must be at least NROW.
!    Output,real*8    A(*), the matrix in CSR Compressed Sparse Row format.
!           integer*4 JA(*), IA(NROW+1)
!    Output,integer*4 IERR, error indicator. 
!                           0 means normal return;
!                           i means that the the code stopped while processing row i, because
!                              there was no space left in A and JA, as defined by NZMAX.
  implicit none
  integer*4  ncol, ndns, nrow
  integer*4  i, ia(nrow+1), ierr
  integer*4  j, next, nzmax
  integer*4  ja(*)
  real*8     a(*)
  real*8     dns(ndns,ncol)
  ierr = 0
  next = 1
  ia(1) = 1

  do i = 1, nrow
    do j = 1, ncol
      if ( dns(i,j) /= 0.0D+00 ) then
        if ( nzmax < next ) then
          ierr = i
          return
        end if
        ja(next) = j
        a(next) = dns(i,j)
        next = next + 1
      end if
    end do
    ia(i+1) = next
  end do
  return
endsubroutine
!complex version of dnscsr_r routine
subroutine dnscsr ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )
  implicit none
  integer*4  ncol, ndns, nrow
  integer*4  i, ia(nrow+1), ierr
  integer*4  j, next, nzmax
  integer*4  ja(*)
  complex*16 a(*)
  complex*16 dns(ndns,ncol)
  ierr = 0
  next = 1
  ia(1) = 1

  do i = 1, nrow
    do j = 1, ncol
      if ( dns(i,j) /= 0.0D+00 ) then
        if ( nzmax < next ) then
          ierr = i
          return
        end if
        ja(next) = j
        a(next) = dns(i,j)
        next = next + 1
      end if
    end do
    ia(i+1) = next
  end do
  return
endsubroutine
!complex version of csrdns_r routine
subroutine csrdns ( nrow, ncol, a, ja, ia, dns, ndns, ierr )
!*****************************************************************************80
!! CSRDNS converts Compressed Sparse Row to Dense format.
!  Discussion:
!    This routine converts a row-stored sparse matrix into a densely stored one.
!  Modified: 07 January 2004
!  Author: Youcef Saad
!  Parameters:
!    Input, integer*4 NROW, the row dimension of the matrix.
!    Input, integer*4 NCOL, the column dimension of the matrix.
!    Input, complex*16 A(*), the matrix in CSR Compressed Sparse Row format.
!           integer*4 JA(*), IA(NROW+1)
!    Output,complex*16 DNS(NDNS,NDNS), the dense array containing a copy of the matrix.
!    Input, integer*4 NDNS, the dimension of the DNS array.
!    Output,integer*4 IERR, error indicator.
!                           0, means normal return
!                           i, means that the code has stopped when processing
!                              row number i, because it found a column number > ncol.
  implicit none
  integer*4  ncol, ndns
  integer*4  i, j, k, nrow, ierr
  complex*16 dns(ndns,ncol)
  complex*16 a(*)
  integer*4  ia(*)
  integer*4  ja(*)

  ierr = 0
  dns(1:nrow,1:ncol) = 0.0D+00
  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( ncol < j ) then
        ierr = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do
  return
endsubroutine

subroutine csrdns_r ( nrow, ncol, a, ja, ia, dns, ndns, ierr )
!*****************************************************************************80
!! CSRDNS converts Compressed Sparse Row to Dense format.
!  Discussion:
!    This routine converts a row-stored sparse matrix into a densely stored one.
!  Modified: 07 January 2004
!  Author: Youcef Saad
!  Parameters:
!    Input, integer*4 NROW, the row dimension of the matrix.
!    Input, integer*4 NCOL, the column dimension of the matrix.
!    Input, real*8    A(*), the matrix in CSR Compressed Sparse Row format.
!           integer*4 JA(*), IA(NROW+1)
!    Output,real*8    DNS(NDNS,NDNS), the dense array containing a copy of the matrix.
!    Input, integer*4 NDNS, the dimension of the DNS array.
!    Output,integer*4 IERR, error indicator.
!                           0, means normal return
!                           i, means that the code has stopped when processing
!                              row number i, because it found a column number > ncol.
  implicit none
  integer*4  ncol, ndns
  integer*4  i, j, k, nrow, ierr
  real*8     dns(ndns,ncol)
  real*8     a(*)
  integer*4  ia(*)
  integer*4  ja(*)

  ierr = 0
  dns(1:nrow,1:ncol) = 0.0D+00
  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( ncol < j ) then
        ierr = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do
  return
endsubroutine
subroutine aplb ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, iw, ierr )
!*****************************************************************************80
!! APLB performs the CSR matrix sum C = A + B.
!  Modified: 07 January 2004
!  Author: Youcef Saad
!  Parameters:
!    Input, integer*4  NROW, the row dimension of A and B.
!    Input, integer*4  NCOL, the column dimension of A and B.
!    Input, integer*4  JOB.  When JOB = 0, only the structure (i.e. the arrays jc, ic) is computed and the
!                            real values are ignored.
!    Input, real*8     A(*), B(*), the matrix in CSR Compressed Sparse Row format.
!           integer*4  JA(*), IA(NROW+1)
!           integer*4  JB(*), IB(NROW+1)
!           integer*4  nzmax, The  length of the arrays c and jc. 
!                             amub will stop if the result matrix C  has a number of elements that 
!                             exceeds exceeds nzmax. See ierr.
! on return:
!           real*8     C(*),  resulting matrix C in compressed sparse row sparse format.
!           integer*4  JC(*), IC(NROW+1)
!           integer*4  ierr,  serving as error message.
!                             ierr = 0 means normal return,
!                             ierr > 0 means that amub stopped while computing the
!                             i-th row  of C with i = ierr, because the number
!                             of elements in C exceeds nzmax.
! work arrays:
!           integer*4  iw,    work array of length equal to the number of columns in A.
!
  implicit none
  integer*4  ncol, nrow
  real*8     a(*)
  real*8     b(*)
  real*8     c(*)
  integer*4   ia(nrow+1), ib(nrow+1), ic(nrow+1)
  integer*4   ierr, ii, iw(ncol)
  integer*4   jcol, job, jpos, k, ka, kb, nzmax
  integer*4   len
  integer*4   ja(*)
  integer*4   jb(*)
  integer*4   jc(*)
  logical     values

  values = ( job /= 0 )
  ierr = 0
  len = 0
  ic(1) = 1
  iw(1:ncol) = 0

  do ii = 1, nrow
!  Row I.
     do ka = ia(ii), ia(ii+1)-1
        len = len + 1
        jcol = ja(ka)
        if ( nzmax < len ) then
          ierr = ii
          return
        end if
        jc(len) = jcol
        if ( values ) then
          c(len) = a(ka)
        end if
        iw(jcol) = len
     end do
     do kb = ib(ii), ib(ii+1)-1
        jcol = jb(kb)
        jpos = iw(jcol)
        if ( jpos == 0 ) then
           len = len + 1
           if ( nzmax < len ) then
             ierr = ii
             return
           end if
           jc(len) = jcol
           if ( values ) then
             c(len) = b(kb)
           end if
           iw(jcol)= len
        else
           if ( values ) then
             c(jpos) = c(jpos) + b(kb)
           end if
        end if
     end do
     do k = ic(ii), len
       iw(jc(k)) = 0
     end do
     ic(ii+1) = len+1
  end do

  return
endsubroutine

endmodule
