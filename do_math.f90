#include "alias.inc"
module do_math
   use mpi_setup
   use print_io

#ifdef F08
   interface matprod
     module procedure :: matproduct_complex
     module procedure :: matproduct_real
     module procedure :: matproduct_complex_normal
     module procedure :: matproduct_real_normal
     module procedure :: matproduct_real_general
     module procedure :: matproduct_real_general_x1
     module procedure :: matproduct_complex_general
   end interface
#endif

contains
function matproduct_real_general_x1(M,K, JOBA, A, K_, JOBB, B) result(C)
   implicit none
   integer*4   M, K
   integer*4   K_,N
   real*8      alpha, beta
   character*1 JOBA, JOBB
   real*8,intent(in) :: A( M,  K), B( K_)
   real*8      C( M)
   character(*), parameter :: func = 'matproduct_real_general'

   ! computes alpha * A * B + beta * C = C
   alpha = 1d0
   beta  = 0d0
!  JOBA(or B) = 'N' or 'n',  op( A ) = A.
!  JOBA(or B) = 'T' or 't',  op( A ) = A**T.
!  JOBA(or B) = 'C' or 'c',  op( A ) = A**H.

   if( K .ne. K_) then
     write(message, '(A,A)')' got an error from function: ',func  ; write_msg
     write(message, '(A,I3,A,I3)')'  K .ne. K_:  K=', K, ' K_=',K_  ; write_msg
     stop
   endif

   call DGEMM(JOBA, JOBB, M, 1, K, alpha, A, M, B, K, beta, C, M)

   return
endfunction

function matproduct_real_general(M,K, JOBA, A, K_, N, JOBB, B) result(C)
   implicit none
   integer*4   M, K
   integer*4   K_,N 
   real*8      alpha, beta
   character*1 JOBA, JOBB
   real*8      A( M,  K), B( K_,  N)
   real*8      C( M,  N)
   character(*), parameter :: func = 'matproduct_real_general'
  
   ! computes alpha * A * B + beta * C = C
   alpha = 1d0
   beta  = 0d0
!  JOBA(or B) = 'N' or 'n',  op( A ) = A.
!  JOBA(or B) = 'T' or 't',  op( A ) = A**T.
!  JOBA(or B) = 'C' or 'c',  op( A ) = A**H.

   if( K .ne. K_) then
     write(message, '(A,A)')' got an error from function: ',func  ; write_msg
     write(message, '(A,I3,A,I3)')'  K .ne. K_:  K=', K, ' K_=',K_  ; write_msg
     stop
   endif

   call DGEMM(JOBA, JOBB, M, N, K, alpha, A, M, B, K, beta, C, M)

   return
endfunction

function matproduct_complex_general(M,K, JOBA, A, K_, N, JOBB, B) result(C)
   implicit none
   integer*4,intent(in) ::  M, K
   integer*4,intent(in) ::  K_,N
   complex*16  alpha, beta
   complex*16  A( M,  K), B( K_,  N)
   complex*16  C( M,  N)
   character*1,intent(in) :: JOBA, JOBB
   character(*), parameter :: func = 'matproduct_real_general'

   ! computes alpha * A * B + beta * C = C
   alpha = 1d0
   beta  = 0d0

   if( K .ne. K_) then
     write(message, '(A,A)')' got an error from function: ',func  ; write_msg
     write(message, '(A,I3,A,I3)')'  K .ne. K_:  K=', K, ' K_=',K_  ; write_msg
     stop
   endif

   call ZGEMM(JOBA, JOBB, M, N, K, alpha, A, M, B, K, beta, C, M)

   return
endfunction


function matproduct_real(msize, JOBA, A, JOBB, B) result(C)
   implicit none
   integer*4   msize
   real*8      A(msize, msize), B(msize, msize)
   real*8      C(msize, msize)
   real*8      alpha, beta
   character*1 JOBA, JOBB

   alpha = 1d0
   beta  = 0d0
   call DGEMM(JOBA, JOBB, msize, msize, msize, alpha , &
               A, msize, B, msize, beta , C, msize)

   return
endfunction

function matproduct_complex(msize,JOBA, A, JOBB, B) result(C)
   implicit none
   integer*4   msize
   complex*16  A(msize, msize), B(msize, msize)
   complex*16  C(msize, msize)
   complex*16  alpha, beta
   character*1 JOBA, JOBB
!ZGEMM  performs one of the matrix-matrix operations
!   C := alpha*op( A )*op( B ) + beta*C,
!where  op( X ) is one of
!   op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
!alpha and beta are scalars, and A, B and C are matrices, with op( A )
!an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!         JOBA or JOBB is CHARACTER*1
!          On entry, TRANSA specifies the form of op( A ) to be used in
!          the matrix multiplication as follows:
!             JOBA = 'N' or 'n',  op( A ) = A.
!             JOBA = 'T' or 't',  op( A ) = A**T.
!             JOBA = 'C' or 'c',  op( A ) = A**H.

   alpha = 1d0
   beta  = 0d0
   call ZGEMM(JOBA, JOBB, msize, msize, msize, alpha , &
               A, msize, B, msize, beta , C, msize)

   return
endfunction

function matproduct_real_normal(msize, A, B) result(C)
   implicit none
   integer*4   msize
   real*8      A(msize, msize), B(msize, msize)
   real*8      C(msize, msize)
   real*8      alpha, beta

   alpha = 1d0
   beta  = 0d0
   call DGEMM('N', 'N', msize, msize, msize, alpha , &
               A, msize, B, msize, beta , C, msize)

   return
endfunction

function matproduct_complex_normal(msize,A,B) result(C)
   implicit none
   integer*4   msize
   complex*16  A(msize, msize), B(msize, msize)
   complex*16  C(msize, msize)
   complex*16  alpha, beta

   alpha = 1d0
   beta  = 0d0
   call ZGEMM('N', 'N', msize, msize, msize, alpha , &
               A, msize, B, msize, beta , C, msize)

   return
endfunction

! generate diagonal msize x msize identity matrix
function diagonal(msize) result(A)
   implicit none
   integer*4,intent(in) :: msize
   integer*4               i
   real*8                  A(msize,msize)

   A = 0d0

   do i = 1, msize
    A(i,i) = 1d0
   enddo

   return
endfunction

!This function is to calculate determinant of the complex matrix
!The source is adoped from : https://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/
complex*16 function determinant(N, mat)
    implicit none
    integer*4, intent(in) :: N
    complex*16  mat(N,N)
    integer*4   i, info
    integer*4   ipiv(N)
    real*8      sgn

    ipiv = 0
    call ZGETRF(N, N, mat, N, ipiv, info)

    determinant = (1d0,0d0)
    do i = 1, N
        determinant = determinant*mat(i, i)
    end do

    sgn = 1d0
    do i = 1, N
        if(ipiv(i) /= i) then
            sgn = -sgn
        end if
    end do
    determinant = sgn*determinant

endfunction

integer*4 function determinant3i(mat)
    implicit none
    integer*4   mat(3,3)

    determinant3i =  mat(1,1) * ( mat(2,2) * mat(3,3) - mat(2,3) * mat(3,2) ) &
                    -mat(1,2) * ( mat(2,1) * mat(3,3) - mat(2,3) * mat(3,1) ) &
                    +mat(1,3) * ( mat(2,1) * mat(3,2) - mat(2,2) * mat(3,1) )

endfunction

integer*4 function tracei(mat,mat_size)
    implicit none
    integer*4   i, mat_size
    integer*4   mat(mat_size,mat_size)

    tracei = 0

    do i = 1, mat_size
      tracei = tracei + mat(i,i)
    enddo

endfunction

!The routine computes the singular value decomposition (SVD) 
!of a rectangular complex matrix M, optionally the left and/or 
!right singular vectors. The SVD is written as:
! M = U*S*WT
subroutine get_svd(M, U, WT, msize)
    implicit none
    integer*4  msize, iflag
    integer*4  lwork
    complex*16 M(msize,msize)
    real*8     S(msize,msize)
    complex*16 U(msize,msize)
    complex*16 WT(msize,msize)
    complex*16, allocatable :: work(:)
    real*8     rwork(12*msize)
    character(*), parameter :: func = 'get_svd'

    allocate(work(12*msize))

    lwork = -1
    CALL ZGESVD( 'A', 'A', msize, msize, M, msize, U, S, msize, WT, msize, work, lwork, rwork, iflag)
    if (iflag .eq. 0 .and. real(work(1)) .gt. 0) then
       lwork= work(1)
       deallocate(work)
       allocate(work(lwork))
    else
       write(message, '(A,A)')' got an error from function: ',func  ; write_msg
       stop
    endif

    CALL ZGESVD( 'A', 'A', msize, msize, M, msize, S, U, msize, WT, msize, work, lwork, rwork, iflag)
    if (iflag .ne. 0) then
       write(message, '(A,A)')' got an error from function: ',func  ; write_msg
       stop
    endif

!   call print_matrix_r(S,msize,msize,'Ssss',0)

   return
endsubroutine

! The routine computes all the eigenvalues, and optionally, the eigenvectors of a 
! complex generalized Hermitian positive-definite eigenproblem, of the form
! A*x = λ*B*x, 
! A*B*x = λ*x, or 
! B*A*x = λ*x.
! Here A and B are assumed to be Hermitian and B is also positive definite.
subroutine cal_gen_eig_hermitian(H, S, msize, E, flag_get_orbital)
    implicit none
    integer*4  itype
    integer*4, intent(in) ::  msize
    integer*4  iflag
    complex*16 H(msize,msize)
    complex*16 S(msize,msize)
    complex*16 H_(msize,msize)
    real*8     rwork(12*msize)
    real*8     E(msize)
    logical    flag_get_orbital
    character(*), parameter :: func = 'cal_eig_hermitian'
    integer*4  lwork
    complex*16, allocatable :: work(:)

   itype = 1

   if(msize .eq. 1) then
     E = H(1,1)
     return
   endif

   lwork = -1
   allocate(work(12*msize))
   call  ZHEGV( itype, 'N', 'U', msize, H, msize, S, msize, E, work, lwork, rwork, iflag )
   if (iflag .eq. 0 .and. real(work(1)) .gt. 0) then
      lwork= work(1)
      deallocate(work)
      allocate(work(lwork))
   else
      write(message, '(A,A)')' got an error from function: ',func  ; write_msg
      stop
   endif

    if(flag_get_orbital) then
      !JOB = 'V'
      CALL ZHEGV( itype, 'V', 'U', msize, H, msize, S, msize, E, work, lwork, rwork, iflag )
    elseif(.not. flag_get_orbital) then
      !JOB = 'N'
      H_ = H
      CALL ZHEGV( itype, 'N', 'U', msize, H_, msize, S, msize, E, work, lwork, rwork, iflag )
    endif

   return
endsubroutine

! The routine computes all the eigenvalues and, optionally, 
! the eigenvectors of a square complex Hermitian matrix H.
! The eigenvector v(j) of A satisfies the following formula:
! H*v(j) = E(j)*v(j) where E(j) is its eigenvalue.
! If "JOB" operator is 'V' corresponding eigenvector v is 
! returned in H matrix.
! The computed eigenvectors are orthonormal.
subroutine cal_eig_hermitian(H, msize, E, flag_get_orbital)
    implicit none
    integer*4, intent(in) ::  msize 
    integer*4  iflag
    complex*16 H(msize,msize)
    complex*16 H_(msize,msize)
    real*8     rwork(12*msize)
    real*8     E(msize)
    logical    flag_get_orbital
    character(*), parameter :: func = 'cal_eig_hermitian'
    integer*4  lwork
    complex*16, allocatable :: work(:)

   if(msize .eq. 1) then
     E = H(1,1)
     return
   endif

   lwork = -1
   allocate(work(12*msize))
   call  ZHEEV( 'N', 'U', msize, H, msize, E, work, lwork, rwork, iflag )
   if (iflag .eq. 0 .and. real(work(1)) .gt. 0) then
      lwork= work(1)
      deallocate(work)
      allocate(work(lwork))
   else
      write(message, '(A,A)')' got an error from function: ',func  ; write_msg
      stop
   endif

    if(flag_get_orbital) then
      !JOB = 'V'
      CALL ZHEEV( 'V', 'U', msize, H, msize, E, work, lwork, rwork, iflag )
    elseif(.not. flag_get_orbital) then
      !JOB = 'N'
      H_ = H
      CALL ZHEEV( 'N', 'U', msize, H_, msize, E, work, lwork, rwork, iflag )
    endif

return
end subroutine


! ZHEGVX Computes selected eigenvalues and, optionally, eigenvectors of a complex 
! generalized Hermitian positive-definite eigenproblem.
! The routine computes selected eigenvalues, and optionally, 
! the eigenvectors of a complex generalized Hermitian positive-definite eigenproblem, of the form
! A * X = lambda * B * X .
! Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range 
! of indices for the desired eigenvalues.
subroutine cal_gen_eig_hermitianx(H,S,msize,iband,nband,E,V,flag_get_orbital)
    implicit none
    integer*4  itype
    integer*4  msize, iflag
    integer*4  lwork
    integer*4  iband,fband,nband,nband_
    integer    ifail(msize), iwork(7*msize)
    real*8     rwork(12*msize)
    real*8     E(nband), E_(msize)
    real*8     vl,vu
    real*8     abstol
    logical    flag_get_orbital
    complex*16 H(msize,msize)
    complex*16 S(msize,msize)
    complex*16 V(msize,nband)
    complex*16, allocatable :: work(:)
    character(*), parameter :: func = 'cal_gen_eig_hermitianx'
    character*1 JOBZ
    real*8      DLAMCH
    external    DLAMCH
  
! itype = 1, the problem type is A*x = λ*B*x
! itype = 2, the problem type is A*B*x = λ*x
! itype = 3, the problem type is B*A*x = λ*x
   itype = 1

   abstol = DLAMCH('U')
   fband = iband + nband - 1
   nband_= nband
   JOBZ  = 'N'
   if(msize .eq. 1) then
     E = H(1,1)
     return
   endif

   lwork = -1
   allocate(work(12*msize))
   call ZHEGVX(itype,'N','I','U',msize,H,msize,S,msize,vl,vu,iband,fband,abstol,nband_,E_,V,msize,work,lwork,rwork,iwork,ifail,iflag)
   if (iflag .eq. 0 .and. real(work(1)) .gt. 0) then
      lwork= work(1)
      deallocate(work)
      allocate(work(lwork))
   else
      write(message, '(A,A)')' got an error from function: ',func  ; write_msg
      stop
   endif
   if(flag_get_orbital) JOBZ='V'
   call ZHEGVX(itype,JOBZ,'I','U',msize,H,msize,S,msize,vl,vu,iband,fband,abstol,nband,E_,V,msize,work,lwork,rwork,iwork,ifail,iflag)
   E(1:nband) = E_(1:nband)

return
end subroutine

! ZHEEVX computes selected eigenvalues and, optionally, eigenvectors
! of a complex Hermitian matrix H.  Eigenvalues and eigenvectors can
! be selected by specifying either a range of values or a range of
! indices for the desired eigenvalues.
subroutine cal_eig_hermitianx(H,msize,iband,nband,E,V,flag_get_orbital)
    implicit none
    integer*4  msize, iflag
    integer*4  lwork
    integer*4  iband,fband,nband,nband_
    integer    ifail(msize), iwork(7*msize)
    real*8     rwork(12*msize)
    real*8     E(nband), E_(msize)
    real*8     vl,vu
    real*8     abstol
    logical    flag_get_orbital
    complex*16 H(msize,msize)
    complex*16 V(msize,nband)
    complex*16, allocatable :: work(:)
    character(*), parameter :: func = 'cal_eig_hermitianx'
    character*1 JOBZ
    real*8      DLAMCH
    external    DLAMCH

!   ZOBZ   'N':  Compute eigenvalues only
!          'V':  Compute eigenvalues and eigenvectors
!   RANGE  'A': all eigenvalues will be found
!          'V': all eigenvalues in the half-open interval (VL,VU] will be found
!          'I': the IL-th through IU-th eigenvalues will be found
!   UPLO   'U':  Upper triangle of A is stored
!          'L':  Lower triangle of A is stored
!    N     The order of the matrix A.
!    A     On entry, the Hermitian matrix A
!          On exit, the lower triangle (if UPLO='L') or the upper triangle (if UPLO='U') of A, 
!          including the diagonal, is destroyed
!   LDA    The leading dimension of the array A.
!    VL    If RANGE='V', the lower bound of the interval to be searched for eigenvalues. VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!    VU    If RANGE='V', the upper bound of the interval to be searched for eigenvalues. VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!    IL    If RANGE='I', the index of the smallest eigenvalue to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!    IU    If RANGE='I', the index of the largest eigenvalue to be returned. 
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!  ABSTOL  The absolute error tolerance for the eigenvalues. 
!          An approximate eigenvalue is accepted as converged
!          when it is determined to lie in an interval [a,b]
!          of width less than or equal to ABSTOL + EPS *   max( |a|,|b| ),
!          where EPS is the machine precision.
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!          If this routine returns with INFO>0, indicating that some
!          eigenvectors did not converge, try setting ABSTOL to 2*DLAMCH('S').
!    M     The total number of eigenvalues found. 0 <= M <= N.
!          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!    W     On normal exit, the first M elements contain the selected eigenvalues in ascending order.
!    Z     If JOBZ = 'V', then if INFO = 0, the first M columns of Z contain the orthonormal eigenvectors of the matrix A
!          corresponding to the selected eigenvalues, with the i-th column of Z holding the eigenvector associated with W(i).
!          If an eigenvector fails to converge, then that column of Z contains the latest approximation to the eigenvector, 
!          and the index of the eigenvector is returned in IFAIL. If JOBZ = 'N', then Z is not referenced.
!          Note: the user must ensure that at least max(1,M) columns are supplied in the array Z; if RANGE = 'V', 
!          the exact value of M is not known in advance and an upper bound must be used.
!   LDZ    The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
!   WORK   On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!  LWORK   The length of the array WORK.  LWORK >= 1, when N <= 1; otherwise 2*N.
!          For optimal efficiency, LWORK >= (NB+1)*N, where NB is the max of the blocksize for ZHETRD and for
!          ZUNMTR as returned by ILAENV. 
!          If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, 
!          returns this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
!  RWORK   RWORK is DOUBLE PRECISION array, dimension (7*N)
!  IWORK   IWORK is INTEGER array, dimension (5*N)
!  IFAIL   If JOBZ = 'V', then if INFO = 0, the first M elements of IFAIL are zero.  If INFO > 0, then IFAIL contains the
!          indices of the eigenvectors that failed to converge. If JOBZ = 'N', then IFAIL is not referenced.
!  INFO    = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, then i eigenvectors failed to converge. Their indices are stored in array IFAIL.
!  abstol=1.d-10
!  abstol=2.0d0*tiny(abstol)
   abstol = DLAMCH('U')
   fband = iband + nband - 1
   nband_= nband
   JOBZ  = 'N'
   if(msize .eq. 1) then
     E = H(1,1)
     return
   endif
   
   lwork = -1
   allocate(work(12*msize))
   call ZHEEVX('N','I','U',msize,H,msize,vl,vu,iband,fband,abstol,nband_,E_,V,msize,work,lwork,rwork,iwork,ifail,iflag)
   if (iflag .eq. 0 .and. real(work(1)) .gt. 0) then
      lwork= work(1)
      deallocate(work)
      allocate(work(lwork))
   else
      write(message, '(A,A)')' got an error from function: ',func  ; write_msg
      stop
   endif
   if(flag_get_orbital) JOBZ='V'
   call ZHEEVX(JOBZ,'I','U',msize,H,msize,vl,vu,iband,fband,abstol,nband,E_,V,msize,work,lwork,rwork,iwork,ifail,iflag)
   E(1:nband) = E_(1:nband)
return
endsubroutine
#ifdef SCALAPACK
! subroutine pcal_eig_hermitianx(H,msize,iband,nband,E,V,flag_get_orbital)
!   use mpi_setup
!   implicit none
!   integer*4  msize, iflag
!   integer*4  lwork, lrwork, liwork
!   integer*4  iband,fband,nband,nband_, nband_found, nv
!   integer    ifail(msize)
!   real*8     E(nband)
!   real*8     vl,vu
!   real*8     abstol, orfac
!   logical    flag_get_orbital
!   complex*16 H(msize,msize)
!   complex*16 V(msize,nband)
!   complex*16 myV(msize,nband)
!   complex*16, allocatable :: work(:)
!   real*8,     allocatable :: rwork(:), gap(:)
!   integer*4,  allocatable :: iwork(:), iclustr(:)
!   character(*), parameter :: func = 'pcal_eig_hermitianx'
!   character*1 JOBZ
!   real*8      PDLAMCH
!   external    PDLAMCH

!   fband = iband + nband - 1
!   nband_= nband
!   JOBZ  = 'N'
!   if(msize .eq. 1) then
!     E = H(1,1)
!     return
!   endif

!   !this setting yields the most orthogonal eigenvectors
!   abstol = PDLAMCH(CONTEXT, 'U')
!   orfac = 0.001d0 ! default

!   !allocate arrays used in PZHEEVX
!   allocate(iclustr(2*NPROW*NPCOL), gap(NPROW*NPCOL))

!   ! first carry out a workspace query
!   lwork = -1
!   lrwork = -1
!   liwork = -1
!   ALLOCATE(work(1))
!   ALLOCATE(rwork(1))
!   ALLOCATE(iwork(1))
!   CALL PZHEEVX(JOBZ,'I','U', N, myH, 1, 1, DESC_myH, vl, vu, &
!                iband, fband, abstol, nband_found, nv, &
!                E, orfac, myV, 1, 1, DESC_myV, &
!                work, lwork, rwork, lrwork, iwork, liwork, &
!                ifail, iclustr, gap, iflag)
!   lwork = INT(ABS(work(1))) ; lrwork = INT(ABS(rwork(1))) ; liwork =INT (ABS(iwork(1)))
!   DEALLOCATE(work)          ; DEALLOCATE(rwork)           ; DEALLOCATE(iwork)
!   ALLOCATE(work(lwork))     ; ALLOCATE(rwork(lrwork))     ; ALLOCATE(iwork(liwork))

!   ! actual calculation
!   if(flag_get_orbital) JOBZ='V'
!   call PZHEEVX(JOBZ,'I','U', N, myH, 1, 1, DESC_myH, vl, vu, &
!                iband, fband, abstol, nband_found, nv, &
!                E, orfac, myV, 1, 1, DESC_myV, &
!                work, lwork, rwork, lrwork, iwork, liwork, &
!                ifail, iclustr, gap, iflag)
!   if(nband_found .ne. nband) then
!     write(message,'(A)') ' got an error from function: (nband .ne. nband_found)', func  ; write_msg
!     stop
!   endif

! endsubroutine
#endif
#ifdef MKL_SPARSE
subroutine cal_gen_eig_hermitianx_sparse(SHk,SSk,emin,emax,nemax,ne_found,ne_guess,E,V,flag_vector,fpm,iflag,ne_prev)
    use parameters, only: spmat
    use time
    implicit none
    type(spmat) :: SHk, SSk
    integer*4   msize, iflag
    integer*4   ne_found, nemax, ne_guess, ne_prev
    integer*4   loop, fpm(128)
    integer*4   iter, max_iter
    real*8      emin, emax
    real*8      epsout
    logical     flag_vector
    real*8      E_(nemax)
    real*8      E(nemax)
    real*8      res(nemax)
    complex*16  V_(SHk%msize,nemax)
    complex*16  V(SHk%msize,nemax)
    character(*), parameter :: func = 'cal_eig_hermitianx_sparse'
    character*1 UPLO
    logical     flag_success
    real*8      t1, t0

    flag_success = .false.
    max_iter = 5; iter = 1

    msize = SHk%msize

    if(fpm(5) .eq. 1) then
      V_ = 0d0
      if(ne_prev .ge. 1) V_(:,1:ne_prev) = V(:, 1:ne_prev)
    endif

    UPLO = 'F'
    ! FEAST eigensolver for Complex and Hermitian Sparse 3-array matrix
    ! Extended Eigensolver interface for generalized eigenvalue problem with sparse matrices
    do while (.not. flag_success .and. iter .le. 5)
      !call zfeast_hcsrgv(uplo, n, a, ia, ja, b, ib, jb, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
      !call zfeast_hcsrev(uplo, n, a, ia, ja,            fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)

      call zfeast_hcsrgv(UPLO, msize, SHk%H, SHk%I, SHk%J, SSk%H, SSk%I, SSk%J, fpm, epsout, loop, emin, emax, ne_guess, &
                         E_, V_, ne_found, res, iflag)
      call report_error_feast_scsrev(iflag, fpm, flag_success, iter, max_iter, emin, emax, ne_guess, ne_found, nemax)
    enddo

    if(ne_found .ge. 1) then
      E = 0d0
      E(1:ne_found) = E_(1:ne_found)

      if(fpm(5) .eq. 1) then
        V = 0d0
        V(:, 1:ne_found) = V_(:, 1:ne_found)
      else
        if(flag_vector) then
          V = 0d0
          V(:, 1:ne_found) = V_(:, 1:ne_found)
        endif
      endif

    elseif(ne_found .eq. 0) then
      E = 0d0

      if(fpm(5) .eq. 1) then
        V = V_ !return to V_ if no result obtained
      else
        if(flag_vector) then
          V = 0d0
        endif
      endif

    endif


    return
endsubroutine
subroutine cal_eig_hermitianx_sparse(SHk, emin,emax,nemax,ne_found,ne_guess,E,V,flag_vector, fpm, iflag, ne_prev)
    use parameters, only: spmat
    use time
    implicit none
    type(spmat) :: SHk
    integer*4   msize, iflag
    integer*4   ne_found, nemax, ne_guess, ne_prev
    integer*4   loop, fpm(128)
    integer*4   iter, max_iter
    real*8      emin, emax
    real*8      epsout
    logical     flag_vector
    real*8      E_(nemax)
    real*8      E(nemax)
    real*8      res(nemax)
    complex*16  V_(SHk%msize,nemax)
    complex*16  V(SHk%msize,nemax)
    character(*), parameter :: func = 'cal_eig_hermitianx_sparse'
    character*1 UPLO
    logical     flag_success
    real*8      t1, t0

    flag_success = .false.
    max_iter = 5; iter = 1

    msize = SHk%msize

    if(fpm(5) .eq. 1) then
      V_ = 0d0
      if(ne_prev .ge. 1) V_(:,1:ne_prev) = V(:, 1:ne_prev)
    endif

    UPLO = 'F'
    ! FEAST eigensolver for Complex and Hermitian Sparse 3-array matrix
    do while (.not. flag_success .and. iter .le. 5)
      call zfeast_hcsrev(UPLO, msize, SHk%H, SHk%I, SHk%J, fpm, epsout, loop, emin, emax, ne_guess, &
                         E_, V_, ne_found, res, iflag)
      call report_error_feast_scsrev(iflag, fpm, flag_success, iter, max_iter, emin, emax, ne_guess, ne_found, nemax)
    enddo

    if(ne_found .ge. 1) then
      E = 0d0
      E(1:ne_found) = E_(1:ne_found)

      if(fpm(5) .eq. 1) then
        V = 0d0
        V(:, 1:ne_found) = V_(:, 1:ne_found)
      else
        if(flag_vector) then
          V = 0d0
          V(:, 1:ne_found) = V_(:, 1:ne_found)
        endif
      endif

    elseif(ne_found .eq. 0) then
      E = 0d0

      if(fpm(5) .eq. 1) then
        V = V_ !return to V_ if no result obtained
      else
        if(flag_vector) then
          V = 0d0
        endif
      endif

    endif

return
endsubroutine
subroutine report_error_feast_scsrev(iflag, fpm, flag_success, iter, max_iter, emin, emax, ne_guess, ne_found, nemax)
    implicit none
    integer*4    iflag, idummy
    integer*4    iter, max_iter
    logical      flag_success
    integer*4    fpm(128)
    character*8  argument(16)
    real*8       emin, emax
    integer*4    ne_guess, ne_found   
    integer*4    nemax

    argument(1)  = 'UPLO    '
    argument(2)  = 'msize   '
    argument(3)  = 'SHk     '
    argument(4)  = 'I       '
    argument(5)  = 'J       '
    argument(6)  = 'fpm     '
    argument(7)  = 'epsout  '
    argument(8)  = 'loop    '
    argument(9)  = 'emin    '
    argument(10) = 'emax    '
    argument(11) = 'E_      '
    argument(12) = 'V_      '
    argument(13) = 'ne_guess'
    argument(14) = 'ne_found'
    argument(15) = 'res     '
    argument(16) = 'iflag   '

    select case(iflag)

      ! Error messages. Program will be stop.
      case(202)
        write(message,'(A)'               )'   !ERROR! feast_scsrev: IFLAG=202,   Problem with size of the system "msize"'  ; write_msg
        stop
      case(201)                    
        write(message,'(A)'               )'   !ERROR! feast_scsrev: IFLAG=201,   Problem with size of subspace   "ne_guess"'  ; write_msg
        stop
      case(200)                    
        write(message,'(A)'               )'   !ERROR! feast_scsrev: IFLAG=200,   Problem with "emin", "emax"'  ; write_msg
        stop
      case(100:199)
        idummy = iflag - 100
        write(message,'(A,I0,A,I0,A     )')'   !ERROR! feast_scsrev: IFLAG= 100+',idummy,', Problem with ',idummy,'-th value'  ; write_msg
        write(message,'(A,I0,A)'          )'           of the input FEAST parameter (i.e. fpm(',idummy,'))'  ; write_msg
        stop

      ! Warning messages. Program will be continue.
      case(6  )                    
        write(message,'(A)'               )'   !WARN!  feast_scsrev: IFLAG=6  ,   FEAST converges but subspace is not bi-orthogonal'  ; write_msg
        flag_success = .true.
      case(5  )                    
        write(message,'(A)'               )'   !WARN!  feast_scsrev: IFLAG=5  ,   Only stochastic estimation of #eigenvalues'  ; write_msg
        write(message,'(A)'               )'                                      returned fpm(14)=2'  ; write_msg
        flag_success = .true.
      case(4  )                    
        write(message,'(A)'               )'   !WARN!  feast_scsrev: IFLAG=4  ,   Only the subspace has been returned using fpm(14)=1'  ; write_msg
        flag_success = .true.
      case(3  )                    
        write(message,'(A)'               )'   !WARN!  feast_scsrev: IFLAG=3  ,   Size of the subspace "NE_GUESS" is too small'  ; write_msg
        write(message,'(A)'               )'                                      The proper condition is: 0 <= NE * 1.5 < NE_MAX <= NEIG'  ; write_msg
        write(message,'(A)'               )'                                      The eigenvalues less than NE_MAX will be stored.'  ; write_msg
        write(message,'(A)'               )'                                      Please increase NE_GUESS.'  ; write_msg
        write(message,'(A,I0)'            )'                                      NE_FOUND = ', ne_found  ; write_msg
        write(message,'(A,I0)'            )'                                      NE_GUESS = ', ne_guess  ; write_msg
        write(message,'(A   )'            )'                                       ==> NE_GUES_new = ne_guess + ceiling(0.2*NE_MAX)'  ; write_msg
        write(message,'(A,I0)'            )'                                                       = ',ne_guess + ceiling(0.2*nemax)  ; write_msg
        ne_guess = ne_guess + ceiling(0.1*nemax)
        flag_success = .false.
        iter = iter + 1
      case(2  )                    
        write(message,'(A)'               )'   !WARN!  feast_scsrev: IFLAG=2  ,   No Convergence (#iteration loops > fpm(4))'  ; write_msg
        write(message,'(A, I0,A, I0)'     )'           -> increase refinement loops from ', fpm(4),' to ', fpm(4) + 10  ; write_msg
        fpm(4) = fpm(4) + 10
        flag_success = .false.
        iter = iter + 1
      case(1  )                    
!       write(message,'(A)'               )'   !WARN!  feast_scsrev: IFLAG=1  ,   No Eigenvalues found in the search interval:'  ; write_msg
!       write(message,'(2(A,F10.5),A)'    )'                                      [EMIN:EMAX] = [',emin,':',emax,']'  ; write_msg
        flag_success = .true. 
      ! Sucessful exit
      case(0  )                    
        flag_success = .true.
!       write(message,'(A)'               )'   !WARN!  feast_scsrev: IFLAG=1  ,   No Eigenvalues found in the search interval'  ; write_msg
      
      ! Error messages. Program will be stop.
      case(-1 )                    
        write(message,'(A)'               )'   !ERROR! feast_scsrev: IFLAG=-1 ,   Internal error for allocation memory'  ; write_msg
        stop
      case(-2 )                    
        write(message,'(A)'               )'   !ERROR! feast_scsrev: IFLAG=-2 ,   Internal error of the inner system solver in'  ; write_msg
        write(message,'(A)'               )'                                      FEAST predefinded interfaces'  ; write_msg
        stop
      case(-3 )                    
        write(message,'(A)'               )'   !ERROR! feast_scsrev: IFLAG=-3 ,   Internal error of the reduced eigenvalue solver'  ; write_msg
        write(message,'(A)'               )'                                      Possible cause for Hermitian problem:'  ; write_msg
        write(message,'(A)'               )'                                          -->  matrix H may not be positive definite'  ; write_msg
        stop
      case(-120:-101 )                    
        idummy = (iflag + 100) * -1
        write(message,'(A,I2,A,I0,A)'     )'   !ERROR! feast_scsrev: IFLAG=-(100+',idummy,'),  Problem with the ',idummy,'-th argument'  ; write_msg
        write(message,'(A,A)'             )'                                           of the FEAST interface:',trim(argument(idummy))  ; write_msg
        stop
    end select

    if(.not. flag_success .and. iter .gt. max_iter) then
        write(message,'(A          )'     )'   !ERROR! feast_scsrev: Calculation has not been successfully finished until "max_iter".'  ; write_msg
        write(message,'(A          )'     )'           Please increase "max_iter" parameter of your source code or check your input again.'  ; write_msg
        write(message,'(A          )'     )'           Exit program...'  ; write_msg
        stop
    endif

return
endsubroutine
#endif
!  ZGEEV computes for an N-by-N complex nonsymmetric matrix H, the
!  eigenvalues and, optionally, the left and/or right eigenvectors.
!  The right eigenvector VR of H satisfies
!  H*VR = E*VR where E is its eigenvalue.
!  The left eigenvector VL of H satisfies
!  VLT*H = E*VLT where VLT denotes the conjugate transpose (T) of VL. 
!  The computed eigenvectors are normalized to have Euclidean norm 
!  equal to 1 and largest component real.
subroutine cal_eig_nonsymm(H, msize, E)
   implicit none
   integer*4    msize
   integer*4    lwork
   complex*16   H(msize,msize)
   complex*16   E(msize)
   complex*16   VL(msize,msize), VR(msize,msize)
   complex*16, allocatable :: work(:)
   real*8       rwork(12*msize)
   integer*4    iflag
   character(*), parameter :: func = 'cal_eig_nonsym'

   allocate(work(12*msize))

   if(msize .eq. 1) then
     E = H(1,1)
     return
   endif

   lwork = -1
   call  ZGEEV( 'N', 'N', msize, H, msize, E, VL, msize, VR, msize, work, lwork, rwork, iflag )
   if (iflag .eq. 0 .and. real(work(1)) .gt. 0) then
      lwork= work(1)
      deallocate(work)
      allocate(work(lwork))
   else
      write(message, '(A,A)')' got an error from function: ',func  ; write_msg
      stop
   endif

   call  ZGEEV( 'N', 'N', msize, H, msize, E, VL, msize, VR, msize, work, lwork, rwork, iflag )
   if (iflag .ne. 0) then
      write(message, '(A,A)')' got an error from function: ',func  ; write_msg
      stop
   endif

   return
endsubroutine

subroutine set_identity_mat_c(msize, A)
   implicit none
   integer*4    msize
   integer*4    i
   complex*16   A(msize,msize)

   A = (0.d0,0.d0)

   do i = 1, msize

     A(i,i) = (1.d0,0.d0)

   enddo

   return
endsubroutine

subroutine set_identity_mat_r(msize, A)
   implicit none
   integer*4    msize
   integer*4    i
   real*8       A(msize,msize)

   A = 0d0

   do i = 1, msize

     A(i,i) = 1d0

   enddo
   
   return
endsubroutine

elemental real*8 function arg(z)

   use parameters, only : pi
   complex*16,intent(in) :: z

   arg = atan2(aimag(z),real(z))
!  arg = aimag(log(z))

   return
endfunction

subroutine rotate_vector(rot_m, local_vector_radial_coord)
   use parameters
   implicit none
   real*8   local_vector_radial_coord(3) ! local_vector(1) magnitude of the local_vector (negative sign indicates spin minority)
                                         ! local_vector(2) theta for in (deg) unit,
                                         ! local_vector(3) phi   for in (deg) unit. (2) and (3) is used for the rotation
   real*8   theta_rad, phi_rad
   real*8   rot_m(3)        ! 1: x, 2:y, 3:z ! vector returned with rotated

   theta_rad = local_vector_radial_coord(2) * pi2 / 360d0
   phi_rad   = local_vector_radial_coord(3) * pi2 / 360d0

   rot_m(1)= local_vector_radial_coord(1) * cos(theta_rad)*sin(phi_rad)
   rot_m(2)= local_vector_radial_coord(1) * sin(theta_rad)*sin(phi_rad)
   rot_m(3)= local_vector_radial_coord(1) * cos(phi_rad)

return
endsubroutine

! Returns the inverse of a matrix calculated by finding the LU
! decomposition.
! this routine is copied from: http://fortranwiki.org/fortran/show/Matrix+inversion
function inv(A) result(Ainv)
   implicit none
   real*8, dimension(:,:), intent(in) :: A
   real*8, dimension(size(A,1),size(A,2)) :: Ainv
   real*8, dimension(size(A,1)) :: work ! work array for LAPACK
   integer*4, dimension(size(A,1)) :: ipiv ! pipov indices
   integer*4 :: n, info

   ! External procedures defined in LAPACK
   external DGETRF
   external DGETRI

   ! Store A in Ainv to prevent it from being overwritten by LAPACK
   Ainv = A
   n = size(A,1)

   ! DGETRF computes an LU factorization of a general M-by-N matrix A
   ! using partial pivoting with row interchanges.
   call DGETRF(n, n, Ainv, n, ipiv, info)

!  if (info /= 0) then
!     stop 'Matrix is numerically singular!'
!  end if

   ! DGETRI computes the inverse of a matrix using the LU factorization
   ! computed by DGETRF.
   call DGETRI(n, Ainv, n, ipiv, work, n, info)

!  if (info /= 0) then
!     stop 'Matrix inversion failed!'
!  end if

endfunction

function invc(A) result(Ainv)
   implicit none
   complex*16, dimension(:,:), intent(in) :: A
   complex*16, dimension(size(A,1), size(A,2)) :: Ainv
   complex*16, dimension(size(A,1)) :: work
   integer*4,  dimension(size(A,1)) :: ipiv
   integer*4                        :: n, info
   external ZGETRF
   external ZGETRI

   Ainv = A
   n = size(A,1)

   call ZGETRF(n, n, Ainv, n, ipiv, info)

   call ZGETRI(n, Ainv, n, ipiv, work, n, info)

return
endfunction

function area(a,b) result(S)
   implicit none
   real*8   S
   real*8   a(3)
   real*8   b(3)
   real*8   axb(3)

   call vcross(axb,a,b)

   S = dsqrt( dot_product(axb,axb) )

   return
endfunction

! caculate the degree <abc (theta) between two vector ca and cb
function angle(a,b,c) result (theta)
   use parameters,  only: pi
   implicit none
   integer*4     i
   real*8        a(3), b(3), c(3)
   real*8        ba(3), cb(3)
   real*8        theta 
   real*8        enorm
   external      enorm

   ba = a - b
   cb = c - b

   theta = acos ( dot_product(ba,cb)/ ( enorm(3,ba) * enorm(3,cb) ) ) * 180d0 / pi

   return 
endfunction

function degen(E) result(D)
   implicit none
   integer*4   n
   real*8, dimension(:), intent(in) :: E
   real*8, dimension(3,size(E))     :: D
!  real*8, dimension(  size(E))     :: delta_above
!  real*8, dimension(  size(E))     :: delta_below

   n                  = size(E)
   D                  = 0d0

   ! delta_above
   D(2,1:n-1)         = E(2:n  ) - E(1:n-1)
   D(2,  n  )         = 1d0

   ! delta_below
   D(3,2:n  )         = D(2,1:n-1)
   D(3,1    )         = 1d0
 
   ! delta_below * delta_above
   D(1, :   )         = D(2,:) * D(3,:)

   return
endfunction

! gaussian distribution function
elemental real*8 function fgauss(sigma, x)
   use parameters, only : pi, pi2
   implicit none
   real*8,intent(in) :: sigma
   real*8,intent(in) :: x
   real*8               sigma2, xx

   xx = x**2 ; sigma2 = sigma**2

  !fgauss= exp(-0.5d0*xx/sigma2)/(sigma*sqrt(pi2))
   fgauss= 1d0 / (sigma*sqrt(pi2)) * exp( -xx / (2d0*sigma2) )

return
end function

! fermi-dirac distribution function
! k_B=8.6173303d-5  !! Boltzmann constant = 2m/hbar**2 [1/eV Ang^2]   
elemental real*8 function fdirac(x, u, T)
   use parameters, only : k_B, eta
   real*8,intent(in) ::  x ! = e - ef
   real*8,intent(in) ::  u ! chemical potential 
   real*8,intent(in) ::  T
   real*8                kT
   
   kT = k_B * (T + eta) ! the small value eta is introduced to avoid the numerical error

   fdirac = 1d0 / ( 1d0 + exp( (x-u) / kT) )

return
endfunction

! gradient of fermi-dirac distribution function w.r.t. the chemical potential
elemental real*8 function grad_u_fdirac(x, u, T)
   use parameters, only : k_B, eta
   real*8,intent(in) ::  x ! = e - ef
   real*8,intent(in) ::  u ! chemical potential 
   real*8,intent(in) ::  T
   real*8                kT

   ! Note: using the relation cosh(2a) + 1 = 2*cosh(a))**2 and 
   !       derivative rule with f(x) = 1/(1+exp(x)) -> f' = exp(x)*x' / (1+exp(x))**2,
   !       grad_u fdirac = grad_u (1 + exp( (x - u) / kT))**-1
   !                     =        (2*cosh( (e - u)/2kT)**-2
   !

   kT = k_B * (T + eta) + eta ! the small value eta is introduced to avoid the numerical error
   grad_u_fdirac = 1d0/ kT / (2d0 * cosh((x - u)/(2d0*kT)))**2

return
endfunction


! return index of original 2D array from reshaped 1D array index from 2D array
elemental integer*4 function idx2Dj(idx1D, ni)
   implicit none
   integer*4,intent(in) :: ni
   integer*4,intent(in) :: idx1D

   idx2Dj = floor( real(idx1D - 1)/real(ni) ) + 1
return
endfunction
elemental integer*4 function idx2Di(idx1D, ni)
   implicit none
   integer*4,intent(in) :: ni
   integer*4,intent(in) :: idx1D

   idx2Di = mod( idx1D - 1, ni ) + 1
return
endfunction

endmodule
