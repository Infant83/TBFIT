module do_math

   interface matprod
     module procedure :: matproduct_complex
     module procedure :: matproduct_real
     module procedure :: matproduct_complex_normal
     module procedure :: matproduct_real_normal
   end interface

contains

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
       write(6, '(A,A)')' got an error from function: ',func
       stop
    endif

    CALL ZGESVD( 'A', 'A', msize, msize, M, msize, S, U, msize, WT, msize, work, lwork, rwork, iflag)
    if (iflag .ne. 0) then
       write(6, '(A,A)')' got an error from function: ',func
       stop
    endif

!   call print_matrix_r(S,msize,msize,'Ssss',0)

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
    integer*4  msize, iflag
    complex*16 H(msize,msize)
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
      write(6, '(A,A)')' got an error from function: ',func
      stop
   endif

    if(flag_get_orbital) then
      !JOB = 'V'
      CALL ZHEEV( 'V', 'U', msize, H, msize, E, work, 16*msize, rwork, iflag )
    elseif(.not. flag_get_orbital) then
      !JOB = 'N'
      CALL ZHEEV( 'N', 'U', msize, H, msize, E, work, 16*msize, rwork, iflag )
    endif

return
end subroutine

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
      write(6, '(A,A)')' got an error from function: ',func
      stop
   endif

   call  ZGEEV( 'N', 'N', msize, H, msize, E, VL, msize, VR, msize, work, lwork, rwork, iflag )
   if (iflag .ne. 0) then
      write(6, '(A,A)')' got an error from function: ',func
      stop
   endif

   return
endsubroutine

subroutine set_identity_mat_c(msize, A)
   implicit none
   integer*4    msize
   integer*4    i
   complex*16   A(msize,msize)

   A = (0d0,0d0)

   do i = 1, msize

     A(i,i) = (1d0,0d0)

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
   integer*8, dimension(size(A,1)) :: ipiv ! pipov indices
   integer*8 :: n, info

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

endmodule
