#include "alias.inc"
subroutine test()
#ifdef SCALAPACK
   use mpi_setup
   implicit none
      integer*4          mpierr
      integer*4          IA, JA, I, J, K
      complex*16         alpha
      integer*4          msize, iflag
      integer*4          lwork, lrwork, liwork
      integer*4          iband, fband, nband, nband_, nband_found, nv
      real*8             vl, vu, abstol, orfac
      logical            flag_get_orbital
      character*1        JOBZ, TOP
      real*8, external:: PDLAMCH
      complex*16,allocatable :: myH(:, :)
      complex*16,allocatable :: myV(:, :)
      complex*16,allocatable :: work(:)
      real*8,    allocatable :: E(:), rwork(:), gap(:)
      integer*4, allocatable :: iwork(:), ifail(:), iclustr(:)
      integer*4, allocatable :: imap(:, :)
      character(*), parameter :: func = 'pcal_eig_hermitianx'
      integer*4          DESC_myH(DLEN_), DESC_myV(DLEN_)
      integer*4          CONTEXT, CONTEXT_, CONTEXT1 
      integer*4          MY_COL, MY_ROW, NB, NP_COL, NP_ROW
      integer*4          NPROW_, NPCOL_, MYROW_, MYCOL_
      integer*4          BLACS_PNUM
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, &
                         BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO, &
                         BLACS_SETUP, DESCINIT, PZHEEVX,&
                         PZLAPRNT, PZLAMODHILB, PZELGET, PZELSET, BLACS_PNUM

      msize = 5
      NB = 1
      call MPI_BARRIER(mpi_comm_earth, mpierr)
      do I = 0, nprow-1
        if(COMM_ASIA%mycoord(1) .eq. I) then
          call BLACS_GRIDINIT(COMM_ASIA%mpi_comm, 'R', 1, npcol)
          call BLACS_GRIDINFO(COMM_ASIA%mpi_comm, nprow, npcol, MY_ROW, MY_COL)
          write(6,'(4(A,I0))')"XXX MYID:",myid, ' COMM:', COMM_ASIA%mpi_comm, ' MY_ROW:', MY_ROW, ' MY_COL:', MY_COL
        endif
      enddo
    
      abstol = PDLAMCH(COMM_ASIA%mpi_comm, 'U')
      orfac = 0.001d0 ! default
      TOP = ' '
      kill_job

      if(CONTEXT_ .eq. 2) then
        call DESCINIT( DESC_myH, msize, msize, NB, NB, 0, 0, CONTEXT_, msize, iflag)
        CALL DESCINIT( DESC_myV, msize, msize, NB, NB, 0, 0, CONTEXT_, msize, iflag )

        CALL PZLAMODHILB( msize, myH, 1, 1, DESC_myH, iflag )
     
      endif
      kill_job

!     CALL DESCINIT( DESC_myH, msize, msize, NB, NB, 0, 0, CONTEXT, msize, iflag )
!     CALL DESCINIT( DESC_myV, msize, msize, NB, NB, 0, 0, CONTEXT, msize, iflag )

!     allocate(myH(ld_myH,ld_myH))
!     allocate(myV(ld_myV,ld_myV))

      CALL PZLAMODHILB( msize, myH, 1, 1, DESC_myH, iflag )
      CALL PZLAPRNT( msize, msize, myH, 1, 1, DESC_myH, 0, 0, 'myH', 6, work )

!     CALL PZELGET ('A', top, alpha, A, 2, 1, DESCA)

      allocate(iclustr(2*NPROW*NPCOL), gap(NPROW*NPCOL))
    ! first carry out a workspace query
    lwork = -1
    lrwork = -1
    liwork = -1
    ALLOCATE(work(1))
    ALLOCATE(rwork(1))
    ALLOCATE(iwork(1))
    CALL PZHEEVX(JOBZ,'I','U', msize, myH, 1, 1, DESC_myH, vl, vu, &
                 iband, fband, abstol, nband_found, nv, &
                 E, orfac, myV, 1, 1, DESC_myV, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 ifail, iclustr, gap, iflag)
    lwork = INT(ABS(work(1))) ; lrwork = INT(ABS(rwork(1))) ; liwork =INT (ABS(iwork(1)))
    DEALLOCATE(work)          ; DEALLOCATE(rwork)           ; DEALLOCATE(iwork)
    ALLOCATE(work(lwork))     ; ALLOCATE(rwork(lrwork))     ; ALLOCATE(iwork(liwork))

    ! actual calculation
    if(flag_get_orbital) JOBZ='V'
    call PZHEEVX(JOBZ,'I','U', msize, myH, 1, 1, DESC_myH, vl, vu, &
                 iband, fband, abstol, nband_found, nv, &
                 E, orfac, myV, 1, 1, DESC_myV, &
                 work, lwork, rwork, lrwork, iwork, liwork, &
                 ifail, iclustr, gap, iflag)
    if(nband_found .ne. nband) then
      write(6,'(A)') ' got an error from function: (nband .ne. nband_found)', func
      stop
    endif

!              JOBZ RANGE UPLO  N  A IA JA  DESC_A   VL     VU   IL IU  ABSTOL  M  NZ  W   ORFAC
!               |     |     |   |  |  |  |    |       |      |    |  |     |    |   |  |     |
!CALL PZHEEVX( 'V',  'A',  'L', 4, A, 1, 1, DESC_A, 0.0D0, 0.0D0, 0, 0, -1.0D0, M, NZ, W, -1.0D0, &
!               Z, 1, 1, DESC_Z, WORK , 0,  RWORK,  0,   IWORK ,  0,  IFAIL, ICLUSTR, GAP, INFO)
!               |  |  |    |      |     |     |     |      |      |     |       |      |     |
!               Z IZ JZ  DESC_Z  WORK LWORK RWORK LRWORK IWORK LIWORK IFAIL  ICLUSTR  GAP  INFO


20 CONTINUE

 kill_job
!call blacs_exit(1)
return
#endif
endsubroutine

#ifdef SCALAPACK
subroutine make_split()
use mpi_setup
implicit none
integer*4 CONTEXT_, mpierr, CONTEXT

!     call MPI_COMM_CREATE(COMM_EARTH%mpi_comm, COMM_ASIA%mpi_comm, CONTEXT_, mpierr)
!       write(6,*)"KKKK ", myid, CONTEXT_, COMM_EARTH%dims(2)
      if(COMM_EARTH%mycoord(2) .eq. 1) then
        call BLACS_GRIDINIT( COMM_ASIA%mpi_comm, 'Row', 1, COMM_EARTH%dims(2))
        write(6,*)"GGGG ", myid, CONTEXT
      endif

kill_job

return
endsubroutine
      SUBROUTINE PZLAMODHILB( N, A, IA, JA, DESCA, INFO )
      use mpi_setup
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
      INTEGER            I, J, MY_COL, MY_ROW, NP_COL, NP_ROW
      EXTERNAL           BLACS_GRIDINFO, PZELSET
      INTRINSIC          DBLE, DCMPLX
!     The matlab code for a real matrix is:
!         hilb(n) + diag( [ 1:-1/n:1/n ] )
!     The matlab code for a complex matrix is:
!         hilb(N) + toeplitz( [ 1 (1:(N-1))*i ] )
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DT_*LLD_*MB_*M_*NB_*N_* RSRC_.LT.0 )RETURN
      INFO = 0
      NP_ROW = nprow
      NP_COL = npcol
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NP_ROW, NP_COL, MY_ROW, MY_COL )
      IF( IA.NE.1 ) THEN
         INFO = -3
      ELSE IF( JA.NE.1 ) THEN
         INFO = -4
      END IF
      DO 20 J = 1, N
         DO 10 I = 1, N
            IF( I.EQ.J ) THEN
               CALL PZELSET( A, I, J, DESCA, DCMPLX( 1d0 / ( DBLE( I+J )-1d0 ) )+ DCMPLX( 1d0 ) )
            ELSE
               CALL PZELSET( A, I, J, DESCA, DCMPLX( 1d0 / ( DBLE( I+J )-1d0 ), DBLE( J-I ) ) )
            END IF
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
#endif
