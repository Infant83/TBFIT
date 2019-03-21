#include "alias.inc"
subroutine scalapack_initialize( SUMMRY, NOUT, N, NRHS, NB, NPROW, NPCOL, WORK, myid, nprocs_)
   use mpi_setup
      CHARACTER*( * )    SUMMRY
      INTEGER            IAM, N, NRHS, NB, NOUT, NPCOL, NPROW, nprocs_
      INTEGER            WORK( * )
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
      CHARACTER*79       USRINFO
      INTEGER            ICTXT

      EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT, BLACS_GRIDINIT, BLACS_SETUP, IGEBR2D, IGEBS2D

      IF( IAM.EQ.0 ) THEN
         OPEN( NIN, FILE='SCAEX.dat', STATUS='OLD' )
         READ( NIN, FMT = * ) SUMMRY
         SUMMRY = ' '
         READ( NIN, FMT = 9999 ) USRINFO
         READ( NIN, FMT = * ) SUMMRY
         READ( NIN, FMT = * ) NOUT
         IF( NOUT.NE.0 .AND. NOUT.NE.6 ) OPEN( NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )

         READ( NIN, FMT = * ) N
         READ( NIN, FMT = * ) NRHS
         READ( NIN, FMT = * ) NB
         READ( NIN, FMT = * ) NPROW
         READ( NIN, FMT = * ) NPCOL
         CLOSE( NIN )


         IF( NPROCS_.LT.1 ) THEN
            NPROCS_ = NPROW * NPCOL
            CALL BLACS_SETUP( IAM, NPROCS_ )
         END IF

         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS_ )

         WORK( 1 ) = N
         WORK( 2 ) = NRHS
         WORK( 3 ) = NB
         WORK( 4 ) = NPROW
         WORK( 5 ) = NPCOL
         CALL IGEBS2D( ICTXT, 'All', ' ', 5, 1, WORK, 5 )

         WRITE( NOUT, FMT = 9999 ) 'SCALAPACK example driver.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )'The matrices A and B are read from a file.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )'An explanation of the input/output parameters follows:'
         WRITE( NOUT, FMT = 9999 ) 'N       : The order of the matrix A.'
         WRITE( NOUT, FMT = 9999 ) 'NRHS    : The number of right and sides.'
         WRITE( NOUT, FMT = 9999 ) 'NB      : The size of the square blocks the matrices A and B are split into.'
         WRITE( NOUT, FMT = 9999 ) 'P       : The number of process rows.'
         WRITE( NOUT, FMT = 9999 ) 'Q       : The number of process columns.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 ) 'The following parameter values will be used:'
         WRITE( NOUT, FMT = 9998 ) 'N    ', N
         WRITE( NOUT, FMT = 9998 ) 'NRHS ', NRHS
         WRITE( NOUT, FMT = 9998 ) 'NB   ', NB
         WRITE( NOUT, FMT = 9998 ) 'P    ', NPROW
         WRITE( NOUT, FMT = 9998 ) 'Q    ', NPCOL
         WRITE( NOUT, FMT = * )

      ELSE

         IF( NPROCS_.LT.1 ) CALL BLACS_SETUP( IAM, NPROCS_ )

         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS_ )
         CALL IGEBR2D( ICTXT, 'All', ' ', 5, 1, WORK, 5, 0, 0 )
         N     = WORK( 1 )
         NRHS  = WORK( 2 )
         NB    = WORK( 3 )
         NPROW = WORK( 4 )
         NPCOL = WORK( 5 )

      END IF

      CALL BLACS_GRIDEXIT( ICTXT )
      RETURN

   20 WRITE( NOUT, FMT = 9997 )
      CLOSE( NIN )
      IF( NOUT.NE.6 .AND. NOUT.NE.0 ) CLOSE( NOUT )
      CALL BLACS_ABORT( ICTXT, 1 )
      STOP
 9999 FORMAT( A )
 9998 FORMAT( 2X, A5, '   :        ', I6 )
 9997 FORMAT( ' Illegal input in file ',40A,'.  Aborting run.' )
endsubroutine
