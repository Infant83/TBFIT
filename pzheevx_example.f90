subroutine pzheevx_example()
   implicit none
   integer*4 
   character*1 ORDER
   

   ORDER = 'R'
   NPROW = 2
   NPCOL = 2
   CALL BLACS_GET(0, 0, ICONTXT)
   CALL BLACS_GRIDINIT(ICONTXT, ORDER, NPROW, NPCOL)
   CALL BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYROW, MYCOL)
 
!              JOBZ RANGE UPLO  N  A IA JA  DESC_A   VL     VU   IL IU  ABSTOL  M  NZ  W   ORFAC
!               |     |     |   |  |  |  |    |       |      |    |  |     |    |   |  |     |
 CALL PZHEEVX( 'V',  'A',  'L', 4, A, 1, 1, DESC_A, 0.0D0, 0.0D0, 0, 0, -1.0D0, M, NZ, W, -1.0D0, &
                Z, 1, 1, DESC_Z, WORK , 0,  RWORK,  0,   IWORK ,  0,  IFAIL, ICLUSTR, GAP, INFO)
!               |  |  |    |      |     |     |     |      |      |     |       |      |     |
!               Z IZ JZ  DESC_Z  WORK LWORK RWORK LRWORK IWORK LIWORK IFAIL  ICLUSTR  GAP  INFO
return
endsubroutine
