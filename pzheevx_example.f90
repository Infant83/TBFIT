PROGRAM PZHEEVX

USE useful_defs
IMPLICIT NONE
INCLUDE 'mpif.h'
!___________________________________________________________________________________

  INTEGER :: LWORK, LRWORK,LIWORK 
  INTEGER :: MAXN  
  INTEGER :: LDA 
  INTEGER :: LDB
  REAL(float) :: diff
!_______________________________________________________________________________________
  INTEGER :: i1,i2,i_loc,j_loc
  INTEGER :: nx,ny,N
  REAL (float) :: t(3)
  COMPLEX (float), ALLOCATABLE :: H(:,:),U(:,:)
  REAL (float) :: E(:)

  ! for density of states
  TYPE (dos_defs) :: dos
 
  !local arrays
   COMPLEX (float),ALLOCATABLE:: LH(:,:), LZ(:,:)
   

   INTEGER  NUMROC,ICEIL
   REAL PDLAMCH
!_____________________________________________________________________________________
!     .. Local Scalars ..
      INTEGER            CONTEXT, IAM, INFO, MYCOL, MYROW, NBROW,NBCOL, NPCOL, NPROCS, NPROW 

!     .. Local Arrays ..
      INTEGER            DESCH( 50 ), DESCZ( 50 ), IWORK(:), IFAIL(:), ICLUSTR(:)
      COMPLEX (FLOAT), ALLOCATABLE ::   A(:,:), WORK(: ),Z( :,: )
      REAL (FLOAT), ALLOCATABLE :: RWORK(:), W(:), GAP(:)

      REAL (FLOAT) :: ABSTOL, ORFAC, VL, VU
      INTEGER :: IL, IU, M, NZ, NNP,NN,NPO,MQO

!     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,BLACS_SETUP, DESCINIT, PDLAMODHILB, PDLAPRNT,PZHEEVX, NUMROC, PDLAMCH,ICEIL


      ! network topology
      CHARACTER (1):: TOP =' '

      INTEGER :: NB,MB,NBROWS,NBCOLUMNS,IBROWS,IBCOLUMNS,COL_MAX,ROW_MAX,ROWS,COLUMNS,ICOL,IROW,LOCr,LOCc
      INTEGER:: ISTATUS(MPI_STATUS_SIZE),REQUEST, IERR,RANK,NODE

      INTRINSIC  CEILING
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     .. Executable Statements ..
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Set up the problem

      NBROW = 64         ! number of rows in a block of A(global array)
      NBCOL = 64         ! number of columns in a block of A(global array)
      NPROW = 1          ! number of rows in the processor grid
      NPCOL = 3          ! number of columns in the processor grid

!     Initialize the BLACS

      CALL SL_init(CONTEXT,NPROW,NPCOL)
      CALL mpi_comm_size(mpi_comm_world,nprocs,info)

!     check that the number of process passed as the argument matches the number of 
!     processes in the processor grid
      IF (nprocs /= nprow * npcol) THEN
         WRITE(*,*) 'Error! Number of processors passed does not match with processors in the grid'
         STOP
      END IF 

      !each processor gets its rank
      CALL  MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERR)

!     Initialize a single BLACS context
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )


      !   allocate memory for the global arrays
      CALL init_vals(nx,ny,N,t,dos)
      ALLOCATE(H(N,N),U(N,N),E(N),dos%total_dos(dos%nbin))
  

!    every process should initialize LDA and LDB
     LDA = NUMROC(N,NBROW,MYROW,0,NPROW)
     LDB = NUMROC(N,NBCOL,MYCOL,0,NPCOL)

     !    check that each processor has the right matrix

      ALLOCATE(A(N,N), W(N), Z(N,N))

      !   allocate memory for the local arrays
      ALLOCATE(LH(LDA,LDB), LZ(LDA,LDB))
  
!     These are basic array descriptors

      CALL DESCINIT( DESCH, N, N, NBROW, NBCOL, 0, 0, CONTEXT, LDA, INFO )
      CALL DESCINIT( DESCZ, N, N, NBROW, NBCOL, 0, 0, CONTEXT, LDA, INFO )
______________________________________________________________________________________________________

! only one processor with coordinates(0,0) executes this part of the code

      IF(myrow.EQ.0.AND.mycol.EQ.0) THEN

         ! Make the Hamiltonian matrix
         CALL make_H(nx,ny,N,t,H)

!        now distribute the hamiltonian matrix to each processor (send the matrix)
         CALL ZGEBS2D(CONTEXT,'A',TOP,N,N,H,N)

      ELSE                              ! other processors can recieve the matrix
         CALL ZGEBR2D(CONTEXT,'A',TOP,N,N,H,N,0,0)

      END IF
____________________________________________________________________________________________________   

      ! let each processor initialize its local array elements with the correct 
      ! values from the distributed matrix
      DO i_loc = 1, N
         DO j_loc = 1,N
            CALL pzelset(LH,i_loc,j_loc,DESCH, H(i_loc,j_loc))
         END DO
      END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Ask PCHEEVX to compute the entire eigendecomposition
       
       !this setting yields the most orthogonal eigenvectors
       ABSTOL = PDLAMCH(CONTEXT, 'U')

       !ORFAC specifies which eigenvectors should be reorthogonalized
       !set to the default value
       ORFAC = 0.001

       !allocate arrays used in PCHEEVX
        ALLOCATE(IFAIL(N), ICLUSTR(2*NPROW*NPCOL), GAP(NPROW*NPCOL))

       ! variables used for allocating different workspaces 
       NNP = MAX(N, NPROW*NPCOL + 1, 4)
       NN = MAX(N,NBROW, 2)
       MQO = NUMROC(NN, NBCOL, 0, 0, NPCOL)
       NPO = NUMROC(NN, NBROW, 0, 0,NPROW)

       ! first carry out a workspace query
       lwork = -1
       lrwork = -1
       liwork = -1
       ALLOCATE(work(1))
       ALLOCATE(rwork(1))
       ALLOCATE(iwork(1))
       CALL PZHEEVX( 'V', 'A', 'U', N, LH, 1, 1, DESCH, VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, LZ, 1, 1, DESCZ,&
            WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO ) 
 
       lwork = INT(ABS(work(1)))
       DEALLOCATE(work)
       ALLOCATE(work(lwork))
      
       lrwork = INT(ABS(rwork(1)))
       DEALLOCATE(rwork)
       ALLOCATE(rwork(lrwork))

       liwork =INT (ABS(iwork(1)))
       DEALLOCATE(iwork)
       ALLOCATE(iwork(liwork))

      CALL PZHEEVX( 'V', 'A', 'U', N, LH, 1, 1, DESCH, VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, LZ, 1, 1, DESCZ,&
            WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
!_________________________________________________________________________________________

! collect all the orthonormal eigenvectors of the matrix

!Begin master part

IF (rank.EQ.0) THEN
   !obtain dimensions of the local array on master node
   LOCr = NUMROC(N,NBROW,MYROW,0,NPROW)
   LOCc = NUMROC(N,NBCOL,MYCOL,0,NPCOL)

   !k is the rank of the processor that sent the data
   !nprocs is the number of processors in the communicator
   DO node = 0, nprocs -1 
      IF (node.EQ.0) THEN
         irow = myrow
         icol = mycol
      ELSE
         
         !node -- is the rank of the processor that sent the data
         !status -- message status array(of type integer)
         CALL MPI_Recv(irow, 1, MPI_Integer, NODE, 10, MPI_COMM_WORLD, ISTATUS, IERR)
         CALL MPI_Recv(icol, 1, MPI_Integer, NODE, 20, MPI_COMM_WORLD, ISTATUS, IERR)

         CALL MPI_Recv(LOCr, 1, MPI_Integer, NODE, 40, MPI_COMM_WORLD, ISTATUS, IERR)
         CALL MPI_Recv(LOCc, 1, MPI_Integer, NODE, 50, MPI_COMM_WORLD, ISTATUS,IERR )
         
         !create a local matrix with the size of the matrix passed
          DEALLOCATE(LZ)
          ALLOCATE( LZ(locr,locc))
       
          !recieve the local matrix sent from node
          CALL MPI_Recv(LZ, locr*locc, MPI_DOUBLE_COMPLEX, NODE, 30,MPI_COMM_WORLD, ISTATUS, IERR)
      ENDIF

      !compute the number of blocks in each local array
      nbrows = CEILING(REAL(LOCr)/REAL(NBROW))    !number of blocks in the row
      nbcolumns = CEILING(REAL(LOCc)/REAL(NBCOL))    !number of blocks in the columns

     !loop over each block in the column
      DO ibcolumns = 1, nbcolumns  
         !special case - number of columns is less than NBCOL     
         IF(ibcolumns.EQ.nbcolumns) THEN
            col_max = locc-(nbcolumns-1)*NBCOL 
         ELSE
            col_max = NBCOL !number of columns in a block        
         ENDIF

         !loop over each block in the row
         DO ibrows = 1, nbrows
            !special case - number of columns is less than NBCO
            IF (ibrows.EQ.nbrows) THEN
               row_max = locr-(nbrows-1)*NBROW
            ELSE
               row_max = NBROW !number of rows in a block
            ENDIF
            !for each column in the block, loop over each row in the block
            DO columns = 1, col_max
               DO rows = 1, row_max
                  Z((irow+(ibrows-1)*nprow)*NBROW + rows, (icol+(ibcolumns-1)*npcol)*NBCOL + columns)&
                        = LZ((ibrows-1)*NBROW + rows, (ibcolumns-1)*NBCOL + columns)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   !End master part.

   ! Begin slave part.
ELSE 
   !send the processor coordinates in grid
   CALL MPI_Send(myrow, 1, MPI_Integer, 0, 10, MPI_COMM_WORLD, IERR)
   CALL MPI_Send(mycol, 1, MPI_Integer, 0, 20, MPI_COMM_WORLD, IERR)
   !send number rows and columns
   CALL MPI_Send(LDA, 1, MPI_Integer, 0, 40, MPI_COMM_WORLD, IERR)
   CALL MPI_Send(LDB, 1, MPI_Integer, 0, 50, MPI_COMM_WORLD, IERR)
   !send the local matrix
   CALL MPI_Send(LZ,LDA*LDB , MPI_DOUBLE_COMPLEX, 0, 30, MPI_COMM_WORLD, IERR)

ENDIF
! End slave part.
!------------------------------------------------------------------------------------------

OPEN(17,file='differences.dat')

IF(myrow.EQ.0.AND.mycol.EQ.0) THEN
!OPEN(17,file='differences.dat')
 DO i1=1,N
  diff = SUM(ABS(MATMUL(H,Z(:,i1)) - W(i1)*(Z(:,i1))))
 WRITE(17,101) 'difference = ', diff
END DO
 
    ! changed H to A and E to W
    ! Now get the density of states
    CALL total_dos(N,Z,W,dos)

    OPEN(15,file='test.dat')
    DO i1 = 1,dos%nbin
       WRITE(15,101) dos%wmin+(i1-1)*dos%dw , dos%total_dos(i1)
    END DO
   101 FORMAT(2(1x,g13.6))
  END IF 

  DEALLOCATE(H,U,E,dos%total_dos)
  DEALLOCATE(LH,LZ)
!_________________________________________________________________________________________

	CALL BLACS_GRIDEXIT( CONTEXT )
	CALL BLACS_EXIT( 0 )

END PROGRAM PZHEEVX
