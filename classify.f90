module classify

contains

! Following kmeans subroutine is from https://rosettacode.org/wiki/K-means%2B%2B_clustering#Fortran 
!**********************************************************************
! KMPP - K-Means++ - Traditional data clustering with a special initialization
! Public Domain - This program may be used by any person for any purpose.
!
! Origin:
!    Hugo Steinhaus, 1956
!
! Refer to:
!    "kmeans++: the advantages of careful seeding"
!    David Arthur and Sergei Vassilvitskii
!    Proceedings of the eighteenth annual ACM-SIAM symposium 
!      on Discrete algorithms, 2007
!
!____Variable_______I/O_______Description___________________Type_______
!    X(P,N)         In        Data points                   Real
!    P              In        Dimension of the data         Integer
!    N              In        Number of points              Integer
!    K              In        # clusters                    Integer
!    C(P,K)         Out       Center points of clusters     Real
!    Z(N)           Out       What cluster a point is in    Integer
!    WORK(N)        Neither                                 Real
!    IFAULT         Out       Error code                    Integer
!***********************************************************************
subroutine kmeanspp (X, P, N, K, C, Z, WORK, IFAULT)
       IMPLICIT NONE
       INTEGER P, N, K, Z(N), IFAULT
       REAL*8 X(P,N), C(P,K), WORK(N)
 
       INTEGER ITER                 ! maximum iterations
       REAL BIG                     ! arbitrary large number
       PARAMETER (ITER = 1000,BIG = 1E33)
       INTEGER H           ! count iterations
       INTEGER I           ! count points
       INTEGER I1          ! point marked as initial center
       INTEGER J           ! count dimensions
       INTEGER L           ! count clusters
       INTEGER L0          ! present cluster ID
       INTEGER L1          ! new cluster ID
       REAL*8  BEST        ! shortest distance to a center
       REAL*8  D2          ! squared distance
       REAL*8  TOT         ! a total
       REAL*8  W           ! temp scalar
       LOGICAL CHANGE      ! whether any points have been reassigned

       IF (K < 1 .OR. K > N) THEN       ! K out of bounds
         IFAULT = 3
         RETURN
       END IF
 
!      Begin.
       IFAULT = 0 ; Z = 0 ; WORK = BIG
 
       ! choose initial center from data points
       CALL RANDOM_NUMBER (W)
       I1 = MIN(INT(W * FLOAT(N)) + 1, N)  ! choose first center at random
       C(:,1) = X(:,I1)
 
       ! initialize other centers
       DO L = 2, K              
         TOT = 0.
         DO I = 1, N                     ! measure from each point
           BEST = WORK(I)
           D2 = 0.                         ! to prior center
           DO J = 1, P
             D2 = D2 + (X(J,I) - C(J,L-1)) **2  ! Squared Euclidean distance
             IF (D2 .GE. BEST) GO TO 10               ! needless to add to D2
           END DO                          ! next J
           IF (D2 < BEST) BEST = D2          ! shortest squared distance 
           WORK(I) = BEST 
  10       TOT = TOT + BEST             ! cumulative squared distance
         END DO                      ! next data point
 
         ! Choose center with probability proportional to its squared distance from existing centers.
         CALL RANDOM_NUMBER (W)
         W = W * TOT    ! uniform at random over cumulative distance
         TOT = 0.
         DO I = 1, N
           I1 = I
           TOT = TOT + WORK(I)
           IF (TOT > W) GO TO 20
         END DO                ! next I
  20     CONTINUE
         C(:,L) = X(:,I1)
       END DO   

!      main loop for k-means
       DO H = 1, ITER
         CHANGE = .FALSE.
 
!        find nearest center for each point 
         DO I = 1, N
           L0 = Z(I)
           L1 = 0
           BEST = BIG
           DO L = 1, K
             D2 = 0.
             DO J = 1, P
               D2 = D2 + (X(J,I) - C(J,L)) **2
               IF (D2 .GE. BEST) GO TO 30
             END DO
  30         CONTINUE
             IF (D2 < BEST) THEN           ! new nearest center
               BEST = D2
               L1 = L
             END IF             
           END DO        ! next L
 
           IF (L0 .NE. L1) THEN
             Z(I) = L1                   !  reassign point 
             CHANGE = .TRUE.
           END IF
         END DO         ! next I
         IF (.NOT. CHANGE) RETURN      ! success
 
!        find cluster centers
         WORK(1:K) = 0.
         C(1:P,1:K) = 0.
 
         DO I = 1, N
           L = Z(I)
           WORK(L) = WORK(L) + 1.             ! count
           DO J = 1, P
             C(J,L) = C(J,L) + X(J,I)         ! add
           END DO
         END DO
 
         DO L = 1, K
           IF (WORK(L) < 0.5) THEN          ! empty cluster check
             IFAULT = 1                     ! fatal error
             RETURN
           END IF
           W = 1. / WORK(L)
           DO J = 1, P
             C(J,L) = C(J,L) * W     ! multiplication is faster than division
           END DO
         END DO
 
       END DO                   ! next H
       IFAULT = 2                ! too many iterations

       RETURN
endsubroutine  ! of KMPP

! Following kmeans subroutine is from https://people.sc.fsu.edu/~jburkardt/f_src/kmeans/kmeans.html
! KMEANS, a FORTRAN90 code which handles the K-Means problem, which organizes a set of N points in M dimensions into K clusters;
subroutine kmeans ( dim_num, point_num, cluster_num, it_max, it_num, point, &
                    cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_01 applies the K-Means algorithm.
!
!  Discussion:
!
!    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
!    observations are to be allocated to CLUSTER_NUM clusters in such 
!    a way that the within-cluster sum of squares is minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    Original FORTRAN77 version by David Sparks.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58: 
!    Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates which cluster
!    each point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    cluster energies.
!

  implicit none

  integer ( kind = 4 ) cluster_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) cluster(point_num)
  real ( kind = 8 ) cluster_center(dim_num,cluster_num)
  real ( kind = 8 ) cluster_energy(cluster_num)
  integer ( kind = 4 ) cluster_population(cluster_num)
  real ( kind = 8 ) dc
  real ( kind = 8 ) de
  real ( kind = 8 ) f(point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(1)
  real ( kind = 8 ) point(dim_num,point_num)
  integer ( kind = 4 ) swap

  it_num = 0
  !  Idiot checks.
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) 'KMEANS - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop 1
  end if
  if ( dim_num < 1 ) then
    write ( *, '(a)' ) 'KMEANS - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop 1
  end if
  if ( point_num < 1 ) then
    write ( *, '(a)' ) 'KMEANS - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop 1
  end if

  !  For each observation, calculate the distance from each cluster center, and assign to the nearest.
  do i = 1, point_num
    do j = 1, cluster_num
      cluster_energy(j) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do
    list = minloc ( cluster_energy(1:cluster_num) )
    cluster(i) = list(1)
  end do

  !  Determine the cluster population counts.
  cluster_population(1:cluster_num) = 0
  do i = 1, point_num
    j = cluster(i)
    cluster_population(j) = cluster_population(j) + 1
  end do

  !  Calculate the mean and sum of squares for each cluster.
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00
  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) + point(1:dim_num,i)
  end do
  do i = 1, cluster_num
    if ( 0 < cluster_population(i) ) then
      cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) / real ( cluster_population(i), kind = 8 )
    end if
  end do

  !  Set the point energies.
  f(1:point_num) = 0.0D+00
  do i = 1, point_num
    j = cluster(i)
    f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
  end do

  !  Set the cluster energies.
  cluster_energy(1:cluster_num) = 0.0D+00
  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + f(i)
  end do

  !  Adjust the point energies by a weight factor.
  do i = 1, point_num
    j = cluster(i)
    if ( 1 < cluster_population(j) ) then
      f(i) = f(i) * real ( cluster_population(j))/ real ( cluster_population(j) - 1)
    end if
  end do

  !  Examine each observation in turn to see if it should be reassigned to a different cluster.
  it_num = 0
  do while ( it_num < it_max )
    it_num = it_num + 1
    swap = 0
    do i = 1, point_num
      il = cluster(i)
      ir = il
      if ( cluster_population(il) <= 1 ) then
        cycle
      end if
      dc = f(i)
      do j = 1, cluster_num
        if ( j /= il ) then
          de = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) * real(cluster_population(j)) / real ( cluster_population(j) + 1)
          if ( de < dc ) then
            dc = de
            ir = j
          end if
        end if
      end do

      if ( ir == il ) cycle !  If the lowest value was obtained by staying in the current cluster, then cycle.

      !  Reassign the point from cluster IL to cluster IR.
      cluster_center(1:dim_num,il) =( cluster_center(1:dim_num,il) * real ( cluster_population(il)) &
                                    - point(1:dim_num,i) ) / real ( cluster_population(il) - 1)

      cluster_center(1:dim_num,ir) = ( cluster_center(1:dim_num,ir) * real ( cluster_population(ir)) &
                                    + point(1:dim_num,i) ) / real ( cluster_population(ir) + 1)

      cluster_energy(il) = cluster_energy(il) - f(i)
      cluster_energy(ir) = cluster_energy(ir) + dc
      cluster_population(ir) = cluster_population(ir) + 1
      cluster_population(il) = cluster_population(il) - 1

      cluster(i) = ir

      !  Adjust the value of F for points in clusters IL and IR.
      do j = 1, point_num
        k = cluster(j)
        if ( k == il .or. k == ir ) then
          f(j) = sum ( ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )
          if ( 1 < cluster_population(k) ) then
            f(j) = f(j) * real ( cluster_population(k))/ ( real ( cluster_population(k) - 1) )
          end if
        end if
      end do
      swap = swap + 1
    end do
    !  Exit if no reassignments were made during this iteration.
    if ( swap == 0 ) exit
  end do

  !  Compute the cluster energies.
  cluster_energy(1:cluster_num) = 0.0D+00
  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
  end do

  return
endsubroutine

endmodule

