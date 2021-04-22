! Following kmeans_01 subroutine is from https://people.sc.fsu.edu/~jburkardt/f_src/kmeans/kmeans.html
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
!
!  Idiot checks.
!
  if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop 1
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop 1
  end if

  if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop 1
  end if
!
!  For each observation, calculate the distance from each cluster
!  center, and assign to the nearest.
!
  do i = 1, point_num

    do j = 1, cluster_num
      cluster_energy(j) = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do

    list = minloc ( cluster_energy(1:cluster_num) )
    cluster(i) = list(1)

  end do
!

!
!  Determine the cluster population counts.
!
  cluster_population(1:cluster_num) = 0

  do i = 1, point_num
    j = cluster(i)
    cluster_population(j) = cluster_population(j) + 1
  end do
!
!  Calculate the mean and sum of squares for each cluster.
!
  cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
      + point(1:dim_num,i)
  end do

  do i = 1, cluster_num
    if ( 0 < cluster_population(i) ) then
      cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
        / real ( cluster_population(i), kind = 8 )
    end if
  end do
!
!  Set the point energies.
!
  f(1:point_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
  end do
!
!  Set the cluster energies.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + f(i)
  end do
!
!  Adjust the point energies by a weight factor.
!
  do i = 1, point_num
    j = cluster(i)
    if ( 1 < cluster_population(j) ) then
      f(i) = f(i) * real ( cluster_population(j), kind = 8 ) &
        / real ( cluster_population(j) - 1, kind = 8 )
    end if
  end do
!
!  Examine each observation in turn to see if it should be
!  reassigned to a different cluster.
!

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

          de = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) &
            * real ( cluster_population(j), kind = 8 ) &
            / real ( cluster_population(j) + 1, kind = 8 )

          if ( de < dc ) then
            dc = de
            ir = j
          end if

        end if

      end do
!
!  If the lowest value was obtained by staying in the current cluster,
!  then cycle.
!
      if ( ir == il ) then
        cycle
      end if
!
!  Reassign the point from cluster IL to cluster IR.
!
      cluster_center(1:dim_num,il) = &
        ( cluster_center(1:dim_num,il) &
        * real ( cluster_population(il), kind = 8 ) &
        - point(1:dim_num,i) ) / real ( cluster_population(il) - 1, kind = 8 )

      cluster_center(1:dim_num,ir) = &
        ( cluster_center(1:dim_num,ir) &
        * real ( cluster_population(ir), kind = 8 ) &
        + point(1:dim_num,i) ) / real ( cluster_population(ir) + 1, kind = 8 )

      cluster_energy(il) = cluster_energy(il) - f(i)
      cluster_energy(ir) = cluster_energy(ir) + dc
      cluster_population(ir) = cluster_population(ir) + 1
      cluster_population(il) = cluster_population(il) - 1

      cluster(i) = ir
!
!  Adjust the value of F for points in clusters IL and IR.
!

      do j = 1, point_num

        k = cluster(j)

        if ( k == il .or. k == ir ) then

          f(j) = sum ( &
           ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )

          if ( 1 < cluster_population(k) ) then
            f(j) = f(j) * real ( cluster_population(k), kind = 8 ) &
              / ( real ( cluster_population(k) - 1, kind = 8 ) )
          end if

        end if

      end do

      swap = swap + 1

    end do
!
!  Exit if no reassignments were made during this iteration.
!
    if ( swap == 0 ) then
      exit
    end if

  end do
!
!  Compute the cluster energies.
!
  cluster_energy(1:cluster_num) = 0.0D+00

  do i = 1, point_num

    j = cluster(i)

    cluster_energy(j) = cluster_energy(j) + sum ( &
      ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

  end do

  return
end

