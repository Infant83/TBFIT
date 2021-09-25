! mpif90 -o para_feast para_feast.f90 -L../LIB/FEAST/4.0/lib/x64 -lpfeast -mkl=parallel
program para_feast
  implicit none
  include 'mpif.h'
  integer :: n,nnz
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa,c
  integer,dimension(:),allocatable :: isa,jsa,ic,jc

  integer, dimension(64)     :: fpm
  integer                    :: M0, M, i, info ! ne_max
  double precision           :: Emin = -0.35d0, Emax = 0.23d0
  double precision           :: rea,img
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: X
  double precision, dimension(:), allocatable :: E, res
  double precision               :: epsout
  integer                        :: nL3, rank, loop
  real t1, t0
  integer :: code
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

  nL3=1 ; M0 = 40
  open(10,file='system2.mtx',status='old')
  read(10,*) n,n,nnz
  allocate(ic(nnz))
  allocate(jc(nnz))
  allocate(c(nnz))
     do i=1,nnz
        read(10,*) ic(i),jc(i),rea,img
        c(i)=rea*(1.0d0,0.0d0)+img*(0.0d0,1.0d0)
     end do
     close(10)

  allocate(isa(1:n+1))
  allocate(jsa(1:nnz))
  allocate(sa(1:nnz))
! if(rank .eq. 0) then
  call zcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
! endif
  deallocate(ic,jc,c)
  allocate(E(M0))
  allocate( res(M0))
! if(rank .eq. 0) then
  allocate(X(n, M0))
! endif

  call cpu_time(t0)
  call pfeastinit(fpm, MPI_COMM_WORLD, nL3)
  fpm(1) = 1
  write(6,*)"L3 comm", fpm(49), rank
  write(6,*)"L2 comm", fpm( 9), rank
  call pzfeast_hcsrev('F',  n,    sa,    isa,  jsa,   fpm,epsout, loop,  Emin,Emax ,M0 ,E ,X ,M, res ,info)
  call cpu_time(t1)
  if(rank == 2) then
    write(6,*)"mpi time: ", t1 -t0
    write(6,*)' MPI Solutions(E, EV, res) at rank L2'
    do i = 1, M
      write(6,*)' E=',E(i), 'res=', res(i)
    enddo
  endif

! if(rank == 0) then
!   call cpu_time(t0)
!   call feastinit(fpm)
!   fpm(1) = 0
!   call zfeast_hcsrev('F',  n,    sa,    isa,  jsa,   fpm,epsout, loop,  Emin,Emax ,M0 ,E ,X ,M, res ,info)
!   call cpu_time(t1)
!   write(6,*)"serial time: ", t1 -t0

!   write(6,*)' Serial Solutions(E, EV, res) at rank L2'
!   do i = 1, M
!     write(6,*)' E=',E(i), 'res=', res(i)
!   enddo
! endif

  call MPI_FINALIZE(code)
end program para_feast
