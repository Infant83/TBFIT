#include "alias.inc"
subroutine test()
   use mpi_setup
   implicit none
   integer*4    mpierr

   integer*4          pid, msize_sub
   external           blacs_pnum  ! external blacs function
   external           blacs_get, blacs_gridinit, blacs_gridinfo   ! external blacs subroutine
   integer*4          blacs_pnum
   integer*4          ourjob(0)
!  integer*4          myid_blacs
!  integer*4          nprow, npcol, myrow, mycol, icntxt

#ifdef MPI
   call MPI_BARRIER(mpi_comm_earth, mpierr)
   call MPI_BARRIER(mpi_comm_earth, mpierr)
#ifdef SCALAPACK

!  ourjob(0) = 999 
!  write(6,*)"VVVV "

!  write(6,*)"VVVV ", ourjob(0)
  
!  stop
!  call blacs_initialize()

!  pid = 66
!  msize_sub = 6
!  nprow = 2
!  npcol = 2



!  call blacs_get(-1,0,mpi_comm_blacs)
!  call blacs_gridinit(mpi_comm_blacs, 'Row-major', nprow, npcol)
!  call blacs_gridinfo(mpi_comm_blacs, nprow, npcol, myrow, mycol); if(mpi_comm_blacs .eq. -1) go to 99
!  myid_blacs= blacs_pnum(mpi_comm_blacs, myrow, mycol)
!  call blacs_pcoord(icntxt, myid, myrow, mycol)
   write(6,*)"XXX ", myid, myid_blacs,  myrow, mycol
#endif

   ourjob(0)=987
   write(6,*)"BBBB ", ourjob(0)

99 call MPI_BARRIER(mpi_comm_earth, mpierr)
   call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif
   kill_job
return
endsubroutine
