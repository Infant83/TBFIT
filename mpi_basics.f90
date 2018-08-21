module mpi_basics
   implicit none
   include 'mpif.h'
   logical :: flag_use_mpi
   integer :: nprocs
   integer :: myid
   integer :: mpi_comm_earth
!  integer, parameter :: mpi_int4 = MPI_INTEGER4
!  integer, parameter :: mpi_real8  = MPI_REAL8
!  integer, parameter :: mpi_comp16 = MPI_COMPLEX16

endmodule mpi_basics
