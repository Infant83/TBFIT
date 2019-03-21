module mpi_basics
   implicit none
   include 'mpif.h'
   logical :: flag_use_mpi
   integer :: npar
   integer :: kpar
   integer :: nproc_per_band
   integer :: mpi_comm_earth
   integer :: myid
   integer :: nprocs

   type mpicomm 
        integer :: mpi_comm
        integer :: myid
        integer :: nprocs
   endtype

endmodule mpi_basics
