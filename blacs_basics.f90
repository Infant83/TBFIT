#include "alias.inc"
module blacs_basics
    implicit none
    logical :: flag_use_blacs
    integer :: mpi_comm_blacs
    integer :: myid_blacs
    integer :: nprow, npcol, myrow, mycol

endmodule blacs_basics
