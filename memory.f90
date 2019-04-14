#include "alias.inc"
module memory
   implicit none

contains

  subroutine report_memory(sizeof_, length, arg_name)
    integer*4    kind
    integer*4    sizeof_, length
    character(len=*), optional :: arg_name

    if(present(arg_name)) then
      write(6,'(3A,F10.3,A)')'  MEMORY USAGE (',arg_name,'): ', &
                        real(sizeof_)*real(length) / 1024d0 / 1024d0 / 1024d0,' GB'
    else
      write(6,'( A,F10.3,A)')'  MEMORY USAGE : ',&
                        real(sizeof_)*real(length) / 1024d0 / 1024d0 / 1024d0,' GB'
    endif

    return
  endsubroutine

  subroutine report_memory_total(ispinor, ispin, nspin, neig, nband, nkp, &
                                 flag_stat, flag_sparse, flag_use_mpi, ncpus)
    integer*4   ispinor, ispin, nspin, neig, nband, nkp
    integer*4   size_E, size_V, size_Hk, size_Hm, size_Hs
    integer*4   ncpus, impi
    logical     flag_stat, flag_sparse, flag_use_mpi

    impi = 1 

    if(flag_use_mpi) impi = 2 ! due to MPI_REDUCE requires additional memory

!   size_V = neig*ispin  *nband*nspin  * nkp * impi
    size_E =              nband*nspin  * nkp * ncpus
    size_V = neig*ispin  *nband*nspin  * nkp + &  ! EE%V 
             neig*ispin  *nband*nspin  * nkp      !    V in root node
    size_Hk= neig*ispinor*neig *ispinor* ncpus  
    size_Hm= neig*ispinor*neig *ispinor* ncpus
    size_Hs= neig*ispinor*neig *ispinor* ncpus

    if(flag_stat) then
      call report_memory(size_E ,8 , 'Eigen values')
      call report_memory(size_V ,16, 'Eigen vector')
      call report_memory(size_Hk,16, 'H(total)    ')
      if(.not. flag_sparse) then
        call report_memory(size_Hm,16, 'H(magnetic) ')
        call report_memory(size_Hs,16, 'H(spin-orb) ')
      endif
    endif
  
    return
  endsubroutine

endmodule
