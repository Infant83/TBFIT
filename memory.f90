#include "alias.inc"
module memory
   use print_io
   use mpi_setup
   implicit none

contains

  subroutine report_memory(sizeof_, length, arg_name)
    integer*8    sizeof_
    integer*4    length
    character(len=*), optional :: arg_name

    if(present(arg_name)) then
      write(message,'(3A,F10.3,A)')'   MEMORY USAGE (',arg_name,'): ', real(sizeof_)*real(length)/1024d0/1024d0/1024d0,' GB'
    else
      write(message,'( A,F10.3,A)')'   MEMORY USAGE : ',               real(sizeof_)*real(length)/1024d0 /1024d0 /1024d0,' GB'
    endif
    write_msg

    return
  endsubroutine

  subroutine report_memory_total(ispinor, ispin, nspin, neig, nband, nkp, &
                                 flag_stat, flag_sparse, ncpus)
    integer*4   ispinor, ispin, nspin, neig, nband, nkp
    integer*8   size_E, size_V, size_Hk, size_Hm, size_Hs
    integer*4   ncpus, impi
    logical     flag_stat, flag_sparse

#ifdef MPI
    impi = 2 ! due to MPI_REDUCE requires additional memory
#else
    impi = 1 
#endif

!   size_V = neig*ispin  *nband*nspin  * nkp * impi
    size_E =              int8(nband)*int8(nspin)  * int8(nkp) * int8(ncpus)
    size_V =int8(neig)*int8(ispin  )*int8(nband)*int8(nspin  )* int8(nkp  ) + &  ! EE%V 
            int8(neig)*int8(ispin  )*int8(nband)*int8(nspin  )* int8(nkp  )    !    V in root node
    size_Hk=int8(neig)*int8(ispinor)*int8(neig )*int8(ispinor)* int8(ncpus)  
    size_Hm=int8(neig)*int8(ispinor)*int8(neig )*int8(ispinor)* int8(ncpus)
    size_Hs=int8(neig)*int8(ispinor)*int8(neig )*int8(ispinor)* int8(ncpus)

    if(flag_stat) then
      call report_memory(size_E ,8 , 'Eigen values')
      call report_memory(size_V ,16, 'Eigen vector')
      if(.not. flag_sparse) then
        call report_memory(size_Hk,16, 'H(total)    ')
        call report_memory(size_Hm,16, 'H(magnetic) ')
        call report_memory(size_Hs,16, 'H(spin-orb) ')
      endif
    endif
  
    return
  endsubroutine

endmodule
