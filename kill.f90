#include "alias.inc"
module kill
  use mpi_setup
  use print_io
  implicit none

  contains

    subroutine check_kill_tbfit(PINPT)
      use parameters, only: incar
      logical        flag_exist
      integer*4      mpierr
      type(incar) :: PINPT

      inquire(file="KILLFIT",exist=flag_exist)

#ifdef MPI
      call MPI_BCAST(flag_exist, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
#endif

      if(flag_exist) then

        if( PINPT%slater_koster_type .gt. 10) then
          call update_param_nrl( PINPT )
        else
          call update_param( PINPT )
        endif

        if_main call execute_command_line('\rm -f KILLFIT')

        if_main call print_param (PINPT, PINPT%pfileoutnm, PINPT%flag_print_param)
        write(message,'(A)') ' '  ; write_msg
        write(message,'(A)') ' Termination of job is requested by providing KILLFIT file.' ; write_msg
        write(message,'(2A)') ' The latest updates of PARAMETERS will be written in ', trim(PINPT%pfileoutnm) ; write_msg
        write(message,'(A)') ' Kill the job now...' ; write_msg
        write(message,'(A)') ' '  ; write_msg

        kill_job

      endif
      return
    endsubroutine

endmodule
