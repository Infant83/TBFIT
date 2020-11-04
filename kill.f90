#include "alias.inc"
module kill
  use mpi_setup
  use print_io
  implicit none

  contains

    subroutine check_kill_tbfit(PINPT,PPRAM, PWGHT)
      use parameters, only: incar, weight, params
      logical        flag_exist
      integer*4      mpierr
      type(incar)                            :: PINPT
      type(params)                           :: PPRAM
      type(weight), dimension(PINPT%nsystem) :: PWGHT

      inquire(file="KILLFIT",exist=flag_exist)

#ifdef MPI
      call MPI_BCAST(flag_exist, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
#endif

      if(flag_exist) then

        if( PPRAM%slater_koster_type .gt. 10) then
          call update_param_nrl( PPRAM )
        else
          call update_param( PPRAM )
        endif

        if_main call execute_command_line('\rm -f KILLFIT')

        ! NOTE: only PWGHT(1) info is printed along with
        !       even though PINPT%nsystem > 1, this is because we just consider PPRAM to be applied
        !       for entire systems in fitting process. 
        !       To avoid any confusion, in the future release, this could be 
        !       corrected or other approach in printing PPRAM can be considered, but in this version
        !       we just keep this stratege for the convenience... 31.Oct.2020 HJK
        if_main call print_param (PINPT, PPRAM, PWGHT(1), PPRAM%pfileoutnm, .TRUE.)
        write(message,'(A)') ' '  ; write_msg
        write(message,'(A)') ' Termination of job is requested by providing KILLFIT file.' ; write_msg
        write(message,'(2A)') ' The latest updates of PARAMETERS will be written in ', trim(PPRAM%pfileoutnm) ; write_msg
        write(message,'(A)') ' Kill the job now...' ; write_msg
        write(message,'(A)') ' '  ; write_msg

        kill_job

      endif
      return
    endsubroutine

endmodule
