#include "alias.inc"
module print_io
   character*2048    message
   character*2048    message1
   character*2048    message2
   character*2048    message3
   character*2048    message4
   character*2048    message5
   character*2048    message6
   character*2048    message7
   character*2048    message8
   character*2048    message9
   character*2048    message_pack(1024) ! just assumed that maximum number of cpus = 1024 (one can increase according to your system))
   integer*4, public, parameter :: pid_log = 17
   integer*4, public            :: iverbose     ! verbosity 1: full, 2: no
   integer*4, public            :: print_mode   ! this is used to control verbosity when 
                                                ! calling write_log subroutine
                                                ! In some subroutines, one needs to write only
                                                ! to the output log file, not to screen.
                                                ! In this case, one use "print_mode" instead of
                                                ! "iverbose" parameter when calling write_log
                                                ! subroutine. 
contains        
   subroutine open_log(fnamelog, myid)
      implicit none
      character(len=*)     :: fnamelog
      integer*4,intent(in) :: myid

      if_main open(pid_log, file=trim(fnamelog), status='unknown')

      return
   endsubroutine

   subroutine write_log(string, imode, myid)
      implicit none
      character(len=*)      :: string
      integer*4                imode
      integer*4, intent(in) :: myid

      if(iverbose .eq. 2) return

      select case (imode)

        case(1)
          if_main write(pid_log, '(A)') trim(string)
    
        case(2)
          if_main write(6,       '(A)') trim(string)
   
        case(3)
          if_main write(pid_log, '(A)') trim(string)
          if_main write(6,       '(A)') trim(string)
       
        case(11)
                  write(pid_log, '(A)') trim(string)

        case(12)
                  write(6,       '(A)') trim(string)

        case(13)
                  write(pid_log, '(A)') trim(string)
                  write(6,       '(A)') trim(string)
        
        case(21)  ! write with advance=no to file
          if_main write(pid_log, '(A)', ADVANCE='NO') trim(string)

        case(22)  ! write with advance=no to file
          if_main write(6      , '(A)', ADVANCE='NO') trim(string)

        case(23)  ! write with advance=no to file
          if_main write(pid_log, '(A)', ADVANCE='NO') trim(string)
          if_main write(6      , '(A)', ADVANCE='NO') trim(string)

        case(99) ! return without writing anything
          return

        case(97) 
          return
         !if_main write(pid_log, '(A)') trim(string)

      end select

      return
   endsubroutine

   subroutine close_log(myid)
      implicit none
      integer*4, intent(in) :: myid

      if_main close(pid_log)

      return
   endsubroutine

endmodule
