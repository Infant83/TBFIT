#include "alias.inc"
module version
  use time
  use mpi_setup
  use print_io 

contains 

  function version_tag() result(ver_tag)
    implicit none
    character*132  ver_tag
    character*132  dummy, version, inputline
    integer*4      i_continue, idx1
    integer*4      mpierr

    i_continue = 0

    if(myid .eq. 0) then

        call execute_command_line('grep TBFIT_VERSION makefile > _foo_version')
        open(11, file='_foo_version', status='unknown')
        do while (i_continue .ge. 0)
            read(11, '(A)', iostat=i_continue) inputline
            if(i_continue .lt. 0) exit
            inputline = adjustl(trim(inputline))

            if(index(inputline,'#') .ne. 1 .and. &
               index(inputline, 'TBFIT_VERSION') .eq. 1) then
                idx1 = index(inputline, '=')
                version = inputline(idx1+1:132)
                read(version,*) version
            endif
        enddo
        close(11)
        call execute_command_line('rm -f _foo_version')
    endif
    ver_tag='# TBFIT version '//adjustl(trim(version))//' (build: ' // __DATE__// ' ' //__TIME__// ') '
    
#ifdef MPI
    call MPI_BCAST(ver_tag, 132, MPI_CHARACTER, 0 , mpi_comm_earth, mpierr)
#endif

    return
  endfunction

  subroutine version_stamp(t_start)
    implicit none
   !character*132  ver !, version_tag
    real*8         t_start
    
   !ver_tag='# TBFIT 2019. Apr. 23. (build ' // __DATE__// ' ' //__TIME__// ') '
   !ver_tag='# TBFIT 2020. Jun. 25. (build ' // __DATE__// ' ' //__TIME__// ') '
   !ver_tag='# TBFIT version 0.4.1  (build ' // __DATE__// ' ' //__TIME__// ') '
   !ver_tag='# TBFIT version 0.5.1  (build ' // __DATE__// ' ' //__TIME__// ') '
   !ver = verion_tag()

    call write_log(trim( version_tag()          ),3,myid)
    call write_log(' ',3,myid)
    call write_log("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",3,myid)
    call write_log("#! Copyright 2018. Hyun-Jung Kim All rights reserved.!",3,myid)
    call write_log("#!           (FJZ,  h.kim@fz-juelich.de)             !",3,myid)
    call write_log("#!           (KIAS, Infant@kias.re.kr)               !",3,myid)
    call write_log("#!           (angpangmokjang@hanmail.net)            !",3,myid)
    call write_log("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",3,myid)
    call write_log(' ',3,myid)

   !if(t_start .gt. 0d0) then
      call write_log(" -------------------------------------------------------------------------",3,myid)
      call timestamp('| Program start on', t_start)
      call report_hostname()
      call write_log(" -------------------------------------------------------------------------",3,myid)
   !endif

    return
  endsubroutine

endmodule
