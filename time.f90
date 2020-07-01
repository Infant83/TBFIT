#include "alias.inc"
module time
   use print_io
contains

  subroutine time_check(time_elapsed, time_old, init)
    use mpi_setup
    implicit none
    real*8   time_elapsed
    real*8   time_now
    real*8   time_old
    logical  flag_init
    character(len=*), optional :: init

    ! How much time has been elapsed since time_old?
    ! If time_now  = time_old : report current time
    ! If time_now /= time_old : report time_now - time_old 
    !                           report time_old as current time

    if(present(init) .and. trim(init) .eq. 'init') then
      time_elapsed = 0d0
      time_old     = 0d0
    endif

    if(time_elapsed .eq. time_old) then
#ifdef MPI
      time_now=MPI_Wtime()
#else
      call cpu_time(time_now)
#endif
      time_elapsed = 0d0
      time_old = time_now
    elseif(time_elapsed .ne. time_old) then
#ifdef MPI
      time_now=MPI_Wtime()
#else
      call cpu_time(time_now)
#endif
      time_elapsed = time_now - time_old
      time_old = time_now
    endif

    return
  endsubroutine

subroutine timestamp (time_step,t)
! Note: This subroutine "timestamp" is copied from MINPACK library of John Burkardt
!  Licensing: This code is distributed under the GNU LGPL license.
!  Modified: 18 May 2013
!  Author: John Burkardt
!  Parameters: None
  use mpi_setup
  implicit none
  character ( len = 8 ) ampm
  integer*4 d
  integer*4 h
  integer*4 m
  integer*4 mm
  character(*) time_step
  character*5, parameter, dimension(12) :: month = (/ &
    'Jan. ', 'Feb. ', 'Mar. ', 'Apr. ', &
    'May. ', 'Jun. ', 'Jul. ', 'Aug. ', &
    'Sep. ', 'Oct. ', 'Nov. ', 'Dec. ' /)
  integer*4 n
  integer*4 s
  integer*4 values(8)
  integer*4 y
  real*8    t

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write (message, '(a,2x,i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim(time_step), d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
  call write_log(message,3,myid)

    call cpu_time(t)
  return
endsubroutine

endmodule

