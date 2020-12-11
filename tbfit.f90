#include "alias.inc"
program tbfit
!*****************************************************************************80
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Written:
!    13 Dec. 2017 ~ under development..
!  Author:
!    Hyun-Jung Kim (Korea Institute for Advanced Study), Infant@kias.re.kr
!
  use parameters
  use mpi_setup
  use time
  use version
  use print_io
  implicit none
  external  get_eig
  real*8    t_start,t_end
  integer*4         mpierr
  character*132     gnu_command
  type(incar  )                           :: PINPT       ! parameters for input arguments
  type(params )                           :: PPRAM_FIT   ! parameters fitted
  type(params ),allocatable, dimension(:) :: PPRAM       ! parameters for input TB parameters
  type(kpoints),allocatable, dimension(:) :: PKPTS
  type(energy ),allocatable, dimension(:) :: EDFT
  type(energy ),allocatable, dimension(:) :: ETBA
  type(weight ),allocatable, dimension(:) :: PWGHT
  type(poscar ),allocatable, dimension(:) :: PGEOM
  type(hopping),allocatable, dimension(:) :: NN_TABLE
  type(berry  ),allocatable, dimension(:) :: PINPT_BERRY
  type(dos    ),allocatable, dimension(:) :: PINPT_DOS
  type(replot ),allocatable, dimension(:) :: PRPLT

  call parse_very_init(PINPT)
  allocate(PPRAM(PINPT%nsystem))
  allocate(PKPTS(PINPT%nsystem))
  allocate(EDFT(PINPT%nsystem))
  allocate(ETBA(PINPT%nsystem))
  allocate(PWGHT(PINPT%nsystem))
  allocate(PGEOM(PINPT%nsystem))
  allocate(NN_TABLE(PINPT%nsystem))
  allocate(PINPT_BERRY(PINPT%nsystem))
  allocate(PINPT_DOS(PINPT%nsystem))
  allocate(PRPLT(PINPT%nsystem))

#ifdef MPI
  call mpi_initialize(PINPT%fnamelog)
#else
  call open_log(PINPT%fnamelog,myid)
#endif

  call version_stamp(t_start)
  call parse(PINPT)
  if_test call test()

  call get_fit(PINPT, PPRAM_FIT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PINPT_BERRY, PINPT_DOS)
    
  call post_process(PINPT, PPRAM, PPRAM_FIT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PINPT_BERRY, PINPT_DOS, PRPLT)
  
  call get_replot(PINPT, PGEOM, PKPTS, PRPLT)

  write(message,*)''  ; write_msg
  write(message,'(A)')' -------------------------------------------------------------------------' ; write_msg
  call timestamp ('| Program ends on',t_end)
  write(message,'(A)')' -------------------------------------------------------------------------' ; write_msg
  write(message,*)''; write_msg
  write(message,'(A,F13.3)')'Time elapsed (total, sec): ',t_end - t_start; write_msg
  write(message,*)''; write_msg

  if(PINPT%flag_plot) then
    write(gnu_command, '(A,A)')'gnuplot ', trim(PINPT%filenm_gnuplot)
    if_main call execute_command_line(gnu_command)
  endif

#ifdef MPI
  call mpi_finish()
#endif
  call close_log(myid)
  stop
end program
