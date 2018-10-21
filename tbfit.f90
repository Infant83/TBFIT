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
  implicit none
  external  get_eig
  real*8    t_start,t_end
  integer*4         mpierr
  type (incar)   :: PINPT       ! parameters for input arguments
  type (dos)     :: PINPT_DOS   ! parameters for density of states calculation
  type (berry)   :: PINPT_BERRY ! parameters for berry phase related quantity calculation
  type (poscar)  :: PGEOM       ! parameters for geometry info
  type (kpoints) :: PKPTS       ! parameters for kpoints
  type (energy)  :: EDFT        ! target energy to be fitted to
  type (energy)  :: ETBA        ! calculated energy with given TB-parameters
  type (weight)  :: PWGHT       ! weight factor for the fitting to target energy
  type (hopping) :: NN_TABLE    ! table for hopping index
  type (gainp)   :: PKAIA       ! input/control parameters for genetic algorithm

  !call test(PINPT, PKPTS, PGEOM)

#ifdef MPI
  call mpi_initialize()
#endif
  if_main call timestamp ('Program start on',t_start)
          call parse(PINPT) ; if_test call test()
          call read_input(PINPT,PINPT_DOS,PINPT_BERRY,PKPTS,PGEOM,PWGHT,EDFT,NN_TABLE,PKAIA)
  if(PINPT%flag_tbfit) call get_fit(PINPT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PKAIA)

  if(PINPT%flag_get_band .or. PINPT%flag_get_berry_curvature) then
    call allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA)
    call get_eig(NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, PINPT, ETBA%E, ETBA%V, PGEOM%neig, &
                 PINPT%init_erange, PINPT%nband, PINPT%flag_get_orbital, PINPT%flag_sparse, .true., .true.)
    if(PINPT%flag_get_band .and. myid .eq. 0) call print_energy(PKPTS, ETBA%E, ETBA%V, PGEOM, PINPT) 
    if(PINPT%flag_get_berry_curvature) call get_berry_curvature(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS, ETBA)
    if(PINPT%flag_plot_stm_image)      call plot_stm_image(PINPT,PGEOM,PKPTS, ETBA)
    if(PINPT%flag_plot_eigen_state)    call plot_eigen_state(PINPT,PGEOM,PKPTS, ETBA)
  endif

#ifdef MPI
  call MPI_Barrier(mpi_comm_earth,mpierr)
#endif

  if(PINPT%flag_get_dos)             call get_dos(NN_TABLE, PINPT, PINPT_DOS, PGEOM, PKPTS)

  if(PINPT%flag_get_zak_phase)       call get_zak_phase(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
  if(PINPT%flag_get_wcc)             call get_wcc(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
  if(PINPT%flag_get_z2)              call get_z2(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
  if(PINPT%flag_get_parity)          call get_parity(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)


! call MPI_Barrier(mpi_comm_earth, mpierr)
  if_main call timestamp ('Program ends on',t_end)
  if_main write(6,'(A,F13.3)')'  Time elapsed:',t_end - t_start
  
  deallocate(ETBA%E)
  if(PINPT%flag_get_orbital) deallocate(ETBA%V)

#ifdef MPI
  call mpi_finish()
#endif
  stop
end program
