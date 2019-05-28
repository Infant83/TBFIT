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
  type (replot)  :: PRPLT
#ifdef MPI
  call mpi_initialize()
! call test()
#endif
  if_main call version_stamp(t_start)
          call parse(PINPT, PKPTS)
  if_test call test()
          call read_input(PINPT,PINPT_DOS,PINPT_BERRY,PKPTS,PGEOM,PWGHT,EDFT,NN_TABLE,PKAIA,PRPLT)

  if( (.not. PRPLT%flag_replot_dos) .and. (.not. PRPLT%flag_replot_ldos) .and. (.not. PRPLT%flag_replot_sldos) .and. &
      (.not. PRPLT%flag_replot_proj_band) .and. (.not. PRPLT%flag_replot_band) ) then
    if(PINPT%flag_tbfit) call get_fit(PINPT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PKAIA)
    if(PINPT%flag_get_band .or. PINPT%flag_get_berry_curvature ) then
      call allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA)
      call get_eig(NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, PINPT, ETBA%E, ETBA%V, PGEOM%neig, &
                   PINPT%init_erange, PINPT%nband, PINPT%flag_get_orbital, PINPT%flag_sparse, .true., .true.)
      if(PINPT%flag_get_band   .and. myid .eq. 0) call print_energy(PKPTS, ETBA%E, ETBA%V, PGEOM, PINPT) 
      if(PINPT%flag_print_proj .and. myid .eq. 0) call print_energy_proj(PKPTS, ETBA%E, ETBA%V, PGEOM, PINPT) 
      if(PINPT%flag_get_berry_curvature) call get_berry_curvature(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS, ETBA)
      if(PINPT%flag_plot_stm_image)      call plot_stm_image(PINPT,PGEOM,PKPTS, ETBA)
      if(PINPT%flag_plot_eigen_state)    call plot_eigen_state(PINPT,PGEOM,PKPTS, ETBA)
      if(PINPT%flag_get_effective_ham)   call get_eig_downfold(PINPT,PKPTS,PGEOM,NN_TABLE) ! NOTE: THIS IS EXPERIMENTAL, BUT WORKS ANYWAY.
    endif
#ifdef MPI
    call MPI_Barrier(mpi_comm_earth,mpierr)
#endif

    ! In these routin, they call "get_eig" to get eigenvalues & eigenvectors 
    if(PINPT%flag_get_dos)             call get_dos(NN_TABLE, PINPT, PINPT_DOS, PGEOM, PKPTS)
    if(PINPT%flag_get_zak_phase)       call get_zak_phase(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
    if(PINPT%flag_get_wcc)             call get_wcc(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
    if(PINPT%flag_get_z2)              call get_z2(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
    if(PINPT%flag_get_parity)          call get_parity(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)


  elseif( (PRPLT%flag_replot_dos .or. PRPLT%flag_replot_ldos .or. PRPLT%flag_replot_sldos &
      .or. PRPLT%flag_replot_proj_band .or.PRPLT%flag_replot_band ) ) then  ! program run in replot mode

    call replot_dos_band(PINPT, PGEOM, PKPTS, PRPLT)

  endif

! call MPI_Barrier(mpi_comm_earth, mpierr)
  if_main write(6,*)''
  if_main write(6,'(A)')' -------------------------------------------------------------------------'
  if_main call timestamp ('| Program ends on',t_end)
  if_main write(6,'(A)')' -------------------------------------------------------------------------'
  if_main write(6,*)''
  if_main write(6,'(A,F13.3)')'Time elapsed (total, sec): ',t_end - t_start
  if_main write(6,*)''
  
  if(allocated(ETBA%E))      deallocate(ETBA%E)
  if(PINPT%flag_get_orbital .and. allocated(ETBA%V) ) deallocate(ETBA%V)

#ifdef MPI
  call mpi_finish()
#endif
  stop
end program
