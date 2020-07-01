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
  use reorder_band
  use print_io
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

  call parse_log(PINPT)
#ifdef MPI
  call mpi_initialize(PINPT%fnamelog)
#else
  call open_log(PINPT%fnamelog,myid)
#endif

  call version_stamp(t_start)
  call parse(PINPT, PKPTS)
  if_test call test()
          call read_input(PINPT,PINPT_DOS,PINPT_BERRY,PKPTS,PGEOM,PWGHT,EDFT,NN_TABLE,PKAIA,PRPLT)
  if( (.not. PRPLT%flag_replot_dos) .and. (.not. PRPLT%flag_replot_ldos) .and. (.not. PRPLT%flag_replot_sldos) .and. &
      (.not. PRPLT%flag_replot_proj_band) .and. (.not. PRPLT%flag_replot_band) .and. (.not. PRPLT%flag_replot_didv) ) then
    if(PINPT%flag_tbfit) call get_fit(PINPT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PKAIA)
    if(PINPT%flag_get_band .or. PINPT%flag_get_berry_curvature ) then
      call allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA)
      call get_eig(NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, PINPT, ETBA%E, ETBA%V, PGEOM%neig, &
                   PINPT%init_erange, PINPT%nband, PINPT%flag_get_orbital, PINPT%flag_sparse, .true., PINPT%flag_phase) !, PINPT%flag_get_band_order)
      if(PINPT%flag_get_band_order .or. PINPT%flag_get_band_order_print_only) then
         call get_ordered_band(ETBA, PKPTS%nkpoint, PGEOM%neig, PINPT%init_erange, PINPT%nband, PINPT, .false., PWGHT)
      endif
      if(PINPT%flag_print_energy_singlek) then
        write(message,'(A)')'    !WARN! Band structure information is printed into separate file "band_structure_TBA.kp_*.dat" by request.' ; write_msg
        write(message,'(A)')'           However, due to some technical things, program will stop at this point. In the future release it will be updated.' ; write_msg
        write(message,'(A)')'           Program stops...' ; write_msg
      else

        if(PINPT%flag_print_energy_diff) then
          if(PINPT%flag_get_band   .and. myid .eq. 0) call print_energy(PKPTS, ETBA%E, EDFT%E, ETBA%V, PGEOM%neig, PINPT, PWGHT, .TRUE.,'') 
        else
          if(PINPT%flag_get_band   .and. myid .eq. 0) then 
            call print_energy(PKPTS, ETBA%E, ETBA%E, ETBA%V, PGEOM%neig, PINPT, PWGHT, .TRUE., '') 
            if(PINPT%flag_get_band_order) then
              call print_energy(PKPTS, ETBA%E_ORD, ETBA%E_ORD, ETBA%V_ORD, PGEOM%neig, PINPT, PWGHT, .TRUE., '_ordered' ) 
            endif
          endif
        endif

        if(PINPT%flag_print_proj .and. myid .eq. 0) call print_energy_proj(PKPTS, ETBA%E, ETBA%V, PGEOM, PINPT) 
        if(PINPT%flag_get_berry_curvature) call get_berry_curvature(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS, ETBA)
        if(PINPT%flag_plot_stm_image)      call plot_stm_image(PINPT,PGEOM,PKPTS, ETBA)
        if(PINPT%flag_plot_eigen_state)    call plot_eigen_state(PINPT,PGEOM,PKPTS, ETBA)
        if(PINPT%flag_get_effective_ham)   call get_eig_downfold(PINPT,PKPTS,PGEOM,NN_TABLE) ! NOTE: THIS IS EXPERIMENTAL, BUT WORKS ANYWAY.
      endif
    endif
#ifdef MPI
    call MPI_Barrier(mpi_comm_earth,mpierr)
#endif

    ! In the following routines, they call "get_eig" to get eigenvalues & eigenvectors 
    if(PINPT%flag_get_dos)             call get_dos(NN_TABLE, PINPT, PINPT_DOS, PGEOM, PKPTS)
    if(PINPT%flag_get_zak_phase)       call get_zak_phase(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
    if(PINPT%flag_get_wcc)             call get_wcc(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
    if(PINPT%flag_get_z2)              call get_z2(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
    if(PINPT%flag_get_parity)          call get_parity(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
    if(PINPT%flag_get_symmetry)        call get_symmetry_eig(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)


  elseif( (PRPLT%flag_replot_dos .or. PRPLT%flag_replot_ldos .or. PRPLT%flag_replot_sldos &
      .or. PRPLT%flag_replot_proj_band .or.PRPLT%flag_replot_band .or. PRPLT%flag_replot_didv) ) then  ! program run in replot mode

    call replot_dos_band(PINPT, PGEOM, PKPTS, PRPLT)

  endif

  write(message,*)''  ; write_msg
  write(message,'(A)')' -------------------------------------------------------------------------' ; write_msg
  call timestamp ('| Program ends on',t_end)
  write(message,'(A)')' -------------------------------------------------------------------------' ; write_msg
  write(message,*)''; write_msg
  write(message,'(A,F13.3)')'Time elapsed (total, sec): ',t_end - t_start; write_msg
  write(message,*)''; write_msg

  if(allocated(ETBA%E))      deallocate(ETBA%E)
  if(PINPT%flag_get_orbital .and. allocated(ETBA%V) ) deallocate(ETBA%V)

#ifdef MPI
  call mpi_finish()
#endif
  call close_log(myid)
  stop
end program
