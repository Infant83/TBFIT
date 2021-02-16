#include "alias.inc"
subroutine post_process(PINPT, PPRAM, PPRAM_FIT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PINPT_BERRY, PINPT_DOS, PRPLT)
  use parameters, only: hopping, incar, energy, spmat, params, poscar, kpoints, weight, dos, berry, gainp, replot
  use set_default, only: init_params
  use mpi_setup
  use reorder_band
  use total_energy
  implicit none
  type(incar)                             :: PINPT
  type(params )                           :: PPRAM_FIT
  type(params ), dimension(PINPT%nsystem) :: PPRAM
  type(kpoints), dimension(PINPT%nsystem) :: PKPTS
  type(energy ), dimension(PINPT%nsystem) :: EDFT
  type(energy ), dimension(PINPT%nsystem) :: ETBA
  type(weight ), dimension(PINPT%nsystem) :: PWGHT
  type(poscar ), dimension(PINPT%nsystem) :: PGEOM
  type(hopping), dimension(PINPT%nsystem) :: NN_TABLE
  type(berry  ), dimension(PINPT%nsystem) :: PINPT_BERRY
  type(dos    ), dimension(PINPT%nsystem) :: PINPT_DOS
  type(replot ), dimension(PINPT%nsystem) :: PRPLT
  type(gainp  )                           :: PKAIA ! temp
  integer*4                                  i
  integer*4                                  mpierr
! integer*4                                  print_mode

 
  if(PINPT%flag_tbfit_finish) then
    ! this print_mode applies to the read_input routine, and make not to print input settings as it reads 
    ! INCAR-TB again. This constraint make the output report redundunt
    print_mode = 99
  endif

  do i = 1, PINPT%nsystem
    write(message,'( A)')' '  ; write_msg
    write(message,'( A)')' #======================================================='  ; write_msg
    if(PINPT%flag_tbfit_finish) then
      write(message,'(2A)')'   START POST-PROCESSING PROCEDURE: ',trim(PGEOM(i)%gfilenm) ; write_msg
    else
      write(message,'(1A)')'   START POST-PROCESSING PROCEDURE: '                        ; write_msg
    endif
    write(message,'( A)')' #======================================================='  ; write_msg
   !write(message,'( A)')' '  ; write_msg

!   call read_input(PINPT,PPRAM(i),PINPT_DOS(i),PINPT_BERRY(i),PKPTS(i),PGEOM(i),PWGHT(i), &
!                   EDFT(i), NN_TABLE(i),PKAIA,PRPLT(i), i)
    call read_input(PINPT,PPRAM(i),PKPTS(i), PGEOM(i), PWGHT(i), EDFT(i), NN_TABLE(i), &
                    PINPT_DOS(i), PINPT_BERRY(i), PKAIA, PRPLT(i), i)


    if(PINPT%flag_tbfit_finish) then
      call init_params(PPRAM(i), PINPT)
      PPRAM(i) = PPRAM_FIT ! just use fitted parameter if flag_tbfit_finish
    endif

    if(PRPLT(i)%flag_replot) cycle

    if(PINPT%flag_get_band .or. PINPT%flag_get_berry_curvature) then
      call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i))
      call get_eig(NN_TABLE(i), PKPTS(i)%kpoint, PKPTS(i)%nkpoint, PINPT, PPRAM, ETBA(i)%E, ETBA(i)%V, ETBA(i)%SV, PGEOM(i)%neig, &
                    PGEOM(i)%init_erange, PGEOM(i)%nband, PINPT%flag_get_orbital, PINPT%flag_sparse, .true., PINPT%flag_phase)

      call get_total_energy(ETBA(i), PINPT, PKPTS(i), PGEOM(i), .true., .true.)

      if(PINPT%flag_get_band_order .or. PINPT%flag_get_band_order_print_only) then
         call get_ordered_band(ETBA(i), PKPTS(i), PGEOM(i), PWGHT(i), PINPT, .false., PPRAM(i)%flag_use_overlap)
      endif
      call print_band_structure(PKPTS(i), ETBA(i), EDFT(i), PGEOM(i), PINPT, PPRAM(i), PWGHT(i))

      ! POST PROCESSING
      if(PINPT%flag_get_berry_curvature) call get_berry_curvature(NN_TABLE(i), PINPT, PINPT_BERRY(i), PGEOM(i), PKPTS(i), ETBA(i))
      if(PINPT%flag_get_circ_dichroism)  call get_circular_dichroism(NN_TABLE(i), PINPT, PGEOM(i), PKPTS(i), ETBA(i))
      if(PINPT%flag_plot_stm_image)      call plot_stm_image(PINPT,PGEOM(i),PKPTS(i), ETBA(i))
      if(PINPT%flag_plot_eigen_state)    call plot_eigen_state(PINPT,PGEOM(i),PKPTS(i), ETBA(i))
      if(PINPT%flag_get_effective_ham)   call get_eig_downfold(PINPT,PPRAM(i),PKPTS(i),PGEOM(i),NN_TABLE(i)) ! NOTE: THIS IS EXPERIMENTAL, BUT WORKS ANYWAY.
    endif

#ifdef MPI
    call MPI_Barrier(mpi_comm_earth,mpierr)
#endif

    ! In the following routines, they call "get_eig" to get eigenvalues & eigenvectors in their routine
    if(PINPT%flag_get_dos)             call get_dos(NN_TABLE(i), PINPT, PPRAM(i), PINPT_DOS(i), PGEOM(i), PKPTS(i))
    if(PINPT%flag_get_zak_phase)       call get_zak_phase(NN_TABLE(i), PINPT, PPRAM(i), PINPT_BERRY(i), PGEOM(i), PKPTS(i))
    if(PINPT%flag_get_wcc)             call get_wcc(NN_TABLE(i), PINPT, PPRAM(i), PINPT_BERRY(i), PGEOM(i), PKPTS(i))
    if(PINPT%flag_get_z2)              call get_z2(NN_TABLE(i), PINPT, PPRAM(i), PINPT_BERRY(i), PGEOM(i), PKPTS(i))
    if(PINPT%flag_get_parity)          call get_parity(NN_TABLE(i), PINPT, PPRAM(i), PINPT_BERRY(i), PGEOM(i), PKPTS(i))
    if(PINPT%flag_get_symmetry)        call get_symmetry_eig(NN_TABLE(i), PINPT, PPRAM(i), PINPT_BERRY(i), PGEOM(i), PKPTS(i))

  enddo

  return
endsubroutine
