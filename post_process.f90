#include "alias.inc"
subroutine post_process(PINPT, PPRAM, PPRAM_FIT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PINPT_BERRY, PINPT_DOS, PRPLT)
  use parameters, only: hopping, incar, energy, spmat, params, poscar, kpoints, weight, dos, berry, gainp, replot
  use set_default, only: init_params
  use mpi_setup
  use reorder_band
  use total_energy
#ifdef MKL_SPARSE
  use sparse_tool, only: feast_initialize
#endif
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
  integer*4                                  i, ik, my_ik
  integer*4                                  mpierr
  integer*8                                  buffer, buffer_limit, buffer2
  integer*4, allocatable                  :: feast_ne(:,:), ne_found(:)
  integer*4, allocatable                  :: ourjob(:), ourjob_disp(:)
  real*8,    allocatable                  :: E(:,:)
  if(PINPT%flag_tbfit_finish) then
    ! this print_mode applies to the read_input routine, and make not to print input settings as it reads 
    ! INCAR-TB again. This constraint make the output report redundunt
    print_mode = 99
  endif
  buffer_limit = int8(2)**int8(31) - int8(1)

  do i = 1, PINPT%nsystem
    write(message,'( A)')' '  ; write_msg
    write(message,'( A)')' #======================================================='  ; write_msg
    if(PINPT%flag_tbfit_finish) then
      write(message,'(2A)')'   START POST-PROCESSING PROCEDURE: ',trim(PGEOM(i)%gfilenm) ; write_msg
    else
      write(message,'(1A)')'   START POST-PROCESSING PROCEDURE: '                        ; write_msg
    endif
    write(message,'( A)')' #======================================================='  ; write_msg

    call read_input(PINPT,PPRAM(i),PKPTS(i), PGEOM(i), PWGHT(i), EDFT(i), NN_TABLE(i), &
                    PINPT_DOS(i), PINPT_BERRY(i), PKAIA, PRPLT(i), i)

    if(PINPT%flag_tbfit_finish) then
      call init_params(PPRAM(i), PINPT)
      PPRAM(i) = PPRAM_FIT ! just use fitted parameter if flag_tbfit_finish
    endif

    if(PRPLT(i)%flag_replot) cycle

    if(PINPT%flag_get_band .or. PINPT%flag_get_berry_curvature) then

      if(.not. PINPT%flag_print_energy_singlek) then
            if(.not. PINPT%flag_distribute_nkp) then
              call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i), .FALSE., 0)
              call get_eig(NN_TABLE(i), PKPTS(i)%kpoint, PKPTS(i)%nkpoint, PINPT, PPRAM(i), ETBA(i)%E, ETBA(i)%V, ETBA(i)%SV, PGEOM(i)%neig, &
                            PGEOM(i)%init_erange, PGEOM(i)%nband, PINPT%flag_get_orbital, PINPT%flag_sparse, .true., PINPT%flag_phase)


            ! Single K mode  ############
            elseif(PINPT%flag_distribute_nkp) then
#ifdef MKL_SPARSE
              if(PINPT%flag_sparse) then
                call feast_initialize(PINPT, PKPTS(i)%nkpoint)
#ifdef MPI
                if(allocated(feast_ne)) deallocate(feast_ne) ; allocate(feast_ne(PINPT%nspin, PKPTS(i)%nkpoint)); feast_ne = 0
#endif

                if(PINPT%nsystem .gt. 1) then
                    write(message, '(A)')'    !WARN! in somehow, PINPT%feast_ne variable should be stored in ETBA as ' ; write_msg
                    write(message, '(A)')'           PINPT is unique and disapear as do loop runs over "i"'  ; write_msg
                    write(message, '(A)')'           Please modify the source code! (message for myself.. as I am lazy now..' ; write_msg
                    write(message, '(A)')'           Left this tasks for the future works... :)'  ; write_msg
                    write(message, '(A)')'           HJ Kim. 19. Mar. 2021.'  ; write_msg
                    kill_job
                endif
              endif
#endif
              if(allocated(E))         deallocate(E)         ; allocate(E(PGEOM(i)%nband*PINPT%nspin,PKPTS(i)%nkpoint))
              if(allocated(ne_found))  deallocate(ne_found)  ; allocate(ne_found(PINPT%nspin)) ; ne_found = 0
              if(allocated(ourjob)) deallocate(ourjob) ; allocate(ourjob(nprocs))
              if(allocated(ourjob_disp)) deallocate(ourjob_disp) ; allocate(ourjob_disp(0:nprocs-1))
              call mpi_job_distribution_chain(PKPTS(i)%nkpoint, nprocs, ourjob, ourjob_disp)
              call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i), .TRUE., ourjob(myid+1))
              if(PINPT%flag_get_orbital) ETBA(i)%V = (0.d0,0.d0)
              if( (PINPT%flag_get_orbital .and. PPRAM(i)%flag_use_overlap) ) ETBA(i)%SV = (0.d0,0.d0)
#ifdef MKL_SPARSE
              PINPT%feast_neguess = PINPT%feast_nemax !initialize to nemax but ne_guess will be adjusted using the ne_found in the previous step.
#endif
              do ik= sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
                my_ik = ik - sum(ourjob(1:myid))
#ifdef MKL_SPARSE
                if(PINPT%flag_sparse) then
                  if(PINPT%flag_get_orbital)then
                    if(PINPT%feast_fpm(5) .eq. 1 .and. my_ik .gt. 1) then
                      ETBA(i)%V(:,:,my_ik) = ETBA(i)%V(:,:,my_ik-1)
                    endif
                  endif
                endif
#endif
                call get_eig_singlek(NN_TABLE(i), PKPTS(i)%kpoint, ik, PINPT, PPRAM(i), &
                                     ETBA(i)%E(:,ik), ETBA(i)%V(:,:,my_ik), ETBA(i)%SV(:,:,my_ik), PGEOM(i)%neig, &
                                     PGEOM(i)%init_erange, PGEOM(i)%nband, ne_found, &
                                     PINPT%flag_get_orbital, PINPT%flag_sparse, .true., PINPT%flag_phase)

#ifdef MKL_SPARSE
                PINPT%feast_ne(:,ik) = ne_found
#endif
              enddo
#ifdef MPI
              call MPI_ALLREDUCE(ETBA(i)%E, E, size(E), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr) ; ETBA(i)%E = E
#endif

#ifdef MKL_SPARSE
              if(PINPT%flag_sparse) then
#ifdef MPI
                call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
                PINPT%feast_ne = feast_ne
#endif        
                write(message,'(A,I0)')'   MAX_NE_FOUND (NE_MAX): ',maxval(PINPT%feast_ne) ; write_msg
              endif
#endif
            !##################################
            endif

      elseif(PINPT%flag_print_energy_singlek) then
        call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i), .FALSE., 0)
        call get_eig_sepk(NN_TABLE(i), PKPTS(i)%kpoint, PKPTS(i)%nkpoint, PINPT, PPRAM(i), ETBA(i)%E, PGEOM(i)%neig, &
                      PGEOM(i)%init_erange, PGEOM(i)%nband, PINPT%flag_get_orbital, PINPT%flag_sparse, .true., PINPT%flag_phase)
      endif

      call get_total_energy(ETBA(i), PINPT, PKPTS(i), PGEOM(i), .true., .true.)

      if(PINPT%flag_get_band_order .or. PINPT%flag_get_band_order_print_only) then
         call get_ordered_band(ETBA(i), PKPTS(i), PGEOM(i), PWGHT(i), PINPT, .false., PPRAM(i)%flag_use_overlap)
      endif

      call print_band_structure(PKPTS(i), ETBA(i), EDFT(i), PGEOM(i), PINPT, PPRAM(i), PWGHT(i))
      ! POST PROCESSING

      if(PINPT%flag_get_berry_curvature) call get_berry_curvature(NN_TABLE(i), PINPT, PINPT_BERRY(i), PGEOM(i), PKPTS(i), ETBA(i))
      if(PINPT%flag_get_circ_dichroism)  call get_circular_dichroism(NN_TABLE(i), PINPT, PGEOM(i), PKPTS(i), ETBA(i))
      if(PINPT%flag_plot_stm_image)      call plot_stm_image(PINPT,PGEOM(i),PKPTS(i), ETBA(i), PPRAM(i)%flag_use_overlap)
      if(PINPT%flag_plot_eigen_state)    call plot_eigen_state(PINPT,PGEOM(i),PKPTS(i), ETBA(i), PPRAM(i)%flag_use_overlap)
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
