#include "alias.inc"
subroutine get_eig(NN_TABLE, kp, nkp, PINPT, PPRAM, E, V, SV, neig, iband, nband, flag_vector, flag_sparse, flag_stat, flag_phase)
  use parameters, only: hopping, incar, energy, spmat, params
  use mpi_setup
  use time
  use memory
  use reorder_band
  use print_io
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (params ) :: PPRAM
  type (energy ) :: EE
  type (spmat  ) :: SHm, SHs
  integer*4  mpierr, iadd, ii
  integer*4  ik, my_ik, neig
  integer*4  iband,nband ! if sparse: iband = 1, nband = feast_nemax
  integer*4  nkp, is, ie,fe, im, fm
  real*8     percent, kp(3,nkp)
  real*8     E(nband*PINPT%nspin,nkp)
  complex*16 V(neig*PINPT%ispin,nband*PINPT%nspin,nkp)
  complex*16 SV(neig*PINPT%ispin,nband*PINPT%nspin,nkp) ! overlap integeral multiplied S * V
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  logical, intent(in) :: flag_vector, flag_stat
  real*8     t0, t1
  character*4 timer
  logical    flag_phase, flag_init, flag_sparse
  logical    flag_sparse_SHm, flag_sparse_SHs ! if .false. sparse Hamiltonian SHm(collinear magnetic) and SHs(SOC) will not
                                              ! be added up and constructed along with the k_loop due to the total number 
                                              ! of non-zero element is zero. This is determined in the get_ham_mag(soc)_sparse 
  integer*4  feast_ne(PINPT%nspin, nkp)
  integer*4  ncpu, id
  integer*4, allocatable :: ourjob(:), ourjob_disp(:), ourgroup(:)
  integer*4                 mygroup(0:nprocs-1),groupid

  if(flag_stat) then
    write(message,'(A)') '  ' ; write_msg
    write(message,'(A)') ' #--- START: BAND STRUCTURE EVALUATION -----------' ; write_msg
  endif
  timer = 'init'

  if(COMM_KOREA%flag_split) then
    ncpu = COMM_KOREA%nprocs
    id   = COMM_KOREA%myid
  else
    ncpu = nprocs
    id   = myid
  endif
#ifdef PSPARSE
  if(flag_sparse) then
#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    call get_npar(trim(PINPT%ifilenm(1)))
    allocate(ourgroup(npar))
    allocate(ourjob(npar))
    allocate(ourjob_disp(0:npar-1))
    call mpi_job_distribution_group(npar, nkp, ourgroup, mygroup, ourjob, ourjob_disp)
    call mpi_comm_anmeldung(COMM_JUELICH, npar, mygroup, COMM_JUELICH_RATHAUS)
    groupid = COMM_JUELICH%color
    id      = groupid
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#else
    write(message,'(A)')' DPSARSE and DMPI option should be activated together. Currently, only DPSPARSE is activated.' ; write_msg
    write(message,'(A)')' Please set DPSPARSE option in your makefile OPTIONS and recompile the code. Exit...' ; write_msg
    kill_job
#endif
  else
    allocate(ourjob(ncpu))
    allocate(ourjob_disp(0:ncpu-1))
  endif
#else
  allocate(ourjob(ncpu))
  allocate(ourjob_disp(0:ncpu-1))
#endif

#ifdef PSPARSE
  if(flag_sparse) then
    call report_job_distribution_group(flag_stat, ourjob, ourgroup)
  else
    call mpi_job_distribution_chain(nkp, ncpu, ourjob, ourjob_disp)
    if(.not. COMM_KOREA%flag_split) call report_job_distribution(flag_stat, ourjob)
  endif
#else
  call mpi_job_distribution_chain(nkp, ncpu, ourjob, ourjob_disp)
  if(.not. COMM_KOREA%flag_split) call report_job_distribution(flag_stat, ourjob)
#endif

  call initialize_all (EE, neig, nband, nkp, ourjob(id+1), &
                       PINPT, flag_vector, PPRAM%flag_use_overlap, flag_sparse, flag_stat, &
                       ii, iadd, t1, t0, flag_init)

  if_main call report_memory_total(PINPT%ispinor, PINPT%ispin, PINPT%nspin, neig, nband, nkp, &
                                   flag_stat, flag_sparse, ncpu)

  if(flag_stat) then 
    write(message,'(A)'             ) ' ' ; write_msg
    write(message,'(A)') '   STATUS: ' 
    call write_log(trim(message), 1, myid)
    call write_log(trim(message),22, myid)
  endif

 k_loop:do ik= sum(ourjob(1:id))+1, sum(ourjob(1:id+1))
    my_ik = ik - sum(ourjob(1:id))
    if(flag_sparse) then
#ifdef MKL_SPARSE
      !NOTE: SHm is k-independent -> keep unchankged from first call
      !      SHs is k-independent if flag_slater_koster
      !             k-dependent   if .not. slater_koster 
      call cal_eig_Hk_sparse(SHm, SHs, EE%E(:,ik), EE%V(:,:,my_ik), PINPT, PPRAM, NN_TABLE, kp(:,ik), &
                             neig, nband, flag_vector, flag_init, flag_phase, &
                             PINPT%feast_ne(1:PINPT%nspin,ik),&
                             flag_sparse_SHm, flag_sparse_SHs, timer)
#else

      write(message,'(A)')'    !WARN! The EWINDOW tag is only available if you have put -DMKL_SPARSE option'  ; write_msg
      write(message,'(A)')'           in your make file. Please find the details in the instruction. Exit program...'  ; write_msg
      kill_job
#endif
    elseif(.not.flag_sparse) then
      call cal_eig_Hk_dense ( Hm,  Hs, EE%E(:,ik), EE%V(:,:,my_ik), EE%SV(:,:,my_ik), PINPT, PPRAM, NN_TABLE, kp(:,ik), &
                             neig, iband, nband, flag_vector, flag_init, flag_phase,ik)
    endif

#ifdef PSPARSE
    if(flag_sparse) then
      if(flag_stat .and. myid .eq. 0) call print_eig_status(ik, ii, iadd, ourjob)
    elseif(.not. flag_sparse) then
      if(flag_stat .and. id .eq. 0) call print_eig_status(ik, ii, iadd, ourjob)
    endif
#else
    if(flag_stat .and. id .eq. 0) call print_eig_status(ik, ii, iadd, ourjob)
#endif

  enddo k_loop

#ifdef MPI
  if(flag_stat .and. myid .eq. 0) then
    write(message,'(A)')'  ' ; write_msg
    write(message,'(A)')'   Gathering all results to main node 0 ...' 
    call write_log(trim(message), 23, myid)
  endif

  if(.not. COMM_KOREA%flag_split) then
    if(.not. flag_sparse) then
      call MPI_ALLREDUCE(EE%E, E, size(E,kind=4), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)  ! share all results 
      if(flag_vector .and. nkp .ne. 1) then
        call MPI_GATHERV(EE%V,size(EE%V, kind=4), MPI_COMPLEX16, V, &
                         ourjob      *neig*PINPT%ispin*nband*PINPT%nspin, &
                         ourjob_disp *neig*PINPT%ispin*nband*PINPT%nspin, &
                         MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)  ! only main node keep wave vector information
        if(PPRAM%flag_use_overlap) then
          call MPI_GATHERV(EE%SV,size(EE%SV,kind=4), MPI_COMPLEX16, SV, &
                           ourjob      *     neig*PINPT%ispin*nband*PINPT%nspin , &
                           ourjob_disp *     neig*PINPT%ispin*nband*PINPT%nspin , &
                           MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)  ! only main node keep overlap matrix information
        endif
      endif                     
    elseif(flag_sparse) then
#ifdef PSPARSE      
      if(COMM_JUELICH_RATHAUS%mpi_comm .ne. MPI_COMM_NULL) then
        call MPI_ALLREDUCE(EE%E, E, size(E,kind=4), MPI_REAL8, MPI_SUM, COMM_JUELICH_RATHAUS%mpi_comm, mpierr)  ! share all results 
      endif
      if(flag_vector .and. nkp .ne. 1 .and. COMM_JUELICH_RATHAUS%mpi_comm .ne. MPI_COMM_NULL) then
        call MPI_GATHERV(EE%V,size(EE%V, kind=4), MPI_COMPLEX16, V, &
                         ourjob      *neig*PINPT%ispin*nband*PINPT%nspin, &
                         ourjob_disp *neig*PINPT%ispin*nband*PINPT%nspin, &
                         MPI_COMPLEX16, 0, COMM_JUELICH_RATHAUS%mpi_comm, mpierr)  ! only main node keep wave vector information
        if(PPRAM%flag_use_overlap) then
          call MPI_GATHERV(EE%SV,size(EE%SV,kind=4), MPI_COMPLEX16, SV, &
                           ourjob      *     neig*PINPT%ispin*nband*PINPT%nspin , &
                           ourjob_disp *     neig*PINPT%ispin*nband*PINPT%nspin , &
                           MPI_COMPLEX16, 0, COMM_JUELICH_RATHAUS%mpi_comm, mpierr)  ! only main node keep overlap matrix information
        endif
      endif
#else
      call MPI_ALLREDUCE(EE%E, E, size(E,kind=4), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)  ! share all results 
      if(flag_vector .and. nkp .ne. 1) then
        call MPI_GATHERV(EE%V,size(EE%V, kind=4), MPI_COMPLEX16, V, &
                         ourjob      *neig*PINPT%ispin*nband*PINPT%nspin, &
                         ourjob_disp *neig*PINPT%ispin*nband*PINPT%nspin, &
                         MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)  ! only main node keep wave vector information
        if(PPRAM%flag_use_overlap) then
          call MPI_GATHERV(EE%SV,size(EE%SV,kind=4), MPI_COMPLEX16, SV, &
                           ourjob      *     neig*PINPT%ispin*nband*PINPT%nspin , &
                           ourjob_disp *     neig*PINPT%ispin*nband*PINPT%nspin , &
                           MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)  ! only main node keep overlap matrix information
        endif
      endif
#endif
    endif
  elseif(COMM_KOREA%flag_split) then
    call MPI_ALLREDUCE(EE%E, E, size(E,kind=4), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)  ! share all results 
    if(flag_vector .and. nkp .ne. 1) then
      call MPI_GATHERV(EE%V,size(EE%V,kind=4), MPI_COMPLEX16, V, &
                       ourjob      *neig*PINPT%ispin*nband*PINPT%nspin, &
                       ourjob_disp *neig*PINPT%ispin*nband*PINPT%nspin, &
                       MPI_COMPLEX16, 0, COMM_KOREA%mpi_comm, mpierr)  ! only main node keep wave vector information
      if(PPRAM%flag_use_overlap) then
        call MPI_GATHERV(EE%SV,size(EE%SV,kind=4), MPI_COMPLEX16, SV, &
                         ourjob      *neig*PINPT%ispin*nband*PINPT%nspin, &
                         ourjob_disp *neig*PINPT%ispin*nband*PINPT%nspin, &
                         MPI_COMPLEX16, 0, COMM_KOREA%mpi_comm, mpierr)  ! only main node keep overlap matrix information
      endif
    endif
  endif
  if(flag_stat) then
    write(message,'(A)')'  done!' ; write_msg
    write(message,'(A  )')' ' ; write_msg
  endif
#ifdef MKL_SPARSE
  if(flag_sparse) then 
    if(.not.COMM_KOREA%flag_split) then
#ifdef PSPARSE
      call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, COMM_JUELICH%mpi_comm, mpierr)
      if(COMM_JUELICH_RATHAUS%mpi_comm .ne. MPI_COMM_NULL) then
        call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, COMM_JUELICH_RATHAUS%mpi_comm, mpierr)
      endif
#else
      call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
#endif
    elseif(COMM_KOREA%flag_split) then
      call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
    endif
    PINPT%feast_ne = feast_ne
    if(flag_stat) then
      write(message,'(A,I0)')'   MAX_NE_FOUND (NE_MAX): ',maxval(PINPT%feast_ne) ; write_msg
    endif
  endif
#endif
#else
  E = EE%E ; if(flag_vector) V = EE%V
#ifdef MKL_SPARSE
  if(flag_stat) then
    write(message,'(A,I0)')'   MAX_NE_FOUND (NE_MAX): ',maxval(PINPT%feast_ne) ; write_msg
  endif
#endif
#endif

  call finalize_all(EE, SHm, SHs, t1, t0, PINPT, flag_stat, flag_vector, flag_sparse)

  if(flag_stat) then
    write(message,'(A)')' #--- END: BAND STRUCTURE EVALUATION -----------' ; write_msg
  endif

return
endsubroutine

subroutine get_eig_sepk(NN_TABLE, kp, nkp, PINPT, PPRAM, E, neig, iband, nband, flag_vector, flag_sparse, flag_stat, flag_phase)
  use parameters, only: hopping, incar, energy, spmat, params
  use mpi_setup
  use time
  use memory
  use reorder_band
  use print_io
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (params ) :: PPRAM
  type (energy ) :: EE
  type (spmat  ) :: SHm, SHs
  integer*4  mpierr, iadd, ii
  integer*4  ik, neig, ijob
  integer*4  iband,nband ! if sparse: iband = 1, nband = feast_nemax
  integer*4  nkp, is, ie,fe, im, fm
  real*8     percent, kp(3,nkp)
  real*8     E(nband*PINPT%nspin, nkp)
  complex*16 V(neig*PINPT%ispin,nband*PINPT%nspin)
  complex*16 SV(neig*PINPT%ispin,nband*PINPT%nspin) ! overlap integeral multiplied S * V
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  logical, intent(in) :: flag_vector, flag_stat
  real*8     t0, t1
  character*4 timer
  character*80 fname_header, fname
  character*20,external ::   int2str
  logical    flag_file_exist
  logical    flag_phase, flag_init, flag_sparse
  logical    flag_sparse_SHm, flag_sparse_SHs ! if .false. sparse Hamiltonian SHm(collinear magnetic) and SHs(SOC) will not
                                              ! be added up and constructed along with the k_loop due to the total number 
                                              ! of non-zero element is zero. This is determined in the get_ham_mag(soc)_sparse 
  integer*4  feast_ne(PINPT%nspin, nkp)
  integer*4  ncpu, id
  integer*4, allocatable :: ourjob(:), ourjob_disp(:), ourgroup(:)
  integer*4, allocatable :: non_skip_list(:)
  integer*4                 mygroup(0:nprocs-1), groupid
  integer*4                 n_non_skip

  if(flag_stat) then
    write(message,'(A)') '  ' ; write_msg
    write(message,'(A)') ' #--- START: BAND STRUCTURE EVALUATION -----------' ; write_msg
  endif
  timer = 'init'

  EE%E = 0d0
  if(flag_vector) V = (0.d0,0.d0)
  if(flag_vector .and. PPRAM%flag_use_overlap) SV = (0.d0,0.d0)

  if(COMM_KOREA%flag_split) then
    ncpu = COMM_KOREA%nprocs
    id   = COMM_KOREA%myid
  else
    ncpu = nprocs
    id   = myid
  endif

  ! check non-skip list
  allocate(non_skip_list(nkp)) ; non_skip_list = 0 ; n_non_skip = 0
  do ik = 1, nkp
    ! check if the calculation is already done.. check restart
    fname_header = 'band_structure_TBA'//trim(PINPT%title(NN_TABLE%mysystem))//'.kp_'//trim(ADJUSTL(int2str(ik)))
    if(.not. PINPT%flag_write_unformatted) then
      call get_fname(fname_header, fname, 1, PINPT%flag_collinear, PINPT%flag_noncollinear)
    elseif(PINPT%flag_write_unformatted) then
      call get_fname_bin(fname_header, fname, 1, PINPT%flag_collinear, PINPT%flag_noncollinear)
    endif
    inquire(file=trim(fname),exist=flag_file_exist)
    if(flag_file_exist) then
        write(message,'(A,A,A,I0)')'    Band structure file:', trim(fname),' is already exist => skip to calculate KPT: ',ik ; write_msg
    else
        n_non_skip = n_non_skip + 1
        non_skip_list(n_non_skip) = ik
        cycle
    endif
    
    if(PINPT%nspin .eq. 2 .and. flag_file_exist) then
      if(.not. PINPT%flag_write_unformatted) then
        call get_fname(fname_header, fname, 2, PINPT%flag_collinear, PINPT%flag_noncollinear)
      elseif(PINPT%flag_write_unformatted) then
        call get_fname_bin(fname_header, fname, 2, PINPT%flag_collinear, PINPT%flag_noncollinear)
      endif
      inquire(file=trim(fname),exist=flag_file_exist)
      if(flag_file_exist) then
        write(message,'(A,A,A,I0)')'    Band structure file:', trim(fname),' is already exist => skip to calculate KPT: ',ik ; write_msg
      else
        n_non_skip = n_non_skip + 1
        non_skip_list(n_non_skip) = ik
        cycle
      endif
    endif
  enddo

#ifdef PSPARSE
  if(flag_sparse) then
#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    call get_npar(trim(PINPT%ifilenm(1)))
    allocate(ourgroup(npar))
    allocate(ourjob(npar))
    allocate(ourjob_disp(0:npar-1))
    call mpi_job_distribution_group(npar, n_non_skip, ourgroup, mygroup, ourjob, ourjob_disp)
    call mpi_comm_anmeldung(COMM_JUELICH, npar, mygroup, COMM_JUELICH_RATHAUS)
    groupid = COMM_JUELICH%color
    id      = groupid
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#else
    write(message,'(A)')' DPSARSE and DMPI option should be activated together. Currently, only DPSPARSE is activated.' ; write_msg
    write(message,'(A)')' Please set DPSPARSE option in your makefile OPTIONS and recompile the code. Exit...' ; write_msg
    kill_job
#endif
  else
    allocate(ourjob(ncpu))
    allocate(ourjob_disp(0:ncpu-1))
  endif
#else
  allocate(ourjob(ncpu))
  allocate(ourjob_disp(0:ncpu-1))
#endif

#ifdef PSPARSE
  if(flag_sparse) then
    call report_job_distribution_group(flag_stat, ourjob, ourgroup)
  else
    call mpi_job_distribution_chain(n_non_skip, ncpu, ourjob, ourjob_disp)
    if(.not. COMM_KOREA%flag_split) call report_job_distribution(flag_stat, ourjob)
  endif
#else
  call mpi_job_distribution_chain(n_non_skip, ncpu, ourjob, ourjob_disp)
  if(.not. COMM_KOREA%flag_split) call report_job_distribution(flag_stat, ourjob)
#endif

  call initialize_all (EE, neig, nband, nkp, ourjob(id+1), PINPT, flag_vector, PPRAM%flag_use_overlap, flag_sparse, flag_stat, &
                       ii, iadd, t1, t0, flag_init)

  if_main call report_memory_singlek(PINPT%ispinor, PINPT%ispin, PINPT%nspin, neig, nband, nkp, &
                                   flag_stat, flag_sparse, ncpu)

  if(flag_stat) then
    write(message,'(A)'             ) ' ' ; write_msg
    write(message,'(A)') '   STATUS: '
    call write_log(trim(message), 1, myid)
    call write_log(trim(message),22, myid)
  endif

 k_loop:do ijob= sum(ourjob(1:id))+1, sum(ourjob(1:id+1))
    ik = non_skip_list(ijob)

    if(flag_sparse) then
#ifdef MKL_SPARSE
      !NOTE: SHm is k-independent -> keep unchankged from first call
      !      SHs is k-independent if flag_slater_koster
      !             k-dependent   if .not. slater_koster 
      call cal_eig_Hk_sparse(SHm, SHs, EE%E(:,ik), V, PINPT, PPRAM, NN_TABLE, kp(:,ik), &
                             neig, nband, flag_vector, flag_init, flag_phase, &
                             PINPT%feast_ne(1:PINPT%nspin,ik), &
                             flag_sparse_SHm, flag_sparse_SHs, timer)

#else
      write(message,'(A)')'    !WARN! The EWINDOW tag is only available if you have put -DMKL_SPARSE option'  ; write_msg
      write(message,'(A)')'           in your make file. Please find the details in the instruction. Exit program...'  ; write_msg
      kill_job
#endif
    elseif(.not.flag_sparse) then
      call cal_eig_Hk_dense ( Hm,  Hs, EE%E(:,ik), V, SV, PINPT, PPRAM, NN_TABLE, kp(:,ik), &
                             neig, iband, nband, flag_vector, flag_init, flag_phase,ik)
    endif

    if(flag_stat .and. id .eq. 0) call print_eig_status(ijob, ii, iadd, ourjob)
    call print_energy_singlek(EE%E(:,ik), V, SV, neig, iband, nband, ik, kp(:,ik), &
                              PINPT, PPRAM%flag_use_overlap, NN_TABLE%mysystem)
  enddo k_loop

#ifdef MPI
! call MPI_BARRIER(mpi_comm_earth, mpierr)
  if(flag_stat .and. myid .eq. 0) then
    write(message,'(A)')'  ' ; write_msg
    write(message,'(A)')'   Gathering all results to main node 0 ...'
    call write_log(trim(message), 23, myid)
  endif

  if(.not. COMM_KOREA%flag_split) then
    if(.not. flag_sparse) then
      call MPI_ALLREDUCE(EE%E, E, size(E), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)  ! share all results 
    elseif(flag_sparse) then
#ifdef PSPARSE
      if(COMM_JUELICH_RATHAUS%mpi_comm .ne. MPI_COMM_NULL) then
        call MPI_ALLREDUCE(EE%E, E, size(E), MPI_REAL8, MPI_SUM, COMM_JUELICH_RATHAUS%mpi_comm, mpierr)
      endif
#else
      call MPI_ALLREDUCE(EE%E, E, size(E), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)  ! share all results 
#endif
    endif
  elseif(COMM_KOREA%flag_split) then
    call MPI_ALLREDUCE(EE%E, E, size(E), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)  ! share all results 
  endif

  if(flag_stat) then
    write(message,'(A)')'  done!' ; write_msg
    write(message,'(A  )')' ' ; write_msg
  endif

#ifdef MKL_SPARSE
  if(flag_sparse) then
    if(.not.COMM_KOREA%flag_split) then
#ifdef PSPARSE
      if(COMM_JUELICH_RATHAUS%mpi_comm .ne. MPI_COMM_NULL) then
        call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, COMM_JUELICH_RATHAUS%mpi_comm, mpierr)
      endif
#else
      call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
#endif
    elseif(COMM_KOREA%flag_split) then
      call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
    endif
    PINPT%feast_ne = feast_ne
    if(flag_stat) then
      write(message,'(A,I0)')'   MAX_NE_FOUND (NE_MAX): ',maxval(PINPT%feast_ne) ; write_msg
    endif
  endif
#endif

#else
  E = EE%E 
#ifdef MKL_SPARSE
  if(flag_stat) then
    write(message,'(A,I0)')'   MAX_NE_FOUND (NE_MAX): ',maxval(PINPT%feast_ne) ; write_msg
  endif
#endif

#endif

  call finalize_all(EE, SHm, SHs, t1, t0, PINPT, flag_stat, flag_vector, flag_sparse)

  if(flag_stat) then
    write(message,'(A)')' #--- END: BAND STRUCTURE EVALUATION -----------' ; write_msg
  endif

  ! NEED TO BE UPDATED HERE !!! HJ KIM 21.Oct 2020
  if(PINPT%flag_print_energy_singlek) then
    write(message,'(A)')'    !WARN! The calculation finished sucessfully. ' ; write_msg
    write(message,'(A)')'           Band structure information is printed into separate file "band_structure_TBA.kp_*'//trim(PINPT%title(NN_TABLE%mysystem))//'.dat" by request.' ; write_msg
    write(message,'(A)')'           However, due to some technical things, program will stop at this point. In the future release it will be updated.' ; write_msg
    write(message,'(A)')'          ' ; write_msg
    write(message,'(A)')'           Program stops...' ; write_msg
    kill_job
  endif

return
endsubroutine

#ifdef MKL_SPARSE
subroutine cal_eig_Hk_sparse(SHm, SHs, E, V, PINPT, PPRAM, NN_TABLE, kp, neig, &
                             nemax, flag_vector, flag_init, flag_phase, ne_found,&
                             flag_sparse_SHm, flag_sparse_SHs, timer)
  use parameters, only : incar, hopping, spmat, params
  use mpi_setup
  use sparse_tool
  use time
  use print_io
  implicit none
  type(hopping) :: NN_TABLE
  type(incar  ) :: PINPT
  type(params ) :: PPRAM
  type(spmat  ) :: SHk, SSk, SH0, SS0, SHm, SHs
  integer*4        mpierr
  integer*4        neig
  integer*4        nemax  !nemax <= neig * nspin. Choosing optimal value is critical for performance.
  integer*4        ne_found(PINPT%nspin)
  integer*4        ik
  integer*4        ie, fe, im, fm, is
  integer*4        feast_info
  real*8           kp(3)
  real*8           emin, emax
  real*8           E(nemax*PINPT%nspin)                   ! will store all the energy eigenvalue for each spin
  complex*16       V(neig*PINPT%ispin,nemax*PINPT%nspin)  ! will store all the spin block at once in the first dimension
  logical          flag_vector, flag_init, flag_phase, flag_sparse_SHm, flag_sparse_SHs
  real*8           t1, t0 
  character*4      timer

  emin = PINPT%feast_emin ; emax = PINPT%feast_emax

#ifdef MPI
#ifdef PSPARSE
  call MPI_BARRIER(COMM_JUELICH%mpi_comm, mpierr)
#endif
#endif

 sp:do is = 1, PINPT%nspin
      ! if(flag_init) SHm, SHs will be kept during ik-run. (SHs will be modified if flag_slater_koster=.false.)
      call time_check(t1,t0,timer)
      if(.not. PPRAM%flag_use_overlap) then
        call         get_hamk_sparse(SHk,      SH0,      SHm, SHs, is, kp, &
                                     PINPT, PPRAM, neig, NN_TABLE, flag_init, flag_phase, flag_sparse_SHm, flag_sparse_SHs)
      elseif(  PPRAM%flag_use_overlap) then
        call get_hamk_sparse_overlap(SHk, SSk, SH0, SS0, SHm, SHs, is, kp, &
                                     PINPT, PPRAM, neig, NN_TABLE, flag_init, flag_phase, flag_sparse_SHm, flag_sparse_SHs)
      endif
      call time_check(t1,t0) 

      if(timer .eq. 'init' .and. myid .eq. 0) then 
        write(message,'(A,F10.4,A)')'   TIME for SPARSE MATRIX CONSTRUCTION: ',t1, ' (sec)' ; write_msg
        timer = 'off'
      else
        timer = 'off'
      endif
      call get_matrix_index(ie, fe, im, fm, is, nemax, neig, PINPT%ispinor)
      if(.not. PPRAM%flag_use_overlap) then
#ifdef PSPARSE
        call cal_eig_hermitianx_psparse(SHk, emin, emax, nemax, ne_found(is), &
                                        PINPT%feast_neguess, E(ie:fe), V(im:fm,ie:fe), flag_vector,&
                                        PINPT%feast_fpm, feast_info)
#else
        call cal_eig_hermitianx_sparse( SHk, emin, emax, nemax, ne_found(is), &
                                        PINPT%feast_neguess, E(ie:fe), V(im:fm,ie:fe), flag_vector,&
                                        PINPT%feast_fpm, feast_info)
#endif
      elseif(  PPRAM%flag_use_overlap) then
#ifdef PSPARSE
        call cal_gen_eig_hermitianx_psparse(SHk, SSk, emin, emax, nemax, ne_found(is), &
                                            PINPT%feast_neguess, E(ie:fe),V(im:fm,ie:fe), flag_vector, &
                                            PINPT%feast_fpm, feast_info)
#else
        call cal_gen_eig_hermitianx_sparse( SHk, SSk, emin, emax, nemax, ne_found(is), &
                                            PINPT%feast_neguess, E(ie:fe),V(im:fm,ie:fe), flag_vector, &
                                            PINPT%feast_fpm, feast_info)
#endif
      endif

      call adjust_ne_guess(feast_info, is, ne_found(is), kp, neig, nemax, PINPT)
      if(allocated(SHk%H)) deallocate(SHk%H)
      if(allocated(SHk%I)) deallocate(SHk%I)
      if(allocated(SHk%J)) deallocate(SHk%J) ! SHk should be deallocated for each run
      if(allocated(SSk%H)) deallocate(SSk%H)
      if(allocated(SSk%I)) deallocate(SSk%I)
      if(allocated(SSk%J)) deallocate(SSk%J)
    enddo sp

  if(PINPT%flag_soc .and. .not. PPRAM%flag_slater_koster) then
    if(allocated(SHs%H)) deallocate(SHs%H)
    if(allocated(SHs%I)) deallocate(SHs%I)
    if(allocated(SHs%J)) deallocate(SHs%J)
  endif
  deallocate(SH0%H)
  deallocate(SH0%I)
  deallocate(SH0%J)
  if(allocated(SS0%H)) deallocate(SS0%H)
  if(allocated(SS0%I)) deallocate(SS0%I)
  if(allocated(SS0%J)) deallocate(SS0%J)

#ifdef MPI
#ifdef PSPARSE
  call MPI_BARRIER(COMM_JUELICH%mpi_comm, mpierr)
#endif
#endif

return
endsubroutine
#endif
subroutine cal_eig_Hk_dense(Hm, Hs, E, V, SV, PINPT, PPRAM, NN_TABLE, kp, neig, iband, nband, flag_vector, flag_init, flag_phase,ik)
  use parameters, only : incar, hopping, pauli_0, params
  use kronecker_prod, only: kproduct
  use mpi_setup
  use print_matrix
  use time
  use do_math
  use phase_factor
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (params ) :: PPRAM
  integer*4  neig, iband, nband
  integer*4  nkpoint, is, ie,fe, im, fm
  integer*4  ik
  real*8     kp(3)
  real*8     E(nband*PINPT%nspin)                      ! will store all the energy eigenvalue for each spin
  complex*16 V(neig*PINPT%ispin,nband*PINPT%nspin)     ! will store all the spin block at once in the first dimension
  complex*16 SV(neig*PINPT%ispin,nband*PINPT%nspin)     ! will store all the spin block at once in the first dimension
  complex*16 H0(neig,neig),S0(neig,neig)               ! slater-koster hopping (k-dependent) hamiltonian and overlap matrix
  complex*16 Hk(neig*PINPT%ispinor,neig*PINPT%ispinor) ! total hamiltonian (k-dependent)
  complex*16 Sk(neig*PINPT%ispinor,neig*PINPT%ispinor) ! total overlap matrix (k-dependent)
  complex*16 Sk_(neig*PINPT%ispinor,neig*PINPT%ispinor) ! temp. overlap matrix
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  logical    flag_vector, flag_init, flag_phase, flag_load_nntable
  real*8     t1, t0
  real*8     E_(4)
  character*20,external ::   int2str

  ! This routine calculates all the eigenvalues within [iband:iband+nband-1] using Hamiltonian Hk with dense matrix format.
   E_ = 0d0
   flag_load_nntable = PINPT%flag_load_nntable

 sp:do is = 1, PINPT%nspin
      ! if(flag_init) Hm, Hs will be kept during ik-run. (Hs will be modified if flag_slater_koster=.false.)
         call get_hamk_dense(Hk, H0, Hm, Hs, is, kp, PINPT, PPRAM, neig, NN_TABLE, flag_init, flag_phase) 
         if(PINPT%flag_print_hamk) then
           call print_matrix_c(Hk,neig*PINPT%ispinor, neig*PINPT%ispinor, &
                              'Hk_K'//trim(ADJUSTL(int2str(ik)))//'_SP'//trim(ADJUSTL(int2str(is))),1, 'F15.10')
         endif

         if(PPRAM%flag_use_overlap .and. is .eq. 1) then
           call set_ham0(S0, kp, PPRAM, neig, NN_TABLE, F_IJ, flag_phase, .true., flag_load_nntable) ! .true. : flag_set_overlap
           if(PINPT%flag_noncollinear) then 
             Sk = kproduct(pauli_0, S0, 2, 2, neig, neig)
           else
             Sk = S0
           endif
           if(PINPT%flag_print_hamk) then
             call print_matrix_c(Sk,neig*PINPT%ispinor, neig*PINPT%ispinor, &
                                 'Sk_K'//trim(ADJUSTL(int2str(ik)))//'_SP'//trim(ADJUSTL(int2str(is))),1, 'F15.10')
           endif
         endif

         call get_matrix_index(ie, fe, im, fm, is, nband, neig, PINPT%ispinor)
         
         if(.not. PPRAM%flag_use_overlap) then
           call cal_eig(Hk, neig, PINPT%ispinor, PINPT%ispin, iband, nband, E(ie:fe), V(im:fm,ie:fe), flag_vector)
         elseif(  PPRAM%flag_use_overlap) then
           if(flag_vector) then
             Sk_ = Sk ! save temp (since Sk is modified after calling ZHEGE or ZHEGEX)
           endif
           call cal_gen_eig(Hk, Sk, neig, PINPT%ispinor, PINPT%ispin, iband, nband, E(ie:fe), V(im:fm,ie:fe), flag_vector)

           if(flag_vector) then ! calculate V' = Sk|V> which will be used to calculate Mulliken charge rho = <V|Sk|V>
             SV = matprod(fm-im+1, fm-im+1, 'N', Sk_, fm-im+1, fe-ie+1, 'N', V(im:fm,ie:fe))
           endif
         endif

    enddo sp

return
endsubroutine
subroutine initialize_all(EE, neig, nband, nkp, my_nkp, PINPT, flag_vector, flag_overlap, flag_sparse, flag_stat, &
                          ii, iadd, t1, t0, flag_init)
  use parameters, only: incar, energy
  use mpi_setup
  use time
  implicit none
  type(incar) :: PINPT
  type(energy):: EE
  integer*4  neig, nband, nkp, my_nkp
  integer*4  iadd, ii
  logical    flag_vector, flag_stat, flag_init, flag_sparse, flag_overlap
  real*8     t1, t0
  integer*4  nL3

  flag_init = .true.
  call time_check(t1,t0,'init')

#ifdef MKL_SPARSE
  if(flag_sparse) then 
#ifdef PSPARSE
    nL3 = 1
#ifdef MPI
    call pfeastinit(PINPT%feast_fpm, COMM_JUELICH%mpi_comm, nL3)
#else
    call  feastinit(PINPT%feast_fpm)
#endif
#else
    call  feastinit(PINPT%feast_fpm)
#endif
    if(allocated(PINPT%feast_ne)) deallocate(PINPT%feast_ne)
    allocate(PINPT%feast_ne(PINPT%nspin, nkp))
    PINPT%feast_ne = 0 !initialize to zero
    PINPT%feast_neguess = PINPT%feast_nemax ! initialize to nemax. During calculations ne_guess will be adjusted using the ne_found in
                                            ! the previous step.

    PINPT%feast_fpm(1) = 0  ! Specifies whether Extended Eigensolver routines print runtime status (0:F, 1:T)
    PINPT%feast_fpm(2) = 4  ! The number of contour points N_e (see the description of FEAST algorithm
                            ! Ref: E. Polizzi, Phys. Rev. B 79, 115112 (2009) 
    PINPT%feast_fpm(3) = 11 ! Error trace double precisiion stopping criteria e ( e = 10^-feast_fpm(3) )
    PINPT%feast_fpm(4) = 20 ! Maximum number of Extended Eigensolver refinement loops allowed. 
                            ! If no convergence is reached within fpm(4) refinement loops, 
                            ! Extended Eigensolver routines return info=2.
    PINPT%feast_fpm(5) = 0  ! User initial subspace. If fpm(5)=0 then Extended Eigensolver routines generate
                            ! initial subspace, if fpm(5)=1 the user supplied initial subspace is used.
    PINPT%feast_fpm(6) = 1  ! Extended Eigensolver stopping test.
                            ! fpm(6)=0 : Extended Eigensolvers are stopped if the residual stopping test is satisfied.
                            ! fpm(6)=1 : Extended Eigensolvers are stopped if this trace stopping test is satisfied.
    PINPT%feast_fpm(14)= 0  ! If 1, return the computed eigenvectors subspace after one single contour integration.
  endif
#endif

  if(flag_stat) then 
    call initialize_eig_status(ii, iadd, nkp)
  endif

  allocate(EE%E(nband*PINPT%nspin, nkp))
  EE%E = 0d0 

  if(.not. PINPT%flag_print_energy_singlek) then
    if(allocated(EE%V)) deallocate(EE%V)
    if(allocated(EE%SV)) deallocate(EE%SV)
    allocate(EE%V(neig*PINPT%ispin, nband*PINPT%nspin, my_nkp))
    allocate(EE%SV(neig*PINPT%ispin, nband*PINPT%nspin, my_nkp))
    if(flag_vector) EE%V = (0.d0,0.d0)
    if(flag_vector .and. flag_overlap) EE%SV = (0.d0,0.d0)
  endif

return
endsubroutine
subroutine finalize_all(EE, SHm, SHs, t1, t0, PINPT, flag_stat, flag_vector, flag_sparse)
  use parameters, only : energy, spmat, incar
  use mpi_setup
  use time
  use print_io
  implicit none
  type(energy) :: EE
  type(spmat ) :: SHm, SHs
  type(incar ) :: PINPT
  logical    flag_vector, flag_stat, flag_sparse
  real*8     t1, t0

  call time_check(t1, t0)
  if(flag_stat) then
    write(message,*)' ' ; write_msg
    write(message,'(A,F20.6)')"   TIME for EIGENVALUE SOLVE (s)", t1 ; write_msg
  endif
  if(allocated(EE%E)) deallocate(EE%E)
  if(allocated(EE%V)) deallocate(EE%V)
  if(allocated(EE%SV)) deallocate(EE%SV)

  if(flag_sparse) then
    if(allocated(SHm%H)) deallocate(SHm%H)
    if(allocated(SHm%I)) deallocate(SHm%I)
    if(allocated(SHm%J)) deallocate(SHm%J)
    if(allocated(SHs%H)) deallocate(SHs%H)
    if(allocated(SHs%I)) deallocate(SHs%I)
    if(allocated(SHs%J)) deallocate(SHs%J)
  endif
return
endsubroutine
subroutine get_hamk_dense(Hk, H0, Hm, Hs, is, kpoint, PINPT, PPRAM, neig, NN_TABLE, flag_init, flag_phase)
  use parameters, only: hopping, incar, pauli_0, pauli_x, pauli_y, pauli_z, params
  use kronecker_prod, only: kproduct
  use mpi_setup
  use phase_factor
  use do_math
  use print_matrix
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (params ) :: PPRAM
  logical    flag_init, flag_phase, flag_load_nntable
  integer*4  neig, is
  real*8     kpoint(3)
  complex*16 H0(neig,neig)                             ! slater-koster hopping (k-dependent)
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  complex*16 Hk(neig*PINPT%ispinor,neig*PINPT%ispinor) ! total hamiltonian (k-dependent)
  logical    flag_set_overlap

  flag_load_nntable = PINPT%flag_load_nntable
  flag_set_overlap  = .false. ! in the first run set H

  if(is .eq. 1) call set_ham0(H0, kpoint, PPRAM, neig, NN_TABLE, F_IJ, flag_phase, flag_set_overlap, flag_load_nntable)

  if(flag_init) then
    if(PINPT%flag_collinear) then
      call set_ham_mag(Hm, NN_TABLE, PPRAM, neig, PINPT%ispinor, PINPT%flag_collinear, PINPT%flag_noncollinear)
    elseif(PINPT%flag_noncollinear) then
      call set_ham_mag(Hm, NN_TABLE, PPRAM, neig, PINPT%ispinor, PINPT%flag_collinear, PINPT%flag_noncollinear)
      if(PINPT%flag_soc .and. PPRAM%flag_slater_koster) &
        call set_ham_soc(Hs, 0d0, PPRAM, neig, NN_TABLE, F_IJ, flag_phase)
    endif
    flag_init = .false.
  endif

  if(PINPT%flag_collinear) then
    Hk = H0 + ((-1d0)**(is+1)) * Hm
  elseif(PINPT%flag_noncollinear) then
    if(PINPT%flag_soc) then
      !set up k-dependent SOC in the case of 'cc' orbitals
      if(.not. PPRAM%flag_slater_koster) &
        call set_ham_soc(Hs, kpoint, PPRAM, neig, NN_TABLE, F_IJ, flag_phase)
      Hk = kproduct(pauli_0, H0, 2, 2, neig, neig) + Hm + Hs
    else
      Hk = kproduct(pauli_0, H0, 2, 2, neig, neig) + Hm
    endif
  elseif(.not. PINPT%flag_collinear .and. .not. PINPT%flag_noncollinear) then
    Hk = H0
  endif

  return
endsubroutine
subroutine get_matrix_index(ie, fe, im, fm, is, nband, neig, ispinor)
  implicit none
  integer*4, intent(out)::  ie, fe, im, fm
  integer*4, intent(in) ::  is, nband, neig, ispinor

  ! initial and final index for row index of E(ie:fe) and column index of V(:,ie:fe)
  ie = 1 + (is-1)*nband
  fe = nband + (is-1)*nband

  ! initial and final index for column index of V(im:fm,:)
  im = 1 + (is-1)*neig 
  fm = neig*ispinor  + (is-1)*neig
 
return
endsubroutine

subroutine set_ham0(H, kpoint, PPRAM, neig, NN_TABLE, FIJ, flag_phase, flag_set_overlap, flag_load_nntable)
  use parameters, only : zi, hopping, params
  use phase_factor
  implicit none
  interface 
    function FIJ(k,R)
      complex*16   FIJ
      real*8, intent(in) :: k(3)
      real*8, intent(in) :: R(3)
    endfunction
  end interface
  type (hopping) :: NN_TABLE
  type (params ) :: PPRAM
  integer*4         neig , i, j, nn
  integer*4         ivel_axis
  real*8            kpoint(3), tol, tij_sk, tij_cc
  complex*16        H(neig ,neig)
  complex*16        Eij
  external          tij_sk, tij_cc
  logical           flag_phase  
  logical           flag_set_overlap
  logical           flag_load_nntable

  H=0.0d0
  tol=NN_TABLE%onsite_tolerance
  do nn=1,NN_TABLE%n_neighbor
    i=NN_TABLE%i_matrix(nn)
    j=NN_TABLE%j_matrix(nn)

    call get_hopping_integral(Eij, NN_TABLE, nn, PPRAM, tol, kpoint, FIJ, flag_phase, flag_set_overlap, flag_load_nntable)

    if(.not. flag_set_overlap) then
      if(i .eq. j .and. NN_TABLE%Dij(nn) <= tol) then
        if(PPRAM%slater_koster_type .gt. 10) then
          if(nint(PPRAM%param_const_nrl(4,1,NN_TABLE%local_U_param_index(i))) .ge. 1) then
            H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PPRAM%param_const_nrl(5,1,(NN_TABLE%local_U_param_index(i)))
          else
            H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PPRAM%param_nrl(1,(NN_TABLE%local_U_param_index(i)))
          endif
        else
          if(nint(PPRAM%param_const(4,NN_TABLE%local_U_param_index(i))) .ge. 1) then
            H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PPRAM%param_const(5,(NN_TABLE%local_U_param_index(i)))
          else
            H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PPRAM%param((NN_TABLE%local_U_param_index(i)))

          endif
        endif

      elseif(i .eq. j .and. NN_TABLE%Dij(nn) > tol) then
        H(i,j) = H(i,j) + Eij
      else
        H(i,j) = H(i,j) + Eij
        H(j,i) = H(j,i) + conjg(Eij)
      endif

    elseif( flag_set_overlap) then
      if(i .eq. j .and. NN_TABLE%Dij(nn) <= tol) then 
        H(i,j) = H(i,j) + 1d0
      elseif(i .eq. j .and. NN_TABLE%Dij(nn) > tol) then
        H(i,j) = H(i,j) + Eij
      else
        H(i,j) = H(i,j) + Eij
        H(j,i) = H(j,i) + conjg(Eij)
      endif
    endif

  enddo

return
endsubroutine

subroutine set_ham_mag(H, NN_TABLE, PPRAM, neig, ispinor, flag_collinear, flag_noncollinear)
    use parameters, only: zi, hopping, incar, params
    implicit none
    type (hopping) :: NN_TABLE
    type (incar  ) :: PINPT
    type (params ) :: PPRAM
    integer*4    neig, ispinor
    integer*4    i, ii
    complex*16   H(neig*ispinor,neig*ispinor)
    logical      flag_collinear, flag_noncollinear

    H=0d0
    if(flag_collinear) then
      do i = 1, neig
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(PPRAM%slater_koster_type .gt. 10) then
            if(nint(PPRAM%param_const_nrl(4,1,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i) = -0.5d0 * NN_TABLE%local_moment(1,i) * PPRAM%param_const_nrl(5,1,NN_TABLE%stoner_I_param_index(i))
            else
              H(i,i) = -0.5d0 * NN_TABLE%local_moment(1,i) * PPRAM%param_nrl(1,NN_TABLE%stoner_I_param_index(i))
            endif
          else
            if(nint(PPRAM%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i) = -0.5d0 * NN_TABLE%local_moment(1,i) * PPRAM%param_const(5,NN_TABLE%stoner_I_param_index(i))
            else
              H(i,i) = -0.5d0 * NN_TABLE%local_moment(1,i) * PPRAM%param(NN_TABLE%stoner_I_param_index(i))
            endif
          endif
        endif
      enddo

    elseif(flag_noncollinear) then
      do i = 1, neig  ! Hx
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(PPRAM%slater_koster_type .gt. 10) then
            if(nint(PPRAM%param_const_nrl(4,1,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i+neig) = H(i,i+neig) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param_const_nrl(5,1,NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i) = H(i+neig,i) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param_const_nrl(5,1,NN_TABLE%stoner_I_param_index(i))
            else
              H(i,i+neig) = H(i,i+neig) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param_nrl(1,NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i) = H(i+neig,i) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param_nrl(1,NN_TABLE%stoner_I_param_index(i))
            endif
          else
            if(nint(PPRAM%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i+neig) = H(i,i+neig) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param_const(5,NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i) = H(i+neig,i) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param_const(5,NN_TABLE%stoner_I_param_index(i))
            else
              H(i,i+neig) = H(i,i+neig) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param(NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i) = H(i+neig,i) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                                * PPRAM%param(NN_TABLE%stoner_I_param_index(i))
            endif
          endif
        endif
      enddo

      do i = 1, neig  ! Hy
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(PPRAM%slater_koster_type .gt. 10) then
            if(nint(PPRAM%param_const_nrl(4,1,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i+neig) = H(i,i+neig)  + 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param_const_nrl(5,1,NN_TABLE%stoner_I_param_index(i)) * zi
              H(i+neig,i) = H(i+neig,i)  - 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param_const_nrl(5,1,NN_TABLE%stoner_I_param_index(i)) * zi
            else
              H(i,i+neig) = H(i,i+neig)  + 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param_nrl(1,NN_TABLE%stoner_I_param_index(i)) * zi
              H(i,i+neig) = H(i,i+neig)  - 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param_nrl(1,NN_TABLE%stoner_I_param_index(i)) * zi
            endif
          else
            if(nint(PPRAM%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i+neig) = H(i,i+neig)  + 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param_const(5,NN_TABLE%stoner_I_param_index(i)) * zi
              H(i+neig,i) = H(i+neig,i)  - 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param_const(5,NN_TABLE%stoner_I_param_index(i)) * zi
            else
              H(i,i+neig) = H(i,i+neig)  + 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param(NN_TABLE%stoner_I_param_index(i)) * zi
              H(i,i+neig) = H(i,i+neig)  - 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                                 * PPRAM%param(NN_TABLE%stoner_I_param_index(i)) * zi
            endif
          endif
        endif
      enddo

      do i = 1, neig ! Hz
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(PPRAM%slater_koster_type .gt. 10) then
            if(nint(PPRAM%param_const_nrl(4,1,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i) = H(i,i) - 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                      * PPRAM%param_const_nrl(5,1,NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i+neig) = H(i+neig,i+neig) + 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                                          * PPRAM%param_const_nrl(5,1,NN_TABLE%stoner_I_param_index(i))
            else
              H(i,i) = H(i,i) - 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                      * PPRAM%param_nrl(1,NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i+neig) = H(i+neig,i+neig) + 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                                          * PPRAM%param_nrl(1,NN_TABLE%stoner_I_param_index(i))
            endif
          else
            if(nint(PPRAM%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
              H(i,i) = H(i,i) - 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                      * PPRAM%param_const(5,NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i+neig) = H(i+neig,i+neig) + 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                                          * PPRAM%param_const(5,NN_TABLE%stoner_I_param_index(i))
            else
              H(i,i) = H(i,i) - 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                      * PPRAM%param(NN_TABLE%stoner_I_param_index(i))
              H(i+neig,i+neig) = H(i+neig,i+neig) + 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                                          * PPRAM%param(NN_TABLE%stoner_I_param_index(i))
            endif
          endif
        endif
      enddo

    endif

return
endsubroutine

subroutine allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA, flag_distribute_vector, my_nkp)
   use parameters
   implicit none
   type (incar)   :: PINPT       ! parameters for input arguments
   type (energy)  :: ETBA        ! target energy to be fitted to
   type (poscar)  :: PGEOM       ! parameters for geometry info
   type (kpoints) :: PKPTS       ! parameters for kpoints
   logical           flag_distribute_vector
   integer*4         my_nkp, nkp_vector
   ! nband : number of eigenvalues to be stored for each spin
   !         = PGEOM%neig*PINPT%ispinor  (if .not. PINPT%flag_erange, default)
   !         = PINPT%feast_nemax (if EWINDOW tag is on and nemax is smaller than PGEOM%neig*PINPT%ispinor)
   !         = PGEOM%fina_erange - PGEOM%init_erange + 1 ( if PINPT%flag_erange)
   ! neig  : number of orbital basis
   ! nspin : 2 for collinear 1 for non-collinear
   ! ispin : 2 for collinear 2 for non-collinear

   if(flag_distribute_vector) then
     nkp_vector = my_nkp
   elseif(.not. flag_distribute_vector) then
     nkp_vector = PKPTS%nkpoint
   endif

   if(allocated(ETBA%E)) deallocate(ETBA%E)
   allocate(ETBA%E(PGEOM%nband*PINPT%nspin, PKPTS%nkpoint))

   ! need to find better way to allocate V and SV, since if we don't turn on LORBIT, 
   ! V and SV is not need to be saved. To save memory one need to make it simpler.
   ! But, now, just keep this way, to make my life easier. (H.-J. Kim, 01. Feb. 2021)
   if(.not. PINPT%flag_print_energy_singlek) then
     if(allocated(ETBA%V)) deallocate(ETBA%V)
     allocate(ETBA%V( PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin, nkp_vector))

     if(allocated(ETBA%SV)) deallocate(ETBA%SV)
     allocate(ETBA%SV(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin, nkp_vector))
   endif

   ! This can also be allocated if LROBIT = TRUE in future version (will check later on, H.-J. Kim, 01. Feb. 2021)
   ! In current version, this information is only used in get_orbital_projection routine which is only called by
   ! get_dE routine and is activated when PINPT%flag_fit_orbital = .true.
   if(PINPT%flag_fit_orbital) then
     if(allocated(ETBA%ORB)) deallocate(ETBA%ORB)
     allocate(ETBA%ORB(PINPT%lmmax,PGEOM%nband*PINPT%nspin, PKPTS%nkpoint))
   endif

   if(PINPT%flag_get_band_order) then
     if(allocated(ETBA%IDX)) deallocate(ETBA%IDX)
     if(allocated(ETBA%E_ORD)) deallocate(ETBA%E_ORD)
     if(allocated(ETBA%V_ORD)) deallocate(ETBA%V_ORD)
     if(allocated(ETBA%SV_ORD)) deallocate(ETBA%SV_ORD)
     allocate(ETBA%IDX(PGEOM%nband*PINPT%nspin, PKPTS%nkpoint))
     allocate(ETBA%E_ORD(PGEOM%nband*PINPT%nspin, PKPTS%nkpoint))
     allocate(ETBA%V_ORD(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin, nkp_vector))
     allocate(ETBA%SV_ORD(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin, nkp_vector))
   endif

   if(PINPT%flag_get_total_energy) then
     if(allocated(ETBA%F_OCC)) deallocate(ETBA%F_OCC)
     if(allocated(ETBA%E_BAND)) deallocate(ETBA%E_BAND)
     if(allocated(ETBA%E_TOT)) deallocate(ETBA%E_TOT)
     allocate(ETBA%F_OCC(PGEOM%nband*PINPT%nspin, PKPTS%nkpoint))
     allocate(ETBA%E_BAND(PINPT%nspin))
     allocate(ETBA%E_TOT (PINPT%nspin))
   endif
return
endsubroutine

subroutine initialize_eig_status(ii, iadd, nkpoint)
   use mpi_setup
   use print_io
   implicit none
   integer*4     iadd, ii, nkpoint
  
   if(nkpoint .le. 25) then
     iadd=10
     ii = 1
   else
     iadd=5
     ii = 1
   endif
return
endsubroutine
subroutine print_eig_status(ik, ii, iadd, ourjob)
   use mpi_setup
   use print_io
   implicit none
   integer*4 ik, ii, iadd
   integer*4 ourjob(nprocs)

   if( floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0) .ge. real(iadd*ii) ) then
     if( floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0) .ge. 100) then
       write(6,'(I0,A,$)') floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0),' Done!'
       write(6,'(A)')''
       write(message,'(I0,A)') floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0),' Done!'
       call write_log("            "//trim(message), 1, myid)
     else
       if(ik .gt. 1) then
         write(6,'(I0," ",$)') floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0)
         write(message,'(I3," ")') floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0)
         call write_log("            "//trim(message), 1 , myid)
       else
         write(6,'(4x, I0," ",$)')floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0) 
         write(message,'(I3," ")')floor((ik-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0) 
         call write_log("            "//trim(message), 1 , myid)
       endif
     endif
     ii = ii + 1
   endif

   return
endsubroutine

subroutine get_ham_Hk(Hk, NN_TABLE, PINPT, PPRAM, kpoint, is, neig, flag_phase)
   use parameters, only: hopping, incar, params
   use mpi_setup
   use phase_factor
   use do_math
   implicit none
   type (hopping) :: NN_TABLE
   type (incar  ) :: PINPT
   type (params ) :: PPRAM
   integer*4         neig
   integer*4         mpierr
   integer*4         is
   real*8            kpoint(3)
   complex*16 H0(neig,neig)                             ! slater-koster hopping (k-dependent)
   complex*16 Hk(neig*PINPT%ispinor,neig*PINPT%ispinor) ! total hamiltonian (k-dependent)
   complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
   complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
   logical    flag_phase, flag_init

   flag_init=.true.

   call get_hamk_dense(Hk, H0, Hm, Hs, is, kpoint, PINPT, PPRAM, neig, NN_TABLE, flag_init, flag_phase)

   return
endsubroutine

subroutine stop_get_eig(msize, nband)
   use mpi_setup
   use print_io
   implicit none
   integer*4    msize, nband

   write(message,'(A)')        '    !WARN! Check NERANGE! NBAND should be less equal to matrix size MSIZE (NBAND <= MSIZE = N_ORBIT*ISPINOR),'   ; write_msg
   write(message,'(A)')        '           where ISPINOR = 2/1 if LSORB = .TRUE./.FALSE. and N_ORBIT = total number of orbitals.'   ; write_msg
   write(message,'(A,I0,A,I0)')'           MSIZE = ',msize, ' , NBAND = FINA_E - INIT_E + 1 = ',nband  ; write_msg
   write(message,'(A)')        '           Exit program...'  ; write_msg
   stop

return
endsubroutine
subroutine cal_gen_eig(Hk, Sk, neig, ispinor, ispin, iband, nband, E, V, flag_vector)
   use do_math
   use mpi_setup
   implicit none
   integer*4, intent(in):: iband, nband, neig, ispinor, ispin
   integer*4              msize
   complex*16             Hk(neig*ispinor,neig*ispinor)
   complex*16             Sk(neig*ispinor,neig*ispinor)
   complex*16             V(neig*ispinor,nband)
   real*8                 E(nband)
   logical                flag_vector

   msize = neig * ispinor

   if(msize .eq. nband) then
     ! nband = neig*ispinor if no erange defined
     call cal_gen_eig_hermitian(Hk, Sk, msize, E, flag_vector)
     if(flag_vector) V = Hk
   elseif(msize .ne. nband .and. nband .lt. msize) then
     call cal_gen_eig_hermitianx(Hk, Sk, msize, iband, nband, E, V, flag_vector)
   else
     call stop_get_eig(msize,nband)
   endif

return
endsubroutine
subroutine cal_eig(Hk, neig, ispinor, ispin, iband, nband, E, V, flag_vector)
   use do_math
   use mpi_setup
   implicit none
   integer*4, intent(in):: iband, nband, neig, ispinor, ispin
   integer*4              msize
   complex*16             Hk(neig*ispinor,neig*ispinor)
   complex*16             V(neig*ispinor,nband)
   real*8                 E(nband)
   logical                flag_vector

   msize = neig * ispinor

   if(msize .eq. nband) then
     ! nband = neig*ispinor if no erange defined
     call cal_eig_hermitian(Hk, msize, E, flag_vector)
     if(flag_vector) V = Hk
   elseif(msize .ne. nband .and. nband .lt. msize) then
     call cal_eig_hermitianx(Hk, msize, iband, nband, E, V, flag_vector)
   else
     call stop_get_eig(msize,nband)
   endif

return
endsubroutine

#ifdef MKL_SPARSE
subroutine adjust_ne_guess(feast_info, is, ne_found, kp, neig, nemax, PINPT)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4    neig, feast_info
   real*8       kp(3)
   integer*4    is, ne_found, ik
   integer*4    nemax
 
   if(ne_found .eq. 0) then
     PINPT%feast_neguess = nint(PINPT%feast_nemax/2d0) + nint(PINPT%feast_nemax*0.4d0)
   elseif(ne_found .ge. 1) then
     PINPT%feast_neguess = nint(ne_found * 1.5 + 4)
     if(PINPT%feast_neguess .gt. nemax) PINPT%feast_neguess = nemax
   endif

return
endsubroutine
#endif
