#include "alias.inc"
subroutine get_eig(NN_TABLE, kp, nkp, PINPT, E, V, neig, iband, nband, flag_vector, flag_sparse, flag_stat, flag_phase)
  use parameters, only: hopping, incar, energy, spmat
  use mpi_setup
  use time
  use memory
! use print_matrix
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (energy ) :: EE
  type (spmat  ) :: SHm, SHs
  integer*4  mpierr, iadd, ii
  integer*4  ik, my_ik, neig
  integer*4  iband,nband ! if sparse: iband = 1, nband = feast_nemax
  integer*4  nkp, is, ie,fe, im, fm
  integer*4  ne_prev(PINPT%nspin)
  real*8     percent, kp(3,nkp)
  real*8     E(nband*PINPT%nspin,nkp)
  complex*16 V(neig*PINPT%ispin,nband*PINPT%nspin,nkp)
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  logical, intent(in) :: flag_vector, flag_stat
  character*100 stat
  real*8     t0, t1
  character*4 timer
  logical    flag_phase, flag_init, flag_sparse
  logical    flag_sparse_SHm, flag_sparse_SHs ! if .false. sparse Hamiltonian SHm(collinear magnetic) and SHs(SOC) will not
                                              ! be added up and constructed along with the k_loop due to the total number 
                                              ! of non-zero element is zero. This is determined in the get_ham_mag(soc)_sparse 
  integer*4  feast_ne(PINPT%nspin, nkp)
#ifdef MPI
  integer*4  ourjob(nprocs)
  integer*4  ourjob_disp(0:nprocs-1)
  call mpi_job_distribution_chain(nkp, ourjob, ourjob_disp)
  if(flag_stat .and. myid .eq. 0) write(6,'(A)') 'START: BAND STRUCTURE EVALUATION'
  call report_job_distribution(flag_stat, ourjob)
#else
  integer*4  ourjob(1) 
  integer*4  ourjob_disp(0)
  call mpi_job_distribution_chain(nkp, ourjob, ourjob_disp)
  if(flag_stat) write(6,'(A)') 'START: BAND STRUCTURE EVALUATION'
  call report_job_distribution(flag_stat, ourjob)
#endif

  timer = 'init'
  call initialize_all (EE, neig, nband, nkp, ourjob(myid+1), PINPT, flag_vector, flag_sparse, flag_stat, &
                       ii, iadd, stat, t1, t0, flag_init)
  if_main call report_memory_total(PINPT%ispinor, PINPT%ispin, PINPT%nspin, neig, nband, nkp, &
                                   flag_stat, flag_sparse, flag_use_mpi, nprocs)
 k_loop:do ik= sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
    my_ik = ik - sum(ourjob(1:myid))
    if(flag_sparse) then
      if(PINPT%feast_fpm(5) .eq. 1 .and. .not. flag_init) then 
!       EE%V(:,:,ik) = EE%V(:,:,ik-1)
        EE%V(:,:,my_ik) = EE%V(:,:,my_ik-1)
        ne_prev = PINPT%feast_ne(1:PINPT%nspin,ik-1)
      else
        ne_prev = 0
      endif
#ifdef MKL_SPARSE
      !NOTE: SHm is k-independent -> keep unchankged from first call
      !      SHs is k-independent if flag_slater_koster
      !             k-dependent   if .not. slater_koster 
      call cal_eig_Hk_sparse(SHm, SHs, EE%E(:,ik), EE%V(:,:,my_ik), PINPT, NN_TABLE, kp(:,ik), &
                             neig, nband, flag_vector, flag_init, flag_phase, &
                             PINPT%feast_ne(1:PINPT%nspin,ik),ik, &
                             flag_sparse_SHm, flag_sparse_SHs, ne_prev, timer)
      
#else

      if_main write(6,'(A)')'    !WARN! The EWINDOW tag is only available if you have put -DMKL_SPARSE option'
      if_main write(6,'(A)')'           in your make file. Please find the details in the instruction. Exit program...'
      kill_job
#endif
    elseif(.not.flag_sparse) then
      call cal_eig_Hk_dense ( Hm,  Hs, EE%E(:,ik), EE%V(:,:,my_ik), PINPT, NN_TABLE, kp(:,ik), &
                             neig, iband, nband, flag_vector, flag_init, flag_phase)
    endif

!   call print_eig_status(ik*ourjob(myid+1), ii, iadd, stat, nkp, flag_stat)
    call print_eig_status(ik-sum(ourjob(1:myid)), ii, iadd, stat, ourjob(myid+1), flag_stat)
  enddo k_loop

#ifdef MPI
  call MPI_ALLREDUCE(EE%E, E, size(E), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
  if(flag_vector) then
    call MPI_GATHERV(EE%V,size(EE%V), MPI_COMPLEX16, V, &
                     ourjob     *neig*PINPT%ispin*nband*PINPT%nspin, &
                     ourjob_disp*neig*PINPT%ispin*nband*PINPT%nspin, &
                     MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
  endif
#ifdef MKL_SPARSE
  if(flag_sparse) then 
    call MPI_ALLREDUCE(PINPT%feast_ne, feast_ne, size(feast_ne), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
    PINPT%feast_ne = feast_ne
  endif
#endif
#else
  E = EE%E ; if(flag_vector) V = EE%V
#endif

  call finalize_all(EE, SHm, SHs, t1, t0, PINPT, flag_stat, flag_vector, flag_sparse)
  if(flag_stat) then
    if_main write(6,'(A)')'END: BAND STRUCTURE EVALUATION'
  endif
return
endsubroutine

#ifdef MKL_SPARSE
subroutine cal_eig_Hk_sparse(SHm, SHs, E, V, PINPT, NN_TABLE, kp, neig, &
                             nemax, flag_vector, flag_init, flag_phase, ne_found,ik, &
                             flag_sparse_SHm, flag_sparse_SHs,ne_prev, timer)
  use parameters, only : incar, hopping, spmat
  use mpi_setup
  use sparse_tool
  use time
  implicit none
  type(hopping) :: NN_TABLE
  type(incar  ) :: PINPT
  type(spmat  ) :: SHk, SH0, SHm, SHs
  integer*4        neig
  integer*4        nemax  !nemax <= neig * nspin. Choosing optimal value is critical for performance.
  integer*4        ne_found(PINPT%nspin), ne_prev(PINPT%nspin)
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
   
 sp:do is = 1, PINPT%nspin
      ! if(flag_init) SHm, SHs will be kept during ik-run. (SHs will be modified if flag_slater_koster=.false.)
      call time_check(t1,t0,timer)
      call get_hamk_sparse(SHk, SH0, SHm, SHs, is, kp, PINPT, neig, NN_TABLE, flag_init, flag_phase, flag_sparse_SHm, flag_sparse_SHs)
      call time_check(t1,t0) 
      if(timer .eq. 'init' .and. myid .eq. 0) then 
        write(6,'(A,F10.4,A)')'   TIME for SPARSE MATRIX CONSTRUCTION: ',t1, ' (sec)'
        timer = 'off'
      else
        timer = 'off'
      endif
      call get_matrix_index(ie, fe, im, fm, is, nemax, neig, PINPT%ispinor)
      call cal_eig_sparse(SHk, neig, PINPT%ispinor, PINPT%ispin, nemax, PINPT%feast_neguess, E(ie:fe), V(im:fm,ie:fe), flag_vector, &
                          emin, emax, ne_found(is), PINPT%feast_fpm, feast_info, ne_prev(is))
      call adjust_ne_guess(feast_info, is, ne_found(is), kp, ik, neig, nemax, PINPT)
      if(allocated(SHk%H)) deallocate(SHk%H)
      if(allocated(SHk%I)) deallocate(SHk%I)
      if(allocated(SHk%J)) deallocate(SHk%J) ! SHk should be deallocated for each run
    enddo sp

  if(PINPT%flag_soc .and. .not. PINPT%flag_slater_koster) then
    if(allocated(SHs%H)) deallocate(SHs%H)
    if(allocated(SHs%I)) deallocate(SHs%I)
    if(allocated(SHs%J)) deallocate(SHs%J)
  endif
  deallocate(SH0%H)
  deallocate(SH0%I)
  deallocate(SH0%J)

return
endsubroutine
#endif
subroutine cal_eig_Hk_dense(Hm, Hs, E, V, PINPT, NN_TABLE, kp, neig, iband, nband, flag_vector, flag_init, flag_phase)
  use parameters, only : incar, hopping
  use mpi_setup
  use print_matrix
  use time
  use do_math
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  integer*4  neig, iband, nband
  integer*4  nkpoint, is, ie,fe, im, fm
  real*8     kp(3)
  real*8     E(nband*PINPT%nspin)                      ! will store all the energy eigenvalue for each spin
  complex*16 V(neig*PINPT%ispin,nband*PINPT%nspin)     ! will store all the spin block at once in the first dimension
  complex*16 H0(neig,neig)                             ! slater-koster hopping (k-dependent)
  complex*16 Hk(neig*PINPT%ispinor,neig*PINPT%ispinor) ! total hamiltonian (k-dependent)
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  logical    flag_vector, flag_init, flag_phase
  real*8     t1, t0
  real*8     E_(4)
  ! This routine calculates all the eigenvalues within [iband:iband+nband-1] using Hamiltonian Hk with dense matrix format.
   E_ = 0d0
 sp:do is = 1, PINPT%nspin
      ! if(flag_init) Hm, Hs will be kept during ik-run. (Hs will be modified if flag_slater_koster=.false.)
      call get_hamk_dense(Hk, H0, Hm, Hs, is, kp, PINPT, neig, NN_TABLE, flag_init, flag_phase) 
      call get_matrix_index(ie, fe, im, fm, is, nband, neig, PINPT%ispinor)
      call cal_eig(Hk, neig, PINPT%ispinor, PINPT%ispin, iband, nband, E(ie:fe), V(im:fm,ie:fe), flag_vector)
    enddo sp

return
endsubroutine
subroutine initialize_all(EE, neig, nband, nkp, my_nkp, PINPT, flag_vector, flag_sparse, flag_stat, &
                          ii, iadd, stat, t1, t0, flag_init)
  use parameters, only: incar, energy
  use mpi_setup
  use time
  implicit none
  type(incar) :: PINPT
  type(energy):: EE
  integer*4  neig, nband, nkp, my_nkp
  integer*4  iadd, ii
  logical    flag_vector, flag_stat, flag_init, flag_sparse
  character*100 stat
  real*8     t1, t0

  flag_init = .true.
  call time_check(t1,t0,'init')

#ifdef MKL_SPARSE
  if(flag_sparse) then 
    call feastinit(PINPT%feast_fpm)
    if(allocated(PINPT%feast_ne)) deallocate(PINPT%feast_ne)
    allocate(PINPT%feast_ne(PINPT%nspin, nkp))
    PINPT%feast_ne = 0 !initialize to zero
    PINPT%feast_neguess = PINPT%feast_nemax ! initialize to nemax. During calculations ne_guess will be adjusted using the ne_found in
                                            ! the previous step.

    PINPT%feast_fpm(1) = 0  ! Specifies whether Extended Eigensolver routines print runtime status (0:F, 1:T)
    PINPT%feast_fpm(2) = 4  ! The number of contour points N_e (see the description of FEAST algorithm
                            ! Ref: E. Polizzi, Phys. Rev. B 79, 115112 (2009) 
    PINPT%feast_fpm(3) = 11 ! Error trace double precisiion stopping criteria e ( e = 10^-feast_fpm(3) )
    PINPT%feast_fpm(4) = 3  ! Maximum number of Extended Eigensolver refinement loops allowed. 
                            ! If no convergence is reached within fpm(4) refinement loops, 
                            ! Extended Eigensolver routines return info=2.
    PINPT%feast_fpm(5) = 0  ! User initial subspace. If fpm(5)=0 then Extended Eigensolver routines generate
                            ! initial subspace, if fpm(5)=1 the user supplied initial subspace is used.
    PINPT%feast_fpm(6) = 1  ! Extended Eigensolver stopping test.
                            ! fpm(6)=0 : Extended Eigensolvers are stopped if the residual stopping test is satisfied.
                            ! fpm(6)=1 : Extended Eigensolvers are stopped if this trace stopping test is satisfied.
    PINPT%feast_fpm(7) = 5  ! Error trace single precision stopping criteria (10-fpm(7)).
    PINPT%feast_fpm(14)= 0  ! If 1, return the computed eigenvectors subspace after one single contour integration.
    PINPT%feast_fpm(27)= 0  ! Specifies whether Extended Eigensolver routines check input matrices (applies to CSR format only).
    PINPT%feast_fpm(28)= 0  ! Check if matrix B is positive definite. Set fpm(28) = 1 to check if B is positive definite.
  endif
#endif

  if(flag_stat) then 
    call initialize_eig_status(ii, iadd, stat, nkp)
  endif

  allocate(EE%E(nband*PINPT%nspin, nkp))
! allocate(EE%V(neig*PINPT%ispin, nband*PINPT%nspin, nkp))
  allocate(EE%V(neig*PINPT%ispin, nband*PINPT%nspin, my_nkp))
  EE%E = 0d0 ; if(flag_vector) EE%V = (0.d0,0.d0)

return
endsubroutine
subroutine finalize_all(EE, SHm, SHs, t1, t0, PINPT, flag_stat, flag_vector, flag_sparse)
  use parameters, only : energy, spmat, incar
  use mpi_setup
  use time
  implicit none
  type(energy) :: EE
  type(spmat ) :: SHm, SHs
  type(incar ) :: PINPT
  logical    flag_vector, flag_stat, flag_sparse
  real*8     t1, t0

  call time_check(t1, t0)
  if(flag_stat .and. myid .eq. 0) write(6,'(A,F12.6)')"TIME for EIGENVALUE SOLVE (s)", t1

  if(allocated(EE%E)) deallocate(EE%E)
  if(allocated(EE%V)) deallocate(EE%V)

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
subroutine get_hamk_dense(Hk, H0, Hm, Hs, is, kpoint, PINPT, neig, NN_TABLE, flag_init, flag_phase)
  use parameters, only: hopping, incar, pauli_0, pauli_x, pauli_y, pauli_z
  use kronecker_prod, only: kproduct
  use mpi_setup
  use phase_factor
  use do_math
  use print_matrix
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  logical    flag_init, flag_phase
  integer*4  neig, is
  real*8     kpoint(3)
  complex*16 H0(neig,neig)                             ! slater-koster hopping (k-dependent)
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  complex*16 Hk(neig*PINPT%ispinor,neig*PINPT%ispinor) ! total hamiltonian (k-dependent)

  if(is .eq. 1) call set_ham0(H0, kpoint, PINPT, neig, NN_TABLE, F_IJ, flag_phase)
  if(flag_init) then
    if(PINPT%flag_collinear) then
      call set_ham_mag(Hm, NN_TABLE, PINPT, neig)
    elseif(PINPT%flag_noncollinear) then
      call set_ham_mag(Hm, NN_TABLE, PINPT, neig)
      if(PINPT%flag_soc .and. PINPT%flag_slater_koster) &
        call set_ham_soc(Hs, 0d0, PINPT, neig, NN_TABLE, F_IJ, flag_phase)
    endif
    flag_init = .false.
  endif

  if(PINPT%flag_collinear) then
    Hk = H0 + ((-1d0)**(is+1)) * Hm
  elseif(PINPT%flag_noncollinear) then
    if(PINPT%flag_soc) then
      !set up k-dependent SOC in the case of 'cc' orbitals
      if(.not. PINPT%flag_slater_koster) &
        call set_ham_soc(Hs, kpoint, PINPT, neig, NN_TABLE, F_IJ, flag_phase)
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
subroutine set_ham0_file(H, kpoint, PINPT, neig, NN_TABLE, FIJ, flag_phase)
  use parameters, only : zi, hopping, incar
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
  type (incar  ) :: PINPT
  integer*4         neig , i, j, nn
  integer*4         ivel_axis
  real*8            kpoint(3), tol, tij_sk, tij_cc
  complex*16        H(neig ,neig)
  complex*16        Eij
  external          tij_sk, tij_cc
  logical           flag_phase

  H=0.0d0
  tol=NN_TABLE%onsite_tolerance
  do nn=1,NN_TABLE%n_neighbor
    i=NN_TABLE%i_matrix(nn)
    j=NN_TABLE%j_matrix(nn)

    if(PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = NN_TABLE%tij_file(nn)         * FIJ( kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = NN_TABLE%tij_file(nn)         * FIJ( kpoint, NN_TABLE%R  (1:3,nn))
      endif
    elseif(.not.PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = NN_TABLE%tij_file(nn)         * FIJ( kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = NN_TABLE%tij_file(nn)         * FIJ( kpoint, NN_TABLE%R  (1:3,nn))
      endif
    endif

    if(i .eq. j .and. NN_TABLE%Dij(nn) <= tol) then
      if(nint(PINPT%param_const(4,NN_TABLE%local_U_param_index(i))) .ge. 1) then
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param_const(5,(NN_TABLE%local_U_param_index(i)))
      else
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param((NN_TABLE%local_U_param_index(i)))
      endif
    elseif(i .eq. j .and. NN_TABLE%Dij(nn) >  tol) then
      H(i,j) = H(i,j) + Eij
    else
      H(i,j) = H(i,j) + Eij
      H(j,i) = H(j,i) + conjg(Eij)
    endif
  enddo

return
endsubroutine
subroutine set_ham0(H, kpoint, PINPT, neig, NN_TABLE, FIJ, flag_phase)
  use parameters, only : zi, hopping, incar
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
  type (incar  ) :: PINPT
  integer*4         neig , i, j, nn
  integer*4         ivel_axis
  real*8            kpoint(3), tol, tij_sk, tij_cc
  complex*16        H(neig ,neig)
  complex*16        Eij
  external          tij_sk, tij_cc
  logical           flag_phase  

  H=0.0d0
  tol=NN_TABLE%onsite_tolerance
  do nn=1,NN_TABLE%n_neighbor
    i=NN_TABLE%i_matrix(nn)
    j=NN_TABLE%j_matrix(nn)

    if(PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ( kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ( kpoint, NN_TABLE%R  (1:3,nn))
      endif
    elseif(.not.PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ( kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ( kpoint, NN_TABLE%R  (1:3,nn))
      endif
    endif

!   if(i .eq. j ) then !.and. NN_TABLE%Dij(nn) <= tol) then
    if(i .eq. j .and. NN_TABLE%Dij(nn) <= tol) then
      if(nint(PINPT%param_const(4,NN_TABLE%local_U_param_index(i))) .ge. 1) then
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param_const(5,(NN_TABLE%local_U_param_index(i)))
      else
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param((NN_TABLE%local_U_param_index(i)))
      endif
    elseif(i .eq. j .and. NN_TABLE%Dij(nn) > tol) then
      H(i,j) = H(i,j) + Eij
    else
      H(i,j) = H(i,j) + Eij
      H(j,i) = H(j,i) + conjg(Eij)
    endif
  enddo

return
endsubroutine

subroutine set_ham0_(H, kpoint, PINPT, neig , NN_TABLE, FIJ, flag_phase)
  use parameters, only : zi, hopping, incar
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
  type (incar  ) :: PINPT
  integer*4         neig , i, j, nn
  integer*4         ivel_axis
  real*8            kpoint(3), tol, tij_sk, tij_cc
  complex*16        H(neig ,neig)
  complex*16        Eij
  external          tij_sk, tij_cc
  logical           flag_phase

  H=0.0d0
  tol=NN_TABLE%onsite_tolerance
  do nn=1,NN_TABLE%n_neighbor
    i=NN_TABLE%i_matrix(nn)
    j=NN_TABLE%j_matrix(nn)

    if(PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    elseif(.not.PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    endif

    if(i .eq. j .and. NN_TABLE%Dij(nn) <= tol) then
      if(nint(PINPT%param_const(4,NN_TABLE%local_U_param_index(i))) .ge. 1) then
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param_const(5,(NN_TABLE%local_U_param_index(i)))
      else
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param((NN_TABLE%local_U_param_index(i)))
      endif
    elseif(i .eq. j .and. NN_TABLE%Dij(nn) > tol) then
      H(i,j) = H(i,j) + Eij
    else
      H(i,j) = H(i,j) + Eij
      H(j,i) = H(j,i) + conjg(Eij)
    endif
  enddo

return
endsubroutine



subroutine set_ham0_vel(H, kpoint, PINPT, neig , NN_TABLE, FIJ, flag_phase)
  use parameters, only : zi, hopping, incar
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
  type (incar  ) :: PINPT
  integer*4         neig , i, j, nn
  integer*4         ivel_axis
  real*8            kpoint(3), tol, tij_sk, tij_cc
  complex*16        H(neig ,neig)
  complex*16        Eij
  external          tij_sk, tij_cc
  logical           flag_phase

  H=0.0d0
  tol=NN_TABLE%onsite_tolerance
  do nn=1,NN_TABLE%n_neighbor
    i=NN_TABLE%i_matrix(nn)
    j=NN_TABLE%j_matrix(nn)

    if(PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    elseif(.not.PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    endif

    H(i,j) = H(i,j) + Eij
    if(i .ne. j) H(j,i) = H(j,i) + conjg(Eij)
  enddo

return
endsubroutine

subroutine set_ham_mag(H, NN_TABLE, PINPT, neig)
    use parameters, only: zi, hopping, incar
    implicit none
    type (hopping) :: NN_TABLE
    type (incar  ) :: PINPT
    integer*4    neig
    integer*4    i, ii
    complex*16   H(neig*PINPT%ispinor,neig*PINPT%ispinor)

    H=0d0
    if(PINPT%flag_collinear) then
      do i = 1, neig
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(nint(PINPT%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
            H(i,i) = -0.5d0 * NN_TABLE%local_moment(1,i) * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i))
          else
            H(i,i) = -0.5d0 * NN_TABLE%local_moment(1,i) * PINPT%param(NN_TABLE%stoner_I_param_index(i))
          endif
        endif
      enddo

    elseif(PINPT%flag_noncollinear) then
      do i = 1, neig  ! Hx
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(nint(PINPT%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
            H(i,i+neig) = H(i,i+neig) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                              * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i))
            H(i+neig,i) = H(i+neig,i) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                              * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i))
          else
            H(i,i+neig) = H(i,i+neig) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                              * PINPT%param(NN_TABLE%stoner_I_param_index(i))
            H(i+neig,i) = H(i+neig,i) - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                              * PINPT%param(NN_TABLE%stoner_I_param_index(i))
          endif
        endif
      enddo

      do i = 1, neig  ! Hy
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(nint(PINPT%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
            H(i,i+neig) = H(i,i+neig)  + 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                               * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i)) * zi
            H(i+neig,i) = H(i+neig,i)  - 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                               * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i)) * zi
          else
            H(i,i+neig) = H(i,i+neig)  + 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                               * PINPT%param(NN_TABLE%stoner_I_param_index(i)) * zi
            H(i,i+neig) = H(i,i+neig)  - 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                               * PINPT%param(NN_TABLE%stoner_I_param_index(i)) * zi
          endif
        endif
      enddo

      do i = 1, neig ! Hz
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(nint(PINPT%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
            H(i,i) = H(i,i) - 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                    * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i))
            H(i+neig,i+neig) = H(i+neig,i+neig) + 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                                        * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i))
          else
            H(i,i) = H(i,i) - 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                    * PINPT%param(NN_TABLE%stoner_I_param_index(i))
            H(i+neig,i+neig) = H(i+neig,i+neig) + 0.5d0 * NN_TABLE%local_moment_rot(3,i) &
                                                        * PINPT%param(NN_TABLE%stoner_I_param_index(i))
          endif
        endif
      enddo

    endif

return
endsubroutine

subroutine report_job_distribution(flag_stat, ourjob)
   use mpi_setup
   implicit none
   integer*4    mpierr
   integer*4    ourjob(nprocs)
   logical      flag_stat

   if(flag_stat) then
!    if_main write(6,'(A)')                   'START: BAND STRUCTURE EVALUATION'
     if_main write(6,'(A)')                   '       JOB DISTRUBUTION :'
     call MPI_BARRIER(mpi_comm_earth, mpierr)
             write(6,'(A,I0,A,I0,A)')         '       ->cpuid(',myid,'): ', ourjob(myid+1),' k-points'
     call MPI_BARRIER(mpi_comm_earth, mpierr)
   endif

   return
endsubroutine

subroutine allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA)
   use parameters
   implicit none
   type (incar)   :: PINPT       ! parameters for input arguments
   type (energy)  :: ETBA        ! target energy to be fitted to
   type (poscar)  :: PGEOM       ! parameters for geometry info
   type (kpoints) :: PKPTS       ! parameters for kpoints
   ! nband : number of eigenvalues to be stored for each spin
   ! neig  : number of orbital basis
   ! nspin : 2 for collinear 1 for non-collinear
   ! ispin : 2 for collinear 2 for non-collinear


   allocate(ETBA%E(PINPT%nband*PINPT%nspin, PKPTS%nkpoint))
   allocate(ETBA%V(PGEOM%neig*PINPT%ispin,PINPT%nband*PINPT%nspin, PKPTS%nkpoint))

return
endsubroutine
subroutine initialize_eig_status(ii, iadd, stat, nkpoint)
   use mpi_setup
   implicit none
   character*100 stat
   integer*4     iadd, ii, nkpoint
    
   stat = '****************************************************************************************************'
   if(myid .eq. 0) write(6,'(A)')stat
   if(nkpoint .le. 25) then
!    iadd=floor(real(1)/real(nkpoint)*100d0)
     iadd=10
     ii = 1
   else
     iadd=5
     ii = 1
   endif
!#write(6,*)"XX ", iadd
!stop
return
endsubroutine

subroutine print_eig_status(ik, ii, iadd, stat, nkpoint, flag_stat)
   use mpi_setup
   implicit none
   character*100   stat
   integer*4       ik, ii, iadd, nkpoint
   real*8          percent
   logical         flag_stat

   if(.not.flag_stat) return
   if(myid .ne. 0) return

   percent =  ik / real(nkpoint) * 100d0
   if( floor(percent) .ge. real(iadd*ii) ) then
     if_main write(6,'(A,I3)')stat(1:iadd*ii),floor(percent)
     ii = ii + 1
   endif

return
endsubroutine

subroutine get_ham_Hk(Hk, NN_TABLE, PINPT, kpoint, is, neig, flag_phase)
   use parameters, only: hopping, incar
   use mpi_setup
   use phase_factor
   use do_math
   implicit none
   type (hopping) :: NN_TABLE
   type (incar  ) :: PINPT
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

   call get_hamk_dense(Hk, H0, Hm, Hs, is, kpoint, PINPT, neig, NN_TABLE, flag_init, flag_phase)

   return
endsubroutine

subroutine stop_get_eig(msize, nband)
   use mpi_setup
   implicit none
   integer*4    msize, nband

   if_main write(6,'(A)')        '    !WARN! Check NERANGE! NBAND should be less equal to matrix size MSIZE (NBAND <= MSIZE = N_ORBIT*ISPINOR),' 
   if_main write(6,'(A)')        '           where ISPINOR = 2/1 if LSORB = .TRUE./.FALSE. and N_ORBIT = total number of orbitals.' 
   if_main write(6,'(A,I0,A,I0)')'           MSIZE = ',msize, ' , NBAND = FINA_E - INIT_E + 1 = ',nband
   if_main write(6,'(A)')        '           Exit program...'
   stop

return
endsubroutine

subroutine cal_eig(Hk, neig, ispinor, ispin, iband, nband, E, V, flag_vector)
   use do_math
   use mpi_setup
   implicit none
   integer*4, intent(in):: iband, nband, neig, ispinor, ispin
!  integer*4              neig, iband, nband, ispinor, ispin, msize
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
subroutine cal_eig_sparse(SHk, neig, ispinor, ispin, nemax, ne_guess, E, V, flag_vector, &
                          emin, emax, ne_found, feast_fpm, feast_info, ne_prev)
   use parameters, only : spmat
   use do_math
   use mpi_setup
   implicit none
   type(spmat):: SHk
   integer*4     neig, nemax, ispinor, ispin, ne_prev
   complex*16    V(neig*ispinor,nemax)
   real*8        E(nemax)
   real*8        emin, emax
   integer*4     ne_found, ne_guess, feast_info
   integer*4     feast_fpm(128)
   logical       flag_vector

   call cal_eig_hermitianx_sparse(SHk, emin, emax, nemax, ne_found, ne_guess, E, V, flag_vector,&
                                  feast_fpm, feast_info, ne_prev)

return
endsubroutine
subroutine adjust_ne_guess(feast_info, is, ne_found, kp, ik, neig, nemax, PINPT)
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
