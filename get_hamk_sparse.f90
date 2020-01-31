#include "alias.inc"
#ifdef MKL_SPARSE
subroutine get_hamk_sparse_overlap(SHk, SSk, SH0,SS0, SHm, SHs, is, kp, PINPT, neig, NN_TABLE, flag_init, flag_phase, &
                           flag_sparse_zero_SHm, flag_sparse_zero_SHs)
  use parameters, only : incar, hopping, spmat
  use mpi_setup
  use phase_factor
  use do_math
  use kronecker_prod
  use sparse_tool
#ifdef MKL_SPARSE
  use MKL_SPBLAS
#endif
  use time
  implicit none
  type (incar  ) :: PINPT
  type (hopping) :: NN_TABLE
  type (spmat  ) :: SHk, SH0, SHm, SHs
  type (spmat  ) :: SSk, SS0  ! for overlap matrix
#ifdef MKL_SPARSE
  type (SPARSE_MATRIX_T) ::  Sk,  Sk_  , S0, Sm, Ss
  type (SPARSE_MATRIX_T) ::S_Sk,S_Sk_, S_S0 ! for overlap matrix
#endif
  logical           flag_init, flag_phase
  integer*4         is, neig, ispinor, ispin
  integer*4         istat
  real*8            kp(3)
  complex*16        alpha
  logical           flag_sparse_zero_SHm ! if the sparse matrix has no non-zero element (nnz = 0),
  logical           flag_sparse_zero_SHs ! this flag will be .true. and will skip related construction routine.
  real*8            t1, t0
  integer*4         mpierr

  ! DEALLOCATION of array
  ! Hk: initialized every call (deallocated with cal_eig_Hk_sparse exit)
  ! H0: initialized every call (deallocated with cal_eig_Hk_sparse exit) ! same rule holds for overlap matrix construction
  ! Hm: initialized in the first call (deallocated after get_eig exit)
  ! Hs: initialized in the first call if slater-koster (deallocated after get_eig exit)
  !     initialized every call        if .not. slater-koster (deallocated after cal_eig_Hk_sparse exit)
 

  ! H0 will be constructed for spin-1. For spin-2, copied from spin-1 
  if(PINPT%flag_use_overlap) then
    if(is .eq. 1) call set_ham0_sparse_overlap  (SH0, SS0, kp, PINPT, neig, NN_TABLE, flag_phase)
  elseif(.not. PINPT%flag_use_overlap) then
    if(is .eq. 1) call set_ham0_sparse          (SH0,      kp, PINPT, neig, NN_TABLE, flag_phase)
  endif
! if(PINPT%flag_use_overlap) then
!   if(is .eq. 1) call set_ham0_sparse (SS0, kp, PINPT, neig, NN_TABLE, flag_phase, .true. )
! endif

  if(flag_init) then   ! setup k-independent Hamiltonian: Hm, Hs (if .not. slater_koster)
    if(PINPT%flag_collinear) then
      call set_ham_mag_sparse(SHm, NN_TABLE, PINPT, neig, flag_sparse_zero_SHm)
    elseif(PINPT%flag_noncollinear) then
      call set_ham_mag_sparse(SHm, NN_TABLE, PINPT, neig, flag_sparse_zero_SHm)

      if(PINPT%flag_soc .and. PINPT%flag_slater_koster) then
        call set_ham_soc_sparse(SHs, 0d0, PINPT, neig, NN_TABLE, flag_phase, flag_sparse_zero_SHs, F_IJ)
      endif
    endif
    flag_init = .false.
  endif

  if(PINPT%flag_collinear) then
    if(.not. flag_sparse_zero_SHm) then
      call sparse_create_csr_handle(S0, SH0)
      if(PINPT%flag_use_overlap) then
        call sparse_create_csr_handle(S_Sk, SS0) ! since overlap matrix does not need to add Sm, SS0 is directly handled by S_Sk
      endif
      call sparse_create_csr_handle(Sm, SHm)
      alpha = ((-1d0)**(is+1))

      ! Sk(up) = S0 + Sm , Sk(dn) = S0 - Sm
      istat = MKL_SPARSE_z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0 , alpha, Sm, Sk)
      call sparse_error_report('MKL_SPARSE_z_ADD: S0+Sm=Sk', istat)

      call sparse_export_csr(Sk, SHk)
      if(PINPT%flag_use_overlap) then
        call sparse_export_csr(S_Sk, SSk)
      endif

      istat = MKL_SPARSE_DESTROY(Sm)
      call sparse_error_report('MKL_SPARSE_DESTROY: Sm ', istat)
      istat = MKL_SPARSE_DESTROY(S0)
      call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
!     if(PINPT%flag_use_overlap) then
!       istat = MKL_SPARSE_DESTROY(S_S0)
!       call sparse_error_report('MKL_SPARSE_DESTROY: S_S0 ', istat)
!     endif
    elseif(flag_sparse_zero_SHm) then
      SHk%msize = SH0%msize
      SHk%nnz   = SH0%nnz
      allocate(SHk%H(SH0%nnz    )) ; SHk%H = SH0%H
      allocate(SHk%J(SH0%nnz    )) ; SHk%J = SH0%J
      allocate(SHk%I(SH0%msize+1)) ; SHk%I = SH0%I
      if(PINPT%flag_use_overlap) then
        SSk%msize = SS0%msize
        SSk%nnz   = SS0%nnz
        allocate(SSk%H(SS0%nnz    )) ; SSk%H = SS0%H
        allocate(SSk%J(SS0%nnz    )) ; SSk%J = SS0%J
        allocate(SSk%I(SS0%msize+1)) ; SSk%I = SS0%I
      endif
    endif

  elseif(PINPT%flag_noncollinear) then
    if(PINPT%flag_soc) then
      if(.not. PINPT%flag_slater_koster) then
        !set up k-dependent SOC in the case of 'cc' orbitals
        call set_ham_soc_sparse(SHs, kp, PINPT, neig, NN_TABLE, flag_phase, flag_sparse_zero_SHs, F_IJ)
!       if_main write(6,'(A)')'    !WARN! Current version does not support Sparse Matrix '
!       if_main write(6,'(A)')'           for non-Slater-Koster type Hamiltonian'
!       if_main write(6,'(A)')'           Exit program...'
!       kill_job
      endif

      call kproduct_pauli_0_CSR(SH0)
      call sparse_create_csr_handle(S0, SH0)
      if(PINPT%flag_use_overlap) then
        call kproduct_pauli_0_CSR(SS0)
        call sparse_create_csr_handle(S_Sk, SS0)
      endif
      if(.not. flag_sparse_zero_SHm) call sparse_create_csr_handle(Sm, SHm)
      if(.not. flag_sparse_zero_SHs) call sparse_create_csr_handle(Ss, SHs)
      alpha = 1d0

      if(.not. flag_sparse_zero_SHm) then
        istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0  , alpha, Sm, Sk_)
        call sparse_error_report('MKL_SPARSE_z_ADD: S0+Sm=Sk_', istat)
        istat = MKL_SPARSE_DESTROY(Sm)
        if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sm ', istat)
        istat = MKL_SPARSE_DESTROY(S0)
        if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
      endif

      if(.not. flag_sparse_zero_SHm) then
        if(.not. flag_sparse_zero_SHs) then
          istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, Sk_ , alpha, Ss, Sk)
          call sparse_error_report('MKL_SPARSE_z_ADD: Sk_+Ss=Sk', istat)
          istat = MKL_SPARSE_DESTROY(Sk_)
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sk_', istat)
          istat = MKL_SPARSE_DESTROY(Ss)
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Ss ', istat)
        endif
      else
        if(.not. flag_sparse_zero_SHs) then
          istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0  , alpha, Ss, Sk)
          call sparse_error_report('MKL_SPARSE_z_ADD: S0 +Ss=Sk', istat)
          istat = MKL_SPARSE_DESTROY(S0)
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
          istat = MKL_SPARSE_DESTROY(Ss)
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Ss ', istat)
        endif
      endif

      if(.not. flag_sparse_zero_SHm) then
        if(.not. flag_sparse_zero_SHs) then
          call sparse_export_csr(Sk, SHk)
        else
          call sparse_export_csr(Sk_, SHk)
        endif
      elseif(  flag_sparse_zero_SHm) then
        if(.not. flag_sparse_zero_SHs) then
          call sparse_export_csr(Sk, SHk)
        elseif(  flag_sparse_zero_SHs) then
          call sparse_export_csr(S0, SHk)
        endif
      endif

      if(PINPT%flag_use_overlap) then
        call sparse_export_csr(S_Sk, SSk)
      endif

    else
      call kproduct_pauli_0_CSR(SH0)
      call sparse_create_csr_handle(S0, SH0)
      if(PINPT%flag_use_overlap) then
        call kproduct_pauli_0_CSR(SS0)
        call sparse_create_csr_handle(S_Sk, SS0)
      endif

      if(.not. flag_sparse_zero_SHm) call sparse_create_csr_handle(Sm, SHm)
      alpha = 1d0

      if(.not. flag_sparse_zero_SHm) then
        istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0  , alpha, Sm, Sk)
        call sparse_error_report('MKL_SPARSE_z_ADD: S0+Sm=Sk', istat)
      endif

      if(.not. flag_sparse_zero_SHm) then
        call sparse_export_csr(Sk, SHk)
      else
        call sparse_export_csr(S0, SHk)
      endif

      if(PINPT%flag_use_overlap) then
        call sparse_export_csr(S_Sk, SSk)
      endif

      istat = MKL_SPARSE_DESTROY(Sm)
      if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sm ', istat)
      istat = MKL_SPARSE_DESTROY(S0)
      if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
    endif

  elseif(.not. PINPT%flag_collinear .and. .not. PINPT%flag_noncollinear) then

    allocate(SHk%H(SH0%nnz))
    allocate(SHk%J(SH0%nnz))
    allocate(SHk%I(SH0%msize + 1))
    SHk%H     = SH0%H
    SHk%I     = SH0%I
    SHk%J     = SH0%J
    SHk%nnz   = SH0%nnz
    SHk%msize = SH0%msize

    if(PINPT%flag_use_overlap) then
      allocate(SSk%H(SS0%nnz))
      allocate(SSk%J(SS0%nnz))
      allocate(SSk%I(SS0%msize + 1))
      SSk%H     = SS0%H
      SSk%I     = SS0%I
      SSk%J     = SS0%J
      SSk%nnz   = SS0%nnz
      SSk%msize = SS0%msize
    endif

  endif

  return
endsubroutine

subroutine get_hamk_sparse(SHk, SH0, SHm, SHs, is, kp, PINPT, neig, NN_TABLE, flag_init, flag_phase, &
                           flag_sparse_zero_SHm, flag_sparse_zero_SHs)
  use parameters, only : incar, hopping, spmat
  use mpi_setup
  use phase_factor
  use do_math
  use kronecker_prod
  use sparse_tool
#ifdef MKL_SPARSE
  use MKL_SPBLAS
#endif
  use time
  implicit none
  type (incar  ) :: PINPT
  type (hopping) :: NN_TABLE
  type (spmat  ) :: SHk, SH0, SHm, SHs
#ifdef MKL_SPARSE
  type (SPARSE_MATRIX_T) :: Sk, Sk_, S0, Sm, Ss
#endif
  logical           flag_init, flag_phase  
  integer*4         is, neig, ispinor, ispin
  integer*4         istat
  real*8            kp(3)  
  complex*16        alpha
  logical           flag_sparse_zero_SHm ! if the sparse matrix has no non-zero element (nnz = 0),
  logical           flag_sparse_zero_SHs ! this flag will be .true. and will skip related construction routine.
  real*8            t1, t0
  integer*4         mpierr

  ! DEALLOCATION of array
  ! Hk: initialized every call (deallocated with cal_eig_Hk_sparse exit)
  ! H0: initialized every call (deallocated with cal_eig_Hk_sparse exit)
  ! Hm: initialized in the first call (deallocated after get_eig exit)
  ! Hs: initialized in the first call if slater-koster (deallocated after get_eig exit)
  !     initialized every call        if .not. slater-koster (deallocated after cal_eig_Hk_sparse exit)

  ! H0 will be constructed for spin-1. For spin-2, copied from spin-1 
  if(is .eq. 1) call set_ham0_sparse   (SH0, kp, PINPT, neig, NN_TABLE, flag_phase)

  if(flag_init) then   ! setup k-independent Hamiltonian: Hm, Hs (if .not. slater_koster)
    if(PINPT%flag_collinear) then
      call set_ham_mag_sparse(SHm, NN_TABLE, PINPT, neig, flag_sparse_zero_SHm)
    elseif(PINPT%flag_noncollinear) then
      call set_ham_mag_sparse(SHm, NN_TABLE, PINPT, neig, flag_sparse_zero_SHm)

      if(PINPT%flag_soc .and. PINPT%flag_slater_koster) then
        call set_ham_soc_sparse(SHs, 0d0, PINPT, neig, NN_TABLE, flag_phase, flag_sparse_zero_SHs, F_IJ)
      endif
    endif
    flag_init = .false.
  endif

  if(PINPT%flag_collinear) then
    if(.not. flag_sparse_zero_SHm) then
      call sparse_create_csr_handle(S0, SH0)
      call sparse_create_csr_handle(Sm, SHm)
      alpha = ((-1d0)**(is+1))

      ! Sk(up) = S0 + Sm , Sk(dn) = S0 - Sm
      istat = MKL_SPARSE_z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0 , alpha, Sm, Sk)
      call sparse_error_report('MKL_SPARSE_z_ADD: S0+Sm=Sk', istat)

      call sparse_export_csr(Sk, SHk)

      istat = MKL_SPARSE_DESTROY(Sm)
      call sparse_error_report('MKL_SPARSE_DESTROY: Sm ', istat)
      istat = MKL_SPARSE_DESTROY(S0)
      call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
    elseif(flag_sparse_zero_SHm) then
      SHk%msize = SH0%msize
      SHk%nnz   = SH0%nnz
      allocate(SHk%H(SH0%nnz    )) ; SHk%H = SH0%H
      allocate(SHk%J(SH0%nnz    )) ; SHk%J = SH0%J
      allocate(SHk%I(SH0%msize+1)) ; SHk%I = SH0%I
    endif

  elseif(PINPT%flag_noncollinear) then
    if(PINPT%flag_soc) then
      if(.not. PINPT%flag_slater_koster) then 
        !set up k-dependent SOC in the case of 'cc' orbitals
        call set_ham_soc_sparse(SHs, kp, PINPT, neig, NN_TABLE, flag_phase, flag_sparse_zero_SHs, F_IJ)
!       if_main write(6,'(A)')'    !WARN! Current version does not support Sparse Matrix '
!       if_main write(6,'(A)')'           for non-Slater-Koster type Hamiltonian'
!       if_main write(6,'(A)')'           Exit program...'
!       kill_job
      endif
      
      call kproduct_pauli_0_CSR(SH0)
      call sparse_create_csr_handle(S0, SH0)
      if(.not. flag_sparse_zero_SHm) call sparse_create_csr_handle(Sm, SHm)
      if(.not. flag_sparse_zero_SHs) call sparse_create_csr_handle(Ss, SHs)
      alpha = 1d0
     
      if(.not. flag_sparse_zero_SHm) then
        istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0  , alpha, Sm, Sk_)    
        call sparse_error_report('MKL_SPARSE_z_ADD: S0+Sm=Sk_', istat)
        istat = MKL_SPARSE_DESTROY(Sm)    
        if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sm ', istat)
        istat = MKL_SPARSE_DESTROY(S0)    
        if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
      endif

      if(.not. flag_sparse_zero_SHm) then
        if(.not. flag_sparse_zero_SHs) then
          istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, Sk_ , alpha, Ss, Sk)
          call sparse_error_report('MKL_SPARSE_z_ADD: Sk_+Ss=Sk', istat)
          istat = MKL_SPARSE_DESTROY(Sk_)
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sk_', istat)
          istat = MKL_SPARSE_DESTROY(Ss)    
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Ss ', istat)
        endif
      else
        if(.not. flag_sparse_zero_SHs) then
          istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0  , alpha, Ss, Sk)
          call sparse_error_report('MKL_SPARSE_z_ADD: S0 +Ss=Sk', istat)
          istat = MKL_SPARSE_DESTROY(S0)    
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
          istat = MKL_SPARSE_DESTROY(Ss)    
          if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Ss ', istat)
        endif
      endif

      if(.not. flag_sparse_zero_SHm) then
        if(.not. flag_sparse_zero_SHs) then
          call sparse_export_csr(Sk, SHk)
        else
          call sparse_export_csr(Sk_, SHk)
        endif
      elseif(  flag_sparse_zero_SHm) then
        if(.not. flag_sparse_zero_SHs) then
          call sparse_export_csr(Sk, SHk)
        elseif(  flag_sparse_zero_SHs) then
          call sparse_export_csr(S0, SHk)
        endif
      endif

    else
      call kproduct_pauli_0_CSR(SH0)
      call sparse_create_csr_handle(S0, SH0)
      if(.not. flag_sparse_zero_SHm) call sparse_create_csr_handle(Sm, SHm)
      alpha = 1d0

      if(.not. flag_sparse_zero_SHm) then
        istat = MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, S0  , alpha, Sm, Sk)
        call sparse_error_report('MKL_SPARSE_z_ADD: S0+Sm=Sk', istat)
      endif

      if(.not. flag_sparse_zero_SHm) then
        call sparse_export_csr(Sk, SHk)
      else
        call sparse_export_csr(S0, SHk)
      endif

      istat = MKL_SPARSE_DESTROY(Sm)
      if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sm ', istat)
      istat = MKL_SPARSE_DESTROY(S0)
      if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
    endif

  elseif(.not. PINPT%flag_collinear .and. .not. PINPT%flag_noncollinear) then

    allocate(SHk%H(SH0%nnz))
    allocate(SHk%J(SH0%nnz))
    allocate(SHk%I(SH0%msize + 1))
 
    SHk%H     = SH0%H
    SHk%I     = SH0%I
    SHk%J     = SH0%J
    SHk%nnz   = SH0%nnz
    SHk%msize = SH0%msize

  endif

  return
endsubroutine
subroutine set_ham0_sparse_overlap(SH0, SS0, kpoint, PINPT, neig, NN_TABLE, flag_phase)
  use parameters, only : zi, hopping, incar, spmat, eta
  use phase_factor
  use mpi_setup
  use sparse_tool
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (spmat  ) :: SH0_COO, SH0
  type (spmat  ) :: SS0_COO, SS0
  integer*4         neig , ii, jj, nn, m, mm
  integer*4         ivel_axis
  real*8            kpoint(3), tol
  integer*4         I(NN_TABLE%n_neighbor*2) ! default maximum : n_neighbor
  integer*4         J(NN_TABLE%n_neighbor*2)
  complex*16        H(NN_TABLE%n_neighbor*2)
  complex*16        S(NN_TABLE%n_neighbor*2) ! for overlap
  complex*16        IJ(NN_TABLE%n_neighbor*2)
  complex*16        Eij, Tij
  complex*16        Sij
  logical           flag_phase
  real*8            abstol

  H=0.0d0
  S=0.0d0
  IJ= (-1d0,-1d0)
  mm = 0
  tol=NN_TABLE%onsite_tolerance

nn_:do nn=1,NN_TABLE%n_neighbor
    ii=NN_TABLE%i_matrix(nn)
    jj=NN_TABLE%j_matrix(nn)
    call get_hopping_integral(Eij, NN_TABLE, nn, PINPT, tol, kpoint, F_IJ, flag_phase, .false.)
    call get_hopping_integral(Sij, NN_TABLE, nn, PINPT, tol, kpoint, F_IJ, flag_phase, .true. )

    if(ii .eq. jj .and. NN_TABLE%Dij(nn) <= tol) then
      if(nint(PINPT%param_const(4,NN_TABLE%local_U_param_index(ii))) .ge. 1) then
        Tij =  Eij + NN_TABLE%local_charge(ii)*PINPT%param_const(5,(NN_TABLE%local_U_param_index(ii)))
      else
        Tij =  Eij + NN_TABLE%local_charge(ii)*PINPT%param((NN_TABLE%local_U_param_index(ii)))
      endif

      if(abs(Tij) .ge. eta) then
        do m = 1, mm ! check previously stored array: if exist, overwrite Tij (Tii), otherwise, save new array
          if( ii .eq. nint(real(IJ(m))) .and. jj .eq. nint(aimag(IJ(m))) ) then
            H(m) = H(m) + Tij
            S(m) = S(m) + 1d0 ! if onsite, overlap integral should be 1d0.
            cycle nn_
          endif
        enddo
        mm = mm + 1
        I(mm) = ii
        J(mm) = jj
        H(mm) = Tij
        S(mm) = 1d0 ! if onsite, overlap integral should be 1d0.
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi
      endif
    elseif(ii .eq. jj .and. NN_TABLE%Dij(nn) > tol) then
      Tij =  Eij

      if(abs(Tij) .ge. eta) then
        do m = 1, mm ! check previously stored array: if exist, overwrite Tij (Tii), otherwise, save new array
          if( ii .eq. nint(real(IJ(m))) .and. jj .eq. nint(aimag(IJ(m))) ) then
            H(m) = H(m) + Tij
            S(m) = S(m) + Sij
            cycle nn_
          endif
        enddo
        mm = mm + 1
        I(mm) = ii
        J(mm) = jj
        H(mm) = Tij
        S(mm) = Sij
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi
      endif

    else
      Tij = Eij

      if(abs(Tij) .ge. eta) then
        do m = 1, mm ! check previously stored array: if exist, overwrite Tij and Tji, otherwise, save new array
          if( ii .eq. nint(real(IJ(m))) .and. jj .eq. nint(aimag(IJ(m))) ) then
            H(m)   = H(m) + Tij
            H(m+1) = H(m+1) + conjg(Tij)
            S(m)   = S(m) + Sij
            S(m+1) = S(m+1) + conjg(Sij)
            cycle nn_
          endif
        enddo
        mm = mm + 1
        I(mm) = ii
        J(mm) = jj
        H(mm) = Tij
        S(mm) = Sij
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi
        mm = mm + 1
        I(mm) = jj
        J(mm) = ii
        H(mm) = conjg(Tij)
        S(mm) = conjg(Sij)
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi

        if(mm .gt. NN_TABLE%n_neighbor*2) then
          if_main write(6,'(A)')'    !WARN! Number of non-zero element NNZ > n_neighbor*2.'
          if_main write(6,'(A)')'           Please check "set_ham0_sparse" routine. Exit program...'
          stop
        endif
      endif

    endif

  enddo nn_

  if(mm .eq. 0) then
    if_main write(6,'(A)')'    !WARN! No non-zero element is found in preparing Hamiltonian H0.'
    if_main write(6,'(A)')'           Please check your system, input file etc. Exit program...'
    stop
  endif

  allocate(SH0_COO%H(mm))
  allocate(SH0_COO%I(mm))
  allocate(SH0_COO%J(mm))
  SH0_COO%nnz   = mm
  SH0_COO%msize = neig
  SH0_COO%H     = H(1:mm)
  SH0_COO%I     = I(1:mm)
  SH0_COO%J     = J(1:mm)
  call sparse_convert_coo_csr(SH0_COO, SH0)

  allocate(SS0_COO%H(mm))
  allocate(SS0_COO%I(mm))
  allocate(SS0_COO%J(mm))
  SS0_COO%nnz   = mm
  SS0_COO%msize = neig
  SS0_COO%H     = S(1:mm)
  SS0_COO%I     = I(1:mm)
  SS0_COO%J     = J(1:mm)
  call sparse_convert_coo_csr(SS0_COO, SS0)

  deallocate(SH0_COO%H)
  deallocate(SH0_COO%I)
  deallocate(SH0_COO%J)

return
endsubroutine

subroutine set_ham0_sparse(SH0, kpoint, PINPT, neig, NN_TABLE, flag_phase)
  use parameters, only : zi, hopping, incar, spmat, eta
  use phase_factor
  use mpi_setup
  use sparse_tool
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (spmat  ) :: SH0_COO, SH0
  integer*4         neig , ii, jj, nn, m, mm
  integer*4         ivel_axis
  real*8            kpoint(3), tol
  integer*4         I(NN_TABLE%n_neighbor*2) ! default maximum : n_neighbor
  integer*4         J(NN_TABLE%n_neighbor*2)
  complex*16        H(NN_TABLE%n_neighbor*2)
  complex*16        IJ(NN_TABLE%n_neighbor*2)
  complex*16        Eij, Tij
  logical           flag_phase
  real*8            abstol

  H=0.0d0
  IJ= (-1d0,-1d0)
  mm = 0
  tol=NN_TABLE%onsite_tolerance

nn_:do nn=1,NN_TABLE%n_neighbor
    ii=NN_TABLE%i_matrix(nn)
    jj=NN_TABLE%j_matrix(nn)
    call get_hopping_integral(Eij, NN_TABLE, nn, PINPT, tol, kpoint, F_IJ, flag_phase, .false.)

    if(ii .eq. jj .and. NN_TABLE%Dij(nn) <= tol) then
      if(nint(PINPT%param_const(4,NN_TABLE%local_U_param_index(ii))) .ge. 1) then
        Tij =  Eij + NN_TABLE%local_charge(ii)*PINPT%param_const(5,(NN_TABLE%local_U_param_index(ii)))
      else
        Tij =  Eij + NN_TABLE%local_charge(ii)*PINPT%param((NN_TABLE%local_U_param_index(ii)))
      endif

      if(abs(Tij) .ge. eta) then
        do m = 1, mm ! check previously stored array: if exist, overwrite Tij (Tii), otherwise, save new array
          if( ii .eq. nint(real(IJ(m))) .and. jj .eq. nint(aimag(IJ(m))) ) then
            H(m) = H(m) + Tij
            cycle nn_
          endif
        enddo
        mm = mm + 1
        I(mm) = ii
        J(mm) = jj
        H(mm) = Tij
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi
      endif
    elseif(ii .eq. jj .and. NN_TABLE%Dij(nn) > tol) then
      Tij =  Eij

      if(abs(Tij) .ge. eta) then
        do m = 1, mm ! check previously stored array: if exist, overwrite Tij (Tii), otherwise, save new array
          if( ii .eq. nint(real(IJ(m))) .and. jj .eq. nint(aimag(IJ(m))) ) then
            H(m) = H(m) + Tij
            cycle nn_
          endif
        enddo
        mm = mm + 1
        I(mm) = ii
        J(mm) = jj
        H(mm) = Tij
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi
      endif

    else
      Tij = Eij

      if(abs(Tij) .ge. eta) then
        do m = 1, mm ! check previously stored array: if exist, overwrite Tij and Tji, otherwise, save new array
          if( ii .eq. nint(real(IJ(m))) .and. jj .eq. nint(aimag(IJ(m))) ) then
            H(m)   = H(m) + Tij
            H(m+1) = H(m+1) + conjg(Tij)
            cycle nn_
          endif
        enddo
        mm = mm + 1
        I(mm) = ii
        J(mm) = jj
        H(mm) = Tij
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi
        mm = mm + 1
        I(mm) = jj
        J(mm) = ii
        H(mm) = conjg(Tij)
        IJ(mm) = real(I(mm)) + real(J(mm)) * zi

        if(mm .gt. NN_TABLE%n_neighbor*2) then
          if_main write(6,'(A)')'    !WARN! Number of non-zero element NNZ > n_neighbor*2.'
          if_main write(6,'(A)')'           Please check "set_ham0_sparse" routine. Exit program...'
          stop
        endif
      endif

    endif

  enddo nn_

  if(mm .eq. 0) then
    if_main write(6,'(A)')'    !WARN! No non-zero element is found in preparing Hamiltonian H0.'
    if_main write(6,'(A)')'           Please check your system, input file etc. Exit program...'
    stop
  endif

  allocate(SH0_COO%H(mm))
  allocate(SH0_COO%I(mm))
  allocate(SH0_COO%J(mm))
  SH0_COO%nnz   = mm
  SH0_COO%msize = neig
  SH0_COO%H     = H(1:mm)
  SH0_COO%I     = I(1:mm)
  SH0_COO%J     = J(1:mm)

  call sparse_convert_coo_csr(SH0_COO, SH0)

  deallocate(SH0_COO%H)
  deallocate(SH0_COO%I)
  deallocate(SH0_COO%J)

return
endsubroutine
subroutine set_ham_mag_sparse(SHm, NN_TABLE, PINPT, neig, flag_sparse_zero)
  use parameters, only : zi, hopping, incar, spmat, eta
  use mpi_setup
  use sparse_tool
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (spmat  ) :: SHm_COO, SHm
  integer*4         neig
  integer*4         nn, ii, mm, m
  complex*16        Tij, Tij_x, Tij_y, Tij_z
  integer*4         I(NN_TABLE%n_neighbor*4) ! default maximum : n_neighbor*4
  integer*4         J(NN_TABLE%n_neighbor*4)
  complex*16        H(NN_TABLE%n_neighbor*4)
  complex*16        IJ(NN_TABLE%n_neighbor*4) 
  logical           flag_sparse_zero

    flag_sparse_zero = .false.

    H=0d0
    mm = 0
    IJ = (-1d0,-1d0)

    if(PINPT%flag_collinear) then
 nn_c:do nn = 1, neig
        if(NN_TABLE%stoner_I_param_index(nn) .gt. 0) then   ! if stoner parameter has been set...
          if(nint(PINPT%param_const(4,NN_TABLE%stoner_I_param_index(nn))) .eq. 1) then ! if i-th basis has constraint .true.
            Tij = -0.5d0 * NN_TABLE%local_moment(1,nn) * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(nn))
          else
            Tij = -0.5d0 * NN_TABLE%local_moment(1,nn) * PINPT%param(NN_TABLE%stoner_I_param_index(nn))
          endif
          
          if(abs(Tij) .ge. eta) then
            do m = 1, mm ! check previously stored array: if exist, overwrite Tij (Tii), otherwise, save new array
              if( nn .eq. nint(real(IJ(m))) ) then
                H(m) = H(m) + Tij
                cycle nn_c
              endif
            enddo
            mm = mm + 1
            I(mm) = nn
            J(mm) = nn
            H(mm) = Tij
            IJ(mm) = real(I(mm))

            if(mm .gt. NN_TABLE%n_neighbor*4) then
              if_main write(6,'(A)')'    !WARN! Number of non-zero element NNZ > n_neighbor*4.'
              if_main write(6,'(A)')'           Please check "set_ham_mag_sparse" routine. Exit program...'
              stop
            endif
          endif

        endif
      enddo nn_c

    elseif(PINPT%flag_noncollinear) then
nn_nc:do nn = 1, neig
        if(NN_TABLE%stoner_I_param_index(nn) .gt. 0) then   ! if stoner parameter has been set...
          if(nint(PINPT%param_const(4,NN_TABLE%stoner_I_param_index(nn))) .eq. 1) then ! if i-th basis has constraint .true.
            Tij_x  = - 0.5d0 * NN_TABLE%local_moment_rot(1,nn) * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(nn))
            Tij_y  =   0.5d0 * NN_TABLE%local_moment_rot(2,nn) * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(nn)) * zi
            Tij_z  = - 0.5d0 * NN_TABLE%local_moment_rot(3,nn) * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(nn))
          else
            Tij_x  = - 0.5d0 * NN_TABLE%local_moment_rot(1,nn) * PINPT%param(NN_TABLE%stoner_I_param_index(nn))
            Tij_y  =   0.5d0 * NN_TABLE%local_moment_rot(2,nn) * PINPT%param(NN_TABLE%stoner_I_param_index(nn)) * zi
            Tij_z  = - 0.5d0 * NN_TABLE%local_moment_rot(3,nn) * PINPT%param(NN_TABLE%stoner_I_param_index(nn))
          endif

          if(abs(Tij_x) + abs(Tij_y) + abs(Tij_z) .ge. eta) then
            do m = 1, mm
              if( nn .eq. nint(real(IJ(m))) .and. nn .eq. nint(aimag(IJ(m))) ) then
                H(m)   = H(m  ) + Tij_z
                H(m+1) = H(m+1) - Tij_z
                H(m+2) = H(m+2) + Tij_x + Tij_y
                H(m+3) = H(m+3) + Tij_x - Tij_y
                cycle nn_nc
              endif
            enddo
            mm = mm + 1
             I(mm) = nn
             J(mm) = nn
             H(mm) = Tij_z
             IJ(mm) = real(I(mm)) + zi * real(J(mm))
            mm = mm + 1
             I(mm) = nn + neig
             J(mm) = nn + neig
             H(mm) = -Tij_z
             IJ(mm) = real(I(mm)) + zi * real(J(mm))
            mm = mm + 1
             I(mm) = nn
             J(mm) = nn + neig
             H(mm) = Tij_x + Tij_y
             IJ(mm) = real(I(mm)) + zi * real(J(mm))
            mm = mm + 1
             I(mm) = nn  + neig
             J(mm) = nn
             H(mm) = Tij_x - Tij_y
             IJ(mm) = real(I(mm)) + zi * real(J(mm))
           
            if(mm .gt. NN_TABLE%n_neighbor*4) then
              if_main write(6,'(A)')'    !WARN! Number of non-zero element NNZ > n_neighbor*4.'
              if_main write(6,'(A)')'           Please check "set_ham_mag_sparse" routine. Exit program...'
              stop
            endif
          endif

        endif

      enddo nn_nc

    endif

    if(mm .eq. 0) then 
      flag_sparse_zero = .true.
    elseif(mm .ge. 1) then
      allocate(SHm_COO%H(mm))
      allocate(SHm_COO%I(mm))
      allocate(SHm_COO%J(mm))
      SHm_COO%nnz   = mm
      SHm_COO%msize = neig * PINPT%ispinor
      SHm_COO%H     = H(1:mm)
      SHm_COO%I     = I(1:mm)
      SHm_COO%J     = J(1:mm)

      call sparse_convert_coo_csr(SHm_COO, SHm)

      deallocate(SHm_COO%H)
      deallocate(SHm_COO%I)
      deallocate(SHm_COO%J)
    endif

return
endsubroutine
subroutine set_ham_soc_sparse(SHs, kp , PINPT, neig, NN_TABLE, flag_phase, flag_sparse_zero, FIJ)
  use parameters, only : zi, hopping, incar, spmat, eta, pi, pi2
  use sparse_tool
  use kronecker_prod
  use mpi_setup
  use phase_factor
  use get_parameter
#ifdef MKL_SPARSE
  use MKL_SPBLAS
#endif
  implicit none
  interface
    function FIJ(k,R)
      complex*16 :: FIJ
      real*8, intent(in) :: k(3)
      real*8, intent(in) :: R(3)
    endfunction
  end interface
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (spmat  ) :: COO_x, COO_y, COO_z, CSR_x, CSR_y, CSR_z, SHs
  type (SPARSE_MATRIX_T) :: SHx, SHy, SHz, SHxy, SHxyz
  integer*4         neig
  integer*4         nn, ii,jj, mm, m
  integer*4         istat
  integer*4         soc_index, rashba_index
  real*8            lambda_soc, lambda_rashba
  real*8            kp(3)
  logical           flag_phase, flag_sparse_zero
  complex*16        Tij, Tij_x, Tij_y, Tij_z
  complex*16        F
  integer*4         I(NN_TABLE%n_neighbor*2) ! default maximum : n_neighbor*4
  integer*4         J(NN_TABLE%n_neighbor*2)
  complex*16        Hx(NN_TABLE%n_neighbor*2)
  complex*16        Hy(NN_TABLE%n_neighbor*2)
  complex*16        Hz(NN_TABLE%n_neighbor*2)
  complex*16        IJ(NN_TABLE%n_neighbor*2)
  complex*16        alpha
  complex*16        L_x, L_y, L_z
  external          L_x, L_y, L_z
  character*8       ci_orb, cj_orb
  character*20      ci_atom, cj_atom
  real*8            lsign, hsign
  real*8            hop_signx, hop_signy, hop_signatom
  complex*16        prod

  flag_sparse_zero = .false.
  alpha = 1d0
  mm = 0
  IJ = (-1d0, -1d0)

    if(PINPT%flag_slater_koster) then
      Hx = 0d0
      Hy = 0d0
      Hz = 0d0

nn_sk:do nn = 1, NN_TABLE%n_neighbor
        soc_index = NN_TABLE%soc_param_index(nn)
        ii = NN_TABLE%i_matrix(nn) ; jj = NN_TABLE%j_matrix(nn)

        ! set SOC hamiltonian based on atomic orbitals or
        ! set SOC hamiltonian based on 'xx' type orbitals which are composed by linear combination of atomic orbitals
        if( soc_index .gt. 0 .and. (NN_TABLE%p_class(nn) .eq. 'pp' .or. NN_TABLE%p_class(nn) .eq. 'dd' &
                                                                   .or. NN_TABLE%p_class(nn) .eq. 'xx'  ) ) then
          call get_param(PINPT,    soc_index, 1, lambda_soc   )
         !call get_param(PINPT,    soc_index, lambda_soc   )
          ! CALCULATE  <orb_i|LS|orb_j> 
          Tij_x = lambda_soc * L_x(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))
          Tij_y = lambda_soc * L_y(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))
          Tij_z = lambda_soc * L_z(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))

          if(abs(Tij_x) + abs(Tij_y) + abs(Tij_z) .ge. eta) then
            mm = mm + 1
             I(mm) = ii
             J(mm) = jj
             Hx(mm) = Tij_x
             Hy(mm) = Tij_y
             Hz(mm) = Tij_z
            mm = mm + 1
             I(mm) = jj
             J(mm) = ii
             Hx(mm) = conjg(Tij_x)
             Hy(mm) = conjg(Tij_y)
             Hz(mm) = conjg(Tij_z)

            if(mm .gt. NN_TABLE%n_neighbor*4) then
              if_main write(6,'(A)')'    !WARN! Number of non-zero element NNZ > n_neighbor*4.'
              if_main write(6,'(A)')'           Please check "set_ham_soc_sparse" routine. Exit program...'
              stop
            endif
          endif

        endif

      enddo nn_sk

      if(mm .eq. 0) then
        flag_sparse_zero = .true. 
      elseif(mm .ge. 1) then
        allocate(COO_x%H(mm));allocate(COO_y%H(mm));allocate(COO_z%H(mm))
        allocate(COO_x%I(mm));allocate(COO_y%I(mm));allocate(COO_z%I(mm))
        allocate(COO_x%J(mm));allocate(COO_y%J(mm));allocate(COO_z%J(mm))
        COO_x%nnz = mm ; COO_x%msize = neig
        COO_y%nnz = mm ; COO_y%msize = neig
        COO_z%nnz = mm ; COO_z%msize = neig
        COO_x%H   = Hx(1:mm) ; COO_x%I   = I(1:mm) ; COO_x%J   = J(1:mm)
        COO_y%H   = Hy(1:mm) ; COO_y%I   = I(1:mm) ; COO_y%J   = J(1:mm)
        COO_z%H   = Hz(1:mm) ; COO_z%I   = I(1:mm) ; COO_z%J   = J(1:mm)


        call sparse_convert_coo_csr(COO_x, CSR_x)
        call sparse_convert_coo_csr(COO_y, CSR_y)
        call sparse_convert_coo_csr(COO_z, CSR_z)

        !SET UP Hamiltonian H_soc*sigma   
        call kproduct_pauli_x_CSR(CSR_x)
        call kproduct_pauli_y_CSR(CSR_y)
        call kproduct_pauli_z_CSR(CSR_z)
        call sparse_create_csr_handle(SHx, CSR_x)
        call sparse_create_csr_handle(SHy, CSR_y)
        call sparse_create_csr_handle(SHz, CSR_z)

        istat = MKL_SPARSE_z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, SHx , alpha, SHy, SHxy)
        call sparse_error_report('MKL_SPARSE_z_ADD: SHx+SHy=SHxy', istat)

        istat = MKL_SPARSE_z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, SHxy, alpha, SHz, SHxyz)
        call sparse_error_report('MKL_SPARSE_z_ADD: SHxy+SHz=SHxyz', istat)

        call sparse_export_csr(SHxyz, SHs)

        deallocate(COO_x%H);deallocate(COO_y%H);deallocate(COO_z%H)
        deallocate(COO_x%I);deallocate(COO_y%I);deallocate(COO_z%I)
        deallocate(COO_x%J);deallocate(COO_y%J);deallocate(COO_z%J)
        deallocate(CSR_x%H);deallocate(CSR_y%H);deallocate(CSR_z%H)
        deallocate(CSR_x%I);deallocate(CSR_y%I);deallocate(CSR_z%I)
        deallocate(CSR_x%J);deallocate(CSR_y%J);deallocate(CSR_z%J)

        istat = MKL_SPARSE_DESTROY(SHx)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHx', istat)
        istat = MKL_SPARSE_DESTROY(SHy)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHy', istat)
        istat = MKL_SPARSE_DESTROY(SHz)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHz', istat)
        istat = MKL_SPARSE_DESTROY(SHxy)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHxy', istat)
      endif

    elseif(.not.PINPT%flag_slater_koster) then
      Hx = 0d0
      Hy = 0d0
      Hz = 0d0

     !WARNING!! This setting is only valid for Bi/Si(110) case with certain atomic geometry and lattice vectors,
     !          since the sign convention is only valid and meaningful for this particular case.
     !          If you are dealing with other system, please construct your own hamltonian setup.
nn_cc:do nn = 1, NN_TABLE%n_neighbor
        soc_index    = NN_TABLE%cc_index_set(2,nn)
        rashba_index = NN_TABLE%cc_index_set(3,nn)
        ii = NN_TABLE%i_matrix(nn) ; jj = NN_TABLE%j_matrix(nn)

        if(flag_phase) then
          F = FIJ(kp, NN_TABLE%Rij(:,nn))
        elseif(.not. flag_phase) then
          F = FIJ(kp, NN_TABLE%R  (:,nn))
        endif

        if( soc_index .ge. 1 .and. rashba_index .ge. 1) then
          call get_param(PINPT,    soc_index, 1, lambda_soc   )
          call get_param(PINPT, rashba_index, 1, lambda_rashba)
         !call get_param(PINPT,    soc_index, lambda_soc   )
         !call get_param(PINPT, rashba_index, lambda_rashba)

!         ! set Rashba-SOC between i_orb and j_orb separated by |dij|, originated from E-field normal to surface
          Tij_x = zi * lambda_rashba * NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn) * F  ! sigma_x
          Tij_y = zi * lambda_rashba *-NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn) * F  ! sigma_y

          ! set SOC between i_orb and j_orb separated by |dij|, 
          ! originated from E-field due to neighbor atom nearby the hopping path
          ! H_SOC = sigma_<<ij>> i * lsoc * v_ij * ci' * sigma_z * cj
          ! v_ij = di x dj / (|di x dj|) , di(dj) are vector connecting nearest neigbohor
          ! atom from i (to j).
          hop_signx= NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn)
          hop_signy= NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn)
          ci_atom = NN_TABLE%site_cindex(NN_TABLE%i_atom(nn))
          cj_atom = NN_TABLE%site_cindex(NN_TABLE%j_atom(nn))
          if( ci_atom(1:2) .eq. 'b1') hop_signatom = 1.0
          if( ci_atom(1:2) .eq. 'b2') hop_signatom =-1.0
          Tij_z = zi * lambda_soc * hop_signatom * (hop_signx + hop_signy) * F

          call save_Hsoc_sparse(Tij_x, Tij_y, Tij_z, Hx, Hy, Hz, mm, ii, jj, I, J, IJ, NN_TABLE%n_neighbor)

        elseif( soc_index .ge. 1 .and. rashba_index .eq. 0) then
          call get_param(PINPT,    soc_index, 1, lambda_soc   )
         !call get_param(PINPT,    soc_index, lambda_soc   )

!         ! This model is only for Kane-mele type of SOC. Be careful..
          prod=exp(-2d0*zi * pi * dot_product((/2.45d0,0d0/), NN_TABLE%Rij(1:2,nn)))
          hsign  = sign(1d0,aimag(prod))
          ci_atom = NN_TABLE%site_cindex(NN_TABLE%i_atom(nn))
          cj_atom = NN_TABLE%site_cindex(NN_TABLE%j_atom(nn))
          if( ci_atom(1:1) .eq. 'a') lsign = -1.0d0
          if( ci_atom(1:1) .eq. 'b') lsign =  1.0d0

          Tij_x             = (0d0, 0d0)
          Tij_y             = (0d0, 0d0)
          Tij_z             = zi * lambda_soc * lsign * hsign * F
          call save_Hsoc_sparse(Tij_x, Tij_y, Tij_z, Hx, Hy, Hz, mm, ii, jj, I, J, IJ, NN_TABLE%n_neighbor)

        elseif( soc_index .eq. 0 .and. rashba_index .gt. 1 ) then ! WARN: only the AB-a hopping is considered (for Bi/Si110 case)
          call get_param(PINPT, rashba_index, 1, lambda_rashba)
         !call get_param(PINPT, rashba_index, lambda_rashba)

!         ! set Rashba-SOC between i_orb and j_orb separated by |dij|, originated from E-field normal to surface
          Tij_x = zi * lambda_rashba * NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn) * F  ! sigma_x
          Tij_y = zi * lambda_rashba *-NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn) * F  ! sigma_y
          Tij_z = (0d0,0d0)

          call save_Hsoc_sparse(Tij_x, Tij_y, Tij_z, Hx, Hy, Hz, mm, ii, jj, I, J, IJ, NN_TABLE%n_neighbor)
        endif

      enddo nn_cc

      if(mm .eq. 0) then
        flag_sparse_zero = .true.
      elseif(mm .ge. 1) then

        allocate(COO_x%H(mm));allocate(COO_y%H(mm));allocate(COO_z%H(mm))
        allocate(COO_x%I(mm));allocate(COO_y%I(mm));allocate(COO_z%I(mm))
        allocate(COO_x%J(mm));allocate(COO_y%J(mm));allocate(COO_z%J(mm))
        COO_x%nnz = mm ; COO_x%msize = neig
        COO_y%nnz = mm ; COO_y%msize = neig
        COO_z%nnz = mm ; COO_z%msize = neig
        COO_x%H   = Hx(1:mm) ; COO_x%I   = I(1:mm) ; COO_x%J   = J(1:mm)
        COO_y%H   = Hy(1:mm) ; COO_y%I   = I(1:mm) ; COO_y%J   = J(1:mm)
        COO_z%H   = Hz(1:mm) ; COO_z%I   = I(1:mm) ; COO_z%J   = J(1:mm)


        call sparse_convert_coo_csr(COO_x, CSR_x)
        call sparse_convert_coo_csr(COO_y, CSR_y)
        call sparse_convert_coo_csr(COO_z, CSR_z)

        !SET UP Hamiltonian H_soc*sigma   
        call kproduct_pauli_x_CSR(CSR_x)
        call kproduct_pauli_y_CSR(CSR_y)
        call kproduct_pauli_z_CSR(CSR_z)
        call sparse_create_csr_handle(SHx, CSR_x)
        call sparse_create_csr_handle(SHy, CSR_y)
        call sparse_create_csr_handle(SHz, CSR_z)

        istat = MKL_SPARSE_z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, SHx , alpha, SHy, SHxy)
        call sparse_error_report('MKL_SPARSE_z_ADD: SHx+SHy=SHxy', istat)

        istat = MKL_SPARSE_z_ADD(SPARSE_OPERATION_NON_TRANSPOSE, SHxy, alpha, SHz, SHxyz)
        call sparse_error_report('MKL_SPARSE_z_ADD: SHxy+SHz=SHxyz', istat)

        call sparse_export_csr(SHxyz, SHs)

        deallocate(COO_x%H);deallocate(COO_y%H);deallocate(COO_z%H)
        deallocate(COO_x%I);deallocate(COO_y%I);deallocate(COO_z%I)
        deallocate(COO_x%J);deallocate(COO_y%J);deallocate(COO_z%J)
        deallocate(CSR_x%H);deallocate(CSR_y%H);deallocate(CSR_z%H)
        deallocate(CSR_x%I);deallocate(CSR_y%I);deallocate(CSR_z%I)
        deallocate(CSR_x%J);deallocate(CSR_y%J);deallocate(CSR_z%J)

        istat = MKL_SPARSE_DESTROY(SHx)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHx', istat)
        istat = MKL_SPARSE_DESTROY(SHy)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHy', istat)
        istat = MKL_SPARSE_DESTROY(SHz)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHz', istat)
        istat = MKL_SPARSE_DESTROY(SHxy)
        call sparse_error_report('MKL_SPARSE_DESTROY: SHxy', istat)
      endif
       
    endif
return
endsubroutine
subroutine save_Hsoc_sparse(Tij_x, Tij_y, Tij_z, Hx, Hy, Hz, mm, ii, jj, I, J, IJ, n_neighbor)
   use parameters, only: zi, eta
   use mpi_setup
   implicit none
   complex*16   Tij_x, Tij_y, Tij_z
   complex*16   Hx(n_neighbor*2), Hy(n_neighbor*2), Hz(n_neighbor*2), IJ(n_neighbor*2)
   integer*4    I(n_neighbor*2), J(n_neighbor*2)
   integer*4    n_neighbor, ii, jj, m, mm
   integer*4    mpierr

   if(abs(Tij_x) + abs(Tij_y) + abs(Tij_z) .ge. eta) then
     do m = 1, mm ! check previously stored array: if exist, overwrite Tij and Tji, otherwise, save new array
       if( ii .eq. nint(real(IJ(m))) .and. jj .eq. nint(aimag(IJ(m))) ) then
         Hx(m)   = Hx(m) + Tij_x
         Hy(m)   = Hy(m) + Tij_y
         Hz(m)   = Hz(m) + Tij_z
         Hx(m+1) = conjg(Hx(m))
         Hy(m+1) = conjg(Hy(m))
         if(ii .eq. jj) then ! check if diagonal part
           Hz(m+1) = 0d0
         else
           Hz(m+1) = conjg(Hz(m))
         endif
         return
       endif
     enddo
     mm = mm + 1
      I(mm) = ii
      J(mm) = jj
      Hx(mm) = Tij_x
      Hy(mm) = Tij_y
      Hz(mm) = Tij_z
      IJ(mm) = real(I(mm)) + real(J(mm)) * zi
     mm = mm + 1
      I(mm) = jj
      J(mm) = ii
      Hx(mm) = conjg(Hx(mm-1))
      Hy(mm) = conjg(Hy(mm-1))
      if(ii .eq. jj) then  !check if diagonal part 
        Hz(mm) = 0d0
      else
        Hz(mm) = conjg(Hz(mm-1))
      endif
      IJ(mm) = real(I(mm)) + real(J(mm)) * zi

     if(mm .gt. n_neighbor*4) then
       if_main write(6,'(A)')'    !WARN! Number of non-zero element NNZ > n_neighbor*4.'
       if_main write(6,'(A)')'           Please check "set_ham_soc_sparse" routine. Exit program...'
       kill_job
     endif
   else

     return
   endif

return
endsubroutine
#endif

subroutine get_hopping_integral(Eij, NN_TABLE, nn, PINPT, tol, kpoint, FIJ_, flag_phase, flag_set_overlap)
  use parameters, only : zi, hopping, incar
  use phase_factor
  interface
    function FIJ_(k,R)
      complex*16   FIJ_
      real*8, intent(in) :: k(3)
      real*8, intent(in) :: R(3)
    endfunction
  end interface
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  integer*4         nn
  real*8            kpoint(3), tol, tij_sk, tij_cc
  complex*16        Eij
  external          tij_sk, tij_cc
  logical           flag_phase, flag_set_overlap

  if(PINPT%flag_slater_koster) then
    if(flag_phase) then
      Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false., flag_set_overlap) * FIJ_( kpoint, NN_TABLE%Rij(1:3,nn))
    elseif(.not. flag_phase) then
      Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false., flag_set_overlap) * FIJ_( kpoint, NN_TABLE%R  (1:3,nn))
    endif
  elseif(.not.PINPT%flag_slater_koster) then
    if(flag_phase) then
      Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ_( kpoint, NN_TABLE%Rij(1:3,nn))
    elseif(.not. flag_phase) then
      Eij = tij_cc(NN_TABLE,nn,PINPT,tol,.false.) * FIJ_( kpoint, NN_TABLE%R  (1:3,nn))
    endif
  endif

return
endsubroutine
