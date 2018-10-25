#include "alias.inc"
#ifdef MKL_SPARSE
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

  ! DEALLOCATION of array
  ! Hk: initialized every call (deallocated with cal_eig_Hk_sparse exit)
  ! H0: initialized every call (deallocated with cal_eig_Hk_sparse exit)
  ! Hm: initialized in the first call (deallocated after get_eig exit)
  ! Hs: initialized in the first call if slater-koster (deallocated after get_eig exit)
  !     initialized every call        if .not. slater-koster (deallocated after cal_eig_Hk_sparse exit)

  ! H0 will be constructed for spin-1. For spin-2, copied from spin-1 multiplied by -1
  if(is .eq. 1) call set_ham0_sparse   (SH0, kp, PINPT, neig, NN_TABLE, flag_phase)

  if(flag_init) then   ! setup k-independent Hamiltonian: Hm, Hs (if .not. slater_koster)
    if(PINPT%flag_collinear) then
      call set_ham_mag_sparse(SHm, NN_TABLE, PINPT, neig, flag_sparse_zero_SHm)
    elseif(PINPT%flag_noncollinear) then
      call set_ham_mag_sparse(SHm, NN_TABLE, PINPT, neig, flag_sparse_zero_SHm)

      if(PINPT%flag_soc .and. PINPT%flag_slater_koster) then
        call set_ham_soc_sparse(SHs, 0d0, PINPT, neig, NN_TABLE, flag_phase, flag_sparse_zero_SHs)
      endif
    endif
    flag_init = .false.
  endif

  if(PINPT%flag_collinear) then
    if(.not. flag_sparse_zero_SHm) then
      call sparse_create_csr_handle(S0, SH0)
      call sparse_create_csr_handle(Sm, SHm)
      alpha = ((-1d0)**(is+1))

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
!       call set_ham_soc_sparse(SHs, kp, PINPT, neig, NN_TABLE, flag_phase)
        if_main write(6,'(A)')'    !WARN! Current version does not support Sparse Matrix '
        if_main write(6,'(A)')'           for non-Slater-Koster type Hamiltonian'
        if_main write(6,'(A)')'           Exit program...'
        stop
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
!write(6,*)"CCCCC"
!stop
!     istat = MKL_SPARSE_DESTROY(Sk_)
!     if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sk_', istat)
!     istat = MKL_SPARSE_DESTROY(Ss)    
!     if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Ss ', istat)
!     istat = MKL_SPARSE_DESTROY(Sm)    
!     if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: Sm ', istat)
!     istat = MKL_SPARSE_DESTROY(S0)    
!     if(istat .gt. 1) call sparse_error_report('MKL_SPARSE_DESTROY: S0 ', istat)
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
    call get_hopping_integral(Eij, NN_TABLE, nn, PINPT, tol, kpoint, F_IJ, flag_phase)

    if(ii .eq. jj) then
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
subroutine set_ham_soc_sparse(SHs, kp , PINPT, neig, NN_TABLE, flag_phase, flag_sparse_zero)
  use parameters, only : zi, hopping, incar, spmat, eta
  use sparse_tool
  use kronecker_prod
  use mpi_setup
#ifdef MKL_SPARSE
  use MKL_SPBLAS
#endif
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  type (spmat  ) :: COO_x, COO_y, COO_z, CSR_x, CSR_y, CSR_z, SHs
  type (SPARSE_MATRIX_T) :: SHx, SHy, SHz, SHxy, SHxyz
  integer*4         neig
  integer*4         nn, ii,jj, mm, m
  integer*4         istat
  integer*4         soc_index 
  real*8            lambda_soc
  real*8            kp(3)
  logical           flag_phase, flag_sparse_zero
  complex*16        Tij, Tij_x, Tij_y, Tij_z
  integer*4         I(NN_TABLE%n_neighbor*2) ! default maximum : n_neighbor*4
  integer*4         J(NN_TABLE%n_neighbor*2)
  complex*16        Hx(NN_TABLE%n_neighbor*2)
  complex*16        Hy(NN_TABLE%n_neighbor*2)
  complex*16        Hz(NN_TABLE%n_neighbor*2)
  complex*16        IJ(NN_TABLE%n_neighbor*2)
  complex*16        alpha
  complex*16        L_x, L_y, L_z
  external          L_x, L_y, L_z

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
          call get_param(PINPT,    soc_index, lambda_soc   )
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
   
      if_main write(6,*)'    !WARN! Current version does not support Sparse Matrix for non-Slater-Koster type Hamiltonian'
      if_main write(6,*)'           Exit program...'
      stop

    endif

return
endsubroutine
#endif

subroutine get_hopping_integral(Eij, NN_TABLE, nn, PINPT, tol, kpoint, FIJ_, flag_phase)
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
  logical           flag_phase

  if(PINPT%flag_slater_koster) then
    if(flag_phase) then
      Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ_( kpoint, NN_TABLE%Rij(1:3,nn))
    elseif(.not. flag_phase) then
      Eij = tij_sk(NN_TABLE,nn,PINPT,tol,.false.) * FIJ_( kpoint, NN_TABLE%R  (1:3,nn))
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
