subroutine get_eig(NN_TABLE, kpoint, nkpoint, PINPT, E, V, neig, flag_get_orbital, flag_count_percent, flag_phase)
  use parameters, only: hopping, incar, pauli_0, pauli_x, pauli_y, pauli_z
  use kronecker_prod, only: kproduct
  use mpi_setup
  use phase_factor
  use do_math
  implicit none
  type (hopping) :: NN_TABLE
  type (incar  ) :: PINPT
  integer*4  iadd, ii
  integer*4  ik,neig
  integer*4  nkpoint
  integer*4  mpierr
  integer*4  is, ispin
  integer*4  irow,frow, icol,fcol
  real*8     percent
  real*8     kpoint(3,nkpoint)
  real*8     E(neig*PINPT%ispin,nkpoint)
  complex*16 V(neig*PINPT%ispin,neig*PINPT%ispin,nkpoint) 
#ifdef MPI
  real*8     E_(neig*PINPT%ispin,nkpoint)
  complex*16 V_(neig*PINPT%ispin,neig*PINPT%ispin,nkpoint) 
#endif
  complex*16 H0(neig,neig)                             ! slater-koster hopping (k-dependent)
  complex*16 Hk(neig*PINPT%ispinor,neig*PINPT%ispinor) ! total hamiltonian (k-dependent)
  complex*16 Hm(neig*PINPT%ispinor,neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
  complex*16 Hs(neig*PINPT%ispinor,neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
  logical, intent(in) :: flag_get_orbital, flag_count_percent
  character*100 cpercent
  real*8   time_0, time_1
  logical  flag_phase

#ifdef MPI
  time_0=MPI_Wtime()
#else
  call cpu_time(time_0)
#endif

  E = 0d0
  if(flag_get_orbital) V = (0.d0,0.d0)
#ifdef MPI
  E_= 0d0
  if(flag_get_orbital) V_= (0.d0,0.d0)
#endif

  if(flag_count_percent .and. myid .eq. 0) call initialize_eig_status(ii, iadd, cpercent, nkpoint)

  if(PINPT%flag_collinear) then
    call set_ham_mag(Hm, NN_TABLE, PINPT, neig)
  elseif(PINPT%flag_noncollinear) then
    call set_ham_mag(Hm, NN_TABLE, PINPT, neig)
    if(PINPT%flag_soc .and. PINPT%flag_slater_koster) call set_ham_soc(Hs, 0d0, PINPT, neig, NN_TABLE, F_IJ, flag_phase)
  endif


 k_loop:do ik= 1 + myid, nkpoint, nprocs
    call set_ham0(H0, kpoint(:,ik), PINPT, neig, NN_TABLE, F_IJ, flag_phase)
    if(PINPT%flag_collinear) then
      do is = 1, PINPT%nspin
        Hk = H0 + ((-1d0)**(is+1)) * Hm
        irow = 1 + (is-1)*neig ; frow = neig + (is-1)*neig
        icol = 1 + (is-1)*neig ; fcol = neig + (is-1)*neig
#ifdef MPI
        call cal_eig_hermitian(Hk, neig, E_(irow:frow,ik), flag_get_orbital)
        if(flag_get_orbital) V_(irow:frow,icol:fcol,ik) = Hk(1:neig,1:neig)
#else
        call cal_eig_hermitian(Hk, neig, E(irow:frow,ik), flag_get_orbital)
        if(flag_get_orbital) V(irow:frow,icol:fcol,ik) = Hk(1:neig,1:neig)
#endif
      enddo
 
    elseif(PINPT%flag_noncollinear) then
      if(PINPT%flag_soc) then
        !set up k-dependent SOC in the case of 'cc' orbitals
        if(.not. PINPT%flag_slater_koster) call set_ham_soc(Hs, kpoint(:,ik), PINPT, neig, NN_TABLE, F_IJ, flag_phase) 
        Hk = kproduct(pauli_0, H0, 2, 2, neig, neig) + Hm + Hs
      else
        Hk = kproduct(pauli_0, H0, 2, 2, neig, neig) + Hm 
      endif


#ifdef MPI
      call cal_eig_hermitian(Hk, neig*2, E_(1:neig*2,ik), flag_get_orbital)
      if(flag_get_orbital) V_(1:neig*2,1:neig*2,ik) = Hk(1:neig*2,1:neig*2)
#else
      call cal_eig_hermitian(Hk, neig*2, E(1:neig*2,ik), flag_get_orbital)
      if(flag_get_orbital) V(1:neig*2,1:neig*2,ik) = Hk(1:neig*2,1:neig*2)
#endif
 
    elseif(.not. PINPT%flag_collinear .and. .not. PINPT%flag_noncollinear) then
#ifdef MPI
      call cal_eig_hermitian(H0, neig, E_(1:neig,ik), flag_get_orbital)
      if(flag_get_orbital) V_(1:neig,1:neig,ik) = H0(1:neig,1:neig)
#else
      call cal_eig_hermitian(H0, neig, E(1:neig,ik), flag_get_orbital)
      if(flag_get_orbital) V(1:neig,1:neig,ik) = H0(1:neig,1:neig)
#endif
    endif
 
    if(flag_count_percent .and. myid .eq. 0) call print_eig_status(ik, ii, iadd, cpercent, nkpoint)
  enddo k_loop

#ifdef MPI
  call MPI_ALLREDUCE(E_, E, size(E), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
  if(flag_get_orbital) call MPI_ALLREDUCE(V_, V, size(V), MPI_COMPLEX16, MPI_SUM, mpi_comm_earth, mpierr)
  time_1=MPI_Wtime()
  if(flag_count_percent .and. myid .eq. 0) write(6,'(A,F12.6)')"TIME for EIGENVALUE SOLVE (s)", time_1 - time_0
#else
  call cpu_time(time_1)
  if(flag_count_percent .and. myid .eq. 0) write(6,'(A,F12.6)')"TIME for EIGENVALUE SOLVE (s)", time_1 - time_0
#endif


return
endsubroutine
subroutine set_ham0(H, kpoint, PINPT, neig , NN_TABLE, FIJ, flag_phase)
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
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    elseif(.not.PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    endif

    if(i .eq. j) then
      if(nint(PINPT%param_const(4,NN_TABLE%local_U_param_index(i))) .ge. 1) then
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param_const(5,(NN_TABLE%local_U_param_index(i)))
      else
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param((NN_TABLE%local_U_param_index(i)))
      endif
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
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    elseif(.not.PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    endif

    if(i .eq. j) then
      if(nint(PINPT%param_const(4,NN_TABLE%local_U_param_index(i))) .ge. 1) then
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param_const(5,(NN_TABLE%local_U_param_index(i)))
      else
        H(i,j) = H(i,j) + Eij + NN_TABLE%local_charge(i)*PINPT%param((NN_TABLE%local_U_param_index(i)))
      endif
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
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_sk(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
      endif
    elseif(.not.PINPT%flag_slater_koster) then
      if(flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%Rij(1:3,nn))
      elseif(.not. flag_phase) then
        Eij = tij_cc(NN_TABLE,nn,PINPT,tol) * FIJ(kpoint, NN_TABLE%R  (1:3,nn))
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
            H(i,i+neig) = H(i,i) + zi * 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                              * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i))
          else
            H(i,i+neig) = H(i,i) + zi * 0.5d0 * NN_TABLE%local_moment_rot(2,i) &
                                              * PINPT%param(NN_TABLE%stoner_I_param_index(i))
          endif
        endif
      enddo

      do i = 1, neig  ! Hy
        if(NN_TABLE%stoner_I_param_index(i) .gt. 0) then   ! if stoner parameter has been set...
          if(nint(PINPT%param_const(4,NN_TABLE%stoner_I_param_index(i))) .eq. 1) then ! if i-th basis has constraint .true.
            H(i,i+neig) = H(i,i+neig)  - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                               * PINPT%param_const(5,NN_TABLE%stoner_I_param_index(i))
          else
            H(i,i+neig) = H(i,i+neig)  - 0.5d0 * NN_TABLE%local_moment_rot(1,i) &
                                               * PINPT%param(NN_TABLE%stoner_I_param_index(i))
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

subroutine allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA)
   use parameters
   implicit none
   type (incar)   :: PINPT       ! parameters for input arguments
   type (energy)  :: ETBA        ! target energy to be fitted to
   type (poscar)  :: PGEOM       ! parameters for geometry info
   type (kpoints) :: PKPTS       ! parameters for kpoints


   allocate(ETBA%E(PGEOM%neig*PINPT%ispin, PKPTS%nkpoint))
   allocate(ETBA%V(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin, PKPTS%nkpoint))

return
endsubroutine
subroutine initialize_eig_status(ii, iadd, cpercent, nkpoint)
   use mpi_setup
   implicit none
   character*100 cpercent
   integer*4     iadd, ii, nkpoint
    
   cpercent = '****************************************************************************************************'
   if(myid .eq. 0) write(6,'(A)')cpercent
   if(nkpoint .le. 25) then
     iadd=floor(real(1)/real(nkpoint)*100d0)
     ii = 1
   else
     iadd=4
     ii = 1
   endif

return
endsubroutine

subroutine print_eig_status(ik, ii, iadd, cpercent, nkpoint)
   use mpi_setup
   implicit none
   character*100   cpercent
   integer*4       ik, ii, iadd, nkpoint
   real*8          percent

   percent =  ik / real(nkpoint) * 100d0
   if( floor(percent) .ge. iadd*ii ) then
     write(6,'(A,I3)')cpercent(1:iadd*ii),floor(percent)
     ii = ii + 1
   endif

return
endsubroutine

