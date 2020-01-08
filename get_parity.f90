#include "alias.inc"
subroutine get_parity(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
   use parameters, only : hopping, incar, berry, poscar, kpoints, pi2, pauli_0, pauli_x, pauli_y, pauli_z
   use berry_phase
   use mpi_setup
   use phase_factor
   use kronecker_prod, only: kproduct
   use print_matrix
   implicit none
   type (incar)   :: PINPT       ! parameters for input arguments
   type (berry)   :: PINPT_BERRY ! parameters for berry phase related quantity calculation
   type (poscar)  :: PGEOM       ! parameters for geometry info
   type (kpoints) :: PKPTS       ! parameters for kpoints
   type (hopping) :: NN_TABLE    ! table for hopping index
   integer*4         nkp
   integer*4         i,j, is, ik
   integer*4         ie, je
   integer*4         mpierr
   complex*16        P_OP(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        Hk(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        Hk_(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        V(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        Sij(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor) ! symmetry matrix S_ij = 
   real*8            E(PGEOM%neig*PINPT%ispinor), E_(PGEOM%neig*PINPT%ispinor)
   complex*16        S_eig(PGEOM%neig*PINPT%ispinor)
   logical           flag_phase, flag_phase_shift, flag_print_ham
   integer*4         neig, ispin, nspin, ispinor
   character*2       spin(2)
   real*8            kpoint(3,PINPT_BERRY%parity_nkpoint)
   real*8            kpoint_reci(3,PINPT_BERRY%parity_nkpoint)
   character*128     title
   character*10      kpoint_name(PINPT_BERRY%parity_nkpoint)
   complex*16        phi_i(PGEOM%neig*PINPT%ispinor), phi_j(PGEOM%neig*PINPT%ispinor)
   complex*16        phase_shift(PGEOM%neig*PINPT%ispinor)
   real*8            origin(3)
   real*8            R(3,3)
   real*8            Rcart(3,3) ! rotation operator for crystal structure (cartesian coord)

   neig       = PGEOM%neig
   ispin      = PINPT%ispin
   nspin      = PINPT%nspin
   ispinor    = PINPT%ispinor
   nkp        = PINPT_BERRY%parity_nkpoint
   kpoint     = PINPT_BERRY%parity_kpoint(:,:)
   kpoint_reci= PINPT_BERRY%parity_kpoint_reci(:,:)
   kpoint_name= PINPT_BERRY%parity_kpoint_name
   origin     = PINPT_BERRY%parity_origin
   R          = PINPT_BERRY%parity_operator
  !Rcart(1,:)  = (/ cos(theta/180d0*pi), sin(theta/180d0*pi), 0d0/)
  !Rcart(2,:)  = (/-sin(theta/180d0*pi), cos(theta/180d0*pi), 0d0/)
  !Rcart(3,:)  = (/                 0d0,                 0d0, 1d0/)
   spin       = (/'up','dn'/)
   flag_phase = PINPT_BERRY%flag_parity_phase
   flag_phase_shift = .false.
   flag_print_ham = PINPT_BERRY%flag_print_hamiltonian_parity

   if_main write(6,*)' '
   if_main write(6,'(A)')'START: PARITY EIGENVALUE CALCULATION'

   call set_parity_matrix_op(PGEOM, PINPT, P_OP, R, origin)
   !if_main call print_matrix_i(int(P_OP), neig*ispinor, neig*ispinor, 'P_OP', 0, 'I3')
   if(flag_print_ham .and. myid .eq. 0) call print_matrix_i(int(P_OP), neig*ispinor, neig*ispinor, 'P_OP', 0, 'I3')
   
sp:do is = 1, nspin
      if(nspin .eq. 2 .and. myid .eq. 0) write(6,'(3A)')'     SPIN',spin(is),':'

  kp:do ik = 1, nkp
       if_main write(6,'(A)')' '
       if_main write(6,'(A,I3,1x,A10,A,3F10.5)')'       KPT',ik,kpoint_name(ik),' : ',kpoint_reci(:,ik)

       call get_ham_Hk(Hk, NN_TABLE, PINPT, kpoint(:,ik), is, neig, flag_phase) ; V = Hk
       !call symmetrize_hamk(Hk, V, P_OP, neig, ispinor) ! return V

       if(flag_print_ham .and. myid .eq. 0) call print_matrix_c(V, neig*ispinor, neig*ispinor, 'H_'//trim(PINPT_BERRY%parity_kpoint_name(ik)), 1, 'F6.3')

       call cal_eig_hermitian(V, neig*ispinor, E, .true.)
       !call get_phase_shift(phase_shift,kpoint(:,ik),PGEOM,ispinor)
       phase_shift = 1d0
       call get_symmetry_matrix(Sij, S_eig, V, phase_shift, P_OP, E, PGEOM, neig, ispinor, flag_phase_shift)
       if(flag_print_ham .and. myid .eq. 0) call print_matrix_c(Sij, neig*ispinor, neig*ispinor, 'Sij ', 0, 'F6.3')

       if_main_then
         if(nspin .eq. 2) write(6,'(5A)')'                   E(n,k)      <n,',trim(spin(is)),'|P|n,',trim(spin(is)),'> (parity)'
         if(nspin .eq. 1) write(6,'( A)')'                   E(n,k)      <n|P|n> (parity)'
         do i = 1, neig*ispinor
           !write(6,'(A,I5,A,F10.5,*(F11.2,A,F4.1,A))')'      E(n=',i,')', E(i) , real(Sij(i,i)),'+',aimag(Sij(i,i)),'i'
           write(6,'(A,I5,A,F10.5,*(3x,F7.4,A,F7.4,A))')'      E(n=',i,')', E(i) , real(S_eig(i)),' + ',aimag(S_eig(i)),'i'
           if(i .eq. PINPT_BERRY%noccupied) then
             write(6,'(A)')'     ------------- occupied level --------------------'
           endif
         enddo
       if_main_end

     enddo kp 
  
   enddo sp
 
   if_main write(6,*)' '
   if_main write(6,'(A)')'  END: PARITY EIGENVALUE CALCULATION'

   return
endsubroutine
