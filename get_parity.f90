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
   complex*16        parity_eigenvalue(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        parity_matrix_op(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        Hk(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        Hk_(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16        V(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   real*8            E(PGEOM%neig*PINPT%ispinor), E_(PGEOM%neig*PINPT%ispinor)
   logical           flag_phase
   integer*4         neig, ispin, nspin, ispinor
   character*2       spin(2)
   real*8            kpoint(3,PINPT_BERRY%parity_nkpoint)
   real*8            kpoint_reci(3,PINPT_BERRY%parity_nkpoint)
   character*128     title
   character*10      kpoint_name(PINPT_BERRY%parity_nkpoint)
   complex*16        phi_i(PGEOM%neig*PINPT%ispinor), phi_j(PGEOM%neig*PINPT%ispinor)
   real*8            origin(3)

   neig       = PGEOM%neig
   ispin      = PINPT%ispin
   nspin      = PINPT%nspin
   ispinor    = PINPT%ispinor
   nkp        = PINPT_BERRY%parity_nkpoint
   kpoint     = PINPT_BERRY%parity_kpoint(:,:)
   kpoint_reci= PINPT_BERRY%parity_kpoint_reci(:,:)
   kpoint_name= PINPT_BERRY%parity_kpoint_name
   origin     = PINPT_BERRY%parity_origin
   spin       = (/'up','dn'/)
   flag_phase = .false. 

   if_main write(6,*)' '
   if_main write(6,'(A)')'START: PARITY EIGENVALUE CALCULATION'

   call set_parity_matrix_op(PGEOM, PINPT, parity_matrix_op, PINPT_BERRY%parity_operator, origin)

sp:do is = 1, nspin
     if_main_then
       if(nspin .eq. 2) write(6,'(3A)')'     SPIN',spin(is),':'
     if_main_end

  kp:do ik = 1, nkp
       if_main write(6,'(A)')' '
       if_main write(6,'(A,I3,1x,A10,A,3F10.5)')'       KPT',ik,kpoint_name(ik),' : ',kpoint_reci(:,ik)
       parity_eigenvalue = 0d0

       call get_hamk(Hk, NN_TABLE, PINPT, kpoint(:,ik), is, neig, flag_phase);Hk_ = Hk
       call cal_eig_hermitian(Hk_, neig*ispinor, E, .true.)
       call symmetrize_hamk(Hk, V, parity_matrix_op, neig, ispinor)
       call cal_eig_hermitian(V, neig*ispinor, E_, .true.)
       call get_symmetry_eigenvalue(parity_eigenvalue, V, parity_matrix_op, neig, ispinor)

       if_main_then
!        if(nspin .eq. 2) write(title,'(5A)')'               <i,',trim(spin(is)),'|P|j,',trim(spin(is)),'>'
!        if(nspin .eq. 1) write(title,'( A)')'               <i|P|j>'
!        call print_matrix_r(real(parity_eigenvalue),neig*ispinor, neig*ispinor, trim(title),0,'x,F4.0')
         if(nspin .eq. 2) write(6,'(5A)')'                   E(n,k)      <n,',trim(spin(is)),'|P|n,',trim(spin(is)),'>'
         if(nspin .eq. 1) write(6,'( A)')'                   E(n,k)      <n|P|n>'
         do i = 1, neig*ispinor
           write(6,'(A,I5,A,F10.5,3x,*(F5.2,"+",F5.2,"i"))')'      E(n=',i,')', E(i), parity_eigenvalue(i,i)
         enddo
       if_main_end

     enddo kp 
  
   enddo sp
 
   if_main write(6,*)' '
   if_main write(6,'(A)')'  END: PARITY EIGENVALUE CALCULATION'

   return
endsubroutine
