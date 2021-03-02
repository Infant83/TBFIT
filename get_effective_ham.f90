#include "alias.inc"
subroutine get_eig_downfold(PINPT, PPRAM, PKPTS, PGEOM, NN_TABLE)
   use mpi_setup
   use parameters, only : incar, poscar, hopping, kpoints, energy, params
   use time
   use do_math
   use print_matrix
   use print_io
   implicit none
   type(hopping) :: NN_TABLE
   type(incar  ) :: PINPT
   type(params ) :: PPRAM
   type(energy ) :: ETBA_FULL
   type(energy)  :: ETBA_EFF
   type(kpoints) :: PKPTS
   type(poscar)  :: PGEOM
   integer*4     neig, neig_eff, nband, nkp, ispinor
   integer*4     is, ik, ii, ie, fe, im, fm, iadd
   integer*4     iband
   integer*4     mpierr
   character*100 stat
   logical       flag_vector, flag_sparse, flag_init
   logical       flag_phase, flag_stat
   real*8        t1, t0
   real*8        kp(3,PKPTS%nkpoint)
   complex*16    Hm(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor) ! collinear magnetism hamiltonian (k-independent)
   complex*16    Hs(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor) ! 1st-order SO coupling hamiltonian (k-dependent if .not. SK)
   complex*16    H0(PGEOM%neig,PGEOM%neig)                             ! slater-koster hopping (k-dependent)
   complex*16    Hk(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor) ! total hamiltonian (k-dependent)
   real*8        E(PGEOM%nband*PINPT%nspin)                      ! will store all the energy eigenvalue for each spin
   integer*4     index_a(PGEOM%neig_eff * PINPT%ispin)
   integer*4     index_b((PGEOM%neig - PGEOM%neig_eff) * PINPT%ispin)
   complex*16    V(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin)     ! will store all the spin block at once in the first dimension 
   complex*16    Ha (PGEOM%neig_eff*PINPT%ispinor,PGEOM%neig_eff*PINPT%ispinor)                                 !  Hk = [ Ha   Hab ]  
   complex*16    Hb ((PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor,(PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor)   !       [ Hab'  Hb ]
   complex*16    Hab(PGEOM%neig_eff *PINPT%ispinor,(PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor)                 !
   complex*16    Hef(PGEOM%neig_eff*PINPT%ispinor,PGEOM%neig_eff*PINPT%ispinor)  ! Hef = Ha + Hab * (e - Hb)^-1 * Hab' 
   real*8        EEF(PGEOM%neig_eff*PINPT%ispin,PKPTS%nkpoint)                                ! will store all the energy eigenvalue for each spin
   complex*16    VEF(PGEOM%neig_eff*PINPT%ispin,PGEOM%neig_eff*PINPT%ispin,PKPTS%nkpoint)     ! will store all the spin block at once in the first dimension 
   character*80  fname_header
#ifdef MPI
   integer*4     ourjob(nprocs)
   integer*4     ourjob_disp(0:nprocs-1)
   call mpi_job_distribution_chain(PKPTS%nkpoint, nprocs, ourjob, ourjob_disp) 
#else
   integer*4     ourjob(1)
   integer*4     ourjob_disp(0)
   call mpi_job_distribution_chain(PKPTS%nkpoint, nprocs, ourjob, ourjob_disp)
#endif

!  NOTE: This subroutine has not been MPI parallized in the current version, 
!        since constructing effective hamiltonian 
!        does not take too much resources. However, in the future release, we will
!        address this subroutine and make it viable along with MPI.

   flag_vector = .true.
   flag_sparse = .false.
   flag_phase  = .true.
   flag_stat   = .true.
   write(fname_header,*)'band_structure_TBA_EFF',trim(PINPT%title(NN_TABLE%mysystem))
!  fname_header='band_structure_TBA_EFF'

   ispinor = PINPT%ispinor
   neig    = PGEOM%neig
   neig_eff= PGEOM%neig_eff
   nband   = PGEOM%nband
   iband   = PGEOM%init_erange
   kp      = PKPTS%kpoint   
   nkp     = PKPTS%nkpoint
   index_a = NN_TABLE%i_eff_orb ; call get_index_b(index_a, index_b, PGEOM, PINPT)
   EEF     = 0d0
   VEF     = (0d0,0d0)
   allocate(ETBA_EFF%E(PGEOM%neig_eff*PINPT%ispin,PKPTS%nkpoint))
   allocate(ETBA_EFF%V(PGEOM%neig_eff*PINPT%ispin,PGEOM%neig_eff*PINPT%ispin, PKPTS%nkpoint))

   if(iband .ne. 1 .and. nband .ne. PGEOM%neig*PINPT%ispinor) then
     write(message,'(A)')'    !WARN! In the case of "SET EFFECTIVE ", ERANGE tag cannot not be applied.'  ; write_msg
!    kill_job
   endif

   write(message,'(A)')' '  ; write_msg
   write(message,'(A)')'    * START: GET EFFECTIVE HAMILTONIAN AND ITS EIGENVALUE '  ; write_msg
   write(message,'(A)')' '  ; write_msg

   flag_init = .true.
   call time_check(t1,t0,'init')
   call initialize_eig_status(ii, iadd, stat, nkp)
   E = 0d0 ; V = (0.d0,0.d0) 
   
 k_loop:do ik= sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
     do is = 1, PINPT%nspin
       call get_hamk_dense(Hk, H0, Hm, Hs, is, kp(:,ik), PINPT, PPRAM, neig, NN_TABLE, flag_init, flag_phase)
       call subdivide_H(Hk, Ha, Hb, Hab, is, PINPT, PGEOM, index_a, index_b)

       call get_effective_ham(Hef, Ha, Hb, Hab, PINPT, PGEOM)
       call get_matrix_index(ie, fe, im, fm, is, neig_eff*ispinor, neig_eff, PINPT%ispinor)
       call cal_eig_hermitian(Hef,neig_eff*ispinor, EEF(ie:fe,ik), flag_vector) ; if(flag_vector) VEF(im:fm,ie:fe,ik) = Hef
     enddo
   enddo k_loop

#ifdef MPI
   call MPI_ALLREDUCE(EEF, ETBA_EFF%E, size(ETBA_EFF%E), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   if(flag_vector) call MPI_ALLREDUCE(VEF, ETBA_EFF%V, size(ETBA_EFF%V), MPI_COMPLEX16, MPI_SUM, mpi_comm_earth, mpierr)
#else
   ETBA_EFF%E = EEF
   ETBA_EFF%V = VEF
#endif

   if_main call print_energy_eff(PKPTS, ETBA_EFF%E, ETBA_EFF%V, PGEOM, PINPT, neig_eff, fname_header)

   write(message,'(A)')' '  ; write_msg
   write(message,'(A)')'    * END: GET EFFECTIVE HAMILTONIAN AND ITS EIGENVALUE '  ; write_msg
   write(message,'(A)')' '  ; write_msg

   return
endsubroutine

subroutine get_effective_ham(Hef, Ha, Hb, Hab, PINPT, PGEOM)
   use parameters, only : incar, poscar
   use do_math
   use print_matrix
   implicit none
   type(incar)   ::   PINPT
   type(poscar)  ::   PGEOM
   real*8        e
   complex*16    Ha (PGEOM%neig_eff*PINPT%ispinor,PGEOM%neig_eff*PINPT%ispinor)
   complex*16    Hb ((PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor,(PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor)
   complex*16    Hab(PGEOM%neig_eff*PINPT%ispinor,(PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor)
   complex*16    Hef(PGEOM%neig_eff*PINPT%ispinor,PGEOM%neig_eff*PINPT%ispinor)  

   e = (PINPT%eff_emax + PINPT%eff_emin)/2d0
   ! NOTE: In the current version only the energy dependent downfolded Hamiltonian will be constructed.
   !       In this case, we will restrict the energy e will be set by center of EFF_EWINDOW for the 
   !       convenience. In the near future, energy-independent version will be implemented. HJK.
   !       The detailed descriptiob and procedures can be found this the following paper:
   !       E. Zurek, O. Jepsen, and O. K. Anderson, "Muffin-Tin Orbital Wannier-Like Functions for 
   !       Insulators and Metals", Chem. Phys. Chem. 6, 1934-1942 (2005)
   !       or
   !       R. Norman-Elvenich, "Effective Low-Energy Hamiltonians for Correlated Transition Metal
   !       Compounds", Diploma Thesis, Technische Universitat Graz, (2010)

   ! Hef = Ha + Hab * (e - Hb)^-1 * Hab'
   Hef = Ha + matmul(matmul(Hab, invc((e - Hb))), conjg(transpose(Hab)))

   return
endsubroutine

subroutine subdivide_H(Hk, Ha, Hb, Hab, is, PINPT, PGEOM, index_a, index_b)
   use parameters, only : incar, poscar
   use print_matrix
   implicit none
   type(incar)   ::   PINPT
   type(poscar)  ::   PGEOM
   integer*4     is, na, nb, ia, fa, ib, fb
   integer*4     index_a(PGEOM%neig_eff * PINPT%ispin)
   integer*4     index_b((PGEOM%neig - PGEOM%neig_eff) * PINPT%ispin)
   integer*4     index_a_(PGEOM%neig_eff * PINPT%ispinor)
   integer*4     index_b_((PGEOM%neig-PGEOM%neig_eff) * PINPT%ispinor)
   complex*16    Hk(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16    Ha (PGEOM%neig_eff*PINPT%ispinor,PGEOM%neig_eff*PINPT%ispinor)
   complex*16    Hb ((PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor,(PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor)
   complex*16    Hab(PGEOM%neig_eff *PINPT%ispinor,(PGEOM%neig - PGEOM%neig_eff)*PINPT%ispinor)

   na = PGEOM%neig_eff * PINPT%ispinor
   nb =(PGEOM%neig-PGEOM%neig_eff)* PINPT%ispinor
   ia= 1+PGEOM%neig_eff * (is - 1)
   fa= PGEOM%neig_eff * (is - 1) * PINPT%ispinor
   ib= 1+(PGEOM%neig-PGEOM%neig_eff) * (is - 1)
   fb= (PGEOM%neig-PGEOM%neig_eff) * (is - 1) * PINPT%ispinor
   index_a_ = index_a ( ia:fa ) - PGEOM%neig * (is - 1)
   index_b_ = index_b ( ib:fb ) - PGEOM%neig * (is - 1)

   Ha = Hk(index_a_, index_a_)
   Hb = Hk(index_b_, index_b_)
   Hab= Hk(index_a_, index_b_)

   return
endsubroutine

subroutine get_index_b(index_a, index_b, PGEOM, PINPT)
   use parameters, only : poscar, incar
   implicit none
   type(poscar)  ::   PGEOM
   type(incar )  ::   PINPT
   integer*4          i, ia, ib
   integer*4          index_a(PGEOM%neig_eff * PINPT%ispin)
   integer*4          index_b((PGEOM%neig - PGEOM%neig_eff) * PINPT%ispin)

   ia = 1
   ib = 1
lp:do i = 1, PGEOM%neig * PINPT%ispin
     if(index_a(ia) .eq. i ) then
       ia = ia + 1
       cycle lp
     elseif(index_a(ia) .ne. i) then
       index_b(ib) = i
       ib = ib + 1
       cycle lp
     endif
   enddo lp

   return
endsubroutine

subroutine set_effective_orbital_index(PINPT,PGEOM, NN_TABLE)
   use mpi_setup
   use parameters, only : incar, poscar, hopping
   use sorting
   implicit none
   type(incar)   ::   PINPT
   type(poscar)  ::   PGEOM
   type(hopping) ::   NN_TABLE
   character*40       c_atom_orb
   integer*4          i, iorb, ieff, imatrix
   integer*4          nitems, neff_set, neig_eff
   external           nitems
   character*40, allocatable :: c_eff_orb(:)
   integer*4,    allocatable :: i_eff_orb(:)
   integer*4                    i_eff_orb_(PGEOM%neig)

   neff_set = nitems(PINPT%eff_orb_dummyc)
   allocate(c_eff_orb(neff_set))
   read(PINPT%eff_orb_dummyc, *) c_eff_orb(1:neff_set)
   i_eff_orb_ = 0

   neig_eff = 0
   do ieff = 1, neff_set
     do i = 1, PGEOM%n_atom
       do iorb = 1, PGEOM%n_orbital(i)
         write(c_atom_orb,'(A,A,A)')trim(PGEOM%c_spec(PGEOM%spec(i))), ':', &
                                    trim(PGEOM%c_orbital(iorb,i))
         imatrix= sum( PGEOM%n_orbital(1:i) ) - PGEOM%n_orbital(i) + iorb

         if(trim(c_atom_orb) .eq. trim(c_eff_orb(ieff))) then
           neig_eff = neig_eff + 1
           i_eff_orb_(neig_eff) = imatrix
         endif
       enddo
     enddo
   enddo

   allocate(i_eff_orb(neig_eff))
   i_eff_orb(1:neig_eff) = i_eff_orb_(1:neig_eff)
   call get_sort_variable_1D_int(i_eff_orb, neig_eff, 'inc')

   if(PINPT%ispin .eq. 2) then
     allocate(NN_TABLE%i_eff_orb( neig_eff * 2 ))
     NN_TABLE%i_eff_orb(1:neig_eff) = i_eff_orb
     NN_TABLE%i_eff_orb(1+neig_eff:neig_eff*2) = i_eff_orb + PGEOM%neig
   else
     allocate(NN_TABLE%i_eff_orb( neig_eff ))
     NN_TABLE%i_eff_orb(1:neig_eff) = i_eff_orb
   endif

   PGEOM%neig_eff = neig_eff

   deallocate(i_eff_orb)
   deallocate(c_eff_orb)
   return
endsubroutine
