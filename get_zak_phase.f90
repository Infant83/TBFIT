#include "alias.inc"
subroutine get_zak_phase(NN_TABLE, PINPT, PPRAM, PINPT_BERRY, PGEOM, PKPTS)
   use parameters, only : incar, hopping, poscar, berry, energy, kpoints, pi2, params
   use berry_phase
   use mpi_setup
   use print_io
   implicit none
   type(hopping)          :: NN_TABLE
   type(incar)            :: PINPT
   type(params)           :: PPRAM
   type(berry)            :: PINPT_BERRY
   type(poscar)           :: PGEOM
   type(kpoints)          :: PKPTS
   integer*4                 mpierr
   integer*4                 ikpath, is
   integer*4                 nkdiv, nkpath, nerange
   integer*4                 iband, nband
   integer*4, allocatable :: erange(:)
   real*8                    time1, time2
   real*8                    G1(3), G2(3)
   real*8,    allocatable :: E(:,:)
   complex*16,allocatable :: V(:,:,:)
   complex*16,allocatable :: SV(:,:,:)
   real*8                    zak_phase(PINPT%nspin,PINPT_BERRY%zak_nkpath)
   real*8                    zak_phase_(PINPT%nspin,PINPT_BERRY%zak_nkpath)
   real*8                    kpoint(3,PINPT_BERRY%zak_nkdiv,PINPT_BERRY%zak_nkpath)
   logical                   flag_phase, flag_sparse, flag_order
#ifdef MPI
   if_main time1 = MPI_Wtime()
#else
   call cpu_time(time1)
#endif

   call set_berry_erange(PINPT_BERRY, PGEOM, PINPT, 'zk')
   allocate(erange(PINPT_BERRY%zak_nerange))
   allocate(E(PINPT_BERRY%zak_nerange,PINPT_BERRY%zak_nkdiv))
   allocate(V(PGEOM%neig*PINPT%ispin,PINPT_BERRY%zak_nerange,PINPT_BERRY%zak_nkdiv))
   allocate(SV(PGEOM%neig*PINPT%ispin,PINPT_BERRY%zak_nerange,PINPT_BERRY%zak_nkdiv))

   call set_berry_kpath (PINPT_BERRY, PGEOM, PINPT, 'zk')

   allocate(PINPT_BERRY%zak_phase(PINPT%nspin,PINPT_BERRY%zak_nkpath))
   allocate(PINPT_BERRY%polarization(PINPT%nspin))
   PINPT_BERRY%zak_phase = 0d0
   PINPT_BERRY%polarization = 0d0
   flag_sparse= .false.
   flag_phase = PINPT_BERRY%flag_zak_phase
!  flag_phase= .TRUE.  
   zak_phase = 0d0
   zak_phase_= 0d0
#ifdef MPI
   zak_phase_= 0d0
#endif
   flag_order= .false. ! it should not be .true.

!  NOTE: MPI parallelism need to be improved, and the memory handling as well.

   nkdiv  = PINPT_BERRY%zak_nkdiv
   nkpath = PINPT_BERRY%zak_nkpath
   G1     = PINPT_BERRY%zak_kpoint(:,nkdiv,1) - PINPT_BERRY%zak_kpoint(:,1,1)
   G2     = PINPT_BERRY%zak_kpoint(:,1,nkpath) - PINPT_BERRY%zak_kpoint(:,1,1)
   nerange= PINPT_BERRY%zak_nerange
   nband  = PINPT_BERRY%zak_nerange/PINPT%nspin
   erange = PINPT_BERRY%zak_erange
   kpoint = PINPT_BERRY%zak_kpoint
   iband  = erange(1)

   write(message,*)'' ; write_msg
   write(message,'(A)')'START: ZAK PHASE EVALUATION' ; write_msg
   write(message,'(A,A)')'  BAND INDEX: ',adjustl(trim(PINPT_BERRY%strip_zak_range)) ; write_msg

   do ikpath = 1, nkpath
     call get_eig(NN_TABLE, kpoint(:,:,ikpath), nkdiv, PINPT, PPRAM, E, V, SV, PGEOM%neig, iband, nband, .true.,flag_sparse, .false., flag_phase)
     call set_periodic_gauge(V, G1, PINPT, PGEOM, nkdiv, erange, nerange)
#ifdef F08
     call get_berry_phase(zak_phase(:,ikpath), kpoint(:,:,ikpath), V, PINPT, PGEOM, nkdiv, erange, nerange)
#else
     call get_berry_phase_det(zak_phase(:,ikpath), kpoint(:,:,ikpath), V, PINPT, PGEOM, nkdiv, erange, nerange)
#endif
     write(message,'(A,I0,A,I0)')"  STATUS: ",ikpath,' / ',nkpath ; write_msg
   enddo

#ifdef MPI 
   PINPT_BERRY%zak_phase = zak_phase
   if_main call print_zakphase(PINPT, PINPT_BERRY)
#else
   PINPT_BERRY%zak_phase = zak_phase
   call print_zakphase(PINPT, PINPT_BERRY)
#endif

   do is = 1, PINPT%nspin
     do ikpath = 1, nkpath-1
       PINPT_BERRY%polarization(is) = PINPT_BERRY%polarization(is) + PINPT_BERRY%zak_phase(is,ikpath) / (nkpath-1)
     enddo
     write(message,*)" POLARIZATION:", PINPT_BERRY%polarization(is) * area(G1,G2) ; write_msg
   enddo

#ifdef MPI
   if_main time2 = MPI_Wtime()
#else
   call cpu_time(time2)
#endif
   write(message,'(A,F12.3)')'END: ZAK PHASE EVALUATION. TIME ELAPSED (s) =',time2-time1 ; write_msg

   return
endsubroutine
