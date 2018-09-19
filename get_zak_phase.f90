#include "alias.inc"
subroutine get_zak_phase(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
   use parameters, only : incar, hopping, poscar, berry, energy, kpoints, pi2
   use berry_phase
   use mpi_setup
   implicit none
   type(hopping) :: NN_TABLE
   type(incar)   :: PINPT
   type(berry)   :: PINPT_BERRY
   type(poscar)  :: PGEOM
   type(kpoints) :: PKPTS
   integer*4        mpierr
   integer*4        ikpath, is
   integer*4        nkdiv, nkpath, nerange
   integer*4        erange(PINPT_BERRY%zak_nerange)
   real*8           time1, time2
   real*8           G1(3), G2(3)
   real*8           E(PGEOM%neig*PINPT%ispin,PINPT_BERRY%zak_nkdiv)
   complex*16       V(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin,PINPT_BERRY%zak_nkdiv)
   real*8           zak_phase(PINPT%nspin,PINPT_BERRY%zak_nkpath)
   real*8           zak_phase_(PINPT%nspin,PINPT_BERRY%zak_nkpath)
   real*8           kpoint(3,PINPT_BERRY%zak_nkdiv,PINPT_BERRY%zak_nkpath)
   logical          flag_phase
#ifdef MPI
   if_main time1 = MPI_Wtime()
#else
   call cpu_time(time1)
#endif

   allocate(PINPT_BERRY%zak_phase(PINPT%nspin,PINPT_BERRY%zak_nkpath))
   allocate(PINPT_BERRY%polarization(PINPT%nspin))
   PINPT_BERRY%zak_phase = 0d0
   PINPT_BERRY%polarization = 0d0
   flag_phase = PINPT_BERRY%flag_zak_phase
!  flag_phase= .TRUE.  
   zak_phase = 0d0
   zak_phase_= 0d0
#ifdef MPI
   zak_phase_= 0d0
#endif

   nkdiv  = PINPT_BERRY%zak_nkdiv
   nkpath = PINPT_BERRY%zak_nkpath
   G1     = PINPT_BERRY%zak_kpoint(:,nkdiv,1) - PINPT_BERRY%zak_kpoint(:,1,1)
   G2     = PINPT_BERRY%zak_kpoint(:,1,nkpath) - PINPT_BERRY%zak_kpoint(:,1,1)
   nerange= PINPT_BERRY%zak_nerange
   erange = PINPT_BERRY%zak_erange
   kpoint = PINPT_BERRY%zak_kpoint

   if_main write(6,*)''
   if_main write(6,'(A)')'START: ZAK PHASE EVALUATION'
   if_main write(6,'(A,A)')'  BAND INDEX: ',adjustl(trim(PINPT_BERRY%strip_zak_range))

   do ikpath = 1, nkpath
     call get_eig(NN_TABLE, kpoint(:,:,ikpath), nkdiv, PINPT, E, V, PGEOM%neig,.true.,.false., flag_phase)
     call set_periodic_gauge(V, G1, PINPT, PGEOM, nkdiv, erange, nerange)
#ifdef F08
     call get_berry_phase(zak_phase(:,ikpath), kpoint(:,:,ikpath), V, PINPT, PGEOM, nkdiv, erange, nerange)
#else
     call get_berry_phase_det(zak_phase(:,ikpath), kpoint(:,:,ikpath), V, PINPT, PGEOM, nkdiv, erange, nerange)
#endif
     if_main write(6,'(A,I,A,I)')"  STATUS: ",ikpath,' / ',nkpath
   enddo

#ifdef MPI 
   ! MPI routine is not supported yet...
!  call MPI_Reduce(zak_phase, zak_phase_, size(zak_phase), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
!  PINPT_BERRY%zak_phase = zak_phase_

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
     write(6,*)" POLARIZATION:", PINPT_BERRY%polarization(is) * area(G1,G2)
   enddo

#ifdef MPI
   if_main time2 = MPI_Wtime()
#else
   call cpu_time(time2)
#endif
   if_main write(6,'(A,F12.3)')'END: ZAK PHASE EVALUATION. TIME ELAPSED (s) =',time2-time1

   return
endsubroutine
