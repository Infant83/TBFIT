#include "alias.inc"
subroutine get_wcc(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
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
   integer*4        is
   integer*4        ikpath, nkpath
   integer*4        nkdiv, nerange
   integer*4        erange(PINPT_BERRY%wcc_nerange)
   real*8           G(3)
   real*8           time1, time2
   real*8           E(PGEOM%neig*PINPT%ispin,PINPT_BERRY%wcc_nkdiv)
   complex*16       V(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin,PINPT_BERRY%wcc_nkdiv)
   real*8           wcc(PINPT_BERRY%wcc_nerange/PINPT%nspin,PINPT%nspin,PINPT_BERRY%wcc_nkpath)
   real*8           largest_gap(PINPT%nspin,PINPT_BERRY%wcc_nkpath)
   integer*4        clock_direct (PINPT%nspin,PINPT_BERRY%wcc_nkpath)
   integer*4        z2_index(PINPT%nspin)
#ifdef MPI
   real*8           wcc_(PINPT_BERRY%wcc_nerange/PINPT%nspin,PINPT%nspin,PINPT_BERRY%wcc_nkpath)
#endif
   real*8           kpoint(3,PINPT_BERRY%wcc_nkdiv,PINPT_BERRY%wcc_nkpath)
   logical          flag_phase, flag_get_chern

#ifdef MPI
   if_main time1 = MPI_Wtime()
#else
   call cpu_time(time1)
#endif

   if_main write(6,*)''
   if_main write(6,'(A)')'START: WCC EVALUATION'
   if_main write(6,'(A,A)')'  BAND INDEX: ',adjustl(trim(PINPT_BERRY%strip_wcc_range))

   ! NOTE : The range of WCC will be [0:1] with the unit of lattice vector. 
   !        To get polarization of i direction, you can multiply lattice parameter a_i of i-direction
   !        and electric charge e. See Section II-B of PRB 95, 075146 (2017) for the details.
   allocate(PINPT_BERRY%wcc(PINPT_BERRY%wcc_nerange/PINPT%nspin,PINPT%nspin,PINPT_BERRY%wcc_nkpath))
   if(PINPT_BERRY%flag_wcc_get_chern) then
     allocate(PINPT_BERRY%wcc_chern(PINPT%nspin))
     allocate(PINPT_BERRY%wcc_polarization(PINPT%nspin,PINPT_BERRY%wcc_nkpath))
     PINPT_BERRY%wcc_chern = 0d0
     PINPT_BERRY%wcc_polarization = 0d0
   endif

   PINPT_BERRY%wcc = 0d0 
   flag_phase = .TRUE. 
   flag_get_chern = PINPT_BERRY%flag_wcc_get_chern
   wcc = 0d0
#ifdef MPI
   wcc_= 0d0
#endif

   nkdiv  = PINPT_BERRY%wcc_nkdiv
   nkpath = PINPT_BERRY%wcc_nkpath
   nerange= PINPT_BERRY%wcc_nerange
   erange = PINPT_BERRY%wcc_erange(:)
   kpoint = PINPT_BERRY%wcc_kpoint
   G      = PINPT_BERRY%wcc_kpoint(:,nkdiv,1) - PINPT_BERRY%wcc_kpoint(:,1,1)
   do ikpath = 1,  nkpath
     call get_eig(NN_TABLE, kpoint(:,:,ikpath), nkdiv, PINPT, E, V, PGEOM%neig,.true.,.false., flag_phase)
     call set_periodic_gauge(V, G, PINPT, PGEOM, nkdiv, erange, nerange)
     call get_berry_phase(wcc(:,:,ikpath),kpoint(:,:,ikpath), V, PINPT, PGEOM, nkdiv, erange, nerange)
     if_main write(6,'(A,I,A,I)')"  STATUS: ",ikpath,' / ',nkpath
   enddo

#ifdef MPI
   ! MPI routine is not supported yet...
!  call MPI_Reduce(wcc, wcc_, size(wcc), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
!  PINPT_BERRY%wcc = wcc_

   PINPT_BERRY%wcc = wcc
!  if_main call print_wcc(PINPT, PINPT_BERRY)
   call find_largest_gap(largest_gap, clock_direct, z2_index, PINPT_BERRY%wcc, PINPT%nspin, nkpath, nerange)
   if(flag_get_chern) call get_chern_number(PINPT_BERRY%wcc_chern(:), PINPT_BERRY%wcc_polarization(:,:), &
                                            PINPT_BERRY%wcc(:,:,:), PINPT%nspin, nkpath, nerange)
   if_main call print_wcc(PINPT, PINPT_BERRY%wcc_filenm, PINPT_BERRY%wcc_gap_filenm, &
                          PINPT_BERRY%wcc, PINPT_BERRY%wcc_kpath, nkpath, &
                          nerange, PINPT_BERRY%strip_wcc_range, largest_gap, clock_direct, z2_index, &
                          PINPT_BERRY%wcc_polarization(:,:), PINPT_BERRY%wcc_chern(:), flag_get_chern)

#else
   PINPT_BERRY%wcc = wcc 
!  call print_wcc(PINPT, PINPT_BERRY)
   call find_largest_gap(largest_gap, clock_direct, z2_index, PINPT_BERRY%wcc, PINPT%nspin, nkpath, nerange)
   if(flag_get_chern) call get_chern_number(PINPT_BERRY%wcc_chern(:), PINPT_BERRY%wcc_polarization(:,:), &
                                            PINPT_BERRY%wcc(:,:,:), PINPT%nspin, nkpath, nerange)
   call print_wcc(PINPT, PINPT_BERRY%wcc_filenm, PINPT_BERRY%wcc_gap_filenm, &
                  PINPT_BERRY%wcc, PNPT_BERRY%wcc_kpath, nkpath, &
                  nerange, PINPT_BERRY%strip_wcc_range, largest_gap, clock_direct, z2_index, &
                  PINPT_BERRY%wcc_polarization(:,:), PINPT_BERRY%wcc_chern(:), flag_get_chern)
#endif

#ifdef MPI
   if_main time2 = MPI_Wtime()
#else
   call cpu_time(time2)
#endif
   if_main write(6,'(A,I1)')'  Z2 INDEX     =   ',z2_index(:)
   if_main write(6,'(A,F12.3)')'END: WCC EVALUATION. TIME ELAPSED (s) =',time2-time1

   return
endsubroutine
