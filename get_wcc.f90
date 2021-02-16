#include "alias.inc"
subroutine get_wcc(NN_TABLE, PINPT, PPRAM, PINPT_BERRY, PGEOM, PKPTS)
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
   integer*4                 is
   integer*4                 ikpath, nkpath
   integer*4                 nkdiv, nerange, nband, iband
   integer*4, allocatable :: erange(:)
   real*8                    G(3)
   real*8                    time1, time2
   real*8,    allocatable :: E(:,:)
   complex*16,allocatable :: V(:,:,:)
   complex*16,allocatable :: SV(:,:,:)
   real*8,    allocatable :: wcc(:,:,:)
#ifdef MPI                   
   real*8,    allocatable :: wcc_(:,:,:)
#endif                       
   real*8                    largest_gap(PINPT%nspin,PINPT_BERRY%wcc_nkpath)
   integer*4                 clock_direct (PINPT%nspin,PINPT_BERRY%wcc_nkpath)
   integer*4                 z2_index(PINPT%nspin)
   real*8                    kpoint(3,PINPT_BERRY%wcc_nkdiv,PINPT_BERRY%wcc_nkpath)
   logical                   flag_phase
   logical                   flag_sparse, flag_get_chern, flag_get_chern_spin, flag_order


#ifdef MPI
   if_main time1 = MPI_Wtime()
#else
   call cpu_time(time1)
#endif

   call set_berry_erange(PINPT_BERRY, PGEOM, PINPT, 'wc')
   allocate(erange(PINPT_BERRY%wcc_nerange))
   allocate(E(PINPT_BERRY%wcc_nerange,PINPT_BERRY%wcc_nkdiv))
   allocate(V(PGEOM%neig*PINPT%ispin,PINPT_BERRY%wcc_nerange,PINPT_BERRY%wcc_nkdiv))
   allocate(SV(PGEOM%neig*PINPT%ispin,PINPT_BERRY%wcc_nerange,PINPT_BERRY%wcc_nkdiv))
   allocate(wcc(PINPT_BERRY%wcc_nerange/PINPT%nspin,PINPT%nspin,PINPT_BERRY%wcc_nkpath))
#ifdef MPI
   allocate(wcc_(PINPT_BERRY%wcc_nerange/PINPT%nspin,PINPT%nspin,PINPT_BERRY%wcc_nkpath))
#endif
   call set_berry_kpath (PINPT_BERRY, PGEOM, PINPT, 'wc')

   write(message,*)'' ; write_msg
   write(message,'(A)')      ' #--- START: WCC EVALUATION' ; write_msg
   write(message,'(A,A)')'  BAND INDEX: ',adjustl(trim(PINPT_BERRY%strip_wcc_range)) ; write_msg
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
   flag_sparse= .false. ! current version does not support sparse matrix for wcc evaluation
   flag_phase = PINPT_BERRY%flag_wcc_phase
   flag_get_chern      = PINPT_BERRY%flag_wcc_get_chern
   flag_get_chern_spin = PINPT_BERRY%flag_wcc_get_chern_spin
   wcc = 0d0
#ifdef MPI
   wcc_= 0d0
#endif
   flag_order = .false. ! it should not be .true.

   nkdiv  = PINPT_BERRY%wcc_nkdiv
   nkpath = PINPT_BERRY%wcc_nkpath
   nerange= PINPT_BERRY%wcc_nerange
   nband  = PINPT_BERRY%wcc_nerange/PINPT%nspin
   erange = PINPT_BERRY%wcc_erange(:)
   kpoint = PINPT_BERRY%wcc_kpoint
   G      = PINPT_BERRY%wcc_kpoint(:,nkdiv,1) - PINPT_BERRY%wcc_kpoint(:,1,1)
   iband  = erange(1)

   do ikpath = 1,  nkpath
     call get_eig(NN_TABLE, kpoint(:,:,ikpath), nkdiv, PINPT, PPRAM, E, V, SV, PGEOM%neig, iband, nband,.true., flag_sparse, .false., flag_phase)
     if_main call set_periodic_gauge(V, G, PINPT, PGEOM, nkdiv, erange, nerange)
#ifdef F08
     if_main call get_berry_phase(wcc(:,:,ikpath),kpoint(:,:,ikpath), V, PINPT, PGEOM, nkdiv, erange, nerange)
#else
     if_main call get_berry_phase_svd(wcc(:,:,ikpath),kpoint(:,:,ikpath), V, PINPT, PGEOM, nkdiv, erange, nerange)
#endif
     write(message,'(A,I0,A,I0)')"  STATUS: ",ikpath,' / ',nkpath ; write_msg
   enddo

#ifdef MPI
   ! NOTE: MPI routine is not supported yet... 
   ! However, get_eig is MPI activated, which gives fairly good performance.
   ! Probabliy this routine need to be updated to use full cpu loads in "get_berry_phase" routines, 
   ! but still fast enough since determining bulk topology usually calculated within primitive unitcell 
   ! where hamiltonian size is quite small.

   PINPT_BERRY%wcc = wcc
   call find_largest_gap(largest_gap, clock_direct, z2_index, PINPT_BERRY%wcc, PINPT%nspin, nkpath, nerange)
   if(flag_get_chern) call get_chern_number(PINPT_BERRY%wcc_chern(:), PINPT_BERRY%wcc_polarization(:,:), &
                                            PINPT_BERRY%wcc(:,:,:), PINPT%nspin, nkpath, nerange)
   if_main call print_wcc(PINPT, PINPT_BERRY%wcc_filenm, PINPT_BERRY%wcc_gap_filenm, &
                          PINPT_BERRY%wcc, PINPT_BERRY%wcc_kpath, nkpath, &
                          nerange, PINPT_BERRY%strip_wcc_range, largest_gap, clock_direct, z2_index, &
                          PINPT_BERRY%wcc_polarization(:,:), PINPT_BERRY%wcc_chern(:), flag_get_chern)

#else
   PINPT_BERRY%wcc = wcc 
   call find_largest_gap(largest_gap, clock_direct, z2_index, PINPT_BERRY%wcc, PINPT%nspin, nkpath, nerange)
   if(flag_get_chern) call get_chern_number(PINPT_BERRY%wcc_chern(:), PINPT_BERRY%wcc_polarization(:,:), &
                                            PINPT_BERRY%wcc(:,:,:), PINPT%nspin, nkpath, nerange)
   call print_wcc(PINPT, PINPT_BERRY%wcc_filenm, PINPT_BERRY%wcc_gap_filenm, &
                  PINPT_BERRY%wcc, PINPT_BERRY%wcc_kpath, nkpath, &
                  nerange, PINPT_BERRY%strip_wcc_range, largest_gap, clock_direct, z2_index, &
                  PINPT_BERRY%wcc_polarization(:,:), PINPT_BERRY%wcc_chern(:), flag_get_chern)
#endif

#ifdef MPI
   if_main time2 = MPI_Wtime()
#else
   call cpu_time(time2)
#endif
   if_main_then
      if(.not. flag_get_chern) then 
        write(message,'(A,I1)')'  Z2 INDEX     =   ',z2_index(:) ; write_msg
      endif
      write(message,'(A,F12.3)')' #--- END: WCC EVALUATION. TIME ELAPSED (s) =',time2-time1 ; write_msg
      write(message,'(A      )')' '; write_msg
   if_main_end
   return
endsubroutine
