#include "alias.inc"
! this routine calculates circular dichroism 
! see Ref: W. Yao, D. Xiao, and Q. Niu, PRB 77, 235406 (2008) 
subroutine get_circular_dichroism(NN_TABLE, PINPT, PGEOM, PKPTS, ETBA)
   use parameters, only : incar, hopping, poscar, energy, kpoints, pauli_0, hbar, eta
   use berry_phase
   use phase_factor
   use mpi_setup
   use time
   use memory
   use kronecker_prod, only : kproduct
   use print_io
   implicit none
   type(hopping) :: NN_TABLE
   type(incar)   :: PINPT
   type(poscar)  :: PGEOM
   type(energy)  :: ETBA
   type(kpoints) :: PKPTS
   real*8           kpoint(3,PKPTS%nkpoint)
   real*8           time1, time2
   character*4      timer
   complex*16,allocatable :: myV(:,:,:)
   integer*4        ii, ik, is, ie
   integer*4        nkpoint, neig, nband, ispin, nspin, ispinor, msize
   integer*4        ieig, feig
   integer*4        sizebuff
   integer*4        istat, iadd
   integer*4        ourjob(nprocs), ourjob_disp(0:nprocs-1)
   integer*4        id, my_ik, mpierr
   complex*16       dxH(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16       dyH(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   complex*16       dzH(PGEOM%neig*PINPT%ispinor,PGEOM%neig*PINPT%ispinor)
   real*8           circ_dichroism(PINPT%nspin, PKPTS%nkpoint, PINPT%ncirc_dichroism)
   real*8           circ_dichroism_(PINPT%nspin, PKPTS%nkpoint, PINPT%ncirc_dichroism)
   logical          flag_phase
   character*100    stat
   character*2      spin(2)
   ! NOTE1: this routine does not work properly with sparse matrix and erange tag
   !        i.e., if you provide ERANGE and EWINDOW tag, there should be error...
   !        We will add some warning messages later on. 

   ! NOTE2: this routine does not support non-orthogonal basis case, where "USE_OVERLAP .TRUE."
   !        in the future release we hope to put thoes functionality. 
   !        For the non-orthgonal basis case, one can refer to following paper:
   !        Lee et al., PRB 98, 115115 (2018): "Tight-binding calculations of optical matrix element
   !                                            for conductivity using nonorthogonal atomic orbitals:
   !                                            Anomalous Hall conductivity in bcc Fe"

   call time_check(time2,time1,timer)
   write(message,*)''  ; write_msg
   write(message,'(A)')' ---- START: CIRCULAR DICHROISM -----------'   ; write_msg

   kpoint         = PKPTS%kpoint   ; nkpoint         = PKPTS%nkpoint
   neig           = PGEOM%neig     ; ispinor         = PINPT%ispinor
   nband          = PGEOM%nband    ; nspin           = PINPT%nspin
   msize          = neig * ispinor ; ispin           = PINPT%ispin
   circ_dichroism = 0d0            ; circ_dichroism_ = 0d0
   flag_phase     = .TRUE.         ; spin(:)         = (/'up','dn'/)! default
   sizebuff       = neig*ispin*nband*nspin
   call mpi_job_distribution_chain(nkpoint, ourjob, ourjob_disp)
   call report_job_distribution(.true., ourjob)
   if_main call report_memory( int8(sizebuff*nkpoint*2), 16, 'Eigen vectors')

   if(.not. allocated(myV)) allocate(myV(neig*ispinor,nband*nspin,ourjob(myid+1)))
#ifdef MPI
   sizebuff = neig*ispin*nband*nspin 
   call MPI_SCATTERV(ETBA%V, ourjob*sizebuff, &
                             ourjob_disp*sizebuff, &
                             MPI_COMPLEX16, myV, ourjob(myid+1)*sizebuff, &
                             MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
#else
   myV=ETBA%V
#endif

   do ii = 1, PINPT%ncirc_dichroism
     do is = 1, PINPT%nspin
       ieig = PINPT%circ_dichroism_pair(1,ii) + (is-1)*nband
       feig = PINPT%circ_dichroism_pair(2,ii) + (is-1)*nband
       write(message,*)' '  ; write_msg
       if(nspin .eq. 2) then
         write(message,'(3A,I0,A,I0)')'     -- BAND INDEX (',trim(spin(is)),') : from ', ieig, ' to ', feig ; write_msg
       else
         write(message, '(A,I0,A,I0)')'     -- BAND INDEX : from ', ieig, ' to ', feig ; write_msg
       endif
       
       call initialize_eig_status(istat, iadd, nkpoint)
       write(message,'(A)') '        STATUS: ' ; call write_log(trim(message), 1, myid) ; call write_log(trim(message),22, myid)

       do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
         my_ik = ik - sum(ourjob(1:myid)) 
         call get_velocity_matrix(PINPT, NN_TABLE, kpoint(:,ik), neig, dzH, dzF_IJ, flag_phase)
         call get_velocity_matrix(PINPT, NN_TABLE, kpoint(:,ik), neig, dxH, dxF_IJ, flag_phase)
         call get_velocity_matrix(PINPT, NN_TABLE, kpoint(:,ik), neig, dyH, dyF_IJ, flag_phase)
         call get_optical_transition(circ_dichroism(is,ik,ii), ETBA%E((/ieig,feig/),ik), &
                                                                myV(:,(/ieig,feig/),my_ik), dxH, dyH, dzH, msize, 1)
       
         if_main call print_eig_status(ik, istat, iadd, ourjob)
      
       enddo ! ik

     enddo ! is
   enddo ! ii

#ifdef MPI
   call MPI_ALLREDUCE(circ_dichroism, circ_dichroism_, size(circ_dichroism), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   circ_dichroism = circ_dichroism_
#endif

   if_main  write(message,'(A)')' ' ; write_msg
   if_main  call print_circ_dichroism(circ_dichroism, ETBA%E, PINPT, PKPTS, PGEOM)

   call time_check(time2, time1)
   write(message,*)' ' ; write_msg
   write(message,'(A,F12.3)')'   TIME for CIRCULAR DICHROISM CALCULATION (s) ', time2 ; write_msg
   write(message,'(A      )')' ---- END: CIRCULAR DICHROISM -----------'; write_msg

   return
endsubroutine

subroutine get_optical_transition(optical_selectivity, E, V, dxH, dyH, dzH, msize, imode)
   use parameters, only : eta, zi
   implicit none
   integer*4   n, m
   integer*4   im
   integer*4   msize, imode
   real*8      E(2), de2
   complex*16  V(msize,2)
   complex*16  dxH(msize, msize)
   complex*16  dyH(msize, msize)
   complex*16  dzH(msize, msize)
   complex*16  transition_matrix_left
   complex*16  transition_matrix_right
   real*8      optical_selectivity
   complex*16  vx_nm, vy_nm, vz_nm
   complex*16  vx_mn, vy_mn, vz_mn
   complex*16  psi_n(msize), psi_m(msize)

!  transition_matrix = 0d0

   if(imode .eq. 1) then
     psi_n = V(:,1)
     psi_m = V(:,2)
!    do im = 1, msize
!      psi_m = 
!    enddo
    !write(6,*)"ZZZZ ", psi_n
    !vy_nm = dot_product( psi_n, matmul(dyH,psi_m) )
    !vx_nm = dot_product( psi_n, matmul(dxH,psi_m) )

     vx_mn = dot_product( psi_m, matmul(dxH,psi_n) )
     vy_mn = dot_product( psi_m, matmul(dyH,psi_n) )
  
     transition_matrix_left = vx_mn + zi * vy_mn
     transition_matrix_right= vx_mn - zi * vy_mn
     optical_selectivity    =   ((abs(transition_matrix_left))**2 - (abs(transition_matrix_right))**2) & 
                              / ((abs(transition_matrix_left))**2 + (abs(transition_matrix_right))**2)
  
   endif

   return
endsubroutine
