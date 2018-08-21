#include "alias.inc"
subroutine print_wcc(PINPT, filenm, filenm_gap, wcc, kpath, nkpath, nerange_tot, strip_erange, largest_gap, clock_direct, z2_index, polarization, chern, flag_get_chern)
   use parameters, only : incar, pid_wcc
   implicit none
   type(incar)   :: PINPT
   integer*4     pid_gap 
   integer*4     nerange_tot
   integer*4     nkpath, nerange, nspin
   real*8        wcc(nerange_tot/PINPT%nspin,PINPT%nspin,nkpath)
   real*8        kpath(3,2,nkpath)
   real*8        largest_gap(PINPT%nspin,nkpath)
   real*8        polarization(PINPT%nspin, nkpath)
   real*8        chern(PINPT%nspin)
   integer*4     clock_direct (PINPT%nspin,nkpath)
   integer*4     z2_index(PINPT%nspin)
   character*256 strip_erange
   character*40  filenm, filenm_gap
   logical       flag_get_chern

   pid_gap = pid_wcc + 1
   nspin   = PINPT%nspin
   nerange = nerange_tot / nspin
   call write_wcc_header(pid_wcc, pid_gap, filenm, filenm_gap, strip_erange, nspin, kpath, nkpath, z2_index, clock_direct, flag_get_chern, chern) 
   call write_wcc_main(pid_wcc, pid_gap, nerange, nspin, nkpath, wcc, largest_gap, polarization, flag_get_chern)
   return
endsubroutine
subroutine write_wcc_main(pid_wcc_,pid_gap_, nerange, nspin, nkpath, wcc, largest_gap, polarization, flag_get_chern)
   implicit none
   integer*4     pid_gap_, pid_wcc_
   integer*4     ie
   integer*4     is, ikpath
   integer*4     nerange, nspin, nkpath
   real*8        largest_gap(nspin,nkpath)
   real*8        wcc(nerange,nspin,nkpath)
   real*8        polarization(nspin,nkpath)
   real*8        kk
   logical       flag_get_chern


   do ie = 1, nerange
     if(nspin .eq. 2) then
       write(pid_wcc_,'(A,I,A)')'# k-path, ',ie,'-th wcc (up, dn)'
       if(flag_get_chern) then
         if(ie .eq. 1) write(pid_gap_,'(A    )')'# k-path, largest gap of wannier charge center (up, dn), polarization evolution (up, dn)'
       else
         if(ie .eq. 1) write(pid_gap_,'(A    )')'# k-path, largest gap of wannier charge center (up, dn)'
       endif
     else
       write(pid_wcc_,'(A,I,A)')'# k-path, ',ie,'-th wcc'
       if(flag_get_chern) then
         if(ie .eq. 1) write(pid_gap_,'(A    )')'# k-path, largest gap of wannier charge center, polarization evolution'
       else
         if(ie .eq. 1) write(pid_gap_,'(A    )')'# k-path, largest gap of wannier charge center'
       endif
     endif

     do ikpath = 1, nkpath
       if(nspin .eq. 2) then
         kk= 1d0 / (nkpath-1) * (ikpath - 1)
         write(pid_wcc_,'(F10.6,2F16.8)')kk, (wcc(ie,is,ikpath),is=1,nspin)
         if(flag_get_chern) then
           if(ie .eq. 1) write(pid_gap_,'(F10.6,2F16.8,2x,2F16.8)')kk, (largest_gap(is,ikpath),is=1,nspin), (polarization(is,ikpath),is=1,nspin)
         else
           if(ie .eq. 1) write(pid_gap_,'(F10.6,2F16.8)')kk, (largest_gap(is,ikpath),is=1,nspin)
         endif
       else
         kk= 1d0 / (nkpath-1) * (ikpath - 1)
         write(pid_wcc_,'(F10.6, F16.8)')kk,  wcc(ie,1,ikpath)
         if(flag_get_chern) then
           if(ie .eq. 1) write(pid_gap_,'(F10.6, F16.8, 2x, F16.8)')kk, largest_gap(1,ikpath), polarization(1,ikpath)
         else
           if(ie .eq. 1) write(pid_gap_,'(F10.6, F16.8)')kk, largest_gap(1,ikpath)
         endif
       endif
     enddo

     write(pid_wcc_,*)' '
     write(pid_wcc_,*)' '

   enddo

   close(pid_wcc_)
   close(pid_gap_)

   return
endsubroutine
subroutine write_wcc_header(pid_wcc_, pid_gap_, wcc_filenm, wcc_gap_filenm, strip_wcc_range, nspin, kpath, nkpath, z2_index, clock_direct, flag_get_chern, chern)
   implicit none
   integer*4     pid_gap_, pid_wcc_
   integer*4     is, ikpath
   integer*4     nspin, nkpath
   real*8        largest_gap(nspin,nkpath)
   real*8        kpath(3,2,nkpath)
   integer*4     clock_direct (nspin,nkpath)
   integer*4     z2_index(nspin)
   real*8        chern(nspin)
   logical       flag_get_chern
   character*2   c_spin(2)
   character*40  wcc_filenm, wcc_gap_filenm
   character*256 strip_wcc_range

   c_spin(1) = 'up'; c_spin(2) = 'dn'

   open(pid_wcc_, file=trim(wcc_filenm), status='unknown')
   open(pid_gap_, file=trim(wcc_gap_filenm), status='unknown')
   write(pid_wcc_,'(A,A)')'# BAND INDEX: ',adjustl(trim(strip_wcc_range))
   write(pid_wcc_,'(A)')'# NOTE: (-1)^Z2_INDEX = PI_ik Gap_jumps(ik) (ik/nkpath < 0.5)'
   write(pid_wcc_,'(A)')'#       "Gap_jumps -1 -> odd  number of gap jumps,'
   write(pid_wcc_,'(A)')'#       "Gap_jumps  1 -> even (or zero) number of gap jumps'
   if(nspin .eq. 2) then
     write(pid_wcc_,'(3A,I3)')'# TOPOLOGICAL INDEX:    Z2 INDEX for spin-',c_spin(1 ),' = ', z2_index(1)
     write(pid_wcc_,'(3A,I3)')'#                                for spin-',c_spin(2 ),' = ', z2_index(2)
     if(flag_get_chern) then
       write(pid_wcc_,'(3A,I3)')'#                  : CHERN NUMBR for spin-',c_spin(1 ),' = ', nint(chern(1))
       write(pid_wcc_,'(3A,I3)')'#                                for spin-',c_spin(2 ),' = ', nint(chern(2))
     endif
   elseif(nspin .eq. 1) then
     write(pid_wcc_,'( A,I3)')'# TOPOLOGICAL INDEX:    Z2 INDEX = ', z2_index(1)
     if(flag_get_chern) then
       write(pid_wcc_,'( A,I3)')'#                  : CHERN NUMBR = ', nint(chern(1))
     endif
   endif

   do ikpath = 1, nkpath
     if(nspin .eq. 2) then
       write(pid_wcc_,'(A,I5,A,3F10.5,A,3F10.5,A,2I3)',advance='yes')'# K-PATH',ikpath,': ', &
                                                  kpath(1:3,1,ikpath),' ->', &
                                                  kpath(1:3,2,ikpath),' Gap_jumps(up,dn):',clock_direct(:,ikpath)
     elseif(nspin .eq. 1) then
       write(pid_wcc_,'(A,I5,A,3F10.5,A,3F10.5,A, I3)',advance='yes')'# K-PATH',ikpath,': ', &
                                                  kpath(1:3,1,ikpath),' ->', &
                                                  kpath(1:3,2,ikpath),' Gap_jumps:',clock_direct(1,ikpath)

     endif
   enddo

   return
endsubroutine
