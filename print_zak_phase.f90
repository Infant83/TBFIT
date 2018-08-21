#include "alias.inc"
subroutine print_zakphase(PINPT, PINPT_BERRY)
   use parameters, only : incar, berry, pid_zak
   implicit none
   type(incar)   :: PINPT
   type(berry)   :: PINPT_berry
   integer*4    ikpath
   integer*4    is
   real*8    kk

   open(pid_zak, file=trim(PINPT_BERRY%zak_filenm), status='unknown')
   write(pid_zak,'(A,A)')"# BAND INDEX: ",adjustl(trim(PINPT_BERRY%strip_zak_range))

   do ikpath = 1, PINPT_BERRY%zak_nkpath
     write(pid_zak,'(A,I5,A,3F10.5,A,3F10.5)',advance='yes')'# K-PATH',ikpath,': ', &
                                                PINPT_BERRY%zak_kpath(1:3,1,ikpath),' ->', &
                                                PINPT_BERRY%zak_kpath(1:3,2,ikpath)
   enddo

   if(PINPT%nspin .eq. 2) then
     write(pid_zak,'(A)')'# k-path,  zak_phase (spin-up, spin-dn)'
   else
     write(pid_zak,'(A)')'# k-path,  zak_phase'
   endif

   do ikpath = 1, PINPT_BERRY%zak_nkpath
     if(PINPT%nspin .eq. 2) then
       kk= 1d0 / (PINPT_BERRY%zak_nkpath-1) * (ikpath - 1)
       write(pid_zak,'(F10.6,2F16.8)')kk, (PINPT_BERRY%zak_phase(is,ikpath),is=1,PINPT%nspin)
     else
       kk= 1d0 / (PINPT_BERRY%zak_nkpath-1) * (ikpath - 1)
       write(pid_zak,'(F10.6, F16.8)')kk,  PINPT_BERRY%zak_phase(1,ikpath)
     endif
   enddo

   close(pid_zak)

   return
endsubroutine
