#include "alias.inc"
subroutine set_penalty_orb(NN_TABLE, PINPT, PGEOM, PKPTS, PWGHT, strip_kp, strip_tb, strip_orb, strip_site, strip_pen)
   use parameters
   use mpi_setup
   implicit none
   type(incar )  :: PINPT
   type(poscar)  :: PGEOM
   type(kpoints) :: PKPTS
   type(weight)  :: PWGHT
   type(hopping) :: NN_TABLE
   character(*), parameter :: func = 'set_penalty_orb'
   integer*4        size_krange, size_tbrange, size_orbrange, size_siterange
   integer*4        nrange
   integer*4        krange(PKPTS%nkpoint)
   integer*4        tbrange(PGEOM%neig*PINPT%ispin)
   integer*4        orbrange(PGEOM%max_orb)
   integer*4        siterange(PGEOM%neig*PINPT%ispin)
   real*8           pen_orb
   character(*)     strip_kp, strip_tb, strip_pen, strip_orb, strip_site
   logical          flag_number
   external         flag_number
   integer*4        mpierr

   if_main_then
   write(6,'(10A)')' O_PENALTY: KRANGE=[',trim(strip_kp),'],TB_RANGE=[',trim(strip_tb), &
                   '],ORB_RANGE=[',trim(strip_orb),'],SITE_RANGE=[',trim(strip_site), &
                   '],PENALTY=',trim(strip_pen)
   if_main_end
   !krange set
   call get_krange(strip_kp, PKPTS, krange, size_krange)

   !tbrange set
   call get_tbrange(strip_tb, PGEOM, PINPT, tbrange, size_tbrange)

   !orbrange set
   call get_orbrange(strip_orb, PGEOM, orbrange, size_orbrange)

   !siterange set
   call get_siterange(strip_site, NN_TABLE, PGEOM, PINPT, orbrange(1:size_orbrange), size_orbrange, siterange, size_siterange)

   !weight set
   if(flag_number(trim(strip_pen))) then
     call str2real(trim(strip_pen),pen_orb)
   elseif(.not. flag_number(trim(strip_pen))) then
     if_main write(6,*)'  !WARN! PENALTY setting argument is not properly asigned. '
     if_main write(6,*)'         Check! Exit...',func
     kill_job
   endif

   PWGHT%PENALTY_ORB(siterange(1:size_siterange),tbrange(1:size_tbrange),krange(1:size_krange))=pen_orb

return
endsubroutine
subroutine set_weight(PINPT, PGEOM, PKPTS, PWGHT, EDFT, EDFT_all, strip_kp, strip_tb, strip_df, strip_wt)
   use parameters
   use mpi_setup
   implicit none
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4        nrange
   integer*4        init,fina
   integer*4, allocatable :: irange(:), irange_new(:), krange(:), tbrange(:), dfrange(:)
   type(incar )  ::  PINPT
   type(poscar)  ::  PGEOM
   type(kpoints) ::  PKPTS
   type(weight)  ::  PWGHT
   type(energy)  ::  EDFT
   type(energy)  ::  EDFT_all
   real*8           w
   character(*)     strip_kp, strip_tb, strip_df, strip_wt
   character*20,allocatable :: str(:)
   character*20     dummy
   character(*), parameter :: func = 'set_weight'
   logical          flag_number, flag_collinear, flag_noncollinear
   external         nitems   
   external         flag_number
   integer*4        mpierr

   flag_collinear = PINPT%flag_collinear
   flag_noncollinear = PINPT%flag_noncollinear

   if_main_then
   write(6,'(8A)')' C_WEIGHT: KRANGE=[',trim(strip_kp),'],TB_RANGE=[',trim(strip_tb), &
                     '],DF_RANGE=[',trim(strip_df),'],WEIGHT=',trim(strip_wt)
   if_main_end

   !krange set
   nrange=nitems(strip_kp)
   allocate(str(nrange))
   read(strip_kp,*)str(1:nrange)

   do i=1,nrange

     if(index(str(i),':') .eq. 0) then

       if(flag_number(str(i))) then
         call str2int( trim(str(i)), init )
         fina=init
       elseif(.not. flag_number(str(i))) then
         if(trim(str(i)) .eq. 'NKP' .or. trim(str(i)) .eq. 'KEND') then
           init=PKPTS%nkpoint
           fina=init
         elseif(trim(str(i)) .eq. 'ALL') then
           init=1
           fina=PKPTS%nkpoint
         else
           if_main write(6,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     elseif(index(str(i),':') .gt. 0)then

       pos=index(str(i),':')
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then
         dummy='ALL'
       else
         dummy=dummy(1:pos-1)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NKP' .or. trim(dummy) .eq. 'KEND') then
           init=PKPTS%nkpoint
         elseif(trim(dummy) .eq. 'ALL') then
           init=1
         else
           if_main write(6,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), fina )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NKP' .or. trim(dummy) .eq. 'KEND' .or. trim(dummy) .eq. 'ALL') then 
           fina=PKPTS%nkpoint
         else
           if_main write(6,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     endif

     size_new= size( (/(k, k=init,fina)/) )
     if(.not. allocated(irange)) then
       allocate(irange(size_new))
       irange(:)=( (/(k, k=init,fina)/) )
     elseif(  allocated(irange)) then
       size_now= size(irange)
       allocate(irange_new(size_now+size_new))
!      irange_new=(/irange,init:fina/)
       irange_new=(/irange,(/(k, k=init,fina)/)/)
       deallocate(irange)
       allocate(irange(size_now+size_new))
       irange=irange_new
       deallocate(irange_new)
     endif

   enddo

   allocate(krange(size(irange)))
   krange = irange
   deallocate(str)
   if(allocated(irange))deallocate(irange)
   if(allocated(irange_new))deallocate(irange_new)

   !tbrange set
   nrange=nitems(strip_tb)
   allocate(str(nrange))
   read(strip_tb,*)str(1:nrange)

   do i=1,nrange

     if(index(str(i),':') .eq. 0) then
       if(flag_number(str(i))) then
         call str2int( trim(str(i)), init )
         fina=init
       elseif(.not. flag_number(str(i))) then
         if(trim(str(i)) .eq. 'NEIG' .or. trim(str(i)) .eq. 'EEND') then
           if(flag_noncollinear)then
             init=PGEOM%neig * 2
           else
             init=PGEOM%neig
           endif
           fina=init
         elseif(trim(str(i)) .eq. 'ALL') then
           init=1
           if(flag_noncollinear) then
             fina=PGEOM%neig * 2
           else
             fina=PGEOM%neig
           endif
         else
           if_main write(6,*)' 1!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     elseif(index(str(i),':') .gt. 0)then

       pos=index(str(i),':')
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(1:pos-1)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NEIG' .or. trim(dummy) .eq. 'EEND') then 
           if(flag_noncollinear) then
             init=PGEOM%neig * 2
           else
             init=PGEOM%neig
           endif
         elseif(trim(dummy) .eq. 'ALL')then
           init=1
         else
           if_main write(6,*)' 2!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), fina )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NEIG' .or. trim(dummy) .eq. 'EEND' .or. trim(dummy) .eq. 'ALL') then 
           if(flag_noncollinear) then
             fina=PGEOM%neig * 2
           else
             fina=PGEOM%neig
           endif
         else
           if_main write(6,*)' 3!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     endif

     size_new= size( (/(k, k=init,fina)/) )
     if(.not. allocated(irange)) then
       allocate(irange(size_new))
       irange(:)=( (/(k, k=init,fina)/) )
     elseif(  allocated(irange)) then
       size_now= size(irange)
       allocate(irange_new(size_now+size_new))
       irange_new=(/irange, (/(k, k=init,fina)/) /)
       deallocate(irange)
       allocate(irange(size_now+size_new))
       irange=irange_new
       deallocate(irange_new)
     endif

   enddo

   if(flag_collinear) then
     allocate( tbrange(size(irange)*2) )
     tbrange = (/irange, irange+PGEOM%neig/)
   else
     allocate( tbrange(size(irange)) )
     tbrange = irange
   endif

   deallocate(str)
   if(allocated(irange))    deallocate(irange)
   if(allocated(irange_new))deallocate(irange_new)
 
   !dfrange set
   nrange=nitems(strip_df)
   allocate(str(nrange))
   read(strip_df,*)str(1:nrange)

   do i=1,nrange

     if(index(str(i),':') .eq. 0) then

       if(flag_number(str(i))) then
         call str2int( trim(str(i)), init )
         fina=init
       elseif(.not. flag_number(str(i))) then
         if(trim(str(i)) .eq. 'IBAND') then
           init=PWGHT%iband
           fina=PWGHT%iband
         elseif(trim(str(i)) .eq. 'FBAND') then
           init=PWGHT%fband
           fina=PWGHT%fband
         elseif(trim(str(i)) .eq. 'ALL') then
           init=1
           if(flag_collinear) then
             fina=PGEOM%neig_target/2
           else
             fina=PGEOM%neig_target
           endif
         else
           if_main write(6,*)' 1!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     elseif(index(str(i),':') .gt. 0)then

       pos=index(str(i),':')
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(1:pos-1)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'IBAND') then
           init=PWGHT%iband
         elseif(trim(dummy) .eq. 'FBAND') then
           init=PWGHT%fband
         elseif(trim(dummy) .eq. 'ALL') then
           init=1
         else
           if_main write(6,*)' 2!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), fina )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'IBAND') then
           fina=PWGHT%iband
         elseif(trim(dummy) .eq. 'FBAND') then
           fina=PWGHT%fband
         elseif(trim(dummy) .eq. 'ALL') then
           if(flag_collinear) then
             fina=PGEOM%neig_target/2
           else
             fina=PGEOM%neig_target
           endif
         else 
           if_main write(6,*)' 3!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     endif

     size_new= size( (/(k, k=init,fina)/) )
     if(.not. allocated(irange)) then
       allocate(irange(size_new))
       irange(:)=( (/(k, k=init,fina)/) )
     elseif(  allocated(irange)) then
       size_now= size(irange)
       allocate(irange_new(size_now+size_new))
       irange_new=(/irange, (/(k, k=init,fina)/) /)
       deallocate(irange)
       allocate(irange(size_now+size_new))
       irange=irange_new
       deallocate(irange_new)
     endif

   enddo

   if(flag_collinear) then
     allocate( dfrange(size(irange)*2) )
     dfrange = (/irange, irange + (PGEOM%neig_target/2)/)
   else
     allocate( dfrange(size(irange)) )
     dfrange = irange
   endif

   deallocate(str)
   if(allocated(irange))    deallocate(irange)
   if(allocated(irange_new))deallocate(irange_new)
 
   !weight set
   call str2real(trim(strip_wt),w)

   !reset target energy & weight
   if(size(tbrange) .ne. size(dfrange)) then
     if_main write(6,'(A,A)')'  !WARN! size of TBRANGE and DFRANGE is not same. Please check "SET WEIGHT" tag. Exit...',func
     kill_job
   endif
   EDFT%E(tbrange,krange) = EDFT_all%E(dfrange,krange)
   PWGHT%WT(tbrange,krange) = w

return
endsubroutine

subroutine get_siterange(strip_site, NN_TABLE, PGEOM, PINPT, orbrange, size_orbrange, siterange_dummy, size_siterange)
   use parameters
   use mpi_setup
   implicit none
   type(poscar)  :: PGEOM
   type(incar )  :: PINPT
   type(hopping) :: NN_TABLE
   integer*4        i, j, k, nitems, pos
   integer*4        size_irange
   integer*4, intent(in)  :: size_orbrange
   integer*4, intent(out) :: size_siterange
   integer*4, intent(out) :: siterange_dummy(PGEOM%neig*PINPT%ispin)
   integer*4        range_dummy(PGEOM%n_atom), size_range_dummy
   integer*4        nrange
   integer*4        init, fina
   integer*4, allocatable :: irange(:), irange_new(:), siterange(:)
   integer*4, intent(in) :: orbrange(size_orbrange)
   character(*)     strip_site
   character*20,allocatable :: str(:)
   character*20     dummy, dummy_
   character(*), parameter :: func = 'get_siterange'
   logical          flag_number, flag_collinear, flag_noncollinear
   external         nitems
   external         flag_number
   integer*4        mpierr

   siterange_dummy = 0
   flag_collinear = PINPT%flag_collinear
   flag_noncollinear = PINPT%flag_noncollinear
 
   !siterange set
   nrange=nitems(strip_site)
   allocate(str(nrange))
   read(strip_site,*)str(1:nrange)
 
   do i=1,nrange

     if(index(str(i),':') .eq. 0) then
       dummy = str(i)
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
         if(.not. allocated(irange)) then
           allocate(irange(1))
           irange(1) = init
         elseif(allocated(irange)) then  
           allocate(irange_new(size(irange) + 1))
           irange_new = (/irange, init/)
           deallocate(irange)
           allocate(irange(size(irange_new)))
           irange = irange_new
           deallocate(irange_new)
         endif

       elseif(.not. flag_number(dummy)) then
         call get_site_number(NN_TABLE, PGEOM, dummy, range_dummy, size_range_dummy)
         if(.not. allocated(irange)) then
           allocate(irange(size_range_dummy))
           irange(1:size_range_dummy) = range_dummy(1:size_range_dummy)
         elseif(allocated(irange)) then
           allocate(irange_new(size(irange) + size_range_dummy))
           irange_new = (/irange,range_dummy(1:size_range_dummy)/)
           deallocate(irange)
           allocate(irange(size(irange_new)))
           irange = irange_new
           deallocate(irange_new)
         endif

       endif

     elseif(index(str(i),':') .gt. 0)then
       pos=index(str(i),':')
       dummy = str(i)
       if(len_trim(dummy) .eq. 1) then
         dummy = 'ALL'
       else
         dummy = dummy(1:pos-1)
         dummy_= str(i)
         dummy_= dummy_(pos+1:)
       endif

       if(trim(dummy) .eq. 'ALL') then
         if(allocated(irange)) then
           if_main write(6,*)'  !WARN! SITERANGE setting argument is not properly asigned. '
           if_main write(6,*)'         Check! Exit...',func
           kill_job
         elseif(.not. allocated(irange)) then
           allocate(irange(PGEOM%n_atom))
           irange(1:PGEOM%n_atom) = (/(k, k=1,PGEOM%n_atom)/)
         endif
       else
         if(flag_number(dummy)) then
           call str2int( trim(dummy), init)
         elseif(.not. flag_number(dummy)) then
           call get_site_number(NN_TABLE, PGEOM, dummy, range_dummy, size_range_dummy)
           if(size_range_dummy .ne. 1) then
             if_main write(6,*)'  !WARN! SITERANGE setting argument is not properly asigned. '
             if_main write(6,*)'         Check! Exit...',func
             kill_job
           endif
           init = range_dummy(1)
         endif
         if(flag_number(dummy_)) then
           call str2int( trim(dummy_), fina)
         elseif(.not. flag_number(dummy_)) then
           call get_site_number(NN_TABLE, PGEOM, dummy_, range_dummy, size_range_dummy)
           if(size_range_dummy .ne. 1) then
             if_main write(6,*)'  !WARN! SITERANGE setting argument is not properly asigned. '
             if_main write(6,*)'         Check! Exit...',func
             kill_job
           endif
           fina = range_dummy(1)
         endif
         
         if(.not. allocated(irange)) then
           size_range_dummy = fina - init + 1
           allocate(irange(size_range_dummy))
           irange(1:size_range_dummy) = ((/(k, k=init,fina)/) )
         elseif(allocated(irange)) then
           size_range_dummy = fina - init + 1
           allocate(irange_new(size(irange) + size_range_dummy))
           irange_new = (/irange, (/(k, k=init,fina)/) /)
           deallocate(irange)
           allocate(irange(size(irange_new)))
           irange = irange_new
           deallocate(irange_new)
         endif

       endif

     endif

   enddo

   size_irange = size(irange)

   size_siterange = 0
   do i = 1, size_irange
     do j = 1, size_orbrange
       size_siterange = size_siterange + 1
       siterange_dummy(size_siterange) = orbrange(j) &
                                        +sum(PGEOM%n_orbital(1:irange(i))) - PGEOM%n_orbital(i)
     enddo
   enddo

   if(flag_collinear .or. flag_noncollinear) then
     siterange_dummy = (/siterange_dummy(1:size_siterange), siterange_dummy(1:size_siterange) + PGEOM%neig/)
     size_siterange = size_siterange * 2
   endif

return
endsubroutine
subroutine get_orbrange(strip_orb, PGEOM, orbrange_dummy, size_orbrange)
   use parameters
   use mpi_setup
   type(poscar)  :: PGEOM
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4, intent(out) :: size_orbrange
   integer*4, intent(out) :: orbrange_dummy(PGEOM%max_orb)
   integer*4        nrange
   integer*4        init, fina
   integer*4, allocatable :: irange(:), irange_new(:), orbrange(:)
   character(*)     strip_orb
   character*20,allocatable :: str(:)
   character*20     dummy
   character(*), parameter :: func = 'get_orbrange'
   logical          flag_number, flag_collinear, flag_noncollinear
   external         nitems
   external         flag_number
   integer*4        mpierr

   !orbrange set
   nrange=nitems(strip_orb)
   allocate(str(nrange))
   read(strip_orb,*)str(1:nrange)

   do i=1,nrange

     if(index(str(i),':') .eq. 0) then

       if(flag_number(str(i))) then
         call str2int( trim(str(i)), init )
         fina=init
       elseif(.not. flag_number(str(i))) then
         if_main write(6,*)'  !WARN! ORBRANGE setting argument is not properly asigned. Check! Exit...',func
         kill_job
       endif

     elseif(index(str(i),':') .gt. 0)then

       pos=index(str(i),':')
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then
         dummy='ALL'
       else
         dummy=dummy(1:pos-1)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'ALL') then
           if_main write(6,*)'  !WARN! ORBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), fina )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'ALL') then
           if_main write(6,*)'  !WARN! ORBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     endif

     size_new= size( (/(k, k=init,fina)/) )
     if(.not. allocated(irange)) then
       allocate(irange(size_new))
       irange(:)=( (/(k, k=init,fina)/) )
     elseif(  allocated(irange)) then
       size_now= size(irange)
       allocate(irange_new(size_now+size_new))
       irange_new=(/irange, (/(k, k=init,fina)/) /)
       deallocate(irange)
       allocate(irange(size_now+size_new))
       irange=irange_new
       deallocate(irange_new)
     endif

   enddo

   allocate(orbrange(size(irange)))
   orbrange = irange

   size_orbrange = size(irange)
   orbrange_dummy = 0
   orbrange_dummy(1:size_orbrange) = orbrange

   deallocate(str)
   deallocate(orbrange)
   if(allocated(irange))deallocate(irange)
   if(allocated(irange_new))deallocate(irange_new)


return
endsubroutine
subroutine get_dfrange(strip_df, PWGHT, PINPT, PGEOM, dfrange_dummy, size_dfrange)
   use parameters
   use mpi_setup
   implicit none
   type(incar )  :: PINPT
   type(weight)  :: PWGHT
   type(poscar)  :: PGEOM
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4, intent(out) :: size_dfrange
   integer*4, intent(out) :: dfrange_dummy(PGEOM%neig_target)
   integer*4        nrange
   integer*4        init, fina
   integer*4        mpierr
   integer*4, allocatable :: irange(:), irange_new(:), dfrange(:)
   character(*)     strip_df
   character*20,allocatable :: str(:)
   character*20     dummy
   character(*), parameter :: func = 'get_dfrange'
   logical          flag_number, flag_collinear, flag_noncollinear
   external         nitems
   external         flag_number

   flag_collinear = PINPT%flag_collinear
   flag_noncollinear = PINPT%flag_noncollinear

   !dfrange set
   nrange=nitems(strip_df)
   allocate(str(nrange))
   read(strip_df,*)str(1:nrange)

   do i=1,nrange

     if(index(str(i),':') .eq. 0) then

       if(flag_number(str(i))) then
         call str2int( trim(str(i)), init )
         fina=init
       elseif(.not. flag_number(str(i))) then
         if(trim(str(i)) .eq. 'IBAND') then
           init=PWGHT%iband
           fina=PWGHT%iband
         elseif(trim(str(i)) .eq. 'FBAND') then
           init=PWGHT%fband
           fina=PWGHT%fband
         elseif(trim(str(i)) .eq. 'ALL') then
           init=1
           if(flag_collinear) then
             fina=PGEOM%neig_target/2
           else
             fina=PGEOM%neig_target
           endif
         else
           if_main write(6,*)' 1!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     elseif(index(str(i),':') .gt. 0)then

       pos=index(str(i),':')
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(1:pos-1)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'IBAND') then
           init=PWGHT%iband
         elseif(trim(dummy) .eq. 'FBAND') then
           init=PWGHT%fband
         elseif(trim(dummy) .eq. 'ALL') then
           init=1
         else
           if_main write(6,*)' 2!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), fina )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'IBAND') then
           fina=PWGHT%iband
         elseif(trim(dummy) .eq. 'FBAND') then
           fina=PWGHT%fband
         elseif(trim(dummy) .eq. 'ALL') then
           if(flag_collinear) then
             fina=PGEOM%neig_target/2
           else
             fina=PGEOM%neig_target
           endif
         else
           if_main write(6,*)' 3!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     endif

     size_new= size( (/(k, k=init,fina)/) )
     if(.not. allocated(irange)) then
       allocate(irange(size_new))
       irange(:)=( (/(k, k=init,fina)/) )
     elseif(  allocated(irange)) then
       size_now= size(irange)
       allocate(irange_new(size_now+size_new))
       irange_new=(/irange,(/(k, k=init,fina)/) /)
       deallocate(irange)
       allocate(irange(size_now+size_new))
       irange=irange_new
       deallocate(irange_new)
     endif

   enddo

   if(flag_collinear) then
     allocate( dfrange(size(irange)*2) )
     dfrange = (/irange, irange + (PGEOM%neig_target/2)/)
     
     size_dfrange = size(irange) * 2
     dfrange_dummy = 0
     dfrange_dummy(1:size_dfrange) = (/irange, irange + (PGEOM%neig_target/2)/)

   else
     allocate( dfrange(size(irange)) )
     dfrange = irange

     size_dfrange = size(irange)
     dfrange_dummy = 0
     dfrange_dummy(1:size_dfrange) = (/irange/)


   endif

   deallocate(str)
   if(allocated(irange))    deallocate(irange)
   if(allocated(irange_new))deallocate(irange_new)



return
endsubroutine
subroutine get_tbrange(strip_tb, PGEOM, PINPT, tbrange_dummy, size_tbrange)
   use parameters
   use mpi_setup
   implicit none
   type(incar )  :: PINPT
   type(poscar)  :: PGEOM
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4, intent(out) :: size_tbrange
   integer*4, intent(out) :: tbrange_dummy(PGEOM%neig*PINPT%ispin)
   integer*4        nrange
   integer*4        init, fina
   integer*4, allocatable :: irange(:), irange_new(:), tbrange(:)
   character(*)     strip_tb
   character*20,allocatable :: str(:)
   character*20     dummy
   character(*), parameter :: func = 'get_tbrange'
   logical          flag_number, flag_collinear, flag_noncollinear
   external         nitems
   external         flag_number
   integer*4        mpierr

   flag_collinear = PINPT%flag_collinear
   flag_noncollinear = PINPT%flag_noncollinear

   !tbrange set
   nrange=nitems(strip_tb)
   allocate(str(nrange))
   read(strip_tb,*)str(1:nrange)

   do i=1,nrange

     if(index(str(i),':') .eq. 0) then
       if(flag_number(str(i))) then
         call str2int( trim(str(i)), init )
         fina=init
       elseif(.not. flag_number(str(i))) then
         if(trim(str(i)) .eq. 'NEIG' .or. trim(str(i)) .eq. 'EEND') then
           if(flag_noncollinear)then
             init=PGEOM%neig * 2
           else
             init=PGEOM%neig
           endif
           fina=init
         elseif(trim(str(i)) .eq. 'ALL') then
           init=1
           if(flag_noncollinear) then
             fina=PGEOM%neig * 2
           else
             fina=PGEOM%neig
           endif
         else
           if_main write(6,*)' 1!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     elseif(index(str(i),':') .gt. 0)then

       pos=index(str(i),':')
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(1:pos-1)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NEIG' .or. trim(dummy) .eq. 'EEND') then
           if(flag_noncollinear) then
             init=PGEOM%neig * 2
           else
             init=PGEOM%neig
           endif
         elseif(trim(dummy) .eq. 'ALL')then
           init=1
         else
           if_main write(6,*)' 2!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), fina )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NEIG' .or. trim(dummy) .eq. 'EEND' .or. trim(dummy) .eq. 'ALL') then
           if(flag_noncollinear) then
             fina=PGEOM%neig * 2
           else
             fina=PGEOM%neig
           endif
         else
           if_main write(6,*)' 3!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     endif

     size_new= size( (/(k, k=init,fina)/) )
     if(.not. allocated(irange)) then
       allocate(irange(size_new))
       irange(:)=( (/(k, k=init,fina)/) )
     elseif(  allocated(irange)) then
       size_now= size(irange)
       allocate(irange_new(size_now+size_new))
       irange_new=(/irange, (/(k, k=init,fina)/)/)
       deallocate(irange)
       allocate(irange(size_now+size_new))
       irange=irange_new
       deallocate(irange_new)
     endif

   enddo

   if(flag_collinear) then
     allocate( tbrange(size(irange)*2) )
     tbrange = (/irange, irange+PGEOM%neig/)

     size_tbrange = size(irange) * 2
     tbrange_dummy = 0
     tbrange_dummy(1:size_tbrange) = (/irange, irange + PGEOM%neig/)
   else
     allocate( tbrange(size(irange)) )
     tbrange = irange

     size_tbrange = size(irange)
     tbrange_dummy = 0
     tbrange_dummy(1:size_tbrange) = irange
   endif

   deallocate(str)
   deallocate(tbrange)
   if(allocated(irange))    deallocate(irange)
   if(allocated(irange_new))deallocate(irange_new)

return
endsubroutine

subroutine get_krange(strip_kp, PKPTS, krange_dummy, size_krange)
   use parameters
   use mpi_setup
   implicit none
   type(kpoints) :: PKPTS
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4, intent(out) :: size_krange
   integer*4, intent(out) :: krange_dummy(PKPTS%nkpoint)
   integer*4        nrange
   integer*4        init, fina
   integer*4, allocatable :: irange(:), irange_new(:), krange(:)
   character(*)     strip_kp
   character*20,allocatable :: str(:)
   character*20     dummy
   character(*), parameter :: func = 'get_krange'
   logical          flag_number
   external         nitems
   external         flag_number
   integer*4        mpierr

   !krange set
   nrange=nitems(strip_kp)
   allocate(str(nrange))
   read(strip_kp,*)str(1:nrange)

   do i=1,nrange

     if(index(str(i),':') .eq. 0) then

       if(flag_number(str(i))) then
         call str2int( trim(str(i)), init )
         fina=init
       elseif(.not. flag_number(str(i))) then
         if(trim(str(i)) .eq. 'NKP' .or. trim(str(i)) .eq. 'KEND') then
           init=PKPTS%nkpoint
           fina=init
         elseif(trim(str(i)) .eq. 'ALL') then
           init=1
           fina=PKPTS%nkpoint
         else
           if_main write(6,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     elseif(index(str(i),':') .gt. 0)then

       pos=index(str(i),':')
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then
         dummy='ALL'
       else
         dummy=dummy(1:pos-1)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), init )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NKP' .or. trim(dummy) .eq. 'KEND') then
           init=PKPTS%nkpoint
         elseif(trim(dummy) .eq. 'ALL') then
           init=1
         else
           if_main write(6,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif
       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif
       if(flag_number(dummy)) then
         call str2int( trim(dummy), fina )
       elseif(.not. flag_number(str(i))) then
         if(trim(dummy) .eq. 'NKP' .or. trim(dummy) .eq. 'KEND' .or. trim(dummy) .eq. 'ALL') then
           fina=PKPTS%nkpoint
         else
           if_main write(6,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func
           kill_job
         endif
       endif

     endif
 
     size_new= size( (/(k, k=init,fina)/) )
     if(.not. allocated(irange)) then
       allocate(irange(size_new))
       irange(:)=( (/(k, k=init,fina)/) )
     elseif(  allocated(irange)) then
       size_now= size(irange)
       allocate(irange_new(size_now+size_new))
       irange_new=(/irange, (/(k, k=init,fina)/)/)
       deallocate(irange)
       allocate(irange(size_now+size_new))
       irange=irange_new
       deallocate(irange_new)
     endif

   enddo

   allocate(krange(size(irange)))
   krange = irange

   size_krange = size(irange)
   krange_dummy = 0
   krange_dummy(1:size_krange) = krange

   deallocate(str)
   deallocate(krange)
   if(allocated(irange))deallocate(irange)
   if(allocated(irange_new))deallocate(irange_new)
  
return
endsubroutine

