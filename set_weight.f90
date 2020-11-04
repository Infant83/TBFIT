#include "alias.inc"
subroutine set_penalty_orb(PINPT, PGEOM, PKPTS, PWGHT, strip_kp, strip_tb, strip_orb, strip_site, strip_pen)
   use parameters
   use mpi_setup
   use print_io
   implicit none
   type(incar )  :: PINPT
   type(poscar)  :: PGEOM
   type(kpoints) :: PKPTS
   type(weight)  :: PWGHT
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
   character*40     fm_ob
   character*20,external :: int2str

   if( PWGHT%max_len_strip_ob .gt. 1)then
     write(fm_ob,'( "(A22,A",A,",A16,A",A,",A17,A",A,",A18,A",A,",A13,A10)" )') &
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_kp))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_tb))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_ob))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_st)))
   else
     write(fm_ob,'( "(A22,A10",",A16,A10",",A17,A10",",A18,A10",",A13,A10)" )')
   endif

   if_main_then
     write(message,fm_ob)' O_PENALTY: KRANGE= [ ',trim(strip_kp), &
                                 ' ], TB_RANGE= [ ',trim(strip_tb), &
                                 ' ], ORB_RANGE= [ ',trim(strip_orb), &
                                 ' ], SITE_RANGE= [ ', trim(strip_site), &
                                 ' ], PENALTY= ',trim(strip_pen) ; write_msg
   if_main_end
   !krange set
   call get_krange(strip_kp, PKPTS, krange, size_krange)

   !tbrange set
   call get_tbrange(strip_tb, PWGHT, PGEOM, PINPT, tbrange, size_tbrange)

   !orbrange set
   call get_orbrange(strip_orb, PGEOM, orbrange, size_orbrange)

   !siterange set
   call get_siterange(strip_site, PGEOM, PINPT, orbrange(1:size_orbrange), size_orbrange, siterange, size_siterange)

   !weight set
   if(flag_number(trim(strip_pen))) then
     call str2real(trim(strip_pen),pen_orb)
   elseif(.not. flag_number(trim(strip_pen))) then
     write(message,*)'  !WARN! PENALTY setting argument is not properly asigned. '  ; write_msg
     write(message,*)'         Check! Exit...',func  ; write_msg
     kill_job
   endif

   PWGHT%PENALTY_ORB(siterange(1:size_siterange),tbrange(1:size_tbrange),krange(1:size_krange))=pen_orb

return
endsubroutine

subroutine set_weight(PINPT, PGEOM, PKPTS, PWGHT, EDFT, EDFT_all, strip_kp, strip_tb, strip_df, strip_wt, WT, imode)
   use parameters
   use mpi_setup
   use print_io
   implicit none
   type(incar )           :: PINPT
   type(poscar)           :: PGEOM
   type(kpoints)          :: PKPTS
   type(weight)           :: PWGHT
   type(energy)           :: EDFT
   type(energy)           :: EDFT_all
   integer*4                 i, k, nitems, pos, size_now, size_new
   integer*4                 ikp
   integer*4                 nrange
   integer*4                 init,fina
   integer*4                 imode
   integer*4                 size_krange, size_tbrange, size_dfrange
   integer*4, allocatable :: irange(:), irange_new(:)
   integer*4                 krange(PKPTS%nkpoint)
   integer*4                 tbrange(PGEOM%neig*PINPT%ispin)
   integer*4                 dfrange(PGEOM%neig_target)
   real*8                    w
   character(*)              strip_kp, strip_tb, strip_df, strip_wt
   character*20              dummy, strk
   logical                   flag_number, flag_collinear, flag_noncollinear
   external                  nitems   
   external                  flag_number
   integer*4                 mpierr
   real*8                    WT(PGEOM%neig*PINPT%ispin, PKPTS%nkpoint)
   character*40              fm_wt
   character*20,allocatable :: str(:)
   character(*), parameter  :: func = 'set_weight'
   character*20,external    :: int2str

   flag_collinear = PINPT%flag_collinear
   flag_noncollinear = PINPT%flag_noncollinear
   write(fm_wt,'( "(A21,A",A,",A16,A",A,",A16,A",A,",A12,A10)" )') trim(ADJUSTL(int2str(PWGHT%max_len_strip_kp))),&
                                                               trim(ADJUSTL(int2str(PWGHT%max_len_strip_tb))),&
                                                               trim(ADJUSTL(int2str(PWGHT%max_len_strip_df)))

   if_main_then
     if(imode .eq. 1) then
       write(message,fm_wt) ' C_WEIGHT: KRANGE= [ ',trim(strip_kp), &
                                 ' ], TB_RANGE= [ ',trim(strip_tb), &
                                 ' ], DF_RANGE= [ ',trim(strip_df), &
                                 ' ], WEIGHT= ',trim(strip_wt) ; write_msg
     elseif(imode .eq. 2) then
       write(message,fm_wt) ' C_WEIGHT: KRANGE= [ ',trim(strip_kp), &
                                 ' ], TB_RANGE= [ ',trim(strip_tb), &
                                 ' ], DF_RANGE= [ ',trim(strip_df), &
                                 ' ], DEGENW= ',trim(strip_wt) ; write_msg
     endif
   if_main_end

   !krange set
   call get_krange(strip_kp, PKPTS, krange, size_krange)   

   !tbrange set
   call get_tbrange(strip_tb, PWGHT, PGEOM, PINPT, tbrange, size_tbrange)

   !dfrange set
   call get_dfrange(strip_df, PWGHT, PINPT, PGEOM, dfrange, size_dfrange)
 
   !weight set
   call str2real(trim(strip_wt),w)

   !reset target energy & weight
   if( size_tbrange .ne. size_dfrange ) then
     write(message,'(A,A)')'  !WARN! size of TBRANGE and DFRANGE is not same. Please check "SET WEIGHT" tag. Exit...',func  ; write_msg
     kill_job
   endif
   EDFT%E(tbrange(1:size_tbrange),krange(1:size_krange)) = EDFT_all%E(dfrange(1:size_dfrange),krange(1:size_krange))
   WT(tbrange(1:size_tbrange),krange(1:size_krange)) = w

return
endsubroutine

subroutine get_siterange(strip_site, PGEOM, PINPT, orbrange, size_orbrange, siterange_dummy, size_siterange)
   use parameters
   use mpi_setup
   use print_io
   implicit none
   type(poscar)  :: PGEOM
   type(incar )  :: PINPT
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
         call get_site_number(PGEOM, dummy, range_dummy, size_range_dummy)
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
           write(message,*)'  !WARN! SITERANGE setting argument is not properly asigned. '  ; write_msg
           write(message,*)'         Check! Exit...',func  ; write_msg
           kill_job
         elseif(.not. allocated(irange)) then
           allocate(irange(PGEOM%n_atom))
           irange(1:PGEOM%n_atom) = (/(k, k=1,PGEOM%n_atom)/)
         endif
       else
         if(flag_number(dummy)) then
           call str2int( trim(dummy), init)
         elseif(.not. flag_number(dummy)) then
           call get_site_number(PGEOM, dummy, range_dummy, size_range_dummy)
           if(size_range_dummy .ne. 1) then
             write(message,*)'  !WARN! SITERANGE setting argument is not properly asigned. '  ; write_msg
             write(message,*)'         Check! Exit...',func  ; write_msg
             kill_job
           endif
           init = range_dummy(1)
         endif
         if(flag_number(dummy_)) then
           call str2int( trim(dummy_), fina)
         elseif(.not. flag_number(dummy_)) then
           call get_site_number(PGEOM, dummy_, range_dummy, size_range_dummy)
           if(size_range_dummy .ne. 1) then
             write(message,*)'  !WARN! SITERANGE setting argument is not properly asigned. '  ; write_msg
             write(message,*)'         Check! Exit...',func  ; write_msg
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
   use print_io
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
         write(message,*)'  !WARN! ORBRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
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
           write(message,*)'  !WARN! ORBRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
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
           write(message,*)'  !WARN! ORBRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
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
   use print_io
   implicit none
   type(incar )  :: PINPT
   type(weight)  :: PWGHT
   type(poscar)  :: PGEOM
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4, intent(out) :: size_dfrange
   integer*4, intent(out) :: dfrange_dummy(PGEOM%neig_target)
   integer*4        nrange
   integer*4        init, fina, iband
   integer*4        mpierr
   integer*4, allocatable :: irange(:), irange_new(:), dfrange(:)
   character(*)     strip_df
   character*20,allocatable :: str(:)
   character*20     dummy, strdf
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
         strdf = trim(str(i))
         call find_band_index(iband, PWGHT, strdf)

         if(iband .eq. 0) then
           write(message,*)' 1!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
           kill_job
         elseif(iband .ge. 1) then ! single band specified
           init = iband
           fina = iband
         elseif(iband .le. -1) then ! ALL
           init = 1
           if(flag_collinear) then
             fina=PGEOM%neig_target/2
           else
             fina=PGEOM%neig_target
           endif
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

       strdf = trim(dummy)
       call find_band_index(iband, PWGHT, strdf)
       if(iband .eq. 0) then
         write(message,*)' 2!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
         kill_job
       elseif(iband .ge. 1) then
         init = iband
       elseif(iband .le. -1) then
         init = 1
       endif      

       dummy=str(i)
       if(len_trim(dummy) .eq. 1) then ! if only ':' mark is present
         dummy='ALL'
       else
         dummy=dummy(pos+1:)
       endif

       strdf = trim(dummy)
       call find_band_index(iband, PWGHT, strdf)
       if(iband .eq. 0) then
         write(message,*)' 3!WARN! DFRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
         kill_job
       elseif(iband .ge. 1) then
         fina = iband
       elseif(iband .le. -1) then
         if(flag_collinear) then
           fina=PGEOM%neig_target/2
         else
           fina=PGEOM%neig_target
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
subroutine get_tbrange(strip_tb, PWGHT, PGEOM, PINPT, tbrange_dummy, size_tbrange)
   use parameters
   use mpi_setup
   use print_io
   implicit none
   type(incar )  :: PINPT
   type(poscar)  :: PGEOM
   type(weight)  :: PWGHT
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4, intent(out) :: size_tbrange
   integer*4, intent(out) :: tbrange_dummy(PGEOM%neig*PINPT%ispin)
   integer*4        nrange
   integer*4        init, fina, iband
   integer*4, allocatable :: irange(:), irange_new(:), tbrange(:)
   character(*)     strip_tb
   character*20,allocatable :: str(:)
   character*20     dummy, strtb
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
           strtb = trim(str(i))
           call find_band_index(iband, PWGHT, strtb)

           if(iband .eq. 0) then
             write(message,*)' 1!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
             kill_job
           elseif(iband .ge. 1) then
             init = iband
             fina = iband
           endif
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
           strtb = trim(dummy)
           call find_band_index(iband, PWGHT, strtb)

           if(iband .eq. 0) then
             write(message,*)' 2!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
             kill_job
           elseif(iband .ge. 1) then
             init = iband
           endif

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
           strtb = trim(dummy)
           call find_band_index(iband, PWGHT, strtb)

           if(iband .eq. 0) then
             write(message,*)' 3!WARN! TBRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
             kill_job
           elseif(iband .ge. 1) then
             fina = iband
           endif
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
   use print_io
   implicit none
   type(kpoints) :: PKPTS
   integer*4        i, k, nitems, pos, size_now, size_new
   integer*4, intent(out) :: size_krange
   integer*4, intent(out) :: krange_dummy(PKPTS%nkpoint)
   integer*4        nrange
   integer*4        init, fina, ikp
   integer*4, allocatable :: irange(:), irange_new(:), krange(:)
   character(*)     strip_kp
   character*20,allocatable :: str(:)
   character*20     dummy, strk
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
           strk = trim(str(i))
           call find_kp_index(ikp, PKPTS, strk)
 
           if(ikp .eq. 0) then          
             write(message,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
             kill_job
           else
             init = ikp
             fina = init
           endif
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
           strk = trim(dummy)
           call find_kp_index(ikp, PKPTS, strk)
           
           if(ikp .eq. 0) then
             write(message,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
             kill_job
           else
             init = ikp
           endif
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
           strk = trim(dummy)
           call find_kp_index(ikp, PKPTS, strk)
 
           if(ikp .eq. 0) then          
             write(message,*)'  !WARN! KRANGE setting argument is not properly asigned. Check! Exit...',func  ; write_msg
             kill_job
           else
             fina = ikp
           endif
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

subroutine find_kp_index(ikp, PKPTS, strk)
   use parameters, only: kpoints
   implicit none
   type(kpoints) :: PKPTS
   integer*4        ikp, iline
   integer*4        idummy, pos, iadd
   character*20     strk, dummy
   character*8      kname
   ikp = 0

   ! NOTE: this routine only valid if PKPTS%flag_klinemode = .true.

   do iline = 1, PKPTS%nline+1
     if(iline .ne. 1) then
       kname = trim(PKPTS%k_name( (iline-1)*2 ))

       if( index(trim(strk),trim(kname)) .ge. 1) then
         if(PKPTS%idiv_mode .eq. 3) then ! aims-like
           ikp = sum( PKPTS%ndiv(1:iline-1))+1
         elseif(PKPTS%idiv_mode .eq. 1) then ! vasp-like
           ikp = PKPTS%ndiv(1)*(iline-1)+1
         elseif(PKPTS%idiv_mode .eq. 2) then ! FLEUR-like
           ikp = PKPTS%ndiv(1)*(iline-1) ! NOTE: need to be checked whether it is write
         endif

         if(index(strk,'+') .ge. 1 .or. index(strk,'-') .ge. 1) then
           pos = index(strk,'+')
           if(pos .eq. 0) pos = index(strk,'-')
           dummy = strk(pos:)
           call str2int( trim(dummy), iadd )
           ikp = ikp + iadd
         endif

       endif

     elseif(iline .eq. 1) then
       kname = trim(PKPTS%k_name(iline))

       if( index(trim(strk),trim(kname)) .ge. 1) then
         ikp = 1
       endif

       if(index(strk,'+') .ge. 1 .or. index(strk,'-') .ge. 1) then
         pos = index(strk,'+')
         if(pos .eq. 0) pos = index(strk,'-')
         dummy = strk(pos:)
         call str2int( trim(dummy), iadd )
         ikp = ikp + iadd
       endif
     endif
   enddo

return
endsubroutine

subroutine find_band_index(iband, PWGHT, strband)
   use parameters, only : weight
   implicit none
   type(weight)  :: PWGHT
   integer*4    iband
   integer*4    idummy, pos, iadd
   character*20 strband, str, dummy
   integer*4    mpierr
   logical      flag_number
   external     flag_number

   iband = 0
   iadd  = 0

   if(index(strband,'+') .ge. 1 .or. index(strband,'-') .ge. 1) then
     pos = index(strband,'+')
     if(pos .eq. 0) pos = index(strband,'-')

     dummy = strband(pos:)
     call str2int( trim(dummy), iadd)
     str = strband(1:pos-1)
   else
     str = trim(strband)
   endif

   if(flag_number(str)) then
     call str2int(trim(str),iband)
   else
     select case (trim(str))
     
       case('VBMD')
         iband = PWGHT%vbmd + iadd
       
       case('VBMT')
         iband = PWGHT%vbmt + iadd
       
       case('CBMD')
         iband = PWGHT%cbmd + iadd
     
       case('CBMT')
         iband = PWGHT%cbmt + iadd
     
       case('IBAND')
         iband = PWGHT%iband+ iadd
     
       case('FBAND')
         iband = PWGHT%fband+ iadd
     
       case('ALL') 
         iband = -1
     
     endselect
   endif

return
endsubroutine
