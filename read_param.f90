#include "alias.inc"
subroutine read_param(PINPT, PWGHT, param_const, param_const_nrl)
   use parameters, only : incar, weight, max_nparam, pid_param
   use read_incar, only : set_weight_factor
   use mpi_setup
   use print_io
   implicit none
   type(incar ) :: PINPT
   type(weight) :: PWGHT
   character(*), parameter :: func = 'read_param'
   integer*4       nitems
   integer*4       i, i_dummy, i_dummy2
   integer*4       i_continue, i_index
   logical         flag_index
   real*8          param_const(5,max_nparam), param_const_nrl(5,4,max_nparam)
   character*132   inputline
   character*40    desc_str, dummy, dummy2
   character*40    param_name
   real*8          r_dummy
   logical         flag_skip, flag_number
   external        nitems, flag_number
   integer*4       mpierr

   ! start reading "weight" for fitting if PINPT%flag_use_weight
   if(PINPT%flag_use_weight) then
     open(pid_param,FILE=PINPT%pfilenm,status='old',iostat=i_continue)
     if(i_continue .ne. 0) then
        write(message,'(A,A,A)')'Unknown error opening file:',trim(PINPT%pfilenm),func ; write_msg
        kill_job
     else
        PINPT%flag_pfile=.true.
     endif
     dummy = 'WEIGHT'
     call set_weight_factor(PINPT, PWGHT, dummy)
     close(pid_param)
   endif

   ! start reading basic information
   open(pid_param,FILE=PINPT%pfilenm,status='old',iostat=i_continue)
   if(i_continue .ne. 0) then
      write(message,'(A,A,A)')'Unknown error opening file:',trim(PINPT%pfilenm),func ; write_msg
      kill_job
   else
      PINPT%flag_pfile=.true.
   endif

   !check number of parameters
   i=0
   i_dummy=0
   i_dummy2=9999
   do
     read(pid_param,'(A)',iostat=i_continue) inputline
     if(i_continue<0) exit               ! end of file reached
     if(i_continue>0) then
       write(message,*)'Unknown error reading file:',trim(PINPT%pfilenm),func ; write_msg
     endif
     i=i+1
     call check_comment(inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_empty  (inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_options(inputline, PINPT, flag_skip)
     if(flag_skip) then
       i = i - 1; cycle
     endif

   enddo
   close(pid_param)
   if( i .ne. 0) then
     PINPT%nparam = i
   elseif ( i .eq. 0 ) then
     write(message,'(A,A,A)')'Error in reading ',trim(PINPT%pfilenm),' file. Empty file' ; write_msg
     kill_job
   endif

   ! reading parameters
   open(pid_param,FILE=PINPT%pfilenm,status='old',iostat=i_continue)
   allocate( PINPT%param(PINPT%nparam) )
   allocate( PINPT%param_name(PINPT%nparam) )
   allocate( PINPT%param_nsub(PINPT%nparam))

   param_const(1,1:PINPT%nparam) = 0d0  ! is equal to..    !initialize no matter param_const_ is already provided (PFILE is considered in priori)
   param_const(2,1:PINPT%nparam) = 20d0 ! default upper bound
   param_const(3,1:PINPT%nparam) =-20d0 ! default lower bound
   param_const(4,1:PINPT%nparam) = 0d0  ! fixed (1) or not 
   param_const(5,1:PINPT%nparam) = 0d0  ! if fixed, initial value will be stored here
   
   if(PINPT%slater_koster_type .gt. 10) then
     allocate( PINPT%param_nrl(PINPT%param_nsub_max,PINPT%nparam) )
     param_const_nrl(1,:,1:PINPT%nparam) = 0d0      !initialize no matter param_const_ is already provided (PFILE is considered in priori)
     param_const_nrl(2,:,1:PINPT%nparam) = 1112900d0 ! default upper bound
     param_const_nrl(3,:,1:PINPT%nparam) =-1112900d0 ! default lower bound
     param_const_nrl(4,:,1:PINPT%nparam) = 0d0
     param_const_nrl(5,:,1:PINPT%nparam) = 0d0
   endif

   !read parameters if nparam .ne. 0
   i=0
   do
     read(pid_param,'(A)',iostat=i_continue) inputline
     if(i_continue<0) exit               ! end of file reached
     if(i_continue>0) then 
       write(message,*)'Unknown error reading file:',trim(PINPT%pfilenm),func ; write_msg
     endif
     i=i+1
     call check_comment(inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_empty  (inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_options(inputline, PINPT, flag_skip)
     if(flag_skip) then
       i = i - 1; cycle
     endif

     i_dummy = nitems(inputline) - 1
     call check_param_indexing(flag_index, inputline)

     if( flag_index ) then 
       i_dummy = i_dummy - 1
     else
       i_dummy = i_dummy
     endif

     if(PINPT%slater_koster_type .gt. 10) then ! if NRL type SK parameter

       if(i_dummy .eq. 4) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PINPT%param_name(i), PINPT%param_nrl(1:4,i)
         else
           read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i), PINPT%param_nrl(1:4,i)
         endif
       elseif(i_dummy .eq. 5) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PINPT%param_name(i), PINPT%param_nrl(1:4,i), dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i), PINPT%param_nrl(1:4,i), dummy
         endif
         dummy = trim(dummy)
         if(.not. flag_number(dummy) ) then
           if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
             param_const_nrl(4,:,i) = 1d0
             param_const_nrl(5,:,i) = PINPT%param_nrl(:,i)
           elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
             param_const_nrl(4,:,i) = 0d0
           endif
         elseif( flag_number(dummy) ) then
           call str2real(dummy, r_dummy)
           PINPT%param_nrl(:,i) = PINPT%param_nrl(:,i) * r_dummy ! re_scaled
         endif
       elseif(i_dummy .eq. 6) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PINPT%param_name(i), PINPT%param_nrl(1:4,i), dummy2, dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i), PINPT%param_nrl(1:4,i), dummy2, dummy
         endif
         dummy2= trim(dummy2)
         if(.not.flag_number(dummy2)) then
           write(message,'(A)') "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file." ; write_msg
           kill_job
         endif
         dummy = trim(dummy)
         if(flag_number(dummy)) then 
           write(message,'(A)') "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file." ; write_msg
           kill_job
         endif
         call str2real(dummy2, r_dummy)
         PINPT%param_nrl(:,i) = PINPT%param_nrl(:,i) * r_dummy ! re_scaled
         if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
           param_const_nrl(4,:,i) = 1d0
           param_const_nrl(5,:,i) = PINPT%param_nrl(:,i)
         elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
           param_const_nrl(4,:,i) = 0d0
         endif

       elseif( i_dummy .lt. 4) then

         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) param_name
         else
           read(inputline,*,iostat=i_continue) i_index, param_name
         endif
         param_name = trim(param_name)     

         if(param_name(1:2) .ne. 'e_' .and. param_name(1:2) .ne. 'ss' .and. param_name(1:2) .ne. 'sp' .and. &
            param_name(1:2) .ne. 'pp' .and. param_name(1:2) .ne. 'pd' .and. param_name(1:2) .ne. 'dp' .and. &
            param_name(1:2) .ne. 'dd' .and. param_name(1:2) .ne. 's_' .and. param_name(1:2) .ne. 'os' .and. &
            param_name(1:2) .ne. 'o_' ) then
           if(i_dummy .eq. 1) then
             if(.not. flag_index) then
               read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param_nrl(1,i)
             else
               read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i),PINPT%param_nrl(1,i)
             endif
           elseif(i_dummy .eq. 2) then
             if(.not. flag_index) then
               read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param_nrl(1,i), dummy
             else
               read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i),PINPT%param_nrl(1,i), dummy
             endif
             dummy = trim(dummy)
             if(.not. flag_number(dummy) ) then
               if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then
                 param_const_nrl(4,:,i) = 1d0
                 param_const_nrl(5,:,i) = PINPT%param_nrl(:,i)
               elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r') then
                 param_const_nrl(4,:,i) = 0d0
               endif
             elseif( flag_number(dummy) ) then
               call str2real(dummy, r_dummy)
               PINPT%param_nrl(:,i) = PINPT%param_nrl(:,i) * r_dummy
             endif
           elseif( i_dummy .eq. 3) then
             if(.not. flag_index) then
               read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param_nrl(1,i),dummy2, dummy
             else
               read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i),PINPT%param_nrl(1,i),dummy2, dummy
             endif
             dummy2= trim(dummy2)
             if(.not.flag_number(dummy2)) then
               write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msg
               kill_job
             endif
             dummy = trim(dummy)
             if(flag_number(dummy)) then
               write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msg
               kill_job
             endif
             call str2real(dummy2, r_dummy)
             PINPT%param_nrl(:,i) = PINPT%param_nrl(:,i) * r_dummy ! re_scaled
             if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
               param_const_nrl(4,:,i) = 1d0
               param_const_nrl(5,:,i) = PINPT%param_nrl(:,i)
             elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
               param_const_nrl(4,:,i) = 0d0
             endif
           endif

         else
           write(message,'(A)') "  !!WARN!! Wrong syntax in PFILE, check your PFILE."    ; write_msg
           write(message,'(A)') "           You have set SK_SCALE_MODE <= 10, but more than 3 arguments" ; write_msg
           write(message,'(A)') "           has been asigned." ; write_msg
           kill_job
         endif

       else
         write(message,'(A)') "  !!WARN!! Wrong syntax in PFILE, check your PFILE."   ; write_msg
         write(message,'(A)') "           You have set SK_SCALE_MODE <= 10, but more than 3 arguments" ; write_msg
         write(message,'(A)') "           has been asigned." ; write_msg
         kill_job
       endif

       param_name = trim(PINPT%param_name(i))
       if(param_name(1:2) .ne. 'e_' .and. param_name(1:2) .ne. 'ss' .and. param_name(1:2) .ne. 'sp' .and. &
          param_name(1:2) .ne. 'pp' .and. param_name(1:2) .ne. 'pd' .and. param_name(1:2) .ne. 'dp' .and. &
          param_name(1:2) .ne. 'dd' .and. param_name(1:2) .ne. 's_' .and. param_name(1:2) .ne. 'os' .and. &
          param_name(1:2) .ne. 'o_' ) then

         PINPT%param_nsub(i) = 1
         if(param_name(1:4) .eq. 'lons' ) then
           param_const_nrl(2,1,i) = 2d0 ! upper bound
           param_const_nrl(3,1,i) = 0d0 ! lower bound
           param_const_nrl(2,2:4,i) = 0d0 ! upper bound
           param_const_nrl(3,2:4,i) = 0d0 ! lower bound
         endif
       else
         PINPT%param_nsub(i) = 4  ! number of sub parameters (a, b, c, d for NRL hopping, alpha, beta, gamma, xi for NRL onsite)
         if(param_name(1:2) .ne. 'e_') then  ! set bound for parameter 'd' in expenential decay function ( exp(-d**2 * R)
           param_const_nrl(2,4,i) = 2d0 ! upper bound
           param_const_nrl(3,4,i) = 0d0 ! lower bound
         endif
       endif

     else

       if(i_dummy .eq. 1) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i)
         else
           read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i),PINPT%param(i)
         endif
         call set_scale_default(param_const(3,i), PINPT%param_name(i)) 
       elseif( i_dummy .eq. 2) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i),dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i),PINPT%param(i),dummy
         endif
         call set_scale_default(param_const(3,i), PINPT%param_name(i))
         dummy = trim(dummy)
         if(.not. flag_number(dummy) ) then
           if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
             param_const(4,i) = 1d0
             param_const(5,i) = PINPT%param(i)
           elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
             param_const(4,i) = 0d0
           endif
         elseif( flag_number(dummy) ) then
           call str2real(dummy, r_dummy)
           PINPT%param(i) = PINPT%param(i) * r_dummy ! re_scaled
         endif
       elseif( i_dummy .eq. 3) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i),dummy2, dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PINPT%param_name(i),PINPT%param(i),dummy2, dummy
         endif
         call set_scale_default(param_const(3,i), PINPT%param_name(i))
         dummy2= trim(dummy2)
         if(.not.flag_number(dummy2)) then 
           write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msg
           kill_job
         endif
         dummy = trim(dummy)
         if(flag_number(dummy)) then
           write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msg
           kill_job
         endif
         call str2real(dummy2, r_dummy)
         PINPT%param(i) = PINPT%param(i) * r_dummy ! re_scaled
         if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
           param_const(4,i) = 1d0
           param_const(5,i) = PINPT%param(i)
         elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
           param_const(4,i) = 0d0
         endif

       elseif( i_dummy .ge. 4) then
         write(message,'(A)') "  !!WARN!! Wrong syntax in PFILE, check your PFILE." ; write_msg
         write(message,'(A)') "           You have set SK_SCALE_MODE <= 10, but more than 3 arguments" ; write_msg
         write(message,'(A)') "           has been asigned." ; write_msg
         kill_job
       endif
       param_name = trim(PINPT%param_name(i))
       if(param_name(1:2) .ne. 'e_' .and. param_name(1:2) .ne. 'ss' .and. param_name(1:2) .ne. 'sp' .and. &
          param_name(1:2) .ne. 'pp' .and. param_name(1:2) .ne. 'pd' .and. param_name(1:2) .ne. 'dp' .and. &
          param_name(1:2) .ne. 'dd' .and. param_name(1:2) .ne. 's_' .and. param_name(1:2) .ne. 'os' .and. &
          param_name(1:2) .ne. 'o_' ) then

         PINPT%param_nsub(i) = 1
       else
         PINPT%param_nsub(i) = 1
       endif
     endif

   enddo
   close(pid_param)

return
endsubroutine

subroutine set_scale_default( param_const, param_name)
   implicit none
   character*40  param_name
   real*8        param_const
!  integer*4     slater_koster_type

!  if(slater_koster_type .gt. 10) then

!    if(param_name(1:2) .ne. 'e_' .and.  param_name(1:2) .ne. 'la' .and. &
!       param_name(1:2) .ne. 'st' .and.  param_name(1:2) .ne. 'lo' .and. &
!       param_name(1:2) .ne. 'lr' ) then

!      param_const = 0.033 ! set the lower bound of scaling factor if param_name start with "sss, sps, pps, .. etc" or
!                          ! "o_sss, o_sps, o_pps, .. etc".
!    endif
!  else
     if(param_name(1:2) .eq. 's_') then
       param_const = 0.001 ! set the lower bound of scaling factor
     elseif(param_name(1:3) .eq. 'os_') then
       param_const = 0.001 ! set the lower bound of scaling factor
     endif
!  endif

return
endsubroutine

subroutine set_param_const(PINPT,PGEOM)
   use parameters, only : incar, poscar
   use mpi_setup
   use print_io
   implicit none
   integer*4     i, ii, i_a, i_b
   integer*4     j
   character*40  dummy
   type(poscar)  :: PGEOM
   type(incar)   :: PINPT
   integer*4     mpierr

!PINPT%param_const(i,:) i=1 -> is same as
!                       i=2 -> is lower than (.le.) : set maximum bound  ! functionality is not available yet
!                       i=3 -> is lower than (.ge.) : set minimum bound  ! functionality is not available yet
!                       i=4 -> is fixed : not to be fitted, just stay    ! its original value will be stored in PINPT%param_const(i=5,:)

!  allocate( PINPT%param_const(5,PINPT%nparam) )
!  PINPT%param_const(1,:) = 0d0  ! same as 
!  PINPT%param_const(2,:) = 20d0 ! upper bound
!  PINPT%param_const(3,:) =-20d0 ! lower bound ; if start with s_ (scale factor) then >= 0.001 (default)
!  PINPT%param_const(4,:) = 0d0  ! fixed 1:true 0:no
!  PINPT%param_const(5,:) = 0d0  ! fixed value

!  NOTE: same rule applies for PINPT%param_const_nrl(1:5,:,:) as PINPT%param_const(1:5,:)

   if(PINPT%slater_koster_type .le. 10) then
     do i = 1, PINPT%nparam_const
       if( trim(PINPT%c_const(2,i)) .eq. '=' ) then
         dummy = trim(PINPT%c_const(3,i))
         if( dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f')then  ! check if fixed
           do ii = 1, PINPT%nparam
             if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then 
               i_a = ii
               PINPT%param_const(4,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
               PINPT%param_const(5,i_a) = PINPT%param(i_a) ! save its value to PINPT%param_const(5,i_a)
             endif
           enddo
         else
           i_a = 0
           do ii = 1, PINPT%nparam
             if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) i_a = ii
           enddo
           do ii = 1, PINPT%nparam
             if( trim(PINPT%c_const(3,i)) .eq. trim(PINPT%param_name(ii)) .and. i_a .ne. 0 ) then 
               i_b = ii
               PINPT%param_const(1,i_a) = real(i_b)
             endif
           enddo
         endif
         cycle
       elseif( trim(PINPT%c_const(2,i)) .eq. '<=' ) then
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then 
             i_a = ii
             call str2real( trim(PINPT%c_const(3,i)), PINPT%param_const(2,i_a) )
           endif
         enddo
         cycle
       elseif( trim(PINPT%c_const(2,i)) .eq. '>=' ) then
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then 
             i_a = ii
             call str2real( trim(PINPT%c_const(3,i)), PINPT%param_const(3,i_a) )
           endif
         enddo
         cycle
       elseif( trim(PINPT%c_const(2,i)) .eq. '==' ) then
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then 
             i_a = ii
             PINPT%param_const(4,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
             PINPT%param_const(5,i_a) = PINPT%param(i_a) ! save its value to PINPT%param_const(5,i_a)
           endif
         enddo
         cycle
       else
         write(message,'(A)')'  !WARNING! parameter constraint is not properly defined. Please check again. Exit...' ; write_msg
         kill_job
       endif
     enddo
  
   elseif(PINPT%slater_koster_type .gt. 10) then
     do i = 1, PINPT%nparam_const
       if( trim(PINPT%c_const(2,i)) .eq. '=' ) then
         dummy = trim(PINPT%c_const(3,i))
         if( dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f')then  ! check if fixed
           do ii = 1, PINPT%nparam
             if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then
               i_a = ii
               PINPT%param_const_nrl(4,:,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
               PINPT%param_const_nrl(5,:,i_a) = PINPT%param_nrl(:,i_a) ! save its value to PINPT%param_const(5,i_a)
             endif
           enddo
         else
           i_a = 0
           do ii = 1, PINPT%nparam
             if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) i_a = ii
           enddo
           do ii = 1, PINPT%nparam
             if( trim(PINPT%c_const(3,i)) .eq. trim(PINPT%param_name(ii)) .and. i_a .ne. 0 ) then
               i_b = ii
               PINPT%param_const_nrl(1,:,i_a) = real(i_b)
             endif
           enddo
         endif
         cycle
       elseif( trim(PINPT%c_const(2,i)) .eq. '<=' ) then  ! set upper bound
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then
             i_a = ii
             do j = 1, PINPT%param_nsub_max
               call str2real( trim(PINPT%c_const(3,i)), PINPT%param_const_nrl(2,j, i_a) )
             enddo
           endif
         enddo
         cycle
       elseif( trim(PINPT%c_const(2,i)) .eq. '>=' ) then
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then
             i_a = ii
             do j = 1, PINPT%param_nsub_max
               call str2real( trim(PINPT%c_const(3,i)), PINPT%param_const_nrl(3,j, i_a) )
             enddo
           endif
         enddo
         cycle
       elseif( trim(PINPT%c_const(2,i)) .eq. '==' ) then
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) then
             i_a = ii
             PINPT%param_const_nrl(4,:, i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
             PINPT%param_const_nrl(5,:, i_a) = PINPT%param_nrl(:,i_a) ! save its value to PINPT%param_const(5,i_a)
           endif
         enddo
         cycle
       else
         write(message,'(A)')'  !WARNING! parameter constraint is not properly defined. Please check again. Exit...' ; write_msg
         kill_job
       endif
     enddo
   endif

return
endsubroutine

! update with constraint values which is redirect with (1,i) indexing stored in param_const
subroutine update_param(PINPT)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4      i

   do i = 1, PINPT%nparam
     if( nint(PINPT%param_const(1,i)) .ge. 1) then ! 
       PINPT%param(i) = PINPT%param( nint(PINPT%param_const(1,i)) )
     endif
   enddo

return
endsubroutine
subroutine update_param_nrl(PINPT)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4      i, j, k

   do j = 1, PINPT%nparam
     k = nint(PINPT%param_const_nrl(1,1,j))
     if( k .ge. 1) then ! 
       do i = 1,PINPT%param_nsub(j)
         PINPT%param_nrl(i,j) = PINPT%param_nrl(i, k)
       enddo
     endif
   enddo

return
endsubroutine

subroutine check_param_indexing(flag_index,inputline)
   implicit none
   logical              flag_index
   character*40         dummy
   character*132        inputline
   logical, external :: flag_number

   read(inputline, *) dummy
   
   if( flag_number(trim(dummy)) ) then
     flag_index = .true.
   else
     flag_index = .false.
   endif

endsubroutine

subroutine check_options(inputline, PINPT, flag_skip)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   logical        flag_skip  
   character*132  inputline
   character*40   dummy, dummy1

   read(inputline, *) dummy

   ! read options if it is predefined, otherwise it will be recognized as 'parameters'
   if( trim(dummy) .eq. 'PRINT_INDEX' .or. trim(dummy) .eq. 'print_index') then
     read(inputline, *) dummy1, PINPT%flag_pfile_index
     flag_skip = .true.
   elseif( trim(dummy) .eq. 'USE_OVERLAP' .or. trim(dummy) .eq. 'use_overlap' ) then
     read(inputline, *) dummy1, PINPT%flag_use_overlap
     flag_skip = .true.
   elseif( trim(dummy) .eq. 'SK_SCALE_MODE' .or. trim(dummy) .eq. 'SCALE_MODE' .or. &
           trim(dummy) .eq. 'sk_scale_mode' .or. trim(dummy) .eq. 'scale_mode') then
     read(inputline, *) dummy1, PINPT%slater_koster_type 
     flag_skip = .true. 
     if(PINPT%slater_koster_type .gt. 10) then
       PINPT%param_nsub_max = 4
     else
       PINPT%param_nsub_max = 1
     endif
   elseif( trim(dummy) .eq. 'L_BROADEN' .or. trim(dummy) .eq. 'l_broaden' ) then
     read(inputline, *) dummy1, PINPT%l_broaden
     flag_skip = .true.
   else
     flag_skip = .false.
   endif

   return
endsubroutine

subroutine check_print_tag(inputline, flag_pfile_index, flag_skip)
   implicit none
   logical        flag_pfile_index, flag_skip
   character*132  inputline
   character*40   dummy, dummy1

   read(inputline, *) dummy

   if( trim(dummy) .eq. 'PRINT_INDEX' ) then
     read(inputline, *) dummy1, flag_pfile_index
     flag_skip = .true.
   else
     flag_skip = .false.
   endif

endsubroutine

