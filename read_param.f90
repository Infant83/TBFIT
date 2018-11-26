subroutine read_param(PINPT, param_const)
   use parameters, only : incar, max_nparam, pid_param
   implicit none
   type(incar ) :: PINPT
   character(*), parameter :: func = 'read_param'
   integer*4       nitems
   integer*4       i, i_dummy, i_dummy2
   integer*4       i_continue
   real*8          param_const(5,max_nparam)
   character*132   inputline
   character*40    desc_str, dummy, dummy2
   real*8          r_dummy
   logical         flag_skip, flag_number
   external        nitems, flag_number
  
   ! start reading basic information
   open(pid_param,FILE=PINPT%pfilenm,status='old',iostat=i_continue)
   if(i_continue .ne. 0) then
      write(6,'(A,A,A)')'Unknown error opening file:',trim(PINPT%pfilenm),func
      stop
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
     if(i_continue>0) write(6,*)'Unknown error reading file:',trim(PINPT%pfilenm),func
     i=i+1
     call check_comment(inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_empty  (inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
   enddo
   close(pid_param)
   if( i .ne. 0) then
     PINPT%nparam = i
   elseif ( i .eq. 0 ) then
     write(6,'(A,A,A)')'Error in reading ',trim(PINPT%pfilenm),' file. Empty file'
     stop
   endif

   ! reading parameters
   open(pid_param,FILE=PINPT%pfilenm,status='old',iostat=i_continue)
   allocate( PINPT%param(PINPT%nparam) )
   allocate( PINPT%param_name(PINPT%nparam) )
   param_const(1,1:PINPT%nparam) = 0d0      !initialize no matter param_const_ is already provided (PFILE is considered in priori)
   param_const(2,1:PINPT%nparam) = 20d0 ! default upper bound
   param_const(3,1:PINPT%nparam) =-20d0 ! default lower bound
   param_const(4,1:PINPT%nparam) = 0d0
   param_const(5,1:PINPT%nparam) = 0d0
   

   !read parameters if nparam .ne. 0
   i=0
   do
     read(pid_param,'(A)',iostat=i_continue) inputline
     if(i_continue<0) exit               ! end of file reached
     if(i_continue>0) write(6,*)'Unknown error reading file:',trim(PINPT%pfilenm),func
     i=i+1
     call check_comment(inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_empty  (inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     i_dummy = nitems(inputline) - 1

     if(i_dummy .eq. 1) then
       read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i)
       call set_scale_default(param_const(3,i), PINPT%param_name(i))

     elseif( i_dummy .eq. 2) then
       read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i),dummy
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
       read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i),dummy2, dummy
       call set_scale_default(param_const(3,i), PINPT%param_name(i))

       dummy2= trim(dummy2)
       if(.not.flag_number(dummy2)) stop "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file."
       dummy = trim(dummy)
       if(flag_number(dummy))     stop "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file."

       call str2real(dummy2, r_dummy)
       PINPT%param(i) = PINPT%param(i) * r_dummy ! re_scaled

       if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
         param_const(4,i) = 1d0
         param_const(5,i) = PINPT%param(i)
       elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
         param_const(4,i) = 0d0
       endif

     endif

   enddo
   close(pid_param)

return
endsubroutine

subroutine set_scale_default( param_const, param_name )
   implicit none
   character*40  param_name
   real*8        param_const

   if(param_name(1:2) .eq. 's_') then
     param_const = 0.001 ! set the minimum scaling rule
   endif

return
endsubroutine

subroutine set_param_const(PINPT,PGEOM)
   use parameters, only : incar, poscar
   use mpi_setup
   implicit none
   integer*4     i, ii, i_a, i_b
   character*40  dummy
   type(poscar)  :: PGEOM
   type(incar)   :: PINPT

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
   do i = 1, PINPT%nparam_const

     if( trim(PINPT%c_const(2,i)) .eq. '=' ) then
       dummy = trim(PINPT%c_const(3,i))
       if( dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f')then  ! check if fixed
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) i_a = ii
         enddo
         PINPT%param_const(4,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
         PINPT%param_const(5,i_a) = PINPT%param(i_a) ! save its value to PINPT%param_const(5,i_a)
       else
         do ii = 1, PINPT%nparam
           if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) i_a = ii
           if( trim(PINPT%c_const(3,i)) .eq. trim(PINPT%param_name(ii)) ) i_b = ii
         enddo
         PINPT%param_const(1,i_a) = real(i_b)
       endif
       cycle

     elseif( trim(PINPT%c_const(2,i)) .eq. '<=' ) then

       do ii = 1, PINPT%nparam
         if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) i_a = ii
       enddo
       call str2real( trim(PINPT%c_const(3,i)), PINPT%param_const(2,i_a) )
       cycle

     elseif( trim(PINPT%c_const(2,i)) .eq. '>=' ) then

       do ii = 1, PINPT%nparam
         if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) i_a = ii
       enddo
       call str2real( trim(PINPT%c_const(3,i)), PINPT%param_const(3,i_a) )
       cycle

     elseif( trim(PINPT%c_const(2,i)) .eq. '==' ) then

       do ii = 1, PINPT%nparam
         if( trim(PINPT%c_const(1,i)) .eq. trim(PINPT%param_name(ii)) ) i_a = ii
       enddo
       PINPT%param_const(4,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
       PINPT%param_const(5,i_a) = PINPT%param(i_a) ! save its value to PINPT%param_const(5,i_a)
       cycle

     else
       if(myid .eq. 0) write(6,'(A)')'  !WARNING! parameter constraint is not properly defined. Please check again. Exit...'
       stop
     endif
   enddo

return
endsubroutine

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
