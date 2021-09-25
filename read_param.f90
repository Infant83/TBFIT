#include "alias.inc"
subroutine read_tb_param(PINPT, PPRAM, PWGHT )
   use parameters, only : incar, weight, max_nparam, pid_param, params
   use read_incar, only : set_weight_factor
   use mpi_setup
   use print_io
   use random_mod
   implicit none
   type(incar ) :: PINPT
   type(weight) :: PWGHT
   type(params) :: PPRAM
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
   logical         flag_skip, flag_number, flag_exist
   external        nitems, flag_number
   integer*4       mpierr
   integer*4       iseed

  !iseed = 9342
  !call random_init(iseed)
   call random_seed()

   inquire(file=trim(PPRAM%pfilenm),exist=flag_exist)
   if(.not. flag_exist) then
     write(message,'(A,A,A)')'    !WARN! Parameter file:',trim(PPRAM%pfilenm),' does not exist!! Exit...' ; write_msgi
     kill_job
   endif

   ! start reading "weight" for fitting if PINPT%flag_use_weight
   if(PINPT%flag_use_weight) then
     open(pid_param,FILE=PPRAM%pfilenm,status='old',iostat=i_continue)
     dummy = 'WEIGHT'
     call set_weight_factor(PINPT, PWGHT, dummy)
     close(pid_param)
   endif

   ! start reading basic information
   write(message,*)' '  ; write_msgi
   write(message,*)'#- READING INPUT PARAMER FILE: ',trim(PPRAM%pfilenm)  ; write_msgi

   open(pid_param,FILE=PPRAM%pfilenm,status='old',iostat=i_continue)

   !check number of parameters
   i=0
   i_dummy=0
   i_dummy2=9999
   do
     read(pid_param,'(A)',iostat=i_continue) inputline
     if(i_continue<0) exit               ! end of file reached
     if(i_continue>0) then
       write(message,*)'Unknown error reading file:',trim(PPRAM%pfilenm),func ; write_msgi
     endif
     i=i+1
     call check_comment(inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_empty  (inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_options(inputline, PPRAM, flag_skip)
     if(flag_skip) then
       i = i - 1; cycle
     endif

   enddo
   close(pid_param)
   if( i .ne. 0) then
     PPRAM%nparam = i
   elseif ( i .eq. 0 ) then
     write(message,'(A,A,A)')'Error in reading ',trim(PPRAM%pfilenm),' file. Empty file' ; write_msgi
     kill_job
   endif

   ! reading parameters
   open(pid_param,FILE=PPRAM%pfilenm,status='old',iostat=i_continue)
   allocate( PPRAM%param(PPRAM%nparam) )
   allocate( PPRAM%param_best(PPRAM%nparam) )
   allocate( PPRAM%param_name(PPRAM%nparam) )
   allocate( PPRAM%param_nsub(PPRAM%nparam))
   
   PPRAM%param(1:PPRAM%nparam  ) = 0d0
   PPRAM%param_best(1:PPRAM%nparam  ) = 0d0

  !param_const(:,0             ) = 0d0
   param_const(1,1:PPRAM%nparam) = 0d0  ! is equal to..    !initialize no matter param_const_ is already provided (PFILE is considered in priori)
   param_const(2,1:PPRAM%nparam) = 20d0 ! default upper bound
   param_const(3,1:PPRAM%nparam) =-20d0 ! default lower bound
   param_const(4,1:PPRAM%nparam) = 0d0  ! fixed (1) or not 
   param_const(5,1:PPRAM%nparam) = 0d0  ! if fixed, initial value will be stored here
   
   if(PPRAM%slater_koster_type .gt. 10) then
     allocate( PPRAM%param_nrl(PPRAM%param_nsub_max,PPRAM%nparam) )

     PPRAM%param_nrl(:,1:PPRAM%nparam  ) = 0d0

!    param_const_nrl(:,:,0             ) = 0d0      !initialize no matter param_const_ is already provided (PFILE is considered in priori)
     param_const_nrl(1,:,1:PPRAM%nparam) = 0d0      !initialize no matter param_const_ is already provided (PFILE is considered in priori)
     param_const_nrl(2,:,1:PPRAM%nparam) = 1112900d0 ! default upper bound
     param_const_nrl(3,:,1:PPRAM%nparam) =-1112900d0 ! default lower bound
     param_const_nrl(4,:,1:PPRAM%nparam) = 0d0
     param_const_nrl(5,:,1:PPRAM%nparam) = 0d0
   endif

   !read parameters if nparam .ne. 0
   i=0
   do
     read(pid_param,'(A)',iostat=i_continue) inputline
     if(i_continue<0) exit               ! end of file reached
     if(i_continue>0) then 
       write(message,*)'Unknown error reading file:',trim(PPRAM%pfilenm),func ; write_msgi
     endif
     i=i+1
     call check_comment(inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_empty  (inputline,i_dummy2,i,flag_skip) ; if(flag_skip) cycle
     call check_options(inputline, PPRAM, flag_skip)
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

     if(PPRAM%slater_koster_type .gt. 10) then ! if NRL type SK parameter

       if(i_dummy .eq. 4) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PPRAM%param_name(i), PPRAM%param_nrl(1:4,i)
         else
           read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i), PPRAM%param_nrl(1:4,i)
         endif
       elseif(i_dummy .eq. 5) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PPRAM%param_name(i), PPRAM%param_nrl(1:4,i), dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i), PPRAM%param_nrl(1:4,i), dummy
         endif
         dummy = trim(dummy)
         if(.not. flag_number(dummy) ) then
           if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
             param_const_nrl(4,:,i) = 1d0
             param_const_nrl(5,:,i) = PPRAM%param_nrl(:,i)
           elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'random'
         !   param_const_nrl(4,:,i) = 0d0
             r_dummy = random()*2d0 - 1.d0
#ifdef MPI
             call MPI_BCAST(r_dummy, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif
             PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy
           endif
         elseif( flag_number(dummy) ) then
           call str2real(dummy, r_dummy)
           PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy ! re_scaled
         endif
       elseif(i_dummy .eq. 6) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PPRAM%param_name(i), PPRAM%param_nrl(1:4,i), dummy2, dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i), PPRAM%param_nrl(1:4,i), dummy2, dummy
         endif
         dummy2= trim(dummy2)
         if(.not.flag_number(dummy2) .and. (dummy2(1:1) .ne. 'r' .and. dummy2(1:1) .ne. 'R')) then
           write(message,'(A)') "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file." ; write_msgi
           kill_job
         endif
         dummy = trim(dummy)
         if(flag_number(dummy)) then 
           write(message,'(A)') "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file." ; write_msgi
           kill_job
         endif
         if(dummy2(1:1) .ne. 'r' .and. dummy2(1:1) .ne. 'R') then
           call str2real(dummy2, r_dummy)
           PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy ! re_scaled
         elseif(dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then
           r_dummy = random()*2d0 - 1.d0
#ifdef MPI
           call MPI_BCAST(r_dummy, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif
           PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy + r_dummy
         endif
         if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
           param_const_nrl(4,:,i) = 1d0
           param_const_nrl(5,:,i) = PPRAM%param_nrl(:,i)
!        elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'random'
!          param_const_nrl(4,:,i) = 0d0
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
               read(inputline,*,iostat=i_continue) PPRAM%param_name(i),PPRAM%param_nrl(1,i)
             else
               read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i),PPRAM%param_nrl(1,i)
             endif
           elseif(i_dummy .eq. 2) then
             if(.not. flag_index) then
               read(inputline,*,iostat=i_continue) PPRAM%param_name(i),PPRAM%param_nrl(1,i), dummy
             else
               read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i),PPRAM%param_nrl(1,i), dummy
             endif
             dummy = trim(dummy)
             if(.not. flag_number(dummy) ) then
               if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then
                 param_const_nrl(4,:,i) = 1d0
                 param_const_nrl(5,:,i) = PPRAM%param_nrl(:,i)
               elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r') then
!                param_const_nrl(4,:,i) = 0d0
                 r_dummy = random()*2d0 - 1.d0
#ifdef MPI      
                 call MPI_BCAST(r_dummy, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif          
                 PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy + r_dummy
               endif
             elseif( flag_number(dummy) ) then
               call str2real(dummy, r_dummy)
               PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy
             endif
           elseif( i_dummy .eq. 3) then
             if(.not. flag_index) then
               read(inputline,*,iostat=i_continue) PPRAM%param_name(i),PPRAM%param_nrl(1,i),dummy2, dummy
             else
               read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i),PPRAM%param_nrl(1,i),dummy2, dummy
             endif
             dummy2= trim(dummy2)
             if(.not.flag_number(dummy2) .and. (dummy2(1:1) .ne. 'r' .and. dummy2(1:1) .ne. 'R')) then
               write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msgi
               kill_job
             endif
             dummy = trim(dummy)
             if(flag_number(dummy)) then
               write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msgi
               kill_job
             endif
             if(dummy2(1:1) .ne. 'r' .and. dummy2(1:1) .ne. 'R') then
               call str2real(dummy2, r_dummy)
               PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy ! re_scaled
             elseif( dummy2(1:1) .eq. 'r' .or. dummy2(1:1) .eq. 'R' ) then
               r_dummy = random()*2d0 - 1.d0
#ifdef MPI     
               call MPI_BCAST(r_dummy, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif         
               PPRAM%param_nrl(:,i) = PPRAM%param_nrl(:,i) * r_dummy + r_dummy
             endif
             if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
               param_const_nrl(4,:,i) = 1d0
               param_const_nrl(5,:,i) = PPRAM%param_nrl(:,i)
!            elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
!              param_const_nrl(4,:,i) = 0d0
             endif
           endif

         else
           write(message,'(A)') "  !!WARN!! Wrong syntax in PFILE, check your PFILE."    ; write_msgi
           write(message,'(A)') "           You have set SK_SCALE_MODE <= 10, but more than 3 arguments" ; write_msgi
           write(message,'(A)') "           has been asigned." ; write_msgi
           kill_job
         endif

       else
         write(message,'(A)') "  !!WARN!! Wrong syntax in PFILE, check your PFILE."   ; write_msgi
         write(message,'(A)') "           You have set SK_SCALE_MODE <= 10, but more than 3 arguments" ; write_msgi
         write(message,'(A)') "           has been asigned." ; write_msgi
         kill_job
       endif

       param_name = trim(PPRAM%param_name(i))
       if(param_name(1:2) .ne. 'e_' .and. param_name(1:2) .ne. 'ss' .and. param_name(1:2) .ne. 'sp' .and. &
          param_name(1:2) .ne. 'pp' .and. param_name(1:2) .ne. 'pd' .and. param_name(1:2) .ne. 'dp' .and. &
          param_name(1:2) .ne. 'dd' .and. param_name(1:2) .ne. 's_' .and. param_name(1:2) .ne. 'os' .and. &
          param_name(1:2) .ne. 'o_' ) then

         PPRAM%param_nsub(i) = 1
         if(param_name(1:4) .eq. 'lons' ) then
           param_const_nrl(2,1,i) = 2d0 ! upper bound
           param_const_nrl(3,1,i) = 0d0 ! lower bound
           param_const_nrl(2,2:4,i) = 0d0 ! upper bound
           param_const_nrl(3,2:4,i) = 0d0 ! lower bound
         endif
       else
         PPRAM%param_nsub(i) = 4  ! number of sub parameters (a, b, c, d for NRL hopping, alpha, beta, gamma, xi for NRL onsite)
         if(param_name(1:2) .ne. 'e_') then  ! set bound for parameter 'd' in expenential decay function ( exp(-d**2 * R)
           param_const_nrl(2,4,i) = 2d0 ! upper bound
           param_const_nrl(3,4,i) = 0d0 ! lower bound
         endif
       endif

     else

       if(i_dummy .eq. 1) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PPRAM%param_name(i),PPRAM%param(i)
         else
           read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i),PPRAM%param(i)
         endif
         call set_scale_default(param_const(3,i), PPRAM%param_name(i)) 
       elseif( i_dummy .eq. 2) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PPRAM%param_name(i),PPRAM%param(i),dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i),PPRAM%param(i),dummy
         endif
         call set_scale_default(param_const(3,i), PPRAM%param_name(i))
         dummy = trim(dummy)
         if(.not. flag_number(dummy) ) then
           if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
             param_const(4,i) = 1d0
             param_const(5,i) = PPRAM%param(i)
           elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'random '
!            param_const(4,i) = 0d0
             r_dummy = random()*2d0 - 1.d0
#ifdef MPI
             call MPI_BCAST(r_dummy, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif
             PPRAM%param(i) = PPRAM%param(i) * r_dummy + r_dummy ! re_scaled
           endif
         elseif( flag_number(dummy) ) then
           call str2real(dummy, r_dummy)
           PPRAM%param(i) = PPRAM%param(i) * r_dummy ! re_scaled
         endif
       elseif( i_dummy .eq. 3) then
         if(.not. flag_index) then
           read(inputline,*,iostat=i_continue) PPRAM%param_name(i),PPRAM%param(i),dummy2, dummy
         else
           read(inputline,*,iostat=i_continue) i_index, PPRAM%param_name(i),PPRAM%param(i),dummy2, dummy
         endif
         call set_scale_default(param_const(3,i), PPRAM%param_name(i))
         dummy2= trim(dummy2)
         if(.not.flag_number(dummy2) .and. (dummy2(1:1) .ne. 'r' .and. dummy2(1:1) .ne. 'R')) then 
           write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msgi
           kill_job
         endif
         dummy = trim(dummy)
         if(flag_number(dummy)) then
           write(message,'(A)') "  !!WARN!! wrong syntax in PFILE, check your PFILE." ; write_msgi
           kill_job
         endif
         if(dummy2(1:1) .ne. 'r' .and. dummy2(1:1) .ne. 'R') then
           call str2real(dummy2, r_dummy)
           PPRAM%param(i) = PPRAM%param(i) * r_dummy ! re_scaled
         elseif(dummy2(1:1) .eq. 'r' .or. dummy2(1:1) .eq. 'R') then
           r_dummy = random()*2d0 - 1.d0
#ifdef MPI
           call MPI_BCAST(r_dummy, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif
           PPRAM%param(i) = PPRAM%param(i) * r_dummy +r_dummy ! re_scaled
         endif
         if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
           param_const(4,i) = 1d0
           param_const(5,i) = PPRAM%param(i)
!        elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
!          param_const(4,i) = 0d0
         endif

       elseif( i_dummy .ge. 4) then
         write(message,'(A)') "  !!WARN!! Wrong syntax in PFILE, check your PFILE." ; write_msgi
         write(message,'(A)') "           You have set SK_SCALE_MODE <= 10, but more than 3 arguments" ; write_msgi
         write(message,'(A)') "           has been asigned." ; write_msgi
         kill_job
       endif
       param_name = trim(PPRAM%param_name(i))
       if(param_name(1:2) .ne. 'e_' .and. param_name(1:2) .ne. 'ss' .and. param_name(1:2) .ne. 'sp' .and. &
          param_name(1:2) .ne. 'pp' .and. param_name(1:2) .ne. 'pd' .and. param_name(1:2) .ne. 'dp' .and. &
          param_name(1:2) .ne. 'dd' .and. param_name(1:2) .ne. 's_' .and. param_name(1:2) .ne. 'os' .and. &
          param_name(1:2) .ne. 'o_' ) then

         PPRAM%param_nsub(i) = 1
       else
         PPRAM%param_nsub(i) = 1
       endif
     endif

   enddo
   close(pid_param)

   PPRAM%param_best = PPRAM%param  

   ! SET PARAMETER CONSTRAINT : from PFILE
   allocate( PPRAM%param_const(5,PPRAM%nparam) )
!  allocate( PPRAM%param_const_best(5,PPRAM%nparam) )
   !initialize
!   PPRAM%param_const(:,0) = 0d0 ! this value should be zero
    PPRAM%param_const(1,:) =param_const(1,1:PPRAM%nparam)  ! if gt 0 and it is "i", param is same as i-th parameter 
    PPRAM%param_const(2,:) =param_const(2,1:PPRAM%nparam)  ! default upper bound 
    PPRAM%param_const(3,:) =param_const(3,1:PPRAM%nparam)  ! default lower bound
    PPRAM%param_const(4,:) =param_const(4,1:PPRAM%nparam)  ! if set to 1; fix 
    PPRAM%param_const(5,:) =param_const(5,1:PPRAM%nparam)  ! if set to 1; fix and save the parameter as constant
!   PPRAM%param_const_best = PPRAM%param_const

   if(PPRAM%slater_koster_type .gt. 10) then
     allocate( PPRAM%param_const_nrl(5,4,PPRAM%nparam) )
     !initialize
!    PPRAM%param_const_nrl(:,:,0) = 0d0 ! this value should be zero
     PPRAM%param_const_nrl(1,:,:) =param_const_nrl(1,:,1:PPRAM%nparam)  ! if gt 0 and it is "i", param is same as i-th parameter 
     PPRAM%param_const_nrl(2,:,:) =param_const_nrl(2,:,1:PPRAM%nparam)  ! default upper bound 
     PPRAM%param_const_nrl(3,:,:) =param_const_nrl(3,:,1:PPRAM%nparam)  ! default lower bound
     PPRAM%param_const_nrl(4,:,:) =param_const_nrl(4,:,1:PPRAM%nparam)  ! if set to 1; fix 
     PPRAM%param_const_nrl(5,:,:) =param_const_nrl(5,:,1:PPRAM%nparam)  ! if set to 1; fix and save the parameter as constant
   endif
   
   ! SET PARAMETER CONSTRAINT : from IFILE (INCAR-TB)
   call set_param_const(PPRAM) ! set parameter constraints (SET CONSTRAINT TBPARAM)

   ! REPORT PARAM
   call report_param(PINPT, PPRAM)

   ! Check parameter type and condition
   if(.not. PPRAM%flag_slater_koster .and. PPRAM%flag_use_overlap) then
     write(message,'(A)')'    !WARN! Construction of overlap matrix is only available within Slakter-Koster method turned on.'              ; write_msgi
     write(message,'(A)')'           Set "USE_OVERLAP .FALSE." in your PARAM_FIT.dat file or set "IS_SK .TRUE." and prepare '    ; write_msgi
     write(message,'(A)')'           proper parameter set for overlap integrals, for example, o_pps_1_CC, o_sps_1_CC, and etc., to proceed' ; write_msgi
     write(message,'(A)')'           Exit program...' ; write_msgi
     kill_job
   endif

   write(message,*)'#- END READING PARAMETER FILE ---------------------'  ; write_msgi
   write(message,*)' '  ; write_msgi

return
endsubroutine

subroutine set_scale_default( param_const, param_name)
   implicit none
   character*40  param_name
   real*8        param_const

   if(param_name(1:2) .eq. 's_') then
     param_const = 0.001 ! set the lower bound of scaling factor
   elseif(param_name(1:3) .eq. 'os_') then
     param_const = 0.001 ! set the lower bound of scaling factor
   endif

return
endsubroutine

subroutine set_param_const(PPRAM)
   use parameters, only : params
   use mpi_setup
   use print_io
   implicit none
   integer*4     i, ii, i_a, i_b
   integer*4     j
   character*40  dummy
   type(params)  :: PPRAM
   integer*4     mpierr

!PPRAM%param_const(i,:) i=1 -> is same as
!                       i=2 -> is lower than (.le.) : set maximum bound  ! functionality is not available yet
!                       i=3 -> is lower than (.ge.) : set minimum bound  ! functionality is not available yet
!                       i=4 -> is fixed : not to be fitted, just stay    ! its original value will be stored in PPRAM%param_const(i=5,:)

!  allocate( PPRAM%param_const(5,PPRAM%nparam) )
!  PPRAM%param_const(1,:) = 0d0  ! same as 
!  PPRAM%param_const(2,:) = 20d0 ! upper bound
!  PPRAM%param_const(3,:) =-20d0 ! lower bound ; if start with s_ (scale factor) then >= 0.001 (default)
!  PPRAM%param_const(4,:) = 0d0  ! fixed 1:true 0:no
!  PPRAM%param_const(5,:) = 0d0  ! fixed value

!  NOTE: same rule applies for PPRAM%param_const_nrl(1:5,:,:) as PPRAM%param_const(1:5,:)

   if(.not. PPRAM%flag_set_param_const) return

   if(PPRAM%slater_koster_type .le. 10) then
     do i = 1, PPRAM%nparam_const
       if( trim(PPRAM%c_const(2,i)) .eq. '=' ) then
         dummy = trim(PPRAM%c_const(3,i))
         if( dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f')then  ! check if fixed
           do ii = 1, PPRAM%nparam
             if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then 
               i_a = ii
               PPRAM%param_const(4,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
               PPRAM%param_const(5,i_a) = PPRAM%param(i_a) ! save its value to PPRAM%param_const(5,i_a)
             endif
           enddo
         else
           i_a = 0
           do ii = 1, PPRAM%nparam
             if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) i_a = ii
           enddo
           do ii = 1, PPRAM%nparam
             if( trim(PPRAM%c_const(3,i)) .eq. trim(PPRAM%param_name(ii)) .and. i_a .ne. 0 ) then 
               i_b = ii
               PPRAM%param_const(1,i_a) = real(i_b)
             endif
           enddo
         endif
         cycle
       elseif( trim(PPRAM%c_const(2,i)) .eq. '<=' ) then
         do ii = 1, PPRAM%nparam
           if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then 
             i_a = ii
             call str2real( trim(PPRAM%c_const(3,i)), PPRAM%param_const(2,i_a) )
           endif
         enddo
         cycle
       elseif( trim(PPRAM%c_const(2,i)) .eq. '>=' ) then
         do ii = 1, PPRAM%nparam
           if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then 
             i_a = ii
             call str2real( trim(PPRAM%c_const(3,i)), PPRAM%param_const(3,i_a) )
           endif
         enddo
         cycle
       elseif( trim(PPRAM%c_const(2,i)) .eq. '==' ) then
         do ii = 1, PPRAM%nparam
           if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then 
             i_a = ii
             PPRAM%param_const(4,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
             PPRAM%param_const(5,i_a) = PPRAM%param(i_a) ! save its value to PPRAM%param_const(5,i_a)
           endif
         enddo
         cycle
       else
         write(message,'(A)')'  !WARNING! parameter constraint is not properly defined. Please check again. Exit...' ; write_msgi
         kill_job
       endif
     enddo
  
   elseif(PPRAM%slater_koster_type .gt. 10) then
     do i = 1, PPRAM%nparam_const
       if( trim(PPRAM%c_const(2,i)) .eq. '=' ) then
         dummy = trim(PPRAM%c_const(3,i))
         if( dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f')then  ! check if fixed
           do ii = 1, PPRAM%nparam
             if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then
               i_a = ii
               PPRAM%param_const_nrl(4,:,i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
               PPRAM%param_const_nrl(5,:,i_a) = PPRAM%param_nrl(:,i_a) ! save its value to PPRAM%param_const(5,i_a)
             endif
           enddo
         else
           i_a = 0
           do ii = 1, PPRAM%nparam
             if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) i_a = ii
           enddo
           do ii = 1, PPRAM%nparam
             if( trim(PPRAM%c_const(3,i)) .eq. trim(PPRAM%param_name(ii)) .and. i_a .ne. 0 ) then
               i_b = ii
               PPRAM%param_const_nrl(1,:,i_a) = real(i_b)
             endif
           enddo
         endif
         cycle
       elseif( trim(PPRAM%c_const(2,i)) .eq. '<=' ) then  ! set upper bound
         do ii = 1, PPRAM%nparam
           if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then
             i_a = ii
             do j = 1, PPRAM%param_nsub_max
               call str2real( trim(PPRAM%c_const(3,i)), PPRAM%param_const_nrl(2,j, i_a) )
             enddo
           endif
         enddo
         cycle
       elseif( trim(PPRAM%c_const(2,i)) .eq. '>=' ) then
         do ii = 1, PPRAM%nparam
           if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then
             i_a = ii
             do j = 1, PPRAM%param_nsub_max
               call str2real( trim(PPRAM%c_const(3,i)), PPRAM%param_const_nrl(3,j, i_a) )
             enddo
           endif
         enddo
         cycle
       elseif( trim(PPRAM%c_const(2,i)) .eq. '==' ) then
         do ii = 1, PPRAM%nparam
           if( trim(PPRAM%c_const(1,i)) .eq. trim(PPRAM%param_name(ii)) ) then
             i_a = ii
             PPRAM%param_const_nrl(4,:, i_a) = 1d0 !turn on fix the parameter to be remained as the initial guess
             PPRAM%param_const_nrl(5,:, i_a) = PPRAM%param_nrl(:,i_a) ! save its value to PPRAM%param_const(5,i_a)
           endif
         enddo
         cycle
       else
         write(message,'(A)')'  !WARNING! parameter constraint is not properly defined. Please check again. Exit...' ; write_msgi
         kill_job
       endif
     enddo
   endif

return
endsubroutine

! update with constraint values which is redirect with (1,i) indexing stored in param_const
subroutine update_param(PPRAM)
   use parameters, only : params
   implicit none
   type(params) :: PPRAM
   integer*4       i

   do i = 1, PPRAM%nparam
     if( nint(PPRAM%param_const(1,i)) .ge. 1) then ! 
       PPRAM%param(i) = PPRAM%param( nint(PPRAM%param_const(1,i)) )
     endif
   enddo

return
endsubroutine
subroutine update_param_nrl(PPRAM)
   use parameters, only : params
   implicit none
   type(params) :: PPRAM
   integer*4       i, j, k

   do j = 1, PPRAM%nparam
     k = nint(PPRAM%param_const_nrl(1,1,j))
     if( k .ge. 1) then ! 
       do i = 1,PPRAM%param_nsub(j)
         PPRAM%param_nrl(i,j) = PPRAM%param_nrl(i, k)
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

subroutine check_options(inputline, PPRAM, flag_skip)
   use parameters, only : params
   implicit none
   type(params) :: PPRAM
   logical         flag_skip  
   character*132   inputline
   character*40    dummy, dummy1

   read(inputline, *) dummy

   ! read options if it is predefined, otherwise it will be recognized as 'parameters'
   if( trim(dummy) .eq. 'PRINT_INDEX' .or. trim(dummy) .eq. 'print_index') then
     read(inputline, *) dummy1, PPRAM%flag_pfile_index
     flag_skip = .true.
   elseif( trim(dummy) .eq. 'USE_OVERLAP' .or. trim(dummy) .eq. 'use_overlap' ) then
     read(inputline, *) dummy1, PPRAM%flag_use_overlap
     flag_skip = .true.
   elseif( trim(dummy) .eq. 'USE_SKPARAM' .or. trim(dummy) .eq. 'IS_SK' .or. &
           trim(dummy) .eq. 'SLATER_KOSTER' ) then
     read(inputline, *) dummy1, PPRAM%flag_slater_koster 
     flag_skip = .true.
   elseif( trim(dummy) .eq. 'SK_SCALE_MODE' .or. trim(dummy) .eq. 'SCALE_MODE' .or. &
           trim(dummy) .eq. 'sk_scale_mode' .or. trim(dummy) .eq. 'scale_mode') then
     read(inputline, *) dummy1, PPRAM%slater_koster_type 
     flag_skip = .true. 
     if(PPRAM%slater_koster_type .gt. 10) then
       PPRAM%param_nsub_max = 4
     else
       PPRAM%param_nsub_max = 1
     endif
   elseif( trim(dummy) .eq. 'L_BROADEN' .or. trim(dummy) .eq. 'l_broaden' ) then
     read(inputline, *) dummy1, PPRAM%l_broaden
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

subroutine report_param(PINPT, PPRAM)
   use parameters, only : incar, params
   use print_io
   use mpi_setup
   implicit none
   type(incar)   :: PINPT
   type(params)  :: PPRAM
   integer*4        mpierr
   integer*4        i

   if(PPRAM%flag_slater_koster) then
     write(message,'(A)'   )'TYP_PARAM: (parameter type) Slater-Koster'  ; write_msgi
     write(message,'(A,I0)')'TYP_SCALE: (parameter scaling function) mode = ',PPRAM%slater_koster_type; write_msgi

     select case(PPRAM%slater_koster_type)
       case(1)
         write(message,'(A)'   )'           => f_s = Exp( (R0 - D)/(dda_s*R0) ), Ref: PRB 85.195458 (2012)' ; write_msgi
       case(2)
         write(message,'(A)'   )'           => f_s = Exp( -(D/R0)**dda_s ), Ref: PRB 92.205108 (2015)' ; write_msgi
       case(3)
         write(message,'(A)'   )'           => f_s = (R0/D)**(dda_s), Ref: PRB 51.16772 (1995)' ; write_msgi
       case(4)
         write(message,'(A)'   )'           => f_s = Exp( -abs(dda_s) * (D - R0) ), PRB 93.241407 (2016)' ; write_msgi
       case(11)
         write(message,'(A)'   )'           => e_onsite = e(1) + e(2) * (rho_at**(2/3)) + e(3) * (rho_at**(4/3)) + e(4) * (rho_at**(2))'; write_msgi
         write(message,'(A)'   )'           => f_s      = (dda(1)+dda(2)*D+dda(3)*(D**2)) * Exp(-(dda(4)**2) * D) * f_cut'; write_msgi
         write(message,'(A)'   )'              rho_at   = sum(Exp(-(l_onsite(:)**2) * R_nn(:) ) * f_cut(:))'  ; write_msgi
         write(message,'(A)'   )'              f_cut    = (1 + Exp( (D - R0)/l_broaden + 5 ) )**(-1) ' ; write_msgi
     end select

     write(message,'(A)'   )'              D        = distance to neighbor atom'; write_msgi
     write(message,'(A)'   )'              R0       = reference distance' ;write_msgi
     write(message,'(A)'   )'           ==> Scaled parameter = PARAM * f_s' ; write_msgi
   else
     write(message,'(A)')'TYP_PARAM: (parameter type) user defined'  ; write_msgi
   endif

   write(message,'(A,I8)')'  N_PARAM:',PPRAM%nparam ; write_msgi
   if(PPRAM%slater_koster_type .gt. 10) then
     write(message,'(A)') '         : NRL TB scheme is applied in parameterization' ; write_msgi
     write(message,'(A,F9.4)') '           => L_BROADEN (cutoff function) = ', PPRAM%l_broaden ; write_msgi
   endif

   do i=1,PPRAM%nparam
     if(PPRAM%slater_koster_type .gt. 10) then
       write(message,'(A,2x,A14,1x,*(F10.5))')'  C_PARAM:',PPRAM%param_name(i),PPRAM%param_nrl(1:PPRAM%param_nsub(i),i) ; write_msgi_file
     else
       write(message,'(A,2x,A14,1x,F10.5)')'  C_PARAM:',PPRAM%param_name(i),PPRAM%param(i) ; write_msgi_file
     endif
   enddo


   return
endsubroutine

