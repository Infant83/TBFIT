#include "alias.inc"
subroutine parse(PINPT)
   use parameters, only: incar
   use mpi_setup
   implicit none
   character*20        option, value
   character*20        dummy
   integer*4           narg, iarg
   logical,external :: flag_number
   logical             flag_logical, flag
   type(incar)      :: PINPT
   
   PINPT%flag_parse       = .false.
   PINPT%flag_tbfit_parse = .false.
   PINPT%flag_tbfit_test  = .false.
   PINPT%flag_kfile_parse = .false.
   PINPT%flag_miter_parse = .false.
   PINPT%flag_mxfit_parse = .false.
   PINPT%flag_lorbit_parse= .false.
   PINPT%flag_ldos_parse  = .false.
   PINPT%flag_inputcard_fname_parse = .false.
   narg = iargc()
   

  if_main  write(6,*)' '
  if_main  write(6,*)'---- READING INPUT TAG FROM: COMMAND-LINE ARGUMENTS'

 arg:do iarg = 1, narg
       call getarg(iarg, option)
       if(.not. flag_number(trim(option))) then
         if(trim(option) .eq. '-h') then
           if_main call help()

         elseif(trim(option) .eq. '-input') then
           PINPT%flag_inputcard_fname_parse = .true.
           call getarg(iarg+1, value)
           if_main write(6,'(2A)')'  IN_FILE: ',trim(value)
           inquire(file=trim(value),exist=flag)
           if(.not. flag) then
             if_main write(6,'(3A)')'    !WARN! INPUT_FILE:"',trim(value),'" does not exist! Check again... Exit program.'
             kill_job
           elseif(flag) then
             PINPT%ifilenm = trim(value)
           endif

         elseif(trim(option) .eq. '-fit' .or. trim(option) .eq. '-f') then
           PINPT%flag_tbfit_parse  = .true.
           PINPT%flag_tbfit_parse_ = .TRUE.
           if_main write(6,'(A)')'  L_TBFIT:   .TRUE. (enforce by -fit option)'

         elseif(trim(option) .eq. '-nofit' .or. trim(option) .eq. '-n') then
           PINPT%flag_tbfit_parse  = .true.
           PINPT%flag_tbfit_parse_ = .FALSE.
           if_main write(6,'(A)')'  L_TBFIT:   .FALSE. (enforce by -nofit option)'

         elseif(trim(option) .eq. '-test' .or. trim(option) .eq. '-t') then         
           PINPT%flag_tbfit_parse  = .true.
           PINPT%flag_tbfit_test   = .true.
           if_main write(6,'(A)')'  !WARN! Program will run with TEST mode (-test option is detected).'
           if_main write(6,'(A)')'         After calling test() routine, program will stop immediately.'

         elseif(trim(option) .eq. '-param' .or. trim(option) .eq. '-p') then
           PINPT%flag_parse  = .true.
           PINPT%flag_pfile = .true.
           call getarg(iarg+1, PINPT%pfilenm)
           if_main write(6,'(A,A)')' PARA_FNM:  ',trim(PINPT%pfilenm)

         elseif(trim(option) .eq. '-kpooint' .or. trim(option) .eq. '-kp' .or. trim(option) .eq. '-k') then
           PINPT%flag_parse = .true.
           PINPT%flag_kfile_parse = .true.
           call getarg(iarg+1, PINPT%kfilenm)
           if_main write(6,'(A,A)')' KPTS_FNM:  ',trim(PINPT%kfilenm)

         elseif(trim(option) .eq. '-miter' .or. trim(option) .eq. '-m') then
           PINPT%flag_parse = .true.
           PINPT%flag_miter_parse = .true.
           call getarg(iarg+1, value)
           read(value, *) PINPT%miter
           if_main write(6,'(A,I8)')' MAX_ITER:  ',PINPT%miter

         elseif(trim(option) .eq. '-mxfit' .or. trim(option) .eq. '-mf') then
           PINPT%flag_parse = .true.
           PINPT%flag_mxfit_parse = .true.
           call getarg(iarg+1, value)
           read(value, *) PINPT%mxfit
           if_main write(6,'(A,I8)')'  MAX_FIT:  ',PINPT%mxfit
    
         elseif(trim(option) .eq. '-lorbit') then
           PINPT%flag_parse = .true.
           PINPT%flag_lorbit_parse = .true.
           PINPT%flag_get_orbital = .true.
           PINPT%flag_print_mag = .false.
           call getarg(iarg+1, value)
           if(trim(value) .eq. 're' ) then
             PINPT%flag_print_mag = .TRUE.
             read(value, *) PINPT%axis_print_mag
             if_main write(6,'(2A)')'  L_ORBIT: .TRUE. | print out real part of wavefnc.: ', PINPT%axis_print_mag
           elseif(trim(value) .eq. 'im' ) then
             PINPT%flag_print_mag = .TRUE.
             read(value, *) PINPT%axis_print_mag
             if_main write(6,'(2A)')'  L_ORBIT: .TRUE. | print out imag part of wavefnc.: ', PINPT%axis_print_mag
           elseif(trim(value) .eq. 'mx' .or. trim(value) .eq. 'my' .or. trim(value) .eq. 'mz') then
             PINPT%flag_print_mag = .TRUE.
             read(value, *) PINPT%axis_print_mag
             if_main write(6,'(2A)')'  L_ORBIT: .TRUE. | print out magnetization <sigma>: ', PINPT%axis_print_mag
           elseif(trim(value) .eq. 'no' ) then
             PINPT%flag_lorbit_parse=.true.
             PINPT%flag_get_orbital =.false.
             if_main write(6,'( A)')'  L_ORBIT: .FALSE.'
           elseif(trim(value) .eq. 'yes' ) then
             PINPT%flag_lorbit_parse=.true.
             PINPT%flag_get_orbital =.true. 
             if_main write(6,'( A)')'  L_ORBIT: .TRUE. | print out projected orbital weight'
           else
             PINPT%flag_lorbit_parse=.true.
             PINPT%flag_get_orbital =.true. 
             if_main write(6,'( A)')'  L_ORBIT: .TRUE. | print out projected orbital weight'
           endif

         elseif(trim(option) .eq. '-ldos') then
           PINPT%flag_parse = .true.
           PINPT%flag_lorbit_parse = .true.
           PINPT%flag_get_orbital = .true.
           PINPT%flag_print_mag = .false.
           PINPT%flag_ldos_parse = .true.

           if(iarg + 1 .le. narg) then 
             call getarg(iarg+1, value)
             read(value,*)dummy
             call str2logical(dummy,flag_logical, flag)
             if(flag_logical) then 
               PINPT%flag_print_ldos = flag
             else
               PINPT%flag_print_ldos = .true.
             endif
           elseif(iarg + 1 .gt. narg) then
              PINPT%flag_print_ldos = .true.
           endif
    
           if(PINPT%flag_print_ldos) then
             if_main write(6,'(A)')'  L_LDOS: .TRUE. | print out atom projected orbital weight'
           elseif(.not. PINPT%flag_print_ldos) then
             if_main write(6,'(A)')'  L_LDOS: .FALSE.'
           endif

         endif

       endif

     enddo arg

   if_main  write(6,*)'---- END READING INPUT FROM: COMMAND-LINE ARGUMENTS ---------------------'
   if_main  write(6,*)' '

   return
endsubroutine

subroutine help()
   implicit none

   write(6,'(A)')"          **** TBFIT: COMMAND LINE ARGUMENTS ***"
   write(6,'(A)')" "
   write(6,'(A)')" "
   write(6,'(A)')" ## POSSIBLE OPTIONS ##"
   write(6,'(A)')"                                 "
   write(6,'(A)')"   -h         : print help pages"
   write(6,'(A)')"                                 "
   write(6,'(A)')"   -input INP : enforce to read INP file instead of INCAR-TB for the input card file"
   write(6,'(A)')"   -fit       : enforce to run with 'fitting' mode"
   write(6,'(A)')"   -nofit     : enforce not to run with 'fitting' even if the TBFIT tag is set to .TRUE."
   write(6,'(A)')"   -miter MIT : enforce to set maximum number of iteration for LMDIF to MIT in prior to the MITER tag"
   write(6,'(A)')"   -param  PF : enforce to read parameter file 'PF' in prior to the PFILE tag"
   write(6,'(A)')"   -kpoint KF : enforce to read k-points  file 'KF' in prior to the KFILE tag"
   write(6,'(A)')"   -lorbit    : enforce    to print orbital information"
   write(6,'(A)')"           yes: enforce    to print orbital information"
   write(6,'(A)')"           no : enforce not to print orbital information"
   write(6,'(A)')"           mx : enforce     to print magnetization mx " 
   write(6,'(A)')"           my : enforce     to print magnetization mz " 
   write(6,'(A)')"           mz : enforce     to print magnetization mz " 
   write(6,'(A)')"   -ldos      : enforce    to print orbital information for each atom in separate file"
   stop

   return
endsubroutine
