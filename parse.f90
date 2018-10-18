#include "alias.inc"
subroutine parse(PINPT)
   use parameters, only: incar
   use mpi_setup
   implicit none
   character*20        option, value
   integer*4           narg, iarg
   logical,external :: flag_number
   type(incar)      :: PINPT
   
   PINPT%flag_tbfit_parse = .false.
   PINPT%flag_tbfit_test  = .false.
   narg = iargc()
   


 arg:do iarg = 1, narg
       call getarg(iarg, option)
       if(.not. flag_number(trim(option))) then
         if(trim(option) .eq. '-h') then
           if_main call help()

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
         endif

       endif

     enddo arg

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
   write(6,'(A)')"   -fit       : enforce to run with 'fitting' mode"
   write(6,'(A)')"   -nofit     : enforce not to run with 'fitting' even if the TBFIT tag is set to .TRUE."

   stop

   return
endsubroutine
