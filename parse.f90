#include "alias.inc"

subroutine parse_very_init_py(PINPT, nsystem, ifilenm_dummy )
   use parameters, only: incar, max_dummy
   use mpi_setup
   implicit none
   type(incar)       :: PINPT
   integer*4            narg, iarg, i, nsystem
   character*132        option, value
   character*132        ifilenm_dummy(nsystem)
   character*132        dummy
   logical, external :: flag_number

   PINPT%fnamelog = 'TBFIT.out' ! default
   PINPT%nsystem = nsystem
   PINPT%flag_python_module = .TRUE.

   allocate(PINPT%ifilenm(PINPT%nsystem))
   allocate(PINPT%title(PINPT%nsystem))

   if(PINPT%nsystem .eq. 1) then  
     PINPT%ifilenm(1) = trim(ifilenm_dummy(1))
     PINPT%title(1) = ''
   elseif(PINPT%nsystem .ge. 2) then
     do i = 1, PINPT%nsystem
       PINPT%ifilenm(i) = trim(ifilenm_dummy(i))
       write(PINPT%title(i),'(A,I0)') '.',i
     enddo
   endif


   return
endsubroutine

subroutine parse_very_init(PINPT)
   use parameters, only: incar, max_dummy
   use mpi_setup
   implicit none
   type(incar)       :: PINPT
   integer*4            narg, iarg, i
   character*132        option, value
   character*132        ifilenm_dummy(max_dummy)
   character*132        dummy
   logical, external :: flag_number

   PINPT%fnamelog = 'TBFIT.out' ! default
   PINPT%nsystem = 0 
   narg = iargc()
   PINPT%flag_python_module = .FALSE.

   do iarg = 1, narg
     call getarg(iarg, option)
     if(.not. flag_number(trim(option))) then
       if(trim(option) .eq. '-log' .or. trim(option) .eq. '-o' .or. trim(option) .eq. '-out') then
         call getarg(iarg+1, PINPT%fnamelog) ! set output file name
       elseif(trim(option) .eq. '-input' .or. trim(option) .eq. '-i') then
         PINPT%nsystem = PINPT%nsystem + 1
         call getarg(iarg+1, dummy) ! set input setup file name
         ifilenm_dummy(PINPT%nsystem) = dummy
       endif
     endif
   enddo

   if(PINPT%nsystem .eq. 0) then
     allocate(PINPT%ifilenm(1))    ! default
     PINPT%ifilenm(1) = 'INCAR-TB' ! default
     PINPT%nsystem = 1             ! default
     allocate(PINPT%title(1))
     PINPT%title(1) = '' ! no title by default
   elseif(PINPT%nsystem .eq. 1) then
     allocate(PINPT%ifilenm(PINPT%nsystem))
     allocate(PINPT%title(PINPT%nsystem))
     do i = 1, PINPT%nsystem
       PINPT%ifilenm(i) = trim(ifilenm_dummy(i))
       PINPT%title(i) = '' ! no title by default
     enddo
   elseif(PINPT%nsystem .ge. 2) then
     allocate(PINPT%ifilenm(PINPT%nsystem))
     allocate(PINPT%title(PINPT%nsystem))
     do i = 1, PINPT%nsystem
       PINPT%ifilenm(i) = trim(ifilenm_dummy(i))
       write(PINPT%title(i),'(A,I0)') '.',i  ! numbering by default
     enddo
   endif

   return
endsubroutine

subroutine parse(PINPT)
   use parameters, only: incar
   use mpi_setup
   use print_io
   use version
   implicit none
   character*20        option, value
   character*20        dummy
   integer*4           narg, iarg, mpierr
   integer*4           i, j
   integer*4           iverbose_
   logical,external :: flag_number
   logical             flag_logical, flag, flag_ifile_exist
   logical             flag_exist
   type(incar)      :: PINPT

   PINPT%flag_parse       = .false.
   PINPT%flag_plot        = .false.
   PINPT%flag_tbfit_parse = .false.
   PINPT%flag_ls_type_parse = .false.
   PINPT%flag_tbfit       = .false.
   PINPT%flag_tbfit_finish= .false. ! only activate when exit from fitting routine 
   PINPT%flag_tbfit_test  = .false.
   PINPT%flag_kfile_parse = .false.
   PINPT%flag_miter_parse = .false.
   PINPT%flag_mxfit_parse = .false.
   PINPT%flag_lorbit_parse= .false.
   PINPT%flag_filenm_gnuplot_parse = .false.
   PINPT%flag_inputcard_fname_parse = .false. ! deprecated
   PINPT%flag_ndiv_line_parse = .false.
   PINPT%flag_ndiv_grid_parse = .false.
   PINPT%flag_print_only_target = .false.
   PINPT%flag_pso_verbose_parse = .false.
   PINPT%flag_fit_orbital_parse = .false.
   PINPT%flag_reduce_overlap_parse = .false.
   PINPT%flag_reduce_hopping_parse = .false.
   PINPT%reduce_overlap_parse = 1.0d0
   PINPT%reduce_hopping_parse = 1.0d0

   narg = iargc()
   
   if(PINPT%flag_python_module) then
     iverbose = 2 ! 1: full, 2: no
     iverbose_= 2
     print_mode = 99 ! default verbosity 
   else
     iverbose = 1 ! 1: full, 2: no
     iverbose_= 1
     print_mode = 3  ! default verbosity 
   endif
   write(message,*)' '; write_msg
   write(message,*)'#--- READING INPUT TAG FROM: COMMAND-LINE ARGUMENTS ---------------------'  ; write_msg

 arg:do iarg = 1, narg
       call getarg(iarg, option)
       if(.not. flag_number(trim(option))) then
         if(trim(option) .eq. '-h') then
           if_main call help()
           kill_job
         elseif(trim(option) .eq. '-v') then
           call getarg(iarg+1, value )
           read(value, *) iverbose_ 
           if(iverbose_ .eq. 1) then
             write(message,'(A)')'  VERBOSE: full'
           elseif(iverbose_ .eq. 2) then
             write(message,'(A)')'  VERBOSE: none'
           elseif(iverbose_ .eq. -2) then
             write(message,'(A)')'  VERBOSE: write to file only'
           endif
    
         elseif(trim(option) .eq. '-fit' .or. trim(option) .eq. '-f') then
           PINPT%flag_tbfit_parse  = .true.
           PINPT%flag_tbfit = .true.
           write(message,'(A)')'  L_TBFIT:   .TRUE. (enforce by -fit option)'  ; write_msg

         elseif(trim(option) .eq. '-nofit' .or. trim(option) .eq. '-n') then
           PINPT%flag_tbfit_parse  = .true.
           PINPT%flag_tbfit = .false.
           write(message,'(A)')'  L_TBFIT:   .FALSE. (enforce by -nofit option)'  ; write_msg

         elseif(trim(option) .eq. '-np') then
           PINPT%flag_tbfit_parse  = .true.
           PINPT%flag_tbfit = .false.
           write(message,'(A)')'  L_TBFIT:   .FALSE. (enforce by -np option)'  ; write_msg
           PINPT%flag_plot         = .TRUE.
           write(message,'(A)')'   L_PLOT:   .TRUE. (gnuplot gnuBAND-TB.gpi)'  ; write_msg

         elseif(trim(option) .eq. '-red_ovl') then
           PINPT%flag_reduce_overlap_parse  = .true.
           call getarg(iarg+1, value )
           read(value, *) PINPT%reduce_overlap_parse
           write(message,'(A,F9.4)')'  RED_OVL:  ', PINPT%reduce_overlap_parse  ; write_msg

         elseif(trim(option) .eq. '-red_hop') then
           PINPT%flag_reduce_hopping_parse  = .true.
           call getarg(iarg+1, value )
           read(value, *) PINPT%reduce_hopping_parse
           write(message,'(A,F9.4)')'  RED_HOP:  ', PINPT%reduce_hopping_parse  ; write_msg

         elseif(trim(option) .eq. '-ls' .or. trim(option) .eq. '-lstype') then
           PINPT%flag_ls_type_parse = .true.
           call getarg(iarg+1, PINPT%ls_type)
           write(message,'(2A)')'  LS_TYPE:  ',trim(PINPT%ls_type) ; write_msg

         elseif(trim(option) .eq. '-print_only' .or. trim(option) .eq. '-po' .or. trim(option) .eq. '-print_weight') then
           PINPT%flag_print_only_target = .true.
           PINPT%flag_tbfit_parse  = .true.
           PINPT%flag_tbfit = .true.
           write(message,'(A)')'  L_TBFIT:   .TRUE. (enforce by -print_only or equivalent option)'  ; write_msg
           write(message,'(A)')'        AND  PRINT_ONLY_TARGET requested ...                     '  ; write_msg

         elseif(trim(option) .eq. '-plot' .or. trim(option) .eq. '-pl') then
           PINPT%flag_plot         = .TRUE.
           write(message,'(A)')'   L_PLOT:   .TRUE. (gnuplot gnuBAND-TB.gpi)'  ; write_msg

         elseif(trim(option) .eq. '-gnuplot') then
           PINPT%flag_filenm_gnuplot_parse = .true.
           call getarg(iarg+1, PINPT%filenm_gnuplot_parse)
           PINPT%flag_filenm_gnuplot_parse = .true.
           write(message,'(2A)')'   L_PLOT: GNUPLOT FILE NAME = ', trim(PINPT%filenm_gnuplot_parse); write_msg

         elseif(trim(option) .eq. '-test' .or. trim(option) .eq. '-t') then         
           PINPT%flag_tbfit_test   = .true.
           write(message,'(A)')'  !WARN! Program will run with TEST mode (-test option is detected).'  ; write_msg
           write(message,'(A)')'         After calling test() routine, program will stop immediately.'  ; write_msg

         elseif(trim(option) .eq. '-param' .or. trim(option) .eq. '-p') then
           PINPT%flag_pfile_parse = .true.
           call getarg(iarg+1, PINPT%pfilenm_parse)
           write(message,'(A,A)')' PARA_FNM:  ',trim(PINPT%pfilenm_parse)  ; write_msg
           inquire(file=trim(PINPT%pfilenm_parse),exist=flag_exist)
           if(.not. flag_exist) then
             write(message,'(A,A,A)')'    !WARN! Parameter file:',trim(PINPT%pfilenm_parse),' does not exist!! Exit...' ; write_msg
             kill_job
           endif

         elseif(trim(option) .eq. '-weight' .or. trim(option) .eq. '-w') then
           PINPT%flag_wfile_parse = .true.
           PINPT%flag_set_weight_from_file = .true.
           call getarg(iarg+1, PINPT%wfilenm_parse)
           write(message,'(A,A)')' WGHT_FNM:  ',trim(PINPT%wfilenm_parse)  ; write_msg
           inquire(file=trim(PINPT%wfilenm_parse),exist=flag_exist)
           if(.not. flag_exist) then
             write(message,'(A,A,A)')'    !WARN! Weight file:',trim(PINPT%wfilenm_parse),' does not exist!! Exit...' ; write_msg
             kill_job
           endif

         elseif(trim(option) .eq. '-kpoint' .or. trim(option) .eq. '-kp' .or. trim(option) .eq. '-k') then
           PINPT%flag_parse = .true.
           PINPT%flag_kfile_parse = .true.
           call getarg(iarg+1, PINPT%kfilenm_parse)
           write(message,'(A,A)')' KPTS_FNM:  ',trim(PINPT%kfilenm_parse)  ; write_msg

         elseif(trim(option) .eq. '-miter' .or. trim(option) .eq. '-m') then
           PINPT%flag_parse = .true.
           PINPT%flag_miter_parse = .true.
           call getarg(iarg+1, value)
           read(value, *) PINPT%miter

         elseif(trim(option) .eq. '-mxfit' .or. trim(option) .eq. '-mf') then
           PINPT%flag_parse = .true.
           PINPT%flag_mxfit_parse = .true.
           call getarg(iarg+1, value)
           read(value, *) PINPT%mxfit

         elseif(trim(option) .eq. '-pso_verbose') then
           call getarg(iarg+1, value)
           read(value, *) PINPT%pso_verbose
           PINPT%flag_pso_verbose_parse = .true.

         elseif(trim(option) .eq. '-lorbit') then
           PINPT%flag_parse = .true.
           PINPT%flag_lorbit_parse = .true.
           PINPT%flag_get_orbital = .true.
           PINPT%flag_print_orbital = .true.
           PINPT%flag_print_mag = .false.
           call getarg(iarg+1, value)
           if(trim(value) .eq. 're' ) then
             PINPT%flag_print_mag = .TRUE.
             read(value, *) PINPT%axis_print_mag
           elseif(trim(value) .eq. 'im' ) then
             PINPT%flag_print_mag = .TRUE.
             read(value, *) PINPT%axis_print_mag
           elseif(trim(value) .eq. 'mx' .or. trim(value) .eq. 'my' .or. trim(value) .eq. 'mz') then
             PINPT%flag_print_mag = .TRUE.
             read(value, *) PINPT%axis_print_mag
           elseif(trim(value) .eq. 'no' .or. trim(value) .eq. 'F' .or. trim(value) .eq. 'f' .or. &
                  trim(value) .eq. '.false.' .or. trim(value) .eq. '.FALSE.') then
             PINPT%flag_lorbit_parse=.true.
             PINPT%flag_get_orbital =.false.
           elseif(trim(value) .eq. 'yes' .or. trim(value) .eq. 'T' .or. trim(value) .eq. '.TRUE.' .or. &
                  trim(value) .eq. '.true.' .or. trim(value) .eq. 't') then
             PINPT%flag_lorbit_parse=.true.
             PINPT%flag_get_orbital =.true. 
           else
             PINPT%flag_lorbit_parse=.true.
             PINPT%flag_get_orbital =.true. 
           endif

         endif

       endif

     enddo arg

   ! check whether input file exists
   write(message,'(A,I0)')'  NSYSTEM: ', PINPT%nsystem ; write_msg
   do i = 1, PINPT%nsystem
     write(message,'(A,I0,2A)')'   SYSTEM  ',i, ' from : ',trim(PINPT%ifilenm(i)) ; write_msg
     inquire(file=trim(PINPT%ifilenm(i)),exist=flag_ifile_exist)
     if(.not. flag_ifile_exist) then
       write(message,'(A,A,A)')'    !WARN! Input file:',trim(PINPT%ifilenm(i)),' does not exist!! Exit...' ; write_msg
       kill_job
     endif
   enddo
   ! sanity check for the inputs if PINPT%nsystem > 1
   do i = 1, PINPT%nsystem
     do j = 1, PINPT%nsystem
       if(i .ne. j) then
         if( trim(PINPT%ifilenm(i)) .eq. trim(PINPT%ifilenm(j))) then
           write(message,'(A,I0,3A,I0,3A)')'    !WARN! ',i,'-th input file ',trim(PINPT%ifilenm(i)),' is same with ', &
                                                         j,'-th input file ',trim(PINPT%ifilenm(i)),' -> Exit...'; write_msg
           kill_job
         endif
       endif
     enddo
   enddo

   write(message,*)'#--- END READING INPUT FROM: COMMAND-LINE ARGUMENTS ---------------------' ; write_msg
   write(message,*)' ' ; write_msg

   iverbose = iverbose_

   return
endsubroutine

subroutine help()
   use mpi_setup
   use mykind
   implicit none
   integer(kind=sp) mpierr

   write(6,'(A)')"          **** TBFIT: COMMAND LINE ARGUMENTS ***"
   write(6,'(A)')" "
   write(6,'(A)')" "
   write(6,'(A)')" ## POSSIBLE OPTIONS ##"
   write(6,'(A)')"                                 "
   write(6,'(A)')"   -h                       : print this help page for command line arguments"
   write(6,'(A)')"                                               "
   write(6,'(A)')"   -log (or -o) FNAME       : output log will be written in FNAME. Default: TBFIT.out"
   write(6,'(A)')"   -input INP               : enforce to read INP file instead of INCAR-TB for the input card file"
   write(6,'(A)')"   -fit or -f               : enforce to run with 'fitting' mode"
   write(6,'(A)')"   -nofit or -n             : enforce not to run with 'fitting' even if the TBFIT tag is set to .TRUE."
   write(6,'(A)')"   -plot or -pl             : run gnuplot script after calculation with 'gnuplot gnuBAND-TB.gpi' "
   write(6,'(A)')"   -np                      : same as ' -nofit + -plot ' option applied "
   write(6,'(A)')"   -ls_type (or -ls) METHOD : set fitting METHOD manually. One of following can be selected: LMDIF, PSO, PSO+LMDIF"
   write(6,'(A)')"   -gnuplot xx.gpi          : with -np or -plot tag, the gnuplot command will be run with 'gnuplot xx.gpi' "
   write(6,'(A)')"   -print_only(or -po)      : enforce not to proceed fitting but print out weighting information and stop"
   write(6,'(A)')"   -miter MIT               : enforce to set maximum number of iteration for LMDIF to MIT "
   write(6,'(A)')"                              in prior to the MITER tag"
   write(6,'(A)')"   -mxfit(or -mf) MXF       : enforce to set maximum number of reapeat of iteration for LMDIF to MXF "
   write(6,'(A)')"                              in prior to the MITER tag"
   write(6,'(A)')"   -param(or -p) PF         : enforce to read parameter file 'PF' in prior to the PFILE tag"
   write(6,'(A)')"   -weight(or -w) WF        : enforce to read weight    file 'WF' in prior to the WFILE tag or 'SET WEIGHT'"
   write(6,'(A)')"   -kpoint(or -kp, -k) KF   : enforce to read k-points  file 'KF' in prior to the KFILE tag"
   write(6,'(A)')"   -lorbit                  : enforce to print orbital information"
   write(6,'(A)')"       .true.  or T         : enforce to print orbital information"
   write(6,'(A)')"       .false. or F         : enforce not to print orbital information"
   write(6,'(A)')"           mx               : enforce to print magnetization mx " 
   write(6,'(A)')"           my               : enforce to print magnetization mz " 
   write(6,'(A)')"           mz               : enforce to print magnetization mz " 
!  write(6,'(A)')"   -ldos .true.  or T       : enforce to print orbital information for each atom in separate file"
!  write(6,'(A)')"         .false. or F       : enforce not to print orbital information"
   write(6,'(A)')"   -pso_verbose  IVERBOSE   : determine verbosity in PSO routine. "
   write(6,'(A)')"                            :  -> IVERBOSE = 1 : write all info (incl. cost for each particles)  "
   write(6,'(A)')"                            :  -> IVERBOSE = 2 : write PSO results only"
   write(6,'(A)')"   -v IVERBOSE              : determin verbosity in the whole process."
   write(6,'(A)')"                            :  -> IVERBOSE = 1 : write all info "
   write(6,'(A)')"                            :  -> IVERBOSE = 2 : write no info neither on screen and file"
   write(6,'(A)')"                            :  -> IVERBOSE =-2 : write only to file "
   write(6,'(A)')"   -red_ovl  RED            : Multipy RED to overlap parameters o_*   "
   write(6,'(A)')"   -red_hop  RED            : Multipy RED to hopping parameters sps_, pps_, ... etc. "
   write(6,'(A)')"   -test                    : run test routine, for the development perpose only."

   return
endsubroutine
