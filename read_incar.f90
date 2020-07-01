#include "alias.inc"
module read_incar
   use mpi_setup
   use parameters
   use print_io
   implicit none

contains

   subroutine set_eigplot(PINPT,desc_str)
      type(incar)  ::  PINPT
      integer*4     i
      integer*4     i_continue
      integer*4     i_dummy, i_dummy2
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      character*40  plot_mode
      character*10  stm_window_tag
      character*40  dummy, dummy1, dummy2
      real*8        stm_emin(max_dummy), stm_emax(max_dummy)
      external      nitems 
      character(*), parameter :: func = 'set_eigplot'
      plot_mode = trim(desc_str)
      PINPT%flag_get_orbital = .true.
      PINPT%flag_get_band = .true.

mode: select case ( trim(plot_mode) )

        case('EIGPLOT')
            PINPT%flag_plot_eigen_state = .true.
  plot_eig: do while(trim(desc_str) .ne. 'END')
              read(pid_incar,'(A)',iostat=i_continue) inputline
              read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
              if(i_continue .ne. 0) cycle      ! skip empty line
              if(desc_str(1:1).eq.'#') cycle   ! skip comment
              if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

    eig_mode: select case ( trim(desc_str) )
                case('IEIG')
                  PINPT%n_eig_print = nitems(inputline) - 1
                  allocate(PINPT%i_eig_print(PINPT%n_eig_print))
                  read(inputline,*,iostat=i_continue) desc_str,PINPT%i_eig_print(1:PINPT%n_eig_print)
                  write(message,'(A,*(I6))')' EIG_PLOT:  ',PINPT%i_eig_print  ; write_msg
                  cycle plot_eig

                case('IKPT')
                  PINPT%n_kpt_print = nitems(inputline) - 1
                  allocate(PINPT%i_kpt_print(PINPT%n_kpt_print))
                  read(inputline,*,iostat=i_continue) desc_str,PINPT%i_kpt_print(1:PINPT%n_kpt_print)
                  write(message,'(A,*(I6))')' KPT_PLOT:  ',PINPT%i_kpt_print  ; write_msg
                  cycle plot_eig

                case('NGRID')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 3) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT%ngrid(1:3)
                    PINPT%flag_default_ngrid = .false.
                  else
                    write(message,'(A)')'    !WARN! NGRID tag of "SET EIGPLOT" should be three consequent integer numbers.'  ; write_msg
                    write(message,'(A,A)')'           Please check NGRID tag again. Exit... ',func  ; write_msg
                    stop
                  endif

                case('RORIGIN')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 3) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT%r_origin(1:3)
                    write(message,'(A,3(F15.8))')'   R_ORIG:  ',PINPT%r_origin(1:3) ; write_msg
                    PINPT%flag_default_rorigin = .false.
                  else
                    write(message,'(A)')'    !WARN! RORIGIN tag of "SET EIGPLOT" should be three consequent real values.'  ; write_msg
                    write(message,'(A,A)')'           Please check RORIGIN tag again. Exit... ',func  ; write_msg
                    stop
                  endif

                case('WAVEPLOT','WAV_PLOT')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_plot_wavefunction
                    write(message,'(A,L)')' WAV_PLOT:  ',PINPT%flag_plot_wavefunction  ; write_msg
                  else
                    write(message,'(A)')'    !WARN! WAV_PLOT tag of "SET EIGPLOT" should be .TRUE. or .FALSE.'  ; write_msg
                    write(message,'(A,A)')'           Please check WAV_PLOT tag again. Exit... ',func  ; write_msg
                    stop
                  endif

                case('RCUT', 'RCUT_ORB')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str, PINPT%rcut_orb_plot
                    write(message,'(A,F9.4)')' RCUT_ORB:  ',PINPT%rcut_orb_plot  ; write_msg
                  else
                    write(message,'(A)')'    !WARN! RCUT_ORB tag of "SET EIGPLOT" should be single real type parameter.'  ; write_msg
                    write(message,'(A,A)')'           Please check RCUT_ORB tag again. Exit... ',func  ; write_msg
                  endif
                case('REPEAT_CELL')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .ne. 3) then
                    write(message,'(A)')'   !WARN! REPEAT_CELL only accepts three logical arguments, for example "T T F",'  ; write_msg
                    write(message,'(A)')'          which implies that the orbitals along a1 and a2 direction are repeated and'  ; write_msg
                    write(message,'(A)')'          not for the a3 direction. check manual. Stop program...'  ; write_msg
                    stop
                  endif
                  read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_repeat_cell_orb_plot(1:3)

              end select eig_mode

            enddo plot_eig

        case('STMPLOT')
            PINPT%flag_plot_stm_image = .true.
            PINPT%flag_repeat_cell_orb_plot(1:3) = .true.

            write(message,'(A,L)')' STM_PLOT:  ',PINPT%flag_plot_stm_image  ; write_msg
            PINPT%n_stm= 0
  plot_stm: do while(trim(desc_str) .ne. 'END')
              read(pid_incar,'(A)',iostat=i_continue) inputline
              read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
              if(i_continue .ne. 0) cycle      ! skip empty line
              if(desc_str(1:1).eq.'#') cycle   ! skip comment
              if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

    stm_mode: select case ( trim(desc_str) )

                case('NGRID')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 3) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT%stm_ngrid(1:3)
                    PINPT%flag_default_stm_ngrid = .false.
                  else
                    write(message,'(A)')'    !WARN! NGRID tag of "SET STMPLOT" should be three consequent integer numbers.'  ; write_msg
                    write(message,'(A,A)')'           Please check NGRID tag again. Exit... ',func  ; write_msg
                    stop
                  endif

                case('STM_ERANGE','STM_WINDOW')
                  stm_window_tag=trim(desc_str)
                  PINPT%n_stm=PINPT%n_stm + 1
                  call strip_off (trim(inputline), dummy, stm_window_tag, ' ' , 2)   ! get dos_range
                  i_dummy=index(dummy,':')
                  if( i_dummy .eq. 0) then ! if ':' is not provided
                    i_dummy2 = nitems(dummy)
                    if(i_dummy2 .eq. 2)then ! if 'emax' and 'emin' is both provided
                      read(dummy,*,iostat=i_continue) stm_emin(PINPT%n_stm), stm_emax(PINPT%n_stm)
                      write(message,'(A,I2,A,F9.4,A,F9.4,A)')' STM_ERAN: PLOT(',PINPT%n_stm,')= [',stm_emin(PINPT%n_stm),':',stm_emax(PINPT%n_stm),']'  ; write_msg
                    else
                      write(message,'(A)')'    !WARNING!  STM_ERANGE is not properly set up.'  ; write_msg
                      write(message,'(A)')'    !WARNING!  Proper usage is as follows:'  ; write_msg
                      write(message,'(A)')'    !WARNING!    STM_ERANGE  EMIN:EMAX or EMIN EMAX'  ; write_msg
                      write(message,'(A)')'    !WARNING!  Exit program...'  ; write_msg
                      stop
                    endif
                  elseif(i_dummy .gt. 1) then ! if ':' is provided
                    call strip_off (trim(dummy), dummy1,' ',':',0)
                    call str2real(dummy1,stm_emin(PINPT%n_stm))
                    call strip_off (trim(dummy), dummy2,':',' ',2)
                    call str2real(dummy2,stm_emax(PINPT%n_stm))
                    write(message,'(A,I2,A,F9.4,A,F9.4,A)')' STM_ERAN: PLOT(',PINPT%n_stm,')= [',stm_emin(PINPT%n_stm),':',stm_emax(PINPT%n_stm),']'  ; write_msg
                  endif

                case('RCUT', 'RCUT_ORB')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str, PINPT%rcut_orb_plot
                    write(message,'(A,F9.4)')' RCUT_ORB:  ',PINPT%rcut_orb_plot  ; write_msg
                  else
                    write(message,'(A)')'    !WARN! RCUT_ORB tag of "SET EIGPLOT" should be single real type parameter.'  ; write_msg
                    write(message,'(A,A)')'           Please check RCUT_ORB tag again. Exit... ',func  ; write_msg
                  endif
                case('REPEAT_CELL')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .ne. 3) then
                    write(message,'(A)')'   !WARN! REPEAT_CELL only accepts three logical arguments, for example "T T F",'  ; write_msg
                    write(message,'(A)')'          which implies that the orbitals along a1 and a2 direction are repeated and'  ; write_msg
                    write(message,'(A)')'          not for the a3 direction. check manual. Stop program...'  ; write_msg
                    stop
                  endif
                  read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_repeat_cell_orb_plot(1:3)
                  
              end select stm_mode

            enddo plot_stm
            allocate(PINPT%stm_emax(PINPT%n_stm))
            allocate(PINPT%stm_emin(PINPT%n_stm))
            PINPT%stm_emax = stm_emax(1:PINPT%n_stm)
            PINPT%stm_emin = stm_emin(1:PINPT%n_stm)


      end select mode

      do i = 1, 3
        if(PINPT%flag_repeat_cell_orb_plot(i)) then
          PINPT%repeat_cell_orb_plot(i) = 1
        else
          PINPT%repeat_cell_orb_plot(i) = 0
          write(message,'(A,I1,A)')'   !WARN! You have set the orbitals will not be repeated along a',i,'-direction' ; write_msg
          write(message,'(A)')     '          in the STM or EIGPLOT. Proceed anyway...' ; write_msg
        endif
      enddo


      return    
   endsubroutine

   ! deprecated option .. use set_tbparam_file instead
   subroutine set_tbparam(PINPT,param_const,desc_str)
      type(incar)  ::  PINPT
      integer*4     i, ii, i_continue
      integer*4     i_dummy
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      character*40  dummy, dummy2
      real*8        r_dummy
      logical       flag_number
      external      nitems, flag_number
      real*8        param_const(5,max_nparam)
      character(*), parameter :: func = 'set_tbparam'

      if(PINPT%flag_pfile) then
        write(message,'(A,A,A)')'    !WARNING! ',trim(PINPT%pfilenm),' is alread provided,'  ; write_msg
        write(message,'(A,A,A)')'    !WARNING! ','but the TBPARAM is also provided in the INCAR-TB file.'  ; write_msg
        write(message,'(A,A,A)')'    !WARNING! ','Those values from INCAR-TB will not be read'  ; write_msg
      elseif(.not. PINPT%flag_pfile) then
        i=0
        do
          read(pid_incar,'(A)',iostat=i_continue) inputline
          read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
          if(i_continue .ne. 0) then
           i=i+1
           cycle              ! skip empty line
          elseif(desc_str(1:1).eq.'#') then
           i=i+1
           cycle  ! skip comment
          elseif(desc_str .eq. 'END' ) then
           i=i+1
           exit   ! finish reading TB-parameters
          else
           i=i+1
           PINPT%nparam=PINPT%nparam+1
          endif
        enddo
        do ii=1,i
          backspace(pid_incar)
        enddo
        if(PINPT%nparam .ne. 0) then
          allocate( PINPT%param(PINPT%nparam), PINPT%param_name(PINPT%nparam) )
          if(.not. PINPT%flag_pfile) then
            param_const(1,1:PINPT%nparam) = 0d0
            param_const(2,1:PINPT%nparam) = 20d0
            param_const(3,1:PINPT%nparam) =-20d0
            param_const(4,1:PINPT%nparam) = 0d0
            param_const(5,1:PINPT%nparam) = 0d0

          endif

          i=0
          do
            read(pid_incar,'(A)',iostat=i_continue) inputline
            read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
            if(i_continue .ne. 0) cycle              ! skip empty line
            if (desc_str(1:1).eq.'#') cycle  ! skip comment
            if (desc_str .eq. 'END' ) exit   ! finish reading TB-parameters
            i=i+1
            i_dummy = nitems(inputline) - 1

            if(i_dummy .eq. 1) then
              read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i)

            elseif( i_dummy .eq. 2) then
              read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i),dummy
              dummy = trim(dummy)
              if(.not. flag_number(dummy)) then
                if(.not. PINPT%flag_pfile) then
                  if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
                    param_const(4,i) = 1d0
                    param_const(5,i) = PINPT%param(i)
                  elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
                    param_const(4,i) = 0d0
                  endif
                endif
              elseif(flag_number(dummy)) then
                if(.not. PINPT%flag_pfile) then
                  call str2real(dummy, r_dummy)
                  PINPT%param(i) = PINPT%param(i) * r_dummy ! re_scaled
                endif
              endif

            elseif( i_dummy .eq. 3) then
              read(inputline,*,iostat=i_continue) PINPT%param_name(i),PINPT%param(i),dummy2, dummy
              dummy2= trim(dummy2)
              if(.not.flag_number(dummy2)) stop "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file."
              dummy = trim(dummy)
              if(flag_number(dummy))     stop "  !!WARN!! wrong syntax in PARAM.dat, check your PARAM file."

              call str2real(dummy2, r_dummy)
              PINPT%param(i) = PINPT%param(i) * r_dummy ! re_scaled
              
              if(.not. PINPT%flag_pfile) then
                if(dummy(1:1) .eq. 'F' .or. dummy(1:1) .eq. 'f') then ! if set 'fixed' or 'Fixed'
                  param_const(4,i) = 1d0
                  param_const(5,i) = PINPT%param(i)
                elseif( dummy(1:1) .eq. 'R' .or. dummy(1:1) .eq. 'r' ) then ! if set 'relaxed' or 'Relaxed'
                  param_const(4,i) = 0d0
                endif
              endif

            endif
          enddo
          PINPT%flag_pincar = .true.
        elseif ( PINPT%nparam .eq. 0 ) then
          write(message,'(A,A,A)')'    !WARNING! TBPARAM tag is set in the INCAR-TB but'  ; write_msg
          write(message,'(A,A,A)')'    !WARNING! TB-parameter is not provided.         '  ; write_msg
          write(message,'(A,A,A)')'    !WARNING! TB-parameter will be read from externl file if set by PFILE'  ; write_msg
        endif
      endif
      return
   endsubroutine

   subroutine set_constraint(PINPT,desc_str)
      type(incar)  ::  PINPT
      integer*4     i, ii, i_continue
      integer*4     nitems
      external      nitems
      character(*), parameter :: func = 'set_constraint'
      character*132 inputline
      character*40  desc_str, dummy, dummy1,dummy2,dummy3     

      i=0
      do
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) then
         i=i+1
         cycle              ! skip empty line
        elseif(desc_str(1:1).eq.'#') then
         i=i+1
         cycle  ! skip comment
        elseif(desc_str .eq. 'END' ) then
         i=i+1
         exit   ! finish reading parameter constraint
        else
         i=i+1
         PINPT%nparam_const=PINPT%nparam_const+1
        endif
      enddo
      do ii=1,i
        backspace(pid_incar)
      enddo

      if(PINPT%nparam_const .ne. 0) then
        PINPT%flag_set_param_const = .true.
        allocate( PINPT%c_const(3,PINPT%nparam_const) )
        i=0
        do
          read(pid_incar,'(A)',iostat=i_continue) inputline
          read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
          if(i_continue .ne. 0) cycle              ! skip empty line
          if (desc_str(1:1).eq.'#') cycle  ! skip comment
          if (desc_str .eq. 'END' ) exit   ! finish reading TB-parameters
          i=i+1
          read(inputline,*,iostat=i_continue) dummy1, dummy2, dummy3
          PINPT%c_const(1,i)=dummy1
          PINPT%c_const(2,i)=dummy2
          PINPT%c_const(3,i)=dummy3
        enddo
      endif

      return
   endsubroutine

   subroutine set_nn_class(PGEOM, desc_str)
      type(poscar)  ::  PGEOM
      integer*4     i, ii, i_continue
      integer*4     nitems
      external      nitems
      character(*), parameter :: func = 'set_nn_class'
      character*132 inputline
      character*40  desc_str, dummy, dummy1,dummy2,dummy3
      character*80  strip_nn_pair_(max_pair_type)
      real*8        strip_nn_dist_(max_pair_type), strip_nn_r0_(max_pair_type)

            i=0
            do while(trim(desc_str) .ne. 'END')
              read(pid_incar,'(A)',iostat=i_continue) inputline
              read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
              if(i_continue .ne. 0) cycle      ! skip empty line
              if(desc_str(1:1).eq.'#') cycle   ! skip comment
              if(trim(desc_str).eq.'END') exit ! exit loop if 'END'
              i=i+1
              call strip_off (trim(inputline), strip_nn_pair_(i), ' ', ':', 0)   ! get nn-pair type. header
              call strip_off (trim(inputline), dummy, ':', ' ', 2)   ! get nn-pair dist

              if( nitems(dummy) .eq. 1) then
                call str2real(trim(dummy),strip_nn_dist_(i))
                strip_nn_r0_(i)   = strip_nn_dist_(i)
                strip_nn_dist_(i) = strip_nn_dist_(i)*1.1d0
              elseif( nitems(dummy) .eq. 2) then
                read(dummy,*,iostat=i_continue) dummy1, dummy2
                call str2real(trim(dummy1),strip_nn_dist_(i))
                call str2real(trim(dummy2),strip_nn_r0_(i))
              elseif( nitems(dummy) .eq. 3) then
                read(dummy,*,iostat=i_continue) dummy1, dummy2, dummy3
                call str2real(trim(dummy1),strip_nn_dist_(i))
                call str2real(trim(dummy3),strip_nn_r0_(i))
              endif
              write(message,'(A,A12,2(A,1F8.4))')'  NN_PAIR:  ',strip_nn_pair_(i),' R0_max:',strip_nn_dist_(i),'   R0:',strip_nn_r0_(i)  ; write_msg
            enddo
            allocate( PGEOM%nn_pair(i) )
            allocate( PGEOM%nn_dist(i) )
            allocate( PGEOM%nn_r0(i)   )
            PGEOM%nn_pair(1:i) = strip_nn_pair_(1:i) ! Atom1-Atom2
            PGEOM%nn_dist(1:i) = strip_nn_dist_(1:i) ! maximum nn dist ATom1-Atom2
            PGEOM%nn_r0(1:i)   = strip_nn_r0_(1:i)   ! R0 dist ATom1-Atom2
            PGEOM%n_nn_type    = i


      return
   endsubroutine

   subroutine set_replot(PRPLT, desc_str)
      type(replot )  :: PRPLT
      integer*4     i, ii, k, i_continue
      integer*4     nitems
      external      nitems
      character(*), parameter :: func = 'set_replot'
      character*132 inputline
      logical       flag_number
      external      flag_number
!     character*40  desc_str, dummy, dummy1,dummy2,dummy3
      character*40  desc_str, dummy, dummy2,dummy3
      character*132 dummy1
      integer*4     i_dummy,i_dummy1,i_dummy2,i_dummy3,i_dummy4,i_dummy5
      character*40,  allocatable :: strip_dummy(:)
      integer*4     i_dummyr(max_dummy)
      integer*4     mpierr
      character*2   axis_dummy, axis_print_mag(max_dummy)
      character*4   form_dummy 
      logical       flag_replot_print_single(max_dummy)
      logical       flag_replot_write_unformatted(max_dummy)

      ! setup default values
!     PRPLT%replot_ldos_natom    =  0 
      PRPLT%replot_dos_smearing  =  0.025d0 
      PRPLT%replot_dos_emin      = -10d0
      PRPLT%replot_dos_emax      =  10d0
      PRPLT%replot_dos_nediv     =  1000
      PRPLT%replot_sldos_cell    = (/1,1,1/)
      PRPLT%r_origin             = 0d0
      PRPLT%bond_cut             = 3d0 
      PRPLT%flag_replot_formatted= .true.  ! default: read formatted band_structure_TBA file
!     PRPLT%replot_axis_print_mag= 'rh'



 set_rpl:do while(trim(desc_str) .ne. 'END')
           read(pid_incar,'(A)',iostat=i_continue) inputline
           read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
           if(i_continue .ne. 0) cycle      ! skip empty line
           if(desc_str(1:1).eq.'#') cycle   ! skip comment
           if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

  case_rpl:select case ( trim(desc_str) )
             case('REPLOT_DOS')
               read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_dos 
               if(PRPLT%flag_replot_dos) then
                 write(message,'(A)')'  REPLT_DOS: .TRUE.'  ; write_msg
               elseif(.not. PRPLT%flag_replot_dos) then
                 write(message,'(A)')'  REPLT_DOS: .FALSE.'  ; write_msg
               endif

             case('REPLOT_SLDOS')
               i_dummy = nitems(inputline) - 1
               if(i_dummy .eq. 1) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_sldos
                 if(PRPLT%flag_replot_sldos) then
                   write(message,'(A)')'REPLT_SLDOS: .TRUE. ==> written in ',trim(PRPLT%replot_sldos_fname)  ; write_msg
                 elseif(.not. PRPLT%flag_replot_sldos) then
                   write(message,'(A)')'REPLT_SLDOS: .FALSE.'  ; write_msg
                 endif
               elseif(i_dummy .eq. 2) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_sldos,PRPLT%replot_sldos_fname
                 if(PRPLT%flag_replot_sldos) then
                   write(message,'(2A)')'REPLT_SLDOS: .TRUE. ==> written in ',trim(PRPLT%replot_sldos_fname)  ; write_msg
                 elseif(.not. PRPLT%flag_replot_sldos) then
                   write(message,'(A)')'REPLT_SLDOS: .FALSE.'  ; write_msg
                 endif
               endif

             case('REPLOT_ONLY')
               read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_only
               if(PRPLT%flag_replot_only) then
                 write(message,'(A)')' REPLT_ONLY: .TRUE.'  ; write_msg
               elseif(.not. PRPLT%flag_replot_only) then
                 write(message,'(A)')' REPLT_ONLY: .FALSE.'  ; write_msg
               endif

             case('FILE_FORMAT')
               read(inputline,*,iostat=i_continue) desc_str, desc_str
               if(desc_str(1:3) .eq. 'bin') then
                 write(message,'(A)')' REPLT_FILE_FORM: read unformatted (binary) .bin file'  ; write_msg
                 PRPLT%flag_replot_formatted = .false.
               elseif(desc_str(1:3) .eq. 'asc' .or. desc_str(1:3) .eq. 'dat') then
                 write(message,'(A)')' REPLT_FILE_FORM: read formatted (ascii) .dat file'  ; write_msg
                 PRPLT%flag_replot_formatted = .true.
               endif

             case('REPLOT_BAND') ! convert band~.bin file with 'wf' format into band~replot.dat with 'axis_print_mag' format.
               i_dummy = nitems(inputline) - 1
               if(i_dummy .eq. 1) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_band
                 if(PRPLT%flag_replot_band) then
                   write(message,'(A)')' REPLT_BAND: .TRUE. with <phi_i|psi_nk> (rh) ; phi_i : atomic orbital'  ; write_msg

                   PRPLT%replot_nband = PRPLT%replot_nband + 1
                   axis_print_mag(PRPLT%replot_nband) = 'rh'
                   flag_replot_print_single(PRPLT%replot_nband)= .FALSE.
                   flag_replot_write_unformatted(PRPLT%replot_nband) = .FALSE.
                   
                 elseif(.not. PRPLT%flag_replot_band) then
                   write(message,'(A)')' REPLT_BAND: .FALSE.  rh'  ; write_msg
                 endif
               elseif(i_dummy .eq. 2) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_band, axis_dummy
                 if(PRPLT%flag_replot_band) then
                   if    (axis_dummy .eq. 'mx') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <psi_nk|sigma_x|psi_nk> (mx)'  ; write_msg
                   elseif(axis_dummy .eq. 'my') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <psi_nk|sigma_y|psi_nk> (my)'  ; write_msg
                   elseif(axis_dummy .eq. 'mz') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <psi_nk|sigma_z|psi_nk> (mz)'  ; write_msg
                   elseif(axis_dummy .eq. 'wf') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with total wavefunction (wf)'  ; write_msg
                   elseif(axis_dummy .eq. 'rh') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <phi_i|psi_nk> (rh) ; phi_i : atomic orbital'  ; write_msg
                   elseif(axis_dummy .eq. 'no') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. only with eigenvalues'  ; write_msg
                   endif
                   
                   PRPLT%replot_nband = PRPLT%replot_nband + 1
                   axis_print_mag(PRPLT%replot_nband) = axis_dummy
                   flag_replot_print_single(PRPLT%replot_nband)= .FALSE.
                   flag_replot_write_unformatted(PRPLT%replot_nband) = .FALSE.

                 elseif(.not. PRPLT%flag_replot_band) then
                   write(message,'(2A)')' REPLT_BAND: .FALSE.  ',trim(axis_dummy)  ; write_msg
                 endif
               elseif(i_dummy .eq. 3) then

                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_band, axis_dummy, form_dummy
                 if(PRPLT%flag_replot_band) then
                   if    (axis_dummy .eq. 'mx') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <psi_nk|sigma_x|psi_nk> (mx)'  ; write_msg
                   elseif(axis_dummy .eq. 'my') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <psi_nk|sigma_y|psi_nk> (my)'  ; write_msg
                   elseif(axis_dummy .eq. 'mz') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <psi_nk|sigma_z|psi_nk> (mz)'  ; write_msg
                   elseif(axis_dummy .eq. 'wf') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with total wavefunction (wf)'  ; write_msg
                   elseif(axis_dummy .eq. 'rh') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. with <phi_i|psi_nk> (rh) ; phi_i : atomic orbital'  ; write_msg
                   elseif(axis_dummy .eq. 'no') then
                     write(message,'(2A)')' REPLT_BAND: .TRUE. only with eigenvalues'  ; write_msg
                   endif

                   PRPLT%replot_nband = PRPLT%replot_nband + 1
                   axis_print_mag(PRPLT%replot_nband) = axis_dummy

                   if(form_dummy(1:3) .eq. 'bin') then
                     flag_replot_write_unformatted(PRPLT%replot_nband) = .TRUE. 
                     if(form_dummy(1:4) .eq. 'bin4') then
                       flag_replot_print_single(PRPLT%replot_nband)= .TRUE.
                       write(message,'(2A)')'           : with binary format (single).'  ; write_msg
                     else
                       flag_replot_print_single(PRPLT%replot_nband)= .FALSE.
                       write(message,'(2A)')'           : with binary format.'  ; write_msg
                     endif
                   else
                     flag_replot_write_unformatted(PRPLT%replot_nband) = .FALSE. 
                     flag_replot_print_single(PRPLT%replot_nband)= .FALSE.
                   endif

                 elseif(.not. PRPLT%flag_replot_band) then
                   write(message,'(2A)')' REPLT_BAND: .FALSE.  ',trim(axis_dummy)  ; write_msg
                 endif

               endif

!            case('REPLOT_EIG')
!              i_dummy = nitems(inputline) - 1
!              if(i_dummy .eq. 1) then
!                read(inputline, *, iostate=i_continue) desc_str, PRPLT%flag_replot_eig, 
!                if(PRPLT%flag_replot_eig) then
!                  write(message,'(A)')' REPLT_EIG: .TRUE.'  ; write_msg
!                endif
!              endif

             case('REPLOT_PROJ_BAND', 'REPLOT_PROJ_SUM')
               i_dummy = nitems(inputline) - 1
               if(i_dummy .eq. 1) then
                 read(inputline,*,iostat=i_continue) desc_str, PRPLT%flag_replot_proj_band
                 if(PRPLT%flag_replot_proj_band) then
                   write(message,'(A)')'  !!!!WARN: REPLOT_PROJ_BAND tag is activated but atom indext to be summed up '  ; write_msg
                   write(message,'(A)')'            did not specified.'  ; write_msg
                   write(message,'(A)')'            Correct usage is, for example, if you want to resolve'  ; write_msg
                   write(message,'(A)')'            contribution of atom 1 to 10 and 15, then you can specify as follows'  ; write_msg
                   write(message,'(A)')'            REPLOT_PROJ_BAND .TRUE.  1:10 15'  ; write_msg
                   write(message,'(A)')'  Stop program...'  ; write_msg
                   kill_job
                 endif

                 if( .not. PRPLT%flag_replot_proj_band) then
                   write(message,'(A)')' REPLT_PROJ_BAND: .FALSE.'  ; write_msg
                 endif

               elseif(i_dummy .gt. 1) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_proj_band
                 read(inputline,*,iostat=i_continue) desc_str,dummy
                 if(PRPLT%flag_replot_proj_band) then
                   PRPLT%replot_nproj_sum = PRPLT%replot_nproj_sum + 1
                   call strip_off (trim(inputline), dummy1, trim(dummy), ' ' , 2)   ! get dos_ensurf
                   i_dummy1=index(dummy1,':')
                   if(i_dummy1 .eq. 0) then
                     if(PRPLT%replot_nproj_sum .eq. 1) then 
                       allocate( PRPLT%replot_proj_natom(max_dummy2) )
                       allocate( PRPLT%replot_proj_atom(max_dummy2*10,max_dummy2) )
                       PRPLT%replot_ldos_natom = 0
                       PRPLT%replot_ldos_atom  = 0
                     endif
                     PRPLT%replot_proj_natom(PRPLT%replot_nproj_sum) = nitems(dummy1)
                     read(dummy1,*,iostat=i_continue) PRPLT%replot_proj_atom(1:PRPLT%replot_proj_natom(PRPLT%replot_nproj_sum), PRPLT%replot_nproj_sum)
                     write(message,'(A,I0,2A)')'REPLT_PROJ_BAND: .TRUE. , Atom_index SET ',PRPLT%replot_nproj_sum,' = ',trim(dummy1)  ; write_msg
                   elseif(i_dummy1 .ge. 1)then
                     i_dummy2 = nitems(dummy1)
                     allocate( strip_dummy(i_dummy2) )
                     read(dummy1,*,iostat=i_continue) (strip_dummy(i),i=1,i_dummy2)
                     ii = 0
                     do i = 1, i_dummy2
                       i_dummy3 = index(strip_dummy(i),':')
                       if(i_dummy3 .eq. 0) then
                         ii = ii + 1
                         call str2int(strip_dummy(i),i_dummy4)
                         i_dummyr(ii) = i_dummy4
                       elseif(i_dummy3 .gt. 1) then
                         ii = ii + 1
                         call strip_off(trim(strip_dummy(i)), dummy3, ' ', ':', 0)
                         call str2int(dummy3,i_dummy4)
                         call strip_off(trim(strip_dummy(i)), dummy3, ':', ' ', 2)
                         call str2int(dummy3,i_dummy5)
                         i_dummyr(ii:ii+i_dummy5 - i_dummy4) = (/ (k, k=i_dummy4, i_dummy5) /)
                         ii = ii + i_dummy5 - i_dummy4
                       endif
                     enddo
                     if(PRPLT%replot_nproj_sum .eq. 1) then
                       allocate( PRPLT%replot_proj_natom(max_dummy2) )
                       allocate( PRPLT%replot_proj_atom(max_dummy2*10,max_dummy2) )
                       PRPLT%replot_ldos_natom = 0
                       PRPLT%replot_ldos_atom  = 0
                     endif
                     PRPLT%replot_proj_natom(PRPLT%replot_nproj_sum) = ii
                     deallocate( strip_dummy )
                     PRPLT%replot_proj_atom(1:ii,PRPLT%replot_nproj_sum) = i_dummyr(1:ii)
                     write(message,'(A,I0,2A)')'REPLT_PROJ_BAND: .TRUE. , Atom_index SET ',PRPLT%replot_nproj_sum,' = ',trim(dummy1)  ; write_msg
                   endif
                 endif

               endif


             case('REPLOT_LDOS')
               i_dummy = nitems(inputline) - 1
               if(i_dummy .eq. 1) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_ldos
                 if(PRPLT%flag_replot_ldos) then
                   write(message,'(A)')'  !WARNING! DOS_LDOS -> .TRUE. but the list of target atoms is not specified.'  ; write_msg
                   write(message,'(A)')'            Correct usage is, for example, if you want to resolve'  ; write_msg
                   write(message,'(A)')'            contribution of atom 1 to 10 and 15, then you can specify as follows'  ; write_msg
                   write(message,'(A)')'            REPLOT_LDOS .TRUE.  1:10 15'  ; write_msg
                   write(message,'(A)')'            If you want to resolve all atoms, then you can specify'  ; write_msg
                   write(message,'(A)')'            REPLOT_LDOS .TRUE.  1:NATOM '  ; write_msg
                   write(message,'(A)')'            where "NATOM" must be replaced implicitly according to the provided geometry file'  ; write_msg
                   write(message,'(A)')'  Stop program...'  ; write_msg
                   kill_job
                 endif

               elseif(i_dummy .gt. 1) then

                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_ldos
                 read(inputline,*,iostat=i_continue) desc_str,dummy
                 if(PRPLT%flag_replot_ldos) then
                   PRPLT%replot_nldos_sum = PRPLT%replot_nldos_sum + 1
                   call strip_off (trim(inputline), dummy1, trim(dummy), ' ' , 2)   ! get dos_ensurf
                   i_dummy1=index(dummy1,':')
                   if(i_dummy1 .eq. 0) then
                     if(PRPLT%replot_nldos_sum .eq. 1) then 
                       allocate( PRPLT%replot_ldos_natom(max_dummy2) )
                       allocate( PRPLT%replot_ldos_atom(max_dummy2*10,max_dummy2) )
                       PRPLT%replot_ldos_natom = 0
                       PRPLT%replot_ldos_atom  = 0
                     endif
                     PRPLT%replot_ldos_natom(PRPLT%replot_nldos_sum) = nitems(dummy1)
                     read(dummy1,*,iostat=i_continue) PRPLT%replot_ldos_atom(1:PRPLT%replot_ldos_natom(PRPLT%replot_nldos_sum), PRPLT%replot_nldos_sum)
                     write(message,'(A,I0,2A)')' REPLT_LDOS: .TRUE. , Atom_index SET ',PRPLT%replot_nldos_sum,' = ',trim(dummy1)  ; write_msg
                   elseif(i_dummy1 .ge. 1)then
                     i_dummy2 = nitems(dummy1)
                     allocate( strip_dummy(i_dummy2) )
                     read(dummy1,*,iostat=i_continue) (strip_dummy(i),i=1,i_dummy2)
                     ii = 0
                     do i = 1, i_dummy2
                       i_dummy3 = index(strip_dummy(i),':')
                       if(i_dummy3 .eq. 0) then
                         ii = ii + 1
                         call str2int(strip_dummy(i),i_dummy4)
                         i_dummyr(ii) = i_dummy4
                       elseif(i_dummy3 .gt. 1) then
                         ii = ii + 1
                         call strip_off(trim(strip_dummy(i)), dummy3, ' ', ':', 0)
                         call str2int(dummy3,i_dummy4)
                         call strip_off(trim(strip_dummy(i)), dummy3, ':', ' ', 2)
                         call str2int(dummy3,i_dummy5)
                         i_dummyr(ii:ii+i_dummy5 - i_dummy4) = (/ (k, k=i_dummy4, i_dummy5) /)
                         ii = ii + i_dummy5 - i_dummy4
                       endif
                     enddo
                     if(PRPLT%replot_nldos_sum .eq. 1) then 
                       allocate( PRPLT%replot_ldos_natom(max_dummy2) )
                       allocate( PRPLT%replot_ldos_atom(max_dummy2*10,max_dummy2) )
                       PRPLT%replot_ldos_natom = 0
                       PRPLT%replot_ldos_atom  = 0
                     endif
                     PRPLT%replot_ldos_natom(PRPLT%replot_nldos_sum) = ii
                     deallocate( strip_dummy )
                     PRPLT%replot_ldos_atom(1:ii,PRPLT%replot_nldos_sum) = i_dummyr(1:ii)
                     write(message,'(A,I0,2A)')' REPLT_LDOS: .TRUE. , Atom_index SET ',PRPLT%replot_nldos_sum,' = ', trim(dummy1)  ; write_msg
                   endif
                 endif
               endif

             case('REPLOT_DIDV', 'REPLOT_STS') 
               i_dummy = nitems(inputline) - 1
               if(i_dummy .eq. 1) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_didv
                 if(PRPLT%flag_replot_didv) then
                   write(message,'(2A)')' REPLT_DIDV: .TRUE. ==> written in ',trim(PRPLT%replot_didv_fname)  ; write_msg
                 elseif(.not. PRPLT%flag_replot_didv) then
                   write(message,'(A)')' REPLT_DIDV: .FALSE.'  ; write_msg
                 endif
               elseif(i_dummy .eq. 2) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%flag_replot_didv,PRPLT%replot_didv_fname
                 if(PRPLT%flag_replot_didv) then
                   write(message,'(2A)')' REPLT_DIDV: .TRUE. ==> written in ',trim(PRPLT%replot_didv_fname)  ; write_msg
                 elseif(.not. PRPLT%flag_replot_didv) then
                   write(message,'(A)')' REPLT_DIDV: .FALSE.'  ; write_msg
                 endif
               endif

             case('NEDOS','REPLOT_NEDOS')
               read(inputline,*,iostat=i_continue) desc_str,PRPLT%replot_dos_nediv
               write(message,'(A,I8)')' REPLT_NDIV:', PRPLT%replot_dos_nediv  ; write_msg
               if(PRPLT%replot_dos_nediv .le. 0) then
                 write(message,'(3A,I0)')'    !WARNING! ',trim(desc_str),' should be larger than or equal to 1. current value: ',PRPLT%replot_dos_nediv ; write_msg
                 write(message,'(A)')'               Exit program...'  ; write_msg
                 kill_job
               endif

             case('DOS_ERANGE', 'DOS_EWINDOW', 'REPLOT_ERANGE', 'REPLOT_EWINDOW')
               call strip_off (trim(inputline), dummy, trim(desc_str), ' ' , 2)   ! get dos_range
               i_dummy=index(dummy,':')
               call strip_off (trim(dummy), dummy1,' ',':',0)
               if( i_dummy .eq. 0) then
                 i_dummy2 = nitems(dummy)
                 if(i_dummy2 .eq. 2)then
                   read(dummy,*,iostat=i_continue) PRPLT%replot_dos_emin,PRPLT%replot_dos_emax
                   write(message,'(A,F15.8)')' REPLT_EMIN:  ',PRPLT%replot_dos_emin  ; write_msg
                   write(message,'(A,F15.8)')' REPLT_EMAX:  ',PRPLT%replot_dos_emax  ; write_msg
                 else
                   write(message,'(A)')'    !WARNING!  REPLOT_EWINDOW is not properly set up.'  ; write_msg
                   write(message,'(A)')'    !WARNING!  Proper usage is as follows:'  ; write_msg
                   write(message,'(A)')'    !WARNING!    REPLOT_WINDOW  EMIN:EMAX , :EMAX, EMIN: or :'  ; write_msg
                   write(message,'(A)')'    !WARNING! or REPLOT_EWINDOW  EMIN EMAX'  ; write_msg
                   write(message,'(A)')'    !WARNING!  Exit program...'  ; write_msg
                   kill_job
                 endif
               elseif(i_dummy .ge. 1) then
                 if(len_trim(dummy1) .eq. 0) then
                   PRPLT%replot_dos_emin = -10.0d0 ! default dos_emin
                   write(message,'(A,F15.8)')' REPLT_EMIN:  ',PRPLT%replot_dos_emin  ; write_msg
                 else
                   call str2real(dummy1,PRPLT%replot_dos_emin)
                   write(message,'(A,F15.8)')' REPLT_EMIN:  ',PRPLT%replot_dos_emin  ; write_msg
                 endif
                 call strip_off (trim(dummy), dummy2,':',' ',2)
                 if(len_trim(dummy2) .eq. 0) then
                   PRPLT%replot_dos_emax =  10.0d0 ! default dos_emax
                   write(message,'(A,F15.8)')' REPLT_EMAX:  ',PRPLT%replot_dos_emax  ; write_msg
                 else
                   call str2real(dummy2,PRPLT%replot_dos_emax)
                   write(message,'(A,F15.8)')' REPLT_EMAX:  ',PRPLT%replot_dos_emax  ; write_msg
                 endif
               endif

             case('SMEARING', 'REPLOT_SMEARING') ! gaussian smearing
               read(inputline,*,iostat=i_continue) desc_str,PRPLT%replot_dos_smearing
               write(message,'(A,F8.4)')'REPLT_SIGMA: GAUSSIAN WIDTH = ',PRPLT%replot_dos_smearing  ; write_msg

             case('REPEAT_CELL')
               i_dummy = nitems(inputline) -1
               if(i_dummy .eq. 1) then
                 read(inputline,*,iostat=i_continue) desc_str, i_dummy2 
                 PRPLT%replot_sldos_cell = i_dummy2
               elseif(i_dummy .eq. 2) then
                 read(inputline,*,iostat=i_continue) desc_str, PRPLT%replot_sldos_cell(1:2)
               elseif(i_dummy .eq. 3) then
                 read(inputline,*,iostat=i_continue) desc_str, PRPLT%replot_sldos_cell(1:3) 
               endif
               write(message,'(A,3(I6))')' REPLT_CELL:',PRPLT%replot_sldos_cell(1:3)  ; write_msg

             case('RORIGIN')
               i_dummy = nitems(inputline) - 1
               if(i_dummy .eq. 3) then
                 read(inputline,*,iostat=i_continue) desc_str,PRPLT%r_origin(1:3)
                 write(message,'(A,3(F15.8))')' REPLT_ORIG:  ',PRPLT%r_origin(1:3) ; write_msg
               else
                 write(message,'(A)')'    !WARN! RORIGIN tag of "SET REPLOT" should be three consequent real values.'  ; write_msg
                 write(message,'(A,A)')'           Please check RORIGIN tag again. Exit... ',func  ; write_msg
                 kill_job
               endif

             case('BOND_CUT')
               read(inputline,*,iostat=i_continue) desc_str,PRPLT%bond_cut
               write(message,'(A,3(F15.8))')' REPLT_RCUT:  ',PRPLT%bond_cut ; write_msg

           endselect case_rpl

         enddo set_rpl

!     if(PRPLT%flag_replot_formatted .and. PRPLT%flag_replot_band) then
!       write(message,'(A)') '    !WARN! FILE_FORMAT for band_structure_TBA has been set to '  ; write_msg
!       write(message,'(A)') '           ascii(formatted) and also requested REPLOT_BAND to .TRUE.'  ; write_msg
!       write(message,'(A)') '           which is not accepted offer.'  ; write_msg
!       write(message,'(A)') '           Note that REPLOT_BAND is to convert band_structure_TBA.bin to'  ; write_msg
!       write(message,'(A)') '           band_structure_TBA.dat which is ascii (formatted) format.'  ; write_msg
!       write(message,'(A)') '           Exit program...'  ; write_msg
!       kill_job
!     endif

      if(PRPLT%replot_nproj_sum .ge. 1) PRPLT%flag_replot_proj_band = .true.
      if(PRPLT%replot_nldos_sum .ge. 1) PRPLT%flag_replot_ldos = .true.

      if(PRPLT%replot_nband     .ge. 1) then 
        PRPLT%flag_replot_band = .true.
        allocate(PRPLT%replot_axis_print_mag(PRPLT%replot_nband))
        allocate(PRPLT%flag_replot_print_single(PRPLT%replot_nband))
        allocate(PRPLT%flag_replot_write_unformatted(PRPLT%replot_nband))
        do i=1, PRPLT%replot_nband
          PRPLT%replot_axis_print_mag(i) = axis_print_mag(i)
          PRPLT%flag_replot_print_single(i) = flag_replot_print_single(i)
          PRPLT%flag_replot_write_unformatted(i) = flag_replot_write_unformatted(i)
        enddo
      endif
    
      if(PRPLT%replot_dos_emin .eq. PRPLT%replot_dos_emax .and. PRPLT%replot_dos_nediv .ne. 1) then
         PRPLT%replot_dos_nediv = 1
         write(message,'(A)') '    !WARN! EMIN and EMAX value of "DOS_EWINDOW" tag of "SET REPLOT" has been set to be equal: EMIN=EMAX'  ; write_msg
         write(message,'(A)') '           and REPLOT_NEDOS or (NEDOS) is larger than 2 which is nonsense, hence, REPLOT_NEDOS is '  ; write_msg
         write(message,'(A)') '           enfornced to be 1. Proceed calculations..'  ; write_msg
      endif
      return
   endsubroutine
   subroutine set_density_of_states(PINPT, PINPT_DOS, desc_str)
      type(incar )  ::  PINPT
      type(dos)     :: PINPT_DOS
      integer*4     i, ii, k, i_continue
      integer*4     nitems
      external      nitems
      character(*), parameter :: func = 'set_dos'
      character*132 inputline
      character*40  desc_str, dummy, dummy1,dummy2,dummy3
      integer*4     i_dummy,i_dummy1,i_dummy2,i_dummy3,i_dummy4,i_dummy5
      character*40,  allocatable :: strip_dummy(:)
      integer*4     i_dummyr(max_dummy)
      logical       flag_number
      external      flag_number
    
            PINPT%flag_get_dos = .true.
!           PINPT%flag_get_band = .true.
            write(message,'(A)')'  GET_DOS: .TRUE.'  ; write_msg
            !initialize the default values
            PINPT_DOS%dos_kgrid(1:3) = 0d0
            PINPT_DOS%dos_emin = -10d0
            PINPT_DOS%dos_emax =  10d0
            PINPT_DOS%dos_nediv = 1000
            PINPT_DOS%dos_kshift   = 0d0
            PINPT_DOS%dos_flag_gamma = .false.
            PINPT_DOS%dos_flag_print_kpoint = .false.
            PINPT_DOS%dos_kfilenm    = 'IBZKPT-DOS_TB'
            PINPT_DOS%dos_filenm     = 'DOS_TB_projected.dat'
            PINPT_DOS%ldos_filenm    = 'LDOS_TB_projected'   ! default. atom index will be appended after.
            PINPT_DOS%dos_flag_print_eigen = .false.
            PINPT_DOS%dos_kunit = 'R'
            PINPT_DOS%dos_iband = 1
            PINPT_DOS%dos_fband = 999999
            PINPT_DOS%dos_smearing = 0.025
            PINPT_DOS%dos_flag_sparse = .false.
            PINPT_DOS%dos_flag_print_ldos = .false.
            PINPT_DOS%dos_ldos_natom = 0

   set_dos: do while(trim(desc_str) .ne. 'END')
              read(pid_incar,'(A)',iostat=i_continue) inputline
              read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
              if(i_continue .ne. 0) cycle      ! skip empty line
              if(desc_str(1:1).eq.'#') cycle   ! skip comment
              if(trim(desc_str).eq.'END') exit ! exit loop if 'END'
    case_dos: select case ( trim(desc_str) )
                case('NEDOS')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_nediv
                  write(message,'(A,I8)')' DOS_EDIV:', PINPT_DOS%dos_nediv  ; write_msg
                case('DOS_ERANGE', 'DOS_EWINDOW')
                  call strip_off (trim(inputline), dummy, trim(desc_str), ' ' , 2)   ! get dos_range
                  i_dummy=index(dummy,':')
                  call strip_off (trim(dummy), dummy1,' ',':',0)
                  if( i_dummy .eq. 0) then
                    i_dummy2 = nitems(dummy)
                    if(i_dummy2 .eq. 2)then
                      read(dummy,*,iostat=i_continue) PINPT_DOS%dos_emin,PINPT_DOS%dos_emax
                      write(message,'(A,F15.8)')' DOS_EMIN:  ',PINPT_DOS%dos_emin  ; write_msg
                      write(message,'(A,F15.8)')' DOS_EMAX:  ',PINPT_DOS%dos_emax  ; write_msg
                    else
                      write(message,'(A)')'    !WARNING!  DOS_ERANGE is not properly set up.'  ; write_msg
                      write(message,'(A)')'    !WARNING!  Proper usage is as follows:'  ; write_msg
                      write(message,'(A)')'    !WARNING!    DOS_ERANGE  EMIN:EMAX , :EMAX, EMIN: or :'  ; write_msg
                      write(message,'(A)')'    !WARNING! or DOS_ERANGE  EMIN EMAX'  ; write_msg
                      write(message,'(A)')'    !WARNING!  Exit program...'  ; write_msg
                      stop
                    endif
                  elseif(i_dummy .ge. 1) then
                    if(len_trim(dummy1) .eq. 0) then
                      PINPT_DOS%dos_emin = -10.0d0 ! default dos_emin
                      write(message,'(A,F15.8)')' DOS_EMIN:  ',PINPT_DOS%dos_emin  ; write_msg
                    else
                      call str2real(dummy1,PINPT_DOS%dos_emin)
                      write(message,'(A,F15.8)')' DOS_EMIN:  ',PINPT_DOS%dos_emin  ; write_msg
                    endif
                    call strip_off (trim(dummy), dummy2,':',' ',2)
                    if(len_trim(dummy2) .eq. 0) then
                      PINPT_DOS%dos_emax =  10.0d0 ! default dos_emax
                      write(message,'(A,F15.8)')' DOS_EMAX:  ',PINPT_DOS%dos_emax  ; write_msg
                    else
                      call str2real(dummy2,PINPT_DOS%dos_emax)
                      write(message,'(A,F15.8)')' DOS_EMAX:  ',PINPT_DOS%dos_emax  ; write_msg
                    endif
                  endif
                case('DOS_NRANGE', 'DOS_NE_MAX')
                  call strip_off (trim(inputline), dummy, trim(desc_str), ' ' , 2)   ! get dos_range
                  i_dummy=index(dummy,':')
                  call strip_off (trim(dummy), dummy1,' ',':',0)
                  if( i_dummy .eq. 0) then
                    i_dummy2 = nitems(dummy)
                    if(i_dummy2 .eq. 2)then
                      read(dummy,*,iostat=i_continue) PINPT_DOS%dos_iband,PINPT_DOS%dos_fband
                      write(message,'(A,I8)')' DOS_IBND:',PINPT_DOS%dos_iband  ; write_msg
                      write(message,'(A,I8)')' DOS_FBND:',PINPT_DOS%dos_fband  ; write_msg
                    elseif(i_dummy2 .eq. 1) then
                      read(dummy,*,iostat=i_continue) PINPT_DOS%dos_fband
                      PINPT_DOS%dos_iband = 1
                      write(message,'(A,I8)')' DOS_IBND: (default)',PINPT_DOS%dos_iband  ; write_msg
                      write(message,'(A,I8)')' DOS_FBND:',PINPT_DOS%dos_fband  ; write_msg
                    else
                      write(message,'(A)')'    !WARNING!  DOS_NRANGE is not properly set up.'  ; write_msg
                      write(message,'(A)')'    !WARNING!  Proper usage is as follows:'  ; write_msg
                      write(message,'(A)')'    !WARNING!    DOS_NRANGE  IBAND:FBAND , :FBAND, IBAND: or :'  ; write_msg
                      write(message,'(A)')'    !WARNING! or DOS_NRANGE  IBAND FBAND'  ; write_msg
                      write(message,'(A)')'    !WARNING!  Exit program...'  ; write_msg
                      stop
                    endif
                  elseif(i_dummy .ge. 1 .and. len_trim(dummy) .ne. 1) then
                    if(len_trim(dummy1) .eq. 0) then
                      PINPT_DOS%dos_iband = 1 ! default dos_emin
                      write(message,'(A,I8)')' DOS_IBND:',PINPT_DOS%dos_iband  ; write_msg
                    else
                      call str2int(dummy1,PINPT_DOS%dos_iband)
                      write(message,'(A,I8)')' DOS_IBND:',PINPT_DOS%dos_iband  ; write_msg
                    endif
                    call strip_off (trim(dummy), dummy2,':',' ',2)
                    if(len_trim(dummy2) .eq. 0) then
                      PINPT_DOS%dos_fband =  999999 ! default dos_emax
                      write(message,'(A)')' DOS_FBND: automatically set to NEIG'  ; write_msg
                    else
                      if(flag_number( trim(dummy2) ) ) then
                        call str2int(dummy2,PINPT_DOS%dos_fband)
                        write(message,'(A,I8)')' DOS_FBND:',PINPT_DOS%dos_fband  ; write_msg
                      elseif(.not. flag_number( trim(dummy2) ) .and. trim(dummy2) .eq. 'NEIG' ) then
                        PINPT_DOS%dos_fband =  999999
                        write(message,'(A)')' DOS_FBND: NEIG = N_ORBIT'  ; write_msg
                      else
                        write(message,'(A,I8)')'  !WARNING! DOS_FBND setting is inproper. Please check syntax again. Exit..'  ; write_msg
                        stop
                      endif
                    endif
                  elseif(i_dummy .eq. 1 .and. len_trim(dummy) .eq. 1) then
                    PINPT_DOS%dos_iband = 1
                    write(message,'(A,I8)')' DOS_IBND: 1'  ; write_msg
                    write(message,'(A,I8)')' DOS_FBND: NEIG'  ; write_msg
                  endif
                case('MKGRID')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kgrid(1:3)
                  write(message,'(A,3I6)')' DOS_KDIV: Monkhorst-Pack grid',PINPT_DOS%dos_kgrid(1:3)  ; write_msg
                  PINPT_DOS%dos_flag_gamma=.false.
                case('GKGRID')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kgrid(1:3)
                  write(message,'(A,3I6)')' DOS_KDIV: Gamma-centered grid',PINPT_DOS%dos_kgrid(1:3)  ; write_msg
                  PINPT_DOS%dos_flag_gamma=.true.
                case('KGRID')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kgrid(1:3)
                  write(message,'(A,3I6)')' DOS_KDIV: Monkhorst-Pack grid',PINPT_DOS%dos_kgrid(1:3)  ; write_msg
                  PINPT_DOS%dos_flag_gamma=.false.
                case('KSHIFT') ! optional shift of the mesh (s_1, s_2, s_3) ; the usage is same as VASP 
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kshift(1:3)
                  write(message,'(A,3F6.2)')' DOS_KSFT:',PINPT_DOS%dos_kshift(1:3)  ; write_msg
                case('PRINT_KPTS') ! print kpoint information (reciprocal unit) into file
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_kpoint
                  elseif(i_dummy .eq. 2) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_kpoint, PINPT_DOS%dos_kfilenm
                  endif
                  if(PINPT_DOS%dos_flag_print_kpoint) then
                    write(message,'(A,A)')' DOS_KOUT: KPOINT OUT -> .TRUE. ',trim(PINPT_DOS%dos_kfilenm)  ; write_msg
                  elseif(.not. PINPT_DOS%dos_flag_print_kpoint) then
                    write(message,'(A,A)')' DOS_KOUT: KPOINT OUT -> .FALSE.'  ; write_msg
                  endif
                case('PRINT_EIG') ! print eigenstate energy information into ENSURF.EIG.NEIG.dat file
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_eigen
                    if(PINPT_DOS%dos_flag_print_eigen) then
                      write(message,'(A)')'  !WARNING! DOS_EOUT: EIGEN  OUT -> .TRUE. but the target energy to be plotted is not specified.'  ; write_msg
                      write(message,'(A)')'  !WARNING! Please check PRINT_EIG tag again. The proper use is : PRINT_EIG .TRUE. N1 N2 N3... or N1:N3 '  ; write_msg
                      write(message,'(A)')'  !WARNING! Exit program...'  ; write_msg
                      stop
                    endif
                  elseif(i_dummy .gt. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_eigen
                    read(inputline,*,iostat=i_continue) desc_str,dummy
                    if(PINPT_DOS%dos_flag_print_eigen) then
                      call strip_off (trim(inputline), dummy1, trim(dummy), ' ' , 2)   ! get dos_ensurf
                      i_dummy1=index(dummy1,':')
                      if(i_dummy1 .eq. 0) then
                        PINPT_DOS%dos_n_ensurf = nitems(dummy1)
                        allocate( PINPT_DOS%dos_ensurf(PINPT_DOS%dos_n_ensurf) )
                        read(dummy1,*,iostat=i_continue) PINPT_DOS%dos_ensurf(1:PINPT_DOS%dos_n_ensurf)
                        write(message,'(A,*(I5))')' DOS_EOUT: ENSURF OUT -> .TRUE.',PINPT_DOS%dos_ensurf(1:PINPT_DOS%dos_n_ensurf)  ; write_msg
                      elseif(i_dummy1 .ge. 1)then
                        i_dummy2 = nitems(dummy1)
                        allocate( strip_dummy(i_dummy2) )
                        read(dummy1,*,iostat=i_continue) (strip_dummy(i),i=1,i_dummy2)
                        ii = 0
                        do i = 1, i_dummy2
                          i_dummy3 = index(strip_dummy(i),':')
                          if(i_dummy3 .eq. 0) then
                            ii = ii + 1
                            call str2int(strip_dummy(i),i_dummy4)
                            i_dummyr(ii) = i_dummy4
                          elseif(i_dummy3 .gt. 1) then
                            ii = ii + 1
                            call strip_off(trim(strip_dummy(i)), dummy3, ' ', ':', 0)
                            call str2int(dummy3,i_dummy4)
                            call strip_off(trim(strip_dummy(i)), dummy3, ':', ' ', 2)
                            call str2int(dummy3,i_dummy5)
                            i_dummyr(ii:ii+i_dummy5 - i_dummy4) = (/ (k, k=i_dummy4, i_dummy5) /)
                            ii = ii + i_dummy5 - i_dummy4
                          endif
                        enddo
                        PINPT_DOS%dos_n_ensurf = ii
                        allocate( PINPT_DOS%dos_ensurf(PINPT_DOS%dos_n_ensurf) )
                        deallocate( strip_dummy )
                        PINPT_DOS%dos_ensurf(1:ii) = i_dummyr(1:ii)
                      endif
                    endif
                  endif
                case('PRINT_UNIT') ! kpoint unit : RECIPROCAL (fractional) or ANGSTROM (1/A)
                  read(inputline,*,iostat=i_continue) desc_str,desc_str
                  dummy = trim(desc_str)
                  PINPT_DOS%dos_kunit=dummy(1:1)
                  if(PINPT_DOS%dos_kunit .eq. 'r' .or. PINPT_DOS%dos_kunit .eq. 'R') then
                    PINPT_DOS%dos_kunit = 'R'
                    write(message,'(A)')' DOS_UNIT: KPOINT UNIT = Reciprocal'  ; write_msg
                  elseif(PINPT_DOS%dos_kunit .eq. 'a' .or. PINPT_DOS%dos_kunit .eq. 'A') then
                    PINPT_DOS%dos_kunit = 'A'
                    write(message,'(A)')' DOS_UNIT: KPOINT UNIT = Angstrom unit (1/A)'  ; write_msg
                  endif
                case('DOS_FNAME') ! DOS output file name
                  read(inputline,*,iostat=i_continue) desc_str,desc_str
                  PINPT_DOS%dos_filenm = trim(desc_str)
                  write(message,'(A,A)')' DOS_FNAM: ',trim(PINPT_DOS%dos_filenm)  ; write_msg
                case('LDOS_FNAME') ! LDOS output file name
                  read(inputline,*,iostat=i_continue) desc_str, desc_str
                  PINPT_DOS%ldos_filenm = trim(desc_str)
                  write(message,'(A,A)')'LDOS_FNAM: ',trim(PINPT_DOS%dos_filenm)  ; write_msg
                case('SMEARING') ! gaussian smearing
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_smearing
                  write(message,'(A,F8.4)')'DOS_SIGMA: GAUSSIAN WIDTH = ',PINPT_DOS%dos_smearing  ; write_msg
                case('DOS_SPARSE') ! logical flag for sparse matrix setup, if true, ERANGE will be 
                                   ! the energy window, and DOS_NRANGE will be used as NE_MAX
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_sparse
                  if(PINPT_DOS%dos_flag_sparse) then
                    write(message,'(A)')'DOS_SPARSE: .TRUE.'  ; write_msg
                  elseif(.not. PINPT_DOS%dos_flag_sparse) then
                    write(message,'(A)')'DOS_SPARSE: .FALSE.'  ; write_msg
                  endif

                case('PRINT_LDOS')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_ldos
                    if(PINPT_DOS%dos_flag_print_ldos) then
                      write(message,'(A)')'  !WARNING! DOS_LDOS -> .TRUE. but the list of target atoms is not specified.'  ; write_msg
                      write(message,'(A)')'  !WARNING!          -> LDOS for all atoms of the system will be evaluated by default.'  ; write_msg
                    endif
                  elseif(i_dummy .gt. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_ldos
                    read(inputline,*,iostat=i_continue) desc_str,dummy
                    if(PINPT_DOS%dos_flag_print_ldos) then
                      call strip_off (trim(inputline), dummy1, trim(dummy), ' ' , 2)   ! get dos_ensurf
                      i_dummy1=index(dummy1,':')
                      if(i_dummy1 .eq. 0) then
                        PINPT_DOS%dos_ldos_natom = nitems(dummy1)
                        allocate( PINPT_DOS%dos_ldos_atom(PINPT_DOS%dos_ldos_natom) )
                        read(dummy1,*,iostat=i_continue) PINPT_DOS%dos_ldos_atom(1:PINPT_DOS%dos_ldos_natom)
                        write(message,'(A,A)')' DOS_LDOS: .TRUE. , Atom_index = ',trim(dummy1)  ; write_msg
                      elseif(i_dummy1 .ge. 1)then
                        i_dummy2 = nitems(dummy1)
                        allocate( strip_dummy(i_dummy2) )
                        read(dummy1,*,iostat=i_continue) (strip_dummy(i),i=1,i_dummy2)
                        ii = 0
                        do i = 1, i_dummy2
                          i_dummy3 = index(strip_dummy(i),':')
                          if(i_dummy3 .eq. 0) then
                            ii = ii + 1
                            call str2int(strip_dummy(i),i_dummy4)
                            i_dummyr(ii) = i_dummy4
                          elseif(i_dummy3 .gt. 1) then
                            ii = ii + 1
                            call strip_off(trim(strip_dummy(i)), dummy3, ' ', ':', 0)
                            call str2int(dummy3,i_dummy4)
                            call strip_off(trim(strip_dummy(i)), dummy3, ':', ' ', 2)
                            call str2int(dummy3,i_dummy5)
                            i_dummyr(ii:ii+i_dummy5 - i_dummy4) = (/ (k, k=i_dummy4, i_dummy5) /)
                            ii = ii + i_dummy5 - i_dummy4
                          endif
                        enddo
                        PINPT_DOS%dos_ldos_natom = ii
                        allocate( PINPT_DOS%dos_ldos_atom(PINPT_DOS%dos_ldos_natom) )
                        deallocate( strip_dummy )
                        PINPT_DOS%dos_ldos_atom(1:ii) = i_dummyr(1:ii)
                        write(message,'(A,A)')' DOS_LDOS: .TRUE. , Atom_index = ',trim(dummy1)  ; write_msg
                      endif
                    endif
                  endif

              end select case_dos

            enddo set_dos
      return
   endsubroutine

   subroutine set_ribbon(PINPT, flag_kfile_ribbon, desc_str)
      type(incar)  ::  PINPT
      integer*4     i_continue
      integer*4     i_dummy
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      external      nitems
      character(*), parameter :: func = 'set_ribbon'
      logical       flag_kfile_ribbon

      PINPT%flag_set_ribbon = .true.
      write(message,'(A)')' SET_RIBN: SET UP RIBBON GEOMETRY = .TRUE.'  ; write_msg
      PINPT%ribbon_nslab(1:3) = -1
      PINPT%ribbon_vacuum(1:3) =  0   ! default setting for vacuum size

set_rib: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('NSLAB')
            read(inputline,*,iostat=i_continue) desc_str,PINPT%ribbon_nslab(1:3)

          case('VACUUM')
            read(inputline,*,iostat=i_continue) desc_str,PINPT%ribbon_vacuum(1:3)

          case('PRINT_ONLY_R')
            read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_only_ribbon_geom

          case('KFILE_R','KFILE_RIB','KFILE_RIBBON')
            read(inputline,*,iostat=i_continue) desc_str, PINPT%ribbon_kfilenm
            flag_kfile_ribbon = .true.
            PINPT%kfilenm = PINPT%ribbon_kfilenm
            write(message,'(A,A)')' KPTS_FNM:  (for ribbon) ',trim(PINPT%kfilenm)  ; write_msg

        end select
      enddo set_rib

      if(minval(PINPT%ribbon_nslab(1:3)) .le. 0) then
        write(message,'(A)')'  !WARNING! NSLAB tag for the RIBBON calculation is not properly defined. Please check the tag.'  ; write_msg
        write(message,'(A)')'  !WARNING! Proper usage: NSLAB   N1 N2 N3'  ; write_msg
        write(message,'(A)')'  !WARNING! Exit program...'  ; write_msg
        stop
      endif
      if(minval(PINPT%ribbon_nslab(1:3)) .lt. 0) then
        write(message,'(A)')'  !WARNING! VACUUM tag for the RIBBON calculation is not properly defined. Please check the tag.'  ; write_msg
        write(message,'(A)')'  !WARNING! Proper usage: VACUUM   vac_1 vac_2 vac_3'  ; write_msg
        write(message,'(A)')'  !WARNING! Exit program...'  ; write_msg
        stop
      endif

      write(message,'(A,3I5)')' RIB_SLAB: (N1*A1,N2*A2,N3*A3) => N1,N2,N3 =', PINPT%ribbon_nslab(1:3)  ; write_msg
      write(message,'(A,3F12.6)')' RIB_VACU: (in Angstrom)', PINPT%ribbon_vacuum(1:3)  ; write_msg
      if(PINPT%flag_print_only_ribbon_geom) then
        write(message,'(A,3F12.6)')' RIB_GEOM: PRINT_ONLY = .TRUE.'  ; write_msg
      else
        write(message,'(A,3F12.6)')' RIB_GEOM: PRINT_ONLY = .FALSE.'  ; write_msg
      endif

      return
   endsubroutine

   subroutine set_efield(PINPT, desc_str)
      type(incar)  ::  PINPT
      integer*4     i_continue
      integer*4     i_dummy
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      external      nitems
      character(*), parameter :: func = 'set_efield'

      PINPT%flag_efield = .true.

      do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('EFIELD')
            read(inputline,*,iostat=i_continue) desc_str,PINPT%efield(1:3)
            write(message,'(A,3F12.6)')'  E_FIELD:  ',PINPT%efield(1:3)  ; write_msg
          case('EF_CENTER','EF_ORIGIN','EF_ORIGIN_FRAC','EF_CENTER_FRAC') ! field_origin
            PINPT%flag_efield_frac = .true.
            PINPT%flag_efield_cart = .false.
            read(inputline,*,iostat=i_continue) desc_str,PINPT%efield_origin(1:3)
!           write(message,'(A,3F12.6)')'EF_ORIGIN:  (in factional coord) ',PINPT%efield_origin(1:3)  ; write_msg
          case('EF_CENTER_CART','EF_ORIGIN_CART','EF_CORIGIN') ! field_origin
            PINPT%flag_efield_frac = .false.
            PINPT%flag_efield_cart = .true.
            read(inputline,*,iostat=i_continue) desc_str,PINPT%efield_origin_cart(1:3)
!           write(message,'(A,3F12.6)')'EF_ORIGIN:  (in cartesian coord) ',PINPT%efield_origin_cart(1:3)  ; write_msg
        end select
      enddo 

      return
   endsubroutine

   subroutine set_weight_factor(PINPT, PWGHT, desc_str)
      type(incar)  ::  PINPT
      type(weight)  :: PWGHT
      integer*4     i, i_orb, i_deg
      integer*4     i_continue
      integer*4     i_dummy
      integer*4     nitems
      character*132 inputline, inputline_dummy
      character*40  desc_str
      character*40  str2lowcase
      external      nitems, str2lowcase
      character(*), parameter :: func = 'set_weight_factor'
      character*132 strip_kp_(max_set_weight),strip_tb_(max_set_weight),strip_df_(max_set_weight),strip_wt_(max_set_weight)
      character*132 strip_kp_orb_(max_set_weight),strip_tb_orb_(max_set_weight),strip_orb_(max_set_weight),strip_site_(max_set_weight)
      character*132 strip_pen_orb_(max_set_weight)
      character*132 strip_kp_deg_(max_set_weight), strip_tb_deg_(max_set_weight),strip_df_deg_(max_set_weight),strip_wt_deg_(max_set_weight)

      i=0
      i_orb = 0
      i_deg = 0
 wgt: do while(trim(desc_str) .ne. 'END' )

        if(PINPT%flag_use_weight) then
          read(pid_param,'(A)',iostat=i_continue) inputline
          call take_comment(trim(inputline), inputline_dummy) ! take comment line (get string from comment mark )
          inputline = inputline_dummy
          if(i_continue .lt. 0) exit ! exit loop if end-of-file
        else
          read(pid_incar,'(A)',iostat=i_continue) inputline
          if(i_continue .lt. 0) exit ! exit loop if end-of-file
        endif

        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(index(inputline,'#') .gt. 1) then
          call strip_off (trim(inputline), inputline_dummy, '', '#', 0) ! strip-off '#' comments
          inputline = inputline_dummy
        endif
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        if( index(inputline,'ORBT_I') .eq. 0 .and. index(inputline,'KRANGE') .gt. 0 .and. index(inputline,'DEGENW') .eq. 0 .and. &
                                                                                          index(inputline,'WEIGHT') .gt. 0) then
          i=i+1
          call strip_off (trim(inputline), strip_kp_(i), 'KRANGE', 'TBABND', 1)   ! get KRANGE strip
          call strip_off (trim(inputline), strip_tb_(i), 'TBABND', 'DFTBND', 1)   ! get TBABND strip
          call strip_off (trim(inputline), strip_df_(i), 'DFTBND', 'WEIGHT', 1)   ! get DFTBND strip
          call strip_off (trim(inputline), strip_wt_(i), 'WEIGHT', ' '     , 2)   ! get WEIGHT strip
        elseif( index(inputline,'ORBT_I') .gt. 0 .and. index(inputline,'KRANGE') .gt. 0 .and. index(inputline,'DEGENW') .eq. 0) then
          i_orb = i_orb + 1
          call strip_off (trim(inputline), strip_kp_orb_(i_orb), 'KRANGE', 'TBABND', 1)   ! get KRANGE strip
          call strip_off (trim(inputline), strip_tb_orb_(i_orb), 'TBABND', 'ORBT_I', 1)   ! get TBABND strip
          call strip_off (trim(inputline), strip_orb_(i_orb),    'ORBT_I', 'SITE_I', 1)   ! get ORBT_I strip
          call strip_off (trim(inputline), strip_site_(i_orb),   'SITE_I', 'PENALTY', 1)   ! get SITE_I strip
          strip_site_(i_orb) = str2lowcase(strip_site_(i_orb))
          call strip_off (trim(inputline), strip_pen_orb_(i_orb), 'PENALTY', ' '     , 2)   ! get PENALTY strip
        elseif( index(inputline,'ORBT_I') .eq. 0 .and. index(inputline,'KRANGE') .gt. 0 .and. index(inputline,'DEGENW') .gt. 0) then
          i_deg = i_deg + 1
          call strip_off (trim(inputline), strip_kp_deg_(i_deg), 'KRANGE', 'TBABND', 1)   ! get KRANGE strip
          call strip_off (trim(inputline), strip_tb_deg_(i_deg), 'TBABND', 'DFTBND', 1)   ! get TBABND strip
          call strip_off (trim(inputline), strip_df_deg_(i_deg), 'DFTBND', 'DEGENW', 1)   ! get DFTBND strip
          call strip_off (trim(inputline), strip_wt_deg_(i_deg), 'DEGENW', ' '     , 2)   ! get DEGENERACY WEIGHT strip
        endif
      enddo wgt

      if(nitems(inputline) .eq. 3) then
        read(inputline,*,iostat=i_continue) desc_str, desc_str, desc_str
        if(desc_str .eq. 'PRINT_ONLY') PINPT%flag_print_only_target = .true.
      endif

      PWGHT%nweight=i
      PINPT%nweight=i
      PWGHT%npenalty_orb = i_orb
      PINPT%npenalty_orb = i_orb
      PWGHT%ndegenw=i_deg
      PINPT%ndegenw=i_deg

      ! Set WEIGHT
      if(PWGHT%nweight .eq. 0) then
        PWGHT%flag_weight_default = .true.
      elseif( PWGHT%nweight .ge. 1) then
        PWGHT%flag_weight_default = .false.
        write(message,'(A,I8)')' N_CONSTR:',PWGHT%nweight  ; write_msg
        if(allocated(PINPT%strip_kp) ) deallocate(PINPT%strip_kp) ; allocate( PINPT%strip_kp(PWGHT%nweight) )
        if(allocated(PINPT%strip_tb) ) deallocate(PINPT%strip_tb) ; allocate( PINPT%strip_tb(PWGHT%nweight) )
        if(allocated(PINPT%strip_df) ) deallocate(PINPT%strip_df) ; allocate( PINPT%strip_df(PWGHT%nweight) )
        if(allocated(PINPT%strip_wt) ) deallocate(PINPT%strip_wt) ; allocate( PINPT%strip_wt(PWGHT%nweight) )
        do i = 1, PWGHT%nweight
          PINPT%strip_kp(i)=strip_kp_(i)
          PINPT%strip_tb(i)=strip_tb_(i)
          PINPT%strip_df(i)=strip_df_(i)
          PINPT%strip_wt(i)=strip_wt_(i)
        enddo
      endif

      ! Set DEGENERACY WEIGHT
      if(PWGHT%ndegenw .eq. 0) then
        PINPT%flag_fit_degeneracy = .false.
      elseif( PWGHT%ndegenw .ge. 1) then
        PINPT%flag_fit_degeneracy = .true. 
        write(message,'(A,I8)')' N_DEGENW:',PWGHT%ndegenw  ; write_msg
        if(allocated(PINPT%strip_kp_deg) ) deallocate(PINPT%strip_kp_deg) ; allocate( PINPT%strip_kp_deg(PWGHT%ndegenw) )
        if(allocated(PINPT%strip_tb_deg) ) deallocate(PINPT%strip_tb_deg) ; allocate( PINPT%strip_tb_deg(PWGHT%ndegenw) )
        if(allocated(PINPT%strip_df_deg) ) deallocate(PINPT%strip_df_deg) ; allocate( PINPT%strip_df_deg(PWGHT%ndegenw) )
        if(allocated(PINPT%strip_wt_deg) ) deallocate(PINPT%strip_wt_deg) ; allocate( PINPT%strip_wt_deg(PWGHT%ndegenw) )
        do i = 1, PWGHT%ndegenw
          PINPT%strip_kp_deg(i)=strip_kp_deg_(i)
          PINPT%strip_tb_deg(i)=strip_tb_deg_(i)
          PINPT%strip_df_deg(i)=strip_df_deg_(i)
          PINPT%strip_wt_deg(i)=strip_wt_deg_(i)
        enddo
      endif


      ! Set ORBITAL PENALTY 
      if(PWGHT%npenalty_orb .eq. 0) then
        PWGHT%flag_weight_default_orb = .true.
      elseif( PWGHT%npenalty_orb .ge. 1) then
        PWGHT%flag_weight_default_orb = .false.
        write(message,'(A,I8,A)')' N_CONSTR:',PWGHT%npenalty_orb,' (# of penalty weight for orbital)'  ; write_msg
        if(allocated(PINPT%strip_kp_orb)) deallocate(PINPT%strip_kp_orb) ; allocate( PINPT%strip_kp_orb(PWGHT%npenalty_orb) )
        if(allocated(PINPT%strip_tb_orb)) deallocate(PINPT%strip_tb_orb) ; allocate( PINPT%strip_tb_orb(PWGHT%npenalty_orb) )
        if(allocated(PINPT%strip_pen_orb))deallocate(PINPT%strip_pen_orb); allocate( PINPT%strip_pen_orb(PWGHT%npenalty_orb) )
        if(allocated(PINPT%strip_orb))    deallocate(PINPT%strip_orb)    ; allocate( PINPT%strip_orb(PWGHT%npenalty_orb) )
        if(allocated(PINPT%strip_site))   deallocate(PINPT%strip_site)   ; allocate( PINPT%strip_site(PWGHT%npenalty_orb) )
        do i = 1, PWGHT%npenalty_orb
          PINPT%strip_kp_orb(i)=strip_kp_orb_(i)
          PINPT%strip_tb_orb(i)=strip_tb_orb_(i)
          PINPT%strip_pen_orb(i)=strip_pen_orb_(i)
          PINPT%strip_orb(i)=strip_orb_(i)
          PINPT%strip_site(i)=strip_site_(i)
        enddo

      endif

      return
   endsubroutine

   subroutine set_z2(PINPT, PINPT_BERRY, desc_str)
      type(incar)  ::  PINPT
      type(berry)  ::  PINPT_BERRY
      integer*4     i_continue
      integer*4     i_dummy, i_dummy2
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      character*40  dummy, dummy2
      external      nitems
      character(*), parameter :: func = 'set_z2'

      PINPT%flag_get_z2 = .true.
      PINPT%flag_berry = .true.
      PINPT_BERRY%z2_nkdiv = 11 ! default
      PINPT_BERRY%z2_nkpath=  0 ! default
      PINPT_BERRY%z2_filenm= 'Z2.WCC' ! default
      PINPT_BERRY%z2_gap_filenm= 'Z2.GAP' ! default
      PINPT_BERRY%z2_dimension = 3 ! default
      allocate(PINPT_BERRY%z2_axis(3))
      PINPT_BERRY%z2_axis(:) = (/1,2,3/) ! default
      PINPT_BERRY%flag_z2_get_chern = .false. ! default

      write(message,'(A)')'   GET_Z2: .TRUE. '  ; write_msg
  set_wann: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('Z2_ERANGE')  ! get z2_range
            call strip_off (trim(inputline), PINPT_BERRY%strip_z2_range, 'Z2_ERANGE', ' ' , 2)
            write(message,'(3A)')'  Z2_RANG:  ',adjustl(trim(PINPT_BERRY%strip_z2_range)),' (ENERGY RANGE TO BE CALCULATED)'  ; write_msg

          case('Z2_DIMENSION')  ! which direction? 
            call strip_off (trim(inputline), dummy, 'Z2_DIMENSION', ' ' , 2)
            i_dummy=index(dummy,'D')
            if(i_dummy .eq. 0) then
              write(message,'(A)')' !WARN!: Z2_DIMENSION should be start with one of following strings:'  ; write_msg
              write(message,'(A)')'  3D  => Z2_DIMENSION 3D'  ; write_msg
              write(message,'(A)')'  2D  => Z2_DIMENSION 2D:z # z is normal direction of surface'  ; write_msg
              write(message,'(A)')'  1D  => Z2_DIMENSION 1D:x # x is parallel direction of 1D system'  ; write_msg
              write(message,'(A)')'  For 2D and 1D, "x", "y", or "z" indicates the normal and parallel'         ; write_msg
              write(message,'(A)')'  direction of the system'  ; write_msg
              write(message,'(A)')'  Exit program anyway... Check your "SET Z2" in "INCAR-TB" again.'  ; write_msg
              stop
            endif

            call str2int(dummy(1:i_dummy-1),PINPT_BERRY%z2_dimension)

            if(PINPT_BERRY%z2_dimension .eq. 3) then
              deallocate(PINPT_BERRY%z2_axis)
              allocate(PINPT_BERRY%z2_axis(3))
              PINPT_BERRY%z2_axis(:) = (/1,2,3/)
              write(message,'(3A)')'   Z2_DIM: ', trim(dummy), ' mode (three direction will be checked)'  ; write_msg

            elseif(PINPT_BERRY%z2_dimension .eq. 2) then
              deallocate(PINPT_BERRY%z2_axis)
              allocate(PINPT_BERRY%z2_axis(1))
              i_dummy = index(dummy,':')
              if(i_dummy .eq. 0) then
                PINPT_BERRY%z2_axis(1) = 3 ! default : third reciprocal axis will be chosen for the WCC surface normal axis
                write(message,'(3A)')'   Z2_DIM: ', trim(dummy), ' mode (three direction will be checked)'  ; write_msg
                write(message,'( A)')'  Z2_AXIS: B3 direction of the reciprocal lattice vector (default)'  ; write_msg
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','b1', 'X', 'A', 'B1')
                    PINPT_BERRY%z2_axis(1) = 1
                    write(message,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode'  ; write_msg
                    write(message,'( A)')'  Z2_AXIS: B1 direction of the reciprocal lattice vector (B1-B2 plane, xy)'  ; write_msg
                  case('2','y','b','b2', 'Y', 'B', 'B2')
                    PINPT_BERRY%z2_axis(1) = 2
                    write(message,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode'  ; write_msg
                    write(message,'( A)')'  Z2_AXIS: B2 direction of the reciprocal lattice vector (B2-B3 plane, yz)'  ; write_msg
                  case('3','z','c','b3', 'Z', 'C', 'B3')
                    PINPT_BERRY%z2_axis(1) = 3
                    write(message,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode'  ; write_msg
                    write(message,'( A)')'  Z2_AXIS: B3 direction of the reciprocal lattice vector (B3-B1 plane, zx)'  ; write_msg
                endselect
              endif

            elseif(PINPT_BERRY%z2_dimension .eq. 1) then
              deallocate(PINPT_BERRY%z2_axis)
              allocate(PINPT_BERRY%z2_axis(1))
              i_dummy = index(dummy, ':')
              if(i_dummy .eq. 0) then    
                PINPT_BERRY%z2_axis(1) = 1
                write(message,'(3A)')'   Z2_DIM: ', trim(dummy), ' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                write(message,'( A)')'  Z2_AXIS: B1 direction of the reciprocal lattice vector (default, parallel to B1)'  ; write_msg
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','kx', 'b1', 'X', 'A', 'B1', 'KX')
                    PINPT_BERRY%z2_axis(1) = 1
                    write(message,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                    write(message,'( A)')'  Z2_AXIS: B1 direction of the reciprocal lattice vector (parallel to B1)'  ; write_msg
                  case('2','y','b','ky', 'b2', 'Y', 'B', 'B2', 'KY')
                    PINPT_BERRY%z2_axis(1) = 2
                    write(message,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                    write(message,'( A)')'  Z2_AXIS: B2 direction of the reciprocal lattice vector (parallel to B2)'  ; write_msg
                  case('3','z','c','kz', 'b3', 'Z', 'C', 'B3', 'KZ')
                    PINPT_BERRY%z2_axis(1) = 3
                    write(message,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                    write(message,'( A)')'  Z2_AXIS: B3 direction of the reciprocal lattice vector (parallel to B3)'  ; write_msg
                endselect

              endif
            endif

          case('Z2_NKDIV')
            i_dummy = nitems(inputline) -1
            if(i_dummy .eq. 1) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%z2_nkdiv
            elseif(i_dummy .eq. 2) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%z2_nkdiv,PINPT_BERRY%z2_nkpath
            endif
    
          case('GET_CHERN','Z2_CHERN')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_z2_get_chern
            if(PINPT_BERRY%flag_z2_get_chern) then
              write(message,'(3A)')'GET_CHERN: (Z2 SET) .TRUE. (evaluate Chern number for the bands with ERANGE)'  ; write_msg
            else
              write(message,'(3A)')'GET_CHERN: (Z2 SET) .FALSE.'  ; write_msg
            endif

          case('SET_PHASE', 'GET_PHASE', 'PHASE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_z2_phase

        end select

      enddo set_wann

      if(PINPT_BERRY%z2_dimension .eq. 3) then
        if(PINPT_BERRY%z2_nkpath .eq. 0) then
          write(message,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv,' (division along k-line for wcc))'  ; write_msg
          write(message,'(A)')' WARN: The direction of WCC evolution did not provided,'  ; write_msg
          write(message,'(A)')'       so we will setup with default =>  Z2_NKPATH = 21'  ; write_msg
          PINPT_BERRY%z2_nkpath = 21
        elseif(PINPT_BERRY%z2_nkpath .gt. 0) then
          write(message,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv, ' (division along k-line for wcc))'  ; write_msg
          write(message,'(A, I4,A)')' Z2_NKPATH:',PINPT_BERRY%z2_nkpath,' (division along k-direction for wcc evolution))'  ; write_msg
        endif
        PINPT_BERRY%z2_nplane = 2
      elseif(PINPT_BERRY%z2_dimension .eq. 2) then
        if(PINPT_BERRY%z2_nkpath .eq. 0) then
          write(message,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv,' (division along k-line for wcc))'  ; write_msg
          write(message,'(A)')' WARN: The direction of WCC evolution did not provided,'  ; write_msg
          write(message,'(A)')'       so we will setup with default =>  Z2_NKPATH = 21'  ; write_msg
          PINPT_BERRY%z2_nkpath = 21
        elseif(PINPT_BERRY%z2_nkpath .gt. 0) then
          write(message,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv, ' (division along k-line for wcc))'  ; write_msg
          write(message,'(A, I4,A)')' Z2_NKPATH:',PINPT_BERRY%z2_nkpath,' (division along k-direction for wcc evolution))'  ; write_msg
        endif
        PINPT_BERRY%z2_nplane = 1
      elseif(PINPT_BERRY%z2_dimension .eq. 1) then
        if(PINPT_BERRY%z2_nkpath .eq. 0) then
          write(message,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv,' (division along k-line for wcc))'  ; write_msg
          PINPT_BERRY%z2_nkpath = 1
        elseif(PINPT_BERRY%z2_nkpath .gt. 0) then
          write(message,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv, ' (division along k-line for wcc))'  ; write_msg
          PINPT_BERRY%z2_nkdiv = 1 ! for 1D case, it is forced to be 1. There is no reason to scan other direction.
        endif
        PINPT_BERRY%z2_nplane = 1
      endif

      allocate(PINPT_BERRY%z2_kpoint     (3,PINPT_BERRY%z2_nkdiv,PINPT_BERRY%z2_nkpath, &
                                            PINPT_BERRY%z2_nplane,size(PINPT_BERRY%z2_axis)) )
      allocate(PINPT_BERRY%z2_kpoint_reci(3,PINPT_BERRY%z2_nkdiv,PINPT_BERRY%z2_nkpath, &
                                            PINPT_BERRY%z2_nplane,size(PINPT_BERRY%z2_axis)) )

      return
   endsubroutine
   subroutine set_wcc(PINPT, PINPT_BERRY, desc_str)
      type(incar)  ::  PINPT
      type(berry)  ::  PINPT_BERRY
      integer*4     i_continue
      integer*4     i_dummy
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      external      nitems
      character(*), parameter :: func = 'set_wcc'
      real*8        wcc_kpath_dummy(3,2,max_kpoint)
      character*2   k_direct(3)
      integer*4     ikpath, ik, i

      k_direct(1)='B1'
      k_direct(2)='B2'
      k_direct(3)='B3'
      PINPT%flag_get_wcc = .true.
      PINPT%flag_berry = .true.
      PINPT_BERRY%wcc_nkdiv = 11 ! default
      PINPT_BERRY%wcc_nkdiv2= 21 ! default
      PINPT_BERRY%wcc_filenm= 'WCC.OUT.dat' ! default
      PINPT_BERRY%wcc_gap_filenm= 'WCC.GAP.dat' ! default
      PINPT_BERRY%wcc_kpath_shift = (/0d0,0d0,0d0/) ! default
      PINPT_BERRY%flag_wcc_get_chern = .false. ! default
      PINPT_BERRY%flag_wcc_get_chern_spin = .false. ! default
      ikpath = 0
      ik = 0

      write(message,'(A)')'  GET_WCC: .TRUE. (GET Wannier charge center)'  ; write_msg
  set_wann: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('WCC_ERANGE')  ! get wcc_range
            call strip_off (trim(inputline), PINPT_BERRY%strip_wcc_range, 'WCC_ERANGE', ' ' , 2)
            write(message,'(3A)')' WCC_RANG:  ',adjustl(trim(PINPT_BERRY%strip_wcc_range)),' (ENERGY RANGE TO BE CALCULATED)'  ; write_msg
          case('WCC_FNAME') ! WCC output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%wcc_filenm = trim(desc_str)
            write(message,'(A,A)')' WCC_FNAM:  ',trim(PINPT_BERRY%wcc_filenm)  ; write_msg

          case('WCC_GAP_FNAME','WCC_FNAME_GAP') ! WCC largest gap output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%wcc_gap_filenm = trim(desc_str)
            write(message,'(A,A)')' WCC_GAPF:  (largest gap file)',trim(PINPT_BERRY%wcc_gap_filenm)  ; write_msg

          ! this option is only for test purpose...
          case('WCC_PATH_SHIFT','WCC_KPATH_SHIFT')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%wcc_kpath_shift(1:3)

          case('WCC_PATH','WCC_KPATH')
            ikpath = ikpath + 1
            read(inputline,*,iostat=i_continue) desc_str, wcc_kpath_dummy(1:3,1,ikpath), wcc_kpath_dummy(1:3,2,ikpath)
            write(message,'(A,3F12.8,A,3F12.8,A)')' WCC_PATH: (KPATH) ', wcc_kpath_dummy(1:3,1,ikpath),' --> ', wcc_kpath_dummy(1:3,2,ikpath),' (closed loop for WCC, reci unit)'  ; write_msg

          case('WCC_DIREC','WCC_DIRECT','WCC_DIRECTION')
            PINPT_BERRY%flag_wcc_evolve = .true.
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%wcc_direction
            write(message,'(3A)')' WCC_DIR :  ',k_direct(PINPT_BERRY%wcc_direction),' (k-direction of WCC evolution)'  ; write_msg

          case('WCC_NKDIV')
            i_dummy = nitems(inputline) -1
            if(i_dummy .eq. 1) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%wcc_nkdiv
            elseif(i_dummy .eq. 2) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%wcc_nkdiv,PINPT_BERRY%wcc_nkdiv2
            endif

          case('GET_CHERN','WCC_CHERN')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_wcc_get_chern
            if(PINPT_BERRY%flag_wcc_get_chern) then
              write(message,'(3A)')'GET_CHERN: (WCC SET) .TRUE. (Chern number for the bands with ERANGE)'  ; write_msg
            else
              write(message,'(3A)')'GET_CHERN: (WCC SET) .FALSE.'  ; write_msg
            endif

          case('GET_CHERN_SPIN','WCC_CHERN_SPIN', 'GET_SPIN_CHERN', 'WCC_SPIN_CHERN')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_wcc_get_chern_spin
            if(PINPT_BERRY%flag_wcc_get_chern_spin) then
              write(message,'(3A)')'GET_CHERN:  SPIN_CHERN .TRUE. (WCC SET) (spin Chern number for the bands with ERANGE)'  ; write_msg
            else
              write(message,'(3A)')'GET_CHERN:  SPIN_CHERN .FALSE. (WCC SET)'  ; write_msg
            endif

          case('SET_PHASE', 'GET_PHASE', 'PHASE', 'WCC_PHASE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_wcc_phase

        end select

      enddo set_wann

      if(PINPT_BERRY%flag_wcc_evolve .and. PINPT_BERRY%wcc_nkdiv2 .eq. 1) then
        PINPT_BERRY%wcc_nkdiv2 = 31 ! set default if wcc_evolution along WCC_DIRECT is requested
      endif

      if(PINPT_BERRY%wcc_nkdiv2 .gt. 1) then
        write(message,'(A,2I4,3A)')' WCC_KDIV: ',PINPT_BERRY%wcc_nkdiv, PINPT_BERRY%wcc_nkdiv2,' (division along k-line(wcc) and k-direct(evolution along ', k_direct(PINPT_BERRY%wcc_direction),')'  ; write_msg
      elseif(PINPT_BERRY%wcc_nkdiv2 .eq. 1) then
        write(message,'(A, I4,A)')' WCC_KDIV: ',PINPT_BERRY%wcc_nkdiv, ' (division along k-path(wcc))'  ; write_msg
        write(message,'(A)')' WARN: The direction of WCC evolution did not provided,'  ; write_msg
        write(message,'(A)')'       so we will setup with default => WCC_DIRECT = 2 ,'  ; write_msg
        write(message,'(A)')'                                        WCC_KDIV2  = 31'  ; write_msg
        PINPT_BERRY%wcc_nkdiv2 = 31
        PINPT_BERRY%flag_wcc_evolve = .true.
      endif

      allocate(PINPT_BERRY%wcc_kpath(3,2,PINPT_BERRY%wcc_nkdiv2))
      allocate(PINPT_BERRY%wcc_kpoint(3,PINPT_BERRY%wcc_nkdiv,PINPT_BERRY%wcc_nkdiv2))
      allocate(PINPT_BERRY%wcc_kpoint_reci(3,PINPT_BERRY%wcc_nkdiv,PINPT_BERRY%wcc_nkdiv2))
      do ik = 1, PINPT_BERRY%wcc_nkdiv2
        do i = 1, 3
          if(ik .eq. 1) then
              PINPT_BERRY%wcc_kpath(i,:,1) = wcc_kpath_dummy(i,:,1) + PINPT_BERRY%wcc_kpath_shift(i)
          elseif(ik .gt. 1 .and. i .eq. PINPT_BERRY%wcc_direction) then
              PINPT_BERRY%wcc_kpath(i,1,ik) = wcc_kpath_dummy(i,1,1) + 1d0/(PINPT_BERRY%wcc_nkdiv2-1)*(ik-1)+ PINPT_BERRY%wcc_kpath_shift(i)
              PINPT_BERRY%wcc_kpath(i,2,ik) = wcc_kpath_dummy(i,2,1) + 1d0/(PINPT_BERRY%wcc_nkdiv2-1)*(ik-1)+ PINPT_BERRY%wcc_kpath_shift(i)
          elseif(ik .gt. 1 .and. i .ne. PINPT_BERRY%wcc_direction) then
              PINPT_BERRY%wcc_kpath(i,1,ik) = wcc_kpath_dummy(i,1,1)+ PINPT_BERRY%wcc_kpath_shift(i)
              PINPT_BERRY%wcc_kpath(i,2,ik) = wcc_kpath_dummy(i,2,1)+ PINPT_BERRY%wcc_kpath_shift(i)
          endif
        enddo
      enddo
      PINPT_BERRY%wcc_nkpath     = PINPT_BERRY%wcc_nkdiv2

      return
   endsubroutine

   subroutine set_zak_phase(PINPT, PINPT_BERRY, desc_str)
      type(incar)  ::  PINPT
      type(berry)  ::  PINPT_BERRY
      integer*4     i_continue
      integer*4     i_dummy
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      external      nitems
      character(*), parameter :: func = 'set_zak_phase'
      real*8        zak_kpath_dummy(3,2,max_kpoint)
      character*2   k_direct(3)
      integer*4     ikpath, ik, i

      k_direct(1)='B1'
      k_direct(2)='B2'
      k_direct(3)='B3'
      PINPT%flag_get_zak_phase = .true.
      PINPT%flag_berry = .true.
      PINPT_BERRY%zak_nkdiv = 11 ! default
      PINPT_BERRY%zak_nkdiv2=  1 ! default
      PINPT_BERRY%zak_kpath_shift = (/0d0,0d0,0d0/) ! default
      PINPT_BERRY%zak_filenm='ZAK_PHASE.OUT.dat'
      ikpath = 0
      ik = 0

      write(message,'(A)')' ZAK_PHAS: .TRUE. (GET ZAK PHASE)'  ; write_msg
   set_zak: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('ZAK_ERANGE')  ! get zak_range
            call strip_off (trim(inputline), PINPT_BERRY%strip_zak_range, 'ZAK_ERANGE', ' ' , 2) 
            write(message,'(3A)')' ZAK_RANG:  ',adjustl(trim(PINPT_BERRY%strip_zak_range)), ' (ENERGY RANGE TO BE CALCULATED)'  ; write_msg

          case('ZAK_FNAME') ! ZAK phase output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%zak_filenm = trim(desc_str)
            write(message,'(A,A)')' ZAK_FNAM:  ',trim(PINPT_BERRY%zak_filenm)  ; write_msg

          case('ZAK_SEPARATE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_zak_separate
            if(PINPT%flag_zak_separate) then
              write(message,'(A)')' ZAK_SEPR: .TRUE. (also separate into each eigenstate)'  ; write_msg
            elseif(.not. PINPT%flag_zak_separate) then
              write(message,'(A)')' ZAK_SEPR: .FALSE. (do not separate into each eigenstate)'  ; write_msg
            endif

          ! this option is only for test purpose...
          case('ZAK_PATH_SHIFT','ZAK_KPATH_SHIFT')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%zak_kpath_shift(1:3)

          case('ZAK_PATH','ZAK_KPATH')
            ikpath = ikpath + 1
            read(inputline,*,iostat=i_continue) desc_str, zak_kpath_dummy(1:3,1,ikpath), zak_kpath_dummy(1:3,2,ikpath)
            write(message,'(A,3F12.8,A,3F12.8,A)')' ZAK_PATH: (KPATH) ', zak_kpath_dummy(1:3,1,ikpath),' --> ', zak_kpath_dummy(1:3,2,ikpath),' (closed loop for Zak phase, reci unit)'  ; write_msg
          case('ZAK_DIREC','ZAK_DIRECT','ZAK_DIRECTION')
            PINPT_BERRY%flag_zak_evolve = .true.
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%zak_direction
            write(message,'(3A)')' ZAK_DIR :  ',k_direct(PINPT_BERRY%zak_direction),' (k-direction of Zak phase evolution)'  ; write_msg

          case('ZAK_NKDIV')
            i_dummy = nitems(inputline) -1
            if(i_dummy .eq. 1) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%zak_nkdiv
            elseif(i_dummy .eq. 2) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%zak_nkdiv,PINPT_BERRY%zak_nkdiv2
            endif

          case('SET_PHASE','GET_PHASE', 'PHASE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_zak_phase

        end select
       
      enddo set_zak

      if(PINPT_BERRY%flag_zak_evolve .and. PINPT_BERRY%zak_nkdiv2 .eq. 1) then
        PINPT_BERRY%zak_nkdiv2 = 31 ! set default if zak_evolution along ZAK_DIRECT is requested
      endif

      if(PINPT_BERRY%zak_nkdiv2 .gt. 1) then
        write(message,'(A,2I4,3A)')' ZAK_KDIV: ',PINPT_BERRY%zak_nkdiv,PINPT_BERRY%zak_nkdiv2,' (division along k-path(zak_phase) and k-direct(evolution along ', k_direct(PINPT_BERRY%zak_direction),')'  ; write_msg
      elseif(PINPT_BERRY%zak_nkdiv2 .eq. 1) then
        write(message,'(A, I4,A)')' ZAK_KDIV: ',PINPT_BERRY%zak_nkdiv,' (division along k-path(zak_phase))'  ; write_msg
      endif

      if(PINPT_BERRY%zak_nkdiv2 .gt. 1) then
        allocate(PINPT_BERRY%zak_kpath(3,2,PINPT_BERRY%zak_nkdiv2))
        allocate(PINPT_BERRY%zak_kpoint(3,PINPT_BERRY%zak_nkdiv,PINPT_BERRY%zak_nkdiv2))
        allocate(PINPT_BERRY%zak_kpoint_reci(3,PINPT_BERRY%zak_nkdiv,PINPT_BERRY%zak_nkdiv2))
        do ik = 1, PINPT_BERRY%zak_nkdiv2
          do i = 1, 3
            if(ik .eq. 1) then
                PINPT_BERRY%zak_kpath(i,:,1)  = zak_kpath_dummy(i,:,1)+ PINPT_BERRY%zak_kpath_shift(i)
            elseif(ik .gt. 1 .and. i .eq. PINPT_BERRY%zak_direction) then
                PINPT_BERRY%zak_kpath(i,1,ik) = zak_kpath_dummy(i,1,1) + 1d0/(PINPT_BERRY%zak_nkdiv2-1)*(ik-1)+ PINPT_BERRY%zak_kpath_shift(i)
                PINPT_BERRY%zak_kpath(i,2,ik) = zak_kpath_dummy(i,2,1) + 1d0/(PINPT_BERRY%zak_nkdiv2-1)*(ik-1)+ PINPT_BERRY%zak_kpath_shift(i)
            elseif(ik .gt. 1 .and. i .ne. PINPT_BERRY%zak_direction) then
                PINPT_BERRY%zak_kpath(i,1,ik) = zak_kpath_dummy(i,1,1) + PINPT_BERRY%zak_kpath_shift(i)
                PINPT_BERRY%zak_kpath(i,2,ik) = zak_kpath_dummy(i,2,1) + PINPT_BERRY%zak_kpath_shift(i)
            endif
          enddo
        enddo
        PINPT_BERRY%zak_nkpath     = PINPT_BERRY%zak_nkdiv2

      elseif(PINPT_BERRY%zak_nkdiv2 .eq. 1) then
        if(ikpath .ge. 1) then
          allocate(PINPT_BERRY%zak_kpath(3,2,ikpath))
          allocate(PINPT_BERRY%zak_kpoint(3,PINPT_BERRY%zak_nkdiv,ikpath))
          allocate(PINPT_BERRY%zak_kpoint_reci(3,PINPT_BERRY%zak_nkdiv,ikpath))

          do ik = 1, ikpath
            do i = 1, PINPT_BERRY%zak_nkdiv
              PINPT_BERRY%zak_kpath(:,i,ik) = zak_kpath_dummy(:,i,ik) + PINPT_BERRY%zak_kpath_shift(:)
            enddo
          enddo
          PINPT_BERRY%zak_nkpath     = ikpath

        elseif(ikpath .eq. 0) then
          allocate(PINPT_BERRY%zak_kpath(3,2,1))
          allocate(PINPT_BERRY%zak_kpoint(3,PINPT_BERRY%zak_nkdiv,1))
          allocate(PINPT_BERRY%zak_kpoint_reci(3,PINPT_BERRY%zak_nkdiv,1))
          PINPT_BERRY%zak_nkpath         = 1
          PINPT_BERRY%zak_kpath(1:3,1,1) = (/0.0d0,0d0,0d0/)+ PINPT_BERRY%zak_kpath_shift(:)
          PINPT_BERRY%zak_kpath(1:3,2,1) = (/1.0d0,0d0,0d0/)+ PINPT_BERRY%zak_kpath_shift(:)
         write(message,'(A,3F12.8,A,3F12.8)')' ZAK_PATH:  (KPATH,default) ',PINPT_BERRY%zak_kpath(1:3,1,1),'  --> ',PINPT_BERRY%zak_kpath(1:3,2,1)  ; write_msg
        endif

      endif

      return
   endsubroutine

   subroutine set_berry_curvature(PINPT, PINPT_BERRY, desc_str)
      type(incar)  ::  PINPT
      type(berry)  ::  PINPT_BERRY
      integer*4     i_continue
      integer*4     i_dummy
      integer*4     nitems
      character*132 inputline
      character*40  desc_str
      character*40  dummy, dummy2
      external      nitems
      character(*), parameter :: func = 'set_berry_curvature'

      PINPT%flag_get_berry_curvature = .true.
      PINPT%flag_berry = .true.
      PINPT_BERRY%strip_bc_range=':' ! default, calculate for all eigenstates
      PINPT_BERRY%flag_bc_filenm_provided = .false.
      PINPT_BERRY%flag_bc_method_kubo = .true. ! default
      PINPT_BERRY%bc_dimension = 2  ! default
      allocate(PINPT_BERRY%bc_axis(1)); PINPT_BERRY%bc_axis(1) = 3 ! default
      PINPT_BERRY%bc_nkdiv(1:2)= 21 ! default
      PINPT_BERRY%bc_nkdiv(3)  = 1  ! default

      write(message,'(A)')' '  ; write_msg
      write(message,'(A)')' BERRYCRV: .TRUE. (GET BERRY CURVATURE)'  ; write_msg

 set_berryc: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('BERRYC_ERANGE')
            call strip_off (trim(inputline), PINPT_BERRY%strip_bc_range, 'BERRYC_ERANGE', ' ' , 2) ! get berrycurv_range
            write(message,'(3A)')'   ERANGE:  ',adjustl(trim(PINPT_BERRY%strip_bc_range)),' (ENERGY RANGE TO BE CALCULATED)'  ; write_msg

          case('BERRYC_FNM','BERRYC_FNAME') ! BERRY CURVATURE output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%bc_filenm = trim(desc_str)
            write(message,'(3A)')'   FONAME:  ',trim(PINPT_BERRY%bc_filenm),' (OUTPUT FILE NAME HEADER for BERRY CURVATURE)'  ; write_msg
            PINPT_BERRY%flag_bc_filenm_provided = .true.

          case('BERRYC_METHOD') ! BERRY CURVATURE evalulation method
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            if(trim(desc_str) .eq. 'KUBO') then
              PINPT_BERRY%flag_bc_method_kubo  = .true.
              PINPT_BERRY%flag_bc_method_resta = .false.
              PINPT%flag_get_orbital = .true.
              write(message,'(A)')'   METHOD:  Berry curvature (Omega) is obtained via Kubo formula, i.e., '  ; write_msg
              write(message,'(A)')'             Omega_z(k) = -2Im*SIGMA_{m/=m} v_nm,x * v_mn,y / (e_m - e_n)^2'  ; write_msg
              write(message,'(A)')'                 v_nm,i = <u_nk|dH(k)/dk_i|u_mk> .'  ; write_msg
              write(message,'(A)')'             See Ref. [Wang et. al., PRB 74, 195118 (2006)]'  ; write_msg
            elseif(trim(desc_str) .eq. 'RESTA' .or. trim(desc_str) .eq. 'FUKUI') then
              PINPT_BERRY%flag_bc_method_kubo = .false.
              PINPT_BERRY%flag_bc_method_resta = .true.
              PINPT%flag_get_orbital = .true.
              write(message,'(A)')'   METHOD:  Berry curvature (Omega) is obtained via descritization method,'  ; write_msg
              write(message,'(A)')'            where PI_i runs over closed loop of k normal to kx,ky plane, i.e.,'  ; write_msg
              write(message,'(A)')'             Omega_z(k) = -Im ln PI_i det M(k_i, k_i+1) / ka'  ; write_msg
              write(message,'(A)')'             Here, M is overlap matrix and the matrix element M_mn is as follow'  ; write_msg
              write(message,'(A)')'             M_mn(k_i, k_i+1) = <u_m(k)|u_n(k+1)>, and ka is area of the closed loop'  ; write_msg
              write(message,'(A)')'            See Ref. [Resta, J. Phys.: Condens. Matter 12, R107 (2000)]'  ; write_msg
            endif

          case('BERRYC_NKDIV')
            i_dummy = nitems(inputline) -1
            if(i_dummy .eq. 1) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%bc_nkdiv(1)
              PINPT_BERRY%bc_nkdiv(2)=1
              PINPT_BERRY%bc_nkdiv(3)=1
            elseif(i_dummy .eq. 2) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%bc_nkdiv(1:2)
              PINPT_BERRY%bc_nkdiv(3)=1
            elseif(i_dummy .eq. 3) then
              read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%bc_nkdiv(1:3)
            endif

          case('BERRYC_DIMENSION')  ! which direction to get Chern number integration? 
            call strip_off (trim(inputline), dummy, 'BERRYC_DIMENSION', ' ' , 2)
            i_dummy=index(dummy,'D')
            if(i_dummy .eq. 0) then
              write(message,'(A)')' !WARN!: BERRYC_DIMENSION should be start with one of following strings:'  ; write_msg
              write(message,'(A)')'  3D  => BERRYC_DIMENSION 3D'  ; write_msg
              write(message,'(A)')'  2D  => BERRYC_DIMENSION 2D:z # z is normal direction of surface'  ; write_msg
              write(message,'(A)')'  1D  => BERRYC_DIMENSION 1D:x # x is parallel direction of 1D system'  ; write_msg
              write(message,'(A)')'  For 2D and 1D, "x", "y", or "z" indicates the normal and parallel'  ; write_msg
              write(message,'(A)')'  direction of the system'  ; write_msg
              write(message,'(A)')'  Exit program anyway... Check your "SET BERRYCURV" in "INCAR-TB" again.'  ; write_msg
              stop
            endif

            call str2int(dummy(1:i_dummy-1),PINPT_BERRY%bc_dimension)

            if(PINPT_BERRY%bc_dimension .eq. 3) then
              deallocate(PINPT_BERRY%bc_axis)
              allocate(PINPT_BERRY%bc_axis(3))
              PINPT_BERRY%bc_axis(:) = (/1,2,3/)
              write(message,'(3A)')'   BC_DIM: ', trim(dummy), ' mode (three direction will be checked)'  ; write_msg

            elseif(PINPT_BERRY%bc_dimension .eq. 2) then
              deallocate(PINPT_BERRY%bc_axis)
              allocate(PINPT_BERRY%bc_axis(1))
              i_dummy = index(dummy,':')
              if(i_dummy .eq. 0) then
                PINPT_BERRY%bc_axis(1) = 3 ! default : third reciprocal axis will be chosen for the WCC surface normal axis
                write(message,'(3A)')'   BC_DIM: ', trim(dummy), ' mode (three direction will be checked)'  ; write_msg
                write(message,'( A)')'  BC_AXIS: B3 direction of the reciprocal lattice vector (default)'  ; write_msg
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','b1', 'X', 'A', 'B1')
                    PINPT_BERRY%bc_axis(1) = 1
                    write(message,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode'  ; write_msg
                    write(message,'( A)')'  BC_AXIS: B1 direction of the reciprocal lattice vector (B2-B3 plane, yz)'  ; write_msg
                  case('2','y','b','b2', 'Y', 'B', 'B2')
                    PINPT_BERRY%bc_axis(1) = 2
                    write(message,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode'  ; write_msg
                    write(message,'( A)')'  BC_AXIS: B2 direction of the reciprocal lattice vector (B3-B1 plane, zx)'  ; write_msg
                  case('3','z','c','b3', 'Z', 'C', 'B3')
                    PINPT_BERRY%bc_axis(1) = 3
                    write(message,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode'  ; write_msg
                    write(message,'( A)')'  BC_AXIS: B3 direction of the reciprocal lattice vector (B1-B2 plane, xy)'  ; write_msg
                endselect
              endif

            elseif(PINPT_BERRY%bc_dimension .eq. 1) then
              deallocate(PINPT_BERRY%bc_axis)
              allocate(PINPT_BERRY%bc_axis(1))
              i_dummy = index(dummy, ':')
              if(i_dummy .eq. 0) then
                PINPT_BERRY%bc_axis(1) = 1
                write(message,'(3A)')'   BC_DIM: ', trim(dummy), ' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                write(message,'( A)')'  BC_AXIS: B1 direction of the reciprocal lattice vector (default, parallel to B1)'  ; write_msg
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','kx', 'b1', 'X', 'A', 'B1', 'KX')
                    PINPT_BERRY%bc_axis(1) = 1
                    write(message,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                    write(message,'( A)')'  BC_AXIS: B1 direction of the reciprocal lattice vector (parallel to B1)'  ; write_msg
                  case('2','y','b','ky', 'b2', 'Y', 'B', 'B2', 'KY')
                    PINPT_BERRY%bc_axis(1) = 2
                    write(message,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                    write(message,'( A)')'  BC_AXIS: B2 direction of the reciprocal lattice vector (parallel to B2)'  ; write_msg
                  case('3','z','c','kz', 'b3', 'Z', 'C', 'B3', 'KZ')
                    PINPT_BERRY%bc_axis(1) = 3
                    write(message,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'  ; write_msg
                    write(message,'( A)')'  BC_AXIS: B3 direction of the reciprocal lattice vector (parallel to B3)'  ; write_msg
                endselect

              endif
            endif

          case('SET_PHASE', 'GET_PHASE', 'PHASE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_bc_phase

        end select
      enddo set_berryc

      return
   endsubroutine

   subroutine set_symmetry_check(PINPT, PINPT_BERRY, desc_str)
      type(incar)         :: PINPT
      type(berry)    :: PINPT_BERRY
      character*132          inputline
      character*40           desc_str
      character(*), parameter :: func = 'set_symmetry_check'
      integer*4, external :: nitems
      logical,   external :: flag_number
      integer*4              i_dummy, i_continue
      integer*4              i, nkp
      character*10           kp_name(10)
      real*8                 kp(3,10) !reciprocal unit 
      PINPT%flag_get_symmetry = .true.
      nkp = 0

      ! default rotation matrix for the coordinates
      PINPT_BERRY%symmetry_operator(:,1) = (/-1d0,    0d0,    0d0/)
      PINPT_BERRY%symmetry_operator(:,2) = (/   0d0, -1d0,    0d0/)
      PINPT_BERRY%symmetry_operator(:,3) = (/   0d0,    0d0, -1d0/)


 set_symm:   do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('SYMMETRY_KP', 'SYMM_KP')
            ! KPOINT (usually time-reversal invariant momentum point would be meaningful)
            nkp = nkp + 1

            i_dummy = nitems(inputline)

            if(i_dummy .eq. 4) then
              read(inputline,*,iostat=i_continue) desc_str, kp(1:3, nkp)
              write(kp_name(nkp),'(A,I0)')'KP',nkp
              write(message,'(A,3F10.5)')'  SYMM_KP: ', kp(:,nkp), kp_name(nkp)  ; write_msg
            elseif(i_dummy .gt. 4) then
              read(inputline,*,iostat=i_continue) desc_str, kp(1:3, nkp), kp_name(nkp)
              if(.not. flag_number(kp_name(nkp))) then
                write(message,'(A,3F10.5,2x,A)')'  SYMM_KP: ', kp(:,nkp), trim(kp_name(nkp))  ; write_msg
              else
                write(message,'(4A)')'    !WANR! Check SYMMETRY_EIG   SETTING tags ->',trim(desc_str), ' ', trim(func)  ; write_msg
                stop
              endif
            endif
          
          case('SYMMETRY_ORIGIN', 'ORIGIN', 'ORIGIN_SHIFT')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%symmetry_origin(1:3)
            write(message,'(A,3F10.5,A)')'   ORIGIN: ', PINPT_BERRY%symmetry_origin(1:3), ' (used in SYMMETRY_EIG)'  ; write_msg

          case('SYMM_OP1','SYMMETRY_OP1','ROTATION_MAT1','ROTATION1')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%symmetry_operator(1:3,1)
          case('SYMM_OP2','SYMMETRY_OP2','ROTATION_MAT2','ROTATION2')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%symmetry_operator(1:3,2)
          case('SYMM_OP3','SYMMETRY_OP3','ROTATION_MAT3','ROTATION3')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%symmetry_operator(1:3,3)

          case('NOCC', 'NOCCUPIED','N_OCC','N_OCCUPIED')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%noccupied
          case('ROT_ANGLE', 'ROTATION_ANGLE',  'ANGLE_ROT', 'ANGLE_ROTATION')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%symmetry_theta

          case('PRINT_HAMILTONIAN', 'PRINT_MATRIX', 'PRINT_HAM' )
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%flag_print_hamiltonian_symmetry
          case('SYMMETRY_PHASE' )
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%flag_symmetry_phase

        end select
      enddo set_symm

      if( nkp .eq. 0) then
        PINPT_BERRY%symmetry_nkpoint = 1
        allocate(PINPT_BERRY%symmetry_kpoint(3,1))
        allocate(PINPT_BERRY%symmetry_kpoint_reci(3,1))
        allocate(PINPT_BERRY%symmetry_kpoint_name(1))
        write(message,'(A,3F10.5,2x,A)')'SYMMETRY_KP: ', (/0d0, 0d0, 0d0/),'Gamma'  ; write_msg
        PINPT_BERRY%symmetry_kpoint_reci(:,1) = (/0d0, 0d0, 0d0/) ! set default
        PINPT_BERRY%symmetry_kpoint_name(1) = 'Gamma'

      elseif( nkp .ge. 1) then
        PINPT_BERRY%symmetry_nkpoint = nkp
        allocate(PINPT_BERRY%symmetry_kpoint(3,nkp))
        allocate(PINPT_BERRY%symmetry_kpoint_reci(3,nkp))
        allocate(PINPT_BERRY%symmetry_kpoint_name(nkp))
        PINPT_BERRY%symmetry_kpoint_reci(:,:) = kp(1:3,1:nkp)
        do i = 1, nkp
          PINPT_BERRY%symmetry_kpoint_name(i) = kp_name(i)
        enddo
      endif   
   
      write(message,'(A,3(F8.5,2x),A)')'[SYMMETRY] [  ', PINPT_BERRY%symmetry_operator(:,1),']                          '  ; write_msg
      write(message,'(A,3(F8.5,2x),A)')'[OPERATOR]=[  ', PINPT_BERRY%symmetry_operator(:,2),'] => R(inv) = S * ( R - ORIGIN )'  ; write_msg
      write(message,'(A,3(F8.5,2x),A)')'[    S   ] [  ', PINPT_BERRY%symmetry_operator(:,3),']                          '  ; write_msg

      return
   endsubroutine

   subroutine set_parity_check(PINPT, PINPT_BERRY, desc_str)
      type(incar)         :: PINPT
      type(berry)    :: PINPT_BERRY
      character*132          inputline
      character*40           desc_str
      character(*), parameter :: func = 'set_parity_check'
      integer*4, external :: nitems
      logical,   external :: flag_number
      integer*4              i_dummy, i_continue
      integer*4              i, nkp
      character*10           kp_name(10)
      real*8                 kp(3,10) !reciprocal unit 
      PINPT%flag_get_parity = .true.
      nkp = 0

      ! default rotation matrix for the coordinates
      PINPT_BERRY%parity_operator(:,1) = (/-1d0,    0d0,    0d0/)
      PINPT_BERRY%parity_operator(:,2) = (/   0d0, -1d0,    0d0/)
      PINPT_BERRY%parity_operator(:,3) = (/   0d0,    0d0, -1d0/)


 set_parity: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('PARITY_KP', 'TRIM_POINT', 'TRIM', 'TRIM_KP', 'TRIM_KPOINT', 'PARITY_KPOINT') 
            ! KPOINT (usually time-reversal invariant momentum point would be meaningful)
            nkp = nkp + 1

            i_dummy = nitems(inputline)

            if(i_dummy .eq. 4) then
              read(inputline,*,iostat=i_continue) desc_str, kp(1:3, nkp)
              write(kp_name(nkp),'(A,I0)')'KP',nkp
              write(message,'(A,3F10.5)')'PARITY_KP: ', kp(:,nkp), kp_name(nkp)  ; write_msg
            elseif(i_dummy .gt. 4) then
              read(inputline,*,iostat=i_continue) desc_str, kp(1:3, nkp), kp_name(nkp)
              if(.not. flag_number(kp_name(nkp))) then
                write(message,'(A,3F10.5,2x,A)')'PARITY_KP: ', kp(:,nkp), trim(kp_name(nkp))  ; write_msg
              else
                write(message,'(4A)')'    !WANR! Check PARITY_CHECK SETTING tags ->',trim(desc_str), ' ', trim(func)  ; write_msg
                stop
              endif
            endif
          
          case('PARITY_ORIGIN', 'ORIGIN', 'ORIGIN_SHIFT')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_origin(1:3)
            write(message,'(A,3F10.5,A)')'   ORIGIN: ', PINPT_BERRY%parity_origin(1:3), ' (used in PARITY_CHECK)'  ; write_msg

          case('PARITY_OP1','SYMMETRY_OP1','ROTATION_MAT1','ROTATION1')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_operator(1:3,1)
          case('PARITY_OP2','SYMMETRY_OP2','ROTATION_MAT2','ROTATION2')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_operator(1:3,2)
          case('PARITY_OP3','SYMMETRY_OP3','ROTATION_MAT3','ROTATION3')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_operator(1:3,3)

          case('NOCC', 'NOCCUPIED','N_OCC','N_OCCUPIED')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%noccupied

          case('PRINT_HAMILTONIAN', 'PRINT_MATRIX', 'PRINT_HAM' )
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%flag_print_hamiltonian_parity
          case('PARITY_PHASE' )
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%flag_parity_phase

        end select
      enddo set_parity

      if( nkp .eq. 0) then
        PINPT_BERRY%parity_nkpoint = 1
        allocate(PINPT_BERRY%parity_kpoint(3,1))
        allocate(PINPT_BERRY%parity_kpoint_reci(3,1))
        allocate(PINPT_BERRY%parity_kpoint_name(1))
        write(message,'(A,3F10.5,2x,A)')'PARITY_KP: ', (/0d0, 0d0, 0d0/),'Gamma'  ; write_msg
        PINPT_BERRY%parity_kpoint_reci(:,1) = (/0d0, 0d0, 0d0/) ! set default
        PINPT_BERRY%parity_kpoint_name(1) = 'Gamma'

      elseif( nkp .ge. 1) then
        PINPT_BERRY%parity_nkpoint = nkp
        allocate(PINPT_BERRY%parity_kpoint(3,nkp))
        allocate(PINPT_BERRY%parity_kpoint_reci(3,nkp))
        allocate(PINPT_BERRY%parity_kpoint_name(nkp))
        PINPT_BERRY%parity_kpoint_reci(:,:) = kp(1:3,1:nkp)
        do i = 1, nkp
          PINPT_BERRY%parity_kpoint_name(i) = kp_name(i)
        enddo
      endif   
   
      write(message,'(A,3(F8.5,2x),A)')'  PARITY  [  ', PINPT_BERRY%parity_operator(:,1),']                          '  ; write_msg
      write(message,'(A,3(F8.5,2x),A)')' OPERATOR=[  ', PINPT_BERRY%parity_operator(:,2),'] => R(inv) = P * ( R - ORIGIN )'  ; write_msg
      write(message,'(A,3(F8.5,2x),A)')'     P    [  ', PINPT_BERRY%parity_operator(:,3),']                          '  ; write_msg

      return
   endsubroutine

   subroutine set_magnetism_tag(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_magnetism_tag'
      character*40  str2lowcase
      external      str2lowcase

      read(inputline,*,iostat=i_continue) desc_str,desc_str
      desc_str = str2lowcase(desc_str)
      if(desc_str(1:6) .eq. 'noncol' .or. desc_str(1:7) .eq. 'non-col' .or. desc_str(1:2) .eq. 'nc') then
        PINPT%flag_noncollinear = .true.
        PINPT%flag_collinear    = .false.
        write(message,'(A)')'  LNONCOL: .TRUE. (non-collinear = .TRUE.)'  ; write_msg
        write(message,'(A)')'    ISPIN: 2 (noncollinear = .TRUE.)'  ; write_msg
        write(message,'(A)')'  ISPINOR: 2 (noncollinear = .TRUE.)'  ; write_msg
        PINPT%ispin = 2
        PINPT%ispinor = 2
        PINPT%nspin = 1
      elseif(desc_str(1:6) .eq. 'collin' .or. desc_str(1:3) .eq. 'col' .or. desc_str(1:2) .eq. 'cl') then
        PINPT%flag_collinear    = .true.
        PINPT%flag_noncollinear = .false.
        write(message,'(A)')'    ISPIN: 2 (collinear = .TRUE.)'  ; write_msg
        write(message,'(A)')'  ISPINOR: 1 (collinear = .TRUE.)'  ; write_msg
        PINPT%ispin = 2
        PINPT%ispinor = 1
        PINPT%nspin =2
      elseif(desc_str(1:6) .eq. 'nonmag' .or. desc_str(1:2) .eq. 'nm' .or. desc_str(1:7) .eq. 'non-mag') then
        PINPT%flag_collinear    = .false.
        PINPT%flag_noncollinear = .false.
        write(message,'(A)')'    ISPIN: 1 (non-magnetic = .true.)'  ; write_msg
        write(message,'(A)')'  ISPINOR: 1 (non-magnetic = .true.)'  ; write_msg
        PINPT%ispin = 1
        PINPT%ispinor = 1
        PINPT%nspin =1
      else
        write(message,'(A)')'  !WARNING! TYPMAG is not properly set. '  ; write_msg
        write(message,'(A)')'  !WARNING! usage: for non-magnetic run : "NM", "NONMAGNETIC", "NON-MAGNETIC" '  ; write_msg
        write(message,'(A)')'  !WARNING!        for collinear    run : "COLLINEAR", "COL", "CL" '  ; write_msg
        write(message,'(A)')'  !WARNING!        for noncollinear run : "NONCOLLINEAR", "NON-COLLINEAR", "NC"'  ; write_msg
        write(message,'(A)')'  !WARNING!        are acceptable (case insensitive). Other option is not allowed. Exit program.'  ; write_msg
        stop
      endif

      return
   endsubroutine

   subroutine set_spin_orbit_tag(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_spin_orbit_tag'
      character*40  str2lowcase
      external      str2lowcase

      read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_soc
      if(PINPT%flag_soc .and. .not. PINPT%flag_collinear) then
        write(message,'(A)')'    LSORB: .TRUE. (non-collinear = .TRUE. .AND. spin-orbit = .TRUE.)'  ; write_msg
        PINPT%ispin = 2
        PINPT%ispinor = 2
        PINPT%flag_noncollinear = .true.
      elseif(.not.PINPT%flag_soc)then
        write(message,'(A)')'    LSORB: .FALSE. ( spin-orbit = .FALSE.)'  ; write_msg
      elseif(PINPT%flag_soc .and. PINPT%flag_collinear) then
        write(message,'(A)')'  !WARNING! LSORB and TYPMAG is not properly set. '  ; write_msg
        write(message,'(A)')'  !WARNING! COLLINEAR = .true. and LSORB = .true. which is not proper.'  ; write_msg
        write(message,'(A)')'  !WARNING! Both option are not available simultaneously. Exit program.'  ; write_msg
        stop
      endif

      return
   endsubroutine

   subroutine set_tbparam_file(PINPT, PWGHT, param_const, param_const_nrl, inputline)
      type(incar)  ::  PINPT
      type(weight) ::  PWGHT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str, dummy
      real*8        param_const(5,max_nparam), param_const_nrl(5,4,max_nparam)
      integer*4     nitems, i_dummy
      external      nitems
      character(*), parameter :: func = 'set_tbparam_file'

      i_dummy = nitems(inputline) - 1

      if(.not. PINPT%flag_pfile) then
        if(i_dummy .eq. 1) then
          read(inputline,*,iostat=i_continue) desc_str, PINPT%pfilenm
        elseif(i_dummy .gt. 1) then
          read(inputline,*,iostat=i_continue) desc_str, PINPT%pfilenm, dummy
          if(trim(dummy) .eq. 'USE_WEIGHT' .or. trim(dummy) .eq. 'use_weight') then
            PINPT%flag_use_weight = .true.
          endif
        endif

        if(allocated(PINPT%param) .or. allocated(PINPT%param_name) .or. PINPT%flag_pincar) then
           PINPT%flag_pincar=.false.
           deallocate(PINPT%param)
           deallocate(PINPT%param_name)
           write(message,'(A)')'  !WARN!  TB-parameter is already set by TBPARAM in the INCAR-TB,'  ; write_msg
           write(message,'(A)')'  !WARN!  however, since the external TB-parameter file is provided,'  ; write_msg
           write(message,'(A,A,A)')'  !WARN!  those valeus from ',trim(PINPT%pfilenm),' will be read in priori'  ; write_msg
        endif
      endif
      write(message,'(A,A)')' PARA_FNM:  ',trim(PINPT%pfilenm)  ; write_msg
      call read_param(PINPT, PWGHT, param_const, param_const_nrl)

      return
   endsubroutine

   subroutine set_kpoint_file(PINPT, flag_kfile_ribbon, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_kpoint_file'
      logical       flag_kfile_ribbon 
      integer*4     nitems, i_dummy
      external      nitems
      
      i_dummy = nitems(inputline) -1

      if(i_dummy .eq. 1) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%kfilenm
        if(flag_kfile_ribbon) then
          PINPT%kfilenm = PINPT%ribbon_kfilenm
          write(message,'(A,A)')' KPTS_FNM:  (for ribbon) ',trim(PINPT%kfilenm)  ; write_msg
        else
          read(inputline,*,iostat=i_continue) desc_str, PINPT%kfilenm
          write(message,'(A,A)')' KPTS_FNM:  ',trim(PINPT%kfilenm)  ; write_msg
        endif
      elseif(i_dummy .eq. 2) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%kfilenm, PINPT%kline_type  
        if(flag_kfile_ribbon) then
          PINPT%kfilenm = PINPT%ribbon_kfilenm
          write(message,'(A,A)')' KPTS_FNM:  (for ribbon) ',trim(PINPT%kfilenm)  ; write_msg
        else
          read(inputline,*,iostat=i_continue) desc_str, PINPT%kfilenm
          write(message,'(A,A)')' KPTS_FNM:  ',trim(PINPT%kfilenm)  ; write_msg
        endif
      endif

      return
   endsubroutine

   subroutine set_geom_file(PINPT, inputline, imode)
      type(incar)  ::  PINPT
      integer*4     i_continue
      integer*4     nitems
      integer*4     imode
      character*132 inputline
      character*40  desc_str
      external      nitems
      integer*4     i_dummy
      logical       flag, flag_logical
      character(*), parameter :: func = 'set_geom_file'

      if(imode .eq. 1) then
        i_dummy = nitems(inputline) -1
        if(i_dummy .eq. 1) then
          read(inputline,*,iostat=i_continue) desc_str, PINPT%gfilenm
          write(message,'(A,A)')' GEOM_FNM:  ',trim(PINPT%gfilenm)  ; write_msg

        elseif(i_dummy .ge. 2) then
          read(inputline,*,iostat=i_continue) desc_str, PINPT%gfilenm, desc_str
          call str2logical(trim(desc_str),flag_logical,flag)
          if(flag_logical) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%gfilenm, PINPT%flag_report_geom
            write(message,'(A,A)')' GEOM_FNM:  ',trim(PINPT%gfilenm)  ; write_msg
            write(message,'(A  )')'         :  PRINT_GEOM = .TRUE.'  ; write_msg

          elseif(.not. flag_logical) then
            if(trim(desc_str) .eq. 'PRINT_GEOM') then
              read(inputline,*,iostat=i_continue) desc_str, PINPT%gfilenm, desc_str, desc_str
              call str2logical(trim(desc_str),flag_logical,flag)
              if(flag_logical) then 
                read(inputline,*,iostat=i_continue) desc_str, PINPT%gfilenm, desc_str, PINPT%flag_report_geom
                write(message,'(A,A)')' GEOM_FNM:  ',trim(PINPT%gfilenm)  ; write_msg
                write(message,'(A  )')'         :  PRINT_GEOM = .TRUE.'  ; write_msg
              endif
            endif
          endif
        endif

      elseif(imode .eq. 2) then
        read(inputline, *,iostat=i_continue) desc_str, desc_str
        call str2logical(trim(desc_str),flag_logical,flag)
        if(flag_logical) then
          read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_report_geom
          if(PINPT%flag_report_geom) then
            write(message,'(A  )')'         :  PRINT_GEOM = .TRUE.'  ; write_msg
          elseif(.not.PINPT%flag_report_geom) then
            write(message,'(A  )')'         :  PRINT_GEOM = .FALSE.'  ; write_msg
          endif
        endif
      endif

      return
   endsubroutine

   subroutine set_tbparam_out_file(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_tbparam_out_file'

      read(inputline,*,iostat=i_continue) desc_str, PINPT%pfileoutnm
      write(message,'(A,A40)')' POUT_FNM:  ',PINPT%pfileoutnm  ; write_msg
      PINPT%flag_print_param = .true.

      return
   endsubroutine

   subroutine set_hopping_type(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      integer*4     nitems, i_dummy 
      character(*), parameter :: func = 'set_hopping_type'
      external      nitems

     !i_dummy = nitems(inputline) - 1

     !if(i_dummy .eq. 1) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_slater_koster
        if(PINPT%flag_slater_koster) then
          !write(message,'(A,A)')' TYPE_HOP:  ','SLATER_KOSTER: mode 1'  ; write_msg
          write(message,'(A,A)')' TYPE_HOP:  ','SLATER_KOSTER'  ; write_msg
        elseif(.not. PINPT%flag_slater_koster) then
          write(message,'(A,A)')' TYPE_HOP:  ','USER_DEFINED'  ; write_msg
        endif

    ! elseif(i_dummy .eq. 3) then
 
    !   read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_slater_koster, desc_str, PINPT%slater_koster_type
    !   if(PINPT%flag_slater_koster) then
    !     write(message,'(A,A,I0)')' TYPE_HOP:  ','SLATER_KOSTER: mode ', PINPT%slater_koster_type   ; write_msg
    !     ! parameterization scheme : if 11 : see Mehl & Papaconstantopoulos PRB 54, 4519 (1996)
    !   elseif(.not. PINPT%flag_slater_koster) then
    !     write(message,'(A,A)')' TYPE_HOP:  ','USER_DEFINED'  ; write_msg
    !   endif

    ! endif

      return
   endsubroutine

   subroutine set_local_charge(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_local_charge'

      read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_local_charge
      if(PINPT%flag_local_charge) then
        write(message,'(A)')'   LOCCHG: .TRUE. (read local charge density from GFILE)'  ; write_msg
      elseif(.not. PINPT%flag_local_charge) then
        write(message,'(A)')'   LOCCHG: .FALSE.'  ; write_msg
      endif

      return
   endsubroutine

   subroutine set_plus_U_scheme(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_plus_U_scheme'

      read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_plus_U
      if(PINPT%flag_plus_U) then
        write(message,'(A)')'   PLUS+U: .TRUE. '  ; write_msg
      elseif(.not. PINPT%flag_plus_U) then
        write(message,'(A)')'   LOCCHG: .FALSE.'  ; write_msg
      endif

      return
   endsubroutine

   subroutine set_load_nntable(PINPT, inputline, tag_name)
      type(incar)  :: PINPT
      integer*4       i_continue
      character*132   inputline
      character*40    desc_str, tag_name
      character(*), parameter :: func = 'set_load_nntable'
      logical         flag_load_nntable, flag_exist
      integer*4       nitems, i_dummy
      external        nitems

      read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_load_nntable

      if(PINPT%flag_load_nntable) then
        i_dummy = nitems(inputline) -1

        if(i_dummy .eq. 2) then
          read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_load_nntable, PINPT%nnfilenm
          inquire(file=PINPT%nnfilenm,exist=flag_exist)
          if(.not.flag_exist) then
            write(message,'(5A)')'  !WARN! You requested a HOPPING file via "', trim(tag_name),'" tag, but the file "', trim(PINPT%nnfilenm),'" does not exist! Exit...'  ; write_msg
            stop
          endif

        elseif(i_dummy .eq. 1) then
          write(message,'(4A)')'  !WARN! You requested a HOPPING file via "', trim(tag_name), '" tag, but the file name has not been defined. Exit...'  ; write_msg
          stop
        endif
      endif

      return
   endsubroutine

   subroutine set_band_order(PINPT, inputline)
      type(incar)  :: PINPT
      integer*4       i_continue
      character*132   inputline, inputline_dummy
      character*40    desc_str_
      character*40    desc_str
      character*40    str2lowcase
      integer*4       nitems, i_dummy
      external        nitems, str2lowcase
      real*8          dummy
      character(*), parameter :: func = 'set_band_order'

      call comment_off(trim(inputline), inputline_dummy) ! comment off (take out comment mark)
      inputline = trim(inputline_dummy)

      i_dummy = nitems(inputline) -1

      read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_get_band_order
  
      if(PINPT%flag_get_band_order) then
        if(i_dummy .eq. 3) then
          read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_get_band_order, desc_str_ , dummy
          desc_str = str2lowcase(desc_str_)
          if(trim(desc_str) .eq. 'cutoff'  .or. trim(desc_str) .eq. 'ov_cut') then
            PINPT%band_order_overlap_cutoff = dummy
            write(message,'(A,L,A,F9.4)')'   LORDER: ',PINPT%flag_get_band_order, '    OV_CUT: ',PINPT%band_order_overlap_cutoff  ; write_msg
          else
            write(message,'( A)')' '  ; write_msg
            write(message,'(3A)')'   !!WARN!!  You specified ',trim(desc_str), ' which is not recognized tag.'  ; write_msg
            write(message,'( A)')'             TBFIT will set OV_CUT to default value: sqrt(2)/2. check: LORDER tag'  ; write_msg
            write(message,'(A,L,A,F9.4)')'   LORDER: ',PINPT%flag_get_band_order, ' OV_CUT: ',sqrt(2d0)/2d0  ; write_msg
            PINPT%band_order_overlap_cutoff = sqrt(2d0)/2d0
          endif
        elseif( (i_dummy .eq. 4 .or. i_dummy .eq. 2) .and. index(inputline,'PRINT_ONLY') .ne. 0 ) then
          PINPT%flag_get_band_order_print_only = .TRUE.
          select case(i_dummy)
            case(4)
              if(index(inputline,'PRINT_ONLY') .gt. index(inputline,'OV_CUT')) then
                read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_get_band_order, desc_str_ , dummy
              elseif(index(inputline,'PRINT_ONLY') .lt. index(inputline,'OV_CUT')) then
                read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_get_band_order, desc_str_, desc_str_ , dummy
              endif
              desc_str = str2lowcase(desc_str_)
              if(trim(desc_str) .eq. 'cutoff'  .or. trim(desc_str) .eq. 'ov_cut') then
                PINPT%band_order_overlap_cutoff = dummy
                write(message,'(A,L,A,F9.4)')'   LORDER: ',PINPT%flag_get_band_order, '    OV_CUT: ',PINPT%band_order_overlap_cutoff  ; write_msg
                write(message,'(A)'         )'   LORDER + PRINT_ONLY requested -> re-order will not work in fitting routines (if TBFIT = .TRUE.)'  ; write_msg
                write(message,'(A)'         )'                       but thre re-ordered band will be written at the end of band calculations.'  ; write_msg
              else
                write(message,'( A)')' '  ; write_msg
                write(message,'(3A)')'   !!WARN!!  You specified ',trim(desc_str_), ' which is not recognized tag.'  ; write_msg
                write(message,'( A)')'             TBFIT will set OV_CUT to default value: sqrt(2)/2. check: LORDER tag'  ; write_msg
                write(message,'(A,L,A,F9.4)')'   LORDER: ',PINPT%flag_get_band_order, ' OV_CUT: ',sqrt(2d0)/2d0  ; write_msg
                PINPT%band_order_overlap_cutoff = sqrt(2d0)/2d0
                write(message,'(A)'         )'   LORDER + PRINT_ONLY requested -> re-order will not work in fitting routines (if TBFIT = .TRUE.)'  ; write_msg
                write(message,'(A)'         )'                       but thre re-ordered band will be written at the end of band calculations.'  ; write_msg
              endif
            case(2)
              read(inputline, *,iostat=i_continue) desc_str, PINPT%flag_get_band_order
              write(message,'(A,L,A,F9.4)')'   LORDER: ',PINPT%flag_get_band_order, '    OV_CUT: ',sqrt(2d0)/2d0  ; write_msg
              PINPT%band_order_overlap_cutoff = sqrt(2d0)/2d0
              write(message,'(A)'         )'   LORDER + PRINT_ONLY requested -> re-order will not work in fitting routines (if TBFIT = .TRUE.)'  ; write_msg
              write(message,'(A)'         )'                       but thre re-ordered band will be written at the end of band calculations.'  ; write_msg
          end select
        elseif(i_dummy .eq. 1) then
          write(message,'(A,L,A,F9.4)')'   LORDER: ',PINPT%flag_get_band_order, '    OV_CUT: ',sqrt(2d0)/2d0  ; write_msg
          PINPT%band_order_overlap_cutoff = sqrt(2d0)/2d0
        endif
      else
        write(message,'(A,L,A,F9.4)')'   LORDER: ',PINPT%flag_get_band_order  ; write_msg
      endif

      return
   endsubroutine
   subroutine set_target_file(PINPT, flag_read_energy, inputline, input_tag)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_target_file'
      character*40  input_tag
      logical       flag_read_energy
      integer*4     nitems, i_dummy
      external      nitems

      select case(input_tag)
        case('EFILE')
          i_dummy = nitems(inputline) -1
          if(i_dummy .eq. 1) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmu
          elseif(i_dummy .eq. 2) then
            read(inputline,*,iostat=i_continue) desc_str, desc_str
            if(trim(desc_str) .eq. 'VASP' .or. trim(desc_str) .eq. 'vasp') then
              PINPT%efile_type='vasp'
              read(inputline,*,iostat=i_continue) desc_str, desc_str, PINPT%efilenmu !, PINPT%read_energy_column_index
            elseif(trim(desc_str) .eq. 'USER' .or. trim(desc_str) .eq. 'user') then
              PINPT%efile_type='user'
              read(inputline,*,iostat=i_continue) desc_str, desc_str, PINPT%efilenmu !, PINPT%read_energy_column_index
            else
              PINPT%efile_type='user'
              read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmu, PINPT%read_energy_column_index
            endif
          elseif(i_dummy .eq. 3) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efile_type, desc_str, PINPT%itarget_e_start
            if(trim(PINPT%efile_type) .eq. 'VASP' .or. trim(PINPT%efile_type) .eq. 'vasp') then
              PINPT%efile_type='vasp'
              read(inputline,*,iostat=i_continue) desc_str, desc_str, PINPT%efilenmu, PINPT%itarget_e_start !, PINPT%read_energy_column_index
            elseif(trim(PINPT%efile_type) .eq. 'USER' .or. trim(PINPT%efile_type) .eq. 'user' ) then
              PINPT%efile_type='user'
              read(inputline,*,iostat=i_continue) desc_str, desc_str, PINPT%efilenmu, PINPT%itarget_e_start !, PINPT%read_energy_column_index
            endif

          endif
          write(message,'(A,A)')' EDFT_FNM: ',trim(PINPT%efilenmu)  ; write_msg
          flag_read_energy=.true.

        case('EFILEU')
          i_dummy = nitems(inputline) -1
          if(i_dummy .eq. 1) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmu
          elseif(i_dummy .eq. 2) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmu, PINPT%read_energy_column_index
          endif
          write(message,'(A,A)')' EDFT_FNM: spin-up= ',trim(PINPT%efilenmu)  ; write_msg
          flag_read_energy=.true.

        case('EFILED')
          i_dummy = nitems(inputline) -1
          if(i_dummy .eq. 1) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmd
          elseif(i_dummy .eq. 2) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmd, PINPT%read_energy_column_index_dn
          endif
          write(message,'(A,A)')' EDFT_FNM: spin-dn= ',trim(PINPT%efilenmd)  ; write_msg
          flag_read_energy=.true.

      end select

      return
   endsubroutine

   subroutine set_target_band_init_fina(PWGHT, inputline, input_tag)
      type(weight)  ::  PWGHT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_taget_band_init_fina'
      character*40  input_tag
      logical       flag_read_energy
      integer*4     nitems, i_dummy
      external      nitems

      select case(input_tag)
        case('IBAND')
          read(inputline,*,iostat=i_continue) desc_str, PWGHT%iband
          write(message,'(A,I8)')' INI_BAND:',PWGHT%iband  ; write_msg

         case('FBAND')
          read(inputline,*,iostat=i_continue) desc_str, PWGHT%fband
          write(message,'(A,I8)')' FIN_BAND:',PWGHT%fband  ; write_msg

      end select

      return
   endsubroutine

   subroutine set_target_scissor(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_target_scissor'

      read(inputline,*,iostat=i_continue) desc_str, PINPT%i_scissor, PINPT%r_scissor
      PINPT%flag_scissor = .true.
      write(message,'(A,F8.2,A,I5,A)')'  SCISSOR: .TRUE. => EDFT(n,k) + ',PINPT%r_scissor,' (if n >=',PINPT%i_scissor,')'  ; write_msg

      return
   endsubroutine

   subroutine set_ldos_project_print(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4        i_continue, mpierr
      integer*4        nitems
      integer*4        i, k, ii
      character*132    inputline, dummy_
      character*2      str2lowcase
      character*40     desc_str, dummy, dummy1,dummy3
      integer*4        i_dummy,i_dummy1,i_dummy2,i_dummy3,i_dummy4,i_dummy5
      character*40, allocatable :: strip_dummy(:)
      integer*4        i_dummyr(max_dummy)
      character(*), parameter :: func = 'set_ldos_project_print'
      external         nitems, str2lowcase

      PINPT%flag_print_proj_sum = .false.

      call strip_off (inputline, dummy_, ' ', '#', 0) ! cut off unnecessary comments

      if(index(inputline, 'LDOS_SUM') .ge. 1 .or. index(inputline, 'PROJ_SUM') .ge. 1 .or. index(inputline, 'PROJ') .ge. 1) then
        PINPT%flag_print_proj_sum = .true.
      endif

      i_dummy = nitems(inputline) - 1
      if(i_dummy .eq. 1) then
        if(PINPT%flag_print_proj_sum) then
          write(message,'(A)')'  !!!!WARN: PROJ_SUM tag is activated but atom indext to be summed up '  ; write_msg
          write(message,'(A)')'            did not specified.'  ; write_msg
          write(message,'(A)')'            Correct usage is, for example, if you want to resolve'  ; write_msg
          write(message,'(A)')'            contribution of atom 1 to 10 and 15, then you can specify as follows'  ; write_msg
          write(message,'(A)')'            PROJ_SUM .TRUE.  1:10 15'  ; write_msg
          write(message,'(A)')'  Stop program...'  ; write_msg
          kill_job
        endif
        read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_proj
        if(PINPT%flag_print_proj) then
          write(message,'(A)')'   L_LDOS: .TRUE. | print out projected orbital weight for each atom'  ; write_msg
          PINPT%flag_get_orbital = .true.
        elseif( .not. PINPT%flag_print_proj) then
          write(message,'(A)')'   L_LDOS: .FALSE.'  ; write_msg
          PINPT%flag_get_orbital = .false.
        endif

      elseif(i_dummy .gt. 1) then
        read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_print_proj
        read(inputline,*,iostat=i_continue) desc_str,dummy
        if(PINPT%flag_print_proj) then
          PINPT%nproj_sum = PINPT%nproj_sum + 1
          call strip_off (trim(inputline), dummy1, trim(dummy), ' ' , 2)   ! get dos_ensurf
          i_dummy1=index(dummy1,':')
          if(i_dummy1 .eq. 0) then
            if(PINPT%nproj_sum .eq. 1) allocate( PINPT%proj_natom(max_dummy2) )
            if(PINPT%nproj_sum .eq. 1) allocate( PINPT%proj_atom(max_dummy2*10,max_dummy2) )
            PINPT%proj_natom(PINPT%nproj_sum) = nitems(dummy1)
            read(dummy1,*,iostat=i_continue) PINPT%proj_atom(1:PINPT%proj_natom(PINPT%nproj_sum), PINPT%nproj_sum)
            write(message,'(A,I0,2A)')' LDOS_SUM: .TRUE. , Atom_index SET ',PINPT%nproj_sum,' = ',trim(dummy1)  ; write_msg
          elseif(i_dummy1 .ge. 1)then
            i_dummy2 = nitems(dummy1)
            allocate( strip_dummy(i_dummy2) )
            read(dummy1,*,iostat=i_continue) (strip_dummy(i),i=1,i_dummy2)
            ii = 0   
            do i = 1, i_dummy2 
              i_dummy3 = index(strip_dummy(i),':')
              if(i_dummy3 .eq. 0) then
                ii = ii + 1 
                call str2int(strip_dummy(i),i_dummy4)
                i_dummyr(ii) = i_dummy4 
              elseif(i_dummy3 .gt. 1) then
                ii = ii + 1 
                call strip_off(trim(strip_dummy(i)), dummy3, ' ', ':', 0)
                call str2int(dummy3,i_dummy4)
                call strip_off(trim(strip_dummy(i)), dummy3, ':', ' ', 2)
                call str2int(dummy3,i_dummy5)
                i_dummyr(ii:ii+i_dummy5 - i_dummy4) = (/ (k, k=i_dummy4, i_dummy5) /)
                ii = ii + i_dummy5 - i_dummy4 
              endif    
            enddo    
            if(PINPT%nproj_sum .eq. 1) allocate( PINPT%proj_natom(max_dummy2) )
            if(PINPT%nproj_sum .eq. 1) allocate( PINPT%proj_atom(max_dummy2*10,max_dummy2) )
            PINPT%proj_natom(PINPT%nproj_sum) = ii
            deallocate( strip_dummy )
            PINPT%proj_atom(1:ii,PINPT%nproj_sum) = i_dummyr(1:ii)
            write(message,'(A,I0,2A)')' LDOS_SUM: .TRUE. , Atom_index SET ',PINPT%nproj_sum,' = ',trim(dummy1)  ; write_msg
          endif    
        endif

      endif

      return
   endsubroutine

   subroutine set_local_orbital_print(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      integer*4     nitems
      integer*4     i_dummy, mpierr
      character*132 inputline
      character*2   dummy
      character*40  desc_str
      character*4   str2lowcase
      character*6   c_dummy
      character(*), parameter :: func = 'set_local_orbital_print'
      external      nitems, str2lowcase

      PINPT%axis_print_mag = 'rh' ! write <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j by default
      PINPT%flag_print_single = .false.

      i_dummy = nitems(inputline) - 1
      if(i_dummy .eq. 1) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_orbital
        if(PINPT%flag_print_orbital) then
          PINPT%flag_get_orbital = .true.
          PINPT%axis_print_mag = 'rh' ! <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j
          write(message,'(2A)')'  L_ORBIT: .TRUE. | print out projected orbital weight: ', PINPT%axis_print_mag  ; write_msg
        elseif( .not. PINPT%flag_print_orbital) then
          write(message,'(A)')'  L_ORBIT: .FALSE.'  ; write_msg
          PINPT%flag_get_orbital = .false.
          PINPT%axis_print_mag = 'no' 
        endif
    
      elseif(i_dummy .eq. 2) then
        PINPT%flag_print_mag = .TRUE.
        read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_orbital, desc_str !PINPT%axis_print_mag
        desc_str = str2lowcase(trim(desc_str))
        if(desc_str(1:4) .eq. 'bin4') then
          PINPT%axis_print_mag = 'bi'
          PINPT%flag_print_single = .true.
        else
          PINPT%flag_print_single = .false.
          PINPT%axis_print_mag = desc_str(1:2)
        endif

        if(PINPT%flag_print_orbital) then
          if(PINPT%axis_print_mag .eq. 're') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out real part of wavefnc.: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'im') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out imag part of wavefnc.: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'wf') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out total wavefnc.: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'bi') then ! write binary format
            PINPT%flag_get_orbital = .true.
            PINPT%flag_write_unformatted = .true.
            PINPT%axis_print_mag = 'rh' ! <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out projected orbital weight with binary format (unformatted): ', PINPT%axis_print_mag  ; write_msg
          elseif(PINPT%axis_print_mag(1:1) .eq. 'm') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out magnetization <sigma>: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'rh') then 
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out projected orbital weight: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          endif
        elseif( .not. PINPT%flag_print_orbital) then
          write(message,'(A)')'  L_ORBIT: .FALSE.'  ; write_msg
          PINPT%flag_get_orbital = .false.
          PINPT%axis_print_mag = 'no'
        endif
   
      elseif(i_dummy .eq. 3) then
        PINPT%flag_print_mag = .TRUE.
        read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_orbital, PINPT%axis_print_mag, desc_str !c_dummy
        PINPT%axis_print_mag=str2lowcase(PINPT%axis_print_mag)
        desc_str=str2lowcase(trim(desc_str)) 

        if(desc_str(1:4) .eq. 'bin4') then
          c_dummy = 'bi'
          PINPT%flag_print_single = .true.
        else
          PINPT%flag_print_single = .false.
          c_dummy = desc_str(1:2)
        endif

        if(PINPT%flag_print_orbital) then
          if(PINPT%axis_print_mag .eq. 're') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out real part of wavefnc.: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'im') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out imag part of wavefnc.: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'wf') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out total wavefnc.: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'bi') then ! write binary format
            PINPT%flag_get_orbital = .true.
            PINPT%flag_write_unformatted = .true.
            PINPT%axis_print_mag = 'rh' ! <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out projected orbital weight with binary format (unformatted): ',PINPT%axis_print_mag  ; write_msg
            if(PINPT%flag_print_single) then
              write(message,'( A)')'                  | with "single" precision as requiested by "bin4" tag'  ; write_msg
            endif
          elseif(PINPT%axis_print_mag(1:1) .eq. 'm') then
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out magnetization <sigma>: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'rh') then 
            write(message,'(2A)')'  L_ORBIT: .TRUE. | print out projected orbital weight: ', PINPT%axis_print_mag  ; write_msg
            PINPT%flag_get_orbital = .true.
          endif

          if(c_dummy(1:2) .eq. 'bi') then
            PINPT%flag_get_orbital = .true.
            PINPT%flag_write_unformatted = .true.
            write(message,'(1A)')'                  | with binary format (unformatted)'  ; write_msg
!           write(message,'( A)')'                    with "single" precision as requiested by "bin4" tag'  ; write_msg
          endif

        elseif( .not. PINPT%flag_print_orbital) then
          write(message,'(A)')'  L_ORBIT: .FALSE.'  ; write_msg
          PINPT%flag_get_orbital = .false.
          PINPT%axis_print_mag = 'no'
        endif

      endif

      return
   endsubroutine

   subroutine set_kpoint_unit(PKPTS, inputline)
      type(kpoints)  ::  PKPTS
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_kpoint_unit'
      character*40  dummy

      read(inputline,*,iostat=i_continue) desc_str,desc_str
      dummy = trim(desc_str)
      PKPTS%kunit=dummy(1:1)
      if(PKPTS%kunit .eq. 'r' .or. PKPTS%kunit .eq. 'R') then
        PKPTS%kunit = 'R'
        write(message,'(A)')' KPT_UNIT: KPOINT UNIT = Reciprocal'  ; write_msg
      elseif(PKPTS%kunit .eq. 'a' .or. PKPTS%kunit .eq. 'A') then
        PKPTS%kunit = 'A'
        write(message,'(A)')' KPT_UNIT: KPOINT UNIT = Angstrom unit (1/A)'  ; write_msg
      endif

      return
   endsubroutine

   subroutine set_tbfit(PINPT, inputline, input_tag)
      type(incar  )  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character*40  input_tag
      character(*), parameter :: func = 'set_tbfit'

      select case(input_tag)

        case('TBFIT') !set TBFIT or not
          read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_tbfit
          if(PINPT%flag_tbfit) then
             write(message,'(A)')'  L_TBFIT: .TRUE.'  ; write_msg
          elseif(.not. PINPT%flag_tbfit) then
             write(message,'(A)')'  L_TBFIT: .FALSE.'  ; write_msg
          else
             write(message,'(A)')'  L_TBFIT:  unknown input variable.. exit program'  ; write_msg
             stop
          endif

        case('LSTYPE') !set non-linear regression scheme
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ls_type
          if(PINPT%ls_type .eq. 'LMDIF' .or. PINPT%ls_type .eq. 'lmdif' ) then
            PINPT%ls_type = 'LMDIF'
            write(message,'(A,A8,A)')'  LS_TYPE:  ',PINPT%ls_type,', Levenberg-Marquardt with finite-difference for Jacobian'  ; write_msg
          elseif(PINPT%ls_type .eq. 'PIKAIA' .or. PINPT%ls_type .eq. 'GA' .or. &
                 PINPT%ls_type .eq. 'pikaia' .or. PINPT%ls_type .eq. 'ga' ) then
            PINPT%ls_type = 'GA'
            write(message,'(A,A8,A)')'  LS_TYPE:  ',PINPT%ls_type,', Genetic algorithm based on PIKAIA library'  ; write_msg
          elseif(PINPT%ls_type .eq. 'ga+lmdif' .or. PINPT%ls_type .eq. 'GA+LMDIF' .or. &
                 PINPT%ls_type .eq. 'lmdif+ga' .or. PINPT%ls_type .eq. 'LMDIF+GA') then
            PINPT%flag_ga_with_lmdif=.true.
            write(message,'(A,A8,A)')'  LS_TYPE:  ',PINPT%ls_type,', Genetic algorithm based + Levenberg-Marquardt non-linear regression'  ; write_msg
            PINPT%ls_type = 'GA'
          else
            write(message,'(A,A6,A)')'  LS_TYPE:  ',PINPT%ls_type,' is not defined or not available in the current version. Exit..'  ; write_msg
            stop
          endif

        case('PTOL') !set parameter tolerance for iteration step
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ptol
          write(message,'(A,F15.8)')'    P_TOL:  ',PINPT%ptol  ; write_msg

        case('FTOL') !set function tolerance for iteration step
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ftol
          write(message,'(A,F15.8)')'    F_TOL:  ',PINPT%ftol  ; write_msg

        case('FDIFF') !set function tolerance for fitting step (compared with previous fitting steps)
          read(inputline,*,iostat=i_continue) desc_str, PINPT%fdiff
          write(message,'(A,F15.8)')'   F_DIFF:  ',PINPT%fdiff  ; write_msg

        case('MITER') !set maximum iteration % maximum # of generations for GA
          if(.not. PINPT%flag_miter_parse) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%miter
            write(message,'(A,I8)')' MAX_ITER:  ',PINPT%miter  ; write_msg
          endif     

        case('MXFIT') !set maximum number of repeat for LMDIF after achieving local minima
          if(.not. PINPT%flag_mxfit_parse) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%mxfit
            write(message,'(A,I8)')'  MAX_FIT:  ',PINPT%mxfit  ; write_msg
          endif

        case('NPOP') ! set number of generation for genetic algorithm
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ga_npop
          write(message,'(A,I8)')'  GA_NPOP:  ',PINPT%ga_npop  ; write_msg

      end select

      return
   endsubroutine

   subroutine set_print_only_target(PINPT, inputline)
      type(incar  )  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_print_only_target'

      read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_only_target
      write(message,'(A,L)')' P_TARGET:',PINPT%flag_print_only_target  ; write_msg

      return
   endsubroutine

   subroutine set_onsite_tol(NN_TABLE, inputline)
      type(hopping)  ::  NN_TABLE
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_onsite_tol'

      read(inputline,*,iostat=i_continue) desc_str, desc_str, NN_TABLE%onsite_tolerance
      write(message,'(A,F16.4)')'ONSITETOL:',NN_TABLE%onsite_tolerance  ; write_msg

      return
   endsubroutine

   subroutine set_effective(PINPT, desc_str)
      type(incar)    ::  PINPT     
      integer*4          i_continue,i_dummy
      character*132      inputline, dummy_
      character*40       desc_str
      integer*4          nitems
      external           nitems
      real*8             init,fina
      character(*), parameter :: func = 'set_effective'

 set_eff: do while(trim(desc_str) .ne. 'END')
        read(pid_incar, '(A)', iostat = i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        call strip_off (inputline, dummy_, ' ', '#', 0) ! cut off unnecessary comments
        if(index(dummy_, trim(desc_str)) .ge. 1) then
          inputline = dummy_
        endif

        select case( trim(desc_str) )
          case('EFF_ORB') 
            call strip_off(inputline, PINPT%eff_orb_dummyc, trim(desc_str), ' ',2)
            write(message,'(A,A)')'  EFF_ORB: ',trim(PINPT%eff_orb_dummyc)  ; write_msg

          case('EFF_EWINDOW')
            call get_window(init,fina,inputline, desc_str)
            write(message,'(A,F10.5,A,F10.5,A)')' EFF_WIND: (EFF_EWINDOW) [', init,' : ', fina,']'  ; write_msg
            PINPT%eff_emin = init
            PINPT%eff_emax = fina
        endselect

      enddo set_eff

      PINPT%flag_get_effective_ham = .true.

      return
   endsubroutine

   subroutine set_gainp(PKAIA, desc_str)
      type(gainp)    ::  PKAIA
      integer*4          i_continue
      integer*4          i_dummy
      character*132      inputline, dummy_
      character*40       desc_str
      integer*4          nitems
      external           nitems
      character(*), parameter :: func = 'set_gainp'
     

      write(message,'(A,3F10.5,A)')'   * PARAMETERS FOR GENETIC ALGORITHM *'  ; write_msg
      write(message,'(A,3F10.5,A)')'   -----------------------------------'  ; write_msg

  set_pkaia: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        call strip_off (inputline, dummy_, ' ', '#', 0) ! cut off unnecessary comments
        if(index(dummy_, trim(desc_str)) .ge. 1) then
          inputline = dummy_
        endif

        select case ( trim(desc_str) )
!         case('LMDIF')
!            read(inputline,*,iostat=i_continue) desc_str, PKAIA%flag_ga_with_lmdif
!            write(message,'(A,L       )')'    LMDIF: ', PKAIA%flag_ga_with_lmdif  ; write_msg

          case('MGEN')
             read(inputline,*,iostat=i_continue) desc_str, PKAIA%mgen
             write(message,'(A,I8      )')'     MGEN: ', PKAIA%mgen  ; write_msg

          case('NPOP')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%npop
            write(message,'(A,I8      )')'     NPOP: ', PKAIA%npop  ; write_msg

          case('NGENE')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%ngene
            write(message,'(A,I8      )')'    NGENE: ', PKAIA%ngene  ; write_msg

          case('PCROSS')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pcross
            write(message,'(A, F10.5  )')'   PCROSS: ', PKAIA%pcross  ; write_msg

          case('RMUTMIN')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pmutmn
            write(message,'(A, F10.5  )')' MIN_MUTR: ', PKAIA%pmutmn  ; write_msg

          case('RMUTMAX')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pmutmx
            write(message,'(A, F10.5  )')' MAX_MUTR: ', PKAIA%pmutmx  ; write_msg

          case('RMUTINI')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pmut
            write(message,'(A, F10.5  )')' INI_MUTR: ', PKAIA%pmut  ; write_msg

          case('MUT_MOD')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%imut
            write(message,'(A, I8     )')' MUT_MODE: ', PKAIA%imut  ; write_msg

          case('FDIF')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%fdif
            write(message,'(A, F10.5  )')'     FDIF: ', PKAIA%fdif  ; write_msg

          case('IREP')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%irep
            write(message,'(A, I8     )')'     IREP: ', PKAIA%irep  ; write_msg

          case('IELITE')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%ielite
            write(message,'(A, I8     )')'   IELITE: ', PKAIA%ielite  ; write_msg

          case('VERBOSE')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%ivrb  
            write(message,'(A, I8     )')'  VERBOSE: ', PKAIA%ivrb    ; write_msg

          case('CONVTOL')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%convtol
            write(message,'(A, F10.5  )')'  CONVTOL: ', PKAIA%convtol  ; write_msg

          case('CONVWIN')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%convwin
            write(message,'(A, I8     )')'  CONVWIN: ', PKAIA%convwin  ; write_msg

          case('IGUESSF')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%iguessf
            write(message,'(A, F10.5  )')'  IGUESSF: ', PKAIA%iguessf  ; write_msg

          case('ISEED'  )
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%iseed  
            write(message,'(A, I8     )')' RANDSEED: ', PKAIA%iseed    ; write_msg

          case('UBOUND')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%upper_bound
            write(message,'(A, F10.5  )')'   UBOUND: ', PKAIA%upper_bound  ; write_msg

          case('LBOUND')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%lower_bound
            write(message,'(A, F10.5  )')'   LBOUND: ', PKAIA%lower_bound  ; write_msg

        end select
      enddo set_pkaia

      write(message,'(A,3F10.5,A)')'   * END READING: GENETIC ALGORITHM PARAMETERS*'  ; write_msg
      write(message,'(A,3F10.5,A)')'   ---------------------------------------------'  ; write_msg

      return
   endsubroutine

   subroutine set_energy_window(PINPT, inputline, desc_str)
      type(incar  )  ::  PINPT
      integer*4          i_continue
      integer*4          i_dummy
      character*132      inputline, dummy_
      character*40       desc_str, dummy, dummy1, dummy2
      integer*4          nitems
      external           nitems
      character(*), parameter :: func = 'set_energy_window'
      PINPT%feast_nemax = -1 ! default
      if(PINPT%flag_erange) then
        write(message,'(A)') '    !WARN! ERANGE tag and EWINDOW tag cannot be used simultaneously.'  ; write_msg
        write(message,'(A)') '           EWINDOW tag sets to use extended eigensolver routines based on '  ; write_msg
        write(message,'(A)') '           FEAST algorithm for the sparse matrix eigenvalue problem.'  ; write_msg
        write(message,'(A)') '           ERANGE tag, on the other hand, sets to use ZHEEVX function in LAPACK routines '  ; write_msg
        write(message,'(A)') '           for the dense matrix eigenvalue problem.'  ; write_msg
        write(message,'(A)') '           Please deactivate one of the tag to proceed.'  ; write_msg
        write(message,'(A)') '           If you dealing with very large system for the band structure calculation,'  ; write_msg
        write(message,'(A)') '           it is recommended to use EWINDOW tag rather than ERANGE'  ; write_msg
        write(message,'(A)') '           Note that EWINDOW does not applies to the post processings requested by'  ; write_msg
        write(message,'(A)') '           various "SET" tags, for example, Z2_INDEX, WCC, ZAK_PHASE, PARITY_CHECK, etc.'  ; write_msg
        write(message,'(A)') '           In these routines, ERANGE tags there in will be used and ZHEEVX function will be'  ; write_msg
        write(message,'(A)') '           used instead.'  ; write_msg
        write(message,'(A)') '           Exit program anyway...'  ; write_msg
        stop
      endif

      call strip_off (inputline, dummy_, ' ', '#', 0) ! cut off unnecessary comments
      if(index(dummy_, trim(desc_str)) .ge. 1) then
        inputline = dummy_
      endif

      i_dummy = index(inputline,'NE_MAX')
      if(i_dummy .ge. 1) then 
        call strip_off (inputline, dummy_, 'NE_MAX', ' ', 2)
        call str2int(trim(dummy_), PINPT%feast_nemax)
        call strip_off (inputline, dummy_, ' ', 'NE_MAX', 0)
        inputline = dummy_
      elseif(i_dummy .eq. 0) then
        i_dummy = index(inputline,'NE_MAX')
        if(i_dummy .ge. 1) then
          call strip_off (inputline, dummy_, 'NE_MAX', ' ', 2)
          call str2int(trim(dummy_), PINPT%feast_nemax)
        call strip_off (inputline, dummy_, ' ', 'NE_MAX', 0)
          inputline = dummy_
        endif
      endif

      call strip_off (inputline, dummy, trim(desc_str), ' ', 2)

      i_dummy = index(dummy,':')
      
      if(i_dummy .eq. 0) then
        if( nitems(inputline) -1 .eq. 0) then
          write(message,'(A,A)')'  !WARN! No energy window has been defined! Please check EWINDOW tag atain. Exit... ', func  ; write_msg
          stop
        elseif( nitems(inputline) -1 .eq. 2) then
          read(inputline,*,iostat=i_continue) desc_str,PINPT%feast_emin, PINPT%feast_emax
          write(message,'(A,F12.5,A,F12.5,A)')' E_WINDOW:  [',PINPT%feast_emin,',',PINPT%feast_emax,'], [emin:emax] (eV)'  ; write_msg
          PINPT%flag_sparse = .true. ! call FEAST algorithom
        else
          write(message,'(A,A)')'  !WARN! Please check EWINDOW tag syntax. Exit... Current setting = ',trim(dummy)  ; write_msg
          write(message,'(A,A)')'         Only accept following syntax: EWINDOW EMIN:EMAX  or'  ; write_msg
          write(message,'(A,A)')'                                       EWINDOW EMIN EMAX '  ; write_msg
          write(message,'(A,A)')'         Here, EMIN and EMAX are real*8 values for the eigenvalues to be searched.'  ; write_msg
          stop
        endif

      elseif(i_dummy .ge. 1) then
        call strip_off(dummy, dummy1,' ',':',0)
        call strip_off(dummy, dummy2,':',' ',2)
        if(len_trim(dummy1) .eq. 0 .or. len_trim(dummy2) .eq. 0) then
          write(message,'(A,A)')'  !WARN! Please check ERANGE tag syntax. Exit... Current setting = ',trim(dummy)  ; write_msg
          write(message,'(A,A)')'         Only accept following syntax: ERANGE INIT:FINA  or'  ; write_msg
          write(message,'(A,A)')'                                       ERANGE INIT FINA'  ; write_msg
          write(message,'(A,A)')'         Here, INIT and FINA are integer values for the eigenvalue index.'  ; write_msg
          stop
        endif
        call str2real(dummy1,PINPT%feast_emin)
        call str2real(dummy2,PINPT%feast_emax)
        write(message,'(A,F12.5,A,F12.5,A)')' E_WINDOW:  [',PINPT%feast_emin,',',PINPT%feast_emax,'], [emin:emax] (eV)'  ; write_msg
        PINPT%flag_sparse = .true. ! call FEAST algorithom
      endif

      if(PINPT%feast_nemax .ge. 1) then
        ! feast_nemax is upper bound of states to be stored in your output. 
        ! The states, beyond this criteria will be thrown away.
        write(message,'(A,I0,A)')'   NE_MAX: ', PINPT%feast_nemax, ' (maximum # of states within [emin:emax])'  ; write_msg
      elseif(PINPT%feast_nemax .le. 0) then
        write(message,'(A,I0,A)')'   NE_MAX: N_ORBIT * ISPINOR [N_ORBIT:number of orbitals ]'  ; write_msg
        write(message,'(A)')     '                             [ISPINOR:2 for LSORB=.TRUE. ]'  ; write_msg
        write(message,'(A)')     '                             [       :1 for LSORB=.FALSE.]'  ; write_msg
      endif

      return
   endsubroutine

   subroutine set_energy_range(PINPT, inputline, desc_str)
      type(incar  )  ::  PINPT
      integer*4     i_continue
      integer*4     i_dummy
      character*132 inputline
      character*40  desc_str, dummy, dummy1, dummy2
      integer*4     nitems
      external      nitems
      character(*), parameter :: func = 'set_energy_range'

      if(PINPT%flag_sparse) then
        write(message,'(A)') '    !WARN! ERANGE tag and EWINDOW tag cannot be used simultaneously.'  ; write_msg
        write(message,'(A)') '           EWINDOW tag sets to use extended eigensolver routines based on '  ; write_msg
        write(message,'(A)') '           FEAST algorithm for the sparse matrix eigenvalue problem.'  ; write_msg
        write(message,'(A)') '           ERANGE tag, on the other hand, sets to use ZHEEVX function in LAPACK routines '  ; write_msg
        write(message,'(A)') '           for the dense matrix eigenvalue problem.'  ; write_msg
        write(message,'(A)') '           Please deactivate one of the tag to proceed.'  ; write_msg
        write(message,'(A)') '           If you dealing with very large system for the band structure calculation,'  ; write_msg
        write(message,'(A)') '           it is recommended to use EWINDOW tag rather than ERANGE'  ; write_msg
        write(message,'(A)') '           Note that EWINDOW does not applies to the post processings requested by'  ; write_msg
        write(message,'(A)') '           various "SET" tags, for example, Z2_INDEX, WCC, ZAK_PHASE, PARITY_CHECK, etc.'  ; write_msg
        write(message,'(A)') '           In these routines, ERANGE tags there in will be used and ZHEEVX function will be'  ; write_msg
        write(message,'(A)') '           used instead.'  ; write_msg
        write(message,'(A)') '           Exit program anyway...'  ; write_msg
        stop
      endif

      call strip_off (inputline, dummy,trim(desc_str),' ',2)
      i_dummy=index(dummy,':')
      if(i_dummy .eq. 0) then
        if( nitems(inputline) -1 .eq. 0) then
          write(message,'(A,A)')'  !WARN! No range has been defined! Please check ERANGE tag atain. Exit... ', func  ; write_msg
          stop
        elseif( nitems(inputline) -1 .eq. 1) then
          read(inputline,*,iostat=i_continue) desc_str,PINPT%init_erange
          PINPT%fina_erange = PINPT%init_erange
          write(message,'(A,I6,A,I6,A)')'  E_RANGE:  FROM ',PINPT%init_erange,'-th TO ',PINPT%fina_erange,'-th EIGENSTATES'  ; write_msg
          PINPT%flag_erange = .true.
        elseif( nitems(inputline) -1 .eq. 2) then
          read(inputline,*,iostat=i_continue) desc_str,PINPT%init_erange, PINPT%fina_erange
          write(message,'(A,I6,A,I6,A)')'  E_RANGE:  FROM ',PINPT%init_erange,'-th TO ',PINPT%fina_erange,'-th EIGENSTATES'  ; write_msg
          PINPT%flag_erange = .true.
        else
          write(message,'(A,A)')'  !WARN! Please check ERANGE tag syntax. Exit... Current setting = ',trim(dummy)  ; write_msg
          write(message,'(A,A)')'         Only accept following syntax: ERANGE INIT:FINA  or'  ; write_msg
          write(message,'(A,A)')'                                       ERANGE INIT FINA'  ; write_msg
          write(message,'(A,A)')'         Here, INIT and FINA are integer values for the eigenvalue index.'  ; write_msg
          stop
        endif

      elseif(i_dummy .ge. 1) then
        call strip_off(dummy, dummy1,' ',':',0)
        call strip_off(dummy, dummy2,':',' ',2)
        if(len_trim(dummy1) .eq. 0 .or. len_trim(dummy2) .eq. 0) then
          write(message,'(A,A)')'  !WARN! Please check ERANGE tag syntax. Exit... Current setting = ',trim(dummy)  ; write_msg
          write(message,'(A,A)')'         Only accept following syntax: ERANGE INIT:FINA  or'  ; write_msg
          write(message,'(A,A)')'                                       ERANGE INIT FINA'  ; write_msg
          write(message,'(A,A)')'         Here, INIT and FINA are integer values for the eigenvalue index.'  ; write_msg
          stop
        endif
        call str2int(dummy1,PINPT%init_erange)
        call str2int(dummy2,PINPT%fina_erange)
        write(message,'(A,I6,A,I6,A)')'  E_RANGE:  FROM ',PINPT%init_erange,'-th TO ',PINPT%fina_erange,'-th EIGENSTATES'  ; write_msg
        PINPT%flag_erange = .true.
      endif
      PINPT%nband = PINPT%fina_erange - PINPT%init_erange + 1

      return
   endsubroutine

   subroutine set_get_band(PINPT, inputline, desc_str)
      implicit none
      type(incar) :: PINPT
      character*132  inputline
      character*40   desc_str
      integer*4      i_continue

      read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_get_band
      if(PINPT%flag_get_band) then
        write(message,'(A)')'  GET_BND: .TRUE.'  ; write_msg
      elseif(.not.PINPT%flag_get_band) then
        write(message,'(A)')'  GET_BND: .FALSE.'  ; write_msg
      endif
      return
   endsubroutine

   subroutine set_nn_max(PINPT, inputline,desc_str)
      implicit none
      type(incar) :: PINPT
      character*132  inputline
      character*40   desc_str, dummy_, dummy
      integer*4      i_continue, i_dummy
      integer*4      nitems
      external       nitems


      call strip_off (inputline, dummy_, ' ', '#', 0) ! cut off unnecessary comments
      if(index(dummy_, trim(desc_str)) .ge. 1) then
        inputline = dummy_
      endif

      if(nitems(inputline) -1 .eq. 3) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%nn_max(1:3)
        write(message,'(A,3I5)')'   NN_MAX:', PINPT%nn_max(1:3)  ; write_msg
      elseif(nitems(inputline) -1 .eq. 1) then
        read(inputline,*,iostat=i_continue) desc_str, i_dummy
        PINPT%nn_max(1:3) = i_dummy
        write(message,'(A,3I5)')'   NN_MAX:', PINPT%nn_max(1:3)  ; write_msg
      else
        write(message,'(A)')'    !WARN! Error in reading NN_MAX tag. Exit program...'  ; write_msg
        stop
      endif

      return
   endsubroutine

#ifdef SPGLIB
   subroutine set_spglib_write(PINPT, inputline, desc_str)
      implicit none
      type(incar) :: PINPT
      character*132  inputline
      character*40   desc_str
      integer*4      i_continue

      read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_spglib

      if(PINPT%flag_spglib) then
        write(message,'(A)')'  LSPGLIB: .TRUE.'  ; write_msg
      elseif(.not.PINPT%flag_spglib) then
        write(message,'(A)')'  LSPGLIB: .FALSE.'  ; write_msg
      endif

      return
   endsubroutine
#endif


endmodule
