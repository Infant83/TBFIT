#include "alias.inc"
module read_incar
   use mpi_setup
   use parameters
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
                  if_main write(6,'(A,*(I6))')' EIG_PLOT:  ',PINPT%i_eig_print
                  cycle plot_eig

                case('IKPT')
                  PINPT%n_kpt_print = nitems(inputline) - 1
                  allocate(PINPT%i_kpt_print(PINPT%n_kpt_print))
                  read(inputline,*,iostat=i_continue) desc_str,PINPT%i_kpt_print(1:PINPT%n_kpt_print)
                  if_main write(6,'(A,*(I6))')' KPT_PLOT:  ',PINPT%i_kpt_print
                  cycle plot_eig

                case('NGRID')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 3) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT%ngrid(1:3)
                    PINPT%flag_default_ngrid = .false.
                  else
                    if_main write(6,'(A)')'    !WARN! NGRID tag of "SET EIGPLOT" should be three consequent integer numbers.'
                    if_main write(6,'(A,A)')'           Please check NGRID tag again. Exit... ',func
                    stop
                  endif

                case('RORIGIN')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 3) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT%r_origin(1:3)
                    if_main  write(6,'(A,3(F15.8))')'   R_ORIG:  ',PINPT%r_origin(1:3)
                    PINPT%flag_default_rorigin = .false.
                  else
                    if_main write(6,'(A)')'    !WARN! RORIGIN tag of "SET EIGPLOT" should be three consequent real values.'
                    if_main write(6,'(A,A)')'           Please check RORIGIN tag again. Exit... ',func
                    stop
                  endif

                case('WAVEPLOT','WAV_PLOT')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_plot_wavefunction
                    if_main write(6,'(A,L)')' WAV_PLOT:  ',PINPT%flag_plot_wavefunction
                  else
                    if_main write(6,'(A)')'    !WARN! WAV_PLOT tag of "SET EIGPLOT" should be .TRUE. or .FALSE.'
                    if_main write(6,'(A,A)')'           Please check WAV_PLOT tag again. Exit... ',func
                    stop
                  endif

                case('RCUT', 'RCUT_ORB')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str, PINPT%rcut_orb_plot
                    if_main write(6,'(A,F9.4)')' RCUT_ORB:  ',PINPT%rcut_orb_plot
                  else
                    if_main write(6,'(A)')'    !WARN! RCUT_ORB tag of "SET EIGPLOT" should be single real type parameter.'
                    if_main write(6,'(A,A)')'           Please check RCUT_ORB tag again. Exit... ',func
                  endif

              end select eig_mode

            enddo plot_eig

        case('STMPLOT')
            PINPT%flag_plot_stm_image = .true.
            PINPT%flag_repeat_cell_orb_plot(1:3) = .true.

            if_main write(6,'(A,L)')' STM_PLOT:  ',PINPT%flag_plot_stm_image
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
                    if_main write(6,'(A)')'    !WARN! NGRID tag of "SET STMPLOT" should be three consequent integer numbers.'
                    if_main write(6,'(A,A)')'           Please check NGRID tag again. Exit... ',func
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
                      if_main write(6,'(A,I2,A,F9.4,A,F9.4,A)')' STM_ERAN: PLOT(',PINPT%n_stm,')= [', &
                                                                       stm_emin(PINPT%n_stm),':',stm_emax(PINPT%n_stm),']'
                    else
                      if_main write(6,'(A)')'    !WARNING!  STM_ERANGE is not properly set up.'
                      if_main write(6,'(A)')'    !WARNING!  Proper usage is as follows:'
                      if_main write(6,'(A)')'    !WARNING!    STM_ERANGE  EMIN:EMAX or EMIN EMAX'
                      if_main write(6,'(A)')'    !WARNING!  Exit program...'
                      stop
                    endif
                  elseif(i_dummy .gt. 1) then ! if ':' is provided
                    call strip_off (trim(dummy), dummy1,' ',':',0)
                    call str2real(dummy1,stm_emin(PINPT%n_stm))
                    call strip_off (trim(dummy), dummy2,':',' ',2)
                    call str2real(dummy2,stm_emax(PINPT%n_stm))
                    if_main write(6,'(A,I2,A,F9.4,A,F9.4,A)')' STM_ERAN: PLOT(',PINPT%n_stm,')= [', &
                                                                     stm_emin(PINPT%n_stm),':',stm_emax(PINPT%n_stm),']'
                  endif

                case('RCUT', 'RCUT_ORB')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str, PINPT%rcut_orb_plot
                    if_main write(6,'(A,F9.4)')' RCUT_ORB:  ',PINPT%rcut_orb_plot
                  else
                    if_main write(6,'(A)')'    !WARN! RCUT_ORB tag of "SET EIGPLOT" should be single real type parameter.'
                    if_main write(6,'(A,A)')'           Please check RCUT_ORB tag again. Exit... ',func
                  endif
                case('REPEAT_CELL')
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .ne. 3) then
                    if_main write(6,'(A)')'   !WARN! REPEAT_CELL only accepts three logical arguments, for example "T T F",'
                    if_main write(6,'(A)')'          which implies that the orbitals along a1 and a2 direction are repeated and'
                    if_main write(6,'(A)')'          not for the a3 direction. check manual. Stop program...'
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
          write(6,'(A,I1,A)')'   !WARN! You have set the orbital along a',i,'-direction will not'
          write(6,'(A)')     '          be repeated in the STM or EIGPLOT. Proceed anyway...'
        endif
      enddo


      return    
   endsubroutine

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
        if_main write(6,'(A,A,A)')'    !WARNING! ',trim(PINPT%pfilenm),' is alread provided,'
        if_main write(6,'(A,A,A)')'    !WARNING! ','but the TBPARAM is also provided in the INCAR-TB file.'
        if_main write(6,'(A,A,A)')'    !WARNING! ','Those values from INCAR-TB will not be read'
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
            param_const(2,1:PINPT%nparam) = 99999d0
            param_const(3,1:PINPT%nparam) =-99999d0
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
          if_main write(6,'(A,A,A)')'    !WARNING! TBPARAM tag is set in the INCAR-TB but'
          if_main write(6,'(A,A,A)')'    !WARNING! TB-parameter is not provided.         '
          if_main write(6,'(A,A,A)')'    !WARNING! TB-parameter will be read from externl file if set by PFILE'
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
              if_main write(6,'(A,A12,2(A,1F8.4))')'  NN_PAIR:  ',strip_nn_pair_(i), &
                                                           ' R0_max:',strip_nn_dist_(i), &
                                                           '   R0:',strip_nn_r0_(i)
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

   subroutine set_density_of_states(PINPT, PINPT_DOS, desc_str)
      type(incar )  ::  PINPT
      type(dos)     :: PINPT_DOS
!     integer*4, parameter :: max_dummy = 9999999
      integer*4     i, ii, i_continue
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
            if_main write(6,'(A)')'  GET_DOS: .TRUE.'
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
            PINPT_DOS%dos_flag_print_eigen = .false.
            PINPT_DOS%dos_kunit = 'R'
            PINPT_DOS%dos_iband = 1
            PINPT_DOS%dos_fband = 999999
            PINPT_DOS%dos_smearing = 0.025

   set_dos: do while(trim(desc_str) .ne. 'END')
              read(pid_incar,'(A)',iostat=i_continue) inputline
              read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
              if(i_continue .ne. 0) cycle      ! skip empty line
              if(desc_str(1:1).eq.'#') cycle   ! skip comment
              if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

    case_dos: select case ( trim(desc_str) )
                case('NEDOS')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_nediv
                  if_main write(6,'(A,I8)')' DOS_EDIV:', PINPT_DOS%dos_nediv
                case('DOS_ERANGE')
                  call strip_off (trim(inputline), dummy, 'DOS_ERANGE', ' ' , 2)   ! get dos_range
                  i_dummy=index(dummy,':')
                  call strip_off (trim(dummy), dummy1,' ',':',0)
                  if( i_dummy .eq. 0) then
                    i_dummy2 = nitems(dummy)
                    if(i_dummy2 .eq. 2)then
                      read(dummy,*,iostat=i_continue) PINPT_DOS%dos_emin,PINPT_DOS%dos_emax
                      if_main write(6,'(A,F15.8)')' DOS_EMIN:  ',PINPT_DOS%dos_emin
                      if_main write(6,'(A,F15.8)')' DOS_EMAX:  ',PINPT_DOS%dos_emax
                    else
                      if_main write(6,'(A)')'    !WARNING!  DOS_ERANGE is not properly set up.'
                      if_main write(6,'(A)')'    !WARNING!  Proper usage is as follows:'
                      if_main write(6,'(A)')'    !WARNING!    DOS_ERANGE  EMIN:EMAX , :EMAX, EMIN: or :'
                      if_main write(6,'(A)')'    !WARNING! or DOS_ERANGE  EMIN EMAX'
                      if_main write(6,'(A)')'    !WARNING!  Exit program...'
                      stop
                    endif
                  elseif(i_dummy .ge. 1) then
                    if(len_trim(dummy1) .eq. 0) then
                      PINPT_DOS%dos_emin = -10.0d0 ! default dos_emin
                      if_main write(6,'(A,F15.8)')' DOS_EMIN:  ',PINPT_DOS%dos_emin
                    else
                      call str2real(dummy1,PINPT_DOS%dos_emin)
                      if_main write(6,'(A,F15.8)')' DOS_EMIN:  ',PINPT_DOS%dos_emin
                    endif
                    call strip_off (trim(dummy), dummy2,':',' ',2)
                    if(len_trim(dummy2) .eq. 0) then
                      PINPT_DOS%dos_emax =  10.0d0 ! default dos_emax
                      if_main write(6,'(A,F15.8)')' DOS_EMAX:  ',PINPT_DOS%dos_emax
                    else
                      call str2real(dummy2,PINPT_DOS%dos_emax)
                      if_main write(6,'(A,F15.8)')' DOS_EMAX:  ',PINPT_DOS%dos_emax
                    endif
                  endif
                case('DOS_NRANGE')
                  call strip_off (trim(inputline), dummy, 'DOS_NRANGE', ' ' , 2)   ! get dos_range
                  i_dummy=index(dummy,':')
                  call strip_off (trim(dummy), dummy1,' ',':',0)
                  if( i_dummy .eq. 0) then
                    i_dummy2 = nitems(dummy)
                    if(i_dummy2 .eq. 2)then
                      read(dummy,*,iostat=i_continue) PINPT_DOS%dos_iband,PINPT_DOS%dos_fband
                      if_main write(6,'(A,I8)')' DOS_IBND:',PINPT_DOS%dos_iband
                      if_main write(6,'(A,I8)')' DOS_FBND:',PINPT_DOS%dos_fband
                    else
                      if_main write(6,'(A)')'    !WARNING!  DOS_NRANGE is not properly set up.'
                      if_main write(6,'(A)')'    !WARNING!  Proper usage is as follows:'
                      if_main write(6,'(A)')'    !WARNING!    DOS_NRANGE  IBAND:FBAND , :FBAND, IBAND: or :'
                      if_main write(6,'(A)')'    !WARNING! or DOS_NRANGE  IBAND FBAND'
                      if_main write(6,'(A)')'    !WARNING!  Exit program...'
                      stop
                    endif
                  elseif(i_dummy .ge. 1 .and. len_trim(dummy) .ne. 1) then
                    if(len_trim(dummy1) .eq. 0) then
                      PINPT_DOS%dos_iband = 1 ! default dos_emin
                      if_main write(6,'(A,I8)')' DOS_IBND:',PINPT_DOS%dos_iband
                    else
                      call str2int(dummy1,PINPT_DOS%dos_iband)
                      if_main write(6,'(A,I8)')' DOS_IBND:',PINPT_DOS%dos_iband
                    endif
                    call strip_off (trim(dummy), dummy2,':',' ',2)
                    if(len_trim(dummy2) .eq. 0) then
                      PINPT_DOS%dos_fband =  999999 ! default dos_emax
                      if_main write(6,'(A)')' DOS_FBND: automatically set to NEIG'
                    else
                      if(flag_number( trim(dummy2) ) ) then
                        call str2int(dummy2,PINPT_DOS%dos_fband)
                        if_main write(6,'(A,I8)')' DOS_FBND:',PINPT_DOS%dos_fband
                      elseif(.not. flag_number( trim(dummy2) ) .and. trim(dummy2) .eq. 'NEIG' ) then
                        PINPT_DOS%dos_fband =  999999
                        if_main write(6,'(A)')' DOS_FBND: NEIG = N_ORBIT'
                      else
                        if_main write(6,'(A,I8)')'  !WARNING! DOS_FBND setting is inproper. Please check syntax again. Exit..'
                        stop
                      endif
                    endif
                  elseif(i_dummy .eq. 1 .and. len_trim(dummy) .eq. 1) then
                    PINPT_DOS%dos_iband = 1
                    if_main write(6,'(A,I8)')' DOS_IBND: 1'
                    if_main write(6,'(A,I8)')' DOS_FBND: NEIG'
                  endif
                case('MKGRID')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kgrid(1:3)
                  if_main write(6,'(A,3I6)')' DOS_KDIV: Monkhorst-Pack grid',PINPT_DOS%dos_kgrid(1:3)
                  PINPT_DOS%dos_flag_gamma=.false.
                case('GKGRID')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kgrid(1:3)
                  if_main write(6,'(A,3I6)')' DOS_KDIV: Gamma-centered grid',PINPT_DOS%dos_kgrid(1:3)
                  PINPT_DOS%dos_flag_gamma=.true.
                case('KGRID')
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kgrid(1:3)
                  if_main write(6,'(A,3I6)')' DOS_KDIV: Monkhorst-Pack grid',PINPT_DOS%dos_kgrid(1:3)
                  PINPT_DOS%dos_flag_gamma=.false.
                case('KSHIFT') ! optional shift of the mesh (s_1, s_2, s_3) ; the usage is same as VASP 
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_kshift(1:3)
                  if_main write(6,'(A,3F6.2)')' DOS_KSFT:',PINPT_DOS%dos_kshift(1:3)
                case('PRINT_KPTS') ! print kpoint information (reciprocal unit) into file
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_kpoint
                  elseif(i_dummy .eq. 2) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_kpoint, PINPT_DOS%dos_kfilenm
                  endif
                  if(PINPT_DOS%dos_flag_print_kpoint) then
                    if_main write(6,'(A,A)')' DOS_KOUT: KPOINT OUT -> .TRUE. ',trim(PINPT_DOS%dos_kfilenm)
                  elseif(.not. PINPT_DOS%dos_flag_print_kpoint) then
                    if_main write(6,'(A,A)')' DOS_KOUT: KPOINT OUT -> .FALSE.'
                  endif
                case('PRINT_EIG') ! print eigenstate energy information into ENSURF.EIG.NEIG.dat file
                  i_dummy = nitems(inputline) - 1
                  if(i_dummy .eq. 1) then
                    read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_flag_print_eigen
                    if(PINPT_DOS%dos_flag_print_eigen) then
                      if_main write(6,'(A)')'  !WARNING! DOS_EOUT: EIGEN  OUT -> .TRUE. but the target energy to be plotted is not specified.'
                      if_main write(6,'(A)')'  !WARNING! Please check PRINT_EIG tag again. The proper use is : PRINT_EIG .TRUE. N1 N2 N3... or N1:N3 '
                      if_main write(6,'(A)')'  !WARNING! Exit program...'
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
                        if_main write(6,'(A,*(I5))')' DOS_EOUT: ENSURF OUT -> .TRUE.',PINPT_DOS%dos_ensurf(1:PINPT_DOS%dos_n_ensurf)
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
                            i_dummyr(ii:ii+i_dummy5 - i_dummy4) = (/i_dummy4:i_dummy5/)
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
                    if_main write(6,'(A)')' DOS_UNIT: KPOINT UNIT = Reciprocal'
                  elseif(PINPT_DOS%dos_kunit .eq. 'a' .or. PINPT_DOS%dos_kunit .eq. 'A') then
                    PINPT_DOS%dos_kunit = 'A'
                    if_main write(6,'(A)')' DOS_UNIT: KPOINT UNIT = Angstrom unit (1/A)'
                  endif
                case('DOS_FNAME') ! DOS output file name
                  read(inputline,*,iostat=i_continue) desc_str,desc_str
                  PINPT_DOS%dos_filenm = trim(desc_str)
                  if_main write(6,'(A,A)')' DOS_FNAM: ',trim(PINPT_DOS%dos_filenm)
                case('SMEARING') ! gaussian smearing
                  read(inputline,*,iostat=i_continue) desc_str,PINPT_DOS%dos_smearing
                  if_main write(6,'(A,F8.4)')'DOS_SIGMA: GAUSSIAN WIDTH = ',PINPT_DOS%dos_smearing

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
      if_main write(6,'(A)')' SET_RIBN: SET UP RIBBON GEOMETRY = .TRUE.'
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
            if_main write(6,'(A,A)')' KPTS_FNM:  (for ribbon) ',trim(PINPT%kfilenm)

        end select
      enddo set_rib

      if(minval(PINPT%ribbon_nslab(1:3)) .le. 0) then
        if_main write(6,'(A)')'  !WARNING! NSLAB tag for the RIBBON calculation is not properly defined. Please check the tag.'
        if_main write(6,'(A)')'  !WARNING! Proper usage: NSLAB   N1 N2 N3'
        if_main write(6,'(A)')'  !WARNING! Exit program...'
        stop
      endif
      if(minval(PINPT%ribbon_nslab(1:3)) .lt. 0) then
        if_main write(6,'(A)')'  !WARNING! VACUUM tag for the RIBBON calculation is not properly defined. Please check the tag.'
        if_main write(6,'(A)')'  !WARNING! Proper usage: VACUUM   vac_1 vac_2 vac_3'
        if_main write(6,'(A)')'  !WARNING! Exit program...'
        stop
      endif

      if_main write(6,'(A,3I5)')' RIB_SLAB: (N1*A1,N2*A2,N3*A3) => N1,N2,N3 =', PINPT%ribbon_nslab(1:3)
      if_main write(6,'(A,3F12.6)')' RIB_VACU: (in Angstrom)', PINPT%ribbon_vacuum(1:3)
      if(PINPT%flag_print_only_ribbon_geom) then
        if_main write(6,'(A,3F12.6)')' RIB_GEOM: PRINT_ONLY = .TRUE.'
      else
        if_main write(6,'(A,3F12.6)')' RIB_GEOM: PRINT_ONLY = .FALSE.'
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
            if_main write(6,'(A,3F12.6)')'  E_FIELD:  ',PINPT%efield(1:3)
          case('EF_CENTER','EF_ORIGIN','EF_ORIGIN_FRAC','EF_CENTER_FRAC') ! field_origin
            PINPT%flag_efield_frac = .true.
            PINPT%flag_efield_cart = .false.
            read(inputline,*,iostat=i_continue) desc_str,PINPT%efield_origin(1:3)
!           if_main write(6,'(A,3F12.6)')'EF_ORIGIN:  (in factional coord) ',PINPT%efield_origin(1:3)
          case('EF_CENTER_CART','EF_ORIGIN_CART','EF_CORIGIN') ! field_origin
            PINPT%flag_efield_frac = .false.
            PINPT%flag_efield_cart = .true.
            read(inputline,*,iostat=i_continue) desc_str,PINPT%efield_origin_cart(1:3)
!           if_main write(6,'(A,3F12.6)')'EF_ORIGIN:  (in cartesian coord) ',PINPT%efield_origin_cart(1:3)
        end select
      enddo 

      return
   endsubroutine

   subroutine set_weight_factor(PINPT, PWGHT, desc_str)
      type(incar)  ::  PINPT
      type(weight)  :: PWGHT
      integer*4     i, i_orb
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

      i=0
      i_orb = 0
      do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(index(inputline,'#') .gt. 1) then
          call strip_off (trim(inputline), inputline_dummy, '', '#', 0) ! strip-off '#' comments
          inputline = inputline_dummy
        endif
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'
        if( index(inputline,'ORBT_I') .eq. 0 ) then
          i=i+1
          call strip_off (trim(inputline), strip_kp_(i), 'KRANGE', 'TBABND', 1)   ! get KRANGE strip
          call strip_off (trim(inputline), strip_tb_(i), 'TBABND', 'DFTBND', 1)   ! get TBABND strip
          call strip_off (trim(inputline), strip_df_(i), 'DFTBND', 'WEIGHT', 1)   ! get DFTBND strip
          call strip_off (trim(inputline), strip_wt_(i), 'WEIGHT', ' '     , 2)   ! get WEIGHT strip
        elseif( index(inputline,'ORBT_I') .gt. 0 ) then
          i_orb = i_orb + 1
          call strip_off (trim(inputline), strip_kp_orb_(i_orb), 'KRANGE', 'TBABND', 1)   ! get KRANGE strip
          call strip_off (trim(inputline), strip_tb_orb_(i_orb), 'TBABND', 'ORBT_I', 1)   ! get TBABND strip
          call strip_off (trim(inputline), strip_orb_(i_orb),    'ORBT_I', 'SITE_I', 1)   ! get ORBT_I strip
          call strip_off (trim(inputline), strip_site_(i_orb),   'SITE_I', 'PENALTY', 1)   ! get SITE_I strip
          strip_site_(i_orb) = str2lowcase(strip_site_(i_orb))
          call strip_off (trim(inputline), strip_pen_orb_(i_orb), 'PENALTY', ' '     , 2)   ! get PENALTY strip
        endif
      enddo
      if(nitems(inputline) .eq. 3) then
        read(inputline,*,iostat=i_continue) desc_str, desc_str, desc_str
        if(desc_str .eq. 'PRINT_ONLY') PINPT%flag_print_only_target = .true.
      endif
      PWGHT%nweight=i
      PINPT%nweight=i
      PWGHT%npenalty_orb = i_orb
      PINPT%npenalty_orb = i_orb
      if(PWGHT%nweight .eq. 0) then
        PWGHT%flag_weight_default = .true.
      elseif( PWGHT%nweight .ge. 1) then
        PWGHT%flag_weight_default = .false.
        if_main write(6,'(A,I8)')' N_CONSTR:',PWGHT%nweight
        allocate( PINPT%strip_kp(PWGHT%nweight) )
        allocate( PINPT%strip_tb(PWGHT%nweight) )
        allocate( PINPT%strip_df(PWGHT%nweight) )
        allocate( PINPT%strip_wt(PWGHT%nweight) )
        do i = 1, PWGHT%nweight
          PINPT%strip_kp(i)=strip_kp_(i)
          PINPT%strip_tb(i)=strip_tb_(i)
          PINPT%strip_df(i)=strip_df_(i)
          PINPT%strip_wt(i)=strip_wt_(i)
        enddo
      endif

      if(PWGHT%npenalty_orb .eq. 0) then
        PWGHT%flag_weight_default_orb = .true.
      elseif( PWGHT%npenalty_orb .ge. 1) then
        PWGHT%flag_weight_default_orb = .false.
        if_main write(6,'(A,I8,A)')' N_CONSTR:',PWGHT%npenalty_orb,' (# of penalty weight for orbital)'
        allocate( PINPT%strip_kp_orb(PWGHT%npenalty_orb) )
        allocate( PINPT%strip_tb_orb(PWGHT%npenalty_orb) )
        allocate( PINPT%strip_pen_orb(PWGHT%npenalty_orb) )
        allocate( PINPT%strip_orb(PWGHT%npenalty_orb) )
        allocate( PINPT%strip_site(PWGHT%npenalty_orb) )
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

      if_main write(6,'(A)')'   GET_Z2: .TRUE. '
  set_wann: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('Z2_ERANGE')  ! get z2_range
            call strip_off (trim(inputline), PINPT_BERRY%strip_z2_range, 'Z2_ERANGE', ' ' , 2)
            if_main write(6,'(3A)')'  Z2_RANG:  ',adjustl(trim(PINPT_BERRY%strip_z2_range)), &
                                           ' (ENERGY RANGE TO BE CALCULATED)'

          case('Z2_DIMENSION')  ! which direction? 
            call strip_off (trim(inputline), dummy, 'Z2_DIMENSION', ' ' , 2)
            i_dummy=index(dummy,'D')
            if(i_dummy .eq. 0) then
              if_main write(6,'(A)')' !WARN!: Z2_DIMENSION should be start with one of following strings:'
              if_main write(6,'(A)')'  3D  => Z2_DIMENSION 3D'
              if_main write(6,'(A)')'  2D  => Z2_DIMENSION 2D:z # z is normal direction of surface'
              if_main write(6,'(A)')'  1D  => Z2_DIMENSION 1D:x # x is parallel direction of 1D system'
              if_main write(6,'(A)')'  For 2D and 1D, "x", "y", or "z" indicates the normal and parallel'       
              if_main write(6,'(A)')'  direction of the system'
              if_main write(6,'(A)')'  Exit program anyway... Check your "SET Z2" in "INCAR-TB" again.'
              stop
            endif

            call str2int(dummy(1:i_dummy-1),PINPT_BERRY%z2_dimension)

            if(PINPT_BERRY%z2_dimension .eq. 3) then
              deallocate(PINPT_BERRY%z2_axis)
              allocate(PINPT_BERRY%z2_axis(3))
              PINPT_BERRY%z2_axis(:) = (/1,2,3/)
              if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy), ' mode (three direction will be checked)'

            elseif(PINPT_BERRY%z2_dimension .eq. 2) then
              deallocate(PINPT_BERRY%z2_axis)
              allocate(PINPT_BERRY%z2_axis(1))
              i_dummy = index(dummy,':')
              if(i_dummy .eq. 0) then
                PINPT_BERRY%z2_axis(1) = 3 ! default : third reciprocal axis will be chosen for the WCC surface normal axis
                if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy), ' mode (three direction will be checked)'
                if_main write(6,'( A)')'  Z2_AXIS: B3 direction of the reciprocal lattice vector (default)'
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','b1', 'X', 'A', 'B1')
                    PINPT_BERRY%z2_axis(1) = 1
                    if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode'
                    if_main write(6,'( A)')'  Z2_AXIS: B1 direction of the reciprocal lattice vector (B1-B2 plane, xy)'
                  case('2','y','b','b2', 'Y', 'B', 'B2')
                    PINPT_BERRY%z2_axis(1) = 2
                    if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode'
                    if_main write(6,'( A)')'  Z2_AXIS: B2 direction of the reciprocal lattice vector (B2-B3 plane, yz)'
                  case('3','z','c','b3', 'Z', 'C', 'B3')
                    PINPT_BERRY%z2_axis(1) = 3
                    if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode'
                    if_main write(6,'( A)')'  Z2_AXIS: B3 direction of the reciprocal lattice vector (B3-B1 plane, zx)'
                endselect
              endif

            elseif(PINPT_BERRY%z2_dimension .eq. 1) then
              deallocate(PINPT_BERRY%z2_axis)
              allocate(PINPT_BERRY%z2_axis(1))
              i_dummy = index(dummy, ':')
              if(i_dummy .eq. 0) then    
                PINPT_BERRY%z2_axis(1) = 1
                if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy), ' mode (Only Berry (Zak) phase will be calculated)'
                if_main write(6,'( A)')'  Z2_AXIS: B1 direction of the reciprocal lattice vector (default, parallel to B1)'
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','kx', 'b1', 'X', 'A', 'B1', 'KX')
                    PINPT_BERRY%z2_axis(1) = 1
                    if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'
                    if_main write(6,'( A)')'  Z2_AXIS: B1 direction of the reciprocal lattice vector (parallel to B1)'
                  case('2','y','b','ky', 'b2', 'Y', 'B', 'B2', 'KY')
                    PINPT_BERRY%z2_axis(1) = 2
                    if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'
                    if_main write(6,'( A)')'  Z2_AXIS: B2 direction of the reciprocal lattice vector (parallel to B2)'
                  case('3','z','c','kz', 'b3', 'Z', 'C', 'B3', 'KZ')
                    PINPT_BERRY%z2_axis(1) = 3
                    if_main write(6,'(3A)')'   Z2_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'
                    if_main write(6,'( A)')'  Z2_AXIS: B3 direction of the reciprocal lattice vector (parallel to B3)'
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
              if_main write(6,'(3A)')'GET_CHERN: (Z2 SET) .TRUE. (evaluate Chern number for the bands with ERANGE)'
            else
              if_main write(6,'(3A)')'GET_CHERN: (Z2 SET) .FALSE.'
            endif

          case('SET_PHASE', 'GET_PHASE', 'PHASE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_z2_phase

        end select

      enddo set_wann

      if(PINPT_BERRY%z2_dimension .eq. 3) then
        if(PINPT_BERRY%z2_nkpath .eq. 0) then
          if_main write(6,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv,' (division along k-line for wcc))'
          if_main write(6,'(A)')' WARN: The direction of WCC evolution did not provided,'
          if_main write(6,'(A)')'       so we will setup with default =>  Z2_NKPATH = 21'
          PINPT_BERRY%z2_nkpath = 21
        elseif(PINPT_BERRY%z2_nkpath .gt. 0) then
          if_main write(6,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv, ' (division along k-line for wcc))'
          if_main write(6,'(A, I4,A)')' Z2_NKPATH:',PINPT_BERRY%z2_nkpath,' (division along k-direction for wcc evolution))'
        endif
        PINPT_BERRY%z2_nplane = 2
      elseif(PINPT_BERRY%z2_dimension .eq. 2) then
        if(PINPT_BERRY%z2_nkpath .eq. 0) then
          if_main write(6,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv,' (division along k-line for wcc))'
          if_main write(6,'(A)')' WARN: The direction of WCC evolution did not provided,'
          if_main write(6,'(A)')'       so we will setup with default =>  Z2_NKPATH = 21'
          PINPT_BERRY%z2_nkpath = 21
        elseif(PINPT_BERRY%z2_nkpath .gt. 0) then
          if_main write(6,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv, ' (division along k-line for wcc))'
          if_main write(6,'(A, I4,A)')' Z2_NKPATH:',PINPT_BERRY%z2_nkpath,' (division along k-direction for wcc evolution))'
        endif
        PINPT_BERRY%z2_nplane = 1
      elseif(PINPT_BERRY%z2_dimension .eq. 1) then
        if(PINPT_BERRY%z2_nkpath .eq. 0) then
          if_main write(6,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv,' (division along k-line for wcc))'
          PINPT_BERRY%z2_nkpath = 1
        elseif(PINPT_BERRY%z2_nkpath .gt. 0) then
          if_main write(6,'(A, I4,A)')' Z2_NKDIV: ',PINPT_BERRY%z2_nkdiv, ' (division along k-line for wcc))'
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

      if_main write(6,'(A)')'  GET_WCC: .TRUE. (GET Wannier charge center)'
  set_wann: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('WCC_ERANGE')  ! get wcc_range
            call strip_off (trim(inputline), PINPT_BERRY%strip_wcc_range, 'WCC_ERANGE', ' ' , 2)
            if_main write(6,'(3A)')' WCC_RANG:  ',adjustl(trim(PINPT_BERRY%strip_wcc_range)), &
                                           ' (ENERGY RANGE TO BE CALCULATED)'
          case('WCC_FNAME') ! WCC output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%wcc_filenm = trim(desc_str)
            if_main write(6,'(A,A)')' WCC_FNAM:  ',trim(PINPT_BERRY%wcc_filenm)

          case('WCC_GAP_FNAME','WCC_FNAME_GAP') ! WCC largest gap output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%wcc_gap_filenm = trim(desc_str)
            if_main write(6,'(A,A)')' WCC_GAPF:  (largest gap file)',trim(PINPT_BERRY%wcc_gap_filenm)

          ! this option is only for test purpose...
          case('WCC_PATH_SHIFT','WCC_KPATH_SHIFT')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%wcc_kpath_shift(1:3)

          case('WCC_PATH','WCC_KPATH')
            ikpath = ikpath + 1
            read(inputline,*,iostat=i_continue) desc_str, wcc_kpath_dummy(1:3,1,ikpath), wcc_kpath_dummy(1:3,2,ikpath)
            if_main write(6,'(A,3F12.8,A,3F12.8,A)')' WCC_PATH: (KPATH) ', &
                                                              wcc_kpath_dummy(1:3,1,ikpath),' --> ', &
                                                              wcc_kpath_dummy(1:3,2,ikpath),' (closed loop for WCC, reci unit)'

          case('WCC_DIREC','WCC_DIRECT','WCC_DIRECTION')
            PINPT_BERRY%flag_wcc_evolve = .true.
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%wcc_direction
            if_main write(6,'(3A)')' WCC_DIR :  ',k_direct(PINPT_BERRY%wcc_direction),' (k-direction of WCC evolution)'

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
              if_main write(6,'(3A)')'GET_CHERN: (WCC SET) .TRUE. (Chern number for the bands with ERANGE)'
            else
              if_main write(6,'(3A)')'GET_CHERN: (WCC SET) .FALSE.'
            endif

          case('GET_CHERN_SPIN','WCC_CHERN_SPIN', 'GET_SPIN_CHERN', 'WCC_SPIN_CHERN')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_wcc_get_chern_spin
            if(PINPT_BERRY%flag_wcc_get_chern_spin) then
              if_main write(6,'(3A)')'GET_CHERN:  SPIN_CHERN .TRUE. (WCC SET) (spin Chern number for the bands with ERANGE)'
            else
              if_main write(6,'(3A)')'GET_CHERN:  SPIN_CHERN .FALSE. (WCC SET)'
            endif

          case('SET_PHASE', 'GET_PHASE', 'PHASE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_wcc_phase

        end select

      enddo set_wann

      if(PINPT_BERRY%flag_wcc_evolve .and. PINPT_BERRY%wcc_nkdiv2 .eq. 1) then
        PINPT_BERRY%wcc_nkdiv2 = 31 ! set default if wcc_evolution along WCC_DIRECT is requested
      endif

      if(PINPT_BERRY%wcc_nkdiv2 .gt. 1) then
        if_main write(6,'(A,2I4,3A)')' WCC_KDIV: ',PINPT_BERRY%wcc_nkdiv, &
                                                          PINPT_BERRY%wcc_nkdiv2, &
                                            ' (division along k-line(wcc) and k-direct(evolution along ', &
                                            k_direct(PINPT_BERRY%wcc_direction),')'
      elseif(PINPT_BERRY%wcc_nkdiv2 .eq. 1) then
        if_main write(6,'(A, I4,A)')' WCC_KDIV: ',PINPT_BERRY%wcc_nkdiv, &
                                            ' (division along k-path(wcc))'
        if_main write(6,'(A)')' WARN: The direction of WCC evolution did not provided,'
        if_main write(6,'(A)')'       so we will setup with default => WCC_DIRECT = 2 ,'
        if_main write(6,'(A)')'                                        WCC_KDIV2  = 31'
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

      if_main write(6,'(A)')' ZAK_PHAS: .TRUE. (GET ZAK PHASE)'
   set_zak: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('ZAK_ERANGE')  ! get zak_range
            call strip_off (trim(inputline), PINPT_BERRY%strip_zak_range, 'ZAK_ERANGE', ' ' , 2) 
            if_main write(6,'(3A)')' ZAK_RANG:  ',adjustl(trim(PINPT_BERRY%strip_zak_range)), &
                                           ' (ENERGY RANGE TO BE CALCULATED)'

          case('ZAK_FNAME') ! ZAK phase output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%zak_filenm = trim(desc_str)
            if_main write(6,'(A,A)')' ZAK_FNAM:  ',trim(PINPT_BERRY%zak_filenm)

          case('ZAK_SEPARATE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT%flag_zak_separate
            if(PINPT%flag_zak_separate) then
              if_main write(6,'(A)')' ZAK_SEPR: .TRUE. (also separate into each eigenstate)'
            elseif(.not. PINPT%flag_zak_separate) then
              if_main write(6,'(A)')' ZAK_SEPR: .FALSE. (do not separate into each eigenstate)'
            endif

          ! this option is only for test purpose...
          case('ZAK_PATH_SHIFT','ZAK_KPATH_SHIFT')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%zak_kpath_shift(1:3)

          case('ZAK_PATH','ZAK_KPATH')
            ikpath = ikpath + 1
            read(inputline,*,iostat=i_continue) desc_str, zak_kpath_dummy(1:3,1,ikpath), zak_kpath_dummy(1:3,2,ikpath)
            if_main write(6,'(A,3F12.8,A,3F12.8,A)')' ZAK_PATH: (KPATH) ', &
                                                              zak_kpath_dummy(1:3,1,ikpath),' --> ', &
                                                              zak_kpath_dummy(1:3,2,ikpath),' (closed loop for Zak phase, reci unit)'
          case('ZAK_DIREC','ZAK_DIRECT','ZAK_DIRECTION')
            PINPT_BERRY%flag_zak_evolve = .true.
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%zak_direction
            if_main write(6,'(3A)')' ZAK_DIR :  ',k_direct(PINPT_BERRY%zak_direction),' (k-direction of Zak phase evolution)'

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
        if_main write(6,'(A,2I4,3A)')' ZAK_KDIV: ',PINPT_BERRY%zak_nkdiv, &
                                                  PINPT_BERRY%zak_nkdiv2, &
                                    ' (division along k-path(zak_phase) and k-direct(evolution along ', &
                                    k_direct(PINPT_BERRY%zak_direction),')'
      elseif(PINPT_BERRY%zak_nkdiv2 .eq. 1) then
        if_main write(6,'(A, I4,A)')' ZAK_KDIV: ',PINPT_BERRY%zak_nkdiv, &
                                    ' (division along k-path(zak_phase))'
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
         if_main write(6,'(A,3F12.8,A,3F12.8)')' ZAK_PATH:  (KPATH,default) ',PINPT_BERRY%zak_kpath(1:3,1,1), &
                                                                             '  --> ',PINPT_BERRY%zak_kpath(1:3,2,1)
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

      if_main write(6,'(A)')' '
      if_main write(6,'(A)')' BERRYCRV: .TRUE. (GET BERRY CURVATURE)'

 set_berryc: do while(trim(desc_str) .ne. 'END')
        read(pid_incar,'(A)',iostat=i_continue) inputline
        read(inputline,*,iostat=i_continue) desc_str  ! check INPUT tag
        if(i_continue .ne. 0) cycle      ! skip empty line
        if(desc_str(1:1).eq.'#') cycle   ! skip comment
        if(trim(desc_str).eq.'END') exit ! exit loop if 'END'

        select case ( trim(desc_str) )
          case('BERRYC_ERANGE')
            call strip_off (trim(inputline), PINPT_BERRY%strip_bc_range, 'BERRYC_ERANGE', ' ' , 2) ! get berrycurv_range
            if_main write(6,'(3A)')'   ERANGE:  ',adjustl(trim(PINPT_BERRY%strip_bc_range)),' (ENERGY RANGE TO BE CALCULATED)'

          case('BERRYC_FNM','BERRYC_FNAME') ! BERRY CURVATURE output file name
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            PINPT_BERRY%bc_filenm = trim(desc_str)
            if_main write(6,'(3A)')'   FONAME:  ',trim(PINPT_BERRY%bc_filenm),' (OUTPUT FILE NAME HEADER for BERRY CURVATURE)'
            PINPT_BERRY%flag_bc_filenm_provided = .true.

          case('BERRYC_METHOD') ! BERRY CURVATURE evalulation method
            read(inputline,*,iostat=i_continue) desc_str,desc_str
            if(trim(desc_str) .eq. 'KUBO') then
              PINPT_BERRY%flag_bc_method_kubo  = .true.
              PINPT_BERRY%flag_bc_method_resta = .false.
              PINPT%flag_get_orbital = .true.
              if_main write(6,'(A)')'   METHOD:  Berry curvature (Omega) is obtained via Kubo formula, i.e., '
              if_main write(6,'(A)')'             Omega_z(k) = -2Im*SIGMA_{m/=m} v_nm,x * v_mn,y / (e_m - e_n)^2'
              if_main write(6,'(A)')'                 v_nm,i = <u_nk|dH(k)/dk_i|u_mk> .'
              if_main write(6,'(A)')'             See Ref. [Wang et. al., PRB 74, 195118 (2006)]'
            elseif(trim(desc_str) .eq. 'RESTA' .or. trim(desc_str) .eq. 'FUKUI') then
              PINPT_BERRY%flag_bc_method_kubo = .false.
              PINPT_BERRY%flag_bc_method_resta = .true.
              PINPT%flag_get_orbital = .true.
              if_main write(6,'(A)')'   METHOD:  Berry curvature (Omega) is obtained via descritization method,'
              if_main write(6,'(A)')'            where PI_i runs over closed loop of k normal to kx,ky plane, i.e.,'
              if_main write(6,'(A)')'             Omega_z(k) = -Im ln PI_i det M(k_i, k_i+1) / ka'
              if_main write(6,'(A)')'             Here, M is overlap matrix and the matrix element M_mn is as follow'
              if_main write(6,'(A)')'             M_mn(k_i, k_i+1) = <u_m(k)|u_n(k+1)>, and ka is area of the closed loop'
              if_main write(6,'(A)')'            See Ref. [Resta, J. Phys.: Condens. Matter 12, R107 (2000)]'
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
              if_main write(6,'(A)')' !WARN!: BERRYC_DIMENSION should be start with one of following strings:'
              if_main write(6,'(A)')'  3D  => BERRYC_DIMENSION 3D'
              if_main write(6,'(A)')'  2D  => BERRYC_DIMENSION 2D:z # z is normal direction of surface'
              if_main write(6,'(A)')'  1D  => BERRYC_DIMENSION 1D:x # x is parallel direction of 1D system'
              if_main write(6,'(A)')'  For 2D and 1D, "x", "y", or "z" indicates the normal and parallel'
              if_main write(6,'(A)')'  direction of the system'
              if_main write(6,'(A)')'  Exit program anyway... Check your "SET BERRYCURV" in "INCAR-TB" again.'
              stop
            endif

            call str2int(dummy(1:i_dummy-1),PINPT_BERRY%bc_dimension)

            if(PINPT_BERRY%bc_dimension .eq. 3) then
              deallocate(PINPT_BERRY%bc_axis)
              allocate(PINPT_BERRY%bc_axis(3))
              PINPT_BERRY%bc_axis(:) = (/1,2,3/)
              if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy), ' mode (three direction will be checked)'

            elseif(PINPT_BERRY%bc_dimension .eq. 2) then
              deallocate(PINPT_BERRY%bc_axis)
              allocate(PINPT_BERRY%bc_axis(1))
              i_dummy = index(dummy,':')
              if(i_dummy .eq. 0) then
                PINPT_BERRY%bc_axis(1) = 3 ! default : third reciprocal axis will be chosen for the WCC surface normal axis
                if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy), ' mode (three direction will be checked)'
                if_main write(6,'( A)')'  BC_AXIS: B3 direction of the reciprocal lattice vector (default)'
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','b1', 'X', 'A', 'B1')
                    PINPT_BERRY%bc_axis(1) = 1
                    if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode'
                    if_main write(6,'( A)')'  BC_AXIS: B1 direction of the reciprocal lattice vector (B2-B3 plane, yz)'
                  case('2','y','b','b2', 'Y', 'B', 'B2')
                    PINPT_BERRY%bc_axis(1) = 2
                    if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode'
                    if_main write(6,'( A)')'  BC_AXIS: B2 direction of the reciprocal lattice vector (B3-B1 plane, zx)'
                  case('3','z','c','b3', 'Z', 'C', 'B3')
                    PINPT_BERRY%bc_axis(1) = 3
                    if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode'
                    if_main write(6,'( A)')'  BC_AXIS: B3 direction of the reciprocal lattice vector (B1-B2 plane, xy)'
                endselect
              endif

            elseif(PINPT_BERRY%bc_dimension .eq. 1) then
              deallocate(PINPT_BERRY%bc_axis)
              allocate(PINPT_BERRY%bc_axis(1))
              i_dummy = index(dummy, ':')
              if(i_dummy .eq. 0) then
                PINPT_BERRY%bc_axis(1) = 1
                if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy), ' mode (Only Berry (Zak) phase will be calculated)'
                if_main write(6,'( A)')'  BC_AXIS: B1 direction of the reciprocal lattice vector (default, parallel to B1)'
              elseif(i_dummy .ne. 0) then
                call strip_off (trim(dummy), dummy2, ':', ' ' , 2)
                select case (trim(dummy2))
                  case('1','x','a','kx', 'b1', 'X', 'A', 'B1', 'KX')
                    PINPT_BERRY%bc_axis(1) = 1
                    if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'
                    if_main write(6,'( A)')'  BC_AXIS: B1 direction of the reciprocal lattice vector (parallel to B1)'
                  case('2','y','b','ky', 'b2', 'Y', 'B', 'B2', 'KY')
                    PINPT_BERRY%bc_axis(1) = 2
                    if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'
                    if_main write(6,'( A)')'  BC_AXIS: B2 direction of the reciprocal lattice vector (parallel to B2)'
                  case('3','z','c','kz', 'b3', 'Z', 'C', 'B3', 'KZ')
                    PINPT_BERRY%bc_axis(1) = 3
                    if_main write(6,'(3A)')'   BC_DIM: ', trim(dummy(1:i_dummy-1)),' mode (Only Berry (Zak) phase will be calculated)'
                    if_main write(6,'( A)')'  BC_AXIS: B3 direction of the reciprocal lattice vector (parallel to B3)'
                endselect

              endif
            endif

          case('SET_PHASE', 'GET_PHASE', 'PHASE')
            read(inputline,*,iostat=i_continue) desc_str,PINPT_BERRY%flag_bc_phase

        end select
      enddo set_berryc

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
              if_main write(6,'(A,3F10.5)')'PARITY_KP: ', kp(:,nkp), kp_name(nkp)
            elseif(i_dummy .gt. 4) then
              read(inputline,*,iostat=i_continue) desc_str, kp(1:3, nkp), kp_name(nkp)
              if(.not. flag_number(kp_name(nkp))) then
                if_main write(6,'(A,3F10.5,2x,A)')'PARITY_KP: ', kp(:,nkp), trim(kp_name(nkp))
              else
                if_main write(6,'(4A)')'    !WANR! Check PARITY_CHECK SETTING tags ->',trim(desc_str), ' ', trim(func)
                stop
              endif
            endif
          
          case('PARITY_ORIGIN', 'ORIGIN', 'ORIGIN_SHIFT')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_origin(1:3)
            if_main write(6,'(A,3F10.5,A)')'   ORIGIN: ', PINPT_BERRY%parity_origin(1:3), ' (used in PARITY_CHECK)'

          case('PARITY_OP1','SYMMETRY_OP1','ROTATION_MAT1','ROTATION1')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_operator(1:3,1)
          case('PARITY_OP2','SYMMETRY_OP2','ROTATION_MAT2','ROTATION2')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_operator(1:3,2)
          case('PARITY_OP3','SYMMETRY_OP3','ROTATION_MAT3','ROTATION3')
            read(inputline,*,iostat=i_continue) desc_str, PINPT_BERRY%parity_operator(1:3,3)

        end select
      enddo set_parity

      if( nkp .eq. 0) then
        PINPT_BERRY%parity_nkpoint = 1
        allocate(PINPT_BERRY%parity_kpoint(3,1))
        allocate(PINPT_BERRY%parity_kpoint_reci(3,1))
        allocate(PINPT_BERRY%parity_kpoint_name(1))
        if_main write(6,'(A,3F10.5,2x,A)')'PARITY_KP: ', (/0d0, 0d0, 0d0/),'Gamma'
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
   
      if_main write(6,'(A,3(F8.5,2x),A)')'  PARITY  [  ', PINPT_BERRY%parity_operator(:,1),']                          '
      if_main write(6,'(A,3(F8.5,2x),A)')' OPERATOR=[  ', PINPT_BERRY%parity_operator(:,2),'] => R(inv) = P * ( R - ORIGIN )'
      if_main write(6,'(A,3(F8.5,2x),A)')'     P    [  ', PINPT_BERRY%parity_operator(:,3),']                          '

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
        if_main write(6,'(A)')'  LNONCOL: .TRUE. (non-collinear = .TRUE.)'
        if_main write(6,'(A)')'    ISPIN: 2 (noncollinear = .TRUE.)'
        if_main write(6,'(A)')'  ISPINOR: 2 (noncollinear = .TRUE.)'
        PINPT%ispin = 2
        PINPT%ispinor = 2
        PINPT%nspin = 1
      elseif(desc_str(1:6) .eq. 'collin' .or. desc_str(1:3) .eq. 'col' .or. desc_str(1:2) .eq. 'cl') then
        PINPT%flag_collinear    = .true.
        PINPT%flag_noncollinear = .false.
        if_main write(6,'(A)')'    ISPIN: 2 (collinear = .TRUE.)'
        if_main write(6,'(A)')'  ISPINOR: 1 (collinear = .TRUE.)'
        PINPT%ispin = 2
        PINPT%ispinor = 1
        PINPT%nspin =2
      elseif(desc_str(1:6) .eq. 'nonmag' .or. desc_str(1:2) .eq. 'nm' .or. desc_str(1:7) .eq. 'non-mag') then
        PINPT%flag_collinear    = .false.
        PINPT%flag_noncollinear = .false.
        if_main write(6,'(A)')'    ISPIN: 1 (non-magnetic = .true.)'
        if_main write(6,'(A)')'  ISPINOR: 1 (non-magnetic = .true.)'
        PINPT%ispin = 1
        PINPT%ispinor = 1
        PINPT%nspin =1
      else
        if_main write(6,'(A)')'  !WARNING! TYPMAG is not properly set. '
        if_main write(6,'(A)')'  !WARNING! usage: for non-magnetic run : "NM", "NONMAGNETIC", "NON-MAGNETIC" '
        if_main write(6,'(A)')'  !WARNING!        for collinear    run : "COLLINEAR", "COL", "CL" '
        if_main write(6,'(A)')'  !WARNING!        for noncollinear run : "NONCOLLINEAR", "NON-COLLINEAR", "NC"'
        if_main write(6,'(A)')'  !WARNING!        are acceptable (case insensitive). Other option is not allowed. Exit program.'
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
        if_main write(6,'(A)')'    LSORB: .TRUE. (non-collinear = .TRUE. .AND. spin-orbit = .TRUE.)'
        PINPT%ispin = 2
        PINPT%ispinor = 2
        PINPT%flag_noncollinear = .true.
      elseif(.not.PINPT%flag_soc)then
        if_main write(6,'(A)')'    LSORB: .FALSE. ( spin-orbit = .FALSE.)'
      elseif(PINPT%flag_soc .and. PINPT%flag_collinear) then
        if_main write(6,'(A)')'  !WARNING! LSORB and TYPMAG is not properly set. '
        if_main write(6,'(A)')'  !WARNING! COLLINEAR = .true. and LSORB = .true. which is not proper.'
        if_main write(6,'(A)')'  !WARNING! Both option are not available simultaneously. Exit program.'
        stop
      endif

      return
   endsubroutine

   subroutine set_tbparam_file(PINPT, param_const, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      real*8        param_const(5,max_nparam)
      character(*), parameter :: func = 'set_tbparam_file'

      if(.not. PINPT%flag_pfile) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%pfilenm
        if(allocated(PINPT%param) .or. allocated(PINPT%param_name) .or. PINPT%flag_pincar) then
           PINPT%flag_pincar=.false.
           deallocate(PINPT%param)
           deallocate(PINPT%param_name)
           if_main write(6,'(A)')'  !WARN!  TB-parameter is alread set by TBPARAM in the INCAR-TB,'
           if_main write(6,'(A)')'  !WARN!  however, since the external TB-parameter file is provided,'
           if_main write(6,'(A,A,A)')'  !WARN!  those valeus from ',trim(PINPT%pfilenm),' will be read in priori'
        endif
      endif
      if_main write(6,'(A,A)')' PARA_FNM:  ',trim(PINPT%pfilenm)
      call read_param(PINPT, param_const)

      return
   endsubroutine

   subroutine set_kpoint_file(PINPT, flag_kfile_ribbon, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_kpoint_file'
      logical       flag_kfile_ribbon 

      read(inputline,*,iostat=i_continue) desc_str, PINPT%kfilenm
      if(flag_kfile_ribbon) then
        PINPT%kfilenm = PINPT%ribbon_kfilenm
        if_main write(6,'(A,A)')' KPTS_FNM:  (for ribbon) ',trim(PINPT%kfilenm)
      else
        read(inputline,*,iostat=i_continue) desc_str, PINPT%kfilenm
        if_main write(6,'(A,A)')' KPTS_FNM:  ',trim(PINPT%kfilenm)
      endif

      return
   endsubroutine

   subroutine set_geom_file(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_geom_file'

      read(inputline,*,iostat=i_continue) desc_str, PINPT%gfilenm
      if_main write(6,'(A,A)')' GEOM_FNM:  ',trim(PINPT%gfilenm)

      return
   endsubroutine

   subroutine set_tbparam_out_file(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_tbparam_out_file'

      read(inputline,*,iostat=i_continue) desc_str, PINPT%pfileoutnm
      if_main write(6,'(A,A40)')' POUT_FNM:  ',PINPT%pfileoutnm
      PINPT%flag_print_param = .true.

      return
   endsubroutine

   subroutine set_hopping_type(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_hopping_type'

      read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_slater_koster
      if(PINPT%flag_slater_koster) then
        if_main write(6,'(A,A)')' TYPE_HOP:  ','SLATER_KOSTER'
      elseif(.not. PINPT%flag_slater_koster) then
        if_main write(6,'(A,A)')' TYPE_HOP:  ','USER_DEFINED'
      endif

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
        if_main write(6,'(A)')'   LOCCHG: .TRUE. (read local charge density from GFILE)'
      elseif(.not. PINPT%flag_local_charge) then
        if_main write(6,'(A)')'   LOCCHG: .FALSE.'
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
        if_main write(6,'(A)')'   PLUS+U: .TRUE. '
      elseif(.not. PINPT%flag_plus_U) then
        if_main write(6,'(A)')'   LOCCHG: .FALSE.'
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
            if_main write(6,'(5A)')'  !WARN! You requested a HOPPING file via "', trim(tag_name), &
                                   '" tag, but the file "', trim(PINPT%nnfilenm),'" does not exist! Exit...'
            stop
          endif

        elseif(i_dummy .eq. 1) then
          if_main write(6,'(4A)')'  !WARN! You requested a HOPPING file via "', trim(tag_name), &
                                   '" tag, but the file name has not been defined. Exit...'
          stop
        endif
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
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmu, PINPT%read_energy_column_index
          endif
          if_main write(6,'(A,A)')' EDFT_FNM: ',trim(PINPT%efilenmu)
          flag_read_energy=.true.

        case('EFILEU')
          i_dummy = nitems(inputline) -1
          if(i_dummy .eq. 1) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmu
          elseif(i_dummy .eq. 2) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmu, PINPT%read_energy_column_index
          endif
          if_main write(6,'(A,A)')' EDFT_FNM: spin-up= ',trim(PINPT%efilenmu)
          flag_read_energy=.true.

        case('EFILED')
          i_dummy = nitems(inputline) -1
          if(i_dummy .eq. 1) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmd
          elseif(i_dummy .eq. 2) then
            read(inputline,*,iostat=i_continue) desc_str, PINPT%efilenmd, PINPT%read_energy_column_index_dn
          endif
          if_main write(6,'(A,A)')' EDFT_FNM: spin-dn= ',trim(PINPT%efilenmd)
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
          if_main write(6,'(A,I8)')' INI_BAND:',PWGHT%iband

         case('FBAND')
          read(inputline,*,iostat=i_continue) desc_str, PWGHT%fband
          if_main write(6,'(A,I8)')' FIN_BAND:',PWGHT%fband

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
      if_main write(6,'(A,F8.2,A,I5,A)')'  SCISSOR: .TRUE. => EDFT(n,k) + ',PINPT%r_scissor,' (if n >=',PINPT%i_scissor,')'

      return
   endsubroutine

   subroutine set_local_orbital_print(PINPT, inputline)
      type(incar)  ::  PINPT
      integer*4     i_continue
      integer*4     nitems
      integer*4     i_dummy
      character*132 inputline
      character*2   dummy
      character*40  desc_str
      character*2   str2lowcase
      character(*), parameter :: func = 'set_local_orbital_plot'
      external      nitems, str2lowcase

      i_dummy = nitems(inputline) - 1
      if(i_dummy .eq. 1) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_orbital
        if(PINPT%flag_print_orbital) then
          if_main write(6,'(A)')'  L_ORBIT: .TRUE. | print out projected orbital weight'
          PINPT%flag_get_orbital = .true.
        elseif( .not. PINPT%flag_print_orbital) then
          if_main write(6,'(A)')'  L_ORBIT: .FALSE.'
          PINPT%flag_get_orbital = .false.
        endif
    
      elseif(i_dummy .eq. 2) then
        PINPT%flag_print_mag = .TRUE.
        read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_print_orbital, PINPT%axis_print_mag
        PINPT%axis_print_mag=str2lowcase(PINPT%axis_print_mag)
        if(PINPT%flag_print_orbital) then
          if(PINPT%axis_print_mag .ne. 're' .or. PINPT%axis_print_mag .ne. 'im') then
            if_main write(6,'(2A)')'  L_ORBIT: .TRUE. | print out magnetization <sigma>: ', PINPT%axis_print_mag
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 're') then
            if_main write(6,'(2A)')'  L_ORBIT: .TRUE. | print out real part of wavefnc.: ', PINPT%axis_print_mag
            PINPT%flag_get_orbital = .true.
          elseif(PINPT%axis_print_mag .eq. 'im') then
            if_main write(6,'(2A)')'  L_ORBIT: .TRUE. | print out imag part of wavefnc.: ', PINPT%axis_print_mag
          endif
        elseif( .not. PINPT%flag_print_orbital) then
          if_main write(6,'(A)')'  L_ORBIT: .FALSE.'
          PINPT%flag_get_orbital = .false.
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
        if_main write(6,'(A)')' KPT_UNIT: KPOINT UNIT = Reciprocal'
      elseif(PKPTS%kunit .eq. 'a' .or. PKPTS%kunit .eq. 'A') then
        PKPTS%kunit = 'A'
        if_main write(6,'(A)')' KPT_UNIT: KPOINT UNIT = Angstrom unit (1/A)'
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
             if_main write(6,'(A)')'  L_TBFIT: .TRUE.'
          elseif(.not. PINPT%flag_tbfit) then
             if_main write(6,'(A)')'  L_TBFIT: .FALSE.'
          else
             if_main write(6,'(A)')'  L_TBFIT:  unknown input variable.. exit program'
             stop
          endif

        case('LSTYPE') !set non-linear regression scheme
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ls_type
          if(PINPT%ls_type .eq. 'LMDIF' .or. PINPT%ls_type .eq. 'lmdif' ) then
            if_main write(6,'(A,A6,A)')'  LS_TYPE:  ',PINPT%ls_type,', Levenberg-Marquardt with finite-difference for Jacobian'
          elseif(PINPT%ls_type .eq. 'PIKAIA' .or. PINPT%ls_type .eq. 'GA' .or. &
                 PINPT%ls_type .eq. 'pikaia' .or. PINPT%ls_type .eq. 'ga' ) then
            if_main write(6,'(A,A6,A)')'  LS_TYPE:  ',PINPT%ls_type,', Genetic algorithm based on PIKAIA library'
          else
            if_main write(6,'(A,A6,A)')'  LS_TYPE:  ',PINPT%ls_type,' is not defined or not available in the current version. Exit..'
            stop
          endif

        case('PTOL') !set parameter tolerance for iteration step
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ptol
          if_main write(6,'(A,F15.8)')'    P_TOL:  ',PINPT%ptol

        case('FTOL') !set function tolerance for iteration step
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ftol
          if_main write(6,'(A,F15.8)')'    F_TOL:  ',PINPT%ftol

        case('MITER') !set maximum iteration % maximum # of generations for GA
          read(inputline,*,iostat=i_continue) desc_str, PINPT%miter
          if_main write(6,'(A,I8)')' MAX_ITER:  ',PINPT%miter
          
        case('NPOP') ! set number of generation for genetic algorithm
          read(inputline,*,iostat=i_continue) desc_str, PINPT%ga_npop
          if_main write(6,'(A,I8)')'  GA_NPOP:  ',PINPT%ga_npop

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
      if_main write(6,'(A,L)')' P_TARGET:',PINPT%flag_print_only_target

      return
   endsubroutine

   subroutine set_onsite_tol(NN_TABLE, inputline)
      type(hopping)  ::  NN_TABLE
      integer*4     i_continue
      character*132 inputline
      character*40  desc_str
      character(*), parameter :: func = 'set_onsite_tol'

      read(inputline,*,iostat=i_continue) desc_str, desc_str, NN_TABLE%onsite_tolerance
      if_main write(6,'(A,F16.4)')'ONSITETOL:',NN_TABLE%onsite_tolerance

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
     

      if_main write(6,'(A,3F10.5,A)')'   * PARAMETERS FOR GENETIC ALGORITHM *'
      if_main write(6,'(A,3F10.5,A)')'   -----------------------------------'

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

          case('NPOP')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%npop
            if_main write(6,'(A,I8      )')'     NPOP: ', PKAIA%npop

          case('NGENE')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%ngene
            if_main write(6,'(A,I8      )')'    NGENE: ', PKAIA%ngene

          case('PCROSS')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pcross
            if_main write(6,'(A, F10.5  )')'   PCROSS: ', PKAIA%pcross

          case('RMUTMIN')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pmutmn
            if_main write(6,'(A, F10.5  )')' MIN_MUTR: ', PKAIA%pmutmn

          case('RMUTMAX')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pmutmx
            if_main write(6,'(A, F10.5  )')' MAX_MUTR: ', PKAIA%pmutmx

          case('RMUTINI')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%pmut
            if_main write(6,'(A, F10.5  )')' INI_MUTR: ', PKAIA%pmut

          case('MUT_MOD')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%imut
            if_main write(6,'(A, I8     )')' MUT_MODE: ', PKAIA%imut

          case('FDIF')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%fdif
            if_main write(6,'(A, F10.5  )')'     FDIF: ', PKAIA%fdif

          case('IREP')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%irep
            if_main write(6,'(A, I8     )')'     IREP: ', PKAIA%irep

          case('IELITE')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%ielite
            if_main write(6,'(A, I8     )')'   IELITE: ', PKAIA%ielite

          case('VERBOSE')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%ivrb  
            if_main write(6,'(A, I8     )')'  VERBOSE: ', PKAIA%ivrb  

          case('CONVTOL')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%convtol
            if_main write(6,'(A, F10.5  )')'  CONVTOL: ', PKAIA%convtol

          case('CONVWIN')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%convwin
            if_main write(6,'(A, F10.5  )')'  CONVWIN: ', PKAIA%convwin

          case('IGUESSF')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%iguessf
            if_main write(6,'(A, F10.5  )')'  IGUESSF: ', PKAIA%iguessf

          case('ISEED'  )
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%iseed  
            if_main write(6,'(A, F10.5  )')' RANDSEED: ', PKAIA%iseed  

          case('UBOUND')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%upper_bound
            if_main write(6,'(A, F10.5  )')'   UBOUND: ', PKAIA%upper_bound

          case('LBOUND')
            read(inputline,*,iostat=i_continue) desc_str, PKAIA%lower_bound
            if_main write(6,'(A, F10.5  )')'   LBOUND: ', PKAIA%lower_bound

        end select
      enddo set_pkaia

      if_main write(6,'(A,3F10.5,A)')'   * END READING: GENETIC ALGORITHM PARAMETERS*'
      if_main write(6,'(A,3F10.5,A)')'   ---------------------------------------------'

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
        if_main write(6,'(A)') '    !WARN! ERANGE tag and EWINDOW tag cannot be used simultaneously.'
        if_main write(6,'(A)') '           EWINDOW tag sets to use extended eigensolver routines based on '
        if_main write(6,'(A)') '           FEAST algorithm for the sparse matrix eigenvalue problem.'
        if_main write(6,'(A)') '           ERANGE tag, on the other hand, sets to use ZHEEVX function in LAPACK routines '
        if_main write(6,'(A)') '           for the dense matrix eigenvalue problem.'
        if_main write(6,'(A)') '           Please deactivate one of the tag to proceed.'
        if_main write(6,'(A)') '           If you dealing with very large system for the band structure calculation,'
        if_main write(6,'(A)') '           it is recommended to use EWINDOW tag rather than ERANGE'
        if_main write(6,'(A)') '           Note that EWINDOW does not applies to the post processings requested by'
        if_main write(6,'(A)') '           various "SET" tags, for example, Z2_INDEX, WCC, ZAK_PHASE, PARITY_CHECK, etc.'
        if_main write(6,'(A)') '           In these routines, ERANGE tags there in will be used and ZHEEVX function will be'
        if_main write(6,'(A)') '           used instead.'
        if_main write(6,'(A)') '           Exit program anyway...'
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
          if_main write(6,'(A,A)')'  !WARN! No energy window has been defined! Please check EWINDOW tag atain. Exit... ', func
          stop
        elseif( nitems(inputline) -1 .eq. 2) then
          read(inputline,*,iostat=i_continue) desc_str,PINPT%feast_emin, PINPT%feast_emax
          if_main write(6,'(A,F12.5,A,F12.5,A)')' E_WINDOW:  [',PINPT%feast_emin,',',PINPT%feast_emax,'], [emin:emax] (eV)'
          PINPT%flag_sparse = .true. ! call FEAST algorithom
        else
          if_main write(6,'(A,A)')'  !WARN! Please check EWINDOW tag syntax. Exit... Current setting = ',trim(dummy)
          if_main write(6,'(A,A)')'         Only accept following syntax: EWINDOW EMIN:EMAX  or'
          if_main write(6,'(A,A)')'                                       EWINDOW EMIN EMAX '
          if_main write(6,'(A,A)')'         Here, EMIN and EMAX are real*8 values for the eigenvalues to be searched.'
          stop
        endif

      elseif(i_dummy .ge. 1) then
        call strip_off(dummy, dummy1,' ',':',0)
        call strip_off(dummy, dummy2,':',' ',2)
        if(len_trim(dummy1) .eq. 0 .or. len_trim(dummy2) .eq. 0) then
          if_main write(6,'(A,A)')'  !WARN! Please check ERANGE tag syntax. Exit... Current setting = ',trim(dummy)
          if_main write(6,'(A,A)')'         Only accept following syntax: ERANGE INIT:FINA  or'
          if_main write(6,'(A,A)')'                                       ERANGE INIT FINA'
          if_main write(6,'(A,A)')'         Here, INIT and FINA are integer values for the eigenvalue index.'
          stop
        endif
        call str2real(dummy1,PINPT%feast_emin)
        call str2real(dummy2,PINPT%feast_emax)
        if_main write(6,'(A,F12.5,A,F12.5,A)')' E_WINDOW:  [',PINPT%feast_emin,',',PINPT%feast_emax,'], [emin:emax] (eV)'
        PINPT%flag_sparse = .true. ! call FEAST algorithom
      endif

      if(PINPT%feast_nemax .ge. 1) then
        ! feast_nemax is upper bound of states to be stored in your output. 
        ! The states, beyond this criteria will be thrown away.
        if_main write(6,'(A,I0,A)')'   NE_MAX: ', PINPT%feast_nemax, ' (maximum # of states within [emin:emax])'
      elseif(PINPT%feast_nemax .le. 0) then
        if_main write(6,'(A,I0,A)')'   NE_MAX: N_ORBIT * ISPINOR [N_ORBIT:number of orbitals ]'
        if_main write(6,'(A)')     '                             [ISPINOR:2 for LSORB=.TRUE. ]'
        if_main write(6,'(A)')     '                             [       :1 for LSORB=.FALSE.]'
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
        if_main write(6,'(A)') '    !WARN! ERANGE tag and EWINDOW tag cannot be used simultaneously.'
        if_main write(6,'(A)') '           EWINDOW tag sets to use extended eigensolver routines based on '
        if_main write(6,'(A)') '           FEAST algorithm for the sparse matrix eigenvalue problem.'
        if_main write(6,'(A)') '           ERANGE tag, on the other hand, sets to use ZHEEVX function in LAPACK routines '
        if_main write(6,'(A)') '           for the dense matrix eigenvalue problem.'
        if_main write(6,'(A)') '           Please deactivate one of the tag to proceed.'
        if_main write(6,'(A)') '           If you dealing with very large system for the band structure calculation,'
        if_main write(6,'(A)') '           it is recommended to use EWINDOW tag rather than ERANGE'
        if_main write(6,'(A)') '           Note that EWINDOW does not applies to the post processings requested by'
        if_main write(6,'(A)') '           various "SET" tags, for example, Z2_INDEX, WCC, ZAK_PHASE, PARITY_CHECK, etc.'
        if_main write(6,'(A)') '           In these routines, ERANGE tags there in will be used and ZHEEVX function will be'
        if_main write(6,'(A)') '           used instead.'
        if_main write(6,'(A)') '           Exit program anyway...'
        stop
      endif

      call strip_off (inputline, dummy,trim(desc_str),' ',2)
      i_dummy=index(dummy,':')
      if(i_dummy .eq. 0) then
        if( nitems(inputline) -1 .eq. 0) then
          if_main write(6,'(A,A)')'  !WARN! No range has been defined! Please check ERANGE tag atain. Exit... ', func
          stop
        elseif( nitems(inputline) -1 .eq. 1) then
          read(inputline,*,iostat=i_continue) desc_str,PINPT%init_erange
          PINPT%fina_erange = PINPT%init_erange
          if_main write(6,'(A,I6,A,I6,A)')'  E_RANGE:  FROM ',PINPT%init_erange,'-th TO ',PINPT%fina_erange,'-th EIGENSTATES'
          PINPT%flag_erange = .true.
        elseif( nitems(inputline) -1 .eq. 2) then
          read(inputline,*,iostat=i_continue) desc_str,PINPT%init_erange, PINPT%fina_erange
          if_main write(6,'(A,I6,A,I6,A)')'  E_RANGE:  FROM ',PINPT%init_erange,'-th TO ',PINPT%fina_erange,'-th EIGENSTATES'
          PINPT%flag_erange = .true.
        else
          if_main write(6,'(A,A)')'  !WARN! Please check ERANGE tag syntax. Exit... Current setting = ',trim(dummy)
          if_main write(6,'(A,A)')'         Only accept following syntax: ERANGE INIT:FINA  or'
          if_main write(6,'(A,A)')'                                       ERANGE INIT FINA'
          if_main write(6,'(A,A)')'         Here, INIT and FINA are integer values for the eigenvalue index.'
          stop
        endif

      elseif(i_dummy .ge. 1) then
        call strip_off(dummy, dummy1,' ',':',0)
        call strip_off(dummy, dummy2,':',' ',2)
        if(len_trim(dummy1) .eq. 0 .or. len_trim(dummy2) .eq. 0) then
          if_main write(6,'(A,A)')'  !WARN! Please check ERANGE tag syntax. Exit... Current setting = ',trim(dummy)
          if_main write(6,'(A,A)')'         Only accept following syntax: ERANGE INIT:FINA  or'
          if_main write(6,'(A,A)')'                                       ERANGE INIT FINA'
          if_main write(6,'(A,A)')'         Here, INIT and FINA are integer values for the eigenvalue index.'
          stop
        endif
        call str2int(dummy1,PINPT%init_erange)
        call str2int(dummy2,PINPT%fina_erange)
        if_main write(6,'(A,I6,A,I6,A)')'  E_RANGE:  FROM ',PINPT%init_erange,'-th TO ',PINPT%fina_erange,'-th EIGENSTATES'
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
        if_main write(6,'(A)')'  GET_BND: .TRUE.'
      elseif(.not.PINPT%flag_get_band) then
        if_main write(6,'(A)')'  GET_BND: .FALSE.'
      endif
      return
   endsubroutine

   subroutine set_nn_max(PINPT, inputline,desc_str)
      implicit none
      type(incar) :: PINPT
      character*132  inputline
      character*40   desc_str, dummy_, dummy
      integer*4      i_continue
      integer*4      nitems
      external       nitems


      call strip_off (inputline, dummy_, ' ', '#', 0) ! cut off unnecessary comments
      if(index(dummy_, trim(desc_str)) .ge. 1) then
        inputline = dummy_
      endif

      if(nitems(inputline) -1 .eq. 3) then
        read(inputline,*,iostat=i_continue) desc_str, PINPT%nn_max(1:3)
        if_main write(6,'(A,3I5)')'   NN_MAX:', PINPT%nn_max(1:3)
      else
        if_main write(6,'(A)')'    !WARN! Error in reading NN_MAX tag. Exit program...'
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
        if_main write(6,'(A)')'   SPGLIB: .TRUE.'
      elseif(.not.PINPT%flag_spglib) then
        if_main write(6,'(A)')'   SPGLIB: .FALSE.'
      endif

      return
   endsubroutine
#endif

endmodule
