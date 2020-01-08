subroutine print_param (PINPT, iter, title, flag_print_param)
  use parameters, only : incar
  implicit none
  integer*4 i,j,k,iter, pid_param_new
  integer*4 nsub
  character ( len = * ) title
  logical flag_print_param
  type (incar)   :: PINPT

  if(.not. flag_print_param) then
    pid_param_new = 6
  elseif(flag_print_param)then
    pid_param_new = 81
    open(pid_param_new, file = trim(title), status = 'unknown')
  endif

  if( iter .ge. 1 ) then
    
    if(pid_param_new .eq. 6) write ( pid_param_new, '(A,i4,A)',ADVANCE='NO') '   PARAM(iter=',iter,')='
    do i = 1, PINPT%nparam
      if(PINPT%slater_koster_type .gt. 10) then
         k=nint(PINPT%param_const_nrl(1,1,i))
         nsub = PINPT%param_nsub(i)
        if( k .eq. 0) then
          write ( pid_param_new, '(*(F9.3))',ADVANCE='NO' ) PINPT%param_nrl(1:nsub,i)
        elseif( k .ge. 1) then
          write ( pid_param_new, '(*(F9.3))',ADVANCE='NO' ) PINPT%param_nrl(1:nsub,k)
        endif
      else
        if( nint(PINPT%param_const(1,i)) .eq. 0) then
          write ( pid_param_new, '(F9.3)',ADVANCE='NO' ) PINPT%param(i)
        elseif( nint(PINPT%param_const(1,i)) .ge. 1) then
          write ( pid_param_new, '(F9.3)',ADVANCE='NO' ) PINPT%param( nint(PINPT%param_const(1,i)) )
        endif
      endif
    enddo
    write ( *, '(A)' ) ''

  elseif ( iter .le. 0 ) then

    write ( pid_param_new, '(A)' ) ' '
    if(flag_print_param) then
      write ( pid_param_new, '(A,A)' ) '# ',trim ( title )
      if(PINPT%flag_collinear) then
        write(pid_param_new,'(A,A)')'# E_TARGET FILE NAME (up): ',trim(PINPT%efilenmu)
        write(pid_param_new,'(A,A)')'# E_TARGET FILE NAME (dn): ',trim(PINPT%efilenmd)
      else
        write(pid_param_new,'(A,A)')'# E_TARGET FILE NAME : ',trim(PINPT%efilenmu)
      endif
      if(PINPT%flag_scissor) then 
        write(pid_param_new,'(A,F8.2,A,I5,A)')'# SCISSOR: .TRUE. => EDFT(n,k) + ',PINPT%r_scissor,' (if n >=',PINPT%i_scissor,')'
      endif
      if(PINPT%nweight .ge. 1) then
        do i = 1, PINPT%nweight
          write(pid_param_new,'(3(A13,A20),A13,A8)')'# KRANGE ',trim(PINPT%strip_kp(i)), &
                                                    '  TBABND ',trim(PINPT%strip_tb(i)), &
                                                    '  DFTBND ',trim(PINPT%strip_df(i)), &
                                                    '  WEIGHT ',trim(PINPT%strip_wt(i))
        enddo
      endif
      if(PINPT%npenalty_orb .ge. 1) then
        do i = 1, PINPT%npenalty_orb
          write(pid_param_new,'(5(A13,A20))')'# KRANGE ',trim(PINPT%strip_kp_orb(i)), &
                                         '  TBABND ',trim(PINPT%strip_tb_orb(i)), &
                                         '  ORB_RANGE ',trim(PINPT%strip_orb(i)), &
                                         '  SITE_RANGE ',trim(PINPT%strip_site(i)), &
                                         '  PENALTY ',trim(PINPT%strip_pen_orb(i))
        enddo
      endif
      
      if(PINPT%flag_pfile_index) then
        write(pid_param_new, '(A20,A11)') ' PRINT_INDEX  ','.TRUE.'
      else
        write(pid_param_new, '(A20,A11)') ' PRINT_INDEX  ','.FALSE.'
      endif

      if(PINPT%flag_use_overlap) then
        write(pid_param_new, '(A20,A11)') ' USE_OVERLAP  ','.TRUE.'
      else
        write(pid_param_new, '(A20,A11)') ' USE_OVERLAP  ','.FALSE.'
      endif

      if(PINPT%slater_koster_type .gt. 10) then
        write(pid_param_new, '(A20,I11,A)') ' SK_SCALE_MODE  ', PINPT%slater_koster_type, '  # 11 = NRL type scaling'
        write(pid_param_new, '(A20,2x,F9.4,A)') ' L_BROADEN  ', PINPT%l_broaden, '  # angstrom '
      else
        write(pid_param_new, '(A20,I11,A)') ' SK_SCALE_MODE  ', PINPT%slater_koster_type
      endif
    else
      write ( pid_param_new, '(A)' ) trim ( title )
    endif

    write(pid_param_new, '(A)') ' '

    do i = 1, PINPT%nparam
      if(PINPT%slater_koster_type .gt. 10) then
        k = nint(PINPT%param_const_nrl(1,1,i))
        nsub = PINPT%param_nsub(i)
        if( k .eq. 0) then ! if do not have pair
          if( nint(PINPT%param_const_nrl(4,1,i)) .eq. 1) then ! if fixed
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,*(2x,F20.8))', ADVANCE='NO' ) trim(PINPT%param_name(i)), &
                                                                PINPT%param_const_nrl(5,1:nsub,i)
              write ( pid_param_new, '(A)')'  Fixed'
            else
              write ( pid_param_new, '(i6, 2x,A20,*(2x,F20.8))',ADVANCE='NO' ) i, trim(PINPT%param_name(i)), &
                                                                    PINPT%param_const_nrl(5,1:nsub,i)
              write ( pid_param_new, '(A)')'  Fixed'
            endif
          else
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,*(2x,F20.8))', ADVANCE='NO' ) trim(PINPT%param_name(i)), &
                                                                PINPT%param_nrl(1:nsub,i)
              write ( pid_param_new, '(A)')'  '
            else
              write ( pid_param_new, '(i6, 2x,A20,*(2x,F20.8))', ADVANCE='NO' ) i, trim(PINPT%param_name(i)), &
                                                                    PINPT%param_nrl(1:nsub,i)
              write ( pid_param_new, '(A)')'  '
            endif
          endif
        elseif( k .ge. 1) then
          if( nint( PINPT%param_const_nrl(4,1,k) ) .eq. 1) then
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,*(2x,F20.8))', ADVANCE='NO' ) trim(PINPT%param_name(i)), &
                           PINPT%param_const_nrl(5,1:nsub,k)
              write ( pid_param_new, '(A)')'  Fixed'
            else
              write ( pid_param_new, '(i6, 2x,A20,*(2x,F20.8))', ADVANCE='NO' ) i, trim(PINPT%param_name(i)), &
                           PINPT%param_const_nrl(5,1:nsub,k)
              write ( pid_param_new, '(A)')'  Fixed'
            endif
          else
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,*(2x,F20.8))',ADVANCE='NO' ) trim(PINPT%param_name(i)), &
                           PINPT%param_nrl(1:nsub,k)
              write ( pid_param_new, '(A)')'  '
            else
              write ( pid_param_new, '(i6, 2x,A20,*(2x,F20.8))',ADVANCE='NO' ) i, trim(PINPT%param_name(i)), &
                           PINPT%param_nrl(1:nsub,k)
              write ( pid_param_new, '(A)')'  '
            endif
          endif
        endif

      else ! if not NRL type hopping 

        if( nint(PINPT%param_const(1,i)) .eq. 0) then ! if do not have pair
          if( nint(PINPT%param_const(4,i)) .eq. 1) then ! if fixed
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), PINPT%param_const(5,i), '  Fixed'
            else
              write ( pid_param_new, '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), PINPT%param_const(5,i), '  Fixed'
            endif
          else
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), PINPT%param(i),'  '
            else
              write ( pid_param_new, '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), PINPT%param(i),'  '
            endif
          endif
        elseif( nint(PINPT%param_const(1,i)) .ge. 1) then
          if( nint( PINPT%param_const( 4,nint(PINPT%param_const(1,i)) ) ) .eq. 1) then
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), &
                                                             PINPT%param_const(5,nint(PINPT%param_const(1,i))),'  Fixed'
            else
              write ( pid_param_new, '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), &
                                                                    PINPT%param_const(5,nint(PINPT%param_const(1,i))),'  Fixed'
            endif
          else
            if( .not. PINPT%flag_pfile_index) then
              write ( pid_param_new, '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), &
                                                             PINPT%param( nint(PINPT%param_const(1,i)) ),'  '
            else
              write ( pid_param_new, '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), &
                                                                    PINPT%param( nint(PINPT%param_const(1,i)) ),'  '
            endif
          endif
        endif
      endif
    enddo

    if(.not. flag_print_param) write ( *, '(a)' ) ' '

  endif

  if(flag_print_param) close(pid_param_new)

return
endsubroutine
