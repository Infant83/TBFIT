subroutine print_param (PINPT, title, flag_print_param)
  use parameters, only : incar
  use print_io
  implicit none
  integer*4 i,j,k,pid_param_new
  integer*4 nsub, im
  character ( len = * ) title
  logical flag_print_param
  type (incar)   :: PINPT
  character*40   fm_wt, fm_ob
  character*20,external ::   int2str
  write(fm_wt,'( "(A9,A",A,",A9,A",A,",A9,A",A,",A9,A10)" )') trim(ADJUSTL(int2str(PINPT%max_len_strip_kp))),&
                                                             trim(ADJUSTL(int2str(PINPT%max_len_strip_tb))),&
                                                             trim(ADJUSTL(int2str(PINPT%max_len_strip_df)))
  if( PINPT%max_len_strip_ob .gt. 1)then
    write(fm_ob,'( "(A9,A",A,",A9,A",A,",A9,A",A,",A9,A",A,",A9,A10)" )') &
                                                             trim(ADJUSTL(int2str(PINPT%max_len_strip_kp))),&
                                                             trim(ADJUSTL(int2str(PINPT%max_len_strip_tb))),&
                                                             trim(ADJUSTL(int2str(PINPT%max_len_strip_ob))),&
                                                             trim(ADJUSTL(int2str(PINPT%max_len_strip_st)))
  else
    write(fm_ob,'( "(A9,A10",",A9,A10",",A9,A10",",A9,A10",",A9,A10)" )')
  endif

  im = 0

  if(.not. flag_print_param) then
    pid_param_new = 6
  elseif(flag_print_param)then
    pid_param_new = 81
    open(pid_param_new, file = trim(title), status = 'unknown')
  endif

  im=im+1; write(message_pack(im), '(A)' ) ' '
  if(flag_print_param) then
    im=im+1; write(message_pack(im), '(A,A)' ) '# ',trim ( title )
    if(PINPT%flag_collinear) then
      im=im+1; write(message_pack(im),'(A,A)')'# E_TARGET FILE NAME (up): ',trim(PINPT%efilenmu)
      im=im+1; write(message_pack(im),'(A,A)')'# E_TARGET FILE NAME (dn): ',trim(PINPT%efilenmd)
    else
      im=im+1; write(message_pack(im),'(A,A)')'# E_TARGET FILE NAME : ',trim(PINPT%efilenmu)
    endif
    if(PINPT%flag_scissor) then 
      im=im+1; write(message_pack(im),'(A,F8.2,A,I5,A)')'# SCISSOR: .TRUE. => EDFT(n,k) + ',PINPT%r_scissor,' (if n >=',PINPT%i_scissor,')'
    endif
    if(PINPT%nweight .ge. 1) then
      do i = 1, PINPT%nweight
        im=im+1; write(message_pack(im),fm_wt)'# KRANGE ',trim(PINPT%strip_kp(i)), '  TBABND ',trim(PINPT%strip_tb(i)), '  DFTBND ',trim(PINPT%strip_df(i)), '  WEIGHT ',trim(PINPT%strip_wt(i))
      enddo
    endif
    if(PINPT%npenalty_orb .ge. 1) then
      do i = 1, PINPT%npenalty_orb
        im=im+1; write(message_pack(im),fm_ob)'# KRANGE ',trim(PINPT%strip_kp_orb(i)), '  TBABND ',trim(PINPT%strip_tb_orb(i)), '  ORBT_I ',trim(PINPT%strip_orb(i)), '  SITE_I ',trim(PINPT%strip_site(i)), '  PENALTY ',trim(PINPT%strip_pen_orb(i))
      enddo
    endif
    if(PINPT%ndegenw .ge. 1) then
      do i = 1, PINPT%ndegenw
        im=im+1; write(message_pack(im),fm_wt)'# KRANGE ',trim(PINPT%strip_kp_deg(i)), '  TBABND ',trim(PINPT%strip_tb_deg(i)), '  DFTBND ',trim(PINPT%strip_df_deg(i)), '  DEGENW ',trim(PINPT%strip_wt_deg(i))
      enddo
    endif

    
    if(PINPT%flag_pfile_index) then
      im=im+1; write(message_pack(im), '(A20,A11)') ' PRINT_INDEX  ','.TRUE.'
    else
      im=im+1; write(message_pack(im), '(A20,A11)') ' PRINT_INDEX  ','.FALSE.'
    endif

    if(PINPT%flag_use_overlap) then
      im=im+1; write(message_pack(im), '(A20,A11)') ' USE_OVERLAP  ','.TRUE.'
    else
      im=im+1; write(message_pack(im), '(A20,A11)') ' USE_OVERLAP  ','.FALSE.'
    endif

    if(PINPT%slater_koster_type .gt. 10) then
      im=im+1; write(message_pack(im), '(A20,I11,A)') ' SK_SCALE_MODE  ', PINPT%slater_koster_type, '  # 11 = NRL type scaling'
      im=im+1; write(message_pack(im), '(A20,2x,F9.4,A)') ' L_BROADEN  ', PINPT%l_broaden, '  # angstrom '
    else
      im=im+1; write(message_pack(im), '(A20,I11,A)') ' SK_SCALE_MODE  ', PINPT%slater_koster_type
    endif
  else
    im=im+1; write(message_pack(im), '(A)' ) trim ( title )
  endif

  im=im+1; write(message_pack(im), '(A)') ' '

  do i = 1, PINPT%nparam
    if(PINPT%slater_koster_type .gt. 10) then
      k = nint(PINPT%param_const_nrl(1,1,i))
      nsub = PINPT%param_nsub(i)
      if( k .eq. 0) then ! if do not have pair
        if( nint(PINPT%param_const_nrl(4,1,i)) .eq. 1) then ! if fixed
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,*(2x,F20.8))') trim(PINPT%param_name(i)), PINPT%param_const_nrl(5,1:nsub,i)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  Fixed'
          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,*(2x,F20.8))') i, trim(PINPT%param_name(i)), PINPT%param_const_nrl(5,1:nsub,i)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  Fixed'
          endif
        else
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,*(2x,F20.8))') trim(PINPT%param_name(i)), PINPT%param_nrl(1:nsub,i)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  '
          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,*(2x,F20.8))') i, trim(PINPT%param_name(i)), PINPT%param_nrl(1:nsub,i)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  '
          endif
        endif
      elseif( k .ge. 1) then
        if( nint( PINPT%param_const_nrl(4,1,k) ) .eq. 1) then
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,*(2x,F20.8))') trim(PINPT%param_name(i)), PINPT%param_const_nrl(5,1:nsub,k)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  Fixed'

          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,*(2x,F20.8))') i, trim(PINPT%param_name(i)), PINPT%param_const_nrl(5,1:nsub,k)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  Fixed'

          endif
        else
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,*(2x,F20.8))') trim(PINPT%param_name(i)), PINPT%param_nrl(1:nsub,k)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  '
          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,*(2x,F20.8))') i, trim(PINPT%param_name(i)), PINPT%param_nrl(1:nsub,k)
                     write(message_pack(im), '(2A)')trim(message_pack(im)),'  '
          endif
        endif
      endif

    else ! if not NRL type hopping 

      if( nint(PINPT%param_const(1,i)) .eq. 0) then ! if do not have pair
        if( nint(PINPT%param_const(4,i)) .eq. 1) then ! if fixed
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), PINPT%param_const(5,i), '  Fixed'
          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), PINPT%param_const(5,i), '  Fixed'
          endif
        else
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), PINPT%param(i),'  '
          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), PINPT%param(i),'  '
          endif
        endif
      elseif( nint(PINPT%param_const(1,i)) .ge. 1) then
        if( nint( PINPT%param_const( 4,nint(PINPT%param_const(1,i)) ) ) .eq. 1) then
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), PINPT%param_const(5,nint(PINPT%param_const(1,i))),'  Fixed'
          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), PINPT%param_const(5,nint(PINPT%param_const(1,i))),'  Fixed'
          endif
        else
          if( .not. PINPT%flag_pfile_index) then
            im=im+1; write(message_pack(im), '(2x,A20,2x,F16.8,A)' ) trim(PINPT%param_name(i)), PINPT%param( nint(PINPT%param_const(1,i)) ),'  '
          else
            im=im+1; write(message_pack(im), '(i6, 2x,A20,2x,F16.8,A)' ) i, trim(PINPT%param_name(i)), PINPT%param( nint(PINPT%param_const(1,i)) ),'  '
          endif
        endif
      endif
    endif
  enddo

  if(.not. flag_print_param) then 
    im=im+1; write(message_pack(im), '(a)' ) ' '
  endif

  do i = 1, im
    
    if(.not. flag_print_param) then 
      call write_log(trim(message_pack(i)), 3, 0)
    elseif(flag_print_param) then
      write(pid_param_new,'(A)')trim(message_pack(i))
    endif
  enddo

  if(flag_print_param) close(pid_param_new)

return
endsubroutine
