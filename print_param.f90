subroutine print_param (PINPT, PPRAM, PWGHT, title, flag_print_file)
  use parameters, only : incar, weight, params
  use print_io
  implicit none
  type (incar)          :: PINPT
  type (weight)         :: PWGHT
  type (params)         :: PPRAM
  integer*4                i,j,k,pid_param_new
  integer*4                nsub, im, i_fix
  character ( len = * )    title
  character ( len = 40)    param_name
  real*8                   param(PPRAM%param_nsub_max)
  logical                  flag_print_file
  character*40             fm_wt, fm_ob
  character*20,external :: int2str
  write(fm_wt,'( "(A9,A",A,",A9,A",A,",A9,A",A,",A9,A10)" )') trim(ADJUSTL(int2str(PWGHT%max_len_strip_kp))),&
                                                             trim(ADJUSTL(int2str(PWGHT%max_len_strip_tb))),&
                                                             trim(ADJUSTL(int2str(PWGHT%max_len_strip_df)))
  if( PWGHT%max_len_strip_ob .gt. 1)then
    write(fm_ob,'( "(A9,A",A,",A9,A",A,",A9,A",A,",A9,A",A,",A9,A10)" )') &
                                                             trim(ADJUSTL(int2str(PWGHT%max_len_strip_kp))),&
                                                             trim(ADJUSTL(int2str(PWGHT%max_len_strip_tb))),&
                                                             trim(ADJUSTL(int2str(PWGHT%max_len_strip_ob))),&
                                                             trim(ADJUSTL(int2str(PWGHT%max_len_strip_st)))
  else
    write(fm_ob,'( "(A9,A10",",A9,A10",",A9,A10",",A9,A10",",A9,A10)" )')
  endif

  im = 0

  if(.not. flag_print_file) then
    pid_param_new = 6
  elseif(flag_print_file)then
    pid_param_new = 81
    open(pid_param_new, file = trim(title), status = 'unknown')
  endif

  im=im+1; write(message_pack(im), '(A)' ) ' '
  if(flag_print_file) then
    im=im+1; write(message_pack(im), '(A,A)' ) '# ',trim ( title )
    if(PINPT%flag_collinear) then
      im=im+1; write(message_pack(im),'(A,A)')'# E_TARGET FILE NAME (up): ',trim(PWGHT%efilenmu)
      im=im+1; write(message_pack(im),'(A,A)')'# E_TARGET FILE NAME (dn): ',trim(PWGHT%efilenmd)
    else
      im=im+1; write(message_pack(im),'(A,A)')'# E_TARGET FILE NAME : ',trim(PWGHT%efilenmu)
    endif
    if(PINPT%flag_scissor) then 
      im=im+1; write(message_pack(im),'(A,F8.2,A,I5,A)')'# SCISSOR: .TRUE. => EDFT(n,k) + ',PINPT%r_scissor, &
                                                                                ' (if n >=',PINPT%i_scissor,')'
    endif

    if( .not. PINPT%flag_set_weight_from_file) then
      if(PWGHT%nweight .ge. 1) then
        do i = 1, PWGHT%nweight
          im=im+1; write(message_pack(im),fm_wt)'# KRANGE ',trim(PWGHT%strip_kp(i)), &
                                                '  TBABND ',trim(PWGHT%strip_tb(i)), &
                                                '  DFTBND ',trim(PWGHT%strip_df(i)), &
                                                '  WEIGHT ',trim(PWGHT%strip_wt(i))
        enddo
      endif
      if(PWGHT%npenalty_orb .ge. 1) then
        do i = 1, PWGHT%npenalty_orb
          im=im+1; write(message_pack(im),fm_ob)'# KRANGE ',trim(PWGHT%strip_kp_orb(i)), &
                                                '  TBABND ',trim(PWGHT%strip_tb_orb(i)), &
                                                '  ORBT_I ',trim(PWGHT%strip_orb(i)), &
                                                '  SITE_I ',trim(PWGHT%strip_site(i)), &
                                                '  PENALTY ',trim(PWGHT%strip_pen_orb(i))
        enddo
      endif
      if(PWGHT%ndegenw .ge. 1) then
        do i = 1, PWGHT%ndegenw
          im=im+1; write(message_pack(im),fm_wt)'# KRANGE ',trim(PWGHT%strip_kp_deg(i)), &
                                                '  TBABND ',trim(PWGHT%strip_tb_deg(i)), &
                                                '  DFTBND ',trim(PWGHT%strip_df_deg(i)), &
                                                '  DEGENW ',trim(PWGHT%strip_wt_deg(i))
        enddo
      endif

    elseif(PINPT%flag_set_weight_from_file) then
      im=im+1; write(message_pack(im),'(A,A)')'# WEIGHT FILE NAME : ', trim(PWGHT%wfilenm)
    endif

    
    if(PPRAM%flag_pfile_index) then
      im=im+1; write(message_pack(im), '(A20,A11)') ' PRINT_INDEX  ','.TRUE.'
    else
      im=im+1; write(message_pack(im), '(A20,A11)') ' PRINT_INDEX  ','.FALSE.'
    endif

    if(PPRAM%flag_use_overlap) then
      im=im+1; write(message_pack(im), '(A20,A11)') ' USE_OVERLAP  ','.TRUE.'
      im=im+1; write(message_pack(im), '(A20,F12.8)')' REDUCE_OVERLAP  ',PPRAM%reduce_overlap
      im=im+1; write(message_pack(im), '(A20,F12.8)')' REDUCE_HOPPING  ',PPRAM%reduce_hopping
    else
      im=im+1; write(message_pack(im), '(A20,A11)') ' USE_OVERLAP  ','.FALSE.'
      im=im+1; write(message_pack(im), '(A20,F12.8)')' REDUCE_HOPPING  ',PPRAM%reduce_hopping
    endif

    if(PPRAM%flag_slater_koster) then
      im=im+1; write(message_pack(im), '(A20,A11)') '       IS_SK  ','.TRUE.'
    else
      im=im+1; write(message_pack(im), '(A20,A11)') ' USE_OVERLAP  ','.FALSE.'
    endif

    if(PPRAM%slater_koster_type .gt. 10) then
      im=im+1; write(message_pack(im), '(A20,I11,A)') ' SK_SCALE_MODE  ', PPRAM%slater_koster_type, '  # 11 = NRL type scaling'
      im=im+1; write(message_pack(im), '(A20,2x,F9.4,A)') ' L_BROADEN  ', PPRAM%l_broaden, '  # angstrom '
    else
      im=im+1; write(message_pack(im), '(A20,I11,A)') ' SK_SCALE_MODE  ', PPRAM%slater_koster_type
    endif
  else
    im=im+1; write(message_pack(im), '(A)' ) trim ( title )
  endif

  im=im+1; write(message_pack(im), '(A)') ' '

  do i = 1, PPRAM%nparam
    call param_select(PPRAM, i, i_fix, param_name, param, nsub)
    call set_print_message(  i, i_fix, param_name, param, nsub, PPRAM%flag_pfile_index, message)
    im=im+1 ; message_pack(im) = message
  enddo

  if(.not. flag_print_file) then 
    im=im+1; write(message_pack(im), '(a)' ) ' '
  endif

  do i = 1, im
    
    if(.not. flag_print_file) then 
      call write_log(trim(message_pack(i)), print_mode, 0)
    elseif(flag_print_file) then
      write(pid_param_new,'(A)')trim(message_pack(i))//'  '
    endif
  enddo

  if(flag_print_file) close(pid_param_new)

return
endsubroutine

subroutine print_param_pso (PINPT, PPRAM, PWGHT)
  use parameters, only : incar, weight, params
  use print_io
  use sorting
  use directory
  implicit none
  type (incar)          :: PINPT
  type (weight)         :: PWGHT
  type (params)         :: PPRAM
  integer*4                i,j,k,ik, pid_param_pso
  integer*4                iparam, nsub, im, i_fix, bestn
  character ( len = 80)    title, title_pid
  character ( len = 40)    param_name
  real*8                   param(PPRAM%param_nsub_max)
  real*8                   last_costs(PPRAM%pso_nparticles)
  real*8                   cost
  integer*4                ilast_costs(PPRAM%pso_nparticles)
  integer*4, allocatable:: ibest_costs(:)
  logical                  flag_print_file
  character*40             fm_wt, fm_ob
  character*20,external :: int2str
  real*8                   param_temp(PPRAM%nparam)
  real*8                   enorm, dist
  real*8, allocatable   :: dist_min(:,:)
  integer*4,allocatable :: dist_nn(:,:)
  external                 enorm
  character*13          :: pso_folder

  title = 'PARAM_FIT.PSO.best_particles.dat'
  pso_folder= './PSO_params/'
  pid_param_pso = 181
  k = 3
  last_costs = 0d0
  param_temp = 0d0 ; param_temp = PPRAM%param
  bestn = int(real(PPRAM%pso_nparticles) * PPRAM%pso_report_ratio)
  if(bestn .le. 0) bestn = PPRAM%pso_nparticles

  allocate(ibest_costs(bestn))
  allocate(dist_nn(bestn,k))
  allocate(dist_min(bestn,k))

  last_costs = PPRAM%pso_cost_history_i(PPRAM%pso_miter, :)

  call get_sort_index_1D(ilast_costs, last_costs, PPRAM%pso_nparticles, 'ascending ')
  ibest_costs = ilast_costs(1:bestn)

  call makedir(pso_folder)
 !open(pid_param_pso, file = trim(pso_folder)//trim(title), status = 'unknown')
  open(pid_param_pso, file = trim(title), status = 'unknown')
  do i=1, bestn
    j = ibest_costs(i)
    cost = last_costs(j)
    write(pid_param_pso,'(A,I0,A,I0,A,F20.9)')'# Best particle : rank= ',i, ', ID= ',j, ', COST= ', cost

    PPRAM%param(PPRAM%iparam_free(:)) = PPRAM%pso_pbest_history(PPRAM%pso_miter, j, :)

    do iparam = 1, PPRAM%nparam
      call param_select(PPRAM, iparam, i_fix, param_name, param, nsub)
      if(i_fix .eq. 1) then
        write(pid_param_pso, '(I6,F20.8,A,I0,F20.9)') iparam, param(1:nsub), '   #   Fixed  # ', j, cost 
      else
        write(pid_param_pso, '(I6,F20.8,A,I0,F20.9)') iparam, param(1:nsub), '   #          # ', j, cost
      endif
    enddo
  
    write(pid_param_pso,'(A)')'  '
    write(pid_param_pso,'(A)')'  '
 
   !write(title_pid,'(A,A,I0,A,I0,A)')trim(pso_folder),'PARAM_FIT.PSO.RANK.',i,'.PID.',j,'.dat'
    write(title_pid,'(A,A,I0,A)')trim(pso_folder),'PARAM_FIT.PSO.RANK.',i,'.dat'
    call print_param (PINPT, PPRAM, PWGHT, trim(title_pid), .TRUE.)
  enddo

  ! calculate parameter distance between each parameter sets
ii:do i=1, bestn
    dist_min(i,:) = huge(dist)
    dist_nn(i,:)  = -1
 jj:do j=1, bestn
      if (i == j) cycle
      dist = enorm(PPRAM%nparam_free, PPRAM%pso_pbest_history(PPRAM%pso_miter, ibest_costs(i), :)-PPRAM%pso_pbest_history(PPRAM%pso_miter, ibest_costs(j), :))
   kk:do ik=1, k
        if(dist .lt. dist_min(i,ik) .and. ik .lt. k) then 
          dist_min(i,ik+1:k)= dist_min(i,ik:k-1)
          dist_nn(i,ik+1:k) = dist_nn(i,ik:k-1)
          dist_min(i,ik)= dist
          dist_nn(i,ik) = ibest_costs(j)
          cycle jj
        elseif(dist .lt. dist_min(i,ik) .and. ik .eq. k) then 
          dist_min(i,ik)= dist
          dist_nn(i,ik) = ibest_costs(j)
          cycle jj
        endif
      enddo kk
    enddo jj
  ! write(6,'(A,I0,3(F16.8),3(I6))')"  DISTANCE BETWEEN PARAMETERS: (PARAM_I, DIST(IJ), PARAM_J)", ibest_costs(i), dist_min(i,:), dist_nn(i,:)
  enddo ii

!write(6,*)"ZZZ", PPRAM%pso_pbest_history(PINPT%miter, ibest_costs(1), :), PPRAM%pso_pbest_history(PINPT%miter, ibest_costs(3), :)
!write(6,*)"ZZZ", PPRAM%pso_pbest_history(PINPT%miter, ibest_costs(3), :)
  close(pid_param_pso)
  PPRAM%param = param_temp

  return
endsubroutine

subroutine param_select (PPRAM, i, i_fix, param_name, param, nsub)

   use parameters, only : params
   implicit none
   type (params)  :: PPRAM
   real*8            param(PPRAM%param_nsub_max)
   character*40      param_name
   integer*4         k, i, nsub
   integer*4         i_fix
   character*80      fmt
   character*2048    msg

    if(PPRAM%slater_koster_type .gt. 10) then
      k = nint(PPRAM%param_const_nrl(1,1,i))
      nsub = PPRAM%param_nsub(i)
      param_name = trim(PPRAM%param_name(i))
      if( k .eq. 0) then ! if do not have pair
        i_fix = nint(PPRAM%param_const_nrl(4,1,i))
        if( i_fix .eq. 1) then ! if fixed
          param      = PPRAM%param_const_nrl(5,1:nsub,i)
        else
          param      = PPRAM%param_nrl(1:nsub,i)
        endif
      elseif( k .ge. 1) then
        i_fix = nint( PPRAM%param_const_nrl(4,1,k) )  ! if fixed
        if( nint( PPRAM%param_const_nrl(4,1,k) ) .eq. 1) then
          param      = PPRAM%param_const_nrl(5,1:nsub,k)
        else
          param      = PPRAM%param_nrl(1:nsub,k)
        endif
      endif
    else ! if not NRL type hopping
      k = nint(PPRAM%param_const(1,i)) ! ! if do not have pair
      nsub = 1
      param_name = trim(PPRAM%param_name(i))
      if( k .eq. 0) then ! if do not have pair
        i_fix = nint(PPRAM%param_const(4,i))
        if( i_fix .eq. 1) then ! if fixed
          param      = PPRAM%param_const(5,i)
        else
          param      = PPRAM%param(i)
        endif
      elseif( k .ge. 1) then ! if have pair
        i_fix = nint(PPRAM%param_const(4,k))
        if( i_fix .eq. 1) then
          param      = PPRAM%param_const(5,k)
        else
          param      = PPRAM%param( nint(PPRAM%param_const(1,i)) )
        endif
      endif
    endif

    ! applying constraint that monotonic decrease rule
    if( (param_name(1:2) .eq. 's_' .or. param_name(1:3) .eq. 'os_') .and. &
        (PPRAM%slater_koster_type .eq. 4 .or. PPRAM%slater_koster_type .eq. 6) ) param = abs(param)

return
endsubroutine

subroutine set_print_param_format(fmt, nsub, flag_index)
   implicit none
   character*80  fmt
   character*10  fmt_
   character*6   c_index
   character*11  c_param
   character*3   c_fix  
   logical       flag_index
   integer*4     nsub

   if(flag_index) then
     c_index = 'i6,2x,'
   else
     c_index = '   2x,'
   endif
    
   if(nsub .lt. 1) then
     c_param = '*(2x,E20.8)'
   else
     c_param = '  2x,E20.8 '
   endif

   write(fmt, '(5A)')'(',trim(c_index),' A20, ', c_param, ' ,A)'

return
endsubroutine

subroutine set_print_message(i, i_fix, param_name, param, nsub, flag_index, msg)

   implicit none
   integer*4      i, i_fix, nsub
   logical        flag_index, flag_fix
   character*80   fmt
   character*40   param_name
   real*8         param(nsub)
   character*2048 msg 

   if(i_fix .eq. 1) then
     flag_fix = .true.
   else
     flag_fix = .false.
   endif

   call set_print_param_format(fmt, nsub, flag_index)

   if(flag_index) then
     if(flag_fix) then
       write(msg, fmt) i, trim(param_name), param(1:nsub), '       Fixed'
     else
       write(msg, fmt) i, trim(param_name), param(1:nsub), '       '
     endif 
   else
     if(flag_fix) then
       write(msg, fmt)    trim(param_name), param(1:nsub), '       Fixed'
     else
       write(msg, fmt)    trim(param_name), param(1:nsub), '       '
     endif
  endif
return
endsubroutine
