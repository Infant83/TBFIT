#include "alias.inc"
subroutine read_input(PINPT, PPRAM, PKPTS, PGEOM, PWGHT, EDFT, NN_TABLE, PINPT_DOS, PINPT_BERRY, PKAIA, PRPLT, PUFLD, mysystem)
  use parameters
  use read_incar
  use berry_phase
  use mpi_setup
  use print_io
  use projected_band
  use set_default
  implicit none
  integer*4                 mpierr
  integer*4                 i_continue
  integer*4                 linecount
  integer*4                 mysystem
  character*132             inputline
  character*40              desc_str
  integer*4, external    :: nitems
  type(incar)            :: PINPT
  type(params)           :: PPRAM
  type(dos)              :: PINPT_DOS
  type(berry)            :: PINPT_BERRY
  type(poscar)           :: PGEOM
  type(kpoints)          :: PKPTS
  type(energy)           :: EDFT
  type(weight)           :: PWGHT 
  type(hopping)          :: NN_TABLE
  type(gainp)            :: PKAIA
  type(replot)           :: PRPLT
  type(unfold)           :: PUFLD
  character(*),parameter :: func = 'read_input'
  logical                   flag_write_nntable
  
  flag_write_nntable = (.not. PINPT%flag_python_module)

  call init_incar(PINPT)
  call init_params(PPRAM, PINPT)
  call init_berry(PINPT_BERRY, PINPT)
  call init_kpoints(PKPTS, PINPT)
  call init_weight(PWGHT)
  call init_poscar(PGEOM)
  call init_hopping(NN_TABLE)
  call init_gainp(PKAIA)
  call init_replot(PRPLT)
  call init_dos(PINPT_DOS, PINPT)
  call init_unfold(PUFLD, PINPT)
  ! READ INCAR-TB
  call read_input_tags(PINPT, PPRAM, PKPTS, PGEOM, NN_TABLE, PWGHT, &
                       EDFT, PKAIA, PINPT_BERRY, PINPT_DOS, PRPLT, PUFLD, mysystem)
  ! SET BASIC SYSTEM
  call read_tb_param(PINPT, PPRAM, PWGHT)
  call read_geometry(PGEOM, PINPT, NN_TABLE, PPRAM)
  call read_kpoint(PKPTS, PINPT, PGEOM)
  call set_nn_table(NN_TABLE, PINPT, PPRAM, PGEOM, flag_write_nntable, PRPLT%flag_replot)
  call set_target_and_weight(EDFT, PWGHT, PINPT, PGEOM, PKPTS)

  ! SET input arguments for post-processings
  call set_ngrid(PINPT, PGEOM)
  call set_ldos_atom(PINPT, PGEOM)

  if(PINPT%flag_tbfit) PINPT%flag_distribute_nkp = .FALSE.

  if(PINPT%flag_tbfit .and. PINPT%flag_print_energy_diff ) then 
    if(PINPT%flag_print_orbital) then 
      PINPT%flag_print_energy_diff = .false.
      write(message,'(A)')'    !WARN! You have requested both (LORBIT .TRUE.) and (PRTDIFF .TRUE.)' ; write_msgi
      write(message,'(A)')'    !WARN! This two function does not work togeter. We deactivate PRTDIFF by default.' ; write_msgi
    endif
  else
    PINPT%flag_print_energy_diff = .false.
  endif

  write(message, '(A)' )' ' ; write_msgi
  write(message, '(A)' )' -------------------------------------------------------------';  write_msgi
  write(message, '(2A)')' #-- END READING INPUT FILE : ', trim(PINPT%ifilenm(mysystem)) ;  write_msgi
  write(message, '(A)' )' -------------------------------------------------------------';  write_msgi
  write(message, '(A)' )' ' ; write_msgi

return
endsubroutine

subroutine rewrite_incar(PINPT, PPRAM, PWGHT)
   use parameters, only : incar, weight, pid_incar, params
   type(incar  ) :: PINPT
   type(params ) :: PPRAM
   type(weight ) :: PWGHT
   integer*4        pid_incar_temp
   integer*4        i_continue, i_continue_
   integer*4        l0, i
   character*256    inputline
   character*256    inputline_temp
   character*40     desc_str
   character*40     fm_wt, fm_ob
   character*20,external ::   int2str
   character*132    gnu_command
   logical          flag_written_weight
   character*132    fname, fname_temp

   i_continue = 0
   i_continue_= 0
   flag_written_weight = .false. 
   
   fname = trim(PINPT%ifilenm(PWGHT%mysystem))
   fname_temp = '__'//trim(fname)

   write(fm_wt,'( "(A11,A",A,",A9,A",A,",A9,A",A,",A9,A10)" )') trim(ADJUSTL(int2str(PWGHT%max_len_strip_kp))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_tb))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_df)))
   if( PWGHT%max_len_strip_ob .gt. 1)then
     write(fm_ob,'( "(A11,A",A,",A9,A",A,",A9,A",A,",A9,A",A,",A9,A10)" )') &
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_kp))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_tb))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_ob))),&
                                                              trim(ADJUSTL(int2str(PWGHT%max_len_strip_st)))
   else
     write(fm_ob,'( "(A11,A10",",A9,A10",",A9,A10",",A9,A10",",A9,A10)" )')
   endif

   pid_incar_temp = pid_incar + 10
   open(pid_incar, file=trim(fname), status = 'unknown')
   open(pid_incar_temp, file=trim(fname_temp), status = 'unknown')
    
   do while (i_continue .ge. 0)
     read(pid_incar, '(A)',iostat=i_continue) inputline
     if(i_continue .lt. 0) exit
     if(index(inputline,'WEIGHT') .gt. 1 .and. index(inputline,'SET') .ge. 1 .and. &
        index(inputline,'WEIGHT') .gt. index(inputline,'SET') .and.  &
        index(adjustl(trim(inputline)),'#') .ne. 1) then 
 
        write(pid_incar_temp,'(A)')trim(inputline)

        ! insert weight info of PFILE to INCAR-TB
        write(pid_incar_temp,'(3A)')' #  Weight info of PFILE(',trim(PPRAM%pfilenm),') is written upon the USE_WEIGHT request'
        if(PWGHT%nweight .ge. 1) then
          do i = 1, PWGHT%nweight
            write(pid_incar_temp,fm_wt)'    KRANGE ',trim(PWGHT%strip_kp(i)), '  TBABND ',trim(PWGHT%strip_tb(i)), &
                                       '  DFTBND ',trim(PWGHT%strip_df(i)), '  WEIGHT ',trim(PWGHT%strip_wt(i))
          enddo
        endif
        if(PWGHT%npenalty_orb .ge. 1) then
          do i = 1, PWGHT%npenalty_orb
            write(pid_incar_temp,fm_ob)'    KRANGE ',trim(PWGHT%strip_kp_orb(i)), '  TBABND ',trim(PWGHT%strip_tb_orb(i)), &
                                       '  ORBT_I ',trim(PWGHT%strip_orb(i)),    '  SITE_I ',trim(PWGHT%strip_site(i)), &
                                       '  PENALTY ',trim(PWGHT%strip_pen_orb(i))
          enddo
        endif
  
        write(pid_incar_temp,'(A)')' '
        write(pid_incar_temp,'(A)')' # The existing weight information was commented out.'
        i_continue_ = 0
        do while (i_continue_ .ge. 0) 
          read(pid_incar,'(A)',iostat=i_continue) inputline
          read(inputline, *,iostat=i_continue_) desc_str
          if(trim(desc_str) .eq. 'END') then 
            i_continue_ = -1
            write(pid_incar_temp,'(A)')trim(inputline)
            flag_written_weight = .true.
          elseif(trim(desc_str) .eq. 'KRANGE') then
            i_continue_ =  0
            write(pid_incar_temp,'(2A)')' # ',trim(inputline)
          endif
        enddo
     else
       write(pid_incar_temp,'(A)')trim(inputline)
     endif
   enddo
   
   if(.not. flag_written_weight) then
     write(pid_incar_temp,'(A)')' '
     write(pid_incar_temp,'(A)')'    SET  WEIGHT'
     ! insert weight info of PFILE to INCAR-TB
     write(pid_incar_temp,'(3A)')' #  Weight info of PFILE(',trim(PPRAM%pfilenm),') is written upon the USE_WEIGHT request'
     if(PWGHT%nweight .ge. 1) then
       do i = 1, PWGHT%nweight
         write(pid_incar_temp,fm_wt)'    KRANGE ',trim(PWGHT%strip_kp(i)), '  TBABND ',trim(PWGHT%strip_tb(i)), &
                                    '  DFTBND ',trim(PWGHT%strip_df(i)), '  WEIGHT ',trim(PWGHT%strip_wt(i))
       enddo
     endif
     if(PWGHT%npenalty_orb .ge. 1) then
       do i = 1, PWGHT%npenalty_orb
         write(pid_incar_temp,fm_ob)'    KRANGE ',trim(PWGHT%strip_kp_orb(i)), '  TBABND ',trim(PWGHT%strip_tb_orb(i)), &
                                    '  ORBT_I ',trim(PWGHT%strip_orb(i)),    '  SITE_I ',trim(PWGHT%strip_site(i)), &
                                    '  PENALTY ',trim(PWGHT%strip_pen_orb(i))
       enddo
     endif
     write(pid_incar_temp,'(A)')'    END  WEIGHT'
   endif

   close(pid_incar)
   close(pid_incar_temp)

   write(gnu_command, '(4A)')'yes | cp ',trim(fname), ' ', trim(fname)//'_pristine'
   call execute_command_line(gnu_command, .false.)
   write(gnu_command, '(4A)')'yes | mv -f ',trim(fname_temp),' ', trim(fname)
   call execute_command_line(gnu_command, .false.)

return
endsubroutine
