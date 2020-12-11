#include "alias.inc"
subroutine read_kpoint(PKPTS, PINPT, PGEOM)
  use parameters, only : kpoints, incar, poscar
  use mpi_setup
  use print_io
  implicit none
  type(incar)   :: PINPT
  type(kpoints) :: PKPTS
  type(poscar)  :: PGEOM
  integer*4        mpierr
  logical          flag_kfile_exist

  ! read info: kpoint 
  inquire(file=PKPTS%kfilenm,exist=flag_kfile_exist)
  if(flag_kfile_exist) then
    call read_kpoint_file(PKPTS,PGEOM,PINPT%flag_ndiv_line_parse,PINPT%flag_ndiv_grid_parse)
  elseif( (.not. flag_kfile_exist .and. PINPT%flag_get_band) .or. &
          (.not. flag_kfile_exist .and. PINPT%flag_get_berry_curvature) ) then
    write(message,'(A,A,A)')'  !WARN! ',trim(PKPTS%kfilenm),' does not exist!! Exit...' ; write_msgi
    kill_job
  endif

  return
endsubroutine
subroutine read_kpoint_file(PKPTS, PGEOM, flag_ndiv_line_parse, flag_ndiv_grid_parse)
  use parameters, only : kpoints, poscar, pid_kpoint, incar
  use mpi_setup
  use print_io
  implicit none
  integer*4, parameter :: max_kline=100
  integer*4     i_continue, nitems,ndiv_temp
  integer*4     i,linecount, i_dummy
  integer*4     idiv_mode, ik
  integer*4     mpierr
  integer*4     iline
  real*8        kline_dummy(3,max_kline)
  real*8, allocatable :: kpts_cart(:,:), kpts_reci(:,:)
  character*132 inputline,fname
  character*40  desc_str,dummy, k_name_dummy(max_kline)
  character*20  int2str
  character(*), parameter :: func = 'read_kpoint'
  logical       flag_skip
  external      nitems,int2str
  logical       flag_ndiv_line_parse, flag_ndiv_grid_parse

  type(kpoints) :: PKPTS
  type(poscar)  :: PGEOM
  PKPTS%flag_cartesianK = .false.
  PKPTS%flag_reciprocal = .false.
  PKPTS%flag_kgridmode  = .false.
  PKPTS%flag_klinemode  = .false.
  fname                 = trim(PKPTS%kfilenm)

  write(message,*)' '  ; write_msgi
  write(message,*)'*- READING KPOINTS FILE: ',trim(fname)  ; write_msgi
  open (pid_kpoint, FILE=fname,iostat=i_continue)
  linecount = 0 

 line: do
        read(pid_kpoint,'(A)',iostat=i_continue) inputline
        if(i_continue<0) exit               ! end of file reached
        if(i_continue>0) then
          write(message,*)'Unknown error reading file:',trim(fname),func  ; write_msgi
          kill_job
        endif
        linecount = linecount + 1
        call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle
        call check_empty(inputline,linecount,i,flag_skip) ; if(flag_skip) cycle

        if(i_continue .ne. 0) cycle              ! skip empty line

        ! head
         if(linecount .eq. 1) then
           cycle

        ! dividing factor
         elseif(linecount .eq. 2) then
           PKPTS%n_ndiv = nitems(inputline)
           read(inputline,*,iostat=i_continue) ndiv_temp
           if(ndiv_temp .ne. 0) then
             if(.not. allocated(PKPTS%ndiv)) allocate( PKPTS%ndiv(PKPTS%n_ndiv) )
             read(inputline,*,iostat=i_continue) PKPTS%ndiv(1:PKPTS%n_ndiv)
             if(PKPTS%n_ndiv .gt. 1) PKPTS%kline_type = 'FHI-AIMS' ! enforce 
           endif
           cycle

        ! k-grid or -line mode
         elseif(linecount .eq. 3) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'L' .or. desc_str(1:1) .eq. 'l') then
             write(message,'(A)')'   K_MODE: Line-mode'  ; write_msgi
             PKPTS%flag_klinemode=.true.
             if(.not. flag_ndiv_line_parse) then
               write(message,'(A,*(I8))')'   #N_DIV: (number of ndiv)',PKPTS%n_ndiv  ; write_msgi
               write(message,'(A,*(I8))')'    N_DIV:',PKPTS%ndiv  ; write_msgi
             else
               write(message,'(A,*(I8))')'    N_DIV: (set by -nkp_line option) ',PKPTS%ndiv  ; write_msgi
             endif
           elseif(desc_str(1:1) .eq. 'G' .or. desc_str(1:1) .eq. 'g') then
             linecount = linecount + 1
             write(message,'(A)')'   K_MODE: Gamma-centered'  ; write_msgi
             PKPTS%flag_gamma=.true.
             if(.not. allocated(PKPTS%ndiv)) allocate( PKPTS%ndiv(3) )
           elseif(desc_str(1:1) .eq. 'M' .or. desc_str(1:1) .eq. 'm') then
             linecount = linecount + 1
             write(message,'(A)')'   K_MODE: non Gamma-centered'  ; write_msgi
             PKPTS%flag_gamma=.false.
             if(.not. allocated(PKPTS%ndiv)) allocate( PKPTS%ndiv(3) )
           endif

         ! k-vector type
         elseif(linecount .eq. 4) then
           read(inputline,*,iostat=i_continue) desc_str
           if( (desc_str(1:1) .eq. 'R' .or. desc_str(1:1) .eq. 'r') .and. &
                PKPTS%flag_klinemode ) then 
             PKPTS%flag_reciprocal=.true.
             PKPTS%flag_cartesianK=.false.
             write(message,'(A)')'   K_TYPE: Reciprocal unit'  ; write_msgi
           elseif( (desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c'  .or.  &
                    desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') .and. &
                    PKPTS%flag_klinemode ) then
             PKPTS%flag_reciprocal=.false.
             PKPTS%flag_cartesianK=.true.
             write(message,'(A)')'   K_TYPE: Cartesian unit (1/A)'  ; write_msgi
           endif

         ! k-grid if .not. 'linemode' .and. 'kgridmode)
         elseif(linecount .eq. 5 .and. .not. PKPTS%flag_klinemode) then
           if(.not. flag_ndiv_grid_parse) then
             read(inputline,*,iostat=i_continue) PKPTS%ndiv(1:3)
             write(message,'(A,4x,3I4)')'   K_GRID:',PKPTS%ndiv(1:3)  ; write_msgi
           else
             read(inputline,*,iostat=i_continue) i_dummy, i_dummy, i_dummy
             write(message,'(A,4x,3I4)')'   K_GRID: (set by -nkp_grid option)',PKPTS%ndiv(1:3)  ; write_msgi
           endif
           PKPTS%flag_kgridmode=.true.
         elseif(linecount .eq. 6 .and. .not. PKPTS%flag_klinemode) then
           read(inputline,*,iostat=i_continue) PKPTS%k_shift(1:3)
           write(message,'(A,4x,3F9.5)')'  K_SHIFT:',PKPTS%k_shift(1:3)  ; write_msgi

         ! k-line if 'linemode'
         elseif(linecount .ge. 5 .and. PKPTS%flag_klinemode) then
           backspace(pid_kpoint)
           PKPTS%nline=0;i=0
    kline: do 
             read(pid_kpoint,'(A)',iostat=i_continue) inputline
             if(i_continue<0) exit               ! end of file reached
             if(i_continue>0) then
               write(message,*)'Unknown error reading file:',trim(fname),func  ; write_msgi
               kill_job
             endif
             linecount = linecount+1 ; i=i+1
             call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle kline
             call check_empty(inputline,linecount,i,flag_skip) ; if (flag_skip) cycle kline
             read(inputline,*) kline_dummy(1:3,i),k_name_dummy(i)
             if( mod(i,2) .eq. 1 .and. i .ge. 1) then
               PKPTS%nline = PKPTS%nline + 1
               write(message,'(A,I2,A,4x,3F12.8,1x,A2)')' K_LINE',PKPTS%nline,': ',kline_dummy(1:3,i),trim(k_name_dummy(i))  !; write_msgi
             elseif( mod(i,2) .eq. 0 .and. i .ge. 1 ) then
               write(message,'(2A,3F12.8,1x,A2)')trim(message), '  --> ', kline_dummy(1:3,i),trim(k_name_dummy(i))  ; write_msgi
             endif
           enddo kline
           write(message,'(A,I8)')'   N_LINE:',PKPTS%nline  ; write_msgi
           if(PKPTS%nline .ne. PKPTS%n_ndiv .and. (trim(PKPTS%kline_type) .eq. 'FHI-AIMS' .or. trim(PKPTS%kline_type) .eq. 'FHI-aims' )) then
             write(message,'(A,I0,A,I0,A)')'   !WARN! You specified ',PKPTS%nline, &
                                           ' k-path in your KFILE, but it is mismatch with the variable (#N_DIV= ', PKPTS%n_ndiv, &
                                           ') in second line of your KFILE'  ; write_msgi
             write(message,'(2A)')         '          Please check your KFILE: ',trim(fname)   ; write_msgi
             kill_job
           endif
           allocate( PKPTS%kline(3,PKPTS%nline * 2) )
           allocate( PKPTS%k_name(PKPTS%nline * 2) )
           PKPTS%kline(1:3,1:PKPTS%nline*2) = kline_dummy(1:3,1:PKPTS%nline * 2)
           PKPTS%k_name(1:PKPTS%nline*2) = k_name_dummy(1:PKPTS%nline * 2)
         endif
         if(i_continue<0) exit ! I put this line for the gfortran issue, but I'm
!        not  sure why this should be here. However, if this is not here,
!        program get some error... Check later on. 11.Dec. 2020, HJK
      enddo line

  if (linecount == 0) then
    write(message,*)'Attention - empty input file: INCAR-TB ',func  ; write_msgi
    stop
  endif
  close(pid_kpoint)

  if(PKPTS%flag_klinemode .and. .not. PKPTS%flag_kgridmode) then
     if(trim(PKPTS%kline_type) .eq. 'FLEUR' .or. trim(PKPTS%kline_type) .eq. 'fleur') then
       idiv_mode = 2 ! division type: fleur-like. n division between kpoint A and B and total n+1 points
     elseif(trim(PKPTS%kline_type) .eq. 'FHI-AIMS' .or. trim(PKPTS%kline_type) .eq. 'FHI-aims') then
       idiv_mode = 3 ! division type: vasp-like with n-1 division between each segments. In this mode, however, 
                     ! every path has different division. This is same as FHI-AIMS code does.
     else
       idiv_mode = 1 ! division type: vasp-like. n-1 division between kpoint A and B and total n points
     endif
     PKPTS%idiv_mode = idiv_mode
     call get_kpath(PKPTS, PGEOM, PKPTS%kunit, idiv_mode)
     write(message,'(A,I8)')'  NKPOINT:',PKPTS%nkpoint  ; write_msgi
     if(idiv_mode .eq. 1) then ! VASP type
       write(message,'(A,A8,*(A8))')'   K-PATH:       ', PKPTS%k_name(1), (PKPTS%k_name((ik-1)*2),ik=2,PKPTS%nline+1); write_msgi
       write(message,'(A,*(I8))')   '  (index):', (PKPTS%ndiv*(ik-1)+1,ik=1,PKPTS%nline) , PKPTS%ndiv*PKPTS%nline ; write_msgi
     elseif(idiv_mode .eq. 3) then ! AIMS type
       write(message,'(A,A8,*(A8))')'   K-PATH:       ', PKPTS%k_name(1), (PKPTS%k_name((ik-1)*2),ik=2,PKPTS%nline+1); write_msgi
       write(message,'(A,*(I8))')   '  (index):', 1,  (sum(PKPTS%ndiv(1:ik))+1,ik=1,PKPTS%nline-1) , sum(PKPTS%ndiv) ; write_msgi
     elseif(idiv_mode .eq. 2) then ! FLEUR type
       write(message,'(A,A8,*(A8))')'   K-PATH:       ', PKPTS%k_name(1), (PKPTS%k_name((ik-1)*2),ik=2,PKPTS%nline+1); write_msgi
       write(message,'(A,*(I8))')   '  (index):', 1,  (PKPTS%ndiv*ik,ik=1,PKPTS%nline) ; write_msgi
     endif
  elseif(PKPTS%flag_kgridmode .and. .not. PKPTS%flag_klinemode) then
     PKPTS%nkpoint = PKPTS%ndiv(1)*PKPTS%ndiv(2)*PKPTS%ndiv(3)
     write(message,'(A,I8)')'  NKPOINT:',PKPTS%nkpoint  ; write_msgi
     allocate( kpts_cart(3,PKPTS%nkpoint) )
     allocate( kpts_reci(3,PKPTS%nkpoint) )
     allocate( PKPTS%kpoint(3,PKPTS%nkpoint) )
     allocate( PKPTS%kpoint_reci(3,PKPTS%nkpoint) )
     call get_kgrid(kpts_cart, kpts_reci, PKPTS%ndiv(1), PKPTS%ndiv(2), PKPTS%ndiv(3), PKPTS%k_shift(1:3), PGEOM, PKPTS%flag_gamma)
     PKPTS%kpoint(:,:)=kpts_cart(:,:)
     PKPTS%kpoint_reci(:,:)=kpts_reci(:,:)
     deallocate( kpts_cart )
     deallocate( kpts_reci )
  elseif(PKPTS%flag_kgridmode .and. PKPTS%flag_klinemode) then
     write(message,'(A)')'   !WARN! Check KPOINT file. Both linemode and MP mode set simulatneously. Exit..'  ; write_msgi
     stop
  elseif( .not. PKPTS%flag_kgridmode .and. .not. PKPTS%flag_klinemode) then
     write(message,'(A)')'   !WARN! Check KPOINT file. Both linemode & MP mode does not set. Exit..'  ; write_msgi
     stop
  endif

  write(message,*)'*- END READING KPOINT FILE ---------------------'  ; write_msgi
  write(message,*)' '  ; write_msgi
return
endsubroutine
