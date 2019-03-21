#include "alias.inc"
subroutine read_kpoint(fname, PKPTS, PGEOM, PINPT)
  use parameters, only : kpoints, poscar, pid_kpoint, incar
  use mpi_setup
  implicit none
  integer*4, parameter :: max_kline=100
  integer*4     i_continue, nitems,ndiv_temp
  integer*4     i,linecount, i_dummy
  real*8        kline_dummy(3,max_kline)
  real*8, allocatable :: kpts_cart(:,:), kpts_reci(:,:)
  character*132 inputline,fname
  character*40  desc_str,dummy, k_name_dummy(max_kline)
  character(*), parameter :: func = 'read_kpoint'
  logical       flag_skip
  external      nitems
  type(kpoints) :: PKPTS
  type(poscar)  :: PGEOM
  type(incar )  :: PINPT
  PKPTS%flag_cartesianK = .false.
  PKPTS%flag_reciprocal = .false.
  PKPTS%flag_kgridmode  = .false.
  PKPTS%flag_klinemode  = .false.

  if_main write(6,*)' '
  if_main write(6,*)'*- READING KPOINTS FILE: ',trim(fname)
  open (pid_kpoint, FILE=fname,iostat=i_continue)
  linecount = 0 

 line: do
        read(pid_kpoint,'(A)',iostat=i_continue) inputline
        if(i_continue<0) exit               ! end of file reached
        if(i_continue>0) then
          if_main write(6,*)'Unknown error reading file:',trim(fname),func
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
           read(inputline,*,iostat=i_continue) ndiv_temp
           cycle

        ! k-grid mode
         elseif(linecount .eq. 3) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'L' .or. desc_str(1:1) .eq. 'l') then
             if_main write(6,'(A)')'   K_MODE: Line-mode'
             PKPTS%flag_klinemode=.true.
             if(.not. PINPT%flag_ndiv_line_parse) then
               if(.not. allocated(PKPTS%ndiv)) allocate( PKPTS%ndiv(1) )
               PKPTS%ndiv(1) = ndiv_temp
               if_main write(6,'(A,I8)')'    N_DIV:',PKPTS%ndiv
             else
               if_main write(6,'(A,I8)')'    N_DIV: (set by -nkp_line option) ',PKPTS%ndiv
             endif
           elseif(desc_str(1:1) .eq. 'G' .or. desc_str(1:1) .eq. 'g') then
             linecount = linecount + 1
             if_main write(6,'(A)')'   K_MODE: Gamma-centered'
             PKPTS%flag_gamma=.true.
             if(.not. allocated(PKPTS%ndiv)) allocate( PKPTS%ndiv(3) )
           elseif(desc_str(1:1) .eq. 'M' .or. desc_str(1:1) .eq. 'm') then
             linecount = linecount + 1
             if_main write(6,'(A)')'   K_MODE: non Gamma-centered'
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
             if_main write(6,'(A)')'   K_TYPE: Reciprocal unit'
           elseif( (desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c'  .or.  &
                    desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') .and. &
                    PKPTS%flag_klinemode ) then
             PKPTS%flag_reciprocal=.false.
             PKPTS%flag_cartesianK=.true.
             if_main write(6,'(A)')'   K_TYPE: Cartesian unit (1/A)'
           endif

         ! k-grid if .not. 'linemode' .and. 'kgridmode)
         elseif(linecount .eq. 5 .and. .not. PKPTS%flag_klinemode) then
           if(.not. PINPT%flag_ndiv_grid_parse) then
             read(inputline,*,iostat=i_continue) PKPTS%ndiv(1:3)
             if_main write(6,'(A,4x,3I4)')'   K_GRID:',PKPTS%ndiv(1:3)
           else
             read(inputline,*,iostat=i_continue) i_dummy, i_dummy, i_dummy
             if_main write(6,'(A,4x,3I4)')'   K_GRID: (set by -nkp_grid option)',PKPTS%ndiv(1:3)
           endif
           PKPTS%flag_kgridmode=.true.
         elseif(linecount .eq. 6 .and. .not. PKPTS%flag_klinemode) then
           read(inputline,*,iostat=i_continue) PKPTS%k_shift(1:3)
           if(myid .eq. 0) write(6,'(A,4x,3F9.5)')'  K_SHIFT:',PKPTS%k_shift(1:3)

         ! k-line if 'linemode'
         elseif(linecount .ge. 5 .and. PKPTS%flag_klinemode) then
           backspace(pid_kpoint)
           PKPTS%nline=0;i=0
    kline: do 
             read(pid_kpoint,'(A)',iostat=i_continue) inputline
             if(i_continue<0) exit               ! end of file reached
             if(i_continue>0) then
               if(myid .eq. 0) write(6,*)'Unknown error reading file:',trim(fname),func
             endif
             linecount = linecount+1 ; i=i+1
             call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle kline
             call check_empty(inputline,linecount,i,flag_skip) ; if (flag_skip) cycle kline
             read(inputline,*,iostat=i_continue) kline_dummy(1:3,i),k_name_dummy(i)
             if( mod(i,2) .eq. 1 .and. i .ge. 1) then
               PKPTS%nline = PKPTS%nline + 1
               if(myid .eq. 0) write(6,'(A,I2,A,4x,3F12.8,1x,A2)',ADVANCE='NO')' K_LINE',PKPTS%nline,': ',kline_dummy(1:3,i),trim(k_name_dummy(i))
             elseif( mod(i,2) .eq. 0 .and. i .ge. 1 ) then
               if(myid .eq. 0) write(6,'(A,3F12.8,1x,A2)')'  --> ', kline_dummy(1:3,i),trim(k_name_dummy(i))
             endif
           enddo kline
           if(myid .eq. 0) write(6,'(A,I8)')'   N_LINE:',PKPTS%nline
           allocate( PKPTS%kline(3,PKPTS%nline * 2) )
           allocate( PKPTS%k_name(PKPTS%nline * 2) )
           PKPTS%kline(1:3,1:PKPTS%nline*2) = kline_dummy(1:3,1:PKPTS%nline * 2)
           PKPTS%k_name(1:PKPTS%nline*2) = k_name_dummy(1:PKPTS%nline * 2)
         endif

      enddo line

  if (linecount == 0) then
    if(myid .eq. 0) write(6,*)'Attention - empty input file: INCAR-TB ',func
    stop
  endif
  close(pid_kpoint)

  if(PKPTS%flag_klinemode .and. .not. PKPTS%flag_kgridmode) then
     call get_kpath(PKPTS, PGEOM, PKPTS%kunit)
     if(myid .eq. 0) write(6,'(A,I8)')'  NKPOINT:',PKPTS%nkpoint
  elseif(PKPTS%flag_kgridmode .and. .not. PKPTS%flag_klinemode) then
     PKPTS%nkpoint = PKPTS%ndiv(1)*PKPTS%ndiv(2)*PKPTS%ndiv(3)
     if(myid .eq. 0) write(6,'(A,I8)')'  NKPOINT:',PKPTS%nkpoint
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
     if(myid .eq. 0) write(6,'(A)')'   !WARN! Check KPOINT file. Both linemode & MP mode set simulatneously. Exit..'
     stop
  elseif( .not. PKPTS%flag_kgridmode .and. .not. PKPTS%flag_klinemode) then
     if(myid .eq. 0) write(6,'(A)')'   !WARN! Check KPOINT file. Both linemode & MP mode does not set. Exit..'
     stop
  endif

  if(myid .eq. 0) write(6,*)'*- END READING KPOINT FILE ---------------------'
  if(myid .eq. 0) write(6,*)' '
return
endsubroutine
