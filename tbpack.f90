#include "alias.inc"
subroutine get_kpath(PKPTS, PGEOM, kunit, idiv_mode)
  use parameters, only : kpoints, poscar
  use mpi_setup
  use print_io
  implicit none
  type(kpoints) :: PKPTS
  type(poscar)  :: PGEOM
  real*8           a1(3),a2(3),a3(3), b1(3),b2(3),b3(3)
  real*8           dk(3),dk_(3), dk_reci(3)
  real*8           kdist, enorm
  real*8, allocatable :: PK(:,:), PK_reci(:,:)
  integer*4        i,ii,ik,iline,ndivk, iline_
  integer*4        idiv_mode
  character*1      kunit
  external         enorm
  a1=PGEOM%a_latt(1:3,1)
  a2=PGEOM%a_latt(1:3,2)
  a3=PGEOM%a_latt(1:3,3)
  call get_reci(b1,b2,b3, a1,a2,a3)

  ! k-line mode
  if(PKPTS%flag_klinemode) then
    if(idiv_mode .eq. 1) then ! VASP type
      ndivk=PKPTS%ndiv(1)
      PKPTS%nkpoint = ndivk * PKPTS%nline
    elseif(idiv_mode .eq. 2) then ! FLEUR type
      ndivk=PKPTS%ndiv(1)
      PKPTS%nkpoint = ndivk * PKPTS%nline + 1
    elseif(idiv_mode .eq. 3) then ! AIMS type
      PKPTS%nkpoint = sum( PKPTS%ndiv(:) )
    endif
    allocate( PK(3,PKPTS%nline*2) )
    allocate( PK_reci(3,PKPTS%nline*2) )
    allocate( PKPTS%kpoint(3,PKPTS%nkpoint) )
    allocate( PKPTS%kpoint_reci(3,PKPTS%nkpoint) )
    
    do iline=1,PKPTS%nline*2
      PK(1:3,iline)=PKPTS%kline(1,iline)*b1(1:3)+ &
                    PKPTS%kline(2,iline)*b2(1:3)+ &
                    PKPTS%kline(3,iline)*b3(1:3)
      PK_reci(1:3,iline)= PKPTS%kline(1:3,iline)
    enddo

    ii=0
    kdist = 0.0
    if(trim(PKPTS%k_name(1)) .eq. "G" .or. trim(PKPTS%k_name(1)) .eq. 'g') then
      write(message,'(A,A,4x,F10.6,2A,4A)')'    KINIT','= ',kdist,' ; ','KNAME_INIT','=','"','{/Symbol G}','"'  ; write_msg
    else
      write(message,'(A,A,4x,F10.6,2A,4A)')'    KINIT','= ',kdist,' ; ','KNAME_INIT','=','"',trim(PKPTS%k_name(1)),'"'  ; write_msg
    endif
    iline_ = 0
    do iline=1,PKPTS%nline*2, 2
      if(idiv_mode .eq. 1) then
        dk_(:)= ( PKPTS%kline(:,iline + 1) - PKPTS%kline(:,iline) ) / real(ndivk-1)
      elseif(idiv_mode .eq. 2) then
        dk_(:)= ( PKPTS%kline(:,iline + 1) - PKPTS%kline(:,iline) ) / real(ndivk)
      elseif(idiv_mode .eq. 3) then
        iline_ = iline_ + 1
        ndivk = PKPTS%ndiv(iline_)
        dk_(:)= ( PKPTS%kline(:,iline + 1) - PKPTS%kline(:,iline) ) / real(ndivk-1)
      endif
      dk(1:3)=dk_(1) * b1(1:3) + dk_(2) * b2(1:3) + dk_(3) * b3(1:3)
      dk_reci(1:3)=dk_(1:3)

      if(iline .ge. 1 .and. iline .ne. PKPTS%nline*2-1) then
        if(idiv_mode .eq. 1 .or. idiv_mode .eq. 3) then
          kdist= kdist + enorm(3,dk) * real(ndivk-1)
        elseif(idiv_mode .eq. 2) then
          kdist= kdist + enorm(3,dk) * real(ndivk)
        endif
        if( (iline-1)/2 + 2 .lt. 10) then
          if(trim(PKPTS%k_name(iline+1)) .eq. "G" .or. trim(PKPTS%k_name(iline+1)) .eq. 'g') then
            write(message,'(A,I1,A,4x,F10.6,2A,I1,4A)')'       K',(iline-1)/2 + 2,'= ',kdist, ' ; ','KNAME_',(iline-1)/2 + 2,'   =','"','{/Symbol G}','"'  ; write_msg
          else
            write(message,'(A,I1,A,4x,F10.6,2A,I1,4A)')'       K',(iline-1)/2 + 2,'= ',kdist, ' ; ','KNAME_',(iline-1)/2 + 2,'   =','"',trim(PKPTS%k_name(iline+1)),'"'  ; write_msg
          endif
        elseif(iline .lt. 100) then
          if(trim(PKPTS%k_name(iline+1)) .eq. "G" .or. trim(PKPTS%k_name(iline+1)) .eq. 'g') then
            write(message,'(A,I2,A,4x,F10.6,2A,I2,4A)')'      K',(iline-1)/2 + 2,'= ',kdist, ' ; ','KNAME_',(iline-1)/2 + 2,'  =','"','{/Symbol G}','"'  ; write_msg
          else
            write(message,'(A,I2,A,4x,F10.6,2A,I2,4A)')'      K',(iline-1)/2 + 2,'= ',kdist, ' ; ','KNAME_',(iline-1)/2 + 2,'  =','"',trim(PKPTS%k_name(iline+1)),'"'  ; write_msg
          endif
        endif
      elseif(iline .ge. 1 .and. iline .eq. PKPTS%nline*2-1) then
        if(idiv_mode .eq. 1 .or. idiv_mode .eq. 3) then
          kdist=kdist + enorm(3,dk) * real(ndivk-1)
        elseif(idiv_mode .eq. 2) then
          kdist=kdist + enorm(3,dk) * real(ndivk)
        endif
        if(trim(PKPTS%k_name(iline+1)) .eq. "G" .or. trim(PKPTS%k_name(iline+1)) .eq. 'g') then
          write(message,'(A,A,4x,F10.6,2A,4A)')'     KEND','= ',kdist,' ; ','KNAME_END',' =','"','{/Symbol G}','"'  ; write_msg
        else
          write(message,'(A,A,4x,F10.6,2A,4A)')'     KEND','= ',kdist,' ; ','KNAME_END',' =','"',trim(PKPTS%k_name(iline+1)),'"'  ; write_msg
        endif
      endif

      do ik=1,ndivk
        ii=ii + 1
        PKPTS%kpoint(1:3, ii ) = PK(1:3,iline) + dk(1:3)*(ik-1)
        PKPTS%kpoint_reci(1:3, ii ) = PK_reci(1:3,iline) + dk_reci(1:3)*(ik-1)
      enddo    
      if (idiv_mode .eq. 2 .and. iline .eq. PKPTS%nline*2-1) then
        ii=ii + 1
        PKPTS%kpoint(1:3, ii ) = PK(1:3,iline) + dk(1:3)*ndivk
        PKPTS%kpoint_reci(1:3, ii ) = PK_reci(1:3,iline) + dk_reci(1:3)*ndivk
      endif
    enddo

    do iline=1, PKPTS%nline
      if(iline .eq. 1) then
        write(message,'(A)')' set xtics (KNAME_INIT KINIT,' 
      elseif(iline .ge. 2 .and. iline .lt. PKPTS%nline) then
        write(message,'(2A,i0,A,i0,A)')trim(message),' KNAME_',iline,' K',iline,','
      endif

      if(iline .eq. PKPTS%nline .and. iline .eq. 1) then
        write(message,'(2A)')trim(message), ' KNAME_END KEND) nomirror'
      elseif(iline .eq. PKPTS%nline .and. iline .gt. 1) then
        write(message,'(2A,i0,A,i0,A)')trim(message),' KNAME_',iline,' K',iline,', KNAME_END KEND) nomirror'
      endif
    enddo
    write_msg

  ! k-grid mode
  elseif(PKPTS%flag_kgridmode) then

    call get_kgrid(PKPTS%kpoint, PKPTS%kpoint_reci, PKPTS%ndiv(1), PKPTS%ndiv(2), PKPTS%ndiv(3),&
                   PKPTS%k_shift, PGEOM, PKPTS%flag_gamma)

  endif

return
endsubroutine
subroutine get_kline(kp,kp_reci,nk,PGEOM,kinit_reci,kend_reci)
   use parameters, only : poscar
   implicit none
   type(poscar) :: PGEOM
   integer*4       nk
   integer*4       ik
   real*8          kinit(3), kend(3)
   real*8          kinit_reci(3), kend_reci(3)
   real*8          kp(3,nk), kp_reci(3,nk) 
   real*8          dk(3), dk_reci(3)
   real*8          a1(3),a2(3),a3(3), b1(3),b2(3),b3(3)

   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   call get_reci(b1,b2,b3, a1,a2,a3)

   kinit = kinit_reci(1) * b1(:) + kinit_reci(2) * b2(:) + kinit_reci(3) * b3(:)
   kend  =  kend_reci(1) * b1(:) +  kend_reci(2) * b2(:) +  kend_reci(3) * b3(:)
   dk_reci(:)= ( kend_reci - kinit_reci ) / (nk-1)
   dk(:)     =dk_reci(1) * b1(1:3) + dk_reci(2) * b2(1:3) + dk_reci(3) * b3(1:3)

   do ik=1,nk
     kp(1:3, ik)      = kinit(1:3)      + dk(1:3)     *(ik-1)
     kp_reci(1:3, ik) = kinit_reci(1:3) + dk_reci(1:3)*(ik-1)
   enddo

   return
endsubroutine
subroutine get_kpoint(kp,kp_,nk,PGEOM)
   use parameters, only : poscar
   type(poscar)  :: PGEOM
   integer*4        nk, ik
   real*8           kp(3,nk), kp_(3,nk) ! kp_=reciprocal unit
   real*8           a1(3),a2(3),a3(3), b1(3), b2(3), b3(3)
  
   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   call get_reci(b1,b2,b3, a1, a2, a3)
   
   do ik = 1, nk
     kp(:,ik) =  b1(:)*kp_(1,ik) + b2(:)*kp_(2,ik) + b3(:)*kp_(3,ik)
   enddo

   return
endsubroutine
subroutine get_kgrid(kp,kp_,nk1,nk2,nk3,kshift,PGEOM, flag_gamma)
  use parameters, only : poscar
  implicit none
  integer*4     nk1,nk2,nk3, ii, ik1,ik2,ik3
  real*8        kp(3,nk1*nk2*nk3), kp_(3,nk1*nk2*nk3)  ! kp_ reciprocal unit
  real*8        a1(3),a2(3),a3(3), b1(3),b2(3),b3(3)
  real*8        kshift(3),r_dummy1,r_dummy2,r_dummy3
  logical       flag_gamma
  type(poscar)  :: PGEOM

  a1=PGEOM%a_latt(1:3,1)
  a2=PGEOM%a_latt(1:3,2)
  a3=PGEOM%a_latt(1:3,3)
  call get_reci(b1,b2,b3, a1,a2,a3)
  if(.not. flag_gamma) then 
    kshift = kshift + (/0.5d0, 0.5d0, 0.5d0/)
  endif

  ii = 0
  do ik3 = 0, nk3-1
    do ik2 = 0, nk2-1
      do ik1 = 0, nk1-1
        ii = ii + 1
        r_dummy1 = (dble(ik1)+ kshift(1)) / dble(nk1)
        r_dummy2 = (dble(ik2)+ kshift(2)) / dble(nk2)
        r_dummy3 = (dble(ik3)+ kshift(3)) / dble(nk3)

        if(r_dummy1 .gt. 0.5d0+tiny(0.5d0)) then 
          kp_(1,ii) = r_dummy1- 1.0d0
        else
          kp_(1,ii) = r_dummy1
        endif

        if(r_dummy2 .gt. 0.5d0+tiny(0.5d0)) then 
          kp_(2,ii) = r_dummy2- 1.0d0
        else
          kp_(2,ii) = r_dummy2
        endif

        if(r_dummy3 .gt. 0.5d0+tiny(0.5d0)) then 
          kp_(3,ii) = r_dummy3- 1.0d0
        else
          kp_(3,ii) = r_dummy3
        endif

        kp(:,ii) =  b1(:)*kp_(1,ii) + b2(:)*kp_(2,ii) + b3(:)*kp_(3,ii) 

      enddo
    enddo
  enddo
  
return
endsubroutine
subroutine print_kpoint (kpoint_reci, nkpoint, fname)
  use parameters, only : poscar, pid_ibzkpt
  implicit none
  integer*4               i, nkpoint
  real*8                  kpoint_reci(3,nkpoint) 
  character(*) fname
  type(poscar) :: PGEOM

  open(pid_ibzkpt, file = trim(fname) , status = 'unknown')
  write(pid_ibzkpt,'(A)')'Automatically generated mesh'
  write(pid_ibzkpt,'(I8)') nkpoint
  write(pid_ibzkpt,'(A)')'Reciprocal lattice'
  do i = 1, nkpoint
    write(pid_ibzkpt,'(3F20.14,I6)') kpoint_reci(:,i), 1
  enddo

  close(pid_ibzkpt)

return
endsubroutine

! subroutine for computing reciprocal lattice vector
subroutine get_reci(b1,b2,b3, a1,a2,a3)
  use parameters, only : pi
  implicit real*8 (a-h,o-z)
  real*8 a1(3),a2(3),a3(3), a2xa3(3)
  real*8 b1(3),b2(3),b3(3), b1xb2(3)

  call vcross(a2xa3,a2,a3)
  Vcell=dot_product(a1,a2xa3)

  call vcross(b1,a2,a3)
  call vcross(b2,a3,a1)
  call vcross(b3,a1,a2)

  b1=2.d0*pi*b1/Vcell
  b2=2.d0*pi*b2/Vcell
  b3=2.d0*pi*b3/Vcell
return
endsubroutine get_reci
subroutine vcross(a,b,c)
  ! subroutine for computing vector cross-product
  implicit none
  real*8 a(3),b(3),c(3)

  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
return
end subroutine vcross
subroutine check_here_c(check_flag)
  implicit none
  character(*) check_flag

  write(6,*)'  !!**!! ',trim(check_flag)

return
endsubroutine
subroutine check_here_i(check_flag)
  implicit none
  integer*4   check_flag

  write(6,*)'  !!**!! ',check_flag

return
endsubroutine
subroutine check_here_r(check_flag)
  implicit none
  real*8   check_flag

  write(6,*)'  !!**!! ',check_flag

return
endsubroutine
subroutine get_window(init,fina,inputline,desc_str)
  use mpi_setup
  use print_io
  implicit none
  integer*4      nitems, i_dummy, mpierr
  external       nitems
  character*132  inputline, dummy, dummy_
  character*40   desc_str, dummy1, dummy2
  real*8         init,fina

  call strip_off(inputline, dummy, ' ', '#',0) ! cut off unnecessary comments
! call strip_off(dummy, dummy_, ' ', '!', 0) ! cut off unecessary comments   
! dummy_ = dummy 
  if(index(dummy, trim(desc_str)) .ge. 1) then
    inputline = dummy
  endif

  call strip_off(inputline, dummy, trim(desc_str), ' ', 2) ! cut off description strip

  if(index(dummy, ':') .gt. 1) then
    call strip_off(dummy, dummy1, ' ', ':', 0) ! cut off init
    call str2real(dummy1, init)
    call strip_off(dummy, dummy1, ':', ' ', 2) ! cut off fina
    call str2real(dummy1, fina)
  elseif(index(dummy,':') .eq. 0) then
    if(nitems(dummy) .ne. 2) then
      write(message,'(A)')'    !WARN! Number of items to be read is not equal to 2 in the EWINDOW'  ; write_msg
      write(message,'(A)')'           The correct syntax is -> INIT_E:FINA_E (ex, -2:2) '  ; write_msg
      write(message,'(A)')'           or                       INIT_E FINA_E (ex, -2 2) '  ; write_msg
      write(message,'(A)')'           Exit program...'  ; write_msg
      kill_job
    elseif(nitems(dummy) .eq. 2) then
      read(dummy,*) init, fina
    endif
  elseif(index(dummy,':') .eq. 1) then
      write(message,'(A)')'    !WARN! INIT_ENERGY and FINAL_ENERGY has not been declaired properly.'  ; write_msg
      write(message,'(A)')'           The correct syntax is -> INIT_E:FINA_E (ex, -2:2) '  ; write_msg
      write(message,'(A)')'           Exit program...'  ; write_msg
      kill_job
  endif

return
endsubroutine
subroutine take_comment(string, strip)
! comment off from the "string" and return as "strip".
! take out comment mark only and return remaining part
   implicit none
   logical blank
   character(*) string, strip
   integer*4    l0, init

   init = -1
   strip = ''
   l0 = len_trim(string)

   init = index(string,'#',.FALSE.)
   if(init .eq. 0) then
     strip = string
   elseif(init .ge. 1) then
     strip = string(init+1:l0)
   endif

return
endsubroutine
subroutine comment_off(string, strip) 
! comment off from the "string" and return as "strip".
! take out comment mark only and return remaining part
   implicit none
   logical blank
   character(*) string, strip
   integer*4    l0, init
    
   init = -1
   strip = ''
   l0 = len_trim(string)

   init = index(string,'#',.FALSE.)
   if(init .eq. 0) then
     strip = string
   elseif(init .ge. 1) then
     strip = string(1:init-1)
   endif

return 
endsubroutine
subroutine strip_off (string, strip, strip_a, strip_b, mode)
!strip   : strip to be extract out of string
!strip_a : strip_index a   
!strip_b : strip_index b
!ex) string = hello world ! -> strip_a='hello', strip_b='!', strip =' world ' and mode = 1

!mode 0: strip-off where              strip   < strip_b of string 
!mode 1: strip-off where   strip_a <  strip   < strip_b of string only if strip_a =/ strip b
!mode 2: strip-off where   strip_a <  strip 
!mode 3: strip-off where   strip_a <  strip   < strip_b of string only if strip_a = strip_b 
!mode 4: strip-off where   strip_a <  strip_b < strip   of string only if strip_a = strip_b 
!mode 5: same sa mode 2 but ignore '#' tag.
  implicit none
  logical blank
  character(*) string, strip_a, strip_b, strip
  integer*4 mode,mode_check,l0,la,lb,ls, i, j, ii, init,fini
  strip=''

  l0 = len_trim(string)
  mode_check= mode
  la=len_trim(strip_a)
  lb=len_trim(strip_b)
  if (la .eq. lb .and. mode .eq. 1) then
    mode_check = 3
  else
    mode_check = mode
  endif

  if(mode_check .eq. 0) then
    do i = 1, l0
      if(string(i:i+lb-1) .eq. trim(strip_b))then
        strip=adjustl(trim(string(1:i-1)))
        exit
      endif
    enddo

  elseif(mode_check .eq. 1) then
    do i = 1, l0
      if(string(i:i+la-1) .eq. trim(strip_a)) then 
        init=i+la
        do j = init + 1 , l0
          if(string(j:j+lb-1) .eq. trim(strip_b)) then
            fini=j-1
            strip=adjustl(trim(string(init:fini)))
            exit
          endif
        enddo
      endif
!     if(string(i:i+lb-1) .eq. trim(strip_b)) then
!       fini=i-1
!       strip=adjustl(trim(string(init:fini)))
!     endif
    enddo

  elseif(mode_check .eq. 2) then
    do i = 1, l0
      if(string(i:i) .eq. '#') then
        fini=i - 1
        exit
      endif
      fini = l0
    enddo
    do i = 1, l0
      if(string(i:i+la-1) .eq. trim(strip_a)) then
        init=i+la
        strip=adjustl(trim(string(init:fini)))
        exit
      endif
    enddo

  elseif(mode_check .eq. 3) then
    ii = 0
    do i = 1, l0
      ii = ii + 1
      if(string(ii:ii+la-1) .eq. trim(strip_a)) then 
        init=ii+la
        ii = ii + la
      endif
      if(string(ii:ii+lb-1) .eq. trim(strip_b)) then
        fini=ii-1
        strip=adjustl(trim(string(init:fini)))
      endif
    enddo

  elseif(mode_check .eq. 4) then
    init = index(string,':',.TRUE.) + 1
    do i = ii+1, l0
      if(string(i:i) .eq. '#') then
        fini=i - 1
        exit
      endif
      fini = l0
    enddo
    strip=adjustl(trim(string(init:fini)))
  elseif(mode_check .eq. 5) then
    do i = 1, l0
      fini = l0
    enddo
    do i = 1, l0
      if(string(i:i+la-1) .eq. trim(strip_a)) then
        init=i+la
        strip=adjustl(trim(string(init:fini)))
        exit
      endif
    enddo

  endif
return
endsubroutine
subroutine check_comment(inputline,linecount,i,flag_skip)
  implicit none
  integer*4 i,linecount,i_continue
  character(*)inputline
  character*40 desc_str
  logical flag_skip
  read(inputline,*,iostat=i_continue) desc_str
  if (linecount .ne. 1 .and. desc_str(1:1).eq.'#') then
    linecount=linecount - 1
    flag_skip = .true.
    i=i - 1
  else
    flag_skip = .false.
    i=i
  endif

return
endsubroutine
subroutine check_empty(inputline,linecount,i,flag_skip)
  implicit none
  integer*4 i,linecount,i_continue
  character(*)inputline
  character*40 desc_str
  logical flag_skip

  read(inputline,*,iostat=i_continue) desc_str
  if(i_continue .ne. 0) then
   flag_skip = .true.
   i = i - 1
   linecount = linecount - 1
  else
   flag_skip = .false.
   i = i
  endif

return
endsubroutine
! check number of items in "string" except comments # or !
function nitems(string)
  implicit none
  logical blank
  integer*4 nitems,l,i
  character(*),intent(in) :: string
  nitems=0
  l=len_trim(string)
  blank = .true.
  do i=1,l
   if(string(i:i) .eq. '#' .or. string(i:i) .eq. '!') exit

   if (blank .and. string(i:i) .ne. ' ' ) then
     blank=.false.
     nitems=nitems + 1
   elseif( .not. blank .and. string(i:i) .eq. ' ') then
     blank=.true.
   endif
  enddo
  return
endfunction
subroutine str2logical(string,flag_logical,flag)
  use mpi_setup
  implicit none
  character(*),intent(in) :: string
  logical, intent(out) :: flag, flag_logical
  integer*4 iflag

  read(string, *, iostat=iflag) flag

  if(iflag .ge. 1) then
    flag_logical = .false.
    flag = .true.
  elseif(iflag .eq. 0) then
    flag_logical = .true.
    flag = .true.
  elseif(iflag .lt. 0) then
    flag_logical = .true.
    flag = .true.
  endif

  return
endsubroutine
function int2str(w) result(string)
  implicit none
! character(*), intent(out) :: string
  character*20  string
  integer*4,    intent(in)  :: w

  write(string,*) w

  return
endfunction
!subroutine int2str(w,string)
! implicit none
! character(*), intent(out) :: string
! integer*4,    intent(in)  :: w

! write(string,*) w

! return
!endsubroutine
subroutine str2int(string,w)
  implicit none
  character(*),intent(in) :: string
  integer*4,intent(out)   :: w

  read(string,*) w

return
end subroutine
subroutine str2real(string,w)
  use mpi_setup
  use print_io
  implicit none
  character(*),intent(in) :: string
  real*8,intent(out)   :: w
  integer*4 i_number, mpierr

  read(string,*,iostat=i_number) w

  if(i_number .ne. 0) then
    write(message,'(A)')' !!!! Error !!! "string" cannot be stored in "real" value. str2real. Exit.'  ; write_msg
    kill_job
  endif

return
end subroutine
function flag_number(string)
  implicit none
  character(*), intent(in) :: string
  logical    flag_number
  real*8     x
  integer*4  i_number

  read(string,*,iostat=i_number) x
  flag_number = i_number == 0

endfunction
function enorm ( n, x )
  implicit none

  integer*4 n
  real*8 x(n),enorm

 !enorm = sqrt ( sum ( x(1:n) ** 2 ))
  enorm = sqrt ( sum ( x(1:n) ** 2 ))
! enorm = abs  ( sum ( abs (x(1:n)) ))
  return
end
subroutine kronprod(A,B,row_a,col_a,row_b,col_b,AB) 
 implicit none
 integer*4, intent(in)   :: row_a,col_a, row_b,col_b
 complex*16, intent(in)  :: A(row_a,col_a)
 complex*16, intent(in)  :: B(row_b,col_b)
 complex*16, intent(out) :: AB(row_a*row_b,col_a*col_b)
 integer i, j , k , l 
 integer m, n , p , q 

 AB = 0d0
 do i = 1,row_a
  do j = 1,col_a
   n=(i-1)*row_b + 1
   m=n+row_b -1
   p=(j-1)*col_b + 1
   q=p+col_b -1
   AB(n:m,p:q) = A(i,j)*B
  enddo
 enddo

return
endsubroutine
subroutine show(pid,A,row_a,col_a) !Write array A in row,column order.
 integer*4  pid 
 integer*4  i_row, row_a, col_a
 real*8     A(row_a,col_a) 

  do i_row = 1,row_a
    write (pid,'(*(F10.4))') A(i_row,1:col_a) 
  enddo

return
endsubroutine
subroutine show_c(pid,A,row_a,col_a) !Write array A in row,column order.
 integer*4  pid
 integer*4  i_row, row_a, col_a
 complex*16 A(row_a,col_a)

  do i_row = 1,row_a
    write (pid,'(*(F10.4,F10.4))') A(i_row,1:col_a)
  enddo

return
endsubroutine
!Below the StrUpCase and StrLowCase is adopted from 'computer-programming-forum.com' and the author is Paul van Delst. 
!And it is slightly modified for the purpose...
FUNCTION str2upcase ( string )
  CHARACTER(*), INTENT(IN) :: string 
  CHARACTER(LEN(string))   :: str2upcase
  INTEGER :: i, n 
  CHARACTER*26  LOWER_CASE, UPPER_CASE
  LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  str2upcase = string
  DO i = 1, LEN( string ) 
    n = INDEX( LOWER_CASE, string(i:i) ) 
    IF ( n /= 0 ) str2upcase(i:i) = UPPER_CASE(n:n) 
  END DO 
END FUNCTION str2upcase 
FUNCTION str2lowcase( string )
  CHARACTER(*), INTENT(IN) :: string
  CHARACTER(LEN(string))   :: str2lowcase
  INTEGER :: i, n 
  CHARACTER*26  LOWER_CASE, UPPER_CASE
  LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  str2lowcase = string
  DO i = 1, LEN(string) 
    n = INDEX( UPPER_CASE, string(i:i) ) 
    IF ( n /= 0 ) str2lowcase(i:i) = LOWER_CASE(n:n) 
  END DO 
END FUNCTION str2lowcase

! return index of original 2D array from reshaped 1D array index from 2D array
!function idx2Di(idx1D, ni) result(i)
!   implicit none
!   integer*4    i, ni
!   integer*4    idx1D
!   
!   i = floor( real(idx1D)/real(ni) ) + 1
!!  j = mod( idx1D, ni )
!return
!endfunction
!function idx2Dj(idx1D, ni) result(j)
!  implicit none
!  integer*4    j, ni 
!  integer*4    idx1D

!!  i = floor( real(idx1D)/real(ni) ) + 1
!  j = mod( idx1D, ni )
!return
!endfunction


