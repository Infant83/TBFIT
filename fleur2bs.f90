! conversion program for "struct.xsf" file of FLEUR into "struct.bs" format for xbs.
! by Hyun-Jung Kim. PGI/IAS, Forschungszentrum Juelich, h.kim@fz-juelich.de
! last update : 2020.08.04 by KHJ
! last update : 2020.12.07 by KHJ

program fleur2bs
  implicit none
  integer,parameter::matom=129024
  integer,parameter::mtype=16
  character(len=80)::title,fname,atom_name, dummy1, dummy2, title_
  integer::nx,ny,nz, ishift
  real*8   shift_nx,shift_ny,shift_nz
  character(len=132)::line, desc_str
  integer::lline,column,tend,i,j,k,ibound
  real(8),dimension(3,3)::lattice,r_lattice
  real(8)::vol
  integer::ntype, natom_total
  integer,dimension(mtype)::natom
  character(len=40)::TCHAR_
  character(len=1) ::TCHAR
  logical::selective,cartesian
  integer::itype,iatom,iclass
  real(8),dimension(3,matom,mtype)::coord
  logical,dimension(3,matom,mtype)::movable
  real(8),dimension(3)::c,tvec
  logical,dimension(3)::m

  character(len=5),dimension(mtype)::aname
  real(8),dimension(0:mtype)::radius
  real(8),dimension(0:mtype)::radi_R
  real(8),dimension(3,mtype)::color
  real(8),dimension(3)::bcolor
  integer::idum, pid_atoms, pid_geom
  real(8)::x_min,x_max,y_min,y_max,z_min,z_max,c_min,c_max
 
  logical              flag_done 
  logical, external :: is_string
  real(8)           :: bcolor_type, radi_R_min, radi_RR
  integer, external :: nitems

  ishift = 0
  pid_atoms = 77
  pid_geom  = 78
  flag_done = .false.
  bcolor_type=0.0 ! default
  ntype = 0
  ! read atoms_fleur.d header info.
  open(pid_atoms,file='atoms.d',status='old',form='formatted')
  rewind(pid_atoms)
  read(pid_atoms,*)fname
  read(pid_atoms,*)nx,ny,nz, desc_str
  if (desc_str(1:1).ne.'#') then
   rewind(pid_atoms)
   read(pid_atoms,*)fname
   read(pid_atoms,*)nx,ny,nz, shift_nx,shift_ny,shift_nz
   ishift = 1
  endif
! read(pid_atoms,*)ntype
! read(pid_atoms,*)natom(1:ntype)

  
  do while (.not. flag_done)
    read(pid_atoms,*)dummy1, dummy2
    if(is_string(trim(dummy2))) then 
      ntype = ntype + 1
      backspace(pid_atoms)
      read(pid_atoms,*)natom(ntype),aname(ntype),radius(ntype),color(:,ntype),radi_R(ntype)
      flag_done = .false.
    else
      flag_done = .true.
      backspace(pid_atoms)
    endif
  enddo
! do i=1,ntype
!  write(6,*)"ZZZ ", natom(i), aname(i),radius(i),color(:,i),radi_R(i)
! end do
  
! stop
  ! read POSCAR file
  open(pid_geom,file=fname,status='old',form='formatted')
  rewind(pid_geom)
  ! title
  read(pid_geom,'(a)')title
  write(0,'(a)')title(1:len_trim(title))
  read(pid_geom,'(a)')title
  write(0,'(a)')title(1:len_trim(title))
  ! lattice vectors
  do i=1,3; read(pid_geom,*)lattice(1:3,i); end do
  r_lattice(1,1)=lattice(2,2)*lattice(3,3)-lattice(3,2)*lattice(2,3)
  r_lattice(2,1)=lattice(3,2)*lattice(1,3)-lattice(1,2)*lattice(3,3)
  r_lattice(3,1)=lattice(1,2)*lattice(2,3)-lattice(2,2)*lattice(1,3)
  r_lattice(1,2)=lattice(2,3)*lattice(3,1)-lattice(3,3)*lattice(2,1)
  r_lattice(2,2)=lattice(3,3)*lattice(1,1)-lattice(1,3)*lattice(3,1)
  r_lattice(3,2)=lattice(1,3)*lattice(2,1)-lattice(2,3)*lattice(1,1)
  r_lattice(1,3)=lattice(2,1)*lattice(3,2)-lattice(3,1)*lattice(2,2)
  r_lattice(2,3)=lattice(3,1)*lattice(1,2)-lattice(1,1)*lattice(3,2)
  r_lattice(3,3)=lattice(1,1)*lattice(2,2)-lattice(2,1)*lattice(1,2)
  vol=dot_product(lattice(1:3,1),r_lattice(1:3,1))
  r_lattice(1:3,1:3)=r_lattice(1:3,1:3)/vol
  do i=1,3; write(0,*)lattice(1:3,i),r_lattice(1:3,i); end do
  write(0,*)vol
  ! PRIMECOORD
  read(pid_geom,*)title
  read(pid_geom,*) natom_total
  if(natom_total .ne. sum(natom(1:ntype))) call fatal
  write(0,*)natom(1:ntype)
  ! format of coord
  cartesian=.true. 

  ! read atoms_fleur.d
! do i=1,ntype
!  read(pid_atoms,*)idum,aname(i),radius(i),color(:,i),radi_R(i)
! end do
  read(pid_atoms,'(A)')title ; backspace(pid_atoms)
  call take_before_comment(title,title_) ; title = title_
  if(nitems(title) .eq. 4) then
    read(pid_atoms,*)radius(0),bcolor(1:3)              !bond width
  elseif(nitems(title) .eq. 5) then
    read(pid_atoms,*)radius(0),bcolor(1:3), bcolor_type !bond width
  endif
  
! read(pid_atoms,*)radius(0),bcolor(1:3) !bond width
  read(pid_atoms,*)ibound                !bounding box factor: 0(off), 1(on)
  if((ibound .gt. 0) .AND. (ibound .lt. 4))then
  read(pid_atoms,*)c_min,c_max
  elseif(ibound .eq. 4)then
  read(pid_atoms,*)x_min,x_max
  read(pid_atoms,*)y_min,y_max
  read(pid_atoms,*)z_min,z_max
  endif
  close(pid_atoms)

  ! coord info
  do itype=1,ntype
    do iatom=1,natom(itype)
      read(pid_geom,*)idum, c(1:3)
      if(cartesian)then
        do i=1,3
          tvec(i)=dot_product(r_lattice(1:3,i),c(1:3))
        end do
      else
        tvec(1:3)=c(1:3)
      end if
      tvec(:)=tvec(:)+(/shift_nx,shift_ny,shift_nz/)
      coord(1:3,iatom,itype)=tvec(1:3)-nint(tvec(1:3))
    end do
  end do
  close(pid_geom)


  ! generate bs file
  open(9,file='struct.bs',status='replace',form='formatted')
  do i=1,nx
  do j=1,ny
  do k=1,nz
    do itype=1,ntype
      do iatom=1,natom(itype)
        c(1)=dot_product((/i-1,j-1,k-1/)+coord(1:3,iatom,itype),lattice(1,1:3))
        c(2)=dot_product((/i-1,j-1,k-1/)+coord(1:3,iatom,itype),lattice(2,1:3))
        c(3)=dot_product((/i-1,j-1,k-1/)+coord(1:3,iatom,itype),lattice(3,1:3))
        if((c(3) .gt. z_min .AND. c(3) .lt. z_max) .AND. &
           (c(2) .gt. y_min .AND. c(2) .lt. y_max) .AND. &
           (c(1) .gt. x_min .AND. c(1) .lt. x_max) .AND. (ibound .eq. 4))then
          write(*,'("atom ",a6,3f20.11," # ",i5)')aname(itype),c(1:3),iatom
          write(9,'("atom ",a6,3f20.11," # ",i5)')aname(itype),c(1:3),iatom
        ELSEIF(ibound .eq. 0)then
          write(*,'("atom ",a6,3f20.11," # ",i5)')aname(itype),c(1:3),iatom
          write(9,'("atom ",a6,3f20.11," # ",i5)')aname(itype),c(1:3),iatom
        ELSEIF(ibound .gt. 0 .AND. ibound .lt. 4)then
          if(c(ibound) .gt. c_min .AND. c(ibound) .lt. c_max)then
            write(*,'("atom ",a6,3f20.11," # ",i5)')aname(itype),c(1:3),iatom
            write(9,'("atom ",a6,3f20.11," # ",i5)')aname(itype),c(1:3),iatom
          endif
        endif
      end do
    end do
  end do
  end do
  end do



  do itype=1,ntype
    write(*,'("spec ",a,1x,f9.5,1x,3(f8.4,1x))')&
     aname(itype)(1:len_trim(aname(itype))),radi_R(itype)*0.30d0,&
     color(:,itype)
    write(9,'("spec ",a,1x,f9.5,1x,3(f8.4,1x))')&
     aname(itype)(1:len_trim(aname(itype))),radi_R(itype)*0.30d0,&
     color(:,itype)
  end do
  do itype=1,ntype
    do i=itype,ntype

      if(bcolor_type .gt. 0.0d0) then
        bcolor(1:3) = (color(:,itype) + color(:,i)) * 0.5d0
        radi_R_min = min(radi_R(itype)*0.30d0, radi_R(i)*0.30d0) ! minimum radius of atom ball
        if( radius(0) .gt. radi_R_min) then  ! if radi_R_min is narrower then cylinder radii
          radi_RR = radi_R_min               !    set cylinder radii as radi_R_min
        else                                 ! if radi_R_min is wider than cylinder radii
          radi_RR = radius(0)                !    use cylinder radii
        endif
      else
        radi_R_min = min(radi_R(itype)*0.30d0, radi_R(i)*0.30d0) ! minimum radius of atom ball
        if( radius(0) .gt. radi_R_min) then  ! if radi_R_min is narrower then cylinder radii
          radi_RR = radi_R_min               !    set cylinder radii as radi_R_min
        else                                 ! if radi_R_min is wider than cylinder radii
          radi_RR = radius(0)                !    use cylinder radii
        endif
      endif

      write(*,'("bonds ",2a6," 0.0 ",2(f9.5,1x),3(f8.4,1x))')&
       aname(itype),aname(i),(radius(itype)+radius(i)),& 
       radi_RR,bcolor(1:3)
      write(9,'("bonds ",2a6," 0.0 ",2(f9.5,1x),3(f8.4,1x))')&
       aname(itype),aname(i),(radius(itype)+radius(i)),&
       radi_RR,bcolor(1:3)
    end do
  end do
contains
  subroutine fatal; implicit none
    write(0,*)'please check your poscar file'; stop
  end subroutine fatal
end program fleur2bs
subroutine skipwhite(line,lline,column)
  implicit none
  character(len=*),intent(in)::line
  integer,intent(in)::lline
  integer,intent(inout)::column
  integer::i
  do i=column,lline
    if(ichar(line(i:i))>32)then
      column=i
      return
    end if
  end do
  column=0
  return
end subroutine skipwhite
subroutine gettoken(line,lline,column,tend)
  implicit none
  character(len=*),intent(in)::line
  integer,intent(in)::lline
  integer,intent(in)::column
  integer,intent(out)::tend
  integer::i
  do i=column,lline
    if(ichar(line(i:i))<=32)then
      tend = i - 1
      return
    end if
  end do
  tend=lline
  return
end subroutine gettoken

FUNCTION is_string(string)
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: string
LOGICAL :: is_string
REAL :: x
INTEGER :: e
is_string = .false.
READ(string,*,IOSTAT=e ) x
is_string = e .ne. 0
END FUNCTION is_string

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

subroutine take_before_comment(string, strip)
! comment off from the "string" and return as "strip".
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

! vim:shiftwidth=2:smarttab:autoindent:expandtab
