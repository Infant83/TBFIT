!Modified by Hyun-Jung Kim. Dept. of Physics Hanyang Univ.2008,1,21
!minor Modified by Hyun-Jung Kim. Dept. of Physics Hanyang Univ.2013,5,1 : radi_R added, bounding box added
!include shift option, by Hyun-Jung Kim, KIAS 2018.04.09 : shift_nx, shift_ny, shift_nz
!last update : 2018.05.16 by KHJ

program poscar2bs
  implicit none
  integer,parameter::matom=129024
  integer,parameter::mtype=16
  character(len=80)::title,fname,atom_name
  integer::nx,ny,nz, ishift
  real*8   shift_nx,shift_ny,shift_nz
  character(len=132)::line, desc_str
  integer::lline,column,tend,i,j,k,ibound
  real(8),dimension(3)::sc
  real(8),dimension(3,3)::lattice,r_lattice
  real(8)::vol
  integer::ntype
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

  ishift = 0
  pid_atoms = 77
  pid_geom  = 78

  ! read atoms.d header info.
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

  ! read POSCAR file
  open(pid_geom,file=fname,status='old',form='formatted')
  rewind(pid_geom)
  ! title
  read(pid_geom,'(a)')title
  write(0,'(a)')title(1:len_trim(title))
  ! scale factor(s)
  read(pid_geom,'(a)')line
  lline=len_trim(line);                        if(lline==0)call fatal
  write(0,*)line(1:len_trim(line)),lline
  !
  column=1; call skipwhite(line,lline,column); if(column==0)call fatal
  call gettoken(line,lline,column,tend); read(line(column:tend),*)sc(1)
  write(0,*)line(column:tend),column,tend,sc(1)
  !
  column=tend+1; call skipwhite(line,lline,column)
  if(column==0)then; sc(2)=sc(1); sc(3)=sc(1)
  else
    call gettoken(line,lline,column,tend); read(line(column:tend),*)sc(2)
    write(0,*)line(column:tend),column,tend,sc(2)
    !
    column=tend+1; call skipwhite(line,lline,column); if(column==0)call fatal
    call gettoken(line,lline,column,tend)
    write(0,*)line(column:tend),column,tend
    read(line(column:lline),*)sc(3)
  end if
  write(0,*)sc(1:3)
  ! lattice vectors
  do i=1,3; read(pid_geom,*)lattice(1:3,i); end do
  do i=1,3; lattice(i,1:3)=lattice(i,1:3)*sc(i); end do
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
  ! atom name : Infant added
  read(pid_geom,*)atom_name
  ! natom(s)
  read(pid_geom,'(a)')line
  lline=len_trim(line); if(lline==0)call fatal
  ntype=0
  column=1
  do
    call skipwhite(line,lline,column); if(column==0)exit
    call gettoken(line,lline,column,tend)
    ntype=ntype+1; read(line(column:tend),*) natom(ntype)
    column=tend+1
  end do
  if(ntype==0)call fatal
  write(0,*)natom(1:ntype)
  ! format of coord
  read(pid_geom,'(a)')tchar_
  tchar_ = adjustl(trim(tchar_))
  tchar = tchar_(1:1)
  selective=.false.
  if(adjustl(trim(tchar))=='s'.or.adjustl(trim(tchar))=='S')then
    selective=.true.
    read(pid_geom,'(a)')tchar_
    tchar_ = adjustl(trim(tchar_))
    tchar = tchar_(1:1)
  end if
  ! type of coord
  cartesian=.false.
  if(adjustl(trim(tchar))=='c'.or.adjustl(trim(tchar))=='C'.or. &
     adjustl(trim(tchar))=='k'.or.adjustl(trim(tchar))=='K')cartesian=.true.

  ! read atoms.d
  do i=1,ntype
   read(pid_atoms,*)idum,aname(i),radius(i),color(:,i),radi_R(i)
  end do
  read(pid_atoms,*)radius(0),bcolor(1:3) !bond width
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
      if(selective)then
        read(pid_geom,*)c(1:3),m(1:3)
      else
        read(pid_geom,*)c(1:3)
      end if
      if(cartesian)then
        do i=1,3
          tvec(i)=dot_product(r_lattice(1:3,i),c(1:3)*sc(1:3))
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
  open(9,file='poscar.bs',status='replace',form='formatted')
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
!    aname(itype)(1:len_trim(aname(itype))),radius(itype)*0.30d0,&
     aname(itype)(1:len_trim(aname(itype))),radi_R(itype)*0.30d0,&
     color(:,itype)
    write(9,'("spec ",a,1x,f9.5,1x,3(f8.4,1x))')&
!    aname(itype)(1:len_trim(aname(itype))),radius(itype)*0.30d0,&
     aname(itype)(1:len_trim(aname(itype))),radi_R(itype)*0.30d0,&
     color(:,itype)
  end do
  do itype=1,ntype
    do i=itype,ntype
      !write(*,'("bonds ",2a6," 0.0 3.0 0.1 1.0 1.0 1.0")')aname(itype),aname(i)
      write(*,'("bonds ",2a6," 0.0 ",2(f9.5,1x),3(f8.4,1x))')&
       aname(itype),aname(i),(radius(itype)+radius(i)),& 
       radius(0),bcolor(1:3)
      !0.5d0*(color(:,itype)+color(:,i))
      write(9,'("bonds ",2a6," 0.0 ",2(f9.5,1x),3(f8.4,1x))')&
       aname(itype),aname(i),(radius(itype)+radius(i)),&
       radius(0),bcolor(1:3)
    end do
  end do
contains
  subroutine fatal; implicit none
    write(0,*)'please check your poscar file'; stop
  end subroutine fatal
end program poscar2bs
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
! vim:shiftwidth=2:smarttab:autoindent:expandtab
