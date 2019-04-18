#include "alias.inc"
subroutine read_energy(PINPT, PGEOM,PKPTS,EDFT, EDFT_all, PWGHT)
  use parameters
  use mpi_setup
  implicit none
  character*132     fnameu,fnamed
  character*132    inputline
  character*40  desc_str,dummy, dummy_(1000)
  character(*), parameter :: func = 'read_energy'
  integer*4, parameter :: max_eig = 1000000
  integer*4        ie, ik, k, i_continue, line_tot
  integer*4                 size_range
  integer*4, allocatable :: irange(:),irange_up(:), irange_dn(:)

  type(incar)   :: PINPT
  type(poscar)  :: PGEOM
  type(kpoints) :: PKPTS
  type(energy)  :: EDFT_
  type(energy)  :: EDFT_all
  type(energy)  :: EDFT
  type(weight)  :: PWGHT

  if(PINPT%flag_collinear) then
   fnameu = trim(PINPT%efilenmu)
   fnamed = trim(PINPT%efilenmd)
   allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(PGEOM%neig*2 .gt. max_eig) then
     if_main write(6,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'
     if_main write(6,'(A,A)')'  !WARNING!  Exit program... ',func
     stop
   endif
  elseif(PINPT%flag_noncollinear) then
   fnameu = trim(PINPT%efilenmu)
   allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(PGEOM%neig .gt. max_eig) then
     if_main write(6,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'
     if_main write(6,'(A,A)')'  !WARNING!  Exit program... ',func
     stop
   endif
  else
   fnameu = trim(PINPT%efilenmu)
   allocate( EDFT%E(PGEOM%neig,  PKPTS%nkpoint) )
   if(PGEOM%neig .gt. max_eig) then
     if_main write(6,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'
     if_main write(6,'(A,A)')'  !WARNING!  Exit program... ',func
     stop
   endif
  endif

  allocate( EDFT_%E(max_eig, PKPTS%nkpoint) )
  if(PINPT%ispin .eq. 2)then
    open(pid_energy,  file=trim(fnameu), status='old', iostat=i_continue)
    open(pid_energy+1,file=trim(fnamed), status='old', iostat=i_continue)
  else
    open(pid_energy,file=trim(fnameu), status='old', iostat=i_continue)
  endif

  if_main write(6,*)' '
  if(PINPT%flag_collinear) then
    if_main write(6,*)'*- READING TARGET ENERGY FILE: ',trim(fnameu), ' and ', trim(fnamed)
  else
    if_main write(6,*)'*- READING TARGET ENERGY FILE: ',trim(fnameu)
  endif

  ie=0
  line_tot = 0

 loop1:do
         ik=0
         ie=ie + 1
         do while (ik .lt. PKPTS%nkpoint)
           read(pid_energy,'(A)',iostat=i_continue) inputline
           if(i_continue<0) exit loop1         ! end of file reached
           if(i_continue>0) then 
             if_main write(6,*)'Unknown error reading file:',trim(fnameu),func
             stop
           endif
           ! check comment and emtpy line
           read(inputline,*,iostat=i_continue) desc_str
           if (desc_str(1:1).eq.'#') cycle  ! skip comment
           if (i_continue .ne. 0)    cycle  ! skip empty

           ik = ik + 1
           line_tot = line_tot + 1
           if(PINPT%read_energy_column_index .eq. 2) then
             read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik)
           elseif(PINPT%read_energy_column_index .eq. 3) then
             read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik)
           elseif(PINPT%read_energy_column_index .eq. 4) then
             read(inputline,*,iostat=i_continue) dummy, dummy, dummy, EDFT_%E(ie, ik)
           elseif(PINPT%read_energy_column_index .eq. 5) then
             read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
           elseif(PINPT%read_energy_column_index .eq. 6) then
             read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
           elseif(PINPT%read_energy_column_index .eq. 7) then
             read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
           elseif(PINPT%read_energy_column_index .eq. 8) then
             read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
           elseif(PINPT%read_energy_column_index .ge. 9) then
             read(inputline,*,iostat=i_continue) dummy_(1:PINPT%read_energy_column_index-1), EDFT_%E(ie, ik)
           endif
         enddo

         if(ik .ne. PKPTS%nkpoint) then
           if_main write(6,'(A,A )')'   !WARN! error in reading ',trim(fnameu)
           if_main write(6,'(A,I8)')'   Exit.. NKPTS !=',ik
           stop
         endif
       enddo loop1

       if(PINPT%flag_collinear) then
         ie = ie - 1
         line_tot = 0
 loop2:  do
           ik=0
           ie=ie + 1
           do while (ik .lt. PKPTS%nkpoint)
             read(pid_energy+1,'(A)',iostat=i_continue) inputline
             if(i_continue<0) exit loop2         ! end of file reached
             if(i_continue>0) then 
               if_main write(6,*)'Unknown error reading file:',trim(fnamed),func
               stop
             endif
             ! check comment and emtpy line
             read(inputline,*,iostat=i_continue) desc_str
             if (desc_str(1:1).eq.'#') cycle  ! skip comment
             if (i_continue .ne. 0)    cycle  ! skip empty
        
             ik = ik + 1
             line_tot = line_tot + 1
             if(PINPT%read_energy_column_index_dn .eq. 2) then
               read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik)
             elseif(PINPT%read_energy_column_index_dn .eq. 3) then
               read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik)
             endif
           enddo
           if(ik .ne. PKPTS%nkpoint) then
             if_main write(6,'(A,A )')'   !WARN! error in reading ',trim(fnamed)
             if_main write(6,'(A,I8)')'   Exit.. NKPTS !=',ik
             stop
           endif
         enddo loop2
       endif

       PGEOM%neig_target=ie - 1
       if_main write(6,'(A,I8)')' N_TARGET:',PGEOM%neig_target

  allocate( EDFT_all%E(PGEOM%neig_target, PKPTS%nkpoint) )

  if(PINPT%flag_collinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig
  elseif(PINPT%flag_noncollinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig*2
  else
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig
  endif

  size_range = size( (/(k, k=PWGHT%iband,PWGHT%fband)/) )
  if(PINPT%flag_collinear) then
    allocate( irange_up(size_range) )
    allocate( irange_dn(size_range) )
    allocate( irange(size_range*PINPT%ispin) )
  elseif(PINPT%flag_noncollinear) then
    allocate( irange(size_range*PINPT%ispin) )
  else
    allocate( irange(size_range) )
  endif
  
  if(PINPT%flag_collinear) then
    irange_up = (/(k, k=PWGHT%iband,PWGHT%fband)/)
    irange_dn = irange_up + PGEOM%neig_target / 2
    irange = (/irange_up, irange_dn/)
  elseif(PINPT%flag_noncollinear) then
    irange = (/(k, k=PWGHT%iband,PWGHT%fband)/)
  else
    irange = (/(k, k=PWGHT%iband,PWGHT%fband)/)
  endif

  if(PINPT%flag_collinear) then
    if(PINPT%flag_scissor) then
      EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
            EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
      EDFT_%E(PGEOM%neig+PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
            EDFT_%E(PGEOM%neig+PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
    endif
    EDFT%E(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
  elseif(PINPT%flag_noncollinear) then
    if(PINPT%flag_scissor) then
      EDFT_%E(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
            EDFT_%E(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
    endif
    EDFT%E(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
  else
    if(PINPT%flag_scissor) then
      EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
            EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
    endif
    EDFT%E(1:PGEOM%neig,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
  endif
  EDFT_all%E(1:PGEOM%neig_target,1:PKPTS%nkpoint) = EDFT_%E(1:PGEOM%neig_target,1:PKPTS%nkpoint)

  if_main write(6,*)'*- END READING TARGET ENERGY FILE --------------'
  if_main write(6,*)' '
  close(pid_energy)
  if(PINPT%flag_collinear) close(pid_energy+1)
  deallocate(EDFT_%E)
return
endsubroutine

! read band_structure_TBA.dat generated by TBFIT.
subroutine read_energy_tbfit(E, V, nspin, norb, nband, nk, kmode, flag_vector)
   use parameters, only : pid_energy
   implicit none
   integer*4              nspin, norb, nband, nk
   integer*4              pid, i_continue
   integer*4              ie,ik,is
   integer*4              ii, fi, im, fm
   integer*4              linecount
   real*8                 E(nband*nspin,nk)
   real*8                 V(norb,nband*nspin,nk)
   real*8                 kpoint
   character*132          inputline
   character*40           desc_str
   character*4            kmode
   external               nitems
   integer*4              nitems, idummy
   logical                flag_go, flag_vector
   character*80           fnameu, fnamed

   pid = pid_energy
   linecount = 0
   ie= 0
   flag_go = .false.

   if(nspin .eq. 2) then ! collinear
     fnameu = 'band_structure_TBA.up.dat'
     fnamed = 'band_structure_TBA.dn.dat'
   elseif(nspin .eq. 1) then ! nm or non-collinear
     fnameu = 'band_structure_TBA.dat'
   endif
  

   do is = 1, nspin
     if(is .eq. 1) then
       open(pid,file = trim(fnameu), iostat=i_continue)
       write(6,'(2A)')'     ... loading: ', trim(fnameu)
       if(i_continue .ne. 0) then 
         write(6,'(A)')'  !!! ERROR IN READING (up) ', trim(fnameu)
         write(6,'(A)')'      EXIT PROGRAM : read_energy_tbfit '
         stop
       endif
     elseif(is .eq. 2) then
       open(pid,file = trim(fnamed), iostat=i_continue)
       write(6,'(2A)')'  ... loading: ', trim(fnameu)
       if(i_continue .ne. 0) then 
         write(6,'(A)')'  !!! ERROR IN READING (dn) ', trim(fnamed)
         write(6,'(A)')'      EXIT PROGRAM : read_energy_tbfit '
         stop
       endif
     endif

     ii = 1+nband*(is-1) ; fi = nband+nband*(is-1)
     im = 1              ; fm = norb            
     call load_band_structure(E(ii:fi,:), V(im:fm,ii:fi,:), nband, norb, nk, kmode, pid, flag_vector)

     close(pid)
   enddo 
   return
endsubroutine

subroutine load_band_structure(E, V, nband, norb, nk, kmode, pid, flag_vector)
   implicit none
   integer*4              pid
   integer*4              ie, ik
   integer*4              nband, norb, nk
   integer*4              nskip
   integer*4              idummy
   integer*4              nitems
   external               nitems
   character*132          inputline
   real*8                 E(nband, nk)
   real*8                 V(norb, nband, nk)
   logical                flag_go, flag_vector
   character*4            kmode
   real*8, allocatable :: kpoint(:)

   if(trim(kmode) .eq. 'grid') then   
     allocate(kpoint(3))
     nskip = 3
   elseif(trim(kmode) .eq. 'line') then
     allocate(kpoint(1))
     nskip = 1
   endif

   do ie = 1, nband
     flag_go = .false.
     do while (.not.flag_go)
       read(pid,'(A)')inputline
       idummy = nitems(inputline)
       if(idummy .le. 0) then
         flag_go = .false.
       elseif(idummy .gt. 0) then
         flag_go = .true.
         backspace(pid)
       endif
     enddo

     do ik = 1, nk
       read(pid,'(A)')inputline
       idummy = nitems(inputline) - nskip
       if(idummy .gt. 1) then
         backspace(pid)
         if(flag_vector) then
           read(pid,*)kpoint(:), E(ie,ik),V(1:norb,ie,ik)
         elseif(.not. flag_vector) then
           read(pid,*)kpoint(:), E(ie,ik)
         endif
       elseif(idummy .eq. 0) then
         backspace(pid)
         read(pid,*)kpoint(:)
         E(ie,ik) = -999d0
         if(flag_vector) V(1:norb,ie,ik) = 0d0
       endif
     enddo
   enddo

   deallocate(kpoint)

   return
endsubroutine
