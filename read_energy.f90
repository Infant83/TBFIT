subroutine read_energy(PINPT, PGEOM,PKPTS,EDFT, EDFT_all, PWGHT)
  use parameters
  implicit none
  character*132     fnameu,fnamed
  character*132    inputline
  character*40  desc_str,dummy
  character(*), parameter :: func = 'read_energy'
  integer*4, parameter :: max_eig = 1000000
  integer*4        ie, ik, i_continue, line_tot
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
     write(6,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'
     write(6,'(A,A)')'  !WARNING!  Exit program... ',func
     stop
   endif
  elseif(PINPT%flag_noncollinear) then
   fnameu = trim(PINPT%efilenmu)
   allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(PGEOM%neig .gt. max_eig) then
     write(6,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'
     write(6,'(A,A)')'  !WARNING!  Exit program... ',func
     stop
   endif
  else
   fnameu = trim(PINPT%efilenmu)
   allocate( EDFT%E(PGEOM%neig,  PKPTS%nkpoint) )
   if(PGEOM%neig .gt. max_eig) then
     write(6,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'
     write(6,'(A,A)')'  !WARNING!  Exit program... ',func
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

  write(6,*)' '
  if(PINPT%flag_collinear) then
    write(6,*)'*- READING TARGET ENERGY FILE: ',trim(fnameu), ' and ', trim(fnamed)
  else
    write(6,*)'*- READING TARGET ENERGY FILE: ',trim(fnameu)
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
             write(6,*)'Unknown error reading file:',trim(fnameu),func
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
           endif
         enddo
         if(ik .ne. PKPTS%nkpoint) then
           write(6,'(A,A )')'   !WARN! error in reading ',trim(fnameu)
           write(6,'(A,I8)')'   Exit.. NKPTS !=',ik
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
               write(6,*)'Unknown error reading file:',trim(fnamed),func
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
             write(6,'(A,A )')'   !WARN! error in reading ',trim(fnamed)
             write(6,'(A,I8)')'   Exit.. NKPTS !=',ik
             stop
           endif
         enddo loop2
       endif

       PGEOM%neig_target=ie - 1
       write(6,'(A,I8)')' N_TARGET:',PGEOM%neig_target

  allocate( EDFT_all%E(PGEOM%neig_target, PKPTS%nkpoint) )

  if(PINPT%flag_collinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig
  elseif(PINPT%flag_noncollinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig*2
  else
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig
  endif

  size_range = size( (/PWGHT%iband:PWGHT%fband/) )
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
    irange_up = (/PWGHT%iband:PWGHT%fband/)
    irange_dn = irange_up + PGEOM%neig_target / 2
    irange = (/irange_up, irange_dn/)
  elseif(PINPT%flag_noncollinear) then
    irange = (/PWGHT%iband:PWGHT%fband/)
  else
    irange = (/PWGHT%iband:PWGHT%fband/)
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

  write(6,*)'*- END READING TARGET ENERGY FILE --------------'
  write(6,*)' '
  close(pid_energy)
  if(PINPT%flag_collinear) close(pid_energy+1)
  deallocate(EDFT_%E)
return
endsubroutine
