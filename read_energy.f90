#include "alias.inc"
subroutine read_energy_vasp(PINPT, PGEOM,PKPTS,EDFT, EDFT_all, PWGHT) ! read EIGENVAL file
  use parameters
  use mpi_setup
  use print_io
  implicit none
  character*132     fname, fnameu,fnamed
  character*132    inputline
  character*40  desc_str,dummy, dummy_(1000)
  character(*), parameter :: func = 'read_energy'
  integer*4, parameter :: max_eig = 1000000
  integer*4        ie, ik, k, i_continue, line_tot
  integer*4                 size_range
  integer*4        idummy, ispin, nkp_target, neig_target, nelect, ie_start
  integer*4        iread_occ
  integer*4, allocatable :: irange(:),irange_up(:), irange_dn(:), occ(:,:), occ_(:,:)

  type(incar)   :: PINPT
  type(poscar)  :: PGEOM
  type(kpoints) :: PKPTS
  type(energy)  :: EDFT_, EDFT__
  type(energy)  :: EDFT_all
  type(energy)  :: EDFT
  type(weight)  :: PWGHT

  iread_occ=1 ! default ; if LORBIT = 0 and the EIGENVAL does not contain occupation info, 
              !           set iread_occ = 0 when energy surface extraction (-scf)

  if(PINPT%flag_collinear) then
    allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
    if(PINPT%flag_fit_degeneracy) then
      allocate( EDFT%D(3,PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
    endif
    if(PGEOM%neig*2 .gt. max_eig) then
      write(message,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'  ; write_msg
      write(message,'(A,A)')'  !WARNING!  Exit program... ',func  ; write_msg
      stop
    endif
  elseif(PINPT%flag_noncollinear) then
    allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
    if(PINPT%flag_fit_degeneracy) then
      allocate( EDFT%D(3,PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
    endif
    if(PGEOM%neig .gt. max_eig) then
      write(message,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'  ; write_msg
      write(message,'(A,A)')'  !WARNING!  Exit program... ',func  ; write_msg
      stop
    endif
  else
    allocate( EDFT%E(PGEOM%neig,  PKPTS%nkpoint) )
    if(PINPT%flag_fit_degeneracy) then
      allocate( EDFT%D(3,PGEOM%neig,  PKPTS%nkpoint) )
    endif
    if(PGEOM%neig .gt. max_eig) then
      write(message,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'  ; write_msg
      write(message,'(A,A)')'  !WARNING!  Exit program... ',func  ; write_msg
      stop
    endif
  endif
  fname = trim(PINPT%efilenmu)

  allocate( EDFT_%E(max_eig, PKPTS%nkpoint) )
  allocate( EDFT__%E(max_eig, PKPTS%nkpoint) )
  allocate( occ(max_eig, PKPTS%nkpoint) )
  allocate( occ_(max_eig, PKPTS%nkpoint) )
  open(pid_energy,file=trim(fname), status='old', iostat=i_continue)

  write(message,*)' '  ; write_msg
  write(message,*)'*- READING TARGET ENERGY FILE: ',trim(fname)  ; write_msg

  ie_start=PINPT%itarget_e_start
  line_tot = 0

  ! reading header for EIGENVAL
  read(pid_energy,*) idummy, idummy, idummy, ispin
  do k=1, 4
    read(pid_energy, *) dummy
  enddo
  read(pid_energy,*) nelect, nkp_target, neig_target
  ! read energy for each k-point
  do ik=1,nkp_target
    read(pid_energy,*) dummy
    do ie = 1, neig_target
      if( iread_occ .eq. 1) then
        if(ispin .eq. 2) then
          read(pid_energy,*) dummy, EDFT__%E(ie, ik), EDFT__%E(ie+neig_target, ik), occ_(ie,ik), occ_(ie+neig_target,ik)
        else
          read(pid_energy,*) dummy, EDFT__%E(ie, ik), occ_(ie,ik)
        endif
      else
        if(ispin .eq. 2) then
          read(pid_energy,*) dummy, EDFT__%E(ie, ik), EDFT__%E(ie+neig_target, ik)
        else
          read(pid_energy,*) dummy, EDFT__%E(ie, ik)
        endif
      endif
    enddo
  enddo

  EDFT__%E = EDFT__%E - PINPT%efile_ef

  if(PINPT%flag_collinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig
  elseif(PINPT%flag_noncollinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig*2
  else
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig
  endif

  size_range = size( (/(k, k=PWGHT%iband,PWGHT%fband)/) )

  ! reset according to the initial target
  if(ispin .eq. 2) then
    if( (neig_target - PINPT%itarget_e_start + 1) .ge. size_range ) then
      EDFT_%E(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PINPT%itarget_e_start:neig_target,1:nkp_target)
      EDFT_%E(1+neig_target-PINPT%itarget_e_start+1: (neig_target-PINPT%itarget_e_start+1)*2,1:nkp_target) = EDFT__%E(1+neig_target:neig_target*2,1:nkp_target)
      occ(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target) = occ_(PINPT%itarget_e_start:neig_target,1:nkp_target)
      occ(1+neig_target-PINPT%itarget_e_start+1: (neig_target-PINPT%itarget_e_start+1)*2,1:nkp_target) = occ_(1+neig_target:neig_target*2,1:nkp_target)
    else ! if neig_target - PINPT%itarget_e_start + 1 is less than fband - iband + 1
      idummy = size_range - (neig_target - PINPT%itarget_e_start + 1)
      EDFT_%E(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PINPT%itarget_e_start:neig_target,1:nkp_target)
      EDFT_%E(neig_target-PINPT%itarget_e_start+1+1:size_range,1:nkp_target) = 100d0 ! dummy value is added
      EDFT_%E(size_range + 1: size_range + neig_target-PINPT%itarget_e_start+1,1:nkp_target) = EDFT__%E(1+neig_target:neig_target*2,1:nkp_target)
      EDFT_%E(size_range + neig_target-PINPT%itarget_e_start+1+1:size_range*2,1:nkp_target) = 100d0

      occ(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target)                            = occ_(PINPT%itarget_e_start:neig_target,1:nkp_target)
      occ(neig_target-PINPT%itarget_e_start+1+1:size_range,1:nkp_target)                 = 100d0 ! dummy value is added
      occ(size_range + 1: size_range + neig_target-PINPT%itarget_e_start+1,1:nkp_target) = occ_(1+neig_target:neig_target*2,1:nkp_target)
      occ(size_range + neig_target-PINPT%itarget_e_start+1+1:size_range*2,1:nkp_target)  = 100d0
    endif
  else
    if( (neig_target - PINPT%itarget_e_start + 1) .ge. size_range ) then
      EDFT_%E(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PINPT%itarget_e_start:neig_target,1:nkp_target)
      occ(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target) = occ_(PINPT%itarget_e_start:neig_target,1:nkp_target)
    else
      EDFT_%E(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PINPT%itarget_e_start:neig_target,1:nkp_target)
      EDFT_%E(neig_target-PINPT%itarget_e_start+1+1:size_range,1:nkp_target) = 100d0
 
      occ(1:neig_target-PINPT%itarget_e_start+1,1:nkp_target)            = occ_(PINPT%itarget_e_start:neig_target,1:nkp_target)
      occ(neig_target-PINPT%itarget_e_start+1+1:size_range,1:nkp_target) = 100d0
    endif
  endif

  if(neig_target - PINPT%itarget_e_start + 1 .ge. size_range) then
    PGEOM%neig_target=neig_target * ispin
    write(message,'(A,I8)')' N_TARGET:',PGEOM%neig_target  ; write_msg
  else
    PGEOM%neig_target=size_range * ispin
    write(message,'(A,I8)')' N_TARGET: (adjusted)',PGEOM%neig_target  ; write_msg
  endif

  allocate( EDFT_all%E(PGEOM%neig_target, nkp_target) )

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

  write(message,*)'*- END READING TARGET ENERGY FILE --------------'  ; write_msg
  write(message,*)' '  ; write_msg
  close(pid_energy)
! if(PINPT%flag_collinear) close(pid_energy+1)
  deallocate(EDFT_%E)
  if(allocated(EDFT__%E)) deallocate(EDFT__%E)
  if(allocated(occ)) deallocate(occ)
  if(allocated(occ_)) deallocate(occ_)

return
endsubroutine

subroutine read_energy(PINPT, PGEOM,PKPTS,EDFT, EDFT_all, PWGHT)
  use parameters
  use mpi_setup
  use print_io
  implicit none
  character*132             fnameu,fnamed
  character*132             inputline
  character*40              desc_str,dummy, dummy_(1000)
  character(*), parameter:: func = 'read_energy'
  integer*4, parameter   :: max_eig = 1000000
  integer*4                 ie, ik, k, i_continue
  integer*4                 size_range
  integer*4                 mpierr
  integer*4                 ik_, kreduce
  integer*4, allocatable :: irange(:),irange_up(:), irange_dn(:)
  type(incar)   :: PINPT
  type(poscar)  :: PGEOM
  type(kpoints) :: PKPTS
  type(energy)  :: EDFT_
  type(energy)  :: EDFT_all
  type(energy)  :: EDFT
  type(weight)  :: PWGHT
  logical          flag_order

  flag_order = PINPT%flag_get_band_order 

  ! NOTE: kreduce with greater than 1 is only applicaple with PINPT%
  if(PKPTS%kreduce .gt. 1) then
    kreduce = PKPTS%kreduce
    if( trim(PINPT%kline_type) .eq. 'vasp' .or. trim(PINPT%kline_type) .eq. 'VASP') then
      write(message,'(A)')'  !WARNING!  KREDUCE should be 1 if you read VASP type of target energy file.'  ; write_msg
      write(message,'(A)')'             Please read user manual with KREDUCE tag and KFILE tag.'  ; write_msg
      write(message,'(A)')'             Exit program...'  ; write_msg
      kill_job
    endif
  else
    kreduce = 1
  endif

  if(PINPT%flag_collinear) then
   fnameu = trim(PINPT%efilenmu)
   fnamed = trim(PINPT%efilenmd)
   allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(flag_order) allocate( EDFT%E_ORD(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) ) 
   if(PINPT%flag_fit_degeneracy) then
     allocate( EDFT%D(3,PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   endif
   if(PGEOM%neig*2 .gt. max_eig) then
     write(message,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'  ; write_msg
     write(message,'(A,A)')'  !WARNING!  Exit program... ',func  ; write_msg
     stop
   endif
  elseif(PINPT%flag_noncollinear) then
   fnameu = trim(PINPT%efilenmu)
   allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(flag_order) allocate( EDFT%E_ORD(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(PINPT%flag_fit_degeneracy) then
     allocate( EDFT%D(3,PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   endif
   if(PGEOM%neig .gt. max_eig) then
     write(message,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'  ; write_msg
     write(message,'(A,A)')'  !WARNING!  Exit program... ',func  ; write_msg
     stop
   endif
  else
   fnameu = trim(PINPT%efilenmu)
   allocate( EDFT%E(PGEOM%neig,  PKPTS%nkpoint) )
   if(flag_order) allocate( EDFT%E_ORD(PGEOM%neig,  PKPTS%nkpoint) )
   if(PINPT%flag_fit_degeneracy) then
     allocate( EDFT%D(3,PGEOM%neig,  PKPTS%nkpoint) )
   endif
   if(PGEOM%neig .gt. max_eig) then
     write(message,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'  ; write_msg
     write(message,'(A,A)')'  !WARNING!  Exit program... ',func  ; write_msg
     stop
   endif
  endif

  allocate( EDFT_%E(max_eig, PKPTS%nkpoint) )
  if(flag_order) allocate( EDFT_%E_ORD(max_eig, PKPTS%nkpoint) )

  if(PINPT%ispin .eq. 2)then
    open(pid_energy,  file=trim(fnameu), status='old', iostat=i_continue)
    open(pid_energy+1,file=trim(fnamed), status='old', iostat=i_continue)
  else
    open(pid_energy,file=trim(fnameu), status='old', iostat=i_continue)
  endif

!# ORDERED ENERGY => E_ORD
!# PRISITINE ENERGY => E 

  write(message,*)' '  ; write_msg
  if(PINPT%flag_collinear) then
    write(message,*)'*- READING TARGET ENERGY FILE: ',trim(fnameu), ' and ', trim(fnamed)  ; write_msg
  else
    write(message,*)'*- READING TARGET ENERGY FILE: ',trim(fnameu)  ; write_msg
  endif

  ie=0
 loop1:do
         ik=0; ik_ = 0
         ie=ie + 1
         do while (ik .lt. PKPTS%nkpoint)
           read(pid_energy,'(A)',iostat=i_continue) inputline
           if(i_continue<0) exit loop1         ! end of file reached
           if(i_continue>0) then 
             write(message,*)'Unknown error reading file:',trim(fnameu),func  ; write_msg
             stop
           endif
           ! check comment and emtpy line
           read(inputline,*,iostat=i_continue) desc_str
           if (desc_str(1:1).eq.'#') cycle  ! skip comment
           if (i_continue .ne. 0)    cycle  ! skip empty

           ik_ = ik_ + 1
           if( mod(ik_ - 1, kreduce) .eq. 0)  then
             ik = ik + 1
             if(PINPT%read_energy_column_index .eq. 2) then
               if(.not. flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik)
               elseif(flag_order) then 
                 ! this is for read re-ordered target band and save it to E_ORD, 
                 ! hence, be sure that file should contain this information too!
                 read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik), EDFT_%E_ORD(ie, ik)
               endif
             elseif(PINPT%read_energy_column_index .eq. 3) then
               if(.not. flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik)
               elseif(flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik), EDFT_%E_ORD(ie, ik)
               endif
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
           else
             read(inputline,*,iostat=i_continue) dummy
           endif
         enddo

         if(ik .ne. PKPTS%nkpoint) then
           write(message,'(A,A )')'   !WARN! error in reading ',trim(fnameu)  ; write_msg
           write(message,'(A,I8)')'   Exit.. NKPTS !=',ik  ; write_msg
           stop
         endif

       enddo loop1
       if(PINPT%flag_collinear) then
         ie = ie - 1
 loop2:  do
           ik=0
           ie=ie + 1
           do while (ik .lt. PKPTS%nkpoint)
             read(pid_energy+1,'(A)',iostat=i_continue) inputline
             if(i_continue<0) exit loop2         ! end of file reached
             if(i_continue>0) then 
               write(message,*)'Unknown error reading file:',trim(fnamed),func  ; write_msg
               stop
             endif
             ! check comment and emtpy line
             read(inputline,*,iostat=i_continue) desc_str
             if (desc_str(1:1).eq.'#') cycle  ! skip comment
             if (i_continue .ne. 0)    cycle  ! skip empty
        
             ik = ik + 1
             if(PINPT%read_energy_column_index_dn .eq. 2) then
               if(.not.flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik)
               elseif(flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik), EDFT_%E_ORD(ie, ik)
               endif
             elseif(PINPT%read_energy_column_index_dn .eq. 3) then
               if(.not. flag_order)then
                 read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik)
               elseif(flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik), EDFT_%E_ORD(ie, ik)
               endif
             endif
           enddo
           if(ik .ne. PKPTS%nkpoint) then
             write(message,'(A,A )')'   !WARN! error in reading ',trim(fnamed)  ; write_msg
             write(message,'(A,I8)')'   Exit.. NKPTS !=',ik  ; write_msg
             stop
           endif
         enddo loop2
       endif

       PGEOM%neig_target=ie - 1
       write(message,'(A,I8)')' N_TARGET:',PGEOM%neig_target  ; write_msg

  allocate( EDFT_all%E(PGEOM%neig_target, PKPTS%nkpoint) )
  if(flag_order) allocate( EDFT_all%E_ORD(PGEOM%neig_target, PKPTS%nkpoint) )

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

  EDFT_%E = EDFT_%E - PINPT%efile_ef
  if(flag_order) EDFT_%E_ORD = EDFT_%E_ORD - PINPT%efile_ef

  if(PINPT%flag_collinear) then
    if(PINPT%flag_scissor) then
      EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
            EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
      EDFT_%E(PGEOM%neig+PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
            EDFT_%E(PGEOM%neig+PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
      if(flag_order) then
        EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
              EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
        EDFT_%E_ORD(PGEOM%neig+PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
              EDFT_%E_ORD(PGEOM%neig+PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
      endif
    endif
    EDFT%E(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
    if(flag_order) EDFT%E_ORD(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E_ORD(irange,1:PKPTS%nkpoint)
  elseif(PINPT%flag_noncollinear) then
    if(PINPT%flag_scissor) then
      EDFT_%E(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
            EDFT_%E(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
      if(flag_order) then
        EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
              EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
      endif
    endif
    EDFT%E(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
    if(flag_order) EDFT%E_ORD(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E_ORD(irange,1:PKPTS%nkpoint)
  else
    if(PINPT%flag_scissor) then
      EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
            EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
      if(flag_order) then
        EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
              EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
      endif
    endif
    EDFT%E(1:PGEOM%neig,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
    if(flag_order) EDFT%E_ORD(1:PGEOM%neig,1:PKPTS%nkpoint) = EDFT_%E_ORD(irange,1:PKPTS%nkpoint)
  endif
  EDFT_all%E(1:PGEOM%neig_target,1:PKPTS%nkpoint) = EDFT_%E(1:PGEOM%neig_target,1:PKPTS%nkpoint)
  if(flag_order) EDFT_all%E_ORD(1:PGEOM%neig_target,1:PKPTS%nkpoint) = EDFT_%E_ORD(1:PGEOM%neig_target,1:PKPTS%nkpoint)

  write(message,*)'*- END READING TARGET ENERGY FILE --------------'  ; write_msg
  write(message,*)' '  ; write_msg
  close(pid_energy)
  if(PINPT%flag_collinear) close(pid_energy+1)
  deallocate(EDFT_%E)
  if(flag_order) deallocate(EDFT_%E_ORD)
return
endsubroutine

subroutine read_energy_tbfit_V(E, V, ne_found, ispin, nspin, nbasis, nband, nk, kmode, flag_vector, flag_wf, flag_formatted, flag_exit)
! read band_structure_TBA.dat generated by TBFIT.
   use parameters, only : pid_energy, t1, t0
   use time
   use print_io
   use mpi_setup
   implicit none
   integer*4              ispin, nspin, nbasis, nband, nk
   integer*4              ispinor
   integer*4              pid, i_continue
   integer*4              ie,ik,is
   integer*4              ii, fi, im, fm
   integer*4              linecount
   integer*4              ne_found(nspin, nk)
   real*8                 E(nband*nspin,nk)
   complex*16             V(nbasis*ispin,nband*nspin,nk)
   real*8                 kpoint
   character*132          inputline
   character*40           desc_str
   character*4            kmode
   external               nitems
   integer*4              nitems, idummy
   logical                flag_go, flag_vector
   logical                flag_wf, flag_formatted
   character*80           fnameu, fnamed
   character*11           form_
   logical                flag_exit

   call time_check(t1,t0,'init')

   flag_exit = .false.

   if(flag_formatted)      form_ = 'formatted'
   if(.not.flag_formatted) form_ = 'unformatted'

   pid = pid_energy
   linecount = 0
   ie= 0
   flag_go = .false.
   if(ispin .eq. 2 .and. nspin .eq. 1) then 
     ispinor = 2
   else
     ispinor = 1
   endif   

   if(nspin .eq. 2) then ! collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA.up.dat'
       fnamed = 'band_structure_TBA.dn.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA.up.bin'
       fnamed = 'band_structure_TBA.dn.bin'
     endif
   elseif(nspin .eq. 1) then ! nm or non-collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA.bin'
     endif
   endif
  

   do is = 1, nspin
     if(is .eq. 1) then
       if(flag_formatted)       open(pid,file = trim(fnameu),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnameu), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then 
         write(message,'(A)')'  !!! ERROR IN READING (up) ', trim(fnameu)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit '  ; write_msg
         flag_exit = .true.
       endif
     elseif(is .eq. 2) then
       if(flag_formatted)       open(pid,file = trim(fnamed),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnamed), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then 
         write(message,'(A)')'  !!! ERROR IN READING (dn) ', trim(fnamed)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit '  ; write_msg
         flag_exit = .true.
       endif
     endif

     if(.not. flag_exit) then
       ii = 1+nband*(is-1) ; fi = nband+nband*(is-1)
       im = 1+nbasis*(is-1)  ; fm = nbasis*ispinor +nbasis*(is-1)

       if(flag_formatted) then
         call load_band_structure_V(E(ii:fi,:), V(im:fm,ii:fi,:), ne_found(is,:), &
                                    ispinor, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector)
!                                   ispinor, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector, flag_exit)
       elseif(.not. flag_formatted) then
         call load_band_structure_bin(E(ii:fi,:), V(im:fm,ii:fi,:), ne_found(is,:), &
                                      ispinor, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector, flag_exit)
       endif   
      
       close(pid)
     endif
   enddo 

   call time_check(t1,t0,'end')
   if(nspin .eq. 1) then
     write(message,'(A,F12.6)')"   TIME for LOADING band_structure_TBA (s): ", t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   elseif(nspin .eq. 2) then
     write(message,'(A,F12.6)')"   TIME for LOADING band_structure_TBA.up(dn) (s): ", t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   endif

   return
endsubroutine
subroutine read_energy_tbfit_V2(E, V2, ne_found, ispin, nspin, nbasis, nband, nk, kmode, flag_vector, flag_wf, flag_formatted, flag_exit)
!read band_structure_TBA.dat with c_mode = 'rh'
   use parameters, only : pid_energy, t1, t0
   use time
   use mpi_setup
   use print_io
   implicit none
   integer*4              ispin, nspin, nbasis, nband, nk
   integer*4              ispinor
   integer*4              pid, i_continue
   integer*4              ie,ik,is
   integer*4              ii, fi, im, fm
   integer*4              linecount
   integer*4              ne_found(nspin, nk)
   real*8                 E(nband*nspin,nk)
   real*8                 V2(nbasis,nband*nspin,nk)
   real*8                 kpoint
   character*132          inputline
   character*40           desc_str
   character*4            kmode
   external               nitems
   integer*4              nitems, idummy
   logical                flag_go, flag_vector
   logical                flag_wf
   logical                flag_formatted
   logical                flag_exit
   character*80           fnameu, fnamed
   character*11           form_

   call time_check(t1,t0,'init')

   flag_exit = .false.

   pid = pid_energy
   linecount = 0
   ie= 0
   flag_go = .false.
   if(ispin .eq. 2 .and. nspin .eq. 1) then
     ispinor = 2
   else
     ispinor = 1
   endif

   if(nspin .eq. 2) then ! collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA.up.dat'
       fnamed = 'band_structure_TBA.dn.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA.up.bin'
       fnamed = 'band_structure_TBA.dn.bin'
     endif
   elseif(nspin .eq. 1) then ! nm or non-collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA.bin'
     endif
   endif

   do is = 1, nspin
     if(is .eq. 1) then
       if(flag_formatted)       open(pid,file = trim(fnameu),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnameu), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then
         write(message,'(A)')'  !!! ERROR IN READING (up) ', trim(fnameu)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit '  ; write_msg
         flag_exit = .true.
       endif
     elseif(is .eq. 2) then
       if(flag_formatted)       open(pid,file = trim(fnamed),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnamed), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then
         write(message,'(A)')'  !!! ERROR IN READING (dn) ', trim(fnamed)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit '  ; write_msg
         flag_exit = .true.
       endif
     endif


     if(.not. flag_exit) then
       ii = 1+nband*(is-1) ; fi = nband+nband*(is-1)
       im = 1              ; fm = nbasis

       if(flag_formatted) then
         call load_band_structure_V2(E(ii:fi,:), V2(im:fm,ii:fi,:), ne_found(is,:), ispin, nspin, nband, nbasis, nk, kmode, pid, &
                                  flag_vector)
       elseif(.not. flag_formatted) then
         call load_band_structure_V2_bin(E(ii:fi,:), V2(im:fm,ii:fi,:), ne_found(is,:), &
                                      ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector, flag_exit)

       endif
       close(pid)
     endif

   enddo

   call time_check(t1,t0,'end')
   if(nspin .eq. 1) then
     write(message,'(A,F12.6)')"   TIME for LOADING band_structure_TBA (s): ", t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   elseif(nspin .eq. 2) then
     write(message,'(A,F12.6)')"   TIME for LOADING band_structure_TBA.up(dn) (s): ", t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   endif

   return
endsubroutine
subroutine load_band_structure_bin(E, V, ne_found, ispinor, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector, flag_exit)
   use mpi_setup
   use print_io
   implicit none
   integer*4              pid
   integer*4              ie, ik, im
   integer*4              ispin,nband, nbasis, nk, ispinor, nspin
   integer*4              nskip
   integer*4              idummy
   integer*4              nitems
   integer*4              ne_found(nk)
   external               nitems
   character*132          inputline
   real*8                 E(nband, nk)
   real*4                 E4(nband, nk)  ! kind=4
   complex*16             V(nbasis*ispinor,nband,nk)
   complex*8,allocatable::V4(:,:,:) ! kind=4
   logical                flag_go, flag_vector
   character*4            kmode
   character*12           c_emin,c_emax
   real*8                 emin, emax
   real*8, allocatable :: kpoint_(:,:)
   real*4, allocatable :: kpoint4_(:,:) !kind=4
   integer*4              ikmode, nbasis_, nk_, nband_, ispin_, nspin_, ispinor_
   integer*4              init_erange_, fina_erange_
   integer*4              nemax_
   real*8                 emin_, emax_
   logical                flag_vector_, flag_erange_, flag_sparse_
   logical                flag_single_
   logical                flag_exit
   character*2            c_mode_
   flag_exit = .false.

   ne_found = 0
   if(flag_vector) then
     V        = 0d0
   endif

   ! read header
   read(pid) ikmode, flag_vector_, flag_single_, flag_erange_, flag_sparse_, nbasis_, nk_, nband_, &
                     ispin_, nspin_, ispinor_, c_mode_
   ! allocate temporal kpoint_ array which would not be used but (should be) same as PKPTS%kpoint
   allocate(kpoint_(ikmode,nk))

   if(flag_vector .and. flag_single_) then
     allocate(V4(nbasis*ispinor,nband,nk))
     V4       = 0
     allocate(kpoint4_(ikmode,nk))
   endif
  
   ! set E as -999d0 for all (ie,ik) if flag_sparse_ by default. 
   ! if ne_found(ik) is not zero, E(1:ne_found(ik),ik) will be filled in the following step
   ! otherwise, remain as -999d0
   ! this setting is necessary to avoid unexpected DOS which is originated from E(.not. 1:ne_found(ik),ik)
   if(flag_sparse_) then
     E = -999d0
   endif

   if(.not. (c_mode_ .eq. 'wf' .or. c_mode_ .eq. 'no' .or. c_mode_ .eq. 'rh') ) then 
     write(message,'(A)') '    !WARN current version does not support to read C_MODE =/ "wf" or "no" '  ; write_msg
     write(message,'(A)') '     Exit program...'  ; write_msg
     flag_exit = .true.
   endif
 
   if(.not. flag_exit) then

     if(flag_erange_) then
       read(pid) flag_erange_, init_erange_, fina_erange_
     else
       read(pid) flag_erange_
     endif
     if(flag_sparse_) then
       read(pid) flag_sparse_, emin_, emax_, nemax_
     else
       read(pid) flag_sparse_
     endif

   ! read main wavefunction information
     if(flag_vector) then
       if(ispinor .eq. 2) then
         if(c_mode_ .eq. 'wf') then
           if(.not.flag_single_) then
            !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik)), (((V(im,ie,ik),V(im+nbasis,ie,ik)),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk )
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik) 
                do im = 1, nbasis
                  read(pid)  V(im,ie,ik),V(im+nbasis,ie,ik)
                enddo
               enddo
             enddo
           elseif(flag_single_) then
            !read(pid) ((ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik)), (((V4(im,ie,ik),V4(im+nbasis,ie,ik)),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk)
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik)
                do im = 1, nbasis
                  read(pid)  V4(im,ie,ik),V4(im+nbasis,ie,ik)
                enddo
               enddo
             enddo

             kpoint_ = kpoint4_
             E = E4
             V = V4
           endif
         elseif(c_mode_ .eq. 'rh') then ! current version does not support reading 'rh' format but 
                                        ! trying to keep thi line for the future.
                                        ! Since array V is for complex numbers, reading 'rh' data which is 'real' is nonsense.
           write(message,'(A)')'    !ERROR Reading "rh" data is not available in the current version. Check again...'  ; write_msg
           write(message,'(A)')'           Proceed anaway...'  ; write_msg
           if(.not. flag_single_) then
            !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik)), (( V(im,ie,ik),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk)
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik)
                do im = 1, nbasis
                  read(pid)  V(im,ie,ik)
                enddo
               enddo
             enddo

           elseif(flag_single_) then
            !read(pid) ((ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik)), ((V4(im,ie,ik),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk)
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik)
                do im = 1, nbasis
                  read(pid)  V4(im,ie,ik)
                enddo
               enddo
             enddo

             kpoint_ = kpoint4_
             E = E4
             V = V4
           endif
         endif
       elseif(ispinor .eq. 1) then
         if(c_mode_ .eq. 'wf') then
           if(.not. flag_single_) then
            !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik)), ((V(im,ie,ik),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk)
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik)
                do im = 1, nbasis
                  read(pid)  V(im,ie,ik)
                enddo
               enddo
             enddo
             
           elseif(flag_single_) then
            !read(pid) ((ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik)), ((V4(im,ie,ik),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk)
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik)
                do im = 1, nbasis
                  read(pid)  V4(im,ie,ik)
                enddo
               enddo
             enddo
             
             kpoint_ = kpoint4_
             E = E4
             V = V4
           endif
         elseif(c_mode_ .eq. 'rh') then ! current version does not support reading 'rh' format but 
                                        ! trying to keep thi line for the future.
                                        ! Since array V is for complex numbers, reading 'rh' data which is 'real' is nonsense.
           write(message,'(A)')'    !ERROR Reading "rh" data is not available in the current version. Check again...'  ; write_msg
           write(message,'(A)')'           Proceed anaway...'  ; write_msg
           if(.not. flag_single_) then
            !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik)), (( V(im,ie,ik),im=1,nbasis), ie=1,ne_found(ik))),ik=1,nk)
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik)
                do im = 1, nbasis
                  read(pid)  V(im,ie,ik)
                enddo
               enddo
             enddo
             
           elseif(flag_single_) then
            !read(pid) ((ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik)), ((V4(im,ie,ik),im=1,nbasis), ie=1,ne_found(ik))),ik=1,nk)
             do ik = 1, nk
               read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
               do ie = 1, ne_found(ik)
                do im = 1, nbasis
                  read(pid)  V4(im,ie,ik)
                enddo
               enddo
             enddo

             kpoint_ = kpoint4_
             E = E4
             V = V4
           endif
         endif
       endif
     elseif(.not. flag_vector) then
      !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))),ik=1,nk)
       do ik = 1, nk
         read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
       enddo

     endif

   endif ! flag_exit

   deallocate(kpoint_)
   if(allocated(kpoint4_)) deallocate(kpoint_)
   if(allocated(V4))      deallocate(V4)

   return
endsubroutine
subroutine load_band_structure_V(E,V,ne_found, ispinor, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector)
!read band_structure_TBA.dat with c_mode = 'wf' (LORBIT = .TRUE. wf)
   implicit none
   integer*4              pid
   integer*4              ie, ik, im
   integer*4              ispin,nband, nbasis, nk, ispinor, nspin
   integer*4              nskip
   integer*4              idummy
   integer*4              nitems
   integer*4              ne_found(nk)
   external               nitems
   character*132          inputline
   real*8                 E(nband, nk)
   complex*16             V(nbasis*ispinor, nband , nk)
   real*8                 V_(nbasis*ispinor*2)
   logical                flag_go, flag_vector
   character*4            kmode
   character*12           c_emin,c_emax
   real*8                 emin, emax
   real*8, allocatable :: kpoint(:)

   ne_found = 0
   V = 0d0
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
         ne_found(ik) = ne_found(ik) + 1
         backspace(pid)
         if(flag_vector) then
           read(pid,*)kpoint(:), E(ie,ik),V_(:)
         elseif(.not. flag_vector) then
           read(pid,*)kpoint(:), E(ie,ik)
         endif

         do im = 1, nbasis
           if(ispinor .eq. 2) then
             V(im, ie, ik)      = cmplx( V_(im*4-3), V_(im*4-2) )
             V(im+nbasis, ie, ik) = cmplx( V_(im*4-1), V_(im*4) )
           elseif(ispinor .eq. 1) then
             V(im, ie, ik) = cmplx( V_(im*2-1), V_(im*2) )
           endif
         enddo  

       elseif(idummy .eq. 0) then
         backspace(pid)
         read(pid,*)kpoint(:)
         E(ie,ik) = -999d0
         if(flag_vector) V(:,ie,ik) = 0d0
       endif
     enddo
   enddo

   deallocate(kpoint)
   return
endsubroutine

subroutine load_band_structure_V2(E, V2, ne_found, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector)
!read band_structure_TBA.dat with c_mode = 'rh' (LORBIT = .TRUE. rh or LORBIT = .TRUE.)
   implicit none
   integer*4              pid
   integer*4              ie, ik
   integer*4              ispin,nband, nbasis, nk, ispinor, nspin
   integer*4              nskip
   integer*4              idummy
   integer*4              nitems
   integer*4              ne_found(nk)
   external               nitems
   character*132          inputline
   real*8                 E(nband, nk)
   real*8                 V2(nbasis, nband , nk)
   logical                flag_go, flag_vector
   character*4            kmode
   character*12           c_emin,c_emax
   real*8                 emin, emax
   real*8, allocatable :: kpoint(:)

   ne_found = 0
   V2 = 0d0

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
       if(idummy .ge. 1) then
         ne_found(ik) = ne_found(ik) + 1
         backspace(pid)
         if(flag_vector) then
           read(pid,*)kpoint(:), E(ie,ik),V2(1:nbasis,ie,ik)
         elseif(.not. flag_vector) then
           read(pid,*)kpoint(:), E(ie,ik)
         endif
       elseif(idummy .eq. 0) then
         backspace(pid)
         read(pid,*)kpoint(:)
         E(ie,ik) = -999d0
         if(flag_vector) V2(1:nbasis,ie,ik) = 0d0
       endif
     enddo
   enddo

   deallocate(kpoint)

   return
endsubroutine

subroutine load_band_structure_V2_bin(E, V2, ne_found, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector, flag_exit)
   use print_io
   use mpi_setup
   implicit none
   integer*4              pid
   integer*4              ie, ik, im
   integer*4              ispin,nband, nbasis, nk, nspin
   integer*4              nskip
   integer*4              idummy
   integer*4              nitems
   integer*4              ne_found(nk)
   external               nitems
   character*132          inputline
   real*8                 E(nband, nk)
   real*4                 E4(nband, nk)  ! kind4
   real*8                 V2(nbasis,nband,nk)
   real*4, allocatable :: V4(:,:,:) ! kind4
   logical                flag_go, flag_vector
   character*4            kmode
   character*12           c_emin,c_emax
   real*8                 emin, emax
   real*8, allocatable :: kpoint_(:,:)
   real*4, allocatable :: kpoint4_(:,:) !kind4
   integer*4              ikmode, nbasis_, nk_, nband_, ispin_, nspin_, ispinor_
   integer*4              init_erange_, fina_erange_
   integer*4              nemax_
   real*8                 emin_, emax_
   logical                flag_vector_, flag_erange_, flag_sparse_
   logical                flag_single_
   logical                flag_exit
   character*2            c_mode_
   flag_exit = .false.
   ne_found = 0
   if(flag_vector) then
     V2       = 0d0
   endif

   ! read header
   read(pid) ikmode, flag_vector_, flag_single_, flag_erange_, flag_sparse_, nbasis_, nk_, nband_, &
                     ispin_, nspin_, ispinor_, c_mode_

   ! allocate temporal kpoint_ array which would not be used but (should be) same as PKPTS%kpoint
   allocate(kpoint_(ikmode,nk))

   if(flag_vector .and. flag_single_) then
     allocate(V4(nbasis,nband,nk))
     V4       = 0
     allocate(kpoint4_(ikmode,nk))
   endif
  
   ! set E as -999d0 for all (ie,ik) if flag_sparse_ by default. 
   ! if ne_found(ik) is not zero, E(1:ne_found(ik),ik) will be filled in the following step
   ! otherwise, remain as -999d0
   ! this setting is necessary to avoid unexpected DOS which is originated from E(.not. 1:ne_found(ik),ik)
   if(flag_sparse_) then
     E = -999d0
   endif

   if(.not. (c_mode_ .eq. 'wf' .or. c_mode_ .eq. 'no' .or. c_mode_ .eq. 'rh') ) then 
     write(message,'(A)') '    !WARN current version does not support to read C_MODE =/ "wf" or "no" '  ; write_msg
     write(message,'(A)') '     Exit program...'  ; write_msg
     flag_exit = .true.
   endif
 
   if(.not. flag_exit) then

     if(flag_erange_) then
       read(pid) flag_erange_, init_erange_, fina_erange_
     else
       read(pid) flag_erange_
     endif
     if(flag_sparse_) then
       read(pid) flag_sparse_, emin_, emax_, nemax_
     else
       read(pid) flag_sparse_
     endif

   ! read main wavefunction information
     if(flag_vector) then
       if(c_mode_ .eq. 'wf') then ! current routine does not support reading 'wf' format but 
                                  ! trying to keep thi line for the future.
                                  ! Since array V2 is for real numbers, reading 'wf' data which is 'complex' is nonsense.
         if(.not. flag_single_) then
          !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik)), ((V2(im,ie,ik),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk)
           do ik = 1, nk
             read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
             do ie = 1, ne_found(ik)
              do im = 1, nbasis
                read(pid)  V2(im,ie,ik)
              enddo
             enddo
           enddo

         elseif(flag_single_) then
          !read(pid) ((ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik)), ((V4(im,ie,ik),im=1,nbasis),ie=1,ne_found(ik))),ik=1,nk)
           do ik = 1, nk
             read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
             do ie = 1, ne_found(ik)
              do im = 1, nbasis
                read(pid)  V4(im,ie,ik)
              enddo
             enddo
           enddo
           
           kpoint_ = kpoint4_
           E  = E4
           V2 = V4
         endif
       elseif(c_mode_ .eq. 'rh') then
         if(.not. flag_single_) then
          !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik)), (( V2(im,ie,ik),im=1,nbasis), ie=1,ne_found(ik))),ik=1,nk)
           do ik = 1, nk
             read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
             do ie = 1, ne_found(ik)
              do im = 1, nbasis
                read(pid)  V2(im,ie,ik)
              enddo
             enddo
           enddo
          
         elseif(flag_single_) then
          !read(pid) ((ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik)), ((V4(im,ie,ik),im=1,nbasis), ie=1,ne_found(ik))),ik=1,nk)
           do ik = 1, nk
             read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
             do ie = 1, ne_found(ik)
              do im = 1, nbasis
                read(pid)  V4(im,ie,ik)
              enddo
             enddo
           enddo
           
           kpoint_ = kpoint4_
           E  = E4
           V2 = V4
         endif
       endif
     elseif(.not. flag_vector) then
      !read(pid) ((ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))),ik=1,nk)
       do ik = 1, nk
         read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
       enddo
       
     endif

   endif ! flag_exit

   deallocate(kpoint_)
   if(allocated(kpoint4_)) deallocate(kpoint4_)
   if(allocated(V4)) deallocate(V4)

   return
endsubroutine
