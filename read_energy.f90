#include "alias.inc"
subroutine set_target_and_weight(EDFT, PWGHT, PINPT, PGEOM, PKPTS)
   use parameters, only : energy, incar, weight, poscar, kpoints
   use mpi_setup
   use print_io
   implicit none
   type(incar)   :: PINPT
   type(energy)  :: EDFT, EDFT_ALL
   type(weight)  :: PWGHT
   type(poscar)  :: PGEOM
   type(kpoints) :: PKPTS
   integer*4        mpierr
   integer*4        i
   logical          flag_efile_exist, flag_efileu_exist, flag_efiled_exist
   character*80     fname

   if(PINPT%flag_tbfit_finish) return   
   if(.not. PINPT%flag_tbfit) return

   ! read info: target energy to be fitted with
   if(PINPT%flag_tbfit) then
     if(PINPT%flag_collinear .and. PWGHT%efile_type .eq. 'user') then
       if(len_trim(PWGHT%efilenmu) .eq. 0 .or. len_trim(PWGHT%efilenmd) .eq. 0) then
         write(message,'(A)')'  !WARN!  EFILE has not been set properly. Check EFILE or EFILEU, EFILED. Exit program.' ; write_msg
         kill_job
       endif
       inquire(file=PWGHT%efilenmu,exist=flag_efileu_exist)
       inquire(file=PWGHT%efilenmd,exist=flag_efiled_exist)
     else
       inquire(file=PWGHT%efilenmu,exist=flag_efileu_exist)
     endif
   endif

   if(PINPT%flag_collinear .and. PINPT%flag_tbfit) then
     if(flag_efileu_exist .and. flag_efiled_exist .and. PWGHT%efile_type .eq. 'user') then
       flag_efile_exist = .true.
     elseif(flag_efileu_exist .and. PWGHT%efile_type .eq. 'vasp' ) then
       flag_efile_exist = .true.
     else
       flag_efile_exist = .false.
     endif
   elseif(.not. PINPT%flag_collinear .and. PINPT%flag_tbfit) then
     if(flag_efileu_exist) then
       flag_efile_exist = .true.
     else
       flag_efile_exist = .false.
     endif
   endif
 
   if(PINPT%flag_tbfit .and. flag_efile_exist) then
     
     if(PWGHT%flag_weight_default) then
       allocate(PWGHT%WT(PGEOM%neig*PINPT%ispin, PKPTS%nkpoint))
       PWGHT%WT(:,:)=0.00001d0 !initialize
       if(PINPT%flag_fit_degeneracy) allocate(PWGHT%DEGENERACY_WT(PGEOM%neig*PINPT%ispin, PKPTS%nkpoint))
       if(PINPT%flag_fit_degeneracy) PWGHT%DEGENERACY_WT(:,:)=0d0       !initialize
 
 
       if(PWGHT%efile_type .eq. 'vasp') then
         call read_energy_vasp(PINPT,PGEOM,PKPTS,EDFT,EDFT_all,PWGHT)
       else
         call read_energy(PINPT,PGEOM,PKPTS,EDFT,EDFT_all,PWGHT)
       endif
       
       write(message,*)' '  ; write_msg
       write(message,*)'#-START SET WEIGHT FOR FITTING  --------------------'; write_msg
       if(PINPT%flag_fit_degeneracy) then
         call get_degeneracy(EDFT, PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint, PINPT)
         do i=1, PWGHT%ndegenw
           PWGHT%max_len_strip_kp =max(len_trim(PWGHT%strip_kp_deg(i)), PWGHT%max_len_strip_kp)
           PWGHT%max_len_strip_tb =max(len_trim(PWGHT%strip_tb_deg(i)), PWGHT%max_len_strip_tb)
           PWGHT%max_len_strip_df =max(len_trim(PWGHT%strip_df_deg(i)), PWGHT%max_len_strip_df)
           call set_weight(PINPT, PGEOM, PKPTS, PWGHT, EDFT, EDFT_all, PWGHT%strip_kp_deg(i), PWGHT%strip_tb_deg(i), &
                           PWGHT%strip_df_deg(i), PWGHT%strip_wt_deg(i), PWGHT%DEGENERACY_WT, 2)
         enddo
       endif
       if(PINPT%flag_print_only_target ) then
         fname = 'band_structure_DFT'//trim(PINPT%title(EDFT%mysystem))//'.dat'
         if_main call print_energy_weight( PKPTS%kpoint, PKPTS%nkpoint, EDFT, PWGHT, PGEOM%neig, PINPT, &
                                           fname,PINPT%flag_get_band_order)
         write(message,'(A    )')'  !WARN! PRINT_ONLY_TARGET requested..'         ; write_msg
         write(message,'(A,A,A)')'  !WARN! check ', trim(fname), ' Exit..'  ; write_msg
         kill_job
       endif
 
     elseif(.not. PWGHT%flag_weight_default) then
       allocate(PWGHT%WT(PGEOM%neig*PINPT%ispin, PKPTS%nkpoint))
       if(PINPT%flag_fit_degeneracy) allocate( PWGHT%DEGENERACY_WT(PGEOM%neig*PINPT%ispin, PKPTS%nkpoint) )
 
       if(PWGHT%efile_type .eq. 'vasp') then
         call read_energy_vasp(PINPT,PGEOM,PKPTS,EDFT,EDFT_all,PWGHT)
       else
         call read_energy(PINPT,PGEOM,PKPTS,EDFT,EDFT_all,PWGHT)
       endif
 
       write(message,*)' '  ; write_msg
       write(message,*)'#-START SET WEIGHT FOR FITTING  --------------------'; write_msg
       PWGHT%WT(:,:)=0.00001d0 !initialize
       if(PINPT%flag_fit_degeneracy) PWGHT%DEGENERACY_WT(:,:)=0d0       !initialize
 
       do i=1, PWGHT%nweight
         PWGHT%max_len_strip_kp =max(len_trim(PWGHT%strip_kp(i)), PWGHT%max_len_strip_kp)
         PWGHT%max_len_strip_tb =max(len_trim(PWGHT%strip_tb(i)), PWGHT%max_len_strip_tb)
         PWGHT%max_len_strip_df =max(len_trim(PWGHT%strip_df(i)), PWGHT%max_len_strip_df)
       enddo
       do i=1, PWGHT%nweight
         call set_weight(PINPT, PGEOM, PKPTS, PWGHT, EDFT, EDFT_all, PWGHT%strip_kp(i), PWGHT%strip_tb(i), &
                         PWGHT%strip_df(i), PWGHT%strip_wt(i), PWGHT%WT, 1)
       enddo
       
       if(PINPT%flag_fit_degeneracy) then
       ! get degeneracy information for DFT target
         call get_degeneracy(EDFT, PGEOM%neig*PINPT%ispin, PKPTS%nkpoint, PINPT)
         do i=1, PWGHT%ndegenw
           PWGHT%max_len_strip_kp =max(len_trim(PWGHT%strip_kp_deg(i)), PWGHT%max_len_strip_kp)
           PWGHT%max_len_strip_tb =max(len_trim(PWGHT%strip_tb_deg(i)), PWGHT%max_len_strip_tb)
           PWGHT%max_len_strip_df =max(len_trim(PWGHT%strip_df_deg(i)), PWGHT%max_len_strip_df)
         enddo
         do i=1, PWGHT%ndegenw
           call set_weight(PINPT, PGEOM, PKPTS, PWGHT, EDFT, EDFT_all, PWGHT%strip_kp_deg(i), PWGHT%strip_tb_deg(i), &
                           PWGHT%strip_df_deg(i), PWGHT%strip_wt_deg(i), PWGHT%DEGENERACY_WT, 2)
         enddo
       endif
 
       ! normalize weight so that their sum to be 1
       !PWGHT%WT = PWGHT%WT / sum(PWGHT%WT)
 
       if(PINPT%flag_print_only_target ) then
         fname = 'band_structure_DFT'//trim(PINPT%title(EDFT%mysystem))//'.dat'
         if_main call print_energy_weight( PKPTS%kpoint, PKPTS%nkpoint, EDFT, PWGHT, PGEOM%neig, PINPT, &
                                           fname,PINPT%flag_get_band_order)
 
         write(message,'(A    )')'  !WARN! PRINT_ONLY_TARGET requested..'        ; write_msg
         write(message,'(A,A,A)')'  !WARN! check ', trim(fname), ' Exit..' ; write_msg
         kill_job
       endif
 
     endif
 
     if(.not. PWGHT%flag_weight_orb) then
       allocate(PWGHT%PENALTY_ORB(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin, PKPTS%nkpoint))
       PWGHT%PENALTY_ORB(:,:,:) = 0.d0 ! initialize & set default value (zeros for all)
     elseif(PWGHT%flag_weight_orb) then
       allocate(PWGHT%PENALTY_ORB(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin, PKPTS%nkpoint))
       PWGHT%PENALTY_ORB(:,:,:) = 0.d0 ! initialize
       do i = 1, PWGHT%npenalty_orb
         PWGHT%max_len_strip_kp =max(len_trim(PWGHT%strip_kp_orb(i)), PWGHT%max_len_strip_kp)
         PWGHT%max_len_strip_tb =max(len_trim(PWGHT%strip_tb_orb(i)), PWGHT%max_len_strip_tb)
         PWGHT%max_len_strip_ob =max(len_trim(PWGHT%strip_pen_orb(i)),PWGHT%max_len_strip_ob)
         PWGHT%max_len_strip_st =max(len_trim(PWGHT%strip_site(i)),   PWGHT%max_len_strip_st)
       enddo
       do i = 1, PWGHT%npenalty_orb
         call set_penalty_orb(PINPT, PGEOM, PKPTS, PWGHT, PWGHT%strip_kp_orb(i), PWGHT%strip_tb_orb(i), &
                              PWGHT%strip_orb(i), PWGHT%strip_site(i), PWGHT%strip_pen_orb(i) )
       enddo
     endif
 
     write(message,*)'#- END SET WEIGHT FOR FITTING --------------------'  ; write_msg
     write(message,*)' '  ; write_msg

   elseif(PINPT%flag_tbfit .and. .not. flag_efile_exist ) then
     write(message,'(A,A,A)')'  !WARN!  The target file ',trim(PWGHT%efilenmu),' does not exist or not specified!! '; write_msg
     write(message,'(A    )')'          Please check EFILE tag again. Exit...' ; write_msg
     kill_job
   endif
 
   return
endsubroutine
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
  integer*4        idummy, ispin, nkp_target, neig_target, nelect_vasp, ie_start
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
  fname = trim(PWGHT%efilenmu)

  allocate( EDFT_%E(max_eig, PKPTS%nkpoint) )
  allocate( EDFT__%E(max_eig, PKPTS%nkpoint) )
  allocate( occ(max_eig, PKPTS%nkpoint) )
  allocate( occ_(max_eig, PKPTS%nkpoint) )
  open(pid_energy,file=trim(fname), status='old', iostat=i_continue)

  write(message,*)' '  ; write_msg
  write(message,*)'#- READING TARGET ENERGY FILE: ',trim(fname)  ; write_msg

  ie_start=PWGHT%itarget_e_start
  line_tot = 0

  ! reading header for EIGENVAL
  read(pid_energy,*) idummy, idummy, idummy, ispin
  do k=1, 4
    read(pid_energy, *) dummy
  enddo
  read(pid_energy,*) nelect_vasp, nkp_target, neig_target
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

  EDFT__%E = EDFT__%E - PWGHT%efile_ef

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
    if( (neig_target - PWGHT%itarget_e_start + 1) .ge. size_range ) then
      EDFT_%E(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PWGHT%itarget_e_start:neig_target,1:nkp_target)
      EDFT_%E(1+neig_target-PWGHT%itarget_e_start+1: (neig_target-PWGHT%itarget_e_start+1)*2,1:nkp_target) = EDFT__%E(1+neig_target:neig_target*2,1:nkp_target)
      occ(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = occ_(PWGHT%itarget_e_start:neig_target,1:nkp_target)
      occ(1+neig_target-PWGHT%itarget_e_start+1: (neig_target-PWGHT%itarget_e_start+1)*2,1:nkp_target) = occ_(1+neig_target:neig_target*2,1:nkp_target)
    else ! if neig_target - PWGHT%itarget_e_start + 1 is less than fband - iband + 1
      idummy = size_range - (neig_target - PWGHT%itarget_e_start + 1)
      EDFT_%E(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PWGHT%itarget_e_start:neig_target,1:nkp_target)
      EDFT_%E(neig_target-PWGHT%itarget_e_start+1+1:size_range,1:nkp_target) = 100d0 ! dummy value is added
      EDFT_%E(size_range + 1: size_range + neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = EDFT__%E(1+neig_target:neig_target*2,1:nkp_target)
      EDFT_%E(size_range + neig_target-PWGHT%itarget_e_start+1+1:size_range*2,1:nkp_target) = 100d0

      occ(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target)                            = occ_(PWGHT%itarget_e_start:neig_target,1:nkp_target)
      occ(neig_target-PWGHT%itarget_e_start+1+1:size_range,1:nkp_target)                 = 100d0 ! dummy value is added
      occ(size_range + 1: size_range + neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = occ_(1+neig_target:neig_target*2,1:nkp_target)
      occ(size_range + neig_target-PWGHT%itarget_e_start+1+1:size_range*2,1:nkp_target)  = 100d0
    endif
  else
    if( (neig_target - PWGHT%itarget_e_start + 1) .ge. size_range ) then
      EDFT_%E(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PWGHT%itarget_e_start:neig_target,1:nkp_target)
      occ(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = occ_(PWGHT%itarget_e_start:neig_target,1:nkp_target)
    else
      EDFT_%E(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target) = EDFT__%E(PWGHT%itarget_e_start:neig_target,1:nkp_target)
      EDFT_%E(neig_target-PWGHT%itarget_e_start+1+1:size_range,1:nkp_target) = 100d0
 
      occ(1:neig_target-PWGHT%itarget_e_start+1,1:nkp_target)            = occ_(PWGHT%itarget_e_start:neig_target,1:nkp_target)
      occ(neig_target-PWGHT%itarget_e_start+1+1:size_range,1:nkp_target) = 100d0
    endif
  endif

  if(neig_target - PWGHT%itarget_e_start + 1 .ge. size_range) then
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

  write(message,*)'#- END READING TARGET ENERGY FILE --------------'  ; write_msg
  write(message,*)' '  ; write_msg
  close(pid_energy)
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
  type(incar)            :: PINPT
  type(poscar)           :: PGEOM
  type(kpoints)          :: PKPTS
  type(energy)           :: EDFT_
  type(energy)           :: EDFT_all
  type(energy)           :: EDFT
  type(weight)           :: PWGHT
  character*132             fnameu,fnamed
  character*132             inputline
  character*40              desc_str,dummy, dummy_(1000)
  character(*), parameter:: func = 'read_energy'
  integer*4, parameter   :: max_eig = 100000
  integer*4                 ie, ik, k, i_continue
  integer*4                 size_range
  integer*4                 mpierr
  integer*4                 ik_, kreduce
  integer*4, allocatable :: irange(:),irange_up(:), irange_dn(:)
  logical          flag_order, flag_fit_orbital
  integer*4                 vasp_orb_order(PINPT%lmmax)
  real*8,    allocatable :: TOT(:,:)

  flag_order = PINPT%flag_get_band_order 
  flag_fit_orbital = PINPT%flag_fit_orbital ! this is only and should available with spd orbital resolved data file
                                            ! one can obtain such data, for example via VASP PROCAR file
                                            ! The target DFT file format should be as follows
                                            ! 1st-column    2nd-column  3rd  4th  5th  6th  7th  8th  9th  10th 11th 12th  13th  14th  15th 
                                            !  k-path(k1)    energy(e1)  s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k2            e1   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k3            e1   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k4            e1   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k5            e1   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k6            .....
                                            !         k7            .....
                                            !         k8            .....

                                            !         k1            e2   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k2            e2   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k3            e2   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k4            .....
                                            !         k5            .....
                                            !         k6            .....
                                            
                                            !         k1            e3   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k2            e3   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k3            e3   s    p    d   s    py    pz   px  dxy  dyz   dz2  dxz    dx2   tot
                                            !         k4            .....
                                            !         k5            .....
                                            !         k6            .....

                                            ! NOTE: Here, the data should be normalized so that the "tot" to be unity after reading the data.
                                            
  ! NOTE: kreduce with greater than 1 is only applicaple with PINPT%
  if(PKPTS%kreduce .gt. 1) then
    kreduce = PKPTS%kreduce
    if( trim(PKPTS%kline_type) .eq. 'vasp' .or. trim(PKPTS%kline_type) .eq. 'VASP') then
      write(message,'(A)')'  !WARNING!  KREDUCE should be 1 if you read VASP type of target energy file.'  ; write_msg
      write(message,'(A)')'             Please read user manual with KREDUCE tag and KFILE tag.'  ; write_msg
      write(message,'(A)')'             Exit program...'  ; write_msg
      kill_job
    endif
  else
    kreduce = 1
  endif

  if(PINPT%flag_collinear) then
   fnameu = trim(PWGHT%efilenmu)
   fnamed = trim(PWGHT%efilenmd)
   allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(flag_order) allocate( EDFT%E_ORD(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) ) 
   if(PINPT%flag_fit_degeneracy) then
     allocate( EDFT%D(3,PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   endif
   if(flag_fit_orbital) allocate( EDFT%ORB(PINPT%lmmax,PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(PGEOM%neig*2 .gt. max_eig) then
     write(message,'(A)')'  !WARNING!  NEIG is exceeding predefined max_eig=1000000. Please increase "max_eig"'  ; write_msg
     write(message,'(A,A)')'  !WARNING!  Exit program... ',func  ; write_msg
     stop
   endif
  elseif(PINPT%flag_noncollinear) then
   fnameu = trim(PWGHT%efilenmu)
   allocate( EDFT%E(PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
   if(flag_fit_orbital) allocate( EDFT%ORB(PINPT%lmmax,PGEOM%neig*PINPT%ispin,  PKPTS%nkpoint) )
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
   fnameu = trim(PWGHT%efilenmu)
   allocate( EDFT%E(PGEOM%neig,  PKPTS%nkpoint) )
   if(flag_fit_orbital) allocate( EDFT%ORB(PINPT%lmmax,PGEOM%neig,  PKPTS%nkpoint) )
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
  if(flag_fit_orbital) then 
    allocate( EDFT_%ORB(PINPT%lmmax,max_eig,PKPTS%nkpoint))
    allocate(       TOT(            max_eig,PKPTS%nkpoint)) ! temporal, for the normalization purpose
    TOT = 0d0
    if(PINPT%lmmax .eq. 9) then
      vasp_orb_order=(/1,3,4,2,7,9,5,8,6/)
    elseif(PINPT%lmmax .eq. 3) then
      vasp_orb_order=(/1,2,3/)
    endif
  endif
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
    write(message,*)'#- READING TARGET ENERGY FILE: ',trim(fnameu), ' and ', trim(fnamed)  ; write_msg
  else
    write(message,*)'#- READING TARGET ENERGY FILE: ',trim(fnameu)  ; write_msg
  endif

  ie=0
  !loop1 for spin up 
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
             if(PWGHT%read_energy_column_index .eq. 2) then
               if(.not. flag_order) then
                 if(flag_fit_orbital) then
                   if(PINPT%lmmax .eq.9) then
                     read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik), dummy,dummy,dummy,EDFT_%ORB(vasp_orb_order,ie,ik)
                   elseif(PINPT%lmmax .eq. 3) then
                     read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik),                   EDFT_%ORB(vasp_orb_order,ie,ik)
                   endif
                     TOT(ie,ik) = sum(EDFT_%ORB(vasp_orb_order,ie,ik))! NOTE: not work for f-electron case
                 else
                   read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik)
                 endif
               elseif(flag_order) then 
                 ! this is for read re-ordered target band and save it to E_ORD, 
                 ! hence, be sure that file should contain this information too!
                 read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik), EDFT_%E_ORD(ie, ik)
               endif
             elseif(PWGHT%read_energy_column_index .eq. 3) then
               if(.not. flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik)
               elseif(flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, dummy, EDFT_%E(ie, ik), EDFT_%E_ORD(ie, ik)
               endif
             elseif(PWGHT%read_energy_column_index .eq. 4) then
               read(inputline,*,iostat=i_continue) dummy, dummy, dummy, EDFT_%E(ie, ik)
             elseif(PWGHT%read_energy_column_index .eq. 5) then
               read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
             elseif(PWGHT%read_energy_column_index .eq. 6) then
               read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
             elseif(PWGHT%read_energy_column_index .eq. 7) then
               read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
             elseif(PWGHT%read_energy_column_index .eq. 8) then
               read(inputline,*,iostat=i_continue) dummy, dummy, dummy, dummy, dummy, dummy, dummy, EDFT_%E(ie, ik)
             elseif(PWGHT%read_energy_column_index .ge. 9) then
               read(inputline,*,iostat=i_continue) dummy_(1:PWGHT%read_energy_column_index-1), EDFT_%E(ie, ik)
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
         !loop2 for spin dn
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
             if(PWGHT%read_energy_column_index_dn .eq. 2) then
               if(.not.flag_order) then
                 if(flag_fit_orbital) then
                   if(PINPT%lmmax .eq. 9) then
                     read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik), dummy,dummy,dummy,EDFT_%ORB(vasp_orb_order,ie,ik)
                   elseif(PINPT%lmmax .eq. 3) then
                     read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik),                   EDFT_%ORB(vasp_orb_order,ie,ik)             ! NOTE: not work for f-electron case
                   endif
                   TOT(ie,ik) = sum(EDFT_%ORB(vasp_orb_order,ie,ik))
                 else
                   read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik)
                 endif
               elseif(flag_order) then
                 read(inputline,*,iostat=i_continue) dummy, EDFT_%E(ie, ik), EDFT_%E_ORD(ie, ik)
               endif
             elseif(PWGHT%read_energy_column_index_dn .eq. 3) then
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
       endif ! flag_collinear

       PGEOM%neig_target=ie - 1
       write(message,'(A,I8)')' N_TARGET:',PGEOM%neig_target  ; write_msg

  allocate( EDFT_all%E(PGEOM%neig_target, PKPTS%nkpoint) )
  if(flag_fit_orbital) allocate( EDFT_all%ORB(PINPT%lmmax, PGEOM%neig_target, PKPTS%nkpoint))
  if(flag_order) allocate( EDFT_all%E_ORD(PGEOM%neig_target, PKPTS%nkpoint) )

  if(PINPT%flag_collinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig
  elseif(PINPT%flag_noncollinear) then
    if(PWGHT%fband .eq. -9999) PWGHT%fband = PGEOM%neig*2 !(multiplying ispinor)
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

  EDFT_%E = EDFT_%E - PWGHT%efile_ef
  if(flag_order) EDFT_%E_ORD = EDFT_%E_ORD - PWGHT%efile_ef

  if(PINPT%flag_collinear) then
    if(PINPT%flag_scissor) then ! setup scissor
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

    ! store target data to EDFT%E
    EDFT%E(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
    if(flag_fit_orbital) EDFT%ORB(:,1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%ORB(:,irange,1:PKPTS%nkpoint)
    if(flag_order) EDFT%E_ORD(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E_ORD(irange,1:PKPTS%nkpoint)

  elseif(PINPT%flag_noncollinear) then
    if(PINPT%flag_scissor) then  ! setup scissor
      EDFT_%E(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
      EDFT_%E(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
      if(flag_order) then
        EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) = &
        EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig*2,1:PKPTS%nkpoint) + PINPT%r_scissor
      endif
    endif

    ! store target data to EDFT%E
    EDFT%E(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
    if(flag_order) EDFT%E_ORD(1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%E_ORD(irange,1:PKPTS%nkpoint)
    if(flag_fit_orbital) EDFT%ORB(:,1:PGEOM%neig*PINPT%ispin,1:PKPTS%nkpoint) = EDFT_%ORB(:,irange,1:PKPTS%nkpoint)

  else ! if nonmagnetic
    if(PINPT%flag_scissor) then  ! setup scissor
      EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
      EDFT_%E(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
      if(flag_order) then
        EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) = &
        EDFT_%E_ORD(PINPT%i_scissor:PGEOM%neig,1:PKPTS%nkpoint) + PINPT%r_scissor
      endif
    endif

    ! store target data to EDFT%E
    EDFT%E(1:PGEOM%neig,1:PKPTS%nkpoint) = EDFT_%E(irange,1:PKPTS%nkpoint)
    if(flag_order) EDFT%E_ORD(1:PGEOM%neig,1:PKPTS%nkpoint) = EDFT_%E_ORD(irange,1:PKPTS%nkpoint)
    if(flag_fit_orbital) then 
      do ik = 1, PKPTS%nkpoint
        do ie = 1, PGEOM%neig
          EDFT%ORB(:,ie,ik) = EDFT_%ORB(:,irange(ie),ik)/TOT(irange(ie),ik)
        enddo
      enddo
    endif

  endif

  EDFT_all%E(1:PGEOM%neig_target,1:PKPTS%nkpoint) = EDFT_%E(1:PGEOM%neig_target,1:PKPTS%nkpoint)
  if(flag_fit_orbital) EDFT_all%ORB(:,1:PGEOM%neig_target,1:PKPTS%nkpoint) = EDFT_%ORB(:,1:PGEOM%neig_target,1:PKPTS%nkpoint)
  if(flag_order) EDFT_all%E_ORD(1:PGEOM%neig_target,1:PKPTS%nkpoint) = EDFT_%E_ORD(1:PGEOM%neig_target,1:PKPTS%nkpoint)

  write(message,*)'#- END READING TARGET ENERGY FILE --------------'  ; write_msg
  write(message,*)' '  ; write_msg
  close(pid_energy)
  if(PINPT%flag_collinear) close(pid_energy+1)
  deallocate(EDFT_%E)
  if(flag_order) deallocate(EDFT_%E_ORD)
return
endsubroutine

subroutine read_energy_tbfit_V(PRPLT, PINPT, E, V, ne_found, ispin, nspin, nbasis, nband, nk, kmode, flag_vector, flag_wf, flag_formatted, flag_exit)
! read band_structure_TBA.xxx generated by TBFIT.
   use parameters, only : pid_energy, t1, t0, replot, incar
   use time
   use print_io
   use mpi_setup
   implicit none
   type(incar)         :: PINPT
   type(replot)        :: PRPLT
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
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.dat'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.bin'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.bin'
     endif
   elseif(nspin .eq. 1) then ! nm or non-collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.bin'
     endif
   endif
  

   do is = 1, nspin
     if(is .eq. 1) then
       if(flag_formatted)       open(pid,file = trim(fnameu),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnameu), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then 
         write(message,'(A)')'  !!! ERROR IN READING (up) ', trim(fnameu)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit_V '  ; write_msg
         flag_exit = .true.
       endif
     elseif(is .eq. 2) then
       if(flag_formatted)       open(pid,file = trim(fnamed),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnamed), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then 
         write(message,'(A)')'  !!! ERROR IN READING (dn) ', trim(fnamed)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit_V '  ; write_msg
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
         call load_band_structure_V_bin(E(ii:fi,:), V(im:fm,ii:fi,:), ne_found(is,:), &
                                      ispinor, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector, flag_exit)
       endif   
      
       close(pid)
     endif
   enddo 

   call time_check(t1,t0,'end')
   if(nspin .eq. 1) then
     write(message,'(3A,F12.6)')"   TIME for LOADING (s) ",trim(fnameu)," : ", t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   elseif(nspin .eq. 2) then
     write(message,'(5A,F12.6)')"   TIME for LOADING (s) ",trim(fnameu)," and ", trim(fnamed), " : " , t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   endif

   return
endsubroutine
subroutine read_energy_tbfit_V2(PRPLT, PINPT, E, V2, ne_found, ispin, nspin, nbasis, nband, nk, kmode, flag_vector, flag_wf, flag_formatted, flag_exit)
!read band_structure_TBA.xxx with c_mode = 'rh'
   use parameters, only : pid_energy, t1, t0, replot, incar
   use time
   use mpi_setup
   use print_io
   implicit none
   type(incar)         :: PINPT
   type(replot)        :: PRPLT
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
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.dat'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.bin'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.bin'
     endif
   elseif(nspin .eq. 1) then ! nm or non-collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.bin'
     endif
   endif

   do is = 1, nspin
     if(is .eq. 1) then
       if(flag_formatted)       open(pid,file = trim(fnameu),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnameu), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then
         write(message,'(A)')'  !!! ERROR IN READING (up) ', trim(fnameu)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit_V2 '  ; write_msg
         flag_exit = .true.
       endif
     elseif(is .eq. 2) then
       if(flag_formatted)       open(pid,file = trim(fnamed),                    status='old', iostat=i_continue)
       if(.not. flag_formatted) open(pid,file = trim(fnamed), form='unformatted',status='old', iostat=i_continue)

       if(i_continue .ne. 0) then
         write(message,'(A)')'  !!! ERROR IN READING (dn) ', trim(fnamed)  ; write_msg
         write(message,'(A)')'      EXIT PROGRAM : read_energy_tbfit_V2 '  ; write_msg
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
     write(message,'(3A,F12.6)')"   TIME for LOADING (s) ",trim(fnameu)," : ", t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   elseif(nspin .eq. 2) then
     write(message,'(5A,F12.6)')"   TIME for LOADING (s) ",trim(fnameu)," and ", trim(fnamed), " : " , t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   endif

   return
endsubroutine
subroutine read_energy_tbfit(PRPLT, PINPT, E, ne_found, ispin, nspin, nbasis, nband, nk, kmode, flag_formatted, flag_exit)
   !read band_structure_TBA.xxx with c_mode = 'no'
   use parameters, only : pid_energy, t1, t0, replot, incar
   use time
   use mpi_setup
   use print_io
   implicit none
   type(incar)         :: PINPT
   type(replot)        :: PRPLT
   integer*4              ispin, nspin, nbasis, nband, nk
   integer*4              ispinor
   integer*4              pid, i_continue
   integer*4              ie,ik,is
   integer*4              ii, fi, im, fm
   integer*4              linecount
   integer*4              ne_found(nspin, nk)
   real*8                 E(nband*nspin,nk)
   real*8                 kpoint
   character*132          inputline
   character*40           desc_str
   character*4            kmode
   external               nitems
   integer*4              nitems, idummy
   logical                flag_go
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
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.dat'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.bin'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.bin'
     endif
   elseif(nspin .eq. 1) then ! nm or non-collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dat'
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.bin'
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

       if(flag_formatted) then
         call load_band_structure(E(ii:fi,:), ne_found(is,:), ispin, nspin, nband, nbasis, nk, kmode, pid)
       elseif(.not. flag_formatted) then
         call load_band_structure_bin(E(ii:fi,:), ne_found(is,:), ispin, nspin, nband, nbasis, nk, kmode, pid, flag_exit)

       endif
       close(pid)
     endif

   enddo

   call time_check(t1,t0,'end')
   if(nspin .eq. 1) then
     write(message,'(3A,F12.6)')"   TIME for LOADING (s) ",trim(fnameu)," : ", t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   elseif(nspin .eq. 2) then
     write(message,'(5A,F12.6)')"   TIME for LOADING (s) ",trim(fnameu)," and ", trim(fnamed), " : " , t1  ; write_msg
     write(message,'(A)')''  ; write_msg
   endif

   return
endsubroutine
subroutine load_band_structure_V_bin(E, V, ne_found, ispinor, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector, flag_exit)
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
!read band_structure_TBA.xxx with c_mode = 'wf' (LORBIT = .TRUE. wf)
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
subroutine load_band_structure(E,ne_found, ispin, nspin, nband, nbasis, nk, kmode, pid)
!read band_structure_TBA.xxx with c_mode = 'no'
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
   logical                flag_go
   character*4            kmode
   character*12           c_emin,c_emax
   real*8                 emin, emax
   real*8, allocatable :: kpoint(:)

   ne_found = 0
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
         read(pid,*)kpoint(:), E(ie,ik)

       elseif(idummy .eq. 0) then
         backspace(pid)
         read(pid,*)kpoint(:)
         E(ie,ik) = -999d0
       endif
     enddo
   enddo

   deallocate(kpoint)
   return
endsubroutine

subroutine load_band_structure_V2(E, V2, ne_found, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_vector)
!read band_structure_TBA.xxt with c_mode = 'rh' (LORBIT = .TRUE. rh or LORBIT = .TRUE.)
   implicit none
   integer*4              pid
   integer*4              ie, ik
   integer*4              ispin,nband, nbasis, nk, nspin
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

! NOTE: 1. This routine is not spin dependent. In fname, ".up.dat/.up.bin or .dn.dat/.dn.bin" is specified
! and read through "fname" and store E and V, C accordingly.
! 2. It is necesarry to make the routine toward flag_overlap case if c_mode_ == 'wf'.
! For this, one need to load Overlap matrix multiplied wavevector S|psi> should be read together. 
! This is the future work. In this case, the orbital projection information stored in "V" variable,
! can be evaluated from wavefunction coefficient C and overlap matrix S multiplied C, SC. => V = C' * SC
! 3. For c_mode_ == 'rh' and flag_overlap case, the stored orbital projection is Mulliken charge, C'*SC
! HJ Kim. 18. March 2021
subroutine load_band_singlek(PGEOM, fname, nemax, nbasis, ispinor, E, C, V, ne_found, kp, gg, &
                             flag_get_ldos, flag_wf, flag_sparse, flag_phase, flag_formatted)
    use parameters, only : poscar
    use print_io
    use mpi_setup
    use phase_factor
    implicit none
    type(poscar) :: PGEOM
    character*80    fname
    character*256   inputline
    integer*4       my_pid, idummy
    integer*4       mpierr, i_continue
    integer*4,external :: nitems
    integer*4       nemax, nbasis, ispinor, ne_found
    integer*4       ie, im
    real*8          kp(3), gg(3)
    real*4          kp4(3)
    logical         flag_go, flag_get_ldos, flag_sparse, flag_formatted
    logical         flag_phase, flag_wf
    integer*4       ikmode, nbasis_, nk_, nband_, ispin_, nspin_, ispinor_
    integer*4       init_erange_, fina_erange_, nemax_
    real*8          emin_, emax_
    logical         flag_vector_, flag_erange_, flag_sparse_, flag_single_
    character*2     c_mode_
    real*8          E(nemax)
    real*4          E4(nemax)
    real*8          C_(nbasis*ispinor * 2)
    real*8          V(nbasis,nemax) ! Be aware that V here is orbital projection
    real*4          V4(nbasis,nemax) ! Be aware that V here is orbital projection
    complex*16      C(nbasis*ispinor,nemax) ! Here C is the wavefunction coefficient (wave vector) multiplied by phase 
    complex*16      c_up, c_dn
    complex*8       c4_up, c4_dn
    complex*16      phase

    ne_found  = 0
    E = -999d0 ; V = 0d0 ; C = 0d0;
    my_pid = 100 + myid
    phase = (1d0,0d0)

    if(flag_wf) C = 0d0
    if(flag_formatted) then
      open(my_pid,file=trim(fname), status='old', iostat=i_continue)
      if(flag_sparse) then
        read(my_pid, '(A)') inputline
        read(my_pid, '(A)') inputline
        read(my_pid, *    ) inputline, inputline, ne_found
      else
        read(my_pid, '(A)') inputline
        ne_found = nemax
      endif

      if(ne_found .eq. 0) then
        flag_go = .false. ! check whether blank is found : if found, it means that a new eigen block is found
        do while (.not.flag_go)
          read(my_pid,'(A)')inputline
          idummy = nitems(inputline)
          if(idummy .le. 0) then
            flag_go = .false.
          elseif(idummy .gt. 0) then
            flag_go = .true.
            backspace(my_pid)
          endif
        enddo
        read(my_pid,*) kp(:)
        close(my_pid)
        return
      elseif(ne_found .gt. 0) then
    
        do ie = 1, ne_found
          flag_go = .false. 
          do while (.not.flag_go)
            read(my_pid,'(A)')inputline
            idummy = nitems(inputline)
            if(idummy .le. 0) then
              flag_go = .false.
            elseif(idummy .gt. 0) then
              flag_go = .true.
              backspace(my_pid)
            endif
          enddo
       
          if(flag_get_ldos) then ! make sure that your band_strucgture file contains orbital/wavefunction information
            if(flag_wf) then
              read(my_pid,*) kp(:), E(ie), C_(:)
              do im = 1, nbasis
                if(flag_phase) phase = F_IJ(-(kp(:)+gg(:)), PGEOM%o_coord_cart(:,im))
                if(ispinor .eq. 2) then
                  C(im,ie) = cmplx(C_(im*4-3), C_(im*4-2)) * phase
                  C(im+nbasis,ie) = cmplx(C_(im*4-1), C_(im*4  )) * phase
                  c_up     = cmplx(C_(im*4-3), C_(im*4-2))
                  c_dn     = cmplx(C_(im*4-1), C_(im*4  ))
                  V(im,ie) = real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! orbital projection
                elseif(ispinor .eq. 1) then
                  C(im,ie) = cmplx(C_(im*2-1), C_(im*2)) * phase
                  c_up     = cmplx(C_(im*2-1), C_(im*2))
                  V(im,ie) = real( conjg(c_up)*c_up )
                endif
              enddo
            elseif(.not. flag_wf) then
              read(my_pid,*) kp(:), E(ie), V(1:nbasis,ie)
            endif
          else
            read(my_pid,*) kp(:), E(ie)
          endif

        enddo
        close(my_pid)
      endif

    elseif(.not. flag_formatted) then
      open(my_pid,file = trim(fname), form='unformatted',status='old', iostat=i_continue)

      ! read header
      read(my_pid) ikmode, flag_vector_, flag_single_, flag_erange_, flag_sparse_, nbasis_, nk_, nband_, &
                ispin_, nspin_, ispinor_, c_mode_
      if(flag_vector_ .and. flag_single_ .and. c_mode_ .eq. 'rh') then
        V4 = 0d0
      endif
      ! assume ikmode = 3 (by default if PRTSEPK = .TRUE.)

      if(flag_sparse_) E = -999d0

      if(flag_erange_) then 
        read(my_pid) flag_erange_, init_erange_, fina_erange_
      elseif(.not. flag_erange_) then
        read(my_pid) flag_erange_
      endif

      if(flag_sparse_) then
        read(my_pid) flag_sparse_, emin_, emax_, nemax_
      elseif(.not. flag_sparse_) then
        read(my_pid) flag_sparse_
      endif

      ! read main energy/wavefunction or orbital projection
      if(.not. flag_get_ldos) then ! save as asking flag_vector
        if(.not. flag_single_) then
          read(my_pid) ne_found, kp(:), (E(ie),ie=1,ne_found)
        elseif(flag_single_) then
          read(my_pid) ne_found, kp4(:), (E4(ie),ie=1,ne_found)
          kp = kp4
          E  = E4
        endif
      elseif(flag_get_ldos) then

        if(c_mode_ .eq. 'wf') then
          if(flag_phase) phase = F_IJ(-(kp(:)+gg(:)), PGEOM%o_coord_cart(:,im))
          if(.not. flag_single_) then
            read(my_pid) ne_found, kp(:), (E(ie),ie=1, ne_found)
          elseif(flag_single_) then
            read(my_pid) ne_found, kp4(:), (E4(ie),ie=1, ne_found)
            kp = kp4
            E  = E4
          endif
          
          do ie = 1, ne_found
            do im = 1, nbasis
              if(ispinor .eq. 2) then
                if(.not. flag_single_) then
                 !read(my_pid) C(im, ie), C(im + nbasis, ie)
                  read(my_pid) c_up, c_dn
                elseif(flag_single_) then
                  read(my_pid) c4_up, c4_dn
                  c_up = c4_up ; c_dn = c4_dn
                endif
                V(im,        ie) = real( conjg(c_up)*c_up + conjg(c_dn)*c_dn)
                C(im,        ie) = c_up * phase
                C(im+nbasis, ie) = c_dn * phase
              elseif(ispinor .eq. 1) then
                if(.not. flag_single_) then
                  read(my_pid) c_up
                elseif(flag_single_) then
                  read(my_pid) c4_up
                  c_up = c4_up
                endif
                V(im,        ie) = real( conjg(c_up)*c_up )
                C(im,        ie) = c_up * phase
              endif
            enddo
          enddo
        
        elseif(c_mode_ .eq. 'rh') then
          if(.not. flag_single_) then
            read(my_pid) ne_found, kp(:), (E(ie),ie=1, ne_found)
          elseif(flag_single_) then
            read(my_pid) ne_found, kp4(:), (E4(ie),ie=1, ne_found)
            kp = kp4
            E  = E4
          endif

          do ie = 1, ne_found
            do im = 1, nbasis
              if(.not. flag_single_) then
                read(my_pid) V(im, ie) 
              elseif(flag_single_) then
                read(my_pid) V4(im, ie) 
              endif
            enddo
          enddo
          V = V4
        endif

      endif ! get_ldos
      close(my_pid)
    endif  ! flag_unformatted

    return
endsubroutine
subroutine load_band_structure_bin(E, ne_found, ispin, nspin, nband, nbasis, nk, kmode, pid, flag_exit)
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
   logical                flag_go
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

   ! read header
   read(pid) ikmode, flag_vector_, flag_single_, flag_erange_, flag_sparse_, nbasis_, nk_, nband_, &
                     ispin_, nspin_, ispinor_, c_mode_
   
   ! allocate temporal kpoint_ array which would not be used but (should be) same as PKPTS%kpoint
   allocate(kpoint_(ikmode,nk))

   if(flag_single_) then
     allocate(kpoint4_(ikmode,nk))
   endif
  
   ! set E as -999d0 for all (ie,ik) if flag_sparse_ by default. 
   ! if ne_found(ik) is not zero, E(1:ne_found(ik),ik) will be filled in the following step
   ! otherwise, remain as -999d0
   ! this setting is necessary to avoid unexpected DOS which is originated from E(.not. 1:ne_found(ik),ik)
   if(flag_sparse_) then
     E = -999d0
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
     if(.not. flag_single_) then
       do ik = 1, nk
         read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
       enddo
     elseif(flag_single_) then
       do ik = 1, nk
         read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
       enddo
       kpoint_ = kpoint4_
       E  = E4
     endif
   endif ! flag_exit

   deallocate(kpoint_)
   if(allocated(kpoint4_)) deallocate(kpoint4_)

   return
endsubroutine
