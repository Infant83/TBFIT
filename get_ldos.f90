#include "alias.inc"
program get_ldos
   use parameters, only : incar, poscar, hopping,  pid_geom,  zi 
   use mpi_setup
   use time
   use print_io
   use do_math, only : fgauss
   use set_default, only : init_poscar
   implicit none
   real*8                    t_start, t_end  
   integer*4                 mpierr
   integer*4                 nediv, nkp, nemax
   integer*4                        nkp_, nkp_x2
   integer*4                 nkp_sldos, ne_sldos
   integer*4                 nkdiv
   integer*4                 ii, i, j, k, ie, my_k, im, ik
   character*40              ie_str, fe_str
   integer*4                 d, iatom
   integer*4                 igx, igy, igz, ig(3), ik_shift
   integer*4                 ispinor
   integer*4                 natom_ldos,  nspec, norb
   integer*4                 my_pid
   integer*4                 narg,iarg
   integer*4                 lmmax
   integer*4, allocatable :: ne_found(:)
   integer*4, allocatable :: ne_found_(:)
   integer*4, allocatable :: iatom_ldos(:)
   integer*4, allocatable :: kidx_sldos(:), eidx_sldos(:)
   integer*4                 nx, ny, nz, ix, iy, iz
   real*8                    erange_i, erange_f, sigma
   real*8                    dos_, ldos_, ldos2_
   real*8                    weight
   real*8                    dk
   real*8                    T(3,3), coord(3), r_origin(3)
   real*8,    allocatable :: coord_cart(:,:)
   real*8,    allocatable :: E(:,:), V(:,:)
   complex*16,allocatable ::         C(:,:)
   complex*16                c_up, c_dn 
   complex*16                sw_up, sw_dn, sw
   real*8                    sw_real
   real*8,    allocatable :: E_(:,:)
   real*8,    allocatable :: e_range(:)
   real*8,    allocatable :: ldos(:,:), dos(:)
   real*8,    allocatable :: ldos__(:,:), dos__(:)
   real*8,    allocatable :: ldos_k(  :,:), ldos2_k(:,:), dos_k(:,:)
   real*8,    allocatable :: ldos_k_(  :,:), ldos2_k_(:,:), dos_k_(:,:)
   real*8,    allocatable :: ldos3_k(:,:,:), ldos3_k_(:,:,:)
   real*8,    allocatable :: ldos4_k(:,:,:), ldos4_k_(:,:,:)
   real*8,    allocatable :: ldos5_k(:,:,:), ldos5_k_(:,:,:)
   real*8,    allocatable :: ldos4_sum(:), ldos4_sum_(:)
   real*8,    allocatable :: ldos_kk( :,:)
   real*8,    allocatable :: sw_k(:,:), sw2_k(:,:)
   real*8,    allocatable :: sw_k_(:,:), sw2_k_(:,:)
   real*8,    allocatable :: sw3_k(:,:), sw3_k_(:,:)
   real*8,    allocatable :: kk(:,:), kline(:), kline_(:)
   real*8,    allocatable :: kk_(:,:)
   real*8,    allocatable :: kp(:,:), kp_(:,:), kline_p(:), kline_p_(:)
   real*8,    allocatable :: sldos_sum(:,:), sldos_sum_(:,:)
   real*8                    gg(3)
   real*8,    allocatable :: qpi(:,:), q(:)
   logical                   flag_read_orbital
   logical                   flag_get_band, flag_get_dos
   logical                   flag_get_ldos, flag_get_qpi, flag_ewindow, flag_phase, flag_get_unfold
   logical                   flag_wf, flag_spectral_band, flag_formatted, flag_get_sldos
   logical                   flag_kselect, flag_enselect, flag_sldos_select
   character*80              fname_in, fname, header, fname_geom
   character*40              ftype
   character*20,external  :: int2str
   integer*4, allocatable :: ourjob(:), ourjob_disp(:)
   type (incar )            :: PINPT
   type (poscar)            :: PGEOM
   type (hopping)           :: NN_TABLE
   character*40              fnamelog
   character*136              path, filenm_sldos
   character*136              path_sldos
   fnamelog = 'TBFIT.GET_LDOS.log'
   path     = './'        ! default path where the bands are stored
   path_sldos= './didv/'  ! default path where the sldos to be save
#ifdef MPI 
   call mpi_initialize(fnamelog)
#else
   call open_log(fnamelog)
#endif
   flag_get_ldos   = .false.
   flag_get_qpi    = .false.
   flag_get_unfold = .false.
   flag_phase      = .false.
   flag_wf         = .false.
   flag_formatted  = .true.  ! true: dat, false: bin
   flag_kselect    = .false.
   flag_enselect   = .false.
   flag_sldos_select = .false.
   flag_read_orbital = .false.
   flag_get_band   = .true.
   flag_get_dos    = .false.

   ig = 1 ! default (should not be less than 1 and should be 1 if flag_get_unfold=False
   nediv = 1
   allocate(ourjob(nprocs))
   allocate(ourjob_disp(0:nprocs-1))
!  parse basic information
   ispinor = 2 ! SOC case
   fname_in = 'ldos.in'
  !fname_geom = 'POSCAR-TB'
   narg = iargc()
   nkp_sldos = 0

! IMPORTANT NOTE: The spin collinear case has not been implemented. Please be aware and need ot be implemented 
!                 in the near future. HJ Kim. 17. March. 2021


   if(narg .eq. 1) call getarg(1, fname_in)

   open(999, file=trim(fname_in), status='old')
     write(message,'(A)')"-- title of the calculation (without blank)"  ; write_msg
       read(999,*)header
       write(message,'(A,A   )')"   CASE = ",trim(header)  ; write_msg
     write(message,'(A)')"-- type of data storage: binary-> bin, ascii-> dat"  ; write_msg
       read(999,*)ftype
       if(ftype(1:1) .eq. 'b' .or. ftype(1:1) .eq. 'B') then
         write(message,'(A,A   )')"   FORMAT = ",trim(ftype)  ; write_msg
         flag_formatted = .FALSE.
       elseif(ftype(1:1) .eq. 'd' .or. ftype(1:1) .eq. 'D') then
         write(message,'(A,A   )')" FORMAT = ",trim(ftype)  ; write_msg
         flag_formatted = .TRUE. 
       endif
     write(message,'(A)')"-- path of data storage (excluding data file name): ex) ./PATH" ; write_msg
       read(999,*) path
       write(message,'(A,A)')"   PATH = ",trim(path) ; write_msg
     write(message,'(A)')"-- name of the geometry file"  ; write_msg
       read(999,*)fname_geom
       write(message,'(A,A   )')"   GEOMETRY = ",trim(fname_geom)  ; write_msg
       iverbose = 2
       call init_poscar(PGEOM)
       PGEOM%gfilenm = trim(fname_geom)
       call read_poscar(PINPT, PGEOM, NN_TABLE)
       iverbose = 1
       ! NONTE: probably this will works fine with flag_slater_koster
       !        but did not check with non-slater-koseter type orbitals. Need some check
       !        16.March.2021 HJ Kim
       PINPT%lmmax = 9
       call set_orbital_index(PGEOM, PINPT%lmmax)
       lmmax = PINPT%lmmax
       allocate(coord_cart(3,PGEOM%n_atom))
       do iatom=1, PGEOM%n_atom
          coord(:) = PGEOM%a_coord(:,iatom) + r_origin(:) ; coord(:) = coord(:) - nint(coord(:))
          coord_cart(:,iatom) = coord(1)*PGEOM%a_latt(:,1)+coord(2)*PGEOM%a_latt(:,2)+coord(3)*PGEOM%a_latt(:,3)
       enddo
     write(message,'(A)')"-- SOC (yes:2, no:1)?"   ; write_msg
       read(999,*)ispinor
       write(message,'(A,(I6))')"   ISPINOR=",ispinor  ; write_msg
     write(message,'(A)')"-- NEIG? (how many eigenvalues?, write ne_max written in band_structure_TBA.kp_xx.dat only if EWINDOW mode)"  ; write_msg
       read(999,*)nemax
       write(message,'(A,(I6))')"   NEMAX=",nemax  ; write_msg
     write(message,'(A)')"-- NSPEC ? (ex, if the system is consist of carbon and oxygen: 2)"  ; write_msg
       write(message,'(A,(I6))')"   NSPEC=",PGEOM%n_spec  ; write_msg
     write(message,'(A)')"-- NATOM for each species ? (ex, 2 carbon 2 oxygen: 2 2)"  ; write_msg
       write(message,'(A,*(I6))')"   NATOM_SPEC=" ,PGEOM%i_spec  ; write_msg
     write(message,'(A)')"-- NKP ?"  ; write_msg
       read(999,*)nkp
       write(message,'(A,*(I6))')"   NKP = " ,nkp  ; write_msg
       allocate(kk(3,nkp)) ; kk = 0d0
       allocate(kk_(3,nkp)); kk_= 0d0
       allocate(kline(nkp)); kline = 0d0
       allocate(E(nemax,nkp)) ; E = .0d0
       allocate(E_(nemax,nkp)) ; E_ = .0d0
       allocate(ne_found(nkp)); ne_found = 0
       allocate(ne_found_(nkp)); ne_found_ = 0
     write(message,'(A)')"-- NKDIV ?"  ; write_msg
       read(999,*)nkdiv
       write(message,'(A,*(I5))')"   NKDIV = ",nkdiv  ; write_msg
     write(message,'(A)')"-- EWINDOW? (SPARSE ?) "  ; write_msg
       read(999,*)flag_ewindow
       write(message,'(A,(L1))')"   FLAG EWINDOW= ", flag_ewindow  ; write_msg
     write(message,'(A)')"-- Plot band structure?"  ; write_msg
       read(999,*)flag_get_band
       write(message,'(A,(L1))')"   FLAG GET BAND= ", flag_get_band ; write_msg
     write(message,'(A)')"-- Plot DOS ? "  ; write_msg
       read(999,*)flag_get_dos 
       write(message,'(A,(L1))')"   FLAG GET DOS= ", flag_get_dos  ; write_msg
     write(message,'(A)')"-- Gaussian smearing?"  ; write_msg
       read(999,*)sigma
       write(message,'(A,(F10.4),A)')"   SMEARING = ",sigma, ' (eV)'  ; write_msg
     write(message,'(A)')"-- Number of division? (NEDIV)"  ; write_msg
       read(999,*)nediv
      !if(nediv < 1) nediv = 1
       write(message,'(A,(I6))')"   NEDIV = ",nediv  ; write_msg
       allocate( dos_k(nediv,nkp)) ; dos_k = 0d0
       allocate( dos_k_(nediv,nkp)) ; dos_k_ = 0d0
       allocate( dos(nediv)) ; dos = 0d0
       allocate( dos__(nediv)) ; dos__ = 0d0
     write(message,'(A)')"-- Energy range for dos plot (erange_i  erange_f)"  ; write_msg
       read(999,*)erange_i, erange_f
       write(message,'(A,F10.4,A,F10.4)')"   ERANGE = ",erange_i,' : ',erange_f  ; write_msg

     write(message,'(A)')"-- Read wavefunction (wf, .true.) or density (rh, .false.)? Logical value (.true. or .false.)"  ; write_msg
       read(999,*)flag_wf
       if(flag_wf) then
         write(message,'(A)')'   WAVEFUNCTION |psi_nk> READ = .TRUE.'  ; write_msg
       elseif(.not. flag_wf) then
         write(message,'(A)')'   RHO READ = .TRUE.'  ; write_msg
       endif

     write(message,'(A)')"-- GET LDOS ? (evaluate LDOS from band_structure_TBA.dat with specified energy range)"  ; write_msg
       read(999,*)flag_get_ldos
       write(message,'(A,L1,I0)')"   GET_LDOS = ", flag_get_ldos ; write_msg
       if(flag_get_ldos) then
         allocate(ldos_k(      nediv,nkp)) ; ldos_k = 0d0
         allocate(ldos_k_(      nediv,nkp)) ; ldos_k_ = 0d0
         allocate(ldos2_k(nemax,nkp)) ; ldos2_k = 0d0
         allocate(ldos2_k_(nemax,nkp)) ; ldos2_k_ = 0d0
         allocate(ldos3_k(nemax,lmmax,nkp)) ; ldos3_k = 0d0
         allocate(ldos3_k_(nemax,lmmax,nkp)) ; ldos3_k_ = 0d0
         allocate(ldos(PGEOM%n_atom,nediv)) ; ldos = 0d0
         allocate(ldos__(PGEOM%n_atom,nediv)) ; ldos__ = 0d0
       endif
       if(flag_get_ldos) then ; write(message,'(A)')"-- Number of atoms for LDOS plot "  ; write_msg ; endif
         read(999,*)natom_ldos
       if(flag_get_ldos) then ; write(message,'(A,(I6))')"   NATOM_LDOS = ",natom_ldos  ; write_msg ; endif
         allocate(iatom_ldos(natom_ldos))
       if(flag_get_ldos) then ; write(message,'(A)')"-- Atom index for ldos plot"  ; write_msg ;endif
         read(999,*)iatom_ldos(1:natom_ldos)
       if(flag_get_ldos) then ; write(message,'(A,*(I5))')"   LDOS ATOM INDEX = ",iatom_ldos  ; write_msg ;endif

     write(message,'(A)')"-- GET UNFOLD? "  ; write_msg
       read(999,*)flag_get_unfold
       write(message,'(A,L)')"   GET UNFOLD = ", flag_get_unfold  ; write_msg
       if(flag_get_unfold) then; write(message,'(A)')"-- Multiply phase factor to each wavefunction coefficient?"  ; write_msg ;endif
         read(999,*)flag_phase
       if(flag_get_unfold) then ; write(message,'(A,L)')"   SET PHASE = ", flag_phase  ; write_msg ; endif
       if(flag_get_unfold) then ; write(message,'(A)')"-- G-vector direction: ngx, ngy, ngz "  ; write_msg ; endif
         read(999,*)ig(1:3)
       if(flag_get_unfold)  then ; write(message,'(A,3(I6))')"   G_DIRECT= ",ig(:)  ; write_msg ; endif
       if(flag_get_unfold) then
         allocate(sw_k  (       nediv,nkp*ig(1)*ig(2)*ig(3))) ;    sw_k = 0d0
         allocate(sw_k_ (       nediv,nkp*ig(1)*ig(2)*ig(3))) ;    sw_k_= 0d0
         allocate(sw2_k (       nemax,nkp*ig(1)*ig(2)*ig(3))) ;   sw2_k = 0d0
         allocate(sw2_k_(       nemax,nkp*ig(1)*ig(2)*ig(3))) ;   sw2_k_= 0d0
       endif

     write(message,'(A)')"-- GET QPI ? "  ; write_msg
       read(999,*)flag_get_qpi 
       write(message,'(A,(L1))')"   GET_QPI = ", flag_get_qpi  ; write_msg
       if(flag_get_qpi) then
         allocate(ldos_kk(     nediv,nkp*ig(1)*ig(2)*ig(3))) ; ldos_kk= 0d0
         allocate(kline_(nkp*ig(1)*ig(2)*ig(3))); kline_ = 0d0
         allocate(   qpi(      nediv,nkp*2)) ;    qpi = 0d0
         allocate(     q(            nkp*2)) ;      q = 0d0
       endif
       if(flag_get_qpi) then ;write(message,'(A)')"-- Weight "  ; write_msg ;endif
         read(999,*)weight
       if(flag_get_qpi) then;  write(message,'(A,F10.4        )')"   WEIGHT = ",weight  ; write_msg ; endif

     write(message,'(A)')"-- Spectral band plot ? "  ; write_msg
       read(999,*)flag_spectral_band
       write(message,'(A,(L1))')"   SPECTRAL_BAND = ", flag_spectral_band  ; write_msg
       if(flag_spectral_band) then
         allocate(sw3_k (       nemax,nkp                  )) ;   sw3_k = 0d0
         allocate(sw3_k_(       nemax,nkp                  )) ;   sw3_k_= 0d0
       endif

     write(message,'(A,F9.4,A,F9.4)')"-- Spatial LDOS within erange= ",erange_i,' : ',erange_f; write_msg
       read(999,*)flag_get_sldos
       write(message,'(A,(L1))')"    SLDOS = ", flag_get_sldos
       if(flag_get_sldos) then
         allocate(sldos_sum(PGEOM%n_atom, nediv)) ; sldos_sum = 0d0
         allocate(sldos_sum_(PGEOM%n_atom, nediv)) ; sldos_sum_ = 0d0
       endif
       if(flag_get_sldos)  then ; write(message,'(A)')'-- cell repeat?: nx, ny, nz' ; write_msg ; endif
         T = PGEOM%a_latt
         read(999,*)nx, ny, nz
         if(flag_get_sldos)  then ;write(message,'(A, I3,I3,I3)')'   CELL_REPEAT (nx, ny, nz) = ', nx, ny, nz ; write_msg ; endif
       if(flag_get_sldos) then ; write(message,'(A)')'-- Shift of origin of atomic coord.(fractional coord)' ; write_msg ; endif
         read(999,*)r_origin(1), r_origin(2), r_origin(3)
         if(flag_get_sldos)  then;  write(message,'(A,3(F9.4))')'   R_ORIGIN: ',r_origin(:); write_msg ; endif
       if(flag_get_sldos)  then;  write(message,'(A)')"-- path to SLDOS be stored: ex) './didv/' "; write_msg; endif
         read(999,*)path_sldos
         if(flag_get_sldos) then ; write(message,'(A,A)')"   SLDOS SAVE PATH = ",trim(path_sldos) ; write_msg ; endif

     write(message,'(A)')'-- Spatial LDOS with k-point selective ?'
       read(999,*)flag_kselect
       write(message,'(A,L)')"    SELECT K = ",flag_kselect
       if(flag_kselect) then ; write(message,'(A)')'-- number of k-points to be resolved' ; write_msg; endif
       read(999,*)nkp_sldos
       if(flag_kselect) then ; write(message,'(A,I0)')'   NKP RESOLVE = ', nkp_sldos ; write_msg; endif
       if(flag_kselect) then ; write(message,'(A)')'-- indices of k-points to be resolved' ; write_msg; endif
       allocate(kidx_sldos(nkp_sldos)) ; kidx_sldos = 0
       read(999,*)kidx_sldos(1:nkp_sldos)
       if(flag_kselect) then ; write(message,'(A,*(I5))')'   IKP RESOLVE = ', kidx_sldos ; write_msg; endif
     write(message,'(A)')'-- Spatial LDOS with energy level selective ?'
       read(999,*)flag_enselect
       write(message,'(A,L)')"    SELECT EN = ",flag_enselect
       if(flag_enselect) then ; write(message,'(A)')'-- number of eigenvalues to be resolved per each k-point' ; write_msg; endif
       read(999,*)ne_sldos
       if(flag_enselect) then ; write(message,'(A,I0)')'   N-EN RESOLVE = ', ne_sldos ; write_msg; endif
       if(flag_enselect) then ; write(message,'(A)')'-- indices of energy level to be resolved' ; write_msg; endif
       allocate(eidx_sldos(ne_sldos)) ; eidx_sldos = 0
       read(999,*)eidx_sldos(1:ne_sldos)
       if(flag_enselect) then ; write(message,'(A,*(I6))')'   IEN RESOLVE = ', eidx_sldos ; write_msg; endif
     if(flag_kselect .and. flag_enselect) flag_sldos_select = .True.
     if(flag_sldos_select) then
       allocate(ldos4_k(PGEOM%n_atom,ne_sldos,nkp_sldos))  ; ldos4_k = 0d0
       allocate(ldos4_k_(PGEOM%n_atom,ne_sldos,nkp_sldos)) ; ldos4_k_ = 0d0
       allocate(ldos4_sum(PGEOM%n_atom))  ; ldos4_sum = 0d0
       allocate(ldos4_sum_(PGEOM%n_atom)) ; ldos4_sum_ = 0d0
     endif
     if(flag_kselect .and. .not. flag_enselect) then
       allocate(ldos5_k( PGEOM%n_atom, nediv, nkp_sldos)) ; ldos5_k = 0d0
       allocate(ldos5_k_(PGEOM%n_atom, nediv, nkp_sldos)) ; ldos5_k_= 0d0
     endif

     write(message,'(A)')" END READING ldos.in"  ; write_msg
     write(message,'(A)')" " ; write_msg
     close(999)

!  setup system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(nediv .gt. 1) then
     e_range =  erange_i + (/(k, k=0,nediv-1)/) * (erange_f - erange_i)/dble(nediv - 1)
     e_range(0) = e_range(1) - (e_range(2) - e_range(1))
   elseif(nediv .eq. 1) then
     allocate(e_range(0:1))
     e_range =  erange_i
     e_range(0) = e_range(1)
   endif

   norb = sum(PGEOM%n_orbital)
   allocate(V(norb,nemax)) ; V = 0d0
   allocate(C(norb*ispinor,nemax)) ; C = 0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if( flag_kselect .or. flag_sldos_select ) then
     do k = 1, nkp_sldos
       if(flag_formatted) then ! need to extend toward spin polarized case in the future. 18. March. 2021, HJ Kim
         fname = trim(path)//'./band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(kidx_sldos(k))))//'.dat'
       else
         fname = trim(path)//'./band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(kidx_sldos(k))))//'.bin'
       endif
       E(:,kidx_sldos(k)) = -999d0 ; gg = 0d0
       call load_band_singlek(PGEOM, fname, nemax, norb, ispinor, E(:,kidx_sldos(k)), C, V, ne_found(kidx_sldos(k)), kk(:,kidx_sldos(k)), gg, &
                              .TRUE., flag_wf, flag_ewindow, flag_phase, flag_formatted)
       if(flag_sldos_select) then
         write(message,'(A    )')" "; write_msg
         write(message,'(A    )')" SLDOS (K & EN)-select mode"; write_msg
         call mpi_job_distribution_chain(ne_sldos, nprocs, ourjob, ourjob_disp)
         do i = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
           write(message,'(A,I0,A,I0    )')" STATUS KP: ", kidx_sldos(k),' EN: ', eidx_sldos(i) ; write_msg_all
           do j = 1,  PGEOM%n_atom
             ii = sum(PGEOM%n_orbital(1:j)) - PGEOM%n_orbital(j) + 1 ! starting orbital index for atom j
             ldos4_k(j, i, k) = ldos4_k(j,i,k) + sum(V(ii:ii+PGEOM%n_orbital(j)-1,eidx_sldos(i)))
             ldos4_sum(j) = ldos4_sum(j) + sum(V(ii:ii+PGEOM%n_orbital(j)-1,eidx_sldos(i)))
           enddo
         enddo
       elseif(flag_kselect .and. .not. flag_enselect) then
         write(message,'(A    )')" "; write_msg
         write(message,'(A    )')" SLDOS K-select mode"; write_msg
         write(message,'(A,I0 )')" STATUS KP: ", kidx_sldos(k); write_msg
         call mpi_job_distribution_chain(nediv, nprocs, ourjob, ourjob_disp)
         do ie = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
           write(message,'(A,F10.4,A)')"    STAT EN: ", real(ie - sum(ourjob(1:myid)))/real(sum(ourjob(1:myid+1)))*100.d0 , ' %' ; write_msg
           do i = 1, ne_found(kidx_sldos(k))
             dos_ = fgauss(sigma, e_range(ie) - E(i,kidx_sldos(k))) / real(nkp)
             do j = 1, PGEOM%n_atom
               ii = sum(PGEOM%n_orbital(1:j)) - PGEOM%n_orbital(j) + 1 ! starting orbital index for atom j
               ldos5_k(j,ie, k) = ldos5_k(j,ie, k) + dos_ * sum(V(ii:ii+PGEOM%n_orbital(j)-1,i))
             enddo
           enddo
         enddo
       elseif(.not. flag_kselect .and. flag_enselect) then
         write(message,'(A    )')" "; write_msg
         write(message,'(A    )')" !! WARN !! FLAG_KSELECT=False and FLAG_ENSELECT=True"; write_msg
         write(message,'(A    )')"            Such combination is not allowed. " ; write_msg
         write(message,'(A    )')"            Instead, you can set erange with same values (0.3 0.3 for example) " ; write_msg
         write(message,'(A    )')"            and set nediv 1 so that certain energy level will be calculated" ; write_msg
         write(message,'(A    )')"            with SLDOS plot mode."  ; write_msg
         kill_job
       endif

     enddo !k
   endif

#ifdef MPI
call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif 

   if(flag_get_band .or. flag_get_dos .or. flag_get_unfold .or. flag_get_sldos .or. flag_spectral_band .or. flag_get_ldos) then
   write(message,'(A)')" "; write_msg
   do igz=0,ig(3)-1
   do igy=0,ig(2)-1
   do igx=0,ig(1)-1
     ik_shift=(igx+igy+igz)*nkp
     call mpi_job_distribution_chain(nkp, nprocs, ourjob, ourjob_disp)
     gg = PGEOM%b_latt(:,1) * real(igx) + PGEOM%b_latt(:,2) * real(igy) + PGEOM%b_latt(:,3) * real(igz) ! G vector 
     do k = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
       my_k = k - sum(ourjob(1:myid))
       if(flag_formatted) then ! need to extend toward spin polarized case in the future. 18. March. 2021, HJ Kim
         fname = trim(path)//'./band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(k)))//'.dat'
       else
         fname = trim(path)//'./band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(k)))//'.bin'
       endif
       E(:,k) = -999d0

       flag_read_orbital = (flag_get_ldos .or. flag_get_unfold .or. flag_get_sldos .or. flag_spectral_band)
       call load_band_singlek(PGEOM, fname, nemax, norb, ispinor, E(:,k), C, V, ne_found(k), kk(:,k), gg, &
                              flag_read_orbital, flag_wf, flag_ewindow, flag_phase, flag_formatted)

       if(flag_get_unfold) then
         write(message,'(A,F10.4,A,3I3)')" STAT KP: ", real(my_k)/real(sum(ourjob(1:myid+1)))*100,' %, ig(1),ig(2),ig(3) = ',igx+1, igy+1, igz+1 ; write_msg
       else
         write(message,'(A,F10.4,A    )')" STAT KP: ", real(my_k)/real(sum(ourjob(1:myid+1)))*100,' %' ; write_msg
       endif
       do ie = 1, nediv
         if_main write(6,'(A,F10.4,A)')"    STAT EN: ", real(ie)/real(nediv)*100.d0 , ' %'
         do i = 1, ne_found(k)
           dos_ = fgauss(sigma, e_range(ie) - E(i,k)) / real(nkp)
           if( (igx+igy+igz) .eq. 0 ) then
             dos(ie) = dos(ie) + dos_
             dos_k(ie,k) = dos_k(ie,k) + dos_
           endif
       
           if( (flag_get_ldos .or. flag_get_sldos)  .and. (igx+igy+igz) .eq. 0) then
             if(flag_get_ldos ) then
               do j = 1, natom_ldos
                 ii = sum(PGEOM%n_orbital(1:iatom_ldos(j))) - PGEOM%n_orbital(iatom_ldos(j)) + 1 ! starting orbital index for atom iatom_ldos(j)
                 ldos_ =  dos_ * sum(V(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1,i))
                 ldos(j,ie) = ldos(j,ie) + ldos_
                 ldos_k(ie,k) = ldos_k(ie,k) + ldos_
       
                 if(ie .eq. 1) then ! for atom projected band
                   ldos2_ = sum(V(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1,i))
                   ldos2_k(i,k) = ldos2_k(i,k) + ldos2_
                   ldos3_k(i,1:PGEOM%n_orbital(iatom_ldos(j)),k) = ldos3_k(i,1:PGEOM%n_orbital(iatom_ldos(j)),k) + V(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1,i)
                 endif
               enddo
             endif
       
             if(flag_get_sldos) then
               do j = 1, PGEOM%n_atom
                 ii = sum(PGEOM%n_orbital(1:j)) - PGEOM%n_orbital(j) + 1 ! starting orbital index for atom j
                 sldos_sum(j, ie) = sldos_sum(j, ie) + dos_ * sum(V(ii:ii+PGEOM%n_orbital(j)-1,i))
               enddo
             endif
           endif
       
           if(flag_get_unfold) then
             sw_up = 0d0; sw_dn = 0d0; sw_real = 0d0
             do j = 1, natom_ldos
               ii = sum(PGEOM%n_orbital(1:iatom_ldos(j))) - PGEOM%n_orbital(iatom_ldos(j)) + 1 ! starting orbital index for atom iatom_ldos(j)
               ! NOTE: to get unfolding band flag_wf should be TRUE
               if(ispinor .eq. 2) then
                 c_up = sum(C(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1, i)) 
                 c_dn = sum(C(norb+ii:norb+ii+PGEOM%n_orbital(iatom_ldos(j))-1, i))
                 sw_up= sw_up + c_up ; sw_dn = sw_dn + c_dn
               else
                 c_up = sum(C(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1, i)) 
                 sw_up= sw_up + c_up
               endif
             enddo
             if(ispinor .eq. 2) then
               sw_k(ie,k+ik_shift) = sw_k(ie,k+ik_shift) + real(conjg(sw_up)*sw_up + conjg(sw_dn)*sw_dn) * dos_ 
             else
               sw_k(ie,k+ik_shift) = sw_k(ie,k+ik_shift) + real(conjg(sw_up)*sw_up                     ) * dos_ 
             endif
             if(ie .eq. 1 .and. ispinor .eq. 2) sw2_k(i,k+ik_shift) = real(conjg(sw_up)*sw_up + conjg(sw_dn)*sw_dn) 
             if(ie .eq. 1 .and. ispinor .eq. 1) sw2_k(i,k+ik_shift) = real(conjg(sw_up)*sw_up                     ) 
           endif
       
           if(flag_spectral_band .and. .not. flag_get_unfold) then
             sw_up = 0d0; sw_dn=0d0 ; sw_real = 0d0
             do j = 1, natom_ldos
               ii = sum(PGEOM%n_orbital(1:iatom_ldos(j))) - PGEOM%n_orbital(iatom_ldos(j)) + 1
               if(flag_wf) then
                 if(ispinor .eq. 2) then
                   c_up = sum(C(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1, i))
                   c_dn = sum(C(norb+ii:norb+ii+PGEOM%n_orbital(iatom_ldos(j))-1, i))
                   sw_up= sw_up + c_up ; sw_dn = sw_dn + c_dn
                 else
                   c_up = sum(C(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1, i))
                   sw_up= sw_up + c_up
                 endif
               elseif(.not. flag_wf) then
                 sw_real= sw_real + sum(V(ii:ii+PGEOM%n_orbital(iatom_ldos(j))-1, i))
               endif
             enddo
       
             if(flag_wf) then
               if(ispinor .eq. 2) then
                 sw3_k(ie, k) = sw3_k(ie, k) + real(conjg(sw_up)*sw_up + conjg(sw_dn)*sw_dn) * dos_
               else
                 sw3_k(ie, k) = sw3_k(ie, k) + real(conjg(sw_up)*sw_up                     ) * dos_
               endif
             elseif(.not. flag_wf) then
               sw3_k(ie, k) = sw3_k(ie, k) + sw_real * dos_
             endif
       
           elseif(flag_spectral_band .and. flag_get_unfold) then
             write(message,'(A)')' !WARN! You have requested band unfolding and spectral band plot simulatneously.'
             write(message,'(A)')' !WARN! Current version does not support it. Please operate seperately. Exit program...'
             kill_job
           endif
       
         enddo ! ie
       enddo ! nediv

     enddo ! ik
   enddo !ig
   enddo
   enddo
   endif ! flags..

#ifdef MPI
   call MPI_BARRIER(mpi_comm_earth, mpierr)

   call MPI_ALLREDUCE(E, E_,size(E),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   E = E_
   call MPI_ALLREDUCE(ne_found, ne_found_, size(ne_found), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
   ne_found = ne_found_
   call MPI_ALLREDUCE(kk, kk_, size(kk), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   kk = kk_
   if(flag_get_dos) then
     call MPI_ALLREDUCE(dos, dos__,size(dos), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     dos = dos__
     call MPI_ALLREDUCE(dos_k,dos_k_,size(dos_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     dos_k = dos_k_
   endif
   if(flag_get_ldos) then
     call MPI_ALLREDUCE(ldos, ldos__,size(ldos), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos = ldos__
     call MPI_ALLREDUCE(ldos_k, ldos_k_, size(ldos_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos_k = ldos_k_
     call MPI_ALLREDUCE(ldos2_k, ldos2_k_, size(ldos2_k),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos2_k = ldos2_k_
     call MPI_ALLREDUCE(ldos3_k, ldos3_k_, size(ldos3_k),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos3_k = ldos3_k_
   endif
   if(flag_get_unfold) then
     call MPI_ALLREDUCE(  sw_k,   sw_k_, size(  sw_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
       sw_k =   sw_k_
     call MPI_ALLREDUCE(  sw2_k,   sw2_k_, size(  sw2_k),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
       sw2_k =   sw2_k_
   endif
   if(flag_spectral_band) then
     call MPI_ALLREDUCE(  sw3_k,   sw3_k_, size(  sw3_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     sw3_k =   sw3_k_
   endif
   if(flag_get_sldos) then
     call MPI_ALLREDUCE(sldos_sum, sldos_sum_, size(sldos_sum), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   endif
   if(flag_sldos_select) then
     call MPI_ALLREDUCE(ldos4_sum, ldos4_sum_, size(ldos4_sum), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos4_sum = ldos4_sum_
     call MPI_ALLREDUCE(ldos4_k,   ldos4_k_,   size(ldos4_k),   MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos4_k = ldos4_k_
   endif
   if(.not. flag_enselect .and. flag_kselect) then
     call MPI_ALLREDUCE(ldos5_k, ldos5_k_, size(ldos5_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos5_k = ldos5_k_
   endif
#endif
   write(message,'(A)')" Replot DOS ... Done!"  ; write_msg
   write(message,'(A)')" " ; write_msg

   call get_kline_dist(kk, nkp, kline)

   !########   write band structure #########################
   if(myid .eq. 0 .and. flag_get_band) then
     open(999,file='band_structure_TBA.'//trim(header)//'.dat',status='unknown')
     if(flag_get_ldos) write(999,'(A,*(1x,I0))')'# LDOS for atoms: ',iatom_ldos(:)
     do i = 1, nemax
       write(999,'(A,I0,A)')'# ',i, ' -th eigen'
       write(999,'(A     )',ADVANCE='NO')'#       k-dist          E(eV)   '
       do im = 1, PGEOM%n_orbital(iatom_ldos(1))
         write(999,'(2X, A8)', ADVANCE='NO')     trim(PGEOM%c_orbital(im,iatom_ldos(1)))
       enddo
       write(999,'(A)',ADVANCE='YES')'  tot(atom_sum)'
       do k = 1,nkp
         if(flag_get_ldos) then
           write(999,'(2F16.6, *(F10.4))')kline(k), E(i,k), ldos3_k(i,1:PGEOM%n_orbital(iatom_ldos(1)),k), ldos2_k(i,k)
         else
           write(999,'(2F16.6          )')kline(k), E(i,k)
         endif
       enddo
       write(999,*)' '
       write(999,*)' '
     enddo
     do i = 1, nemax
       write(999,'(A,I0,A)')'# ',i, ' -th eigen'
       write(999,'(A     )',ADVANCE='NO')'#       k-dist          E(eV)   '
       do im = 1, PGEOM%n_orbital(iatom_ldos(1))
         write(999,'(2X, A8)', ADVANCE='NO')     trim(PGEOM%c_orbital(im,iatom_ldos(1)))
       enddo
       write(999,'(A)',ADVANCE='YES')'  tot(atom_sum)'
       do k = 1,nkp
         if(flag_get_ldos) then
           write(999,'(2F16.6, *(F10.4))')kline(k), E(i,k), ldos3_k(i,1:PGEOM%n_orbital(iatom_ldos(1)),k), ldos2_k(i,k)
         else
           write(999,'(2F16.6          )')kline(k), E(i,k)
         endif
       enddo
       write(999,*)' '
       write(999,*)' '
     enddo

     write(message,'(A)')' -> Band structure (+atom projected) written: band_structure_TBA.'//trim(header)//'.dat' ; write_msg
     close(999)
   endif
   !######################################################

   !###### write spectral band ###########################
   if(myid .eq. 0 .and. flag_spectral_band) then
     open(999,file='band_structure_TBA.SPECTRAL.'//trim(header)//'.dat',status='unknown')
     write(999,'(A,*(1x,I0))')'# Spectral band arround local atoms: ',iatom_ldos(:)
     do k = 1, nkp
       do ie = 1, nediv
         write(999, '(F16.8, F16.8, F16.8)')kline(k), e_range(ie), sw3_k(ie,k)
       enddo
       write(999,*)' '
     enddo
     write(message,'(A)')' -> Spectral band is written: band_structure_TBA.SPECTRAL.'//trim(header)//'.dat ' ; write_msg        
   endif
   !######################################################
    
   !#  write unfold band  ################################
   if(myid .eq. 0 .and. flag_get_unfold) then  
     ! IMPORTANT NOTE:  the result of this routine is only valid if you just draw 1D k-path so that applying gg to kk is just valid 
     open(999,file='band_structure_TBA.unfold.dat',status='unknown')
     write(999,'(A,*(1x,I0))')'# LDOS for atoms: ',iatom_ldos(:)
     do i = 1, nemax
       write(999,'(A,I0,A)')'# ',i, ' -th eigen'
       do igz=0,ig(3) - 1
       do igy=0,ig(2) - 1
       do igx=0,ig(1) - 1
         do k = 1,nkp
           write(999,'(2F16.6, *(F10.4))')kline(k)+kline(nkp)*(igx+igy+igz), E(i,k),   sw2_k(i,k+nkp*(igx+igy+igz))
         enddo
       enddo
       enddo
       enddo
       write(999,*)' '
       write(999,*)' '
     enddo
     write(message,'(A)')' -> Band structure (+ spectral weight, + atom projected) written: band_structure_TBA.unfold.dat' ; write_msg

     close(999)
   endif
   !######################################################

   !# write DOS ##########################################
   if(myid .eq. 0 .and. flag_get_dos) then
     open(999,file='DOS.total.dat',status='unknown')
     if(flag_get_ldos) then
       open(998,file='LDOS.'//trim(header)//'.dat',status='unknown')
       write(998,'(A,*(1x,I0))')'# LDOS for atoms: ',iatom_ldos(:)
     endif
     do i =1, nediv
       write(999, '(2F16.4)') e_range(i), dos(i)
       if(flag_get_ldos) then
         write(998, '(3F16.4)') e_range(i), dos(i), sum(ldos(:,i))
       endif
     enddo
     write(message,'(A)')' -> DOS written: DOS.total.dat and LDOS.'//trim(header)//'.dat' ; write_msg
     close(999) 
     if(flag_get_ldos) close(998)
   endif
   !######################################################

   !####### calculate spectral weight (band unfolding) ###
   if(myid .eq. 0 .and. flag_get_unfold) then 
     open(999,file='SW.'//trim(header)//'.dat',status='unknown')
     write(999,'(A,*(1x,I0))')'# Spectral weight arround local atoms: ',iatom_ldos(:)

     do igz=0,ig(3)-1
     do igy=0,ig(2)-1
     do igx=0,ig(1)-1
       do k = 1, nkp
         do ie = 1, nediv
           write(999, '(F16.8, F16.8, F16.8)')kline(k)+kline(nkp)*(igx+igy+igz), e_range(ie), sw_k(ie,k+nkp*(igx+igy+igz))
         enddo
         write(999,*)' '
       enddo
     enddo
     enddo
     enddo

     close(999)
     write(message,'(A)')' -> Spectral weight into extended (unfolded) BZ is written: SW.'//trim(header)//'.dat ' ; write_msg
   endif
   !######################################################

   !### calculate QPI-like pattern using unfolded spectral weight and write QPI result
   if(myid .eq. 0 .and. flag_get_qpi .and. flag_get_ldos .and. flag_get_unfold) then
     j=0 ; i = 0
     do igz=0, ig(3)-1
     do igy=0, ig(2)-1
     do igx=0, ig(1)-1
       do k = 1, nkp
         j = j + 1
         if( mod(j-1, nkdiv) .ne. 0 .or. j .eq. 1) then
           i = i + 1
           ldos_kk(:,i) = sw_k(:,j)
           kline_(i) = kline(k) + kline(nkp)*(igx+igy+igz)
         endif
       enddo
     enddo 
     enddo
     enddo
     nkp_ = i
     nkp_x2 = nkdiv * 4 - 3

     dk=kline(2) - kline(1)
     do i=1, nkp_x2      !nkp-1
         do j=1, nkp_x2      !nkp-1
           d = nint(kline_(j) / dk) + 1
           q(d) = kline_(j)
           do ie=1, nediv
             qpi(ie,d) = qpi(ie,d) + ldos_kk(ie,i) * ldos_kk(ie,i+d-1) * weight
           enddo
         enddo
     enddo

     ! write QPI result 
     open(999,file='QPI.'//trim(header)//'.dat',status='unknown')
     write(999,'(A,*(1x,I0))')'# QPI arround local atoms (calculated from unfolded band): ',iatom_ldos(:)
       do d =1, nkp_x2
         do ie=1, nediv
           write(999, '(F16.8, F16.8, F16.8)')q(d), e_range(ie), qpi(ie,d)
         enddo
         write(999,*)' '
       enddo
     close(999)
     write(message,'(A)')' -> QPI written: QPI.'//trim(header)//'.dat ' ; write_msg
   endif 
   !######################################################

#ifdef MPI
   call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

   !# write SLDOS with k-select (no en-select) ###########
   if(.not. flag_enselect .and. flag_kselect) then
     do ik = 1, nkp_sldos
       call system('mkdir -p  '//trim(path_sldos))
       call mpi_job_distribution_chain( nediv, nprocs, ourjob, ourjob_disp)
       do ie = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
          write(ie_str, '(F9.4)')e_range(ie)
          write(filenm_sldos,'(3A,I0,3A,I0,A)') trim(path_sldos)//'/DIDV.SLDOS.',trim(header),'.',ie,'.',trim(adjustl(ie_str)),'.KP',kidx_sldos(ik),'.dat'
          open(999+myid, file=trim(filenm_sldos), status='unknown')
          write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5)*10:($5) w p lc palette ps vari pt 7"
          write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5) w p lc palette pt 7"
          write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:($6==1?$2:1/0):3:($5) w p lc palette pt 7 # conditional plot"
          write(999+myid,'(A,I0)')'# Kpoint index: ', kidx_sldos(ik)
          write(999+myid,'(A,2(F16.8,A))')'# Spatial local density of states : EWINDOW = [',erange_i,':',erange_f,']'
          write(999+myid,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom
          write(999+myid,'(2(A,F16.8))'  )'# dE(eV)= ',e_range(1)-e_range(0), ' ,SMEARING(eV)= ', sigma
       
          ! note collinear case should be modified...
          write(999+myid,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)        ENERGY(eV)             SLDOS        ATOM_SPECIES      ATOM_NUMBER    SITE_INDEX'
          write(999+myid,'(A,F16.8)')'# ENERGY (eV): ', e_range(ie)
          do iatom = 1, PGEOM%n_atom
            do ix = 1, nx
              do iy = 1, ny
                do iz = 1, nz
                  coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                                  T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                                  T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
                 !if(PINPT%nspin .eq. 2) then
                 !  write(mypid_ldos,'(4(F16.8,1x),2(F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), PRPLT%replot_sldos_sum(iatom, :, ie), &
                 !                                                                                     PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
                 !else
                    write(999+myid,'(4(F16.8,1x), (F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), ldos5_k(iatom,ie,ik), PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
                 !endif
                enddo
              enddo
            enddo
          enddo
          write(message,'(A,I0,A,A)')' -> SLDOS at ',kidx_sldos(ik),'-th K point is written: ', trim(filenm_sldos); write_msg_all
          close(999+myid)
       enddo ! ie
#ifdef MPI
       call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif
     enddo
   endif
   !######################################################


#ifdef MPI
   call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

   ! SLDOS SELECT K and EIGEN 
   if(flag_sldos_select) then
     call system('mkdir -p  '//trim(path_sldos))
     call mpi_job_distribution_chain(nkp_sldos, nprocs, ourjob, ourjob_disp)
     do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
     do ie = 1, ne_sldos
        write(ie_str, '(F9.4)')e_range(ie)
        write(filenm_sldos,'(3A,I0,A,I0,A)') trim(path_sldos)//'/EIGEN.',trim(header),'.EN',eidx_sldos(ie),'.KP',kidx_sldos(ik),'.dat'
        open(999+myid, file=trim(filenm_sldos), status='unknown')
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5)*10:($5) w p lc palette ps vari pt 7"
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5) w p lc palette pt 7"
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:($6==1?$2:1/0):3:($5) w p lc palette pt 7 # conditional plot"
        write(999+myid,'(A,I0,A,F16.5,A,I0)')'# eigenvalue index: ',eidx_sldos(ie), ', energy level: ', E(eidx_sldos(ie),kidx_sldos(ik)),' Kpoint index: ', kidx_sldos(ik)
        write(999+myid,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom

        ! note collinear case should be modified...
        write(999+myid,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)        ENERGY(eV)             SLDOS        ATOM_SPECIES      ATOM_NUMBER    SITE_INDEX'
        write(999+myid,'(A,F16.8)')'# ENERGY (eV): ', e_range(ie)
        do iatom = 1, PGEOM%n_atom
          do ix = 1, nx
            do iy = 1, ny
              do iz = 1, nz
                coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                                T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                                T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
               !if(PINPT%nspin .eq. 2) then
               !  write(mypid_ldos,'(4(F16.8,1x),2(F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), PRPLT%replot_sldos_sum(iatom, :, ie), &
               !                                                                                     PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
               !else
                  write(999+myid,'(4(F16.8,1x), (F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), ldos4_k(iatom, ie, ik), PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
               !endif
              enddo
            enddo
          enddo
        enddo
        write(message,'(A,A)')' -> EIGEN state written: ', trim(filenm_sldos); write_msg_all
        close(999+myid)
     enddo
     enddo

     if(myid .eq. 0) then
       write(ie_str, '(F9.4)')e_range(ie)
       write(filenm_sldos,'(3A)') trim(path_sldos)//'/EIGEN.',trim(header),'.sum.dat'
       open(999+myid, file=trim(filenm_sldos), status='unknown')
       write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5)*10:($5) w p lc palette ps vari pt 7"
       write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5) w p lc palette pt 7"
       write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:($6==1?$2:1/0):3:($5) w p lc palette pt 7 # conditional plot"
       write(999+myid,'(A,*(I6))')'# k-point    index: ',kidx_sldos(:)
       write(999+myid,'(A,*(I6))')'# eigenvalue index: ',eidx_sldos(:)
       write(999+myid,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom
  
       ! note collinear case should be modified...
       write(999+myid,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)        ENERGY(eV)             SLDOS        ATOM_SPECIES      ATOM_NUMBER    SITE_INDEX'
       write(999+myid,'(A,F16.8)')'# ENERGY (eV): ', e_range(ie)
       do iatom = 1, PGEOM%n_atom
         do ix = 1, nx
           do iy = 1, ny
             do iz = 1, nz
               coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                               T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                               T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
              !if(PINPT%nspin .eq. 2) then
              !  write(mypid_ldos,'(4(F16.8,1x),2(F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), PRPLT%replot_sldos_sum(iatom, :, ie), &
              !                                                                                     PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
              !else
                 write(999+myid,'(4(F16.8,1x), (F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), ldos4_sum(iatom), PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
              !endif
             enddo
           enddo
         enddo
       enddo
       write(message,'(A,A)')' -> EIGEN state written: ', trim(filenm_sldos); write_msg_all
       close(999+myid)
     endif
   endif
#ifdef MPI
   call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

   ! write SLDOS to didv folder
   if(flag_get_sldos) then
     call system('mkdir -p  '//trim(path_sldos)) ! ./didv')
     call mpi_job_distribution_chain( nediv, nprocs, ourjob, ourjob_disp)

     do ie = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        write(ie_str, '(F9.4)')e_range(ie)
        write(filenm_sldos,'(3A,I0,3A)') trim(path_sldos)//'/DIDV.SLDOS.',trim(header),'.',ie,'.',trim(adjustl(ie_str)),'.dat'
        open(999+myid, file=trim(filenm_sldos), status='unknown')
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5)*10:($5) w p lc palette ps vari pt 7"
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5) w p lc palette pt 7"
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:($6==1?$2:1/0):3:($5) w p lc palette pt 7 # conditional plot"
        write(999+myid,'(A,2(F16.8,A))')'# Spatial local density of states : EWINDOW = [',erange_i,':',erange_f,']'
        write(999+myid,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom
        write(999+myid,'(2(A,F16.8))'  )'# dE(eV)= ',e_range(1)-e_range(0), ' ,SMEARING(eV)= ', sigma
        
        ! note collinear case should be modified...
        write(999+myid,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)        ENERGY(eV)             SLDOS        ATOM_SPECIES      ATOM_NUMBER    SITE_INDEX'
        write(999+myid,'(A,F16.8)')'# ENERGY (eV): ', e_range(ie)
        do iatom = 1, PGEOM%n_atom
          do ix = 1, nx
            do iy = 1, ny
              do iz = 1, nz
                coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                                T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                                T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
               !if(PINPT%nspin .eq. 2) then
               !  write(mypid_ldos,'(4(F16.8,1x),2(F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), PRPLT%replot_sldos_sum(iatom, :, ie), &
               !                                                                                     PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
               !else
                  write(999+myid,'(4(F16.8,1x), (F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), sldos_sum(iatom, ie), PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
               !endif
              enddo
            enddo
          enddo
        enddo
        write(message,'(A,A)')' -> SLDOS written: ', trim(filenm_sldos); write_msg_all
        close(999+myid)
     enddo
#ifdef MPI
     call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

     if(myid .eq. 0) then
        write(ie_str, '(F9.4)')e_range(0)
        write(fe_str, '(F9.4)')e_range(nediv)
        write(filenm_sldos,'(7A)') trim(path_sldos)//'/DIDV.SLDOS.sum.',trim(header),'.',trim(adjustl(ie_str)),'_',trim(adjustl(fe_str)),'.dat'
        open(999+myid, file=trim(filenm_sldos), status='unknown')
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5)*10:($5) w p lc palette ps vari pt 7"
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:2:3:($5) w p lc palette pt 7"
        write(999+myid,'(3A)')"# simple gnuplot command: splot '",trim(filenm_sldos),"' u 1:($6==1?$2:1/0):3:($5) w p lc palette pt 7 # conditional plot"
        write(999+myid,'(A,2(F16.8,A))')'# Integrated spatial local density of states : EWINDOW = [',erange_i,':',erange_f,']'
        write(999+myid,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom
        write(999+myid,'(2(A,F16.8))'  )'# dE(eV)= ',e_range(1)-e_range(0), ' ,SMEARING(eV)= ', sigma

        ! note collinear case should be modified...
        write(999+myid,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)             SLDOS        ATOM_SPECIES'
        do iatom = 1, PGEOM%n_atom
          do ix = 1, nx
            do iy = 1, ny
              do iz = 1, nz
                coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                                T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                                T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
                if(PINPT%nspin .eq. 2) then
                  if_main write(999+myid,'(3(F16.8,1x),2(F16.8,1x),I12)')coord, sum(sldos_sum(iatom, :)), PGEOM%spec(iatom)
                else
                  if_main write(999+myid,'(3(F16.8,1x), (F16.8,1x),I12)')coord, sum(sldos_sum(iatom, :)), PGEOM%spec(iatom)
                endif
              enddo
            enddo
          enddo
        enddo
        write(message,'(2A)')' -> SLDOS sum over erange is written: ', filenm_sldos  ; write_msg
        close(999+myid)
     endif
        
   endif ! sldos


#ifdef MPI
  call MPI_BARRIER(mpi_comm_earth, mpierr)
  call mpi_finish()
#endif
  stop
endprogram
