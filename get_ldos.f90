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
   integer*4                 nkdiv
   integer*4                 ii, i, j, k, ie, my_k, im
   integer*4                 d
   integer*4                 igx, igy, igz, ig(3), ik_shift
   integer*4                 ispinor
   integer*4                 natom_ldos,  nspec, norb
   integer*4                 my_pid
   integer*4                 narg,iarg
   integer*4                 lmmax
   integer*4, allocatable :: ne_found(:)
   integer*4, allocatable :: ne_found_(:)
   integer*4, allocatable :: iatom_ldos(:)
   real*8                    erange_i, erange_f, sigma
   real*8                    dos_, ldos_, ldos2_
   real*8                    weight
   real*8                    dk
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
   real*8,    allocatable :: ldos_kk( :,:)
   real*8,    allocatable :: sw_k(:,:), sw2_k(:,:)
   real*8,    allocatable :: sw_k_(:,:), sw2_k_(:,:)
   real*8,    allocatable :: sw3_k(:,:), sw3_k_(:,:)
   real*8,    allocatable :: kk(:,:), kline(:), kline_(:)
   real*8,    allocatable :: kk_(:,:)
   real*8,    allocatable :: kp(:,:), kp_(:,:), kline_p(:), kline_p_(:)
   real*8                    gg(3)
   real*8,    allocatable :: qpi(:,:), q(:)
   logical                   flag_get_ldos, flag_get_qpi, flag_ewindow, flag_phase, flag_get_unfold
   logical                   flag_wf, flag_spectral_band, flag_formatted
   character*80              fname_in, fname, header, fname_geom
   character*40              ftype
   character*20,external  :: int2str
   integer*4, allocatable :: ourjob(:), ourjob_disp(:)
   type (incar )            :: PINPT
   type (poscar)            :: PGEOM
   type (hopping)           :: NN_TABLE
   character*40              fnamelog
   fnamelog = 'TBFIT_LDOS2.log'
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

   ig = 1

   allocate(ourjob(nprocs))
   allocate(ourjob_disp(0:nprocs-1))
!  parse basic information
   ispinor = 2 ! SOC case
   fname_in = 'ldos.in'
  !fname_geom = 'POSCAR-TB'
   narg = iargc()

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
       allocate(E(nemax,nkp)) ; E = 0d0
       allocate(E_(nemax,nkp)) ; E_ = 0d0
       allocate(ne_found(nkp)); ne_found = 0
       allocate(ne_found_(nkp)); ne_found_ = 0
     write(message,'(A)')"-- NKDIV ?"  ; write_msg
       read(999,*)nkdiv
       write(message,'(A,*(I5))')"   NKDIV = ",nkdiv  ; write_msg
     write(message,'(A)')"-- EWINDOW ?"  ; write_msg
       read(999,*)flag_ewindow
       write(message,'(A,(L1))')"   FLAG EWINDOW= ", flag_get_ldos  ; write_msg
     write(message,'(A)')"-- Gaussian smearing?"  ; write_msg
       read(999,*)sigma
       write(message,'(A,(F10.4),A)')"   SMEARING = ",sigma, ' (eV)'  ; write_msg
     write(message,'(A)')"-- Number of division? (NEDIV)"  ; write_msg
       read(999,*)nediv
       write(message,'(A,(I6))')"   NEDIV = ",nediv  ; write_msg
       allocate( dos_k(nediv,nkp)) ; dos_k = 0d0
       allocate( dos_k_(nediv,nkp)) ; dos_k_ = 0d0
       allocate( dos(nediv)) ; dos = 0d0
       allocate( dos__(nediv)) ; dos__ = 0d0
     write(message,'(A)')"-- Energy range for dos plot (erange_i  erange_f)"  ; write_msg
       read(999,*)erange_i, erange_f
       write(message,'(A,F10.4,A,F10.4)')"   ERANGE = ",erange_i,' : ',erange_f  ; write_msg
     write(message,'(A)')"-- GET LDOS ? (evaluate LDOS from band_structure_TBA.dat with specified energy range)"  ; write_msg
       read(999,*)flag_get_ldos
       write(message,'(A,(L1))')"   GET_LDOS = ", flag_get_ldos  ; write_msg
     if(flag_get_ldos) then
       allocate(ldos_k(      nediv,nkp)) ; ldos_k = 0d0
       allocate(ldos_k_(      nediv,nkp)) ; ldos_k_ = 0d0
       allocate(ldos2_k(nemax,nkp)) ; ldos2_k = 0d0
       allocate(ldos2_k_(nemax,nkp)) ; ldos2_k_ = 0d0
       allocate(ldos3_k(nemax,lmmax,nkp)) ; ldos3_k = 0d0
       allocate(ldos3_k_(nemax,lmmax,nkp)) ; ldos3_k_ = 0d0
       allocate(ldos(PGEOM%n_atom,nediv)) ; ldos = 0d0
       allocate(ldos__(PGEOM%n_atom,nediv)) ; ldos__ = 0d0
       write(message,'(A)')"-- Read wavefunction (wf, .true.) or density (rh, .false.)? Logical value (.true. or .false.)"  ; write_msg
         read(999,*)flag_wf
         if(flag_wf) then
           write(message,'(A)')'   WAVEFUNCTION |psi_nk> READ = .TRUE.'  ; write_msg
         elseif(.not. flag_wf) then
           write(message,'(A)')'   RHO READ = .TRUE.'  ; write_msg
         endif
       write(message,'(A)')"-- Number of atoms for LDOS plot "  ; write_msg
         read(999,*)natom_ldos
         write(message,'(A,(I6))')"   NATOM_LDOS = ",natom_ldos  ; write_msg
         allocate(iatom_ldos(natom_ldos))
       write(message,'(A)')"-- Atom index for ldos plot"  ; write_msg
         read(999,*)iatom_ldos(1:natom_ldos)
         write(message,'(A,*(I5))')"   LDOS ATOM INDEX = ",iatom_ldos  ; write_msg
     endif
     write(message,'(A)')"-- GET UNFOLD? "  ; write_msg
       read(999,*)flag_get_unfold
       write(message,'(A,L)')"   GET UNFOLD = ", flag_get_unfold  ; write_msg
     if(flag_get_unfold) then
       write(message,'(A)')"-- Multiply phase factor to each wavefunction coefficient?"  ; write_msg
         read(999,*)flag_phase
         write(message,'(A,L)')"   SET PHASE = ", flag_phase  ; write_msg
       write(message,'(A)')"-- G-vector direction: ngx, ngy, ngz "  ; write_msg
         read(999,*)ig(1:3)
         write(message,'(A,3(I6))')"   G_DIRECT= ",ig(:)  ; write_msg
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
       write(message,'(A)')"-- Weight "  ; write_msg
         read(999,*)weight
         write(message,'(A,F10.4        )')"   WEIGHT = ",weight  ; write_msg
     endif
     write(message,'(A)')"-- Spectral band plot ? "  ; write_msg
       read(999,*)flag_spectral_band
       write(message,'(A,(L1))')"   SPECTRAL_BAND = ", flag_spectral_band  ; write_msg
     if(flag_spectral_band) then
       allocate(sw3_k (       nemax,nkp                  )) ;   sw3_k = 0d0
       allocate(sw3_k_(       nemax,nkp                  )) ;   sw3_k_= 0d0
     endif
     close(999)

!  setup system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! set erange
   e_range =  erange_i + (/(k, k=0,nediv-1)/) * (erange_f - erange_i)/dble(nediv - 1)

   norb = sum(PGEOM%n_orbital)
   allocate(V(norb,nemax)) ; V = 0d0
   allocate(C(norb*ispinor,nemax)) ; C = 0d0
!  if(flag_wf) then 
!    allocate(C(norb*ispinor,nemax)); C = 0d0
!  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do igz=0,ig(3)-1
   do igy=0,ig(2)-1
   do igx=0,ig(1)-1
     ik_shift=(igx+igy+igz)*nkp
     call mpi_job_distribution_chain(nkp, nprocs, ourjob, ourjob_disp)
     gg = PGEOM%b_latt(:,1) * real(igx) + PGEOM%b_latt(:,2) * real(igy) + PGEOM%b_latt(:,3) * real(igz) ! G vector 
     do k = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
       my_k = k - sum(ourjob(1:myid))
       if(flag_formatted) then ! need to extend toward spin polarized case in the future. 18. March. 2021, HJ Kim
         fname = './band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(k)))//'.dat'
       else
         fname = './band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(k)))//'.bin'
       endif
       E(:,k) = -999d9
       call load_band_singlek(PGEOM, fname, nemax, norb, ispinor, E(:,k), C, V, ne_found(k), kk(:,k), gg, &
                              flag_get_ldos, flag_wf, flag_ewindow, flag_phase, flag_formatted)

!      if(.not. flag_wf) then
!        call load_band(PGEOM, fname, nemax, norb, ne_found(k), E(:,k), V, kk(:,k), flag_get_ldos, flag_ewindow, flag_formatted)
!      elseif(flag_wf) then
!        call load_band_C(PGEOM, fname, nemax, norb, ispinor, ne_found(k), E(:,k), C, V, kk(:,k), &
!                         flag_get_ldos,flag_ewindow, flag_phase, gg, flag_formatted)
!      endif

       do ie = 1, nediv
         do i = 1, ne_found(k)
           dos_ = fgauss(sigma, e_range(ie) - E(i,k)) / real(nkp)
           if( (igx+igy+igz) .eq. 0 ) then
             dos(ie) = dos(ie) + dos_
             dos_k(ie,k) = dos_k(ie,k) + dos_
           endif

           if(flag_get_ldos .and. (igx+igy+igz) .eq. 0) then
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

           if(flag_get_unfold) then
             sw_up = 0d0; sw_dn = 0d0; sw_real = 0d0
             do j = 1, natom_ldos
               ii = sum(PGEOM%n_orbital(1:iatom_ldos(j))) - PGEOM%n_orbital(iatom_ldos(j)) + 1 ! starting orbital index for atom iatom_ldos(j)
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

           elseif(flag_spectral_band .and. .not. flag_get_unfold) then
             write(message,'(A)')' !WARN! You have requested band unfolding and spectral band plot simulatneously.'
             write(message,'(A)')' !WARN! Current version does not support it. Please operate seperately. Exit program...'
             kill_job
           endif

         enddo
       enddo

       if(flag_get_unfold) then
         write(message,'(A,F10.4,A,3I3)')" STAT KP: ", real(my_k)/real(sum(ourjob(1:myid+1)))*100,' %, ig(1),ig(2),ig(3) = ',igx+1, igy+1, igz+1 ; write_msg
       else
         write(message,'(A,F10.4,A    )')" STAT KP: ", real(my_k)/real(sum(ourjob(1:myid+1)))*100,' %' ; write_msg
       endif
     enddo
   enddo
   enddo
   enddo
#ifdef MPI
   call MPI_ALLREDUCE(E, E_,size(E),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   E = E_
   call MPI_ALLREDUCE(ne_found, ne_found_, size(ne_found), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
   ne_found = ne_found_
   call MPI_ALLREDUCE(kk, kk_, size(kk), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   kk = kk_
   if(flag_get_ldos) then
     call MPI_ALLREDUCE(dos, dos__,size(dos), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     dos = dos__
     call MPI_ALLREDUCE(dos_k,dos_k_,size(dos_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     dos_k = dos_k_
     call MPI_ALLREDUCE(ldos, ldos__,size(ldos), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos = ldos__
     call MPI_ALLREDUCE(ldos_k, ldos_k_, size(ldos_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos_k = ldos_k_
     call MPI_ALLREDUCE(ldos2_k, ldos2_k_, size(ldos2_k),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos2_k = ldos2_k_
     call MPI_ALLREDUCE(ldos3_k, ldos3_k_, size(ldos3_k),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos3_k = ldos3_k_
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
   endif
#endif
   write(message,'(A)')" Replot DOS ... Done!"  ; write_msg

   call get_kline_dist(kk, nkp, kline)
   ! write band structure
   if(myid .eq. 0) then
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

     if(flag_spectral_band) then
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
    
     if(flag_get_unfold) then  ! the result of this routine is only valid if you just draw 1D k-path so that applying gg to kk is just valid 
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

     ! write DOS
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

     ! calculate spectral weight (band unfolding)
     if(flag_get_unfold) then 
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

     ! calculate QPI-like pattern using unfolded spectral weight
     if(flag_get_qpi .and. flag_get_ldos .and. flag_get_unfold) then
      !j=0
      !do k=1,nkp*ig(1)*ig(2)*ig(3)
      !  if(mod(k-1,nkdiv) .ne. 0 .or. k .eq. 1) then
      !    j = j+1
      !    ldos_kk(:,j) = sw_k(:,k)
      !    kline_(j) = kline(k)
      !  endif
      !enddo
      !nkp_ = j

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

       open(999,file='QPI.'//trim(header)//'.dat',status='unknown')
       write(999,'(A,*(1x,I0))')'# QPI arround local atoms (calculated from unfolded band): ',iatom_ldos(:)
         do d =1, nkp_x2
           do ie=1, nediv
             write(999, '(F16.8, F16.8, F16.8)')q(d), e_range(ie), qpi(ie,d)
           enddo
           write(999,*)' '
        !  write(999,*)' '
         enddo
       close(999)
       write(message,'(A)')' -> QPI written: QPI.'//trim(header)//'.dat ' ; write_msg

     endif


   endif


#ifdef MPI
  call MPI_BARRIER(mpi_comm_earth, mpierr)
  call mpi_finish()
#endif
  stop
endprogram
! only applicable routine for PRTSEPK .TRUE. 
!ubroutine load_band_C(PGEOM, fname, nemax, norb, ispinor, ne_found, E, C,V, kk, flag_get_ldos, flag_ewindow, flag_phase, gg, flag_formatted)
!  use parameters, only : poscar
!  use mpi_setup
!  use phase_factor
!  use print_io
!  implicit none
!  character*80    fname
!  character*256   inputline
!  integer*4       my_pid !, myid
!  integer*4       mpierr, i_continue
!  integer*4       norb, nemax, ne_found, ispinor
!  integer*4       nskip
!  integer*4       ie, ik, nk, im
!  integer*4       idummy
!  real*8          kk(3), gg(3)
!  real*8          E(nemax)
!  real*8          V(norb,nemax)
!  complex*16      C(norb*ispinor,nemax)
!  complex*16      c_up, c_dn, sw_up, sw_dn
!  complex*16      phase
!  real*8          C_(norb*ispinor * 2)
!  logical         flag_go, flag_get_ldos, flag_ewindow, flag_phase
!  logical         flag_formatted
!  complex*16      F_IJ
!  integer*4,external ::  nitems
!  type (poscar)            :: PGEOM
!  E = -999d0 ; V = 0d0 ; C = 0d0
!  nk = 1  ! for each file only one k-points has been written
!  nskip = 3 ! for k-grid:3, kline:1
!  my_pid = 100 + myid
!  phase = (1d0,0)

!  open(my_pid,file=trim(fname), status='old', iostat=i_continue)
!  if(flag_ewindow) then
!    read(my_pid, '(A)') inputline
!    read(my_pid, '(A)') inputline
!    read(my_pid, *    ) inputline, inputline, ne_found
!  else
!    read(my_pid, '(A)') inputline
!    ne_found = nemax
!  endif

!  if(ne_found .eq. 0) then
!    flag_go = .false.
!    do while (.not.flag_go)
!      read(my_pid,'(A)')inputline
!      idummy = nitems(inputline)
!      if(idummy .le. 0) then
!        flag_go = .false.
!      elseif(idummy .gt. 0) then
!        flag_go = .true.
!        backspace(my_pid)
!      endif
!    enddo
!    read(my_pid,*) kk(:)
!    return
!  else
!    do ie = 1, ne_found
!      flag_go = .false.
!      do while (.not.flag_go)
!        read(my_pid,'(A)')inputline
!        idummy = nitems(inputline)
!        if(idummy .le. 0) then
!          flag_go = .false.
!        elseif(idummy .gt. 0) then
!          flag_go = .true.
!          backspace(my_pid)
!        endif
!      enddo

!      do ik = 1, nk
!        read(my_pid,*) kk(:), E(ie), C_(:)
!        do im = 1, norb
!          if(flag_phase) phase = F_IJ(-(kk(:)+gg(:)), PGEOM%o_coord_cart(:,im))
!          if(ispinor .eq. 2) then
!            C(im,ie)      = cmplx(C_(im*4-3), C_(im*4-2)) * phase
!            C(im+norb,ie) = cmplx(C_(im*4-1), C_(im*4  )) * phase
!            c_up = cmplx(C_(im*4-3), C_(im*4-2))
!            c_dn = cmplx(C_(im*4-1), C_(im*4  ))
!            V(im,ie) = real( conjg(c_up)*c_up + conjg(c_dn)*c_dn)
!          elseif(ispinor .eq. 1) then
!            C(im,ie)      = cmplx(C_(im*2-1), C_(im*2)) * phase
!            c_up = cmplx(C_(im*2-1), C_(im*2))
!            V(im,ie) = real( conjg(c_up)*c_up )
!          endif
!        enddo
!      enddo
!    enddo
!    close(my_pid)
!  endif
!eturn
!ndsubroutine

!ubroutine load_band(PGEOM, fname, nemax, norb, ne_found, E, V, kk, flag_get_ldos, flag_ewindow, flag_formatted)
!  use parameters, only : poscar
!  use mpi_setup
!  use print_io
!  implicit none
!  character*80    fname
!  character*256   inputline
!  integer*4       my_pid !, myid
!  integer*4       mpierr, i_continue
!  integer*4       norb, nemax, ne_found
!  integer*4       nskip
!  integer*4       ie, ik, nk
!  integer*4       idummy
!  real*8          kk(3)
!  real*8          E(nemax)
!  real*8          V(norb,nemax)
!  logical         flag_go, flag_get_ldos, flag_ewindow, flag_formatted
!  integer*4,external ::  nitems
!  type (poscar)            :: PGEOM
!
!  E = -999d0 ; V = 0d0   
!  nk = 1  ! for each file only one k-points has been written
!  nskip = 3 ! for k-grid:3, kline:1
!  my_pid = 100 + myid

!  if(flag_formatted) then
!    open(my_pid,file=trim(fname), status='old', iostat=i_continue)
!    if(flag_ewindow) then
!      read(my_pid, '(A)') inputline
!      read(my_pid, '(A)') inputline
!      read(my_pid, *    ) inputline, inputline, ne_found
!    else
!      read(my_pid, '(A)') inputline
!      ne_found = nemax
!    endif
! 
!    if(ne_found .eq. 0) then 
!      flag_go = .false.
!      do while (.not.flag_go)
!        read(my_pid,'(A)')inputline
!        idummy = nitems(inputline)
!        if(idummy .le. 0) then
!          flag_go = .false.
!        elseif(idummy .gt. 0) then
!          flag_go = .true.
!          backspace(my_pid)
!        endif
!      enddo
!      read(my_pid,*) kk(:)    
!      close(my_pid)
!      return
!    else
!      do ie = 1, ne_found
!        flag_go = .false.
!        do while (.not.flag_go)
!          read(my_pid,'(A)')inputline
!          idummy = nitems(inputline)
!          if(idummy .le. 0) then
!            flag_go = .false.
!          elseif(idummy .gt. 0) then
!            flag_go = .true.
!            backspace(my_pid)
!          endif
!        enddo
!      
!        do ik = 1, nk
!          if(flag_get_ldos) then
!            read(my_pid,*) kk(:), E(ie), V(1:norb,ie)
!          else
!            read(my_pid,*) kk(:), E(ie)
!          endif
!        enddo
!      enddo
!      close(my_pid)
!    endif

!  elseif(.not. flag_formatted) then
!!!!! CHECK FROM HERE!!
!    open(my_pid,file=trim(fname), status='old', form='unformatted', iostat=i_continue)
!    read(my_pid) ikmode, flag_vector_, flag_single_, flag_erange_, flag_sparse_, nbasis_, nk_, nband_, &
!                 ispin_, nspin_, ispinor_, c_mode_
!    if(.not. flag_exit) then
!      if(flag_erange_) then
!        read(my_pid) flag_erange_, init_erange_, fina_erange_
!      else
!        read(pid) flag_erange_
!      endif
!      if(flag_sparse_) then
!        read(pid) flag_sparse_, emin_, emax_, nemax_
!      else
!        read(pid) flag_sparse_
!      endif

!      ! read main wavefunction information
!      if(.not. flag_single_) then
!        do ik = 1, nk
!          read(pid) ne_found(ik), kpoint_(1:ikmode,ik), (E(ie,ik),ie=1,ne_found(ik))
!        enddo
!      elseif(flag_single_) then
!        do ik = 1, nk
!          read(pid) ne_found(ik), kpoint4_(1:ikmode,ik), (E4(ie,ik),ie=1,ne_found(ik))
!        enddo
!        kpoint_ = kpoint4_
!        E  = E4
!      endif

!    endif
!  endif
!eturn
!ndsubroutine

