#include "alias.inc"
program get_ldos
   use parameters, only : incar, poscar, pid_geom,  zi 
   use mpi_setup
   use time
   use print_io
   implicit none
   real*8                    t_start, t_end  
   integer*4                 mpierr
   integer*4                 nediv, nkp, nemax
   integer*4                        nkp_, nkp_x2
   integer*4                 nkdiv
   integer*4                 ii, i, j, k, ie, my_k
   integer*4                 d
   integer*4                 igx, igy, igz, ig(3), ik_shift
   integer*4                 ispinor
   integer*4                 natom_ldos,  nspec, norb
   integer*4                 my_pid
   integer*4                 narg,iarg
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
   complex*16                sw_up, sw_dn
   real*8,    allocatable :: E_(:,:)
   real*8,    allocatable :: e_range(:)
   real*8,    allocatable :: ldos(:,:), dos(:)
   real*8,    allocatable :: ldos__(:,:), dos__(:)
   real*8,    allocatable :: ldos_k(  :,:), ldos2_k(:,:), dos_k(:,:)
   real*8,    allocatable :: ldos_k_(  :,:), ldos2_k_(:,:), dos_k_(:,:)
   real*8,    allocatable :: ldos_kk( :,:)
   real*8,    allocatable :: sw_k(:,:), sw2_k(:,:)
   real*8,    allocatable :: sw_k_(:,:), sw2_k_(:,:)
   real*8,    allocatable :: kk(:,:), kline(:), kline_(:)
   real*8,    allocatable :: kk_(:,:)
   real*8,    allocatable :: kp(:,:), kp_(:,:), kline_p(:), kline_p_(:)
   real*8                    gg(3)
   real*8,    allocatable :: qpi(:,:), q(:)
   real*8,    external    :: fgauss
   logical                   flag_get_ldos, flag_get_qpi, flag_ewindow, flag_phase, flag_get_unfold
   logical                   flag_wf
   character*80              fname_in, fname, header, fname_geom
   character*20,external  :: int2str
   integer*4, allocatable :: ourjob(:), ourjob_disp(:)
   type (incar )            :: PINPT
   type (poscar)            :: PGEOM
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
 
   ig = 1

   allocate(ourjob(nprocs))
   allocate(ourjob_disp(0:nprocs-1))
!  parse basic information
   ispinor = 2 ! SOC case
   fname_in = 'ldos.in'
   fname_geom = 'POSCAR-TB'
   narg = iargc()

   call read_geom(PINPT,PGEOM, fname_geom)

   if(narg .eq. 1) call getarg(1, fname_in)

   open(999, file=trim(fname_in), status='old')
     write(message,'(A)')"-- title of the calculation (without blank)"  ; write_msg
       read(999,*)header
       write(message,'(A,A   )')"   CASE = ",trim(header)  ; write_msg
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
     

     close(999)
!  setup system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! set erange
   e_range =  erange_i + (/(k, k=0,nediv-1)/) * (erange_f - erange_i)/dble(nediv - 1)

   norb = sum(PGEOM%n_orbital)
   allocate(V(norb,nemax)) ; V = 0d0
   if(flag_wf) then 
     allocate(C(norb*ispinor,nemax)); C = 0d0
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   do igz=0,ig(3)-1
   do igy=0,ig(2)-1
   do igx=0,ig(1)-1
     ik_shift=(igx+igy+igz)*nkp
     call mpi_job_distribution_chain(nkp, ourjob, ourjob_disp)
     gg = PGEOM%b_latt(:,1) * real(igx) + PGEOM%b_latt(:,2) * real(igy) + PGEOM%b_latt(:,3) * real(igz) ! G vector 
     do k = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
       my_k = k - sum(ourjob(1:myid))
       fname = './band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(k)))//'.dat'
       E(:,k) = -999d9
       if(.not. flag_wf) then
         call load_band(PGEOM, fname, nemax, norb, ne_found(k), E(:,k), V, kk(:,k), flag_get_ldos, flag_ewindow)
       elseif(flag_wf) then
         call load_band_C(PGEOM, fname, nemax, norb, ispinor, ne_found(k), E(:,k), C, V, kk(:,k), &
                          flag_get_ldos,flag_ewindow, flag_phase, gg)
       endif
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
               endif
             enddo
           endif
   
           if(flag_get_unfold) then
             sw_up = 0d0; sw_dn = 0d0
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
     if(flag_get_unfold) then
       call MPI_ALLREDUCE(  sw_k,   sw_k_, size(  sw_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
         sw_k =   sw_k_
       call MPI_ALLREDUCE(  sw2_k,   sw2_k_, size(  sw2_k),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
         sw2_k =   sw2_k_
     endif
   endif
#endif
   write(message,'(A)')" Replot DOS ... Done!"  ; write_msg

   call get_kline_dist(kk, nkp, kline)
   ! write band structure
   if(myid .eq. 0) then
     open(999,file='band_structure_TBA.total.dat',status='unknown')
     if(flag_get_ldos) write(999,'(A,*(1x,I0))')'# LDOS for atoms: ',iatom_ldos(:)
     do i = 1, nemax
       write(999,'(A,I0,A)')'# ',i, ' -th eigen'
       do k = 1,nkp
         if(flag_get_ldos) then
           write(999,'(2F16.6, *(F10.4))')kline(k), E(i,k), ldos2_k(i,k)
         else
           write(999,'(2F16.6          )')kline(k), E(i,k)
         endif
       enddo
       write(999,*)' '
       write(999,*)' '
     enddo
     write(message,'(A)')' -> Band structure (+atom projected) written: band_structure_TBA.total.dat' ; write_msg
     close(999)

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
subroutine load_band_C(PGEOM, fname, nemax, norb, ispinor, ne_found, E, C,V, kk, flag_get_ldos, flag_ewindow, flag_phase, gg)
   use parameters, only : poscar
   use mpi_setup
   use phase_factor
   use print_io
   implicit none
   character*80    fname
   character*256   inputline
   integer*4       my_pid !, myid
   integer*4       mpierr, i_continue
   integer*4       norb, nemax, ne_found, ispinor
   integer*4       nskip
   integer*4       ie, ik, nk, im
   integer*4       idummy
   real*8          kk(3), gg(3)
   real*8          E(nemax)
   real*8          V(norb,nemax)
   complex*16      C(norb*ispinor,nemax)
   complex*16      c_up, c_dn, sw_up, sw_dn
   complex*16      phase
   real*8          C_(norb*ispinor * 2)
   logical         flag_go, flag_get_ldos, flag_ewindow, flag_phase
!  complex*16      F_IJ
   integer*4,external ::  nitems
   type (poscar)            :: PGEOM
   E = -999d0 ; V = 0d0 ; C = 0d0
   nk = 1  ! for each file only one k-points has been written
   nskip = 3 ! for k-grid:3, kline:1
   my_pid = 100 + myid
   phase = (1d0,0)

   open(my_pid,file=trim(fname), status='old', iostat=i_continue)
   if(flag_ewindow) then
     read(my_pid, '(A)') inputline
     read(my_pid, '(A)') inputline
     read(my_pid, *    ) inputline, inputline, ne_found
   else
     read(my_pid, '(A)') inputline
     ne_found = nemax
   endif

   if(ne_found .eq. 0) then
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
     read(my_pid,*) kk(:)
     return
   else
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

       do ik = 1, nk
         read(my_pid,*) kk(:), E(ie), C_(:)
         do im = 1, norb
           if(flag_phase) phase = F_IJ(-(kk(:)+gg(:)), PGEOM%o_coord_cart(:,im))
           if(ispinor .eq. 2) then
             C(im,ie)      = cmplx(C_(im*4-3), C_(im*4-2)) * phase
             C(im+norb,ie) = cmplx(C_(im*4-1), C_(im*4  )) * phase
             c_up = cmplx(C_(im*4-3), C_(im*4-2))
             c_dn = cmplx(C_(im*4-1), C_(im*4  ))
             V(im,ie) = real( conjg(c_up)*c_up + conjg(c_dn)*c_dn)
           elseif(ispinor .eq. 1) then
             C(im,ie)      = cmplx(C_(im*2-1), C_(im*2)) * phase
             c_up = cmplx(C_(im*2-1), C_(im*2))
             V(im,ie) = real( conjg(c_up)*c_up )
           endif
         enddo
       enddo
     enddo
     close(my_pid)
   endif
return
endsubroutine

subroutine load_band(PGEOM, fname, nemax, norb, ne_found, E, V, kk, flag_get_ldos, flag_ewindow)
   use parameters, only : poscar
   use mpi_setup
   use print_io
   implicit none
   character*80    fname
   character*256   inputline
   integer*4       my_pid !, myid
   integer*4       mpierr, i_continue
   integer*4       norb, nemax, ne_found
   integer*4       nskip
   integer*4       ie, ik, nk
   integer*4       idummy
   real*8          kk(3)
   real*8          E(nemax)
   real*8          V(norb,nemax)
   logical         flag_go, flag_get_ldos, flag_ewindow
   integer*4,external ::  nitems
   type (poscar)            :: PGEOM
 
   E = -999d0 ; V = 0d0   
   nk = 1  ! for each file only one k-points has been written
   nskip = 3 ! for k-grid:3, kline:1
!  myid = 0
   my_pid = 100 + myid

   open(my_pid,file=trim(fname), status='old', iostat=i_continue)
   if(flag_ewindow) then
     read(my_pid, '(A)') inputline
     read(my_pid, '(A)') inputline
     read(my_pid, *    ) inputline, inputline, ne_found
   else
     read(my_pid, '(A)') inputline
     ne_found = nemax
   endif

   if(ne_found .eq. 0) then 
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
     read(my_pid,*) kk(:)    
     return
   else
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
     
       do ik = 1, nk
         if(flag_get_ldos) then
           read(my_pid,*) kk(:), E(ie), V(1:norb,ie)
         else
           read(my_pid,*) kk(:), E(ie)
         endif
       enddo
     enddo
     close(my_pid)
   endif
return
endsubroutine
function fgauss(sigma, x)
   implicit none
   real*8   sigma,sigma2
   real*8   x,xx
   real*8   fgauss
   real*8 pi, pi2
   pi = 4.d0*atan(1.d0)
   pi2=pi*2d0

   xx = x**2
   sigma2 = sigma**2

   fgauss= exp(-0.5d0*xx/sigma2)/(sigma*sqrt(pi2))

return
end function

function int2str(w) result(string)
  implicit none
! character(*), intent(out) :: string
  character*20  string
  integer*4,    intent(in)  :: w

  write(string,*) w

  return
endfunction

function nitems(string)
  implicit none
  logical blank
  integer*4 nitems,l,i
  character(*),intent(in) :: string
  nitems=0
  l=len_trim(string)
  blank = .true.
  do i=1,l
   if(string(i:i) .eq. '#') exit

   if (blank .and. string(i:i) .ne. ' ' ) then
     blank=.false.
     nitems=nitems + 1
   elseif( .not. blank .and. string(i:i) .eq. ' ') then
     blank=.true.
   endif
  enddo
  return
endfunction

subroutine get_kline_dist(kp, nkp, kline)
   implicit none
   integer*4    ik, nkp
   real*8       kline(nkp),k0(3), enorm
   real*8       kp(3,nkp)
   external     enorm

   do ik=1,nkp
     if(ik .eq. 1) then
      k0=kp(:,1)
      kline(1)=0
     else
      k0=kp(:,ik-1)
      kline(ik)=kline(ik-1)
     endif
     kline(ik)=kline(ik)+ enorm(3, kp(:,ik)-k0(:) )
   enddo

return
endsubroutine

function enorm ( n, x )
  implicit none

  integer*4 n
  real*8 x(n),enorm

  enorm = sqrt ( sum ( x(1:n) ** 2 ))
  return
end

subroutine read_geom(PINPT,PGEOM,fname)
  use parameters,  only : incar, poscar, pid_geom
! use inverse_mat, only : inv
  use do_math, only : rotate_vector, inv
  use mpi_setup
  use print_io
  implicit none
  integer*4, parameter     :: max_orb_temp = 20
  integer*4                   mpierr
  integer*4                   i_continue,nitems
  integer*4                   i,j,ii,linecount, i_dummy, i_dummy1, i_dummy2, i_dummy3
  integer*4                   iorb
  integer*4                   pos_index(20), i_index, min_pos, n_flag_sel
  real*8, allocatable      :: local_charge_(:), local_moment_(:,:)
  real*8                      t_latt_inv(3,3)
  real*8                      t_coord(3)
  character*264               inputline, str2lowcase
  character*132               inputline_dummy, dummy
  character*80                fname
  character*40                desc_str,dummy1,dummy2,dummy3
  real*8, allocatable      :: temp_orbital_sign(:,:)
  character*8, allocatable :: temp_orbital(:,:)
  character*8                 temp
  character*20,allocatable :: site_c_index_(:)
  logical,     allocatable :: flag_site_c_index_(:)
  character*10                site_index
  character*20                locpot_index
  character(*), parameter  :: func = 'read_poscar'
  logical                     flag_skip, flag_read_moment, flag_moment_cart
  external                    nitems, str2lowcase
  type (incar )            :: PINPT  
  type (poscar)            :: PGEOM

! fname        = 'POSCAR-TB'
  PGEOM%n_spec = 0
  PGEOM%flag_selective = .false.
  flag_read_moment = .false.
  flag_moment_cart = .false.
  pos_index = 0

  write(message,*)' ' ; write_msg
  write(message,*)'*- READING INPUT GEOMETRY FILE: ',trim(fname) ; write_msg
  open (pid_geom, FILE=fname,iostat=i_continue)
  linecount = 0
  ii = 0
line: do
        read(pid_geom,'(A)',iostat=i_continue) inputline
        if(i_continue<0) exit               ! end of file reached
        if(i_continue>0) then 
          write(message,*)'Unknown error reading file:',trim(fname),func ; write_msg
        endif

        if(linecount .eq. 0) then 
          call check_comment(inputline,linecount,i,flag_skip)
          if(flag_skip) linecount = linecount + 1
        else
          call check_comment(inputline,linecount,i,flag_skip)
          if(flag_skip) linecount = linecount + 1
          if (flag_skip) cycle
        endif
        linecount = linecount + 1

        ! head
         if(linecount .eq. 1) then
           PGEOM%system_name = trim(inputline)
           write(message,'(A,A)')'   SYSTEM:  ', trim(PGEOM%system_name) ; write_msg
           cycle

        ! scaling factor
         elseif(linecount .eq. 2) then
           read(inputline,*,iostat=i_continue) PGEOM%a_scale
           write(message,'(A,F15.8)')'  A_SCALE:  ',PGEOM%a_scale ; write_msg
           cycle

        ! lattice parameter
         elseif(linecount .eq. 3 ) then
           backspace(pid_geom)
           do i=1,3
             read(pid_geom,'(A)',iostat=i_continue) inputline
             read(inputline,*,iostat=i_continue) PGEOM%a_latt(1:3,i)
             write(message,'(A,i1,A,3F15.8)')'  A_LATT',i,':  ',PGEOM%a_latt(1:3,i) ; write_msg
           enddo
           call get_reci(PGEOM%b_latt(:,1), PGEOM%b_latt(:,2), PGEOM%b_latt(:,3), &
                         PGEOM%a_latt(:,1), PGEOM%a_latt(:,2), PGEOM%a_latt(:,3))
           do i=1,3
             write(message,'(A,i1,A,3F15.8)')'  B_RECI',i,':  ',PGEOM%b_latt(1:3,i) ; write_msg
           enddo
           linecount = linecount + 2
           cycle

        ! species name and number of atoms
         elseif(linecount .eq. 6 ) then
           PGEOM%n_spec=nitems(inputline)
           allocate( PGEOM%c_spec(PGEOM%n_spec), PGEOM%i_spec(PGEOM%n_spec) )
           write(message,'(A,i8)')'   N_SPEC:',PGEOM%n_spec ; write_msg
           read(inputline,*,iostat=i_continue) PGEOM%c_spec(1:PGEOM%n_spec)
           read(pid_geom,'(A)',iostat=i_continue) inputline
           linecount = linecount + 1
           call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle
           read(inputline,*,iostat=i_continue) PGEOM%i_spec(1:PGEOM%n_spec)
           PGEOM%n_atom=sum ( PGEOM%i_spec(1:PGEOM%n_spec) )

           allocate( PGEOM%spec(PGEOM%n_atom) )
           do i=1,PGEOM%n_spec
             PGEOM%spec( sum(PGEOM%i_spec(1:i)) -PGEOM%i_spec(i)+1 : sum(PGEOM%i_spec(1:i)) ) = i
           enddo

           write(message,'(A,i8)')'   N_ATOM:',PGEOM%n_atom ; write_msg
           allocate( PGEOM%a_coord(3,PGEOM%n_atom), &
                     PGEOM%a_coord_cart(3,PGEOM%n_atom), &
                     PGEOM%n_orbital(PGEOM%n_atom), &
                     local_charge_(PGEOM%n_atom*max_orb_temp), &
                     local_moment_(3,PGEOM%n_atom*max_orb_temp), &
                     site_c_index_(PGEOM%n_atom), flag_site_c_index_(PGEOM%n_atom), &
                     temp_orbital(max_orb_temp, PGEOM%n_atom), &
                     temp_orbital_sign(max_orb_temp, PGEOM%n_atom ) )
                     local_charge_ = 0d0 ! initialize as zero
                     local_moment_ = 0d0 ! initialize as zero
                     flag_site_c_index_ = .false. ! initialize as .false.
                     temp_orbital_sign = 1d0 ! initialize as unity.
           do i=1,PGEOM%n_spec
             if(i .eq. 1)then
               write(message,'(A,I2,A,A4,1x,i8)')'  SPEC',i,':',trim(PGEOM%c_spec(i)),PGEOM%i_spec(i) ; write_msg
             else
               write(message,'(A,I2,A,A4,1x,i8)')'   SPEC',i,':',trim(PGEOM%c_spec(i)),PGEOM%i_spec(i) ; write_msg
             endif
           enddo

        ! constraint and coordinate type
         elseif(linecount .eq. 8 ) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'S' .or. desc_str(1:1) .eq. 's') then 
             PGEOM%flag_selective = .true.
             write(message,'(A)')' L_CONSTR:  .TRUE.' ; write_msg
           elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then 
             PGEOM%flag_selective = .false.
             PGEOM%flag_direct=.true.
             PGEOM%flag_cartesian=.false.
             linecount = linecount + 1
             write(message,'(A)')' L_CONSTR:  .TRUE.' ; write_msg
             write(message,'(A)')' C_CRDTYP:  DIRECT' ; write_msg
           elseif(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
                  desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then 
             PGEOM%flag_selective = .false.
             PGEOM%flag_direct=.false.
             PGEOM%flag_cartesian=.true.
             linecount = linecount + 1
             write(message,'(A)')' L_CONSTR:  .FALSE.' ; write_msg
             write(message,'(A)')' C_CRDTYP:  CARTESIAN' ; write_msg

           endif
         elseif(linecount .eq. 9 ) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
              desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then 
             PGEOM%flag_direct=.false.
             PGEOM%flag_cartesian=.true.
             write(message,'(A)')' C_CRDTYP:  CARTESIAN' ; write_msg
           elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then 
             PGEOM%flag_direct=.true.
             PGEOM%flag_cartesian=.false.
             write(message,'(A)')' C_CRDTYP:  DIRECT' ; write_msg
           endif

         ! atomic coordinate & atomic orbital information
         elseif(linecount .eq. 10 ) then
           backspace(pid_geom)

           if(PGEOM%flag_selective) then
             n_flag_sel = 3
           elseif(.not. PGEOM%flag_selective) then
             n_flag_sel = 0
           endif

           do i=1,PGEOM%n_atom
             read(pid_geom,'(A)',iostat=i_continue) inputline
             linecount = linecount + 1
             call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle
             i_dummy1=index(inputline,'#')
             if(i_dummy1 .ne. 0) inputline = inputline(1:i_dummy1-1) !check comment
             inputline = str2lowcase(inputline)  

             !check 'charge', 'moment', 'site_c_index'
             pos_index(1)=index(trim(inputline),'charge') !check whether 'charge' has been set up
             locpot_index='charge'
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local_pot')
               locpot_index='local_pot' 
             endif
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local.pot')
               locpot_index='local.pot'
             endif
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local_potential')
               locpot_index='local_potential'
             endif
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local.potential')
               locpot_index='local.potential'
             endif

             pos_index(2)=index(trim(inputline),'moment') !check whether 'moment' has been set up
             site_index='site_index'
             pos_index(3)=index(trim(inputline),'site_index') !check whether 'site_index' has been set up
             if(pos_index(3) .eq. 0) then 
               site_index='site_indx'
               pos_index(3)=index(trim(inputline),'site_indx')
             endif
             if(pos_index(3) .eq. 0) then 
               site_index='site_idx'
               pos_index(3)=index(trim(inputline),'site_idx')
             endif
             if(pos_index(3) .eq. 0) then 
               site_index='site_name'
               pos_index(3)=index(trim(inputline),'site_name')
             endif
             min_pos = 999999
             do i_index = 1, 3
               if(pos_index(i_index) .ne. 0) then
                 if(pos_index(i_index) .lt. min_pos) then
                   min_pos = pos_index(i_index)
                 endif
               endif
             enddo
            
             if(min_pos .eq. 999999) then  ! if no other values are defined.
               inputline_dummy = inputline
             else
               inputline_dummy = inputline(1:min_pos-1)
             endif
             ii = ii + 1
             PGEOM%n_orbital(i) = nitems(inputline_dummy) - 3 - n_flag_sel

             if(pos_index(1) .ne. 0) then ! charge
               if(pos_index(2) .ne. 0) then 
                 call strip_off(trim(inputline), dummy, trim(locpot_index), 'moment', 1)
               elseif(pos_index(2) .eq. 0 .and. pos_index(3) .ne. 0) then
                 call strip_off(trim(inputline), dummy, trim(locpot_index), trim(site_index), 1)
               else
                 call strip_off(trim(inputline), dummy, trim(locpot_index), '', 2)
               endif
               i_dummy = nitems(dummy)
               if(i_dummy .ne. PGEOM%n_orbital(i)) then
                 write(message,'(A,I6)')'  !WARNING! Charge setting is inproper. Number of items should be same as N_ORBITAL(i). iatom=',i ; write_msg
                 write(message,'(A)')   '  !WARNING! Please check GFILE. Exit program...' ; write_msg
                 stop
               else
                 read(dummy,*)local_charge_(ii:ii+PGEOM%n_orbital(i)-1)
               endif
             endif ! charge

             if(pos_index(2) .ne. 0) then ! moment
               if(pos_index(3) .ne. 0) then
                 if(index(trim(inputline),'moment.r') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.r', trim(site_index), 1)
                 elseif(index(trim(inputline),'moment.c') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.c', trim(site_index), 1)
                   flag_moment_cart = .true.
                 else
                   call strip_off(trim(inputline), dummy, 'moment', trim(site_index), 1)
                 endif
               else
                 if(index(trim(inputline),'moment.r') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.r', '', 2)
                 elseif(index(trim(inputline),'moment.c') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.c', '', 2)
                   flag_moment_cart = .true.
                 else
                   call strip_off(trim(inputline), dummy, 'moment', '', 2)
                 endif
               endif

               i_dummy = nitems(dummy)
               if(PINPT%flag_collinear) then
                 if(i_dummy .ne. PGEOM%n_orbital(i)) then
                   write(message,'(A)')'  !WARNING! Moment setting is inproper. Number of items should be same as N_ORBITAL(i) in the collinear setting.' ; write_msg
                   write(message,'(A,I6)')'  !WARNING! Please check GFILE. Exit program...  iatom=',i ; write_msg
                   stop
                 else 
                   read(dummy,*)local_moment_(1,ii:ii+PGEOM%n_orbital(i)-1)
                 endif
               endif
               if(PINPT%flag_noncollinear) then
                 if(i_dummy .ne. PGEOM%n_orbital(i)*3) then
                   write(message,'(A)')'  !WARNING! Moment setting is inproper. Number of items should be same as N_ORBITAL(i)*3 in the non-collinear setting.' ; write_msg 
                   write(message,'(A,I6)')'  !WARNING! N_ORBITAL(i) * 3 = number of items (M, theta, phi) or (Mx,My,Mz). Please check GFILE. Exit program...  iatom=',i ; write_msg
                   stop
                 else
                   read(dummy,*)((local_moment_(j,i_dummy),j=1,3),i_dummy = ii, ii+PGEOM%n_orbital(i)-1)
                 endif
               endif ! moment
             endif

             if(pos_index(3) .ne. 0) then ! site_index 
               !if site_index is predefined in the POSCAR with Site_index tag: flag_site_cindex = .true. and use it as site_cindex
               call strip_off(trim(inputline), dummy, trim(site_index), '', 2)
               i_dummy = nitems(dummy)
               if(i_dummy .ne. 1) then
                 write(message,'(A)')'  !WARNING! Site_index setting is inproper. Number of items should be one and the data type should be character(20)' ; write_msg
                 write(message,'(A,I6)')'  !WARNING! Please check GFILE. Exit program...  iatom=',i ; write_msg
                 stop
               else
                 read(dummy,*)site_c_index_(i)
                 flag_site_c_index_(i) = .true.
               endif
             else
               !if site_index is not predefined in the POSCAR with Site_index tag: flag_site_cindex = .false. and use atom_name+atom_number as site_cindex
               dummy2 = PGEOM%c_spec(PGEOM%spec(i))
               if( i .lt. 10) then
                 write(dummy3,'(I1)') i
               elseif( i .ge. 10 .and. i .lt. 100) then
                 write(dummy3,'(I2)') i
               elseif( i .ge. 100 .and. i .lt. 1000) then
                 write(dummy3,'(I3)') i
               elseif( i .ge. 1000 .and. i .lt. 10000) then
                 write(dummy3,'(I4)') i
               elseif( i .ge. 10000 .and. i .lt. 100000) then
                 write(dummy3,'(I5)') i
               elseif( i .ge. 100000 .and. i .lt. 1000000) then
                 write(dummy3,'(I6)') i
               endif
               i_dummy2 = len_trim(dummy2)
               i_dummy3 = len_trim(dummy3)
               write(site_c_index_(i),'(2A)') dummy2(1:i_dummy2),dummy3(1:i_dummy3)
               site_c_index_(i) = str2lowcase(site_c_index_(i))
             endif !site_index

             ii = ii+PGEOM%n_orbital(i)-1

             ! read atomic coordinate & atomic orbital & magnetic moment info (if 'moment' tag is provided)              
             if(PGEOM%n_orbital(i) .eq. 0) then
               read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i)
               temp_orbital(1,i)='_na_'
             elseif(PGEOM%n_orbital(i) .ge. 1) then
               if(PGEOM%flag_selective) then
                 read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i), desc_str, desc_str, desc_str, &
                                                     temp_orbital(1:PGEOM%n_orbital(i),i)
               else
                 read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i), &
                                                     temp_orbital(1:PGEOM%n_orbital(i),i)
                 do i_dummy = 1, PGEOM%n_orbital(i)
                   temp = trim(temp_orbital(i_dummy,i))
                   if(temp(1:1) .eq. '-') then
                     temp_orbital(i_dummy,i) = temp(2:)
                     temp_orbital_sign(i_dummy,i) = -1d0
                   endif
                 enddo
               endif
             endif

           enddo !n_atom
           PGEOM%neig=sum(PGEOM%n_orbital(1:PGEOM%n_atom))
           PGEOM%neig_total = PGEOM%neig * PINPT%ispin
           PGEOM%nbasis = PGEOM%neig 

           if(PGEOM%neig  == 0) then
             write(message,'(A)')'  !! Check geometry input file. atomic orbital is not asigned!' ; write_msg
           elseif(PGEOM%neig >= 1) then
             write(message,'(A,i8)')'  N_ORBIT:',PGEOM%neig ; write_msg
             PGEOM%max_orb=maxval( PGEOM%n_orbital(:) )
             allocate( PGEOM%c_orbital(PGEOM%max_orb,PGEOM%n_atom) ) ; PGEOM%c_orbital = '_na_'
             allocate( PGEOM%orb_sign(PGEOM%max_orb,PGEOM%n_atom) )  ; PGEOM%orb_sign = 1d0
             PGEOM%c_orbital(1:PGEOM%max_orb,1:PGEOM%n_atom) = temp_orbital(1:PGEOM%max_orb,1:PGEOM%n_atom)
             PGEOM%orb_sign(1:PGEOM%max_orb,1:PGEOM%n_atom)  = temp_orbital_sign(1:PGEOM%max_orb,1:PGEOM%n_atom)
             do i=1,PGEOM%n_atom
               if(PGEOM%n_orbital(i) .eq. 0) then
                 if(PINPT%flag_report_geom) then
                   write(message,'(A,I4,A,I3,2x,10A7)')' ATOM',i,': ',PGEOM%n_orbital(i), PGEOM%c_orbital(1,i) ; write_msg
                 endif

               elseif(PGEOM%n_orbital(i) .gt. 0) then
                 if(PINPT%flag_report_geom) then
                   write(message,'(A,I4,A,I3,2x,10A7)')' ATOM',i,': ',PGEOM%n_orbital(i), PGEOM%c_orbital(1:PGEOM%n_orbital(i),i) ; write_msg
                   write(message,'(A,A20)'            )' SITE_IDX:   ',site_c_index_(i) ; write_msg
                 endif

                 if(PINPT%flag_local_charge) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(PINPT%flag_report_geom) then
                     write(message,'(A,*(F10.4))')'   CHARGE:   ',local_charge_(i_dummy:i_dummy1) ; write_msg
                   endif
                 endif

                 if(PINPT%flag_collinear) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(PINPT%flag_report_geom) then
                     write(message,'(A,*(F10.4))')'   MAGMOM:   ',local_moment_(1,i_dummy:i_dummy1) ; write_msg
                   endif
                 elseif(PINPT%flag_noncollinear) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(flag_moment_cart)then
                     if(PINPT%flag_report_geom) then
                       write(message,'(A,*(3F7.3,2x))')'   MAGMOM: (Mx,My,Mz) ',(local_moment_(1:3,i_dummy2),i_dummy2=i_dummy,i_dummy1) ; write_msg
                     endif
                   else
                     if(PINPT%flag_report_geom) then
                       write(message,'(A,*(3F7.3,2x))')'   MAGMOM: (M,theta,phi) ',(local_moment_(1:3,i_dummy2),i_dummy2=i_dummy,i_dummy1) ; write_msg
                     endif
                   endif
                 endif

               endif
             enddo
           elseif(PGEOM%neig < 0)then
             write(message,'(A)')'  !! Check geometry input file. negative number of atomic orbitals ??' ; write_msg
           endif
         endif ! linecount

         if (  i_continue .ne. 0 ) cycle  ! skip empty line 

      enddo line
  
! call check_sanity(PINOT,PGEOM) !!!!! for future work
  if(PGEOM%flag_cartesian) then ! cartesian
    t_latt_inv = inv(PGEOM%a_latt)
    ii = 0
    allocate(PGEOM%o_coord(3,PGEOM%neig))
    allocate(PGEOM%o_coord_cart(3,PGEOM%neig))
    do i = 1, PGEOM%n_atom
      t_coord=PGEOM%a_coord(:,i)
      PGEOM%a_coord(1,i) = dot_product(t_latt_inv(1,:), t_coord(:)) 
      PGEOM%a_coord(2,i) = dot_product(t_latt_inv(2,:), t_coord(:))
      PGEOM%a_coord(3,i) = dot_product(t_latt_inv(3,:), t_coord(:))
      PGEOM%a_coord(:,i) = PGEOM%a_coord(:,i) - int(PGEOM%a_coord(:,i))
      PGEOM%a_coord_cart(:,i) = PGEOM%a_coord(1,i) * PGEOM%a_latt(:,1) + &
                                PGEOM%a_coord(2,i) * PGEOM%a_latt(:,2) + &
                                PGEOM%a_coord(3,i) * PGEOM%a_latt(:,3)
      do iorb=1,PGEOM%n_orbital(i)
        ii = ii + 1
        PGEOM%o_coord(:,ii) = PGEOM%a_coord(:,i)
        PGEOM%o_coord_cart(:,ii) = PGEOM%o_coord(1,ii) * PGEOM%a_latt(:,1) + &
                                   PGEOM%o_coord(2,ii) * PGEOM%a_latt(:,2) + &
                                   PGEOM%o_coord(3,ii) * PGEOM%a_latt(:,3)
      enddo
    enddo
  else ! direct
    allocate(PGEOM%o_coord(3,PGEOM%neig))
    allocate(PGEOM%o_coord_cart(3,PGEOM%neig))
    ii = 0
    do i = 1, PGEOM%n_atom
      PGEOM%a_coord_cart(:,i) = PGEOM%a_coord(1,i) * PGEOM%a_latt(:,1) + &
                                PGEOM%a_coord(2,i) * PGEOM%a_latt(:,2) + &
                                PGEOM%a_coord(3,i) * PGEOM%a_latt(:,3)
      do iorb=1,PGEOM%n_orbital(i)
        ii = ii + 1
        PGEOM%o_coord(:,ii) = PGEOM%a_coord(:,i)
        PGEOM%o_coord_cart(:,ii) = PGEOM%o_coord(1,ii) * PGEOM%a_latt(:,1) + &
                                   PGEOM%o_coord(2,ii) * PGEOM%a_latt(:,2) + &
                                   PGEOM%o_coord(3,ii) * PGEOM%a_latt(:,3) 
      enddo
    enddo
  endif

  ! store local charge
! allocate(NN_TABLE%local_charge(PGEOM%neig)) ! local charge for each atomic orbital basis
! NN_TABLE%local_charge = 0d0
! if(PINPT%flag_local_charge) then
!   NN_TABLE%local_charge(:) = local_charge_(1:PGEOM%neig)
! endif

! ! store site_index
! allocate(NN_TABLE%site_cindex(PGEOM%n_atom)) ! site_index for each atomic site
! allocate(NN_TABLE%flag_site_cindex(PGEOM%n_atom)) ! flag for site_index for each atomic site
! NN_TABLE%site_cindex = site_c_index_
! NN_TABLE%flag_site_cindex =flag_site_c_index_
! allocate(PGEOM%site_cindex(PGEOM%n_atom))
! PGEOM%site_cindex = site_c_index_

  if (linecount == 0) then
    write(message,*)'Attention - empty input file: ',trim(fname),' , ',func ; write_msg
#ifdef MPI
    call MPI_Abort(mpi_comm_earth, 0, mpierr)
#else
    stop
#endif
  endif
  close(pid_geom)


  write(message,*)'*- END READING GEOMETRY FILE ---------------------' ; write_msg
  write(message,*)' ' ; write_msg
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

