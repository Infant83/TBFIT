#include "alias.inc"
subroutine get_dos(NN_TABLE, PINPT, PINPT_DOS, PGEOM, PKPTS)
   use parameters  
   use mpi_setup
   use time
   use memory
   use print_io
   implicit none
   type(hopping)           :: NN_TABLE
   type(dos    )           :: PINPT_DOS
   type(incar  )           :: PINPT
   type(poscar )           :: PGEOM
   type(kpoints)           :: PKPTS
   integer*4                  i,k,ie,nkpoint,nparam,neig, nediv, ispin, nspin, ispinor
   integer*4                  my_ie
   integer*4                  iadd, iaddk, inc, inck
   integer*4                  ik,nk1,nk2,nk3
   integer*4                  iband, fband, nband, nband_found
   integer*4                  is
   integer*4                  im, imatrix, iatom, ia, ii, mm
   integer*4                  mpierr
   real*8                     kshift(3)
   real*8                     e_range(PINPT_DOS%dos_nediv)
   real*8                     emax, emin
   real*8                     sigma, x
   integer*4,  allocatable :: ne_found(:,:) 
   real*8,     allocatable :: kpoint(:,:), kpoint_reci(:,:)
   real*8,     allocatable :: param(:), param_const(:,:)
   real*8,     allocatable :: dos_tot(:,:), dos_ ! nspin,nediv
   real*8,     allocatable :: ldos_tot(:,:,:,:)                  ! nbasis+tot, dos_natom_ldos, nediv,spin
   real*8,     allocatable :: E(:,:), E_(:,:)
   complex*16, allocatable :: V(:,:,:), V_(:,:,:)
   complex*16, allocatable :: myV(:,:) !nbasis,nband
   real*8                     rho
   real*8                     a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
   real*8                     b2xb3(3),bzvol,dkv
   real*8                     fgauss
   external                   fgauss  
   character*40               fname_header
   real*8                     time1, time2, time3, time4
   logical                    flag_sparse, flag_exit_dsum
   logical                    flag_order
   integer*4                  feast_nemax_save
#ifdef PSPARSE
   integer*4                  feast_fpm_save(64)
#else
   integer*4                  feast_fpm_save(128)
#endif
   integer*4, allocatable  :: feast_ne_save(:,:)
   real*8                     feast_emin_save, feast_emax_save
#ifdef MPI                  
  integer*4                   ourjob(nprocs)
  integer*4                   ourjob_disp(0:nprocs-1)
#else                        
  integer*4                   ourjob(1)
  integer*4                   ourjob_disp(0)
#endif

#ifdef MPI
   if_main time1 = MPI_Wtime()
#else
   call cpu_time(time1)
#endif

   write(message,*)''  ; write_msg
   write(message,'(A)')'START: DOS EVALUATION'  ; write_msg

   neig    = PGEOM%neig
   ispin   = PINPT%ispin
   nspin   = PINPT%nspin
   ispinor = PINPT%ispinor
   nk1     = PINPT_DOS%dos_kgrid(1)  
   nk2     = PINPT_DOS%dos_kgrid(2)  
   nk3     = PINPT_DOS%dos_kgrid(3)  
   nkpoint = nk1 * nk2 * nk3
   nparam  = PINPT%nparam
   kshift  = PINPT_DOS%dos_kshift
   nediv   = PINPT_DOS%dos_nediv
   emax    = PINPT_DOS%dos_emax
   emin    = PINPT_DOS%dos_emin
   e_range = emin + (/(k, k=0,nediv-1)/) * (emax - emin)/dble(nediv - 1)
   sigma   = PINPT_DOS%dos_smearing
   iband   = PINPT_DOS%dos_iband
   fband   = PINPT_DOS%dos_fband 

   flag_order = PINPT%flag_get_band_order

   iadd       = 10 ; iaddk = 8  ; inck = 1

   call mpi_job_distribution_chain(nediv, ourjob, ourjob_disp)

   if(PINPT%flag_noncollinear) then ! set default fband if fband has not been pre-defined
     if(fband .eq. 999999) fband = neig * 2
   else
     if(fband .eq. 999999) fband = neig
   endif
   nband   = fband - iband + 1

   if(PINPT_DOS%dos_flag_sparse) then ! set for FEAST eigen solver if DOS_SPARSE=.TRUE.
     flag_sparse = PINPT_DOS%dos_flag_sparse

     if(PINPT%flag_sparse) then
       feast_emin_save = PINPT%feast_emin
       feast_emax_save = PINPT%feast_emax
       feast_fpm_save  = PINPT%feast_fpm
       feast_nemax_save= PINPT%feast_nemax
       if(allocated(PINPT%feast_ne)) then
         allocate(feast_ne_save(PINPT%nspin, PKPTS%nkpoint))
         feast_ne_save = PINPT%feast_ne
         deallocate(PINPT%feast_ne)
       endif
     endif

     PINPT%feast_emin  = emin
     PINPT%feast_emax  = emax
     PINPT%feast_nemax = nband

     if(PINPT%feast_nemax .gt. neig * ispinor) then
       write(message,'(A,I0,A)')'    !WARN! The NE_MAX (',PINPT%feast_nemax,') of DOS_EWINDOW tag is larger than the eigenvalues (NEIG)'  ; write_msg
       write(message,'(A,I0,A)')'           of the system (',PGEOM%neig * PINPT%ispinor,'). Hence, we enforce NEMAX = NEIG.'  ; write_msg
       write(message,'(A,I0,A)')'           Otherwise, you can reduce the expected NE_MAX within the EWINDOW with a proper guess.'  ; write_msg
       PINPT%feast_nemax = PINPT%nband
     endif
   elseif(.not. PINPT_DOS%dos_flag_sparse) then
     flag_sparse = PINPT_DOS%dos_flag_sparse
   endif

   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   call get_reci(b1,b2,b3, a1,a2,a3)
   call vcross(b2xb3,b2,b3)
   bzvol=dot_product(b1,b2xb3)
   dkv  = bzvol / dble(nkpoint)

   allocate( param(nparam) )
   allocate( param_const(5,nparam) )
   param   = PINPT%param
   param_const = PINPT%param_const
   allocate( E(nband*nspin,nkpoint) )
   if_main allocate( V(neig*ispin,nband*nspin,nkpoint) )
   allocate( myV(neig*ispin,nband*nspin) )
   allocate( kpoint(3,nkpoint) )
   allocate( kpoint_reci(3,nkpoint) )
   allocate( PINPT_DOS%dos_kpoint(3,nkpoint) )
   allocate( PINPT_DOS%dos_erange(nediv) )

   allocate( PINPT_DOS%dos_tot(nspin,nediv), dos_tot(nspin,nediv) )
   PINPT_DOS%dos_tot = 0d0
             dos_tot = 0d0
   
   if(PINPT_DOS%dos_flag_print_ldos) then
     allocate(PINPT_DOS%ldos_tot(PGEOM%max_orb,PINPT_DOS%dos_ldos_natom,nspin,nediv) )
     allocate(          ldos_tot(PGEOM%max_orb,PINPT_DOS%dos_ldos_natom,nspin,nediv) )
     PINPT_DOS%ldos_tot= 0d0
               ldos_tot= 0d0
   endif

   PINPT_DOS%dos_erange = e_range
   call get_kgrid(kpoint,kpoint_reci,nk1,nk2,nk3,kshift, PGEOM, PINPT_DOS%dos_flag_gamma)
   if(PINPT_DOS%dos_flag_print_kpoint) then 
     if_main call print_kpoint(kpoint_reci, nkpoint, PINPT_DOS%dos_kfilenm)
   endif

   call get_eig(NN_TABLE,kpoint,nkpoint,PINPT, E, V, neig, iband, nband, &
                PINPT%flag_get_orbital, flag_sparse, .true., .true.) !, flag_order)

   if(flag_sparse) then
     allocate(ne_found(PINPT%nspin, nkpoint))
     ne_found = PINPT%feast_ne
   else
     allocate(ne_found(PINPT%nspin, nkpoint))
     ne_found = PINPT%nband
   endif

#ifdef MPI
   call MPI_Barrier(mpi_comm_earth, mpierr)
#endif

   ! main routine for DOS evaluation
!kp:do ik = 1 + myid, nkpoint, nprocs
 write(message,'(A)')' ... calculating DOS ...'  ; write_msg
 if_main call time_check(time4,time3,'init')
 if(PINPT_DOS%dos_flag_print_ldos) then
   if_main call report_memory(int8(size(ldos_tot)) * nprocs * 2, 8, 'LDOS(total)   ')
 endif

kp:do ik = 1,  nkpoint
     if(PINPT_DOS%dos_flag_print_ldos) then
#ifdef MPI
       if_main myV = V(:,:,ik)
       call MPI_BCAST(myV, size(myV), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
#else
       myV = V(:,:,ik)
#endif
     endif
     if(nkpoint .lt. 10) then
       write(message,'(A,I0,A,I0)')  '         STAT KP: ', ik,'/',nkpoint  ; write_msg
     else
       if( ik/real(nkpoint)*100d0 .ge. real(iaddk*inck) ) then
         write(message,'(A,F10.3,A)')'         STAT KP: ', ik/real(nkpoint)*100d0, ' %'  ; write_msg
         inck = inck + 1
       endif
     endif
!eig:do ie = 1 + myid, nediv, nprocs
     inc = 1
 eig:do ie = sum(ourjob(1:myid)) + 1, sum(ourjob(1:myid+1))
       if(nkpoint .lt. 10) then
         if( (ie-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0 .ge. real(iadd*inc) ) then
           write(message,'(A,F10.3,A)')'            STAT EN: ', (ie-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0, ' %'  ; write_msg
           inc = inc + 1
         endif
       endif

    sp:do is = 1, nspin
    dsum:do i = 1, ne_found(is,ik) ! init_e, fina_e
!          dos_ = fgauss(sigma, e_range(ie) - E(i+nband*(is-1),ik) ) * dkv
           dos_ = fgauss(sigma, e_range(ie) - E(i+nband*(is-1),ik) ) / nkpoint
           dos_tot(is,ie) = dos_tot(is,ie) + dos_

           if(PINPT_DOS%dos_flag_print_ldos) then
             ii = i+nband*(is-1) ! get band index according to the spin index
             do ia = 1, PINPT_DOS%dos_ldos_natom
!              iatom = PINPT_DOS%dos_ldos_atom(ia) ! which atom will be resolved?
               ! which matrix index is the starting point for "iatom"
               imatrix = sum(PGEOM%n_orbital(1:PINPT_DOS%dos_ldos_atom(ia))) - &
                             PGEOM%n_orbital(PINPT_DOS%dos_ldos_atom(ia)) + 1 
               do im = 1, PGEOM%n_orbital(PINPT_DOS%dos_ldos_atom(ia))
                 mm = im+imatrix-1 + PGEOM%neig*(is-1) ; rho =       real(conjg(myV(mm     ,ii))*myV(mm     ,ii))
                 if(ispinor .eq. 2)                      rho = rho + real(conjg(myV(mm+neig,ii))*myV(mm+neig,ii))
                 ldos_tot(im,ia,is,ie) = ldos_tot(im,ia,is,ie) + rho * dos_
               enddo
             enddo
           endif

         enddo dsum
       enddo sp
     enddo eig
   enddo kp
#ifdef MPI
   call MPI_REDUCE(dos_tot, PINPT_DOS%dos_tot, nediv*nspin, MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
   if_main call print_dos(PINPT_DOS, PINPT)
   if(PINPT_DOS%dos_flag_print_ldos) then
     call MPI_ALLREDUCE(ldos_tot, PINPT_DOS%ldos_tot, size(ldos_tot), MPI_REAL8, MPI_SUM,  mpi_comm_earth, mpierr)
     call print_ldos(PINPT_DOS, PINPT, PGEOM)
   endif
#else
   PINPT_DOS%dos_tot = dos_tot
   call print_dos(PINPT_DOS, PINPT)
   if(PINPT_DOS%dos_flag_print_ldos) then
     PINPT_DOS%ldos_tot = ldos_tot
     call print_ldos(PINPT_DOS, PINPT, PGEOM)
   endif
#endif
   if_main call time_check(time4,time3)
   write(message,'(A,F10.4,A)')' ... calculating DOS ... DONE  : ',time4, ' (sec)'  ; write_msg


   ! NOTE: if flag_sparse = .true. dos_flag_print_eigen will not be activated due to the eigenvalue 
   !       ordering is not well defined in this case.
   if(PINPT_DOS%dos_flag_print_eigen .and. myid .eq. 0 .and. .not. flag_sparse) then
     allocate(E_(nspin,nkpoint))
     allocate(V_(neig*ispin,nspin,nkpoint))

     do i = 1, PINPT_DOS%dos_n_ensurf
       call get_ensurf_fname_header(PINPT_DOS%dos_ensurf(i), fname_header)
       do is = 1, nspin
         E_(is,:)   = E(  PINPT_DOS%dos_ensurf(i) - iband + 1 - (is-1)*(neig-PINPT_DOS%dos_n_ensurf), :)
         V_(:,is,:) = V(:,PINPT_DOS%dos_ensurf(i) - iband + 1 - (is-1)*(neig-PINPT_DOS%dos_n_ensurf), :)
       enddo
      call print_energy_ensurf(kpoint, nkpoint,PINPT_DOS%dos_ensurf(i), nspin, E_, V_, PGEOM, PINPT, fname_header, PINPT_DOS%dos_kunit)
     enddo
   endif

   deallocate( param )
   deallocate( param_const)
   deallocate( E )
   if(allocated(V )) deallocate( V )
   if(allocated(E_)) deallocate( E_)
   if(allocated(V_)) deallocate( V_)
   if(allocated(myV))deallocate(myV)
   deallocate( kpoint )
   deallocate( kpoint_reci )
   deallocate( ne_found )
   deallocate( dos_tot )
   if(allocated(ldos_tot)) deallocate(ldos_tot)

   if(PINPT_DOS%dos_flag_sparse) then ! set for FEAST eigen solver if DOS_SPARSE=.TRUE.
     if(PINPT%flag_sparse) then
       PINPT%feast_emin = feast_emin_save
       PINPT%feast_emax = feast_emax_save
       PINPT%feast_fpm  = feast_fpm_save 
       PINPT%feast_nemax= feast_nemax_save
       if(allocated(PINPT%feast_ne)) then
         deallocate(PINPT%feast_ne)
         allocate(PINPT%feast_ne(PINPT%nspin, PKPTS%nkpoint))
         PINPT%feast_ne = feast_ne_save
         deallocate(feast_ne_save)
       endif
     endif
   endif

#ifdef MPI
   if_main time2 = MPI_Wtime()
#else
   call cpu_time(time2)
#endif
   write(message,'(A,F12.3)')'END: DOS EVALUATION. TIME ELAPSED (s) =',time2-time1  ; write_msg

return
endsubroutine

function fgauss(sigma, x)
   use parameters, only : pi, pi2
   implicit none
   real*8   sigma,sigma2
   real*8   x,xx
   real*8   fgauss
   xx = x**2
   sigma2 = sigma**2

   fgauss= exp(-0.5d0*xx/sigma2)/(sigma*sqrt(pi2))

return
end function
subroutine get_ensurf_fname_format(i, format_string)
   implicit none
   integer*4    i
   character*40 format_string

   if(i .gt. 0 .and. i .lt. 10) then
     write(format_string, '( "(A,",A,",I1,A)" )')'"0000"'
   elseif(i .ge. 10 .and. i .lt. 100) then
     write(format_string, '( "(A,",A,",I2,A)" )')'"000"'
   elseif(i .ge. 100 .and. i .lt. 1000) then
     write(format_string, '( "(A,",A,",I3,A)" )')'"00"'
   elseif(i .ge. 1000 .and. i .lt. 10000) then
     write(format_string, '( "(A,",A,",I4,A)" )')'"0"'
   elseif(i .ge. 10000 .and. i .lt. 100000) then
     write(format_string, '( "(A,      I5,A)" )')
   endif

return
endsubroutine
subroutine get_ensurf_fname_header(i, fname_header)
   implicit none
   integer*4    i
   character*40 format_string
   character*40 fname_header

   call get_ensurf_fname_format(i, format_string)

   write(fname_header, format_string)'ENSURF.EIG.',i

return
endsubroutine

subroutine print_dos(PINPT_DOS, PINPT)
  use parameters, only : dos, pid_dos, incar
  implicit none
  type(dos)    :: PINPT_DOS
  type(incar)  :: PINPT
  integer*4       i, ie
  real*8          dos_data(PINPT_DOS%dos_nediv,PINPT%nspin)
  real*8          e_range(PINPT_DOS%dos_nediv)
  character*40    filenm
  logical         flag_collinear

  filenm = PINPT_DOS%dos_filenm
  e_range = PINPT_DOS%dos_erange

  open(pid_dos, file=trim(filenm), status = 'unknown')

  write(pid_dos,'(A,I8,A,F16.8)')'# NDIV = ',PINPT_DOS%dos_nediv,' dE=',e_range(2)-e_range(1)
  if(.not.PINPT%flag_collinear) then
    write(pid_dos,'(A)')'#    energy (ev)           dos '
  elseif(PINPT%flag_collinear) then
    write(pid_dos,'(A)')'#    energy (ev)          dos-up           dos-dn'
  endif
  do ie = 1, PINPT_DOS%dos_nediv
    if(.not.PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,1x,F16.8)')e_range(ie), PINPT_DOS%dos_tot(1,ie)
    elseif(PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,1x,F16.8)')e_range(ie), PINPT_DOS%dos_tot(1:2,ie)
    endif
  enddo

  close(pid_dos)

return
endsubroutine

subroutine print_ldos(PINPT_DOS, PINPT, PGEOM)
   use parameters, only: dos, pid_ldos, incar, poscar
   use mpi_setup
   implicit none
   type(dos)    :: PINPT_DOS
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   integer*4       im, ia, iaa, ie, iatom
   integer*4       my_pid_ldos
   real*8          e_range(PINPT_DOS%dos_nediv)
   character*40    filenm
#ifdef MPI
   integer*4       ourjob(nprocs)
   integer*4       ourjob_disp(0:nprocs-1)
   call mpi_job_distribution_chain(PINPT_DOS%dos_ldos_natom, ourjob, ourjob_disp)   
#else
   integer*4  ourjob(1)
   integer*4  ourjob_disp(0)
   call mpi_job_distribution_chain(PINPT_DOS%dos_ldos_natom, ourjob, ourjob_disp)
#endif

   e_range = PINPT_DOS%dos_erange

   do ia = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
     iatom = PINPT_DOS%dos_ldos_atom(ia)
     my_pid_ldos = pid_ldos + myid
     write(filenm,'(A,A,I0,A)') trim(PINPT_DOS%ldos_filenm),'_atom.',iatom,'.dat'
     open(my_pid_ldos, file=trim(filenm), status = 'unknown')

     write(my_pid_ldos,'(A,I0,A,A,A )')'# ATOM = ',iatom,' (spec = ',trim(PGEOM%c_spec(PGEOM%spec(iatom))),' )'
     write(my_pid_ldos,'(A,I8,A,F16.8)')'# NDIV = ',PINPT_DOS%dos_nediv,' dE=',e_range(2)-e_range(1)
   
     if(.not. PINPT%flag_collinear) then
       write(my_pid_ldos,'(A)',ADVANCE='NO')          '#    energy (ev)        dos_total'
       do im=1, PGEOM%n_orbital(ia)-1
         write(my_pid_ldos,'(A15,I1,1x)',ADVANCE='NO')'      orb-',im
       enddo
       write(my_pid_ldos,'(A15,I1,1x)',ADVANCE='YES')  '      orb-',im
     elseif(PINPT%flag_collinear) then
       write(my_pid_ldos,'(A)',ADVANCE='NO')          '#    energy (ev)     dos_total-up     dos_total-dn'
       do im=1, PGEOM%n_orbital(ia)
         write(my_pid_ldos,'(A15,I1,1x)',ADVANCE='NO')' up-orb-',im
       enddo
       do im=1, PGEOM%n_orbital(ia)-1
         write(my_pid_ldos,'(A15,I1,1x)',ADVANCE='NO')' dn-orb-',im
       enddo
       write(my_pid_ldos,'(A15,I1,1x)',ADVANCE='YES')  ' dn-orb-',im
     endif  

     do ie = 1, PINPT_DOS%dos_nediv
       if(.not.PINPT%flag_collinear) then
         write(my_pid_ldos,'(F16.8,1x,  F16.8    , *(F16.8,1x))')e_range(ie), sum(PINPT_DOS%ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie)), & 
                                                                                  PINPT_DOS%ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie)
       elseif(PINPT%flag_collinear) then                                                                                          
         write(my_pid_ldos,'(F16.8,1x,2(F16.8,1x), *(F16.8,1x))')e_range(ie), sum(PINPT_DOS%ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie)), &
                                                                              sum(PINPT_DOS%ldos_tot(1:PGEOM%n_orbital(iatom),ia,2,ie)), &
                                                                                  PINPT_DOS%ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie) , &
                                                                                  PINPT_DOS%ldos_tot(1:PGEOM%n_orbital(iatom),ia,2,ie)   
       endif
     enddo
     
     close(my_pid_ldos)
   
   enddo

   return
endsubroutine
