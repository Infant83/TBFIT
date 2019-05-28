#include "alias.inc"
subroutine replot_dos_band(PINPT, PGEOM, PKPTS, PRPLT)
   use parameters, only : incar, poscar, replot, kpoints, pi, pi2
   use mpi_setup
   use time
   use memory
   implicit none
   type(incar  )            :: PINPT
   type(poscar )            :: PGEOM
   type(replot )            :: PRPLT
   type(kpoints)            :: PKPTS
   integer*4                   nbasis, nband, nspin
   integer*4                   ispinor, ispin, nkp
   integer*4                   nediv
   integer*4                   k
   integer*4                   ne_found(PINPT%nspin, PKPTS%nkpoint)
   integer*4                   init_erange, fina_erange
   real*8                      emin, emax, e_range(PRPLT%replot_dos_nediv)
   real*8                      emin_band, emax_band
   real*8                      sigma
   integer*4                   mpierr
   character*4                 kmode
   logical                     flag_vector
   logical                     flag_replot_dos, flag_replot_ldos, flag_replot_sldos
   logical                     flag_replot_proj_band, flag_replot_band
   logical                     flag_erange
   logical                     flag_wf
   logical                     flag_formatted
   logical                     flag_exit
   real*8,     allocatable  :: E(:,:) ! (nband * nspin, nkpoint)
   ! be careful that PINPT%nband should be same as nband which will be found by check_band_header routine...
   ! this routine should be modified in this sense...  KHJ May. 15 2019.
   complex*16, allocatable  :: V(:,:,:), V2(:,:,:)
!  complex*16                  V(PGEOM%neig*PINPT%ispin,PINPT%nband*PINPT%nspin,PKPTS%nkpoint) ! wavevector
!  real*8                      V2(PGEOM%neig,PINPT%nband*PINPT%nspin,PKPTS%nkpoint) ! V2(m,n,k) = <phi_m|, psi_nk>, m = orbital, n = band, k = kpoint index
   real*8                      time2, time1
   character*80                cjob
   character*2                 c_mode,c_mode_print
   logical                     flag_sparse

   flag_exit = .false.

   if_main call time_check(time2,time1,'init')

   flag_replot_dos       = PRPLT%flag_replot_dos
   flag_replot_ldos      = PRPLT%flag_replot_ldos
   flag_replot_sldos     = PRPLT%flag_replot_sldos
   flag_replot_proj_band = PRPLT%flag_replot_proj_band
   flag_replot_band      = PRPLT%flag_replot_band
   flag_formatted        = PRPLT%flag_replot_formatted
   c_mode_print          = PRPLT%replot_axis_print_mag

   if(PKPTS%flag_klinemode) kmode = 'line'
   if(PKPTS%flag_kgridmode) kmode = 'grid'

   if_main write(6,*)''
   if_main write(6,*)''
   if_main write(6,'(A)')' ** Program run in REPLOT mode '
   if_main write(6,*)''
   cjob = ' '
   if(flag_replot_dos)      write(cjob,'(A,A)') trim(cjob), ' + DOS'
   if(flag_replot_ldos)     write(cjob,'(A,A)') trim(cjob), ' + LDOS'
   if(flag_replot_sldos)    write(cjob,'(A,A)') trim(cjob), ' + SLDOS'
   if(flag_replot_proj_band)write(cjob,'(A,A)') trim(cjob), ' + PROJ_BAND'
   if(flag_replot_band) then
     if(c_mode_print .eq. 'mx' .or. c_mode_print .eq. 'my' .or. c_mode_print .eq. 'mz') then
                            write(cjob,'(A,4A)')trim(cjob), ' + BAND (<sigma_i>, i=',c_mode_print,')'
     elseif(c_mode_print .eq. 'wf') then
                            write(cjob,'(A,2A)')trim(cjob), ' + BAND (wavefunction)'
     elseif(c_mode_print .eq. 'rh') then
                            write(cjob,'(A,2A)')trim(cjob), ' + BAND (<phi_i|psi_nk>, i=orbital)'
     elseif(c_mode_print .eq. 'no') then
                            write(cjob,'(A,2A)')trim(cjob), ' + BAND '
     endif
   endif

   if_main write(6,'(A,A,A)')'  START REPLOT: ',trim(cjob),' EVALUATION'
   if_main write(6,*)''

   nspin                 = PINPT%nspin
   nbasis                = PGEOM%neig
   call check_band_header(flag_vector, flag_wf, flag_sparse, nband, emin_band, emax_band, nspin, &
                          flag_erange, init_erange, fina_erange, c_mode, flag_formatted)

   if(.not.flag_sparse) nband = PINPT%nband
   ispinor               = PINPT%ispinor
   ispin                 = PINPT%ispin
   nkp                   = PKPTS%nkpoint
   nediv                 = PRPLT%replot_dos_nediv
   emax                  = PRPLT%replot_dos_emax
   emin                  = PRPLT%replot_dos_emin
   e_range               = emin + (/(k, k=0,nediv-1)/) * (emax - emin)/dble(nediv - 1)
   sigma                 = PRPLT%replot_dos_smearing
   
   allocate(E(nband*nspin,nkp))
   if(flag_vector) then
     allocate(V2(nbasis,nband*nspin,nkp))
     if(flag_wf) then
       allocate(V(nbasis*ispin,nband*nspin,nkp))
     endif
   endif

   if(nband .ne. PINPT%nband) then
     write(6,'(A)')'    ! ERROR: nband read from band structure file is differ from the setup PINPT%nband. This discrepancy may'
     write(6,'(A)')'           : be originated from NE_MAX tag of your band structure that has been calculated under EWINDOW'
     write(6,'(A)')'           : tag which uses sparse matrix solver. Exit program...'
     kill_job
   endif

   ! reading energy + wave vector (if requested)
   if(flag_wf) then
     if(flag_vector) then
       if_main call report_memory(int8(size(V)), 16, 'Wave vector   ') 
     endif
     if_main call read_energy_tbfit_V(E, V, ne_found, ispin, nspin, nbasis, nband, nkp, kmode, &
                                      flag_vector, flag_wf, flag_formatted, flag_exit)

#ifdef MPI
     call MPI_BCAST(flag_exit, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
#endif
     if(flag_exit) then
       kill_job
     endif

   elseif(.not. flag_wf) then
     if(.not. flag_formatted) then
       if_main write(6,'(A)')'    ! WARN: reading with C_MODE:"rh" is not supported in the current version.'
       if_main write(6,'(A)')'      Exit program... '
       kill_job
     endif
     if(flag_vector) then
       if_main call report_memory(int8(size(V2)), 8, 'Wave vector   ') 
     endif
     if_main call read_energy_tbfit_V2(E, V2, ne_found, ispin, nspin, nbasis, nband, nkp, kmode, &
                                       flag_vector, flag_wf)
   endif

   ! convert wave vector into its squared form 
   if(flag_vector .and. flag_wf .and. myid .eq. 0) then
     if_main call V_to_V2(V, V2, ispin, nspin, nbasis, nband, nkp)
   endif

   allocate( PRPLT%replot_dos_erange(nediv) )
   PRPLT%replot_dos_erange = e_range

#ifdef MPI
   call MPI_BCAST(E, size(E), MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

   ! replot dos & ldos
   if(flag_replot_dos .or. flag_replot_ldos .or. flag_replot_sldos) then
     call get_dos_ldos(PINPT, PGEOM, PRPLT, E, V2, nspin, nbasis, nband, nkp, e_range, nediv, sigma, &
                       flag_replot_dos, flag_replot_ldos, flag_replot_sldos)
   endif   

   ! replot pband
   if(flag_replot_proj_band) then
     if_main call print_replot_energy_proj(PKPTS, E, V2, PGEOM, PRPLT, ne_found, emin_band, emax_band, &
                                           ispinor, nband, nspin, ispin, nbasis, &
                                           flag_erange, flag_vector, flag_sparse, &
                                           init_erange, fina_erange, c_mode)
   endif

   ! replot band
   if(flag_replot_band) then
     if_main call print_replot_energy(PKPTS, E, V, PGEOM, PRPLT, ne_found, emin_band, emax_band, &
                                           ispinor, nband, nspin, ispin, nbasis, &
                                           flag_erange, flag_vector, flag_sparse, &
                                           init_erange, fina_erange, c_mode_print)
   endif


   if_main call time_check(time2,time1)
   if_main write(6,*)''
   if_main write(6,'(3A,F10.4,A)')'  END REPLOT: ',trim(cjob),' EVALUATION : ',time2, ' (sec)'
   if_main write(6,*)''


!  if(allocated(V)) deallocate(V)
!  if(allocated(V2)) deallocate(V2)

   return
endsubroutine

subroutine get_dos_ldos(PINPT, PGEOM, PRPLT, E, V2, nspin, nbasis, nband, nkp, e_range, nediv, sigma, &
                        flag_replot_dos, flag_replot_ldos, flag_replot_sldos)
   use parameters, only : incar, poscar, replot
   use mpi_setup
   use memory
   implicit none
   type(incar  )            :: PINPT
   type(poscar )            :: PGEOM
   type(replot )            :: PRPLT
   integer*4                   ik, ie, ia, is, i
   integer*4                   iadd, inc, iaddk, inck
   integer*4                   iatom, ii, fi
   integer*4                   nbasis, nband, nspin, nkp
   integer*4                   nediv
   integer*4                   mpierr
   real*8                      e_range(nediv), de, sigma
   real*8                      E(nband*nspin,nkp)
   real*8                      V2(nbasis,nband*nspin,nkp) ! not wavevector, 
   real*8                      myV(nbasis,nband*nspin)
   real*8                      replot_dos_tot(nspin,nediv)
   real*8                      replot_dos_ntot(nspin,0:nediv)
   real*8, allocatable      :: replot_ldos_tot(:,:,:,:)
   real*8, allocatable      :: replot_sldos_sum(:,:,:)
   logical                     flag_replot_dos, flag_replot_ldos, flag_replot_sldos
   real*8                      dos_
   real*8                      fgauss
   external                    fgauss
   real*8                      a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
   real*8                      b2xb3(3),bzvol,dkv
   integer*4                   ldos_natom, ldos_atom(PRPLT%replot_ldos_natom)
   integer*4                   n_orbital(PGEOM%n_atom)
   real*8                      coord_cart(3,PGEOM%n_atom)
#ifdef MPI                  
   integer*4                   ourjob(nprocs)
   integer*4                   ourjob_disp(0:nprocs-1)
#else                        
   integer*4                   ourjob(1)
   integer*4                   ourjob_disp(0)
#endif

   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   call get_reci(b1,b2,b3, a1,a2,a3)
   call vcross(b2xb3,b2,b3)
   bzvol=dot_product(b1,b2xb3)
   dkv  = bzvol / dble(nkp)
   de   = e_range(2)-e_range(1)

   ldos_natom = PRPLT%replot_ldos_natom
   ldos_atom  = PRPLT%replot_ldos_atom
   n_orbital  = PGEOM%n_orbital
   iadd       = 10 ; iaddk = 8  ; inck = 1

   allocate(PRPLT%replot_dos_tot(nspin, nediv))
   PRPLT%replot_dos_tot  = 0d0
         replot_dos_tot  = 0d0
   allocate(PRPLT%replot_dos_ntot(nspin,0:nediv))
   PRPLT%replot_dos_ntot = 0d0
         replot_dos_ntot = 0d0 
   if_main write(6,*)''
   if_main call report_memory(int8(size(replot_dos_tot))*2*nprocs, 8, 'DOS(total)    ')

   if(flag_replot_ldos) then
     allocate(PRPLT%replot_ldos_tot(PGEOM%max_orb, PRPLT%replot_ldos_natom, nspin, nediv))
     allocate(      replot_ldos_tot(PGEOM%max_orb, PRPLT%replot_ldos_natom, nspin, nediv))
     PRPLT%replot_ldos_tot = 0d0
           replot_ldos_tot = 0d0
     if_main call report_memory(int8(size(replot_ldos_tot))*2*nprocs, 8, 'LDOS(total)   ')
   endif

   if(flag_replot_sldos) then
     allocate(PRPLT%replot_sldos_sum(PGEOM%n_atom, nspin, nediv))
     allocate(      replot_sldos_sum(PGEOM%n_atom, nspin, nediv))
     PRPLT%replot_sldos_sum = 0d0
           replot_sldos_sum = 0d0
     if_main call report_memory(int8(size(replot_sldos_sum))*2*nprocs, 8, 'SLDOS(total)  ')
   endif

   if(flag_replot_ldos .or. flag_replot_sldos) then
     if_main call report_memory(int8(size(myV))*nprocs + size(V2), 8, 'Wave vector2  ')
   endif

   call mpi_job_distribution_chain(nediv, ourjob, ourjob_disp)
   if_main write(6,*)''
   do ik = 1, nkp
     if(nkp .lt. 10) then
       if_main write(6,'(A,I0,A,I0)')  '     STAT KP: ', ik,'/',nkp
     else
       if( ik/real(nkp)*100d0 .ge. real(iaddk*inck) ) then
         if_main write(6,'(A,F10.3,A)')'     STAT KP: ', ik/real(nkp)*100d0, ' %'
         inck = inck + 1
       endif
     endif

     if(flag_replot_ldos .or. flag_replot_sldos) then
#ifdef MPI
       myV = V2(:,:,ik)
       call MPI_BCAST(myV, size(myV), MPI_REAL8, 0, mpi_comm_earth, mpierr)
#else
       myV = V2(:,:,ik)
#endif
     endif

     inc = 1
     do ie = sum(ourjob(1:myid)) + 1, sum(ourjob(1:myid+1))
       if(nkp .lt. 10) then
         if( (ie-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0 .ge. real(iadd*inc) ) then
           if_main write(6,'(A,F10.3,A)')'    STAT EN: ', (ie-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0, ' %'
           inc = inc + 1
         endif
       endif

       do is = 1, nspin
         do i = 1, nband

           dos_ = fgauss(sigma, e_range(ie) - E(i+nband*(is-1),ik)) / nkp
           replot_dos_tot(is, ie) = replot_dos_tot(is, ie) + dos_

           if(flag_replot_ldos) then
             do ia = 1, ldos_natom
               iatom = ldos_atom(ia)
               ii = sum(n_orbital(1:iatom)) - n_orbital(iatom) + 1
               fi = ii + n_orbital(iatom)-1
               replot_ldos_tot(1:n_orbital(iatom),ia,is,ie) = replot_ldos_tot(1:n_orbital(iatom),ia,is,ie) + dos_ * myV(ii:fi,i+nband*(is-1))
             enddo
           endif

           if(flag_replot_sldos) then
             do iatom = 1, PGEOM%n_atom
               ii = sum(n_orbital(1:iatom)) - n_orbital(iatom) + 1
               fi = ii + n_orbital(iatom)-1
               replot_sldos_sum(iatom,is,ie) = replot_sldos_sum(iatom,is,ie) + dos_ * sum(myV(ii:fi,i+nband*(is-1)))
             enddo
           endif


         enddo
       enddo
     enddo
   enddo

#ifdef MPI
   call MPI_REDUCE(replot_dos_tot, PRPLT%replot_dos_tot, size(replot_dos_tot), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
   do ie = 1, nediv
     do is = 1, PINPT%nspin
       PRPLT%replot_dos_ntot(is,ie) = PRPLT%replot_dos_ntot(is,ie-1) + PRPLT%replot_dos_tot(is, ie) * de
     enddo
   enddo
   if(flag_replot_dos) then
     if_main call print_replot_dos(PRPLT, PINPT)   
     if_main write(6,'(A)')' Done!'
   endif
   if(flag_replot_ldos) then
     call MPI_ALLREDUCE(replot_ldos_tot, PRPLT%replot_ldos_tot, size(replot_ldos_tot), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     call print_replot_ldos(PRPLT, PINPT, PGEOM)
     if_main write(6,'(A)')' Done!'
   endif
   if(flag_replot_sldos) then
     call MPI_REDUCE(replot_sldos_sum, PRPLT%replot_sldos_sum, size(replot_sldos_sum), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
     if_main call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart)
     call MPI_BCAST(coord_cart, size(coord_cart), MPI_REAL8, 0, mpi_comm_earth, mpierr)
     call print_bond(PRPLT, PINPT, PGEOM, coord_cart)
     if_main write(6,'(A)')' Done!'
   endif

#else
   PRPLT%replot_dos_tot = replot_dos_tot
   do ie = 1, nediv
     do is = 1, PINPT%nspin
       PRPLT%replot_dos_ntot(is,ie) = PRPLT%replot_dos_ntot(is,ie-1) + PRPLT%replot_dos_tot(is, ie) * de
     enddo
   enddo
   if(flag_replot_dos) then
     call print_replot_dos(PRPLT, PINPT)
     write(6,'(A)')'    write DOS ... Done!'
   endif
   if(flag_replot_ldos) then
     PRPLT%replot_ldos_tot = replot_ldos_tot
     call print_replot_ldos(PRPLT, PINPT, PGEOM)
     write(6,'(A)')'    write LDOS ... Done!'
   endif
   if(flag_replot_sldos) then
     PRPLT%replot_sldos_sum = replot_sldos_sum
     call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart)
     write(6,'(A)')'    write SLDOS ... Done!'
   endif

#endif

   if(allocated(replot_ldos_tot))  deallocate(replot_ldos_tot)
   if(allocated(replot_sldos_sum)) deallocate(replot_sldos_sum)
   return
endsubroutine  

subroutine print_replot_dos(PRPLT, PINPT)
  use parameters, only : replot, pid_dos, incar
  implicit none
  type(replot) :: PRPLT 
  type(incar)  :: PINPT
  integer*4       i, ie
  real*8          e_range(PRPLT%replot_dos_nediv)
  character*40    filenm
  logical         flag_collinear

  filenm = 'DOS.replot.dat'
  e_range = PRPLT%replot_dos_erange

  open(pid_dos, file=trim(filenm), status = 'unknown')
  write(6,'(A)')' '
  write(6,'(3A)',ADVANCE='no')'   WRITING.... density of states (DOS): ',trim(filenm),' ... '

  write(pid_dos,'(A,I8,A,F16.8)')'# NDIV = ',PRPLT%replot_dos_nediv,' dE=',e_range(2)-e_range(1)
  if(.not.PINPT%flag_collinear) then
    write(pid_dos,'(A)')'#    energy (ev)           dos             nelect'
  elseif(PINPT%flag_collinear) then
    write(pid_dos,'(A)')'#    energy (ev)          dos-up           dos-dn           nelect-up        nelect-dn'
  endif
  do ie = 1, PRPLT%replot_dos_nediv
    if(.not.PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,*(1x,F16.8))')e_range(ie), PRPLT%replot_dos_tot(1,ie), PRPLT%replot_dos_ntot(1,ie)
    elseif(PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,*(1x,F16.8))')e_range(ie), PRPLT%replot_dos_tot(1:2,ie), PRPLT%replot_dos_ntot(1:2,ie)
    endif
  enddo

  close(pid_dos)

return
endsubroutine
subroutine print_replot_energy(PKPTS,E,V,PGEOM,PRPLT, ne_found, emin, emax, ispinor, nband, nspin, ispin, nbasis, &
                                    flag_erange, flag_vector, flag_sparse, init_erange, fina_erange, c_mode)
   use parameters, only: pid_energy, incar, poscar, kpoints, replot, zi
   use mpi_setup
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   type(kpoints):: PKPTS
   integer*4       ie, is, ik, im
   integer*4       nbasis, ispinor, ispin
   integer*4       nband, nspin
   integer*4       nkpoint
   integer*4       init_erange, fina_erange
   integer*4       init_e, fina_e
   integer*4       ne_found(nspin, PKPTS%nkpoint)
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   logical         flag_erange, flag_vector
   logical         flag_sparse
   logical         flag_collinear, flag_noncollinear
   real*8          E(nband*nspin,PKPTS%nkpoint)
   complex*16      V(nbasis*ispin,nband*nspin,PKPTS%nkpoint) ! complex wavefunction (total)
   complex*16      c_up, c_dn
   real*8          emin, emax
   character*80    fname_header
   character*80    fname
   character*6     kunit_
   character*28    kmode
   character*8     sigma
   character*2     c_mode

   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = PGEOM%neig
   sigma='sigma_0 '

   flag_collinear = .false.
   flag_noncollinear = .false.
   if(ispin .eq. 2 .and. nspin .eq. 2) then
     flag_collinear = .true.
   elseif(ispin .eq. 2 .and. nspin .eq. 1) then
     flag_noncollinear = .true.
   endif

   if(flag_vector .and. c_mode .eq. 'no') then
     flag_print_orbital = .false.
   else
     flag_print_orbital = flag_vector
   endif

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   call get_e_range(init_e, fina_e, nbasis, .false., ispinor, flag_erange, init_erange, fina_erange)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

spin:do is = 1, nspin
       fname_header = 'band_structure_TBA.replot'
       call get_fname(fname_header, fname, is, flag_collinear, flag_noncollinear)
       open(pid_energy, file=trim(fname), status = 'unknown')
       if(trim(c_mode) .eq. 'rh') then
         write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(c_mode), &
               ' ] -> <phi_ij|psi_nk> ; i,j => orbital j in atom i;  n,k => band index (n) and kpoint index (k)'
       elseif(trim(c_mode) .eq. 'wf') then
         write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(c_mode), ' ] -> wavefunction coefficients '
       elseif(trim(c_mode) .eq. 'no') then
         write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(c_mode), ' ]'
       elseif(trim(c_mode) .eq. 'mx' .or. &
              trim(c_mode) .eq. 'my' .or. &
              trim(c_mode) .eq. 'mz') then
         write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(c_mode), ' ] -> magnetization <sigma_i>'
       endif
       if(flag_sparse) then
         write(pid_energy, '(A,2(F10.4,A),I0)')'# The EWINDOW mode: energy window [EMIN:EMAX]=[ ',emin, &
                           ' : ',emax,' ], NE_MAX= ',nband
         do ik = 1, nkpoint
           write(pid_energy, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
         enddo
       elseif(flag_erange) then
          write(pid_energy, '(A,I0,A,I0,A)')'#   ERANGE=[ ',init_erange,' : ',fina_erange,' ]'
       endif

   eig:do ie =1, nband !init_e, fina_e
         write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', init_e + ie - 1,' -th eigen'
         if(.not. flag_print_orbital) then
           write(pid_energy,'(A)',ADVANCE='NO')''
         elseif(  flag_print_orbital) then
           if(c_mode .eq. 'mz') sigma='sigma_z '
           if(c_mode .eq. 'mx') sigma='sigma_x '
           if(c_mode .eq. 'my') sigma='sigma_y '
           if(c_mode .ne. 'wf') then
             write(pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
             write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
           elseif(c_mode .eq. 'wf') then
             write(pid_energy, '(1A)',ADVANCE='YES') '# wavefunction coeff.:          |ci>               '
             write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
           endif
        mm:do im=1,nbasis
             if(c_mode .ne. 'wf') then
               write(pid_energy, '(I9)',ADVANCE='NO')im
             elseif(c_mode .eq. 'wf') then
               if(ispinor .eq. 2) then
                 write(pid_energy, '(I38)',ADVANCE='NO')im
               elseif(ispinor .eq. 1) then
                 write(pid_energy, '(I19)',ADVANCE='NO')im
               endif
             endif
             if(im .ge. 30 .and. im .lt. nbasis) then
               write(pid_energy, '(A)',ADVANCE='NO')' ... '
               exit mm
             endif
           enddo mm
           write(pid_energy,'(A)')''
         endif

      kp:do ik = 1, nkpoint
           if(flag_klinemode) then
             if( ie .le. ne_found(is, ik) ) then
               write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+nband*(is-1),ik)
             elseif( ie .gt. ne_found(is, ik)) then
               write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
             endif
           elseif(flag_kgridmode) then
             if( ie .le. ne_found(is, ik) ) then
               write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+nband*(is-1),ik)
             elseif(ie .gt. ne_found(is, ik)) then
               write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
             endif
           endif

           if( ie .le. ne_found(is, ik) ) then
             if(flag_print_orbital) then
         basis:do im=1,nbasis-1
                 if(ispinor .eq. 2) then
                   c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)

                   if    (c_mode .eq. 'mz') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
                   elseif(c_mode .eq. 'mx') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
                   elseif(c_mode .eq. 'my') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
                   elseif(c_mode .eq. 'wf') then
                     write(pid_energy,'(2(F9.4,F9.4," "))',ADVANCE='NO') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                   else
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
                   endif
                 elseif(ispinor .eq. 1) then
                   c_up = V(im+nbasis*(is-1),ie+PINPT%nband*(is-1),ik)
                   if(c_mode .eq. 'wf') then
                     write(pid_energy,'(1(F9.4,F9.4," "))',ADVANCE='NO') c_up
                   else
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*c_up)
                   endif
                 endif
               enddo basis
               if(ispinor .eq. 2) then
                 c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
                 if    (c_mode .eq. 'mz') then
                   write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
                 elseif(c_mode .eq. 'mx') then
                   write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
                 elseif(c_mode .eq. 'my') then
                   write(pid_energy,'(*(F9.4))',ADVANCE='YES') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
                 elseif(c_mode .eq. 'wf') then
                   write(pid_energy,'(2(F9.4,F9.4," "))',ADVANCE='YES') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                 else
                   write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
                 endif
               elseif(ispinor .eq. 1) then
                 c_up = V(nbasis+nbasis*(is-1),ie+nband*(is-1),ik)
                 if(c_mode .eq. 'wf') then
                   write(pid_energy,'(1(F9.4,F9.4," "))',ADVANCE='YES') c_up
                 else
                   write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*c_up)
                 endif
               endif
             endif
             if(.not.flag_print_orbital) write(pid_energy,*)''
           elseif(ie .gt. ne_found(is, ik)) then
             write(pid_energy,*)''
           endif
         enddo kp
         write(pid_energy,*)''
         write(pid_energy,*)''
       enddo eig
       close(pid_energy)
     enddo spin


return
endsubroutine

subroutine print_replot_energy_proj(PKPTS,E,V,PGEOM,PRPLT, ne_found, emin, emax, ispinor, nband, nspin, ispin, nbasis, &
                                    flag_erange, flag_vector, flag_sparse, init_erange, fina_erange, c_mode)
   use parameters, only: pid_energy, incar, poscar, kpoints, replot, zi
   use mpi_setup
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   type(kpoints):: PKPTS
   integer*4       ie, is, ik, im, ia
   integer*4       isum, iatom
   integer*4       nbasis, ispinor, ispin
   integer*4       nband, nspin
   integer*4       nkpoint
   integer*4       proj_natom, proj_atom(maxval(PRPLT%replot_proj_natom(1:PRPLT%replot_nproj_sum)))
   integer*4       init_erange, fina_erange
   integer*4       init_e, fina_e
   integer*4       ne_found(nspin, PKPTS%nkpoint)
   integer*4       imatrix
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   logical         flag_proj_sum
   logical         flag_erange, flag_vector
   logical         flag_sparse
   logical         flag_collinear, flag_noncollinear
   real*8          E(nband*nspin,PKPTS%nkpoint)
   real*8          V(nbasis,nband*nspin,PKPTS%nkpoint) ! rho_nk = <psi_nk|psi_nk> is stored, not complex wavefunction
   real*8          c_up, c_dn, c_tot
   real*8          c_sum(nband,PKPTS%nkpoint)
   real*8          emin, emax
   character*80    fname_header, fname_header_sum
   character*80    fname, fname_sum
   character*6     kunit_
   character*28    kmode
   character*8     sigma
   character*2     c_mode

   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = flag_vector
   flag_proj_sum = PRPLT%flag_replot_proj_band
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = PGEOM%neig
   sigma='sigma_0 '

   flag_collinear = .false.
   flag_noncollinear = .false.
   if(ispin .eq. 2 .and. nspin .eq. 2) then
     flag_collinear = .true.
   elseif(ispin .eq. 2 .and. nspin .eq. 1) then
     flag_noncollinear = .true.
   endif

   do isum=1, PRPLT%replot_nproj_sum
     proj_natom = PRPLT%replot_proj_natom(isum)
     proj_atom  = PRPLT%replot_proj_atom(1:proj_natom, isum)
 
     call get_kunit(PKPTS%kunit, kunit_)
     call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
     call get_e_range(init_e, fina_e, PGEOM%neig, .false., ispinor, flag_erange, init_erange, fina_erange)
     if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

  spin:do is = 1, nspin
         if(flag_proj_sum) then
           c_sum = 0d0
           write(fname_header_sum,'(A,I0)')'band_structure_TBA_atom.sum',isum
           call get_fname(fname_header_sum, fname_sum, is, flag_collinear, flag_noncollinear)
           write(6,'(A)')' '
           write(6,'(3A)',ADVANCE='no')'   WRITING.... projected band structure: ',trim(fname_sum),' ... '
           open(pid_energy+100, file = trim(fname_sum), status = 'unknown')
           if(flag_sparse) then
             write(pid_energy+100, '(A,2(F10.4,A))')'# The EWINDOW mode: energy window [EMIN:EMAX]=[', &
                                                 emin,':', emax,']'
             do ik = 1, nkpoint
               write(pid_energy+100, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
             enddo
           endif
           write(pid_energy+100, '(A, *(I0,1x))'),'#  ATOM_INDEX to be sum up: ', proj_atom(1:proj_natom)
         endif

     atom:do iatom = 1, proj_natom
         ia = proj_atom(iatom)
         imatrix = sum( PGEOM%n_orbital(1:ia) ) - PGEOM%n_orbital(ia) + 1
         write(fname_header,'(A,I0)')'band_structure_TBA_atom.',ia
         call get_fname(fname_header, fname, is, flag_collinear, flag_noncollinear)
         open(pid_energy, file = trim(fname), status = 'unknown')

         if(flag_sparse) then
           write(pid_energy, '(A,2(F10.4,A))')'# The EWINDOW mode: energy window [EMIN:EMAX]=[', &
                                               emin,':', emax,']'
           do ik = 1, nkpoint
             write(pid_energy, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
           enddo
         endif

     eig:do ie = 1, nband ! init_e, fina_e
           write(pid_energy,'(2A,I8,A,I8,3A)',ADVANCE='yes')kmode,'  energy(eV) :',init_e+ie-1,' -th eigen | ',ia, &
                                                      ' -th atom (spec= ',trim(PGEOM%c_spec(PGEOM%spec(ia))),' )'
           if(trim(c_mode) .eq. 'mz') sigma='sigma_z '
           if(trim(c_mode) .eq. 'mx') sigma='sigma_x '
           if(trim(c_mode) .eq. 'my') sigma='sigma_y '

           write(pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
           write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
           do im=imatrix, imatrix + PGEOM%n_orbital(ia) - 1
             write(pid_energy, '(I9)',ADVANCE='NO')im
           enddo
           write(pid_energy,'(A9)',ADVANCE='YES') ' tot'

           if(iatom .eq. proj_natom .and. flag_proj_sum) then
             write(pid_energy+100,'(2A,I8,A      )',ADVANCE='yes')kmode,'  energy(eV) :',init_e+ie-1,' -th eigen '
             if(trim(c_mode) .eq. 'mz') sigma='sigma_z '
             if(trim(c_mode) .eq. 'mx') sigma='sigma_x '
             if(trim(c_mode) .eq. 'my') sigma='sigma_y '

             write(pid_energy+100, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
             write(pid_energy+100, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)     E(ev), '
             write(pid_energy+100,'(A)',ADVANCE='YES') '  tot(atom_sum)'
           endif

        kp:do ik = 1, nkpoint
             if(flag_klinemode) then
               if( ie .le. ne_found(is, ik) ) then
                 write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+nband*(is-1),ik)
               elseif( ie .gt. ne_found(is, ik)) then
                 write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
               endif
             elseif(flag_kgridmode) then
               if( ie .le. ne_found(is, ik) ) then
                 write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+nband*(is-1),ik)
               elseif(ie .gt. ne_found(is, ik)) then
                 write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
               endif
             endif

             if(flag_proj_sum .and. iatom .eq. proj_natom) then
               if(flag_klinemode) then
                 if( ie .le. ne_found(is, ik) ) then
                   write(pid_energy+100,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+nband*(is-1),ik)
                 elseif( ie .gt. ne_found(is, ik)) then
                   write(pid_energy+100,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
                 endif
               elseif(flag_kgridmode) then
                 if( ie .le. ne_found(is, ik) ) then
                   write(pid_energy+100,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+nband*(is-1),ik)
                 elseif(ie .gt. ne_found(is, ik)) then
                   write(pid_energy+100,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
                 endif
               endif
             endif

             if( ie .le. ne_found(is, ik) ) then
               if(flag_print_orbital) then
                 c_tot = 0d0 !initialize
           basis:do im=imatrix, imatrix+PGEOM%n_orbital(ia) - 1
                   if(ispinor .eq. 2) then
                     c_up = V(im,ie,ik)
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') c_up
                     c_tot = c_tot + c_up
                   elseif(ispinor .eq. 1) then
                     c_up = V(im+PGEOM%neig*(is-1),ie+nband*(is-1),ik)
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') c_up
                     c_tot = c_tot + c_up
                   endif
                 enddo basis
                 write(pid_energy,'(*(F9.4))',ADVANCE='YES') c_tot

                 if(flag_proj_sum) c_sum(ie,ik) = c_sum(ie,ik) + c_tot
                 if(flag_proj_sum .and. iatom .eq. proj_natom) then
                   write(pid_energy+100,'(*(F9.4))',ADVANCE='YES') c_sum(ie,ik)
                 endif
               endif
               if(.not.flag_print_orbital) write(pid_energy,*)''
               if(.not.flag_print_orbital .and. flag_proj_sum) write(pid_energy+100,*)'' ! maybe do not need.. but how knows?
             elseif(ie .gt. ne_found(is, ik)) then
               write(pid_energy,*)''
               if(iatom .eq. proj_natom .and. flag_proj_sum) write(pid_energy+100,*)''
             endif
           enddo kp

           write(pid_energy,*)''
           write(pid_energy,*)''
           if(iatom .eq. proj_natom .and. flag_proj_sum) write(pid_energy+100,*)''
           if(iatom .eq. proj_natom .and. flag_proj_sum) write(pid_energy+100,*)''

         enddo eig

         close(pid_energy)
       enddo atom
       if(iatom-1 .eq. proj_natom .and. flag_proj_sum) close(pid_energy+100)
       write(6,'(A)',ADVANCE='YES')' DONE!'
     enddo spin

   enddo


return
endsubroutine
subroutine print_replot_ldos(PRPLT, PINPT, PGEOM)
   use parameters, only: replot, pid_ldos, incar, poscar
   use mpi_setup
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   integer*4       im, ia, iaa, ie, iatom
   integer*4       my_pid_ldos
   real*8          e_range(PRPLT%replot_dos_nediv)
   character*40    filenm
#ifdef MPI
   integer*4       ourjob(nprocs)
   integer*4       ourjob_disp(0:nprocs-1)
   call mpi_job_distribution_chain(PRPLT%replot_ldos_natom, ourjob, ourjob_disp)
#else
   integer*4       ourjob(1)
   integer*4       ourjob_disp(0)
   call mpi_job_distribution_chain(PRPLT%replot_ldos_natom, ourjob, ourjob_disp)
#endif

   e_range = PRPLT%replot_dos_erange

   do ia = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
     iatom = PRPLT%replot_ldos_atom(ia)
     my_pid_ldos = pid_ldos + myid
     write(filenm,'(A,I0,A)') 'LDOS.replot.',iatom,'.dat'
     open(my_pid_ldos, file=trim(filenm), status = 'unknown')

     write(6,'(A)')' '
     write(6,'(3A)',ADVANCE='no')'   WRITING.... local density of states (LDOS): ',trim(filenm),' ... '

     write(my_pid_ldos,'(A,I0,A,A,A )')'# ATOM = ',iatom,' (spec = ',trim(PGEOM%c_spec(PGEOM%spec(iatom))),' )'
     write(my_pid_ldos,'(A,I8,A,F16.8)')'# NDIV = ',PRPLT%replot_dos_nediv,' dE=',e_range(2)-e_range(1)

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

     do ie = 1, PRPLT%replot_dos_nediv
       if(.not.PINPT%flag_collinear) then
         write(my_pid_ldos,'(F16.8,1x,  F16.8    , *(F16.8,1x))')e_range(ie), sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie)), &
                                                                                  PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie)
       elseif(PINPT%flag_collinear) then
         write(my_pid_ldos,'(F16.8,1x,2(F16.8,1x), *(F16.8,1x))')e_range(ie), sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie)), &
                                                                              sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,2,ie)), &
                                                                                  PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie) , &
                                                                                  PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,2,ie)
       endif
     enddo

     close(my_pid_ldos)

   enddo

   return
endsubroutine

subroutine print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart)
   use parameters, only: replot, pid_ldos, incar, poscar
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   integer*4       iatom, is, i
   character*40    filenm
   integer*4       nx, ny, nz
   integer*4       ix, iy, iz
   integer*4       nkp
   real*8          sldos(PGEOM%n_atom,PINPT%nspin) ! integrated spatial LDOS
   real*8          coord(3), coord_cart(3,PGEOM%n_atom), T(3,3), de

   filenm = 'SLDOS.replot.dat' 

   nx = PRPLT%replot_sldos_cell(1)
   ny = PRPLT%replot_sldos_cell(2)
   nz = PRPLT%replot_sldos_cell(3)
   T(:,1) = PGEOM%a_latt(:,1)
   T(:,2) = PGEOM%a_latt(:,2)
   T(:,3) = PGEOM%a_latt(:,3)
   de     = PRPLT%replot_dos_erange(2) - PRPLT%replot_dos_erange(1)

  ! coord info
   do iatom=1,PGEOM%n_atom
     coord(:) = PGEOM%a_coord(:,iatom) + PRPLT%r_origin(:)
     coord(:) = coord(:) - nint(coord(:))
     coord_cart(:,iatom) = coord(1)*PGEOM%a_latt(:,1)+coord(2)*PGEOM%a_latt(:,2)+coord(3)*PGEOM%a_latt(:,3)
   end do

   do is = 1, PINPT%nspin
     do iatom = 1, PGEOM%n_atom
       sldos(iatom,is) = sum(PRPLT%replot_sldos_sum(iatom, is, :)) * de
     enddo
   enddo

   open(pid_ldos, file=trim(filenm), status = 'unknown')
   write(6,'(A)')' '
   write(6,'(3A)',ADVANCE='no')'   WRITING.... spatial local density of states (SLDOS) with energy window integrated: ',trim(filenm),' ... '

   write(pid_ldos,'(A,2(F16.8,A))')'# Integrated spatial local density of states : EWINDOW = [',PRPLT%replot_dos_emin,':',PRPLT%replot_dos_emax,']'
   write(pid_ldos,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom
   
   if(PINPT%flag_collinear) then
     write(pid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)            SLDOS(UP)        SLDOS(DN)        ATOM_SPECIES'
   else
     write(pid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)              SLDOS        ATOM_SPECIES'
   endif

   do iatom = 1, PGEOM%n_atom
     do ix = 1, nx
       do iy = 1, ny
         do iz = 1, nz
           coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                           T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                           T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
           if(PINPT%nspin .eq. 2) then
             write(pid_ldos,'(3(F16.8,1x),2(F16.8,1x),I12)')coord, sldos(iatom, :), PGEOM%spec(iatom)
           else
             write(pid_ldos,'(3(F16.8,1x), (F16.8,1x),I12)')coord, sldos(iatom, :), PGEOM%spec(iatom)
           endif
         enddo
       enddo
     enddo
   enddo

   close(pid_ldos)

   return
endsubroutine

subroutine print_bond(PRPLT, PINPT, PGEOM, coord_cart)
   use parameters, only: replot, pid_ldos, incar, poscar, onsite_tolerance
   use mpi_setup
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   integer*4       nx, ny, nz
   integer*4       ix, iy, iz, jx, jy, jz
   integer*4       ia, ja, iaa, i
   integer*4       natom
   real*8          T(3,3), ti(3), tj(3), DIJ
   character*40    filenm
   real*8          enorm
   external        enorm
   real*8          coord_cart(3,PGEOM%n_atom)
   real*8, allocatable :: bond(:,:), bond_(:,:)
   integer*4       nbond, ibond
   integer*4       nn_class
   real*8          r0
   integer*4       mpierr
#ifdef MPI
   integer*4       ourjob(nprocs)
   integer*4       ourjob_disp(0:nprocs-1)
   integer*4       ourbond(nprocs)
   call mpi_job_distribution_chain(PGEOM%n_atom, ourjob, ourjob_disp)
#else
   integer*4       ourjob(1)
   integer*4       ourjob_disp(0)
   call mpi_job_distribution_chain(PGEOM%n_atom, ourjob, ourjob_disp)
#endif

   filenm = 'BOND.replot.dat'

   nx = PRPLT%replot_sldos_cell(1)
   ny = PRPLT%replot_sldos_cell(2)
   nz = PRPLT%replot_sldos_cell(3)
   T(:,1) = PGEOM%a_latt(:,1)
   T(:,2) = PGEOM%a_latt(:,2)
   T(:,3) = PGEOM%a_latt(:,3)
   natom  = PGEOM%n_atom
   allocate(bond_(8,PGEOM%n_atom * nx * ny * nz * 10))
   bond_ = 0d0
   ibond = 0

   do ia = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
     do ix = 1, nx
       do iy = 1, ny
         do iz = 1, nz
           do ja = ia, natom
             do jx = 1, nx
               do jy = 1, ny
                 do jz = 1, nz

                   ti(:)=(/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                           T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                           T(3,1)*(ix-1) + T(3,2)*(iy-1) + T(3,3)*(iz-1)/)
                   tj(:)=(/T(1,1)*(jx-1) + T(1,2)*(jy-1) + T(1,3)*(jz-1), &
                           T(2,1)*(jx-1) + T(2,2)*(jy-1) + T(2,3)*(jz-1), &
                           T(3,1)*(jx-1) + T(3,2)*(jy-1) + T(3,3)*(jz-1)/)
                   DIJ = enorm(3, (coord_cart(:,ia) + ti(:)) - (coord_cart(:,ja) + tj(:)) )
                   call get_nn_class(PGEOM, ia, ja, DIJ, onsite_tolerance, nn_class, r0)
                   if(nn_class .eq. 1 .and. DIJ .le. PRPLT%bond_cut) then
                     ibond = ibond + 1
                     bond_(1:3, ibond) = coord_cart(1:3, ia) + ti(:)
                     bond_(4:6, ibond) = coord_cart(1:3, ja) + tj(:)
                     bond_(7:8, ibond) = (/real(ia),real(ja)/)
                   endif
                 enddo !jz
               enddo !jy
             enddo !jx
           enddo !j
         enddo !iz
       enddo !iy
     enddo !ix
   enddo  !i

#ifdef MPI
   call MPI_ALLREDUCE(ibond, nbond, 1, MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
   allocate(bond(8,nbond))
   call MPI_ALLGATHER(ibond, 1, MPI_INTEGER4, ourbond,1,MPI_INTEGER4, mpi_comm_earth, mpierr)  
   if_main_then
     do i = 1, nprocs - 1
       ourjob_disp(i) = ourjob_disp(i-1) + 8 * ourbond(i)
     enddo
   if_main_end
   call MPI_BCAST(ourjob_disp, size(ourjob_disp), MPI_INTEGER4, 0, mpi_comm_earth, mpierr)
   call MPI_GATHERV(bond_(1:8,1:ibond), size(bond_(1:8,1:ibond)), MPI_REAL8, bond, &
                    ourbond * 8, ourjob_disp,  MPI_REAL8, 0, &
                    mpi_comm_earth, mpierr)
#else
   nbond = ibond
   allocate(bond(8,nbond))
   bond(1:8,1:nbond) = bond_(1:8, 1:nbond)
#endif

   if_main_then
     open(pid_ldos+777, file = trim(filenm), status = 'unknown')
     do ibond = 1, nbond
       write(pid_ldos+777,'(3(F16.8,1x),I8)') bond(1:3, ibond), nint(bond(7,ibond))
       write(pid_ldos+777,'(3(F16.8,1x),I8)') bond(4:6, ibond), nint(bond(8,ibond))
       write(pid_ldos+777,'(A)')'  '
       write(pid_ldos+777,'(A)')'  '
     enddo
     close(pid_ldos+777)
   if_main_end

   return
endsubroutine
subroutine check_band_header(flag_vector, flag_wf, flag_sparse, nemax, emin, emax, nspin, &
                             flag_erange, init_erange, fina_erange, c_mode, flag_formatted)
   use parameters, only : pid_energy
   use mpi_setup
   implicit none
   integer*4              pid
   integer*4              nspin
   integer*4              i_continue, mpierr
   integer*4              nemax
   integer*4              init_erange, fina_erange
   integer*4              ikmode_, neig_, nkpoint_, nband_  ! read from band file but not stored, 
   integer*4              ispin_, nspin_, ispinor_          ! if not matched with defined settings, calculation will be crashed.
   character*132          inputline
   character*80           fnameu, fnamed
   character*12           c_emin, c_emax
   character*12           c_nemax
   character*2            c_mode
   real*8                 emin, emax
   logical                flag_sparse, flag_erange
   logical                flag_wf
   logical                flag_vector
   logical                flag_formatted
   logical                flag_go
   logical                flag_exist
   logical                flag_single_
   character*11           form_
   flag_sparse = .false.
   flag_erange = .false.

   flag_go     = .false.

   if(flag_formatted)      form_ = 'formatted'
   if(.not.flag_formatted) form_ = 'unformatted'

   pid = pid_energy
   if(nspin .eq. 2) then ! collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA.up.dat'
       fnamed = 'band_structure_TBA.dn.dat'
       inquire(file=trim(fnameu),exist=flag_exist)
       if(flag_exist) then 
         inquire(file=trim(fnamed),exist=flag_exist)
         if(.not. flag_exist) then
           if_main write(6,'(3A)')'   !ERROR file: ', trim(fnamed),' does not exist!!'
           if_main write(6,'(A)') '    Exit program... : check_band_header'
           kill_job
         endif
       elseif(.not. flag_exist) then
         if_main write(6,'(3A)')'   !ERROR file: ', trim(fnameu),' does not exist!!'
         if_main write(6,'(A)') '    Exit program... : check_band_header'
         kill_job
       endif
       
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA.up.bin'
       fnamed = 'band_structure_TBA.dn.bin'
       inquire(file=trim(fnameu),exist=flag_exist)
       if(flag_exist) then
         inquire(file=trim(fnamed),exist=flag_exist)
         if(.not. flag_exist) then
           if_main write(6,'(3A)')'   !ERROR file: ', trim(fnamed),' does not exist!!'
           if_main write(6,'(A)') '    Exit program... : check_band_header'
           kill_job
         endif
       elseif(.not. flag_exist) then
         if_main write(6,'(3A)')'   !ERROR file: ', trim(fnameu),' does not exist!!'
         if_main write(6,'(A)') '    Exit program... : check_band_header'
         kill_job
       endif

     endif
   elseif(nspin .eq. 1) then ! nm or non-collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA.dat'
       inquire(file=trim(fnameu),exist=flag_exist)
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA.bin'
       inquire(file=trim(fnameu),exist=flag_exist)
     endif
     if(.not. flag_exist) then
       if_main write(6,'(3A)')'   !ERROR file: ', trim(fnameu),' does not exist!!'
       if_main write(6,'(A)') '    Exit program... : check_band_header'
       kill_job
     endif
   endif


   ! it is sufficient to check with is = 1 case only.
   if(nspin .eq. 1) then
     if_main write(6,'(2A)')'   LOADING.... ', trim(fnameu)
   elseif(nspin .eq. 2) then
     if_main write(6,'(3A)')'   LOADING.... ', trim(fnameu), ' and ', trim(fnamed)
   endif

   open(pid, file=trim(fnameu), form=trim(form_), status='old', iostat=i_continue)
   if(i_continue .ne. 0) then
     write(6,'(2A)')'  !!! ERROR IN READING ', trim(fnameu)
     write(6,'(A)')'      Exit program... : check_band_header '
     kill_job
   endif

   if(flag_formatted) then
     ! check whether sparse eigen solver was used... If yes, read emin, emax
     ! Read first two lines for the basic informations.
     do while (.not. flag_go)
       read(pid,'(A)')inputline
       if(index(inputline,'EWINDOW') .ge. 1) then
         flag_sparse = .true.
         if(flag_sparse) then
           call strip_off(trim(inputline), c_emin, '[EMIN:EMAX]=[',':',1)
           call strip_off(trim(inputline), c_emax, ':',']',1)
           call str2real(c_emin,emin)
           call str2real(c_emax,emax)
           call strip_off(trim(inputline), c_nemax,'NE_MAX=',' ',5)
           call str2int(c_nemax, nemax)
           if_main write(6,'(A,2(F12.4,A),I0)')'   => EWINDOW [EMIN:EMAX]=[ ',emin,' : ',emax,' ], NE_MAX= ',nemax
         endif
       elseif(index(inputline, 'ERANGE') .ge. 1) then
         flag_erange = .true.
         if(flag_erange) then
           call strip_off(trim(inputline), c_emin, 'ERANGE=[',':',1)
           call strip_off(trim(inputline), c_emax, ':',']',1)
           call str2int(c_emin, init_erange)
           call str2int(c_emax, fina_erange)
           nemax = fina_erange - init_erange + 1
           if_main write(6,'(A,I0,A,I0,A,I0)')'   => ERANGE=[ ',init_erange,' : ',fina_erange,' ], NERANGE= ',nemax
         endif
       elseif(index(inputline, 'LORBIT=') .ge. 1) then
         call strip_off(trim(inputline), c_mode, 'LORBIT=[',']',1)
         if(trim(c_mode) .eq. 'rh') then
           flag_wf = .false.
           flag_vector = .true.
           if_main write(6,'(A,I0,A,I0,A,I0)')'   => C_MODE= rh ; <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j'
         elseif(trim(c_mode) .eq. 'wf') then
           flag_wf = .true.
           flag_vector = .true.
           if_main write(6,'(A,I0,A,I0,A,I0)')'   => C_MODE= wf ; wavefunction coefficent will be read.'
         else
           flag_wf = .false.
           flag_vector = .false.
           if_main write(6,'(A,I0,A,I0,A,I0)')'   => C_MODE= no ; only read eigenvalues.'
         endif
       else
         flag_go = .true.
       endif

     enddo

   elseif(.not. flag_formatted) then
     read(pid) ikmode_, flag_vector, flag_single_, flag_erange, flag_sparse, &
               neig_, nkpoint_, nband_, ispin_, nspin_, ispinor_, c_mode
     if(flag_erange) then
       read(pid) flag_erange, init_erange, fina_erange
       nemax = fina_erange - init_erange + 1
       if_main write(6,'(A,I0,A,I0,A,I0)')'   => ERANGE=[ ',init_erange,' : ',fina_erange,' ], NERANGE= ',nemax  
       if(nemax .ne. nband_) then 
         write(6,'(3A)')'         !ERROR in reading ',trim(fnameu),' : nemax .ne. nband_ '
         write(6,'(A)') '          Exit program...'
         kill_job
       endif
     else
       read(pid) flag_erange ! .FALSE.
     endif
     if(flag_sparse) then
       read(pid) flag_sparse, emin, emax, nemax
       if(nemax .ne. nband_) then
         write(6,'(3A)')'         !ERROR in reading ',trim(fnameu),' : nemax .ne. nband_ '
         write(6,'(A)') '          Exit program...'
         kill_job
       endif
       if_main write(6,'(A,2(F12.4,A),I0)')'   => EWINDOW [EMIN:EMAX]=[ ',emin,' : ',emax,' ], NE_MAX= ',nemax
     else
       read(pid) flag_sparse ! .FALSE.
     endif

     if(c_mode .eq. 'wf') then
       flag_wf = .true.  ! current version only support reading band_....bin file with wavefuction information is written.
       if_main write(6,'(A)')'   => C_MODE= wf ; wavefunction coefficent will be read.'
     elseif(c_mode .eq. 'rh') then
       flag_wf = .false.
       if_main write(6,'(A)')'   => C_MODE= rh ; <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j'
     elseif(c_mode .eq. 'no') then
       flag_wf = .false. ! current version only support reading band_....bin file with wavefuction information is written.
       if_main write(6,'(A)')'   => C_MODE= no ; only read eigenvalues.'
     else
       if_main write(6,'(5A)')'    !ERROR in reading ',trim(fnameu),' : cannot recognize C_MODE=',c_mode, &
                                          ' in current version. Exit program...'
       kill_job
     endif
   endif

   close(pid)

   if(flag_vector .and. trim(c_mode) .eq. 'no') then 
     if(nspin .eq. 2) then
       if_main write(6,'(4A)')'    !ERROR in reading ',trim(fnameu), ' and ', trim(fnamed)
     elseif(nspin .eq. 1) then
       if_main write(6,'(2A)')'    !ERROR in reading ',trim(fnameu)
     endif
     if_main write(6,'(A)') '     cannot read orbital information from avove files. Check again...'
     if_main write(6,'(A)') '     Exit program...'
     kill_job
   endif

   return
endsubroutine
subroutine V_to_V2(V, V2, ispin, nspin, nbasis, nband, nkp)
   implicit none
   integer*4    ik, is, ie, im
   integer*4    ispin, nspin, nbasis, nband, nkp
   complex*16   V(nbasis*ispin, nband*nspin, nkp)   
   real*8       V2(nbasis, nband*nspin, nkp)  
   complex*16   c_up, c_dn

   do ik = 1, nkp
     do is = 1, nspin
       do ie = 1, nband 

         if(nspin .eq. 1 .and. ispin .eq. 2) then ! non-collinear
           do im = 1, nbasis
             c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
             V2(im, ie, ik) = real(conjg(c_up)*c_up + conjg(c_dn)*c_dn)
           enddo
         elseif(nspin .eq. 1 .and. ispin .eq. 1) then  !nonmagnetic
           do im = 1, nbasis
             c_up = V(im,ie,ik)
             V2(im, ie, ik) = real(conjg(c_up)*c_up)
           enddo
         elseif(nspin .eq. 2 .and. ispin .eq. 2) then ! collinear
           do im = 1, nbasis
             c_up = V(im+(is-1)*nbasis,ie+(is-1)*nband,ik)
             V2(im, ie+(is-1)*nband, ik) = real(conjg(c_up)*c_up)
           enddo
         endif

       enddo 
     enddo
   enddo

   return
endsubroutine
