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
   integer*4                   norb, nband, nspin
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
   logical                     flag_replot_proj_band
   logical                     flag_erange
   logical                     flag_wf
   real*8,     allocatable  :: E(:,:)
!  real*8                      E(PINPT%nband*PINPT%nspin,PKPTS%nkpoint)
   ! be careful that PINPT%nband should be same as nband which will be found by check_band_header routine...
   ! this routine should be modified in this sense...  KHJ May. 15 2019.
   real*8                      V(PGEOM%neig,PINPT%nband*PINPT%nspin,PKPTS%nkpoint) ! not wavevector, 
   complex*16                  V2(PGEOM%neig*PINPT%ispin,PINPT%nband*PINPT%nspin,PKPTS%nkpoint) ! not wavevector, 
!  real*8,     allocatable  :: V2(:,:,:)  ! not wavevector, but its sqare 
!  complex*16, allocatable  :: V(:,:,:)  ! wavevector
   real*8                      time2, time1
   character*80                cjob
   character*12                c_mode
   logical                     flag_sparse

   if_main call time_check(time2,time1,'init')

   flag_replot_dos       = PRPLT%flag_replot_dos
   flag_replot_ldos      = PRPLT%flag_replot_ldos
   flag_replot_sldos     = PRPLT%flag_replot_sldos
   flag_replot_proj_band = PRPLT%flag_replot_proj_band
   if(PKPTS%flag_klinemode) kmode = 'line'
   if(PKPTS%flag_kgridmode) kmode = 'grid'
   if(flag_replot_ldos .or. flag_replot_sldos .or. flag_replot_proj_band) then
     flag_vector = .true.
   elseif(.not. flag_replot_ldos .and. .not. flag_replot_sldos .and. .not. flag_replot_proj_band) then
     flag_vector = .false.
   endif

   if_main write(6,*)''
   if_main write(6,*)''
   if_main write(6,'(A)')' ** Program run in REPLOT mode '
   if_main write(6,*)''
   cjob = ' '
   if(flag_replot_dos)      write(cjob,'(A,A)') trim(cjob), ' + DOS'
   if(flag_replot_ldos)     write(cjob,'(A,A)') trim(cjob), ' + LDOS'
   if(flag_replot_sldos)    write(cjob,'(A,A)') trim(cjob), ' + SLDOS'
   if(flag_replot_proj_band)write(cjob,'(A,A)') trim(cjob), ' + PROJ_BAND'

   if_main write(6,'(A,A,A)')'  START REPLOT: ',trim(cjob),' EVALUATION'
   if_main write(6,*)''

   nspin                 = PINPT%nspin
   norb                  = PGEOM%neig
   call check_band_header(flag_vector, flag_wf, flag_sparse, nband, emin_band, emax_band, nspin, flag_erange, init_erange, fina_erange, c_mode)

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

   if(nband .ne. PINPT%nband) then
     write(6,'(A)')'    ! ERROR: nband read from band structure file is differ from the setup PINPT%nband. This discrepancy may'
     write(6,'(A)')'           : be originated from NE_MAX tag of your band structure that has been calculated under EWINDOW'
     write(6,'(A)')'           : tag which uses sparse matrix solver. Exit program...'
     kill_job
   endif
!  if_main allocate(V2(norb,nband*nspin,nkp))
!  if(flag_wf) then
!    if_main allocate(V(norb*ispin,nband*nspin,nkp))
!  endif
!  if(flag_vector) then
!    if_main allocate(V2(norb,nband*nspin,nkp))
!  endif
   if(flag_wf) then
     if(flag_vector) then
!      if_main allocate(V(norb*ispin,nband*nspin,nkp))
       if_main call report_memory(int8(size(V)), 16, 'Wave vector  ') 
     endif
     if_main call read_energy_tbfit_V(E, V, ne_found, ispin, nspin, norb, nband, nkp, kmode, flag_vector, flag_wf)
   elseif(.not. flag_wf) then
     if_main call read_energy_tbfit_V2(E, V2, ne_found, ispin, nspin, norb, nband, nkp, kmode, flag_vector, flag_wf)
   endif
   if(flag_vector .and. flag_wf .and. myid .eq. 0) then
     if_main call V_to_V2(V, V2, ispin, nspin, norb, nband, nkp)
   endif

   allocate( PRPLT%replot_dos_erange(nediv) )
   PRPLT%replot_dos_erange = e_range

#ifdef MPI
   call MPI_BCAST(E, size(E), MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

   if(flag_replot_dos .or. flag_replot_ldos .or. flag_replot_sldos) then
     call get_dos_ldos(PINPT, PGEOM, PRPLT, E, V2, nspin, norb, nband, nkp, e_range, nediv, sigma, &
                       flag_replot_dos, flag_replot_ldos, flag_replot_sldos)
   endif   

   if(flag_replot_proj_band) then
     if_main call print_replot_energy_proj(PKPTS, E, V2, PGEOM, PRPLT, ne_found, emin_band, emax_band, ispinor, nband, nspin, ispin, norb, &
                                   flag_erange, flag_vector, flag_sparse, init_erange, fina_erange,c_mode)
   endif

   if_main call time_check(time2,time1)
   if_main write(6,*)''
   if_main write(6,'(3A,F10.4,A)')'  END REPLOT: ',trim(cjob),' EVALUATION : ',time2, ' (sec)'
   if_main write(6,*)''


!  if(allocated(V)) deallocate(V)
!  if(allocated(V2)) deallocate(V2)

   return
endsubroutine

subroutine get_dos_ldos(PINPT, PGEOM, PRPLT, E, V2, nspin, norb, nband, nkp, e_range, nediv, sigma, &
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
   integer*4                   norb, nband, nspin, nkp
   integer*4                   nediv
   integer*4                   mpierr
   real*8                      e_range(nediv), de, sigma
   real*8                      E(nband*nspin,nkp)
   real*8                      V2(norb,nband*nspin,nkp) ! not wavevector, 
   real*8                      myV(norb,nband*nspin)
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
     if_main call report_memory(int8(size(myV))*nprocs + size(V2), 8, 'Wave vector2 ')
   endif

   call mpi_job_distribution_chain(nediv, ourjob, ourjob_disp)
   if_main write(6,*)''
   do ik = 1, nkp
     if(nkp .lt. 10) then
       if_main write(6,'(A,I0,A,I0)')  '         STAT KP: ', ik,'/',nkp
     else
       if( ik/real(nkp)*100d0 .ge. real(iaddk*inck) ) then
         if_main write(6,'(A,F10.3,A)')'         STAT KP: ', ik/real(nkp)*100d0, ' %'
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
           if_main write(6,'(A,F10.3,A)')'            STAT EN: ', (ie-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0, ' %'
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
     if_main write(6,'(A)')'      write DOS ... Done!'
   endif
   if(flag_replot_ldos) then
     call MPI_ALLREDUCE(replot_ldos_tot, PRPLT%replot_ldos_tot, size(replot_ldos_tot), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
     call print_replot_ldos(PRPLT, PINPT, PGEOM)
     if_main write(6,'(A)')'      write LDOS ... Done!'
   endif
   if(flag_replot_sldos) then
     call MPI_REDUCE(replot_sldos_sum, PRPLT%replot_sldos_sum, size(replot_sldos_sum), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
     if_main call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart)
     call MPI_BCAST(coord_cart, size(coord_cart), MPI_REAL8, 0, mpi_comm_earth, mpierr)
     call print_bond(PRPLT, PINPT, PGEOM, coord_cart)
     if_main write(6,'(A)')'      write SLDOS ... Done!'
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
     write(6,'(A)')'      write DOS ... Done!'
   endif
   if(flag_replot_ldos) then
     PRPLT%replot_ldos_tot = replot_ldos_tot
     call print_replot_ldos(PRPLT, PINPT, PGEOM)
     write(6,'(A)')'      write LDOS ... Done!'
   endif
   if(flag_replot_sldos) then
     PRPLT%replot_sldos_sum = replot_sldos_sum
     call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart)
     write(6,'(A)')'      write SLDOS ... Done!'
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
subroutine print_replot_energy_proj(PKPTS,E,V,PGEOM,PRPLT, ne_found, emin, emax, ispinor, nband, nspin, ispin, norb, &
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
   integer*4       ispin_print, nbasis, ispinor, ispin
   integer*4       nband, nspin, norb
   integer*4       nkpoint
   integer*4       proj_natom, proj_atom(maxval(PRPLT%replot_proj_natom(1:PRPLT%replot_nproj_sum)))
   integer*4       init_erange, fina_erange
   integer*4       init_e, fina_e
   integer*4       emin, emax
   integer*4       ne_found(nspin, PKPTS%nkpoint)
   integer*4       imatrix
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   logical         flag_proj_sum
   logical         flag_erange, flag_vector
   logical         flag_sparse
   logical         flag_collinear, flag_noncollinear
   real*8          E(nband*nspin,PKPTS%nkpoint)
   real*8          V(norb,nband*nspin,PKPTS%nkpoint) ! rho_nk = <psi_nk|psi_nk> is stored, not complex wavefunction
   real*8          c_up, c_dn, c_tot
   real*8          c_sum(nband,PKPTS%nkpoint)
   character*80    fname_header, fname_header_sum
   character*80    fname, fname_sum
   character*6     kunit_
   character*28    kmode
   character*8     sigma
   character*12    c_mode

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
     call get_ispin_print(flag_collinear, ispin_print)
     if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

  spin:do is = 1, ispin_print
         if(flag_proj_sum) then
           c_sum = 0d0
           write(fname_header_sum,'(A,I0)')'band_structure_TBA_atom.sum',isum
           call get_fname(fname_header_sum, fname_sum, is, flag_collinear, flag_noncollinear)
           write(6,'(3A)',ADVANCE='no')'      writing projected band structure... ',trim(fname_sum),' ... '
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
subroutine check_band_header(flag_vector, flag_wf, flag_sparse, nemax, emin, emax, nspin, flag_erange, init_erange, fina_erange, c_mode)
   use parameters, only : pid_energy
   use mpi_setup
   implicit none
   integer*4              pid
   integer*4              nspin
   integer*4              i_continue, mpierr
   integer*4              nemax
   integer*4              init_erange, fina_erange
   character*132          inputline
   character*80           fnameu, fnamed
   character*12           c_emin, c_emax
   character*12           c_nemax
   character*12           c_mode
   real*8                 emin, emax
   logical                flag_sparse, flag_erange
   logical                flag_wf
   logical                flag_vector

   flag_sparse = .false.
   flag_erange = .false.

   pid = pid_energy
   if(nspin .eq. 2) then
     fnameu = 'band_structure_TBA.up.dat'
     fnamed = 'band_structure_TBA.dn.dat'
   elseif(nspin .eq. 1) then
     fnameu = 'band_structure_TBA.dat'
   endif

   ! it is sufficient to check with is = 1 case only.
   if(nspin .eq. 1) then
     if_main write(6,'(2A)')'  LOADING.... ', trim(fnameu)
   elseif(nspin .eq. 2) then
     if_main write(6,'(3A)')'  LOADING.... ', trim(fnameu), ' and ', trim(fnamed)
   endif

   open(pid, file=trim(fnameu), iostat=i_continue)
   if(i_continue .ne. 0) then
     write(6,'(A)')'  !!! ERROR IN READING (up) ', trim(fnameu)
     write(6,'(A)')'      EXIT PROGRAM : check_band_header '
     stop
   endif

   ! check whether sparse eigen solver was used... If yes, read emin, emax
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
       if_main write(6,'(A,2(F12.4,A),I0)')'         EWINDOW [EMIN:EMAX]=[ ',emin,' : ',emax,' ], NE_MAX= ',nemax
     endif
   elseif(index(inputline, 'ERANGE') .ge. 1) then
     flag_erange = .true.
     if(flag_erange) then
       call strip_off(trim(inputline), c_emin, 'ERANGE=[',':',1)
       call strip_off(trim(inputline), c_emax, ':',']',1)
       call str2int(c_emin, init_erange)
       call str2int(c_emax, fina_erange)
       nemax = fina_erange - init_erange + 1
       if_main write(6,'(A,I0,A,I0,A,I0)')'         ERANGE=[ ',init_erange,' : ',fina_erange,' ], NERANGE= ',nemax
     endif
   elseif(index(inputline, 'LORBIT=') .ge. 1) then
     call strip_off(trim(inputline), c_mode, 'LORBIT=[',']',1)
     if(trim(c_mode) .eq. 'rh') then
       flag_wf = .false.
     elseif(trim(c_mode) .eq. 'wf') then
       flag_wf = .true.
     endif
   endif

   read(pid,'(A)') inputline
   if(index(inputline,'LORBIT=') .ge. 1) then
     call strip_off(trim(inputline), c_mode, 'LORBIT=[',']',1)
     if(trim(c_mode) .eq. 'rh') then
       flag_wf = .false.
     elseif(trim(c_mode) .eq. 'wf') then
       flag_wf = .true.
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
subroutine V_to_V2(V, V2, ispin, nspin, norb, nband, nkp)
   implicit none
   integer*4    ik, is, ie, im
   integer*4    ispin, nspin, norb, nband, nkp
   complex*16   V(norb*ispin, nband*nspin, nkp)   
   real*8       V2(norb, nband*nspin, nkp)  
   complex*16   c_up, c_dn

   do ik = 1, nkp
     do is = 1, nspin
       do ie = 1, nband 

         if(nspin .eq. 1 .and. ispin .eq. 2) then ! non-collinear
           do im = 1, norb
             c_up = V(im,ie,ik); c_dn = V(im + norb,ie,ik)
             V2(im, ie, ik) = real(conjg(c_up)*c_up + conjg(c_dn)*c_dn)
           enddo
         elseif(nspin .eq. 1 .and. ispin .eq. 1) then  !nonmagnetic
           do im = 1, norb
             c_up = V(im,ie,ik)
             V2(im, ie, ik) = real(conjg(c_up)*c_up)
           enddo
         elseif(nspin .eq. 2 .and. ispin .eq. 2) then ! collinear
           do im = 1, norb
             c_up = V(im+(is-1)*norb,ie+(is-1)*nband,ik)
             V2(im, ie+(is-1)*nband, ik) = real(conjg(c_up)*c_up)
           enddo
         endif

       enddo 
     enddo
   enddo

   return
endsubroutine
