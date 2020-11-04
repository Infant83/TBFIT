#include "alias.inc"
subroutine get_replot(PINPT, PGEOM, PKPTS, PRPLT)
   use parameters, only : incar, poscar, replot, kpoints
   use mpi_setup
   implicit none
   type(incar)   :: PINPT
   type(poscar ), dimension(PINPT%nsystem) :: PGEOM
   type(kpoints), dimension(PINPT%nsystem) :: PKPTS
   type(replot ), dimension(PINPT%nsystem) :: PRPLT
   integer*4                                  i
   integer*4                                  mpierr

   do i = 1, PINPT%nsystem
     if(.not. PRPLT(i)%flag_replot) cycle

     call replot_dos_band(PINPT, PGEOM(i), PKPTS(i), PRPLT(i))

   enddo

   return
endsubroutine

subroutine replot_dos_band(PINPT, PGEOM, PKPTS, PRPLT)
   use parameters, only : incar, poscar, replot, kpoints, pi, pi2
   use mpi_setup
   use time
   use memory
   use print_io
   implicit none
   type(incar  )            :: PINPT
   type(poscar )            :: PGEOM
   type(replot )            :: PRPLT
   type(kpoints)            :: PKPTS
   integer*4                   nbasis, nband, nspin
   integer*4                   ispinor, ispin, nkp
   integer*4                   nediv
   integer*4                   i, k
   integer*4                   ne_found(PINPT%nspin, PKPTS%nkpoint)
   integer*4                   init_erange, fina_erange
   real*8                      emin, emax, e_range(PRPLT%replot_dos_nediv)
   real*8                      emin_band, emax_band
   real*8                      sigma
   integer*4                   mpierr
   character*4                 kmode
   logical                     flag_vector
   logical                     flag_replot_dos, flag_replot_ldos, flag_replot_sldos, flag_replot_didv
   logical                     flag_replot_proj_band, flag_replot_band
   logical                     flag_erange
   logical                     flag_wf
   logical                     flag_formatted
   logical                     flag_exit
   real*8,     allocatable  :: E(:,:) ! (nband * nspin, nkpoint)
   ! be careful that PGEOM%nband should be same as nband which will be found by check_band_header routine...
   ! this routine should be modified in this sense...  KHJ May. 15 2019.
   complex*16, allocatable  :: V(:,:,:)
   real*8,     allocatable  :: V2(:,:,:)
   real*8                      time2, time1
   character*132               cjob
   character*2                 c_mode, c_mode_
   character*2                 c_mode_print(PRPLT%replot_nband)
   logical                     flag_sparse

   flag_exit = .false.

   if_main call time_check(time2,time1,'init')

   flag_replot_dos       = PRPLT%flag_replot_dos
   flag_replot_ldos      = PRPLT%flag_replot_ldos
   flag_replot_sldos     = PRPLT%flag_replot_sldos
   flag_replot_didv      = PRPLT%flag_replot_didv
   flag_replot_proj_band = PRPLT%flag_replot_proj_band
   flag_replot_band      = PRPLT%flag_replot_band
   flag_formatted        = PRPLT%flag_replot_formatted
   c_mode_print          = PRPLT%replot_axis_print_mag

   if(PKPTS%flag_klinemode) kmode = 'line'
   if(PKPTS%flag_kgridmode) kmode = 'grid'

   write(message,*)''  ; write_msg
   write(message,*)''  ; write_msg
   write(message,'(A)')' ** Program run in REPLOT mode '  ; write_msg
   write(message,*)''  ; write_msg
   cjob = ' '

   nspin                 = PINPT%nspin
   nbasis                = PGEOM%neig
   call check_band_header(PRPLT, PINPT, flag_vector, flag_wf, flag_sparse, nband, emin_band, emax_band, nspin, &
                          flag_erange, init_erange, fina_erange, c_mode, flag_formatted)

   if(flag_replot_dos)      write(cjob,'(A,A)') trim(cjob), ' + DOS'
   if(flag_vector) then
     if(flag_replot_ldos)     write(cjob,'(A,A)') trim(cjob), ' + LDOS'
     if(flag_replot_sldos)    write(cjob,'(A,A)') trim(cjob), ' + SLDOS'
     if(flag_replot_didv)     write(cjob,'(A,A)') trim(cjob), ' + dI/dV'
     if(flag_replot_proj_band)write(cjob,'(A,A)') trim(cjob), ' + PROJ_BAND'
   endif

   if(flag_replot_band) then
     do i=1, PRPLT%replot_nband
       if(c_mode_print(i) .eq. 'mx' .or. c_mode_print(i) .eq. 'my' .or. c_mode_print(i) .eq. 'mz') then
                              write(cjob,'(A,4A)')trim(cjob), ' + BAND (<sigma_i>, i=',c_mode_print(i),')'
       elseif(c_mode_print(i) .eq. 'wf') then
                              write(cjob,'(A,2A)')trim(cjob), ' + BAND (wavefunction)'
       elseif(c_mode_print(i) .eq. 'rh') then
                              write(cjob,'(A,2A)')trim(cjob), ' + BAND (<phi_i|psi_nk>, i=orbital)'
       elseif(c_mode_print(i) .eq. 'no') then
                              write(cjob,'(A,2A)')trim(cjob), ' + BAND (eig only)'
       endif
     enddo
   endif

   write(message,*)''  ; write_msg
   write(message,'(A,A,A)')'   START REPLOT: ',trim(cjob),' EVALUATION'  ; write_msg
   write(message,*)''  ; write_msg


   if(.not.flag_sparse) nband = PGEOM%nband
   ispinor               = PINPT%ispinor
   ispin                 = PINPT%ispin
   nkp                   = PKPTS%nkpoint
   nediv                 = PRPLT%replot_dos_nediv
   emax                  = PRPLT%replot_dos_emax
   emin                  = PRPLT%replot_dos_emin
   if(nediv .gt. 1) then
     e_range               = emin + (/(k, k=0,nediv-1)/) * (emax - emin)/dble(nediv - 1)
   elseif(nediv .eq. 1) then
     e_range               = (/emin,emax/)
   endif
   sigma                 = PRPLT%replot_dos_smearing
   
   allocate(E(nband*nspin,nkp))

   if(flag_vector) then
     allocate(V2(nbasis,nband*nspin,nkp))
     if(flag_wf) then
       allocate(V(nbasis*ispin,nband*nspin,nkp))
     endif
   endif

   if(nband .ne. PGEOM%nband) then
     write(message,'(A)')'    ! ERROR: nband read from band structure file is differ from the setup PGEOM%nband. This discrepancy may' ; write_msg
     write(message,'(A)')'           : be originated from NE_MAX tag of your band structure that has been calculated under EWINDOW' ; write_msg
     write(message,'(A)')'           : tag which uses sparse matrix solver. Exit program...' ; write_msg
     kill_job
   endif

   ! reading energy + wave vector (if requested)
   if(flag_wf) then
     if_main call read_energy_tbfit_V(PRPLT, PINPT, E, V, ne_found, ispin, nspin, nbasis, nband, nkp, kmode, &
                                      flag_vector, flag_wf, flag_formatted, flag_exit)
#ifdef MPI
     call MPI_BCAST(flag_exit, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
     call MPI_BCAST(ne_found,size(ne_found), MPI_INTEGER, 0, mpi_comm_earth, mpierr)
#endif
     if(flag_exit) then
       kill_job
     endif

   elseif(.not. flag_wf) then
   
     if(.not. flag_formatted .and. PRPLT%flag_replot_band) then
       do i = 1, PRPLT%replot_nband
         c_mode_ = c_mode_print(i)
         if(c_mode_(1:1) .eq. 'm' .or. c_mode_(1:1) .eq. 'w') then
           write(message,'(A)') '    ! WARN: reading with C_MODE:"rh" is not available since you have '  ; write_msg
           write(message,'(3A)')'            requested to replot ', c_mode_print(i), ' in your REPLOT_BAND tag'   ; write_msg
           write(message,'(A)') '            which requires to read band_structure_TBA file with wavefunction "wf" information.'  ; write_msg
           write(message,'(A)') '      Exit program... '  ; write_msg
           kill_job
         endif
       enddo
     endif
     if(flag_vector) then
       if_main call read_energy_tbfit_V2(PRPLT, PINPT, E, V2, ne_found, ispin, nspin, nbasis, nband, nkp, kmode, &
                                         flag_vector, flag_wf, flag_formatted, flag_exit)
     elseif(.not. flag_vector) then
       if_main call read_energy_tbfit(PRPLT, PINPT, E, ne_found, ispin, nspin, nbasis, nband, nkp, kmode, flag_formatted, flag_exit)
     endif

#ifdef MPI
       call MPI_BCAST(flag_exit, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
       call MPI_BCAST(ne_found,size(ne_found), MPI_INTEGER, 0, mpi_comm_earth, mpierr)
#endif
     if(flag_exit) then
       kill_job
     endif

   endif

   if(flag_wf) then
     if(flag_vector) then
       if_main call report_memory(int8(size(V)), 16, 'Wave vector   ')
     endif
   elseif(.not. flag_wf) then
     if(flag_vector) then
       if_main call report_memory(int8(size(V2)), 8, 'Wave vector   ')
     endif
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
   if(flag_replot_dos .or. flag_replot_ldos .or. flag_replot_sldos .or. flag_replot_didv) then
     if(flag_vector) then
     call get_dos_ldos(PINPT, PGEOM, PRPLT, ne_found, E, V2, nspin, nbasis, nband, nkp, e_range, nediv, sigma, &
                       flag_replot_dos, flag_replot_ldos, flag_replot_sldos, flag_replot_didv, flag_vector)
     elseif(.not. flag_vector) then
       call get_dos_ldos(PINPT, PGEOM, PRPLT, ne_found, E, 0d0, nspin, nbasis, nband, nkp, e_range, nediv, sigma, &
                         flag_replot_dos, flag_replot_ldos, flag_replot_sldos, flag_replot_didv, flag_vector)
     endif
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
     if(flag_wf) then
       if_main call print_replot_energy(PINPT, PKPTS, E, V, PGEOM, PRPLT, ne_found, emin_band, emax_band, &
                                             ispinor, nband, nspin, ispin, nbasis, &
                                             flag_erange, flag_vector, flag_sparse, &
                                             init_erange, fina_erange) !, c_mode_print)
     elseif(.not. flag_wf) then
       if_main call print_replot_energy_V2(PINPT, PKPTS, E, V2, PGEOM, PRPLT, ne_found, emin_band, emax_band, &
                                             ispinor, nband, nspin, ispin, nbasis, &
                                             flag_erange, flag_vector, flag_sparse, &
                                             init_erange, fina_erange) !, c_mode_print)
     endif
   endif


   if_main call time_check(time2,time1)
   write(message,*)''  ; write_msg
   write(message,'(3A,F10.4,A)')'  END REPLOT: ',trim(cjob),' EVALUATION : ',time2, ' (sec)'  ; write_msg
   write(message,*)''  ; write_msg

   return
endsubroutine

subroutine get_dos_ldos(PINPT, PGEOM, PRPLT, ne_found, E, V2, nspin, nbasis, nband, nkp, e_range, nediv, sigma, &
                        flag_replot_dos, flag_replot_ldos, flag_replot_sldos, flag_replot_didv, flag_vector)
   use parameters, only : incar, poscar, replot
   use mpi_setup
   use memory
   use print_io
   use do_math, only : fgauss
   implicit none
   type(incar  )            :: PINPT
   type(poscar )            :: PGEOM
   type(replot )            :: PRPLT
   integer*4                   ik, ie, ia, is, i
   integer*4                   iadd, inc, iaddk, inck
   integer*4                   iatom, ii, fi
   integer*4                   isum
   integer*4                   nbasis, nband, nspin, nkp
   integer*4                   nediv
   integer*4                   mpierr
   integer*4                   ne_found(nspin, nkp)
   real*8                      e_range(nediv), de, sigma
   real*8                      E(nband*nspin,nkp)
   real*8                      V2(nbasis,nband*nspin,nkp) ! not wavevector, 
   real*8                      myV(nbasis,nband*nspin)
   real*8                      replot_dos_tot(nspin,nediv)
   real*8                      replot_dos_ntot(nspin,0:nediv)
   real*8, allocatable      :: replot_ldos_tot(:,:,:,:,:) ! (n_orbital(iatom), maxval(ldos_natom), nspin, nediv,nldos_sum)
   real*8, allocatable      :: replot_sldos_sum(:,:,:)
   logical                     flag_replot_dos, flag_replot_ldos, flag_replot_sldos, flag_replot_didv
   logical                     flag_vector
   logical                     flag_sum
   real*8                      dos_
   real*8                      a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
   real*8                      b2xb3(3),bzvol,dkv
   integer*4                   ldos_natom, ldos_atom(maxval(PRPLT%replot_ldos_natom(1:PRPLT%replot_nldos_sum)))
   integer*4                   n_orbital(PGEOM%n_atom)
   real*8                      coord_cart(3,PGEOM%n_atom)
#ifdef MPI                  
   integer*4                   ourjob(nprocs)
   integer*4                   ourjob_disp(0:nprocs-1)
#else                        
   integer*4                   ourjob(1)
   integer*4                   ourjob_disp(0)
#endif
   logical                     flag_para_band
   integer*4                   ie_init, ie_fina, i_init, i_fina

   flag_para_band = .false. ! default

   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   call get_reci(b1,b2,b3, a1,a2,a3)
   call vcross(b2xb3,b2,b3)
   bzvol=dot_product(b1,b2xb3)
   dkv  = bzvol / dble(nkp)
   if(nediv .ne. 1) then
     de   = e_range(2)-e_range(1)
   elseif(nediv .eq. 1) then
     de   = 1d0
   endif

   n_orbital  = PGEOM%n_orbital
   iadd       = 10 ; iaddk = 8  ; inck = 1

   allocate(PRPLT%replot_dos_tot(nspin, nediv))
   PRPLT%replot_dos_tot  = 0d0
         replot_dos_tot  = 0d0
   allocate(PRPLT%replot_dos_ntot(nspin,0:nediv))
   PRPLT%replot_dos_ntot = 0d0
         replot_dos_ntot = 0d0 
   write(message,*)''  ; write_msg
   if_main call report_memory(int8(size(replot_dos_tot))*2*nprocs, 8, 'DOS(total)    ')

   if(flag_replot_ldos .and. flag_vector) then
     allocate(PRPLT%replot_ldos_tot(PGEOM%max_orb, maxval(PRPLT%replot_ldos_natom(:)), nspin, nediv, PRPLT%replot_nldos_sum))
     allocate(      replot_ldos_tot(PGEOM%max_orb, maxval(PRPLT%replot_ldos_natom(:)), nspin, nediv, PRPLT%replot_nldos_sum))
     PRPLT%replot_ldos_tot = 0d0
           replot_ldos_tot = 0d0
     if_main call report_memory(int8(size(replot_ldos_tot))*2*nprocs, 8, 'LDOS(total)   ')
   endif

   if((flag_replot_sldos .or. flag_replot_didv) .and. flag_vector) then
     allocate(PRPLT%replot_sldos_sum(PGEOM%n_atom, nspin, nediv))
     allocate(      replot_sldos_sum(PGEOM%n_atom, nspin, nediv))
     PRPLT%replot_sldos_sum = 0d0
           replot_sldos_sum = 0d0
     if_main call report_memory(int8(size(replot_sldos_sum))*2*nprocs, 8, 'SLDOS(total)  ')
   endif

   if((flag_replot_ldos .or. flag_replot_sldos .or. flag_replot_didv) .and. flag_vector) then
     if_main call report_memory(int8(size(myV))*nprocs + size(V2), 8, 'Wave vector2  ')
   endif

   if(nband .le. nediv) then
     flag_para_band = .false.
   elseif(nband .gt. nediv) then
     flag_para_band = .true.
   endif

   write(message,*)''  ; write_msg
   do ik = 1, nkp
     if(nkp .lt. 10) then
       write(message,'(A,I0,A,I0)')  '     STAT KP: ', ik,'/',nkp  ; write_msg
     else
       if( ik/real(nkp)*100d0 .ge. real(iaddk*inck) ) then
         write(message,'(A,F10.3,A)')'     STAT KP: ', ik/real(nkp)*100d0, ' %'  ; write_msg
         inck = inck + 1
       endif
     endif

     if((flag_replot_ldos .or. flag_replot_sldos .or. flag_replot_didv) .and. flag_vector) then
#ifdef MPI
       if_main myV = V2(:,:,ik)
       call MPI_BCAST(myV, size(myV), MPI_REAL8, 0, mpi_comm_earth, mpierr)
#else
       myV = V2(:,:,ik)
#endif
     endif

     inc = 1

     if(.not. flag_para_band) then
       call mpi_job_distribution_chain(nediv, ourjob, ourjob_disp)
       ie_init = sum(ourjob(1:myid)) + 1 ; ie_fina  = sum(ourjob(1:myid+1))
     elseif(flag_para_band) then
       ie_init = 1                       ; ie_fina  = nediv
     endif

 eig:do ie = ie_init, ie_fina

       ! This routine should be checked later on, from here to
       if(nkp .lt. 10) then
         if(.not. flag_para_band) then
           if( (ie-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0 .ge. real(iadd*inc) ) then
             write(message,'(A,F10.3,A)')'    STAT EN: ', (ie-sum(ourjob(1:myid)))/real(ourjob(myid+1))*100d0, ' %'  ; write_msg
             inc = inc + 1
           endif
         elseif(flag_para_band) then
           if( (ie-(ie_init - 1)      )/real(ie_fina       )*100d0 .ge. real(iadd*inc) ) then
             write(message,'(A,F10.3,A)')'    STAT EN: ', (ie- (ie_init - 1)     )/real( ie_fina      )*100d0, ' %'  ; write_msg
             inc = inc + 1
           endif
         endif
       endif
      ! here! I'm not quite sure that "STAT EN:" in the case of "flag_para_band" reports correct status.
      ! Please check again.. but the main result does not altered whether above "if ~~~ endif" routine is correct or not. KHJ
    

       do is = 1, nspin

         if(.not. flag_para_band) then
           i_init  = 1                       ; i_fina   = ne_found(is,ik)
         elseif(flag_para_band) then
           call mpi_job_distribution_chain(ne_found(is,ik), ourjob, ourjob_disp)
           i_init  = sum(ourjob(1:myid)) + 1 ; i_fina   = sum(ourjob(1:myid+1))
         endif

         do i = i_init, i_fina

           dos_ = fgauss(sigma, e_range(ie) - E(i+nband*(is-1),ik)) / real(nkp)

           replot_dos_tot(is, ie) = replot_dos_tot(is, ie) + dos_

           if(flag_replot_ldos .and. flag_vector) then
             do isum = 1, PRPLT%replot_nldos_sum
               ldos_natom = PRPLT%replot_ldos_natom(isum)
               ldos_atom  = PRPLT%replot_ldos_atom(1:ldos_natom,isum)

               do ia = 1, ldos_natom
                 iatom = ldos_atom(ia)
                 ii = sum(n_orbital(1:iatom)) - n_orbital(iatom) + 1
                 fi = ii + n_orbital(iatom)-1
                 replot_ldos_tot(1:n_orbital(iatom),ia,is,ie,isum) = replot_ldos_tot(1:n_orbital(iatom),ia,is,ie,isum) + &
                                                                     dos_ * myV(ii:fi,i+nband*(is-1))
               enddo
             enddo
           endif

           if((flag_replot_sldos .or. flag_replot_didv) .and. flag_vector) then
             do iatom = 1, PGEOM%n_atom
               ii = sum(n_orbital(1:iatom)) - n_orbital(iatom) + 1
               fi = ii + n_orbital(iatom)-1
               replot_sldos_sum(iatom,is,ie) = replot_sldos_sum(iatom,is,ie) + dos_ * sum(myV(ii:fi,i+nband*(is-1)))
             enddo
           endif


         enddo
       enddo
     enddo eig
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
     write(message,'(A)')'    DONE!'  ; write_msg
   endif

   if(flag_vector) then
     if(flag_replot_ldos) then
       call MPI_ALLREDUCE(replot_ldos_tot, PRPLT%replot_ldos_tot, size(replot_ldos_tot), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
       call print_replot_ldos(PRPLT, PINPT, PGEOM)
     endif
     if(flag_replot_sldos .or. flag_replot_didv) then
       call MPI_REDUCE(replot_sldos_sum, PRPLT%replot_sldos_sum, size(replot_sldos_sum), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
       if(flag_replot_sldos) then
         flag_sum = .true.
         if_main call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart, flag_sum)
       endif
       if(flag_replot_didv) then
         flag_sum = .false.
         call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart, flag_sum)
       endif
       call MPI_BCAST(coord_cart, size(coord_cart), MPI_REAL8, 0, mpi_comm_earth, mpierr)
       call print_bond(PRPLT, PINPT, PGEOM, coord_cart)
     endif
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
     write(message,'(A)')'    write DOS ... Done!' ; write_msg
   endif

   if(flag_vector) then
     if(flag_replot_ldos) then
       PRPLT%replot_ldos_tot = replot_ldos_tot
       call print_replot_ldos(PRPLT, PINPT, PGEOM)
     endif
     if(flag_replot_sldos .or. flag_replot_didv) then
       PRPLT%replot_sldos_sum = replot_sldos_sum
       if(flag_replot_sldos) then
         flag_sum = .true.
         call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart, flag_sum)
       endif
       if(flag_replot_didv) then
         flag_sum = .false.
         call print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart, flag_sum)
       endif
       call print_bond(PRPLT, PINPT, PGEOM, coord_cart)
     endif
   endif

#endif

   if(allocated(replot_ldos_tot))  deallocate(replot_ldos_tot)
   if(allocated(replot_sldos_sum)) deallocate(replot_sldos_sum)
   return
endsubroutine  

subroutine print_replot_dos(PRPLT, PINPT)
  use parameters, only : replot, pid_dos, incar
  use print_io
  use mpi_setup
  implicit none
  type(replot) :: PRPLT 
  type(incar)  :: PINPT
  integer*4       i, ie
  real*8          e_range(PRPLT%replot_dos_nediv)
  character*40    filenm
  logical         flag_collinear
! filenm = 'DOS.replot.dat'
! filenm = 'DOS.replot'//trim(PINPT%title(PRPLT%mysystem))//'.dat'
  filenm = trim(PRPLT%replot_dos_fname)//trim(PINPT%title(PRPLT%mysystem))//'.dat'
  e_range = PRPLT%replot_dos_erange

  open(pid_dos, file=trim(filenm), status = 'unknown')
  write(message,'(A)')' ' ; write_msg
  write(message,'(3A)')'   WRITING.... density of states (DOS): ',trim(filenm),' ... ' ;write_msg

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
subroutine print_replot_energy(PINPT,PKPTS,E,V,PGEOM,PRPLT, ne_found, emin, emax, ispinor, nband, nspin, ispin, nbasis, &
                                    flag_erange, flag_vector, flag_sparse, init_erange, fina_erange) !, c_mode)
   use parameters, only: pid_energy, incar, poscar, kpoints, replot, zi
   use mpi_setup
   use print_io
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   type(kpoints):: PKPTS
   integer*4       ie, is, ik, im, ib
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
   character*2     c_mode, c_mode_(PRPLT%replot_nband)

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

   c_mode_ = PRPLT%replot_axis_print_mag

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   call get_e_range(init_e, fina_e, nbasis, .false., ispinor, flag_erange, init_erange, fina_erange)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

nb:do ib = 1, PRPLT%replot_nband
spin:do is = 1, nspin

       c_mode = c_mode_(ib)
       if(flag_vector .and. c_mode .eq. 'no') then
         flag_print_orbital = .false.
       else
         flag_print_orbital = flag_vector
       endif

       write(fname_header,'(4A)')'band_structure_TBA',trim(PINPT%title(PRPLT%mysystem)),'.replot_',trim(c_mode)
       
       call get_fname(fname_header, fname, is, flag_collinear, flag_noncollinear)
       open(pid_energy, file=trim(fname), status = 'unknown')
       write(message,'(A)')' ' ; write_msg
       write(message,'(5A)')'   WRITING.... band structure with ',trim(c_mode),' : ', trim(fname),' ... ' ; write_msg

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
                   c_up = V(im+nbasis*(is-1),ie+PGEOM%nband*(is-1),ik)
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
       write(message,'(A)')'    DONE!' ; write_msg
     enddo spin
   enddo nb

return
endsubroutine

subroutine print_replot_energy_proj(PINPT, PKPTS,E,V2,PGEOM,PRPLT, ne_found, emin, emax, ispinor, nband, nspin, ispin, nbasis, &
                                    flag_erange, flag_vector, flag_sparse, init_erange, fina_erange, c_mode)
   use parameters, only: pid_energy, incar, poscar, kpoints, replot, zi
   use mpi_setup
   use print_io
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
   real*8          V2(nbasis,nband*nspin,PKPTS%nkpoint) ! rho_nk = <psi_nk|psi_nk> is stored, not complex wavefunction
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
           write(fname_header_sum,'(2A,I0)')'band_structure_TBA_atom',trim(PINPT%title(PRPLT%mysystem)),'.sum',isum
  
           call get_fname(fname_header_sum, fname_sum, is, flag_collinear, flag_noncollinear)
           write(message,'(A)')' ' ; write_msg
           write(message,'(3A)')'   WRITING.... projected band structure: ',trim(fname_sum),' ... ' ; write_msg
           open(pid_energy+100, file = trim(fname_sum), status = 'unknown')
           if(flag_sparse) then
             write(pid_energy+100, '(A,2(F10.4,A))')'# The EWINDOW mode: energy window [EMIN:EMAX]=[', &
                                                 emin,':', emax,']'
             do ik = 1, nkpoint
               write(pid_energy+100, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
             enddo
           endif
           write(pid_energy+100, '(A, *(I0,1x))') '#  ATOM_INDEX to be sum up: ', proj_atom(1:proj_natom)
         endif

     atom:do iatom = 1, proj_natom
         ia = proj_atom(iatom)
         imatrix = sum( PGEOM%n_orbital(1:ia) ) - PGEOM%n_orbital(ia) + 1
         write(fname_header,'(2A,I0)')'band_structure_TBA_atom',trim(PINPT%title(PRPLT%mysystem)),ia
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
                     c_up = V2(im,ie,ik)
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') c_up
                     c_tot = c_tot + c_up
                   elseif(ispinor .eq. 1) then
!                    c_up = V(im+PGEOM%neig*(is-1),ie+nband*(is-1),ik)
                     c_up = V2(im,ie+nband*(is-1),ik)
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
       write(message,'(A)')'    DONE!' ; write_msg
     enddo spin

   enddo


return
endsubroutine
subroutine print_replot_ldos(PRPLT, PINPT, PGEOM)
   use parameters, only: replot, pid_ldos, incar, poscar
   use mpi_setup
   use print_io
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   integer*4       mpierr
   integer*4       im, ia, iaa, ie, iatom
   integer*4       isum
   integer*4       my_pid_ldos
   integer*4       ldos_natom
   integer*4       ldos_atom(maxval(PRPLT%replot_ldos_natom(1:PRPLT%replot_nldos_sum)))
   real*8          e_range(PRPLT%replot_dos_nediv), de
   character*40    filenm
   character*40    filenm_sum
   real*8          ldos_sum(PINPT%nspin,PRPLT%replot_dos_nediv,PRPLT%replot_nldos_sum)
   real*8          ldos_sum_ntot(PINPT%nspin)
#ifdef MPI
   real*8          ldos_sum_(PINPT%nspin,PRPLT%replot_dos_nediv,PRPLT%replot_nldos_sum)
   integer*4       ourjob(nprocs)
   integer*4       ourjob_disp(0:nprocs-1)
!  call mpi_job_distribution_chain(PRPLT%replot_ldos_natom, ourjob, ourjob_disp)
#else
   integer*4       ourjob(1)
   integer*4       ourjob_disp(0)
!  call mpi_job_distribution_chain(PRPLT%replot_ldos_natom, ourjob, ourjob_disp)
#endif

   e_range = PRPLT%replot_dos_erange
   if(PRPLT%replot_dos_nediv .gt. 1) then
     de      = e_range(2) - e_range(1)
   elseif(PRPLT%replot_dos_nediv .eq. 1) then
     de      = 1d0
   endif
   ldos_sum= 0d0

 sum1:do isum = 1, PRPLT%replot_nldos_sum

        ldos_natom = PRPLT%replot_ldos_natom(isum)
        ldos_atom  = PRPLT%replot_ldos_atom(1:ldos_natom, isum)

        write(filenm_sum,'(A,I0,2A)')'LDOS.replot.sum',isum,trim(PINPT%title(PRPLT%mysystem)),'.dat'
        if_main open(pid_ldos, file=trim(filenm_sum), status = 'unknown') 
        if_main write(pid_ldos,'(A, *(I0,1x))')'# ATOM_INDEX to be sum up: ', ldos_atom(1:ldos_natom)
        if_main write(pid_ldos,'(A,I8,A,F16.8)')'# NDIV = ',PRPLT%replot_dos_nediv,' dE=',e_range(2)-e_range(1)
        if(.not. PINPT%flag_collinear) then ! orbital information would not be stored since each atom would have different orbital sets in general
          if_main write(pid_ldos,'(A)',ADVANCE='YES')          '#    energy (ev)        dos_total           nelect'
        elseif(PINPT%flag_collinear) then
          if_main write(pid_ldos,'(A)',ADVANCE='YES')          '#    energy (ev)     dos_total-up     dos_total-dn           nelect-up        nelect-dn'
        endif

        write(message,'(A)')' '  ; write_msg
        write(message,'(3A)')'   WRITING.... local density of states (LDOS): ',trim(filenm_sum),' ... '  ; write_msg
        
        call mpi_job_distribution_chain(ldos_natom, ourjob, ourjob_disp)
   atom:do ia = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
          iatom = ldos_atom(ia)
          my_pid_ldos = pid_ldos + myid + 1
          write(filenm,'(A,I0,2A)') 'LDOS.replot.atom_',iatom,trim(PINPT%title(PRPLT%mysystem)),'.dat'
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
      
       en:do ie = 1, PRPLT%replot_dos_nediv
            if(.not.PINPT%flag_collinear) then
              ldos_sum(1,ie,isum) = ldos_sum(1,ie,isum) + sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie,isum))
              write(my_pid_ldos,'(F16.8,1x,  F16.8    , *(F16.8,1x))')e_range(ie), sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie,isum)), &
                                                                                       PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie,isum)
            elseif(PINPT%flag_collinear) then
              ldos_sum(1,ie,isum) = ldos_sum(1,ie,isum) + sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie,isum))
              ldos_sum(2,ie,isum) = ldos_sum(2,ie,isum) + sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,2,ie,isum))
              write(my_pid_ldos,'(F16.8,1x,2(F16.8,1x), *(F16.8,1x))')e_range(ie), sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie,isum)), &
                                                                                   sum(PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,2,ie,isum)), &
                                                                                       PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,1,ie,isum) , &
                                                                                       PRPLT%replot_ldos_tot(1:PGEOM%n_orbital(iatom),ia,2,ie,isum)
            endif
          enddo en
      
          close(my_pid_ldos)
      
        enddo atom
#ifdef MPI
          call MPI_REDUCE(ldos_sum, ldos_sum_, size(ldos_sum), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
          if_main ldos_sum = ldos_sum_
#endif          
        ldos_sum_ntot = 0d0
        do ie = 1, PRPLT%replot_dos_nediv
          if(.not.PINPT%flag_collinear) then
            ldos_sum_ntot=ldos_sum_ntot + ldos_sum(1, ie, isum) * de
            if_main write(pid_ldos,'(F16.8,1x, 2F16.8    )')e_range(ie), ldos_sum(1, ie, isum), ldos_sum_ntot
          elseif(PINPT%flag_collinear) then
            ldos_sum_ntot(:)=ldos_sum_ntot(:) + ldos_sum(:, ie, isum)
            if_main write(pid_ldos,'(F16.8,1x,4(F16.8,1x))')e_range(ie), ldos_sum(1, ie, isum), ldos_sum(2, ie, isum), ldos_sum_ntot(:)
          endif
        enddo

        if_main close(pid_ldos) 
        write(message,'(A)')'    DONE!'  ; write_msg
      enddo sum1

   return
endsubroutine

subroutine print_replot_sldos(PRPLT, PINPT, PGEOM, nkp, coord_cart, flag_sum)
   use parameters, only: replot, pid_ldos, incar, poscar
   use mpi_setup
   use print_io
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   integer*4       iatom, is, i
   integer*4       ie
   character*40    filenm
   integer*4       nx, ny, nz
   integer*4       ix, iy, iz
   integer*4       nkp
   real*8          sldos(PGEOM%n_atom,PINPT%nspin) ! integrated spatial LDOS
   real*8          e_range(PRPLT%replot_dos_nediv)
   real*8          coord(3), coord_cart(3,PGEOM%n_atom), T(3,3), de
   logical         flag_sum
   integer*4       mypid_ldos, mpierr
!#ifdef MPI
   integer*4       ourjob(nprocs)
   integer*4       ourjob_disp(0:nprocs-1)
!#else

!#endif
   e_range = PRPLT%replot_dos_erange

   nx = PRPLT%replot_sldos_cell(1)
   ny = PRPLT%replot_sldos_cell(2)
   nz = PRPLT%replot_sldos_cell(3)
   T(:,1) = PGEOM%a_latt(:,1)
   T(:,2) = PGEOM%a_latt(:,2)
   T(:,3) = PGEOM%a_latt(:,3)
   if(PRPLT%replot_dos_nediv .ne. 1) then
     de     = PRPLT%replot_dos_erange(2) - PRPLT%replot_dos_erange(1)
   elseif(PRPLT%replot_dos_nediv .eq. 1) then
     de     = 1d0
   endif

  ! coord info
   do iatom=1,PGEOM%n_atom
     coord(:) = PGEOM%a_coord(:,iatom) + PRPLT%r_origin(:)
     coord(:) = coord(:) - nint(coord(:))
     coord_cart(:,iatom) = coord(1)*PGEOM%a_latt(:,1)+coord(2)*PGEOM%a_latt(:,2)+coord(3)*PGEOM%a_latt(:,3)
   end do

   if(flag_sum) then
     do is = 1, PINPT%nspin
       do iatom = 1, PGEOM%n_atom
         sldos(iatom,is) = sum(PRPLT%replot_sldos_sum(iatom, is, :)) * de
       enddo
     enddo
   elseif(.not. flag_sum) then
#ifdef MPI
     call MPI_BCAST(PRPLT%replot_sldos_sum, size(PRPLT%replot_sldos_sum), MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

!    do is = 1, PINPT%nspin
!      do iatom = 1, PGEOM%n_atom
!        
!      enddo
!    enddo
   endif

   if(flag_sum) then
     filenm = trim(PRPLT%replot_sldos_fname)//trim(PINPT%title(PRPLT%mysystem))//'.dat'
     if_main open(pid_ldos, file=trim(filenm), status = 'unknown')
     write(message,'(A)')' '  ; write_msg
     write(message,'(3A)')'   WRITING.... spatial local density of states (SLDOS) with energy window integrated: ',trim(filenm),' ... '  ; write_msg

     if_main write(pid_ldos,'(A,2(F16.8,A))')'# Integrated spatial local density of states : EWINDOW = [',PRPLT%replot_dos_emin,':',PRPLT%replot_dos_emax,']'
     if_main write(pid_ldos,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom
    !if_main write(pid_ldos,'(A,F16.8)'     )'# dE(eV)= ',de
     if_main write(pid_ldos,'(2(A,F16.8))'  )'# dE(eV)= ',de, ' ,SMEARING(eV)= ', PRPLT%replot_dos_smearing

     if(PINPT%flag_collinear) then
       if_main write(pid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)           SLDOS(UP)        SLDOS(DN)        ATOM_SPECIES'
     else
       if_main write(pid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)             SLDOS        ATOM_SPECIES'
     endif

     do iatom = 1, PGEOM%n_atom
       do ix = 1, nx
         do iy = 1, ny
           do iz = 1, nz
             coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                             T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                             T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
             if(PINPT%nspin .eq. 2) then
               if_main write(pid_ldos,'(3(F16.8,1x),2(F16.8,1x),I12)')coord, sldos(iatom, :), PGEOM%spec(iatom)
             else
               if_main write(pid_ldos,'(3(F16.8,1x), (F16.8,1x),I12)')coord, sldos(iatom, :), PGEOM%spec(iatom)
             endif
           enddo
         enddo
       enddo
     enddo
     write(message,'(A)')'    write SLDOS ... Done!'  ; write_msg
     if_main close(pid_ldos)
   elseif(.not. flag_sum) then
     write(message,'(A)')' '  ; write_msg
     write(message,'(4A)')'   WRITING.... spatial local density of states (SLDOS .or. dI/dV) within energy window : ./didv/',&
                          trim(PRPLT%replot_didv_fname),trim(PINPT%title(PRPLT%mysystem)),'.??? ...'  ; write_msg
     call system('mkdir -p ./didv')

     call mpi_job_distribution_chain(PRPLT%replot_dos_nediv, ourjob, ourjob_disp)

     !do ie = 1, PRPLT%replot_dos_nediv
     do ie = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
       filenm = trim(PRPLT%replot_didv_fname)//trim(PINPT%title(PRPLT%mysystem))
       mypid_ldos = pid_ldos + myid
       write(filenm,'(2A,I0)')trim(PRPLT%replot_didv_fname),'.',ie
       open(mypid_ldos, file='./didv/'//trim(filenm), status = 'unknown')

       write(mypid_ldos,'(A,2(F16.8,A))')'# Spatial local density of states : EWINDOW = [',PRPLT%replot_dos_emin,':',PRPLT%replot_dos_emax,']'
       write(mypid_ldos,'(A,I0)'        )'# NATOM = ',PGEOM%n_atom
       write(mypid_ldos,'(2(A,F16.8))'  )'# dE(eV)= ',de, ' ,SMEARING(eV)= ', PRPLT%replot_dos_smearing

       if(PINPT%flag_collinear) then
         write(mypid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)        ENERGY(eV)           SLDOS(UP)        SLDOS(DN)        ATOM_SPECIES      ATOM_NUMBER    SITE_INDEX'
       else
         write(mypid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)        ENERGY(eV)             SLDOS        ATOM_SPECIES      ATOM_NUMBER    SITE_INDEX'
       endif

       write(mypid_ldos,'(A,F16.8)')'# ENERGY (eV): ', e_range(ie)
       do iatom = 1, PGEOM%n_atom
         do ix = 1, nx
           do iy = 1, ny
             do iz = 1, nz
               coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                               T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                               T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
               if(PINPT%nspin .eq. 2) then
                 write(mypid_ldos,'(4(F16.8,1x),2(F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), PRPLT%replot_sldos_sum(iatom, :, ie), &
                                                                                                    PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
               else
                 write(mypid_ldos,'(4(F16.8,1x), (F16.8,1x),I12,4x,I12,12x,2A)')coord, e_range(ie), PRPLT%replot_sldos_sum(iatom, :, ie), &
                                                                                                    PGEOM%spec(iatom), iatom, '#', trim(PGEOM%site_cindex(iatom))
               endif
             enddo
           enddo
         enddo
       enddo
      !write(pid_ldos,'(A)')' '
      !write(pid_ldos,'(A)')' '
       close(mypid_ldos)
     enddo
#ifdef MPI
     call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif
     write(message,'(A)')'    write dI/dV ... Done!'  ; write_msg
   endif


   return
endsubroutine

subroutine print_bond(PRPLT, PINPT, PGEOM, coord_cart)
   use parameters, only: replot, pid_ldos, incar, poscar, onsite_tolerance
   use mpi_setup
   use print_io
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
subroutine check_band_header(PRPLT, PINPT, flag_vector, flag_wf, flag_sparse, nemax, emin, emax, nspin, &
                             flag_erange, init_erange, fina_erange, c_mode, flag_formatted)
   use parameters, only : pid_energy, replot, incar
   use mpi_setup
   use print_io
   implicit none
   type(incar)         :: PINPT
   type(replot)        :: PRPLT
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
   flag_wf     = .false.
   flag_vector = .false.
   flag_go     = .false.

   if(flag_formatted)      form_ = 'formatted'
   if(.not.flag_formatted) form_ = 'unformatted'

   pid = pid_energy
   if(nspin .eq. 2) then ! collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.dat'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.dat'
       inquire(file=trim(fnameu),exist=flag_exist)
       if(flag_exist) then 
         inquire(file=trim(fnamed),exist=flag_exist)
         if(.not. flag_exist) then
           write(message,'(3A)')'   !ERROR file: ', trim(fnamed),' does not exist!!'  ; write_msg
           write(message,'(A)') '    Exit program... : check_band_header'  ; write_msg
           kill_job
         endif
       elseif(.not. flag_exist) then
         write(message,'(3A)')'   !ERROR file: ', trim(fnameu),' does not exist!!'  ; write_msg
         write(message,'(A)') '    Exit program... : check_band_header'  ; write_msg
         kill_job
       endif
       
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.up.bin'
       fnamed = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dn.bin'
       inquire(file=trim(fnameu),exist=flag_exist)
       if(flag_exist) then
         inquire(file=trim(fnamed),exist=flag_exist)
         if(.not. flag_exist) then
           write(message,'(3A)')'   !ERROR file: ', trim(fnamed),' does not exist!!'  ; write_msg
           write(message,'(A)') '    Exit program... : check_band_header'  ; write_msg
           kill_job
         endif
       elseif(.not. flag_exist) then
         write(message,'(3A)')'   !ERROR file: ', trim(fnameu),' does not exist!!'  ; write_msg
         write(message,'(A)') '    Exit program... : check_band_header'  ; write_msg
         kill_job
       endif

     endif
   elseif(nspin .eq. 1) then ! nm or non-collinear
     if(flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.dat'
       inquire(file=trim(fnameu),exist=flag_exist)
     elseif(.not. flag_formatted) then
       fnameu = 'band_structure_TBA'//trim(PINPT%title(PRPLT%mysystem))//'.bin'
       inquire(file=trim(fnameu),exist=flag_exist)
     endif
     if(.not. flag_exist) then
       write(message,'(3A)')'   !ERROR file: ', trim(fnameu),' does not exist!!'  ; write_msg
       write(message,'(A)') '    Exit program... : check_band_header'  ; write_msg
       kill_job
     endif
   endif


   ! it is sufficient to check with is = 1 case only.
   if(nspin .eq. 1) then
     write(message,'(2A)')'   LOADING.... ', trim(fnameu)  ; write_msg
   elseif(nspin .eq. 2) then
     write(message,'(3A)')'   LOADING.... ', trim(fnameu), ' and ', trim(fnamed)  ; write_msg
   endif

   open(pid, file=trim(fnameu), form=trim(form_), status='old', iostat=i_continue)
   if(i_continue .ne. 0) then
     write(message,'(2A)')'  !!! ERROR IN READING ', trim(fnameu) ; write_msg
     write(message,'(A)')'      Exit program... : check_band_header ' ; write_msg
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
           write(message,'(A,2(F12.4,A),I0)')'   => EWINDOW [EMIN:EMAX]=[ ',emin,' : ',emax,' ], NE_MAX= ',nemax  ; write_msg
         endif
       elseif(index(inputline, 'ERANGE') .ge. 1) then
         flag_erange = .true.
         if(flag_erange) then
           call strip_off(trim(inputline), c_emin, 'ERANGE=[',':',1)
           call strip_off(trim(inputline), c_emax, ':',']',1)
           call str2int(c_emin, init_erange)
           call str2int(c_emax, fina_erange)
           nemax = fina_erange - init_erange + 1
           write(message,'(A,I0,A,I0,A,I0)')'   => ERANGE=[ ',init_erange,' : ',fina_erange,' ], NERANGE= ',nemax  ; write_msg
         endif
       elseif(index(inputline, 'LORBIT=') .ge. 1) then
         call strip_off(trim(inputline), c_mode, 'LORBIT=[',']',1)
         if(trim(c_mode) .eq. 'rh') then
           flag_wf = .false.
           flag_vector = .true.
           write(message,'(A,I0,A,I0,A,I0)')'   => C_MODE= rh ; <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j'  ; write_msg
         elseif(trim(c_mode) .eq. 'wf') then
           flag_wf = .true.
           flag_vector = .true.
           write(message,'(A,I0,A,I0,A,I0)')'   => C_MODE= wf ; wavefunction coefficent will be read.'  ; write_msg
         else
           flag_wf = .false.
           flag_vector = .false.
           write(message,'(A,I0,A,I0,A,I0)')'   => C_MODE= no ; only read eigenvalues.'  ; write_msg
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
       write(message,'(A,I0,A,I0,A,I0)')'   => ERANGE=[ ',init_erange,' : ',fina_erange,' ], NERANGE= ',nemax    ; write_msg
       if(nemax .ne. nband_) then 
         write(message,'(3A)')'         !ERROR in reading ',trim(fnameu),' : nemax .ne. nband_ ' ; write_msg
         write(message,'(A)') '          Exit program...' ; write_msg
         kill_job
       endif
     else
       read(pid) flag_erange ! .FALSE.
     endif
     if(flag_sparse) then
       read(pid) flag_sparse, emin, emax, nemax
       if(nemax .ne. nband_) then
         write(message,'(3A)')'         !ERROR in reading ',trim(fnameu),' : nemax .ne. nband_ ' ; write_msg
         write(message,'(A)') '          Exit program...' ; write_msg
         kill_job
       endif
       write(message,'(A,2(F12.4,A),I0)')'   => EWINDOW [EMIN:EMAX]=[ ',emin,' : ',emax,' ], NE_MAX= ',nemax  ; write_msg
     else
       read(pid) flag_sparse ! .FALSE.
     endif

     if(c_mode .eq. 'wf') then
       flag_wf = .true.  ! current version only support reading band_....bin file with wavefuction information is written.
       write(message,'(A)')'   => C_MODE= wf ; wavefunction coefficent will be read.'  ; write_msg
     elseif(c_mode .eq. 'rh') then
       flag_wf = .false.
       write(message,'(A)')'   => C_MODE= rh ; <phi_ij|psi_nk>, phi_ij : atomic orbital of atom i orbital j'  ; write_msg
     elseif(c_mode .eq. 'no') then
       flag_wf = .false. ! current version only support reading band_....bin file with wavefuction information is written.
       write(message,'(A)')'   => C_MODE= no ; only read eigenvalues.'  ; write_msg
     else
       write(message,'(5A)')'    !ERROR in reading ',trim(fnameu),' : cannot recognize C_MODE=',c_mode, ' in current version. Exit program...' ; write_msg
       kill_job
     endif
   endif

   close(pid)

   if(flag_vector .and. trim(c_mode) .eq. 'no') then 
     if(nspin .eq. 2) then
       write(message,'(4A)')'    !ERROR in reading ',trim(fnameu), ' and ', trim(fnamed)  ; write_msg
     elseif(nspin .eq. 1) then
       write(message,'(2A)')'    !ERROR in reading ',trim(fnameu)  ; write_msg
     endif
     write(message,'(A)') '     cannot read orbital information from avove files. Check again...'  ; write_msg
     write(message,'(A)') '     Exit program...'  ; write_msg
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
subroutine print_replot_energy_V2(PINPT, PKPTS,E,V2,PGEOM,PRPLT, ne_found, emin, emax, ispinor, nband, nspin, ispin, nbasis, &
                                    flag_erange, flag_vector, flag_sparse, init_erange, fina_erange) !, c_mode)
   use parameters, only: pid_energy, incar, poscar, kpoints, replot, zi
   use mpi_setup
   use print_io
   implicit none
   type(replot) :: PRPLT
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   type(kpoints):: PKPTS
   integer*4       ie, is, ik, im, ib
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
   real*8          V2(nbasis,nband*nspin,PKPTS%nkpoint) ! complex wavefunction (total)
   real*8          emin, emax
   character*80    fname_header
   character*80    fname
   character*6     kunit_
   character*28    kmode
   character*8     sigma
   character*2     c_mode, c_mode_(PRPLT%replot_nband)
   integer*4       ikmode
   real*8, allocatable :: kpoint_(:,:)

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

   c_mode_ = PRPLT%replot_axis_print_mag

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   call get_e_range(init_e, fina_e, nbasis, .false., ispinor, flag_erange, init_erange, fina_erange)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)


nb:do ib = 1, PRPLT%replot_nband

   if(.not. PRPLT%flag_replot_write_unformatted(ib)) then

  spin:do is = 1, nspin

         c_mode = c_mode_(ib)
         if(flag_vector .and. c_mode .eq. 'no') then
           flag_print_orbital = .false.
         else
           flag_print_orbital = flag_vector
         endif

         write(fname_header,'(4A)')'band_structure_TBA',trim(PINPT%title(PRPLT%mysystem)),'.replot_',trim(c_mode)
         call get_fname(fname_header, fname, is, flag_collinear, flag_noncollinear)
         open(pid_energy, file=trim(fname), status = 'unknown')
         write(message,'(A)')' ' ; write_msg
         write(message,'(5A)')'   WRITING.... band structure with ',trim(c_mode),' : ', trim(fname),' ... ' ; write_msg

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
             if(c_mode .eq. 'rh') then
               write(pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
               write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
             endif
          mm:do im=1,nbasis
               if(c_mode .eq. 'rh') then
                 write(pid_energy, '(I9)',ADVANCE='NO')im
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
                     if    (c_mode .eq. 'rh') then
                       write(pid_energy,'(*(F9.4))',ADVANCE='NO') V2(im,ie,ik) ! up + dn : total
                     endif
                   elseif(ispinor .eq. 1) then
                     if(c_mode .eq. 'rh') then
                       write(pid_energy,'(*(F9.4))',ADVANCE='NO') V2(im,ie+PGEOM%nband*(is-1),ik)
                     endif
                   endif
                 enddo basis
                 if(ispinor .eq. 2) then
                   if(c_mode .eq. 'rh') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='YES') V2(im,ie,ik) ! up + dn : total
                   endif
                 elseif(ispinor .eq. 1) then
                   if(c_mode .eq. 'rh') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='YES') V2(im,ie+PGEOM%nband*(is-1),ik)
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
         write(message,'(A)')'    DONE!' ; write_msg
       enddo spin

     elseif(PRPLT%flag_replot_write_unformatted(ib)) then
       c_mode = c_mode_(ib)
       if(flag_klinemode) then
         ikmode = 1
         allocate(kpoint_(1,nkpoint))
       elseif(flag_kgridmode) then
         ikmode = 3
         allocate(kpoint_(3,nkpoint))
       endif

       if(flag_vector .and. c_mode .eq. 'no') then
         flag_print_orbital = .false.
       else
         flag_print_orbital = flag_vector
       endif

 spinb:do is = 1, nspin

         write(fname_header,'(4A)')'band_structure_TBA',trim(PINPT%title(PRPLT%mysystem)),'.replot_',trim(c_mode)
         call get_fname_bin(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)

         ! write header
         open(pid_energy, file=trim(fname), form='unformatted', status='unknown')
         write(message,'(A)')' ' ; write_msg
         write(message,'(5A)')'   WRITING.... band structure with ',trim(c_mode),' : ', trim(fname),' ... ' ; write_msg
         write(pid_energy) ikmode, flag_print_orbital, PRPLT%flag_replot_print_single(ib), flag_erange, &
                           flag_sparse, nbasis, nkpoint, nband, &
                           ispin, nspin, ispinor, trim(c_mode)
         if(flag_erange) then
           write(pid_energy) flag_erange, init_erange, fina_erange
         else
           write(pid_energy) flag_erange ! .FALSE.
         endif
         if(flag_sparse) then
           write(pid_energy) flag_sparse, emin, emax, nband
         else
           write(pid_energy) flag_sparse ! .FALSE.
         endif

         ! write main wavefunction information
           if(flag_print_orbital) then
             if(ispinor .eq. 2) then
               if(c_mode .eq. 'wf') then
                 write(message,'(A)')'   ! ERROR: cannot write wavefunction information from "rh". Please check your REPLOT_BAND tag...' ; write_msg
                 stop
               elseif(c_mode .eq. 'rh') then
                 if(.not.PRPLT%flag_replot_print_single(ib)) then
                   do ik = 1, nkpoint
                     write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                     do ie = 1, ne_found(is,ik)
                       do im = 1, nbasis
                         write(pid_energy) V2(im,ie,ik)
                       enddo
                     enddo
                   enddo
                 elseif(PRPLT%flag_replot_print_single(ib)) then
                   do ik = 1, nkpoint
                     write(pid_energy) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
                     do ie = 1, ne_found(is,ik)
                       do im = 1, nbasis
                         write(pid_energy) real(V2(im,ie,ik),kind=4)
                       enddo
                     enddo
                   enddo
                 endif
               endif
             elseif(ispinor .eq. 1) then
               if(c_mode .eq. 'wf') then
                 write(message,'(A)')'   ! ERROR: cannot write wavefunction information from "rh". Please check your REPLOT_BAND tag...' ; write_msg
                 stop
               elseif(c_mode .eq. 'rh') then
                 if(.not.PRPLT%flag_replot_print_single(ib)) then
                   do ik = 1, nkpoint
                     write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                     do ie = 1+nband*(is-1),nband*(is-1)+ne_found(is,ik)
                       do im = 1+nbasis*(is-1),nbasis*is
                         write(pid_energy) V2(im,ie+PGEOM%nband*(is-1),ik)
                       enddo
                     enddo
                   enddo

                 elseif(PRPLT%flag_replot_print_single(ib)) then
                   do ik = 1, nkpoint
                     write(pid_energy) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
                     do ie = 1+nband*(is-1),nband*(is-1)+ne_found(is,ik)
                       do im = 1+nbasis*(is-1),nbasis*is
                         write(pid_energy) real(V2(im,ie+PGEOM%nband*(is-1),ik),kind=4)
                       enddo
                     enddo
                   enddo

                 endif
               endif
             endif
           else
             do ik = 1, nkpoint
               write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
             enddo

           endif

         close(pid_energy)
         write(message,'(A)')'    DONE!' ; write_msg
       enddo spinb
       

     endif

   enddo nb


   if(allocated(kpoint_)) deallocate(kpoint_)

return
endsubroutine
