#include "alias.inc"
subroutine print_band_structure(PKPTS, ETBA, EDFT, PGEOM, PINPT, PPRAM, PWGHT)
   use parameters, only : kpoints, energy, poscar, incar, params, weight
   use mpi_setup
   use print_io
   implicit none
   type(incar)   :: PINPT
   type(params)  :: PPRAM
   type(energy)  :: ETBA, EDFT
   type(poscar)  :: PGEOM
   type(weight)  :: PWGHT
   type(kpoints) :: PKPTS
   integer*4        mpierr

!  if(.not. PINPT%flag_distribute_nkp) then
     if( myid .ne. 0) return
!  endif

   if(PINPT%flag_print_energy_diff ) then

     if(PINPT%flag_get_band) call print_energy(PKPTS, ETBA%E, EDFT%E, ETBA%V, ETBA%SV, PGEOM%neig, PGEOM%init_erange, PGEOM%nband, &
                                               PINPT, PWGHT, PPRAM%flag_use_overlap, .TRUE.,'')

   elseif(.not. PINPT%flag_print_energy_diff ) then

     if(PINPT%flag_get_band) then
!      if(.not. PINPT%flag_distribute_nkp) then
         call print_energy(PKPTS, ETBA%E, ETBA%E, ETBA%V, ETBA%SV, PGEOM%neig, PGEOM%init_erange, PGEOM%nband, &
                           PINPT, PWGHT, PPRAM%flag_use_overlap, .TRUE., '')
         if(PINPT%flag_get_band_order) then
           call print_energy(PKPTS, ETBA%E_ORD, ETBA%E_ORD, ETBA%V_ORD, ETBA%SV_ORD, PGEOM%neig, PGEOM%init_erange, PGEOM%nband, &
                           PINPT, PWGHT, PPRAM%flag_use_overlap, .TRUE., '_ordered')
         endif
!      elseif(PINPT%flag_distribute_nkp) then
!       !call print_energy_distributed(ETBA%E, ETBA%E, ETBA%V, ETBA%SV, PGEOM%neig, PGEOM%init_erange, PGEOM%nband, &
!       !                              PINPT, PWGHT, PPRAM%flag_use_overlap, .TRUE., '')
!        write(message,'(A)')'    !WARN! The current version does not support PINPT%flag_distribute_nkp (LDISTRK) = .TRUE. ' ; write_msg
!        write(message,'(A)')'           for band structure storing. The subroutine for this "print_energy_distributed" ' ; write_msg
!        write(message,'(A)')'           will be implemented in the near future.' ; write_msg
!        write(message,'(A)')'           Exit program... 19. Mar. 2021. HJ Kim' ; write_msg
!        kill_job
!      endif
     endif

   endif

   if(PINPT%flag_print_proj ) then
     call print_energy_proj(PKPTS, ETBA%E, ETBA%V, ETBA%SV, PGEOM, PINPT, PPRAM)
   endif
   
    
   return
endsubroutine
subroutine print_energy_ensurf (kpoint, nkpoint,ie, nspin, E, V, SV, PGEOM, PINPT, fname_header, kunit, flag_use_overlap)
   use parameters, only : pid_energy, incar, poscar, zi
   implicit none
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   integer*4       is, ie, ik, im
   integer*4       nkpoint
   integer*4       nspin, nbasis
   real*8          kpoint(3,nkpoint)
   logical         flag_print_orbital
   real*8          E(nspin,nkpoint)
   complex*16      V(PGEOM%neig*PINPT%ispin,nspin,nkpoint)  
   complex*16      SV(PGEOM%neig*PINPT%ispin,nspin,nkpoint)  
   complex*16      c_up,  c_dn
   complex*16      sc_up, sc_dn
   character(*)    fname_header
   character*80    fname
   character*1     kunit
   character*6     kunit_
   character*8     sigma
   character*28    kmode
   logical         flag_use_overlap

   ! NOTE: this routine write eigenvalues for certain energy level "ie" (and their opposite spin state if collinear case)
   !       hence, variable size of E is (nspin,nkpoint)

   sigma='sigma_0 '
   nbasis = PGEOM%neig
   call get_kunit(kunit, kunit_)
   call get_plotmode(.false., .true., kunit_, kmode)
 spin:do is = 1, nspin
        call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
        open(pid_energy, file=trim(fname), status = 'unknown')

            write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', ie,' -th eigen'
            if(.not. PINPT%flag_print_orbital) then
              write(pid_energy,'(A)',ADVANCE='yes')''
            elseif(  PINPT%flag_print_orbital) then
              if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
              if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
              if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
              write(pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma= ',sigma
              write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
          mm: do im=1,nbasis
                write(pid_energy, '(I9)',ADVANCE='NO')im
                if(im .ge. 30 .and. im .lt. nbasis) then
                  write(pid_energy, '(A)',ADVANCE='NO')' ... '
                  exit mm
                endif
              enddo mm
              write(pid_energy,'(A)',ADVANCE='yes')''
            endif

         kp:do ik = 1, nkpoint
              write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(1+(is-1),ik)
              if(PINPT%flag_print_orbital) then
          basis:do im=1,nbasis-1
                  if(PINPT%ispinor .eq. 2) then
                    c_up = V(im,is,ik); c_dn = V(im + nbasis,is,ik)
                    if(flag_use_overlap) then
                      sc_up = SV(im,is,ik); sc_dn = SV(im + nbasis,is,ik)
                    else
                      sc_up = V(im,is,ik); sc_dn = V(im + nbasis,is,ik)
                    endif
                    if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*sc_up - conjg(c_dn)*sc_dn) ! up - dn : mz 
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*sc_up + conjg(c_up)*sc_dn) ! up*dn + dn*up : mx
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*sc_up - conjg(c_up)*sc_dn)*zi) ! (up*dn - dn*up)*i : my
                    else
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn) ! up + dn : total
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    c_up = V(im+PGEOM%neig*(is-1),is,ik)
                    if(flag_use_overlap) then
                      sc_up = SV(im+PGEOM%neig*(is-1),is,ik)
                    else
                      sc_up = V(im+PGEOM%neig*(is-1),is,ik)
                    endif
                    write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*sc_up)
                  endif
                enddo basis
                if(PINPT%ispinor .eq. 2) then
                  c_up = V(im,is,ik); c_dn = V(im + nbasis,is,ik)
                  if(flag_use_overlap) then
                    sc_up = SV(im,is,ik); sc_dn = SV(im + nbasis,is,ik)
                  else
                    sc_up =  V(im,is,ik); sc_dn =  V(im + nbasis,is,ik)
                  endif
                  if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                    write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*sc_up - conjg(c_dn)*sc_dn) ! up - dn : mz
                  elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                    write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_dn)*sc_up + conjg(c_up)*sc_dn) ! up*dn + dn*up : mx
                  elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                    write(pid_energy,'(*(F9.4))',ADVANCE='YES') real((conjg(c_dn)*sc_up - conjg(c_up)*sc_dn)*zi) ! (up*dn - dn*up)*i : my
                  else
                    write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn) ! up + dn : total
                  endif
                elseif(PINPT%ispinor .eq. 1) then
                  c_up = V(nbasis+PGEOM%neig*(is-1),is,ik)
                  if(flag_use_overlap) then
                    sc_up = SV(nbasis+PGEOM%neig*(is-1),is,ik)
                  else
                    sc_up = V(nbasis+PGEOM%neig*(is-1),is,ik)
                  endif
                  write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*sc_up)
                endif
              endif
              if(.not.PINPT%flag_print_orbital) write(pid_energy,*)''
            enddo kp
            write(pid_energy,*)''
            write(pid_energy,*)''

      close(pid_energy)

      enddo spin
   

return
endsubroutine
subroutine print_energy_eff( PKPTS, E, V, PGEOM, PINPT, neig, fname_header )
   use parameters, only : pid_energy, incar, poscar, kpoints, zi
   implicit none
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   type(kpoints):: PKPTS
   integer*4       neig
   integer*4       ie,is,ik,im
   integer*4       nspin, nbasis, nkpoint
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   real*8          E(neig*PINPT%ispin,PKPTS%nkpoint)
   complex*16      V(neig*PINPT%ispin,neig*PINPT%ispin,PKPTS%nkpoint)
   complex*16      c_up, c_dn
   character*80    fname_header
   character*80    fname
   character*6     kunit_
   character*28    kmode
   character*8     sigma

!  NOTE: this subroutine is not for SPARSE matrix
!        and probably need to be updated to accept ERANGE tag where init_erange is not 1

   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = PINPT%flag_print_orbital
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = neig
   nspin  = PINPT%nspin
   sigma='sigma_0 '

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

spin:do is = 1, nspin
     call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear) 
     open(pid_energy, file=trim(fname), status = 'unknown')
 eig:do ie =1, neig*PINPT%ispinor
       write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', ie,' -th eigen'
       if(.not. flag_print_orbital) then
         write(pid_energy,'(A)',ADVANCE='NO')''
       elseif(  flag_print_orbital) then
         if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
         if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
         if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
         write(pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
         write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
     mm: do im=1,nbasis
           write(pid_energy, '(I9)',ADVANCE='NO')im
           if(im .ge. 30 .and. im .lt. nbasis) then
             write(pid_energy, '(A)',ADVANCE='NO')' ... '
             exit mm
           endif
         enddo mm
         write(pid_energy,'(A)')''
       endif
    kp:do ik = 1, nkpoint
         if(flag_klinemode) then
           write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+neig*(is-1),ik)
         elseif(flag_kgridmode) then
           write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+neig*(is-1),ik)
         endif

         if(flag_print_orbital) then
     basis:do im = 1, nbasis-1
             if(PINPT%ispinor .eq. 2) then
               c_up = V(im,ie,ik); c_dn = V(im+nbasis,ie,ik)
               if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                 write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
               elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                 write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
               elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                 write(pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
               else
                 write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
               endif
             elseif(PINPT%ispinor .eq. 1) then
               c_up = V(im+neig*(is-1), ie+neig*(is-1),ik)
               write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*c_up) 
             endif
           enddo basis
           if(PINPT%ispinor .eq. 2) then
             c_up = V(im,ie,ik); c_dn = V(im+nbasis,ie,ik)
             if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
               write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
             elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
               write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
             elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
               write(pid_energy,'(*(F9.4))',ADVANCE='YES') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
             else
               write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
             endif
           elseif(PINPT%ispinor .eq. 1) then
             c_up = V(im+neig*(is-1), ie+neig*(is-1),ik)
             write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*c_up)
           endif
         endif
         if(.not. flag_print_orbital) write(pid_energy,*)''
       enddo kp
       write(pid_energy,*)''
       write(pid_energy,*)''
     enddo eig
     close(pid_energy)
   enddo spin

return
endsubroutine
subroutine print_energy_proj(PKPTS,E,V,SV,PGEOM,PINPT,PPRAM)
   use parameters, only: pid_energy, incar, poscar, kpoints, params, zi
   implicit none
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   type(kpoints):: PKPTS
   type(params ):: PPRAM
   integer*4       mysystem
   integer*4       ie,is,ik,im,ia
   integer*4       isum, iatom, iorb
   integer*4       nspin, nbasis
   integer*4       nkpoint
   integer*4       proj_natom, proj_atom(maxval(PINPT%proj_natom(1:PINPT%nproj_sum)))
   integer*4       init_e, fina_e
   integer*4       ne_found(PINPT%nspin, PKPTS%nkpoint)
   integer*4       imatrix
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital, flag_use_overlap
   logical         flag_proj_sum, flag_proj_atom
   real*8          E(PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      V(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      SV(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      c_up, c_dn, c_tot
   complex*16      sc_up, sc_dn
   complex*16      c_sum(PGEOM%nband,PKPTS%nkpoint)
   complex*16      c_sum_orb(PGEOM%max_orb,PGEOM%nband,PKPTS%nkpoint)
   complex*16      c_tot_orb(PGEOM%max_orb)
   character*80    fname_header, fname_header_sum
   character*80    fname, fname_sum
   character*6     kunit_
   character*28    kmode
   character*8     sigma

   mysystem       = PGEOM%mysystem
   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = PINPT%flag_get_orbital
   flag_proj_sum = PINPT%flag_print_proj_sum
   flag_proj_atom= PINPT%flag_print_proj_atom
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = PGEOM%neig
   nspin  = PINPT%nspin
   sigma='sigma_0 '
   flag_use_overlap = PPRAM%flag_use_overlap

   do isum = 1, PINPT%nproj_sum
     proj_natom = PINPT%proj_natom(isum)
     proj_atom  = PINPT%proj_atom(1:proj_natom,isum)
     
     if(PINPT%flag_sparse) then
       ne_found = PINPT%feast_ne
     else
       ne_found = PGEOM%nband
     endif

     call get_kunit(PKPTS%kunit, kunit_)
     call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
     init_e = PGEOM%init_erange ; fina_e = init_e + PGEOM%nband - 1

     if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)
   
  spin:do is = 1, nspin

         if(flag_proj_sum) then
           c_sum = 0d0
           c_sum_orb = 0d0
           write(fname_header_sum,'(3A,I0)')'band_structure_TBA_atom',trim(PINPT%title(mysystem)),'.sum',isum
           call get_fname(fname_header_sum, fname_sum, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
           open(pid_energy+100, file = trim(fname_sum), status = 'unknown')
           if(PINPT%flag_sparse) then
             write(pid_energy+100, '(A,2(F10.4,A),I0)')'# The EWINDOW mode: energy window [EMIN:EMAX]=[ ', &
                                                 PINPT%feast_emin,' : ', PINPT%feast_emax,' ], NE_MAX= ',PINPT%feast_nemax
             do ik = 1, nkpoint
               write(pid_energy+100, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
             enddo
           elseif(PINPT%flag_erange) then
             write(pid_energy+100, '(A,I0,A,I0,A)')'#   ERANGE=[ ',init_e,' : ',fina_e,' ]'
           endif
           write(pid_energy+100, '(A, *(I0,1x))') '#  ATOM_INDEX to be sum up: ', proj_atom(1:proj_natom)
         endif

     atom:do iatom = 1, proj_natom
         ia = proj_atom(iatom)
         imatrix = sum( PGEOM%n_orbital(1:ia) ) - PGEOM%n_orbital(ia) + 1
         if(flag_proj_atom) then
           write(fname_header,'(3A,I0)')'band_structure_TBA_atom',trim(PINPT%title(mysystem)),'.',ia 
           call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
           open(pid_energy, file = trim(fname), status = 'unknown')
           
           if(PINPT%flag_sparse) then
             write(pid_energy, '(A,2(F10.4,A),I0)')'# The EWINDOW mode: energy window [EMIN:EMAX]=[ ', &
                                                 PINPT%feast_emin,' : ', PINPT%feast_emax,'], NE_MAX= ',PINPT%feast_nemax
             do ik = 1, nkpoint
               write(pid_energy, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
             enddo
           elseif(PINPT%flag_erange) then
             write(pid_energy, '(A,I0,A,I0,A)')'#   ERANGE=[ ',init_e,' : ',fina_e,' ]'
           endif
         endif
     eig:do ie = 1, PGEOM%nband ! init_e, fina_e
           if(flag_proj_atom) then
             write(pid_energy,'(2A,I8,A,I8,3A)',ADVANCE='yes')kmode,'  energy(eV) :',init_e+ie-1,' -th eigen | ',ia, &
                                                        ' -th atom (spec= ',trim(PGEOM%c_spec(PGEOM%spec(ia))),' )'
           endif

           if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
           if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
           if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
           
           if(flag_proj_atom) then
             write(pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma     
             write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
             do im=imatrix, imatrix + PGEOM%n_orbital(ia) - 1
               write(pid_energy, '(I9)',ADVANCE='NO')im
             enddo
             write(pid_energy,'(A9)',ADVANCE='YES') ' tot'
           endif

           if(iatom .eq. proj_natom .and. flag_proj_sum) then
             write(pid_energy+100,'(2A,I8,A      )',ADVANCE='yes')kmode,'  energy(eV) :',init_e+ie-1,' -th eigen '
             if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
             if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
             if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '

             write(pid_energy+100, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
             write(pid_energy+100, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)     E(eV), '
             do im=1      , PGEOM%n_orbital(ia) ! note: PROJ_BAND should be summed up with those atoms having same number of basis and same orbitals
              !write(pid_energy+100, '(I9)',ADVANCE='NO')im
               write(pid_energy+100, '(1X,A8)',ADVANCE='NO')trim(PGEOM%c_orbital(im,ia))
             enddo
             write(pid_energy+100,'(A)',ADVANCE='YES') '  tot(atom_sum)'
           endif

        kp:do ik = 1, nkpoint
             if(flag_klinemode .and. flag_proj_atom) then
               if( ie .le. ne_found(is, ik) ) then
                 write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+PGEOM%nband*(is-1),ik)
               elseif( ie .gt. ne_found(is, ik)) then
                 write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
               endif
             elseif(flag_kgridmode .and. flag_proj_atom) then
               if( ie .le. ne_found(is, ik) ) then
                 write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+PGEOM%nband*(is-1),ik)
               elseif(ie .gt. ne_found(is, ik)) then
                 write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
               endif
             endif

             if(flag_proj_sum .and. iatom .eq. proj_natom) then
               if(flag_klinemode) then
                 if( ie .le. ne_found(is, ik) ) then
                   write(pid_energy+100,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+PGEOM%nband*(is-1),ik)
                 elseif( ie .gt. ne_found(is, ik)) then
                   write(pid_energy+100,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
                 endif
               elseif(flag_kgridmode) then
                 if( ie .le. ne_found(is, ik) ) then
                   write(pid_energy+100,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+PGEOM%nband*(is-1),ik)
                 elseif(ie .gt. ne_found(is, ik)) then
                   write(pid_energy+100,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
                 endif
               endif
             endif

             if( ie .le. ne_found(is, ik) ) then
               if(flag_print_orbital) then
                 c_tot = 0d0 !initialize
                 c_tot_orb = 0d0 ! initialize
                 iorb=0
           basis:do im=imatrix, imatrix+PGEOM%n_orbital(ia) - 1
                    iorb=iorb + 1
                   if(PINPT%ispinor .eq. 2) then
                     c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
                     if(flag_use_overlap) then
                       sc_up = SV(im,ie,ik); sc_dn = SV(im + nbasis,ie,ik)
                     else
                       sc_up =  V(im,ie,ik); sc_dn =  V(im + nbasis,ie,ik)
                     endif
                     if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                       c_tot_orb(iorb) = real( conjg(c_up)*sc_up - conjg(c_dn)*sc_dn) ! ! up - dn : mz
                     elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                       c_tot_orb(iorb) = real( conjg(c_dn)*sc_up + conjg(c_up)*sc_dn) ! ! up*dn + dn*up : mx
                     elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                       c_tot_orb(iorb) = real((conjg(c_dn)*sc_up - conjg(c_up)*sc_dn)*zi) ! (up*dn - dn*up)*i : my
                     else
                       c_tot_orb(iorb) = real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn) ! up + dn : total
                     endif
                   elseif(PINPT%ispinor .eq. 1) then
                     c_up = V(im+PGEOM%neig*(is-1),ie+PGEOM%nband*(is-1),ik)
                     if(flag_use_overlap) then
                       sc_up = SV(im+PGEOM%neig*(is-1),ie+PGEOM%nband*(is-1),ik)
                     else
                       sc_up =  V(im+PGEOM%neig*(is-1),ie+PGEOM%nband*(is-1),ik)
                     endif
                     c_tot_orb(iorb) = real(conjg(c_up)*sc_up)
                   endif

                   if(flag_proj_atom) then
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(c_tot_orb(iorb))
                   endif
                   c_tot = c_tot + real(c_tot_orb(iorb))

                 enddo basis
                 if(flag_proj_atom) then
                   write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(c_tot)
                 endif

                 if(flag_proj_sum) then 
                   c_sum(ie,ik) = c_sum(ie,ik) + c_tot
                   c_sum_orb(1:PGEOM%n_orbital(ia),ie,ik) = c_sum_orb(1:PGEOM%n_orbital(ia),ie,ik) + c_tot_orb(1:PGEOM%n_orbital(ia))
                 endif
                 if(flag_proj_sum .and. iatom .eq. proj_natom) then
                   write(pid_energy+100,'(*(F9.4))',ADVANCE='YES') real(c_sum_orb(1:PGEOM%n_orbital(ia),ie,ik)), real(c_sum(ie, ik))
                 endif

               endif

               if(.not.flag_print_orbital .and. flag_proj_atom) write(pid_energy,*)''
               if(.not.flag_print_orbital .and. flag_proj_sum) write(pid_energy+100,*)'' ! maybe do not need.. but who knows?
             elseif(ie .gt. ne_found(is, ik)) then
               if(flag_proj_atom) write(pid_energy,*)''
               if(iatom .eq. proj_natom .and. flag_proj_sum) write(pid_energy+100,*)'' 
             endif
           enddo kp

           if(flag_proj_atom) write(pid_energy,*)''
           if(flag_proj_atom) write(pid_energy,*)''
           if(iatom .eq. proj_natom .and. flag_proj_sum) write(pid_energy+100,*)'' 
           if(iatom .eq. proj_natom .and. flag_proj_sum) write(pid_energy+100,*)'' 

         enddo eig

         if(flag_proj_atom) close(pid_energy)
       enddo atom
       if(iatom-1 .eq. proj_natom .and. flag_proj_sum) close(pid_energy+100)
     enddo spin

   enddo
return
endsubroutine
subroutine print_energy_singlek( E, V, SV, neig, iband, nband, iik , kp, PINPT, flag_use_overlap, mysystem)
   use parameters, only : pid_energy, incar, poscar, zi
   use mpi_setup
   type(incar)  :: PINPT
   integer*4       mysystem
   integer*4       neig, iband, nband
   integer*4       iik, nkp
   integer*4       ie, is ,ik, im
   integer*4       nspin, nbasis
   integer*4       ikmode
   integer*4       init_e, fina_e
   integer*4       ne_found(PINPT%nspin)
   integer*4       my_pid_energy, mpierr
   integer*4       my_pid_energyS
   real*8          E(nband*PINPT%nspin)
   real*8          kp(3)
   complex*16      V(neig*PINPT%ispin,nband*PINPT%nspin)
   complex*16      SV(neig*PINPT%ispin,nband*PINPT%nspin) ! overlap matrix multiplied wavefunction
   complex*16      c_up, c_dn
   complex*16      sc_up, sc_dn ! for overlap matrix multiplied wavefunction
   character*80    fname_header
   character*80    fname
   character*80    fname_headerS, fnameS ! for overlap matrix multiplied wavefunction
   character*28    kmode
   character*6     kunit_
   character*20,external ::   int2str
   logical         flag_print_orbital, flag_use_overlap, flag_file_exist
   character*8     sigma

   my_pid_energy = pid_energy + myid
   my_pid_energyS= pid_energy + myid + nprocs
   fname_header = 'band_structure_TBA'//trim(PINPT%title(mysystem))//'.kp_'//trim(ADJUSTL(int2str(iik)))
   fname_headerS = 'band_structure_TBA_S'//trim(PINPT%title(mysystem))//'.kp_'//trim(ADJUSTL(int2str(iik)))

   flag_print_orbital = PINPT%flag_print_orbital
   nbasis = neig
   nspin  = PINPT%nspin
   sigma='sigma_0 '
   ikmode = 3
   flag_file_exist = .FALSE.

   if(PINPT%flag_sparse) then
     ne_found = PINPT%feast_ne(1:PINPT%nspin,iik)
   else
     ne_found = nband
   endif

   call get_kunit('A', kunit_) ! enforce as default
   call get_plotmode(.false., .true., kunit_, kmode)
   init_e = iband ; fina_e = init_e + nband - 1

   if(.not. PINPT%flag_write_unformatted) then
 spin:do is = 1, nspin
        call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
        inquire(file=trim(fname),exist=flag_file_exist)
        if(flag_file_exist) return

        open(my_pid_energy, file=trim(fname), status = 'unknown')
        if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
          call get_fname(fname_headerS, fnameS, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
          open(my_pid_energyS, file=trim(fnameS), status = 'unknown')
        endif
          if(trim(PINPT%axis_print_mag) .eq. 'rh') then
            write(my_pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), &
                  ' ] -> <phi_ij|psi_nk> ; i,j => orbital j in atom i;  n,k => band index (n) and kpoint index (k)'
          elseif(trim(PINPT%axis_print_mag) .eq. 'wf') then
            write(my_pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ] -> wavefunction coefficients '
            if(flag_use_overlap) then
              write(my_pid_energyS,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ] -> wavefunction coefficients multiplied with overlap matrix S(k)'
            endif
          elseif(trim(PINPT%axis_print_mag) .eq. 'no') then
            write(my_pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ]'
          elseif(trim(PINPT%axis_print_mag) .eq. 'mx' .or. &
                 trim(PINPT%axis_print_mag) .eq. 'my' .or. &
                 trim(PINPT%axis_print_mag) .eq. 'mz') then
            write(my_pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ] -> magnetization <sigma_i>'
          endif
          if(PINPT%flag_sparse) then
            write(my_pid_energy, '(A,2(F10.4,A),I0)')'# The EWINDOW mode: energy window [EMIN:EMAX]=[ ',PINPT%feast_emin, &
                              ' : ',PINPT%feast_emax,' ], NE_MAX= ',PINPT%feast_nemax
            write(my_pid_energy, '(A,I0,A,I0)')'#   NE_FOUND(ik=',iik,')= ',ne_found(is)
            if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
              write(my_pid_energyS, '(A,2(F10.4,A),I0)')'# The EWINDOW mode: energy window [EMIN:EMAX]=[ ',PINPT%feast_emin, &
                                ' : ',PINPT%feast_emax,' ], NE_MAX= ',PINPT%feast_nemax
              write(my_pid_energyS, '(A,I0,A,I0)')'#   NE_FOUND(ik=',iik,')= ',ne_found(is)
            endif
          elseif(PINPT%flag_erange) then
             write(my_pid_energy, '(A,I0,A,I0,A)')'#   ERANGE=[ ',init_e,' : ',fina_e,' ]'
             if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
               write(my_pid_energy, '(A,I0,A,I0,A)')'#   ERANGE=[ ',init_e,' : ',fina_e,' ]'
             endif
          endif

      eig:do ie =1, nband !init_e, fina_e
            write(my_pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', init_e + ie - 1,' -th eigen'
            if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
              write(my_pid_energyS, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', init_e + ie - 1,' -th eigen'
            endif

            if(.not. flag_print_orbital) then
              write(my_pid_energy,'(A)',ADVANCE='NO')''
            elseif(  flag_print_orbital) then
              if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
              if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
              if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
              if(PINPT%axis_print_mag .ne. 'wf') then
                if(PINPT%axis_print_mag .ne. 'wf') then
                  write(my_pid_energy, '(2A)',ADVANCE='YES') '# overlap S multiplied wavefunction coeff.: <ci|sigma*S|ci>,sigma=',sigma
                else
                  write(my_pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
                endif
                write(my_pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
              elseif(PINPT%axis_print_mag .eq. 'wf') then
                write(my_pid_energy, '(1A)',ADVANCE='YES') '# wavefunction coeff.:          |ci>               '
                write(my_pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
                if(flag_use_overlap) then
                  write(my_pid_energyS, '(1A)',ADVANCE='YES') '# wavefunction coeff.:          |ci>               '
                  write(my_pid_energyS, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
                endif
              endif
           mm:do im=1,nbasis
                if(PINPT%axis_print_mag .ne. 'wf') then
                  write(my_pid_energy, '(I9)',ADVANCE='NO')im
                elseif(PINPT%axis_print_mag .eq. 'wf') then
                  if(PINPT%ispinor .eq. 2) then
                    write(my_pid_energy, '(I38)',ADVANCE='NO')im
                    if(flag_use_overlap) then
                      write(my_pid_energyS, '(I38)',ADVANCE='NO')im
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    write(my_pid_energy, '(I19)',ADVANCE='NO')im
                    if(flag_use_overlap) then
                      write(my_pid_energyS, '(I19)',ADVANCE='NO')im
                    endif
                  endif
                endif
                if(im .ge. 30 .and. im .lt. nbasis) then
                  write(my_pid_energy, '(A)',ADVANCE='NO')' ... '
                  if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                    write(my_pid_energyS, '(A)',ADVANCE='NO')' ... '
                  endif
                  exit mm
                endif
              enddo mm
              write(my_pid_energy,'(A)')''
              if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                write(my_pid_energyS,'(A)')''
              endif
            endif

            if( ie .le. ne_found(is    ) ) then
              write(my_pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kp(:), E(ie+nband*(is-1))
              if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                write(my_pid_energyS,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kp(:), E(ie+nband*(is-1))
              endif
            elseif(ie .gt. ne_found(is    )) then
              write(my_pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kp(:)
              if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                write(my_pid_energyS,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kp(:)
              endif
            endif
              if( ie .le. ne_found(is    ) ) then
                if(flag_print_orbital) then
            basis:do im=1,nbasis-1
                    if(PINPT%ispinor .eq. 2) then
                      c_up = V(im,ie); c_dn = V(im + nbasis,ie)
                      if(flag_use_overlap) then
                        sc_up = SV(im,ie); sc_dn = SV(im + nbasis,ie)
                      else
                        sc_up =  V(im,ie); sc_dn =  V(im + nbasis,ie)
                      endif
                      if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                        write(my_pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*sc_up - conjg(c_dn)*sc_dn) ! up - dn : mz
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                        write(my_pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*sc_up + conjg(c_up)*sc_dn) ! up*dn + dn*up : mx
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                        write(my_pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*sc_up - conjg(c_up)*sc_dn)*zi) ! (up*dn - dn*up)*i : my
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                        write(my_pid_energy,'(2(F9.4,F9.4," "))',ADVANCE='NO') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                        if(flag_use_overlap) then
                          write(my_pid_energyS,'(2(F9.4,F9.4," "))',ADVANCE='NO') sc_up, sc_dn ! sc_up and sc_dn (real,imag) S multiplied wavefunction coefficient
                        endif
                      else
                        write(my_pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn) ! up + dn : total
                      endif
                    elseif(PINPT%ispinor .eq. 1) then
                      c_up = V(im+nbasis*(is-1),ie+nband*(is-1))
                      if(flag_use_overlap) then
                        sc_up = SV(im+nbasis*(is-1),ie+nband*(is-1))
                      else
                        sc_up =  V(im+nbasis*(is-1),ie+nband*(is-1))
                      endif
                      if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                        write(my_pid_energy,'(1(F9.4,F9.4," "))',ADVANCE='NO') c_up
                      else
                        write(my_pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*sc_up)
                      endif
                    endif
                  enddo basis
                  if(PINPT%ispinor .eq. 2) then
                    c_up = V(im,ie); c_dn = V(im + nbasis,ie)
                    if(flag_use_overlap) then
                      sc_up = SV(im,ie); sc_dn = SV(im + nbasis,ie)
                    else
                      sc_up =  V(im,ie); sc_dn =  V(im + nbasis,ie)
                    endif
                    if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                      write(my_pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*sc_up - conjg(c_dn)*sc_dn) ! up - dn : mz
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                      write(my_pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_dn)*sc_up + conjg(c_up)*sc_dn) ! up*dn + dn*up : mx
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                      write(my_pid_energy,'(*(F9.4))',ADVANCE='YES') real((conjg(c_dn)*sc_up - conjg(c_up)*sc_dn)*zi) ! (up*dn - dn*up)*i : my
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                      write(my_pid_energy,'(2(F9.4,F9.4," "))',ADVANCE='YES') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                      if(flag_use_overlap) then
                        write(my_pid_energyS,'(2(F9.4,F9.4," "))',ADVANCE='YES') sc_up, sc_dn ! c_up and c_dn (real,imag) S multiploed wavefunction coefficient
                      endif
                    else
                      write(my_pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn) ! up + dn : total
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    c_up = V(nbasis+nbasis*(is-1),ie+nband*(is-1))
                    if(flag_use_overlap) then
                      sc_up = SV(nbasis+nbasis*(is-1),ie+nband*(is-1))
                    else
                      sc_up =  V(nbasis+nbasis*(is-1),ie+nband*(is-1))
                    endif
                    if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                      write(my_pid_energy,'(1(F9.4,F9.4," "))',ADVANCE='YES') c_up
                      if(flag_use_overlap) then
                        write(my_pid_energyS,'(1(F9.4,F9.4," "))',ADVANCE='YES') sc_up
                      endif
                    else
                      write(my_pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*sc_up)
                    endif
                  endif
                endif
                if(.not.flag_print_orbital) then
                  write(my_pid_energy,*)''
                endif
              elseif(ie .gt. ne_found(is)) then
                write(my_pid_energy,*)''
                if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                  write(my_pid_energyS,*)''
                endif
              endif
            write(my_pid_energy,*)''
            write(my_pid_energy,*)''
            if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
              write(my_pid_energyS,*)''
              write(my_pid_energyS,*)''
            endif
          enddo eig

        close(my_pid_energy)
        if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
          close(my_pid_energyS)
        endif
      enddo spin

   elseif(PINPT%flag_write_unformatted) then
     spinb:do is = 1, nspin
       call get_fname_bin(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
       inquire(file=trim(fname),exist=flag_file_exist)
       if(flag_file_exist) return

       if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
         call get_fname_bin(fname_headerS, fnameS, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
       endif

       open(my_pid_energy, file=trim(fname), form='unformatted', status='unknown')
       write(my_pid_energy) ikmode, flag_print_orbital, PINPT%flag_print_single, PINPT%flag_erange, &
                         PINPT%flag_sparse, neig, 1, nband, &
                         PINPT%ispin, PINPT%nspin, PINPT%ispinor, PINPT%axis_print_mag
       if(PINPT%flag_erange) then
         write(my_pid_energy) PINPT%flag_erange, init_e , fina_e
       else
         write(my_pid_energy) PINPT%flag_erange ! .FALSE.
       endif
       if(PINPT%flag_sparse) then
         write(my_pid_energy) PINPT%flag_sparse, PINPT%feast_emin, PINPT%feast_emax, PINPT%feast_nemax
       else
         write(my_pid_energy) PINPT%flag_sparse ! .FALSE.
       endif

       if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
         open(my_pid_energyS, file=trim(fnameS), form='unformatted', status='unknown')
         write(my_pid_energyS) ikmode, flag_print_orbital, PINPT%flag_print_single, PINPT%flag_erange, &
                           PINPT%flag_sparse, neig, 1, nband, &
                           PINPT%ispin, PINPT%nspin, PINPT%ispinor, PINPT%axis_print_mag
         if(PINPT%flag_erange) then
           write(my_pid_energyS) PINPT%flag_erange, init_e , fina_e
         else
           write(my_pid_energyS) PINPT%flag_erange ! .FALSE.
         endif
         if(PINPT%flag_sparse) then
           write(my_pid_energyS) PINPT%flag_sparse, PINPT%feast_emin, PINPT%feast_emax, PINPT%feast_nemax
         else
           write(my_pid_energyS) PINPT%flag_sparse ! .FALSE.
         endif
       endif

     ! write main wavefunction information
       if(flag_print_orbital) then
         if(PINPT%ispinor .eq. 2) then
           if(PINPT%axis_print_mag .eq. 'wf') then

             if(.not.PINPT%flag_print_single) then
               write(my_pid_energy) ne_found(is), kp(:), (E(ie+nband*(is-1)),ie=1,ne_found(is))
               if(flag_use_overlap) then
                 write(my_pid_energyS) ne_found(is), kp(:), (E(ie+nband*(is-1)),ie=1,ne_found(is))
               endif
               do ie = 1, ne_found(is)
                 do im = 1, nbasis
                   write(my_pid_energy) V(im,ie), V(im+nbasis,ie)
                   if(flag_use_overlap) write(my_pid_energyS) SV(im,ie), SV(im+nbasis,ie)
                 enddo
               enddo
             elseif(PINPT%flag_print_single) then
                 write(my_pid_energy) ne_found(is), real(kp(:),kind=4), (real(E(ie+nband*(is-1)),kind=4),ie=1,ne_found(is))
                 do ie = 1, ne_found(is)
                   do im = 1, nbasis
                     if(.not.flag_use_overlap) then
                       write(my_pid_energy) cmplx((/V(im,ie), V(im+nbasis,ie)/),kind=4)
                     else
                       write(my_pid_energy) cmplx((/V(im,ie), SV(im+nbasis,ie)/),kind=4)
                     endif
                   enddo
                 enddo
             endif

           elseif(PINPT%axis_print_mag .eq. 'rh') then
             if(.not.PINPT%flag_print_single) then
                 write(my_pid_energy) ne_found(is), kp(:), (E(ie+nband*(is-1)),ie=1,ne_found(is))
                 do ie = 1, ne_found(is)
                   do im = 1, nbasis
                     if(.not.flag_use_overlap) then
                       write(my_pid_energy) real( conjg(V(im,ie))*V(im,ie) + conjg(V(im+nbasis,ie))*V(im+nbasis,ie) )
                     else
                       write(my_pid_energy) real( conjg(V(im,ie))*SV(im,ie) + conjg(V(im+nbasis,ie))*SV(im+nbasis,ie) )
                     endif
                   enddo
                 enddo

             elseif(PINPT%flag_print_single) then
                 write(my_pid_energy) ne_found(is), real(kp(:),kind=4), (real(E(ie+nband*(is-1)),kind=4),ie=1,ne_found(is))
                 do ie = 1, ne_found(is)
                   do im = 1, nbasis
                     if(.not. flag_use_overlap) then
                       write(my_pid_energy) real( conjg(V(im,ie))*V(im,ie)+conjg(V(im+nbasis,ie))*V(im+nbasis,ie), kind=4 )
                     else
                       write(my_pid_energy) real( conjg(V(im,ie))*SV(im,ie)+conjg(V(im+nbasis,ie))*SV(im+nbasis,ie), kind=4 )
                     endif
                   enddo
                 enddo

             endif
           endif

         elseif(PINPT%ispinor .eq. 1) then
           if(PINPT%axis_print_mag .eq. 'wf') then
             if(.not.PINPT%flag_print_single) then
                 write(my_pid_energy) ne_found(is), kp(:), (E(ie+nband*(is-1)),ie=1,ne_found(is))
                 if(flag_use_overlap) then
                   write(my_pid_energyS) ne_found(is), kp(:), (E(ie+nband*(is-1)),ie=1,ne_found(is))
                 endif
                 do ie = 1+nband*(is-1), nband*(is-1)+ne_found(is)
                   do im = 1+nbasis*(is-1),nbasis*is
                     write(my_pid_energy) V(im,ie)
                     if(flag_use_overlap) then
                       write(my_pid_energyS) SV(im,ie)
                     endif
                   enddo
                 enddo
             elseif(PINPT%flag_print_single) then
                 write(my_pid_energy) ne_found(is), real(kp(:),kind=4), (real(E(ie+nband*(is-1)),kind=4),ie=1,ne_found(is))
                 if(flag_use_overlap) then
                   write(my_pid_energyS) ne_found(is), real(kp(:),kind=4), (real(E(ie+nband*(is-1)),kind=4),ie=1,ne_found(is))
                 endif
                 do ie = 1+nband*(is-1), nband*(is-1)+ne_found(is)
                   do im = 1+nbasis*(is-1),nbasis*is
                     write(my_pid_energy) cmplx(V(im,ie), kind=4)
                     if(flag_use_overlap) then
                       write(my_pid_energyS) cmplx(SV(im,ie), kind=4)
                     endif
                   enddo
                 enddo
             endif

           elseif(PINPT%axis_print_mag .eq. 'rh') then
             if(.not.PINPT%flag_print_single) then
                 write(my_pid_energy) ne_found(is), kp(:), (E(ie+nband*(is-1)),ie=1,ne_found(is))
                 do ie = 1+nband*(is-1),nband*(is-1)+ne_found(is)
                   do im = 1+nbasis*(is-1),nbasis*is
                     if(.not. flag_use_overlap) then
                       write(my_pid_energy) real(conjg(V(im,ie))*V(im,ie))
                     else
                       write(my_pid_energy) real(conjg(V(im,ie))*SV(im,ie))
                     endif
                   enddo
                 enddo
             elseif(PINPT%flag_print_single) then
                 write(my_pid_energy) ne_found(is), real(kp(:),kind=4), (real(E(ie+nband*(is-1)),kind=4),ie=1,ne_found(is))
                 do ie = 1+nband*(is-1),nband*(is-1)+ne_found(is)
                   do im = 1+nbasis*(is-1),nbasis*is
                     if(.not. flag_use_overlap) then
                       write(my_pid_energy) real(conjg(V(im,ie))*V(im,ie),kind=4)
                     else
                       write(my_pid_energy) real(conjg(V(im,ie))*SV(im,ie),kind=4)
                     endif
                   enddo
                 enddo
             endif

           endif
         endif

       elseif(.not. flag_print_orbital) then
         if(.not. PINPT%flag_print_single) then
           write(my_pid_energy) ne_found(is), kp(:), (E(ie+nband*(is-1)),ie=1,ne_found(is))
         elseif(PINPT%flag_print_single) then
           write(my_pid_energy) ne_found(is), real(kp(:),kind=4), (real(E(ie+nband*(is-1)),kind=4),ie=1,ne_found(is))
         endif
       endif

       close(my_pid_energy)
       if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
         close(my_pid_energyS)
       endif

     enddo spinb
   endif 
      
   return
endsubroutine
subroutine print_energy( PKPTS, E, E2, V, SV, neig, iband, nband, PINPT, PWGHT, flag_use_overlap, flag_final_step, suffix)
   use parameters, only : pid_energy, incar, poscar, weight, kpoints, zi
   use print_io
   use mpi_setup
   type(incar)  :: PINPT
   type(kpoints):: PKPTS 
   type(weight ):: PWGHT 
   integer*4       mysystem
   integer*4       neig, iband, nband
   integer*4       ie,is,ik,im
   integer*4       nspin, nbasis
   integer*4       ikmode
   integer*4       init_e, fina_e
   integer*4       ne_found(PINPT%nspin, PKPTS%nkpoint)
   integer*4       irecl, i_continue, irec
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   real*8, allocatable :: kpoint_(:,:)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   real*8          E(nband*PINPT%nspin,PKPTS%nkpoint) ! TB energy
   real*8          D(3,nband*PINPT%nspin,PKPTS%nkpoint) ! TB degeneracy info
   real*8          E2(nband*PINPT%nspin,PKPTS%nkpoint) ! Target DFT energy
   real*8          D2(3,nband*PINPT%nspin,PKPTS%nkpoint) ! Target DFT degeneracy info
   complex*16      V(neig*PINPT%ispin,nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      SV(neig*PINPT%ispin,nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      c_up, c_dn
   complex*16      sc_up, sc_dn
   character*80    fname_header, fname_headerS !, suffix
   character(*)    suffix
   character*80    fname, fnameS
   integer*4       pid_energyS
   character*6     kunit_
   character*28    kmode
   character*8     sigma
   logical         flag_print_energy_diff, flag_fit_degeneracy
   logical         flag_final_step ! if true, print D2 - D as well with D (only if flag_print_energy_diff and flag_fit_degeneracy)
                                   ! if false, print D only (only if flag_print_energy_diff and flag_fit_degeneracy)
                                   ! The purpose of this tag is to enable one can track the changes of D and its deviation from D2 (of DFT) 
                                   ! during fitting, it does not save D2 - D during iteration step but after finishing iterations, it can save 
                                   ! D2 - D. Which is only relevant to plot it, in my opinion.
   logical         flag_use_overlap                                  
   real*8          max_wt

   mysystem = PKPTS%mysystem
   max_wt=maxval(PWGHT%WT(:,:))
   if( max_wt .eq. 0) max_wt = 1

   fname_header = trim('band_structure_TBA'//trim(suffix))//trim(PINPT%title(mysystem))
   fname_headerS= trim('band_structure_TBA_S'//trim(suffix))//trim(PINPT%title(mysystem))
   pid_energyS  = pid_energy + 100

   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = PINPT%flag_print_orbital
   flag_print_energy_diff = PINPT%flag_print_energy_diff ! only valid if tbfit = .true. and lorbit = .false. (flag_print_orbital) 
   flag_fit_degeneracy = PINPT%flag_fit_degeneracy

   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = neig
   nspin  = PINPT%nspin
   sigma='sigma_0 '

   if(flag_print_energy_diff .and. flag_fit_degeneracy) then
     call get_degeneracy_serial(E, D, nband*PINPT%nspin, PKPTS%nkpoint, PINPT)
     if(flag_final_step) then
       call get_degeneracy_serial(E2, D2, nband*PINPT%nspin, PKPTS%nkpoint, PINPT)
     endif
   endif

  !if(flag_print_energy_diff .and. flag_print_orbital) then
  !  write(message,'(A)')'    !WARN! You have requested both (LORBIT .TRUE.) and (PRTDIFF .TRUE.)' ; write_msg
  !  write(message,'(A)')'    !WARN! This two function does not work togeter. We deactivate PRTDIFF by default.' ; write_msg
  !endif

   if(flag_klinemode) then
     ikmode = 1
     allocate(kpoint_(1,nkpoint))
   elseif(flag_kgridmode) then
     ikmode = 3
     allocate(kpoint_(3,nkpoint))
   endif

   if(PINPT%flag_sparse) then
     ne_found = PINPT%feast_ne
   else
     ne_found = nband
   endif

   call get_kunit(PKPTS%kunit, kunit_) ! default: angstrom (be carefull)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   init_e = iband ; fina_e = init_e + nband - 1

   if(flag_klinemode) then 
     call get_kline_dist(kpoint, nkpoint, kline)
     kpoint_(1,:) = kline(:)
   elseif(flag_kgridmode) then
     kpoint_ = kpoint
   endif


   if(.not. PINPT%flag_write_unformatted) then
 spin:do is = 1, nspin
        call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear) 
        if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
          call get_fname(fname_headerS, fnameS, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
        endif
        if(flag_final_step) then
          write(message,'(A)') ' '  ; write_msg
          write(message,'(2A)')' #- Writing band structure file: ', trim(fname) ; write_msg
        endif
        open(pid_energy, file=trim(fname), status = 'unknown')
        if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
          open(pid_energyS, file=trim(fnameS), status = 'unknown')
        endif
          if(trim(PINPT%axis_print_mag) .eq. 'rh') then
            write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), &
                  ' ] -> <phi_ij|psi_nk> ; i,j => orbital j in atom i;  n,k => band index (n) and kpoint index (k)'
          elseif(trim(PINPT%axis_print_mag) .eq. 'wf') then
            write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ] -> wavefunction coefficients '
            if(flag_use_overlap) then
              write(pid_energyS,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ] -> overlap matrix multiplied wavefunction coefficients S|V>'
            endif
          elseif(trim(PINPT%axis_print_mag) .eq. 'no') then
            write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ]'
          elseif(trim(PINPT%axis_print_mag) .eq. 'mx' .or. &
                 trim(PINPT%axis_print_mag) .eq. 'my' .or. &
                 trim(PINPT%axis_print_mag) .eq. 'mz') then
            write(pid_energy,'(3A)')'#   MODE LORBIT=[ ',trim(PINPT%axis_print_mag), ' ] -> magnetization <sigma_i>'
          endif
          if(PINPT%flag_sparse) then
            write(pid_energy, '(A,2(F10.4,A),I0)')'# The EWINDOW mode: energy window [EMIN:EMAX]=[ ',PINPT%feast_emin, &
                              ' : ',PINPT%feast_emax,' ], NE_MAX= ',PINPT%feast_nemax
            if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
              write(pid_energyS, '(A,2(F10.4,A),I0)')'# The EWINDOW mode: energy window [EMIN:EMAX]=[ ',PINPT%feast_emin, &
                                ' : ',PINPT%feast_emax,' ], NE_MAX= ',PINPT%feast_nemax
            endif                  
            do ik = 1, nkpoint
              write(pid_energy, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
              if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                write(pid_energyS, '(A,I0,A,I0)')'#   NE_FOUND(ik=',ik,')= ',ne_found(is,ik)
              endif
            enddo
          elseif(PINPT%flag_erange) then
             write(pid_energy, '(A,I0,A,I0,A)')'#   ERANGE=[ ',init_e,' : ',fina_e,' ]'
             if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
               write(pid_energyS, '(A,I0,A,I0,A)')'#   ERANGE=[ ',init_e,' : ',fina_e,' ]'
             endif
          endif
      eig:do ie =1, nband !init_e, fina_e
            if(flag_print_energy_diff) then
              if(.not. flag_fit_degeneracy) write(pid_energy, '(A,I8,A)', ADVANCE = 'yes') trim(kmode), init_e + ie - 1,' -th eigen,    energy(eV),     EDFT-ETBA     kp    eig'
              if(      flag_fit_degeneracy) write(pid_energy, '(A,I8,A)', ADVANCE = 'yes') trim(kmode), init_e + ie - 1,' -th eigen,    energy(eV),   EDFT-ETBA, DTBA_DEGENERACY, DDFT-DTBA'
              if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                if(.not. flag_fit_degeneracy) write(pid_energyS, '(A,I8,A)', ADVANCE = 'yes') trim(kmode), init_e + ie - 1,' -th eigen,    energy(eV),   EDFT-ETBA'
                if(      flag_fit_degeneracy) write(pid_energyS, '(A,I8,A)', ADVANCE = 'yes') trim(kmode), init_e + ie - 1,' -th eigen,    energy(eV),   EDFT-ETBA, DTBA_DEGENERACY, DDFT-DTBA'
              endif
            else
              write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', init_e + ie - 1,' -th eigen'     
              if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                write(pid_energyS, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', init_e + ie - 1,' -th eigen'
              endif
            endif
            if(.not. flag_print_orbital) then
              write(pid_energy,'(A)',ADVANCE='NO')''
            elseif(  flag_print_orbital) then
              if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
              if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
              if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
              if(PINPT%axis_print_mag .ne. 'wf') then
                if(flag_use_overlap) then
                  write(pid_energy, '(2A)',ADVANCE='YES') '# overlap S multiplied wavefunction coeff.: <ci|sigma*S|ci>,sigma=',sigma
                else
                  write(pid_energy, '(2A)',ADVANCE='YES') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
                endif
                write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
              elseif(PINPT%axis_print_mag .eq. 'wf') then
                write(pid_energy, '(1A)',ADVANCE='YES') '# wavefunction coeff.:          |ci>               '
                write(pid_energy, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
                if(flag_use_overlap) then
                  write(pid_energyS, '(1A)',ADVANCE='YES') '# S multiploed wavefunction coeff.:         S|ci>               '
                  write(pid_energyS, '( A)',ADVANCE='NO')  '# k-dist   (ci: wfn coeff for i-th orb)   E(eV), i='
                endif
              endif
           mm:do im=1,nbasis
                if(PINPT%axis_print_mag .ne. 'wf') then
                  write(pid_energy, '(I9)',ADVANCE='NO')im
                elseif(PINPT%axis_print_mag .eq. 'wf') then
                  if(PINPT%ispinor .eq. 2) then
                    write(pid_energy, '(I38)',ADVANCE='NO')im
                    if(flag_use_overlap) then
                      write(pid_energyS, '(I38)',ADVANCE='NO')im
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    write(pid_energy, '(I19)',ADVANCE='NO')im
                    if(flag_use_overlap) then
                      write(pid_energyS, '(I19)',ADVANCE='NO')im
                    endif
                  endif
                endif
                if(im .ge. 30 .and. im .lt. nbasis) then
                  write(pid_energy, '(A)',ADVANCE='NO')' ... '
                  if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                    write(pid_energyS, '(A)',ADVANCE='NO')' ... '
                  endif
                  exit mm
                endif
              enddo mm
              write(pid_energy,'(A)')''
              if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                write(pid_energyS,'(A)')''
              endif
            endif

         kp:do ik = 1, nkpoint
              if(flag_klinemode) then
                if( ie .le. ne_found(is, ik) ) then
                  if(flag_print_energy_diff) then
                    if(flag_fit_degeneracy .and. flag_final_step) then 
                      write(pid_energy,'(1x,F12.6,24x,*(F12.6,1x))',ADVANCE='NO')kline(ik), &
                                                                                 E(ie+nband*(is-1),ik), &
                                                                                 E2(ie+nband*(is-1),ik)-E(ie+nband*(is-1),ik), &
                                                                                 D(1,ie+nband*(is-1),ik), &
                                                                                 D2(1,ie+nband*(is-1),ik) - D(1,ie+nband*(is-1),ik)
                    elseif(flag_fit_degeneracy .and. .not. flag_final_step) then
                      write(pid_energy,'(1x,F12.6,24x,3(F12.6,1x))',ADVANCE='NO')kline(ik), &
                                                                                 E(ie+nband*(is-1),ik), &
                                                                                 E2(ie+nband*(is-1),ik)-E(ie+nband*(is-1),ik), &
                                                                                 D(1,ie+nband*(is-1),ik)
                    elseif(.not. flag_fit_degeneracy) then
                      write(pid_energy,'(1x,F12.6,24x,2(F14.6,1x),2I6)',ADVANCE='NO')kline(ik), E(ie+nband*(is-1),ik), &
                                                                               E2(ie+nband*(is-1),ik)-E(ie+nband*(is-1),ik), &
                                                                               ik,ie+nband*(is-1) ! k, eig index
                    endif
                  else
                    write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+nband*(is-1),ik)
                    if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                      write(pid_energyS,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+nband*(is-1),ik)
                    endif
                  endif
                elseif( ie .gt. ne_found(is, ik)) then
                  write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
                  if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                    write(pid_energyS,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
                  endif
                endif
              elseif(flag_kgridmode) then
                if( ie .le. ne_found(is, ik) ) then
                  if(flag_print_energy_diff) then
                    write(pid_energy,'(1x,3F12.6,2(F12.6,1x),2I6)',ADVANCE='NO')kpoint(:,ik), E(ie+nband*(is-1),ik), &
                                                                            E2(ie+nband*(is-1),ik) - E(ie+nband*(is-1),ik) , &
                                                                            ik,ie+nband*(is-1)
                  else
                    write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+nband*(is-1),ik)
                    if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                      write(pid_energyS,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+nband*(is-1),ik)
                    endif
                  endif
                elseif(ie .gt. ne_found(is, ik)) then
                  write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
                  if(flag_use_overlap .and. trim(PINPT%axis_print_mag) .eq. 'wf') then
                    write(pid_energyS,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
                  endif
                endif
              endif
              if( ie .le. ne_found(is, ik) ) then
                if(flag_print_orbital) then
            basis:do im=1,nbasis-1
                    if(PINPT%ispinor .eq. 2) then
                      c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
                      if(flag_use_overlap) then
                        sc_up = SV(im,ie,ik); sc_dn = SV(im + nbasis,ie,ik)
                      else
                        sc_up =  V(im,ie,ik); sc_dn =  V(im + nbasis,ie,ik)
                      endif
                      if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*sc_up - conjg(c_dn)*sc_dn) ! up - dn : mz
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*sc_up + conjg(c_up)*sc_dn) ! up*dn + dn*up : mx
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*sc_up - conjg(c_up)*sc_dn)*zi) ! (up*dn - dn*up)*i : my
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                        write(pid_energy,'(2(F9.4,F9.4," "))',ADVANCE='NO') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                        if(flag_use_overlap) then
                          write(pid_energyS,'(2(F9.4,F9.4," "))',ADVANCE='NO') sc_up, sc_dn ! sc_up and sc_dn (real,imag) overlap matrix multiplied wavefunction coefficient
                        endif
                      else ! case 'rh', need to be 1 if USE_OVERLAP = .FALSE. otherwise does not need to be.
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn) ! up + dn : total
                      endif
                    elseif(PINPT%ispinor .eq. 1) then
                      c_up = V(im+nbasis*(is-1),ie+nband*(is-1),ik)
                      if(flag_use_overlap) then
                        sc_up = SV(im+nbasis*(is-1),ie+nband*(is-1),ik)
                      else
                        sc_up =  V(im+nbasis*(is-1),ie+nband*(is-1),ik)
                      endif
                      if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                        write(pid_energy,'(1(F9.4,F9.4," "))',ADVANCE='NO') c_up
                        if(flag_use_overlap) then
                          write(pid_energyS,'(1(F9.4,F9.4," "))',ADVANCE='NO') sc_up
                        endif
                      else
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*sc_up)
                      endif
                    endif
                  enddo basis
                  if(PINPT%ispinor .eq. 2) then
                    c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
                    if(flag_use_overlap) then
                      sc_up = SV(im,ie,ik); sc_dn = SV(im + nbasis,ie,ik)
                    else
                      sc_up =  V(im,ie,ik); sc_dn =  V(im + nbasis,ie,ik)
                    endif
                    if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*sc_up - conjg(c_dn)*sc_dn) ! up - dn : mz
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_dn)*sc_up + conjg(c_up)*sc_dn) ! up*dn + dn*up : mx
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real((conjg(c_dn)*sc_up - conjg(c_up)*sc_dn)*zi) ! (up*dn - dn*up)*i : my
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                      write(pid_energy,'(2(F9.4,F9.4," "))',ADVANCE='YES') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                      if(flag_use_overlap) then
                        write(pid_energyS,'(2(F9.4,F9.4," "))',ADVANCE='YES') sc_up, sc_dn ! sc_up and sc_dn (real,imag) overlap matrix multiplied wavefunction coefficient
                      endif
                    else ! case 'rh', need to be 1 if USE_OVERLAP = .FALSE. otherwise does not need to be.
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn) ! up + dn : total
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    c_up = V(nbasis+nbasis*(is-1),ie+nband*(is-1),ik)
                    if(flag_use_overlap) then
                      sc_up = SV(nbasis+nbasis*(is-1),ie+nband*(is-1),ik)
                    else
                      sc_up =  V(nbasis+nbasis*(is-1),ie+nband*(is-1),ik)
                    endif
                    if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                      write(pid_energy,'(1(F9.4,F9.4," "))',ADVANCE='YES') c_up
                      if(flag_use_overlap) then
                        write(pid_energyS,'(1(F9.4,F9.4," "))',ADVANCE='YES') sc_up
                      endif
                    else
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*sc_up)
                    endif
                  endif
                else
                  write(pid_energy,'(A)',ADVANCE='YES')''
                  if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
                    write(pid_energyS,'(A)',ADVANCE='YES')''
                  endif
                endif
              elseif(ie .gt. ne_found(is, ik)) then
                write(pid_energy,*)''
                if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
                  write(pid_energyS,*)''
                endif
              endif
            enddo kp
            write(pid_energy,'(A)',ADVANCE='YES')''
            write(pid_energy,'(A)',ADVANCE='YES')''
            if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
              write(pid_energyS,'(A)',ADVANCE='YES')''
              write(pid_energyS,'(A)',ADVANCE='YES')''
            endif
          enddo eig

        close(pid_energy)
        if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
          close(pid_energyS)
        endif
      enddo spin

    elseif(PINPT%flag_write_unformatted) then ! only work with flag_tbfit = .false.
spinb:do is = 1, nspin
        call get_fname_bin(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
        if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
          call get_fname_bin(fname_headerS, fnameS, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
        endif
        ! write header
        write(message,'(A)') ' '  ; write_msg
        write(message,'(2A)')' #- Writing band structure file: ', trim(fname) ; write_msg
        open(pid_energy, file=trim(fname), form='unformatted', status='unknown')
        write(pid_energy) ikmode, flag_print_orbital, PINPT%flag_print_single, PINPT%flag_erange, & 
                          PINPT%flag_sparse, neig, PKPTS%nkpoint, nband, &
                          PINPT%ispin, PINPT%nspin, PINPT%ispinor, PINPT%axis_print_mag
        if(PINPT%flag_erange) then
          write(pid_energy) PINPT%flag_erange, init_e , fina_e 
        else
          write(pid_energy) PINPT%flag_erange ! .FALSE.
        endif
        if(PINPT%flag_sparse) then
          write(pid_energy) PINPT%flag_sparse, PINPT%feast_emin, PINPT%feast_emax, PINPT%feast_nemax
        else
          write(pid_energy) PINPT%flag_sparse ! .FALSE.
        endif

        if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
          open(pid_energyS, file=trim(fnameS), form='unformatted', status='unknown')
          write(pid_energyS) ikmode, flag_print_orbital, PINPT%flag_print_single, PINPT%flag_erange, & 
                            PINPT%flag_sparse, neig, PKPTS%nkpoint, nband, &
                            PINPT%ispin, PINPT%nspin, PINPT%ispinor, PINPT%axis_print_mag
          if(PINPT%flag_erange) then
            write(pid_energyS) PINPT%flag_erange, init_e , fina_e
          else
            write(pid_energyS) PINPT%flag_erange ! .FALSE.
          endif
          if(PINPT%flag_sparse) then
            write(pid_energyS) PINPT%flag_sparse, PINPT%feast_emin, PINPT%feast_emax, PINPT%feast_nemax
          else
            write(pid_energyS) PINPT%flag_sparse ! .FALSE.
          endif
        endif

        ! write main wavefunction information
          if(flag_print_orbital) then
            if(PINPT%ispinor .eq. 2) then
              if(PINPT%axis_print_mag .eq. 'wf') then
                if(.not.PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                    if(flag_use_overlap) then
                      write(pid_energyS) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                    endif
                    do ie = 1, ne_found(is, ik)
                      do im = 1, nbasis
                        write(pid_energy) V(im,ie,ik), V(im+nbasis,ie,ik)
                        if(flag_use_overlap) write(pid_energyS) SV(im,ie,ik), SV(im+nbasis,ie,ik)
                      enddo
                    enddo
                  enddo
                elseif(PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
                    do ie = 1, ne_found(is, ik)
                      do im = 1, nbasis
                        if(.not.flag_use_overlap) then
                          write(pid_energy) cmplx((/V(im,ie,ik), V(im+nbasis,ie,ik)/),kind=4)
                        else
                          write(pid_energy) cmplx((/V(im,ie,ik), SV(im+nbasis,ie,ik)/),kind=4)
                        endif
                      enddo
                    enddo
                  enddo

                endif
              elseif(PINPT%axis_print_mag .eq. 'rh') then
                if(.not.PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                    do ie = 1, ne_found(is,ik)
                      do im = 1, nbasis
                        if(.not.flag_use_overlap) then
                          write(pid_energy) real( conjg(V(im,ie,ik))*V(im,ie,ik) + conjg(V(im+nbasis,ie,ik))*V(im+nbasis,ie,ik) )
                        else
                          write(pid_energy) real( conjg(V(im,ie,ik))*SV(im,ie,ik) + conjg(V(im+nbasis,ie,ik))*SV(im+nbasis,ie,ik) )
                        endif
                      enddo
                    enddo
                  enddo
                  
                elseif(PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
                    do ie = 1, ne_found(is,ik)
                      do im = 1, nbasis
                        if(.not. flag_use_overlap) then
                          write(pid_energy) real( conjg(V(im,ie,ik))*V(im,ie,ik)+conjg(V(im+nbasis,ie,ik))*V(im+nbasis,ie,ik), kind=4 )
                        else
                          write(pid_energy) real( conjg(V(im,ie,ik))*SV(im,ie,ik)+conjg(V(im+nbasis,ie,ik))*SV(im+nbasis,ie,ik), kind=4 )
                        endif
                      enddo
                    enddo
                  enddo
                 
                endif
              endif
            elseif(PINPT%ispinor .eq. 1) then
              if(PINPT%axis_print_mag .eq. 'wf') then
                if(.not.PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                    if(flag_use_overlap) then
                      write(pid_energyS) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                    endif
                    do ie = 1+nband*(is-1), nband*(is-1)+ne_found(is,ik)
                      do im = 1+nbasis*(is-1),nbasis*is
                        write(pid_energy) V(im,ie,ik)
                        if(flag_use_overlap) then
                          write(pid_energyS) SV(im,ie,ik)
                        endif
                      enddo
                    enddo
                  enddo
                 
                elseif(PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
                    if(flag_use_overlap) then
                      write(pid_energyS) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
                    endif
                    do ie = 1+nband*(is-1), nband*(is-1)+ne_found(is,ik)
                      do im = 1+nbasis*(is-1),nbasis*is
                        write(pid_energy) cmplx(V(im,ie,ik), kind=4)
                        if(flag_use_overlap) then
                          write(pid_energyS) cmplx(SV(im,ie,ik), kind=4)
                        endif
                      enddo
                    enddo
                  enddo
                  
                endif
              elseif(PINPT%axis_print_mag .eq. 'rh') then
                if(.not.PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
                    do ie = 1+nband*(is-1),nband*(is-1)+ne_found(is,ik) 
                      do im = 1+nbasis*(is-1),nbasis*is
                        if(.not. flag_use_overlap) then
                          write(pid_energy) real(conjg(V(im,ie,ik))*V(im,ie,ik))
                        else
                          write(pid_energy) real(conjg(V(im,ie,ik))*SV(im,ie,ik))
                        endif
                      enddo
                    enddo
                  enddo

                elseif(PINPT%flag_print_single) then
                  do ik = 1, nkpoint
                    write(pid_energy) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
                    do ie = 1+nband*(is-1),nband*(is-1)+ne_found(is,ik) 
                      do im = 1+nbasis*(is-1),nbasis*is
                        if(.not. flag_use_overlap) then
                          write(pid_energy) real(conjg(V(im,ie,ik))*V(im,ie,ik),kind=4)
                        else
                          write(pid_energy) real(conjg(V(im,ie,ik))*SV(im,ie,ik),kind=4)
                        endif
                      enddo
                    enddo
                  enddo

                endif
              endif
            endif
          else
            if(.not. PINPT%flag_print_single) then
              do ik = 1, nkpoint
                write(pid_energy) ne_found(is,ik), kpoint_(:,ik), (E(ie+nband*(is-1),ik),ie=1,ne_found(is,ik))
              enddo
            elseif(PINPT%flag_print_single) then
              do ik = 1, nkpoint
                write(pid_energy) ne_found(is,ik), real(kpoint_(:,ik),kind=4), (real(E(ie+nband*(is-1),ik),kind=4),ie=1,ne_found(is,ik))
              enddo
            endif
          endif

        close(pid_energy)
        if(flag_use_overlap .and. PINPT%axis_print_mag .eq. 'wf') then
          close(pid_energyS)
        endif
      enddo spinb

    endif

    deallocate(kpoint_)
return
endsubroutine

subroutine get_kline_dist(kpoint, nkpoint, kline)
   implicit none
   integer*4    ik, nkpoint
   real*8       kline(nkpoint),k0(3), enorm
   real*8       kpoint(3,nkpoint)
   external     enorm

   do ik=1,nkpoint
     if(ik .eq. 1) then
      k0=kpoint(:,1)
      kline(1)=0
     else
      k0=kpoint(:,ik-1)
      kline(ik)=kline(ik-1)
     endif
     kline(ik)=kline(ik)+ enorm(3, kpoint(:,ik)-k0(:) )
   enddo

return
endsubroutine
subroutine get_fname_bin(fname_header, fname, is, flag_collinear, flag_noncollinear)
   implicit none
   integer*4    is
   character(*) fname_header
   character*80 fname
   logical      flag_noncollinear, flag_collinear

   fname = ' '

   if(flag_noncollinear) then
     write(fname, '(A,A4)')trim(fname_header),'.bin'
   elseif(flag_collinear) then
     if(is .eq. 1) then
       write(fname,    '(A,A7)')trim(fname_header),'.up.bin'
     elseif(is .eq. 2) then
       write(fname, '(A,A7)')trim(fname_header),'.dn.bin'
     endif
   elseif(.not. flag_noncollinear .and. .not. flag_collinear) then
     write(fname, '(A,A4)')trim(fname_header),'.bin'
   endif

return
endsubroutine
subroutine get_fname(fname_header, fname, is, flag_collinear, flag_noncollinear)
   implicit none
   integer*4    is
   character(*) fname_header
   character*80 fname
   logical      flag_noncollinear, flag_collinear

   fname = ' '

   if(flag_noncollinear) then
     write(fname, '(A,A4)')trim(fname_header),'.dat'
   elseif(flag_collinear) then
     if(is .eq. 1) then
       write(fname,    '(A,A7)')trim(fname_header),'.up.dat'
     elseif(is .eq. 2) then
       write(fname, '(A,A7)')trim(fname_header),'.dn.dat'
     endif
   elseif(.not. flag_noncollinear .and. .not. flag_collinear) then
     write(fname, '(A,A4)')trim(fname_header),'.dat'  
   endif

return
endsubroutine

subroutine get_e_range(init_e, fina_e, neig, flag_ensurf, ispinor, flag_erange, init_erange, fina_erange)
   implicit none
   integer*4      init_e, fina_e, neig, ispinor
   integer*4      init_erange, fina_erange
   logical        flag_ensurf, flag_erange

   if(.not. flag_ensurf) then
     if(flag_erange) then
       init_e = init_erange
       fina_e = fina_erange
     elseif(.not. flag_erange) then
       init_e = 1
       fina_e = neig * ispinor
     endif

   elseif(flag_ensurf) then
     init_e   = 1
     fina_e   = 1
   endif
return
endsubroutine

subroutine get_kunit(kunit, kunit_)
   implicit none
   character*1  kunit
   character*6  kunit_

   if(kunit .eq. 'R') then
     kunit_ = '(reci)'
   elseif(kunit .eq. 'A') then
     kunit_ = '(A^-1)'
   endif

return
endsubroutine
subroutine get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   implicit none
   logical      flag_klinemode, flag_kgridmode
   character*6  kunit_
   character*28 kmode

   if(flag_klinemode) then
     write(kmode,'(3A)')'#    k-dist ',trim(kunit_),' '
   elseif(flag_kgridmode) then
     write(kmode,'(3A)')'#    k-point',trim(kunit_),' '
   endif

return
endsubroutine
subroutine print_energy_weight (kpoint, nkpoint, EDFT, PWGHT, neig, PINPT, fname, flag_print_order)
  use parameters, only : energy, weight, incar, pid_energy
  implicit none
  type(energy)  :: EDFT
  type(weight)  :: PWGHT
  type(incar )  :: PINPT
  integer*4, intent(in) :: neig, nkpoint
  integer*4 i,ie,ik, ispin, neig_
  real*8 kpoint(3,nkpoint), kline(nkpoint),k0(3), enorm
  real*8 max_wt
  character(*) fname
  logical  flag_collinear, flag_noncollinear
  logical  flag_print_order
  external enorm
! pid_energy=32

  max_wt=maxval(PWGHT%WT(:,:))
  if( max_wt .eq. 0) max_wt = 1

  ispin = PINPT%ispin
  flag_collinear = PINPT%flag_collinear
  flag_noncollinear = PINPT%flag_noncollinear

  if(flag_noncollinear) then
    neig_ = neig * 2
  else
    neig_ = neig
  endif

  do ik=1,nkpoint
    if(ik .eq. 1) then
     k0=kpoint(:,1)
     kline(1)=0
    else
     k0=kpoint(:,ik-1)
     kline(ik)=kline(ik-1)
    endif
    kline(ik)=kline(ik)+ enorm(3, kpoint(:,ik)-k0(:) )
  enddo

  open(pid_energy, file=trim(fname))
  do ie=1,neig_
   if(flag_collinear) then
     if(.not.PINPT%flag_fit_degeneracy) then
       write(pid_energy, '(A,I8,A)') '#k-dist(A^-1)UPenergy(eV)   weight DNenergy(eV)   weight, ORDER_Eu(eV) ORDER_Ed(eV),', ie,' -th eigen'
     elseif(PINPT%flag_fit_degeneracy) then
       write(pid_energy, '(2A,I8,A)')'#k-dist(A^-1)UPenergy(eV)   weight DNenergy(eV)   weight,',&
                                                  'UPdegeneracy   UPdegeneracy ', ie,' -th eigen'
     endif
   else
     if(.not.PINPT%flag_fit_degeneracy) then
       write(pid_energy, '(A,I8,A)') '#k-dist(A^-1)  energy(eV)   weight   en_ord(eV),', ie,' -th eigen'
     elseif(PINPT%flag_fit_degeneracy) then
       write(pid_energy, '(A,I8,A)') '#k-dist(A^-1)  energy(eV)   weight   degeneracy,', ie,' -th eigen'
     endif
   endif
   do ik=1,nkpoint
    if(flag_collinear) then
      if(.not.PINPT%flag_fit_degeneracy) then
        if(.not. flag_print_order) then
          write(pid_energy,'(2x,F9.6, 2(F12.6, F11.3) )',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt, &
                                                                                 EDFT%E(ie+neig,ik), PWGHT%WT(ie+neig,ik)/max_wt
        elseif(flag_print_order) then ! note: printing re-ordered target energy is only possible without "fit_degeneracy" option
          write(pid_energy,'(2x,F9.6, 3(F12.6, F11.3), F12.6, F12.6 )',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt, &
                                                                                 EDFT%E(ie+neig,ik), PWGHT%WT(ie+neig,ik)/max_wt, &
                                                                                 EDFT%E_ORD(ie,ik), &
                                                                                 EDFT%E_ORD(ie+neig,ik)
        endif
      elseif(PINPT%flag_fit_degeneracy) then
        write(pid_energy,'(2x,F9.6, 2(F12.6, F11.3), 2F11.3 )',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt, &
                                                                               EDFT%E(ie+neig,ik), PWGHT%WT(ie+neig,ik)/max_wt, &
                                                                               EDFT%D(1,ie,ik), EDFT%D(1,ie+neig,ik)
      endif
    else
      if(.not.PINPT%flag_fit_degeneracy) then
        if(.not. flag_print_order) then
          write(pid_energy,'(2x,F9.6, F12.6, F12.6)',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt
        elseif(flag_print_order) then
          write(pid_energy,'(2x,F9.6, F12.6, F12.6, F12.6)',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt, EDFT%E_ORD(ie,ik)
        endif
      elseif(PINPT%flag_fit_degeneracy) then
        write(pid_energy,'(2x,F9.6, F12.6,2F11.3)',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt, EDFT%D(1,ie,ik)
      endif
    endif
    write(pid_energy,*)''
   enddo !ik
   write(pid_energy,*)''
   write(pid_energy,*)''
  enddo !ie
  close(pid_energy)

return
endsubroutine
