subroutine print_energy_ensurf (kpoint, nkpoint,ie, ispin_print, E, V, PGEOM, PINPT, fname_header, kunit)
   use parameters, only : pid_energy, incar, poscar, zi
   type(incar) :: PINPT
   type(poscar):: PGEOM
   integer*4      is, ie, ik
   integer*4      ispin_print, nbasis
   real*8         kpoint(3,nkpoint)
   logical        flag_print_orbital
   real*8         E(ispin_print,nkpoint)
   complex*16     V(PGEOM%neig*PINPT%ispin,ispin_print,nkpoint)  
   complex*16     c_up, c_dn
   character(*)   fname_header
   character*80   fname
   character*1    kunit
   character*6    kunit_
   character*8    sigma

   sigma='sigma_0 '
   nbasis = PGEOM%neig
   call get_kunit(kunit, kunit_)
   call get_plotmode(.false., .true., kunit_, kmode)
 spin:do is = 1, ispin_print
        call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
        open(pid_energy, file=trim(fname), status = 'unknown')

            write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', ie,' -th eigen'
            if(.not. PINPT%flag_print_orbital) then
              write(pid_energy,'(A)',ADVANCE='yes')''
            elseif(  PINPT%flag_print_orbital) then
              if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
              if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
              if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
              write(pid_energy, '(2A)',ADVANCE='NO') '# wavefunction coeff.: <ci|sigma|ci>,sigma= ',sigma
              do im=1,nbasis-1
                write(pid_energy, '(I9)',ADVANCE='NO')im
              enddo
              write(pid_energy,'(I9)',ADVANCE='yes')im
            endif

         kp:do ik = 1, nkpoint
              write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(1+(is-1),ik)
              if(PINPT%flag_print_orbital) then
          basis:do im=1,nbasis-1
                  if(PINPT%ispinor .eq. 2) then
                    c_up = V(im,is,ik); c_dn = V(im + nbasis,is,ik)
                    if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz 
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
                    else
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    c_up = V(im+PGEOM%neig*(is-1),is,ik)
                    write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*c_up)
                  endif
                enddo basis
                if(PINPT%ispinor .eq. 2) then
                  c_up = V(im,is,ik); c_dn = V(im + nbasis,is,ik)
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
                  c_up = V(nbasis+PGEOM%neig*(is-1),is,ik)
                  write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*c_up)
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
   integer*4       ispin_print, nbasis, nkpoint
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
   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = PINPT%flag_print_orbital
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = neig
   sigma='sigma_0 '

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   call get_ispin_print(PINPT%flag_collinear, ispin_print)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

spin:do is = 1, ispin_print
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
         write(pid_energy, '(2A)',ADVANCE='NO') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
         do im=1,nbasis
           write(pid_energy, '(I9)',ADVANCE='NO')im
         enddo
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
subroutine print_energy_ldos(PKPTS,E,V,PGEOM,PINPT)
   use parameters, only: pid_energy, incar, poscar, kpoints, zi
   type(incar)  :: PINPT
   type(poscar) :: PGEOM
   type(kpoints):: PKPTS
   integer*4       ie,is,ik,im,ia
   integer*4       ispin_print, nbasis, natom
   integer*4       init_e, fina_e
   integer*4       ne_found(PINPT%nspin, PKPTS%nkpoint)
   integer*4       imatrix
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   real*8          E(PINPT%nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      V(PGEOM%neig*PINPT%ispin,PINPT%nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      c_up, c_dn, c_tot
   character*80    fname_header
   character*80    fname
   character*6     kunit_
   character*28    kmode
   character*8     sigma

   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = PINPT%flag_get_orbital
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = PGEOM%neig
   sigma='sigma_0 '
   natom= PGEOM%n_atom
   
 
   if(PINPT%flag_sparse) then
     ne_found = PINPT%feast_ne
   else
     ne_found = PINPT%nband
   endif

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   call get_e_range(init_e, fina_e, PGEOM%neig, .false., PINPT)
   call get_ispin_print(PINPT%flag_collinear, ispin_print)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)
 
spin:do is = 1, ispin_print
   atom:do ia = 1, natom
       imatrix = sum( PGEOM%n_orbital(1:ia) ) - PGEOM%n_orbital(ia) + 1
       write(fname_header,'(A,I0)')'band_structure_TBA_atom.',ia 
       call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
       open(pid_energy, file = trim(fname), status = 'unknown')
       
       if(PINPT%flag_sparse) then
         write(pid_energy, '(A,2(F10.4,A))')'# The EWINDOW mode: energy window [EMIN:EMAX]=[', &
                                             PINPT%feast_emin,':', PINPT%feast_emax,']'
       endif
   eig:do ie = 1, PINPT%nband ! init_e, fina_e
         write(pid_energy,'(2A,I8,A,I8,3A)',ADVANCE='yes')kmode,'  energy(eV) :',init_e+ie-1,' -th eigen | ',ia, &
                                                    ' -th atom (spec= ',trim(PGEOM%c_spec(PGEOM%spec(ia))),' )'
         if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
         if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
         if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
         
         write(pid_energy, '(2A)',ADVANCE='NO') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma     
         do im=imatrix, imatrix + PGEOM%n_orbital(ia) - 1
           write(pid_energy, '(I9)',ADVANCE='NO')im
         enddo
         write(pid_energy,'(A9)',ADVANCE='YES') ' tot'
!        write(pid_energy,'(A)')''     

      kp:do ik = 1, nkpoint
           if(flag_klinemode) then
             if( ie .le. ne_found(is, ik) ) then
               write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+PINPT%nband*(is-1),ik)
             elseif( ie .gt. ne_found(is, ik)) then
               write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
             endif
           elseif(flag_kgridmode) then
             if( ie .le. ne_found(is, ik) ) then
               write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+PINPT%nband*(is-1),ik)
             elseif(ie .gt. ne_found(is, ik)) then
               write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
             endif
           endif
           if( ie .le. ne_found(is, ik) ) then
             if(flag_print_orbital) then
               c_tot = 0d0 !initialize
         basis:do im=imatrix, imatrix+PGEOM%n_orbital(ia) - 1
                 if(PINPT%ispinor .eq. 2) then
                   c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
                   if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
                     c_tot = c_tot + real( conjg(c_up)*c_up - conjg(c_dn)*c_dn)
                   elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
                     c_tot = c_tot + real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) 
                   elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
                     c_tot = c_tot + real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi)
                   else
                     write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
                     c_tot = c_tot + real( conjg(c_up)*c_up + conjg(c_dn)*c_dn)
                   endif
                 elseif(PINPT%ispinor .eq. 1) then
                   c_up = V(im+PGEOM%neig*(is-1),ie+PINPT%nband*(is-1),ik)
                   write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*c_up)
                   c_tot = c_tot + real(conjg(c_up)*c_up)
                 endif
               enddo basis
               write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(c_tot)

!              if(PINPT%ispinor .eq. 2) then
!                c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
!                if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
!                  write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
!                elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
!                  write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
!                elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
!                  write(pid_energy,'(*(F9.4))',ADVANCE='YES') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
!                else
!                  write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
!                endif
!              elseif(PINPT%ispinor .eq. 1) then
!                c_up = V(nbasis+PGEOM%neig*(is-1),ie+PINPT%nband*(is-1),ik)
!                write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*c_up)
!              endif
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
     enddo atom
   enddo spin

return
endsubroutine
subroutine print_energy( PKPTS, E, V, PGEOM, PINPT)
   use parameters, only : pid_energy, incar, poscar, kpoints, zi
   type(incar)  :: PINPT
   type(poscar) :: PGEOM 
   type(kpoints):: PKPTS 
   integer*4       ie,is,ik,im
   integer*4       ispin_print, nbasis
   integer*4       init_e, fina_e
   integer*4       ne_found(PINPT%nspin, PKPTS%nkpoint)
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   real*8          E(PINPT%nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      V(PGEOM%neig*PINPT%ispin,PINPT%nband*PINPT%nspin,PKPTS%nkpoint)
   complex*16      c_up, c_dn
   character*80    fname_header
   character*80    fname
   character*6     kunit_
   character*28    kmode
   character*8     sigma
   fname_header = 'band_structure_TBA'
   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = PINPT%flag_print_orbital
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = PGEOM%neig
   sigma='sigma_0 '


   if(PINPT%flag_sparse) then
     ne_found = PINPT%feast_ne
   else
     ne_found = PINPT%nband
   endif

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   call get_e_range(init_e, fina_e, PGEOM%neig, .false., PINPT)
   call get_ispin_print(PINPT%flag_collinear, ispin_print)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

 spin:do is = 1, ispin_print
        call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear) 
        open(pid_energy, file=trim(fname), status = 'unknown')
          if(PINPT%flag_sparse) then
            write(pid_energy, '(A,2(F10.4,A))')'# The EWINDOW mode: energy window [EMIN:EMAX]=[',PINPT%feast_emin,':',PINPT%feast_emax,']'
          endif
      eig:do ie =1, PINPT%nband !init_e, fina_e
            write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', init_e + ie - 1,' -th eigen'     
            if(.not. flag_print_orbital) then
              write(pid_energy,'(A)',ADVANCE='NO')''
            elseif(  flag_print_orbital) then
              if(PINPT%axis_print_mag .eq. 'mz') sigma='sigma_z '
              if(PINPT%axis_print_mag .eq. 'mx') sigma='sigma_x '
              if(PINPT%axis_print_mag .eq. 'my') sigma='sigma_y '
              if(PINPT%axis_print_mag .ne. 'wf') then
                write(pid_energy, '(2A)',ADVANCE='NO') '# wavefunction coeff.: <ci|sigma|ci>,sigma=',sigma
              elseif(PINPT%axis_print_mag .eq. 'wf') then
                write(pid_energy, '(1A)',ADVANCE='NO') '# wavefunction coeff.:          |ci>               '
              endif
              do im=1,nbasis
                if(PINPT%axis_print_mag .ne. 'wf') then
                  write(pid_energy, '(I9)',ADVANCE='NO')im
                elseif(PINPT%axis_print_mag .eq. 'wf') then
                  write(pid_energy, '(I38)',ADVANCE='NO')im
                endif
              enddo
              write(pid_energy,'(A)')''

!             if(PINPT%ispinor .eq. 2 .and. PINPT%axis_print_mag .eq. 'wf') then
!               do im=1,nbasis
!                 write(pid_energy, '(51x, 2A19)',ADVANCE='NO')'up','dn'
!               enddo
!               write(pid_energy,'(A)')''
!             endif
            endif

         kp:do ik = 1, nkpoint
              if(flag_klinemode) then
                if( ie .le. ne_found(is, ik) ) then
                  write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik), E(ie+PINPT%nband*(is-1),ik)
                elseif( ie .gt. ne_found(is, ik)) then
                  write(pid_energy,'(1x,F12.6,24x,F14.6,1x)',ADVANCE='NO')kline(ik)
                endif
              elseif(flag_kgridmode) then
                if( ie .le. ne_found(is, ik) ) then
                  write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik), E(ie+PINPT%nband*(is-1),ik)
                elseif(ie .gt. ne_found(is, ik)) then
                  write(pid_energy,'(1x,3F12.6,F14.6,1x)',ADVANCE='NO')kpoint(:,ik)
                endif
              endif
              if( ie .le. ne_found(is, ik) ) then
                if(flag_print_orbital) then
            basis:do im=1,nbasis-1
                    if(PINPT%ispinor .eq. 2) then
                      c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)

                      if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
                      elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                        write(pid_energy,'(2(F9.4,F9.4,"i"))',ADVANCE='NO') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                      else
                        write(pid_energy,'(*(F9.4))',ADVANCE='NO') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
                      endif
                    elseif(PINPT%ispinor .eq. 1) then
                      c_up = V(im+PGEOM%neig*(is-1),ie+PINPT%nband*(is-1),ik)
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO') real(conjg(c_up)*c_up)
                    endif
                  enddo basis
                  if(PINPT%ispinor .eq. 2) then
                    c_up = V(im,ie,ik); c_dn = V(im + nbasis,ie,ik)
                    if    (PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up - conjg(c_dn)*c_dn) ! up - dn : mz
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mx') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_dn)*c_up + conjg(c_up)*c_dn) ! up*dn + dn*up : mx
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'my') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real((conjg(c_dn)*c_up - conjg(c_up)*c_dn)*zi) ! (up*dn - dn*up)*i : my
                    elseif(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'wf') then
                      write(pid_energy,'(2(F9.4,F9.4,"i"))',ADVANCE='YES') c_up, c_dn ! c_up and c_dn (real,imag) wavefunction coefficient
                    else
                      write(pid_energy,'(*(F9.4))',ADVANCE='YES') real( conjg(c_up)*c_up + conjg(c_dn)*c_dn) ! up + dn : total
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    c_up = V(nbasis+PGEOM%neig*(is-1),ie+PINPT%nband*(is-1),ik)
                    write(pid_energy,'(*(F9.4))',ADVANCE='YES') real(conjg(c_up)*c_up)
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

subroutine get_e_range(init_e, fina_e, neig, flag_ensurf, PINPT)
   use parameters, only : incar, poscar
   implicit none
   type(incar) :: PINPT
   integer*4      init_e, fina_e, neig
   logical        flag_ensurf 

   if(.not. flag_ensurf) then
     if(PINPT%flag_erange) then
       init_e = PINPT%init_erange
       fina_e = PINPT%fina_erange
     elseif(.not. PINPT%flag_erange) then
       init_e = 1
       fina_e = neig * PINPT%ispinor
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
     write(kmode,'(3A)')'#        k-dist ',kunit_,'      '
   elseif(flag_kgridmode) then
     write(kmode,'(3A)')'#        k-point',kunit_,'      '
   endif

return
endsubroutine
subroutine get_ispin_print(flag_collinear, ispin_print)
   implicit none
   logical       flag_collinear
   integer*4     ispin_print

   if(flag_collinear) then

     ispin_print = 2
  
   else

     ispin_print = 1

   endif

return
endsubroutine
subroutine print_energy_weight (kpoint, nkpoint, EDFT, PWGHT, neig, PINPT, fname)
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
     write(pid_energy, '(A,I8,A)') '#k-dist(A^-1)UPenergy(eV)   weight DNenergy(eV)   weight,', ie,' -th eigen'
   else
     write(pid_energy, '(A,I8,A)') '#k-dist(A^-1)  energy(eV)   weight,', ie,' -th eigen'
   endif
   do ik=1,nkpoint
    if(flag_collinear) then
      write(pid_energy,'(2x,F9.6, 2(F12.6, F11.3) )',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt, &
                                                                             EDFT%E(ie+neig,ik), PWGHT%WT(ie+neig,ik)/max_wt
    else
      write(pid_energy,'(2x,F9.6, F12.6, F11.3)',ADVANCE='NO')kline(ik), EDFT%E(ie,ik), PWGHT%WT(ie,ik)/max_wt
    endif
    write(pid_energy,*)''
   enddo !ik
   write(pid_energy,*)''
   write(pid_energy,*)''
  enddo !ie
  close(pid_energy)

return
endsubroutine


