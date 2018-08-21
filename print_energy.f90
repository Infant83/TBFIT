subroutine print_energy_ensurf (kpoint, nkpoint, ispin_print, E, V, PGEOM, PINPT, fname_header, kunit)
   use parameters, only : pid_energy, incar, poscar
   type(incar) :: PINPT
   type(poscar):: PGEOM
   integer*4      is, ie, ik
   integer*4      ispin_print, nbasis
   integer*4      init_e, fina_e
   real*8         kpoint(3,nkpoint)
   logical        flag_print_orbital
   real*8         E(ispin_print,nkpoint)
   complex*16     V(PGEOM%neig*PINPT%ispin,ispin_print,nkpoint)  
   character(*)   fname_header
   character*80   fname
   character*1    kunit
   character*6    kunit_

   nbasis = PGEOM%neig
   call get_kunit(kunit, kunit_)
   call get_plotmode(.false., .true., kunit_, kmode)
   call get_e_range(init_e, fina_e, PGEOM, .true., PINPT)
  
 spin:do is = 1, ispin_print
        call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
        open(pid_energy, file=trim(fname), status = 'unknown')

      eig:do ie =init_e, fina_e
            write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', ie,' -th eigen'
            if(.not. PINPT%flag_print_orbital) then
              write(pid_energy,'(A)',ADVANCE='yes')' '
            elseif(  PINPT%flag_print_orbital) then
              write(pid_energy, '(A)',ADVANCE='NO') '# wavefunction coeff.:                              '
              do im=1,nbasis-1
                write(pid_energy, '(I9)',ADVANCE='NO')im
              enddo
              write(pid_energy,'(I9)',ADVANCE='yes')im
            endif

         kp:do ik = 1, nkpoint
              write(pid_energy,'(1x,3F12.6,F12.6,3x)',ADVANCE='NO')kpoint(:,ik), E(ie+(is-1),ik)
              if(PINPT%flag_print_orbital) then
          basis:do im=1,nbasis-1
                  if(PINPT%ispinor .eq. 2) then
                    write(pid_energy,'(*(F9.4))',ADVANCE='NO')abs(V(im,ie,ik)) + abs(V(im + nbasis, ie, ik))
                  elseif(PINPT%ispinor .eq. 1) then
                    write(pid_energy,'(*(F9.4))',ADVANCE='NO')abs(V(im+PGEOM%neig*(is-1),ie+(is-1),ik))
                  endif
                enddo basis
                if(PINPT%ispinor .eq. 2) then
                  write(pid_energy,'(*(F9.4))',ADVANCE='YES')abs(V(nbasis,ie,ik)) + abs(V(nbasis*2, ie, ik))
                elseif(PINPT%ispinor .eq. 1) then
                  write(pid_energy,'(*(F9.4))',ADVANCE='YES')abs(V(nbasis+PGEOM%neig*(is-1),ie+(is-1),ik))
                endif
              endif
              if(.not.PINPT%flag_print_orbital) write(pid_energy,*)''
            enddo kp
            write(pid_energy,*)''
            write(pid_energy,*)''
          enddo eig

      close(pid_energy)

      enddo spin
   

return
endsubroutine
subroutine print_energy( PKPTS, E, V, PGEOM, PINPT)
   use parameters, only : pid_energy, incar, poscar, kpoints
   type(incar)  :: PINPT
   type(poscar) :: PGEOM 
   type(kpoints):: PKPTS 
   integer*4       ie,is,ik
   integer*4       ispin_print, nbasis
   integer*4       init_e, fina_e
   real*8          kline(PKPTS%nkpoint),kpoint(3,PKPTS%nkpoint)
   logical         flag_klinemode, flag_kgridmode, flag_print_orbital
   real*8          E(PGEOM%neig*PINPT%ispin,PKPTS%nkpoint)
   complex*16      V(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin,PKPTS%nkpoint)
   character*80    fname_header
   character*80    fname
   character*6     kunit_
   character*28    kmode

   fname_header = 'band_structure_TBA'
   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   flag_print_orbital = PINPT%flag_get_orbital
   kpoint = PKPTS%kpoint
   nkpoint= PKPTS%nkpoint
   nbasis = PGEOM%neig
 
   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   call get_e_range(init_e, fina_e, PGEOM, .false., PINPT)
   call get_ispin_print(PINPT%flag_collinear, ispin_print)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

 spin:do is = 1, ispin_print
        call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear) 
        open(pid_energy, file=trim(fname), status = 'unknown')

      eig:do ie =init_e, fina_e
            write(pid_energy, '(2A,I8,A)', ADVANCE = 'yes') kmode,'  energy(eV) :', ie,' -th eigen'     
            if(.not. flag_print_orbital) then
              write(pid_energy,'(A)',ADVANCE='NO')' '
            elseif(  flag_print_orbital) then
              write(pid_energy, '(A)',ADVANCE='NO') '# wavefunction coeff.:                              '
              do im=1,nbasis
                write(pid_energy, '(I9)',ADVANCE='NO')im
              enddo
              write(pid_energy,'(A)')''
            endif

         kp:do ik = 1, nkpoint
              if(flag_klinemode) then
                write(pid_energy,'(1x,F12.6,24x,F12.6,3x)',ADVANCE='NO')kline(ik), E(ie+PGEOM%neig*(is-1),ik)
              elseif(flag_kgridmode) then
                write(pid_energy,'(1x,3F12.6,F12.6,3x)',ADVANCE='NO')kpoint(:,ik), E(ie+PGEOM%neig*(is-1),ik)
              endif
              if(flag_print_orbital) then
          basis:do im=1,nbasis-1
                  if(PINPT%ispinor .eq. 2) then
                    if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO')abs(V(im,ie,ik)) - abs(V(im + nbasis, ie, ik)) ! up - dn : mz
                    else
                      write(pid_energy,'(*(F9.4))',ADVANCE='NO')abs(V(im,ie,ik)) + abs(V(im + nbasis, ie, ik)) ! up + dn : total
                    endif
                  elseif(PINPT%ispinor .eq. 1) then
                    write(pid_energy,'(*(F9.4))',ADVANCE='NO')abs(V(im+PGEOM%neig*(is-1),ie+PGEOM%neig*(is-1),ik))
                  endif
                enddo basis
                if(PINPT%ispinor .eq. 2) then
                  if(PINPT%flag_print_mag .and. PINPT%axis_print_mag .eq. 'mz') then
                    write(pid_energy,'(*(F9.4))',ADVANCE='YES')abs(V(nbasis,ie,ik)) - abs(V(nbasis*2, ie, ik)) ! up - dn : mz
                  else
                    write(pid_energy,'(*(F9.4))',ADVANCE='YES')abs(V(nbasis,ie,ik)) + abs(V(nbasis*2, ie, ik)) ! up + dn : total
                  endif
                elseif(PINPT%ispinor .eq. 1) then
                  write(pid_energy,'(*(F9.4))',ADVANCE='YES')abs(V(nbasis+PGEOM%neig*(is-1),ie+PGEOM%neig*(is-1),ik))
                endif
              endif
              if(.not.flag_print_orbital) write(pid_energy,*)''
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

subroutine get_e_range(init_e, fina_e, PGEOM, flag_ensurf, PINPT)
   use parameters, only : incar, poscar
   implicit none
   type(incar) :: PINPT
   type(poscar):: PGEOM
   integer*4      init_e, fina_e
   logical        flag_ensurf 

   if(.not. flag_ensurf) then
     if(PINPT%flag_erange) then
       init_e = PINPT%init_erange
       fina_e = PINPT%fina_erange
     elseif(.not. PINPT%flag_erange) then
       init_e = 1
       fina_e = PGEOM%neig * PINPT%ispinor
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


