#include "alias.inc"
subroutine replot_dos(PINPT, PGEOM, PKPTS, PRPLT)
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
   real*8                      emin, emax, e_range(PRPLT%replot_dos_nediv)
   real*8                      sigma
   integer*4                   mpierr
   character*4                 kmode
   logical                     flag_vector
   logical                     flag_replot_dos, flag_replot_ldos, flag_replot_sldos
   real*8                      E(PINPT%nband*PINPT%nspin,PKPTS%nkpoint)
   real*8                      V(PGEOM%neig,PINPT%nband*PINPT%nspin,PKPTS%nkpoint) ! not wavevector, 
   real*8                      time2, time1
   character*40                cjob

   if_main call time_check(time2,time1,'init')

   norb             = PGEOM%neig
   nband            = PINPT%nband
   nspin            = PINPT%nspin
   ispinor          = PINPT%ispinor
   ispin            = PINPT%ispin
   nkp              = PKPTS%nkpoint
   flag_replot_dos  = PRPLT%flag_replot_dos
   flag_replot_ldos = PRPLT%flag_replot_ldos
   flag_replot_sldos= PRPLT%flag_replot_sldos
   nediv            = PRPLT%replot_dos_nediv
   emax             = PRPLT%replot_dos_emax
   emin             = PRPLT%replot_dos_emin
   e_range          = emin + (/(k, k=0,nediv-1)/) * (emax - emin)/dble(nediv - 1)
   sigma            = PRPLT%replot_dos_smearing

   allocate( PRPLT%replot_dos_erange(nediv) )
   PRPLT%replot_dos_erange = e_range

   if(PKPTS%flag_klinemode) kmode = 'line'
   if(PKPTS%flag_kgridmode) kmode = 'grid'
   if(flag_replot_ldos .or. flag_replot_sldos) then
     flag_vector = .true.
   elseif(.not. flag_replot_ldos .and. .not. flag_replot_ldos) then
     flag_vector = .false.
   endif
   if_main write(6,*)''
   if_main write(6,*)''
   if_main write(6,'(A)')' ** Program run in REPLOT mode '
   cjob = ' '
   if(flag_replot_dos)   write(cjob,'(A,A)') trim(cjob), ' DOS'
   if(flag_replot_ldos)  write(cjob,'(A,A)') trim(cjob), ' + LDOS'
   if(flag_replot_sldos) write(cjob,'(A,A)') trim(cjob), ' + SLDOS'
   if_main write(6,'(A,A,A)')'    START REPLOT: ',trim(cjob),' EVALUATION'

   if_main call read_energy_tbfit(E, V, nspin, norb, nband, nkp, kmode, flag_vector)
 
#ifdef MPI
   call MPI_BCAST(E, size(E), MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

   if(flag_replot_dos .or. flag_replot_ldos .or. flag_replot_sldos) then
     call get_dos_ldos(PINPT, PGEOM, PRPLT, E, V, nspin, norb, nband, nkp, e_range, nediv, sigma, &
                       flag_replot_dos, flag_replot_ldos, flag_replot_sldos)
   endif   

   if_main call time_check(time2,time1)
   if_main write(6,'(3A,F10.4,A)')'    END REPLOT: ',trim(cjob),' EVALUATION : ',time2, ' (sec)'

   return
endsubroutine

subroutine get_dos_ldos(PINPT, PGEOM, PRPLT, E, V, nspin, norb, nband, nkp, e_range, nediv, sigma, &
                        flag_replot_dos, flag_replot_ldos, flag_replot_sldos)
   use parameters, only : incar, poscar, replot
   use mpi_setup
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
   real*8                      e_range(nediv), sigma
   real*8                      E(PINPT%nband*PINPT%nspin,nkp)
   real*8                      V(PGEOM%neig,PINPT%nband*PINPT%nspin,nkp) ! not wavevector, 
   real*8                      myV(PGEOM%neig,PINPT%nband*PINPT%nspin)
   real*8                      replot_dos_tot(nspin,nediv)
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
  
   ldos_natom = PRPLT%replot_ldos_natom
   ldos_atom  = PRPLT%replot_ldos_atom
   n_orbital  = PGEOM%n_orbital
   iadd       = 10 ; iaddk = 8  ; inck = 1

   allocate(PRPLT%replot_dos_tot(nspin, nediv))
   PRPLT%replot_dos_tot  = 0d0
         replot_dos_tot  = 0d0
   if(flag_replot_ldos) then
     allocate(PRPLT%replot_ldos_tot(PGEOM%max_orb, PRPLT%replot_ldos_natom, nspin, nediv))
     allocate(      replot_ldos_tot(PGEOM%max_orb, PRPLT%replot_ldos_natom, nspin, nediv))
     PRPLT%replot_ldos_tot = 0d0
           replot_ldos_tot = 0d0
   endif
   if(flag_replot_sldos) then
     allocate(PRPLT%replot_sldos_sum(PGEOM%n_atom, nspin, nediv))
     allocate(      replot_sldos_sum(PGEOM%n_atom, nspin, nediv))
     PRPLT%replot_sldos_sum = 0d0
           replot_sldos_sum = 0d0
   endif
   call mpi_job_distribution_chain(nediv, ourjob, ourjob_disp)
   do ik = 1, nkp
     if(nkp .lt. 10) then
       if_main write(6,'(A,I0,A,I0)')  '         STAT KP: ', ik,'/',nkp
     else
       if( ik/real(nkp)*100d0 .ge. real(iaddk*inck) ) then
         if_main write(6,'(A,F10.3,A)')'         STAT KP: ', ik/real(nkp)*100d0, ' %'
         inck = inck + 1
       endif
     endif
#ifdef MPI
     myV = V(:,:,ik)
     call MPI_BCAST(myV, size(myV), MPI_REAL8, 0, mpi_comm_earth, mpierr)
#else
     myV = V(:,:,ik)
#endif
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
!          dos_ = fgauss(sigma, e_range(ie) - E(i+nband*(is-1),ik)) * dkv
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
    write(pid_dos,'(A)')'#    energy (ev)           dos '
  elseif(PINPT%flag_collinear) then
    write(pid_dos,'(A)')'#    energy (ev)          dos-up           dos-dn'
  endif
  do ie = 1, PRPLT%replot_dos_nediv
    if(.not.PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,1x,F16.8)')e_range(ie), PRPLT%replot_dos_tot(1,ie)
    elseif(PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,1x,F16.8)')e_range(ie), PRPLT%replot_dos_tot(1:2,ie)
    endif
  enddo

  close(pid_dos)

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
     write(pid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)            SLDOS(UP)        SLDOS(DN)'
   else
     write(pid_ldos,'(A)'           )'#           CARTESIAN COORDINATE Rx,Ry,Rz(Ang)              SLDOS'
   endif

   do iatom = 1, PGEOM%n_atom
     do ix = 1, nx
       do iy = 1, ny
         do iz = 1, nz
           coord = coord_cart(:,iatom) + (/T(1,1)*(ix-1) + T(1,2)*(iy-1) + T(1,3)*(iz-1), &
                                           T(2,1)*(ix-1) + T(2,2)*(iy-1) + T(2,3)*(iz-1), &
                                           T(3,1)*(iz-1) + T(3,2)*(iz-1) + T(3,3)*(iz-1)/)
           write(pid_ldos,'(3(F16.8,1x),*(F16.8,1x))') coord, sldos(iatom, :)

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
