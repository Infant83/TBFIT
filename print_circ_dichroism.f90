#include "alias.inc"
subroutine print_circ_dichroism(circ_dichroism, E, PINPT, PKPTS, PGEOM)
   use parameters, only: incar, kpoints, pid_circ, poscar
   use sorting
   implicit none
   type(incar)   :: PINPT
   type(kpoints) :: PKPTS
   type(poscar)  :: PGEOM
   integer*4        ik, ii, ieig, feig
   integer*4        ik_sort
   integer*4        nkpoint, iband, fband, is
   real*8           kpoint(3,PKPTS%nkpoint), kpoint_reci(3,PKPTS%nkpoint)
   real*8           kline(PKPTS%nkpoint)
   integer*4        isort_index(PKPTS%nkpoint)
   character*2      spin(2)
   real*8           circ_dichroism(PINPT%nspin, PKPTS%nkpoint, PINPT%ncirc_dichroism)
   real*8           circ_max, circ_min
   integer*4        icirc_max(1), icirc_min(1)
   logical          flag_klinemode, flag_kgridmode
   character*80     fname_header, fname
   character*5      plot_mode
   character*6      kunit
   real*8           E(PGEOM%nband*PINPT%nspin, PKPTS%nkpoint)

   kpoint         = PKPTS%kpoint
   kpoint_reci    = PKPTS%kpoint_reci
   nkpoint        = PKPTS%nkpoint
   spin(1)        = 'no' ; if(PINPT%flag_collinear) spin(:) = (/'up','dn'/)

   if(PKPTS%flag_klinemode) then 
     isort_index = (/(ik, ik=1,nkpoint)/)
   else
     call get_sort_index(isort_index, kpoint_reci, size(kpoint_reci(:,1)),nkpoint)
   endif

  
   call get_kunit(PKPTS%kunit, kunit)
   if(PKPTS%flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

   do ii = 1, PINPT%ncirc_dichroism
     
     do is = 1, PINPT%nspin
       ieig = PINPT%circ_dichroism_pair(1,ii) + (is-1)*PGEOM%nband
       feig = PINPT%circ_dichroism_pair(2,ii) + (is-1)*PGEOM%nband
       write(fname_header,'(A,I0,A,I0)')'CIRC_DICHROISM.EIG_',ieig,'-',feig
       call get_fname(fname_header, fname, is, PINPT%flag_collinear, PINPT%flag_noncollinear)
       icirc_max = maxloc(circ_dichroism(is,:,ii)) ; icirc_min = minloc(circ_dichroism(is,:,ii)) 
        circ_max = maxval(circ_dichroism(is,:,ii)) ;  circ_min = minval(circ_dichroism(is,:,ii))
       call write_header_circ_dichroism(fname, ieig, feig, spin(is),PKPTS%flag_klinemode, kunit, &
                                        kpoint_reci(:,(/icirc_max, icirc_min/)), circ_max, circ_min)

       do ik = 1, nkpoint
         ik_sort = isort_index(ik)
         if(PKPTS%flag_klinemode) then
           write(pid_circ, '(F11.6,22X,3F20.4,5x,3F11.6)') kline(ik_sort)   ,circ_dichroism(is,ik_sort,ii), &
                                                             E(ieig,ik_sort),E(feig,ik_sort), kpoint_reci(:,ik_sort)
         elseif(.not. PKPTS%flag_klinemode) then
           write(pid_circ, '(3F11.6   ,3F20.4,5x,3F11.6)') kpoint(:,ik_sort),circ_dichroism(is,ik_sort,ii), &
                                                             E(ieig,ik_sort),E(feig,ik_sort), kpoint_reci(:,ik_sort)
         endif
       enddo

     enddo ! is

   enddo ! ii
   return
endsubroutine

subroutine write_header_circ_dichroism(fname, ieig, feig, spin, flag_klinemode, kunit, kpoint_reci, circ_max, circ_min)
   use parameters, only: pid_circ
   character*80    fname
   integer*4       ieig, feig
   real*8          circ_max, circ_min
   character*2     spin
   character*6     kunit
   logical         flag_klinemode
   real*8          kpoint_reci(3,2)
   open(pid_circ, file=trim(fname), status='unknown')

   if(spin .eq. 'no') then
     write(pid_circ,'( A,I0,A,I0)')'# OPTICAL SELECTIVITY BETWEEN BANDS V and C: ', ieig,' and ', feig
   else
     write(pid_circ,'(3A,I0,A,I0)')'# OPTICAL SELECTIVITY BETWEEN BANDS V and C (spin-',spin,') : ', ieig, ' and ', feig
   endif

   write(pid_circ,'(A)')"# eta(k,w_cv)= |P(k,s,cv,+)|^2 - |P(k,s,cv,-)|^2"
   write(pid_circ,'(A)')"#              ---------------------------------"
   write(pid_circ,'(A)')"#              |P(k,s,cv,+)|^2 + |P(k,s,cv,-)|^2"
   write(pid_circ,'(A)')"#  The TRANSITION MATRIX ELEMENT P ="
   write(pid_circ,'(A)')"#   P(k,s,cv,+ or -) = 1/sqrt(2)[p_x(k,cv,s) + (or -) i*p_y(k,cv,s)]"
   write(pid_circ,'(A)')"#  THE INTERBAND TRANSITION MATRIX p_x,y ="
   write(pid_circ,'(A)')"#   p_x,y(k,cv,s)=<psi(k,c,s)|-i*hbar*1/dx(y)|psi(k,v,s>"
   write(pid_circ,'(A,3F16.6,A,F16.6)')"# MAXVAL of SELECTIVITY at kx,ky,kz (in reci, ",kpoint_reci(:,1),") = ", circ_max
   write(pid_circ,'(A,3F16.6,A,F16.6)')"# MINVAL of SELECTIVITY at kx,ky,kz (in reci, ",kpoint_reci(:,2),") = ", circ_min

   if(flag_klinemode) then
     write(pid_circ,'(3A)')"# (cart) k-dist ",kunit,"                selectivity(eta(k)),       energ(v)           energy(c)(eV), (recip)kx         ky         kz"
   else                                                                                                        
     write(pid_circ,'(3A)')"# (cart) kx        ky        kz",kunit," selectivity(eta(k)),       energ(v)           energy(c)(eV), (recip)kx         ky         kz"
   endif

   return
endsubroutine
