#include "alias.inc"
subroutine print_berrycurvature(PINPT, PINPT_BERRY, PKPTS, PGEOM, E0, plot_mode)
   use parameters, only  : incar, kpoints, poscar, berry, pid_berrycurv
   use sorting
   type(incar)   :: PINPT
   type(berry)   :: PINPT_berry
   type(kpoints) :: PKPTS
   type(poscar ) :: PGEOM
   integer*4        nerange
   integer*4        is, ie, iee, ik, ik_sort, k
   integer*4        ieig, feig
   integer*4        nkpoint
   integer*4        erange(PGEOM%neig*PINPT%ispin)
   integer*4        isort_index(PKPTS%nkpoint)
   real*8           nerange_
   real*8           omega(PGEOM%neig*PINPT%ispin,3,PKPTS%nkpoint)
   real*8           kpoint(3,PKPTS%nkpoint)
   real*8           kpoint_reci(3,PKPTS%nkpoint)
   real*8           kline(PKPTS%nkpoint)
   real*8           E0(PGEOM%neig*PINPT%ispin,PKPTS%nkpoint)
   character*80     fname_header, fname
   character*6      kunit_
   character*28     kmode  
   character*5      plot_mode
   character*2      spin(2)
   character*1      axis(3)
   logical          flag_klinemode, flag_kgridmode
   axis(1) = 'x';   axis(2) = 'y';   axis(3) = 'z'

   if(PINPT%flag_collinear) then
     spin(1) = 'up'
     spin(2) = 'dn'
   else
     spin(1) = ' '
   endif

   if(plot_mode .eq. 'total') then
     nerange = PGEOM%neig * PINPT%ispin
     erange  = (/(k, k=1,nerange)/)
     fname_header   = 'BERRYCURV_TBA.total'
   elseif(plot_mode .eq. 'range') then
     nerange= PINPT_BERRY%bc_nerange / PINPT%nspin
     erange(1:nerange)  = PINPT_BERRY%bc_erange(1:nerange)
     if(PINPT_BERRY%flag_bc_filenm_provided) then
       write(fname_header,'(A)')trim(PINPT_BERRY%bc_filenm)
     else
       write(fname_header,'(A)')'BERRYCURV_TBA.range'
     endif
   endif

   flag_klinemode = PKPTS%flag_klinemode
   flag_kgridmode = PKPTS%flag_kgridmode
   kpoint = PKPTS%kpoint
   kpoint_reci = PKPTS%kpoint_reci
   nkpoint= PKPTS%nkpoint

   if(flag_kgridmode) call get_sort_index(isort_index,kpoint_reci,size(kpoint_reci(:,1)),nkpoint)

   call get_kunit(PKPTS%kunit, kunit_)
   call get_plotmode(flag_klinemode, flag_kgridmode, kunit_, kmode)
   if(flag_klinemode) call get_kline_dist(kpoint, nkpoint, kline)

   if(plot_mode .eq. 'total') then

  sp:do is = 1, PINPT%nspin
       call get_fname(fname_header, fname,is, PINPT%flag_collinear, PINPT%flag_noncollinear) 
       open(pid_berrycurv, file=trim(fname), status='unknown')
   eig:do ie = 1, nerange ! need to be improved so that ENRANGE is properly applied...
         iee = erange(ie) + (is -1 )*PGEOM%neig
         write(pid_berrycurv, '(2A,I8,2A)') kmode,'  Berrycurvature (A-1):', iee,' -th eigen ',spin(is)
         if(PINPT_BERRY%bc_dimension .eq. 3) then
           write(pid_berrycurv, '(A,36x,A)') '#','  Energy(eV)     Omega(x)     Omega(y)     Omega(z)       kpoint (reciprocal unit)'
         elseif(PINPT_BERRY%bc_dimension .le. 2) then
           write(pid_berrycurv, '(A,36x,3A)') '#','  Energy(eV)     Omega(',axis(PINPT_BERRY%bc_axis(1)),&
                                             ')       kpoint (reciprocal unit)        '
         endif

      kp:do ik = 1, nkpoint
           if(flag_klinemode) then
             if(PINPT_BERRY%bc_dimension .eq. 3) then
               write(pid_berrycurv,'(1x,F12.6,24x,4(E12.4,1x),3F12.4)')kline(ik), E0(iee,ik), &
                                                                       PINPT_BERRY%omega(erange(ie),:,is,ik), &
                                                                       kpoint_reci(:,ik)
             elseif(PINPT_BERRY%bc_dimension .le. 2) then
               write(pid_berrycurv,'(1x,F12.6,24x,2(E12.4,1x),3F12.4)')kline(ik), E0(iee,ik), &
                                                                           PINPT_BERRY%omega(erange(ie),PINPT_BERRY%bc_axis(1),is,ik), &
                                                                           kpoint_reci(:,ik)
             endif
           elseif(flag_kgridmode) then
             ik_sort = isort_index(ik)
             if(PINPT_BERRY%bc_dimension .eq. 3) then
               write(pid_berrycurv,'(1x,3F12.6,4(E12.4,1x),3F12.4)')kpoint(:,ik_sort), E0(iee,ik_sort), &
                                                                    PINPT_BERRY%omega(erange(ie),:,is,ik_sort), &
                                                                    kpoint_reci(:,ik_sort)
             elseif(PINPT_BERRY%bc_dimension .le. 2) then
               write(pid_berrycurv,'(1x,3F12.6,2(E12.4,1x),3F12.4)')kpoint(:,ik_sort), E0(iee,ik_sort), &
                                                                        PINPT_BERRY%omega(erange(ie),PINPT_BERRY%bc_axis(1),is,ik_sort), &
                                                                        kpoint_reci(:,ik_sort)
             endif
           endif
         enddo kp
         write(pid_berrycurv,*)''
         write(pid_berrycurv,*)''
       enddo eig
       close(pid_berrycurv)
     enddo sp

   elseif(plot_mode .eq. 'range') then

 sp_:do is = 1, PINPT%nspin
       call get_fname(fname_header, fname,is, PINPT%flag_collinear, PINPT%flag_noncollinear)
       open(pid_berrycurv, file=trim(fname), status='unknown')
         write(pid_berrycurv, '(5A)') kmode,'  Berrycurvature (A^-1): sum over [', trim(PINPT_BERRY%strip_bc_range),']-th eigen states ',spin(is)
         if(PINPT_BERRY%bc_dimension .eq. 3) then
           write(pid_berrycurv, '(A,36x,A)') '#','    Omega(x)     Omega(y)     Omega(z)       kpoint (reciprocal unit)'
         elseif(PINPT_BERRY%bc_dimension .le. 2) then
           write(pid_berrycurv, '(A,36x,3A)') '#','   Omega(',axis(PINPT_BERRY%bc_axis(1)),')       kpoint (reciprocal unit)        '
         endif
      kp_:do ik = 1, nkpoint
           if(flag_klinemode) then
             if(PINPT_BERRY%bc_dimension .eq. 3) then
               write(pid_berrycurv,'(1x,F12.6,24x,3(E12.4,1x),3F12.4)')kline(ik), &
                                                                       sum(PINPT_BERRY%omega(erange(1:nerange),1,is,ik)), &
                                                                       sum(PINPT_BERRY%omega(erange(1:nerange),2,is,ik)), &
                                                                       sum(PINPT_BERRY%omega(erange(1:nerange),3,is,ik)), &
                                                                       kpoint_reci(:,ik)
             elseif(PINPT_BERRY%bc_dimension .le. 2) then
               write(pid_berrycurv,'(1x,F12.6,24x,1(E12.4,1x),3F12.4)')kline(ik), &
                                                                       sum(PINPT_BERRY%omega(erange(1:nerange),PINPT_BERRY%bc_axis(1),is,ik)), &
                                                                       kpoint_reci(:,ik)
             endif
           elseif(flag_kgridmode) then
             ik_sort = isort_index(ik)
             if(PINPT_BERRY%bc_dimension .eq. 3) then
               write(pid_berrycurv,'(1x,3F12.6,3(E12.4,1x),3F12.4)')kpoint(:,ik_sort), &
                                                                    sum(PINPT_BERRY%omega(erange(1:nerange),1,is,ik_sort)), &
                                                                    sum(PINPT_BERRY%omega(erange(1:nerange),2,is,ik_sort)), &
                                                                    sum(PINPT_BERRY%omega(erange(1:nerange),3,is,ik_sort)), &
                                                                    kpoint_reci(:,ik_sort)
             elseif(PINPT_BERRY%bc_dimension .le. 2) then
               write(pid_berrycurv,'(1x,3F12.6,1(E12.4,1x),3F12.4)')kpoint(:,ik_sort), &
                                                                    sum(PINPT_BERRY%omega(erange(1:nerange),PINPT_BERRY%bc_axis(1),is,ik_sort)), &
                                                                    kpoint_reci(:,ik_sort)
             endif
           endif
         enddo kp_
         write(pid_berrycurv,*)''
         write(pid_berrycurv,*)''
       close(pid_berrycurv)
     enddo sp_

   endif

   return
endsubroutine
