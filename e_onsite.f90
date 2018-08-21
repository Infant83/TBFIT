! this function is for onsite energy of TaS2 system especially...
function e_onsite_xx(ci_orb,cj_orb,PINPT,site_index)
  use parameters, only : rt3, incar
  implicit none
  type (incar  ) :: PINPT
  integer*4         i
  integer*4         e_xp1_center, e_xp2_center
  integer*4         e_xp1_wing1, e_xp2_wing1, e_xp3_wing1, e_xp12_wing1
  integer*4         e_xp1_wing2, e_xp2_wing2, e_xp3_wing2, e_xp12_wing2
  real*8            ptmp1, ptmp2
  real*8            e_onsite_xx
  character*8       ci_orb, cj_orb
  character*20      site_index
  real*8            H(4)

  e_onsite_xx  = 0.d0
  ptmp1        = 0.d0
  ptmp2        = 0.d0

  ! user defined onsite parameter order as defined in PARAM_FIT.dat file.
  ! IMPORTANT: The order should be same as in PARAM_FIT.dat !!!
  e_xp1_center =  1
  e_xp2_center =  2

  e_xp1_wing1  =  3
  e_xp2_wing1  =  4
  e_xp3_wing1  =  5
  e_xp12_wing1 =  6

  e_xp1_wing2  =  7
  e_xp2_wing2  =  8
  e_xp3_wing2  =  9
  e_xp12_wing2 =  10

site:  select case( trim(site_index) )
 
         case ('center')
           ! Be careful, "tuned_onsite" routine only properly works if you have already provide e_tuning parameter in your
           ! PARAMETER file with 19-th variable.
           call tuned_onsite(PINPT, 0, H)
           if( trim(ci_orb) .eq. trim(cj_orb) .and. trim(ci_orb) .eq. 'xp1' ) then
!            ptmp1= PINPT%param(e_xp1_center)
!            if( nint(PINPT%param_const(4, e_xp1_center)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp1_center)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(1)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) )  .and. &
                   ((trim(ci_orb) .eq. 'xp2') .or. (trim(ci_orb) .eq. 'xp3')) ) then
!            ptmp1= PINPT%param(e_xp2_center)
!            if( nint(PINPT%param_const(4, e_xp2_center)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_center)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(2)
           endif

         case ('wing1_a')
           ! Be careful, "tuned_onsite" routine only properly works if you have already provide e_tuning parameter in your
           ! PARAMETER file with 19-th variable.
           call tuned_onsite(PINPT, 1, H)
           if( trim(ci_orb) .eq. trim(cj_orb) .and. trim(ci_orb) .eq. 'xp1' ) then
!            ptmp1= PINPT%param(e_xp1_wing1)
!            if( nint(PINPT%param_const(4, e_xp1_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp1_wing1)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(1)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp2') ) then
!            ptmp1= PINPT%param(e_xp2_wing1)
!            if( nint(PINPT%param_const(4, e_xp2_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing1)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(2)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp3') ) then
!            ptmp1= PINPT%param(e_xp3_wing1)
!            if( nint(PINPT%param_const(4, e_xp3_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp3_wing1)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(3)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp2') .or. &
                    (trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing1)
!            if( nint(PINPT%param_const(4, e_xp12_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing1)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(4)
           endif

         case ('wing1_b')
           call tuned_onsite(PINPT, 1, H)
           if( trim(ci_orb) .eq. trim(cj_orb) .and. trim(ci_orb) .eq. 'xp1' ) then
!            ptmp1= PINPT%param(e_xp1_wing1)
!            if( nint(PINPT%param_const(4, e_xp1_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp1_wing1)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(1)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp2') ) then
!            ptmp1= PINPT%param(e_xp2_wing1)
!            ptmp2= PINPT%param(e_xp3_wing1)
!            if( nint(PINPT%param_const(4, e_xp2_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing1)
!            if( nint(PINPT%param_const(4, e_xp3_wing1)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing1)
!            e_onsite_xx = (ptmp1 + 3.d0*ptmp2)/4.d0
             e_onsite_xx = (H(2) + 3.d0*H(3))/4.d0
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp3') ) then
!            ptmp1= PINPT%param(e_xp2_wing1)
!            ptmp2= PINPT%param(e_xp3_wing1)
!            if( nint(PINPT%param_const(4, e_xp2_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing1)
!            if( nint(PINPT%param_const(4, e_xp3_wing1)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing1)
!            e_onsite_xx = (3.d0*ptmp1 + ptmp2)/4.d0
             e_onsite_xx = (3.d0*H(2) + H(3))/4.d0
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp2') .or. &
                    (trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing1)
!            if( nint(PINPT%param_const(4, e_xp12_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing1)
!            e_onsite_xx = -0.5d0*ptmp1
             e_onsite_xx = -0.5d0*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing1)
!            if( nint(PINPT%param_const(4, e_xp12_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing1)
!            e_onsite_xx = -0.5d0*rt3*ptmp1
             e_onsite_xx = -0.5d0*rt3*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp2')) ) then
!            ptmp1= PINPT%param(e_xp2_wing1)
!            ptmp2= PINPT%param(e_xp3_wing1)
!            if( nint(PINPT%param_const(4, e_xp2_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing1)
!            if( nint(PINPT%param_const(4, e_xp3_wing1)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing1)
!            e_onsite_xx =  rt3/4.d0*(ptmp1 - ptmp2)
             e_onsite_xx =  rt3/4.d0*(H(2) - H(3))
           endif

         case ('wing1_c')
           call tuned_onsite(PINPT, 1, H)
           if( trim(ci_orb) .eq. trim(cj_orb) .and. trim(ci_orb) .eq. 'xp1' ) then
!            ptmp1= PINPT%param(e_xp1_wing1)
!            if( nint(PINPT%param_const(4, e_xp1_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp1_wing1)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(1)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp2') ) then
!            ptmp1= PINPT%param(e_xp2_wing1)
!            ptmp2= PINPT%param(e_xp3_wing1)
!            if( nint(PINPT%param_const(4, e_xp2_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing1)
!            if( nint(PINPT%param_const(4, e_xp3_wing1)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing1)
!            e_onsite_xx = (ptmp1 + 3.d0*ptmp2)/4.d0
             e_onsite_xx = (H(2) + 3.d0*H(3))/4.d0
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp3') ) then
!            ptmp1= PINPT%param(e_xp2_wing1)
!            ptmp2= PINPT%param(e_xp3_wing1)
!            if( nint(PINPT%param_const(4, e_xp2_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing1)
!            if( nint(PINPT%param_const(4, e_xp3_wing1)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing1)
!            e_onsite_xx = (3.d0*ptmp1 + ptmp2)/4.d0
             e_onsite_xx = (3.d0*H(2) + H(3))/4.d0
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp2') .or. &
                    (trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing1)
!            if( nint(PINPT%param_const(4, e_xp12_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing1)
!            e_onsite_xx = -0.5d0*ptmp1
             e_onsite_xx = -0.5d0*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing1)
!            if( nint(PINPT%param_const(4, e_xp12_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing1)
!            e_onsite_xx =  0.5d0*rt3*ptmp1
             e_onsite_xx =  0.5d0*rt3*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp2')) ) then
!            ptmp1= PINPT%param(e_xp2_wing1)
!            ptmp2= PINPT%param(e_xp3_wing1)
!            if( nint(PINPT%param_const(4, e_xp2_wing1)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing1)
!            if( nint(PINPT%param_const(4, e_xp3_wing1)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing1)
!            e_onsite_xx = -rt3/4.d0*(ptmp1 - ptmp2)
             e_onsite_xx = -rt3/4.d0*(H(2) - H(3))
           endif

         case ('wing2_a')
           call tuned_onsite(PINPT, 2, H)
           if( trim(ci_orb) .eq. trim(cj_orb) .and. trim(ci_orb) .eq. 'xp1' ) then
!            ptmp1= PINPT%param(e_xp1_wing2)
!            if( nint(PINPT%param_const(4, e_xp1_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp1_wing2)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(1)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp2') ) then
!            ptmp1= PINPT%param(e_xp2_wing2)
!            if( nint(PINPT%param_const(4, e_xp2_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing2)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(2)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp3') ) then
!            ptmp1= PINPT%param(e_xp3_wing2)
!            if( nint(PINPT%param_const(4, e_xp3_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp3_wing2)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(3)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp2') .or. &
                    (trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing2)
!            if( nint(PINPT%param_const(4, e_xp12_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing2)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(4)
           endif

         case ('wing2_b')
           call tuned_onsite(PINPT, 2, H)
           if( trim(ci_orb) .eq. trim(cj_orb) .and. trim(ci_orb) .eq. 'xp1' ) then
!            ptmp1= PINPT%param(e_xp1_wing2)
!            if( nint(PINPT%param_const(4, e_xp1_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp1_wing2)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(1)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp2') ) then
!            ptmp1= PINPT%param(e_xp2_wing2)
!            ptmp2= PINPT%param(e_xp3_wing2)
!            if( nint(PINPT%param_const(4, e_xp2_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing2)
!            if( nint(PINPT%param_const(4, e_xp3_wing2)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing2)
!            e_onsite_xx = (ptmp1 + 3.d0*ptmp2)/4.d0
             e_onsite_xx = (H(2) + 3.d0*H(3))/4.d0
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp3') ) then
!            ptmp1= PINPT%param(e_xp2_wing2)
!            ptmp2= PINPT%param(e_xp3_wing2)
!            if( nint(PINPT%param_const(4, e_xp2_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing2)
!            if( nint(PINPT%param_const(4, e_xp3_wing2)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing2)
!            e_onsite_xx = (3.d0*ptmp1 + ptmp2)/4.d0
             e_onsite_xx = (3.d0*H(2) + H(3))/4.d0
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp2') .or. &
                    (trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing2)
!            if( nint(PINPT%param_const(4, e_xp12_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing2)
!            e_onsite_xx = -0.5d0*ptmp1
             e_onsite_xx = -0.5d0*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing2)
!            if( nint(PINPT%param_const(4, e_xp12_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing2)
!            e_onsite_xx = -0.5d0*rt3*ptmp1
             e_onsite_xx = -0.5d0*rt3*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp2')) ) then
!            ptmp1= PINPT%param(e_xp2_wing2)
!            ptmp2= PINPT%param(e_xp3_wing2)
!            if( nint(PINPT%param_const(4, e_xp2_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing2)
!            if( nint(PINPT%param_const(4, e_xp3_wing2)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing2)
!            e_onsite_xx =  rt3/4.d0*(ptmp1 - ptmp2)
             e_onsite_xx =  rt3/4.d0*(H(2) - H(3))
           endif

         case ('wing2_c')
           call tuned_onsite(PINPT, 2, H)
           if( trim(ci_orb) .eq. trim(cj_orb) .and. trim(ci_orb) .eq. 'xp1' ) then
!            ptmp1= PINPT%param(e_xp1_wing2)
!            if( nint(PINPT%param_const(4, e_xp1_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp1_wing2)
!            e_onsite_xx = ptmp1
             e_onsite_xx = H(1)
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp2') ) then
!            ptmp1= PINPT%param(e_xp2_wing2)
!            ptmp2= PINPT%param(e_xp3_wing2)
!            if( nint(PINPT%param_const(4, e_xp2_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing2)
!            if( nint(PINPT%param_const(4, e_xp3_wing2)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing2)
!            e_onsite_xx = (ptmp1 + 3.d0*ptmp2)/4.d0
             e_onsite_xx = (H(2) + 3.d0*H(3))/4.d0
           elseif( ( trim(ci_orb) .eq. trim(cj_orb) ) .and. (trim(ci_orb) .eq. 'xp3') ) then
!            ptmp1= PINPT%param(e_xp2_wing2)
!            ptmp2= PINPT%param(e_xp3_wing2)
!            if( nint(PINPT%param_const(4, e_xp2_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing2)
!            if( nint(PINPT%param_const(4, e_xp3_wing2)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing2)
!            e_onsite_xx = (3.d0*ptmp1 + ptmp2)/4.d0
             e_onsite_xx = (3.d0*H(2) + H(3))/4.d0
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp2') .or. &
                    (trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing2)
!            if( nint(PINPT%param_const(4, e_xp12_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing2)
!            e_onsite_xx = -0.5d0*ptmp1
             e_onsite_xx = -0.5d0*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp1') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp1')) ) then
!            ptmp1= PINPT%param(e_xp12_wing2)
!            if( nint(PINPT%param_const(4, e_xp12_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp12_wing2)
!            e_onsite_xx =  0.5d0*rt3*ptmp1
             e_onsite_xx =  0.5d0*rt3*H(4)
           elseif( ( trim(ci_orb) .ne. trim(cj_orb) ) .and. &
                   ((trim(ci_orb) .eq. 'xp2') .and. (trim(cj_orb) .eq. 'xp3') .or. &
                    (trim(ci_orb) .eq. 'xp3') .and. (trim(cj_orb) .eq. 'xp2')) ) then
!            ptmp1= PINPT%param(e_xp2_wing2)
!            ptmp2= PINPT%param(e_xp3_wing2)
!            if( nint(PINPT%param_const(4, e_xp2_wing2)) .eq. 1) ptmp1=PINPT%param_const(5, e_xp2_wing2)
!            if( nint(PINPT%param_const(4, e_xp3_wing2)) .eq. 1) ptmp2=PINPT%param_const(5, e_xp3_wing2)
!            e_onsite_xx = -rt3/4.d0*(ptmp1 - ptmp2)
             e_onsite_xx = -rt3/4.d0*(H(2) - H(3))
           endif

       end select site

return
endfunction

! this particular subroutine is written only for TaS2 system
subroutine tuned_onsite(PINPT, iwing, H)
   use parameters, only : rt2, incar
   implicit none
   type (incar  ) :: PINPT
   real*8            e0, e1, m0, m1
   integer*4         i
   integer*4         iwing
   real*8            H(4)

   select case (iwing)   
     case(0)
       do i = 1, 2
         if(PINPT%param_const(4, i) .eq. 1) then
           H(i) = PINPT%param_const(5, i)
         else
           H(i) = PINPT%param(i)
         endif
       enddo
       e0 = (H(1) + 2d0*H(2))/3d0
       m0 = (H(1) -     H(2))/3d0

       if(PINPT%param_const(4, 19) .eq. 1) then
         e0 = e0 * PINPT%param_const(5,19)
       else
         e0 = e0 * PINPT%param(19)
       endif

       H(1) = e0 + 2d0 * m0
       H(2) = e0 -       m0
       H(3) = 0d0
       H(4) = 0d0

     case(1)
       do i = 3, 6
         if(nint(PINPT%param_const(4, i)) .eq. 1) then
           H(i - 2) = PINPT%param_const(5, i) 
         else
           H(i - 2) = PINPT%param(i)
         endif
       enddo
       e0 = 1.d0/3.d0 * H(1) + 2.d0/3.d0 * H(2)                    + 2.d0*rt2/3.d0 * H(4)
       e1 = 1.d0/3.d0 * H(1) + 1.d0/6.d0 * H(2) + 1.d0/2.d0 * H(3) -      rt2/3.d0 * H(4)
       m0 = 1.d0/3.d0 * H(1) - 1.d0/3.d0 * H(2)                    + 1.d0/rt2/3.d0 * H(4)
       m1 = 1.d0/3.d0 * H(1) + 1.d0/6.d0 * H(2) - 1.d0/2.d0 * H(3) -      rt2/3.d0 * H(4)
       ! apply tunning parameter lambda which is defined in PARAMETER file (19-th parameter)
       if(nint(PINPT%param_const(4, 19)) .eq. 1) then
         e0 = e0 * PINPT%param_const(5,19)
         e1 = e1 * PINPT%param_const(5,19)
       else
         e0 = e0 * PINPT%param(19)
         e1 = e1 * PINPT%param(19)
       endif
       H(1) = 1.d0/3.d0 * (        e0 + 2.d0 * e1 + 4.d0 * m0 + 2.d0 * m1)
       H(2) = 1.d0/3.d0 * ( 2.d0 * e0 +        e1 - 4.d0 * m0 +        m1)
       H(3) =                                  e1 -                    m1
       H(4) = rt2/3.d0  * (        e0 -        e1 +        m0 -        m1)

     case(2)
       do i = 7, 10
         if(nint(PINPT%param_const(4, i)) .eq. 1) then
           H(i - 6) = PINPT%param_const(5, i)
         else
           H(i - 6) = PINPT%param(i)
         endif
       enddo
       e0 = 1d0/3d0 * H(1) + 2d0/3d0 * H(2)                  + 2d0*rt2/3d0 * H(4)
       e1 = 1d0/3d0 * H(1) + 1d0/6d0 * H(2) + 1d0/2d0 * H(3) -     rt2/3d0 * H(4)
       m0 = 1d0/3d0 * H(1) - 1d0/3d0 * H(2)                  + 1d0/rt2/3d0 * H(4)
       m1 = 1d0/3d0 * H(1) + 1d0/6d0 * H(2) - 1d0/2d0 * H(3) -     rt2/3d0 * H(4)
       ! apply tunning parameter lambda which is defined in PARAMETER file (19-th parameter)
       if(nint(PINPT%param_const(4, 19)) .eq. 1) then
         e0 = e0 * PINPT%param_const(5,19)
         e1 = e1 * PINPT%param_const(5,19)
       else
         e0 = e0 * PINPT%param(19)
         e1 = e1 * PINPT%param(19)
       endif
       H(1) = 1d0/3d0 * (       e0 + 2d0 * e1 + 4d0 * m0 + 2d0 * m1)
       H(2) = 1d0/3d0 * ( 2d0 * e0 +       e1 - 4d0 * m0 +       m1)
       H(3) =                              e1 -                  m1
       H(4) = rt2/3d0 * (       e0 -       e1 +       m0 -       m1)

   endselect

   return
endsubroutine
