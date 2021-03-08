#include "alias.inc"
function tij_sk(NN_TABLE,ii,PPRAM,tol,flag_set_overlap)
  use parameters, only : pi, rt2, rt3, hopping, params
  use get_parameter
  use print_io
  use mpi_setup
  implicit none
  type (hopping)  :: NN_TABLE
  type (params )  :: PPRAM
  integer*4 i,ii, iscale_mode
  integer*4 n_nn, i_atom
  integer*4 i_o, i_s, i_p, i_d, i_ovl 
  real*8   l, m, n, ll, mm, nn, lm, ln, mn, lmp, lmm
  real*8   rij(3),dij, d0
  real*8   tij_sk,tol
  real*8   e_(4), s_(4),   p_(4),   d_(4)  ! for NRL SK parameterization
  real*8   e, s,   p,   d                  ! for normal SK parameterization
  real*8      s_s, p_s, d_s
  real*8   i_sign, j_sign
  character*8  ci_orb,cj_orb
  character*20 site_index
  real*8             l_onsite(NN_TABLE%n_nn(NN_TABLE%i_atom(ii)))
  real*8, external:: f_s
  real*8, external:: f_s2
  real*8, external:: f_s_nrl
  real*8, external:: e_onsite_xx
  real*8, external:: e_nrl
!!logical  flag_init, flag_set_overlap
  logical  flag_set_overlap

  ! NN_TABLE%sk_set_index(0:6(+6;if use_overlap),nn)
  !    0        1       2      3       4        5       6
  ! e_onsite  sigma     pi   delta   sigma_s   pi_s   delta_s    ! _s indicates scaling factor
  !             1+6     2+6    3+6     4+6      5+6     6+6     
  !         o_sigma   o_pi o_delta o_sigma_s o_pi_s o_delta_s    ! _s indicates scaling factor, o_ indicates overlap integral

  ! NOTE: if param_class = "xx", param_name for onsite_energy should look like below,
  !       e_'orb_name'_'site_index'
  !       else if param_class = 'ss', 'pp', 'dd', 'sp' etc., i.e., normal slater-koster type,
  !       then the param_name for onsite_energy should look like below,
  !       e_'orb_name'_'species_name'
  !       For more details, see get_param_name routine to see how the 'param_name' is constructed 
  !       in the case of param_class = 'ss', 'pp', 'dd', 'sp' etc, or,
  !       see e_onsite_xx routine to see how the 'param_name' is constructed
  !       in the case of param_class = 'xx' which is based on the 'site_index' subtag.

  i_ovl = 0
  if(flag_set_overlap) i_ovl = 6

  iscale_mode= PPRAM%slater_koster_type ! default = 1, see 'f_s' function for the detail.
  i_atom     = NN_TABLE%i_atom(ii)
  rij(1:3)   = NN_TABLE%Rij(1:3,ii)
  dij        = NN_TABLE%Dij(    ii)
  ci_orb     = NN_TABLE%ci_orb( ii)
  cj_orb     = NN_TABLE%cj_orb( ii)
   i_sign    = NN_TABLE%i_sign( ii)
   j_sign    = NN_TABLE%j_sign( ii)
  site_index = NN_TABLE%site_cindex( i_atom)
  if( NN_TABLE%n_class(ii) .gt. 0 ) then
    l        = rij(1)/dij
    m        = rij(2)/dij
    n        = rij(3)/dij
    ll       = l*l
    mm       = m*m
    nn       = n*n
    lm       = l*m
    ln       = l*n
    mn       = m*n
    lmp      = ll + mm
    lmm      = ll - mm
  endif


  if(flag_set_overlap .and. NN_TABLE%n_class(ii)  .eq.                    0    .and.  &
                            NN_TABLE%i_matrix(ii) .eq. NN_TABLE%j_matrix(ii)  ) then
    e  = 1d0
    if(iscale_mode .gt. 10) e_(1:4)= (/1d0,0d0,0d0,0d0/)
  else
    e  = 0d0
    if(iscale_mode .gt. 10) e_ = 0d0
  endif
  s  = 0d0
  p  = 0d0
  d  = 0d0
  s_s= 1d0 !default if not provided
  p_s= 1d0 !default if not provided
  d_s= 1d0 !default if not provided
  d0  = NN_TABLE%Dij0(ii)
  if(iscale_mode .gt. 10 ) then
    n_nn      = NN_TABLE%n_nn(i_atom)
  endif

  if( NN_TABLE%n_class(ii) .eq. 0 ) then

      if( NN_TABLE%p_class(ii) .eq. 'xx' ) then  
        ! set user defined onsite energy modifications in 'e_onsite_xx' function
        e = e_onsite_xx(ci_orb,cj_orb, PPRAM, site_index)

      else
        ! set onsite energy species by species and orbital by orbital
        i_o = NN_TABLE%sk_index_set(0,ii) ! onsite energy parameter index
        if(i_o .gt. 0 .and. .not. flag_set_overlap) then ! if onsite energy for Hk
          if(iscale_mode .le. 10) call get_param(PPRAM, i_o, 1, e)
          if(iscale_mode .gt. 10) then 
            do i = 1, 4
              call get_param(PPRAM,i_o, i, e_(i)) ! get onsite energy parameters
            enddo
            do i = 1, n_nn
              call get_param(PPRAM,NN_TABLE%l_onsite_param_index(NN_TABLE%j_nn(i,i_atom)),1, l_onsite(i)) ! get l_onsite for each j_nn
            enddo
          endif
        else
          l_onsite = 0d0
        endif
      endif ! set_onsite

  elseif(NN_TABLE%n_class(ii) .gt. 0) then
      if(iscale_mode .gt. 10) then
        i_s = NN_TABLE%sk_index_set(1+i_ovl,ii)
        i_p = NN_TABLE%sk_index_set(2+i_ovl,ii)
        i_d = NN_TABLE%sk_index_set(3+i_ovl,ii)
        if( i_s .ne. 0) then
          do i = 1, 4
            call get_param(PPRAM, i_s, i , s_(i) ) ! sigma
          enddo
        endif
        if( i_p .ne. 0) then
          do i = 1, 4
            call get_param(PPRAM, i_p, i, p_(i) ) ! pi
          enddo
        endif
        if( i_d .ne. 0) then 
          do i = 1, 4
            call get_param(PPRAM, i_d, i, d_(i) ) ! delta
          enddo
        endif
      else
        i_s = NN_TABLE%sk_index_set(1+i_ovl,ii)
        i_p = NN_TABLE%sk_index_set(2+i_ovl,ii)
        i_d = NN_TABLE%sk_index_set(3+i_ovl,ii)
        if( i_s .ne. 0) call get_param(PPRAM, i_s, 1, s  ) ! sigma
        if( i_p .ne. 0) call get_param(PPRAM, i_p, 1, p  ) ! pi
        if( i_d .ne. 0) call get_param(PPRAM, i_d, 1, d  ) ! delta
        i_s = NN_TABLE%sk_index_set(4+i_ovl,ii) 
        i_p = NN_TABLE%sk_index_set(5+i_ovl,ii)
        i_d = NN_TABLE%sk_index_set(6+i_ovl,ii)
        if( i_s .ne. 0) call get_param(PPRAM, i_s, 1, s_s) ! sigma_scale
        if( i_p .ne. 0) call get_param(PPRAM, i_p, 1, p_s) ! pi_scale
        if( i_d .ne. 0) call get_param(PPRAM, i_d, 1, d_s) ! delta_scale
      endif

  endif !check n_class


  if( iscale_mode .gt. 10 .and. NN_TABLE%n_class(ii) .gt. 0) then
    s = f_s_nrl ( s_, d0, dij, iscale_mode, PPRAM%l_broaden) * i_sign * j_sign
    p = f_s_nrl ( p_, d0, dij, iscale_mode, PPRAM%l_broaden) * i_sign * j_sign
    d = f_s_nrl ( d_, d0, dij, iscale_mode, PPRAM%l_broaden) * i_sign * j_sign

  elseif(iscale_mode .le. 10 .and. NN_TABLE%n_class(ii) .gt. 0) then

    s = s * f_s ( s_s, d0, dij, iscale_mode) * i_sign * j_sign
    p = p * f_s ( p_s, d0, dij, iscale_mode) * i_sign * j_sign
    d = d * f_s ( d_s, d0, dij, iscale_mode) * i_sign * j_sign

  endif

!endif
  ! SK-energy integral if nn_class > 0
  if( NN_TABLE%n_class(ii) .ne. 0) then
sk: select case ( NN_TABLE%p_class(ii) )

      case ('ss')
        if    (ci_orb(1:1) .eq. 's' .and. cj_orb(1:1) .eq. 's') then
          tij_sk = s
        else
          write(message,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.' ; write_msg
          write(message,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb) ; write_msg
          write(message,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...' ; write_msg
          stop
        endif

      case ('pp') 
        if    (ci_orb(1:2) .eq. 'px' .and. cj_orb(1:2) .eq. 'px') then
          tij_sk = ll*s + (1d0-ll)*p 

        elseif(ci_orb(1:2) .eq. 'py' .and. cj_orb(1:2) .eq. 'py') then
          tij_sk = mm*s + (1d0-mm)*p 

        elseif(ci_orb(1:2) .eq. 'pz' .and. cj_orb(1:2) .eq. 'pz') then
          tij_sk = nn*s + (1d0-nn)*p 

        elseif( (ci_orb(1:2) .eq. 'px' .and. cj_orb(1:2) .eq. 'py') .or. &
                (ci_orb(1:2) .eq. 'py' .and. cj_orb(1:2) .eq. 'px') ) then
          tij_sk = lm*(s - p)

        elseif( (ci_orb(1:2) .eq. 'px' .and. cj_orb(1:2) .eq. 'pz') .or. &
                (ci_orb(1:2) .eq. 'pz' .and. cj_orb(1:2) .eq. 'px') ) then
          tij_sk = ln*(s - p)

        elseif( (ci_orb(1:2) .eq. 'py' .and. cj_orb(1:2) .eq. 'pz') .or. &
                (ci_orb(1:2) .eq. 'pz' .and. cj_orb(1:2) .eq. 'py') ) then
          tij_sk = mn*(s - p)
     
        else
          write(message,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.' ; write_msg
          write(message,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb) ; write_msg
          write(message,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...' ; write_msg
          stop
        endif
  
      case ('dd')
        if    (ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:3) .eq. 'dz2') then
          tij_sk = ((nn-0.5d0*lmp)**2)*s + 3d0*nn*lmp*p + 0.75d0*(lmp**2)*d

        elseif(ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dx2') then
          tij_sk = 0.75d0*(lmm**2)*s + (lmp-(lmm**2))*p + (nn+0.25d0*(lmm**2))*d

        elseif(ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dxy') then
          tij_sk = 3d0*(lm**2)*s + (lmp-4d0*(lm**2))*p + (nn+(lm**2))*d

        elseif(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dxz') then
          tij_sk = 3d0*ll*nn*s + (ll+nn - 4d0*ll*nn)*p + (mm+ll*nn)*d 

        elseif(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dyz') then
          tij_sk = 3d0*mm*nn*s + (mm+nn - 4d0*mm*nn)*p + (ll+mm*nn)*d

        elseif( (ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:3) .eq. 'dx2') .or. &
                (ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dz2') ) then
          tij_sk = sin(pi/3d0)*lmm*(nn-0.5d0*lmp)*s - 2d0*sin(pi/3d0)*nn*lmm*p + sin(pi/3d0)/2d0*(1d0+nn)*lmm*d

        elseif( (ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:3) .eq. 'dxy') .or. &
                (ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dz2') ) then
          tij_sk = 2d0*sin(pi/3d0)*lm*(nn-0.5d0*lmp)*s - 4d0*sin(pi/3d0)*lm*nn*p + sin(pi/3d0)*lm*(1d0+nn)*d

        elseif( (ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:3) .eq. 'dxz') .or. &
                (ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dz2') ) then
          tij_sk = 2d0*sin(pi/3d0)*ln*(nn-0.5d0*lmp)*s + 2d0*sin(pi/3d0)*ln*(lmp-nn)*p - sin(pi/3d0)*ln*lmp*d

        elseif( (ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:3) .eq. 'dyz') .or. &
                (ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dz2') ) then
          tij_sk = 2d0*sin(pi/3)*mn*(nn-0.5d0*lmp)*s + 2d0*sin(pi/3d0)*mn*(lmp-nn)*p - sin(pi/3d0)*mn*lmp*d

        elseif( (ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dxy') .or. &
                (ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dx2') ) then
          tij_sk = 1.5d0*lm*lmm*s - 2d0*lm*lmm*p + 0.5d0*lm*lmm*d

        elseif( (ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dxz') .or. &
                (ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dx2') ) then
          tij_sk = 1.5d0*ln*lmm*s + ln*(1d0-2d0*lmm)*p - ln*(1d0-0.5d0*lmm)*d

        elseif( (ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dyz') .or. &
                (ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dx2') ) then
          tij_sk = 1.5d0*mn*lmm*s - mn*(1d0+2d0*lmm)*p + mn*(1d0+0.5d0*lmm)*d

        elseif( (ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dxz') .or. &
                (ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dxy') ) then
          tij_sk = 3d0*ll*mn*s + mn*(1d0-4d0*ll)*p + mn*(ll-1d0)*d

        elseif( (ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dyz') .or. &
                (ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dxy') ) then
          tij_sk = 3d0*ln*mm*s + ln*(1d0-4d0*mm)*p + ln*(mm-1d0)*d

        elseif( (ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dyz') .or. &
                (ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dxz') ) then
          tij_sk = 3d0*lm*nn*s + lm*(1d0-4d0*nn)*p + lm*(nn-1d0)*d

        else
          write(message,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.' ; write_msg
          write(message,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb) ; write_msg
          write(message,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...' ; write_msg
          stop
        endif

      case ('sp', 'ps')
        if    (ci_orb(1:1) .eq. 's' .and. cj_orb(1:2) .eq. 'px') then
          tij_sk = l*s

        elseif(ci_orb(1:2) .eq. 'px' .and. cj_orb(1:1) .eq. 's') then
          tij_sk =-l*s

        elseif(ci_orb(1:1) .eq. 's' .and. cj_orb(1:2) .eq. 'py') then
          tij_sk = m*s

        elseif(ci_orb(1:2) .eq. 'py' .and. cj_orb(1:1) .eq. 's') then
          tij_sk =-m*s

        elseif(ci_orb(1:1) .eq. 's' .and. cj_orb(1:2) .eq. 'pz') then
          tij_sk = n*s

        elseif(ci_orb(1:2) .eq. 'pz' .and. cj_orb(1:1) .eq. 's') then
          tij_sk =-n*s

        else
          write(message,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.' ; write_msg
          write(message,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb) ; write_msg
          write(message,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...' ; write_msg
          stop
        endif

      case ('sd', 'ds')
        if    ( (ci_orb(1:1) .eq. 's'   .and. cj_orb(1:3) .eq. 'dz2') .or. &
                (ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:1) .eq.  's' ) ) then
          tij_sk = (nn-0.5d0*lmp)*s

        elseif( (ci_orb(1:1) .eq. 's'   .and. cj_orb(1:3) .eq. 'dx2') .or. &
                (ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:1) .eq.  's' ) ) then
          tij_sk = sin(pi/3d0)*lmm*s

        elseif( (ci_orb(1:1) .eq. 's'   .and. cj_orb(1:3) .eq. 'dxy') .or. &
                (ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:1) .eq.  's' ) ) then
          tij_sk = 2d0*sin(pi/3d0)*lm*s

        elseif( (ci_orb(1:1) .eq. 's'   .and. cj_orb(1:3) .eq. 'dxz') .or. &
                (ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:1) .eq.  's' ) ) then
          tij_sk = 2d0*sin(pi/3d0)*ln*s

        elseif( (ci_orb(1:1) .eq. 's'   .and. cj_orb(1:3) .eq. 'dyz') .or. &
                (ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:1) .eq.  's' ) ) then
          tij_sk = 2d0*sin(pi/3d0)*mn*s

        else
          write(message,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.' ; write_msg
          write(message,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb) ; write_msg
          write(message,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...' ; write_msg
          stop
        endif

      case ('pd', 'dp')
        if    (ci_orb(1:2) .eq. 'px'  .and. cj_orb(1:3) .eq. 'dz2' ) then
          tij_sk = l*(nn-0.5d0*lmp)*s - 2d0*sin(pi/3d0)*l*nn*p
        elseif(ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:2) .eq. 'px'  ) then
          tij_sk =-l*(nn-0.5d0*lmp)*s + 2d0*sin(pi/3d0)*l*nn*p

        elseif(ci_orb(1:2) .eq. 'px'  .and. cj_orb(1:3) .eq. 'dx2' ) then
          tij_sk = sin(pi/3d0)*l*lmm*s + l*(1d0-lmm)*p
        elseif(ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:2) .eq. 'px'  ) then
          tij_sk =-sin(pi/3d0)*l*lmm*s - l*(1d0-lmm)*p

        elseif(ci_orb(1:2) .eq. 'px'  .and. cj_orb(1:3) .eq. 'dxy' ) then
          tij_sk = 2d0*sin(pi/3d0)*ll*m*s + m*(1d0-2d0*ll)*p 
        elseif(ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:2) .eq. 'px'  ) then
          tij_sk =-2d0*sin(pi/3d0)*ll*m*s - m*(1d0-2d0*ll)*p

        elseif(ci_orb(1:2) .eq. 'px'  .and. cj_orb(1:3) .eq. 'dxz' ) then
          tij_sk = 2d0*sin(pi/3d0)*ll*n*s + n*(1d0-2d0*ll)*p
        elseif(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:2) .eq. 'px'  ) then
          tij_sk =-2d0*sin(pi/3d0)*ll*n*s - n*(1d0-2d0*ll)*p

        elseif(ci_orb(1:2) .eq. 'px'  .and. cj_orb(1:3) .eq. 'dyz' ) then
          tij_sk = 2d0*sin(pi/3d0)*lm*n*s - 2d0*lm*n*p
        elseif(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:2) .eq. 'px'  ) then
          tij_sk =-2d0*sin(pi/3d0)*lm*n*s + 2d0*lm*n*p

        elseif(ci_orb(1:2) .eq. 'py'  .and. cj_orb(1:3) .eq. 'dz2' ) then
          tij_sk = m*(nn-0.5d0*lmp)*s - 2d0*sin(pi/3d0)*m*nn*p
        elseif(ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:2) .eq. 'py'  ) then
          tij_sk =-m*(nn-0.5d0*lmp)*s + 2d0*sin(pi/3d0)*m*nn*p

        elseif(ci_orb(1:2) .eq. 'py'  .and. cj_orb(1:3) .eq. 'dx2' ) then
          tij_sk = sin(pi/3d0)*m*lmm*s - m*(1d0+lmm)*p 
        elseif(ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:2) .eq. 'py'  ) then
          tij_sk =-sin(pi/3d0)*m*lmm*s + m*(1d0+lmm)*p

        elseif(ci_orb(1:2) .eq. 'py'  .and. cj_orb(1:3) .eq. 'dxy' ) then
          tij_sk = 2d0*sin(pi/3d0)*mm*l*s + l*(1d0-2d0*mm)*p 
        elseif(ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:2) .eq. 'py'  ) then
          tij_sk =-2d0*sin(pi/3d0)*mm*l*s - l*(1d0-2d0*mm)*p

        elseif(ci_orb(1:2) .eq. 'py'  .and. cj_orb(1:3) .eq. 'dxz' ) then
          tij_sk = 2d0*sin(pi/3d0)*mn*l*s - 2d0*mn*l*p
        elseif(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:2) .eq. 'py'  ) then
          tij_sk =-2d0*sin(pi/3d0)*mn*l*s + 2d0*mn*l*p

        elseif(ci_orb(1:2) .eq. 'py'  .and. cj_orb(1:3) .eq. 'dyz' ) then
          tij_sk = 2d0*sin(pi/3d0)*mm*n*s + n*(1d0-2d0*mm)*p 
        elseif(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:2) .eq. 'py'  ) then
          tij_sk =-2d0*sin(pi/3d0)*mm*n*s - n*(1d0-2d0*mm)*p

        elseif(ci_orb(1:2) .eq. 'pz'  .and. cj_orb(1:3) .eq. 'dz2' ) then
          tij_sk = n*(nn-0.5d0*lmp)*s + 2d0*sin(pi/3d0)*n*lmp*p
        elseif(ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:2) .eq. 'pz'  ) then
          tij_sk =-n*(nn-0.5d0*lmp)*s - 2d0*sin(pi/3d0)*n*lmp*p

        elseif(ci_orb(1:2) .eq. 'pz'  .and. cj_orb(1:3) .eq. 'dx2' ) then
          tij_sk = sin(pi/3d0)*n*lmm*s - n*lmm*p
        elseif(ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:2) .eq. 'pz'  ) then
          tij_sk =-sin(pi/3d0)*n*lmm*s + n*lmm*p

        elseif(ci_orb(1:2) .eq. 'pz'  .and. cj_orb(1:3) .eq. 'dxy' ) then
          tij_sk = 2d0*sin(pi/3d0)*ln*m*s - 2d0*ln*m*p
        elseif(ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:2) .eq. 'pz'  ) then
          tij_sk =-2d0*sin(pi/3d0)*ln*m*s + 2d0*ln*m*p

        elseif(ci_orb(1:2) .eq. 'pz'  .and. cj_orb(1:3) .eq. 'dxz' ) then
          tij_sk = 2d0*sin(pi/3d0)*nn*l*s + l*(1d0-2d0*nn)*p
        elseif(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:2) .eq. 'pz'  ) then
          tij_sk =-2d0*sin(pi/3d0)*nn*l*s - l*(1d0-2d0*nn)*p

        elseif(ci_orb(1:2) .eq. 'pz'  .and. cj_orb(1:3) .eq. 'dyz' ) then
          tij_sk = 2d0*sin(pi/3d0)*nn*m*s + m*(1d0-2d0*nn)*p
        elseif(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:2) .eq. 'pz'  ) then
          tij_sk =-2d0*sin(pi/3d0)*nn*m*s - m*(1d0-2d0*nn)*p

        else
          write(message,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.' ; write_msg
          write(message,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb) ; write_msg
          write(message,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...' ; write_msg
          stop
        endif

      ! set up interatomic hopping parameter for the customized (user defined) function
      case ('xx')
        ! xp1 = phi1 = dz2
        if    (ci_orb(1:3) .eq. 'xp1' .and. cj_orb(1:3) .eq. 'xp1') then
          tij_sk = ((nn-0.5d0*lmp)**2)*s + 3d0*nn*lmp*p + 0.75d0*(lmp**2)*d

        ! xp2 = -2/rt6 dx2 + 1/rt3 dyz
        elseif(ci_orb(1:3) .eq. 'xp2' .and. cj_orb(1:3) .eq. 'xp2') then
          tij_sk = 2.d0/3.d0 * ( 0.75d0*(lmm**2)*s + (lmp-(lmm**2))*p + (nn+0.25d0*(lmm**2))*d ) &
               +1.d0/3.d0 * ( 3d0*mm*nn*s + (mm+nn - 4d0*mm*nn)*p + (ll+mm*nn)*d ) &
               -2.d0/3.d0 * rt2 * ( 1.5d0*mn*lmm*s - mn*(1d0+2d0*lmm)*p + mn*(1d0+0.5d0*lmm)*d )

        ! xp3 = -2/rt6 dxy + 1/rt3 dxz
        elseif(ci_orb(1:3) .eq. 'xp3' .and. cj_orb(1:3) .eq. 'xp3') then
          tij_sk = 2.d0/3.d0 * ( 3d0*(lm**2)*s + (lmp-4d0*(lm**2))*p + (nn+(lm**2))*d ) &
               +1.d0/3.d0 * ( 3d0*ll*nn*s + (ll+nn - 4d0*ll*nn)*p + (mm+ll*nn)*d ) &
               -2.d0/3.d0 * rt2 * ( 3d0*ll*mn*s + mn*(1d0-4d0*ll)*p + mn*(ll-1d0)*d )

        elseif( (ci_orb(1:3) .eq. 'xp1' .and. cj_orb(1:3) .eq. 'xp2') .or. &
                (ci_orb(1:3) .eq. 'xp2' .and. cj_orb(1:3) .eq. 'xp1') ) then
          tij_sk =-1.d0/rt3 * ( rt2 * ( sin(pi/3d0)*lmm*(nn-0.5d0*lmp)*s - 2d0*sin(pi/3d0)*nn*lmm*p + sin(pi/3d0)/2d0*(1d0+nn)*lmm*d ) &
                            -      ( 2d0*sin(pi/3)*mn*(nn-0.5d0*lmp)*s + 2d0*sin(pi/3d0)*mn*(lmp-nn)*p - sin(pi/3d0)*mn*lmp*d ) )

        elseif( (ci_orb(1:3) .eq. 'xp1' .and. cj_orb(1:3) .eq. 'xp3') .or. &
                (ci_orb(1:3) .eq. 'xp3' .and. cj_orb(1:3) .eq. 'xp1') ) then
          tij_sk = 1.d0/rt3 * (-rt2 * ( 2d0*sin(pi/3d0)*lm*(nn-0.5d0*lmp)*s - 4d0*sin(pi/3d0)*lm*nn*p + sin(pi/3d0)*lm*(1d0+nn)*d ) &
                            +      ( 2d0*sin(pi/3d0)*ln*(nn-0.5d0*lmp)*s + 2d0*sin(pi/3d0)*ln*(lmp-nn)*p - sin(pi/3d0)*ln*lmp*d ) )

        elseif( (ci_orb(1:3) .eq. 'xp2' .and. cj_orb(1:3) .eq. 'xp3') .or. &
                (ci_orb(1:3) .eq. 'xp3' .and. cj_orb(1:3) .eq. 'xp2') ) then
          tij_sk = 2.d0/3.d0 * ( 1.5d0*lm*lmm*s - 2d0*lm*lmm*p + 0.5d0*lm*lmm*d ) &
               -1.d0/3.d0 * rt2 * ( 1.5d0*ln*lmm*s + ln*(1d0-2d0*lmm)*p - ln*(1d0-0.5d0*lmm)*d ) &
               -1.d0/3.d0 * rt2 * ( 3d0*ln*mm*s + ln*(1d0-4d0*mm)*p + ln*(mm-1d0)*d ) &
               +1.d0/3.d0 *       ( 3d0*lm*nn*s + lm*(1d0-4d0*nn)*p + lm*(nn-1d0)*d )

        else
          write(message,'(A)')' !WARNING! SK_C energy integral is not properly defined or orbital name is improper.' ; write_msg
          write(message,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb) ; write_msg
          write(message,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...' ; write_msg
          stop
        endif

    end select sk

  ! onsite energy if nn_class == 0
  elseif( NN_TABLE%n_class(ii) .eq. 0 ) then 
    if(NN_TABLE%flag_efield .and. (NN_TABLE%i_matrix(ii) .eq. NN_TABLE%j_matrix(ii)) ) then ! if E-field
      if(.not.flag_set_overlap) then

        if(iscale_mode .gt. 10) then
          if(n_nn .gt. 0) then
            tij_sk = e_nrl(e_, NN_TABLE%R0_nn(1:n_nn,i_atom), NN_TABLE%R_nn(1:n_nn,i_atom), n_nn, l_onsite, PPRAM%l_broaden) &
                    -dot_product(NN_TABLE%efield(1:3),NN_TABLE%i_coord(1:3,ii)-NN_TABLE%efield_origin_cart(1:3))
          elseif(n_nn .eq. 0) then
            tij_sk = e_(1)
          endif
        else
          tij_sk = e - dot_product(NN_TABLE%efield(1:3),NN_TABLE%i_coord(1:3,ii)-NN_TABLE%efield_origin_cart(1:3))
        endif
      elseif(flag_set_overlap) then
        tij_sk = e 
      endif
    else ! if not E-field
      if(iscale_mode .gt. 10) then
        if(n_nn .gt. 0 .and. .not. flag_set_overlap) then
          tij_sk = e_nrl(e_, NN_TABLE%R0_nn(1:n_nn,i_atom), NN_TABLE%R_nn(1:n_nn,i_atom), n_nn, l_onsite, PPRAM%l_broaden) 
        elseif(n_nn .eq. 0) then
          tij_sk = e_(1)
        endif
      else
        tij_sk = e 
      endif
    endif

  endif

return
endfunction
function e_nrl(e, R0, R_nn, n_nn, l_onsite, l_broaden)
   implicit none
   integer*4    i, n_nn
   real*8       e(4)
   real*8       l_onsite(n_nn), l_broaden
   real*8       R_nn(n_nn), R0(n_nn)
   real*8       f_cut(n_nn)
   real*8       e_nrl, rho_at

   ! rho : local atomic density at atom i with additional parameter ldensity
   rho_at = 0d0
   f_cut = 1d0/(1d0 + Exp( (R_nn(:) - R0(:))/l_broaden + 5d0 ) )
   rho_at= sum(Exp(-(l_onsite(:)**2d0) * R_nn(:) ) * f_cut(:))

   e_nrl = e(1) + e(2) * (rho_at**(2d0/3d0)) + &
                  e(3) * (rho_at**(4d0/3d0)) + &
                  e(4) * (rho_at**(2d0    ))
return
endfunction
function f_s_nrl(dda, d0, d, mode, l_broaden)
   implicit none
   integer*4    mode
   real*8       f_s_nrl, l_broaden, f_cut
   real*8       dda(4), d0, d
   
!    mode : scaling function
!     11 = see Ref. PRB 54, 4519 (1996) : Naval Resarch Laboratory (NRL) TB scheme
   if(mode .eq. 11) then
     f_cut   = (1d0 + Exp( (d - d0)/l_broaden + 5d0 ) )**(-1d0)
     f_s_nrl = ( dda(1) + dda(2) * d + dda(3) * (d**2d0) ) * Exp(-(dda(4)**2d0) * d) * f_cut
   endif

return
endfunction
function f_s(dda_s,d0,d, mode)
   implicit none
   integer*4 mode
   real*8    f_s
   real*8    dda_s,d0,d

! mode : scaling function
!  1 = see Ref. PRB 85.195458 (2012): for interlayer pz-pz interaction of twisted BL graphene
!    f_s=exp( (-D+D0) / (SFACTOR(1)) )
!  2 = see Ref. PRB 92.205108 (2015): for interlayer p-p interaction of layered TMDC material
!    f_s=exp( -(d/d0)^(SFACTOR(1)) )
!  3 = see PRB 51.16772 (1995): for s-p or p-p interaction of Silicon or Germanium crystal
!    f_s=(d0/d)^(SFACTOR(1))
!  4 = see PRB 93.241407 (2016) : In-Si case
!    f_s=exp( (d0-d) * SFACTOR(1) )
!  6 = see Europhys.Lett.9.701 (1989): GSP parameterization for carbon system
!      SFACTOR(1) = m ; SFACTOR(2) = d_c , critical distance ; SFACTOR(3) = m_c
!    f_s=(d0/d)^(SFACTOR(1))*exp(SFACTOR(1)*( -(d/SFACTOR(2))^(SFACTOR(3)) + -(d0/SFACTOR(2))^(SFACTOR(3)) ))

   if(mode .eq. 1) then
     f_s = Exp( (d0 - d)/(dda_s*d0) )
!    f_s = Exp( (d0 - d)/(dda_s) )
   elseif(mode .eq. 2) then
     f_s = Exp( -(d/d0)**dda_s )
   elseif(mode .eq. 3) then
     f_s = (d0/d)**(dda_s)
   elseif(mode .eq. 4) then
    !f_s = Exp( (d0 - d)*dda_s )
     f_s = Exp( -abs(dda_s) * (d - d0) )
   elseif(mode .eq. 5) then
     f_s = (1 - dda_s * (d - d0) )

!  elseif(mode .eq. 6)
!    f_s = (d0/d)**( dda_s(1) ) * Exp( dda_s(1) * ( -(d/dda_s(2))**dda_s(3) - (d0/dda_s(2))**dda_s(3) ) )
   endif
return
endfunction

! deprecated.. HJK, 03.Dec.2020
function f_s2(dda_s, dda_s2, d0, d, mode)
   implicit none
   integer*4 mode
   real*8    f_s2
   real*8    dda_s, dda_s2, d0,d

! mode : scaling function
! 11 = see Ref. PRB 85.195458 (2012): for interlayer pz-pz interaction of twisted BL graphene
!    f_s=exp( (-D+D0) / (SFACTOR(1)) )
! 12 = see Ref. PRB 92.205108 (2015): for interlayer p-p interaction of layered TMDC material
!    f_s=exp( -(d/d0)^(SFACTOR(1)) )
! 13 = see PRB 51.16772 (1995): for s-p or p-p interaction of Silicon or Germanium crystal
!    f_s=(d0/d)^(SFACTOR(1))
! 14 = see Europhys.Lett.9.701 (1989): GSP parameterization for carbon system
!      SFACTOR(1) = m ; SFACTOR(2) = d_c , critical distance ; SFACTOR(3) = m_c
!    f_s=(d0/d)^(SFACTOR(1))*exp(SFACTOR(1)*( -(d/SFACTOR(2))^(SFACTOR(3)) + -(d0/SFACTOR(2))^(SFACTOR(3)) ))

   if(mode .eq. 11) then
     f_s2 = Exp( (d0 - d)/(dda_s*d0) ) * dda_s2
   elseif(mode .eq. 12) then
     f_s2 = Exp( -(d/d0)**dda_s )
   elseif(mode .eq. 13) then
     f_s2 = (d0/d)**(dda_s)
!  elseif(mode .eq. 14)
!    f_s2 = (d0/d)**( dda_s(1) ) * Exp( dda_s(1) * ( -(d/dda_s(2))**dda_s(3) - (d0/dda_s(2))**dda_s(3) ) )
   endif
return
endfunction

subroutine get_hopping_integral(Eij, NN_TABLE, nn, PPRAM, tol, kpoint, FIJ_, &
                                flag_phase, flag_set_overlap, flag_load_nntable)
  use parameters, only : zi, hopping, params
  use phase_factor
  implicit none
  interface
    function FIJ_(k,R)
      complex*16   FIJ_
      real*8, intent(in) :: k(3)
      real*8, intent(in) :: R(3)
    endfunction
  end interface
  type (hopping) :: NN_TABLE
  type (params ) :: PPRAM
  integer*4         nn
  real*8            kpoint(3), tol, tij_sk, tij_cc
  real*8            Rij(3)
  complex*16        Eij, phase_ij, tij
  external          tij_sk, tij_cc
  logical           flag_phase, flag_set_overlap, flag_load_nntable

  if(flag_phase) then
    Rij = NN_TABLE%Rij(1:3,nn)
  else
    Rij = NN_TABLE%R  (1:3,nn)
  endif

  phase_ij = FIJ_(kpoint, Rij)

  if(flag_load_nntable) then
    tij = NN_TABLE%tij_file(nn)
    Eij = tij * phase_ij
    ! NOTE: if flag_set_overlap = .true. and flag_load_nntable = .true. this will not work properly
    !       since NN_TABLE%sij_file(nn) is not properly defined.
    !       Need to be updated in the future release. 29.Oct.2020: H.-J. Kim
    return
  endif

  if(PPRAM%flag_slater_koster) then
    tij = tij_sk(NN_TABLE,nn,PPRAM,tol, flag_set_overlap)
    Eij = tij * phase_ij
    return
  else
    tij = tij_cc(NN_TABLE,nn,PPRAM,tol)
    Eij = tij * phase_ij
    return
  endif

return
endsubroutine

