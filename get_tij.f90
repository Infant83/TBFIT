function tij_sk(NN_TABLE,ii,PINPT,tol,flag_init)
  use parameters, only : pi, rt2, rt3, hopping, incar
  implicit none
  type (hopping)  :: NN_TABLE
  type (incar  )  :: PINPT
  integer*4 i,ii, iscale_mode
  real*8   l, m, n, ll, mm, nn, lm, ln, mn, lmp, lmm
  real*8   rij(3),dij, d0
  real*8   tij_sk,tol
  real*8   e, s,   p,   d
  real*8      s_s, p_s, d_s
  character*8  ci_orb,cj_orb
  character*20 site_index
  real*8, external:: f_s
  real*8, external:: f_s2
  real*8, external:: e_onsite_xx
  logical  flag_init

  ! NN_TABLE%sk_set_index(0:6,nn)
  !    0        1     2    3     4      5     6
  ! e_onsite  sigma   pi delta sigma_s pi_s delta_s    ! _s indicates scaling factor

  ! NOTE: if param_class = "xx", param_name for onsite_energy should look like below,
  !       e_'orb_name'_'site_index'
  !       else if param_class = 'ss', 'pp', 'dd', 'sp' etc., i.e., normal slater-koster type,
  !       then the param_name for onsite_energy should look like below,
  !       e_'orb_name'_'species_name'
  !       For more details, see get_param_name routine to see how the 'param_name' is constructed 
  !       in the case of param_class = 'ss', 'pp', 'dd', 'sp' etc, or,
  !       see e_onsite_xx routine to see how the 'param_name' is constructed
  !       in the case of param_class = 'xx' which is based on the 'site_index' subtag.

  if(PINPT%flag_load_nntable .and. .not. flag_init) then
    tij_sk = NN_TABLE%tij_file(ii)
    return
  endif

  iscale_mode = 1 ! default = 1, see 'f_s' function for the detail.

  rij(1:3)   = NN_TABLE%Rij(1:3,ii)
  dij        = NN_TABLE%Dij(    ii)
  ci_orb     = NN_TABLE%ci_orb( ii)
  cj_orb     = NN_TABLE%cj_orb( ii)
  site_index = NN_TABLE%site_cindex( NN_TABLE%i_atom(ii) )

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

  e   = 0d0
  s   = 0d0
  p   = 0d0
  d   = 0d0
  s_s = 1d0 !default if not provided
  p_s = 1d0 !default if not provided
  d_s = 1d0 !default if not provided
  d0  = NN_TABLE%Dij0(ii)

  if( NN_TABLE%n_class(ii) .eq. 0 ) then

    if( NN_TABLE%p_class(ii) .eq. 'xx' ) then  
      ! set user defined onsite energy modifications in 'e_onsite_xx' function
      e = e_onsite_xx(ci_orb,cj_orb, PINPT, site_index)
    else
      ! set onsite energy species by species and orbital by orbital
      if(NN_TABLE%sk_index_set(0,ii) .gt. 0) then
        call get_param(PINPT, NN_TABLE%sk_index_set(0,ii), e)
        
      endif
    endif ! set_onsite

  elseif(NN_TABLE%n_class(ii) .gt. 0) then

    if(NN_TABLE%sk_index_set(1,ii) .ne. 0) then
      call get_param(PINPT, NN_TABLE%sk_index_set(1,ii), s)
    endif
   
    if(NN_TABLE%sk_index_set(2,ii) .ne. 0) then
      call get_param(PINPT, NN_TABLE%sk_index_set(2,ii), p)
    endif
   
    if(NN_TABLE%sk_index_set(3,ii) .ne. 0) then 
      call get_param(PINPT, NN_TABLE%sk_index_set(3,ii), d)
    endif
   
    if(NN_TABLE%sk_index_set(4,ii) .gt. 0) then
      call get_param(PINPT, NN_TABLE%sk_index_set(4,ii), s_s)
    endif
   
    if(NN_TABLE%sk_index_set(5,ii) .ne. 0) then 
      call get_param(PINPT, NN_TABLE%sk_index_set(5,ii), p_s)
    endif
   
    if(NN_TABLE%sk_index_set(6,ii) .ne. 0) then
      call get_param(PINPT, NN_TABLE%sk_index_set(6,ii), d_s)
    endif

  endif !check n_class

  s = s * f_s ( s_s, d0, dij, iscale_mode)
  p = p * f_s ( p_s, d0, dij, iscale_mode)
  d = d * f_s ( d_s, d0, dij, iscale_mode)

  ! SK-energy integral if nn_class > 0
  if( NN_TABLE%n_class(ii) .ne. 0) then
sk: select case ( NN_TABLE%p_class(ii) )

      case ('ss')
        if    (ci_orb(1:1) .eq. 's' .and. cj_orb(1:1) .eq. 's') then
          tij_sk = s
        else
          write(6,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.'
          write(6,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb)
          write(6,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...'
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
          write(6,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.'
          write(6,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb)
          write(6,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...'
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
          write(6,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.'
          write(6,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb)
          write(6,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...'
          stop
        endif

      case ('sp')
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
          write(6,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.'
          write(6,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb)
          write(6,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...'
          stop
        endif

      case ('sd')
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
          write(6,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.'
          write(6,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb)
          write(6,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...'
          stop
        endif

      case ('pd')
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
          write(6,'(A)')' !WARNING! SK energy integral is not properly defined or orbital name is improper.'
          write(6,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb)
          write(6,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...'
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
          write(6,'(A)')' !WARNING! SK_C energy integral is not properly defined or orbital name is improper.'
          write(6,'(A,A)')' !WARNING! CI_ORB = ',trim(ci_orb)
          write(6,'(A,A)')' !WARNING! CJ_ORB = ',trim(cj_orb),' Exit...'
          stop
        endif

    end select sk

  ! onsite energy if nn_class == 0
  elseif( NN_TABLE%n_class(ii) .eq. 0 ) then 
    if(PINPT%flag_efield .and. (NN_TABLE%i_matrix(ii) .eq. NN_TABLE%j_matrix(ii)) ) then
      tij_sk = e - dot_product(PINPT%efield(1:3), NN_TABLE%i_coord(1:3,ii) - PINPT%efield_origin_cart(1:3))
    else
      tij_sk = e 
    endif
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
!  4 = see Europhys.Lett.9.701 (1989): GSP parameterization for carbon system
!      SFACTOR(1) = m ; SFACTOR(2) = d_c , critical distance ; SFACTOR(3) = m_c
!    f_s=(d0/d)^(SFACTOR(1))*exp(SFACTOR(1)*( -(d/SFACTOR(2))^(SFACTOR(3)) + -(d0/SFACTOR(2))^(SFACTOR(3)) ))

   if(mode .eq. 1) then
     f_s = Exp( (d0 - d)/(dda_s*d0) )
!    f_s = Exp( (d0 - d)/(dda_s) )
   elseif(mode .eq. 2) then
     f_s = Exp( -(d/d0)**dda_s )
   elseif(mode .eq. 3) then
     f_s = (d0/d)**(dda_s)
!  elseif(mode .eq. 4)
!    f_s = (d0/d)**( dda_s(1) ) * Exp( dda_s(1) * ( -(d/dda_s(2))**dda_s(3) - (d0/dda_s(2))**dda_s(3) ) )
   endif
return
endfunction
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

