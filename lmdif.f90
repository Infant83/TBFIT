#include "alias.inc"
subroutine leasqr_lm (get_eig, NN_TABLE, EDFT, PWGHT, PINPT, PPRAM, PKPTS, PGEOM, fnorm)
  use parameters
  use mpi_setup
  use print_io
  implicit none
  type (incar)                             :: PINPT 
  type (params)                            :: PPRAM 
  type (hopping), dimension(PINPT%nsystem) :: NN_TABLE
  type (weight) , dimension(PINPT%nsystem) :: PWGHT
  type (kpoints), dimension(PINPT%nsystem) :: PKPTS
  type (poscar) , dimension(PINPT%nsystem) :: PGEOM
  type (energy) , dimension(PINPT%nsystem) :: EDFT
  integer*4                                   info, maxfev
  integer*4                                   nparam_free, mpierr
  integer*4                                   ldjac, imode ! ldjac > nparam_free
  integer*4                                   i
  real*8                                      epsfcn,factor, tol, xtol, ftol, gtol
  external                                    get_eig
  logical                                     flag_write_info
  real*8                                      fnorm  ! fnorm of last step
  
  if(PPRAM%slater_koster_type .gt. 10) then
    nparam_free = PPRAM%nparam_nrl_free ! total number of free parameters
                                        ! PPRAM%nparam total number of sk parameter set
                                        ! PPRAM%nparam_free total number of sk free parameter set
  else
    nparam_free = PPRAM%nparam_free ! total number of free parameters
  endif

  if( sum(PKPTS(:)%nkpoint) .lt. nparam_free ) then
    imode = 1
    ldjac = 0
    do i = 1, PINPT%nsystem
      ldjac = ldjac + PKPTS(i)%nkpoint * PGEOM(i)%nband * PINPT%nspin
    enddo
  elseif( sum(PKPTS(:)%nkpoint) .gt. nparam_free ) then
    imode = 2
    ldjac = sum(PKPTS(:)%nkpoint)
  endif

  if( PINPT%ls_type == 'LMDIF' ) then
   write(message,*)' Start: fitting procedures with ',PINPT%ls_type,' method.'  ; write_msg
    if ( nparam_free <= 0 ) then
      write(message,*)'    !!!! WARN !!!! nparam_free (total number of free parameters) <= 0'  ; write_msg
      write(message,*)'                   Check parameter settings. Exit program.'  ; write_msg
      kill_job
    elseif( ldjac < nparam_free ) then 
      write(message,*)'    !!!! WARN !!!! number of kpoints * number of eigenvalues <= number of free parameters'  ; write_msg
      write(message,*)'                   => Increase kpoints, and try again. Exit program.'  ; write_msg
      kill_job
    endif
    factor = 100.0D+00
    maxfev = PINPT%miter * ( nparam_free + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00;epsfcn = 0.000D+00
    flag_write_info = .true.
    call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm,maxfev, epsfcn, factor, info, flag_write_info)
   if_main  call infostamp(info,PINPT%ls_type)
   write(message,*)" End: fitting procedures"  ; write_msg

  elseif( PINPT%ls_type == 'GA' ) then
    if ( nparam_free <= 0 ) then
      write(message,*)'    !!!! WARN !!!! nparam_free (number of free parameters) <= 0'  ; write_msg
      write(message,*)'                   Check parameter settings. Exit program.'  ; write_msg
      kill_job
    elseif( ldjac < nparam_free ) then
      write(message,*)'    !!!! WARN !!!! number of kpoints * number of eigenvalues <= number of free parameters'  ; write_msg
      write(message,*)'                   => Increase kpoints, and try again. Exit program.'  ; write_msg
      kill_job
    endif
    factor = 100.0D+00
    maxfev = PINPT%miter * ( nparam_free + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00;epsfcn = 0.0D+00
    flag_write_info = .false. 
    call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm, maxfev, epsfcn, factor, info, flag_write_info)
  endif

  return
endsubroutine
subroutine lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
                 ftol, xtol, gtol, fnorm, maxfev, epsfcn, factor, info, flag_write_info)

!*****************************************************************************80
!
!! LMDIF minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Modified: by Hyun-Jung Kim (KIAS, Infant@kias.re.kr)
!    14 December 2017
!    Modified for the TBFIT purpose.
!    The original source code can be found in: https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html
!  Reference:
!
  use parameters
  use cost_function
  use mpi_setup
  use kill
  use reorder_band
  use print_io
  use projected_band
  implicit none
  type(incar)                              :: PINPT
  type(params)                             :: PPRAM
  type(weight)  , dimension(PINPT%nsystem) :: PWGHT
  type(poscar)  , dimension(PINPT%nsystem) :: PGEOM
  type(hopping) , dimension(PINPT%nsystem) :: NN_TABLE
  type(energy)  , dimension(PINPT%nsystem) :: ETBA_FIT, EDFT
  type(kpoints) , dimension(PINPT%nsystem) :: PKPTS
  integer*4     nparam_free ! number of free parameters (PPRAM%nparam_nrl_free if slater_koster_type >= 10, PPRAM%nparam_free if < 10)
  integer*4     k, nsub
  integer*4     i,j,l,iter,info,ipvt(nparam_free),maxfev,nfev
  integer*4     irange(PPRAM%param_nsub_max), irange_(PPRAM%param_nsub_max)
  integer*4     ldjac ! leading demension of jacobean fjac
  integer*4     imode 
  real*8        actred,delta,diag(nparam_free),dirder,enorm,epsfcn,epsmch,factor
  real*8        fjac(ldjac, nparam_free)
  real*8        fvec(ldjac)
  real*8        fvec_plain(ldjac)
  real*8        wa4(ldjac)
  real*8        fnorm,fnorm1,fnorm_, ftol,gnorm,gtol,par
  real*8        pnorm,prered,qtf(nparam_free),ratio,sum2,temp,temp1,temp2,xnorm,xtol
  real*8        wa1(nparam_free),wa2(nparam_free),wa3(nparam_free)
  real*8        wa2_temp(nparam_free)
  character*132 pfileoutnm_temp
  logical       flag_write_info
  external      get_eig
  integer*4     i_dummy
  character*132 gnu_command
  logical       flag_wait_plot
  logical       flag_order, flag_order_weight
  integer*4     mpierr


  flag_order          = PINPT%flag_get_band_order .and. (.not. PINPT%flag_get_band_order_print_only)
  flag_order_weight   = .false. ! experimental feature

  do i = 1, PINPT%nsystem
    allocate(ETBA_FIT(i)%E(PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
    allocate(ETBA_FIT(i)%V(PGEOM(i)%neig*PINPT%ispin,PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
    allocate(ETBA_FIT(i)%ORB(PINPT%lmmax,PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
    allocate(ETBA_FIT(i)%SV(PGEOM(i)%neig*PINPT%ispin,PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
    if(flag_order) then
      allocate(ETBA_FIT(i)%E_ORD(PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
      allocate(ETBA_FIT(i)%V_ORD(PGEOM(i)%neig*PINPT%ispin,PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
      allocate(ETBA_FIT(i)%SV_ORD(PGEOM(i)%neig*PINPT%ispin,PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
      allocate(ETBA_FIT(i)%IDX(PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint))
    endif
    if(PINPT%flag_plot_fit .or. PINPT%flag_print_energy_diff) then
      flag_wait_plot = .false.
      write(gnu_command, '(A,A)')'gnuplot ', trim(PINPT%filenm_gnuplot)
    endif
  enddo

  fnorm_ = 0d0
  i_dummy = 0
  epsmch = epsilon ( epsmch )
  info = 0 ; nfev = 0
  
  if (ftol < 0.0D+00 .or. xtol < 0.0D+00 .or. gtol < 0.0D+00 .or. maxfev <= 0) go to 300

! Evaluate degeneracy information for DFT band : in the beginning
  do i = 1, PINPT%nsystem
    call get_degeneracy(EDFT(i), PGEOM(i)%neig*PINPT%ispin, PKPTS(i)%nkpoint, PINPT)
  enddo

  call get_dE(fvec, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA_FIT, PWGHT, PGEOM, PKPTS) 
  nfev = 1
  fnorm = enorm ( ldjac , fvec )

!  Initialize Levenberg-Marquardt parameter and iteration counter.
  iter = 1 ; par = 0.0D+00

  if(flag_write_info) then
    write(message,'(A)')' '  ; write_msg
    write(message,'(A,I8, 2(A,F16.6))')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ',fnorm, &
                                                      ', (EDFT-ETBA) = ', enorm ( ldjac , fvec_plain )   ; write_msg
  endif

30 continue   !  Beginning of the outer loop.
!  Calculate the jacobian matrix.
  if(flag_order) then
    do i = 1, PINPT%nsystem
      call get_ordered_band(ETBA_FIT(i), PKPTS(i), PGEOM(i), PWGHT(i), PINPT, flag_order_weight, PPRAM%flag_use_overlap)
    enddo
  endif

  call fdjac2 (get_eig,NN_TABLE,ldjac,imode,PINPT,PPRAM,PGEOM,fvec,ETBA_FIT,EDFT,nparam_free,PWGHT,fjac,epsfcn,PKPTS)

  nfev = nfev + nparam_free

!  Compute the QR factorization of the jacobian.
  call qrfac ( ldjac, nparam_free, fjac, ipvt, wa1, wa2 )
!  On the first iteration, scale according to the norms of the columns of the initial jacobian.
     if ( iter == 1 ) then
       do j = 1, PPRAM%nparam_free
         if(PPRAM%slater_koster_type .gt. 10) then
           nsub = PPRAM%param_nsub(PPRAM%iparam_free(j))
           i    = PPRAM%iparam_free_nrl(j)
           diag(i:i+nsub-1) = wa2(i:i+nsub-1)
         else
           diag(j) = wa2(j)
         endif
       enddo
       do j = 1, nparam_free
         if ( wa2(j) == 0.0D+00 ) diag(j) = 1.0D+00
       enddo

!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
       do j = 1, PPRAM%nparam_free
         if(PPRAM%slater_koster_type .gt. 10) then
           nsub = PPRAM%param_nsub(PPRAM%iparam_free(j))
           i    = PPRAM%iparam_free_nrl(j)
           wa3(i:i+nsub-1) = diag(i:i+nsub-1)*PPRAM%param_nrl(1:nsub,PPRAM%iparam_free(j))
         else
           wa3(j) = diag(j)*PPRAM%param(PPRAM%iparam_free(j))
         endif
       enddo

       xnorm = enorm ( nparam_free, wa3 )
       delta = factor * xnorm
       if ( delta == 0.0D+00 ) delta = factor
     endif  ! if iter = 1

!  Form Q' * FVEC and store the first N components in QTF.
     wa4(:) = fvec(:)
     do j = 1, nparam_free
       if ( fjac(j,j) /= 0.0D+00 ) then
         sum2 = dot_product ( wa4(j:ldjac), fjac(j:ldjac,j) )
         temp = -sum2 / fjac(j,j)
         wa4(j:ldjac) = wa4(j:ldjac) + fjac(j:ldjac,j) * temp
       endif
       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)
     enddo

!  Compute the norm of the scaled gradient.
     gnorm = 0.0D+00
     if ( fnorm /= 0.0D+00 ) then
       do j = 1, nparam_free
         l = ipvt(j)
         if (wa2(l) /= 0.0D+00) then
           sum2  = sum( fjac(:,j)*(qtf(:)/fnorm) )
           gnorm = max( gnorm, abs(sum2/wa2(l)) )
         endif
       enddo
     endif

!  Test for convergence of the gradient norm.
     if (gnorm <= gtol) then
       info = 4
       go to 300
     endif

!  Rescale if necessary.
     diag(1:nparam_free) = max ( diag(1:nparam_free), wa2(1:nparam_free) )

200  continue !  Beginning of the inner loop.

!  Determine the Levenberg-Marquardt parameter.
        call lmpar ( nparam_free, fjac, ldjac, ipvt, diag, qtf, delta, par, wa1, wa2 )

!  Store the direction P and X + P. 
!  Calculate the norm of P.
        wa1(1:nparam_free) = -wa1(1:nparam_free)

        do j = 1, PPRAM%nparam_free
          if(PPRAM%slater_koster_type .gt. 10) then
            nsub = PPRAM%param_nsub(PPRAM%iparam_free(j)) 
            i    = PPRAM%iparam_free_nrl(j)
            wa2(i:i+nsub-1) = PPRAM%param_nrl(1:nsub,PPRAM%iparam_free(j)) + wa1(i:i+nsub-1)
          else
            wa2(j) = PPRAM%param(PPRAM%iparam_free(j)) + wa1(j)
          endif
        enddo
        wa3(1:nparam_free) = diag(1:nparam_free) * wa1(1:nparam_free)
        pnorm = enorm (nparam_free, wa3)
!  On the first iteration, adjust the initial step bound.
        if ( iter == 1 ) delta = min (delta, pnorm)

!  Evaluate the function at X + P and calculate its norm.
        if(PPRAM%slater_koster_type .gt. 10) then
          do j = 1, PPRAM%nparam_free  
            nsub = PPRAM%param_nsub(PPRAM%iparam_free(j)) ; i = PPRAM%iparam_free_nrl(j)
            wa2_temp(i:i+nsub-1) = PPRAM%param_nrl(1:nsub,PPRAM%iparam_free(j))  ! store temp
            PPRAM%param_nrl(1:nsub,PPRAM%iparam_free(j)) = wa2(i:i+nsub-1)       ! update param
          enddo
        else
          wa2_temp(1:nparam_free) = PPRAM%param(PPRAM%iparam_free) ! store temp
          PPRAM%param(PPRAM%iparam_free) = wa2(1:nparam_free)      ! update param
        endif
        
        call get_dE(wa4, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA_FIT, PWGHT, PGEOM, PKPTS)

        nfev = nfev + 1
        fnorm1 = enorm ( ldjac, wa4 )

        if(PPRAM%slater_koster_type .gt. 10) then
          do j = 1, PPRAM%nparam_free
            nsub = PPRAM%param_nsub(PPRAM%iparam_free(j)) ; i = PPRAM%iparam_free_nrl(j)
            PPRAM%param_nrl(1:nsub,PPRAM%iparam_free(j)) = wa2_temp(i:i+nsub-1) ! restore param
          enddo
        else
          PPRAM%param(PPRAM%iparam_free) = wa2_temp(1:nparam_free) ! restore param
        endif

!  Compute the scaled actual reduction.
        if ( 0.1D+00 * fnorm1 < fnorm ) then
          actred = 1.0D+00 - (fnorm1/fnorm)**2
        else
          actred = -1.0D+00
        endif
!  Compute the scaled predicted reduction and the scaled directional derivative.
        do j = 1, nparam_free
          wa3(j) = 0.0D+00
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        enddo

        temp1 = enorm ( nparam_free, wa3 ) / fnorm
        temp2 = ( sqrt ( par ) * pnorm ) / fnorm
        prered = temp1 ** 2 + temp2 ** 2 / 0.5D+00
        dirder = - ( temp1 ** 2 + temp2 ** 2 )

!  Compute the ratio of the actual to the predicted reduction.
        ratio = 0.0D+00
        if ( prered /= 0.0D+00 ) ratio = actred / prered


!  Update the step bound.
        if ( ratio <= 0.25D+00 ) then
           if (actred >= 0.0D+00) temp = 0.5D+00
           if (actred < 0.0D+00) temp = 0.5D+00*dirder/(dirder+0.5D+00*actred)
           if (0.1D+00*fnorm1 >= fnorm .or. temp < 0.1D+00) temp = 0.1D+00

           delta = temp * min ( delta, pnorm / 0.1D+00  )
           par = par / temp
        else
           if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
             delta = 2.0D+00 * pnorm
             par = 0.5D+00 * par
           endif
        endif

!  Test for successful iteration.
!  Successful iteration. update X, FVEC, and their norms.
        if ( 0.0001D+00 <= ratio ) then
          do j = 1, PPRAM%nparam_free
            if(PPRAM%slater_koster_type .gt. 10) then
              nsub = PPRAM%param_nsub(PPRAM%iparam_free(j)) ; i = PPRAM%iparam_free_nrl(j)
              PPRAM%param_nrl(1:nsub,PPRAM%iparam_free(j)) = wa2(i:i+nsub-1) ! save
              wa2(i:i+nsub-1) = diag(i:i+nsub-1)*PPRAM%param_nrl(1:nsub,PPRAM%iparam_free(j)) ! update
 
            else
              PPRAM%param(PPRAM%iparam_free(j)) = wa2(j)
              wa2(j) = diag(j) * PPRAM%param(PPRAM%iparam_free(j))
            endif
          enddo
          fvec(1:ldjac) = wa4(1:ldjac)
          xnorm = enorm ( nparam_free, wa2 )
          fnorm = fnorm1
          iter = iter + 1
          if(flag_write_info) then
            write(message,'(A)')' '  ; write_msg
            write(message,'(A,I8, 2(A,F16.6))')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ',fnorm, &
                                                              ', (EDFT-ETBA) = ', enorm ( ldjac , fvec_plain )   ; write_msg
            if_main write(pfileoutnm_temp,'(A,A)')trim(PPRAM%pfileoutnm),'_temp'
            if_main call print_param(PINPT,PPRAM,PWGHT(1),pfileoutnm_temp,.TRUE.) ! only main system will be printed..
            fnorm_ = fnorm ! fnorm of previous step

            if(PINPT%flag_plot_fit .or. PINPT%flag_print_energy_diff) then

              do i = 1, PINPT%nsystem
                call get_eig(NN_TABLE(i), PKPTS(i)%kpoint, PKPTS(i)%nkpoint, PINPT, PPRAM, ETBA_FIT(i)%E, ETBA_FIT(i)%V, ETBA_FIT(i)%SV, &
                             PGEOM(i)%neig, PGEOM(i)%init_erange, PGEOM(i)%nband, PINPT%flag_get_orbital, .false., .false., PINPT%flag_phase)
                if(PINPT%flag_get_band_order) then 
                  call get_ordered_band(ETBA_FIT(i), PKPTS(i), PGEOM(i), PWGHT(i), PINPT, flag_order_weight, PPRAM%flag_use_overlap)
                  
                  if_main call print_energy(PKPTS(i), ETBA_FIT(i)%E_ORD, ETBA_FIT(i)%E_ORD, ETBA_FIT(i)%V_ORD, ETBA_FIT(i)%SV_ORD, PGEOM(i)%neig, &
                                            PINPT, PWGHT(i), PPRAM%flag_use_overlap, PINPT%flag_plot_fit,'_ordered')
                endif
                if_main call print_energy(PKPTS(i), ETBA_FIT(i)%E, ETBA_FIT(i)%E, ETBA_FIT(i)%V, ETBA_FIT(i)%SV, PGEOM(i)%neig, &
                                          PINPT, PWGHT(i), PPRAM%flag_use_overlap, PINPT%flag_plot_fit,'')
              enddo

            endif

          endif
        endif

!  Tests for convergence.
        if (abs(actred) <= ftol .and. prered <= ftol .and. 0.5D+00*ratio <= 1.0D+00) info=1
        if (delta <= xtol*xnorm ) info=2
        if (abs(actred) <= ftol .and. prered <= ftol .and. 0.5D+00*ratio <= 1.0D+00 .and. info == 2) info=3
        if ( info /= 0 ) then
          go to 300
        endif

!  Tests for termination and stringent tolerances.
        if ( maxfev <= nfev ) info=5
        if (abs(actred) <= epsmch .and. prered <= epsmch .and. 0.5D+00*ratio <= 1.0D+00) info=6
        if (delta <= epsmch*xnorm) info=7
        if (gnorm <= epsmch) info=4
        if (info /= 0) then 
          go to 300
        endif

!  End of the inner loop.  Repeat if iteration unsuccessful.
        if ( ratio < 0.0001D+00 ) go to 200

!  check kill
        call check_kill_tbfit(PINPT,PPRAM, PWGHT)

!  End of the outer loop.
     go to 30
300 continue

!  Termination, either normal or user imposed.
    if(flag_write_info) then
      if_main_then
        if(info .eq. 1) then
          write(message,*)' '  ; write_msg
          write(message,101)"  Termination INFO=",info,' , condition: |actred|,prered <= ftol, ratio <= 2 :',' actred=',actred, &
                                                                                                       ' prered=',prered, &
                                                                                                       ' ratio=',ratio  ; write_msg
        elseif(info .eq. 2) then
          write(message,*)' '  ; write_msg
          write(message,102)"  Termination INFO=",info,' , condition: delta <= xtol*xnorm :',' delta=',delta, &
                                                                                       ' xtol*xnorm=',xtol*xnorm  ; write_msg
        elseif(info .eq. 3) then
          write(message,*)' '  ; write_msg
          write(message,103)"  Termination INFO=",info,' , condition: |actred|,prered <= ftol, delta<=xtol*xnorm :',' actred=',actred, &
                                                                                                              ' prered=',prered, &
                                                                                                              ' delta=',delta, &
                                                                                                              ' xtol*xnorm=',xtol*xnorm  ; write_msg
        elseif(info .eq. 5) then
          write(message,*)' '  ; write_msg
          write(message,104)"  Termination INFO=",info,' , condition: maxfev(=miter*(nparam_free+1)) <= nfev :',' miter=',PINPT%miter, & 
                                                                                            ' nparam_free=',nparam_free, &
                                                                                            ' nfev=',nfev ; write_msg
        elseif(info .eq. 6) then
          write(message,*)' '  ; write_msg
          write(message,105)"  Termination INFO=",info,' , condition: |actred|,prered <= epsmch, ratio <= 2 :',' actred=',actred, &  
                                                                                                         ' prered=',prered, &
                                                                                                         ' epsmch=',epsmch, &
                                                                                                         ' ratio=',ratio ; write_msg
        elseif(info .eq. 7) then
          write(message,*)' '  ; write_msg
          write(message,106)"  Termination INFO=",info,' , condition: delta <= epsmch*xnorm :',' delta=',delta, &  
                                                                                         ' epsmch*xnorm=',epsmch*xnorm; write_msg
        elseif(info .eq. 4) then
          write(message,*)' '  ; write_msg
          write(message,107)"  Termination INFO=",info,' , condition: gnorm <= epsmch :',' gnorm=',gnorm, &  
                                                                                   ' epsmch=',epsmch; write_msg
        else
          write(message,*)' '  ; write_msg
          write(message,108)"  Termination INFO=",info  ; write_msg
        endif
      if_main_end 
    endif

101 format(A,I2,A,3(A,E14.7))
102 format(A,I2,A,2(A,E14.7))
103 format(A,I2,A,4(A,E14.7))
104 format(A,I2,A,3(A,I6   ))
105 format(A,I2,A,4(A,E14.7))
106 format(A,I2,A,2(A,E14.7))
107 format(A,I2,A,2(A,E14.7))
108 format(A,I2)

    if( PPRAM%slater_koster_type .gt. 10) then
      call update_param_nrl( PPRAM )
    else
      call update_param( PPRAM )
    endif

  return
endsubroutine
subroutine fdjac2 (get_eig,NN_TABLE,ldjac,imode,PINPT,PPRAM,PGEOM,fvec,ETBA_FIT,EDFT,nparam_free,PWGHT,fjac,epsfcn,PKPTS)
!*****************************************************************************80
!
!! FDJAC2 estimates an M by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This subroutine computes a forward-difference approximation
!    to the M by N jacobian matrix associated with a specified
!    problem of M functions in N variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Modified: by Hyun-Jung Kim (KIAS, Infant@kias.re.kr)
!    14 December 2017
!    Modified for the TBFIT purpose.
!    The original source code can be found in: https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
  use parameters
  use cost_function
  use mpi_setup
  use reorder_band
  use projected_band
  implicit none
  type (incar  )                           :: PINPT
  type (params )                           :: PPRAM
  type (hopping), dimension(PINPT%nsystem) :: NN_TABLE
  type (weight) , dimension(PINPT%nsystem) :: PWGHT
  type (energy) , dimension(PINPT%nsystem) :: ETBA_FIT, EDFT
  type (kpoints), dimension(PINPT%nsystem) :: PKPTS
  type (poscar) , dimension(PINPT%nsystem) :: PGEOM
  integer*4  i,j,ii, nparam_free
  integer*4  ldjac, imode
  real*8     eps,epsfcn,epsmch,h,temp,fjac(ldjac,nparam_free)
  real*8     wa(ldjac), fvec(ldjac), fvec_plain(ldjac)
  logical    flag_order_weight, flag_order
  external   get_eig
  character*20, external  :: int2str

  flag_order          = PINPT%flag_get_band_order .and. (.not. PINPT%flag_get_band_order_print_only)
  flag_order_weight   = .false.

  if(PPRAM%slater_koster_type .gt. 10) ii = 0
  epsmch = epsilon(epsmch)
  eps = sqrt( max(epsfcn,epsmch) )
  if(PPRAM%slater_koster_type .gt. 10) then
    do j = 1, PPRAM%nparam_free
        do i = 1, PPRAM%param_nsub(PPRAM%iparam_free(j))
          ! update param with h added, original param will be saved to temp, and calculate wa
          temp = PPRAM%param_nrl(i,PPRAM%iparam_free(j))
          h = eps*abs(temp)
          if (h == 0.0D+00 ) h=eps
          PPRAM%param_nrl(i,PPRAM%iparam_free(j)) = temp+h

          call get_dE(wa, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA_FIT, PWGHT, PGEOM, PKPTS)

          ! restore param from temp, and calculate derivation fjac from wa and fvec
          PPRAM%param_nrl(i,PPRAM%iparam_free(j)) = temp
          ii = ii + 1
          fjac(:,ii) = ( wa(:) - fvec(:) ) / h
        enddo
     enddo
  else

    do j = 1, PPRAM%nparam_free
      temp = PPRAM%param(PPRAM%iparam_free(j))
      h = eps*abs(temp)
      if (h == 0.0D+00 ) h=eps
      PPRAM%param(PPRAM%iparam_free(j)) = temp+h
      call get_dE(wa, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA_FIT, PWGHT, PGEOM, PKPTS)

      PPRAM%param(PPRAM%iparam_free(j)) = temp
      fjac(:,j) = ( wa(:) - fvec(:) ) / h
    enddo
  endif

  return
endsubroutine
