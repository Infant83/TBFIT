#include "alias.inc"
subroutine leasqr_lm (get_eig, NN_TABLE, EDFT, neig, PWGHT, PINPT, PKPTS, fnorm)
  use parameters
  use mpi_setup
  use print_io
  implicit none
  type (incar)           :: PINPT 
  type (hopping)         :: NN_TABLE
  type (weight)          :: PWGHT
  type (kpoints)         :: PKPTS
  type (energy)          :: EDFT
  integer*4                 nkpoint, neig, iband, nband, info, maxfev
  integer*4                 nparam, nparam_free, mpierr
  integer*4                 ldjac, imode
  integer*4                 i, j, k, ii
  real*8                    epsfcn,factor, tol, xtol, ftol, gtol, kpoint(3,PKPTS%nkpoint)
  external                  get_eig
  logical                   flag_write_info
  real*8                    fnorm(2) ! fnorm(1) fnorm of last step, fnorm(2) delta fnorm of last step
  
  nkpoint = PKPTS%nkpoint
   kpoint = PKPTS%kpoint
  iband   = PINPT%init_erange
  nband   = PINPT%nband

  if(PINPT%slater_koster_type .gt. 10) then
    nparam = PINPT%nparam_nrl ! total number of parameters
    nparam_free = PINPT%nparam_nrl_free ! total number of free parameters
                                        ! PINPT%nparam total number of sk parameter set
                                        ! PINPT%nparam_free total number of sk free parameter set
  else
    nparam = PINPT%nparam ! total number of parameters 
    nparam_free = PINPT%nparam_free ! total number of free parameters
  endif
  if( nkpoint * nband * PINPT%nspin .le. nparam) then
    ldjac = nkpoint * nband * PINPT%nspin
    imode = 1
  else
    ldjac = nkpoint
    imode = 2
  endif

  if( PINPT%ls_type == 'LMDIF' ) then
   write(message,*)' Start: fitting procedures with ',PINPT%ls_type,' method.'  ; write_msg
    if ( nparam_free <= 0 ) then
      write(message,*)'    !!!! WARN !!!! nparam_free (total number of free parameters) <= 0'  ; write_msg
      write(message,*)'                   Check parameter settings. Exit program.'  ; write_msg
      kill_job
    elseif(nkpoint * nband * PINPT%nspin < nparam_free ) then 
      write(message,*)'    !!!! WARN !!!! number of kpoints * number of eigenvalues <= number of free parameters'  ; write_msg
      write(message,*)'                   => Increase kpoints, and try again. Exit program.'  ; write_msg
      kill_job
    endif
    factor = 100.0D+00
    maxfev = PINPT%miter * ( nparam_free + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00;epsfcn = 0.000D+00
    flag_write_info = .true.
    call lmdif(get_eig, NN_TABLE, kpoint, nkpoint, ldjac, imode, PINPT, PKPTS, EDFT, nparam, nparam_free, neig, iband, nband, PWGHT, &
               ftol, xtol, gtol, fnorm(1),fnorm(2), maxfev, epsfcn, factor, info, flag_write_info)
   if_main  call infostamp(info,PINPT%ls_type)
   write(message,*)" End: fitting procedures"  ; write_msg

  elseif( PINPT%ls_type == 'GA' ) then
    if ( nparam_free <= 0 ) then
      write(message,*)'    !!!! WARN !!!! nparam_free (number of free parameters) <= 0'  ; write_msg
      write(message,*)'                   Check parameter settings. Exit program.'  ; write_msg
      kill_job
    elseif(nkpoint * nband * PINPT%nspin < nparam_free ) then
      write(message,*)'    !!!! WARN !!!! number of kpoints * number of eigenvalues <= number of free parameters'  ; write_msg
      write(message,*)'                   => Increase kpoints, and try again. Exit program.'  ; write_msg
      kill_job
    endif
    factor = 100.0D+00
    maxfev = PINPT%miter * ( nparam_free + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00;epsfcn = 0.0D+00
    flag_write_info = .false. 
    call lmdif(get_eig, NN_TABLE, kpoint, nkpoint, ldjac, imode, PINPT, PKPTS, EDFT, nparam, nparam_free, neig, iband, nband, PWGHT, &
               ftol, xtol, gtol, fnorm(1),fnorm(2), maxfev, epsfcn, factor, info, flag_write_info)

  endif

  return
endsubroutine
subroutine lmdif(get_eig, NN_TABLE, kpoint, nkpoint, ldjac, imode, PINPT, PKPTS, EDFT, nparam, nparam_free, neig, iband, nband, PWGHT, &
                 ftol, xtol, gtol, fnorm,fnorm2, maxfev, epsfcn, factor, info, flag_write_info)
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
  implicit none
  type(incar)   :: PINPT
  type(weight)  :: PWGHT
  type(hopping) :: NN_TABLE
  type(energy)  :: ETBA_FIT, EDFT
  type(kpoints) :: PKPTS
  integer*4     nkpoint,neig,iband,nband
  integer*4     nparam ! number of parameters (differ from PINPT%nparam if slater_koster_type > 10)
  integer*4     nparam_free ! number of free parameters (differ from PINPT%nparam_free if slater_koster_type > 10)
  integer*4     k, nsub
  integer*4     i,j,l,iter,info,ipvt(nparam_free),maxfev,nfev
  integer*4     irange(PINPT%param_nsub_max), irange_(PINPT%param_nsub_max)
  integer*4     ldjac ! leading demension of jacobean fjac
  integer*4     imode 
  real*8        kpoint(3,nkpoint)
  real*8        actred,delta,diag(nparam_free),dirder,enorm,epsfcn,epsmch,factor
  real*8        fjac(ldjac, nparam_free)
  real*8        fvec(ldjac)
  real*8        fvec_plain(ldjac)
  real*8        dvec(ldjac)
  real*8        wa4(ldjac)
  real*8        fnorm,fnorm1,fnorm2,fnorm_, ftol,gnorm,gtol,par
  real*8        pnorm,prered,qtf(nparam_free),ratio,sum2,temp,temp1,temp2,xnorm,xtol
  real*8        wa1(nparam_free),wa2(nparam_free),wa3(nparam_free)
  real*8        wa2_temp(nparam_free)
  character*132 pfileoutnm_temp
  logical       flag_get_orbital, flag_collinear, flag_noncollinear, flag_soc
  logical       flag_write_info
  external      get_eig
  integer*4     i_dummy
  character*132 gnu_command
  logical       flag_wait_plot
  logical       flag_fit_degeneracy, flag_order, flag_order_weight
  integer*4     mpierr


  flag_order          = PINPT%flag_get_band_order .and. (.not. PINPT%flag_get_band_order_print_only)
  flag_get_orbital    = ((.not. PWGHT%flag_weight_default_orb) .or. flag_order)
  flag_fit_degeneracy = PINPT%flag_fit_degeneracy
  flag_order_weight   = .false. ! experimental feature

  ! ETBA_FIT: this variable is temporal and used only in this routine 
  allocate(ETBA_FIT%E(nband*PINPT%nspin, nkpoint))
  allocate(ETBA_FIT%V(neig*PINPT%ispin,PINPT%nband*PINPT%nspin, nkpoint))
  if(flag_order) then
    allocate(ETBA_FIT%E_ORD(nband*PINPT%nspin, nkpoint))
    allocate(ETBA_FIT%V_ORD(neig*PINPT%ispin,PINPT%nband*PINPT%nspin, nkpoint))
    allocate(ETBA_FIT%IDX(nband*PINPT%nspin, nkpoint))
  endif
  if(PINPT%flag_plot_fit .or. PINPT%flag_print_energy_diff) then
    flag_wait_plot = .false.
    write(gnu_command, '(A,A)')'gnuplot ', trim(PINPT%filenm_gnuplot)
  endif


  fnorm_ = 0d0
  i_dummy = 0
  epsmch = epsilon ( epsmch )
  info = 0 ; nfev = 0
  
  if (ftol < 0.0D+00 .or. xtol < 0.0D+00 .or. gtol < 0.0D+00 .or. maxfev <= 0) go to 300
!  Evaluate the function at the starting point and calculate its norm.
  call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, ETBA_FIT%E, ETBA_FIT%V,  &
               neig, iband, nband, flag_get_orbital, .false., .false., PINPT%flag_phase) !, flag_order)
! Evaluate degeneracy information for DFT band : in the beginning
  if(flag_fit_degeneracy) call get_degeneracy(EDFT%E, EDFT%D, neig*PINPT%ispin, nkpoint, PINPT)
! Evaluate degeneracy information for DFT band : after calling get_eig
  if(flag_fit_degeneracy) call get_degeneracy(ETBA_FIT%E, ETBA_FIT%D, nband*PINPT%nspin,nkpoint, PINPT)

  if(flag_order) call get_ordered_band(ETBA_FIT, nkpoint, neig, iband, nband, PINPT, flag_order_weight, PWGHT) 
  call get_cost_function(fvec  , ETBA_FIT, EDFT, neig, iband, nband, PINPT, nkpoint, PWGHT, ldjac, imode, flag_order, fvec_plain)
  nfev = 1
  fnorm = enorm ( ldjac , fvec )

!  Initialize Levenberg-Marquardt parameter and iteration counter.
  iter = 1 ; par = 0.0D+00

  if(flag_write_info) then
    write(message,'(A)')' '  ; write_msg
   !write(message,'(A,I4,A,F16.6)')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ',fnorm  ; write_msg
    write(message,'(A,I8, 2(A,F16.6))')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ',fnorm, ', (EDFT-ETBA) = ', enorm ( ldjac , fvec_plain )   ; write_msg
  endif

30 continue   !  Beginning of the outer loop.
!  Calculate the jacobian matrix.
  if(flag_order) call get_ordered_band(ETBA_FIT, nkpoint, neig, iband, nband, PINPT, flag_order_weight, PWGHT) 
  call fdjac2 ( get_eig, NN_TABLE, ldjac, imode, kpoint, nkpoint, PINPT, fvec, ETBA_FIT, EDFT, nparam_free, &
                neig, iband, nband, PWGHT, fjac, epsfcn, flag_get_orbital,PKPTS)

  nfev = nfev + nparam_free

!  Compute the QR factorization of the jacobian.
  call qrfac ( ldjac, nparam_free, fjac, ipvt, wa1, wa2 )
!  On the first iteration, scale according to the norms of the columns of the initial jacobian.
     if ( iter == 1 ) then
       do j = 1, PINPT%nparam_free
         if(PINPT%slater_koster_type .gt. 10) then
           nsub = PINPT%param_nsub(PINPT%iparam_free(j))
           i    = PINPT%iparam_free_nrl(j)
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
       do j = 1, PINPT%nparam_free
         if(PINPT%slater_koster_type .gt. 10) then
           nsub = PINPT%param_nsub(PINPT%iparam_free(j))
           i    = PINPT%iparam_free_nrl(j)
           wa3(i:i+nsub-1) = diag(i:i+nsub-1)*PINPT%param_nrl(1:nsub,PINPT%iparam_free(j))
         else
           wa3(j) = diag(j)*PINPT%param(PINPT%iparam_free(j))
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

        do j = 1, PINPT%nparam_free
          if(PINPT%slater_koster_type .gt. 10) then
            nsub = PINPT%param_nsub(PINPT%iparam_free(j)) 
            i    = PINPT%iparam_free_nrl(j)
            wa2(i:i+nsub-1) = PINPT%param_nrl(1:nsub,PINPT%iparam_free(j)) + wa1(i:i+nsub-1)
          else
            wa2(j) = PINPT%param(PINPT%iparam_free(j)) + wa1(j)
          endif
        enddo
        wa3(1:nparam_free) = diag(1:nparam_free) * wa1(1:nparam_free)
        pnorm = enorm (nparam_free, wa3)
!  On the first iteration, adjust the initial step bound.
        if ( iter == 1 ) delta = min (delta, pnorm)

!  Evaluate the function at X + P and calculate its norm.
        if(PINPT%slater_koster_type .gt. 10) then
          do j = 1, PINPT%nparam_free  
            nsub = PINPT%param_nsub(PINPT%iparam_free(j)) ; i = PINPT%iparam_free_nrl(j)
            wa2_temp(i:i+nsub-1) = PINPT%param_nrl(1:nsub,PINPT%iparam_free(j))  ! store temp
            PINPT%param_nrl(1:nsub,PINPT%iparam_free(j)) = wa2(i:i+nsub-1)       ! update param
          enddo
          call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, ETBA_FIT%E, ETBA_FIT%V, neig, iband, nband, flag_get_orbital, .false., .false., PINPT%flag_phase)
          do j = 1, PINPT%nparam_free
            nsub = PINPT%param_nsub(PINPT%iparam_free(j)) ; i = PINPT%iparam_free_nrl(j)
            PINPT%param_nrl(1:nsub,PINPT%iparam_free(j)) = wa2_temp(i:i+nsub-1) ! restore param
          enddo
        else
          wa2_temp(1:nparam_free) = PINPT%param(PINPT%iparam_free) ! store temp
          PINPT%param(PINPT%iparam_free) = wa2(1:nparam_free)      ! update param
          call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, ETBA_FIT%E, ETBA_FIT%V, neig, iband, nband, flag_get_orbital, .false., .false., PINPT%flag_phase)
          PINPT%param(PINPT%iparam_free) = wa2_temp(1:nparam_free) ! restore param
        endif

        if(flag_fit_degeneracy) call get_degeneracy(ETBA_FIT%E, ETBA_FIT%D, nband*PINPT%nspin, nkpoint, PINPT)
        if(flag_order) call get_ordered_band(ETBA_FIT, nkpoint, neig, iband, nband, PINPT, flag_order_weight, PWGHT) 
        call get_cost_function(wa4  , ETBA_FIT, EDFT, neig, iband, nband, PINPT, nkpoint, PWGHT, ldjac, imode, flag_order, fvec_plain)
        nfev = nfev + 1
        fnorm1 = enorm ( ldjac, wa4 )

!if_main write(6,'(A, I6,I6, 2(F20.10))')"VVV ", iter, nfev, sum(wa4), fnorm1

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
          do j = 1, PINPT%nparam_free
            if(PINPT%slater_koster_type .gt. 10) then
              nsub = PINPT%param_nsub(PINPT%iparam_free(j)) ; i = PINPT%iparam_free_nrl(j)
              PINPT%param_nrl(1:nsub,PINPT%iparam_free(j)) = wa2(i:i+nsub-1) ! save
              wa2(i:i+nsub-1) = diag(i:i+nsub-1)*PINPT%param_nrl(1:nsub,PINPT%iparam_free(j)) ! update
 
            else
              PINPT%param(PINPT%iparam_free(j)) = wa2(j)
              wa2(j) = diag(j) * PINPT%param(PINPT%iparam_free(j))
            endif
          enddo
          fvec(1:ldjac) = wa4(1:ldjac)
          xnorm = enorm ( nparam_free, wa2 )
          fnorm = fnorm1
          iter = iter + 1
          if(flag_write_info) then
            fnorm2 = abs(fnorm_ - fnorm)/fnorm ! delta fnorm
            write(message,'(A)')' '  ; write_msg
           !write(message,'(A,I4,A,F16.6)')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ',fnorm  ; write_msg
            write(message,'(A,I8, 2(A,F16.6))')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ',fnorm, ', (EDFT-ETBA) = ', enorm ( ldjac , fvec_plain )   ; write_msg
            if_main write(pfileoutnm_temp,'(A,A)')trim(PINPT%pfileoutnm),'_temp'
            if_main call print_param(PINPT,pfileoutnm_temp,.TRUE.)
            fnorm_ = fnorm ! fnorm of previous step

            if(PINPT%flag_plot_fit .or. PINPT%flag_print_energy_diff) then
              call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, ETBA_FIT%E, ETBA_FIT%V, &
                                     neig, PINPT%init_erange, PINPT%nband, PINPT%flag_get_orbital, .false., .false., PINPT%flag_phase)
              if(PINPT%flag_get_band_order) then 
                call get_ordered_band(ETBA_FIT, nkpoint, neig, PINPT%init_erange, PINPT%nband, PINPT, flag_order_weight, PWGHT) 
                if_main call print_energy(PKPTS, ETBA_FIT%E_ORD, ETBA_FIT%E_ORD, ETBA_FIT%V_ORD, neig, PINPT, PWGHT, PINPT%flag_plot_fit,'_ordered')
              endif
              if_main call print_energy(PKPTS, ETBA_FIT%E, ETBA_FIT%E, ETBA_FIT%V, neig, PINPT, PWGHT, PINPT%flag_plot_fit,'')
              
             !if(PINPT%flag_plot_fit) then
             !  if_main call execute_command_line(gnu_command, flag_wait_plot)
             !endif
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
        call check_kill_tbfit(PINPT)

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

    if( PINPT%slater_koster_type .gt. 10) then
      call update_param_nrl( PINPT )
    else
      call update_param( PINPT )
    endif
   !if(flag_write_info) then
   !  if_main call print_param (PINPT, '  Fitted param(i):', .FALSE.)
   !endif

  return
endsubroutine
subroutine param_to_param_nrl(param,param_nsub, nparam_tot, nparam, param_nrl)
  use mpi_setup
  integer*4    i,j, nsub
  integer*4    nparam, nparam_tot
  integer*4    param_nsub(nparam)
  real*8       param_nrl(4,nparam), param(nparam_tot)

  do j = 1, nparam
    nsub = param_nsub(j)
    i    = sum(param_nsub(1:j-1)) + 1
    param_nrl(1:nsub,j) = param(i:i+nsub-1)
  enddo

  return
endsubroutine
subroutine param_nrl_to_param(param_nrl,param_nsub, iparam_free_nrl, nparam_tot, nparam, param)
  use mpi_setup
  integer*4    i,j, nsub
  integer*4    nparam, nparam_tot
  integer*4    param_nsub(nparam)
  integer*4    iparam_free_nrl(nparam)
  real*8       param_nrl(4,nparam), param(nparam_tot)
  
  do j = 1, nparam
    nsub = param_nsub(j)
   !i    = sum(param_nsub(1:j-1)) + 1
    i    = iparam_free_nrl(j)
    param(i:i+nsub-1) = param_nrl(1:nsub,j)
  enddo

  return
endsubroutine
subroutine fdjac2 (get_eig, NN_TABLE, ldjac, imode, kpoint, nkpoint, PINPT, fvec, ETBA_FIT, EDFT, nparam_free, &
                   neig, iband, nband, PWGHT, fjac, epsfcn, flag_get_orbital,PKPTS)

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
  implicit none
  type (incar  ) :: PINPT
  type (hopping) :: NN_TABLE
  type (weight)  :: PWGHT
  type (energy)  :: ETBA_FIT, EDFT
  type (kpoints) :: PKPTS
  integer*4  nkpoint,neig,iband,nband
  integer*4  i,j,ii, nparam_free
  integer*4  ldjac, imode
  real*8     eps,epsfcn,epsmch,h,temp,fjac(ldjac,nparam_free)
  real*8     wa(ldjac), fvec(ldjac)
  real*8                fvec_plain(ldjac)
  real*8     kpoint(3,nkpoint)
  logical    flag_get_orbital, flag_order_weight, flag_order
  external   get_eig
  character*20, external  :: int2str

  flag_order_weight   = .false.
  flag_order          = PINPT%flag_get_band_order .and. (.not. PINPT%flag_get_band_order_print_only)

  if(PINPT%slater_koster_type .gt. 10) ii = 0
  epsmch = epsilon(epsmch)
  eps = sqrt( max(epsfcn,epsmch) )
  if(PINPT%slater_koster_type .gt. 10) then
    do j = 1, PINPT%nparam_free
        do i = 1, PINPT%param_nsub(PINPT%iparam_free(j))
          ! update param with h added, original param will be saved to temp, and calculate wa
          temp = PINPT%param_nrl(i,PINPT%iparam_free(j))
          h = eps*abs(temp)
          if (h == 0.0D+00 ) h=eps
          PINPT%param_nrl(i,PINPT%iparam_free(j)) = temp+h
          call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, ETBA_FIT%E, ETBA_FIT%V, neig, &
                       iband, nband, flag_get_orbital, .false., .false., PINPT%flag_phase)
  
          if(PINPT%flag_fit_degeneracy) call get_degeneracy(ETBA_FIT%E, ETBA_FIT%D, nband*PINPT%nspin, nkpoint, PINPT)
  
          if(flag_order) call get_ordered_band(ETBA_FIT, nkpoint, neig, iband, nband, PINPT, flag_order_weight, PWGHT) 
          call get_cost_function(wa, ETBA_FIT, EDFT, neig, iband, nband, PINPT, nkpoint, PWGHT, ldjac, imode, flag_order, fvec_plain)
          ! restore param from temp, and calculate derivation fjac from wa and fvec
          PINPT%param_nrl(i,PINPT%iparam_free(j)) = temp
          ii = ii + 1
          fjac(:,ii) = ( wa(:) - fvec(:) ) / h
        enddo
     enddo
  else
!write(6,'(A,F9.4, *(I3))')"KKK ", sum(wa), 0,     ETBA_FIT%IDX(1:11,5)
    do j = 1, PINPT%nparam_free
      temp = PINPT%param(PINPT%iparam_free(j))
      h = eps*abs(temp)
      if (h == 0.0D+00 ) h=eps
      PINPT%param(PINPT%iparam_free(j)) = temp+h
      call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, ETBA_FIT%E, ETBA_FIT%V, neig, &
                   iband, nband, flag_get_orbital, .false., .false., .true.)
      if(PINPT%flag_fit_degeneracy) call get_degeneracy(ETBA_FIT%E, ETBA_FIT%D, nband*PINPT%nspin, nkpoint, PINPT)

      !if(flag_order) call get_ordered_band(ETBA_FIT, nkpoint, neig, iband, nband, PINPT, flag_order_weight, PWGHT) 
      if(flag_order) call update_order(ETBA_FIT%IDX, ETBA_FIT%E, ETBA_FIT%E_ORD, nkpoint, nband, PINPT)
      call get_cost_function(wa, ETBA_FIT, EDFT, neig, iband, nband, PINPT, nkpoint, PWGHT, ldjac, imode, flag_order, fvec_plain)

!do ii = 1, 10
!write(6,'(A,F9.4, *(I3))')"KKK ", sum(wa), j,     ETBA_FIT%IDX(1:11,5)
!call print_energy(PKPTS, ETBA_FIT%E_ORD, ETBA_FIT%E_ORD, ETBA_FIT%V_ORD, neig, PINPT, PWGHT, .true., '_'//adjustl(trim(int2str(j))))
!enddo
      PINPT%param(PINPT%iparam_free(j)) = temp
      fjac(:,j) = ( wa(:) - fvec(:) ) / h
!write(6,'(A, 72(F19.4), I3     )')"ZZZ ", fjac(:,j), j
!write(6,'(A, I3, 72(F19.12)  )')"ZZZ ", j, h
!write(6,'(A, I3, 72(F19.8)  )')"JJJ ", j, temp+h
    enddo
  endif
!stop
  return
endsubroutine
