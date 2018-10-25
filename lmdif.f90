#include "alias.inc"
subroutine leasqr_lm (get_eig, NN_TABLE, kpoint, nkpoint, E_DFT, neig, iband, nband, PWGHT, PINPT)
  use parameters
  use mpi_setup
  implicit none
  type (incar)   :: PINPT 
  type (hopping) :: NN_TABLE
  type (weight)  :: PWGHT
  integer*4 nkpoint, neig, iband, nband, info, maxfev
  real*8  epsfcn,factor, tol, xtol, ftol, gtol, kpoint(3,nkpoint)
  real*8  E_DFT(neig*PINPT%ispin,nkpoint)
  external get_eig
  if( PINPT%ls_type == 'LMDIF' ) then
   if_main write(6,*)' Start: fitting procedures with ',PINPT%ls_type,' method.'
    if ( PINPT%nparam <= 0 .or. nkpoint < PINPT%nparam .or. tol < 0.0D+00) return
    factor = 100.0D+00
    maxfev = PINPT%miter * ( PINPT%nparam + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00;epsfcn = 0.0D+00
    call lmdif(get_eig, NN_TABLE, kpoint, nkpoint, PINPT, E_DFT, neig, iband, nband, PWGHT, &
               ftol, xtol, gtol, maxfev, epsfcn, factor, info)
   if_main  call infostamp(info,PINPT%ls_type)
   if_main write(6,*)" End: fitting procedures"
  endif
  return
endsubroutine
subroutine lmdif(get_eig, NN_TABLE, kpoint, nkpoint, PINPT, E_DFT, neig, iband, nband, PWGHT, &
                 ftol, xtol, gtol, maxfev, epsfcn, factor, info)
  use parameters
  use mpi_setup
  implicit none
  type(incar)   :: PINPT
  type(weight)  :: PWGHT
  type(hopping) :: NN_TABLE
  integer*4     nkpoint,neig,iband,nband
  integer*4     i,j,l,iter,info,ipvt(PINPT%nparam),maxfev,nfev
  real*8        kpoint(3,nkpoint)
  real*8        actred,delta,diag(PINPT%nparam),dirder,enorm,epsfcn,epsmch,factor
  real*8        fjac(nkpoint,PINPT%nparam),fvec(nkpoint)
  real*8        E_DFT(neig*PINPT%ispin,nkpoint),E_TBA(nband*PINPT%nspin,nkpoint),wa4(nkpoint)
  complex*16    V(neig*PINPT%ispin,nband*PINPT%nspin,nkpoint)
  real*8        fnorm,fnorm1,ftol,gnorm,gtol,par
  real*8        pnorm,prered,qtf(PINPT%nparam),ratio,sum2,temp,temp1,temp2,xnorm,xtol
  real*8        wa1(PINPT%nparam),wa2(PINPT%nparam),wa3(PINPT%nparam),param(PINPT%nparam)
  real*8        wa2_temp(PINPT%nparam)
  character*132 pfileoutnm_temp
  logical       flag_get_orbital, flag_collinear, flag_noncollinear, flag_soc
  external      get_eig
! flag_collinear = PINPT%flag_collinear
! flag_noncollinear = PINPT%flag_noncollinear
! flag_soc = PINPT%flag_soc 
  flag_get_orbital = (.not. PWGHT%flag_weight_default_orb)

  epsmch = epsilon ( epsmch )
  info = 0 ; nfev = 0
  if (ftol < 0.0D+00 .or. xtol < 0.0D+00 .or. gtol < 0.0D+00 .or. maxfev <= 0) go to 300

!  Evaluate the function at the starting point and calculate its norm.
  call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, E_TBA, V,  neig, iband, nband, flag_get_orbital, .false., .false., .true.)
  if(.not. flag_get_orbital) then
    call get_fvec(fvec, E_TBA, E_DFT, 0, neig, iband, nband, PINPT, nkpoint, PWGHT)
  else
    call get_fvec(fvec, E_TBA, E_DFT, V, neig, iband, nband, PINPT, nkpoint, PWGHT)
  endif
  nfev = 1
  fnorm = enorm ( nkpoint, fvec )

!  Initialize Levenberg-Marquardt parameter and iteration counter.
  iter = 1 ; par = 0.0D+00
30 continue   !  Beginning of the outer loop.
!  Calculate the jacobian matrix.
  call fdjac2 ( get_eig, NN_TABLE, kpoint, nkpoint, PINPT, fvec, E_DFT, neig, iband, nband, PWGHT, fjac, epsfcn, flag_get_orbital)
  nfev = nfev + PINPT%nparam
!  Compute the QR factorization of the jacobian.
  call qrfac ( nkpoint, PINPT%nparam, fjac, ipvt, wa1, wa2 )

!  On the first iteration, scale according to the norms of the columns of the initial jacobian.
     if ( iter == 1 ) then
       do j = 1, PINPT%nparam
         if( nint(PINPT%param_const(1,j)) .eq. 0 ) then
           diag(j) = wa2(j)
         elseif( nint(PINPT%param_const(1,j)) .ge. 1) then
           diag(j) = wa2( nint(PINPT%param_const(1,j)) )
         endif
       enddo
       do j = 1, PINPT%nparam
         if ( wa2(j) == 0.0D+00 ) diag(j) = 1.0D+00
       enddo
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
       do j = 1, PINPT%nparam
         if( nint(PINPT%param_const(1,j)) .eq. 0 ) then
           wa3(j) = diag(j)*PINPT%param(j)
         elseif( nint(PINPT%param_const(1,j)) .ge. 1) then
           wa3(j) = diag( nint(PINPT%param_const(1,j)) )*PINPT%param( nint(PINPT%param_const(1,j)) )
         endif
       enddo
       xnorm = enorm ( PINPT%nparam, wa3 )
       delta = factor * xnorm
       if ( delta == 0.0D+00 ) delta = factor
     endif

!  Form Q' * FVEC and store the first N components in QTF.
     wa4(1:nkpoint) = fvec(1:nkpoint)
     do j = 1, PINPT%nparam
       if ( fjac(j,j) /= 0.0D+00 ) then
         sum2 = dot_product ( wa4(j:nkpoint), fjac(j:nkpoint,j) )
         temp = -sum2 / fjac(j,j)
         wa4(j:nkpoint) = wa4(j:nkpoint) + fjac(j:nkpoint,j) * temp
       endif
       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)
     enddo

!  Compute the norm of the scaled gradient.
     gnorm = 0.0D+00
     if ( fnorm /= 0.0D+00 ) then
       do j = 1, PINPT%nparam
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
     diag(1:PINPT%nparam) = max ( diag(1:PINPT%nparam), wa2(1:PINPT%nparam) )

200  continue !  Beginning of the inner loop.
!  Determine the Levenberg-Marquardt parameter.
        call lmpar ( PINPT%nparam, fjac, nkpoint, ipvt, diag, qtf, delta, par, wa1, wa2 )

!  Store the direction P and X + P. 
!  Calculate the norm of P.
        wa1(1:PINPT%nparam) = -wa1(1:PINPT%nparam)
        do j = 1, PINPT%nparam
          if( nint(PINPT%param_const(1,j)) .eq. 0 ) then
            wa2(j) = PINPT%param(j) + wa1(j)
if(j .eq. 12) write(6,*)"ZZZZZ", wa2(12), nint(PINPT%param_const(1,12))
          elseif( nint(PINPT%param_const(1,j)) .ge. 1) then
            wa2(j) = PINPT%param( nint(PINPT%param_const(1,j)) ) + wa1(j)
          endif
        enddo
        wa3(1:PINPT%nparam) = diag(1:PINPT%nparam) * wa1(1:PINPT%nparam)
        pnorm = enorm (PINPT%nparam, wa3)

!  On the first iteration, adjust the initial step bound.
        if ( iter == 1 ) delta = min (delta, pnorm)

!  Evaluate the function at X + P and calculate its norm.
        wa2_temp = PINPT%param 
        PINPT%param = wa2
        call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, E_TBA, V, neig, iband, nband, flag_get_orbital, .false., .false., .true.)
        PINPT%param = wa2
        wa2 = PINPT%param
        PINPT%param = wa2_temp
        if(.not. flag_get_orbital) then
          call get_fvec(wa4, E_TBA, E_DFT, 0, neig, iband, nband, PINPT, nkpoint, PWGHT)
        else
          call get_fvec(wa4, E_TBA, E_DFT, V, neig, iband, nband, PINPT, nkpoint, PWGHT)
        endif
        nfev = nfev + 1
        fnorm1 = enorm ( nkpoint, wa4 )

!  Compute the scaled actual reduction.
        if ( 0.1D+00 * fnorm1 < fnorm ) then
          actred = 1.0D+00 - (fnorm1/fnorm)**2
        else
          actred = -1.0D+00
        endif

!  Compute the scaled predicted reduction and the scaled directional derivative.
        do j = 1, PINPT%nparam
          wa3(j) = 0.0D+00
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        enddo

        temp1 = enorm ( PINPT%nparam, wa3 ) / fnorm
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
          do j = 1, PINPT%nparam
            if( nint(PINPT%param_const(1,j)) .eq. 0 ) then
              PINPT%param(j) = wa2(j)
              wa2(j) = diag(j) * PINPT%param(j)
            elseif( nint(PINPT%param_const(1,j)) .ge. 1 ) then
              PINPT%param(j) = wa2( nint(PINPT%param_const(1,j)) )
              wa2(j) = diag( nint(PINPT%param_const(1,j)) ) * PINPT%param( nint(PINPT%param_const(1,j)) )
            endif
          enddo
          fvec(1:nkpoint) = wa4(1:nkpoint)
          xnorm = enorm ( PINPT%nparam, wa2 )
          fnorm = fnorm1
          iter = iter + 1
          if_main write(6,'(A)')' '
          if_main write(6,'(A,I4,A,F16.6)')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ',fnorm
          if_main write(pfileoutnm_temp,'(A,A)')trim(PINPT%pfileoutnm),'_temp'
          if_main call print_param(PINPT,0,pfileoutnm_temp,.TRUE.)
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

!  End of the outer loop.
     go to 30
300 continue

!  Termination, either normal or user imposed.
  if_main_then
  if(info .eq. 1) then
    write(6,*)' '
    write(6,101)"  Termination INFO=",info,' , condition: |actred|,prered <= ftol, ratio <= 2 :',' actred=',actred, &
                                                                                                 ' prered=',prered, &
                                                                                                 ' ratio=',ratio
  elseif(info .eq. 2) then
    write(6,*)' '
    write(6,102)"  Termination INFO=",info,' , condition: delta <= xtol*xnorm :',' delta=',delta, &
                                                                                 ' xtol*xnorm=',xtol*xnorm
  elseif(info .eq. 3) then
    write(6,*)' '
    write(6,103)"  Termination INFO=",info,' , condition: |actred|,prered <= ftol, delta<=xtol*xnorm :',' actred=',actred, &
                                                                                                        ' prered=',prered, &
                                                                                                        ' delta=',delta, &
                                                                                                        ' xtol*xnorm=',xtol*xnorm
  elseif(info .eq. 5) then
    write(6,*)' '
    write(6,104)"  Termination INFO=",info,' , condition: miter*(nparam+1) <= nfev :',' miter=',PINPT%miter, &
                                                                                      ' nparam=',PINPT%nparam, &
                                                                                      ' nfev=',nfev
  elseif(info .eq. 6) then
    write(6,*)' '
    write(6,105)"  Termination INFO=",info,' , condition: |actred|,prered <= epsmch, ratio <= 2 :',' actred=',actred, & 
                                                                                                   ' prered=',prered, &
                                                                                                   ' epsmch=',epsmch, &
                                                                                                   ' ratio=',ratio
  elseif(info .eq. 7) then
    write(6,*)' '
    write(6,106)"  Termination INFO=",info,' , condition: delta <= epsmch*xnorm :',' delta=',delta, &
                                                                                   ' epsmch*xnorm=',epsmch*xnorm
  elseif(info .eq. 4) then
    write(6,*)' '
    write(6,107)"  Termination INFO=",info,' , condition: gnorm <= epsmch :',' gnorm=',gnorm, &
                                                                             ' epsmch=',epsmch
  else
    write(6,*)' '
    write(6,108)"  Termination INFO=",info
  endif
  if_main_end 

101 format(A,I2,A,3(A,E14.7))
102 format(A,I2,A,2(A,E14.7))
103 format(A,I2,A,4(A,E14.7))
104 format(A,I2,A,3(A,I6   ))
105 format(A,I2,A,4(A,E14.7))
106 format(A,I2,A,2(A,E14.7))
107 format(A,I2,A,2(A,E14.7))
108 format(A,I2)
  if_main call print_param (PINPT, 0, '  Fitted param(i):', .FALSE.)
  return
endsubroutine
subroutine fdjac2 (get_eig, NN_TABLE, kpoint, nkpoint, PINPT, fvec, E_DFT, neig, iband, nband, PWGHT, fjac, epsfcn, flag_get_orbital)
  use parameters
  use mpi_setup
  implicit none
  type (incar  ) :: PINPT
  type (hopping) :: NN_TABLE
  type (weight)  :: PWGHT
  integer*4  nkpoint,neig,iband,nband
  integer*4  i,j
  real*8     eps,epsfcn,epsmch,h,temp,fjac(nkpoint,PINPT%nparam)
  real*8     wa(nkpoint),fvec(nkpoint),E_TBA(nband*PINPT%nspin,nkpoint),E_DFT(neig*PINPT%ispin,nkpoint),kpoint(3,nkpoint)
  complex*16 V(neig*PINPT%ispin,nband*PINPT%nspin,nkpoint)
  logical    flag_get_orbital
  external   get_eig

  epsmch = epsilon(epsmch)
  eps = sqrt( max(epsfcn,epsmch) )
  do j = 1, PINPT%nparam
    temp = PINPT%param(j)
    h = eps*abs(temp)
    if (h == 0.0D+00 ) h=eps
    PINPT%param(j) = temp+h

    call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, E_TBA, V, neig, iband, nband, flag_get_orbital, .false., .false., .true.)
    if(.not. flag_get_orbital) then
      call get_fvec(wa, E_TBA, E_DFT, 0, neig, iband, nband, PINPT, nkpoint, PWGHT)
    else 
      call get_fvec(wa, E_TBA, E_DFT, V, neig, iband, nband, PINPT, nkpoint, PWGHT)
    endif

    PINPT%param(j) = temp
    fjac(1:nkpoint,j) = ( wa(1:nkpoint) - fvec(1:nkpoint) ) / h
  enddo

  return
endsubroutine
