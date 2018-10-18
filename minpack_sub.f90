subroutine lmpar ( nparam, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )
  implicit none
  integer*4 ldr,nparam
  integer*4 i,j,k,l,iter,ipvt(nparam),nsing
  real*8 delta,diag(nparam),dwarf,dxnorm,enorm,gnorm,fp
  real*8 par,parc,parl,paru,sum2,temp
  real*8 qtb(nparam),r(ldr,nparam),sdiag(nparam),wa1(nparam),wa2(nparam),x(nparam)
 
!  DWARF is the smallest positive magnitude.
  dwarf = tiny ( dwarf )
 
!  Compute and store in X the Gauss-Newton direction.
!  If the jacobian is rank-deficient, obtain a least squares solution.
  nsing = nparam

  do j = 1, nparam
    wa1(j) = qtb(j)
    if ( r(j,j) == 0.0D+00 .and. nsing == nparam ) then
      nsing = j - 1
    end if
    if ( nsing < nparam ) then
      wa1(j) = 0.0D+00
    end if
  end do

  do k = 1, nsing
    j = nsing - k + 1
    wa1(j) = wa1(j) / r(j,j)
    temp = wa1(j)
    wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  end do

  do j = 1, nparam
    l = ipvt(j)
    x(l) = wa1(j)
  end do
 
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
  iter = 0
  wa2(1:nparam) = diag(1:nparam) * x(1:nparam)
  dxnorm = enorm ( nparam, wa2 )
  fp = dxnorm - delta

  if ( fp <= 0.1D+00 * delta ) then
    if ( iter == 0 ) then
      par = 0.0D+00
    end if
    return
  end if
 
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function. Otherwise set this bound to zero.
  parl = 0.0D+00

  if ( nparam <= nsing ) then
    do j = 1, nparam
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, nparam
      sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
      wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
    end do

    temp = enorm ( nparam, wa1 )
    parl = ( ( fp / delta ) / temp ) / temp
  end if
 
!  Calculate an upper bound, PARU, for the zero of the function.
  do j = 1, nparam
    sum2 = dot_product ( qtb(1:j), r(1:j,j) )
    l = ipvt(j)
    wa1(j) = sum2 / diag(l)
  end do

  gnorm = enorm ( nparam, wa1 )
  paru = gnorm / delta

  if ( paru == 0.0D+00 ) then
    paru = dwarf / min ( delta, 0.1D+00 )
  end if
 
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
  par = max ( par, parl )
  par = min ( par, paru )
  if ( par == 0.0D+00 ) then
    par = gnorm / dxnorm
  end if
 
!  Beginning of an iteration.
  do
    iter = iter + 1
 
!  Evaluate the function at the current value of PAR.
    if ( par == 0.0D+00 ) then
      par = max ( dwarf, 0.001D+00 * paru )
    end if

    wa1(1:nparam) = sqrt ( par ) * diag(1:nparam)

    call qrsolv ( nparam, r, ldr, ipvt, wa1, qtb, x, sdiag )

    wa2(1:nparam) = diag(1:nparam) * x(1:nparam)
    dxnorm = enorm ( nparam, wa2 )
    temp = fp
    fp = dxnorm - delta
 
!  If the function is small enough, accept the current value of PAR.
    if ( abs ( fp ) <= 0.1D+00 * delta ) then
      exit
    end if
 
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
    if ( parl == 0.0D+00 .and. fp <= temp .and. temp < 0.0D+00 ) then
      exit
    else if ( iter == 10 ) then
      exit
    end if
 
!  Compute the Newton correction.
    do j = 1, nparam
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, nparam
      wa1(j) = wa1(j) / sdiag(j)
      temp = wa1(j)
      wa1(j+1:nparam) = wa1(j+1:nparam) - r(j+1:nparam,j) * temp
    end do

    temp = enorm ( nparam, wa1 )
    parc = ( ( fp / delta ) / temp ) / temp
 
!  Depending on the sign of the function, update PARL or PARU.
    if ( 0.0D+00 < fp ) then
      parl = max ( parl, par )
    else if ( fp < 0.0D+00 ) then
      paru = min ( paru, par )
    end if
 
!  Compute an improved estimate for PAR.
    par = max ( parl, par + parc )
 
!  End of an iteration.
  end do
 
!  Termination.
  if ( iter == 0 ) then
    par = 0.0D+00
  end if

  return
end
subroutine qrfac ( m, n, fjac, ipvt, rdiag, acnorm )
  implicit none

  integer*4 m,n
  integer*4 i,j,k,i4_temp,ipvt(n),kmax,minmn
  real*8 fjac(m,n),acnorm(n),ajnorm,enorm,epsmch
  real*8 r8_temp(m),rdiag(n),temp,wa(n)
  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, fjac(1:m,j) )
  enddo

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  do j = 1, n
    ipvt(j) = j
  end do
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!  Bring the column of largest norm into the pivot position.
      kmax = j
      do k = j, n
        if ( rdiag(kmax) < rdiag(k) ) then
          kmax = k
        end if
      end do

      if ( kmax /= j ) then

        r8_temp(1:m) = fjac(1:m,j)
        fjac(1:m,j)     = fjac(1:m,kmax)
        fjac(1:m,kmax)  = r8_temp(1:m)

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        i4_temp    = ipvt(j)
        ipvt(j)    = ipvt(kmax)
        ipvt(kmax) = i4_temp

      end if
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1, fjac(j,j) )

    if ( ajnorm /= 0.0D+00 ) then

      if ( fjac(j,j) < 0.0D+00 ) then
        ajnorm = -ajnorm
      end if

      fjac(j:m,j) = fjac(j:m,j) / ajnorm
      fjac(j,j) = fjac(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      do k = j + 1, n

        temp = dot_product ( fjac(j:m,j), fjac(j:m,k) ) / fjac(j,j)

        fjac(j:m,k) = fjac(j:m,k) - temp * fjac(j:m,j)

        if (rdiag(k) /= 0.0D+00 ) then

          temp = fjac(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00-temp ** 2 ) )

          if ( 0.05D+00 * ( rdiag(k) / wa(k) ) ** 2 <= epsmch ) then
            rdiag(k) = enorm ( m-j, fjac(j+1,k) )
            wa(k) = rdiag(k)
          end if

        end if

      end do

    end if

    rdiag(j) = - ajnorm

  end do

  return
end
subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )
  implicit none

  integer*4 ldr,n
  integer*4 i,ipvt(n),j,k,l,nsing
  real*8 c,cotan,diag(n)
  real*8 qtb(n),qtbpj,r(ldr,n),s,sdiag(n),sum2,t,temp,wa(n),x(n)
!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
  do j = 1, n
    r(j:n,j) = r(j,j:n)
    x(j) = r(j,j)
  end do

  wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
  do j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
    l = ipvt(j)

    if ( diag(l) /= 0.0D+00 ) then

      sdiag(j:n) = 0.0D+00
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
      qtbpj = 0.0D+00

      do k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
        if ( sdiag(k) /= 0.0D+00 ) then

          if ( abs ( r(k,k) ) < abs ( sdiag(k) ) ) then
            cotan = r(k,k) / sdiag(k)
            s = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
            c = s * cotan
          else
            t = sdiag(k) / r(k,k)
            c = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * t ** 2 )
            s = c * t
          end if
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          do i = k+1, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          end do

        end if

      end do

    end if
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  end do
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  do j = 1, n

    if ( sdiag(j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if

    if ( nsing < n ) then
      wa(j) = 0.0D+00
    end if

  end do

  do j = nsing, 1, -1
    sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
    wa(j) = ( wa(j) - sum2 ) / sdiag(j)
  end do
!
!  Permute the components of Z back to components of X.
!
  do j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
  end do

  return
end
subroutine infostamp (info, lm_method)
  implicit none
  integer*4 info
  character ( len = * ) lm_method

  if ( trim(lm_method) == 'lmdif' ) then

   if(info .eq. 0) then
    write(6,'(A)')"  INFO: 0, improper input parameters"
   elseif(info .eq. 1) then
    write(6,'(A)')"  INFO: 1, both actual and predicted relative reductions in the sum of squares are at most FTOL."
   elseif(info .eq. 2) then
    write(6,'(A)')"  INFO: 2, relative error between two consecutive iterates is at most XTOL."
   elseif(info .eq. 3) then
    write(6,'(A)')"  INFO: 3, conditions for INFO = 1 and INFO = 2 both hold. normal ends"
   elseif(info .eq. 4) then
    write(6,'(A)')"  INFO: 4, the cosine of the angle between FVEC and any column of the jacobian is at most GTOL in absolute value."
   elseif(info .eq. 5) then
    write(6,'(A)')"  INFO: 5, number of calls to FCN has reached or exceeded MAXFEV."
   elseif(info .eq. 6) then
    write(6,'(A)')"  INFO: 6, FTOL is too small.  No further reduction in the sum of squares is possible."
   elseif(info .eq. 7) then
    write(6,'(A)')"  INFO: 7, XTOL is too small.  No further improvement in the approximate solution X is possible."
   elseif(info .eq. 8) then
    write(6,'(A)')"  INFO: 8, GTOL is too small.  FVEC is orthogonal to the columns of the jacobian to machine precision."
   elseif(info .lt. 0) then
    write(6,'(A)')"  INFO: <0, If the user has terminated execution, INFO is set to the (negative) value of IFLAG."
   endif

  endif

end
