! This it the README for the geodesic Levenberg-Marquardt algorithm v1.0
  
! Geodesic Levenberg-Marquardt is a variant of Levenberg-Marquardt that 
! adds third order corrections to the proposed step from the directional 
! second derivative (either by a finite difference estimate or an analytic 
! evaluation).  The main routine is in the file geodesiclm.f90

! The method makes use of the BLAS and LAPACK subroutines for matrix manipulation.
! You should link to these libraries when compiling.

! Although we have used this routine successfully in our own research, 
! we do not guarantee that it is bug free.  
! If you encounter a bug (or an unexpected behavior) please let us know.  
! Send details about the bug to Mark Transtrum: mktranstrum@byu.edu

! If you use this code, please acknowledge such by referencing one one of 
! the following papers in any published work:

! Transtrum M.K., Machta B.B., and Sethna J.P, 
! Why are nonlinear fits to data so challenging?  
! Phys. Rev. Lett. 104, 060201 (2010)

! Transtrum M.K., Machta B.B., and Sethna J.P., 
! The geometry of nonlinear least squares with applications to sloppy model and optimization.  
! Phys. Rev. E. 80, 036701 (2011)

! Disclaimer: No guarantee whatsoever is provided. No liability whatsoever is accepted for 
! any loss or damage of any kind resulting from any defect or inaccuracy in this code.

module lm_geodesic

! MIT License
  
! Copyright (c) 2016 Mark K. Transtrum
 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

contains

SUBROUTINE Acceptance(n,C, Cnew, Cbest, ibold, accepted,  dtd, v, vold)
  IMPLICIT NONE
  INTEGER n,accepted, ibold
  REAL (KIND=8) C, Cnew, Cbest, beta
  REAL (KIND=8) dtd(n,n), v(n), vold(n)

  IF( Cnew .LE. C) THEN !! Accept all downhill steps
     accepted = MAX(accepted + 1, 1)
  ELSE
     !! Calculate beta
     IF (DOT_PRODUCT(vold,vold) .EQ. 0.0d+0) THEN
        beta = 1.0d+0
     ELSE
        beta = DOT_PRODUCT(v,MATMUL(dtd, vold))
        beta = beta/SQRT( DOT_PRODUCT(v,MATMUL(dtd,v)) * DOT_PRODUCT(vold,MATMUL(dtd,vold) ))
        beta = min(1.0d+0,1.0d+0-beta)
     END IF
     SELECT CASE (ibold)
     CASE(0) !! Only downhill steps 
        IF( Cnew .LE. C) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted - 1, -1)
        END IF
     CASE(1)
        IF(beta*Cnew .LE. Cbest) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     CASE(2)
        IF(beta*beta*Cnew .LE. Cbest) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     CASE(3)
        IF(beta*Cnew .LE. C) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     CASE(4)
        IF(beta*beta*Cnew .LE. C) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     END SELECT
  END IF
END SUBROUTINE Acceptance

SUBROUTINE convergence_check(m, n, converged, accepted, counter, C, Cnew, x, fvec, fjac, lam, xnew, &
       nfev, maxfev, njev, maxjev, naev, maxaev, maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,cos_alpha)
  IMPLICIT NONE
  INTEGER m,n, converged, accepted, counter, nfev, maxfev, njev, maxjev, naev, maxaev
  REAL (KIND=8) C, Cnew, x(n), fvec(m), fjac(m,n), xnew(n), grad(n), lam, rpar(m)
  REAL (KIND=8) maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol, cos_alpha
  INTEGER i

!  The first few criteria should be checked every iteration, since
!  they depend on counts and the Jacobian but not the proposed step.

!  nfev
  IF(maxfev .GT. 0) THEN
     IF(nfev.GE.maxfev) THEN
        converged = -2
        counter = 0
        RETURN
     ENDIF
  ENDIF

!  njev
  IF(maxjev .GT. 0) THEN
     IF(njev.GE.maxjev) THEN
        converged = -3
        RETURN
     ENDIF
  ENDIF


!  naev
  IF(maxaev .GT. 0) THEN
     IF(naev.GE.maxaev) THEN
        converged = -4
        RETURN
     END IF
  ENDIF

!  maxlam
  IF(maxlam .GT. 0.0d+0) THEN
     IF(lam .GE. maxlam) THEN
        converged = -5
        RETURN
     END IF
  END IF

!  minlam
  IF(minlam .GT. 0.0d+0 .AND. lam .GT. 0.0d+0) THEN
     IF(lam .LE. minlam) THEN
        counter = counter + 1
        IF(counter .GE. 3) THEN
           converged = -6
           RETURN
        END IF
        RETURN
     END IF
  END IF

! artol -- angle between residual vector and tangent plane
  IF( artol .GT. 0.0d+0) THEN
     !! Only calculate the projection if artol > 0
     !! CALL projection(m,n,fvec, fjac, rpar,eps)
     !! cos_alpha = SQRT(DOT_PRODUCT(rpar,rpar)/DOT_PRODUCT(fvec,fvec))
     IF( cos_alpha .LE. artol) THEN
        converged = 1
        RETURN
     END IF
  END IF

! If gradient is small
  grad = -1.0d+0*MATMUL(fvec, fjac)
  IF(SQRT(DOT_PRODUCT(grad,grad)).LE.gtol) THEN
     converged = 3
     RETURN
  ENDIF

! If cost is sufficiently small
  IF (C .LT. Cgoal) THEN !! Check every iteration in order to catch a cost small on the first iteration
     converged = 2
     RETURN
  END IF


!  If step is not accepted, then don't check remaining criteria
  IF(accepted.LT.0) THEN
     counter = 0
     converged = 0
     RETURN
  ENDIF

! If step size is small
  if(SQRT(DOT_PRODUCT(x-xnew,x-xnew)).LT. xtol) THEN
     converged = 4
     RETURN
  ENDIF

! If each parameter is moving relatively small
  xrtolcheck: DO i = 1,n
     converged = 5
     IF(  ABS(x(i) - xnew(i)) .GT. xrtol*ABS(x(i)) .OR.  (xnew(i) .NE. xnew(i)) ) converged = 0 !! continue if big step or nan in xnew
     IF( converged .EQ. 0) EXIT xrtolcheck
  END DO xrtolcheck
  IF(converged .EQ. 5) RETURN

! If cost is not decreasing -- this can happen by accident, so we require that it occur three times in a row
  IF( (C - Cnew).LE. ftol .AND.(C-Cnew).GE.0.) THEN
     counter = counter + 1
     IF(counter .GE. 3) THEN
        converged = 6
        RETURN
     END IF
     RETURN
  ENDIF

! If cost is not decreasing relatively -- again can happen by accident so require three times in a row
  IF( (C - Cnew).LE.(frtol*C).AND.(C-Cnew).GE.0.) THEN
     counter = counter + 1
     IF(counter .GE. 3) THEN
        converged = 7
        RETURN
     ENDIF
     RETURN
  ENDIF

! If none of the above: continue
  counter = 0
  converged = 0
END SUBROUTINE convergence_check

SUBROUTINE FDAvv(m,n,x,v,fvec,fjac, func,acc, jac_uptodate, h2)
  IMPLICIT NONE
  INTEGER m, n
  REAL (KIND=8) x(n), v(n), fvec(m), fjac(m,n), acc(m), xtmp(n), ftmp(m), h2
  LOGICAL jac_uptodate
  EXTERNAL func

  IF( jac_uptodate) THEN
     xtmp = x + h2*v
     CALL func(m,n,xtmp,ftmp)
     acc = (2.0d+0/h2)*( (ftmp - fvec)/h2 - MATMUL(fjac,v) )
  ELSE !if jacobian not up to date, do not use jacobian in F.D. (needs one more function call)
     xtmp = x + h2*v
     CALL func(m,n,xtmp,ftmp)
     xtmp = x - h2*v
     CALL func(m,n,xtmp,acc)
     acc = (ftmp - 2*fvec + acc)/(h2*h2)
  ENDIF
END SUBROUTINE FDAvv
! -*- f90 -*-            
! ****************************************
! Routine for calculating finite-difference jacobian

SUBROUTINE FDJAC(m,n,x,fvec,fjac,func,eps,center_diff)
  IMPLICIT NONE
  INTEGER m,n, i
  REAL (KIND=8) x(n), dx(n), fvec(m), fjac(m,n), eps, epsmach, dpmpar
  LOGICAL center_diff
  REAL (KIND=8) h, temp1(m), temp2(m) 
  
  epsmach = dpmpar(1)
  IF(center_diff) THEN
     DO i = 1, n
        h = eps*ABS(x(i))
        IF (h < epsmach) h = eps
        dx(:) = 0.0D+00
        dx(i) = 0.5d+0*h
        CALL func(m,n,x+dx,temp1,0)
        CALL func(m,n,x-dx,temp2,0)
        fjac(:,i) = (temp1 - temp2)/h
     END DO
  ELSE
     DO i = 1, n
        h = eps*ABS(x(i))
        IF (h < epsmach) h = eps
        dx(:) = 0.0D+00
        dx(i) = h
        CALL func(m,n,x+dx,temp1,0)
        fjac(:,i) = (temp1 - fvec)/h
     END DO
  END IF
 
END SUBROUTINE FDJAC
! -*- f90 -*-
! file leastsq.f90
! Main Geodesic-Bold-BroydenUpdate-Levenberg-Marquardt routine
! version 1.0.2

SUBROUTINE geodesiclm(func, jacobian, Avv, &
       x, fvec, fjac, n, m, &
       callback, info, &
       analytic_jac, analytic_Avv, &
       center_diff, h1, h2,&
       dtd, damp_mode, &
       niters, nfev, njev, naev, &
       maxiter, maxfev, maxjev, maxaev, maxlam, minlam, &
       artol, Cgoal, gtol, xtol, xrtol, ftol, frtol, &
       converged, &
       print_level, print_unit, &
       imethod, iaccel, ibold, ibroyden, &
       initialfactor, factoraccept, factorreject, avmax)

!*****************************************************************
! 
!    subroutine geodesicLM
!    
!    The purpose of geolevmar is to minimize the sum of the squares
!    of m nonlinear functions of n variables by a modification of
!    the Levenberg-Marquardt algorithm that utilizes the geodesic
!    acceleration step correction, bold acceptance criterion, and
!    a Broyden update of the jacobian matrix.  The method employs one
!    of several possible schemes for updating the Levenberg-Marquardt
!    parameter.  The user must provide a subroutine which calcualtes
!    the functions, and optionally the jacobian and a directional 
!    derivative of the functions.  The latter two will be estimated
!    by finite differences if not supplied.
!
!    If you use this code, please acknowledge such by referencing one
!    one of the following papers in any published work:
!    
!    Transtrum M.K., Machta B.B., and Sethna J.P, Why are nonlinear
!    fits to data so challenging?  Phys. Rev. Lett. 104, 060201 (2010)
!
!    Transtrum M.K., Machta B.B., and Sethna J.P., The geometry of
!    nonlinear least squares with applications to sloppy model and
!    optimization.  Phys. Rev. E. 80, 036701 (2011)
!
!
!    The subroutine statement is:
!
!    geodesicLM(func, jacobian, Avv, x, fvec, fjac, n, m, callback, info,
!              analytic_jac, analytic_Avv, center_diff, h1, h2,
!              dtd, damp_mode, niteres, nfev, njev, naev,
!              maxiters, maxfev, maxjev, maxaev, maxlam, minlam,
!              artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,
!              converged, print_level, print_unit,
!              imethod, iaccel, ibold, ibroyden,
!              initialfactor, factoraccept, factorreject, avmax)
!
!    where
!
!    func is a user supplied subroutine which calculates the functions and
!    should be written as follows:
!
!      subroutine func(m, n, x, fvec)
!      integer m, n
!      double precision x(n), fvec(m)
!      --------------------------------------------------------------------
!      calculates the function at x and returns their values in fvec
!      x, m, and n should be left unchanged
!      --------------------------------------------------------------------
!      end subroutine func
!
!    jacobian is a user supplied subroutine which calculates the jacobian of
!    of the functions if analytic_jac is .TRUE.
!    jacobian should be writen as follows     
!
!      subroutine jacobian(m, n, x, fjac)
!      integer m, n
!      double precision x(n), fjac(m,n)
!      --------------------------------------------------------------------
!      calculates the jacobian at x and returns their values in fjac
!      x, m, and n should be left unchanged
!      --------------------------------------------------------------------
!      end subroutine jacobian
!
!    Avv is a user supplied subroutine which calculates the directional
!    second derivative of the functions if analytic_Avv is .TRUE.
!    Avv should be writen as follows     
!
!      subroutine Avv(m, n, x, v, acc)
!      integer m, n
!      double precision x(n), v(n), acc(m)
!      --------------------------------------------------------------------
!      calculates the directional second derivative at x in the direction 
!      of v and returns the values in acc
!      x, v, m, and n should be left unchanged
!      --------------------------------------------------------------------
!      end subroutine Avv
!
!    x is an array of length n.  On input it contains an initial estimate of
!    the solution.  On exit, it contains the final estimate of the solution.
!
!    fvec is an output array of length m containing the funtion evaluation at
!    the final solution
!
!    fjac is an output array of dimension(m,n) containing the jacobian evaluation
!    the final solution.  The array MATMUL( TRANSPOSE(fjac), fjac) is an estimate
!    of the covariance matrix of parameters at the final solution
!
!    n an input integer set to the number of parameters
!
!    m an input integer set to the number of functions
!
!    callback a user supplied subroutine to be called after each iteration of the
!    algorithm.  
!    callback should be written as follows
!
!      subroutine callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
!      integer m, n, accepted, info
!      double precision x(n), v(n), a(n), fvec(m), fjac(m,n), acc(m), lam, dtd(n,n), fvec_new(m)
!      --------------------------------------------------------------------
!      m, n, x, v, a, fvec, fjac, acc, lam, dtd, fvec_new, accepted, should be left unchanged
!      On input, info = 0 and should be changed to a nonzero value if the user
!      wishes to terminate calculation
!      --------------------------------------------------------------------
!      end subroutine callback
!
!    info an output integer set to a nonzero value if the user terminated the routine
!    (see callback).
!
!    analytic_jac an input boolean set to .TRUE. if the subroutine jacobian calculates
!    the jacobian.  If .FALSE. then a finite difference estimate will be used.
!
!    analytic_Avv an input boolean set to .TRUE. if the subroutine Avv calculates
!    the directional second derivative.  If .FALSE. then a finite difference estimate
!    will be used.
!
!    center_diff an input boolean.  If finite differences are used to estimate the jacobian
!    then center differences will used if center_diff is .TRUE., otherwise, forward
!    differences will be used.  Note that center differences are more accurate by require
!    more function evaluations.
!
!    h1 an input double precision specifying the step size for the finite difference estimates
!    of the jacobian.
!
!    h2 an input double precision specifying the steps ize for the finite difference estiamtes
!    of the directional second derivative.
!
!    dtd a double precision array of dimension(n,n).  dtd is used as the damping matrix in the 
!    Levenberg-Marquardt routine.  It's exact treatment is specified by the damp_mode input.
!
!    damp_mode an input integer specifying the details of the LM damping as follows:
!      damp_mode = 0: dtd is set to the identity.
!      damp_mode = 1: dtd should be a positive definite, diagonal matrix whose entries are dynamically
!                updated based on the elements of the jacobian.
!
!    niters an output integer specifying the number of iterations of the algorithm.
!
!    nfev an output integer specifying the number of calls to func.  
!
!    njev an output integer specifying the number of calls to jacobian.
!
!    naev an output integer specifying the number of calls to Avv.
!
!    maxiter an input integer specifying the maximum number of routine iterations.
!
!    maxfev an input integer specifying the maximum number of function calls
!    if maxfev = 0, then there is no limit to the number of function calls.
!
!    maxjev an input integer specifying the maximum number of jacobian calls
!    if maxjev = 0, then there is no limit to the number of jacobian calls.
!
!    maxaev an input integer specifying the maximum number of Avv calls
!    if maxaev = 0, then there is no limit to the number of Avv calls.
!
!    maxlam an input double precision specifying the maximum allowed value of 
!    the damping term lambda. If this is negative, then there is no limit.
!
!    minlam an input double precision specifying the minimum allowed value of 
!    the damping term lambda. If lambda is smaller than this value for three consecutive steps
!    the routine terminates.  If this is negative, then there is no limit.
!
!    artol an input double precision.  The method will terminate when the cosine of the
!    angle between the residual vector and the range of the jacobian is less than artol.
!
!    Cgoal an input double precision.  The method will terminate when the cost (one half
!    the sum of squares of the function) falls below Cgoal.
!
!    gtol an input double precision.  The method will terminate when norm of Cost gradient 
!    falls below gtol.
!    
!    xtol an input double precision.  The method will terminate when parameters change by
!    less than xtol.
!
!    xrtol an input double precision.  The method will terminate if the relative change in
!    each of the parameters is less than xrtol.
!
!    ftol an input double precision.  The method will termiante if the Cost fails to decrease
!    by more than ftol for 3 consecutive iterations.
!
!    frtol an input double precision.  The method will terminate if the relative decrease in
!    Cost is less than frtol 3 consecutive iterations.
!
!    converged an output integer indicated the reason for termination:
!      converged = 1: artol 
!      converged = 2: Cgoal
!      converged = 3: gtol
!      converged = 4: xtol
!      converged = 5: xrtol
!      converged = 6: ftol
!      converged = 7: frtol
!      converged = -1: maxiters exeeded
!      converged = -2: maxfev exceeded
!      converged = -3: maxjev exceeded
!      converged = -4: maxaev exceeded
!      converged = -10: user requested termination in callback via info
!      converged = -11: Either the initial function evalaution or subsequent jacobian
!                       evaluations produced Nans.
!
!    print_level an input integer specifying the amount of details to be printed.
!    acceptable values range from 0 to 5, with larger number printing more details.
!
!    print_unit an input integer specifying the unit number details should be written to.
!
!    imethod an input integer specifying the method for updating the LM parameter
!      imethod = 0: adjusted by fixed factors after accepted/rejected steps
!      imethod = 1: adjusted as described in Nielson
!      imethod = 2: adjusted according to an unpublished method due to Cyrus Umrigar and Peter Nightingal
!      imethod = 10: step size Delta adjusted by fixed factors after accepted/rejected steps
!      imethod = 11: step size adjusted as described in More'
!
!    initialfactor an input double precision for specifying either the initial LM parameter
!    of the initial step size.
!
!    factoraccept an input double precision (larger than 1.0) specifying the factor by which
!    either the LM parameter or the step size will be adjusted after an accepted step if
!    imethod = 0 or 10
!
!    factorreject an input double precision (larger than 1.0) specifying the factor by which
!    either the LM parameter of the step size will be adjusted after a rejected step if
!    imethod = 0 or 10
!
!    avmax an input double precision specifying the maximum norm of the geodesic acceleration 
!    relative to the velocity vector.
!
!*****************************************************************

  IMPLICIT NONE
  !! Passed parameters
  EXTERNAL func, jacobian, Avv, callback

  REAL (KIND=8) x(n), fvec(m), fjac(m,n)
  INTEGER n, m
  LOGICAL analytic_jac, analytic_Avv, center_diff
  REAL (KIND=8) h1, h2
  REAL (KIND=8) dtd(n,n)
  INTEGER damp_mode, info
  INTEGER niters, nfev, njev, naev
  INTEGER maxiter, maxfev, maxaev, maxjev
  REAL (KIND=8) maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol
  INTEGER converged
  INTEGER print_level, print_unit
  INTEGER iaccel, ibroyden, ibold, imethod
  REAL (KIND=8) avmax, initialfactor, factoraccept, factorreject

  !! Internal parameters

  REAL (KIND=8) acc(m), v(n), vold(n), a(n), lam, delta, cos_alpha, av
  REAL (KIND=8) fvec_new(m), fvec_best(m), C, Cnew, Cbest, Cold
  REAL (KIND=8) x_new(n), x_best(n)
  REAL (KIND=8) jtj(n,n), g(n,n)
  REAL (KIND=8) temp1, temp2, pred_red, dirder, actred, rho, a_param
  INTEGER i, j, istep, accepted, counter
  
  character(16) :: converged_info(-11:7)
  LOGICAL jac_uptodate, jac_force_update, valid_result

  ! strings for concluding print statement
  converged_info = '????????'
  converged_info(1) = 'artol reached'
  converged_info(2) = 'Cgoal reached'
  converged_info(3) = 'gtol reached'
  converged_info(4) = 'xtol reached'
  converged_info(5) = 'xrtol reached'
  converged_info(6) = 'ftol reached'
  converged_info(7) = 'frtol reached'
  converged_info(-1) = 'maxiters exeeded'
  converged_info(-2) = 'maxfev exceeded'
  converged_info(-3) = 'maxjev exceeded'
  converged_info(-4) = 'maxaev exceeded'
  converged_info(-10) = 'User Termination '
  converged_info(-11) = 'NaN Produced'

  IF(print_level .GE. 1) THEN
     WRITE(print_unit, *) "Optimizing with Geodesic-Levenberg-Marquardt algorithm, version 1.0.2"
     WRITE(print_unit, *) "Method Details:"
     WRITE(print_unit, *) "  Update method:   ", imethod
     WRITE(print_unit, *) "  acceleration:    ", iaccel
     WRITE(print_unit, *) "  Bold method:     ", ibold
     WRITE(print_unit, *) "  Broyden updates: ", ibroyden
     FLUSH(print_unit)
  ENDIF

  !! Initialize variables
  niters = 0
  nfev = 0
  naev = 0
  njev = 0
  converged = 0
  v(:) = 0.0d+0
  vold(:) = 0.0d+0
  a(:) = 0.0d+0
  cos_alpha = 1.0d+0
  av = 0.0d+0
  a_param = 0.5

  accepted = 0
  counter = 0
  CALL func(m,n,x,fvec)
  nfev = nfev + 1
  C = 0.5d+0*DOT_PRODUCT(fvec,fvec)
  IF(print_level .GE. 1) THEN
     WRITE(print_unit, *) "  Initial Cost:    ", C
     FLUSH(print_unit)
  ENDIF
  valid_result = .TRUE.
  !! Check for nans in fvec
  checkfvec: DO i = 1,m     
     IF(fvec(i) /= fvec(i)) THEN
        valid_result = .FALSE.
        EXIT checkfvec
     END IF
  END DO checkfvec
  IF (.NOT. valid_result) THEN
     converged = -11
     maxiter = 0
  ENDIF
  Cbest = C
  fvec_best = fvec
  x_best = x
  IF(analytic_jac) THEN
     CALL jacobian(m,n,x,fjac)
     njev = njev + 1
  ELSE 
     CALL fdjac(m,n,x,fvec,fjac,func,h1,center_diff)
     IF (center_diff) THEN
        nfev = nfev + 2*n
     ELSE
        nfev = nfev + n
     ENDIF
  ENDIF
  jac_uptodate = .TRUE.
  jac_force_update = .FALSE.
  jtj = MATMUL(TRANSPOSE(fjac), fjac)

  !! Check fjac for nans
  valid_result = .TRUE.
  checkfjac_initial: DO i = 1,m     
     DO j = 1,n
        IF(fjac(i,j) /= fjac(i,j)) THEN
           valid_result = .FALSE.
           EXIT checkfjac_initial
        END IF
     END DO
  END DO checkfjac_initial
  IF( .NOT. valid_result) THEN
     converged = -11
     maxiter = 0
  ENDIF

  acc(:) = 0.0d+0
  a(:) = 0.0d+0

  !! Initialize scaling matrix
  IF(damp_mode.EQ.0) THEN
     dtd(:,:) = 0.0d+0
     DO i=1,n
        dtd(i,i) = 1.0d+0
     END DO
  ELSEIF(damp_mode.EQ.1) THEN
     DO i = 1,n
        dtd(i,i) = MAX(jtj(i,i),dtd(i,i))
     END DO
  ENDIF

  !! Initialize lambda
  IF(imethod .LT. 10) THEN
     lam = jtj(1,1)
     DO i = 2,n
        lam = MAX(jtj(i,i),lam)
     END DO
     lam = lam * initialfactor
  !! Initialize step bound if using trust region method
  ELSEIF(imethod .GE. 10) THEN
     delta = initialfactor*SQRT(DOT_PRODUCT(x,MATMUL(dtd,x)))
     lam = 1.0d+0
     IF(delta .EQ. 0.0d+0) delta = 100d+0
     IF( converged .EQ. 0) CALL TrustRegion(n,m,fvec,fjac,dtd,delta,lam) !! Do not call this if there were nans in either fvec or fjac
  ENDIF

  !! Main Loop
  main: DO istep=1, maxiter
     
     info = 0
     CALL callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
     IF( info .NE. 0) THEN
        converged = -10
        exit main
     ENDIF
     !! Update Functions

     !! Full or partial Jacobian Update?
     IF (accepted .GT. 0 .AND. ibroyden .LE. 0) jac_force_update = .TRUE.
     IF (accepted + ibroyden .LE. 0 .AND. .NOT. jac_uptodate) jac_force_update = .TRUE.  !Force jac update after too many failed attempts

     IF (accepted .GT. 0 .AND. ibroyden .GT. 0 .AND. .NOT. jac_force_update) THEN !! Rank deficient update of jacobian matrix
        CALL UPDATEJAC(m,n,fjac, fvec, fvec_new, acc, v, a)
        jac_uptodate = .FALSE.
     ENDIF

     IF( accepted .GT. 0) THEN !! Accepted step
        fvec = fvec_new
        x = x_new
        vold = v
        C = Cnew
        IF( C .LE. Cbest) THEN
           x_best = x
           Cbest = C
           fvec_best = fvec
        ENDIF
     ENDIF
     

     IF( jac_force_update ) THEN !! Full rank update of jacobian
        IF(analytic_jac) THEN
           CALL jacobian(m,n,x,fjac)
           njev = njev + 1
        ELSE
           CALL fdjac(m,n,x,fvec,fjac,func,h1,center_diff)
           IF (center_diff) THEN
              nfev = nfev + 2*n
           ELSE
              nfev = nfev + n
           ENDIF
        ENDIF
        jac_uptodate = .TRUE.
        jac_force_update = .FALSE.
     ENDIF

     !! Check fjac for nans
     valid_result = .TRUE.
     checkfjac: DO i = 1,m     
        DO j = 1,n
           IF(fjac(i,j) /= fjac(i,j)) THEN
              valid_result = .FALSE.
              EXIT checkfjac
           END IF
        END DO
     END DO checkfjac

     IF (valid_result) THEN !! If no nans in jacobian

        jtj = MATMUL(TRANSPOSE(fjac), fjac)
     
        !! Update Scaling/lam/TrustRegion
        IF(istep .GT. 1) THEN !! Only necessary after first step
           IF( damp_mode .EQ. 1) THEN
              DO i = 1,n
                 dtd(i,i) = MAX(jtj(i,i),dtd(i,i))
              END DO
           ENDIF
           !! Can add other lam-delta update methods 
           SELECT CASE( imethod )
           CASE( 0 )
              !! Update lam directly by fixed factors
              CALL Updatelam_factor(lam, accepted, factoraccept, factorreject)
           CASE( 1 )
              !! Update lam directly based on Gain Factor rho (see Nielson reference)
              CALL Updatelam_nelson(lam, accepted, factoraccept, factorreject, rho)
           CASE( 2 )
              !! Update lam directly using method of Umrigar and Nightingale [unpublished]
              CALL Updatelam_Umrigar(m,n,lam,accepted, v, vold, fvec, fjac, dtd, a_param, Cold, Cnew) 
           CASE(10)
              !! Update delta by fixed factors
              CALL UpdateDelta_factor(delta, accepted, factoraccept, factorreject)
              CALL TrustRegion(n,m,fvec, fjac, dtd, delta, lam)
           CASE(11)
              !! Update delta as described in More' reference
              CALL UpdateDelta_more(delta, lam, n, v, dtd, rho, C, Cnew, dirder, actred, av, avmax)
              CALL TrustRegion(n,m,fvec, fjac, dtd, delta, lam)
           END SELECT
        ENDIF

        !! Propose Step
        !! metric aray
        g = jtj + lam*dtd
        !! Cholesky decomposition
        CALL DPOTRF('U', n, g, n, info)
        !! CALL inv(n, g, info)
     ELSE !! If nans in jacobian
        converged = -11
        exit main
     ENDIF


     IF(info .EQ. 0) THEN  !! If matrix decomposition successful:
        !! v = -1.0d+0*MATMUL(g,MATMUL(fvec,fjac)) ! velocity
        v = -1.0d+0*MATMUL(fvec, fjac)
        CALL DPOTRS('U', n, 1, g, n, v, n, info)
        
        ! Calcualte the predicted reduction and the directional derivative -- useful for updating lam methods
        temp1 = 0.5d+0*DOT_PRODUCT(v,MATMUL(jtj, v))/C
        temp2 = 0.5d+0*lam*DOT_PRODUCT(v,MATMUL(dtd,v))/C
        pred_red = temp1 + 2.0d+0*temp2
        dirder = -1.0d+0*(temp1 + temp2)
        ! calculate cos_alpha -- cos of angle between step direction (in data space) and residual vector
        cos_alpha = ABS(DOT_PRODUCT(fvec, MATMUL(fjac, v)))
        cos_alpha = cos_alpha/SQRT( DOT_PRODUCT(fvec, fvec)*DOT_PRODUCT(MATMUL(fjac,v), MATMUL(fjac, v)) )
        IF ( imethod .LT. 10) delta = SQRT(DOT_PRODUCT(v, MATMUL(dtd, v)))  !! Update delta if not set directly
        !! update acceleration
        IF(iaccel .GT. 0) THEN
           IF( analytic_Avv ) THEN 
              CALL Avv(m,n,x,v,acc)
              naev = naev + 1
           ELSE
              CALL FDAvv(m,n,x,v,fvec, fjac, func, acc, jac_uptodate, h2)
              IF(jac_uptodate) THEN
                 nfev = nfev + 1
              ELSE 
                 nfev = nfev + 2 !! we don't use the jacobian if it is not up to date
              ENDIF
           ENDIF
           !! Check accel for nans
           valid_result = .TRUE.
           checkAccel: DO i = 1,m     
              IF(acc(i) /= acc(i)) THEN
                 valid_result = .FALSE.
                 EXIT checkAccel
              END IF
           END DO checkAccel
           IF (valid_result ) THEN
              a = -1.0d+0*MATMUL(acc, fjac)
              CALL DPOTRS('U', n, 1, g, n, a, n, info)
              !!a = -1.0d+0*MATMUL(g,MATMUL(acc,fjac))
           ELSE 
              a(:) = 0.0d+0 !! If nans in acc, we will ignore the acceleration term
           ENDIF
        ENDIF

        !! Evaluate at proposed step -- only necessary if av <= avmax
        av = SQRT(DOT_PRODUCT(a,MATMUL(dtd, a))/DOT_PRODUCT(v,MATMUL(dtd,v)))
        IF( av .LE. avmax) THEN
           x_new = x + v + 0.5d+0*a
           CALL func(m,n,x_new,fvec_new)
           nfev = nfev + 1
           Cnew = 0.5d+0*DOT_PRODUCT(fvec_new,fvec_new)
           Cold = C
           valid_result = .TRUE.
           !! Check for nans in fvec_new
           checkfvec_new: DO i = 1,m     
              IF(fvec_new(i) /= fvec_new(i)) THEN
                 valid_result = .FALSE.
                 EXIT checkfvec_new
              END IF
           END DO checkfvec_new
           IF (valid_result) THEN  !! If no nans, proceed as normal
              ! update rho and actred
              actred = 1.0d+0 - Cnew/C
              rho = 0.0d+0
              IF(pred_red .NE. 0.0d+0) rho = (1.0d+0 - Cnew/C)/pred_red
              !! Accept or Reject proposed step
              CALL Acceptance(n,C, Cnew, Cbest, ibold, accepted, dtd, v, vold)
           ELSE !! If nans in fvec_new, reject step
              actred = 0.0d+0
              rho = 0.0d+0
              accepted = MIN(accepted - 1, -1)
           ENDIF
        ELSE !! If acceleration too large, then reject
           accepted = MIN(accepted - 1, -1)
        ENDIF
     ELSE !! If matrix factorization fails, reject the proposed step
        accepted = MIN(accepted -1, -1)
     ENDIF

     !! Check Convergence
     IF (converged .EQ. 0) THEN
        CALL convergence_check(m, n, converged, accepted, counter, &
             & C, Cnew, x, fvec, fjac, lam, x_new, &
             & nfev, maxfev, njev, maxjev, naev, maxaev, maxlam, minlam, &
             & artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,cos_alpha)
        IF (converged .EQ. 1 .AND. .NOT. jac_uptodate) THEN  
           !! If converged by artol with an out of date jacobian, update the jacoban to confirm true convergence
           converged = 0
           jac_force_update = .TRUE.
        END IF
     ENDIF
     

     !! Print status
     IF (print_level .EQ. 2 .AND. accepted .GT. 0) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        FLUSH(print_unit)
     ELSEIF(print_level .EQ. 3) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        FLUSH(print_unit)
     ENDIF
     IF (print_level .EQ. 4 .AND. accepted .GT. 0) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        WRITE(print_unit, *) "  x = ", x
        WRITE(print_unit, *) "  v = ", v
        WRITE(print_unit, *) "  a = ", a
        FLUSH(print_unit)
     ELSEIF (print_level .EQ. 5) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        WRITE(print_unit, *) "  x = ", x
        WRITE(print_unit, *) "  v = ", v
        WRITE(print_unit, *) "  a = ", a
        FLUSH(print_unit)
     ENDIF

     ! If converged -- return
     IF(converged .NE. 0) THEN
        exit main
     ENDIF

     IF (accepted .GE. 0) jac_uptodate = .FALSE. !jacobian is now out of date
     
  END DO main
  ! end main loop
  
  ! If not converged
  IF(converged .EQ. 0) converged = -1
  niters = istep

  ! Return best fit found
  ! If the method converged, but final x is different from x_best -- what to do?
  x = x_best
  fvec = fvec_best

  IF(print_level .GE. 1) THEN
     WRITE(print_unit,*) "Optimization finished"
     WRITE(print_unit,*) "Results:"
     WRITE(print_unit,*) "  Converged:    ", converged_info(converged), converged
     WRITE(print_unit,*) "  Final Cost: ", 0.5d+0*DOT_PRODUCT(fvec,fvec)
     WRITE(print_unit,*) "  Cost/DOF: ", 0.5d+0*DOT_PRODUCT(fvec,fvec)/(m-n)
     WRITE(print_unit,*) "  niters:     ", istep
     WRITE(print_unit,*) "  nfev:       ", nfev
     WRITE(print_unit,*) "  njev:       ", njev
     WRITE(print_unit,*) "  naev:       ", naev
     FLUSH(print_unit)
  ENDIF

END SUBROUTINE geodesiclm
     

     



! -*- f90 -*-
!****************************************
! Routines for updating lam

SUBROUTINE TrustRegion(n,m, fvec, fjac, dtd, delta, lam)
  !! Calls dgqt supplied by minpack to calculate the step and Lagrange multiplier
  IMPLICIT NONE
  INTEGER n, m, i, itmax, info
  REAL (KIND=8) fvec(m), fjac(m,n), dtd(n,n), delta, lam, v(n)
  REAL (KIND=8) z(n), wa1(n), wa2(n) !! Work arrays for dgqt
  REAL (KIND=8) rtol, atol, f
  REAL (KIND=8) jtilde(m,n), gradCtilde(n), g(n,n)

  !! Parameters for dgqt
  rtol = 1.0d-03
  atol = 1.0d-03 
  itmax = 10

  DO i = 1,n
     jtilde(:,i) = fjac(:,i)/SQRT(dtd(i,i)) !! This assumes that dtd is diagonal...
  END DO
  gradCtilde = MATMUL(fvec, jtilde)
  g = MATMUL(TRANSPOSE(jtilde), jtilde)
  CALL dgqt(n, g, n, gradCtilde, delta, rtol, atol, itmax, lam, f, v, info, i, z, wa1, wa2)
  !! Transform v back to non-dtd units
!!$  DO i = 1,n
!!$     v(i) = v(i)/SQRT(dtd(i,i))
!!$  END DO
  RETURN
END SUBROUTINE TrustRegion
!! traditional update methods

SUBROUTINE Updatelam_factor(lam, accepted, factoraccept, factorreject)
  !! Update lam based on accepted/rejected step
  IMPLICIT NONE
  INTEGER accepted
  REAL (KIND=8) lam, factoraccept, factorreject

  IF(accepted.GE.0) THEN
     lam = lam / factoraccept
  ELSE
     lam = lam * factorreject
  ENDIF
END SUBROUTINE Updatelam_factor

SUBROUTINE Updatelam_nelson(lam, accepted, factoraccept, factorreject, rho)
  !! Update method due to Nelson [ref]
  IMPLICIT NONE
  INTEGER accepted, i
  DOUBLE PRECISION lam, factoraccept, factorreject, nu, rho
  IF(accepted.GE.0) THEN
     lam = lam * MAX( 1.0d+0/factoraccept, 1.0d+0 - (factorreject - 1.0d+0)*(2.0d+0*rho - 1.0d+0)**3 )
  ELSE
     nu = factorreject
     DO i = 2,-1*accepted !! double nu for each rejection
        nu = nu*2.0d+0
     END DO
     lam = lam * nu
  ENDIF
END SUBROUTINE Updatelam_nelson


SUBROUTINE Updatelam_Umrigar(m,n,lam, accepted, v, vold, fvec, fjac, dtd, a_param,C,Cnew)
  !! Method due to Umrigar and Nightingale [unpublished]
  IMPLICIT NONE
  INTEGER n,m,accepted, info
  DOUBLE PRECISION lam,v(n),vold(n),fvec(m), fjac(m,n), g(n,n), dtd(n,n), C, Cnew
  DOUBLE PRECISION lamold, d1, d2, grad(n), factor, a_param
  DOUBLE PRECISION amemory, cos_on

  amemory = EXP(-1.0d+0/5.0d+0)
  cos_on = DOT_PRODUCT(v,MATMUL(dtd, vold))
  cos_on = cos_on/SQRT( DOT_PRODUCT(v, MATMUL(dtd,v))*DOT_PRODUCT(vold, MATMUL(dtd,vold)))
  IF( accepted.GE.0) THEN

     IF( Cnew .LE. C ) THEN
        IF( cos_on .GT. 0) THEN
           a_param = amemory*a_param + 1.0 - amemory
        ELSE
           a_param =  amemory*a_param + 0.5*(1.0 - amemory)
        END IF
     ELSE
        a_param =  amemory*a_param + 0.5*(1.0 - amemory)
     END IF

     factor = MIN( 100.0d+0, MAX(1.1d+0, 1.0d+0/(2.2D-16 + 1.0d+0-ABS(2.0d+0*a_param - 1.0d+0))**2))
     IF (Cnew .LE. C .AND. cos_on .GE. 0) THEN
        lam = lam/factor
     ELSEIF(Cnew .GT. C) THEN
        lam = lam*SQRT(factor)
     END IF

  ELSE
     a_param =  amemory*a_param
     factor = MIN( 100.0d+0, MAX(1.1d+0, 1.0d+0/(2.2D-16 + 1.0d+0-ABS(2.0d+0*a_param - 1.0d+0))**2))
     lamold = lam
     IF( cos_on .GT. 0 ) THEN
        lam = lam * SQRT(factor)
     ELSE
        lam = lam * factor
     END IF
     ! Check for a 10% change in drift
     ! Umrigar and Nightingal suggest a check that the the proposed change in lam actually produces a meaningful change in the step.
     ! But this code produces strange results in a few cases.  -MKT
!!$     IF(accepted .EQ. -1) THEN
!!$        d1 = SQRT(DOT_PRODUCT(v,v))
!!$        g = MATMUL(TRANSPOSE(fjac),fjac) + lam*dtd
!!$        CALL DPOTRF('U', n, g, n, info)
!!$        grad = MATMUL(fvec, fjac)
!!$        CALL DPOTRS('U', n, 1, g, n, grad, n, info)
!!$        d2 = SQRT(DOT_PRODUCT(grad, grad) )
!!$        IF( 10.0d+0*ABS(d2-d1) .LT. d2 ) lam = lam - 0.1*d2*(lamold - lam)/(d1 - d2)
!!$     END IF
  ENDIF
END SUBROUTINE Updatelam_Umrigar


!! Trust region update methods

SUBROUTINE Updatedelta_factor(delta, accepted, factoraccept, factorreject)
  !! Update lam based on accepted/rejected step
  IMPLICIT NONE
  INTEGER accepted
  REAL (KIND=8) delta, factoraccept, factorreject

  IF(accepted.GE.0) THEN
     delta = delta * factoraccept
  ELSE
     delta = delta / factorreject
  ENDIF
END SUBROUTINE Updatedelta_factor

SUBROUTINE Updatedelta_more(delta, lam, n, v, dtd, rho, C, Cnew, dirder, actred, av, avmax)
  IMPLICIT NONE
  INTEGER n
  DOUBLE PRECISION delta, lam, v(n), dtd(n,n), rho, C, Cnew, dirder, actred, av, avmax
  DOUBLE PRECISION pnorm, temp
  pnorm = SQRT(DOT_PRODUCT(v,MATMUL(dtd, v)))
  IF (rho .GT. 0.25d+0) THEN
     IF (lam .GT. 0.0d+0 .AND. rho .LT. 0.75d+0) THEN
        temp = 1.0d+0
     ELSE
        temp = 2.0d+0*pnorm/delta
     END IF
  ELSE
     IF ( actred .GE. 0.0d+0) THEN
        temp = 0.5d+0
     ELSE
        temp = 0.5d+0*dirder/(dirder + 0.5d+0*actred)
     END IF
     IF ( 0.01*Cnew .GE. C .OR. temp .LT. 0.1d+0) temp = 0.1d+0
  END IF
  !! We need to make sure that if acceleration is too big, we decrease the step size
  IF (av .GT. avmax) THEN
     temp = MIN(temp,MAX(avmax/av,0.1d+0))
  END IF

  delta = temp*MIN(delta,10.0d+0*pnorm)
  lam = lam/temp
END SUBROUTINE Updatedelta_more


! -*- f90 -*-
!****************************************
! Routine for rank-deficient jacobian update

SUBROUTINE UPDATEJAC(m,n,fjac, fvec, fvec_new, acc, v, a)
  IMPLICIT NONE
  INTEGER m, n, i, j
  REAL (KIND=8) fjac(m,n), fvec(m), fvec_new(m), acc(m)
  REAL (KIND=8) v(n), a(n), djac(m), v2(n), r1(m)

  r1 = fvec + 0.5*MATMUL(fjac,v) + 0.125d+0*acc
  djac = 2.0*(r1 - fvec - 0.5*MATMUL(fjac,v))/DOT_PRODUCT(v,v)
  DO i = 1,m
     DO j = 1,n
        fjac(i,j) = fjac(i,j) + djac(i)*0.5d+0*v(j)
     END DO
  END DO
  v2 = 0.5d+0*(v + a)
  djac = 0.5*(fvec_new - r1 - MATMUL(fjac,v2))/DOT_PRODUCT(v2,v2)
  DO i = 1,m
     DO j = 1,n
        fjac(i,j) = fjac(i,j) + djac(i)*v2(j)
     END DO
  END DO
END SUBROUTINE UPDATEJAC

 



endmodule
! -*- f90 -*-            
! ****************************************
! Routine for calculating finite-difference second directional derivative

      subroutine destsv(n,r,ldr,svmin,z)
      integer ldr, n
      double precision svmin
      double precision r(ldr,n), z(n)
!     **********
!
!     Subroutine destsv
!
!     Given an n by n upper triangular matrix R, this subroutine
!     estimates the smallest singular value and the associated
!     singular vector of R.
!
!     In the algorithm a vector e is selected so that the solution
!     y to the system R'*y = e is large. The choice of sign for the
!     components of e cause maximal local growth in the components
!     of y as the forward substitution proceeds. The vector z is
!     the solution of the system R*z = y, and the estimate svmin
!     is norm(y)/norm(z) in the Euclidean norm.
!
!     The subroutine statement is
!
!       subroutine estsv(n,r,ldr,svmin,z)
!
!     where
!
!       n is an integer variable.
!         On entry n is the order of R.
!         On exit n is unchanged.
!
!       r is a double precision array of dimension (ldr,n)
!         On entry the full upper triangle must contain the full
!            upper triangle of the matrix R.
!         On exit r is unchanged.
!
!       ldr is an integer variable.
!         On entry ldr is the leading dimension of r.
!         On exit ldr is unchanged.
!
!       svmin is a double precision variable.
!         On entry svmin need not be specified.
!         On exit svmin contains an estimate for the smallest
!            singular value of R.
!
!       z is a double precision array of dimension n.
!         On entry z need not be specified.
!         On exit z contains a singular vector associated with the
!            estimate svmin such that norm(R*z) = svmin and
!            norm(z) = 1 in the Euclidean norm.
!
!     Subprograms called
!
!       Level 1 BLAS ... dasum, daxpy, dnrm2, dscal
!
!     MINPACK-2 Project. October 1993.
!     Argonne National Laboratory
!     Brett M. Averick and Jorge J. More'.
!
!     **********
      double precision one, p01, zero
      parameter (zero=0.d0,p01=1.0d-2,one=1.0d0)

      integer i, j
      double precision e, s, sm, temp, w, wm, ynorm, znorm

      double precision dasum, dnrm2
      external dasum, daxpy, dnrm2, dscal

      do 10 i = 1, n
         z(i) = zero
   10 continue

!     This choice of e makes the algorithm scale invariant.

      e = abs(r(1,1))
      if (e .eq. zero) then
         svmin = zero
         z(1) = one

         return

      end if

!     Solve R'*y = e.

      do 30 i = 1, n
         e = sign(e,-z(i))

!        Scale y. The factor of 0.01 reduces the number of scalings.

         if (abs(e-z(i)) .gt. abs(r(i,i))) then
            temp = min(p01,abs(r(i,i))/abs(e-z(i)))
            call dscal(n,temp,z,1)
            e = temp*e
         end if

!        Determine the two possible choices of y(i).

         if (r(i,i) .eq. zero) then
            w = one
            wm = one
         else
            w = (e-z(i))/r(i,i)
            wm = -(e+z(i))/r(i,i)
         end if

!        Choose y(i) based on the predicted value of y(j) for j > i.

         s = abs(e-z(i))
         sm = abs(e+z(i))
         do 20 j = i + 1, n
            sm = sm + abs(z(j)+wm*r(i,j))
   20    continue
         if (i .lt. n) then
            call daxpy(n-i,w,r(i,i+1),ldr,z(i+1),1)
            s = s + dasum(n-i,z(i+1),1)
         end if
         if (s .lt. sm) then
            temp = wm - w
            w = wm
            if (i .lt. n) call daxpy(n-i,temp,r(i,i+1),ldr,z(i+1),1)
         end if
         z(i) = w
   30 continue
      ynorm = dnrm2(n,z,1)

!     Solve R*z = y.

      do 40 j = n, 1, -1

!        Scale z.

         if (abs(z(j)) .gt. abs(r(j,j))) then
            temp = min(p01,abs(r(j,j))/abs(z(j)))
            call dscal(n,temp,z,1)
            ynorm = temp*ynorm
         end if
         if (r(j,j) .eq. zero) then
            z(j) = one
         else
            z(j) = z(j)/r(j,j)
         end if
         temp = -z(j)
         call daxpy(j-1,temp,r(1,j),1,z,1)
   40 continue

!     Compute svmin and normalize z.

      znorm = one/dnrm2(n,z,1)
      svmin = ynorm*znorm
      call dscal(n,znorm,z,1)

      end
      subroutine dgqt(n,a,lda,b,delta,rtol,atol,itmax,par,f,x,info,iter, z,wa1,wa2)
      integer n, lda, itmax, info
      double precision delta, rtol, atol, par, f
      double precision a(lda,n), b(n), x(n), z(n), wa1(n), wa2(n)
!     ***********
!
!     Subroutine dgqt
!
!     Given an n by n symmetric matrix A, an n-vector b, and a
!     positive number delta, this subroutine determines a vector
!     x which approximately minimizes the quadratic function
!
!           f(x) = (1/2)*x'*A*x + b'*x
!
!     subject to the Euclidean norm constraint
!
!           norm(x) <= delta.
!
!     This subroutine computes an approximation x and a Lagrange
!     multiplier par such that either par is zero and
!
!            norm(x) <= (1+rtol)*delta,
!
!     or par is positive and
!
!            abs(norm(x) - delta) <= rtol*delta.
!
!     If xsol is the solution to the problem, the approximation x
!     satisfies
!
!            f(x) <= ((1 - rtol)**2)*f(xsol)
!
!     The subroutine statement is
!
!       subroutine dgqt(n,a,lda,b,delta,rtol,atol,itmax,
!                        par,f,x,info,z,wa1,wa2)
!
!     where
!
!       n is an integer variable.
!         On entry n is the order of A.
!         On exit n is unchanged.
!
!       a is a double precision array of dimension (lda,n).
!         On entry the full upper triangle of a must contain the
!            full upper triangle of the symmetric matrix A.
!         On exit the array contains the matrix A.
!
!       lda is an integer variable.
!         On entry lda is the leading dimension of the array a.
!         On exit lda is unchanged.
!
!       b is an double precision array of dimension n.
!         On entry b specifies the linear term in the quadratic.
!         On exit b is unchanged.
!
!       delta is a double precision variable.
!         On entry delta is a bound on the Euclidean norm of x.
!         On exit delta is unchanged.
!
!       rtol is a double precision variable.
!         On entry rtol is the relative accuracy desired in the
!            solution. Convergence occurs if
!
!              f(x) <= ((1 - rtol)**2)*f(xsol)
!
!         On exit rtol is unchanged.
!
!       atol is a double precision variable.
!         On entry atol is the absolute accuracy desired in the
!            solution. Convergence occurs when
!
!              norm(x) <= (1 + rtol)*delta
!
!              max(-f(x),-f(xsol)) <= atol
!
!         On exit atol is unchanged.
!
!       itmax is an integer variable.
!         On entry itmax specifies the maximum number of iterations.
!         On exit itmax is unchanged.
!
!       par is a double precision variable.
!         On entry par is an initial estimate of the Lagrange
!            multiplier for the constraint norm(x) <= delta.
!         On exit par contains the final estimate of the multiplier.
!
!       f is a double precision variable.
!         On entry f need not be specified.
!         On exit f is set to f(x) at the output x.
!
!       x is a double precision array of dimension n.
!         On entry x need not be specified.
!         On exit x is set to the final estimate of the solution.
!
!       info is an integer variable.
!         On entry info need not be specified.
!         On exit info is set as follows:
!
!            info = 1  The function value f(x) has the relative
!                      accuracy specified by rtol.
!
!            info = 2  The function value f(x) has the absolute
!                      accuracy specified by atol.
!
!            info = 3  Rounding errors prevent further progress.
!                      On exit x is the best available approximation.
!
!            info = 4  Failure to converge after itmax iterations.
!                      On exit x is the best available approximation.
!
!       z is a double precision work array of dimension n.
!
!       wa1 is a double precision work array of dimension n.
!
!       wa2 is a double precision work array of dimension n.
!
!     Subprograms called
!
!       MINPACK-2  ......  destsv
!
!       LAPACK  .........  dpotrf
!
!       Level 1 BLAS  ...  dasum, daxpy, dcopy, ddot, dnrm2, dscal
!
!       Level 2 BLAS  ...  dtrmv, dtrsv
!
!     MINPACK-2 Project. July 1994.
!     Argonne National Laboratory and University of Minnesota.
!     Brett M. Averick, Richard Carter, and Jorge J. More'
!
!     ***********
      double precision one, p001, p5, zero
      parameter (zero=0.0d0,p001=1.0d-3,p5=0.5d0,one=1.0d0)

      logical rednc
      integer indef, iter, j
      double precision alpha, anorm, bnorm, parc, parf, parl, pars, paru, prod, rxnorm, rznorm, temp, xnorm

      double precision dasum, ddot, dnrm2
      external daxpy, dcopy, ddot, destsv, dnrm2, dpotrf, dscal, dtrsv

!     Initialization.

      parf = zero
      xnorm = zero
      rxnorm = zero
      rednc = .false.
      do 10 j = 1, n
         x(j) = zero
         z(j) = zero
   10 continue

!     Copy the diagonal and save A in its lower triangle.

      call dcopy(n,a,lda+1,wa1,1)
      do 20 j = 1, n - 1
         call dcopy(n-j,a(j,j+1),lda,a(j+1,j),1)
   20 continue

!     Calculate the l1-norm of A, the Gershgorin row sums,
!     and the l2-norm of b.

      anorm = zero
      do 30 j = 1, n
         wa2(j) = dasum(n,a(1,j),1)
         anorm = max(anorm,wa2(j))
   30 continue
      do 40 j = 1, n
         wa2(j) = wa2(j) - abs(wa1(j))
   40 continue
      bnorm = dnrm2(n,b,1)

!     Calculate a lower bound, pars, for the domain of the problem.
!     Also calculate an upper bound, paru, and a lower bound, parl,
!     for the Lagrange multiplier.

      pars = -anorm
      parl = -anorm
      paru = -anorm
      do 50 j = 1, n
         pars = max(pars,-wa1(j))
         parl = max(parl,wa1(j)+wa2(j))
         paru = max(paru,-wa1(j)+wa2(j))
   50 continue
      parl = max(zero,bnorm/delta-parl,pars)
      paru = max(zero,bnorm/delta+paru)

!     If the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.

      par = max(par,parl)
      par = min(par,paru)

!     Special case: parl = paru.

      paru = max(paru,(one+rtol)*parl)

!     Beginning of an iteration.

      info = 0
      do 90 iter = 1, itmax

!        Safeguard par.

         if (par .le. pars .and. paru .gt. zero) par = max(p001, sqrt(parl/paru))*paru
!        write (2,1000) parl,paru,pars,par
!1000    format(4d12.3)

!        Copy the lower triangle of A into its upper triangle and
!        compute A + par*I.

         do 60 j = 1, n - 1
            call dcopy(n-j,a(j+1,j),1,a(j,j+1),lda)
   60    continue
         do 70 j = 1, n
            a(j,j) = wa1(j) + par
   70    continue

!        Attempt the  Cholesky factorization of A without referencing
!        the lower triangular part.

         call dpotrf('U',n,a,lda,indef)

!        Case 1: A + par*I is positive definite.

         if (indef .eq. 0) then

!           Compute an approximate solution x and save the
!           last value of par with A + par*I positive definite.

            parf = par
            call dcopy(n,b,1,wa2,1)
            call dtrsv('U','T','N',n,a,lda,wa2,1)
            rxnorm = dnrm2(n,wa2,1)
            call dtrsv('U','N','N',n,a,lda,wa2,1)
            call dcopy(n,wa2,1,x,1)
            call dscal(n,-one,x,1)
            xnorm = dnrm2(n,x,1)

!           Test for convergence.

            if (abs(xnorm-delta) .le. rtol*delta .or. (par .eq. zero .and. xnorm .le. (one+rtol)*delta)) info = 1

!           Compute a direction of negative curvature and use this
!           information to improve pars.

            call destsv(n,a,lda,rznorm,z)
            pars = max(pars,par-rznorm**2)

!           Compute a negative curvature solution of the form
!           x + alpha*z where norm(x+alpha*z) = delta.

            rednc = .false.
            if (xnorm .lt. delta) then

!              Compute alpha

               prod = ddot(n,z,1,x,1)/delta
               temp = (delta-xnorm)*((delta+xnorm)/delta)
               alpha = temp/(abs(prod)+sqrt(prod**2+temp/delta))
               alpha = sign(alpha,prod)

!              Test to decide if the negative curvature step
!              produces a larger reduction than with z = 0.

               rznorm = abs(alpha)*rznorm
               if ((rznorm/delta)**2+par*(xnorm/delta)**2 .le. par) rednc = .true.

!              Test for convergence.

               if (p5*(rznorm/delta)**2 .le. rtol*(one-p5*rtol)*(par+(rxnorm/delta)**2)) then 
                  info = 1
               else if (p5*(par+(rxnorm/delta)**2) .le. (atol/delta)/delta .and. info .eq. 0) then
                  info = 2
               else if (xnorm .eq. zero) then
                  info = 1
               end if
            end if

!           Compute the Newton correction parc to par.

            if (xnorm .eq. zero) then
               parc = -par
            else
               call dcopy(n,x,1,wa2,1)
               temp = one/xnorm
               call dscal(n,temp,wa2,1)
               call dtrsv('U','T','N',n,a,lda,wa2,1)
               temp = dnrm2(n,wa2,1)
               parc = (((xnorm-delta)/delta)/temp)/temp
            end if

!           Update parl or paru.

            if (xnorm .gt. delta) parl = max(parl,par)
            if (xnorm .lt. delta) paru = min(paru,par)
         else

!           Case 2: A + par*I is not positive definite.

!           Use the rank information from the Cholesky
!           decomposition to update par.

            if (indef .gt. 1) then

!              Restore column indef to A + par*I.

               call dcopy(indef-1,a(indef,1),lda,a(1,indef),1)
               a(indef,indef) = wa1(indef) + par

!              Compute parc.

               call dcopy(indef-1,a(1,indef),1,wa2,1)
               call dtrsv('U','T','N',indef-1,a,lda,wa2,1)
               a(indef,indef) = a(indef,indef) - dnrm2(indef-1,wa2,1)**2
               call dtrsv('U','N','N',indef-1,a,lda,wa2,1)
            end if
            wa2(indef) = -one
            temp = dnrm2(indef,wa2,1)
            parc = -(a(indef,indef)/temp)/temp
            pars = max(pars,par,par+parc)

!           If necessary, increase paru slightly.
!           This is needed because in some exceptional situations
!           paru is the optimal value of par.

            paru = max(paru,(one+rtol)*pars)
         end if

!        Use pars to update parl.

         parl = max(parl,pars)

!        Test for termination.

         if (info .eq. 0) then
            if (iter .eq. itmax) info = 4
            if (paru .le. (one+p5*rtol)*pars) info = 3
            if (paru .eq. zero) info = 2
         end if

!        If exiting, store the best approximation and restore
!        the upper triangle of A.

         if (info .ne. 0) then

!           Compute the best current estimates for x and f.

            par = parf
            f = -p5*(rxnorm**2+par*xnorm**2)
            if (rednc) then
               f = -p5*((rxnorm**2+par*delta**2)-rznorm**2)
               call daxpy(n,alpha,z,1,x,1)
            end if

!           Restore the upper triangle of A.

            do 80 j = 1, n - 1
               call dcopy(n-j,a(j+1,j),1,a(j,j+1),lda)
   80       continue
            call dcopy(n,wa1,1,a,lda+1)

            return

         end if

!        Compute an improved estimate for par.

         par = max(parl,par+parc)

!        End of an iteration.

   90 continue

      end
      double precision function dpmpar(i)
      integer i
!     **********
!
!     Function dpmpar
!
!     This function provides double precision machine parameters
!     when the appropriate set of data statements is activated (by
!     removing the c from column 1) and all other data statements are
!     rendered inactive. Most of the parameter values were obtained
!     from the corresponding Bell Laboratories Port Library function.
!
!     The function statement is
!
!       double precision function dpmpar(i)
!
!     where
!
!       i is an integer input variable set to 1, 2, or 3 which
!         selects the desired machine parameter. If the machine has
!         t base b digits and its smallest and largest exponents are
!         emin and emax, respectively, then these parameters are
!
!         dpmpar(1) = b**(1 - t), the machine precision,
!
!         dpmpar(2) = b**(emin - 1), the smallest magnitude,
!
!         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
!
!     Argonne National Laboratory. MINPACK Project. November 1996.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
!
!     **********
      integer mcheps(4)
      integer minmag(4)
      integer maxmag(4)
      double precision dmach(3)
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
!
!     Machine constants for the IBM 360/370 series,
!     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
!     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
!
!     data mcheps(1),mcheps(2) / z34100000, z00000000 /
!     data minmag(1),minmag(2) / z00100000, z00000000 /
!     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
!
!     Machine constants for the Honeywell 600/6000 series.
!
!     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
!     data minmag(1),minmag(2) / o402400000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
!
!     Machine constants for the CDC 6000/7000 series.
!
!     data mcheps(1) / 15614000000000000000b /
!     data mcheps(2) / 15010000000000000000b /
!
!     data minmag(1) / 00604000000000000000b /
!     data minmag(2) / 00000000000000000000b /
!
!     data maxmag(1) / 37767777777777777777b /
!     data maxmag(2) / 37167777777777777777b /
!
!     Machine constants for the PDP-10 (KA processor).
!
!     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
!     data minmag(1),minmag(2) / "033400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
!
!     Machine constants for the PDP-10 (KI processor).
!
!     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
!     data minmag(1),minmag(2) / "000400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
!
!     Machine constants for the PDP-11. 
!
!     data mcheps(1),mcheps(2) /   9472,      0 /
!     data mcheps(3),mcheps(4) /      0,      0 /
!
!     data minmag(1),minmag(2) /    128,      0 /
!     data minmag(3),minmag(4) /      0,      0 /
!
!     data maxmag(1),maxmag(2) /  32767,     -1 /
!     data maxmag(3),maxmag(4) /     -1,     -1 /
!
!     Machine constants for the Burroughs 6700/7700 systems.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o7770000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o7777777777777777 /
!
!     Machine constants for the Burroughs 5700 system.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o0000000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o0007777777777777 /
!
!     Machine constants for the Burroughs 1700 system.
!
!     data mcheps(1) / zcc6800000 /
!     data mcheps(2) / z000000000 /
!
!     data minmag(1) / zc00800000 /
!     data minmag(2) / z000000000 /
!
!     data maxmag(1) / zdffffffff /
!     data maxmag(2) / zfffffffff /
!
!     Machine constants for the Univac 1100 series.
!
!     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
!     data minmag(1),minmag(2) / o000040000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
!
!     Machine constants for the Data General Eclipse S/200.
!
!     Note - it may be appropriate to include the following card -
!     static dmach(3)
!
!     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
!     data mcheps/32020k,3*0/
!
!     Machine constants for the Harris 220.
!
!     data mcheps(1),mcheps(2) / '20000000, '00000334 /
!     data minmag(1),minmag(2) / '20000000, '00000201 /
!     data maxmag(1),maxmag(2) / '37777777, '37777577 /
!
!     Machine constants for the Cray-1.
!
!     data mcheps(1) / 0376424000000000000000b /
!     data mcheps(2) / 0000000000000000000000b /
!
!     data minmag(1) / 0200034000000000000000b /
!     data minmag(2) / 0000000000000000000000b /
!
!     data maxmag(1) / 0577777777777777777777b /
!     data maxmag(2) / 0000007777777777777776b /
!
!     Machine constants for the Prime 400.
!
!     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
!     data minmag(1),minmag(2) / :10000000000, :00000100000 /
!     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
!
!     Machine constants for the VAX-11.
!
!     data mcheps(1),mcheps(2) /   9472,  0 /
!     data minmag(1),minmag(2) /    128,  0 /
!     data maxmag(1),maxmag(2) / -32769, -1 /
!
!     Machine constants for IEEE machines.
!
      data dmach(1) /2.22044604926d-16/
      data dmach(2) /2.22507385852d-308/
      data dmach(3) /1.79769313485d+308/
!
      dpmpar = dmach(i)
      return
!
!     Last card of function dpmpar.
!
      end
