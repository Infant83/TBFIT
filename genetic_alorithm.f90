#include "alias.inc"
subroutine gen_algo( get_eig, NN_TABLE, E_DFT, PWGHT, PINPT, PPRAM, PKPTS, PGEOM, PKAIA)
    use parameters
    use mpi_setup
    use time
    use pikaia_module, only: pikaia_class
    use print_io
    implicit none
    type(incar)        :: PINPT
    type(params)       :: PPRAM
    type(kpoints)      :: PKPTS
    type(poscar)       :: PGEOM 
    type(gainp)        :: PKAIA
    type(hopping)      :: NN_TABLE
    type(weight)       :: PWGHT
    type(pikaia_class) :: ga_pikaia
    integer*4             nkpoint, neig, iband, nband, info
    integer*4                      ne  , ib   , nb
    real*8                kpoint(3,PKPTS%nkpoint)
    real*8                E_DFT(PGEOM%neig*PINPT%ispin, PKPTS%nkpoint)
    external              get_eig
    integer*4             nparam                  
    integer*4             istat
    real*8                param(PPRAM%nparam) ! parameters
    real*8                gofit               ! goodness of fit as a function of param gofit(param) = [0:1]
    real*8                xl(PPRAM%nparam) ! lower bound for x
    real*8                xu(PPRAM%nparam) ! upper bound for x
    logical               flag_header_write
    external              report_iter, get_gofit
    write(message,'(A)')' Start: fitting procedures with Genetic Algorithm '  ; write_msg
    write(message,'(A)')'        based on PIKAIA library.'  ; write_msg
 
    nkpoint= PKPTS%nkpoint
    kpoint = PKPTS%kpoint
    neig   = PGEOM%neig
    nband  = PGEOM%nband
    iband  = PGEOM%init_erange  
    nparam = PPRAM%nparam
    param  = PPRAM%param(:)
    xl     = PPRAM%param_const(3, :)
    xu     = PPRAM%param_const(2, :)

    ! Initialize the class:
    call ga_pikaia%init(nparam, xl, xu, get_gofit, istat, nkpoint, kpoint, &
                        ne                 = neig         ,&
                        ib                 = iband        ,&
                        nb                 = nband        ,&
                        iter_f             = report_iter  ,&
                        np                 = PKAIA%npop   ,&
                        ngen               = PKAIA%mgen   ,&
                        nd                 = PKAIA%ngene  ,&
                        pcross             = PKAIA%pcross ,&
                        pmutmn             = PKAIA%pmutmn ,&
                        pmutmx             = PKAIA%pmutmx ,&
                        pmut               = PKAIA%pmut   ,&
                        imut               = PKAIA%imut   ,&
                        fdif               = PKAIA%fdif   ,&
                        irep               = PKAIA%irep   ,&
                        ielite             = PKAIA%ielite ,&
                        ivrb               = PKAIA%ivrb   ,&
                        convergence_tol    = PKAIA%convtol,&
                        convergence_window = PKAIA%convwin,&
                        initial_guess_frac = PKAIA%iguessf,&
                        iseed              = PKAIA%iseed)


    ! Run fitting by calling pikaia: 
    call ga_pikaia%solve(param,gofit,istat,PINPT,PPRAM,PKPTS,PGEOM,NN_TABLE,E_DFT,PWGHT)

    PPRAM%param(:) = param
    write(message,'(A)')'   End: fitting procedures.'    ; write_msg

    return
endsubroutine

subroutine report_iter(me,iter,param,gofit,nparam,header_written)
    use pikaia_module, only : pikaia_class
    use mpi_setup
    use print_io
    implicit none
    class(pikaia_class),intent(inout) :: me
    integer*4,intent(in)              :: iter, nparam
    real*8,dimension(:),intent(in)    :: param
    real*8,intent(in)                 :: gofit
    character*10                      :: xheader(nparam)
    integer                           :: i
    logical                              header_written

    write(message,'(A)')' '  ; write_msg
    write(message,'(A,I4,A,F16.6)')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ', -gofit  ; write_msg

    return
endsubroutine

subroutine get_gofit(ga_pikaia, param, gofit,kpoint,nkpoint,neig, iband, nband, NN_TABLE, E_DFT, PWGHT, PINPT, PPRAM, PKPTS, PGEOM, flag_lmdif)
  use parameters, only : hopping, weight, incar, params, kpoints, poscar
  use pikaia_module, only: pikaia_class
  implicit none
  class(pikaia_class),intent(inout) :: ga_pikaia
  type(incar)                       :: PINPT
  type(params)                      :: PPRAM
  type(kpoints)                     :: PKPTS
  type(poscar)                      :: PGEOM
  type(weight)                      :: PWGHT
  type(hopping)                     :: NN_TABLE
  integer*4,intent(in)              :: nkpoint, neig, iband, nband
  real*8,intent(in)                 :: E_DFT(neig*PINPT%ispin,nkpoint)
  real*8,intent(out)                :: gofit   ! fitness value
  real*8,intent(inout)                 :: param(PPRAM%nparam)    !unscaled x vector: [xu,xl]
  real*8,intent(in)                 :: kpoint(3,nkpoint)
  logical,intent(in)                :: flag_lmdif
   
  logical                              flag_get_orbital, flag_order
  complex*16                           V(neig*PINPT%ispin,nband*PINPT%nspin,nkpoint)
  complex*16                           SV(neig*PINPT%ispin,nband*PINPT%nspin,nkpoint)
  real*8                               E_TBA(nband*PINPT%nspin,nkpoint) 
  real*8                               fvec(nkpoint)
  real*8                               fnorm, fnorm_(2)
  real*8, external                  :: enorm
  external                             get_eig

  ! store param which is generated by PIKAIA to PINPT%param before calling get_eig subroutine.
  PPRAM%param = param

  flag_order       = PINPT%flag_get_band_order
  flag_get_orbital = ( PWGHT%flag_weight_orb .or. flag_order)

  if(flag_lmdif) then
    call leasqr_lm ( get_eig, NN_TABLE, E_DFT, PWGHT, PINPT, PPRAM, PKPTS, PGEOM, fnorm_)
  endif

  call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, PPRAM, E_TBA, V, SV, neig, iband, nband, &
               flag_get_orbital, .false., .false., PINPT%flag_phase)

  if(.not. flag_get_orbital) then
    call get_fvec(fvec, E_TBA, E_DFT, (0d0,0d0), neig, iband, nband, PINPT, nkpoint, PWGHT)
  else
    call get_fvec(fvec, E_TBA, E_DFT, V, neig, iband, nband, PINPT, nkpoint, PWGHT)
  endif

  fnorm = enorm ( nkpoint, fvec )
  gofit = -fnorm

  return
endsubroutine

subroutine get_fvec (fvec, E_TBA, E_DFT, V, neig, iband, nband, PINPT, nkpoint, PWGHT)
  use parameters, only: weight, incar
  implicit none
  type (weight)  :: PWGHT
  type (incar )  :: PINPT
  integer*4  ik, neig, nkpoint
  integer*4  is, ie, ie_, iband, nband
  real*8     E_TBA(nband*PINPT%nspin,nkpoint)
  real*8     E_DFT(neig*PINPT%ispin,nkpoint), dE(nband*PINPT%nspin), fvec(nkpoint)
  real*8     orbital_weight_penalty(nband*PINPT%nspin)
  complex*16, intent(in) :: V(neig*PINPT%ispin,nband*PINPT%nspin,nkpoint)
  real*8     enorm
  external   enorm

  if(.not. PWGHT%flag_weight_orb) then

    do ik=1,nkpoint
      do is = 1, PINPT%nspin
        do ie = 1, nband
          ie_ = ie + iband - 1 + (is-1) * neig
          dE(ie+(is-1)*nband) =    (E_TBA(ie+(is-1)*nband,ik) - E_DFT(ie_,ik))    * PWGHT%WT(ie_,ik)
        enddo
      enddo
      fvec(ik) = enorm ( nband*PINPT%nspin, dE )
    enddo

  elseif(PWGHT%flag_weight_orb) then

    do ik=1,nkpoint
      do is = 1, PINPT%nspin
        do ie = 1, nband
          ie_ = ie + iband - 1 + (is-1) * neig
          orbital_weight_penalty(ie+(is-1)*nband) = sum( PWGHT%PENALTY_ORB(:,ie_,ik)*abs(V(:,ie+(is-1)*nband,ik)) )
          dE(ie+(is-1)*nband) = (E_TBA(ie+(is-1)*nband,ik) - E_DFT(ie_,ik)) * PWGHT%WT(ie_,ik) + orbital_weight_penalty(ie_)
        enddo
      enddo
      fvec(ik) = enorm ( nband*PINPT%nspin, dE )
    enddo

  endif
return
endsubroutine

