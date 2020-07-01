#include "alias.inc"
subroutine gen_algo( get_eig, NN_TABLE, kpoint, nkpoint, E_DFT, neig, &
                     iband, nband, PWGHT, PINPT, PKAIA)
    use parameters
    use mpi_setup
    use time
    use pikaia_module, only: pikaia_class
    use print_io
    implicit none
    type(incar)        :: PINPT
    type(gainp)        :: PKAIA
    type(hopping)      :: NN_TABLE
    type(weight)       :: PWGHT
    type(pikaia_class) :: ga_pikaia
    integer*4             nkpoint, neig, iband, nband, info
    integer*4                      ne  , ib   , nb
    real*8                kpoint(3,nkpoint)
    real*8                E_DFT(neig*PINPT%ispin, nkpoint)
    external              get_eig
    integer*4             nparam                  
    integer*4             istat
    real*8                param(PINPT%nparam) ! parameters
    real*8                gofit               ! goodness of fit as a function of param gofit(param) = [0:1]
    real*8                xl(PINPT%nparam) ! lower bound for x
    real*8                xu(PINPT%nparam) ! upper bound for x
    logical               flag_header_write
    external              report_iter, get_gofit
    write(message,'(A)')' Start: fitting procedures with Genetic Algorithm '  ; write_msg
    write(message,'(A)')'        based on PIKAIA library.'  ; write_msg
   
    nparam = PINPT%nparam
    param  = PINPT%param(:)
    xl     = PINPT%param_const(3, :)
    xu     = PINPT%param_const(2, :)

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
    call ga_pikaia%solve(param,gofit,istat,PINPT,NN_TABLE,E_DFT,PWGHT)

    PINPT%param(:) = param
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

    !the first time it is called, also write a header:
!   if (.not. header_written) then
!     do i=1,nparam
!       write(xheader(i),'(I10)') i
!       xheader(i) = 'X'//trim(adjustl(xheader(i)))
!       xheader(i) = repeat(' ',10-len_trim(xheader(i)))//xheader(i)
!     end do
!     write(message,'(A5,1X,*(A10,1X))') 'ITER',xheader,'F'  ; write_msg
!     header_written = .true.
!   end if

!   write(message,'(I5,1X,*(F10.6,1X))') iter,param,gofit  ; write_msg
    write(message,'(A)')' '  ; write_msg
    write(message,'(A,I4,A,F16.6)')'   ITER=',iter,',(EDFT-ETBA)*WEIGHT = ', -gofit  ; write_msg

    return
endsubroutine

subroutine get_gofit(ga_pikaia, param, gofit,kpoint,nkpoint,neig, iband, nband, NN_TABLE, E_DFT, PWGHT, PINPT, flag_lmdif)
  use parameters, only : hopping, weight, incar
  use pikaia_module, only: pikaia_class
  implicit none
  class(pikaia_class),intent(inout) :: ga_pikaia
  type(incar)                       :: PINPT
  type(weight)                      :: PWGHT
  type(hopping)                     :: NN_TABLE
  integer*4,intent(in)              :: nkpoint, neig, iband, nband
  real*8,intent(in)                 :: E_DFT(neig*PINPT%ispin,nkpoint)
  real*8,intent(out)                :: gofit   ! fitness value
  real*8,intent(inout)                 :: param(PINPT%nparam)    !unscaled x vector: [xu,xl]
  real*8,intent(in)                 :: kpoint(3,nkpoint)
  logical,intent(in)                :: flag_lmdif
   
  logical                              flag_get_orbital, flag_order
  complex*16                           V(neig*PINPT%ispin,nband*PINPT%nspin,nkpoint)
  real*8                               E_TBA(nband*PINPT%nspin,nkpoint) 
  real*8                               fnorm, fvec(nkpoint)
  real*8, external                  :: enorm
  external                             get_eig

  ! store param which is generated by PIKAIA to PINPT%param before calling get_eig subroutine.
  PINPT%param = param

  flag_order       = PINPT%flag_get_band_order
  flag_get_orbital = ((.not. PWGHT%flag_weight_default_orb) .or. flag_order)

  if(flag_lmdif) then
    call leasqr_lm ( get_eig, NN_TABLE, kpoint, nkpoint, E_DFT, neig, iband, nband, PWGHT, PINPT)
  endif

  call get_eig(NN_TABLE, kpoint, nkpoint, PINPT, E_TBA, V, neig, iband, nband, flag_get_orbital, .false., .false., PINPT%flag_phase, flag_order)

  if(.not. flag_get_orbital) then
    call get_fvec(fvec, E_TBA, E_DFT, 0, neig, iband, nband, PINPT, nkpoint, PWGHT)
  else
    call get_fvec(fvec, E_TBA, E_DFT, V, neig, iband, nband, PINPT, nkpoint, PWGHT)
  endif

  fnorm = enorm ( nkpoint, fvec )
  gofit = -fnorm
! gofit = -exp(-1d0/(fnorm**2d0)) + 1d0

  return
endsubroutine

