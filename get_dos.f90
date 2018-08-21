subroutine get_dos(NN_TABLE, PINPT, PINPT_DOS, PGEOM, PKPTS)
   use parameters  
   use mpi_setup
   ! some part of this routine is adopted from dos.f90 file of WannierTools originally written by QuanSheng Wu 
   !  (Ref. https://doi.org/10.1016/j.cpc.2017.09.)
   implicit none
   type(hopping) :: NN_TABLE
   type(dos    ) :: PINPT_DOS
   type(incar  ) :: PINPT
   type(poscar ) :: PGEOM
   type(kpoints) :: PKPTS
   integer*4        i,ie,nkpoint,nparam,neig, nediv, ispin
   integer*4        ik,nk1,nk2,nk3
   integer*4        iband, fband, nband
   integer*4        is, ispin_print
   integer*4        mpierr
   real*8           kshift(3)
   real*8           e_range(PINPT_DOS%dos_nediv)
   real*8           emax, emin
   real*8           sigma, x, g_smear
   real*8,     allocatable :: kpoint(:,:), kpoint_reci(:,:)
   real*8,     allocatable :: param(:), param_const(:,:)
   real*8,     allocatable :: dos_up(:), dos_dn(:)
   real*8,     allocatable :: E(:,:), E_(:,:)
   complex*16, allocatable :: V(:,:,:), V_(:,:,:)
   real*8           a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
   real*8           b2xb3(3),bzvol,dkv
   real*8           fgauss
   external         fgauss  
   character*40     fname_header
   
   neig    = PGEOM%neig
   ispin   = PINPT%ispin
   nk1     = PINPT_DOS%dos_kgrid(1)  
   nk2     = PINPT_DOS%dos_kgrid(2)  
   nk3     = PINPT_DOS%dos_kgrid(3)  
   nkpoint = nk1 * nk2 * nk3
   nparam  = PINPT%nparam
   kshift  = PINPT_DOS%dos_kshift
   nediv   = PINPT_DOS%dos_nediv
   emax    = PINPT_DOS%dos_emax
   emin    = PINPT_DOS%dos_emin
   e_range = emin + (/0:nediv-1/) * (emax - emin)/dble(nediv - 1)
   g_smear = PINPT_DOS%dos_smearing
   iband   = PINPT_DOS%dos_iband
   fband   = PINPT_DOS%dos_fband 
   if(PINPT%flag_noncollinear) then
     if(fband .eq. 999999) fband = neig * 2
   else
     if(fband .eq. 999999) fband = neig
   endif
   nband   = PINPT_DOS%dos_fband - PINPT_DOS%dos_iband + 1

   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   call get_reci(b1,b2,b3, a1,a2,a3)
   call vcross(b2xb3,b2,b3)
   bzvol=dot_product(b1,b2xb3)
   dkv  = bzvol / dble(nkpoint)

   allocate( param(nparam) )
   allocate( param_const(5,nparam) )
   param   = PINPT%param
   param_const = PINPT%param_const
   allocate( E(neig*ispin,nkpoint) )
   allocate( V(neig*ispin,neig*ispin,nkpoint) )
   allocate( kpoint(3,nkpoint) )
   allocate( kpoint_reci(3,nkpoint) )
   allocate( PINPT_DOS%dos_kpoint(3,nkpoint) )
   allocate( PINPT_DOS%dos_erange(nediv) )
   if(PINPT%flag_collinear) then
     allocate( PINPT_DOS%dos_up(nediv), dos_up(nediv) )
     allocate( PINPT_DOS%dos_dn(nediv), dos_dn(nediv) )
     
     PINPT_DOS%dos_up = 0d0
     PINPT_DOS%dos_dn = 0d0
               dos_up = 0d0
               dos_dn = 0d0
   else
     allocate( PINPT_DOS%dos_up(nediv), dos_up(nediv) )
     PINPT_DOS%dos_up = 0d0
               dos_up = 0d0
   endif

   PINPT_DOS%dos_erange = e_range
   call get_kgrid(kpoint,kpoint_reci,nk1,nk2,nk3,kshift, PGEOM, PINPT_DOS%dos_flag_gamma)
   if(PINPT_DOS%dos_flag_print_kpoint) then 
     if(myid .eq. 0) call print_kpoint(kpoint_reci, nkpoint, PINPT_DOS%dos_kfilenm)
   endif

   call get_eig(NN_TABLE,kpoint,nkpoint,PINPT, E, V, neig, PINPT%flag_get_orbital, .true., .true.)
   sigma = g_smear

   call MPI_Barrier(mpi_comm_earth, mpierr)

kp:do ik = 1 + myid, nkpoint, nprocs
 eig:do ie = 1, nediv
 dosum:do i = iband, fband
         if(PINPT%flag_collinear) then
           x = e_range(ie) - E(i,ik)
           dos_up(ie) = dos_up(ie) + fgauss(sigma,x) * dkv
           x = e_range(ie) - E(i+neig,ik)
           dos_dn(ie) = dos_dn(ie) + fgauss(sigma,x) * dkv
         else
           x = e_range(ie) - E(i,ik)
           dos_up(ie) = dos_up(ie) + fgauss(sigma,x) * dkv
         endif
       enddo dosum
     enddo eig
   enddo kp

   if(flag_use_mpi)  call MPI_REDUCE(dos_up, PINPT_DOS%dos_up, nediv, MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
   if(flag_use_mpi .and. PINPT%flag_collinear) call MPI_REDUCE(dos_dn, PINPT_DOS%dos_dn, nediv, MPI_REAL8, MPI_SUM, 0, &
                                                               mpi_comm_earth, mpierr)
   if(myid .eq. 0)  call print_dos(PINPT_DOS, PINPT)

   if(PINPT_DOS%dos_flag_print_eigen .and. myid .eq. 0) then
     call get_ispin_print(PINPT%flag_collinear, ispin_print)
     allocate(E_(ispin_print,nkpoint))
     allocate(V_(neig*ispin,ispin_print,nkpoint))

     do i = 1, PINPT_DOS%dos_n_ensurf
       call get_ensurf_fname_header(PINPT_DOS%dos_ensurf(i), fname_header)
       do is = 1, ispin_print
         E_(is,:)   = E(  PINPT_DOS%dos_ensurf(i) + (is-1)*neig, :)
         V_(:,is,:) = V(:,PINPT_DOS%dos_ensurf(i) + (is-1)*neig, :)
       enddo
      call print_energy_ensurf(kpoint, nkpoint, ispin_print, E_, V_, PGEOM, PINPT, fname_header, PINPT_DOS%dos_kunit)
     enddo
   endif

   deallocate( param )
   deallocate( param_const)
   deallocate( E )
   deallocate( V )
   if(allocated(E_)) deallocate( E_)
   if(allocated(V_)) deallocate( V_)
   deallocate( kpoint )
   deallocate( kpoint_reci )

return
endsubroutine
function fgauss(sigma, x)
   use parameters, only : pi, pi2
   implicit none
   real*8   sigma,sigma2
   real*8   x,xx
   real*8   fgauss
   xx = x**2
   sigma2 = sigma**2

   fgauss= exp(-0.5d0*xx/sigma2)/(sigma*sqrt(pi2))

return
end function
subroutine get_ensurf_fname_format(i, format_string)
   implicit none
   integer*4    i
   character*40 format_string

   if(i .gt. 0 .and. i .lt. 10) then
     write(format_string, '( "(A,",A,",I1,A)" )')'"0000"'
   elseif(i .ge. 10 .and. i .lt. 100) then
     write(format_string, '( "(A,",A,",I2,A)" )')'"000"'
   elseif(i .ge. 100 .and. i .lt. 1000) then
     write(format_string, '( "(A,",A,",I3,A)" )')'"00"'
   elseif(i .ge. 1000 .and. i .lt. 10000) then
     write(format_string, '( "(A,",A,",I4,A)" )')'"0"'
   elseif(i .ge. 10000 .and. i .lt. 100000) then
     write(format_string, '( "(A,      I5,A)" )')
   endif

return
endsubroutine
subroutine get_ensurf_fname_header(i, fname_header)
   implicit none
   integer*4    i
   character*40 format_string
   character*80 fname_header

   call get_ensurf_fname_format(i, format_string)

   write(fname_header, format_string)'ENSURF.EIG.',i

return
endsubroutine

subroutine print_dos(PINPT_DOS, PINPT)
  use parameters, only : dos, pid_dos, incar
  implicit none
  type(dos)    :: PINPT_DOS
  type(incar)  :: PINPT
  integer*4       i, ie
  real*8          dos_data(PINPT_DOS%dos_nediv,PINPT%nspin)
  real*8          e_range(PINPT_DOS%dos_nediv)
  character*40    filenm
  logical         flag_collinear


  filenm = PINPT_DOS%dos_filenm
  e_range = PINPT_DOS%dos_erange

  if(.not.PINPT%flag_collinear) then
     dos_data(1:PINPT_DOS%dos_nediv, 1) = PINPT_DOS%dos_up(1:PINPT_DOS%dos_nediv)
  elseif(PINPT%flag_collinear) then
     dos_data(1:PINPT_DOS%dos_nediv, 1) = PINPT_DOS%dos_up(1:PINPT_DOS%dos_nediv)
     dos_data(1:PINPT_DOS%dos_nediv, 2) = PINPT_DOS%dos_dn(1:PINPT_DOS%dos_nediv)
  endif

  open(pid_dos, file=trim(filenm), status = 'unknown')

  write(pid_dos,'(A,I8,A,F16.8)')'# NDIV = ',PINPT_DOS%dos_nediv,' dE=',e_range(2)-e_range(1)
  if(.not.PINPT%flag_collinear) then
    write(pid_dos,'(A)')'#    energy (ev)           dos '
  elseif(PINPT%flag_collinear) then
    write(pid_dos,'(A)')'#    energy (ev)          dos-up           dos-dn'
  endif
  do ie = 1, PINPT_DOS%dos_nediv
    if(.not.PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,1x,F16.8)')e_range(ie), dos_data(ie,1)
    elseif(.not.PINPT%flag_collinear) then
      write(pid_dos,'(F16.8,1x,2F16.8)')e_range(ie), dos_data(ie,1:2)
    endif
  enddo

  close(pid_dos)

return
endsubroutine

