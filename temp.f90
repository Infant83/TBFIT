#include "alias.inc"
program get_ldos2
   use parameters, only : zi 
   use mpi_setup
   use time
   implicit none
   real*8                    t_start, t_end  
   integer*4                 mpierr
   integer*4                 nediv, nkp, nemax
   integer*4                 nkdiv
   integer*4                 ii, i, j, k, ie, my_k
   integer*4                 d
   integer*4                 ispinor
   integer*4                 natom, natom_ldos,  nspec, norb
   integer*4                 my_pid
   integer*4                 narg,iarg
   integer*4, allocatable :: ne_found(:)
   integer*4, allocatable :: ne_found_(:)
   integer*4, allocatable :: natom_spec(:), norb_spec(:), norb_atom(:), iatom_ldos(:)
   real*8                    erange_i, erange_f, sigma
   real*8                    dos_, ldos_, ldos2_
   real*8                    weight
   real*8                    dk
   real*8,    allocatable :: E(:,:), V(:,:)
   real*8,    allocatable :: E_(:,:)
   real*8,    allocatable :: e_range(:)
   real*8,    allocatable :: ldos(:,:), dos(:)
   real*8,    allocatable :: ldos__(:,:), dos__(:)
   real*8,    allocatable :: ldos_k(  :,:), ldos2_k(:,:), dos_k(:,:)
   real*8,    allocatable :: ldos_k_(  :,:), ldos2_k_(:,:), dos_k_(:,:)
   real*8,    allocatable :: kp(:,:), kline(:)
   real*8,    allocatable :: kp_(:,:)
   real*8,    allocatable :: qpi(:,:), q(:)
   real*8,    external    :: fgauss
   logical                   flag_get_ldos, flag_get_qpi
   character*80              fname_in, fname, header
   character*20,external  :: int2str
!  integer*4, allocatable :: ourjob(nprocs), ourjob_disp(0:nprocs-1)
   integer*4, allocatable :: ourjob(:), ourjob_disp(:)
#ifdef MPI 
   call mpi_initialize()
#endif
   allocate(ourjob(nprocs))
   allocate(ourjob_disp(0:nprocs-1))
!  parse basic information
   ispinor = 2 ! SOC case
   fname_in = 'ldos2.in'
   narg = iargc()
   if(narg .eq. 1) call getarg(1, fname_in)

   open(999, file=trim(fname_in), status='old')
     if_main write(6,'(A)')"-- title of the calculation (without blank)"
       read(999,*)header
       if_main write(6,'(A,A   )')"   CASE = ",trim(header)
     if_main write(6,'(A)')"-- SOC (yes:2, no:1)?" 
       read(999,*)ispinor
       if_main write(6,'(A,(I6))')"   ISPINOR=",ispinor
     if_main write(6,'(A)')"-- NEIG? (how many eigenvalues?, write ne_max written in band_structure_TBA.kp_xx.dat)"
       read(999,*)nemax
       if_main write(6,'(A,(I6))')"   NEMAX=",nemax
     if_main write(6,'(A)')"-- NSPEC ? (ex, if the system is consist of carbon and oxygen: 2)"
       read(999,*)nspec
       if_main write(6,'(A,(I6))')"   NSPEC=",nspec
       allocate(natom_spec(nspec))
       allocate(norb_spec(nspec))
     if_main write(6,'(A)')"-- NATOM for each species ? (ex, 2 carbon 2 oxygen: 2 2)"
       read(999,*)natom_spec(1:nspec)
       if_main write(6,'(A,*(I6))')"   NATOM_SPEC=" ,natom_spec
       natom = sum(natom_spec)
       allocate(norb_atom(natom))
     if_main write(6,'(A)')"-- NORB for each species ? (ex, carbon (px, py, pz)  oxygen (s px py pz) : 3 4)"
       read(999,*)norb_spec(1:nspec)
       if_main write(6,'(A,*(I6))')"   NORB_SPEC=" ,norb_spec
     if_main write(6,'(A)')"-- NKP ?"
       read(999,*)nkp
       if_main write(6,'(A,*(I6))')"   NKP = " ,nkp
       allocate(kp(3,nkp)) ;kp = 0d0
       allocate(kp_(3,nkp));kp_= 0d0
       allocate(kline(nkp));kline = 0d0
       allocate(E(nemax,nkp)) ; E = 0d0
       allocate(E_(nemax,nkp)) ; E_ = 0d0
       allocate(ne_found(nkp)); ne_found = 0
       allocate(ne_found_(nkp)); ne_found_ = 0
     if_main write(6,'(A)')"-- NKDIV ?"
       read(999,*)nkdiv
       if_main write(6,'(A,*(I5))')"   NKDIV = ",nkdiv
     if_main write(6,'(A)')"-- Gaussian smearing?"
       read(999,*)sigma
       if_main write(6,'(A,(F10.4),A)')"   SMEARING = ",sigma, ' (eV)'
     if_main write(6,'(A)')"-- Number of division? (NEDIV)"
       read(999,*)nediv
       if_main write(6,'(A,(I6))')"   NEDIV = ",nediv
       allocate( dos_k(nediv,nkp)) ; dos_k = 0d0
       allocate( dos_k_(nediv,nkp)) ; dos_k_ = 0d0
       allocate( dos(nediv)) ; dos = 0d0
       allocate( dos__(nediv)) ; dos__ = 0d0
     if_main write(6,'(A)')"-- Energy range for dos plot (erange_i  erange_f)"
       read(999,*)erange_i, erange_f
       if_main write(6,'(A,F10.4,A,F10.4)')"   ERANGE = ",erange_i,' : ',erange_f
     if_main write(6,'(A)')"-- GET LDOS ? (evaluate LDOS from band_structure_TBA.dat with specified energy range)"
       read(999,*)flag_get_ldos
       if_main write(6,'(A,(L1))')"   GET_LDOS = ", flag_get_ldos
     if(flag_get_ldos) then
       allocate(ldos_k(      nediv,nkp)) ; ldos_k = 0d0
       allocate(ldos_k_(      nediv,nkp)) ; ldos_k_ = 0d0
       allocate(   qpi(      nediv,nkp)) ;    qpi = 0d0
       allocate(     q(            nkp)) ;      q = 0d0
       allocate(ldos2_k(nemax,nkp)) ; ldos2_k = 0d0
       allocate(ldos2_k_(nemax,nkp)) ; ldos2_k_ = 0d0
       allocate(ldos(natom,nediv)) ; ldos = 0d0
       allocate(ldos__(natom,nediv)) ; ldos__ = 0d0
       if_main write(6,'(A)')"-- Number of atoms for LDOS plot "
         read(999,*)natom_ldos
         if_main write(6,'(A,(I6))')"   NATOM_LDOS = ",natom_ldos
         allocate(iatom_ldos(natom_ldos))
       if_main write(6,'(A)')"-- Atom index for ldos plot"
         read(999,*)iatom_ldos(1:natom_ldos)
         if_main write(6,'(A,*(I5))')"   LDOS ATOM INDEX = ",iatom_ldos
     endif
     if_main write(6,'(A)')"-- GET QPI ? "
       read(999,*)flag_get_qpi 
       if_main write(6,'(A,(L1))')"   GET_QPI = ", flag_get_qpi
     if(flag_get_qpi) then
       if_main write(6,'(A)')"-- Weight "
         read(999,*)weight
         if_main write(6,'(A,F10.4        )')"   WEIGHT = ",weight
     endif
     
     close(999)
!  setup system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! set erange
   e_range =  erange_i + (/(k, k=0,nediv-1)/) * (erange_f - erange_i)/dble(nediv - 1)

   ! define number of orbitals basis for each atom
   do j = 1, nspec
     do i = 1, natom_spec(j)
       k = sum(natom_spec(1:j)) - natom_spec(j) + i
       norb_atom(k) = norb_spec(j)
     enddo
   enddo
   norb = sum(norb_atom)
   allocate(V(norb,nemax)) ; V = 0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call mpi_job_distribution_chain(nkp, ourjob, ourjob_disp)
  !do k = 1, nkp
   do k = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
     my_k = k - sum(ourjob(1:myid))
     if(myid .eq. 0) write(6,'(A,F10.4,A)')" STAT KP: ", real(my_k)/real(sum(ourjob(1:myid+1)))*100,' %'
     fname = 'band_structure_TBA'//'.kp_'//trim(ADJUSTL(int2str(k)))//'.dat'
     E(:,k) = -999d9
     call load_band(fname, nemax, norb, ne_found(k), E(:,k), V, kp(:,k), flag_get_ldos)
     do ie = 1, nediv
       do i = 1, ne_found(k)
         dos_ = fgauss(sigma, e_range(ie) - E(i,k)) / real(nkp)
         dos(ie) = dos(ie) + dos_
         dos_k(ie,k) = dos_k(ie,k) + dos_

         if(flag_get_ldos) then
           do j = 1, natom_ldos
             ii = sum(norb_atom(1:iatom_ldos(j))) - norb_atom(iatom_ldos(j)) + 1
             ldos_ =  dos_ * sum(V(ii:ii+norb_atom(iatom_ldos(j))-1,i))
             ldos(j,ie) = ldos(j,ie) + ldos_
             ldos_k(ie,k) = ldos_k(ie,k) + ldos_

             if(ie .eq. 1) then
               ldos2_ = sum(V(ii:ii+norb_atom(iatom_ldos(j))-1,i))
               ldos2_k(i,k) = ldos2_k(i,k) + ldos2_
             endif
           enddo
         endif

       enddo
     enddo
   enddo

#ifdef MPI
   call MPI_ALLREDUCE(E, E_,size(E),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   E = E_
   call MPI_ALLREDUCE(ne_found, ne_found_, size(ne_found), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
   ne_found = ne_found_
   call MPI_ALLREDUCE(kp, kp_, size(kp), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   kp = kp_
   if(flag_get_ldos) then
     call MPI_ALLREDUCE(dos, dos__,size(dos), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     dos = dos__
     call MPI_ALLREDUCE(dos_k,dos_k_,size(dos_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     dos_k = dos_k_
     call MPI_ALLREDUCE(ldos, ldos__,size(ldos), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos = ldos__
     call MPI_ALLREDUCE(ldos_k, ldos_k_, size(ldos_k), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos_k = ldos_k_
     call MPI_ALLREDUCE(ldos2_k, ldos2_k_, size(ldos2_k),MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
     ldos2_k = ldos2_k_
   endif
#endif
   if_main write(6,'(A)')" Replot DOS ... Done!"

   call get_kline_dist(kp, nkp, kline)
   ! write band structure
   if(myid .eq. 0) then
     open(999,file='band_structure_TBA.total.dat',status='unknown')
     do i = 1, nemax
       write(999,'(A,I0,A)')'# ',i, ' -th eigen'
       do k = 1,nkp
         if(flag_get_ldos) then
           write(999,'(2F16.6, *(F10.4))')kline(k), E(i,k), ldos2_k(i,k)
         else
           write(999,'(2F16.6          )')kline(k), E(i,k)
         endif
       enddo
       write(999,*)' '
       write(999,*)' '
     enddo
     write(6,'(A)')' -> Band structure (+atom projected) written: band_structure_TBA.total.dat'
     close(999)

     ! write DOS
     open(999,file='DOS.total.dat',status='unknown')
     if(flag_get_ldos) then
       open(998,file='LDOS.'//trim(header)//'.dat',status='unknown')
       write(998,'(A,*(1x,I0))')'# LDOS for atoms: ',iatom_ldos(:)
     endif
     do i =1, nediv
       write(999, '(2F16.4)') e_range(i), dos(i)
       if(flag_get_ldos) then
         write(998, '(3F16.4)') e_range(i), dos(i), sum(ldos(:,i))
       endif
     enddo
     write(6,'(A)')' -> DOS written: DOS.total.dat and LDOS.'//trim(header)//'.dat'
     close(999) 
     if(flag_get_ldos) close(998)

     dk=kline(2) - kline(1)
     if(flag_get_qpi .and. flag_get_ldos) then
       ! calculate QPI-like pattern
       do i=1, nkp
         if(i .ne. nkdiv+1) then
           do j=1,nkp
             if(i+j .ne. nkdiv+1 .and. i+j .ne. nkdiv*2+1 .and. i+j .ne. nkdiv*3+1) then
               if(i+j .gt. nkp) then
                 d = nint(kline(nkp)+kline((i+j - nkp)) - kline(i))/dk + 1
                 q(d) = kline(nkp)+kline((i+j - nkp)) - kline(i)
                 do ie=1, nediv
                   qpi(ie,d) = qpi(ie,d) + ldos_k(ie,i) * ldos_k(ie,i+j-nkp) * weight
                 enddo
              write(6,*)"XXX ", d,q(d)
               else
                 d = nint(kline(j)/dk) + 1
                 q(d) = kline(j)
                 do ie=1, nediv
                   qpi(ie,d) = qpi(ie,d) + ldos_k(ie,i) * ldos_k(ie,i+j) * weight
                 enddo
               endif
             endif
           enddo
         endif
       enddo

!      do i=1, nkp
!        do j=1, nkp
!          if( i .ne. nkdiv+1 .and. j .ne. nkdiv+1) then
!            if(j .ge. i) then
!              if(j .ge. nkdiv + 1 .and. i .lt. nkdiv+1) then
!                d = j - i + 1 - 1
!              else
!                d = j - i + 1
!              endif
!              q(d) = kline(j) - kline(i)
!              do ie = 1, nediv
!                qpi(ie,d)=qpi(ie,d) + ldos_k(ie,i) * ldos_k(ie,j) * weight
!              enddo
!            endif
!          endif
!        enddo
!      enddo

       open(999,file='QPI.'//trim(header)//'.dat',status='unknown')
       write(999,'(A,*(1x,I0))')'# QPI arround local atoms: ',iatom_ldos(:)
         do d =1, nkp-1
           do ie=1, nediv
             write(999, '(F16.8, F16.8, F16.8)')q(d), e_range(ie), qpi(ie,d)
           enddo
           write(999,*)' '
        !  write(999,*)' '
         enddo
       close(999)
       write(6,'(A)')' -> QPI written: QPI.'//trim(header)//'.dat '

     endif  
   endif

#ifdef MPI
  call MPI_BARRIER(mpi_comm_earth, mpierr)
  call mpi_finish()
#endif
  stop
endprogram

subroutine load_band(fname, nemax, norb, ne_found, E, V, kp, flag_get_ldos)
   use mpi_setup
   implicit none
   character*80    fname
   character*256   inputline
   integer*4       my_pid !, myid
   integer*4       mpierr, i_continue
   integer*4       norb, nemax, ne_found
   integer*4       nskip
   integer*4       ie, ik, nk
   integer*4       idummy
   real*8          kp(3)
   real*8          E(nemax)
   real*8          V(norb,nemax)
   logical         flag_go, flag_get_ldos
   integer*4,external ::  nitems
   E = -999d0 ; V = 0d0   
   nk = 1  ! for each file only one k-points has been written
   nskip = 3 ! for k-grid:3, kline:1
!  myid = 0
   my_pid = 100 + myid

   open(my_pid,file=trim(fname), status='old', iostat=i_continue)
   read(my_pid, '(A)') inputline
   read(my_pid, '(A)') inputline
   read(my_pid, *    ) inputline, inputline, ne_found

   if(ne_found .eq. 0) then 
     flag_go = .false.
     do while (.not.flag_go)
       read(my_pid,'(A)')inputline
       idummy = nitems(inputline)
       if(idummy .le. 0) then
         flag_go = .false.
       elseif(idummy .gt. 0) then
         flag_go = .true.
         backspace(my_pid)
       endif
     enddo
     read(my_pid,*) kp(:)    
     return
   else
     do ie = 1, ne_found
       flag_go = .false.
       do while (.not.flag_go)
         read(my_pid,'(A)')inputline
         idummy = nitems(inputline)
         if(idummy .le. 0) then
           flag_go = .false.
         elseif(idummy .gt. 0) then
           flag_go = .true.
           backspace(my_pid)
         endif
       enddo
     
       do ik = 1, nk
         if(flag_get_ldos) then
           read(my_pid,*) kp(:), E(ie), V(1:norb,ie)
         else
           read(my_pid,*) kp(:), E(ie)
         endif
       enddo
     enddo
     close(my_pid)
   endif
return
endsubroutine
function fgauss(sigma, x)
   implicit none
   real*8   sigma,sigma2
   real*8   x,xx
   real*8   fgauss
   real*8 pi, pi2
   pi = 4.d0*atan(1.d0)
   pi2=pi*2d0

   xx = x**2
   sigma2 = sigma**2

   fgauss= exp(-0.5d0*xx/sigma2)/(sigma*sqrt(pi2))

return
end function

function int2str(w) result(string)
  implicit none
! character(*), intent(out) :: string
  character*20  string
  integer*4,    intent(in)  :: w

  write(string,*) w

  return
endfunction

function nitems(string)
  implicit none
  logical blank
  integer*4 nitems,l,i
  character(*),intent(in) :: string
  nitems=0
  l=len_trim(string)
  blank = .true.
  do i=1,l
   if(string(i:i) .eq. '#') exit

   if (blank .and. string(i:i) .ne. ' ' ) then
     blank=.false.
     nitems=nitems + 1
   elseif( .not. blank .and. string(i:i) .eq. ' ') then
     blank=.true.
   endif
  enddo
  return
endfunction

subroutine get_kline_dist(kp, nkp, kline)
   implicit none
   integer*4    ik, nkp
   real*8       kline(nkp),k0(3), enorm
   real*8       kp(3,nkp)
   external     enorm

   do ik=1,nkp
     if(ik .eq. 1) then
      k0=kp(:,1)
      kline(1)=0
     else
      k0=kp(:,ik-1)
      kline(ik)=kline(ik-1)
     endif
     kline(ik)=kline(ik)+ enorm(3, kp(:,ik)-k0(:) )
   enddo

return
endsubroutine

function enorm ( n, x )
  implicit none

  integer*4 n
  real*8 x(n),enorm

  enorm = sqrt ( sum ( x(1:n) ** 2 ))
  return
end

