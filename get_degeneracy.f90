#include "alias.inc"
subroutine get_D(E, D, n, imode)
   use do_math, only : degen
   implicit none
   integer*4                n
   integer*4, intent(in) :: imode
   real*8                   E(n), D(3,n)

 ! for nonmag
   if(imode .eq. 1) then
     D(1:3,1:n)     = degen(E(1    :n  ))

 ! for collinear
   elseif(imode .eq. 2) then
     D(1:3,1:n/2)   = degen(E(1    :n/2)) ! up
     D(1:3,n/2+1:n) = degen(E(n/2+1:n  )) ! dn

 ! for noncollinear
   elseif(imode .eq. 4) then
     D(1:3,1:n)     = degen(E(1    :n  ))

   endif

   return
endsubroutine

subroutine get_degeneracy(EN, n, nkp, PINPT)
   use parameters, only : incar, energy
   use mpi_setup
   use do_math, only : degen
   implicit none
   type (incar)           :: PINPT
   type (energy)          :: EN   
   integer*4                 n, nkp
   real*8                    E(n,nkp)
   real*8                    D(3,n,nkp) ! (1,:,:)  : dE(above) * dE(below)
                                        ! (2,:,:)  : dE(above) = E(n+1) - E(n  )
                                        ! (3,:,:)  : dE(below) = E(n  ) - E(n-1)
!  real*8                    D_(3,n,nkp)
   integer*4                 ik
   integer*4                 mpierr, ncpu, id
  !integer*4                 ourjob(nprocs), ourjob_disp(0:nprocs-1)
   integer*4, allocatable :: ourjob(:), ourjob_disp(:)

   if(.not. PINPT%flag_fit_degeneracy) return
   E = EN%E
   D = 0d0
  !D_= 0d0
 
   if(COMM_KOREA%flag_split) then
     ncpu = COMM_KOREA%nprocs
     id   = COMM_KOREA%myid
   else
     ncpu = nprocs
     id   = myid
   endif
   allocate(ourjob(ncpu))
   allocate(ourjob_disp(0:ncpu-1))

  !call mpi_job_ourjob(nkp, ourjob)
   call mpi_job_distribution_chain(nkp, ncpu, ourjob, ourjob_disp)
   
   do ik = sum(ourjob(1:myid)) + 1, sum(ourjob(1:myid+1))
     call get_D(E(:,ik), D(1:3,:,ik), n, PINPT%ispin*PINPT%ispinor)
   enddo

#ifdef MPI
   if(COMM_KOREA%flag_split) then
     call MPI_ALLREDUCE(D, EN%D, size(D), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
   elseif(.not. COMM_KOREA%flag_split) then
     call MPI_ALLREDUCE(D, EN%D, size(D), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
   endif
#else
   EN%D = D
#endif

   return
endsubroutine

subroutine get_degeneracy_serial(E,D, n, nkp, PINPT)
   use parameters, only : incar
   use do_math, only : degen
   implicit none
   type (incar)           :: PINPT
   integer*4                 n, nkp
   real*8                    E(n,nkp)
   real*8                    D(3,n,nkp)
   integer*4                 ik

   D = 0d0

   do ik = 1, nkp
     call get_D(E(:,ik), D(1:3,:,ik), n, PINPT%ispin*PINPT%ispinor)
   enddo

   return
endsubroutine
