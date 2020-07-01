#include "alias.inc"
subroutine get_D(E, D, n, imode)
   use do_math, only : degen
   implicit none
   integer*4                n
   integer*4, intent(in) :: imode
   real*8                   E(n), D(3,n)

   if(imode .eq. 1) then
     ! for nonmag
     D(1:3,1:n)     = degen(E(1    :n  ))
   elseif(imode .eq. 2) then
     ! for collinear
     D(1:3,1:n/2)   = degen(E(1    :n/2)) ! up
     D(1:3,n/2+1:n) = degen(E(n/2+1:n  )) ! dn
   elseif(imode .eq. 4) then
     ! for noncollinear
     D(1:3,1:n)     = degen(E(1    :n  ))
   endif

   return
endsubroutine

subroutine get_degeneracy(E,D, n, nkp, PINPT)
   use parameters, only : incar
   use mpi_setup
   use do_math, only : degen
   implicit none
   type (incar)           :: PINPT
   integer*4                 n, nkp
   real*8                    E(n,nkp)
   real*8                    D(3,n,nkp) ! (1,:,:)  : dE(above) * dE(below)
                                        ! (2,:,:)  : dE(above) = E(n+1) - E(n  )
                                        ! (3,:,:)  : dE(below) = E(n  ) - E(n-1)
   real*8                    D_(3,n,nkp)
   integer*4                 ik
   integer*4                 ourjob(nprocs), ourjob_disp(0:nprocs-1)
   integer*4                 mpierr

   D = 0d0
   D_= 0d0
 
   call mpi_job_ourjob(nkp, ourjob)
   do ik = sum(ourjob(1:myid)) + 1, sum(ourjob(1:myid+1))
     call get_D(E(:,ik), D_(1:3,:,ik), n, PINPT%ispin*PINPT%ispinor)
   enddo

#ifdef MPI
   call MPI_ALLREDUCE(D_, D, size(D), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
#else
   D = D_
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
