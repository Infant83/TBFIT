module mpi_setup

#ifdef MPI
   use mpi_basics

   implicit none
   integer, allocatable :: task_list(:)

#else
   integer*4, public, parameter :: myid   = 0
   integer*4, public, parameter :: nprocs = 1
#endif

   contains

#ifdef MPI
   subroutine mpi_initialize()
     integer*4  mpierr

     call MPI_INIT(mpierr)

     flag_use_mpi = .false.

     if(mpierr .eq. -1) then
 
       flag_use_mpi = .false.
       nprocs = 0

     elseif(mpierr .ne. MPI_SUCCESS) then
       write(6,*) "Error in MPI initialization."
       write(6,*) "Exiting..."
       stop
     
     else
       flag_use_mpi = .true.
       myid = 0
     endif


     mpi_comm_earth = MPI_COMM_WORLD

     call get_my_task()

     if(flag_use_mpi .and. myid .eq. 0) then
       write(6,'(A,I0,A)') " MPI-parallelism will be employed. Running on ",nprocs," total cores."
     endif

   endsubroutine

   subroutine get_my_task()
     integer*4  mpierr
     
     call MPI_COMM_SIZE(mpi_comm_earth, nprocs, mpierr)
     call MPI_COMM_RANK(mpi_comm_earth, myid, mpierr)

     if (mpierr.ne.MPI_SUCCESS) then
       write(6,'(1X,A)') "* get_my_task() failed"
       write(6,'(1X,A)') "* Exiting..."
       stop
     endif
   endsubroutine

   subroutine get_my_procs(name, length)
     character(LEN=MPI_MAX_PROCESSOR_NAME) :: name
     integer :: length
     integer :: mpierr
     
     call MPI_GET_PROCESSOR_NAME(name, length, mpierr)

     if(mpierr .ne. MPI_SUCCESS) then
       write(6,'(1X,A)') "* get_my_processor() failed"
       write(6,'(1X,A)') "* Exiting..."
       stop
     endif
   endsubroutine

   subroutine mpi_finish()
     integer*4  mpierr

     call MPI_FINALIZE(mpierr)
     stop
   endsubroutine 

#endif

   subroutine mpi_job_distribution_chain(njob, ourjob)
     implicit none
     integer*4    njob
     integer*4    mynjob
     integer*4    cpuid, mpierr
     integer*4    nresidue
#ifdef MPI
     integer*4    ourjob(nprocs)
#else
     integer*4    ourjob(1)
!    integer*4    myid, nprocs
!    myid = 0
!    nprocs = 1
#endif

     mynjob = floor ( real(njob)/real(nprocs) )
     nresidue = nint (real(njob) - real(mynjob) * real(nprocs))

     do cpuid = 1, nprocs
       if( cpuid .le. nresidue ) then
          ourjob(cpuid) = mynjob + 1
       else
          ourjob(cpuid) = mynjob
       endif
     enddo

   endsubroutine

endmodule
