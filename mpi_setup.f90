module mpi_setup

   use mpi_basics

   implicit none
   integer, allocatable :: task_list(:)
   contains

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
       write(6,*) " MPI-parallelism will be employed. myid=",myid
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
   endsubroutine 

   subroutine mpi_job_distribution(njob_total, myjob_init, myjob_fina)
     implicit none
     integer ::   njob_total, myjob_init, myjob_fina
     integer ::   njob_each    

     njob_each = njob_total / nprocs


   endsubroutine
endmodule
