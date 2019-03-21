#include "alias.inc"
module mpi_setup

#ifdef MPI
   use mpi_basics
#if SCALAPACK 
   use blacs_basics
#endif
   implicit none
   integer, allocatable :: task_list(:)
   type(mpicomm), target:: COMM_EARTH
   type(mpicomm), target:: COMM_ASIA
   type(mpicomm), target:: COMM_KOREA

#else
   integer*4, public            :: myid   = 0
   integer*4, public            :: nprocs = 1
   integer*4, public            :: npar   = 1
   integer*4, public            :: kpar   = 1
   integer*4, public            :: nproc_per_band = 1

   integer*4, public            :: myid_blacs = 0
   integer*4, public            :: nprow = 1
   integer*4, public            :: npcol = 1
   integer*4, public            :: myrow = 0
   integer*4, public            :: mycol = 0

   type mpicomm
        integer, public            :: myid = 0
        integer, public            :: nprocs = 1
   endtype

!  integer*4, public, parameter :: myid   = 0
!  integer*4, public, parameter :: nprocs = 1
!  integer*4, public, parameter :: npar   = 1
!  integer*4, public, parameter :: kpar   = 1
!  integer*4, public, parameter :: nproc_per_band = 1

!  integer*4, public, parameter :: myid_blacs = 0
!  integer*4, public, parameter :: nprow = 1
!  integer*4, public, parameter :: npcol = 1
!  integer*4, public, parameter :: myrow = 0
!  integer*4, public, parameter :: mycol = 0

!  type mpicomm
!       integer, public, parameter :: myid = 0
!       integer, public, parameter :: nprocs = 1
!  endtype

   type(mpicomm), target:: COMM_EARTH
   type(mpicomm), target:: COMM_ASIA
   type(mpicomm), target:: COMM_KOREA

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
       write(6,'(A,I0,A)') " MPI-parallelism will be employed. Running on ", &
                           nprocs," total cores."
     endif

#ifdef SCALAPACK
     call mpi_division()
#else
     ! NOTE: in the current version, we did not consider eigenvalue parallization
     !       if -DSCALAPACK is not activated, i.e., only k-point parallization will
     !       be performed unless -DMPI is activated.
     npar = 1
     kpar = 1
#endif

     return
   endsubroutine

   subroutine mpi_division()
     integer*4      mpierr

     call get_npar_kpar()
     call blacs_initialize()

     return
   endsubroutine

   subroutine blacs_initialize()
     integer*4   mpierr
     integer*4   npar_, kpar_     

     COMM_EARTH%mpi_comm = mpi_comm_earth
     COMM_EARTH%nprocs   = nprocs
     COMM_EARTH%myid     = myid
     npar_               = npar
     kpar_               = kpar

     call mpi_divide(COMM_EARTH, COMM_ASIA, COMM_KOREA, npar_, kpar_)

     if_main write(6,'(A,2(I0,A))')' Each k-point on ', COMM_ASIA%nprocs, &
                                   ' cores, ',npar,' groups.'

kill_job

     return
   endsubroutine

! NOTE: this subroutine is copied and modified from "M_initc" subroutine of VASP code
   subroutine mpi_init_comm(COMM)
     integer*4   mpierr
     type(mpicomm) :: COMM

     call MPI_COMM_RANK(COMM%mpi_comm, COMM%myid, mpierr)
     if(mpierr .ne. MPI_SUCCESS) then
       write(6,'(A)')' Error in MPI_COMM_RANK : mpi_init_comm'
       kill_job
     endif
     call MPI_COMM_SIZE(COMM%mpi_comm, COMM%nprocs, mpierr)
     if(mpierr .ne. MPI_SUCCESS) then
       write(6,'(A)')' Error in MPI_COMM_SIZE : mpi_init_comm'
       kill_job
     endif

     call MPI_BARRIER(COMM%mpi_comm, mpierr)
     if(mpierr .ne. MPI_SUCCESS) then
       write(6,'(A)')' Error in MPI_BARRIER : mpi_init_comm'
       kill_job
     endif

     return
   endsubroutine

! NOTE: this subroutine is copied and modified from "M_divide" subroutine of VASP code
   subroutine mpi_divide(COMM_EARTH, COMM_ASIA, COMM_KOREA, npar_, kpar_)
     integer*4   mpierr
     integer*4   npar_, kpar_
     integer*4   dims(2), dims_adam(2), dims_eve(2)
     logical     flag_period(2), flag_dims_adam(2), flag_dims_eve(2), reorder
     logical     grp_row, grp_col
     type(mpicomm) :: COMM_EARTH, COMM_MARS, COMM_ASIA, COMM_KOREA

     if(npar_ > COMM_EARTH%nprocs) then 
       if_main write(6,'(A)')' Error in mpi_divide: NPAR >= NPROCS'
       kill_job
     endif

     ! create new communicator for parents with cartesian topology [dim1, dim2]
     flag_period = .false.
     reorder     = .false.
     dims(1) = npar_
     dims(2) = COMM_EARTH%nprocs / npar_

     if(dims(1) * dims(2) .ne. COMM_EARTH%nprocs) then
       write(6,'(A,I0,A,I0)') &
                 " mpi_divide: can't subdivide ",COMM_EARTH%nprocs,' cpus by ', npar_
       write(6,'(A,I0,I0)')' Please check your NPAR or KPAR setting. Exit...', dims(1), dims(2)
       kill_job
     endif
     call MPI_CART_CREATE(COMM_EARTH%mpi_comm, 2,(/2, 2/),(/.false.,.false./),.false., &
                          COMM_MARS%mpi_comm, mpierr)
     call mpi_init_comm(COMM_MARS); COMM_EARTH = COMM_MARS

     ! create new communicator for "ASIA"
     grp_row = .false.  
     grp_col = .true.   
     call MPI_CART_SUB(COMM_EARTH%mpi_comm, (/grp_row, grp_col/), &
                       COMM_ASIA%mpi_comm, mpierr)
     call mpi_init_comm(COMM_ASIA)

     checkAAA , COMM_ASIA%mpi_comm, COMM_ASIA%nprocs, COMM_ASIA%myid, myid
kill_job

     return
   endsubroutine

   subroutine get_npar_kpar()
     integer*4      pid
     integer*4      i_continue, linecount, mpierr
     character*132  inputline
     character*40   desc_str
     logical        flag_fail

     flag_fail = .false.
     pid = 78

     ! set default values
     npar = 1
     kpar = 1

     if(myid .eq. 0) then
       open (pid, file='INCAR-TB', iostat=i_continue)
       do
         read(pid, '(A)', iostat=i_continue) inputline
         if(i_continue < 0) exit
         if(i_continue > 0) then
           write(6,*)'Unknown error reading file: mpi_division'
           flag_fail = .true. ; exit
         endif

         read(inputline,*,iostat=i_continue) desc_str
         if(i_continue .ne. 0) cycle              ! skip empty line
         if (desc_str(1:1).eq.'#') cycle  ! skip comment

         select case (desc_str)
           case('NPAR')
             read(inputline,*,iostat=i_continue) desc_str, npar
           case('KPAR')
             read(inputline,*,iostat=i_continue) desc_str, kpar
         endselect

       enddo
       close(pid)
     endif

     call MPI_BCAST(flag_fail, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
     call MPI_BCAST(npar     , 1, MPI_INTEGER, 0, mpi_comm_earth, mpierr)
     call MPI_BCAST(kpar     , 1, MPI_INTEGER, 0, mpi_comm_earth, mpierr)
     if(flag_fail) then
       kill_job
     endif

     return
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

     if_main write(6,'(A)') ' MPI-parallelism will be finished. Good luck.'

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
