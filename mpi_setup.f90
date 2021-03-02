#include "alias.inc"
module mpi_setup
   use print_io
#ifdef MPI
   use mpi_basics
#if SCALAPACK 
   use blacs_basics
#endif
   implicit none
   integer, allocatable :: task_list(:)
   type(mpicomm), target:: COMM_EARTH
   type(mpicomm), target:: COMM_ASIA
   type(mpicomm), public, target:: COMM_KOREA

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
        integer, public         :: myid = 0
        integer, public         :: nprocs = 1
        logical, public         :: flag_split = .FALSE.
    
        integer, public         :: key = 0
        integer, public         :: color = 0
        integer, public, allocatable :: group_main(:)
   endtype

   type(mpicomm), target:: COMM_EARTH
   type(mpicomm), target:: COMM_ASIA
   type(mpicomm), public, target:: COMM_KOREA

#endif

   contains

!!!!!!! start if_def MPI
#ifdef MPI
   subroutine mpi_initialize(fnamelog)
     implicit none
     integer*4  mpierr
#ifdef SCALAPACK
     integer*4  CONTEXT
#endif
     character*40 fnamelog
     yourid = 99
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

     call open_log(fnamelog, myid)

     write(message,'(A,I0,A)')"#MPI-parallelism will be employed. Running on ", nprocs," total cores." 
     call write_log(message,3,myid)

#ifdef SCALAPACK
     call mpi_division()
     call proc_map() !need to be updated...
#else
     ! NOTE: in the current version, we did not consider eigenvalue parallization
     !       if -DSCALAPACK is not activated, i.e., only k-point parallization will
     !       be performed unless -DMPI is activated.
     npar = 1
     kpar = 1
#endif

     return
   endsubroutine

#ifdef SCALAPACK
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
     COMM_EARTH%npar     = npar
     COMM_EARTH%kpar     = kpar

     call mpi_divide(COMM_EARTH, COMM_ASIA)

     write(message,'(A,2(I0,A))')' Each k-point on ', COMM_ASIA%nprocs, ' cores, ',COMM_EARTh%npar,' groups.' ; write_msg
     return
   endsubroutine

! NOTE: this subroutine is copied and modified from "M_initc" subroutine of VASP code
   subroutine mpi_init_comm(COMM)
     integer*4   mpierr
     type(mpicomm) :: COMM

     call MPI_COMM_RANK(COMM%mpi_comm, COMM%myid, mpierr)
     if(mpierr .ne. MPI_SUCCESS) then
       write(message,'(A)')' Error in MPI_COMM_RANK : mpi_init_comm' ; write_msg
       call MPI_BARRIER(mpi_comm_earth, mpierr)
       call mpi_finish()
     endif
     call MPI_COMM_SIZE(COMM%mpi_comm, COMM%nprocs, mpierr)
     if(mpierr .ne. MPI_SUCCESS) then
       write(message,'(A)')' Error in MPI_COMM_SIZE : mpi_init_comm' ; write_msg
       call MPI_BARRIER(mpi_comm_earth, mpierr)
       call mpi_finish()
     endif

     call MPI_BARRIER(COMM%mpi_comm, mpierr)
     if(mpierr .ne. MPI_SUCCESS) then
       write(message,'(A)')' Error in MPI_BARRIER : mpi_init_comm' ; write_msg
       call MPI_BARRIER(mpi_comm_earth, mpierr)
       call mpi_finish()
     endif

     return
   endsubroutine

   subroutine proc_map()
     implicit none
     integer*4, allocatable :: imap(:,:)
     integer*4  NP_COL, MY_ROW, MY_COL
     integer*4  CONTEXT
     integer*4  CONTEXT_, NPROW_, NPCOL_, MYROW_, MYCOL_
     integer*4  i, j, k, mpierr
     integer*4  BLACS_PNUM
     external   BLACS_PNUM

     nprow = COMM_EARTH%dims(1)
     npcol = COMM_EARTH%dims(2)

     allocate(imap(nprow,npcol))
     allocate(id_blacs(COMM_EARTH%nprocs))
    
     CONTEXT_ = COMM_EARTH%mpi_comm

     call BLACS_GRIDINIT( CONTEXT_, 'Row', 1, COMM_EARTH%nprocs)
     call BLACS_GRIDINFO( CONTEXT_, NPROW_, NPCOL_, MYROW_, MYCOL_)
     myid_blacs = BLACS_PNUM(CONTEXT_, MYROW_, MYCOL_)
     call MPI_ALLGATHER(myid_blacs, 1, MPI_INTEGER4, id_blacs(1), 1, &
                        MPI_INTEGER4, COMM_EARTH%mpi_comm, mpierr)

     if(mpierr .ne. 0) then
       write(message,'(A,I0,A)')' ERROR ALLGATHER in "proc_map" (error code= ', mpierr,' )' ; write_msg
     endif

     k = 1
     do j = 1, npcol
       do i = 1, nprow
         imap(i,j) = id_blacs(k)
         k = k + 1
       enddo
     enddo

     call MPI_BARRIER(COMM_EARTH%mpi_comm, mpierr)
     call BLACS_GRIDEXIT(CONTEXT_)

!    call BLACS_GET(0,0,CONTEXT)
     CONTEXT = COMM_EARTH%mpi_comm
     call BLACS_GRIDMAP( CONTEXT, imap, nprow, nprow, npcol)
     call BLACS_GRIDINFO( CONTEXT, nprow, npcol, myrow, mycol)
     return
   endsubroutine

! NOTE: this subroutine is copied and modified from "M_divide" subroutine of VASP code
   subroutine mpi_divide(COMM_EARTH, COMM_ASIA)
     integer*4   mpierr
     integer*4   ndims, dims(2) ! two dimensional cartesian topology
     logical     flag_period(2), reorder ! topology without periodic boundary condition
     logical     grp_row, grp_col
     type(mpicomm) :: COMM_EARTH, COMM_MARS, COMM_ASIA

     if(COMM_EARTH%npar > COMM_EARTH%nprocs) then 
       write(message, '(A)')' Error in mpi_divide: NPAR >= NPROCS' ; write_msg
       call MPI_BARRIER(mpi_comm_earth, mpierr)
       call mpi_finish()
     endif

     ! create new communicator for parents with cartesian topology [dim1, dim2]
     flag_period = .false.
     reorder     = .false.
     ndims   = 2 ! 2 dimensional topology
     COMM_EARTH%dims(1) = COMM_EARTH%npar
     COMM_EARTH%dims(2) = COMM_EARTH%nprocs / COMM_EARTH%npar

     if(COMM_EARTH%dims(1) * COMM_EARTH%dims(2) .ne. COMM_EARTH%nprocs) then
       write(message,'(A,I0,A,I0)') " mpi_divide: can't subdivide ",COMM_EARTH%nprocs,' cpus by ', COMM_EARTH%npar    ; write_msg
       write(message,'(A,A,I0,A,I0,A)')' Please check your NPAR setting. Exit...', ' DIM(1:2)= (',COMM_EARTH%dims(1),',', COMM_EARTH%dims(2),')' ; write_msg
       call MPI_BARRIER(mpi_comm_earth, mpierr)
       call mpi_finish()
     endif
     call MPI_CART_CREATE(COMM_EARTH%mpi_comm, ndims,(/COMM_EARTH%dims(1), COMM_EARTH%dims(2)/), &
                          flag_period,reorder, COMM_MARS%mpi_comm, mpierr)
     call mpi_init_comm(COMM_MARS)
     COMM_EARTH%mpi_comm = COMM_MARS%mpi_comm
     COMM_EARTH%nprocs   = COMM_MARS%nprocs  
     COMM_EARTH%myid     = COMM_MARS%myid    

     ! Create cartesian coordinate for processors
     call MPI_CART_COORDS(COMM_EARTH%mpi_comm,COMM_EARTH%myid, ndims, COMM_EARTH%mycoord, mpierr)

     ! Create new communicator for "ASIA" subdividing "EARTH" with "column" major.
     grp_row = .false.  
     grp_col = .true.   
     call MPI_CART_SUB(COMM_EARTH%mpi_comm, (/grp_row, grp_col/), &
                       COMM_ASIA%mpi_comm, mpierr)
     call mpi_init_comm(COMM_ASIA) ! generate ASIA%myid, ASIA%nprocs
     COMM_ASIA%dims = COMM_EARTH%dims
     COMM_ASIA%mycoord = COMM_EARTH%mycoord
     return
   endsubroutine

#endif

   subroutine get_my_task()
     integer*4  mpierr
     
     call MPI_COMM_SIZE(mpi_comm_earth, nprocs, mpierr)
     call MPI_COMM_RANK(mpi_comm_earth, myid, mpierr)

     if (mpierr.ne.MPI_SUCCESS) then
       write(message,'(1X,A)') "* get_my_task() failed"  ; write_msg_all
       write(message,'(1X,A)') "* Exiting..."  ; write_msg_all
       stop
     endif
   endsubroutine

   subroutine get_my_procs(name, length)
     character(LEN=MPI_MAX_PROCESSOR_NAME) :: name
     integer :: length
     integer :: mpierr
     
     call MPI_GET_PROCESSOR_NAME(name, length, mpierr)

     if(mpierr .ne. MPI_SUCCESS) then
       write(message,'(1X,A)') "* get_my_processor() failed" ; write_msg_all
       write(message,'(1X,A)') "* Exiting..." ; write_msg_all
       stop
     endif
   endsubroutine

   subroutine mpi_finish()
     integer*4  mpierr

     write(message,'(A)') ' MPI-parallelism will be finished. Good luck.' ; write_msg
     call MPI_FINALIZE(mpierr)
     stop
   endsubroutine 

#endif  
!!!!!!! end if_def MPI

   ! NOTE: Anmeldung (Deutsch) = registration 
   ! This subroutine split current world communicator mpi_comm_earth into several group.
   ! As I'm come from Korea and working at Geermany, I need to registrate my color of eye
   ! and get the id. This is just for fun, kind of joke but inspiring my identity and
   ! refreshing wonderful working environment in Germany :) H.-J. Kim (FZJ, 25. Feb. 2021)
   subroutine mpi_comm_anmeldung(COMM_KOREA, ngroup, mygroup)
     implicit none
     type(mpicomm)::COMM_KOREA
     integer*4      ngroup, mpierr
     integer*4      ourgroup(ngroup)
     integer*4      mygroup(0:nprocs-1)
     integer*4      i, group_main(ngroup)

     COMM_KOREA%color = mygroup(myid)
     COMM_KOREA%key   = myid
     COMM_KOREA%flag_split = .TRUE.

     call MPI_COMM_SPLIT(mpi_comm_earth, COMM_KOREA%color, COMM_KOREA%key, COMM_KOREA%mpi_comm, mpierr)
     call MPI_COMM_RANK(COMM_KOREA%mpi_comm, COMM_KOREA%myid, mpierr)
     call MPI_COMM_SIZE(COMM_KOREA%mpi_comm, COMM_KOREA%nprocs, mpierr)

     if(allocated(COMM_KOREA%group_main)) deallocate(COMM_KOREA%group_main)
     allocate(COMM_KOREA%group_main(npar)) 
     COMM_KOREA%group_main = 0
     do i = 0, ngroup - 1
       if(COMM_KOREA%color .eq. i .and. COMM_KOREA%myid .eq. 0) then
         COMM_KOREA%group_main(i+1) = myid
       endif
     enddo
     
     call MPI_ALLREDUCE(COMM_KOREA%group_main, group_main, ngroup, MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
     COMM_KOREA%group_main = group_main
     return
   endsubroutine

   subroutine mpi_job_ourjob(njob, ourjob)
     implicit none 
     integer*4     njob, mynjob
     integer*4     ourjob(nprocs)
     integer*4     cpuid, mpierr
     integer*4     nresidue
 
     mynjob = floor ( real(njob)/real(nprocs) )
     nresidue = nint ( real(njob) - real(mynjob) * real(nprocs) )   
     ourjob = 0
     do cpuid = 1, nprocs   
       if( cpuid .le. nresidue ) then
         ourjob(cpuid) = mynjob + 1
       else
         ourjob(cpuid) = mynjob
       endif
     enddo

   endsubroutine

   subroutine mpi_job_distribution_group(ngroup, njob, ourgroup, mygroup, ourjob)
     implicit none
     integer*4      ngroup, nmember, nresidue
     integer*4      njob,   ngroupjob, nresidue_
     integer*4      groupid, cpuid, id
     integer*4      ourgroup(ngroup), mygroup(0:nprocs-1)
     integer*4      ourjob(ngroup)
     integer*4      mpierr

     ! ncpus per each group
     nmember   = floor( real(nprocs)/real(ngroup) )
     nresidue  = nint ( real(nprocs) - real(nmember) * real(ngroup) )

     ! njobs per each group
     ngroupjob = floor( real(njob)/real(ngroup) )
     nresidue_ = nint ( real(njob) - real(ngroupjob) * real(ngroup) )

     ! each group have nmember + alpha, alpha is the distributed from nresidue over groups
     do groupid = 1, ngroup
       if(groupid .le. nresidue) then
         ourgroup(groupid) = nmember + 1
       else
         ourgroup(groupid) = nmember
       endif
     enddo

     do groupid = 1, ngroup
       if(groupid .le. nresidue_) then
         ourjob(groupid) = ngroupjob + 1
       else
         ourjob(groupid) = ngroupjob
       endif
     enddo 
    
     cpuid = 0
     do groupid = 0, ngroup - 1
       do id = 1, ourgroup(groupid+1)
         mygroup(cpuid) = groupid
         cpuid = cpuid + 1
       enddo
     enddo

     return
   endsubroutine
   subroutine mpi_job_distribution_chain(njob, ncpu, ourjob, ourjob_disp)
     implicit none
     integer*4    njob, ncpu
     integer*4    mynjob
     integer*4    cpuid, mpierr
     integer*4    nresidue
     integer*4    ourjob(ncpu)
     integer*4    ourjob_disp(0:ncpu-1)

     mynjob = floor ( real(njob)/real(ncpu) )
     nresidue = nint (real(njob) - real(mynjob) * real(ncpu))
     ourjob = 0

     do cpuid = 1, ncpu
       if( cpuid .le. nresidue ) then
          ourjob(cpuid) = mynjob + 1
       else
          ourjob(cpuid) = mynjob
       endif
     enddo

    ourjob_disp(0) = 0
#ifdef MPI
    do cpuid = 1, ncpu-1
      ourjob_disp(cpuid)= ourjob_disp(cpuid - 1) + ourjob(cpuid)
    enddo
#endif

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
           write(message,'(A)') 'Unknown error reading file: get_npar_kpar' ; write_msg_all
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

     if(nprocs .eq. 1 .and. npar .gt. 1) then
       write(message,'(A)') ' !WARN! NPROCS = 1 and NPAR > NPROCS --> enforce NPAR = 1 ' ; write_msg_all
       npar = 1
     endif

#ifdef MPI
     call MPI_BCAST(flag_fail, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
     call MPI_BCAST(npar     , 1, MPI_INTEGER, 0, mpi_comm_earth, mpierr)
     call MPI_BCAST(kpar     , 1, MPI_INTEGER, 0, mpi_comm_earth, mpierr)
#endif
     if(flag_fail) then
       call MPI_BARRIER(mpi_comm_earth, mpierr)
       call mpi_finish()
     endif

     return
   endsubroutine

subroutine report_job_distribution(flag_stat, ourjob, jobname)
   implicit none
   integer*4    mpierr, i, id
   integer*4    ourjob(nprocs)
   logical      flag_stat
   character(len=80), optional, intent(in) :: jobname


   if(flag_stat) then
     if(present(jobname)) then
       write(message,'(A,A)')               '       JOB DISTRUBUTION for ',trim(jobname),' :' ; write_msg
     else
       write(message,'(A)')                 '       JOB DISTRUBUTION :' ; write_msg
     endif

     write(message,'(A,I0,A,I0,A)')         '       -> cpuid( ',myid,' ): ', ourjob(myid+1),' k-points'
#ifdef MPI
     call MPI_GATHER(message, 2048, MPI_CHARACTER, message_pack, 2048, MPI_CHARACTER, 0, mpi_comm_earth, mpierr)
#else
     message_pack(1) = message
#endif
     do i = 1, nprocs
       call write_log(trim(message_pack(i)), 3, myid)
     enddo

   endif

   return
endsubroutine

subroutine report_hostname()
   implicit none
   character(len=80) :: myhost
   integer*4            i, mpierr, len_host
   character*80         hosts(nprocs), hosts_
   character*20, external::int2str

   call HOSTNM(myhost)   

#ifdef MPI
   call MPI_GATHER(myhost, 80, MPI_CHARACTER, hosts, 80, MPI_CHARACTER, 0, mpi_comm_earth, mpierr)
#else
   hosts(1) = myhost  
#endif

   do i = 1, nprocs
     call write_log('| Executed on  '//trim(myhost)//' : NODE = '//trim(int2str(i-1)),3,myid)
   enddo

   return
endsubroutine
endmodule
