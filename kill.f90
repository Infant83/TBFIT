!module kill
! use mpi_setup
! use mykind
! implicit none

! contains

    subroutine kill_job()
        use mpi_setup
        use mykind
        implicit none
        integer(kind=sp)    mpierr

#ifdef MPI
        if(COMM_KOREA%flag_split) then
            call MPI_BARRIER(COMM_KOREA%mpi_comm, mpierr)
            call mpi_finish()
        elseif(.not. COMM_KOREA%flag_split) then
            call MPI_BARRIER(mpi_comm_earth, mpierr)
            call mpi_finish()
        endif
#else

        stop

#endif        

    endsubroutine

!endmodule
