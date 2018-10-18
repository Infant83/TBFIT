#define if_main            if(myid .eq. 0)
#define if_main_then       if(myid .eq. 0) then
#define if_main_end        endif
#define if_nsp2            if(PINPT%nspin .eq. 2)
#define if_ncol            if(PINPT%ispinor .eq. 2)
#define mpi_kill           call MPI_Abort(mpi_comm_earth, 0, mpierr)
#define if_test            if(PINPT%flag_tbfit_test)
