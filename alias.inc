#define if_main            if(myid .eq. 0)
#define if_main_then       if(myid .eq. 0) then
#define if_main_end        endif
#define if_nsp2            if(PINPT%nspin .eq. 2)
#define if_ncol            if(PINPT%ispinor .eq. 2)
#define mpi_kill           call MPI_Abort(mpi_comm_earth, 0, mpierr)
#define if_test            if(PINPT%flag_tbfit_test)
#define checkXXX           write(6,*)'XXX '
#define checkZZZ           write(6,*)'ZZZ '
#define checkVVV           write(6,*)'VVV '
#define checkJJJ           write(6,*)'JJJ '
#define checkKKK           write(6,*)'KKK '
#define checkAAA           write(6,*)'AAA '
#define write_msg          call write_log(message,3,myid)
#define write_msg_file     call write_log(message,1,myid)
#define write_msg_screen   call write_log(message,2,myid)
#define write_msg_all      call write_log(message,13,myid)
#ifdef MPI
#define kill_job           call MPI_BARRIER(mpi_comm_earth, mpierr); call mpi_finish()
#else
#define kill_job           stop
#endif
