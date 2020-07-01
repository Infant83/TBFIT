#include "alias.inc"
module reorder_band
   use do_math
   use mpi_setup
   use print_io

contains
   subroutine update_order(IDX, E, E_ORD, nkp, nband, PINPT)
      use parameters, only: incar
      implicit none
      type(incar  ) :: PINPT
      integer*4        nband, nkp
      integer*4        ik
      integer*4        mpierr 
      real*8           E(nband*PINPT%nspin,nkp)
      real*8           E_ORD(nband*PINPT%nspin,nkp)
      integer*4        IDX(nband*PINPT%nspin, nkp)
  
      do ik = 1, nkp
        E_ORD(:,ik)   = E(IDX(:,ik),ik)
      enddo

      return
   endsubroutine

   subroutine get_ordered_band(ETBA, nkp, neig, iband, nband, PINPT, flag_weight, PWGHT)
      use parameters, only: energy, incar, poscar, weight, pi, pi2, zi, rt2
      use print_matrix
      implicit none
      type(incar  ) :: PINPT
      type(weight ) :: PWGHT
      type(energy ) :: ETBA
      integer*4        neig, nband, iband, nkp
      integer*4        ik, ii, jj, i, msize
      integer*4        max_jj
      integer*4        mpierr
      real*8           E(nband*PINPT%nspin,nkp) ! original band
      real*8           E_ORD(nband*PINPT%nspin,nkp) ! re-ordered band
      complex*16       V(neig*PINPT%ispin,nband*PINPT%nspin,nkp)
      complex*16       V_ORD(neig*PINPT%ispin,nband*PINPT%nspin,nkp)
      real*8           MM(nband*PINPT%nspin,nband*PINPT%nspin,nkp)
      real*8           MM_(nband*PINPT%nspin,nband*PINPT%nspin,nkp)
      real*8           R(nband*PINPT%nspin,nband*PINPT%nspin)
      real*8           Mij, Mij_
      real*8           ov_cut
      integer*4        IDX(nband*PINPT%nspin, nkp) ! band index
      real*8           IDXR(nband*PINPT%nspin) ! band index
      character*20,external :: int2str
      logical          flag_order
      logical          flag_weight ! always .false. in the current version... HJK (19.06.2020)
      integer*4        ourjob(nprocs), ourjob_disp(0:nprocs-1)
     
      call mpi_job_distribution_chain(nkp-1, ourjob, ourjob_disp)

      msize = nband*PINPT%nspin
      R = diagonal(msize)
     !IDXR = real((/1:msize/))
      do ii=1, msize
        IDXR(ii) = real(ii)
      enddo
      ov_cut = PINPT%band_order_overlap_cutoff
#ifdef MPI
      call MPI_BCAST(ETBA%V, size(ETBA%V), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
#endif
      MM = 0d0
      Mij = 0d0

      if(.not. PINPT%flag_print_orbital) then
        write(message,'(A)')' '                                                                                ; write_msg
        write(message,'(A)')'    !WARN! LORBIT = .false. and LORDER = .TRUE., which is incompatible request.'  ; write_msg
        write(message,'(A)')'           Try again with LORBIT = .TRUE.'                                        ; write_msg
        write(message,'(A)')'           Exit program anyway...'                                                ; write_msg
        write(message,'(A)')' '                                                                                ; write_msg
        kill_job
      endif
   
! MPI does not applied, since it make it slow down
   kp:do ik = 2, nkp
      m:do ii = 1, msize
           Mij_ = 0d0
        n:do jj = 1, msize
             Mij = real(abs( dot_product(ETBA%V(:,ii, ik  ), ETBA%V(:,jj,ik-1)) ))
             if( Mij .gt. ov_cut .and. Mij .gt. Mij_) then
               Mij_ = Mij ! update if larger than sqrt(2)/2 (default)
               max_jj = jj
             endif
          enddo n
          if(Mij_ .eq. 0d0) then
            MM(ii,ii,ik) = 1d0
          else
            MM(ii,max_jj,ik) = real( nint(Mij_) )
          endif
        enddo m
      enddo kp

      ETBA%IDX(:,1) = int(IDXR)
      ETBA%E_ORD(:,1) = ETBA%E(:,1)
      do ik = 2, nkp
        R       = matprod (msize, msize, 'N', R, msize, msize, 'N', MM(:,:, ik))
        ! band swap by multiplying R_k by E(:,ik), where R_k = PI_k(MM(:,:,ik)) act as rotational matrix for E
        ETBA%IDX(:,ik) = int( matprod(msize, msize, 'N', R, msize, 'N', IDXR) )
        ETBA%E_ORD(:,ik)   = ETBA%E(ETBA%IDX(:,ik),ik)
        if(flag_weight) PWGHT%WT(:,ik)   = PWGHT%WT(ETBA%IDX(:,ik),ik)
      enddo

     !call print_matrix_r(MM(:,:,2),msize, msize, 'MM', 0, 'F6.3')

      if_main_then 
        ETBA%V_ORD(:,:,1) = ETBA%V(:,:,1)
        do ik = 2, nkp 
          ETBA%V_ORD(:,:,ik) = ETBA%V(:,ETBA%IDX(:,ik),ik) 
        enddo
      if_main_end

      return
   endsubroutine

endmodule
