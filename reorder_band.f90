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

   subroutine get_ordered_band(ETBA, PKPTS, PGEOM, PWGHT, PINPT, flag_order_weight, flag_use_overlap)
      use parameters, only: energy, incar, poscar, kpoints, weight, pi, pi2, zi, rt2
      use print_matrix
      implicit none
      type(incar  ) :: PINPT
      type(weight ) :: PWGHT
      type(energy ) :: ETBA
      type(kpoints) :: PKPTS
      type(poscar ) :: PGEOM
      integer*4        ik, ii, jj, i, msize
      integer*4        max_jj
      integer*4        mpierr
      real*8           E(PGEOM%nband*PINPT%nspin,PKPTS%nkpoint) ! original band
      real*8           E_ORD(PGEOM%nband*PINPT%nspin,PKPTS%nkpoint) ! re-ordered band
      complex*16       V(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
      complex*16       SV(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
      complex*16       V_ORD(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
      complex*16       SV_ORD(PGEOM%neig*PINPT%ispin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
      real*8           MM(PGEOM%nband*PINPT%nspin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
      real*8           MM_(PGEOM%nband*PINPT%nspin,PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
      real*8           R(PGEOM%nband*PINPT%nspin,PGEOM%nband*PINPT%nspin)
      real*8           Mij, Mij_
      real*8           ov_cut
      integer*4        IDX(PGEOM%nband*PINPT%nspin, PKPTS%nkpoint) ! band index
      real*8           IDXR(PGEOM%nband*PINPT%nspin) ! band index
      character*20,external :: int2str
      logical          flag_order
      logical          flag_order_weight ! always .false. in the current version... HJK (19.06.2020)
      logical          flag_use_overlap
      integer*4        ourjob(nprocs), ourjob_disp(0:nprocs-1)
     
      call mpi_job_distribution_chain(PKPTS%nkpoint-1, ourjob, ourjob_disp)

      msize = PGEOM%nband*PINPT%nspin
      R = diagonal(msize)
     !IDXR = real((/1:msize/))
      do ii=1, msize
        IDXR(ii) = real(ii)
      enddo
      ov_cut = PINPT%band_order_overlap_cutoff
#ifdef MPI
      call MPI_BCAST(ETBA%V, size(ETBA%V), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
      if(flag_use_overlap) then
        call MPI_BCAST(ETBA%SV, size(ETBA%SV), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
      endif
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
   kp:do ik = 2, PKPTS%nkpoint
      m:do ii = 1, msize
           Mij_ = 0d0
        n:do jj = 1, msize
             if(flag_use_overlap) then
               Mij = real(abs( dot_product(ETBA%V(:,ii, ik  ), ETBA%SV(:,jj,ik-1)) ))
             else
               Mij = real(abs( dot_product(ETBA%V(:,ii, ik  ), ETBA%V(:,jj,ik-1)) ))
             endif
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
      do ik = 2, PKPTS%nkpoint
        R       = matprod (msize, msize, 'N', R, msize, msize, 'N', MM(:,:, ik))
        ! band swap by multiplying R_k by E(:,ik), where R_k = PI_k(MM(:,:,ik)) act as rotational matrix for E
        ETBA%IDX(:,ik) = int( matprod(msize, msize, 'N', R, msize, 'N', IDXR) )
        ETBA%E_ORD(:,ik)   = ETBA%E(ETBA%IDX(:,ik),ik)
        if(flag_order_weight) PWGHT%WT(:,ik)   = PWGHT%WT(ETBA%IDX(:,ik),ik)
      enddo

      if_main_then 
        ETBA%V_ORD(:,:,1) = ETBA%V(:,:,1)
        if(flag_use_overlap) then
          ETBA%SV_ORD(:,:,1) = ETBA%SV(:,:,1)
        endif
        do ik = 2, PKPTS%nkpoint 
          ETBA%V_ORD(:,:,ik) = ETBA%V(:,ETBA%IDX(:,ik),ik) 
          if(flag_use_overlap) then
            ETBA%SV_ORD(:,:,ik) = ETBA%SV(:,ETBA%IDX(:,ik),ik) 
          endif
        enddo
      if_main_end

      return
   endsubroutine

endmodule
