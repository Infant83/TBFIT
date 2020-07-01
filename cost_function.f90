#include "alias.inc"
module cost_function
   use mpi_setup

contains 

  subroutine get_cost_function(fvec, ETBA_FIT, EDFT, neig, iband, nband, PINPT, nkpoint, PWGHT, ldjac, imode, flag_order)
    use parameters, only: weight, incar, energy
    use mpi_setup
    implicit none
    type (weight)          ::  PWGHT
    type (incar )          ::  PINPT
    type (energy)          ::  ETBA_FIT, EDFT
    integer*4                  i, ik, neig, nkpoint
    integer*4                  is, ie, ie_, iband, nband
    integer*4                  mpierr
    integer*4                  ldjac, imode
    real*8                     fvec(ldjac)
    real*8                     fvec_(ldjac)
    real*8                     OW
    real*8                     dE,         dD(3)
    real*8                     dE2(nband), dD2(3,nband)
    real*8                     cost_func(nband*PINPT%nspin)
    real*8                     enorm
    external                   enorm
    integer*4                  ourjob(nprocs), ourjob_disp(0:nprocs-1), sumk, sumk1, my_i
    logical                    flag_get_orbital, flag_fit_degeneracy
    logical                    flag_order
    real*8                     E_TBA(nband*PINPT%nspin,nkpoint)
    real*8                     E_DFT(neig*PINPT%ispin,nkpoint)
    complex*16                 V(neig*PINPT%ispin,nband*PINPT%nspin,nkpoint)
    ! imode : 1, if ldjac < nparam (total number of parameters) -> unusual cases. try to avoid.
    ! imode : 2, if ldjac > nparam -> usual cases
  
    flag_get_orbital = (.not. PWGHT%flag_weight_default_orb)
    flag_fit_degeneracy = PINPT%flag_fit_degeneracy

    OW = 0d0
    dD = 0d0
    dE = 0d0
    cost_func = 0d0
    if(flag_fit_degeneracy) dD2= 0d0
    if(imode .eq. 1) then
      if(.not. flag_get_orbital ) then
        call mpi_job_ourjob(nkpoint, ourjob)
        sumk = sum(ourjob(1:myid)) ; sumk1 = sum(ourjob(1:myid+1)) ; fvec_ = 0d0
        do ik= sumk + 1, sumk1
          my_i  = (ik - 1) * nband * PINPT%nspin
          do is = 1, PINPT%nspin
            ie_ = iband - 1 + (is-1) * neig
            if(flag_order) then
              dE2 =         (ETBA_FIT%E_ORD(1+(is-1)*nband:is*nband,ik) - EDFT%E_ORD(ie_+1:ie_+nband,ik))**2
            else
              dE2 =         (ETBA_FIT%E(1+(is-1)*nband:is*nband,ik) - EDFT%E(ie_+1:ie_+nband,ik))**2
            endif
            if(flag_fit_degeneracy) then ! would not work with flag_order 
              dD2(1,:)= (ETBA_FIT%D(1,1+(is-1)*nband:is*nband,ik) - EDFT%D(1,ie_+1:ie_+nband,ik))**2
              fvec_(my_i+1 +nband*(is-1):my_i+nband*is) = dE2(:)   *            PWGHT%WT(ie_+1:ie_+nband,ik) + &
                                                          dD2(1,:) * PWGHT%DEGENERACY_WT(ie_+1:ie_+nband,ik)
            else
              fvec_(my_i+1 +nband*(is-1):my_i+nband*is) = dE2(:)   *            PWGHT%WT(ie_+1:ie_+nband,ik)
            endif

          enddo
        enddo
#ifdef MPI
        call MPI_ALLREDUCE(fvec_, fvec, size(fvec), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
#else
        fvec = fvec_
#endif
      elseif(flag_get_orbital) then
! NOTE: this routine is not been parallelized yet. In the future it is highly recommended to parallelize it.
!       However, at this moment, with a quick test with parallization, it seems does not much improved with 
!       parallelized routines. Thereby I decided left it to be serial routine. One can check it later on.
!       The parallism can be achieved with the similar scheme in the above routine with (.not. flag_get_orbital) case. HJ Kim. 04.June.2020
        do ik=1,nkpoint
          do is = 1, PINPT%nspin
            do ie = 1, nband
              ie_ = ie + iband - 1 + (is-1) * neig
              if(flag_order) then
                OW  = sum( PWGHT%PENALTY_ORB(:,ie_,ik)*abs(ETBA_FIT%V_ORD(:,ie+(is-1)*nband,ik)) ) ! orbital weight
                dE  = (ETBA_FIT%E_ORD(ie+(is-1)*nband,ik) - EDFT%E_ORD(ie_,ik)) * PWGHT%WT(ie_,ik)
              else
                OW  = sum( PWGHT%PENALTY_ORB(:,ie_,ik)*abs(ETBA_FIT%V(:,ie+(is-1)*nband,ik)) ) ! orbital weight
                dE  = (ETBA_FIT%E(ie+(is-1)*nband,ik) - EDFT%E(ie_,ik)) * PWGHT%WT(ie_,ik)
              endif
              if(flag_fit_degeneracy) dD(1)  = (ETBA_FIT%D(1,ie+(is-1)*nband,ik) - EDFT%D(1,ie_,ik)) * PWGHT%DEGENERACY_WT(ie_,ik) ! degeneracy 

              cost_func(ie+(is-1)*nband) = dE + OW + dD(1)

            enddo
          enddo
          fvec(ik) = enorm ( nband*PINPT%nspin, cost_func )
        enddo
      endif
  
    elseif(imode .eq. 2) then
      if(.not. flag_get_orbital) then
        do ik=1,nkpoint
          do is = 1, PINPT%nspin
            do ie = 1, nband
              ie_ = ie + iband - 1 + (is-1) * neig
              if(flag_order) then
                dE  = (ETBA_FIT%E_ORD(ie+(is-1)*nband,ik) - EDFT%E_ORD(ie_,ik))    * PWGHT%WT(ie_,ik)
              else
                dE  = (ETBA_FIT%E(ie+(is-1)*nband,ik) - EDFT%E(ie_,ik))    * PWGHT%WT(ie_,ik)
              endif
              if(flag_fit_degeneracy) dD(1:3)  = (ETBA_FIT%D(1:3,ie+(is-1)*nband,ik) - EDFT%D(1:3,ie_,ik))       * PWGHT%DEGENERACY_WT(ie_,ik)

              cost_func(ie+(is-1)*nband) = abs(dE)  + sum(abs(dD(1:3))) !!/3d0
            enddo
          enddo
          fvec(ik) = enorm ( nband*PINPT%nspin, cost_func )
        enddo

      elseif(flag_get_orbital) then
#ifdef MPI
        if(flag_order) then
          call MPI_BCAST(ETBA_FIT%V_ORD, size(ETBA_FIT%V_ORD), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
        else
          call MPI_BCAST(ETBA_FIT%V, size(ETBA_FIT%V), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
        endif
#endif
        do ik=1,nkpoint
          do is = 1, PINPT%nspin
            do ie = 1, nband
              ie_ = ie + iband - 1 + (is-1) * neig
              if(flag_order) then
                OW  = sum( PWGHT%PENALTY_ORB(:,ie_,ik)*abs(ETBA_FIT%V_ORD(:,ie+(is-1)*nband,ik)) )
                dE  = (ETBA_FIT%E_ORD(ie+(is-1)*nband,ik) - EDFT%E_ORD(ie_,ik)) * PWGHT%WT(ie_,ik)
              else
                OW  = sum( PWGHT%PENALTY_ORB(:,ie_,ik)*abs(ETBA_FIT%V(:,ie+(is-1)*nband,ik)) )
                dE  = (ETBA_FIT%E(ie+(is-1)*nband,ik) - EDFT%E(ie_,ik)) * PWGHT%WT(ie_,ik)
              endif
              if(flag_fit_degeneracy) dD(1) = (ETBA_FIT%D(1,ie+(is-1)*nband,ik) - EDFT%D(1,ie_,ik)) * PWGHT%DEGENERACY_WT(ie_,ik)

              cost_func(ie+(is-1)*nband) = dE + OW + dD(1)

            enddo
          enddo
          fvec(ik) = enorm ( nband*PINPT%nspin, cost_func )
        enddo
      endif
    endif

    return

  endsubroutine


endmodule
