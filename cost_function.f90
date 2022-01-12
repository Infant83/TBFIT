#include "alias.inc"
module cost_function
   use mpi_setup

contains 
  subroutine get_dE(fvec, fvec_plain, fvec_orb, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA_FIT, PWGHT, PGEOM, PKPTS, flag_fit_orbital)
    use parameters, only: weight, incar, energy, params, poscar, hopping, kpoints
    use projected_band
    use reorder_band
    implicit none
    type(incar)                            :: PINPT
    type(params)                           :: PPRAM
    type(hopping),dimension(PINPT%nsystem) :: NN_TABLE
    type(poscar) ,dimension(PINPT%nsystem) :: PGEOM
    type(weight) ,dimension(PINPT%nsystem) :: PWGHT
    type(energy) ,dimension(PINPT%nsystem) :: ETBA_FIT, EDFT
    type(kpoints),dimension(PINPT%nsystem) :: PKPTS
    integer*4                                 mpierr
    integer*4                                 i
    integer*4                                 ldjac, my_ldjac
    integer*4                                 imode
    logical                                   flag_order
    logical                                   flag_get_orbital
    logical                                   flag_order_weight
    logical                                   flag_fit_orbital
    real*8                                    fvec(ldjac)
    real*8                                    fvec_plain(ldjac)
    real*8                                    fvec_orb(ldjac)
    integer*4                                 ildjac, fldjac

    flag_order_weight   = .false. ! experimental feature
    flag_order          = PINPT%flag_get_band_order .and. (.not. PINPT%flag_get_band_order_print_only)
    ildjac              = 0
    fldjac              = 0

    do i = 1, PINPT%nsystem
      flag_get_orbital = (PWGHT(i)%flag_weight_orb .or. flag_order .or. flag_fit_orbital)
      if(imode .eq. 1 .or. imode .eq. 11 .or. imode .eq. 13) then
        my_ldjac = PKPTS(i)%nkpoint * PGEOM(i)%nband * PINPT%nspin
      elseif(imode .eq. 2 .or. imode .eq. 12) then
        my_ldjac = PKPTS(i)%nkpoint
      endif
      ildjac = fldjac + 1
      fldjac = ildjac + my_ldjac - 1
        
      ! Evaluate the function at the starting point and calculate its norm.
      call get_eig(NN_TABLE(i), PKPTS(i)%kpoint, PKPTS(i)%nkpoint, PINPT, PPRAM, &
                   ETBA_FIT(i)%E, ETBA_FIT(i)%V, ETBA_FIT(i)%SV, &
                   PGEOM(i)%neig, PGEOM(i)%init_erange, PGEOM(i)%nband, &
                   flag_get_orbital, .false., .false., PINPT%flag_phase)
      ! Evaluate degeneracy information for TBA band : after calling get_eig
      call get_degeneracy(ETBA_FIT(i), PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint, PINPT)

      ! Evaluate ordered band for TBA band
      if(flag_order) then 
        call get_ordered_band(ETBA_FIT(i), PKPTS(i), PGEOM(i), PWGHT(i), PINPT, flag_order_weight, PPRAM%flag_use_overlap)
      endif

      ! Evaluate orbital projection for TBA band
      if(flag_fit_orbital) then
        call get_orbital_projection(ETBA_FIT(i), PKPTS(i), PINPT, PGEOM(i), PPRAM%flag_use_overlap)
      endif

      call get_cost_function(fvec(ildjac:fldjac), fvec_plain(ildjac:fldjac), fvec_orb(ildjac:fldjac), &
                             ETBA_FIT(i), EDFT(i), PGEOM(i), PINPT, PKPTS(i), PWGHT(i), my_ldjac, imode, flag_order, flag_fit_orbital)
    enddo

    return
  endsubroutine

  subroutine get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA_FIT, PWGHT, PGEOM, PKPTS, flag_fit_orbital)
    use parameters, only: weight, incar, energy, params, poscar, hopping, kpoints
    use projected_band
    use reorder_band
    implicit none
    type(incar)                            :: PINPT
    type(params)                           :: PPRAM
    type(hopping),dimension(PINPT%nsystem) :: NN_TABLE
    type(poscar) ,dimension(PINPT%nsystem) :: PGEOM
    type(weight) ,dimension(PINPT%nsystem) :: PWGHT
    type(energy) ,dimension(PINPT%nsystem) :: ETBA_FIT, EDFT
    type(kpoints),dimension(PINPT%nsystem) :: PKPTS
    integer*4                                 mpierr
    integer*4                                 i
    integer*4                                 ldjac, my_ldjac
    integer*4                                 imode
    logical                                   flag_order
    logical                                   flag_get_orbital
    logical                                   flag_order_weight
    logical                                   flag_fit_orbital
    real*8                                    fvec(1)
    real*8                                    fvec_plain(1)
    real*8                                    fvec_orb(1)
    integer*4                                 ildjac, fldjac

    flag_order_weight   = .false. ! experimental feature
    flag_order          = PINPT%flag_get_band_order .and. (.not. PINPT%flag_get_band_order_print_only)
    imode = 12

    do i = 1, PINPT%nsystem
      flag_get_orbital = (PWGHT(i)%flag_weight_orb .or. flag_order .or. flag_fit_orbital)

      ! Evaluate the function at the starting point and calculate its norm.
      call get_eig(NN_TABLE(i), PKPTS(i)%kpoint, PKPTS(i)%nkpoint, PINPT, PPRAM, &
                   ETBA_FIT(i)%E, ETBA_FIT(i)%V, ETBA_FIT(i)%SV, &
                   PGEOM(i)%neig, PGEOM(i)%init_erange, PGEOM(i)%nband, &
                   flag_get_orbital, .false., .false., PINPT%flag_phase)

      ! Evaluate degeneracy information for TBA band : after calling get_eig
      call get_degeneracy(ETBA_FIT(i), PGEOM(i)%nband*PINPT%nspin, PKPTS(i)%nkpoint, PINPT)

      ! Evaluate ordered band for TBA band
      if(flag_order) then
        call get_ordered_band(ETBA_FIT(i), PKPTS(i), PGEOM(i), PWGHT(i), PINPT, flag_order_weight, PPRAM%flag_use_overlap)
      endif

      ! Evaluate orbital projection for TBA band
      if(flag_fit_orbital) then
        call get_orbital_projection(ETBA_FIT(i), PKPTS(i), PINPT, PGEOM(i), PPRAM%flag_use_overlap)
      endif
      call get_cost_function(fvec(1), fvec_plain(1), fvec_orb(1), ETBA_FIT(i), EDFT(i), PGEOM(i), PINPT, PKPTS(i), PWGHT(i), 1, imode, flag_order, flag_fit_orbital)
    enddo

    return
  endsubroutine

  subroutine get_cost_function(fvec, fvec_plain, fvec_orb, ETBA_FIT, EDFT, PGEOM, PINPT, PKPTS, PWGHT, ldjac, imode, flag_order, flag_fit_orbital)
    use parameters, only: weight, incar, energy, pi2, poscar, kpoints
    use do_math,    only: fgauss 
    use mpi_setup
    implicit none
    type (weight)          ::  PWGHT
    type (incar )          ::  PINPT
    type (poscar)          ::  PGEOM
    type (kpoints)         ::  PKPTS
    type (energy)          ::  ETBA_FIT, EDFT
    integer*4                  i, ik
    integer*4                  is, ie, ie_, iband
    integer*4                  je, je_
    integer*4                  my_ik, sizebuff
    integer*4                  mpierr
    integer*4                  ldjac, imode
    real*8                     fvec(ldjac)
    real*8                     fvec_(ldjac)
    real*8                     fvec_plain(ldjac)
    real*8                     fvec_plain_(ldjac)
    real*8                     fvec_orb(ldjac)
    real*8                     fvec_orb_(ldjac)
    real*8                     OW, dD(3)
    real*8                     dE, dE_plain, dE_smooth, dORB_smooth
    real*8                     cost_func(PGEOM%nband*PINPT%nspin)
    real*8                     cost_func_plain(PGEOM%nband*PINPT%nspin)
    real*8                     sigma, maxgauss, max_gauss
    real*8                     enorm
    external                   enorm
    logical                    flag_weight_orbital, flag_fit_degeneracy
    logical                    flag_order
    real*8                     E_TBA(PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
    real*8                     dE_TBA(PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
    real*8                     dORB_TBA(PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
    real*8                     E_DFT(PGEOM%neig*PINPT%ispin,PKPTS%nkpoint)
    complex*16, allocatable :: myV(:,:,:)
    real*8,     allocatable :: myORB_TBA(:,:,:), myORB_DFT(:,:,:)
    real*8,     allocatable :: orb_TBA(:), orb_DFT(:,:)
    real*8,     allocatable :: dE_gauss(:) , dE_ref(:), dORB(:), dORB_ref(:)
    logical                    flag_fit_orbital
    integer*4                  ie_cutoff
   !integer*4                  ourjob(nprocs), ourjob_disp(0:nprocs-1), my_i
    integer*4                  ncpu, id
    integer*4,  allocatable :: ourjob(:), ourjob_disp(:), my_i
    
    ! imode : 1, if ldjac < nparam (total number of parameters) -> unusual cases. try to avoid.
    !               ldjac = nkpoint * nband * PINPT%nspin
    ! imode : 2, if ldjac > nparam -> usual cases
    !               ldjac = nkpoint

    ! Important NOTE:
    ! imode = 1       ldjac != 1 : return fvec, fvec_plain (ldjac = nkpoint * nband)
    ! imode = 2       ldjac != 1 : return fvec, fvec_plain (ldjac = nkpoint)
    ! imode = 11 with ldjac != 1 : return ETBA%dE, fvec, fvec_plain (ldjac = nkpoint * nband)
    ! imode = 12 with ldjac != 1 : return ETBA%dE, fvec, fvec_plain (ldjac = nkpoint)
    ! imode = 13 with ldjac  = 1 : return ETBA%dE 
    ! imode = 13 with ldjac != 1 : return          fvec, fvec_plain (ldjac = nkpoint * nband)

    flag_fit_degeneracy = PINPT%flag_fit_degeneracy
    flag_weight_orbital = PWGHT%flag_weight_orb

    dE        = 0d0 ;  sigma            = PINPT%orbital_fit_smearing   ! default 
    cost_func = 0d0 ;  cost_func_plain  = 0d0 
    if(ldjac .ne. 1) then
      fvec      = 0d0 ;  fvec_plain       = 0d0  ; fvec_orb  = 0d0
      fvec_     = 0d0 ;  fvec_plain_      = 0d0  ; fvec_orb_ = 0d0
    endif
    maxgauss  = 1d0/(sigma*sqrt(2d0*pi2)) *real(PGEOM%neig)
    max_gauss = fgauss(sigma, 0d0)
    iband     = PGEOM%init_erange
    ie_cutoff = PWGHT%ie_cutoff
    if(flag_fit_degeneracy) dD= 0d0
    if(flag_weight_orbital) OW= 0d0
    if(imode .gt. 10) dE_TBA = 0d0
    if(imode .gt. 10) dORB_TBA = 0d0

    if(COMM_KOREA%flag_split) then
      ncpu = COMM_KOREA%nprocs
      id   = COMM_KOREA%myid
    else
      ncpu = nprocs
      id   = myid
    endif
    allocate(ourjob(ncpu))
    allocate( ourjob_disp(0:ncpu-1))

    call mpi_job_distribution_chain(PKPTS%nkpoint, ncpu, ourjob, ourjob_disp)
    if(flag_order) then
      E_TBA = ETBA_FIT%E_ORD
      E_DFT = EDFT%E_ORD
      if(flag_weight_orbital) then
        sizebuff = PGEOM%neig*PINPT%ispinor*PGEOM%nband*PINPT%nspin
        allocate(myV(PGEOM%neig*PINPT%ispinor, PGEOM%nband*PINPT%nspin, ourjob(id+1)))
#ifdef MPI
        ! Distribute wavefunction data to neighboring nodes
        if(COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%V_ORD, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_COMPLEX16, myV, ourjob(id+1)*sizebuff, &
                            MPI_COMPLEX16, 0, COMM_KOREA%mpi_comm, mpierr)
        elseif(.not. COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%V_ORD, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_COMPLEX16, myV, ourjob(id+1)*sizebuff, &
                            MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
        endif
#else
        myV = ETBA_FIT%V_ORD
#endif
      endif

      if(flag_fit_orbital) then
        sizebuff = PINPT%lmmax * PGEOM%nband*PINPT%nspin
        allocate(dE_gauss(PGEOM%nband))
        allocate(dE_ref(  PGEOM%nband))
        allocate(dORB_ref(PGEOM%nband))
        allocate(myORB_TBA(PINPT%lmmax, PGEOM%nband*PINPT%nspin, ourjob(id+1)))
        allocate(myORB_DFT(PINPT%lmmax, PGEOM%nband*PINPT%nspin, ourjob(id+1))) 
        allocate(orb_TBA(PINPT%lmmax))
        allocate(orb_DFT(PINPT%lmmax, PGEOM%nband))
        ! note: the size of second column, nband*PINPT%nspin is differ from read_energy routine
        !       but should work anyway, since ERANGE or EWINDOW tag should not be applied if
        !       fitting is requested by TBFIT tag of INCAR
#ifdef MPI
        ! Distribute orbital info to neighboring nodes
        if(COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%ORB, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_TBA, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, COMM_KOREA%mpi_comm, mpierr)
          call MPI_SCATTERV(EDFT%ORB    , ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_DFT, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, COMM_KOREA%mpi_comm, mpierr)
        elseif(.not. COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%ORB, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_TBA, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, mpi_comm_earth, mpierr)
          call MPI_SCATTERV(EDFT%ORB    , ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_DFT, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, mpi_comm_earth, mpierr)
        endif
#else
        myORB_TBA = ETBA_FIT%ORB
        myORB_DFT = EDFT%ORB
#endif
      endif

    elseif(.not. flag_order) then
      E_TBA = ETBA_FIT%E
      E_DFT = EDFT%E
      if(flag_weight_orbital) then
        sizebuff = PGEOM%neig*PINPT%ispinor*PGEOM%nband*PINPT%nspin
        allocate(myV(PGEOM%neig*PINPT%ispinor, PGEOM%nband*PINPT%nspin, ourjob(id+1)))
#ifdef MPI
        ! Distribute wavefunction data to neighboring nodes
        if(COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%V, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_COMPLEX16, myV, ourjob(id+1)*sizebuff, &
                            MPI_COMPLEX16, 0, COMM_KOREA%mpi_comm, mpierr)
        elseif(.not. COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%V, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_COMPLEX16, myV, ourjob(id+1)*sizebuff, &
                            MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
        endif

#else    
        myV = ETBA_FIT%V
#endif
      endif

      if(flag_fit_orbital) then
        sizebuff = PINPT%lmmax * PGEOM%nband*PINPT%nspin
        allocate(dE_gauss(PGEOM%nband))
        allocate(dE_ref(  PGEOM%nband))
        allocate(dORB_ref(PGEOM%nband))
        allocate(dORB(PGEOM%nband))
        allocate(myORB_TBA(PINPT%lmmax, PGEOM%nband*PINPT%nspin, ourjob(id+1)))
        allocate(myORB_DFT(PINPT%lmmax, PGEOM%nband*PINPT%nspin, ourjob(id+1))) 
        allocate(orb_TBA(PINPT%lmmax))
        allocate(orb_DFT(PINPT%lmmax, PGEOM%nband))
        ! note: the size of second column, nband*PINPT%nspin is differ from read_energy routine
        !       but should work anyway, since ERANGE or EWINDOW tag should not be applied if
        !       fitting is requested by TBFIT tag of INCAR
#ifdef MPI
        ! Distribute orbital info to neighboring nodes
        if(COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%ORB, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_TBA, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, COMM_KOREA%mpi_comm, mpierr)
          call MPI_SCATTERV(EDFT%ORB    , ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_DFT, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, COMM_KOREA%mpi_comm, mpierr)
        elseif(.not. COMM_KOREA%flag_split) then
          call MPI_SCATTERV(ETBA_FIT%ORB, ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_TBA, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, mpi_comm_earth, mpierr)
          call MPI_SCATTERV(EDFT%ORB    , ourjob*sizebuff, ourjob_disp*sizebuff, &
                            MPI_REAL8, myORB_DFT, ourjob(id+1)*sizebuff, &
                            MPI_REAL8, 0, mpi_comm_earth, mpierr)
        endif
#else
        myORB_TBA = ETBA_FIT%ORB
        myORB_DFT = EDFT%ORB
#endif
      endif
    endif

!#### IMODE = 1 : if nkpoint is less than number of free parameters to be fitted, this routine works
    if(imode .eq. 1 .or. imode .eq. 11) then

      do ik = sum(ourjob(1:id))+1, sum(ourjob(1:id+1))
        my_i = (ik - 1) * PGEOM%nband * PINPT%nspin
        my_ik= ik - sum(ourjob(1:id))
        do is = 1, PINPT%nspin
          do ie = 1, PGEOM%nband
            ie_ = ie + iband - 1 + (is-1)*PGEOM%neig

            dE_plain  = E_TBA(ie+(is-1)*PGEOM%nband,ik) - E_DFT(ie_,ik)
            dE        = dE_plain * PWGHT%WT(ie_,ik)
            if(flag_fit_degeneracy) dD= (ETBA_FIT%D(:,ie+(is-1)*PGEOM%nband,ik) - EDFT%D(:,ie_,ik)) * PWGHT%DEGENERACY_WT(ie_,ik)
            if(flag_weight_orbital) OW = sum( PWGHT%PENALTY_ORB(:,ie_,ik)*abs(myV(:,ie+(is-1)*PGEOM%nband,my_ik)) )
            if(imode .eq. 11) then 
              dE_TBA(ie+(is-1)*PGEOM%nband,ik) = dE_plain
            else
              fvec(my_i + ie+(is-1)*PGEOM%nband) = abs(dE) + (sum(abs(dD(:)))) + OW
              if(ie_cutoff .gt. 0) then
                if(ie .le. ie_cutoff) then
                  fvec_plain(my_i + ie+(is-1)*PGEOM%nband) = abs(dE_plain)
                endif
              else
                fvec_plain(my_i + ie+(is-1)*PGEOM%nband) = abs(dE_plain)
              endif
            endif
          enddo ! ie
        enddo ! is 
      enddo ! ik

!#### IMODE = 2 ; ldjac > nparam -> usual cases
    elseif(imode .eq. 2 .or. imode .eq. 12) then
      do ik = sum(ourjob(1:id))+1, sum(ourjob(1:id+1))
        my_ik = ik - sum(ourjob(1:id))
        do is = 1, PINPT%nspin

          do ie = 1, PGEOM%nband
            ie_      = ie + iband - 1 + (is-1) * PGEOM%neig
            dE_plain = E_TBA(ie+(is-1)*PGEOM%nband,ik) - E_DFT(ie_,ik)
            
            if(flag_fit_degeneracy) dD       = (ETBA_FIT%D(:,ie+(is-1)*PGEOM%nband,ik) - EDFT%D(:,ie_,ik)) * PWGHT%DEGENERACY_WT(ie_,ik)
            if(flag_weight_orbital) OW       = sum( PWGHT%PENALTY_ORB(:,ie_,ik)*abs(myV(:,ie+(is-1)*PGEOM%nband,my_ik)) )
            dE       = abs(dE_plain * PWGHT%WT(ie_,ik))
            if(flag_fit_orbital) then
              orb_TBA    = myORB_TBA(:,ie+(is-1)*PGEOM%nband,my_ik)
              orb_DFT    = myORB_DFT(:,1+(is-1)*PGEOM%nband:PGEOM%nband+(is-1)*PGEOM%nband,my_ik)
              dORB       = fdORB(orb_TBA, orb_DFT, PINPT%lmmax, PGEOM%nband)
              dORB_ref   = fdORB(orb_DFT(:,ie), orb_DFT, PINPT%lmmax, PGEOM%nband)

            endif
            if(imode .eq. 12) then
              dE_TBA(ie+(is-1)*PGEOM%nband,ik) = dE_plain
              if(flag_fit_orbital) dORB_TBA(ie+(is-1)*PGEOM%nband,ik)= sum(abs(dORB - dORB_ref))
            else
              cost_func(ie+(is-1)*PGEOM%nband)       = dE + sum(abs(dD)) + OW
              if(ie_cutoff .gt. 0) then
                if(ie .le. ie_cutoff) then
                  cost_func_plain(ie+(is-1)*PGEOM%nband) = abs(dE_plain)
                endif
              else
                cost_func_plain(ie+(is-1)*PGEOM%nband) = abs(dE_plain)
              endif
            endif
          enddo

        enddo ! is
        if(imode .eq. 2) then
          fvec(ik) = enorm(PGEOM%nband*PINPT%nspin, cost_func)
          fvec_plain(ik) = enorm(PGEOM%nband*PINPT%nspin, cost_func_plain)
        endif
      enddo ! ik

    elseif(imode .eq. 13) then
      do ik = sum(ourjob(1:id))+1, sum(ourjob(1:id+1))
        my_i = (ik - 1) * PGEOM%nband * PINPT%nspin
        my_ik= ik - sum(ourjob(1:id))
        do is = 1, PINPT%nspin
          do ie = 1, PGEOM%nband
            ie_      = ie + iband - 1 + (is-1) * PGEOM%neig
            dE_plain = E_TBA(ie+(is-1)*PGEOM%nband,ik) - E_DFT(ie_,ik)
            dE       = dE_plain * PWGHT%WT(ie_,ik)
            dE_TBA(ie+(is-1)*PGEOM%nband,ik) = dE_plain

            if(flag_fit_orbital) then
              ! NOTE: what is the best way to make the orbital different to be a "smooth" function?
              ! orbital    :           s          px          py          pz         dz2         dx2         dxy         dxz         dyz"
              orb_TBA    = myORB_TBA(:,ie+(is-1)*PGEOM%nband,my_ik)
              orb_DFT    = myORB_DFT(:,1+(is-1)*PGEOM%nband:PGEOM%nband+(is-1)*PGEOM%nband,my_ik)
              dORB       = fdORB(orb_TBA, orb_DFT, PINPT%lmmax, PGEOM%nband)
              dORB_ref   = fdORB(orb_DFT(:,ie), orb_DFT, PINPT%lmmax, PGEOM%nband)
              dORB_smooth= sum(abs(dORB - dORB_ref)) * PWGHT%WT(ie_,ik)
            endif

            if(ldjac .ne. 1) then
              if(ie_cutoff .gt. 0) then
                if(ie .le. ie_cutoff) then
                  fvec(my_i + ie+(is-1)*PGEOM%nband) = abs(dE) 
                  fvec_plain(my_i + ie+(is-1)*PGEOM%nband) = abs(dE_plain)
                  if(flag_fit_orbital) then
                    fvec_orb(my_i + ie+(is-1)*PGEOM%nband) = dORB_smooth
                  endif
                endif
              else
                fvec(my_i + ie+(is-1)*PGEOM%nband) = abs(dE) 
                fvec_plain(my_i + ie+(is-1)*PGEOM%nband) = abs(dE_plain)
                if(flag_fit_orbital) then
                    fvec_orb(my_i + ie+(is-1)*PGEOM%nband) = dORB_smooth
                endif
              endif
            endif
          enddo ! ie
        enddo ! is
      enddo ! ik
    endif ! IMODE ?

    if(imode .lt. 10) then
#ifdef MPI
      if(COMM_KOREA%flag_split) then
        call MPI_ALLREDUCE(fvec, fvec_, size(fvec), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
        call MPI_ALLREDUCE(fvec_plain, fvec_plain_, size(fvec_plain), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
      elseif(.not. COMM_KOREA%flag_split) then
        call MPI_ALLREDUCE(fvec, fvec_, size(fvec), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
        call MPI_ALLREDUCE(fvec_plain, fvec_plain_, size(fvec_plain), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
      endif
      fvec = fvec_
      fvec_plain = fvec_plain_
#endif

    elseif(imode .gt. 10 .and. imode .le. 12) then
#ifdef MPI
      if(COMM_KOREA%flag_split) then
        call MPI_ALLREDUCE(dE_TBA, ETBA_FIT%dE, size(dE_TBA), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
        if(flag_fit_orbital) then
          call MPI_ALLREDUCE(dORB_TBA, ETBA_FIT%dORB, size(dORB_TBA), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
        endif
        if(ldjac .ne. 1) then
          call MPI_ALLREDUCE(fvec, fvec_, size(fvec), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
          fvec = fvec_
          call MPI_ALLREDUCE(fvec_plain, fvec_plain_, size(fvec_plain), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
          fvec_plain = fvec_plain_
          if(flag_fit_orbital) then
            call MPI_ALLREDUCE(fvec_orb, fvec_orb_, size(fvec_orb), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
            fvec_orb = fvec_orb_
          endif
        endif
      elseif(.not. COMM_KOREA%flag_split) then
        call MPI_ALLREDUCE(dE_TBA, ETBA_FIT%dE, size(dE_TBA), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
        if(flag_fit_orbital) then
          call MPI_ALLREDUCE(dORB_TBA, ETBA_FIT%dORB, size(dORB_TBA), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
        endif
        if(ldjac .ne. 1) then
          call MPI_ALLREDUCE(fvec, fvec_, size(fvec), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
          fvec = fvec_
          call MPI_ALLREDUCE(fvec_plain, fvec_plain_, size(fvec_plain), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
          fvec_plain = fvec_plain_
          if(flag_fit_orbital) then
            call MPI_ALLREDUCE(fvec_orb, fvec_orb_, size(fvec_orb), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
            fvec_orb = fvec_orb_
          endif
        endif
      endif
#else
      ETBA_FIT%dE = dE_TBA
      if(flag_fit_orbital) then
        ETBA_FIT%dORB = dORB_TBA
      endif
#endif
    elseif(imode .ge. 13) then
#ifdef MPI
      if(COMM_KOREA%flag_split) then
        if(ldjac .ne. 1) then
          call MPI_ALLREDUCE(fvec, fvec_, size(fvec), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
          fvec = fvec_
          call MPI_ALLREDUCE(fvec_plain, fvec_plain_, size(fvec_plain), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
          fvec_plain = fvec_plain_
          if(flag_fit_orbital) then
            call MPI_ALLREDUCE(fvec_orb, fvec_orb_, size(fvec_orb), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
            fvec_orb = fvec_orb_
          endif
        elseif(ldjac .eq. 1) then
          call MPI_ALLREDUCE(dE_TBA, ETBA_FIT%dE, size(dE_TBA), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
        endif
      elseif(.not. COMM_KOREA%flag_split) then
        if(ldjac .ne. 1) then
          call MPI_ALLREDUCE(fvec, fvec_, size(fvec), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
          fvec = fvec_
          call MPI_ALLREDUCE(fvec_plain, fvec_plain_, size(fvec_plain), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
          fvec_plain = fvec_plain_
          if(flag_fit_orbital) then
            call MPI_ALLREDUCE(fvec_orb, fvec_orb_, size(fvec_orb), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
            fvec_orb = fvec_orb_
          endif
        elseif(ldjac .eq. 1) then
          call MPI_ALLREDUCE(dE_TBA, ETBA_FIT%dE, size(dE_TBA), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
        endif
      endif
#else
      if(ldjac .eq. 1) then
        ETBA_FIT%dE = dE_TBA
      endif
#endif
    endif

    return
  endsubroutine

  function fdORB(myORB_TBA, myORB_DFT, lmmax, nband)
    use do_math,    only: fgauss
    implicit none
    integer*4    ie
    integer*4    lmmax
    integer*4    nband
    real*8       myORB_TBA(lmmax)
    real*8       myORB_DFT(lmmax, nband)
    real*8       fdORB(nband)
    external     enorm
    real*8       enorm

    do ie = 1, nband
      fdORB(ie) = dot_product(myORB_TBA, myORB_DFT(:,ie))/(sqrt(dot_product(myORB_TBA, myORB_TBA)) * sqrt(dot_product(myORB_DFT(:,ie),myORB_DFT(:,ie))))
    enddo

    return
  endfunction


! function fdORB2(E_DFT, myORB_TBA, myORB_DFT, lmmax, nband)
!   use do_math,    only: fgauss
!   implicit none
!   integer*4    ie
!   integer*4    lmmax
!   integer*4    nband
!   real*8       myORB_TBA(lmmax)
!   real*8       myORB_DFT(lmmax, nband)
!   real*8       fdORB(nband)
!   real*8       E_DFT(nband)
!   external     enorm
!   real*8       enorm

!   do ie = 1, nband
!     fdORB(ie) = dot_product(myORB_TBA, myORB_DFT(:,ie))/(sqrt(dot_product(myORB_TBA, myORB_TBA)) * sqrt(dot_product(myORB_DFT(:,ie),myORB_DFT(:,ie))))
!   enddo

!   return
! endfunction

endmodule
