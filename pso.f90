#include "alias.inc"
subroutine pso_fit (PINPT, PPRAM, PKPTS, PWGHT, PGEOM, NN_TABLE, EDFT)
    use mykind
    use mpi_setup
    use parameters, only: incar, params, kpoints, weight, poscar, hopping, energy
    use cost_function, only: get_dE, get_cost_function, get_dE12
    use print_io
    use random_mod
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(kpoints ), dimension(PINPT%nsystem)    :: PKPTS
    type(weight  ), dimension(PINPT%nsystem)    :: PWGHT
    type(poscar  ), dimension(PINPT%nsystem)    :: PGEOM
    type(hopping ), dimension(PINPT%nsystem)    :: NN_TABLE
    type(energy  ), dimension(PINPT%nsystem)    :: EDFT, ETBA
    logical                                        flag_stat
    integer(kind=sp)                               n_particles, nparam_free
    integer(kind=sp)                               iparam, iter, iptcl
    real(kind=dp), allocatable                  :: fvec(:), fvec_plain(:)
    real(kind=dp), allocatable                  :: vel(:,:), pos(:,:)
    real(kind=dp), allocatable                  :: pbest(:,:), gbest(:)
    real(kind=dp), allocatable                  :: cpbest(:)
    real(kind=dp)                               :: cgbest, cost
    real(kind=dp), allocatable                  :: costs(:)
    integer(kind=sp)                               imode, ldjac
    integer(kind=sp)                               min_loc(1)
    logical                                        flag_order, flag_go
    real(kind=dp)                                  c1, c2, w, r1, r2
    integer(kind=sp)                               iseed
    integer(kind=sp)                               mpierr
    real(kind=dp), external                     :: enorm
    real(kind=dp)                                  fnorm
    external                                       get_eig
    integer(kind=sp)                               info, maxfev
    real(kind=dp)                                  epsfcn, factor, xtol, ftol, gtol
    logical                                        flag_with_lmdif
    
    flag_with_lmdif = PINPT%flag_pso_with_lmdif

    iseed      = 123 ! set by default 
    call random_init(iseed)
    flag_order = .FALSE.
    imode      = 13
    n_particles= PPRAM%pso_nparticles
    nparam_free= PPRAM%nparam_free
    c1    = PPRAM%pso_c1
    c2    = PPRAM%pso_c2
    w     = PPRAM%pso_w 
    r1    = 1d0  ; r2    = 1d0

    factor = 100.0D+00
    maxfev = 100 * ( PPRAM%nparam_free + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00; epsfcn = 0.000D+00

    if(allocated(PPRAM%pso_cost_history)) deallocate(PPRAM%pso_cost_history)
    allocate(PPRAM%pso_cost_history(PINPT%miter))
    PPRAM%pso_cost_history = 0d0

    do i = 1, PINPT%nsystem
      call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i))
      if(allocated(ETBA(i)%dE)) deallocate(ETBA(i)%dE)
      allocate(ETBA(i)%dE( size(ETBA(i)%E(:,1)) , size(ETBA(i)%E(1,:)) ))
      ETBA(i)%E = 0d0 ; ETBA(i)%V = 0d0 ; ETBA(i)%SV = 0d0 ; ETBA(i)%dE=0d0
    enddo

    allocate(vel(  n_particles, nparam_free))
    allocate(pos(  n_particles, nparam_free))
    allocate(pbest(n_particles, nparam_free))
    allocate(gbest(              nparam_free))
    allocate(cpbest(n_particles             ))
    allocate(costs(n_particles              ))
    costs = 0.d0
    vel   = 0.d0 ; pos   = 0.d0 ; pbest = 0.d0 ; gbest = 0.d0
                                 cpbest = 0.d0 ;cgbest = 0.d0 
                                
    ldjac      = 0
    do i = 1, PINPT%nsystem
        ldjac = ldjac + PKPTS(i)%nkpoint * PGEOM(i)%nband * PINPT%nspin
    enddo
    if(ldjac .ne. 1) then
      allocate(fvec(ldjac))
      allocate(fvec_plain(ldjac))
    endif
    fvec = 0d0 ; fvec_plain = 0d0

    ! initialize parameters with random noise
    do iptcl = 1, n_particles
      do iparam = 1, nparam_free
         flag_go = .FALSE.
         do while (.not. flag_go)
#ifdef MPI
           if_main r1 = (random()*2d0 - 1d0) * 5d0
           call MPI_BCAST(r1, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#else
           r1 = (random()*2d0 - 1d0) * 5d0
#endif
           pos(iptcl, iparam) = PPRAM%param( PPRAM%iparam_free(iparam) ) + r1
           vel(iptcl, iparam) = r1
           if( pos(iptcl, iparam) .ge. PPRAM%param_const(2,iparam) .or. &
               pos(iptcl, iparam) .le. PPRAM%param_const(3,iparam) ) then
               flag_go = .FALSE.
           else
               flag_go = .TRUE. 
           endif

         enddo
      enddo
    enddo   

    ! prepare initial costs with initial parameters
    do iptcl = 1, n_particles
      PPRAM%param(PPRAM%iparam_free(:)) = pos(iptcl, :)
      if(.not. flag_with_lmdif) then
        call get_dE(fvec, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
!       call get_eig(NN_TABLE(1), PKPTS(1)%kpoint, PKPTS(1)%nkpoint, PINPT, PPRAM, &
!                    ETBA(1)%E, ETBA(1)%V, ETBA(1)%SV, PGEOM(1)%neig, PGEOM(1)%init_erange, PGEOM(1)%nband, &
!                    PINPT%flag_get_orbital, PINPT%flag_sparse, .FALSE., PINPT%flag_phase)
!       call get_cost_function(fvec, fvec_plain, ETBA(1), EDFT(1), PGEOM(1), &
!                              PINPT, PKPTS(1), PWGHT(1), ldjac, imode, flag_order)
        costs(iptcl) = enorm(ldjac, fvec)
      elseif(flag_with_lmdif) then
        call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
                 ftol, xtol, gtol, fnorm, maxfev, epsfcn, factor, info, .FALSE.)
        costs(iptcl) = fnorm
      endif
    enddo
    min_loc = minloc(costs(:))
    gbest = pos(min_loc(1),:)
    cgbest= costs(min_loc(1))
    pbest = pos
    cpbest= costs

    ! main loop
    do iter = 1, PINPT%miter

      ! update particle velocity and position
      do iptcl = 1, n_particles
        do iparam = 1, nparam_free
#ifdef MPI
          if_main r1 = random() 
          if_main r2 = random()
          call MPI_BCAST(r1, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
          call MPI_BCAST(r2, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#else
          r1 = random() ; r2 = random()
#endif
          vel(iptcl, iparam) = w * vel(iptcl, iparam) + c1 * r1 * (pbest(iptcl, iparam) - pos(iptcl, iparam)) &
                                                      + c2 * r2 * (gbest(iparam)        - pos(iptcl, iparam))
        enddo
      enddo
      pos = pos + vel

      ! get cost for each particle 
      costs = 0d0
      do iptcl = 1, n_particles
        PPRAM%param(PPRAM%iparam_free(:)) = pos(iptcl, :)
        if(.not. flag_with_lmdif) then
          call get_dE(fvec, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
         !call get_eig(NN_TABLE(1), PKPTS(1)%kpoint, PKPTS(1)%nkpoint, PINPT, PPRAM, &
         !             ETBA(1)%E, ETBA(1)%V, ETBA(1)%SV, PGEOM(1)%neig, PGEOM(1)%init_erange, PGEOM(1)%nband, &
         !             PINPT%flag_get_orbital, PINPT%flag_sparse, .FALSE., PINPT%flag_phase)
         !call get_cost_function(fvec, fvec_plain, ETBA(1), EDFT(1), PGEOM(1), &
         !                       PINPT, PKPTS(1), PWGHT(1), ldjac, imode, flag_order)
          
          costs(iptcl) = enorm(ldjac, fvec)
        elseif(flag_with_lmdif) then
          call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
                   ftol, xtol, gtol, fnorm, maxfev, epsfcn, factor, info, .FALSE.)
          costs(iptcl) = fnorm
        endif
      enddo

      ! update gbest, pbest
      min_loc = minloc(costs(:))
      if( costs(min_loc(1)) .lt. cgbest) then
        cgbest = costs(min_loc(1))
        gbest  = pos(min_loc(1),:)
      endif
      do iptcl = 1, n_particles
        if( costs(iptcl) .lt. cpbest(iptcl) ) then
          cpbest(iptcl)   = costs(iptcl)
           pbest(iptcl,:) = pos(iptcl,:)
        endif
      enddo

      PPRAM%pso_cost_history(iter) = cgbest
      PPRAM%niter                  = iter

      write(message,'(A, i5, A, F20.8)')" PSO STATUS: iter=",iter, " GBest COST = ", cgbest ; write_msg

    enddo

    ! calculate with final parameters and save costs ETBA%dE
    PPRAM%param(PPRAM%iparam_free(:)) = gbest(:)
    if(.not. flag_with_lmdif) then
      ldjac = 1
      call get_dE(fvec, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
     !call get_eig(NN_TABLE(1), PKPTS(1)%kpoint, PKPTS(1)%nkpoint, PINPT, PPRAM, &
     !             ETBA(1)%E, ETBA(1)%V, ETBA(1)%SV, PGEOM(1)%neig, PGEOM(1)%init_erange, PGEOM(1)%nband, &
     !             PINPT%flag_get_orbital, PINPT%flag_sparse, .FALSE., PINPT%flag_phase)
     !call get_cost_function(fvec, fvec_plain, ETBA(1), EDFT(1), PGEOM(1), &
     !                       PINPT, PKPTS(1), PWGHT(1), ldjac, imode, flag_order)
    elseif(flag_with_lmdif) then
      call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm, maxfev, epsfcn, factor, info, .FALSE.)
      call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
    endif

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

    return
endsubroutine
