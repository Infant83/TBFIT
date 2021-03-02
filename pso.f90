#include "alias.inc"
subroutine pso_fit (PINPT, PPRAM, PKPTS, PWGHT, PGEOM, NN_TABLE, EDFT, iseed_, pso_miter)
    use mykind
    use mpi_setup
    use parameters, only: incar, params, kpoints, weight, poscar, hopping, energy
    use cost_function, only: get_dE, get_cost_function, get_dE12
    use print_io
    use random_mod
    implicit none
    type(incar   )                              :: PINPT ! only INCAR from first system is used
    type(params  )                              :: PPRAM ! asume each system shares same parameter
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
    real(kind=dp)                               :: cgbest, cgbest_plain, cost
    real(kind=dp), allocatable                  :: costs(:), costs_plain(:)
    real(kind=dp), allocatable                  :: costs_(:), costs_plain_(:)
    integer(kind=sp)                               i, min_id
    integer(kind=sp)                               imode, ldjac
    integer(kind=sp)                               min_loc(1), pso_miter
    logical                                        flag_order, flag_go
    real(kind=dp)                                  c1, c2, w, r1, r2
    integer(kind=sp)                               iseed, iseed_
    integer(kind=sp)                               mpierr
    real(kind=dp), external                     :: enorm
    real(kind=dp)                                  fnorm, fnorm_plain
    external                                       get_eig
    integer(kind=sp)                               info, maxfev
    real(kind=dp)                                  epsfcn, factor, xtol, ftol, gtol
    logical                                        flag_with_lmdif
    integer(kind=sp), allocatable               :: ourgroup(:), ourjob(:)
    integer(kind=sp)                               mygroup(0:nprocs-1), groupid
#ifdef MPI
    integer(kind=sp)                               mpistat(MPI_STATUS_SIZE)
#endif    

    flag_with_lmdif = PINPT%flag_pso_with_lmdif
    flag_order      = .FALSE.
    imode           = 13
    n_particles     = PPRAM%pso_nparticles
    nparam_free     = PPRAM%nparam_free
    c1              = PPRAM%pso_c1
    c2              = PPRAM%pso_c2
    w               = PPRAM%pso_w 
    iseed           = iseed_ ! set by default 
    r1              = 1d0
    r2              = 1d0
    call random_init(iseed)

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    !Note: n_particles will be solved in n_groups of cpu family.
    !      Each group will have n_member ~ nprocs/ngroup
    !      The n_member will be further parallized to solve eigenvalue problem in get_eig routine
    call get_npar_kpar() ; allocate(ourgroup(npar)); allocate(ourjob(npar))
    call mpi_job_distribution_group(npar, n_particles, ourgroup, mygroup, ourjob)
    call mpi_comm_anmeldung(COMM_KOREA, npar, mygroup)
    groupid = COMM_KOREA%color
    write(message,'(A,I0,A)')      '  JOB DISTRUBUTION for PSO particles over NPAR (NPAR=',npar,'):' ; write_msg
    do i = 1, npar
      write(message,'(A,3(I0,A))') '       -> Group_id(',i-1,'): ',ourgroup(i), &
                                   ' cpus asigned and ',ourjob(i),' jobs are distributed' ; write_msg
    enddo
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#else
    allocate(ourjob(1))
    ourjob = n_particles
    groupid= 0
#endif

    factor = 100.0D+00
    maxfev = 100 * ( PPRAM%nparam_free + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00; epsfcn = 0.000D+00

    if(allocated(PPRAM%pso_cost_history)) deallocate(PPRAM%pso_cost_history)
    allocate(PPRAM%pso_cost_history(pso_miter))
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
    allocate(costs_plain(n_particles        ))
    allocate(costs_(n_particles              ))
    allocate(costs_plain_(n_particles        ))
    costs = 0.d0 ; costs_plain = 0.d0 ; costs_ = 0.d0 ; costs_plain_ = 0.d0
    vel   = 0.d0 ; pos   = 0.d0 
    pbest = 0.d0 ; gbest = 0.d0
    cpbest = 0.d0 ;cgbest = 0.d0 ; cgbest_plain = 0.d0
                                
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
           if_main r1 = (random()*2d0 - 1d0) * PPRAM%pso_max_noise_amplitude
           if_main r2 = (random()*2d0 - 1d0) + 1d-10  ! sign
           if_main r2 = r2 / abs(r2)
           call MPI_BCAST(r1, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
           call MPI_BCAST(r2, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
           ! NOTE: all processors share random number over "iptcl" loop
           !       so that pos and vel update is to be mpi_comm_earth wide
           !       => up to this point, every processor knows vel and pos for all particles
#else
           r1 = (random()*2d0 - 1d0) * 5d0
           r2 = (random()*2d0 - 1d0) + 1d-10  ! sign
           r2 = r2 / abs(r2)
#endif
           if(iptcl .eq. 1) then ! use initial parameter without noise for particle-1 only
            r1 = 0d0
            r2 = 1d0
           endif
           pos(iptcl, iparam) = PPRAM%param( PPRAM%iparam_free(iparam) )*r2 + r1
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

    ! initial costs with initial parameters
    write(message,'(A)')" "; write_msg
    write(message,'(A)')"  Preparing initial costs" ; write_msg
    do iptcl = sum(ourjob(1:groupid)) + 1, sum(ourjob(1:groupid+1))
      PPRAM%param(PPRAM%iparam_free(:)) = pos(iptcl, :)
      if(.not. flag_with_lmdif) then
        call get_dE(fvec, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
        costs(iptcl) = enorm(ldjac, fvec)
        costs_plain(iptcl) = enorm(ldjac, fvec_plain)
      elseif(flag_with_lmdif) then
        call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
                 ftol, xtol, gtol, fnorm, fnorm_plain, maxfev, epsfcn, factor, info, .FALSE.)
        costs(iptcl) = fnorm
        costs_plain(iptcl) = fnorm_plain
      endif

    enddo

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    if(COMM_KOREA%flag_split) then
      if(COMM_KOREA%myid .eq. 0 .and. myid .ne. 0) then
        call MPI_SEND(costs(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1))), ourjob(groupid+1), &
                      MPI_REAL8, 0, 111, mpi_comm_earth, mpierr)
        call MPI_SEND(costs_plain(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1))), ourjob(groupid+1), &
                      MPI_REAL8, 0, 222, mpi_comm_earth, mpierr)
      elseif( myid .eq. 0 ) then
        do i = 1, npar - 1
            call MPI_RECV(costs(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                          MPI_REAL8, COMM_KOREA%group_main(i+1), 111, mpi_comm_earth, mpistat, mpierr)
            call MPI_RECV(costs_plain(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                          MPI_REAL8, COMM_KOREA%group_main(i+1), 222, mpi_comm_earth, mpistat, mpierr)
        enddo
      endif
    endif
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif
    ! report initial cost
    do iptcl = 1, n_particles
      write(message,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl); write_msg
    enddo

    if(myid .eq. 0) then
      if(.not. PPRAM%flag_fit_plain) then
        min_loc = minloc(costs(:))
      elseif(PPRAM%flag_fit_plain) then
        min_loc = minloc(costs_plain(:))
      endif
      min_id= min_loc(1)
      gbest = pos(min_loc(1),:)
      cgbest= costs(min_loc(1))
      cgbest_plain = costs_plain(min_loc(1))
      pbest = pos
      if(.not. PPRAM%flag_fit_plain) then
        cpbest= costs
      elseif(PPRAM%flag_fit_plain) then
        cpbest= costs_plain
      endif
    endif

#ifdef MPI
    ! share cost information with other members
    call MPI_BCAST(gbest       , size(gbest) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    call MPI_BCAST(cgbest      , 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    call MPI_BCAST(cgbest_plain, 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    call MPI_BCAST(pbest       , size(pbest) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    call MPI_BCAST(cpbest      , size(cpbest),  MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

    ! main PSO loop
    write(message,'(A)')" "; write_msg
    if(PPRAM%flag_fit_plain) then
        write(message,'(A)')"  Main PSO Loop start: minimizing COST w/o WEIGHT" ; write_msg
    elseif(.not. PPRAM%flag_fit_plain) then
        write(message,'(A)')"  Main PSO Loop start: minimizing COST " ; write_msg
    endif
 it:do iter = 1, pso_miter
      write(message,'(A)')" "; write_msg
      write(message,'(A, i5)')"  *Iter: ",iter ; write_msg

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
      costs = 0d0 ; costs_plain = 0d0
 ptcl:do iptcl = sum(ourjob(1:groupid)) + 1, sum(ourjob(1:groupid+1))
        PPRAM%param(PPRAM%iparam_free(:)) = pos(iptcl, :)
        if(.not. flag_with_lmdif) then
          call get_dE(fvec, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
          costs(iptcl) = enorm(ldjac, fvec)
          costs_plain(iptcl) = enorm(ldjac, fvec_plain)
        elseif(flag_with_lmdif) then
          call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
                   ftol, xtol, gtol, fnorm, fnorm_plain, maxfev, epsfcn, factor, info, .FALSE.)
          costs(iptcl) = fnorm
          costs_plain(iptcl) = fnorm_plain
        endif
      enddo ptcl

#ifdef MPI
      ! merge costs
      call MPI_BARRIER(mpi_comm_earth, mpierr)
      if(COMM_KOREA%flag_split) then
        if(COMM_KOREA%myid .eq. 0 .and. myid .ne. 0) then
          call MPI_SEND(costs(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1))), ourjob(groupid+1), &
                        MPI_REAL8, 0, 333, mpi_comm_earth, mpierr)
          call MPI_SEND(costs_plain(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1))), ourjob(groupid+1), &
                        MPI_REAL8, 0, 444, mpi_comm_earth, mpierr)
        elseif( myid .eq. 0 ) then
          do i = 1, npar - 1
              call MPI_RECV(costs(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                            MPI_REAL8, COMM_KOREA%group_main(i+1), 333, mpi_comm_earth, mpistat, mpierr)
              call MPI_RECV(costs_plain(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                            MPI_REAL8, COMM_KOREA%group_main(i+1), 444, mpi_comm_earth, mpistat, mpierr)
          enddo
        endif
      endif
      call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif
      ! report initial cost
      do iptcl = 1, n_particles
        write(message,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl); write_msg
      enddo

      ! update gbest, pbest
      if(myid .eq. 0) then
        if(.not. PPRAM%flag_fit_plain) then
          min_loc = minloc(costs(:))
          if( costs(min_loc(1)) .lt. cgbest) then
            cgbest = costs(min_loc(1))
            cgbest_plain = costs_plain(min_loc(1))
            gbest  = pos(min_loc(1),:)
            min_id = min_loc(1)
          endif
          do iptcl = 1, n_particles
            if( costs(iptcl) .lt. cpbest(iptcl) ) then
              cpbest(iptcl)   = costs(iptcl)
               pbest(iptcl,:) = pos(iptcl,:)
            endif
          enddo
        elseif(PPRAM%flag_fit_plain) then
          min_loc = minloc(costs_plain(:))
          if( costs_plain(min_loc(1)) .lt. cgbest_plain) then
            cgbest = costs(min_loc(1))
            cgbest_plain = costs_plain(min_loc(1))
            gbest  = pos(min_loc(1),:)
          endif
          do iptcl = 1, n_particles
            if( costs_plain(iptcl) .lt. cpbest(iptcl) ) then
              cpbest(iptcl)   = costs_plain(iptcl)
               pbest(iptcl,:) = pos(iptcl,:)
            endif
          enddo
        endif
      endif
#ifdef MPI
      ! share cost information with other members
      call MPI_BCAST(gbest       , size(gbest) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(cgbest      , 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(cgbest_plain, 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(pbest       , size(pbest) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(cpbest      , size(cpbest),  MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

      if(.not. PPRAM%flag_fit_plain) then
        PPRAM%pso_cost_history(iter) = cgbest
      elseif(PPRAM%flag_fit_plain) then
        PPRAM%pso_cost_history(iter) = cgbest_plain
      endif
      PPRAM%niter                  = iter

      write(message,'(A, i5, 2(A, F24.12), A,I0)')"  PSO RESULT: iter =",iter, &
                                                  ", GBest COST = ", cgbest, " ,COST(w/o weight) = ", cgbest_plain, &
                                                  ", Best Particle ID = ", min_id ; write_msg

    enddo it ! iter

    ! calculate with final parameters and save costs ETBA%dE
    PPRAM%param(PPRAM%iparam_free(:)) = gbest(:)
    if(.not. flag_with_lmdif) then
      ldjac = 1
      call get_dE(fvec, fvec_plain, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
    elseif(flag_with_lmdif) then
      call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm, fnorm_plain, maxfev, epsfcn, factor, info, .FALSE.)
      call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS)
    endif

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    COMM_KOREA%flag_split = .FALSE.
#endif
   
    return
endsubroutine
