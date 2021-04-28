subroutine pso_fit_best (PINPT, PPRAM, PKPTS, PWGHT, PGEOM, NN_TABLE, EDFT, iseed_, pso_miter)
    use mykind
    use mpi_setup
    use parameters, only: incar, params, kpoints, weight, poscar, hopping, energy
    use cost_function, only: get_dE, get_cost_function, get_dE12
    use print_io
    use random_mod
    use sorting
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
    real(kind=dp), allocatable                  :: fvec(:), fvec_plain(:), fvec_orb(:)
    real(kind=dp), allocatable                  :: vel(:,:), pos(:,:), pos_(:,:)
    real(kind=dp), allocatable                  :: pbest(:,:), gbest(:)
    real(kind=dp), allocatable                  :: pbest_history(:,:,:)
    real(kind=dp), allocatable                  :: cpbest(:), cpbest_plain(:), cpbest_orb(:)
    real(kind=dp)                               :: cgbest, cgbest_plain, cost
    real(kind=dp)                               ::         cgbest_orb
    real(kind=dp), allocatable                  :: costs(:), costs_plain(:), costs_orb(:)
    real(kind=dp), allocatable                  :: costs_(:), costs_plain_(:), costs_orb_(:)
    integer(kind=sp)                               i, min_id
    integer(kind=sp)                               imode, ldjac
    integer(kind=sp)                               min_loc(1), pso_miter
    logical                                        flag_order, flag_go
    real(kind=dp)                                  c1, c2, w, r1, r2, r3
    integer(kind=sp)                               iseed, iseed_
    integer(kind=sp)                               mpierr
    integer(kind=sp)                               iselect_mode
    real(kind=dp), external                     :: enorm
    real(kind=dp)                                  fnorm, fnorm_plain, fnorm_orb
    real(kind=dp)                                  pmax, pmin
    external                                       get_eig
    integer(kind=sp)                               info, maxfev, bestn
    real(kind=dp)                                  epsfcn, factor, xtol, ftol, gtol
    logical                                        flag_with_lmdif
    logical                                        flag_fit_plain, flag_fit_orbital
    integer(kind=sp), allocatable               :: ourgroup(:), ourjob(:)
    integer(kind=sp)                               mygroup(0:nprocs-1), groupid
#ifdef MPI
    integer(kind=sp)                               mpistat(MPI_STATUS_SIZE)
#endif    

   !PPRAM%flag_fit_plain = .TRUE.
    flag_fit_orbital= PINPT%flag_fit_orbital
    flag_fit_plain  = PPRAM%flag_fit_plain
    flag_with_lmdif = PINPT%flag_pso_with_lmdif
    flag_order      = .FALSE.
    imode           = 13
    n_particles     = PPRAM%pso_nparticles
    nparam_free     = PPRAM%nparam_free
    c1              = PPRAM%pso_c1
    c2              = PPRAM%pso_c2
    w               = PPRAM%pso_w 
    iseed           = iseed_ 
    r1              = 1d0
    r2              = 1d0
    r3              = 0d0
    bestn           = int(real(PPRAM%pso_nparticles) * PPRAM%pso_report_ratio)
    call random_init(iseed)
    if(flag_fit_orbital .and. flag_fit_plain) then
        flag_fit_plain = .FALSE. ! we consider fit_orbital instead of fit_plain
    endif

    if(flag_fit_orbital) then
      iselect_mode       = 3
    elseif(flag_fit_plain) then        
      iselect_mode       = 2
    else
      iselect_mode       = 1
    endif

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    !Note: n_particles will be solved in n_groups of cpu family.
    !      Each group will have n_member ~ nprocs/ngroup
    !      The n_member will be further parallized to solve eigenvalue problem in get_eig routine
    call get_npar_kpar(trim(PINPT%ifilenm(1)))
    allocate(ourgroup(npar)); allocate(ourjob(npar))
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
    maxfev = PINPT%miter * ( PPRAM%nparam_free + 1 )
    ftol = PINPT%ftol    ;xtol = PINPT%ptol ; gtol = 0.0D+00; epsfcn = 0.000D+00

    if(bestn .lt. npar) bestn = npar
    if(bestn .gt. PPRAM%pso_nparticles) bestn = PPRAM%pso_nparticles
    if(allocated(PPRAM%pso_cost_history)) deallocate(PPRAM%pso_cost_history)
    if(allocated(PPRAM%pso_cost_history_i)) deallocate(PPRAM%pso_cost_history_i)
    allocate(PPRAM%pso_cost_history(pso_miter))
    allocate(PPRAM%pso_cost_history_i(pso_miter,n_particles))
    PPRAM%pso_cost_history = 0d0
    PPRAM%pso_cost_history_i = 0d0

    do i = 1, PINPT%nsystem
      call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i), .FALSE., 0)
      if(allocated(ETBA(i)%dE)) deallocate(ETBA(i)%dE)
      allocate(ETBA(i)%dE( size(ETBA(i)%E(:,1)) , size(ETBA(i)%E(1,:)) ))
      ETBA(i)%E = 0d0 ; ETBA(i)%V = 0d0 ; ETBA(i)%SV = 0d0 ; ETBA(i)%dE=0d0
    enddo

    allocate(vel(  n_particles, nparam_free))
    allocate(pos(  n_particles, nparam_free))
    allocate(pbest(n_particles, nparam_free))
    if(allocated(PPRAM%pso_pbest_history)) deallocate(PPRAM%pso_pbest_history)
    allocate(PPRAM%pso_pbest_history(pso_miter, n_particles, nparam_free))
    allocate(gbest(              nparam_free))
    allocate(cpbest(n_particles             ))
    allocate(cpbest_plain(n_particles             ))
    allocate(cpbest_orb(n_particles             ))
    allocate(costs(n_particles              ))
    allocate(costs_plain(n_particles        ))
    allocate(costs_orb(n_particles        ))
    allocate(costs_(n_particles              ))
    allocate(costs_plain_(n_particles        ))
    allocate(costs_orb_(n_particles        ))
    costs = 0.d0 ; costs_plain = 0.d0 ; costs_ = 0.d0 ; costs_plain_ = 0.d0
    costs_orb = 0d0 ; costs_orb_ = 0d0
    vel   = 0.d0 ; pos   = 0.d0 
    pbest = 0.d0 ; gbest = 0.d0
    cpbest = 0.d0 ; cpbest_plain = 0.d0 ; cpbest_orb = 0d0
    cgbest = 0.d0 ; cgbest_plain = 0.d0 ; cgbest_orb = 0d0
    PPRAM%pso_pbest_history = 0d0
                            
    ldjac      = 0
    do i = 1, PINPT%nsystem
        ldjac = ldjac + PKPTS(i)%nkpoint * PGEOM(i)%nband * PINPT%nspin
    enddo
    if(ldjac .ne. 1) then
      allocate(fvec(ldjac))
      allocate(fvec_plain(ldjac))
    ! if(flag_fit_orbital) allocate(fvec_orb(ldjac))
      allocate(fvec_orb(ldjac))
    endif
    fvec = 0d0 ; fvec_plain = 0d0 
    if(flag_fit_orbital) fvec_orb = 0d0

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
           if( PPRAM%param( PPRAM%iparam_free(iparam) ) .ge. PPRAM%param_const(2, PPRAM%iparam_free(iparam) ) .or. &
               PPRAM%param( PPRAM%iparam_free(iparam) ) .le. PPRAM%param_const(3, PPRAM%iparam_free(iparam) )) then
             pmax = PPRAM%param_const(2, PPRAM%iparam_free(iparam) ) 
             pmin = PPRAM%param_const(3, PPRAM%iparam_free(iparam) )
             pos(iptcl, iparam) = (pmax+pmin)*0.5d0*r2 + r1
           else
             pos(iptcl, iparam) = PPRAM%param( PPRAM%iparam_free(iparam) )*r2 + r1
           endif
           vel(iptcl, iparam) = r1

           if( pos(iptcl, iparam) .ge. PPRAM%param_const(2, PPRAM%iparam_free(iparam) ) .or. &
               pos(iptcl, iparam) .le. PPRAM%param_const(3, PPRAM%iparam_free(iparam) )) then
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
        call get_dE(fvec, fvec_plain, fvec_orb, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, flag_fit_orbital)
        costs(iptcl) = enorm(ldjac, fvec)
        costs_plain(iptcl) = enorm(ldjac, fvec_plain)
        if(flag_fit_orbital) costs_orb(iptcl)   = enorm(ldjac, fvec_orb)
      elseif(flag_with_lmdif) then
        call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
                 ftol, xtol, gtol, fnorm, fnorm_plain, fnorm_orb, maxfev, epsfcn, factor, info, .FALSE., flag_fit_orbital)
        costs(iptcl) = fnorm
        costs_plain(iptcl) = fnorm_plain
        pos(iptcl, :) = PPRAM%param(PPRAM%iparam_free(:))
        if(flag_fit_orbital) costs_orb(iptcl)  = fnorm_orb
      endif
      if(COMM_KOREA%flag_split) then
        if(iselect_mode .eq. 3) then
          if(COMM_KOREA%myid .eq. 0) write(6,'(A, i5, 3(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl), &
                                                                                                                     " ,COST(orbital   ) = ", costs_orb(iptcl)
        else
          if(COMM_KOREA%myid .eq. 0) write(6,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl)
        endif
      else
        if(iselect_mode .eq. 3) then
          if_main                    write(6,'(A, i5, 3(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl), &
                                                                                                                     " ,COST(orbital   ) = ", costs_orb(iptcl)
        else
          if_main                    write(6,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl)
        endif
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
        if(flag_fit_orbital) then
          call MPI_SEND(costs_orb(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1))), ourjob(groupid+1), &
                        MPI_REAL8, 0, 222, mpi_comm_earth, mpierr)
        endif
       !if(flag_with_lmdif) then
        call MPI_SEND(pos(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1)), :), ourjob(groupid+1)*PPRAM%nparam_free, &
                      MPI_REAL8, 0, 110, mpi_comm_earth, mpierr)
       !endif
        call MPI_SEND(vel(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1)), :), ourjob(groupid+1)*PPRAM%nparam_free, &
                      MPI_REAL8, 0, 110, mpi_comm_earth, mpierr)
       
      elseif( myid .eq. 0 ) then
        do i = 1, npar - 1
            call MPI_RECV(costs(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                          MPI_REAL8, COMM_KOREA%group_main(i+1), 111, mpi_comm_earth, mpistat, mpierr)
            call MPI_RECV(costs_plain(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                          MPI_REAL8, COMM_KOREA%group_main(i+1), 222, mpi_comm_earth, mpistat, mpierr)
            if(flag_fit_orbital) then
              call MPI_RECV(costs_orb(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                            MPI_REAL8, COMM_KOREA%group_main(i+1), 222, mpi_comm_earth, mpistat, mpierr)
            endif
           !if(flag_with_lmdif) then
            call MPI_RECV(pos(sum(ourjob(1:i))+1:sum(ourjob(1:i+1)),:), ourjob(i+1)*PPRAM%nparam_free, &
            MPI_REAL8, COMM_KOREA%group_main(i+1), 110, mpi_comm_earth, mpistat, mpierr)
           !endif
            call MPI_RECV(vel(sum(ourjob(1:i))+1:sum(ourjob(1:i+1)),:), ourjob(i+1)*PPRAM%nparam_free, &
            MPI_REAL8, COMM_KOREA%group_main(i+1), 110, mpi_comm_earth, mpistat, mpierr)
        enddo
      endif
    endif
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif
    ! report initial cost
    do iptcl = 1, n_particles
      if(iselect_mode .eq. 3) then
        write(message,'(A, i5, 3(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl), &
                                                                                              " ,COST(orbital   ) = ", costs_orb(iptcl); write_msg_file
      else
        write(message,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl); write_msg_file
      endif
    enddo

    ! sort costs and pick bestn particles only and save
    call get_bestn_particles(bestn, n_particles, nparam_free, pos, vel, costs, costs_plain, costs_orb)

    ! report sorted cost
    write(message,'(A)')" "; write_msg
    do iptcl = 1, bestn
      if(iselect_mode .eq. 3) then
        write(message,'(A, i5, 3(A, F20.8))')"  Best PTCL: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl), &
                                                                                              " ,COST(orbital   ) = ", costs_orb(iptcl); write_msg
      else
        write(message,'(A, i5, 2(A, F20.8))')"  Best PTCL: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl); write_msg
      endif
    enddo

    ! find best particles and costs
    if(myid .eq. 0) then
      if(iselect_mode .eq. 1) then
        min_loc = minloc(costs(:))
      elseif(iselect_mode .eq. 2 ) then
        min_loc = minloc(costs_plain(:))
      elseif(iselect_mode .eq. 3) then
        min_loc = minloc(costs_orb(:))
      endif
      min_id= min_loc(1)
      gbest = pos(min_loc(1),:)
      cgbest= costs(min_loc(1))
      cgbest_plain = costs_plain(min_loc(1))
      cgbest_orb = costs_orb(min_loc(1))
      cpbest = costs
      cpbest_plain = costs_plain
      cpbest_orb   = costs_orb
      pbest = pos
    endif

#ifdef MPI
    ! share cost information with other members
    call MPI_BCAST(gbest       , size(gbest) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    call MPI_BCAST(cgbest      , 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    call MPI_BCAST(cgbest_plain, 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    if(iselect_mode .eq. 3) then
      call MPI_BCAST(cgbest_orb, 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    endif
    call MPI_BCAST(pbest       , size(pbest) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    call MPI_BCAST(cpbest      , size(cpbest),  MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

    ! report best one
    if(iselect_mode .eq. 1 .or. iselect_mode .eq. 2) then
      write(message,'(A, i5, 2(A, F24.12)      )')"  PSO INITIAL RESULT: iter =",0, &
                                                  ", GBest COST = ", cgbest, " ,COST(w/o weight) = ", cgbest_plain ; write_msg
!                                                 ", Best Particle ID = ", min_id ; write_msg
    elseif(iselect_mode .eq. 3) then
      write(message,'(A, i5, 3(A, F24.12)      )')"  PSO INITIAL RESULT: iter =",0, &
                                                  ", GBest COST = ", cgbest, " ,COST(w/o weight) = ", cgbest_plain, &
                                                  ", COST(orbital) = ", cgbest_orb ; write_msg
!                                                 ", Best Particle ID = ", min_id ; write_msg
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! main PSO loop !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(message,'(A)')" "; write_msg
    if(iselect_mode .eq. 2) then
        write(message,'(A)')"  Main PSO Loop start: minimizing COST w/o WEIGHT" ; write_msg
    elseif(iselect_mode .eq. 1) then
        write(message,'(A)')"  Main PSO Loop start: minimizing COST " ; write_msg
    elseif(iselect_mode .eq. 3) then
        write(message,'(A)')"  Main PSO Loop start: minimizing COST w orbital" ; write_msg
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
          if_main r3 = (random()*2d0 - 1d0) * PPRAM%pso_max_noise_amplitude
          call MPI_BCAST(r1, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
          call MPI_BCAST(r2, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
          call MPI_BCAST(r3, 1, MPI_REAL8, 0, mpi_comm_earth, mpierr)
#else
          r1 = random() ; r2 = random() ; r3 = (random()*2d0 - 1d0) * PPRAM%pso_max_noise_amplitude
#endif
          vel(iptcl, iparam) = w * vel(iptcl, iparam) + c1 * r1 * (pbest(iptcl, iparam) - pos(iptcl, iparam)) &
                                                      + c2 * r2 * (gbest(iparam)        - pos(iptcl, iparam)) + r3*0.3d0
        enddo
      enddo
      pos = pos + vel

      ! get cost for each particle 
      costs = 0d0 ; costs_plain = 0d0
 ptcl:do iptcl = sum(ourjob(1:groupid)) + 1, sum(ourjob(1:groupid+1))
        PPRAM%param(PPRAM%iparam_free(:)) = pos(iptcl, :)

        if(.not. flag_with_lmdif) then
          call get_dE(fvec, fvec_plain, fvec_orb, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, flag_fit_orbital)
          costs(iptcl) = enorm(ldjac, fvec)
          costs_plain(iptcl) = enorm(ldjac, fvec_plain)
          if(flag_fit_orbital) costs_orb(iptcl)   = enorm(ldjac, fvec_orb)
        elseif(flag_with_lmdif) then
          call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
                   ftol, xtol, gtol, fnorm, fnorm_plain, fnorm_orb, maxfev, epsfcn, factor, info, .FALSE., flag_fit_orbital)
          costs(iptcl) = fnorm
          costs_plain(iptcl) = fnorm_plain
          if(flag_fit_orbital) costs_orb(iptcl)   = fnorm_orb
          pos(iptcl, :) = PPRAM%param(PPRAM%iparam_free(:))
        endif
        if(COMM_KOREA%flag_split) then
          if(iselect_mode .eq. 3) then
            if(COMM_KOREA%myid .eq. 0) write(6,'(A, i5, 3(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl), &
                                                                                                                       " ,COST(orbital   ) = ", costs_orb(iptcl);
          else
            if(COMM_KOREA%myid .eq. 0) write(6,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl)
          endif
        else
          if(iselect_mode .eq. 3) then
            if_main                    write(6,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl),  &
                                                                                                                       " ,COST(orbital   ) = ", costs_plain(iptcl)
          else
            if_main                    write(6,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl)
          endif
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
          if(iselect_mode .eq. 3) then
            call MPI_SEND(costs_orb(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1))), ourjob(groupid+1), &
                          MPI_REAL8, 0, 444, mpi_comm_earth, mpierr)
          endif
          call MPI_SEND(pos(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1)),:), ourjob(groupid+1)*PPRAM%nparam_free, &
                        MPI_REAL8, 0, 330, mpi_comm_earth, mpierr)
          call MPI_SEND(vel(sum(ourjob(1:groupid))+1:sum(ourjob(1:groupid+1)),:), ourjob(groupid+1)*PPRAM%nparam_free, &
                        MPI_REAL8, 0, 330, mpi_comm_earth, mpierr)
        elseif( myid .eq. 0 ) then
          do i = 1, npar - 1
              call MPI_RECV(costs(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                            MPI_REAL8, COMM_KOREA%group_main(i+1), 333, mpi_comm_earth, mpistat, mpierr)
              call MPI_RECV(costs_plain(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                            MPI_REAL8, COMM_KOREA%group_main(i+1), 444, mpi_comm_earth, mpistat, mpierr)
              if(iselect_mode .eq. 3) then
                call MPI_RECV(costs_orb(sum(ourjob(1:i))+1:sum(ourjob(1:i+1))), ourjob(i+1), &
                              MPI_REAL8, COMM_KOREA%group_main(i+1), 444, mpi_comm_earth, mpistat, mpierr)
              endif
              call MPI_RECV(pos(sum(ourjob(1:i))+1:sum(ourjob(1:i+1)),:), ourjob(i+1)*PPRAM%nparam_free, &
                            MPI_REAL8, COMM_KOREA%group_main(i+1), 330, mpi_comm_earth, mpistat, mpierr)
              call MPI_RECV(vel(sum(ourjob(1:i))+1:sum(ourjob(1:i+1)),:), ourjob(i+1)*PPRAM%nparam_free, &
                            MPI_REAL8, COMM_KOREA%group_main(i+1), 330, mpi_comm_earth, mpistat, mpierr)
          enddo
        endif
      endif
      call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif


      ! report cost
      do iptcl = 1, n_particles
        if(iselect_mode .eq. 3) then
          write(message,'(A, i5, 3(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl), &
                                                                                                " ,COST(orbital   ) = ", costs_orb(iptcl); write_msg_file
        else
          write(message,'(A, i5, 2(A, F20.8))')"   Particle: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl); write_msg_file
        endif
      enddo

      ! sort costs and pick bestn particles only and save
      call get_bestn_particles(bestn, n_particles, nparam_free, pos, vel, costs, costs_plain, costs_orb)

      ! report sorted cost
      write(message,'(A)')" "; write_msg
      do iptcl = 1, bestn
        if(iselect_mode .eq. 3) then
          write(message,'(A, i5, 3(A, F20.8))')"  Best PTCL: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl), &
                                                                                                " ,COST(orbital   ) = ", costs_orb(iptcl); write_msg
        else
          write(message,'(A, i5, 2(A, F20.8))')"  Best PTCL: ",iptcl, " COST = ", costs(iptcl), " ,COST(w/o weight) = ", costs_plain(iptcl); write_msg
        endif
      enddo

      ! update gbest, pbest
      if(myid .eq. 0) then
        if(iselect_mode .eq. 1) then
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
              cpbest_plain(iptcl)   = costs_plain(iptcl)
               pbest(iptcl,:) = pos(iptcl,:)
            endif
          enddo
        elseif(iselect_mode .eq. 2) then
          min_loc = minloc(costs_plain(:))
          if( costs_plain(min_loc(1)) .lt. cgbest_plain) then
            cgbest = costs(min_loc(1))
            cgbest_plain = costs_plain(min_loc(1))
            gbest  = pos(min_loc(1),:)
            min_id = min_loc(1)
          endif
          do iptcl = 1, n_particles
            if( costs_plain(iptcl) .lt. cpbest_plain(iptcl) ) then
              cpbest(iptcl)         = costs(iptcl)
              cpbest_plain(iptcl)   = costs_plain(iptcl)
              pbest(iptcl,:) = pos(iptcl,:)
            endif
          enddo
        elseif(iselect_mode .eq. 3) then
          min_loc = minloc(costs_orb(:))
          if( costs_orb(min_loc(1)) .lt. cgbest_orb) then
            cgbest = costs(min_loc(1))
            cgbest_plain = costs_plain(min_loc(1))
            cgbest_orb = costs_orb(min_loc(1))
            gbest  = pos(min_loc(1),:)
            min_id = min_loc(1)
          endif
          do iptcl = 1, n_particles
            if( costs_orb(iptcl) .lt. cpbest_orb(iptcl) ) then
              cpbest(iptcl)         = costs(iptcl)
              cpbest_plain(iptcl)   = costs_plain(iptcl)
              cpbest_orb(iptcl)   = costs_orb(iptcl)
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
      if(iselect_mode .eq. 3) then
        call MPI_BCAST(cgbest_orb, 1           ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      endif
      call MPI_BCAST(pbest       , size(pbest) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(cpbest      , size(cpbest),  MPI_REAL8, 0, mpi_comm_earth, mpierr)
#endif

      if(iselect_mode .eq. 1) then
        PPRAM%pso_cost_history(iter) = cgbest
        PPRAM%pso_cost_history_i(iter,:) = cpbest
      elseif(iselect_mode .eq. 2) then
        PPRAM%pso_cost_history(iter) = cgbest_plain
        PPRAM%pso_cost_history_i(iter,:) = cpbest_plain
      elseif(iselect_mode .eq. 3) then
        PPRAM%pso_cost_history(iter) = cgbest_orb
        PPRAM%pso_cost_history_i(iter,:) = cpbest_orb
      endif
      PPRAM%pso_pbest_history(iter,:,:) = pbest
      PPRAM%niter                  = iter

      if(iselect_mode .eq. 1 .or. iselect_mode .eq. 2) then
        write(message,'(A, i5, 2(A, F24.12)      )')"  PSO RESULT: iter =",iter, &
                                                    ", GBest COST = ", cgbest, " ,COST(w/o weight) = ", cgbest_plain ; write_msg
!                                                   ", Best Particle ID = ", min_id ; write_msg
      elseif(iselect_mode .eq. 3) then
        write(message,'(A, i5, 3(A, F24.12)      )')"  PSO RESULT: iter =",iter, &
                                                    ", GBest COST = ", cgbest, " ,COST(w/o weight) = ", cgbest_plain, &
                                                    ", COST(orbital) = ", cgbest_orb  ; write_msg
!                                                   ", Best Particle ID = ", min_id ; write_msg
      endif
    enddo it ! iter

    ! calculate with final parameters and save costs ETBA%dE
    PPRAM%param(PPRAM%iparam_free(:)) = gbest(:)
    if(.not. flag_with_lmdif) then
      ldjac = 1
      call get_dE(fvec, fvec_plain, fvec_orb, ldjac, imode, PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, flag_fit_orbital)
    elseif(flag_with_lmdif) then
      call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm, fnorm_plain, fnorm_orb, maxfev, epsfcn, factor, info, .FALSE., flag_fit_orbital)
      call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, flag_fit_orbital)
    endif

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    COMM_KOREA%flag_split = .FALSE.
#endif
   
    return
endsubroutine

subroutine get_bestn_particles(bestn, n_particles, nparam, pos, vel, costs, costs_plain, costs_orb)
    use sorting
    use mpi_setup
    implicit none
    integer*4   iptcl
    integer*4   bestn
    integer*4   n_particles
    integer*4   nparam
    integer*4   icosts(n_particles)
    integer*4   icosts_bestn(bestn)
    real*8      pos(n_particles,nparam)
    real*8      pos_(n_particles, nparam)
    real*8      vel(n_particles,nparam)
    real*8      vel_(n_particles, nparam)
    real*8      costs(n_particles)
    real*8      costs_(n_particles)
    real*8      costs_plain(n_particles)
    real*8      costs_plain_(n_particles)
    real*8      costs_orb(n_particles)
    real*8      costs_orb_(n_particles)
    integer*4   mpierr

    call get_sort_index_1D(icosts, costs, n_particles, 'ascending ')
    icosts_bestn = icosts(1:bestn)

    do iptcl = 1, n_particles
      pos_(iptcl,:)       = pos(         icosts_bestn( mod(iptcl-1,bestn)+1) ,:)
      vel_(iptcl,:)       = vel(         icosts_bestn( mod(iptcl-1,bestn)+1) ,:)
      costs_(iptcl)       = costs(       icosts_bestn( mod(iptcl-1,bestn)+1) )
      costs_plain_(iptcl) = costs_plain( icosts_bestn( mod(iptcl-1,bestn)+1) )
      costs_orb_(iptcl)   = costs_orb(   icosts_bestn( mod(iptcl-1,bestn)+1) )
    enddo

    pos         = pos_
    vel         = vel_
    costs       = costs_
    costs_plain = costs_plain_
    costs_orb   = costs_orb_
 
#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
    if(COMM_KOREA%flag_split) then
      call MPI_BCAST(pos         , size(pos)         ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(vel         , size(vel)         ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(costs       , size(costs)       ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(costs_plain , size(costs_plain) ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
      call MPI_BCAST(costs_orb   , size(costs_orb)   ,  MPI_REAL8, 0, mpi_comm_earth, mpierr)
    endif
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

    return
endsubroutine
