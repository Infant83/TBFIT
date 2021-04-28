#include "alias.inc"
! Note: if flag_tbfit = .true.,  program stops at the end of the routine get_fit with the band structure with fitted parameters
subroutine get_fit(PINPT, PPRAM_FIT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PINPT_BERRY, PINPT_DOS)
   use parameters
   use mpi_setup
!  use kill
   use print_io
   use time
   implicit none
   integer*4         mpierr
   external          get_eig
   type(incar)                             :: PINPT 
   type(params )                           :: PPRAM_FIT
   type(params ), dimension(PINPT%nsystem) :: PPRAM
   type(kpoints), dimension(PINPT%nsystem) :: PKPTS
   type(energy ), dimension(PINPT%nsystem) :: EDFT
   type(weight ), dimension(PINPT%nsystem) :: PWGHT
   type(poscar ), dimension(PINPT%nsystem) :: PGEOM
   type(hopping), dimension(PINPT%nsystem) :: NN_TABLE
   type(gainp  ), dimension(PINPT%nsystem) :: PKAIA
   type(berry  ), dimension(PINPT%nsystem) :: PINPT_BERRY
   type(dos    ), dimension(PINPT%nsystem) :: PINPT_DOS
   logical                                    flag_exit
   integer*4                                  ifit
   integer*4                                  i
   integer*4                                  iseed   ! random seed, default = 123
   real*8                                     fnorm, fnorm_
   
   ! Important note: In this routine, we assume that PPRAM applies to all system,
   !                 so we just use only PPRAM(1).
   !                 This policy can be changed afterwards on the purpose if needed but at this moment
   !                 we keep in this way.

   if(PINPT%flag_tbfit_parse) then 
     if(.not. PINPT%flag_tbfit) return  ! if -nofit, -np is requested from command line, exit this routine
   else
     call check_tbfit_true(PINPT)
     if(.not. PINPT%flag_tbfit) return  ! if flag_tbfit = .false. from ifile, just exit this routine
   endif

   fnorm_ = 0d0        ; flag_exit = .false. ;    ifit = 0
   call time_check(t1, t0, 'init')

   call initialize_fit(PINPT, PPRAM_FIT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PKAIA, PINPT_BERRY, PINPT_DOS)

   call set_free_parameters(PPRAM_FIT)

   call report_init(PINPT, PPRAM_FIT, PKPTS, EDFT, PWGHT, PGEOM)


   if(trim(PINPT%ls_type) .eq. 'LMDIF' .or. trim(PINPT%ls_type) .eq. 'lmdif') then

 fit:do while( .not. flag_exit)
       write(message,'(A)')' '  ; write_msg
       write(message,'(A,I0,A)')' #-START ',ifit+1,'-th LMDIF run'  ; write_msg
       write(message,'(A)')' '  ; write_msg

       if(allocated(PPRAM_FIT%cost_history)) deallocate(PPRAM_FIT%cost_history)
       allocate(PPRAM_FIT%cost_history(PINPT%miter))
       PPRAM_FIT%cost_history = 0d0

       call leasqr_lm ( get_eig, NN_TABLE, EDFT, PWGHT, PINPT, PPRAM_FIT, PKPTS, PGEOM, fnorm)

       call check_conv_and_constraint(PPRAM_FIT, PINPT, flag_exit, ifit, fnorm, fnorm_)

     enddo fit
    
   elseif(trim(PINPT%ls_type) .eq. 'PSO'    .or. trim(PINPT%ls_type) .eq. 'pso') then
       write(message,'(A)')' '  ; write_msg
       write(message,'(A,I0,A)')' #-START ',ifit+1,'-th PSO run'  ; write_msg
       write(message,'(A)')' '  ; write_msg

       if(trim(PPRAM_FIT%pso_mode) .eq. 'pso' .or. trim(PPRAM_FIT%pso_mode) .eq. 'PSO') then
         call pso_fit ( PINPT, PPRAM_FIT, PKPTS, PWGHT, PGEOM, NN_TABLE, EDFT, PPRAM_FIT%pso_iseed, PPRAM_FIT%pso_miter  )
       elseif(trim(PPRAM_FIT%pso_mode) .eq. 'pso_bestn' .or. trim(PPRAM_FIT%pso_mode) .eq. 'PSO_BESTN') then
         call pso_fit_best ( PINPT, PPRAM_FIT, PKPTS, PWGHT, PGEOM, NN_TABLE, EDFT, PPRAM_FIT%pso_iseed, PPRAM_FIT%pso_miter  )
       endif

      !call check_conv_and_constraint(PPRAM_FIT, PINPT, flag_exit, ifit, fnorm, fnorm_)

   elseif(trim(PINPT%ls_type) .eq. 'PIKAIA' .or. trim(PINPT%ls_type) .eq. 'GA' .or. &
          trim(PINPT%ls_type) .eq. 'pikaia' .or. trim(PINPT%ls_type) .eq. 'ga' ) then

     ! Important NOTE: Current version does not support multiple systems to be fitted with genetic algorithm
     !                 => only main system is passed to the gen_algo subroutine.
     if(PPRAM_FIT%slater_koster_type .gt. 10) then
       write(message,'(A)')'    !WARN! Current version does not support to use SK_SCALE_TYPE >= 10 along with LSTYPE=GA or PIKAIYA'  ; write_msg
       write(message,'(A)')'           Please set LSTYPE to LMDIF. The option will be available in the near future. Good luck...'  ; write_msg
       kill_job
     else
       call gen_algo  ( get_eig, NN_TABLE(1), EDFT(1)%E, PWGHT(1), PINPT, PPRAM_FIT, PKPTS(1), PGEOM(1), PKAIA(1))
     endif

   endif

!  ! print fitted parameters to file
        ! NOTE: only PWGHT(1) info is printed along with
        !       even though PINPT%nsystem > 1, this is because we just consider PPRAM to be applied
        !       for entire systems in fitting process.
        !       To avoid any confusion, in the future release, this could be
        !       corrected or other approach in printing PPRAM can be considered, but in this version
        !       we just keep this stratege for the convenience... 31.Oct.2020 HJK
   if_main call print_param (PINPT, PPRAM_FIT, PWGHT(1), PPRAM_FIT%pfileoutnm, .TRUE.)
   if_main call print_param (PINPT, PPRAM_FIT, PWGHT(1),  '   Fitted param(i):', .FALSE.)
   if(trim(PINPT%ls_type) .eq. 'PSO'    .or. trim(PINPT%ls_type) .eq. 'pso') then
     if(PINPT%flag_pso_report_particles) then
       if_main call print_param_pso(PINPT, PPRAM_FIT, PWGHT(1))
     endif

   endif

   write(message,'(A,A)')'  Fitted parameters will be written in ',PPRAM_FIT%pfileoutnm  ; write_msg

   call time_check(t1, t0)
   write(message,*)' ' ; write_msg
   write(message,'(A,F12.6)')"  TIME ELAPSED for FITTING (s)", t1 ; write_msg
   write(message,'( A)')' '  ; write_msg
   write(message,'( A)')' #===================================================='  ; write_msg
   write(message,'(2A)')'   END FITTING PROCEDURE: ', trim(PINPT%ls_type)  ; write_msg
   write(message,'( A)')' #===================================================='  ; write_msg
   write(message,'( A)')' '  ; write_msg

   PINPT%flag_tbfit_finish = .true.

   return
endsubroutine

subroutine initialize_fit(PINPT, PPRAM, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PKAIA, PINPT_BERRY, PINPT_DOS)
   use parameters, only : incar, params, kpoints, energy, weight, poscar, hopping, gainp
   use parameters, only : dos, berry, replot
   use print_io
   use mpi_setup
   implicit none
   type(incar  )                           :: PINPT, PINPT_
   type(params )                           :: PPRAM, PPRAM1, PPRAM_
   type(kpoints), dimension(PINPT%nsystem) :: PKPTS
   type(energy ), dimension(PINPT%nsystem) :: EDFT 
   type(weight ), dimension(PINPT%nsystem) :: PWGHT
   type(poscar ), dimension(PINPT%nsystem) :: PGEOM 
   type(hopping), dimension(PINPT%nsystem) :: NN_TABLE
   type(gainp  ), dimension(PINPT%nsystem) :: PKAIA   
   type(berry  ), dimension(PINPT%nsystem) :: PINPT_BERRY 
   type(dos    ), dimension(PINPT%nsystem) :: PINPT_DOS  

   type(replot)                             :: PRPLT ! dummy
   integer*4                                   i

   if(PINPT%flag_python_module) then
     print_mode = 99
   else
     print_mode = 3
   endif

   do i = 1, PINPT%nsystem
     call read_input(PINPT,PPRAM_,PKPTS(i), PGEOM(i),PWGHT(i),EDFT(i),NN_TABLE(i), &
                     PINPT_DOS(i),PINPT_BERRY(i),PKAIA(i),PRPLT, i)

     ! save main input PINPT and PPRAM for other systems
     if(i .eq. 1) then
       PINPT_ = PINPT  ! main input tags are read from system 1
       PPRAM1 = PPRAM_  ! main parameters are read from system 1
     elseif(i .gt. 1) then
       PINPT_%title(i) = PINPT%title(i)
       PINPT_%ifilenm(i) = PINPT%ifilenm(i)
     endif

   enddo
   
   PINPT = PINPT_
   PPRAM = PPRAM1

   return
endsubroutine
subroutine report_init(PINPT, PPRAM, PKPTS, EDFT, PWGHT, PGEOM)
   use parameters, only : params, incar, kpoints, energy, weight, poscar
   use print_io
   use mpi_setup
   implicit none
   type (incar)                             :: PINPT       ! parameters for input arguments
   type (params)                            :: PPRAM       ! parameters for input arguments
   type (poscar) , dimension(PINPT%nsystem) :: PGEOM       ! parameters for geometry info
   type (kpoints), dimension(PINPT%nsystem) :: PKPTS       ! parameters for kpoints
   type (energy) , dimension(PINPT%nsystem) :: EDFT        ! target energy to be fitted to
   type (weight) , dimension(PINPT%nsystem) :: PWGHT       ! weight factor for the fitting to target energy
   integer*4                                   mpierr
   integer*4                                   i
   character*80                                fname, add_LMDIF
   logical                                     flag_with_lmdif

   flag_with_lmdif = .false.
   add_LMDIF = ''
   if(trim(PINPT%ls_type) .eq. 'PSO') then
     if(PINPT%flag_pso_with_lmdif)  flag_with_lmdif = .true.
     add_LMDIF='+LMDIF'
   elseif(trim(PINPT%ls_type) .eq. 'GA') then
     if(PINPT%flag_ga_with_lmdif)  flag_with_lmdif = .true.
     add_LMDIF='+LMDIF'
   endif

   write(message,'( A)')' '  ; write_msg
   write(message,'( A)')' #===================================================='  ; write_msg
   write(message,'(2A)')'   START FITTING PROCEDURE: ', trim(PINPT%ls_type)//trim(add_LMDIF)  ; write_msg
   write(message,'( A)')' #===================================================='  ; write_msg
   write(message,'( A)')' '  ; write_msg
   write(message,'( A,I0)')'  N_PARAM (free)     : ', PPRAM%nparam_free  ; write_msg
   if(trim(PINPT%ls_type) .eq. 'PSO') then
     write(message,'( A,I0)')'  PSO_MITER          : ', PPRAM%pso_miter           ; write_msg
     write(message,'( A,I0)')'  P_PARTICLE(PSO_NP) : ', PPRAM%pso_nparticles      ; write_msg
     write(message,'( A,3F10.4)')'  PSO_OPT(C1,C2,W)   : ', PPRAM%pso_c1,PPRAM%pso_c2,PPRAM%pso_w ; write_msg
     write(message,'( A, F10.4)')'  PSO_NOISE          : ', PPRAM%pso_max_noise_amplitude ; write_msg
     write(message,'( A,A )')'  PSO_MODE           : ', trim(PPRAM%pso_mode)       ; write_msg
     if(trim(PPRAM%pso_mode) .eq. 'pso_bestn' .or. trim(PPRAM%pso_mode) .eq. 'PSO_BESTN') then
       write(message,'( A,I0,A,F5.1,A)')'  PSO_BESTN          : ', int(real(PPRAM%pso_nparticles) * PPRAM%pso_report_ratio), ' (= ',PPRAM%pso_report_ratio*100d0,' %)'; write_msg
     endif

   endif
   write(message,'( A,I0)')'  N_SYSTEM           : ', PINPT%nsystem      ; write_msg
   do i = 1, PINPT%nsystem
     write(message,'( A,I0,2A)')'  GFILE (system-',i,')   : ', trim(PGEOM(i)%gfilenm); write_msg
   enddo
   write(message,'( A)')' '  ; write_msg

   !print target energy to file
   do i = 1, PINPT%nsystem
     if_main_then 
       fname = 'band_structure_DFT'//trim(PINPT%title(i))//'.dat'

       call print_energy_weight( PKPTS(i)%kpoint, PKPTS(i)%nkpoint, EDFT(i), PWGHT(i), PGEOM(i)%neig, &
                                 PINPT, trim(fname), PINPT%flag_get_band_order)
     if_main_end
   enddo

!  !replace WEIGHT info of INCAR with those in PFILE if USE_WEIGHT is activated.
   if(PINPT%flag_use_weight) then ! NOTE: only active if nsystem = 1
     if_main call rewrite_incar(PINPT, PPRAM, PWGHT(1))
   endif

   return
endsubroutine
subroutine check_conv_and_constraint(PPRAM, PINPT, flag_exit, ifit, fnorm, fnorm_) 
   use parameters, only : params, incar
   use print_io
   use mpi_setup  
   implicit none
   type (incar)   :: PINPT
   type (params)  :: PPRAM
   integer*4         i, j, k, ifit
   real*8            fnorm, fnorm_, fdiff
   logical           flag_conv, flag_exit

   flag_conv = .true.
   ifit = ifit + 1

   write(message,'(A)')' '  ; write_msg
   write(message,'(A,I0,A)')'           Check constraint: upper and lower bounds'  ; write_msg
   write(message,'(A,I0,A)')'    #========================================================='  ; write_msg

   do j = 1, PPRAM%nparam

     if(PPRAM%slater_koster_type .gt. 10) then
       do i = 1, PPRAM%param_nsub(j)
         k = nint(PPRAM%param_const_nrl(1,i,j))
         if( k .eq. 0) then
           if(PPRAM%param_nrl(i,j) .lt. PPRAM%param_const_nrl(3,i,j)) then ! lower bound
             if( PPRAM%param_nrl(i,j) .lt. 0d0 .and. i .eq. 4 .and. &
                -PPRAM%param_nrl(i,j) .ge. PPRAM%param_const_nrl(2,i,j)) then
               PPRAM%param_nrl(i,j) = PPRAM%param_const_nrl(2,i,j)
             elseif( PPRAM%param_nrl(i,j) .le. 0d0 .and. i .eq. 4 .and. &
                -PPRAM%param_nrl(i,j) .lt. PPRAM%param_const_nrl(2,i,j)) then
               PPRAM%param_nrl(i,j) = -PPRAM%param_nrl(i,j)
             else
               PPRAM%param_nrl(i,j) = PPRAM%param_const_nrl(3,i,j)
             endif
             flag_conv = .false.
             write(message,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PPRAM%param_name(j)),'(',i,') < ', &
                                                      PPRAM%param_const_nrl(3,i,j),' --> ', &
                                                      trim(PPRAM%param_name(j)),'(',i,') = ',PPRAM%param_nrl(i,j) ; write_msg
           elseif(PPRAM%param_nrl(i,j) .gt. PPRAM%param_const_nrl(2,i,j)) then !upper bound
             PPRAM%param_nrl(i,j) = PPRAM%param_const_nrl(2,i,j)
             flag_conv = .false.
             write(message,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PPRAM%param_name(j)),'(',i,') > ', &
                                                      PPRAM%param_const_nrl(2,i,j),' --> ', &
                                                      trim(PPRAM%param_name(j)),'(',i,') = ',PPRAM%param_nrl(i,j) ; write_msg
           endif
         elseif( k .ge. 1) then ! if same as k-th parameter
           if(PPRAM%param_nrl(i,j) .lt. PPRAM%param_const_nrl(3,i, k )) then !lower bound
             if(PPRAM%param_nrl(i,j) .lt. 0d0 .and. i .eq. 4 .and. &
                -PPRAM%param_nrl(i,j) .ge. PPRAM%param_const_nrl(2,i,k)) then
               PPRAM%param_nrl(i,j) = PPRAM%param_const_nrl(2,i,k)
             elseif( PPRAM%param_nrl(i,j) .le. 0d0 .and. i .eq. 4 .and. &
                -PPRAM%param_nrl(i,j) .lt. PPRAM%param_const_nrl(2,i,k)) then
               PPRAM%param_nrl(i,j) = -PPRAM%param_nrl(i,j)
             else
               PPRAM%param_nrl(i,j) = PPRAM%param_const_nrl(3,i, k )
             endif
             flag_conv = .false.
             write(message,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PPRAM%param_name(j)),'(',i,') < ', &
                                                      PPRAM%param_const_nrl(3,i,k),' --> ', &
                                                      trim(PPRAM%param_name(j)),'(',i,') = ',PPRAM%param_nrl(i,j) ; write_msg
           elseif(PPRAM%param_nrl(i,j) .gt. PPRAM%param_const_nrl(2,i, k)) then !upper bound
             PPRAM%param_nrl(i,j) = PPRAM%param_const_nrl(2,i, k)
             flag_conv = .false.
             write(message,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PPRAM%param_name(j)),'(',i,') > ', &
                                                      PPRAM%param_const_nrl(2,i,k),' --> ', &
                                                      trim(PPRAM%param_name(j)),'(',i,') = ',PPRAM%param_nrl(i,j) ; write_msg
           endif
         endif
       enddo
     else
       if( nint(PPRAM%param_const(1,j)) .eq. 0) then
         if(PPRAM%param(j) .lt. PPRAM%param_const(3,j)) then ! lower bound
           PPRAM%param(j) = PPRAM%param_const(3,j)
           flag_conv = .false.
           write(message,'(6x,2A,F10.4,3A,F10.4)')trim(PPRAM%param_name(j)),' < ',PPRAM%param_const(3,j),' --> ', &
                                                    trim(PPRAM%param_name(j)),' = ',PPRAM%param(j) ; write_msg
         elseif(PPRAM%param(j) .gt. PPRAM%param_const(2,j)) then !upper bound
           PPRAM%param(j) = PPRAM%param_const(2,j)
           flag_conv = .false.
           write(message,'(6x,2A,F10.4,3A,F10.4)')trim(PPRAM%param_name(j)),' > ',PPRAM%param_const(2,j),' --> ', &
                                                    trim(PPRAM%param_name(j)),' = ',PPRAM%param(j) ; write_msg
         endif
       elseif(nint(PPRAM%param_const(1,j)) .ge. 1) then
         if(PPRAM%param(j) .lt. PPRAM%param_const(3,nint(PPRAM%param_const(1,j)))) then !lower bound
           PPRAM%param(j) = PPRAM%param_const(3,nint(PPRAM%param_const(1,j)))
           flag_conv = .false.
           write(message,'(6x,2A,F10.4,3A,F10.4)')trim(PPRAM%param_name(j)),' < ',PPRAM%param_const(3,nint(PPRAM%param_const(1,j))),' --> ', &
                                                    trim(PPRAM%param_name(j)),' = ',PPRAM%param(j) ; write_msg
         elseif(PPRAM%param(j) .gt. PPRAM%param_const(2,nint(PPRAM%param_const(1,j)))) then !upper bound
           PPRAM%param(j) = PPRAM%param_const(2,nint(PPRAM%param_const(1,j)))
           flag_conv = .false.
           write(message,'(6x,2A,F10.4,3A,F10.4)')trim(PPRAM%param_name(j)),' > ',PPRAM%param_const(2,nint(PPRAM%param_const(1,j))),' --> ', &
                                                    trim(PPRAM%param_name(j)),' = ',PPRAM%param(j) ; write_msg
         endif
       endif
     endif
   enddo

   if(flag_conv) then
   write(message,'(A)')' '  ; write_msg
     write(message,'(A)')'           Check fitness function updates '  ; write_msg
     write(message,'(A)')'    #============================================'  ; write_msg
     ! 2. check fitness function updates
     fdiff = sqrt((fnorm_ - fnorm)**2)/fnorm*100
     if( fdiff .le. PINPT%fdiff) then
       flag_conv = .true.
       write(message,'(A,F16.8,A,F16.8)')'  FDIFF(',PINPT%fdiff,')>=', fdiff  ; write_msg
     elseif( (fdiff .gt. PINPT%fdiff ) ) then
       flag_conv = .false.
       write(message,'(A,F16.8,A,F16.8)')'  FDIFF(',PINPT%fdiff,')<=', fdiff  ; write_msg
     endif
   endif
   fnorm_ = fnorm

   if(PINPT%flag_plot_fit) then
     if_main call execute_command_line('gnuplot '//trim(PINPT%filenm_gnuplot), .false.)
   endif

   if(ifit .ge. PINPT%mxfit .or. flag_conv) then
     flag_exit = .true.
   else
     flag_exit = .false.
   endif

   return
endsubroutine

subroutine set_free_parameters(PPRAM)
   use parameters, only : params
   implicit none
   type (params)          :: PPRAM
   integer*4                 i, j, k
   integer*4, allocatable :: iparam_free(:)
   integer*4, allocatable :: iparam_free_nrl(:)
   
   if(PPRAM%slater_koster_type .gt. 10) then
     PPRAM%nparam_nrl = sum(PPRAM%param_nsub(1:PPRAM%nparam))
     allocate(iparam_free(PPRAM%nparam))
     allocate(iparam_free_nrl(PPRAM%nparam))
     PPRAM%nparam_nrl_free = 0
     PPRAM%nparam_free = 0
     do i = 1, PPRAM%nparam
       k = nint(PPRAM%param_const_nrl(1,1, i)) ! is equal to
       j = nint(PPRAM%param_const_nrl(4,1, i)) ! is fixed
       if( k .eq. 0 .and. j .eq. 0) then
         PPRAM%nparam_free     = PPRAM%nparam_free     + 1
         iparam_free_nrl(PPRAM%nparam_free) = PPRAM%nparam_nrl_free + 1
         iparam_free(PPRAM%nparam_free) = i
         PPRAM%nparam_nrl_free = PPRAM%nparam_nrl_free + PPRAM%param_nsub(i)
       endif
     enddo
     allocate(PPRAM%iparam_free(PPRAM%nparam_free))
     allocate(PPRAM%iparam_free_nrl(PPRAM%nparam_free))
     PPRAM%iparam_free(1:PPRAM%nparam_free) = iparam_free(1:PPRAM%nparam_free)
     PPRAM%iparam_free_nrl(1:PPRAM%nparam_free) = iparam_free_nrl(1:PPRAM%nparam_free)
     deallocate(iparam_free)
     deallocate(iparam_free_nrl)

   else
     allocate(iparam_free(PPRAM%nparam))
     PPRAM%nparam_free = 0
     do i = 1, PPRAM%nparam
       k = nint(PPRAM%param_const(1, i)) ! is equal to
       j = nint(PPRAM%param_const(4, i)) ! is fixed
       if( k .eq. 0 .and. j .eq. 0) then
         PPRAM%nparam_free = PPRAM%nparam_free + 1
         iparam_free(PPRAM%nparam_free) = i
       endif
     enddo
     allocate(PPRAM%iparam_free(PPRAM%nparam_free))
     PPRAM%iparam_free(1:PPRAM%nparam_free) = iparam_free(1:PPRAM%nparam_free)
     deallocate(iparam_free)
   endif

   return
endsubroutine

subroutine check_tbfit_true(PINPT)
   use parameters, only : incar, pid_incar
   use mpi_setup
   implicit none
   type(incar)   :: PINPT
   logical          flag_tbfit
   integer*4        i_continue, linecount
   integer*4        mpierr
   character*132    inputline
   character*40     desc_str
   character(*),parameter :: func = 'check_tbfit_true'

   ! assume that iflenm(1) is the main input file
   open(pid_incar, file=trim(PINPT%ifilenm(1)), iostat = i_continue) ; linecount = 0
   
   line: do
           read(pid_incar,'(A)', iostat=i_continue) inputline
           if(i_continue < 0) exit  ! end of file reached
           if(i_continue > 0) then
             call write_log('Unknown error reading file:'//trim(PINPT%ifilenm(1))//' '//trim(func),3,myid) ! check only master input ifilenm(1)
             kill_job
           endif

           linecount = linecount + 1

           ! check INPUT tag with TBFIT
           read(inputline, *, iostat=i_continue) desc_str
           if(i_continue .ne. 0) cycle              ! skip empty line
           if (desc_str(1:1).eq.'#') cycle  ! skip comment

           ! head
           select case (desc_str)
            
             case('TBFIT')
               read(inputline,*,iostat=i_continue) desc_str, PINPT%flag_tbfit
             case('SET')
               read(inputline,*,iostat=i_continue) desc_str, desc_str

               if(trim(desc_str) .eq. 'REPLOT') then
                  PINPT%flag_tbfit = .false.
                  close(pid_incar)
                  return
               endif
           end select
            
         enddo line

   close(pid_incar)

   return
endsubroutine

    subroutine check_kill_tbfit(PINPT,PPRAM, PWGHT)
      use parameters, only: incar, weight, params
      use print_io
      use mpi_setup
      logical        flag_exist
      integer*4      mpierr
      type(incar)                            :: PINPT
      type(params)                           :: PPRAM
      type(weight), dimension(PINPT%nsystem) :: PWGHT

      inquire(file="KILLFIT",exist=flag_exist)

#ifdef MPI
      if(COMM_KOREA%flag_split) then
        call MPI_BCAST(flag_exist, 1, MPI_LOGICAL, 0, COMM_KOREA%mpi_comm, mpierr)
      elseif(.not. COMM_KOREA%flag_split) then
        call MPI_BCAST(flag_exist, 1, MPI_LOGICAL, 0, mpi_comm_earth, mpierr)
      endif
#endif

      if(flag_exist) then

        if( PPRAM%slater_koster_type .gt. 10) then
          call update_param_nrl( PPRAM )
        else
          call update_param( PPRAM )
        endif

        if_main call execute_command_line('\rm -f KILLFIT')

        ! NOTE: only PWGHT(1) info is printed along with
        !       even though PINPT%nsystem > 1, this is because we just consider PPRAM to be applied
        !       for entire systems in fitting process. 
        !       To avoid any confusion, in the future release, this could be 
        !       corrected or other approach in printing PPRAM can be considered, but in this version
        !       we just keep this stratege for the convenience... 31.Oct.2020 HJK
        if_main call print_param (PINPT, PPRAM, PWGHT(1), PPRAM%pfileoutnm, .TRUE.)
        write(message,'(A)') ' '  ; write_msg
        write(message,'(A)') ' Termination of job is requested by providing KILLFIT file.' ; write_msg
        write(message,'(2A)') ' The latest updates of PARAMETERS will be written in ', trim(PPRAM%pfileoutnm) ; write_msg
        write(message,'(A)') ' Kill the job now...' ; write_msg
        write(message,'(A)') ' '  ; write_msg

        kill_job

      endif
      return
    endsubroutine

