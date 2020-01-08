#include "alias.inc"
subroutine get_fit(PINPT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PKAIA)
   use parameters
   use mpi_setup
   use kill
   implicit none
   integer*4         mpierr
   external          get_eig
   type (incar)   :: PINPT       ! parameters for input arguments
   type (dos)     :: PINPT_DOS   ! parameters for density of states calculation
   type (poscar)  :: PGEOM       ! parameters for geometry info
   type (kpoints) :: PKPTS       ! parameters for kpoints
   type (energy)  :: EDFT        ! target energy to be fitted to
   type (energy)  :: ETBA_FIT    ! final fitted energies
   type (weight)  :: PWGHT       ! weight factor for the fitting to target energy
   type (hopping) :: NN_TABLE    ! table for hopping index
   type (gainp)   :: PKAIA       ! input/control parameters for genetic algorithm
   logical           flag_conv
   integer*4         i,j,k, ifit,mxfit
   real*8            fdiff, fnorm(2), fnorm_
   character*132     gnu_command
   logical           flag_wait_plot 
   integer*4, allocatable :: iparam_free(:)
   integer*4, allocatable :: iparam_free_nrl(:)
   
   if(PINPT%flag_plot_fit) then
     allocate(ETBA_FIT%E(PINPT%nband*PINPT%nspin, PKPTS%nkpoint))
     allocate(ETBA_FIT%V(PGEOM%neig*PINPT%ispin,PINPT%nband*PINPT%nspin, PKPTS%nkpoint))
     flag_wait_plot = .false.
     write(gnu_command, '(A,A)')'gnuplot ', trim(PINPT%filenm_gnuplot)
   endif

   if(PINPT%flag_sparse) then
     if_main write(6,'(A)') '    !WARN! EWINDOW tag cannot be used in FITTING procedures'
     if_main write(6,'(A)') '           Please deactivate EWINDOW for the further process.'
     if_main write(6,'(A)') '           Exit program...'
     kill_job
   endif
   fnorm_ = 0d0
   flag_conv = .false.
   mxfit = PINPT%mxfit
   if(mxfit .le. 0) then
     if_main write(6,'(A)')'    !WARN! The current MXFIT value, maximum number of LMDIF attempt, is <= 0 .'
     if_main write(6,'(A)')'           Please make sure that MXFIT > 0 in your input tag.'
     kill_job
   endif
   ifit = 0

   if(PINPT%slater_koster_type .gt. 10) then
     PINPT%nparam_nrl = sum(PINPT%param_nsub(1:PINPT%nparam))
     allocate(iparam_free(PINPT%nparam))
     allocate(iparam_free_nrl(PINPT%nparam))
     PINPT%nparam_nrl_free = 0
     PINPT%nparam_free = 0
     do i = 1, PINPT%nparam
       k = nint(PINPT%param_const_nrl(1,1, i)) ! is equal to
       j = nint(PINPT%param_const_nrl(4,1, i)) ! is fixed
       if( k .eq. 0 .and. j .eq. 0) then
         PINPT%nparam_free     = PINPT%nparam_free     + 1
         iparam_free_nrl(PINPT%nparam_free) = PINPT%nparam_nrl_free + 1
         iparam_free(PINPT%nparam_free) = i
         PINPT%nparam_nrl_free = PINPT%nparam_nrl_free + PINPT%param_nsub(i)
       endif
     enddo
     allocate(PINPT%iparam_free(PINPT%nparam_free))
     allocate(PINPT%iparam_free_nrl(PINPT%nparam_free))
     PINPT%iparam_free(1:PINPT%nparam_free) = iparam_free(1:PINPT%nparam_free)
     PINPT%iparam_free_nrl(1:PINPT%nparam_free) = iparam_free_nrl(1:PINPT%nparam_free)
     deallocate(iparam_free)
     deallocate(iparam_free_nrl)

   else
     allocate(iparam_free(PINPT%nparam))
     PINPT%nparam_free = 0
     do i = 1, PINPT%nparam
       k = nint(PINPT%param_const(1, i)) ! is equal to
       j = nint(PINPT%param_const(4, i)) ! is fixed
       if( k .eq. 0 .and. j .eq. 0) then
         PINPT%nparam_free = PINPT%nparam_free + 1
         iparam_free(PINPT%nparam_free) = i
       endif
     enddo
     allocate(PINPT%iparam_free(PINPT%nparam_free))
     PINPT%iparam_free(1:PINPT%nparam_free) = iparam_free(1:PINPT%nparam_free)
     deallocate(iparam_free)

   endif
   if_main write(6,'( A)')' '
   if_main write(6,'( A)')' ====================================================='
   if_main write(6,'(2A)')'   START FITTING PROCEDURE: ', trim(PINPT%ls_type)
   if_main write(6,'( A)')' ====================================================='
   if_main write(6,'( A)')' '

   if_main write(6,'( A,I0)')'  N_PARAM (free):  ', PINPT%nparam_free
   if_main write(6,'( A)')' '

   !print target energy to file
   if_main call print_energy_weight( PKPTS%kpoint, PKPTS%nkpoint, EDFT, PWGHT, PGEOM%neig, &
                                     PINPT, 'band_structure_DFT.dat')

   !initial parameter
   if_main call print_param ( PINPT, 0, '  Initial param(i):', .FALSE.)

   ! fitting parameter
   if(trim(PINPT%ls_type) .eq. 'LMDIF' .or. trim(PINPT%ls_type) .eq. 'lmdif') then
 fit:do while( ifit .le. mxfit .or. .not. flag_conv)
       ifit = ifit + 1

       if_main write(6,'(A)')' '
       if_main write(6,'(A,I0,A)')' **START ',ifit,'-th LMDIF run'
       if_main write(6,'(A)')' '

       call leasqr_lm ( get_eig, NN_TABLE, EDFT%E, PGEOM%neig, PWGHT, PINPT, PKPTS, fnorm)

       if_main write(6,'(A,I0,A)')'           Check constraint: upper and lower bounds'
       if_main write(6,'(A,I0,A)')'    =========================================================='

       flag_conv = .true.
       ! 1. check parameter constraints 
       do j = 1, PINPT%nparam

         if(PINPT%slater_koster_type .gt. 10) then
           do i = 1, PINPT%param_nsub(j)
             k = nint(PINPT%param_const_nrl(1,i,j))
             if( k .eq. 0) then
               if(PINPT%param_nrl(i,j) .lt. PINPT%param_const_nrl(3,i,j)) then ! lower bound  
                 if( PINPT%param_nrl(i,j) .lt. 0d0 .and. i .eq. 4 .and. &
                    -PINPT%param_nrl(i,j) .ge. PINPT%param_const_nrl(2,i,j)) then 
                   PINPT%param_nrl(i,j) = PINPT%param_const_nrl(2,i,j)
                 elseif( PINPT%param_nrl(i,j) .le. 0d0 .and. i .eq. 4 .and. &
                    -PINPT%param_nrl(i,j) .lt. PINPT%param_const_nrl(2,i,j)) then
                   PINPT%param_nrl(i,j) = -PINPT%param_nrl(i,j)
                 else
                   PINPT%param_nrl(i,j) = PINPT%param_const_nrl(3,i,j)
                 endif
                 flag_conv = .false.
                 if_main write(6,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PINPT%param_name(j)),'(',i,') < ', &
                                                          PINPT%param_const_nrl(3,i,j),' --> ', &
                                                          trim(PINPT%param_name(j)),'(',i,') = ',PINPT%param_nrl(i,j)
               elseif(PINPT%param_nrl(i,j) .gt. PINPT%param_const_nrl(2,i,j)) then !upper bound
                 PINPT%param_nrl(i,j) = PINPT%param_const_nrl(2,i,j)
                 flag_conv = .false.
                 if_main write(6,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PINPT%param_name(j)),'(',i,') > ', &
                                                          PINPT%param_const_nrl(2,i,j),' --> ', &
                                                          trim(PINPT%param_name(j)),'(',i,') = ',PINPT%param_nrl(i,j)
               endif
             elseif( k .ge. 1) then ! if same as k-th parameter
               if(PINPT%param_nrl(i,j) .lt. PINPT%param_const_nrl(3,i, k )) then !lower bound
                 if(PINPT%param_nrl(i,j) .lt. 0d0 .and. i .eq. 4 .and. &
                    -PINPT%param_nrl(i,j) .ge. PINPT%param_const_nrl(2,i,k)) then
                   PINPT%param_nrl(i,j) = PINPT%param_const_nrl(2,i,k)
                 elseif( PINPT%param_nrl(i,j) .le. 0d0 .and. i .eq. 4 .and. &
                    -PINPT%param_nrl(i,j) .lt. PINPT%param_const_nrl(2,i,k)) then
                   PINPT%param_nrl(i,j) = -PINPT%param_nrl(i,j)
                 else
                   PINPT%param_nrl(i,j) = PINPT%param_const_nrl(3,i, k )
                 endif
                 flag_conv = .false.
                 if_main write(6,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PINPT%param_name(j)),'(',i,') < ', &
                                                          PINPT%param_const_nrl(3,i,k),' --> ', &
                                                          trim(PINPT%param_name(j)),'(',i,') = ',PINPT%param_nrl(i,j)
               elseif(PINPT%param_nrl(i,j) .gt. PINPT%param_const_nrl(2,i, k)) then !upper bound
                 PINPT%param_nrl(i,j) = PINPT%param_const_nrl(2,i, k)
                 flag_conv = .false.
                 if_main write(6,'(6x,2A,I0,A,F10.4,3A,I0,A,F10.4)')trim(PINPT%param_name(j)),'(',i,') > ', &
                                                          PINPT%param_const_nrl(2,i,k),' --> ', &
                                                          trim(PINPT%param_name(j)),'(',i,') = ',PINPT%param_nrl(i,j)
               endif
             endif
           enddo
         else
           if( nint(PINPT%param_const(1,j)) .eq. 0) then
             if(PINPT%param(j) .lt. PINPT%param_const(3,j)) then ! lower bound  
               PINPT%param(j) = PINPT%param_const(3,j)
               flag_conv = .false.
               if_main write(6,'(6x,2A,F10.4,3A,F10.4)')trim(PINPT%param_name(j)),' < ',PINPT%param_const(3,j),' --> ', &
                                                        trim(PINPT%param_name(j)),' = ',PINPT%param(j)
             elseif(PINPT%param(j) .gt. PINPT%param_const(2,j)) then !upper bound
               PINPT%param(j) = PINPT%param_const(2,j)
               flag_conv = .false.
               if_main write(6,'(6x,2A,F10.4,3A,F10.4)')trim(PINPT%param_name(j)),' > ',PINPT%param_const(2,j),' --> ', &
                                                        trim(PINPT%param_name(j)),' = ',PINPT%param(j)
             endif
           elseif(nint(PINPT%param_const(1,j)) .ge. 1) then
             if(PINPT%param(j) .lt. PINPT%param_const(3,nint(PINPT%param_const(1,j)))) then !lower bound
               PINPT%param(j) = PINPT%param_const(3,nint(PINPT%param_const(1,j)))
               flag_conv = .false.
               if_main write(6,'(6x,2A,F10.4,3A,F10.4)')trim(PINPT%param_name(j)),' < ',PINPT%param_const(3,nint(PINPT%param_const(1,j))),' --> ', &
                                                        trim(PINPT%param_name(j)),' = ',PINPT%param(j)
             elseif(PINPT%param(j) .gt. PINPT%param_const(2,nint(PINPT%param_const(1,j)))) then !upper bound
               PINPT%param(j) = PINPT%param_const(2,nint(PINPT%param_const(1,j)))
               flag_conv = .false.
               if_main write(6,'(6x,2A,F10.4,3A,F10.4)')trim(PINPT%param_name(j)),' > ',PINPT%param_const(2,nint(PINPT%param_const(1,j))),' --> ', &
                                                        trim(PINPT%param_name(j)),' = ',PINPT%param(j)
             endif
           endif
         endif
       enddo

       if(flag_conv) then
         if_main write(6,'(A)')'           Check fitness function updates '
         if_main write(6,'(A)')'    ============================================='
         ! 2. check fitness function updates
         fdiff = sqrt((fnorm_ - fnorm(1))**2)/fnorm(1)*100
         if( fdiff .le. PINPT%fdiff) then
           flag_conv = .true. 
           if_main write(6,'(A,F16.8,A,F16.8)')'  FDIFF(',PINPT%fdiff,')>=', fdiff
         elseif( (fdiff .gt. PINPT%fdiff ) ) then
           flag_conv = .false.
           if_main write(6,'(A,F16.8,A,F16.8)')'  FDIFF(',PINPT%fdiff,')<=', fdiff
         endif
       endif
       fnorm_ = fnorm(1)
 
       if_main write(6,'(A,I0,A)')'    =========================================================='
       if_main write(6,'(A)')' '
       if_main write(6,'(A,I0,A)')' **  END ',ifit,'-th LMDIF run'
       if_main write(6,'(A)')' '
       if(flag_conv .or. ifit .ge. mxfit) exit fit
     enddo fit

   ! run Genetic algorithm
   elseif(trim(PINPT%ls_type) .eq. 'PIKAIA' .or. trim(PINPT%ls_type) .eq. 'GA' .or. &
          trim(PINPT%ls_type) .eq. 'pikaia' .or. trim(PINPT%ls_type) .eq. 'ga' ) then
     if(PINPT%slater_koster_type .gt. 10) then
       if_main write(6,'(A)')'    !WARN! Current version does not support to use SK_SCALE_TYPE >= 10 along with LSTYPE=GA or PIKAIYA'
       if_main write(6,'(A)')'           Please set LSTYPE to LMDIF. The option will be available in the near future. Good luck...'
       kill_job
     else
       call gen_algo  ( get_eig, NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, EDFT%E, PGEOM%neig, &
                                 PINPT%init_erange, PINPT%nband, PWGHT, PINPT, PKAIA)
     endif
   endif
   ! print fitted parameters to file
   if(PINPT%flag_print_param) then
     if_main call print_param (PINPT, 0,PINPT%pfileoutnm, PINPT%flag_print_param)
     if_main call print_param (PINPT, 0, '   Fitted param(i):', .FALSE.)
     if_main write(6,'(A,A)')'  Fitted parameters will be written in ',PINPT%pfileoutnm
   endif

   if(allocated(ETBA_FIT%E)) deallocate(ETBA_FIT%E)
   if(allocated(ETBA_FIT%V)) deallocate(ETBA_FIT%V)

return
endsubroutine
