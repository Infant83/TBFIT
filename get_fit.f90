#include "alias.inc"
subroutine get_fit(PINPT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE, PKAIA)
   use parameters
   use mpi_setup
   implicit none
   integer*4         mpierr
   external          get_eig
   type (incar)   :: PINPT       ! parameters for input arguments
   type (dos)     :: PINPT_DOS   ! parameters for density of states calculation
   type (poscar)  :: PGEOM       ! parameters for geometry info
   type (kpoints) :: PKPTS       ! parameters for kpoints
   type (energy)  :: EDFT        ! target energy to be fitted to
   type (weight)  :: PWGHT       ! weight factor for the fitting to target energy
   type (hopping) :: NN_TABLE    ! table for hopping index
   type (gainp)   :: PKAIA       ! input/control parameters for genetic algorithm
   logical           flag_conv
   integer*4         j,ifit,mxfit

   if(PINPT%flag_sparse) then
     if_main write(6,'(A)') '    !WARN! EWINDOW tag cannot be used in FITTING procedures'
     if_main write(6,'(A)') '           Please deactivate EWINDOW for the further process.'
     if_main write(6,'(A)') '           Exit program...'
     stop
   endif

   flag_conv = .false.
   mxfit = PINPT%mxfit
   if(mxfit .le. 0) then
     if_main write(6,'(A)')'    !WARN! The current MXFIT value, maximum number of LMDIF attempt, is <= 0 .'
     if_main write(6,'(A)')'           Please make sure that MXFIT > 0 in your input tag.'
     kill_job
   endif
   ifit = 0

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

       call leasqr_lm ( get_eig, NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, EDFT%E, PGEOM%neig, &
                                 PINPT%init_erange, PINPT%nband, PWGHT, PINPT)

       if_main write(6,'(A,I0,A)')'           Check constraint: upper and lower bounds'
       if_main write(6,'(A,I0,A)')'    =========================================================='

       flag_conv = .true.
       do j = 1, PINPT%nparam
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
       enddo

       if_main write(6,'(A,I0,A)')'    =========================================================='
       if_main write(6,'(A)')' '
       if_main write(6,'(A,I0,A)')' **  END ',ifit,'-th LMDIF run'
       if_main write(6,'(A)')' '
       if(flag_conv .or. ifit .eq. mxfit) exit fit
     enddo fit

   elseif(trim(PINPT%ls_type) .eq. 'PIKAIA' .or. trim(PINPT%ls_type) .eq. 'GA' .or. &
          trim(PINPT%ls_type) .eq. 'pikaia' .or. trim(PINPT%ls_type) .eq. 'ga' ) then
     call gen_algo  ( get_eig, NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, EDFT%E, PGEOM%neig, &
                               PINPT%init_erange, PINPT%nband, PWGHT, PINPT, PKAIA)
   endif
   ! print fitted parameters to file
   if(PINPT%flag_print_param) then
     if_main call print_param (PINPT, 0,PINPT%pfileoutnm, PINPT%flag_print_param)
     if_main call print_param (PINPT, 0, '   Fitted param(i):', .FALSE.)
     if_main write(6,'(A,A)')'  Fitted parameters will be written in ',PINPT%pfileoutnm
   endif

return
endsubroutine
