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

   if(PINPT%flag_sparse) then
     if_main write(6,'(A)') '    !WARN! EWINDOW tag cannot be used in FITTING procedures'
     if_main write(6,'(A)') '           Please deactivate EWINDOW for the further process.'
     if_main write(6,'(A)') '           Exit program...'
     stop
   endif

   !print target energy to file
   if_main call print_energy_weight( PKPTS%kpoint, PKPTS%nkpoint, EDFT, PWGHT, PGEOM%neig, &
                                     PINPT, 'band_structure_DFT.dat')

   !initial parameter
   if_main call print_param ( PINPT, 0, '  Initial param(i):', .FALSE.)

   ! fitting parameter
   if(trim(PINPT%ls_type) .eq. 'LMDIF' .or. trim(PINPT%ls_type) .eq. 'lmdif') then
     call leasqr_lm ( get_eig, NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, EDFT%E, PGEOM%neig, &
                               PINPT%init_erange, PINPT%nband, PWGHT, PINPT)
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
