subroutine get_fit(PINPT, PKPTS, EDFT, PWGHT, PGEOM, NN_TABLE)
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

   !print target energy to file
   if(myid .eq. 0) call print_energy_weight( PKPTS%kpoint, PKPTS%nkpoint, EDFT, PWGHT, PGEOM%neig, PINPT, 'band_structure_DFT.dat')

   !initial parameter
   if(myid .eq. 0) call print_param ( PINPT, 0, '  Initial param(i):', .FALSE.)

   ! fitting parameter
   call leasqr_lm ( get_eig, NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, EDFT%E, PGEOM%neig, PWGHT, PINPT)

   ! print fitted parameters to file
   if(PINPT%flag_print_param) then
     if(myid .eq. 0) call print_param (PINPT, 0,PINPT%pfileoutnm, PINPT%flag_print_param)
     if(myid .eq. 0) write(6,'(A,A)')'  Fitted parameters will be written in ',PINPT%pfileoutnm
   endif

return
endsubroutine
