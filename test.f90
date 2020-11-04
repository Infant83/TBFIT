#include "alias.inc"
subroutine test(PINPT)
  use print_io
  use parameters
  implicit none
  type (incar)   :: PINPT       ! parameters for input arguments
  type (params)  :: PPRAM       ! parameters for input arguments
  type (dos)     :: PINPT_DOS   ! parameters for density of states calculation
  type (berry)   :: PINPT_BERRY ! parameters for berry phase related quantity calculation
  type (poscar)  :: PGEOM       ! parameters for geometry info
  type (kpoints) :: PKPTS       ! parameters for kpoints
  type (energy)  :: EDFT        ! target energy to be fitted to
  type (energy)  :: ETBA        ! calculated energy with given TB-parameters
  type (weight)  :: PWGHT       ! weight factor for the fitting to target energy
  type (hopping) :: NN_TABLE    ! table for hopping index
  type (gainp)   :: PKAIA       ! input/control parameters for genetic algorithm
  type (replot)  :: PRPLT
  character*132     filenmin

  call read_input(PINPT,PPRAM,PINPT_DOS,PINPT_BERRY,PKPTS,PGEOM,PWGHT,EDFT,NN_TABLE,PKAIA,PRPLT, 1)

      call allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA)
      call get_eig(NN_TABLE, PKPTS%kpoint, PKPTS%nkpoint, PINPT, PPRAM, ETBA%E, ETBA%V, ETBA%SV, PGEOM%neig, &
                   PGEOM%init_erange, PGEOM%nband, PINPT%flag_get_orbital, PINPT%flag_sparse, .true., PINPT%flag_phase)

      call print_band_structure(PKPTS, ETBA, EDFT, PGEOM, PINPT, PPRAM, PWGHT)

  stop 
endsubroutine
