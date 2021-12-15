#include "alias.inc"
! Module for unfolding band structure
! written by: 
!   Hyun-Jung Kim (FZJ, h.kim@fz-juelich.de)
!   Jejune Park   (KIAS, jejunepark@kias.re.kr)
! Ref: C.-C. Lee, Y. Yamada-Takamura, and T. Ozaki, 
!      J. Phys.: Condens. Matter 25, 345501 (2013)
!


! NOTE: This module is under developing now...
module unfold
    use do_math
    use print_io
    use mpi_setup
    use mykind

contains

    subroutine set_kpoints_PBZ(PINPT, PKPTS)
        use parameters, only: incar, kpoints
        implicit none
        type(incar   )  :: PINPT
        type(kpoints )  :: PKPTS
        real(kind=dp)      b_PBZ(3,3)
        real(kind=dp)      b_SBZ(3,3)

        
        
        return
    endsubroutine

endmodule
