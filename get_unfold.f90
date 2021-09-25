#include "alias.inc"
subroutine get_unfold(PINPT, PPRAM, PUFLD)
    use parameters
    use mpi_setup
    use print_io
    use time
    use mykind
    use set_default
    use projected_band
    implicit none
    integer(kind=sp)                           mpierr
    external                                   get_eig
    type(incar)                             :: PINPT, PINPT_UF
    type(kpoints)                           :: PKPTS_PBZ, PKPTS_SBZ
    type(poscar )                           :: PGEOM_PC, PGEOM_SC
    type(hopping)                           :: NN_TABLE_PC, NN_TABLE_SC
    type(params )                           :: PPRAM  
    type(unfold )                           :: PUFLD
    ! TEMPORARY TYPES-
    type(unfold )                           :: PUFLD_TMP
    type(weight)                            :: PWGHT
    type(energy)                            :: EDFT
    type(dos)                               :: PINPT_DOS
    type(berry)                             :: PINPT_BERRY
    type(gainp)                             :: PKAIA
    type(replot)                            :: PRPLT
    !------------------
    logical                                    flag_kline_mode
    real(kind=dp)                              b_SBZ(3,3)
    real(kind=dp)                              b_PBZ(3,3)
    real(kind=dp)                              kpoint_PBZ(3,PUFLD%unfold_nkpoint) ! PBZ k-points (cartesian)
    real(kind=dp)                              kline(PUFLD%unfold_nkpoint)        ! kline if kline_mode

    write(message,'( A)')' '  ; write_msg
    write(message,'( A)')' #===================================================='  ; write_msg
    write(message,'( A)')'   START BAND UNFOLDING ' ; write_msg
    write(message,'( A)')' #===================================================='  ; write_msg
    write(message,'( A)')' '  ; write_msg

    call time_check(t1, t0, 'init')
    call init_incar(PINPT_UF)
    PINPT_UF = PINPT
    PINPT_UF%nsystem = 2
    if(allocated (PINPT_UF%ifilenm)) deallocate(PINPT_UF%ifilenm)
    if(allocated (PINPT_UF%title  )) deallocate(PINPT_UF%title  )
    allocate(PINPT_UF%ifilenm(2))
    allocate(PINPT_UF%title  (2))
    PINPT_UF%ifilenm(1) = trim(PINPT%ifilenm(1))
    PINPT_UF%ifilenm(2) = trim(PUFLD%unfold_ifilenm_PBZ)
    PINPT_UF%title(1)   = 'SC'
    PINPT_UF%title(2)   = 'PC'
    
    ! initialize Primitive cell (PC) and Supercell (SC) info
    write(message,'(A)')' #---- SETUP SUPERCELL -----------'; write_msg
    call read_input(PINPT_UF,PPRAM,PKPTS_SBZ, PGEOM_SC,PWGHT,EDFT,NN_TABLE_SC, &
                    PINPT_DOS,PINPT_BERRY,PKAIA,PRPLT, PUFLD_TMP, 1)
    write(message,'(A)')' #---- SETUP PRIMITIVE CELL -----------'; write_msg
    call read_input(PINPT_UF,PPRAM,PKPTS_PBZ, PGEOM_PC,PWGHT,EDFT,NN_TABLE_PC, &
                    PINPT_DOS,PINPT_BERRY,PKAIA,PRPLT, PUFLD_TMP, 2)

    
!   call initialize_unfold(PINPT, PPRAM, PKPTS, PGEOM, NN_TABLE)
    
    

    write(6,*)"BBB UNFOLD"    
      
    return
endsubroutine
