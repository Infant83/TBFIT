#include "alias.inc"
subroutine get_symmetry_info(PGEOM)
   use parameters, only: poscar, onsite_tolerance, max_nsym, alphabet
   use spglib_interface
   use mpi_setup
   use do_math
   use iso_c_binding, only:  c_char, c_int, c_double, c_ptr, c_f_pointer
   implicit none
   type (poscar)       :: PGEOM       ! parameters for geometry info
   type(SpglibDataset) :: dset
   integer*4              i_dummy
   integer*4              i, isym
   real*8                 r(3)
   character*26        :: wyckoff_symbol
   character*128       :: fmt
   real*8                 a_latt_primitive(3,3)
   real*8                 a_coord_primitive(3,PGEOM%n_atom) 
   integer*4              spec_primitive(PGEOM%n_atom)
   integer*4              n_atom_primitive
   integer*4              std_mapping_to_primitive
   integer*4              trans_mat(3,3)
   integer*4              type_W
!  real*8                 R_WT(3)
   external               type_W !, R_WT
!  real*8, external   ::  R_WT(3)

   wyckoff_symbol(1:26)=alphabet

   if_main  write(6,'(A)')' '
   if_main  write(6,'(A)')'  :--------------------------:'
   if_main  write(6,'(A)')'  : SPACE GROUP INFORMATIONS :'
   if_main  write(6,'(A)')'  :----------------------------------------------------------------------------------------------:'
   if_main  write(6,'(A)')' '

   dset = spg_get_dataset(transpose(PGEOM%a_latt), PGEOM%a_coord, &
                          PGEOM%spec, PGEOM%n_atom, onsite_tolerance)

   a_latt_primitive  = transpose(PGEOM%a_latt)
   a_coord_primitive = PGEOM%a_coord
   spec_primitive    = PGEOM%spec

   n_atom_primitive  = spg_find_primitive(a_latt_primitive ,a_coord_primitive, &
                                          spec_primitive, &
                                          PGEOM%n_atom, onsite_tolerance)
   call get_crystal_system(PGEOM%spg_crystal_system,dset%spacegroup_number)
   PGEOM%spg_Hermann_Mauguin_number = spg_get_pointgroup( PGEOM%spg_point_group, &
                                      trans_mat, dset%rotations, dset%n_operations)
    ! Note: Returned point group number are mapped to the international table symbol
    !       and number as follows.
    !      1 :  1           2 :  -1          3 :  2           4 :  m     
    !      5 :  2/m         6 :  222         7 :  mm2         8 :  mmm   
    !      9 :  4           10:  -4          11:  4/m         12:  422   
    !      13:  4mm         14:  -42m        15:  4/mmm       16:  3     
    !      17:  -3          18:  32          19:  3m          20:  -3m   
    !      21:  6           22:  -6          23:  6/m         24:  622   
    !      25:  6mm         26:  -62m        27:  6/mmm       28:  23    
    !      29:  m-3         30:  432         31:  -43m        32:  m-3m  
     
   allocate(PGEOM%spg_a_coord_primitive(3,n_atom_primitive))
   allocate(PGEOM%spg_spec_primitive   (  n_atom_primitive))
   allocate(PGEOM%spg_rotations(3,3,dset%n_operations))   
   allocate(PGEOM%spg_translations(3,dset%n_operations))
   allocate(PGEOM%spg_wyckoffs(PGEOM%n_atom))
   allocate(PGEOM%spg_equivalent_atoms(PGEOM%n_atom))
   allocate(PGEOM%spg_det_W(dset%n_operations))
   allocate(PGEOM%spg_tr_W(dset%n_operations))
   allocate(PGEOM%spg_type_W(dset%n_operations))
   allocate(PGEOM%spg_a_coord_operated(3,PGEOM%n_atom,dset%n_operations))

   PGEOM%spg_error                 = dset%spglib_error
   PGEOM%spg_space_group           = dset%spacegroup_number
   PGEOM%spg_hall_number           = dset%hall_number
   PGEOM%spg_transformation_matrix = dset%transformation_matrix
   PGEOM%spg_origin_shift          = dset%origin_shift
   PGEOM%spg_n_operations          = dset%n_operations
   PGEOM%spg_international         = trim(dset%international_symbol)
   PGEOM%spg_hall_symbol           = trim(dset%hall_symbol)
   PGEOM%spg_choice                = trim(dset%choice)
   PGEOM%spg_rotations             = dset%rotations
   PGEOM%spg_translations          = dset%translations
   PGEOM%spg_wyckoffs              = dset%wyckoffs
   PGEOM%spg_equivalent_atoms      = dset%equivalent_atoms
   PGEOM%spg_a_latt_primitive      = transpose(a_latt_primitive)
   PGEOM%spg_a_coord_primitive     = a_coord_primitive(:,1:n_atom_primitive)
   PGEOM%spg_spec_primitive        = spec_primitive(1:n_atom_primitive)
   PGEOM%spg_n_atom_primitive      = n_atom_primitive

   do i = 1, PGEOM%spg_n_operations
      PGEOM%spg_det_W(i) = determinant3i(PGEOM%spg_rotations(:,:,i))
      PGEOM%spg_tr_W(i) = tracei(PGEOM%spg_rotations(:,:,i),3)
      PGEOM%spg_type_W(i) = type_W(PGEOM%spg_det_W(i), PGEOM%spg_tr_W(i))
   enddo

   if (PGEOM%spg_space_group .eq. 0) then
     if_main stop "  Space group could not be found. Exit..."
   else
     PGEOM%spg_schoenflies = ' '
     i_dummy = spg_get_schoenflies(PGEOM%spg_schoenflies, transpose(PGEOM%a_latt), &
                                   PGEOM%a_coord, PGEOM%spec, PGEOM%n_atom, onsite_tolerance)
   endif
   
   if_main_then
   ! write general information
   write(6,'(A,A )'   )'   - Crystal system             : ', trim(PGEOM%spg_crystal_system)
   write(6,'(A,I3)'   )'   - Space group number         : ', PGEOM%spg_space_group
   write(6,'(A,A )'   )'   - Crystal choice             : ', trim(PGEOM%spg_choice)
   write(6,'(A,A)'    )'   - International short symbol : ', trim(PGEOM%spg_international)
   write(6,'(A,A )'   )'   - Point group symbol         : ', trim(PGEOM%spg_point_group)
   write(6,'(A,A)'    )'   - Schoenflies symbol         : ', trim(PGEOM%spg_schoenflies)
   write(6,'(3A,I3,A)')'   - Hall symbol [number]       : ', trim(PGEOM%spg_hall_symbol), &
                                                        ' [', PGEOM%spg_hall_number,']'
   write(6,'(A)'      )'   '

   ! write Wyckoff symbol and equivalent atoms
   ! Note: The Wyckoff symbol is determined using "Coordinates" in the Wyckoff position data set,
   !       and the "Coordinates (direct/fractional)" are listed in the 
   !       'T. Hahn, International Tables for Crystallography Volume A (Wiley, 2011)',
   !       with respect to the "Primitive" lattice vector which is found by "spg_find_primitive" 
   !       function of SPGLIB library.
   !       The details can be found in following reference: https://arxiv.org/pdf/1808.01590.pdf
   !        [Ref] A. Togo and I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv.1808.01590 (2018)
   write(6,'(A       )')"    :--------------* Original lattice vector A --------------------------------------------------:" 
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A,3F15.8)')"    :   Ai = [Ai1,Ai2,Ai3]         [A11 A12 A13] ",PGEOM%a_latt(:,1)
   write(6,'(A,3F15.8)')"   -:   Ai'= transpose(Ai)  => A = [A21 A22 A23]=",PGEOM%a_latt(:,2)
   write(6,'(A,3F15.8)')"    : =>A' = [A1',A2',A3']         [A31 A32 A33] ",PGEOM%a_latt(:,3)
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :--------------* Primitive lattice vector B -------------------------------------------------:" 
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A,3F15.8)')"    :   Bi = [Bi1,Bi2,Bi3]         [B11 B12 B13] ",PGEOM%spg_a_latt_primitive(:,1)
   write(6,'(A,3F15.8)')"   -:   Bi'= transpose(Bi)  => B = [B21 B22 B23]=",PGEOM%spg_a_latt_primitive(:,2)
   write(6,'(A,3F15.8)')"    : =>B' = [B1',B2',B3']         [B31 B32 B33] ",PGEOM%spg_a_latt_primitive(:,3)
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :--------------* Transformation matrix P ----------------------------------------------------:" 
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A,3F15.8)')"    :  B'*P = A' or P'*B = A       [P11 P12 P13] ",PGEOM%spg_transformation_matrix(:,1)
   write(6,'(A,3F15.8)')"   -:  or A'*inv(P) = B'    => P = [P21 P22 P23]=",PGEOM%spg_transformation_matrix(:,2)
   write(6,'(A,3F15.8)')"    :  or inv(P')*A = B            [P31 P32 P33] ",PGEOM%spg_transformation_matrix(:,3)
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :--------------* Origin shift O -------------------------------------------------------------:" 
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A,3F15.8)')"   -: O = Origin_A-Origin_B => O = [O1  O2  O3 ]=",PGEOM%spg_origin_shift(:)
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :--------------* Coordinate transformaton by P and O ----------------------------------------:" 
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :    Rs = atomic position w.r.t. primitive lattice vector (with origin shift)"  
   write(6,'(A       )')"   -:    R  = atomic position w.r.t. current (A) lattice vector "
   write(6,'(A       )')"    :    Rs'= P * R' + O'   "
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :--------------* Wyckoff symbol (Ws) and Equivalent atom (Eqv.) -----------------------------:" 
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A)'       )"    : ATOM Eqv.   : W :            R (input geom)           :           Rs (primitive geom)       "
   write(6,'(A       )')"   -: -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -:" 
   write(fmt,'(A)')"(A,I5,1x,A,I5,1x,A,A,A,3F12.8,A,3F12.8)"
   do i = 1, PGEOM%n_atom
     write(6,fmt)"    :", i,'=',PGEOM%spg_equivalent_atoms(i)+1, ': ', &
                          wyckoff_symbol(PGEOM%spg_wyckoffs(i)+1:PGEOM%spg_wyckoffs(i)+1), &
                    ' :', PGEOM%a_coord(:,i),' :', PGEOM%spg_a_coord_primitive(:,i)
   enddo
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :--------------------------------------------------------------------------------------------:" 

   ! write symmetry operations 
   write(6,'(A)'       )'    '
   write(6,'(A       )')"    :--------------* Symmetry operations {W|T} --------------------------------------------------:" 
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A,I3)'    )"    :  Total number of symmetry operation (nsym) = ",PGEOM%spg_n_operations                        
   write(6,'(A)'       )"    :                                                                                            :"
   write(6,'(A       )')"    :--------------------------------------------------------------------------------------------:" 
   write(6,'(A)'       )'    '
   if_main_end
   do isym = 1, PGEOM%spg_n_operations
     if_main write(6,'(A,I3,3(A,I3)  )')'       - operation # ',isym, ' : det(W) = ',PGEOM%spg_det_W(isym), &
                                                           ' , tr(W) = ', PGEOM%spg_tr_W(isym), &
                                                           ' => type(W) = ', PGEOM%spg_type_W(isym)
                                                      ! /* Look-up table */
                                                      ! /* Operation   -6 -4 -3 -2 -1  1  2  3  4  6 */
                                                      ! /* Trace     -  2 -1  0  1 -3  3 -1  0  1  2 */
                                                      ! /* Determinant -1 -1 -1 -1 -1  1  1  1  1  1 */
     if_main write(6,97)PGEOM%spg_rotations(:,1,isym), mod(PGEOM%spg_translations(1,isym)+10d0,1d0)
     if_main write(6,98)PGEOM%spg_rotations(:,2,isym), mod(PGEOM%spg_translations(2,isym)+10d0,1d0)
     if_main write(6,99)PGEOM%spg_rotations(:,3,isym), mod(PGEOM%spg_translations(3,isym)+10d0,1d0)

     call get_operated_coord(PGEOM%spg_rotations(:,:,isym), PGEOM%spg_translations(:,isym), &
                             PGEOM%a_coord, PGEOM%spg_a_coord_operated(:,:,isym), &
                             PGEOM%spg_n_operations, PGEOM%n_atom)
     if_main_then
     do i = 1, PGEOM%n_atom
       write(6,'(A,I5,A,3F10.5,A,3F10.5)')"           atom",i,':',mod(PGEOM%a_coord(:,i)+10d0,1d0), &
                                                           ' -> ',PGEOM%spg_a_coord_operated(:,i,isym)
     enddo
     if_main_end
   enddo
   if_main write(6,'(A       )')"    :--------------------------------------------------------------------------------------------:" 
   if_main write(6,'(A)')' '


97 format('           [',I3,1x,I3,1x,I3,']   [x]','   [',F8.4,']')
98 format('           [',I3,1x,I3,1x,I3,'] * [y]',' + [',F8.4,']')
99 format('           [',I3,1x,I3,1x,I3,']   [z]','   [',F8.4,']')
return
endsubroutine

subroutine get_operated_coord(rot, t, R, R_, nsym, natom)
   implicit none
   integer*4   nsym, natom
   integer*4   rot(3,3)
   real*8      t(3)
   real*8      R(3,natom)
   real*8      R_(3,natom)
   integer*4   i

     do i = 1, natom
       R_(:,i) = mod( matmul( transpose(rot(:,:)), R(:,i) ) + t(:) + 10d0,1d0)
     enddo

return
endsubroutine

function type_W(det,tr)
   implicit none
   integer*4   type_W
   integer*4   det
   integer*4   tr

   select case (det)

     case(1)
       select case(tr)
         case( 3)
           type_W = 1
         case(-1)
           type_W = 2
         case( 0)
           type_W = 3
         case( 1)
           type_W = 4
         case( 2)
           type_W = 6
       end select

     case(-1)
       select case(tr)
         case(-3)
           type_W =-1
         case( 1)
           type_W =-2
         case( 0)
           type_W =-3
         case(-1)
           type_W =-4
         case(-2)
           type_W =-6
       end select

   end select

return
endfunction


subroutine  get_crystal_system(crystal_system,space_group) !result(crystal_system)
   implicit none
   integer*4    space_group
   character*11 international
   character*12 crystal_system

   crystal_system(1:12) = '            '   

   select case (space_group)

     case(1:2)
       crystal_system='triclinic'

     case(3:15)
       crystal_system='monoclinic'

     case(16:74)
       crystal_system='orthorhombic'

     case(75:142)
       crystal_system='tetragonal'  

     case(143:167)
!      if(international(1:1) .eq. 'P')
!        lattice_system='tetragonal'
       crystal_system='trigonal'    

     case(168:194)
       crystal_system='hexagonal'  
   
     case(195:230)
       crystal_system='cubic' 
   
   end select

   return
endsubroutine
