#include "alias.inc"
#ifdef SPGLIB
subroutine get_symmetry_info(PGEOM)
   use parameters, only: poscar, onsite_tolerance, max_nsym, alphabet
   use spglib_interface
   use mpi_setup
   use do_math
   use iso_c_binding, only:  c_char, c_int, c_double, c_ptr, c_f_pointer
   use print_io
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

   write(message,'(A)')' '  ; write_msg
   write(message,'(A)')'  :--------------------------:'  ; write_msg
   write(message,'(A)')'  : SPACE GROUP INFORMATIONS :'  ; write_msg
   write(message,'(A)')'  :----------------------------------------------------------------------------------------------:'  ; write_msg
   write(message,'(A)')' '  ; write_msg

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
     if_main_then
       write(message,'(A)') "  Space group could not be found. Exit..."  ; write_msg
     if_main_end
   else
     PGEOM%spg_schoenflies = ' '
     i_dummy = spg_get_schoenflies(PGEOM%spg_schoenflies, transpose(PGEOM%a_latt), &
                                   PGEOM%a_coord, PGEOM%spec, PGEOM%n_atom, onsite_tolerance)
   endif

   if_main_then
   ! write general information
   write(message,'(A,A )'   )'   - Crystal system             : ', trim(PGEOM%spg_crystal_system)  ; write_msg
   write(message,'(A,I3)'   )'   - Space group number         : ', PGEOM%spg_space_group  ; write_msg
   write(message,'(A,A )'   )'   - Crystal choice             : ', trim(PGEOM%spg_choice)  ; write_msg
   write(message,'(A,A)'    )'   - International short symbol : ', trim(PGEOM%spg_international)  ; write_msg
!  write(message,'(A,A)'    )'   - Point group symbol         : ', trim(PGEOM%spg_point_group)  ; write_msg
!  write(message,'(A,A)'    )'   - Schoenflies symbol         : ', trim(PGEOM%spg_schoenflies)  ; write_msg
   write(message,'(3A,I3,A)')'   - Hall symbol [number]       : ', trim(PGEOM%spg_hall_symbol), & 
                                                        ' [', PGEOM%spg_hall_number,']' ; write_msg
   write(message,'(A)'      )'   '  ; write_msg

   ! write Wyckoff symbol and equivalent atoms
   ! Note: The Wyckoff symbol is determined using "Coordinates" in the Wyckoff position data set,
   !       and the "Coordinates (direct/fractional)" are listed in the 
   !       'T. Hahn, International Tables for Crystallography Volume A (Wiley, 2011)',
   !       with respect to the "Primitive" lattice vector which is found by "spg_find_primitive" 
   !       function of SPGLIB library.
   !       The details can be found in following reference: https://arxiv.org/pdf/1808.01590.pdf
   !        [Ref] A. Togo and I. Tanaka, "Spglib: a software library for crystal symmetry search", arXiv.1808.01590 (2018)
   write(message,'(A       )')"    :--------------* Original lattice vector A --------------------------------------------------:"   ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A,3F15.8)')"    :   Ai = [Ai1,Ai2,Ai3]         [A11 A12 A13] ",PGEOM%a_latt(:,1)  ; write_msg
   write(message,'(A,3F15.8)')"   -:   Ai'= transpose(Ai)  => A = [A21 A22 A23]=",PGEOM%a_latt(:,2)  ; write_msg
   write(message,'(A,3F15.8)')"    : =>A' = [A1',A2',A3']         [A31 A32 A33] ",PGEOM%a_latt(:,3)  ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :--------------* Primitive lattice vector B -------------------------------------------------:"   ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A,3F15.8)')"    :   Bi = [Bi1,Bi2,Bi3]         [B11 B12 B13] ",PGEOM%spg_a_latt_primitive(:,1)  ; write_msg
   write(message,'(A,3F15.8)')"   -:   Bi'= transpose(Bi)  => B = [B21 B22 B23]=",PGEOM%spg_a_latt_primitive(:,2)  ; write_msg
   write(message,'(A,3F15.8)')"    : =>B' = [B1',B2',B3']         [B31 B32 B33] ",PGEOM%spg_a_latt_primitive(:,3)  ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :--------------* Transformation matrix P ----------------------------------------------------:"   ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A,3F15.8)')"    :  B'*P = A' or P'*B = A       [P11 P12 P13] ",PGEOM%spg_transformation_matrix(:,1)  ; write_msg
   write(message,'(A,3F15.8)')"   -:  or A'*inv(P) = B'    => P = [P21 P22 P23]=",PGEOM%spg_transformation_matrix(:,2)  ; write_msg
   write(message,'(A,3F15.8)')"    :  or inv(P')*A = B            [P31 P32 P33] ",PGEOM%spg_transformation_matrix(:,3)  ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :--------------* Origin shift O -------------------------------------------------------------:"   ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A,3F15.8)')"   -: O = Origin_A-Origin_B => O = [O1  O2  O3 ]=",PGEOM%spg_origin_shift(:)  ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :--------------* Coordinate transformaton by P and O ----------------------------------------:"   ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :    Rs = atomic position w.r.t. primitive lattice vector (with origin shift)"    ; write_msg
   write(message,'(A       )')"   -:    R  = atomic position w.r.t. current (A) lattice vector "  ; write_msg
   write(message,'(A       )')"    :    Rs'= P * R' + O'   "  ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :--------------* Wyckoff symbol (Ws) and Equivalent atom (Eqv.) -----------------------------:"   ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A)'       )"    : ATOM Eqv.   : W :            R (input geom)           :           Rs (primitive geom)       "  ; write_msg
   write(message,'(A       )')"   -: -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -:"   ; write_msg
   write(fmt,'(A)')"(A,I5,1x,A,I5,1x,A,A,A,3F12.8,A,3F12.8)"
   do i = 1, PGEOM%n_atom
     write(message,fmt)"    :", i,'=',PGEOM%spg_equivalent_atoms(i)+1, ': ', wyckoff_symbol(PGEOM%spg_wyckoffs(i)+1:PGEOM%spg_wyckoffs(i)+1), ' :', PGEOM%a_coord(:,i),' :', PGEOM%spg_a_coord_primitive(:,i)  ; write_msg
   enddo
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :--------------------------------------------------------------------------------------------:"   ; write_msg

   ! write symmetry operations 
   write(message,'(A)'       )'    '  ; write_msg
   write(message,'(A       )')"    :--------------* Symmetry operations {W|T} --------------------------------------------------:"   ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A,I3)'    )"    :  Total number of symmetry operation (nsym) = ",PGEOM%spg_n_operations                          ; write_msg
   write(message,'(A)'       )"    :                                                                                            :"  ; write_msg
   write(message,'(A       )')"    :--------------------------------------------------------------------------------------------:"   ; write_msg
   write(message,'(A)'       )'    '  ; write_msg
   if_main_end
   do isym = 1, PGEOM%spg_n_operations
     write(message,'(A,I3,3(A,I3)  )')'       - operation # ',isym, ' : det(W) = ',PGEOM%spg_det_W(isym), ' , tr(W) = ', PGEOM%spg_tr_W(isym), ' => type(W) = ', PGEOM%spg_type_W(isym)  ; write_msg
                                                      ! /* Look-up table */
                                                      ! /* Operation   -6 -4 -3 -2 -1  1  2  3  4  6 */
                                                      ! /* Trace     -  2 -1  0  1 -3  3 -1  0  1  2 */
                                                      ! /* Determinant -1 -1 -1 -1 -1  1  1  1  1  1 */
     write(message,97)PGEOM%spg_rotations(:,1,isym), mod(PGEOM%spg_translations(1,isym)+10d0,1d0)  ; write_msg
     write(message,98)PGEOM%spg_rotations(:,2,isym), mod(PGEOM%spg_translations(2,isym)+10d0,1d0)  ; write_msg
     write(message,99)PGEOM%spg_rotations(:,3,isym), mod(PGEOM%spg_translations(3,isym)+10d0,1d0)  ; write_msg

     call get_operated_coord(PGEOM%spg_rotations(:,:,isym), PGEOM%spg_translations(:,isym), &
                             PGEOM%a_coord, PGEOM%spg_a_coord_operated(:,:,isym), &
                             PGEOM%spg_n_operations, PGEOM%n_atom)
     if_main_then
     do i = 1, PGEOM%n_atom
       write(message,'(A,I5,A,3F10.5,A,3F10.5)')"           atom",i,':',mod(PGEOM%a_coord(:,i)+10d0,1d0), ' -> ',PGEOM%spg_a_coord_operated(:,i,isym)  ; write_msg
     enddo
     if_main_end
   enddo
   write(message,'(A       )')"    :--------------------------------------------------------------------------------------------:"   ; write_msg
   write(message,'(A)')' '  ; write_msg


97 format('           [',I3,1x,I3,1x,I3,']   [x]','   [',F8.4,']')
98 format('           [',I3,1x,I3,1x,I3,'] * [y]',' + [',F8.4,']')
99 format('           [',I3,1x,I3,1x,I3,']   [z]','   [',F8.4,']')
return
endsubroutine
#endif

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

subroutine set_symmetry_operator(S, ROT, O, theta, PGEOM, PINPT)
   use parameters, only: poscar, incar, onsite_tolerance, zi, pi
   use print_matrix
   use mpi_setup
   use do_math
   use element_info
   implicit none
   type (incar)   :: PINPT       ! parameters for input arguments
   type (poscar)  :: PGEOM       ! parameters for geometry info
   real*8            ROT(3,3)    ! rotation matrix (rotation in fractional coord)
   real*8            O(3)        ! origin
   complex*16        S(PGEOM%neig*PINPT%ispinor, PGEOM%neig*PINPT%ispinor) ! symmetry operator for each orbital
   complex*16,allocatable :: S_block(:,:)
   real*8            a1(3), a2(3), a3(3)
   real*8            R(3), R_(3) ! R and R' where R' = ROT * R
   real*8            pos_i(3), pos_j(3) ! cartesian position for i(j) orbital
   real*8            shift(3)
   real*8            theta, r0(3)
   integer*4         iorb, jorb, imatrix, jmatrix
   integer*4         iinit, iend, jinit, jend
   integer*4         neig, niorb, njorb
   integer*4         ix, iy, iz, i, j, is
   integer*4         ispec, jspec
   character*8       orb_name
   integer*4         l, mpierr
   real*8,external:: enorm

   neig = PGEOM%neig
   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   r0=(/0d0,0d0,0d0/)
   S = 0d0

ii:do i = 1, PGEOM%n_atom
     niorb = PGEOM%n_orbital(i)
     if(niorb .eq. 0) cycle ii
     R = PGEOM%a_coord(:,i) - O
     pos_i = R(1)*a1 + R(2)*a2 + R(3)*a3
     ispec = atomic_number(PGEOM%c_spec(PGEOM%spec(i)))
     imatrix = sum( PGEOM%n_orbital(1:i) ) - PGEOM%n_orbital(i) + 1
     do ix = -2,2
     do iy = -2,2
     do iz = -2,2
    jj:do j = 1, PGEOM%n_atom
         njorb = PGEOM%n_orbital(j)
         jspec = atomic_number(PGEOM%c_spec(PGEOM%spec(j)))
         if(njorb .eq. 0 .and. PGEOM%spec(i) .eq. PGEOM%spec(j) ) cycle jj

         R_ = PGEOM%a_coord(:,j) - O
         R_ = matmul( transpose(ROT), R_ )
         pos_j= R_(1)*a1 + R_(2)*a2 + R_(3)*a3
         shift= ix*a1+iy*a2+iz*a3
         
         if( enorm(3, pos_i - pos_j - shift ) .lt. onsite_tolerance) then
           ! equivalent (atomic position) under rotation ROT
           allocate( S_block( njorb, njorb ) )
           call set_symmetry_block(S_block, njorb, PGEOM%c_orbital(:,1: njorb), theta)
           jmatrix = sum( PGEOM%n_orbital(1:j) ) - PGEOM%n_orbital(j) + 1
           if(PINPT%ispinor .eq. 2) then
             iinit=imatrix ; iend=imatrix+niorb-1 ; jinit=jmatrix ; jend=jmatrix+njorb-1
             S(iinit:iend, jinit:jend) = S_block * exp( zi*theta/180d0*pi/2d0) !spin-up

             iinit=imatrix+neig ; iend=imatrix+niorb-1+neig ; jinit=jmatrix+neig ; jend=jmatrix+njorb-1+neig
             S(iinit:iend, jinit:jend) = S_block * exp(-zi*theta/180d0*pi/2d0) !spin-dn
           elseif(PINPT%ispinor .eq. 1) then
             iinit=imatrix ; iend=imatrix+niorb-1 ; jinit=jmatrix ; jend=jmatrix+njorb-1
             S(iinit:iend, jinit:jend) = S_block                               !spin-up
           endif
  
          deallocate(S_block)
         endif
       enddo jj
     enddo
     enddo
     enddo     
   enddo ii  

   ! call print_matrix_c(S, PGEOM%neig*PINPT%ispinor, PGEOM%neig*PINPT%ispinor, 'SYMM', 0, 'F6.3')

   return
endsubroutine
subroutine set_symmetry_block(S_block, norb, corb, theta)
   use parameters, only: pi
   use mpi_setup
   use print_matrix
   implicit none
   integer*4    i
   integer*4    iorb
   integer*4    iporb(3), idorb(5)
   integer*4    norb, mpierr
   integer*4    p_order, d_order
   integer*4    ispinor
   real*8       theta
   complex*16   S_block(norb, norb)
   character*8  corb(norb), corb_, porb_
   character*8  orb_s, orb_px, orb_py, orb_pz
   character*8  orb_dx2, orb_dxy, orb_dz2, orb_dxz, orb_dyz
   real*8       c, s
   logical      flag_skip_p, flag_skip_s, flag_skip_d

   S_block = 0d0
   p_order = 0
   d_order = 0
   c = cos(theta / 180d0 * pi)
   s = sin(theta / 180d0 * pi)
   flag_skip_s = .false.
   flag_skip_p = .false.
   flag_skip_d = .false.
   
   do iorb=1,norb
     corb_ = trim(corb(iorb))
     if( corb_(1:1) .eq. 'p' .and. .not. flag_skip_p) then
       do i = 1, norb
         if(corb(i) .eq. 'px') then
           iporb(1) = i
         elseif(corb(i) .eq. 'py') then
           iporb(2) = i
         elseif(corb(i) .eq. 'pz') then
           iporb(3) = i
         endif
       enddo     
       flag_skip_p = .true.
     elseif( corb_(1:1) .eq. 'd' .and. .not. flag_skip_d) then
       do i = 1, norb
         if(corb(i) .eq. 'dx2') then
           idorb(1) = i
         elseif(corb(i) .eq. 'dxy') then
           idorb(2) = i
         elseif(corb(i) .eq. 'dz2') then
           idorb(3) = i
         elseif(corb(i) .eq. 'dxz') then
           idorb(4) = i
         elseif(corb(i) .eq. 'dyz') then
           idorb(5) = i
         endif
       enddo
       flag_skip_d = .true.
     endif
   enddo
orb:do iorb=1, norb
      corb_ = trim(corb(iorb))
      if( corb_(1:1) .eq. 's') then
        S_block(iorb, iorb) = 1d0
      elseif( corb_(1:2) .eq. 'px') then
        S_block(iporb(1), iporb) = (/   c,   s, 0d0/)
      elseif( corb_(1:2) .eq. 'py') then
        S_block(iporb(2), iporb) = (/  -s,   c, 0d0/)
      elseif( corb_(1:2) .eq. 'pz') then
        S_block(iporb(3), iporb) = (/ 0d0, 0d0, 1d0/)
      endif

    enddo orb

 !call print_matrix_c(S_block, norb, norb, 'SSSS', 0, 'F9.4')

!kill_job

   return
endsubroutine
subroutine get_symmetry_matrix(Sij, S_eig, V, phase_shift, S_OP, E, PGEOM,  neig, ispinor, flag_phase_shift)
   use parameters, only: poscar, eta
   use do_math
   use print_matrix
   use berry_phase
   implicit none
   type (poscar)  :: PGEOM       ! parameters for geometry info
   integer*4    neig, ispinor
   integer*4    ie, je, init,iend
   complex*16   V(neig*ispinor,neig*ispinor)
   complex*16   S_(neig*ispinor,neig*ispinor)
   complex*16   S_OP(neig*ispinor,neig*ispinor)
   complex*16   Sij(neig*ispinor,neig*ispinor)
   complex*16   S_eig(neig*ispinor)
   complex*16   S_eig_(neig*ispinor)
   complex*16   phase_shift(neig*ispinor)
   real*8       E(neig*ispinor)
   logical      flag_phase_shift 
   real*8       very_small

   very_small = 0.0000001d0
   init   = 1
   Sij = 0d0

   do ie = 1, neig*ispinor
     do je = 1, neig*ispinor
      Sij(ie,je) = dot_product( V(:,ie), matmul(S_OP, phase_shift*V(:,je) ) )
      !if(flag_phase_shift) then
      !  Sij(ie,je) = dot_product( V(:,ie), matmul(S_OP,      phase_shift  *   V(:,je) ) )
      !else
      !  Sij(ie,je) = dot_product( V(:,ie), matmul(S_OP, V(:,je) ) )
      !endif
     enddo
   enddo

   S_ = Sij  ! save
 
   ! check degeneracy & indexing eigenvalue and subspace
   do ie=1,neig*ispinor - 1
     if(abs(E(ie)-E(ie+1)) .gt. very_small .and. ie .ne. neig*ispinor-1) then
       iend = ie
       call cal_eig_nonsymm(S_(init:iend,init:iend), iend-init+1, S_eig(init:iend))
       init = ie+1
     elseif(abs(E(ie)-E(ie+1)) .gt. very_small .and. ie .eq. neig*ispinor-1) then
       iend = ie
       call cal_eig_nonsymm(S_(init:iend,init:iend), iend-init+1, S_eig(init:iend))
       init = ie+1
       iend = ie+1
       call cal_eig_nonsymm(S_(init:iend,init:iend), iend-init+1, S_eig(init:iend))
     elseif(abs(E(ie)-E(ie+1)) .lt. very_small .and. ie .eq. neig*ispinor-1) then
       iend = ie + 1
       call cal_eig_nonsymm(S_(init:iend,init:iend), iend-init+1, S_eig(init:iend))

     endif
   enddo

   return
endsubroutine
