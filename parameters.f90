module parameters
  use mpi_setup
  use iso_c_binding, only:  c_char, c_int, c_double, c_ptr, c_f_pointer
#ifdef MKL_SPARSE
  use MKL_SPBLAS 
#endif
!#ifndef MPI
!  integer*4,    public, parameter :: nprocs = 1
!  integer*4,    public, parameter :: myid   = 0
!#endif
  character*26, public, parameter ::alphabet='abcdefghijklmnopqrstuvwxyz'
  real*8    , public, parameter   ::      pi=4.d0*atan(1.d0) ! Leibniz's formula for Pi
  real*8    , public, parameter   ::     pi2=2.d0*pi
  real*8    , public, parameter   ::    bohr=0.52917721067d0 ! meter Bohr radius
  real*8    , public, parameter   :: hartree=27.21138602d0 ! eV   Hartree energy
  real*8    , public, parameter   ::    hbar=4.135667662d-15/pi2  ! eV*s Plank constant
  real*8    , public, parameter   ::       c=0.262465831d0 !! constant c = 2m/hbar**2 [1/eV Ang^2] 
  real*8    , public, parameter   :: g_elect=2.0023193043617 !! g-factor[see Mol.Phys. 98, 1597(2000) for sign]
  real*8    , public, parameter   ::     rt2=sin( 4.d0*atan(1.d0)/4.d0 ) * 2.d0 ! sqrt(2)
  real*8    , public, parameter   ::     rt3=sin( 4.d0*atan(1.d0)/3.d0 ) * 2.d0 ! sqrt(3)
  real*8    , public, parameter   :: onsite_tolerance= 0.0001 !! symmetry precision
  real*8    , public, parameter   ::     eta=epsilon(1d0) ! tiny value
  real*8    , public              :: t1, t0 ! time
  complex*16, public, parameter   ::      zi=(0.d0,1.d0)
  complex*16, public, parameter   ::     pzi= pi*zi
  complex*16, public, parameter   ::    pzi2=2*pzi
  complex*16, public, dimension(2,2),   parameter :: pauli_0 = reshape((/  1,  0,  0,  1 /), (/2,2/))
  complex*16, public, dimension(2,2),   parameter :: pauli_x = reshape((/  0,  1,  1,  0 /), (/2,2/))
  complex*16, public, dimension(2,2),   parameter :: pauli_y = reshape((/  0,  1, -1,  0 /) * zi, (/2,2/))
  complex*16, public, dimension(2,2),   parameter :: pauli_z = reshape((/  1,  0,  0, -1 /), (/2,2/))
  integer*4,  public, dimension(3,3,3), parameter :: levi_civita = reshape((/0,0,0, 0,0,-1, 0,1,0, & 
                                                                     0,0,1, 0,0,0, -1,0,0, & 
                                                                     0,-1,0, 1,0,0, 0,0,0/), (/3,3,3/))
  integer*4,  public, dimension(3,2),   parameter :: cyclic_axis = reshape((/2,3,1,3,1,2/), (/3,2/))

  integer*4,  public, parameter   :: max_nparam = 500    !! maximum number of onsite and hopping parameters
  integer*4,  public, parameter   :: max_kpoint = 100000 !! maximum number of kpoints
  integer*4,  public, parameter   :: max_set_weight = 100000  !! maximum number of "SET WEIGHT" tag
  integer*4,  public, parameter   :: max_pair_type  = 1000    !! maximum number of NN_CLASS pair type
  integer*4,  public, parameter   :: max_dummy      = 9999999 !! maximun number of dummy index for arbitral purpose   
  integer*4,  public, parameter   :: max_nsym       = 1000000 !! maximun number of symmetry operation for SPGLIB
  integer*4,  public, parameter   :: pid_energy      = 30 
  integer*4,  public, parameter   :: pid_nntable     = 33
  integer*4,  public, parameter   :: pid_incar       = 78
  integer*4,  public, parameter   :: pid_kpoint      = 79 
  integer*4,  public, parameter   :: pid_stm         = 81 
  integer*4,  public, parameter   :: pid_chg         = 82 
  integer*4,  public, parameter   :: pid_ibzkpt      = 83 
  integer*4,  public, parameter   :: pid_dos         = 84 
  integer*4,  public, parameter   :: pid_param       = 88 
  integer*4,  public, parameter   :: pid_matrix      = 99
  integer*4,  public, parameter   :: pid_berrycurv   = 100
  integer*4,  public, parameter   :: pid_zak         = 101
  integer*4,  public, parameter   :: pid_wcc         = 102 ! pid_wcc + 1 = pid_gap
  integer*4,  public, parameter   :: pid_geom        = 177 
  integer*4,  public, parameter   :: pid_geom_ribbon = 178

  type incar !PINPT
       logical                       flag_get_band  ! default = .true.

       logical                       flag_tbfit_parse, flag_tbfit_parse_
       logical                       flag_tbfit_test

       real*8                        ptol
       real*8                        ftol
       integer*4                     miter,nparam,nparam_const
       integer*4                     read_energy_column_index, read_energy_column_index_dn
       logical                       flag_tbfit, flag_pfile, flag_pincar
       logical                       flag_print_only_target, flag_print_param
       logical                       flag_print_orbital, flag_get_orbital
       logical                       flag_set_param_const
       logical                       flag_slater_koster ! default .true.
       logical                       flag_print_mag
       logical                       flag_load_nntable ! default .false.
       character*2                   axis_print_mag
       real*8,       allocatable  :: param(:)
       character*20, allocatable  :: param_name(:)
       character*40, allocatable  :: c_const(:,:)
       real*8,       allocatable  :: param_const(:,:) ! i=1 -> 'is same as'
                                                      ! i=2 -> 'is lower than' (.le.) : maximum bound  ! not available yet
                                                      ! i=3 -> 'is lower than' (.ge.) : minimum bound  ! not available yet
                                                      ! i=4 -> 'is fixed' : fixed during fitting       ! original value will be 
                                                                                                       ! copied to PINPT%param_const(i=5,:)
       character*8                   ls_type   ! fitting method
       character*132                 kfilenm,gfilenm    ! kpoint file, geometry file
       character*132                 ribbon_kfilenm     ! kpoint file for ribbon geometry defined in 'SET RIBBON'
       character*132                 pfilenm,pfileoutnm ! SK parameter file input & output
       character*132                 efilenmu,efilenmd  ! target energy file (spin up & dn)
       character*132                 nnfilenm           ! hopping integral file name (default = hopping.dat)     

       integer*8                     nweight, npenalty_orb
       character*132, allocatable :: strip_kp(:), strip_tb(:), strip_df(:), strip_wt(:)
       character*132, allocatable :: strip_kp_orb(:), strip_tb_orb(:), strip_orb(:), strip_site(:), strip_pen_orb(:)

       ! plot_eig mode
       logical                       flag_plot_eigen_state  !default : .false.
       logical                       flag_plot_wavefunction !default : .true. => WAV_PLOT = .true. if .not. CHG_PLOT = .true.
       logical                       flag_default_ngrid, flag_default_rorigin
       integer*4                     n_eig_print,n_kpt_print
       integer*4,    allocatable  :: i_eig_print(:)
       integer*4,    allocatable  :: i_kpt_print(:)
       integer*4                     ngrid(3)
       real*8                        r_origin(3)   ! direct coordinate for shift of the origin of chg file
       real*8                        rcut_orb_plot ! cutoff radius for orbital plot (applies STM mode too if provied in SET STMPLOT)
                                                   ! default = 5 Ang (see PRB 64, 235111)
       ! plot_stm mode
       logical                       flag_repeat_cell_orb_plot(3) ! default = .true. true. true.
       integer*4                     repeat_cell_orb_plot(3) ! default = 1
       logical                       flag_plot_stm_image    !default : .false.
       logical                       flag_default_stm_ngrid, flag_default_stm_rorigin
       integer*4                     n_stm        ! how many stm images will be plotted
       integer*4                     stm_ngrid(3) 
       real*8,       allocatable  :: stm_emax(:), stm_emin(:)

       logical                       flag_get_dos

       logical                       flag_get_wcc
       logical                       flag_get_zak_phase, flag_zak_separate, flag_zak_kfile_read
       logical                       flag_get_z2
       logical                       flag_get_berry_curvature, flag_berryc_separate
       logical                       flag_berry
       logical                       flag_get_parity

       integer*4                     ispin   ! nonmag: 1, collinear: 2, non-collinear: 2
       integer*4                     ispinor ! nonmag: 1, collinear: 1, non-collinear: 2
       integer*4                     nspin   ! nonmag: 1, collinear: 2, non-collinear: 1
       logical                       flag_local_charge, flag_collinear, flag_noncollinear, flag_soc

       logical                       flag_erange
       integer*4                     init_erange, fina_erange  ! ie:fe
       integer*4                     nband ! nonmag: ie:fe, collinear: ie:fe for up or dn, non-collinear: ie:fe
       logical                       flag_sparse ! if EWINDOW tag has been set up, flag_sparse is forced to be .true.
       integer*4                     feast_nemax ! maximum number of eigenvalues (=nband_guess)
       real*8                        feast_emin, feast_emax ! energy window [emin:emax] for FEAST algorithm
       integer*4                     feast_fpm(128) ! FEAST parameters
       integer*4,    allocatable  :: feast_ne(:,:)  ! Number of states found in erange [emin:emax], (:,:) = (nspin,nkpoint)
       integer*4                     feast_neguess  ! initial guess for the number of states to be found in interval

       logical                       flag_set_ribbon, flag_print_only_ribbon_geom
       integer*4                     ribbon_nslab(3)
       real*8                        ribbon_vacuum(3)

       logical                       flag_plus_U

       logical                       flag_scissor
       integer*4                     i_scissor ! i_scissor level will be operated by scissor operator 
       real*8                        r_scissor ! EDFT(n,k) + r_scissor if n >= i_scissor

       logical                       flag_efield, flag_efield_frac, flag_efield_cart
       real*8                        efield(3)
       real*8                        efield_origin(3)
       real*8                        efield_origin_cart(3)

  endtype incar

  type poscar !PGEOM
       integer*4                     neig       ! either neig_up and neig_dn
       integer*4                     neig_total ! neig_up + neig_dn
       integer*4                     neig_target
       integer*4                     nbasis  ! normally neig = nbasis

       integer*4                     n_spec, n_atom
       integer*4,   allocatable   :: i_spec(:) ! number of atoms per each species (1:n_spec)
       character*8, allocatable   :: c_spec(:) ! character of species for each species (1:n_spec)
       integer*4,   allocatable   :: spec(:)   ! species information for each atom (1:n_atom)
       real*8,      allocatable   :: a_coord(:,:) ! atomic  coordinate (1:3, 1:n_atom) (direct, fractional)
       real*8,      allocatable   :: o_coord(:,:) ! orbital coordinate (1:3, 1:neig) (direct, fractional)
       real*8,      allocatable   :: o_coord_cart(:,:) ! orbital coordinate (1:3, 1:neig) (cartesian)

       integer*4                     max_orb  ! maximum number of orbitals asigned in each atomic site
       integer*4,   allocatable   :: n_orbital(:) ! number of orbitals per atomic site
       character*8, allocatable   :: c_orbital(:,:) ! name of atomic orbitals for each atomic sites

       character*40                  system_name
       real*8                        a_scale, a_latt(3,3) ! lattice vector (unit of Ang)
       real*8                        b_latt(3,3) ! reciprocal lattice vector(unit of A^-1)
       logical                       flag_selective, flag_direct, flag_cartesian

       integer*4                     n_nn_type
       character*80,allocatable   :: nn_pair(:)
       real*8,      allocatable   :: nn_dist(:), nn_r0(:)

       real*8,      allocatable   :: local_charge(:,:)   ! local net charge (rho_up - rho_dn)
                                          
       !SPGLIB related variables
       integer*4                     spg_error
       integer*4                     spg_space_group !space group index
       integer*4                     spg_hall_number
       integer*4                     spg_Hermann_Mauguin_number
       real*8                        spg_transformation_matrix(3,3)
       real*8                        spg_origin_shift(3)
       integer*4                     spg_n_operations ! nsym, number of symmetry operations
       character*12                  spg_international
       character*18                  spg_hall_symbol
       character*6                   spg_choice
       character*7                   spg_point_group
       character*12                  spg_crystal_system
       integer*4,   allocatable   :: spg_rotations(:,:,:)   ! {->w,   t} (3,3,spg_n_operations)
       real*8,      allocatable   :: spg_translations(:,:) !  {  w, ->t} (3,spg_n_operations)
       integer*4,   allocatable   :: spg_wyckoffs(:)
       integer*4,   allocatable   :: spg_equivalent_atoms(:) 
       real*8,      allocatable   :: spg_a_coord_operated(:,:,:) ! (3,n_atom,spg_n_operations)
       character*7                   spg_schoenflies
       real*8                        spg_a_latt_primitive(3,3)
       real*8,      allocatable   :: spg_a_coord_primitive(:,:) !(3,spg_n_atom_primitive)
       integer*4,   allocatable   :: spg_spec_primitive(:)      !(  spg_n_atom_primitive)
       integer*4                     spg_n_atom_primitive
       integer*4,   allocatable   :: spg_det_w(:) ! (spg_n_operation) ! determinant of rotation w (1=proper, -1=improper)
       integer*4,   allocatable   :: spg_tr_w(:) ! (spg_n_operation) ! trace of rotation w 
       integer*4,   allocatable   :: spg_type_w(:) ! types of rotation operation of space group
  endtype poscar

  type kpoints !PKPTS
       integer*4                     nkpoint,nline
       integer*4,   allocatable   :: ndiv(:)
       real*8,      allocatable   :: kpoint(:,:),kline(:,:)
       real*8,      allocatable   :: kpoint_reci(:,:)
       real*8                        k_shift(3)
       character*8, allocatable   :: k_name(:)
       logical                       flag_klinemode
       logical                       flag_kgridmode, flag_gamma
       logical                       flag_reciprocal, flag_cartesianK
       character*1                   kunit
  endtype kpoints

  type energy !EDFT / ETBA / ETBA_DOS
       real*8,      allocatable   :: E(:,:)   !E(neig*ispin,nkpoint) (for nspin=2, up=1:neig, dn=neig+1:neig*2)
       complex*16,  allocatable   :: V(:,:,:) !V(nbasis=neig*ispin,neig*ispin,nkpoint) order is same as E above
  endtype energy

  type weight !PWGHT
       integer*4                     nweight, iband, fband
       integer*4                     npenalty_orb
       real*8,      allocatable   :: WT(:,:)
       real*8,      allocatable   :: PENALTY_ORB(:,:,:)
       logical                       flag_weight_default
       logical                       flag_weight_default_orb
  endtype weight

  type hopping !NN_TABLE (nearest neighbor table; but not restricted to nearest)
       integer*4,   allocatable   :: i_atom(:)     !(n) n = nn_neighbor index
       integer*4,   allocatable   :: j_atom(:)
       real*8,      allocatable   :: i_coord(:,:)  ! cartesian coordinate of atom_i (1:3,n) = nn_neighbor index
       real*8,      allocatable   :: j_coord(:,:)  ! cartesian coordinate of atom_j
       real*8,      allocatable   :: Rij(:,:) ! vector connecting two orbital i and j
       real*8,      allocatable   :: R(:,:)   ! cell periodicity where orbital j sits on
       real*8,      allocatable   :: Dij(:)   ! distance between orbital i and j
       real*8,      allocatable   :: Dij0(:)
       integer*4,   allocatable   :: i_matrix(:) 
       character*8, allocatable   :: ci_orb(:)
       integer*4,   allocatable   :: j_matrix(:)
       character*8, allocatable   :: cj_orb(:)
       character*2, allocatable   :: p_class(:)
       integer*4,   allocatable   :: n_class(:)
       integer*4,   allocatable   :: sk_index_set(:,:) ! (i,nn), i=onsite(0),sigma(1),  pi(2),  delta(3)
                                                       !           scaled->s_sigma(4),s_pi(5),s_delta(6)
       integer*4,   allocatable   :: cc_index_set(:,:) ! (i,nn), i=onsite(0),t_typ(1),t_typ(2),t_typ(3), t_typ(...)
       real*8,      allocatable   :: tij(:) ! hopping amplitude (except SOC, magnetic coupling)
       real*8                        onsite_tolerance
       integer*4                     n_neighbor

       character*20,allocatable   :: site_cindex(:)  ! site indicator
       real*8,      allocatable   :: local_charge(:) !i=ham_index(neig)
       real*8,      allocatable   :: local_moment(:,:) !(1:3,i) (1:3)=(mi, theta, phi), i=ham_idx(neig)
       real*8,      allocatable   :: local_moment_rot(:,:) !(1:3,i) (1:3)=(mx, my, mz)=-I*mi_dot_sigma, i=ham_idx
       integer*4,   allocatable   :: stoner_I_param_index(:) ! array size = i, i=ham_index(neig), 
       integer*4,   allocatable   :: local_U_param_index(:) ! array size = i, i=ham_index(neig), 
       integer*4,   allocatable   :: plus_U_param_index(:) ! array size = i, i=ham_index(neig), 
       integer*4,   allocatable   :: soc_param_index(:)    ! soc parameter index for each orbital-orbital pair (nn)
!      real*8,      allocatable   :: orbital_moment(:,:)   ! (1:3,i) (1:3)=(Lx, Ly, Lz), i=nn_index (nn-class=0)

       real*8,      allocatable   :: tij_file(:) ! hopping amplitude read from file


  endtype hopping

  type dos ! PINPT_DOS
       character*40                  dos_kfilenm,dos_filenm
       integer*4                     dos_kgrid(3)
       integer*4                     dos_nediv 
       integer*4                     dos_iband, dos_fband
       real*8                        dos_emin, dos_emax
       real*8                        dos_smearing
       real*8                        dos_kshift(1:3)
       real*8,      allocatable   :: dos_kpoint(:,:)       
       real*8,      allocatable   :: dos(:), dos_erange(:)   ! size=nediv
       real*8,      allocatable   :: dos_up(:), dos_dn(:)    ! size=nediv
       logical                       dos_flag_gamma, dos_flag_print_kpoint
   
       logical                       dos_flag_print_eigen
       integer*4,   allocatable   :: dos_ensurf(:) !integer array of band index. size=n_ensurf
       integer*4                     dos_n_ensurf
       character*1                   dos_kunit
  endtype dos

  type berry ! PINPT_BERRY
       logical                       flag_wcc_evolve
       character*40                  wcc_filenm, wcc_gap_filenm
       integer*4                     wcc_nerange
       integer*4,   allocatable   :: wcc_erange(:)
       real*8,      allocatable   :: wcc(:,:,:) ! (nerange/nspin,nspin,nkpath)
       real*8,      allocatable   :: wcc_kpath(:,:,:) ! k-path (1:3, 2, nkline) (1:3, init-fina, nkline)
       real*8                        wcc_kpath_shift(3)
       real*8,      allocatable   :: wcc_kpoint(:,:,:) ! k-points along the k-path (1:3,1:nk,nkline)
       real*8,      allocatable   :: wcc_kpoint_reci(:,:,:) ! k-points along the k-path (1:3,1:nk,nkline) (reci unit)
       integer*4                     wcc_nkdiv, wcc_nkdiv2         ! nkdiv along k-path, nkdiv2 along k-path wcc_evolve
       integer*4                     wcc_nkpath
       integer*4                     wcc_direction
       character*256                 strip_wcc_range
       logical                       flag_wcc_phase
       logical                       flag_wcc_get_chern
       logical                       flag_wcc_get_chern_spin
       real*8,      allocatable   :: wcc_polarization(:,:)   ! (nspin, nkpath)
       real*8,      allocatable   :: wcc_chern(:)           ! chern number (nspin)

       logical                       flag_zak_evolve
       character*40                  zak_filenm 
       integer*4                     zak_nerange
       integer*4,   allocatable   :: zak_erange(:)          ! default size = neig * ispin
       real*8,      allocatable   :: zak_phase(:,:) !(nspin,nkpath) = (erange)
       real*8,      allocatable   :: zak_kpath(:,:,:) ! k-path (1:3, 2, nkline) (1:3, init-fina, nkpath)
       real*8                        zak_kpath_shift(3)
       real*8,      allocatable   :: zak_kpoint(:,:,:) ! k-points along the k-path (1:3,1:nk,nkpath)
       real*8,      allocatable   :: zak_kpoint_reci(:,:,:) ! k-points along the k-path (1:3,1:nk,nkpath) (reci unit)
       integer*4                     zak_nkdiv, zak_nkdiv2  ! nkdiv along one k-path, nkdiv2 how many k-path (zak_evolve) along zak_dir
       integer*4                     zak_nkpath
       integer*4                     zak_direction
       character*256                 strip_zak_range
       real*8,      allocatable   :: polarization(:)
       logical                       flag_zak_phase
 

!      logical                       
       character*256                 strip_z2_range
       character*40                  z2_filenm, z2_gap_filenm
       integer*4                     z2_nkdiv, z2_nkpath
       integer*4                     z2_nerange
       integer*4,   allocatable   :: z2_erange(:)          ! default size = neig * ispin
       integer*4                     z2_dimension, z2_nplane
       integer*4,   allocatable   :: z2_axis(:)
       real*8,      allocatable   :: z2_kpoint(:,:,:,:,:)      ! kp along k-path (3,nkdiv,nkpath,nplane,naxis) (nxis =3(3D), 1(2D,1D))
       real*8,      allocatable   :: z2_kpoint_reci(:,:,:,:,:) ! kp along k-path (3,nkdiv,nkpath,nplane,naxis) (nplane = 2(3D), 1(2D,1D))
       real*8,      allocatable   :: z2_wcc(:,:,:,:,:)         ! wcc (nerange/nspin,nspin,nkpath,nplane,naxis)
       logical                       flag_z2_get_chern
       real*8,      allocatable   :: z2_chern(:,:,:)           ! chern number (nspin, nplane, naxis)
       real*8,      allocatable   :: z2_polarization(:,:,:,:)           ! chern number (nspin, nkpath, nplane, naxis)
       logical                       flag_z2_phase

       character*40                  bc_filenm 
       logical                       flag_bc_filenm_provided
       logical                       flag_bc_method_kubo, flag_bc_method_resta
       integer*4                     bc_dimension
       integer*4                     bc_nkdiv(3)
       integer*4,   allocatable   :: bc_axis(:)
       integer*4                     bc_nerange
       integer*4,   allocatable   :: bc_erange(:) ! default size = neig * ispin
       real*8,      allocatable   :: berrycurv(:,:,:) ! (:,:,:) = (1:3 ->Omega_x:z, nkpoints, erange)
       real*8,      allocatable   :: omega(:,:,:,:) ! (neig * ispinor, 3, nspin, nkpoint)
       character*256                 strip_bc_range
       logical                       flag_bc_phase

       real*8                        parity_origin(3) !direct coordinate of the orgiin of the system in the given unit cell
       real*8                        parity_operator(3,3) 
       integer*4                     parity_nkpoint
       real*8,      allocatable   :: parity_kpoint(:,:) !(3,parity_nkpoint)
       real*8,      allocatable   :: parity_kpoint_reci(:,:) !(3,parity_nkpoint)
       character*10,allocatable   :: parity_kpoint_name(:)

  endtype berry

  type spmat ! sparse matrix with Compressed Sparse Row format
       integer*4                     nnz    ! number of non-zero elements
       integer*4                     msize  ! matrix size
       complex*16  ,allocatable   :: H(:)   ! sparse array of square matrix H_square(msize,msize), (n_neighbor)
       integer*4   ,allocatable   :: I(:)   ! Row    array : if CSR format : I(msize + 1), I(msize+1)-1 = nnz
                                            !                if COO format : I(nnz)
       integer*4   ,allocatable   :: J(:)   ! Column array : J(nnz), same size with H
!#ifdef MKL_SPARSE
!       type(c_ptr)                :: cptrH   , cptrB   , cptrE   , cptrJ
!       integer*4,       pointer   ::            ptrB(:),  ptrE(:),  ptrJ(:)
!       complex*16,      pointer   ::  ptrH(:)
!#endif
  endtype spmat
endmodule 
