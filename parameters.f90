module parameters
  use mpi_setup
  use iso_c_binding, only:  c_char, c_int, c_double, c_ptr, c_f_pointer
#ifdef MKL_SPARSE
  use MKL_SPBLAS 
#endif
  character*26, public, parameter ::alphabet='abcdefghijklmnopqrstuvwxyz'
  real*8    , public, parameter   ::      pi=4.d0*atan(1.d0) ! Leibniz's formula for Pi
  real*8    , public, parameter   ::     pi2=2.d0*pi
  real*8    , public, parameter   ::    bohr=0.52917721067d0 ! meter Bohr radius
  real*8    , public, parameter   :: hartree=27.21138602d0 ! eV   Hartree energy
  real*8    , public, parameter   ::    hbar=4.135667662d-15/pi2  ! eV*s Plank constant
  real*8    , public, parameter   ::       c=0.262465831d0 !! constant c = 2m/hbar**2 [1/eV Ang^2] 
  real*8    , public, parameter   ::     k_B=8.6173303d-5  !! Boltzmann constant = 2m/hbar**2 [1/eV Ang^2] 
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

  integer*4,  public, parameter   :: max_nparam      = 500    !! maximum number of onsite and hopping parameters in total
  integer*4,  public, parameter   :: max_kpoint      = 100000 !! maximum number of kpoints
  integer*4,  public, parameter   :: max_set_weight  = 100000  !! maximum number of "SET WEIGHT" tag
  integer*4,  public, parameter   :: max_pair_type   = 1000    !! maximum number of NN_CLASS pair type
  integer*4,  public, parameter   :: max_dummy       = 9999999 !! maximun number of dummy index for arbitral purpose   
  integer*4,  public, parameter   :: max_dummy2      = 1000    !! maximun number of dummy index for arbitral purpose   
  integer*4,  public, parameter   :: max_nsym        = 1000000 !! maximun number of symmetry operation for SPGLIB
  integer*4,  public, parameter   :: pid_energy      = 30 
  integer*4,  public, parameter   :: pid_nntable     = 33
  integer*4,  public, parameter   :: pid_incar       = 78
  integer*4,  public, parameter   :: pid_kpoint      = 79 
  integer*4,  public, parameter   :: pid_stm         = 81 
  integer*4,  public, parameter   :: pid_chg         = 82 
  integer*4,  public, parameter   :: pid_ibzkpt      = 83 
  integer*4,  public, parameter   :: pid_dos         = 84 
  integer*4,  public, parameter   :: pid_ldos        = 85 
  integer*4,  public, parameter   :: pid_param       = 88 
  integer*4,  public, parameter   :: pid_matrix      = 99
  integer*4,  public, parameter   :: pid_berrycurv   = 100
  integer*4,  public, parameter   :: pid_zak         = 101
  integer*4,  public, parameter   :: pid_wcc         = 102 ! pid_wcc + 1 = pid_gap
  integer*4,  public, parameter   :: pid_geom        = 177 
  integer*4,  public, parameter   :: pid_geom_ribbon = 178
  integer*4,  public, parameter   :: pid_circ        = 179
 
  type incar !PINPT
       ! set in parsing step
       integer*4                     nsystem        ! number of geometry
       character*132, allocatable :: ifilenm(:)     ! input tag file (default = INCAR-TB, specified by -input ifilenm in command-line)
       character*132, allocatable :: title(:)       ! title of the system
       character*132                 pfilenm_parse  ! parsed parameter file name
       character*132                 kfilenm_parse  ! parsed kpoints   file name

       character*40                  fnamelog       ! log file name, default = TBFIT.log
       logical                       flag_get_band  ! default = .true.
       logical                       flag_spglib    ! default = .true. ! write space group information
       logical                       flag_tbfit_parse, flag_tbfit_parse_
       logical                       flag_kfile_parse
       logical                       flag_pfile_parse
       logical                       flag_ndiv_line_parse, flag_ndiv_grid_parse
       logical                       flag_fit_degeneracy ! fitting is helped by fitting degeneracy as well.
       logical                       flag_miter_parse, flag_mxfit_parse
       logical                       flag_lorbit_parse, flag_proj_parse
       logical                       flag_parse
       logical                       flag_tbfit_test
       logical                       flag_inputcard_fname_parse
       logical                       flag_ga_with_lmdif ! default = PKAIA%flag_ga_with_lmdif
       logical                       flag_write_unformatted ! default = .false.
       logical                       flag_report_geom
   
       logical                       flag_phase  ! default = .true.  ! apply phase factor to atomic orbital

       real*8                        ptol
       real*8                        ftol
       real*8                        fdiff
       integer*4                     ga_npop
       integer*4                     miter,mxfit
       logical                       flag_tbfit, flag_pincar
       logical                       flag_tbfit_finish ! when exiting fitting routine, it will be turn to .true.
       logical                       flag_print_only_target
       logical                       flag_print_energy_diff ! also print energy different between target and tight-binding energies
       logical                       flag_print_orbital  ! activated with LORBIT .TRUE. tag
       logical                       flag_get_orbital    ! whether request wavevector in diagonalize routine
       logical                       flag_print_proj
       logical                       flag_print_mag
       logical                       flag_print_single ! default .false. (write single precision for wavefunction)
       logical                       flag_load_nntable ! default .false.
       character*2                   axis_print_mag ! mx, my, mz (pauli matrices), 
                                                    ! re, im (real or imag part), 
                                                    ! wf (full wf), bi (enforce to write wf with binary format)

       character*8                   ls_type   ! fitting method
       integer*4                     nn_max(3)          ! cell reapeat for the nearest neighbor finding in find_nn routine (default:3 3 3)
       logical                       flag_use_weight    ! if true use "weight" information written in PFILE and 
                                                        ! replace it with SET WEIGHT info of INCAR-TB file after fitting procedures

       ! plot_eig mode
       logical                       flag_plot_eigen_state  !default : .false.
       logical                       flag_plot_wavefunction !default : .true. => WAV_PLOT = .true. if .not. CHG_PLOT = .true.
       integer*4                     n_eig_print,n_kpt_print
       integer*4,    allocatable  :: i_eig_print(:)
       integer*4,    allocatable  :: i_kpt_print(:)
       real*8                        rcut_orb_plot ! cutoff radius for orbital plot (applies STM mode too if provied in SET STMPLOT)
                                                   ! default = 5 Ang (see PRB 64, 235111)
       ! plot_stm mode
       logical                       flag_repeat_cell_orb_plot(3) ! default = .true. true. true.
       integer*4                     repeat_cell_orb_plot(3) ! default = 1
       logical                       flag_plot_stm_image    !default : .false.
       integer*4                     n_stm        ! how many stm images will be plotted
       real*8,       allocatable  :: stm_emax(:), stm_emin(:)

       logical                       flag_get_dos

       logical                       flag_get_wcc
       logical                       flag_get_zak_phase, flag_zak_separate, flag_zak_kfile_read
       logical                       flag_get_z2
       logical                       flag_get_berry_curvature, flag_berryc_separate
       logical                       flag_berry
       logical                       flag_get_parity, flag_get_symmetry

       !                             nonmag: ispin = 1, nspin = 1, ispinor = 1
       !                             noncol: ispin = 2, nspin = 1, ispinor = 2
       !                             collin: ispin = 2, nspin = 2, ispinor = 1
       integer*4                     ispin   ! nonmag: 1, collinear: 2, non-collinear: 2
       integer*4                     ispinor ! nonmag: 1, collinear: 1, non-collinear: 2
       integer*4                     nspin   ! nonmag: 1, collinear: 2, non-collinear: 1
       logical                       flag_local_charge, flag_collinear, flag_noncollinear, flag_soc

       logical                       flag_erange
       logical                       flag_sparse ! if EWINDOW tag has been set up, flag_sparse is forced to be .true.
       integer*4                     feast_nemax ! maximum number of eigenvalues (=nband_guess)
       real*8                        feast_emin, feast_emax ! energy window [emin:emax] for FEAST algorithm
#ifdef PSPARSE
       integer*4                     feast_fpm(64)  ! FEAST parameters set by original feastinit
#else
       integer*4                     feast_fpm(128) ! FEAST parameters set by MKL extended eigensolver feast support routine
#endif
       integer*4,    allocatable  :: feast_ne(:,:)  ! Number of states found in erange [emin:emax], (:,:) = (nspin,nkpoint)
       integer*4                     feast_neguess  ! initial guess for the number of states to be found in interval

       logical                       flag_plus_U

       logical                       flag_scissor
       integer*4                     i_scissor ! i_scissor level will be operated by scissor operator 
       real*8                        r_scissor ! EDFT(n,k) + r_scissor if n >= i_scissor

       ! construct effective hamiltonian
       logical                       flag_get_effective_ham
       character*132                 eff_orb_dummyc
       real*8                        eff_emin, eff_emax

       ! projected band
       integer*4                     nproj_sum
       integer*4,   allocatable   :: proj_atom(:,:) ! integer array of atom index. maxsize=n_atom
       integer*4,   allocatable   :: proj_natom(:) ! how many atoms to be plotted for projected band
       logical                       flag_print_proj_sum

       logical                       flag_plot_fit
       logical                       flag_plot
       character*132                 filenm_gnuplot, filenm_gnuplot_parse
       logical                       flag_filenm_gnuplot_parse

       logical                       flag_print_energy_singlek
       logical                       flag_print_hamk

       logical                       flag_get_band_order   ! flag whether perform band re-ordering by calculating overlap integral
       logical                       flag_get_band_order_print_only ! if true, band re-order will not be performed in the fitting routines
       real*8                        band_order_overlap_cutoff ! cutoff of overlap integral to perform eigenvalue swap : sqrt(2)/2 by default

       logical                       flag_get_total_energy 
       real*8                        electronic_temperature ! Temperature, default = 0 (K)

       ! circular dichroism 
       logical                       flag_get_circ_dichroism ! flag whether perform circular dichroism calculation
       integer*4                     ncirc_dichroism         
       integer*4,   allocatable   :: circ_dichroism_pair(:,:)

       integer*4                     lmmax ! 9 (default) =>  1:s, 2~4:px,py,pz, 5~9:dz2,dx2,dxy,dxz,dyz
                                           ! 3           =>  1:s, 2  :px,py,pz, 3  :dz2,dx2,dxy,dxz,dyz
       logical                       flag_fit_orbital ! fit orbital character as well?
       real*8                        orbital_fit_smearing ! gaussian smearing for orbital fitting                                           

  endtype incar

  type params !PPRAM
     ! integer*4                     mysystem   ! my system index (deprecated...)
       logical                       flag_set_param_const     ! whether set constrain on parameters
       logical                       flag_pfile_index         ! whether put numbering in PARAM_FIT.dat after fitting
       logical                       flag_use_overlap         ! whether setup overlap hamiltonian. Automatically activated if the overlap integral parameters
                                                              ! are provided in the PARAM_FIT.dat file in priori (start with o_ )
       real*8                        l_broaden                ! broadening of the cutoff-function for NRL type SK parameters (used if SK_SCALE_MODE > 10 in PFILE)
       logical                       flag_slater_koster       ! default .true.
       logical                       flag_nrl_slater_koster   ! default .false.
       integer*4                     slater_koster_type       ! 1 ~ 3: exponential scaling, 11: Mehl & Papaconstantopoulos NRL method (PRB 54, 4519 (1996))
       character*132                 pfilenm,pfileoutnm       ! SK parameter file input & output
       integer*4                     nparam                   ! total number of parameters 
       integer*4                     nparam_const             ! total number of parameters constraints
       integer*4                     nparam_free              ! total number of free parameters size(iparam_free(:))
       real*8,       allocatable  :: param(:)                 ! TB parameters
       real*8,       allocatable  :: param_nrl(:,:)           ! NRL TB parameters
       integer*4,    allocatable  :: iparam_free(:)           ! parameter index for free parameters only (nparam_free) 
       integer*4,    allocatable  :: iparam_free_nrl(:)       ! parameter index for free parameters only (nparam_free), sum(param_nsub(iparam_free(1:j-1)))+1
       character*40, allocatable  :: param_name(:)            
       character*40, allocatable  :: c_const(:,:)             
       real*8,       allocatable  :: param_const(:,:)         ! i=1 -> 'is same as'
                                                              ! i=2 -> 'is lower than' (.le.) : maximum bound  ! <=  20 
                                                              ! i=3 -> 'is lower than' (.ge.) : minimum bound  ! >= -20 or >= 0.001 (if scale factor)
                                                              ! i=4 -> 'is fixed' : fixed during fitting       ! original value will be 
                                                                                                               ! copied to PINPT%param_const(i=5,:)
       real*8,       allocatable  :: param_const_nrl(:,:,:)   ! (i,j,nparam) ! j=1:4 -> a, b, c, d  for sk-parameters 
                                                              !                      -> alpha, beta, gamma, xi for onsite parameters
                                                              !                if size(,:,) = 1, then lonsite_ or lambda_, stoner_ , etc....
                                                              ! note: this parameter is only activated if slater_koster_type = 11
       integer*4                     nparam_nrl               ! total number of parameters sum(param_nsub(:))
       integer*4                     nparam_nrl_free          ! total number of free parameters sum(param_nsub(iparam_free(:)))
       integer*4                     param_nsub_max           ! 4 if SK_SCALE_MODE >10, 1 if SK_SCALE_MODE <= 10
       integer*4,    allocatable  :: param_nsub(:)            ! (nparam) number of sub-parameters for each parameter
                                                              ! 4: NRL type of e_, ssp_, pds_,...(with s_, o_, os_ )-> onsite, hopping, hopping_scale, 
                                                              !                                               hopping(overlap),hopping_scale(overlap)
                                                              ! 1: other parameters -> local_U, lrashba_, lambda_, lsoc_, stoner_I_ , ..., etc.
  endtype params

  type poscar !PGEOM
       integer*4                     mysystem   ! my system index

       character*132                 title      ! system name (should have no blank)
       character*132                 gfilenm    ! geometry file
       real*8,      allocatable   :: nelect(:)  ! number of total electrons in the system  (nspin, for up & dn if nspin=2), set by NELECT tag
                                                ! this parameter is used in total energy calculations to determine Fermi level of the system.
                                                ! NOTE(06.Oct.2020): The default will be automatically set in the future as following:
                                                ! by default this will be determined based on the valence electron configuration of the atom.
                                                ! For example, if you have two Si(3s1_3p3) atom and four H(1s1) atom, 
                                                ! then you will have 4*2 + 1*4 = 12 electrons in total with non-magnetic or noncollinear system 
                                                ! and 6 (up) and 6 (dn) electrons, respectively, with collinear magnetic systems.
       integer*4                     neig       ! either neig_up and neig_dn, total number of atomic orbitals
                                                ! = sum(PGEOM%n_orbital(1:PGEOM%n_atom)) (defined in read_poscar)
                                                ! The naming is somewhat confusion since it can be misreading as number of eigen states but
                                                ! we just let it be for the legacy. It represents TOTAL NUMBER OF ORBITAL BASIS.
       integer*4                     neig_total ! neig_up + neig_dn
       integer*4                     neig_target
       integer*4                     nbasis  ! normally neig = nbasis
       integer*4                     neig_eff   ! number of orbital basis for the effective ham (for up and dn)

       integer*4                     init_erange, fina_erange  ! ie:fe, set by ERANGE tag, otherwise, 1 and neig*ispinor 
       integer*4                     nband ! nonmag: ie:fe, collinear: ie:fe for up or dn, non-collinear: ie:fe
                                           ! = PGEOM%neig*PINPT%ispinor  (if .not. PINPT%flag_erange, default)
                                           ! = PINPT%feast_nemax (if EWINDOW tag is on and nemax is smaller than PGEOM%neig*PINPT%ispinor)
                                           ! = PGEOM%fina_erange - PGEOM%init_erange + 1 ( if PINPT%flag_erange)

       integer*4                     n_spec, n_atom
       integer*4,   allocatable   :: i_spec(:) ! number of atoms per each species (1:n_spec)
       character*8, allocatable   :: c_spec(:) ! character of species for each species (1:n_spec)
       logical,     allocatable   :: flag_site_cindex(:) ! flag site indicator has been defined or not
       character*20,allocatable   :: site_cindex(:)  ! site indicator. NOTE: same as site_cindex of NN_TABLE 
       integer*4,   allocatable   :: spec(:)        ! species information for each atom (1:n_atom). The order of appearance in the POSCAR
       integer*4,   allocatable   :: spec_equiv(:)  ! species information specifying equivalent atom species
       real*8,      allocatable   :: a_coord(:,:) ! atomic  coordinate (1:3, 1:n_atom) (direct, fractional)
       real*8,      allocatable   :: a_coord_cart(:,:) ! atomic  coordinate (1:3, 1:n_atom) (cartesian)
       real*8,      allocatable   :: o_coord(:,:) ! orbital coordinate (1:3, 1:neig) (direct, fractional)
       real*8,      allocatable   :: o_coord_cart(:,:) ! orbital coordinate (1:3, 1:neig) (cartesian)
       integer*4,   allocatable   :: i_eff_orb(:) ! matrix index (diagonal) for effective hamiltonian setup (1:neig_eff * nspin)
                                                  ! this is only used for constructing NN_TABLE array

       integer*4                     max_orb  ! maximum number of orbitals asigned in each atomic site
       integer*4,   allocatable   :: n_orbital(:) ! number of orbitals per atomic site
       character*8, allocatable   :: c_orbital(:,:) ! name of atomic orbitals for each atomic sites
       real*8,      allocatable   :: orb_sign(:,:)  ! sign of orbital

       integer*4,   allocatable   :: ispec(:)     ! atomic number. As you can find in the periodic table.
       integer*4,   allocatable   :: l_quantum(:,:) ! angular momentum quantum number: s p d ..?
       integer*4,   allocatable   :: orb_n_quantum(:,:) ! principal quantum number for the orbital of l_quantum: 1s, 2s ?? ...
       real*8,      allocatable   :: n_quantum(:) ! principal quantum number for the atomic species
       real*8,      allocatable   :: z_eff_nuc(:,:) ! effective nuclear charge for each atomic orbital


       character*40                  system_name
       real*8                        a_scale, a_latt(3,3) ! lattice vector (unit of Ang)
       real*8                        b_latt(3,3) ! reciprocal lattice vector(unit of A^-1)
       logical                       flag_selective, flag_direct, flag_cartesian

       integer*4                     n_nn_type
       character*80,allocatable   :: nn_pair(:)
       real*8,      allocatable   :: nn_dist(:), nn_r0(:)

       integer*4,   allocatable   :: orb_index(:) ! store lm(atomic orbital, s, px,py...), from 1 to 9, information for each basis (neig)

       ! setup for ribbon geometry
       logical                       flag_set_ribbon, flag_print_only_ribbon_geom
       integer*4                     ribbon_nslab(3)
       real*8                        ribbon_vacuum(3)

       ! setup grid for charge plot
       integer*4                     ngrid(3)
       integer*4                     stm_ngrid(3) 
       real*8                        r_origin(3)   ! direct coordinate for shift of the origin of chg file

       ! local charge/moment (same as NN_TABLE)
       real*8,      allocatable   :: local_charge(:) !i=ham_index(neig)
       real*8,      allocatable   :: local_moment(:,:) !(1:3,i) (1:3)=(mi, theta, phi), i=ham_idx(neig)
       real*8,      allocatable   :: local_moment_rot(:,:) !(1:3,i) (1:3)=(mx, my, mz)=-I*mi_dot_sigma, i=ham_idx

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
       integer*4                     mysystem   ! my system index

       character*132                 kfilenm            ! kpoint file
       character*132                 ribbon_kfilenm     ! kpoint file for ribbon geometry defined in 'SET RIBBON'
       character*8                   kline_type ! FHI-AIMS, VASP, FLEUR
       integer*4                     nkpoint,nline
       integer*4                     n_ndiv
       integer*4                     idiv_mode ! kpath division mode
                                               ! 1 -> division type: vasp-like. n-1 division between kpoint A and B and total n points
                                               ! 2 -> division type: fleur-like. n division between kpoint A and B and total n+1 points 
                                               ! 3 -> division type: vasp-like with n-1 division between each segments. In this mode, however,
                                               !                    every path has different division. This is same as FHI-AIMS code does.

       integer*4,   allocatable   :: ndiv(:)
       integer*4                     kreduce ! Should be 1 if kline_type is not FLEUR, in the current version.
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
       integer*4                     mysystem   ! my system index
       ! energy/wavefunction etc setup
       real*8,      allocatable   :: E(:,:)   !E(neig*ispin,nkpoint) (for nspin=2, up=1:neig, dn=neig+1:neig*2), eigenstate
       complex*16,  allocatable   :: V(:,:,:) !V(nbasis=neig*ispin,neig*ispin,nkpoint) wavevector (basis,eigen,kpooint)
       real*8,      allocatable   :: ORB(:,:,:) ! orbital projected density of state (lm, neig*ispin (or nband*nspin), nkpoint), lm=s,px,py....,dyz
       complex*16,  allocatable   :: SV(:,:,:) !SV(nbasis=neig*ispin,neig*ispin,nkpoint) overlap matrix S multiplied V, S_ij = <i|j> (i,j,kpoint)
       real*8,      allocatable   :: D(:,:,:)   !D(3,neig*ispin,nkpoint) Information for degeneracy is stored.
                                                !D(1,:,:) -> degeneracy D_above * D_below 
                                                !D(2,:,:) -> degeneracy D_above = E_n+1 - E_n  (exception, E_neig = 1)
                                                !D(3,:,:) -> degeneracy D_below = E_n   - E_n-1(exception, E_1    = 1)
       integer*4,   allocatable   :: IDX(:,:) ! band index after band re-ordering (valid if LORDER = .TRUE.)
       real*8,      allocatable   :: E_ORD(:,:) ! re-ordered band
       complex*16,  allocatable   :: V_ORD(:,:,:) ! re-ordered eigenvector
       complex*16,  allocatable   :: SV_ORD(:,:,:) ! re-ordered eigenvector multiplied with S

       real*8,      allocatable   :: E_BAND(:)  ! band  energy of the system for each spin (nspin) 
       real*8,      allocatable   :: E_TOT(:)   ! total energy of the system for each spin (nspin)
       real*8,      allocatable   :: F_OCC(:,:) ! Fermi-dirac occupation function (ispin*nband, nkp)
       real*8                        E_F        ! Fermi level
  endtype energy

  type weight !PWGHT
       integer*4                     mysystem   ! my system index
       character*132                 efilenmu,efilenmd  ! target energy file (spin up & dn)
       character*16                  efile_type         ! type of target_energy file: "VASP" (for future release we will add AIMS, QE ...) or  "user"
       real*8                        efile_ef           ! fermi level for efile (energy shift: energy - efile_ef will be applied when reading efile)
       integer*4                     itarget_e_start    ! from which energy level TBFIT read as target energy from target file?
       integer*4                     read_energy_column_index, read_energy_column_index_dn ! which column of the target file should be read?
       
       integer*4                     iband, fband
       integer*4                     vbmd, cbmd   ! valence and conduction band maximum and minimum: DFT
       integer*4                     vbmt, cbmt   ! valence and conduction band maximum and minimum: TBA
       real*8,      allocatable   :: WT(:,:)
       real*8,      allocatable   :: DEGENERACY_WT(:,:)
       real*8,      allocatable   :: PENALTY_ORB(:,:,:)
       logical                       flag_weight_default
       logical                       flag_weight_orb

       ! weight setting
       integer*4                     nweight      ! number of           weight  (WEIGHT)  strip in "SET WEIGHT" tag
       integer*4                     ndegenw      ! number of dgeneracy weight  (DEGENW)  strip in "SET WEIGHT" tag
       integer*4                     npenalty_orb ! number of orbital   penalty (PENALTY) strip in "SET WEIGHT" tag
       character*132, allocatable :: strip_kp(:), strip_tb(:), strip_df(:), strip_wt(:)
       character*132, allocatable :: strip_kp_orb(:), strip_tb_orb(:), strip_orb(:), strip_site(:), strip_pen_orb(:)
       character*132, allocatable :: strip_kp_deg(:), strip_tb_deg(:), strip_df_deg(:), strip_wt_deg(:)
       integer*4                     max_len_strip_kp, max_len_strip_tb, max_len_strip_df
       integer*4                     max_len_strip_ob, max_len_strip_st

  endtype weight

  type hopping !NN_TABLE (nearest neighbor table; but not restricted to nearest)
       integer*4                     mysystem   ! my system index

       character*132                 nnfilenm      ! hopping integral file name (default = hopping.dat)     
       character*132                 nnfilenmo     ! hopping integral file name (default = hopping.dat)     
       integer*4,   allocatable   :: i_atom(:)     !(n) n = nn_neighbor index
       integer*4,   allocatable   :: j_atom(:)
       real*8,      allocatable   :: i_coord(:,:)  ! cartesian coordinate of atom_i (1:3,n) = nn_neighbor index
       real*8,      allocatable   :: j_coord(:,:)  ! cartesian coordinate of atom_j
       real*8,      allocatable   :: Rij(:,:) ! vector connecting two orbital i and j
       real*8,      allocatable   :: R(:,:)   ! cell periodicity where orbital j sits on
       real*8,      allocatable   :: Dij(:)   ! distance between orbital i and j
       real*8,      allocatable   :: Dij0(:)
       integer*4,   allocatable   :: i_matrix(:) ! matrix index (1:neig)
       character*8, allocatable   :: ci_orb(:)   ! orbital character for n-th nn_neighbor hopping
       real*8,      allocatable   :: i_sign(:)   ! orbital sign      for n-th nn_neighbor hopping
       integer*4,   allocatable   :: j_matrix(:) ! matrix index (1:neig) 
       character*8, allocatable   :: cj_orb(:)   ! orbital character for n-th nn_neighbor hopping
       real*8,      allocatable   :: j_sign(:)   ! orbital sign      for n-th nn_neighbor hopping
       character*2, allocatable   :: p_class(:)
       integer*4,   allocatable   :: n_class(:)
       integer*4,   allocatable   :: sk_index_set(:,:) ! (i,nn), i=onsite(0),sigma(1),  pi(2),  delta(3)
                                                       !           scaled->s_sigma(4),s_pi(5),s_delta(6)
       integer*4,   allocatable   :: cc_index_set(:,:) ! (i,nn), i=onsite(0),t_typ(1),t_typ(2),t_typ(3), t_typ(...)
       real*8,      allocatable   :: tij(:) ! hopping amplitude (except SOC, magnetic coupling)
       real*8,      allocatable   :: sij(:) ! overlap integral if USE_OVERLAP activated in PARAM_FIT.dat
       real*8                        onsite_tolerance
       integer*4                     n_neighbor ! total number of pairwise two-center hoppings
       integer*4,   allocatable   :: n_nn(:) !(natom) total number of nearest neighbor pairs for each atoms within cutoff radius (valid if SK_SCALE_TYPE >= 11)
       real*8,      allocatable   :: R_nn(:,:)  ! (n_nn, n_atom) distance information for nearest neighbor pairs. cutoff distance d0_cut
       real*8,      allocatable   :: R0_nn(:,:) ! (n_nn, n_atom) distance information for nearest neighbor pairs. reference cutoff distance d0
       integer*4,   allocatable   :: j_nn(:,:) ! (n_nn, n_atom) atom index for nearest neighbor pairs
       integer*4                     max_nn_pair  ! max(n_nn)
       integer*4,   allocatable   :: l_onsite_param_index(:) ! array (n_atom)

       character*20,allocatable   :: site_cindex(:)  ! site indicator
       logical,     allocatable   :: flag_site_cindex(:) ! flag site indicator has been defined or not
       real*8,      allocatable   :: local_charge(:) !i=ham_index(neig)
       real*8,      allocatable   :: local_moment(:,:) !(1:3,i) (1:3)=(mi, theta, phi), i=ham_idx(neig)
       real*8,      allocatable   :: local_moment_rot(:,:) !(1:3,i) (1:3)=(mx, my, mz)=-I*mi_dot_sigma, i=ham_idx
       integer*4,   allocatable   :: stoner_I_param_index(:) ! array size = i, i=ham_index(neig), 
       integer*4,   allocatable   :: local_U_param_index(:) ! array size = i, i=ham_index(neig), 
       integer*4,   allocatable   :: plus_U_param_index(:) ! array size = i, i=ham_index(neig), 
       integer*4,   allocatable   :: soc_param_index(:)    ! soc parameter index for each orbital-orbital pair (nn)

       real*8,      allocatable   :: tij_file(:) ! hopping amplitude read from file (hopping.dat)
       real*8,      allocatable   :: sij_file(:) ! overlap integral  read from file (overlap.dat)

       integer*4,   allocatable   :: i_eff_orb(:) ! effective orbital index (1:neig * ispin) ! only meaningful if ERANGE=full

       ! E-field
       logical                       flag_efield, flag_efield_frac, flag_efield_cart
       real*8                        efield(3)
       real*8                        efield_origin(3)
       real*8                        efield_origin_cart(3)
  endtype hopping

  type dos ! PINPT_DOS
       integer*4                     mysystem   ! my system index

       character*40                  dos_kfilenm
       character*40                  dos_filenm, ldos_filenm
       integer*4                     dos_kgrid(3)
       integer*4                     dos_nediv 
       integer*4                     dos_iband, dos_fband
       real*8                        dos_emin, dos_emax
       real*8                        dos_smearing
       real*8                        dos_kshift(1:3)
       real*8,      allocatable   :: dos_kpoint(:,:)       
       real*8,      allocatable   :: dos_erange(:)   ! size=nediv
       real*8,      allocatable   :: dos_tot(:,:)    ! size=nspin,nediv
       logical                       dos_flag_gamma, dos_flag_print_kpoint
       logical                       dos_flag_print_ldos
       logical                       dos_flag_print_eigen

       integer*4,   allocatable   :: dos_ensurf(:) !integer array of band index. size=n_ensurf
       integer*4                     dos_n_ensurf
       integer*4,   allocatable   :: dos_ldos_atom(:) ! integer array of atom index. maxsize=n_atom
       integer*4                     dos_ldos_natom ! how many atoms to be plotted for ldos
       real*8,      allocatable   :: ldos_tot(:,:,:,:) !LDOS for each (ORB, ATOM, ENERGY, SPIN)

       character*1                   dos_kunit

       logical                       dos_flag_sparse
  endtype dos

  type berry ! PINPT_BERRY
       integer*4                     mysystem   ! my system index

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
       integer*4                     noccupied
       logical                       flag_print_hamiltonian_parity, flag_parity_phase

       real*8                        symmetry_origin(3)
       real*8                        symmetry_operator(3,3)
       integer*4                     symmetry_nkpoint
       real*8,      allocatable   :: symmetry_kpoint(:,:)
       real*8,      allocatable   :: symmetry_kpoint_reci(:,:)
       character*10,allocatable   :: symmetry_kpoint_name(:)
       real*8                        symmetry_theta ! rotation angle along z-axis (definite specification)
       logical                       flag_print_hamiltonian_symmetry, flag_symmetry_phase
 

  endtype berry

  type spmat ! sparse matrix with Compressed Sparse Row format
       integer*4                     nnz     ! number of non-zero elements
       integer*4                     msize   ! matrix size
       complex*16,  allocatable   :: H(:)    ! sparse array of square matrix H_square(msize,msize), (n_neighbor)
       integer*4,   allocatable   :: I(:)    ! Row    array : if CSR format : I(msize + 1), I(msize+1)-1 = nnz
                                             !                if COO format : I(nnz)
       integer*4,   allocatable   :: J(:)    ! Column array : J(nnz), same size with H
  endtype spmat

  type gainp ! input/control parameters for genetic algorithm
       integer*4                     mysystem   ! my system index

       integer*4                     mgen    ! maximum number of GA iterations (generation) (default : 500)
       integer*4                     npop    ! number of individuals in a population (default is 100)
       integer*4                     ngene   ! number of significant digits (number of genes) retained 
                                             !           in chromosomal encoding (default is 6).
       real*8                        pcross  ! crossover probability ; must be  <= 1.0 (default : 0.85)
       real*8                        pmutmn  ! minimum mutation rate; must be >= 0.0 (default is 0.0005)
       real*8                        pmutmx  ! mmaximum mutation rate; must be <= 1.0 (default is 0.25)
       real*8                        pmut    ! initial mutation rate; should be small (default is 0.005)
       integer*4                     imut    ! mutation mode; 1/2/3/4/5 (default is 2).
                                             ! 1=one-point mutation, fixed rate.
                                             ! 2=one-point, adjustable rate based on fitness.
                                             ! 3=one-point, adjustable rate based on distance.
                                             ! 4=one-point+creep, fixed rate.
                                             ! 5=one-point+creep, adjustable rate based on fitness.
                                             ! 6=one-point+creep, adjustable rate based on distance.
       real*8                        fdif    ! relative fitness differential; range from 0(none) to 1(maximum). (default 1)
       integer*4                     irep    ! reproduction plan; 1 = Full generational replacement
                                             !                    2 = Steady-state-replace-random
                                             !                    3 = Steady-state-replace-worst (default)
       integer*4                     ielite  ! elitism flag; 0/1=off/on (default is 0)
                                             !               (Applies only to reproduction plans 1 and 2)
       integer*4                     ivrb    ! printed output 0/1/2=None/Minimal/Verbose  (default : 0)
       real*8                        convtol ! convergence tolerance; must be > 0.0 (default is 0.0001)
       integer*4                     convwin ! convergence window; must be >= 0
                                             ! This is the number of consecutive solutions
                                             ! within the tolerance for convergence to be declared (default is 20)
       real*8                        iguessf ! fraction of the initial population to set equal to the initial guess.
                                             ! (none) to 1.0 (all). (default is 0.1 or 10%).
       integer*4                     iseed   ! random seed value; must be > 0 (default is 999)
       real*8                        lower_bound, upper_bound
  endtype gainp

  type replot ! PRPLT, required parameters for replot dos/ldos with given bandstructure file
       integer*4                     mysystem   ! my system index

       logical                       flag_replot            ! whether turn on replot mode

       logical                       flag_replot_dos        ! flag for density of states
       logical                       flag_replot_ldos       ! flag for local density of states
       logical                       flag_replot_sldos      ! flag for spatial local density of states
       logical                       flag_replot_band       ! flag for replot band structure
       logical                       flag_replot_didv       ! flag for replot spatial density of states (LDOS) for given energy window
       logical                       flag_replot_proj_band  ! flag for replot projected band

       logical                       flag_replot_formatted ! flag whether band_structure_TBA file is formatted (.dat) or unformatted (.bin)

       character*80                  fname_band_up, fname_band_dn  ! file name the be read (set to default)
       real*8                        replot_dos_smearing ! gaussian smearing
       real*8,      allocatable   :: replot_dos_erange(:)   ! size=nediv, division of energy 
       integer*4                     replot_dos_nediv     ! number of division
       real*8                        replot_dos_emin, replot_dos_emax  ! from emin to emax
       real*8,      allocatable   :: replot_dos_tot(:,:)     ! density of states, (nspin, nediv)
       real*8,      allocatable   :: replot_dos_ntot(:,:) ! integrated DOS from initial to the energy level (nspin,0:nediv)

       integer*4                     replot_nldos_sum  ! number of ldos plot with set of atoms (how many REPLOT_LDOS tag given?)
       integer*4,   allocatable   :: replot_ldos_natom(:) ! number of atoms to be resolved for each REPLOT_LDOS request (nldos_sum)
       integer*4,   allocatable   :: replot_ldos_atom(:,:) ! atom index
       real*8,      allocatable   :: replot_ldos_tot(:,:,:,:,:)   ! local DOS (n_orbital(iatom), maxval(ldos_natom), nspin, nediv,nldos_sum)

       character*40                  replot_sldos_fname
       character*40                  replot_dos_fname
       character*40                  replot_didv_fname
       real*8,      allocatable   :: replot_sldos_sum(:,:,:) ! spatial local density of states sum over( natom, nspin, nediv )
       integer*4                     replot_sldos_cell(3)  ! repeat cell along a1, a2, a3 direction
       real*8                        r_origin(3) ! direct coordinate for shift of the origin of atomic coordinate of SLDOS
       real*8                        bond_cut    ! bond length <= bond_cut will not be written in BOND.replot.dat Default: 3.0 (ang)


       integer*4                     replot_nband ! number of REPLOT_BAND request
       character*2, allocatable   :: replot_axis_print_mag(:) ! (nband)
       logical,     allocatable   :: flag_replot_print_single(:) ! (nband)
       logical,     allocatable   :: flag_replot_write_unformatted(:) ! (nband)

       ! projected band
       integer*4                     replot_nproj_sum
       integer*4,   allocatable   :: replot_proj_natom(:) ! how many atoms to be plotted for projected band
       integer*4,   allocatable   :: replot_proj_atom(:,:) ! integer array of atom index. maxsize=n_atom
     
  endtype replot
endmodule 
