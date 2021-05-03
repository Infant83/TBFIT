module parameters
  use mpi_setup
  use iso_c_binding, only:  c_char, c_int, c_double, c_ptr, c_f_pointer
  use mykind
#ifdef MKL_SPARSE
  use MKL_SPBLAS 
#endif
  character(len=26), public, parameter ::alphabet='abcdefghijklmnopqrstuvwxyz'
  real(kind=dp)   , public, parameter   ::      pi=4.d0*atan(1.d0) ! Leibniz's formula for Pi
  real(kind=dp)   , public, parameter   ::     pi2=2.d0*pi
  real(kind=dp)   , public, parameter   ::    bohr=0.52917721067d0 ! meter Bohr radius
  real(kind=dp)   , public, parameter   :: hartree=27.21138602d0 ! eV   Hartree energy
  real(kind=dp)   , public, parameter   ::    hbar=4.135667662d-15/pi2  ! eV*s Plank constant
  real(kind=dp)   , public, parameter   ::       c=0.262465831d0 !! constant c = 2m/hbar**2 [1/eV Ang^2] 
  real(kind=dp)   , public, parameter   ::     k_B=8.6173303d-5  !! Boltzmann constant = 2m/hbar**2 [1/eV Ang^2] 
  real(kind=dp)   , public, parameter   :: g_elect=2.0023193043617 !! g-factor[see Mol.Phys. 98, 1597(2000) for sign]
  real(kind=dp)   , public, parameter   ::     rt2=sin( 4.d0*atan(1.d0)/4.d0 ) * 2.d0 ! sqrt(2)
  real(kind=dp)   , public, parameter   ::     rt3=sin( 4.d0*atan(1.d0)/3.d0 ) * 2.d0 ! sqrt(3)
  real(kind=dp)   , public, parameter   :: onsite_tolerance= 0.0001 !! symmetry precision
  real(kind=dp)   , public, parameter   ::     eta=epsilon(1d0) ! tiny value
  real(kind=dp)   , public              :: t1, t0 ! time
  complex(kind=dp), public, parameter   ::      zi=(0.d0,1.d0)
  complex(kind=dp), public, parameter   ::     pzi= pi*zi
  complex(kind=dp), public, parameter   ::    pzi2=2*pzi
  complex(kind=dp), public, dimension(2,2),   parameter :: pauli_0 = reshape((/  1,  0,  0,  1 /), (/2,2/))
  complex(kind=dp), public, dimension(2,2),   parameter :: pauli_x = reshape((/  0,  1,  1,  0 /), (/2,2/))
  complex(kind=dp), public, dimension(2,2),   parameter :: pauli_y = reshape((/  0,  1, -1,  0 /) * zi, (/2,2/))
  complex(kind=dp), public, dimension(2,2),   parameter :: pauli_z = reshape((/  1,  0,  0, -1 /), (/2,2/))
  integer(kind=sp), public, dimension(3,3,3), parameter :: levi_civita = reshape((/0,0,0, 0,0,-1, 0,1,0, & 
                                                                     0,0,1, 0,0,0, -1,0,0, & 
                                                                     0,-1,0, 1,0,0, 0,0,0/), (/3,3,3/))
  integer(kind=sp),  public, dimension(3,2),   parameter :: cyclic_axis = reshape((/2,3,1,3,1,2/), (/3,2/))

! integer(kind=sp),  public              :: iverbose 
  integer(kind=sp),  public, parameter   :: max_nparam      = 500    !! maximum number of onsite and hopping parameters in total
  integer(kind=sp),  public, parameter   :: max_kpoint      = 100000 !! maximum number of kpoints
  integer(kind=sp),  public, parameter   :: max_set_weight  = 100000  !! maximum number of "SET WEIGHT" tag
  integer(kind=sp),  public, parameter   :: max_pair_type   = 1000    !! maximum number of NN_CLASS pair type
  integer(kind=sp),  public, parameter   :: max_dummy       = 9999999 !! maximun number of dummy index for arbitral purpose   
  integer(kind=sp),  public, parameter   :: max_dummy2      = 1000    !! maximun number of dummy index for arbitral purpose   
  integer(kind=sp),  public, parameter   :: max_nsym        = 1000000 !! maximun number of symmetry operation for SPGLIB
  integer(kind=sp),  public, parameter   :: pid_energy      = 30 
  integer(kind=sp),  public, parameter   :: pid_nntable     = 33
  integer(kind=sp),  public, parameter   :: pid_incar       = 78
  integer(kind=sp),  public, parameter   :: pid_kpoint      = 79 
  integer(kind=sp),  public, parameter   :: pid_stm         = 81 
  integer(kind=sp),  public, parameter   :: pid_chg         = 82 
  integer(kind=sp),  public, parameter   :: pid_ibzkpt      = 83 
  integer(kind=sp),  public, parameter   :: pid_dos         = 84 
  integer(kind=sp),  public, parameter   :: pid_ldos        = 85 
  integer(kind=sp),  public, parameter   :: pid_param       = 88 
  integer(kind=sp),  public, parameter   :: pid_matrix      = 99
  integer(kind=sp),  public, parameter   :: pid_berrycurv   = 100
  integer(kind=sp),  public, parameter   :: pid_zak         = 101
  integer(kind=sp),  public, parameter   :: pid_wcc         = 102 ! pid_wcc + 1 = pid_gap
  integer(kind=sp),  public, parameter   :: pid_geom        = 177 
  integer(kind=sp),  public, parameter   :: pid_geom_ribbon = 178
  integer(kind=sp),  public, parameter   :: pid_circ        = 179
 
  type incar !PINPT
       ! set in parsing step
       integer(kind=sp)                     nsystem        ! number of geometry
       character(len=132), allocatable :: ifilenm(:)     ! input tag file (default = INCAR-TB, specified by -input ifilenm in command-line)
       character(len=132), allocatable :: title(:)       ! title of the system
       character(len=132)                 pfilenm_parse  ! parsed parameter file name
       character(len=132)                 kfilenm_parse  ! parsed kpoints   file name
       logical                       flag_python_module ! .true. if call from python module tbfitpy

       character(len=40)                  fnamelog       ! log file name, default = TBFIT.log
       logical                       flag_get_band  ! default = .true.
       logical                       flag_spglib    ! default = .true. ! write space group information
       logical                       flag_tbfit_parse, flag_tbfit_parse_
       logical                       flag_kfile_parse
       logical                       flag_pfile_parse
       logical                       flag_ndiv_line_parse, flag_ndiv_grid_parse
       logical                       flag_fit_degeneracy ! fitting is helped by fitting degeneracy as well.
       logical                       flag_miter_parse, flag_mxfit_parse
       logical                       flag_lorbit_parse
!      logical                       flag_proj_parse
       logical                       flag_parse
       logical                       flag_tbfit_test
       logical                       flag_inputcard_fname_parse
       logical                       flag_ga_with_lmdif ! use lmdif method for local optimization in GA method, default = .false.
       logical                       flag_pso_with_lmdif ! use lmdif method for local optimization in PSO method, default = .false.
       logical                       flag_pso_report_particles
       logical                       flag_write_unformatted ! default = .false.
       logical                       flag_report_geom
   
       logical                       flag_phase  ! default = .true.  ! apply phase factor to atomic orbital

       real(kind=dp)                        ptol
       real(kind=dp)                        ftol
       real(kind=dp)                        fdiff
       integer(kind=sp)                     ga_npop
       integer(kind=sp)                     miter,mxfit
       logical                       flag_tbfit, flag_pincar
       logical                       flag_tbfit_finish ! when exiting fitting routine, it will be turn to .true.
       logical                       flag_print_only_target
       logical                       flag_print_energy_diff ! also print energy different between target and tight-binding energies
       logical                       flag_print_orbital  ! activated with LORBIT .TRUE. tag
       logical                       flag_get_orbital    ! whether request wavevector in diagonalize routine
       logical                       flag_print_mag
       logical                       flag_print_single ! default .false. (write single precision for wavefunction)
       logical                       flag_load_nntable ! default .false.
       character(len=2)                   axis_print_mag ! mx, my, mz (pauli matrices), 
                                                    ! re, im (real or imag part), 
                                                    ! wf (full wf), bi (enforce to write wf with binary format)

       integer(kind=sp)                    npar      ! number of group for k-point parallelization
       character(len=10)                   ls_type   ! fitting method
       integer(kind=sp)                     nn_max(3)          ! cell reapeat for the nearest neighbor finding in find_nn routine (default:3 3 3)
       logical                       flag_use_weight    ! if true use "weight" information written in PFILE and 
                                                        ! replace it with SET WEIGHT info of INCAR-TB file after fitting procedures

       ! plot_eig mode
       logical                       flag_plot_eigen_state  !default : .false.
       logical                       flag_plot_wavefunction !default : .true. => WAV_PLOT = .true. if .not. CHG_PLOT = .true.
       integer(kind=sp)                     n_eig_print,n_kpt_print
       integer(kind=sp),    allocatable  :: i_eig_print(:)
       integer(kind=sp),    allocatable  :: i_kpt_print(:)
       real(kind=dp)                        rcut_orb_plot ! cutoff radius for orbital plot (applies STM mode too if provied in SET STMPLOT)
                                                   ! default = 5 Ang (see PRB 64, 235111)
       ! plot_stm mode
       logical                       flag_repeat_cell_orb_plot(3) ! default = .true. true. true.
       integer(kind=sp)                     repeat_cell_orb_plot(3) ! default = 1
       logical                       flag_plot_stm_image    !default : .false.
       integer(kind=sp)                     n_stm        ! how many stm images will be plotted
       real(kind=dp),       allocatable  :: stm_emax(:), stm_emin(:)

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
       integer(kind=sp)              ispin   ! nonmag: 1, collinear: 2, non-collinear: 2
       integer(kind=sp)              ispinor ! nonmag: 1, collinear: 1, non-collinear: 2
       integer(kind=sp)              nspin   ! nonmag: 1, collinear: 2, non-collinear: 1
       logical                       flag_local_charge, flag_collinear, flag_noncollinear, flag_soc

       logical                       flag_erange
       logical                       flag_sparse ! if EWINDOW tag has been set up, flag_sparse is forced to be .true.
       integer(kind=sp)              feast_nemax ! maximum number of eigenvalues (=nband_guess)
       real(kind=dp)                 feast_emin, feast_emax ! energy window [emin:emax] for FEAST algorithm

       ! NOTE: new version of FEAST has fpm array with 64 but old version have 128. Should be properly chosen according to the FEAST LIB
       integer(kind=sp)              feast_fpm(128)   ! FEAST parameters set by original feastinit
       integer(kind=sp)              feast_fpm64(64)    ! FEAST parameters set by original feastinit
       integer(kind=sp)              feast_fpm128(128)  ! FEAST parameters set by MKL extended eigensolver feast support routine
!#ifdef PSPARSE
!       integer(kind=sp)              feast_fpm(64)  ! FEAST parameters set by original feastinit
!#else
!       integer(kind=sp)              feast_fpm(128) ! FEAST parameters set by MKL extended eigensolver feast support routine
!#endif
       integer(kind=sp),    allocatable  :: feast_ne(:,:)  ! Number of states found in erange [emin:emax], (:,:) = (nspin,nkpoint)
       integer(kind=sp)              feast_neguess  ! initial guess for the number of states to be found in interval

       logical                       flag_plus_U

       logical                       flag_scissor
       integer(kind=sp)              i_scissor ! i_scissor level will be operated by scissor operator 
       real(kind=dp)                 r_scissor ! EDFT(n,k) + r_scissor if n >= i_scissor

       ! construct effective hamiltonian
       logical                       flag_get_effective_ham
       character(len=132)            eff_orb_dummyc
       real(kind=dp)                 eff_emin, eff_emax

       ! projected band
       integer(kind=sp)                     nproj_sum
       integer(kind=sp),   allocatable   :: proj_atom(:,:) ! integer array of atom index. maxsize=n_atom
       integer(kind=sp),   allocatable   :: proj_natom(:) ! how many atoms to be plotted for projected band
       logical                       flag_print_proj
       logical                       flag_print_proj_sum
       logical                       flag_print_proj_atom

       logical                       flag_plot_fit
       logical                       flag_plot
       character(len=132)            filenm_gnuplot, filenm_gnuplot_parse
       logical                       flag_filenm_gnuplot_parse

       logical                       flag_print_energy_singlek
       logical                       flag_print_hamk

       logical                       flag_get_band_order   ! flag whether perform band re-ordering by calculating overlap integral
       logical                       flag_get_band_order_print_only ! if true, band re-order will not be performed in the fitting routines
       real(kind=dp)                 band_order_overlap_cutoff ! cutoff of overlap integral to perform eigenvalue swap : sqrt(2)/2 by default

       logical                       flag_get_total_energy 
       real(kind=dp)                 electronic_temperature ! Temperature, default = 0 (K)

       ! circular dichroism 
       logical                       flag_get_circ_dichroism ! flag whether perform circular dichroism calculation
       integer(kind=sp)                     ncirc_dichroism         
       integer(kind=sp),   allocatable   :: circ_dichroism_pair(:,:)

       integer(kind=sp)                     lmmax ! 9 (default) =>  1:s, 2~4:px,py,pz, 5~9:dz2,dx2,dxy,dxz,dyz
                                           ! 3           =>  1:s, 2  :px,py,pz, 3  :dz2,dx2,dxy,dxz,dyz
       logical                       flag_fit_orbital ! fit orbital character as well?
       real(kind=dp)                        orbital_fit_smearing ! gaussian smearing for orbital fitting                                           

       logical                              flag_distribute_nkp ! distribute V and SV over cpu nodes (my_nkp). Only valid in post-processing routine

       integer(kind=sp)              iseed ! random seed
       logical                       flag_pso_verbose_parse
       integer(kind=sp)              pso_verbose 
       
  endtype incar

  type params !PPRAM
     ! integer(kind=sp)                     mysystem   ! my system index (deprecated...)
       logical                       flag_set_param_const     ! whether set constrain on parameters
       logical                       flag_pfile_index         ! whether put numbering in PARAM_FIT.dat after fitting
       logical                       flag_use_overlap         ! whether setup overlap hamiltonian. Automatically activated if the overlap integral parameters
                                                              ! are provided in the PARAM_FIT.dat file in priori (start with o_ )
       logical                       flag_fit_plain           ! take costfunction without multiplying weight factor for PSO scheme
       real(kind=dp)                 l_broaden                ! broadening of the cutoff-function for NRL type SK parameters (used if SK_SCALE_MODE > 10 in PFILE)
       logical                       flag_slater_koster       ! default .true.
       logical                       flag_nrl_slater_koster   ! default .false.
       integer(kind=sp)              slater_koster_type       ! 1 ~ 3: exponential scaling, 11: Mehl & Papaconstantopoulos NRL method (PRB 54, 4519 (1996))
       character(len=132)            pfilenm,pfileoutnm       ! SK parameter file input & output
       integer(kind=sp)                     nparam                   ! total number of parameters 
       integer(kind=sp)                     nparam_const             ! total number of parameters constraints
       integer(kind=sp)                     nparam_free              ! total number of free parameters size(iparam_free(:))
       real(kind=dp),       allocatable  :: param(:)                 ! TB parameters
       real(kind=dp),       allocatable  :: param_nrl(:,:)           ! NRL TB parameters
       integer(kind=sp),    allocatable  :: iparam_free(:)           ! parameter index for free parameters only (nparam_free) 
       integer(kind=sp),    allocatable  :: iparam_free_nrl(:)       ! parameter index for free parameters only (nparam_free), sum(param_nsub(iparam_free(1:j-1)))+1
       character(len=40), allocatable  :: param_name(:)            
       character(len=40), allocatable  :: c_const(:,:)             
       real(kind=dp),       allocatable  :: param_const(:,:)         ! i=1 -> 'is same as'
                                                              ! i=2 -> 'is lower than' (.le.) : maximum bound  ! <=  20 
                                                              ! i=3 -> 'is lower than' (.ge.) : minimum bound  ! >= -20 or >= 0.001 (if scale factor)
                                                              ! i=4 -> 'is fixed' : fixed during fitting       ! original value will be 
                                                                                                               ! copied to PINPT%param_const(i=5,:)
       real(kind=dp),       allocatable  :: param_const_nrl(:,:,:)   ! (i,j,nparam) ! j=1:4 -> a, b, c, d  for sk-parameters 
                                                              !                      -> alpha, beta, gamma, xi for onsite parameters
                                                              !                if size(,:,) = 1, then lonsite_ or lambda_, stoner_ , etc....
                                                              ! note: this parameter is only activated if slater_koster_type = 11
       integer(kind=sp)                     nparam_nrl               ! total number of parameters sum(param_nsub(:))
       integer(kind=sp)                     nparam_nrl_free          ! total number of free parameters sum(param_nsub(iparam_free(:)))
       integer(kind=sp)                     param_nsub_max           ! 4 if SK_SCALE_MODE >10, 1 if SK_SCALE_MODE <= 10
       integer(kind=sp),    allocatable  :: param_nsub(:)            ! (nparam) number of sub-parameters for each parameter
                                                              ! 4: NRL type of e_, ssp_, pds_,...(with s_, o_, os_ )-> onsite, hopping, hopping_scale, 
                                                              !                                               hopping(overlap),hopping_scale(overlap)
                                                              ! 1: other parameters -> local_U, lrashba_, lambda_, lsoc_, stoner_I_ , ..., etc.
       integer(kind=sp)                     niter ! number of iteractions performed

       integer(kind=sp)                     pso_miter  !  maximum number of iteration in PSO method
       character(len=40)                    pso_mode          ! 'pso':default method, 'pso_bestn': select bestN particles for next iteration
       integer(kind=sp)                     pso_nparticles    ! number of particles (parameters) used in particle swarm optimization method
       real(kind=dp)                        pso_c1            ! c1, c2, w defines cognitive, social, and inertia factor of PSO method  
       real(kind=dp)                        pso_c2            ! vel = w * vel + c1 * rand() * (pbest - x) + c2 * rand() * (gbest - x)
       real(kind=dp)                        pso_w 
       real(kind=dp)                        pso_report_ratio  ! (pso_report_ratio*100) percent of best particles would be reported
       integer(kind=sp)                     pso_iseed         ! random seed
       real(kind=dp)                        pso_max_noise_amplitude ! random noize amplitude
       real(kind=dp),      allocatable   :: pso_cost_history(:)  ! save cost function history w.r.t. the iteration in PSO method. size:(miter)
       real(kind=dp),      allocatable   :: pso_cost_history_i(:,:)  ! save cost function history w.r.t. the iteration in PSO method. size:(miter, pso_nparticles)
       real(kind=dp),      allocatable   :: pso_pbest_history(:,:,:) ! parameter history in PSO method (miter, pso_nparticles, nparam)
       real(kind=dp),      allocatable   :: cost_history(:)  ! save cost function history w.r.t. the iteration in LMDIF method. size:(miter)

  endtype params

  type poscar !PGEOM
       integer(kind=sp)                     mysystem   ! my system index

       character(len=132)                 title      ! system name (should have no blank)
       character(len=132)                 gfilenm    ! geometry file
       real(kind=dp),      allocatable   :: nelect(:)  ! number of total electrons in the system  (nspin, for up & dn if nspin=2), set by NELECT tag
                                                ! this parameter is used in total energy calculations to determine Fermi level of the system.
                                                ! NOTE(06.Oct.2020): The default will be automatically set in the future as following:
                                                ! by default this will be determined based on the valence electron configuration of the atom.
                                                ! For example, if you have two Si(3s1_3p3) atom and four H(1s1) atom, 
                                                ! then you will have 4*2 + 1*4 = 12 electrons in total with non-magnetic or noncollinear system 
                                                ! and 6 (up) and 6 (dn) electrons, respectively, with collinear magnetic systems.
       integer(kind=sp)                     neig       ! either neig_up and neig_dn, total number of atomic orbitals
                                                ! = sum(PGEOM%n_orbital(1:PGEOM%n_atom)) (defined in read_poscar)
                                                ! The naming is somewhat confusion since it can be misreading as number of eigen states but
                                                ! we just let it be for the legacy. It represents TOTAL NUMBER OF ORBITAL BASIS.
       integer(kind=sp)                     neig_total ! neig_up + neig_dn
       integer(kind=sp)                     neig_target
       integer(kind=sp)                     nbasis  ! normally neig = nbasis
       integer(kind=sp)                     neig_eff   ! number of orbital basis for the effective ham (for up and dn)

       integer(kind=sp)                     init_erange, fina_erange  ! ie:fe, set by ERANGE tag, otherwise, 1 and neig*ispinor 
       integer(kind=sp)                     nband ! nonmag: ie:fe, collinear: ie:fe for up or dn, non-collinear: ie:fe
                                           ! = PGEOM%neig*PINPT%ispinor  (if .not. PINPT%flag_erange, default)
                                           ! = PINPT%feast_nemax (if EWINDOW tag is on and nemax is smaller than PGEOM%neig*PINPT%ispinor)
                                           ! = PGEOM%fina_erange - PGEOM%init_erange + 1 ( if PINPT%flag_erange)

       integer(kind=sp)                     n_spec, n_atom
       integer(kind=sp),   allocatable   :: i_spec(:) ! number of atoms per each species (1:n_spec)
       character(len=8), allocatable   :: c_spec(:) ! character of species for each species (1:n_spec)
       logical,     allocatable   :: flag_site_cindex(:) ! flag site indicator has been defined or not
       character(len=20),allocatable   :: site_cindex(:)  ! site indicator. NOTE: same as site_cindex of NN_TABLE 
       integer(kind=sp),   allocatable   :: spec(:)        ! species information for each atom (1:n_atom). The order of appearance in the POSCAR
       integer(kind=sp),   allocatable   :: spec_equiv(:)  ! species information specifying equivalent atom species in the system (n_atom)
!      integer(kind=sp),   allocatable   :: spec_define(:)  ! species information specifying equivalent atom species in periodic table (n_atom)
       real(kind=dp),      allocatable   :: a_coord(:,:) ! atomic  coordinate (1:3, 1:n_atom) (direct, fractional)
       real(kind=dp),      allocatable   :: a_coord_cart(:,:) ! atomic  coordinate (1:3, 1:n_atom) (cartesian)
       real(kind=dp),      allocatable   :: o_coord(:,:) ! orbital coordinate (1:3, 1:neig) (direct, fractional)
       real(kind=dp),      allocatable   :: o_coord_cart(:,:) ! orbital coordinate (1:3, 1:neig) (cartesian)
       integer(kind=sp),   allocatable   :: i_eff_orb(:) ! matrix index (diagonal) for effective hamiltonian setup (1:neig_eff * nspin)
                                                  ! this is only used for constructing NN_TABLE array

       integer(kind=sp)                     max_orb  ! maximum number of orbitals asigned in each atomic site
       integer(kind=sp),   allocatable   :: n_orbital(:) ! number of orbitals per atomic site
       character(len=8), allocatable   :: c_orbital(:,:) ! name of atomic orbitals for each atomic sites
       real(kind=dp),      allocatable   :: orb_sign(:,:)  ! sign of orbital

       integer(kind=sp),   allocatable   :: ispec(:)     ! atomic number. As you can find in the periodic table.
       integer(kind=sp),   allocatable   :: l_quantum(:,:) ! angular momentum quantum number: s p d ..?
       integer(kind=sp),   allocatable   :: orb_n_quantum(:,:) ! principal quantum number for the orbital of l_quantum: 1s, 2s ?? ...
       real(kind=dp),      allocatable   :: n_quantum(:) ! principal quantum number for the atomic species
       real(kind=dp),      allocatable   :: z_eff_nuc(:,:) ! effective nuclear charge for each atomic orbital


       character(len=40)                  system_name
       real(kind=dp)                        a_scale, a_latt(3,3) ! lattice vector (unit of Ang)
       real(kind=dp)                        b_latt(3,3) ! reciprocal lattice vector(unit of A^-1)
       logical                       flag_selective, flag_direct, flag_cartesian

       integer(kind=sp)                     n_nn_type
       character(len=80),allocatable   :: nn_pair(:)
       real(kind=dp),      allocatable   :: nn_dist(:), nn_r0(:)

       integer(kind=sp),   allocatable   :: orb_index(:) ! store lm(atomic orbital, s, px,py...), from 1 to 9, information for each basis (neig)

       ! setup for ribbon geometry
       logical                       flag_set_ribbon, flag_print_only_ribbon_geom
       integer(kind=sp)                     ribbon_nslab(3)
       real(kind=dp)                        ribbon_vacuum(3)

       ! setup grid for charge plot
       integer(kind=sp)                     ngrid(3)
       integer(kind=sp)                     stm_ngrid(3) 
       real(kind=dp)                        r_origin(3)   ! direct coordinate for shift of the origin of chg file

       ! local charge/moment (same as NN_TABLE)
       real(kind=dp),      allocatable   :: local_charge(:) !i=ham_index(neig)
       real(kind=dp),      allocatable   :: local_moment(:,:) !(1:3,i) (1:3)=(mi, theta, phi), i=ham_idx(neig)
       real(kind=dp),      allocatable   :: local_moment_rot(:,:) !(1:3,i) (1:3)=(mx, my, mz)=-I*mi_dot_sigma, i=ham_idx

       !SPGLIB related variables
       integer(kind=sp)                     spg_error
       integer(kind=sp)                     spg_space_group !space group index
       integer(kind=sp)                     spg_hall_number
       integer(kind=sp)                     spg_Hermann_Mauguin_number
       real(kind=dp)                        spg_transformation_matrix(3,3)
       real(kind=dp)                        spg_origin_shift(3)
       integer(kind=sp)                     spg_n_operations ! nsym, number of symmetry operations
       character(len=12)            spg_international
       character(len=18)            spg_hall_symbol
       character(len=6)             spg_choice
       character(len=7)             spg_point_group
       character(len=12)            spg_crystal_system
       integer(kind=sp),   allocatable   :: spg_rotations(:,:,:)   ! {->w,   t} (3,3,spg_n_operations)
       real(kind=dp),      allocatable   :: spg_translations(:,:) !  {  w, ->t} (3,spg_n_operations)
       integer(kind=sp),   allocatable   :: spg_wyckoffs(:)
       integer(kind=sp),   allocatable   :: spg_equivalent_atoms(:) 
       real(kind=dp),      allocatable   :: spg_a_coord_operated(:,:,:) ! (3,n_atom,spg_n_operations)
       character(len=7)                     spg_schoenflies
       real(kind=dp)                        spg_a_latt_primitive(3,3)
       real(kind=dp),      allocatable   :: spg_a_coord_primitive(:,:) !(3,spg_n_atom_primitive)
       integer(kind=sp),   allocatable   :: spg_spec_primitive(:)      !(  spg_n_atom_primitive)
       integer(kind=sp)                     spg_n_atom_primitive
       integer(kind=sp),   allocatable   :: spg_det_w(:) ! (spg_n_operation) ! determinant of rotation w (1=proper, -1=improper)
       integer(kind=sp),   allocatable   :: spg_tr_w(:) ! (spg_n_operation) ! trace of rotation w 
       integer(kind=sp),   allocatable   :: spg_type_w(:) ! types of rotation operation of space group
  endtype poscar

  type kpoints !PKPTS
       integer(kind=sp)                     mysystem   ! my system index

       character(len=132)            kfilenm            ! kpoint file
       character(len=132)            ribbon_kfilenm     ! kpoint file for ribbon geometry defined in 'SET RIBBON'
       character(len=8)              kline_type ! FHI-AIMS, VASP, FLEUR
       integer(kind=sp)                     nkpoint,nline
       integer(kind=sp)                     n_ndiv
       integer(kind=sp)                     idiv_mode ! kpath division mode
                                               ! 1 -> division type: vasp-like. n-1 division between kpoint A and B and total n points
                                               ! 2 -> division type: fleur-like. n division between kpoint A and B and total n+1 points 
                                               ! 3 -> division type: vasp-like with n-1 division between each segments. In this mode, however,
                                               !                    every path has different division. This is same as FHI-AIMS code does.

       integer(kind=sp),   allocatable   :: ndiv(:)
       integer(kind=sp)                     kreduce ! Should be 1 if kline_type is not FLEUR, in the current version.
       real(kind=dp),      allocatable   :: kpoint(:,:),kline(:,:)
       real(kind=dp),      allocatable   :: kpoint_reci(:,:)
       real(kind=dp)                        k_shift(3)
       character(len=8),   allocatable   :: k_name(:)
       character(len=8),   allocatable   :: k_name2(:)
       integer(kind=sp),   allocatable   :: k_name_index(:)
       logical                       flag_klinemode
       logical                       flag_kgridmode, flag_gamma
       logical                       flag_reciprocal, flag_cartesianK
       character(len=1)            kunit
  endtype kpoints

  type energy !EDFT / ETBA / ETBA_DOS
       integer(kind=sp)                     mysystem   ! my system index
       ! energy/wavefunction etc setup
       real(kind=dp),      allocatable   :: E(:,:)   !E(neig*ispin,nkpoint) (for nspin=2, up=1:neig, dn=neig+1:neig*2), eigenstate
       real(kind=dp),      allocatable   :: dE(:,:)   !dE(neig*ispin,nkpoint) (for nspin=2, up=1:neig, dn=neig+1:neig*2), EDFT - ETBA, if flag_python_module
       complex(kind=dp),  allocatable   :: V(:,:,:) !V(nbasis=neig*ispin,neig*ispin,nkpoint) wavevector (basis,eigen,kpooint)
       real(kind=dp),      allocatable   :: ORB(:,:,:) ! orbital projected density of state (lm, neig*ispin (or nband*nspin), nkpoint), lm=s,px,py....,dyz
       complex(kind=dp),  allocatable   :: SV(:,:,:) !SV(nbasis=neig*ispin,neig*ispin,nkpoint) overlap matrix S multiplied V, S_ij = <i|j> (i,j,kpoint)
       real(kind=dp),      allocatable   :: D(:,:,:)   !D(3,neig*ispin,nkpoint) Information for degeneracy is stored.
                                                !D(1,:,:) -> degeneracy D_above * D_below 
                                                !D(2,:,:) -> degeneracy D_above = E_n+1 - E_n  (exception, E_neig = 1)
                                                !D(3,:,:) -> degeneracy D_below = E_n   - E_n-1(exception, E_1    = 1)
       integer(kind=sp),   allocatable   :: IDX(:,:) ! band index after band re-ordering (valid if LORDER = .TRUE.)
       real(kind=dp),      allocatable   :: E_ORD(:,:) ! re-ordered band
       complex(kind=dp),  allocatable   :: V_ORD(:,:,:) ! re-ordered eigenvector
       complex(kind=dp),  allocatable   :: SV_ORD(:,:,:) ! re-ordered eigenvector multiplied with S

       real(kind=dp),      allocatable   :: E_BAND(:)  ! band  energy of the system for each spin (nspin) 
       real(kind=dp),      allocatable   :: E_TOT(:)   ! total energy of the system for each spin (nspin)
       real(kind=dp),      allocatable   :: F_OCC(:,:) ! Fermi-dirac occupation function (ispin*nband, nkp)
       real(kind=dp)                        E_F        ! Fermi level
!      real(kind=dp),      allocatable   :: pso_cost_history(:)  ! save cost function history w.r.t. the iteration in PSO method. size:(miter)
  endtype energy

  type weight !PWGHT
       integer(kind=sp)                     mysystem   ! my system index
       character(len=132)            efilenmu,efilenmd  ! target energy file (spin up & dn)
       character(len=16)             efile_type         ! type of target_energy file: "VASP" (for future release we will add AIMS, QE ...) or  "user"
       real(kind=dp)                        efile_ef           ! fermi level for efile (energy shift: energy - efile_ef will be applied when reading efile)
       integer(kind=sp)                     itarget_e_start    ! from which energy level TBFIT read as target energy from target file?
       integer(kind=sp)                     read_energy_column_index, read_energy_column_index_dn ! which column of the target file should be read?
       
       integer(kind=sp)                     iband, fband
       integer(kind=sp)                     vbmd, cbmd   ! valence and conduction band maximum and minimum: DFT
       integer(kind=sp)                     vbmt, cbmt   ! valence and conduction band maximum and minimum: TBA
       integer(kind=sp)                     ie_cutoff    ! if(ie_cutoff .gt. 0) sum_n=1,ie_cutoff (EDFT_n - ETBA_n) , default = 0
       real(kind=dp),      allocatable   :: WT(:,:)
       real(kind=dp),      allocatable   :: DEGENERACY_WT(:,:)
       real(kind=dp),      allocatable   :: PENALTY_ORB(:,:,:)
       logical                       flag_weight_default
       logical                       flag_weight_orb

       ! weight setting
       integer(kind=sp)                     nweight      ! number of           weight  (WEIGHT)  strip in "SET WEIGHT" tag
       integer(kind=sp)                     ndegenw      ! number of dgeneracy weight  (DEGENW)  strip in "SET WEIGHT" tag
       integer(kind=sp)                     npenalty_orb ! number of orbital   penalty (PENALTY) strip in "SET WEIGHT" tag
       character(len=132), allocatable :: strip_kp(:), strip_tb(:), strip_df(:), strip_wt(:)
       character(len=132), allocatable :: strip_kp_orb(:), strip_tb_orb(:), strip_orb(:), strip_site(:), strip_pen_orb(:)
       character(len=132), allocatable :: strip_kp_deg(:), strip_tb_deg(:), strip_df_deg(:), strip_wt_deg(:)
       integer(kind=sp)                     max_len_strip_kp, max_len_strip_tb, max_len_strip_df
       integer(kind=sp)                     max_len_strip_ob, max_len_strip_st

  endtype weight

  type hopping !NN_TABLE (nearest neighbor table; but not restricted to nearest)
       integer(kind=sp)                     mysystem   ! my system index

       character(len=132)                 nnfilenm      ! hopping integral file name (default = hopping.dat)     
       character(len=132)                 nnfilenmo     ! hopping integral file name (default = hopping.dat)     
       integer(kind=sp),   allocatable   :: i_atom(:)     !(n) n = nn_neighbor index
       integer(kind=sp),   allocatable   :: j_atom(:)
       real(kind=dp),      allocatable   :: i_coord(:,:)  ! cartesian coordinate of atom_i (1:3,n) = nn_neighbor index
       real(kind=dp),      allocatable   :: j_coord(:,:)  ! cartesian coordinate of atom_j
       real(kind=dp),      allocatable   :: Rij(:,:) ! vector connecting two orbital i and j
       real(kind=dp),      allocatable   :: R(:,:)   ! cell periodicity where orbital j sits on
       real(kind=dp),      allocatable   :: Dij(:)   ! distance between orbital i and j
       real(kind=dp),      allocatable   :: Dij0(:)
       integer(kind=sp),   allocatable   :: i_matrix(:) ! matrix index (1:neig)
       character(len=8), allocatable   :: ci_orb(:)   ! orbital character for n-th nn_neighbor hopping
       real(kind=dp),      allocatable   :: i_sign(:)   ! orbital sign      for n-th nn_neighbor hopping
       integer(kind=sp),   allocatable   :: j_matrix(:) ! matrix index (1:neig) 
       character(len=8), allocatable   :: cj_orb(:)   ! orbital character for n-th nn_neighbor hopping
       real(kind=dp),      allocatable   :: j_sign(:)   ! orbital sign      for n-th nn_neighbor hopping
       character(len=2), allocatable   :: p_class(:)
       integer(kind=sp),   allocatable   :: n_class(:)
       integer(kind=sp),   allocatable   :: sk_index_set(:,:) ! (i,nn), i=onsite(0),sigma(1),  pi(2),  delta(3)
                                                       !           scaled->s_sigma(4),s_pi(5),s_delta(6)
       integer(kind=sp),   allocatable   :: cc_index_set(:,:) ! (i,nn), i=onsite(0),t_typ(1),t_typ(2),t_typ(3), t_typ(...)
       real(kind=dp),      allocatable   :: tij(:) ! hopping amplitude (except SOC, magnetic coupling)
       real(kind=dp),      allocatable   :: sij(:) ! overlap integral if USE_OVERLAP activated in PARAM_FIT.dat
       real(kind=dp)                        onsite_tolerance
       integer(kind=sp)                     n_neighbor ! total number of pairwise two-center hoppings
       integer(kind=sp),   allocatable   :: n_nn(:) !(natom) total number of nearest neighbor pairs for each atoms within cutoff radius (valid if SK_SCALE_TYPE >= 11)
       real(kind=dp),      allocatable   :: R_nn(:,:)  ! (n_nn, n_atom) distance information for nearest neighbor pairs. cutoff distance d0_cut
       real(kind=dp),      allocatable   :: R0_nn(:,:) ! (n_nn, n_atom) distance information for nearest neighbor pairs. reference cutoff distance d0
       integer(kind=sp),   allocatable   :: j_nn(:,:) ! (n_nn, n_atom) atom index for nearest neighbor pairs
       integer(kind=sp)                     max_nn_pair  ! max(n_nn)
       integer(kind=sp),   allocatable   :: l_onsite_param_index(:) ! array (n_atom)

       character(len=20),allocatable   :: site_cindex(:)  ! site indicator
       logical,     allocatable   :: flag_site_cindex(:) ! flag site indicator has been defined or not
       real(kind=dp),      allocatable   :: local_charge(:) !i=ham_index(neig)
       real(kind=dp),      allocatable   :: local_moment(:,:) !(1:3,i) (1:3)=(mi, theta, phi), i=ham_idx(neig)
       real(kind=dp),      allocatable   :: local_moment_rot(:,:) !(1:3,i) (1:3)=(mx, my, mz)=-I*mi_dot_sigma, i=ham_idx
       integer(kind=sp),   allocatable   :: stoner_I_param_index(:) ! array size = i, i=ham_index(neig), 
       integer(kind=sp),   allocatable   :: local_U_param_index(:) ! array size = i, i=ham_index(neig), 
       integer(kind=sp),   allocatable   :: plus_U_param_index(:) ! array size = i, i=ham_index(neig), 
       integer(kind=sp),   allocatable   :: soc_param_index(:)    ! soc parameter index for each orbital-orbital pair (nn)

       real(kind=dp),      allocatable   :: tij_file(:) ! hopping amplitude read from file (hopping.dat)
       real(kind=dp),      allocatable   :: sij_file(:) ! overlap integral  read from file (overlap.dat)

       integer(kind=sp),   allocatable   :: i_eff_orb(:) ! effective orbital index (1:neig * ispin) ! only meaningful if ERANGE=full

       ! E-field
       logical                       flag_efield, flag_efield_frac, flag_efield_cart
       real(kind=dp)                        efield(3)
       real(kind=dp)                        efield_origin(3)
       real(kind=dp)                        efield_origin_cart(3)
  endtype hopping

  type dos ! PINPT_DOS
       integer(kind=sp)                     mysystem   ! my system index

       character(len=40 )            dos_kfilenm
       character(len=40 )            dos_filenm, ldos_filenm
       integer(kind=sp)                     dos_kgrid(3)
       integer(kind=sp)                     dos_nediv 
       integer(kind=sp)                     dos_iband, dos_fband
       real(kind=dp)                        dos_emin, dos_emax
       real(kind=dp)                        dos_smearing
       real(kind=dp)                        dos_kshift(1:3)
       real(kind=dp),      allocatable   :: dos_kpoint(:,:)       
       real(kind=dp),      allocatable   :: dos_erange(:)   ! size=nediv
       real(kind=dp),      allocatable   :: dos_tot(:,:)    ! size=nspin,nediv
       logical                       dos_flag_gamma, dos_flag_print_kpoint
       logical                       dos_flag_print_ldos
       logical                       dos_flag_print_eigen

       integer(kind=sp),   allocatable   :: dos_ensurf(:) !integer array of band index. size=n_ensurf
       integer(kind=sp)                     dos_n_ensurf
       integer(kind=sp),   allocatable   :: dos_ldos_atom(:) ! integer array of atom index. maxsize=n_atom
       integer(kind=sp)                     dos_ldos_natom ! how many atoms to be plotted for ldos
       real(kind=dp),      allocatable   :: ldos_tot(:,:,:,:) !LDOS for each (ORB, ATOM, ENERGY, SPIN)

       character(len=1  )            dos_kunit

       logical                       dos_flag_sparse
  endtype dos

  type berry ! PINPT_BERRY
       integer(kind=sp)                     mysystem   ! my system index

       logical                       flag_wcc_evolve
       character(len=40)                  wcc_filenm, wcc_gap_filenm
       integer(kind=sp)                     wcc_nerange
       integer(kind=sp),   allocatable   :: wcc_erange(:)
       real(kind=dp),      allocatable   :: wcc(:,:,:) ! (nerange/nspin,nspin,nkpath)
       real(kind=dp),      allocatable   :: wcc_kpath(:,:,:) ! k-path (1:3, 2, nkline) (1:3, init-fina, nkline)
       real(kind=dp)                        wcc_kpath_shift(3)
       real(kind=dp),      allocatable   :: wcc_kpoint(:,:,:) ! k-points along the k-path (1:3,1:nk,nkline)
       real(kind=dp),      allocatable   :: wcc_kpoint_reci(:,:,:) ! k-points along the k-path (1:3,1:nk,nkline) (reci unit)
       integer(kind=sp)                     wcc_nkdiv, wcc_nkdiv2         ! nkdiv along k-path, nkdiv2 along k-path wcc_evolve
       integer(kind=sp)                     wcc_nkpath
       integer(kind=sp)                     wcc_direction
       character(len=256)                 strip_wcc_range
       logical                       flag_wcc_phase
       logical                       flag_wcc_get_chern
       logical                       flag_wcc_get_chern_spin
       real(kind=dp),      allocatable   :: wcc_polarization(:,:)   ! (nspin, nkpath)
       real(kind=dp),      allocatable   :: wcc_chern(:)           ! chern number (nspin)

       logical                       flag_zak_evolve
       character(len=40)                  zak_filenm 
       integer(kind=sp)                     zak_nerange
       integer(kind=sp),   allocatable   :: zak_erange(:)          ! default size = neig * ispin
       real(kind=dp),      allocatable   :: zak_phase(:,:) !(nspin,nkpath) = (erange)
       real(kind=dp),      allocatable   :: zak_kpath(:,:,:) ! k-path (1:3, 2, nkline) (1:3, init-fina, nkpath)
       real(kind=dp)                        zak_kpath_shift(3)
       real(kind=dp),      allocatable   :: zak_kpoint(:,:,:) ! k-points along the k-path (1:3,1:nk,nkpath)
       real(kind=dp),      allocatable   :: zak_kpoint_reci(:,:,:) ! k-points along the k-path (1:3,1:nk,nkpath) (reci unit)
       integer(kind=sp)                     zak_nkdiv, zak_nkdiv2  ! nkdiv along one k-path, nkdiv2 how many k-path (zak_evolve) along zak_dir
       integer(kind=sp)                     zak_nkpath
       integer(kind=sp)                     zak_direction
       character(len=256)                 strip_zak_range
       real(kind=dp),      allocatable   :: polarization(:)
       logical                       flag_zak_phase
 

!      logical                       
       character(len=256)                 strip_z2_range
       character(len=40)                  z2_filenm, z2_gap_filenm
       integer(kind=sp)                     z2_nkdiv, z2_nkpath
       integer(kind=sp)                     z2_nerange
       integer(kind=sp),   allocatable   :: z2_erange(:)          ! default size = neig * ispin
       integer(kind=sp)                     z2_dimension, z2_nplane
       integer(kind=sp),   allocatable   :: z2_axis(:)
       real(kind=dp),      allocatable   :: z2_kpoint(:,:,:,:,:)      ! kp along k-path (3,nkdiv,nkpath,nplane,naxis) (nxis =3(3D), 1(2D,1D))
       real(kind=dp),      allocatable   :: z2_kpoint_reci(:,:,:,:,:) ! kp along k-path (3,nkdiv,nkpath,nplane,naxis) (nplane = 2(3D), 1(2D,1D))
       real(kind=dp),      allocatable   :: z2_wcc(:,:,:,:,:)         ! wcc (nerange/nspin,nspin,nkpath,nplane,naxis)
       logical                       flag_z2_get_chern
       real(kind=dp),      allocatable   :: z2_chern(:,:,:)           ! chern number (nspin, nplane, naxis)
       real(kind=dp),      allocatable   :: z2_polarization(:,:,:,:)           ! chern number (nspin, nkpath, nplane, naxis)
       logical                       flag_z2_phase

       character(len=40)                  bc_filenm 
       logical                       flag_bc_filenm_provided
       logical                       flag_bc_method_kubo, flag_bc_method_resta
       integer(kind=sp)                     bc_dimension
       integer(kind=sp)                     bc_nkdiv(3)
       integer(kind=sp),   allocatable   :: bc_axis(:)
       integer(kind=sp)                     bc_nerange
       integer(kind=sp),   allocatable   :: bc_erange(:) ! default size = neig * ispin
       real(kind=dp),      allocatable   :: berrycurv(:,:,:) ! (:,:,:) = (1:3 ->Omega_x:z, nkpoints, erange)
       real(kind=dp),      allocatable   :: omega(:,:,:,:) ! (neig * ispinor, 3, nspin, nkpoint)
       character(len=256)                 strip_bc_range
       logical                       flag_bc_phase

       real(kind=dp)                        parity_origin(3) !direct coordinate of the orgiin of the system in the given unit cell
       real(kind=dp)                        parity_operator(3,3) 
       integer(kind=sp)                     parity_nkpoint
       real(kind=dp),      allocatable   :: parity_kpoint(:,:) !(3,parity_nkpoint)
       real(kind=dp),      allocatable   :: parity_kpoint_reci(:,:) !(3,parity_nkpoint)
       character(len=10),allocatable   :: parity_kpoint_name(:)
       integer(kind=sp)                     noccupied
       logical                       flag_print_hamiltonian_parity, flag_parity_phase

       real(kind=dp)                        symmetry_origin(3)
       real(kind=dp)                        symmetry_operator(3,3)
       integer(kind=sp)                     symmetry_nkpoint
       real(kind=dp),      allocatable   :: symmetry_kpoint(:,:)
       real(kind=dp),      allocatable   :: symmetry_kpoint_reci(:,:)
       character(len=10),allocatable   :: symmetry_kpoint_name(:)
       real(kind=dp)                        symmetry_theta ! rotation angle along z-axis (definite specification)
       logical                       flag_print_hamiltonian_symmetry, flag_symmetry_phase
 

  endtype berry

  type spmat ! sparse matrix with Compressed Sparse Row format
       integer(kind=sp)                     nnz     ! number of non-zero elements
       integer(kind=sp)                     msize   ! matrix size
       complex(kind=dp),  allocatable   :: H(:)    ! sparse array of square matrix H_square(msize,msize), (n_neighbor)
       integer(kind=sp),   allocatable   :: I(:)    ! Row    array : if CSR format : I(msize + 1), I(msize+1)-1 = nnz
                                             !                if COO format : I(nnz)
       integer(kind=sp),   allocatable   :: J(:)    ! Column array : J(nnz), same size with H
  endtype spmat

  type gainp ! input/control parameters for genetic algorithm
       integer(kind=sp)                     mysystem   ! my system index

       integer(kind=sp)                     mgen    ! maximum number of GA iterations (generation) (default : 500)
       integer(kind=sp)                     npop    ! number of individuals in a population (default is 100)
       integer(kind=sp)                     ngene   ! number of significant digits (number of genes) retained 
                                             !           in chromosomal encoding (default is 6).
       real(kind=dp)                        pcross  ! crossover probability ; must be  <= 1.0 (default : 0.85)
       real(kind=dp)                        pmutmn  ! minimum mutation rate; must be >= 0.0 (default is 0.0005)
       real(kind=dp)                        pmutmx  ! mmaximum mutation rate; must be <= 1.0 (default is 0.25)
       real(kind=dp)                        pmut    ! initial mutation rate; should be small (default is 0.005)
       integer(kind=sp)                     imut    ! mutation mode; 1/2/3/4/5 (default is 2).
                                             ! 1=one-point mutation, fixed rate.
                                             ! 2=one-point, adjustable rate based on fitness.
                                             ! 3=one-point, adjustable rate based on distance.
                                             ! 4=one-point+creep, fixed rate.
                                             ! 5=one-point+creep, adjustable rate based on fitness.
                                             ! 6=one-point+creep, adjustable rate based on distance.
       real(kind=dp)                        fdif    ! relative fitness differential; range from 0(none) to 1(maximum). (default 1)
       integer(kind=sp)                     irep    ! reproduction plan; 1 = Full generational replacement
                                             !                    2 = Steady-state-replace-random
                                             !                    3 = Steady-state-replace-worst (default)
       integer(kind=sp)                     ielite  ! elitism flag; 0/1=off/on (default is 0)
                                             !               (Applies only to reproduction plans 1 and 2)
       integer(kind=sp)                     ivrb    ! printed output 0/1/2=None/Minimal/Verbose  (default : 0)
       real(kind=dp)                        convtol ! convergence tolerance; must be > 0.0 (default is 0.0001)
       integer(kind=sp)                     convwin ! convergence window; must be >= 0
                                             ! This is the number of consecutive solutions
                                             ! within the tolerance for convergence to be declared (default is 20)
       real(kind=dp)                        iguessf ! fraction of the initial population to set equal to the initial guess.
                                             ! (none) to 1.0 (all). (default is 0.1 or 10%).
       integer(kind=sp)                     iseed   ! random seed value; must be > 0 (default is 999)
       real(kind=dp)                        lower_bound, upper_bound
  endtype gainp

  type replot ! PRPLT, required parameters for replot dos/ldos with given bandstructure file
       integer(kind=sp)                     mysystem   ! my system index

       logical                       flag_replot            ! whether turn on replot mode

       logical                       flag_replot_dos        ! flag for density of states
       logical                       flag_replot_ldos       ! flag for local density of states
       logical                       flag_replot_sldos      ! flag for spatial local density of states
       logical                       flag_replot_band       ! flag for replot band structure
       logical                       flag_replot_didv       ! flag for replot spatial density of states (LDOS) for given energy window
       logical                       flag_replot_proj_band  ! flag for replot projected band

       logical                       flag_replot_formatted ! flag whether band_structure_TBA file is formatted (.dat) or unformatted (.bin)

       character(len=80)                  fname_band_up, fname_band_dn  ! file name the be read (set to default)
       real(kind=dp)                        replot_dos_smearing ! gaussian smearing
       real(kind=dp),      allocatable   :: replot_dos_erange(:)   ! size=nediv, division of energy 
       integer(kind=sp)                     replot_dos_nediv     ! number of division
       real(kind=dp)                        replot_dos_emin, replot_dos_emax  ! from emin to emax
       real(kind=dp),      allocatable   :: replot_dos_tot(:,:)     ! density of states, (nspin, nediv)
       real(kind=dp),      allocatable   :: replot_dos_ntot(:,:) ! integrated DOS from initial to the energy level (nspin,0:nediv)

       integer(kind=sp)                     replot_nldos_sum  ! number of ldos plot with set of atoms (how many REPLOT_LDOS tag given?)
       integer(kind=sp),   allocatable   :: replot_ldos_natom(:) ! number of atoms to be resolved for each REPLOT_LDOS request (nldos_sum)
       integer(kind=sp),   allocatable   :: replot_ldos_atom(:,:) ! atom index
       real(kind=dp),      allocatable   :: replot_ldos_tot(:,:,:,:,:)   ! local DOS (n_orbital(iatom), maxval(ldos_natom), nspin, nediv,nldos_sum)

       character(len=40)                  replot_sldos_fname
       character(len=40)                  replot_dos_fname
       character(len=40)                  replot_didv_fname
       real(kind=dp),      allocatable   :: replot_sldos_sum(:,:,:) ! spatial local density of states sum over( natom, nspin, nediv )
       integer(kind=sp)                     replot_sldos_cell(3)  ! repeat cell along a1, a2, a3 direction
       real(kind=dp)                        r_origin(3) ! direct coordinate for shift of the origin of atomic coordinate of SLDOS
       real(kind=dp)                        bond_cut    ! bond length <= bond_cut will not be written in BOND.replot.dat Default: 3.0 (ang)


       integer(kind=sp)                     replot_nband ! number of REPLOT_BAND request
       character(len=2), allocatable   :: replot_axis_print_mag(:) ! (nband)
       logical,     allocatable   :: flag_replot_print_single(:) ! (nband)
       logical,     allocatable   :: flag_replot_write_unformatted(:) ! (nband)

       ! projected band
       integer(kind=sp)                     replot_nproj_sum
       integer(kind=sp),   allocatable   :: replot_proj_natom(:) ! how many atoms to be plotted for projected band
       integer(kind=sp),   allocatable   :: replot_proj_atom(:,:) ! integer array of atom index. maxsize=n_atom
     
  endtype replot
endmodule 
