#include "alias.inc"
module pyfit
    use parameters
    use mpi_setup
    use print_io
    use mykind
    use do_math
    use projected_band
    use,intrinsic :: iso_fortran_env
    implicit none

    type incar_py  ! copied from parameters.f90: type incar
       ! set in parsing step
       logical                            flag_python_module 
       logical                            flag_get_band  
       logical                            flag_fit_degeneracy 
       logical                            flag_report_geom
       logical                            flag_phase  
       logical                            flag_tbfit
       logical                            flag_tbfit_finish 
       logical                            flag_print_only_target
       logical                            flag_print_energy_diff 
       logical                            flag_print_orbital
       logical                            flag_get_orbital  
       logical                            flag_print_mag
       logical                            flag_print_single 
       logical                            flag_local_charge
       logical                            flag_collinear
       logical                            flag_noncollinear
       logical                            flag_soc
       logical                            flag_plus_U
       logical                            flag_scissor
       logical                            flag_print_proj
       logical                            flag_print_proj_sum
       logical                            flag_print_proj_atom
       logical                            flag_plot_fit
       logical                            flag_plot
       logical                            flag_get_band_order   
       logical                            flag_get_band_order_print_only 
       logical                            flag_use_weight 
       logical                            flag_fit_orbital
       logical                            flag_fit_orbital_parse
       logical                            flag_pso_with_lmdif
       logical                            flag_get_unfold     
       real(kind=dp)                      ptol
       real(kind=dp)                      ftol
       real(kind=dp)                      fdiff
       real(kind=dp)                      r_scissor
       integer(kind=sp)                   nsystem        
       integer(kind=sp)                   lmmax 
       integer(kind=sp)                   miter
       integer(kind=sp)                   mxfit
       integer(kind=sp)                   nn_max(3) 
       integer(kind=sp)                   ispin   
       integer(kind=sp)                   ispinor 
       integer(kind=sp)                   nspin   
       integer(kind=sp)                   nproj_sum
       integer(kind=sp)                   iverbose
       integer(kind=sp)                   i_scissor
       character(len=10)                  ls_type 
       character(len=40)                  fnamelog       

       integer(kind=sp),   allocatable :: proj_atom(:,:) 
       integer(kind=sp),   allocatable :: proj_natom(:) 
       character(len=132), allocatable :: ifilenm(:)     
       character(len=132), allocatable :: title(:)       
    
       logical                            flag_get_total_energy
       real(kind=dp)                      electronic_temperature

       real(kind=dp)                      orbital_fit_smearing
    endtype incar_py

  type params_py !PPRAM_PY
       logical                            flag_set_param_const     
       logical                            flag_pfile_index         
       logical                            flag_use_overlap         
       logical                            flag_slater_koster       
       logical                            flag_nrl_slater_koster   
       logical                            flag_fit_plain
       real(kind=dp)                      l_broaden                
       real(kind=dp)                      pso_c1
       real(kind=dp)                      pso_c2
       real(kind=dp)                      pso_w 
       real(kind=dp)                      pso_max_noise_amplitude
       integer(kind=sp)                   pso_miter
       integer(kind=sp)                   slater_koster_type       
       integer(kind=sp)                   nparam                   
       integer(kind=sp)                   nparam_const             
       integer(kind=sp)                   nparam_nrl               
       integer(kind=sp)                   nparam_nrl_free          
       integer(kind=sp)                   param_nsub_max           
       integer(kind=sp)                   nparam_free              
       integer(kind=sp)                   pso_nparticles
       character(len=132)                 pfilenm
       character(len=132)                 pfileoutnm       

       real(kind=dp),     allocatable  :: param(:)
       real(kind=dp),     allocatable  :: param_best(:)
       real(kind=dp),     allocatable  :: param_nrl(:,:)           
       real(kind=dp),     allocatable  :: param_const(:,:)         
!      real(kind=dp),     allocatable  :: param_const_best(:,:)         
       real(kind=dp),     allocatable  :: param_const_nrl(:,:,:)   
       integer(kind=sp),  allocatable  :: param_nsub(:)            
       integer(kind=sp),  allocatable  :: iparam_free(:)           
       integer(kind=sp),  allocatable  :: iparam_free_nrl(:)       
       character(len=40), allocatable  :: param_name(:)
       character(len=40), allocatable  :: c_const(:,:)

       real(kind=dp),    allocatable   :: pso_cost_history(:)
       real(kind=dp),    allocatable   :: pso_cost_history_i(:,:)
       real(kind=dp),    allocatable   :: pso_pbest_history(:,:,:)
       real(kind=dp),    allocatable   ::     cost_history(:)
       integer(kind=sp)                   niter
  endtype params_py
    
  type kpoints_py !PKPTS_PY
       logical                            flag_klinemode
       logical                            flag_kgridmode
       logical                            flag_gamma
       logical                            flag_reciprocal
       logical                            flag_cartesianK
       real(kind=dp)                      k_shift(3)
       integer(kind=sp)                   mysystem   
       integer(kind=sp)                   nkpoint,nline
       integer(kind=sp)                   n_ndiv
       integer(kind=sp)                   idiv_mode 
       integer(kind=sp)                   kreduce 
       character(len=132)                 kfilenm            
       character(len=132)                 ribbon_kfilenm     
       character(len=8)                   kline_type 
       character(len=1)                   kunit

       real(kind=dp),      allocatable :: kpoint(:,:)
       real(kind=dp),      allocatable :: kline(:,:)
       real(kind=dp),      allocatable :: kpoint_reci(:,:)
       integer(kind=sp),   allocatable :: ndiv(:)
       character(len=8),   allocatable :: k_name(:)
       character(len=8),   allocatable :: k_name2(:)
       integer(kind=sp),   allocatable :: k_name_index(:)

       real(kind=dp),      allocatable :: kdist(:) ! only for python module
  endtype kpoints_py

  type weight_py 
       logical                            flag_weight_default
       logical                            flag_weight_orb
       real(kind=dp)                      efile_ef           
       integer(kind=sp)                   mysystem   
       integer(kind=sp)                   itarget_e_start    
       integer(kind=sp)                   read_energy_column_index
       integer(kind=sp)                   read_energy_column_index_dn 
       integer(kind=sp)                   ie_cutoff    
       character(len=132)                 efilenmu
       character(len=132)                 efilenmd  
       character(len=16)                  efile_type         

       real(kind=dp),      allocatable :: WT(:,:)
  endtype weight_py

  type poscar_py 
       logical                            flag_selective
       logical                            flag_direct
       logical                            flag_cartesian
       integer(kind=sp)                   mysystem   
       integer(kind=sp)                   neig       
       integer(kind=sp)                   neig_total 
       integer(kind=sp)                   neig_target
       integer(kind=sp)                   nbasis  
       integer(kind=sp)                   neig_eff   
       integer(kind=sp)                   init_erange
       integer(kind=sp)                   fina_erange  
       integer(kind=sp)                   nband 
       integer(kind=sp)                   n_spec
       integer(kind=sp)                   n_atom
       integer(kind=sp)                   max_orb  
       character(len=132)                 title      
       character(len=132)                 gfilenm    
       character(len=40)                  system_name
       character(len=8), allocatable   :: c_orbital(:,:)
       character(len=8), allocatable   :: c_spec(:) ! character of species for each species (1:n_spec)
       real(kind=dp)                      a_scale
       real(kind=dp)                      a_latt(3,3) 
       real(kind=dp)                      b_latt(3,3) 

       real(kind=dp),    allocatable   :: nelect(:)  
       integer(kind=sp), allocatable   :: n_orbital(:) 
       integer(kind=sp), allocatable   :: orb_index(:) 
       integer(kind=sp), allocatable   :: spec(:) ! species information for each atom (1:n_atom). The order of appearance in the POSCAR
  endtype poscar_py

  type hopping_py !NN_TABLE_PY (nearest neighbor table; but not restricted to nearest)
       logical                            flag_efield
       logical                            flag_efield_frac
       logical                            flag_efield_cart
       real(kind=dp)                      onsite_tolerance
       real(kind=dp)                      efield(3)
       real(kind=dp)                      efield_origin(3)
       real(kind=dp)                      efield_origin_cart(3)
       integer(kind=sp)                   mysystem   
       integer(kind=sp)                   n_neighbor 
       integer(kind=sp)                   max_nn_pair  
       character(len=132)                 nnfilenm      
       character(len=132)                 nnfilenmo     

       logical,          allocatable   :: flag_site_cindex(:) 
       real(kind=dp),    allocatable   :: i_coord(:,:)  
       real(kind=dp),    allocatable   :: j_coord(:,:)  
       real(kind=dp),    allocatable   :: Rij(:,:) 
       real(kind=dp),    allocatable   :: R(:,:)   
       real(kind=dp),    allocatable   :: Dij(:)   
       real(kind=dp),    allocatable   :: Dij0(:)
       real(kind=dp),    allocatable   :: i_sign(:)   
       real(kind=dp),    allocatable   :: j_sign(:)   
       real(kind=dp),    allocatable   :: R_nn(:,:)  
       real(kind=dp),    allocatable   :: R0_nn(:,:) 
       real(kind=dp),    allocatable   :: local_charge(:) 
       real(kind=dp),    allocatable   :: local_moment(:,:) 
       real(kind=dp),    allocatable   :: local_moment_rot(:,:) 
       integer(kind=sp), allocatable   :: i_atom(:)     
       integer(kind=sp), allocatable   :: j_atom(:)
       integer(kind=sp), allocatable   :: i_matrix(:) 
       integer(kind=sp), allocatable   :: j_matrix(:) 
       integer(kind=sp), allocatable   :: n_class(:)
       integer(kind=sp), allocatable   :: sk_index_set(:,:)  ! (0:6 + ovrl, nn)
       integer(kind=sp), allocatable   :: cc_index_set(:,:)  ! (0:3, nn)
       integer(kind=sp), allocatable   :: n_nn(:)
       integer(kind=sp), allocatable   :: j_nn(:,:) 
       integer(kind=sp), allocatable   :: l_onsite_param_index(:) 
       integer(kind=sp), allocatable   :: stoner_I_param_index(:) 
       integer(kind=sp), allocatable   :: local_U_param_index(:) 
       integer(kind=sp), allocatable   :: plus_U_param_index(:) 
       integer(kind=sp), allocatable   :: soc_param_index(:)    
       character(len=8), allocatable   :: ci_orb(:)   
       character(len=8), allocatable   :: cj_orb(:)   
       character(len=2), allocatable   :: p_class(:)
       character(len=20),allocatable   :: site_cindex(:)  
  endtype hopping_py

  type energy_py
       integer(kind=sp)                   mysystem

       real(kind=dp),    allocatable   :: E(:,:)
       real(kind=dp),    allocatable   :: dE(:,:)
       real(kind=dp),    allocatable   :: dORB(:,:)
       real(kind=dp),    allocatable   :: ORB(:,:,:)
       complex(kind=dp), allocatable   :: V(:,:,:)
       complex(kind=dp), allocatable   :: SV(:,:,:)
       real(kind=dp),    allocatable   :: V2(:,:,:) ! projected orbital |psi_i><psi_i|
       real(kind=dp),    allocatable   :: E_BAND(:)
       real(kind=dp),    allocatable   :: E_TOT(:)
       real(kind=dp),    allocatable   :: F_OCC(:,:)
       real(kind=dp)                   :: E_F       
       
  endtype 

contains

subroutine init3(comm, PINPT_PY, PPRAM_PY, PKPTS_PY1, PWGHT_PY1, PGEOM_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                                           PKPTS_PY2, PWGHT_PY2, PGEOM_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                                           PKPTS_PY3, PWGHT_PY3, PGEOM_PY3, NN_TABLE_PY3, EDFT_PY3, ETBA_PY3)
    use read_incar
    use set_default
    type(incar_py),       intent(inout) :: PINPT_PY
    type(params_py),      intent(inout) :: PPRAM_PY
    type(kpoints_py),     intent(inout) :: PKPTS_PY1   , PKPTS_PY2    , PKPTS_PY3
    type(weight_py),      intent(inout) :: PWGHT_PY1   , PWGHT_PY2    , PWGHT_PY3
    type(poscar_py),      intent(inout) :: PGEOM_PY1   , PGEOM_PY2    , PGEOM_PY3
    type(hopping_py),     intent(inout) :: NN_TABLE_PY1, NN_TABLE_PY2 , NN_TABLE_PY3
    type(energy_py),      intent(inout) :: EDFT_PY1    , EDFT_PY2     , EDFT_PY3
    type(energy_py),      intent(inout) :: ETBA_PY1    , ETBA_PY2     , ETBA_PY3
    integer,              intent(in)    :: comm
    integer(kind=sp)                       i, j, k
    type(incar   )                      :: PINPT, PINPT_
    type(params  )                      :: PPRAM, PPRAM1, PPRAM_
    type(kpoints ),dimension(3)         :: PKPTS        ! assume nsystem = 3
    type(weight  ),dimension(3)         :: PWGHT        ! assume nsystem = 3
    type(poscar  ),dimension(3)         :: PGEOM        ! assume nsystem = 3
    type(hopping ),dimension(3)         :: NN_TABLE     ! assume nsystem = 3
    type(energy  ),dimension(3)         :: EDFT, ETBA   ! assume nsystem = 3
    type(dos     )                      :: TMP1 ! temp
    type(berry   )                      :: TMP2 ! temp
    type(gainp   )                      :: TMP3 ! temp
    type(replot  )                      :: TMP4 ! temp
    type(unfold  )                      :: PUFLD ! unfold
    integer(kind=sp)                       nsystem
    real(kind=dp), allocatable          :: param_best(:)
#ifdef MPI
   ! initialize MPI setup
    mpi_comm_earth = comm
    call get_my_task()
#endif

   ! init system
    call parse_very_init_py(PINPT, PINPT_PY%nsystem, PINPT_PY%ifilenm(1))
    call parse(PINPT)
    PINPT%flag_fit_orbital_parse = PINPT_PY%flag_fit_orbital_parse
    PINPT%flag_fit_orbital       = PINPT%flag_fit_orbital_parse
    nsystem = PINPT%nsystem

    call set_verbose(PINPT_PY%iverbose)
    if(allocated(PPRAM_PY%param_best)) then
        param_best = PPRAM_PY%param_best
    endif
    call read_input(PINPT,PPRAM_,PKPTS(1),PGEOM(1),PWGHT(1),EDFT(1),NN_TABLE(1), TMP1,TMP2,TMP3,TMP4, PUFLD, 1)
    PINPT_ = PINPT
    PPRAM1 = PPRAM_
    PINPT  = PINPT_
    PPRAM  = PPRAM1
    do i = 2,nsystem
      call read_input(PINPT,PPRAM_,PKPTS(i),PGEOM(i),PWGHT(i),EDFT(i),NN_TABLE(i), TMP1,TMP2,TMP3,TMP4, PUFLD, i)
      PINPT_%title(i) = PINPT%title(i)
      PINPT_%ifilenm(i)= PINPT%ifilenm(i)
      PINPT = PINPT_
      PPRAM = PPRAM1
    enddo

    call set_free_parameters(PPRAM)

    do i = 1, nsystem
      call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i), .FALSE., 0)
      if(allocated(ETBA(i)%dE))   deallocate(ETBA(i)%dE)
      if(allocated(ETBA(i)%dORB)) deallocate(ETBA(i)%dORB)
      allocate(ETBA(i)%dE(   size(ETBA(i)%E(:,1)),size(ETBA(i)%E(1,:)) ))
      allocate(ETBA(i)%dORB( size(ETBA(i)%E(:,1)),size(ETBA(i)%E(1,:)) ))

      ETBA(i)%E = 0d0 ; ETBA(i)%V = 0d0 ; ETBA(i)%SV = 0d0 ; ETBA(i)%dE=0d0 ; ETBA(i)%dORB = 0d0
      if(PINPT%flag_fit_orbital) then
        ETBA(i)%ORB  = 0d0
        ETBA(i)%V2   = 0d0
      endif

    enddo

    if(PKPTS(1)%flag_klinemode) then
      if(allocated(PKPTS_PY1%kdist)) deallocate(PKPTS_PY1%kdist)
      allocate(PKPTS_PY1%kdist(PKPTS(1)%nkpoint))
      call get_kline_dist(PKPTS(1)%kpoint, PKPTS(1)%nkpoint, PKPTS_PY1%kdist) ! assume that tbfit uses line mode for kpoints
    endif

    if(PKPTS(2)%flag_klinemode) then
      if(allocated(PKPTS_PY2%kdist)) deallocate(PKPTS_PY2%kdist)
      allocate(PKPTS_PY2%kdist(PKPTS(2)%nkpoint))
      call get_kline_dist(PKPTS(2)%kpoint, PKPTS(2)%nkpoint, PKPTS_PY2%kdist) ! assume that tbfit uses line mode for kpoints
    endif

    if(PKPTS(3)%flag_klinemode) then
      if(allocated(PKPTS_PY3%kdist)) deallocate(PKPTS_PY3%kdist)
      allocate(PKPTS_PY3%kdist(PKPTS(3)%nkpoint))
      call get_kline_dist(PKPTS(3)%kpoint, PKPTS(3)%nkpoint, PKPTS_PY3%kdist) ! assume that tbfit uses line mode for kpoints
    endif


    call copy_incar(PINPT_PY,      PINPT, 1)
    call copy_params(PPRAM_PY,     PPRAM, 1)

    call copy_all(PKPTS_PY1, PGEOM_PY1, PWGHT_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                  PKPTS(1) , PGEOM(1) , PWGHT(1) , NN_TABLE(1) , EDFT(1) , ETBA(1) , 1)
    call copy_all(PKPTS_PY2, PGEOM_PY2, PWGHT_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                  PKPTS(2) , PGEOM(2) , PWGHT(2) , NN_TABLE(2) , EDFT(2) , ETBA(2) , 1)
    call copy_all(PKPTS_PY3, PGEOM_PY3, PWGHT_PY3, NN_TABLE_PY3, EDFT_PY3, ETBA_PY3, &
                  PKPTS(3) , PGEOM(3) , PWGHT(3) , NN_TABLE(3) , EDFT(3) , ETBA(3) , 1)

    if(allocated(param_best)) then
        PPRAM_PY%param_best = param_best
    endif
!   write(6,*)"ZZZZ ", PINPT%flag_fit_orbital, PINPT_PY%flag_fit_orbital
!stop

    return
endsubroutine

subroutine init2(comm, PINPT_PY, PPRAM_PY, PKPTS_PY1, PWGHT_PY1, PGEOM_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                                           PKPTS_PY2, PWGHT_PY2, PGEOM_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2)
    use read_incar
    use set_default
    type(incar_py),       intent(inout) :: PINPT_PY
    type(params_py),      intent(inout) :: PPRAM_PY
    type(kpoints_py),     intent(inout) :: PKPTS_PY1   , PKPTS_PY2
    type(weight_py),      intent(inout) :: PWGHT_PY1   , PWGHT_PY2
    type(poscar_py),      intent(inout) :: PGEOM_PY1   , PGEOM_PY2
    type(hopping_py),     intent(inout) :: NN_TABLE_PY1, NN_TABLE_PY2
    type(energy_py),      intent(inout) :: EDFT_PY1    , EDFT_PY2
    type(energy_py),      intent(inout) :: ETBA_PY1    , ETBA_PY2
    integer,              intent(in)    :: comm
    integer(kind=sp)                       i, j, k
    type(incar   )                      :: PINPT, PINPT_
    type(params  )                      :: PPRAM, PPRAM1, PPRAM_
    type(kpoints ),dimension(2)         :: PKPTS ! assume nsystem = 2
    type(weight  ),dimension(2)         :: PWGHT ! assume nsystem = 2
    type(poscar  ),dimension(2)         :: PGEOM ! assume nsystem = 2
    type(hopping ),dimension(2)         :: NN_TABLE ! assume nsystem = 2
    type(energy  ),dimension(2)         :: EDFT, ETBA ! assume nsystem = 2
    type(dos     )                      :: TMP1 ! temp
    type(berry   )                      :: TMP2 ! temp
    type(gainp   )                      :: TMP3 ! temp
    type(replot  )                      :: TMP4 ! temp
    type(unfold  )                      :: PUFLD ! unfold
    real(kind=dp), allocatable          :: param_best(:)
#ifdef MPI
   ! initialize MPI setup
    mpi_comm_earth = comm
    call get_my_task()
#endif

   ! init system
    call parse_very_init_py(PINPT, PINPT_PY%nsystem, PINPT_PY%ifilenm(1))
    call parse(PINPT)
    PINPT%flag_fit_orbital_parse = PINPT_PY%flag_fit_orbital_parse
    PINPT%flag_fit_orbital       = PINPT%flag_fit_orbital_parse

    call set_verbose(PINPT_PY%iverbose)
    if(allocated(PPRAM_PY%param_best)) then
        param_best = PPRAM_PY%param_best
    endif
    call read_input(PINPT,PPRAM_,PKPTS(1),PGEOM(1),PWGHT(1),EDFT(1),NN_TABLE(1), TMP1,TMP2,TMP3,TMP4, PUFLD, 1)
    PINPT_ = PINPT
    PPRAM1 = PPRAM_
    PINPT  = PINPT_
    PPRAM  = PPRAM1
    call read_input(PINPT,PPRAM_,PKPTS(2),PGEOM(2),PWGHT(2),EDFT(2),NN_TABLE(2), TMP1,TMP2,TMP3,TMP4, PUFLD, 2)
    PINPT_%title(2) = PINPT%title(2)
    PINPT_%ifilenm(2)= PINPT%ifilenm(2)
    PINPT = PINPT_
    PPRAM = PPRAM1

    call set_free_parameters(PPRAM)

    do i = 1, 2
      call allocate_ETBA(PGEOM(i), PINPT, PKPTS(i), ETBA(i), .FALSE., 0)
      if(allocated(ETBA(i)%dE))   deallocate(ETBA(i)%dE)
      if(allocated(ETBA(i)%dORB)) deallocate(ETBA(i)%dORB)
      allocate(ETBA(i)%dE(   size(ETBA(i)%E(:,1)),size(ETBA(i)%E(1,:)) ))
      allocate(ETBA(i)%dORB( size(ETBA(i)%E(:,1)),size(ETBA(i)%E(1,:)) ))

      ETBA(i)%E = 0d0 ; ETBA(i)%V = 0d0 ; ETBA(i)%SV = 0d0 ; ETBA(i)%dE=0d0 ; ETBA(i)%dORB = 0d0
      if(PINPT%flag_fit_orbital) then
        ETBA(i)%ORB  = 0d0
        ETBA(i)%V2   = 0d0
      endif

    enddo
    if(PKPTS(1)%flag_klinemode) then
      if(allocated(PKPTS_PY1%kdist)) deallocate(PKPTS_PY1%kdist)
      allocate(PKPTS_PY1%kdist(PKPTS(1)%nkpoint))
      call get_kline_dist(PKPTS(1)%kpoint, PKPTS(1)%nkpoint, PKPTS_PY1%kdist) ! assume that tbfit uses line mode for kpoints
    endif

    if(PKPTS(2)%flag_klinemode) then
      if(allocated(PKPTS_PY2%kdist)) deallocate(PKPTS_PY2%kdist)
      allocate(PKPTS_PY2%kdist(PKPTS(2)%nkpoint))
      call get_kline_dist(PKPTS(2)%kpoint, PKPTS(2)%nkpoint, PKPTS_PY2%kdist) ! assume that tbfit uses line mode for kpoints
    endif

    call copy_incar(PINPT_PY,      PINPT, 1)
    call copy_params(PPRAM_PY,     PPRAM, 1)

    call copy_all(PKPTS_PY1, PGEOM_PY1, PWGHT_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                  PKPTS(1) , PGEOM(1) , PWGHT(1) , NN_TABLE(1) , EDFT(1) , ETBA(1) , 1)

    call copy_all(PKPTS_PY2, PGEOM_PY2, PWGHT_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                  PKPTS(2) , PGEOM(2) , PWGHT(2) , NN_TABLE(2) , EDFT(2) , ETBA(2) , 1)

    if(allocated(param_best)) then
        PPRAM_PY%param_best = param_best
    endif

    return
endsubroutine
subroutine init(comm, PINPT_PY, PPRAM_PY, PKPTS_PY, PWGHT_PY, PGEOM_PY, NN_TABLE_PY, EDFT_PY, ETBA_PY)
    use read_incar
    use set_default
    type(incar_py),       intent(inout) :: PINPT_PY
    type(params_py),      intent(inout) :: PPRAM_PY
    type(kpoints_py),     intent(inout) :: PKPTS_PY
    type(weight_py),      intent(inout) :: PWGHT_PY
    type(poscar_py),      intent(inout) :: PGEOM_PY
    type(hopping_py),     intent(inout) :: NN_TABLE_PY
    type(energy_py),      intent(inout) :: EDFT_PY
    type(energy_py),      intent(inout) :: ETBA_PY
    integer,              intent(in)    :: comm
    integer(kind=sp)                       i, j, k
    type(incar   )                      :: PINPT
    type(params  )                      :: PPRAM
    type(kpoints )                      :: PKPTS ! assume nsystem = 1
    type(weight  )                      :: PWGHT ! assume nsystem = 1
    type(poscar  )                      :: PGEOM ! assume nsystem = 1
    type(hopping )                      :: NN_TABLE ! assume nsystem = 1
    type(energy  )                      :: EDFT, ETBA ! assume nsystem = 1
    type(dos     )                      :: TMP1 ! temp
    type(berry   )                      :: TMP2 ! temp
    type(gainp   )                      :: TMP3 ! temp
    type(replot  )                      :: TMP4 ! temp
    type(unfold  )                      :: PUFLD ! unfold
    real(kind=dp), allocatable          :: param_best(:)
#ifdef MPI
   ! initialize MPI setup
    mpi_comm_earth = comm
    call get_my_task()
#endif

   ! init system
    call parse_very_init_py(PINPT, PINPT_PY%nsystem, PINPT_PY%ifilenm(1))
    call parse(PINPT)

    PINPT%flag_fit_orbital_parse = PINPT_PY%flag_fit_orbital_parse
    PINPT%flag_fit_orbital       = PINPT%flag_fit_orbital_parse

    call set_verbose(PINPT_PY%iverbose)
    if(allocated(PPRAM_PY%param_best)) then
        param_best = PPRAM_PY%param_best
    endif
    call read_input(PINPT,PPRAM,PKPTS,PGEOM,PWGHT,EDFT,NN_TABLE, TMP1,TMP2,TMP3,TMP4, PUFLD, 1)
    call set_free_parameters(PPRAM)

    call allocate_ETBA(PGEOM, PINPT, PKPTS, ETBA, .FALSE., 0)
    if(allocated(ETBA%dE))   deallocate(ETBA%dE)
    if(allocated(ETBA%dORB)) deallocate(ETBA%dORB)
    allocate(ETBA%dE(   size(ETBA%E(:,1)),size(ETBA%E(1,:)) )) 
    allocate(ETBA%dORB( size(ETBA%E(:,1)),size(ETBA%E(1,:)) ))

    ETBA%E = 0d0 ; ETBA%V = 0d0 ; ETBA%SV = 0d0 ; ETBA%dE=0d0 ; ETBA%dORB = 0d0
    if(PINPT%flag_fit_orbital) then 
      ETBA%ORB  = 0d0
      ETBA%V2   = 0d0
    endif

    if(PKPTS%flag_klinemode) then
      if(allocated(PKPTS_PY%kdist)) deallocate(PKPTS_PY%kdist)
      allocate(PKPTS_PY%kdist(PKPTS%nkpoint))
      call get_kline_dist(PKPTS%kpoint, PKPTS%nkpoint, PKPTS_PY%kdist) ! assume that tbfit uses line mode for kpoints
    endif

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_kpoints(PKPTS_PY, PKPTS, 1)
    call copy_poscar(PGEOM_PY, PGEOM, 1)
    call copy_weight(PWGHT_PY, PWGHT, 1)
    call copy_hopping(NN_TABLE_PY, NN_TABLE, 1)
    call copy_energy(EDFT_PY, EDFT, 1)
    call copy_energy(ETBA_PY, ETBA, 1)

    if(allocated(param_best)) then
        PPRAM_PY%param_best = param_best
    endif
    return
end subroutine init

subroutine pso3(comm, PINPT_PY, PPRAM_PY, PKPTS_PY1, PWGHT_PY1, PGEOM_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                                          PKPTS_PY2, PWGHT_PY2, PGEOM_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                                          PKPTS_PY3, PWGHT_PY3, PGEOM_PY3, NN_TABLE_PY3, EDFT_PY3, ETBA_PY3, iseed, pso_miter)
    use cost_function, only: get_dE12
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY1   ,PKPTS_PY2    ,PKPTS_PY3
    type(poscar_py),      intent(inout)         :: PGEOM_PY1   ,PGEOM_PY2    ,PGEOM_PY3
    type(weight_py),      intent(inout)         :: PWGHT_PY1   ,PWGHT_PY2    ,PWGHT_PY3
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY1,NN_TABLE_PY2 ,NN_TABLE_PY3
    type(energy_py),      intent(inout)         :: ETBA_PY1    ,ETBA_PY2     ,ETBA_PY3
    type(energy_py),      intent(inout)         :: EDFT_PY1    ,EDFT_PY2     ,EDFT_PY3
    integer,              intent(inout)         :: comm
    integer(kind=sp),     intent(inout)         :: iseed, pso_miter
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(kpoints ), dimension(PINPT_PY%nsystem) :: PKPTS
    type(weight  ), dimension(PINPT_PY%nsystem) :: PWGHT
    type(poscar  ), dimension(PINPT_PY%nsystem) :: PGEOM
    type(hopping ), dimension(PINPT_PY%nsystem) :: NN_TABLE
    type(energy  ), dimension(PINPT_PY%nsystem) :: EDFT, ETBA
    integer(kind=sp)                               mpierr
    logical                                        flag_with_lmdif

#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif

    call set_verbose(PINPT_PY%iverbose)

    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_all(PKPTS_PY1, PGEOM_PY1, PWGHT_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                    PKPTS(1) , PGEOM(1) , PWGHT(1) , NN_TABLE(1) , EDFT(1) , ETBA(1) , 2)
    call copy_all(PKPTS_PY2, PGEOM_PY2, PWGHT_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                    PKPTS(2) , PGEOM(2) , PWGHT(2) , NN_TABLE(2) , EDFT(2) , ETBA(2) , 2)
    call copy_all(PKPTS_PY3, PGEOM_PY3, PWGHT_PY3, NN_TABLE_PY3, EDFT_PY3, ETBA_PY3, &
                    PKPTS(3) , PGEOM(3) , PWGHT(3) , NN_TABLE(3) , EDFT(3) , ETBA(3) , 2)

    call pso_fit(PINPT, PPRAM, PKPTS, PWGHT, PGEOM, NN_tABLE, EDFT, iseed, pso_miter)
    call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, PINPT%flag_fit_orbital)

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_energy(ETBA_PY1, ETBA(1), 1)
    call copy_energy(ETBA_PY2, ETBA(2), 1)
    call copy_energy(ETBA_PY3, ETBA(3), 1)

    return
endsubroutine

subroutine pso2(comm, PINPT_PY, PPRAM_PY, PKPTS_PY1, PWGHT_PY1, PGEOM_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                                          PKPTS_PY2, PWGHT_PY2, PGEOM_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, iseed, pso_miter)
    use cost_function, only: get_dE12
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY1   ,PKPTS_PY2
    type(poscar_py),      intent(inout)         :: PGEOM_PY1   ,PGEOM_PY2
    type(weight_py),      intent(inout)         :: PWGHT_PY1   ,PWGHT_PY2
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY1,NN_TABLE_PY2
    type(energy_py),      intent(inout)         :: ETBA_PY1    ,ETBA_PY2
    type(energy_py),      intent(inout)         :: EDFT_PY1    ,EDFT_PY2
    integer,              intent(inout)         :: comm
    integer(kind=sp),     intent(inout)         :: iseed, pso_miter
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(kpoints ), dimension(PINPT_PY%nsystem) :: PKPTS
    type(weight  ), dimension(PINPT_PY%nsystem) :: PWGHT
    type(poscar  ), dimension(PINPT_PY%nsystem) :: PGEOM
    type(hopping ), dimension(PINPT_PY%nsystem) :: NN_TABLE
    type(energy  ), dimension(PINPT_PY%nsystem) :: EDFT, ETBA
    integer(kind=sp)                               mpierr
    logical                                        flag_with_lmdif

#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif

    call set_verbose(PINPT_PY%iverbose)

    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_all(PKPTS_PY1, PGEOM_PY1, PWGHT_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                    PKPTS(1) , PGEOM(1) , PWGHT(1) , NN_TABLE(1) , EDFT(1) , ETBA(1) , 2)
    call copy_all(PKPTS_PY2, PGEOM_PY2, PWGHT_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                    PKPTS(2) , PGEOM(2) , PWGHT(2) , NN_TABLE(2) , EDFT(2) , ETBA(2) , 2)

    call pso_fit(PINPT, PPRAM, PKPTS, PWGHT, PGEOM, NN_tABLE, EDFT, iseed, pso_miter)
    call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, PINPT%flag_fit_orbital)

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_energy(ETBA_PY1, ETBA(1), 1)
    call copy_energy(ETBA_PY2, ETBA(2), 1)

    return
endsubroutine

subroutine pso(comm, PINPT_PY, PPRAM_PY, PKPTS_PY, PWGHT_PY, PGEOM_PY, NN_TABLE_PY, EDFT_PY, ETBA_PY, iseed, pso_miter)
    use cost_function, only: get_dE12
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY
    type(poscar_py),      intent(inout)         :: PGEOM_PY
    type(weight_py),      intent(inout)         :: PWGHT_PY
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY
    type(energy_py),      intent(inout)         :: ETBA_PY
    type(energy_py),      intent(inout)         :: EDFT_PY
    integer,              intent(inout)         :: comm
    integer(kind=sp),     intent(inout)         :: iseed, pso_miter
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(kpoints ), dimension(PINPT_PY%nsystem) :: PKPTS
    type(weight  ), dimension(PINPT_PY%nsystem) :: PWGHT
    type(poscar  ), dimension(PINPT_PY%nsystem) :: PGEOM
    type(hopping ), dimension(PINPT_PY%nsystem) :: NN_TABLE
    type(energy  ), dimension(PINPT_PY%nsystem) :: EDFT, ETBA
    integer(kind=sp)                               mpierr
    logical                                        flag_with_lmdif
     
#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif

    call set_verbose(PINPT_PY%iverbose)

    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_all(PKPTS_PY , PGEOM_PY , PWGHT_PY , NN_TABLE_PY , EDFT_PY , ETBA_PY , &
                    PKPTS(1) , PGEOM(1) , PWGHT(1) , NN_TABLE(1) , EDFT(1) , ETBA(1) , 2)
    
    call pso_fit(PINPT, PPRAM, PKPTS, PWGHT, PGEOM, NN_tABLE, EDFT, iseed, pso_miter)
    call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, PINPT%flag_fit_orbital)

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_energy(ETBA_PY, ETBA(1), 1)

    return
endsubroutine

subroutine eig3(comm, PINPT_PY, PPRAM_PY, PKPTS_PY, PGEOM_PY, NN_TABLE_PY, ETBA_PY, &
                                          PKPTS_PY2, PGEOM_PY2, NN_TABLE_PY2, ETBA_PY2, &
                                          PKPTS_PY3, PGEOM_PY3, NN_TABLE_PY3, ETBA_PY3)
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY   ,PKPTS_PY2   ,PKPTS_PY3
    type(poscar_py),      intent(inout)         :: PGEOM_PY   ,PGEOM_PY2   ,PGEOM_PY3
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY,NN_TABLE_PY2,NN_TABLE_PY3
    type(energy_py),      intent(inout)         :: ETBA_PY    ,ETBA_PY2    ,ETBA_PY3
    integer,              intent(inout)         :: comm
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(hopping ), dimension(PINPT_PY%nsystem) :: NN_TABLE
    type(energy  ), dimension(PINPT_PY%nsystem) :: ETBA
    type(kpoints ), dimension(PINPT_PY%nsystem) :: PKPTS
    type(poscar  ), dimension(PINPT_PY%nsystem) :: PGEOM
    logical                                        flag_stat

#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif
    call set_verbose(PINPT_PY%iverbose)
    if(iverbose .eq. 1) then
        flag_stat = .TRUE.
    else
        flag_stat = .FALSE.
    endif

    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_hopping(NN_TABLE_PY , NN_TABLE(1), 2)
    call copy_hopping(NN_TABLE_PY2, NN_TABLE(2), 2)
    call copy_hopping(NN_TABLE_PY3, NN_TABLE(3), 2)

    call get_eig(NN_TABLE(1), PKPTS_PY%kpoint, PKPTS_PY%nkpoint, PINPT, PPRAM, &
                 ETBA_PY%E, ETBA_PY%V, ETBA_PY%SV, PGEOM_PY%neig, PGEOM_PY%init_erange, PGEOM_PY%nband, &
                 PINPT%flag_get_orbital, PINPT%flag_sparse, flag_stat, PINPT%flag_phase)
    call get_eig(NN_TABLE(2), PKPTS_PY2%kpoint, PKPTS_PY2%nkpoint, PINPT, PPRAM, &
                 ETBA_PY2%E, ETBA_PY2%V, ETBA_PY2%SV, PGEOM_PY2%neig, PGEOM_PY2%init_erange, PGEOM_PY2%nband, &
                 PINPT%flag_get_orbital, PINPT%flag_sparse, flag_stat, PINPT%flag_phase)
    call get_eig(NN_TABLE(3), PKPTS_PY3%kpoint, PKPTS_PY3%nkpoint, PINPT, PPRAM, &
                 ETBA_PY3%E, ETBA_PY3%V, ETBA_PY3%SV, PGEOM_PY3%neig, PGEOM_PY3%init_erange, PGEOM_PY3%nband, &
                 PINPT%flag_get_orbital, PINPT%flag_sparse, flag_stat, PINPT%flag_phase)

    if(PINPT%flag_get_orbital) then
      call copy_kpoints(PKPTS_PY,    PKPTS(1), 2)
      call copy_poscar(PGEOM_PY,     PGEOM(1), 2)
      call copy_energy(ETBA_PY,      ETBA(1),  2)
      call get_orbital_projection(ETBA(1), PKPTS(1), PINPT, PGEOM(1), PPRAM%flag_use_overlap)
      call copy_energy(ETBA_PY,      ETBA(1),  1)

      call copy_kpoints(PKPTS_PY2,    PKPTS(2), 2)
      call copy_poscar(PGEOM_PY2,     PGEOM(2), 2)
      call copy_energy(ETBA_PY2,      ETBA(2),  2)
      call get_orbital_projection(ETBA(2), PKPTS(2), PINPT, PGEOM(2), PPRAM%flag_use_overlap)
      call copy_energy(ETBA_PY2,      ETBA(2),  1)

      call copy_kpoints(PKPTS_PY3,    PKPTS(3), 2)
      call copy_poscar(PGEOM_PY3,     PGEOM(3), 2)
      call copy_energy(ETBA_PY3,      ETBA(3),  2)
      call get_orbital_projection(ETBA(3), PKPTS(3), PINPT, PGEOM(3), PPRAM%flag_use_overlap)
      call copy_energy(ETBA_PY3,      ETBA(3),  1)
    endif

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_hopping(NN_TABLE_PY, NN_TABLE(1), 1)
    call copy_hopping(NN_TABLE_PY2, NN_TABLE(2), 1)
    call copy_hopping(NN_TABLE_PY3, NN_TABLE(3), 1)
    return
endsubroutine

subroutine eig2(comm, PINPT_PY, PPRAM_PY, PKPTS_PY, PGEOM_PY, NN_TABLE_PY, ETBA_PY, &
                                          PKPTS_PY2, PGEOM_PY2, NN_TABLE_PY2, ETBA_PY2)
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY   ,PKPTS_PY2
    type(poscar_py),      intent(inout)         :: PGEOM_PY   ,PGEOM_PY2
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY,NN_TABLE_PY2
    type(energy_py),      intent(inout)         :: ETBA_PY    ,ETBA_PY2
    integer,              intent(inout)         :: comm
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(hopping ), dimension(PINPT_PY%nsystem) :: NN_TABLE
    type(energy  ), dimension(PINPT_PY%nsystem) :: ETBA
    type(kpoints ), dimension(PINPT_PY%nsystem) :: PKPTS
    type(poscar  ), dimension(PINPT_PY%nsystem) :: PGEOM
    logical                                        flag_stat

#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif
    call set_verbose(PINPT_PY%iverbose)
    if(iverbose .eq. 1) then
        flag_stat = .TRUE.
    else
        flag_stat = .FALSE.
    endif

    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_hopping(NN_TABLE_PY , NN_TABLE(1), 2)
    call copy_hopping(NN_TABLE_PY2, NN_TABLE(2), 2)

    call get_eig(NN_TABLE(1), PKPTS_PY%kpoint, PKPTS_PY%nkpoint, PINPT, PPRAM, &
                 ETBA_PY%E, ETBA_PY%V, ETBA_PY%SV, PGEOM_PY%neig, PGEOM_PY%init_erange, PGEOM_PY%nband, &
                 PINPT%flag_get_orbital, PINPT%flag_sparse, flag_stat, PINPT%flag_phase)
    call get_eig(NN_TABLE(2), PKPTS_PY2%kpoint, PKPTS_PY2%nkpoint, PINPT, PPRAM, &
                 ETBA_PY2%E, ETBA_PY2%V, ETBA_PY2%SV, PGEOM_PY2%neig, PGEOM_PY2%init_erange, PGEOM_PY2%nband, &
                 PINPT%flag_get_orbital, PINPT%flag_sparse, flag_stat, PINPT%flag_phase)
    
    if(PINPT%flag_get_orbital) then
      call copy_kpoints(PKPTS_PY,    PKPTS(1), 2)
      call copy_poscar(PGEOM_PY,     PGEOM(1), 2)
      call copy_energy(ETBA_PY,      ETBA(1),  2)
      call get_orbital_projection(ETBA(1), PKPTS(1), PINPT, PGEOM(1), PPRAM%flag_use_overlap)
      call copy_energy(ETBA_PY,      ETBA(1),  1)

      call copy_kpoints(PKPTS_PY2,    PKPTS(2), 2)
      call copy_poscar(PGEOM_PY2,     PGEOM(2), 2)
      call copy_energy(ETBA_PY2,      ETBA(2),  2)
      call get_orbital_projection(ETBA(2), PKPTS(2), PINPT, PGEOM(2), PPRAM%flag_use_overlap)
      call copy_energy(ETBA_PY2,      ETBA(2),  1)
    endif

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_hopping(NN_TABLE_PY, NN_TABLE(1), 1)
    call copy_hopping(NN_TABLE_PY2, NN_TABLE(2), 1)
    return
endsubroutine

subroutine eig(comm, PINPT_PY, PPRAM_PY, PKPTS_PY, PGEOM_PY, NN_TABLE_PY, ETBA_PY)
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY
    type(poscar_py),      intent(inout)         :: PGEOM_PY
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY
    type(energy_py),      intent(inout)         :: ETBA_PY
    integer,              intent(inout)         :: comm
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(hopping )                              :: NN_TABLE
    type(energy  )                              :: ETBA     
    type(kpoints )                              :: PKPTS
    type(poscar  )                              :: PGEOM 
    logical                                        flag_stat

#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif
    call set_verbose(PINPT_PY%iverbose)
    if(iverbose .eq. 1) then
        flag_stat = .TRUE. 
    else
        flag_stat = .FALSE.
    endif

    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_hopping(NN_TABLE_PY, NN_TABLE, 2)
    call get_eig(NN_TABLE, PKPTS_PY%kpoint, PKPTS_PY%nkpoint, PINPT, PPRAM, &
                 ETBA_PY%E, ETBA_PY%V, ETBA_PY%SV, PGEOM_PY%neig, PGEOM_PY%init_erange, PGEOM_PY%nband, &
                 PINPT%flag_get_orbital, PINPT%flag_sparse, flag_stat, PINPT%flag_phase)

    if(PINPT%flag_get_orbital) then 
      call copy_kpoints(PKPTS_PY,    PKPTS, 2)
      call copy_poscar(PGEOM_PY,     PGEOM, 2)
      call copy_energy(ETBA_PY,      ETBA,  2)
      call get_orbital_projection(ETBA, PKPTS, PINPT, PGEOM, PPRAM%flag_use_overlap)
      call copy_energy(ETBA_PY,      ETBA,  1)
    endif

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_hopping(NN_TABLE_PY, NN_TABLE, 1)
    return
endsubroutine

subroutine toten(comm, PINPT_PY, PKPTS_PY, PGEOM_PY, ETBA_PY)
    use total_energy
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY
    type(poscar_py),      intent(inout)         :: PGEOM_PY
    type(energy_py),      intent(inout)         :: ETBA_PY
    integer,              intent(in)            :: comm
    type(incar   )                              :: PINPT
    type(kpoints )                              :: PKPTS
    type(poscar  )                              :: PGEOM
    type(energy  )                              :: ETBA 
    logical                                        flag_report
    logical                                        flag_opt_ef
    integer(kind=sp)                               nelect(2)
    integer(kind=sp)                               mpierr
#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif

    call set_verbose(PINPT_PY%iverbose)

    if(PGEOM_PY%nelect(1) .le. 0.d0) then
       write(6,'(A)')'    !WARN! The total nergy calculation is requested but number of electons are not explicitly specified. '
       write(6,'(A)')'           Please check "nelect" tag for the information. Exit...'
       kill_job
    else
       if(.not. PINPT_PY%flag_get_total_energy) then
         if(allocated(PGEOM_PY%nelect)) then
            nelect = PGEOM_PY%nelect
            deallocate(PGEOM_PY%nelect)
         endif
         if(PINPT_PY%nspin .eq. 2) then
            allocate(PGEOM_PY%nelect(2))
            PGEOM_PY%nelect = nelect
         else
            allocate(PGEOM_PY%nelect(1))
            PGEOM_PY%nelect = nelect(1)
         endif
         PINPT_PY%flag_get_total_energy = .TRUE.
         if(.not. allocated(ETBA_PY%F_OCC)) allocate(ETBA_PY%F_OCC(PGEOM_PY%nband*PINPT_PY%nspin, PKPTS_PY%nkpoint))
         if(.not. allocated(ETBA_PY%E_BAND))allocate(ETBA_PY%E_BAND(PINPT_PY%nspin))
         if(.not. allocated(ETBA_PY%E_TOT)) allocate(ETBA_PY%E_TOT (PINPT_PY%nspin))
       endif
    endif

    call copy_incar(  PINPT_PY, PINPT, 2)
    call copy_kpoints(PKPTS_PY, PKPTS, 2)
    call copy_poscar( PGEOM_PY, PGEOM, 2)
    call copy_energy(  ETBA_PY, ETBA , 2)
    if(iverbose .eq. 1) then
        flag_report = .TRUE.  ! but not be written if iverbose = 2
    else
        flag_report = .FALSE. ! but not be written if iverbose = 2
    endif

    flag_opt_ef = .TRUE.
    call get_total_energy(ETBA, PINPT, PKPTS, PGEOM, flag_opt_ef, flag_report)

    call copy_energy(  ETBA_PY, ETBA , 1)
    
    return
endsubroutine

subroutine fit3(comm, PINPT_PY, PPRAM_PY, PKPTS_PY1, PWGHT_PY1, PGEOM_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1,&
                                          PKPTS_PY2, PWGHT_PY2, PGEOM_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2,&
                                          PKPTS_PY3, PWGHT_PY3, PGEOM_PY3, NN_TABLE_PY3, EDFT_PY3, ETBA_PY3)
    use read_incar
    use cost_function, only: get_dE12
    use set_default
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY1   , PKPTS_PY2   , PKPTS_PY3    
    type(weight_py),      intent(inout)         :: PWGHT_PY1   , PWGHT_PY2   , PWGHT_PY3
    type(poscar_py),      intent(inout)         :: PGEOM_PY1   , PGEOM_PY2   , PGEOM_PY3
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY1, NN_TABLE_PY2, NN_TABLE_PY3
    type(energy_py),      intent(inout)         :: EDFT_PY1    , EDFT_PY2    , EDFT_PY3
    type(energy_py),      intent(inout)         :: ETBA_PY1    , ETBA_PY2    , ETBA_PY3
    integer,              intent(in)            :: comm
    integer(kind=sp)                               i, j, k
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(kpoints ), dimension(3               ) :: PKPTS
    type(weight  ), dimension(3               ) :: PWGHT
    type(poscar  ), dimension(3               ) :: PGEOM
    type(hopping ), dimension(3               ) :: NN_TABLE
    type(energy  ), dimension(3               ) :: EDFT, ETBA
    integer(kind=sp)                               ifit
    logical                                        flag_exit
    real(kind=dp)                                  fnorm, fnorm_plain, fnorm_orb,  fnorm_
    external                                       get_eig
    logical                                        flag_order, flag_get_orbital, flag_write_info
    integer(kind=sp)                               mpierr
    integer(kind=sp)                               info, maxfev, ldjac, imode, nparam_free
    real(kind=dp)                                  epsfcn, factor, xtol, ftol, gtol
#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif
    fnorm_ = 0d0        ; flag_exit = .false. ;    ifit = 0
    fnorm_plain = 0d0   ; fnorm_orb = 0d0
    flag_order       = PINPT_PY%flag_get_band_order .and. (.not. PINPT_PY%flag_get_band_order_print_only)
    flag_get_orbital = (PWGHT_PY1%flag_weight_orb .or. flag_order .or. PINPT_PY%flag_fit_orbital)

    if(allocated(PPRAM_PY%cost_history)) deallocate(PPRAM_PY%cost_history)
    allocate(PPRAM_PY%cost_history(PINPT_PY%miter))
    PPRAM_PY%cost_history = 0d0
    flag_write_info = .false.
    call set_verbose(PINPT_PY%iverbose)
    if(iverbose .eq. 1) flag_write_info = .true.

    call copy_incar(PINPT_PY, PINPT, 2) ; PINPT%nsystem=3
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_all(PKPTS_PY1, PGEOM_PY1, PWGHT_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                    PKPTS(1) , PGEOM(1) , PWGHT(1) , NN_TABLE(1) , EDFT(1) , ETBA(1) , 2) ! system-1
    call copy_all(PKPTS_PY2, PGEOM_PY2, PWGHT_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                    PKPTS(2) , PGEOM(2) , PWGHT(2) , NN_TABLE(2) , EDFT(2) , ETBA(2) , 2) ! system-2
    call copy_all(PKPTS_PY3, PGEOM_PY3, PWGHT_PY3, NN_TABLE_PY3, EDFT_PY3, ETBA_PY3, &
                    PKPTS(3) , PGEOM(3) , PWGHT(3) , NN_TABLE(3) , EDFT(3) , ETBA(3) , 2) ! system-3
    
    nparam_free = PPRAM%nparam_free ! total number of free parameters. Note: nrl parameterization is not supported yet.
    imode  = 13 !
    factor = 100.0D+00
    maxfev = PINPT_PY%miter * ( PPRAM_PY%nparam_free + 1 )
    ftol = PINPT_PY%ftol    ;xtol = PINPT_PY%ptol ; gtol = 0.0D+00; epsfcn = 0.000D+00
    ldjac = 0

    if(imode .eq. 13) then
      do i = 1, PINPT%nsystem
          ldjac = ldjac + PKPTS(i)%nkpoint * PGEOM(i)%nband * PINPT%nspin
      enddo
    elseif(imode .eq. 2) then
      ldjac = sum(PKPTS(:)%nkpoint)
    endif

    call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm, fnorm_plain, fnorm_orb, maxfev, epsfcn, factor, info, flag_write_info, .false. )
    call check_conv_and_constraint(PPRAM, PINPT, flag_exit, ifit, fnorm, fnorm_)

    ! before return calculate again with final parameter to get get energy and dE
    call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, PINPT%flag_fit_orbital )

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_energy(ETBA_PY1, ETBA(1), 1)
    call copy_energy(ETBA_PY2, ETBA(2), 1)
    call copy_energy(ETBA_PY3, ETBA(3), 1)

    return
endsubroutine

subroutine fit2(comm, PINPT_PY, PPRAM_PY, PKPTS_PY1, PWGHT_PY1, PGEOM_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1,&
                                          PKPTS_PY2, PWGHT_PY2, PGEOM_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2)
    use read_incar
    use cost_function, only: get_dE12
    use set_default
    type(incar_py),       intent(inout)         :: PINPT_PY      
    type(params_py),      intent(inout)         :: PPRAM_PY     
    type(kpoints_py),     intent(inout)         :: PKPTS_PY1   , PKPTS_PY2
    type(weight_py),      intent(inout)         :: PWGHT_PY1   , PWGHT_PY2    
    type(poscar_py),      intent(inout)         :: PGEOM_PY1   , PGEOM_PY2    
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY1, NN_TABLE_PY2 
    type(energy_py),      intent(inout)         :: EDFT_PY1    , EDFT_PY2     
    type(energy_py),      intent(inout)         :: ETBA_PY1    , ETBA_PY2     
    integer,              intent(in)            :: comm         
    integer(kind=sp)                               i, j, k
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(kpoints ), dimension(2               ) :: PKPTS
    type(weight  ), dimension(2               ) :: PWGHT
    type(poscar  ), dimension(2               ) :: PGEOM
    type(hopping ), dimension(2               ) :: NN_TABLE
    type(energy  ), dimension(2               ) :: EDFT, ETBA
    integer(kind=sp)                               ifit
    logical                                        flag_exit
    real(kind=dp)                                  fnorm, fnorm_plain, fnorm_orb,  fnorm_
    external                                       get_eig
    logical                                        flag_order, flag_get_orbital, flag_write_info
    integer(kind=sp)                               mpierr
    integer(kind=sp)                               info, maxfev, ldjac, imode, nparam_free
    real(kind=dp)                                  epsfcn, factor, xtol, ftol, gtol
#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif
    fnorm_ = 0d0        ; flag_exit = .false. ;    ifit = 0
    fnorm_plain = 0d0   ; fnorm_orb = 0d0
    flag_order       = PINPT_PY%flag_get_band_order .and. (.not. PINPT_PY%flag_get_band_order_print_only)
    flag_get_orbital = (PWGHT_PY1%flag_weight_orb .or. flag_order .or. PINPT_PY%flag_fit_orbital)

    if(allocated(PPRAM_PY%cost_history)) deallocate(PPRAM_PY%cost_history)
    allocate(PPRAM_PY%cost_history(PINPT_PY%miter))
    PPRAM_PY%cost_history = 0d0
    flag_write_info = .false.
    call set_verbose(PINPT_PY%iverbose)
    if(iverbose .eq. 1) flag_write_info = .true.

    call copy_incar(PINPT_PY, PINPT, 2) ; PINPT%nsystem=2
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_all(PKPTS_PY1, PGEOM_PY1, PWGHT_PY1, NN_TABLE_PY1, EDFT_PY1, ETBA_PY1, &
                    PKPTS(1) , PGEOM(1) , PWGHT(1) , NN_TABLE(1) , EDFT(1) , ETBA(1) , 2) ! system-1
    call copy_all(PKPTS_PY2, PGEOM_PY2, PWGHT_PY2, NN_TABLE_PY2, EDFT_PY2, ETBA_PY2, &
                    PKPTS(2) , PGEOM(2) , PWGHT(2) , NN_TABLE(2) , EDFT(2) , ETBA(2) , 2) ! system-2
    
    nparam_free = PPRAM%nparam_free ! total number of free parameters. Note: nrl parameterization is not supported yet.
    imode  = 13 !
    factor = 100.0D+00
    maxfev = PINPT_PY%miter * ( PPRAM_PY%nparam_free + 1 )
    ftol = PINPT_PY%ftol    ;xtol = PINPT_PY%ptol ; gtol = 0.0D+00; epsfcn = 0.000D+00
    ldjac = 0

    if(imode .eq. 13) then
      do i = 1, PINPT%nsystem
          ldjac = ldjac + PKPTS(i)%nkpoint * PGEOM(i)%nband * PINPT%nspin
      enddo
    elseif(imode .eq. 2) then
      ldjac = sum(PKPTS(:)%nkpoint)
    endif

    call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm, fnorm_plain, fnorm_orb, maxfev, epsfcn, factor, info, flag_write_info, .false. )
    call check_conv_and_constraint(PPRAM, PINPT, flag_exit, ifit, fnorm, fnorm_)

    ! before return calculate again with final parameter to get get energy and dE
    call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, PINPT%flag_fit_orbital )

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_energy(ETBA_PY1, ETBA(1), 1)
    call copy_energy(ETBA_PY2, ETBA(2), 1)

    return
endsubroutine

subroutine fit(comm, PINPT_PY, PPRAM_PY, PKPTS_PY, PWGHT_PY, PGEOM_PY, NN_TABLE_PY, &
               EDFT_PY, ETBA_PY)
    use read_incar
    use cost_function, only: get_dE12
    use set_default
    type(incar_py),       intent(inout)         :: PINPT_PY
    type(params_py),      intent(inout)         :: PPRAM_PY
    type(kpoints_py),     intent(inout)         :: PKPTS_PY
    type(weight_py),      intent(inout)         :: PWGHT_PY
    type(poscar_py),      intent(inout)         :: PGEOM_PY
    type(hopping_py),     intent(inout)         :: NN_TABLE_PY
    type(energy_py),      intent(inout)         :: EDFT_PY
    type(energy_py),      intent(inout)         :: ETBA_PY
    integer,              intent(in)            :: comm
    integer(kind=sp)                               i, j, k
    type(incar   )                              :: PINPT
    type(params  )                              :: PPRAM
    type(kpoints ), dimension(PINPT_PY%nsystem) :: PKPTS
    type(weight  ), dimension(PINPT_PY%nsystem) :: PWGHT
    type(poscar  ), dimension(PINPT_PY%nsystem) :: PGEOM
    type(hopping ), dimension(PINPT_PY%nsystem) :: NN_TABLE
    type(energy  ), dimension(PINPT_PY%nsystem) :: EDFT, ETBA
    integer(kind=sp)                               ifit
    logical                                        flag_exit
    real(kind=dp)                                  fnorm, fnorm_plain, fnorm_orb,  fnorm_
    external                                       get_eig
    logical                                        flag_order, flag_get_orbital, flag_write_info
    integer(kind=sp)                               mpierr 
    integer(kind=sp)                               info, maxfev, ldjac, imode, nparam_free
    real(kind=dp)                                  epsfcn, factor, xtol, ftol, gtol
#ifdef MPI
    mpi_comm_earth = comm
    call get_my_task()
#endif

    fnorm_ = 0d0        ; flag_exit = .false. ;    ifit = 0
    fnorm_plain = 0d0   ; fnorm_orb = 0d0
    flag_order       = PINPT_PY%flag_get_band_order .and. (.not. PINPT_PY%flag_get_band_order_print_only)
    flag_get_orbital = (PWGHT_PY%flag_weight_orb .or. flag_order .or. PINPT_PY%flag_fit_orbital)

    if(allocated(PPRAM_PY%cost_history)) deallocate(PPRAM_PY%cost_history)
    allocate(PPRAM_PY%cost_history(PINPT_PY%miter))
    PPRAM_PY%cost_history = 0d0
    flag_write_info = .false.
    call set_verbose(PINPT_PY%iverbose)
    if(iverbose .eq. 1) flag_write_info = .true.

    ! init system
    ! NOTE: in this subroutine, we assume that nsystem=1. nsystem>1 case, it will be updated in near future. H.-J. Kim. (01. Feb. 2021)
    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_params(PPRAM_PY, PPRAM, 2)
    call copy_kpoints(PKPTS_PY, PKPTS(1), 2)
    call copy_poscar(PGEOM_PY, PGEOM(1), 2)
    call copy_weight(PWGHT_PY, PWGHT(1), 2)
    call copy_hopping(NN_TABLE_PY, NN_TABLE(1), 2)
    call copy_energy(EDFT_PY, EDFT(1), 2)
    call copy_energy(ETBA_PY, ETBA(1), 2)
    
    nparam_free = PPRAM%nparam_free ! total number of free parameters. Note: nrl parameterization is not supported yet.
    imode  = 13 !
    factor = 100.0D+00
    maxfev = PINPT_PY%miter * ( PPRAM_PY%nparam_free + 1 )
    ftol = PINPT_PY%ftol    ;xtol = PINPT_PY%ptol ; gtol = 0.0D+00; epsfcn = 0.000D+00    
    ldjac = 0

    if(imode .eq. 13) then
      do i = 1, PINPT_PY%nsystem
          ldjac = ldjac + PKPTS(i)%nkpoint * PGEOM(i)%nband * PINPT%nspin
      enddo
    elseif(imode .eq. 2) then
      ldjac = sum(PKPTS(:)%nkpoint)
    endif

    call lmdif(get_eig, NN_TABLE, ldjac, imode, PINPT, PPRAM, PKPTS, PGEOM, EDFT, nparam_free, PWGHT, &
               ftol, xtol, gtol, fnorm, fnorm_plain, fnorm_orb, maxfev, epsfcn, factor, info, flag_write_info, .false. )
    call check_conv_and_constraint(PPRAM, PINPT, flag_exit, ifit, fnorm, fnorm_)

    ! before return calculate again with final parameter to get get energy and dE
    call get_dE12(PINPT, PPRAM, NN_TABLE, EDFT, ETBA, PWGHT, PGEOM, PKPTS, PINPT%flag_fit_orbital )

    call copy_incar(PINPT_PY, PINPT, 1)
    call copy_params(PPRAM_PY, PPRAM, 1)
    call copy_energy(ETBA_PY, ETBA(1), 1)

#ifdef MPI
    call MPI_BARRIER(mpi_comm_earth, mpierr)
#endif

    return
end subroutine fit

function init_incar_py(ifilenm, nsystem) result(PINPT_PY)
    type(incar_py)                         PINPT_PY
    type(incar   )                      :: PINPT
    integer(kind=sp),     intent(in)    :: nsystem
    character(len=132),   intent(in)    :: ifilenm(nsystem)
    integer(kind=sp)                       i, j, mpierr
    logical                                flag_ifile_exist

    iverbose = 2 ! 1: full, 2: no
    print_mode = 99  ! default verbosity 
 
    PINPT_PY%nsystem = nsystem
    PINPT_PY%flag_python_module = .TRUE.

    if(PINPT_PY%nsystem .eq. 1) then
      allocate(PINPT_PY%ifilenm(PINPT_PY%nsystem))
      allocate(PINPT_PY%title(PINPT_PY%nsystem))
      PINPT_PY%ifilenm(1) = trim(ifilenm(1))
      PINPT_PY%title(1) = ''
    elseif(PINPT_PY%nsystem .ge. 2) then
      allocate(PINPT_PY%ifilenm(PINPT_PY%nsystem))
      allocate(PINPT_PY%title(PINPT_PY%nsystem))
      do i = 1, PINPT_PY%nsystem
        PINPT_PY%ifilenm(i) = trim(ifilenm(i))
        write(PINPT_PY%title(i),'(A,I0)') '.',i
      enddo
    endif

    do i = 1, PINPT_PY%nsystem
      inquire(file=trim(PINPT_PY%ifilenm(i)),exist=flag_ifile_exist)
      if(.not. flag_ifile_exist) then
        write(6,'(A,A,A)')'    !WARN! Input file:',trim(PINPT_PY%ifilenm(i)),' does not exist!! Exit...' ; write_msg
        kill_job
      endif
    enddo

    PINPT_PY%flag_python_module                      = .TRUE. 
    PINPT_PY%flag_get_band                           = .FALSE.
    PINPT_PY%flag_fit_degeneracy                     = .FALSE.
    PINPT_PY%flag_report_geom                        = .FALSE.
    PINPT_PY%flag_phase                              = .TRUE. 
    PINPT_PY%flag_tbfit                              = .TRUE. 
    PINPT_PY%flag_tbfit_finish                       = .FALSE.
    PINPT_PY%flag_print_only_target                  = .FALSE.
    PINPT_PY%flag_print_energy_diff                  = .TRUE. 
    PINPT_PY%flag_print_orbital                      = .FALSE.
    PINPT_PY%flag_get_orbital                        = .FALSE.
    PINPT_PY%flag_print_mag                          = .FALSE.
    PINPT_PY%flag_print_single                       = .FALSE.
    PINPT_PY%flag_local_charge                       = .FALSE.
    PINPT_PY%flag_collinear                          = .FALSE.
    PINPT_PY%flag_noncollinear                       = .FALSE.
    PINPT_PY%flag_soc                                = .FALSE.
    PINPT_PY%flag_plus_U                             = .FALSE.
    PINPT_PY%flag_scissor                            = .FALSE.
    PINPT_PY%flag_print_proj                         = .FALSE.
    PINPT_PY%flag_print_proj_sum                     = .FALSE.
    PINPT_PY%flag_print_proj_atom                    = .FALSE.
    PINPT_PY%flag_plot_fit                           = .FALSE.
    PINPT_PY%flag_plot                               = .FALSE.
    PINPT_PY%flag_get_band_order                     = .FALSE.
    PINPT_PY%flag_get_band_order_print_only          = .FALSE.
    PINPT_PY%flag_use_weight                         = .FALSE.
    PINPT_PY%flag_fit_orbital                        = .FALSE.
    PINPT_PY%flag_fit_orbital_parse                  = .FALSE.
    PINPT_PY%flag_get_total_energy                   = .FALSE.
    PINPT_PY%flag_pso_with_lmdif                     = .FALSE.
    PINPT_PY%flag_get_unfold                         = .FALSE.
    PINPT_PY%iverbose                                = 2 ! do not print : 2, print all: 1
    PINPT_PY%ptol                                    = 0.00001d0
    PINPT_PY%ftol                                    = 0.00001d0
    PINPT_PY%fdiff                                   = 0.001d0
    PINPT_PY%electronic_temperature                  = 11.60452251451435653502d0 ! electronic temp T; default, so that k_B * T = 0.001 eV
    PINPT_PY%lmmax                                   = 9
    PINPT_PY%miter                                   = 30
    PINPT_PY%mxfit                                   = 1
    PINPT_PY%nn_max(3)                               = 3
    PINPT_PY%ispin                                   = 1
    PINPT_PY%ispinor                                 = 1
    PINPT_PY%nspin                                   = 1
    PINPT_PY%fnamelog                                = 'TBFIT.out'
    PINPT_PY%orbital_fit_smearing                    = 0.1d0 ! default
    return
endfunction

subroutine copy_incar(PINPT_PY, PINPT, imode)
    type(incar    )                        PINPT
    type(incar_py)                         PINPT_PY
    integer(kind=sp)                       imode, n1, n2

    if(imode .eq. 1) then
       PINPT_PY%flag_python_module                 =      PINPT%flag_python_module             
       PINPT_PY%flag_get_band                      =      PINPT%flag_get_band
       PINPT_PY%flag_fit_degeneracy                =      PINPT%flag_fit_degeneracy
       PINPT_PY%flag_report_geom                   =      PINPT%flag_report_geom
       PINPT_PY%flag_phase                         =      PINPT%flag_phase
       PINPT_PY%flag_tbfit                         =      PINPT%flag_tbfit
       PINPT_PY%flag_tbfit_finish                  =      PINPT%flag_tbfit_finish
       PINPT_PY%flag_print_only_target             =      PINPT%flag_print_only_target
       PINPT_PY%flag_print_energy_diff             =      PINPT%flag_print_energy_diff
       PINPT_PY%flag_print_orbital                 =      PINPT%flag_print_orbital
       PINPT_PY%flag_get_orbital                   =      PINPT%flag_get_orbital
       PINPT_PY%flag_print_mag                     =      PINPT%flag_print_mag
       PINPT_PY%flag_print_single                  =      PINPT%flag_print_single
       PINPT_PY%flag_local_charge                  =      PINPT%flag_local_charge
       PINPT_PY%flag_collinear                     =      PINPT%flag_collinear
       PINPT_PY%flag_noncollinear                  =      PINPT%flag_noncollinear
       PINPT_PY%flag_soc                           =      PINPT%flag_soc
       PINPT_PY%flag_plus_U                        =      PINPT%flag_plus_U
       PINPT_PY%flag_scissor                       =      PINPT%flag_scissor
       PINPT_PY%flag_print_proj                    =      PINPT%flag_print_proj
       PINPT_PY%flag_print_proj_sum                =      PINPT%flag_print_proj_sum
       PINPT_PY%flag_print_proj_atom               =      PINPT%flag_print_proj_atom
       PINPT_PY%flag_plot_fit                      =      PINPT%flag_plot_fit
       PINPT_PY%flag_plot                          =      PINPT%flag_plot
       PINPT_PY%flag_get_band_order                =      PINPT%flag_get_band_order
       PINPT_PY%flag_get_band_order_print_only     =      PINPT%flag_get_band_order_print_only
       PINPT_PY%flag_use_weight                    =      PINPT%flag_use_weight
       PINPT_PY%flag_fit_orbital                   =      PINPT%flag_fit_orbital
       PINPT_PY%flag_fit_orbital_parse             =      PINPT%flag_fit_orbital_parse
       PINPT_PY%flag_get_total_energy              =      PINPT%flag_get_total_energy
       PINPT_PY%flag_pso_with_lmdif                =      PINPT%flag_pso_with_lmdif
       PINPT_PY%flag_get_unfold                    =      PINPT%flag_get_unfold    
       PINPT_PY%ptol                               =      PINPT%ptol
       PINPT_PY%ftol                               =      PINPT%ftol
       PINPT_PY%fdiff                              =      PINPT%fdiff
       PINPT_PY%r_scissor                          =      PINPT%r_scissor
       PINPT_PY%electronic_temperature             =      PINPT%electronic_temperature
       PINPT_PY%i_scissor                          =      PINPT%i_scissor
       PINPT_PY%nsystem                            =      PINPT%nsystem
       PINPT_PY%lmmax                              =      PINPT%lmmax
       PINPT_PY%miter                              =      PINPT%miter
       PINPT_PY%mxfit                              =      PINPT%mxfit
       PINPT_PY%nn_max                             =      PINPT%nn_max
       PINPT_PY%ispin                              =      PINPT%ispin
       PINPT_PY%ispinor                            =      PINPT%ispinor
       PINPT_PY%nspin                              =      PINPT%nspin
       PINPT_PY%ls_type                            =      PINPT%ls_type
       PINPT_PY%fnamelog                           =      PINPT%fnamelog
       PINPT_PY%nproj_sum                          =      PINPT%nproj_sum
       PINPT_PY%orbital_fit_smearing               =      PINPT%orbital_fit_smearing

       if(allocated( PINPT_PY%ifilenm       ))         deallocate(PINPT_PY%ifilenm)
       if(allocated( PINPT_PY%title         ))         deallocate(PINPT_PY%title  )
       if(allocated( PINPT_PY%proj_natom    ))         deallocate(PINPT_PY%proj_natom)
       if(allocated( PINPT_PY%proj_atom     ))         deallocate(PINPT_PY%proj_atom)

       if(allocated( PINPT%ifilenm )) then
         n1 = size(PINPT%ifilenm) 
         allocate( PINPT_PY%ifilenm(n1)) ; PINPT_PY%ifilenm = PINPT%ifilenm
       endif
       if(allocated( PINPT%title )) then
         n1 = size(PINPT%title)  
         allocate( PINPT_PY%title(n1)) ; PINPT_PY%title = PINPT%title
       endif
       if(allocated( PINPT%proj_atom )) then
         n1 = size(PINPT%proj_atom(:,1)) ; n2 = size(PINPT%proj_atom(1,:))
         allocate( PINPT_PY%proj_atom(n1,n2)) ; PINPT_PY%proj_atom = PINPT%proj_atom
       endif
       if(allocated( PINPT%proj_natom )) then
         n1 = size(PINPT%proj_natom)
         allocate( PINPT_PY%proj_natom(n1)) ; PINPT_PY%proj_natom = PINPT%proj_natom
       endif


    elseif( imode .eq. 2) then
       PINPT%flag_python_module                 =      PINPT_PY%flag_python_module             
       PINPT%flag_get_band                      =      PINPT_PY%flag_get_band
       PINPT%flag_fit_degeneracy                =      PINPT_PY%flag_fit_degeneracy
       PINPT%flag_report_geom                   =      PINPT_PY%flag_report_geom
       PINPT%flag_phase                         =      PINPT_PY%flag_phase
       PINPT%flag_tbfit                         =      PINPT_PY%flag_tbfit
       PINPT%flag_tbfit_finish                  =      PINPT_PY%flag_tbfit_finish
       PINPT%flag_print_only_target             =      PINPT_PY%flag_print_only_target
       PINPT%flag_print_energy_diff             =      PINPT_PY%flag_print_energy_diff
       PINPT%flag_print_orbital                 =      PINPT_PY%flag_print_orbital
       PINPT%flag_get_orbital                   =      PINPT_PY%flag_get_orbital
       PINPT%flag_print_mag                     =      PINPT_PY%flag_print_mag
       PINPT%flag_print_single                  =      PINPT_PY%flag_print_single
       PINPT%flag_local_charge                  =      PINPT_PY%flag_local_charge
       PINPT%flag_collinear                     =      PINPT_PY%flag_collinear
       PINPT%flag_noncollinear                  =      PINPT_PY%flag_noncollinear
       PINPT%flag_soc                           =      PINPT_PY%flag_soc
       PINPT%flag_plus_U                        =      PINPT_PY%flag_plus_U
       PINPT%flag_scissor                       =      PINPT_PY%flag_scissor
       PINPT%flag_print_proj                    =      PINPT_PY%flag_print_proj
       PINPT%flag_print_proj_sum                =      PINPT_PY%flag_print_proj_sum
       PINPT%flag_print_proj_atom               =      PINPT_PY%flag_print_proj_atom
       PINPT%flag_plot_fit                      =      PINPT_PY%flag_plot_fit
       PINPT%flag_plot                          =      PINPT_PY%flag_plot
       PINPT%flag_get_band_order                =      PINPT_PY%flag_get_band_order
       PINPT%flag_get_band_order_print_only     =      PINPT_PY%flag_get_band_order_print_only
       PINPT%flag_use_weight                    =      PINPT_PY%flag_use_weight
       PINPT%flag_fit_orbital                   =      PINPT_PY%flag_fit_orbital
       PINPT%flag_fit_orbital_parse             =      PINPT_PY%flag_fit_orbital_parse
       PINPT%flag_get_total_energy              =      PINPT_PY%flag_get_total_energy
       PINPT%flag_pso_with_lmdif                =      PINPT_PY%flag_pso_with_lmdif
       PINPT%flag_get_unfold                    =      PINPT_PY%flag_get_unfold     
       PINPT%ptol                               =      PINPT_PY%ptol
       PINPT%ftol                               =      PINPT_PY%ftol
       PINPT%fdiff                              =      PINPT_PY%fdiff
       PINPT%r_scissor                          =      PINPT_PY%r_scissor
       PINPT%electronic_temperature             =      PINPT_PY%electronic_temperature
       PINPT%i_scissor                          =      PINPT_PY%i_scissor
       PINPT%nsystem                            =      PINPT_PY%nsystem
       PINPT%lmmax                              =      PINPT_PY%lmmax
       PINPT%miter                              =      PINPT_PY%miter
       PINPT%mxfit                              =      PINPT_PY%mxfit
       PINPT%nn_max                             =      PINPT_PY%nn_max
       PINPT%ispin                              =      PINPT_PY%ispin
       PINPT%ispinor                            =      PINPT_PY%ispinor
       PINPT%nspin                              =      PINPT_PY%nspin
       PINPT%ls_type                            =      PINPT_PY%ls_type
       PINPT%fnamelog                           =      PINPT_PY%fnamelog
       PINPT%nproj_sum                          =      PINPT_PY%nproj_sum
       PINPT%orbital_fit_smearing               =      PINPT_PY%orbital_fit_smearing

       if(allocated( PINPT%ifilenm          ))         deallocate(PINPT%ifilenm)
       if(allocated( PINPT%title            ))         deallocate(PINPT%title  )
       if(allocated( PINPT%proj_natom       ))         deallocate(PINPT%proj_natom)
       if(allocated( PINPT%proj_atom        ))         deallocate(PINPT%proj_atom)

       if(allocated( PINPT_PY%ifilenm )) then
         n1 = size(PINPT_PY%ifilenm) 
         allocate( PINPT%ifilenm(n1)) ; PINPT%ifilenm = PINPT_PY%ifilenm
       endif
       if(allocated( PINPT_PY%title )) then
         n1 = size(PINPT_PY%title) 
         allocate( PINPT%title(n1)) ; PINPT%title = PINPT_PY%title
       endif
       if(allocated( PINPT_PY%proj_atom )) then
         n1 = size(PINPT_PY%proj_atom(:,1)) ; n2 = size(PINPT_PY%proj_atom(1,:))
         allocate( PINPT%proj_atom(n1,n2)) ; PINPT%proj_atom = PINPT_PY%proj_atom
       endif
       if(allocated( PINPT_PY%proj_natom )) then
         n1 = size(PINPT_PY%proj_natom)  
         allocate( PINPT%proj_natom(n1)) ; PINPT%proj_natom = PINPT_PY%proj_natom
       endif

    endif
    
    return
endsubroutine

function init_params_py() result(PPRAM_PY)
    type(params_py)                        PPRAM_PY

    PPRAM_PY%flag_slater_koster                      = .TRUE.
    PPRAM_PY%flag_use_overlap                        = .FALSE.
    PPRAM_PY%flag_pfile_index                        = .FALSE.
    PPRAM_PY%flag_set_param_const                    = .FALSE.
    PPRAM_PY%flag_fit_plain                          = .FALSE.
    PPRAM_PY%pfilenm                                 = 'PARAM_FIT.dat' 
    PPRAM_PY%pfileoutnm                              = 'PARAM_FIT.new.dat' 
    PPRAM_PY%l_broaden                               = 0.15d0 
    PPRAM_PY%nparam                                  = 0 
    PPRAM_PY%nparam_const                            = 0 
    PPRAM_PY%param_nsub_max                          = 1 
    PPRAM_PY%slater_koster_type                      = 1 
    PPRAM_PY%nparam_nrl                              = 0
    PPRAM_PY%nparam_nrl_free                         = 0
    PPRAM_PY%pso_nparticles                          = 50
    PPRAM_PY%pso_c1                                  = 2.0d0
    PPRAM_PY%pso_c2                                  = 2.0d0
    PPRAM_PY%pso_w                                   = 0.9d0
    PPRAM_PY%pso_max_noise_amplitude                 = 5.0d0
    PPRAM_PY%niter                                   = 0
    PPRAM_PY%pso_miter                               = 10

    if(allocated(PPRAM_PY%param))           deallocate(PPRAM_PY%param)
    if(allocated(PPRAM_PY%param_best))      deallocate(PPRAM_PY%param_best)
    if(allocated(PPRAM_PY%param_nrl))       deallocate(PPRAM_PY%param_nrl)
    if(allocated(PPRAM_PY%iparam_free))     deallocate(PPRAM_PY%iparam_free)
    if(allocated(PPRAM_PY%iparam_free_nrl)) deallocate(PPRAM_PY%iparam_free_nrl)
    if(allocated(PPRAM_PY%param_name))      deallocate(PPRAM_PY%param_name)
    if(allocated(PPRAM_PY%c_const))         deallocate(PPRAM_PY%c_const)
    if(allocated(PPRAM_PY%param_const))     deallocate(PPRAM_PY%param_const)
   !     if(allocated(PPRAM_PY%param_const_best)) deallocate(PPRAM_PY%param_const_best)
    if(allocated(PPRAM_PY%param_const_nrl)) deallocate(PPRAM_PY%param_const_nrl)
    if(allocated(PPRAM_PY%param_nsub))      deallocate(PPRAM_PY%param_nsub)

    if(allocated(PPRAM_PY%pso_cost_history)) deallocate(PPRAM_PY%pso_cost_history )
    if(allocated(PPRAM_PY%pso_cost_history_i)) deallocate(PPRAM_PY%pso_cost_history_i )
    if(allocated(PPRAM_PY%pso_pbest_history)) deallocate(PPRAM_PY%pso_pbest_history )
    if(allocated(PPRAM_PY%cost_history))     deallocate(PPRAM_PY%cost_history )
    return
endfunction

subroutine copy_params_best(PPRAM_PY, imode)
    type(params_py)    PPRAM_PY
    integer(kind=sp)         n1, imode

    if(imode .eq. 1) then ! param -> best
      PPRAM_PY%param_best = PPRAM_PY%param
    elseif(imode .eq. 2) then ! best -> param
      PPRAM_PY%param = PPRAM_PY%param_best
    endif

    return
endsubroutine

subroutine copy_params(PPRAM_PY, PPRAM, imode)
    type(params    )                       PPRAM   
    type(params_py)                        PPRAM_PY
    integer(kind=sp)                       imode, n1, n2, n3

    if(imode .eq. 1) then
       PPRAM_PY%flag_set_param_const          =      PPRAM%flag_set_param_const
       PPRAM_PY%flag_pfile_index              =      PPRAM%flag_pfile_index
       PPRAM_PY%flag_use_overlap              =      PPRAM%flag_use_overlap
       PPRAM_PY%flag_slater_koster            =      PPRAM%flag_slater_koster
       PPRAM_PY%flag_nrl_slater_koster        =      PPRAM%flag_nrl_slater_koster
       PPRAM_PY%flag_fit_plain                =      PPRAM%flag_fit_plain
       PPRAM_PY%l_broaden                     =      PPRAM%l_broaden
       PPRAM_PY%slater_koster_type            =      PPRAM%slater_koster_type
       PPRAM_PY%nparam                        =      PPRAM%nparam
       PPRAM_PY%nparam_const                  =      PPRAM%nparam_const
       PPRAM_PY%nparam_nrl                    =      PPRAM%nparam_nrl
       PPRAM_PY%nparam_nrl_free               =      PPRAM%nparam_nrl_free
       PPRAM_PY%param_nsub_max                =      PPRAM%param_nsub_max
       PPRAM_PY%nparam_free                   =      PPRAM%nparam_free
       PPRAM_PY%pfilenm                       =      PPRAM%pfilenm
       PPRAM_PY%pfileoutnm                    =      PPRAM%pfileoutnm
       PPRAM_PY%pso_nparticles                =      PPRAM%pso_nparticles
       PPRAM_PY%pso_c1                        =      PPRAM%pso_c1        
       PPRAM_PY%pso_c2                        =      PPRAM%pso_c2        
       PPRAM_PY%pso_w                         =      PPRAM%pso_w         
       PPRAM_PY%pso_max_noise_amplitude       =      PPRAM%pso_max_noise_amplitude
       PPRAM_PY%niter                         =      PPRAM%niter 
       PPRAM_PY%pso_miter                     =      PPRAM%pso_miter

       if(allocated( PPRAM_PY%param         ))         deallocate(PPRAM_PY%param         )
       if(allocated( PPRAM_PY%param_best    ))         deallocate(PPRAM_PY%param_best    )
       if(allocated( PPRAM_PY%param_nrl     ))         deallocate(PPRAM_PY%param_nrl     )
       if(allocated( PPRAM_PY%param_const   ))         deallocate(PPRAM_PY%param_const   )
       if(allocated( PPRAM_PY%param_const_nrl))        deallocate(PPRAM_PY%param_const_nrl)
       if(allocated( PPRAM_PY%param_nsub    ))         deallocate(PPRAM_PY%param_nsub    )
       if(allocated( PPRAM_PY%iparam_free   ))         deallocate(PPRAM_PY%iparam_free   )
       if(allocated( PPRAM_PY%iparam_free_nrl))        deallocate(PPRAM_PY%iparam_free_nrl)
       if(allocated( PPRAM_PY%param_name    ))         deallocate(PPRAM_PY%param_name    )
       if(allocated( PPRAM_PY%c_const       ))         deallocate(PPRAM_PY%c_const       )
       if(allocated( PPRAM_PY%pso_cost_history  ))     deallocate(PPRAM_PY%pso_cost_history )
       if(allocated( PPRAM_PY%pso_cost_history_i))     deallocate(PPRAM_PY%pso_cost_history_i )
       if(allocated( PPRAM_PY%pso_pbest_history))     deallocate(PPRAM_PY%pso_pbest_history )
       if(allocated( PPRAM_PY%cost_history  ))         deallocate(PPRAM_PY%cost_history )

       if(allocated(PPRAM%param                 )) then
          n1 = size(PPRAM%param) 
          allocate( PPRAM_PY%param(n1) )
                    PPRAM_PY%param = PPRAM%param
       endif

       if(allocated(PPRAM%param_best            )) then
          n1 = size(PPRAM%param_best)
          allocate( PPRAM_PY%param_best(n1) )
                    PPRAM_PY%param_best = PPRAM%param_best
       endif

       if(allocated(PPRAM%param_nsub            ))  then
          n1 = size(PPRAM%param_nsub)
          allocate( PPRAM_PY%param_nsub(n1)         )
                    PPRAM_PY%param_nsub = PPRAM%param_nsub
       endif
       if(allocated(PPRAM%iparam_free           )) then
          n1 = size(PPRAM%iparam_free)
          allocate(PPRAM_PY%iparam_free(n1)        )
                   PPRAM_PY%iparam_free = PPRAM%iparam_free
       endif

       if(allocated(PPRAM%iparam_free_nrl       )) then
          n1 = size(PPRAM%iparam_free_nrl)
          allocate(PPRAM_PY%iparam_free_nrl(n1)    )
                   PPRAM_PY%iparam_free_nrl = PPRAM%iparam_free_nrl
       endif

       if(allocated(PPRAM%param_name             ))  then
          n1 = size(PPRAM%param_name)
          allocate(PPRAM_PY%param_name(n1)         )
                   PPRAM_PY%param_name = PPRAM%param_name
       endif

       if(allocated(PPRAM%pso_cost_history  )) then
          n1 = size(PPRAM%pso_cost_history(:))
           allocate(PPRAM_PY%pso_cost_history(n1) ) ; PPRAM_PY%pso_cost_history =  PPRAM%pso_cost_history
       endif

       if(allocated(PPRAM%pso_cost_history_i  )) then
          n1 = size(PPRAM%pso_cost_history_i(:,1))   ; n2 = size(PPRAM%pso_cost_history_i(1,:))
           allocate(PPRAM_PY%pso_cost_history_i(n1,n2) ) 
                    PPRAM_PY%pso_cost_history_i =  PPRAM%pso_cost_history_i
       endif

       if(allocated(PPRAM%cost_history  )) then
          n1 = size(PPRAM%cost_history(:))
           allocate(PPRAM_PY%cost_history(n1) ) ; PPRAM_PY%cost_history =  PPRAM%cost_history
       endif

       if(allocated(PPRAM%pso_pbest_history  )) then
          n1 = size(PPRAM%pso_pbest_history(:,1,1))   ; n2 = size(PPRAM%pso_pbest_history(1,:,1)) ; n3 = size(PPRAM%pso_pbest_history(1,1,:))
           allocate(PPRAM_PY%pso_pbest_history(n1,n2,n3) )
                    PPRAM_PY%pso_pbest_history =  PPRAM%pso_pbest_history
       endif

       if(allocated(PPRAM%param_nrl             )) then
          n1 = size(PPRAM%param_nrl(:,1)) ; n2 = size(PPRAM%param_nrl(1,:))
          allocate(PPRAM_PY%param_nrl(n1,n2)        )
                   PPRAM_PY%param_nrl = PPRAM%param_nrl
       endif
       if(allocated(PPRAM%param_const           )) then
          n1 = size(PPRAM%param_const(:,1)) ; n2 = size(PPRAM%param_const(1,:))
          allocate(PPRAM_PY%param_const(n1,n2)      )
                   PPRAM_PY%param_const = PPRAM%param_const
       endif

!      if(allocated(PPRAM%param_const_best           )) then
!         n1 = size(PPRAM%param_const_best(:,1)) ; n2 = size(PPRAM%param_const_best(1,:))
!         allocate(PPRAM_PY%param_const_best(n1,n2)      )
!                  PPRAM_PY%param_const_best = PPRAM%param_const_best
!      endif

       if(allocated(PPRAM%c_const               )) then
          n1 = size(PPRAM%c_const(:,1)) ; n2 = size(PPRAM%c_const(1,:))
          allocate(PPRAM_PY%c_const(n1,n2)          )
                   PPRAM_PY%c_const = PPRAM%c_const
       endif
       if(allocated(PPRAM%param_const_nrl       )) then
          n1 = size(PPRAM%param_const_nrl(:,1,1)) ; n2 = size(PPRAM%param_const_nrl(1,:,1)) ; n3 = size(PPRAM%param_const_nrl(1,1,:))
          allocate(PPRAM_PY%param_const_nrl(n1,n2,n3))
                   PPRAM_PY%param_const_nrl = PPRAM%param_const_nrl
       endif

    elseif(imode .eq. 2) then
       PPRAM%flag_set_param_const          =      PPRAM_PY%flag_set_param_const
       PPRAM%flag_pfile_index              =      PPRAM_PY%flag_pfile_index
       PPRAM%flag_use_overlap              =      PPRAM_PY%flag_use_overlap
       PPRAM%flag_slater_koster            =      PPRAM_PY%flag_slater_koster
       PPRAM%flag_nrl_slater_koster        =      PPRAM_PY%flag_nrl_slater_koster
       PPRAM%flag_fit_plain                =      PPRAM_PY%flag_fit_plain
       PPRAM%l_broaden                     =      PPRAM_PY%l_broaden
       PPRAM%slater_koster_type            =      PPRAM_PY%slater_koster_type
       PPRAM%nparam                        =      PPRAM_PY%nparam
       PPRAM%nparam_const                  =      PPRAM_PY%nparam_const
       PPRAM%nparam_nrl                    =      PPRAM_PY%nparam_nrl
       PPRAM%nparam_nrl_free               =      PPRAM_PY%nparam_nrl_free
       PPRAM%param_nsub_max                =      PPRAM_PY%param_nsub_max
       PPRAM%nparam_free                   =      PPRAM_PY%nparam_free
       PPRAM%pfilenm                       =      PPRAM_PY%pfilenm
       PPRAM%pfileoutnm                    =      PPRAM_PY%pfileoutnm
       PPRAM%pso_nparticles                =      PPRAM_PY%pso_nparticles
       PPRAM%pso_c1                        =      PPRAM_PY%pso_c1      
       PPRAM%pso_c2                        =      PPRAM_PY%pso_c2      
       PPRAM%pso_w                         =      PPRAM_PY%pso_w       
       PPRAM%pso_max_noise_amplitude       =      PPRAM_PY%pso_max_noise_amplitude
       PPRAM%niter                         =      PPRAM_PY%niter 

       if(allocated( PPRAM%param         ))         deallocate(PPRAM%param         )
       if(allocated( PPRAM%param_best    ))         deallocate(PPRAM%param_best    )
       if(allocated( PPRAM%param_nrl     ))         deallocate(PPRAM%param_nrl     )
       if(allocated( PPRAM%param_const   ))         deallocate(PPRAM%param_const   )
!      if(allocated( PPRAM%param_const_best))       deallocate(PPRAM%param_const_best   )
       if(allocated( PPRAM%param_const_nrl))        deallocate(PPRAM%param_const_nrl)
       if(allocated( PPRAM%param_nsub    ))         deallocate(PPRAM%param_nsub    )
       if(allocated( PPRAM%iparam_free   ))         deallocate(PPRAM%iparam_free   )
       if(allocated( PPRAM%iparam_free_nrl))        deallocate(PPRAM%iparam_free_nrl)
       if(allocated( PPRAM%param_name    ))         deallocate(PPRAM%param_name    )
       if(allocated( PPRAM%c_const       ))         deallocate(PPRAM%c_const       )
       if(allocated( PPRAM%pso_cost_history  ))     deallocate(PPRAM%pso_cost_history )
       if(allocated( PPRAM%pso_cost_history_i))     deallocate(PPRAM%pso_cost_history_i )
       if(allocated( PPRAM%pso_pbest_history))     deallocate(PPRAM%pso_pbest_history )
       if(allocated( PPRAM%cost_history  ))         deallocate(PPRAM%cost_history )

       if(allocated(PPRAM_PY%param                 )) then
          n1 = size(PPRAM_PY%param)
          allocate( PPRAM%param(n1) )
                    PPRAM%param = PPRAM_PY%param
       endif
       if(allocated(PPRAM_PY%param_best            )) then
          n1 = size(PPRAM_PY%param_best)
          allocate( PPRAM%param_best(n1) )
                    PPRAM%param_best = PPRAM_PY%param_best
       endif


       if(allocated(PPRAM_PY%param_nsub            ))  then
          n1 = size(PPRAM_PY%param_nsub)
          allocate( PPRAM%param_nsub(n1)         )
                    PPRAM%param_nsub = PPRAM_PY%param_nsub
       endif
       if(allocated(PPRAM_PY%iparam_free           )) then
          n1 = size(PPRAM_PY%iparam_free)
          allocate(PPRAM%iparam_free(n1)        )
                   PPRAM%iparam_free = PPRAM_PY%iparam_free
       endif

       if(allocated(PPRAM_PY%iparam_free_nrl       )) then
          n1 = size(PPRAM_PY%iparam_free_nrl)
          allocate(PPRAM%iparam_free_nrl(n1)    )
                   PPRAM%iparam_free_nrl = PPRAM_PY%iparam_free_nrl
       endif

       if(allocated(PPRAM_PY%param_name             ))  then
          n1 = size(PPRAM_PY%param_name)
          allocate(PPRAM%param_name(n1)         )
                   PPRAM%param_name = PPRAM_PY%param_name
       endif

       if(allocated(PPRAM_PY%pso_cost_history  )) then
          n1 = size(PPRAM_PY%pso_cost_history(:))
           allocate(PPRAM%pso_cost_history(n1) ) ; PPRAM%pso_cost_history =  PPRAM_PY%pso_cost_history
       endif

       if(allocated(PPRAM_PY%pso_cost_history_i  )) then
          n1 = size(PPRAM_PY%pso_cost_history_i(:,1)) ; n2 = size(PPRAM_PY%pso_cost_history_i(1,:))
           allocate(PPRAM%pso_cost_history_i(n1,n2) ) ; PPRAM%pso_cost_history_i =  PPRAM_PY%pso_cost_history_i
       endif

       if(allocated(PPRAM_PY%cost_history  )) then
          n1 = size(PPRAM_PY%cost_history(:))
           allocate(PPRAM%cost_history(n1) ) ; PPRAM%cost_history =  PPRAM_PY%cost_history
       endif

       if(allocated(PPRAM_PY%pso_pbest_history  )) then
          n1 = size(PPRAM_PY%pso_pbest_history(:,1,1))   ; n2 = size(PPRAM_PY%pso_pbest_history(1,:,1)) ; n3 = size(PPRAM_PY%pso_pbest_history(1,1,:))
           allocate(PPRAM%pso_pbest_history(n1,n2,n3) )
                    PPRAM%pso_pbest_history =  PPRAM_PY%pso_pbest_history
       endif

       if(allocated(PPRAM_PY%param_nrl             )) then       
          n1 = size(PPRAM_PY%param_nrl(:,1)) ; n2 = size(PPRAM_PY%param_nrl(1,:))
          allocate(PPRAM%param_nrl(n1,n2)        )
                   PPRAM%param_nrl = PPRAM_PY%param_nrl
       endif
       if(allocated(PPRAM_PY%param_const           )) then
          n1 = size(PPRAM_PY%param_const(:,1)) ; n2 = size(PPRAM_PY%param_const(1,:))
          allocate(PPRAM%param_const(n1,n2)      )
                   PPRAM%param_const = PPRAM_PY%param_const
       endif
!      if(allocated(PPRAM_PY%param_const_best           )) then
!         n1 = size(PPRAM_PY%param_const_best(:,1)) ; n2 = size(PPRAM_PY%param_const_best(1,:))
!         allocate(PPRAM%param_const_best(n1,n2)      )
!                  PPRAM%param_const_best = PPRAM_PY%param_const_best
!      endif

       if(allocated(PPRAM_PY%c_const               )) then
          n1 = size(PPRAM_PY%c_const(:,1)) ; n2 = size(PPRAM_PY%c_const(1,:))
          allocate(PPRAM%c_const(n1,n2)          )
                   PPRAM%c_const = PPRAM_PY%c_const
       endif
       if(allocated(PPRAM_PY%param_const_nrl       )) then
          n1 = size(PPRAM_PY%param_const_nrl(:,1,1)) ; n2 = size(PPRAM_PY%param_const_nrl(1,:,1)) ; n3 = size(PPRAM_PY%param_const_nrl(1,1,:))
          allocate(PPRAM%param_const_nrl(n1,n2,n3))
                   PPRAM%param_const_nrl = PPRAM_PY%param_const_nrl
       endif

    endif      

endsubroutine

function init_kpoints_py() result(PKPTS_PY)
    type(kpoints_py)                       PKPTS_PY
   
    PKPTS_PY%kfilenm                                  = 'KPOINTS_BAND'
    PKPTS_PY%kline_type                               = 'vasp'
    PKPTS_PY%kunit                                    = 'A' 
    PKPTS_PY%kreduce                                  = 1
    PKPTS_PY%n_ndiv                                   = 1
    PKPTS_PY%idiv_mode                                = 1 

    if(allocated(PKPTS_PY%ndiv       ))  deallocate(PKPTS_PY%ndiv       )
    if(allocated(PKPTS_PY%kpoint     ))  deallocate(PKPTS_PY%kpoint     )
    if(allocated(PKPTS_PY%kline      ))  deallocate(PKPTS_PY%kline      )
    if(allocated(PKPTS_PY%kpoint_reci))  deallocate(PKPTS_PY%kpoint_reci)
    if(allocated(PKPTS_PY%k_name     ))  deallocate(PKPTS_PY%k_name     )
    if(allocated(PKPTS_PY%k_name2    ))  deallocate(PKPTS_PY%k_name2    )
    if(allocated(PKPTS_PY%k_name_index)) deallocate(PKPTS_PY%k_name_index)

    if(allocated(PKPTS_PY%kdist      ))  deallocate(PKPTS_PY%kdist      ) ! only for python module (not to be copied to PKPTS)

    return
endfunction

subroutine copy_kpoints(PKPTS_PY, PKPTS, imode)
    type(kpoints   )                       PKPTS
    type(kpoints_py)                       PKPTS_PY
    integer(kind=sp)                       imode, n1, n2, n3

    if(imode .eq. 1) then
       PKPTS_PY%flag_klinemode                =      PKPTS%flag_klinemode  
       PKPTS_PY%flag_kgridmode                =      PKPTS%flag_kgridmode
       PKPTS_PY%flag_gamma                    =      PKPTS%flag_gamma
       PKPTS_PY%flag_reciprocal               =      PKPTS%flag_reciprocal
       PKPTS_PY%flag_cartesianK               =      PKPTS%flag_cartesianK
       PKPTS_PY%k_shift                       =      PKPTS%k_shift   
       PKPTS_PY%mysystem                      =      PKPTS%mysystem
       PKPTS_PY%nkpoint                       =      PKPTS%nkpoint
       PKPTS_PY%nline                         =      PKPTS%nline
       PKPTS_PY%n_ndiv                        =      PKPTS%n_ndiv
       PKPTS_PY%idiv_mode                     =      PKPTS%idiv_mode
       PKPTS_PY%kreduce                       =      PKPTS%kreduce
       PKPTS_PY%kfilenm                       =      PKPTS%kfilenm
       PKPTS_PY%ribbon_kfilenm                =      PKPTS%ribbon_kfilenm
       PKPTS_PY%kline_type                    =      PKPTS%kline_type
       PKPTS_PY%kunit                         =      PKPTS%kunit

       if(allocated( PKPTS_PY%ndiv              ))      deallocate(PKPTS_PY%ndiv              )
       if(allocated( PKPTS_PY%k_name            ))      deallocate(PKPTS_PY%k_name            )
       if(allocated( PKPTS_PY%k_name2           ))      deallocate(PKPTS_PY%k_name2           )
       if(allocated( PKPTS_PY%k_name_index      ))      deallocate(PKPTS_PY%k_name_index      )
       if(allocated( PKPTS_PY%kpoint            ))      deallocate(PKPTS_PY%kpoint            )
       if(allocated( PKPTS_PY%kline             ))      deallocate(PKPTS_PY%kline             )
       if(allocated( PKPTS_PY%kpoint_reci       ))      deallocate(PKPTS_PY%kpoint_reci       )

       if(allocated( PKPTS%ndiv              )) then
           n1 = size(PKPTS%ndiv)
            allocate(PKPTS_PY%ndiv(n1)          )
                     PKPTS_PY%ndiv = PKPTS%ndiv
       endif

       if(allocated( PKPTS%k_name            )) then
           n1 = size(PKPTS%k_name)
            allocate(PKPTS_PY%k_name(n1)        )
                     PKPTS_PY%k_name = PKPTS%k_name
       endif

       if(allocated( PKPTS%k_name2           )) then
           n1 = size(PKPTS%k_name2)
            allocate(PKPTS_PY%k_name2(n1)        )
                     PKPTS_PY%k_name2 = PKPTS%k_name2
       endif

       if(allocated( PKPTS%k_name_index      )) then
           n1 = size(PKPTS%k_name_index)
            allocate(PKPTS_PY%k_name_index(n1)        )
                     PKPTS_PY%k_name_index = PKPTS%k_name_index
       endif

       if(allocated( PKPTS%kpoint            )) then
           n1 = size(PKPTS%kpoint(:,1)) ; n2 = size(PKPTS%kpoint(1,:))
            allocate(PKPTS_PY%kpoint(n1,n2)     )
                     PKPTS_PY%kpoint = PKPTS%kpoint
       endif

       if(allocated( PKPTS%kline             )) then
           n1 = size(PKPTS%kline(:,1)) ; n2 = size(PKPTS%kline(1,:))
            allocate(PKPTS_PY%kline(n1,n2)     )
                     PKPTS_PY%kline = PKPTS%kline
       endif

       if(allocated( PKPTS%kpoint_reci       )) then
           n1 = size(PKPTS%kpoint_reci(:,1)) ; n2 = size(PKPTS%kpoint_reci(1,:))
            allocate(PKPTS_PY%kpoint_reci(n1,n2)     )
                     PKPTS_PY%kpoint_reci = PKPTS%kpoint_reci
       endif

    elseif(imode .eq. 2) then
       PKPTS%flag_klinemode                =      PKPTS_PY%flag_klinemode
       PKPTS%flag_kgridmode                =      PKPTS_PY%flag_kgridmode
       PKPTS%flag_gamma                    =      PKPTS_PY%flag_gamma
       PKPTS%flag_reciprocal               =      PKPTS_PY%flag_reciprocal
       PKPTS%flag_cartesianK               =      PKPTS_PY%flag_cartesianK
       PKPTS%k_shift                       =      PKPTS_PY%k_shift
       PKPTS%mysystem                      =      PKPTS_PY%mysystem
       PKPTS%nkpoint                       =      PKPTS_PY%nkpoint
       PKPTS%nline                         =      PKPTS_PY%nline
       PKPTS%n_ndiv                        =      PKPTS_PY%n_ndiv
       PKPTS%idiv_mode                     =      PKPTS_PY%idiv_mode
       PKPTS%kreduce                       =      PKPTS_PY%kreduce
       PKPTS%kfilenm                       =      PKPTS_PY%kfilenm
       PKPTS%ribbon_kfilenm                =      PKPTS_PY%ribbon_kfilenm
       PKPTS%kline_type                    =      PKPTS_PY%kline_type
       PKPTS%kunit                         =      PKPTS_PY%kunit

       if(allocated( PKPTS%ndiv              ))      deallocate(PKPTS%ndiv              )
       if(allocated( PKPTS%k_name            ))      deallocate(PKPTS%k_name            )
       if(allocated( PKPTS%k_name2           ))      deallocate(PKPTS%k_name2           )
       if(allocated( PKPTS%k_name_index      ))      deallocate(PKPTS%k_name_index      )
       if(allocated( PKPTS%kpoint            ))      deallocate(PKPTS%kpoint            )
       if(allocated( PKPTS%kline             ))      deallocate(PKPTS%kline             )
       if(allocated( PKPTS%kpoint_reci       ))      deallocate(PKPTS%kpoint_reci       )

       if(allocated( PKPTS_py%ndiv              )) then
           n1 = size(PKPTS_py%ndiv)
            allocate(PKPTS%ndiv(n1)          )
                     PKPTS%ndiv = PKPTS_PY%ndiv
       endif

       if(allocated( PKPTS_PY%k_name            )) then
           n1 = size(PKPTS_PY%k_name)
            allocate(PKPTS%k_name(n1)        )
                     PKPTS%k_name = PKPTS_PY%k_name
       endif

       if(allocated( PKPTS_PY%k_name2           )) then
           n1 = size(PKPTS_PY%k_name2)
            allocate(PKPTS%k_name2(n1)        )
                     PKPTS%k_name2 = PKPTS_PY%k_name2
       endif

       if(allocated( PKPTS_PY%k_name_index            )) then
           n1 = size(PKPTS_PY%k_name_index)
            allocate(PKPTS%k_name_index(n1)        )
                     PKPTS%k_name_index = PKPTS_PY%k_name_index
       endif

       if(allocated( PKPTS_PY%kpoint            )) then
           n1 = size(PKPTS_PY%kpoint(:,1)) ; n2 = size(PKPTS_PY%kpoint(1,:))
            allocate(PKPTS%kpoint(n1,n2)     )
                     PKPTS%kpoint = PKPTS_PY%kpoint
       endif

       if(allocated( PKPTS_PY%kline             )) then
           n1 = size(PKPTS_PY%kline(:,1)) ; n2 = size(PKPTS_PY%kline(1,:))
            allocate(PKPTS%kline(n1,n2)     )
                     PKPTS%kline = PKPTS_PY%kline
       endif

       if(allocated( PKPTS_PY%kpoint_reci       )) then
           n1 = size(PKPTS_PY%kpoint_reci(:,1)) ; n2 = size(PKPTS_PY%kpoint_reci(1,:))
            allocate(PKPTS%kpoint_reci(n1,n2)     )
                     PKPTS%kpoint_reci = PKPTS_PY%kpoint_reci
       endif

    endif

    return
endsubroutine

function init_weight_py() result(PWGHT_PY)
    type(weight_py)                        PWGHT_PY

    PWGHT_PY%flag_weight_default                      = .TRUE.
    PWGHT_PY%flag_weight_orb                          = .FALSE.
    PWGHT_PY%itarget_e_start                          = 1 ! default
    PWGHT_PY%efile_ef                                 = 0d0
    PWGHT_PY%efile_type                               = 'user'
    PWGHT_PY%efilenmu                                 = ' '
    PWGHT_PY%efilenmd                                 = ' '
    PWGHT_PY%read_energy_column_index                 = 2 
    PWGHT_PY%read_energy_column_index_dn              = 2 
    PWGHT_PY%ie_cutoff                                = 0

    if(allocated(PWGHT_PY%WT))   deallocate(PWGHT_PY%WT)

    return
endfunction

subroutine copy_weight(PWGHT_PY, PWGHT, imode)
    type(weight    )                       PWGHT
    type(weight_py )                       PWGHT_PY
    integer(kind=sp)                       imode, n1, n2, n3

    if(imode .eq. 1) then
       PWGHT_PY%flag_weight_default           =      PWGHT%flag_weight_default
       PWGHT_PY%flag_weight_orb               =      PWGHT%flag_weight_orb
       PWGHT_PY%efile_ef                      =      PWGHT%efile_ef
       PWGHT_PY%mysystem                      =      PWGHT%mysystem
       PWGHT_PY%itarget_e_start               =      PWGHT%itarget_e_start
       PWGHT_PY%read_energy_column_index      =      PWGHT%read_energy_column_index
       PWGHT_PY%read_energy_column_index_dn   =      PWGHT%read_energy_column_index_dn
       PWGHT_PY%ie_cutoff                     =      PWGHT%ie_cutoff
       PWGHT_PY%efilenmu                      =      PWGHT%efilenmu
       PWGHT_PY%efilenmd                      =      PWGHT%efilenmd
       PWGHT_PY%efile_type                    =      PWGHT%efile_type

       if(allocated( PWGHT_PY%WT        ))      deallocate(PWGHT_PY%WT)

       if(allocated( PWGHT%WT           ))  then
           n1 = size(PWGHT%WT(:,1)) ; n2 = size(PWGHT%WT(1,:))
            allocate(PWGHT_PY%WT(n1,n2)    )
                     PWGHT_PY%WT =  PWGHT%WT
       endif

    elseif(imode .eq. 2) then

       PWGHT%flag_weight_default           =      PWGHT_PY%flag_weight_default
       PWGHT%flag_weight_orb               =      PWGHT_PY%flag_weight_orb
       PWGHT%efile_ef                      =      PWGHT_PY%efile_ef
       PWGHT%mysystem                      =      PWGHT_PY%mysystem
       PWGHT%itarget_e_start               =      PWGHT_PY%itarget_e_start
       PWGHT%read_energy_column_index      =      PWGHT_PY%read_energy_column_index
       PWGHT%read_energy_column_index_dn   =      PWGHT_PY%read_energy_column_index_dn
       PWGHT%ie_cutoff                     =      PWGHT_PY%ie_cutoff
       PWGHT%efilenmu                      =      PWGHT_PY%efilenmu
       PWGHT%efilenmd                      =      PWGHT_PY%efilenmd
       PWGHT%efile_type                    =      PWGHT_PY%efile_type

       if(allocated( PWGHT%WT        ))      deallocate(PWGHT%WT)

       if(allocated( PWGHT_PY%WT           ))  then
           n1 = size(PWGHT_PY%WT(:,1)) ; n2 = size(PWGHT_PY%WT(1,:))
            allocate(PWGHT%WT(n1,n2)    )
                     PWGHT%WT =  PWGHT_PY%WT
       endif

    endif

endsubroutine

function init_poscar_py() result(PGEOM_PY)
    type(poscar_py)                        PGEOM_PY

    PGEOM_PY%title                                    = '' ! default
    PGEOM_PY%gfilenm                                  = 'POSCAR-TB' !default

    if(allocated(PGEOM_PY%nelect               )) deallocate(PGEOM_PY%nelect)
    if(allocated(PGEOM_PY%n_orbital            )) deallocate(PGEOM_PY%n_orbital            )
    if(allocated(PGEOM_PY%orb_index            )) deallocate(PGEOM_PY%orb_index            )
    if(allocated(PGEOM_PY%spec                 )) deallocate(PGEOM_PY%spec                 )
    if(allocated(PGEOM_PY%c_orbital            )) deallocate(PGEOM_PY%c_orbital            )
    if(allocated(PGEOM_PY%c_spec               )) deallocate(PGEOM_PY%c_spec               )

    return
endfunction

subroutine copy_poscar(PGEOM_PY, PGEOM, imode)
    type(poscar    )                       PGEOM
    type(poscar_py )                       PGEOM_PY
    integer(kind=sp)                       imode, n1, n2, n3

    if(imode .eq. 1) then
       PGEOM_PY%flag_selective                =      PGEOM%flag_selective
       PGEOM_PY%flag_direct                   =      PGEOM%flag_direct
       PGEOM_PY%flag_cartesian                =      PGEOM%flag_cartesian
       PGEOM_PY%mysystem                      =      PGEOM%mysystem
       PGEOM_PY%neig                          =      PGEOM%neig
       PGEOM_PY%neig_total                    =      PGEOM%neig_total
       PGEOM_PY%neig_target                   =      PGEOM%neig_target
       PGEOM_PY%nbasis                        =      PGEOM%nbasis
       PGEOM_PY%neig_eff                      =      PGEOM%neig_eff
       PGEOM_PY%init_erange                   =      PGEOM%init_erange
       PGEOM_PY%fina_erange                   =      PGEOM%fina_erange
       PGEOM_PY%nband                         =      PGEOM%nband
       PGEOM_PY%n_spec                        =      PGEOM%n_spec
       PGEOM_PY%n_atom                        =      PGEOM%n_atom
       PGEOM_PY%max_orb                       =      PGEOM%max_orb
       PGEOM_PY%title                         =      PGEOM%title
       PGEOM_PY%gfilenm                       =      PGEOM%gfilenm
       PGEOM_PY%system_name                   =      PGEOM%system_name
       PGEOM_PY%a_scale                       =      PGEOM%a_scale
       PGEOM_PY%a_latt                        =      PGEOM%a_latt
       PGEOM_PY%b_latt                        =      PGEOM%b_latt

       if(allocated( PGEOM_PY%nelect           ))   deallocate(PGEOM_PY%nelect)
       if(allocated( PGEOM_PY%n_orbital        ))   deallocate(PGEOM_PY%n_orbital)
       if(allocated( PGEOM_PY%orb_index        ))   deallocate(PGEOM_PY%orb_index)
       if(allocated( PGEOM_PY%spec             ))   deallocate(PGEOM_PY%spec     )
       if(allocated( PGEOM_PY%c_orbital        ))   deallocate(PGEOM_PY%c_orbital)
       if(allocated( PGEOM_PY%c_spec           ))   deallocate(PGEOM_PY%c_spec   )

       if(allocated( PGEOM%nelect           ))  then
           n1 = size(PGEOM%nelect)
            allocate(PGEOM_PY%nelect(n1)    )
                     PGEOM_PY%nelect = PGEOM%nelect
       endif

       if(allocated( PGEOM%n_orbital        ))  then
           n1 = size(PGEOM%n_orbital)
            allocate(PGEOM_PY%n_orbital(n1) )
                     PGEOM_PY%n_orbital = PGEOM%n_orbital
       endif

       if(allocated( PGEOM%orb_index        ))  then
           n1 = size(PGEOM%orb_index )
            allocate(PGEOM_PY%orb_index(n1) )
                     PGEOM_PY%orb_index = PGEOM%orb_index
       endif

       if(allocated( PGEOM%spec             ))  then
           n1 = size(PGEOM%spec      )
            allocate(PGEOM_PY%spec(n1) )
                     PGEOM_PY%spec      = PGEOM%spec      
       endif

       if(allocated( PGEOM%c_orbital             )) then
           n1 = size(PGEOM%c_orbital(:,1)) ; n2 = size(PGEOM%c_orbital(1,:))
            allocate(PGEOM_PY%c_orbital(n1,n2)          )
                     PGEOM_PY%c_orbital = PGEOM%c_orbital
       endif

       if(allocated( PGEOM%c_spec           ))  then
           n1 = size(PGEOM%c_spec    )
            allocate(PGEOM_PY%c_spec(n1) )
                     PGEOM_PY%c_spec    = PGEOM%c_spec   
       endif


    elseif(imode .eq. 2) then
       PGEOM%flag_selective                =      PGEOM_PY%flag_selective
       PGEOM%flag_direct                   =      PGEOM_PY%flag_direct
       PGEOM%flag_cartesian                =      PGEOM_PY%flag_cartesian
       PGEOM%mysystem                      =      PGEOM_PY%mysystem
       PGEOM%neig                          =      PGEOM_PY%neig
       PGEOM%neig_total                    =      PGEOM_PY%neig_total
       PGEOM%neig_target                   =      PGEOM_PY%neig_target
       PGEOM%nbasis                        =      PGEOM_PY%nbasis
       PGEOM%neig_eff                      =      PGEOM_PY%neig_eff
       PGEOM%init_erange                   =      PGEOM_PY%init_erange
       PGEOM%fina_erange                   =      PGEOM_PY%fina_erange
       PGEOM%nband                         =      PGEOM_PY%nband
       PGEOM%n_spec                        =      PGEOM_PY%n_spec
       PGEOM%n_atom                        =      PGEOM_PY%n_atom
       PGEOM%max_orb                       =      PGEOM_PY%max_orb
       PGEOM%title                         =      PGEOM_PY%title
       PGEOM%gfilenm                       =      PGEOM_PY%gfilenm
       PGEOM%system_name                   =      PGEOM_PY%system_name
       PGEOM%a_scale                       =      PGEOM_PY%a_scale
       PGEOM%a_latt                        =      PGEOM_PY%a_latt
       PGEOM%b_latt                        =      PGEOM_PY%b_latt

       if(allocated( PGEOM%nelect           ))   deallocate(PGEOM%nelect)
       if(allocated( PGEOM%n_orbital        ))   deallocate(PGEOM%n_orbital)
       if(allocated( PGEOM%orb_index        ))   deallocate(PGEOM%orb_index)
       if(allocated( PGEOM%spec             ))   deallocate(PGEOM%spec     )
       if(allocated( PGEOM%c_orbital        ))   deallocate(PGEOM%c_orbital)
       if(allocated( PGEOM%c_spec           ))   deallocate(PGEOM%c_spec   )

       if(allocated( PGEOM_PY%nelect           ))  then
           n1 = size(PGEOM_PY%nelect)
            allocate(PGEOM%nelect(n1)    )
                     PGEOM%nelect = PGEOM_PY%nelect
       endif

       if(allocated( PGEOM_PY%n_orbital        ))  then
           n1 = size(PGEOM_PY%n_orbital)
            allocate(PGEOM%n_orbital(n1) )
                     PGEOM%n_orbital = PGEOM_PY%n_orbital
       endif

       if(allocated( PGEOM_PY%orb_index        ))  then
           n1 = size(PGEOM_PY%orb_index )
            allocate(PGEOM%orb_index(n1) )
                     PGEOM%orb_index = PGEOM_PY%orb_index
       endif

       if(allocated( PGEOM_PY%spec             ))  then
           n1 = size(PGEOM_PY%spec      )
            allocate(PGEOM%spec(n1) )
                     PGEOM%spec      = PGEOM_PY%spec     
       endif

       if(allocated( PGEOM_PY%c_orbital             )) then
           n1 = size(PGEOM_PY%c_orbital(:,1)) ; n2 = size(PGEOM_PY%c_orbital(1,:))
            allocate(PGEOM%c_orbital(n1,n2)          )
                     PGEOM%c_orbital = PGEOM_PY%c_orbital
       endif

       if(allocated( PGEOM_PY%c_spec           ))  then
           n1 = size(PGEOM_PY%c_spec    )
            allocate(PGEOM%c_spec(n1) )
                     PGEOM%c_spec    = PGEOM_PY%c_spec   
       endif

    endif

    return
endsubroutine

function init_hopping_py() result(NN_TABLE_PY)
    type(hopping_py)                       NN_TABLE_PY

    NN_TABLE_PY%nnfilenmo = 'hopping.dat' ! default
    NN_TABLE_PY%nnfilenm  = 'hopping.dat'
    NN_TABLE_PY%onsite_tolerance = onsite_tolerance ! default defined in parameters.f90
    NN_TABLE_PY%flag_efield = .false.

    if(allocated(NN_TABLE_PY%flag_site_cindex      ))    deallocate(NN_TABLE_PY%flag_site_cindex      )
    if(allocated(NN_TABLE_PY%i_coord               ))    deallocate(NN_TABLE_PY%i_coord               )
    if(allocated(NN_TABLE_PY%j_coord               ))    deallocate(NN_TABLE_PY%j_coord               )
    if(allocated(NN_TABLE_PY%Rij                   ))    deallocate(NN_TABLE_PY%Rij                   )
    if(allocated(NN_TABLE_PY%R                     ))    deallocate(NN_TABLE_PY%R                     )
    if(allocated(NN_TABLE_PY%Dij                   ))    deallocate(NN_TABLE_PY%Dij                   )
    if(allocated(NN_TABLE_PY%Dij0                  ))    deallocate(NN_TABLE_PY%Dij0                  )
    if(allocated(NN_TABLE_PY%i_sign                ))    deallocate(NN_TABLE_PY%i_sign                )
    if(allocated(NN_TABLE_PY%j_sign                ))    deallocate(NN_TABLE_PY%j_sign                )
    if(allocated(NN_TABLE_PY%R_nn                  ))    deallocate(NN_TABLE_PY%R_nn                  )
    if(allocated(NN_TABLE_PY%R0_nn                 ))    deallocate(NN_TABLE_PY%R0_nn                 )
    if(allocated(NN_TABLE_PY%local_charge          ))    deallocate(NN_TABLE_PY%local_charge          )
    if(allocated(NN_TABLE_PY%local_moment          ))    deallocate(NN_TABLE_PY%local_moment          )
    if(allocated(NN_TABLE_PY%local_moment_rot      ))    deallocate(NN_TABLE_PY%local_moment_rot      )
    if(allocated(NN_TABLE_PY%i_atom                ))    deallocate(NN_TABLE_PY%i_atom                )
    if(allocated(NN_TABLE_PY%j_atom                ))    deallocate(NN_TABLE_PY%j_atom                )
    if(allocated(NN_TABLE_PY%i_matrix              ))    deallocate(NN_TABLE_PY%i_matrix              )
    if(allocated(NN_TABLE_PY%j_matrix              ))    deallocate(NN_TABLE_PY%j_matrix              )
    if(allocated(NN_TABLE_PY%n_class               ))    deallocate(NN_TABLE_PY%n_class               )
    if(allocated(NN_TABLE_PY%sk_index_set          ))    deallocate(NN_TABLE_PY%sk_index_set          )
    if(allocated(NN_TABLE_PY%cc_index_set          ))    deallocate(NN_TABLE_PY%cc_index_set          )
    if(allocated(NN_TABLE_PY%n_nn                  ))    deallocate(NN_TABLE_PY%n_nn                  )
    if(allocated(NN_TABLE_PY%j_nn                  ))    deallocate(NN_TABLE_PY%j_nn                  )
    if(allocated(NN_TABLE_PY%l_onsite_param_index  ))    deallocate(NN_TABLE_PY%l_onsite_param_index  )
    if(allocated(NN_TABLE_PY%stoner_I_param_index  ))    deallocate(NN_TABLE_PY%stoner_I_param_index  )
    if(allocated(NN_TABLE_PY%local_U_param_index   ))    deallocate(NN_TABLE_PY%local_U_param_index   )
    if(allocated(NN_TABLE_PY%plus_U_param_index    ))    deallocate(NN_TABLE_PY%plus_U_param_index    )
    if(allocated(NN_TABLE_PY%soc_param_index       ))    deallocate(NN_TABLE_PY%soc_param_index       )
    if(allocated(NN_TABLE_PY%ci_orb                ))    deallocate(NN_TABLE_PY%ci_orb                )
    if(allocated(NN_TABLE_PY%cj_orb                ))    deallocate(NN_TABLE_PY%cj_orb                )
    if(allocated(NN_TABLE_PY%p_class               ))    deallocate(NN_TABLE_PY%p_class               )
    if(allocated(NN_TABLE_PY%site_cindex           ))    deallocate(NN_TABLE_PY%site_cindex           )

    return
endfunction

subroutine copy_hopping(NN_TABLE_PY, NN_TABLE, imode)
    type(hopping   )                       NN_TABLE
    type(hopping_py)                       NN_TABLE_PY
    integer(kind=sp)                       imode, n1, n2, n3

    if(imode .eq. 1) then
       NN_TABLE_PY%flag_efield                =      NN_TABLE%flag_efield
       NN_TABLE_PY%flag_efield_frac           =      NN_TABLE%flag_efield_frac
       NN_TABLE_PY%flag_efield_cart           =      NN_TABLE%flag_efield_cart
       NN_TABLE_PY%onsite_tolerance           =      NN_TABLE%onsite_tolerance
       NN_TABLE_PY%efield                     =      NN_TABLE%efield
       NN_TABLE_PY%efield_origin              =      NN_TABLE%efield_origin
       NN_TABLE_PY%efield_origin_cart         =      NN_TABLE%efield_origin_cart
       NN_TABLE_PY%mysystem                   =      NN_TABLE%mysystem
       NN_TABLE_PY%n_neighbor                 =      NN_TABLE%n_neighbor
       NN_TABLE_PY%max_nn_pair                =      NN_TABLE%max_nn_pair
       NN_TABLE_PY%nnfilenm                   =      NN_TABLE%nnfilenm
       NN_TABLE_PY%nnfilenmo                  =      NN_TABLE%nnfilenmo
    
       if(allocated( NN_TABLE_PY%flag_site_cindex        ))   deallocate( NN_TABLE_PY%flag_site_cindex)    
       if(allocated( NN_TABLE_PY%Dij                     ))   deallocate( NN_TABLE_PY%Dij)                 
       if(allocated( NN_TABLE_PY%Dij0                    ))   deallocate( NN_TABLE_PY%Dij0)                
       if(allocated( NN_TABLE_PY%i_sign                  ))   deallocate( NN_TABLE_PY%i_sign)              
       if(allocated( NN_TABLE_PY%j_sign                  ))   deallocate( NN_TABLE_PY%j_sign)              
       if(allocated( NN_TABLE_PY%local_charge            ))   deallocate( NN_TABLE_PY%local_charge)        
       if(allocated( NN_TABLE_PY%i_atom                  ))   deallocate( NN_TABLE_PY%i_atom)              
       if(allocated( NN_TABLE_PY%j_atom                  ))   deallocate( NN_TABLE_PY%j_atom)              
       if(allocated( NN_TABLE_PY%i_matrix                ))   deallocate( NN_TABLE_PY%i_matrix)            
       if(allocated( NN_TABLE_PY%j_matrix                ))   deallocate( NN_TABLE_PY%j_matrix)            
       if(allocated( NN_TABLE_PY%n_class                 ))   deallocate( NN_TABLE_PY%n_class)             
       if(allocated( NN_TABLE_PY%n_nn                    ))   deallocate( NN_TABLE_PY%n_nn)                
       if(allocated( NN_TABLE_PY%l_onsite_param_index    ))   deallocate( NN_TABLE_PY%l_onsite_param_index)
       if(allocated( NN_TABLE_PY%stoner_I_param_index    ))   deallocate( NN_TABLE_PY%stoner_I_param_index)
       if(allocated( NN_TABLE_PY%local_U_param_index     ))   deallocate( NN_TABLE_PY%local_U_param_index) 
       if(allocated( NN_TABLE_PY%plus_U_param_index      ))   deallocate( NN_TABLE_PY%plus_U_param_index)  
       if(allocated( NN_TABLE_PY%soc_param_index         ))   deallocate( NN_TABLE_PY%soc_param_index)     
       if(allocated( NN_TABLE_PY%ci_orb                  ))   deallocate( NN_TABLE_PY%ci_orb)              
       if(allocated( NN_TABLE_PY%cj_orb                  ))   deallocate( NN_TABLE_PY%cj_orb)              
       if(allocated( NN_TABLE_PY%p_class                 ))   deallocate( NN_TABLE_PY%p_class)             
       if(allocated( NN_TABLE_PY%site_cindex             ))   deallocate( NN_TABLE_PY%site_cindex)         


       if(allocated( NN_TABLE%flag_site_cindex        ))   then
          n1 = size( NN_TABLE%flag_site_cindex        )  
           allocate( NN_TABLE_PY%flag_site_cindex(n1)     )  ;  NN_TABLE_PY%flag_site_cindex          = NN_TABLE%flag_site_cindex
       endif

       if(allocated( NN_TABLE%Dij                     ))   then
          n1 = size( NN_TABLE%Dij                     )  
           allocate( NN_TABLE_PY%Dij(n1)                  )  ;  NN_TABLE_PY%Dij                       = NN_TABLE%Dij
       endif

       if(allocated( NN_TABLE%Dij0                    ))   then
          n1 = size( NN_TABLE%Dij0                    )  
           allocate( NN_TABLE_PY%Dij0(n1)                 )  ;  NN_TABLE_PY%Dij0                      = NN_TABLE%Dij0
       endif

       if(allocated( NN_TABLE%i_sign                  ))   then
          n1 = size( NN_TABLE%i_sign                  )  
           allocate( NN_TABLE_PY%i_sign(n1)               )  ;  NN_TABLE_PY%i_sign                    = NN_TABLE%i_sign
       endif

       if(allocated( NN_TABLE%j_sign                  ))   then
          n1 = size( NN_TABLE%j_sign                  )  
           allocate( NN_TABLE_PY%j_sign(n1)               )  ;  NN_TABLE_PY%j_sign                    = NN_TABLE%j_sign
       endif

       if(allocated( NN_TABLE%local_charge            ))   then
          n1 = size( NN_TABLE%local_charge            )  
           allocate( NN_TABLE_PY%local_charge(n1)         )  ;  NN_TABLE_PY%local_charge              = NN_TABLE%local_charge
       endif

       if(allocated( NN_TABLE%i_atom                  ))   then
          n1 = size( NN_TABLE%i_atom                  )  
           allocate( NN_TABLE_PY%i_atom(n1)               )  ;  NN_TABLE_PY%i_atom                    = NN_TABLE%i_atom
       endif

       if(allocated( NN_TABLE%j_atom                  ))   then
          n1 = size( NN_TABLE%j_atom                  )  
           allocate( NN_TABLE_PY%j_atom(n1)               )  ;  NN_TABLE_PY%j_atom                    = NN_TABLE%j_atom
       endif

       if(allocated( NN_TABLE%i_matrix                ))   then
          n1 = size( NN_TABLE%i_matrix                )  
           allocate( NN_TABLE_PY%i_matrix(n1)             )  ;  NN_TABLE_PY%i_matrix                  = NN_TABLE%i_matrix
       endif

       if(allocated( NN_TABLE%j_matrix                ))   then
          n1 = size( NN_TABLE%j_matrix                )  
           allocate( NN_TABLE_PY%j_matrix(n1)             )  ;  NN_TABLE_PY%j_matrix                  = NN_TABLE%j_matrix
       endif

       if(allocated( NN_TABLE%n_class                 ))   then
          n1 = size( NN_TABLE%n_class                 )  
           allocate( NN_TABLE_PY%n_class(n1)              )  ;  NN_TABLE_PY%n_class                   = NN_TABLE%n_class
       endif

       if(allocated( NN_TABLE%n_nn                    ))   then
          n1 = size( NN_TABLE%n_nn                    )  
           allocate( NN_TABLE_PY%n_nn(n1)                 )  ;  NN_TABLE_PY%n_nn                      = NN_TABLE%n_nn
       endif

       if(allocated( NN_TABLE%l_onsite_param_index    ))   then
          n1 = size( NN_TABLE%l_onsite_param_index    )  
           allocate( NN_TABLE_PY%l_onsite_param_index(n1) )  ;  NN_TABLE_PY%l_onsite_param_index      = NN_TABLE%l_onsite_param_index
       endif

       if(allocated( NN_TABLE%stoner_I_param_index    ))   then
          n1 = size( NN_TABLE%stoner_I_param_index    )  
           allocate( NN_TABLE_PY%stoner_I_param_index(n1) )  ;  NN_TABLE_PY%stoner_I_param_index      = NN_TABLE%stoner_I_param_index
       endif

       if(allocated( NN_TABLE%local_U_param_index     ))   then
          n1 = size( NN_TABLE%local_U_param_index     )  
           allocate( NN_TABLE_PY%local_U_param_index(n1)  )  ;  NN_TABLE_PY%local_U_param_index       = NN_TABLE%local_U_param_index
       endif

       if(allocated( NN_TABLE%plus_U_param_index      ))   then
          n1 = size( NN_TABLE%plus_U_param_index      )  
           allocate( NN_TABLE_PY%plus_U_param_index(n1)   )  ;  NN_TABLE_PY%plus_U_param_index        = NN_TABLE%plus_U_param_index
       endif

       if(allocated( NN_TABLE%soc_param_index         ))   then
          n1 = size( NN_TABLE%soc_param_index         )  
           allocate( NN_TABLE_PY%soc_param_index(n1)      )  ;  NN_TABLE_PY%soc_param_index           = NN_TABLE%soc_param_index
       endif

       if(allocated( NN_TABLE%ci_orb                  ))   then
          n1 = size( NN_TABLE%ci_orb                  )  
           allocate( NN_TABLE_PY%ci_orb(n1)               )  ;  NN_TABLE_PY%ci_orb                    = NN_TABLE%ci_orb
       endif

       if(allocated( NN_TABLE%cj_orb                  ))   then
          n1 = size( NN_TABLE%cj_orb                  )  
           allocate( NN_TABLE_PY%cj_orb(n1)               )  ;  NN_TABLE_PY%cj_orb                    = NN_TABLE%cj_orb
       endif

       if(allocated( NN_TABLE%p_class                 ))   then
          n1 = size( NN_TABLE%p_class                 )  
           allocate( NN_TABLE_PY%p_class(n1)              )  ;  NN_TABLE_PY%p_class                   = NN_TABLE%p_class
       endif

       if(allocated( NN_TABLE%site_cindex             ))   then
          n1 = size( NN_TABLE%site_cindex             )  
           allocate( NN_TABLE_PY%site_cindex(n1)          )  ;  NN_TABLE_PY%site_cindex               = NN_TABLE%site_cindex
       endif


       if(allocated( NN_TABLE_PY%i_coord                 ))   deallocate( NN_TABLE_PY%i_coord)         
       if(allocated( NN_TABLE_PY%j_coord                 ))   deallocate( NN_TABLE_PY%j_coord)         
       if(allocated( NN_TABLE_PY%Rij                     ))   deallocate( NN_TABLE_PY%Rij)             
       if(allocated( NN_TABLE_PY%R                       ))   deallocate( NN_TABLE_PY%R)               
       if(allocated( NN_TABLE_PY%R_nn                    ))   deallocate( NN_TABLE_PY%R_nn)            
       if(allocated( NN_TABLE_PY%R0_nn                   ))   deallocate( NN_TABLE_PY%R0_nn)           
       if(allocated( NN_TABLE_PY%local_moment            ))   deallocate( NN_TABLE_PY%local_moment)    
       if(allocated( NN_TABLE_PY%local_moment_rot        ))   deallocate( NN_TABLE_PY%local_moment_rot)
       if(allocated( NN_TABLE_PY%sk_index_set            ))   deallocate( NN_TABLE_PY%sk_index_set)    
       if(allocated( NN_TABLE_PY%cc_index_set            ))   deallocate( NN_TABLE_PY%cc_index_set)    
       if(allocated( NN_TABLE_PY%j_nn                    ))   deallocate( NN_TABLE_PY%j_nn)            

       if(allocated( NN_TABLE%i_coord                 ))  then 
          n1 = size( NN_TABLE%i_coord(:,1)            ) ;  n2 = size( NN_TABLE%i_coord(1,:)            ) ;
           allocate( NN_TABLE_PY%i_coord(n1,n2)          )  ;    NN_TABLE_PY%i_coord           = NN_TABLE%i_coord   
       endif

       if(allocated( NN_TABLE%j_coord                 ))  then 
          n1 = size( NN_TABLE%j_coord(:,1)            ) ;  n2 = size( NN_TABLE%j_coord(1,:)            ) ;
           allocate( NN_TABLE_PY%j_coord(n1,n2)          )  ;    NN_TABLE_PY%j_coord           = NN_TABLE%j_coord   
       endif

       if(allocated( NN_TABLE%Rij                     ))  then 
          n1 = size( NN_TABLE%Rij(:,1)                ) ;  n2 = size( NN_TABLE%Rij(1,:)                ) ;
           allocate( NN_TABLE_PY%Rij(n1,n2)              )  ;    NN_TABLE_PY%Rij               = NN_TABLE%Rij   
       endif

       if(allocated( NN_TABLE%R                       ))  then 
          n1 = size( NN_TABLE%R(:,1)                  ) ;  n2 = size( NN_TABLE%R(1,:)                  ) ;
           allocate( NN_TABLE_PY%R(n1,n2)                )  ;    NN_TABLE_PY%R                 = NN_TABLE%R   
       endif

       if(allocated( NN_TABLE%R_nn                    ))  then 
          n1 = size( NN_TABLE%R_nn(:,1)               ) ;  n2 = size( NN_TABLE%R_nn(1,:)               ) ;
           allocate( NN_TABLE_PY%R_nn(n1,n2)             )  ;    NN_TABLE_PY%R_nn              = NN_TABLE%R_nn   
       endif

       if(allocated( NN_TABLE%R0_nn                   ))  then 
          n1 = size( NN_TABLE%R0_nn(:,1)              ) ;  n2 = size( NN_TABLE%R0_nn(1,:)              ) ;
           allocate( NN_TABLE_PY%R0_nn(n1,n2)            )  ;    NN_TABLE_PY%R0_nn             = NN_TABLE%R0_nn
       endif

       if(allocated( NN_TABLE%local_moment            ))  then 
          n1 = size( NN_TABLE%local_moment(:,1)       ) ;  n2 = size( NN_TABLE%local_moment(1,:)       ) ;
           allocate( NN_TABLE_PY%local_moment(n1,n2)     )  ;    NN_TABLE_PY%local_moment      = NN_TABLE%local_moment
       endif

       if(allocated( NN_TABLE%local_moment_rot        ))  then 
          n1 = size( NN_TABLE%local_moment_rot(:,1)   ) ;  n2 = size( NN_TABLE%local_moment_rot(1,:)   ) ;
           allocate( NN_TABLE_PY%local_moment_rot(n1,n2) )  ;    NN_TABLE_PY%local_moment_rot  = NN_TABLE%local_moment_rot
       endif

       if(allocated( NN_TABLE%sk_index_set            ))  then 
          n1 = size( NN_TABLE%sk_index_set(:,1)       ) ;  n2 = size( NN_TABLE%sk_index_set(1,:)       ) ;
           allocate( NN_TABLE_PY%sk_index_set(0:n1-1,n2)     )  ;    NN_TABLE_PY%sk_index_set      = NN_TABLE%sk_index_set
       endif

       if(allocated( NN_TABLE%cc_index_set            ))  then 
          n1 = size( NN_TABLE%cc_index_set(:,1)       ) ;  n2 = size( NN_TABLE%cc_index_set(1,:)       ) ;
           allocate( NN_TABLE_PY%cc_index_set(0:n1-1,n2)     )  ;    NN_TABLE_PY%cc_index_set      = NN_TABLE%cc_index_set
       endif

       if(allocated( NN_TABLE%j_nn                    ))  then 
          n1 = size( NN_TABLE%j_nn(:,1)               ) ;  n2 = size( NN_TABLE%j_nn(1,:)               ) ;
           allocate( NN_TABLE_PY%j_nn(n1,n2)             )  ;    NN_TABLE_PY%j_nn              = NN_TABLE%j_nn
       endif

    elseif(imode .eq. 2) then
    
       NN_TABLE%flag_efield                =      NN_TABLE_PY%flag_efield
       NN_TABLE%flag_efield_frac           =      NN_TABLE_PY%flag_efield_frac
       NN_TABLE%flag_efield_cart           =      NN_TABLE_PY%flag_efield_cart
       NN_TABLE%onsite_tolerance           =      NN_TABLE_PY%onsite_tolerance
       NN_TABLE%efield                     =      NN_TABLE_PY%efield
       NN_TABLE%efield_origin              =      NN_TABLE_PY%efield_origin
       NN_TABLE%efield_origin_cart         =      NN_TABLE_PY%efield_origin_cart
       NN_TABLE%mysystem                   =      NN_TABLE_PY%mysystem
       NN_TABLE%n_neighbor                 =      NN_TABLE_PY%n_neighbor
       NN_TABLE%max_nn_pair                =      NN_TABLE_PY%max_nn_pair
       NN_TABLE%nnfilenm                   =      NN_TABLE_PY%nnfilenm
       NN_TABLE%nnfilenmo                  =      NN_TABLE_PY%nnfilenmo

       if(allocated( NN_TABLE%flag_site_cindex        ))   deallocate( NN_TABLE%flag_site_cindex)
       if(allocated( NN_TABLE%Dij                     ))   deallocate( NN_TABLE%Dij)
       if(allocated( NN_TABLE%Dij0                    ))   deallocate( NN_TABLE%Dij0)
       if(allocated( NN_TABLE%i_sign                  ))   deallocate( NN_TABLE%i_sign)
       if(allocated( NN_TABLE%j_sign                  ))   deallocate( NN_TABLE%j_sign)
       if(allocated( NN_TABLE%local_charge            ))   deallocate( NN_TABLE%local_charge)
       if(allocated( NN_TABLE%i_atom                  ))   deallocate( NN_TABLE%i_atom)
       if(allocated( NN_TABLE%j_atom                  ))   deallocate( NN_TABLE%j_atom)
       if(allocated( NN_TABLE%i_matrix                ))   deallocate( NN_TABLE%i_matrix)
       if(allocated( NN_TABLE%j_matrix                ))   deallocate( NN_TABLE%j_matrix)
       if(allocated( NN_TABLE%n_class                 ))   deallocate( NN_TABLE%n_class)
       if(allocated( NN_TABLE%n_nn                    ))   deallocate( NN_TABLE%n_nn)
       if(allocated( NN_TABLE%l_onsite_param_index    ))   deallocate( NN_TABLE%l_onsite_param_index)
       if(allocated( NN_TABLE%stoner_I_param_index    ))   deallocate( NN_TABLE%stoner_I_param_index)
       if(allocated( NN_TABLE%local_U_param_index     ))   deallocate( NN_TABLE%local_U_param_index)
       if(allocated( NN_TABLE%plus_U_param_index      ))   deallocate( NN_TABLE%plus_U_param_index)
       if(allocated( NN_TABLE%soc_param_index         ))   deallocate( NN_TABLE%soc_param_index)
       if(allocated( NN_TABLE%ci_orb                  ))   deallocate( NN_TABLE%ci_orb)
       if(allocated( NN_TABLE%cj_orb                  ))   deallocate( NN_TABLE%cj_orb)
       if(allocated( NN_TABLE%p_class                 ))   deallocate( NN_TABLE%p_class)
       if(allocated( NN_TABLE%site_cindex             ))   deallocate( NN_TABLE%site_cindex)

       if(allocated( NN_TABLE_PY%flag_site_cindex        ))   then
          n1 = size( NN_TABLE_PY%flag_site_cindex        )
           allocate( NN_TABLE%flag_site_cindex(n1)     )  ;  NN_TABLE%flag_site_cindex          = NN_TABLE_PY%flag_site_cindex
       endif

       if(allocated( NN_TABLE_PY%Dij                     ))   then
          n1 = size( NN_TABLE_PY%Dij                     )
           allocate( NN_TABLE%Dij(n1)                  )  ;  NN_TABLE%Dij                       = NN_TABLE_PY%Dij
       endif

       if(allocated( NN_TABLE_PY%Dij0                    ))   then
          n1 = size( NN_TABLE_PY%Dij0                    )
           allocate( NN_TABLE%Dij0(n1)                 )  ;  NN_TABLE%Dij0                      = NN_TABLE_PY%Dij0
       endif

       if(allocated( NN_TABLE_PY%i_sign                  ))   then
          n1 = size( NN_TABLE_PY%i_sign                  )
           allocate( NN_TABLE%i_sign(n1)               )  ;  NN_TABLE%i_sign                    = NN_TABLE_PY%i_sign
       endif

       if(allocated( NN_TABLE_PY%j_sign                  ))   then
          n1 = size( NN_TABLE_PY%j_sign                  )
           allocate( NN_TABLE%j_sign(n1)               )  ;  NN_TABLE%j_sign                    = NN_TABLE_PY%j_sign
       endif

       if(allocated( NN_TABLE_PY%local_charge            ))   then
          n1 = size( NN_TABLE_PY%local_charge            )
           allocate( NN_TABLE%local_charge(n1)         )  ;  NN_TABLE%local_charge              = NN_TABLE_PY%local_charge
       endif

       if(allocated( NN_TABLE_PY%i_atom                  ))   then
          n1 = size( NN_TABLE_PY%i_atom                  )
           allocate( NN_TABLE%i_atom(n1)               )  ;  NN_TABLE%i_atom                    = NN_TABLE_PY%i_atom
       endif

       if(allocated( NN_TABLE_PY%j_atom                  ))   then
          n1 = size( NN_TABLE_PY%j_atom                  )
           allocate( NN_TABLE%j_atom(n1)               )  ;  NN_TABLE%j_atom                    = NN_TABLE_PY%j_atom
       endif

       if(allocated( NN_TABLE_PY%i_matrix                ))   then
          n1 = size( NN_TABLE_PY%i_matrix                )
           allocate( NN_TABLE%i_matrix(n1)             )  ;  NN_TABLE%i_matrix                  = NN_TABLE_PY%i_matrix
       endif

       if(allocated( NN_TABLE_PY%j_matrix                ))   then
          n1 = size( NN_TABLE_PY%j_matrix                )
           allocate( NN_TABLE%j_matrix(n1)             )  ;  NN_TABLE%j_matrix                  = NN_TABLE_PY%j_matrix
       endif

       if(allocated( NN_TABLE_PY%n_class                 ))   then
          n1 = size( NN_TABLE_PY%n_class                 )
           allocate( NN_TABLE%n_class(n1)              )  ;  NN_TABLE%n_class                   = NN_TABLE_PY%n_class
       endif

       if(allocated( NN_TABLE_PY%n_nn                    ))   then
          n1 = size( NN_TABLE_PY%n_nn                    )
           allocate( NN_TABLE%n_nn(n1)                 )  ;  NN_TABLE%n_nn                      = NN_TABLE_PY%n_nn
       endif

       if(allocated( NN_TABLE_PY%l_onsite_param_index    ))   then
          n1 = size( NN_TABLE_PY%l_onsite_param_index    )
           allocate( NN_TABLE%l_onsite_param_index(n1) )  ;  NN_TABLE%l_onsite_param_index      = NN_TABLE_PY%l_onsite_param_index
       endif

       if(allocated( NN_TABLE_PY%stoner_I_param_index    ))   then
          n1 = size( NN_TABLE_PY%stoner_I_param_index    )
           allocate( NN_TABLE%stoner_I_param_index(n1) )  ;  NN_TABLE%stoner_I_param_index      = NN_TABLE_PY%stoner_I_param_index
       endif

       if(allocated( NN_TABLE_PY%local_U_param_index     ))   then
          n1 = size( NN_TABLE_PY%local_U_param_index     )
           allocate( NN_TABLE%local_U_param_index(n1)  )  ;  NN_TABLE%local_U_param_index       = NN_TABLE_PY%local_U_param_index
       endif

       if(allocated( NN_TABLE_PY%plus_U_param_index      ))   then
          n1 = size( NN_TABLE_PY%plus_U_param_index      )
           allocate( NN_TABLE%plus_U_param_index(n1)   )  ;  NN_TABLE%plus_U_param_index        = NN_TABLE_PY%plus_U_param_index
       endif

       if(allocated( NN_TABLE_PY%soc_param_index         ))   then
          n1 = size( NN_TABLE_PY%soc_param_index         )
           allocate( NN_TABLE%soc_param_index(n1)      )  ;  NN_TABLE%soc_param_index           = NN_TABLE_PY%soc_param_index
       endif

       if(allocated( NN_TABLE_PY%ci_orb                  ))   then
          n1 = size( NN_TABLE_PY%ci_orb                  )
           allocate( NN_TABLE%ci_orb(n1)               )  ;  NN_TABLE%ci_orb                    = NN_TABLE_PY%ci_orb
       endif

       if(allocated( NN_TABLE_PY%cj_orb                  ))   then
          n1 = size( NN_TABLE_PY%cj_orb                  )
           allocate( NN_TABLE%cj_orb(n1)               )  ;  NN_TABLE%cj_orb                    = NN_TABLE_PY%cj_orb
       endif

       if(allocated( NN_TABLE_PY%p_class                 ))   then
          n1 = size( NN_TABLE_PY%p_class                 )
           allocate( NN_TABLE%p_class(n1)              )  ;  NN_TABLE%p_class                   = NN_TABLE_PY%p_class
       endif

       if(allocated( NN_TABLE_PY%site_cindex             ))   then
          n1 = size( NN_TABLE_PY%site_cindex             )
           allocate( NN_TABLE%site_cindex(n1)          )  ;  NN_TABLE%site_cindex               = NN_TABLE_PY%site_cindex
       endif

       if(allocated( NN_TABLE%i_coord                 ))   deallocate( NN_TABLE%i_coord)
       if(allocated( NN_TABLE%j_coord                 ))   deallocate( NN_TABLE%j_coord)
       if(allocated( NN_TABLE%Rij                     ))   deallocate( NN_TABLE%Rij)
       if(allocated( NN_TABLE%R                       ))   deallocate( NN_TABLE%R)
       if(allocated( NN_TABLE%R_nn                    ))   deallocate( NN_TABLE%R_nn)
       if(allocated( NN_TABLE%R0_nn                   ))   deallocate( NN_TABLE%R0_nn)
       if(allocated( NN_TABLE%local_moment            ))   deallocate( NN_TABLE%local_moment)
       if(allocated( NN_TABLE%local_moment_rot        ))   deallocate( NN_TABLE%local_moment_rot)
       if(allocated( NN_TABLE%sk_index_set            ))   deallocate( NN_TABLE%sk_index_set)
       if(allocated( NN_TABLE%cc_index_set            ))   deallocate( NN_TABLE%cc_index_set)
       if(allocated( NN_TABLE%j_nn                    ))   deallocate( NN_TABLE%j_nn)

       if(allocated( NN_TABLE_PY%i_coord                 ))  then
          n1 = size( NN_TABLE_PY%i_coord(:,1)            ) ;  n2 = size( NN_TABLE_PY%i_coord(1,:)            ) ;
           allocate( NN_TABLE%i_coord(n1,n2)          )  ;    NN_TABLE%i_coord           = NN_TABLE_PY%i_coord
       endif

       if(allocated( NN_TABLE_PY%j_coord                 ))  then
          n1 = size( NN_TABLE_PY%j_coord(:,1)            ) ;  n2 = size( NN_TABLE_PY%j_coord(1,:)            ) ;
           allocate( NN_TABLE%j_coord(n1,n2)          )  ;    NN_TABLE%j_coord           = NN_TABLE_PY%j_coord
       endif

       if(allocated( NN_TABLE_PY%Rij                     ))  then
          n1 = size( NN_TABLE_PY%Rij(:,1)                ) ;  n2 = size( NN_TABLE_PY%Rij(1,:)                ) ;
           allocate( NN_TABLE%Rij(n1,n2)              )  ;    NN_TABLE%Rij               = NN_TABLE_PY%Rij
       endif

       if(allocated( NN_TABLE_PY%R                       ))  then
          n1 = size( NN_TABLE_PY%R(:,1)                  ) ;  n2 = size( NN_TABLE_PY%R(1,:)                  ) ;
           allocate( NN_TABLE%R(n1,n2)                )  ;    NN_TABLE%R                 = NN_TABLE_PY%R
       endif

       if(allocated( NN_TABLE_PY%R_nn                    ))  then
          n1 = size( NN_TABLE_PY%R_nn(:,1)               ) ;  n2 = size( NN_TABLE_PY%R_nn(1,:)               ) ;
           allocate( NN_TABLE%R_nn(n1,n2)             )  ;    NN_TABLE%R_nn              = NN_TABLE_PY%R_nn
       endif

       if(allocated( NN_TABLE_PY%R0_nn                   ))  then
          n1 = size( NN_TABLE_PY%R0_nn(:,1)              ) ;  n2 = size( NN_TABLE_PY%R0_nn(1,:)              ) ;
           allocate( NN_TABLE%R0_nn(n1,n2)            )  ;    NN_TABLE%R0_nn             = NN_TABLE_PY%R0_nn
       endif

       if(allocated( NN_TABLE_PY%local_moment            ))  then
          n1 = size( NN_TABLE_PY%local_moment(:,1)       ) ;  n2 = size( NN_TABLE_PY%local_moment(1,:)       ) ;
           allocate( NN_TABLE%local_moment(n1,n2)     )  ;    NN_TABLE%local_moment      = NN_TABLE_PY%local_moment
       endif

       if(allocated( NN_TABLE_PY%local_moment_rot        ))  then
          n1 = size( NN_TABLE_PY%local_moment_rot(:,1)   ) ;  n2 = size( NN_TABLE_PY%local_moment_rot(1,:)   ) ;
           allocate( NN_TABLE%local_moment_rot(n1,n2) )  ;    NN_TABLE%local_moment_rot  = NN_TABLE_PY%local_moment_rot
       endif

       if(allocated( NN_TABLE_PY%sk_index_set            ))  then
          n1 = size( NN_TABLE_PY%sk_index_set(:,1)       ) ;  n2 = size( NN_TABLE_PY%sk_index_set(1,:)       ) ;
           allocate( NN_TABLE%sk_index_set(0:n1-1,n2)     )  ;    NN_TABLE%sk_index_set      = NN_TABLE_PY%sk_index_set
       endif

       if(allocated( NN_TABLE_PY%cc_index_set            ))  then
          n1 = size( NN_TABLE_PY%cc_index_set(:,1)       ) ;  n2 = size( NN_TABLE_PY%cc_index_set(1,:)       ) ;
           allocate( NN_TABLE%cc_index_set(0:n1-1,n2)     )  ;    NN_TABLE%cc_index_set      = NN_TABLE_PY%cc_index_set
       endif

       if(allocated( NN_TABLE_PY%j_nn                    ))  then
          n1 = size( NN_TABLE_PY%j_nn(:,1)               ) ;  n2 = size( NN_TABLE_PY%j_nn(1,:)               ) ;
           allocate( NN_TABLE%j_nn(n1,n2)             )  ;    NN_TABLE%j_nn              = NN_TABLE_PY%j_nn
       endif

    endif

    return
endsubroutine

function init_energy_py() result(E_PY)
    type(energy_py)                        E_PY 

    E_PY%mysystem = 1

    if(allocated( E_PY%E    ))    deallocate( E_PY%E    )
    if(allocated( E_PY%dE   ))    deallocate( E_PY%dE   )
    if(allocated( E_PY%dORB ))    deallocate( E_PY%dORB )
    if(allocated( E_PY%ORB  ))    deallocate( E_PY%ORB  )
    if(allocated( E_PY%V    ))    deallocate( E_PY%V    )
    if(allocated( E_PY%SV   ))    deallocate( E_PY%SV   )
    if(allocated( E_PY%V2   ))    deallocate( E_PY%V2   )

    if(allocated( E_PY%E_BAND))   deallocate( E_PY%E_BAND )
    if(allocated( E_PY%E_TOT ))   deallocate( E_PY%E_TOT  )
    if(allocated( E_PY%F_OCC ))   deallocate( E_PY%F_OCC  )

    return
endfunction

subroutine copy_energy(E_PY, E, imode)
    type(energy    )                       E
    type(energy_py )                       E_PY
    integer(kind=sp)                       imode, n1, n2, n3

    if(imode .eq. 1) then
       E_PY%mysystem                       =   E%mysystem
       E_PY%E_F                            =   E%E_F     

       if(allocated( E_PY%E_BAND            ))   deallocate( E_PY%E_BAND   )
       if(allocated( E_PY%E_TOT             ))   deallocate( E_PY%E_TOT    )
       if(allocated( E_PY%F_OCC             ))   deallocate( E_PY%F_OCC    )
       if(allocated( E_PY%E                 ))   deallocate( E_PY%E    )
       if(allocated( E_PY%dE                ))   deallocate( E_PY%dE   )
       if(allocated( E_PY%dORB              ))   deallocate( E_PY%dORB )
       if(allocated( E_PY%ORB               ))   deallocate( E_PY%ORB  )
       if(allocated( E_PY%V                 ))   deallocate( E_PY%V    )
       if(allocated( E_PY%SV                ))   deallocate( E_PY%SV   )
       if(allocated( E_PY%V2                ))   deallocate( E_PY%V2   )

       if(allocated( E%E_BAND            )) then
          n1 = size( E%E_BAND(:)          ) 
           allocate( E_PY%E_BAND(n1)      ) ; E_PY%E_BAND=  E%E_BAND
       endif

       if(allocated( E%E_TOT             )) then
          n1 = size( E%E_TOT(:)           )
           allocate( E_PY%E_TOT(n1)       ) ; E_PY%E_TOT =  E%E_TOT 
       endif

       if(allocated( E%F_OCC             )) then
          n1 = size( E%F_OCC(:,1)         ) ; n2 = size( E%F_OCC(1,:)       )
           allocate( E_PY%F_OCC(n1,n2)    ) ; E_PY%F_OCC =  E%F_OCC
       endif

       if(allocated( E%E                 )) then
          n1 = size( E%E(:,1)             ) ; n2 = size( E%E(1,:)       )
           allocate( E_PY%E(n1,n2)        ) ; E_PY%E     =  E%E
       endif

       if(allocated( E%dE                )) then
          n1 = size( E%dE(:,1)            ) ; n2 = size( E%dE(1,:)      )
           allocate( E_PY%dE(n1,n2)       ) ; E_PY%dE     =  E%dE
       endif

       if(allocated( E%dORB              )) then
          n1 = size( E%dORB(:,1)          ) ; n2 = size( E%dORB(1,:)      )
           allocate( E_PY%dORB(n1,n2)     ) ; E_PY%dORB   =  E%dORB
       endif

       if(allocated( E%ORB               )) then
          n1 = size( E%ORB(:,1,1)         ) ; n2 = size( E%ORB(1,:,1)   );n3 = size( E%ORB(1,1,:)   ) 
           allocate( E_PY%ORB(n1,n2,n3)   ) ; E_PY%ORB   =  E%ORB
       endif

       if(allocated( E%V                 )) then
          n1 = size( E%V(:,1,1)           ) ; n2 = size( E%V(1,:,1)     );n3 = size( E%V(1,1,:)     ) 
           allocate( E_PY%V(n1,n2,n3)     ) ; E_PY%V     =  E%V
       endif

       if(allocated( E%SV                )) then
          n1 = size( E%SV(:,1,1)          ) ; n2 = size( E%SV(1,:,1)    );n3 = size( E%SV(1,1,:)    ) 
           allocate( E_PY%SV(n1,n2,n3)    ) ; E_PY%SV    =  E%SV
       endif

       if(allocated( E%V2                )) then
          n1 = size( E%V2(:,1,1)          ) ; n2 = size( E%V2(1,:,1)    );n3 = size( E%V2(1,1,:)    )
           allocate( E_PY%V2(n1,n2,n3)    ) ; E_PY%V2    =  E%V2
       endif


    elseif(imode .eq. 2) then

       E%mysystem                       =       E_PY%mysystem
       E%E_F                            =       E_PY%E_F     

       if(allocated( E%E_BAND            ))   deallocate( E%E_BAND)
       if(allocated( E%E_TOT             ))   deallocate( E%E_TOT )
       if(allocated( E%F_OCC             ))   deallocate( E%F_OCC )
       if(allocated( E%E                 ))   deallocate( E%E    )
       if(allocated( E%dE                ))   deallocate( E%dE   )
       if(allocated( E%dORB              ))   deallocate( E%dORB )
       if(allocated( E%ORB               ))   deallocate( E%ORB  )
       if(allocated( E%V                 ))   deallocate( E%V    )
       if(allocated( E%SV                ))   deallocate( E%SV   )
       if(allocated( E%V2                ))   deallocate( E%V2   )

       if(allocated( E_PY%E_BAND            )) then
          n1 = size( E_PY%E_BAND(:)          )
           allocate( E%E_BAND(n1)      ) ; E%E_BAND=  E_PY%E_BAND
       endif

       if(allocated( E_PY%E_TOT             )) then
          n1 = size( E_PY%E_TOT(:)           )
           allocate( E%E_TOT(n1)       ) ; E%E_TOT =  E_PY%E_TOT
       endif


       if(allocated( E_PY%F_OCC             )) then
          n1 = size( E_PY%F_OCC(:,1)         ) ; n2 = size( E_PY%F_OCC(1,:)       )
           allocate( E%F_OCC(n1,n2)    ) ; E%F_OCC =  E_PY%F_OCC
       endif

       if(allocated( E_PY%E                 )) then
          n1 = size( E_PY%E(:,1)             ) ; n2 = size( E_PY%E(1,:)       )
           allocate( E%E(n1,n2)              ) ; E%E     = E_PY%E  
       endif

       if(allocated( E_PY%dE                )) then
          n1 = size( E_PY%dE(:,1)            ) ; n2 = size( E_PY%dE(1,:)      )
           allocate( E%dE(n1,n2)       ) ; E%dE     =  E_PY%dE
       endif

       if(allocated( E_PY%dORB              )) then
          n1 = size( E_PY%dORB(:,1)            ) ; n2 = size( E_PY%dORB(1,:)      )
           allocate( E%dORB(n1,n2)       ) ; E%dORB     =  E_PY%dORB
       endif

       if(allocated( E_PY%ORB               )) then
          n1 = size( E_PY%ORB(:,1,1)         ) ; n2 = size( E_PY%ORB(1,:,1)   );n3 = size( E_PY%ORB(1,1,:)   ) 
           allocate( E%ORB(n1,n2,n3)   ) ; E%ORB   =  E_PY%ORB
       endif

       if(allocated( E_PY%V                 )) then
          n1 = size( E_PY%V(:,1,1)           ) ; n2 = size( E_PY%V(1,:,1)     );n3 = size( E_PY%V(1,1,:)     ) 
           allocate( E%V(n1,n2,n3)     ) ; E%V     =  E_PY%V  
       endif

       if(allocated( E_PY%SV                )) then
          n1 = size( E_PY%SV(:,1,1)          ) ; n2 = size( E%SV(1,:,1)    );n3 = size( E%SV(1,1,:)    ) 
           allocate( E%SV(n1,n2,n3)    ) ; E%SV    =  E_PY%SV 
       endif

       if(allocated( E_PY%V2                )) then
          n1 = size( E_PY%V2(:,1,1)          ) ; n2 = size( E%V2(1,:,1)    );n3 = size( E%V2(1,1,:)    )
           allocate( E%V2(n1,n2,n3)    ) ; E%V2    =  E_PY%V2
       endif

    endif

    return
endsubroutine

subroutine print_param_py(PINPT_PY, PPRAM_PY, pfileoutnm)
    use print_io
    implicit none
    type(incar_py), intent(inout) :: PINPT_PY
    type(params_py),intent(inout) :: PPRAM_PY
    type(params   )               :: PPRAM
    integer(kind=sp)                 i, j, k, pid_param_new
    integer(kind=sp)                 nsub, i_fix
    character(len=40)                param_name
    real(kind=dp)                    param(PPRAM_PY%param_nsub_max)
    character(len=20), external   :: int2str
    character(len=132),intent(inout):: pfileoutnm

    pid_param_new = 81
    open(pid_param_new, file = trim(pfileoutnm), status = 'unknown')

    write(pid_param_new, '(A,A)')'# ', trim(pfileoutnm)
    write(pid_param_new, '(A)')'# FITTED PARAMETERS FROM TBFITPY_MOD'

    if(PINPT_PY%flag_scissor) then
       write(pid_param_new,'(A,F8.2,A,I5,A)')'# SCISSOR: .TRUE. => EDFT(n,k) + ',PINPT_PY%r_scissor, &
                                                                     ' (if n >=',PINPT_PY%i_scissor,')' 
    endif

    !write(pid_param_new,'(A)')'# NOTE: Weight informations to obtain this parameter is written in WEIGHT.dat file'

    if(PPRAM_PY%flag_pfile_index) then
      write(pid_param_new, '(A20,A11)') ' PRINT_INDEX  ','.TRUE.'
    else
      write(pid_param_new, '(A20,A11)') ' PRINT_INDEX  ','.FALSE.'
    endif

    if(PPRAM_PY%flag_use_overlap) then
      write(pid_param_new, '(A20,A11)') ' USE_OVERLAP  ','.TRUE.'
    else
      write(pid_param_new, '(A20,A11)') ' USE_OVERLAP  ','.FALSE.'
    endif

    if(PPRAM_PY%flag_slater_koster) then
      write(pid_param_new, '(A20,A11)') '       IS_SK  ','.TRUE.'
    else
      write(pid_param_new, '(A20,A11)') ' USE_OVERLAP  ','.FALSE.'
    endif

    if(PPRAM_PY%slater_koster_type .gt. 10) then
      write(pid_param_new, '(A20,I11,A)') ' SK_SCALE_MODE  ', PPRAM_PY%slater_koster_type, '  # 11 = NRL type scaling'
      write(pid_param_new, '(A20,2x,F9.4,A)') ' L_BROADEN  ', PPRAM_PY%l_broaden, '  # angstrom '
    else
      write(pid_param_new, '(A20,I11,A)') ' SK_SCALE_MODE  ', PPRAM_PY%slater_koster_type
    endif

    write(pid_param_new,'(A)')' '   
 
    call copy_params(PPRAM_PY, PPRAM, 2)   
    do i = 1, PPRAM_PY%nparam
        call param_select(PPRAM, i, i_fix, param_name, param, nsub)
        call set_print_message(  i, i_fix, param_name, param, nsub, PPRAM_PY%flag_pfile_index, message)
        write(pid_param_new, '(A)')trim(message)
    enddo

    close(pid_param_new)

    return
endsubroutine

subroutine print_weight(PWGHT_PY, wfileoutnm)
    use print_io
    implicit none
    type(weight_py),intent(inout)    :: PWGHT_PY
    integer(kind=sp)                    neig, nk, pid_weight
    integer(kind=sp)                    ie, ik
    character(len=132),intent(inout) :: wfileoutnm

    pid_weight = 201
    open(pid_weight, file = trim(wfileoutnm), status = 'unknown')

    write(pid_weight, '(A,A)')'# ', trim(wfileoutnm)
    write(pid_weight, '(A)')'# FITTED WEIGHT FROM TBFITPY_MOD'
    write(pid_weight, '(A)')'#   do ik = 1, nk '
    write(pid_weight, '(A)')'#      print ( WEIGHT(1:neig, ik) )'
    write(pid_weight, '(A)')'#   enddo '

    neig = size(PWGHT_PY%WT(:,1))
    nk   = size(PWGHT_PY%WT(1,:))
    
    do ik = 1, nk   
        write(pid_weight, '(*(F14.5))') PWGHT_PY%WT(:,ik)
    enddo
    
    close(pid_weight)

    return
endsubroutine

subroutine print_target(PINPT_PY, PKPTS_PY, EDFT_PY, PWGHT_PY, PGEOM_PY, tfileoutnm)
    use print_io
    implicit none
    type(incar_py)  ,intent(inout)    :: PINPT_PY 
    type(kpoints_py),intent(inout)    :: PKPTS_PY 
    type(energy_py) ,intent(inout)    :: EDFT_PY  
    type(weight_py) ,intent(inout)    :: PWGHT_PY
    type(poscar_py) ,intent(inout)    :: PGEOM_PY
    type(incar     )                     PINPT
    type(energy    )                     EDFT
    type(weight    )                     PWGHT
    character(len=132),intent(inout)  :: tfileoutnm

    PINPT_PY%fnamelog = 'TEST'
    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_weight(PWGHT_PY, PWGHT, 2)
    call copy_energy(EDFT_PY, EDFT, 2)
    
    call print_energy_weight(PKPTS_PY%kpoint, PKPTS_PY%nkpoint, EDFT, PWGHT, PGEOM_PY%neig, &
                             PINPT, trim(tfileoutnm), PINPT%flag_get_band_order)

    return
endsubroutine

subroutine print_etba(PINPT_PY, PKPTS_PY, ETBA_PY, EDFT_PY, PWGHT_PY, PGEOM_PY, suffix, flag_use_overlap)
    use print_io
    implicit none
    type(incar_py)  ,intent(inout)    :: PINPT_PY
    type(kpoints_py),intent(inout)    :: PKPTS_PY
    type(energy_py) ,intent(inout)    :: ETBA_PY
    type(energy_py) ,intent(inout)    :: EDFT_PY
    type(weight_py) ,intent(inout)    :: PWGHT_PY
    type(poscar_py) ,intent(inout)    :: PGEOM_PY
    type(incar     )                     PINPT
    type(kpoints   )                     PKPTS
    type(weight    )                     PWGHT
    logical                              flag_use_overlap
    character(len=132), intent(inout) :: suffix

    call copy_incar(PINPT_PY, PINPT, 2)
    call copy_weight(PWGHT_PY, PWGHT, 2)
    call copy_kpoints(PKPTS_PY, PKPTS, 2)

    if(PINPT_PY%flag_print_energy_diff) then
    
      call print_energy(PKPTS, ETBA_PY%E, EDFT_PY%E, ETBA_PY%V, ETBA_PY%SV, &
                        PGEOM_PY%neig, PGEOM_PY%init_erange, PGEOM_PY%nband, &
                        PINPT, PWGHT, flag_use_overlap, .TRUE.,trim(suffix))
    else
 
      call print_energy(PKPTS, ETBA_PY%E, ETBA_PY%E, ETBA_PY%V, ETBA_PY%SV, &
                        PGEOM_PY%neig, PGEOM_PY%init_erange, PGEOM_PY%nband, &
                        PINPT, PWGHT, flag_use_overlap, .TRUE.,trim(suffix))
    endif

    return
endsubroutine

subroutine set_verbose(iverbose_mod)
    implicit none
    integer(kind=sp)    iverbose_mod

    if(iverbose_mod .eq. 1) then
        iverbose = 1 ; print_mode = 3  ! for full report
    elseif(iverbose_mod .eq. 2) then
        iverbose = 2 ; print_mode = 99 ! no report
    else
        iverbose = 2 ; print_mode = 99 ! no report
    endif
   
    return
endsubroutine

subroutine update_free_param(param, nparam_free, PPRAM)
    implicit none
    type(params)                                :: PPRAM
    integer(kind=sp)                               nparam_free
    real(kind=dp)                                  param(nparam_free)
    integer(kind=sp)                               iparam

    do iparam = 1, nparam_free
        PPRAM%param(PPRAM%iparam_free(iparam))  = param(iparam)
    enddo

    return
endsubroutine

subroutine copy_all2(PKPTS_PY,PGEOM_PY,NN_TABLE_PY,ETBA_PY, PKPTS,PGEOM,NN_TABLE,ETBA,imode)
    implicit none
    type(kpoints_py)                    :: PKPTS_PY
    type(poscar_py)                     :: PGEOM_PY
    type(hopping_py)                    :: NN_TABLE_PY
    type(energy_py)                     :: ETBA_PY
    type(kpoints )                      :: PKPTS
    type(poscar  )                      :: PGEOM
    type(hopping )                      :: NN_TABLE
    type(energy  )                      :: ETBA
    integer(kind=sp)                       imode

    call copy_kpoints(PKPTS_PY,    PKPTS, imode)
    call copy_poscar(PGEOM_PY,     PGEOM, imode)
    call copy_hopping(NN_TABLE_PY, NN_TABLE, imode)
    call copy_energy(ETBA_PY,      ETBA, imode)

    return
endsubroutine

subroutine copy_all(PKPTS_PY, PGEOM_PY, PWGHT_PY, NN_TABLE_PY, EDFT_PY, ETBA_PY, &
                     PKPTS, PGEOM, PWGHT, NN_TABLE , EDFT, ETBA, imode)
    implicit none
    type(kpoints_py)                    :: PKPTS_PY
    type(weight_py)                     :: PWGHT_PY
    type(poscar_py)                     :: PGEOM_PY
    type(hopping_py)                    :: NN_TABLE_PY
    type(energy_py)                     :: EDFT_PY, ETBA_PY
    type(kpoints )                      :: PKPTS
    type(weight  )                      :: PWGHT
    type(poscar  )                      :: PGEOM
    type(hopping )                      :: NN_TABLE
    type(energy  )                      :: EDFT, ETBA
    integer(kind=sp)                       imode

    call copy_kpoints(PKPTS_PY,    PKPTS, imode)
    call copy_poscar(PGEOM_PY,     PGEOM, imode)
    call copy_weight(PWGHT_PY,     PWGHT, imode)
    call copy_hopping(NN_TABLE_PY, NN_TABLE, imode)
    call copy_energy(EDFT_PY,      EDFT, imode)
    call copy_energy(ETBA_PY,      ETBA, imode)

    return
endsubroutine

!function random() result(r)
!    implicit none
!    real(kind=dp) :: r
!    call random_number(r)
!endfunction

!subroutine random_init(iseed)
!    implicit none
!    integer(kind=sp), intent(in) :: iseed
!    integer(kind=sp),allocatable :: seed(:)
!    integer(kind=sp)             :: n
!    call random_seed(size=n)
!    allocate(seed(n))
!    seed = iseed
!    call random_seed(put=seed)
!    deallocate(seed)
!endsubroutine

endmodule
