module set_default
    use parameters, only : incar, params, dos, berry, kpoints, weight, poscar, energy, &
                           hopping, gainp, replot, unfold, onsite_tolerance
    implicit none

    contains

    subroutine init_incar(PINPT)
        type(incar   )      :: PINPT

        PINPT%flag_get_band=.true.
        PINPT%flag_get_band_order=.false.
        PINPT%flag_get_band_order_print_only =.false.
        PINPT%flag_erange=.false.
        PINPT%flag_plot_fit=.false.
        PINPT%flag_print_energy_diff = .false.
        if(.not. PINPT%flag_filenm_gnuplot_parse) then
          PINPT%filenm_gnuplot = 'gnuBAND-TB.gpi' ! default
        else
          PINPT%filenm_gnuplot = trim(PINPT%filenm_gnuplot_parse)
        endif
        PINPT%flag_plot_stm_image = .false.
        PINPT%flag_plot_eigen_state=.false.
        PINPT%flag_plot_wavefunction = .true.
        if(.not. PINPT%flag_lorbit_parse) PINPT%flag_print_orbital=.false.
        PINPT%flag_print_single=.false.
        PINPT%flag_print_energy_singlek=.false.
        PINPT%flag_phase = .true.
        PINPT%flag_fit_degeneracy = .false.
        if(.not. PINPT%flag_lorbit_parse) PINPT%flag_get_orbital=.false.
        if(.not. PINPT%flag_lorbit_parse) PINPT%flag_print_mag=.false.

        PINPT%flag_print_proj=.false.
        PINPT%flag_print_proj_sum=.false.
        PINPT%flag_print_proj_atom=.false.
        PINPT%nproj_sum = 0
        if(allocated(PINPT%proj_atom))  deallocate(PINPT%proj_atom)
        if(allocated(PINPT%proj_natom)) deallocate(PINPT%proj_natom)

        PINPT%ncirc_dichroism = 0
        PINPT%flag_use_weight = .false.  ! whether read weight factor for fit from PFILE or not
        if(.not. PINPT%flag_wfile_parse) PINPT%flag_set_weight_from_file = .false.
        PINPT%flag_get_dos=.false.
        PINPT%flag_get_z2=.false.
        PINPT%flag_get_zak_phase=.false.
        PINPT%flag_zak_separate=.false.
        PINPT%flag_get_parity=.false.
        PINPT%flag_get_symmetry=.false.
        PINPT%flag_berryc_separate=.false.
        PINPT%flag_zak_kfile_read = .false.
        PINPT%flag_get_berry_curvature = .false.
        PINPT%flag_get_circ_dichroism = .false.
        PINPT%flag_berry = .false.
        PINPT%flag_collinear=.false.
        PINPT%flag_noncollinear=.false.
        PINPT%flag_local_charge=.false.
        PINPT%flag_plus_U=.false.
        PINPT%flag_scissor = .false.
        PINPT%flag_load_nntable= .false.
        PINPT%flag_sparse = .false.
        PINPT%flag_get_effective_ham=.false.
        PINPT%flag_write_unformatted=.false.
        PINPT%flag_fit_orbital = .false.
        PINPT%orbital_fit_smearing = 0.1d0 ! default
        PINPT%electronic_temperature = 11.60452251451435653502d0 ! electronic temp T; default, so that k_B * T = 0.001 eV
        PINPT%flag_get_total_energy = .false.
#ifdef SPGLIB
        PINPT%flag_spglib = .false.
#endif
        if(.not. PINPT%flag_miter_parse) PINPT%miter = 30      ! default 
        if(.not. PINPT%flag_mxfit_parse) PINPT%mxfit = 1       ! default 
        PINPT%ftol  = 0.00001d0 ! default 
        PINPT%ptol  = 0.00001d0 ! default 
        PINPT%fdiff = 0.001d0   ! default
        PINPT%ispin = 1 ! default (1 : nonmag, 2: collinear & noncollinear)
        PINPT%ispinor = 1 ! default (1 : nonmag, collinear 2: noncollinear)
        PINPT%nspin = 1 ! default (1 : nonmag, 2: collinear, 1: noncollinear)
        PINPT%rcut_orb_plot = 5 ! default (unit = angstrom)
        PINPT%nn_max(1:3) = 3 ! default
        PINPT%lmmax=9 ! default 9 => s(1),p(2-4),d(5-9), 3 => s(1), p(2), d(3)

        PINPT%flag_ga_with_lmdif=.false.
        PINPT%flag_pso_with_lmdif=.false.
        PINPT%flag_pso_report_particles=.true.
        PINPT%flag_report_geom = .true.

        PINPT%flag_distribute_nkp = .false.

        PINPT%iseed = 123 ! random seed
        if(.not. PINPT%flag_pso_verbose_parse) PINPT%pso_verbose = 1 

        PINPT%flag_get_unfold = .false.

        return
    endsubroutine

    subroutine init_params(PPRAM, PINPT)
        type(params  )      :: PPRAM
        type(incar   )      :: PINPT

        PPRAM%param_nsub_max = 1 ! default if SK_SCALE_TYPE <= 10
        PPRAM%flag_slater_koster = .true.
        PPRAM%slater_koster_type = 1 ! default scaling method
        PPRAM%flag_use_overlap=.false.
        PPRAM%flag_pfile_index=.false.
        PPRAM%flag_set_param_const=.false.
        PPRAM%flag_fit_plain = .FALSE.
        if( .not. PINPT%flag_pfile_parse) then
          PPRAM%pfilenm='PARAM_FIT.dat' !default
        elseif(PINPT%flag_pfile_parse) then
          PPRAM%pfilenm=PINPT%pfilenm_parse
        endif
        PPRAM%pfileoutnm='PARAM_FIT.new.dat' !default
        PPRAM%l_broaden = 0.15 ! angstrong unit, default for cutoff-function broadening
        PPRAM%nparam= 0 ! number of TB parameters
        PPRAM%nparam_const = 0 ! number of constratint acting on parameters
        PPRAM%niter = 0 ! number of interations performed in fitting procedures
        PPRAM%pso_nparticles                          = 50
        PPRAM%pso_c1                                  = 0.3d0
        PPRAM%pso_c2                                  = 0.4d0
        PPRAM%pso_w                                   = 0.2d0
        PPRAM%pso_max_noise_amplitude                 = 5.0d0
        PPRAM%pso_report_ratio                        = 0.2d0 ! 20 percent of best particles would be reported
        PPRAM%pso_iseed                               = 123 ! random seed
        PPRAM%pso_mode                                = 'pso'
        PPRAM%pso_miter                               = 10

        if(allocated(PPRAM%param))           deallocate(PPRAM%param)
        if(allocated(PPRAM%param_best))      deallocate(PPRAM%param_best)
        if(allocated(PPRAM%param_nrl))       deallocate(PPRAM%param_nrl)
        if(allocated(PPRAM%iparam_free))     deallocate(PPRAM%iparam_free)
        if(allocated(PPRAM%iparam_free_nrl)) deallocate(PPRAM%iparam_free_nrl)
        if(allocated(PPRAM%param_name))      deallocate(PPRAM%param_name)
        if(allocated(PPRAM%c_const))         deallocate(PPRAM%c_const)
        if(allocated(PPRAM%param_const))     deallocate(PPRAM%param_const)
!       if(allocated(PPRAM%param_const_best)) deallocate(PPRAM%param_const_best)
        if(allocated(PPRAM%param_const_nrl)) deallocate(PPRAM%param_const_nrl)
        if(allocated(PPRAM%param_nsub))      deallocate(PPRAM%param_nsub)

        if(allocated(PPRAM%pso_cost_history)) deallocate(PPRAM%pso_cost_history)
        if(allocated(PPRAM%pso_cost_history_i)) deallocate(PPRAM%pso_cost_history_i)
        if(allocated(PPRAM%pso_pbest_history)) deallocate(PPRAM%pso_pbest_history)
        if(allocated(PPRAM%cost_history))     deallocate(PPRAM%cost_historY)

        return
    endsubroutine

    subroutine init_dos(PINPT_DOS, PINPT)
        type(dos     )      :: PINPT_DOS
        type(incar   )      :: PINPT

        if(allocated(PINPT_DOS%dos_kpoint))     deallocate(PINPT_DOS%dos_kpoint)
        if(allocated(PINPT_DOS%dos_erange))     deallocate(PINPT_DOS%dos_erange)
        if(allocated(PINPT_DOS%dos_tot))        deallocate(PINPT_DOS%dos_tot)
        if(allocated(PINPT_DOS%dos_ensurf))     deallocate(PINPT_DOS%dos_ensurf)
        if(allocated(PINPT_DOS%dos_ldos_atom))  deallocate(PINPT_DOS%dos_ldos_atom)
        if(allocated(PINPT_DOS%ldos_tot))       deallocate(PINPT_DOS%ldos_tot)

        return
    endsubroutine

    subroutine init_berry(PINPT_BERRY, PINPT)
        type(berry   )      :: PINPT_BERRY
        type(incar   )      :: PINPT

        PINPT_BERRY%flag_wcc_phase = .false.
        PINPT_BERRY%flag_bc_phase  = .true.
        PINPT_BERRY%flag_z2_phase  = .false.
        PINPT_BERRY%flag_zak_phase = .false.
        PINPT_BERRY%noccupied      = 0
        PINPT_BERRY%flag_print_hamiltonian_parity = .false.
        PINPT_BERRY%flag_print_hamiltonian_symmetry = .false.
        PINPT_BERRY%flag_parity_phase = .false.
        PINPT_BERRY%flag_symmetry_phase = .false.

        if(allocated(PINPT_BERRY%wcc_erange))           deallocate(PINPT_BERRY%wcc_erange)
        if(allocated(PINPT_BERRY%wcc       ))           deallocate(PINPT_BERRY%wcc       )
        if(allocated(PINPT_BERRY%wcc_kpath ))           deallocate(PINPT_BERRY%wcc_kpath )
        if(allocated(PINPT_BERRY%wcc_kpoint))           deallocate(PINPT_BERRY%wcc_kpoint)
        if(allocated(PINPT_BERRY%wcc_kpoint_reci))      deallocate(PINPT_BERRY%wcc_kpoint_reci)
        if(allocated(PINPT_BERRY%wcc_polarization))     deallocate(PINPT_BERRY%wcc_polarization)
        if(allocated(PINPT_BERRY%wcc_chern ))           deallocate(PINPT_BERRY%wcc_chern )     
        if(allocated(PINPT_BERRY%zak_erange))           deallocate(PINPT_BERRY%zak_erange)     
        if(allocated(PINPT_BERRY%zak_phase ))           deallocate(PINPT_BERRY%zak_phase )     
        if(allocated(PINPT_BERRY%zak_kpath ))           deallocate(PINPT_BERRY%zak_kpath )     
        if(allocated(PINPT_BERRY%zak_kpoint))           deallocate(PINPT_BERRY%zak_kpoint)     
        if(allocated(PINPT_BERRY%zak_kpoint_reci))      deallocate(PINPT_BERRY%zak_kpoint_reci)
        if(allocated(PINPT_BERRY%polarization))         deallocate(PINPT_BERRY%polarization)
        if(allocated(PINPT_BERRY%z2_erange ))           deallocate(PINPT_BERRY%z2_erange )
        if(allocated(PINPT_BERRY%z2_axis   ))           deallocate(PINPT_BERRY%z2_axis   )
        if(allocated(PINPT_BERRY%z2_kpoint ))           deallocate(PINPT_BERRY%z2_kpoint )
        if(allocated(PINPT_BERRY%z2_kpoint_reci))       deallocate(PINPT_BERRY%z2_kpoint_reci )
        if(allocated(PINPT_BERRY%z2_wcc    ))           deallocate(PINPT_BERRY%z2_wcc    )
        if(allocated(PINPT_BERRY%z2_chern  ))           deallocate(PINPT_BERRY%z2_chern  )
        if(allocated(PINPT_BERRY%z2_polarization))      deallocate(PINPT_BERRY%z2_polarization)
        if(allocated(PINPT_BERRY%bc_axis   ))           deallocate(PINPT_BERRY%bc_axis   )
        if(allocated(PINPT_BERRY%bc_erange ))           deallocate(PINPT_BERRY%bc_erange )
        if(allocated(PINPT_BERRY%berrycurv ))           deallocate(PINPT_BERRY%berrycurv )
        if(allocated(PINPT_BERRY%omega     ))           deallocate(PINPT_BERRY%omega     )
        if(allocated(PINPT_BERRY%parity_kpoint))        deallocate(PINPT_BERRY%parity_kpoint)
        if(allocated(PINPT_BERRY%parity_kpoint_reci))   deallocate(PINPT_BERRY%parity_kpoint_reci)
        if(allocated(PINPT_BERRY%parity_kpoint_name))   deallocate(PINPT_BERRY%parity_kpoint_name)
        if(allocated(PINPT_BERRY%symmetry_kpoint))      deallocate(PINPT_BERRY%symmetry_kpoint)
        if(allocated(PINPT_BERRY%symmetry_kpoint_reci)) deallocate(PINPT_BERRY%symmetry_kpoint_reci)
        if(allocated(PINPT_BERRY%symmetry_kpoint_name)) deallocate(PINPT_BERRY%symmetry_kpoint_name)

        return
    endsubroutine

    subroutine init_kpoints(PKPTS, PINPT)
        type(kpoints )      :: PKPTS  
        type(incar   )      :: PINPT  

        if(.not. PINPT%flag_kfile_parse) then
          PKPTS%kfilenm='KPOINTS_BAND' !default
        else
          PKPTS%kfilenm=PINPT%kfilenm_parse
        endif
        PKPTS%kline_type = 'vasp' ! default. 'vasp' or 'fleur'
        PKPTS%kunit = 'A' !default 'A' : angstrom or 'R' : reciprocal unit is available
        PKPTS%kreduce = 1 ! default
        PKPTS%n_ndiv = 1 ! default
        PKPTS%idiv_mode= 1 ! default

        if(allocated(PKPTS%ndiv           ))  deallocate(PKPTS%ndiv           )
        if(allocated(PKPTS%kpoint         ))  deallocate(PKPTS%kpoint         )
        if(allocated(PKPTS%kline          ))  deallocate(PKPTS%kline          )
        if(allocated(PKPTS%kpoint_reci    ))  deallocate(PKPTS%kpoint_reci    )
        if(allocated(PKPTS%k_name         ))  deallocate(PKPTS%k_name         )
        if(allocated(PKPTS%k_name2        ))  deallocate(PKPTS%k_name2        )
        if(allocated(PKPTS%k_name_index   ))  deallocate(PKPTS%k_name_index   )

        return
    endsubroutine

    subroutine init_weight(PWGHT)
        type(weight  )      :: PWGHT  

        PWGHT%itarget_e_start = 1 ! default
        PWGHT%efile_ef   = 0d0
        PWGHT%efile_type = 'user'
        PWGHT%efilenmu=' '
        PWGHT%efilenmd=' '
        PWGHT%read_energy_column_index = 2 ! default, read 2nd column as target energy (NM or spin_up)
        PWGHT%read_energy_column_index_dn = 2 ! default, read 2nd column as target energy (spin_dn)
        PWGHT%flag_weight_default = .true.
        PWGHT%flag_weight_orb = .false.
        PWGHT%iband = 1
        PWGHT%fband = -9999
        PWGHT%vbmd  = 0
        PWGHT%cbmd  = 0
        PWGHT%vbmt  = 0
        PWGHT%cbmt  = 0
        PWGHT%nweight = 0
        PWGHT%ie_cutoff = 0
        PWGHT%max_len_strip_kp = 0
        PWGHT%max_len_strip_tb = 0
        PWGHT%max_len_strip_df = 0
        PWGHT%max_len_strip_ob = 0
        PWGHT%max_len_strip_st = 0

        if(allocated(PWGHT%WT            ))   deallocate(PWGHT%WT            )
        if(allocated(PWGHT%DEGENERACY_WT ))   deallocate(PWGHT%DEGENERACY_WT )
        if(allocated(PWGHT%PENALTY_ORB   ))   deallocate(PWGHT%PENALTY_ORB   )
        if(allocated(PWGHT%strip_kp      ))   deallocate(PWGHT%strip_kp      )
        if(allocated(PWGHT%strip_tb      ))   deallocate(PWGHT%strip_tb      )
        if(allocated(PWGHT%strip_df      ))   deallocate(PWGHT%strip_df      )
        if(allocated(PWGHT%strip_wt      ))   deallocate(PWGHT%strip_wt      )
        if(allocated(PWGHT%strip_kp_orb  ))   deallocate(PWGHT%strip_kp_orb  )
        if(allocated(PWGHT%strip_tb_orb  ))   deallocate(PWGHT%strip_tb_orb  )
        if(allocated(PWGHT%strip_tb      ))   deallocate(PWGHT%strip_tb      )
        if(allocated(PWGHT%strip_site    ))   deallocate(PWGHT%strip_site    )
        if(allocated(PWGHT%strip_pen_orb ))   deallocate(PWGHT%strip_pen_orb )
        if(allocated(PWGHT%strip_kp_deg  ))   deallocate(PWGHT%strip_kp_deg  )
        if(allocated(PWGHT%strip_tb_deg  ))   deallocate(PWGHT%strip_tb_deg  )
        if(allocated(PWGHT%strip_df_deg  ))   deallocate(PWGHT%strip_df_deg  )
        if(allocated(PWGHT%strip_wt_deg  ))   deallocate(PWGHT%strip_wt_deg  )

        return
    endsubroutine

    subroutine init_poscar(PGEOM)
        type(poscar  )      :: PGEOM

        PGEOM%title='' ! default
        PGEOM%gfilenm='POSCAR-TB' !default
        PGEOM%flag_set_ribbon=.false.
        PGEOM%flag_print_only_ribbon_geom=.false.
        PGEOM%r_origin(1:3) = 0d0
        PGEOM%ngrid = -1
        PGEOM%stm_ngrid = -1
        
        if(allocated(PGEOM%nelect               )) deallocate(PGEOM%nelect) 
        if(allocated(PGEOM%i_spec               )) deallocate(PGEOM%i_spec               )
        if(allocated(PGEOM%c_spec               )) deallocate(PGEOM%c_spec               )
        if(allocated(PGEOM%flag_site_cindex     )) deallocate(PGEOM%flag_site_cindex     )
        if(allocated(PGEOM%site_cindex          )) deallocate(PGEOM%site_cindex          )
        if(allocated(PGEOM%spec                 )) deallocate(PGEOM%spec                 )
        if(allocated(PGEOM%spec_equiv           )) deallocate(PGEOM%spec_equiv           )
!       if(allocated(PGEOM%spec_define          )) deallocate(PGEOM%spec_define          )
        if(allocated(PGEOM%a_coord              )) deallocate(PGEOM%a_coord              )
        if(allocated(PGEOM%a_coord_cart         )) deallocate(PGEOM%a_coord_cart         )
        if(allocated(PGEOM%o_coord              )) deallocate(PGEOM%o_coord              )
        if(allocated(PGEOM%o_coord_cart         )) deallocate(PGEOM%o_coord_cart         )
        if(allocated(PGEOM%i_eff_orb            )) deallocate(PGEOM%i_eff_orb            )
        if(allocated(PGEOM%n_orbital            )) deallocate(PGEOM%n_orbital            )
        if(allocated(PGEOM%c_orbital            )) deallocate(PGEOM%c_orbital            )
        if(allocated(PGEOM%orb_sign             )) deallocate(PGEOM%orb_sign             )
        if(allocated(PGEOM%ispec                )) deallocate(PGEOM%ispec                )
        if(allocated(PGEOM%l_quantum            )) deallocate(PGEOM%l_quantum            )
        if(allocated(PGEOM%orb_n_quantum        )) deallocate(PGEOM%orb_n_quantum        )
        if(allocated(PGEOM%n_quantum            )) deallocate(PGEOM%n_quantum            )
        if(allocated(PGEOM%z_eff_nuc            )) deallocate(PGEOM%z_eff_nuc            )
        if(allocated(PGEOM%nn_pair              )) deallocate(PGEOM%nn_pair              )
        if(allocated(PGEOM%nn_dist              )) deallocate(PGEOM%nn_dist              )
        if(allocated(PGEOM%nn_r0                )) deallocate(PGEOM%nn_r0                )
        if(allocated(PGEOM%orb_index            )) deallocate(PGEOM%orb_index            )
        if(allocated(PGEOM%local_charge         )) deallocate(PGEOM%local_charge         )
        if(allocated(PGEOM%local_moment         )) deallocate(PGEOM%local_moment         )
        if(allocated(PGEOM%local_moment_rot     )) deallocate(PGEOM%local_moment_rot     )
        if(allocated(PGEOM%spg_rotations        )) deallocate(PGEOM%spg_rotations        )
        if(allocated(PGEOM%spg_translations     )) deallocate(PGEOM%spg_translations     )
        if(allocated(PGEOM%spg_wyckoffs         )) deallocate(PGEOM%spg_wyckoffs         )
        if(allocated(PGEOM%spg_equivalent_atoms )) deallocate(PGEOM%spg_equivalent_atoms )
        if(allocated(PGEOM%spg_a_coord_operated )) deallocate(PGEOM%spg_a_coord_operated )
        if(allocated(PGEOM%spg_a_coord_primitive)) deallocate(PGEOM%spg_a_coord_primitive)
        if(allocated(PGEOM%spg_spec_primitive   )) deallocate(PGEOM%spg_spec_primitive   )
        if(allocated(PGEOM%spg_det_w            )) deallocate(PGEOM%spg_det_w            )
        if(allocated(PGEOM%spg_tr_w             )) deallocate(PGEOM%spg_tr_w             )
        if(allocated(PGEOM%spg_type_w           )) deallocate(PGEOM%spg_type_w           )

        allocate(PGEOM%nelect(2)) ; PGEOM%nelect = -1.d0

        return
    endsubroutine

!   type(energy  )      :: EN

    subroutine init_hopping(NN_TABLE)
       type(hopping )      :: NN_TABLE

       NN_TABLE%nnfilenmo = 'hopping.dat' ! default
       NN_TABLE%nnfilenm  = 'hopping.dat' 
       NN_TABLE%onsite_tolerance = onsite_tolerance ! default defined in parameters.f90
       NN_TABLE%efield_origin(1:3) = 0d0 ! default
       NN_TABLE%flag_efield = .false.
       NN_TABLE%flag_efield_frac = .false.
       NN_TABLE%flag_efield_cart = .false.
    
       if(allocated(NN_TABLE%i_atom                ))    deallocate(NN_TABLE%i_atom                )
       if(allocated(NN_TABLE%j_atom                ))    deallocate(NN_TABLE%j_atom                )
       if(allocated(NN_TABLE%i_coord               ))    deallocate(NN_TABLE%i_coord               )
       if(allocated(NN_TABLE%j_coord               ))    deallocate(NN_TABLE%j_coord               )
       if(allocated(NN_TABLE%Rij                   ))    deallocate(NN_TABLE%Rij                   )
       if(allocated(NN_TABLE%R                     ))    deallocate(NN_TABLE%R                     )
       if(allocated(NN_TABLE%Dij                   ))    deallocate(NN_TABLE%Dij                   )
       if(allocated(NN_TABLE%Dij0                  ))    deallocate(NN_TABLE%Dij0                  )
       if(allocated(NN_TABLE%i_matrix              ))    deallocate(NN_TABLE%i_matrix              )
       if(allocated(NN_TABLE%ci_orb                ))    deallocate(NN_TABLE%ci_orb                )
       if(allocated(NN_TABLE%i_sign                ))    deallocate(NN_TABLE%i_sign                )
       if(allocated(NN_TABLE%j_matrix              ))    deallocate(NN_TABLE%j_matrix              )
       if(allocated(NN_TABLE%cj_orb                ))    deallocate(NN_TABLE%cj_orb                )
       if(allocated(NN_TABLE%j_sign                ))    deallocate(NN_TABLE%j_sign                )
       if(allocated(NN_TABLE%p_class               ))    deallocate(NN_TABLE%p_class               )
       if(allocated(NN_TABLE%n_class               ))    deallocate(NN_TABLE%n_class               )
       if(allocated(NN_TABLE%sk_index_set          ))    deallocate(NN_TABLE%sk_index_set          )
       if(allocated(NN_TABLE%cc_index_set          ))    deallocate(NN_TABLE%cc_index_set          )
       if(allocated(NN_TABLE%tij                   ))    deallocate(NN_TABLE%tij                   )
       if(allocated(NN_TABLE%sij                   ))    deallocate(NN_TABLE%sij                   )
       if(allocated(NN_TABLE%n_nn                  ))    deallocate(NN_TABLE%n_nn                  )
       if(allocated(NN_TABLE%R_nn                  ))    deallocate(NN_TABLE%R_nn                  )
       if(allocated(NN_TABLE%R0_nn                 ))    deallocate(NN_TABLE%R0_nn                 )
       if(allocated(NN_TABLE%j_nn                  ))    deallocate(NN_TABLE%j_nn                  )
       if(allocated(NN_TABLE%l_onsite_param_index  ))    deallocate(NN_TABLE%l_onsite_param_index  )
       if(allocated(NN_TABLE%site_cindex           ))    deallocate(NN_TABLE%site_cindex           )
       if(allocated(NN_TABLE%flag_site_cindex      ))    deallocate(NN_TABLE%flag_site_cindex      )
       if(allocated(NN_TABLE%local_charge          ))    deallocate(NN_TABLE%local_charge          )
       if(allocated(NN_TABLE%local_moment          ))    deallocate(NN_TABLE%local_moment          )
       if(allocated(NN_TABLE%local_moment_rot      ))    deallocate(NN_TABLE%local_moment_rot      )
       if(allocated(NN_TABLE%stoner_I_param_index  ))    deallocate(NN_TABLE%stoner_I_param_index  )
       if(allocated(NN_TABLE%local_U_param_index   ))    deallocate(NN_TABLE%local_U_param_index   )
       if(allocated(NN_TABLE%plus_U_param_index    ))    deallocate(NN_TABLE%plus_U_param_index    )
       if(allocated(NN_TABLE%soc_param_index       ))    deallocate(NN_TABLE%soc_param_index       )
       if(allocated(NN_TABLE%tij_file              ))    deallocate(NN_TABLE%tij_file              )
       if(allocated(NN_TABLE%sij_file              ))    deallocate(NN_TABLE%sij_file              )
       if(allocated(NN_TABLE%i_eff_orb             ))    deallocate(NN_TABLE%i_eff_orb             )
        
       return
    endsubroutine

    subroutine init_gainp(PKAIA)
       type(gainp   )      :: PKAIA   

       PKAIA%mgen    = 500
       PKAIA%npop    = 100
       PKAIA%ngene   = 6
       PKAIA%pcross  = 0.85d0
       PKAIA%pmutmn  = 0.0005d0
       PKAIA%pmutmx  = 0.25d0
       PKAIA%pmut    = 0.005d0
       PKAIA%imut    = 2
       PKAIA%fdif    = 1.0d0
       PKAIA%irep    = 3
       PKAIA%ielite  = 0
       PKAIA%ivrb    = 0
       PKAIA%convtol = 0.0001d0
       PKAIA%convwin = 20.0d0
       PKAIA%iguessf = 0.1d0
       PKAIA%iseed   = 999

       return
    endsubroutine

    subroutine init_replot(PRPLT)
       type(replot  )      :: PRPLT   

       PRPLT%flag_replot       = .false.
       PRPLT%flag_replot_dos   = .false.
       PRPLT%flag_replot_ldos  = .false.
       PRPLT%flag_replot_sldos = .false.
       PRPLT%flag_replot_didv  = .false.
       PRPLT%flag_replot_proj_band  = .false.
       PRPLT%replot_nproj_sum  = 0
       PRPLT%replot_nldos_sum  = 0
       PRPLT%replot_nband      = 0
       PRPLT%replot_sldos_fname= 'SLDOS.replot'
       PRPLT%replot_didv_fname = 'DIDV.replot'
       PRPLT%replot_dos_fname  = 'DOS.replot'
       PRPLT%replot_dos_smearing  =  0.025d0
       PRPLT%replot_dos_emin      = -10d0
       PRPLT%replot_dos_emax      =  10d0
       PRPLT%replot_dos_nediv     =  1000
       PRPLT%replot_sldos_cell    = (/1,1,1/)
       PRPLT%r_origin             = 0d0
       PRPLT%bond_cut             = 3d0
       PRPLT%flag_replot_formatted= .true.  ! default: read formatted band_structure_TBA file

       if(allocated(PRPLT%replot_dos_erange             ))  deallocate(PRPLT%replot_dos_erange             )
       if(allocated(PRPLT%replot_dos_tot                ))  deallocate(PRPLT%replot_dos_tot                )
       if(allocated(PRPLT%replot_dos_ntot               ))  deallocate(PRPLT%replot_dos_ntot               )
       if(allocated(PRPLT%replot_ldos_natom             ))  deallocate(PRPLT%replot_ldos_natom             )
       if(allocated(PRPLT%replot_ldos_atom              ))  deallocate(PRPLT%replot_ldos_atom              )
       if(allocated(PRPLT%replot_ldos_tot               ))  deallocate(PRPLT%replot_ldos_tot               )
       if(allocated(PRPLT%replot_sldos_sum              ))  deallocate(PRPLT%replot_sldos_sum              )
       if(allocated(PRPLT%replot_axis_print_mag         ))  deallocate(PRPLT%replot_axis_print_mag         )
       if(allocated(PRPLT%flag_replot_print_single      ))  deallocate(PRPLT%flag_replot_print_single      )
       if(allocated(PRPLT%flag_replot_write_unformatted ))  deallocate(PRPLT%flag_replot_write_unformatted )
       if(allocated(PRPLT%replot_proj_natom             ))  deallocate(PRPLT%replot_proj_natom             )
       if(allocated(PRPLT%replot_proj_atom              ))  deallocate(PRPLT%replot_proj_atom              )

       return
    endsubroutine

    subroutine init_unfold(PUFLD, PINPT)
       type(unfold) :: PUFLD
       type(incar)  :: PINPT

       PUFLD%unfold_flag_sparse = .false.
       PUFLD%unfold_smearing    = 0.025d0       ! gaussian smearing in spectral function
       PUFLD%unfold_kshift      = 0d0
       PUFLD%unfold_emin        = -10d0
       PUFLD%unfold_emax        =  10d0          
       PUFLD%unfold_nediv       =  1000         ! number of division between emin:emax
       PUFLD%unfold_nemax       =  0            ! number of maximum eigenvalues between emin:emax (only sparse true)
       PUFLD%unfold_kfilenm_PBZ = 'KPOINTS-PBZ' ! defined in primitive cell BZ
       PUFLD%unfold_kfilenm_SBZ = 'KPOINTS-SBZ' ! defined in supercell BZ
       PUFLD%unfold_gfilenm_PC  = 'POSCAR-PC'   ! geometry with primitive cell
       PUFLD%unfold_gfilenm_SC  = 'POSCAR-SC'   ! geometry with supercell
       PUFLD%unfold_ifilenm_PBZ = 'INCAR_PC'    ! input file for primitive cell
       PUFLD%unfold_kgrid       = 0
       if(allocated(PUFLD%unfold_kpoint_PBZ))       deallocate(PUFLD%unfold_kpoint_PBZ)

       return
    endsubroutine
endmodule
