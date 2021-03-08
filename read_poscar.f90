#include "alias.inc"
subroutine read_geometry(PGEOM, PINPT, NN_TABLE, PPRAM)
  use parameters, only : incar, poscar, hopping, params, k_B
  use mpi_setup
  use print_io
  implicit none
  type (incar )            :: PINPT
  type (poscar)            :: PGEOM
  type (hopping)           :: NN_TABLE
  type (params)            :: PPRAM   
  integer*4                   mpierr
  logical                     flag_gfile_exist
  integer*4                   nelect(2)

  inquire(file=PGEOM%gfilenm,exist=flag_gfile_exist)

  if(.not. flag_gfile_exist) then
    write(message,'(A,A,A)')'  !WARN! ',trim(PGEOM%gfilenm),' does not exist!! Exit...' ; write_msgi
    kill_job
  endif

  ! set ribbon if required
  if(PGEOM%flag_set_ribbon) call set_ribbon_geom(PGEOM) 

  ! read geometry
  write(message,*)' '  ; write_msgi
  write(message,*)'#- READING INPUT GEOMETRY FILE: ',trim(PGEOM%gfilenm)  ; write_msgi
  call read_poscar(PINPT, PGEOM, NN_TABLE)

  !set equivalent atom if defined in CONSTRAINT set
  call set_equiv_spec(PPRAM, PGEOM)

  !set nband and initial band if flag_erange = .false.
  if(.not.PINPT%flag_erange) then
    PGEOM%init_erange = 1 ! default
    PGEOM%fina_erange = PGEOM%neig*PINPT%ispinor ! default
    PGEOM%nband = PGEOM%neig*PINPT%ispinor ! default
    if(PINPT%flag_sparse) then
      ! initial number of ncsr will not exeed n_neighbor.
      ! Note that this will be adjusted later in the "get_eig" routine by
      ! constructing the Hamiltonian with Compressed Sparse array format.
      if(PINPT%feast_nemax .le. 0) then
        PGEOM%nband = PGEOM%neig*PINPT%ispinor ! default (do not save memory...)
        PGEOM%init_erange = 1
        PGEOM%fina_erange = PGEOM%neig*PINPT%ispinor
      elseif(PINPT%feast_nemax .ge. 1) then
        if(PINPT%feast_nemax .gt. PGEOM%neig * PINPT%ispinor) then
          write(message,'(A,I0,A)')'    !WARN! The NE_MAX (',PINPT%feast_nemax,') of EWINDOW tag is larger than the eigenvalues (NEIG)'   ; write_msgi
          write(message,'(A,I0,A)')'           of the system (',PGEOM%neig * PINPT%ispinor,'). Hence, we enforce NEMAX = NEIG.'           ; write_msgi
          write(message,'(A,I0,A)')'           Otherwise, you can reduce the expected NE_MAX within the EWINDOW with a proper guess.'     ; write_msgi
          PGEOM%nband = PGEOM%neig * PINPT%ispinor
          PINPT%feast_nemax = PGEOM%nband
        else
          PGEOM%nband = PINPT%feast_nemax
        endif
        PGEOM%init_erange = 1
        PGEOM%fina_erange = PGEOM%nband
      endif
    endif
  endif

  write(message,'(A,i8)')'    NBAND:',PGEOM%nband  ; write_msgi

  ! set orbital index (s,p,d...) according to their name specified in PFILE (only work with flag_slater_koster = .true.)
  if(PPRAM%flag_slater_koster) then
    call set_orbital_index(PGEOM, PINPT%lmmax)
  endif

  ! set effective nunclear charge which will be used in eigen state plot or stm plot
  if(PPRAM%flag_slater_koster .and. (PINPT%flag_plot_stm_image .or. PINPT%flag_plot_eigen_state) ) then
    call set_effective_nuclear_charge(PPRAM, PGEOM)
  endif

  ! setup for effective hamiltonian calculation
  if(PINPT%flag_sparse .and. PINPT%flag_get_effective_ham) then
    write(message,'(A)')'    !WARN! The EWINDOW tag and LOWDIN cannot be used simulatneously in the current version.'  ; write_msgi
    write(message,'(A)')'           Exit program...' ; write_msgi
    kill_job
  elseif(.not. PINPT%flag_sparse .and. PINPT%flag_get_effective_ham) then
    call set_effective_orbital_index(PINPT, PGEOM, NN_TABLE)
  endif

  ! set nelect for total energy calculation (if flag_get_total_energy = .true.)
   if(PINPT%flag_get_total_energy) then
     if(allocated(PGEOM%nelect)) then
       nelect = PGEOM%nelect
       deallocate(PGEOM%nelect)
     endif
     if(PINPT%nspin .eq. 2) then
       allocate(PGEOM%nelect(2))
       PGEOM%nelect = nelect
       write(message,'(A)')'  L_TOTEN: .TRUE.' ; write_msgi
       write(message,'(A,F12.5,A)')'   ELTEMP: ', PINPT%electronic_temperature, ' (electronic temperature (in K))' ; write_msgi
       write(message,'(A,F12.5,A)')'           ', PINPT%electronic_temperature*k_B, ' (gaussian broadening = kB*T)' ; write_msgi

       if(nelect(1) .lt. 0d0) then
         write(message,'(A)')'    !WARN! The total nergy calculation is requested but number of electons are not explicitly specified. ' ; write_msgi
         write(message,'(A)')'           Please check NELECT tag for the information. Exit...' ; write_msgi
         kill_job
       else
         write(message,'(A,F12.5)')'   NELECT: (up)', PGEOM%nelect(1) ; write_msgi
         write(message,'(A,F12.5)')'   NELECT: (dn)', PGEOM%nelect(2) ; write_msgi
       endif
     else
       allocate(PGEOM%nelect(1))
       PGEOM%nelect = nelect(1)
       write(message,'(A)')'  L_TOTEN: .TRUE.' ; write_msgi
       write(message,'(A,F12.5,A)')'   ELTEMP: ', PINPT%electronic_temperature, ' (electronic temperature (in K))' ; write_msgi
       write(message,'(A,F12.5,A)')'           ', PINPT%electronic_temperature*k_B, ' (gaussian broadening = kB*T)' ; write_msgi
       if(nelect(1) .lt. 0d0) then
         write(message,'(A)')'    !WARN! The total nergy calculation is requested but number of electons are not explicitly specified. ' ; write_msgi
         write(message,'(A)')'           Please check NELECT tag for the information. Exit...' ; write_msgi
         kill_job
       else
         write(message,'(A,F12.5)')'   NELECT: ', PGEOM%nelect(1) ; write_msgi
       endif
     endif
   else
     write(message,'(A)')'  L_TOTEN: .FALSE.' ; write_msgi
   endif

  write(message,*)'#- END READING GEOMETRY FILE ---------------------'  ; write_msgi
  write(message,*)' '  ; write_msgi

  return
endsubroutine
subroutine read_poscar(PINPT,PGEOM,NN_TABLE)
  use parameters,  only : incar, poscar, hopping, pid_geom
  use do_math, only : rotate_vector, inv
  use print_io
  use mpi_setup
  implicit none
  integer*4, parameter     :: max_orb_temp = 20
  integer*4                   mpierr
  integer*4                   i_continue,nitems
  integer*4                   i,j,ii,linecount, i_dummy, i_dummy1, i_dummy2, i_dummy3
  integer*4                   iorb
  integer*4                   pos_index(20), i_index, min_pos, n_flag_sel
  real*8, allocatable      :: local_charge_(:), local_moment_(:,:)
  real*8                      t_latt_inv(3,3)
  real*8                      t_coord(3)
  character*264               inputline, str2lowcase
  character*132               fname, inputline_dummy, dummy
  character*40                desc_str,dummy1,dummy2,dummy3
  real*8, allocatable      :: temp_orbital_sign(:,:)
  character*8, allocatable :: temp_orbital(:,:)
  character*8                 temp
  character*20,allocatable :: site_c_index_(:)
  logical,     allocatable :: flag_site_c_index_(:)
  character*10                site_index
  character*20                locpot_index
  character(*), parameter  :: func = 'read_poscar'
  logical                     flag_skip, flag_read_moment, flag_moment_cart
  external                    nitems, str2lowcase
  type (incar )            :: PINPT  
  type (poscar)            :: PGEOM  
  type (hopping)           :: NN_TABLE

  fname        = PGEOM%gfilenm
  PGEOM%n_spec = 0
  PGEOM%flag_selective = .false.
  flag_read_moment = .false.
  flag_moment_cart = .false.
  pos_index = 0

  open (pid_geom, FILE=fname,iostat=i_continue)
  linecount = 0
  ii = 0
line: do
        read(pid_geom,'(A)',iostat=i_continue) inputline
        if(i_continue<0) exit               ! end of file reached
        if(i_continue>0) then 
          write(message,*)'Unknown error reading file:',trim(fname),func  ; write_msgi
        endif

        if(linecount .eq. 0) then 
          call check_comment(inputline,linecount,i,flag_skip)
          if(flag_skip) linecount = linecount + 1
        else
          call check_comment(inputline,linecount,i,flag_skip)
          if(flag_skip) linecount = linecount + 1
          if (flag_skip) cycle
        endif
        linecount = linecount + 1

        ! head
         if(linecount .eq. 1) then
           PGEOM%system_name = trim(inputline)
           write(message,'(A,A)')'   SYSTEM:  ', trim(PGEOM%system_name)  ; write_msgi
           cycle

        ! scaling factor
         elseif(linecount .eq. 2) then
           read(inputline,*,iostat=i_continue) PGEOM%a_scale
           write(message,'(A,F15.8)')'  A_SCALE:  ',PGEOM%a_scale  ; write_msgi
           cycle

        ! lattice parameter
         elseif(linecount .eq. 3 ) then
           backspace(pid_geom)
           write(message,'(A            )')'          _______________________________________________'; write_msgi
           do i=1,3
             read(pid_geom,'(A)',iostat=i_continue) inputline
             read(inputline,*,iostat=i_continue) PGEOM%a_latt(1:3,i)
             write(message,'(A,i1,A,3F15.8)')'  LATT A',i,':  ',PGEOM%a_latt(1:3,i)  ; write_msgi
           enddo
           call get_reci(PGEOM%b_latt(:,1), PGEOM%b_latt(:,2), PGEOM%b_latt(:,3), &
                         PGEOM%a_latt(:,1), PGEOM%a_latt(:,2), PGEOM%a_latt(:,3))
           write(message,'(A            )')'          -----------------------------------------------'; write_msgi
           do i=1,3
             write(message,'(A,i1,A,3F15.8)')'  RECI B',i,':  ',PGEOM%b_latt(1:3,i)  ; write_msgi
           enddo
           write(message,'(A            )')'          _______________________________________________'; write_msgi
           linecount = linecount + 2
           cycle

        ! species name and number of atoms
         elseif(linecount .eq. 6 ) then
           PGEOM%n_spec=nitems(inputline)
           allocate( PGEOM%c_spec(PGEOM%n_spec), PGEOM%i_spec(PGEOM%n_spec) )
           write(message,'(A,i8)')'   N_SPEC:',PGEOM%n_spec  ; write_msgi
           read(inputline,*,iostat=i_continue) PGEOM%c_spec(1:PGEOM%n_spec)
           read(pid_geom,'(A)',iostat=i_continue) inputline
           linecount = linecount + 1
           call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle
           read(inputline,*,iostat=i_continue) PGEOM%i_spec(1:PGEOM%n_spec)
           PGEOM%n_atom=sum ( PGEOM%i_spec(1:PGEOM%n_spec) )

           allocate( PGEOM%spec(PGEOM%n_atom) )
           allocate( PGEOM%spec_equiv(PGEOM%n_atom) )
!          allocate( PGEOM%spec_define(PGEOM%n_atom) )
           do i=1,PGEOM%n_spec
             PGEOM%spec( sum(PGEOM%i_spec(1:i)) -PGEOM%i_spec(i)+1 : sum(PGEOM%i_spec(1:i)) ) = i
           enddo
           PGEOM%spec_equiv = PGEOM%spec ! initialize. 

           write(message,'(A,i8)')'   N_ATOM:',PGEOM%n_atom  ; write_msgi
           allocate( PGEOM%a_coord(3,PGEOM%n_atom), &
                     PGEOM%a_coord_cart(3,PGEOM%n_atom), &
                     PGEOM%n_orbital(PGEOM%n_atom), &
                     local_charge_(PGEOM%n_atom*max_orb_temp), &
                     local_moment_(3,PGEOM%n_atom*max_orb_temp), &
                     site_c_index_(PGEOM%n_atom), flag_site_c_index_(PGEOM%n_atom), &
                     temp_orbital(max_orb_temp, PGEOM%n_atom), &
                     temp_orbital_sign(max_orb_temp, PGEOM%n_atom ) )
                     local_charge_ = 0d0 ! initialize as zero
                     local_moment_ = 0d0 ! initialize as zero
                     flag_site_c_index_ = .false. ! initialize as .false.
                     temp_orbital_sign = 1d0 ! initialize as unity.
           do i=1,PGEOM%n_spec
             if(i .eq. 1)then
               write(message,'(A,I2,A,A4,1x,i8)')'   SPEC',i,':',trim(PGEOM%c_spec(i)),PGEOM%i_spec(i)  ; write_msgi
             else
               write(message,'(A,I2,A,A4,1x,i8)')'   SPEC',i,':',trim(PGEOM%c_spec(i)),PGEOM%i_spec(i)  ; write_msgi
             endif
           enddo

        ! constraint and coordinate type
         elseif(linecount .eq. 8 ) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'S' .or. desc_str(1:1) .eq. 's') then 
             PGEOM%flag_selective = .true.
            !write(message,'(A)')' L_CONSTR:  .TRUE.'  ; write_msgi
           elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then 
             PGEOM%flag_selective = .false.
             PGEOM%flag_direct=.true.
             PGEOM%flag_cartesian=.false.
             linecount = linecount + 1
            !write(message,'(A)')' L_CONSTR:  .TRUE.'  ; write_msgi
             write(message,'(A)')' C_CRDTYP:  DIRECT'  ; write_msgi
           elseif(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
                  desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then 
             PGEOM%flag_selective = .false.
             PGEOM%flag_direct=.false.
             PGEOM%flag_cartesian=.true.
             linecount = linecount + 1
            !write(message,'(A)')' L_CONSTR:  .FALSE.'  ; write_msgi
             write(message,'(A)')' C_CRDTYP:  CARTESIAN'  ; write_msgi
           endif
         elseif(linecount .eq. 9 ) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
              desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then 
             PGEOM%flag_direct=.false.
             PGEOM%flag_cartesian=.true.
             write(message,'(A)')' C_CRDTYP:  CARTESIAN'  ; write_msgi
           elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then 
             PGEOM%flag_direct=.true.
             PGEOM%flag_cartesian=.false.
             write(message,'(A)')' C_CRDTYP:  DIRECT'  ; write_msgi
           endif

         ! atomic coordinate & atomic orbital information
         elseif(linecount .eq. 10 ) then
           backspace(pid_geom)

           if(PGEOM%flag_selective) then
             n_flag_sel = 3
           elseif(.not. PGEOM%flag_selective) then
             n_flag_sel = 0
           endif

           do i=1,PGEOM%n_atom
             read(pid_geom,'(A)',iostat=i_continue) inputline
             linecount = linecount + 1
             call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle
             i_dummy1=index(inputline,'#')
             if(i_dummy1 .ne. 0) inputline = inputline(1:i_dummy1-1) !check comment
             inputline = str2lowcase(inputline)  

             !check 'charge', 'moment', 'site_c_index'
             pos_index(1)=index(trim(inputline),'charge') !check whether 'charge' has been set up
             locpot_index='charge'
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local_pot')
               locpot_index='local_pot' 
             endif
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local.pot')
               locpot_index='local.pot'
             endif
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local_potential')
               locpot_index='local_potential'
             endif
             if(pos_index(1) .eq. 0) then
               pos_index(1)=index(trim(inputline),'local.potential')
               locpot_index='local.potential'
             endif

             pos_index(2)=index(trim(inputline),'moment') !check whether 'moment' has been set up
             site_index='site_index'
             pos_index(3)=index(trim(inputline),'site_index') !check whether 'site_index' has been set up
             if(pos_index(3) .eq. 0) then 
               site_index='site_indx'
               pos_index(3)=index(trim(inputline),'site_indx')
             endif
             if(pos_index(3) .eq. 0) then 
               site_index='site_idx'
               pos_index(3)=index(trim(inputline),'site_idx')
             endif
             if(pos_index(3) .eq. 0) then 
               site_index='site_name'
               pos_index(3)=index(trim(inputline),'site_name')
             endif
             min_pos = 999999
             do i_index = 1, 3
               if(pos_index(i_index) .ne. 0) then
                 if(pos_index(i_index) .lt. min_pos) then
                   min_pos = pos_index(i_index)
                 endif
               endif
             enddo
            
             if(min_pos .eq. 999999) then  ! if no other values are defined.
               inputline_dummy = inputline
             else
               inputline_dummy = inputline(1:min_pos-1)
             endif
             ii = ii + 1
             PGEOM%n_orbital(i) = nitems(inputline_dummy) - 3 - n_flag_sel

             if(pos_index(1) .ne. 0) then ! charge
               if(pos_index(2) .ne. 0) then 
                 call strip_off(trim(inputline), dummy, trim(locpot_index), 'moment', 1)
               elseif(pos_index(2) .eq. 0 .and. pos_index(3) .ne. 0) then
                 call strip_off(trim(inputline), dummy, trim(locpot_index), trim(site_index), 1)
               else
                 call strip_off(trim(inputline), dummy, trim(locpot_index), '', 2)
               endif
               i_dummy = nitems(dummy)
               if(i_dummy .ne. PGEOM%n_orbital(i)) then
                 write(message,'(A,I6)')'  !WARNING! Charge setting is inproper. Number of items should be same as N_ORBITAL(i). iatom=',i  ; write_msgi
                 write(message,'(A)')   '  !WARNING! Please check GFILE. Exit program...'  ; write_msgi
                 stop
               else
                 read(dummy,*)local_charge_(ii:ii+PGEOM%n_orbital(i)-1)
               endif
             endif ! charge

             if(pos_index(2) .ne. 0) then ! moment
               if(pos_index(3) .ne. 0) then
                 if(index(trim(inputline),'moment.r') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.r', trim(site_index), 1)
                 elseif(index(trim(inputline),'moment.c') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.c', trim(site_index), 1)
                   flag_moment_cart = .true.
                 else
                   call strip_off(trim(inputline), dummy, 'moment', trim(site_index), 1)
                 endif
               else
                 if(index(trim(inputline),'moment.r') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.r', '', 2)
                 elseif(index(trim(inputline),'moment.c') .gt. 0) then
                   call strip_off(trim(inputline), dummy, 'moment.c', '', 2)
                   flag_moment_cart = .true.
                 else
                   call strip_off(trim(inputline), dummy, 'moment', '', 2)
                 endif
               endif

               i_dummy = nitems(dummy)
               if(PINPT%flag_collinear) then
                 if(i_dummy .ne. PGEOM%n_orbital(i)) then
                   write(message,'(A)')'  !WARNING! Moment setting is inproper. Number of items should be same as N_ORBITAL(i) in the collinear setting.'  ; write_msgi
                   write(message,'(A,I6)')'  !WARNING! Please check GFILE. Exit program...  iatom=',i  ; write_msgi
                   stop
                 else 
                   read(dummy,*)local_moment_(1,ii:ii+PGEOM%n_orbital(i)-1)
                 endif
               endif
               if(PINPT%flag_noncollinear) then
                 if(i_dummy .ne. PGEOM%n_orbital(i)*3) then
                   write(message,'(A)')'  !WARNING! Moment setting is inproper. Number of items should be same as N_ORBITAL(i)*3 in the non-collinear setting.'  ; write_msgi
                   write(message,'(A,I6)')'  !WARNING! N_ORBITAL(i) * 3 = number of items (M, theta, phi) or (Mx,My,Mz). Please check GFILE. Exit program...  iatom=',i  ; write_msgi
                   stop
                 else
                   read(dummy,*)((local_moment_(j,i_dummy),j=1,3),i_dummy = ii, ii+PGEOM%n_orbital(i)-1)
                 endif
               endif ! moment
             endif

             if(pos_index(3) .ne. 0) then ! site_index 
               !if site_index is predefined in the POSCAR with Site_index tag: flag_site_cindex = .true. and use it as site_cindex
               call strip_off(trim(inputline), dummy, trim(site_index), '', 2)
               i_dummy = nitems(dummy)
               if(i_dummy .ne. 1) then
                 write(message,'(A)')'  !WARNING! Site_index setting is inproper. Number of items should be one and the data type should be character(20)'  ; write_msgi
                 write(message,'(A,I6)')'  !WARNING! Please check GFILE. Exit program...  iatom=',i  ; write_msgi
                 stop
               else
                 read(dummy,*)site_c_index_(i)
                 flag_site_c_index_(i) = .true.
               endif
             else
               !if site_index is not predefined in the POSCAR with Site_index tag: flag_site_cindex = .false. and use atom_name+atom_number as site_cindex
               dummy2 = PGEOM%c_spec(PGEOM%spec(i))
               if( i .lt. 10) then
                 write(dummy3,'(I1)') i
               elseif( i .ge. 10 .and. i .lt. 100) then
                 write(dummy3,'(I2)') i
               elseif( i .ge. 100 .and. i .lt. 1000) then
                 write(dummy3,'(I3)') i
               elseif( i .ge. 1000 .and. i .lt. 10000) then
                 write(dummy3,'(I4)') i
               elseif( i .ge. 10000 .and. i .lt. 100000) then
                 write(dummy3,'(I5)') i
               elseif( i .ge. 100000 .and. i .lt. 1000000) then
                 write(dummy3,'(I6)') i
               endif
               i_dummy2 = len_trim(dummy2)
               i_dummy3 = len_trim(dummy3)
               write(site_c_index_(i),'(2A)') dummy2(1:i_dummy2),dummy3(1:i_dummy3)
               site_c_index_(i) = str2lowcase(site_c_index_(i))
             endif !site_index

             ii = ii+PGEOM%n_orbital(i)-1

             ! read atomic coordinate & atomic orbital & magnetic moment info (if 'moment' tag is provided)              
             if(PGEOM%n_orbital(i) .eq. 0) then
               read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i)
               temp_orbital(1,i)='_na_'
             elseif(PGEOM%n_orbital(i) .ge. 1) then
               if(PGEOM%flag_selective) then
                 read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i), desc_str, desc_str, desc_str, &
                                                     temp_orbital(1:PGEOM%n_orbital(i),i)
               else
                 read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i), &
                                                     temp_orbital(1:PGEOM%n_orbital(i),i)
                 do i_dummy = 1, PGEOM%n_orbital(i)
                   temp = trim(temp_orbital(i_dummy,i))
                   if(temp(1:1) .eq. '-') then
                     temp_orbital(i_dummy,i) = temp(2:)
                     temp_orbital_sign(i_dummy,i) = -1d0
                   endif
                 enddo
               endif
             endif

           enddo !n_atom
           
           ! number of atomic orbitals (Hamiltonian matrix size)
           PGEOM%neig=sum(PGEOM%n_orbital(1:PGEOM%n_atom))
           PGEOM%neig_total = PGEOM%neig * PINPT%ispin
           PGEOM%nbasis = PGEOM%neig 

           if(PGEOM%neig .eq. 0) then
             write(message,'(A)')'  !! Check geometry input file. atomic orbital is not asigned!'  ; write_msgi
           elseif(PGEOM%neig .gt. 1) then
             write(message,'(A,i8)')'  N_ORBIT:',PGEOM%neig  ; write_msgi
             PGEOM%max_orb=maxval( PGEOM%n_orbital(:) )
             allocate( PGEOM%c_orbital(PGEOM%max_orb,PGEOM%n_atom) ) ; PGEOM%c_orbital = '_na_'
             allocate( PGEOM%orb_sign(PGEOM%max_orb,PGEOM%n_atom) )  ; PGEOM%orb_sign = 1d0
             PGEOM%c_orbital(1:PGEOM%max_orb,1:PGEOM%n_atom) = temp_orbital(1:PGEOM%max_orb,1:PGEOM%n_atom)
             PGEOM%orb_sign(1:PGEOM%max_orb,1:PGEOM%n_atom)  = temp_orbital_sign(1:PGEOM%max_orb,1:PGEOM%n_atom)
             do i=1,PGEOM%n_atom
               if(PGEOM%n_orbital(i) .eq. 0) then
                 if(PINPT%flag_report_geom) then
                   write(message,'(A,I4,A,I3,2x,10A7)')' ATOM',i,': ',PGEOM%n_orbital(i), PGEOM%c_orbital(1,i) ; write_msgi_file
                 endif
               elseif(PGEOM%n_orbital(i) .gt. 0) then
                 if(PINPT%flag_report_geom) then
                   write(message,'(A,I4,A,I3,2x,10A7)')' ATOM',i,': ',PGEOM%n_orbital(i), PGEOM%c_orbital(1:PGEOM%n_orbital(i),i) ; write_msgi_file
                   write(message,'(A,A20)'            )' SITE_IDX:   ',site_c_index_(i) ; write_msgi_file
                 endif
                 if(PINPT%flag_local_charge) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(PINPT%flag_report_geom) then
                     write(message,'(A,*(F10.4))')'   CHARGE:   ',local_charge_(i_dummy:i_dummy1) ; write_msgi_file
                   endif
                 endif

                 if(PINPT%flag_collinear) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(PINPT%flag_report_geom) then
                     write(message,'(A,*(F10.4))')'   MAGMOM:   ',local_moment_(1,i_dummy:i_dummy1) ; write_msgi_file
                   endif
                 elseif(PINPT%flag_noncollinear) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(flag_moment_cart)then
                     if(PINPT%flag_report_geom) then
                       write(message,'(A,*(3F7.3,2x))')'   MAGMOM: (Mx,My,Mz) ',(local_moment_(1:3,i_dummy2),i_dummy2=i_dummy,i_dummy1) ; write_msgi_file
                     endif
                   else
                     if(PINPT%flag_report_geom) then
                       write(message,'(A,*(3F7.3,2x))')'   MAGMOM: (M,theta,phi) ',(local_moment_(1:3,i_dummy2),i_dummy2=i_dummy,i_dummy1) ; write_msgi_file
                     endif
                   endif
                 endif

               endif
             enddo
           elseif(PGEOM%neig < 0)then
             write(message,'(A)')'  !! Check geometry input file. negative number of atomic orbitals ??'  ; write_msgi
           endif
         endif ! linecount

         if (  i_continue .ne. 0 ) cycle  ! skip empty line 

      enddo line

! call check_sanity(PINOT,PGEOM) !!!!! for future work
  if(PGEOM%flag_cartesian) then ! cartesian
    t_latt_inv = inv(PGEOM%a_latt)
    ii = 0
    allocate(PGEOM%o_coord(3,PGEOM%neig))
    allocate(PGEOM%o_coord_cart(3,PGEOM%neig))
    do i = 1, PGEOM%n_atom
      t_coord=PGEOM%a_coord(:,i)
      PGEOM%a_coord(1,i) = dot_product(t_latt_inv(1,:), t_coord(:)) 
      PGEOM%a_coord(2,i) = dot_product(t_latt_inv(2,:), t_coord(:))
      PGEOM%a_coord(3,i) = dot_product(t_latt_inv(3,:), t_coord(:))
      PGEOM%a_coord(:,i) = PGEOM%a_coord(:,i) - real(int(PGEOM%a_coord(:,i)))
      PGEOM%a_coord_cart(:,i) = PGEOM%a_coord(1,i) * PGEOM%a_latt(:,1) + &
                                PGEOM%a_coord(2,i) * PGEOM%a_latt(:,2) + &
                                PGEOM%a_coord(3,i) * PGEOM%a_latt(:,3)
      do iorb=1,PGEOM%n_orbital(i)
        ii = ii + 1
        PGEOM%o_coord(:,ii) = PGEOM%a_coord(:,i)
        PGEOM%o_coord_cart(:,ii) = PGEOM%o_coord(1,ii) * PGEOM%a_latt(:,1) + &
                                   PGEOM%o_coord(2,ii) * PGEOM%a_latt(:,2) + &
                                   PGEOM%o_coord(3,ii) * PGEOM%a_latt(:,3)
      enddo
    enddo
  else ! direct
    allocate(PGEOM%o_coord(3,PGEOM%neig))
    allocate(PGEOM%o_coord_cart(3,PGEOM%neig))
    ii = 0
    do i = 1, PGEOM%n_atom
      PGEOM%a_coord_cart(:,i) = PGEOM%a_coord(1,i) * PGEOM%a_latt(:,1) + &
                                PGEOM%a_coord(2,i) * PGEOM%a_latt(:,2) + &
                                PGEOM%a_coord(3,i) * PGEOM%a_latt(:,3)
      do iorb=1,PGEOM%n_orbital(i)
        ii = ii + 1
        PGEOM%o_coord(:,ii) = PGEOM%a_coord(:,i)
        PGEOM%o_coord_cart(:,ii) = PGEOM%o_coord(1,ii) * PGEOM%a_latt(:,1) + &
                                   PGEOM%o_coord(2,ii) * PGEOM%a_latt(:,2) + &
                                   PGEOM%o_coord(3,ii) * PGEOM%a_latt(:,3) 
      enddo
    enddo
  endif

  ! set efield_origin (if "efield" requested)
  if(NN_TABLE%flag_efield .and. NN_TABLE%flag_efield_frac) then
    NN_TABLE%efield_origin_cart(1:3) = NN_TABLE%efield_origin(1) * PGEOM%a_latt(1:3,1) + &
                                       NN_TABLE%efield_origin(2) * PGEOM%a_latt(1:3,2) + &
                                       NN_TABLE%efield_origin(3) * PGEOM%a_latt(1:3,3)
    write(message,'(A,3F12.6)')'EF_ORIGIN:  (in cartesian coord) ',NN_TABLE%efield_origin_cart(1:3)  ; write_msgi
  elseif(NN_TABLE%flag_efield .and. .not. NN_TABLE%flag_efield_frac .and. .not. NN_TABLE%flag_efield_cart) then
    NN_TABLE%efield_origin(1) = 0
    NN_TABLE%efield_origin(2) = 0
    NN_TABLE%efield_origin(3) = ( maxval(PGEOM%a_coord(3,:)) - minval(PGEOM%a_coord(3,:)) ) * 0.5d0

    NN_TABLE%efield_origin_cart(1:3) = NN_TABLE%efield_origin(1) * PGEOM%a_latt(1:3,1) + &
                                       NN_TABLE%efield_origin(2) * PGEOM%a_latt(1:3,2) + &
                                       NN_TABLE%efield_origin(3) * PGEOM%a_latt(1:3,3)
    write(message,'(A,3F12.6)')'EF_ORIGIN:  (in cartesian coord) ',NN_TABLE%efield_origin_cart(1:3)  ; write_msgi
  endif


  ! store site_index
  allocate(PGEOM%site_cindex(PGEOM%n_atom))
  allocate(PGEOM%flag_site_cindex(PGEOM%n_atom))
  PGEOM%site_cindex = site_c_index_
  PGEOM%flag_site_cindex = flag_site_c_index_

  ! store local moment
  allocate(PGEOM%local_moment(3,PGEOM%neig)) ! local net moment for each atomic orbital basis
  allocate(PGEOM%local_moment_rot(3,PGEOM%neig)) ! local moment rotated
  PGEOM%local_moment     = 0d0 ;   PGEOM%local_moment_rot = 0d0
  if(PINPT%flag_collinear) then
    PGEOM%local_moment(1,:) = local_moment_(1,1:PGEOM%neig)
  elseif(PINPT%flag_noncollinear) then
    PGEOM%local_moment(1,:) = local_moment_(1,1:PGEOM%neig)
    PGEOM%local_moment(2,:) = local_moment_(2,1:PGEOM%neig)
    PGEOM%local_moment(3,:) = local_moment_(3,1:PGEOM%neig)
    if(flag_moment_cart) then
      do i = 1, PGEOM%neig
        PGEOM%local_moment_rot(1:3,i) = PGEOM%local_moment(1:3,i)
      enddo
    else
      do i = 1, PGEOM%neig
        call rotate_vector( PGEOM%local_moment_rot(1:3,i), PGEOM%local_moment(1:3,i) )
      enddo
    endif
  endif

  ! store local charge
  allocate(PGEOM%local_charge(PGEOM%neig)) ! local charge for each atomic orbital basis
  PGEOM%local_charge = 0d0
  if(PINPT%flag_local_charge) then
    PGEOM%local_charge(:) = local_charge_(1:PGEOM%neig)
  endif

  if (linecount == 0) then
    write(message,*)'Attention - empty input file: ',trim(fname),' , ',func  ; write_msgi
#ifdef MPI
    call MPI_Abort(mpi_comm_earth, 0, mpierr)
#else
    stop
#endif
  endif
  close(pid_geom)

#ifdef SPGLIB
  ! get space group information by SPGLIB
  if(PINPT%flag_spglib) call get_symmetry_info(PGEOM)
#endif



return
endsubroutine
subroutine set_orbital_index(PGEOM, lmmax)
   use parameters, only: poscar
   use mpi_setup
   implicit none
   type (poscar)            :: PGEOM
   integer*4                   iatom, iorb, imatrix, lmmax
   integer*4, allocatable   :: orb_index(:)
   character*8                 orb
   integer*4                   ourjob(nprocs), ourjob_disp(0:nprocs-1), mpierr

   call mpi_job_distribution_chain(PGEOM%n_atom, nprocs, ourjob, ourjob_disp)
   allocate(PGEOM%orb_index(PGEOM%neig))
   allocate(orb_index(PGEOM%neig))
   PGEOM%orb_index = 0 ; orb_index = 0

   do iatom = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
     do iorb = 1, PGEOM%n_orbital(iatom)
       orb = adjustl(trim(PGEOM%c_orbital(iorb,iatom)))
       imatrix= sum( PGEOM%n_orbital(1:iatom) ) - PGEOM%n_orbital(iatom) + iorb
       
       if(lmmax .eq. 9) then
         select case (trim(orb))
           case ('s')
             PGEOM%orb_index(imatrix) = 1
  
           case ('px')
             PGEOM%orb_index(imatrix) = 2
           case ('py')
             PGEOM%orb_index(imatrix) = 3
           case ('pz')
             PGEOM%orb_index(imatrix) = 4
              
           case ('dz2')
             PGEOM%orb_index(imatrix) = 5
           case ('dx2')
             PGEOM%orb_index(imatrix) = 6
           case ('dxy')
             PGEOM%orb_index(imatrix) = 7
           case ('dxz')
             PGEOM%orb_index(imatrix) = 8
           case ('dyz')
             PGEOM%orb_index(imatrix) = 9
  
         endselect
       elseif(lmmax .eq. 3) then
         select case (trim(orb))
           case ('s')
             PGEOM%orb_index(imatrix) = 1

           case('px', 'py', 'pz') 
             PGEOM%orb_index(imatrix) = 2

           case('dz2', 'dx2', 'dxy', 'dxz', 'dyz')
             PGEOM%orb_index(imatrix) = 3
         endselect
       endif

     enddo ! iorb
   enddo ! iatom

#ifdef MPI
   call MPI_ALLREDUCE(PGEOM%orb_index,orb_index,size(PGEOM%orb_index),MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
   PGEOM%orb_index = orb_index
#endif

   return
endsubroutine

subroutine set_effective_nuclear_charge(PPRAM, PGEOM)
  use parameters,  only : poscar, params
  use mpi_setup
  use element_info, only : z_eff
  use print_io
  implicit none
  type (poscar)            :: PGEOM
  type (params)            :: PPRAM
  integer*4                   iatom, jorb
  character*8                 spec, spec_
  character*8                 orb
  integer*4                   ii, ispec, n, l, orb_n
  integer*4                   mpierr, ierror
  integer*4, allocatable ::   ispec_(:)
  integer*4, allocatable ::   l_quantum_(:,:)
  integer*4, allocatable ::   orb_n_quantum_(:,:)
  real*8,    allocatable ::   n_quantum_(:) ! just set to "real" for the convinience.
  real*8,    allocatable ::   z_eff_nuc_(:,:)
  character*40                dummy1, dummy2
  integer*4                   ourjob(nprocs)
  integer*4                   ourjob_disp(0:nprocs-1)
  call mpi_job_distribution_chain(PGEOM%n_atom, nprocs, ourjob, ourjob_disp)

  allocate(PGEOM%n_quantum(PGEOM%n_atom), &
                n_quantum_(PGEOM%n_atom) ) ! n=1, 2, ...?
  allocate(PGEOM%l_quantum(PGEOM%max_orb, PGEOM%n_atom), &
                l_quantum_(PGEOM%max_orb, PGEOM%n_atom)) ! s, p, d?? for each orb
  allocate(PGEOM%orb_n_quantum(PGEOM%max_orb, PGEOM%n_atom), &
                orb_n_quantum_(PGEOM%max_orb, PGEOM%n_atom)) ! 1s, 2p, 3d?? for each orb
  allocate(PGEOM%z_eff_nuc(PGEOM%max_orb, PGEOM%n_atom), &
                z_eff_nuc_(PGEOM%max_orb, PGEOM%n_atom)) 
  allocate(PGEOM%ispec(PGEOM%n_atom), &
                ispec_(PGEOM%n_atom))
  l_quantum_ = 0; orb_n_quantum_=0;z_eff_nuc_=0d0 ;n_quantum_ = 0d0 
  ispec_ = 0

lp1:do iatom = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
      spec_=adjustl(trim(PGEOM%c_spec( PGEOM%spec(iatom) )))
      call find_spec_const(PPRAM, PGEOM, spec_, spec)
      call find_spec(spec, ispec,n)
      if(ispec .eq. -1) then
        write(6,'(A,I0,3A)')'   !WARN! The ',iatom,'-th atom(',trim(spec_),') cannot be identified in the periodic table.'
        write(6,'(A)')     '          You should define the species information in the "SET CONSTRAINT" in your INCAR-TB.'
        write(6,'(A)')     '          For example, if your atom "Nb1" is actually meant to be Niobium (Nb),'
        write(6,'(A)')     '          then you can add following line in the "CONSTRAINT" section as follows:'
        write(6,'(A)')     '          SET CONSTRAINT '
        write(6,'(A)')     '            ... other constraint ...'
        write(6,'(A)')     '            ... ... '
        write(6,'(A)')     '            atom_Nb1 = spec_Nb '
        write(6,'(A)')     '          END CONSTRAINT '
        write(6,'(A)')     ' '
        write(6,'(A)')     '          Exit program...'
        stop
      endif
      ispec_(iatom) = ispec
      n_quantum_(iatom) = real(n)
      do jorb = 1, PGEOM%n_orbital(iatom)
        orb = adjustl(trim(PGEOM%c_orbital(jorb,iatom)))
        call find_lmom(ispec, orb, l, orb_n)
        l_quantum_(jorb,iatom)     = l
        orb_n_quantum_(jorb,iatom) = orb_n
        z_eff_nuc_(jorb,iatom) = z_eff(ispec, orb_n, l)
        if( z_eff_nuc_(jorb,iatom) .lt. 0.d0) then
          write(6,'(A,2A)')   '   ! WARN ! Be sure that the effective nuclear charge information for species "',trim(spec_),'" has not been stored'
          write(6,'(A)')      '            in "z_eff" function in element_info.f90 module. Please update the database based on'
          write(6,'(A)')      '            the information tabulated in following link: http://www.knowledgedoor.com'
          write(6,'(A)')      '            For example, if you want to incode "Helium" information, you can use following link:'
          write(6,'(A)')      '            http://www.knowledgedoor.com/2/elements_handbook/clementi-raimondi_effective_nuclear_charge_part_2.html#helium '
          write(6,'(A)')      '            Ref.1: Clementi, E., and D. L. Raimondi. "Atomic Screening Constants from SCF Functions.'
          write(6,'(A)')      '                   Journal of Chemical Physics, volume 38, number 11, 1963, pp. 2686–2689. doi:10.1063/1.1733573'
          write(6,'(A)')      '            Ref.2: Clementi, E., D. L. Raimondi, and W. P. Reinhardt. "Atomic Screening Constants from SCF Functions. II. Atoms with 37 to 86 Electrons.'
          write(6,'(A)')      '                   Journal of Chemical Physics, volume 47, number 4, 1967, pp. 1300–1307. doi:10.1063/1.1712084'
          write(6,'(A)')      ' '
          write(6,'(A)')      '            Exit program...'
          stop
        endif
      enddo
    enddo lp1

#ifdef MPI
  call MPI_ALLREDUCE(ispec_,         PGEOM%ispec,         size(PGEOM%ispec),         MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
  call MPI_ALLREDUCE(l_quantum_,     PGEOM%l_quantum,     size(PGEOM%l_quantum),     MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
  call MPI_ALLREDUCE(orb_n_quantum_, PGEOM%orb_n_quantum, size(PGEOM%orb_n_quantum), MPI_INTEGER4, MPI_SUM, mpi_comm_earth, mpierr)
  call MPI_ALLREDUCE(n_quantum_,     PGEOM%n_quantum,     size(PGEOM%n_quantum),     MPI_REAL8,    MPI_SUM, mpi_comm_earth, mpierr)
  call MPI_ALLREDUCE(z_eff_nuc_,     PGEOM%z_eff_nuc,     size(PGEOM%z_eff_nuc),     MPI_REAL8,    MPI_SUM, mpi_comm_earth, mpierr)
#else
  PGEOM%ispec         = ispec_
  PGEOM%l_quantum     = l_quantum_
  PGEOM%orb_n_quantum = orb_n_quantum_
  PGEOM%n_quantum     = n_quantum_
  PGEOM%z_eff_nuc     = z_eff_nuc_
#endif
  
   deallocate(ispec_)
   deallocate(l_quantum_)
   deallocate(orb_n_quantum_)
   deallocate(n_quantum_)
   deallocate(z_eff_nuc_)

   return
endsubroutine

subroutine find_lmom(ispec,orb, l, orb_n)
   use element_info
   implicit none
   integer*4    ispec, n, l, orb_n
   integer*4    i, lo
   character*8  orb
 
   !lo = len_trim(orb)

   select case( orb(1:1) )
     case('s')
       l = 1
     case('p')
       l = 2
     case('d')
       l = 3
     case('f')
       l = 4
   endselect

   orb_n = l_qnumb(l,ispec)

   return
endsubroutine
subroutine find_spec(spec,ispec,n)
  use element_info
  implicit none
  character*8   spec
  character*2   ele
  integer*4     n, ispec
  integer*4     i, le, ls, mpierr
  
   ispec = -1
  l1:do i = 1, 112
       ele = adjustl(trim(element(i)))
       le = len_trim(ele)
       ls = len_trim(spec)
       if( spec(1:ls)  .eq.  ele(1:le) ) then
         ispec = i
         n     = n_qnumb(i)
         exit l1   
       endif
     enddo l1
  
  return
endsubroutine

subroutine set_equiv_spec(PPRAM, PGEOM)
  use parameters, only: params, poscar
  implicit none
  type (params)            :: PPRAM
  type (poscar)            :: PGEOM
  integer*4                   i, ia, ja, ispec, jspec, spec_i, spec_j
  character*40                dummy1, dummy2
  character*8                 c_spec_i, c_spec_j

lp1:do i = 1, PPRAM%nparam_const
      if ( trim(PPRAM%c_const(2,i)) .eq. '=' ) then
        dummy1 = trim(PPRAM%c_const(3,i))
        dummy2 = trim(PPRAM%c_const(1,i))
        if(dummy1(1:4) .eq. 'atom' .and. dummy2(1:4) .eq. 'atom') then
          call strip_off(dummy1, c_spec_i, '_',' ',2)
          call strip_off(dummy2, c_spec_j, '_',' ',2)
          spec_i = 0
      lp2:do ispec = 1, PGEOM%n_spec
            if( trim(PGEOM%c_spec(ispec)) .eq. trim(c_spec_i) ) then
              spec_i = ispec
              exit lp2
            endif
          enddo lp2

          if(spec_i .ne. 0) then
            spec_j = 0
        lp3:do jspec = 1, PGEOM%n_spec
              if( trim(PGEOM%c_spec(jspec)) .eq. trim(c_spec_j)) then
                spec_j = jspec
                exit lp3
              endif
            enddo lp3
            if(spec_j .ne. 0) then
          lp4:do ja = 1, PGEOM%n_atom
                if( PGEOM%spec(ja) .eq. spec_j ) then
                  PGEOM%spec_equiv(ja) = spec_i
                endif
              enddo lp4
            else 
              cycle lp1
            endif ! if spec_j
          else
            cycle lp1
          endif ! if spec_i
        else
          cycle lp1
        endif ! if atom

      endif ! if '='
    enddo lp1

  return
endsubroutine

subroutine find_spec_const(PPRAM, PGEOM, spec, target_spec)
  use parameters, only: params, poscar
  implicit none
  type (params)            :: PPRAM
  type (poscar)            :: PGEOM
  integer*4                   i, ia, ja, ispec, jspec, spec_i, spec_j
  character*40                dummy1, dummy2
  character*8                 c_spec_i, c_spec_j
  character*8                 spec, target_spec

  target_spec = spec
lp1:do i = 1, PPRAM%nparam_const
        if ( trim(PPRAM%c_const(2,i)) .eq. '=' ) then
          dummy1 = trim(PPRAM%c_const(3,i))
          dummy2 = trim(PPRAM%c_const(1,i))
          if(dummy1(1:4) .eq. 'spec' .and. dummy2(1:4) .eq. 'atom') then
            call strip_off(dummy1, c_spec_i, '_',' ',2)
            call strip_off(dummy2, c_spec_j, '_',' ',2)
            if(trim(c_spec_j) .eq. trim(spec)) then
              target_spec = c_spec_i
              exit lp1
            endif
          endif
       endif
    enddo lp1

  return
endsubroutine
