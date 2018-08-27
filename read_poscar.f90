subroutine read_poscar(PINPT,PGEOM,NN_TABLE)
  use parameters,  only : incar, poscar, hopping, pid_geom
! use inverse_mat, only : inv
  use do_math, only : rotate_vector, inv
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
  character*8, allocatable :: temp_orbital(:,:)
  character*20,allocatable :: site_c_index_(:)
  character*10                site_index
  character*20                locpot_index
  character(*), parameter  :: func = 'read_poscar'
  logical                     flag_skip, flag_read_moment, flag_moment_cart
! real*8                      rot_m(3)
  external                    nitems, str2lowcase
  type (incar )            :: PINPT  
  type (poscar)            :: PGEOM  
  type (hopping)           :: NN_TABLE

  fname        = PINPT%gfilenm
  PGEOM%n_spec = 0
  PGEOM%flag_selective = .false.
  flag_read_moment = .false.
  flag_moment_cart = .false.
  pos_index = 0

  if(myid .eq. 0) write(6,*)' '
  if(myid .eq. 0) write(6,*)'*- READING INPUT GEOMETRY FILE: ',trim(fname)
  open (pid_geom, FILE=fname,iostat=i_continue)
  linecount = 0
  ii = 0
line: do
        read(pid_geom,'(A)',iostat=i_continue) inputline
        if(i_continue<0) exit               ! end of file reached
        if(i_continue>0) then 
          if(myid .eq. 0) write(6,*)'Unknown error reading file:',trim(fname),func
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
           if(myid .eq. 0) write(6,'(A,A)')'   SYSTEM:  ', trim(PGEOM%system_name)
           cycle

        ! scaling factor
         elseif(linecount .eq. 2) then
           read(inputline,*,iostat=i_continue) PGEOM%a_scale
           if(myid .eq. 0) write(6,'(A,F15.8)')'  A_SCALE:  ',PGEOM%a_scale
           cycle

        ! lattice parameter
         elseif(linecount .eq. 3 ) then
           backspace(pid_geom)
           do i=1,3
             read(pid_geom,'(A)',iostat=i_continue) inputline
             read(inputline,*,iostat=i_continue) PGEOM%a_latt(1:3,i)
             if(myid .eq. 0) write(6,'(A,i1,A,3F15.8)')'  A_LATT',i,':  ',PGEOM%a_latt(1:3,i)
           enddo
           call get_reci(PGEOM%b_latt(:,1), PGEOM%b_latt(:,2), PGEOM%b_latt(:,3), &
                         PGEOM%a_latt(:,1), PGEOM%a_latt(:,2), PGEOM%a_latt(:,3))
           do i=1,3
             if(myid .eq. 0) write(6,'(A,i1,A,3F15.8)')'  B_RECI',i,':  ',PGEOM%b_latt(1:3,i)
           enddo
           linecount = linecount + 2
           cycle

        ! species name and number of atoms
         elseif(linecount .eq. 6 ) then
           PGEOM%n_spec=nitems(inputline)
           allocate( PGEOM%c_spec(PGEOM%n_spec), PGEOM%i_spec(PGEOM%n_spec) )
           if(myid .eq. 0) write(6,'(A,i8)',ADVANCE='YES')'   N_SPEC:',PGEOM%n_spec
           read(inputline,*,iostat=i_continue) PGEOM%c_spec(1:PGEOM%n_spec)
           read(pid_geom,'(A)',iostat=i_continue) inputline
           linecount = linecount + 1
           call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle
           read(inputline,*,iostat=i_continue) PGEOM%i_spec(1:PGEOM%n_spec)
           PGEOM%n_atom=sum ( PGEOM%i_spec(1:PGEOM%n_spec) )

           allocate( PGEOM%spec(PGEOM%n_atom) )
           do i=1,PGEOM%n_spec
             PGEOM%spec( sum(PGEOM%i_spec(1:i)) -PGEOM%i_spec(i)+1 : sum(PGEOM%i_spec(1:i)) ) = i
           enddo

           if(myid .eq. 0) write(6,'(A,i8)')'   N_ATOM:',PGEOM%n_atom
           allocate( PGEOM%a_coord(3,PGEOM%n_atom), &
                     PGEOM%n_orbital(PGEOM%n_atom), &
                     local_charge_(PGEOM%n_atom*max_orb_temp), &
                     local_moment_(3,PGEOM%n_atom*max_orb_temp), &
                     site_c_index_(PGEOM%n_atom), &
                     temp_orbital(max_orb_temp, PGEOM%n_atom) )
                     local_charge_ = 0d0 ! initialize as zero
                     local_moment_ = 0d0 ! initialize as zero
           if(myid .eq. 0) write(6,'(A)',ADVANCE='NO')' '
           do i=1,PGEOM%n_spec
             if(i .eq. 1)then
               if(myid .eq. 0) write(6,'(A,I2,A,A4,1x,i8)')'  SPEC',i,':',trim(PGEOM%c_spec(i)),PGEOM%i_spec(i)
             else
               if(myid .eq. 0) write(6,'(A,I2,A,A4,1x,i8)')'   SPEC',i,':',trim(PGEOM%c_spec(i)),PGEOM%i_spec(i)
             endif
           enddo

        ! constraint and coordinate type
         elseif(linecount .eq. 8 ) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'S' .or. desc_str(1:1) .eq. 's') then 
             PGEOM%flag_selective = .true.
             if(myid .eq. 0) write(6,'(A)')' L_CONSTR:  .TRUE.'
           elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then 
             PGEOM%flag_selective = .false.
             PGEOM%flag_direct=.true.
             PGEOM%flag_cartesian=.false.
             linecount = linecount + 1
             if(myid .eq. 0) write(6,'(A)')' L_CONSTR:  .TRUE.'
             if(myid .eq. 0) write(6,'(A)')' C_CRDTYP:  DIRECT'
           elseif(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
                  desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then 
             PGEOM%flag_selective = .false.
             PGEOM%flag_direct=.false.
             PGEOM%flag_cartesian=.true.
             linecount = linecount + 1
             if(myid .eq. 0) write(6,'(A)')' L_CONSTR:  .FALSE.'
             if(myid .eq. 0) write(6,'(A)')' C_CRDTYP:  CARTESIAN'
           endif
         elseif(linecount .eq. 9 ) then
           read(inputline,*,iostat=i_continue) desc_str
           if(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
              desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then 
             PGEOM%flag_direct=.false.
             PGEOM%flag_cartesian=.true.
             if(myid .eq. 0) write(6,'(A)')' C_CRDTYP:  CARTESIAN'
           elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then 
             PGEOM%flag_direct=.true.
             PGEOM%flag_cartesian=.false.
             if(myid .eq. 0) write(6,'(A)')' C_CRDTYP:  DIRECT'
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
                 if(myid .eq. 0) write(6,'(A,I6)')'  !WARNING! Charge setting is inproper. Number of items should be same as N_ORBITAL(i). iatom=',i
                 if(myid .eq. 0) write(6,'(A)')   '  !WARNING! Please check GFILE. Exit program...'
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
                   if(myid .eq. 0) write(6,'(A)')'  !WARNING! Moment setting is inproper. Number of items should be same as N_ORBITAL(i) in the collinear setting.'
                   if(myid .eq. 0) write(6,'(A,I6)')'  !WARNING! Please check GFILE. Exit program...  iatom=',i
                   stop
                 else 
                   read(dummy,*)local_moment_(1,ii:ii+PGEOM%n_orbital(i)-1)
                 endif
               endif
               if(PINPT%flag_noncollinear) then
                 if(i_dummy .ne. PGEOM%n_orbital(i)*3) then
                   if(myid .eq. 0) write(6,'(A)')'  !WARNING! Moment setting is inproper. Number of items should be same as N_ORBITAL(i)*3 in the non-collinear setting.'
                   if(myid .eq. 0) write(6,'(A,I6)')'  !WARNING! N_ORBITAL(i) * 3 = number of items (M, theta, phi) or (Mx,My,Mz). Please check GFILE. Exit program...  iatom=',i
                   stop
                 else
                   read(dummy,*)((local_moment_(j,i_dummy),j=1,3),i_dummy = ii, ii+PGEOM%n_orbital(i)-1)
                 endif
               endif ! moment
             endif

             if(pos_index(3) .ne. 0) then ! site_index
               call strip_off(trim(inputline), dummy, trim(site_index), '', 2)
               i_dummy = nitems(dummy)
               if(i_dummy .ne. 1) then
                 if(myid .eq. 0) write(6,'(A)')'  !WARNING! Site_index setting is inproper. Number of items should be one and the data type should be character(20)'
                 if(myid .eq. 0) write(6,'(A,I6)')'  !WARNING! Please check GFILE. Exit program...  iatom=',i
                 stop
               else
                 read(dummy,*)site_c_index_(i)
               endif
             else
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
                 read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i), desc_str, desc_str, desc_str, temp_orbital(1:PGEOM%n_orbital(i),i)
               else
                 read(inputline,*,iostat=i_continue) PGEOM%a_coord(1:3,i), temp_orbital(1:PGEOM%n_orbital(i),i)
               endif
             endif

           enddo !n_atom

           PGEOM%neig=sum(PGEOM%n_orbital(1:PGEOM%n_atom))
           PGEOM%neig_total = PGEOM%neig * PINPT%ispin
           PGEOM%nbasis = PGEOM%neig 

           if(PGEOM%neig  == 0) then
             if(myid .eq. 0) write(6,'(A)')'  !! Check geometry input file. atomic orbital is not asigned!'
           elseif(PGEOM%neig >= 1) then
             if(myid .eq. 0) write(6,'(A,i8)')'  N_ORBIT:',PGEOM%neig
             PGEOM%max_orb=maxval( PGEOM%n_orbital(:) )
             allocate( PGEOM%c_orbital(PGEOM%max_orb,PGEOM%n_atom) ) ; PGEOM%c_orbital = '_na_'
             PGEOM%c_orbital(1:PGEOM%max_orb,1:PGEOM%n_atom) = temp_orbital(1:PGEOM%max_orb,1:PGEOM%n_atom)
             do i=1,PGEOM%n_atom
               if(PGEOM%n_orbital(i) .eq. 0) then
                 if(myid .eq. 0) write(6,'(A,I4,A,I3,2x,10A7)')' ATOM',i,': ',PGEOM%n_orbital(i), &
                                                              PGEOM%c_orbital(1,i)
               elseif(PGEOM%n_orbital(i) .gt. 0) then
                 if(myid .eq. 0) write(6,'(A,I4,A,I3,2x,10A7)')' ATOM',i,': ',PGEOM%n_orbital(i), &
                                                              PGEOM%c_orbital(1:PGEOM%n_orbital(i),i)
                 if(myid .eq. 0) write(6,'(A,A20)'            )' SITE_IDX:   ',site_c_index_(i)

                 if(PINPT%flag_local_charge) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(myid .eq. 0) write(6,'(A)',ADVANCE='NO')'   CHARGE:   '
                   do i_dummy2 = i_dummy, i_dummy1
                     if(myid .eq. 0) write(6,'(F7.3)',ADVANCE='NO')local_charge_(i_dummy2)
                   enddo
                   if(myid .eq. 0) write(6,*)''
                 endif

                 if(PINPT%flag_collinear) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(myid .eq. 0) write(6,'(A)',ADVANCE='NO')'   MAGMOM:   '
                   do i_dummy2 = i_dummy, i_dummy1
                     if(myid .eq. 0) write(6,'(F7.3)',ADVANCE='NO')local_moment_(1,i_dummy2)
                   enddo
                   if(myid .eq. 0) write(6,*)'' ! # note (m_i) for each orbital, m_i = magnetization of atom i'
                 elseif(PINPT%flag_noncollinear) then
                   i_dummy = sum(PGEOM%n_orbital(1:i)) - PGEOM%n_orbital(i) + 1
                   i_dummy1= sum(PGEOM%n_orbital(1:i))
                   if(flag_moment_cart)then
                     if(myid .eq. 0) write(6,'(A)',ADVANCE='NO')'   MAGMOM: (Mx,My,Mz) '
                   else
                     if(myid .eq. 0) write(6,'(A)',ADVANCE='NO')'   MAGMOM: (M,theta,phi) '
                   endif
                   do i_dummy2 = i_dummy, i_dummy1
                     if(myid .eq. 0) write(6,'(3F7.3,2x)',ADVANCE='NO') local_moment_(1:3,i_dummy2)
                   enddo
                   if(myid .eq. 0) write(6,*)'' ! # note (m_i,theta,phi) for each orbital, m_i = magnetization of atom i'
                 endif

               endif
             enddo
           elseif(PGEOM%neig < 0)then
             if(myid .eq. 0) write(6,'(A)')'  !! Check geometry input file. negative number of atomic orbitals ??'
           endif
         endif ! linecount

         if (  i_continue .ne. 0 ) cycle  ! skip empty line 

      enddo line
  
! call check_sanity(PINOT,PGEOM) !!!!! for future work
  if(PGEOM%flag_cartesian) then
    t_latt_inv = inv(PGEOM%a_latt)
    ii = 0
    allocate(PGEOM%o_coord(3,PGEOM%neig))
    allocate(PGEOM%o_coord_cart(3,PGEOM%neig))
    do i = 1, PGEOM%n_atom
      t_coord=PGEOM%a_coord(:,i)
      PGEOM%a_coord(1,i) = dot_product(t_latt_inv(1,:), t_coord(:)) 
      PGEOM%a_coord(2,i) = dot_product(t_latt_inv(2,:), t_coord(:))
      PGEOM%a_coord(3,i) = dot_product(t_latt_inv(3,:), t_coord(:))
      PGEOM%a_coord(:,i) = PGEOM%a_coord(:,i) - int(PGEOM%a_coord(:,i))
      do iorb=1,PGEOM%n_orbital(i)
        ii = ii + 1
        PGEOM%o_coord(:,ii) = PGEOM%a_coord(:,i)
        PGEOM%o_coord_cart(:,ii) = PGEOM%o_coord(1,ii) * PGEOM%a_latt(:,1) + &
                                   PGEOM%o_coord(2,ii) * PGEOM%a_latt(:,2) + &
                                   PGEOM%o_coord(3,ii) * PGEOM%a_latt(:,3)
      enddo
    enddo
  else
    allocate(PGEOM%o_coord(3,PGEOM%neig))
    allocate(PGEOM%o_coord_cart(3,PGEOM%neig))
    ii = 0
    do i = 1, PGEOM%n_atom
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
  if(PINPT%flag_efield .and. PINPT%flag_efield_frac) then
    PINPT%efield_origin_cart(1:3) = PINPT%efield_origin(1) * PGEOM%a_latt(1:3,1) + &
                                    PINPT%efield_origin(2) * PGEOM%a_latt(1:3,2) + &
                                    PINPT%efield_origin(3) * PGEOM%a_latt(1:3,3)
    if(myid .eq. 0) write(6,'(A,3F12.6)')'EF_ORIGIN:  (in cartesian coord) ',PINPT%efield_origin_cart(1:3)
  elseif(PINPT%flag_efield .and. .not. PINPT%flag_efield_frac .and. .not. PINPT%flag_efield_cart) then
    PINPT%efield_origin(1) = 0
    PINPT%efield_origin(2) = 0
    PINPT%efield_origin(3) = ( maxval(PGEOM%a_coord(3,:)) - minval(PGEOM%a_coord(3,:)) ) * 0.5d0

    PINPT%efield_origin_cart(1:3) = PINPT%efield_origin(1) * PGEOM%a_latt(1:3,1) + &
                                    PINPT%efield_origin(2) * PGEOM%a_latt(1:3,2) + &
                                    PINPT%efield_origin(3) * PGEOM%a_latt(1:3,3)
    if(myid .eq. 0) write(6,'(A,3F12.6)')'EF_ORIGIN:  (in cartesian coord) ',PINPT%efield_origin_cart(1:3)
  endif

  ! store local moment
  allocate(NN_TABLE%local_moment(3,PGEOM%neig)) ! local net moment for each atomic orbital basis
  allocate(NN_TABLE%local_moment_rot(3,PGEOM%neig)) ! local moment rotated
  NN_TABLE%local_moment     = 0d0
  NN_TABLE%local_moment_rot = 0d0
  if(PINPT%flag_collinear) then
    NN_TABLE%local_moment(1,:) = local_moment_(1,1:PGEOM%neig)
  elseif(PINPT%flag_noncollinear) then
    NN_TABLE%local_moment(1,:) = local_moment_(1,1:PGEOM%neig)
    NN_TABLE%local_moment(2,:) = local_moment_(2,1:PGEOM%neig)
    NN_TABLE%local_moment(3,:) = local_moment_(3,1:PGEOM%neig)
    if(flag_moment_cart) then
      do i = 1, PGEOM%neig
        NN_TABLE%local_moment_rot(1:3,i) = NN_TABLE%local_moment(1:3,i)
      enddo
    else
      do i = 1, PGEOM%neig
        call rotate_vector( NN_TABLE%local_moment_rot(1:3,i), NN_TABLE%local_moment(1:3,i) )
      enddo
    endif
  endif

  ! store local charge
  allocate(NN_TABLE%local_charge(PGEOM%neig)) ! local charge for each atomic orbital basis
  NN_TABLE%local_charge = 0d0
  if(PINPT%flag_local_charge) then
    NN_TABLE%local_charge(:) = local_charge_(1:PGEOM%neig)
  endif

  ! store site_index
  allocate(NN_TABLE%site_cindex(PGEOM%n_atom)) ! site_index for each atomic site
  NN_TABLE%site_cindex = site_c_index_

  if (linecount == 0) then
    if(myid .eq. 0) write(6,*)'Attention - empty input file: ',trim(fname),' , ',func
#ifdef MPI
    call MPI_Abort(mpi_comm_earth, 0, mpierr)
#else
    stop
#endif
  endif
  close(pid_geom)

  if(myid .eq. 0) write(6,*)'*- END READING GEOMETRY FILE ---------------------'
  if(myid .eq. 0) write(6,*)' '

return
endsubroutine
