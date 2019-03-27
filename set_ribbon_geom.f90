#include "alias.inc"
subroutine set_ribbon_geom(PINPT)
   use parameters 
   use mpi_setup
   use do_math, only : inv
   implicit none
   integer*4                 i,i_continue, linecount
   integer*4                 ispec, nspec, natom_spec(100),n_atom_unit, iatom, ld
   integer*4                 ix,iy,iz
   real*8                    t_latt_inv(3,3), t_latt(3,3)
   character*132             filenm_unitcell,filenm_ribbon
   character*264             inputline, c_dummy
   character*40              desc_str
   character(*), parameter  :: func = 'set_ribbon_geom'
   real*8                    a_latt(3,3),a_latt_ribbon(3,3)
   real*8                    enorm
   real*8                    coord_unit(3), coord_ribbon(3)
   real*8                    coord_unit_(3)
   integer*8                 i_dummy, nitems
   logical                   flag_skip
   logical                   flag_selective, flag_direct, flag_cartesian
   external                  enorm, nitems
   type(incar)           ::  PINPT   
   integer*4                 mpierr

   filenm_unitcell = PINPT%gfilenm 
   write(filenm_ribbon,*)trim(filenm_unitcell),'-ribbon'
   if_main write(6,'(A,A,A,A)')' GEOM_FNM: ',trim(PINPT%gfilenm),' ==> ',trim(filenm_ribbon)
   PINPT%gfilenm = filenm_ribbon

!  pid_geom_ribbon = pid_geom + 1
   
   open (pid_geom, FILE=filenm_unitcell,iostat=i_continue, status='old')
   open (pid_geom_ribbon, FILE=filenm_ribbon,iostat=i_continue, status='unknown')

   do 
     read(pid_geom,'(A)',iostat=i_continue) inputline
     if(i_continue<0) exit               ! end of file reached
     if(i_continue>0) write(6,*)'Unknown error reading file:',trim(filenm_unitcell),func
     linecount = linecount + 1
     call check_comment(inputline,linecount,i,flag_skip) ; if (flag_skip ) cycle


     ! head
      if(linecount .le. 2) then
        write(pid_geom_ribbon,'(A)')trim(inputline)
        cycle

     ! lattice parameter
      elseif(linecount .eq. 3) then
        backspace(pid_geom)
        do i=1,3
          read(pid_geom,'(A)',iostat=i_continue) inputline
          read(inputline,*,iostat=i_continue) a_latt(1:3,i)
          write(pid_geom_ribbon,'(3F20.12)') a_latt(1:3,i) * ( PINPT%ribbon_nslab(i) + PINPT%ribbon_vacuum(i)/enorm(3,a_latt(1:3,i)) )
          a_latt_ribbon(1:3,i) = a_latt(1:3,i) * ( PINPT%ribbon_nslab(i) + PINPT%ribbon_vacuum(i)/enorm(3,a_latt(1:3,i)) )
        enddo
        t_latt_inv = inv(a_latt)
        linecount = linecount + 2
        cycle
      elseif(linecount .eq. 6 ) then
        write(pid_geom_ribbon,'(A)')trim(inputline)
        nspec = nitems(trim(inputline))
        cycle

      elseif(linecount .eq. 7 ) then
        read(inputline,*,iostat=i_continue) natom_spec(1:nspec)
        write(pid_geom_ribbon,'(*(I8))') natom_spec(1:nspec) * PINPT%ribbon_nslab(1) * PINPT%ribbon_nslab(2) * PINPT%ribbon_nslab(3)
        n_atom_unit = sum(natom_spec(1:nspec))

     ! constraint and coordinate type
      elseif(linecount .eq. 8 ) then
        read(inputline,*,iostat=i_continue) desc_str
          if(desc_str(1:1) .eq. 'S' .or. desc_str(1:1) .eq. 's') then
            flag_selective = .true.
            write(pid_geom_ribbon,'(A)')'Selective'
          elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then
            flag_selective = .false.
            flag_direct=.true.
            flag_cartesian=.false.
            linecount = linecount + 1
            write(pid_geom_ribbon,'(A)')'Direct'
          elseif(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
                 desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then
            flag_selective = .false.
            flag_direct=.false.
            flag_cartesian=.true.
            linecount = linecount + 1
!           write(pid_geom_ribbon,'(A)')'Cartesian'
            write(pid_geom_ribbon,'(A)')'Direct' ! we enforce the output to be 'direct' coordinates
          endif
      elseif(linecount .eq. 9 ) then
        read(inputline,*,iostat=i_continue) desc_str
        if(desc_str(1:1) .eq. 'C' .or. desc_str(1:1) .eq. 'c' .or. &
           desc_str(1:1) .eq. 'K' .or. desc_str(1:1) .eq. 'k') then
          flag_direct=.false.
          flag_cartesian=.true.
!         write(pid_geom_ribbon,'(A)')'Cartesian'
          write(pid_geom_ribbon,'(A)')'Direct' ! we enforce the output to be 'direct' coordinates
        elseif(desc_str(1:1) .eq. 'D' .or. desc_str(1:1) .eq. 'd') then
          flag_direct=.true.
          flag_cartesian=.false.
          write(pid_geom_ribbon,'(A)')'Direct'
        endif

      ! atomic coordinate & atomic orbital information
      elseif(linecount .eq. 10 ) then
        backspace(pid_geom)

        do ispec = 1, nspec
          do iatom = 1, natom_spec(ispec)
            read(pid_geom,'(A)',iostat=i_continue) inputline
            i_dummy = nitems(inputline)
            !read coordinate
            if(i_dummy .ge. 4) then 
              read(inputline,*) coord_unit(1:3), desc_str
              if(flag_cartesian) then
!               t_latt = transpose(a_latt)
                coord_unit_ = coord_unit
                coord_unit(1) = dot_product(t_latt_inv(1,:), coord_unit_(:))
                coord_unit(2) = dot_product(t_latt_inv(2,:), coord_unit_(:))
                coord_unit(3) = dot_product(t_latt_inv(3,:), coord_unit_(:))
                coord_unit(:) = coord_unit(:) - int(coord_unit(:))

              endif
!             do i = 1,3
!               if(coord_unit(i) .lt. 0) coord_unit(i) = coord_unit(i) + 1d0
!             enddo
              ld=len_trim(desc_str)
            endif
            if(i_dummy .eq. 3) then 
              read(inputline,*) coord_unit(1:3)
              if(flag_cartesian) then
!               t_latt_inv = inv(a_latt)
                coord_unit_ = coord_unit
                coord_unit(1) = dot_product(t_latt_inv(1,:), coord_unit_(:))
                coord_unit(2) = dot_product(t_latt_inv(2,:), coord_unit_(:))
                coord_unit(3) = dot_product(t_latt_inv(3,:), coord_unit_(:))
                coord_unit(:) = coord_unit(:) - int(coord_unit(:))
              endif
!             do i = 1,3
!               if(coord_unit(i) .lt. 0) coord_unit(i) = coord_unit(i) + 1d0
!             enddo
            endif

            !set ribbon coordinate
            if(i_dummy .ge. 5) call strip_off(trim(inputline), c_dummy, desc_str(1:ld), ' ', 2)
            if(i_dummy .eq. 3) then 
              do iz = 1, PINPT%ribbon_nslab(3)
                do iy = 1, PINPT%ribbon_nslab(2)
                  do ix = 1, PINPT%ribbon_nslab(1)
                    coord_ribbon(1) = coord_unit(1) / real(PINPT%ribbon_nslab(1))+ (ix - 1)/real(PINPT%ribbon_nslab(1))
                    coord_ribbon(2) = coord_unit(2) / real(PINPT%ribbon_nslab(2))+ (iy - 1)/real(PINPT%ribbon_nslab(2))
                    coord_ribbon(3) = coord_unit(3) / real(PINPT%ribbon_nslab(3))+ (iz - 1)/real(PINPT%ribbon_nslab(3))
                    coord_ribbon(1) = coord_ribbon(1) - coord_ribbon(1) * PINPT%ribbon_vacuum(1)/enorm(3,a_latt_ribbon(1:3,1))
                    coord_ribbon(2) = coord_ribbon(2) - coord_ribbon(2) * PINPT%ribbon_vacuum(2)/enorm(3,a_latt_ribbon(1:3,2))
                    coord_ribbon(3) = coord_ribbon(3) - coord_ribbon(3) * PINPT%ribbon_vacuum(3)/enorm(3,a_latt_ribbon(1:3,3))
                    write(pid_geom_ribbon,'(3F20.12)') coord_ribbon(1:3) 
                  enddo
                enddo
              enddo
            elseif(i_dummy .eq. 4) then
              do iz = 1, PINPT%ribbon_nslab(3)
                do iy = 1, PINPT%ribbon_nslab(2)
                  do ix = 1, PINPT%ribbon_nslab(1)
                    coord_ribbon(1) = coord_unit(1) / real(PINPT%ribbon_nslab(1))+ (ix - 1)/real(PINPT%ribbon_nslab(1))
                    coord_ribbon(2) = coord_unit(2) / real(PINPT%ribbon_nslab(2))+ (iy - 1)/real(PINPT%ribbon_nslab(2))
                    coord_ribbon(3) = coord_unit(3) / real(PINPT%ribbon_nslab(3))+ (iz - 1)/real(PINPT%ribbon_nslab(3))
                    coord_ribbon(1) = coord_ribbon(1) - coord_ribbon(1) * PINPT%ribbon_vacuum(1)/enorm(3,a_latt_ribbon(1:3,1))
                    coord_ribbon(2) = coord_ribbon(2) - coord_ribbon(2) * PINPT%ribbon_vacuum(2)/enorm(3,a_latt_ribbon(1:3,2))
                    coord_ribbon(3) = coord_ribbon(3) - coord_ribbon(3) * PINPT%ribbon_vacuum(3)/enorm(3,a_latt_ribbon(1:3,3))
                    write(pid_geom_ribbon,'(3F20.12,2x,A)') coord_ribbon(1:3), trim(desc_str)
                  enddo 
                enddo
              enddo
            elseif(i_dummy .ge. 5) then
              do iz = 1, PINPT%ribbon_nslab(3)
                do iy = 1, PINPT%ribbon_nslab(2)
                  do ix = 1, PINPT%ribbon_nslab(1)
                    coord_ribbon(1) = coord_unit(1) / real(PINPT%ribbon_nslab(1))+ (ix - 1)/real(PINPT%ribbon_nslab(1))
                    coord_ribbon(2) = coord_unit(2) / real(PINPT%ribbon_nslab(2))+ (iy - 1)/real(PINPT%ribbon_nslab(2))
                    coord_ribbon(3) = coord_unit(3) / real(PINPT%ribbon_nslab(3))+ (iz - 1)/real(PINPT%ribbon_nslab(3))
                    coord_ribbon(1) = coord_ribbon(1) - coord_ribbon(1) * PINPT%ribbon_vacuum(1)/enorm(3,a_latt_ribbon(1:3,1))
                    coord_ribbon(2) = coord_ribbon(2) - coord_ribbon(2) * PINPT%ribbon_vacuum(2)/enorm(3,a_latt_ribbon(1:3,2))
                    coord_ribbon(3) = coord_ribbon(3) - coord_ribbon(3) * PINPT%ribbon_vacuum(3)/enorm(3,a_latt_ribbon(1:3,3))
                    write(pid_geom_ribbon,'(3F20.12,2x,A,2x,A)') coord_ribbon(1:3), trim(desc_str), trim(c_dummy)
                  enddo 
                enddo
              enddo
            endif
          enddo
        enddo
      endif
   enddo

   close(pid_geom)
   close(pid_geom_ribbon)

   if(PINPT%flag_print_only_ribbon_geom) then
     if_main write(6,'(A,A,A)')'  !WARN! PRINT_ONLY_RIBBON_GEOM requested..'
     if_main write(6,'(A,A,A)')'  !WARN! The ribbon geometry will be wrote down in /GFILE/-ribbon'
     if_main write(6,'(A,A,A)')'  !WARN! Exit program...'
     kill_job
   endif
return
endsubroutine
