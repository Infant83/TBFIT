subroutine get_site_number(PGEOM, site_cindex, range_dummy, size_range_dummy)
   use parameters, only : poscar
   implicit none
   type(poscar)  :: PGEOM
   integer*4        i, j, pos
   integer*4        lsite
   integer*4        range_dummy(PGEOM%n_atom), size_range_dummy
   character*20     site_cindex, site_name
   character*20     site_dummy
   logical          flag_all

   flag_all = .false.
   range_dummy = 0
   pos = 0

   if(index(site_cindex,'*') .gt. 0) then 
     flag_all = .true.
     pos = index(site_cindex,'*') 
   endif

   if(flag_all) then
     site_name = site_cindex(1:pos-1)
     lsite=len_trim(site_name)
     size_range_dummy = 0

     do i = 1, PGEOM%n_atom
       site_dummy = trim(PGEOM%site_cindex(i))
       if(site_dummy(1:lsite) .eq. site_name) then
         size_range_dummy = size_range_dummy + 1
         range_dummy(size_range_dummy) = i
       endif
     enddo
 
   elseif(.not. flag_all) then
     size_range_dummy = 0
     do i = 1, PGEOM%n_atom
       site_dummy = trim(PGEOM%site_cindex(i))
       if(trim(site_dummy) .eq. trim(site_cindex)) then
         size_range_dummy = size_range_dummy + 1
         range_dummy(size_range_dummy) = i
       endif
     enddo
 
   endif

return
endsubroutine
