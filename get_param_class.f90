subroutine get_nn_class(PGEOM, iatom,jatom,dij,onsite_tol, nn_class, r0)
   use parameters, only : poscar
   real*8      r0  ! reference distance
   real*8      dij, onsite_tol
   integer*4   n, iatom, jatom, li, lj
   integer*4   ii, nn, max_nn, nn_class, i_dummy
   character*8 c_atom_i, c_atom_j
   real*8      d_range
   real*8      dummy
   real*8      dist_nn(0:100)
   real*8      r0_nn(0:100)
   type(poscar) :: PGEOM

   c_atom_i=adjustl(trim(PGEOM%c_spec( PGEOM%spec(iatom) )))
   li=len_trim(c_atom_i)
   c_atom_j=adjustl(trim(PGEOM%c_spec( PGEOM%spec(jatom) )))
   lj=len_trim(c_atom_j)
   
   max_nn = 0
   dist_nn(0)=onsite_tol
   do n = 1, PGEOM%n_nn_type
     index_iatom=index( PGEOM%nn_pair(n),c_atom_i(1:li) )
     index_jatom=index( PGEOM%nn_pair(n),c_atom_j(1:lj) ,.TRUE. )

     if(index_iatom .eq. 0 .or. index_jatom .eq. 0) then
       cycle
     elseif(index_iatom .ge. 1 .and. index_jatom .ge. 1 .and. index_iatom .ne. index_jatom) then
       nn = len_trim(PGEOM%nn_pair(n)) - li - lj
     else
       cycle
     endif

     if( nn .ge. 1) then
       max_nn = max(nn, max_nn)
       dist_nn(nn) = PGEOM%nn_dist(n)
       r0_nn(nn)   = PGEOM%nn_r0(n)
     endif
   enddo

   nn_class = -9999
   do n = 0, max_nn
     if( n .eq. 0 .and. dij .le. dist_nn(n) )then  ! check onsite
       nn_class = 0
       r0       = 0
       exit
     elseif ( n .ge. 1 .and. dij .gt. dist_nn(n - 1) .and. dij .le. dist_nn(n) ) then
       nn_class = n
       r0       = r0_nn(n)
       exit
     endif
   enddo

return
endsubroutine
subroutine get_param_class(PGEOM,iorb,jorb,iatom,jatom,param_class)
   use parameters, only : poscar
   implicit none
   character*8 c_orb_i, c_orb_j
   real*8      r0, dij, onsite_tol
   integer*4   n,iorb,jorb,iatom,jatom,li, lj
   integer*4   index_iatom,index_jatom
   integer*4   ii, nn, max_nn, nn_class,i_dummy
   character*2 param_class
   character*2 c_orb_ij
   character*8 c_atom_i, c_atom_j
   real*8      d_range
   real*8      dummy
   real*8      dist_nn(0:100)
   real*8      r0_nn(0:100)
   type(poscar) :: PGEOM

   ! set hopping property: ss, sp, pp, dd ... etc. ?
   c_orb_i=PGEOM%c_orbital(iorb,iatom)
   c_orb_j=PGEOM%c_orbital(jorb,jatom)
   write(c_orb_ij,'(A1,A1)')c_orb_i(1:1),c_orb_j(1:1)
   
   select case ( c_orb_ij )

     case('ss')
       param_class = 'ss'
     case('sp', 'ps')
       param_class = 'sp'
     case('sd', 'ds')
       param_class = 'sd'
     case('sf', 'fs')
       param_class = 'sf'
     case('pp')
       param_class = 'pp'
     case('pd', 'dp')
       param_class = 'pd'
     case('pf', 'fp')
       param_class = 'pf'
     case('dd')
       param_class = 'dd'
     case('df', 'fd')
       param_class = 'df'
     case('ff')
       param_class = 'ff'
     case('cc')
       param_class = 'cc' ! user defined (for lattice model, BiSi110 example)
     case('xx')
       param_class = 'xx' ! user defined (for sk-type model, TaS2 example..)
   end select
 
return
endsubroutine
subroutine get_onsite_param_index(ionsite_param_index, PINPT, ci_orb, cj_orb, c_atom)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4   i, lio, ljo, la
   integer*4   ionsite_param_index
   character*8  ci_orb, cj_orb, c_atom
   character*20 ee_name

   lio = len_trim(ci_orb)
   ljo = len_trim(cj_orb)
   la = len_trim(c_atom)
   write(ee_name,*)'e_',ci_orb(1:lio),'_',c_atom(1:la)
   ionsite_param_index = 0

   if(ci_orb(1:lio) .eq. cj_orb(1:ljo)) then
    call get_param_index(PINPT, ee_name, ionsite_param_index)
   endif

   if ( ionsite_param_index .eq. -1) then
     write(6,'(A,A,A)')'    !WARNING! Onsite energy for ', adjustl(trim(ee_name)), &
                       ' is not asigned. Please check "SET TBPARAM" tag. Exit...'
     stop
   endif
   
return
endsubroutine
subroutine get_sk_index_set(index_sigma,index_pi,index_delta, &
                            index_sigma_scale,index_pi_scale,index_delta_scale, &
                            PINPT, param_class, nn_class, ci_atom, cj_atom, i_atom, j_atom, &
                            ci_site, cj_site, flag_use_site_cindex)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4                   nn_class
   integer*4                   i, lia, lja, lp
   integer*4                   i_atom, j_atom !species order as defined in GFILE
   integer*4                   index_sigma,index_pi,index_delta
   integer*4                   index_sigma_scale,index_pi_scale,index_delta_scale
   character*16                cij_pair
   character*2                 param_class
   character*8, intent(in) ::  ci_atom , cj_atom
   character*20,intent(in) ::  ci_site , cj_site 
   character*28            ::  ci_atom_, cj_atom_
   logical      flag_scale, flag_use_site_cindex

   ! initialize  
   index_sigma       =  0
   index_pi          =  0
   index_delta       =  0
   index_sigma_scale =  0
   index_pi_scale    =  0
   index_delta_scale =  0

   if(.not.flag_use_site_cindex) then

     flag_scale = .false.
     call get_param_name_index(PINPT, param_class, 'sigma', nn_class, ci_atom, cj_atom, flag_scale, index_sigma)
     call get_param_name_index(PINPT, param_class, 'pi'   , nn_class, ci_atom, cj_atom, flag_scale, index_pi   )
     call get_param_name_index(PINPT, param_class, 'delta', nn_class, ci_atom, cj_atom, flag_scale, index_delta)

     flag_scale = .true. 
     call get_param_name_index(PINPT, param_class, 'sigma', nn_class, ci_atom, cj_atom, flag_scale, index_sigma_scale)
     call get_param_name_index(PINPT, param_class, 'pi'   , nn_class, ci_atom, cj_atom, flag_scale, index_pi_scale   )
     call get_param_name_index(PINPT, param_class, 'delta', nn_class, ci_atom, cj_atom, flag_scale, index_delta_scale)

   elseif(flag_use_site_cindex) then

     write(ci_atom_,'(A,A)')trim(ci_atom),trim(ci_site)
     write(cj_atom_,'(A,A)')trim(cj_atom),trim(cj_site)

     flag_scale = .false.
     call get_param_name_index(PINPT, param_class, 'sigma', nn_class, ci_atom_, cj_atom_, flag_scale, index_sigma)
     call get_param_name_index(PINPT, param_class, 'pi'   , nn_class, ci_atom_, cj_atom_, flag_scale, index_pi   )
     call get_param_name_index(PINPT, param_class, 'delta', nn_class, ci_atom_, cj_atom_, flag_scale, index_delta)

     flag_scale = .true.
     call get_param_name_index(PINPT, param_class, 'sigma', nn_class, ci_atom_, cj_atom_, flag_scale, index_sigma_scale)
     call get_param_name_index(PINPT, param_class, 'pi'   , nn_class, ci_atom_, cj_atom_, flag_scale, index_pi_scale   )
     call get_param_name_index(PINPT, param_class, 'delta', nn_class, ci_atom_, cj_atom_, flag_scale, index_delta_scale)

   endif

return
endsubroutine
subroutine get_local_U_param_index(local_U_param_index, PINPT, nn_class, param_class, ci_atom)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4     nn_class
   integer*4     i, ii, liatom, lparam
   integer*4     local_U_param_index
   character*8   ci_atom
   character*2   param_class, param_class_
   character*40  local_U_param_name 


   !initialize
   local_U_param_index = 0

   liatom = len_trim(ci_atom)
   lparam = len_trim(param_class)
   param_class_ = adjustl(trim(param_class))

   if( nn_class .eq. 0) then
     if( adjustl(trim(param_class_)) .eq. 'ss' ) write(local_U_param_name ,  99)'local_U_s_',ci_atom(1:liatom)
     if( adjustl(trim(param_class_)) .eq. 'pp' ) write(local_U_param_name ,  99)'local_U_p_',ci_atom(1:liatom)
     if( adjustl(trim(param_class_)) .eq. 'dd' ) write(local_U_param_name ,  99)'local_U_d_',ci_atom(1:liatom)

     ! user defined param class and its corresponding parameter type. Here, in this case (TaS2), xp1, xp2, xp3 
     ! orbitals are used so that the local_U index follows local_U_'x'_ someting.
     if( adjustl(trim(param_class_)) .eq. 'xx' ) write(local_U_param_name ,  99)'local_U_x_',ci_atom(1:liatom)

     ! user defined param class and its corresponding parameter type. Here, in this case (BiSi110), cp1 orbital
     ! orbitals are used so that the local_U index follows local_U_'c'_ someting.
     if( adjustl(trim(param_class_)) .eq. 'cc' ) write(local_U_param_name ,  99)'local_U_c_',ci_atom(1:liatom)
   endif

   call get_param_index(PINPT, local_U_param_name, local_U_param_index)
99 format(A,A)   

return
endsubroutine

subroutine get_plus_U_param_index(plus_U_param_index, PINPT, nn_class, param_class, ci_atom)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4      nn_class
   integer*4      i, ii, liatom, lparam
   integer*4      plus_U_param_index
   character*8    ci_atom
   character*2   param_class, param_class_
   character*40  plus_U_param_name

   !initialize
   plus_U_param_index = 0
   
   liatom = len_trim(ci_atom)
   lparam = len_trim(param_class)
   param_class_ = adjustl(trim(param_class))

   if( nn_class .eq. 0 ) then
     if( adjustl(trim(param_class_)) .eq. 'ss' ) write(plus_U_param_name ,  99)'plus_U_s_',ci_atom(1:liatom)
     if( adjustl(trim(param_class_)) .eq. 'pp' ) write(plus_U_param_name ,  99)'plus_U_p_',ci_atom(1:liatom)
     if( adjustl(trim(param_class_)) .eq. 'dd' ) write(plus_U_param_name ,  99)'plus_U_d_',ci_atom(1:liatom)
   endif

   call get_param_index(PINPT, plus_U_param_name, plus_U_param_index)

99 format(A,A)

return
endsubroutine
subroutine get_stoner_I_param_index(stoner_I_param_index, PINPT, nn_class, param_class, ci_atom)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4     nn_class
   integer*4     i, ii, liatom, lparam
   integer*4     stoner_I_param_index
   character*8   ci_atom
   character*2   param_class, param_class_
   character*40  stoner_I_param_name


  !initialize
  stoner_I_param_index = 0

   liatom = len_trim(ci_atom)
   lparam = len_trim(param_class)
   param_class_ = adjustl(trim(param_class))

   if( nn_class .eq. 0) then
     if( adjustl(trim(param_class_)) .eq. 'ss' ) write(stoner_I_param_name ,  99)'stoner_I_s_',ci_atom(1:liatom)
     if( adjustl(trim(param_class_)) .eq. 'pp' ) write(stoner_I_param_name ,  99)'stoner_I_p_',ci_atom(1:liatom)
     if( adjustl(trim(param_class_)) .eq. 'dd' ) write(stoner_I_param_name ,  99)'stoner_I_d_',ci_atom(1:liatom)

     ! for TaS2, user defined xx-class
     if( adjustl(trim(param_class_)) .eq. 'xx' ) write(stoner_I_param_name ,  99)'stoner_I_x_',ci_atom(1:liatom)

     ! for Bi/Si110, user defined cc-class
     if( adjustl(trim(param_class_)) .eq. 'cc' ) write(stoner_I_param_name ,  99)'stoner_I_c_',ci_atom(1:liatom)

   endif

   call get_param_index(PINPT, stoner_I_param_name, stoner_I_param_index)

99 format(A,A)

return
endsubroutine
subroutine get_param_name(param_name, param_class, param_type, nn_class, ci_atom, cj_atom, flag_scale)
   implicit none
   integer*4    lia, lja, lp
   integer*4    nn_class
   character*2  param_class
   character*8  param_type
   character*40 param_name
   logical      flag_scale
   character(*), intent(in) ::  ci_atom, cj_atom

   lia = len_trim(ci_atom)
   lja = len_trim(cj_atom)
   lp  = len_trim(param_class)

   if(.not. flag_scale) then
     if(nn_class .lt. 10) then
       write(param_name,99)      param_class(1:lp),param_type(1:1),'_',nn_class,'_',ci_atom(1:lia),cj_atom(1:lja)
     elseif(nn_class .ge. 10) then
       write(param_name,98)      param_class(1:lp),param_type(1:1),'_',nn_class,'_',ci_atom(1:lia),cj_atom(1:lja)
     endif
   elseif(flag_scale) then
     if(nn_class .lt. 10) then
       write(param_name,89) 's_',param_class(1:lp),param_type(1:1),'_',nn_class,'_',ci_atom(1:lia),cj_atom(1:lja)
     elseif(nn_class .ge. 10) then
       write(param_name,88) 's_',param_class(1:lp),param_type(1:1),'_',nn_class,'_',ci_atom(1:lia),cj_atom(1:lja)
     endif
   endif

99 format(3A,I1,3A)
89 format(4A,I1,3A)
98 format(3A,I2,3A)
88 format(4A,I2,3A)

return
endsubroutine

subroutine get_param_name_index(PINPT, param_class, param_type, nn_class, ci_atom, cj_atom, &
                                flag_scale, param_index)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4      i_atempt
   integer*4      nn_class
   integer*4      param_index
   character(*),intent(in) ::    ci_atom, cj_atom
   character*2    param_class
   character*8    param_type
   character*40   param_name
   logical        flag_scale  

loop:do i_atempt = 0, 1
       if(i_atempt .eq. 0) call get_param_name(param_name, param_class, trim(param_type), nn_class, &
                                               ci_atom, cj_atom, flag_scale)
       if(i_atempt .eq. 1) call get_param_name(param_name, param_class, trim(param_type), nn_class, &
                                               cj_atom, ci_atom, flag_scale)
       call get_param_index(PINPT, param_name, param_index)
       if(param_index .gt. 0) exit loop 
     enddo loop
   return
endsubroutine
subroutine get_param_index(PINPT, param_name, param_index)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4      i
   integer*4      param_index
   character(*), intent(in) ::  param_name
   character*40   pname_file, pname 
   character*40   str2lowcase
   external    :: str2lowcase

   pname = adjustl(trim(param_name))
   pname = str2lowcase(pname)
   param_index = 0
   ! find parameter index with given parameter name
   do i = 1, PINPT%nparam
     pname_file=adjustl(trim(PINPT%param_name(i)))
     pname_file=str2lowcase(pname_file)
!    if( adjustl(trim(PINPT%param_name(i))) .eq. adjustl(trim(param_name)) ) then

     if( trim(pname_file) .eq. adjustl(trim(pname)) ) then
       if( nint(PINPT%param_const(1,i)) .ge. 1 ) then
         param_index = nint(PINPT%param_const(1,i), 4) ! set constraint condition: if equal to
!  write(6,*)"XXX ", pname, pname_file, i, param_index
       else
         param_index = i
       endif
       exit
     endif
   enddo

return
endsubroutine
subroutine get_param(PINPT, param_index, param)
   use parameters, only : incar
   implicit none
   type(incar) :: PINPT
   integer*4      param_index
   real*8         param


   if(param_index .ne. 0) then
     ! get parameter value with given parameter index
     if( nint(PINPT%param_const(4, param_index )) .ge. 1 ) then
       param = PINPT%param_const(5, param_index ) ! set constraint condition: if fixed
     else
       param = PINPT%param( param_index )
     endif

   elseif(param_index .eq. 0) then
     param = 0.d0
   endif

return
endsubroutine
