module element_info

! character*2  element(112)
! integer*4    n_qnumb(112)
 !character*2, dimension(112), parameter :: element = &
  character*2, parameter::  element(112) = &
   (/'H '   ,&
     'He'  ,&
     'Li'  ,&
     'Be'  ,&
     'B '   ,&
     'C '   ,&
     'N '   ,&
     'O '   ,&
     'F '   ,&
     'Ne'  ,&
     'Na'  ,&
     'Mg'  ,&
     'Al'  ,&
     'Si'  ,&
     'P '   ,&
     'S '   ,&
     'Cl'  ,&
     'Ar'  ,&
     'K '   ,&
     'Ca'  ,&
     'Sc'  ,&
     'Ti'  ,&
     'V '   ,&
     'Cr'  ,&
     'Mn'  ,&
     'Fe'  ,& 
     'Co'  ,&
     'Ni'  ,&
     'Cu'  ,&
     'Zn'  ,&
     'Ga'  ,&
     'Ge'  ,&
     'As'  ,&
     'Se'  ,&
     'Br'  ,&
     'Kr'  ,&
     'Rb'  ,&
     'Sr'  ,&
     'Y '   ,&
     'Zr'  ,&
     'Nb'  ,&
     'Mo'  ,&
     'Tc'  ,&
     'Ru'  ,&
     'Rh'  ,&
     'Pd'  ,&
     'Ag'  ,&
     'Cd'  ,&
     'In'  ,&
     'Sn'  ,&
     'Sb'  ,&
     'Te'  ,&
     'I '   ,&
     'Xe'  ,&
     'Cs'  ,&
     'Ba'  ,&
     'La'  ,&
     'Ce'  ,&
     'Pr'  ,&
     'Nd'  ,&
     'Pm'  ,&
     'Sm'  ,&
     'Eu'  ,&
     'Gd'  ,&
     'Tb'  ,&
     'Dy'  ,&
     'Ho'  ,&
     'Er'  ,&
     'Tm'  ,&
     'Yb'  ,&
     'Lu'  ,&
     'Hf'  ,&
     'Ta'  ,&
     'W '   ,&
     'Re'  ,&
     'Os'  ,&  
     'Ir'  ,&
     'Pt'  ,&
     'Au'  ,&
     'Hg'  ,&
     'Tl'  ,&
     'Pb'  ,&
     'Bi'  ,&
     'Po'  ,&
     'At'  ,&
     'Rn'  ,&
     'Fr'  ,&
     'Ra'  ,&
     'Ac'  ,&
     'Th'  ,&
     'Pa'  ,&
     'U '   ,&
     'Np'  ,&
     'Pu'  ,&
     'Am'  ,&
     'Cm'  ,&
     'Bk'  ,&
     'Cf'  ,&
     'Es'  ,&
     'Fm'  ,&
     'Md'  ,&
     'No'  ,&
     'Lr'  ,&
     'Rf'  ,&
     'Db'  ,&
     'Sg'  ,&
     'Bh'  ,&
     'Hs'  ,&
     'Mt'  ,&
     'Ds'  ,&
     'Rg'  ,&
     'Cn'  /)

  integer*4,   dimension(112)            :: n_qnumb = &
   (/1  ,&
     1  ,&
     2  ,&
     2  ,&
     2  ,&
     2  ,&
     2  ,&
     2  ,&
     2  ,&
     2  ,&
     3  ,&
     3  ,&
     3  ,&
     3  ,&
     3  ,&
     3  ,&
     3  ,&
     3  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     4  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     5  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     6  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  ,&
     7  /)


  ! l_qnumb( l ,1:112), s=1, p=2, d=3, f=4
  ! if not zero, the number represent orbital (valence) configuration, for example,
  ! for Silicon with 3s3p for their valence configuration,
  ! l_qnumb( 1 , index_Si) = 3
  ! l_qnumb( 2 , index_Si) = 3
  ! l_qnumb( 3 , index_Si) = 0
  ! l_qnumb( 4 , index_Si) = 0
  integer*4,   dimension(4, 112) :: l_qnumb = reshape( &
           ! s p d f
       (/    1,0,0,0       ,&
             1,0,0,0       ,&
             2,0,0,0       ,&
             2,0,0,0       ,&
             2,2,0,0       ,&
             2,2,0,0       ,&
             2,2,0,0       ,&
             2,2,0,0       ,&
             2,2,0,0       ,&
             2,2,0,0       ,&
             3,0,0,0       ,&
             3,0,0,0       ,&
             3,3,0,0       ,&
             3,3,0,0       ,&
             3,3,0,0       ,&
             3,3,0,0       ,&
             3,3,0,0       ,&
             3,3,0,0       ,&
             4,0,0,0       ,&
             4,0,0,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,0,3,0       ,&
             4,4,3,0       ,&
             4,4,3,0       ,&
             4,4,3,0       ,&
             4,4,3,0       ,&
             4,4,3,0       ,&
             4,4,3,0       ,&
             5,0,0,0       ,&
             5,0,0,0       ,&
             5,0,4,0       ,&
             5,0,4,0       ,&
             5,0,4,0       ,&
             5,0,4,0       ,&
             5,0,4,0       ,&
             5,0,4,0       ,&
             5,0,4,0       ,&
             0,0,4,0       ,&
             5,0,4,0       ,&
             5,0,4,0       ,&
             5,5,4,0       ,&
             5,5,4,0       ,&
             5,5,4,0       ,&
             5,5,4,0       ,&
             5,5,4,0       ,&
             5,5,4,0       ,&
             6,0,0,0       ,&
             6,0,0,0       ,&
             6,0,5,0       ,&
             6,0,5,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,5,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,0,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             6,0,5,4       ,&
             0,6,0,0       ,&
             0,6,0,0       ,&
             5,6,0,0       ,& 
             0,6,0,0       ,&
             0,6,0,0       ,&
             0,6,0,0       ,&
             7,0,0,0       ,&
             7,0,0,0       ,&
             7,0,6,0       ,&
             7,0,6,0       ,&
             7,0,6,5       ,&
             7,0,6,5       ,&
             7,0,6,5       ,&
             7,0,0,5       ,&
             7,0,0,5       ,&
             7,0,6,5       ,&
             7,0,0,5       ,&
             7,0,0,5       ,&
             7,0,0,5       ,&
             7,0,0,5       ,&
             7,0,0,5       ,&
             7,0,0,5       ,&
             7,7,0,5       ,&
             7,0,5,5       ,&
             0,0,0,0       ,&
             0,0,0,0       ,&
             0,0,0,0       ,&
             0,0,0,0       ,&
             0,0,0,0       ,&
             0,0,0,0       ,&
             0,0,0,0       ,&
             0,0,0,0       /), (/4,112/) )

 character*1,parameter ::  angular(4)=(/'s','p','d','f'/)

contains

 function  z_eff(ispec, orb_n, l) result(Zeff)
    implicit none
    integer*4    ispec
    integer*4    orb_n
    integer*4    l
    character*2  orb
    real*8       Zeff

    write(orb,'(I0,A1)')orb_n, angular(l)

    select case (ispec)

      case(1) !'H'   
        if(orb .eq. '1s') Zeff = 1d0

      case(2) !'He'
        if(orb .eq. '1s') Zeff = 1.6875d0
        
      ! http://www.knowledgedoor.com/2/elements_handbook/clementi-raimondi_effective_nuclear_charge_part_4.html#sulfur
      case(16) !'S'
        if(orb .eq. '1s') Zeff = 15.5409d0
        if(orb .eq. '2s') Zeff = 10.629d0  
        if(orb .eq. '2p') Zeff = 11.977d0
        if(orb .eq. '3s') Zeff =  6.3669d0
        if(orb .eq. '3p') Zeff =  5.4819d0

      case(83) !'Bi'
        if(orb .eq. '1s') Zeff = 81.3982d0
        if(orb .eq. '2s') Zeff = 61.1760d0
        if(orb .eq. '2p') Zeff = 78.4670d0
        if(orb .eq. '3s') Zeff = 58.8855d0
        if(orb .eq. '3p') Zeff = 59.9322d0
        if(orb .eq. '4s') Zeff = 47.7072d0
        if(orb .eq. '3d') Zeff = 69.5415d0
        if(orb .eq. '4p') Zeff = 46.8504d0
        if(orb .eq. '5s') Zeff = 31.0290d0
        if(orb .eq. '4d') Zeff = 45.2392d0
        if(orb .eq. '5p') Zeff = 29.0210d0
        if(orb .eq. '6s') Zeff = 15.2400d0
        if(orb .eq. '4f') Zeff = 45.0692d0
        if(orb .eq. '5d') Zeff = 24.2440d0
        if(orb .eq. '6p') Zeff = 13.3400d0

      case(42) !'Mo'
        if(orb .eq. '1s') Zeff = 41.1256d0
        if(orb .eq. '2s') Zeff = 30.8768d0
        if(orb .eq. '2p') Zeff = 37.9718d0
        if(orb .eq. '3s') Zeff = 25.9820d0
        if(orb .eq. '3p') Zeff = 25.4740d0
        if(orb .eq. '4s') Zeff = 16.0960d0
        if(orb .eq. '3d') Zeff = 27.2280d0
        if(orb .eq. '4p') Zeff = 14.9770d0
        if(orb .eq. '5s') Zeff =  6.1060d0
        if(orb .eq. '4d') Zeff = 11.3920d0


    endselect

    return
 endfunction

 function atomic_number(spec) result (ispec)
    implicit none
    character*8    spec
    character*2    ele
    integer*4      ispec
    integer*4      i, le, ls

 ll:do i = 1, 112
      ele = adjustl(trim(element(i)))
      le  = len_trim(ele)
      ls  = len_trim(spec)

      if( spec(1:ls) .eq. ele(1:le) ) then
        ispec = i
        return
      endif 

    enddo ll
  
    return
 endfunction
endmodule
