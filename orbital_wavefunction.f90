module orbital_wavefunction
   implicit none

contains

!function psi_rho(phi_r, nbasis, iee, ikk, ETBA, flag_plot_wavefunction, spin_index)
function psi_rho(phi_r, nbasis, ispin, V,SV, flag_plot_wavefunction, spin_index, flag_use_overlap)
   use parameters, only : incar, energy
   type(energy)  :: ETBA
   integer*4        nbasis, ispin
   complex*16       phi_r(nbasis)
   complex*16       psi_rho
   complex*16       V(nbasis*ispin)
   complex*16       SV(nbasis*ispin)
   character*2      spin_index
   logical          flag_plot_wavefunction, flag_use_overlap

   select case(spin_index)
     case('up')
       if(flag_plot_wavefunction) then
         psi_rho = sum( phi_r(:)*V(1:nbasis) )
       elseif(.not. flag_plot_wavefunction) then
         if(.not. flag_use_overlap) then
           psi_rho = sum( phi_r(:)*V(1:nbasis)*conjg(V(1:nbasis)) )
         else
           psi_rho = sum( phi_r(:)*V(1:nbasis)*conjg(SV(1:nbasis)) )
         endif
       endif

     case('dn') 
       if(flag_plot_wavefunction) then
         psi_rho = sum( phi_r(:)*V(1+nbasis:nbasis*2) )
       elseif(.not. flag_plot_wavefunction) then
         if(.not. flag_use_overlap) then
           psi_rho = sum( phi_r(:)*V(1+nbasis:nbasis*2)*conjg(V(1+nbasis:nbasis*2)) )
         else
           psi_rho = sum( phi_r(:)*V(1+nbasis:nbasis*2)*conjg(SV(1+nbasis:nbasis*2)) )
         endif
       endif

   end select

   return
endfunction

!elemental complex*16 function get_phi_r(xx,yy,zz,c_orb) !,Zeff,l)
elemental complex*16 function get_phi_r(xx,yy,zz,c_orb,Zeff,orb,n,l) !,Zeff,l)
   use parameters, only : bohr, pi, pi2, zi, pzi, pzi2, rt2, rt3
   use element_info, only: angular
   real*8,intent(in) ::   xx, yy, zz
   real*8                 x_, y_, z_
   real*8                 rr, r2 , rho, R
   complex*16             Y
   complex*16             Y1, Y2
   real*8                 Z, Z1, Z2, ZTa  ! effective nuclear charge
   real*8,intent(in) ::   Zeff ! effective nuclear charge
   real*8,intent(in) ::   n   ! principal quantum number 
   integer*4,intent(in) :: l   ! angular momentum
   character*2,intent(in) :: orb ! 1s, 2s, 3p, 5d ???
   character(*), intent(in) :: c_orb

!  for Ta_5d orbitals, Z_ = 16.37
!  (see https://www.webelements.com/periodicity/eff_nuc_chg_clem_6s/ for 6s orbitals)
!  S_3p orbitals, Z_ =  5.48
!  http://www.knowledgedoor.com/2/elements_handbook/clementi-raimondi_effective_nuclear_charge_part_2.html#molybdenum
!  https://winter.group.shef.ac.uk/orbitron/AOs/3p/equations.html
   !Z  = 16.37d0*1.5d0 / bohr  ! for Ta 5d
   !Z  = 11.392d0*1.5d0 / bohr  ! for Mo 4d
   ! Z  = 10.808d0*1.5d0 / bohr  ! for Te 5p
!  Z1  = 11.392d0 ! for Mo 4d
!  Z2  = 10.808d0 ! for Te 5p
!  ZTa = 16.37d0  ! for Ta 5d

   x_  = xx * bohr
   y_  = yy * bohr
   z_  = zz * bohr
   rr = sqrt( x_**2 + y_**2 + z_**2 )
   r2 = rr**2

   Z   = Zeff
   rho = 2d0 * rr * Z / real(l)

!  select case ( c_orb(1:1) )
!    case('d')
!      Z = Z1
!    case('p')
!      Z = Z2
!    case('x')
!      Z = ZTa
!  endselect

   select case ( orb )
     case('1s')
       R = 2d0 * (Z**1.5d0) * Exp(-rho/2d0)
     case('2s')
       R = 1d0/(2d0*sqrt(2d0)) * ( 2d0 - rho ) * (Z**1.5d0) * Exp(-rho/2d0)
     case('2p')
       R = 1d0/(2d0*sqrt(6d0)) * rho * (Z**1.5d0) * Exp(-rho/2d0)
     case('3s')
       R = 1d0/(9d0*sqrt(3d0)) * ( 6d0 - 6d0*rho + rho**2) * (Z**1.5d0) * Exp(-rho/2d0)
     case('3p')
       R = 1d0/(9d0*sqrt(6d0)) * rho * ( 4d0 - rho ) * (Z**1.5d0) * Exp(-rho/2d0)
     case('3d')
       R = 1d0/(9d0*sqrt(30d0)) * (rho**2) * (Z**1.5d0) * Exp(-rho/2d0)
     case('4s')
       R = 1d0/(96d0) * (24d0 - 36d0*rho + 12d0*(rho**2) - rho**3) * (Z**1.5d0) * Exp(-rho/2d0)
     case('4p')
       R = 1d0/(32d0*sqrt(15d0)) * rho* (20d0 - 10d0*rho + rho**2 ) * (Z**1.5d0) * Exp(-rho/2d0)
     case('4d')
       R = 1d0/(96d0*sqrt(5d0)) * (rho**2) * (6d0 - rho) * (Z**1.5d0) * Exp(-rho/2d0)
   !  case('4f')
      ! not available in the current version.
     case('5s')
       R = 1d0/(300d0*sqrt(5d0)) * (120d0 - 240d0*rho + 120d0*(rho**2) - 20d0*(rho**3)+rho**4) * (Z**1.5d0) * Exp(-rho/2d0)
     case('5p')
       R = 1d0/(150d0*sqrt(30d0)) * rho * (120d0 - 90d0*rho + 18d0*(rho**2) - rho**3) * (Z**1.5d0) * Exp(-rho/2d0)
     case('5d')
       R = 1d0/(150d0*sqrt(70d0)) * (rho**2) * (42d0 - 14d0*rho + rho**2) * (Z**1.5d0) * Exp(-rho/2d0)
     case('6s')
       R = 1d0/(2160d0*sqrt(6d0)) * (720d0 - 1800d0*rho + 1200d0*(rho**2) - 300d0*(rho**3) + 30d0*(rho**4) - rho**5) * (Z**1.5d0) * Exp(-rho/2d0)
     case('6p')
       R = 1d0/(432d0*sqrt(210d0)) * rho * (840d0 - 840d0*rho + 252d0*(rho**2) - 28d0*(rho**3) + (rho**4)) * (Z**1.5d0) * Exp(-rho/2d0)
     case('6d')
       R = 1d0/(864d0*sqrt(105d0)) * (rho**2) * (336d0 - 168d0*rho + 24d0*(rho**2) - (rho**3) ) * (Z**1.5d0) * Exp(-rho/2d0)

   endselect

   select case (trim(c_orb))
     case('s') 
       Y = 1d0 * sqrt(1d0 / (4d0*pi))
     case('px')
       Y   = sqrt(3d0) * x_ / rr * sqrt( 1d0/(4d0*pi) )
     case('py')
       Y   = sqrt(3d0) * y_ / rr * sqrt( 1d0/(4d0*pi) )
     case('pz')
       Y   = sqrt(3d0) * z_ / rr * sqrt( 1d0/(4d0*pi) )
     case('dz2')
       Y   = sqrt(5d0/4d0) * (3d0*(z_**2) - r2)/r2 * sqrt( 1d0/(4d0*pi) )
     case('dx2')
       Y   = sqrt(15d0/4d0) * ( x_**2 - y_**2 ) / r2 * sqrt( 1d0/(4d0*pi) )
     case('dxy')
       Y   = sqrt(15d0) * x_ * y_ / r2 * sqrt( 1d0/(4d0*pi) )
     case('dxz')
       Y   = sqrt(15d0) * x_ * z_ / r2 * sqrt( 1d0/(4d0*pi) )
     case('dyz')
       Y   = sqrt(15d0) * y_ * z_ / r2 * sqrt( 1d0/(4d0*pi) )

     case('xp1') ! phi1 = dz2
      Y   = sqrt(5d0/4d0) * (3d0*(z_**2) - r2)/r2 * sqrt( 1d0/(4d0*pi) )

     case('xp2') ! phi2 = -2/r6*dx2 + 1/r3*dyz
      Y1  = sqrt(15d0/4d0) * ( x_**2 - y_**2 ) / r2 * sqrt( 1d0/(4d0*pi) )  ! dx2
      Y2  = sqrt(15d0) * y_ * z_ / r2 * sqrt( 1d0/(4d0*pi) )  ! dyz
      Y   = -2d0/rt2/rt3 *  Y1 + 1d0/rt3 * Y2

     case('xp3') ! phi2 = 2i/r6*dxy - i/r3*dxz ==> phi2_ = phi2 * zi
      Y1  = sqrt(15d0) * x_ * y_ / r2 * sqrt( 1d0/(4d0*pi) )  ! dxy
      Y2  = sqrt(15d0) * x_ * z_ / r2 * sqrt( 1d0/(4d0*pi) )  ! dzx
      Y   = -2d0/rt2/rt3 * Y1 + 1d0/rt3 * Y2

   endselect


   get_phi_r = R * Y

return
endfunction

end module

