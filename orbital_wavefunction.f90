module orbital_wavefunction
   implicit none

contains

function psi_rho(phi_r, nbasis, iee, ikk, ETBA, flag_plot_wavefunction, spin_index)
   use parameters, only : incar, energy
   type(energy)  :: ETBA
   integer*4        nbasis, iee, ikk
   complex*16       phi_r(nbasis)
   complex*16       psi_rho
   character*2      spin_index
   logical          flag_plot_wavefunction

   select case(spin_index)
     case('up')
       if(flag_plot_wavefunction) then
         psi_rho = sum( phi_r(:)*ETBA%V(1:nbasis,iee,ikk) )
       elseif(.not. flag_plot_wavefunction) then
         psi_rho = sum( phi_r(:)*ETBA%V(1:nbasis,iee,ikk)*conjg(ETBA%V(1:nbasis,iee,ikk)) )
       endif

     case('dn') 
       if(flag_plot_wavefunction) then
         psi_rho = sum( phi_r(:)*ETBA%V(1+nbasis:nbasis*2,iee,ikk) )
       elseif(.not. flag_plot_wavefunction) then
         psi_rho = sum( phi_r(:)*ETBA%V(1+nbasis:nbasis*2,iee,ikk)*conjg(ETBA%V(1+nbasis:nbasis*2,iee,ikk)) )
       endif

   end select

   return
endfunction

elemental complex*16 function get_phi_r(xx,yy,zz,c_orb)
   use parameters, only : bohr, pi, pi2, zi, pzi, pzi2, rt2, rt3
   real*8,intent(in) :: xx, yy, zz
   real*8               rr, r2 , rho, R
   complex*16           Y
   complex*16           Y1, Y2
   real*8               Z  ! effective nuclear charge
!  real*8               Z_
   character(*), intent(in) :: c_orb

!  for Ta_5d orbitals, Z_ = 16.37
!  (see https://www.webelements.com/periodicity/eff_nuc_chg_clem_6s/ for 6s orbitals)
!  S_3p orbitals, Z_ =  5.48

   Z  = 16.37d0*1.5d0 / bohr 
!  Z  = Z_ / bohr 
   rr = sqrt( xx**2 + yy**2 + zz**2 )
   r2 = rr**2
 

   select case ( trim(c_orb) )
     case('dz2')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y   = sqrt(5d0/4d0) * (3d0*(zz**2) - r2)/r2 * sqrt( 1d0/(4d0*pi) )
     case('dx2')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y   = sqrt(15d0/4d0) * ( xx**2 - yy**2 ) / r2 * sqrt( 1d0/(4d0*pi) )
     case('dxy')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y   = sqrt(15d0) * xx * yy / r2 * sqrt( 1d0/(4d0*pi) )
     case('dxz')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y   = sqrt(15d0) * xx * zz / r2 * sqrt( 1d0/(4d0*pi) )
     case('dyz')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y   = sqrt(15d0) * yy * zz / r2 * sqrt( 1d0/(4d0*pi) )

     case('px')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0 / 9d0 / sqrt(6d0) * Exp(-rho/2d0) * rho * (-rho + 4d0) * (Z**1.5d0)
      Y   = sqrt(3d0) * xx / rr * sqrt( 1d0/(4d0*pi) )
     case('py')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0 / 9d0 / sqrt(6d0) * Exp(-rho/2d0) * rho * (-rho + 4d0) * (Z**1.5d0)
      Y   = sqrt(3d0) * yy / rr * sqrt( 1d0/(4d0*pi) )
     case('pz')
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0 / 9d0 / sqrt(6d0) * Exp(-rho/2d0) * rho * (-rho + 4d0) * (Z**1.5d0)
      Y   = sqrt(3d0) * zz / rr * sqrt( 1d0/(4d0*pi) )

     !user defined atomic orbital
     case('xp1') ! phi1 = dz2
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y   = sqrt(5d0/4d0) * (3d0*(zz**2) - r2)/r2 * sqrt( 1d0/(4d0*pi) )
 
     case('xp2') ! phi2 = -2/r6*dx2 + 1/r3*dyz
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y1  = sqrt(15d0/4d0) * ( xx**2 - yy**2 ) / r2 * sqrt( 1d0/(4d0*pi) )  ! dx2
      Y2  = sqrt(15d0) * yy * zz / r2 * sqrt( 1d0/(4d0*pi) )  ! dyz
      Y   = -2d0/rt2/rt3 *  Y1 + 1d0/rt3 * Y2

     case('xp3') ! phi2 = 2i/r6*dxy - i/r3*dxz ==> phi2_ = phi2 * zi
      rho = 2d0 * rr * Z / 5d0
      R   = 1d0/(150d0*sqrt(70d0)) * Exp(-rho/2d0) * ( 42d0 - rho * 14d0 + rho**2 )* (rho**2) * (Z**1.5d0)
      Y1  = sqrt(15d0) * xx * yy / r2 * sqrt( 1d0/(4d0*pi) )  ! dxy
      Y2  = sqrt(15d0) * xx * zz / r2 * sqrt( 1d0/(4d0*pi) )  ! dzx
      Y   = -2d0/rt2/rt3 * Y1 + 1d0/rt3 * Y2

   end select

   get_phi_r = R * Y

return
endfunction

end module

