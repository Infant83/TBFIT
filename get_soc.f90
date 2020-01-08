subroutine set_ham_soc(H, k, PINPT, neig, NN_TABLE, FIJ, flag_phase)
    use parameters, only: zi, pi, pi2, pauli_x, pauli_y, pauli_z, hopping, incar
    use kronecker_prod, only: kproduct
    use phase_factor
    use print_matrix
    use get_parameter
    implicit none
    interface
      function FIJ(kp,R)
        complex*16 :: FIJ
        real*8, intent(in) :: kp(3)
        real*8, intent(in) :: R(3)
      endfunction
    end interface
    type(hopping)           :: NN_TABLE
    type(incar  )           :: PINPT
    integer*4, intent(in)   :: neig
    integer*4                  nn, i, j
    integer*4                  soc_index, rashba_index
    real*8                     lambda_soc, lambda_rashba
    real*8                     hop_signx, hop_signy, hop_signatom
    real*8                     k(3)
    complex*16                 H(neig*2,neig*2) 
    complex*16                 Hx(neig,neig), Hy(neig,neig), Hz(neig,neig)
    complex*16                 F
    character*8                ci_orb, cj_orb
    character*20               ci_atom, cj_atom
    complex*16                 L_x, L_y, L_z
    external                   L_x, L_y, L_z
    logical                    flag_phase
    complex*16                 prod
    real*8                     lsign, hsign

    if(PINPT%flag_slater_koster) then
      Hx = 0d0
      Hy = 0d0
      Hz = 0d0
      H  = 0d0

      do nn = 1, NN_TABLE%n_neighbor
        soc_index = NN_TABLE%soc_param_index(nn)
        i = NN_TABLE%i_matrix(nn) ; j = NN_TABLE%j_matrix(nn)

        ! set SOC hamiltonian based on atomic orbitals
        if( soc_index .gt. 0 .and. (NN_TABLE%p_class(nn) .eq. 'pp' .or. NN_TABLE%p_class(nn) .eq. 'dd') ) then
          call get_param(PINPT, soc_index, 1, lambda_soc)
         !if(PINPT%slater_koster_type .gt. 10) then
         !  call get_param_nrl(PINPT,    soc_index, lambda_soc   )
         !else
         !  call get_param(PINPT,    soc_index, lambda_soc   )
         !endif
          ! CALCULATE  <orb_i|LS|orb_j> 
          Hx(i,j) = lambda_soc * L_x(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))
          Hy(i,j) = lambda_soc * L_y(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))
          Hz(i,j) = lambda_soc * L_z(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))
       
          Hx(j,i) = conjg(Hx(i,j))
          Hy(j,i) = conjg(Hy(i,j))
          Hz(j,i) = conjg(Hz(i,j))

        ! set SOC hamiltonian based on 'xx' type orbitals which are composed by linear combination of atomic orbitals
        elseif( soc_index .gt. 0 .and. NN_TABLE%p_class(nn) .eq. 'xx' ) then
          call get_param(PINPT, soc_index, 1, lambda_soc)
         !if(PINPT%slater_koster_type .gt. 10) then
         !  call get_param_nrl(PINPT,    soc_index, lambda_soc   )
         !else  
         !  call get_param(PINPT,    soc_index, lambda_soc   )
         !endif

          ! CALCULATE  <orb_i|LS|orb_j> 
          Hx(i,j) = lambda_soc * L_x(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))
          Hy(i,j) = lambda_soc * L_y(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))
          Hz(i,j) = lambda_soc * L_z(NN_TABLE%ci_orb(nn), NN_TABLE%cj_orb(nn), NN_TABLE%p_class(nn))

          Hx(j,i) = conjg(Hx(i,j)) 
          Hy(j,i) = conjg(Hy(i,j)) 
          Hz(j,i) = conjg(Hz(i,j)) 

        endif

      enddo

      !SET UP Hamiltonian H_soc*sigma   
      H = kproduct(pauli_x, Hx, 2, 2, neig, neig) &
         +kproduct(pauli_y, Hy, 2, 2, neig, neig) &
         +kproduct(pauli_z, Hz, 2, 2, neig, neig)

    elseif(.not.PINPT%flag_slater_koster) then
      Hx = 0d0
      Hy = 0d0
      Hz = 0d0
      H  = 0d0
      
     !WARNING!! This setting is only valid for Bi/Si(110) case with certain atomic geometry and lattice vectors,
     !          since the sign convention is only valid and meaningful for this particular case.
     !          If you are dealing with other system, please construct your own hamltonian setup.
nn_cc:do nn = 1, NN_TABLE%n_neighbor
        soc_index    = NN_TABLE%cc_index_set(2,nn)
        rashba_index = NN_TABLE%cc_index_set(3,nn)
        i = NN_TABLE%i_matrix(nn) ; j = NN_TABLE%j_matrix(nn)

        if(flag_phase) then
          F = FIJ(k, NN_TABLE%Rij(:,nn))
        elseif(.not. flag_phase) then
          F = FIJ(k, NN_TABLE%R  (:,nn))
        endif

        if( soc_index .ge. 1 .and. rashba_index .ge. 1) then
          call get_param(PINPT,    soc_index, 1, lambda_soc   )
          call get_param(PINPT, rashba_index, 1, lambda_rashba)
         !call get_param(PINPT,    soc_index, lambda_soc   )
         !call get_param(PINPT, rashba_index, lambda_rashba)

          ! set Rashba-SOC between i_orb and j_orb separated by |dij|, originated from E-field normal to surface
          H(i,j+neig) = H(i,j+neig) + zi*lambda_rashba * NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn) * F &               ! sigma_x
                                    + zi*lambda_rashba *-NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn) * F * -zi           ! sigma_y
          H(j,i+neig) = H(j,i+neig) + zi*lambda_rashba *-NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn) * conjg(F) &        ! sigma_x
                                    + zi*lambda_rashba * NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn) * conjg(F) * -zi    ! sigma_y 
          H(j+neig,i) = conjg(H(i,j+neig))
          H(i+neig,j) = conjg(H(j,i+neig))


          ! set SOC between i_orb and j_orb separated by |dij|, 
          ! originated from E-field due to neighbor atom nearby the hopping path
          ! H_SOC = sigma_<<ij>> i * lsoc * v_ij * ci' * sigma_z * cj
          ! v_ij = di x dj / (|di x dj|) , di(dj) are vector connecting nearest neigbohor
          ! atom from i (to j).
          hop_signx= NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn)
          hop_signy= NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn)
          ci_atom = NN_TABLE%site_cindex(NN_TABLE%i_atom(nn))
          cj_atom = NN_TABLE%site_cindex(NN_TABLE%j_atom(nn))
          if( ci_atom(1:2) .eq. 'b1') hop_signatom = 1.0
          if( ci_atom(1:2) .eq. 'b2') hop_signatom =-1.0
          H(i,j)           = H(i     ,j     ) + zi*lambda_soc * hop_signatom * hop_signx * F & ! for spin-up sigma_z
                                              + zi*lambda_soc * hop_signatom * hop_signy * F
          H(i+neig,j+neig) = H(i+neig,j+neig) - zi*lambda_soc * hop_signatom * hop_signx * F & ! for spin-dn sigma_z
                                              - zi*lambda_soc * hop_signatom * hop_signy * F
          if(i .ne. j) then
            H(j,i)            = conjg(H(i,j))
            H(j+neig, i+neig) = conjg(H(i+neig, j+neig))
          endif
        elseif( soc_index .ge. 1 .and. rashba_index .eq. 0) then
          call get_param(PINPT,    soc_index, 1, lambda_soc   )
         !call get_param(PINPT,    soc_index, lambda_soc   )

          ! This model is only for Kane-mele type of SOC. Be careful..
          prod=exp(-2d0*zi * pi * dot_product((/2.45d0,0d0/), NN_TABLE%Rij(1:2,nn)))
          hsign  = sign(1d0,aimag(prod))
          ci_atom = NN_TABLE%site_cindex(NN_TABLE%i_atom(nn))
          cj_atom = NN_TABLE%site_cindex(NN_TABLE%j_atom(nn))
          if( ci_atom(1:1) .eq. 'a') lsign = -1.0d0
          if( ci_atom(1:1) .eq. 'b') lsign =  1.0d0
          H(i,j)            = H(i,j)           + zi * lambda_soc * lsign * hsign * F 
          H(i+neig,j+neig)  = H(i+neig,j+neig) - zi * lambda_soc * lsign * hsign * F 
          if(i .ne. j) then
            H(j,i)            = conjg(H(i,j))
            H(j+neig, i+neig) = conjg(H(i+neig, j+neig))
          endif
        elseif( soc_index .eq. 0 .and. rashba_index .gt. 1 ) then ! WARN: only the AB-a hopping is considered
          call get_param(PINPT, rashba_index, 1, lambda_rashba)
         !call get_param(PINPT, rashba_index, lambda_rashba)

          ! set Rashba-SOC between i_orb and j_orb separated by |dij|, originated from E-field normal to surface
          H(i,j+neig) = H(i,j+neig) + zi*lambda_rashba * NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn) * F &               ! sigma_x
                                    + zi*lambda_rashba *-NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn) * F * -zi           ! sigma_y
          H(j,i+neig) = H(j,i+neig) + zi*lambda_rashba *-NN_TABLE%Rij(2,nn)/NN_TABLE%Dij(nn) * conjg(F) &        ! sigma_x
                                    + zi*lambda_rashba * NN_TABLE%Rij(1,nn)/NN_TABLE%Dij(nn) * conjg(F) * -zi    ! sigma_y 
          H(j+neig,i) = conjg(H(i,j+neig))
          H(i+neig,j) = conjg(H(j,i+neig))
        endif

      enddo nn_cc 

    endif
 

return
endsubroutine
subroutine get_soc_param_index(index_lambda,ci_orb, cj_orb, c_atom, PINPT, param_class)
    use parameters, only : incar
    implicit none
    type(incar) :: PINPT
    integer*4     index_lambda
    integer*4     i, lio, ljo, la
    character*8   ci_orb, cj_orb, c_atom
    character*20  lambda_name
    character*2   param_class

    lio = len_trim(ci_orb)
    ljo = len_trim(cj_orb)
    la = len_trim(c_atom)
    index_lambda = 0 !initialize

    write(lambda_name,*)'lambda_',param_class(1:1),'_',c_atom(1:la)

    if(ci_orb(1:lio) .ne. cj_orb(1:ljo)) then 
     call get_param_index(PINPT, lambda_name, index_lambda)
    endif

return
endsubroutine
function L_x(ci_orb, cj_orb, param_class)
   use parameters, only : zi
   implicit none
   integer*4   li,lj
   character*8 ci_orb,cj_orb
   character*2 param_class
   complex*16  L_x

   L_x = 0d0

   select case (param_class)
     case('pp')
       if(ci_orb(1:2) .eq. 'py' .and. cj_orb(1:2) .eq. 'pz') then
         L_x = -zi * 0.5d0
       elseif(ci_orb(1:2) .eq. 'pz' .and. cj_orb(1:2) .eq. 'py') then
         L_x =  zi * 0.5d0
       endif

     case('dd')
       if(ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:3) .eq. 'dyz' ) then
         L_x =  zi * sqrt(3d0) * 0.5d0
       elseif(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dz2' ) then
         L_x = -zi * sqrt(3d0) * 0.5d0
       endif

       if(ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dyz' ) then
         L_x =  zi * 0.5d0
       elseif(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dx2' ) then
         L_x = -zi * 0.5d0
       endif

       if(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dxy' ) then
         L_x =  zi * 0.5d0
       elseif(ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dxz' ) then
         L_x = -zi * 0.5d0
       endif

     case('xx')
       if(ci_orb(1:3) .eq. 'xp1' .and. cj_orb(1:3) .eq. 'xp2' ) then
!        L_x = zi * sqrt(3d0) * 0.5d0 / sqrt(3d0)
         L_x = zi             * 0.5d0 
       elseif(ci_orb(1:3) .eq. 'xp2' .and. cj_orb(1:3) .eq. 'xp1' ) then
!        L_x =-zi * sqrt(3d0) * 0.5d0 / sqrt(3d0)
         L_x =-zi             * 0.5d0
       endif
       

   end select 
return
endfunction
function L_y(ci_orb, cj_orb, param_class)
   use parameters, only : zi
   implicit none
   character*8 ci_orb,cj_orb
   character*2 param_class
   complex*16  L_y

   L_y = 0d0

   select case (param_class)
     case('pp')
       if(ci_orb(1:2) .eq. 'px' .and. cj_orb(1:2) .eq. 'pz') then
         L_y =  zi * 0.5d0
       elseif(ci_orb(1:2) .eq. 'pz' .and. cj_orb(1:2) .eq. 'px') then
         L_y = -zi * 0.5d0
       endif

     case('dd')
       if(ci_orb(1:3) .eq. 'dz2' .and. cj_orb(1:3) .eq. 'dxz' ) then
         L_y = -zi * sqrt(3d0) * 0.5d0
       elseif(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dz2' ) then
         L_y =  zi * sqrt(3d0) * 0.5d0
       endif
   
       if(ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dxz' ) then
         L_y =  zi * 0.5d0
       elseif(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dx2' ) then
         L_y = -zi * 0.5d0
       endif

       if(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dxy' ) then
         L_y = -zi * 0.5d0
       elseif(ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dyz' ) then
         L_y =  zi * 0.5d0
       endif

     case('xx')
       if(ci_orb(1:3) .eq. 'xp1' .and. cj_orb(1:3) .eq. 'xp3' ) then
!        L_y = -zi * sqrt(3d0) * 0.5d0 / sqrt(3d0)
         L_y = -zi             * 0.5d0
       elseif(ci_orb(1:3) .eq. 'xp3' .and. cj_orb(1:3) .eq. 'xp1' ) then
!        L_y =  zi * sqrt(3d0) * 0.5d0 / sqrt(3d0)
         L_y =  zi             * 0.5d0
       endif

   end select
return
endfunction
function L_z(ci_orb, cj_orb, param_class)
   use parameters, only : zi
   implicit none
   character*8 ci_orb,cj_orb
   character*2 param_class
   complex*16  L_z

   L_z = 0d0

   select case (param_class)
     case('pp')
       if(ci_orb(1:2) .eq. 'px' .and. cj_orb(1:2) .eq. 'py') then
         L_z = -zi * 0.5d0
       elseif(ci_orb(1:2) .eq. 'py' .and. cj_orb(1:2) .eq. 'px') then
         L_z =  zi * 0.5d0
       endif

     case('dd')
       if(ci_orb(1:3) .eq. 'dx2' .and. cj_orb(1:3) .eq. 'dxy' ) then
         L_z = -zi 
       elseif(ci_orb(1:3) .eq. 'dxy' .and. cj_orb(1:3) .eq. 'dx2' ) then
         L_z =  zi
       endif

       if(ci_orb(1:3) .eq. 'dyz' .and. cj_orb(1:3) .eq. 'dxz' ) then
         L_z =  zi * 0.5d0
       elseif(ci_orb(1:3) .eq. 'dxz' .and. cj_orb(1:3) .eq. 'dyz' ) then
         L_z = -zi * 0.5d0
       endif

     case('xx')
       if(ci_orb(1:3) .eq. 'xp2' .and. cj_orb(1:3) .eq. 'xp3' ) then
!        L_z =  zi * ( 0.5d0 / 3d0 - 2d0 / 3d0 )
         L_z = -zi * 0.5d0
       elseif(ci_orb(1:3) .eq. 'xp3' .and. cj_orb(1:3) .eq. 'xp2' ) then
!        L_z = -zi * ( 0.5d0 / 3d0 - 2d0 / 3d0 )
         L_z =  zi * 0.5d0
       endif


   end select
return
endfunction
