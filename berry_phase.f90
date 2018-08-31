#include "alias.inc"
module berry_phase
   use do_math
   use mpi_setup
   use sorting

   interface get_berry_phase
    module procedure :: get_berry_phase_det
    module procedure :: get_berry_phase_svd
   end interface

contains

subroutine get_berry_phase_svd(wcc, kpoint, V, PINPT, PGEOM, nkdiv, erange_tot, nerange_tot)
   use parameters, only: incar, poscar, zi, pi2, pi
   implicit none
   type(incar  ) :: PINPT
   type(poscar)  :: PGEOM
   integer*4        icol, fcol
   integer*4        ik, is, i, nkdiv, neig
   integer*4        nerange_tot, nerange
   integer*4        erange_tot(nerange_tot)
   integer*4        erange(nerange_tot/PINPT%nspin)
   real*8           kpoint(3,nkdiv), dk(3)
   complex*16       V(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin,nkdiv)
   complex*16       M(nerange_tot/PINPT%nspin,nerange_tot/PINPT%nspin)
   complex*16       L(nerange_tot/PINPT%nspin,nerange_tot/PINPT%nspin)
   complex*16       L_eig(nerange_tot/PINPT%nspin,PINPT%nspin)
   complex*16       U (nerange_tot/PINPT%nspin,nerange_tot/PINPT%nspin)
   complex*16       A (nerange_tot/PINPT%nspin,nerange_tot/PINPT%nspin)
   complex*16       WT(nerange_tot/PINPT%nspin,nerange_tot/PINPT%nspin)
   complex*16       phase_shift(PGEOM%neig*PINPT%ispinor)
   real*8           wcc(nerange_tot/PINPT%nspin,PINPT%nspin)

   nerange = nerange_tot/PINPT%nspin
   neig    = PGEOM%neig*PINPT%ispinor
   dk      = kpoint(:,2)-kpoint(:,1)

sp:do is = 1, PINPT%nspin
     call set_identity_mat_c(nerange, L)
     erange = erange_tot(1+(is-1)*nerange:nerange*is)

  kp:do ik = 1, nkdiv-1
       icol=1+neig*(is-1) ; fcol=neig*is 
       call get_phase_shift(phase_shift, dk, PGEOM, PINPT%ispinor)
       call get_overlap_matrix(M,V(icol:fcol,erange,ik:ik+1),phase_shift,neig,nerange)
       call get_svd(M,U,WT,nerange) ! perform singular valued decomposition for unitary rotation of M
       L = matprod(nerange,L,matprod(nerange,U,WT)) 
     enddo kp

     call cal_eig_nonsymm(L, nerange, L_eig(:,is))
     wcc(:,is) = dmod(-arg(L_eig(:,is))/pi2+10d0,1d0) ! get unimodular of -arg(L_eig) where L_eig = exp(-i*phi)
                                                      ! For details See: [D. Gresch et al., PRB 95 075146 (2017)]
     call get_sort_variable_1D(wcc(:,is),nerange,'ascending')
   enddo sp

   return
endsubroutine

subroutine get_berry_phase_det(berryp, kpoint, V, PINPT, PGEOM, nkdiv, erange_tot, nerange_tot)
   use parameters, only: incar, poscar, zi, pi2, pi
   implicit none
   type(incar  ) :: PINPT
   type(poscar)  :: PGEOM
   integer*4        icol, fcol
   integer*4        ik, is,  nkdiv, neig
   integer*4        nerange_tot, nerange
   integer*4        erange_tot(nerange_tot)
   integer*4        erange(nerange_tot/PINPT%nspin)
   real*8           kpoint(3,nkdiv), dk(3)
   complex*16       V(PGEOM%neig*PINPT%ispin,PGEOM%neig*PINPT%ispin,nkdiv)
   complex*16       M(nerange_tot/PINPT%nspin,nerange_tot/PINPT%nspin)
   complex*16       phase_shift(PGEOM%neig*PINPT%ispinor)
   real*8           berryp(PINPT%nspin)
   complex*16       det_M(nkdiv-1,PINPT%nspin)

   nerange = nerange_tot/PINPT%nspin
   neig    = PGEOM%neig*PINPT%ispinor
   dk      = kpoint(:,2)-kpoint(:,1)

sp:do is = 1, PINPT%nspin
     erange = erange_tot(1+(is-1)*nerange:nerange*is)
  kp:do ik = 1, nkdiv-1
       icol=1+neig*(is-1) ; fcol=neig*is
       call get_phase_shift(phase_shift, dk, PGEOM, PINPT%ispinor)
       call get_overlap_matrix(M,V(icol:fcol,erange,ik:ik+1),phase_shift,neig,nerange)
       det_M(ik,is) = determinant(nerange,M)
     enddo kp
     berryp(is) = -arg(product(det_M(:,is)))/pi2

   enddo sp

   return
endsubroutine

subroutine get_phase_shift(phase_shift, dk, PGEOM, ispinor)
   use parameters, only: poscar, incar
   use phase_factor, only : F_IJ
   implicit none
   type(poscar) :: PGEOM
   integer*4       i, is
   integer*4       ibasis
   integer*4       ispinor
   real*8          dk(3)
   complex*16      phase_shift(PGEOM%neig*ispinor)

   do is = 1, ispinor
     do ibasis = 1, PGEOM%neig
       phase_shift(ibasis + (is-1)*PGEOM%neig) = F_IJ(-dk, PGEOM%o_coord_cart(:, ibasis))
     enddo
   enddo

   return
endsubroutine

subroutine get_overlap_matrix(M, V, phase_shift, neig, nerange)
   use parameters, only : eta
   use phase_factor, only : F_IJ
   implicit none
   integer*4       i, j, k
   integer*4       neig, nerange
   complex*16      V(neig, nerange, 2)
   complex*16      M(nerange,nerange)
   complex*16      phase_shift(neig)
   M = (0d0,0d0)
   do j = 1, nerange
     do i = 1, nerange
        M(i,j) = dot_product(V(:,i,1), V(:,j,2) )
!       M(i,j) = dot_product(V(:,i,1), V(:,j,2) *phase_shift)
     enddo
   enddo
   return
endsubroutine

subroutine set_periodic_gauge(V, G, PINPT, PGEOM, nkdiv, erange, nerange)
   use parameters, only : poscar, incar, zi
   use phase_factor
   implicit none
   type(poscar)            :: PGEOM
   type(incar )            :: PINPT
   integer*4                  ie, is, im
   integer*4                  nkdiv, nerange
   integer*4                  erange(nerange)
   real*8    , intent(in)  :: G(3)
   complex*16                 V(PGEOM%neig*PINPT%ispin, PGEOM%neig*PINPT%ispin, nkdiv)

   do ie = 1, nerange
     do is = 1, PINPT%ispin
       do im = 1, PGEOM%neig
         V(im+(is-1)*PGEOM%neig,erange(ie),nkdiv) = V(im+(is-1)*PGEOM%neig,erange(ie),1) 
!        V(im+(is-1)*PGEOM%neig,erange(ie),nkdiv) = V(im+(is-1)*PGEOM%neig,erange(ie),1)  * F_IJ(-G, PGEOM%o_coord_cart(:,im))
       enddo
     enddo
   enddo

   return
endsubroutine

subroutine get_velocity_matrix(PINPT, NN_TABLE, kpoint, neig, dHk, dF_IJ, flag_phase)
   use parameters, only: incar, hopping, pauli_0
   use kronecker_prod, only: kproduct
   use phase_factor
   implicit none
   interface
     function dF_IJ(k,R)
       complex*16 :: dF_IJ
       real*8, intent(in) :: k(3)
       real*8, intent(in) :: R(3)
     endfunction
   end interface
   type(incar  ) :: PINPT
   type(hopping) :: NN_TABLE
   integer*4        neig
   integer*4        i
   real*8           kpoint(3)
   complex*16       dH0(neig,neig)
   complex*16       dHs(neig*PINPT%ispinor,neig*PINPT%ispinor) 
   complex*16       dHk(neig*PINPT%ispinor,neig*PINPT%ispinor) 
   logical          flag_phase

   call set_ham0_vel(dH0, kpoint, PINPT, neig, NN_TABLE, dF_IJ, flag_phase)

   if(PINPT%flag_noncollinear) then

     if(PINPT%flag_slater_koster) then 
       dHk = kproduct(pauli_0, dH0, 2, 2, neig, neig)
     else 
       call set_ham_soc(dHs,kpoint(:),PINPT,neig,NN_TABLE,dF_IJ,flag_phase)
       dHk = kproduct(pauli_0, dH0, 2, 2, neig, neig) !+ dHs
     endif

   else

     dHk = dH0

   endif  

   return
endsubroutine  

subroutine get_velocity_matrix_direct(PINPT, NN_TABLE, kpoint, neig, dHk, dk)
   use parameters, only: incar, hopping, pauli_0, eta
   use kronecker_prod, only: kproduct
   use phase_factor
   implicit none
   type(incar  ) :: PINPT
   type(hopping) :: NN_TABLE
   integer*4        neig
   complex*16       dH0(neig,neig),dH1(neig,neig),dH2(neig,neig)
   complex*16       dHk(neig*PINPT%ispinor,neig*PINPT%ispinor)
   complex*16       Hs1(neig*PINPT%ispinor,neig*PINPT%ispinor)
   complex*16       Hs2(neig*PINPT%ispinor,neig*PINPT%ispinor)
   complex*16       dHs(neig*PINPT%ispinor,neig*PINPT%ispinor)
   real*8           k1(3), k2(3), kpoint(3), dk(3), enorm
   external         enorm
   logical          flag_phase

   flag_phase = .true.
   k1 = kpoint - dk; k2 = kpoint + dk
   call set_ham0_(dH2, k2 , PINPT, neig, NN_TABLE, F_IJ, flag_phase)
   call set_ham0_(dH1, k1 , PINPT, neig, NN_TABLE, F_IJ, flag_phase)
   dH0= (dH2 - dH1) / (2d0*eta)

   if(PINPT%ispinor .eq. 2) then
     if(PINPT%flag_soc .and. .not. PINPT%flag_slater_koster) then 
       call set_ham_soc(Hs2, k2, PINPT, neig, NN_TABLE, F_IJ, flag_phase)
       call set_ham_soc(Hs1, k1, PINPT, neig, NN_TABLE, F_IJ, flag_phase)
       dHs= (Hs2 - Hs1) / (2d0*eta)
       dHk= kproduct(pauli_0, dH0, 2, 2, neig, neig) + dHs
     else
       dHk= kproduct(pauli_0, dH0, 2, 2, neig, neig)
     endif
   else
     dHk= dH0
   endif

   return
endsubroutine

subroutine get_omega(omega_, E_, V_, dxH, dyH,dzH, msize)
   use parameters, only : eta, zi
   implicit none
   integer*4    n, m
   integer*4    msize
   real*8       E_(msize)
   complex*16   V_(msize,msize)
   complex*16   dxH(msize,msize)
   complex*16   dyH(msize,msize)
   complex*16   dzH(msize,msize)
   complex*16   omega(msize,3) 
   real*8       omega_(msize,3) 
   complex*16   vx_nm,vy_nm,vz_nm
   complex*16   vx_mn,vy_mn,vz_mn
   complex*16   psi_n(msize), psi_m(msize)
   real*8       de2
   complex*16   zeta

!  zeta  = zi * eta
   omega = (0d0,0d0)
   omega_= (0d0,0d0)
   do n = 1, msize
     psi_n = V_(:,n)
     do m = 1, msize
       if(m .ne. n) then
         psi_m = V_(:,m)
         de2= (E_(m) - E_(n) )**2 + eta
         vy_nm = dot_product( psi_n, matmul(dyH,psi_m) ) ! < psi_n| dH/dky | psi_m>
         vy_mn = dot_product( psi_m, matmul(dyH,psi_n) ) ! < psi_m| dH/dky | psi_n>
         vz_nm = dot_product( psi_n, matmul(dzH,psi_m) ) 
         vz_mn = dot_product( psi_m, matmul(dzH,psi_n) ) 
         vx_nm = dot_product( psi_n, matmul(dxH,psi_m) ) 
         vx_mn = dot_product( psi_m, matmul(dxH,psi_n) ) 

         omega(n,1) = omega(n,1) + vy_nm * vz_mn / de2
         omega(n,2) = omega(n,2) + vz_nm * vx_mn / de2
         omega(n,3) = omega(n,3) + vx_nm * vy_mn / de2
       endif

     enddo
   enddo

   omega_= -2d0*aimag(omega)
!  omega_= real(omega)
   return
endsubroutine

subroutine set_ham0_sym(H0, neig)
   implicit none
   integer*4    i
   integer*4    neig
   complex*16   diag(neig)
   complex*16   H0(neig,neig)

   do i = 1, neig
     diag(i) = H0(i,i)
     H0(i,i)= (0d0,0d0)
   enddo
   H0 = H0 + conjg(transpose(H0))
   do i = 1, neig
     H0(i,i)= diag(i)
   enddo

   return
endsubroutine

! This routine is to get 1st Chern number of given bands defined with ERANGE
! by tracking the evolution of the polarization over the brillouin zone, as described
! in the eq(A4) of APPENDIX A in Ref.[PRB 95, 075146 (2017)].
subroutine get_chern_number(chern, polarization, wcc, nspin, nkpath, nerange_tot)
   use parameters, only : pi2
   implicit none
   integer*4    ik, is
   integer*4    imin
   integer*4    nerange_tot, nerange
   integer*4    nspin, nkpath
   real*8       wcc(nerange_tot/nspin, nspin, nkpath)
   real*8       polarization_(nspin, nkpath)
   real*8       polarization (nspin, nkpath)
   real*8       chern(nspin)
   real*8       delta_p(5)
   real*8       p(5)

   polarization_= 0d0
   polarization = 0d0
   chern        = 0d0
   nerange      = nerange_tot/nspin

   do is = 1, nspin
     do ik = 1, nkpath
       polarization_(is, ik) = dmod(sum(wcc(:,is,ik)),1d0)
     enddo
   enddo

   polarization(:,1) = polarization_(:,1)
   do is = 1, nspin
     do ik = 1, nkpath - 1
       delta_p(:) =  polarization_(is, ik+1) - polarization_(is, ik)+ (/-2:2/)
       chern(is) = chern(is) + delta_p(minloc( abs(delta_p),1 ))
       polarization(is,ik+1) = polarization(is,ik) + delta_p(minloc( abs(delta_p),1 ))
     enddo    
   enddo

   if(nspin .eq. 2) then
     if_main write(6,'(A,2I3)')'  Chern number (up, dn) = ',nint(chern(:))
   elseif(nspin .eq. 1) then
     if_main write(6,'(A,I3)')'  Chern number = ',nint(chern(1))
   endif

   return
endsubroutine
subroutine find_largest_gap(largest_gap, clock_direct, z2_index, wcc, nspin, nkpath, nerange_tot)
   use parameters, only : pi2
   implicit none
   integer*4    ie, is
   integer*4    ie_old, ie_now
   integer*4    ikpath
   integer*4    nspin
   integer*4    nkpath
   integer*4    nerange_tot, nerange
   real*8       gap_now, gap_old
   real*8       e_now, e_old
   real*8       largest_gap(nspin,nkpath)
   real*8       wcc(nerange_tot/nspin,nspin,nkpath)
   real*8       wcc_(0:nerange_tot/nspin,nspin,nkpath)
   real*8       kk
   integer*4    n_gap_jump(nspin,nkpath)
   integer*4    z2_index(nspin)
   integer*4    clock_direct(nspin,nkpath)
   real*8       phi1, phi2, phi3, directed_area

   nerange             = nerange_tot/nspin
   wcc_(0,:,:)         = wcc(nerange,:,:) - 1d0
   wcc_(1:nerange,:,:) = wcc(1:nerange,:,:)
   largest_gap(:,:)    = (wcc(2,:,:) + wcc(1,:,:))/2d0
   z2_index            = 1
   e_now               = 0d0
   e_old               = 0d0

   do ikpath = 1, nkpath
     kk= 1d0 / (nkpath-1) * (ikpath - 1)
     do is = 1, nspin
       gap_now = 0d0
       gap_old = 0d0
       do ie = 1, nerange
         gap_now = wcc_(ie,is,ikpath) - wcc_(ie-1,is,ikpath)
         if(gap_now .ge. gap_old) then
           gap_old = gap_now
           largest_gap(is, ikpath) = dmod( (wcc_(ie,is,ikpath)+wcc_(ie-1,is,ikpath))/2d0 + 1d0, 1d0 )
         endif
       enddo

     enddo
   enddo


   ! evaluation of topological index: see Fig.2 and Eq.(17,18) of PRB 83, 235401 (2011) for the details.
   do ikpath = 1, nkpath-1
     kk= 1d0 / (nkpath-1) * (ikpath - 1)
     do is = 1, nspin
       clock_direct(is,ikpath) = 1
       phi1 = pi2 * largest_gap(is,ikpath)
       phi2 = pi2 * largest_gap(is,ikpath+1)
       do ie = 1, nerange
        phi3 = pi2 * wcc(ie,is,ikpath+1)
        directed_area = sin(phi2-phi1) + sin(phi3-phi2) + sin(phi1-phi3)
        clock_direct(is,ikpath)  = int4(sign(1d0,directed_area)) * clock_direct(is,ikpath)
       enddo
       if(kk .le. 0.5d0) z2_index(is) = z2_index(is) * clock_direct(is,ikpath)
     enddo
   enddo
   clock_direct(:,nkpath) = clock_direct(:,1)
   do is = 1, nspin
     if(z2_index(is) .eq. -1) then
       z2_index(is) = 1
     elseif(z2_index(is) .eq.  1) then
       z2_index(is) = 0
     else
       write(6,*)'  !WARN! Z2 index is not well defined. stop anyway..'
       stop
     endif
   enddo

   return
endsubroutine
subroutine set_berry_erange(PINPT_BERRY, PGEOM, PINPT, mode)
   use parameters, only : incar, berry, poscar
   implicit none
   type(incar)   :: PINPT
   type(berry)   :: PINPT_BERRY
   type(poscar)  :: PGEOM
   character*2      mode
   integer*4        erange_(PGEOM%neig * PINPT%ispin)


   if(mode .eq. 'zk') then ! zak_phase

     call get_tbrange(PINPT_BERRY%strip_zak_range, PGEOM, PINPT, erange_, PINPT_BERRY%zak_nerange)
     allocate(PINPT_BERRY%zak_erange(PINPT_BERRY%zak_nerange))
     PINPT_BERRY%zak_erange(1:PINPT_BERRY%zak_nerange) = erange_(1:PINPT_BERRY%zak_nerange)

   elseif(mode .eq. 'wc') then ! wannier charge center

     call get_tbrange(PINPT_BERRY%strip_wcc_range, PGEOM, PINPT, erange_, PINPT_BERRY%wcc_nerange)
     allocate(PINPT_BERRY%wcc_erange(PINPT_BERRY%wcc_nerange))
     PINPT_BERRY%wcc_erange(1:PINPT_BERRY%wcc_nerange) = erange_(1:PINPT_BERRY%wcc_nerange)

   elseif(mode .eq. 'bc') then ! berry curvature

     call get_tbrange(PINPT_BERRY%strip_bc_range, PGEOM, PINPT, erange_, PINPT_BERRY%bc_nerange)
     allocate(PINPT_BERRY%bc_erange(PINPT_BERRY%bc_nerange))
     PINPT_BERRY%bc_erange(1:PINPT_BERRY%bc_nerange) = erange_(1:PINPT_BERRY%bc_nerange)

   elseif(mode .eq. 'z2') then ! z2 index

     call get_tbrange(PINPT_BERRY%strip_z2_range, PGEOM, PINPT, erange_, PINPT_BERRY%z2_nerange)
     allocate(PINPT_BERRY%z2_erange(PINPT_BERRY%z2_nerange))
     PINPT_BERRY%z2_erange(1:PINPT_BERRY%z2_nerange) = erange_(1:PINPT_BERRY%z2_nerange)

   endif

   return
endsubroutine

subroutine set_berry_kpath(PINPT_BERRY, PGEOM, PINPT, mode)
   use parameters, only : incar, berry, poscar
   implicit none
   type(incar)   :: PINPT
   type(berry)   :: PINPT_BERRY
   type(poscar)  :: PGEOM
   integer*4        ikpath
   character*2      mode

   if(mode .eq. 'zk') then

     do ikpath = 1, PINPT_BERRY%zak_nkpath
       call get_kline(PINPT_BERRY%zak_kpoint(:,:,ikpath), &
                      PINPT_BERRY%zak_kpoint_reci(:,:,ikpath), PINPT_BERRY%zak_nkdiv, PGEOM,&
                      PINPT_BERRY%zak_kpath(:,1,ikpath), &
                      PINPT_BERRY%zak_kpath(:,2,ikpath) )
     enddo

   elseif(mode .eq. 'wc') then

     do ikpath = 1, PINPT_BERRY%wcc_nkpath
       call get_kline(PINPT_BERRY%wcc_kpoint(:,:,ikpath), &
                      PINPT_BERRY%wcc_kpoint_reci(:,:,ikpath), PINPT_BERRY%wcc_nkdiv, PGEOM,&
                      PINPT_BERRY%wcc_kpath(:,1,ikpath), &
                      PINPT_BERRY%wcc_kpath(:,2,ikpath) )
     enddo

   endif

   return
endsubroutine
subroutine set_kpath_plane(kpoint, kpoint_reci, kpath, nkdiv, nkpath, iaxis, kplane, PGEOM)
   use parameters, only : poscar, cyclic_axis
   implicit none
   type(poscar)  :: PGEOM
   integer*4        nkdiv, nkpath, iaxis
   integer*4        ikpath
   real*8           kplane
   real*8           k_init(3), k_end(3)
   real*8           shiftk(3)
   real*8           kaxis(3,3)
   real*8           kpath(3,2,nkpath)
   real*8           kpoint(3,nkdiv,nkpath)
   real*8           kpoint_reci(3,nkdiv,nkpath)

   kaxis(:,1)     = (/1d0,0d0,0d0/)
   kaxis(:,2)     = (/0d0,1d0,0d0/)
   kaxis(:,3)     = (/0d0,0d0,1d0/)

   shiftk = kplane * kaxis(:,iaxis)
   do ikpath = 1, nkpath
     k_init = shiftk + 1d0/(nkpath-1)*(ikpath-1) * kaxis(:,cyclic_axis(iaxis,1))
     k_end  = shiftk + 1d0/(nkpath-1)*(ikpath-1) * kaxis(:,cyclic_axis(iaxis,1)) + kaxis(:,cyclic_axis(iaxis,2))
     kpath(:,1,ikpath) = k_init;kpath(:,2,ikpath) = k_end
     call get_kline(kpoint(:,:,ikpath), kpoint_reci(:,:,ikpath), nkdiv, PGEOM, k_init, k_end)
   enddo

   return
endsubroutine

subroutine get_1st_chern_number(PINPT, PINPT_BERRY, PKPTS, PGEOM)
   use parameters, only : kpoints, incar, poscar, berry, pi2, cyclic_axis
   use do_math, only : area
   implicit none
   type(incar)   :: PINPT
   type(berry)   :: PINPT_BERRY
   type(kpoints) :: PKPTS
   type(poscar ) :: PGEOM
   integer*4        ie, is, ip, ik, iaxis, iplane
   integer*4        nerange
   integer*4        erange(PINPT_BERRY%bc_nerange/PINPT%nspin)
!  real*8           chern(PINPT_BERRY%bc_nerange/PINPT%nspin, PINPT%nspin, 2, size(PINPT_BERRY%bc_axis))
   real*8           chern(PGEOM%neig*PINPT%ispinor, 3, PINPT%nspin)
   real*8           dk2(3)

   dk2(1) = area(PGEOM%b_latt(:,2), PGEOM%b_latt(:,3)) / real(PKPTS%ndiv(2)*PKPTS%ndiv(3),8)
   dk2(2) = area(PGEOM%b_latt(:,3), PGEOM%b_latt(:,1)) / real(PKPTS%ndiv(3)*PKPTS%ndiv(1),8)
   dk2(3) = area(PGEOM%b_latt(:,1), PGEOM%b_latt(:,2)) / real(PKPTS%ndiv(1)*PKPTS%ndiv(2),8)
   nerange = PINPT_BERRY%bc_nerange/PINPT%nspin

   chern = 0d0

   do is = 1, PINPT%nspin
   do iaxis = 1, 3
   do ie = 1, PGEOM%neig*PINPT%ispinor
     chern(ie, iaxis, is) = chern(ie, iaxis, is) + sum(PINPT_BERRY%omega(ie, iaxis, is, :)) / pi2 * dk2(3)
   enddo
   enddo
   enddo
   write(6,*)"XXXX", chern(:,3,1)
   write(6,*)"XXXX", sum(chern(PINPT_BERRY%bc_erange,3,1))
!stop
!ax:do iaxis = 1, size(PINPT_BERRY%bc_axis)
!  sp:do is = 1, PINPT%nspin
!       erange = PINPT_BERRY%bc_erange(1+(is-1)*nerange:nerange + (is-1)*nerange)
!   eig:do ie = 1, nerange
!      pl:do ip = 1, PKPTS%ndiv(PINPT_BERRY%bc_axis(iaxis))
!        kp:do ik = 1, PKPTS%nkpoint
!             if(PKPTS%kpoint_reci(PINPT_BERRY%bc_axis(iaxis),ik) .eq. 0d0) then
!               iplane = 1
!               chern_number(ie, is, iplane, iaxis) = chern_number(ie, is, iplane, iaxis) + PINPT_BERRY%omega(ie,PINPT_BERRY%bc_axis(iaxis),is,ik)
!             elseif(PKPTS%kpoint_reci(PINPT_BERRY%bc_axis(iaxis),ik) .eq. 0.5d0) then                                                         
!               iplane = 2                                                                                                                     
!               chern_number(ie, is, iplane, iaxis) = chern_number(ie, is, iplane, iaxis) + PINPT_BERRY%omega(ie,PINPT_BERRY%bc_axis(iaxis),is,ik)
!             endif
!           enddo kp
!         enddo pl
!       enddo eig
!     enddo sp
!   enddo ax
!   write(6,*)"XXX", sum(chern_number(:,1,1,1))/pi2 * dk2(PINPT_BERRY%bc_axis(1))


!write(6,*)"XXXXX", nerange, PINPT_BERRY%bc_axis(1)
!   do ie = 1, nerange
!     write(6,*)"XXXX", sum(PINPT_BERRY%omega(ie,1,3,:))/pi2*dk2(PINPT_BERRY%bc_axis(1))
!   enddo

!  stop

   return
endsubroutine

end module
