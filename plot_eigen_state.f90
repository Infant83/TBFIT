#include "alias.inc"
subroutine plot_eigen_state(PINPT, PGEOM, PKPTS, ETBA, flag_use_overlap)
   use parameters, only : incar, poscar, kpoints, energy
   use orbital_wavefunction, only : psi_rho
   use mpi_setup
   use time
   use memory
   use print_io
   implicit none
   type (incar)   :: PINPT
   type (poscar)  :: PGEOM   
   type (kpoints) :: PKPTS   
   type (energy)  :: ETBA
   character(*), parameter :: func = 'plot_eigen_state'
   integer*4    ii, i, ie, ik, iatom, iorb, iorbital, icoeff
   integer*4    ikk, iee, neig, igrid
   integer*4    ix,iy,iz,i1,i2,i3, ngrid, ng1, ng2, ng3
   integer*4    nwrite, iline, nline, nresi
   integer*4    pid_chg_up, pid_chg_dn
   integer*4    mpierr
   real*8       a1(3),a2(3),a3(3), vol
   character*8  corb(PGEOM%neig)
   real*8       origin(3,PGEOM%neig), origin_reset(3,PGEOM%neig)
   real*8       grid_a1(0:PGEOM%ngrid(1)-1), rx(PGEOM%neig)
   real*8       grid_a2(0:PGEOM%ngrid(2)-1), ry(PGEOM%neig)
   real*8       grid_a3(0:PGEOM%ngrid(3)-1), rz(PGEOM%neig)
   complex*16   phi_r(PGEOM%neig)
   complex*16   psi_r_up(PGEOM%ngrid(1)*PGEOM%ngrid(2)*PGEOM%ngrid(3))
   complex*16   psi_r_dn(PGEOM%ngrid(1)*PGEOM%ngrid(2)*PGEOM%ngrid(3))
   complex*16   psi_r_up_(PGEOM%ngrid(1)*PGEOM%ngrid(2)*PGEOM%ngrid(3))
   complex*16   psi_r_dn_(PGEOM%ngrid(1)*PGEOM%ngrid(2)*PGEOM%ngrid(3))
   complex*16   V(PGEOM%neig*PINPT%ispin,PINPT%nspin)
   complex*16   SV(PGEOM%neig*PINPT%ispin,PINPT%nspin)
   real*8       time1, time2
   logical      flag_exist_up, flag_exist_dn
   logical      flag_use_overlap
   real*8        t0, t1
   character*4   timer
   real*8        zeff(PGEOM%neig)
   real*8        nqnum(PGEOM%neig)
   integer*4     lqnum(PGEOM%neig)
   character*2   orb(PGEOM%neig)


   timer = 'init'
   call time_check(t1,t0,timer)
   call set_variable_plot_eig(PINPT, PGEOM, neig, ngrid, nwrite, nline, nresi, pid_chg_up, pid_chg_dn, vol, &
                              ng1, ng2, ng3, a1, a2, a3, origin, corb, grid_a1, grid_a2, grid_a3, &
                              zeff, nqnum, lqnum, orb )
   call write_info_plot_eig(PINPT)

en:do ie = 1, PINPT%n_eig_print
     write(message,'(A,I5)')            "   # EIGEN STATE No. : ",PINPT%i_eig_print(ie)  ; write_msg
     iee = PINPT%i_eig_print(ie) - PGEOM%init_erange + 1
 kp: do ik = 1, PINPT%n_kpt_print
       write(message,'(A,*(I5))'       )"    #   K-POINT No. : ",PINPT%i_kpt_print(ik)  ; write_msg
       ikk = PINPT%i_kpt_print(ik)    

#ifdef MPI
       if(PINPT%nspin .eq. 1) then
         if_main  V(:,1)=ETBA%V( :,iee,ikk)
         if(flag_use_overlap) then
           if_main SV(:,1)=ETBA%SV(:,iee,ikk)
         endif
       elseif(PINPT%nspin .eq. 2) then
         if_main  V(:,:)=ETBA%V( :,(/iee,iee+PGEOM%nband/),ikk)
         if(flag_use_overlap) then
           if_main SV(:,:)=ETBA%SV(:,(/iee,iee+PGEOM%nband/),ikk)
         endif
       endif
       call MPI_BCAST(V, size(V), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
       if(flag_use_overlap) then
         call MPI_BCAST(SV, size(SV), MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
       endif
#else
       if(PINPT%nspin .eq. 1) then
         V(:,1)=ETBA%V(:,iee,ikk)
         if(flag_use_overlap) then
           SV(:,1)=ETBA%SV(:,iee,ikk)
         endif
       elseif(PINPT%nspin .eq. 2) then
         V(:,:)=ETBA%V(:,(/iee,iee+PGEOM%nband/),ikk)
         if(flag_use_overlap) then
           SV(:,:)=ETBA%SV(:,(/iee,iee+PGEOM%nband/),ikk)
         endif
       endif
#endif

       call print_PARCHG_head(PINPT, PGEOM, ik, ie, pid_chg_up, pid_chg_dn, flag_exist_up, flag_exist_dn)
       if(.not. flag_exist_up .and. .not. flag_exist_dn) then
         write(message,'(A)')         "      !WARN! Eigen state No.", PINPT%i_eig_print(ie), " has not been found"  ; write_msg
         write(message,'(A)')         "             at this K-POINT within [EMIN:EMAX] specified in your EWINDOW tag. "  ; write_msg
         write(message,'(A)')         "             No result will be written for this K-POINT. Skip.."  ; write_msg
         cycle kp
       elseif(.not. flag_exist_up .and. flag_exist_dn) then
         write(message,'(A)')         "      !WARN! Eigen state No.", PINPT%i_eig_print(ie), " for spin-up has not been found"  ; write_msg
         write(message,'(A)')         "             at this K-POINT within [EMIN:EMAX] specified in your EWINDOW tag. "  ; write_msg
         write(message,'(A)')         "             Only spin-dn result will be written for this K-POINT. "  ; write_msg
       elseif(flag_exist_up .and. .not. flag_exist_dn) then
         write(message,'(A)')         "      !WARN! Eigen state No.", PINPT%i_eig_print(ie), " for spin-dn has not been found"  ; write_msg
         write(message,'(A)')         "             at this K-POINT within [EMIN:EMAX] specified in your EWINDOW tag. "  ; write_msg
         write(message,'(A)')         "             Only spin-up result will be written for this K-POINT. "  ; write_msg
       endif
       call initialize_psi_r(psi_r_up,psi_r_dn, ngrid,PINPT%ispin, flag_exist_up, flag_exist_dn)

#ifdef MPI
       call MPI_BARRIER(mpi_comm_earth,mpierr)
#endif
cell_z:do iz = -1,1
cell_y:do iy = -1,1
cell_x:do ix = -1,1
         call reset_orbital_origin(origin_reset, origin, neig, a1, a2, a3, ix, iy, iz)
 !grid_z:do i3=0,ng3-1
  grid_z:do i3=0+myid,ng3-1, nprocs
  grid_y:do i2=0,ng2-1
  grid_x:do i1=0,ng1-1
           igrid = i1+1+i2*ng1+i3*ng1*ng2
           call get_rxyz(rx,ry,rz, grid_a1, grid_a2, grid_a3, origin_reset, neig, PGEOM%ngrid, a1, a2, a3, i1,i2,i3)
           call get_orbital_wavefunction_phi_r(phi_r, rx,ry,rz, corb, neig, PINPT%rcut_orb_plot, PINPT%flag_plot_wavefunction, zeff, nqnum, lqnum, orb)
           call get_psi_r(psi_r_up,psi_r_dn,igrid,ngrid,neig,phi_r,V,SV,PINPT%ispin,PINPT%nspin,PINPT%flag_plot_wavefunction, &
                          flag_exist_up, flag_exist_dn, flag_use_overlap)
         enddo grid_x
         enddo grid_y
         enddo grid_z
       enddo cell_x 
       enddo cell_y
       enddo cell_z

#ifdef MPI
       call MPI_BARRIER(mpi_comm_earth,mpierr)
       if(PINPT%nspin .eq. 1) then
         call MPI_ALLREDUCE(psi_r_up, psi_r_up_, size(psi_r_up_), MPI_COMPLEX16, MPI_SUM, mpi_comm_earth, mpierr)
         psi_r_up = psi_r_up_
         if(PINPT%ispinor .eq. 2) then
           call MPI_ALLREDUCE(psi_r_dn, psi_r_dn_, size(psi_r_dn_), MPI_COMPLEX16, MPI_SUM, mpi_comm_earth, mpierr)
           psi_r_dn = psi_r_dn_ 
         endif
       elseif(PINPT%nspin .eq. 2) then
         call MPI_ALLREDUCE(psi_r_up, psi_r_up_, size(psi_r_up_), MPI_COMPLEX16, MPI_SUM, mpi_comm_earth, mpierr)
         call MPI_ALLREDUCE(psi_r_dn, psi_r_dn_, size(psi_r_dn_), MPI_COMPLEX16, MPI_SUM, mpi_comm_earth, mpierr)
         psi_r_up = psi_r_up_
         psi_r_dn = psi_r_dn_
       endif
       call MPI_BARRIER(mpi_comm_earth,mpierr)
#endif

       if_main call write_rho_main(pid_chg_up, pid_chg_dn, ngrid, nline, nwrite, nresi, psi_r_up, psi_r_dn, &
                                   PINPT%ispin, PINPT%nspin, PINPT%ispinor, PINPT%flag_plot_wavefunction, &
                                   flag_exist_up, flag_exist_dn)

     enddo kp
!    write(message,'(A)',ADVANCE='yes')" "  ; write_msg
   enddo en

   call time_check(t1,t0)
   write(message,*)' '  ; write_msg
   write(message,'(A,F10.4,A)')'   TIME for EIGENSTATE PLOT : ',t1, ' (sec)'  ; write_msg
   if(PINPT%flag_plot_wavefunction) then
     write(message,'(A)')' #- END WRITING: EIGENSTATE WAVEFUNCTION'  ; write_msg
   elseif(.not. PINPT%flag_plot_wavefunction) then
     write(message,'(A)')' #- END WRITING: EIGENSTATE CHARGE DENSITY'  ; write_msg
   endif


   return
endsubroutine
subroutine initialize_psi_r(psi_r_up,psi_r_dn, ngrid, ispin, flag_exist_up, flag_exist_dn)
   implicit none
   integer*4    ispin, ngrid
   complex*16   psi_r_up(ngrid), psi_r_dn(ngrid)
   logical      flag_plot_wavefunction
   logical      flag_exist_up, flag_exist_dn

   if(flag_exist_up) psi_r_up = (0d0,0d0)
   
   if( ispin .eq. 2 .and. flag_exist_dn) psi_r_dn = (0d0,0d0)

   return
endsubroutine
subroutine get_psi_r(psi_r_up,psi_r_dn,igrid,ngrid,nbasis,phi_r,V,SV,ispin,nspin,flag_plot_wavefunction, &
                     flag_exist_up, flag_exist_dn,flag_use_overlap)
   use parameters, only : incar, energy
   use orbital_wavefunction, only: psi_rho
   type(incar) :: PINPT
   type(energy):: ETBA 
   integer*4    ikk, iee, nbasis, ispin
   integer*4    igrid, ngrid
   complex*16   phi_r(nbasis)
   complex*16   psi_r_up(ngrid)
   complex*16   psi_r_dn(ngrid)
   complex*16   V(nbasis*ispin,nspin)
   complex*16   SV(nbasis*ispin,nspin)
   logical      flag_plot_wavefunction
   logical      flag_exist_up, flag_exist_dn
   logical      flag_use_overlap

   if(flag_exist_up) then
     psi_r_up(igrid) = psi_r_up(igrid) + psi_rho(phi_r, nbasis, ispin, V(:,1),SV(:,1), flag_plot_wavefunction, 'up',flag_use_overlap) 
   endif

   if(flag_exist_dn) then
     if(ispin .eq. 2) then ! we calculate dn(or beta)-spin part if coll. or noncol. case
       if(nspin .eq. 1) then  ! noncollinear
         psi_r_dn(igrid) = psi_r_dn(igrid) + psi_rho(phi_r, nbasis, ispin, V(:,1),SV(:,1), flag_plot_wavefunction, 'dn', flag_use_overlap) 
       elseif(nspin .eq. 2) then ! collinear
         psi_r_dn(igrid) = psi_r_dn(igrid) + psi_rho(phi_r, nbasis, ispin, V(:,2),SV(:,2), flag_plot_wavefunction, 'dn', flag_use_overlap) 
       endif
     endif
   endif

   return
endsubroutine
subroutine write_rho_main(pid_chg_up, pid_chg_dn, ngrid, nline, nwrite, nresi, psi_r_up, psi_r_dn, &
                          ispin, nspin, ispinor, flag_plot_wavefunction, flag_exist_up, flag_exist_dn)
   implicit none
   integer*4      ispin, ispinor, nspin
   integer*4      pid_chg_up, pid_chg_dn, ngrid, nline, nwrite, nresi
   complex*16     psi_r_up(ngrid), psi_r_dn(ngrid)
   logical        flag_plot_wavefunction
   logical        flag_exist_up, flag_exist_dn

   ! write rho (nonmag, noncol), rho_up (collinear), psi_up.majority (collin,noncollin)
   if(    ispinor .eq. 2 .and. .not. flag_plot_wavefunction) then   ! for noncol (rho)
     if(flag_exist_up) then
       call print_PARCHG_main(pid_chg_up, ngrid, nline, nwrite, nresi, psi_r_up+psi_r_dn, .false.)
     endif
   elseif(ispinor .eq. 2 .and.       flag_plot_wavefunction) then   ! for noncol (psi)
     if(flag_exist_up) then
       call print_PARCHG_main(pid_chg_up, ngrid, nline, nwrite, nresi, psi_r_up         , .true. )
       call print_PARCHG_main(pid_chg_dn, ngrid, nline, nwrite, nresi, psi_r_dn         , .true. )
     endif
   elseif(nspin   .eq. 2                                   ) then   ! for collin (rho or psi)
     if(flag_exist_up) then
       call print_PARCHG_main(pid_chg_up, ngrid, nline, nwrite, nresi, psi_r_up,flag_plot_wavefunction)
     endif
     if(flag_exist_dn) then
       call print_PARCHG_main(pid_chg_dn, ngrid, nline, nwrite, nresi, psi_r_dn,flag_plot_wavefunction)
     endif
   elseif(ispin   .eq. 1                                   ) then   ! for nonmag (rho or psi)
     if(flag_exist_up) then
       call print_PARCHG_main(pid_chg_up, ngrid, nline, nwrite, nresi, psi_r_up,flag_plot_wavefunction)
     endif
   endif


   return
endsubroutine
subroutine print_PARCHG_main(pid_chg_, ngrid, nline, nwrite, nresi, psi_r, flag_plot_wavefunction)
   implicit none
   integer*4      i, iline
   integer*4      ngrid, nline, nwrite, nresi
   integer*4      pid_chg_
   complex*16     psi_r(ngrid)
   logical        flag_plot_wavefunction

     do iline = 1, nline
       if(flag_plot_wavefunction) then
         write(pid_chg_,    '(5(1X,E17.10))') (real(psi_r((iline-1)*nwrite+i)),i=1,nwrite) ! real part
         write(pid_chg_+100,'(5(1X,E17.10))') (imag(psi_r((iline-1)*nwrite+i)),i=1,nwrite) ! imag part
       elseif(.not. flag_plot_wavefunction) then
         write(pid_chg_,    '(5(1X,E17.10))') (real(psi_r((iline-1)*nwrite+i)),i=1,nwrite) ! real part
       endif
     enddo

     if(nresi .ge. 1) then
       if(flag_plot_wavefunction) then
         write(pid_chg_,    '(5(1X,E17.10))') (real(psi_r((iline-1)*nwrite+i)),i=1,nresi) ! real part
         write(pid_chg_+100,'(5(1X,E17.10))') (imag(psi_r((iline-1)*nwrite+i)),i=1,nresi) ! imag part
         close(pid_chg_)
         close(pid_chg_+100)
       elseif(.not. flag_plot_wavefunction) then
         write(pid_chg_,    '(5(1X,E17.10))') (real(psi_r((iline-1)*nwrite+i)),i=1,nresi) ! real part
         close(pid_chg_)
       endif 
     endif

   return
endsubroutine
subroutine print_PARCHG_head(PINPT, PGEOM, ik, ie, pid_chg_up, pid_chg_dn, flag_sparse_exist_up, flag_sparse_exist_dn)
   use parameters, only : incar, poscar
   implicit none
   type(incar) :: PINPT
   type(poscar):: PGEOM
   integer*4      ik, ie, pid_chg_up, pid_chg_dn
   character*8    c_extension
   logical        flag_sparse_exist_up, flag_sparse_exist_dn
 
   flag_sparse_exist_up = .true. ! default
   flag_sparse_exist_dn = .true. ! default
   if(PINPT%flag_sparse) then
     if( PINPT%i_eig_print(ie) .le. PINPT%feast_ne(1, PINPT%i_kpt_print(ik))) then
       flag_sparse_exist_up = .true.
       if(PINPT%ispinor .eq. 2) flag_sparse_exist_dn = .true.
     else
       flag_sparse_exist_up = .false.
       if(PINPT%ispinor .eq. 2) flag_sparse_exist_dn = .false.
     endif
     if( PINPT%nspin .eq. 2) then
       if( PINPT%i_eig_print(ie) .le. PINPT%feast_ne(2, PINPT%i_kpt_print(ik))) then 
         flag_sparse_exist_dn = .true.
       else
         flag_sparse_exist_dn = .false.
       endif
     endif
   endif


   if(PINPT%flag_plot_wavefunction) then
     if(PINPT%flag_collinear .or. PINPT%flag_noncollinear) then
       if(flag_sparse_exist_up) then
         c_extension = '-real-up'; call PARCHG_head(pid_chg_up    ,ik, ie, PINPT,PGEOM, c_extension)
         c_extension = '-imag-up'; call PARCHG_head(pid_chg_up+100,ik, ie, PINPT,PGEOM, c_extension)
       endif
       if(flag_sparse_exist_dn) then
         c_extension = '-real-dn'; call PARCHG_head(pid_chg_dn    ,ik, ie, PINPT,PGEOM, c_extension)
         c_extension = '-imag-dn'; call PARCHG_head(pid_chg_dn+100,ik, ie, PINPT,PGEOM, c_extension)
       endif
     elseif(.not. PINPT%flag_collinear .and. .not. PINPT%flag_noncollinear) then
       if(flag_sparse_exist_up) then
         c_extension = '-real'; call PARCHG_head(pid_chg_up    ,ik, ie, PINPT,PGEOM, c_extension)
         c_extension = '-imag'; call PARCHG_head(pid_chg_up+100,ik, ie, PINPT,PGEOM, c_extension)
       endif
     endif

   elseif(.not. PINPT%flag_plot_wavefunction) then
     if(PINPT%flag_collinear) then
       if(flag_sparse_exist_up) then
         c_extension = '-up'; call PARCHG_head(pid_chg_up,ik, ie, PINPT, PGEOM, c_extension)
       endif
       if(flag_sparse_exist_dn) then
         c_extension = '-dn'; call PARCHG_head(pid_chg_dn,ik, ie, PINPT, PGEOM, c_extension)
       endif
     elseif(.not. PINPT%flag_collinear) then
       if(flag_sparse_exist_up) then
         c_extension = ''; call PARCHG_head(pid_chg_up,ik, ie, PINPT, PGEOM, c_extension)
       endif
     endif

   endif

   return
endsubroutine
subroutine PARCHG_head(pid_chg_,ik,ie,PINPT, PGEOM, c_extension)
  use parameters
  implicit none
  integer*4    i, ik, ie
  integer*4    pid_chg_
  character*40 fname
  character*8  c_extension
  type(incar)  :: PINPT
  type(poscar) :: PGEOM


  if(PINPT%i_eig_print(ie) .lt. 10) then

    if(PINPT%i_kpt_print(ik) .lt. 10) then
      write(fname,'(A,I1,A,I1,A)')"PARCHG-W-K000",PINPT%i_kpt_print(ik),"-E000",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 10 .and. PINPT%i_kpt_print(ik) .lt. 100)then
      write(fname,'(A,I2,A,I1,A)')"PARCHG-W-K00",PINPT%i_kpt_print(ik),"-E000",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 100 .and. PINPT%i_kpt_print(ik) .lt. 1000 ) then
      write(fname,'(A,I3,A,I1,A)')"PARCHG-W-K0",PINPT%i_kpt_print(ik),"-E000",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 1000 .and. PINPT%i_kpt_print(ik) .lt. 10000 ) then
      write(fname,'(A,I3,A,I1,A)')"PARCHG-W-K",PINPT%i_kpt_print(ik),"-E000",PINPT%i_eig_print(ie),trim(c_extension)
   
    endif

  elseif(PINPT%i_eig_print(ie) .ge. 10 .and. PINPT%i_eig_print(ie) .lt. 100)then

    if(PINPT%i_kpt_print(ik) .lt. 10) then
      write(fname,'(A,I1,A,I2,A)')"PARCHG-W-K000",PINPT%i_kpt_print(ik),"-E00",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 10 .and. PINPT%i_kpt_print(ik) .lt. 100)then
      write(fname,'(A,I2,A,I2,A)')"PARCHG-W-K00",PINPT%i_kpt_print(ik),"-E00",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 100 .and. PINPT%i_kpt_print(ik) .lt. 1000) then
      write(fname,'(A,I3,A,I2,A)')"PARCHG-W-K0",PINPT%i_kpt_print(ik),"-E00",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 1000 .and. PINPT%i_kpt_print(ik) .lt. 10000) then
      write(fname,'(A,I3,A,I2,A)')"PARCHG-W-K",PINPT%i_kpt_print(ik),"-E00",PINPT%i_eig_print(ie),trim(c_extension)
   
    endif

  elseif(PINPT%i_eig_print(ie) .ge. 100 .and. PINPT%i_eig_print(ie) .lt. 1000 ) then

    if(PINPT%i_kpt_print(ik) .lt. 10) then
      write(fname,'(A,I1,A,I3,A)')"PARCHG-W-K000",PINPT%i_kpt_print(ik),"-E0",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 10 .and. PINPT%i_kpt_print(ik) .lt. 100)then
      write(fname,'(A,I2,A,I3,A)')"PARCHG-W-K00",PINPT%i_kpt_print(ik),"-E0",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 100 .and. PINPT%i_kpt_print(ik) .lt. 1000) then
      write(fname,'(A,I3,A,I3,A)')"PARCHG-W-K0",PINPT%i_kpt_print(ik),"-E0",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 1000 .and. PINPT%i_kpt_print(ik) .lt. 10000) then
      write(fname,'(A,I3,A,I3,A)')"PARCHG-W-K",PINPT%i_kpt_print(ik),"-E0",PINPT%i_eig_print(ie),trim(c_extension)
   
    endif

  elseif(PINPT%i_eig_print(ie) .ge. 1000 .and. PINPT%i_eig_print(ie) .lt. 10000 ) then

    if(PINPT%i_kpt_print(ik) .lt. 10) then
      write(fname,'(A,I1,A,I3,A)')"PARCHG-W-K000",PINPT%i_kpt_print(ik),"-E",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 10 .and. PINPT%i_kpt_print(ik) .lt. 100)then
      write(fname,'(A,I2,A,I3,A)')"PARCHG-W-K00",PINPT%i_kpt_print(ik),"-E",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 100 .and. PINPT%i_kpt_print(ik) .lt. 1000) then
      write(fname,'(A,I3,A,I3,A)')"PARCHG-W-K0",PINPT%i_kpt_print(ik),"-E",PINPT%i_eig_print(ie),trim(c_extension)
   
    elseif(PINPT%i_kpt_print(ik) .ge. 1000 .and. PINPT%i_kpt_print(ik) .lt. 10000) then
      write(fname,'(A,I3,A,I3,A)')"PARCHG-W-K",PINPT%i_kpt_print(ik),"-E",PINPT%i_eig_print(ie),trim(c_extension)
   
    endif

  endif

  open(pid_chg_, file=fname, status='unknown')
  write(pid_chg_,'(A,I5,A,I5,A)')"WAVEFUNCTION: BAND= ",PINPT%i_eig_print(ie), &
                                            " ,KPOINT= ",PINPT%i_kpt_print(ik),trim(c_extension)
  write(pid_chg_,*)PGEOM%a_scale
  write(pid_chg_,'(3F20.16)')PGEOM%a_latt(1:3,1)
  write(pid_chg_,'(3F20.16)')PGEOM%a_latt(1:3,2)
  write(pid_chg_,'(3F20.16)')PGEOM%a_latt(1:3,3)
  write(pid_chg_,*)PGEOM%c_spec(:)
  write(pid_chg_,*)PGEOM%i_spec(:)
   write(pid_chg_,'(A)') "Direct"
  do i = 1, PGEOM%n_atom
    write(pid_chg_,'(3F20.16)') PGEOM%a_coord(:,i)+PGEOM%r_origin(:)
  enddo
  write(pid_chg_,*)" "
  write(pid_chg_,'(1x,3I6)')PGEOM%ngrid(1:3)

 return
endsubroutine

subroutine get_rxyz(rx,ry,rz, grid_a1, grid_a2, grid_a3, origin, nbasis, ngrid, a1, a2, a3, i1,i2,i3)
   implicit none
   integer*4    nbasis
   integer*4    ngrid(3),i1,i2,i3
   real*8       a1(3), a2(3), a3(3)
   real*8       grid_a1(0:ngrid(1)-1), rx(nbasis)
   real*8       grid_a2(0:ngrid(2)-1), ry(nbasis)
   real*8       grid_a3(0:ngrid(3)-1), rz(nbasis)
   real*8       origin(3,nbasis)  

   ! this routine calculates the relative distance to the orbital basis (at 'origin')

   rx(:)= grid_a1(i1)*a1(1) + grid_a2(i2)*a2(1) + grid_a3(i3)*a3(1) - origin(1,1:nbasis)
   ry(:)= grid_a1(i1)*a1(2) + grid_a2(i2)*a2(2) + grid_a3(i3)*a3(2) - origin(2,1:nbasis)
   rz(:)= grid_a1(i1)*a1(3) + grid_a2(i2)*a2(3) + grid_a3(i3)*a3(3) - origin(3,1:nbasis)

   return
endsubroutine

subroutine get_orbital_wavefunction_phi_r(phi_r, rx,ry,rz, corb, nbasis, rcut_orb_plot, flag_plot_wavefunction, &
                                          zeff, nqnum, lqnum, orb)
   use orbital_wavefunction, only : get_phi_r
   implicit none 
   integer*4    nbasis, iorb
   real*8       rx(nbasis)
   real*8       ry(nbasis)
   real*8       rz(nbasis)
   real*8       rcut_orb_plot
   character*8  corb(nbasis)  
   real*8        zeff(nbasis)
   real*8        nqnum(nbasis)
   integer*4     lqnum(nbasis)
   character*2   orb(nbasis)
   complex*16   phi_r(nbasis)
   logical      flag_plot_wavefunction

orb_: do iorb = 1, nbasis
        if( sqrt(rx(iorb)**2+ry(iorb)**2+rz(iorb)**2) .gt. rcut_orb_plot) then 
          phi_r(iorb) = (0d0,0d0)
          cycle orb_
        else
          phi_r(iorb) = get_phi_r(rx(iorb),ry(iorb),rz(iorb),corb(iorb), &
                                  zeff(iorb), orb(iorb), nqnum(iorb), lqnum(iorb) )
          if(.not. flag_plot_wavefunction) phi_r(iorb)= phi_r(iorb) * conjg(phi_r(iorb)) ! <phi(r)|phi(r)> if density plot is requested
        endif
      enddo orb_

   return
endsubroutine

subroutine set_variable_plot_eig(PINPT, PGEOM, neig, ngrid, nwrite, nline, nresi, pid_chg_up, pid_chg_dn, vol, &
                                 ng1, ng2, ng3, a1, a2, a3, origin, corb, grid_a1, grid_a2, grid_a3, &
                                 zeff, nqnum, lqnum, orb)
   use parameters, only : incar, poscar, pid_chg
   use element_info, only : angular
   use mpi_setup
   use print_io
   implicit none
   type(incar)  ::  PINPT
   type(poscar) ::  PGEOM
   integer*4        i, iorbital, iatom, iorb, mpierr
   integer*4        neig, ngrid, nwrite, nline, nresi
   integer*4        ng1, ng2, ng3
   integer*4        pid_chg_up, pid_chg_dn
   real*8           a1(3),a2(3),a3(3), a2xa3(3), vol
   real*8           rshift(3)
   character(*), parameter :: func = 'set_variable_plot_eig'
   character*8      corb(PGEOM%neig)
   real*8           zeff(PGEOM%neig)
   real*8           nqnum(PGEOM%neig)
   integer*4        lqnum(PGEOM%neig)
   character*2      orb(PGEOM%neig), orb_
   real*8           origin(3,PGEOM%neig)
   real*8           grid_d1, grid_a1(0:PGEOM%ngrid(1)-1)
   real*8           grid_d2, grid_a2(0:PGEOM%ngrid(2)-1)
   real*8           grid_d3, grid_a3(0:PGEOM%ngrid(3)-1)

   write(message,*)' '  ; write_msg
   if(PINPT%flag_plot_wavefunction) then
     write(message,'(A)')' #- START WRITING: EIGENSTATE WAVEFUNCTION'  ; write_msg
   elseif(.not. PINPT%flag_plot_wavefunction) then
     write(message,'(A)')' #- START WRITING: EIGENSTATE CHARGE DENSITY'  ; write_msg
   endif
   neig   = PGEOM%neig
   ngrid  = PGEOM%ngrid(1)*PGEOM%ngrid(2)*PGEOM%ngrid(3)
   ng1    = PGEOM%ngrid(1) ; ng2     = PGEOM%ngrid(2) ; ng3     = PGEOM%ngrid(3)
   nwrite = 5
   nline=int(ngrid/nwrite)
   nresi=mod(ngrid,nwrite)

   pid_chg_up = pid_chg
   pid_chg_dn = pid_chg + 10

   rshift(1:3)=PGEOM%r_origin(1:3)
   a1=PGEOM%a_latt(1:3,1)
   a2=PGEOM%a_latt(1:3,2)
   a3=PGEOM%a_latt(1:3,3)
   call vcross(a2xa3,a2,a3)
   vol=dot_product(a1,a2xa3)


   iorbital = 0
   do iatom = 1, PGEOM%n_atom
     do iorb = 1, PGEOM%n_orbital(iatom)
       iorbital = iorbital + 1
       origin(:,iorbital) =( PGEOM%a_coord(1,iatom) + rshift(1) )*a1(:) + &
                           ( PGEOM%a_coord(2,iatom) + rshift(2) )*a2(:) + &
                           ( PGEOM%a_coord(3,iatom) + rshift(3) )*a3(:)
       corb(iorbital)=trim(PGEOM%c_orbital(iorb,iatom))

       zeff(iorbital) = PGEOM%z_eff_nuc(iorb, iatom)
       nqnum(iorbital)= PGEOM%n_quantum(iatom)
       lqnum(iorbital)= PGEOM%l_quantum(iorb,iatom)
       write(orb(iorbital),'(I0,A)') PGEOM%orb_n_quantum(iorb,iatom), angular(PGEOM%l_quantum(iorb,iatom))
       !write(orb_,'(I0,A)') PGEOM%orb_n_quantum(iorb,iatom), PGEOM%l_quantum(iorb,iatom)
       !orb(iorbital) = orb_

     enddo
   enddo

   if (iorbital .ne. PGEOM%neig) then
     write(message,'(A,A)')'  !WARNING! iorbital is not same as neig!, please check again. ',func  ; write_msg
     kill_job
   endif

   grid_d1 = 1d0/dble(ng1)
   grid_d2 = 1d0/dble(ng2)
   grid_d3 = 1d0/dble(ng3)

   grid_a1 = (/(dble(i)*grid_d1, i=0, ng1-1)/)
   grid_a2 = (/(dble(i)*grid_d2, i=0, ng2-1)/)
   grid_a3 = (/(dble(i)*grid_d3, i=0, ng3-1)/)

   return
endsubroutine

subroutine write_info_plot_eig(PINPT)
   use parameters, only : incar
   use mpi_setup
   use print_io
   implicit none
   type(incar)  ::  PINPT

   if(PINPT%flag_plot_wavefunction) then
     write(message,'(A)')      "  PLOT MODE: WAV_PLOT= .TRUE."  ; write_msg
     write(message,'(A)')      "             The real and imaginary part of wavefunction will be plotted in "  ; write_msg
     write(message,'(A)')      "             separate file with '-real' and '-imag' extension, respectively."  ; write_msg
     write(message,'(A,*(I5))')"             Selected band index: ", PINPT%i_eig_print(1:PINPT%n_eig_print)  ; write_msg
     write(message,'(A,*(I5))')"             Selected   k-points: ", PINPT%i_kpt_print(1:PINPT%n_kpt_print)  ; write_msg
   elseif(.not.PINPT%flag_plot_wavefunction) then
     write(message,'(A)')      "  PLOT MODE: CHG_PLOT= .TRUE."  ; write_msg
     write(message,'(A)')      "             The charge density plot is requested."  ; write_msg
     write(message,'(A,*(I5))')"             Selected band index: ", PINPT%i_eig_print(1:PINPT%n_eig_print)  ; write_msg
     write(message,'(A,*(I5))')"             Selected   k-points: ", PINPT%i_kpt_print(1:PINPT%n_kpt_print)  ; write_msg
   endif

   return
endsubroutine

subroutine reset_orbital_origin(origin_, origin, neig, a1, a2, a3, ix, iy, iz)
   implicit none
   integer*4    neig
   integer*4    ix,iy,iz
   real*8       origin(3,neig), origin_(3,neig)
   real*8       rshift_orbital(3)  
   real*8       a1(3),a2(3),a3(3)

   rshift_orbital(:) = ix * a1(:) + iy * a2(:) + iz * a3(:)
   origin_(1,1:neig) = origin(1,1:neig) + rshift_orbital(1)
   origin_(2,1:neig) = origin(2,1:neig) + rshift_orbital(2)
   origin_(3,1:neig) = origin(3,1:neig) + rshift_orbital(3)

   return
endsubroutine

subroutine check_isnan(psi_r_up, psi_r_dn, ngrid, nspin)
   use parameters, only : zi
   implicit none
   integer*4    i
   integer*4    ngrid, nspin
   complex*16   psi_r_up(ngrid)
   complex*16   psi_r_dn(ngrid)

   ! note that 'isnan' operator is called to check whether there are singlular point.
   ! note also that 'isnan' function is valid only for 'intel' compiler.
   ! Therefore, one need to be imporoved to use another simpler or general operator to check it.

   do i = 1, ngrid
     if(isnan(real(psi_r_up(i))) .and. .not. isnan(imag(psi_r_up(i)))) then 
      psi_r_up(i) = 0d0 + imag(psi_r_up(i)) * zi
     elseif(.not. isnan(real(psi_r_up(i)))  .and. isnan(imag(psi_r_up(i)))) then
      psi_r_up(i) = real(psi_r_up(i)) + 0d0 * zi
     elseif(isnan(real(psi_r_up(i))) .and. isnan(imag(psi_r_up(i)))) then
      psi_r_up(i) = (0d0,0d0)
     endif
   enddo

   if(nspin .eq. 2) then
     do i = 1, ngrid
       if(isnan(real(psi_r_dn(i))) .and. .not. isnan(imag(psi_r_dn(i)))) then
        psi_r_dn(i) = 0d0 + imag(psi_r_dn(i)) * zi
       elseif(.not. isnan(real(psi_r_dn(i)))  .and. isnan(imag(psi_r_dn(i)))) then
        psi_r_dn(i) = real(psi_r_dn(i)) + 0d0 * zi
       elseif(isnan(real(psi_r_dn(i))) .and. isnan(imag(psi_r_dn(i)))) then
        psi_r_dn(i) = (0d0,0d0)
       endif
     enddo
   endif

   return
endsubroutine
