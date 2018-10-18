#include "alias.inc"
subroutine get_z2(NN_TABLE, PINPT, PINPT_BERRY, PGEOM, PKPTS)
   use parameters, only : incar, hopping, poscar, berry, energy, kpoints, pi2, cyclic_axis
   use berry_phase
   use mpi_setup
   use time
   implicit none
   type(hopping) :: NN_TABLE
   type(incar)   :: PINPT
   type(berry)   :: PINPT_BERRY
   type(poscar)  :: PGEOM
   type(kpoints) :: PKPTS
   integer*4        mpierr
   logical          flag_phase, flag_sparse, flag_get_chern
   integer*4        i, is
   integer*4        nerange, nkdiv, nkpath, nplane
   integer*4        ix, ip, ikpath
   integer*4        iband, nband
   integer*4        erange(PINPT_BERRY%z2_nerange)
   real*8           E(PINPT_BERRY%z2_nerange,PINPT_BERRY%z2_nkdiv)
   complex*16       V(PGEOM%neig*PINPT%ispin,PINPT_BERRY%z2_nerange,PINPT_BERRY%z2_nkdiv)
   real*8           shift(2)
   real*8           k_init_reci(3), k_end_reci(3), kpath(3,2,PINPT_BERRY%z2_nkpath)
   real*8           G(3)
   real*8           time1, time2
   character*2      caxis(3)
   character*3      cplane(2)
   character*40     fnm_header, fnm_gap_header
   character*40     z2_filenm, z2_gap_filenm
   integer*4        z2_axis(size(PINPT_BERRY%z2_axis))
   real*8           largest_gap(PINPT%nspin,PINPT_BERRY%z2_nkpath)
   integer*4        clock_direct (PINPT%nspin,PINPT_BERRY%z2_nkpath)
   integer*4        z2_index(PINPT%nspin), z2_bulk(0:3,PINPT%nspin)
   integer*4        z2_dimension

   allocate(PINPT_BERRY%z2_wcc(PINPT_BERRY%z2_nerange/PINPT%nspin,PINPT%nspin,PINPT_BERRY%z2_nkpath, &
                               PINPT_BERRY%z2_nplane,size(PINPT_BERRY%z2_axis)))
   PINPT_BERRY%z2_wcc = 0d0
   if(PINPT_BERRY%flag_z2_get_chern) then
     allocate(PINPT_BERRY%z2_chern(PINPT%nspin,PINPT_BERRY%z2_nplane,size(PINPT_BERRY%z2_axis)))
     allocate(PINPT_BERRY%z2_polarization(PINPT%nspin,PINPT_BERRY%z2_nkpath,PINPT_BERRY%z2_nplane,size(PINPT_BERRY%z2_axis)))
     PINPT_BERRY%z2_chern = 0d0
     PINPT_BERRY%z2_polarization = 0d0
   endif

   call time_check(time1, time2, 'init')

   if_main write(6,*)''
   if_main write(6,'(A)')'START: Z2 EVALUATION'
   if_main write(6,'(A,A)')'  BAND INDEX: ',adjustl(trim(PINPT_BERRY%strip_z2_range))
   flag_sparse= .false.
   flag_phase = PINPT_BERRY%flag_z2_phase ! default = .false.
!  flag_phase = .TRUE. 
   flag_get_chern = PINPT_BERRY%flag_z2_get_chern
   if(size(PINPT_BERRY%z2_axis) .eq. 3) nplane = 2
   if(size(PINPT_BERRY%z2_axis) .eq. 1) nplane = 1
   caxis(1) ='B1' ;caxis(2)='B2';caxis(3)='B3'
   cplane(1)='0.0';cplane(2)='0.5'
   erange         = PINPT_BERRY%z2_erange(:)
   nerange        = PINPT_BERRY%z2_nerange
   nband          = PINPT_BERRY%z2_nerange/PINPT%nspin 
   nkdiv          = PINPT_BERRY%z2_nkdiv
   nkpath         = PINPT_BERRY%z2_nkpath
   z2_axis        = PINPT_BERRY%z2_axis
   z2_dimension   = PINPT_BERRY%z2_dimension
   shift(:)       = (/0d0,0.5d0/)
   fnm_header     = PINPT_BERRY%z2_filenm
   fnm_gap_header = PINPT_BERRY%z2_gap_filenm
   z2_bulk        = 0
   iband          = erange(1)
   if(z2_dimension .ge. 2) then
axis:do ix = 1, size(z2_axis)
 plane:do ip = 1, nplane
         call set_kpath_plane(PINPT_BERRY%z2_kpoint(:,:,:,ip,ix), PINPT_BERRY%z2_kpoint_reci(:,:,:,ip,ix), kpath, nkdiv, nkpath, z2_axis(ix), shift(ip), PGEOM)
         G = PINPT_BERRY%z2_kpoint(:,nkdiv,1,ip,ix) - PINPT_BERRY%z2_kpoint(:,1,1,ip,ix)
    path:do ikpath = 1, nkpath
           call get_eig(NN_TABLE, PINPT_BERRY%z2_kpoint(:,:,ikpath,ip,ix), nkdiv, PINPT, E, V, PGEOM%neig, iband, nband, .true., flag_sparse, .false., flag_phase)
           call set_periodic_gauge(V, G, PINPT, PGEOM, nkdiv, erange, nerange)
#ifdef F08
           call get_berry_phase    (PINPT_BERRY%z2_wcc(:,:,ikpath,ip,ix), PINPT_BERRY%z2_kpoint(:,:,ikpath,ip,ix), V, PINPT, PGEOM, nkdiv, erange, nerange)
#else
           call get_berry_phase_svd(PINPT_BERRY%z2_wcc(:,:,ikpath,ip,ix), PINPT_BERRY%z2_kpoint(:,:,ikpath,ip,ix), V, PINPT, PGEOM, nkdiv, erange, nerange)
#endif
           if_main call write_status(ikpath, nkpath, z2_axis(ix), ip)
         enddo path
         call get_z2_fname(z2_filenm    , PINPT_BERRY%z2_filenm    , ip, z2_axis(ix))
         call get_z2_fname(z2_gap_filenm, PINPT_BERRY%z2_gap_filenm, ip, z2_axis(ix))

#ifdef MPI
         ! MPI routine is not supported yet...
         !  call MPI_Reduce(wcc, wcc_, size(wcc), MPI_REAL8, MPI_SUM, 0, mpi_comm_earth, mpierr)
         !  PINPT_BERRY%wcc = wcc_
         call find_largest_gap(largest_gap, clock_direct, z2_index, PINPT_BERRY%z2_wcc(:,:,:,ip,ix), PINPT%nspin, nkpath, nerange)
         call store_z2_bulk_index(z2_bulk, z2_index, z2_axis(ix), ip, PINPT%nspin, z2_dimension)
         if(flag_get_chern) call get_chern_number(PINPT_BERRY%z2_chern(:,ip,ix), PINPT_BERRY%z2_polarization(:,:,ip,ix), &
                                                  PINPT_BERRY%z2_wcc(:,:,:,ip,ix), PINPT%nspin, nkpath, nerange)
         if_main call print_wcc(PINPT, z2_filenm, z2_gap_filenm, PINPT_BERRY%z2_wcc(:,:,:,ip,ix), kpath, nkpath, &
                                nerange, PINPT_BERRY%strip_z2_range, largest_gap, clock_direct, z2_index, &
                                PINPT_BERRY%z2_polarization(:,:,ip,ix), PINPT_BERRY%z2_chern(:,ip,ix), flag_get_chern)
#else
!        PINPT_BERRY%wcc = wcc
         call find_largest_gap(largest_gap, clock_direct, z2_index, PINPT_BERRY%z2_wcc(:,:,:,ip,ix), PINPT%nspin, nkpath, nerange)
         call store_z2_bulk_index(z2_bulk, z2_index, z2_axis(ix), ip, PINPT%nspin, z2_dimension)
         if(flag_get_chern) call get_chern_number(PINPT_BERRY%z2_chern(:,ip,ix), PINPT_BERRY%z2_polarization(:,:,ip,ix), &
                                                  PINPT_BERRY%z2_wcc(:,:,:,ip,ix), PINPT%nspin, nkpath, nerange)
         call print_wcc(PINPT, z2_filenm, z2_gap_filenm, PINPT_BERRY%z2_wcc(:,:,:,ip,ix), kpath, nkpath, &
                        nerange, PINPT_BERRY%strip_z2_range, largest_gap, clock_direct, z2_index, &
                        PINPT_BERRY%z2_polarization(:,:,ip,ix), PINPT_BERRY%z2_chern(:,ip,ix), flag_get_chern)
#endif
       enddo plane
     enddo axis

   else
     if_main write(6,'(A)')'  !WARN! Current version does not support 1D case... Exit anyway.'
     stop
   endif

   call time_check(time1, time2)
   
   if_main_then
     if(.not. flag_get_chern) call write_z2_index(PINPT%nspin, z2_dimension, z2_bulk)
     write(6,'(A,F12.3)')'END: Z2 INDEX CALCULATION. TIME ELAPSED (s) =',time1
   if_main_end

   return
endsubroutine
subroutine write_z2_index(nspin, z2_dimension, z2_bulk)
   implicit none
   integer*4    nspin
   integer*4    z2_dimension
   integer*4    z2_bulk(0:3,nspin)
   integer*4    is
   character*2  cspin(2)
   cspin(1) = 'up'
   cspin(2) = 'dn'

   do is = 1, nspin
     if(nspin .eq. 2) then
       write(6,'(A,A)',ADVANCE='no')'  For spin-',cspin(is)
     endif
     if(z2_dimension .eq. 3) then
       z2_bulk(0,is) = mod(z2_bulk(0,is),2)
       write(6,'(A,4(I1,A))')'  BULK TOPOLOGICAL INDEX: [v0; v1, v2, v3] = [',z2_bulk(0,is),'; ',z2_bulk(1,is),', ', &
                                                                     z2_bulk(2,is),', ',z2_bulk(3,is),']'
     elseif(z2_dimension .le. 2) then
       write(6,'(A,I1)')'  Z2 INDEX:  v = ',z2_bulk(0,is)
     endif
   enddo

   return
endsubroutine
subroutine write_status(ikpath, nkpath, iaxis, ip)
   use parameters, only : cyclic_axis
   implicit none
   integer*4    ip, iaxis
   integer*4    ikpath, nkpath
   character*2  caxis(3)
   character*3  cplane(2)
   caxis(1) ='B1' ;caxis(2)='B2';caxis(3)='B3'
   cplane(1)='0.0';cplane(2)='0.5'
 
   write(6,'(2(A,I3),8A)')'  STATUS: ',ikpath,'/',nkpath,' KPATH, [',  &
                          caxis(cyclic_axis(iaxis,1)),'-',caxis(cyclic_axis(iaxis,2)), &
                          '] K-PLANE with ', cplane(ip),' ',caxis(iaxis)
  
   return
endsubroutine
subroutine get_z2_fname(fname, fnm_header, ip, iaxis)
   use parameters, only : cyclic_axis
   implicit none
   integer*4    ip, iaxis
   character*2  caxis(3)
   character*3  cplane(2)
   character*40 fname, fnm_header

   caxis(1) ='B1' ;caxis(2)='B2';caxis(3)='B3'
   cplane(1)='0.0';cplane(2)='0.5'

!  write(fname,'(10A)')trim(fnm_header),'.',cplane(ip),'-',caxis(iaxis),&
!        '.',caxis(cyclic_axis(iaxis,1)),'_',caxis(cyclic_axis(iaxis,2)),'-PLANE.dat'
!  write(fname,'(6A)')trim(fnm_header),'.',cplane(ip),'-',caxis(iaxis),'-PLANE.dat'
   write(fname,'(6A)')trim(fnm_header),'.',cplane(ip),'-',caxis(iaxis),'.dat'

   return
endsubroutine
subroutine store_z2_bulk_index(z2_bulk, z2_index, iaxis, iplane, nspin, z2_dimension)
   implicit none
   integer*4    nspin, iplane, iaxis
   integer*4    z2_index(nspin)
   integer*4    z2_bulk(0:4, nspin)
   integer*4    z2_dimension

   if(    z2_dimension .eq. 3 .and. iplane .eq. 2 .and. iaxis .ne. 3) then ! store kx1, ky1
     z2_bulk(iaxis,:) =                z2_index(:)
   elseif(z2_dimension .eq. 3 .and. iplane .eq. 2 .and. iaxis .eq. 3) then ! store kz1 and add to kz0
     z2_bulk(iaxis,:) =                z2_index(:)
     z2_bulk(0    ,:) = z2_bulk(0,:) + z2_index(:)
   elseif(z2_dimension .eq. 3 .and. iplane .eq. 1 .and. iaxis .eq. 3) then ! store kz0 and add to kz0
     z2_bulk(0    ,:) = z2_bulk(0,:) + z2_index(:)
   elseif(z2_dimension .le. 2) then
     z2_bulk(0    ,:) =                z2_index(:)
   endif

   return
endsubroutine
