#include "alias.inc"
module total_energy
   use parameters, only: energy, incar, kpoints, poscar
   use mpi_setup
   use sorting, only : get_sort_index_1D
   use do_math, only : idx2Di, idx2Dj
   use print_io
contains

   subroutine get_total_energy(ETBA, PINPT, PKPTS, PGEOM, flag_opt_ef, flag_report, E_F_)
      implicit none
      type(energy)    :: ETBA
      type(incar)     :: PINPT
      type(poscar)    :: PGEOM
      type(kpoints)   :: PKPTS
      integer*4          nkp, nband, nspin  
      integer*4          mpierr
      real*8             eltemp
      real*8             nelect_ref
      real*8             degen ! degeneracy
      real*8             E_F, E0
      real*8             OCC(PGEOM%nband*PINPT%nspin,PKPTS%nkpoint)
      logical            flag_quit, flag_opt_ef, flag_report
      real*8,intent(inout),optional :: E_F_

      if(.not. PINPT%flag_get_total_energy) return

      if(flag_report) then
        write(message,'(A)')' ' ; write_msg
        write(message,'(A)')' #--- START: TOTAL ENERGY EVALUATION -----------' ; write_msg
      endif

      ! initialize
      ETBA%E_BAND = 0d0
      ETBA%E_TOT  = 0d0
      degen       = 1d0
      ETBA%E_F    = 0d0 ! temporary

      ! if present E_F_, use it and return it (if flag_opt_ef -> update) with E_F_, 
      ! otherwise, use ETBA%E_F and (if flag_opt_ef -> update) save it to ETBA%E_F
      if(present(E_F_)) then
        E_F     = E_F_
      else
        E_F     = ETBA%E_F
      endif

      eltemp      = PINPT%electronic_temperature ! default = 0k
      nkp         = PKPTS%nkpoint
      nband       = PGEOM%nband
      nspin       = PINPT%nspin
      nelect_ref  = sum(PGEOM%nelect)   ! reference number of electrons
      if(PINPT%ispin .eq. 1) degen = 2d0


      ! optimize Fermi level with given number of electrons nelect_ref
      if(flag_opt_ef) then
        call get_fermi_level(ETBA%E, E_F, nelect_ref, nband, nspin, nkp, eltemp, degen)
        if(present(E_F_)) then 
          E_F_     = E_F
        else
          ETBA%E_F = E_F
        endif
      endif

      ! evaluate occupation with Fermi level E_F, F_OCC
      call get_occupation(ETBA%E, E_F, ETBA%F_OCC, nband, nspin, nkp, degen, eltemp)
      ! evaluate band energy based on F_OCC
      ETBA%E_BAND = sum(band_energy(ETBA%E, ETBA%F_OCC))/nkp

      ! evaluate band energy based on F_OCC with eltemp = 0d0
      call get_occupation(ETBA%E, E_F, OCC, nband, nspin, nkp, degen, 0d0)
      E0 = sum(band_energy(ETBA%E, OCC))/nkp
        
      if(flag_report) then
        write(message,'(A,F20.8, A)')      '  -Fermi level E_F (in eV)                    : ', E_F ; write_msg
        write(message,'(A,F20.8, A)')      '  -Number of electrons (with E_F)             : ', sum(ETBA%F_OCC)/nkp ; write_msg
        write(message,'(A,F20.8, A)')      '  -Number of electrons (NELECT)               : ', nelect_ref          ; write_msg
        write(message,'(A)')               '  -Total energy components (in eV)              '              ; write_msg
        write(message,'(A,F10.4, A,F20.8)')'     *- Band energy (ELTEMP =  ',eltemp,')    : ', ETBA%E_BAND     ; write_msg
        write(message,'(A,F10.4, A,F20.8)')'     *- Band energy (ELTEMP -> ',0.0d0 ,')    : ', E0          ; write_msg
        write(message,'(A      )')         '     __________________________________________' ; write_msg
        write(message,'(A,F20.8, A)')      '     #- Total energy (in eV)                  : ', ETBA%E_BAND     ; write_msg
        
        write(message,'(A)')' #--- END: TOTAL ENERGY EVALUATION -----------' ; write_msg
        write(message,'(A)')' ' ; write_msg
      endif

      return
   endsubroutine

   subroutine get_fermi_level(E, E_F, nelect_ref, nband, nspin, nkp, eltemp, degen)
      implicit none
      integer*4          ii, ik, ie
      integer*4          nkp, nband, nspin
      integer*4          mpierr
      real*8             eltemp, de_cut
      real*8             nelect
      real*8             nelect_ref
      real*8             grad_nelect
      real*8             degen ! degeneracy
      real*8             E(nband*nspin,nkp)
      real*8             E_F
      logical            flag_quit

      nelect      = 0d0
      flag_quit   = .false.
      grad_nelect = 0d0
      de_cut      = 0.0000001 ! cutoff for the convergence of Fermi level evaluation

      !step0. Rough estimation of Fermi level with mid-gap position
      call get_midgap(E, E_F, nelect_ref, nband, nspin, nkp, degen)

      !step1. get nelect with initial E_F = 0d0
      call get_nelect(E, E_F, nelect, nband, nspin, nkp, eltemp, degen)
      if( abs(nelect - nelect_ref) .lt. de_cut) then
        flag_quit = .true.
      else
        flag_quit = .false.
        call get_grad_nelect(E, E_F, grad_nelect, nband, nspin, nkp, eltemp, degen)
        E_F = E_F - (nelect - nelect_ref)/grad_nelect 
        call get_nelect(E, E_F, nelect, nband, nspin, nkp, eltemp, degen)
      endif

      !step2. if nelect = nelec_ref -> quit, otherwise, recalculate with updated E_F
      !       based on gradient of nelect w.r.t chemical potential
      ii = 0
      do while (.not. flag_quit)
         ii  = ii + 1
        if( abs(nelect - nelect_ref) .lt. de_cut ) then
          flag_quit = .true.
        else
          flag_quit = .false.
          call get_grad_nelect(E, E_F, grad_nelect, nband, nspin, nkp, eltemp, degen)
          E_F = E_F - (nelect - nelect_ref)/grad_nelect
          call get_nelect(E, E_F, nelect, nband, nspin, nkp, eltemp, degen)
        endif
      enddo

      return
   endsubroutine

   elemental real*8 function band_energy(E, occ)
      implicit none
      real*8, intent(in) :: E, occ

      band_energy = E * occ

      return
   endfunction

   subroutine get_occupation(E, E_F, occ, nband, nspin, nkp, degen, eltemp)
      use do_math, only : fdirac
      implicit none
      integer*4    nband, nspin, nkp
      integer*4    ik, ie, is
      integer*4    mpierr
      real*8       E_F
      real*8       eltemp, degen
      real*8       E(nband*nspin, nkp)
      real*8       occ(nband * nspin, nkp), occ_(nband * nspin, nkp)
      integer*4    ourjob(nprocs), ourjob_disp(0:nprocs-1)

      occ_ = 0d0

      call mpi_job_distribution_chain(nkp, nprocs, ourjob, ourjob_disp)

      do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        do is = 1, nspin
          occ_(1+nband*(is-1):nband*is,ik) = fdirac( E(1+nband*(is-1):nband*is,ik), E_F, eltemp ) * degen
        enddo
      enddo
 
#ifdef MPI
      call MPI_ALLREDUCE(occ_, occ, size(occ), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
#else
      occ = occ_
#endif

      return
   endsubroutine


   subroutine get_nelect(E, E_F, nelect, nband, nspin, nkp, eltemp, degen)
      use do_math, only : fdirac
      implicit none
      integer*4    nband, nspin, nkp
      integer*4    ik, ie, is
      integer*4    mpierr
      real*8       dos_tot
      real*8       eltemp, degen
      real*8       E(nband*nspin, nkp)
      real*8       nelect
      real*8       E_F
      integer*4    ourjob(nprocs), ourjob_disp(0:nprocs-1)

      call mpi_job_distribution_chain(nkp, nprocs, ourjob, ourjob_disp)

      dos_tot = 0d0

      do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        do is = 1, nspin
          dos_tot     = dos_tot     + sum( fdirac( E(1+nband*(is-1):nband*is,ik), E_F, eltemp ) ) * degen
        enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE( dos_tot/real(nkp), nelect, 1, MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
#else
      nelect = dos_tot/real(nkp)
#endif

      return
   endsubroutine

   subroutine get_grad_nelect(E, E_F, grad_nelect, nband, nspin, nkp, eltemp, degen)
      use do_math, only : grad_u_fdirac
      implicit none
      integer*4    nband, nspin, nkp
      integer*4    ik, ie, is
      integer*4    mpierr
      real*8       eltemp, degen
      real*8       E(nband*nspin, nkp)
      real*8       grad_nelect, grad_nelect_
      real*8       E_F
      integer*4    ourjob(nprocs), ourjob_disp(0:nprocs-1)

      call mpi_job_distribution_chain(nkp, nprocs, ourjob, ourjob_disp)

      grad_nelect  = 0d0

      do ik = sum(ourjob(1:myid))+1, sum(ourjob(1:myid+1))
        do is = 1, nspin
          grad_nelect = grad_nelect + sum( grad_u_fdirac( E(1+nband*(is-1):nband*is,ik), E_F, eltemp ) ) * degen
        enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(grad_nelect, grad_nelect_, 1, MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
      grad_nelect = grad_nelect_ / real(nkp)
#else
      grad_nelect = grad_nelect / real(nkp)
#endif

      return
   endsubroutine

   subroutine get_midgap(E, E_F, nelect_ref, nband, nspin, nkp, degen)
      use sorting, only : get_sort_index_1D
      use do_math, only : idx2Di, idx2Dj
      implicit none
      integer*4             ii, ie, ik
      integer*4             nband, nspin, nkp, neig_tot
      integer*4             idx_eig1D(nband*nspin*nkp)
      real*8, intent(in) :: E(nband*nspin, nkp)
      real*8                nelect_ref
      real*8                nelect
      real*8                E1, E2
      real*8                E_F
      real*8                degen

      neig_tot    = nband*nspin*nkp
      nelect      = 0d0
      ii          = 0

      ! sort eigenvalue in ascending order
      call get_sort_index_1D(idx_eig1D, reshape(E,[neig_tot]), neig_tot, 'ascending ') 

      do while (nelect .lt. nelect_ref)
        ii = ii + 1
        nelect = nelect + 1d0/real(nkp) * degen
        ie = idx2Di(idx_eig1D(ii),nband*nspin)
        ik = idx2Dj(idx_eig1D(ii),nband*nspin)
      enddo

      ! VBM     
      E1 = E(ie, ik)

      ! CBM
      ie = idx2Di(idx_eig1D(ii),nband*nspin)
      ik = idx2Dj(idx_eig1D(ii),nband*nspin)
      E2 = E(ie, ik)

      E_F = (E1 + E2)*0.5d0
 
      return
   endsubroutine
endmodule
