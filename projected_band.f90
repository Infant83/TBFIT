#include "alias.inc"
module projected_band
   use do_math
   use mpi_setup
   use print_io

contains

   subroutine get_orbital_projection(ETBA, PKPTS, PINPT, PGEOM, flag_use_overlap)
     use parameters, only: energy, incar, poscar, kpoints
     implicit none
     type(incar  ) :: PINPT
     type(poscar ) :: PGEOM
     type(kpoints) :: PKPTS
     type(energy ) :: ETBA 
     integer*4        neig, iband, nband
     integer*4        mpierr
     integer*4        sizebuff
     integer*4        ik, is, my_ik
     integer*4        ie, im
     real*8           ORB(PINPT%lmmax,PGEOM%nband*PINPT%nspin, PKPTS%nkpoint)
     complex*16       c_up, c_dn
     complex*16       sc_up, sc_dn
     complex*16,allocatable :: myV(:,:,:)
     complex*16,allocatable :: mySV(:,:,:)
     integer*4, allocatable :: ourjob(:), ourjob_disp(:)
     integer*4                 ncpu, id
     logical                   flag_use_overlap

     if(.not. PINPT%flag_fit_orbital) return

     sizebuff = PGEOM%neig*PINPT%ispinor*PGEOM%nband*PINPT%nspin
     ETBA%ORB = 0d0

     if(COMM_KOREA%flag_split) then
       ncpu = COMM_KOREA%nprocs
       id   = COMM_KOREA%myid
     else
       ncpu = nprocs
       id   = myid
     endif
     allocate(ourjob(ncpu))
     allocate(ourjob_disp(0:ncpu-1))
     call mpi_job_distribution_chain(PKPTS%nkpoint, ncpu, ourjob, ourjob_disp)

     allocate(myV(PGEOM%neig*PINPT%ispinor,PGEOM%nband*PINPT%nspin,ourjob(id+1)))
     if(flag_use_overlap) then
       allocate(mySV(PGEOM%neig*PINPT%ispinor,PGEOM%nband*PINPT%nspin,ourjob(id+1)))
     endif
#ifdef MPI
     if(.not.COMM_KOREA%flag_split) then
       call MPI_SCATTERV(ETBA%V, ourjob*sizebuff, ourjob_disp*sizebuff, &
                         MPI_COMPLEX16, myV, ourjob(id+1)*sizebuff, &
                         MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
       if(flag_use_overlap) then
         call MPI_SCATTERV(ETBA%SV, ourjob*sizebuff, ourjob_disp*sizebuff, &
                           MPI_COMPLEX16, mySV, ourjob(id+1)*sizebuff, &
                           MPI_COMPLEX16, 0, mpi_comm_earth, mpierr)
       endif
     elseif(COMM_KOREA%flag_split) then
       call MPI_SCATTERV(ETBA%V, ourjob*sizebuff, ourjob_disp*sizebuff, &
                         MPI_COMPLEX16, myV, ourjob(id+1)*sizebuff, &
                         MPI_COMPLEX16, 0, COMM_KOREA%mpi_comm, mpierr)
       if(flag_use_overlap) then
         call MPI_SCATTERV(ETBA%SV, ourjob*sizebuff, ourjob_disp*sizebuff, &
                           MPI_COMPLEX16, mySV, ourjob(id+1)*sizebuff, &
                           MPI_COMPLEX16, 0, COMM_KOREA%mpi_comm, mpierr)
       endif
     endif
#else
     myV=ETBA%V
     if(flag_use_overlap) then
       mySV=ETBA%SV
     endif
#endif

     do is = 1, PINPT%nspin
     do ik = sum(ourjob(1:id))+1, sum(ourjob(1:id+1))
       my_ik = ik - sum(ourjob(1:id))
       do ie=1, PGEOM%nband

         ! NOTE: here we use Mulliken charge if flag_use_overlap = .true., not net charge
         do im=1, PGEOM%neig
           if(PINPT%ispinor .eq. 2) then
             c_up = myV(im,ie,my_ik); c_dn = myV(im+PGEOM%neig,ie,my_ik)
             if(flag_use_overlap) then
               sc_up = mySV(im,ie,my_ik); sc_dn = mySV(im+PGEOM%neig,ie,my_ik)
             else
               sc_up = c_up             ; sc_dn = c_dn
             endif

             ETBA%ORB(PGEOM%orb_index(im), ie, ik) = &
                ETBA%ORB(PGEOM%orb_index(im), ie, ik) + &
                real( conjg(c_up)*sc_up + conjg(c_dn)*sc_dn)

           else
             c_up = myV(im+PGEOM%neig*(is-1),ie+PGEOM%nband*(is-1),my_ik)
             if(flag_use_overlap) then
               sc_up = myV(im+PGEOM%neig*(is-1),ie+PGEOM%nband*(is-1),my_ik)
             else
               sc_up = c_up
             endif

             ETBA%ORB(PGEOM%orb_index(im), ie+PGEOM%nband*(is-1), ik) = &
                ETBA%ORB(PGEOM%orb_index(im), ie+PGEOM%nband*(is-1), ik) + &
                real(conjg(c_up)*sc_up)

           endif
         enddo

       enddo
     enddo
     enddo

#ifdef MPI
!    call MPI_ALLREDUCE(ETBA%ORB, ORB, size(ORB), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
!    ETBA%ORB = ORB

     if(.not. COMM_KOREA%flag_split) then
       call MPI_ALLREDUCE(ETBA%ORB, ORB, size(ORB), MPI_REAL8, MPI_SUM, mpi_comm_earth, mpierr)
       ETBA%ORB = ORB
     elseif(COMM_KOREA%flag_split) then
       call MPI_ALLREDUCE(ETBA%ORB, ORB, size(ORB), MPI_REAL8, MPI_SUM, COMM_KOREA%mpi_comm, mpierr)
       ETBA%ORB = ORB
     endif
#endif

     return
   endsubroutine

   subroutine set_ldos_atom(PINPT, PGEOM)
     use parameters, only : incar, poscar
     implicit none
     type(incar)            :: PINPT
     type(poscar)           :: PGEOM
     integer*4, allocatable :: proj_atom_dummy(:,:), proj_natom_dummy(:)
     integer*4                 i
     integer*4                 mpierr

     ! setup atom index for projected band if not allocated
     if(PINPT%flag_print_proj_sum) then
       allocate(proj_natom_dummy(PINPT%nproj_sum))
       proj_natom_dummy = 0
       proj_natom_dummy(1:PINPT%nproj_sum) = PINPT%proj_natom(1:PINPT%nproj_sum)
       deallocate(PINPT%proj_natom); allocate(PINPT%proj_natom(PINPT%nproj_sum))
       PINPT%proj_natom = proj_natom_dummy
   
       allocate(proj_atom_dummy(maxval(PINPT%proj_natom(1:PINPT%nproj_sum)),PINPT%nproj_sum))
       proj_atom_dummy = 0
       do i = 1, PINPT%nproj_sum
         proj_atom_dummy(1:PINPT%proj_natom(i),i) = PINPT%proj_atom(1:PINPT%proj_natom(i),i)
       enddo
       deallocate(PINPT%proj_atom);
       allocate(PINPT%proj_atom(maxval(PINPT%proj_natom(1:PINPT%nproj_sum)),PINPT%nproj_sum))
       do i = 1, PINPT%nproj_sum
         PINPT%proj_atom(1:PINPT%proj_natom(i),i) = proj_atom_dummy(1:PINPT%proj_natom(i),i)
       enddo
   
       deallocate(proj_natom_dummy)
       deallocate(proj_atom_dummy)
     endif
        
     return
   endsubroutine

endmodule
