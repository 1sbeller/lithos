program reconstruct_3D_wavefield

  use precision_mod
  use constants_mod

  use mpi_mod
  use inputs_outputs_mod, only: read_all_inputs, read_input_box_file
  use interp_mesh_mod,    only: determine_connectivity
  use post_process_wavefields_mod, only: reconstruct_velocity, reconstruct_stress
  use global_parameters_mod, only:ntime, energy

  implicit none

  integer(kind=si) :: isim, ipart

  !*** Begin program (mpi)
  call begin_program

  !*** Prepare reconstruction
  isim=1
  if (myid == 0) then
     call read_all_inputs(isim,nsim)  ! Rotations are defined here
  end if
  call MPI_bcast(npart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
  call MPI_bcast(tbuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

!  do ipart = 1, npart

     !*** Read each input_box file
     if (myid == 0) call read_input_box_file(isim) !,ipart)
     call MPI_bcast(nbrec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

     !*** Avoid empty subdomains
!     if (nbrec < 1) then
!        cycle
!     end if   

     !*** Find points for current subdomain
     if (myid == 0) call determine_connectivity

     !*** Distribute on MPI process
     call alloc_all_mpi
     call bcast_all_mpi
     call distribute_mpi
     
     !*** Reconstruct wavefields
     call reconstruct_velocity() !ipart)  !(isim)
     call reconstruct_stress()   !ipart)    !(isim)

!  end do

  
  !*** Save energy for stalta and infos for subdomains
  if (myid == 0) then
     write(6,*)'Reconstruction done !'
     open(17,file='energy_for_stalta.bin',status='replace',access='direct',recl=cp*ntime)
     write(17,rec=1)energy
     close(17)
     open(17,file='infos_split2subdom',status='replace',action='write')
     do ipart=1,npart
        write(17,'(i8.8,i8.8,i8.8)')tab_box_rec(1:3,ipart)
     end do
     close(17)
  end if

  !*** End the program (mpi)
  call finish_program

contains

  subroutine begin_program

    call init_mpi

  end subroutine begin_program

  subroutine finish_program
    
    call finalize_mpi

  end subroutine finish_program

end program reconstruct_3D_wavefield
