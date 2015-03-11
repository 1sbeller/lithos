program reconstruct_3D_wavefield

  use precision_mod
  use constants_mod

  use mpi_mod
  use inputs_outputs_mod, only: read_all_inputs
  use interp_mesh_mod,    only: determine_connectivity

  implicit none

  integer(kind=si) :: isim

  !*** Begin program (mpi)
  call begin_program

  !*** Prepare reconstruction
  if (myid == 0) then
     call read_all_inputs(isim,nsim)  ! Rotations are defined here
     call determine_connectivity
  end if

  !*** Distribute on MPI process
  call distribute_mpi

  !*** Reconstruct wavefields
  do isim=1,nsim
     call reconstruct_velocity(isim)
     call reconstruct_stress(isim)
  end do

  !*** End the program (mpi)
  call finish_program

contains

  subroutine begin_program

    call init_mpi

  end subroutine begin_program

  subroutine finish_program
    
    call finalize_mpi

  end subroutine finish_program

end program
