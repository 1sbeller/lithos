program reconstruct_3D_wavefield

  use precision_mod
  use constants_mod
  use mpi_mod

  implicit none

  !*** Begin program (mpi)
  call begin_program




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
