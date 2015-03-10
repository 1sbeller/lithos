module mpi_mod
  
  use mpi
  use precision_mod

  !*** Mpi variables
  integer(kind=si) :: myid, nbproc, ierr_mpi

  !*** MPI types
  integer :: mpidi, mpisi, mpisp, mpidp, mpicp, mpihp

contains

!================================================================================
! Inititalize mpi
  subroutine init_mpi

    integer(kind=di) :: dummydi
    integer(kind=si) :: dummysi
    real(kind=sp) :: dummysp
    real(kind=dp) :: dummydp
    real(kind=cp) :: dummycp
    real(kind=hp) :: dummyhp
    
    !*** Get infos from communicator
    call MPI_init(ierr_mpi)
    call MPI_comm_rank(MPI_COMM_WORLD, myid,   ierr_mpi)
    call MPI_comm_size(MPI_COMM_WORLD, nbproc, ierr_mpi)
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)

    !*** Define mpi kinds
    call MPI_Type_create_f90_integer(range(dummydi), mpidi, ierr_mpi)
    call MPI_Type_create_f90_integer(range(dummysi), mpisi, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummysp), range(dummysp), mpisp, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummydp), range(dummydp), mpidp, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummycp), range(dummycp), mpicp, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummyhp), range(dummyhp), mpihp, ierr_mpi)    

  end subroutine init_mpi
!--------------------------------------------------------------------------------

!================================================================================
! Finalize mpi
  subroutine finalize_mpi
    
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
    call MPI_finalize(ierr_mpi)

  end subroutine finalize_mpi
!--------------------------------------------------------------------------------








end module mpi_mod
