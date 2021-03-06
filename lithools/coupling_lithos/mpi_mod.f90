module mpi_mod
  
  use mpi
  use precision_mod
  use global_parameters_mod
  use rotation_matrix_mod

  !*** Mpi variables
  integer(kind=si) :: myid, nb_proc, ierr_mpi

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
    call MPI_comm_size(MPI_COMM_WORLD, nb_proc, ierr_mpi)
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)

    !*** Define mpi kinds
    call MPI_Type_create_f90_integer(range(dummydi), mpidi, ierr_mpi)
    call MPI_Type_create_f90_integer(range(dummysi), mpisi, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummysp), range(dummysp), mpisp, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummydp), range(dummydp), mpidp, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummycp), range(dummycp), mpicp, ierr_mpi)
    call MPI_Type_create_f90_real(precision(dummyhp), range(dummyhp), mpihp, ierr_mpi)    

    call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)

    !*** Define operation
!    call MPI_Op_create(mysum, .true., mympisum, ierr_mpi)   
!    call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)

  end subroutine init_mpi
!--------------------------------------------------------------------------------

!================================================================================
! Finalize mpi
  subroutine finalize_mpi
    
!    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
    call MPI_finalize(ierr_mpi)

  end subroutine finalize_mpi
!--------------------------------------------------------------------------------


!================================================================================
! Allocate mpi arrays
  subroutine alloc_all_mpi

   call mpi_bcast(nsim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(ntime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(nbrec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(nbproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi) 
   call mpi_bcast(ibeg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(nel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(rot_mat,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(trans_rot_mat,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(rot_mat_mesh,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
   call mpi_bcast(trans_rot_mat_mesh,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)

   if(allocated(Mij_scale)) deallocate(Mij_scale)
   allocate(Mij_scale(6))
   Mij_scale(:) = 0.

   if (myid >  0) then 

      if (allocated(tab_box_rec)) deallocate(tab_box_rec)
      allocate(tab_box_rec(1:3,1:npart))

      if (allocated(reciever_geogr)) deallocate(reciever_geogr)
      if (allocated(reciever_sph)) deallocate(reciever_sph)
      if (allocated(reciever_cyl)) deallocate(reciever_cyl)
      if (allocated(reciever_interp_value)) deallocate(reciever_interp_value)
      allocate(reciever_geogr(3,nbrec),reciever_sph(3,nbrec),reciever_cyl(3,nbrec),reciever_interp_value(nbrec))

      if (allocated(data_rec)) deallocate(data_rec)
      if (allocated(stress_rec)) deallocate(stress_rec)
      if (allocated(stress_to_write)) deallocate(stress_to_write)
      if (allocated(strain_rec)) deallocate(strain_rec)
      allocate(data_rec(nbrec,3),stress_rec(nbrec,6),stress_to_write(nbrec,6),strain_rec(nbrec,6))

      if (allocated(data_rec_all)) deallocate(data_rec_all)
      if (allocated(stress_rec_all)) deallocate(stress_rec_all)
      allocate(data_rec_all(nbrec,3),stress_rec_all(nbrec,6))

      if (allocated(f1)) deallocate(f1)
      if (allocated(f2)) deallocate(f2)
      allocate(f1(nsim,nbrec),f2(nsim,nbrec))

      if (allocated(phi)) deallocate(phi)
      allocate(phi(nbrec))
      
      if (allocated(magnitude)) deallocate(magnitude)
      allocate(magnitude(nsim))
      
      if (allocated(Mij)) deallocate(Mij)
      allocate(Mij(nsim,6))

      if (allocated(scoor)) deallocate(scoor)
      if (allocated(zcoor)) deallocate(zcoor)
      allocate(scoor(ibeg:iend,ibeg:iend,nel),zcoor(ibeg:iend,ibeg:iend,nel))

      if (allocated(data_read)) deallocate(data_read)
      allocate(data_read(ibeg:iend,ibeg:iend,nel))

      if (allocated(xi_rec)) deallocate(xi_rec)
      if (allocated(eta_rec)) deallocate(eta_rec)
      allocate(xi_rec(nbrec),eta_rec(nbrec))

      if (allocated(rec2elm)) deallocate(rec2elm)
      allocate(rec2elm(nbrec))

      if (allocated(src_type)) deallocate(src_type)
      allocate(src_type(nsim,2))

   end if

   if (allocated(stress_reduce)) deallocate(stress_reduce)
   if (allocated(data_reduce)) deallocate(data_reduce)
   allocate(stress_reduce(nbrec,6),data_reduce(nbrec,3))
   call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)


  end subroutine alloc_all_mpi
!--------------------------------------------------------------------------------

!================================================================================
! Broacdact arrays
  subroutine bcast_all_mpi
    
    ! single
    call mpi_bcast(data_rec,3*nbrec,MPI_REAL4,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(data_rec_all,3*nbrec,MPI_REAL4,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(stress_to_write,6*nbrec,MPI_REAL4,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(stress_rec,6*nbrec,MPI_REAL4,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(stress_rec_all,6*nbrec,MPI_REAL4,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(strain_rec,6*nbrec,MPI_REAL4,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(f1,nbrec*nsim,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(f2,nbrec*nsim,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(Mij,nsim*6,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(Mij_scale,6,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(magnitude,nsim,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)

    call mpi_bcast(phi,nbrec,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)

    ! double 
    call mpi_bcast(scoor,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(zcoor,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(reciever_geogr,3*nbrec,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)    
    call mpi_bcast(reciever_sph,3*nbrec,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(reciever_cyl,3*nbrec,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(xi_rec,nbrec,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(eta_rec,nbrec,MPI_REAL8,0,MPI_COMM_WORLD,ierr_mpi)

    ! integer
    call mpi_bcast(rec2elm,nbrec,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(tab_box_rec,3*npart,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

    ! character 
     call mpi_bcast(src_type,10*2*nsim,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)
     call mpi_bcast(coup_tool,3,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)

  end subroutine bcast_all_mpi
!--------------------------------------------------------------------------------

!================================================================================
! Distribute array on procs
  subroutine distribute_mpi

    integer(kind=si) :: nb_by_pc, left   
    
    nb_by_pc=(nbrec/nb_proc)
    left=nbrec-nb_by_pc*nb_proc
    
    if (myid < left) then
      irecmin = (myid*(nb_by_pc+1)) + 1 
      irecmax = irecmin + nb_by_pc
    else
      irecmin = myid*(nb_by_pc) + left + 1
      irecmax = irecmin + nb_by_pc - 1   
    end if
 
  end subroutine distribute_mpi
!--------------------------------------------------------------------------------

!================================================================================
! Reduce velocity
  subroutine reduce_mpi_veloc

    call mpi_reduce(data_rec_all,data_reduce,nbrec*3,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
    data_rec_all(:,:)=data_reduce(:,:)

  end subroutine reduce_mpi_veloc
!--------------------------------------------------------------------------------

!================================================================================
! Reduce stress
  subroutine reduce_mpi_stress
    
    call mpi_reduce(stress_rec_all,stress_reduce,nbrec*6,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
    stress_rec_all(:,:)=stress_reduce(:,:)
  
  end subroutine reduce_mpi_stress
!--------------------------------------------------------------------------------

end module mpi_mod
