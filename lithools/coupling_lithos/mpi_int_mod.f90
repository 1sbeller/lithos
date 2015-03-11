module mpi_int_mod

  use precision_mod
  use constants_mod
  use mpi
  use interpolation_parameters_mod

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

  end subroutine init_mpi
!--------------------------------------------------------------------------------

!================================================================================
! Finalize mpi
  subroutine finalize_mpi
    
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
    call MPI_finalize(ierr_mpi)

  end subroutine finalize_mpi
!--------------------------------------------------------------------------------

!================================================================================
! Determine recmin recmax
  subroutine scatter_data

    nb_rec_by_proc=nbrec/nb_proc
    nb_remain_proc=mod(nbrec,nb_proc)
    write(*,*) ' decomposition ', nb_rec_by_proc,nb_remain_proc
    if (myid < nb_remain_proc) then
       irecmin=myid*(nb_rec_by_proc+1) + 1
       irecmax=irecmin+nb_rec_by_proc 
    else
       irecmin=myid*(nb_rec_by_proc) + 1 + nb_remain_proc              
       irecmax=irecmin+nb_rec_by_proc-1
    end if

    do iproc=0,nb_proc-1
       if (iproc < nb_remain_proc) then
                    
          i_inf(iproc+1)=iproc*(nb_rec_by_proc+1) + 1
          i_sup(iproc+1)=i_inf(iproc+1) + nb_rec_by_proc
          
          shift(iproc+1)=iproc*(nb_rec_by_proc+1) + 1 - 1     
          shift_sv(iproc+1)=shift(iproc+1)*ntime 
          
          nb_received(iproc+1)=nb_rec_by_proc+1
          nb_received_sv(iproc+1)=nb_received(iproc+1) 
          
       else
          
          i_inf(iproc+1)=iproc*(nb_rec_by_proc) + 1 + nb_remain_proc
          i_sup(iproc+1)=i_inf(iproc+1) + nb_rec_by_proc - 1
                    
          shift(iproc+1)=iproc*(nb_rec_by_proc) + 1 + nb_remain_proc - 1 
          shift_sv(iproc+1)=shift(iproc+1)*ntime
          
          nb_received(iproc+1)=nb_rec_by_proc
          nb_received_sv(iproc+1)=nb_received(iproc+1) !*ntime       !nb_rec_by_proc*ntime
          
       end if
    end do

    nrec_to_store = irecmax-irecmin+1

  end subroutine scatter_data
!--------------------------------------------------------------------------------

!================================================================================
! Broadcast data
  subroutine broadcast_all_data

    if (myid == 0) write(*,*) 'READING OK',9*ntime*nbrec

    !** integers
    call mpi_bcast(ntime,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(nbrec,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(npts,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(nsim,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ntold,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(nlta,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(nsta,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(istap,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ntnew,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(num_bnd_faces,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ngllsquare,1,MPISI,0,MPI_COMM_WORLD,ierr_mpi)

    call MPI_bcast(abs_bnd_ielem,num_bnd_faces,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(abs_bnd_ijk,3*ngllsquare*num_bnd_faces,MPISI,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(abs_bnd_normal,3*ngllsquare*num_bnd_faces,MPICP,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(abs_bnd_jacobian2Dw,ngllsquare*num_bnd_faces,MPICP,0,MPI_COMM_WORLD,ierr_mpi)

    call MPI_bcast(xcoord,npts,MPICP,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(ycoord,npts,MPICP,0,MPI_COMM_WORLD,ierr_mpi) 
    call MPI_bcast(zcoord,npts,MPICP,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(loc2glob,ngllx*nglly*ngllz*nelem,MPISI,0,MPI_COMM_WORLD,ierr_mpi)

    !** Reals
    call mpi_bcast(dtold,1,MPICP,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(dtnew,1,MPICP,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(thres,1,MPICP,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(alpha,1,MPICP,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(fmax,1,MPICP,0,MPI_COMM_WORLD,ierr_mpi)


    !*** Characters
    call MPI_bcast(fwdtool,3,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)


    if (myid > 0) then
       allocate(working_axisem_dir(nsim))
       
       if (nsim == 1) then
          working_axisem_dir(1)="./"
       else
          working_axisem_dir(1) = "MZZ/"
          working_axisem_dir(2) = "MXX_P_MYY/"
          working_axisem_dir(3) = "MXZ_MYZ/"
          working_axisem_dir(4) = "MXY_MXX_M_MYY/"
       end if
    end if
    
    allocate(nb_received(nb_proc),shift(nb_proc),nb_received_sv(nb_proc),shift_sv(nb_proc))
    allocate(i_inf(nb_proc),i_sup(nb_proc))
    
  end subroutine broadcast_all_data
!--------------------------------------------------------------------------------




end module mpi_int_mod










end module mpi_int_mod
