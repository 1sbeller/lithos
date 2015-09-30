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
    
    if (myid == 0) then !!! More data
       irecmin = 1 
       irecmax = nb_rec_by_proc + nb_remain_proc
    else
       irecmin = myid * nb_rec_by_proc + nb_remain_proc +1
       irecmax = (myid + 1) * nb_rec_by_proc + nb_remain_proc
    end if

!    if (myid < nb_remain_proc) then
!       irecmin=myid*(nb_rec_by_proc+1) + 1
!       irecmax=irecmin+nb_rec_by_proc 
!    else
!       irecmin=myid*(nb_rec_by_proc) + 1 + nb_remain_proc              
!       irecmax=irecmin+nb_rec_by_proc-1
!    end if

    do iproc=0,nb_proc-1
       if (iproc == 0) then

          i_inf(iproc+1) = 1
          i_sup(iproc+1) = nb_rec_by_proc + nb_remain_proc

          nb_received(iproc+1)    = nb_rec_by_proc + nb_remain_proc
          nb_received_sv(iproc+1) = nb_received(iproc+1)
        
       else

          i_inf(iproc+1) = iproc * nb_rec_by_proc + nb_remain_proc + 1
          i_sup(iproc+1) = (iproc + 1) * nb_rec_by_proc + nb_remain_proc

          nb_received(iproc+1)    = nb_rec_by_proc 
          nb_received_sv(iproc+1) = nb_received(iproc+1)
     
       end if

!       if (iproc < nb_remain_proc) then
!                    
!          i_inf(iproc+1)=iproc*(nb_rec_by_proc+1) + 1
!          i_sup(iproc+1)=i_inf(iproc+1) + nb_rec_by_proc
!          
!          shift(iproc+1)=iproc*(nb_rec_by_proc+1) + 1 - 1     
!          shift_sv(iproc+1)=shift(iproc+1)*ntime 
!          
!          nb_received(iproc+1)=nb_rec_by_proc+1
!          nb_received_sv(iproc+1)=nb_received(iproc+1) 
!          
!       else
!          
!          i_inf(iproc+1)=iproc*(nb_rec_by_proc) + 1 + nb_remain_proc
!          i_sup(iproc+1)=i_inf(iproc+1) + nb_rec_by_proc - 1
!                    
!          shift(iproc+1)=iproc*(nb_rec_by_proc) + 1 + nb_remain_proc - 1 
!          shift_sv(iproc+1)=shift(iproc+1)*ntime
!          
!          nb_received(iproc+1)=nb_rec_by_proc
!          nb_received_sv(iproc+1)=nb_received(iproc+1) !*ntime       !nb_rec_by_proc*ntime
!          
!       end if
    end do

    nrec_to_store = irecmax-irecmin+1

  end subroutine scatter_data
!--------------------------------------------------------------------------------

!================================================================================
! Broadcast data
  subroutine broadcast_all_data

    if (myid == 0) write(*,*) 'READING OK',9*ntime*nbrec

    call mpi_bcast(ngll,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

    ngllx = ngll
    nglly = ngll
    ngllz = ngll

    !** integers
    call mpi_bcast(ntime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(nbrec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(tbuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(npts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
!    call mpi_bcast(nptsa,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(nsim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ntold,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
!    call mpi_bcast(nlta,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
!    call mpi_bcast(nsta,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(istap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ntnew,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
!    call mpi_bcast(num_bnd_faces,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ngllx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(nglly,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ngllz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(isconv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(istap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ntstf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(ind,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

    if (myid>0) then
       if (allocated(loc2glob)) deallocate(loc2glob)
       if (allocated(xcoord)) deallocate(xcoord)
       if (allocated(ycoord)) deallocate(ycoord)
       if (allocated(zcoord)) deallocate(zcoord)
       if (allocated(abs_bnd_tag)) deallocate(abs_bnd_tag)
       if (allocated(abs_bnd_ielem)) deallocate(abs_bnd_ielem)
       if (allocated(abs_bnd_ijk)) deallocate(abs_bnd_ijk)
       if (allocated(abs_bnd_jacobian2Dw)) deallocate(abs_bnd_jacobian2Dw)
       if (allocated(abs_bnd_normal)) deallocate(abs_bnd_normal)

       if(.not.allocated(loc2glob)) allocate(loc2glob(ngllx,nglly,ngllz,nelem))
       if(.not.allocated(xcoord)) allocate(xcoord(nptsa))
       if(.not.allocated(ycoord)) allocate(ycoord(nptsa))
       if(.not.allocated(zcoord)) allocate(zcoord(nptsa))
       if(.not.allocated(abs_bnd_tag))       allocate(abs_bnd_tag(num_bnd_faces))
       if(.not.allocated(abs_bnd_ielem))       allocate(abs_bnd_ielem(num_bnd_faces))
       if(.not.allocated(abs_bnd_ijk))         allocate(abs_bnd_ijk(3,ngllsquare,num_bnd_faces))
       if(.not.allocated(abs_bnd_jacobian2Dw)) allocate(abs_bnd_jacobian2Dw(ngllsquare,num_bnd_faces))
       if(.not.allocated(abs_bnd_normal))      allocate(abs_bnd_normal(3,ngllsquare,num_bnd_faces))
       if (isconv == 1) then
            if(.not.allocated(stf)) allocate(stf(ntstf))
       end if
    end if
    if (isconv == 1) then
            call MPI_bcast(stf,ntstf,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
    end if
    call MPI_bcast(abs_bnd_tag,num_bnd_faces,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(abs_bnd_ielem,num_bnd_faces,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(abs_bnd_ijk,3*ngllsquare*num_bnd_faces,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(abs_bnd_normal,3*ngllsquare*num_bnd_faces,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(abs_bnd_jacobian2Dw,ngllsquare*num_bnd_faces,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)

    call MPI_bcast(xcoord,nptsa,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(ycoord,nptsa,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi) 
    call MPI_bcast(zcoord,nptsa,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
    call MPI_bcast(loc2glob,ngllx*nglly*ngllz*nelem,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

    !** Reals
    call mpi_bcast(dtold,1,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(dtnew,1,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
!    call mpi_bcast(thres,1,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(alpha,1,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(fmax,1,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)


    !*** Characters
    call MPI_bcast(fwdtool,3,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)


    if (myid > 0) then
       if (allocated(working_axisem_dir)) deallocate(working_axisem_dir)
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
    
    if (allocated(nb_received)) deallocate(nb_received)
    if (allocated(shift)) deallocate(shift)
    if(allocated(nb_received_sv)) deallocate(nb_received_sv)
    if (allocated(shift_sv)) deallocate(shift_sv)
    if (allocated(i_inf)) deallocate(i_inf)
    if (allocated(i_sup)) deallocate(i_sup)
    allocate(nb_received(nb_proc),shift(nb_proc),nb_received_sv(nb_proc),shift_sv(nb_proc))
    allocate(i_inf(nb_proc),i_sup(nb_proc))
    
  end subroutine broadcast_all_data
!--------------------------------------------------------------------------------




end module mpi_int_mod


