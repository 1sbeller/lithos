program interpolate_3D_wavefield

  use precision_mod
  use constants_mod
  use mpi_int_mod
  use interpolation_parameters_mod
  use interp_process_mod
  use inputs_outputs_interp_mod

  implicit none

  integer(kind=si) :: warning=0, ipart, ind2, lit
  character(len=80) :: ficii
  integer(kind=si), parameter :: IIN=11
  character(len=5) :: myfileend 
  character(len=9)  :: myfileend1
  character(len=80) :: srcfile


  !================================================================================
  ! Prepare reconstruction
  !--------------------------------------------------
  !*** Begin program and mpi
  call begin_program
  
  !*** Read inut files parameters
  call read_all_inputs
  call mpi_bcast(npart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

  !*** Read infos on partitionning
  if(.not.allocated(tab_box_rec)) allocate(tab_box_rec(3,npart))
  if (.not.allocated(part_info)) allocate(part_info(npart,6))
  if (myid == 0) then
     open(155,file='partition_info.txt',status='old',action='read')
     do ipart = 1,npart
        read(155,'(6i9.8)')part_info(ipart,1:6)  !npart, ipart, myne, myndof, myntag, myntag*25
     end do
     close(155)
     open(155,file='infos_split2subdom',status='old',action='read')
     do ipart = 1,npart
        read(155,'(3i8.8)')tab_box_rec(1:3,ipart)  !nbrecloc,istart,iend
     end do
     close(155)
  end if
  print *,'done 1'

  call MPI_bcast(part_info,6*npart,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
  call MPI_bcast(tab_box_rec,3*npart,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
  do ipart = 1, npart

     !*** Get infos
     nelem         = part_info(ipart,3)
     nptsa         = part_info(ipart,4)
     num_bnd_faces = part_info(ipart,5)
     
     if (num_bnd_faces < 1) then 
        cycle
     end if
     call mpi_bcast(ngllsquare,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
     call mpi_bcast(ngll,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

     !*** Read arrays
     write(myfileend,'(a,i4.4)')'_',ipart

     if (myid == 0) then
        if (allocated(loc2glob)) deallocate(loc2glob)
        if (allocated(xcoord)) deallocate(xcoord)
        if (allocated(ycoord)) deallocate(ycoord)
        if (allocated(zcoord)) deallocate(zcoord)
        if (allocated(abs_bnd_tag)) deallocate(abs_bnd_tag)
        if (allocated(abs_bnd_ielem)) deallocate(abs_bnd_ielem)
        if (allocated(abs_bnd_ijk)) deallocate(abs_bnd_ijk)
        if (allocated(abs_bnd_jacobian2Dw)) deallocate(abs_bnd_jacobian2Dw)
        if (allocated(abs_bnd_normal)) deallocate(abs_bnd_normal)

        if(.not.allocated(loc2glob)) allocate(loc2glob(ngll,ngll,ngll,nelem))
        if(.not.allocated(xcoord)) allocate(xcoord(nptsa))
        if(.not.allocated(ycoord)) allocate(ycoord(nptsa))
        if(.not.allocated(zcoord)) allocate(zcoord(nptsa))
        if(.not.allocated(abs_bnd_tag))       allocate(abs_bnd_tag(num_bnd_faces))
        if(.not.allocated(abs_bnd_ielem))       allocate(abs_bnd_ielem(num_bnd_faces))
        if(.not.allocated(abs_bnd_ijk))         allocate(abs_bnd_ijk(3,ngllsquare,num_bnd_faces))
        if(.not.allocated(abs_bnd_jacobian2Dw)) allocate(abs_bnd_jacobian2Dw(ngllsquare,num_bnd_faces))
        if(.not.allocated(abs_bnd_normal))      allocate(abs_bnd_normal(3,ngllsquare,num_bnd_faces))
        
        !*** X-coord
        open(unit=IIN,file=rep(1:len_trim(rep))//'/MESH/'//'coordinates'//myfileend//'.bin',status='old',form='unformatted')
        read(IIN) xcoord(:)
        read(IIN) ycoord(:)
        read(IIN) zcoord(:)
        close(IIN)
        
        !*** Global indexing
        open(unit=IIN,file=rep(1:len_trim(rep))//'/MESH/'//'loc2glob'//myfileend//'.bin',status='old',form='unformatted')
        read(IIN) loc2glob(:,:,:,:)
        close(IIN)
        
        !*** Boundary conditions
        open(unit=IIN,file=rep(1:len_trim(rep))//'/MESH/'//'boundary'//myfileend//'.bin',status='old',form='unformatted')
        read(IIN) num_bnd_faces
        read(IIN) abs_bnd_tag
        read(IIN) abs_bnd_ielem
        read(IIN) abs_bnd_ijk
        read(IIN) abs_bnd_normal
        read(IIN) abs_bnd_jacobian2Dw
        close(IIN)
     end if

     !*** Inputs files
     if (ipart ==1) then
        if(myid == 0)  call define_axisem_dir !(ipart)
     end if

     nbrec = maxval(tab_box_rec(3,:))
     npts = tab_box_rec(1,ipart)
     call broadcast_all_data

     !*** Prepare reading of reconstructed AxiSEM outputs
     if (ipart == 1) then
        ntold = ntime   
        if (myid==0) write(6,*)'Must read ',nbrec,' points for ',ntime,' time steps.'
        if (myid==0) write(6,*)'Check ntold and ntime : ',ntold,ntime
        if (itend > ntold) then
           write(6,*)'WARNING itend > olden will STOP'
           stop 'itend > ntold'
        end if
     end if
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     
!!$     if (myid == 0) then
!!$        if(.not.allocated(data_tmp)) allocate(data_tmp(nbrec,9))
!!$        if (.not.allocated(vxold2)) allocate(vxold2(nbrec,oldlen))
!!$        if (.not.allocated(vyold2)) allocate(vyold2(nbrec,oldlen))
!!$        if (.not.allocated(vzold2)) allocate(vzold2(nbrec,oldlen))
!!$        if (.not.allocated(sxxold2)) allocate(sxxold2(nbrec,oldlen))
!!$        if (.not.allocated(syyold2)) allocate(syyold2(nbrec,oldlen))
!!$        if (.not.allocated(szzold2)) allocate(szzold2(nbrec,oldlen))
!!$        if (.not.allocated(syzold2)) allocate(syzold2(nbrec,oldlen))
!!$        if (.not.allocated(sxzold2)) allocate(sxzold2(nbrec,oldlen))
!!$        if (.not.allocated(sxyold2)) allocate(sxyold2(nbrec,oldlen))
!!$        if (.not.allocated(vxold2)) allocate(vxold2(npts,1))
!!$        if (.not.allocated(vyold2)) allocate(vyold2(npts,1))
!!$        if (.not.allocated(vzold2)) allocate(vzold2(npts,1))
!!$        if (.not.allocated(sxxold2)) allocate(sxxold2(npts,1))
!!$        if (.not.allocated(syyold2)) allocate(syyold2(npts,1))
!!$        if (.not.allocated(szzold2)) allocate(szzold2(npts,1))
!!$        if (.not.allocated(syzold2)) allocate(syzold2(npts,1))
!!$        if (.not.allocated(sxzold2)) allocate(sxzold2(npts,1))
!!$        if (.not.allocated(sxyold2)) allocate(sxyold2(npts,1))
!!$        
!!$        vxold2(:,:) = 0.
!!$        vyold2(:,:) = 0.
!!$        vzold2(:,:) = 0.
!!$        sxxold2(:,:) = 0.
!!$        syyold2(:,:) = 0.
!!$        szzold2(:,:) = 0.
!!$        syzold2(:,:) = 0.
!!$        sxzold2(:,:) = 0.
!!$        sxyold2(:,:) = 0.
!!$     end if
     
     
     write(6,*)'veriflala',nbrec,npts

     call scatter_data


     itbeg  = ind !-ceiling(2.5/fmax)
     itend  = itbeg + ceiling(ntnew * dtnew / dtold)
     tbeg = itbeg * dtold      !*** Starting time of cut signal
     tend = itend * dtold
     oldlen = itend-itbeg+1

     
     if (myid==0) write(6,*)'Infos : itbeg, itend, tbeg, tend, oldlen, ntold, dtold'
     if (myid==0) write(6,*)itbeg, itend, tbeg, tend, oldlen, ntold, dtold

     
     !print *,'DEBUG npts,nptsa,nrec_to_store : ',npts,nptsa,nrec_to_store
     if (allocated(vxold1)) deallocate(vxold1)
     if (allocated(vyold1)) deallocate(vyold1)
     if (allocated(vzold1)) deallocate(vzold1)
     if (allocated(sxxold1)) deallocate(sxxold1)
     if (allocated(syyold1)) deallocate(syyold1)
     if (allocated(szzold1)) deallocate(szzold1)
     if (allocated(sxyold1)) deallocate(sxyold1)
     if (allocated(syzold1)) deallocate(syzold1)
     if (allocated(sxzold1)) deallocate(sxzold1)
     if (.not.allocated(vxold1)) allocate(vxold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(vyold1)) allocate(vyold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(vzold1)) allocate(vzold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(sxxold1)) allocate(sxxold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(syyold1)) allocate(syyold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(szzold1)) allocate(szzold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(syzold1)) allocate(syzold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(sxzold1)) allocate(sxzold1(nrec_to_store,oldlen)) !ntold))
     if (.not.allocated(sxyold1)) allocate(sxyold1(nrec_to_store,oldlen)) !ntold))
     
     vxold1(:,:) = 0.
     vyold1(:,:) = 0.
     vzold1(:,:) = 0.
     sxxold1(:,:) = 0.
     syyold1(:,:) = 0.
     szzold1(:,:) = 0.
     syzold1(:,:) = 0.
     sxzold1(:,:) = 0.
     sxyold1(:,:) = 0.
    
 
     !================================================================================
     ! Read AxiSEM reconstructed files
     !--------------------------------------------------
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)

     if (ipart == 1) then  ! Read only some timestep of all points by master procs
        if (myid == 0) then
           if(.not.allocated(data_tmp)) allocate(data_tmp(nbrec,9))
           if (.not.allocated(vxold2)) allocate(vxold2(nbrec,oldlen))
           if (.not.allocated(vyold2)) allocate(vyold2(nbrec,oldlen))
           if (.not.allocated(vzold2)) allocate(vzold2(nbrec,oldlen))
           if (.not.allocated(sxxold2)) allocate(sxxold2(nbrec,oldlen))
           if (.not.allocated(syyold2)) allocate(syyold2(nbrec,oldlen))
           if (.not.allocated(szzold2)) allocate(szzold2(nbrec,oldlen))
           if (.not.allocated(syzold2)) allocate(syzold2(nbrec,oldlen))
           if (.not.allocated(sxzold2)) allocate(sxzold2(nbrec,oldlen))
           if (.not.allocated(sxyold2)) allocate(sxyold2(nbrec,oldlen))
!!$        if (.not.allocated(vxold2)) allocate(vxold2(npts,1))
!!$        if (.not.allocated(vyold2)) allocate(vyold2(npts,1))
!!$        if (.not.allocated(vzold2)) allocate(vzold2(npts,1))
!!$        if (.not.allocated(sxxold2)) allocate(sxxold2(npts,1))
!!$        if (.not.allocated(syyold2)) allocate(syyold2(npts,1))
!!$        if (.not.allocated(szzold2)) allocate(szzold2(npts,1))
!!$        if (.not.allocated(syzold2)) allocate(syzold2(npts,1))
!!$        if (.not.allocated(sxzold2)) allocate(sxzold2(npts,1))
!!$        if (.not.allocated(sxyold2)) allocate(sxyold2(npts,1))
        
           vxold2(:,:) = 0.
           vyold2(:,:) = 0.
           vzold2(:,:) = 0.
           sxxold2(:,:) = 0.
           syyold2(:,:) = 0.
           szzold2(:,:) = 0.
           syzold2(:,:) = 0.
           sxzold2(:,:) = 0.
           sxyold2(:,:) = 0.
        end if
        
        if (myid==0) write(6,*)'Read AxisEM files...'
        lit = 0
        
        do itime=1,ntime
           
           if (myid == 0) then ! Read on master proc
                      
              if (modulo(itime,tbuff)==0) then
                 write(6,*)'time ',100.*itime/ntime,'%'
              end if
              
              do isim=1,nsim
                 
                 read(ivx(isim))  data_tmp(:,1)
                 read(ivy(isim))  data_tmp(:,2)
                 read(ivz(isim))  data_tmp(:,3)
                 read(isxx(isim)) data_tmp(:,4)
                 read(isyy(isim)) data_tmp(:,5)
                 read(iszz(isim)) data_tmp(:,6)
                 read(isxy(isim)) data_tmp(:,7)
                 read(isxz(isim)) data_tmp(:,8)
                 read(isyz(isim)) data_tmp(:,9)
                 
                 !              vxold2(:,1)  = data_tmp(:,1) ! vxold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),1))
                 !              vyold2(:,1)  = data_tmp(:,2) !vyold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),2))
                 !              vzold2(:,1)  = data_tmp(:,3) !vzold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),3))
                 !              sxxold2(:,1) = data_tmp(:,4) !sxxold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),4))
                 !              syyold2(:,1) = data_tmp(:,5) !syyold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),5))
                 !              szzold2(:,1) = data_tmp(:,6) !szzold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),6))
                 !              sxyold2(:,1) = data_tmp(:,7) !sxyold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),7))
                 !              sxzold2(:,1) = data_tmp(:,8) !sxzold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),8))
                 !              syzold2(:,1) = data_tmp(:,9) !syzold2(:,1) + real(data_tmp(tab_box_rec(2,ipart):tab_box_rec(3,ipart),9))
                 
              end do
              
              if (itime >= itbeg .and. itime <= itend) then
                 
                 lit = lit + 1
                 
                 vxold2(:,lit)  =  data_tmp(:,1) ! vxold2(i_inf(1):i_sup(1),1)
                 vyold2(:,lit)  =  data_tmp(:,2) !vyold2(i_inf(1):i_sup(1),1)
                 vzold2(:,lit)  =  data_tmp(:,3) !vzold2(i_inf(1):i_sup(1),1)
                 sxxold2(:,lit) =  data_tmp(:,4) !sxxold2(i_inf(1):i_sup(1),1)
                 syyold2(:,lit) =  data_tmp(:,5) !syyold2(i_inf(1):i_sup(1),1)
                 szzold2(:,lit) =  data_tmp(:,6) !szzold2(i_inf(1):i_sup(1),1)
                 sxyold2(:,lit) =  data_tmp(:,7) !sxyold2(i_inf(1):i_sup(1),1)
                 sxzold2(:,lit) =  data_tmp(:,8) !sxzold2(i_inf(1):i_sup(1),1)
                 syzold2(:,lit) =  data_tmp(:,9) !syzold2(i_inf(1):i_sup(1),1)
                 
              end if
              
           end if
        end do

     end if
     
     i_inf = i_inf + tab_box_rec(2,ipart) -1
     i_sup = i_sup + tab_box_rec(2,ipart) -1

     !*** Then send/receiv
     !*** Send-receive to scatter data
!     do itime = 1,oldlen
        if (myid ==0) then ! SEND
           
           vxold1(:,:)  =  vxold2(i_inf(1):i_sup(1),:)
           vyold1(:,:)  =  vyold2(i_inf(1):i_sup(1),:)
           vzold1(:,:)  =  vzold2(i_inf(1):i_sup(1),:)
           sxxold1(:,:) =  sxxold2(i_inf(1):i_sup(1),:)
           syyold1(:,:) =  syyold2(i_inf(1):i_sup(1),:)
           szzold1(:,:) =  szzold2(i_inf(1):i_sup(1),:)
           sxyold1(:,:) =  sxyold2(i_inf(1):i_sup(1),:)
           sxzold1(:,:) =  sxzold2(i_inf(1):i_sup(1),:)
           syzold1(:,:) =  syzold2(i_inf(1):i_sup(1),:)
           
           do iproc=1,nb_proc-1
              
              call mpi_send( vxold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq1,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send( vyold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq2,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send( vzold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq3,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send(sxxold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq4,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send(syyold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq5,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send(szzold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq6,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send(sxyold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq7,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send(sxzold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq8,MPI_COMM_WORLD,ierr_mpi)
              call mpi_send(syzold2(i_inf(iproc+1):i_sup(iproc+1),:),nb_received_sv(iproc+1)*oldlen,MPI_REAL,iproc,etq9,MPI_COMM_WORLD,ierr_mpi)
              
           
!!$                 
!!$                 call mpi_send( vxold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq1,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send( vyold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq2,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send( vzold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq3,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send(sxxold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq4,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send(syyold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq5,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send(szzold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq6,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send(sxyold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq7,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send(sxzold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq8,MPI_COMM_WORLD,ierr_mpi)
!!$                 call mpi_send(syzold2(i_inf(iproc+1):i_sup(iproc+1),1),nb_received_sv(iproc+1),MPI_REAL,iproc,etq9,MPI_COMM_WORLD,ierr_mpi)
                 
           end do
           
        else ! RECEIVE
           
           call mpi_recv( vxold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq1,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv( vyold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq2,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv( vzold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq3,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv(sxxold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq4,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv(syyold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq5,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv(szzold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq6,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv(sxyold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq7,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv(sxzold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq8,MPI_COMM_WORLD,statut,ierr_mpi)
           call mpi_recv(syzold1(:,:),nb_received_sv(myid+1)*oldlen,MPI_REAL,0,etq9,MPI_COMM_WORLD,statut,ierr_mpi)
           
        end if

 !    end do
     
     if (ipart == 1) then
        if (myid==0) then
           do isim=1,nsim
              close(ivx(isim))
              close(ivy(isim))
              close(ivz(isim))
              close(isxx(isim))
              close(isyy(isim))
              close(iszz(isim))
              close(isxy(isim))
              close(isxz(isim))
              close(isyz(isim))
           end do
           deallocate(data_tmp)
           write(6,*)'Done'
        end if
     end if
!!$     if (myid ==0) then
!!$        deallocate(vxold2)
!!$        deallocate(vyold2)
!!$        deallocate(vzold2)
!!$        deallocate(sxxold2)
!!$        deallocate(syyold2)
!!$        deallocate(szzold2)
!!$        deallocate(syzold2)
!!$        deallocate(sxzold2)
!!$        deallocate(sxyold2)
!!$     end if

     ntime = oldlen
     ntold = oldlen

!     write(ficii,'(a,i4.4,a)')'sig_',myid,'.dat' 
!     open(20000+myid,file=trim(ficii),action='write')
!     write(20000+myid,*)vzold1(1,:)
!     close(20000+myid)


     !================================================================================
     ! Compute velocity energy and STA/LTA this could be done on many procs)
     !--------------------------------------------------
!!$     !*** Compute energy
!!$     if (myid==0) write(6,*)'Compute STA/LTA...'
!!$     if (.not.allocated(vpow)) allocate(vpow(ntold))
!!$     if (.not.allocated(stalta)) allocate(stalta(ntold))
!!$     vpow = 0.
!!$     stalta = 0.
!!$
!!$     do itold = 1, ntold
!!$        sumx = 0.
!!$        sumy = 0.
!!$        sumz = 0.
!!$        do ipt =  1, nrec_to_store  !irecmin, irecmax
!!$           sumx = sumx + vxold1(ipt,itold)**2
!!$           sumy = sumy + vyold1(ipt,itold)**2
!!$           sumz = sumz + vzold1(ipt,itold)**2
!!$        end do
!!$        call MPI_allreduce(MPI_IN_PLACE,sumx,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
!!$        call MPI_allreduce(MPI_IN_PLACE,sumy,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
!!$        call MPI_allreduce(MPI_IN_PLACE,sumz,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
!!$        vpow(itold) = sqrt(sumx + sumy +sumz)
!!$     end do
!!$
!!$     if (myid == 0) then
!!$        open(26,file='verif.bin',access='direct',recl=cp*ntold)
!!$        write(26,rec=1)vpow
!!$        close(26)
!!$
!!$        if (isconv ==1) then         
!!$           open(26,file='stfint.bin',access='direct',recl=cp*ntold)
!!$           write(26,rec=1)stf
!!$           close(26)
!!$        end if
!!$     end if
!!$
!!$
!!$     !*** Compute STA/LTA
!!$     if (isconv == 1) then
!!$        if(.not.allocated(conv)) allocate(conv(ntold))
!!$        call myconvolution(vpow,stf,ntold,ntstf,conv)
!!$        vpow = conv
!!$        if (allocated(conv)) deallocate(conv)
!!$     end if
!!$     if (myid ==0) then
!!$        open(26,file='verifconv.bin',access='direct',recl=cp*ntold)
!!$        write(26,rec=1)vpow
!!$        close(26)
!!$     end if
!!$     call substalta(vpow, ntold,  nsta, nlta, thres, stalta, ind)
!!$     if (myid ==0) then
!!$        open(26,file='stalta.bin',access='direct',recl=cp*ntold)
!!$        write(26,rec=1)stalta
!!$        close(26)
!!$     end if
!!$
!!$     if (istap == 0) then
!!$	write(6,*)'verif stalta, indice is: '
!!$	write(6,*)ind - 1 
!!$        write(6,*)'corresponding to time step:'
!!$	write(6,*)(ind - 1)* dtold
!!$	write(6,*)'dont quit but beware that stalta may have failed'
!!$        warning =1
!!$	call finalize_mpi
!!$	stop 
!!$     else 
!!$	ind = alpha
!!$     end if
!!$
!!$     !  do i=1,ntold
!!$     !     if (stalta(i) >= thres) then
!!$     !        ind = i;
!!$     !        exit
!!$     !     end if
!!$     !  end do
!!$     !  if (myid ==0) then
!!$     ! open(26,file='verifconv.bin',access='direct',recl=cp*ntold)
!!$     ! write(26,rec=1)vpow
!!$     ! close(26)
!!$     ! end if
!!$
!!$     deallocate(vpow)
!!$     deallocate(stalta)
!!$
!!$
!!$     if (myid==0) write(6,*)'Done'

     !================================================================================
     ! Cut signal and convolve with stf
     !--------------------------------------------------
     !*** Convolve with stf (could also be a filter)
     if (isconv == 1) then  
        !* Oldlen
        oldlen = ntime 
        if(.not.allocated(taper)) allocate(taper(ntold))
        alph1 = 10/(dtold*oldlen)          ! Warning 2s taper
        taper = tuckeywin(oldlen,alph1)

        !* Allocate convtmp
        if (.not.allocated(convtmpvx))  allocate(convtmpvx(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpvy))  allocate(convtmpvy(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpvz))  allocate(convtmpvz(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpsxx)) allocate(convtmpsxx(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpsyy)) allocate(convtmpsyy(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpszz)) allocate(convtmpszz(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpsyz)) allocate(convtmpsyz(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpsxz)) allocate(convtmpsxz(nrec_to_store,oldlen+ntstf-1))
        if (.not.allocated(convtmpsxy)) allocate(convtmpsxy(nrec_to_store,oldlen+ntstf-1))
        convtmpvx = 0.
        convtmpvy = 0.
        convtmpvz = 0.
        convtmpsxx = 0.
        convtmpsyy = 0.
        convtmpszz = 0.
        convtmpsyz = 0.
        convtmpsxz = 0.
        convtmpsxy = 0.

        !* convolve
        do it = 1, ntold
           do itstf = 1, ntstf
              do ipt = 1, nrec_to_store
                 convtmpvx( ipt,it+itstf-1) = convtmpvx( ipt,it+itstf-1) + vxold1(ipt,it)  * stf(itstf) * taper(it)
                 convtmpvy( ipt,it+itstf-1) = convtmpvy( ipt,it+itstf-1) + vyold1(ipt,it)  * stf(itstf) * taper(it)
                 convtmpvz( ipt,it+itstf-1) = convtmpvz( ipt,it+itstf-1) + vzold1(ipt,it)  * stf(itstf) * taper(it)
                 convtmpsxx(ipt,it+itstf-1) = convtmpsxx(ipt,it+itstf-1) + sxxold1(ipt,it) * stf(itstf) * taper(it)
                 convtmpsyy(ipt,it+itstf-1) = convtmpsyy(ipt,it+itstf-1) + syyold1(ipt,it) * stf(itstf) * taper(it)
                 convtmpszz(ipt,it+itstf-1) = convtmpszz(ipt,it+itstf-1) + szzold1(ipt,it) * stf(itstf) * taper(it)
                 convtmpsyz(ipt,it+itstf-1) = convtmpsyz(ipt,it+itstf-1) + syzold1(ipt,it) * stf(itstf) * taper(it)
                 convtmpsxz(ipt,it+itstf-1) = convtmpsxz(ipt,it+itstf-1) + sxzold1(ipt,it) * stf(itstf) * taper(it)
                 convtmpsxy(ipt,it+itstf-1) = convtmpsxy(ipt,it+itstf-1) + sxyold1(ipt,it) * stf(itstf) * taper(it)
              end do
           end do
           if(myid==0) write(6,*)'conv ',it,ntold
        end do

        !* take good part (wrong)
        if (modulo(ntstf,2) == 0) then
           ind2 = ntstf/2 + 1
        else
           ind2 = ceiling(real(ntstf/2,kind=cp))
        end if
        !   ind2 = 1

        vxold1(:,:)  = convtmpvx(:,ind2:ind2+ntold-1)  * dtold
        vyold1(:,:)  = convtmpvy(:,ind2:ind2+ntold-1)  * dtold
        vzold1(:,:)  = convtmpvz(:,ind2:ind2+ntold-1)  * dtold
        sxxold1(:,:) = convtmpsxx(:,ind2:ind2+ntold-1) * dtold
        syyold1(:,:) = convtmpsyy(:,ind2:ind2+ntold-1) * dtold
        szzold1(:,:) = convtmpszz(:,ind2:ind2+ntold-1) * dtold
        syzold1(:,:) = convtmpsyz(:,ind2:ind2+ntold-1) * dtold
        sxzold1(:,:) = convtmpsxz(:,ind2:ind2+ntold-1) * dtold
        sxyold1(:,:) = convtmpsxy(:,ind2:ind2+ntold-1) * dtold

        !* Deallocate
        deallocate(convtmpvx)
        deallocate(convtmpvy)
        deallocate(convtmpvz)
        deallocate(convtmpsxx)
        deallocate(convtmpsyy)
        deallocate(convtmpszz)
        deallocate(convtmpsyz)
        deallocate(convtmpsxz)
        deallocate(convtmpsxy)

        if(allocated(taper)) deallocate(taper)

     end if


     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)

     if (myid==0) write(6,*)'Done.'



     !================================================================================
     !*** Filter seismograms
     if (istap == 1) then
        if(myid==0) write(6,*)'Filtering...'       
        if(.not.allocated(taper)) allocate(taper(ntime))
        if(.not.allocated(convfilt)) allocate(convfilt(ntime))
        convfilt = 0._cp
        alph1 = 10/(dtold*ntime)          ! Warning 2s taper
        taper = tuckeywin(ntime,alph1)
        

        if (myid == 0) print *,myid,fmax,dtold,0.5/dtold 
        do ipt = 1, nrec_to_store
           
           !*** Taper on residuals
           vxold1(ipt,:) = vxold1(ipt,:) * taper(:)
           call bwfilt(vxold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           vxold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           vyold1(ipt,:) = vyold1(ipt,:) * taper(:)
           call bwfilt(vyold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           vyold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           vzold1(ipt,:) = vzold1(ipt,:) * taper(:)
           call bwfilt(vzold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           vzold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           sxxold1(ipt,:) = sxxold1(ipt,:) * taper(:)
           call bwfilt(sxxold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           sxxold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           syyold1(ipt,:) = syyold1(ipt,:) * taper(:)
           call bwfilt(syyold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           syyold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           szzold1(ipt,:) = szzold1(ipt,:) * taper(:)
           call bwfilt(szzold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           szzold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           syzold1(ipt,:) = syzold1(ipt,:) * taper(:)
           call bwfilt(syzold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           syzold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           sxzold1(ipt,:) = sxzold1(ipt,:) * taper(:)
           call bwfilt(sxzold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           sxzold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           sxyold1(ipt,:) = sxyold1(ipt,:) * taper(:)
           call bwfilt(sxyold1(ipt,:),convfilt,dtold,ntime,1,4,1e-3_cp,fmax)
           sxyold1(ipt,:) = convfilt(:)
           convfilt = 0._cp   
           
        end do

        if (allocated(convfilt)) deallocate(convfilt)
        if (allocated(taper)) deallocate(taper)
        if(myid==0) write(6,*)'Done!'
     end if


     !*** Cut old signal
!!$     itbeg  = ind !-ceiling(2.5/fmax)
!!$     itend  = itbeg + ceiling(ntnew * dtnew / dtold)
!!$     tbeg = itbeg * dtold      !*** Starting time of cut signal
!!$     tend = itend * dtold
!!$     oldlen = itend-itbeg+1
!!$
!!$
!!$     if (myid==0) write(6,*)'Infos : itbeg, itend, tbeg, tend, oldlen, ntold, dtold'
!!$     if (myid==0) write(6,*)itbeg, itend, tbeg, tend, oldlen, ntold, dtold
!!$
!!$!!!!!!!! Interpolate stf again
!!$!    if (isconv==1) then
!!$!    write(6,*)'Interpolate source time function...'
!!$!    if(.not.allocated(stfn)) allocate(stfn(oldlen))
!!$!    stf(:) = 0._cp
!!$!    festf = 1. / dtold
!!$! 
!!$!    do itnew = 1, oldlen
!!$!       do itold = 1, ntstf
!!$!          stfn(itnew) = stfn(itnew) + stf(itold) * mysinc(real(festf * itnew * dtold - itold))
!!$!       end do
!!$!    end do
!!$!    if(allocated(tmpstf)) deallocate(tmpstf)
!!$!
!!$!    ntstf = ntold
!!$!
!!$!    write(6,*)'Done !'
!!$!   end if
!!$!$     !*** Interpolate to axisem 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(vxold)) allocate(vxold(nrec_to_store,oldlen))
     vxold(:,:) = vxold1(:,:)
     deallocate(vxold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(vyold)) allocate(vyold(nrec_to_store,oldlen))
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     vyold(:,:) = vyold1(:,:)
     deallocate(vyold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(vzold)) allocate(vzold(nrec_to_store,oldlen))
     vzold(:,:) = vzold1(:,:)
     deallocate(vzold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(sxxold)) allocate(sxxold(nrec_to_store,oldlen))
     sxxold(:,:) = sxxold1(:,:)
     deallocate(sxxold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(syyold)) allocate(syyold(nrec_to_store,oldlen))
     syyold(:,:) = syyold1(:,:)
     deallocate(syyold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(szzold)) allocate(szzold(nrec_to_store,oldlen))
     szzold(:,:) = szzold1(:,:)
     deallocate(szzold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(syzold)) allocate(syzold(nrec_to_store,oldlen))
     syzold(:,:) = syzold1(:,:)
     deallocate(syzold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(sxzold)) allocate(sxzold(nrec_to_store,oldlen))
     sxzold(:,:) = sxzold1(:,:)
     deallocate(sxzold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (.not.allocated(sxyold)) allocate(sxyold(nrec_to_store,oldlen))
     sxyold(:,:) = sxyold1(:,:)
     deallocate(sxyold1)
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)

     ntold = oldlen

     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)



!     write(ficii,'(a,i4.4,a)')'convstf_',myid,'.dat' 
!     open(20000+myid,file=trim(ficii),action='write')
!     write(20000+myid,*)vzold(1,:)
!     close(20000+myid)

!     write(ficii,'(a,i4.4,a)')'filtconvstf_',myid,'.dat' 
!     open(20000+myid,file=trim(ficii),action='write')
!     write(20000+myid,*)vzold(1,:)
!     close(20000+myid)


     !================================================================================
     ! Resample with sinc interpolation
     !--------------------------------------------------
     if (myid==0) write(6,*)'Resample with sinc interpolation....'

     !*** Open file and allocate
     if (myid==0) then
        if (fwdtool == 'SEM') then
           write(srcfile,*)'incident_field'
           write(myfileend1,'(a,i4.4)')'_vel_',ipart
           open(73000,file=trim(srcfile)//myfileend1,status='replace',access='direct',recl=cp*3*ngllsquare*num_bnd_faces*tbuff)
           write(myfileend1,'(a,i4.4)')'_tra_',ipart
           open(73001,file=trim(srcfile)//myfileend1,status='replace',access='direct',recl=cp*3*ngllsquare*num_bnd_faces*tbuff)
        else
           open(20,file='incident_field.bin',status='unknown',form='unformatted')
        end if
     end if
     if (fwdtool == 'SEM') then
        if (allocated(vel_inc)) deallocate(vel_inc)
        if (allocated(trac_inc)) deallocate(trac_inc)
        if (.not.allocated(vel_inc))    allocate(vel_inc(3,ngllsquare,num_bnd_faces,tbuff))
        if (.not.allocated(trac_inc)) allocate(trac_inc(3,ngllsquare,num_bnd_faces,tbuff))
        vel_inc  = 0.
        trac_inc = 0.

!!! Add local versions
        if (allocated(mapipt)) deallocate(mapipt)
        if (.not.allocated(mapipt)) allocate(mapipt(2,ngllsquare*num_bnd_faces))
        mapipt = 0    

        !*** Map igll, iface to ipt
        ipt = 0
        do iface = 1, num_bnd_faces
           do igll = 1, ngllsquare
              ipt = ipt + 1
              mapipt(1,ipt) = igll
              mapipt(2,ipt) = iface
           end do
        end do
     else
        if (allocated(vel_inc2)) deallocate(vel_inc2)
        if (allocated(stress_inc)) deallocate(stress_inc)
        if (.not.allocated(vel_inc2))    allocate(vel_inc2(3,npts,tbuff))     ! warning npts may be wrong
        if (.not.allocated(stress_inc)) allocate(stress_inc(6,npts,tbuff))    ! warning "
        vel_inc2   = 0.
        stress_inc = 0.
     end if
     if (allocated(tab_sinc)) deallocate(tab_sinc)
     if (.not.allocated(tab_sinc)) allocate(tab_sinc(ntold))
     feold = 1./dtold

     if (myid==0) write(6,*)'Infos : feold, ntold, itbeg, itend, tbeg, tend, dtnew, dtold,ntnew'
     if (myid==0) write(6,*)feold,ntold,itbeg,itend,tbeg,tend,dtnew,dtold,ntnew

!     print *,myid,fwdtool
     
     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)

     !*** Loop over new time steps
     do itnew = 1, ntnew

        !*** Locate in buffer
        ibuf = modulo(itnew-1,tbuff) +1    ! Indice in buffer
        nbuf = ((itnew-1)/tbuff) +1        ! Current number of buffer

 !       if (fwdtool == 'SEM') then
 !         vel_inc  = 0.
 !         trac_inc = 0.
 !       else 
 !          vel_inc2   = 0.
 !          stress_inc = 0.
 !       end if

        !*** Compute sinc kernel
        call comp_tab_sinc(itnew,dtnew,feold,ntold,tab_sinc)
        
        do ipt = 1, nrec_to_store

           !*** Loop over old time steps
           vx  = sum( vxold(ipt,:)  * tab_sinc(:) )
           vy  = sum( vyold(ipt,:)  * tab_sinc(:) )
           vz  = sum( vzold(ipt,:)  * tab_sinc(:) )
           sxx = sum( sxxold(ipt,:) * tab_sinc(:) )
           syy = sum( syyold(ipt,:) * tab_sinc(:) )
           szz = sum( szzold(ipt,:) * tab_sinc(:) )
           syz = sum( syzold(ipt,:) * tab_sinc(:) )
           sxz = sum( sxzold(ipt,:) * tab_sinc(:) )
           sxy = sum( sxyold(ipt,:) * tab_sinc(:) )

           !*** Compute traction for sem
           select case(fwdtool)
           case ('SEM') 

              !* 1. Indices
              iptglob = ipt + i_inf(myid+1) - tab_box_rec(2,ipart)  !ipt + i_inf(myid+1) - 1 !irecmin + ipt - 1 
              igll    = mapipt(1,iptglob)
              iface   = mapipt(2,iptglob)

              !Get local indices for GLL point
              ielem = abs_bnd_ielem(iface)
              i = abs_bnd_ijk(1,igll,iface)
              j = abs_bnd_ijk(2,igll,iface)
              k = abs_bnd_ijk(3,igll,iface)
              iglob = loc2glob(i,j,k,ielem)

              !coord
              x=xcoord(iglob)
              y=ycoord(iglob)
              z=zcoord(iglob)

              !* 2. Normals
              nx = sign(abs_bnd_normal(1,igll,iface),x)
              ny = sign(abs_bnd_normal(2,igll,iface),y)
              nz = sign(abs_bnd_normal(3,igll,iface),z)
!              nx = abs_bnd_normal(1,igll,iface)
!              ny = abs_bnd_normal(2,igll,iface)
!              nz = abs_bnd_normal(3,igll,iface)

              !* 3. Tractions
              tx = sxx*nx + sxy*ny + sxz*nz
              ty = sxy*nx + syy*ny + syz*nz
              tz = sxz*nx + syz*ny + szz*nz

              !* 4. Stock
              vel_inc(1,igll,iface,ibuf) = vx
              vel_inc(2,igll,iface,ibuf) = vy
              vel_inc(3,igll,iface,ibuf) = vz

              trac_inc(1,igll,iface,ibuf) = tx
              trac_inc(2,igll,iface,ibuf) = ty
              trac_inc(3,igll,iface,ibuf) = tz

           case('FD')

              iptglob = ipt + i_inf(myid+1) -1

              !*** Stock everything else
              vel_inc2(1,iptglob,ibuf) = vx
              vel_inc2(2,iptglob,ibuf) = vy
              vel_inc2(3,iptglob,ibuf) = vz

              stress_inc(1,iptglob,ibuf) = sxx
              stress_inc(2,iptglob,ibuf) = syy
              stress_inc(3,iptglob,ibuf) = szz
              stress_inc(4,iptglob,ibuf) = syz
              stress_inc(5,iptglob,ibuf) = sxz
              stress_inc(6,iptglob,ibuf) = sxy

           end select

        end do

        
        if (ibuf == tbuff) then
           
           !*** Write everything (check engine)
           select case (fwdtool)
           case('SEM')
              call MPI_allreduce(MPI_IN_PLACE,vel_inc,3*ngllsquare*num_bnd_faces*tbuff,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
              call MPI_allreduce(MPI_IN_PLACE,trac_inc,3*ngllsquare*num_bnd_faces*tbuff,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)   
              if (myid ==0) then 
                 
                 write(6,*)'Progress : ',100.*itnew/ntnew,'%'   
                 
                 !* Save buffer to disk
                 write(73000,rec=nbuf)vel_inc(:,:,:,:)
                 write(73001,rec=nbuf)trac_inc(:,:,:,:)
                 
              end if

             !* Reinit to zero
             vel_inc  = 0.
             trac_inc = 0.
                 
                            

           case('FD')        
              call MPI_allreduce(MPI_IN_PLACE,vel_inc2,3*npts*tbuff,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
              call MPI_allreduce(MPI_IN_PLACE,stress_inc,6*npts*tbuff,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr_mpi) 
              if (myid == 0) write(20)vel_inc2,stress_inc
           end select
           
        end if

     end do

     call MPI_barrier(MPI_COMM_WORLD,ierr_mpi)
     if (myid==0) then
        if (fwdtool == 'SEM') then
           close(73000)
           close(73001)
        end if
     else
        close(20)
     end if
     
     if (myid==0) write(6,*)'End of incident field computation for ',ipart
     

     if(allocated(vxold)) deallocate(vxold)
     if(allocated(vyold)) deallocate(vyold)
     if(allocated(vzold)) deallocate(vzold)
     if(allocated(sxxold)) deallocate(sxxold)
     if(allocated(syyold)) deallocate(syyold)
     if(allocated(szzold)) deallocate(szzold)
     if(allocated(sxyold)) deallocate(sxyold)
     if(allocated(sxzold)) deallocate(sxzold)
     if(allocated(syzold)) deallocate(syzold)
     
     
  end do
  !  if (warning == 1) then      
!  if (myid==0) write(*,*)'WARNING automatic picking has been used it may be wrong...' 
!  if (myid==0) write(6,*)'WARNING please verify.'
!  end if
        

  call finish_program

contains

  subroutine begin_program

    call init_mpi

  end subroutine begin_program

  subroutine finish_program

    call finalize_mpi

  end subroutine finish_program

end program
