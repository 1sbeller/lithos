program interpolate_3D_wavefield

  use precision_mod
  use constants_mod
  use mpi_int_mod
  use interpolation_parameters_mod
  use interp_process_mod
  use inputs_outputs_interp_mod

  implicit none

  !================================================================================
  ! Prepare reconstruction
  !--------------------------------------------------
  !*** Begin program and mpi
  call begin_program
  
  !*** Read inut files parameters
  call read_all_inputs
  call broadcast_all_data

  !*** Prepare reading of reconstructed AxiSEM outputs
  if (myid==0) write(6,*)'Must read ',nbrec,' points for ',ntime,' time steps.'
  if (myid==0) write(6,*)'Check ntold and ntime : ',ntold,ntime
  npts = nbrec

  if (myid == 0) then
     if(.not.allocated(data_tmp)) allocate(data_tmp(nbrec,9))
     if (.not.allocated(vxold2)) allocate(vxold2(npts,1))
     if (.not.allocated(vyold2)) allocate(vyold2(npts,1))
     if (.not.allocated(vzold2)) allocate(vzold2(npts,1))
     if (.not.allocated(sxxold2)) allocate(sxxold2(npts,1))
     if (.not.allocated(syyold2)) allocate(syyold2(npts,1))
     if (.not.allocated(szzold2)) allocate(szzold2(npts,1))
     if (.not.allocated(syzold2)) allocate(syzold2(npts,1))
     if (.not.allocated(sxzold2)) allocate(sxzold2(npts,1))
     if (.not.allocated(sxyold2)) allocate(sxyold2(npts,1))
     
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

  call scatter_data

  if (.not.allocated(vxold1)) allocate(vxold1(nrec_to_store,ntold))
  if (.not.allocated(vyold1)) allocate(vyold1(nrec_to_store,ntold))
  if (.not.allocated(vzold1)) allocate(vzold1(nrec_to_store,ntold))
  if (.not.allocated(sxxold1)) allocate(sxxold1(nrec_to_store,ntold))
  if (.not.allocated(syyold1)) allocate(syyold1(nrec_to_store,ntold))
  if (.not.allocated(szzold1)) allocate(szzold1(nrec_to_store,ntold))
  if (.not.allocated(syzold1)) allocate(syzold1(nrec_to_store,ntold))
  if (.not.allocated(sxzold1)) allocate(sxzold1(nrec_to_store,ntold))
  if (.not.allocated(sxyold1)) allocate(sxyold1(nrec_to_store,ntold))
  
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
  if (myid==0) write(6,*)'Read AxisEM files...'
  do itime=1,ntime

     if (myid == 0) then ! Read on master proc
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
           
           vxold2(:,1)  = vxold2(:,1) + real(data_tmp(:,1))
           vyold2(:,1)  = vyold2(:,1) + real(data_tmp(:,2))
           vzold2(:,1)  = vzold2(:,1) + real(data_tmp(:,3))
           sxxold2(:,1) = sxxold2(:,1) + real(data_tmp(:,4))
           syyold2(:,1) = syyold2(:,1) + real(data_tmp(:,5))
           szzold2(:,1) = szzold2(:,1) + real(data_tmp(:,6))
           sxyold2(:,1) = sxyold2(:,1) + real(data_tmp(:,7))
           sxzold2(:,1) = sxzold2(:,1) + real(data_tmp(:,8))
           syzold2(:,1) = syzold2(:,1) + real(data_tmp(:,9))

        end do
     end if

     !*** Send-receive to scatter data
     if (myid ==0) then ! SEND

        vxold1(:,itime)  =  vxold2(i_inf(1):i_sup(1),1)
        vyold1(:,itime)  =  vyold2(i_inf(1):i_sup(1),1)
        vzold1(:,itime)  =  vzold2(i_inf(1):i_sup(1),1)
        sxxold1(:,itime) = sxxold2(i_inf(1):i_sup(1),1)
        syyold1(:,itime) = syyold2(i_inf(1):i_sup(1),1)
        szzold1(:,itime) = szzold2(i_inf(1):i_sup(1),1)
        sxyold1(:,itime) = sxyold2(i_inf(1):i_sup(1),1)
        sxzold1(:,itime) = sxzold2(i_inf(1):i_sup(1),1)
        syzold1(:,itime) = syzold2(i_inf(1):i_sup(1),1)

        do iproc=1,nb_proc-1
        
           call mpi_send( vxold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send( vyold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send( vzold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send(sxxold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send(syyold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send(szzold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send(sxyold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send(sxzold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
           call mpi_send(syzold2(i_inf(myid+1):i_sup(myid+1),1),nb_received_sv(myid+1),MPICP,myid,etq,MPI_COMM_WORLD,ierr_mpi)
     
        end do
        
     else ! RECEIVE
     
        call mpi_send( vxold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send( vyold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send( vzold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send(sxxold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send(syyold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send(szzold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send(sxyold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send(sxzold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)
        call mpi_send(syzold1(:,itime),nb_received_sv(myid+1),MPICP,0,etq,MPI_COMM_WORLD,statut,ierr_mpi)

     end if

  end do
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
  
  if (myid ==0) then
     deallocate(vxold2)
     deallocate(vyold2)
     deallocate(vzold2)
     deallocate(sxxold2)
     deallocate(syyold2)
     deallocate(szzold2)
     deallocate(syzold2)
     deallocate(sxzold2)
     deallocate(sxyold2)
  end if

  !================================================================================
  ! Compute velocity energy and STA/LTA this could be done on many procs)
  !--------------------------------------------------
  !*** Compute energy
  if (myid==0) write(6,*)'Compute STA/LTA...'
  if (.not.allocated(vpow)) allocate(vpow(ntold))
  if (.not.allocated(stalta)) allocate(stalta(ntold))
  vpow = 0.
  stalta = 0.
  
  do itold = 1, ntold
     sumx = 0.
     sumy = 0.
     sumz = 0.
     do ipt =  1, nrec_to_store  !irecmin, irecmax
        sumx = sumx + vxold1(ipt,itold)**2
        sumy = sumy + vyold1(ipt,itold)**2
        sumz = sumz + vzold1(ipt,itold)**2
     end do
    call MPI_allreduce(MPI_IN_PLACE,sumx,1,MPICP,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
    call MPI_allreduce(MPI_IN_PLACE,sumy,1,MPICP,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
    call MPI_allreduce(MPI_IN_PLACE,sumz,1,MPICP,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
    vpow(itold) = sqrt(sumx + sumy +sumz)
  end do

  !*** Compute STA/LTA
  call substalta(vpow, nsta, nlta, stalta)
  do i=1,ntold
     if (stalta(i) >= thres) then
        ind = i;
        exit
     end if
  end do
  deallocate(vpow)
  deallocate(stalta)
  if (myid==0) write(6,*)'Done'

  !================================================================================
  ! Cut signal and convolve with stf
  !--------------------------------------------------
  !*** Cut old signal
  itbeg  = i-1
  itend  = itbeg + ceiling(ntnew * dtnew / dtold)
  tbeg = itbeg * dtold      !*** Starting time of cut signal
  tend = itend * dtold
  oldlen = itend-itbeg+1
  
  if (myid==0) write(6,*)'Infos : itbeg, itend, tbeg, tend, oldlen'
  if (myid==0) write(6,*)itbeg, itend, tbeg, tend, oldlen

  if (.not.allocated(vxold)) allocate(vxold(nrec_to_store,oldlen))
  if (.not.allocated(vyold)) allocate(vyold(nrec_to_store,oldlen))
  if (.not.allocated(vzold)) allocate(vzold(nrec_to_store,oldlen))
  if (.not.allocated(sxxold)) allocate(sxxold(nrec_to_store,oldlen))
  if (.not.allocated(syyold)) allocate(syyold(nrec_to_store,oldlen))
  if (.not.allocated(szzold)) allocate(szzold(nrec_to_store,oldlen))
  if (.not.allocated(syzold)) allocate(syzold(nrec_to_store,oldlen))
  if (.not.allocated(sxzold)) allocate(sxzold(nrec_to_store,oldlen))
  if (.not.allocated(sxyold)) allocate(sxyold(nrec_to_store,oldlen))

  vxold(:,:) = vxold1(:,itbeg:itend)
  vyold(:,:) = vyold1(:,itbeg:itend)
  vzold(:,:) = vzold1(:,itbeg:itend)
  sxxold(:,:) = sxxold1(:,itbeg:itend)
  syyold(:,:) = syyold1(:,itbeg:itend)
  szzold(:,:) = szzold1(:,itbeg:itend)
  syzold(:,:) = syzold1(:,itbeg:itend)
  sxzold(:,:) = sxzold1(:,itbeg:itend)
  sxyold(:,:) = sxyold1(:,itbeg:itend)

  deallocate(vxold1)
  deallocate(vyold1)
  deallocate(vzold1)
  deallocate(sxxold1)
  deallocate(syyold1)
  deallocate(szzold1)
  deallocate(syzold1)
  deallocate(sxzold1)
  deallocate(sxyold1)

  ntold = oldlen

  !*** Convolve with stf (could also be a filter)
  !!!!! Allocate stf
  !!!!! Read stf
  !!!!! Convolve (can be done on multi procs also)
  !  do ipt=1,npts
  !          
  !     call myconvolve(dble(vxold(ipt,:)),ker,nt,n1,conv)
  !     vxold(ipt,:) = real(conv(ibeg:iend))
  !     call myconvolve(dble(vyold(ipt,:)),ker,nt,n1,conv)
  !     vyold(ipt,:) = real(conv(ibeg:iend))
  !     call myconvolve(dble(vzold(ipt,:)),ker,nt,n1,conv)
  !     vzold(ipt,:) = real(conv(ibeg:iend))
  !     
  !     call myconvolve(dble(sxxold(ipt,:)),ker,nt,n1,conv)
  !     sxxold(ipt,:) = real(conv(ibeg:iend))
  !     call myconvolve(dble(syyold(ipt,:)),ker,nt,n1,conv)
  !     syyold(ipt,:) = real(conv(ibeg:iend))
  !     call myconvolve(dble(szzold(ipt,:)),ker,nt,n1,conv)
  !     szzold(ipt,:) = real(conv(ibeg:iend))
  !     call myconvolve(dble(syzold(ipt,:)),ker,nt,n1,conv)
  !     syzold(ipt,:) = real(conv(ibeg:iend))
  !     call myconvolve(dble(sxzold(ipt,:)),ker,nt,n1,conv)
  !     sxzold(ipt,:) = real(conv(ibeg:iend))
  !     call myconvolve(dble(sxyold(ipt,:)),ker,nt,n1,conv)
  !     sxyold(ipt,:) = real(conv(ibeg:iend))
  !     
  !  end do
  if (myid==0) write(6,*)'Done.'


  !================================================================================
  ! Resample with sinc interpolation
  !--------------------------------------------------
  if (myid==0) write(6,*)'Resample with sinc interpolation....'

  !*** Open file and allocate
  open(20,file='incident_field.bin',status='unknown',form='unformatted')
  if (fwdtool == 'SEM') then
     if (.not.allocated(vel_inc))    allocate(vel_inc(3,ngllsquare,num_bnd_faces))
     if (.not.allocated(trac_inc)) allocate(trac_inc(3,ngllsquare,num_bnd_faces))
     vel_inc  = 0.
     trac_inc = 0.
     
     !!! Add local versions

     
     if (.not.allocated(mapipt)) allocate(mapipt(2,npts))
     
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
     if (.not.allocated(vel_inc2))    allocate(vel_inc2(3,npts))
     if (.not.allocated(stress_inc)) allocate(stress_inc(6,npts))
     vel_inc2   = 0.
     stress_inc = 0.
  end if
  if (.not.allocated(tab_sinc)) allocate(tab_sinc(ntold))
  feold = 1./dtold

  if (myid==0) write(6,*)'Infos : feold, ntold, itbeg, itend, tbeg, tend, dtnew, dtold'
  if (myid==0) write(6,*)feold,ntold,itbeg,itend,tbeg,tend,dtnew,dtold

  !*** Loop over new time steps
  do itnew = 1, ntnew

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
           iptglob = ipt + i_inf(myid+1) -1
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
           !           nx = abs_bnd_normal(1,igll,iface)
           !           ny = abs_bnd_normal(2,igll,iface)
           !           nz = abs_bnd_normal(3,igll,iface)
           
           !* 3. Tractions
           tx = sxx*nx + sxy*ny + sxz*nz
           ty = sxy*nx + syy*ny + syz*nz
           tz = sxz*nx + syz*ny + szz*nz
           
           !* 4. Stock
           vel_inc(1,igll,iface) = vx
           vel_inc(2,igll,iface) = vy
           vel_inc(3,igll,iface) = vz
           
           trac_inc(1,igll,iface) = tx
           trac_inc(2,igll,iface) = ty
           trac_inc(3,igll,iface) = tz
           
        case('FD')
           
           iptglob = ipt + i_inf(myid+1) -1

           !*** Stock everything else
           vel_inc2(1,iptglob) = vx
           vel_inc2(2,iptglob) = vy
           vel_inc2(3,iptglob) = vz
           
           stress_inc(1,iptglob) = sxx
           stress_inc(2,iptglob) = syy
           stress_inc(3,iptglob) = szz
           stress_inc(4,iptglob) = syz
           stress_inc(5,iptglob) = sxz
           stress_inc(6,iptglob) = sxy
           
        end select

     end do
     
     if(myid ==0 .and. mod(100*itnew/ntnew,5)==0) write(6,*)'Progress : ',100*itnew/ntnew,'%'

     !*** Write everything (check engine)
     select case (fwdtool)
     case('SEM')
        call MPI_allreduce(MPI_IN_PLACE,vel_inc,3*ngllsquare*num_bnd_faces,MPICP,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
        call MPI_allreduce(MPI_IN_PLACE,trac_inc,3*ngllsquare*num_bnd_faces,MPICP,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)   
        if (myid ==0) write(20)vel_inc,trac_inc
     case('FD')        
        call MPI_allreduce(MPI_IN_PLACE,vel_inc2,3*npts,MPICP,MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
        call MPI_allreduce(MPI_IN_PLACE,stress_inc,6*npts,MPICP,MPI_SUM,MPI_COMM_WORLD,ierr_mpi) 
        if (myid == 0) write(20)vel_inc2,stress_inc
     end select

  end do

  close(20)

  if (myid==0) write(*,*)'End of incident field computation.'


  call finish_program

contains

  subroutine begin_program

    call init_mpi

  end subroutine begin_program

  subroutine finish_program

    call finalize_mpi

  end subroutine finish_program

end program
