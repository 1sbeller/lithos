module inputs_outputs_mod

  use precision_mod

  implicit none

contains

!################################################################################
!############################# INPUTS FILES #####################################
!################################################################################

!================================================================================
! Routine read_all_inputs used in the main program
  subroutine read_all_inputs(isim, nsim_to_send)

    integer(kind=si), intent(inout) :: isim
    integer(kind=si), intent(inout) :: nsim_to_send

    !*** Read simulations infos
    call read_info_simu(nsim_to_send)

    !*** Read coupling information
    call read_reconstruction_inputs(isim)

    !*** Read mesh infos
    call read_mesh(isim)

  end subroutine read_all_inputs
!--------------------------------------------------------------------------------


!================================================================================
! Read informations for each source of the axisem simulation
  subroutine read_info_simu(nsim_to_send)
  
    use global_parameters_mod, only :nsim,simdir,src_type,ntime,working_axisem_dir,magnitude,Mij
    
    integer(kind=si), intent(inout) :: nsim_to_send

    logical :: use_netcdf

    integer(kind=si), parameter  :: i_param_post=500
    integer(kind=si)             :: iostat, ioerr, isim, iinparam_source
    real(kind=cp)                ::  tshift
    character(len=12)            :: src_file_type
    character(len=256)  :: keyword, keyvalue, line, junk
    
    character(len=100), allocatable, dimension(:) :: bkgrndmodel, rot_rec
    character(len=7),   allocatable, dimension(:) :: stf_type
    
    real(kind=cp), allocatable, dimension(:) :: dt_tmp, period_tmp, dt_seis_tmp
    real(kind=cp), allocatable, dimension(:) :: dt_strain, dt_snap, nt_snap
    real(kind=cp), allocatable, dimension(:) :: srccolat_tmp, srclon_tmp, src_depth_tmp
    real(kind=cp), allocatable, dimension(:) :: shift_fact_tmp, mag

    integer(kind=si), allocatable, dimension(:) :: ishift_deltat, ishift_seisdt, ishift_straindt
    integer(kind=si), allocatable, dimension(:) :: ibeg, iend
    integer(kind=si), allocatable, dimension(:) :: nt, nt_strain, nrec_tmp, nt_seis_tmp

    real(kind=cp) :: amplitude


    !*** Define working directory
    working_axisem_dir='./'

    !*** Open and read inparam_basic
    open(unit=i_param_post, file='inparam_basic', status='old', action='read', &
         iostat=ioerr)
    do
       read(i_param_post, fmt='(a256)', iostat=ioerr) line
       if (ioerr < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
       
       read(line,*) keyword, keyvalue 
       select case(trim(keyword))
       case('SIMULATION_TYPE')
          if (keyvalue == 'single') then
             nsim = 1
             allocate(simdir(nsim))
             simdir(1) = "./"
          elseif (keyvalue == 'moment') then
             nsim = 4
             allocate(simdir(nsim))
             simdir(1) = "MZZ/"
             simdir(2) = "MXX_P_MYY/"
             simdir(3) = "MXZ_MYZ/"
             simdir(4) = "MXY_MXX_M_MYY/"
          elseif (keyvalue == 'force') then
             write(6,*) 'postprocessing for "force" simulation needs work!'
             stop 2
          endif
       end select
    end do
    nsim_to_send=nsim
    close(i_param_post)
    
    !*** Define tshif
    tshift = 0.
    
    !*** Allocate arrays for each simulation
    allocate(bkgrndmodel(nsim), stf_type(nsim))
    allocate(src_type(nsim,2))
    allocate(dt_tmp(nsim), period_tmp(nsim), mag(nsim), magnitude(nsim), dt_seis_tmp(nsim))
    allocate(dt_strain(nsim), dt_snap(nsim))
    allocate(rot_rec(nsim))
    allocate(nt(nsim), nrec_tmp(nsim), nt_seis_tmp(nsim), nt_strain(nsim), nt_snap(nsim))
    allocate(ibeg(nsim), iend(nsim), srccolat_tmp(nsim), srclon_tmp(nsim), src_depth_tmp(nsim))
    allocate(shift_fact_tmp(nsim))
    allocate(ishift_deltat(nsim), ishift_seisdt(nsim), ishift_straindt(nsim))
    allocate(Mij(nsim,1:6)) 
    Mij(:,:) = 0.
    
    !*** For each simulation, read the file SIMULATION.INFO and the the source parameter info 
    do isim = 1,nsim
       open(unit=99,file=trim(simdir(isim))//'/simulation.info')
       read(99,*) bkgrndmodel(isim)
       read(99,*) dt_tmp(isim)
       read(99,*) nt(isim)
       read(99,*) src_type(isim,1)
       read(99,*) src_type(isim,2)
       read(99,*) stf_type(isim)
       read(99,*) src_file_type
       read(99,*) period_tmp(isim)
       read(99,*) src_depth_tmp(isim)
       read(99,*) srccolat_tmp(isim)
       read(99,*) srclon_tmp(isim)
       read(99,*) mag(isim)
       read(99,*) nrec_tmp(isim)
       read(99,*) nt_seis_tmp(isim)
       read(99,*) dt_seis_tmp(isim)
       read(99,*) nt_strain(isim)
       read(99,*) dt_strain(isim)
       read(99,*) nt_snap(isim)
       read(99,*) dt_snap(isim)
       read(99,*) rot_rec(isim)
       read(99,*) ibeg(isim)
       read(99,*) iend(isim)
       read(99,*) shift_fact_tmp(isim)
       read(99,*) ishift_deltat(isim)
       read(99,*) ishift_seisdt(isim)
       read(99,*) ishift_straindt(isim)
       read(99,*) 
       read(99,*)
       read(99,*) use_netcdf
       close(99)


       !*** Tensor moment
       magnitude(isim) = mag(isim)

       select case(src_file_type)
       case('moment')
          write(6,*)'  reading CMTSOLUTION file....'
          open(unit=20000,file='CMTSOLUTION',POSITION='REWIND',status='old')
          read(20000,*) junk
          read(20000,*) junk
          read(20000,*) junk
          read(20000,*) junk
          read(20000,*) junk
          read(20000,*) junk
          read(20000,*) junk
          read(20000,*) junk, Mij(isim,1) !Mrr
          read(20000,*) junk, Mij(isim,2) !Mtt
          read(20000,*) junk, Mij(isim,3) !Mpp
          read(20000,*) junk, Mij(isim,4) !Mrt
          read(20000,*) junk, Mij(isim,5) !Mrp
          read(20000,*) junk, Mij(isim,6) !Mtp
          close(20000)
          
          Mij(isim,:) = Mij(isim,:) / 1.e7 ! CMTSOLUTION given in dyn-cm
          
       case('single')
          iinparam_source = 1132
          open(unit=iinparam_source, file='inparam_source', status='old', action='read', iostat=ioerr)
          if (ioerr /= 0) stop 'Check input file ''inparam_source''! Is it still there?'
          
          do
             read(iinparam_source, fmt='(a256)', iostat=ioerr) line
             if (ioerr < 0) exit
             if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
             
             read(line,*) keyword, keyvalue
             
             if (trim(keyword) == 'SOURCE_AMPLITUDE') then
                read(keyvalue,*) amplitude
             end if
          end do
          Mij(isim,:) = 0.0
          select case(src_type(isim,2))
          case('mrr')
             Mij(isim,1) =  amplitude
          case('mtt_p_mpp')
             Mij(isim,2) =  amplitude
             Mij(isim,3) =  amplitude
          case('mtr', 'mrt')
             Mij(isim,4) =  amplitude
          case('mpr', 'mrp')
             Mij(isim,5) =  amplitude
          case('mtp', 'mpt')
             Mij(isim,6) =  amplitude
          case('mtt_m_mpp')
             Mij(isim,2) =  amplitude
             Mij(isim,3) = -amplitude
          case('explosion')
             Mij(isim,1) =  amplitude
             Mij(isim,2) =  amplitude
             Mij(isim,3) =  amplitude
          case default
             write(6,*) 'unknown source type: ', src_type(isim,2)
          end select
          
       case default
          write(6,*)'unknown simulation type!', src_file_type !type(isim,1)
          stop
       end select
  
	write(6,*)'Simu : ',isim
	write(6,*)'Mij = ',Mij(isim,:)
 
    end do
  

    
  end subroutine read_info_simu
!--------------------------------------------------------------------------------

!================================================================================
! Read informations for reconstruction and from coupling
  subroutine read_reconstruction_inputs(isim)
      
    use global_parameters_mod
    use rotation_matrix_mod

    integer(kind=si), intent(in) :: isim
    
    integer(kind=si) :: i
    
    real(kind=dp) :: srclon,  srclat,  srccolat, az_mesh
    real(kind=dp) :: meshlon, meshlat, meshcolat, meshalpha

    integer(kind=si)   :: iin_file=12, ioerr
    character(len=256) :: line
    character(len=256) :: keyword, keyvalue

    write(6, '(A)', advance='no')' Reading lithos.par for recon infos...'
    open(unit=iin_file,file='lithos.par', status='old', action='read', iostat=ioerr)
    if (ioerr /= 0) stop 'Check input file ''lithos.par''! Is it still there?' 
       
    do 
          
       read(iin_file,fmt='(a256)',iostat=ioerr) line
       if (ioerr < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
       read(line,*) keyword, keyvalue
          
       parameters_recon : select case(trim(keyword))
       case('modeling_tool')
          read(keyvalue, *) coup_tool

       case('simulation_nb_buffer')
          read(keyvalue, *) tbuff       

       case('coupling_box_file')
          read(keyvalue, *) input_point_file

       case('coupling_nb_proc_axisem')
          read(keyvalue, *) nbproc
          
       case('coupling_src_latlon')
          read(line, *) keyword, lat_src, lon_src

       case('coupling_axisem_nsim')
          read(keyvalue, *) nsim

       case('nb_partitions')
          read(keyvalue, *) npart

       case('mesh_coordinates')
          read(line, *) keyword, lat_mesh, lon_mesh, az_mesh

       end select parameters_recon
    end do
    close(iin_file)
    write(6,*)'Done!'
    !*** Read recontruction parameter file
!    open(10,file='reconstruction.par',status='old') 
!    read(10,'(a)')coup_tool
!    read(10,'(a)')input_point_file
!    read(10,*)nbproc
!    read(10,*)lat_src,lon_src
!    read(10,*)lat_mesh,lon_mesh,az_mesh
!    read(10,*)nsim
!    close(10)

    !*** Define AxiSEM solutions files
    input_veloc_name(1)  = 'velocityfiel_us'
    input_veloc_name(2)  = 'velocityfiel_up'
    input_veloc_name(3)  = 'velocityfiel_uz'
    input_stress_name(1) = 'stress_Sg11_sol'
    input_stress_name(2) = 'stress_Sg22_sol'
    input_stress_name(3) = 'stress_Sg33_sol'
    input_stress_name(4) = 'stress_Sg12_sol'
    input_stress_name(5) = 'stress_Sg13_sol'
    input_stress_name(6) = 'stress_Sg23_sol'
    
    !*** Define reconstructed output files
    output_veloc_name(1)  = 'velocityoutp_u1'
    output_veloc_name(2)  = 'velocityoutp_u2'
    output_veloc_name(3)  = 'velocityoutp_u3'
    output_stress_name(1) = 'stress_Sg11_out'
    output_stress_name(2) = 'stress_Sg22_out'
    output_stress_name(3) = 'stress_Sg33_out'
    output_stress_name(4) = 'stress_Sg12_out'
    output_stress_name(5) = 'stress_Sg13_out'
    output_stress_name(6) = 'stress_Sg23_out'
    
    !*** Read number of snapshots from axisem simulation
    open(10,file=trim(working_axisem_dir)//trim(simdir(isim))//'/nb_rec_to_read.par')
    read(10,*) ntime
    close(10)

    !*** Alocate energy vector
    if (.not.allocated(energy)) allocate(energy(ntime))
    energy = 0.

    !*** Compute rotation matrices 
    !* For the source
    srclon = lon_src * pi / 180._dp   
    srccolat =  (90._dp - lat_src) * pi / 180._dp 
    call def_rot_matrix(srccolat,srclon,rot_mat,trans_rot_mat)
    
    !* For the mesh
    meshlon = lon_mesh * pi / 180._dp   
    meshlat = lat_mesh * pi / 180._dp
    meshcolat =  (90._dp - lat_mesh) * pi / 180._dp
    meshalpha = az_mesh * pi/ 180._dp 
    select case (coup_tool)
    case ('DG')
       call def_rot_matrix_DG(meshcolat,meshlon,rot_mat_mesh,trans_rot_mat_mesh)
    case ('SEM')
        call def_rot_matrix_SEM(meshlat,meshlon,meshalpha,rot_mat_mesh,trans_rot_mat_mesh)
    case ('FD')
        call def_rot_matrix_SEM(meshlat,meshlon,meshalpha,rot_mat_mesh,trans_rot_mat_mesh)
    end select
    
    !*** Define permuation matrix
    mat(1,1)=0.;    mat(1,2)=1.;    mat(1,3)=0.
    mat(2,1)=1.;    mat(2,2)=0.;    mat(2,3)=0.
    mat(3,1)=0.;    mat(3,2)=0.;    mat(3,3)=1.
    
    tmat(1,1)=mat(1,1);    tmat(1,2)=mat(2,1);    tmat(1,3)=mat(3,1)
    tmat(2,1)=mat(1,2);    tmat(2,2)=mat(2,2);    tmat(2,3)=mat(3,2)
    tmat(3,1)=mat(1,3);    tmat(3,2)=mat(2,3);    tmat(3,3)=mat(3,3)
    

  end subroutine read_reconstruction_inputs
!--------------------------------------------------------------------------------


!================================================================================
! Read inputs_box files
  subroutine  read_input_box_file(isim) !,ipart)

    use global_parameters_mod
    use rotation_matrix_mod

    integer(kind=si), intent(in) :: isim !, ipart
    
    integer(kind=si) :: i, nbreclocc, ipart, istartrec, iendrec
    
    real(kind=dp) :: srclon,  srclat,  srccolat, az_mesh
    real(kind=dp) :: meshlon, meshlat, meshcolat, meshalpha

    character(len=5) :: myfileend

    if (allocated(tab_box_rec)) deallocate(tab_box_rec)
    allocate(tab_box_rec(1:3,1:npart))

    nbrec     = 0
    istartrec = 0
    iendrec   = 0
    do ipart = 1, npart

       write(6,*)'process subdom ',ipart,'over ',npart   

       !*** Read input box files (containing coordinates of the mesh boundaries)
       !*** Find the file
       write(myfileend,'(a,i4.4)')'_',ipart

       !* Number of points to be read
       open(10,file=trim(working_axisem_dir)//trim(simdir(isim))//'/'//trim(input_point_file(1:len(trim(input_point_file))-4))//myfileend//'.txt')
       read(10,*) nbreclocc
       close(10)
       
       if (nbreclocc > 0) then
          istartrec = nbrec + 1 
          nbrec     = nbrec + nbreclocc
          iendrec   = nbrec
       end if

       tab_box_rec(1:3,ipart) = (/ nbreclocc, istartrec, iendrec /)
!       if (nbreclocc > 0) then
!          istartrec = nbrec + 1
!       end if
       write(6,*)'nbrecloc, istart, iend ',tab_box_rec(:,ipart)   

    end do

    !* Allocate relevant arrays
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
    stress_rec=0.
    data_rec=0.
    
    if (allocated(data_rec_all)) deallocate(data_rec_all)
    if (allocated(stress_rec_all)) deallocate(stress_rec_all)
    allocate(data_rec_all(nbrec,3),stress_rec_all(nbrec,6))
    stress_rec_all = 0.
    data_rec_all   = 0.   
    
    if (allocated(f1)) deallocate(f1)
    if (allocated(f2)) deallocate(f2)
    if (allocated(phi)) deallocate(phi)
    allocate(f1(nsim,nbrec),f2(nsim,nbrec),phi(nbrec))
    f1 = 0.
    f2 = 0.

!    nbrec     = 0
!    istartrec = 1
    do ipart = 1, npart
       write(6,*)'process subdom ',ipart,'over ',npart   

       if (tab_box_rec(1,ipart) < 1) then
          cycle
       end if
       
       write(myfileend,'(a,i4.4)')'_',ipart

       !* Number of points to be read
       open(10,file=trim(working_axisem_dir)//trim(simdir(isim))//'/'//trim(input_point_file(1:len(trim(input_point_file))-4))//myfileend//'.txt')
       read(10,*) nbreclocc

!       if (nbreclocc > 0) then
!          nbrec     = nbrec + nbreclocc
!          iendrec   = nbrec
!       end if

       !* Read and computes coordinates for each point
       do i=tab_box_rec(2,ipart),tab_box_rec(3,ipart) !! radius, latitude, longitude
          read(10,*) reciever_geogr(1,i),reciever_geogr(2,i),reciever_geogr(3,i)
          if (reciever_geogr(1,i) > 6371.00000000_dp) then
             reciever_geogr(1,i) = 6371.00000000
          end if
          reciever_sph(1,i) = reciever_geogr(1,i)*1000._dp
          reciever_sph(2,i) = (90._dp - reciever_geogr(2,i)) * pi / 180._dp  
          reciever_sph(3,i) = reciever_geogr(3,i)  * pi / 180._dp

          call rotate_box(reciever_sph(1,i), reciever_sph(2,i), reciever_sph(3,i),trans_rot_mat)

          reciever_cyl(1,i) = reciever_sph(1,i)  * sin( reciever_sph(2,i))
          reciever_cyl(2,i) = reciever_sph(3,i)
          reciever_cyl(3,i) = reciever_sph(1,i)  * cos( reciever_sph(2,i))

          phi(i) = reciever_cyl(2,i) 

       end do
!       if (nbreclocc > 0) then
!          istartrec = nbrec + 1
!       end if
       close(10)
    end do
    

  end subroutine read_input_box_file
!--------------------------------------------------------------------------------


!================================================================================
! Read axisem mesh informations
  subroutine read_mesh(isim)
    
    use global_parameters_mod
    
    integer(kind=si), intent(in) :: isim

    integer(kind=si)   :: iproc, iel, indx_stored
    character(len=250) :: file_to_read
    character(len=4)   :: appmynum
    
    real(kind=dp), allocatable, dimension(:,:,:) :: ssol, zsol
    
    !*** Number of elements for each procs
    allocate(nb_stored(0:nbproc-1))
    nel=0

    !*** Count total number of elements 
    do iproc=0,nbproc-1

       write(file_to_read,'(a10,i5.5,a4)')'parameters',iproc,'.par' 
       open(10,file=trim(working_axisem_dir)//trim(simdir(isim))//'/'//trim(file_to_read))
       read(10,*) nb_stored(iproc),ibeg,iend  
       close(10)
       nel = nel + nb_stored(iproc)
       write(*,*) iproc,nb_stored(iproc),ibeg,iend  
       write(*,*) nel 
       write(*,*) 

    end do
    allocate(scoor(ibeg:iend,ibeg:iend,nel),zcoor(ibeg:iend,ibeg:iend,nel))
    allocate(data_read(ibeg:iend,ibeg:iend,nel))

    !*** Read mesh (path depend on source type)
    indx_stored=1
    do iproc=0,nbproc-1
       if (nb_stored(iproc) > 0) then
          call define_io_appendix(appmynum,iproc)
          allocate(ssol(ibeg:iend,ibeg:iend, nb_stored(iproc)),zsol(ibeg:iend,ibeg:iend, nb_stored(iproc)))
          file_to_read=trim(working_axisem_dir)//trim(simdir(isim))//'Data/mesh_sol_'//appmynum//'.dat'
          open(10,file=trim(file_to_read),FORM="UNFORMATTED")
          read(10) ssol(ibeg:iend,ibeg:iend,:),zsol(ibeg:iend,ibeg:iend,:)
          close(10)
          scoor(ibeg:iend,ibeg:iend,indx_stored:indx_stored+nb_stored(iproc)-1)= ssol(ibeg:iend,ibeg:iend,:)
          zcoor(ibeg:iend,ibeg:iend,indx_stored:indx_stored+nb_stored(iproc)-1)= zsol(ibeg:iend,ibeg:iend,:)
          indx_stored=indx_stored+nb_stored(iproc)
          deallocate(ssol,zsol)
          
       end if
    end do
        
  contains
    
    !================================================================================
    ! Defines the 4 digit character string appended to any data or io file 
    ! related to process myid
    subroutine define_io_appendix(app,iproc)
      
      implicit none
      
      integer(kind=si), intent(in)  :: iproc
      character(len=4), intent(out) :: app
      
      write(app,"(I4.4)") iproc
      
    end subroutine define_io_appendix
    !--------------------------------------------------------------------------------
    
  end subroutine read_mesh
!--------------------------------------------------------------------------------




!################################################################################
!############################## OUTPUTS FILES ###################################
!################################################################################

!================================================================================
! Write output velocity files
  subroutine write_veloc3D(ivx,ivy,ivz)
    
    use global_parameters_mod, only: data_rec_all
    
    integer(kind=si), intent(in) :: ivx, ivy, ivz

    write(ivx) data_rec_all(:,1)
    write(ivy) data_rec_all(:,2)
    write(ivz) data_rec_all(:,3)

  end subroutine write_veloc3D
!--------------------------------------------------------------------------------


!================================================================================
! Write output stress files 
  subroutine write_stress3D(isxx,isyy,iszz,isxy,isxz,isyz)

    use global_parameters_mod, only: stress_rec_all

    integer(kind=si), intent(in) :: isxx, isyy, iszz, isxy, isxz, isyz
 
    write(isxx) stress_rec_all(:,1)
    write(isyy) stress_rec_all(:,2)
    write(iszz) stress_rec_all(:,3)
    write(isxy) stress_rec_all(:,4)
    write(isxz) stress_rec_all(:,5)
    write(isyz) stress_rec_all(:,6)
    
  end subroutine write_stress3D
!--------------------------------------------------------------------------------


end module inputs_outputs_mod
