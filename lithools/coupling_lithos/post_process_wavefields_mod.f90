module post_process_wavefields_mod

  use precision_mod
  use constants_mod

  implicit none

contains

!================================================================================
! Routine reconstructing the 3D velocity wavefield
  subroutine reconstruct_velocity(isim)

    use mpi_mod
    use global_parameters_mod
    use rotation_matrix_mod
    use inputs_outputs_mod
    
    integer(kind=si), intent(in) :: isim

    integer(kind=si) :: itime, iproc, ifield, i, ivx, ivy, ivz    
    integer(kind=si), dimension(3) :: indx_stored
    integer(kind=si), dimension(:,:), allocatable :: iunit

    real(kind=sp), allocatable, dimension(:,:,:) :: data_to_read

    character(len=4)   :: appmynum
    character(len=256) :: fichier

    allocate(iunit(0:nbproc-1,3))
        

    !*** Unit file
    i=150
    do ifield=1,3
       do iproc=0, nbproc-1
          i=i+1
          iunit(iproc,ifield)=i
       end do
    end do
    
    i=i+1
    ivx=i
    i=i+1
    ivy=i
    i=i+1
    ivz=i
    
    !*** Compute prefactor for wavefield reconstruction
    write(*,*) src_type
    call compute_prefactor(src_type(isim,1),src_type(isim,2))
    
    !*** Read AxiSEM solutions
    if (myid == 0) then
       write(fichier,'(a6,a15)') '/Data/',output_veloc_name(1)
       open(ivx,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_veloc_name(2)
       open(ivy,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_veloc_name(3)
       open(ivz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       
       write(*,*) 'nbrec to write', nbrec
       write(*,*) 'nbrec to read ',sum(nb_stored(:))
       write(*,*) 'nt to write ', ntime
        
       write(ivx) nbrec,ntime
       write(ivy) nbrec,ntime
       write(ivz) nbrec,ntime
       
       data_rec=0.
       
       !*** For each velocity component
       do ifield=1,3 
          write(*,*) ifield
          if (trim(src_type(isim,1))=='monopole' .and. ifield==2) then
             write(*,*) 'monopole => not up'
             cycle
          end if
          
          !* Open files
          do iproc=0, nbproc-1
             call define_io_appendix(appmynum,iproc)
             write(fichier,'(a6,a15,a1)') '/Data/',input_veloc_name(ifield),'_'
             open(unit=iunit(iproc,ifield), file=trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier) &
                  //appmynum//'.bindat', &
                  FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
          end do
       end do
    end if
    
    !*** Read those files
    do itime=1,ntime
       data_rec=0.
       indx_stored=1
       !** us comp
       ifield=1
       if (myid ==0) then 
          do iproc=0, nbproc-1
             if (nb_stored(iproc) > 0) then
                
                allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                     = data_to_read(ibeg:iend,ibeg:iend,:)
                indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                
                deallocate(data_to_read)
             end if
          end do
       end if
       call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi) 
       !* Interpolate us comp
       call interpol_field(ifield)
 
       if (.not.(trim(src_type(isim,1))=='monopole')) then
          ifield=2
          if (myid == 0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                end if
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi)
          call interpol_field(ifield)
       end if

       !** uz comp
       ifield=3
       if (myid ==0) then
          do iproc=0, nbproc-1
             if (nb_stored(iproc) > 0) then
                allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                     = data_to_read(ibeg:iend,ibeg:iend,:)
                indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                deallocate(data_to_read)
             end if
          end do
       end if
       call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi)
       call interpol_field(ifield)           

       !*** Compute cylindrical coordinates
       call compute_3D_cyl

       !*** Perform rotations
       call rotate2cartesian_with_source_in_pole
       call rotate_back_source                    ! global earth cartesian 
       call rotate_back_to_local_cart             ! local cartesian
       call reduce_mpi_veloc
       
       !*** Write wavefield
       if (myid == 0) call write_veloc3D(ivx,ivy,ivz)

    end do ! time step

    if (myid ==0) then
       !*** Close files 
       do ifield=1,3
          do iproc=0, nbproc-1
             close(iunit(iproc,ifield))
          end do
       end do
       close(ivx)
       close(ivy)
       close(ivz)
    end if
    
  end subroutine reconstruct_velocity
!--------------------------------------------------------------------------------


!================================================================================
! Routine reconstructing the 3D velocity wavefield
  subroutine read_stress_field_and_interpol(isim)
  
    use mpi_mod
    use global_parameters_mod
    use rotation_matrix_mod
    use inputs_outputs_mod
    
    integer(kind=si), intent(in) :: isim

    integer(kind=si) :: itime, iproc, ifield, i
    integer(kind=si) :: isxx, isyy, iszz, isxy, isxz, isyz
    integer(kind=si), dimension(6) :: indx_stored
    integer(kind=si), allocatable, dimension(:,:) :: iunit
    
    real(kind=sp), allocatable, dimension(:,:,:)  :: data_to_read
    
    character(len=4)   :: appmynum
    character(len=256) :: fichier
        
    allocate(iunit(0:nbproc-1,6))
        
    ! unit file
    i=150
    do ifield=1,6
       do iproc=0, nbproc-1
          i=i+1
          iunit(iproc,ifield)=i
       end do
    end do
    
    i=i+1
    isxx=i
    i=i+1
    isyy=i
    i=i+1
    iszz=i
    i=i+1
    isxy=i
    i=i+1
    isxz=i
    i=i+1
    isyz=i
    
    call compute_prefactor(src_type(isim,1),src_type(isim,2))
    if (myid == 0) then  
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(1)
       open(isxx,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(2)
       open(isyy,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(3)
       open(iszz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(4)
       open(isxy,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(5)
       open(isxz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(6)
       open(isyz,file= trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier), FORM="UNFORMATTED")
       
       write(isxx) nbrec,ntime
       write(isyy) nbrec,ntime
       write(iszz) nbrec,ntime
       write(isxy) nbrec,ntime
       write(isxz) nbrec,ntime
       write(isyz) nbrec,ntime
       
       data_rec=0.
       do ifield=1,6 
          write(*,*) ifield
          if (trim(src_type(isim,1))=='monopole' .and. (ifield==4 .or. ifield==6)) then
             write(*,*) 'monopole => not up'
             cycle
          end if
          ! open files
          do iproc=0, nbproc-1
             call define_io_appendix(appmynum,iproc)
             write(fichier,'(a6,a15,a1)') '/Data/',input_stress_name(ifield),'_'
             open(unit=iunit(iproc,ifield), file=trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier) &
                  //appmynum//'.bindat', &
                  FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
          end do
       end do
    end if
    
    ! read files
    do itime=1,ntime
       stress_rec=0. 
       indx_stored=1
       !! s11
       ifield=1
       if (myid == 0) then
          do iproc=0, nbproc-1
             if (nb_stored(iproc) > 0) then
                
                allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                     = data_to_read(ibeg:iend,ibeg:iend,:)
                
                indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                deallocate(data_to_read)
                
             end if
          end do
       end if
       call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi)
       call interpol_stress(ifield)
       
       
       !! s22
       ifield=2
       if (myid==0) then
          do iproc=0, nbproc-1
             if (nb_stored(iproc) > 0) then
                
                allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                     = data_to_read(ibeg:iend,ibeg:iend,:)
                indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                deallocate(data_to_read)
                
             end if
          end do
       end if
       call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi)
       call interpol_stress(ifield)
       
       
       !! s33
       ifield=3
       if (myid==0) then
          do iproc=0, nbproc-1
             if (nb_stored(iproc) > 0) then
                
                allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                     = data_to_read(ibeg:iend,ibeg:iend,:)
                indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                deallocate(data_to_read)
                
             end if
             
          end do
       end if
       call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi)
       call interpol_stress(ifield)
       
       !! s12
       if (.not.(trim(src_type(isim,1))=='monopole')) then
          ifield=4
          if (myid==0) then 
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                   
                end if
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi)
          call interpol_stress(ifield)
       end if
       
       
       !! s13
       ifield=5
       if (myid == 0) then
          do iproc=0, nbproc-1
             if (nb_stored(iproc) > 0) then
                
                allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                     = data_to_read(ibeg:iend,ibeg:iend,:)
                indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                deallocate(data_to_read)
                
             end if
             
          end do
       end if
       call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi) 
       call interpol_stress(ifield)
       
       
       !! s23
       if (.not.(trim(src_type(isim,1))=='monopole')) then
          ifield=6
          if (myid ==0) then 
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read(iunit(iproc,ifield))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                   
                end if
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPISP,0,MPI_COMM_WORLD,ierr_mpi)
          call interpol_stress(ifield)
       end if
       
       !*** Compute cylindrical coordinates
       call compute_stress_3D_cyl

       !*** Perform rotations
       call rotate2cartesian_with_source_in_pole
       call rotate_back_source                    ! global earth cartesian 
       call rotate_back_to_local_cart             ! local cartesian
       call reduce_mpi_stress
       
       !*** Write wavefield
       if (myid == 0) call  write_stress3D(isxx,isyy,iszz,isxy,isxz,isyz)
       
    end do
    
    ! close files 
    do ifield=1,6
       do iproc=0, nbproc-1
          close(iunit(iproc,ifield))
       end do
    end do
    
    close(isxx)
    close(isyy)
    close(iszz)
    close(isxy)
    close(isxz)
    close(isyz)
    
  end subroutine read_stress_field_and_interpol
!--------------------------------------------------------------------------------  

!================================================================================
! Interpolate velocity wavefield
  subroutine interpol_field(ifield)
        
    use global_parameters_mod
    use interp_mesh_mod

    integer(kind=si), intent(in) :: ifield

    integer(kind=si) :: irec, iel
    real(kind=cp) :: xi,eta,interp_value
    real(kind=cp), dimension(NGLLX,NGLLY) :: field
    
    data_rec(:,ifield) = 0.
    
    do irec = irecmin, irecmax
       
       xi=xi_rec(irec)
       eta=eta_rec(irec)
       iel=rec2elm(irec)
       field(:,:)=data_read(:,:,iel)
       
       call interpole_field(xi,eta,field,interp_value) 
       
       data_rec(irec,ifield)=interp_value
       
    end do
    
  end subroutine interpol_field
!--------------------------------------------------------------------------------

!================================================================================
! Interpolate stress wavefields
  subroutine interpol_stress(ifield)
        
    use global_parameters_mod
    use interp_mesh_mod
    
    integer(kind=si), intent(in) :: ifield

    integer(kind=si) :: irec, iel
    real(kind=cp) :: xi,eta,interp_value
    real(kind=cp), dimension(NGLLX,NGLLY) :: field
    
    stress_rec(:,ifield)=0.
    do irec = irecmin, irecmax
       
       xi=xi_rec(irec)
       eta=eta_rec(irec)
       iel=rec2elm(irec)
       field(:,:)=data_read(:,:,iel)
       
       call interpole_field(xi,eta,field,interp_value) 
       
       stress_rec(irec,ifield)=interp_value

    end do
    
  end subroutine interpol_stress
!--------------------------------------------------------------------------------

!================================================================================
! Routine computing prefactor for wavefield extrapolation
  subroutine compute_prefactor(src_type,Mcomp)
    
    use global_parameters_mod, only: f1, f2, phi, nbrec
    
    character(len=10), intent(in)  :: src_type
    character(len=10), intent(out) :: Mcomp
    integer(kind=si) :: irec

    select case (trim(src_type))
    case('monopole')
       f1=1.
       f2=0.
    case('dipole')
       select case (trim(Mcomp))
       case('mtr')
          do irec=1,nbrec
             f1(irec)=cos(phi(irec))
             f2(irec)=-sin(phi(irec))
          end do
       case ('thetaforce')
          do irec=1,nbrec
             f1(irec)=cos(phi(irec))
             f2(irec)=-sin(phi(irec))
          end do
       case('mpr')
          do irec=1,nbrec
             f1(irec)=sin(phi(irec))
             f2(irec)=cos(phi(irec))
          end do
       case('phiforce')
          do irec=1,nbrec
             f1(irec)=sin(phi(irec))
             f2(irec)=cos(phi(irec))
          end do
       end select
    case('quadpole')
       select case (trim(Mcomp))
       case ('mtt_m_mpp')
          do irec=1,nbrec
             f1(irec)=cos(2.*phi(irec))
             f2(irec)=-sin(2.*phi(irec))
          end do
       case('mtp')
          do irec=1,nbrec
             f1(irec)=sin(2*phi(irec))
             f2(irec)=cos(2*phi(irec))
          end do
       end select
       
    end select
    
  end subroutine compute_prefactor
!--------------------------------------------------------------------------------


!================================================================================
! Compute cylinfrical coordinates 3D (vel)
  subroutine compute_3D_cyl

    use global_parameters_mod

    integer(kind=si) :: irec

    do irec=irecmin,irecmax
       data_rec(irec,1)=f1(irec)*data_rec(irec,1)
       data_rec(irec,2)=f2(irec)*data_rec(irec,2)
       data_rec(irec,3)=f1(irec)*data_rec(irec,3)
    end do

  end subroutine compute_3D_cyl
!--------------------------------------------------------------------------------

!================================================================================
! Compute cylindrical coordinates 3D (stress)
  subroutine compute_stress_3D_cyl
    
    use global_parameters_mod

    integer(kind=si) :: irec

    do irec=irecmin,irecmax
       stress_rec(irec,1)=f1(irec)*stress_rec(irec,1)
       stress_rec(irec,2)=f1(irec)*stress_rec(irec,2)
       stress_rec(irec,3)=f1(irec)*stress_rec(irec,3)
       stress_rec(irec,4)=f2(irec)*stress_rec(irec,4)
       stress_rec(irec,5)=f1(irec)*stress_rec(irec,5)
       stress_rec(irec,6)=f2(irec)*stress_rec(irec,6)
    end do
    
  end subroutine compute_stress_3D_cyl
!--------------------------------------------------------------------------------

end module post_process_wavefields_mod
