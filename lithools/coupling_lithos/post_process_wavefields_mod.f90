module post_process_wavefields_mod

  use precision_mod
  use constants_mod
!  use mpi_mod
  implicit none

contains

!================================================================================
! Routine reconstructing the 3D velocity wavefield
  subroutine reconstruct_velocity(ipart)  !(isim)

    use global_parameters_mod
    use rotation_matrix_mod
    use inputs_outputs_mod
    use mpi_mod

    integer(kind=si) :: isim, ipart    
    !, intent(in) :: isim

    integer(kind=si) :: itime, iproc, ifield, i, ivx, ivy, ivz    
    integer(kind=si), dimension(3) :: indx_stored
    integer(kind=si), dimension(:,:,:), allocatable :: iunit

    real(kind=sp), allocatable, dimension(:,:,:) :: data_to_read

    character(len=4)   :: appmynum
    character(len=5)   :: myfileend
    character(len=256) :: fichier

    allocate(iunit(0:nbproc-1,3,nsim))
    
    !*** Unit file
    i=150
    do isim = 1,nsim    
       do ifield=1,3
          do iproc=0, nbproc-1
             i=i+1
             iunit(iproc,ifield,isim)=i
          end do
       end do
    end do

    !*** For output
    i=i+1
    ivx=i
    i=i+1
    ivy=i
    i=i+1
    ivz=i

    !*** Compute prefactor for wavefield reconstruction
    do isim = 1, nsim
       call compute_prefactor(src_type(isim,1),src_type(isim,2),isim)
    end do

    !*** Open AxiSEM solutions files and outputs
    if (myid == 0) then

       !*** Outputs velocity
       write(myfileend,'(a,i4.4)')'_',ipart

       write(fichier,'(a15)') output_veloc_name(1)
       open(ivx,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       write(fichier,'(a15)') output_veloc_name(2)
       open(ivy,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       write(fichier,'(a15)') output_veloc_name(3)
       open(ivz,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       
       write(*,*) 'nbrec to write', nbrec
       write(*,*) 'nbrec to read ',sum(nb_stored(:))
       write(*,*) 'nt to write ', ntime
       
       write(ivx) nbrec,ntime
       write(ivy) nbrec,ntime
       write(ivz) nbrec,ntime
       
       data_rec=0.
       
       !*** For each velocity component input
       do isim = 1,nsim
          do ifield=1,3 
             write(*,*) isim, ifield
             if (trim(src_type(isim,1))=='monopole' .and. ifield==2) then
                write(*,*) 'monopole => not up'
                cycle
             end if
             
             !* Open files
             do iproc=0, nbproc-1
                call define_io_appendix(appmynum,iproc)
                write(fichier,'(a6,a15,a1)') '/Data/',input_veloc_name(ifield),'_'
                open(unit= iunit(iproc,ifield,isim), file=trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier) &
                     //appmynum//'.bindat', &
                     FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
             end do
          end do
       end do
    end if
    
        
    if (myid==0) write(6,*)'Percentage : '
    
    !*** Read those files at each time step
    do itime=1,ntime

       data_rec_all = 0.   ! sum over source simulations
       
       do isim=1,nsim
          
          data_rec=0.
          indx_stored=1
          
          !** us comp
          ifield=1
          if (myid ==0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   
                   deallocate(data_to_read)
                end if
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi) 
          call interpol_field(ifield)
          

          !** up comp
          if (.not.(trim(src_type(isim,1))=='monopole')) then
             ifield=2
             if (myid == 0) then
                do iproc=0, nbproc-1
                   if (nb_stored(iproc) > 0) then
                      allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                      read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                      data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                           = data_to_read(ibeg:iend,ibeg:iend,:)
                      indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                      deallocate(data_to_read)
                   end if
                end do
             end if
             call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
             call interpol_field(ifield)
          end if
          
          !** uz comp
          ifield=3
          if (myid ==0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                end if
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
          call interpol_field(ifield)           
       
          !*** Compute cylindrical components (prefactors) and sum 
          call compute_3D_cyl(isim)
          !data_rec_all = data_rec_all + data_rec
          
       end do !isim

       !*** Perform rotations
       call rotate2cartesian_with_source_in_pole
       call rotate_back_source                    ! global earth cartesian 
       call rotate_back_to_local_cart             ! local cartesian
       call reduce_mpi_veloc
       
       !*** Write wavefield
       if (myid == 0) then
          !* Compute energy
          energy(itime) = energy(itime) + sum(sqrt(data_rec_all(:,1)**2 + data_rec_all(:,2)**2 + data_rec_all(:,3)**2))
          
          !* Write to disk
          call write_veloc3D(ivx,ivy,ivz)
          write(6,*)itime,ntime                  
       end if
    end do ! time step
    

    !*** Close files
    if (myid ==0) then
       !*** Close files read
       do isim = 1, nsim
          do ifield=1,3
             do iproc=0, nbproc-1
                close( iunit(iproc,ifield,isim))
             end do
          end do
       end do
       !*** Close output files
       close(ivx)
       close(ivy)
       close(ivz)
    end if
    
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
    


  end subroutine reconstruct_velocity
!--------------------------------------------------------------------------------


!================================================================================
! Routine reconstructing the 3D velocity wavefield
  subroutine reconstruct_stress(ipart)   !(isim)
  
    use mpi_mod
    use global_parameters_mod
    use rotation_matrix_mod
    use inputs_outputs_mod
    
    integer(kind=si) :: isim, ipart   
    !, intent(in) :: isim

    integer(kind=si) :: itime, iproc, ifield, i
    integer(kind=si) :: isxx, isyy, iszz, isxy, isxz, isyz
    integer(kind=si), dimension(6) :: indx_stored
    integer(kind=si), allocatable, dimension(:,:,:) :: iunit
    
    real(kind=sp), allocatable, dimension(:,:,:)  :: data_to_read
    
    character(len=4)   :: appmynum
    character(len=5)   :: myfileend
    character(len=256) :: fichier
        
    allocate(iunit(0:nbproc-1,6,nsim))
        
    !*** Unit file to read
    i=150
    do isim = 1, nsim
       do ifield=1,6
          do iproc=0, nbproc-1
             i=i+1
             iunit(iproc,ifield,isim)=i
          end do
       end do
    end do

    !*** For outputs
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


    !*** Compute prefactor for wavefield reconstruction
    do isim = 1, nsim
       call compute_prefactor(src_type(isim,1),src_type(isim,2),isim)
    end do

    !*** Open AxiSEM solutions files and outputs
    if (myid == 0) then  

       !*** Outputs stress
       write(myfileend,'(a,i4.4)')'_',ipart
       write(fichier,'(a15)') output_stress_name(1)
       open(isxx,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(2)
       open(isyy,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(3)
       open(iszz,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(4)
       open(isxy,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(5)
       open(isxz,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(6)
       open(isyz,file= trim(working_axisem_dir)//trim(fichier)//myfileend, FORM="UNFORMATTED")
       
       write(isxx) nbrec,ntime
       write(isyy) nbrec,ntime
       write(iszz) nbrec,ntime
       write(isxy) nbrec,ntime
       write(isxz) nbrec,ntime
       write(isyz) nbrec,ntime
       
       stress_rec=0.

       !*** For each inpuut stress components
       do isim = 1, nsim
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
                open(unit= iunit(iproc,ifield,isim), file=trim(working_axisem_dir)//trim(simdir(isim))//trim(fichier) &
                     //appmynum//'.bindat', &
                     FORM="UNFORMATTED", STATUS="UNKNOWN", POSITION="REWIND")
             end do
          end do
       end do
    end if
    
    
    if (myid==0) write(6,*)'Percentage : '


    !*** Read those files at each time step
    do itime=1,ntime

       stress_rec_all = 0.    ! sum over source simulations
       
       do isim = 1, nsim
          

          stress_rec=0. 
          indx_stored=1


          !*** s11 comp
          ifield=1
          if (myid == 0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                   
                end if
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
          call interpol_stress(ifield)
       
          
          !*** s22 comp
          ifield=2
          if (myid==0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                   
                end if
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
          call interpol_stress(ifield)
          
       
          !*** s33 comp
          ifield=3
          if (myid==0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                   
                end if
                
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
          call interpol_stress(ifield)
       
          !*** s12 comp
          if (.not.(trim(src_type(isim,1))=='monopole')) then
             ifield=4
             if (myid==0) then 
                do iproc=0, nbproc-1
                   if (nb_stored(iproc) > 0) then
                      
                      allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                      read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                      data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                           = data_to_read(ibeg:iend,ibeg:iend,:)
                      indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                      deallocate(data_to_read)
                      
                   end if
                end do
             end if
             call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
             call interpol_stress(ifield)
          end if
       
       
          !*** s13 comp
          ifield=5
          if (myid == 0) then
             do iproc=0, nbproc-1
                if (nb_stored(iproc) > 0) then
                   
                   allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                   read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                   data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                        = data_to_read(ibeg:iend,ibeg:iend,:)
                   indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                   deallocate(data_to_read)
                   
                end if
                
             end do
          end if
          call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi) 
          call interpol_stress(ifield)
          
       
          !*** s23 comp
          if (.not.(trim(src_type(isim,1))=='monopole')) then
             ifield=6
             if (myid ==0) then 
                do iproc=0, nbproc-1
                   if (nb_stored(iproc) > 0) then
                      
                      allocate(data_to_read(ibeg:iend,ibeg:iend, nb_stored(iproc)))
                      read( iunit(iproc,ifield,isim))  data_to_read(ibeg:iend,ibeg:iend,:)
                      data_read(ibeg:iend, ibeg:iend, indx_stored(ifield):indx_stored(ifield)+nb_stored(iproc)-1) &
                           = data_to_read(ibeg:iend,ibeg:iend,:)
                      indx_stored(ifield)=indx_stored(ifield)+nb_stored(iproc)
                      deallocate(data_to_read)
                      
                   end if
                end do
             end if
             call mpi_bcast(data_read,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_REAL,0,MPI_COMM_WORLD,ierr_mpi)
             call interpol_stress(ifield)
          end if
       
       
          !*** compute cylindrical componenets (with prefactors) and sum
          call compute_stress_3D_cyl(isim)
          !stress_rec_all = stress_rec_all + stress_rec
          
       end do ! isim
       
       !*** Perform rotations
       call rotate2cartesian_with_source_in_pole_stress
       call rotate_back_source_stress                    ! global earth cartesian 
       call rotate_back_to_local_cart_stress             ! local cartesian
       call reduce_mpi_stress
       
       !*** Write wavefield
       if (myid == 0) then
          call  write_stress3D(isxx,isyy,iszz,isxy,isxz,isyz)
          write(6,*)itime,ntime                   
       end if
       
    end do
    
    !*** Close files 
    if (myid == 0) then
       !*** Close inputs files
       do isim = 1, nsim
          do ifield=1,6
             do iproc=0, nbproc-1
                close( iunit(iproc,ifield,isim))
             end do
          end do
       end do
       !*** Close output files
       close(isxx)
       close(isyy)
       close(iszz)
       close(isxy)
       close(isxz)
       close(isyz)
    end if

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
    
  end subroutine reconstruct_stress
!--------------------------------------------------------------------------------  

!================================================================================
! Interpolate velocity wavefield
  subroutine interpol_field(ifield)
        
    use global_parameters_mod
    use interp_mesh_mod

    integer(kind=si), intent(in) :: ifield

    integer(kind=si) :: irec, iel
    real(kind=dp) :: xi,eta,interp_value
    real(kind=dp), dimension(NGLLX,NGLLY) :: field
    
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
    real(kind=dp) :: xi,eta,interp_value
    real(kind=dp), dimension(NGLLX,NGLLY) :: field
    
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
  subroutine compute_prefactor(src_type,Mcomp,isim)
    
    
    use mpi_mod, only: myid, ierr_mpi, MPI_COMM_WORLD 
    use global_parameters_mod, only: Mij, phi, nbrec, nsim, magnitude, mij_scale, f1, f2
    
    implicit none
    
    integer(kind=si), intent(in) :: isim
    character(len=10), intent(in)  :: src_type
    character(len=10), intent(in) :: Mcomp
    integer(kind=si) :: irec

    !*** Scale tensor
    Mij_scale(:) = 0.

    Mij_scale(1) = Mij(isim,1) / magnitude(isim) 
    Mij_scale(2) = Mij(isim,2) / magnitude(isim) 
    Mij_scale(3) = Mij(isim,3) / magnitude(isim) 
    Mij_scale(4) = Mij(isim,4) / magnitude(isim) 
    Mij_scale(5) = Mij(isim,5) / magnitude(isim) 
    Mij_scale(6) = Mij(isim,6) / magnitude(isim) 

    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
    if (myid == 0) then
        print *,'src_type, Mcomp, isim '
        print *,src_type, Mcomp, isim

       write(6,*)'Mij scaled: on proc 0'
       write(6,*) Mij_scale
       write(6,*)'-----------------'
       write(6,*) Mij(isim,:)
       write(6,*)'-----------------'
       write(6,*) magnitude(isim)
       write(6,*)'-----------------'
    end if
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
    if (myid == 1) then
       write(6,*)'Mij scaled: on proc 1'
       write(6,*) Mij_scale
       write(6,*)'-----------------'
       write(6,*) Mij(isim,:)
       write(6,*)'-----------------'
       write(6,*) magnitude(isim)
       write(6,*)'-----------------'
    end if
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)

    !*** Prefactors
    select case (trim(src_type))
    case('monopole')

       if(myid==0) print *,'i fond monopole'

       select case(Mcomp)
       case('mrr')
          f1(isim,:) = Mij_scale(1)
          f2(isim,:) = 0.
       case('mtt_p_mpp')
          f1(isim,:) = Mij_scale(2) + Mij_scale(3)
          f2(isim,:) = 0.
       case('explosion')
          f1(isim,:) = (Mij_scale(1) + Mij_scale(2) + Mij_scale(3)) / 3.
          f2(isim,:) = 0.
       end select

    case('dipole')

       if(myid==0) print *,'i fond diple'

       select case (trim(Mcomp))
       case('mtr','mrt')
          do irec=1,nbrec
             f1(isim,irec) =   Mij_scale(4) * cos(phi(irec)) + Mij_scale(5) * sin(phi(irec))   !   Mij_scale(4) * cos(phi(irec))
             f2(isim,irec) = - Mij_scale(4) * sin(phi(irec)) + Mij_scale(5) * cos(phi(irec))  ! - Mij_scale(4) * sin(phi(irec))
          end do
       case ('thetaforce')   !!! TO verify
          do irec=1,nbrec
             f1(isim,irec) =   cos(phi(irec))
             f2(isim,irec) = - sin(phi(irec))
          end do
       case('mpr','mrp')
          do irec=1,nbrec
             f1(isim,irec) =   Mij_scale(4) * cos(phi(irec)) + Mij_scale(5) * sin(phi(irec))   !   Mij_scale(4) * cos(phi(irec))
             f2(isim,irec) = - Mij_scale(4) * sin(phi(irec)) + Mij_scale(5) * cos(phi(irec))  ! - Mij_scale(4) * sin(phi(irec))
!             f1(isim,irec) = Mij_scale(5) * sin(phi(irec))
!             f2(isim,irec) = Mij_scale(5) * cos(phi(irec))
          end do
       case('phiforce')      !!! TO verify
          do irec=1,nbrec
             f1(isim,irec) = sin(phi(irec))
             f2(isim,irec) = cos(phi(irec))
          end do
       end select

    case('quadpole')

       if(myid==0) print *,'i fond quadpole'

       select case (trim(Mcomp))
       case ('mtt_m_mpp')
          do irec=1,nbrec
             f1(isim,irec) = (Mij_scale(2) - Mij_scale(3)) * cos(2.*phi(irec)) + 2. * Mij_scale(6) * sin(2.*phi(irec)) !!! Check -1/2 
             f2(isim,irec) = (Mij_scale(3) - Mij_scale(2)) * sin(2.*phi(irec)) + 2. * Mij_scale(6) * cos(2.*phi(irec)) !!! Check -1/2
!             f1(isim,irec)=   (Mij_scale(2) - Mij_scale(3)) * cos(2.*phi(irec)) !!! Check -1/2 
!             f2(isim,irec)= - (Mij_scale(3) - Mij_scale(2)) * sin(2.*phi(irec)) !!! Check -1/2
          end do
       case('mtp','mpt')
          do irec=1,nbrec
!             f1(isim,irec) =  2. * Mij_scale(6) * sin(2*phi(irec))    !!! Check -1/2  .. 2.*
!             f2(isim,irec) =  2. * Mij_scale(6) * cos(2*phi(irec))    !!! Check -1/2  .. 2./*
             f1(isim,irec) = (Mij_scale(2) - Mij_scale(3)) * cos(2.*phi(irec)) + 2. * Mij_scale(6) * sin(2.*phi(irec)) !!! Check -1/2 
             f2(isim,irec) = (Mij_scale(3) - Mij_scale(2)) * sin(2.*phi(irec)) + 2. * Mij_scale(6) * cos(2.*phi(irec)) !!! Check -1/2
          end do
       end select
       
    end select


    ! check for debug
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
    if (myid == 0) then
       write(6,*) 'F1,F2 prefactor max proc0:', maxval(f1(isim,:)), maxval(f2(isim,:))
       write(6,*) 'F1,F2 prefactor min proc0:', minval(f1(isim,:)), minval(f2(isim,:))
       write(6,*) 'Phi min and max, proc 0 ', minval(phi(:)), maxval(phi(:))
       write(6,*) 'Nbrec: ',nbrec
    end if
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
    if (myid == 1) then
       write(6,*) 'F1,F2 prefactor max proc1:', maxval(f1(isim,:)), maxval(f2(isim,:))
       write(6,*) 'F1,F2 prefactor min proc1:', minval(f1(isim,:)), minval(f2(isim,:))
       write(6,*) 'Phi min and max, proc 0 ', minval(phi(:)), maxval(phi(:))
       write(6,*) 'Nbrec: ',nbrec
    end if
    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)


!!$    !*** For moment tensor
!!$!    if (.not.allocated(mij_prefact)) allocate(mij_prefact(nbrec,nsim,3))
!!$!    if (.not.allocated(Mij_scale)) allocate(Mij_scale(nsim,3))
!!$    Mij_scale = 0.
!!$    Mij_scale(isim,1) = Mij(isim,1) / magnitude(isim)
!!$    Mij_scale(isim,2) = Mij(isim,2) / magnitude(isim)
!!$    Mij_scale(isim,3) = Mij(isim,3) / magnitude(isim)
!!$    Mij_scale(isim,4) = Mij(isim,4) / magnitude(isim)
!!$    Mij_scale(isim,5) = Mij(isim,5) / magnitude(isim)
!!$    Mij_scale(isim,6) = Mij(isim,6) / magnitude(isim)
!!$    mij_prefact = 0.
!!$    
!!$    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
!!$    if (myid == 0) then
!!$       write(6,*)'Mij scaled: on proc 0'
!!$       write(6,*) Mij_scale
!!$       write(6,*)'-----------------'
!!$       write(6,*) Mij
!!$       write(6,*)'-----------------'
!!$       write(6,*) magnitude
!!$       write(6,*)'-----------------'
!!$    end if
!!$    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
!!$    if (myid == 1) then
!!$       write(6,*)'Mij scaled: on proc 1'
!!$       write(6,*) Mij_scale
!!$       write(6,*)'-----------------'
!!$       write(6,*) Mij
!!$       write(6,*)'-----------------'
!!$       write(6,*) magnitude
!!$       write(6,*)'-----------------'
!!$    end if
!!$    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
!!$
!!$    select case(Mcomp)
!!$    case('mrr')
!!$       do irec=1,nbrec
!!$          mij_prefact(irec,isim,:) = Mij_scale(isim,1)
!!$          mij_prefact(irec,isim,2) = 0.
!!$       end do
!!$       write(6,*) isim, 'Simulation is mrr, prefact:', &
!!$            mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
!!$            mij_prefact(1,isim,3)
!!$    case('mtt_p_mpp')
!!$       do irec=1,nbrec
!!$          mij_prefact(irec,isim,:) = Mij_scale(isim,2) + Mij_scale(isim,3)
!!$          mij_prefact(irec,isim,2) = 0.
!!$       end do
!!$       write(6,*) isim, 'Simulation is mpp, prefact:', &
!!$            mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
!!$            mij_prefact(1,isim,3)
!!$       
!!$    case('mtr', 'mrt', 'mpr', 'mrp')
!!$       do irec=1,nbrec
!!$          mij_prefact(irec,isim,1) =   Mij_scale(isim,4) * cos(phi(irec)) &
!!$               + Mij_scale(isim,5) * sin(phi(irec))
!!$          mij_prefact(irec,isim,2) = - Mij_scale(isim,4) * sin(phi(irec)) &
!!$               + Mij_scale(isim,5) * cos(phi(irec))
!!$          mij_prefact(irec,isim,3) =   Mij_scale(isim,4) * cos(phi(irec)) &
!!$               + Mij_scale(isim,5) * sin(phi(irec))
!!$       end do
!!$       write(6,*) isim, 'Simulation is mtr, prefact:', &
!!$            mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
!!$            mij_prefact(1,isim,3)
!!$       
!!$    case('mtp', 'mpt', 'mtt_m_mpp')
!!$       do irec=1,nbrec
!!$          mij_prefact(irec,isim,1) = (Mij_scale(isim,2) - Mij_scale(isim,3)) * cos(2. * phi(irec))  &
!!$               + 2. * Mij_scale(isim,6)  * sin(2. * phi(irec))
!!$          mij_prefact(irec,isim,2) = (Mij_scale(isim,3) - Mij_scale(isim,2)) * sin(2. * phi(irec)) &
!!$               + 2. * Mij_scale(isim,6)  * cos(2. * phi(irec))
!!$          mij_prefact(irec,isim,3) = (Mij_scale(isim,2) - Mij_scale(isim,3)) * cos(2. * phi(irec))  &
!!$               + 2. * Mij_scale(isim,6)  * sin(2. * phi(irec))
!!$       end do
!!$       write(6,*) isim, 'Simulation is mtp, prefact:', &
!!$            mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
!!$            mij_prefact(1,isim,3)
!!$       
!!$    case('explosion')
!!$       do irec=1,nbrec
!!$          mij_prefact(irec,isim,:) = (Mij_scale(isim,1) + Mij_scale(isim,2) + Mij_scale(isim,3)) / 3.
!!$       end do
!!$       write(6,*) isim, 'Simulation is explosion, prefact:', &
!!$            mij_prefact(1,isim,1), mij_prefact(1,isim,2), &
!!$            mij_prefact(1,isim,3)
!!$       
!!$    case default
!!$       write(6,*) 'unknown source type ', Mcomp
!!$       stop
!!$    end select
!!$
!!$    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
!!$    if (myid == 0) then
!!$       write(6,*) 'Mij phi prefactor proc0:', maxval(mij_prefact(:,isim,1)), &
!!$            maxval(mij_prefact(:,isim,2)), maxval(mij_prefact(:,isim,3))
!!$    end if
!!$    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)
!!$    if (myid == 1) then
!!$       write(6,*) 'Mij phi prefactor proc1:', maxval(mij_prefact(:,isim,1)), &
!!$            maxval(mij_prefact(:,isim,2)), maxval(mij_prefact(:,isim,3))
!!$    end if
!!$    call MPI_barrier(MPI_COMM_WORLD, ierr_mpi)

  end subroutine compute_prefactor
!--------------------------------------------------------------------------------


!================================================================================
! Compute cylinfrical coordinates 3D (vel)
  subroutine compute_3D_cyl(isim)

    use global_parameters_mod

    integer(kind=si) :: irec
    integer(kind=si), intent(in) :: isim

    do irec=irecmin,irecmax
       data_rec_all(irec,1) = data_rec_all(irec,1) + f1(isim,irec) * data_rec(irec,1) !f1
       data_rec_all(irec,2) = data_rec_all(irec,2) + f2(isim,irec) * data_rec(irec,2)          !f2
       data_rec_all(irec,3) = data_rec_all(irec,3) + f1(isim,irec) * data_rec(irec,3)          !f1
    end do

  end subroutine compute_3D_cyl
!--------------------------------------------------------------------------------


!================================================================================
! Compute cylindrical coordinates 3D (stress)
  subroutine compute_stress_3D_cyl(isim)
    
    use global_parameters_mod

    integer(kind=si) :: irec
    integer(kind=si), intent(in) :: isim

    do irec=irecmin,irecmax
       stress_rec_all(irec,1) = stress_rec_all(irec,1) + f1(isim,irec) * stress_rec(irec,1) !f1
       stress_rec_all(irec,2) = stress_rec_all(irec,2) + f1(isim,irec) * stress_rec(irec,2) !f1
       stress_rec_all(irec,3) = stress_rec_all(irec,3) + f1(isim,irec) * stress_rec(irec,3) !f1
       stress_rec_all(irec,4) = stress_rec_all(irec,4) + f2(isim,irec) * stress_rec(irec,4) !f2
       stress_rec_all(irec,5) = stress_rec_all(irec,5) + f1(isim,irec) * stress_rec(irec,5) !f1
       stress_rec_all(irec,6) = stress_rec_all(irec,6) + f2(isim,irec) * stress_rec(irec,6) !f2
    end do
    
  end subroutine compute_stress_3D_cyl
!--------------------------------------------------------------------------------

end module post_process_wavefields_mod
