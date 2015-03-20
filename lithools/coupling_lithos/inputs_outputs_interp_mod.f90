module inputs_outputs_interp_mod

  use precision_mod
  use mpi_int_mod
  use interpolation_parameters_mod
  use interp_process_mod

  implicit none

contains

!================================================================================
! Read all inputs for interpolation
  subroutine read_all_inputs

    if (myid ==0) then

       !*** Read parameters files
       call read_interpolation_param    !** Add here something with respect to convolution
       call read_reconstruction_param
       
       !*** Read modeling tools param
       call read_modeling_params
       if (isconv == 1) then
           call read_stf
       end if

    end if
       
  end subroutine read_all_inputs
!--------------------------------------------------------------------------------

!================================================================================
! Read reconstruction parameters file and prepare the reading of data
  subroutine read_reconstruction_param

    open(10,file='../reconstruction.par',status='old') 
    read(10,'(a)')fwdtool
    read(10,'(a)')input_point_file
    read(10,*)dummy
    read(10,*)lat_src,lon_src
    read(10,*) lat_mesh,lon_mesh,alpha_mesh
    read(10,*) nsim
    close(10)

    if(.not.allocated(working_axisem_dir)) allocate(working_axisem_dir(nsim))

    if (nsim == 1) then
       working_axisem_dir(1)="./"
    else
       working_axisem_dir(1) = "MZZ/"
       working_axisem_dir(2) = "MXX_P_MYY/"
       working_axisem_dir(3) = "MXZ_MYZ/"
       working_axisem_dir(4) = "MXY_MXX_M_MYY/"
    end if
    
    open(10,file=trim( working_axisem_dir(1))//'Data/strain_info.dat0000')
    read(10,*) ntold
    read(10,*) dtold,i
    close(10)
    print *,ntold, dtold,i
  
    output_veloc_name(1)  = 'velocityoutp_u1'
    output_veloc_name(2)  = 'velocityoutp_u2'
    output_veloc_name(3)  = 'velocityoutp_u3'
    
    output_stress_name(1) = 'stress_Sg11_out'
    output_stress_name(2) = 'stress_Sg22_out'
    output_stress_name(3) = 'stress_Sg33_out'
    output_stress_name(4) = 'stress_Sg12_out'
    output_stress_name(5) = 'stress_Sg13_out'
    output_stress_name(6) = 'stress_Sg23_out'
    
    iunit=1000
    allocate(ivx(nsim),ivy(nsim),ivz(nsim))
    allocate(isxx(nsim),isyy(nsim),iszz(nsim))
    allocate(isxy(nsim),isxz(nsim),isyz(nsim))
    
    do isim=1,nsim
       ivx(isim)=next_iunit(iunit)
       ivy(isim)=next_iunit(iunit)
       ivz(isim)=next_iunit(iunit)
       isxx(isim)=next_iunit(iunit)
       isyy(isim)=next_iunit(iunit)
       iszz(isim)=next_iunit(iunit)
       isxy(isim)=next_iunit(iunit)
       isxz(isim)=next_iunit(iunit)
       isyz(isim)=next_iunit(iunit)
       
       write(fichier,'(a6,a15)') '/Data/',output_veloc_name(1)
       open(ivx(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_veloc_name(2)
       open(ivy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_veloc_name(3)
       open(ivz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(1)
       open(isxx(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(2)
       open(isyy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(3)
       open(iszz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(4)
       open(isxy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(5)
       open(isxz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a6,a15)') '/Data/',output_stress_name(6)
       open(isyz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
    end do
    
    do isim=1,nsim
       read(ivx(isim))  nbrec,ntime
       read(ivy(isim))  nbrec,ntime
       read(ivz(isim))  nbrec,ntime
       read(isxx(isim)) nbrec,ntime
       read(isyy(isim)) nbrec,ntime
       read(iszz(isim)) nbrec,ntime
       read(isxy(isim)) nbrec,ntime
       read(isxz(isim)) nbrec,ntime
       read(isyz(isim)) nbrec,ntime
    end do
    
    print *,'nbrec ntime',nbrec,ntime
    print *,'if problem read one integer...'
    
    npts = nbrec
    
  end subroutine read_reconstruction_param
!--------------------------------------------------------------------------------

!================================================================================
! Interpolation    parameter file
  subroutine read_interpolation_param

    open(10,file='../interpolation.par')
    read(10,*)fwdtool
    read(10,'(a)')rep
    read(10,*)nsta,nlta,thres
    read(10,*)istap,alpha
    read(10,*)isconv,ntstf
    read(10,*)stf_file
    read(10,*)dummy
    read(10,*)fmax
    close(10)

  end subroutine read_interpolation_param
!--------------------------------------------------------------------------------
 
!================================================================================
! Read source time function to convolve
  subroutine read_stf

     !*** Read
     if(.not.allocated(stf)) allocate(stf(ntstf))     

     open(10,file=trim(stf_file),access='direct',recl=cp*ntstf)
     read(10,rec=1)stf
     close(10)

     !*** Interpolate to axisem 


     


  end subroutine read_stf
!--------------------------------------------------------------------------------

!================================================================================
! Read modeling tools parameter
  subroutine read_modeling_params
    
    select case (fwdtool)
    case('SEM')
       print *,trim(rep)//'/sem.par'
       open(10,file=trim(rep)//'/sem.par',action='read')
       read(10,*)dummy
       read(10,*)dummy
       read(10,*)dummy
       read(10,*)nelx,nely,nelz,dummy
       read(10,*)ngllx,nglly,ngllz,dummy
       read(10,*)dummy
       read(10,*)dtnew,ntnew,dummy
       read(10,*)dummy
       read(10,*)num_bnd_faces,dummy
       close(10)
       
       nelem = nelx * nely * nelz
       ngllsquare = ngllx * nglly
       nptsa = (nelx * (ngllx -1) +1)*(nely * (nglly -1) +1)*(nelz * (ngllz -1) +1)
       
       print *,dtnew, ntnew
       
       if(.not.allocated(loc2glob)) allocate(loc2glob(ngllx,nglly,ngllz,nelem))
       if(.not.allocated(xcoord)) allocate(xcoord(nptsa))
       if(.not.allocated(ycoord)) allocate(ycoord(nptsa))
       if(.not.allocated(zcoord)) allocate(zcoord(nptsa))
       if(.not.allocated(abs_bnd_ielem))       allocate(abs_bnd_ielem(num_bnd_faces))
       if(.not.allocated(abs_bnd_ijk))         allocate(abs_bnd_ijk(3,ngllsquare,num_bnd_faces))
       if(.not.allocated(abs_bnd_jacobian2Dw)) allocate(abs_bnd_jacobian2Dw(ngllsquare,num_bnd_faces))
       if(.not.allocated(abs_bnd_normal))      allocate(abs_bnd_normal(3,ngllsquare,num_bnd_faces))
       
       open(unit=27,file=rep(1:len_trim(rep))//'/MESH/x.bin',status='old',form='unformatted')
       read(27) xcoord
       close(27)
       open(unit=27,file=rep(1:len_trim(rep))//'/MESH/y.bin',status='old',form='unformatted')
       read(27) ycoord
       close(27)
       open(unit=27,file=rep(1:len_trim(rep))//'/MESH/z.bin',status='old',form='unformatted')
       read(27) zcoord
       close(27)
       open(unit=27,file=rep(1:len_trim(rep))//'/MESH/ibool.bin',status='old',form='unformatted')
       read(27) loc2glob
       close(27)
       open(unit=27,file=rep(1:len_trim(rep))//'/MESH/boundary.bin',status='old',form='unformatted')
       read(27) abs_bnd_ielem
       read(27) abs_bnd_ijk
       read(27) abs_bnd_normal
       read(27) abs_bnd_jacobian2Dw
       close(27)
       
    case('FD')
       stop 'Not implemented yet....'
    end select
    

  end subroutine read_modeling_params
!--------------------------------------------------------------------------------

end module inputs_outputs_interp_mod






