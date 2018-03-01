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

    integer(kind=si)   :: iin_file=12, ioerr
    character(len=256) :: line
    character(len=256) :: keyword, keyvalue

    
print *,'lalalalalalalalalala'
        write(6, '(A)', advance='no')' Reading lithos.par for recon infos...'
    open(unit=iin_file,file='lithos.par', status='old', action='read', iostat=ioerr)
    if (ioerr /= 0) stop 'Check input file ''lithos.par''! Is it still there?' 
    
    do 
       
       read(iin_file,fmt='(a256)',iostat=ioerr) line
print *,line
       if (ioerr < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
       read(line,*) keyword, keyvalue
       
       parameters_recon : select case(trim(keyword))
       case('modeling_tool')
          read(keyvalue, *) fwdtool
          
       case('coupling_box_file')
          read(keyvalue, *) input_point_file
          
       case('coupling_src_latlon')
          read(line, *) keyword, lat_src, lon_src
          
       case('coupling_axisem_nsim')
          read(keyvalue, *) nsim
          
       case('nb_partitions')
          read(keyvalue, *) npart
          
       case('mesh_coordinates')
          read(line, *) keyword, lat_mesh, lon_mesh, alpha_mesh
          
       end select parameters_recon
    end do
    close(iin_file)
    write(6,*)'Done!'
    
!    open(10,file='reconstruction.par',status='old') 
!    read(10,'(a)')fwdtool
!    read(10,'(a)')input_point_file
!    read(10,*)dummy
!    read(10,*)lat_src,lon_src
!    read(10,*) lat_mesh,lon_mesh,alpha_mesh
!    read(10,*) nsim
!    close(10)

    if (allocated(working_axisem_dir)) deallocate(working_axisem_dir)
    if(.not.allocated(working_axisem_dir)) allocate(working_axisem_dir(nsim))

    nsim = 1
    
    if (nsim == 1) then
       working_axisem_dir(1)="./"
    else
       working_axisem_dir(1) = "MZZ/"
       working_axisem_dir(2) = "MXX_P_MYY/"
       working_axisem_dir(3) = "MXZ_MYZ/"
       working_axisem_dir(4) = "MXY_MXX_M_MYY/"
    end if
    
    !open(10,file=trim( working_axisem_dir(1))//'Data/strain_info.dat0000')
    open(10,file=trim( working_axisem_dir(1))//'simulation.info')
    do i=1,15
	read(10,*)dummy
    end do
    read(10,*)ntold,dummy
    read(10,*)dtold,dummy
    close(10)
    !read(10,*) ntold
    !read(10,*) dtold,i
    !close(10)
    print *,ntold, dtold,i

  end subroutine read_reconstruction_param


  subroutine define_axisem_dir !(ipart)

!    integer(kind=si), intent(in) :: ipart
 !   character(len=5)   :: myfileend

 !   write(myfileend,'(a,i4.4)')'_',ipart

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
    if (allocated(ivx)) deallocate(ivx)
    if (allocated(ivy)) deallocate(ivy)
    if (allocated(ivz)) deallocate(ivz)
    if (allocated(isxx)) deallocate(isxx)
    if (allocated(isyy)) deallocate(isyy)
    if (allocated(iszz)) deallocate(iszz)
    if (allocated(isxy)) deallocate(isxy)
    if (allocated(isyz)) deallocate(isyz)
    if (allocated(isxz)) deallocate(isxz)
    allocate(ivx(nsim),ivy(nsim),ivz(nsim))
    allocate(isxx(nsim),isyy(nsim),iszz(nsim))
    allocate(isxy(nsim),isxz(nsim),isyz(nsim))
    
    !*** new
    nsim = 1
    working_axisem_dir(1)="./"
       
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
       
       write(fichier,'(a15)') output_veloc_name(1)
       open(ivx(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_veloc_name(2)
       open(ivy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_veloc_name(3)
       open(ivz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(1)
       open(isxx(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(2)
       open(isyy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(3)
       open(iszz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(4)
       open(isxy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(5)
       open(isxz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
       write(fichier,'(a15)') output_stress_name(6)
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
    
  end subroutine define_axisem_dir
!--------------------------------------------------------------------------------

!================================================================================
! Interpolation    parameter file
  subroutine read_interpolation_param

    integer(kind=si)   :: iin_file=12, ioerr
    character(len=256) :: line
    character(len=256) :: keyword, keyvalue

    write(6, '(A)', advance='no')' Reading lithos.par for interp infos...'
    open(unit=iin_file,file='lithos.par', status='old', action='read', iostat=ioerr)
    if (ioerr /= 0) stop 'Check input file ''lithos.par''! Is it still there?' 
    
    do 
       
       read(iin_file,fmt='(a256)',iostat=ioerr) line
       if (ioerr < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
       read(line,*) keyword, keyvalue
       
       parameters_to_read : select case(trim(keyword))
       case('coupling_interp_rep')
          read(keyvalue, '(a)') rep

       case('coupling_starttime_index')
          read(keyvalue, *) ind

       case('coupling_interp_stf')
          read(keyvalue, *) isconv
        if (isconv == 1) then
             read(line, *) keyword, isconv, ntstf, dtstf, stf_file
          end if
          
       case('coupling_interp_filter')
          read(keyvalue, *) istap
          if (istap == 1) then
             read(line, *) keyword, istap, alpha, fmax
          end if
          
       end select parameters_to_read
       
    end do
    close(iin_file)
    write(6,*)'Done!'
    
 
!    open(10,file='interpolation.par')
!    read(10,*)fwdtool
!    read(10,'(a)')rep
!    read(10,*)nsta,nlta,thres
!    read(10,*)istap,alpha
!    read(10,*)isconv,ntstf,dtstf
!    read(10,*)stf_file
!    read(10,*)dummy
!    read(10,*)fmax
!    close(10)

  end subroutine read_interpolation_param
!--------------------------------------------------------------------------------
 
!================================================================================
! Read source time function to convolve
  subroutine read_stf

    real(kind=cp), dimension(:), allocatable :: tmpstf
    real(kind=cp) :: festf, fenew, told

    integer(kind=si) :: itnew, itold, ntwe

    !*** Read STF
    if(.not.allocated(tmpstf)) allocate(tmpstf(ntstf))     
    
    open(10,file=trim(stf_file),access='direct',recl=cp*ntstf)
    read(10,rec=1)tmpstf
    close(10)
    
    !*** Interpolate
    write(6,*)'Interpolate source time function...'

    festf = 1. / dtstf
    fenew = 1. / dtold
    told  = ntstf * dtstf

    ntwe = ceiling(told * fenew)

    if(.not.allocated(stf)) allocate(stf(ntwe))
    stf(:) = 0._cp
   
    do itnew = 1, ntwe
       do itold = 1, ntstf
          stf(itnew) = stf(itnew) + tmpstf(itold) * mysinc(festf * (real(itnew-1) * dtold - real(itold -1)*dtstf))
       end do
    end do
    if(allocated(tmpstf)) deallocate(tmpstf)

    ntstf = ntwe

    write(6,*)'Done !'
    
  contains

    real(kind=cp) function mysinc(x)

      real(kind=cp) :: x
      real(kind=cp), parameter :: pipi=3.141592653589793
    
      if (abs(x) >= 1e-13) then
         mysinc = sin(pipi*x)/(pipi*x)
      else
         mysinc = 1
      end if
      
    end function mysinc

  end subroutine read_stf
!--------------------------------------------------------------------------------

!================================================================================
! Read modeling tools parameter
  subroutine read_modeling_params
    
    integer(kind=si)   :: iin_file=12, ioerr
    character(len=256) :: line
    character(len=256) :: keyword, keyvalue


    select case (fwdtool)
    case('SEM')

       write(6, '(A)', advance='no')' Reading lithos.par for general infos...'
       open(unit=iin_file,file='lithos.par', status='old', action='read', iostat=ioerr)
       if (ioerr /= 0) stop 'Check input file ''lithos.par''! Is it still there?' 
       
       do 
          
          read(iin_file,fmt='(a256)',iostat=ioerr) line
          if (ioerr < 0) exit
          if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
          read(line,*) keyword, keyvalue
          
          parameters_to_read : select case(trim(keyword))
          case('modeling_order')
             read(keyvalue, *) ord
             ngll = ord+1
             ngllsquare = ngll * ngll

          case('simulation_duration')
             read(keyvalue, *) tsimu

          case('simulation_time_step')
             read(keyvalue, *) dtnew

          case('simulation_nb_time_step')
             read(keyvalue, *) ntnew

          case('simulation_nb_buffer')
             read(keyvalue, *) tbuff

          end select parameters_to_read

       end do
       close(iin_file)

       print *,dtnew, ntnew

    case('FD')
       stop 'Not implemented yet....'
    end select
    

  end subroutine read_modeling_params
!--------------------------------------------------------------------------------

end module inputs_outputs_interp_mod






