module rotation_matrix_mod

  use precision_mod
  use constants_mod

  implicit none

contains

!================================================================================
! Rotate in global cartesian with source still in pole
  subroutine rotate2cartesian_with_source_in_pole_stress

    use global_parameters_mod
    
    integer(kind=si) :: irec, i, j, k
    real(kind=cp), dimension(3,3) :: tmp1, B, st, Bt
    real(kind=cp), dimension(6,6) :: tmp

    ! compute B*st*Bt
    do irec=irecmin,irecmax
          
       ! rotation matrix 
       Bt(1,1)=  cos(phi(irec))
       Bt(1,2)= - sin(phi(irec))
       Bt(1,3)= 0. 
       
       Bt(2,1)=  sin(phi(irec))
       Bt(2,2)=  cos(phi(irec))
       Bt(2,3)=  0.

       Bt(3,1)= 0. 
       Bt(3,2)= 0. 
       Bt(3,3)= 1. 
       
       B=transpose(Bt)

       st(1,1)=stress_rec_all(irec,1)
       st(1,2)=stress_rec_all(irec,4)
       st(1,3)=stress_rec_all(irec,5)
       
       st(2,1)=stress_rec_all(irec,4)
       st(2,2)=stress_rec_all(irec,2)
       st(2,3)=stress_rec_all(irec,6)

       st(3,1)=stress_rec_all(irec,5)
       st(3,2)=stress_rec_all(irec,6)
       st(3,3)=stress_rec_all(irec,3)

       ! st*Bt
       tmp=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp(i,j)=tmp(i,j)+st(i,k)*B(k,j)   !*B(j,k)
             end do
          end do
       end do

       ! B*st*Bt
       tmp1=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp1(i,j)=tmp1(i,j) + Bt(i,k)*tmp(k,j) !B(i,k)*tmp(k,j)
             end do
          end do
       end do
       
       ! stress in cartesian coordinates 
       stress_rec_all(irec,1)=tmp1(1,1)
       stress_rec_all(irec,2)=tmp1(2,2)
       stress_rec_all(irec,3)=tmp1(3,3)
       stress_rec_all(irec,4)=tmp1(1,2) 
       stress_rec_all(irec,5)=tmp1(1,3)
       stress_rec_all(irec,6)=tmp1(2,3)
       
    end do
    
  end subroutine rotate2cartesian_with_source_in_pole_stress
!--------------------------------------------------------------------------------

!================================================================================
! Put the source from pole to right location
  subroutine rotate_back_source

    use global_parameters_mod
    
    integer(kind=si) :: irec, i, j, k
    real(kind=sp), dimension(3) :: tmp, veloc
    
    do irec=irecmin,irecmax
       
       ! veloc in cylindical coordinates
       veloc(1)=data_rec_all(irec,1)
       veloc(2)=data_rec_all(irec,2)
       veloc(3)=data_rec_all(irec,3)
       
       ! R*veloc
       tmp=0.  
       do i=1,3
          do k=1,3
             tmp(i)=tmp(i)+veloc(k)*rot_mat(i,k) !trans_rot_mat(i,k)
          end do
       end do
       
       ! valocity in cartesian
       data_rec_all(irec,1)=tmp(1)
       data_rec_all(irec,2)=tmp(2)
       data_rec_all(irec,3)=tmp(3)
       
    end do
      
  end subroutine rotate_back_source
!--------------------------------------------------------------------------------

!================================================================================
! Rotate to local cartesian 
  subroutine rotate_back_to_local_cart

    use global_parameters_mod

    integer(kind=si) :: irec, i, j, k
    real(kind=sp), dimension(3) :: tmp, veloc

    do irec=irecmin,irecmax
           
       ! veloc in global coordinates
       veloc(1)=data_rec_all(irec,1)
       veloc(2)=data_rec_all(irec,2)
       veloc(3)=data_rec_all(irec,3)
          
       ! Rt*veloc
       tmp=0.  
       do i=1,3
          do k=1,3
             tmp(i)=tmp(i)+veloc(k)*rot_mat_mesh(i,k)  
          end do
       end do
       

       select case (coup_tool)
       case('DG')
          ! valocity in cartesian
          data_rec_all(irec,1)=tmp(1)
          data_rec_all(irec,2)=tmp(2)
          data_rec_all(irec,3)=tmp(3)
       case('SEM')
          data_rec_all(irec,1)=tmp(2)
          data_rec_all(irec,2)=tmp(3)
          data_rec_all(irec,3)=tmp(1)
       case('FD')
          data_rec_all(irec,1)=tmp(2)
          data_rec_all(irec,2)=tmp(3)
          data_rec_all(irec,3)=tmp(1)
       end select
    end do

  end subroutine rotate_back_to_local_cart
!--------------------------------------------------------------------------------

!================================================================================
! Rotate source back for stress
  subroutine rotate_back_source_stress

    use global_parameters_mod

    integer(kind=si) :: irec, i, j, k
    real(kind=sp), dimension(3,3) :: tmp, tmp1, st

    do irec=irecmin,irecmax
       
       ! stress in cylindical coordinates
       st(1,1)=stress_rec_all(irec,1)
       st(1,2)=stress_rec_all(irec,4)
       st(1,3)=stress_rec_all(irec,5)
       
       st(2,1)=stress_rec_all(irec,4)
       st(2,2)=stress_rec_all(irec,2)
       st(2,3)=stress_rec_all(irec,6)
       
       st(3,1)=stress_rec_all(irec,5)
       st(3,2)=stress_rec_all(irec,6)
       st(3,3)=stress_rec_all(irec,3)
       
       ! st*Rt
       tmp=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp(i,j)=tmp(i,j)+st(i,k) * trans_rot_mat(k,j)   !*trans_rot_mat(k,j)
             end do
          end do
       end do
       
       ! R*st*Rt
       tmp1=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp1(i,j)=tmp1(i,j)+tmp(k,j) * rot_mat(i,k) !*rot_mat(i,k)
             end do
          end do
       end do
       
       ! stress in cartesian
       stress_rec_all(irec,1)=tmp1(1,1)
       stress_rec_all(irec,2)=tmp1(2,2)
       stress_rec_all(irec,3)=tmp1(3,3)
       stress_rec_all(irec,4)=tmp1(1,2) 
       stress_rec_all(irec,5)=tmp1(1,3)
       stress_rec_all(irec,6)=tmp1(2,3)
       
    end do
    
  end subroutine rotate_back_source_stress
!--------------------------------------------------------------------------------

!================================================================================
! Rotate local cart (stress)
  subroutine rotate_back_to_local_cart_stress

    use global_parameters_mod

    integer(kind=si) :: irec, i, j, k
    real(kind=sp), dimension(3,3) :: tmp, tmp1, st

    do irec=irecmin,irecmax
           
       ! stress in cylindical coordinates
       st(1,1)=stress_rec_all(irec,1)
       st(1,2)=stress_rec_all(irec,4)
       st(1,3)=stress_rec_all(irec,5)
       
       st(2,1)=stress_rec_all(irec,4)
       st(2,2)=stress_rec_all(irec,2)
       st(2,3)=stress_rec_all(irec,6)
       
       st(3,1)=stress_rec_all(irec,5)
       st(3,2)=stress_rec_all(irec,6)
       st(3,3)=stress_rec_all(irec,3)
       
       ! st*R
       tmp=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp(i,j)=tmp(i,j)+st(i,k)*trans_rot_mat_mesh(k,j)
             end do
          end do
       end do
       
       ! Rt*st*R
       tmp1=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp1(i,j)=tmp1(i,j)+rot_mat_mesh(i,k)*tmp(k,j) 
             end do
          end do
       end do


       select case (coup_tool)
       case('DG')
          ! stress in cartesian
          stress_rec_all(irec,1)=tmp1(1,1)
          stress_rec_all(irec,2)=tmp1(2,2)
          stress_rec_all(irec,3)=tmp1(3,3)
          stress_rec_all(irec,4)=tmp1(1,2) 
          stress_rec_all(irec,5)=tmp1(1,3)
          stress_rec_all(irec,6)=tmp1(2,3)
       case('SEM')
          stress_rec_all(irec,1)=tmp1(2,2)
          stress_rec_all(irec,2)=tmp1(3,3)
          stress_rec_all(irec,3)=tmp1(1,1)
          stress_rec_all(irec,4)=tmp1(2,3)
          stress_rec_all(irec,5)=tmp1(2,1)
          stress_rec_all(irec,6)=tmp1(3,1)
       case('FD')
          stress_rec_all(irec,1)=tmp1(2,2)
          stress_rec_all(irec,2)=tmp1(3,3)
          stress_rec_all(irec,3)=tmp1(1,1)
          stress_rec_all(irec,4)=tmp1(2,3)
          stress_rec_all(irec,5)=tmp1(2,1)
          stress_rec_all(irec,6)=tmp1(3,1)
       end select
    end do
    
  end subroutine rotate_back_to_local_cart_stress
!--------------------------------------------------------------------------------

!================================================================================
! Cartesian source pole rotation
  subroutine rotate2cartesian_with_source_in_pole
    
    use global_parameters_mod
    
    integer(kind=si) :: irec, i, k
    real(kind=sp), dimension(3)   :: tmp, veloc
    real(kind=sp), dimension(3,3) :: B, Bt

    do irec=irecmin,irecmax
           
       ! rotation matrix
       Bt(1,1)=  cos(phi(irec))
       Bt(1,2)= - sin(phi(irec))
       Bt(1,3)= 0. 
       
       Bt(2,1)=  sin(phi(irec))
       Bt(2,2)=  cos(phi(irec))
       Bt(2,3)=  0.
       
       Bt(3,1)= 0. 
       Bt(3,2)= 0. 
       Bt(3,3)= 1. 
       
       B = transpose(Bt)

       ! veloc in cylindical coordinates
       veloc(1)=data_rec_all(irec,1)
       veloc(2)=data_rec_all(irec,2)
       veloc(3)=data_rec_all(irec,3)

       ! B*veloc
       tmp=0.  
       do i=1,3
          do k=1,3
             tmp(i)=tmp(i)+veloc(k)*Bt(i,k)
          end do
       end do
       
       ! valocity in cartesian
       data_rec_all(irec,1)=tmp(1)
       data_rec_all(irec,2)=tmp(2)
       data_rec_all(irec,3)=tmp(3)
       
    end do
    
  end subroutine rotate2cartesian_with_source_in_pole
!--------------------------------------------------------------------------------


!================================================================================
! Rotation matrix for source
  subroutine def_rot_matrix(srccolat,srclon,rot_mat,trans_rot_mat)
      
    real(kind=dp), dimension(3,3), intent(out)  :: rot_mat,trans_rot_mat
    integer(kind=si)               :: i, j
    real(kind=dp)                  :: smallval
    real(kind=dp), intent(in)      :: srccolat,srclon
      
    smallval=1e-11_dp
      
    ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007.
    rot_mat(1,1) = dcos(srccolat) * dcos(srclon)
    rot_mat(2,2) = dcos(srclon)
    rot_mat(3,3) = dcos(srccolat)
    rot_mat(2,1) = dcos(srccolat) * dsin(srclon)
    rot_mat(3,1) = -dsin(srccolat)
    rot_mat(3,2) = 0._dp
    rot_mat(1,2) = -dsin(srclon)
    rot_mat(1,3) = dsin(srccolat) * dcos(srclon)
    rot_mat(2,3) = dsin(srccolat) * dsin(srclon)
    
    where (abs(rot_mat)<smallval) rot_mat = 0._dp
    
    trans_rot_mat = transpose(rot_mat)
    
    write(*,*) ' ROTATION MATRIX '
    write(*,*) rot_mat(1,1),rot_mat(1,2),rot_mat(1,3)
    write(*,*) rot_mat(2,1),rot_mat(2,2),rot_mat(2,3)
    write(*,*) rot_mat(3,1),rot_mat(3,2),rot_mat(3,3)
    write(*,*) 
  
  end subroutine def_rot_matrix
!--------------------------------------------------------------------------------


!================================================================================
! Rotation matrix for DG mesh
  subroutine def_rot_matrix_DG(srccolat,srclon,rot_mat,trans_rot_mat)
      
    real(kind=dp), dimension(3,3), intent(out)  :: rot_mat,trans_rot_mat
    integer(kind=si)               :: i, j
    real(kind=dp)                  :: smallval
    real(kind=dp), intent(in)      :: srccolat,srclon
      
    smallval=1e-11_dp
    
    rot_mat(1,2) = -dcos(srccolat)*dcos(srclon)
    rot_mat(2,1) =  dcos(srclon)
    rot_mat(3,3) =  dcos(srccolat)

    rot_mat(2,3) = dsin(srccolat) * dsin(srclon)
    rot_mat(1,3) = dsin(srccolat) * dcos(srclon)
    rot_mat(1,1) = -dsin(srclon)

    rot_mat(2,2) = -dcos(srccolat)*dsin(srclon)
    rot_mat(3,1) = 0._dp
    rot_mat(3,2) = dsin(srccolat) !sin(srccolat) * sin(srclon)

 
    where (abs(rot_mat)<smallval) rot_mat = 0.0_dp
    
    write(*,*) rot_mat(1,1),rot_mat(1,2),rot_mat(1,3)
    write(*,*) rot_mat(2,1),rot_mat(2,2),rot_mat(2,3)
    write(*,*) rot_mat(3,1),rot_mat(3,2),rot_mat(3,3)
    
    trans_rot_mat = transpose(rot_mat)
    
  end subroutine def_rot_matrix_DG
!--------------------------------------------------------------------------------


!================================================================================
! Roation matrix for SEM mesh
  subroutine def_rot_matrix_SEM(lat,lon,alpha,rotmat,transrotmat)

    real(kind=dp), intent(in)                  :: lon, lat, alpha
    real(kind=dp), dimension(3,3), intent(out) :: rotmat, transrotmat
    real(kind=dp)                              :: smallval

    smallval=1e-11_dp

    rotmat(1,1) =  dcos(lat) * dcos(lon)
    rotmat(1,2) =  dcos(lat) * dsin(lon)
    rotmat(1,3) =  dsin(lat)
    rotmat(2,1) = -dsin(lon) * dcos(alpha) - dsin(alpha) * dsin(lat) * dcos(lon)
    rotmat(2,2) =  dcos(lon) * dcos(alpha) - dsin(alpha) * dsin(lat) * dsin(lon)
    rotmat(2,3) =  dsin(alpha) * dcos(lat)
    rotmat(3,1) =  dsin(lon) * dsin(alpha) - dcos(alpha) * dsin(lat) * dcos(lon)
    rotmat(3,2) = -dcos(lon) * dsin(alpha) - dcos(alpha) * dsin(lat) * dsin(lon)
    rotmat(3,3) =  dcos(alpha) * dcos(lat)

    where (abs(rotmat)<smallval) rotmat = 0.0_dp

    transrotmat = transpose(rotmat)

    write(*,*) 'MESH ROT SEM'
    write(*,*) rotmat(1,1),rotmat(1,2),rotmat(1,3)
    write(*,*) rotmat(2,1),rotmat(2,2),rotmat(2,3)
    write(*,*) rotmat(3,1),rotmat(3,2),rotmat(3,3)


  end subroutine def_rot_matrix_SEM
!--------------------------------------------------------------------------------

!================================================================================
! Rotate the box points
  subroutine rotate_box(r,th,ph,trans_rot_mat)

    real(kind=dp), dimension(3,3), intent(in) :: trans_rot_mat
    real(kind=dp), intent(inout) :: r,th,ph
    real(kind=dp), dimension(3)  :: x_vec, x_vec_rot
    real(kind=dp)                :: r_r, smallval_dble 
    
    smallval_dble=0.0_dp   !! 1e-11  ! a quoi sert ce truc? ! VM VM 
    !! verifier l'effet que çà peut avoir de le mettre a zero
    
    x_vec(1) = r * dsin(th) * dcos(ph)
    x_vec(2) = r * dsin(th) * dsin(ph)
    x_vec(3) = r * dcos(th) 
    
    x_vec_rot = matmul(trans_rot_mat,x_vec)
    
    r_r = dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2 )
    th = dacos((x_vec_rot(3)  + smallval_dble )/ ( r_r + smallval_dble) )
    ph = datan2(x_vec_rot(2),x_vec_rot(1))
    
  end subroutine rotate_box
!--------------------------------------------------------------------------------

end module rotation_matrix_mod
