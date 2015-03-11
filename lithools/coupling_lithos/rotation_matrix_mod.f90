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
    real(kind=sp), dimension(3,3) :: tmp1, B, st
    real(kind=sp), dimension(6,6) :: tmp

    ! compute B*st*Bt
    do irec=irecmin,irecmax
          
       ! rotation matrix 
       B(1,1)=  cos(phi(irec))
       B(1,2)= - sin(phi(irec))
       B(1,3)= 0. 
       
       B(2,1)=  sin(phi(irec))
       B(2,2)=  cos(phi(irec))
       B(2,3)=  0.

       B(3,1)= 0. 
       B(3,2)= 0. 
       B(3,3)= 1. 
       
       st(1,1)=stress_rec(irec,1)
       st(1,2)=stress_rec(irec,4)
       st(1,3)=stress_rec(irec,5)
       
       st(2,1)=stress_rec(irec,4)
       st(2,2)=stress_rec(irec,2)
       st(2,3)=stress_rec(irec,6)

       st(3,1)=stress_rec(irec,5)
       st(3,2)=stress_rec(irec,6)
       st(3,3)=stress_rec(irec,3)

       ! st*Bt
       tmp=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp(i,j)=tmp(i,j)+st(i,k)*B(j,k)
             end do
          end do
       end do

       ! B*st*Bt
       tmp1=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp1(i,j)=tmp1(i,j)+B(i,k)*tmp(k,j)
             end do
          end do
       end do
       
       ! stress in cartesian coordinates 
       stress_rec(irec,1)=tmp1(1,1)
       stress_rec(irec,2)=tmp1(2,2)
       stress_rec(irec,3)=tmp1(3,3)
       stress_rec(irec,4)=tmp1(1,2) 
       stress_rec(irec,5)=tmp1(1,3)
       stress_rec(irec,6)=tmp1(2,3)
       
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
       veloc(1)=data_rec(irec,1)
       veloc(2)=data_rec(irec,2)
       veloc(3)=data_rec(irec,3)
       
       ! R*veloc
       tmp=0.  
       do i=1,3
          do k=1,3
             tmp(i)=tmp(i)+veloc(k)*rot_mat(i,k)
          end do
       end do
       
       ! valocity in cartesian
       data_rec(irec,1)=tmp(1)
       data_rec(irec,2)=tmp(2)
       data_rec(irec,3)=tmp(3)
       
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
       veloc(1)=data_rec(irec,1)
       veloc(2)=data_rec(irec,2)
       veloc(3)=data_rec(irec,3)
          
       ! Rt*veloc
       tmp=0.  
       do i=1,3
          do k=1,3
             tmp(i)=tmp(i)+veloc(k)*trans_rot_mat_mesh(i,k)  
          end do
       end do
       
       ! valocity in cartesian
       data_rec(irec,1)=tmp(1)
       data_rec(irec,2)=tmp(2)
       data_rec(irec,3)=tmp(3)
              
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
       st(1,1)=stress_rec(irec,1)
       st(1,2)=stress_rec(irec,4)
       st(1,3)=stress_rec(irec,5)
       
       st(2,1)=stress_rec(irec,4)
       st(2,2)=stress_rec(irec,2)
       st(2,3)=stress_rec(irec,6)
       
       st(3,1)=stress_rec(irec,5)
       st(3,2)=stress_rec(irec,6)
       st(3,3)=stress_rec(irec,3)
       
       ! st*Rt
       tmp=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp(i,j)=tmp(i,j)+st(i,k)*trans_rot_mat(k,j)
             end do
          end do
       end do
       
       ! R*st*Rt
       tmp1=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp1(i,j)=tmp1(i,j)+tmp(k,j)*rot_mat(i,k)
             end do
          end do
       end do
       
       ! stress in cartesian
       stress_rec(irec,1)=tmp1(1,1)
       stress_rec(irec,2)=tmp1(2,2)
       stress_rec(irec,3)=tmp1(3,3)
       stress_rec(irec,4)=tmp1(1,2) 
       stress_rec(irec,5)=tmp1(1,3)
       stress_rec(irec,6)=tmp1(2,3)
       
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
       st(1,1)=stress_rec(irec,1)
       st(1,2)=stress_rec(irec,4)
       st(1,3)=stress_rec(irec,5)
       
       st(2,1)=stress_rec(irec,4)
       st(2,2)=stress_rec(irec,2)
       st(2,3)=stress_rec(irec,6)
       
       st(3,1)=stress_rec(irec,5)
       st(3,2)=stress_rec(irec,6)
       st(3,3)=stress_rec(irec,3)
       
       ! st*R
       tmp=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp(i,j)=tmp(i,j)+st(i,k)*rot_mat_mesh(k,j)
             end do
          end do
       end do
       
       ! Rt*st*R
       tmp1=0.
       do j=1,3
          do i=1,3
             do k=1,3
                tmp1(i,j)=tmp1(i,j)+trans_rot_mat_mesh(i,k)*tmp(k,j) 
             end do
          end do
       end do
       
       ! stress in cartesian
       stress_rec(irec,1)=tmp1(1,1)
       stress_rec(irec,2)=tmp1(2,2)
       stress_rec(irec,3)=tmp1(3,3)
       stress_rec(irec,4)=tmp1(1,2) 
       stress_rec(irec,5)=tmp1(1,3)
       stress_rec(irec,6)=tmp1(2,3)
       
    end do
    
  end subroutine rotate_back_to_local_cart_stress
!--------------------------------------------------------------------------------

!================================================================================
! Cartesian source pole rotation
  subroutine rotate2cartesian_with_source_in_pole
    
    use global_parameters_mod
    
    integer(kind=si) :: irec, i, k
    real(kind=sp), dimension(3)   :: tmp, veloc
    real(kind=sp), dimension(3,3) :: B

    do irec=irecmin,irecmax
           
       ! rotation matrix
       B(1,1)=  cos(phi(irec))
       B(1,2)= - sin(phi(irec))
       B(1,3)= 0. 
       
       B(2,1)=  sin(phi(irec))
       B(2,2)=  cos(phi(irec))
       B(2,3)=  0.
       
       B(3,1)= 0. 
       B(3,2)= 0. 
       B(3,3)= 1. 
       
       ! veloc in cylindical coordinates
       veloc(1)=data_rec(irec,1)
       veloc(2)=data_rec(irec,2)
       veloc(3)=data_rec(irec,3)

       ! B*veloc
       tmp=0.  
       do i=1,3
          do k=1,3
             tmp(i)=tmp(i)+veloc(k)*B(i,k)
          end do
       end do
       
       ! valocity in cartesian
       data_rec(irec,1)=tmp(1)
       data_rec(irec,2)=tmp(2)
       data_rec(irec,3)=tmp(3)
       
    end do
    
  end subroutine rotate2cartesian_with_source_in_pole
!--------------------------------------------------------------------------------


!================================================================================
! Rotation matrix for source
  subroutine def_rot_matrix(srccolat,srclon,rot_mat,trans_rot_mat)
      
    real(kind=cp), dimension(3,3), intent(out)  :: rot_mat,trans_rot_mat
    integer(kind=si)               :: i, j
    real(kind=cp)                  :: smallval
    real(kind=cp), intent(in)      :: srccolat,srclon
      
    smallval=1e-11
      
    ! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007.
    rot_mat(1,1) = dcos(srccolat) * dcos(srclon)
    rot_mat(2,2) = dcos(srclon)
    rot_mat(3,3) = dcos(srccolat)
    rot_mat(2,1) = dcos(srccolat) * dsin(srclon)
    rot_mat(3,1) = -dsin(srccolat)
    rot_mat(3,2) = 0.d0
    rot_mat(1,2) = -dsin(srclon)
    rot_mat(1,3) = dsin(srccolat) * dcos(srclon)
    rot_mat(2,3) = dsin(srccolat) * dsin(srclon)
    
    where (dabs(rot_mat)<smallval) rot_mat = 0.0
    
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
      
    real(kind=cp), dimension(3,3), intent(out)  :: rot_mat,trans_rot_mat
    integer(kind=si)               :: i, j
    real(kind=cp)                  :: smallval
    real(kind=cp), intent(in)      :: srccolat,srclon
      
    smallval=1e-11
    
    rot_mat(1,2) = -cos(srccolat)*cos(srclon)
    rot_mat(2,1) =  cos(srclon)
    rot_mat(3,3) =  cos(srccolat)

    rot_mat(2,3) = sin(srccolat) * sin(srclon)
    rot_mat(1,3) = sin(srccolat) * cos(srclon)
    rot_mat(1,1) = -sin(srclon)

    rot_mat(2,2) = -cos(srccolat)*sin(srclon)
    rot_mat(3,1) = 0.
    rot_mat(3,2) = sin(srccolat) !sin(srccolat) * sin(srclon)

 
    where (abs(rot_mat)<smallval) rot_mat = 0.0
    
    write(*,*) rot_mat(1,1),rot_mat(1,2),rot_mat(1,3)
    write(*,*) rot_mat(2,1),rot_mat(2,2),rot_mat(2,3)
    write(*,*) rot_mat(3,1),rot_mat(3,2),rot_mat(3,3)
    
    trans_rot_mat = transpose(rot_mat)
    
  end subroutine def_rot_matrix_DG
!--------------------------------------------------------------------------------


!================================================================================
! Rotate the box points
  subroutine rotate_box(r,th,ph,trans_rot_mat)

    real(kind=cp), dimension(3,3), intent(in) :: trans_rot_mat
    real(kind=cp), intent(inout) :: r,th,ph
    real(kind=cp), dimension(3)  :: x_vec, x_vec_rot,
    real(kind=cp)                :: r_r, smallval_dble 
    
    smallval_dble=0._dp !! 1e-11  ! a quoi sert ce truc? ! VM VM 
    !! verifier l'effet que çà peut avoir de le mettre a zero
    
    x_vec(1) = r * dsin(th) * dcos(ph)
    x_vec(2) = r * dsin(th) * dsin(ph)
    x_vec(3) = r * dcos(th) 
    
    x_vec_rot = matmul(trans_rot_mat,x_vec)
    
    r_r = dsqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2 )
    th = dacos((x_vec_rot(3)  + smallval_dble )/ ( r_r + smallval_dble) )
    ph = atan2(x_vec_rot(2),x_vec_rot(1))
    
  end subroutine rotate_box
!--------------------------------------------------------------------------------

end module rotation_matrix_mod
