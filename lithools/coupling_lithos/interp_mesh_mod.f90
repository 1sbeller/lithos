module interp_mesh_mod

  use precision_mod
  use global_parameters_mod, only: NGNOD, NGLLX, NGLLY

  real(kind=dp) :: a, b, c, d, det_jacobian

  real(kind=dp), dimension(NGLLX) :: hxir, hpxir, xigll, wxgll
  real(kind=dp), dimension(NGLLY) :: hpetar, hetar, yigll, wygll

  real(kind=dp), parameter :: GAUSSALPHA=0._dp, GAUSSBETA=0._dp
  real(kind=dp), parameter :: zero=0._dp, one=1._dp
  real(kind=dp), parameter :: three=3._dp, quart=0.25_dp, half=0.5_dp

contains 

!================================================================================
! Interpole wavefield in current element
  subroutine interpole_field(xi,eta,field,interp_value)
  
    real(kind=dp), intent(in)  :: xi, eta
    real(kind=dp), intent(out) :: interp_value

    integer(kind=si) :: igll, jgll
    real(kind=dp)    :: hlagrange
    real(kind=dp), dimension(NGLLX,NGLLY), intent(in) :: field

    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)

    call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
    call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)

    interp_value=0._dp
 
    do jgll = 1,NGLLY
       do igll = 1,NGLLX
                     
          !*** Lagrange polynomial
          hlagrange = hxir(igll)*hetar(jgll) 

          !*** Inteprolate value 
          interp_value = interp_value + field(igll,jgll)*hlagrange
                    
       end do
    end do
   
  end subroutine interpole_field
!--------------------------------------------------------------------------------

!================================================================================
! Find xi eta coordinate of point  
  subroutine find_xix_eta(nodes_crd,xi,eta,s_target,z_target)

    real(kind=dp), dimension(ngnod,2), intent(inout) :: nodes_crd
    real(kind=dp), intent(in)                       :: s_target, z_target
    real(kind=dp), intent(out)                    :: xi, eta

    integer(kind=si) :: inode, iguess, iter_newton, niter_newton

    real(kind=dp), dimension(ngnod) :: sph
    real(kind=dp) :: distmin, dist
    real(kind=dp) :: dxi, deta, s, z 
    
    niter_newton=6

    !*** Find the closest node 
    distmin=1e30_dp
    iguess=1
    do inode = 1, ngnod
       dist=(nodes_crd(inode,1)-s_target)**2 + (nodes_crd(inode,2)-z_target)**2
       if (dist < distmin) then
          distmin=dist
          iguess=inode
       end if
    end do

    !*** Convert to xi,eta initial guess 
    if (iguess==1) then
       xi=-1.
       eta=-1.
    end if
    if (iguess==2) then
       xi=0.
       eta=-1.
    end if
    if (iguess==3) then
       xi=1.
       eta=-1.
    end if
    if (iguess==4) then
       xi=1.
       eta=0.
    end if
    if (iguess==5) then
       xi=1.
       eta=1.
    end if
    if (iguess==6) then
       xi=0.
       eta=1.
    end if
    if (iguess==7) then
       xi=-1.
       eta=1.
    end if
    if (iguess==8) then
       xi=-1.
       eta=0.
    end if

    do iter_newton=1,niter_newton
       det_jacobian =  det_jacobian_shape(xi, eta, nodes_crd)
       call shp8(xi,eta,sph)
       s=0.
       z=0.
       do inode=1,ngnod
          s=s+sph(inode)*nodes_crd(inode,1)
          z=z+sph(inode)*nodes_crd(inode,2)
       end do
       
       dxi=(d*(s-s_target)  -b  *(z-z_target))/det_jacobian
       deta=(-c*(s-s_target) +a  *(z-z_target))/det_jacobian
       
       xi  = xi  - dxi
       eta = eta - deta
       
    end do

    call shp8(xi,eta,sph)
    s=0.
    z=0.
    do inode=1,ngnod
       s=s+sph(inode)*nodes_crd(inode,1)
       z=z+sph(inode)*nodes_crd(inode,2)
    end do
   
  end subroutine find_xix_eta
!--------------------------------------------------------------------------------

!================================================================================
! Shape function for quad8
  subroutine shp8(xil,etal,shp)
    !
    ! This routine computes and returns the quadratic
    ! shape functions axixiociated with a 8-nodes serendip
    ! element for a given point of coordinates (xi,eta).
    !
    ! Topology is defined as follows 
    !
    ! 7 - - - 6 - - - 5
    ! |       ^       |
    ! |   eta |       |
    ! |       |       |
    ! 8        --->   4
    ! |        xi     |
    ! |               |
    ! |               |
    ! 1 - - - 2 - - - 3
    !
    
    real(kind=dp), intent(in) :: xil, etal

    real(kind=dp), dimension(8), intent(out) :: shp

    real(kind=dp) :: xip, xim, etap, etam, xixi, etaeta
    
    shp(:) = zero

    xip    = one +  xil
    xim    = one -  xil 
    etap   = one + etal
    etam   = one - etal
    xixi   =  xil *  xil 
    etaeta = etal * etal 
    
    ! Corners first:
    shp(1) = quart * xim * etam * (xim + etam - three)
    shp(3) = quart * xip * etam * (xip + etam - three)
    shp(5) = quart * xip * etap * (xip + etap - three)
    shp(7) = quart * xim * etap * (xim + etap - three)
    
    ! Then midpoints:
    shp(2) = half  * etam * (one -   xixi)
    shp(4) = half  *  xip * (one - etaeta)
    shp(6) = half  * etap * (one -   xixi)
    shp(8) = half  *  xim * (one - etaeta)      
    
  end subroutine shp8
!-----------------------------------------------------------------------------------------
  
!================================================================================
! Derivative of shape functions for quad8
  subroutine shp8der(xil,etal,shpder)
    !
    ! This routine computes and returns the derivatives
    ! of the shape functions axixiociated with a 8-nodes serendip
    ! element for a given point of coordinates (xi,eta).
    !
    ! Topology is defined as follows
    !
    ! 7 - - - 6 - - - 5
    ! |       ^       |
    ! |   eta |       |
    ! |       |       |
    ! 8        --->   4
    ! |        xi     |
    ! |               |
    ! |               |
    ! 1 - - - 2 - - - 3
    !
    !
    ! shpder(:,1) : derivative wrt xi
    ! shpder(:,2) : derivative wrt eta
    
    real(kind=dp), intent(in)                  :: xil, etal
    real(kind=dp), dimension(8,2), intent(out) :: shpder
    real(kind=dp) :: xip, xim, etap, etam, xixi, etaeta
    
    shpder(:,:) = zero
    
    xip    = one +  xil
    xim    = one -  xil
    etap   = one + etal
    etam   = one - etal
    xixi   =  xil *  xil
    etaeta = etal * etal
    
    ! Corners first:
    shpder(1,1) = -quart * etam * ( xim + xim + etam - three)
    shpder(1,2) = -quart *  xim * (etam + xim + etam - three)
    shpder(3,1) =  quart * etam * ( xip + xip + etam - three)
    shpder(3,2) = -quart *  xip * (etam + xip + etam - three)
    shpder(5,1) =  quart * etap * ( xip + xip + etap - three)
    shpder(5,2) =  quart *  xip * (etap + xip + etap - three)
    shpder(7,1) = -quart * etap * ( xim + xim + etap - three)
    shpder(7,2) =  quart *  xim * (etap + xim + etap - three)
    
    ! Then midside points :
    shpder(2,1) = -one  * xil * etam
    shpder(2,2) = -half * (one - xixi)
    shpder(4,1) =  half * (one - etaeta)
    shpder(4,2) = -one  * etal * xip
    shpder(6,1) = -one  * xil * etap
    shpder(6,2) =  half * (one - xixi)
    shpder(8,1) = -half * (one - etaeta)
    shpder(8,2) = -one  * etal * xim
    
  end subroutine shp8der
!-----------------------------------------------------------------------------------------
  
!================================================================================
! Get determinant of jacobian
  real(kind=dp) function det_jacobian_shape(xil, etal, nodes_crd)
    ! This routines the value of the Jacobian (that is, 
    ! the determinant of the Jacobian matrix), for any point
    ! inside a given element. IT ASSUMES 8 nodes 2D isoparametric
    ! formulation of the geometrical transformation and therefore
    ! requires the knowledge of the coordinated of the 8 control
    ! points, which are defined as follows :
    !
    !     7 - - - 6 - - - 5
    !     |       ^       |
    !     |   eta |       |
    !     |       |       |
    !     8        --->   4
    !     |        xi     |
    !     |               |
    !     |               |
    !     1 - - - 2 - - - 3 .
    !
    real(kind=dp), intent(in) ::xil, etal
    real(kind=dp), dimension(8,2), intent(in) :: nodes_crd

    real(kind=dp), dimension(8,2) :: shpder
    integer(kind=si) :: inode
    
    ! Compute the appropriate derivatives of the shape
    ! functions
    call shp8der(xil,etal,shpder)
    
    a = zero
    b = zero
    c = zero
    d = zero
    
    do inode = 1, 8
       a = a + nodes_crd(inode,1)*shpder(inode,1)          
       d = d + nodes_crd(inode,2)*shpder(inode,2)          
       b = b + nodes_crd(inode,1)*shpder(inode,2)          
       c = c + nodes_crd(inode,2)*shpder(inode,1)          
    end do
    
    det_jacobian_shape = a*d - b*c 
    
  end function det_jacobian_shape
!--------------------------------------------------------------------------------  


!================================================================================
! Lagrange interpolants
  subroutine lagrange_any(xi,NGLL,xigll,h,hprime)

    ! subroutine to compute the Lagrange interpolants based upon the GLL points
    ! and their first derivatives at any point xi in [-1,1]
    
    implicit none
    
    integer(kind=si), intent(in) :: NGLL
    real(kind=dp),    intent(in) :: xi
    real(kind=dp), dimension(NGLL), intent(in)  :: xigll
    real(kind=dp), dimension(NGLL), intent(out) :: h, hprime
    
    integer(kind=si) :: dgr, i, j
    real(kind=cp)    :: prod1, prod2
    
    do dgr=1,NGLL
       
       prod1 = 1.0_dp
       prod2 = 1.0_dp
       do i=1,NGLL
          if(i /= dgr) then
             prod1 = prod1*(xi-xigll(i))
             prod2 = prod2*(xigll(dgr)-xigll(i))
          endif
       enddo
       h(dgr)=prod1/prod2
       
       hprime(dgr)=0.0_dp
       do i=1,NGLL
          if(i /= dgr) then
             prod1=1.0_dp
             do j=1,NGLL
                if(j /= dgr .and. j /= i) prod1 = prod1*(xi-xigll(j))
             enddo
             hprime(dgr) = hprime(dgr)+prod1
          endif
       enddo
       hprime(dgr) = hprime(dgr)/prod2
       
    enddo
    
  end subroutine lagrange_any
!--------------------------------------------------------------------------------


!================================================================================
! Find connections between points and elements/gll
  subroutine determine_connectivity

    use global_parameters_mod
    
    real(kind=dp) :: xi, eta, scur, zcur
    real(kind=dp) :: smin, smax, zmin, zmax

    real(kind=dp), dimension(NGNOD,2) :: nodes_crd
    real(kind=dp), parameter          :: eps=1e-3_dp  !!-3

    integer(kind=si), dimension(8) :: IGRIDs, IGRIDz
    integer(kind=si) :: irec, iel, inode, icheck
      
    ! conversion GLL point to 8-node control elements
    IGRIDs(1)=0
    IGRIDs(2)=2
    IGRIDs(3)=4
    IGRIDs(4)=4
    IGRIDs(5)=4
    IGRIDs(6)=2
    IGRIDs(7)=0
    IGRIDs(8)=0
    
    IGRIDz(1)=0
    IGRIDz(2)=0
    IGRIDz(3)=0
    IGRIDz(4)=2
    IGRIDz(5)=4
    IGRIDz(6)=4
    IGRIDz(7)=4
    IGRIDz(8)=2
   
    if (allocated(xi_rec)) deallocate(xi_rec)
    if (allocated(eta_rec)) deallocate(eta_rec)
    if (allocated(rec2elm)) deallocate(rec2elm)         
    if (allocated(rec2elm2)) deallocate(rec2elm2)         
    allocate(xi_rec(nbrec),eta_rec(nbrec))
    allocate(rec2elm(nbrec))
    allocate(rec2elm2(nbrec))
    rec2elm=-1
    rec2elm2=-1

    PRINT *,'THERIS IS :',NBREC
        print *,'nelem : ',nel    
    !*** CONNECTION POINT <-> MESH------------
    do irec=1,nbrec
       scur=reciever_cyl(1,irec)
       zcur=reciever_cyl(3,irec)
       do iel = 1, NEL
          !*** Element
       icheck = 0
          smin=1e40_dp
          smax=-1e40_dp    !-1e25_cp
          zmin=1e40_dp  !smin
          zmax=-1e40_dp !smax
          do inode=1,NGNOD
             nodes_crd(inode,1)=scoor(IGRIDs(inode),IGRIDz(inode),iel)
             nodes_crd(inode,2)=zcoor(IGRIDs(inode),IGRIDz(inode),iel)
             smin=min(smin, nodes_crd(inode,1))
             smax=max(smax, nodes_crd(inode,1))
             zmin=min(zmin, nodes_crd(inode,2))
             zmax=max(zmax, nodes_crd(inode,2))
          end do

          if ( scur > smin-eps .and. scur < smax + eps .and. zcur > zmin-eps .and. zcur < zmax + eps) then
             call find_xix_eta(nodes_crd,xi,eta,scur,zcur)
             rec2elm2(irec)=iel
             if (xi > -1.05 .and. xi < 1.05 .and. eta > -1.05 .and. eta < 1.05) then
                rec2elm(irec)=iel
                xi_rec(irec)=xi
                eta_rec(irec)=eta
                exit
             end if
          else
                icheck = 1
          end if
       end do
       if (icheck == 1) then 
           write(6,*)scur,zcur
       end if
    end do
    call check_rec2elm
    ! END CONNECTION POINT <-> MESH ----
    
  contains
    
    !==================================================
    ! Check if points have been ommited
    subroutine check_rec2elm()
            
      integer :: irec,forgot_point,forgot_point2
      
      forgot_point=0
      FORGOT_POINT2=0
      do irec=1,nbrec
         if(rec2elm(irec)==-1) then
            forgot_point= forgot_point+1
         end if
         if(rec2elm2(irec)==-1) then
            forgot_point2= forgot_point2+1
         end if
      end do
      if (forgot_point > 0) write(*,*) 'forgot ', forgot_point2,' points'
      if (forgot_point > 0) write(*,*) 'forgot ', forgot_point,' points'
      if (forgot_point > 0) stop
    end subroutine check_rec2elm
    !--------------------------------------------------
    
  end subroutine determine_connectivity
!--------------------------------------------------------------------------------

end module interp_mesh_mod
