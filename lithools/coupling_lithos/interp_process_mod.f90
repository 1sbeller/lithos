module interp_process_mod

  use precision_mod
  use constants_mod
  use interpolatation_parameters_mod

contains

  
!================================================================================
! Compute STA/LTA for picking
  subroutine substalta(sig, nsta, nlta, stalta)

    real(kind=cp), dimension(:), allocatable, intent(in)  :: sig
    real(kind=cp), dimension(:), allocatable, intent(inout) :: stalta

    integer(kind=si), intent(in) :: nsta, nlta

    integer(kind=si) :: m, nsta_1, nlta_1, i

    real(kind=cp), dimension(:), allocatable :: sta, lta, pad_sta, pad_lta, tmp1, tmp2

    m = size(sig)

    nsta_1 = nsta - 1
    nlta_1 = nlta - 1


    allocate(sta(m))
    allocate(lta(m))
    allocate(tmp1(m))
    allocate(tmp2(m))
    allocate(pad_sta(nsta_1))
    allocate(pad_lta(nlta_1))
    sta = 0.
    lta = 0.
    pad_sta = 0.
    pad_lta = 1.

    !*** compute the short time average (STA)
    do i=1,nsta
       tmp1(1:nsta_1) = pad_sta(:)
       tmp1(nsta_1+1:m) = sig(i:m - nsta_1 + i-1)**2
       sta = sta + tmp1
    end do
    sta = sta / nsta

    !*** compute the long time average (LTA)
    do i =1,nlta
       tmp2(1:nlta_1) = pad_lta(:)
       tmp2(nlta_1+1:m) = sig(i:m - nlta_1 + i-1)**2
       lta = lta + tmp2
    end do
    lta = lta / nlta

!    sta(1:nsta_1) = 0.
!    lta(1:nlta_1) = 1.

    do i=1,m
       if (lta(i) < 1e-20) then
          lta(i) = 1.
          sta(i) = 0.
       end if
    end do
    stalta = sta / lta

  end subroutine substalta
!--------------------------------------------------------------------------------

!================================================================================
! Tukey tapering windows
!--------------------------------------------------
! N     : number of samples
! alpha : percentage of signal to taper
! tuk   : tapered window
  function tuckeywin(N,alpha) result(tuk)

    integer(kind=si), intent(in)   :: N
    real(kind=cp), intent(in)      :: alpha

    integer(kind=si) :: i
    real(kind=cp), parameter :: pipi=3.141592653589793
    real(kind=cp), dimension(N) :: tuk

    !*** Central part
    tuk(:) = 1.

    !*** Left part
    do i=0,int(0.5*alpha*(N-1))
       tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-1.)))
    end do


    !*** Right part
    do i=int((N-1)*(1-alpha/2.)),N-1
       tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-(2./alpha)+1.)))
    end do

  end function tuckeywin
!--------------------------------------------------------------------------------

!================================================================================
! Sinc functions
  real(kind=cp) function mysinc(x)

    real(kind=cp) :: x
    real(kind=cp), parameter :: pipi=3.141592653589793

    if (abs(x) >= 1e-13) then
       mysinc = sin(pipi*x)/(pipi*x)
    else
       mysinc = 1
    end if

  end function mysinc
!--------------------------------------------------------------------------------

  subroutine comp_tab_sinc(itn,dtn,fo,nto,tab_sinc)

    integer(kind=si), intent(in) :: itn, nto
    real(kind=cp) , intent(in) :: fo, dtn
    real(kind=cp), dimension(nto), intent(out) :: tab_sinc  
    integer(kind=si) :: j

    do j=1,nto
       tab_sinc(j) = mysinc(real(fo * itn * dtn - j))
    end do
    
  end subroutine comp_tab_sinc


!================================================================================
! Next open unit
  function next_iunit(i)

    integer(kind=si) ::  i,next_iunit

    i=i+1
    next_iunit=i

  end function next_iunit
!--------------------------------------------------------------------------------

  subroutine prepare_adjoint_filter

  !  real(kind=cp), intent(in) :: fmax, dt
  !  real(kind=cp), dimension(:), allocatable, intent(out) :: Ker, conv, taper
  !  integer(kind=si), intent(out) :: ibeg, iend

    !*** Define bandwidth
    n1 = ceiling(4./(fmax * dt))        !before it was 4./fmax...  !bandwidth              !Let's say the frequency
    if (modulo(n1,2) == 0) then
       n1 = n1+1
    end if

    !*** Allocate
    if (.not. allocated(ker)) allocate(ker(n1))
    if (.not. allocated(taptap)) allocate(taptap(nt))
    if (.not. allocated(conv)) allocate(conv(n1+nt-1))
    ker = 0.
    conv = 0.

    !*** Create filter kernel (may be done elsewhere)
    call sinc_kernel(dt, dble(1./fmax), n1, ker)
    
    !*** Get taper
    taptap = tuckeywin(nt,0.05)

    !*** Index (to get good part of conoluion)
    ibeg = ceiling(n1/2.)-1
    iend = ibeg + nt - 1

  end subroutine prepare_adjoint_filter

!!$  subroutine filter_adjoint_source
!!$    
!!$!    integer(kind=si), intent(in) :: irec
!!$!
!!$!    !*** Loop over receivers
!!$!    do irec=1,nrec
!!$!           
!!$       
!!$    !*** Taper on residuals
!!$    toconv(irec,:) = residu_vx(irec,:) * taptap(:)
!!$    residu_vy(irec,:) = residu_vy(irec,:) * taptap(:)
!!$    residu_vz(irec,:) = residu_vz(irec,:) * taptap(:)
!!$
!!$      
!!$    !*** Filter on residuals
!!$    call myconvolve(dble(residu_vx(irec,:)),ker,nt,n1,conv)
!!$    residu_vx(irec,:) = real(conv(ibeg:iend))
!!$    call myconvolve(dble(residu_vy(irec,:)),ker,nt,n1,conv)
!!$    residu_vy(irec,:) = real(conv(ibeg:iend))
!!$    call myconvolve(dble(residu_vz(irec,:)),ker,nt,n1,conv)
!!$    residu_vz(irec,:) = real(conv(ibeg:iend))
!!$
!!$    end do
!!$
!!$  end subroutine filter_adjoint_source

 subroutine sinc_kernel(dx1in, dx1out, n1, ker)

    real(kind=cp), intent(in) :: dx1in, dx1out
    integer(kind=si), intent(in) :: n1

    real(kind=cp), dimension(n1), intent(out) :: ker

    real(kind=cp), dimension(n1) :: taperx

    real(kind=cp) :: nf1, nf2, cf1, cf2, sum
    integer(kind=si) :: i, j

    !*** Original frequency
    nf1 = 1./dx1in
    
    !*** Cut-off frquancies
    cf1 = 1./dx1out

    !*** Creates sinc kernels
    do i=1,n1
       ker(i) = mysinc(real(2.*(cf1/nf1))*real((i-(n1-1)/2)))
    end do

    !*** Get taper
    taperx = dble(tuckeywin(n1,0.5))

    do i=1,n1
       ker(i)= ker(i) * taperx(i)
    end do
       
    !*** Normalize by sum
    sum = 0.
    do i=1,n1
       sum = sum +  ker(i)
    end do
    
     ker(:) = ker(:) / sum

     print *,'WARNBIIING',cf1,nf1,dx1in,dx1out
     print *,size(ker)
     print *,ker

   end subroutine sinc_kernel


  subroutine myconvolve(a,b,na,nb,conv)

    integer(kind=si), intent(in) :: na, nb
    
    real(kind=cp), dimension(na), intent(in) :: a
    real(kind=cp), dimension(nb), intent(in) :: b

    real(kind=cp), dimension(na+nb-1), intent(out) :: conv
    integer(kind=si) :: ia, ib

    !*** Put to zero
    conv = 0.

    !*** Convolve
    do ia=1,na
      do ib=1,nb
         conv(ia+ib-1) = conv(ia+ib-1) + a(ia) * b(ib)
      end do
   end do

  end subroutine myconvolve


end module interp_process_mod
