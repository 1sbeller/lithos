module interp_process_mod

  use precision_mod
  use constants_mod
  use interpolation_parameters_mod

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
!========================================================================
! Convolution routine
  subroutine myconvolution(sig1,sig2,n1,n2,conv)

    integer(kind=si), intent(in) :: n1, n2

    real(kind=cp), dimension(n1), intent(in) :: sig1
    real(kind=cp), dimension(n2), intent(in) :: sig2

    real(kind=cp), dimension(n1), intent(out) ::conv

    real(kind=cp), dimension(n1+n2-1) :: convtmp !, intent(out) :: conv
    integer(kind=si) :: i1, i2, ind

    !*** Put to zero
    convtmp = zero

    !*** Convolve
    do i1=1,n1
      do i2=1,n2
         convtmp(i1+i2-1) = convtmp(i1+i2-1) + sig1(i1) * sig2(i2)
      end do
    end do

    !*** Take good parts
    if (modulo(n2,2) == 0) then
        ind = n2/2+1
    else
        ind = ceiling(real(n2/2,kind=cp))
    end if
    conv(:) = convtmp(ind:ind+n1-1)
  
  end subroutine myconvolution
!------------------------------------------------------------------------

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


end module interp_process_mod
