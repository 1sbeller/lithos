module interp_process_mod

  use precision_mod
  use constants_mod
  use interpolation_parameters_mod

  implicit none

contains

!================================================================================
! Compute stalta ratio and give first pick
  subroutine substalta(sig,n,nsta,nlta,crit,stalta,tpick)

    integer(kind=si), intent(in)            :: n, nsta, nlta
    real(kind=cp), intent(in)               :: crit
    real(kind=cp), dimension(n), intent(in) :: sig

    integer(kind=si), intent(out) :: tpick

    real(kind=cp), dimension(n)    :: sta, lta
    real(kind=cp), dimension(n), intent(out) :: stalta
    real(kind=cp), dimension(n+2*nsta) :: tmpsta
    real(kind=cp), dimension(n+2*nlta) :: tmplta

    integer(kind=si) :: i

    !*** Compute the short time average (STA)
    tmpsta(1:nsta) = sig(1)
    tmpsta(nsta+1:nsta+n) = sig(:) 
    tmpsta(nsta+n+1:n+2*nsta) = sig(n)
    sta = zero
    do i=1+nsta,n+nsta
       sta(i-nsta) = sum(tmpsta(i-nsta:i+nsta)**2)
    end do
    sta = 0.5 * sta / nsta

    !*** Compute the long time average (LTA)
    tmplta(1:nlta) = sig(1)
    tmplta(nlta+1:nlta+n) = sig(:) 
    tmplta(nlta+n+1:n+2*nlta) = sig(n)
    lta = zero
    do i=1+nlta,n+nlta
        lta(i-nlta) = sum(tmplta(i-nlta:i+nlta)**2)
    end do
    lta = 0.5 * lta / nlta


    !*** Compute ratio and gives first pick
    stalta = sta / lta
    do i=1,n
      if (stalta(i) >= crit) then
        tpick = i
        exit
      end if
    end do

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

    !*** Take good parts (this is wrong...)
    if (modulo(n2,2) == 0) then
        ind = n2/2+1
    else
        ind = ceiling(real(n2/2,kind=cp))
    end if
    ind = 1
    conv(1:n1) = convtmp(ind:ind+n1-1)

!	print *,'Indice ',ind,'size ',n1+n2-1
  
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
  subroutine bwfilt (x, y, dt, n, irek, norder, f1, f2)

    ! recursive filtering of data with butterworth filter
    ! x: input array
    ! y: output array
    ! dt: time increment
    ! n: number of data points
    
    ! irek=0: forward filtering only
    ! irek=1: forward and backward filtering
    
    ! norder: order of butterworth filter
    ! norder=0: only filtering, no determination of coefficients
    ! norder<0: no starplots of transfer function and impulse response
    
    ! f1: low cutoff frequency (Hz)
    ! f1=0: low pass filter
    
    ! f2: high cutoff frequency (Hz)
    ! f2>0.5/dt: high pass filter
    
    implicit none

    integer(kind=si) :: iunit, npoles,n,lx    
    integer(kind=si) :: irek,norder    
    real(kind=cp), dimension(n) ::x,y
    real(kind=cp), dimension (10) ::  a, b1, b2
    real(kind=cp) :: dt,f1,f2

    !real(kind(0d0)) :: x(n),y(n)
    
    iunit = 3
    
    if(norder/=0) then
       npoles=abs(norder)
       !determination of filter coefficients
       call bpcoeff(f1,f2,npoles, dt, a,b1, b2)
       if(norder>=0) then
          !plot of transfer function and impuulse response
          lx = 100
          !filtering
       endif
    endif
    
    
    if(n/=0) then
       call rekurs(x,y,n,a,b1,b2,npoles,irek)
    endif
    return
  end subroutine bwfilt
  
  !---------------------------------------------------------------
  
  
  subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
    ! performs recursive filtering of data in array x of length ndat
    ! filtered output in y
    ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
    ! npoles is the number of poles
    ! iflag=0: forward filtering only
    ! iflag/=0: forward and backward filtering
    
    implicit none
    
    real(kind=cp), dimension(10) :: z,z1,z2 ,a,b1,b2
    real(kind=cp)  ::  x1,x2
    integer(kind=si) :: ndat, npoles,n,i
    integer(kind=si) :: iflag
    real(kind=cp), dimension(ndat) :: x, y
    
    !forward
    
    x1 = 0.d0
    x2 = 0.d0
    
    do i = 1, npoles
       z1(i) = 0.d0
       z2(i) = 0.d0
    enddo
    
    do n = 1, ndat
       z(1) = a(1)*(x(n)-x2) -b1(1)*z1(1) -b2(1)*z2(1)
       do i = 2, npoles
          z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
       enddo
       x2=x1
       x1=x(n)
       do i = 1, npoles
          z2(i) =z1(i)
          z1(i) =z(i)
       enddo
       y(n) = z(npoles)
    enddo
    
    if(iflag==0) then
       return
    endif
    
    !backward
    
    x1 =0.d0
    x2 =0.d0
    
    do i = 1, npoles
       z1(i) = 0.d0
       z2(i) = 0.d0
    enddo

    do n = ndat, 1, -1
       z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
       do i =2, npoles
          z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
       enddo
       x2=x1
       x1=y(n)
       do i = 1,npoles
          z2(i)=z1(i)
          z1(i)=z(i)
       enddo
       y(n) = z(npoles)
    enddo
    return
  end subroutine rekurs

  subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
    !determines filtercoefficients for recursive bandpassfilter
    
    real(kind=cp),dimension(10) :: a,b1,b2
    complex(kind=4) :: s(20), t1,t2,p
    real(kind=cp), parameter :: pii = 3.141592653589793d0
    real(kind=cp) :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
    integer(kind=si) :: i,npol2,n,npoles
    
    
    if(npoles>10) then
       stop ' npoles greater than 10: STOP '
    endif
    
    d2= 2.d0/dt
    w1=d2*tan(2.d0*pii*f1/d2)
    w2=d2*tan(2.d0*pii*f2/d2)
    w0=0.5*(w2-w1)
    
    i=1
    npol2=npoles/2+1
    do n =1,npoles
       p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pii/dble(2*npoles)))
       t1 = p*cmplx(w0,0.d0)
       t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
       s(i)=t1+t2
       s(i+1)=t1-t2
       i=i+2
    enddo
    
    do n=1,npoles
       ssum=2*real(s(n))
       sprod=dble(s(n)*conjg(s(n)))
       fact1=d2*d2-d2*ssum+sprod
       fact2=2.d0*(sprod-d2*d2)
       fact3=d2*d2+d2*ssum+sprod
       a(n)=2.d0*d2*w0/fact1
       b1(n)=fact2/fact1
       b2(n)=fact3/fact1
    enddo
    return
  end subroutine bpcoeff


end module interp_process_mod
