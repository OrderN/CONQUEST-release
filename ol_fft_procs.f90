!!****h* Conquest/fft_procedures
!!  NAME
!!   fft_procedures
!!  PURPOSE
!!   Contains FFT algorithms for Bessel Transforms
!!  AUTHOR
!!   Numerical Recipes (R Choudhury)
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2007/01/10 12:18 dave
!!    Tidied and enforced f90 conventions
!!  SOURCE
!!

module fft_procedures

  implicit none

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"
!!***

contains

!!****f* fft_procedures/four1 *
!!
!!  NAME 
!!   four1
!!  USAGE
!!   four1(data,nn,isign)
!!  PURPOSE
!!   Evaluates complex valued FFT of input data array
!!  INPUTS
!!   data - array of complex function values
!!   nn - no of real values in data array (half size of array)
!!   isign - sign of FFT (1 for forward, -1 for reverse)
!!  USES
!!   datatypes
!!  AUTHOR
!!   Numerical Recipes (R Choudhury)
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine four1(data,nn,isign)

    use datatypes
    use numbers, ONLY: zero, half, one, two, twopi

    implicit none

    ! Danielson-Lanczos algorithm
    ! for FFT of complex array 'data' dimension nn (equivalent to
    ! a real array of dimension 2*nn

    integer :: n,i,j,m,mmax,istep,isign,nn
    real(double) ::  wr,wi,wpr,wpi,wtemp,theta
    real(double) :: tempr,tempi
    real(double), dimension(2*nn) :: data

    n = 2*nn
    j = 1
    do i=1,n,2
       if(j>i) then
          tempr = data(j)
          tempi = data(j+1)
          data(j) = data(i)
          data(j+1) = data(i+1)
          data(i) = tempr
          data(i+1) = tempi
       endif
       m = n/2
       do while((m>=2).and.(j>m))
          j = j-m
          m = m/2
       enddo
       j = j+m
    end do
    mmax = 2

    do while(n>mmax)
       istep = 2*mmax
       theta = twopi/(real(isign*mmax,double))
       wpr = -two*sin(half*theta)**2
       wpi = sin(theta)
       wr = one
       wi = zero
       do m=1,mmax,2
          do i=m,n,istep
             j = i+mmax
             !tempr = sngl(wr)*data(j)-sngl(wi)*data(j+1)
             !tempi = sngl(wr)*data(j+1)+sngl(wi)*data(j)
             tempr = wr*data(j)-wi*data(j+1)
             tempi = wr*data(j+1)+wi*data(j)
             data(j) = data(i) - tempr
             data(j+1) = data(i+1)-tempi
             data(i) = data(i) + tempr
             data(i+1) = data(i+1) + tempi
          end do
          wtemp = wr
          wr = wr*wpr-wi*wpi+wr
          wi = wi*wpr+wtemp*wpi+wi
       enddo
       mmax = istep
    end do
    return

  end subroutine four1
!!***

!!****f* fft_procedures/realft *
!!
!!  NAME 
!!   realft
!!  USAGE
!!   realft(data,n,isign)
!!  PURPOSE
!!   creates FFT of 2n real valued data points
!!  INPUTS
!!   data - array of real valued data points
!!   n - data is of dimension (2n)
!!   isign - +1 for forward FFT, -1 for reverse
!!  USES
!!   datatypes
!!  AUTHOR
!!   Numerical Recipes ( R Choudhury)
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine realft(data,n,isign)

    use datatypes
    use numbers, ONLY: zero, half, one, two, pi
    implicit none

    ! Fourier transforms set of 2n real valued data points
    ! N must be a power of 2

    real(double) :: wr,wi,wpr,wpi,wtemp,theta
    real(double) :: wrs,wis,h1r,h2r,h1i,h2i,c1,c2
    integer isign,n2p3,i,n,i1,i2,i3,i4
    real(double), dimension(2*n) :: data

    theta = pi/real(n,double)
    c1 = half
    if(isign.eq.1) then
       c2 = -half
       call four1(data,n,+1)
    else
       c2 = half
       theta = -theta
    endif
    wpr = -two*sin(half*theta)**2
    wpi = sin(theta)
    wr = one+wpr
    wi = wpi
    n2p3 = 2*n+3
    do i=2,n/2
       i1 = 2*i-1
       i2 = i1+1
       i3 = n2p3-i2
       i4 = i3+1
       wrs = wr!sngl(wr)
       wis = wi!sngl(wi)
       h1r = c1*(data(i1)+data(i3))
       h1i = c1*(data(i2)-data(i4))
       h2r = -c2*(data(i2)+data(i4))
       h2i = c2*(data(i1)-data(i3))
       data(i1) = h1r+wrs*h2r-wis*h2i
       data(i2) = h1i+wrs*h2i+wis*h2r
       data(i3) = h1r-wrs*h2r+wis*h2i
       data(i4) = -h1i+wrs*h2i+wis*h2r
       wtemp = wr
       wr = wr*wpr-wi*wpi+wr
       wi = wi*wpr+wtemp*wpi+wi
    end do
    if(isign==1) then
       h1r = data(1)
       data(1) = h1r+data(2)
       data(2) = h1r-data(2)
    else
       h1r = data(1)
       data(1) = c1*(h1r+data(2))
       data(2) = c1*(h1r-data(2))
       call four1(data,n,-1)
    endif
    return
  end subroutine realft
!!***

!!****f* fft_procedures/sinft *
!!
!!  NAME 
!!   sinft
!!  USAGE
!!   sinft(y,n)
!!  PURPOSE
!!   Evaluate sine transform of data
!!  INPUTS
!!   y - array of data
!!   n - dimension of y
!!  USES
!!   datatypes
!!  AUTHOR
!!   Numerical Recipes (R Choudhury)
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!

  subroutine sinft(y,n)

    use datatypes
    use numbers, ONLY: zero, half, one, two, pi

    implicit none

    integer, intent(in) :: n
    real(double), dimension(n), intent(out) :: y

    integer :: j,m
    real(double) wr,wi,wpr,wpi,wtemp,theta
    real(double) y1,y2,sum

    theta = pi/real(n,double)
    wr = one
    wi = zero
    wpr = -two*sin(half*theta)**2
    wpi = sin(theta)
    y(1) = zero
    m = n/2 !not in new num recps
    do j=1,m
       wtemp = wr
       wr = wr*wpr-wi*wpi+wr
       wi = wi*wpr+wtemp*wpi+wi
       y1 = wi*(y(j+1)+y(n-j+1))
       y2 = half*(y(j+1)-y(n-j+1))
       y(j+1) = y1+y2
       y(n-j+1) = y1-y2
    end do
    call realft(y,m,+1)
    sum = zero
    y(1) = half*y(1)
    y(2) = zero
    do j=1,n-1,2
       sum = sum+y(j)
       y(j) = y(j+1)
       y(j+1) = sum
    end do
    return
  end subroutine sinft
!!***
       
!!****f* fft_procedures/cosft *
!!
!!  NAME 
!!   cosft
!!  USAGE
!!   cosft(y,n,isign)
!!  PURPOSE
!!   Evaluate cosine transform of data array
!!  INPUTS
!!   y - data array
!!   n - dimension of data array
!!   isign : +1 for forward transform, -1 for reverse
!!  USES
!!   datatypes
!!  AUTHOR
!!   Numerical REcipes, (R Choudhury)
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2008/03/03 18:47 dave
!!    Changed float to real
!!   2012/06/17 L.Tong
!!   - Changed y to intent(inout), it was intent(out) which was incorrect
!!  SOURCE
!!
  subroutine cosft(y,n,isign)

    use datatypes
    use numbers, ONLY: zero, half, one, two, pi

    implicit none

    integer, intent(in) :: n,isign
    real(double), dimension(n), intent(inout) :: y

    real(double) wr,wi,wpr,wpi,wtemp,theta
    real(double) sum,y1,y2,even,odd,enfo,sumo,sume
    integer i,m,j

    theta = pi/real(n,double)
    wr = one
    wi = zero
    wpr = -two*sin(half*theta)**2
    wpi = sin(theta)
    sum = y(1)
    m = n/2
    do j=1,m-1
       wtemp = wr
       wr = wr*wpr-wi*wpi+wr
       wi = wi*wpr+wtemp*wpi+wi
       y1 = half*(y(j+1)+y(n-j+1))
       y2 = (y(j+1)-y(n-j+1))
       y(j+1) = y1-wi*y2
       y(n-j+1) = y1+wi*y2
       sum = sum+wr*y2
    enddo
    call realft(y,m,+1)
    y(2) = sum
    do j=4,n,2
       sum = sum+y(j)
       y(j) = sum
    end do
    if(isign.eq.-1) then
       even = y(1)
       odd = y(2)
       do i=3,n-1,2
          even = even+y(i)
          odd = odd+y(i+1)
       end do
       enfo = two*(even-odd)
       sumo = y(1)-enfo
       sume = (two*odd/real(n,double))-sumo
       y(1) = half*enfo
       y(2) = y(2)-sume
       do i=3,n-1,2
          y(i) = y(i) -sumo
          y(i+1) = y(i+1)-sume
       end do
    endif
    return
  end subroutine cosft
!!***


 end module fft_procedures
!!***
 
