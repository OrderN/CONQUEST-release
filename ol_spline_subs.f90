! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: ol_spline_subs.f90,v 1.4 2005/05/26 08:36:27 drb Exp $
! ------------------------------------------------------------------------------
! Module cubic_spline_routines
! ------------------------------------------------------------------------------
! Code area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/cubic_spline_routines
!!  NAME
!!   cubic_spline_routines
!!  PURPOSE
!!   contains routines to do cubic spline interpolations
!!   and also to reproject an array of given grid spacing onto
!!   a different grid spacing
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!
module cubic_spline_routines
   
  implicit none
  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id: ol_spline_subs.f90,v 1.4 2005/05/26 08:36:27 drb Exp $"

!!***
contains

!!****f* cubic_spline_routines/spline_new *
!!
!!  NAME 
!!   spline_new
!!  USAGE
!!   spline_new(x,y,n,yp1,ypn,y2)
!!  PURPOSE
!!   creates array of second derivatives of the input argument function
!!  INPUTS
!!   x - abscissae of argument function
!!   y - values of function at abscissae
!!   n - no of points
!!   yp1 - 1st derivative at first x
!!   ypn - 1st derivative at last x
!!   y2 - array of second derivatives
!!  USES
!!   datatypes
!!  AUTHOR
!!   Numerical Recipes
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine spline_new(x,y,n,yp1,ypn,y2)

    use datatypes

    implicit none

    integer, intent(in) ::  n
    real(double), intent(in) ::  yp1,ypn,x(n),y(n)
    real(double), intent(out) :: y2(n)

    !takes arrays x(1:n) and y(1:n) and outputs
    !y2(1:n), the array of 2nd derivatives at
    !the specified abscissae.  yp1,ypn and values of 2nd
    !derivative at 1st and last point respectively
    integer i,k,stat
    real(double) :: p,qn,sig,un
    real(double), allocatable :: u(:)

    allocate(u(n), STAT=stat)
    if (yp1.gt..99e30_double) then
       y2(1) = 0.0_double
       u(1) = 0.0_double
    else
       y2(1) = -0.5_double
       u(1) = (3.0_double/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do  i=2,n-1
       sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p = sig*y2(i-1)+2._double
       y2(i) = (sig-1.0_double)/p
       u(i) = (6.0_double*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
            &/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if(ypn.gt..99e30) then
       qn = 0.0_double
       un = 0.0_double
    else
       qn = 0.5_double
       un = (3.0_double/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_double)
    do  k=n-1,1,-1
       y2(k) = y2(k)*y2(k+1)+u(k)
    enddo
    deallocate(u)
    return

  end subroutine spline_new
  !!***

!!****f* cubic_spline_routines/splint_new2 *
!!
!!  NAME 
!!   splint_new2
!!  USAGE
!!   splint_new2(xa,ya,y2a,n,x,y)
!!  PURPOSE
!!   outputs spline interpolated value of a function
!!   at some specified argument x - uses fact that grid is uniform
!!  INPUTS
!!   xa - array of abscissae
!!   ya - array of function values
!!   y2a - array of second derivatives
!!   n - no of points in arrays
!!   x - value at which function is to be interpolated
!!   y - value of interpolate function at x.
!!  USES
!!   datatypes
!!  AUTHOR
!!   Numerical Recipes\R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2007/01/13 09:14 dave
!!    Tweaked and tidied
!!  SOURCE
!!
  subroutine splint_new2(xa,ya,y2a,n,x,y)

    use datatypes
    use numbers, ONLY: very_small, one, six
    use GenComms, ONLY: cq_abort

    implicit none

    integer, intent(in) :: n
    real(double), intent(in) :: x,xa(n),y2a(n),ya(n)
    real(double), intent(out) :: y

    integer k,khi,klo
    real(double) a,b,h,dx

    dx = xa(n)/real(n-1,double)
    klo = aint(x/dx)+1
    khi = klo+1

    h = xa(khi)-xa(klo)
    if (h<very_small) call cq_abort("bad xa input in splint_new")
    a=(xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)+(a*(a*a-one)*y2a(klo)+b*(b*b-one)*y2a(khi) )*(h*h)/six
    return

  end subroutine splint_new2
!!***

!!****f* cubic_spline_routines/reproject_array *
!!
!!  NAME 
!!   reproject_array
!!  USAGE
!!   reproject_array(yin,del_old,n_old,yout,del_new,n_new)
!!  PURPOSE
!!   Reprojects an array from a given grid spacing onto new grid spacing
!!  INPUTS
!!   yin - array to be reprojected
!!   del_old - old grid spacing
!!   n_old - old no of grid points
!!   yout - reprojected array values
!!   del_new - new array grid spacing
!!   n_new - new no of grid points
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine reproject_array(yin,del_old,n_old,yout,del_new,n_new)

    use datatypes

    implicit none

    integer, intent(in) :: n_old,n_new
    real(double), intent(in) :: del_old,del_new
    real(double), intent(in), dimension(n_old) :: yin
    real(double), intent(out), dimension(n_new) :: yout

    integer :: i,l,m
    real(double), allocatable, dimension(:) :: y2,xin
    real(double) :: yp1,ypn, x,y

    allocate(xin(n_old))
    allocate(y2(n_old))
    !setting up x grid for spline_new subroutine
    do i=1,n_old
       xin(i) = (i-1)*del_old
    enddo

    !setting first derivatives of endpoints for spline
    !nb this DOES introduce an inaccuracy (depending on grid fineness)
    yp1 = (yin(2)-yin(1))/del_old
    ypn = (yin(n_old)-yin(n_old-1))/del_old

    call spline_new(xin,yin,n_old,yp1,ypn,y2)

    do i=1,n_new
       x=(i-1)*del_new
       call splint_new2(xin,yin,y2,n_old,x,y)
       yout(i) = y
    enddo
    deallocate(xin)
    deallocate(y2)

  end subroutine reproject_array
!!***

!!****f* cubic_spline_routines/matcharrays *
!!
!!  NAME 
!!   matcharrays
!!  USAGE
!!   matcharrays(y1,n1,d1,y2,n2,d2,d12,n1_new,n2_new,n12)
!!
!!  PURPOSE
!!   Takes two arrays and projects them onto (same) specified grid spacing 
!!  INPUTS
!!   y1, y2 - (values of) the two arrays
!!   n1, n2 - no of points in each array
!!   d1, d2 - spacing of original arrays
!!   d12 - new grid spacing
!!   n1_new, n2_new  - no of points in reprojected arrays
!!   n12 - size of reprojected arrays - power of two for FFT algorithms
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine matcharrays(y1,n1,d1,y2,n2,d2,d12,n1_new,n2_new,n12)

    use datatypes

    implicit none

    integer, intent(in) :: n1,n2,n12
    real(double), intent(inout), dimension(n12) :: y1,y2
    real(double), allocatable, dimension(:) :: ytemp
    real(double), intent(in) :: d1,d2,d12
    integer,intent(inout) :: n1_new,n2_new
    integer i

    n1_new = ((n1-1)*d1/d12)+1
    allocate(ytemp(n1_new))
    call reproject_array(y1(1:n1),d1,n1,ytemp,d12,n1_new)
    y1(1:n1_new) = ytemp(1:n1_new)
    y1(n1_new+1:n12) = 0.0_double
    deallocate(ytemp)

    n2_new = ((n2-1)*d2/d12)+1
    allocate(ytemp(n2_new))
    call reproject_array(y2(1:n2),d2,n2,ytemp,d12,n2_new)
    y2(1:n2_new) = ytemp(1:n2_new)
    y2(n2_new+1:n12) = 0.0_double
    deallocate(ytemp)

  end subroutine matcharrays
!!***

!!****f* cubic_spline_routines/spline_ol_intval_new2 *
!!
!!  NAME 
!!   spline_ol_intval_new2
!!  USAGE
!!   spline_ol_intval_new2(r,radial_table,radial_table2,npts,del_x,ol_intval)
!!  PURPOSE
!!   Interpolate from array of overlap integral values
!!   different to routine above since looks up table of 2nd derivatives
!!   in radial_table2
!!   a value at a given radius.
!!  INPUTS
!!   r : radius at which value is required
!!   radial_table : array of overlap integral values at different r
!!   radial_table2 - array of 2nd derivatives of radial table
!!   ol_intval : (interpolated) overlap integral value at r
!!   npts,del_x, no of tabulated data points, spacing between points.
!!  USES
!!    datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine spline_ol_intval_new2(r,radial_table,radial_table2,npts,del_x,ol_intval)
    use datatypes
    implicit none
    !routine to perform cubic spline interpolation of value
    !of radial_table_out array for pre-specified R
    !RC 18/09/03 taking into account that y2 tables are initialised already
    real(double), intent(in) :: r
    real(double), intent(out) :: ol_intval
    real(double) :: yp1,ypn,xj1,xj2,a,b,c,d
    real(double) :: del_x
    integer :: i,j,npts,n,p,q,flag_spline_direct,n1,n2
    real(double), intent(in), dimension(npts) :: radial_table,radial_table2
    real(double), allocatable, dimension(:) :: y_points
    real(double), allocatable, dimension(:) :: x_vals
    real(double), allocatable, dimension(:) :: y2
    !RC recoding this routine to do a direct spline interpolation on the 
    !input radial table and its' second derivatives - 04/11/03
    !have to make sure I have enough space to take n-3 points
    flag_spline_direct = 1
    if(flag_spline_direct.eq.1) then
       !          write(*,*) 'using direct spline interpolation'
       !locate points nearest to r
       n1 = aint(r/del_x)
       n2 = n1+1
       if(n2>npts-1) then
          !write(*,*) 'OVERFLOW ERROR ',r,del_x,npts,n2
          ol_intval = 0.0_double
          return
       end if
       xj1 = n1*del_x
       xj2 = xj1+del_x

       a = (xj2-r)/del_x
       b = 1.0_double-a
       c = a*(a*a-1.0_double)*del_x*del_x/6.0_double
       d = b*(b*b-1.0_double)*del_x*del_x/6.0_double

       ol_intval = a*radial_table(n1+1)+b*radial_table(n2+1)+&
            &c*radial_table2(n1+1) + d*radial_table2(n2+1)
       !write(*,*) a*radial_table(n1+1)+b*radial_table(n2+1), 'spline first part'
       !          write(*,*) ol_intval, 'spline returned value'
    else
       return
    endif
  end subroutine spline_ol_intval_new2
 !!***
 
end module cubic_spline_routines

       
       
  
