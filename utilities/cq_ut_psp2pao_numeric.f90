! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_numeric
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_numeric *
!!
!!NAME
!! cq_ut_psp2pao_numeric
!!PURPOSE
!! Numerical calculations (splines, integrals...)
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 24/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_numeric
  
  use cq_ut_psp2pao_types

  implicit none

  ! Constants
  real(double), parameter :: very_small = 1.0e-30_double
  real(double), parameter :: big = 1.0e12_double

  real(double), parameter :: zero  = 0.0_double
  real(double), parameter :: one   = 1.0_double
  real(double), parameter :: two = 2.000000_double
  real(double), parameter :: three = 3.0_double
  real(double), parameter :: four = 4.000000_double
  real(double), parameter :: six   = 6.0_double
  real(double), parameter :: seven = 7.000000_double
  real(double), parameter :: half = one/two
  real(double), parameter :: third = one/three                 !0.333333333_double        ! 1/3
  real(double), parameter :: seven_thirds = seven/three        !2.333333333_double ! 7/3
  real(double), parameter :: four_thirds = four/three          !4.000000_double / 3.000000_double


contains

! -----------------------------------------------------------------------------
! Subroutine spline
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_numeric/spline *
!!
!!NAME
!! spline
!!USAGE
!!
!!PURPOSE
!! Calculates the second dervatives (y2) of a discrete function (x, y),
!!   i.e. tabulated data, according the the spline formula
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 24/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine spline(x, y, length, d2y,yp1,ypn)

     use cq_ut_psp2pao_types

     implicit none

     ! Passed variables

     integer :: length
     real(double) :: x(length)
     real(double) :: y(length)
     real(double) :: d2y(length)
     real(double), OPTIONAL :: yp1, ypn

     ! Local variables

     integer :: stat, i
     real(double) :: dy1, dyn
     real(double) :: incr, product
     real(double) :: tx1, tx2, ty1, ty2, tmp1, tmp2       ! Temporal variables
     real(double), dimension(:), allocatable :: tmpv      ! Temporal vector


     allocate(tmpv(length), STAT=stat)

     ! All this is basically the solution of a tridiagonal
     !  set of algebraic equations for the second derivatives

     !! Left extreme (boundary condition)
     if(length < 3) then
        ! Use natural spline, i.e. derivatives at extremes = 0

        d2y(1) = 0.0_double
        tmpv(1) = 0.0_double
     else
        ! Use spline with defined slopes at extremes

        if(PRESENT(yp1)) then
           dy1 = yp1
        else
           ! Calculate slopes using quadratic approximation

           tx1 = x(2)-x(1)
           tx2 = x(3)-x(1)
           ty1 = y(2)-y(1)
           ty2 = y(3)-y(1)
           dy1 = (ty1/(tx1*tx1) - ty2/(tx2*tx2))*tx1*tx2/(tx2-tx1)
        end if
        !write(*,*) 'Starting slope: ',dy1
        d2y(1) = -0.5_double
        tmpv(1) = (3.0_double/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1)) - dy1)
     end if

     !! Middle points
     do i=2,length-1
        incr = (x(i) - x(i-1))/(x(i+1) - x(i-1))
        product = incr*d2y(i-1) + 2.0_double

        d2y(i) = (incr - 1.0_double)/product

        tmpv(i) = (6.0_double &
                * ((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1)))&
                / (x(i+1)-x(i-1)) &
                - incr*tmpv(i-1))/product
     end do

     !! Right extreme (boundary condition)
     if(length < 3) then
        ! Use natural spline, i.e. derivatives at extremes = 0

        tmp1 = 0.0_double
        tmp2 = 0.0_double
     else
        ! Use spline with defined slopes at extremes

        ! Calculate slopes using quadratic approximation

        if(PRESENT(ypn)) then
           dyn = ypn
        else
           tx1 = x(length)-x(length-1)
           tx2 = x(length)-x(length-2)
           ty1 = y(length)-y(length-1)
           ty2 = y(length)-y(length-2)
           dyn = (ty1/(tx1*tx1) - ty2/(tx2*tx2))*tx1*tx2/(tx2-tx1)
        endif
        !write(*,*) 'Slope at end is ',dyn,-4/(x(length)*x(length))
        tmp1 = 0.5_double
        tmp2 = (3.0_double/(x(length) - x(length-1))) &
             * (dyn - (y(length) - y(length-1)) / (x(length) - x(length-1)))
     end if

     d2y(length) = (tmp2 - tmp1*tmpv(length-1)) / (tmp1*d2y(length-1) + 1.0_double)

     !! Backsubstitution of the tridiagonal solver (i.e. actual solution)
     do i=length-1, 1, -1
        d2y(i) = d2y(i)*d2y(i+1) + tmpv(i)
     end do

     if(allocated(tmpv)) deallocate(tmpv)

  end subroutine spline
!!***

! -----------------------------------------------------------------------------
! Subroutine spline_derivative
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_numeric/spline_derivative *
!!
!!NAME
!! spline_derivative
!!USAGE
!!
!!PURPOSE
!! Calculates a derivative using a spline
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 06/02/2007
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine spline_derivative(x, y, d2y, length, x_der, y_der)

     use cq_ut_psp2pao_types

     implicit none

     ! Passed variables

     integer :: length
     real(double) :: x(length)
     real(double) :: y(length)
     real(double) :: d2y(length)
     real(double) :: x_der, y_der

     ! Local variables
     integer :: xmin, xhalf, xmax
     real(double) :: delta
     real(double) :: tmp1, tmp2

     if(length == 0) y_der = zero

     if(x_der <= x(1)) then
        delta = (x(2) - x(1))
        y_der = (y(2) - y(1)) / delta - (2.0_double * d2y(1) + d2y(2)) * delta / six;
return
     end if

     if(x_der > x(length)) then
        y_der = zero
     end if

     xmin = 1
     xmax = length
     do while(xmax - xmin > 1)
        xhalf = (xmax + xmin) / 2;
        if(x(xhalf) > x_der) then
          xmax = xhalf;
        else
          xmin = xhalf;
        end if
     end do

     delta = x(xmax) - x(xmin)
     if(delta < very_small) then
       write(*,*) 'Error in spline derivative: Increment = 0'
     end if

     tmp1  = (x(xmax) - x_der) / delta;
     tmp2  = (x_der - x(xmin)) / delta;
     y_der = (y(xmax) - y(xmin)) / delta - ((three * tmp1 * tmp1 - one) * d2y(xmin) &
           - (three * tmp2 * tmp2 - one) * d2y(xmax)) * delta / six;

!print *,"AST-splineder", xmin, xmax, tmp1, tmp2, delta, x_der, y_der

  end subroutine spline_derivative
!!***

! -----------------------------------------------------------------------------
! Subroutine spline_interpolation
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_numeric/spline_interpolation *
!!
!!NAME
!! spline_interpolation
!!USAGE
!!
!!PURPOSE
!! Interpolates a point using a spline
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 27/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine spline_interpolation(x, y, length, d2y, int_x, int_y)

     use cq_ut_psp2pao_types

     implicit none

     ! Passed variables

     integer :: length
     real(double) :: x(length)
     real(double) :: y(length)
     real(double) :: d2y(length)
     real(double) :: int_x
     real(double) :: int_y

     ! Local variables

     integer :: min, mid, max
     real(double) :: incr, a, b

     if(length==0) then
       int_y = 0.0_double
       return
     end if

     if(int_x <= x(1)) then
       int_y=y(1)
       return
     end if

     if(int_x > x(length)) then
       int_y=0.0_double
       return
     end if

     min=1
     max=length

     do while (max - min > 1)
        mid = (min + max)/2

        if(x(mid) > int_x) then
          max = mid
        else
          min = mid
        end if
     end do

     incr = x(max) - x(min)
     if (incr == 0) then
        write(*,*) 'Error in spline interpolation: Increment = 0'
        stop
     end if

     a = (x(max) - int_x)/incr
     b = (int_x - x(min))/incr
     int_y = a * y(min) + b * y(max) &
           + ((a*a*a-a)*d2y(min) + (b*b*b-b)*d2y(max)) * incr *incr / 6.0_double

!print *,"AST-splint", min, max, x(min), x(max), y(min), y(max)

  end subroutine spline_interpolation
!!***

! -----------------------------------------------------------------------------
! Subroutine simpson_integral_prod
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_numeric/simpson_integral_prod *
!!
!!NAME
!! simpson_integral_prod
!!USAGE
!!
!!PURPOSE
!! Calculate the integral of the product of two functions using Simpson's rule
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 28/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine simpson_integral_prod(f1, f2, length, log_x, result)

     use cq_ut_psp2pao_types

     implicit none

     ! Passed variables

     integer :: length
     real(double) :: f1(length)
     real(double) :: f2(length)
     real(double) :: log_x
     real(double) :: result

     ! Local variables

     integer :: i
     real(double) :: sum

     sum = f1(1) * f2(1) - f1(length) * f2(length)

     do i=2,length-1,2
        sum = sum + (4.0_double* f1(i)*f2(i) + 2.0_double*f1(i+1)*f2(i+1))
     end do
     result = sum * (log_x / 3.0_double);

  end subroutine simpson_integral_prod
!!***

! -----------------------------------------------------------------------------
! Subroutine finite_differentiation
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_numeric/finite_differentiation *
!!
!!NAME
!! finite_differentiation
!!USAGE
!!
!!PURPOSE
!! Calculate the derivative of a function (5-point formula)
!! Note that the internal global scale (gl_r, gl_points_mesh, gl_dnu) is used
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine finite_differentiation(func, dfunc)

     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_global

     implicit none

     ! Passed variables

     real(double) :: func(gl_points_mesh)
     real(double) :: dfunc(gl_points_mesh)

     ! Local variables

     integer :: i

     ! Parameters
     real(double), parameter :: two = 2.0_double
     real(double), parameter :: three = 3.0_double
     real(double), parameter :: four = 4.0_double
     real(double), parameter :: five = 5.0_double
     real(double), parameter :: six = 6.0_double
     real(double), parameter :: twelve = 12.0_double
     real(double), parameter :: twentyfive = 25.0_double

     ! We really need at least five points
     if(gl_points_mesh < 5) then
        write(*,*) 'Error: The density grid must be made of at least five points'
        stop
     end if

     dfunc(1) = -twentyfive * func(1) / twelve &
              + four * func(2) &
              - three * func(3) &
              + four * func(4) / three &
              - func(5) / four

     dfunc(2) = -func(1) / four &
              - five * func(2) / six &
              + three * func(3) / two &
              - func(4) / two &
              + func(5) / twelve

     do i = 3, gl_points_mesh - 2
        dfunc(i) = func(i-2) / twelve &
                 - two * func(i-1) / three &
                 + two * func(i+1) / three &
                 - func(i+2) / twelve
     end do

     dfunc(gl_points_mesh - 1) = -func(gl_points_mesh - 4) / twelve &
                               + func(gl_points_mesh - 3) / two &
                               - three * func(gl_points_mesh - 2) / two &
                               + five * func(gl_points_mesh - 1) / six &
                               + func(gl_points_mesh - 1) / four

     dfunc(gl_points_mesh) = func(gl_points_mesh - 4) / four &
                           - four * func(gl_points_mesh - 3) / three &
                           + three * func(gl_points_mesh - 2) &
                           - four * func(gl_points_mesh - 1) &
                           + twentyfive * func(gl_points_mesh)


     do i=1, gl_points_mesh
        dfunc(i) = dfunc(i) / (gl_r(i) * gl_dnu)
     end do

  end subroutine finite_differentiation
!!***

! -----------------------------------------------------------------------------
! Subroutine integral_cubic
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_numeric/integral_cubic *
!!
!!NAME
!! integral_cubic
!!USAGE
!!
!!PURPOSE
!! Calculate the integral of (part of) a function, using a method by
!! Mike Teter
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine integral_cubic(x, y, dy, ini, end, result)

     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_global

     implicit none

     ! Passed variables
     real(double) :: x(gl_points_mesh)
     real(double) :: y(gl_points_mesh)
     real(double) :: dy(gl_points_mesh)
     integer :: ini, end
     real(double) :: result

     ! Local variables
     integer :: i
     real(double) :: incr

     result = 0.0_double
     do i=ini,end-1
      incr = x(i+1) - x(i)
      result = result &
             + incr * (6.0_double * (y(i) + y(i+1)) &
             + incr * (dy(i) - dy(i+1))) / 12.0_double
     end do

  end subroutine integral_cubic
!!***

end module cq_ut_psp2pao_numeric
