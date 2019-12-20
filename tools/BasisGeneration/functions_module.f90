! -*- mode: F90; mode: font-lock; column-number-mode: true -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module functions
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/functions *
!!  NAME
!!   functions
!!  PURPOSE
!!   Collects mathematical functions in one module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2016/02/09
!!  MODIFICATION HISTORY
!!  SOURCE
!!
module functions

  implicit none

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"
  !!***

contains

  ! -----------------------------------------------------------
  ! Function erfc
  ! -----------------------------------------------------------

  !!****f* DiagModle/erfc *
  !!
  !!  NAME
  !!   erfc
  !!  USAGE
  !!   erfc(x)
  !!  PURPOSE
  !!   Calculated the complementary error function to rounding-error
  !!   based on erfc() in ewald_module
  !!   accuracy
  !!  INPUTS
  !!   real(double) :: x, argument of complementary error function
  !!  AUTHOR
  !!   Iain Chapman/Lianeng Tong
  !!  CREATION DATE
  !!   2001 sometime/2010/07/26
  !!  MODIFICATION HISTORY
  !!   2016/02/09 08:14 dave
  !!    Moved to functions module from DiagModule
  !!  SOURCE
  !!
  real(double) function erfc(x)

    use datatypes
    use numbers,  only: RD_ERR, one, zero, half, two
    use GenComms, only: cq_abort

    implicit none
    
    real(double), parameter :: erfc_delta = 1.0e-12_double, &
         erfc_gln = 0.5723649429247447e0_double, &
         erfc_fpmax = 1.e30_double
    integer, parameter:: erfc_iterations = 10000

    real(double), intent(in) :: x

    ! local variables
    real(double) :: y, y2
    real(double) :: ap, sum, del
    real(double) :: an, b, c, d, h
    integer :: i

    if(x < zero) then
       y = -x
    else
       y = x
    end if
    ! This expects y^2
    y2 = y*y
    if(y<RD_ERR) then
       erfc = one
       return
    end if
    if (y2 < 2.25_double) then
       ap = half
       sum = two
       del = sum
       do i = 1, erfc_iterations
          ap = ap + 1.0_double
          del = del * y2 / ap
          sum = sum + del
          if (abs(del) < abs(sum) * erfc_delta) exit
       end do
       erfc = one - sum * exp(-y2 + half * log(y2) - erfc_gln)
    else
       b = y2 + half
       c = erfc_fpmax
       d = one / b
       sum = d
       do i = 1, erfc_iterations
          an = - i * (i - half)
          b = b + two
          d = an * d + b
          c = b + an / c
          d = one / d
          del = d * c
          sum = sum * del
          if (abs(del - one) < erfc_delta) exit
       end do
       erfc = sum * exp(-y2 + half * log(y2) - erfc_gln)
    end if
    if (x < zero) erfc = two - erfc
    return
  end function erfc
  !!***

!!****f* functions/j0 *
!!
!!  NAME
!!   j0
!!  USAGE
!!
!!  PURPOSE
!!   Calculates 0th-order Bessel function
!! INPUTS
!!   x
!! OUTPUTS
!!
!!  USES
!!
!!  AUTHOR
!!   N. Watanabe (Mizuho) with TM, DRB
!!  CREATION DATE
!!   2014
!!  MODIFICATION HISTORY
!!   2015/11/09 17:28 dave
!!    - Moved into pseudo_tm_module
!!   2016/02/09 08:23 dave
!!    Moved into functions module
!!  SOURCE
!!
  function j0( x )

    use datatypes
    use numbers, only: very_small, one_sixth, one

    implicit none

    real(double) :: x
    real(double) :: j0

    if( x<very_small ) then
       j0 = one - one_sixth*x*x
    else
       j0 = sin(x)/x
    endif
  end function j0
!!***

!!****f* functions/j1 *
!!
!!  NAME 
!!   j1
!!  USAGE
!!   
!!  PURPOSE
!!   1st order bessel function with provision for very small numbers
!!  INPUTS
!!   
!!   
!!  USES
!!   
!!  AUTHOR
!!   NW (Mizuho) with TM and DRB
!!  CREATION DATE
!!   Mid 2014
!!  MODIFICATION HISTORY
!!   2016/02/09 08:24 dave
!!    Moved to functions module
!!  SOURCE
!!  
  function j1x( x )

    use datatypes
    use numbers

    implicit none

    real(double) :: x
    real(double) :: j1x

    if( x<very_small ) then
       j1x = one_third - one/30.0_double*x*x
    else
       j1x = (sin(x)-x*cos(x))/(x*x*x)
    endif
  end function j1x
!!***

  
end module functions
