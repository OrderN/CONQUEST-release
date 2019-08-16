! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module pao_array_utility
! ------------------------------------------------------------------------------
! Code area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/pao_array_utility
!!  NAME
!!   pao_array_utility
!!  PURPOSE
!!   contains routines to reproject an array of given grid spacing onto
!!   a different grid spacing
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2008/02/06 08:30 dave
!!    Changed for output to file not stdout
!!   2019/08/16 12:12 dave
!!    Now uses splines, also changed name
!!  SOURCE
!!
module pao_array_utility
   
  use global_module, ONLY: io_lun

  implicit none
  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"

!!***
contains

!!****f* pao_array_utility/reproject_array *
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
    use numbers
    use splines, ONLY: spline_nonU

    implicit none

    integer, intent(in) :: n_old,n_new
    real(double), intent(in) :: del_old,del_new
    real(double), intent(in), dimension(n_old) :: yin
    real(double), intent(out), dimension(n_new) :: yout

    integer :: i,j,l,m
    real(double), allocatable, dimension(:) :: y2,xin
    real(double) :: yp1,ypn, x,y, xx
    real(double) :: a, b, c, d, r1, r2, r3, r4

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

    call spline_nonU(n_old,xin,yin,yp1,ypn,y2)

    do i=1,n_new
       x=real(i-1,double)*del_new
       j = floor(x/del_old) + 1
       if(j+1<=n_old) then
          xx = real(j,double)*del_old
          a = (xx - x)/del_old
          b = one - a
          c = a * ( a * a - one ) * del_old * del_old / six
          d = b * ( b * b - one ) * del_old * del_old / six
          r1 = yin(j)
          r2 = yin(j+1)
          r3 = y2(j)
          r4 = y2(j+1)
          yout(i) = a*r1 + b*r2 + c*r3 + d*r4
       end if
    enddo
    deallocate(xin)
    deallocate(y2)

  end subroutine reproject_array
!!***

!!****f* pao_array_utility/matcharrays *
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
 
end module pao_array_utility

       
       
  
