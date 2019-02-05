! $Id$
! -----------------------------------------------------------
! Module spline_module
! -----------------------------------------------------------
! Code area 9: general
! -----------------------------------------------------------

!!****h* Conquest/spline_module
!!  NAME
!! 
!!  PURPOSE
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   23/05/01
!!  MODIFICATION HISTORY
!!   31/05/2001 dave
!!    Added RCS Id tag to header, increased ROBODoc documentation
!!   20/06/2001 dave
!!    Added RCS Log tag and cq_abort
!!   15:48, 08/04/2003 drb 
!!    Bug fixes in splint and dsplint
!!***
module spline_module
  
  use global_module, only: area_general

  implicit none

  ! RCS tag for object file identification 
  character(len=80), save, private :: RCSid = "$Id$"

contains

! -----------------------------------------------------------
! Subroutine spline
! -----------------------------------------------------------

!!****f* spline_module/spline *
!!
!!  NAME 
!!   spline
!!  USAGE
!! 
!!  PURPOSE
!!   Given an array y(1:n) containing a function tabulated at 
!!   REGULAR intervals with spacing delta_x, and given values
!!   yp1 and ypn for the first derivative of the function at 
!!   the end points, this subroutine returns an array d2y(1:n) 
!!   containing the second derivatives of the interpolating 
!!   function at the tabulated points. If yp1 and/or ypn are 
!!   equal or larger than 1.0d30, the routine is signaled to set 
!!   the corresponding boundary condition for a natural spline, 
!!   with zero second derivative at that point. This subroutine
!!   is based on subroutine spline(), from Numerical Recipes (2nd edition). 
!!  INPUTS
!!   integer :: n - size of table to be splined
!!   real(double) :: delta_x - spacing of values in table
!!   real(double) :: yp1, ypn - end points giving boundaries
!!   real(double), dimension(n) :: y, y2 - table and second derivative
!!  USES
!!   datatypes, numbers
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   29/02/96
!!  MODIFICATION HISTORY
!!   23/05/2001 dave
!!    F90, ROBODoc, indentation
!!   31/05/2001 dave
!!    More ROBODoc comments
!!  SOURCE
!!
  subroutine spline( n, delta_x, y, yp1, ypn, y2 )

    use datatypes
    use numbers
    use GenComms, only: cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: n

    real(double) :: delta_x, yp1, ypn
    real(double), dimension(n) :: y, y2

    ! Local variables
    integer :: i, k, stat
    real(double) :: p, qn, sig, un
    real(double), parameter ::huge = 1.0e30_double
    real(double), dimension(:), allocatable :: u, x

    allocate(u(n), x(n), STAT=stat)
    if (stat /= 0) call cq_abort("spline: Error alloc mem: ", n)
    call reg_alloc_mem(area_general, 2*n, type_dbl)

    ! generate the array of x entries
    do i=1, n
       x(i) = delta_x * real(i-1,double)
    end do
    if ( yp1 .gt. huge ) then ! lower boundary condition set to 'natural'
       y2(1) = zero
       u(1) = zero
    else ! or else to have a specified first derivative
       y2(1) = -half
       u(1) = ( three/(x(2) - x(1))) *  &
            ((y(2) - y(1))/(x(2) - x(1)) - yp1 )
    end if

    ! this is the decomposition loop for the tridiagonal algorithm; y2 and
    ! u are used for temporary storage of the decomposed factors
    do i=2, n-1
       sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
       p = sig * y2(i-1) + two
       y2(i) = (sig - one) / p
       u(i) = (six * ((y(i+1) - y(i))/(x(i+1) - x(i)) -  &
            (y(i) - y(i-1))/ &
            (x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig * u(i-1))/p 
    end do

    if (ypn>huge) then ! the upper boundary condition is set to be 'natural'...
       qn = zero
       un = zero
    else ! or else to have a specified first derivative
       qn = half
       un = ( three/(x(n) - x(n-1))) *  &
            (ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
    end if
    y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + one )

    ! this is the backsubstitution loop of the tridiagonal algorithm
    do k=n-1, 1, -1
       y2(k) = y2(k) * y2(k+1) + u(k)
    end do

    deallocate(u, x, STAT=stat)
    if (stat /= 0) call cq_abort("spline: Error dealloc mem")
    call reg_dealloc_mem(area_general, 2*n, type_dbl)

    return
  end subroutine spline
!!***

! -----------------------------------------------------------
! Subroutine spline_nonU
! -----------------------------------------------------------

!!****f* spline_module/spline_nonU *
!!
!!  NAME 
!!   spline_nonU - makes a spline table for non-uniform x array
!!  USAGE
!!   spline_nonU(number, x array, y array,start dy/dx,final dy/dx, d2y/dx2)
!!   spline_nonU(n, y, y, yp1, ypn, y2)
!!  PURPOSE
!!   Given an array y(1:n) containing a function tabulated at 
!!   IRREGULAR intervals given by array x(1:n), and given values
!!   yp1 and ypn for the first derivative of the function at 
!!   the end points, this subroutine returns an array d2y(1:n) 
!!   containing the second derivatives of the interpolating 
!!   function at the tabulated points. If yp1 and/or ypn are 
!!   equal or larger than 1.0d30, the routine is signaled to set 
!!   the corresponding boundary condition for a natural spline, 
!!   with zero second derivative at that point. This subroutine
!!   is based on subroutine spline(), from Numerical Recipes (2nd edition). 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   20/11/97
!!  MODIFICATION HISTORY
!!   23/05/2001 dave
!!    F90, ROBODoc, indentation
!!  SOURCE
!!
  subroutine spline_nonU(n, x, y, yp1, ypn, y2)

    use datatypes
    use numbers
    use GenComms, only: cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables

    integer :: n
    real(double) :: yp1, ypn
    real(double), dimension(n) :: x, y, y2

    !     Local variables
    integer i, k, stat
    real(double) :: p, qn, sig, un
    real(double), parameter :: huge = 1.0e30_double
    real(double), dimension(:), allocatable :: u

    allocate(u(n), STAT=stat)
    if (stat /= 0) call cq_abort("spline_nonU: Error alloc mem: ", n)
    call reg_alloc_mem(area_general, n, type_dbl)

    ! set up boundary conditions
    if (yp1>huge) then ! lower boundary condition set to 'natural'
       y2(1) = zero
       u(1) = zero
    else ! or else to have a specified first derivative
       y2(1) = -half
       u(1) = ( three/(x(2) - x(1))) * ((y(2) - y(1))/(x(2) - x(1)) - yp1 )
    end if

    ! this is the decomposition loop for the tridiagonal algorithm; y2 and
    ! u are used for temporary storage of the decomposed factors
    do i=2, n-1
       sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
       p = sig * y2(i-1) + two
       y2(i) = (sig - one) / p
       u(i) = (six * ((y(i+1) - y(i))/(x(i+1) - x(i)) -  &
            (y(i) - y(i-1))/ &
            (x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig * u(i-1))/p 
    end do

    if (ypn>huge) then ! the upper boundary condition is set to be 'natural'...
       qn = zero
       un = zero
    else ! or else to have a specified first derivative
       qn = half
       un = ( three/(x(n) - x(n-1))) *  &
            (ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
    end if
    y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + one )

    ! this is the backsubstitution loop of the tridiagonal algorithm
    do k=n-1, 1, -1
       y2(k) = y2(k) * y2(k+1) + u(k)
    end do

    deallocate(u, STAT=stat)
    if (stat /= 0) call cq_abort("spline_nonU: Error dealloc mem")
    call reg_dealloc_mem(area_general, n, type_dbl)

    return
  end subroutine spline_nonU
!!***

! -----------------------------------------------------------
! Subroutine splint_nonU
! -----------------------------------------------------------

!!****f* spline_module/splint_nonU *
!!
!!  NAME 
!!   splint_nonU
!!  USAGE
!!   splint_nonU(xarray, yarray, d2y/dx2 array, size of array, x value, y, flag)
!!   splint_nonU(xa,ya,y2a,n,x,y, flag)
!!  PURPOSE
!!   NR routine for interpolation using spline table set up by spline
!!   Takes a spline table with non-uniform x spacings and, for the 
!!   specified value x, finds the y value.
!!  INPUTS
!!   integer :: n (size of table)
!!   real(double) :: x,y (x value specified, y value to be found)
!!   real(double), dimension(n) :: xa, ya, y2a (x, y, d2y/dx2 arrays)
!!   logical :: flag (flag to warn if given x is out of tabulated range)
!!  USES
!!   datatypes
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   20/11/97
!!  MODIFICATION HISTORY
!!   23/05/2001 dave
!!    F90, ROBODoc, indentation
!!   20/06/2001 dave
!!    Added cq_abort
!!   15/04/02 lkd
!!    Changed title from splint to splint_nonU for consistency with rest of module
!!   15/04/02 lkd
!!    Added "use numbers" statement
!!   17/04/2002 lkd
!!    Added logical flag to return false if all OK, true if the given x is out of range
!!   25/01/2019 tsuyoshi
!!    added RD_ERR term for the round-error problem
!! 
!!  SOURCE
!!
  subroutine splint_nonU(xa,ya,y2a,n,x,y,flag)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort

    ! Passed variables
    integer, intent(in) :: n
    real(double), intent(in) :: x
    real(double), intent(out) :: y
    real(double), dimension(n), intent(in) :: xa,y2a,ya
    logical, intent(out) :: flag

    ! Local variables
    integer :: k,khi,klo
    real(double) :: a,b,h

    ! initialize warning flag
    flag = .false.
    !ORI if (x .gt. xa(n) .or. x .lt. xa(1)) then
    if (x .gt. xa(n)+RD_ERR .or. x .lt. xa(1)-RD_ERR) then
       flag = .true. ! the given x lies outside the tabulated values, warn the calling routine
    end if

    ! Search through the table and bracket x value
    klo = 1
    khi = n
    do while((khi-klo)>1)
       k=(khi+klo)/2
       if(xa(k).gt.x) then ! Increment high point
          khi=k
       else                ! Increment low point
          klo=k
       end if
    end do
    h=xa(khi)-xa(klo)      ! Find x spacing
   !2019/Jan/29 TM
    if(khi == klo) then
      if (khi == n .or. khi == 1) then
         y = ya(khi)
         return
      else
         call cq_abort('splint_nonU: bisection failed, h=0')
      endif
    endif
    !ori if(h==0) then
    !ori    call cq_abort('splint_nonU: bisection failed, h=0')
    !ori end if
   !2019/Jan/29 TM
    ! Create value of y
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+ &
         ( (a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi) )*(h*h)/six
    return
  end subroutine splint_nonU
!!***

! -----------------------------------------------------------
! Subroutine splint
! -----------------------------------------------------------

!!****f* spline_module/splint *
!!
!!  NAME 
!!   splint
!!  USAGE
!!   splint(xarray, yarray, d2y/dx2 array, size of array, x value, y)
!!   splint(xa,ya,y2a,n,x,y)
!!  PURPOSE
!!   NR routine for interpolation using spline table set up by spline
!!   Takes a spline table with uniform x spacing delta_x and, for the 
!!   specified value x, finds the y value.
!!  INPUTS
!!   integer :: n (size of table)
!!   real(double) :: x, y, delta_x (x value specified, y value to be found, x-value spacing)
!!   real(double), dimension(n) :: ya, y2a (y, d2y/dx2 arrays)
!!   logical :: flag (flag to warn if given x is out of tabulated range)
!!  USES
!!   datatypes, numbers, GenComms
!!  AUTHOR
!!   L.K.Dash
!!  CREATION DATE
!!   15/04/02
!!  MODIFICATION HISTORY
!!   15:47, 08/04/2003 drb 
!!    Fixed major bug: array references were all one element off
!!   10:00, 28/08/2003 drb 
!!    Fixed ANOTHER major bug: I'd forgotten that the array reference was used to create the x position so needs one taken off it !
!!   08:56, 2003/11/20 dave
!!    Added return to prevent array overrun
!!   25/01/2019 tsuyoshi
!!    added RD_ERR term for the round-error problem
!!  SOURCE
!!
  subroutine splint(delta_x,ya,y2a,n,x,y,flag)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort

    ! Passed variables
    integer, intent(in) :: n
    real(double), intent(in) :: x
    real(double), intent(out) :: y
    real(double), intent(in) :: delta_x
    real(double), dimension(n), intent(in) :: y2a,ya
    logical, intent(out) :: flag

    ! Local variables
    integer :: khi,klo
    real(double) :: a,b,xlo,xhi

    ! initialize warning flag
    flag = .false.
    !ORI if (x>real(n-1,double)*delta_x .or. x<0) then
    if (x>real(n-1,double)*delta_x+RD_ERR .or. x<-RD_ERR) then
       flag = .true. ! the given x lies outside the tabulated values, warn the calling routine
       y = 0.0_double ! Avoid array overruns
       return
    end if

    ! bracket x value   
    klo = aint(x/delta_x) + 1 ! integer division gives lower bound, add 1 for array value
    khi = klo + 1
    !2019/Jan/29 TM
    !ori if(khi>n) call cq_abort("Error in splint: index beyond range ",khi,n)
     if(khi>n+1) then
        call cq_abort("Error in splint: index beyond range ",khi,n)
     elseif(khi==n+1) then
        y = 0.0_double
        return
     endif
     
    
    ! Remember that we added to get array value and subtract it off again !
    xlo =delta_x*real(klo-1,double)
    xhi =delta_x*real(khi-1,double)
    
    ! Create value of y
    a = (xhi-x)/delta_x
    b = (x-xlo)/delta_x
    y=a*ya(klo)+b*ya(khi)+ &
         ( (a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi) )*(delta_x*delta_x)/six
    return
  end subroutine splint
!!***

! -----------------------------------------------------------
! Subroutine dsplint
! -----------------------------------------------------------

!!****f* spline_module/dsplint *
!!
!!  NAME 
!!   dsplint
!!  USAGE
!!   dsplint(xarray, yarray, d2y/dx2 array, size of array, x value, dy)
!!   dsplint(xa,ya,y2a,n,x,dy)
!!  PURPOSE
!!   NR routine for interpolation using spline table set up by spline
!!   Takes a spline table with uniform x spacing delta_x and, for the 
!!   specified value x, finds the derivative value.
!!  INPUTS
!!   integer :: n (size of table)
!!   real(double) :: x, y, delta_x (x value specified, y value to be found, x-value spacing)
!!   real(double), dimension(n) :: ya, y2a (y, d2y/dx2 arrays)
!!   logical :: flag (flag to warn if given x is out of tabulated range)
!!  USES
!!   datatypes, numbers, GenComms
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07:50, 2003/03/19 dave
!!  MODIFICATION HISTORY
!!   15:47, 08/04/2003 drb 
!!    Fixed major bug: array references were all one element off (carried over from splint)
!!   11:04, 03/09/2003 drb 
!!    Added xlo/xhi correction that was added to splint (should have been done earlier !)
!!   08:56, 2003/11/20 dave
!!    Added return to prevent array overrun
!!   08:01, 2003/12/04 dave
!!    Added return of spline-interpolated value as well as derivative
!!   25/01/2019 tsuyoshi
!!    added RD_ERR term for the round-error problem
!!  SOURCE
!!
  subroutine dsplint(delta_x,ya,y2a,n,x,y,dy,flag)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort

    ! Passed variables
    integer, intent(in) :: n
    real(double), intent(in) :: x
    real(double), intent(out) :: y,dy
    real(double), intent(in) :: delta_x
    real(double), dimension(n), intent(in) :: y2a,ya
    logical, intent(out) :: flag

    ! Local variables
    integer :: khi,klo
    real(double) :: a,b,da,db,dc,dd, xlo, xhi

    ! initialize warning flag
    flag = .false.
    !ORI if (x>(n-1)*delta_x.OR.x<0) then
    if (x>(n-1)*delta_x+RD_ERR.OR.x<-RD_ERR) then
       flag = .true. ! the given x lies outside the tabulated values, warn the calling routine
       y = 0.0_double
       dy = 0.0_double ! Avoid array overruns
       return
    end if

    ! bracket x value   
    klo = aint(x/delta_x) + 1 ! integer division gives lower bound, add 1 for array value
    khi = klo + 1
    ! Remember that we added to get array value and subtract it off again !
    xlo =delta_x*real(klo-1,double)
    xhi =delta_x*real(khi-1,double)

    !2019/Jan/29 TM
     if(klo == n) then
       y = 0.0_double
       dy = 0.0_double ! Avoid array overruns
       return
     endif
    !2019/Jan/29 TM

    ! Create value of y
    a = (xhi-x)/delta_x
    b = (x-xlo)/delta_x
    da = -one/delta_x
    db =  one/delta_x
    dc = -delta_x*(three*a*a-one)/six
    dd =  delta_x*(three*b*b-one)/six
    y=a*ya(klo)+b*ya(khi)+ &
         ( (a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi) )*(delta_x*delta_x)/six
    dy =  da*ya(klo)+db*ya(khi)+dc*y2a(klo)+dd*y2a(khi)
    return
  end subroutine dsplint
!!***

end module spline_module



