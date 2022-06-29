! -*- mode: F90; mode: font-lock -*-
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
!!   Collects mathematical and other useful functions in one module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2016/02/09
!!  MODIFICATION HISTORY
!!   2019/11/12 08:13 dave
!!    Adding sort routines
!!  SOURCE
!!
module functions

  implicit none

  !!***

contains


  ! -----------------------------------------------------------
  ! Function is_prime
  ! -----------------------------------------------------------

  !!****f* functions/is_prime *
  !!
  !!  NAME
  !!   prime
  !!  USAGE
  !!   prime(n)
  !!  PURPOSE
  !!   Return if n is a prime number
  !!  INPUTS
  !!   integer :: n
  !!  AUTHOR
  !!   LAT
  !!  CREATION DATE
  !!   2021/06/14
  !!  SOURCE
  !!

  logical function is_prime(n)

    use datatypes
    
    implicit none

    integer, intent(in) :: n
    !
    integer :: i, inf, sup
    !
    is_prime = .true.
    !
    inf = 2
    sup = int( sqrt( real(n,double) ) + 1 )
    !
    do i = inf, sup
       if ( mod(n,i) == 0 ) then
          is_prime = .false.
       end if
    end do
    !
    return
    !
  end function is_prime
    
  ! -----------------------------------------------------------
  ! Function erfc
  ! -----------------------------------------------------------

  !!****f* functions/erfc *
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
  real(double) function erfc_cq(x)

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
    real(double) :: an, b, c, d
    integer :: i

    if(x < zero) then
       y = -x
    else
       y = x
    end if
    ! This expects y^2
    y2 = y*y
    if(y<RD_ERR) then
       erfc_cq = one
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
       erfc_cq = one - sum * exp(-y2 + half * log(y2) - erfc_gln)
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
       erfc_cq = sum * exp(-y2 + half * log(y2) - erfc_gln)
    end if
    if (x < zero) erfc_cq = two - erfc_cq
    return
  end function erfc_cq
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

  !!****f* functions/heapsort_integer_index *
  !!
  !!  NAME 
  !!   heapsort_integer_index
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Heap sort for an integer array returning sorted index (not in-place)
  !!   Note that this may appear subtly different to standard implementations
  !!   Arrays starting at 1 have different parent/child formulae to arrays
  !!   starting at 0.  Contains the sift_down routine.
  !!  INPUTS
  !!   
  !!   
  !!  USES
  !!   
  !!  AUTHOR
  !!   David Bowler
  !!  CREATION DATE
  !!   2019/11/12
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!  
  subroutine heapsort_integer_index(n_arr,array,arr_index,reverse)

    ! Passed variables
    integer :: n_arr
    integer, OPTIONAL :: reverse
    integer, dimension(n_arr) :: array
    integer, dimension(n_arr) :: arr_index, tmp_index

    ! Local variables
    integer :: i, temp

    ! Set up index
    do i=1,n_arr
       arr_index(i) = i
    end do
    if(n_arr==1) return
    ! Create heap
    do i = n_arr/2, 1, -1
       call sift_down_integer_index(i, n_arr)
    end do
    ! Sort array
    do i=n_arr, 2, -1
       temp = arr_index(1)
       arr_index(1) = arr_index(i)
       arr_index(i) = temp
       call sift_down_integer_index(1, i-1)
    end do
    if(present(reverse)) then
       if(reverse==1) then
          do i=1,n_arr
             tmp_index(n_arr-i+1) = arr_index(i)
          end do
          arr_index = tmp_index
       end if
    end if
    return
    
  contains

    ! Note that we inherit access to array and arr_index while this is contained
    ! in the original routine
    subroutine sift_down_integer_index(start,end)

      ! Passed variables
      integer :: start, end

      ! Local variables
      integer :: child, root, temp

      root = start
      ! Left child is i*2
      child = root*2
      do while(child<=end)
         ! Right child is i*2 + 1
         if(child+1<=end) then
            if(array(arr_index(child))<array(arr_index(child+1))) child = child+1
         end if
         if(array(arr_index(root))<array(arr_index(child))) then
            temp = arr_index(child)
            arr_index(child) = arr_index(root)
            arr_index(root) = temp
            root = child
            child = root*2
         else
            return
         end if
      end do
      return
    end subroutine sift_down_integer_index
    
  end subroutine heapsort_integer_index

  !!****f* functions/heapsort_real_index *
  !!
  !!  NAME 
  !!   heapsort_real_index
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Heap sort for an integer array returning sorted index (not in-place)
  !!   Note that this may appear subtly different to standard implementations
  !!   Arrays starting at 1 have different parent/child formulae to arrays
  !!   starting at 0.  Contains the sift_down routine.
  !!
  !!   Contains optional flag to signal descending rather than ascending order
  !!  INPUTS
  !!   
  !!   
  !!  USES
  !!   
  !!  AUTHOR
  !!   David Bowler
  !!  CREATION DATE
  !!   2019/11/12
  !!  MODIFICATION HISTORY
  !!   2020/01/24 08:05 dave
  !!    Corrected type of temp to integer
  !!  SOURCE
  !!  
  subroutine heapsort_real_index(n_arr,array,arr_index,reverse)

    use datatypes

    implicit none
    
    ! Passed variables
    integer :: n_arr
    integer, OPTIONAL :: reverse
    real(double), dimension(n_arr) :: array
    integer, dimension(n_arr) :: arr_index, tmp_index

    ! Local variables
    integer :: i, temp

    ! Set up index
    do i=1,n_arr
       arr_index(i) = i
    end do
    if(n_arr==1) return
    ! Create heap
    do i = n_arr/2, 1, -1
       call sift_down_real_index(i, n_arr)
    end do
    ! Sort array
    do i=n_arr, 2, -1
       temp = arr_index(1)
       arr_index(1) = arr_index(i)
       arr_index(i) = temp
       call sift_down_real_index(1, i-1)
    end do
    ! Reverse order (descending) if required
    if(present(reverse)) then
       if(reverse==1) then
          do i=1,n_arr
             tmp_index(n_arr-i+1) = arr_index(i)
          end do
          arr_index = tmp_index
       end if
    end if
    return
    
  contains

    ! Note that we inherit access to array and arr_index while this is contained
    ! in the original routine
    subroutine sift_down_real_index(start,end)

      use datatypes

      implicit none
    
      ! Passed variables
      integer :: start, end

      ! Local variables
      integer :: child, root, temp

      root = start
      ! Left child is i*2
      child = root*2
      do while(child<=end)
         ! Right child is i*2 + 1
         if(child+1<=end) then
            if(array(arr_index(child))<array(arr_index(child+1))) child = child+1
         end if
         if(array(arr_index(root))<array(arr_index(child))) then
            temp = arr_index(child)
            arr_index(child) = arr_index(root)
            arr_index(root) = temp
            root = child
            child = root*2
         else
            return
         end if
      end do
      return
    end subroutine sift_down_real_index
    
  end subroutine heapsort_real_index
  
end module functions
