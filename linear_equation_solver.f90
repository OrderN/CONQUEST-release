! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module linear_equation_solver
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/linear_equation_solver *
!!  NAME
!!   linear_equation_solver
!!  PURPOSE
!!   Sole purpose is to contain the three subroutines linsolv, 
!!   ludcmp and lubksb 
!!  USES
!!
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
module linear_equation_solver
  implicit none
  save

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"

!!***

contains

! -----------------------------------------------------------
! Subroutine linsolv
! -----------------------------------------------------------

!!****f* linear_equation_solver/linsolv *
!!
!!  NAME 
!!   sbrt_name
!!  USAGE
!! 
!!  PURPOSE
!!   solves the system of linear eqns A . x = b using 
!!   LU decomposition. The code is based on Numerical Recipes, 
!!   2nd ed. Sec. 2.3
!!   sbrt linsolv has as sole purpose the calling of
!!   of the two routines ludcmp and lubksb, which
!!   do the actual work.
!! 
!!  INPUTS
!!   a: the elements of the matrix A
!!   n: the number of linear equations
!!   b: the vector b on the right-hand side
!!   On termination, the solution vector x is in
!!   the array b. 
!! 
!!  USES
!! 
!!  AUTHOR
!!
!!  CREATION DATE
!!
!!
!!  MODIFICATION HISTORY
!!
!!
!!  SOURCE
  subroutine linsolv(a,n,b)
    use datatypes
    implicit none

    ! input-output variables
    integer, intent(in) :: n
    real(double), intent(inout), dimension(:) :: b
    real(double), intent(inout), dimension(:,:) :: a
    ! internal variables
    integer, dimension(n) :: indx
    real(double) :: d

    call ludcmp(a,n,indx,d)
    call lubksb(a,n,indx,b)

    return

  end subroutine linsolv

  ! ----------------------------------------------------------------
!!*** 

! -----------------------------------------------------------
! Subroutine ludcmp
! -----------------------------------------------------------

!!****f* linear_equation_solver/ludcmp *
!!
!!  NAME 
!!   ludcmp
!!  USAGE
!! 
!!  PURPOSE
!!   See Numerical Recipes
!! 
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   Numerical Recipes, transcribed with minor changes by Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!
!!
!!  SOURCE
  subroutine ludcmp(a,n,indx,d)
    use datatypes
    use GenComms, ONLY: cq_abort
    use numbers

    implicit none

    ! input-output variables
    integer, intent(in) :: n
    integer, intent(out), dimension(:) :: indx
    real(double), intent(out) :: d
    real(double), intent(inout), dimension(:,:) :: a
    ! internal variables
    integer :: i, j, k, imax
    real(double) :: aamax, dum, sum
    real(double), dimension(n) :: vv

    d = one
    do i = 1, n
       aamax = zero
       do j = 1, n
          if(abs(a(i,j)) > aamax) aamax = abs(a(i,j))
       end do
       if(aamax == 0.) then
	  call cq_abort('ludcmp: matrix is singular')
       endif
       vv(i) = one/aamax
    end do

    do j = 1, n
       if(j > 1) then
          do i = 1, j-1
             sum = a(i,j)
             if(i > 1) then
                do k = 1, i-1
                   sum = sum - a(i,k)*a(k,j)
                enddo
             end if
             a(i,j) = sum
          end do
       end if

       aamax = zero
       do i = j, n
          sum = a(i,j)
          if(j > 1) then
             do k = 1, j-1
                sum = sum - a(i,k)*a(k,j)
             enddo
          endif
          a(i,j) = sum
          dum = vv(i)*abs(sum)
          if(dum >= aamax) then
             imax = i
             aamax = dum
          end if
       end do

       if(j /= imax) then
          do k = 1, n
             dum = a(imax,k)
             a(imax,k) = a(j,k)
             a(j,k) = dum
          enddo
          d = -d
          vv(imax) = vv(j)
       end if

       indx(j) = imax
       if(a(j,j) == zero) then
	  call cq_abort('ludcmp: zero pivot')
       endif

       if(j /= n) then
          dum = one/a(j,j)
          if(j < n) then
             do i = j+1, n
                a(i,j) = a(i,j)*dum
             end do
          end if
       end if

    end do

    return

  end subroutine ludcmp

  ! -------------------------------------------------------------------
!!*** 

! -----------------------------------------------------------
! Subroutine lubksb
! -----------------------------------------------------------

!!****f* linear_equation_solver/lubksb *
!!
!!  NAME 
!!   lubksb
!!  USAGE
!! 
!!  PURPOSE
!!   See Numerical Recipes
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!!
!!  SOURCE
  subroutine lubksb(a,n,indx,b)
    use datatypes
    use numbers
    implicit none

    ! input-output variables
    integer, intent(in) :: n
    integer, intent(in), dimension(:) :: indx
    real(double), intent(inout), dimension(:) :: b
    real(double), intent(in), dimension(:,:) :: a
    ! internal variables
    integer :: i, ii, j, kk
    real(double) :: sum

    ii = 0
    do i = 1, n
       kk = indx(i)
       sum = b(kk)
       b(kk) = b(i)
       if(ii /= 0) then
          do j = ii, i-1
             sum = sum - a(i,j)*b(j)
          end do
       else if(sum /= zero) then
          ii = i
       end if
       b(i) = sum
    end do

    do i = n, 1, -1
       sum = b(i)
       do j = i+1, n
          sum = sum - a(i,j)*b(j)
       end do
       b(i) = sum/a(i,i)
    enddo

    return

  end subroutine lubksb
!!***

end module linear_equation_solver
