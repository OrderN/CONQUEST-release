! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module vdWMesh_module
! ------------------------------------------------------------------------------
! Code area 12: Integration
! ------------------------------------------------------------------------------

!!****h* Conquest/vdWMesh_module
!! NAME
!!   vdWMesh_module
!! PURPOSE
!!   Contains utility routines for meshes used by Van der
!!   Waals functional
!! METHODS
!!   derivative         : Returns function derivatives the same mesh points
!!   get_mesh           : Returns a preciously set 1D mesh
!!   get_n              : Returns the number of mesh points
!!   integral           : Returns the integral of a function defined in a mesh
!!   interpolation      : Returns interpolated values at arbitary points
!!   locate             : Given x0, it returns real value i0 such that x(i0) = x0
!!   numerov            : Solves d2y/dx2 = f(x,y) = f0(x) + f1(x)*y
!!   polint             : Lagrange polynomial interoplation
!!   set_interpolation  : Sets interpolation method (lagrange | spline)
!!   set_mesh           : Sets a uniform, logarithmic, or arbitary 1D mesh
!! AUTHOR
!!   L.Tong
!! CREATION DATE
!!   2012/02/26
!! MODIFICATION HISTORY
!! SOURCE
!!
module vdWMesh_module

  use datatypes
  use global_module, only: area_integn

  implicit none

  ! public procedures
  public :: &
       derivative,        &! Returns function derivatives the same mesh points
       get_mesh,          &! Returns a preciously set 1D mesh
       get_n,             &! Returns the number of mesh points
       integral,          &! Returns the integral of a function defined in a mesh
       interpolation,     &! Returns interpolated values at arbitary points
       locate,            &! Given x0, it returns real value i0 such that x(i0) = x0
       numerov,           &! Solves d2y/dx2 = f(x,y) = f0(x) + f1(x)*y
       polint,            &! Lagrange polynomial interoplation
       set_interpolation, &! Sets interpolation method (lagrange | spline)
       set_mesh            ! Sets a uniform, logarithmic, or arbitary 1D mesh

  private ! nothings is declared public beyond this point

  character(len=80), save :: &
       RCSid = "$Id$"

  ! Internal parameters
  ! Max. exponent coef. in log. mesh
  real(double), parameter :: amax = 1.e-6_double
  ! Tol. for mesh type identification
  real(double), parameter :: xtol = 1.e-12_double
  ! Tolerance for interpolation of x(i)
  real(double), parameter :: itol = 1.e-12_double

  ! Saved internal variables and arrays
  character(len=11), save:: mesh_type = 'unknown'
  character(len=10), save:: interpolation_method = 'spline'
  logical,  save:: defined_mesh = .false.  ! Has a mesh been defined?
  real(double), save:: aa, b ! Log-mesh coefs: x(i)=x(1)+b*(exp(aa*(i-1))-1)
  real(double), save:: d     ! Uniform mesh interval: x(i)=x(1)+d*(i-1)
  real(double), save:: yp1 = huge(yp1) ! dy/dx at x(1) for spline interp.
  real(double), save:: ypn = huge(ypn) ! dy/dx at x(n) for spline interp.
  real(double), save, dimension(:), allocatable :: &
       ivec,  &! Auxiliary vector to call polint: ivec(i)=i
       sqrxp, &! Sqrt(dx/di) at mesh points x(i). Used by numerov.
       s0,    &! (dx/di)**(3/2) at mesh points. Used by numerov.
       s1,    &! (dx/di)**2 at mesh points. Used by numerov.
       s2,    &! (3*xp2**2 - 2*xp1*xp3) / (4*xp1**2). Used by numerov.
       xi,    &! Mesh points x(i)
       xp1,   &! dx/di at mesh points
       xp2,   &! d2x/di2 at mesh points
       xp3,   &! d3x/di3 at mesh points
       xp4     ! d4x/di4 at mesh points
!!*****

contains

  !!****f* vdWMesh_module/get_n
  !! PURPOSE
  !!   Returns the number of mesh points
  !! USAGE
  !!   n = get_n( xmin, xmax, dxmin, dxmax )
  !! INPUTS
  !!   real(double) xmin  : First mesh point
  !!   real(double) xmax  : Last mesh point
  !!   real(double) dxmin : First mesh interval: x(2)-x(1)
  !!   real(double) dxmax : Next to last mesh interval: x(n+1)-x(n)
  !! RETURN VALUE
  !!   integer n : Number of logarithmic-mesh points
  !! NOTES
  !!   - All arguments are input and required
  !!   - No mesh needs to have been previously set. Instead, this
  !!     function is to help to find the argument n needed to set it.
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/27
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function get_n(xmin, xmax, dxmin, dxmax)

    use datatypes
    use numbers
    use GenComms, only: cq_abort

    implicit none

    ! input parameters
    real(double), intent(in) :: xmin
    real(double), intent(in) :: xmax
    real(double), intent(in) :: dxmin
    real(double), intent(in) :: dxmax
    ! result
    integer :: get_n

    real(double) :: a, b

    ! Check signs of arguments
    if (dxmax * dxmin <= zero .or. (xmax - xmin) * dxmax <= zero) &
         call cq_abort('get_n: ERROR: Bad arguments')

    ! Find required number of points
    if (dxmax == dxmin) then  ! Linear mesh
       get_n = nint((xmax-xmin) / dxmax)
    else
       b = (xmax - xmin) * dxmin / (dxmax - dxmin)
       a = log(1 + dxmin / b)
       get_n = 1 + nint(log(1 + (xmax - xmin) / b) / a)
    end if

  end function get_n
  !!*****


  !!****f* vdWMesh_module/set_mesh
  !! PURPOSE
  !!   Returns a preciously set 1D mesh
  !! USAGE
  !!   call set_mesh( n, x, xmin, xmax, a, dxndx1 )
  !! INPUTS
  !!   integer  n          : Number of mesh points
  !!   real(double) x(n)   : Mesh point values
  !!   real(double) xmin   : Value of first mesh point
  !!   real(double) xmax   : Value of last mesh point
  !!   real(double) a      : Fractional increment of succesive mesh
  !!                         intervals of a logarithmic mesh:
  !!                         x(i) = b * (exp(a*(i-1)) - 1)
  !!   real(double) dxndx1 : Ratio between last and first intervals
  !!                         in a logarithmic mesh
  !! OUTPUT
  !!   Allocates (reallocate if nessary) and sets the module mesh
  !!   variables
  !!   ivec(n), sqrxp(n), s0(n), s1(n), s2(n), xi(n), xp1(n),
  !!   xp2(n), xp3(n), xp4(n) and
  !!   aa, b if log mesh
  !!   d     if uniform mesh
  !! NOTES
  !!   - All arguments are input
  !!   - All the arguments, except n, are optional.
  !!   - If x is present, no other optional arguments may be present
  !!   - If x is not present, xmax must be present, and only one of
  !!     a or dxndx1 may be present
  !!   - If xmin is not present, xmin=0 is assumed
  !!   - If only xmax (and optionally xmin) is present, a regular mesh
  !!     is generated
  !!   - Stops with an error message if x is present but x(i) is not monotonic
  !!   - If x is present and it is a uniform or a logarithmic mesh, it will be
  !!     identified as such, and analytic formulas will be used accordingly.
  !!     However, this takes extra time, so it is preferable to use the
  !!     other arguments to set a uniform or a logarithmic mesh.
  !! EXAMPLE
  !!   - To initialize a previously generated mesh:
  !!     call set_mesh(n, x=x_mesh(1:n))
  !!   - To generate a regular mesh:
  !!     call set_mesh(n, xmax=x_max)
  !!   - To generate a logarithmic mesh:
  !!     call set_mesh(n, dxndx1=delta_x_max/delta_x_min)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/27
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine set_mesh(n, x, xmin, xmax, a, dxndx1)

    use datatypes
    use numbers
    use GenComms, only: cq_abort

    implicit none
    ! passed parameters
    integer,                intent(in) :: n
    real(double), optional, intent(in) :: x(n)
    real(double), optional, intent(in) :: xmin
    real(double), optional, intent(in) :: xmax
    real(double), optional, intent(in) :: a
    real(double), optional, intent(in) :: dxndx1
    ! local variables
    integer      :: i
    real(double) :: e, dx, dd, x1

    ! Check consistency of optional arguments
    if (present(x) .and. &
        (present(xmin) .or. present(xmax) .or. present(a) .or. present(dxndx1))) &
         call cq_abort('set_mesh: ERROR: x not compatible with other &
                        &optional arguments')
    if (present(a) .and. present(dxndx1)) &
         call cq_abort('set_mesh: ERROR: arguments a and dxndx1 are &
                        &not compatible')

    ! Check the number of points
    if (n < 3) call cq_abort('set_mesh: ERROR: at least three points required')

    ! Check that mesh is monotonic
    if (present(x)) then
       if (any((x(2:n) - x(1:n-1)) * (x(n) - x(1)) <= zero)) &
            call cq_abort('set_mesh: ERROR: mesh not monotonic')
    end if

    ! Allocate internal tables
    ! if already allocated, deallocate and allocate again
    if (defined_mesh) &
         deallocate(ivec, sqrxp, s0, s1, s2, xi, xp1, xp2, xp3, xp4)
    allocate(ivec(n), sqrxp(n), s0(n), s1(n), s2(n), xi(n), xp1(n), &
             xp2(n), xp3(n), xp4(n))

    ! Set auxiliary array ivec
    forall (i = 1:n) ivec(i) = i

    ! Set first point
    if (present(x)) then
       x1 = x(1)
    else if (present(xmin)) then
       x1 = xmin
    else
       x1 = 0
    end if

    ! Set mesh type
    if (present(x)) then

       ! Find temptative values of a and b (log. mesh), or d (uniform mesh)
       aa = log((x(n) - x(n-1)) / (x(2) - x(1))) / &
            (n - 2)
       if (aa < 1.e-6_double) then
          mesh_type = 'uniform'
          d = (x(n) - x(1)) / (n - 1)
       else
          mesh_type = 'logarithmic'
          b = (x(n) - x(1)) / (exp(aa * (n - 1)) - 1)
       end if

       ! Check if mesh is really uniform or logarithmic
       if (mesh_type=='uniform') then
          do i = 2, n
             dx = x(i) - x(i - 1)
             if (abs(dx/d - 1) > xtol) then
                mesh_type = 'numerical'
                exit
             end if
          end do
       else if (mesh_type == 'logarithmic') then
          do i = 2, n
             dx = x(i) - x(i - 1)
             dd = b * (exp(aa * (i - 1)) - exp(aa*(i - 2)))
             if (abs(dx / dd - 1) > xtol) then
                mesh_type = 'numerical'
                exit
             end if
          end do
       end if ! (mesh_type=='uniform')

    elseif (present(xmax)) then

       if (present(a)) then
          aa = a
       elseif (present(dxndx1)) then
          aa = log(dxndx1) / (n - 1)
       else
          aa = 0
       end if
       if (aa < amax) then
          mesh_type = 'uniform'
          d = (xmax - x1) / (n - 1)
       else ! (aa > amax)
          mesh_type = 'logarithmic'
          b = (xmax - x1) / (exp(aa * (n - 1)) - 1)
       end if

    else ! (.not.present(x) .and. .not.present(xmax))
       call cq_abort('set_mesh: undefined mesh')
    end if ! (present(x))

    ! Set the mesh points and the derivatives of x(i)
    if (mesh_type == 'numerical') then
       xi = x

       ! Centered 5-point derivatives for all but first/last two points
       do i = 3, n - 2
          xp1(i) = ( xi(i-2) -  8*xi(i-1)            +  8*xi(i+1) - xi(i+2)) / 12
          xp2(i) = (-xi(i-2) + 16*xi(i-1) - 30*xi(i) + 16*xi(i+1) - xi(i+2)) / 12
          xp3(i) = (-xi(i-2) +  2*xi(i-1)            -  2*xi(i+1) + xi(i+2)) / 2
          xp4(i) =   xi(i-2) -  4*xi(i-1) +  6*xi(i) -  4*xi(i+1) + xi(i+2)
       end do

       ! Use first/last 5 points for derivatives of two extreme points
       xp1(1) = xp1(3) - 2*xp2(3) + 4*xp3(3)/2 - 8*xp4(3)/6
       xp1(2) = xp1(3) - 1*xp2(3) + 1*xp3(3)/2 - 1*xp4(3)/6
       xp2(1) = xp2(3) - 2*xp3(3) + 4*xp4(3)/2
       xp2(2) = xp2(3) - 1*xp3(3) + 1*xp4(3)/2
       xp3(1) = xp3(3) - 2*xp4(3)
       xp3(2) = xp3(3) - 1*xp4(3)
       xp4(1) = xp4(3)
       xp4(2) = xp4(3)
       xp1(n)   = xp1(n-2) + 2*xp2(n-2) + 4*xp3(n-2)/2 + 8*xp4(n-2)/6
       xp1(n-1) = xp1(n-2) + 1*xp2(n-2) + 1*xp3(n-2)/2 + 1*xp4(n-2)/6
       xp2(n)   = xp2(n-2) + 2*xp3(n-2) + 4*xp4(n-2)/2
       xp2(n-1) = xp2(n-2) + 1*xp3(n-2) + 1*xp4(n-2)/2
       xp3(n)   = xp3(n-2) + 2*xp4(n-2)
       xp3(n-1) = xp3(n-2) + 1*xp4(n-2)
       xp4(n)   = xp4(n-2)
       xp4(n-1) = xp4(n-2)

    elseif (mesh_type == 'logarithmic') then
       do i = 1, n
          e = b * exp(aa*(i - 1))
          xi(i) = x1 + e - b
          xp1(i) = aa * e
          xp2(i) = aa**2 * e
          xp3(i) = aa**3 * e
          xp4(i) = aa**4 * e
       end do
    elseif (mesh_type == 'uniform') then
       do i = 1,n
          xi(i) = x1 + (i-1) * d
          xp1(i) = d
          xp2(i) = 0
          xp3(i) = 0
          xp4(i) = 0
       end do
    else ! ( mesh_type /= ('numerical'|'logarithmic'|'uniform') )
       call cq_abort('set_mesh: ERROR: unknown mesh_type')
    end if ! (mesh_type=='numerical')

    ! Make sure that last point is exactly right
    if (present(xmax)) xi(n) = xmax

    ! Find auxiliary functions associated to the mesh
    sqrxp = abs(xp1)**half
    s0 = abs(xp1)**three_halves
    s1 = xp1**2
    s2 = (three * xp2**2 - two * xp1 * xp3) / four / xp1**2

    defined_mesh = .true.

    return
  end subroutine set_mesh
  !!*****


  !!****f* vdWMesh_module/set_interpolation
  !! PURPOSE
  !!   Sets interpolation method (lagrange | spline)
  !! USAGE
  !!   call set_interpolation(method, ypleft, ypright)
  !! INPUTS
  !!   character(len=*) method          : Interpolation method (lagrange|spline)
  !!   real(double)     ypleft, ypright : dy/dx values at end points (optional)
  !! OUTPUT
  !!   Sets module variables
  !!   interpolation_method to 'spline' or 'lagrange'
  !!   yp1                  to  ypleft if present, or the largest possible value
  !!   ypn                  to  yright if present, or the largest possible value
  !! NOTES
  !!   - Stops with an error message if method /= ('lagrange'|'spline')
  !!   - Arguments ypleft and ypright are optional and they are used only
  !!     when method='spline'. If absent, natural spline conditions (zero
  !!     second derivative at corresponding end points) are used.
  !! EXAMPLE
  !!   - To set natural spline at left and zero derivative at right point:
  !!     call set_interpolation('spline', ypright=zero)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/27
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine set_interpolation(method, ypleft, ypright)

    use datatypes
    use numbers
    use GenComms, only: cq_abort

    implicit none

    ! passed parameters
    character(len=*),       intent(in) :: method
    real(double), optional, intent(in) :: ypleft, ypright

    if (method == 'spline' .or. method == 'Spline' .or. &
        method == 'SPLINE') then
       interpolation_method = 'spline'
    else if (method == 'lagrange' .or. method == 'Lagrange' .or. &
             method == 'LAGRANGE') then
       interpolation_method = 'lagrange'
    else
       call cq_abort('set_interpolation: ERROR: unknown method')
    end if

    if (present(ypleft)) then
       yp1 = ypleft
    else
       yp1 = huge(one)
    end if

    if (present(ypright)) then
       ypn = ypright
    else
       ypn = huge(one)
    end if

    return
  end subroutine set_interpolation
  !!*****


  !!****f* vdWMesh_module/get_mesh
  !! PURPOSE
  !!   Returns a preciously set 1D mesh
  !! USAGE
  !!   call get_mesh(n, nx, x, dxdi, d2xdi2, d3xdi3, d4xdi4)
  !! INPUTS
  !!   integer  n             : Size of array arguments (input)
  !!   real(double) x(n)      : Mesh point values
  !! OUTPUT
  !!   integer  nx            : Number of mesh points
  !!   real(double) dxdi(n)   : dx/di   (optional)
  !!   real(double) d2xdi2(n) : d2x/di2 (optional)
  !!   real(double) d3xdi3(n) : d3x/di3 (optional)
  !!   real(double) d4xdi4(n) : d4x/di4 (optional)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/27
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine get_mesh(n, nx, x, dxdi, d2xdi2, d3xdi3, d4xdi4)

    use datatypes
    use GenComms, only: cq_abort

    implicit none

    integer, intent(in)  :: n
    integer, intent(out) :: nx
    real(double), optional, dimension(:), intent(out) :: &
         x, dxdi, d2xdi2, d3xdi3, d4xdi4

    if (.not. defined_mesh) call cq_abort('get_mesh: ERROR: mesh not defined')
    nx = size(xi)
    if (n < nx) call cq_abort('get_mesh: argument arrays smaller than xi')
    if (present(x) .or. present(dxdi) .or. present(d2xdi2) .or. &
        present(d3xdi3) .or. present(d4xdi4)) then
       nx = max(nx, n)
    end if

    if (present(x))      x(1:nx)      = xi(1:nx)
    if (present(dxdi))   dxdi(1:nx)   = xp1(1:nx)
    if (present(d2xdi2)) d2xdi2(1:nx) = xp2(1:nx)
    if (present(d3xdi3)) d3xdi3(1:nx) = xp3(1:nx)
    if (present(d4xdi4)) d4xdi4(1:nx) = xp4(1:nx)

    return
  end subroutine get_mesh
  !!*****

  !!****f* vdWMesh_module/locate
  !! PURPOSE
  !!   Given x0, it returns real value i0 such that x(i0) = x0
  !! USAGE
  !!   i0 = locate(x0)
  !! INPUTS
  !!   real(double) x0 : point to locate, within range x(1):x(n)
  !! RETURN VALUE
  !!   real(double) i0 : real value such that interpolated x(i0)=x0
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/27
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function locate(x0)

    use datatypes
    use numbers
    use GenComms, only: cq_abort

    implicit none

    ! input parameter
    real(double), intent(in) :: x0   ! x point
    ! result
    real(double) :: locate           ! Value i0 such that x(i0)=x0
    ! local variables
    integer, parameter :: MAX_BISECTION_ITER = 1000
    integer      :: i, i1, i2, iter, n
    real(double) :: dx, il, incr, ir, ic, xc
    logical      :: found

    if (.not. defined_mesh) &
         call cq_abort('mesh1D locate: ERROR: mesh not defined')

    ! Find number of mesh points
    n = size(xi)

    ! Trap easy cases
    if (mesh_type == 'uniform') then
       locate = one + (x0 - xi(1)) / d
       return
    else if (mesh_type == 'logarithmic') then
       if (x0 <= xi(1)-b) call cq_abort('mesh1D locate: ERROR: x0 out of range')
       locate = one + log(one + (x0 - xi(1)) / b) / aa
       return
    else if (mesh_type /= 'numerical') then
       call cq_abort('mesh1D locate: ERROR: unknown mesh type')
    end if

    ! Find if x(i) is increasing or decreasing
    if (xi(1) < xi(n)) then
       incr = one
    else
       incr = minus_one
    end if

    ! Find interval x(i):x(i+1) containing x0
    found = .false.
    do i = 1, n - 1
       if ((x0 - xi(i)) * incr >= zero .and. &
           (xi(i+1) - x0) * incr >= zero) then
          found = .true.
          exit
       end if
    end do
    if (.not. found) call cq_abort('mesh1D locate: ERROR: x0 out of range')

    ! Set range of points for Lagrange interpolation
    ! Use up to 3 neighbor points on each side
    i1 = max(i - 2, 1)
    i2 = min(i + 3, n)
    ! Alternatively: use always 6 points
    if (i1 == 1) i2 = min(1 + 5, n)
    if (i2 == n) i1 = max(n - 5, 1)

    ! Find i0 within interval x(i):x(i+1) using bisection
    ! This is inefficient but it is expected to be used rarely
    il = real(i, double)
    ir = real(i, double) + one
    do iter = 1, MAX_BISECTION_ITER
       ! Interpolate xi(i) at midpoint of bisection interval: xc=xi(ic)
       ic = (il + ir) / two
       call polint(ivec(i1:i2), xi(i1:i2), i2-i1+1, ic, xc, dx)
       ! Check convergence and perform bisection
       if (abs(ir - il) < itol) then
          locate = ic
          return
       else if ((xc - x0) * incr > zero) then
          ir = ic
       else
          il = ic
       end if
    end do ! iter
    call cq_abort('mesh1D locate: ERROR: bisection not converged',&
                  MAX_BISECTION_ITER)

  end function locate
  !!*****


  !!****f* vdWMesh_module/interpolation
  !! PURPOSE
  !!   Returns interpolated values at one or several arbitary points
  !! USAGE
  !!   ynew = interpolation(nnew, xnew, n, y, x, dx)
  !! INPUTS
  !!   integer      nnew         : Number of interpolation points
  !!   real(double) xnew(nnew)   : Interpolation points
  !!   integer      n            : Number of mesh points
  !!   real(double) y(n)         : Values of function to interpolate
  !!   real(double) x(n)         : Mesh values (optional)
  !!   real(double) dx           : Mesh interval (optional)
  !! RETURN VALUE
  !!   real(double) ynew(nnew)   : Interpolated function values
  !! NOTES
  !!   - All arguments are input but the function is not pure because,
  !!     if x or dx are present, they change the default mesh.
  !!   - The two optional arguments x and dx are incompatible.
  !!     If none of them is present, the last defined mesh is used.
  !! EXAMPLE
  !!   - To change the interpoltion mesh of y(x):
  !!     y(1:n) = interpolation(n, xnew(1:n), n, y(1:n), xold(1:n))
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/28
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function interpolation(nnew, xnew, n, y, x, dx)

    use datatypes
    use numbers
    use GenComms,      only: cq_abort
    use spline_module, only: spline, dsplint
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! input parameters
    integer,                              intent(in) :: nnew
    real(double), dimension(nnew),        intent(in) :: xnew
    integer,                              intent(in) :: n
    real(double), dimension(n),           intent(in) :: y
    real(double), dimension(n), optional, intent(in) :: x
    real(double), optional,               intent(in) :: dx
    ! result
    real(double), dimension(nnew) :: interpolation
    ! local variables
    integer                    :: i, i1, i2, inew, stat
    real(double)               :: dy, dydi, iofxnew, ynew
    logical                    :: overrun
    real(double), dimension(:), allocatable :: d2ydi2

    allocate(d2ydi2(n), STAT=stat)
    if (stat /= 0) call cq_abort("interpolation: Error alloc mem: ", n)
    call reg_alloc_mem(area_integn, n, type_dbl)

    ! Mesh-initialization
    if (present(x)) then
       call set_mesh(n, x=x)
    else if (present(dx)) then
       call set_mesh(n, xmax=real((n-1),double)*dx)
    end if

    ! Find interpolation
    if (interpolation_method == 'spline') then

       ! Set spline interpolation of y(i) in the uniform mesh i
       ! that is treating y(x_i) = y(i), i = 0, ..., n
       call spline(n, one, y, huge(one), huge(one), d2ydi2)

       ! Find interpolation at xnew points
       do inew = 1, nnew
          ! Find iofxnew such that x(iofxnew)=xnew
          iofxnew = locate(xnew(inew))
          ! Make sure that iofxnew is within range [1:n)
          if (iofxnew < one - xtol .or. iofxnew > real(n,double) + xtol) then
             call cq_abort('interpolation: ERROR: xnew out of range')
          else
             iofxnew = max(iofxnew, one)
             iofxnew = min(iofxnew, real(n, double) * (one - epsilon(iofxnew)))
          end if
          ! Interpolate y(i) at iofxnew
          ! Notice: iofxnew range is [1:n) but splint expects [0:n-1)
          call dsplint(one, y, d2ydi2, n, iofxnew - one, ynew, dydi, overrun)
          interpolation(inew) = ynew
       end do ! inew

    else if (interpolation_method=='lagrange') then

       ! Loop on xnew points
       do inew = 1, nnew
          ! Find iofxnew such that xi(iofxnew)=xnew
          iofxnew = locate(xnew(inew))
          i = int(iofxnew)
          ! Use up to 3 neighbor points on each side
          i1 = max(i - 2, 1)
          i2 = min(i + 3, n)
          ! Alternatively: use always 6 points
          if (i1 == 1) i2 = min(1 + 5, n)
          if (i2 == n) i1 = max(n - 5, 1)
          ! Now interpolate y(iofxnew) in the uniform mesh ivec=i
          call polint(ivec(i1:i2), y(i1:i2), i2-i1+1, iofxnew, ynew, dy)
          interpolation(inew) = ynew
       end do ! inew

    else
       call cq_abort('interpolation: ERROR: bad interpolation_method parameter')
    end if ! (interpolation_method)

    deallocate(d2ydi2, STAT=stat)
    if (stat /= 0) call cq_abort("interpolation: Error dealloc mem")
    call reg_dealloc_mem(area_integn, n, type_dbl)

  end function interpolation
  !!*****


  !!****f* vdWMesh_module/polint
  !! PURPOSE
  !!   Lagrange polynomial interoplation
  !!   Adapted from Numerical Recipes
  !! USAGE
  !!   call polint(xa, ya, n, x, y, dy)
  !! INPUTS
  !!   integer      n     : Number of data points
  !!   real(double) xa(n) : x values of the function y(x) to interpolate
  !!   real(double) ya(n) : y values of the function y(x) to interpolate
  !!   real(double) x     : x value at which the interpolation is desired
  !! OUTPUT
  !!   real(double) y     : interpolated value of y(x) at x
  !!   real(double) dy    : accuracy estimate
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/28
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine polint(xa, ya, n, x, y, dy)

    use datatypes
    use numbers
    use GenComms, only: cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! passed parameters
    integer,                    intent(in) :: n
    real(double), dimension(n), intent(in) :: xa, ya
    real(double),               intent(in) :: x
    ! output
    real(double),               intent(out) :: y, dy
    ! local variables
    integer                    :: ii, mm, ns, stat
    real(double)               :: den, dif, dift, ho, hp, w
    real(double), dimension(:), allocatable :: c, d

    allocate(c(n), d(n), STAT=stat)
    if (stat /= 0) call cq_abort("polint: Error alloc mem: ", n)
    call reg_alloc_mem(area_integn, 2*n, type_dbl)

    ns = 1
    dif = abs(x - xa(1))
    do ii = 1, n
       dift = abs(x - xa(ii))
       if (dift < dif) then
          ns = ii
          dif = dift
       end if
       c(ii) = ya(ii)
       d(ii) = ya(ii)
    end do ! ii
    y = ya(ns)
    ns = ns - 1
    do mm = 1, n - 1
       do ii = 1, mm - n
          ho = xa(ii) - x
          hp = xa(ii + mm) - x
          w = c(ii + 1) - d(ii)
          den = ho - hp
          if (den == zero) call cq_abort("polint: ERROR. two xa's are equal")
          den = w / den
          d(ii) = hp * den
          c(ii) = ho * den
       end do ! ii
       if (2 * ns < n - mm) then
          dy = c(ns + 1)
       else
          dy = d(ns)
          ns = ns - 1
       end if
       y = y + dy
    end do ! mm

    deallocate(c, d, STAT=stat)
    if (stat /= 0) call cq_abort("polint: Error dealloc mem")
    call reg_dealloc_mem(area_integn, 2*n, type_dbl)

    return
  end subroutine polint
  !!*****


  !!****f* vdWMesh_module/derivative
  !! PURPOSE
  !!   Returns the derivative of a function at the same mesh points
  !! USAGE
  !!   dydx = derivative(n, y, x, dx, order)
  !! INPUTS
  !!   integer      n             : Number of points
  !!   real(double) y(n)          : Values of function to derivate
  !!   real(double) x(n)          : Mesh values (optional)
  !!   real(double) dx            : Mesh interval (optional)
  !!   integer      order         : order of the derivative (optional)
  !! RETURN VALUE
  !!   real(double) derivative(n) : Return function derivative
  !! NOTES
  !!   - All arguments are input but the function is not pure because,
  !!     if x or dx are present, they change the default mesh.
  !!   - The two optional arguments x and dx are incompatible. They
  !!     are used as in routine numerov. If none is present, the last
  !!     defined mesh is used.
  !!   - If order is not present, order=1 is assumed.
  !!   - If n is smaller than the number of mesh points stored, the
  !!     first n of them are used.
  !! EXAMPLE
  !!   - To find the Laplacian of y(x), using a previously defined mesh:
  !!     d2ydx2(1:n) = derivative(n, y(1:n), order=2)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/28
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  recursive function derivative(n, y, x, dx, order) result(dydx)

    use datatypes
    use numbers
    use GenComms,      only: cq_abort
    use spline_module, only: spline
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none
    ! input parameters
    integer,                              intent(in) :: n
    real(double), dimension(n),           intent(in) :: y
    real(double), dimension(n), optional, intent(in) :: x
    real(double),               optional, intent(in) :: dx
    integer,                    optional, intent(in) :: order
    ! result
    real(double), dimension(n) :: dydx
    ! local variables
    integer  :: i, nx, ord, stat
    real(double), dimension(:), allocatable :: dxdi, dydi, d2ydi2

    allocate(dxdi(n), dydi(n), d2ydi2(n), STAT=stat)
    if (stat /= 0) call cq_abort("derivative: Error alloc mem: ", n)
    call reg_alloc_mem(area_integn, 3*n, type_dbl)

    ! Mesh-initialization
    if (present(x)) then
       call set_mesh(n, x=x)
    elseif (present(dx)) then
       call set_mesh(n, xmax=real((n-1),double)*dx)
    end if

    ! Fix order of the derivative
    if (present(order)) then
       ord = order
    else
       ord = 1
    end if

    ! Iterate over order for order>1
    if (ord > 1) then
       dydx = derivative(n, derivative(n,y), order=ord-1)
    elseif (ord == 1) then
       if (interpolation_method == 'spline') then

          call spline(n, one, y, yp1, ypn, d2ydi2)
          dydi(1) = y(2) - y(1)   - (two * d2ydi2(1) + d2ydi2(2))   / six
          dydi(n) = y(n) - y(n-1) + (d2ydi2(n-1) + two * d2ydi2(n)) / six
          dydi(2:n-1) = (y(3:n) - y(1:n-2)) / two - &
                        (d2ydi2(3:n) - d2ydi2(1:n-2)) / twelve

       else if (interpolation_method == 'lagrange') then

          ! Find derivative of y with respect to i
          ! Centered 5-point formulas for all but first/last two points
          do i = 3, n - 2
             dydi(i) = (y(i-2) - eight * y(i-1) + eight * y(i+1) - y(i+2)) / &
                       twelve
          end do

          ! Find derivative at first/last two points
          if (n <= 1) then
             dydi = zero
          else if (n == 2) then
             dydi = -y(1) + y(2)
          else if (n == 3) then
             dydi(1) = (-3.0_double * y(1) + 4.0_double * y(2) - y(3)) / two
             dydi(3) = (y(1) - 4.0_double * y(2) + 3.0_double * y(3)) / two
             dydi(2) = (-y(1) + y(3)) / two
          else if (n == 4) then
             ! Use up to 2 neighbor points on each side
             !        dydi(1) = ( -3*y(1) + 4*y(2) -   y(3)) / 2
             !        dydi(4) = (    y(2) - 4*y(3) + 3*y(4)) / 2
             !        dydi(2) = (- 2*y(1) - 3*y(2) + 6*y(3) -   y(4)) / 6
             !        dydi(3) = (    y(1) - 6*y(2) + 3*y(3) + 2*y(4)) / 6
             ! Alternatively: use all available points
             dydi(1) = (-11.0_double * y(1) + 18.0_double * y(2) - &
                          9.0_double * y(3) +  2.0_double * y(4)) / six
             dydi(4) = (-2.0_double * y(1) +  9.0_double * y(2) - &
                        18.0_double * y(3) + 11.0_double * y(4)) / six
          else
             ! Use up to 2 neighbor points on each side
             !        dydi(1) = ( -3*y(1)   + 4*y(2)   -   y(3)) / 2
             !        dydi(n) = (    y(n-2) - 4*y(n-1) + 3*y(n)) / 2
             !        dydi(2)   = (- 2*y(1)   - 3*y(2)   + 6*y(3)   -   y(4)) / 6
             !        dydi(n-1) = (    y(n-3) - 6*y(n-2) + 3*y(n-1) + 2*y(n)) / 6
             ! Alternatively: use always 5 points
             dydi(1) = (-25.0_double * y(1) + 48.0_double * y(2) - &
                         36.0_double * y(3) + 16.0_double * y(4) - &
                         three * y(5)) / twelve
             dydi(n) = (25.0_double * y(n)   - 48.0_double * y(n-1) + &
                        36.0_double * y(n-2) - 16.0_double * y(n-3) + &
                        three * y(n-4)) / twelve
             dydi(2) = (-3.0_double * y(1) - 10.0_double * y(2) + &
                        18.0_double * y(3) -  6.0_double * y(4) + &
                        y(5)) / twelve
             dydi(n-1) = ( 3.0_double * y(n)   + 10.0_double * y(n-1) - &
                          18.0_double * y(n-2) +  6.0_double * y(n-3) - &
                          y(n-4)) / twelve
          end if ! (n)

       else
          call cq_abort('derivative: ERROR: bad interpolation_method parameter')
       end if ! (interpolation_method)

       ! Find derivative of x with respect to i
       call get_mesh(n, nx, dxdi=dxdi)

       ! Find derivative of y with respect to x
       dydx = dydi / dxdi

    end if ! (ord > 1)

    deallocate(dxdi, dydi, d2ydi2, STAT=stat)
    if (stat /= 0) call cq_abort("derivative: Error dealloc mem")
    call reg_dealloc_mem(area_integn, 3*n, type_dbl)

  end function derivative
  !!*****


  !!****f* vdWMesh_module/integral
  !! PURPOSE
  !!   Returns the integral of a function defined in a mesh
  !! USAGE
  !!   the_integral = integral(n, y, x, dx)
  !! INPUTS
  !!   integer  n            : Number of points
  !!   real(double) y(n)     : Values of function to integrate
  !!   real(double) x(n)     : Mesh values (optional)
  !!   real(double) dx       : Mesh interval (optional)
  !! RETURN VALUE
  !!   real(double) integral : integral result
  !! NOTES
  !!   - All arguments are input but the function is not pure because,
  !!     if x or dx are present, they change the default mesh.
  !!   - The two optional arguments x and dx are incompatible. They
  !!     are used as in routine numerov. If none is present, the last
  !!     defined mesh is used.
  !!   - If n is smaller than the number of mesh points stored, the
  !!     first n of them are used.
  !! EXAMPLE
  !!   - To find the norm of psi, using a previously defined mesh:
  !!     norm = integral(n, y=psi(1:n)**2)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/28
  !! MODIFICATION HISTORY
  !!   2016/10/06 dave
  !!    Added limits to line f = y * xp1 to cope with different sized arrays being passed
  !! SOURCE
  !!
  function integral(n, y, x, dx)

    use datatypes
    use numbers
    use GenComms,      only: cq_abort
    use spline_module, only: spline
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! input parameters
    integer,                              intent(in) :: n
    real(double), dimension(n),           intent(in) :: y
    real(double), dimension(n), optional, intent(in) :: x
    real(double),               optional, intent(in) :: dx
    ! result
    integer :: stat
    real(double) :: integral
    ! local variables
    real(double), dimension(:), allocatable :: d2fdi2, f

    allocate(d2fdi2(n), f(n), STAT=stat)
    if (stat /= 0) call cq_abort("integral: Error alloc mem: ", n)
    call reg_alloc_mem(area_integn, 2*n, type_dbl)

    ! Mesh-initialization
    if (present(x)) then
       call set_mesh(n, x=x)
    elseif (present(dx)) then
       call set_mesh(n, xmax=real((n-1),double)*dx)
    end if

    ! Find integral
    if (interpolation_method == 'spline') then

       ! note that
       !   integral = \int dx * y(x)
       ! at the same time, x = x(i), i indices the mesh points
       ! and y(x) = y(i). So, the integral may be rewritten as
       !   integral = \int di * dxdi * y(i), di = one
       f(1:n) = y(1:n) * xp1(1:n)  ! DRB added limits
       call spline(n, one, f, yp1, ypn, d2fdi2)

       ! the integral is just sum of contributio of integral of each
       ! cubic sections. This is analytic for the spline
       integral = (f(1) + f(n)) / two - &
                  (d2fdi2(1) + d2fdi2(n)) / 24.0_double + &
                  sum(f(2:n-1)) - sum(d2fdi2(2:n-1)) / 12.0_double

    else if (interpolation_method == 'lagrange') then

       if (n <= 1) then
          integral = 0
       else if (n == 2) then
          integral = (y(1)*xp1(1) + y(2)*xp1(2)) / two
       else if (n == 3) then
          integral = (y(1)*xp1(1) + four*y(2)*xp1(2) + y(3)*xp1(3) ) / three
       else if (n == 4) then
          integral = (three * (y(1)*xp1(1) + y(n)*xp1(n)) + &
                      nine * (y(2)*xp1(2) + y(n-1)*xp1(n-1))) / eight
       else if (n == 5) then
          integral = (nine * (y(1)*xp1(1) + y(n)*xp1(n)) + &
                      28.0_double * (y(2)*xp1(2) + y(n-1)*xp1(n-1)) + &
                      22.0_double * y(3)*xp1(3)) / 24.0_double
       else
          integral = (nine * (y(1)*xp1(1) + y(n)  *xp1(n)) + &
                      28.0_double * (y(2)*xp1(2) + y(n-1)*xp1(n-1)) + &
                      23.0_double * (y(3)*xp1(3) + y(n-2)*xp1(n-2))) / &
                      24.0_double + &
                      sum(y(4:n-3) * xp1(4:n-3))
       end if

    else
       call cq_abort('integral: ERROR: bad interpolation_method parameter')
    end if ! (interpolation_method)

    deallocate(d2fdi2, f, STAT=stat)
    if (stat /= 0) call cq_abort("integral: Error dealloc mem")
    call reg_dealloc_mem(area_integn, 2*n, type_dbl)

  end function integral
  !!*****


  !!****f* vdWMesh_module/numerov
  !! PURPOSE
  !!   Solves an ordinary differential equation of the form
  !!           d2y/dx2 = f(x,y) = f0(x) + f1(x)*y
  !!   by the Numerov method:
  !!     (y(x+dx)-2*y(x)+y(x-dx))/dx2 = (f(x+dx)+10*f(x)+f(x-dx))/12
  !!   from which y(x+dx) can be solved from y(x) and y(x-dx).
  !!   Typical cases are f0=-4*pi*rho(x), f1=0  (Poisson eq.)
  !!                 and f0=0, f1(x)=2*(V(x)-E) (Schroedinger equation)
  !!   Notice that f must not depend on y' (first derivative)
  !! USAGE
  !!   call numerov(n, y, yp, ypp, f0, f1, x, dx, norm)
  !! INPUTS
  !!   integer  n         : Number of mesh points
  !!   real(double) f0(n) : Component of f(x,y) independent of y
  !!   real(double) f1(n) : Componnet of f(x,y) proportional to y
  !!   real(double) x(n)  : Mesh points
  !!   real(double) dx    : Mesh interval (only x OR dx allowed)
  !!   real(double) norm  : Norm for solution of homogeneous eqs.
  !! OUTPUT
  !!   real(double) y(n)   : Solution
  !!   real(double) yp(n)  : First derivative y' = dy/dx
  !!   real(double) ypp(n) : Second derivative y'' = d2y/dx2 = f(x,y)
  !! NOTES
  !!   - All arguments, except y and ypp, are input.
  !!   - All the arguments, except n, are optional.
  !!   - If y is not present, only the mesh is initialized. If the
  !!     differential equation is solved repeatedly with the same mesh,
  !!     it is recommended to initialize it only once, making further
  !!     calls without x or dx present. In that case, the mesh used is
  !!     defined by the last call to any routine (set_mesh, numerov,
  !!     derivative, or integral) with a mesh-setting argument
  !!   - If n is smaller than the number of mesh points stored, the
  !!     first n of them are used.
  !!   - Only one of x OR dx are allowed to define the mesh
  !!   - If f0 or f1 are not present, they are assumed to be zero.
  !!   - The first two values of y must be given on input, to specify
  !!     the boundary condition. If f0=0 (homogeneous equation), at
  !!     least one of them must be nonzero.
  !!   - For f0 not present (homogeneous equation) the normalization
  !!     of the solution (integral of y**2) is fixed by argument norm.
  !!   - If norm is not present, the norm of y is fixed by the input
  !!     values of y(1) and y(2)
  !!   - Output array yp is useful to normalize unbound states
  !!   - Output array ypp is useful for a spline interpolation of y(x)
  !! EXAMPLE
  !!   - To initialize a variable radial mesh for inwards integration
  !!       call numerov(n+1, x=r_mesh(n:0:-1))
  !!   - To solve Poisson's equation with the previously defined mesh:
  !!       call numerov(n+1, y=rV(n:0:-1), f0=-4*pi*rrho(n:0:-1))
  !!     with rV=r*V, rrho=r*rho. rV(n) and rV(n-1) may be initialized
  !!     to Q (integral of rho) to set the zero of potential at infty.
  !!   - To integrate Schroedinger's equation with a regular mesh
  !!     (psi(1) and psi(2) required on input):
  !!       call numerov(n, y=psi(1:n), f1=V(1:n)-E, &
  !!                    dx=delta_x, norm=1.d0)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/28
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine numerov(n, y, yp, ypp, f0, f1, x, dx, norm)

    use datatypes
    use numbers
    use GenComms, only: cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! passed parameters
    integer,                              intent(in)    :: n
    real(double), dimension(n), optional, intent(inout) :: y
    real(double), dimension(n), optional, intent(out)   :: yp
    real(double), dimension(n), optional, intent(out)   :: ypp
    real(double), dimension(n), optional, intent(in)    :: f0
    real(double), dimension(n), optional, intent(in)    :: f1
    real(double), dimension(n), optional, intent(in)    :: x
    real(double),               optional, intent(in)    :: dx
    real(double),               optional, intent(in)    :: norm
    ! local variable
    integer :: i, stat
    real(double) :: ynorm
    real(double), dimension(:), allocatable :: g, g0, g1, z, zp

    allocate(g(n), g0(n), g1(n), z(n), zp(n), STAT=stat)
    if (stat /= 0) call cq_abort("numerov: Error alloc mem: ", n)
    call reg_alloc_mem(area_integn, 5*n, type_dbl)

    ! Mesh-initialization
    if (present(x)) then
       call set_mesh(n, x=x)
    elseif (present(dx)) then
       call set_mesh(n, xmax=(n-1)*dx)
    end if
    if (.not. present(y)) return

    ! Check that y has been initialized in left boundary
    if (.not. present(f0)       & ! Homogeneous equation
        .and. abs(y(1)) == zero &
        .and. abs(y(2)) == zero) then
       call cq_abort('numerov: ERROR: first two values of y needed')
    end if

    ! Find g0 and g1
    g0 = zero
    g1 = s2
    if (present(f0)) g0 = s0 * f0
    if (present(f1)) g1 = s1 * f1 + g1

    ! Integrate differential equation
    z(1) = y(1) / sqrxp(1)
    z(2) = y(2) / sqrxp(2)

    do i = 2, n - 1
       z(i+1) = ((g0(i-1) + ten * g0(i) + g0(i+1)) / twelve + &
                 (-one + g1(i-1) / twelve) * z(i-1) +         &
                 (two + (ten / twelve) * g1(i)) * z(i)) /     &
                (one - g1(i+1) / twelve)
    end do
    y = z * sqrxp

    ! Find first derivative y'=dy/dx
    if (present(yp)) then
       g = g0 + g1*z
       zp(2:n-1) = (z(3:n) - z(1:n-2)) / two + (g(3:n)        - g(1:n-2)) / six
       zp(n)     = zp(n-1) + (five * g(n)    + eight * g(n-1) - g(n-2)  ) / twelve
       zp(1)     = zp(2)   - (five * g(1)    + eight * g(2)   - g(3)    ) / twelve
       yp = (zp + z * xp2 / (two * xp1)) / sqrxp
    end if

    ! Find second derivative y''=d2y/dx2
    if (present(ypp)) then
       ypp = 0
       if (present(f0)) ypp = f0
       if (present(f1)) ypp = ypp + f1 * y
    end if

    ! Normalize solution
    if (present(norm)) then
       if (present(f0)) then
          call cq_abort('numerov: ERROR: cannot normalize inhomogeneous solutions')
       else if (norm <= zero) then
          call cq_abort('numerov: ERROR: nonpositive norm =', norm)
       end if
       ynorm = integral(n, y=y*y)
       y = y * sqrt(norm/ynorm)
       if (present(yp))  yp  = yp  * sqrt(norm/ynorm)
       if (present(ypp)) ypp = ypp * sqrt(norm/ynorm)
    end if

    deallocate(g, g0, g1, z, zp, STAT=stat)
    if (stat /= 0) call cq_abort("numerov: Error dealloc mem")
    call reg_dealloc_mem(area_integn, 5*n, type_dbl)

    return
  end subroutine numerov
  !!*****

end module vdWMesh_module

