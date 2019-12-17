module mesh

  use datatypes
  
  implicit none
  
  integer :: nmesh

  ! Default values
  real(double) :: alpha = log(1.012_double)
  real(double) :: beta = 0.0005_double
  real(double) :: mesh_z
  real(double) :: delta_r_reg
  real(double), allocatable, dimension(:) :: rr, rr_squared, drdi, drdi_squared, sqrt_drdi, rmesh_reg
  integer, parameter :: hamann = 1
  integer, parameter :: siesta = 2
  integer :: mesh_type
  ! DRB 2019/06/04 Changing this so that the polynomial is of order 2*n_poly_lag+1
  integer :: n_poly_lag = 3 ! User adjustable in the long term
  integer :: n_poly = 3 ! Should this be user adjustable ? Actually gives 2*n_poly+1

  logical, save :: flag_allocated_mesh = .false.
  
contains

  subroutine make_mesh(species)

    use datatypes
    use numbers
    use pseudo_tm_info, ONLY: pseudo
    use pseudo_atom_info, ONLY: local_and_vkb
    use global_module, ONLY: iprint
    
    implicit none

    ! Passed variables
    integer :: species
    
    ! Local variables
    integer :: i
    real(double) :: alpha_hamann

    !
    ! Initialise and allocate
    !
    if(flag_allocated_mesh) then
       deallocate(rr,rr_squared,drdi,drdi_squared,sqrt_drdi)
    end if
    allocate(rr(nmesh),rr_squared(nmesh),drdi(nmesh),drdi_squared(nmesh),sqrt_drdi(nmesh))
    flag_allocated_mesh = .true.
    rr = zero
    rr_squared = zero
    drdi = zero
    drdi_squared = zero
    sqrt_drdi = zero
    !
    ! Set mesh parameters
    !
    mesh_z = real(pseudo(species)%z, double)
    if(mesh_type==hamann) then
       ! Check consistency of alpha and grid that has been read
       alpha_hamann = 0.01_double*log(local_and_vkb%rr(101)/local_and_vkb%rr(1))
       write(*,fmt='(2x,"Mesh parameters"/4x,"Alpha from Hamann table: ",f12.9)') alpha_hamann
       write(*,fmt='(4x,"Default alpha:           ",f12.9)') alpha
       write(*,fmt='(4x,"Default beta:            ",f12.9)') beta
       if(abs(alpha-alpha_hamann)>1e-6_double) then ! Arbitrary?
          alpha = alpha_hamann
       end if
       do i=1,nmesh
          rr(i) = (beta/mesh_z)*exp(alpha*real(i-1,double))
          rr_squared(i) = rr(i)*rr(i)
          drdi(i) = alpha*rr(i)
          drdi_squared(i) = drdi(i)*drdi(i)
          sqrt_drdi(i) = sqrt(drdi(i))
       end do
    else if(mesh_type==siesta) then
       if(iprint>1) write(*,fmt='("# Siesta type mesh with alpha, beta: ",2f15.8)') alpha,beta
       do i=1,nmesh
          rr(i) = beta*(exp(alpha*real(i,double))-one) ! Normally in Siesta they have i-1 to get r=0
          rr_squared(i) = rr(i)*rr(i)
          drdi(i) = alpha*(rr(i) + beta)
          drdi_squared(i) = drdi(i)*drdi(i)
          sqrt_drdi(i) = sqrt(drdi(i))
          !write(15,*) i,rr(i),rr_squared(i),drdi(i)
       end do
    end if
  end subroutine make_mesh

  ! Create the regular mesh in passed array from cutoff and number of points
  subroutine make_mesh_reg(x_reg,nmesh_reg,cutoff)

    use numbers
    use datatypes
    
    implicit none

    ! Passed variables
    integer :: nmesh_reg
    real(double), dimension(nmesh_reg) :: x_reg
    real(double) :: cutoff

    ! Local variables
    integer :: i
    real(double) :: dr
    
    dr = cutoff/real(nmesh_reg-1,double)
    x_reg = zero
    do i=1,nmesh_reg
       x_reg(i) = real(i-1,double)*dr
    end do
    return
  end subroutine make_mesh_reg
  
  ! interpolate from a logarithmic to even mesh using Lagrange polynomials
  ! Generates all entries for uniform mesh, starting from zero (uses this
  ! for efficiency in finding starting points in array)
  subroutine interpolate(x_reg,y_reg,n_reg,x_irreg,y_irreg,n_irreg,max_bc)

    use numbers
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: iprint
    
    implicit none

    integer :: n_reg
    real(double), dimension(n_reg) :: x_reg, y_reg
    integer :: n_irreg
    real(double), dimension(n_irreg) :: x_irreg, y_irreg
    real(double), OPTIONAL :: max_bc ! Apply BC at cutoff

    integer :: n, i, j, m, n_start, imin, imax, ii
    integer :: nmax, n_poly_degree, npts, n_end
    real(double) :: interp, this_term
    
    ! We have function y on logarithmic mesh x_irreg (this module) as input
    ! We want to interpolate values y_reg on regular mesh x_reg (passed)

    if(x_reg(2)<x_irreg(1)) call cq_abort("Interpolation out of range: ",x_reg(2),x_irreg(1))
    ! If BC specified at rmax, then only loop to n_reg-1
    if(PRESENT(max_bc)) then
       if(iprint>5) write(*,fmt='("Interpolating with boundary condition: ",f18.10)') max_bc
       nmax = n_reg-1
       y_reg(n_reg) = max_bc
    else
       nmax = n_reg
    end if
    ! Set limits for search for the appropriate point in irregular grid
    imin = 1
    imax = n_irreg
    n_poly_degree = 2*n_poly_lag + 1
    ! Loop over output mesh points x, leaving x=0 to end
    do n=2,nmax
       ! Check for end of table
       if(x_reg(n)>x_irreg(n_irreg)) then
          if(abs(x_reg(n)-x_irreg(n_irreg))<1e-6_double) then
             y_reg(n) = y_irreg(n_irreg)
             cycle
          else
             call cq_abort("Interpolation out of range: ",x_reg(n),x_irreg(n_irreg))
          end if
       end if
       ! Find first point r which exceeds x
       do i=imin,n_irreg
          if(x_irreg(i)>x_reg(n)) exit
       end do
       ! Set starting point for interpolation, taking care to avoid an over-run
       imin = i-1
       n_start = max(1,imin - n_poly_lag)
       n_end = min(n_irreg,imin + n_poly_lag)
       interp = zero
       do j=n_start, n_end !n_start + n_poly_degree
          if(abs(y_irreg(j))<RD_ERR) cycle  ! No point in interpolating zero
          this_term = y_irreg(j)
          do m=n_start, n_end !n_start + n_poly_degree
             if(m==j) cycle
             this_term = this_term*(x_reg(n) - x_irreg(m))/(x_irreg(j) - x_irreg(m))
          end do
          interp = interp + this_term
       end do
       !call polint(x_irreg(n_start:n_end),y_irreg(n_start:n_end),n_end-n_start+1,x_reg(n),this_term)
       !write(50,*) x_reg(n),interp,this_term
       y_reg(n) = interp
    end do
    !write(50,*) '&'
    ! Now extrapolate to zero using the regular grid
    ! Idea from Don Hamann's dpnint routine - the extrapolation is
    ! reasonable for the regular (linear) grid, where for the log
    ! grid is it a huge jump and unreliable
    this_term = zero
    interp = zero
    do j= 2, 2 + n_poly_degree
       this_term = y_reg(j)
       do m= 2, 2 + n_poly_degree
          if(m==j) cycle
          this_term = this_term*(zero - x_reg(m))/(x_reg(j) - x_reg(m))
       end do
       interp = interp + this_term
    end do
    y_reg(1) = interp
    return
  end subroutine interpolate

  ! Return the logarithmic mesh point after r
  subroutine convert_r_to_i(r,i)

    use datatypes
    use numbers
    
    implicit none

    integer :: i
    real(double) :: r
    
    if(mesh_type==hamann) then
       if(r<rr(1)) then
          i=1
       else
          ! Original - found the nearest mesh point
          i = nint(log(r*mesh_z/beta)/alpha  + one)+1
          ! Now find the mesh point after r
          ! i = floor(log(r*mesh_z/beta)/alpha  + one) + 1
       end if
    else if(mesh_type==siesta) then
       i = nint(log(r/beta+one)/alpha + one)+1
    end if
    return
  end subroutine convert_r_to_i

  ! Initially just find gradient (radial)
  subroutine make_derivatives(func,dfdr,sigma) 

    use datatypes
    use numbers
    
    implicit none

    ! Passed variables
    real(double), dimension(nmesh) :: func, dfdr
    real(double), dimension(nmesh), optional :: sigma

    ! Local variables
    integer :: i
    real(double) :: c1, c2, c3, c4, c5 ! FD coefficients for mid-points
    real(double) :: d1, d2, d3, d4, d5 ! FD coefficients for end points
    real(double) :: dfdn

    c1 = one/twelve
    c2 = -eight/twelve
    c3 = zero
    c4 = eight/twelve
    c5 = -one/twelve

    d1 = -five*five/twelve
    d2 = four
    d3 = -three
    d4 = four/three
    d5 = -one/four

    ! Find df/dn then df/dr
    ! i=1 and i=2 done separately
    i=1
    dfdn = d1*func(i) + d2*func(i+1) + d3*func(i+2) + d4*func(i+3) + d5*func(i+4)
    dfdr(i) = dfdn/(drdi(i))
    i=2
    dfdn = d1*func(i) + d2*func(i+1) + d3*func(i+2) + d4*func(i+3) + d5*func(i+4)
    dfdr(i) = dfdn/(drdi(i))
    ! Mid points
    do i=3,nmesh-2
       ! Omit redundant func(i) as it's scaled by c3=0
       dfdn = c1*func(i-2) + c2*func(i-1) + c4*func(i+1) + c5*func(i+2)
       dfdr(i) = dfdn/(drdi(i))
    end do
    ! For the last two points use d1->d5 with signs reversed
    i=nmesh-1
    dfdn = -d1*func(i) -d2*func(i-1) -d3*func(i-2) -d4*func(i-3) -d5*func(i-4)
    dfdr(i) = dfdn/(drdi(i))
    i=nmesh
    dfdn = -d1*func(i) -d2*func(i-1) -d3*func(i-2) -d4*func(i-3) -d5*func(i-4)
    dfdr(i) = dfdn/(drdi(i))
    if(present(sigma)) then
       do i=1,nmesh
          sigma(i) = dfdr(i)*dfdr(i)
       end do
    end if
  end subroutine make_derivatives
end module mesh
