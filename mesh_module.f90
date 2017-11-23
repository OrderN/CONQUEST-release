module mesh

  use datatypes
  
  implicit none
  
  integer :: nmesh
  integer :: nmesh_reg

  ! Default values
  real(double) :: alpha = log(1.012_double)
  real(double) :: beta = 0.0005_double
  real(double) :: mesh_z
  real(double) :: delta_r_reg
  real(double), allocatable, dimension(:) :: rr, rr_squared, drdi, sqrt_rr, drdi_squared, sqrt_drdi, rmesh_reg
  integer, parameter :: hamann = 1
  integer, parameter :: siesta = 2
  integer :: mesh_type
  integer :: n_poly_lag = 5 ! User adjustable ?
  integer :: n_poly = 5 ! Should this be user adjustable ? Actually gives 2*n_poly+1

  
  
contains

  subroutine make_mesh(species)

    use datatypes
    use numbers
    use pseudo_tm_info, ONLY: pseudo
    use global_module, ONLY: iprint
    
    implicit none

    ! Passed variables
    integer :: species
    
    ! Local variables
    integer :: i
    
    if(.NOT.allocated(rr)) then
       allocate(rr(nmesh),rr_squared(nmesh),drdi(nmesh),sqrt_rr(nmesh),drdi_squared(nmesh),sqrt_drdi(nmesh))
       rr = zero
       rr_squared = zero
       drdi = zero
       sqrt_rr = zero
       drdi_squared = zero
       sqrt_drdi = zero
    end if
    !write(*,*) '# Mesh: ',mesh_type,alpha,beta
    mesh_z = pseudo(species)%z
    if(mesh_type==hamann) then
       if(iprint>3) write(*,fmt='("# Hamann type mesh with alpha, beta, Z: ",3f15.8)') alpha,beta, mesh_z
       do i=1,nmesh
          rr(i) = (beta/mesh_z)*exp(alpha*real(i-1,double))
          rr_squared(i) = rr(i)*rr(i)
          sqrt_rr(i) = sqrt(rr(i))
          drdi(i) = alpha*rr(i)
          drdi_squared(i) = drdi(i)*drdi(i)
          sqrt_drdi(i) = sqrt(drdi(i))
       end do
    else if(mesh_type==siesta) then
       if(iprint>1) write(*,fmt='("# Siesta type mesh with alpha, beta: ",2f15.8)') alpha,beta
       do i=1,nmesh
          rr(i) = beta*(exp(alpha*real(i,double))-one) ! Normally in Siesta they have i-1 to get r=0
          rr_squared(i) = rr(i)*rr(i)
          sqrt_rr(i) = sqrt(rr(i))
          drdi(i) = alpha*(rr(i) + beta)
          drdi_squared(i) = drdi(i)*drdi(i)
          sqrt_drdi(i) = sqrt(drdi(i))
          !write(15,*) i,rr(i),rr_squared(i),drdi(i)
       end do
    end if
  end subroutine make_mesh

  ! Create the regular mesh
  subroutine make_mesh_reg(cutoff,n_points,reg_cutoff)

    use numbers
    use datatypes
    
    implicit none

    integer, OPTIONAL :: n_points
    real(double) :: cutoff
    real(double), OPTIONAL :: reg_cutoff

    integer :: i, n_use
    
    if(allocated(rmesh_reg)) deallocate(rmesh_reg)
    if(PRESENT(n_points)) then
       n_use = n_points
       if(PRESENT(reg_cutoff)) then
          reg_cutoff = delta_r_reg*real(n_points-1,double)
       !else
       !   write(*,*) '# WARNING: using regular mesh without changing final cutoff'
       end if
    else
       n_use = cutoff/delta_r_reg + 1 !nmesh_reg
    end if
    allocate(rmesh_reg(n_use))
    rmesh_reg = zero
    !delta_r_reg = cutoff/real(n_use-1)
    do i=1,n_use
       rmesh_reg(i) = real(i-1,double)*delta_r_reg
    end do
    return
  end subroutine make_mesh_reg

  ! New polynomial interpolation following NR
  subroutine new_interpolate(x_reg,y_reg,n_reg,x_irreg,y_irreg,n_irreg,max_bc)

    use numbers
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: iprint
    
    implicit none

    integer :: n_reg
    real(double), dimension(n_reg) :: x_reg, y_reg
    integer :: n_irreg
    real(double), dimension(n_irreg) :: x_irreg, y_irreg
    real(double), OPTIONAL :: max_bc ! Apply BC at cutoff

    integer :: nr, i
    integer :: nmax, nmin, npts
    real(double) :: interp, this_term, r

    ! Loop over regular mesh
    do i=1,n_reg-1
       r = x_reg(i)
       ! Convert r into grid point
       call convert_r_to_i(r,nr)
       if(nr<1) write(*,*) 'Error: ',r,nr
       nmin = max(1,nr - n_poly)
       nmax = min(n_irreg,nr + n_poly)
       npts = nmax - nmin + 1
       call polint(x_irreg(nmin:nmax),y_irreg(nmin:nmax),npts,r,interp)
       y_reg(i) = interp
    end do
    if(PRESENT(max_bc)) then
       y_reg(n_reg) = max_bc
    else
       y_reg(n_reg) = zero
    end if
    return
  end subroutine new_interpolate


  ! Polynomial interpolation. 
  subroutine polint(xa,ya,n,x,y)

    use numbers
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer          :: n
    real(double) :: xa(n),ya(n), x, y, dy

    ! Local variables
    integer          :: i, m, ns
    real(double), dimension(n) :: c, d
    real(double) :: den, diff, dift, ho, hp, w

    ns=1
    diff=abs(x-xa(1))
    do i=1,n 
       dift=abs(x-xa(i))
       if (dift<diff) then
          ns=i
          diff=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    end do ! i
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if (den.eq.zero) call cq_abort('Error in polint: two xas are equal')
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       end do ! i
       if ((2*ns)<(n-m)) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    end do ! m
    return
  end subroutine polint
  
  ! interpolate from a logarithmic to even mesh
  ! Uses Lagrange polynomials for now
  subroutine interpolate(x,y_reg,n_reg,y,n_irreg,max_bc)

    use numbers
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: iprint
    
    implicit none

    integer :: n_reg
    real(double), dimension(n_reg) :: x, y_reg
    integer :: n_irreg
    real(double), dimension(n_irreg) :: y
    real(double), OPTIONAL :: max_bc ! Apply BC at cutoff

    integer :: n, i, j, m, n_start, imin, imax, ii
    integer :: nmax
    real(double) :: interp, this_term
    
    ! We have function y on mesh r as input (r will be the logarithmic mesh defined in this module)
    ! We want to output y_reg on mesh x (likely to be regular)

    if(x(2)<rr(1)) call cq_abort("Interpolation out of range: ",x(2),rr(1))
    ! If BC specified at rmax, then only loop to n_reg-1
    if(PRESENT(max_bc)) then
       if(iprint>5) write(*,fmt='("Interpolating with boundary condition: ",f18.10)') max_bc
       nmax = n_reg-1
       y_reg(n_reg) = max_bc
    else
       nmax = n_reg
    end if
    ! Loop over output mesh points x
    ! NB n=1 corresponding to x=0 is a special case
    imin = 1
    imax = n_irreg
    do n=2,nmax
       if(x(n)>rr(n_irreg)) then
          if(abs(x(n)-rr(n_irreg))<1e-6_double) then
             y_reg(n) = y(n_irreg)
             cycle
          else
             call cq_abort("Interpolation out of range: ",x(n),rr(n_irreg))
          end if
       end if
       ! Find points r which bracket x
       ! Start point
       do i=imin,n_irreg
          if(rr(i)>x(n)) exit
       end do
       imin = i-1
       imax = i
       !%%!imin = 1
       !%%!imax = n_irreg
       !%%!do i=1,n_irreg
       !%%!   ii=(imin+imax)/2
       !%%!   if(x(n)>rr(ii)) then
       !%%!      imin = ii
       !%%!   else
       !%%!      imax = ii
       !%%!   end if
       !%%!   if(imax-imin==1) exit
       !%%!end do
       ! Lagrangian interpolation
       if(mod(n_poly_lag,2)==1) then
          n_start = imin - n_poly_lag/2
       else if(x(n)-rr(imin)<rr(imax)-x(n)) then
          n_start = imin - n_poly_lag/2
       else
          n_start = imax - n_poly_lag/2
       end if
       ! Ensure that the start point doesn't cause an over-run
       n_start = max(1,min(n_start,n_irreg - n_poly_lag))
       interp = zero
       do j=n_start,n_start + n_poly_lag
          if(abs(y(j))<RD_ERR) cycle  ! No point in interpolating zero
          this_term = y(j)
          do m=n_start,n_start + n_poly_lag
             if(m==j) cycle
             this_term = this_term*(x(n) - rr(m))/(rr(j) - rr(m))
          end do
          interp = interp + this_term
       end do
       y_reg(n) = interp
    end do
    ! Now interpolate to zero using the regular grid
    this_term = zero
    interp = zero
    do j= 2, 2 + n_poly_lag
       this_term = y_reg(j)
       do m= 2, 2 + n_poly_lag
          if(m==j) cycle
          this_term = this_term*(zero - x(m))/(x(j) - x(m))
       end do
       interp = interp + this_term
    end do
    y_reg(1) = interp
    !if(PRESENT(max_bc)) y_reg(n_reg) = max_bc
    return
  end subroutine interpolate

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
          i = nint(log(r*mesh_z/beta)/alpha  + one)+1
       end if
    else if(mesh_type==siesta) then
       i = nint(log(r/beta+one)/alpha + one)+1
    end if
    return
  end subroutine convert_r_to_i
end module mesh
