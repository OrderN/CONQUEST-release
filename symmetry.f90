! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module symmetry
! ------------------------------------------------------------------------------
! Code area 11: Basis operations
! ------------------------------------------------------------------------------

!!****h* Conquest/symmetry *
!!  NAME
!!   symmetry
!!  PURPOSE
!!   Sole purpose is to contain three routines that perform
!!   symmetry-related tasks in the calculation of blip
!!   coefficients for support functions initially represented
!!   in terms of PAO's.
!!  USES
!!
!!  AUTHOR
!!   Mike Gillan
!!  CREATION DATE
!!   22/6/02
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
module symmetry

  implicit none

  save

!!***

contains

! -----------------------------------------------------------
! Subroutine star_coeffs
! -----------------------------------------------------------

!!****f* symmetry/star_coeffs *
!!
!!  NAME 
!!   star_coeffs
!!  USAGE
!! 
!!  PURPOSE
!!   For a symmetry type specified by sym_type
!! (currently, can be one of the representations a1g, t1u, t2g, eg
!!   of the cubic point group), makes list of stars of blip-grid points,
!!   and calculates components of symmetry-adapted basis vector 
!! (reciprocal of sqrt of multiplicity) for each.
!! 
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
  subroutine star_coeffs(sym_type,nbh,region_bound_single,&
       &region_bound_double,n_star_in_sphere, &
       n_cube2sphere,coeff,xsite,ysite,zsite)
    use datatypes
    use GenComms, ONLY: cq_abort
    use numbers
    implicit none

    character(len=*), intent(in) :: sym_type
    integer, intent(in) :: nbh
    integer, intent(in), dimension(0:) :: region_bound_single
    integer, intent(in), dimension(0:,0:) :: region_bound_double
    integer, intent(out) :: n_star_in_sphere
    integer, intent(out), dimension(:) :: n_cube2sphere
    integer, intent(out), dimension(:) :: xsite, ysite, zsite
    real(double), intent(out), dimension(:) :: coeff
    integer :: i, ia, ib, ic, ix, iy, iz, i1, is, m, nbh_min, n_star_in_cube
    real(double) :: c, f1, recip_sqrt_two

    ! check validity of input parameters
    if((sym_type /= "a1g") .and. (sym_type /= "t1u") .and. &
         (sym_type /= "t2g") .and. (sym_type /= "eg")) then
       call cq_abort('star_coeffs: unrecognised symmetry type')
    end if
    if(sym_type == "a1g") then
       nbh_min = 0
    else
       nbh_min = 1
    end if
    if(nbh < nbh_min) then
       call cq_abort('star_coeffs: nbh out of range',nbh)
    end if

    ! constants
    recip_sqrt_two = one/sqrt(2.0_double)

    ! for each different symmetry type, make star of blip-grid sites

    select case(sym_type)

       ! ===============================================================
       ! a1g symmetry (fully symmetric, angular momentum l = 0)

    case("a1g")

       n_star_in_cube = ((nbh+1)*(nbh+2)*(nbh+3))/6
       do i = 1, n_star_in_cube
          n_cube2sphere(i) = 0
       end do

       is = 0

       do ix = 0, nbh
          do iy = 0, min(ix,region_bound_single(ix))
             do iz = 0, min(iy,region_bound_double(ix,iy))
                i = (ix*(ix+1)*(ix+2))/6 + (iy*(iy+1))/2 + iz + 1
                is = is + 1
                if(ix == 0) then
                   coeff(is) = 1.0_double
                else if(iy == 0) then
                   coeff(is) = 1.0_double/sqrt(6.0_double)
                else if(iz == 0) then
                   if(ix == iy) then
                      coeff(is) = 1.0_double/sqrt(12.0_double)
                   else
                      coeff(is) = 1.0_double/sqrt(24.0_double)
                   end if
                else
                   if((ix == iy) .and. (iy == iz)) then
                      coeff(is) = 1.0_double/sqrt(8.0_double)
                   else if((ix == iy) .or. (iy == iz)) then
                      coeff(is) = 1.0_double/sqrt(24.0_double)
                   else
                      coeff(is) = 1.0_double/sqrt(48.0_double)
                   end if
                end if
                n_cube2sphere(i) = is
             end do
          end do
       end do

       n_star_in_sphere = is

       ! store vector position of first site in each star 
       do ix = 0, nbh
          do iy = 0, ix
             do iz = 0, iy
                i = (ix*(ix+1)*(ix+2))/6 + (iy*(iy+1))/2 + iz + 1
                is = n_cube2sphere(i)
                if(is /= 0) then
                   xsite(is) = ix
                   ysite(is) = iy
                   zsite(is) = iz
                end if
             end do
          end do
       end do

       ! ===============================================================
       ! t1u symmetry (angular momentum l = 1, z function)

    case("t1u")

       n_star_in_cube = (nbh*(nbh+1)*(nbh+2))/2
       do i = 1, n_star_in_cube
          n_cube2sphere(i) = 0
       end do

       is = 0

       ! loop over x and y values
       do ix = 0, nbh
          ia = (ix*(ix+1))/2
          do iy = 0, min(ix,region_bound_single(ix))
             ib = nbh*(ia+iy)
             if((ix == 0) .and. (iy == 0)) then
                m = 2
                c = 1.0_double/sqrt(2.0_double)
             else if(iy == 0) then
                m = 8
                c = 1.0_double/sqrt(8.0_double)
             else if(ix == iy) then
                m = 8
                c = 1.0_double/sqrt(8.0_double)
             else
                m = 16
                c = 1.0_double/sqrt(16.0_double)
             end if
             if(region_bound_double(ix,iy) > 0) then
                do iz = 1, region_bound_double(ix,iy)
                   i = ib + iz
                   is = is + 1
                   coeff(is) = c
                   n_cube2sphere(i) = is
                end do
             end if
          end do
       end do

       n_star_in_sphere = is

       ! store vector position of first site in each star
       do ix = 0, nbh
          do iy = 0, ix
             do iz = 1, nbh
                i = nbh*((ix*(ix+1))/2+iy) + iz
                is = n_cube2sphere(i)
                if(is /= 0) then
                   xsite(is) = ix
                   ysite(is) = iy
                   zsite(is) = iz
                end if
             end do
          end do
       end do

       ! ===============================================================
       ! t2g symmetry (angular momentum l = 2, xy function)

    case("t2g")

       n_star_in_cube = (nbh*(nbh+1)*(nbh+1))/2
       do i = 1, n_star_in_cube
          n_cube2sphere(i) = 0
       end do

       is = 0

       ! loop over x and y values
       do ix = 1, nbh
          if(region_bound_single(ix) > 0) then
             ia = ((ix-1)*ix)/2
             do iy = 1, min(ix,region_bound_single(ix))
                ib = (nbh+1)*(ia+iy-1) + 1
                if(ix == iy) then
                   m = 4
                   c = 1.0_double/sqrt(4.0_double)
                else
                   m = 8
                   c = 1.0_double/sqrt(8.0_double)
                end if
                is = is + 1
                coeff(is) = c
                n_cube2sphere(ib) = is
                if(region_bound_double(ix,iy) > 0) then
                   do iz = 1, region_bound_double(ix,iy)
                      is = is + 1
                      coeff(is) = c*recip_sqrt_two
                      n_cube2sphere(ib+iz) = is
                   end do
                end if
             end do
          end if
       end do

       n_star_in_sphere = is

       ! store vector position of first site in each star
       do ix = 1, nbh
          do iy = 1, ix
             do iz = 0, nbh 
                i = (nbh+1)*(((ix-1)*ix)/2 + iy - 1) + iz + 1
                is = n_cube2sphere(i)
                if(is /= 0) then
                   xsite(is) = ix
                   ysite(is) = iy
                   zsite(is) = iz
                end if
             end do
          end do
       end do

       ! ===============================================================
       ! eg symmetry (angular momentum l = 2, x^2 - y^2 function)

    case("eg")

       n_star_in_cube = (nbh*(nbh+1)*(nbh+1))/2
       do i = 1, n_star_in_cube
          n_cube2sphere(i) = 0
       end do

       is = 0

       ! loop over x and y values
       do ix = 1, nbh
          ia = ((ix-1)*ix)/2
          do iy = 0, min(ix-1,region_bound_single(ix))
             ib = (nbh+1)*(ia+iy) + 1
             if(iy == 0) then
                m = 4
                c = 1.0_double/sqrt(4.0_double)
             else
                m = 8
                c = 1.0_double/sqrt(8.0_double)
             end if
             is = is + 1
             coeff(is) = c
             n_cube2sphere(ib) = is
             if(region_bound_double(ix,iy) > 0) then
                do iz = 1, region_bound_double(ix,iy)
                   is = is + 1
                   coeff(is) = c*recip_sqrt_two
                   n_cube2sphere(ib+iz) = is
                end do
             end if
          end do
       end do

       n_star_in_sphere = is

       ! store vector position of site in each star
       do ix = 1, nbh
          do iy = 0, ix-1
             do iz = 0, nbh
                i = (nbh+1)*(((ix-1)*ix)/2 + iy) + iz + 1
                is = n_cube2sphere(i)
                if(is /= 0) then
                   xsite(is) = ix
                   ysite(is) = iy
                   zsite(is) = iz
                end if
             end do
          end do
       end do

    end select

    return

  end subroutine star_coeffs
! -------------------------------------------------------------
!!*** 

! -----------------------------------------------------------
! Subroutine symmat
! -----------------------------------------------------------

!!****f* symmetry/symmat *
!!
!!  NAME 
!!   symmat
!!  USAGE
!! 
!!  PURPOSE
!!   Constructs block of blip overlap matrix having symmetry
!!   specified by sym_type.
!! 
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
!!
!!  SOURCE
  subroutine symmat(sym_type,nbh,n_star_in_sphere,n_cube2sphere, &
       xsite,ysite,zsite,coeff,amat)
    use datatypes
    use numbers
    implicit none

    real(double), parameter :: p0 = 151.0_double/140.0_double, &
         p1 = 1191.0_double/2240.0_double, p2 = 3.0_double/56.0_double, &
         p3 = 1.0_double/2240.0_double
    integer, intent(in) :: nbh, n_star_in_sphere
    integer, intent(in), dimension(:) :: n_cube2sphere
    integer, intent(in), dimension(:) :: xsite, ysite, zsite
    character(len=*), intent(in) :: sym_type
    real(double), intent(in), dimension(:) :: coeff
    real(double), intent(out), dimension(:,:) :: amat
    integer :: is, i1, i2, ix1, iy1, iz1, ix2, iy2, iz2, &
         ix2min, ix2max, iy2min, iy2max, iz2min,iz2max, &
         abs_ix2, abs_iy2, abs_iz2, min_xy, max_xy, min_xyz, mid_xyz, max_xyz, &
         sign_x, sign_y, sign_z, sign_xy
    real(double) :: c, xfactor, yfactor, zfactor
    real(double), dimension(0:3) :: overlap

    ! constants
    overlap(0) = p0
    overlap(1) = p1
    overlap(2) = p2
    overlap(3) = p3

    ! initialise accumulator to zero
    do i1 = 1, n_star_in_sphere
       do i2 = 1, n_star_in_sphere
          amat(i1,i2) = zero
       end do
    end do

    ! loop over stars i1
    do i1 = 1, n_star_in_sphere
       ix1 = xsite(i1)
       iy1 = ysite(i1)
       iz1 = zsite(i1)
       ix2max = min(ix1+3,nbh)
       ix2min = max(ix1-3,-nbh)
       iy2max = min(iy1+3,nbh)
       iy2min = max(iy1-3,-nbh)
       iz2max = min(iz1+3,nbh)
       iz2min = max(iz1-3,-nbh)

       ! loop over offsets from first site in star i1

       do ix2 = ix2min, ix2max
          abs_ix2 = abs(ix2)
          xfactor = overlap(abs(ix2-ix1))
          if(ix2 /= 0) sign_x = sign(1,ix2)

          do iy2 = iy2min, iy2max
             abs_iy2 = abs(iy2)
             yfactor = xfactor*overlap(abs(iy2-iy1))
             if(iy2 /= 0) sign_y = sign(1,iy2)
             sign_xy = 1
             if(abs_iy2 > abs_ix2) sign_xy = -1
             min_xy = min(abs_ix2,abs_iy2)
             max_xy = max(abs_ix2,abs_iy2)

             do iz2 = iz2min, iz2max
                abs_iz2 = abs(iz2)
                zfactor = yfactor*overlap(abs(iz2-iz1))
                if(iz2 /= 0) sign_z = sign(1,iz2)

                select case(sym_type)

                case("a1g")

                   max_xyz = max(max_xy,abs_iz2)
                   min_xyz = min(min_xy,abs_iz2)
                   mid_xyz = abs_ix2 + abs_iy2 + abs_iz2 - max_xyz - min_xyz
                   i2 = (max_xyz*(max_xyz+1)*(max_xyz+2))/6 + &
                        (mid_xyz*(mid_xyz+1))/2 + min_xyz + 1
                   is = n_cube2sphere(i2)
                   if(is /= 0) amat(i1,is) = amat(i1,is) + zfactor*coeff(is)

                case("t1u")

                   if(iz2 /= 0) then
                      i2 = nbh*((max_xy*(max_xy+1))/2 + min_xy) + abs_iz2
                      is = n_cube2sphere(i2)
                      if(is /= 0) amat(i1,is) = amat(i1,is) + &
                           sign_z*zfactor*coeff(is)
                   end if

                case("t2g")

                   if(min_xy /= 0) then
                      i2 = (nbh+1)*((max_xy*(max_xy-1))/2 + min_xy - 1) + abs_iz2 + 1
                      is = n_cube2sphere(i2)
                      if(is /= 0) amat(i1,is) = amat(i1,is) + &
                           sign_x*sign_y*zfactor*coeff(is)
                   end if

                case("eg")

                   if(min_xy /= max_xy) then
                      i2 = (nbh+1)*((max_xy*(max_xy-1))/2 + min_xy) + abs_iz2 + 1
                      is = n_cube2sphere(i2)
                      if(is /= 0) amat(i1,is) = amat(i1,is) + &
                           sign_xy*zfactor*coeff(is)
                   end if

                end select

             end do
          end do
       end do

    end do

    do i1 = 1, n_star_in_sphere
       c = one/coeff(i1)
       do i2 = 1, n_star_in_sphere
          amat(i1,i2) = c*amat(i1,i2)
       end do
    end do

    return

  end subroutine symmat
!!*** 

! -----------------------------------------------------------
! Subroutine star2sphere
! -----------------------------------------------------------

!!****f* symmetry/star2sphere *
!!
!!  NAME 
!!   star2sphere
!!  USAGE
!! 
!!  PURPOSE
!!   Takes blip coefficients scal_prod_sym() for stars of sites
!!   associated with a specified symmetry sym_type, and uses
!!   them to make explicit values of blip coefficients for
!!   all sites of blips wholly contained within the support region,
!!   the latter explicit values being put into the array blip_contrib().
!! 
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
  subroutine star2sphere(sym_type,n_angcomp,nbh,n_blips_region,&
       &blip_location,n_cube2sphere,scal_prod_sym,blip_contrib)
    use datatypes
    use numbers
    implicit none

    character(len=*), intent(in) :: sym_type
    integer, intent(in) :: n_angcomp, nbh, n_blips_region
    integer, intent(in), dimension(:,:) :: blip_location
    integer, intent(in), dimension(:) :: n_cube2sphere
    real(double), intent(in), dimension(:) :: scal_prod_sym
    real(double), intent(out), dimension(:) :: blip_contrib
    integer :: abs_nx, abs_ny, abs_nz, i, k, max_xy, min_xy, &
         max_xyz, min_xyz, mid_xyz, n, nx, ny, nz, ndum, &
         sign_x, sign_y, sign_z, sign_xy
    real(double) :: over_sign

    ! --- loop over blips wholly contained within support region ------
    do n = 1 , n_blips_region
       nx = blip_location(1,n)
       ny = blip_location(2,n)
       nz = blip_location(3,n)

       ! select requested symmetry ---------------------------------------
       select case(sym_type)

          ! symmetry: angular momentum l = 0, full spherical symmetry -------
       case("a1g")

          abs_nx = abs(nx) ; abs_ny = abs(ny) ; abs_nz = abs(nz)
          max_xyz = max(abs_nx,abs_ny,abs_nz)
          min_xyz = min(abs_nx,abs_ny,abs_nz)
          mid_xyz = abs_nx + abs_ny + abs_nz - max_xyz - min_xyz
          i = (max_xyz*(max_xyz+1)*(max_xyz+2))/6 + & 
               (mid_xyz*(mid_xyz+1))/2 + min_xyz + 1
          blip_contrib(n) = scal_prod_sym(n_cube2sphere(i))

          ! symmetry: angular momentum l = 1 --------------------------------
          ! n_angcomp = 2, 3, 4, corresponding to x, y, z respectively
       case("t1u")

          if(n_angcomp == 2) then
             ndum = nz ; nz = nx ; nx = ny ; ny = ndum
          else if(n_angcomp == 3) then
             ndum = nz ; nz = ny ; ny = nx ; nx = ndum
          end if
          if(nz == 0) then
             blip_contrib(n) = zero
          else
             sign_z = sign(1,nz)
             abs_nx = abs(nx) ; abs_ny = abs(ny) ; abs_nz = abs(nz)
             max_xy = max(abs_nx,abs_ny)
             min_xy = min(abs_nx,abs_ny)
             i = nbh*((max_xy*(max_xy+1))/2 + min_xy) + abs_nz
             blip_contrib(n) = sign_z*scal_prod_sym(n_cube2sphere(i))
          end if

          ! symmetry: angular momentum l = 2 --------------------------------
          ! n_angcomp = 5, 6, 7, corresponding to yz, zx, xy respectively
       case("t2g")

          if(n_angcomp == 5) then
             ndum = nz ; nz = nx ; nx = ny ; ny = ndum
          else if(n_angcomp == 6) then
             ndum = nz ; nz = ny ; ny = nx ; nx = ndum
          end if
          sign_x = sign(1,nx)
          sign_y = sign(1,ny)
          abs_nx = abs(nx) ; abs_ny = abs(ny) ; abs_nz = abs(nz)
          if(min(abs_nx,abs_ny) == 0) then
             blip_contrib(n) = zero
          else
             max_xy = max(abs_nx,abs_ny)
             min_xy = min(abs_nx,abs_ny)
             i = (nbh+1)*((max_xy*(max_xy-1))/2 + min_xy - 1) + abs_nz + 1
             blip_contrib(n) = sign_x*sign_y*scal_prod_sym(n_cube2sphere(i))
          end if

          ! symmetry: angular momentum l = 2 -------------------------------
          ! n_angcomp = 8, 9, corresponding to x^2 - y^2, 3z^2 - r^2, respectively
       case("eg")

          if(n_angcomp == 8) then
             abs_nx = abs(nx) ; abs_ny = abs(ny) ; abs_nz = abs(nz)
             if(abs_nx == abs_ny) then
                blip_contrib(n) = zero
             else 
                sign_xy = 1
                if(abs_ny > abs_nx) sign_xy = -1
                max_xy = max(abs_nx,abs_ny)
                min_xy = min(abs_nx,abs_ny)
                i = (nbh+1)*((max_xy*(max_xy-1))/2 + min_xy) + abs_nz + 1
                blip_contrib(n) = sign_xy*scal_prod_sym(n_cube2sphere(i))
             end if
          else
             blip_contrib(n) = zero
             over_sign = -1.0_double/sqrt(3.0_double)
             do k = 1, 2
                ndum = nz ; nz = nx ; nx = ny ; ny = ndum
                abs_nx = abs(nx) ; abs_ny = abs(ny) ; abs_nz = abs(nz)
                if(abs_nx /= abs_ny) then
                   sign_xy = 1
                   if(abs_ny > abs_nx) sign_xy = -1
                   max_xy = max(abs_nx,abs_ny)
                   min_xy = min(abs_nx,abs_ny)
                   i = (nbh+1)*((max_xy*(max_xy-1))/2 + min_xy) + abs_nz + 1
                   blip_contrib(n) = blip_contrib(n) + &
                        &over_sign*sign_xy*scal_prod_sym(n_cube2sphere(i))
                end if
                over_sign = -over_sign
             end do
          end if

       end select
       ! end select of requested symmetry ---------------------------------

    end do
    ! --- end loop over blips wholly contained within support region ---

    return

  end subroutine star2sphere
!!*** 

end module symmetry
