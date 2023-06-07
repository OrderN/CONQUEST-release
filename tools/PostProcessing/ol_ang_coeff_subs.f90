! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module angular_coeff_routines
! ------------------------------------------------------------------------------
! Code area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/angular_coeff_routines
!!  NAME
!!   angular_coeff_routines
!!  PURPOSE
!!   Holds routines to calculate angular coefficients relevant to 
!!   evaluation of the radial tables.
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2005/07/11 dave
!!    changed inout specification to in for ifort compilation
!!   2007/01/05 17:30 dave
!!    included prefac and fact arrays
!!   2008/02/06 08:27 dave
!!    Changed for output to file not stdout
!!   2008/06/10 ast
!!    Added timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2018/01/22 12:42 JST dave
!!    Changes to make prefactor arrays allocatable (for NA projectors)
!!   2019/10/30 16:53 dave
!!    Remove calculation of (-1)**m throughout (replace with test for m even/odd
!!    and scaling factor of 1 or -1; also added numbers for module use.  Also
!!    introduced a universal epsilon parameter for small angles
!!  SOURCE
!!

module angular_coeff_routines

  use datatypes
  use numbers,                only: zero, quarter, half, one, three_halves, two, &
       three, pi, eight, very_small

  implicit none

  ! Parameter to avoid small angle errors throughout code
  real(double), parameter :: epsilon = 1.0e-4_double 
  real(double), allocatable, dimension(:,:) :: prefac
  real(double), allocatable, dimension(:,:) :: prefac_real
  ! Store factorial (n!) and double factorial (n!!=n*(n-2)*...)
  real(double), allocatable, dimension(:) :: fact, doublefact
!!***

contains

  
!!****f* angular_coeff_routines/convert_basis *
!!
!!  NAME 
!!   convert_basis
!!  USAGE
!!   convert_basis(x,y,z,r,thet,phi)
!!  PURPOSE
!!   Converts a displacement given in Cartesian coords
!!   to spherical polar coords.
!!  INPUTS
!!   x, y, z Cartesian components of displacement vector 
!!   r, thet, phi - spherical polar components of displacement vector
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   13:16, 21/10/2005 drb 
!!    Removed potential divide by zero (thanks to Andy Gormanly)
!!  SOURCE
!!

  subroutine convert_basis(x,y,z,r,thet,phi)
    use datatypes

    implicit none
    !routine to convert Cartesian displacements into displacements
    !in spherical polar coordinates
    real(double), intent(in) :: x,y,z
    real(double), intent(inout) :: r,thet,phi
    real(double) :: x2,y2,z2,mod_r,mod_r_plane,num

    x2 = x*x
    y2 = y*y
    z2 = z*z
    mod_r = sqrt(x2+y2+z2)
    mod_r_plane = sqrt(x2+y2)
    if(mod_r<very_small) then
       r = zero
       thet = zero
       phi = zero
       return
    end if

    if(abs(z)<very_small) then
       thet = half*pi
    else
       thet = acos(z/mod_r) !need to make sure this ratio doesn't disappear as well..
    endif

    if(abs(x2)<very_small.and.abs(y2)<very_small) then
       phi = zero
    else
       if(x<zero) then
          if(y<zero) then
             !3rd quadrant
             phi = pi+acos(-x/mod_r_plane)
          else 
             !2nd quadrant
             phi = (half*pi)+acos(y/mod_r_plane)
          endif
       else
          if(y<zero) then
             !4th quadrant
             phi = (three_halves*pi)+acos(-y/mod_r_plane)
          else
             !1st quadrant
             phi = acos(x/mod_r_plane)
          endif
       endif
    endif
   
    r = mod_r

  end subroutine convert_basis
  !!***

!!****f* angular_coeff_routines/evaluate_pao *
!!
!!  NAME 
!!   evaluate_pao
!!  USAGE
!!   evaluate_pao(sp,l,nz,m,x,y,z,pao_val)
!!  PURPOSE
!!   Returns value of specified PAO at given point in space
!!  INPUTS
!!   x, y, z Cartesian components of displacement vector 
!!   sp,l,nz,m - species no., l ang momentum value, zeta value
!!   and azimuthal quantum number of specified pao
!!  USES
!!   datatypes, pao_format, pao_array_utility
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   25/09/03
!!  MODIFICATION HISTORY
!!   2019/08/16 15:57 dave
!!    Replace spline_ol_intval_new2
!!  SOURCE
!!

  subroutine evaluate_pao(sp,l,nz,m,x,y,z,pao_val)

    use datatypes
    use numbers
    use pao_format, ONLY : pao
    use GenComms, ONLY: cq_abort

    implicit none

    integer, intent(in) :: sp,l,nz,m
    real(double), intent(in) :: x,y,z
    real(double), intent(out) :: pao_val
    real(double) :: r,theta,phi,del_r,y_val,val
    integer :: npts, j
    real(double) :: a, b, c, d, r1, r2, r3, r4, rr
    
    !convert Cartesians into spherical polars
    call convert_basis(x,y,z,r,theta,phi)
    !interpolate for required value of radial function
    npts = pao(sp)%angmom(l)%zeta(nz)%length
    del_r = (pao(sp)%angmom(l)%zeta(nz)%cutoff/&
         &(pao(sp)%angmom(l)%zeta(nz)%length-1))
    j = floor(r/del_r) + 1
    pao_val = zero
    if(j+1<=npts) then
       rr = real(j,double)*del_r
       a = (rr - r)/del_r
       b = one - a
       c = a * ( a * a - one ) * del_r * del_r / six
       d = b * ( b * b - one ) * del_r * del_r / six
       r1 = pao(sp)%angmom(l)%zeta(nz)%table(j)
       r2 = pao(sp)%angmom(l)%zeta(nz)%table(j+1)
       r3 = pao(sp)%angmom(l)%zeta(nz)%table2(j)
       r4 = pao(sp)%angmom(l)%zeta(nz)%table2(j+1)
       pao_val = a*r1 + b*r2 + c*r3 + d*r4
    end if
    !now multiply by value of spherical harmonic
    y_val = re_sph_hrmnc(l,m,theta,phi)
    
    pao_val = pao_val*y_val
    ! Scale by r**l to remove Siesta normalisation
    if(l==1) then ! p
       pao_val = pao_val*r
    else if(l==2) then ! d
       pao_val = pao_val*r*r
    else if(l==3) then ! f
       pao_val = pao_val*r*r*r
    else if(l>3) then
       call cq_abort("Angular momentum error: l>3 ",l)
    end if
  end subroutine evaluate_pao
  !!***

!!****f* angular_coeff_routines/pp_elem_derivative *
!!
!!  NAME 
!!   pp_elem_derivative
!!  USAGE
!!   pp_elem(i_vector,f_r,df_r,l,x_n,y_n,z_n,r,l,d_val)
!!  PURPOSE
!!   Evaluates value of gradient of f(r)*Re(Y_lm) along
!!   a specified Cartesian unit vector
!!  INPUTS
!!   x_n, y_n, z_n, r : direction cosines, value of radial distance  
!!   l -  ang momentum value
!!   f_r,df_r : value of radial function and first derivative
!!   i_vector : direction of unit vector (1:x, 2:y, 3:z) 
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   30/09/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine pao_elem_derivative_2(i_vector,spec,l,nzeta,m,x_i,y_i,z_i,drvtv_val)

    use datatypes
    use numbers

    implicit none
    
    !RC 09/11/03 using (now debugged) routine pp_elem_derivative (see 
    ! above) as template for this sbrt pao_elem_derivative
    real(double), intent(inout) :: x_i,y_i,z_i
    real(double), intent(out) :: drvtv_val
    integer, intent(in) :: i_vector, l, m, spec, nzeta
    integer :: n1,n2
    real(double) :: del_x,xj1,xj2,f_r,df_r,a,b,c,d,alpha,beta,gamma,delta
    real(double) :: x,y,z,r_my,theta,phi,c_r,c_theta,c_phi,x_n,y_n,z_n
    real(double) :: fm_phi, dfm_phi, val1, val2, rl, rl1, tmpr, tmpdr

    !routine to calculate value of gradient of f(r)*Re(Y_lm(x,y,z))
    !along specified Cartesian unit vector
    r_my = sqrt((x_i*x_i)+(y_i*y_i)+(z_i*z_i))
    call pao_rad_vals_fetch(spec,l,nzeta,m,r_my,f_r,df_r)
    if(r_my>very_small) then 
       x_n = x_i/r_my; y_n = y_i/r_my; z_n = z_i/r_my
    else
       x_n = zero; y_n = zero; z_n = zero
    end if
    
    call pp_gradient(i_vector,f_r,df_r,x_n,y_n,z_n,r_my,l,m,drvtv_val)
    
  end subroutine pao_elem_derivative_2
!!***    

!!****f* angular_coeff_routines/pao_rad_vals_fetch *
!!
!!  NAME 
!!   pao_rad_vals_fetch
!!  USAGE
!! 
!!  PURPOSE
!!   Spline routine (probably redundant)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   R. Choudhury
!!  CREATION DATE
!!   2003 sometime
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine pao_rad_vals_fetch(spec,l,nzeta,m,r,f_r,df_r)

    use datatypes
    use numbers
    use pao_format !,ONLY : pao

    implicit none
    !routine to spline out the values of the radial table of a specified PAO
    integer, intent(in) :: spec,l,m,nzeta
    integer :: n1,n2
    real(double), intent(in) :: r
    real(double), intent(inout) :: f_r,df_r
    real(double) :: del_x,a,b,c,d,xj1,xj2
    real(double) :: alpha,beta,gamma,delta,rl,rl1,tmpr,tmpdr
    
    !RC here have to recall that the SIESTA PAO's are normalised by 1/r**l
    del_x = (pao(spec)%angmom(l)%zeta(nzeta)%cutoff/&
         &(pao(spec)%angmom(l)%zeta(nzeta)%length-1))
    !spline out radial table(r) and d1_radial_table(r)
    !initialise to zero
    n1 = 0; n2 = 0; xj1 = zero; xj2 = zero;a=zero
    b = zero; c = zero; d = zero
    
    n1 = floor(r/del_x)
    n2 = n1+1
    xj1 = n1*del_x
    xj2 = xj1+del_x
    
    if(n2<pao(spec)%angmom(l)%zeta(nzeta)%length) then
       a = (xj2-r)/del_x
       b = (r-xj1)/del_x
       c = (a*a*a-a)*del_x*del_x/six
       d = (b*b*b-b)*del_x*del_x/six
       
       alpha = -one/del_x
       beta = -alpha
       gamma = -del_x*(three*a*a-one)/six
        delta = del_x*(three*b*b-one)/six
       
       f_r = a*pao(spec)%angmom(l)%zeta(nzeta)%table(n1+1)+&
            &b*pao(spec)%angmom(l)%zeta(nzeta)%table(n2+1)+&
            &c*pao(spec)%angmom(l)%zeta(nzeta)%table2(n1+1)+&
            &d*pao(spec)%angmom(l)%zeta(nzeta)%table2(n2+1)
       
       df_r = alpha*pao(spec)%angmom(l)%zeta(nzeta)%table(n1+1)& 
            &+beta*pao(spec)%angmom(l)%zeta(nzeta)%table(n2+1)&
            &+gamma*pao(spec)%angmom(l)%zeta(nzeta)%table2(n1+1)&
            &+delta*pao(spec)%angmom(l)%zeta(nzeta)%table2(n2+1)
       
    else
       f_r = zero
       df_r = zero
    end if
    !RC ended splining..
    !Now changing normalised PAO and PAO derivative into unnormalised values
    if(l>0) then
       rl = r**l
       if(l == 1) then
          rl1 = one
       else
          rl1 = r**(l-1)
       endif
       !if(l==1.AND.r_my<1.0e-8) rl1 = one
       tmpr = rl*f_r
       tmpdr = rl*df_r + l*rl1*f_r
       f_r = tmpr
       df_r = tmpdr
    end if
    
  end subroutine pao_rad_vals_fetch
!!***

  subroutine pp_gradient(i_vector,f_r,df_r,x_n,y_n,z_n,r,l,m,d_val)
    use datatypes
    use GenComms, ONLY: inode, ionode, myid !for debugging purposes
    implicit none
    !calculate pp gradients for arbitrary directions,
    !so long as they're Cartesian
    integer, intent(in) :: i_vector,l,m
    real(double), intent(in) :: f_r,df_r,x_n,y_n,z_n,r
    real(double), intent(out) :: d_val
    real(double) :: x,y,z,out_val
    integer :: arg,dir
    !also testing make_xy_pao_gradient_comps subroutine

    !setting up variables for make_xy..
    if(m.lt.0) then 
       arg = -1
    else
       arg = 1
    endif
    x = x_n*r; y = y_n*r; z = z_n*r
    if(i_vector.eq.1) then
       dir = 1
       call pp_gradient_x(f_r,df_r,x_n,y_n,z_n,r,l,m,d_val)
    else if(i_vector.eq.2) then
       dir = -1
       call pp_gradient_y(f_r,df_r,x_n,y_n,z_n,r,l,m,d_val)
    else !i_vector must be 3
       call pp_gradient_z(f_r,df_r,x_n,y_n,z_n,r,l,m,d_val)
    endif
    
  end subroutine pp_gradient

  subroutine pp_gradient_x(f_r,df_r,x_n,y_n,z_n,r,l,m,d_val)
    use datatypes
    implicit none
    !routine to evaluate gradient of f(r)*Re(Ylm(thet,phi)) in the x
    !direction
    real(double), intent(in) :: f_r,df_r,x_n,y_n,z_n,r
    real(double), intent(out) :: d_val
    integer, intent(in) :: l,m
    integer :: dir
    real(double) :: x,y,z,r_my,theta,phi,c_r,c_theta,c_phi
    real(double) :: fm_phi, dfm_phi, val1, val2,denom1,denom2
    real(double) x_comp_1,x_comp_2,x_comp_3,x_comp_4,coeff1,coeff2

    !unnormalize the direction cosines
    x = x_n*r
    y = y_n*r
    z = z_n*r
    !convert from Cartesian system to spherical polars
    call convert_basis(x,y,z,r_my,theta,phi)
    !if zero radius then lets get out of here..
    if(abs(r).lt.very_small) then
       d_val = zero
       return
    else
       continue
    endif
    dir = 1 !telling fix_dval which gradient we are taking
    call make_grad_prefacs(x_comp_1,x_comp_2,x_comp_3,x_comp_4,coeff1,coeff2,df_r,f_r,r,l,m)
    !write(io_lun,*) x_comp_1,x_comp_2,x_comp_3,x_comp_4,'x_comps'
    !now actually calculate and sum the gradient components 
    !RC need to fix the equations for the case where l=0
    if(abs(theta).lt.epsilon) then
       !write(io_lun,*) 'theta small condition satisfied'
       x_comp_1 = x_comp_1*coeff1*sph_hrmnc_z(l+1,m+1,theta)
       x_comp_3 = -x_comp_3*coeff1*sph_hrmnc_z(l+1,m-1,theta)
       x_comp_2 = -x_comp_2*coeff2*sph_hrmnc_z(l-1,m+1,theta)
       x_comp_4 = x_comp_4*coeff2*sph_hrmnc_z(l-1,m-1,theta)
       d_val =  x_comp_1+x_comp_2+x_comp_3+x_comp_4
       call fix_dval(m,d_val,dir)
    else if(theta.lt.pi+epsilon.and.theta.gt.pi-epsilon) then !theta near pi
       x_comp_1 = x_comp_1*coeff1*sph_hrmnc_z(l+1,m+1,theta)
       x_comp_2 = -x_comp_2*coeff2*sph_hrmnc_z(l-1,m+1,theta)
       x_comp_3 = -x_comp_3*coeff1*sph_hrmnc_z(l+1,m-1,theta)
       x_comp_4 = x_comp_4*coeff2*sph_hrmnc_z(l-1,m-1,theta)
       d_val = x_comp_1+x_comp_2+x_comp_3+x_comp_4
       call fix_dval(m,d_val,dir)
    else
       x_comp_1 =  x_comp_1*coeff1*prefac(l+1,m+1)*assoclegpol(l+1,m+1,theta)
       x_comp_2 = -x_comp_2*coeff2*prefac(l-1,m+1)*assoclegpol(l-1,m+1,theta)
       x_comp_3 = -x_comp_3*coeff1*prefac(l+1,m-1)*assoclegpol(l+1,m-1,theta)
       x_comp_4 =  x_comp_4*coeff2*prefac(l-1,m-1)*assoclegpol(l-1,m-1,theta)
       call pp_gradient_x_multiply(x_comp_1,x_comp_2,x_comp_3,x_comp_4,m,phi,d_val)
    endif

  end subroutine pp_gradient_x
  
  subroutine pp_gradient_x_multiply(x_comp_1,x_comp_2,x_comp_3,x_comp_4,m,phi,d_val)
    use datatypes
    implicit none
    !routine to choose between post-mutiplying by sines or cosines depending
    !on the sign of m
    integer, intent(in) :: m
    real(double), intent(in) :: x_comp_1,x_comp_2,x_comp_3,x_comp_4,phi
    real(double), intent(inout) :: d_val
    
    if(m.ge.0) then
       d_val = (cos((m+1)*phi)*(x_comp_1+x_comp_2))+(cos((m-1)*phi))*(x_comp_3+x_comp_4)
       if(m.eq.0) then
          d_val = (cos((m+1)*phi)*(x_comp_1+x_comp_2))+(cos((m-1)*phi))*(x_comp_3+x_comp_4)
          !write(io_lun,*) cos((m+1)*phi),x_comp_1+x_comp_2,cos((m-1)*phi),x_comp_3+x_comp_4
          d_val = d_val/sqrt(two)
          !write(io_lun,*) d_val, 'd_val'
       else
          continue
       endif
    else !m.lt.0
       d_val = (sin((m+1)*phi)*(x_comp_1+x_comp_2))+((sin((m-1)*phi))*(x_comp_3+x_comp_4))
    endif
  end subroutine pp_gradient_x_multiply

  subroutine pp_gradient_y(f_r,df_r,x_n,y_n,z_n,r,l,m,d_val)
    use datatypes
    implicit none
    !routine to evaluate gradient of f(r)*Re(Ylm(thet,phi)) in the y
    !direction
    real(double), intent(in) :: f_r,df_r,x_n,y_n,z_n,r
    real(double), intent(out) :: d_val
    integer, intent(in) :: l,m
    integer :: dir
    real(double) :: x,y,z,r_my,theta,phi,c_r,c_theta,c_phi
    real(double) :: fm_phi, dfm_phi, val1, val2,denom1,denom2
    real(double) y_comp_1,y_comp_2,y_comp_3,y_comp_4,coeff1,coeff2

    !unnormalise the direction cosines
    x = x_n*r
    y = y_n*r
    z = z_n*r
    !find the relevant spherical polar coordinates
    call convert_basis(x,y,z,r_my,theta,phi)
    ! Set d_val to zero for r=0
    if(abs(r).lt.very_small) then
       d_val = zero
       return
    else
       continue
    endif
    call make_grad_prefacs(y_comp_1,y_comp_2,y_comp_3,y_comp_4,coeff1,coeff2,df_r,f_r,r,l,m)
    dir = 2 !telling fix_dval which direction we are taking the gradient in
    if(abs(theta).lt.epsilon) then
       !write(io_lun,*) 'theta small condition satisfied'
       y_comp_1 =  y_comp_1*coeff1*sph_hrmnc_z(l+1,m+1,theta)
       y_comp_2 = -y_comp_2*coeff2*sph_hrmnc_z(l-1,m+1,theta)
       y_comp_3 =  y_comp_3*coeff1*sph_hrmnc_z(l+1,m-1,theta)
       y_comp_4 = -y_comp_4*coeff2*sph_hrmnc_z(l-1,m-1,theta)
       d_val = y_comp_1+y_comp_2+y_comp_3+y_comp_4
       call fix_dval(m,d_val,dir)
    else if(theta.lt.pi+epsilon.and.theta.gt.pi-epsilon) then !theta near pi
       y_comp_1 = y_comp_1*coeff1*sph_hrmnc_z(l+1,m+1,theta)
       y_comp_2 = -y_comp_2*coeff2*sph_hrmnc_z(l-1,m+1,theta)
       y_comp_3 = y_comp_3*coeff1*sph_hrmnc_z(l+1,m-1,theta)
       y_comp_4 = -y_comp_4*coeff2*sph_hrmnc_z(l-1,m-1,theta)
       d_val = y_comp_1+y_comp_2+y_comp_3+y_comp_4
       call fix_dval(m,d_val,dir)
    else
       y_comp_1 =  y_comp_1*coeff1*prefac(l+1,m+1)*assoclegpol(l+1,m+1,theta)
       y_comp_2 = -y_comp_2*coeff2*prefac(l-1,m+1)*assoclegpol(l-1,m+1,theta)
       y_comp_3 =  y_comp_3*coeff1*prefac(l+1,m-1)*assoclegpol(l+1,m-1,theta)
       y_comp_4 = -y_comp_4*coeff2*prefac(l-1,m-1)*assoclegpol(l-1,m-1,theta)
       call pp_gradient_y_multiply(y_comp_1,y_comp_2,y_comp_3,y_comp_4,m,phi,d_val)
    endif
  end subroutine pp_gradient_y

  subroutine pp_gradient_y_multiply(y_comp_1,y_comp_2,y_comp_3,y_comp_4,m,phi,d_val)
    use datatypes
    implicit none
    !routine to decide whether to post-multiply y derivative with cos(phi) or
    !with sin(phi)
    integer, intent(in) :: m
    real(double), intent(in) :: y_comp_1,y_comp_2,y_comp_3,y_comp_4,phi
    real(double), intent(inout) :: d_val
    
    if(m.ge.0) then
       d_val = (sin((m+1)*phi)*(y_comp_1+y_comp_2))+(sin((m-1)*phi))*(y_comp_3+y_comp_4)
       if(m.eq.0) then
          d_val = (sin((m+1)*phi)*(y_comp_1+y_comp_2))+(sin((m-1)*phi))*(y_comp_3+y_comp_4)
          d_val = d_val/sqrt(two)
       else
          continue
       endif
    else !m.lt.0
       d_val = -(cos((m+1)*phi)*(y_comp_1+y_comp_2))-((cos((m-1)*phi))*(y_comp_3+y_comp_4))
    endif
    
  end subroutine pp_gradient_y_multiply

  !! 2016/11/2 10:29 dave Fixing bug in forces for zero distance
  subroutine pp_gradient_z(f_r,df_r,x_n,y_n,z_n,r,l,m,d_val)
    use datatypes
    implicit none
    !code to calculate PAO gradients in the z direction..
    integer, intent(in) :: l,m
    real(double), intent(in) :: f_r,df_r,x_n,y_n,z_n,r
    real(double), intent(out) :: d_val
    real(double) :: x,y,z,pref1,pref2,coeff1,coeff2,z_comp_1,z_comp_2
    real(double) :: r_my,theta,phi
    real(double) :: denom1,denom2
    
    !unnormalize direction cosines and convert to spherical coordinates
    x = x_n*r; y = y_n*r; z = z_n*r
    call convert_basis(x,y,z,r_my,theta,phi)
    ! Set d_val to zero for r = zero
    if(abs(r).lt.very_small) then
       d_val = zero
       return
    endif
    denom1 = (two*l+one)*(two*l+three)
    denom2 = (two*l-one)*(two*l+one)
    pref1 = sqrt((((l+one)**2)-(m**2))/denom1)
    pref2 = sqrt(((l**2)-(m**2))/denom2)
    ! The commented lines below were present when we didn't exit on small r
    !if(r.gt.epsilon) then
    coeff1 = (df_r-(l*f_r/r)); coeff2 = (df_r+((l+one)*f_r/r))
    !else
    !   coeff1 = df_r; coeff2 = df_r
    !end if
    z_comp_1 = pref1*coeff1; z_comp_2 = pref2*coeff2
    if(abs(theta).lt.epsilon) then
       !write(io_lun,*) 'theta small condition satisfied'
       if(m.ge.0) then
          z_comp_1 = z_comp_1*sph_hrmnc_z(l+1,m,theta)
          z_comp_2 = z_comp_2*sph_hrmnc_z(l-1,m,theta)
          d_val = z_comp_1+z_comp_2
          if(m.gt.0) then
             d_val = d_val*sqrt(two)
          else
             d_val = d_val
          endif
       else !m.lt.0
          d_val = zero
       endif
    else
       z_comp_1 = z_comp_1*prefac(l+1,m)*assoclegpol(l+1,m,theta)
       z_comp_2 = z_comp_2*prefac(l-1,m)*assoclegpol(l-1,m,theta)
       if(m.gt.0) then
          d_val = sqrt(two)*cos(m*phi)*(z_comp_1+z_comp_2)
       else if(m.eq.0) then 
          d_val = (z_comp_1+z_comp_2)
       else !m.lt.0
          d_val = sqrt(two)*sin(m*phi)*(z_comp_1+z_comp_2)
       endif
    endif
    
  end subroutine pp_gradient_z
  
  subroutine make_grad_prefacs(pref1,pref2,pref3,pref4,coeff1,coeff2,df_r,f_r,r,l,m)
    use datatypes
    implicit none
    !just code to make simple gradient prefactors
    real(double), intent(in) :: df_r,f_r,r
    real(double), intent(out) :: pref1,pref2,pref3,pref4,coeff1,coeff2
    integer, intent(in) :: l,m
    real(double) :: denom1,denom2

    !need to take care of the situation at the origin
    if(r.gt.epsilon) then
       coeff1 = df_r-(l*f_r/r)
       coeff2 = df_r+(((l+1)/r)*f_r)
    else
       coeff1 = df_r; coeff2 = df_r
    endif

    denom1 = two*(two*l+one)*(two*l+three)
    denom2 = two*(two*l-one)*(two*l+one)
    
    pref1 = sqrt((l+m+one)*(l+m+two)/denom1)
    pref2 = sqrt((l-m-one)*(l-m)/denom2)
    pref3 = sqrt((l-m+one)*(l-m+two)/denom1)
    pref4 = sqrt((l+m-one)*(l+m)/denom2)   

  end subroutine make_grad_prefacs
  
  subroutine fix_dval(m,d_val,dir)
    use datatypes
    implicit none
    !routine to shorten the pp_gradient_i routines
    integer, intent(in) :: m,dir
    real(double), intent(inout) :: d_val
    if(dir.eq.1) then
       if(m.gt.0) then
          d_val = d_val
       else
          if(m.eq.0) then
             d_val = d_val/sqrt(two)
          else
             d_val = zero
             continue
          endif
       endif
    else if(dir.eq.2) then
       if(m.lt.0) then
          d_val = -d_val
       else
          if(m.eq.0) then
             d_val = d_val/sqrt(two)
          else
             d_val = zero
             continue
          endif
       endif
    endif
    
  end subroutine fix_dval
    
!!****f* angular_coeff_routines/assoclegpol *
!!
!!  NAME 
!!   assoclegpol
!!  USAGE
!!   assoclegpol(l,m,x)
!!  PURPOSE
!!   Function to return value of associated Legendre Polynomial 
!!   having parameters l,m.  Uses standard recursion relations:
!!    (l-m)P^{m}_{l} = x(2l-1)P^{m}_{l-1} - (l+m-1)P^{m}_{l-2}
!!    with
!!    P^{m}_{m} = (-1)^{m}(2m-1)!!(1-x^2)^{m/2} with x = cos(theta)
!!    P^{m}_{m+1} = x(2m+1)P^{m}_{m}
!!  INPUTS
!!   l,m, parameters describing Associated Legendre Polynomial
!!   x - argument of associated Legendre polynomial
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2007/01/25 08:21 dave
!!    Removed two assoc. leg. pol routines and combined into one
!!   2019/10/28 13:03 dave
!!    Tidied to clarify recursion and removed (-1)**m
!!  SOURCE
!!
  function assoclegpol(l,m1,theta)
    use datatypes
    use numbers
    implicit none
    !function to generate the value of the l,mth associated
    !Legendre polynomial with argument x
    real(double) :: assoclegpol
    real(double) :: val_mm,val_mm1,val_mm2
    real(double) :: sto_x,factor,factor2,sto_x_test
    real(double) :: x,x_dashed,theta
    integer i,j,k
    integer l,m1,m

    if(l<0) then
       assoclegpol = zero
       return
    endif
    
    m = abs(m1)
    if(m>l) then
       assoclegpol = zero
       return
    endif
    ! Starting value m=l
    val_mm = one
    x = cos(theta)
    sto_x = sin(theta)
    if(m>0) then
       do i=1,2*m-1,2 ! This is (2m-1)!!
          val_mm = val_mm*i
       enddo
       do i=1,m
          val_mm = val_mm*sto_x
       enddo
    end if ! m>0
    if(m==l) then
       assoclegpol = val_mm
    else
       val_mm1 = x*(two*m+1)*val_mm ! l=m+1
       if(m+1==l) then
          assoclegpol = val_mm1
       else
          !now use recursion relation to find val_ml(x)
          if(l>m+1) then
             do i = m+2,l
                val_mm2 = one/(i-m)*((x*(two*i-1)*val_mm1) - ((i+m-1)*val_mm))
                val_mm = val_mm1
                val_mm1 = val_mm2
             enddo
             assoclegpol = val_mm2
          else
             assoclegpol = one
          end if
       end if ! m+1==l
    end if ! m==l

    ! Conversion from P_{lm} with positive m to the required value with negative m

    if(m/=m1) then !input m1 must be negative
       if(2*(m/2)/=m) then ! Odd so factor -1
          factor2 = -(fact(l-m))/(fact(l+m))
       else
          factor2 = (fact(l-m))/(fact(l+m))
       end if
       assoclegpol = assoclegpol*factor2
    else
       continue
    endif
  end function assoclegpol
!!***

!!****s* angular_coeff_routines/set_prefac *
!!
!!  NAME 
!!   set_prefac
!!  USAGE
!!   set_prefac
!!  PURPOSE
!!   Calculates prefactor of Wigner 3-j coefficients and stores
!!  INPUTS
!!   l,m angular momenta
!! 
!!  USES
!!   datatypes
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   2005 sometime
!!  MODIFICATION HISTORY
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine set_prefac(n)

    use datatypes
    use numbers, ONLY: one, four, pi

    implicit none

    integer :: n
    
    integer :: ll, mm
    real(double) :: g, h

    allocate(prefac(-1:n,-n:n))
    prefac(:,:) = one
    do ll = 0, n!lmax_prefac
       do mm = -ll, ll
          g = fact(ll-mm)/fact(ll+mm)
          h = (2*ll+1)/(four*pi)
          prefac(ll,mm) = sqrt(g*h)
       enddo
    enddo
    return
  end subroutine set_prefac
!!***

!!****s* angular_coeff_routines/set_prefac_real *
!!
!!  NAME 
!!   set_prefac_real
!!  USAGE
!!   set_prefac_real
!!  PURPOSE
!!   Calculates prefactor of Wigner 3-j coefficients and stores
!!  INPUTS
!!   l,m angular momenta
!! 
!!  USES
!!   datatypes
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   2005 sometime
!!  MODIFICATION HISTORY
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine set_prefac_real(n)

    use datatypes
    use numbers, ONLY: one, two, four, twopi

    implicit none

    integer :: n
    
    integer :: ll, mm
    real(double) :: g, h

    allocate(prefac_real(-1:n,-n:n))
    prefac_real(:,:) = one
    do ll = 0, n!lmax_prefac
       do mm = -ll, ll
          g = fact(ll-mm)/fact(ll+mm)
          if(mm/=0) then
             h = (2*ll+1)/(twopi)
          else
             h = (2*ll+1)/(two*twopi)
          endif

          prefac_real(ll,mm) = sqrt(g*h)
          !write(io_lun,*) 'pfr: ',ll,mm,ll-mm,ll+mm,2*ll+1,prefac_real(ll,mm)
       enddo
    enddo
    return
  end subroutine set_prefac_real
!!***

!!****f* angular_coeff_routines/re_sph_hrmnc *
!!
!!  NAME 
!!   re_sph_hrmnc
!!  USAGE
!!   re_sph_hrmnc(l,m,theta,phi)
!!  PURPOSE
!!   calculates value of real spherical harmonic
!!  INPUTS
!!   l,m angular momenta
!!   theta, phi - spherical polar coords, z = rcos(theta)
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   25/09/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!

  function re_sph_hrmnc(l,m,theta,phi)
    use datatypes
    use numbers
    implicit none
    !function to return the value of real spherical harmonic
    !at given value of the solid angle
    integer :: i
    integer, intent(in) :: l,m
    real(double), intent(in) :: theta,phi
    real(double) :: re_sph_hrmnc,val,val_chk,beta,del,val_chk2
    
    val = zero
    val = prefac_real(l,m)*assoclegpol(l,m,theta)
    if(m.lt.0) then
       val = val*sin(m*phi)
    else
       val = val*cos(m*phi)
    endif
    re_sph_hrmnc = val

  end function re_sph_hrmnc
!!***


!!****s* angular_coeff_routines/set_fact *
!!
!!  NAME 
!!   set_fact
!!  USAGE
!!   set_fact
!!  PURPOSE
!!   Calculates factorial and stores
!!  INPUTS
!! 
!!  USES
!!   datatypes
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   2005 sometime
!!  MODIFICATION HISTORY
!!   2008/06/10 ast
!!    Added timers
!!   2018/01/22 12:43 dave
!!    Added allocation of arrays and double factorial (n*(n-2)*...)
!!  SOURCE
!!
  subroutine set_fact(lmax)
    use datatypes
    use numbers, ONLY: one

    implicit none

    integer :: lmax
    
    real(double) :: xx
    integer :: ii, max_fact

    max_fact = 4*lmax+1
    !if(max_fact<2*lmax_prefac) max_fact = 2*lmax_prefac
    if(max_fact<(20+2*lmax+1))  max_fact = 20+2*lmax+1
    !write(*,*) 'max_fact is ',max_fact,lmax
    allocate(fact(-1:max_fact))
    fact(-1:max_fact) = one
    do ii=2, max_fact
       xx = real(ii,double)
       fact(ii)=fact(ii-1)* xx
    enddo
    ! Double factorial
    allocate(doublefact(-1:max_fact))
    doublefact(-1:max_fact) = one
    do ii=2,max_fact
       xx = real(ii,double)
       doublefact(ii) = doublefact(ii-2)*xx
    end do
    return
  end subroutine set_fact
!!***

  function sph_hrmnc_z(l,m,theta)
    use datatypes
    use numbers
    implicit none
    !function to return value of spherical harmonic on the z-axis
    integer, intent(in) :: l,m
    real(double), intent(in) :: theta
    real(double) :: sph_hrmnc_z,part,signfac

    sph_hrmnc_z = zero
    if(l<0) return
    if(abs(m)<=l) then
       if(m==0) then
          part = sqrt((two*l+one)/(four*pi))
          if(abs(theta).lt.epsilon) then !theta roughly zero
             sph_hrmnc_z = part
          else
             !theta is in vicinity of pi
             if(theta < pi+epsilon .and. theta > pi-epsilon) then
                signfac = one
                if(2*(l/2)/=l) signfac = -one
                sph_hrmnc_z = signfac*part
             endif
          endif
       endif
    endif
    
  end function sph_hrmnc_z

end module angular_coeff_routines



