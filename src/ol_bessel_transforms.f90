! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module bessel_integrals
! ------------------------------------------------------------------------------
! Code area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/bessel_integrals *
!!  NAME
!!   bessel_integrals
!!  PURPOSE
!!   Holds routines to calculate spherical Bessel transforms of data arrays.
!!   using mixture of FFT subroutines and radial integrals
!!   in real space (FFT's diverge near origin) to evaluate 
!!   spherical Bessel transforms of basis functions.
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2007/01/09 08:28 dave
!!    Tidying, incorporating changes from TM
!!   2007/08/15 12:04 dave
!!    Changed lmax_fact to 22 to accomodate f functions
!!   2008/02/06 08:28 dave
!!    Changed for output to file not stdout
!!   2018/01/22 12:46 dave
!!    Completely new approach for Bessel transform and changes to calculate
!!    general Bessel function for NA projector functions
!!   2019/06/19 15:32 dave
!!    Changes to prepare for removal of NR sine/cosine transform routines
!!  SOURCE
!!
module bessel_integrals

  use datatypes
  use global_module, ONLY: io_lun, area_basis

  implicit none

  save

  ! Store factorial (n!) and double factorial (n!!=n*(n-2)*...)
  real(double), allocatable, dimension(:) :: fact, doublefact
  real(double), allocatable, dimension(:,:) :: bess_coeff

!!***
   
contains

!!****f* bessel_integrals/bessloop *
!!
!!  NAME 
!!   bessloop
!!  USAGE
!!   bessloop(dummyin,l,npts,npts_2,deltar,rcut,dummyout)
!!  PURPOSE
!!   Choose Bessel Transform subroutine that matches l-value of input array
!!  INPUTS
!!   dummyin : input dat array
!!   l : angular momentum of input data array
!!   npts : actual data points in dummyin
!!   npts_2 : no of actual+padded data points (= 2^n)
!!   deltar : grid spacing of data points
!!   rcut : real space cut off of dummyin(npts_2)
!!   dummyout : bessel transform of original data 
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2007/01/09 08:09 dave
!!    Dimension of dummyout fixed
!!   2018/01/22 12:51 JST dave
!!    Now calls new Bessel transform routines (based on integral representation of Bessel function)
!!  SOURCE
!!
   subroutine bessloop(dummyin,l,npts,npts_2,deltar,rcut,dummyout,sign)

     use datatypes

     implicit none

     !choosing spherical Bessel transform subroutine 
     !for incoming angular momentum value l
     integer, intent(in) :: l,npts,npts_2,sign
     real(double), intent(in), dimension(npts) :: dummyin
     real(double), intent(out), dimension(npts_2/2) :: dummyout 
     real(double), intent(in) :: deltar,rcut

     if(mod(l,2)==0) then
        call new_bessel_transform_evenFFTW(l,dummyin,npts,npts_2,rcut,deltar,dummyout)
     else 
        call new_bessel_transform_oddFFTW(l,dummyin,npts,npts_2,rcut,deltar,dummyout)
     end if
     return
   end subroutine bessloop
!!***

!!****f* bessel_integrals/twon *
!!
!!  NAME 
!!   twon
!!  USAGE
!!   twon(npts,npts_2)
!!  PURPOSE
!!   Calculate next highest no to npts which is an integer power of 2
!!  INPUTS
!!   npts : original no of points
!!   npts_2 : padded up to 2^n no of points
!!  USES
!!   none
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
   subroutine twon(npts,npts_2)
     implicit none

     integer, intent(in) :: npts
     integer, intent(out) :: npts_2
     integer :: i

     do i=1,30
        if(npts.lt.2**i) then
           npts_2=2**i
           exit ! Leave do loop
        endif
     enddo
   end subroutine twon
!!***

!!****f* bessel_integrals/complx_fctr *
!!
!!  NAME 
!!   complx_fctr
!!  USAGE
!!   complx_fctr(l1,l2,l,factor)
!!  PURPOSE
!!   calculate real prefactors for radial tables
!!  INPUTS
!!   l1,l2,l - angular momentum values of triplets in radial table
!!   factor : real prefactor output
!!  USES
!!   none
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
   subroutine complx_fctr(l1,l2,l,factor)

     use datatypes 
     use numbers, ONLY: zero, one

     implicit none

     complex(double_cplx) :: z1,factor1
     real(double), intent(out) :: factor
     integer, intent(in) :: l1,l2,l

     z1 = cmplx(zero,one,double_cplx)

     factor1 = (z1**l1)*(conjg(z1**l2))*(conjg(z1**l))
     factor = factor1-aimag(factor1)
     !write(io_lun,*) factor1,'factor1',factor,'factor'

   end subroutine complx_fctr
!!***

!!****f* bessel_integrals/maxtwon *
!!
!!  NAME 
!!   maxtwon
!!  USAGE
!!   maxtwon(n1,del1,n2,del2,n12,del12,delk,kcut)
!!  PURPOSE
!!   Calculate size of array required for spherical Bessel transforms
!!  INPUTS
!!   n1,n2 : sizes of input arrays 
!!   del1, del2 : grid spacing of original arrays
!!   n12, del12 : size and grid spacing of required arrays
!!   delk, kcut : k space grid spacing and cut-off
!!  USES
!!   datatypes
!!  AUTHOR
!!   R Choudhury
!!  CREATION DATE
!!   24/07/03
!!  MODIFICATION HISTORY
!!   2019/10/31 14:11 dave
!!    Replaced call to twon with simple log calculation
!!  SOURCE
!!
   subroutine maxtwon(n1,del1,n2,del2,n12,del12,delk,kcut)

     use datatypes
     use numbers, ONLY: two, twopi

     implicit none
     !code to get upper bound on size of the two
     !arrays as integral power of 2
     real(double),intent(in):: del1,del2,delk,kcut
     real(double), intent(inout) :: del12
     integer :: n1,n2,n12,s3_max,s3i

     real(double) :: s3

     !choosing according to restrictions (of user)
     !on minimum delk, min kcut...
     del12 = twopi/(two*kcut)
     s3 = two*(kcut)/delk
     s3i = floor(log(s3)/log(two)) + 1
     n12 = 2**s3i
     return
   end subroutine maxtwon
!!***

!!****f* bessel_integrals/multiply_ksq *
!!
!!  NAME 
!!   multiply_ksq
!!  USAGE
!!   multiply_ksq(y,n,dk)
!!  PURPOSE
!!   Multiply input array by k squared
!!  INPUTS
!!   y : input data array
!!   n : no of points in input array
!!   dk : grid spacing of input array
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
   subroutine multiply_ksq(y,n,dk)

     use datatypes

     implicit none

     integer, intent(in) :: n
     real(double),intent(inout), dimension(n) :: y
     real(double), intent(in) ::dk
     real(double) :: k
     integer i
     
     do i=1,n
        k = real(i-1,double)*dk
        y(i) = k*k*y(i)
     enddo

   end subroutine multiply_ksq
!!***

!!****f* bessel_integrals/general_bessel *
!!
!!  NAME 
!!   general_bessel
!!  USAGE
!!   general_bessel(x)
!!  PURPOSE
!!   Calculates a general Bessel function using either recursion or
!!   (for small r) series expansion
!!  INPUTS
!!   r : value of r
!!   n : order of Bessel function
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   22/01/18
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
   function general_bessel(r,n)

     use datatypes
     use numbers
     use GenComms, only: cq_abort
     
     implicit none

     real(double) :: general_bessel
     ! Passed
     real(double) :: r
     integer :: n
     
     ! Local
     real(double), dimension(0:n) :: sph_bess
     integer :: i, s
     real(double) :: term
     logical :: flag_series

     ! Do we need series expansion ? These numbers should have a more quantitative basis !
     flag_series = .false.
     if(n>2.AND.r<0.05_double) flag_series = .true.
     if(n>3.AND.r<0.2_double) flag_series = .true.
     if(n>4.AND.r<0.4_double) flag_series = .true.
     if(n>5.AND.r<0.7_double) flag_series = .true.
     if(n>6.AND.r<1.0_double) flag_series = .true.
     if(n>7.AND.r<1.6_double) flag_series = .true.
     if(n>8.AND.r<2.5_double) flag_series = .true.
     if(n>9.AND.r<3.0_double) flag_series = .true.
     if(n>10.AND.r<3.4_double) flag_series = .true.
     if(n<0) call cq_abort("Error: Can't have spherical bessel with order less than zero ",n)
     ! Find Bessel function based on r and need for series expansion
     if(abs(r)<very_small) then
        if(n>0) then
           general_bessel = zero
        else
           general_bessel = one - r*r/six ! Yes, I know... It should be accurate
        end if
     else if(flag_series) then
        general_bessel = zero
        do s=0,5
           term = ((-1)**s)*fact(s+n)*(r**(2*s)) /(fact(s)*fact(2*s+2*n+1))
           general_bessel = general_bessel + term
        enddo
        do i=1,n
           general_bessel = general_bessel * two * r
        end do
     else
        sph_bess(0) = sin(r)/r
        if(n==0) then
           general_bessel = sph_bess(0)
        else
           sph_bess(1) = (sph_bess(0) - cos(r))/r
           if(n==1) then
              general_bessel = sph_bess(1)
           else
              do i=2,n
                 sph_bess(i) = sph_bess(i-1)*real(2*i-1,double)/r - sph_bess(i-2)
              end do
              general_bessel = sph_bess(n)
           end if
        end if
     end if
   end function general_bessel
!!***
   
!!****f* bessel_integrals/new_bessel_transform_even *
!!
!!  NAME 
!!   new_bessel_transform_even
!!  USAGE
!!   new_bessel_transform_even
!!  PURPOSE
!!   Calculates a general Bessel transform for even order
!!  INPUTS
!!   n : order of Bessel function
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   22/01/18
!!  MODIFICATION HISTORY
!!   2019/06/19 15:29 dave
!!    Changed to use FFTW
!!  SOURCE
!!
   subroutine new_bessel_transform_evenFFTW(n,function_in,npts,npts2,rcut,delta_r,function_out)

     use datatypes
     use numbers
     use fft_interface_module, ONLY: cosft_init_wrapper, cosft_exec_wrapper, &
          sinft_init_wrapper, sinft_exec_wrapper, cosft_dest_wrapper, sinft_dest_wrapper

     implicit none

     ! Passed variables
     integer :: n, npts, npts2
     real(double) :: rcut, delta_r
     real(double), dimension(npts) :: function_in !function_out
     real(double), dimension(npts2/2) :: function_out

     ! Local variables
     integer :: i, j, en, l, poly_order, xi
     real(double) :: k, fac, dk, dk2, dk3, x1,x0,x0_in, x1_in, pmj, r, tee_nm, rcut_large
     real(double), dimension(:), allocatable :: dummy1, dummy2, ess, coeff_poly, coeff_poly1, prefac

     poly_order = 3
     ! en is found such that n = 2*en
     en = floor(half*n)
     ! Allocate
     allocate(dummy1(npts2+1),dummy2(npts2),ess(0:en),coeff_poly(0:poly_order), &
          coeff_poly1(0:poly_order),prefac(0:en))
     dummy1 = zero
     dummy2 = zero
     ess = zero
     ! We start the loop from i=2 because i=1 gives zero.
     do i=2,npts
        r = real(i-1,double) * delta_r
        dummy1(i) = r*r*function_in(i)*delta_r
        dummy2(i-1) = r*r*r*function_in(i)*delta_r
     end do
     call sinft_init_wrapper(npts2-1)
     call cosft_init_wrapper(npts2+1)
     ! Cosine transform for function
     call cosft_exec_wrapper(dummy1,npts2+1,-1)
     ! Sine transform for first derivative (for polynomial fitting)
     call sinft_exec_wrapper(dummy2,npts2-1,-1)
     ! Adjust the array back and set to zero
     do j=npts2,2,-1
        dummy2(j) = -half*dummy2(j-1)
     end do
     dummy2(1) = zero
     dummy1 = half*dummy1
     ! Set k-space interval - the k-space grid is pi/rcut_large for compatibility
     ! with FFT routines which implicitly double the number of points
     dk = twopi/(rcut+delta_r)
     dk2 = dk*dk
     dk3 = dk2*dk
     ! Calculate prefactors for Legendre polynomial
     pmj = one
     do j=0,en
        prefac(j) = pmj*doublefact(2*en+2*j-1)/(fact(2*j)*doublefact(2*en-2*j))
        pmj = -pmj
     end do
     ! Loop to perform transform
     do i=3,npts2,2
        k = real((i-1)/2,double)*dk
        fac = one/k
        x0_in = (k-dk)
        x1_in = k
        ! Cubic coefficients in powers of (x-x0)
        coeff_poly1(0) = dummy1(i-2)
        coeff_poly1(1) = dummy2(i-2)
        coeff_poly1(2) = three*(dummy1(i)-dummy1(i-2))/dk2 - (dummy2(i) + two*dummy2(i-2))/dk
        coeff_poly1(3) =  -two*(dummy1(i)-dummy1(i-2))/dk3 + (dummy2(i) + dummy2(i-2))/dk2
        ! Convert to x
        coeff_poly(3) = coeff_poly1(3)
        coeff_poly(2) = coeff_poly1(2) - three*coeff_poly1(3)*(k-dk)
        coeff_poly(1) = coeff_poly1(1) - two*coeff_poly1(2)*(k-dk) + &
             three*(k-dk)*(k-dk)*coeff_poly1(3)
        coeff_poly(0) = coeff_poly1(0) - coeff_poly1(1)*(k-dk) + &
             (k-dk)*(k-dk)*coeff_poly1(2) - (k-dk)*(k-dk)*(k-dk)*coeff_poly1(3)
        function_out((i+1)/2) = zero
        ! Build Bessel transform
        do j=0,en
           ! k^{2j+1} 
           x0 = x0_in
           x1 = x1_in
           ! Interpolate to find T_nm
           tee_nm = zero
           do xi=0,poly_order 
              tee_nm = tee_nm + coeff_poly(xi)*(x1-x0)/real(2*j+xi+1,double)
              x1 = x1*k
              x0 = x0*(k-dk)
           end do
           x0_in = x0_in*(k-dk)*(k-dk)
           x1_in = x1_in*k*k
           ! Accumulate into S_nm
           ess(j) = ess(j) + tee_nm
           ! Sum over n to get transform
           function_out((i+1)/2) = function_out((i+1)/2) + prefac(j)*ess(j)*fac
           fac = fac/(k*k) ! 1/k^{2j+1}
        end do
     end do
     function_out(1) = zero
     if(n==0) then ! Quadrature for k=0
        do i=1,npts
           r = (i-1)*delta_r
           function_out(1) = function_out(1)+(r*r*delta_r*function_in(i))
        enddo
     end if
     call cosft_dest_wrapper
     call sinft_dest_wrapper
   end subroutine new_bessel_transform_evenFFTW
!!***

!!****f* bessel_integrals/new_bessel_transform_odd *
!!
!!  NAME 
!!   new_bessel_transform_odd
!!  USAGE
!!   new_bessel_transform_odd
!!  PURPOSE
!!   Calculates a general Bessel transform for odd order
!!  INPUTS
!!   n : order of Bessel function
!!  USES
!!   
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   22/01/18
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
   subroutine new_bessel_transform_oddFFTW(n,function_in,npts,npts2,rcut,delta_r,function_out)

     use datatypes
     use numbers
     use fft_interface_module, ONLY: cosft_init_wrapper, cosft_exec_wrapper, &
          sinft_init_wrapper, sinft_exec_wrapper, cosft_dest_wrapper, sinft_dest_wrapper

     implicit none

     ! Passed variables
     integer :: n, npts, npts2
     real(double) :: rcut, delta_r
     real(double), dimension(npts) :: function_in  !, function_out
     real(double), dimension(npts2/2) :: function_out

     ! Local variables
     integer :: i, j, en, l, poly_order, xi
     real(double) :: k, fac, dk, dk2, dk3, x1,x0,x0_in, x1_in, pmj, r, tee_nm, rcut_large
     real(double), dimension(:), allocatable :: dummy1, dummy2, ess, coeff_poly, prefac, coeff_poly1

     poly_order = 3
     ! en is found such that 2*en+1
     en = floor(half*n)
     ! Allocate
     allocate(dummy1(npts2),dummy2(npts2+1),ess(0:en),coeff_poly(0:poly_order), &
          coeff_poly1(0:poly_order),prefac(0:en))
     dummy1 = zero
     dummy2 = zero
     ess = zero
     do i=2,npts
        r = real(i-1,double) * delta_r
        dummy1(i-1) = r*r*delta_r*function_in(i)
        dummy2(i) = r*r*r*delta_r*function_in(i)
     end do
     call sinft_init_wrapper(npts2-1)
     call cosft_init_wrapper(npts2+1)
     ! Sine transform for function
     call sinft_exec_wrapper(dummy1,npts2-1,-1)
     do j=npts2,2,-1
        dummy1(j) = half*dummy1(j-1)
     end do
     dummy1(1) = zero
     ! Cosine transform for first derivative (for polynomial fitting)
     call cosft_exec_wrapper(dummy2,npts2+1,-1)
     dummy2 = half*dummy2
     ! Set k-space interval
     dk = twopi/(rcut+delta_r)
     dk2 = dk*dk
     dk3 = dk2*dk
     ! Calculate prefactors for Legendre polynomial
     pmj = one
     do j=0,en
        prefac(j) = pmj*doublefact(2*en+2*j+1)/(fact(2*j+1)*doublefact(2*en-2*j))
        pmj = -pmj
     end do
     ! Loop to perform transform
     do i=3,npts2, 2
        k = real((i-1)/2,double)*dk
        fac = one/(k*k)
        x0_in = (k-dk)*(k-dk)
        x1_in = k*k
        ! Cubic coefficients in powers of (x-x0)
        coeff_poly1(0) = dummy1(i-2)
        coeff_poly1(1) = dummy2(i-2)
        coeff_poly1(2) = three*(dummy1(i)-dummy1(i-2))/dk2 - (dummy2(i) + two*dummy2(i-2))/dk
        coeff_poly1(3) =  -two*(dummy1(i)-dummy1(i-2))/dk3 + (dummy2(i) + dummy2(i-2))/dk2
        ! Convert to x
        coeff_poly(3) = coeff_poly1(3)
        coeff_poly(2) = coeff_poly1(2) - three*coeff_poly1(3)*(k-dk)
        coeff_poly(1) = coeff_poly1(1) - two*coeff_poly1(2)*(k-dk) + &
             three*(k-dk)*(k-dk)*coeff_poly1(3)
        coeff_poly(0) = coeff_poly1(0) - coeff_poly1(1)*(k-dk) + &
             (k-dk)*(k-dk)*coeff_poly1(2) - (k-dk)*(k-dk)*(k-dk)*coeff_poly1(3)
        function_out((i+1)/2) = zero
        ! Build Bessel transform
        do j=0,en
           ! k^{2j+1 +1}
           x0 = x0_in
           x1 = x1_in
           ! Interpolate to find T_nm
           tee_nm = zero
           do xi=0,poly_order ! Interpolation to find T_nm
              tee_nm = tee_nm + coeff_poly(xi)*(x1-x0)/real(2*j+1+xi+1,double)
              x1 = x1*k
              x0 = x0*(k-dk)
           end do
           x0_in = x0_in*(k-dk)*(k-dk)
           x1_in = x1_in*k*k
           ! Accumulate into S_nm
           ess(j) = ess(j) + tee_nm
           ! Sum over n to get transform
           function_out((i+1)/2) = function_out((i+1)/2) + prefac(j)*ess(j)*fac
           fac = fac/(k*k) ! 1/k^{2j+1+1}
        end do
     end do
     function_out(1) = zero
     call cosft_dest_wrapper
     call sinft_dest_wrapper
   end subroutine new_bessel_transform_oddFFTW
   
 end module bessel_integrals
