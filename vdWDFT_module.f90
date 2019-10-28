! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module vdwDFT_module
! ------------------------------------------------------------------------------
! Code area 3: Operators
! ------------------------------------------------------------------------------

!!****h* Conquest/vdwDFT_module
!! NAME
!!   vdwDFT_module
!! PURPOSE
!!   Collects together all routines used for ab initio Van der Waals
!!   fucntional, introduced by
!!     Dion et al. Phys. Rev. Lett., 2004, 92, 246401
!!   And with practical implemtation introduced by
!!     Roman-Perez and Soler, Phys. Rev. Lett., 2009, 103, 096102
!! METHODS
!!   vdW_decusp    : Energy due to the softening of the VdW kernel cusp
!!   vdW_theta     : Finds function theta_q(rho,grad_rho)
!!   vdW_get_qmesh : Returns size and values of q-mesh
!!   vdW_phi       : Finds and interpolates phi(q1,q2,k)
!!   vdW_set_kcut  : Sets the planewave cutoff kc of the integration grid
!! AUTHOR
!!   L.Tong
!! CREATION DATE
!!   2012/02/21
!! MODIFICATION HISTORY
!! SOURCE
!!
module vdWDFT_module

  use datatypes
  use numbers
  use global_module, only: area_ops

  implicit none

  save

  ! methods
  public ::           &
       vdW_decusp,    &! Energy due to the softening of the VdW kernel cusp
       vdW_theta,     &! Finds function theta_q(rho,grad_rho)
       vdW_get_qmesh, &! Returns size and values of q-mesh
       vdW_phi,       &! Finds and interpolates phi(q1,q2,k)
       vdW_set_kcut,  &! Sets the planewave cutoff kc of the integration grid
       vdWXC_energy,  &! Calculates the XC energy under Dion et al.'s vdW scheme
       vdwxc_energy_slow

  private  ! all module parameters are private beyond this point

  character(len=80) :: &
       RCSid = "$Id$"

  ! Precision parameters for the integral defining phi in routine phi_val
  real(double), parameter :: acutmin = 10.0_double  ! Upper integration limit
  real(double), parameter :: acutbyd = 30.0_double  ! Upper integration limit / d
  real(double), parameter :: damax = 0.5_double     ! Max. integration interval
  real(double), parameter :: damin = 1.e-2_double   ! Min. integration interval
  real(double), parameter :: dabyd = 0.1_double     ! Min. integration interval / d

  ! Precision parameter for the integral in routine dphi
  ! Shorter integration limit / d
  real(double), parameter :: ashortbyd = 2.5_double

  ! Parameters for phi_soft. Some recommended pairs of values are
  ! (dsoft,phi0_soft)=(0.5,0.8)|(0.7|0.6)|(1.0|0.4)|(1.5,0.22)|(2.0,0.12)
  ! Softening matching radius
  real(double), parameter :: dsoft = 1.0_double
  ! phi_soft(0,0) (depends on dsoft)
  real(double), parameter :: phi0_soft = 0.40_double
  ! Delta(d) for derivatives
  real(double), parameter :: dd = 0.01_double

  ! Mesh parameters for table phi(d1,d2)
  ! Number of d mesh points
  integer,      parameter :: nd = 20
  ! Max. value of d mesh
  real(double), parameter :: dcut = 30.0_double
  ! Last d mesh interval / first one
  real(double), parameter :: ddmaxddmin = 20.0_double

  ! Use routine dphi for better efficiency in setting phi table?
  ! .true. => efficiency | .false. => accuracy
  logical, parameter :: use_dphi = .true.

  ! Set derivation methods to use for interpolation table
  character(len=*), parameter :: deriv_method = 'numeric' !('numeric'|'interp')
  character(len=*), parameter :: interp_method= 'spline'  !('lagrange'|'spline')

  ! Parameters to find numerical derivatives of phi, to set interpolation
  real(double), parameter :: ddbyd = 0.01_double  ! Delta to find phi derivs / d
  real(double), parameter :: ddmin = 0.001_double ! Min. delta to find phi derivs

  ! Mesh parameters for table of phi(q1,q2,r) and its Fourier transform
  ! Radial mesh points (power of 2)
  integer,      parameter :: nr = 1024
  ! Total number of q mesh points
  integer,      parameter :: mq = 30
  ! Effective number of q mesh points
  integer,      parameter :: nq = mq - 1
  ! Max. value of q mesh
  real(double), parameter :: qcut = 5.0_double
  ! Last q mesh interval / first one
  real(double), parameter :: dqmaxdqmin = 20.0_double
  ! Radial cutoff: r>rcut => phi = 0
  real(double), parameter :: rcut = 100.0_double
  ! Min. radius as denominator
  real(double), parameter :: rmin = 1.e-6_double
  ! Soften kernel in r<rsoft
  real(double), parameter :: rsoft = zero

  ! Parameters for cutoff function, used in radial Fourier transforms of phi
  ! cutoff(x) = (1 - x**ncut1)**ncut2
  integer, parameter :: ncut1 = 8
  integer, parameter :: ncut2 = 4

  ! Parameters for saturate function, used to enforce that q < qcut
  integer, parameter :: nsat  = 12  ! xsat(x) = 1 - exp(-sum_n=1:nsat x**n/n)

  ! Parameters for saturate_inverse function
  real(double), parameter :: xmaxbyxc = 1.5_double ! qmax/qcut
  real(double), parameter :: ytol = 1.e-15_double  ! Tol. for saturated q

  ! Private module variables and arrays
  ! Mesh points for phi(d1,d2) table
  real(double), dimension(nd) :: dmesh
  ! Mesh points for phi(q1,q2,r)
  real(double), dimension(mq) :: qmesh
  ! Coefs. for bicubic interpolation
  real(double), dimension(0:3, 0:3, nd, nd) :: phi_table
  ! Has phi_table been set?
  logical :: phi_table_set = .false.
  ! Has qmesh been set?
  logical :: qmesh_set = .false.
  ! Has kcut been set?
  logical :: kcut_set = .false.
  ! Table of phi(r), r is radial coordinate
  real(double), dimension(1:nr, mq, mq) :: phir
  ! Table of d2_phi/dr2
  real(double), dimension(1:nr, mq, mq) :: d2phidr2
  ! r-mesh interval
  real(double) :: dr
  ! Table of phi(k), fourier transform of phir w,r,t r
  real(double), dimension(1:nr, mq, mq) :: phik
  ! Table of d2_phi/dk2
  real(double), dimension(1:nr, mq, mq) :: d2phidk2
  ! k-mesh interval
  real(double) :: dk
  ! Planewave cutoff: k > kcut => phi = 0
  real(double) :: kcut
  ! # k points within kcut
  integer :: nk

!!*****

contains


  ! !!****f* vdWDFT_module/bessph
  ! !! PURPOSE
  ! !!   Returns the spherical bessel function J_l(x)
  ! !!   Ref. Abramowitz and Stegun, Formulas 10.1.2 and 10.1.19
  ! !! USAGE
  ! !!   J_l(x) = bessph(l,x)
  ! !! INPUTS
  ! !!   integer      l : angular momentum number, l
  ! !!   real(double) x : variable of the bessel function
  ! !! RETURN VALUE
  ! !!   real(double) bessph : J_l(x)
  ! !! AUTHOR
  ! !!   L.Tong
  ! !! CREATION DATE
  ! !!   2012/04/07
  ! !! MODIFICATION HISTORY
  ! !! SOURCE
  ! !!
  ! function bessph(l, x)

  !   use datatypes
  !   use numbers
  !   use global_module, only: io_lun
  !   use GenComms,      only: cq_abort, inode, ionode

  !   implicit none

  !   ! passed parameters
  !   integer,      intent(in) :: l
  !   real(double), intent(in) :: x
  !   ! result
  !   real(double) :: bessph
  !   ! local variables
  !   integer,      parameter :: nterms = 100
  !   real(double), parameter :: tiny = 1.0e-15_double
  !   integer      :: ii, nn
  !   real(double) :: l_d, ii_d, nn_d
  !   real(double) :: switch, term, x2, sum_t, y, fnm1, fn, fnp1
  !   logical      :: done

  !   l_d = real(l, double)

  !   switch = max(one, two * l_d - one)
  !   if (abs(x) < switch) then
  !      ! use power series
  !      term = one
  !      do ii = 1, l
  !         ii_d = real(ii, double)
  !         term = term * x / (two * ii_d + 1)
  !      end do
  !      x2 = x * x
  !      sum_t = zero
  !      done = .false.
  !      do ii = 1, nterms
  !         ii_d = real(ii, double)
  !         sum_t = sum_t + term
  !         term = - term * x2 / (two * ii_d * (two * ii_d + two * l_d + one))
  !         if (abs(term) < tiny) then
  !            bessph = sum_t
  !            done = .true.
  !            exit
  !         end if
  !      end do
  !      if (.not. done) then
  !         write (io_lun, '(1x,a,i2,a,f15.7)') 'l = ',l,' x = ',x
  !         call cq_abort('bessph: series has not converged, see line above.')
  !      end if
  !   else
  !      ! use explicit expressions or recurrence relation
  !      if (l == 0) then
  !         bessph = sin(x) / x
  !      else if (l == 1) then
  !         bessph = (sin(x) / x - cos(x)) / x
  !      else
  !         y = one / x
  !         fnm1 = sin(x) * y
  !         fn = (fnm1 - cos(x)) * y
  !         do nn = 1, l-1
  !            nn_d = real(nn, double)
  !            fnp1 = (two * nn_d + one) * y * fn - fnm1
  !            fnm1 = fn
  !            fn = fnp1
  !         end do
  !         bessph = fn
  !      end if
  !   end if

  ! end function bessph
  ! !!*****


  ! !!****f* vdWDFT_module/four1
  ! !! PURPOSE
  ! !!   Discrete 1D Fourier transform. Modified and converted to
  ! !!   complex double precision from the same routine in Numerical
  ! !!   Recipes
  ! !! USAGE
  ! !!   call four1(data, n, isign)
  ! !! INPUTS
  ! !!   integer      n         : number of points for function data,
  ! !!                            must be a power of two
  ! !!   real(double) data(2*n) : function to transform, complex values stored in
  ! !!                            the format: (\ 1r, 1i, 2r, 2i, 3r, 3i, ... \)
  ! !!   integer      isign     : isign = +1 | -1 ==> direct/inverse transform
  ! !! OUTPUT
  ! !!   real(double) data(2*n) : the direct Fourier transform (isign = +1), or
  ! !!                            n times the inverse Fourier transform (isign = -1)
  ! !! RETURN VALUE
  ! !! AUTHOR
  ! !!   L.Tong
  ! !! CREATION DATE
  ! !!   2012/04/07
  ! !! MODIFICATION HISTORY
  ! !! SOURCE
  ! !!
  ! subroutine four1(data, n, isign)

  !   use datatypes
  !   use numbers

  !   implicit none

  !   ! passed parameters
  !   integer, intent(in) :: n, isign
  !   real(double), dimension(2*n) :: data

  !   integer      :: ii, istep, jj, mm, mmax, nn
  !   real(double) :: tempi, tempr, theta, wi, wpi, wpr, wr, wtemp

  !   nn = 2 * n
  !   jj = 1
  !   do ii = 1, nn, 2
  !      if (jj > ii) then
  !         tempr = data(jj)
  !         tempi = data(jj+1)
  !         data(jj) = data(ii)
  !         data(jj+1) = data(ii+1)
  !         data(ii) = tempr
  !         data(ii+1) = tempi
  !      end if
  !      mm = nn / 2
  !      do
  !         if ((mm < 2) .or. (jj <= mm)) exit
  !         jj = jj - mm
  !         mm = mm / 2
  !      end do
  !      jj = jj + mm
  !   end do ! ii
  !   mmax = 2
  !   do ! until (n <= mmax)
  !      if (nn <= mmax) exit
  !      istep = 2 * mmax
  !      theta = twopi / real((isign * mmax), double)
  !      wpr = -two * sin(half * theta)**2
  !      wpi = sin(theta)
  !      wr = one
  !      wi = zero
  !      do mm = 1, mmax, 2
  !         do ii = mm, nn, istep
  !            jj = ii + mmax
  !            tempr = wr * data(jj) - wi * data(jj+1)
  !            tempi = wr * data(jj + 1) + wi * data(jj)
  !            data(jj) = data(ii) - tempr
  !            data(jj+1) = data(ii+1) - tempi
  !            data(ii) = data(ii) + tempr
  !            data(ii+1) = data(ii+1) + tempi
  !         end do ! ii
  !         wtemp = wr
  !         wr = wr * wpr - wi * wpi + wr
  !         wi = wi * wpr + wtemp * wpi + wi
  !      end do ! mm
  !      mmax = istep
  !   end do ! until (n <= mmax)

  ! end subroutine four1
  ! !!*****


  !!****f* vdWDFT_module/radfft
  !! PURPOSE
  !!   Calculates the discrete approximation to the continuous radial
  !!   Fourier transform of a given data array. The transform is
  !!   defined as follows:
  !!
  !!   If f(vec_r) = f(r), where r = |vec_r|, and if the Fourier
  !!   transform is defined as:
  !!
  !!              /+infty
  !!   F(vec_k) = |  d^3vec_r exp(-i*2*pi*vec_r.vec_k) f(vec_r)
  !!              /-infty
  !!
  !!   Then it can be shown that:  F(vec_k) = F(k), k = |vec_k| and
  !!
  !!           /+infty   2*sin(2*pi*k*r)
  !!   F(k) =  |      dr --------------- * r * f(r)
  !!           /0               k
  !!
  !!   The inverse transfrom of radial Fourier transform is exactly
  !!   the transform it-self.
  !!
  !!   The subroutine assumes the input data to contain n_r elements
  !!   corresponding to data samples taken from f(r) at
  !!
  !!   r(1:nr) = 0, dr, ..., r_cut - dr.
  !!
  !!   In other words, dr = r_max / n_r.
  !!
  !!   The Nyquist critical frequency (fc) for the discrete data
  !!   sample is defined as fc = 1 / 2*dr. The output (transformed)
  !!   data is thus defined to be n_r data points taken from F(k) at
  !!
  !!   k(1:nr) = 0, dk, ..., fc - dk
  !!
  !!   In other words, kmax = fc - dk, and dk = fc / n_r = 1 / (2*dr*
  !!   n_r). And F(fc) is NOT included in the output. This property
  !!   means the problem when discretised according to the above
  !!   scheme can be written as a discerete sin transform.
  !!
  !!   The k = 0 case is treated separately by doing real-space
  !!   integral directly, due to the 1/k singularity if it is to be
  !!   done with fast sin transfroms.
  !! USAGE
  !!   - Forward transform: f(1:n_r) -> g(1:n_r), f(n_r) = r_cut - dr
  !!       dr = rmax / (n_r - 1)
  !!       call radfft(f, g, n_r, dr)
  !!   - Backward tranform: g(1:n_r) -> f(1:n_r), g(n_r) = fc - dk
  !!       dk = 1/(2*dr) - 1/(2*dr*n_r)
  !!       call radfft(g, f, n_r, dk)
  !! INPUTS
  !!   integer      n_r : number of data points
  !!   real(double) dr  : the uniform interval r(i+1) - r(i), in radial grid
  !!                      or, for reverse transformation, dk = k(i+1) - k(i)
  !!   real(double) func(n_r) : data points to be transformed
  !! OUTPUT
  !!   real(double) func_out(n_r) : transformed data points
  !! NOTES
  !!   - Note although uses discrete transforms (sin FT), this is
  !!     equivalent to the continuous version, and hence has a factor
  !!     of real-space grid spacing multiplied to the final result.
  !!   - Because the subroutine uses discrete transform methods, the
  !!     Nyquist critical frequency (fc) limit still applies. That
  !!     means if the transformed data do not tail to 0 at the end
  !!     (the data is cut off at fc), it means aliasing have probably
  !!     occurred and the result will be unreliable. See Numerical
  !!     Recipes on aliasing and critical frequency.
  !!   - The radial FFT differs from ordinary FFT in the sense that
  !!     for ordinary FFT, k should range from -fc to fc (with fc and
  !!     -fc correspond to aliasing points, and with -(fc - dk), ...,
  !!     -dk elements stored in the n_r/2 + 1 to n_r-th element of the
  !!     output); however for radial FFT, k range from 0 to fc. This
  !!     means dk = 1/(2*dr*n_r) instead of the usual 1/(dr*n_r), and
  !!     there is no 'folding back' of the -fc results as the case
  !!     with oridinary FFTs.
  !!   - For reverse transform, dk should be fc/n_r, with fc the
  !!     Nyquist critical frequency corresponding to dr given by the
  !!     forward transform (fc = 1/(2*dr)). Note that k here is
  !!     defined WITHOUT the factor of 2*pi. So if one works with
  !!     reciprocal vectors q = 2*pi*k, one needs to remember to
  !!     divide by 2*pi at the input. I.e. dk = dq / (2*pi)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/05/08
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine radfft(func, func_out, n_r, dr)
    use datatypes
    use numbers
    use dimens
    use fft_procedures, only: sinft
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use GenComms,       only: cq_abort
    implicit none
    ! passed parameters
    integer,      intent(in) :: n_r
    real(double), intent(in) :: dr
    real(double), dimension(n_r), intent(in)  :: func
    real(double), dimension(n_r), intent(out) :: func_out
    ! local variables
    real(double), dimension(:), allocatable :: data
    integer      :: ir, n_q, ik, stat
    real(double) :: rr, kk, dk

    allocate(data(n_r), STAT=stat)
    if (stat /= 0) call cq_abort("radfft: Error alloc mem: ", nr)
    call reg_alloc_mem(area_ops, nr, type_dbl)

    ! see Numerical Recipies Section 12.1: eq (12.1.2), but this time
    ! k only goes from 0 to fc, fc = 1/(2*dr) is Nyquist critical
    ! frequency. So dk = 1/(2*dr*n_r)
    dk = half / real(dr * n_r, double)
    ! set the real data to 2 * r * func(r)
    do ir = 1, n_r
       rr = (ir - 1) * dr
       data(ir) = two * rr * func(ir)
    end do
    ! do sin FT
    call sinft(data, n_r)
    ! leave the k = 0 point aside, this is treated separately
    do ik = 2, n_r, 1
       kk = (ik - 1) * dk
       func_out(ik) = data(ik) / kk
    end do
    ! doing the k = 0 point
    ! the limit (k -> 0) of 2 * r * sin(2*pi*k*r) / k is 4*pi*r**2
    func_out(1) = zero
    do ir = 1, n_r
       rr = (ir - 1) * dr
       func_out(1) = func_out(1) + (four * pi * rr**2 * func(ir))
    end do
    ! scale func_out by to give approximation to non-discrete form
    func_out = dr * func_out

    deallocate(data, STAT=stat)
    if (stat /= 0) call cq_abort("radfft: Error dealloc mem")
    call reg_dealloc_mem(area_ops, nr, type_dbl)

  end subroutine radfft
  !!*****


  ! !!****f* vdWDFT_module/radfft
  ! !! PURPOSE
  ! !!   Implementation by Soler et al.
  ! !!
  ! !!   Makes a fast Fourier transform of a radial function.
  ! !!   If function f is of the form
  ! !!     f(r_vec) = F(r_mod) * Y_lm(theta,phi)
  ! !!   where Ylm is a spherical harmonic with l = argument L, and
  ! !!   argument F contains on input the real_function F(r_mod), in a
  ! !!   uniform radial grid:
  ! !!     r_mod = ir * rmax / nr,  ir = 0,1,...,nr,
  ! !!   and if g is the 3-dimensional Fourier transform of f:
  ! !!     g(k_vec) = 1/(2*pi)**(3/2) *
  ! !!                Integral(d3_r_vec * exp(-i * k_vec * r_vec) * f(r_vec))
  ! !!   then g has the form
  ! !!     g(k_vec) = (-i)**l * G(k_mod) * Ylm(theta,phi)
  ! !!   where argument G contains on output the real_function G(k_mod)
  ! !!   in a uniform radial grid:
  ! !!     k_mod = ik * k_max / nr, ik = 0,1,...,nr,  k_max = nr * pi / rmax
  ! !!
  ! !!   Ref: J.M.Soler notes of 16/08/95.
  ! !! USAGE
  ! !!   call radfft(l, nr, rmax, F, G)
  ! !! INPUTS
  ! !!   integer l   : Angular momentum quantum number
  ! !!   integer n_r : Number of radial inervals.
  ! !!                2*nr must be an acceptable number of points for
  ! !!                the FFT routine used.
  ! !!   real(double) rmax : maximum radius
  ! !!   real(double) F(0:n_r) : function to be transformed, in radial mesh
  ! !! OUTPUT
  ! !!   real(double) G(0:n_r) : Fourier transform of F (but see point 5 below)
  ! !! NOTES
  ! !! - Units
  ! !!   - Units of rmax and F are arbitrary
  ! !!   - Units of k_max and G are related with those of rmax and F in
  ! !!     the obvious way (see above)
  ! !! - Behaviour
  ! !!   1) F and G may be the same physical array, i.e. it is allowed
  ! !!        call radfft(l, n_r, rmax, F, F)
  ! !!   2) It also works in the opposite direction, but then the factor
  ! !!      multiplying the output is (+i)**l. Thus, the following two
  ! !!      calls
  ! !!        call radfft(l, n_r, rmax, F, G)
  ! !!        call radfft(l, n_r, n_r * pi / rmax, G, H)
  ! !!      make H = F
  ! !!   3) If you will divide the output by q**l, truncation errors may
  ! !!      be quite large for small k's if l and n_r are
  ! !!      large. Therefore, these components are calculated by direct
  ! !!      integration rather than FFT. Parameter ERRFFT is the typical
  ! !!      truncation error in the FFT, and controls which k's are
  ! !!      integrated directly. A good value is 1.0e-8. If you will not
  ! !!      divide by k**l, make ERRFFT = 1.0e-30.
  ! !!   4) The function F is assumed to be zero at and beyond rmax. The
  ! !!      last point F(n_r) is not used to find G. except G(0) for l =
  ! !!      0 (see 5)
  ! !!   5) Because of the 'uncetainty principle', if f(r) is strictly
  ! !!      zero for r > rmax, then g(k) cannot be strictly zero for k >
  ! !!      kmax. Therefore G(n_r), which should be exactly zero, is used
  ! !!      (for l = 0) as a 'reminder' term for the integral of G
  ! !!      beyond kmax, to ensure that F(0) = sum(4 * pi * r**2 * dr *
  ! !!      G(ir)) (this allows to recover F(0) when called again in the
  ! !!      inverse direction). Thus, the last value G(n_r) should be
  ! !!      replaced by zero for any other use. NOTICE: this is
  ! !!      commented out in this version!
  ! !! AUTHOR
  ! !!   L.Tong
  ! !! CREATION DATE
  ! !!   2012/04/08
  ! !! MODIFICATION HISTORY
  ! !! SOURCE
  ! !!
  ! subroutine radfft(l, n_r, rmax, F, G)

  !   use datatypes
  !   use numbers

  !   implicit none

  !   ! passed parameters
  !   integer,      intent(in) :: l
  !   integer,      intent(in) :: n_r
  !   real(double), intent(in) :: rmax
  !   real(double), dimension(0:n_r), intent(in)  :: F
  !   real(double), dimension(0:n_r), intent(out) :: G
  !   ! local variables
  !   real(double), parameter :: ERRFFT = 1.0e-8_double
  !   integer      :: ii, iq, ir, jr, mm, mq, nn, n_q
  !   real(double) :: c, dq, dr, Fr, r, rn, q, qmax
  !   real(double), dimension(0:2*n_r)   :: GG
  !   real(double), dimension(2,0:2*n_r) :: Fn
  !   real(double), dimension(2,0:l,0:l) :: P

  !   ! define some constants
  !   n_q = n_r
  !   dr = rmax / n_r
  !   dq = pi / rmax
  !   qmax = n_q * dq
  !   c = dr / sqrt(two * pi)

  !   ! set up a complex polynomial such that the spherical Bessel function:
  !   !    j_l(x) = Real(sum_n(P(n,l) * x**n) * exp(i*x)) / x**(l+1)
  !   P(1,0,0) = zero
  !   P(2,0,0) = -one
  !   if (l > 0) then
  !      P(1,0,1) = zero
  !      P(2,0,1) = -one
  !      P(1,1,1) = -one
  !      P(2,1,1) = zero
  !      if (l > 1) then
  !         do mm = 2, l
  !            do nn = 0, mm
  !               do ii = 1, 2
  !                  P(ii,nn,mm) = zero
  !                  if (nn < mm) &
  !                       P(ii,nn,mm) = P(ii,nn,mm) + &
  !                                     real(2 * mm - 1, double) * &
  !                                     P(ii,nn,mm-1)
  !                  if (nn >= 2) &
  !                       P(ii,nn,mm) = P(ii,nn,mm) - P(ii,nn-2,mm-2)
  !               end do ! ii
  !            end do ! nn
  !         end do ! mm
  !      end if
  !   end if

  !   ! initialise accumulation array
  !   do iq = 0, n_q
  !      GG(iq) = zero
  !   end do

  !   ! iterate on terms of the j_l(q*r) polynomial
  !   do nn = 0, l
  !      ! set up function to be fast fourier transformed
  !      Fn(1,0) = zero
  !      Fn(2,0) = zero
  !      do jr = 1, 2 * n_r - 1
  !         if (jr < n_r) then
  !            ir = jr
  !            r = ir * dr
  !            Fr = F(ir)
  !         else if (jr == n_r) then
  !            ir = jr
  !            r = ir * dr
  !            Fr = zero
  !         else
  !            ir = 2 * n_r - jr
  !            r = - ir * dr
  !            Fr = F(ir) * (- one)**l
  !         end if
  !         ! find r**2 * r**n / r**(l+1)
  !         rn = r**(nn - l + 1)
  !         Fn(1,jr) = c * Fr * rn * P(1,nn,l)
  !         Fn(2,jr) = c * Fr * rn * P(2,nn,l)
  !      end do ! jr
  !      ! perform one-dimensional complex FFT

  !      ! only the elements from 0 to 2*n_r-1 of Fn are used.
  !      ! (a total of 2*n_r). four1 will receive a one-dimensional array
  !      ! of size 2*n_r

  !      call four1(Fn, 2*n_r, +1)

  !      ! accumulate contribution
  !      do iq = 1, n_q
  !         q = iq * dq
  !         GG(iq) = (GG(iq) + Fn(1,iq)) / q
  !      end do
  !   end do ! nn

  !   ! special case for q = 0
  !   GG(0) = zero
  !   if (l == 0) then
  !      do ir = 1, n_r
  !         r = ir * dr
  !         GG(0) = GG(0) + r * r * F(ir)
  !      end do
  !      GG(0) = GG(0) * two * c
  !   end if

  !   ! Direct integration for the smallest Q's
  !     if (l == 0) then
  !       mq = 0
  !     else
  !       mq = n_q * ERRFFT**(one/real(l, double))
  !     end if
  !     do iq = 1, mq
  !       q = iq * dq
  !       GG(iq) = zero
  !       do ir = 1, n_r
  !         r = ir * dr
  !         GG(iq) = GG(iq) + r * r * F(ir) * bessph(l, q * r)
  !       end do
  !       GG(iq) = GG(iq) * two * c
  !     end do

  !     ! Special case for q = qmax
  !     ! if (l == 0) then
  !     !    G_sum = 0.D0
  !     !    do iq = 1, n_q - 1
  !     !       q = iq * dq
  !     !       G_sum = G_sum + q * q * GG(iq)
  !     !    end do
  !     !    G_sum = G_sum * four * pi * dq
  !     !    GG(n_q) = (two * pi)**1.5_double * F(0) - G_sum
  !     !    GG(n_q) = GG(n_q) / (four * pi * dq * qmax**2)
  !     ! end if

  !     ! Copy from local to output array
  !     do iq = 0, n_q
  !        G(iq) = GG(iq)
  !     end do

  ! end subroutine radfft
  ! !!*****


  !!****f* vdWDFT_module/find_bicubic_coeffs
  !! PURPOSE
  !!   Finds coefficients for bicubic interpolation
  !!   Standard approach (written up in many places,
  !!   including Wikipedia and Numerical Recipes)
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/26
  !! MODIFICATION HISTORY
  !!   2019/10/25 08:17 dave
  !!    Tidying, renaming
  !! SOURCE
  !!
  subroutine find_bicubic_coeffs(n1, n2, x1, x2, y, dydx1, dydx2, d2ydx1dx2, c)

    use datatypes

    implicit none
    
    ! passed variables
    integer,                        intent(in) :: n1, n2
    real(double), dimension(n1),    intent(in) :: x1
    real(double), dimension(n2),    intent(in) :: x2
    real(double), dimension(n1,n2), intent(in) :: y
    real(double), dimension(n1,n2), intent(in) :: dydx1
    real(double), dimension(n1,n2), intent(in) :: dydx2
    real(double), dimension(n1,n2), intent(in) :: d2ydx1dx2
    real(double), dimension(0:3,0:3,n1,n2), intent(out) :: c
    
    ! local variables
    integer      :: i1, i11, i12, i13, i14, i2, i21, i22, i23, i24
    real(double) :: dx1, dx2
    real(double), dimension(16)     :: z
    ! This is the inverse matrix needed to find coefficients for one square
    real(double), dimension(16,16) :: wt
    data wt /1,  0, -3, 2, 4*0, -3, 0, 9, -6, 2, 0, -6, 4, 8*0, 3, 0,    &
             -9, 6, -2, 0, 6, -4, 10*0, 9, -6, 2*0, -6, 4, 2*0, 3, -2,   &
             6*0, -9, 6, 2*0, 6, -4, 4*0, 1, 0, -3, 2, -2, 0, 6, -4,     &
             1, 0, - 3, 2, 8*0, -1, 0, 3, -2, 1, 0, -3, 2, 10*0, -3,     &
             2, 2*0, 3, - 2, 6*0, 3, -2, 2*0, -6, 4, 2*0, 3, -2, 0, 1,   &
             -2, 1, 5*0, -3, 6, -3, 0, 2, -4, 2, 9*0, 3, -6, 3, 0, -2,   &
             4, - 2, 10*0, -3, 3, 2*0, 2, -2, 2*0, -1, 1, 6*0, 3, -3,    &
             2*0, - 2, 2, 5*0, 1, -2, 1, 0, -2, 4, -2, 0, 1, -2, 1, 9*0, &
             -1, 2, -1, 0, 1, -2, 1, 10*0, 1, -1, 2*0, -1, 1, 6*0, -1,   &
             1, 2*0, 2, -2, 2*0, -1, 1/

    ! Set coefs. for i1<n1 and i2<n2
    do i2 = 1, n2 - 1
       do i1 = 1, n1 - 1
          dx1 = x1(i1 + 1) - x1(i1)
          dx2 = x2(i2 + 1) - x2(i2)
          i11 = i1
          i12 = i1 + 1
          i13 = i1 + 1
          i14 = i1
          i21 = i2
          i22 = i2
          i23 = i2 + 1
          i24 = i2 + 1
          z( 1) = y(i11,i21)
          z( 2) = y(i12,i22)
          z( 3) = y(i13,i23)
          z( 4) = y(i14,i24)
          z( 5) = dydx1(i11,i21) * dx1
          z( 6) = dydx1(i12,i22) * dx1
          z( 7) = dydx1(i13,i23) * dx1
          z( 8) = dydx1(i14,i24) * dx1
          z( 9) = dydx2(i11,i21) * dx2
          z(10) = dydx2(i12,i22) * dx2
          z(11) = dydx2(i13,i23) * dx2
          z(12) = dydx2(i14,i24) * dx2
          z(13) = d2ydx1dx2(i11,i21) * dx1 * dx2
          z(14) = d2ydx1dx2(i12,i22) * dx1 * dx2
          z(15) = d2ydx1dx2(i13,i23) * dx1 * dx2
          z(16) = d2ydx1dx2(i14,i24) * dx1 * dx2
          z = matmul(wt,z)
          c(0:3, 0:3, i1, i2) = reshape(z, (/4, 4/), order=(/2, 1/))
       end do ! i1
    end do ! i2

    ! Set c for i1=n1 and i2=n2 (valid only at the border)
    c(:, :, n1, :) = 0
    c(:, :, :, n2) = 0
    c(0, 0, n1, :) = y(n1, :)
    c(0, 0, :, n2) = y(:, n2)

    return
  end subroutine find_bicubic_coeffs
  !!*****


  !!****f* vdWDFT_module/cutoff
  !! PURPOSE
  !!   cutoff(x) = (1 - x**ncut1)**ncut2
  !!   - ncut1 and ncut2 are defined as module variables
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/26
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function cutoff(x)

    use datatypes
    use numbers

    implicit none

    ! passed variables
    real(double), intent(in) :: x
    ! result
    real(double) :: cutoff

    if (x <= zero) then
       cutoff = one
    else if (x >= one) then
       cutoff = zero
    else
       cutoff = (one - x**ncut1)**ncut2
    end if

    return
  end function cutoff
  !!*****


  !!****f* vdWDFT_module/dphi
  !! PURPOSE
  !!   Finds phi(d1, d2) - phi(dmax, dmax), with dmax = max(d1, d2)
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/26
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function dphi(d1, d2)

    use datatypes
    use numbers
    use global_module,  only: area_ops
    use GenComms,       only: cq_abort
    use vdWMesh_module, only: set_mesh, get_mesh, integral, get_n
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! input parameters
    real(double), intent(in) :: d1, d2
    ! result
    real(double) :: dphi
    ! local variables
    integer      :: ia, ib, n, nmesh, nshort, stat
    real(double) :: a, acut, b, da1, dan, deff, dmax, dmin, gamma, t, t0, w
    real(double), dimension(:), allocatable :: amesh, c, dphida, &
                                               dphidb, s, nu0, nu1, &
                                               nu2
    ! Find integration mesh
    dmax = max(abs(d1), abs(d2))
    dmin = min(abs(d1), abs(d2))
    deff = dmax
    acut = max(acutmin, acutbyd * deff)
    da1 = min(damax, dabyd * deff)
    da1 = max(da1, damin)
    dan = damax
    n = get_n(zero, acut, da1, dan)
    allocate(amesh(n), c(n), dphida(n), dphidb(n), s(n), nu0(n), &
             nu1(n), nu2(n), STAT=stat)
    if (stat /= 0) call cq_abort("dphi: failed to allocate memory", n, stat)
    call reg_alloc_mem(area_ops, 8 * n, type_dbl)
    call set_mesh(n, xmax=acut, dxndx1=dan/da1)
    call get_mesh(n, nmesh, amesh)

    ! Find limit of shorter mesh
    nshort = n
    do ia = n - 1, 1, -1
       if (amesh(ia) > ashortbyd * deff) nshort = ia
    end do

    ! Find cos(a), sin(a), nu1(a), and nu2(a)
    ! nu1, nu2, defined in [Dion et al. PRL 92, 246401 (2004)]
    gamma = four * pi / nine
    do ia = 1, n
       a = amesh(ia)
       c(ia) = cos(a)
       s(ia) = sin(a)
       if (a == zero) then
          nu0(ia) = dmax**2 / two / gamma
          nu1(ia) = d1**2 / two / gamma
          nu2(ia) = d2**2 / two / gamma
       else ! (a /= 0)
          if (d1 == zero) then
             nu1(ia) = a * a / two
          else
             nu1(ia) = a * a / two / (one - exp(-gamma * (a / d1)**2))
          end if
          if (d2 == zero) then
             nu2(ia) = a * a / two
          else
             nu2(ia) = a * a / two / (one - exp(-gamma * (a / d2)**2))
          end if
          if (dmax <= zero) then
             nu0(ia) = a * a / two
          else
             nu0(ia) = a * a / two / (one - exp(-gamma * (a / dmax)**2))
          end if
       end if ! (a==0)
    end do

    ! Make integral on variable a
    dphida(1) = zero ! integrand w.r.t a
    do ia = 2, nshort
       a = amesh(ia)
       ! Make integral on variable b
       dphidb(1) = zero ! integrand w.r.t a and b [Dion2004, eq(14)]
       do ib = 2, n
          b = amesh(ib)

          w = two * ((three - a * a) * b * c(ib) * s(ia) +     &
                     (three - b * b) * a * c(ia) * s(ib) +     &
                     (a * a + b * b - three) * s(ia) * s(ib) - &
                     three * a * b * c(ia) * c(ib)) /          &
              (a * b)**3

          t = half * (one / (nu1(ia) + nu1(ib)) +                       &
                      one / (nu2(ia) + nu2(ib))) *                      &
                     (one / (nu1(ia) + nu2(ia)) / (nu1(ib) + nu2(ib)) + &
                      one / (nu1(ia) + nu2(ib)) / (nu2(ia) + nu1(ib)))

          t0 = half * (one / (nu0(ia) + nu0(ib)) +                       &
                       one / (nu0(ia) + nu0(ib))) *                      &
                      (one / (nu0(ia) + nu0(ia)) / (nu0(ib) + nu0(ib)) + &
                       one / (nu0(ia) + nu0(ib)) / (nu0(ia) + nu0(ib)))

          dphidb(ib) = a * a * b * b * w * (t - t0)
    end do ! ib
    dphida(ia) = integral(n, dphidb) - integral(ia, dphidb)
  end do ! ia

  ! note that here, we use the symmetry of ia and ib, so that integral
  ! is only done in a trianble, 0 to max for ia and then ia to max for
  ! ib, and then the result is multiplied by two
  dphi = two / pi**2 * two * integral(nshort, dphida)

  deallocate(amesh, c, dphida, dphidb, s, nu0, nu1, nu2, STAT=stat)
  if (stat /= 0) call cq_abort('dphi: failed to deallocate arrays', stat)
  call reg_dealloc_mem(area_ops, 8 * n, type_dbl)

  return
  end function dphi
  !!*****


  !!****f* vdWDFT_module/dphi_fast
  !! PURPOSE
  !!   Finds phi(d1,d2)-phi(dmax,dmax), with dmax=max(d1,d2), by
  !!   - Direct integration if d < dsoft, where d=sqrt(d1**2+d2**2)
  !!   - Interpolation of phi_table if d > dsoft
  !! USAGE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/28
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function dphi_fast(d1, d2)

    use datatypes
    use numbers

    implicit none

    ! input parameters
    real(double), intent(in) :: d1, d2
    ! results
    real(double) :: dphi_fast
    ! local variables
    real(double) :: d, dmax

    if (.not. phi_table_set) call set_phi_table()

    d = sqrt(d1**2 + d2**2)
    dmax = max(d1, d2)

    if (d < dsoft) then
       dphi_fast = dphi(d1, d2)
    else
       dphi_fast = phi_interp(d1, d2) - phi_interp(dmax, dmax)
    end if

    return
  end function dphi_fast
  !!*****


  !!****f* vdWDFT_module/dphi_soft
  !! PURPOSE
  !!   Finds phi_soft(d1,d2) - phi_soft(dmax,dmax), with dmax = max(d1,d2)
  !! USAGE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/28
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function dphi_soft(d1, d2)

    use datatypes

    implicit none

    ! passed parameters
    real(double), intent(in) :: d1, d2
    ! result
    real(double) :: dphi_soft
    ! local variables
    real(double) :: d, dmax

    d = sqrt(d1**2 + d2**2)
    dmax = max(d1, d2)

    if (d < dsoft) then
       dphi_soft = phi_soft(d1, d2) - phi_soft(dmax, dmax)
    else
       dphi_soft = dphi(d1, d2)
    end if

  end function dphi_soft
  !!*****


  !!****f* vdWDFT_module/iofd
  !! PURPOSE
  !!   Finds index i such that dmesh(i) <= d < dmesh(i+1)
  !! USAGE
  !!   index = iofd(d)
  !! INPUTS
  !!   real(double) d : value of d
  !! RETURN VALUE
  !!   integer iofd   : result
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function iofd(d)

    use datatypes

    implicit none

    ! input parameters
    real(double), intent(in) :: d
    ! result
    integer :: iofd
    ! local variables
    real(double), parameter :: amin = 1.e-12_double
    real(double), save :: a, b
    logical,      save :: first_call = .true.

    if (first_call) then
       a = log((dmesh(nd) - dmesh(nd-1)) / (dmesh(2) - dmesh(1))) / &
           real((nd - 2), double)
       a = max(a, amin)
       b = (dmesh(2) - dmesh(1)) / (exp(a) - one)
       first_call = .false.
    end if

    ! note that int truncates 1.7 to 1, and nint truncates 1.7 to 2
    iofd = int(one + log(one + (d - dmesh(1)) / b) / a)
    iofd = max(1, iofd)
    iofd = min(nd - 1, iofd)

  end function iofd
  !!*****


  !!****f* vdWDFT_module/iofq
  !! PURPOSE
  !!   Finds index i such that qmesh(i) <= q < qmesh(i+1)
  !! USAGE
  !!   index = iofq(q)
  !! INPUTS
  !!   real(double) q : value of q
  !! RETURN VALUE
  !!   integer iofq : result
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function iofq(q)

    use datatypes

    implicit none

    ! input parameter
    real(double), intent(in) :: q
    ! result
    integer :: iofq
    ! local variables
    real(double), parameter :: amin = 1.e-12_double
    real(double), save :: a, b
    logical,      save :: first_call = .true.

    if (first_call) then
       a = log((qmesh(mq) - qmesh(mq-1)) / (qmesh(2) - qmesh(1))) / &
           real((mq - 2), double)
       a = max(a, amin)
       b = (qmesh(2) - qmesh(1)) / (exp(a) - one)
       first_call = .false.
    end if

    ! note that int truncates 1.7 to 1, and nint truncates 1.7 to 2
    iofq = int(one + log(one + (q - qmesh(1)) / b ) / a)
    iofq = max(1, iofq)
    iofq = min(mq-1, iofq)

  end function iofq
  !!*****


  !!****f* vdWDFT_module/phi_fast
  !! PURPOSE
  !!   Finds hard phi(d1, d2) kernel by
  !!   - Direct integration if d < dsoft, where d = sqrt(d1**2 + d2**2)
  !!   - Interpolation of phi_table if d > dsoft
  !! USAGE
  !!   phi(d1,d2) = phi_fast(d1, d2)
  !! INPUTS
  !!   real(double) d1, d2 :: values of d1 and d2
  !! RETURN VALUE
  !!   real(double) phi_fast :: value of phi given d1 and d2
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function phi_fast(d1, d2)

    use datatypes

    implicit none

    ! input parameter
    real(double), intent(in) :: d1, d2
    ! result
    real(double) :: phi_fast
    ! local variable
    real(double) :: d

    if (.not. phi_table_set) call set_phi_table()

    d = sqrt(d1**2 + d2**2)

    if (d < dsoft) then
       phi_fast = phi_val(d1, d2)
    else
       phi_fast = phi_interp(d1, d2)
    end if

  end function phi_fast
  !!*****


  !!****f* vdWDFT_module/phi_interp
  !! PURPOSE
  !!   Finds soft phi(d1, d2) kernel by interpolation of phi_table
  !! USAGE
  !!   phi(d1,d2) = phi_interp(d1, d2)
  !! INPUTS
  !!   real(double) d1, d2 : values of d1 and d2
  !! RETURN VALUE
  !!   real(double) phi_interp : value of phi from interpolation
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function phi_interp(d1, d2)

    use datatypes

    implicit none

    ! input parameters
    real(double), intent(in) :: d1, d2
    ! result
    real(double) :: phi_interp
    ! local variables
    integer :: i1, i2, id1, id2
    real(double), dimension(0:3) :: dd1, dd2

    if (.not. phi_table_set) call set_phi_table()

    if (d1 >= dcut .or. d2 >= dcut) then
       phi_interp = zero
       return
    end if

    id1 = iofd(d1)
    id2 = iofd(d2)

    dd1(0) = one
    dd2(0) = one
    dd1(1) = (d1 - dmesh(id1)) / (dmesh(id1+1) - dmesh(id1))
    dd2(1) = (d2 - dmesh(id2)) / (dmesh(id2+1) - dmesh(id2))
    dd1(2) = dd1(1)**2
    dd2(2) = dd2(1)**2
    dd1(3) = dd1(1)**3
    dd2(3) = dd2(1)**3

    phi_interp = zero
    do i2 = 0, 3
       do i1 = 0, 3
          phi_interp = phi_interp + phi_table(i1,i2,id1,id2) * dd1(i1) * dd2(i2)
       end do
    end do

  end function phi_interp
  !!*****


  !!****f* vdWDFT_module/phiofr
  !! PURPOSE
  !!   Finds phi(q1,q2,r) = phi(q1*r,q2*r) with q1=qmesh(iq1), q2=qmesh(iq2),
  !!   by interpolation from phi_table
  !!   Notice that q is a density parameter, relates to the Fermi wavevector
  !! USAGE
  !!   call phiofr(r, phi)
  !! INPUTS
  !!   real(double) r   : radial value, r of phi(q1, q2, r)
  !! OUTPUT
  !!   real(double) phi : value of phi(1:mq, 1:mq) w.r.t r
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine phiofr(r, phi)

    use datatypes
    use GenComms,      only: cq_abort

    implicit none

    ! input parameter
    real(double), intent(in) :: r
    ! result
    real(double), dimension(:,:), intent(out) :: phi
    ! local variables
    integer      :: iq1, iq2, j
    real(double) :: dphidr, a, b, c, d, r1, r2, r3, r4, rr
    logical      :: stat

    if (size(phi, 1) < mq .or. size(phi, 2) < mq) &
         call cq_abort('phiofr: ERROR: size(phi) too small', size(phi), mq)
    if (.not. qmesh_set) call set_qmesh()

    ! rcut is not included in the radial data set, and hence
    ! interpolation should only be betweem 0 to rcut - dr
    if (r >= rcut - dr) then
       phi(:,:) = zero
    else
       do iq2 = 1, mq
          do iq1 = 1, iq2
             ! Use unfiltered kernel
             ! phi(iq1,iq2) = phi_interp(qmesh(iq1)*r, qmesh(iq2)*r)
             ! Use filtered kernel
             j = floor(r/dr)+1
             if(j+1<=nr) then
                rr = real(j,double)*dr
                a = (rr - r)/dr
                b = one - a
                c = a * ( a * a - one ) * dr * dr / six
                d = b * ( b * b - one ) * dr * dr / six
                r1 = phir(j,iq1,iq2)
                r2 = phir(j+1,iq1,iq2)
                r3 = d2phidr2(j,iq1,iq2)
                r4 = d2phidr2(j+1,iq1,iq2)
                phi(iq1,iq2) = a*r1 + b*r2 + c*r3 + d*r4
                phi(iq2,iq1) = phi(iq1,iq2)
             end if
          end do ! iq1
       end do ! iq2
    end if ! (r>=rcut)

    return
  end subroutine phiofr
  !!*****


  !!****f* vdWDFT_module/phi_soft
  !! PURPOSE
  !!   Finds phi(d1, d2) softened near d1 = d2 = 0
  !! USAGE
  !!   phi(d1,d2) = phi_soft(d1, d2)
  !! INPUTS
  !!   real(double) d1, d2   : values of d1, d2
  !! RETURN VALUE
  !!   real(double) phi_soft : value of phi
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function phi_soft(d1, d2)

    use datatypes

    implicit none

    ! input parameter
    real(double), intent(in) :: d1, d2
    ! result
    real(double) :: phi_soft
    ! local variables
    real(double) :: d, d1m, d1p, d2m, d2p, dphidd, &
                    phi, phi0, phi2, phi4, phim, phip

    d = sqrt(d1**2 + d2**2)

    if (d <= zero) then
       phi_soft = phi0_soft
    else if (d > dsoft) then
       phi_soft = phi_val(d1, d2)
    else ! (0 < d < dsoft)
       d1p = d1 / d * (dsoft + dd)
       d2p = d2 / d * (dsoft + dd)
       d1m = d1 / d * (dsoft - dd)
       d2m = d2 / d * (dsoft - dd)
       phip = phi_val(d1p, d2p)
       phim = phi_val(d1m, d2m)
       phi = (phip + phim) / two
       dphidd = (phip - phim) / (two * dd)
       phi0 = phi0_soft
       phi2 = (four * (phi - phi0) - dphidd * dsoft) / (two * dsoft**2)
       phi4 = (two * (phi0 - phi)  + dphidd * dsoft) / (two * dsoft**4)
       phi_soft = phi0 + phi2 * d**2 + phi4 * d**4
    end if ! (d <= zero)

  end function phi_soft
  !!*****


  !!****f* vdWDFT_module/phi_val
  !! PURPOSE
  !!   Finds kernel phi by direct integration of Eq.(14) of
  !!   Dion et al. PRL 92, 246401 (2004)
  !! USAGE
  !!   phi(d1, d2) = phi_val(d1, d2)
  !! INPUTS
  !!   real(double) d1, d2  : values of d1 and d2
  !! RETURN VALUE
  !!   real(double) phi_val : value of phi
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  function phi_val(d1, d2)

    use datatypes
    use global_module,  only: area_ops
    use GenComms,       only: cq_abort
    use vdWMesh_module, only: set_mesh, get_mesh, integral, get_n
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none
    ! input parameters
    real(double), intent(in) :: d1, d2
    ! result
    real(double) :: phi_val
    ! local variables
    integer      :: ia, ib, n, nmesh, stat
    real(double) :: a, acut, b, da1, dan, deff, dmax, dmin, gamma, t, w
    real(double), dimension(:), allocatable :: amesh, c, dphida, &
                                               dphidb, s, nu1, nu2

    ! Find integration mesh
    dmax = max(abs(d1), abs(d2))
    dmin = min(abs(d1), abs(d2))
    deff = sqrt(d1**2 + d2**2)
    acut = max(acutmin, acutbyd * deff)
    da1 = min(damax, dabyd * deff)
    da1 = max(da1, damin)
    dan = damax
    n = get_n(zero, acut, da1, dan)
    allocate(amesh(n), c(n), dphida(n), dphidb(n), s(n), nu1(n), &
             nu2(n), STAT=stat)
    if (stat /= 0) call cq_abort('phi_val: failed to allocate arrays', stat, n)
    call reg_alloc_mem(area_ops, 7 * n, type_dbl)
    call set_mesh(n, xmax=acut, dxndx1=dan/da1)
    call get_mesh(n, nmesh, amesh)

    ! Find cos(a), sin(a), nu1(a), and nu2(a)
    gamma = four * pi / nine
    do ia = 1, n
       a = amesh(ia)
       c(ia) = cos(a)
       s(ia) = sin(a)
       if (a == zero) then
          nu1(ia) = d1**2 / two / gamma
          nu2(ia) = d2**2 / two / gamma
       else ! (a/=0)
          if (d1 == zero) then
             nu1(ia) = a * a / two
          else
             nu1(ia) = a * a / two / (one - exp(-gamma * (a / d1)**2))
          end if
          if (d2 == zero) then
             nu2(ia) = a * a / two
          else
             nu2(ia) = a * a / two / (one - exp(-gamma * (a / d2)**2))
          end if
       end if ! (a==0)
    end do

    ! Make integral on variable a
    dphida(1) = zero
    do ia = 2, n
       a = amesh(ia)

       ! Make integral on variable b
       dphidb(1) = zero
       do ib = 2, n
          b = amesh(ib)

          w = two * ((3 - a*a) * b * c(ib) * s(ia) + &
                     (three - b*b) * a * c(ia) * s(ib) + &
                     (a*a + b*b - three) * s(ia) * s(ib) - &
                     three * a * b * c(ia) * c(ib)) / (a * b)**3

          t = half * (one / (nu1(ia) + nu1(ib)) + one / (nu2(ia) + nu2(ib))) * &
                     (one / (nu1(ia) + nu2(ia)) / (nu1(ib) + nu2(ib)) + &
                      one / (nu1(ia) + nu2(ib)) / (nu2(ia) + nu1(ib)))

          dphidb(ib) = a*a * b*b * w * t

       end do ! ib
       dphida(ia) = integral(n, dphidb)

    end do ! ia
    phi_val = two / pi**2 * integral(n, dphida)

    deallocate(amesh, c, dphida, dphidb, s, nu1, nu2, STAT=stat)
    if (stat /= 0) call cq_abort('phi_val: failed to deallocate arrays', stat)
    call reg_dealloc_mem(area_ops, 7 * n, type_dbl)

  end function phi_val
  !!*****


  !!****f* vdWDFT_module/pofq
  !! PURPOSE
  !!   Finds the values and derivativs, at q0, of the cubic polynomials
  !!   p_1(q0) such that
  !!            y(q0) = \sum_i p_i(q0) * y_i
  !!   is the cubic spline interpolation at q0 of (any) function y(q) with
  !!   values y_i at mesh points qmesh_i
  !! USAGE
  !!   call pofq(q0, p0, dp0,dq0)
  !! INPUTS
  !!   real(double) q0 : value of q0
  !! OUTPUT
  !!   real(double) p0(mq)     : values of the cubic polynomials on qmesh
  !!   real(double) dp0dq0(mq) : values of the derivative of the polynomials
  !!                             on qmesh
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/02/29
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine pofq(q0, p0, dp0dq0)

    use numbers
    use splines, only: spline_nonU

    implicit none
    ! input parameters
    real(double),                intent(in)  :: q0
    real(double), dimension(mq), intent(out) :: p0
    real(double), dimension(mq), intent(out) :: dp0dq0
    ! local variables
    integer       :: iq, iq0
    real(double)  :: a, b, dq
    logical, save :: first_call = .true.
    real(double), dimension(mq,mq), save :: p, d2pdq2

    ! Set up spline polynomial basis
    if (first_call) then
       p = zero
       do iq = 1, mq
          p(iq,iq) = one
          call spline_nonU(mq, qmesh, p(:,iq), huge(one), huge(one), &
                           d2pdq2(:,iq))
       end do
       first_call = .false.
    end if

    ! Find interval of qmesh in which q0 is included
    if (q0 > qmesh(mq)) then   ! q0 out of range
       p0 = zero
       dp0dq0 = zero
       return
    end if
    iq0 = iofq(q0)

    ! Evaluate polynomials of spline basis
    dq = qmesh(iq0+1) - qmesh(iq0)
    a = (qmesh(iq0+1) - q0) / dq   ! dadq0 = -1/dq
    b = (q0 - qmesh(iq0)) / dq     ! dbdq0 = +1/dq
    do iq = 1, mq
       p0(iq) = a*p(iq0,iq) + b*p(iq0+1,iq) + &
                ((a**3-a)*d2pdq2(iq0,iq) + (b**3-b)*d2pdq2(iq0+1,iq)) * dq**2/6
       dp0dq0(iq) = - (p(iq0,iq) - p(iq0+1,iq)) / dq - &
                   ((3*a**2-1)*d2pdq2(iq0,iq) - (3*b**2-1)*d2pdq2(iq0+1,iq)) * dq/6
    end do

  end subroutine pofq
  !!*****


  !!****f* vdWDFT_module/qofrho
  !! PURPOSE
  !!   Finds the local wavevector parameter q0 defined in Eqs.(11-12) of
  !!   Dion et al. PRL 92, 246401 (2004)
  !! USAGE
  !!   call qofrho(rho_r, grho_r, q_r)
  !! INPUTS
  !!   real(double) rho_r(nspin)     : electron density at r
  !!   real(double) grho_r(3,nspin)  : gradient of electron density at r
  !! OUTPUT
  !!   real(double) q_r              : wave vector parameter q0 at r
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !!   2018/02/13 12:23 dave
  !!    Changed call for Vxc to XC_CQ
  !! SOURCE
  !!
  subroutine qofrho(rho_r, grho_r, q_r)

    use datatypes
    use numbers
    use GenComms,      only: cq_abort
    use global_module, only: nspin
    use XC,     only: Vxc_of_r_LSDA_PW92

    implicit none

    ! input parameters
    real(double), dimension(:),     intent(in)  :: rho_r
    real(double), dimension(:,:),   intent(in)  :: grho_r
    real(double),                   intent(out) :: q_r
    ! local variables
    real(double), parameter    :: zab = -0.8491_double ! See Dion et al
    real(double), dimension(3) :: grho_tot_r
    real(double) :: rho_tot_r, grho2_r, q0_r, kf_r, ex_r, ec_r

    ! Initialise q0
    q_r = zero

    ! Get total density and gradient of density
    rho_tot_r = rho_r(1) + rho_r(nspin)
    grho_tot_r(1:3) = grho_r(1:3,1) + grho_r(1:3,nspin)

    ! loop over grid points
    if (rho_tot_r < RD_ERR) then
       q_r = qcut
    else
       ! Calculate Kf(r)
       kf_r = (three*pi**2 * rho_tot_r)**third

       ! Find exchange and correlation energy densities
       call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_x=ex_r, eps_c=ec_r)

       ! Find q0
       grho2_r = sum(grho_tot_r**2)

       if (abs(ex_r) < RD_ERR) then
          q0_r = qcut
       else
          q0_r = &
               (one + ec_r / ex_r - &
                (zab / nine) * grho2_r / (two * kf_r * rho_tot_r)**2) * kf_r
       end if

       ! Saturate q0 to qcut smoothly
       call saturate(q0_r, qcut, q_r)
    end if

  end subroutine qofrho
  !!*****


  !!****f* vdWDFT_module/saturate
  !! PURPOSE
  !!   Defines a function, given xc (saturation value/cutoff)
  !!              y(x) = xc * (1 - exp(-sum_{n=1}^nsat (x / xc)**n / n))
  !!   It is approx. equal to x for x < xc and it saturates to xc when x->infinity
  !! USAGE
  !!   call saturate(x, xc, y[, dydx])
  !! INPUTS
  !!   real(double)           x    : independent variable
  !!   real(double)           xc   : saturation value
  !! OUTPUT
  !!   real(double)           y    : function value
  !!   real(double), optional dydx : derivative dy/dx
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/03
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine saturate(x, xc, y, dydx)

    use datatypes
    use numbers

    implicit none

    ! input parameters
    real(double), intent(in)  :: x
    real(double), intent(in)  :: xc
    real(double), intent(out) :: y
    real(double), optional, intent(out) :: dydx
    ! local variables
    integer :: n
    real(double):: dpdx, p

    if (nsat >= 100) then
       if (x < xc) then
          y = x
          if (present(dydx)) dydx = one
       else ! (x >= xc)
          y = xc
          if (present(dydx)) dydx = zero
       end if ! (x < xc)
    else ! (nsat < 100)
       !  This is the straightforward polynomial evaluation
       !  p = x/xc
       !  dpdx = 1/xc
       !  do n = 2,nsat
       !    p = p + (x/xc)**n / n
       !    dpdx = dpdx + (x/xc)**(n-1) / xc
       !  end do
       ! And this is a more accurate way to evaluate it
       p = (x / xc) / real(nsat,double)
       dpdx = one / xc
       do n = nsat - 1, 1, -1
          p = (p + one / real(n, double)) * x / xc
          if (present(dydx)) dpdx = (dpdx * x + one) / xc
       end do
       y = xc * (one - exp(-p))
       if (present(dydx)) dydx = xc * dpdx * exp(-p)
    end if ! (nsat >= 100)

  end subroutine saturate
  !!*****


  !!****f* vdWDFT_module/saturate_inverse
  !! PURPOSE
  !!   Finds the inverse of the function defined in saturate subroutine
  !! USAGE
  !!   call saturate_inverse(y, xc, x, dydx)
  !! INPUTS
  !!   real(double) y  : Independent variable
  !!   real(double) xc : Saturation value
  !! OUTPUT
  !!   real(double) x  : Inverse function value
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/03
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine saturate_inverse(y, xc, x)

    use datatypes
    use numbers
    use GenComms, only: cq_abort

    implicit none

    ! passed parameters
    real(double), intent(in) :: y
    real(double), intent(in) :: xc
    real(double), intent(out):: x
    ! local variables
    real(double):: x1, x2, yx

    if (y < zero .or. y > xc) &
         call cq_abort('saturate_inverse: y out of range')
    x1 = zero
    x2 = xmaxbyxc * xc
    do
       x = (x1 + x2) / two
       call saturate(x, xc, yx)
       if (abs(y - yx) < ytol) then
          return
       else if (yx < y) then
          x1 = x
       else
          x2 = x
       end if
    end do

  end subroutine saturate_inverse
  !!*****


  !!****f* vdWDFT_module/set_phi_table
  !! PURPOSE
  !!   Finds and writes in disk the interpolation table (mesh points and
  !!   function values) for the kernel phi(d1,d2). If the table file
  !!   already exists, it is simply read and stored in memory.
  !! USAGE
  !!   call set_phi_table()
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/03
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine set_phi_table()

    use datatypes
    use numbers
    use global_module,  only: io_lun
    use GenComms,       only: cq_abort, gcopy, inode, ionode
    use input_module,   only: io_assign, io_close
    use vdWMesh_module, only: set_mesh, get_mesh, set_interpolation, derivative
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! local variables
    logical      :: file_found
    integer      :: lun, stat
    integer      :: id, id1, id2, nmesh
    real(double) :: d, d1, d1m, d1p, d2, d2m, d2p, dd, phi1, phim, phip,&
                    phimm, phimp, phipm, phipp
    real(double), dimension(:,:), allocatable :: dphidd1, dphidd2, &
                                                 d2phidd1dd2, phi

    allocate(dphidd1(nd,nd), dphidd2(nd,nd), d2phidd1dd2(nd,nd), &
             phi(nd,nd), STAT=stat)
    if (stat /= 0) call cq_abort("set_phi_table: Error alloc mem: ", nd)
    call reg_alloc_mem(area_ops, nd*nd, type_dbl)

    ! Read file with table, if present
    if (inode == ionode) then
       call io_assign(lun)
       inquire(file='vdW_kernel.tab', exist=file_found)
       if (file_found) then
          open(unit=lun, file='vdW_kernel.tab', form='unformatted')
          read (lun, end=1) nmesh
          if (nmesh == nd) then
             read (lun, end=1) dmesh
             read (lun, end=1) phi_table
          end if
          call io_close(lun)
       end if
1      continue ! Come here if end-of-file found
    end if
    ! copy to other nodes
    call gcopy(file_found)
    if (file_found) then
       call gcopy(nmesh)
       call gcopy(dmesh, nd)
       call gcopy(phi_table, 4, 4, nd, nd)
       phi_table_set = .true.
       deallocate(dphidd1, dphidd2, d2phidd1dd2, phi, STAT=stat)
       if (stat /= 0) call cq_abort("set_phi_table: Error dealloc mem")
       call reg_dealloc_mem(area_ops, nd*nd, type_dbl)
       return
    end if

    ! Set d-mesh points
    call set_mesh(nd, xmax=dcut, dxndx1=ddmaxddmin)
    call get_mesh(nd, nmesh, dmesh)

    ! Find function at mesh points
    do id1 = 1, nd
       d1 = dmesh(id1)
       phi1 = phi_soft(d1, d1)
       phi(id1,id1) = phi1
       do id2 = 1, id1 - 1
          d2 = dmesh(id2)
          d = sqrt(d1**2 + d2**2)
          if (d < dsoft) then
             phi(id1,id2) = phi_soft( d1, d2 )
          else
             if (use_dphi) then ! Use dphi for better efficiency
                phi(id1,id2) = phi1 + dphi(d1, d2)
             else ! Use only phi_val, to eliminate uncertainties
                phi(id1,id2) = phi_val(d1, d2)
             end if
          end if
          phi(id2,id1) = phi(id1,id2)
       end do ! id2
    end do ! id1

    if (deriv_method == 'numeric') then
       if (inode == ionode) &
            write (io_lun, '(8x,"set_phi_table: Using numerical derivatives")')

       ! Find derivatives at mesh points
       do id1 = 1, nd
          d1 = dmesh(id1)
          dd = ddbyd * d1
          dd = max(dd, ddmin)
          ! d1 = max(d1, dd)
          d1m = d1 - dd
          d1p = d1 + dd
          phim = phi_soft(d1m, d1m)
          phip = phi_soft(d1p, d1p)
          do id2 = 1,id1
             d2  = dmesh(id2)
             ! d2 = max(d2, dd)
             d = sqrt(d1**2 + d2**2)
             d1m = d1 - dd
             d1p = d1 + dd
             d2m = d2 - dd
             d2p = d2 + dd
             if (d < dsoft) then
                phimm = phi_soft(d1m, d2m)
                phipm = phi_soft(d1p, d2m)
                phipp = phi_soft(d1p, d2p)
                phimp = phi_soft(d1m, d2p)
             else ! (d>dsoft)
                if (use_dphi) then
                   phimm = phim + dphi(d1m, d2m)
                   phipm = phip + dphi(d1p, d2m)
                   phipp = phip + dphi(d1p, d2p)
                   if (id1 == id2) then
                      phimp = phip + dphi(d1m, d2p)
                   else
                      phimp = phim + dphi(d1m, d2p)
                   end if
                else ! (.not.use_dphi)
                   phimm = phi_val(d1m, d2m)
                   phipm = phi_val(d1p, d2m)
                   phipp = phi_val(d1p, d2p)
                   phimp = phi_val(d1m, d2p)
                end if ! (use_dphi)
             end if ! (d<dsoft)

             dphidd1(id1,id2)     = (phipp + phipm - phimp - phimm) / (4*dd)
             dphidd2(id1,id2)     = (phipp - phipm + phimp - phimm) / (4*dd)
             d2phidd1dd2(id1,id2) = (phipp - phipm - phimp + phimm) / (2*dd)**2

             dphidd1(id2,id1)     = dphidd2(id1,id2)
             dphidd2(id2,id1)     = dphidd1(id1,id2)
             d2phidd1dd2(id2,id1) = d2phidd1dd2(id1,id2)

          end do ! id2
       end do ! id1

    else if (deriv_method == 'interp') then
       if (inode == ionode) &
            write (io_lun, '("set_phi_table: Using interpolation for &
                              &derivatives")')

       ! Reset mesh, which has been changed by phi_val
       call set_mesh(nd, xmax=dcut, dxndx1=ddmaxddmin)

       ! Set interpolation method
       if (interp_method == 'Lagrange') then
          call set_interpolation('Lagrange')
       else if (interp_method == 'Spline') then
          !      call set_interpolation('Spline', huge(phi), huge(phi))
          call set_interpolation('Spline', zero, zero)
       else
          call cq_abort('set_phi_val: ERROR: Unknown interp_method')
       end if

       ! Find first partial derivatives d_phi/d_d1 and d_phi/d_d2
       do id = 1, nd
          dphidd1(:,id) = derivative(nd, phi(:,id))
          ! using the fact phi is symmetric in d1 and d2
          dphidd2(id,:) = dphidd1(:,id)
       end do

       ! Find second cross partial derivative d_phi/d_d1/d_d2
       do id = 1, nd
          d2phidd1dd2(id,:) = derivative(nd, dphidd1(id,:))
       end do

       ! Symmetrize d_phi/d_d1/d_d2
       do id2 = 2, nd
          do id1 = 1, id2 - 1
             d2phidd1dd2(id1,id2) = &
                  (d2phidd1dd2(id1,id2) + d2phidd1dd2(id2,id1)) / two
             d2phidd1dd2(id2,id1) = d2phidd1dd2(id1,id2)
          end do
       end do

    else

       call cq_abort('set_phi_table: ERROR: Unknown deriv_method')

    end if ! (deriv_method)

    ! Make values and derivatives strictly zero when d1=dmax or d2=dmax
    phi(:,nd) = zero
    phi(nd,:) = zero
    dphidd1(:,nd) = zero
    dphidd1(nd,:) = zero
    dphidd2(:,nd) = zero
    dphidd2(nd,:) = zero
    d2phidd1dd2(:,nd) = zero
    d2phidd1dd2(nd,:) = zero

    ! Make dphi(d1,d2)/dd1=0 for d1=0 and dphi(d1,d2)/dd2=0 for d2=0
    dphidd1(1,:) = zero
    dphidd2(:,1) = zero
    d2phidd1dd2(:,1) = zero
    d2phidd1dd2(1,:) = zero

    ! Set up bicubic interpolation coefficients
    call find_bicubic_coeffs(nd, nd, dmesh, dmesh, phi, dphidd1, dphidd2, &
                d2phidd1dd2, phi_table)

    ! Save phi_table in file, note the same  mesh and table are on all nodes.
    if (inode == ionode) then
       call io_assign(lun)
       open(unit=lun, file='vdW_kernel.tab', form='unformatted')
       write(lun) nd
       write(lun) dmesh
       write(lun) phi_table
       call io_close(lun)
    end if

    ! Mark table as set
    phi_table_set = .true.

    deallocate(dphidd1, dphidd2, d2phidd1dd2, phi, STAT=stat)
    if (stat /= 0) call cq_abort("set_phi_table: Error dealloc mem")
    call reg_dealloc_mem(area_ops, nd*nd, type_dbl)

  end subroutine set_phi_table
  !!*****


  !!****f* vdWDFT_module/set_qmesh
  !! PURPOSE
  !!   Sets mesh of q values
  !! USAGE
  !!   call set_qmesh()
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/03
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine set_qmesh()

    use datatypes
    use vdWMesh_module, only: set_mesh, get_mesh

    implicit none

    ! local variables
    integer :: nmesh

    call set_mesh(mq, xmax=qcut, dxndx1=dqmaxdqmin)
    call get_mesh(mq, nmesh, qmesh)
    qmesh_set = .true.

  end subroutine set_qmesh
  !!*****


  !!****f* vdWDFT_module/vdW_decusp
  !! PURPOSE
  !!   Finds the local energy correction due to the softening of the
  !!   VdW kernel cusp, defined as
  !!
  !!        eps(rho,grad_rho) = (1/2) * rho * Int 4*pi*r**2*dr *
  !!                            (phi(q1*r,q2*r) - phi_soft(q1*r,q2*r))
  !!
  !!   where q1=q2=q0(rho,grad_rho). Notice that grad_rho is included
  !!   in the value of q0 at the origin but not in the change of
  !!   rho(r) in the integrand.  phi_soft(d1,d2) is a softened version
  !!   of the nonlocal VdW kernel phi(d1,d2) (Eq.(14) of Dion et al)
  !!   in which the logarithmic divergence at d1=d2=0 is substituted
  !!   by a smooth analytic function of the form defined in vdW_phi
  !! USAGE
  !!   call vdW_decusp(rhos, grhos, eps)
  !! INPUTS
  !!   real(double) rho_r(nspin)    : Electron spin density at r
  !!   real(double) grho_r(3,nspin) : Spin density gradient at r
  !! OUTPUT
  !!   real(double) eps_r           : Energy correction density at r
  !! NOTES
  !! - Requires a previous call to vdW_set_kcut. Otherwise stops with
  !!   an error msg.
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/04
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine vdW_decusp(rho_r, grho_r, eps_r)

    use datatypes
    use numbers
    use GenComms,        only: cq_abort
    use global_module,   only: nspin
    use memory_module,   only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    real(double), dimension(:),   intent(in)  :: rho_r
    real(double), dimension(:,:), intent(in)  :: grho_r
    real(double),                 intent(out) :: eps_r

    logical, save :: initialized = .false.
    integer       :: iq1, iq2, ir, stat
    real(double)  :: ptp, phi22, r, rho_tot_r, q_r
    real(double), dimension(mq,mq), save :: table
    real(double), dimension(:),   allocatable :: p, pt, dpdq
    real(double), dimension(:,:), allocatable :: phi

    allocate(p(mq), pt(mq), dpdq(mq), phi(mq,mq), STAT=stat)
    if (stat /= 0) call cq_abort("vdW_decusp: Error alloc mem: ", mq)
    call reg_alloc_mem(area_ops, (3 + mq)*mq, type_dbl)

    if (.not. initialized) then

       if (.not. phi_table_set) call set_phi_table()
       if (.not. kcut_set) &
            call cq_abort('vdW_decusp: ERROR: kcut has not been set')

       table = zero
       do ir = 1, nr
          ! looping over radial mesh
          r = (ir - 1) * dr
          call phiofr(r, phi)
          do iq2 = 1, mq
             do iq1 = 1, iq2
                ! d1 = qmesh(iq1) * r
                ! d2 = qmesh(iq2) * r
                ! d = sqrt(d1**2 + d2**2)
                ! if (d < dsoft) then
                !    table(iq1,iq2) = table(iq1,iq2) + &
                !                     twopi*dr * r**2 * &
                !                     (phi_val(d1, d2) - phi(iq1,iq2))
                ! end if
                table(iq1,iq2) = table(iq1,iq2) - twopi*dr * r**2 * phi(iq1,iq2)
             end do ! iq1
          end do ! iq2
       end do

       do iq2 = 2, mq
          do iq1 = 1, iq2 - 1
             table(iq2,iq1) = table(iq1,iq2)
          end do
       end do

       initialized = .true.

    end if ! (.not.initialized)

    call qofrho(rho_r, grho_r, q_r)
    call pofq(q_r, p, dpdq)
    pt = matmul(p, table)
    ptp = sum(pt * p)
    rho_tot_r = rho_r(1) + rho_r(nspin)
    eps_r = rho_tot_r * ptp

    deallocate(p, pt, dpdq, phi, STAT=stat)
    if (stat /= 0) call cq_abort("vdW_decusp: Error dealloc mem")
    call reg_dealloc_mem(area_ops, (3 + mq)*mq, type_dbl)

  end subroutine vdW_decusp
  !!*****


  !!****f* vdWDFT_module/vdW_get_qmesh
  !! PURPOSE
  !!   Returns size and values of q-mesh
  !! USAGE
  !! OUTPUT
  !!   integer,                n    : Number of q mesh points
  !!   real(double), optional  q(:) : Values of q mesh points
  !! EXAMPLE
  !!   integer                   :: nq
  !!   real(double)              :: kcut
  !!   real(double), allocatable :: qmesh(:)
  !!   kcut = 10.0_double ! 10 Bohr^-1 => 100 Ry (this should be
  !!                                              adapted to your mesh)
  !!   call vdW_set_kcut(kcut)
  !!   call vdW_get_qmesh(nq)
  !!   allocate(qmesh(nq))
  !!   call vdW_get_qmesh(nq, qmesh)
  !! NOTES
  !! - Requires a previous call to vdW_set_kcut. Otherwise stops with
  !!   an error msg.
  !! - If the size of array q is smaller than that of the stored
  !!   qmesh, it is filled with the first size(q) values of qmesh
  !! - The size and values of the logarithmic q mesh are set by
  !!   internal parameters that can be changed only by editing them
  !!   in this module:
  !!     nq             : Number of q mesh points
  !!     qcut=qmesh(nq) : Max. value of q mesh
  !!     dqmaxdqmin     : (qmesh(nq) - qmesh(nq-1)) / (qmesh(2) - qmesh(1))
  !!   Although the presently-set values have been found to yield good
  !!   accuracy in preliminary calculations, more tests may be
  !!   required to guarantee convergence in other systems. The value
  !!   of nq is particularly important: larger nq increases accuracy
  !!   but CPU time increases between nq and nq**2.
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/05
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine vdW_get_qmesh(n, q)

    use datatypes

    implicit none

    integer,                intent(out) :: n
    real(double), optional, intent(out) :: q(:)
    integer:: nmax

    if (.not. qmesh_set) call set_qmesh()
    n = nq
    if (present(q)) then
       nmax = max(nq, size(q))
       q(1:nmax) = qmesh(1:nmax)
    end if
  end subroutine vdW_get_qmesh
  !!*****


  !!****f* vdWDFT_module/vdW_phi
  !! PURPOSE
  !!   Finds phi_soft(q1,q2,k) (Fourier transform of
  !!   phi_soft(q1,q2,r)) for all values of q1 and q2 in qmesh. The
  !!   q's are local wavevectors defined in eqs.(11-12), and phi_soft
  !!   is a smoothed version of Eq.(14) of Dion et al, defined as
  !!   phi_soft(d1,d2) = phi_soft(d,a) = phi0 + phi2*d**2 + phi4*d**4,
  !!   where phi0 is a parameter, d=sqrt(d1**2+d2**2), a=atan(d2/d1),
  !!   and phi2, phi4 are chosen so that phi_soft(d,a) matches
  !!   phi(d,a) in value and slope at d=dsoft (another parameter).
  !! USAGE
  !! INPUTS
  !!   real(double)  k           : Modulus of actual k vector
  !! OUTPUT
  !!   real(double)  phi(:,:)    : phi(q1,q2,k) at given k
  !!                               for all q1,q2 in qmesh
  !!   real(double)  dphidk(:,:) : dphi(q1,q2,k)/dk at given k
  !! EXAMPLE
  !!   integer :: nq
  !!   real(double):: k, kcut
  !!   real(double), allocatable :: phi(:,:), dphidk(:,:)
  !!   kcut = 10.0_double ! 10 Bohr^-1 => 100 Ry (this should be
  !!                                              adapted to your mesh)
  !!   call vdW_set_kcut(kcut)
  !!   call vdW_get_qmesh(nq)
  !!   allocate(phi(nq,nq), dphidk(nq,nq))
  !!   do k points
  !!     call vdW_phi(k, phi, dphidk)
  !! NOTES
  !! - Requires a previous call to vdW_set_kcut. Otherwise stops with
  !!   an error msg.
  !! - Stops with an error message if size of array phi is smaller
  !!   than nq*nq.
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/05
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine vdW_phi(k, phi, dphidk)

    use datatypes
    use numbers
    use GenComms, only: cq_abort

    implicit none

    ! parameters
    real(double),                           intent(in)  :: k
    real(double), dimension(:,:),           intent(out) :: phi
    real(double), dimension(:,:), optional, intent(out) :: dphidk
    ! local variables
    integer      :: ik_lo, ik_hi, iq1, iq2
    real(double) :: k_lo, k_hi, a, b, c, d, dc, dd

    if (.not. kcut_set) &
         call cq_abort('vdW_phi: ERROR: kcut must be previously set')

    ! Check argument sizes
    if (size(phi,1) < nq .or. size(phi,2) < nq) &
         call cq_abort('vdW_phi: ERROR: size(phi) too small')

    ! Find phi values at point k
    if (k >= kcut) then
       phi(:,:) = zero
    else
       ! Expand interpolation inline since this is the hottest point in VdW
       ! See Numerical Recipes for the spline interpolation
       ! implementation use here. Note arrays has dimension (1:nk)
       ik_lo = aint(k / dk) + 1 ! integer division gives lower bound,
                                ! add 1 for correct array index
       ik_hi = ik_lo + 1
       k_lo = dk * real(ik_lo - 1, double)
       k_hi = dk * real(ik_hi - 1, double)
       a = (k_hi - k) / dk
       b = (k - k_lo) / dk
       if (present(dphidk)) then
          dc = - (three * a*a - one) * dk / six
          dd = (three * b*b - one) * dk / six
       end if
       c = (a*a*a - a) * dk*dk / six
       d = (b*b*b - b) * dk*dk / six
       do iq2 = 1, nq
          do iq1 = 1, iq2
             phi(iq1,iq2) = &
                  a * phik(ik_lo,iq1,iq2) + b * phik(ik_hi,iq1,iq2) + &
                  c * d2phidk2(ik_lo,iq1,iq2) + d * d2phidk2(ik_hi,iq1,iq2)
             phi(iq2,iq1) = phi(iq1,iq2)
             if (present(dphidk)) then
                dphidk(iq1,iq2) = &
                     (-phik(ik_lo,iq1,iq2) + phik(ik_hi,iq1,iq2)) / dk + &
                     dc * d2phidk2(ik_lo,iq1,iq2) + dd * d2phidk2(ik_hi,iq1,iq2)
                dphidk(iq2,iq1) = dphidk(iq1,iq2)
             end if
          end do
       end do
    end if

  end subroutine vdW_phi
  !!*****


  ! function vdW_phi_iqjq(k, iq1, iq2)

  !   use datatypes
  !   use numbers
  !   use GenComms, only: cq_abort

  !   implicit none

  !   ! parameters
  !   real(double), intent(in)  :: k
  !   integer,      intent(in)  :: iq1, iq2
  !   ! returned value
  !   real(double) :: vdW_phi_iqjq

  !   if (k >= kcut) then
  !      phi = zero
  !   else
  !      ! Expand interpolation inline since this is the hottest point in VdW
  !      ik = int(k / dk)
  !      a = (real(ik + 1, double) * dk - k) / dk
  !      b = one - a
  !      a2 = (three * a**2 - one) * dk / six
  !      b2 = (three * b**2 - one) * dk / six
  !      a3 = (three**3 - a) * dk**2 / six
  !      b3 = (b**3 - b) * dk**2 / six
  !      vdW_phi_iqjq = &
  !           a * phik(ik,iq1,iq2) + b * phik(ik+1,iq1,iq2) + &
  !           a3 * d2phidk2(ik,iq1,iq2) + b3 * d2phidk2(ik+1,iq1,iq2)
  !   end if

  ! end function vdW_phi_iqjq


  !!****f* vdWDFT_module/vdW_set_kcut
  !! PURPOSE
  !!   Sets the reciprocal planewave cutoff kc of the integration grid,
  !!   and finds the interpolation table to be used by vdW_phi to obtain
  !!   the vdW kernel phi at the reciprocal mesh wavevectors.
  !! USAGE
  !!   call vdW_set_kcut(kc)
  !! INPUTS
  !!   real(double) kc : Planewave cutoff: k>kcut => phi=0
  !! NOTES
  !! - An interpolation table to calculate the VdW kernel phi is read
  !!   from file 'vdW_kernel.tab'. If the file does not exist, the
  !!   table is calculated and the file written when vdW_set_kcut is
  !!   called.
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/05
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine vdW_set_kcut(kc)

    use datatypes
    use numbers
    use splines, only: spline
    use GenComms,      only: cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! input parameters
    real(double), intent(in) :: kc
    ! local variables
    integer      :: ik, iq1, iq2, ir, nrs, stat
    real(double) :: dphids, dphidk0, dphidkmax, dphidr0, dphidrmax, &
                    k, phi0, phi2, phis, q1, q2, rs
    real(double), dimension(:), allocatable :: phi, r

    allocate(phi(1:nr), r(1:nr), STAT=stat)
    if (stat /= 0) call cq_abort("vdW_set_kcut: Error alloc mem: ", nr)
    call reg_alloc_mem(area_ops, 2*nr, type_dbl)

    if (kcut_set) return   ! Work alredy done
    if (.not. qmesh_set) call set_qmesh()

    ! This is the radial grid for phi's, (unrelated to the main 3D FFT grid)
    dr = rcut / nr
    dk = pi / rcut
    ! work out number of array elements just covering kc (i.e. kcut)
    nk = aint(kc/dk) + 1
    nrs = nint(rsoft/dr)
    rs = real(nrs,double) * dr

    ! For each pair of values q1 and q2
    do iq2 = 1, mq
       do iq1 = 1, iq2

          ! Saturated q values
          ! q1 = qmesh(iq1)
          ! q2 = qmesh(iq2)

          ! Find original (unsaturated) q values
          call saturate_inverse(qmesh(iq1), qcut, q1)
          call saturate_inverse(qmesh(iq2), qcut, q2)

          ! Find kernel in real space
          do ir = 1, nr
             r(ir) = real(ir - 1, double) * dr
             phir(ir,iq1,iq2) = phi_interp(q1*r(ir), q2*r(ir))
          end do
          phi(:) = phir(:,iq1,iq2)

          ! Change kernel near origin to a parabola matching at rs
          if (nrs > 0) then
             phis = phi(nrs)
             dphids = (phi(nrs+1) - phi(nrs-1)) / (2*dr)
             phi0 = phis - dphids * rs/2
             phi2 = dphids / (2*rs)
             do ir = 1, nrs
                phir(ir,iq1,iq2) = phi0 + phi2 * r(ir)**2
             end do
          end if ! (nrs>0)

          ! Kill kernel smoothly at r=rcut
          do ir = 1, nr
             phir(ir,iq1,iq2) = phir(ir,iq1,iq2) * cutoff(r(ir)/rcut)
          end do

          ! Optimized filter (warning: inaccurate for very large kcut*rcut)
          !      call filter( 0, nr+1, r(:), phir(:,iq1,iq2), kc, 0 )

          ! Find kernel in reciprocal space
          call radfft(phir(:,iq1,iq2), phik(:,iq1,iq2), nr, dr)

          ! Filter out above kcut (kcut is for the main 3D FFT grid,
          ! this is to ensure the maximum frequency of the radial grid
          ! do not exceed the maximum allowed for the 3D grid)
          phik(nk+1:nr,iq1,iq2) = zero

          ! Soft filter below kcut
          do ik = 1, nk
             k = real(ik - 1, double) * dk
             phik(ik,iq1,iq2) = phik(ik,iq1,iq2) * cutoff(k/kc)
          end do

          ! Find filtered kernel in real space
          ! note that for the radfft implementation k =
          ! reciprocal_vector/(2*pi), so need to divide kmax by
          ! factor of twopi to get the correct dk.
          call radfft(phik(:,iq1,iq2), phir(:,iq1,iq2), nr, dk/twopi)

          ! Set up spline interpolation tables
          dphidr0 = zero
          dphidrmax = zero
          dphidk0 = zero
          dphidkmax = zero
          call spline(nr, dr, phir(:,iq1,iq2), dphidr0, dphidrmax, &
                      d2phidr2(:,iq1,iq2))
          call spline(nk, dk, phik(:,iq1,iq2), dphidk0, dphidkmax, &
                      d2phidk2(:,iq1,iq2))

          ! Fill symmetric elements
          phir(:,iq2,iq1) = phir(:,iq1,iq2)
          phik(:,iq2,iq1) = phik(:,iq1,iq2)
          d2phidr2(:,iq2,iq1) = d2phidr2(:,iq1,iq2)
          d2phidk2(:,iq2,iq1) = d2phidk2(:,iq1,iq2)

          !      if (.false. .and. iq1==iq2) then
          !        print*, 'vdW_set_kcut: iq1,iq2=', iq1, iq2
          !        call window( 0._dp, 5._dp, -1._dp, 4._dp, 0 )
          !        call axes( 0._dp, 1._dp, 0._dp, 1._dp )
          !        call plot( nr+1, r, phi, phir(:,iq1,iq2) )
          !        call window( 0._dp, 10._dp, -0.05_dp, 0.15_dp, 0 )
          !        call axes( 0._dp, 1._dp, 0._dp, 0.05_dp )
          !        call plot( nr+1, r, q1*q2*r**2*phi, q1*q2*r**2*phir(:,iq1,iq2) )
          !        call show()
          !      end if

       end do ! iq1
    end do ! iq2

    !  print'(a,/,(2i6,f12.6))', 'vdW_set_kcut: iq1, iq2, phir(0,iq1,iq2) =', &
    !    ((iq1,iq2,phir(0,iq1,iq2),iq1=2,iq2),iq2=2,mq)

    kcut = kc
    kcut_set = .true.

    deallocate(phi, r, STAT=stat)
    if (stat /= 0) call cq_abort("vdW_set_kcut: Error dealloc mem")
    call reg_dealloc_mem(area_ops, 2*nr, type_dbl)

  end subroutine vdW_set_kcut
  !!*****


  !!****f* vdWDFT_module/vdW_theta
  !! PURPOSE
  !!   Finds the value and derivatives of theta_i(rho,grad_rho) =
  !!   rho*p_i(q0), where q0(rho,grad_rho) is the local wavevector
  !!   defined in eqs.(11-12) of Dion et al, PRL 92, 246401 (2004).
  !!   p_i(q0) are the cubic polynomials such that
  !!     y(q0) = Sum_i p_i(q0) * y_i
  !!   is the cubic spline interpolation at q0 of (any) function y(q)
  !!   with values y_i at mesh points qmesh_i
  !! USAGE
  !!   call vdW_theta(rho_r, grho_r, theta_r)
  !! INPUTS
  !!   real(double) rho_r(nspin)    : electron spin density at r
  !!   real(double) grho_r(3,nspin) : spin density gradient at r
  !! OUTPUT
  !!   real(double) theta_r(nq)     : expansion of rho_r * q(r) in qmesh
  !! NOTES
  !! - Requires a previous call to vdW_set_kcut
  !! - The value of q0(rho,grad_rho) is saturated at qmesh(nq) = qcut
  !!   (a parameter) to avoid incorrect 'interpolation' of
  !!   phi(q1,q2,r12) from qmesh points.
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/08
  !! MODIFICATION HISTORY
  !!   2013/07/10 11:43 dave
  !!    Bug fix for sum over two components of rho even without spin
  !! SOURCE
  !!
  subroutine vdW_theta(rho_r, grho_r, theta_r)

    use datatypes
    use numbers
    use global_module, only: nspin, spin_factor
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use GenComms,      only: cq_abort

    implicit none

    ! passed parameters
    real(double), dimension(:),   intent(in)  :: rho_r
    real(double), dimension(:,:), intent(in)  :: grho_r
    real(double), dimension(:),   intent(out) :: theta_r
    ! local variables
    integer :: stat
    real(double) :: rho_tot_r, q_r
    real(double), dimension(3) :: grho_tot_r
    real(double), dimension(:), allocatable :: dpdq, p

    allocate(dpdq(mq), p(mq), STAT=stat)
    if (stat /= 0) call cq_abort("vdW_theta: Error alloc mem: ", mq)
    call reg_alloc_mem(area_ops, 2*mq, type_dbl)

    theta_r = zero

    rho_tot_r = rho_r(1) + rho_r(nspin)
    if(nspin==1) then
       grho_tot_r(:) = spin_factor * grho_r(:,1)
    else
       grho_tot_r = spin_factor * sum(grho_r, nspin)
    end if
    call qofrho(rho_r, grho_r, q_r)
    call pofq(q_r, p, dpdq)
    theta_r(1:nq) = rho_tot_r * p(1:nq)

    deallocate(dpdq, p, STAT=stat)
    if (stat /= 0) call cq_abort("vdW_theta: Error dealloc mem")
    call reg_dealloc_mem(area_ops, 2*mq, type_dbl)

  end subroutine vdW_theta
  !!*****


  !!****f* vdWDFT_module/vdWXC_energy
  !! PURPOSE
  !!   Calculates the van der Waals xc energy given by Dion et al. PRL
  !!   92, 246401 (2004) and made to order N log2 N using the method
  !!   described by Roman-Perez et al. PRL 103, 096102
  !! USAGE
  !!   call vdWXC_energy(rho, E_vdWXC)
  !! INPUTS
  !!   real(double) rho(maxngrid,nspin) : electron densities in each spin channel
  !! OUTPUT
  !!   real(double) E_vdWXC             : vdW xc energy
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/20
  !! MODIFICATION HISTORY
  !!   2012/05/29 L.Tong
  !!   - Added timers
  !! SOURCE
  !!
  subroutine vdWXC_energy(rho, E_vdWXC)

    use datatypes
    use numbers
    use dimens
    use global_module, only: nspin, spin_factor,       &
                             IPRINT_TIME_THRES2
    use maxima_module, only: maxngrid
    use XC,     only: build_gradient,           &
                             eps_xc_of_r_GGA_PBE,      &
                             Vxc_of_r_LSDA_PW92, functional_gga_pbe96_rev98
    use fft_module,    only: fft3, z_columns_node, recip_vector
    use GenComms,      only: gsum, inode, ionode, cq_abort
    use timer_module
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

! LT_debug 2012/04/30 begin
!    use units
! LT_debug 2012/04/30 end

    implicit none

    ! passed parameters
    real(double), dimension(:,:), intent(in)  :: rho
    real(double),                 intent(out) :: E_vdWXC
    ! lcoal variables
    integer :: n_q, ix, iq, jq, ir, ik, spin, stat
    real(double) :: eps_x, eps_c
    real(double) :: mag_k, phi_k_iq_jq, tk_re, tk_im, eps_c_usp_r,    &
                    rho_tot_r, eps_vdW_nl_r, eps_gga_x, eps_lda_c
    real(double), dimension(3) :: fc
    real(double), dimension(3,nspin) :: grho_r
    real(double), dimension(nspin) :: rho_r
    real(double), dimension(:,:,:), allocatable :: grho
    real(double), dimension(:), allocatable :: theta_r, uk_re, &
                                               uk_im, u_vdW_r
    real(double), dimension(:,:), allocatable :: phi_k
    real(double), dimension(:,:), allocatable :: u_vdW
    complex(double_cplx), dimension(:,:), allocatable :: rcp_u_vdW

    real(double) :: E_vdW_nl, E_gga_x, E_lda_c, E_vdW_dusp
    real(double) :: d3k
    ! timer
    type(cq_timer) :: tmr_get_qmesh, tmr_set_kcut, tmr_build_gradient,&
                      tmr_theta, tmr_fft_theta, tmr_int_E_vdW, &
                      tmr_decusp, tmr_gga_x, tmr_lda_c

    ! note that due to the large amount of memory consumed by arrays
    ! like u_vdW and rcp_u_vdW, we reuse these for thetas too.
    ! However this means theta has to be recalculated point-wise once
    ! u_vdW's are computed and stored in the array u_vdW.

    allocate(grho(maxngrid,3,nspin), theta_r(1:nq), uk_re(1:nq), &
             uk_im(1:nq), u_vdW_r(1:nq), phi_k(1:nq,1:nq), &
             u_vdW(maxngrid,1:nq), rcp_u_vdW(maxngrid,1:nq), STAT=stat)
    if (stat /= 0) call cq_abort("vdWXC_energy: Error alloc mem: ", maxngrid, nq)
    call reg_alloc_mem(area_ops, maxngrid*3*nspin + 4*nq + nq*nq + &
                       maxngrid*nq + 2*maxngrid*nq, type_dbl)

    ! allocate VdW arrays
    call start_timer(tmr_get_qmesh, WITH_LEVEL)
    call vdW_get_qmesh(n_q)
    call stop_print_timer(tmr_get_qmesh, "vdW: setting up q-mesh", &
                          IPRINT_TIME_THRES2)

    ! work out kcut, this should be less or equal to Nyquist critical
    ! frequency, fc = 2*pi/2*dr = pi/dr, dr is the real-space grid
    ! spacing. Also cell in conquest is orthorhombic
    fc(1) = pi / (r_super_x / n_grid_x)
    fc(2) = pi / (r_super_y / n_grid_y)
    fc(3) = pi / (r_super_z / n_grid_z)
    kcut = minval(fc)
    d3k = one / (r_super_x * r_super_y * r_super_z)

    ! set kcut
    call start_timer(tmr_set_kcut, WITH_LEVEL)
    call vdW_set_kcut(kcut)
    call stop_print_timer(tmr_set_kcut, &
                          "vdW: setting up kcut and calc phi &
                          &interpolation table", IPRINT_TIME_THRES2)

    ! get gradient of densities
    call start_timer(tmr_build_gradient, WITH_LEVEL)
    do spin = 1, nspin
       call build_gradient(rho(:,spin), grho(:,:,spin), maxngrid)
    end do
    call stop_print_timer(tmr_build_gradient, "vdW: build grho", &
                          IPRINT_TIME_THRES2)

    ! Find expansion of theta_iq(r) = rho(r) * P_iq(q(r))
    call start_timer(tmr_theta, WITH_LEVEL)
    grho_r = zero
    do ir = 1, n_my_grid_points
       rho_r(1:nspin) = rho(ir,1:nspin)
       do ix = 1, 3
          grho_r(ix,1:nspin) = grho(ir,ix,1:nspin)
       end do
       call vdW_theta(rho_r, grho_r, theta_r)
       ! use u_vdW to store values of theta
       u_vdW(ir,1:nq) = theta_r(1:nq)
    end do
    call stop_print_timer(tmr_theta, "vdW: build theta", &
                          IPRINT_TIME_THRES2)

    ! Fourier transform theta_iq(r), use rcp_u_vdW to store its values
    call start_timer(tmr_fft_theta, WITH_LEVEL)
    do iq = 1, nq
       call fft3(u_vdW(:,iq), rcp_u_vdW(:,iq), maxngrid, -1)
    end do
    ! scale to make the fft consistent with radial fft on phi
    rcp_u_vdW = grid_point_volume * rcp_u_vdW * sqrt(real(maxngrid, double))
    call stop_print_timer(tmr_fft_theta, "vdW: FFT theta", &
                          IPRINT_TIME_THRES2)

    call start_timer(tmr_int_E_vdW, WITH_LEVEL)
    E_vdW_nl = zero
    do ik = 1, z_columns_node(inode) * n_grid_z
       ! find the value of |k|
       mag_k = sqrt(recip_vector(ik,1)**2 + &
                    recip_vector(ik,2)**2 + &
                    recip_vector(ik,3)**2)
       if (mag_k < kcut) then
          call vdW_phi(mag_k, phi_k)
          do jq = 1, nq
             do iq = 1, nq
                E_vdW_nl = E_vdW_nl +                        &
                           phi_k(iq,jq) *                    &
                           (real(rcp_u_vdW(ik,iq), double) * &
                            real(rcp_u_vdW(ik,jq), double) + &
                            aimag(rcp_u_vdW(ik,iq)) *        &
                            aimag(rcp_u_vdW(ik,jq)))
             end do
          end do
       end if
    end do
    E_vdW_nl = half * d3k * E_vdW_nl
    call gsum(E_vdW_nl)
    call stop_print_timer(tmr_int_E_vdW, &
                          "vdW: integrate to get non local vdW energy", &
                          IPRINT_TIME_THRES2)

    ! calculate other energy terms
    E_gga_x = zero
    E_lda_c = zero
    E_vdW_dusp = zero
    do ir = 1, n_my_grid_points
       rho_r(1:nspin) = rho(ir,1:nspin)
       do ix = 1, 3
          grho_r(ix,1:nspin) = grho(ir,ix,1:nspin)
       end do
       rho_tot_r = rho_r(1) + rho_r(nspin)
       ! Get local cusp correction to nonlocal vdW energy integral
       call start_timer(tmr_decusp, WITH_LEVEL)
       call vdW_decusp(rho_r, grho_r, eps_c_usp_r)
       E_vdW_dusp = E_vdW_dusp + rho_tot_r * (eps_C_usp_r)
       call stop_timer(tmr_decusp)
       ! Get exchange from revPBE GGA
       call start_timer(tmr_gga_x, WITH_LEVEL)
       call eps_xc_of_r_GGA_PBE(nspin, functional_gga_pbe96_rev98, &
                                rho_r, grho_r, eps_x=eps_gga_x)
       E_gga_x = E_gga_x + rho_tot_r * eps_gga_x
       call stop_timer(tmr_gga_x)
       ! Get correlation form PW92 LDA
       call start_timer(tmr_lda_c, WITH_LEVEL)
       call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_c=eps_lda_c)
       E_lda_c = E_lda_c + rho_tot_r * eps_lda_c
       call stop_timer(tmr_lda_c)
    end do

    call start_timer(tmr_decusp, WITH_LEVEL)
    call gsum(E_vdW_dusp)
    E_vdW_dusp = grid_point_volume * E_vdW_dusp
    call stop_print_timer(tmr_decusp, "vdW: soft phi correction", &
                          IPRINT_TIME_THRES2)

    call start_timer(tmr_gga_x, WITH_LEVEL)
    call gsum(E_gga_x)
    E_gga_x = grid_point_volume * E_gga_x
    call stop_print_timer(tmr_gga_x, "vdW: exchange from gga", &
                          IPRINT_TIME_THRES2)

    call start_timer(tmr_lda_c, WITH_LEVEL)
    call gsum(E_lda_c)
    E_lda_c = grid_point_volume * E_lda_c
    call stop_print_timer(tmr_lda_c, "vdW: correlation from lda", &
                          IPRINT_TIME_THRES2)

    E_vdW_nl = E_vdW_nl + E_vdW_dusp

! LT_debug 2012/04/30 begin
!     if (inode == ionode) then
!        print *, "LT: E_vdW_nl = ", E_vdW_nl * en_conv, en_units(energy_units)
!        print *, "LT: E_gga_x = ", E_gga_x * en_conv, en_units(energy_units)
!        print *, "LT: E_lda_c = ", E_lda_c * en_conv, en_units(energy_units)
!     end if
! LT_debug 2012/04/30 end

    E_vdWXC = E_gga_x + E_vdW_nl + E_lda_c

    ! ! Get vdW u (eq 11,12 of PRL 103, 096102 (2009))
    ! do ik = 1, z_columns_node(inode) * n_grid_z
    !    ! find the value of |k|
    !    mag_k = sqrt(recip_vector(ik,1)**2 + &
    !                 recip_vector(ik,2)**2 + &
    !                 recip_vector(ik,3)**2)
    !    if (mag_k < kcut) then
    !       ! note phi_k is real
    !       call vdW_phi(mag_k, phi_k)
    !       uk_re = zero
    !       uk_im = zero
    !       do jq = 1, nq
    !          do iq = 1, nq
    !             tk_re = real(rcp_u_vdW(ik,jq), double)
    !             tk_im = aimag(rcp_u_vdW(ik,jq))
    !             uk_re(iq) = uk_re(iq) + phi_k(iq,jq) * tk_re
    !             uk_im(iq) = uk_im(iq) + phi_k(iq,jq) * tk_im
    !          end do ! iq
    !       end do ! jq
    !       ! at this point the ik data slot in rcp_u_vdW is free to be reused
    !       rcp_u_vdW(ik,1:nq) = cmplx(uk_re(1:nq), uk_im(1:nq))
    !    else
    !       rcp_u_vdW(ik,1:nq) = zero
    !    end if
    ! end do ! ik

    ! ! Fourier transform u_vdW_k back to u_vdW_r (this replaces thetas
    ! ! stored in the array)
    ! ! do iq = 1, nq
    ! !    call fft3(u_vdW(:,iq), rcp_u_vdW(:,iq), maxngrid, +1)
    ! ! end do

    ! ! Get vdW non-local correlation energy
    ! ! E_vdWXC = zero
    ! do ir = 1, n_my_grid_points
    !    rho_r(1:nspin) = rho(ir,1:nspin)
    !    do ix = 1, 3
    !       grho_r(ix,1:nspin) = grho(ir,ix,1:nspin)
    !    end do
    !    rho_tot_r = rho_r(1) + rho_r(nspin)
    !    ! make sure it is not zero
    !    rho_tot_r = max(rho_tot_r, RD_ERR)
    !    ! u_vdW_r(1:nq) = u_vdW(ir,1:nq)
    !    do ix = 1, 3
    !       grho_r(ix,1:nspin) = grho(ir,ix,1:nspin)
    !    end do
    !    ! Get local cusp correction to nonlocal vdW energy integral
    !    call vdW_decusp(rho_r, grho_r, eps_c_usp_r)
    !    ! Need to recalculate theta again
    !    call vdW_theta(rho_r, grho_r, theta_r)
    !    ! Add contributions to non-local energy density

    !    ! eps_vdW_nl_r = eps_c_usp_r + &
    !    !                half * grid_point_volume * &
    !    !                sum(u_vdW_r * theta_r) / rho_tot_r

    !    eps_vdW_nl_r = eps_c_usp_r + &
    !                   half * &
    !                   sum(u_vdW_r * theta_r) / rho_tot_r

    !    ! Get exchange from revPBE GGA
    !    call eps_xc_of_r_GGA_PBE(nspin, functional_gga_pbe96_rev98, &
    !                             rho_r, grho_r, eps_x=eps_gga_x)
    !    ! Get correlation form PW92 LDA
    !    call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_c=eps_lda_c)
    !    ! Add contributions to the exchange energy density
    !    eps_x = eps_gga_x
    !    ! Add contributions to the correlation energy density
    !    eps_c = eps_lda_c + eps_vdW_nl_r

    !    ! add correction to nonlocal vdW energy
    !    E_vdWXC = E_vdWXC + rho_tot_r * (eps_x + eps_c)
    ! end do
    ! call gsum(E_vdWXC)
    ! E_vdWXC = grid_point_volume * E_vdWXC

    deallocate(grho, theta_r, uk_re, uk_im, u_vdW_r, phi_k, u_vdW, &
               rcp_u_vdW, STAT=stat)
    if (stat /= 0) call cq_abort("vdWXC_energy: Error dealloc mem")
    call reg_dealloc_mem(area_ops, maxngrid*3*nspin + 4*nq + nq*nq + &
                         maxngrid*nq + 2*maxngrid*nq, type_dbl)

  end subroutine vdWXC_energy
  !!*****


!   subroutine vdWXC_energy(rho, E_vdWXC)

!     use datatypes
!     use numbers
!     use dimens
!     use global_module, only: nspin, spin_factor,       &
!                              functional_gga_pbe96_rev98
!     use maxima_module, only: maxngrid
!     use XC_module,     only: build_gradient,           &
!                              eps_xc_of_r_GGA_PBE,      &
!                              Vxc_of_r_LSDA_PW92,       &
!     use fft_module,    only: fft3, z_columns_node, recip_vector
!     use GenComms,      only: gsum, inode

!     implicit none

!     ! passed parameters
!     real(double), dimension(:,:), intent(in)  :: rho
!     real(double),                 intent(out) :: E_vdWXC
!     ! local variable
!     integer      :: n_q
!     integer      :: iq, jq, ik, ia, ib, ix, spin, rr
!     real(double) :: k, xc_energy, E_vdW_nl, E_x_gga, E_c_lda
!     real(double) :: tk_re, tk_im
!     real(double),         dimension(3) :: reciprocal, g
!     real(double),         dimension(maxngrid,nspin)   :: rho_tmp
!     real(double),         dimension(maxngrid)         :: rho_tot,      &
!                                                          mod_grho_tot, &
!                                                          xc_epsilon,   &
!                                                          xc_potential, &
!                                                          eps_x_gga,    &
!                                                          eps_c_lda,    &
!                                                          eps_c_usp
!     real(double),         dimension(maxngrid,3)       :: grho_tot
!     real(double),         dimension(maxngrid,nspin)   :: mod_grho
!     real(double),         dimension(maxngrid,3,nspin) :: grho
!     real(double),         dimension(1:nq)             :: uk_re, uk_im
!     real(double),         dimension(1:nq,1:nq)        :: phi_k, dphi_kdk
!     real(double),         dimension(maxngrid,1:nq)    :: theta_r, u_vdW_r
!     complex(double_cplx), dimension(maxngrid,1:nq)    :: theta_k, u_vdW_k


! ! LT_debug 2012/04/18 begin
!     real(double) :: E_usp
! ! LT_debug 2012/04/18 end

!     ! allocate VdW arrays
!     call vdW_get_qmesh(n_q)

!     ! work out kcut
!     ! cell in conquest is always orthorhombic
!     reciprocal(1) = (two * pi / r_super_x) * n_grid_x
!     reciprocal(2) = (two * pi / r_super_y) * n_grid_y
!     reciprocal(3) = (two * pi / r_super_z) * n_grid_z
!     kcut = half * minval(reciprocal)

!     print *, 'LT_debug: kcut_set = ', kcut_set, kcut

!     call vdW_set_kcut(kcut)

!     print *, 'LT_debug: kcut_set = ', kcut_set, kcut

!     ! Find vacuum value of theta
!     ! LT: 2012/04/09 may not be necessary for Conquest
!     ! rho_tmp = zero
!     ! grho = zero
!     ! call vdW_theta(rho_tmp, grho, theta_vac)

!     ! get gradient of densities
!     do spin = 1, nspin
!        call build_gradient(rho(:,spin), grho(:,:,spin), maxngrid)
!     end do

!     ! Find expansion of theta_iq(r) = rho(r) * P_iq(q(r))
!     do rr = 1, n_my_grid_points
!        rho_r(1:nspin) = rho(rr,1:nspin)
!        do ix = 1, 3
!           grho_r(ix,1:nspin) = grho(rr,ix,1:nspin)
!        end do
!        call vdW_theta(rho_r, grho_r, theta_r)
!        theta(rr,1:nq) = theta_r(1:nq)
!     end do

!     ! get rho_tot and grho_tot
!     rho_tot = spin_factor * sum(rho, 2)
!     grho_tot = spin_factor * sum(grho, 3)

!     ! Fourier transform theta_iq(r)
!     do iq = 1, nq
!        call fft3(theta(:,iq), rcp_theta(:,iq), maxngrid, -1)
!     end do

!     ! Get vdW u (eq 11,12 of PRL 103, 096102 (2009))
!     do ik = 1, z_columns_node(inode) * n_grid_z
!        ! find the value of |k|
!        k = zero
!        do ix = 1, 3
!           k = k + recip_vector(ik,ix) * recip_vector(ik,ix)
!        end do
!        k = sqrt(k)
!        if (k < kcut) then
!           ! note phi_k is real
!           call vdW_phi(k, phi_k, dphi_kdk)
!           uk_re = zero
!           uk_im = zero
!           do jq = 1, nq
!              do iq = 1, nq
!                 tk_re = real(rcp_theta(ik,jq), double)
!                 tk_im = aimag(rcp_theta(ik,jq))
!                 uk_re(iq) = uk_re(iq) + tk_re * phi_k(iq,jq)
!                 uk_im(iq) = uk_im(iq) + tk_im * phi_k(iq,jq)
!              end do
!           end do
!           u_vdW_k(ik,1:nq) = cmplx(uk_re(1:nq), uk_im(1:nq))
!        else
!           u_vdW_k(ik,1:nq) = zero
!        end if
!     end do

!     ! Fourier transform u_vdW_k back to u_vdW_r
!     do iq = 1, nq
!        call fft3(u_vdW_r(:,iq), u_vdW_k(:,iq), maxngrid, +1)
!     end do

!     ! Get local cusp correction to nonlocal VdW energy integral
!     call vdW_decusp(rho, grho, eps_c_usp)

! ! LT_debug 2012/04/18 begin
!     E_usp = zero
!     do rr = 1, n_my_grid_points
!        E_usp = E_usp + eps_c_usp(rr) * rho_tot(rr) * grid_point_volume
!     end do
!     call gsum(E_usp)
!     print *, 'LT_debug: grid_point_volume, E_usp = ', grid_point_volume, E_usp
! ! LT_debug 2012/04/18 end

!     ! Get vdW non-local correlation energy
!     E_vdW_nl = zero
!     do rr = 1, n_my_grid_points
!        do iq = 1, nq
!           ! scale by grid_point_volume as this is required for u_vdW_r
!           ! (since it is grid_point_volume * discrete convolution)
!           E_vdW_nl = E_vdW_nl + half * theta_r(rr,iq) * &
!                      grid_point_volume * u_vdW_r(rr,iq)
!        end do
!        ! add the local cusp correction
!        E_vdW_nl = E_vdW_nl + rho_tot(rr) * eps_c_usp(rr)
!     end do
!     E_vdW_nl = grid_point_volume * E_vdW_nl
!     call gsum(E_vdW_nl)

!     ! Get exchange part of GGA
!     call get_xc_potential_GGA_PBE(rho_tot, xc_potential, xc_epsilon,  &
!                                   xc_energy, maxngrid,                &
!                                   flavour=functional_gga_pbe96_rev98, &
!                                   x_epsilon=eps_x_gga)

!     print *, 'LT_debug: xc_energy_GGA = ', xc_energy

!     E_x_gga = zero
!     do rr = 1, n_my_grid_points
!        E_x_gga = E_x_gga + rho_tot(rr) * eps_x_gga(rr)
!     end do
!     E_x_gga = grid_point_volume * E_x_gga
!     call gsum(E_x_gga)

!     ! Get correlation part of LDA
!     call get_xc_potential_LDA_PW92(rho_tot, xc_potential, xc_epsilon, &
!                                    xc_energy, maxngrid, c_epsilon=eps_c_lda)

!     print *, 'LT_debug: xc_energy_LDA = ', xc_energy

!     E_c_lda = zero
!     do rr = 1, n_my_grid_points
!        E_c_lda = E_c_lda + rho_tot(rr) * eps_c_lda(rr)
!     end do
!     E_c_lda = grid_point_volume * E_c_lda
!     call gsum(E_c_lda)


!     print *, "LT_debug: E_x_gga, E_c_lda, E_vdW_nl: ", &
!              E_x_gga, E_c_lda, E_vdW_nl

!     ! gather all energies together
!     E_vdWXC = E_x_gga + E_c_lda + E_vdW_nl

!   end subroutine vdWXC_energy


  !!****f* vdWDFT_module/vdWXC_energy_slow
  !! PURPOSE
  !!   Calculates the van der Waals xc energy given by Dion et al. PRL
  !!   92, 246401 (2004) iterally
  !! USAGE
  !!   call vdWXC_energy_slow(rho, vdWXC_energy)
  !! INPUTS
  !!   real(double) rho(maxngrid,nspin) : electron densities in each spin channel
  !! OUTPUT
  !!   real(double) E_vdWXC             : vdW xc energy
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/20
  !! MODIFICATION HISTORY
  !!   Bug fix for sum over two components of rho even without spin
  !! SOURCE
  !!
  subroutine vdWXC_energy_slow(rho, E_vdWXC)

    use datatypes
    use numbers
    use global_module, only: nspin, spin_factor
    use XC,     only: build_gradient,                      &
                             eps_xc_of_r_GGA_PBE,                 &
                             Vxc_of_r_LSDA_PW92, functional_gga_pbe96_rev98
    use maxima_module, only: maxngrid
    use dimens,        only: n_my_grid_points, grid_point_volume, &
                             x_grid, y_grid, z_grid
    use grid_index,    only: grid_point_x, grid_point_y, grid_point_z
    use GenComms,      only: gsum
! LT_debug 2012/04/30 begin
!    use units
! LT_debug 2012/04/30 end

    implicit none

    ! passed parameters
    real(double), dimension(:,:), intent(in)  :: rho
    real(double),                 intent(out) :: E_vdWXC
    ! local variables
    integer :: ix, r1, r2, spin
    real(double) :: q0_r, eps_gga_x, eps_lda_c
    real(double) :: xc_energy, d1, d2, E_vdW_nl, E_x_gga, E_c_lda, &
                    r12, r1_x, r1_y, r1_z, r2_x, r2_y, r2_z, phi_d1d2
    real(double), dimension(nspin) :: rho_r
    real(double), dimension(3,nspin) :: grho_r
    real(double), dimension(:,:,:), allocatable :: grho
    real(double), dimension(:), allocatable :: rho_tot, q0,      &
                                               xc_epsilon,       &
                                               xc_potential,     &
                                               eps_x_gga, eps_c_lda

    allocate(grho(maxngrid,3,nspin), rho_tot(maxngrid), q0(maxngrid), &
             xc_epsilon(maxngrid), xc_potential(maxngrid), &
             eps_x_gga(maxngrid), eps_c_lda(maxngrid))

! LT_debug 2012/04/30 begin
!    real(double) :: E_gga_x, E_lda_c
! LT_debug 2012/04/30 end
    rho_tot = zero
    ! get gradient of densities
    do spin = 1, nspin
       call build_gradient(rho(:,spin), grho(:,:,spin), maxngrid)
       ! get the total density
       rho_tot(:) = rho_tot(:) + spin_factor * rho(:,spin)
    end do


    ! get q0
    do r1 = 1, n_my_grid_points
       rho_r(1:nspin) = rho(r1,1:nspin)
       do ix = 1, 3
          grho_r(ix,1:nspin) = grho(r1,ix,1:nspin)
       end do
       call qofrho(rho_r, grho_r, q0(r1))
    end do

    ! loop over the grid points for r1 and r2
    E_vdW_nl = zero
    do r1 = 1, n_my_grid_points
       if (rho_tot(r1) > RD_ERR) then
          do r2 = 1, n_my_grid_points
             if (rho_tot(r2) > RD_ERR) then
!                print *, "LT: slow progress: r1, r2 = ", r1, r2
                ! get r1, r2
                r1_x = (grid_point_x(r1) - 1) * x_grid
                r1_y = (grid_point_y(r1) - 1) * y_grid
                r1_z = (grid_point_z(r1) - 1) * z_grid
                r2_x = (grid_point_x(r2) - 1) * x_grid
                r2_y = (grid_point_y(r2) - 1) * y_grid
                r2_z = (grid_point_z(r2) - 1) * z_grid
                r12 = sqrt((r1_x - r2_x)**2 + &
                           (r1_y - r2_y)**2 + &
                           (r1_z - r2_z)**2)
                ! get d1, d2
                d1 = q0(r1) * r12
                d2 = q0(r2) * r12
                ! get phi(d1, d2)
                phi_d1d2 = phi_val(d1, d2)
                ! accumulate into the integral
                E_vdW_nl = E_vdW_nl + rho_tot(r1) * rho_tot(r2) * phi_d1d2
             end if
          end do ! r2
       end if
    end do ! r1
    call gsum(E_vdW_nl)
    E_vdW_nl = half * (grid_point_volume**2) * E_vdW_nl

    ! Get exchange part of GGA and correlation part of LDA
    E_vdWXC = zero
! LT_debug 2012/06/11 begin
    ! E_gga_x = zero
    ! E_lda_c = zero
! LT_debug 2012/06/11 end
    do r1 = 1, n_my_grid_points
       rho_r(1:nspin) = rho(r1,1:nspin)
       do ix = 1, 3
          grho_r(ix,1:nspin) = grho(r1,ix,1:nspin)
       end do
       call eps_xc_of_r_GGA_PBE(nspin, functional_gga_pbe96_rev98, &
                                rho_r, grho_r, eps_x=eps_gga_x)
       call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_c=eps_lda_c)
       E_vdWXC = E_vdWXC + rho_tot(r1) * (eps_gga_x + eps_lda_c)
! LT_debug 2012/06/11 begin
       ! E_gga_x = E_gga_x + rho_tot(r1) * eps_gga_x
       ! E_lda_c = E_lda_c + rho_tot(r1) * eps_lda_c
! LT_debug 2012/06/11 end
    end do
    call gsum(E_vdWXC)
    E_vdWXC = E_vdWXC * grid_point_volume
! LT_debug 2012/06/11 begin
    ! call gsum(E_gga_x)
    ! E_gga_x = E_gga_x * grid_point_volume
    ! call gsum(E_lda_c)
    ! E_lda_c = E_lda_c * grid_point_volume
! LT_debug 2012/06/11 end

    ! Add contribution from E_vdW_nl
    E_vdWXC = E_vdWXC + E_vdW_nl

! LT_debug 2012/06/11 begin
!     print *, "LT: E_vdW_nl = ", E_vdW_nl * en_conv, en_units(energy_units)
!     print *, "LT: E_gga_x = ", E_gga_x * en_conv, en_units(energy_units)
!     print *, "LT: E_lda_c = ", E_lda_c * en_conv, en_units(energy_units)
! LT_debug 2012/06/11 end

    deallocate(grho, rho_tot, q0, xc_epsilon, xc_potential, &
               eps_x_gga, eps_c_lda)

  end subroutine vdWXC_energy_slow
  !!*****

end module vdWDFT_module

