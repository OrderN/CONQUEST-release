! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_xc
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_xc *
!!
!!NAME
!! cq_ut_psp2pao_xc
!!PURPOSE
!! Exchange-correlation potential and energy related tasks
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_xc

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine get_v_xc
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_xc/get_v_xc *
!!
!!NAME
!! get_v_xc
!!USAGE
!!
!!PURPOSE
!! Build the xc potential
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine get_v_xc

     use cq_ut_psp2pao_global

     implicit none

     ! Local variables
     integer :: stat

     if(gl_xcflag == 0) then    ! LDA PW92
        ! Note that no spin: only the first gl_points_mesh points returned
        call get_xc_potential_LDA_PW92(gl_v_xc)
     else                       ! GGA
        call get_xc_potential_GGA_PBE(gl_v_xc)
     end if

  end subroutine get_v_xc
!!***


! -----------------------------------------------------------------------------
! Subroutine get_dcE_xc
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_xc/get_dcE_xc *
!!
!!NAME
!! get_dcE_xc
!!USAGE
!!
!!PURPOSE
!! Get the double counting term of the exchange-correlation energy
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 15/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine get_dcE_xc(dcE_xc)

     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_numeric

     implicit none

     !! Passed variables

     ! Result
     real(double) :: dcE_xc

     !! Local variables
     real(double), dimension(gl_points_mesh) :: e_xc
     integer :: i

     !! Parameters
     real(double), parameter :: four_pi = 4.0_double * 3.1415926535897932_double

     if(gl_xcflag == 0) then    ! LDA PW92
        ! Note that no spin: only the first gl_points_mesh points returned
        call get_xc_potential_LDA_PW92(gl_v_xc, e_xc)
     else                       ! GGA
        ! Note that no spin: only the first gl_points_mesh points returned
        call get_xc_potential_GGA_PBE(gl_v_xc, e_xc)
     end if

     do i=1, gl_points_mesh
        gl_tmp1(i) = four_pi * gl_r(i) * gl_r(i) &
                   * ( gl_rho(i) + gl_rho(i + gl_points_mesh) ) &
                   * ( e_xc(i) - gl_v_xc(i))
     end do

     ! Get result, in dcE_xc
     call simpson_integral_prod(gl_r, gl_tmp1, gl_points_mesh, gl_dnu, dcE_xc)

  end subroutine get_dcE_xc
!!***


! -----------------------------------------------------------
! Subroutine get_xc_potential_LDA_PW92
! -----------------------------------------------------------

!!****f* H_matrix_module/get_xc_potential_LDA_PW92 *
!!
!!NAME
!! get_xc_potential_LDA_PW92
!!USAGE
!!
!!PURPOSE
!! Calculates the exchange-correlation potential
!! on the grid within LDA using the Ceperley-Alder
!! interpolation formula. It also calculates the
!! total exchange-correlation energy.
!!
!! Note that this is the Perdew-Wang parameterisation of the Ceperley-Alder
!! results for a homogeneous electron gas, as described in Phys. Rev. B 45, 13244 (1992),
!! with Ceperley-Alder in Phys. Rev. Lett. 45, 566 (1980)
!!INPUTS
!!
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!AUTHOR
!! A.S. Torralba
!!CREATION DATE
!! 06/03/06
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine get_xc_potential_LDA_PW92(xc_potential, xc_energy)

    use cq_ut_psp2pao_types
    use cq_ut_psp2pao_global

    implicit none

    ! Passed variables
    real(double) :: xc_potential(gl_points_mesh)
    real(double), optional :: xc_energy(gl_points_mesh)

    !     Local variables
    integer n
    real(double) :: prefactor, postfactor, denominator, &
                    e_correlation, e_exchange, &
                    rcp_rs, rs, sq_rs, &
                    v_correlation, v_exchange, &
                    delta_prefactor, delta_postfactor
    real(double) :: rho


    !     From Table I, Phys. Rev. B 45, 13244 (1992), for reference
    real(double), parameter :: alpha  = 1.0421234_double
    real(double), parameter :: alpha1 = 0.21370_double
    real(double), parameter :: beta1  = 7.5957_double
    real(double), parameter :: beta2  = 3.5876_double
    real(double), parameter :: beta3  = 1.6382_double
    real(double), parameter :: beta4  = 0.49294_double
    real(double), parameter :: A      = 0.031091_double

    !     Constants
    real(double), parameter :: very_small = 1.0e-30_double
    real(double), parameter :: zero = 0.0_double
    real(double), parameter :: one = 1.0_double
    real(double), parameter :: third = 1.0_double/3.0_double
    real(double), parameter :: four_thirds = 4.0_double/3.0_double
    real(double), parameter :: k00 = 1.611991954_double     ! (4*pi/3)**(1/3)
    real(double), parameter :: k01 = -0.458165347_double    ! -3/(2*pi*alpha)
    real(double), parameter :: k02 = -0.062182_double       ! -2*A
    real(double), parameter :: k03 = -0.0132882934_double   ! -2*A*alpha1
    real(double), parameter :: k04 = 0.4723158174_double    ! 2*A*beta1
    real(double), parameter :: k05 = 0.2230841432_double    ! 2*A*beta2
    real(double), parameter :: k06 = 0.1018665524_double    ! 2*A*beta3
    real(double), parameter :: k07 = 0.03065199508_double   ! 2*A*beta4
    real(double), parameter :: k08 = -0.008858862267_double ! 2*k03/3
    real(double), parameter :: k09 = 0.0787193029_double    ! k04/6
    real(double), parameter :: k10 = 0.074361381067_double  ! k05/3
    real(double), parameter :: k11 = 0.0509332762_double    ! k06/2
    real(double), parameter :: k12 = 0.0204346633867_double ! 2*k07/3


    do n=1,gl_points_mesh ! loop over grid pts and store potl on each
       rho = gl_rho(n) + gl_rho(n + gl_points_mesh)
       if (rho > very_small) then ! Find radius of hole
          rcp_rs = k00 * ( rho**third )
       else
          rcp_rs = zero
       end if

       ! ENERGY

       ! Exchange
       e_exchange = k01 * rcp_rs

       ! Correlation
       if (rcp_rs > zero) then
          rs = one/rcp_rs
       else
          rs = zero
       end if
       sq_rs = dsqrt(rs)

       prefactor = k02 + k03*rs
       denominator = sq_rs * ( k04 + sq_rs * ( k05 + sq_rs * ( k06 + k07 * sq_rs)))
       if (denominator > zero) then
          postfactor = log( 1 + 1/denominator )
       else
          postfactor = 0
       end if

       e_correlation = prefactor * postfactor


       ! Both exchange and correlation

       if(present(xc_energy)) then
          xc_energy(n) = e_exchange + e_correlation
       end if


       ! POTENTIAL

       ! Exchange

       v_exchange = four_thirds * e_exchange

       ! Correlation
       !   (derivative of rho * e_correlation)

       ! NOTE: delta_prefactor is actually the derivative of rho*prefactor
       !       delta_postfactor is rho times the derivative of postfactor
       delta_prefactor  = k02 + k08*rs
       if (sq_rs > zero) then
          delta_postfactor = sq_rs * ( k09 &
                           + sq_rs*(k10 + sq_rs*( k11 + k12 * sq_rs ))) &
                           / ( denominator * ( 1 + denominator ) )
       else
          delta_postfactor = 0
       end if

       v_correlation = delta_prefactor * postfactor + prefactor * delta_postfactor

       xc_potential(n) = v_exchange + v_correlation
    end do

    return
  end subroutine get_xc_potential_LDA_PW92
!!***

! -----------------------------------------------------------
! Subroutine get_xc_potential_GGA_PBE
! -----------------------------------------------------------

!!****f* H_matrix_module/get_xc_potential_GGA_PBE *
!!
!!NAME
!! get_xc_potential_GGA_PBE
!!USAGE
!!
!!PURPOSE
!! Calculates the exchange-correlation potential
!! on the grid within PBE. It also calculates the
!! total exchange-correlation energy.
!!
!!INPUTS
!!
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! A.S. Torralba
!!CREATION DATE
!! 28/01/07
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
subroutine get_xc_potential_GGA_PBE(xc_potential,xc_energy)

  use cq_ut_psp2pao_types
  use cq_ut_psp2pao_global
  use cq_ut_psp2pao_numeric

  implicit none

  ! Passed variables
  real(double) :: xc_potential(gl_points_mesh)
  real(double), optional :: xc_energy(gl_points_mesh)

  !     Local variables
  integer n

!  real(double)  :: grad_density(gl_points_mesh)
  real(double)  :: rho, grad_rho, rho1_3, rho1_6, &
                   ks, s, s2, &
                   t, t2, A, At2, &
                   factor0, factor1, factor2, factor3, factor4, &
                   denominator0, numerator1, denominator1, num_den1, numerator2, &
                   dt2_drho, dA_drho, &
                   dnumerator1_drho, ddenominator1_drho, &
                   dnumerator1_dt2, dnumerator1_dA, &
                   ddenominator1_dt2, ddenominator1_dA, &
                   dfactor2_dt2, dfactor2_dnumerator1, dfactor2_ddenominator1, &
                   dfactor2_drho, dfactor3_dfactor2, dfactor3_drho
  real(double) :: xc_energy_lda_total, xc_energy_lda(gl_points_mesh), &
                  xc_potential_lda(gl_points_mesh), &
                  e_correlation_lda, &
                  de_correlation_lda, &
                  e_exchange, e_correlation, &
                  de_exchange, de_correlation, &
                  dde_exchange, dde_correlation

  real(double) :: rho2(gl_points_mesh)
  real(double) :: df_dgrad_rho(gl_points_mesh), df_dgrad_rho2(gl_points_mesh)
  real(double) :: tmp_der

!*ast*  complex(double_cplx), dimension(3,size) :: rgradient      ! Gradient in reciprocal space

!real(double), parameter :: pi = 3.141592653589763_double 
!real(double), parameter :: alpha = 0.5210617_double

  !     From Phys. Rev. Lett. 77, 3865 (1996)
  real(double), parameter :: mu = 0.21951_double
  real(double), parameter :: kappa = 0.804_double
  real(double), parameter :: beta = 0.066725_double
  real(double), parameter :: gamma = 0.031091_double

  !     Precalculated constants
  real(double), parameter :: mu_kappa = mu/kappa                   !0.27302_double         ! mu/kappa
  real(double), parameter :: beta_gamma = beta/gamma               !2.146119_double        ! beta/gamma
  real(double), parameter :: beta_X_gamma = beta*gamma             !0.002074546_double     ! beta*gamma
  real(double), parameter :: k01 = 0.16162045967_double      ! 1/(2*(3*pi*pi)**(1/3))
  real(double), parameter :: k02 = -0.16212105381_double     ! -3*mu*((4*pi/3)**3)/(2*pi*alpha)=mu*k00*k01(LDA_PW2)=mu*k04
  real(double), parameter :: k03 = 1.98468639_double         ! ((4/pi)*(3*pi*pi)**(1/3))**(1/2)
  real(double), parameter :: k04 = -0.738558852965_double    ! -3*((4*pi/3)**3)/(2*pi*alpha) = k00*k01 in LDA_PW92
  real(double), parameter :: k05 = 0.05240415_double         ! -2*k01*k02


  ! Build the gradient of the density

!*ast* USE SPLINES INSTEAD. ONLY RADIAL PART IS INTERESTING HERE

!*ast*  call build_gradient (density, grad_density, grad_density_xyz, size)

  ! Get the LDA part of the functional
!*ast* THIS CALL NEEDS REWORKING. THE ROUTINE NEEDS OPT ARG FOR xc_energy_lda
!  call get_xc_potential_LDA_PW92(xc_potential_lda, xc_energy_lda_total, xc_energy_lda)
  call get_xc_potential_LDA_PW92(xc_potential_lda, xc_energy_lda)

!*ast* TODO Spline rho here
  ! This is for non-spin polarized calculations
  ! The other spin density is assumed to be identical
  call spline(gl_r, gl_rho, gl_points_mesh, rho2)

!  xc_energy = zero
!  delta_E_xc = zero
  do n=1,gl_points_mesh ! loop over grid pts and store potl on each
     rho = gl_rho(n) + gl_rho(n + gl_points_mesh)
 
     call spline_derivative(gl_r, gl_rho, rho2, gl_points_mesh, gl_r(n), grad_rho)
     
     ! Next line assumes spin unpolarised
     grad_rho = two * grad_rho


     !!!!!! XC GGA ENERGY

     ! Exchange

     if (rho > very_small) then
        rho1_3 = rho ** third
        rho1_6 = sqrt (rho1_3)
        s = k01 * grad_rho / (rho ** four_thirds)
        s2 = s * s
        denominator0 = 1.0 / (1.0 + mu_kappa * s2)
        factor0 = k02 * rho1_3
        factor1 = s2 * denominator0
        ! NOTE: This doesn't look like in Phys. Rev. Lett. 77:18, 3865 (1996)
        !       because the 1 in Fx, has been multiplied by Ex-LDA and is implicit
        !       in xc_energy_lda(n), in the total energy below
        e_exchange = factor0 * factor1
     else
        e_exchange = zero
     end if


     ! Correlation

     if (rho > very_small) then
        e_correlation_lda = xc_energy_lda(n) - k04 * rho1_3

        ! t=grad_rho/(2*rho*ks); ks=sqrt(4*kf/pi); kf=(3*pi*pi*rho)**(1/3); s=grad_rho/(2*rho*kf)
        ks = k03 * rho1_6
        t = grad_rho / ( 2 * ks * rho )
        t2 = t * t

        A = exp(-e_correlation_lda / gamma) - 1.0
        if (A > very_small) then
           A = beta_gamma / A
        else
           A = beta_gamma * BIG
        end if

        At2 = A * t2
        numerator1 = 1.0 + At2
        denominator1 = 1.0 + At2 + At2 * At2
        num_den1 = numerator1 / denominator1

        factor2 = t2 * num_den1
        factor3 = gamma * log(one + beta_gamma * factor2)

        e_correlation = factor3;  !gamma * log( 1.0 + beta_gamma * t2 * num_den1 )
     else
        e_correlation = zero
     end if

     ! Both exchange and correlation
!     xc_energy = xc_energy + (xc_energy_lda(n) + e_exchange + e_correlation)*rho

     if(present(xc_energy)) then
        xc_energy(n) = xc_energy_lda(n) + e_exchange + e_correlation
     end if

     !!!!!! POTENTIAL

     !!! Terms due to df/drho

     ! Exchange

     if (rho > very_small) then
        de_exchange = four_thirds * factor0 * factor1 * ( 1 - 2* denominator0 )
     else
        de_exchange = zero
     end if

     ! Correlation

     if (rho > very_small) then
        de_correlation_lda = ( xc_potential_lda(n) &
                             - four_thirds * k04 * rho1_3 &
                             - e_correlation_lda ) / rho

        dt2_drho = -seven_thirds * t2 / rho
        dA_drho  = A * A * exp(-e_correlation_lda / gamma ) * de_correlation_lda / beta
        dnumerator1_dt2   = A
        dnumerator1_dA    = t2
        factor4 = 1.0 + 2 * At2
        ddenominator1_dt2 = A * factor4
        ddenominator1_dA  = t2 * factor4
        dnumerator1_drho   = dnumerator1_dt2 * dt2_drho &
                           + dnumerator1_dA * dA_drho
        ddenominator1_drho = ddenominator1_dt2 * dt2_drho &
                           + ddenominator1_dA * dA_drho
        dfactor2_dt2 = num_den1
        dfactor2_dnumerator1   = t2 / denominator1
        dfactor2_ddenominator1 = -factor2 / denominator1
        dfactor2_drho = dfactor2_dt2 * dt2_drho &
                      + dfactor2_dnumerator1 * dnumerator1_drho &
                      + dfactor2_ddenominator1 * ddenominator1_drho
        dfactor3_dfactor2 = beta / ( 1.0 + beta_gamma * factor2 )
        dfactor3_drho = dfactor3_dfactor2 * dfactor2_drho
        de_correlation = factor3 + rho * dfactor3_drho
     else
        de_correlation = zero
     end if

     xc_potential(n) = xc_potential_lda(n) + de_exchange + de_correlation


     !!! Terms due to df/d|grad_rho|

     ! Exchange

     if (rho > very_small) then
        dde_exchange = -k05 * s * denominator0 * denominator0!factor1 * denominator0
     else
        dde_exchange = zero
     end if

     ! Correlation

     if (rho > very_small) then
        numerator2 = beta_X_gamma * t * ( 1.0 + 2.0 * At2 ) &
                   / (( gamma * denominator1 + beta * t2 * numerator1 ) * denominator1)
        dde_correlation = numerator2 / ks;
     else
        dde_correlation = zero
     end if

     df_dgrad_rho(n) = dde_exchange + dde_correlation
  end do ! do gl_points_mesh

  ! Gradient times derivative of energy
  call spline(gl_r, df_dgrad_rho, gl_points_mesh, df_dgrad_rho2)

  ! Add divergence to potential
  do n=1,gl_points_mesh
     call spline_derivative(gl_r, df_dgrad_rho, df_dgrad_rho2, gl_points_mesh, gl_r(n), tmp_der)
     xc_potential(n) = xc_potential(n) &
                     - ( two*df_dgrad_rho(n)/gl_r(n) + tmp_der )
  end do

  return
end subroutine get_xc_potential_GGA_PBE
!!***

end module cq_ut_psp2pao_xc
