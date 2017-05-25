! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module XC_module
! ------------------------------------------------------------------------------
! Code area 3: Operators
! ------------------------------------------------------------------------------

!!****h* Conquest/XC_module
!! PURPOSE
!!   Collects together all the routines associated to
!!   exchange-correlation functionals (note the van der Waal
!!   implementations are in its separate module)
!! AUTHOR
!!   L.Tong
!! CREATION DATE
!!   2012/04/16
!! MODIFICATION HISTORY
!!   2015/05/12 08:34 dave
!!    Adding GGA stress terms
!!   2015/09/03 17:27 dave
!!    Added GGA stress term to PBE0
!! SOURCE
!!
module XC_module

  use datatypes
  use global_module, only: area_ops

  implicit none

  real(double), dimension(3), public :: XC_GGA_stress

  ! methods
  public ::                         &
       get_xc_potential,            & ! LDA-PZ81, spin non-polarised
       get_dxc_potential,           &
       get_GTH_xc_potential,        & ! LDA-GTH96, spin non-polarised
       get_GTH_dxc_potential,       &
       Vxc_of_r_LSDA_PW92,          &
       get_xc_potential_LDA_PW92,   & ! LDA-PW92, spin non-polarised
       get_dxc_potential_LDA_PW92,  &
       get_xc_potential_LSDA_PW92,  & ! LDA-PW92, spin polarised
       get_dxc_potential_LSDA_PW92, &
       get_xc_potential_GGA_PBE,    & ! GGA-PBE96 (rev98 or rev99), no spin
       eps_xc_of_r_GGA_PBE,         &
       get_dxc_potential_GGA_PBE,   &
       get_xc_potential_hyb_PBE0,   & ! hyb-PBE0
       build_gradient,              & ! calculates gradient of density
       map_density                    ! map cartisian coordinates to
                                      ! density index, obsolete
  private

    ! RCS tag for object file identification
  character(len=80), save :: &
       RCSid = "$Id$"

!!*****

contains

  !!****f* XC_module/get_xc_potential *
  !!
  !!  NAME
  !!   get_xc_potential
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the exchange-correlation potential
  !!   on the grid within LDA using the Ceperley-Alder
  !!   interpolation formula. It also calculates the
  !!   total exchange-correlation energy.
  !!
  !!   Note that this is the Perdew-Zunger parameterisation of the
  !!   Ceperley-Alder results for a homogeneous electron gas, as
  !!   described in Phys. Rev. B 23, 5048 (1981), with Ceperley-Alder in
  !!   Phys. Rev. Lett. 45, 566 (1980)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   02/03/95
  !!  MODIFICATION HISTORY
  !!   01/11/2005 Antonio
  !!    Bibliography:
  !!       Exchange energy: It can be found in:
  !!                        Phys. Rev. B 45, 13244 (1992)
  !!                      **See Eq. 26
  !!       Correlation energy: The correlation functional is PZ-81
  !!                           Phys. Rev. B 23, 5048 (1981)
  !!                         **See Appendix C, and specially Table XII,
  !!                             Eq. C1 for the xc hole (rs),
  !!                             Eqs. C3 & C4, for rs > 1 and
  !!                             Eqs. C5 & C6, for rs < 1
  !!   17/05/2001 dave
  !!    Converted to F90, added ROBODoc header, shortened call
  !!   08/06/2001 dave
  !!    Changed to use gsum from GenComms and added RCS Id and Log tags
  !!   21/06/2001 dave
  !!    Included in H_matrix_module
  !!   08:00, 2003/03/12 dave
  !!    Added calculation of correction term to band energy
  !!   11:49, 30/09/2003 drb
  !!    Added comments to make separate testing of X and C easier
  !!   2008/03/03 18:32 dave
  !!    Removed dsqrt
  !!   2011/12/12 L.Tong
  !!    Removed third, it is now defined in numbers module
  !!   2012/04/03 L.Tong
  !!   - Added optional parameters x_epsilon and c_epsilon for output of
  !!     the exchange and correlation parts of xc_epsilon respectively.
  !!   2014/09/24 L.Truflandier
  !!   - Added optional x_energy for output
  !!  SOURCE
  !!
  subroutine get_xc_potential(density, xc_potential, xc_epsilon,     &
                              xc_energy, size, x_epsilon, c_epsilon, &
                              x_energy)

    use datatypes
    use numbers
    use GenComms,               only: gsum
    use dimens,                 only: grid_point_volume, n_my_grid_points

    implicit none

    ! Passed variables
    integer,                    intent(in)  :: size
    real(double),               intent(out) :: xc_energy
    real(double), dimension(:), intent(in)  :: density
    real(double), dimension(:), intent(out) :: xc_potential, xc_epsilon
    ! optional
    real(double), dimension(:), intent(out), optional :: x_epsilon, c_epsilon
    real(double),               intent(out), optional :: x_energy

    ! Local variables
    integer      :: n
    real(double) :: denominator, e_correlation, e_exchange, ln_rs, &
                    numerator, rcp_rs, rho, rs, rs_ln_rs, sq_rs,   &
                    v_correlation, v_exchange
    real(double), parameter :: alpha  = -0.45817_double
    real(double), parameter :: beta_1 =  1.0529_double
    real(double), parameter :: beta_2 =  0.3334_double
    real(double), parameter :: gamma  = -0.1423_double
    real(double), parameter :: p =  0.0311_double
    real(double), parameter :: q = -0.048_double
    real(double), parameter :: r =  0.0020_double
    real(double), parameter :: s = -0.0116_double

    xc_energy = zero
    if (present(x_epsilon)) x_epsilon = zero
    if (present(c_epsilon)) c_epsilon = zero
    if (present(x_energy))  x_energy  = zero

    do n = 1, n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       if (rho > RD_ERR) then ! Find radius of hole
          rcp_rs = ( four_thirds * pi * rho )**(third)
       else
          rcp_rs = zero
       end if
       e_exchange = alpha * rcp_rs
       if (present(x_epsilon)) x_epsilon(n) = e_exchange
       v_exchange = four_thirds * e_exchange
       if (rcp_rs>zero) then
          rs = one/rcp_rs
       else
          rs = zero
       end if
       sq_rs = sqrt(rs)
       if (rs>=one) then
          denominator = one / (one + beta_1 * sq_rs + beta_2 * rs)
          numerator   = one + seven_sixths * beta_1 * sq_rs +  &
               four_thirds * beta_2 * rs
          e_correlation = gamma * denominator
          v_correlation = gamma * numerator * denominator * denominator
       else if ((rs<one).and.(rs>RD_ERR)) then
          ln_rs    = log(rs)
          rs_ln_rs = rs * ln_rs
          e_correlation = p * ln_rs + q  + r * rs_ln_rs + s * rs
          v_correlation = e_correlation -  &
               third * (p + s * rs + r * (rs_ln_rs + rs))
       else
          e_correlation = zero
          v_correlation = zero
       end if
       if (present(c_epsilon)) c_epsilon(n) = e_correlation
       if (present(x_energy))  x_energy     = x_energy + e_exchange * density(n)
       ! Both X and C
       xc_energy       = xc_energy  + (e_exchange + e_correlation) * density(n)
       xc_potential(n) = v_exchange + v_correlation
       xc_epsilon(n)   = e_exchange + e_correlation
       ! These two for testing
       ! Just C
       !xc_energy = xc_energy+e_correlation*density(n)
       !xc_potential(n) = v_correlation
       !xc_epsilon(n) = e_correlation
       ! Just X
       !xc_energy = xc_energy+e_exchange*density(n)
       !xc_potential(n) = v_exchange
       !xc_epsilon(n) = e_exchange
    end do ! do n_my_grid_points
    call gsum(xc_energy)
    if (present(x_energy)) call gsum(x_energy)
    ! and 'integrate' the energy over the volume of the grid point
    xc_energy = xc_energy * grid_point_volume
    if (present(x_energy)) x_energy =  x_energy * grid_point_volume

    return
  end subroutine get_xc_potential
  !!***


  !!****f* XC_module/get_GTH_xc_potential *
  !!
  !!  NAME
  !!   get_xc_potential
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the exchange-correlation potential
  !!   on the grid within LDA using the Ceperley-Alder
  !!   interpolation formula. It also calculates the
  !!   total exchange-correlation energy.
  !!
  !!   Note that this is the Goedecker/Teter/Hutter formula which
  !!   involves only ratios of polynomials, and is rather easy to
  !!   differentiate.  See PRB 54, 1703 (1996)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:45, 25/03/2003
  !!  MODIFICATION HISTORY
  !!   2011/12/12 L.Tong
  !!     Removed third, it is now defined in numbers module
  !!  SOURCE
  !!
  subroutine get_GTH_xc_potential(density, xc_potential, xc_epsilon, &
                                  xc_energy, size)

    use datatypes
    use numbers
    use global_module,          only: io_lun
    use GenComms,               only: cq_abort, gsum
    use dimens,                 only: grid_point_volume, n_my_grid_points

    implicit none

    ! Passed variables
    integer,      intent(in)  :: size
    real(double), intent(out) :: xc_energy
    real(double), dimension(:), intent(in)  :: density
    real(double), dimension(:), intent(out) :: xc_potential, xc_epsilon

    ! Local variables
    integer n
    real(double) :: denominator, e_correlation, e_exchange, ln_rs,    &
                    numerator, rcp_rs, rho, rs, rs_ln_rs, sq_rs,      &
                    v_correlation, v_exchange, drs_dRho, t1, t2, dt1, &
                    dt2
    real(double), parameter :: a0=0.4581652932831429_double
    real(double), parameter :: a1=2.217058676663745_double
    real(double), parameter :: a2=0.7405551735357053_double
    real(double), parameter :: a3=0.01968227878617998_double
    real(double), parameter :: b1=1.000000000000000_double
    real(double), parameter :: b2=4.504130959426697_double
    real(double), parameter :: b3=1.110667363742916_double
    real(double), parameter :: b4=0.02359291751427506_double

    xc_energy = zero
    do n = 1, n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       if (rho > RD_ERR) then ! Find radius of hole
          rcp_rs = ( four*third * pi * rho )**(third)
          rs     = one/rcp_rs
          !if (rs < 0.01_double) write (io_lun, *) 'rs out of range ', n
       else
          rcp_rs = zero
          rs     = zero
          !write (io_lun, *) 'rho out of range ', n
       end if
       if (rs > zero) then
          drs_dRho = -rs / (3.0 * rho)
          t1  = a0 + rs*(a1 + rs * (a2 + rs * a3))
          t2  = rs * (b1 + rs * (b2 + rs * (b3 + rs * b4)))
          dt1 = a1 + rs * (2.0 * a2 + rs * 3.0 * a3)
          dt2 = b1 + rs * (2.0 * b2 + rs * (3.0 * b3 + rs * 4.0 * b4))
          xc_energy = xc_energy - (t1/t2)*rho
          xc_potential(n) = &
               -(t1/t2) + rho*drs_dRho * (-dt1 / t2 + t1 * dt2 / (t2 * t2))
          xc_epsilon(n) = -t1/t2   ! 2010.Oct.30 TM
       else
          xc_potential(n) = zero
       end if
    end do ! do n_my_grid_points
    call gsum(xc_energy)
    ! and 'integrate' the energy over the volume of the grid point
    xc_energy = xc_energy * grid_point_volume
    return
  end subroutine get_GTH_xc_potential
  !!***


  !!****f* XC_module/get_xc_potential_LDA_PW92 (Obsolete)*
  !!
  !!  NAME
  !!   get_xc_potential_LDA_PW92 (Obsolete)
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the exchange-correlation potential on the grid within
  !!   LDA using the Ceperley-Alder interpolation formula. It also
  !!   calculates the total exchange-correlation energy.
  !!
  !!   Note that this is the Perdew-Wang parameterisation of the
  !!   Ceperley-Alder results for a homogeneous electron gas, as
  !!   described in Phys. Rev. B 45, 13244 (1992), with Ceperley-Alder
  !!   in Phys. Rev. Lett. 45, 566 (1980)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   01/11/05
  !!  MODIFICATION HISTORY
  !!   2008/03/03 18:32 dave
  !!    Removed dsqrt
  !!   2011/03/22 L.Tong
  !!    Removed local definition of third, it is now defined in numbers
  !!    module Changed 1 to one in the log() functions
  !!   2012/04/03 L.Tong
  !!   - Added optional variables x_epsilon and c_epsilon so that if
  !!     present they will store the output for the exchange and
  !!     correlation parts of xc_epsilon respectively
  !!   2012/04/25 L.Tong
  !!   *** THIS SUBROUTINE IS NOW OBSOLETE, USE THE LSDA VERSION ***
  !!  SOURCE
  !!
  subroutine get_xc_potential_LDA_PW92(density, xc_potential,       &
                                       xc_epsilon, xc_energy_total, &
                                       size, x_epsilon, c_epsilon)

    use datatypes
    use numbers
    use GenComms, only: gsum
    use dimens,   only: grid_point_volume, n_my_grid_points

    implicit none

    ! Passed variables
    integer,                    intent(in)  :: size
    real(double),               intent(out) :: xc_energy_total
    real(double), dimension(:), intent(in)  :: density
    real(double), dimension(:), intent(out) :: xc_epsilon, xc_potential
    ! optional
    real(double), dimension(:), intent(out), optional :: x_epsilon, c_epsilon

    ! Local variables
    integer      :: n
    real(double) :: prefactor, postfactor, denominator, &
                    e_correlation, e_exchange,          &
                    rcp_rs, rho, rs, sq_rs,             &
                    v_correlation, v_exchange,          &
                    delta_prefactor, delta_postfactor


    ! From Table I, Phys. Rev. B 45, 13244 (1992), for reference
    real(double), parameter :: alpha  = 1.0421234_double
    real(double), parameter :: alpha1 = 0.21370_double
    real(double), parameter :: beta1  = 7.5957_double
    real(double), parameter :: beta2  = 3.5876_double
    real(double), parameter :: beta3  = 1.6382_double
    real(double), parameter :: beta4  = 0.49294_double
    real(double), parameter :: A      = 0.031091_double

    ! Precalculated constants
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

    xc_energy_total = zero
    if (present(x_epsilon)) x_epsilon = zero
    if (present(c_epsilon)) c_epsilon = zero

    do n = 1, n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)

       if (rho > RD_ERR) then ! Find radius of hole
          rcp_rs = k00 * ( rho**third )
       else
          rcp_rs = zero
       end if

       ! ENERGY

       ! Exchange
       e_exchange = k01 * rcp_rs
       if (present(x_epsilon)) x_epsilon(n) = e_exchange

       ! Correlation
       if (rcp_rs > zero) then
          rs = one/rcp_rs
       else
          rs = zero
       end if
       sq_rs = sqrt(rs)

       prefactor = k02 + k03*rs
       denominator = sq_rs * ( k04 + sq_rs * ( k05 + sq_rs * ( k06 + k07 * sq_rs)))
       if (denominator > zero) then
          postfactor = log( one + one/denominator )
       else
          postfactor = 0
       end if

       e_correlation = prefactor * postfactor
       if (present(c_epsilon)) c_epsilon(n) = e_correlation

       ! Both exchange and correlation

       xc_epsilon(n)   = e_exchange + e_correlation
       xc_energy_total = xc_energy_total + xc_epsilon(n)*rho

       ! POTENTIAL

       ! Exchange

       v_exchange = four_thirds * e_exchange

       ! Correlation
       !   (derivative of rho * e_correlation)

       ! NOTE: delta_prefactor is actually the derivative of rho*prefactor
       !       delta_postfactor is rho times the derivative of postfactor
       delta_prefactor  = k02 + k08*rs
       if (sq_rs > zero) then
          delta_postfactor = &
               sq_rs * (k09 + sq_rs*(k10 + sq_rs*(k11 + k12 * sq_rs))) / &
                       (denominator * (1 + denominator))
       else
          delta_postfactor = 0
       end if

       v_correlation = delta_prefactor * postfactor + prefactor * delta_postfactor

       xc_potential(n) = v_exchange + v_correlation


    end do ! do n_my_grid_points
    call gsum(xc_energy_total)
    ! and 'integrate' the energy over the volume of the grid point
    xc_energy_total = xc_energy_total * grid_point_volume

    return
  end subroutine get_xc_potential_LDA_PW92
  !!***


  !!****f* XC_modue/Vxc_of_r_LSDA_PW92 *
  !! PURPOSE
  !!   Calculates the exchange-correlation energy density and
  !!   potential as a function of rho(r) (evaluated at single grid
  !!   point r).
  !!
  !!   The LDA Functional is given in Perdew and Wang, PRB 45, 13244
  !!   (1992)
  !! USAGE
  !!   To calculate Vx and Vc
  !!     call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_x, eps_c, Vx, Vc)
  !!   To calculate Vx or Vc only
  !!     call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_x=eps_x, Vx=Vx)
  !!     call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_c=eps_c, Vx=Vc)
  !!   To calculate eps_x and/or eps_c without Vx and Vc
  !!     call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_x=eps_x)
  !!     call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_c=eps_c)
  !!     call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_x=eps_x, eps_c=eps_c)
  !! INPUTS
  !!   integer      nspin          : 1 = spin nonpolarised, 2 = spin polarised
  !!   real(double) rho_r(1:nspin) : value of spin dependent density at r
  !! OUTPUT
  !!   (All are optional)
  !!   real(double) eps_x       : value of exchange energy density at r
  !!   real(double) eps_c       : value of correlation energy density at r
  !!   real(double) Vx(1:nspin) : value of exchange potential at r
  !!   real(double) Vc(1:nspin) : value of correlation potential at r
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/25
  !! MODIFICATION HISTORY
  !! 2012/05/27 L.Tong
  !! - I found that with rounding error, even if rho_r(2) is just tiny
  !!   bit negative, and rho_r(1) positive, this will make rho_tot
  !!   still positive and above RD_ERR, and hence pass the zero rho
  !!   test, but zeta = (rho(1) - rho(2)) / rho_tot will now be
  !!   greater than one. This makes one - zeta < 0 and hence (one -
  !!   zeta)**third undefined.  Therefore to fix this I added a
  !!   constraint that zeta must be in between -one and one.
  !! - Used separate subroutine PW92_G for G(rs)
  !! SOURCE
  !!
  subroutine Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_x, eps_c, Vx, Vc)

    use datatypes
    use numbers
    use GenComms, only: cq_abort

    implicit none

    ! passed variables
    integer,                              intent(in)  :: nspin
    real(double), dimension(:),           intent(in)  :: rho_r
    real(double),               optional, intent(out) :: eps_x, eps_c
    real(double), dimension(:), optional, intent(out) :: Vx, Vc
    ! local variables
    real(double) :: rho_tot_r, rs, sq_rs, rcp_rs, rcp_sq_rs, zeta,     &
                    eps_c0, deps_c0_drs, eps_c1, deps_c1_drs, alpha_c, &
                    dalpha_c_drs, f_function, df_function_dzeta,       &
                    alpha_cDfpp0, dalpha_c_drsDfpp0, deps_c_drs,       &
                    deps_c_dzeta
    real(double) :: onePzeta, oneMzeta, onePzeta_1_3, oneMzeta_1_3,    &
                    onePzeta_4_3, oneMzeta_4_3, zeta_3, zeta_4,        &
                    factor_4_3, factor_1_3
    ! parameters used for convinience of calculation
    real(double), parameter :: K01 =  1.042123522_double ! 2*(9*pi/4)**(-1/3)
    real(double), parameter :: K02 =  1.611991954_double ! (4*pi/3)**(1/3)
    real(double), parameter :: K03 = -0.458165293_double ! -3/(2*pi*K01)

    ! check passed parameters
    if (present(Vx) .and. (.not. present(eps_x))) &
         call cq_abort('Vxc_of_r_LSDA_PW92: Vx and eps_x must present &
                        &at the same time')
    if (present(Vc) .and. (.not. present(eps_c))) &
         call cq_abort('Vxc_of_r_LSDA_PW92: Vc and eps_c must present &
                        &at the same time')

    ! get total density
    rho_tot_r = rho_r(1) + rho_r(nspin);

    ! catch the case when density is zero
    if (rho_tot_r <= RD_ERR) then
       if (present(eps_x)) eps_x = zero
       if (present(eps_c)) eps_c = zero
       if (present(Vx))    Vx = zero
       if (present(Vc))    Vc = zero
       return
    end if

    ! From this point on-wards rho_tot_r is greater than zero

    ! work out rs, 1/rs
    rcp_rs = K02 * (rho_tot_r**third) ! this calculates 1/rs
    rs = one / rcp_rs
    sq_rs = sqrt(rs)
    rcp_sq_rs = sqrt(rcp_rs)

    if (nspin == 2) then ! case for spin polarised

       ! zeta related factors
       zeta = (rho_r(1) - rho_r(nspin)) / rho_tot_r
       ! zeta must be in between -one and one, we need to make sure
       ! of this, rounding errors (when rho_r(spin) are small) can
       ! sometimes make zeta go outside this range.
       zeta = max(-one, zeta)
       zeta = min(one, zeta)
       onePzeta = one + zeta
       oneMzeta = one - zeta
       onePzeta_1_3 = onePzeta**third
       oneMzeta_1_3 = oneMzeta**third
       onePzeta_4_3 = onePzeta * onePzeta_1_3
       oneMzeta_4_3 = oneMzeta * oneMzeta_1_3
       zeta_3 = zeta**3
       zeta_4 = zeta * zeta_3
       factor_4_3 = onePzeta_4_3 + oneMzeta_4_3
       factor_1_3 = onePzeta_1_3 - oneMzeta_1_3

       ! Exchange Part
       ! eps_x: eq (26) PRB 45, 13244 (1992)
       if (present(eps_x)) eps_x = half * K03 * rcp_rs * factor_4_3
       ! Vx (exchange potential)
       if (present(Vx)) then
          Vx(1) = four_thirds * &
                  (eps_x + oneMzeta * K03 * rcp_rs * half * factor_1_3)
          Vx(2) = four_thirds * &
                  (eps_x - onePzeta * K03 * rcp_rs * half * factor_1_3)
       end if

       ! Correlation Part

       if (present(eps_c)) then

          if (present(Vc)) then
             ! Work out eps_c0, and deps_c0/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.031091_double,     &
                         alpha1=0.21370_double, &
                         beta1=7.5957_double,   &
                         beta2=3.5876_double,   &
                         beta3=1.6382_double,   &
                         beta4=0.49294_double,  &
                         p=one,                 &
                         G=eps_c0,              &
                         dG_drs=deps_c0_drs)
             ! Work out eps_c1 and deps_c1/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.015545_double,     &
                         alpha1=0.20548_double, &
                         beta1=14.1189_double,  &
                         beta2=6.1977_double,   &
                         beta3=3.3662_double,   &
                         beta4=0.62517_double,  &
                         p=one,                 &
                         G=eps_c1,              &
                         dG_drs=deps_c1_drs)
             ! Work out alpha_c, and dalpha_c/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.016887_double,     &
                         alpha1=0.11125_double, &
                         beta1=10.357_double,   &
                         beta2=3.6231_double,   &
                         beta3=0.88026_double,  &
                         beta4=0.49671_double,  &
                         p=one,                 &
                         G=alpha_c,             &
                         dG_drs=dalpha_c_drs)
             ! note that the table parameters are for -alpha_c
             alpha_c = -alpha_c
             dalpha_c_drs = -dalpha_c_drs
          else ! calculate eps_c only
             ! Work out eps_c0, and deps_c0/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.031091_double,     &
                         alpha1=0.21370_double, &
                         beta1=7.5957_double,   &
                         beta2=3.5876_double,   &
                         beta3=1.6382_double,   &
                         beta4=0.49294_double,  &
                         p=one,                 &
                         G=eps_c0)
             ! Work out eps_c1 and deps_c1/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.015545_double,     &
                         alpha1=0.20548_double, &
                         beta1=14.1189_double,  &
                         beta2=6.1977_double,   &
                         beta3=3.3662_double,   &
                         beta4=0.62517_double,  &
                         p=one,                 &
                         G=eps_c1)
             ! Work out alpha_c, and dalpha_c/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.016887_double,     &
                         alpha1=0.11125_double, &
                         beta1=10.357_double,   &
                         beta2=3.6231_double,   &
                         beta3=0.88026_double,  &
                         beta4=0.49671_double,  &
                         p=one,                 &
                         G=alpha_c)
             ! note that the table parameters are for -alpha_c
             alpha_c = -alpha_c
          end if

          ! Work out f_function eq.(8) PRB 45, 13244 (1992),
          ! and df_function/dzeta
          f_function = 1.923661051_double * (factor_4_3 - two)
          if (present(Vc)) then
             df_function_dzeta = 2.564881401_double * (factor_1_3)
          end if

          ! Get eps_c, eq.(8) PRB 45, 13244 (1992)
          ! alpha_c / f''(0)
          alpha_cDfpp0 = 0.584822340_double * alpha_c
          eps_c = eps_c0 +                                      &
                  f_function * (alpha_cDfpp0 * (one - zeta_4) + &
                                (eps_c1 - eps_c0) * zeta_4)

          if (present(Vc)) then
             ! Get deps_c_drs and deps_c_dzeta
             ! dalpha_c_drs / f''(0)
             dalpha_c_drsDfpp0 = 0.584822340_double * dalpha_c_drs
             ! eq.(A2) of PRB 45, 13244 (1992)
             deps_c_drs =                                                     &
                  deps_c0_drs + f_function *                                  &
                  (zeta_4 * (deps_c1_drs - deps_c0_drs - dalpha_c_drsDfpp0) + &
                   dalpha_c_drsDfpp0)
             ! eq.(A3) of PRB 45, 13244 (1992)
             deps_c_dzeta = (four * zeta_3 * f_function +      &
                             df_function_dzeta * zeta_4) *     &
                            (eps_c1 - eps_c0 - alpha_cDfpp0) + &
                            df_function_dzeta * alpha_cDfpp0

             ! Get Vc
             Vc(1) = eps_c - third * rs * deps_c_drs + oneMzeta * deps_c_dzeta
             Vc(2) = eps_c - third * rs * deps_c_drs - onePzeta * deps_c_dzeta
          end if

       end if ! present(eps_c)

    else ! case for spin non-polarised

       ! Exchange Part

       ! eps_x: eq (26) PRB 45, 13244 (1992)
       if (present(eps_x)) eps_x = K03 * rcp_rs
       ! Vx (exchange potential)
       if (present(Vx)) Vx(1) = four_thirds * eps_x

       ! Correlation Part

       if (present(eps_c)) then

          if (present(Vc)) then
             ! Work out eps_c0, and deps_c0/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.031091_double,     &
                         alpha1=0.21370_double, &
                         beta1=7.5957_double,   &
                         beta2=3.5876_double,   &
                         beta3=1.6382_double,   &
                         beta4=0.49294_double,  &
                         p=one,                 &
                         G=eps_c0,              &
                         dG_drs=deps_c0_drs)
          else
             ! Work out eps_c0, and deps_c0/drs
             ! fitting parameters taken from Table I,
             ! Phys. Rev. B 45, 13244 (1992)
             call PW92_G(rs,                    &
                         A=0.031091_double,     &
                         alpha1=0.21370_double, &
                         beta1=7.5957_double,   &
                         beta2=3.5876_double,   &
                         beta3=1.6382_double,   &
                         beta4=0.49294_double,  &
                         p=one,                 &
                         G=eps_c0)
          end if

          ! Get eps_c, eq.(8) PRB 45, 13244 (1992)
          eps_c = eps_c0

          if (present(Vc)) then
             ! Get deps_c_drs, eq.(A2) of PRB 45, 13244 (1992)
             deps_c_drs = deps_c0_drs
             ! Get Vc
             Vc(1) = eps_c - third * rs * deps_c_drs
          end if

       end if ! present(eps_c)

    end if ! if (nspin == 2)

    return
  end subroutine Vxc_of_r_LSDA_PW92
  !!*****


  !!****f* XC_module/get_xc_potential_LSDA_PW92 *
  !!
  !!  NAME
  !!   get_xc_potential_LSDA_PW92
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the spin polarized exchange-correlation
  !!   potential on the grid within LDA using the Ceperley-Alder
  !!   interpolation formula. It also calculates the total
  !!   exchange-correlation energy.
  !!
  !!   Note that this is the Perdew-Wang parameterisation of the
  !!   Ceperley-Alder results for a homogeneous electron gas, as
  !!   described in Phys. Rev. B 45, 13244 (1992), with Ceperley-Alder
  !!   in Phys. Rev. Lett. 45, 566 (1980) INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   L. Tong
  !!  CREATION DATE
  !!   22/03/2011
  !!  MODIFICATION HISTORY
  !!   2012/03/14 L.Tong
  !!   - Changed spin implementation
  !!   2012/04/03 L.Tong
  !!   - Added optional variables to output exchange and correlation
  !!     part of epsilon separately, if required
  !!   2012/04/25 L.Tong
  !!   - Removed the optional variables for separate exchange and
  !!     correlation energy densities, these are now covered in
  !!     subroutine Vxc_of_r_LSDA_PW92
  !!   - Moved all the local point-wise calculations of epsilon and
  !!     potential to subroutine Vxc_of_r_LSDA_PW92, this allows the
  !!     point-wise calculations to be done in other parts of the
  !!     program if required---this saves a lot of memory in vdWDFT
  !!     calculations.
  !!   2014/09/24 L.Truflandier
  !!   - Added optional x_energy_total for output
  !!  SOURCE
  !!
  subroutine get_xc_potential_LSDA_PW92(density, xc_potential,       &
                                        xc_epsilon, xc_energy_total, &
                                        size, x_energy_total         )

    use datatypes
    use numbers
    use GenComms,               only: cq_abort, gsum
    use dimens,                 only: grid_point_volume, n_my_grid_points
    use global_module,          only: nspin

    implicit none

    ! Passed variables
    ! size of the real space grid
    integer,                      intent(in)  :: size
    real(double), dimension(:,:), intent(in)  :: density
    real(double), dimension(:,:), intent(out) :: xc_potential
    real(double), dimension(:),   intent(out) :: xc_epsilon
    real(double),                 intent(out) :: xc_energy_total
    real(double), optional,       intent(out) :: x_energy_total

    ! Local variables
    integer      :: rr, spin
    real(double) :: eps_x, eps_c, rho_tot_r
    real(double), dimension(nspin) :: rho_r, Vx, Vc

    ! initialisation
    if (present(x_energy_total)) x_energy_total  = zero
    xc_energy_total = zero

    ! loop over grid points on each node
    do rr = 1, n_my_grid_points
       rho_r(1:nspin) = density(rr,1:nspin)
       rho_tot_r      = rho_r(1) + rho_r(nspin)
       call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_x=eps_x, eps_c=eps_c,&
                               Vx=Vx, Vc=Vc)

       xc_epsilon(rr)           = eps_x       + eps_c
       xc_potential(rr,1:nspin) = Vx(1:nspin) + Vc(1:nspin)
       xc_energy_total = xc_energy_total + xc_epsilon(rr) * rho_tot_r
       if (present(x_energy_total)) x_energy_total  =  x_energy_total + eps_x * rho_tot_r
    end do
    call gsum(xc_energy_total)
    if (present(x_energy_total)) call gsum(x_energy_total)

    xc_energy_total = xc_energy_total * grid_point_volume
    if (present(x_energy_total)) x_energy_total =  x_energy_total * grid_point_volume

    return
  end subroutine get_xc_potential_LSDA_PW92
  !!***


  !!****f* XC_module/eps_xc_of_r_GGA_PBE
  !! PURPOSE
  !!   Calculates the exchange and correlation energy density and its
  !!   derivatives for the PBE functional and its variants at any
  !!   local values of rho(r) and grad rho(r).
  !!
  !!   The exact flavour of the PBE functional is choosen by flavour
  !!   parameter:
  !!
  !!   flavour = functional_gga_pbe96
  !!     use original PBE, PRL 77, 3865 (1996)
  !!   flavour = functional_gga_pbe96_rev98:
  !!     use revPBE, PRL 80, 890 (1998)
  !!   flavour = functional_gga_pbe96_r99:
  !!     use RPBE, PRB 59, 7413 (1999)
  !!   flavour = functional_gga_pbe96_wc:
  !!     use Wu-Cohen, PRB 73, 235116  (2006)
  !! USAGE
  !! INPUTS
  !!   integer      nspin             : nspin = 1 for spin non-polarised,
  !!                                    nspin = 2 for spin polarised
  !!   integer      flavour           : flavour of PBE GGA to use
  !!   real(double) rho_r(1:nspin)    : value of spin density at some point r
  !!   real(double) grho_r(3,1:nspin) : vector gradient of spin density at r
  !! OUTPUT
  !!   All outputs are optional
  !!   real(double) eps_x             : value of exchange energy density at r
  !!   real(double) eps_c             : value of correlation energy density at r
  !!   real(double) drhoEps_x(0:3,1:nspin) : drhoEps_x(0,spin) =
  !!                                         d(rho_tot_r * eps_x) / drho(spin)
  !!                                         drhoEps_x(1:3,spin) =
  !!                                         d(rho_tot_r * eps_x) / dgrho(1:3,spin)
  !!   real(double) drhoEps_c(0:3,1:nspin) : drhoEps_c(0,spin) =
  !!                                         d(rho_tot_r * eps_c) / drho(spin)
  !!                                         drhoEps_x(1:3,spin) =
  !!                                         d(rho_tot_r * eps_c) / dgrho(1:3,spin)
  !!
  !!   (where rho_tot_r = rho_r(1) + rho_r(nspin))
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/25
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine eps_xc_of_r_GGA_PBE(nspin, flavour, rho_r, grho_r, eps_x,&
                                 eps_c, drhoEps_x, drhoEps_c)
    use datatypes
    use numbers
    use global_module, only: functional_gga_pbe96,       &
                             functional_gga_pbe96_rev98, &
                             functional_gga_pbe96_r99,   &
                             functional_gga_pbe96_wc
    use GenComms,      only: cq_abort

    implicit none

    ! passed parameters
    integer,                                intent(in)  :: nspin
    integer,                                intent(in)  :: flavour
    real(double), dimension(:),             intent(in)  :: rho_r
    real(double), dimension(:,:),           intent(in)  :: grho_r
    real(double),                 optional, intent(out) :: eps_x, eps_c
    real(double), dimension(0:3,nspin), optional, intent(out) :: &
         drhoEps_x, drhoEps_c
    ! local variables
    integer      :: spin, Fx_selector, ii
    real(double) :: rho_tot_r, rho_tot_s, mod_grho_tot_r,             &
                    mod_grho_tot_s, eps_x_unif, Fx, grho_tot_s_ii,    &
                    dFx_dgrho
    real(double) :: rs, kF, eps_c_unif, ks, zeta, phi, t, A, H
    real(double) :: drs_drho, dkF_drho, dphi_drho, dphi_dzeta,        &
                    deps_c_unif_drho, dH_drho, dH_dgrho, dks_drho,    &
                    dA_drho
    real(double) :: kappa, mu_kappa, s, s2, kFs, dkFs_drho, F1, F2,   &
                    F3, F4, dFx_drho, ds_drho, ds_dgrho, dt_drho,     &
                    dt_dgrho, dF1_drho, dF2_drho, dF3_drho, dF4_drho, &
                    dF3_dgrho, dF4_dgrho, Xwc, dXwc_ds,  dF1_dgrho
    real(double) :: spin_factor
    real(double), dimension(1)       :: rho_s
    real(double), dimension(1:3,1)   :: grho_s
    real(double), dimension(1:3)     :: grho_tot_r
    real(double), dimension(1:nspin) :: mod_grho_r, Vx_unif, Vc_unif, &
                                        dzeta_drho
    ! constants
    ! From PRL 77, 3865 (1996)
    real(double), parameter :: mu = 0.2195149727645171_double
    real(double), parameter :: beta = 0.066725_double
    real(double), parameter :: gamma = 0.031091_double
    real(double), parameter :: kappa_ori = 0.804_double
    ! From PRL 80, 890 (1998)
    real(double), parameter :: kappa_alt = 1.245_double
    ! From PRB 73, 235116  (2006) - Wu-Cohen
    real(double), parameter :: teneightyone = 0.123456790123_double  ! 10/81
    real(double), parameter :: c_wc = 0.00793746933516_double
    ! Precalculated constants
    real(double), parameter :: mu_kappa_ori = 0.27302_double   ! mu/kappa_ori
    real(double), parameter :: mu_kappa_alt = 0.17631_double   ! mu/kappa_alt
    real(double), parameter :: two_mu = 0.4390299455290342_double ! 2 * mu
    real(double), parameter :: beta_gamma = 2.146119_double    ! beta/gamma
    ! minimum value of rho_tot_r and mod_grho_tot_r allowed
    real(double), parameter :: min_rho = RD_ERR
    real(double), parameter :: min_mod_grho = RD_ERR
    ! Fx_selector options
    integer, parameter :: Fx_ori = 1
    integer, parameter :: Fx_alt = 2
    integer, parameter :: Fx_wc  = 3  ! Wu-Cohen Exchange

    ! check optional outputs
    if (present(drhoEps_x) .and. (.not. present(eps_x))) &
         call cq_abort('eps_xc_of_r_GGA_PBE: both drhoEps_x and eps_x &
                        &must be present')
    if (present(drhoEps_c) .and. (.not. present(eps_c))) &
         call cq_abort('eps_xc_of_r_GGA_PBE: both drhoRps_c and eps_c &
                        &must be present')

    ! choose between PBE and revPBE parameters
    select case (flavour)
    case (functional_gga_pbe96)
       Fx_selector = Fx_ori
       kappa = kappa_ori
       mu_kappa = mu_kappa_ori
    case (functional_gga_pbe96_rev98)
       Fx_selector = Fx_ori
       kappa = kappa_alt
       mu_kappa = mu_kappa_alt
    case (functional_gga_pbe96_r99)
       Fx_selector = Fx_alt
       kappa = kappa_ori
       mu_kappa = mu_kappa_ori
   case (functional_gga_pbe96_wc)
       Fx_selector = Fx_wc
       kappa = kappa_ori
       mu_kappa = mu_kappa_ori
    case default
       call cq_abort('eps_xc_of_r_GGA_PBE: flavour undefined', flavour)
    end select

    ! work out total density and its gradient
    rho_tot_r = rho_r(1) + rho_r(nspin)
    rho_tot_r = max(min_rho, rho_tot_r)
    grho_tot_r(1:3) = grho_r(1:3,1) + grho_r(1:3,nspin)
    do spin = 1, nspin
       mod_grho_r(spin) = &
            sqrt(grho_r(1,spin)**2 + grho_r(2,spin)**2 + grho_r(3,spin)**2)
    end do
    mod_grho_tot_r = &
         sqrt(grho_tot_r(1)**2 + grho_tot_r(2)**2 + grho_tot_r(3)**2)
    mod_grho_tot_r = max(min_mod_grho, mod_grho_tot_r)

    ! work out spin_factor
    if (nspin == 1) then
       spin_factor = two
    else
       spin_factor = one
    end if

    ! Exchange part
    if (present(eps_x)) then

       eps_x = zero

       do spin = 1, nspin
          ! rho_s = 2 * rho_r(spin), grho_s = 2 * |grho(spin)|
          rho_s(1) = rho_r(spin)
          rho_tot_s = max(min_rho, two * rho_s(1))
          mod_grho_tot_s = max(min_mod_grho, two * mod_grho_r(spin))
          kFs = (three * pi**2 * rho_tot_s)**third
          s = mod_grho_tot_s / (two * kFs * rho_tot_s)
          s2 = s*s
          ! choose the Fx flavour
          select case (Fx_selector)
          case (Fx_ori)
             F1 = one + mu_kappa * s2
             Fx = one + kappa - kappa / F1
          case (Fx_alt)
             F1 = exp(-mu_kappa * s2)
             Fx = one + kappa * (one - F1)
          case (Fx_wc)   !  Wu-Cohen exchange
             Xwc = teneightyone * s2 + (mu - teneightyone) *  &
                   s2 * exp(-s2) + log(one + c_wc * s2*s2)

             F1 = one +  Xwc / kappa
             Fx = one + kappa - kappa / F1

             dXwc_ds = two * teneightyone * s +                        &
                 (mu - teneightyone) * exp(-s2) * two*s * (one - s2) + &
                 four * c_wc * s*s2 / (one + c_wc * s2*s2)
          end select

          call Vxc_of_r_LSDA_PW92(1, rho_s, eps_x=eps_x_unif, Vx=Vx_unif)

          ! accumulate eps_x in spin components
          eps_x = eps_x + rho_tot_s * eps_x_unif * Fx

          ! calculate the derivative of rho_tot_r * eps_x
          if (present(drhoEps_x)) then
             dkFs_drho = third * kFs / rho_tot_s
             ds_drho = s * (-(dkFs_drho / kFs) - one / rho_tot_s)
             select case (Fx_selector)
             case (Fx_ori)
                dFx_drho = two_mu * s * ds_drho / (F1 * F1)
             case (Fx_alt)
                dFx_drho = two_mu * s * ds_drho * F1
            case (Fx_wc)
                dF1_drho = dXwc_ds * ds_drho / kappa
                dFx_drho = kappa * dF1_drho / (F1*F1)

             end select
             ! get drhoEps_x / drho(spin)
             drhoEps_x(0,spin) = Vx_unif(1) * Fx + &
                                 rho_tot_s * eps_x_unif * dFx_drho
             ! get drhoEps_x / dgrho(spin)
             do ii = 1, 3
                grho_tot_s_ii = two * grho_r(ii,spin)
                ds_dgrho = (s / mod_grho_tot_s) * grho_tot_s_ii / mod_grho_tot_s
                select case (Fx_selector)
                case (Fx_ori)
                   dFx_dgrho = two_mu * s * ds_dgrho / (F1 * F1)
                case (Fx_alt)
                   dFx_dgrho = two_mu * s * ds_dgrho * F1
                case (Fx_wc)
                   dF1_dgrho = dXwc_ds * ds_dgrho / kappa
                   dFx_dgrho = kappa * dF1_dgrho / (F1 * F1)
                end select
                drhoEps_x(ii,spin) = rho_tot_s * eps_x_unif * dFx_dgrho
             end do
          end if ! present(drhoEps_x)
       end do ! spin

       ! get eps_x
       eps_x = half * spin_factor * eps_x / rho_tot_r

    end if ! present(eps_x)

    ! Correlation part
    if (present(eps_c)) then

       ! Find local correlation energy and potential
       call Vxc_of_r_LSDA_PW92(nspin, rho_r, eps_c=eps_c_unif, Vc=Vc_unif)

       ! Find energy denisty for correlation eps_c
       rs = (three / (four * pi * rho_tot_r))**third
       kF = (three * pi**2 * rho_tot_r)**third
       ks = sqrt(four * kF / pi)

       if (nspin == 2) then

          ! zeta related terms
          zeta = (rho_r(1) - rho_r(2)) / rho_tot_r
          ! bound zeta within -1 and 1
          zeta = max(-one + min_rho, zeta)
          zeta = min(one - min_rho, zeta)
          phi = half * ((one + zeta)**two_thirds + (one - zeta)**two_thirds)
          t = mod_grho_tot_r / (two * phi * ks * rho_tot_r)
          F1 = eps_c_unif / gamma / phi**3
          F2 = exp(-F1)
          ! add RD_ERR in denominator to avoid div by 0
          A = beta_gamma / (F2 - one + RD_ERR)
          F3 = t**2 + A * t**4
          F4 = beta_gamma * F3 / (one + A * F3)
          H = gamma * phi**3 * log(one + F4)
          ! get eps_c
          eps_c = eps_c_unif + H

          ! get drhoEps_c
          if (present(drhoEps_c)) then
             drs_drho = - (third * rs / rho_tot_r)
             dkF_drho = third * kF / rho_tot_r
             dks_drho = half * ks * dkF_drho / kF
             dzeta_drho(1) = (one - zeta) / rho_tot_r
             dzeta_drho(2) = - (one + zeta) / rho_tot_r
             dphi_dzeta = half * two_thirds * &
                  (one / (one + zeta)**third - one / (one - zeta)**third)
             do spin = 1, nspin
                deps_c_unif_drho = (Vc_unif(spin) - eps_c_unif) / rho_tot_r
                dphi_drho = dphi_dzeta * dzeta_drho(spin)
                dt_drho = (- t) * (dphi_drho / phi + dks_drho / ks + &
                          one / rho_tot_r)
                ! add RD_ERR to denominator to avoid div by 0
                dF1_drho = F1 * (deps_c_unif_drho / (eps_c_unif + RD_ERR) - &
                           three * dphi_drho / phi)
                dF2_Drho = (- F2) * dF1_drho
                ! add RD_ERR to denominator to avoid div by 0
                dA_drho = A * dF2_drho / (one - F2 + RD_ERR)
                dF3_drho = (two * t + four * A * t**3) * dt_drho + dA_drho * t**4
                dF4_drho = F4 * (dF3_drho / F3 - &
                           (dA_drho * F3 + A * DF3_drho) / (one + A * F3))
                dH_drho = three * H * dphi_drho / phi + &
                     gamma * phi**3 * DF4_drho / (one + F4)
                ! get drhoEps_c / drho(spin)
                drhoEps_c(0,spin) = Vc_unif(spin) + H + rho_tot_r * dH_drho
                do ii = 1, 3
                   dt_dgrho = (t / mod_grho_tot_r) * &
                              grho_tot_r(ii) / mod_grho_tot_r
                   dF3_dgrho = dt_dgrho * (two * t + four * A * t**3)
                   dF4_dgrho = F4 * dF3_dgrho * (one / F3 - A / (one + A * F3))
                   dH_dgrho = gamma * phi**3 * dF4_dgrho / (one + F4)
                   drhoEps_c(ii,spin) = rho_tot_r * dH_dgrho
                end do
             end do ! spin
          end if ! present(rhoEps_c)

       else ! spin non-polarised case

          phi = one
          t = mod_grho_tot_r / (two * ks * rho_tot_r)
          F1 = eps_c_unif / gamma
          F2 = exp(-F1)
          ! add RD_ERR in denominator to avoid div by 0
          A = beta_gamma / (F2 - one + RD_ERR)
          F3 = t**2 + A * t**4
          F4 = beta_gamma * F3 / (one + A * F3)
          H = gamma * log(one + F4)
          ! get eps_c
          eps_c = eps_c_unif + H

          ! get drhoEps_c / drho
          if (present(drhoEps_c)) then
             drs_drho = - (third * rs / rho_tot_r)
             dkF_drho = third * kF / rho_tot_r
             dks_drho = half * ks * dkF_drho / kF
             deps_c_unif_drho = (Vc_unif(1) - eps_c_unif) / rho_tot_r
             dt_drho = (- t) * (dks_drho / ks + one / rho_tot_r)
             ! add RD_ERR to denominator to avoid div by 0
             dF1_drho = F1 * (deps_c_unif_drho / (eps_c_unif + RD_ERR))
             dF2_Drho = (- F2) * dF1_drho
             ! add RD_ERR to denominator to avoid div by 0
             dA_drho = A * dF2_drho / (one - F2 + RD_ERR)
             dF3_drho = (two * t + four * A * t**3) * dt_drho + dA_drho * t**4
             dF4_drho = F4 * (dF3_drho / F3 - &
                        (dA_drho * F3 + A * DF3_drho) / (one + A * F3))

             dH_drho = gamma * DF4_drho / (one + F4)

             ! get drhoEps_c / drho(spin)
             drhoEps_c(0,1) = Vc_unif(1) + H + rho_tot_r * dH_drho
             do ii = 1, 3
                dt_dgrho = (t / mod_grho_tot_r) * grho_tot_r(ii) / mod_grho_tot_r
                dF3_dgrho = dt_dgrho * (two * t + four * A * t**3)
                dF4_dgrho = F4 * dF3_dgrho * (one / F3 - A / (one + A * F3))
                dH_dgrho = gamma * dF4_dgrho / (one + F4)
                drhoEps_c(ii,1) = rho_tot_r * dH_dgrho
             end do
          end if ! present(drhoEps_c)

       end if ! nspin

    end if ! presnet(eps_c)

    return
  end subroutine eps_xc_of_r_GGA_PBE
  !!*****


  !!****f* XC_module/get_xc_potential_GGA_PBE
  !! PURPOSE
  !!   Calculates the exchange-correlation potential
  !!   on the grid within GGA using three flavours of
  !!   Perdew-Burke-Ernzerhof. It also calculates the
  !!   total exchange-correlation energy.
  !!
  !!   flavour not defined
  !!     use original PBE, PRL 77, 3865 (1996)
  !!   flavour = functional_gga_pbe96_rev98:
  !!     use revPBE, PRL 80, 890 (1998)
  !!   flavour = functional_gga_pbe96_r99:
  !!     use RPBE, PRB 59, 7413 (1999)
  !!
  !! USAGE
  !!   call get_xc_potential_GGA_PBE(density, xc_potential,      &
  !!                                 xc_epsilon, xc_energy, size,&
  !!                                 flavour)
  !! INPUTS
  !!   integer      size                : size of the real space grid
  !!   integer      flavour             : flavour of PBE functional (optional)
  !!   real(double) density(size,nspin) : spin dependent density
  !! OUTPUT
  !!   real(double) xc_energy                : total xc-energy (sum over spin)
  !!   real(double) xc_epsilon(size)         : xc-energy density
  !!   real(double) xc_potential(size,nspin) : xc-potnetial
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/04/27
  !! MODIFICATION HISTORY
  !!   2014/09/24 L.Truflandier
  !!   - Added optional x_energy_total for output
  !!   2015/05/12 08:30 dave
  !!   - Added stress term
  !! SOURCE
  !!
  subroutine get_xc_potential_GGA_PBE(density, xc_potential,            &
                                      xc_epsilon, xc_energy, grid_size, &
                                      flavour, x_energy )
    use datatypes
    use numbers
    use global_module, only: & !rcellx, rcelly, rcellz,          &
                             functional_gga_pbe96,            &
                             functional_gga_pbe96_rev98,      &
                             functional_gga_pbe96_r99, nspin, &
                             spin_factor
    use dimens,        only: grid_point_volume, n_my_grid_points
    use GenComms,      only: gsum, cq_abort
    use fft_module,    only: fft3, recip_vector
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! passed variables
    integer,                      intent(in)  :: grid_size
    real(double), dimension(:,:), intent(in)  :: density
    real(double), dimension(:),   intent(out) :: xc_epsilon
    real(double), dimension(:,:), intent(out) :: xc_potential
    real(double),                 intent(out) :: xc_energy
    real(double), optional,       intent(out) ::  x_energy
    integer,      optional,       intent(in)  :: flavour

    ! local variables
    integer      :: PBE_type
    integer      :: rr, spin, stat, dir
    real(double) :: eps_x, eps_c, rho_tot_r
    real(double),         dimension(nspin)        :: rho_r
    real(double),         dimension(3,nspin)      :: grho_r
    real(double),         dimension(0:3,nspin)    :: drhoEps_x, drhoEps_c
    real(double),         dimension(:,:,:), allocatable :: grad_density
    real(double),         dimension(:),     allocatable :: second_term
    complex(double_cplx), dimension(:,:),   allocatable :: rcp_drhoEps_xc

    allocate(grad_density(grid_size,3,nspin), second_term(grid_size), &
             rcp_drhoEps_xc(grid_size,3), STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_xc_potential_GGA_PBE: Error alloc mem: ", grid_size)
    call reg_alloc_mem(area_ops, grid_size*(1+3*(1+nspin)), type_dbl)

    if (present(flavour)) then
       PBE_type = flavour
    else
       PBE_type = functional_gga_pbe96
    end if

    ! initialisation
    if (present(x_energy)) x_energy = zero
    grad_density = zero
    xc_epsilon   = zero
    xc_potential = zero
    xc_energy    = zero
    XC_GGA_stress = zero

    ! Build the gradient of the density
    do spin = 1, nspin
       call build_gradient(density(:,spin), grad_density(:,:,spin), grid_size)
    end do

    do rr = 1, n_my_grid_points
       rho_tot_r = density(rr,1) + density(rr,nspin)
       do spin = 1, nspin
          rho_r(spin)      = density(rr,spin)
          grho_r(1:3,spin) = grad_density(rr,1:3,spin)
       end do
       call eps_xc_of_r_GGA_PBE(nspin, PBE_type, rho_r, grho_r, &
                                eps_x=eps_x, eps_c=eps_c,      &
                                drhoEps_x=drhoEps_x,           &
                                drhoEps_c=drhoEps_c)
       xc_epsilon(rr) = eps_x + eps_c
       xc_energy      = xc_energy + rho_tot_r * xc_epsilon(rr)
       if (present(x_energy)) x_energy = x_energy + rho_tot_r * eps_x
       ! for potential
       do spin = 1, nspin
          xc_potential(rr,spin) = drhoEps_x(0,spin) + drhoEps_c(0,spin)
          ! note that grad_density(rr,1:3,1:spin) has already been used
          ! at this point, so we can savely reuse this slot to store
          ! d(rho * eps_xc) / dgrho at rr.
          do dir=1,3
             grad_density(rr,dir,spin) = drhoEps_x(dir,spin) + drhoEps_c(dir,spin)
             XC_GGA_stress(dir) = XC_GGA_stress(dir) - grho_r(dir,spin)*grad_density(rr,dir,spin)
          end do
       end do
    end do ! rr
    call gsum(xc_energy)
    if (present(x_energy)) call gsum(x_energy)
    xc_energy = xc_energy * grid_point_volume
    call gsum(XC_GGA_stress,3)
    XC_GGA_stress = XC_GGA_stress*grid_point_volume
    !write(*,*) 'GGA stress term: ',XC_GGA_stress
    if (present(x_energy)) x_energy = x_energy * grid_point_volume

    ! add the second term to potential
    do spin = 1, nspin
       ! initialise for each spin
       rcp_drhoEps_xc = zero
       second_term    = zero
       ! FFT drhoEps_xc / dgrho(spin) to reciprocal space
       call fft3(grad_density(:,1,spin), rcp_drhoEps_xc(:,1), grid_size, -1)
       call fft3(grad_density(:,2,spin), rcp_drhoEps_xc(:,2), grid_size, -1)
       call fft3(grad_density(:,3,spin), rcp_drhoEps_xc(:,3), grid_size, -1)

       ! dot product with i * recip_vectors, and store in the first
       ! component of rcp_drhoEps_xc (used as temp storage)
       rcp_drhoEps_xc(:,1) = &
            rcp_drhoEps_xc(:,1) * minus_i * recip_vector(:,1) + &
            rcp_drhoEps_xc(:,2) * minus_i * recip_vector(:,2) + &
            rcp_drhoEps_xc(:,3) * minus_i * recip_vector(:,3)

       ! FFT back to obtain the convolution
       call fft3(second_term(:), rcp_drhoEps_xc(:,1), grid_size, +1)
       ! accumulate to potential
       do rr = 1, n_my_grid_points
          xc_potential(rr,spin) = xc_potential(rr,spin) + second_term(rr)
       end do
    end do ! spin

    deallocate(grad_density, second_term, rcp_drhoEps_xc, STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_xc_potential_GGA_PBE: Error dealloc mem")
    call reg_dealloc_mem(area_ops, grid_size*(1+3*(1+nspin)), type_dbl)

  end subroutine get_xc_potential_GGA_PBE
  !!*****


  !!****f* XC_module/get_xc_potential_hyb_PBE0
  !! PURPOSE
  !!
  !!   Work routine based on get_xc_potential_GGA_PBE.
  !!   It will be recoded in near futur to handle any hybrid functional.
  !!   For now just handle one coefficient as in PBE0.
  !! USAGE
  !!
  !! INPUTS
  !!
  !!   integer      size                : size of the real space grid
  !!   integer      flavour             : flavour of PBE functional (optional)
  !!   real(double) density(size,nspin) : spin dependent density
  !! OUTPUT
  !!   real(double) xc_energy                : total xc-energy (sum over spin)
  !!   real(double) x_energy                 : exchange energy only (sum over spin)
  !!   real(double) xc_epsilon(size)         : xc-energy density
  !!   real(double) xc_potential(size,nspin) : xc-potential
  !!
  !! AUTHOR
  !!   L.Tong/L.Truflandier
  !! CREATION DATE
  !!   2014/09/24
  !! MODIFICATION HISTORY
  !!   2015/09/03 17:27 dave
  !!   Added GGA XC stress term
  !! SOURCE
  !!
  subroutine get_xc_potential_hyb_PBE0(density, xc_potential, exx_a,     &
                                       xc_epsilon, xc_energy, grid_size, &
                                       flavour, x_energy )
    use datatypes
    use numbers
    !use global_module, only: rcellx, rcelly, rcellz, spin_factor, nspin, &
    use global_module, only: spin_factor, nspin, &
                             functional_gga_pbe96
    use dimens,        only: grid_point_volume, n_my_grid_points
    use GenComms,      only: gsum, cq_abort
    use fft_module,    only: fft3, recip_vector
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! passed variables
    integer,                      intent(in)  :: grid_size
    real(double), dimension(:,:), intent(in)  :: density
    real(double), dimension(:),   intent(out) :: xc_epsilon
    real(double), dimension(:,:), intent(out) :: xc_potential
    real(double),                 intent(out) :: xc_energy
    real(double), optional,       intent(out) :: x_energy
    integer,      optional,       intent(in)  :: flavour
    real(double),                 intent(in)  :: exx_a

    ! local variables
    integer      :: PBE_type
    integer      :: rr, spin, stat
    real(double) :: eps_x, eps_c, rho_tot_r
    real(double),         dimension(nspin)        :: rho_r
    real(double),         dimension(3,nspin)      :: grho_r
    real(double),         dimension(0:3,nspin)    :: drhoEps_x, drhoEps_c
    real(double),         dimension(:,:,:), allocatable :: grad_density
    real(double),         dimension(:),     allocatable :: second_term
    complex(double_cplx), dimension(:,:),   allocatable :: rcp_drhoEps_xc

    allocate(grad_density(grid_size,3,nspin), second_term(grid_size), &
             rcp_drhoEps_xc(grid_size,3), STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_xc_potential_GGA_PBE: Error alloc mem: ", grid_size)
    call reg_alloc_mem(area_ops, grid_size*(1+3*(1+nspin)), type_dbl)

    if (present(flavour)) then
       PBE_type = flavour
    else
       PBE_type = functional_gga_pbe96
    end if

    ! setup exx_a

    ! Initialisation
    if (present(x_energy)) x_energy = zero
    grad_density = zero
    xc_epsilon   = zero
    xc_potential = zero    
    xc_energy    = zero

    ! Build the gradient of the density
    do spin = 1, nspin
       call build_gradient(density(:,spin), grad_density(:,:,spin), grid_size)
    end do

    do rr = 1, n_my_grid_points
       rho_tot_r = density(rr,1) + density(rr,nspin)
       do spin = 1, nspin
          rho_r(spin)      = density(rr,spin)
          grho_r(1:3,spin) = grad_density(rr,1:3,spin)
       end do
       call eps_xc_of_r_GGA_PBE(nspin, PBE_type, rho_r, grho_r, &
                                eps_x    =eps_x, &
                                eps_c    =eps_c, &
                                drhoEps_x=drhoEps_x, &
                                drhoEps_c=drhoEps_c)
       xc_epsilon(rr) = exx_a * eps_x + eps_c
       xc_energy      = xc_energy + rho_tot_r * xc_epsilon(rr)
       if (present(x_energy)) x_energy =  x_energy + exx_a * rho_tot_r * eps_x
       ! for potential
       do spin = 1, nspin
          xc_potential(rr,spin) = exx_a * drhoEps_x(0,spin) + drhoEps_c(0,spin)
          ! note that grad_density(rr,1:3,1:spin) has already been used
          ! at this point, so we can savely reuse this slot to store
          ! d(rho * eps_xc) / dgrho at rr.
          grad_density(rr,1:3,spin) = exx_a * drhoEps_x(1:3,spin) + &
                                      drhoEps_c(1:3,spin)
          XC_GGA_stress(1:3) = XC_GGA_stress(1:3) - grho_r(1:3,spin)*grad_density(rr,1:3,spin)
       end do
    end do ! rr
    call gsum(xc_energy)
    if (present(x_energy)) call gsum(x_energy)
    xc_energy = xc_energy * grid_point_volume
    call gsum(XC_GGA_stress,3)
    XC_GGA_stress = XC_GGA_stress*grid_point_volume
    if (present(x_energy)) x_energy =  x_energy * grid_point_volume

    ! add the second term to potential
    do spin = 1, nspin
       ! initialise for each spin
       rcp_drhoEps_xc = zero
       second_term    = zero
       ! FFT drhoEps_xc / dgrho(spin) to reciprocal space
       call fft3(grad_density(:,1,spin), rcp_drhoEps_xc(:,1), grid_size, -1)
       call fft3(grad_density(:,2,spin), rcp_drhoEps_xc(:,2), grid_size, -1)
       call fft3(grad_density(:,3,spin), rcp_drhoEps_xc(:,3), grid_size, -1)

       ! dot product with i * recip_vectors, and store in the first
       ! component of rcp_drhoEps_xc (used as temp storage)
       rcp_drhoEps_xc(:,1) = &
            rcp_drhoEps_xc(:,1) * minus_i * recip_vector(:,1) + &
            rcp_drhoEps_xc(:,2) * minus_i * recip_vector(:,2) + &
            rcp_drhoEps_xc(:,3) * minus_i * recip_vector(:,3)

       ! FFT back to obtain the convolution
       call fft3(second_term(:), rcp_drhoEps_xc(:,1), grid_size, +1)
       ! accumulate to potential
       do rr = 1, n_my_grid_points
          xc_potential(rr,spin) = xc_potential(rr,spin) + second_term(rr)
       end do
    end do ! spin

    deallocate(grad_density, second_term, rcp_drhoEps_xc, STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_xc_potential_GGA_PBE: Error dealloc mem")
    call reg_dealloc_mem(area_ops, grid_size*(1+3*(1+nspin)), type_dbl)

    return
  end subroutine get_xc_potential_hyb_PBE0
  !!*****


  !!****f* XC_module/get_xc_potential_GGA_PBE_obsolete *
  !!
  !!  NAME
  !!   get_xc_potential_GGA_PBE_obsolete
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the exchange-correlation potential
  !!   on the grid within GGA using the
  !!   Perdew-Burke-Ernzerhof. It also calculates the
  !!   total exchange-correlation energy.
  !!
  !!   Note that this is the functional described in
  !!   Phys. Rev. Lett. 77, 3865 (1996)
  !!
  !!   It is also (depending on an optional parameter)
  !!     either revPBE, Phys. Rev. Lett. 80, 890 (1998),
  !!     or RPBE, Phys. Rev. B 59, 7413 (1999)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   11/11/05
  !!  MODIFICATION HISTORY
  !!   17/07/06 rgradient(1,:) used for storage of final result, before transform,
  !!            instead of direct transform of the sum
  !!            rgradient(1,:) + rgradient(2,:) + rgradient(3,:)
  !!          (this was less efficient)
  !!   15:55, 27/04/2007 drb
  !!     Changed recip_vector, grad_density to (n,3) for speed
  !!   2008/11/13 ast
  !!     Added revPBE and RPBE
  !!   2011/12/12 L.Tong
  !!     Removed third, it is now defined in numbers_module
  !!   2012/04/11 L.Tong
  !!     Added optional parameters x_epsilon and c_epsilon to output
  !!     the exchange and correlation parts of xc_epsilon respectively.
  !!  SOURCE
  !!
  subroutine get_xc_potential_GGA_PBE_obsolete(density, xc_potential, &
                                               xc_epsilon, xc_energy, &
                                               size, flavour,         &
                                               x_epsilon, c_epsilon)

    use datatypes
    use numbers
    use global_module, only: &!rcellx, rcelly, rcellz,               &
                             functional_gga_pbe96_rev98,           &
                             functional_gga_pbe96_r99
    use dimens,        only: grid_point_volume, n_my_grid_points
    use GenComms,      only: gsum, cq_abort
    use fft_module,    only: fft3, recip_vector

    implicit none

    ! Passed variables
    integer,      intent(in)  :: size
    real(double), intent(in)  :: density(size)
    real(double), intent(out) :: xc_potential(size), xc_epsilon(size)
    real(double), intent(out) :: xc_energy
    ! optional
    integer,intent(in), optional :: flavour
    real(double), dimension(size), optional, intent(out) :: x_epsilon, c_epsilon

    !     Local variables
    integer :: n
    integer :: selector

    real(double), allocatable, dimension(:)   :: grad_density
    real(double), allocatable, dimension(:,:) :: grad_density_xyz
    real(double), allocatable, dimension(:)   :: ex_lda, ec_lda
    real(double) :: rho, grad_rho, rho1_3, rho1_6, ks, s, s2, t, t2,   &
                    A, At2, factor0, factor1, factor2, factor3,        &
                    factor4, denominator0, numerator1, denominator1,   &
                    num_den1, numerator2, dt2_drho, dA_drho,           &
                    dnumerator1_drho, ddenominator1_drho,              &
                    dnumerator1_dt2, dnumerator1_dA,                   &
                    ddenominator1_dt2, ddenominator1_dA, dfactor2_dt2, &
                    dfactor2_dnumerator1, dfactor2_ddenominator1,      &
                    dfactor2_drho, dfactor3_dfactor2, dfactor3_drho
    real(double) :: xc_energy_lda_total,                               &
                    e_correlation_lda,                                 &
                    de_correlation_lda,                                &
                    e_exchange, e_correlation,                         &
                    de_exchange, de_correlation,                       &
                    dde_exchange, dde_correlation
    real(double) :: df_dgrad_rho
    real(double) :: kappa, mu_kappa
    real(double) :: rpbe_exp
    ! Gradient in reciprocal space
    complex(double_cplx), allocatable, dimension(:,:) :: rgradient

    !     From Phys. Rev. Lett. 77, 3865 (1996)
    real(double), parameter :: mu = 0.21951_double
    real(double), parameter :: beta = 0.066725_double
    real(double), parameter :: gamma = 0.031091_double
    real(double), parameter :: kappa_ori = 0.804_double

    !     From Phys. Rev. Lett. 80, 890 (1998)
    real(double), parameter :: kappa_alt = 1.245_double

    !     Precalculated constants
    real(double), parameter :: mu_kappa_ori = 0.27302_double     ! mu/kappa_ori
    real(double), parameter :: mu_kappa_alt = 0.17631_double     ! mu/kappa_alt
    real(double), parameter :: two_mu = 0.43902_double           ! 2*mu
    real(double), parameter :: beta_gamma = 2.146119_double      ! beta/gamma
    real(double), parameter :: beta_X_gamma = 0.002074546_double ! beta*gamma
    real(double), parameter :: k01 = 0.16162045967_double        ! 1/(2*(3*pi*pi)**(1/3))
    real(double), parameter :: k02 = -0.16212105381_double       ! -3*mu*((4*pi/3)**(1/3))/(2*pi*alpha)
                                                                 ! =mu*k00*k01(LDA_PW92)=mu*k04
    real(double), parameter :: k03 = 1.98468639_double           ! ((4/pi)*(3*pi*pi)**(1/3))**(1/2)
    real(double), parameter :: k04 = -0.738558852965_double      ! -3*((4*pi/3)**(1/3))/(2*pi*alpha) = k00*k01 in LDA_PW92
    real(double), parameter :: k05 = 0.05240415_double           ! -2*k01*k02
    real(double), parameter :: k06 = -0.593801317784_double      ! k04*kappa_ori
    real(double), parameter :: k07 = -0.984745137287_double      ! 4*k04/3
    real(double), parameter :: seven_thirds = 2.333333333_double ! 7/3

    integer :: stat
    !      Selector options
    integer, parameter :: fx_original    = 1                     ! Used in PBE and revPBE
    integer, parameter :: fx_alternative = 2                     ! Used in RPBE

    ! Choose between PBE or revPBE parameters
    if(PRESENT(flavour)) then
      if(flavour==functional_gga_pbe96_rev98) then
        kappa=kappa_alt
        mu_kappa=mu_kappa_alt
      else
        kappa=kappa_ori
        mu_kappa=mu_kappa_ori
      end if
    else
      kappa=kappa_ori
      mu_kappa=mu_kappa_ori
    end if

    if (present(x_epsilon)) then
       allocate(ex_lda(size), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating ex_lda for PBE functional: ", stat)
    end if

    if (present(c_epsilon)) then
       allocate(ec_lda(size), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating ec_lda for PBE functional ", stat)
    end if

    allocate(grad_density(size), grad_density_xyz(size,3),&
             rgradient(size,3),STAT = stat)
    !initialisation  2010.Oct.30 TM
    grad_density(:) = zero
    grad_density_xyz(:,:) = zero
    rgradient(:,:) = zero
    xc_epsilon(:) = zero
    xc_potential(:) = zero
    if (present(x_epsilon)) then
       x_epsilon = zero
       ex_lda = zero
    end if
    if (present(c_epsilon)) then
       c_epsilon = zero
       ec_lda = zero
    end if
    if (stat /= 0) &
         call cq_abort("Error allocating arrays for PBE functional: ", stat)

    ! Choose functional form
    if (present(flavour)) then
       if (flavour == functional_gga_pbe96_r99) then
          selector = fx_alternative
       else
          selector = fx_original
       end if
    else
       selector = fx_original
    end if

    ! Build the gradient of the density
    call build_gradient(density, grad_density_xyz, size)

    grad_density(:) = sqrt(grad_density_xyz(:,1)**2 + &
                           grad_density_xyz(:,2)**2 + &
                           grad_density_xyz(:,3)**2)

    ! Get the LDA part of the functional
    if (present(x_epsilon) .and. present(c_epsilon)) then
       call get_xc_potential_LDA_PW92(density, xc_potential, xc_epsilon, &
                                      xc_energy_lda_total, size,         &
                                      x_epsilon=ex_lda,                  &
                                      c_epsilon=ec_lda)
    else if (present(x_epsilon)) then
       call get_xc_potential_LDA_PW92(density, xc_potential, xc_epsilon, &
                                      xc_energy_lda_total, size,         &
                                      x_epsilon=ex_lda)
    else if (present(c_epsilon)) then
       call get_xc_potential_LDA_PW92(density, xc_potential, xc_epsilon, &
                                      xc_energy_lda_total, size,         &
                                      c_epsilon=ec_lda)
    else
       call get_xc_potential_LDA_PW92(density, xc_potential, xc_epsilon, &
                                      xc_energy_lda_total, size)
    end if

    xc_energy = zero
    do n = 1, n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       grad_rho = grad_density(n)

       !!!!!!   XC GGA ENERGY

       ! Exchange

       if (rho > RD_ERR) then
          rho1_3 = rho ** third
          rho1_6 = sqrt (rho1_3)
          s = k01 * grad_rho / (rho ** four_thirds)
          s2 = s * s
          if(selector == fx_alternative) then           ! RPBE
            rpbe_exp = exp(-mu_kappa * s2)
            e_exchange = k06*rho1_3*(1.0_double-rpbe_exp)
          else                                          ! PBE, revPBE
            denominator0 = 1.0 / (1.0 + mu_kappa * s2)
            factor0 = k02 * rho1_3
            factor1 = s2 * denominator0
            ! NOTE: This doesn't look like in Phys. Rev. Lett. 77:18,
            !       3865 (1996) because the 1 in Fx, has been
            !       multiplied by Ex-LDA and is implicit in
            !       xc_energy_lda(n), in the total energy below
            e_exchange = factor0 * factor1
          end if
       else
          e_exchange = zero
       end if

       if (present(x_epsilon)) x_epsilon(n) = e_exchange + ex_lda(n)

       ! Correlation

       if (rho > RD_ERR) then
          e_correlation_lda = xc_epsilon(n) - k04 * rho1_3

          ! t=grad_rho/(2*rho*ks); ks=sqrt(4*kf/pi); kf=(3*pi*pi*rho)**(1/3); s=grad_rho/(2*rho*kf)
          ks = k03 * rho1_6
          t = grad_rho / ( 2 * ks * rho )
          t2 = t * t

          A = exp(-e_correlation_lda / gamma) - 1.0
          if (A > RD_ERR) then
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

       if (present(c_epsilon)) c_epsilon(n) = e_correlation + ec_lda(n)

       !*ast* TEST-POINT 1
       ! Both exchange and correlation
       xc_energy = xc_energy + (xc_epsilon(n) + e_exchange + e_correlation)*rho
       ! Only LDA part
       !xc_energy = xc_energy + (xc_epsilon(n))*rho
       ! LDA + exchange
       !xc_energy = xc_energy + (xc_epsilon(n) + e_exchange)*rho
       ! Only exchange
       !xc_energy = xc_energy + (e_exchange)*rho
       ! LDA + correlation
       !xc_energy = xc_energy + (xc_epsilon(n) + e_correlation)*rho
       ! Only correlation
       !xc_energy = xc_energy + (e_correlation)*rho

       !!!!!!   POTENTIAL

       !!!   Terms due to df/drho

       ! Exchange

       if (rho > RD_ERR) then
          if(selector == fx_alternative) then           ! RPBE
            de_exchange = k07 * rho1_3 * (kappa - (two_mu * s2 + kappa) * rpbe_exp)
          else                                          ! PBE, revPBE
            de_exchange = four_thirds * factor0 * factor1 * ( 1 - 2* denominator0 )
          end if
       else
          de_exchange = zero
       end if

       ! Correlation

       if (rho > RD_ERR) then
          de_correlation_lda = ( xc_potential(n) &
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

       !*ast* TEST-POINT 2
       xc_potential(n) = xc_potential(n) + de_exchange + de_correlation
       ! Only LDA part
       !xc_potential(n) = xc_potential(n)
       ! LDA + exchange
       !xc_potential(n) = xc_potential(n) + de_exchange
       ! Only exchange
       !xc_potential(n) = de_exchange
       ! LDA + correlation
       !xc_potential(n) = xc_potential(n) + de_correlation
       ! Only correlation
       !xc_potential(n) = de_correlation

       !!!   Terms due to df/d|grad_rho|

       ! Exchange

       if (rho > RD_ERR) then
          if(selector == fx_alternative) then           ! RPBE
             dde_exchange = -k05 * s * rpbe_exp
          else                                          ! PBE, revPBE
             dde_exchange = -k05 * s * denominator0 * denominator0!factor1 * denominator0
          end if
       else
          dde_exchange = zero
       end if

       ! Correlation

       if (rho > RD_ERR) then
          numerator2 = beta_X_gamma * t * ( 1.0 + 2.0 * At2 ) &
                     / (( gamma * denominator1 + beta * t2 * numerator1 ) * denominator1)
          dde_correlation = numerator2 / ks;
       else
          dde_correlation = zero
       end if

       ! Normalisation (modulus of gradient)

       !*ast* TEST-POINT 3
       if(abs(grad_density(n)) > RD_ERR) then  !DEBUG
       df_dgrad_rho = (dde_exchange + dde_correlation) / grad_density(n)
       else                 !DEBUG
       df_dgrad_rho = zero   !DEBUG
       endif                !DEBUG
       ! Only LDA part
       !df_dgrad_rho = 0.0
       ! LDA + exchange
       !df_dgrad_rho = dde_exchange / grad_density(n)
       ! Only exchange
       !df_dgrad_rho = dde_exchange / grad_density(n)
       ! LDA + correlation
       !df_dgrad_rho = dde_correlation / grad_density(n)
       ! Only correlation
       !df_dgrad_rho = dde_correlation / grad_density(n)


       ! Gradient times derivative of energy

       grad_density_xyz(n,1) = grad_density_xyz(n,1)*df_dgrad_rho
       grad_density_xyz(n,2) = grad_density_xyz(n,2)*df_dgrad_rho
       grad_density_xyz(n,3) = grad_density_xyz(n,3)*df_dgrad_rho

       !xc_epsilon  !2010.Oct.30 TM
       xc_epsilon(n) = xc_epsilon(n) + e_exchange + e_correlation
       !*ast* TEST-POINT 4
       !TM delta_E_xc = delta_E_xc + (xc_epsilon(n) + e_exchange + e_correlation)*rho
       ! Only LDA part
       !delta_E_xc = delta_E_xc + (xc_epsilon(n))*rho
       ! LDA + exchange
       !delta_E_xc = delta_E_xc + (xc_epsilon(n) + e_exchange)*rho
       ! Only exchange
       !delta_E_xc = delta_E_xc + (e_exchange)*rho
       ! LDA + correlation
       !delta_E_xc = delta_E_xc + (xc_epsilon(n) + e_correlation)*rho
       ! Only correlation
       !delta_E_xc = delta_E_xc + (e_correlation)*rho

    end do ! do n_my_grid_points


    !!!   Final steps of the energy calculation

    ! Add the energies and 'integrate' over the volume of the grid points

    call gsum(xc_energy)
    xc_energy = xc_energy * grid_point_volume

    ! Fourier transform the gradient, component by component

    call fft3(grad_density_xyz(:,1), rgradient(:,1), size, -1 )
    call fft3(grad_density_xyz(:,2), rgradient(:,2), size, -1 )
    call fft3(grad_density_xyz(:,3), rgradient(:,3), size, -1 )

    !!!   Get the second term of the potential by taking derivatives

    ! First, get the scalar product of the (normalised) gradient and the wave vector (times i)

    rgradient(:,1) = -rgradient(:,1)*minus_i*recip_vector(:,1)
    rgradient(:,1) = rgradient(:,1) - rgradient(:,2)*minus_i*recip_vector(:,2)
    rgradient(:,1) = rgradient(:,1) - rgradient(:,3)*minus_i*recip_vector(:,3)

    ! Add terms of the scalar product and then Fourier transform back the resultant vector
    ! NOTE: Store the result in the modulus of the gradient (not needed anymore) to save memory

    call fft3( grad_density, rgradient(:,1), size, 1 )

    ! Finally, get the potential
    ! NOTE that here grad_density is NOT the modulus of the gradient,
    !      but a term of the potential (see Fourier transform above)

    do n=1,n_my_grid_points
        xc_potential(n) = xc_potential(n) - grad_density(n)
    end do

    ! I changed the order of deallocation, because I was told that deallocation
    ! of the latest allocated array should be done first.
    ! (though I am not sure whether it is true or not)   2010.Oct.30 TM

    if(allocated(rgradient)) then
       deallocate(rgradient,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating rgradient",stat)
    end if
    if(allocated(grad_density_xyz)) then
       deallocate(grad_density_xyz,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating grad_density_xyz",stat)
    end if
    if(allocated(grad_density)) then
       deallocate(grad_density,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating grad_density",stat)
    end if
    if (allocated(ex_lda)) then
       deallocate(ex_lda,STAT=stat)
       if (stat /= 0) call cq_abort("Error deallocating ex_lda", stat)
    end if
    if (allocated(ec_lda)) then
       deallocate(ec_lda,STAT=stat)
       if (stat /= 0) call cq_abort("Error deallocating ec_lda", stat)
    end if

    return
  end subroutine get_xc_potential_GGA_PBE_obsolete
  !!***


  !!****f* XC_module/build_gradient *
  !!
  !!  NAME
  !!   build_gradient
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the gradient of the density, and its modulus
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   11/11/05
  !!  MODIFICATION HISTORY
  !!   15:55, 27/04/2007 drb
  !!    Changed recip_vector, grad_density to (n,3) for speed
  !!   2012/04/27 L.Tong
  !!   - removed function calculating the modulus of the grad_density,
  !!     since this is now redundant.
  !!   - renamed grad_density_xyz to grad_density
  !!  SOURCE
  !!
  subroutine build_gradient(density, grad_density, size)

    use datatypes
    use numbers
!    use global_module, only: rcellx, rcelly, rcellz
    use dimens,        only: n_grid_x, n_grid_y, n_grid_z
    use fft_module,    only: fft3, recip_vector
    use GenComms,      only: cq_abort

    implicit none

    ! Passed variables
    integer, intent(in) :: size
    real(double), intent(in),  dimension(size) :: density
    !ORI real(double), intent(out), dimension(3, size) :: grad_density_xyz
    real(double), intent(out), dimension(size,3) :: grad_density

    ! Local variables
    ! Density in reciprocal space
    complex(double_cplx), allocatable :: rdensity(:)
    ! Temporal reciprocal density
    complex(double_cplx), allocatable :: rdensity_tmp(:)
    integer :: stat

    !local recip_vec_tm

    allocate(rdensity(size), rdensity_tmp(size), STAT=stat)
    if(stat /= 0) &
         call cq_abort('ERROR in build_gradient : stat,size = ',stat,size)
    grad_density = zero ! to remove SIGFPE   2010.Oct.25 TM
    rdensity = zero     ! to remove SIGFPE   2010.Oct.25 TM
    rdensity_tmp = zero ! to remove SIGFPE   2010.Oct.25 TM

    ! Fourier transform the density
    call fft3(density, rdensity, size, -1)

    ! Compute the derivative with respect to x
    rdensity_tmp = -rdensity * minus_i * recip_vector(:,1)!*rcellx/n_grid_x
    call fft3(grad_density(:,1), rdensity_tmp, size, 1)

    ! Compute the derivative with respect to y
    rdensity_tmp = -rdensity * minus_i * recip_vector(:,2)!*rcelly/n_grid_y
    call fft3(grad_density(:,2), rdensity_tmp, size, 1)

    ! Compute the derivative with respect to z
    rdensity_tmp = -rdensity * minus_i * recip_vector(:,3)!*rcellz/n_grid_z
    call fft3(grad_density(:,3), rdensity_tmp, size, 1)

    return
  end subroutine build_gradient
  !!***


  !!****f* H_matrix_module/map_density *
  !!
  !!  NAME
  !!   map_density
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Map cartesian coordinates to density index
  !! (Only for one processor and old arrangement of the density)
  !! (x,y,z) starts at (1,1,1)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   21/11/05
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine map_density(x, y, z, n)

    use block_module, only: nx_in_block, ny_in_block, nz_in_block
    use group_module, only: blocks

    implicit none

    ! Passed variables
    integer, intent(in) :: x, y, z
    integer, intent(out) :: n

    ! Local variables
    integer :: gx, gy, gz, bx, by, bz, n_block

    if((x > (nx_in_block * blocks%ngcellx)).OR.&
       (y > (ny_in_block * blocks%ngcelly)).OR.&
       (z > (nz_in_block * blocks%ngcellz))) &
    then
        n = -1
        return
    end if

    gx=1+mod(x-1, nx_in_block)
    gy=1+mod(y-1, ny_in_block)
    gz=1+mod(z-1, nz_in_block)

    bx=1+(x-1)/nx_in_block
    by=1+(y-1)/ny_in_block
    bz=1+(z-1)/nz_in_block

    if(mod(by+bz, 2) == 1) then
        bx = blocks%ngcellx - bx + 1
    end if

    if(mod(bz, 2) == 0) then
        by = blocks%ngcelly - by + 1
    end if

    n_block = bx + blocks%ngcellx*((by-1) + blocks%ngcelly*(bz-1))

    n = nx_in_block*ny_in_block*nz_in_block*(n_block-1) + &
        gx + nx_in_block*((gy-1) + ny_in_block*(gz-1))

    return
  end subroutine map_density
  !!***


  !!****f* XC_module/get_dxc_potential *
  !!
  !!  NAME
  !!   get_dxc_potential
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the derivative of the exchange-correlation potential
  !!   on the grid within LDA using the Ceperley-Alder
  !!   interpolation formula. This is only for non self-consistent
  !!   calculations (Harris-Foulkes), for the forces.
  !!
  !!   Note that this is the Perdew-Zunger parameterisation of the
  !!   Ceperley-Alder results for a homogeneous electron gas, as
  !!   described in Phys. Rev. B 23, 5048 (1981), with Ceperley-Alder
  !!   in Phys. Rev. Lett. 45, 566 (1980)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   09:36, 2003/03/20 dave
  !!  MODIFICATION HISTORY
  !!   2007/11/16 12:10 dave
  !!    Bug fix: sign of ninth*(two*s+r)*rs term in dv_correlation was
  !!    wrong
  !!   2008/03/03 18:34 dave
  !!    Removed dsqrt
  !!   2011/12/13 L.Tong
  !!    Removed third, as it is now defined in numbers module
  !!  SOURCE
  !!
  subroutine get_dxc_potential(density, dxc_potential, size)

    use datatypes
    use numbers
    use dimens,   only: grid_point_volume, n_my_grid_points
    use GenComms, only: gsum

    implicit none

    ! Passed variables
    integer size
    real(double), dimension(size) :: density, dxc_potential

    ! Local variables
    integer :: n
    real(double) :: denominator, e_correlation, e_exchange, ln_rs, &
                    vnumerator, rcp_rs, rho, rs, rs_ln_rs, sq_rs,  &
                    dv_correlation, dv_exchange, dfirst, dsecond
    real(double), parameter :: alpha  = -0.45817_double
    real(double), parameter :: beta_1 =  1.0529_double
    real(double), parameter :: beta_2 =  0.3334_double
    real(double), parameter :: gamma  = -0.1423_double
    real(double), parameter :: p =  0.0311_double
    real(double), parameter :: q = -0.048_double
    real(double), parameter :: r =  0.0020_double
    real(double), parameter :: s = -0.0116_double
    real(double), parameter :: ninth = 1.0_double/9.0_double
    real(double), parameter :: sixth = 1.0_double/6.0_double

    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       if (rho>RD_ERR) then ! Find radius of hole
          rcp_rs = ( four_thirds * pi * rho )**(third)
       else
          rcp_rs = zero
       end if
       ! First, the easy part: exchange
       e_exchange = alpha * rcp_rs
       if(rho>RD_ERR) then
          dv_exchange = four*ninth * e_exchange/rho
       else
          dv_exchange = zero
       end if
       if (rcp_rs>zero) then
          rs = one/rcp_rs
       else
          rs = zero
       end if
       sq_rs = sqrt(rs)
       if (rs>=one) then
          denominator = one / (one + beta_1 * sq_rs + beta_2 * rs)
          vnumerator = one + seven_sixths * beta_1 * sq_rs +  &
               four_thirds * beta_2 * rs
          dfirst = (one + beta_1 * sq_rs + beta_2 * rs)*(-seven_thirtysixths*beta_1*sq_rs - four*ninth*beta_2*rs)
          dsecond = -two*vnumerator*(-sixth*beta_1*sq_rs-third*beta_2*rs)
          if(rho>RD_ERR) then
             dv_correlation = gamma * (dfirst+dsecond) * denominator * denominator * denominator / rho
          else
             dv_correlation = zero
          end if
       else if ((rs<one).and.(rs>RD_ERR)) then
          ln_rs = log(rs)
          rs_ln_rs = rs * ln_rs
          if(rho>RD_ERR) then
             ! DRB 2007/11/16 Changed sign of ninth to PLUS to correct error
             dv_correlation = -(third*p+two*ninth*r*rs_ln_rs+ninth*(two*s+r)*rs)/rho
          else
             dv_correlation = zero
          end if
       else
          dv_correlation = zero
       end if
       ! Both X and C
       dxc_potential(n) = (dv_exchange + dv_correlation)
       ! Just C
       !dxc_potential(n) = dv_correlation
       ! Just X
       !dxc_potential(n) = dv_exchange
    end do ! do n_my_grid_points
    return
  end subroutine get_dxc_potential
  !!***


  !!****f* force_module/get_GTH_dxc_potential *
  !!
  !!  NAME
  !!   get_GTH_dxc_potential
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the derivative of the exchange-correlation potential
  !!   on the grid within LDA using the Ceperley-Alder
  !!   interpolation formula. This is only for non self-consistent
  !!   calculations (Harris-Foulkes), for the forces.
  !!
  !!   Note that this is the Goedecker/Teter/Hutter parameterisation -
  !!   see PRB 54, 1703 (1996)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:54, 25/03/2003 drb
  !!  MODIFICATION HISTORY
  !!   2011/12/13 L.Tong
  !!    Removed third, it is now defined in numbers module
  !!  SOURCE
  !!
  subroutine get_GTH_dxc_potential(density, dxc_potential, size)

    use datatypes
    use numbers
    use dimens,   only: grid_point_volume, n_my_grid_points
    use GenComms, only: gsum

    implicit none

    ! Passed variables
    integer size

    real(double) :: density(size), dxc_potential(size)

    !     Local variables
    integer      :: n
    real(double) :: denominator, e_correlation, e_exchange, ln_rs,     &
                    vnumerator, rcp_rs, rho, rs, rs_ln_rs,t1, t2, dt1, &
                    dt2, d2t1, d2t2, sq_rs, dv_correlation,            &
                    dv_exchange, dfirst, dsecond, drs_dRho, d2rs_dRho2
    real(double), parameter :: a0=0.4581652932831429_double
    real(double), parameter :: a1=2.217058676663745_double
    real(double), parameter :: a2=0.7405551735357053_double
    real(double), parameter :: a3=0.01968227878617998_double
    real(double), parameter :: b1=1.000000000000000_double
    real(double), parameter :: b2=4.504130959426697_double
    real(double), parameter :: b3=1.110667363742916_double
    real(double), parameter :: b4=0.02359291751427506_double

    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       if (rho>RD_ERR) then ! Find radius of hole
          rcp_rs = ( 4.0_double*third * pi * rho )**(third)
          rs = one/rcp_rs
       else
          rcp_rs = zero
          rs = zero
          rho = zero
       end if
       if(rs>zero) then
          drs_dRho = -rs / (3.0_double * rho)
          d2rs_dRho2 = 4.0_double * rs / (9.0_double * rho * rho)
          t1 = a0 + rs * (a1 + rs * (a2 + rs * a3))
          t2 = rs * (b1 + rs * (b2 + rs * (b3 + rs * b4)))
          dt1 = a1 + rs * (2.0_double * a2 + rs * 3.0_double * a3)
          dt2 = b1 + rs * (2.0_double * b2 + rs * (3.0_double * b3 + &
                           rs * 4.0_double * b4))
          d2t1 = 2.0_double * a2 + 6.0_double * a3 * rs
          d2t2 = 2.0_double * b2 + rs * (6.0_double * b3 + rs * 12.0_double * b4)
          dxc_potential(n) = 2.0_double*drs_dRho *                &
                             (-dt1 / t2 + t1 * dt2 / (t2 * t2)) + &
                             rho *                                &
                             (drs_dRho * drs_dRho *               &
                             (2.0_double * dt1 * dt2 / (t2*t2) -  &
                              d2t1 / t2 -                         &
                              2.0_double * t1 * dt2 * dt2 /       &
                              (t2*t2*t2) + d2t2 * t1/(t2*t2)) +   &
                             d2rs_dRho2*(-dt1 / t2 + t1 * dt2 / (t2 * t2)))
       else
          dxc_potential(n) = zero
       end if
    end do ! do n_my_grid_points
    return
  end subroutine get_GTH_dxc_potential
  !!***


  !!****f* H_matrix_module/get_dxc_potential_LDA_PW92 *
  !!
  !!  NAME
  !!   get_dxc_potential_LDA_PW92
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the derivative of the exchange-correlation potential
  !!   on the grid within LDA using the Ceperley-Alder
  !!   interpolation formula. This is only for non self-consistent
  !!   calculations (Harris-Foulkes), for the forces.
  !!
  !!   Note that this is the Perdew-Wang parameterisation of the
  !!   Ceperley-Alder results for a homogeneous electron gas, as
  !!   described in Phys. Rev. B 45, 13244 (1992), with Ceperley-Alder
  !!   in Phys. Rev. Lett. 45, 566 (1980)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   30/01/06
  !!  MODIFICATION HISTORY
  !!   2008/03/03 18:34 dave
  !!    Removed dsqrt
  !!   2011/10/17 L.Tong
  !!    Corrected (changed) 1 to one in log (1 + 1/denominator)
  !!   2011/12/13 L.Tong
  !     Removed third, it is now defined in numbers module
  !!  SOURCE
  !!
  subroutine get_dxc_potential_LDA_PW92(density, dxc_potential, size, &
                                        eclda, declda_drho, d2eclda_drho2)

    use datatypes
    use numbers
    use dimens,   only: grid_point_volume, n_my_grid_points
    use GenComms, only: gsum

    implicit none

    ! Passed variables
    integer :: size

    real(double), intent(in)    :: density(size)
    real(double), intent(inout) :: dxc_potential(size)
    real(double), optional   :: eclda(size)
    real(double), optional   :: declda_drho(size)
    real(double), optional   :: d2eclda_drho2(size)

    !     Local variables
    integer n
    real(double) :: prefactor, postfactor, denominator, &
                    e_correlation, e_exchange,          &
                    rcp_rs, rho, rs, sq_rs,             &
                    v_correlation, v_exchange,          &
                    delta_prefactor, delta_postfactor,  &
                    dv_exchange, dv_correlation,        &
                    denominator2, dnumrho_drho,         &
                    dden_drho, rho4_3,                  &
                    dprefactor_drho, d2prefactor_drho2, &
                    dposfactor_drho, d2posfactor_drho2, &
                    dec_drho, d2ec_drho2


    !     From Table I, Phys. Rev. B 45, 13244 (1992), for reference
    real(double), parameter :: alpha  = 1.0421234_double
    real(double), parameter :: alpha1 = 0.21370_double
    real(double), parameter :: beta1  = 7.5957_double
    real(double), parameter :: beta2  = 3.5876_double
    real(double), parameter :: beta3  = 1.6382_double
    real(double), parameter :: beta4  = 0.49294_double
    real(double), parameter :: A      = 0.031091_double

    !     Precalculated constants
    real(double), parameter :: k00 =  1.611991954_double     ! (4*pi/3)**(1/3)
    real(double), parameter :: k01 = -0.458165347_double     ! -3/(2*pi*alpha)
    real(double), parameter :: k02 = -0.062182_double        ! -2*A
    real(double), parameter :: k03 = -0.0132882934_double    ! -2*A*alpha1
    real(double), parameter :: k04 =  0.4723158174_double    ! 2*A*beta1
    real(double), parameter :: k05 =  0.2230841432_double    ! 2*A*beta2
    real(double), parameter :: k06 =  0.1018665524_double    ! 2*A*beta3
    real(double), parameter :: k07 =  0.03065199508_double   ! 2*A*beta4
    real(double), parameter :: k08 = -0.008858862267_double  ! 2*k03/3
    real(double), parameter :: k09 =  0.0787193029_double    ! k04/6
    real(double), parameter :: k10 =  0.074361381067_double  ! k05/3
    real(double), parameter :: k11 =  0.0509332762_double    ! k06/2
    real(double), parameter :: k12 =  0.0204346633867_double ! 2*k07/3
    real(double), parameter :: four_ninths       = 4.0_double/9.0_double
    real(double), parameter :: minus_four_thirds = -4.0_double/3.0_double
    real(double), parameter :: k13 = 0.0027477997778_double ! (k08-k03)/k00

    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)

       !!!!   EXCHANGE
       if (rho > RD_ERR) then ! Find radius of hole
          rcp_rs = k00 * ( rho**third )
          dv_exchange = four_ninths * k01 * rcp_rs / rho    ! 4*e_exchange/(9*rho)
       else
          rcp_rs = zero
          dv_exchange = zero
       end if

       e_exchange = k01*rcp_rs

       !!!!   CORRELATION
       if (rcp_rs > zero) then
          rs = one/rcp_rs
       else
          rs = zero
       end if
       sq_rs = sqrt(rs)

       prefactor = k02 + k03*rs
       denominator = sq_rs * ( k04 + sq_rs * ( k05 + sq_rs * ( k06 + k07 * sq_rs)))
       if (denominator > zero) then
          postfactor = log (one + one/denominator)
       else
          postfactor = zero
       end if

       ! Return correlation energy if requested
       ! NOTE: This "energy" is not the actual integrand, but the integrand divided by rho
       if(present(eclda)) then
         eclda(n) = prefactor*postfactor
       end if


       ! NOTE: delta_prefactor is actually the derivative of rho*prefactor
       !       delta_postfactor is rho times the derivative of postfactor
       delta_prefactor  = k02 + k08*rs
       if (sq_rs > zero) then
          denominator2 = ( denominator * ( 1 + denominator ) )
          delta_postfactor = sq_rs * ( k09 + sq_rs*(k10 + sq_rs*( k11 + k12 * sq_rs ))) &
                           / denominator2

          dnumrho_drho = -sq_rs *( 7.0*k09/6.0 + sq_rs * ( 4.0*k10/3.0 + sq_rs * ( 3.0*k11/2.0 +sq_rs*5.0*k12/3.0)))
          dden_drho = -sq_rs*(k04/6.0 +sq_rs*(third*k05 + sq_rs*(half*k06 + sq_rs*2.0*k07/3.0)))
       else
          delta_postfactor = 0

          dnumrho_drho = zero
          dden_drho = zero
       end if

       if(rho > RD_ERR) then
          rho4_3 = rho**minus_four_thirds
          dprefactor_drho   = k13*rho4_3
          d2prefactor_drho2 = minus_four_thirds*dprefactor_drho/rho

          dposfactor_drho   = delta_postfactor/rho
          d2posfactor_drho2 = (dnumrho_drho &
                            - delta_postfactor * ( 1 + 2*denominator) * dden_drho)/(denominator2 * rho *rho)
       else
          rho4_3 = zero
          dprefactor_drho   = zero
          d2prefactor_drho2 = zero

          dposfactor_drho   = zero
          d2posfactor_drho2 = zero
       end if

       dec_drho = postfactor*dprefactor_drho + prefactor*dposfactor_drho

       ! Return derivative of correlation energy (see note above) if requested
       if(present(declda_drho)) then
          declda_drho(n) = dec_drho
       end if


       d2ec_drho2 = 2*dposfactor_drho*dprefactor_drho &
                  + postfactor*d2prefactor_drho2 + prefactor*d2posfactor_drho2

       ! Return second derivative of correlation energy (see note above) if requested
       if(present(d2eclda_drho2)) then
          d2eclda_drho2(n) = d2ec_drho2
       end if

       dv_correlation = 2*dec_drho + rho*d2ec_drho2

       dxc_potential(n) = dv_exchange + dv_correlation
    end do ! do n_my_grid_points

    return
  end subroutine get_dxc_potential_LDA_PW92
  !!***


  !!****f* H_matrix_module/get_dxc_potential_GGA_PBE *
  !!
  !!  NAME
  !!   get_dxc_potential_GGA_PBE
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the potential needed for the non-self-consistent
  !!   within GGA using the Perdew-Burke-Ernzerhof.
  !!
  !!   Note that this is the functional described in
  !!   Phys. Rev. Lett. 77, 3865 (1996)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A.S. Torralba
  !!  CREATION DATE
  !!   31/01/05
  !!  MODIFICATION HISTORY
  !!   15:54, 27/04/2007 drb
  !!    Changed recip_vector, grad_density and tmp2, tmp3 to (n,3) for speed
  !!   2008/11/13 ast
  !!    Added new PBE functional types
  !!   2011/12/13 L.Tong
  !!    Removed third, and eight, they are now defined in numbers module
  !!  SOURCE
  !!
  subroutine get_dxc_potential_GGA_PBE(density, density_out, &
                                       dxc_potential, size, flavour )

    use datatypes
    use numbers
    use global_module, only: & !rcellx, rcelly, rcellz, &
                             functional_gga_pbe96_rev98, functional_gga_pbe96_r99
    use dimens,        only: grid_point_volume, n_my_grid_points
    use GenComms,      only: gsum
    use fft_module,    only: fft3, recip_vector
!    use energy,        only: delta_E_xc

    implicit none

    ! Passed variables
    integer size
    real(double), intent(in)      :: density(size)
    real(double), intent(in)      :: density_out(size)
    real(double), intent(inout)   :: dxc_potential(size)
    integer, optional, intent(in) :: flavour

    ! Local variables
    integer n, i
    integer selector

    real(double) :: grad_density(size), grad_density_xyz(size,3)
    real(double) :: rho, grad_rho, rho1_3, rho1_6
    real(double) :: e_exchange
    real(double) :: xc_energy_lda_total, dxc_potential_lda(size), &
                    eclda(size), declda_drho(size), d2eclda_drho2(size)
    real(double) :: de_dgrad(size)
    real(double) :: s, s2, factor0, factor1, denominator0, denominator0_2
    real(double) :: d2e_drho2, d2e_dgrad2(size), d2e_dgrad_drho(size)
    real(double) :: dden0_drho, df0f1_drho, ds_grad, dden0_dgrad, ds_dgrad, df1_dgrad
    real(double) :: diff_rho(size)
    complex(double_cplx) :: tmp1(size), tmp2(size,3)
    real(double) :: tmp3(size,3), tmp_factor

    real(double) :: a, a2, t, t2, t3, t4, num1, den1, den1_2, fl
    real(double) :: factor_exp, factor_exp2, factor_exp3
    real(double) :: dt_drho, da_drho, dnum1_drho, dden1_drho, dfl_drho, dec_drho
    real(double) :: dt_dgrad, dnum1_dgrad, dden1_dgrad, dfl_dgrad
    real(double) :: d2t_drho2, d2a_drho2, d2num1_drho2, d2den1_drho2, d2fl_drho2
    real(double) :: d2num1_dgrad2, d2den1_dgrad2, d2fl_dgrad2
    real(double) :: d2t_drho_dgrad, d2num1_drho_dgrad, d2den1_drho_dgrad
    real(double) :: d2fl_drho_dgrad

    real(double) :: kappa, mu_kappa

    !     From Phys. Rev. Lett. 77, 3865 (1996)
    real(double), parameter :: mu = 0.21951_double
    real(double), parameter :: beta = 0.066725_double
    real(double), parameter :: gamma = 0.031091_double
    real(double), parameter :: kappa_ori = 0.804_double

    !     From Phys. Rev. Lett. 80, 890 (1998)
    real(double), parameter :: kappa_alt = 1.245_double

    !     Precalculated constants
    real(double), parameter :: mu_kappa_ori = 0.27302_double     ! mu/kappa_ori
    real(double), parameter :: mu_kappa_alt = 0.17631_double     ! mu/kappa_alt
    real(double), parameter :: beta_gamma = 2.146119_double      ! beta/gamma
    real(double), parameter :: beta_X_gamma = 0.002074546_double ! beta*gamma
    real(double), parameter :: k01 = 0.16162045967_double        ! 1/(2*(3*pi*pi)**(1/3))
    real(double), parameter :: k02 = -0.16212105381_double       ! -3*mu*((4*pi/3)**3)/(2*pi*alpha)
                                                                 ! = mu*k00*k01(LDA_PW92)=mu*k04
    real(double), parameter :: k03 = 1.98468639_double           ! ((4/pi)*(3*pi*pi)**(1/3))**(1/2)
    real(double), parameter :: k04 = -0.738558852965_double      ! -3*((4*pi/3)**3)/(2*pi*alpha) = k00*k01 in LDA_PW92
    real(double), parameter :: k05 = 0.05240415_double           ! -2*k01*k02
    real(double), parameter :: eight_thirds = 8.0_double / 3.0_double

    !      Selector options
    integer, parameter :: fx_original    = 1                     ! Used in PBE and revPBE
    integer, parameter :: fx_alternative = 2                     ! Used in RPBE

    ! Choose between PBE or revPBE parameters
    if(PRESENT(flavour)) then
      if(flavour==functional_gga_pbe96_rev98) then
        kappa=kappa_alt
        mu_kappa=mu_kappa_alt
      else
        kappa=kappa_ori
        mu_kappa=mu_kappa_ori
      end if
    else
      kappa=kappa_ori
      mu_kappa=mu_kappa_ori
    end if

!*ast* AT PRESENT IGNORED
    ! Choose functional form
    if(PRESENT(flavour)) then
      if(flavour==functional_gga_pbe96_r99) then
        selector=fx_alternative
      else
        selector=fx_original
      end if
    else
      selector=fx_original
    end if

!*ast* AT PRESENT, NOT IMPLEMENTED
if(selector == fx_alternative) then
  print *,"!!!!!!!!WARNING!!!!!!!!!!!"
  print *,"WARNING: TO DO !!!!!!!!! RPBE NSC Forces NOT IMPLEMENTED YET. The values will be for standard PBE"
  print *,"!!!!!!!!WARNING!!!!!!!!!!!"
end if

    ! Build the gradient of the density
    call build_gradient(density, grad_density_xyz, size)

    grad_density(:) = sqrt(grad_density_xyz(:,1)**2 + &
                           grad_density_xyz(:,2)**2 + &
                           grad_density_xyz(:,3)**2)

    ! Get the LDA part of the functional

    call get_dxc_potential_LDA_PW92(density, dxc_potential_lda, size, &
                                    eclda, declda_drho, d2eclda_drho2)

    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       grad_rho = grad_density(n)

       diff_rho(n) = density(n) - density_out(n)

       !!!!   EXCHANGE

       !!   Energy factors

       if (rho > RD_ERR) then
          rho1_3 = rho ** third
          rho1_6 = sqrt (rho1_3)
          s = k01 * grad_rho / (rho ** four_thirds)
          s2 = s * s
          denominator0 = one / (one + mu_kappa * s2)
          denominator0_2 = denominator0 * denominator0
          factor0 = k02 * rho1_3
          factor1 = s2 * denominator0
          ! NOTE: This doesn't look like in Phys. Rev. Lett. 77:18, 3865 (1996)
          !       because the 1 in Fx, has been multiplied by Ex-LDA and is implicit
          !       in xc_energy_lda(n), in the total energy below
          e_exchange = factor0 * factor1
       else
          s  = 0.0_double
          s2 = 0.0_double
          denominator0 = 0.0_double
          denominator0_2 = 0.0_double
          factor0 = 0.0_double
          factor1 = 0.0_double
          e_exchange = zero
       end if

       !!   Second derivative of Ex wrt rho
       if (rho > RD_ERR) then
          dden0_drho = mu_kappa * eight_thirds * denominator0_2 * s2 /rho
          df0f1_drho = third * e_exchange * (one - eight * denominator0) / rho
          d2e_drho2 = four_thirds * (df0f1_drho * (one - two * denominator0) &
                    - two * e_exchange * dden0_drho)
       else
          d2e_drho2 = 0.0_double
       end if

       !!   Second derivative of Ex wrt grad
       ds_dgrad = 0.0_double
       if(rho > RD_ERR) ds_dgrad = k01 * rho**(-four_thirds)
       dden0_dgrad = -two * mu_kappa * denominator0_2 * s * ds_dgrad

       d2e_dgrad2(n) = two * rho * factor0 * ds_dgrad * (ds_dgrad * denominator0_2 &
                     + two * s * denominator0 * dden0_dgrad)

       !!   Second derivative of Ex wrt grad and rho
       df1_dgrad = two * s * denominator0_2 * ds_dgrad

       d2e_dgrad_drho(n) = (four_thirds) * factor0 &
                         * (df1_dgrad * (one - two * denominator0) &
                         - two * factor1 * dden0_dgrad)


       !!   First derivative wrt gradient
       de_dgrad(n) = rho * factor0 * df1_dgrad


       !!!!   CORRELATION

       !!   Energy terms
       factor_exp = exp(-eclda(n)/gamma)
       a = factor_exp - one
       if(a > RD_ERR) then
         a = beta_gamma / a
       else
         a = beta_gamma * BIG
       end if
       a2 = a * a

       if(rho > RD_ERR) then
         t = grad_rho/(two*k03*(rho**seven_sixths))
       else
         t = 0.0_double
       end if
       t2 = t * t
       t3 = t2 * t
       t4 = t2 * t2
       num1 = one + a*t2;
       den1 = num1 + a2*t4;
       den1_2 = den1*den1;
       if(den1 > RD_ERR) then
         fl = one + beta_gamma * t2 * num1 / den1;
       else
         fl = one
       end if


       !!   First derivative wrt rho

       if(rho > RD_ERR) then
         dt_drho = -seven_sixths*t/rho
       else
         dt_drho = 0.0_double
       end if
    !TM 2007_10_18
    ! eclda(n) = 0 -> factor_exp = one -> zero**(-2) !!!
    !ORI   factor_exp2 = (factor_exp - one)**(-two)
     if(abs(factor_exp-one) > RD_ERR) then
       factor_exp2 = one/(factor_exp - one)
       factor_exp2 = factor_exp2**2
     else
       factor_exp2 = 0.0_double
     end if

       da_drho = beta_gamma * factor_exp * factor_exp2 * declda_drho(n)/gamma
       dnum1_drho = da_drho * t2 + two*a*t*dt_drho
       dden1_drho = dnum1_drho + two*a*t4*da_drho + four*a2*t3*dt_drho
       if(den1 > RD_ERR) then
         dfl_drho = beta_gamma*((two*t*num1*dt_drho + t2*dnum1_drho)/den1 - t2*num1*dden1_drho/den1_2)
       else
         dfl_drho = 0.0_double
       end if
       dec_drho = gamma * (log(fl) + rho * dfl_drho / fl)

       !!   Second derivative with respect to rho

       if(rho > RD_ERR) then
         d2t_drho2 = seven_sixths*(t/(rho*rho) - dt_drho/rho)
       else
         d2t_drho2 = 0.0_double
       end if
       factor_exp3 = factor_exp * factor_exp2
      if(abs(factor_exp-one) > RD_ERR) then
       d2a_drho2 = (beta_gamma/gamma)*( ( ( (two * factor_exp / (factor_exp - one)) - one )&
                 * factor_exp3 * declda_drho(n)**2) /gamma &
                 + factor_exp3 * d2eclda_drho2(n))
      else
       d2a_drho2 = zero
      end if

       d2num1_drho2 = d2a_drho2*t2 + four*t*da_drho*dt_drho &
                    + two*a*((dt_drho**2) + t*d2t_drho2)
       d2den1_drho2 = d2num1_drho2 + two*((da_drho*t2)**2) &
                    + two*a*d2a_drho2*t4 &
                    + 16.0_double*a*t3*da_drho*dt_drho &
                    + 12.0_double*((a*t*dt_drho)**2) &
                    + four*a2*t3*d2t_drho2
       if(den1 > RD_ERR) then
         d2fl_drho2 = beta_gamma*(two*num1*(dt_drho**2) &
                    + t*(four*(dnum1_drho*dt_drho &
                    - num1*dden1_drho*dt_drho/den1) &
                    + two*num1*d2t_drho2) &
                    + t2*(d2num1_drho2 &
                    + two*(num1*((dden1_drho/den1)**2) &
                    - dnum1_drho*dden1_drho/den1) &
                    - num1*d2den1_drho2/den1))/den1
       else
         d2fl_drho2 = 0.0_double
       end if
       d2e_drho2 = d2e_drho2 + gamma*(two*dfl_drho/fl - rho*((dfl_drho/fl)**2) + rho*d2fl_drho2/fl);


       !!   First derivative with respect to grad

       if(rho > RD_ERR) then
         dt_dgrad = one/(two*k03*(rho**seven_sixths))
       else
         dt_dgrad = 0.0_double
       end if
       dnum1_dgrad = two*a*t*dt_dgrad
       dden1_dgrad = two*a*t*dt_dgrad + four*a2*t3*dt_dgrad
       if(den1 > RD_ERR) then
          dfl_dgrad = two*beta_gamma*t*num1*dt_dgrad/den1 + beta_gamma*t2*dnum1_dgrad/den1  &
                 - beta_gamma*t2*num1*dden1_dgrad/den1_2
       else
          dfl_dgrad = 0.0_double
       end if
       de_dgrad(n) = de_dgrad(n) + gamma * rho * dfl_dgrad / fl


       !!   Second derivative with respect to grad

       d2num1_dgrad2 = two*a*(dt_dgrad**2)
       d2den1_dgrad2 = d2num1_dgrad2 + 12.0*((a*t*dt_dgrad)**2)
       if(den1 > RD_ERR) then
        d2fl_dgrad2 = beta_gamma*(two*num1*(dt_dgrad**2)/den1 &
                    + four*t*dnum1_dgrad*dt_dgrad/den1 &
                    - two*t*num1*dden1_dgrad*dt_dgrad/den1_2 &
                    + t2*d2num1_dgrad2/den1 &
                    - t2*dden1_dgrad*dnum1_dgrad/den1_2 &
                    - two*t*num1*dden1_dgrad*dt_dgrad/den1_2 &
                    - t2*dnum1_dgrad*dden1_dgrad/den1_2 &
                    + two*t2*num1*(dden1_dgrad**2)/(den1_2*den1) &
                    - t2*num1*d2den1_dgrad2/den1_2)
       else
        d2fl_dgrad2 = zero
       end if
       d2e_dgrad2(n) = d2e_dgrad2(n) + gamma*rho*(d2fl_dgrad2/fl - (dfl_dgrad**2)/(fl*fl))

       !!   Second derivative with respect to rho and grad

       ! if((grad_rho > RD_ERR ) .or. (rho > RD_ERR)) then  ! Changed by TM 19Oct2007
       if((grad_rho > RD_ERR ) .and. (rho > RD_ERR)) then
         d2t_drho_dgrad = -seven_sixths*t/(grad_rho*rho)
       else
         d2t_drho_dgrad = 0.0_double
       end if
       d2num1_drho_dgrad = two*(t*da_drho*dt_dgrad + a*dt_drho*dt_dgrad + a*t*d2t_drho_dgrad)
       d2den1_drho_dgrad = d2num1_drho_dgrad + 8.0*a*t3*da_drho*dt_dgrad + 12.0*a2*t2*dt_drho*dt_dgrad &
                         + four*a2*t3*d2t_drho_dgrad
       if(den1 > RD_ERR) then
          d2fl_drho_dgrad = beta_gamma*(two*(dt_drho*dt_dgrad*num1 &
                          + t*(dnum1_drho*dt_dgrad &
                          + num1*(d2t_drho_dgrad &
                          - (dden1_drho*dt_dgrad &
                          + dt_drho*dden1_dgrad)/den1) &
                          + dnum1_dgrad*dt_drho &
                         + t*num1*dden1_drho*dden1_dgrad/den1_2)) &
                         + t2*(d2num1_drho_dgrad &
                         - (dnum1_dgrad*dden1_drho &
                         + dnum1_drho*dden1_dgrad &
                         + num1*d2den1_drho_dgrad)/den1))/den1
       else
         d2fl_drho_dgrad = 0.0_double
       end if
        d2e_dgrad_drho(n) = d2e_dgrad_drho(n) + gamma*(dfl_dgrad + rho*(d2fl_drho_dgrad - dfl_drho*dfl_dgrad/fl))/fl


       !!   Add term L1 to the potential
       dxc_potential(n) = dxc_potential(n) &
                        + diff_rho(n) * (d2e_drho2 + dxc_potential_lda(n))

    end do ! do n_my_grid_points


    ! Fourier transform the difference of densities
     tmp1(:)=cmplx(zero,zero,double_cplx)  !TM
    call fft3(diff_rho, tmp1, size, -1)

    !do n=1, n_my_grid_points  ! debugged 25Oct2007 TM
    do n=1, size
       ! Product by reciprocal vector stored for later use
       tmp2(n,1) = -minus_i*recip_vector(n,1)*tmp1(n)
       tmp2(n,2) = -minus_i*recip_vector(n,2)*tmp1(n)
       tmp2(n,3) = -minus_i*recip_vector(n,3)*tmp1(n)
    end do

    ! Fourier transform the vector back to the grid
    call fft3(tmp3(:,1), tmp2(:,1), size, 1)
    call fft3(tmp3(:,2), tmp2(:,2), size, 1)
    call fft3(tmp3(:,3), tmp2(:,3), size, 1)

    ! Add term L3 to potential
    do n=1, n_my_grid_points
       do i=1,3
          if(grad_density(n) > RD_ERR) then
            dxc_potential(n) = dxc_potential(n) + (tmp3(n,i) &
                                                * d2e_dgrad_drho(n) &
                                                * grad_density_xyz(n,i) )/grad_density(n)
          end if
       end do
    end do

    ! Term L4
    do n=1, n_my_grid_points
       if(grad_density(n) > RD_ERR) then
         tmp_factor =(tmp3(n,1) * grad_density_xyz(n,1) &
                    + tmp3(n,2) * grad_density_xyz(n,2) &
                    + tmp3(n,3) * grad_density_xyz(n,3)) &
                    * (d2e_dgrad2(n) &
                    - de_dgrad(n)/grad_density(n) ) / (grad_density(n)*grad_density(n))
       else
         tmp_factor = zero
       end if
       ! Reuse tmp3
       do i=1,3
          if(grad_density(n) > RD_ERR) then
            tmp3(n,i) = tmp_factor * grad_density_xyz(n,i) &
                      + de_dgrad(n) * tmp3(n,i)/grad_density(n)
          else
            tmp3(n,i) = zero
          end if
       end do
    end do

    ! Terms L2 and L5 (using L4)
    do n=1, n_my_grid_points
       do i=1,3
          if(grad_density(n) > RD_ERR) then
            tmp3(n,i) = tmp3(n,i) &
                      + diff_rho(n) * d2e_dgrad_drho(n)* grad_density_xyz(n,i) &
                      / grad_density(n)
          else
            tmp3(n,i) = zero
          end if
       end do
    end do

    tmp2(:,:) = cmplx(zero,zero,double_cplx) ! 25Oct2007 TM
    call fft3(tmp3(:,1), tmp2(:,1), size, -1)
    call fft3(tmp3(:,2), tmp2(:,2), size, -1)
    call fft3(tmp3(:,3), tmp2(:,3), size, -1)

    !do n=1, n_my_grid_points  ! debugged 25Oct2007 TM
    do n=1, size
       ! Product by reciprocal vector stored for later use
       tmp1(n) = -minus_i &
               *(recip_vector(n,1)*tmp2(n,1) &
               + recip_vector(n,2)*tmp2(n,2) &
               + recip_vector(n,3)*tmp2(n,3))
    end do

    ! Use first component of tmp3 to store final vector
    call fft3(tmp3(:,1), tmp1, size, 1)

    do n=1, n_my_grid_points
       dxc_potential(n) = dxc_potential(n) - tmp3(n,1)
    end do

    return
  end subroutine get_dxc_potential_GGA_PBE
  !!***


  !!****f* XC_module/get_dxc_potential_LSDA_PW92
  !! PURPOSE
  !!   Calculates the derivative of the exchange-correlation potential
  !!   on the grid within LSDA using the Ceperley-Alder interpolation
  !!   formula. This is only for non self-consistent calcuations
  !!   (Harris-Foulkes), for the forces.
  !!
  !!   The functional is given by Phys. Rev. B 45, 13244
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2011/10/18
  !! MODIFICATION HISTORY
  !!   2012/03/26 L.Tong
  !!   - Changed spin implementation
  !!   - density is now array of (maxngrid,nspin)
  !!   - dxc_potential_ddensity is now array of
  !!   - (maxngrid,nspin,nspin), so
  !!     dxc_potential_ddensity(n,spin1,spin2) corresponds to
  !!     dVxc(spin1) / drho(spin2) at grid point n.
  !! SOURCE
  !!
  subroutine get_dxc_potential_LSDA_PW92(density,                &
                                         dxc_potential_ddensity, &
                                         size)
    use datatypes
    use numbers
    use dimens,        only: grid_point_volume, n_my_grid_points
    use GenComms,      only: gsum
    use global_module, only: nspin, spin_factor

    implicit none

    ! passed variables
    integer :: size
    real(double), dimension(:,:)   :: density
    real(double), dimension(:,:,:) :: dxc_potential_ddensity

    ! local variables
    integer      :: n, spin, spin_2
    real(double) :: rho_tot, rs, rcp_rs, zeta, sq_rs, e_c0, e_c1,      &
                    malpha_c, f, de_c0_drs, drs_drho_tot, de_c1_drs,   &
                    dmalpha_c_drs, df_dzeta, d2e_c0_drs2, d2e_c1_drs2, &
                    d2malpha_c_drs2, d2f_dzeta2, de_c_drs, de_c_dzeta, &
                    d2e_c_drs2, d2e_c_drs_dzeta, d2e_c_dzeta2, factor
    real(double), dimension(nspin)       :: rho, drs_drho, dzeta_drho
    real(double), dimension(nspin,nspin) :: dVx_drho, dVc_drho

    ! tabulated parameters
    real(double) :: alpha1, beta1, beta2, beta3, beta4, p, A

    ! precalculated parameters
    real(double), parameter :: k00 =  1.611991954_double ! (4*pi/3)**(1/3)
    real(double), parameter :: k01 = -0.413566994_double ! -(2/pi)**(1/3)/3**(2/3)
    real(double), parameter :: K02 =  1.923661051_double ! 1 / (2**(4/3)-2)
    real(double), parameter :: K03 =  0.584822362_double ! 1 / f''(0)
    real(double), parameter :: K04 = -0.328248341_double ! -(1/pi)**(1/3)/3**(2/3)

    ! loop over grid points
    do n = 1, n_my_grid_points
       rho_tot = zero
       do spin = 1, nspin
          rho(spin) = density(n,spin)
          rho_tot = rho_tot + spin_factor * rho(spin)
       end do
       if (rho_tot > RD_ERR) then
          rcp_rs = k00 * rho_tot**third
          zeta = (rho(1) - rho(nspin)) / rho_tot
       else
          rcp_rs = zero
          ! limit rho(spin) --> 0+ for all spin gives zeta --> zero
          ! (remember rho(spin) >= 0)
          zeta = zero
       end if

       ! exchange (worked out from Mathematica)
       do spin = 1, nspin
          if (rho(spin) > RD_ERR) then
             ! if spin non-polarised, dVx_drho is calculated to be dVx/drho_tot
             ! this is found to be half of dVx_drho(spin=1)
             dVx_drho(spin,spin) = &
                  k01 / (spin_factor * (rho(spin))**two_thirds)
          else
             dVx_drho(spin,spin) = zero
          end if
       end do
       if (nspin == 2) then
          dVx_drho(1,2) = zero ! dVx_drho(1,2) = dVx_up_drho_dn
          dVx_drho(2,1) = zero ! dVx_drho(2,1) = dVx_dn_drho_up
       end if

       ! correlation
       if (rcp_rs > RD_ERR) then
          rs = one / rcp_rs
       else
          ! here rs really should be infty, however, it is a TRICK to
          ! set it to zero here, but it gives the correct end results
          ! for G(rs), G'(rs) and G''(rs), etc.
          rs = zero
       end if
       sq_rs = sqrt(rs)

       ! drs_drho
       if (nspin == 1) then
          if (rho_tot > RD_ERR) then
             drs_drho_tot = - third * rs / rho_tot
          else
             drs_drho_tot = zero
          end if
       else if (nspin == 2) then
          do spin = 1, nspin
             if (rho(spin) > RD_ERR) then
                drs_drho(spin) = - third * rs / rho_tot
                dzeta_drho(spin) = ((-one)**(spin+1) - zeta) / rho_tot
             else
                drs_drho(spin) = zero
                dzeta_drho(spin) = zero
             end if
          end do
       end if

       ! e_c0 and de_c0_drs
       ! from table 1 of Phys. Rev. B 45, 13244
       call PW92_G(rs,                    &
                   A=0.031091_double,     &
                   alpha1=0.21370_double, &
                   beta1=7.5957_double,   &
                   beta2=3.5876_double,   &
                   beta3=1.6382_double,   &
                   beta4=0.49294_double,  &
                   p=one,                 &
                   G=e_c0,                &
                   dG_drs=de_c0_drs,      &
                   d2G_drs2=d2e_c0_drs2)

       ! for spin non-polarised calculations this is enough for
       ! de_c_drs and d2e_c_drs2
       de_c_drs = de_c0_drs
       d2e_c_drs2 = d2e_c0_drs2

       if (nspin == 2) then
          ! e_c1 and de_c1_drs
          ! from table 1 of Phys. Rev. B 45, 13244
          call PW92_G(rs,                    &
                      A=0.015545_double,     &
                      alpha1=0.20548_double, &
                      beta1=14.1189_double,  &
                      beta2=6.1977_double,   &
                      beta3=3.3662_double,   &
                      beta4=0.62517_double,  &
                      p=one,                 &
                      G=e_c1,                &
                      dG_drs=de_c1_drs,      &
                      d2G_drs2=d2e_c1_drs2)
          ! malpha_c and dmalpha_c_drs
          ! from table 1 of Phys. Rev. B 45, 13244
          call PW92_G(rs,                    &
                      A=0.016887_double,     &
                      alpha1=0.11125_double, &
                      beta1=10.357_double,   &
                      beta2=3.6231_double,   &
                      beta3=0.88026_double,  &
                      beta4=0.49671_double,  &
                      p=one,                 &
                      G=malpha_c,            &
                      dG_drs=dmalpha_c_drs,  &
                      d2G_drs2=d2malpha_c_drs2)

          ! f and df_dzeta
          f = K02 * ((one + zeta)**four_thirds + (one - zeta)**four_thirds - two)
          df_dzeta = four_thirds * K02 * ((one + zeta)**one_third - &
                     (one - zeta)**one_third)
          d2f_dzeta2 = one_third * four_thirds * K02 * &
                       ((one + zeta)**(-two_thirds) +  &
                        (one - zeta)**(-two_thirds))

          ! first order derivatives
          de_c_drs = de_c_drs - de_c0_drs * f * zeta**4 + &
                     de_c1_drs * f * zeta**4 + &
                     dmalpha_c_drs * K03 * f * (zeta**4 - one)
          de_c_dzeta = four * zeta**3 * f * (e_c1 - e_c0 + K03 * malpha_c) + &
                       df_dzeta * (zeta**4 * e_c1 - zeta**4 * e_c0 + &
                                   (zeta**4 - one) * K03 * malpha_c)

          ! second order derivatives
          ! d2e_c_drs2
          d2e_c_drs2 = d2e_c_drs2 - d2e_c0_drs2 * f * zeta**4 + &
                       d2e_c1_drs2 * f * zeta**4 + &
                       d2malpha_c_drs2 * K03 * f * (zeta**4 - one)

          ! d2e_c_drs_dzeta
          d2e_c_drs_dzeta = four * zeta**3 * f * &
                            (de_c1_drs - de_c0_drs + K03 * dmalpha_c_drs) + &
                            df_dzeta * (zeta**4 * de_c1_drs -  &
                                        zeta**4 * de_c0_drs + &
                                        (zeta**4 - one) * K03 * dmalpha_c_drs)

          ! d2e_c_dzeta2
          d2e_c_dzeta2 = &
               (twelve * zeta * zeta * f + eight * zeta**3 * df_dzeta) * &
               (e_c1 - e_c0 + K03 * malpha_c) + &
               d2f_dzeta2 * (zeta**4 * e_c1 - zeta**4 * e_c0 + &
                             (zeta**4 - one) * K03 * malpha_c)
       end if

       ! finally get the derivatives of the correlation potentials
       if (nspin == 1) then
          ! for spin non-polarised calculations, we need to calculate
          ! dVc_drho_tot
          dVc_drho(1,1) = (two_thirds * de_c_drs - &
                           one_third * rs * d2e_c_drs2) * &
                          drs_drho_tot
       else
          do spin = 1, nspin
             do spin_2 = 1, nspin
                factor = zeta + (-one)**spin
                dVc_drho(spin,spin_2) = &
                     (two_thirds * de_c_drs - &
                      one_third * rs * d2e_c_drs2 - &
                      factor * d2e_c_drs_dzeta) * &
                     drs_drho(spin_2) - &
                     (one_third * rs * d2e_c_drs_dzeta + &
                      factor * d2e_c_dzeta2) * &
                     dzeta_drho(spin_2)
             end do
          end do
       end if

       ! collect things together
       do spin = 1, nspin
          do spin_2 = 1, nspin
             dxc_potential_ddensity(n,spin,spin_2) = &
                  dVx_drho(spin,spin_2) + dVc_drho(spin,spin_2)
          end do
       end do

    end do ! n = 1, n_my_grid_points

    return
  end subroutine get_dxc_potential_LSDA_PW92
  !!*****


  !!****f* XC_module/PW92_G
  !! PURPOSE
  !!   Calculates the interpolation function G defined in PW92 LSDA
  !!   paper:
  !!   Perdew and Wang, Phys. Rev. B, 1992, 45, 13244-13249
  !! USAGE
  !!   call PW92_G(rs, A, alpha1, beta1, beta2, beta3, beta4, p,
  !!               G, dG_drs, d2G_drs2)
  !! INPUTS
  !!   real(double) rs     : the variable for G(rs) function
  !!   real(double) A      : interpolation parameter
  !!   real(double) alpha1 : interpolation parameter
  !!   real(double) beta1  : interpolation parameter
  !!   real(double) beta2  : interpolation parameter
  !!   real(double) beta3  : interpolation parameter
  !!   real(double) beta4  : interpolation parameter
  !!   real(double) p      : interpolation parameter
  !! OUTPUT
  !!   real(double) G        : G(rs)
  !!   real(double) dG_drs   : dG/drs      (optional)
  !!   real(double) d2G_drs2 : d^2G/drs^2  (optional)
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2013/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine PW92_G(rs, A, alpha1, beta1, beta2, beta3, beta4, p, &
                    G, dG_drs, d2G_drs2)
    use datatypes
    use numbers
    implicit none

    ! passed parameters
    real(double), intent(in) :: rs, A, alpha1, beta1, beta2, beta3, beta4, p
    real(double), intent(out) :: G
    real(double), intent(out), optional :: dG_drs, d2G_drs2
    ! local variables
    real(double) :: sq_rs, Q0, Q1, dQ0_drs, dQ1_drs, d2Q0_drs2, d2Q1_drs2

    if (abs(rs) <= RD_ERR) then
       G = zero
       if (present(dG_drs)) dG_drs = zero
       if (present(d2G_drs2)) d2G_drs2 = zero
       return
    end if

    sq_rs = sqrt(rs)

    Q0 = -two * A * (one + alpha1 * rs)
    Q1 = two * A * (beta1 * sq_rs + beta2 * rs + &
                    beta3 * sq_rs**3 + beta4 * rs**(p + one))
    G = Q0 * log(one + one / Q1)
    if (present(dG_drs) .or. present(d2G_drs2)) then
       dQ0_drs = -two * A * alpha1
       dQ1_drs = A * (beta1 / sq_rs + two * beta2 + &
                      three * beta3 * sq_rs + four * beta4 * rs)
       d2Q1_drs2 = half * A * (-beta1 / (sq_rs**3) + &
                               three * beta3 / sq_rs + &
                               eight * beta4)
    end if
    if (present(dG_drs)) then
       dG_drs = dQ0_drs * log(one + one / Q1) - &
                (Q0 * dQ1_drs) / (Q1 * (Q1 + one))
    end if
    if (present(d2G_drs2)) then
       d2G_drs2 = - (two * dQ0_drs * dQ1_drs + Q0 * d2Q1_drs2) / &
                  (Q1 * (Q1 + one)) + &
                  dQ1_drs * (two * Q1 + one) * Q0 * dQ1_drs / &
                  ((Q1 * (Q1 + one))**2)
    end if

    return
  end subroutine PW92_G
  !*****

end module XC_module
