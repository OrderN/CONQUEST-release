! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module md_control
! ------------------------------------------------------------------------------
! Code area 7: move atoms
! ------------------------------------------------------------------------------

!!****h* Conquest/md_control
!!  NAME
!!   md_control
!!  PURPOSE
!!   Controls temperature and pressure in molecular dynamics. Algorithm
!!   adapted from G. Martyna et al. Mol. Phys. 5, 1117 (1996) ("MTTK Paper")
!!  USES
!!   datatypes, numbers, move_atoms, global_module
!!  AUTHOR
!!   Zamaan Raza 
!!  CREATION DATE
!!   2017/10/24 09:24
!!  MODIFICATION HISTORY
!!   2017/10/24 zamaan
!!    Refactored original NVT code from NVTMD branch
!!   2017/11/09 zamaan
!!    Implemented isotropic MTTK barostat
!!   2017/12/05 zamaan
!!    Added separate NHC chain for cell degrees of freedom
!!   2017/05/29 zamaan
!!    Corrected sign of potential energy contribution in get_thermostat_energy. Added
!!    cell degrees of freedom (cell_ndof) to equations for clarity.
!!   2018/07/12 zamaan
!!    Added SSM integrator, adapted from W. Shinoda et al., PRB 69:134103
!!    (2004)
!!   2018/08/11 zamaan
!!    Added iprint_MD level to all write (Conquest_out) statements to tidy
!!    up output
!!   2018/08/12 zamaan
!!    Removed read/write_thermo/baro_checkpoint, replaced with a unified &
!!    checkpoint that includes ionic velocities
!!   2019/04/23 zamaan
!!    Implemented stochastic velocity scaling MD, Bussi et al., J. Chem.
!!    Phys. 126, 014101 (2007) (NVT) and Bussi et al., J. Chem. Phys.
!!    130, 074101 (2009) (NPT)
!!  SOURCE
!!
module md_control

  use datatypes
  use numbers
  use global_module,    only: ni_in_cell, io_lun, iprint_MD, &
                              temp_ion, flag_MDcontinue, flag_MDdebug
  use move_atoms,       only: fac_Kelvin2Hartree
  use species_module,   only: species, mass
  use GenComms,         only: inode, ionode, cq_abort
  use input_module,     only: leqi

  implicit none

  ! Flags

  character(20) :: md_cell_constraint

  ! Unit conversion factors
  real(double), parameter :: fac_HaBohr32GPa = 29421.02648438959
  real(double), parameter :: fac_fs2atu = 41.3413745758
  real(double), parameter :: fac_invcm2hartree = 4.5563352812122295E-6
  real(double), parameter :: fac_thz2hartree = 0.0001519828500716

  ! Files
  character(20) :: md_position_file = 'md.position'
  character(20) :: md_check_file = "md.checkpoint"
  character(20) :: md_thermo_file = "md.thermostat"
  character(20) :: md_baro_file = "md.barostat"
  character(20) :: md_trajectory_file = "trajectory.xsf"
  character(20) :: md_frames_file = "Frames"
  character(20) :: md_stats_file = "Stats"

  ! Module variables
  character(20) :: md_thermo_type, md_baro_type
  real(double)  :: md_tau_T, md_tau_P, md_target_press, md_bulkmod_est, &
                   md_box_mass, md_ndof_ions, md_omega_t, md_omega_p, &
                   md_tau_T_equil, md_tau_P_equil, md_p_drag, md_t_drag
  integer       :: md_n_nhc, md_n_ys, md_n_mts, md_berendsen_equil
  logical       :: flag_write_xsf, md_cell_nhc, md_calc_xlmass
  logical       :: flag_extended_system = .false.
  real(double), dimension(3,3), target      :: lattice_vec
  real(double), dimension(:), allocatable   :: md_nhc_mass, md_nhc_cell_mass
  real(double), dimension(:,:), allocatable, target :: ion_velocity

  !!****s* md_control/type_thermostat
  !!  NAME
  !!   type_thermostat 
  !!  PURPOSE
  !!   Container for all thermostat-related variables
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 09:26
  !!  SOURCE
  !!  
  type type_thermostat
    ! General thermostat variables
    character(20)       :: thermo_type  ! thermostat type
    character(20)       :: baro_type
    character(3)        :: ensemble
    real(double)        :: T_int        ! instantateous temperature
    real(double)        :: T_ext        ! target temperature
    real(double)        :: ke_ions      ! kinetic energy of ions
    real(double)        :: dt           ! time step
    integer             :: ndof         ! number of degrees of freedom
    integer             :: cell_ndof    ! cell degrees of freedom
    logical             :: append
    logical             :: cell_nhc     ! separate NHC for cell?
    real(double)        :: lambda       ! velocity scaling factor

    ! Weak coupling thermostat variables
    real(double)        :: tau_T        ! temperature coupling time period

    ! Nose-Hoover chain thermostat variables
    integer             :: n_nhc        ! number of Nose-Hoover heat baths
    integer             :: n_ys         ! Yoshida-Suzuki order
    integer             :: n_mts_nhc    ! number of time steps for NHC
    real(double)        :: e_thermostat ! energy of thermostat
    real(double)        :: e_barostat   ! energy of barostat (for NPT SVR)
    real(double)        :: e_nhc_ion    ! energy of ionic NHC thermostats
    real(double)        :: e_nhc_cell   ! energy of cell NHC thermostats
    real(double)        :: ke_nhc_ion   ! ke of ionic NHC thermostats
    real(double)        :: ke_nhc_cell  ! ke of cell NHC thermostats
    real(double)        :: pe_nhc_ion   ! pe of ionic NHC thermostats
    real(double)        :: pe_nhc_cell  ! pe of cell NHC thermostats
    real(double)        :: ke_target    ! for computing NHC force
    real(double)        :: t_drag       ! drag factor for thermostat
    character(40)       :: nhc_fmt      ! format string for printing NHC arrays
    character(40)       :: nhc_fmt2
    real(double), dimension(:), allocatable :: eta    ! thermostat "position"
    real(double), dimension(:), allocatable :: v_eta  ! thermostat "velocity"
    real(double), dimension(:), allocatable :: G_nhc  ! "force" on thermostat
    real(double), dimension(:), allocatable :: m_nhc  ! thermostat mass
    real(double), dimension(:), allocatable :: eta_cell  ! NHC for cell
    real(double), dimension(:), allocatable :: v_eta_cell
    real(double), dimension(:), allocatable :: G_nhc_cell
    real(double), dimension(:), allocatable :: m_nhc_cell
    real(double), dimension(:), allocatable :: dt_ys  ! Yoshida-Suzuki time steps

    contains

      procedure, public   :: init_thermo
      procedure, public   :: init_nhc
      procedure, public   :: init_ys
      procedure, public   :: get_berendsen_thermo_sf
      procedure, public   :: get_svr_thermo_sf
      procedure, public   :: v_rescale
      procedure, public   :: get_thermostat_energy
      procedure, public   :: get_temperature_and_ke
      procedure, public   :: dump_thermo_state

      procedure, private  :: update_G_nhc
      procedure, private  :: propagate_eta
      procedure, private  :: propagate_v_eta_lin
      procedure, private  :: propagate_v_eta_exp
      procedure, private  :: apply_nhc_drag

      procedure, public   :: integrate_nhc
  end type type_thermostat
!!***

  !!****s* md_control/type_barostat
  !!  NAME
  !!   type_barostat
  !!  PURPOSE
  !!   Container for all barostat-related variables
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/08 11:45
  !!  SOURCE
  !!  
  type type_barostat
    ! General barostat variables
    character(20)       :: baro_type    ! thermostat type
    real(double)        :: P_int        ! instantateous pressure
    real(double)        :: P_ext        ! target pressure
    real(double)        :: PV           ! pressure*volume enthalpy term
    real(double)        :: volume, volume_ref, v_old
    real(double)        :: ke_ions      ! kinetic energy of ions
    real(double)        :: dt           ! time step
    integer             :: ndof         ! number of degrees of freedom
    integer             :: cell_ndof    ! cell degrees of freedom
    logical             :: append
    real(double), dimension(3,3)  :: lat       ! lattice vectors
    real(double), dimension(3,3)  :: lat_ref   ! reference lattice vectors
    real(double), dimension(3,3)  :: ke_stress      ! kinetic contrib to stress
    real(double), dimension(3,3)  :: static_stress  ! static contrib to stress
    real(double), dimension(3,3)  :: total_stress   ! total stress

    ! constants for polynomial expansion
    real(double)        :: c2, c4, c6, c8

    ! Weak coupling barostat variables
    real(double)        :: tau_P        ! pressure coupling time period
    real(double)        :: bulkmod      ! estimated bulk modulus

    ! Extended Lagrangian barostat variables
    real(double)        :: box_mass
    real(double)        :: e_barostat
    real(double)        :: p_drag       ! drag factor for barostat

    ! Isotropic variables
    real(double)        :: odnf
    real(double)        :: eps, eps_ref ! 1/3 log(V/V_0)
    real(double)        :: v_eps        ! box velocity
    real(double)        :: G_eps        ! box force
    real(double)        :: mu           ! box scaling factor

    ! Fully flexible cell variables
    real(double), dimension(3,3)  :: h    ! bookkeeping variable for ortho-ssm
    real(double), dimension(3,3)  :: v_h  ! hdot = v_h x h
    real(double), dimension(3,3)  :: G_h
    real(double), dimension(3,3)  :: c_g
    real(double), dimension(3,3)  :: I_e
    real(double), dimension(3,3)  :: I_s
    real(double), dimension(3,3)  :: ident
    real(double), dimension(3,3)  :: onfm
    real(double), dimension(3)    :: lambda

    contains

      procedure, public   :: init_baro
      procedure, public   :: get_pressure_and_stress
      procedure, public   :: get_volume
      procedure, public   :: get_berendsen_baro_sf
      procedure, public   :: propagate_berendsen
      procedure, public   :: get_barostat_energy
      procedure, public   :: propagate_npt_mttk
      procedure, public   :: propagate_r_mttk
      procedure, public   :: propagate_box_mttk
      procedure, public   :: dump_baro_state
      procedure, public   :: update_cell
      procedure, public   :: scale_box_velocity

      procedure, private  :: update_G_box
      procedure, private  :: propagate_eps_lin
      procedure, private  :: propagate_eps_exp
      procedure, private  :: propagate_v_box_lin
      procedure, private  :: propagate_v_box_exp
      procedure, private  :: apply_box_drag
      procedure, private  :: update_vscale_fac
      procedure, private  :: poly_sinhx_x

      procedure, public   :: integrate_box
      procedure, public   :: couple_box_particle_velocity
      procedure, public   :: propagate_box_ssm

  end type type_barostat
!!***

contains

  !!****m* md_control/init_thermo *
  !!  NAME
  !!   init_thermo
  !!  PURPOSE
  !!   initialise thermostat
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2018/05/16 16:45
  !!  SOURCE
  !!  
  subroutine init_thermo(th, thermo_type, baro_type, dt, ndof, tau_T, ke_ions)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    character(*), intent(in)              :: thermo_type, baro_type
    real(double), intent(in)              :: dt
    integer, intent(in)                   :: ndof
    real(double), intent(in)              :: ke_ions, tau_T

    ! local variables
    real(double)                          :: dummy_rn
    integer                               :: i

    if (inode==ionode .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') 'Welcome to init_thermo'
    th%T_int = two*ke_ions/fac_Kelvin2Hartree/md_ndof_ions
    th%T_ext = temp_ion
    th%dt = dt
    th%ndof = ndof
    th%ke_ions = ke_ions
    th%ke_target = half*md_ndof_ions*fac_Kelvin2Hartree*th%T_ext
    th%tau_T = tau_T
    th%thermo_type = thermo_type
    th%baro_type = baro_type
    th%cell_nhc = md_cell_nhc
    th%e_thermostat = zero

    if (leqi(baro_type, 'none')) then
      md_cell_constraint = 'fixed'
      th%cell_nhc = .false.
    end if
    select case(md_cell_constraint)
    case('fixed')
      th%cell_ndof = 0
    case('volume')
      th%cell_ndof = 1
    case('xyz')
      th%cell_ndof = 3
    case default
      call cq_abort('MD.CellConstraint must be "fixed", "volume" or "xyz"')
    end select

    ! For NPT stochastic velocity rescaling we also need to thermostat the box
    if (leqi(th%baro_type, 'svr')) then
      th%ke_target = th%ke_target +  &
                     half*th%cell_ndof*fac_Kelvin2Hartree*th%T_ext
    end if

    select case(th%thermo_type)
    case('none')
    case('berendsen')
      flag_extended_system = .false.
    case('svr')
      ! This will be true for isobaric-isothermal SVR, because we will use
      ! Parrinello-Rahman NPH dynamics with a SVR thermostat (barostat 
      ! is initialised AFTER thermostat)
      flag_extended_system = .false.
    case('nhc')
      call th%init_nhc(dt)
    case default
      call cq_abort("MD.ThermoType must be 'none', 'berendsen', 'svr' or 'nhc'")
    end select

    if (inode==ionode .and. iprint_MD > 1) then
      write(io_lun,'(4x,"Thermostat type: ",a)') th%thermo_type
      write(io_lun,'(4x,"Extended system: ",l)') flag_extended_system
      if (.not. leqi(thermo_type, 'none')) then
        write(io_lun,'(4x,a,f10.2)') 'Target temperature        T_ext = ', &
                                     th%T_ext
        write(io_lun,'(4x,a,f10.2)') 'Coupling time period      tau_T = ', &
                                     th%tau_T
      end if
    end if

  end subroutine init_thermo
  !!***

  !!****m* md_control/init_nhc *
  !!  NAME
  !!   init_nhc
  !!  PURPOSE
  !!   initialise Nose-Hoover chain variables
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 10:30
  !!  SOURCE
  !!  
  subroutine init_nhc(th, dt)

    use memory_module,    only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module,    only: area_moveatoms

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    real(double), intent(in)              :: dt

    ! local variables
    character(40)                         :: fmt1, fmt2
    real(double)                          :: omega_thermo, omega_baro, &
                                             ndof_baro, tauT, tauP
    integer                               :: i

    flag_extended_system = .true.
    th%n_nhc = md_n_nhc
    th%n_ys = md_n_ys
    th%n_mts_nhc = md_n_mts
    th%lambda = one
    th%t_drag = one - (md_t_drag*dt/md_tau_T/md_n_mts/md_n_ys)
    write(th%nhc_fmt,'("(4x,a12,",i4,"e14.6)")') th%n_nhc
    write(th%nhc_fmt2,'("(",i4,"e20.12)")') th%n_nhc

    allocate(th%eta(th%n_nhc))
    allocate(th%v_eta(th%n_nhc))
    allocate(th%G_nhc(th%n_nhc))
    allocate(th%m_nhc(th%n_nhc))
    allocate(th%eta_cell(th%n_nhc))    ! independent thermostat for cell DOFs
    allocate(th%v_eta_cell(th%n_nhc))
    allocate(th%G_nhc_cell(th%n_nhc))
    allocate(th%m_nhc_cell(th%n_nhc))
    call reg_alloc_mem(area_moveatoms, 8*th%n_nhc, type_dbl)

    ! Calculate the masses for extended lagrangian variables?
    if (md_calc_xlmass) then
      ! Guess for thermostat/barostat time constants if not specified during
      ! initial_read
      if (md_tau_T < zero) md_tau_T = dt*ten
      if (md_tau_P < zero) md_tau_P = dt*ten*ten
      th%tau_T = md_tau_T
      ! Seems to work better with the factor twopi (angular frequency)
      omega_thermo = twopi/md_tau_T
      omega_baro = twopi/md_tau_P

      th%m_nhc(1) = md_ndof_ions*th%T_ext*fac_Kelvin2Hartree/omega_thermo**2
      if (th%cell_nhc) then
        select case (md_cell_constraint)
        case('volume')
          ndof_baro = one
        case('xyz')
          ndof_baro = three
        end select
        th%m_nhc_cell(1) = (ndof_baro**2)*th%T_ext*fac_Kelvin2Hartree/omega_baro**2
      end if
      do i=2,th%n_nhc
        th%m_nhc(i) = th%T_ext*fac_Kelvin2Hartree/omega_thermo**2
        if (th%cell_nhc) th%m_nhc_cell(i) = th%T_ext*fac_Kelvin2Hartree/omega_baro**2
      end do
    else
      th%m_nhc = md_nhc_mass
      th%m_nhc_cell = md_nhc_cell_mass
    end if

    ! Defaults for heat bath positions, velocities, masses
    th%append = .false.
    if (flag_MDcontinue) th%append = .true.
    ! initialise thermostat velocities and forces
    th%eta = zero
    th%v_eta = sqrt(two*th%T_ext*fac_Kelvin2Hartree/th%m_nhc(1)) 
!      th%v_eta = zero
    th%G_nhc = zero
    do i=2,th%n_nhc
      th%G_nhc(i) = (th%m_nhc(i-1)*th%v_eta(i-1)**2 - &
                    & th%T_ext*fac_Kelvin2Hartree)/th%m_nhc(i)
    end do
    th%eta_cell = zero
    th%v_eta_cell = zero
    th%G_nhc_cell = zero
    if (th%cell_nhc) then
      th%v_eta_cell = sqrt(two*th%T_ext*fac_Kelvin2Hartree/th%m_nhc(1)) 
      do i=2,th%n_nhc
        th%G_nhc_cell(i) = (th%m_nhc_cell(i-1)*th%v_eta_cell(i-1)**2 - &
                           & th%T_ext*fac_Kelvin2Hartree)/th%m_nhc_cell(i)
      end do
    end if

    ! Yoshida-Suzuki time steps
    call th%init_ys(dt, th%n_ys)

    write(fmt1,'("(4x,a16,",i4,"f12.6)")') th%n_nhc
    write(fmt2,'("(4x,a16,",i4,"f12.2)")') th%n_ys
    if (inode==ionode .and. iprint_MD > 1) then
      write(io_lun,'(2x,a)') 'Welcome to init_nhc'
      write(io_lun,'(4x,a,i10)')   'Number of NHC thermostats n_nhc = ', &
                                   th%n_nhc
      write(io_lun,'(4x,a,f15.8)')   'NHC velocity drag factor t_drag = ', &
                                   th%t_drag
      write(io_lun,'(4x,a,i10)')   'Multiple time step order  n_mts = ', &
                                   th%n_mts_nhc
      write(io_lun,'(4x,a,i10)')   'Yoshida-Suzuki order      n_ys  = ', &
                                   th%n_ys
      write(io_lun,fmt2) 'YS time steps:  ', th%dt_ys
      write(io_lun,fmt1) 'ion  NHC masses:', th%m_nhc
      if (th%cell_nhc) write(io_lun,fmt1) 'cell NHC masses:', th%m_nhc_cell
    end if
    call th%get_thermostat_energy

  end subroutine init_nhc
  !!***

  !!****m* md_control/init_ys *
  !!  NAME
  !!   init_ys
  !!  PURPOSE
  !!   initialise Yoshida-Suzuki time steps
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/12/08 16:24
  !!  SOURCE
  !!  
  subroutine init_ys(th, dt, n_ys)

    use memory_module,    only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module,    only: area_moveatoms

    ! passed variables
    class(type_thermostat), intent(inout)     :: th
    real(double), intent(in)                  :: dt
    integer, intent(in)                       :: n_ys

    ! local variables
    integer                                   :: i, j, k, l, m
    real(double), dimension(:,:), allocatable :: psuz
    real(double)                              :: xnt

    allocate(psuz(n_ys,5))
    allocate(th%dt_ys(n_ys))
    call reg_alloc_mem(area_moveatoms, n_ys, type_dbl)

    do i=2,n_ys
      xnt = one/(two*real(i,double)-one)
      do j=1,5
        if (mod(j,3) == 0) then
          psuz(i,j) = one - four/(four-four**xnt)
        else
          psuz(i,j) = one/(four-four**xnt)
        end if
      end do
    end do

    ! Yoshida-Suzuki time steps
    k = 0
    select case(th%n_ys) ! The YS order
    case(1)
      th%dt_ys(1) = one
    case(3)
      th%dt_ys(1) = one/(two - two**(one/three))
      th%dt_ys(2) = one - two*th%dt_ys(1)
      th%dt_ys(3) = th%dt_ys(1)
    case(5)
      do i=1,5
        k = k+1
        th%dt_ys(k) = psuz(2,i)
      end do
    case(7)
      th%dt_ys(1) = -1.17767998417887_double
      th%dt_ys(2) = 0.235573213359357_double
      th%dt_ys(3) = 0.784513610477560_double
      th%dt_ys(4) = one - two*(th%dt_ys(1)+th%dt_ys(2)+th%dt_ys(3))
      th%dt_ys(5) = th%dt_ys(3)
      th%dt_ys(6) = th%dt_ys(2)
      th%dt_ys(7) = th%dt_ys(1)
    case(15)
      th%dt_ys(1) = 0.914844246229740_double
      th%dt_ys(2) = 0.253693336566229_double
      th%dt_ys(3) = -1.44485223686048_double
      th%dt_ys(4) = -0.158240635368243_double
      th%dt_ys(5) = 1.93813913762276_double
      th%dt_ys(6) = -1.96061023297549_double
      th%dt_ys(7) = 0.102799849391985_double
      th%dt_ys(8) = one - two*(th%dt_ys(1)+th%dt_ys(2)+th%dt_ys(3)+&
                               th%dt_ys(4)+th%dt_ys(5)+th%dt_ys(6)+&
                               th%dt_ys(7))
      th%dt_ys(9) = th%dt_ys(7)
      th%dt_ys(10) = th%dt_ys(6)
      th%dt_ys(11) = th%dt_ys(5)
      th%dt_ys(12) = th%dt_ys(4)
      th%dt_ys(13) = th%dt_ys(3)
      th%dt_ys(14) = th%dt_ys(2)
      th%dt_ys(15) = th%dt_ys(1)
    case(25)
      do j=1,5
        do i=1,5
          k = k+1
          th%dt_ys(k) = psuz(2,i)*psuz(3,j)
        end do
      end do
    case(125)
      do l=1,5
        do j=1,5
          do i=1,5
            k = k+1
            th%dt_ys(k) = psuz(2,i)*psuz(3,j)*psuz(4,l)
          end do
        end do
      end do
    case(625)
      do m=1,5
        do l=1,5
          do j=1,5
            do i=1,5
              k = k+1
              th%dt_ys(k) = psuz(2,i)*psuz(3,j)*psuz(4,l)*psuz(5,m)
            end do
          end do
        end do
      end do
    case default
      call cq_abort("Invalid Yoshida-Suzuki order")
    end select
    th%dt_ys = dt*th%dt_ys/real(th%n_mts_nhc, double)
    deallocate(psuz)

  end subroutine init_ys
  !!***

  !!****m* md_control/get_temperature_and_ke *
  !!  NAME
  !!   get_temperature_and_ke
  !!  PURPOSE
  !!   Compute the kinetic energy and kinetic stress
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2018/07/23 10:12
  !!  MODIFICATION HISTORY
  !!   2018/08/11 zamaan
  !!   Added final_call optional argument so that we only print to
  !!   Conquest_out once per step
  !!  SOURCE
  !!  
  subroutine get_temperature_and_ke(th, baro, v, KE, final_call)

    use move_atoms,       only: fac

    ! Passed variables
    class(type_thermostat), intent(inout)     :: th
    type(type_barostat), intent(inout)        :: baro
    real(double), dimension(3,ni_in_cell), intent(in)  :: v  ! ion velocities
    real(double), intent(out)                 :: KE
    integer, optional                         :: final_call

    ! local variables
    integer                                   :: j, k, iatom
    real(double)                              :: m, trace

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "get_temperature_and_ke"

    ! Update the kinetic energy and stress
    baro%ke_stress = zero
    KE = zero
    do iatom=1,ni_in_cell
      m = mass(species(iatom))*fac
      do j=1,3
        ! Kinetic energy
        KE = KE + m*v(j,iatom)**2
        ! Kinetic contribution to stress tensor
        do k=1,3
          baro%ke_stress(j,k) = baro%ke_stress(j,k) + m*v(j,iatom)*v(k,iatom)
        end do
      end do
    end do 
    th%T_int = KE/fac_Kelvin2Hartree/md_ndof_ions
    KE = half*KE
    th%ke_ions = KE
!    trace = baro%ke_stress(1,1) + baro%ke_stress(2,2) + baro%ke_stress(3,3)

    if (present(final_call)) then
      if (inode==ionode .and. iprint_MD > 1) then
        write(io_lun,'(4x,"KE: ",f12.6," Ha")') KE
  !      write(io_lun,'(4x,"Tr(ke_stress): ",f12.6," Ha")') trace/three
        write(io_lun,'(4x,"T:  ",f12.6," K")') th%T_int
        if (flag_MDdebug .and. iprint_MD > 3) then
          write(io_lun,'(4x,a,3e16.8)') "ke_stress:     ", baro%ke_stress(:,1)
          write(io_lun,'(4x,a,3e16.8)') "               ", baro%ke_stress(:,2)
          write(io_lun,'(4x,a,3e16.8)') "               ", baro%ke_stress(:,3)
          write(io_lun,*)
        end if
      end if
    end if
  end subroutine get_temperature_and_ke

  !!****m* md_control/get_berendsen_thermo_sf *
  !!  NAME
  !!   get_berendsen_thermo_sf
  !!  PURPOSE
  !!   Get velocity scaling factor for Berendsen thermostat
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 14:26
  !!  SOURCE
  !!  
  subroutine get_berendsen_thermo_sf(th, dt)

    ! passed variables
    class(type_thermostat), intent(inout)   :: th
    real(double), intent(in)                :: dt

    th%lambda = sqrt(one + (dt/th%tau_T)* &
                     (th%T_ext/th%T_int - one))

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) then
      write(io_lun,'(2x,a)') "get_berendsen_thermo_sf"
    end if
  end subroutine get_berendsen_thermo_sf
  !!***

  !!****m* md_control/get_svr_thermo_sf *
  !!  NAME
  !!   get_svr_thermo_sf
  !!  PURPOSE
  !!   Get velocity scaling factor for CSVR thermostat: &
  !!   Bussi et al. J. Chem. Phys. 126, 014101 (2007)
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2019/04/19
  !!  MODIFICATION HISTORY
  !!   2019/05/21 zamaan
  !!    Replaced all calls to old RNG with new
  !!  SOURCE
  !!  
  subroutine get_svr_thermo_sf(th, dt, baro)

    use GenComms,      only: gsum, my_barrier
    use rng,           only: type_rng

    ! Passed variables
    class(type_thermostat), intent(inout)   :: th
    real(double), intent(in)                :: dt
    type(type_barostat), intent(inout), optional :: baro

    ! Local variables
    type(type_rng)  :: myrng
    real(double)    :: ke, rn, alpha_sq, exp_dt_tau, r1, sum_ri_sq, temp_fac
    integer         :: i, ndof

    th%lambda = zero
    ! We need to guarantee that the generated random numbers are consistent on
    ! all processes, so only doing this part on ionode, the gsumming the scaling
    ! factor - zamaan
    if (inode==ionode) then
      if (flag_MDdebug .and. iprint_MD > 1) &
        write(io_lun,'(2x,a)') "get_svr_thermo_sf"

      call myrng%init_rng
      call myrng%init_normal(one, zero)
      ndof = th%ndof + th%cell_ndof

      exp_dt_tau = exp(-dt/th%tau_T)
      ke = th%ke_ions
      if (present(baro)) then
        call baro%get_barostat_energy
        th%e_barostat = baro%e_barostat ! save for computing e_thermostat
        ke = ke + th%e_barostat
      end if
      temp_fac = th%ke_target/th%ndof/ke

      ! sum can be replaced by a single number drawn from a gamma distribution
      sum_ri_sq = zero
      r1 = myrng%rng_normal()
      sum_ri_sq = sum_ri_sq + r1*r1
      do i=1,ndof-1
        rn = myrng%rng_normal()
        sum_ri_sq = sum_ri_sq + rn*rn
      end do

      alpha_sq = exp_dt_tau + temp_fac * (one - exp_dt_tau) * sum_ri_sq + &
                 two*sqrt(exp_dt_tau) * sqrt(temp_fac*(one - exp_dt_tau)) * r1
      th%lambda = sqrt(alpha_sq)
    end if
    call gsum(th%lambda)
    call my_barrier

  end subroutine get_svr_thermo_sf
  !!***

  !!****m* md_control/v_rescale *
  !!  NAME
  !!   v_rescale
  !!  PURPOSE
  !!   rescale veolocity for simple velocity rescaling thermostats.
  !!   Note that the scaling factor is computed before the first vVerlet
  !!   velocity update, and the velocities scaled after the second vVerlet
  !!   velocity update.
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 14:26
  !!  SOURCE
  !!  
  subroutine v_rescale(th, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v

    v = v*th%lambda

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) then
      write(io_lun,'(2x,a)') "v_rescale"
      write(io_lun,'(4x,"lambda: ",e16.8)') th%lambda
    end if

  end subroutine v_rescale
  !!***

  !!****m* md_control/update_G_nhc *
  !!  NAME
  !!   update_G_nhc
  !!  PURPOSE
  !!   updates the "force" on thermostat k 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 10:41
  !!  SOURCE
  !!  
  subroutine update_G_nhc(th, k, e_barostat)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k      ! NH thermostat index
    real(double), intent(in)              :: e_barostat ! box ke for pressure coupling

    if (th%cell_nhc) then
      if (k == 1) then
        th%G_nhc(k) = two*th%ke_ions - th%ndof*th%T_ext*fac_Kelvin2Hartree
        th%G_nhc_cell(k) = two*e_barostat - th%T_ext*fac_Kelvin2Hartree
      else
        th%G_nhc(k) = th%m_nhc(k-1)*th%v_eta(k-1)**2 - &
                      th%T_ext*fac_Kelvin2Hartree
        th%G_nhc_cell(k) = th%m_nhc_cell(k-1)*th%v_eta_cell(k-1)**2 - &
                      th%T_ext*fac_Kelvin2Hartree
      end if
      th%G_nhc(k) = th%G_nhc(k)/th%m_nhc(k)
      th%G_nhc_cell(k) = th%G_nhc_cell(k)/th%m_nhc_cell(k)
    else
      if (k == 1) then
        th%G_nhc(k) = two*th%ke_ions - th%ndof*th%T_ext*fac_Kelvin2Hartree + &
                      two*e_barostat
      else
        th%G_nhc(k) = th%m_nhc(k-1)*th%v_eta(k-1)**2 - &
                      th%T_ext*fac_Kelvin2Hartree
      end if
      th%G_nhc(k) = th%G_nhc(k)/th%m_nhc(k)
    end if

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
      write(io_lun,'(2x,a,i2)') "update_G_nhc, k = ", k
      write(io_lun,th%nhc_fmt) "G_nhc:      ", th%G_nhc(k)
      if (th%cell_nhc) write(io_lun,th%nhc_fmt) &
        "G_nhc_cell: ", th%G_nhc_cell(k)
    end if

  end subroutine update_G_nhc
  !!***

  !!****m* md_control/propagate_eta *
  !!  NAME
  !!   propagate_eta
  !!  PURPOSE
  !!   propagate the "position" of thermostat k: linear shift
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 10:55
  !!  SOURCE
  !!  
  subroutine propagate_eta(th, k, dt, dtfac)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k      ! NH thermostat index
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter expansion factor

    th%eta(k) = th%eta(k) + dtfac*dt*th%v_eta(k)
    if (th%cell_nhc) th%eta_cell(k) = th%eta_cell(k) + &
                     dtfac*dt*th%v_eta_cell(k)

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then
      write(io_lun,'(2x,a,i2)') "propagate_eta, k = ", k
      write(io_lun,th%nhc_fmt) "eta:        ", th%eta(k)
      if (th%cell_nhc) write(io_lun,th%nhc_fmt) "eta_cell:   ", th%eta_cell(k)
    end if

  end subroutine propagate_eta
  !!***

  !!****m* md_control/propagate_v_eta_lin *
  !!  NAME
  !!   propagate_v_eta_lin
  !!  PURPOSE
  !!   propagate the "velocity" of thermostat k: linear shift 
  !!  AUTHOR
  !!   Zamaan Raza  
  !!  CREATION DATE
  !!   2017/10/24 11:43
  !!  SOURCE
  !!  
  subroutine propagate_v_eta_lin(th, k, dt, dtfac)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k      ! NH thermostat index
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter expansion factor

    th%v_eta(k) = th%v_eta(k) + dtfac*dt*th%G_nhc(k)
    if (th%cell_nhc) th%v_eta_cell(k) = th%v_eta_cell(k) + &
                                        dtfac*dt*th%G_nhc_cell(k)

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then
      write(io_lun,'(2x,a,i2)') "propagate_v_eta_lin, k = ", k
      write(io_lun,th%nhc_fmt) "v_eta:      ", th%v_eta(k)
      if (th%cell_nhc) write(io_lun,th%nhc_fmt) &
        "v_eta_cell: ", th%v_eta_cell(k)
    end if

  end subroutine propagate_v_eta_lin
  !!***

  !!****m* md_control/propagate_v_eta_exp *
  !!  NAME
  !!   propagate_v_eta_exp
  !!  PURPOSE
  !!   propagate the "velocity" of thermostat k: exponential factor 
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 11:45
  !!  SOURCE
  !!  
  subroutine propagate_v_eta_exp(th, k, dt, dtfac)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k      ! NH thermostat index
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter expansion factor

    th%v_eta(k) = th%v_eta(k)*exp(-dtfac*dt*th%v_eta(k+1))
    if (th%cell_nhc) th%v_eta_cell(k) = &
                          th%v_eta_cell(k)*exp(-dtfac*dt*th%v_eta_cell(k+1))

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then
      write(io_lun,'(2x,a,i2)') "propagate_v_eta_exp, k = ", k
      write(io_lun,th%nhc_fmt) "v_eta:      ", th%v_eta(k)
      if (th%cell_nhc) &
        write(io_lun,th%nhc_fmt) "v_eta_cell: ", th%v_eta_cell(k)
    end if

  end subroutine propagate_v_eta_exp
  !!***

  !!****m* md_control/apply_nhc_drag *
  !!  NAME
  !!   apply_nhc_drag
  !!  PURPOSE
  !!   Scale the NHC velocities by the ad-hoc drag factor
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  CREATION DATE
  !!   2018/10/31 10:00
  !!  SOURCE
  !!  
  subroutine apply_nhc_drag(th, baro, k)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    class(type_barostat), intent(inout)   :: baro
    integer, intent(in)                   :: k

    th%v_eta(k) = th%v_eta(k)*th%t_drag
    if (th%cell_nhc) th%v_eta_cell(k) = th%v_eta_cell(k)*baro%p_drag

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then
      write(io_lun,'(2x,a)') "apply_nhc_drag"
    end if

  end subroutine apply_nhc_drag
  !!***

  !!****m* md_control/integrate_nhc *
  !!  NAME
  !!   integrate_nhc
  !!  PURPOSE
  !!   propagate the particle and cell Nose-Hoover chains.
  !!   G. Martyna et al. Mol. Phys. 5, 1117 (1996)
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 11:44
  !!  MODIFICATION HISTORY
  !!   2018/10/31 zamaan
  !!    Now does both particle and box NHC integration, since they are
  !!    orthogonal. Renamed from propagte_nvt_nhc to integrate_nhc.
  !!  SOURCE
  !!  
  subroutine integrate_nhc(th, baro, v, ke)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: ke
    real(double), dimension(:,:), intent(inout) :: v  ! ion velocities

    ! local variables
    integer         :: i_mts, i_ys, i_nhc, i
    real(double)    :: v_sfac   ! ionic velocity scaling factor
    real(double)    :: box_sfac ! box velocity scaling factor
    real(double)    :: fac

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "integrate_nhc"

    th%ke_ions = ke
    v_sfac = one
    box_sfac = one
    if (th%cell_nhc) call baro%get_barostat_energy
    call th%update_G_nhc(1, baro%e_barostat) ! update force on thermostat 1
    do i_mts=1,th%n_mts_nhc ! MTS loop
      do i_ys=1,th%n_ys     ! Yoshida-Suzuki loop
        if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
          write(io_lun,*) "i_ys  = ", i_ys
          write(io_lun,*) "dt_ys = ", th%dt_ys(i_ys)
        end if
        ! Reverse part of Trotter expansion: update thermostat force/velocity
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
        call th%apply_nhc_drag(baro, th%n_nhc)
        do i_nhc=th%n_nhc-1,1,-1 ! loop over NH thermostats in reverse order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%apply_nhc_drag(baro, i_nhc)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do

        ! adjust the velocity scaling factor, scale the kinetic energy
        fac = exp(-half*th%dt_ys(i_ys)*th%v_eta(1))
        v_sfac = v_sfac*fac
        th%ke_ions = th%ke_ions*v_sfac**2

        if (th%cell_nhc) then
          ! Update the box velocities
          fac = exp(-half*th%dt_ys(i_ys)*th%v_eta_cell(1))
          select case(md_cell_constraint)
          case('volume')
            baro%v_eps = baro%v_eps*fac
          case('xyz')
            do i=1,3
              baro%v_h(i,i) = baro%v_h(i,i)*fac
            end do
          end select
          box_sfac = box_sfac*fac
          call baro%get_barostat_energy
        end if

        call th%update_G_nhc(1, baro%e_barostat) ! update force on thermostat 1
        ! update the thermostat "positions" eta
        do i_nhc=1,th%n_nhc
          call th%propagate_eta(i_nhc, th%dt_ys(i_ys), half)
        end do

        ! Forward part of Trotter expansion: update thermostat force/velocity
        do i_nhc=1,th%n_nhc-1 ! loop over NH thermostats in forward order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%update_G_nhc(i_nhc+1, baro%e_barostat)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
      end do  ! Yoshida-Suzuki loop
    end do    ! MTS loop

    ! scale the ionic velocities
    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) &
      write(io_lun,'(2x,a,e16.8)') 'v_sfac = ', v_sfac
    v = v_sfac*v
    th%lambda = v_sfac

    if (th%cell_nhc) then
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) &
        write(io_lun,'(2x,a,e16.8)') 'box_sfac = ', box_sfac
    end if

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then 
      write(io_lun,th%nhc_fmt) "eta:        ", th%eta
      write(io_lun,th%nhc_fmt) "v_eta:      ", th%v_eta
      write(io_lun,th%nhc_fmt) "G_nhc:      ", th%G_nhc
      if (th%cell_nhc) then
        write(io_lun,th%nhc_fmt) "eta_cell:   ", th%eta_cell
        write(io_lun,th%nhc_fmt) "v_eta_cell: ", th%v_eta_cell
        write(io_lun,th%nhc_fmt) "G_nhc_cell: ", th%G_nhc_cell
      end if
    end if

  end subroutine integrate_nhc
  !!***

  !!****m* md_control/get_thermostat_energy *
  !!  NAME
  !!   get_thermostat_energy
  !!  PURPOSE
  !!   compute the NHC contribution to the conserved quantity
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 13:19
  !!  MODIFICATION HISTORY
  !!   2018/05/29 zamaan
  !!    Corrected sign of potential energy contribution, added cell_ndof
  !!  SOURCE
  !!  
  subroutine get_thermostat_energy(th)

    use input_module,     only: io_assign, io_close

    ! passed variables
    class(type_thermostat), intent(inout) :: th

    ! local variables
    integer                               :: k, lun

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "get_thermostat_energy"

    select case(th%thermo_type)
    case('nhc')
      th%ke_nhc_ion = zero
      th%pe_nhc_ion = zero
      th%ke_nhc_cell = zero
      th%pe_nhc_cell = zero
      th%e_thermostat = zero
      if (th%cell_nhc) then
        th%ke_nhc_ion = half*th%m_nhc(1)*th%v_eta(1)**2
        th%pe_nhc_ion = real(th%ndof, double)*th%T_ext*fac_Kelvin2Hartree*th%eta(1)
        th%ke_nhc_cell = half*th%m_nhc_cell(1)*th%v_eta_cell(1)**2
        th%pe_nhc_cell = &
          real(th%cell_ndof, double)*th%T_ext*fac_Kelvin2Hartree*th%eta_cell(1)
      else
        th%ke_nhc_ion = half*th%m_nhc(1)*th%v_eta(1)**2
        th%pe_nhc_ion = real(th%ndof, double)*th%T_ext*fac_Kelvin2Hartree*th%eta(1)
      end if

      do k=2,th%n_nhc
        th%ke_nhc_ion = th%ke_nhc_ion + half*th%m_nhc(k)*th%v_eta(k)**2
        th%pe_nhc_ion = th%pe_nhc_ion + th%T_ext*fac_Kelvin2Hartree*th%eta(k)
        if (th%cell_nhc) then
          th%ke_nhc_cell = th%ke_nhc_cell + &
                           half*th%m_nhc_cell(k)*th%v_eta_cell(k)**2
          th%pe_nhc_cell = th%pe_nhc_cell + &
                           th%T_ext*fac_Kelvin2Hartree*th%eta_cell(k)
        end if
      end do
      th%e_thermostat = th%ke_nhc_ion + th%pe_nhc_ion + &
                        th%ke_nhc_cell + th%pe_nhc_cell

      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 2) then
        write(io_lun,'(4x,"ke_nhc_ion: ",e16.8)') th%ke_nhc_ion
        write(io_lun,'(4x,"pe_nhc_ion: ",e16.8)') th%pe_nhc_ion
        if (th%cell_nhc) then
          write(io_lun,'(4x,"ke_nhc_cell:",e16.8)') th%ke_nhc_cell
          write(io_lun,'(4x,"pe_nhc_cell:",e16.8)') th%pe_nhc_cell
        end if
        write(io_lun,'(4x,"e_nhc:      ",e16.8)') th%e_thermostat
        call io_assign(lun)
        if (th%append) then
          open(unit=lun,file='nhc.dat',position='append')
        else 
          open(unit=lun,file='nhc.dat',status='replace')
          write(lun,'(5a16)') "KE ion", "PE ion", "KE cell", "PE cell", "total"
          th%append = .true.
        end if
        if (th%cell_nhc) then
          write(lun,'(5e16.8)') th%ke_nhc_ion, th%pe_nhc_ion, &
            th%ke_nhc_cell, th%pe_nhc_cell, th%e_thermostat
        else
          write(lun,'(5e16.8)') th%ke_nhc_ion, th%pe_nhc_ion, &
            zero, zero, th%e_thermostat
        end if
        call io_close(lun)
      end if
    case('svr')
      th%e_thermostat = th%e_thermostat - &
                        th%ke_ions * (th%lambda**2 - one) * th%dt
      if (.not. leqi(th%baro_type, 'none')) then
        th%e_thermostat = th%e_thermostat - &
                          th%e_barostat * (th%lambda**2 - one) * th%dt
      end if
      if (inode==ionode .and. iprint_MD > 2) &
        write(io_lun,'(4x,"ke_svr:     ",e16.8)') th%e_thermostat
    case default
      th%e_thermostat = zero
    end select

  end subroutine get_thermostat_energy
  !!***

  !!****m* md_control/dump_thermo_state *
  !!  NAME
  !!   dump_thermo_state
  !!  PURPOSE
  !!   dump the state of the thermostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 13:36
  !!  SOURCE
  !!  
  subroutine dump_thermo_state(th, step, filename)
    use input_module,     only: io_assign, io_close

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: step
    character(len=*), intent(in)          :: filename

    ! local variables
    integer                               :: lun

    if (inode==ionode) then
      call io_assign(lun)
      if (th%append) then
        open(unit=lun,file=filename,position='append')
      else 
        open(unit=lun,file=filename,status='replace')
        th%append = .true.
      end if
      write(lun,'("step          ",i12)') step
      write(lun,'("T_int         ",f12.4)') th%T_int
      write(lun,'("ke_ions       ",e12.4)') th%ke_ions
      write(lun,'("lambda        ",f12.4)') th%lambda
      if (leqi(th%thermo_type, 'nhc')) then
        write(lun,th%nhc_fmt) "eta:          ", th%eta
        write(lun,th%nhc_fmt) "v_eta:        ", th%v_eta
        write(lun,th%nhc_fmt) "G_nhc:        ", th%G_nhc
        if (th%cell_nhc) then
          write(lun,th%nhc_fmt) "eta_cell:     ", th%eta_cell
          write(lun,th%nhc_fmt) "v_eta_cell:   ", th%v_eta_cell
          write(lun,th%nhc_fmt) "G_nhc_cell:   ", th%G_nhc_cell
          write(lun,'("ke_nhc_ion:   ",e16.4)') th%ke_nhc_ion
          write(lun,'("pe_nhc_ion:   ",e16.4)') th%pe_nhc_ion
          write(lun,'("ke_nhc_cell:  ",e16.4)') th%ke_nhc_cell
          write(lun,'("pe_nhc_cell:  ",e16.4)') th%pe_nhc_cell
        end if
        write(lun,'("e_thermostat: ",e16.4)') th%e_thermostat
      end if
      write(lun,*)
      call io_close(lun)
    end if

  end subroutine dump_thermo_state
  !!***

  !!****m* md_control/init_baro *
  !!  NAME
  !!   init_baro
  !!  PURPOSE
  !!   initialise MTTK barostat
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/17 12:44
  !!  MODIFICATION HISTORY
  !!   2019/04/09 zamaan
  !!    Removed unnecessary stress input argument
  !!  SOURCE
  !!  
  subroutine init_baro(baro, baro_type, dt, ndof, v, tau_P, ke_ions)

    use global_module,    only: rcellx, rcelly, rcellz

    ! passed variables
    class(type_barostat), intent(inout)       :: baro
    character(*), intent(in)                  :: baro_type
    real(double), intent(in)                  :: ke_ions, tau_P, dt
    integer, intent(in)                       :: ndof
    real(double), dimension(:,:), intent(in)  :: v

    ! local variables
    real(double)                              :: tauP, omega_P

    ! Globals
    baro%baro_type = baro_type
    baro%dt = dt
    baro%ndof = ndof
    baro%e_barostat = zero

    select case(md_cell_constraint)
    case('volume')
      baro%cell_ndof = 1
    case('xyz')
      baro%cell_ndof = 3
    end select

    if (md_tau_P < zero) md_tau_P = dt*ten*ten
    if (md_calc_xlmass) then
      omega_P = twopi/md_tau_P

      baro%box_mass = &
        (baro%ndof+baro%cell_ndof)*temp_ion*fac_Kelvin2Hartree/omega_P**2

      if (leqi(baro%baro_type, 'pr')) then
        if (leqi(md_cell_constraint, 'xyz')) &
          baro%box_mass = baro%box_mass/three
      end if
    else
      baro%box_mass = md_box_mass
    end if

    baro%P_ext = md_target_press/fac_HaBohr32GPa
    baro%lat_ref = zero
    baro%lat_ref(1,1) = rcellx
    baro%lat_ref(2,2) = rcelly
    baro%lat_ref(3,3) = rcellz
    baro%lat = baro%lat_ref
    call baro%get_volume
    call baro%get_pressure_and_stress
    baro%volume_ref = baro%volume
    baro%v_old = baro%volume
    baro%tau_P = tau_P
    baro%bulkmod = md_bulkmod_est
    baro%p_drag = one - (md_p_drag*dt/md_tau_P/md_n_mts/md_n_ys)

    if (flag_MDcontinue) baro%append = .true.
    select case(baro%baro_type)
    case('none')
    case('berendsen')
      flag_extended_system = .false.
    case('mttk')
      flag_extended_system = .true.
      baro%append = .false.
      baro%eps_ref = third*log(baro%volume/baro%volume_ref)
      baro%eps = baro%eps_ref
      baro%v_eps = zero
      baro%G_eps = zero
      call baro%get_barostat_energy
      baro%odnf = one + three/baro%ndof
      ! constants for polynomial expanion
      baro%c2 = one/6.0_double
      baro%c4 = baro%c2/20.0_double
      baro%c6 = baro%c4/42.0_double
      baro%c8 = baro%c6/72.0_double
    case('pr')
      flag_extended_system = .true.
      baro%append = .false.
      if (leqi(md_cell_constraint, 'volume')) then
        baro%eps = zero
        baro%v_eps = zero
        baro%G_eps = zero
      else
        baro%h = zero
        baro%v_h = zero
        baro%G_h = zero
      end if
      call baro%get_barostat_energy
      baro%odnf = one + three/baro%ndof
    case default
      call cq_abort("MD.BaroType must be 'none', 'berendsen', 'mttk' or 'pr'")
    end select

    if (inode==ionode .and. iprint_MD > 1) then
      write(io_lun,'(2x,a)') 'Welcome to init_baro'
      write(io_lun,'(4x,"Barostat type:   ",a)') baro%baro_type
      write(io_lun,'(4x,"Extended system: ",l)') flag_extended_system
      if (.not. leqi(baro_type, 'none')) then
        write(io_lun,'(4x,a,f14.2)') 'Target pressure        P_ext = ', &
                                      baro%P_ext
        write(io_lun,'(4x,a,f14.2)') 'Coupling time period   tau_P = ', &
                                      baro%tau_P
        if (flag_extended_system) then
          write(io_lun,'(4x,a,f15.8)') 'Box mass                     = ', &
                                        baro%box_mass
          write(io_lun,'(4x,a,f15.8)') 'Pressure drag factor  p_drag = ', &
                                        baro%p_drag
        end if
      end if
    end if

  end subroutine init_baro
  !!***

  !!****m* md_control/get_pressure_and_stress *
  !!  NAME
  !!   get_pressure_and_stress
  !!  PURPOSE
  !!   Compute the pressure from the stress tensor
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/08 13:27
  !!  MODIFICATION HISTORY
  !!   2018/7/23 zamaan
  !!    changed get_pressure to get_pressure_and_stress, incorporating the
  !!    subroutines get_ke_stress and update_static_stress. Now only need one
  !!    subroutine to update these quantities for less confusion in control
  !!   2018/08/11 zamaan
  !!    Added final_call optional argument so that we only print to
  !!    Conquest_out once per step
  !!   2019/04/09 zamaan
  !!    Minor change since Conquet now computes the off-diagonal stress &
  !!    tensor elements.
  !!  SOURCE
  !!  
  subroutine get_pressure_and_stress(baro, final_call)

    use force_module,     only: stress
    use move_atoms,       only: fac

    ! passed variables
    class(type_barostat), intent(inout)       :: baro
    integer, optional                         :: final_call

    ! local variables
    integer                                   :: i

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "get_pressure_and_stress"

    ! Get the static stress from the global variable
    baro%static_stress = zero
    baro%static_stress = -stress ! note the sign convention!!

    if (present(final_call)) then
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
        write(io_lun,'(4x,a,3e16.8)') "static_stress: ", &
                                      baro%static_stress(:,1)
        write(io_lun,'(4x,a,3e16.8)') "               ", &
                                      baro%static_stress(:,2)
        write(io_lun,'(4x,a,3e16.8)') "               ", &
                                      baro%static_stress(:,3)
        write(io_lun,*)
      end if
    end if

    ! Update the total stress and pressure pressure
    if (.not. (leqi(baro%baro_type, 'none'))) call baro%get_volume
    baro%total_stress = zero
    baro%P_int = zero
    baro%total_stress = (baro%ke_stress + baro%static_stress)/baro%volume
    if (present(final_call)) then
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
        write(io_lun,'(4x,a,3e16.8)') "total_stress:  ", baro%total_stress(:,1)
        write(io_lun,'(4x,a,3e16.8)') "               ", baro%total_stress(:,2)
        write(io_lun,'(4x,a,3e16.8)') "               ", baro%total_stress(:,3)
        write(io_lun,*)
      end if
    end if

    do i=1,3
      baro%P_int = baro%P_int + baro%total_stress(i,i)
    end do
    baro%P_int = baro%P_int*third
    baro%PV = baro%volume*baro%P_ext/fac_HaBohr32GPa
    if (present(final_call)) then
      if (inode==ionode .and. iprint_MD > 1) then
        write(io_lun,'(4x,"P:  ",f12.6," GPa")') baro%P_int*fac_HaBohr32GPa
        write(io_lun,'(4x,"PV: ",f12.6," Ha")') baro%PV
      end if
    end if

  end subroutine get_pressure_and_stress
  !!***

  !!****m* md_control/get_volume *
  !!  NAME
  !!   get_volume
  !!  PURPOSE
  !!   Compute the volume of the cell (orthorhombic only!)
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/09 13:15
  !!  SOURCE
  !!  
  subroutine get_volume(baro)

    ! passed variables
    class(type_barostat), intent(inout)         :: baro

    baro%volume = baro%lat(1,1)*baro%lat(2,2)*baro%lat(3,3)

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) then
      write(io_lun,'(2x,a)') "get_volume"
      write(io_lun,'(4x,a,e16.8)') "volume:     ", baro%volume
    end if

  end subroutine get_volume
  !!***

  !!****m* md_control/get_berendsen_baro_sf *
  !!  NAME
  !!   berendsen_baro_sf
  !!  PURPOSE
  !!   Get cell scaling factor for Berendsen barostat
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/12/01 14:17
  !!  SOURCE
  !!  
  subroutine get_berendsen_baro_sf(baro, dt)

    ! passed variables
    class(type_barostat), intent(inout)         :: baro
    real(double), intent(in)                    :: dt

    baro%mu = (one - (dt/baro%tau_P)*(baro%P_ext-baro%P_int)* &
                fac_HaBohr32GPa/baro%bulkmod)**third

  end subroutine get_berendsen_baro_sf
  !!***

  !!****m* md_control/propagate_berendsen *
  !!  NAME
  !!   propagate_berendsen
  !!  PURPOSE
  !!   Propagate the fractional coordinates and cell parameters for the
  !!   Berendsen barostat. Note that this is equivalent to scaling the,
  !!   fractional coordinates, and is not part of the dynamics
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/12/01 14:21
  !!  SOURCE
  !!  
  subroutine propagate_berendsen(baro, flag_movable)

    use global_module,      only: x_atom_cell, y_atom_cell, z_atom_cell, &
                                  atom_coord_diff, id_glob

    ! passed variables
    class(type_barostat), intent(inout)         :: baro
    logical, dimension(:), intent(in)           :: flag_movable

    ! local variables
    integer                                     :: i, gatom, ibeg_atom
    real(double)                                :: x_old, y_old, z_old
    logical                                     :: flagx, flagy, flagz

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) then
      write(io_lun,'(2x,a)') "propagate_berendsen"
      write(io_lun,'(4x,a,e16.8)') "mu:   :     ", baro%mu
    end if

    baro%lat = baro%lat*baro%mu

    ibeg_atom = 1
    do i=1,ni_in_cell
      gatom = id_glob(i)
      flagx = flag_movable(ibeg_atom)
      flagy = flag_movable(ibeg_atom+1)
      flagz = flag_movable(ibeg_atom+2)

      if (flagx) then
        x_old = x_atom_cell(i)
        x_atom_cell(i) = baro%mu*x_atom_cell(i)
        atom_coord_diff(1,gatom) = x_atom_cell(i) - x_old
      end if
      if (flagy) then
        y_old = y_atom_cell(i)
        y_atom_cell(i) = baro%mu*y_atom_cell(i)
        atom_coord_diff(2,gatom) = y_atom_cell(i) - y_old
      end if
      if (flagz) then
        z_old = z_atom_cell(i)
        z_atom_cell(i) = baro%mu*z_atom_cell(i)
        atom_coord_diff(3,gatom) = z_atom_cell(i) - z_old
      end if
    end do
    call baro%update_cell

  end subroutine propagate_berendsen
  !!***

  !!****m* md_control/get_barostat_energy *
  !!  NAME
  !!   get_barostat_energy
  !!  PURPOSE
  !!   get the box kinetic energy
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:09
  !!  SOURCE
  !!  
  subroutine get_barostat_energy(baro, final_call)

    ! passed variables
    class(type_barostat), intent(inout)         :: baro
    integer, optional                           :: final_call

    ! local variables
    integer                                     :: i

    baro%e_barostat = zero
    if (flag_extended_system) then
      select case(baro%baro_type)
      case('mttk')
        baro%e_barostat = half*baro%box_mass*baro%v_eps**2
      case('pr')
        if (leqi(md_cell_constraint, 'volume')) then
          baro%e_barostat = three*half*baro%box_mass*baro%v_eps**2
        else
          do i=1,3
            baro%e_barostat = baro%e_barostat + half*baro%box_mass*baro%v_h(i,i)**2
          end do
        end if
      end select
    end if

    if (present(final_call)) then
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) then
        write(io_lun,'(2x,a)') "get_barostat_energy"
        write(io_lun,'(4x,a,e16.8)') "e_barostat:     ", baro%e_barostat
      end if
    end if

  end subroutine get_barostat_energy
  !!***

  !!****m* md_control/update_G_box *
  !!  NAME
  !!   update_G_box
  !!  PURPOSE
  !!   update force on isotropic MTTK barostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 13:59
  !!  SOURCE
  !!  
  subroutine update_G_box(baro, th)

    ! passed variables
    class(type_barostat), intent(inout)         :: baro
    type(type_thermostat), intent(in)           :: th

    ! local variables
    integer                                     :: i
    real(double)                                :: tr_ke_stress

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) &
      write(io_lun,'(2x,a)') "update_G_box"

    select case(baro%baro_type)
    case('mttk')
      baro%G_eps = (two*baro%odnf*th%ke_ions + &
                    three*(baro%P_int - baro%P_ext)*baro%volume)/baro%box_mass
    case('pr')
      if (leqi(md_cell_constraint, 'volume')) then
        baro%G_eps = (two*th%ke_ions/md_ndof_ions + &
                     (baro%P_int - baro%P_ext)*baro%volume)/baro%box_mass
      else
        tr_ke_stress = zero
        do i=1,3
          tr_ke_stress = tr_ke_stress + baro%ke_stress(i,i)
        end do
        do i=1,3
          baro%G_h(i,i) = (two*tr_ke_stress/md_ndof_ions + &
                          (baro%total_stress(i,i) - baro%P_ext)*baro%volume)/ &
                           baro%box_mass
        end do
      end if
    end select

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then 
      select case(baro%baro_type)
      case('mttk')
        write(io_lun,'(2x,a,e16.8)') "G_eps: ", baro%G_eps
      case('pr')
        if (leqi(md_cell_constraint, 'volume')) then
          write(io_lun,'(2x,a,e16.8)') "G_eps: ", baro%G_eps
        else
          write(io_lun,'(2x,a,3e16.8)') "G_h:   ", baro%G_h(1,1), &
                                         baro%G_h(2,2), baro%G_h(3,3)
        end if
      end select
    end if

  end subroutine update_G_box
  !!***

  !!****m* md_control/propagate_eps_lin *
  !!  NAME
  !!   propagate_eps_lin
  !!  PURPOSE
  !!   propagate epsilon third*ln(V/V0), linear shift
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:12
  !!  SOURCE
  !!  
  subroutine propagate_eps_lin(baro, dt, dtfac)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter epxansion factor

    baro%eps = baro%eps + dtfac*dt*baro%v_eps

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) &
      write(io_lun,'(2x,a)') "propagate_eps_lin"
      write(io_lun,'(4x,a,e16.8)') "eps:        ", baro%eps

  end subroutine propagate_eps_lin
  !!***

  !!****m* md_control/propagate_eps_exp *
  !!  NAME
  !!   propagate_eps_exp
  !!  PURPOSE
  !!   propagate epsilon third*ln(V/V0), exponential factor
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:14
  !!  SOURCE
  !!  
  subroutine propagate_eps_exp(baro, dt, dtfac, v_eta_1)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter epxansion factor
    real(double), intent(in)              :: v_eta_1  ! v of first NHC thermo

    baro%eps = baro%eps*exp(-dtfac*dt*v_eta_1)

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then
      write(io_lun,'(2x,a)') "propagate_eps_exp"
      write(io_lun,'(4x,a,e16.8)') "eps:        ", baro%eps
    end if

  end subroutine propagate_eps_exp
  !!***

  !!****m* md_control/propagate_v_box_lin *
  !!  NAME
  !!   propagate_v_box_lin
  !!  PURPOSE
  !!   propagate box velocity, linear shift
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:17
  !!  SOURCE
  !!  
  subroutine propagate_v_box_lin(baro, dt, dtfac)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter epxansion factor

    ! local variables
    integer                               :: i

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) &
      write(io_lun,'(2x,a)') "propagate_v_box_lin"

    if (leqi(md_cell_constraint, 'volume')) then
      baro%v_eps = baro%v_eps + dtfac*dt*baro%G_eps
    else
      do i=1,3
        baro%v_h(i,i) = baro%v_h(i,i) + dt*dtfac*baro%G_h(i,i)
      end do
    end if

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then
      write(io_lun,*) "propagate_v_box_lin"
      if (leqi(md_cell_constraint, 'volume')) then
        write(io_lun,'(4x,a,e16.8)') "v_eps       ", baro%v_eps
      else
        write(io_lun,'(4x,a,3e16.8)') "v_h         ", &
          baro%v_h(1,1), baro%v_h(2,2), baro%v_h(3,3)
      end if
    end if

  end subroutine propagate_v_box_lin
  !!***

  !!****m* md_control/propagate_v_box_exp *
  !!  NAME
  !!   propagate_v_box_lin
  !!  PURPOSE
  !!   propagate box velocity, exponential factor
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:18
  !!  SOURCE
  !!  
  subroutine propagate_v_box_exp(baro, dt, dtfac, v_eta_1)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter epxansion factor
    real(double), intent(in)              :: v_eta_1  ! v of first NHC thermo

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) &
      write(io_lun,'(2x,a)') "propagate_v_box_lin"

    if (leqi(md_cell_constraint, 'volume')) then
      baro%v_eps = baro%v_eps*exp(-dtfac*dt*v_eta_1)
    end if

  end subroutine propagate_v_box_exp
  !!***

  !!****m* md_control/apply_box_drag *
  !!  NAME
  !!   apply_box_drag
  !!  PURPOSE
  !!   Scale the box velocities by the ad-hoc drag factor
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  CREATION DATE
  !!   2018/10/31 10:00
  !!  SOURCE
  !!  
  subroutine apply_box_drag(baro)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro

    ! local variables
    integer                               :: i

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) &
      write(io_lun,'(2x,a)') "apply_box_drag"

    if (leqi(md_cell_constraint, 'volume')) then    
      baro%v_eps = baro%v_eps*baro%p_drag
    else
      do i=1,3
        baro%v_h(i,i) = baro%v_h(i,i)*baro%p_drag
      end do
    end if

  end subroutine apply_box_drag
  !!***

  !!****m* md_control/scale_box_velocity *
  !!  NAME
  !!   scale_box_velocity
  !!  PURPOSE
  !!   Scale the box velocity during SVR NPT dynamics
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/05/17
  !!  SOURCE
  !!  
  subroutine scale_box_velocity(baro, th)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    type(type_thermostat), intent(in)     :: th

    ! local variables
    integer :: i

    do i=1,3
      baro%v_h(i,i) = baro%v_h(i,i)*th%lambda
    end do

  end subroutine scale_box_velocity
  !!***

  !!****m* md_control/update_vscale_fac *
  !!  NAME
  !!   update_vscale_fac
  !!  PURPOSE
  !!   update the NHC velocity scaling factor
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/29 17:30
  !!  SOURCE
  !!  
  subroutine update_vscale_fac(baro, th, dt, dtfac, v_eta_1, v_sfac)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    type(type_thermostat), intent(inout)  :: th
    real(double), intent(in)              :: dt      ! time step
    real(double), intent(in)              :: dtfac   ! Trotter epxansion factor
    real(double), intent(in)              :: v_eta_1 ! v of first NHC thermo
    real(double), intent(inout)           :: v_sfac  ! the scaling factor

    ! local variables
    real(double)                          :: expfac

    select case(baro%baro_type)
    case('mttk')
      expfac = exp(-dtfac*dt*(v_eta_1 + baro%odnf*baro%v_eps))
      v_sfac = v_sfac*expfac
      th%ke_ions = th%ke_ions*expfac**2
    end select

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 4) then
      write(io_lun,'(2x,a)') "update_vscale_fac"
      write(io_lun,'(4x,a,e16.8)') "v_sfac:     ", v_sfac
    end if

  end subroutine update_vscale_fac
  !!***

  !!****m* md_control/propagate_r *
  !!  NAME
  !!   propagate_r
  !!  PURPOSE
  !!   Propagate ionic positions for MTTK integrator. This is a modified
  !!   velocity Verlet step, adapted from vVerlet_r_dt in integrators_module
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 15:25
  !!  SOURCE
  !!  
  subroutine propagate_r_mttk(baro, dt, v, flag_movable)

    use global_module,      only: x_atom_cell, y_atom_cell, z_atom_cell, &
                                  atom_coord, atom_coord_diff, id_glob
    use species_module,     only: species, mass
    use move_atoms,         only: fac

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), dimension(:,:), intent(in)    :: v
    logical, dimension(:), intent(in)           :: flag_movable

    ! local variables
    real(double)                          :: exp_v_eps, sinhx_x, fac_r, &
                                             fac_v, massa, x_old, y_old, z_old
    integer                               :: i, speca, gatom, ibeg_atom
    logical                               :: flagx, flagy, flagz

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "propagate_r_mttk"

    if (leqi(md_cell_constraint, 'volume')) then
      exp_v_eps = exp(dt*half*baro%v_eps)
      sinhx_x = baro%poly_sinhx_x(half*dt*baro%v_eps)
      fac_r = exp_v_eps**2
      fac_v = exp_v_eps*sinhx_x*dt

      ibeg_atom = 1
      do i=1,ni_in_cell
        gatom = id_glob(i)
        speca = species(i)
        massa = mass(speca)*fac
        flagx = flag_movable(ibeg_atom)
        flagy = flag_movable(ibeg_atom+1)
        flagz = flag_movable(ibeg_atom+2)
        if (flagx) then
          x_old = x_atom_cell(i)
          x_atom_cell(i) = fac_r*x_atom_cell(i) + fac_v*v(1,i)
          atom_coord_diff(1,gatom) = x_atom_cell(i) - x_old
        end if
        if (flagy) then
          y_old = y_atom_cell(i)
          y_atom_cell(i) = fac_r*y_atom_cell(i) + fac_v*v(2,i)
          atom_coord_diff(2,gatom) = y_atom_cell(i) - y_old
        end if
        if (flagz) then
          z_old = z_atom_cell(i)
          z_atom_cell(i) = fac_r*z_atom_cell(i) + fac_v*v(3,i)
          atom_coord_diff(3,gatom) = z_atom_cell(i) - z_old
        end if
      end do
    end if

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
      write(io_lun,'(4x,a,e16.8)') "fac_v:      ", fac_v
      write(io_lun,'(4x,a,e16.8)') "fac_r:      ", fac_r
    end if

  end subroutine propagate_r_mttk
  !!***

  !!****m* md_control/propagate_box *
  !!  NAME
  !!   propagate_box
  !!  PURPOSE
  !!   Propagate box (lattice vectors) for MTTK integrator
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 15:29
  !!  SOURCE
  !!  
  subroutine propagate_box_mttk(baro, dt)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step

    ! local variables
    real(double)                          :: v_new, v_old, lat_sfac

      write(io_lun,'(2x,a)') "propagate_box_mttk"

    baro%v_old = baro%volume
    if (leqi(md_cell_constraint, 'volume')) then
      call baro%propagate_eps_lin(dt, one)
      v_old = baro%volume
      v_new = baro%volume_ref*exp(three*baro%eps)
      lat_sfac = (v_new/v_old)**third
      baro%mu = lat_sfac
      baro%lat = baro%lat*lat_sfac
    end if

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
      write(io_lun,'(4x,a,e16.8)') "lat_sfac:   ", lat_sfac
    end if

  end subroutine propagate_box_mttk
  !!***

  !!****m* md_control/poly_sinhx_x *
  !!  NAME
  !!   poly_sinhx_x
  !!  PURPOSE
  !!   polynomial expansion of sinh(x)/x
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 15:16
  !!  SOURCE
  !!  
  function poly_sinhx_x(baro, x) result(f)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: x

    ! local variables
    real(double)                          :: f, x2

    x2 = x**2
    f = (((baro%c8*x2 + baro%c6)*x2 + baro%c4)*x2 + baro%c2) + one

  end function poly_sinhx_x
  !!***

  !!****m* md_control/propagate_npt_mttk *
  !!  NAME
  !!   propagate_npt_mttk
  !!  PURPOSE
  !!   Propagate the NHC and barostat variables --- MTTK integrator
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 15:55
  !!  SOURCE
  !!  
  subroutine propagate_npt_mttk(baro, th, ke, v)

    ! passed variables
    class(type_barostat), intent(inout)     :: baro
    type(type_thermostat), intent(inout)    :: th
    real(double), intent(in)                :: ke
    real(double), dimension(3,3), intent(inout) :: v

    ! local variables
    integer                                 :: i_mts, i_ys, i_nhc
    real(double)                            :: v_sfac, v_eta_couple

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "propagate_npt_mttk"

    v_sfac = one
    call baro%get_barostat_energy
    call th%update_G_nhc(1, baro%e_barostat)
    if (leqi(md_cell_constraint, 'volume')) then
      call baro%update_G_box(th)
    end if

    do i_mts=1,th%n_mts_nhc ! MTS loop
      do i_ys=1,th%n_ys     ! Yoshida-Suzuki loop
        if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
          write(io_lun,*) "i_ys  = ", i_ys
          write(io_lun,*) "dt_ys = ", th%dt_ys(i_ys)
        end if
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
        th%v_eta(th%n_nhc) = th%v_eta(th%n_nhc)*th%t_drag
        th%v_eta_cell(th%n_nhc) = th%v_eta_cell(th%n_nhc)*baro%p_drag
        do i_nhc=th%n_nhc-1,1,-1 ! loop over NH thermostats in reverse order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          th%v_eta(i_nhc) = th%v_eta(i_nhc)*th%t_drag
          th%v_eta_cell(i_nhc) = th%v_eta_cell(i_nhc)*baro%p_drag
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do

        ! update box velocities
        if (th%cell_nhc) then
          v_eta_couple = th%v_eta_cell(1) ! separate NHC for the cell
        else
          v_eta_couple = th%v_eta(1)      ! same NHC as ions
        end if
        if (leqi(md_cell_constraint, 'volume')) then
          call baro%propagate_v_box_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
          call baro%propagate_v_box_lin(th%dt_ys(i_ys), quarter)
          baro%v_eps = baro%v_eps*baro%p_drag
          call baro%propagate_v_box_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
        end if

        ! update ionic velocities, scale ion kinetic energy
        if (leqi(md_cell_constraint, 'volume')) then
          call baro%update_vscale_fac(th, th%dt_ys(i_ys), half, v_eta_couple, &
                                      v_sfac)
          call baro%update_G_box(th)
        end if
        call baro%get_barostat_energy
        ! update the thermostat "positions" eta
        do i_nhc=1,th%n_nhc
          call th%propagate_eta(i_nhc, th%dt_ys(i_ys), half)
        end do

        ! update box velocities
        if (th%cell_nhc) then
          v_eta_couple = th%v_eta_cell(1) ! separate NHC for the cell
        else
          v_eta_couple = th%v_eta(1)      ! same NHC as ions
        end if
        if (leqi(md_cell_constraint, 'volume')) then
          call baro%propagate_v_box_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
          call baro%propagate_v_box_lin(th%dt_ys(i_ys), quarter)
          call baro%propagate_v_box_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
        end if

        call baro%get_barostat_energy
        call th%update_G_nhc(1, baro%e_barostat)
        do i_nhc=1,th%n_nhc-1 ! loop over NH thermostats in forward order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%update_G_nhc(i_nhc+1, baro%e_barostat)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
      end do  ! Yoshida-Suzuki loop
    end do    ! MTS loop

    ! scale the velocities
    v = v_sfac*v
    th%lambda = v_sfac

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
      write(io_lun,'(4x,a,e16.8)') "v_sfac:     ", v_sfac
    end if

  end subroutine propagate_npt_mttk
  !!***

  !!****m* md_control/integrate_box *
  !!  NAME
  !!   integrate_box
  !!  PURPOSE
  !!   Integrate the box variables, SSM integrator
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/07/12
  !!  SOURCE
  !!  
  subroutine integrate_box(baro, th)

    ! passed variables
    class(type_barostat), intent(inout)     :: baro
    type(type_thermostat), intent(in)       :: th

    ! local variables
    integer                                 :: i
    real(double)                            :: tr_ke_stress

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "integrate_box"

    if (leqi(md_cell_constraint, 'volume')) then
      baro%G_eps = (two*th%ke_ions/md_ndof_ions + &
                   (baro%P_int - baro%P_ext)*baro%volume)/baro%box_mass
      baro%v_eps = baro%v_eps + baro%G_eps*baro%dt*half
      baro%v_eps = baro%v_eps*baro%p_drag
    else
      tr_ke_stress = zero
      do i=1,3
        tr_ke_stress = tr_ke_stress + baro%ke_stress(i,i)
      end do
      do i=1,3
        baro%G_h(i,i) = (two*tr_ke_stress/md_ndof_ions + &
                        (baro%total_stress(i,i) - baro%P_ext)*baro%volume) / &
                         baro%box_mass
        baro%v_h(i,i) = baro%v_h(i,i) + baro%G_h(i,i)*baro%dt*half
        baro%v_h(i,i) = baro%v_h(i,i)*baro%p_drag
      end do
    end if

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then 
      if (leqi(md_cell_constraint, 'volume')) then
        write(io_lun,'(2x,a,e16.8)') "v_eps: ", baro%v_eps
        write(io_lun,'(2x,a,e16.8)') "G_eps: ", baro%G_eps
      else
        write(io_lun,'(2x,a,3e16.8)') "v_h:   ", baro%v_h(1,1), &
                                       baro%v_h(2,2), baro%v_h(3,3)
        write(io_lun,'(2x,a,3e16.8)') "G_h:   ", baro%G_h(1,1), &
                                       baro%G_h(2,2), baro%G_h(3,3)
      end if
    end if

  end subroutine integrate_box

  !!****m* md_control/couple_box_particle_velocity *
  !!  NAME
  !!   couple_box_particle_velocity
  !!  PURPOSE
  !!   Couple the box and particle velcities, SSM integrator
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/07/12
  !!  SOURCE
  !!  
  subroutine couple_box_particle_velocity(baro, th, v)

    ! passed variables
    class(type_barostat), intent(inout)     :: baro
    type(type_thermostat), intent(inout)    :: th
    real(double), dimension(:,:), intent(inout) :: v

    ! local variables
    real(double)                            :: expfac, const
    integer                                 :: i

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "couple_box_particle_velocity"

    const = zero
    if (leqi(md_cell_constraint, 'volume')) then
      const = three*baro%v_eps/md_ndof_ions
!      expfac = exp(-quarter*baro%dt*baro%odnf*baro%v_eps)
      expfac = exp(-quarter*baro%dt*(baro%v_eps + const))
      v  = v*expfac**2
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) &
        write(io_lun,'(4x,"expfac: ",f16.8)') expfac
    else
      do i=1,3
        const = const + baro%v_h(i,i)
      end do
      const = const/md_ndof_ions
      do i=1,3
        expfac = exp(-quarter*baro%dt*baro%odnf*(baro%v_h(i,i)+const))
        v(i,:) = v(i,:)*expfac**2
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) &
        write(io_lun,'(4x,"expfac ",i2,": ",f16.8)') i, expfac
      end do
    end if

  end subroutine couple_box_particle_velocity

  !!****m* md_control/propagate_box_ssm *
  !!  NAME
  !!   propagate_box_ssm
  !!  PURPOSE
  !!   Propagate the box, SSM integrator
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/07/12
  !!  SOURCE
  !!  
  subroutine propagate_box_ssm(baro)

    ! passed variables
    class(type_barostat), intent(inout)     :: baro

    ! local variables
    real(double)                            :: expfac, tr_vh
    real(double), dimension(3)              :: expfac_h
    integer                                 :: i

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "propagate_box_ssm"

    baro%v_old = baro%volume

    if (leqi(md_cell_constraint, 'volume')) then
      ! Propagate box velocity and compute scaling factor
      baro%eps = baro%eps + half*baro%dt*baro%v_eps
      expfac= exp(half*baro%dt*(baro%v_eps + three*baro%v_eps/ni_in_cell))
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
        write(io_lun,'(4x,"v_eps:    ",f16.8)') baro%v_eps
        write(io_lun,'(4x,"expfac:   ",f16.8)') expfac
      end if
      ! Rescale the box
      baro%lat(1,1) = baro%lat(1,1)*expfac
      baro%lat(2,2) = baro%lat(2,2)*expfac
      baro%lat(3,3) = baro%lat(3,3)*expfac
    else
      tr_vh = zero
      do i=1,3
        tr_vh = tr_vh + baro%v_h(i,i)
      end do
      do i=1,3
        baro%h(i,i) = baro%h(i,i) + half*baro%dt*baro%v_h(i,i)
        expfac_h(i) = exp(half*baro%dt*(baro%v_h(i,i) + tr_vh/md_ndof_ions))
        baro%lat(i,i) = baro%lat(i,i)*expfac_h(i)
      end do
      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
        write(io_lun,'(4x,"v_h:      ",3f16.8)') baro%v_h(1,1), &
                                                 baro%v_h(2,2), baro%v_h(3,3)
        write(io_lun,'(4x,"expfac_h: ",3f16.8)') expfac_h
      end if
    end if

  end subroutine propagate_box_ssm

  !!****m* md_control/dump_baro_state *
  !!  NAME
  !!   dump_baro_state
  !!  PURPOSE
  !!   dump the state of the barostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/09 11:49
  !!  SOURCE
  !!  
  subroutine dump_baro_state(baro, step, filename)
    use input_module,     only: io_assign, io_close

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    integer, intent(in)                   :: step
    character(len=*), intent(in)          :: filename

    ! local variables
    integer                               :: lun

    if (inode==ionode) then
      call io_assign(lun)
      if (baro%append) then
        open(unit=lun,file=filename,position='append')
      else 
        open(unit=lun,file=filename,status='replace')
        baro%append = .true.
      end if
      write(lun,'("step    ",i12)') step
      write(lun,'("P_int   ",f12.4)') baro%P_int*fac_HaBohr32GPa
      write(lun,'("volume  ",f12.4)') baro%volume
      write(lun,'("static_stress: ",3e16.8)') &
            baro%static_stress(1,1)*fac_HaBohr32GPa, &
            baro%static_stress(2,2)*fac_HaBohr32GPa, &
            baro%static_stress(3,3)*fac_HaBohr32GPa
      write(lun,'("ke_stress    : ",3e16.8)') &
            baro%ke_stress(1,1)*fac_HaBohr32GPa, &
            baro%ke_stress(2,2)*fac_HaBohr32GPa, &
            baro%ke_stress(3,3)*fac_HaBohr32GPa
      write(lun,'("cell         : ",3f16.8)') &
            baro%lat(1,1), baro%lat(2,2), baro%lat(3,3)
      if (flag_extended_system) then
        if (leqi(md_cell_constraint, 'volume')) then
          write(lun,'("eps     ",e14.6)') baro%eps
          write(lun,'("v_eps   ",e14.6)') baro%v_eps
          write(lun,'("G_eps   ",e14.6)') baro%G_eps
          write(lun,'("e_barostat  ",e14.6)') baro%e_barostat
          write(lun,'("mu      ",e14.6)') baro%mu
        else
          write(lun,'("h       ",3e14.6)') baro%h(1,1), baro%h(2,2), &
                                           baro%h(3,3)
          write(lun,'("v_h     ",3e14.6)') baro%v_h(1,1), baro%v_h(2,2), &
                                           baro%v_h(3,3)
          write(lun,'("G_h     ",3e14.6)') baro%G_h(1,1), baro%G_h(2,2), &
                                           baro%G_h(3,3)
          write(lun,'("e_barostat  ",3e14.6)') baro%e_barostat
        end if
      else
        write(lun,'("mu      ",e14.6)') baro%mu
      end if
      write(lun,*)
      call io_close(lun)
    end if

  end subroutine dump_baro_state
  !!***


  !!****m* md_control/update_cell *
  !!  NAME
  !!   update_cell
  !!  PURPOSE
  !!   Update the integration grid and density when the cell is rescaled during
  !!   a NPT MD run. Adapted from update_cell_dims in move_atoms_module
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/24 10:43
  !!  SOURCE
  !!  
  subroutine update_cell(baro)

    use units
    use global_module,      only: rcellx, rcelly, rcellz, &
                                  flag_diagonalisation
    use GenComms,           only: inode, ionode
    use density_module,     only: density
    use maxima_module,      only: maxngrid
    use dimens,             only: r_super_x, r_super_y, r_super_z, &
                                  r_super_x_squared, r_super_y_squared, &
                                  r_super_z_squared, volume, &
                                  grid_point_volume, &
                                  one_over_grid_point_volume, n_grid_x, &
                                  n_grid_y, n_grid_z
    use fft_module,         only: recip_vector, hartree_factor, i0
    use DiagModule,         only: kk, nkp

    implicit none

    ! passed variables
    class(type_barostat), intent(inout)   :: baro

    ! local variables
    real(double) :: orcellx, orcelly, orcellz, xvec, yvec, zvec, r2, scale
    integer :: i, j

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "update_cell"

    orcellx = rcellx
    orcelly = rcelly
    orcellz = rcellz

    rcellx = baro%lat(1,1)
    rcelly = baro%lat(2,2)
    rcellz = baro%lat(3,3)

    lattice_vec(1,1) = rcellx
    lattice_vec(2,2) = rcelly
    lattice_vec(3,3) = rcellz

    r_super_x = rcellx
    r_super_y = rcelly
    r_super_z = rcellz

    ! scale the integration grid and volume to the new cell.
    ! Constant number of grid points
    r_super_x_squared = r_super_x * r_super_x
    r_super_y_squared = r_super_y * r_super_y
    r_super_z_squared = r_super_z * r_super_z
    call baro%get_volume
    volume = baro%volume
    ! volume = r_super_x * r_super_y * r_super_z
    ! baro%volume = volume
    grid_point_volume = volume/(n_grid_x*n_grid_y*n_grid_z)
    one_over_grid_point_volume = one / grid_point_volume
    scale = (orcellx*orcelly*orcellz)/volume
    density = density * scale
    if(flag_diagonalisation) then
       do i = 1, nkp
          kk(1,i) = kk(1,i) * orcellx / rcellx
          kk(2,i) = kk(2,i) * orcelly / rcelly
          kk(3,i) = kk(3,i) * orcellz / rcellz
       end do
    end if
    do j = 1, maxngrid
       recip_vector(j,1) = recip_vector(j,1) * orcellx / rcellx
       recip_vector(j,2) = recip_vector(j,2) * orcelly / rcelly
       recip_vector(j,3) = recip_vector(j,3) * orcellz / rcellz
       xvec = recip_vector(j,1)/(two*pi)
       yvec = recip_vector(j,2)/(two*pi)
       zvec = recip_vector(j,3)/(two*pi)
       r2 = xvec*xvec + yvec*yvec + zvec*zvec
       if(j/=i0) hartree_factor(j) = one/r2 ! i0 notates gamma point
    end do
    if (inode == ionode .and. iprint_MD > 1) then
      write(io_lun,'(a,3f12.6)') "cell scaling factors: ", rcellx/orcellx, &
                                 rcelly/orcelly, rcellz/orcellz
      write(io_lun,'(a,3f12.6)') "new cell dimensions:  ", rcellx, rcelly, &
                                 rcellz
    end if

  end subroutine update_cell
  !!***

  !!****m* md_control/write_md_checkpoint *
  !!  NAME
  !!   write_md_checkpoint
  !!  PURPOSE
  !!   Write unified checkpoint file for MD restart (contains ionic 
  !!   velocities and extended system variables)
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/08/12 10:19
  !!  SOURCE
  !!  
  subroutine write_md_checkpoint(th, baro)

    use GenComms,         only: cq_abort
    use input_module,     only: io_assign, io_close
    use global_module,    only: id_glob, id_glob_inv
    use io_module,        only: append_coords, write_atomic_positions, &
                                pdb_template

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    class(type_barostat), intent(inout)   :: baro

    ! local variables
    integer                               :: lun, id_global, ni
    logical                               :: append_coords_bak

    if (inode==ionode) then
      if (iprint_MD > 1) &
        write(io_lun,'(2x,"Writing MD checkpoint: ",a)') md_check_file

      ! Write the ionic positions
      append_coords_bak = append_coords
      append_coords = .false.
      call write_atomic_positions(md_position_file, trim(pdb_template))
      append_coords = append_coords_bak

      call io_assign(lun)
      open(unit=lun,file=md_check_file,status='replace')

      ! Write the ionic velocities (taken from write_velocity in io_module)
      do id_global=1,ni_in_cell
        ni = id_glob_inv(id_global)
        if (id_glob(ni) /= id_global) &
          call cq_abort('Error in global labelling', id_global, id_glob(ni))
          write(lun,'(2i8,3e20.12)') id_global, ni, ion_velocity(1:3,ni)
      end do

      ! Write the extended system variables
      if (flag_extended_system) then
        if (leqi(th%thermo_type, 'nhc')) then
          write(lun,th%nhc_fmt2) th%eta
          write(lun,th%nhc_fmt2) th%v_eta
          write(lun,th%nhc_fmt2) th%G_nhc
          if (th%cell_nhc) then
            write(lun,th%nhc_fmt2) th%eta_cell
            write(lun,th%nhc_fmt2) th%v_eta_cell
            write(lun,th%nhc_fmt2) th%G_nhc_cell
          end if
        end if

        if (leqi(baro%baro_type, 'mttk')) then
          write(lun,'(3e20.12)') baro%lat_ref(1,:)
          write(lun,'(3e20.12)') baro%lat_ref(2,:)
          write(lun,'(3e20.12)') baro%lat_ref(3,:)
          write(lun,'(e20.12)') baro%volume_ref
          write(lun,'(e20.12)') baro%eps_ref
        end if

        if (leqi(md_cell_constraint, 'volume')) then
          write(lun,'(e20.12)') baro%eps
          write(lun,'(e20.12)') baro%v_eps
          write(lun,'(e20.12)') baro%G_eps
        else
          write(lun,'(3e20.12)') baro%h(1,:)
          write(lun,'(3e20.12)') baro%h(2,:)
          write(lun,'(3e20.12)') baro%h(3,:)
          write(lun,'(3e20.12)') baro%v_h(1,:)
          write(lun,'(3e20.12)') baro%v_h(2,:)
          write(lun,'(3e20.12)') baro%v_h(3,:)
          write(lun,'(3e20.12)') baro%G_h(1,:)
          write(lun,'(3e20.12)') baro%G_h(2,:)
          write(lun,'(3e20.12)') baro%G_h(3,:)
        end if
      end if
      call io_close(lun)
    end if

  end subroutine write_md_checkpoint
  !!***

  !!****m* md_control/read_md_checkpoint *
  !!  NAME
  !!   read_md_checkpoint
  !!  PURPOSE
  !!   Write checkpoint files for MD restart
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/08/12 10:19
  !!  SOURCE
  !!  
  subroutine read_md_checkpoint(th, baro)

    use GenComms,         only: cq_abort, gcopy
    use input_module,     only: io_assign, io_close
    use global_module,    only: id_glob, id_glob_inv, flag_read_velocity, &
                                species_glob

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    class(type_barostat), intent(inout)   :: baro

    ! local variables
    integer                               :: lun, id_global, ni, ni2, id_tmp
    real(double), dimension(3,ni_in_cell) :: v

    if (inode==ionode) then
      if (iprint_MD > 1) &
        write(io_lun,'(2x,"Reading MD checkpoint: ",a)') md_check_file
      call io_assign(lun)
      open(unit=lun,file=md_check_file,status='old')

      ! Position restart read in initial_read - zamaan
      ! read the ionic velocities
      do id_global=1,ni_in_cell
         ni=id_glob_inv(id_global)
         if(id_glob(ni) /= id_global) &
           call cq_abort(' ERROR in global labelling ',id_global,id_glob( ni))
         read(lun,'(2i8,3e20.12)') id_tmp, ni2, v(1:3,ni)
         if(ni2 /= ni) then
           write(io_lun,*) &
           ' Order of atom has changed for global id (file_labelling) = ', &
             id_global, &
           ' : corresponding labelling (NOprt labelling) used to be ',ni2, &
           ' : but now it is ',ni
         end if
      end do
      if (flag_read_velocity) ion_velocity = v

      ! Read the extended system variables
      if (flag_extended_system) then
        if (leqi(th%thermo_type, 'nhc')) then
          read(lun,*) th%eta
          read(lun,*) th%v_eta
          read(lun,*) th%G_nhc
          if (th%cell_nhc) then
            read(lun,*) th%eta_cell
            read(lun,*) th%v_eta_cell
            read(lun,*) th%G_nhc_cell
          end if
        end if

        if (leqi(baro%baro_type, 'mttk')) then
          read(lun,*) baro%lat_ref(1,:)
          read(lun,*) baro%lat_ref(2,:)
          read(lun,*) baro%lat_ref(3,:)
          read(lun,*) baro%volume_ref
          read(lun,*) baro%eps_ref
        end if

        if (leqi(md_cell_constraint, 'volume')) then
          read(lun,*) baro%eps
          read(lun,*) baro%v_eps
          read(lun,*) baro%G_eps
        else
          read(lun,*) baro%h(1,:)
          read(lun,*) baro%h(2,:)
          read(lun,*) baro%h(3,:)
          read(lun,*) baro%v_h(1,:)
          read(lun,*) baro%v_h(2,:)
          read(lun,*) baro%v_h(3,:)
          read(lun,*) baro%G_h(1,:)
          read(lun,*) baro%G_h(2,:)
          read(lun,*) baro%G_h(3,:)
        end if
      end if

      if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 3) then
        write(io_lun,'(2x,"flag_extended_system ",l)') flag_extended_system
        if (flag_read_velocity) then
          write(io_lun,'(a)') "Reading velocities from md.checkpoint"
          do ni=1,ni_in_cell
            write(io_lun,'(2i8,3e20.12)') ni, species_glob(ni), &
                                          ion_velocity(:,ni)
          end do
        end if
        if (flag_extended_system) then
          write(io_lun,'(a)') "Reading extended system variables &
                                &from md.checkpoint"
          if (leqi(th%thermo_type, 'nhc')) then
            write(io_lun,th%nhc_fmt) "eta:", th%eta
            write(io_lun,th%nhc_fmt) "v_eta:", th%v_eta
            write(io_lun,th%nhc_fmt) "G_nhc:", th%G_nhc
            if (th%cell_nhc) then
              write(io_lun,th%nhc_fmt) "eta_cell:", th%eta_cell
              write(io_lun,th%nhc_fmt) "v_eta_cell:", th%v_eta_cell
              write(io_lun,th%nhc_fmt) "G_eta_cell:", th%G_nhc_cell
            end if
          end if

          if (leqi(baro%baro_type, 'mttk')) then
            write(io_lun,'(2x,12a,3e20.12)') "lat_ref(1,:):", baro%lat_ref(1,:)
            write(io_lun,'(2x,12a,3e20.12)') "lat_ref(2,:):", baro%lat_ref(2,:)
            write(io_lun,'(2x,12a,3e20.12)') "lat_ref(3,:):", baro%lat_ref(3,:)
            write(io_lun,'(2x,12a,e20.12)') "volume_ref:", baro%volume_ref
            write(io_lun,'(2x,12a,e20.12)') "eps_ref:", baro%eps_ref
          end if

          if (leqi(md_cell_constraint, 'volume')) then
            write(io_lun,'(2x,12a,e20.12)') "eps:", baro%eps
            write(io_lun,'(2x,12a,e20.12)') "v_eps:", baro%v_eps
            write(io_lun,'(2x,12a,e20.12)') "G_eps:", baro%G_eps
          else
            write(io_lun,'(2x,12a,3e20.12)') "h(1,:):", baro%h(1,:)
            write(io_lun,'(2x,12a,3e20.12)') "h(2,:):", baro%h(2,:)
            write(io_lun,'(2x,12a,3e20.12)') "h(3,:):", baro%h(3,:)
            write(io_lun,'(2x,12a,3e20.12)') "v_h(1,:):", baro%v_h(1,:)
            write(io_lun,'(2x,12a,3e20.12)') "v_h(2,:):", baro%v_h(2,:)
            write(io_lun,'(2x,12a,3e20.12)') "v_h(3,:):", baro%v_h(3,:)
            write(io_lun,'(2x,12a,3e20.12)') "G_h(1,:):", baro%G_h(1,:)
            write(io_lun,'(2x,12a,3e20.12)') "G_h(2,:):", baro%G_h(2,:)
            write(io_lun,'(2x,12a,3e20.12)') "G_h(3,:):", baro%G_h(3,:)
          end if
        end if
      end if

      call io_close(lun)
    end if

    call gcopy(ion_velocity,3,ni_in_cell)
    if (flag_extended_system) then
      if (leqi(th%thermo_type, 'nhc')) then
        call gcopy(th%eta,th%n_nhc)
        call gcopy(th%v_eta,th%n_nhc)
        call gcopy(th%G_nhc,th%n_nhc)
        if (md_cell_nhc) then
          call gcopy(th%eta_cell,th%n_nhc)
          call gcopy(th%v_eta_cell,th%n_nhc)
          call gcopy(th%G_nhc_cell,th%n_nhc)
        end if
      end if
      call gcopy(baro%lat_ref,3,3)
      call gcopy(baro%volume_ref)
      call gcopy(baro%eps_ref)
      call gcopy(baro%eps)
      call gcopy(baro%v_eps)
      call gcopy(baro%G_eps)
      call gcopy(baro%h,3,3)
      call gcopy(baro%v_h,3,3)
      call gcopy(baro%G_h,3,3)
    end if

  end subroutine read_md_checkpoint
  !!***

end module md_control
