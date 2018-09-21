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
!!    Corrected sign of potential energy contribution in get_nhc_energy. Added
!!    cell degrees of freedom (cell_ndof) to equations for clarity.
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

  implicit none

  ! Unit conversion factors
  real(double), parameter :: fac_HaBohr32GPa = 29421.02648438959
  real(double), parameter :: fac_fs2atu = 41.3413745758
  real(double), parameter :: fac_invcm2hartree = 4.5563352812122295E-6

  ! Checkpoint files
  character(20) :: thermo_check_file = "cq.thermo"
  character(20) :: baro_check_file = "cq.baro"

  ! Module variables
  character(20) :: md_thermo_type, md_baro_type
  real(double)  :: md_tau_T, md_tau_P, md_target_press, md_bulkmod_est, &
                   md_box_mass, md_ndof_ions, md_omega_t, md_omega_p, &
                   md_tau_T_equil, md_tau_P_equil
  integer       :: md_n_nhc, md_n_ys, md_n_mts, md_berendsen_equil
  logical       :: flag_write_xsf, md_cell_nhc, md_calc_xlmass
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
    real(double)        :: e_nhc        ! energy of NHC thermostats
    real(double)        :: e_nhc_ion    ! energy of ionic NHC thermostats
    real(double)        :: e_nhc_cell   ! energy of cell NHC thermostats
    character(40)       :: nhc_fmt      ! format string for printing NHC arrays
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
      procedure, public   :: update_ke_ions
      procedure, public   :: get_berendsen_thermo_sf
      procedure, public   :: berendsen_v_rescale
      procedure, public   :: get_nhc_energy
      procedure, public   :: get_temperature
      procedure, public   :: dump_thermo_state
      procedure, public   :: propagate_nvt_nhc
      procedure, public   :: write_thermo_checkpoint
      procedure, public   :: read_thermo_checkpoint

      procedure, private  :: update_G_nhc
      procedure, private  :: propagate_eta
      procedure, private  :: propagate_v_eta_lin
      procedure, private  :: propagate_v_eta_exp
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
    real(double)        :: volume, volume_ref
    real(double)        :: ke_ions      ! kinetic energy of ions
    real(double)        :: dt           ! time step
    integer             :: ndof         ! number of degrees of freedom
    logical             :: append
    real(double), dimension(3,3)  :: lat       ! lattice vectors
    real(double), dimension(3,3)  :: lat_ref   ! reference lattice vectors
    real(double), dimension(3,3)  :: ke_stress      ! kinetic contrib to stress
    real(double), dimension(3,3)  :: static_stress  ! static contrib to stress

    ! constants for polynomial expansion
    real(double)        :: c2, c4, c6, c8

    ! Weak coupling barostat variables
    real(double)        :: tau_P        ! pressure coupling time period
    real(double)        :: bulkmod      ! estimated bulk modulus

    ! Extended Lagrangian barostat variables
    real(double)        :: box_mass
    real(double)        :: ke_box

    ! Isotropic variables
    real(double)        :: odnf
    real(double)        :: eps, eps_ref ! 1/3 log(V/V_0)
    real(double)        :: v_eps        ! box velocity
    real(double)        :: G_eps        ! box force
    real(double)        :: mu           ! box scaling factor

    ! Fully flexible cell variables
    real(double), dimension(3,3)  :: v_h
    real(double), dimension(3,3)  :: G_h
    real(double), dimension(3,3)  :: c_g
    real(double), dimension(3,3)  :: I_e
    real(double), dimension(3,3)  :: I_s
    real(double), dimension(3,3)  :: ident
    real(double), dimension(3,3)  :: onfm
    real(double), dimension(3)    :: lambda

    contains

      procedure, public   :: init_baro
      procedure, public   :: update_static_stress
      procedure, public   :: get_pressure
      procedure, public   :: get_volume
      procedure, public   :: get_berendsen_baro_sf
      procedure, public   :: propagate_berendsen
      procedure, public   :: get_box_ke
      procedure, public   :: get_ke_stress
      procedure, public   :: propagate_npt_mttk
      procedure, public   :: propagate_r_mttk
      procedure, public   :: propagate_box_mttk
      procedure, public   :: dump_baro_state
      procedure, public   :: update_cell
      procedure, public   :: write_baro_checkpoint
      procedure, public   :: read_baro_checkpoint

      procedure, private  :: update_G_eps
      procedure, private  :: propagate_eps_lin
      procedure, private  :: propagate_eps_exp
      procedure, private  :: propagate_v_eps_lin
      procedure, private  :: propagate_v_eps_exp
      procedure, private  :: update_vscale_fac
      procedure, private  :: poly_sinhx_x

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

    use input_module,     only: leqi

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    character(*), intent(in)              :: thermo_type, baro_type
    real(double), intent(in)              :: dt
    integer, intent(in)                   :: ndof
    real(double), intent(in)              :: ke_ions, tau_T

    if (inode==ionode) write(io_lun,'(2x,a)') 'Welcome to init_thermo'
    th%T_ext = temp_ion
    th%ndof = ndof
    th%ke_ions = ke_ions
    th%tau_T = tau_T
    th%thermo_type = thermo_type
    th%baro_type = baro_type

    select case(baro_type)
    case('iso-mttk')
      th%cell_ndof = 1
    case('ortho-mttk')
      th%cell_ndof = 3
    case('mttk')
      th%cell_ndof = 3
    end select

    call th%get_temperature
    if (inode==ionode) then
      write(io_lun,'(4x,"Thermostat type: ",a)') th%thermo_type
      write(io_lun,'(4x,a,f10.2)') 'Instantaneous temperature T_int = ', &
                                   th%T_int
      if (.not. leqi(thermo_type, 'none')) then
        write(io_lun,'(4x,a,f10.2)') 'Target temperature        T_ext = ', &
                                     th%T_ext
        write(io_lun,'(4x,a,f10.2)') 'Coupling time period      tau_T = ', &
                                     th%tau_T
      end if
    end if

    select case(th%thermo_type)
    case('none')
    case('berendsen')
    case('nhc')
      call th%init_nhc(dt)
    case default
      call cq_abort("Unknown thermostat type")
    end select

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

    th%n_nhc = md_n_nhc
    th%n_ys = md_n_ys
    th%n_mts_nhc = md_n_mts
    th%lambda = one
    th%cell_nhc  = md_cell_nhc
    call th%get_temperature
    write(th%nhc_fmt,'("(4x,a12,",i4,"f10.4)")') th%n_nhc

    allocate(th%eta(th%n_nhc))
    allocate(th%v_eta(th%n_nhc))
    allocate(th%G_nhc(th%n_nhc))
    allocate(th%m_nhc(th%n_nhc))
    call reg_alloc_mem(area_moveatoms, 4*th%n_nhc, type_dbl)
    if (th%cell_nhc) then
      allocate(th%eta_cell(th%n_nhc))    ! independent thermostat for cell DOFs
      allocate(th%v_eta_cell(th%n_nhc))
      allocate(th%G_nhc_cell(th%n_nhc))
      allocate(th%m_nhc_cell(th%n_nhc))
      call reg_alloc_mem(area_moveatoms, 4*th%n_nhc, type_dbl)
    end if

    ! Calculate the masses for extended lagrangian variables?
    if (md_calc_xlmass) then
      ! convert time scales from fs to atomic units
      tauT = md_tau_T*fac_fs2atu 
      tauP = md_tau_P*fac_fs2atu 
      omega_thermo = twopi/tauT
      omega_baro = twopi/tauP
      th%m_nhc(1) = md_ndof_ions*th%T_ext*fac_Kelvin2Hartree/omega_thermo**2
      if (th%cell_nhc) then
        select case (md_baro_type)
        case('iso-mttk')
          ndof_baro = one
        end select
        th%m_nhc_cell(1) = (ndof_baro**2)*th%T_ext*fac_Kelvin2Hartree/omega_baro**2
      end if
      do i=2,th%n_nhc
        th%m_nhc(i) = th%T_ext*fac_Kelvin2Hartree/omega_thermo**2
        if (th%cell_nhc) th%m_nhc_cell(i) = th%T_ext*fac_Kelvin2Hartree/omega_thermo**2
      end do
    else
      th%m_nhc = md_nhc_mass
      th%m_nhc_cell = md_nhc_cell_mass
    end if

    ! Defaults for heat bath positions, velocities, masses
    if (flag_MDcontinue) then
      th%append = .true.
      call th%read_thermo_checkpoint
    else
      th%append = .false.
      th%eta = zero
      th%v_eta = sqrt(two*th%T_ext*fac_Kelvin2Hartree/th%m_nhc(1)) 
      th%G_nhc = zero
      if (th%cell_nhc) then
        th%eta_cell = zero
        th%v_eta_cell = sqrt(two*th%T_ext*fac_Kelvin2Hartree/th%m_nhc(1)) 
        th%G_nhc_cell = zero
      end if
    end if
    call th%get_nhc_energy

    ! Yoshida-Suzuki time steps
    call th%init_ys(dt, th%n_ys)

    write(fmt1,'("(4x,a16,",i4,"f12.2)")') th%n_nhc
    write(fmt2,'("(4x,a16,",i4,"f12.2)")') th%n_ys
    if (inode==ionode) then
      write(io_lun,'(2x,a)') 'Welcome to init_nhc'
      write(io_lun,'(4x,a,i10)')   'Number of NHC thermostats n_nhc = ', &
                                   th%n_nhc
      write(io_lun,'(4x,a,i10)')   'Multiple time step order  n_mts = ', &
                                   th%n_mts_nhc
      write(io_lun,'(4x,a,i10)')   'Yoshida-Suzuki order      n_ys  = ', &
                                   th%n_ys
      write(io_lun,fmt2) 'YS time steps:  ', th%dt_ys
      write(io_lun,fmt1) 'ion  NHC masses:', th%m_nhc
      if (th%cell_nhc) write(io_lun,fmt1) 'cell NHC masses:', th%m_nhc_cell
    end if

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
    th%dt_ys = dt*th%dt_ys/th%n_mts_nhc
    deallocate(psuz)

  end subroutine init_ys
  !!***

  !!****m* md_control/update_ke_ions *
  !!  NAME
  !!   update_ke_ions
  !!  PURPOSE
  !!   update ionic kinetic energy when thermo_type = 'None'
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/09 10:39
  !!  SOURCE
  !!  
  subroutine update_ke_ions(th, ke_ions)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    real(double), intent(in)              :: ke_ions

    th%ke_ions = ke_ions

  end subroutine update_ke_ions
  !!***

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

  end subroutine get_berendsen_thermo_sf
  !!***

  !!****m* md_control/berendsen_v_rescale *
  !!  NAME
  !!   berendsen_v_rescale
  !!  PURPOSE
  !!   rescale veolocity using Berendsen (weak coupling) velocity rescaling.
  !!   Note that the scaling factor is computed before the first vVerlet
  !!   velocity update, and the velocities scaled after the second vVerlet
  !!   velocity update.
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 14:26
  !!  SOURCE
  !!  
  subroutine berendsen_v_rescale(th, v)

    ! passed variables
    class(type_thermostat), intent(inout)       :: th
    real(double), dimension(:,:), intent(inout) :: v

    v = v*th%lambda

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,'(a,i2)') "thermoDebug: berendsen_v_rescale"
      write(io_lun,*) "lambda:     ", th%lambda
    end if

  end subroutine berendsen_v_rescale
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
  subroutine update_G_nhc(th, k, ke_box)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k      ! NH thermostat index
    real(double), intent(in)              :: ke_box ! box ke for pressure coupling

    if (th%cell_nhc) then
      if (k == 1) then
        th%G_nhc(k) = two*th%ke_ions - th%ndof*th%T_ext*fac_Kelvin2Hartree
        th%G_nhc_cell(k) = two*ke_box - th%T_ext*fac_Kelvin2Hartree
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
                      two*ke_box
      else
        th%G_nhc(k) = th%m_nhc(k-1)*th%v_eta(k-1)**2 - &
                      th%T_ext*fac_Kelvin2Hartree
      end if
      th%G_nhc(k) = th%G_nhc(k)/th%m_nhc(k)
    end if

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,'(a,i2)') "thermoDebug: update_G_nhc, k = ", k
      write(io_lun,*) "G_eta:      ", th%G_nhc(k)
      if (th%cell_nhc) write(io_lun,*) "G_eta_cell: ", th%G_nhc_cell(k)
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

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,'(a,i2)') "thermoDebug: propagate_eta, k = ", k
      write(io_lun,*) "eta:        ", th%eta(k)
      if (th%cell_nhc) write(io_lun,*) "eta_cell:   ", th%eta_cell(k)
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

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,'(a,i2)') "thermoDebug: propagate_v_eta_lin, k = ", k
      write(io_lun,*) "v_eta:      ", th%v_eta(k)
      if (th%cell_nhc) write(io_lun,*) "v_eta_cell: ", th%v_eta_cell(k)
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

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,'(a,i2)') "thermoDebug: propagate_v_eta_exp, k = ", k
      write(io_lun,*) "v_eta:      ", th%v_eta(k)
      if (th%cell_nhc) write(io_lun,*) "v_eta_cell: ", th%v_eta_cell(k)
    end if

  end subroutine propagate_v_eta_exp
  !!***

  !!****m* md_control/propagate_nvt_nhc *
  !!  NAME
  !!   propagate_nvt_nhc
  !!  PURPOSE
  !!   propagate the Nose-Hoover chain. 
  !!   G. Martyna et al. Mol. Phys. 5, 1117 (1996)
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 11:44
  !!  SOURCE
  !!  
  subroutine propagate_nvt_nhc(th, v, ke)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    real(double), intent(in)              :: ke
    real(double), dimension(:,:), intent(inout) :: v  ! ion velocities

    ! local variables
    integer       :: i_mts, i_ys, i_nhc
    real(double)  :: v_sfac   ! ionic velocity scaling factor
    real(double)  :: fac

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "thermoDebug: propagate_nvt_nhc"

    th%ke_ions = ke
    v_sfac = one
    call th%update_G_nhc(1, zero) ! update force on thermostat 1
    do i_mts=1,th%n_mts_nhc ! MTS loop
      do i_ys=1,th%n_ys     ! Yoshida-Suzuki loop
        if (inode==ionode .and. flag_MDdebug) then
          write(io_lun,*) "thermoDebug: i_ys  = ", i_ys
          write(io_lun,*) "thermoDebug: dt_ys = ", th%dt_ys(i_ys)
        end if
        ! Reverse part of Trotter expansion: update thermostat force/velocity
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
        do i_nhc=th%n_nhc-1,1,-1 ! loop over NH thermostats in reverse order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do

        ! adjust the velocity scaling factor, scale the kinetic energy
        fac = exp(-half*th%dt_ys(i_ys)*th%v_eta(1))
        v_sfac = v_sfac*fac
        th%ke_ions = th%ke_ions*v_sfac**2

        call th%update_G_nhc(1, zero) ! update force on thermostat 1
        ! update the thermostat "positions" eta
        do i_nhc=1,th%n_nhc
          call th%propagate_eta(i_nhc, th%dt_ys(i_ys), half)
        end do

        ! Forward part of Trotter expansion: update thermostat force/velocity
        do i_nhc=1,th%n_nhc-1 ! loop over NH thermostats in forward order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%update_G_nhc(i_nhc+1, zero)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
      end do  ! Yoshida-Suzuki loop
    end do    ! MTS loop

    ! scale the ionic velocities
    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) 'thermoDebug: v_sfac = ', v_sfac

    v = v_sfac*v
    th%lambda = v_sfac

  end subroutine propagate_nvt_nhc
  !!***

  !!****m* md_control/write_thermo_checkpoint *
  !!  NAME
  !!   write_thermo_checkpoint
  !!  PURPOSE
  !!   Write checkpoint for restarting NHC thermostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/01/19 13:12
  !!  SOURCE
  !!  
  subroutine write_thermo_checkpoint(th)
    
    use input_module,     only: io_assign, io_close

    ! passed variables
    class(type_thermostat), intent(inout) :: th

    ! local variables
    integer                               :: lun

    call io_assign(lun)
    open(unit=lun,file=thermo_check_file,status='replace')

    write(lun,*) "eta:        ", th%eta
    write(lun,*) "v_eta:      ", th%v_eta
    write(lun,*) "G_nhc:      ", th%G_nhc
    if (th%cell_nhc) then
      write(lun,*) "eta_cell:   ", th%eta_cell
      write(lun,*) "v_eta_cell: ", th%v_eta_cell
      write(lun,*) "G_nhc_cell: ", th%G_nhc_cell
    end if
    call io_close(lun)

  end subroutine write_thermo_checkpoint
  !!***

  !!****m* md_control/read_thermo_checkpoint *
  !!  NAME
  !!   read_thermo_checkpoint
  !!  PURPOSE
  !!   Read checkpoint for restarting NHC thermostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/01/19 13:14
  !!  SOURCE
  !!  
  subroutine read_thermo_checkpoint(th)
    
    use input_module,     only: io_assign, io_close

    ! passed variables
    class(type_thermostat), intent(inout) :: th

    ! local variables
    integer                               :: lun
    character(20)                         :: name

    call io_assign(lun)
    open(unit=lun,file=thermo_check_file,status='old')

    read(lun,*) name, th%eta
    read(lun,*) name, th%v_eta
    read(lun,*) name, th%G_nhc
    if (th%cell_nhc) then
      read(lun,*) name, th%eta_cell
      read(lun,*) name, th%v_eta_cell
      read(lun,*) name, th%G_nhc_cell
    end if
    call io_close(lun)

  end subroutine read_thermo_checkpoint
  !!***

  !!****m* md_control/get_nhc_energy *
  !!  NAME
  !!   get_nhc_energy
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
  subroutine get_nhc_energy(th)

    ! passed variables
    class(type_thermostat), intent(inout) :: th

    ! local variables
    integer                               :: k

    th%e_nhc_cell = zero
    if (th%cell_nhc) then
      th%e_nhc_ion = th%m_nhc(1)*th%v_eta(1)**2 - &
                     real(th%ndof, double)*th%T_ext*fac_Kelvin2Hartree*th%eta(1)
      th%e_nhc_cell = th%e_nhc_cell + &
                      th%m_nhc_cell(1)*th%v_eta_cell(1)**2 - &
                      th%cell_ndof*th%T_ext*fac_Kelvin2Hartree*th%eta_cell(1)
    else
      th%e_nhc_ion = th%m_nhc(1)*th%v_eta(1)**2 - &
                     real((th%ndof+th%cell_ndof), double) * &
                     th%T_ext*fac_Kelvin2Hartree*th%eta(1)
    end if

    do k=2,th%n_nhc
      th%e_nhc_ion = th%e_nhc_ion + th%m_nhc(k)*th%v_eta(k)**2 - &
                     th%T_ext*fac_Kelvin2Hartree*th%eta(k)
      if (th%cell_nhc) then
        th%e_nhc_cell = th%e_nhc_cell + &
                        th%m_nhc_cell(k)*th%v_eta_cell(k)**2 - &
                        th%T_ext*fac_Kelvin2Hartree*th%eta_cell(k)
      end if
    end do
    th%e_nhc = th%e_nhc_ion + th%e_nhc_cell

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,*) "thermoDebug: get_nhc_energy: e_nhc = ", th%e_nhc
      write(io_lun,*) "                        e_nhc_ion  = ", th%e_nhc_ion
      write(io_lun,*) "                        e_nhc_cell = ", th%e_nhc_cell
    end if

  end subroutine get_nhc_energy
  !!***

  !!****m* md_control/get_temperature *
  !!  NAME
  !!   get_temperature
  !!  PURPOSE
  !!   compute the instantaneous temperature
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/27 16:48
  !!  SOURCE
  !!  
  subroutine get_temperature(th)

    ! passed variables
    class(type_thermostat), intent(inout) :: th

    th%T_int = two*th%ke_ions/fac_Kelvin2Hartree/md_ndof_ions
    if (inode==ionode) write(io_lun,'(4x,"T = ", f12.6)') th%T_int

  end subroutine get_temperature
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
      write(lun,'("step        ",i12)') step
      write(lun,'("T_int       ",f12.4)') th%T_int
      write(lun,'("ke_ions     ",e12.4)') th%ke_ions
      write(lun,'("lambda      ",f12.4)') th%lambda
      if (th%thermo_type == 'nhc') then
        write(lun,th%nhc_fmt) "eta:        ", th%eta
        write(lun,th%nhc_fmt) "v_eta:      ", th%v_eta
        write(lun,th%nhc_fmt) "G_nhc:      ", th%G_nhc
        if (th%cell_nhc) then
          write(lun,th%nhc_fmt) "eta_cell:   ", th%eta_cell
          write(lun,th%nhc_fmt) "v_eta_cell: ", th%v_eta_cell
          write(lun,th%nhc_fmt) "G_nhc_cell: ", th%G_nhc_cell
          write(lun,'("e_nhc_ion:  ",e16.4)') th%e_nhc_ion
          write(lun,'("e_nhc_cell: ",e16.4)') th%e_nhc_cell
        end if
        write(lun,'("e_nhc:      ",e16.4)') th%e_nhc
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
  !!  SOURCE
  !!  
  subroutine init_baro(baro, baro_type, ndof, stress, v, tau_P, ke_ions)

    use input_module,     only: leqi
    use global_module,    only: rcellx, rcelly, rcellz

    ! passed variables
    class(type_barostat), intent(inout)       :: baro
    character(*), intent(in)                  :: baro_type
    real(double), intent(in)                  :: ke_ions, tau_P
    integer, intent(in)                       :: ndof
    real(double), dimension(3), intent(in)    :: stress
    real(double), dimension(:,:), intent(in)  :: v

    ! local variables
    real(double)                              :: tauP, omega_P

    ! Globals
    baro%baro_type = baro_type
    if (md_calc_xlmass) then
      select case(baro_type)
      case('iso-mttk')
        tauP = md_tau_P*fac_fs2atu 
        omega_P = twopi/tauP
        baro%box_mass = (md_ndof_ions+one)*temp_ion*fac_Kelvin2Hartree/omega_P**2
      end select
    else
      baro%box_mass = md_box_mass
    end if

    ! constants for polynomial expanion
    baro%c2 = one/6.0_double
    baro%c4 = baro%c2/20.0_double
    baro%c6 = baro%c4/42.0_double
    baro%c8 = baro%c6/72.0_double

    baro%P_ext = md_target_press/fac_HaBohr32GPa
    baro%ke_ions = ke_ions
    baro%ndof = ndof
    baro%lat_ref = zero
    baro%lat_ref(1,1) = rcellx
    baro%lat_ref(2,2) = rcelly
    baro%lat_ref(3,3) = rcellz
    baro%lat = baro%lat_ref
    call baro%get_volume
    call baro%update_static_stress(stress)
    call baro%get_ke_stress(v)
    call baro%get_pressure
    baro%volume_ref = baro%volume
    baro%tau_P = tau_P

    select case(baro%baro_type)
    case('none')
    case('berendsen')
      baro%bulkmod = md_bulkmod_est
    case('iso-mttk')
      if (flag_MDcontinue) then
        baro%append = .true.
        call baro%read_baro_checkpoint
        call baro%get_box_ke
      else
        baro%append = .false.
        baro%eps_ref = third*log(baro%volume/baro%volume_ref)
        baro%eps = baro%eps_ref
        baro%v_eps = zero
        call baro%get_box_ke
      end if
      baro%odnf = one + three/baro%ndof
    case('ortho-mttk')
    case('mttk')
    case default
      call cq_abort("Invalid barostat type")
    end select

    if (inode==ionode) then
      write(io_lun,'(2x,a)') 'Welcome to init_baro'
      write(io_lun,'(4x,"Barostat type: ",a)') baro%baro_type
      write(io_lun,'(4x,a,f14.2)') 'Instantaneous pressure P_int = ', &
                                    baro%P_int
      if (.not. leqi(baro_type, 'none')) then
        write(io_lun,'(4x,a,f14.2)') 'Target pressure        P_ext = ', &
                                      baro%P_ext
        write(io_lun,'(4x,a,f14.2)') 'Coupling time period   tau_P = ', &
                                      baro%tau_P
        if (leqi(baro%baro_type, 'iso-mttk')) then
          write(io_lun,'(4x,a,f14.2)') 'Box mass                     = ', &
                                        baro%box_mass
        end if
      end if
    end if

  end subroutine init_baro
  !!***

  !!****m* md_control/update_static_stress *
  !!  NAME
  !!   update_static_stress
  !!  PURPOSE
  !!   update static stress when baro_type = 'None'
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/09 10:44
  !!  SOURCE
  !!  
  subroutine update_static_stress(baro, stress)

    ! passed variables
    class(type_barostat), intent(inout)     :: baro
    real(double), dimension(3), intent(in)  :: stress

    ! local variables
    integer                                 :: i

    baro%static_stress = zero
    do i=1,3
      baro%static_stress(i,i) = stress(i)
    end do

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,*) "baroDebug: update_static_stress"
      write(io_lun,'(a,3f15.8)') "static_stress: ", baro%static_stress(:,1)
      write(io_lun,'(a,3f15.8)') "static_stress: ", baro%static_stress(:,2)
      write(io_lun,'(a,3f15.8)') "static_stress: ", baro%static_stress(:,3)
      write(io_lun,*)
    end if

  end subroutine update_static_stress
  !!***

  !!****m* md_control/get_ke_stress *
  !!  NAME
  !!   get_ke_stress
  !!  PURPOSE
  !!   Compute the kinetic contribution to the stress tensor
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/08 12:23
  !!  SOURCE
  !!  
  subroutine get_ke_stress(baro, v)

    use move_atoms,       only: fac

    ! passed variables
    class(type_barostat), intent(inout)       :: baro
    real(double), dimension(3,ni_in_cell), intent(in)  :: v

    ! local variables
    integer                                   :: iatom, j, k
    real(double)                              :: m

    baro%ke_stress = zero
    do iatom=1,ni_in_cell
      m = mass(species(iatom))*fac
      do j=1,3
        do k=1,3
          baro%ke_stress(j,k) = baro%ke_stress(j,k) + m*v(j,iatom)*v(k,iatom)
        end do
      end do
    end do 
    baro%ke_stress = baro%ke_stress/baro%volume

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,*) "baroDebug: get_ke_stress"
      write(io_lun,'(a,3f15.8)') "ke_stress:  ", baro%ke_stress(:,1)
      write(io_lun,'(a,3f15.8)') "ke_stress:  ", baro%ke_stress(:,2)
      write(io_lun,'(a,3f15.8)') "ke_stress:  ", baro%ke_stress(:,3)
      write(io_lun,*)
    end if

  end subroutine get_ke_stress
  !!***

  !!****m* md_control/get_pressure *
  !!  NAME
  !!   get_pressure
  !!  PURPOSE
  !!   Compute the pressure from the stress tensor
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/08 13:27
  !!  SOURCE
  !!  
  subroutine get_pressure(baro)

    ! passed variables
    class(type_barostat), intent(inout)         :: baro

    ! local variables
    integer                                     :: i

    baro%P_int = zero
    do i=1,3
      baro%P_int = baro%P_int + baro%ke_stress(i,i) + baro%static_stress(i,i)
    end do
    baro%P_int = -baro%P_int*third/baro%volume

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,*) "baroDebug: get_pressure"
      write(io_lun,*) "P_int:      ", baro%P_int
    end if

  end subroutine get_pressure
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
    baro%PV = baro%volume*baro%P_ext

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,*) "baroDebug: get_volume"
      write(io_lun,*) "volume:     ", baro%volume
      write(io_lun,*) "PV    :     ", baro%PV
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

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,*) "baroDebug: propagate_berendsen"
      write(io_lun,*) "mu:   :     ", baro%mu
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

  end subroutine propagate_berendsen
  !!***

  !!****m* md_control/get_box_ke *
  !!  NAME
  !!   get_box_ke
  !!  PURPOSE
  !!   get the box kinetic energy
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:09
  !!  SOURCE
  !!  
  subroutine get_box_ke(baro)

    ! passed variables
    class(type_barostat), intent(inout)         :: baro

    select case(baro%baro_type)
    case('iso-mttk')
      baro%ke_box = half*baro%box_mass*baro%v_eps**2
    case('ortho-mttk')
    case('mttk')
    end select

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: get_box_ke: ke_box = ", baro%ke_box

  end subroutine get_box_ke
  !!***

  !!****m* md_control/update_G_eps *
  !!  NAME
  !!   update_G_eps
  !!  PURPOSE
  !!   update force on isotropic MTTK barostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 13:59
  !!  SOURCE
  !!  
  subroutine update_G_eps(baro)

    ! passed variables
    class(type_barostat), intent(inout)         :: baro

    baro%G_eps = (two*baro%odnf*baro%ke_ions + three*(baro%P_int - baro%P_ext)*baro%volume)/baro%box_mass
    ! baro%G_eps = (three*baro%ke_ions/baro%ndof + three*(baro%P_int - baro%P_ext)*baro%volume)/baro%box_mass

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: update_G_eps: G_eps = ", baro%G_eps

  end subroutine update_G_eps
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

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_eps_lin: eps = ", baro%eps

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

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_eps_exp: eps = ", baro%eps

  end subroutine propagate_eps_exp
  !!***

  !!****m* md_control/propagate_v_eps_lin *
  !!  NAME
  !!   propagate_v_eps_lin
  !!  PURPOSE
  !!   propagate box velocity, linear shift
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:17
  !!  SOURCE
  !!  
  subroutine propagate_v_eps_lin(baro, dt, dtfac)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter epxansion factor

    baro%v_eps = baro%v_eps + dtfac*dt*baro%G_eps

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_v_eps_lin: v_eps = ", baro%v_eps

  end subroutine propagate_v_eps_lin
  !!***

  !!****m* md_control/propagate_v_eps_exp *
  !!  NAME
  !!   propagate_v_eps_lin
  !!  PURPOSE
  !!   propagate box velocity, exponential factor
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/11/17 14:18
  !!  SOURCE
  !!  
  subroutine propagate_v_eps_exp(baro, dt, dtfac, v_eta_1)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter epxansion factor
    real(double), intent(in)              :: v_eta_1  ! v of first NHC thermo

    baro%v_eps = baro%v_eps*exp(-dtfac*dt*v_eta_1)

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_v_eps_lin: v_eps = ", baro%v_eps

  end subroutine propagate_v_eps_exp
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
  subroutine update_vscale_fac(baro, dt, dtfac, v_eta_1, v_sfac)

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt      ! time step
    real(double), intent(in)              :: dtfac   ! Trotter epxansion factor
    real(double), intent(in)              :: v_eta_1 ! v of first NHC thermo
    real(double), intent(inout)           :: v_sfac  ! the scaling factor

    ! local variables
    real(double)                          :: expfac

    select case(baro%baro_type)
    case('iso-mttk')
      expfac = exp(-dtfac*dt*(v_eta_1 + baro%odnf*baro%v_eps))
      v_sfac = v_sfac*expfac
      baro%ke_ions = baro%ke_ions*expfac**2
    end select

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: update_vscale_fac: v_sfac = ", v_sfac

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

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_r_mttk"

    select case(baro%baro_type)
    case('iso-mttk')
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
    case('ortho-mttk')
    case('mttk')
    end select

    if (inode==ionode .and. flag_MDdebug) then
      write(io_lun,*) "baroDebug: propagate_r_mttk: fac_v = ", fac_v
      write(io_lun,*) "baroDebug: propagate_r_mttk: fac_r = ", fac_r
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

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_box_mttk"

    select case(baro%baro_type)
    case('iso-mttk')
      call baro%propagate_eps_lin(dt, one)
      v_old = baro%volume
      v_new = baro%volume_ref*exp(three*baro%eps)
      lat_sfac = (v_new/v_old)**third
      baro%mu = lat_sfac
      baro%lat = baro%lat*lat_sfac
    case('ortho-mttk')
    case('mttk')
    end select

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_box_mttk: lat_sfac = ", lat_sfac

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

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_npt_mttk"

    baro%ke_ions = ke
    v_sfac = one
    call baro%get_box_ke
    call th%update_G_nhc(1, baro%ke_box)
    select case(baro%baro_type)
    case('iso-mttk')
      call baro%update_G_eps
    end select

    do i_mts=1,th%n_mts_nhc ! MTS loop
      do i_ys=1,th%n_ys     ! Yoshida-Suzuki loop
        if (inode==ionode .and. flag_MDdebug) then
          write(io_lun,*) "baroDebug: i_ys  = ", i_ys
          write(io_lun,*) "baroDebug: dt_ys = ", th%dt_ys(i_ys)
        end if
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
        do i_nhc=th%n_nhc-1,1,-1 ! loop over NH thermostats in reverse order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do

        ! update box velocities
        if (th%cell_nhc) then
          v_eta_couple = th%v_eta_cell(1) ! separate NHC for the cell
        else
          v_eta_couple = th%v_eta(1)      ! same NHC as ions
        end if
        select case(baro%baro_type)
        case('iso-mttk')
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
          call baro%propagate_v_eps_lin(th%dt_ys(i_ys), quarter)
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
        case('ortho-mttk')
        case('mttk')
        end select

        ! update ionic velocities, scale ion kinetic energy
        call baro%update_vscale_fac(th%dt_ys(i_ys), half, v_eta_couple, &
                                    v_sfac)
        th%ke_ions = baro%ke_ions

        select case(baro%baro_type)
        case('iso-mttk')
          call baro%update_G_eps
        case('ortho-mttk')
        case('mttk')
        end select
        call baro%get_box_ke
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
        select case(baro%baro_type)
        case('iso-mttk')
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
          call baro%propagate_v_eps_lin(th%dt_ys(i_ys), quarter)
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, &
                                        v_eta_couple)
        case('ortho-mttk')
        case('mttk')
        end select

        call baro%get_box_ke
        call th%update_G_nhc(1, baro%ke_box)
        do i_nhc=1,th%n_nhc-1 ! loop over NH thermostats in forward order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%update_G_nhc(i_nhc+1, baro%ke_box)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
      end do  ! Yoshida-Suzuki loop
    end do    ! MTS loop

    ! scale the velocities
    v = v_sfac*v
    th%lambda = v_sfac

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: propagate_npt_mttk: v_sfac = ", v_sfac

  end subroutine propagate_npt_mttk
  !!***

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
      select case(baro%baro_type)
      case('iso-mttk')
        write(lun,'("eps     ",f12.6)') baro%eps
        write(lun,'("v_eps   ",f12.6)') baro%v_eps
        write(lun,'("G_eps   ",f12.6)') baro%G_eps
        write(lun,'("ke_box  ",f12.6)') baro%ke_box
        write(lun,'("mu      ",f12.6)') baro%mu
      case('ortho-mttk')
      case('mttk')
      case('berendsen')
        write(lun,'("mu      ",f12.6)') baro%mu
      end select
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
    use input_module,       only: leqi

    implicit none

    ! passed variables
    class(type_barostat), intent(inout)   :: baro

    ! local variables
    real(double) :: orcellx, orcelly, orcellz, xvec, yvec, zvec, r2, scale
    integer :: i, j

    if (inode==ionode .and. flag_MDdebug) &
      write(io_lun,*) "baroDebug: update_cell"

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
    if (inode == ionode) then
      write(io_lun,'(a,3f12.6)') "cell scaling factors: ", rcellx/orcellx, &
                                 rcelly/orcelly, rcellz/orcellz
      write(io_lun,'(a,3f12.6)') "new cell dimensions:  ", rcellx, rcelly, &
                                 rcellz
    end if

  end subroutine update_cell
  !!***

  !!****m* md_control/write_baro_checkpoint *
  !!  NAME
  !!   write_baro_checkpoint
  !!  PURPOSE
  !!   Write checkpoint for restarting MTTK barostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/01/19 13:15
  !!  SOURCE
  !!  
  subroutine write_baro_checkpoint(baro)

    use input_module,     only: io_assign, io_close

    ! passed variables
    class(type_barostat), intent(inout) :: baro

    ! local variables
    integer                               :: lun

    call io_assign(lun)
    open(unit=lun,file=baro_check_file,status='replace')

    write(lun,*) "lat_a_ref", baro%lat_ref(1,:)
    write(lun,*) "lat_b_ref", baro%lat_ref(2,:)
    write(lun,*) "lat_c_ref", baro%lat_ref(3,:)
    write(lun,*) "volume_ref", baro%volume_ref
    write(lun,*) "eps_ref ", baro%eps_ref
    write(lun,*) "eps ", baro%eps
    write(lun,*) "v_eps ", baro%v_eps
    write(lun,*) "G_eps ", baro%G_eps
    call io_close(lun)

  end subroutine write_baro_checkpoint
  !!***

  !!****m* md_control/read_baro_checkpoint *
  !!  NAME
  !!   read_thermo_checkpoint
  !!  PURPOSE
  !!   Read checkpoint for restarting MTTK barostat
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/01/19 13:16
  !!  SOURCE
  !!  
  subroutine read_baro_checkpoint(baro)
    
    use input_module,     only: io_assign, io_close

    ! passed variables
    class(type_barostat), intent(inout) :: baro

    ! local variables
    integer                               :: lun
    character(20)                         :: name

    call io_assign(lun)
    open(unit=lun,file=baro_check_file,status='old')

    select case(baro%baro_type)
    case('iso-mttk')
      read(lun,*) name, baro%lat_ref(1,:)
      read(lun,*) name, baro%lat_ref(2,:)
      read(lun,*) name, baro%lat_ref(3,:)
      read(lun,*) name, baro%volume_ref
      read(lun,*) name, baro%eps_ref
      read(lun,*) name, baro%eps
      read(lun,*) name, baro%v_eps
      read(lun,*) name, baro%G_eps
    case('ortho-mttk')
    case('mttk')
    end select
    call io_close(lun)

  end subroutine read_baro_checkpoint
  !!***

end module md_control
