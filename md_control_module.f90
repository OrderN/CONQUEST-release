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
!!   2017/10/24 09:24 zamaan
!!    Refactored original NVT code from NVTMD branch
!!  SOURCE
!!
module md_control

  use datatypes
  use numbers
  use global_module,    only: ni_in_cell, io_lun, iprint_MD
  use move_atoms,       only: fac_Kelvin2Hartree
  use species_module,   only: species, mass
  use GenComms,         only: inode, ionode

  implicit none

  ! Unit conversion factors
  real(double), parameter :: fac_HaBohr32GPa = 29421.02648438959

  character(20) :: md_thermo_type, md_baro_type
  real(double)  :: md_tau_T, md_tau_P, md_target_press, md_baro_beta, &
                   md_box_mass
  integer       :: md_n_nhc, md_n_ys, md_n_mts
  logical       :: md_write_xsf
  real(double), dimension(3,3), target      :: lattice_vec
  real(double), dimension(:), allocatable   :: md_nhc_mass
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
    real(double)        :: T_int        ! instantateous temperature
    real(double)        :: T_ext        ! target temperature
    real(double)        :: ke_ions      ! kinetic energy of ions
    real(double)        :: dt           ! time step
    integer             :: ndof         ! number of degrees of freedom
    logical             :: append
    real(double)        :: lambda       ! velocity scaling factor

    ! Weak coupling thermostat variables
    real(double)        :: tau_T        ! temperature coupling time period

    ! Nose-Hoover chain thermostat variables
    integer             :: n_nhc        ! number of Nose-Hoover heat baths
    integer             :: n_ys         ! Yoshida-Suzuki order
    integer             :: n_mts_nhc    ! number of time steps for NHC
    real(double)        :: ke_nhc       ! kinetic energy of NHC thermostats
    real(double), dimension(:), allocatable :: eta    ! thermostat "position"
    real(double), dimension(:), allocatable :: v_eta  ! thermostat "velocity"
    real(double), dimension(:), allocatable :: G_nhc  ! "force" on thermostat
    real(double), dimension(:), allocatable :: m_nhc  ! thermostat mass
    real(double), dimension(:), allocatable :: dt_ys  ! Yoshida-Suzuki time steps

    contains

      procedure, public   :: init_thermo_none
      procedure, public   :: init_nhc
      procedure, public   :: init_berendsen
      procedure, public   :: update_ke_ions
      procedure, public   :: get_berendsen_sf
      procedure, public   :: berendsen_v_rescale
      procedure, public   :: get_nhc_energy
      procedure, public   :: get_temperature
      procedure, public   :: dump_thermo_state
      procedure, public   :: propagate_nvt_nhc

      procedure, private  :: update_G_eta
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
    real(double)        :: beta         ! isothermal compressibility

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

      procedure, public   :: init_baro_none
      procedure, public   :: init_baro_mttk
      procedure, public   :: update_static_stress
      procedure, public   :: get_pressure
      procedure, public   :: get_volume
      procedure, public   :: get_box_ke
      procedure, public   :: get_ke_stress
      procedure, public   :: propagate_npt_mttk
      procedure, public   :: propagate_r_mttk
      procedure, public   :: propagate_box_mttk
      procedure, public   :: dump_baro_state
      procedure, public   :: update_cell

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

  !!****m* md_control/init_thermo_none *
  !!  NAME
  !!   init_thermo_none
  !!  PURPOSE
  !!   initialise thermostat just for calculating temperature
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/09 10:11
  !!  SOURCE
  !!  
  subroutine init_thermo_none(th, ndof, ke_ions)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: ndof
    real(double), intent(in)              :: ke_ions

    th%thermo_type = "None"
    th%ndof = ndof
    th%ke_ions = ke_ions
    call th%get_temperature

  end subroutine init_thermo_none
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
  subroutine init_nhc(th, dt, T_ext, ndof, n_nhc, n_ys, n_mts, ke_ions)

    use GenComms,         only: cq_abort

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: ndof, n_nhc, n_ys, n_mts
    real(double), intent(in)              :: T_ext, ke_ions, dt

    ! local variables
    character(40)                         :: fmt

    th%thermo_type = "nhc"
    th%append = .false.
    th%dt = dt
    th%T_ext = T_ext
    th%ndof = ndof
    th%n_nhc = n_nhc
    th%ke_ions = ke_ions
    th%n_ys = n_ys
    th%n_mts_nhc = n_mts
    th%lambda = one
    call th%get_temperature

    allocate(th%eta(n_nhc))
    allocate(th%v_eta(n_nhc))
    allocate(th%G_nhc(n_nhc))
    allocate(th%m_nhc(n_nhc))
    allocate(th%dt_ys(n_ys))

    ! Defaults for heat bath positions, velocities, masses
    th%eta = zero
    th%m_nhc = one
    th%v_eta = sqrt(two*th%T_ext*fac_Kelvin2Hartree/th%m_nhc(1)) 
    th%G_nhc = zero

    ! Yoshida-Suzuki time steps
    select case(th%n_ys) ! The YS order
    case(1)
      th%dt_ys(1) = one
    case(3)              ! n_ys = 3,5 taken from MTTK; higher orders possible
      th%dt_ys(1) = one/(two - two**(one/three))
      th%dt_ys(2) = one - two*th%dt_ys(1)
      th%dt_ys(3) = th%dt_ys(1)
    case(5)
      th%dt_ys(1) = one/(four - four**(one/three))
      th%dt_ys(2) = th%dt_ys(1)
      th%dt_ys(3) = one - four*th%dt_ys(1)
      th%dt_ys(4) = th%dt_ys(1)
      th%dt_ys(5) = th%dt_ys(1)
    case default
      call cq_abort("Invalid Yoshida-Suzuki order")
    end select
    th%dt_ys = th%dt*th%dt_ys/th%n_mts_nhc

    write(fmt,'("(4x,a16,",i4,"f10.4)")') th%n_nhc
    if (inode==ionode) then
      write(io_lun,'(2x,a)') 'Welcome to init_nhc'
      write(io_lun,'(4x,a,f10.2)') 'Target temperature        T_ext = ', &
                                   th%T_ext
      write(io_lun,'(4x,a,f10.2)') 'Instantaneous temperature T_int = ', &
                                   th%T_int
      write(io_lun,'(4x,a,i10)')   'Number of NHC thermostats n_nhc = ', &
                                   th%n_nhc
      write(io_lun,'(4x,a,i10)')   'Multiple time step order  n_mts = ', &
                                   th%n_mts_nhc
      write(io_lun,'(4x,a,i10)')   'Yoshida-Suzuki order      n_ys  = ', &
                                   th%n_ys
      write(io_lun,fmt) 'NHC masses:    ', th%m_nhc
      write(io_lun,fmt) 'YS time steps: ', th%dt_ys
    end if

  end subroutine init_nhc
  !!***

  !!****m* md_control/init_berendsen *
  !!  NAME
  !!   init_berendsen
  !!  PURPOSE
  !!   initialise Berendsen weak coupling thermostat variables
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 10:39
  !!  SOURCE
  !!  
  subroutine init_berendsen(th, dt, T_ext, ndof, tau_T, ke_ions)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: ndof
    real(double), intent(in)              :: T_ext, ke_ions, dt, tau_T

    th%thermo_type = "berendsen"
    th%dt = dt
    th%T_ext = T_ext
    th%ndof = ndof
    th%ke_ions = ke_ions
    th%tau_T = tau_T

  end subroutine init_berendsen
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

  !!****m* md_control/get_berendsen_sf *
  !!  NAME
  !!   get_berendsen_sf
  !!  PURPOSE
  !!   Get velocity scaling factor for Berendsen thermostat
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 14:26
  !!  SOURCE
  !!  
  subroutine get_berendsen_sf(th)

    ! passed variables
    class(type_thermostat), intent(inout)   :: th

    th%lambda = sqrt(one + (th%dt/th%tau_T)*(th%T_ext*fac_Kelvin2Hartree/th%T_int - one))

  end subroutine get_berendsen_sf
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

  end subroutine berendsen_v_rescale
  !!***

  !!****m* md_control/update_G_eta *
  !!  NAME
  !!   update_G_eta
  !!  PURPOSE
  !!   updates the "force" on thermostat k 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  CREATION DATE
  !!   2017/10/24 10:41
  !!  SOURCE
  !!  
  subroutine update_G_eta(th, k, ke_box)

    ! passed variables
    class(type_thermostat), intent(inout) :: th
    integer, intent(in)                   :: k      ! NH thermostat index
    real(double), intent(in)              :: ke_box ! box ke for pressure coupling

    if (k == 1) then
      th%G_nhc(k) = two*th%ke_ions - th%ndof*th%T_ext*fac_Kelvin2Hartree + &
                    two*ke_box
    else
      th%G_nhc(k) = th%m_nhc(k-1)*th%v_eta(k-1)**2 - &
                    th%T_ext*fac_Kelvin2Hartree
    end if
    th%G_nhc(k) = th%G_nhc(k)/th%m_nhc(k)

  end subroutine update_G_eta
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

    if (inode==ionode .and. iprint_MD>0) write(io_lun,*) "Welcome to &
                                                          propagate_nvt_nhc"

    th%ke_ions = ke
    v_sfac = one
    call th%update_G_eta(1, zero)
    do i_mts=1,th%n_mts_nhc ! MTS loop
      do i_ys=1,th%n_ys     ! Yoshida-Suzuki loop
        ! Reverse part of Trotter expansion: update thermostat force/velocity
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
        do i_nhc=th%n_nhc-1,1,-1 ! loop over NH thermostats in reverse order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do

        ! scale the ionic velocities and kinetic energy
        fac = exp(-half*th%dt_ys(i_ys)*th%v_eta(1))
        v_sfac = v_sfac*fac
        if (inode==ionode .and. iprint_MD > 1) write(io_lun,*) 'v_sfac = ', v_sfac
        th%ke_ions = th%ke_ions*v_sfac**2

        call th%update_G_eta(1, zero)
        ! update the thermostat "positions" eta
        do i_nhc=1,th%n_nhc
          call th%propagate_eta(i_nhc, th%dt_ys(i_ys), half)
        end do

        ! Forward part of Trotter expansion: update thermostat force/velocity
        do i_nhc=1,th%n_nhc-1 ! loop over NH thermostats in forward order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          if (i_nhc /= 1) call th%update_G_eta(i_nhc, zero)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do
        call th%update_G_eta(th%n_nhc, zero) ! box ke is zero in NVT ensemble
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
      end do  ! Yoshida-Suzuki loop
    end do    ! MTS loop

    ! scale the ionic velocities
    v = v_sfac*v
    th%lambda = v_sfac

  end subroutine propagate_nvt_nhc
  !!***

  !!****m* md_control/get_nhc_energy *
  !!  NAME
  !!   get_nhc_energy
  !!  PURPOSE
  !!   compute the NHC contribution to the constant of motion
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2017/10/24 13:19
  !!  SOURCE
  !!  
  subroutine get_nhc_energy(th)

    ! passed variables
    class(type_thermostat), intent(inout) :: th

    th%ke_nhc = half*sum(th%m_nhc*th%v_eta**2) ! TODO: check this
    th%ke_nhc = th%ke_nhc + th%ndof*th%T_ext*fac_Kelvin2Hartree*th%eta(1)
    th%ke_nhc = th%ke_nhc + th%T_ext*fac_Kelvin2Hartree*sum(th%eta(2:))

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

    th%T_int = two*th%ke_ions/fac_Kelvin2Hartree/th%ndof
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
    character(40)                         :: fmt
    integer                               :: lun

    if (inode==ionode) then
      call io_assign(lun)
      if (th%append) then
        open(unit=lun,file=filename,position='append')
      else 
        open(unit=lun,file=filename,status='replace')
        th%append = .true.
      end if
      write(fmt,'("(a8,",i4,"e16.4)")') th%n_nhc
      write(lun,'("step    ",i12)') step
      write(lun,'("T_int   ",f12.4)') th%T_int
      write(lun,'("ke_ions ",e12.4)') th%ke_ions
      write(lun,'("lambda  ",f12.4)') th%lambda
      if (th%thermo_type == 'nhc') then
        write(lun,fmt) "eta:    ", th%eta
        write(lun,fmt) "v_eta:  ", th%v_eta
        write(lun,fmt) "G_nhc:  ", th%G_nhc
        write(lun,'("e_nhc: ",e16.4)') th%ke_nhc
      else if (th%thermo_type == 'berendsen') then
        write(lun,'("lambda: ",e16.4)') th%lambda
      end if
      write(lun,*)
      call io_close(lun)
    end if

  end subroutine dump_thermo_state
  !!***

  !!****m* md_control/init_baro_none *
  !!  NAME
  !!   init_baro_none
  !!  PURPOSE
  !!   initialise barostat for constant volume, just to compute pressure
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/08 12:23
  !!  SOURCE
  !!  
  subroutine init_baro_none(baro, stress)

    ! passed variables
    class(type_barostat), intent(inout)     :: baro
    real(double), dimension(3), intent(in)  :: stress

    ! local variables
    integer                                 :: i

    baro%baro_type = 'None'
    baro%static_stress = zero
    do i=1,3
      baro%static_stress(i,i) = stress(i)
    end do

  end subroutine init_baro_none
  !!***

  !!****m* md_control/init_baro_mttk *
  !!  NAME
  !!   init_baro_mttk
  !!  PURPOSE
  !!   initialise MTTK barostat
  !!  AUTHOR
  !!    Zamaan Raza 
  !!  CREATION DATE
  !!   2017/11/17 12:44
  !!  SOURCE
  !!  
  subroutine init_baro_mttk(baro, P_ext, ndof, stress, v, ke_ions)

    use GenComms,         only: cq_abort
    use global_module,    only: rcellx, rcelly, rcellz

    ! passed variables
    class(type_barostat), intent(inout)       :: baro
    real(double), intent(in)                  :: P_ext, ke_ions
    integer, intent(in)                       :: ndof
    real(double), dimension(3), intent(in)    :: stress
    real(double), dimension(:,:), intent(in)  :: v

    ! Globals
    baro%baro_type = md_baro_type
    baro%box_mass = md_box_mass
    baro%tau_P = md_tau_P
    baro%beta = md_baro_beta

    ! constants for polynomial expanion
    baro%c2 = one/6.0_double
    baro%c4 = baro%c2/20.0_double
    baro%c6 = baro%c4/42.0_double
    baro%c8 = baro%c6/72.0_double

    baro%P_ext = P_ext
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

    select case(baro%baro_type)
    case('iso-mttk')
      baro%volume_ref = baro%volume
      baro%eps_ref = third*log(baro%volume/baro%volume_ref)
      baro%eps = baro%eps_ref
      baro%v_eps = zero
      baro%ke_box = zero
      baro%odnf = one + three/baro%ndof
    case default
      call cq_abort("Invalid barostat")
    end select

    if (inode==ionode) then
      write(io_lun,('(2x,a)')) 'Welcome to init_baro_mttk'
      write(io_lun,'(4x,a,f10.2)') 'Target pressure        P_ext = ', &
                                    baro%P_ext
      write(io_lun,'(4x,a,f10.2)') 'Instantaneous pressure P_int = ', &
                                    baro%P_int
      write(io_lun,'(4x,a,f10.2)') 'Box mass                     = ', &
                                    baro%box_mass
    end if

  end subroutine init_baro_mttk
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
    baro%P_int = baro%P_int*third

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

  end subroutine get_volume
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
      baro%ke_box = half*baro%v_eps**2/baro%box_mass 
    end select

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

    baro%G_eps = (baro%odnf*baro%ke_ions + three*(baro%P_int - baro%P_ext)*baro%volume)/baro%box_mass
    ! baro%G_eps = (three*baro%ke_ions/baro%ndof + three*(baro%P_int - baro%P_ext)*baro%volume)/baro%box_mass

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
  subroutine propagate_r_mttk(baro, dt, dtfac, v, flag_movable)

    use global_module,      only: x_atom_cell, y_atom_cell, z_atom_cell, &
                                  atom_coord, atom_coord_diff, id_glob
    use species_module,     only: species, mass
    use move_atoms,         only: fac

    ! passed variables
    class(type_barostat), intent(inout)   :: baro
    real(double), intent(in)              :: dt     ! time step
    real(double), intent(in)              :: dtfac  ! Trotter epxansion factor
    real(double), dimension(:,:), intent(in)    :: v
    logical, dimension(:), intent(in)           :: flag_movable

    ! local variables
    real(double)                          :: exp_v_eps, sinhx_x, fac_r, &
                                             fac_v, massa, x_old, y_old, z_old
    integer                               :: i, speca, gatom, ibeg_atom
    logical                               :: flagx, flagy, flagz

    select case(baro%baro_type)
    case('iso-mttk')
      exp_v_eps = exp(dt*dtfac*baro%v_eps)
      sinhx_x = baro%poly_sinhx_x(dtfac*dt*baro%v_eps)
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
    end select

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

    select case(baro%baro_type)
    case('iso-mttk')
      call baro%propagate_eps_lin(dt, one)
      v_old = baro%volume
      v_new = baro%volume_ref*exp(three*baro%eps)
      lat_sfac = (v_new/v_old)**third
      baro%mu = lat_sfac
      baro%lat = baro%lat*lat_sfac
      if (inode==ionode .and. iprint_MD>1) then
        write(io_lun,*) 'lat_sfac = ', lat_sfac
      end if
    end select

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
  subroutine propagate_npt_mttk(baro, th, stress, ke, v)

    ! passed variables
    class(type_barostat), intent(inout)     :: baro
    type(type_thermostat), intent(inout)    :: th
    real(double), intent(in)                :: ke
    real(double), dimension(3), intent(in)  :: stress
    real(double), dimension(3,3), intent(inout) :: v

    ! local variables
    integer                                 :: i_mts, i_ys, i_nhc
    real(double)                            :: v_sfac

    if (inode==ionode .and. iprint_MD>0) write(io_lun,*) "Welcome to &
                                                          propagate_npt_mttk"

    baro%ke_ions = ke
    v_sfac = one
    call baro%get_box_ke
    call th%update_G_eta(1, baro%ke_box)
    call baro%update_G_eps

    do i_mts=1,th%n_mts_nhc ! MTS loop
      do i_ys=1,th%n_ys     ! Yoshida-Suzuki loop
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
        do i_nhc=th%n_nhc-1,1,-1 ! loop over NH thermostats in reverse order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do

        ! update box velocities
        select case(baro%baro_type)
        case('iso-mttk')
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, th%v_eta(1))
          call baro%propagate_v_eps_lin(th%dt_ys(i_ys), quarter)
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, th%v_eta(1))
        end select

        ! update ionic velocities, scale ion kinetic energy
        call baro%update_vscale_fac(th%dt_ys(i_ys), half, th%v_eta(1), v_sfac)
        th%ke_ions = baro%ke_ions

        call baro%update_G_eps
        call baro%get_box_ke
        call th%update_G_eta(1, baro%ke_box)
        ! update the thermostat "positions" eta
        do i_nhc=1,th%n_nhc
          call th%propagate_eta(i_nhc, th%dt_ys(i_ys), half)
        end do

        ! update box velocities
        select case(baro%baro_type)
        case('iso-mttk')
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, th%v_eta(1))
          call baro%propagate_v_eps_lin(th%dt_ys(i_ys), quarter)
          call baro%propagate_v_eps_exp(th%dt_ys(i_ys), one_eighth, th%v_eta(1))
        end select

        do i_nhc=1,th%n_nhc-1 ! loop over NH thermostats in forward order
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
          if (i_nhc > 1) call th%update_G_eta(i_nhc, baro%ke_box)
          call th%propagate_v_eta_lin(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_exp(i_nhc, th%dt_ys(i_ys), one_eighth)
        end do
        call th%update_G_eta(th%n_nhc, baro%ke_box)
        call th%propagate_v_eta_lin(th%n_nhc, th%dt_ys(i_ys), quarter)
      end do  ! Yoshida-Suzuki loop
    end do    ! MTS loop

    ! scale the velocities
    v = v_sfac*v
    th%lambda = v_sfac

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
      write(lun,'("P_int   ",f12.4)') baro%P_int
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
    use global_module,      only: iprint_MD, rcellx, rcelly, rcellz, &
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

    orcellx = rcellx
    orcelly = rcelly
    orcellz = rcellz

    rcellx = baro%lat(1,1)
    rcelly = baro%lat(2,2)
    rcellz = baro%lat(3,3)

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

end module md_control
