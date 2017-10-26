! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module thermostat
! ------------------------------------------------------------------------------
! Code area 9: General
! ------------------------------------------------------------------------------

!!****h* Conquest/name
!!  NAME
!!   thermostat
!!  PURPOSE
!!   Controls temperature in molecular dynamics. Algorithm adapted from
!!   G. Martyna et al. Mol. Phys. 5, 1117 (1996) ("MTTK Paper")
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
module thermostat

  use datatypes,        only: double
  use numbers
  use global_module,    only: ni_in_cell, io_lun
  use move_atoms,       only: kB

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"

  character(20) :: md_thermo_type
  real(double)  :: md_tau_T
  integer       :: md_n_nhc, md_n_ys, md_n_mts
  real(double), dimension(:), allocatable :: md_nhc_mass

    !!****s* /
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
    real(double)        :: T_ext        ! target temperature
    real(double)        :: T_int        ! instantaneous temperature
    real(double)        :: ke_ions      ! kinetic energy of ions
    real(double)        :: dt           ! time step
    integer             :: ndof         ! number of degrees of freedom
    logical             :: restart

    ! Weak coupling thermostat variables
    real(double)        :: lambda       ! Berendsen scaling factor
    real(double)        :: tau_T        ! temperature coupling time period

    ! Nose-Hoover chain thermostat variables
    integer             :: n_nhc        ! number of Nose-Hoover heat baths
    integer             :: n_ys         ! Yoshida-Suzuki order
    integer             :: n_mts_nhc    ! number of time steps for NHC
    real(double)        :: ke_nhc       ! kinetic energy of NHC thermostats
    real(double), dimension(:), allocatable   :: eta    ! thermostat "position"
    real(double), dimension(:), allocatable   :: v_eta  ! thermostat "velocity"
    real(double), dimension(:), allocatable   :: G_nhc  ! "force" on thermostat
    real(double), dimension(:), allocatable   :: m_nhc  ! thermostat mass
    real(double), dimension(:), allocatable   :: dt_ys  ! Yoshida-Suzuki time steps

    contains

      procedure :: init_nhc
      procedure :: init_berendsen
      procedure :: get_berendsen_sf
      procedure :: berendsen_v_rescale
      procedure :: update_G
      procedure :: propagate_eta
      procedure :: propagate_v_eta_1
      procedure :: propagate_v_eta_2
      procedure :: propagate_nvt_nhc
      procedure :: get_nhc_energy
      procedure :: get_temperature
      procedure :: dump_thermo_state
  end type type_thermostat
!!***

contains

!!****m* / *
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

  th%thermo_type = "nhc"
  th%restart = .false.
  th%dt = dt
  th%T_ext = T_ext
  th%ndof = ndof
  th%n_nhc = n_nhc
  th%ke_ions = ke_ions
  th%n_ys = n_ys
  th%n_mts_nhc = n_mts

  allocate(th%eta(n_nhc))
  allocate(th%v_eta(n_nhc))
  allocate(th%G_nhc(n_nhc))
  allocate(th%m_nhc(n_nhc))
  allocate(th%dt_ys(n_ys))

  ! Defaults for heat bath positions, velocities, masses
  th%eta = zero
  th%m_nhc = one
  th%v_eta = sqrt(two*th%T_ext/th%m_nhc(1)) 
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

  call th%get_temperature

end subroutine init_nhc
!!***

!!****m* / *
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

!!****m* / *
!!  NAME
!!   berendsen_v_rescale
!!  PURPOSE
!!   thermostat using Berendsen (weak coupling) velocity rescaling
!!  AUTHOR
!!    Zamaan Raza 
!!  CREATION DATE
!!   2017/10/24 14:26
!!  SOURCE
!!  
subroutine get_berendsen_sf(th)

  ! passed variables
  class(type_thermostat), intent(inout)   :: th

  th%lambda = sqrt(one + (th%dt/th%tau_T)*(th%T_ext/th%T_int - one))

end subroutine get_berendsen_sf
!!***

!!****m* / *
!!  NAME
!!   berendsen_v_rescale
!!  PURPOSE
!!   thermostat using Berendsen (weak coupling) velocity rescaling. Note that
!!   the scaling factor is computed before the first vVerlet velocity update,
!!   and the velocities scaled after the second vVerlet velocity update.
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

!!****m* / *
!!  NAME
!!   update_G
!!  PURPOSE
!!   updates the "force" on thermostat k 
!!  AUTHOR
!!   Zamaan Raza 
!!  CREATION DATE
!!   2017/10/24 10:41
!!  SOURCE
!!  
subroutine update_G(th, k, ke_box)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: k      ! NH thermostat index
  real(double), intent(in)              :: ke_box ! box ke for pressure coupling

  if (k == 1) then
    th%G_nhc(k) = 2*th%ke_ions - th%ndof*kB*th%T_ext + ke_box
  else
    th%G_nhc(k) = th%m_nhc(k-1)*th%v_eta(k-1)**2 - kB*th%T_ext
  end if
  th%G_nhc(k) = th%G_nhc(k)/th%m_nhc(k)

end subroutine update_G
!!***

!!****m* / *
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

!!****m* / *
!!  NAME
!!   propagate_v_eta_1
!!  PURPOSE
!!   propagate the "velocity" of thermostat k: linear shift 
!!  AUTHOR
!!   Zamaan Raza  
!!  CREATION DATE
!!   2017/10/24 11:43
!!  SOURCE
!!  
subroutine propagate_v_eta_1(th, k, dt, dtfac)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: k      ! NH thermostat index
  real(double), intent(in)              :: dt     ! time step
  real(double), intent(in)              :: dtfac  ! Trotter expansion factor

  th%v_eta(k) = th%v_eta(k) + dtfac*dt*th%G_nhc(k)

end subroutine propagate_v_eta_1
!!***

!!****m* / *
!!  NAME
!!   propagate_v_eta_2
!!  PURPOSE
!!   propagate the "velocity" of thermostat k: exponential factor 
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2017/10/24 11:45
!!  SOURCE
!!  
subroutine propagate_v_eta_2(th, k, dt, dtfac)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: k      ! NH thermostat index
  real(double), intent(in)              :: dt     ! time step
  real(double), intent(in)              :: dtfac  ! Trotter expansion factor

  th%v_eta(k) = th%v_eta(k)*exp(-dtfac*dt*th%v_eta(k+1))

end subroutine propagate_v_eta_2
!!***

!!****m* / *
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
subroutine propagate_nvt_nhc(th, v)

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  real(double), dimension(:,:), intent(inout) :: v  ! ion velocities

  ! local variables
  integer       :: i_mts, i_ys, i_nhc
  real(double)  :: v_sfac   ! ionic velocity scaling factor
  real(double)  :: fac

  v_sfac = one
  do i_mts=1,th%n_mts_nhc ! MTS loop
    do i_ys=1,th%n_ys     ! Yoshida-Suzuki loop
      ! Reverse part of Trotter expansion: update thermostat force/velocity
      do i_nhc=th%n_nhc,1,-1 ! loop over NH thermostats in reverse order
        if (i_nhc==th%n_nhc) then
          call th%update_G(i_nhc, zero) ! box ke is zero in NVT ensemble
          call th%propagate_v_eta_1(i_nhc, th%dt_ys(i_ys), quarter)
        else
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_2(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%propagate_v_eta_1(i_nhc, th%dt_ys(i_ys), quarter)
          call th%update_G(i_nhc, zero)
          call th%propagate_v_eta_2(i_nhc, th%dt_ys(i_ys), one_eighth)
        end if
      end do

      ! scale the ionic velocities and kinetic energy
      fac = exp(-half*th%dt_ys(i_ys)*th%v_eta(1))
      v_sfac = v_sfac*fac
      th%ke_ions = th%ke_ions*fac**2

      ! update the thermostat "positions" eta
      do i_nhc=1,th%n_nhc
        call th%propagate_eta(i_nhc, th%dt_ys(i_ys), half)
      end do

      ! Forward part of Trotter expansion: update thermostat force/velocity
      do i_nhc=1,th%n_nhc ! loop over NH thermostats in forward order
        if (i_nhc<th%n_nhc) then
          ! Trotter expansion to avoid sinh singularity
          call th%propagate_v_eta_2(i_nhc, th%dt_ys(i_ys), one_eighth)
          call th%update_G(i_nhc, zero)
          call th%propagate_v_eta_1(i_nhc, th%dt_ys(i_ys), quarter)
          call th%propagate_v_eta_2(i_nhc, th%dt_ys(i_ys), one_eighth)
        else
          call th%update_G(i_nhc, zero) ! box ke is zero in NVT ensemble
          call th%propagate_v_eta_1(i_nhc, th%dt_ys(i_ys), quarter)
        end if
      end do
    end do  ! Yoshida-Suzuki loop
  end do    ! MTS loop

  ! scale the ionic velocities
  v = v_sfac*v
  call th%get_temperature

end subroutine propagate_nvt_nhc
!!***

!!****m* / *
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
  th%ke_nhc = th%ke_nhc + th%ndof*kB*th%T_ext*th%eta(1)
  th%ke_nhc = th%ke_nhc + kB*th%T_ext*sum(th%eta(2:))

end subroutine get_nhc_energy
!!***

!!****m* / *
!!  NAME
!!   get_temperature
!!  PURPOSE
!!   compute the istantaneous temperature
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2017/10/24 14:48
!!  SOURCE
!!  
subroutine get_temperature(th)

  ! passed variables
  class(type_thermostat), intent(inout) :: th

  th%T_int = two*th%ke_ions/kB/real(th%ndof, double)

end subroutine get_temperature
!!***

!!****m* / *
!!  NAME
!!   dump_thermo_state
!!  PURPOSE
!!   dump the state of the thermostat to unit lun
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2017/10/24 13:36
!!  SOURCE
!!  
subroutine dump_thermo_state(th, step, filename)
  use input_module,     only: io_assign, io_close
  use GenComms,         only: myid

  ! passed variables
  class(type_thermostat), intent(inout) :: th
  integer, intent(in)                   :: step
  character(len=*), intent(in)          :: filename

  ! local variables
  character(40)                         :: fmt
  integer                               :: lun

  if (myid==0) then
    call io_assign(lun)
    if (th%restart) then
      open(unit=lun,file=filename,status='old',position='append')
      write(lun,*)
    else 
      open(unit=lun,file=filename,status='replace')
    end if
    write(fmt,'("(a8,",i4,"e16.4)")') th%n_nhc
    write(lun,'("step   ",i12)') step
    write(lun,'("T_int  ",f12.4)') th%T_int
    if (th%thermo_type == 'nhc') then
      write(lun,fmt) "eta:    ", th%eta
      write(lun,fmt) "v_eta:  ", th%v_eta
      write(lun,fmt) "G_nhc:  ", th%G_nhc
      write(lun,'("e_nhc:  ",e16.4)') th%ke_nhc
    else if (th%thermo_type == 'berendsen') then
      write(lun,'("lambda: ",e16.4)') th%lambda
    end if
    call io_close(lun)
  end if

end subroutine dump_thermo_state
!!***

end module thermostat
