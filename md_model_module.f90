!!****h* Conquest/md_model
!!  NAME
!!   md_model
!!  PURPOSE
!!   Store the current state of the MD run 
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2017/10/31 11:50
!!  MODIFICATION HISTORY
!!   
!!  SOURCE
!!  
module md_model

  use datatypes
  use numbers
  use force_module,     only: tot_force, stress
  use global_module,    only: ni_in_cell, io_lun, atom_coord
  use species_module,   only: species
  use md_control,       only: md_n_nhc, ion_velocity, type_thermostat, &
                              type_barostat, lattice_vec

  implicit none

!!****s* md_model/type_md_model
!!  NAME
!!   type_md_model
!!  PURPOSE
!!   Container for one step of MD data 
!!  AUTHOR
!!   Zamaan Raza
!!  SOURCE
!!  
  type type_md_model

    logical                                 :: append   ! append to file?

    ! Simulation parameters
    integer                                 :: step
    integer                                 :: ndof     ! degrees of freedom
    character(3)                            :: ensemble ! nve, nvt, npt etc
    real(double), dimension(:,:), pointer   :: lattice_vec
    real(double), pointer                   :: volume

    ! ionic variables
    integer                                 :: natoms
    integer, pointer, dimension(:)          :: species
    real(double), pointer, dimension(:,:)   :: atom_coords
    real(double), pointer, dimension(:,:)   :: atom_velocity
    real(double), pointer, dimension(:,:)   :: atom_force

    ! MD variables
    real(double)                            :: ion_kinetic_energy
    real(double)                            :: dft_total_energy
    real(double)                            :: h_prime  ! conserved qty

    ! Thermodynamic variables
    real(double), pointer                   :: T_int    ! internal temperature
    real(double), pointer                   :: T_ext    ! internal temperature
    real(double), pointer                   :: P_int    ! internal pressure
    real(double), pointer                   :: P_ext    ! internal pressure
    real(double), pointer                   :: PV
    real(double)                            :: enthalpy

    ! Thermostat
    character(20), pointer                  :: thermo_type
    real(double), pointer                   :: lambda   ! velocity scaling fac
    real(double), pointer                   :: tau_T    ! T coupling period
    integer, pointer                        :: n_nhc
    real(double), pointer                   :: nhc_energy
    real(double), pointer                   :: nhc_cell_energy
    real(double), pointer                   :: nhc_ion_energy
    real(double), pointer, dimension(:)     :: eta
    real(double), pointer, dimension(:)     :: v_eta
    real(double), pointer, dimension(:)     :: G_nhc
    real(double), pointer, dimension(:)     :: m_nhc
    real(double), pointer, dimension(:)     :: eta_cell
    real(double), pointer, dimension(:)     :: v_eta_cell
    real(double), pointer, dimension(:)     :: G_nhc_cell
    real(double), pointer, dimension(:)     :: m_nhc_cell

    ! Barostat
    character(20), pointer                  :: baro_type
    real(double), pointer, dimension(:)     :: stress
    real(double), pointer, dimension(:,:)   :: static_stress
    real(double), pointer, dimension(:,:)   :: ke_stress
    real(double), pointer                   :: box_kinetic_energy
    real(double), pointer                   :: m_box
    real(double), pointer                   :: eps
    real(double), pointer                   :: v_eps
    real(double), pointer                   :: G_eps
    real(double), dimension(3,3)            :: c_g
    real(double), dimension(3,3)            :: v_g

    contains

      procedure, public   :: init_model
      procedure, public   :: get_cons_qty
      procedure, public   :: dump_stats
      procedure, public   :: dump_frame

      procedure, private  :: dump_mdl_atom_arr

  end type type_md_model
  !!***

contains

  !!****m* md_model/init_model *
  !!  NAME
  !!   init_model 
  !!  PURPOSE
  !!   Initialiase the MD model 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  SOURCE
  !!  
  subroutine init_model(mdl, ensemble, thermo, baro)

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
    character(3), intent(in)              :: ensemble
    type(type_thermostat), intent(in), target :: thermo
    type(type_barostat), intent(in), target   :: baro

    mdl%append = .false.

    ! General MD arrays
    mdl%natoms = ni_in_cell
    mdl%ensemble = ensemble
    mdl%species       => species
    mdl%atom_coords   => atom_coord
    mdl%atom_force    => tot_force
    mdl%atom_velocity => ion_velocity
    mdl%lattice_vec   => lattice_vec
    mdl%stress        => stress

    ! Thermostat
    mdl%T_int         => thermo%T_int
    mdl%T_ext         => thermo%T_ext
    mdl%thermo_type   => thermo%thermo_type
    mdl%lambda        => thermo%lambda
    mdl%tau_T         => thermo%tau_T
    mdl%n_nhc         => thermo%n_nhc
    mdl%nhc_energy    => thermo%e_nhc
    mdl%nhc_ion_energy  => thermo%e_nhc_ion
    mdl%nhc_cell_energy => thermo%e_nhc_cell
    mdl%eta           => thermo%eta
    mdl%v_eta         => thermo%v_eta
    mdl%G_nhc         => thermo%G_nhc
    mdl%m_nhc         => thermo%m_nhc
    mdl%eta_cell      => thermo%eta_cell
    mdl%v_eta_cell    => thermo%v_eta_cell
    mdl%G_nhc_cell    => thermo%G_nhc_cell
    mdl%m_nhc_cell    => thermo%m_nhc_cell

    ! Barostat
    mdl%P_int         => baro%P_int
    mdl%P_ext         => baro%P_ext
    mdl%volume        => baro%volume
    mdl%PV            => baro%PV
    mdl%box_kinetic_energy => baro%ke_box
    mdl%baro_type     => baro%baro_type
    mdl%static_stress => baro%static_stress
    mdl%ke_stress     => baro%ke_stress
    mdl%m_box         => baro%box_mass
    mdl%eps           => baro%eps
    mdl%v_eps         => baro%v_eps
    mdl%G_eps         => baro%G_eps

  end subroutine init_model
  !!***

  !!****m* md_model/get_cons_qty *
  !!  NAME
  !!   get_cons_qty 
  !!  PURPOSE
  !!   Initialiase the MD model 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  SOURCE
  !!  
  subroutine get_cons_qty(mdl)

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl

    select case(mdl%ensemble)
    case("nve")
      mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy
    case("nvt")
      if (mdl%thermo_type == 'nhc') then
        mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + &
                      mdl%nhc_energy
      end if
    case("npt")
      if (mdl%baro_type == 'iso-mttk' .or. mdl%baro_type == 'mttk') then
        mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + &
                      mdl%nhc_energy + mdl%box_kinetic_energy + mdl%PV
      end if
    end select

  end subroutine get_cons_qty
  !!***

  !!****m* md_model/dump_stats *
  !!  NAME
  !!   dump_stats 
  !!  PURPOSE
  !!   dump thermodynamics stats to a file 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  SOURCE
  !!  
  subroutine dump_stats(mdl, filename)

    use input_module,     only: io_assign, io_close
    use GenComms,         only: inode, ionode
    use md_control,       only: fac_HaBohr32GPa

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
    character(len=*), intent(in)          :: filename

    ! local variables
    integer                               :: lun
    real(double)                          :: P_GPa

    ! Convert units if necessary
    P_GPa = mdl%P_int*fac_HaBohr32GPA

    if (inode==ionode) then
      call io_assign(lun)

      if (mdl%append) then
        open(unit=lun,file=filename,position='append')
      else 
        open(unit=lun,file=filename,status='replace')
        select case (mdl%ensemble)
        case ('nve')
          write(lun,'(a10,3a18,2a12)') "step", "pe", "ke", "H'", "T", "P"
        case ('nvt')
          if (mdl%thermo_type == 'nhc') then
            write(lun,'(a10,4a18,2a12)') "step", "pe", "ke", "nhc", "H'", "T", "P"
          else
            write(lun,'(a10,3a16,2a12)') "step", "pe", "ke", "total", "T", "P"
          end if
        case ('npt')
          if (mdl%thermo_type == 'nhc') then
            write(lun,'(a10,6a18,3a12)') "step", "pe", "ke", "nhc", "box", "pV", "H'", "T", "P", "V"
          end if
        case ('nph')
          write(lun,'(a10,5a18,3a12)') "step", "pe", "ke", "box", "pV", "H'", "T", "P", "V"
        end select
        mdl%append = .true.
      end if
      select case (mdl%ensemble)
      case ('nve')
        write(lun,'(i10,3e18.8,2f12.4)') mdl%step, mdl%dft_total_energy, &
          mdl%ion_kinetic_energy, mdl%h_prime, mdl%T_int, P_GPa
      case ('nvt')
        if (mdl%thermo_type == 'nhc') then
          write(lun,'(i10,4e18.8,2f12.4)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%nhc_energy, mdl%h_prime, mdl%T_int, &
            P_GPa
        else
          write(lun,'(i10,3e18.8,2f12.4)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%h_prime, mdl%T_int, P_GPa
        end if
      case ('npt')
        if (mdl%thermo_type == 'nhc') then
          write(lun,'(i10,6e18.8,3f12.4)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%nhc_energy, mdl%box_kinetic_energy, &
            mdl%PV, mdl%h_prime, mdl%T_int, P_GPa, mdl%volume
        end if
      case ('nph')
          write(lun,'(i10,5e18.8,3f12.4)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%box_kinetic_energy, mdl%PV, &
            mdl%h_prime, mdl%T_int, P_GPa, mdl%volume
      end select
      call io_close(lun)
    end if

  end subroutine dump_stats
  !!***

  !!****m* md_model/dump_frame *
  !!  NAME
  !!   dump_frame
  !!  PURPOSE
  !!   Dump all relevant restart/analysis-relevant data to file 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  SOURCE
  !!  
  subroutine dump_frame(mdl, filename)

    use input_module,     only: io_assign, io_close
    use GenComms,         only: inode, ionode

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
    character(len=*), intent(in)          :: filename

    ! local variables
    integer                               :: lun, i


    if (inode==ionode) then
      call io_assign(lun)
      if (mdl%append) then
        open(unit=lun,file=filename,position='append')
      else 
        open(unit=lun,file=filename,status='replace')
        mdl%append = .true.
      end if

      write(lun,'("frame ",i8)') mdl%step
      write(lun,'(a)') "cell_vectors"
      do i=1,3
        write(lun,'(3f12.6)') mdl%lattice_vec(i,:)
      end do
      write(lun,'(a)') "end cell_vectors"
      write(lun,'(a)') "stress_tensor"
      write(lun,'(3e12.6)') mdl%stress(1), zero, zero
      write(lun,'(3e12.6)') zero, mdl%stress(2), zero
      write(lun,'(3e12.6)') zero, zero, mdl%stress(3)
      write(lun,'(a)') "end stress_tensor"
      write(lun,'(a)') "positions"
      call mdl%dump_mdl_atom_arr(lun, mdl%atom_coords)
      write(lun,'(a)') "end positions"
      write(lun,'(a)') "velocities"
      call mdl%dump_mdl_atom_arr(lun, mdl%atom_velocity)
      write(lun,'(a)') "end velocities"
      write(lun,'(a)') "forces"
      call mdl%dump_mdl_atom_arr(lun, mdl%atom_force)
      write(lun,'(a)') "end forces"
      write(lun,'(a)') "end frame"
    end if

    call io_close(lun)

  end subroutine dump_frame
  !!***

  !!****m* md_model/dump_mdl_atom_arr *
  !!  NAME
  !!   dump_mdl_atom_arr 
  !!  PURPOSE
  !!   dump an array of atomic data to file (positions, velocities, forces) 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  SOURCE
  !!  
  subroutine dump_mdl_atom_arr(mdl, lun, arr)

    ! passed variables
    class(type_md_model), intent(inout)       :: mdl
    integer, intent(in)                       :: lun
    real(double), dimension(:,:), intent(in)  :: arr

    ! local variables
    integer                           :: i

    do i=1,mdl%natoms
      write(lun,'(2i5,3e20.10)') i, mdl%species(i), arr(:,i)
    end do

  end subroutine dump_mdl_atom_arr
  !!***

end module md_model
