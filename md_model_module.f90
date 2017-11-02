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
  use force_module,     only: tot_force
  use global_module,    only: ni_in_cell, io_lun, atom_coord
  use species_module,   only: species
  use md_control,       only: md_n_nhc, ion_velocity, type_thermostat

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
    real(double), dimension(3,3)            :: lattice_vec
    real(double)                            :: volume

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
    real(double)                            :: T_int    ! internal temperature
    real(double)                            :: P_int    ! internal pressure
    real(double)                            :: PV
    real(double)                            :: enthalpy

    ! Thermostat
    character(20), pointer                  :: thermo_type
    real(double), pointer                   :: lambda   ! velocity scaling fac
    real(double), pointer                   :: tau_T    ! T coupling period
    integer, pointer                        :: n_nhc
    real(double), pointer                   :: nhc_energy
    real(double), pointer, dimension(:)     :: eta
    real(double), pointer, dimension(:)     :: v_eta
    real(double), pointer, dimension(:)     :: G_nhc
    real(double), pointer, dimension(:)     :: m_nhc

    ! Barostat
    real(double), dimension(3,3)            :: stress
    real(double)                            :: box_kinetic_energy
    real(double)                            :: m_box
    real(double)                            :: eps
    real(double)                            :: v_eps
    real(double)                            :: G_eps
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
  subroutine init_model(mdl, thermo)

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
    type(type_thermostat), intent(in), target  :: thermo

    mdl%append = .false.

    ! General MD arrays
    mdl%natoms = ni_in_cell
    mdl%species       => species
    mdl%atom_coords   => atom_coord
    mdl%atom_force    => tot_force
    mdl%atom_velocity => ion_velocity

    ! NHC thermostat arrays
    mdl%thermo_type   => thermo%thermo_type
    mdl%lambda        => thermo%lambda
    mdl%tau_T         => thermo%tau_T
    mdl%n_nhc         => thermo%n_nhc
    mdl%nhc_energy    => thermo%ke_nhc
    mdl%eta           => thermo%eta
    mdl%v_eta         => thermo%v_eta
    mdl%G_nhc         => thermo%G_nhc
    mdl%m_nhc         => thermo%m_nhc

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
      mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + &
                    mdl%nhc_energy
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

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
    character(len=*), intent(in)          :: filename

    ! local variables
    integer                               :: lun


    if (inode==ionode) then
      call io_assign(lun)
      if (mdl%append) then
        open(unit=lun,file=filename,position='append')
      else 
        open(unit=lun,file=filename,status='replace')
        mdl%append = .true.
      end if
      if (mdl%step == 0) then
        select case (mdl%ensemble)
        case ('nve')
          write(lun,'(a10,5a16)') "step", "pe", "ke", "H'", "T", "P"
        case ('nvt')
          if (mdl%thermo_type == 'nhc') then
            write(lun,'(a10,6a16)') "step", "pe", "ke", "nhc", "H'", "T", "P"
          else
            write(lun,'(a10,5a16)') "step", "pe", "ke", "total", "T", "P"
          end if
        case ('npt')
          if (mdl%thermo_type == 'nhc') then
            write(lun,'(a10,9a16)') "step", "pe", "ke", "nhc", "box", "pV", "H'", "T", "P", "V"
          end if
        case ('nph')
          write(lun,'(a10,8a16)') "step", "pe", "ke", "box", "pV", "H'", "T", "P", "V"
        end select
      end if
      select case (mdl%ensemble)
      case ('nve')
        write(lun,'(i10,5e16.6)') mdl%step, mdl%dft_total_energy, &
          mdl%ion_kinetic_energy, mdl%h_prime, mdl%T_int, mdl%P_int
      case ('nvt')
        if (mdl%thermo_type == 'nhc') then
          write(lun,'(i10,6e16.6)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%nhc_energy, mdl%h_prime, mdl%T_int, &
            mdl%P_int
        else
          write(lun,'(i10,5e16.6)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%h_prime, mdl%T_int, mdl%P_int
        end if
      case ('npt')
        if (mdl%thermo_type == 'nhc') then
          write(lun,'(i10,9e16.6)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%nhc_energy, mdl%box_kinetic_energy, &
            mdl%pV, mdl%h_prime, mdl%T_int, mdl%P_int, mdl%volume
        end if
      case ('nph')
          write(lun,'(i10,8e16.6)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%box_kinetic_energy, mdl%pV, &
            mdl%h_prime, mdl%T_int, mdl%P_int, mdl%volume
      end select
    end if

    call io_close(lun)

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
    character(20), intent(in)             :: filename

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
      do i=1,3
        write(lun,'(3f12.6)') mdl%stress(i,:)
      end do
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