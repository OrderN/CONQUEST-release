!!****h* Conquest/md_model *
!!  NAME
!!   md_model
!!  PURPOSE
!!   Store the current state of the MD run 
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2017/10/31 11:50
!!  MODIFICATION HISTORY
!!   2018/05/30 zamaan
!!    Modified dump_stats to print correct header and columns when using
!!    Berendsen equilibration in NVT ensemble
!!
!!  SOURCE
!!
module md_model

  use datatypes
  use numbers
  use GenComms,         only: inode, ionode
  use input_module,     only: leqi
  use force_module,     only: tot_force, stress
  use global_module,    only: ni_in_cell, io_lun, atom_coord, iprint_MD, &
                              flag_MDcontinue, flag_MDdebug, x_atom_cell, &
                              y_atom_cell, z_atom_cell, rcellx, rcelly, rcellz
  use species_module,   only: species
  use md_control,       only: md_n_nhc, ion_velocity, type_thermostat, &
                              type_barostat, lattice_vec, &
                              flag_extended_system, md_cell_constraint

  implicit none

  logical       :: md_tdep ! Dump output for a TDEP calculation?
  character(20) :: file_meta      = 'infile.meta' ! filenames for TDEP
  character(20) :: file_positions = 'infile.positions'
  character(20) :: file_forces    = 'infile.forces'
  character(20) :: file_stat      = 'infile.stat'
  character(20) :: file_loto      = 'infile.lotosplitting'

  !!****s* md_model/type_md_model *
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
    real(double)                            :: timestep
    integer                                 :: nequil   ! equilibration steps

    ! ionic variables
    integer                                 :: natoms
    integer, pointer, dimension(:)          :: species
    real(double), pointer                   :: lat_a, lat_b, lat_c
    real(double), pointer, dimension(:)     :: pos_x, pos_y, pos_z
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
    real(double), pointer, dimension(:,:)   :: stress
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
      procedure, public   :: print_md_energy
      procedure, public   :: dump_stats
      procedure, public   :: dump_frame
      procedure, public   :: dump_tdep

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
  subroutine init_model(mdl, ensemble, timestep, thermo, baro)

    ! passed variables
    class(type_md_model), intent(inout)       :: mdl
    character(3), intent(in)                  :: ensemble
    real(double), intent(in)                  :: timestep
    type(type_thermostat), intent(in), target :: thermo
    type(type_barostat), intent(in), target   :: baro

    if (inode==ionode .and. iprint_MD > 2) &
      write(io_lun,'(2x,a)') "Initialising model"

    mdl%append = .false.
    if (flag_MDcontinue) mdl%append = .true.

    ! General MD arrays
    mdl%step = 0
    mdl%natoms = ni_in_cell
    mdl%ensemble = ensemble
    mdl%timestep = timestep
    mdl%species       => species
    mdl%lat_a         => rcellx
    mdl%lat_b         => rcelly
    mdl%lat_c         => rcellz
    mdl%pos_x         => x_atom_cell
    mdl%pos_y         => y_atom_cell
    mdl%pos_z         => z_atom_cell
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

    use input_module,     only: leqi

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "get_cons_qty"

    select case(mdl%ensemble)
    case("nve")
      mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy
    case("nvt")
      if (flag_extended_system) then
        mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + &
                      mdl%nhc_energy
      else
        mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy
      end if
    case("nph")
      if (flag_extended_system) then
        mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + &
                      mdl%box_kinetic_energy + mdl%PV
      else
        mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + mdl%PV
      end if
    case("npt")
      if (flag_extended_system) then
        select case(mdl%thermo_type)
        case('nhc')
          mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + &
                        mdl%nhc_energy + mdl%box_kinetic_energy + mdl%PV
        case('svr')
          mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + &
                        mdl%box_kinetic_energy + mdl%PV
        end select
      else
        mdl%h_prime = mdl%ion_kinetic_energy + mdl%dft_total_energy + mdl%PV
      end if
    end select

  end subroutine get_cons_qty
  !!***

  !!****m* md_model/print_md_energy *
  !!  NAME
  !!   print_md_energy 
  !!  PURPOSE
  !!   Print MD output at the end of each ionic step
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  SOURCE
  !!  
  subroutine print_md_energy(mdl)

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl

    ! local variables

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "print_md_energy"

    if (inode==ionode) then
      write (io_lun, '(4x,"Conserved quantity H''  : ",f15.8)') &
        mdl%h_prime
      write (io_lun, '(4x,a)') "Components of conserved quantity"
      write (io_lun, '(4x,"Kinetic energy          : ",f15.8)') &
        mdl%ion_kinetic_energy
      write (io_lun, '(4x,"Potential energy        : ",f15.8)') &
        mdl%dft_total_energy
      select case(mdl%ensemble)
      case('nvt')
        if (flag_extended_system) & 
          write (io_lun, '(4x,"Nose-Hoover energy      : ",f15.8)') &
            mdl%nhc_energy
      case('npt')
        if (flag_extended_system) then
          if (mdl%nequil < 1) then
            if (leqi(mdl%thermo_type, 'nhc')) then
              write (io_lun, '(4x,"Nose-Hoover energy      : ",f15.8)') &
                mdl%nhc_energy
            end if
            write (io_lun, '(4x,"Box kinetic energy      : ",f15.8)') &
              mdl%box_kinetic_energy
            write (io_lun, '(4x,"PV                      : ",f15.8)') &
              mdl%PV
          end if
        end if
      end select
    end if
  end subroutine print_md_energy

  !!****m* md_model/dump_stats *
  !!  NAME
  !!   dump_stats 
  !!  PURPOSE
  !!   dump thermodynamics stats to a file 
  !!  AUTHOR
  !!   Zamaan Raza 
  !!  SOURCE
  !!  
  subroutine dump_stats(mdl, filename, nequil)

    use input_module,     only: io_assign, io_close
    use md_control,       only: fac_HaBohr32GPa, md_berendsen_equil

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
    character(len=*), intent(in)          :: filename
    integer, intent(in)                   :: nequil

    ! local variables
    integer                               :: lun
    real(double)                          :: P_GPa

    if (inode==ionode .and. iprint_MD > 1) &
      write(io_lun,'(2x,"Writing statistics to ",a)') filename

    ! Convert units if necessary
    P_GPa = mdl%P_int*fac_HaBohr32GPa

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
          if (flag_extended_system .or. md_berendsen_equil > 0) then
            write(lun,'(a10,4a18,2a12)') "step", "pe", "ke", "nhc", "H'", "T", "P"
          else
            write(lun,'(a10,3a16,2a12)') "step", "pe", "ke", "H'", "T", "P"
          end if
        case ('nph')
          if (flag_extended_system .or. md_berendsen_equil > 0) then
            write(lun,'(a10,5a18,2a12,a16)') "step", "pe", "ke", "box", &
              "pV", "H'", "T", "P", "V"
          else
            write(lun,'(a10,4a18,2a12,a16)') "step", "pe", "ke", "pV", &
              "H'", "T", "P", "V"
          end if
        case ('npt')
          if (flag_extended_system .or. md_berendsen_equil > 0) then
            if (leqi(mdl%thermo_type, 'nhc')) then
              write(lun,'(a10,6a18,2a12,a16)') "step", "pe", "ke", "nhc", &
                "box", "pV", "H'", "T", "P", "V"
            else
              write(lun,'(a10,5a18,2a12,a16)') "step", "pe", "ke", "box", &
                "pV", "H'", "T", "P", "V"
            end if
          else
            write(lun,'(a10,4a18,2a12,a16)') "step", "pe", "ke", "pV", &
              "H'", "T", "P", "V"
          end if
        end select
      end if
      select case (mdl%ensemble)
      case ('nve')
        write(lun,'(i10,3e18.8,2f12.4)') mdl%step, mdl%dft_total_energy, &
          mdl%ion_kinetic_energy, mdl%h_prime, mdl%T_int, P_GPa
      case ('nvt')
        if (flag_extended_system) then
          if (nequil > 0) then
            write(lun,'(i10,4e18.8,2f12.4)') mdl%step, mdl%dft_total_energy, &
              mdl%ion_kinetic_energy, zero, mdl%h_prime, mdl%T_int, &
              P_GPa
          else
            write(lun,'(i10,4e18.8,2f12.4)') mdl%step, mdl%dft_total_energy, &
              mdl%ion_kinetic_energy, mdl%nhc_energy, mdl%h_prime, mdl%T_int, &
              P_GPa
          end if
        else
          write(lun,'(i10,3e18.8,2f12.4)') mdl%step, mdl%dft_total_energy, &
            mdl%ion_kinetic_energy, mdl%h_prime, mdl%T_int, P_GPa
        end if
      case ('nph')
        if (flag_extended_system) then
          if (nequil > 0) then
            write(lun,'(i10,5e18.8,2f12.4,e16.8)') mdl%step, &
              mdl%dft_total_energy, mdl%ion_kinetic_energy, zero, &
              mdl%PV, mdl%h_prime, mdl%T_int, P_GPa, mdl%volume
          else
            write(lun,'(i10,5e18.8,2f12.4,e16.8)') mdl%step, &
              mdl%dft_total_energy, mdl%ion_kinetic_energy, &
              mdl%box_kinetic_energy, mdl%PV, mdl%h_prime, mdl%T_int, P_GPa, &
              mdl%volume
          end if
        else
          write(lun,'(i10,4e18.8,2f12.4,e16.8)') mdl%step, &
            mdl%dft_total_energy, mdl%ion_kinetic_energy, mdl%PV, &
            mdl%h_prime, mdl%T_int, P_GPa, mdl%volume
        end if
      case ('npt')
        if (flag_extended_system) then
          if (leqi(mdl%thermo_type, 'nhc')) then
            if (nequil > 0) then
              write(lun,'(i10,6e18.8,2f12.4,e16.8)') mdl%step, &
                mdl%dft_total_energy, mdl%ion_kinetic_energy, zero, zero, &
                mdl%PV, mdl%h_prime, mdl%T_int, P_GPa, mdl%volume
            else
              write(lun,'(i10,6e18.8,2f12.4,e16.8)') mdl%step, &
                mdl%dft_total_energy, mdl%ion_kinetic_energy, mdl%nhc_energy, &
                mdl%box_kinetic_energy, mdl%PV, mdl%h_prime, mdl%T_int, &
                P_GPa, mdl%volume
            end if
          else
            write(lun,'(i10,5e18.8,2f12.4,e16.8)') mdl%step, &
              mdl%dft_total_energy, mdl%ion_kinetic_energy, &
              mdl%box_kinetic_energy, mdl%PV, mdl%h_prime, mdl%T_int, P_GPa, &
              mdl%volume
          end if
        else
          write(lun,'(i10,4e18.8,2f12.4,e16.8)') mdl%step, &
            mdl%dft_total_energy, mdl%ion_kinetic_energy, mdl%PV, &
            mdl%h_prime, mdl%T_int, P_GPa, mdl%volume
        end if
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

    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
    character(len=*), intent(in)          :: filename

    ! local variables
    integer                               :: lun, i

    if (inode==ionode) then
      if (iprint_MD > 1) write(io_lun,'(2x,"Writing frame to ",a)') filename
      call io_assign(lun)
      if (mdl%append) then
        open(unit=lun,file=filename,position='append')
      else 
        open(unit=lun,file=filename,status='replace')
      end if

      write(lun,'("frame ",i8)') mdl%step
      write(lun,'(a)') "cell_vectors"
      do i=1,3
        write(lun,'(3f12.6)') mdl%lattice_vec(i,:)
      end do
      write(lun,'(a)') "end cell_vectors"
      write(lun,'(a)') "stress_tensor"
      write(lun,'(3f14.6)') mdl%stress(1,:)
      write(lun,'(3f14.6)') mdl%stress(2,:)
      write(lun,'(3f14.6)') mdl%stress(3,:)
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

  !!****m* md_model/dump_tdep *
  !!  NAME
  !!   dump_tdep_frame
  !!  PURPOSE
  !!   Dump a MD step for TDEP postprocessing
  !!  AUTHOR
  !!   Zamaan Raza
  !!  SOURCE
  !!
  subroutine dump_tdep(mdl)
 
    use input_module,     only: io_assign, io_close
    use units,            only: HaToEv, BohrToAng
    use md_control,       only: fac_HaBohr32GPa
 
    ! passed variables
    class(type_md_model), intent(inout)   :: mdl
 
    ! local variables
    integer                               :: lun1, lun2, i
 
    if (inode==ionode) then
      if (flag_MDdebug .and. iprint_MD > 1) &
        write(io_lun,*) "Writing TDEP output"
      call io_assign(lun1)
      open(unit=lun1,file=file_meta,status='replace')
      write(lun1,'(i12,a)') mdl%natoms, " # N atoms"
      write(lun1,'(i12,a)') mdl%step, " # N time steps"
      write(lun1,'(f12.2,a)') mdl%timestep, " # time step in fs"
      write(lun1,'(f12.2,a)') mdl%T_ext, " # temperature in K"
      call io_close(lun1)

      call io_assign(lun1)
      call io_assign(lun2)
 
      if (mdl%append) then
        open(unit=lun1,file=file_positions,status='old',position='append')
        open(unit=lun2,file=file_forces,status='old',position='append')
      else 
        open(unit=lun1,file=file_positions,status='replace')
        open(unit=lun2,file=file_forces,status='replace')
      end if
 
      do i=1,mdl%natoms
        write(lun1,'(3e20.12)') mdl%pos_x(i)/mdl%lat_a, &
                                mdl%pos_y(i)/mdl%lat_b, &
                                mdl%pos_z(i)/mdl%lat_c
        write(lun2,'(3e20.12)') mdl%atom_force(:,i)*HaToeV/BohrToAng
      end do
      call io_close(lun1)
      call io_close(lun2)

      call io_assign(lun1)
      if (mdl%append) then
        open(unit=lun1,file=file_stat,status='old',position='append')
      else 
        open(unit=lun1,file=file_stat,status='replace')
      end if

      write(lun1,'(i7,f10.1,3e16.8,8e14.6)') mdl%step, mdl%step*mdl%timestep, &
        (mdl%dft_total_energy+mdl%ion_kinetic_energy)*HaToeV, &
        mdl%dft_total_energy*HaToeV, &
        mdl%ion_kinetic_energy*HaToeV, mdl%T_int, mdl%P_int, &
        mdl%stress(1,1)*fac_HaBohr32GPa, mdl%stress(2,2)*fac_HaBohr32GPa, &
        mdl%stress(3,3)*fac_HaBohr32GPa, zero, zero, zero
 
      call io_close(lun1)
    end if

  end subroutine dump_tdep
  !!***

end module md_model
