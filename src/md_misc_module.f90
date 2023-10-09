! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module md_misc
! ------------------------------------------------------------------------------
! Code area 7: move atoms
! ------------------------------------------------------------------------------

!!****h* Conquest/md_misc
!!  NAME
!!   md_misc
!!  PURPOSE
!!   Miscellaneous MD routines, here in part to avoid circular dependence
!!  USES
!!   
!!  AUTHOR
!!   David Bowler
!!  CREATION DATE
!!   2022/09/30 08:43
!!  MODIFICATION HISTORY
!!  SOURCE
!!
module md_misc

  use datatypes
  use numbers
  use global_module,    only: ni_in_cell, io_lun, iprint_MD, &
                              temp_ion, flag_MDcontinue, flag_MDdebug
  use move_atoms,       only: fac_Kelvin2Hartree
  use species_module,   only: species, mass
  use GenComms,         only: inode, ionode, cq_abort
  use input_module,     only: leqi
  use rng,              only: type_rng
  use units,            only: HaBohr3ToGPa
  use md_control,       only: ion_velocity, md_nhc_mass, md_nhc_cell_mass, &
       flag_write_extxyz, flag_write_xsf, md_baro_type, md_ensemble, MDtimestep, &
       flag_nhc, md_tau_P, md_tau_P_equil, md_tau_T, md_tau_T_equil, md_thermo_type

  implicit none

  character(20) :: md_thermo_file = "md.thermostat"
  character(20) :: md_baro_file = "md.barostat"
  character(20) :: md_trajectory_file = "trajectory.xsf"
  character(20) :: md_frames_file = "md.frames"
  character(20) :: md_stats_file = "md.stats"
  character(20) :: md_heat_flux_file = "md.heatflux"

contains

  !!****m* md_misc/init_md *
  !!  NAME
  !!   integrate_ensemble
  !!  PURPOSE
  !!   Initialise all ensemble-related MD variables
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/04/32  16:52
  !!  MODIFICATION HISTORY
  !!   2018/05/30 zamaan
  !!    Made BerendsenEquil work with NVT ensemble.
  !!   2019/04/09 zamaan
  !!    Removed unnecessary references to stress tensor
  !!   2019/05/22 14:44 dave & tsuyoshi
  !!    Tweak for initial velocities
  !!   2020/10/07 tsuyoshi
  !!    allocation of atom_vels has been moved from read_atomic_positions
  !!    (introduced atom_vels in July, 2020.)
  !!   2022/10/04 17:27 dave
  !!    Reworking to set initial KE of ions correctly
  !!  
  subroutine init_md(baro, thermo, mdl, md_ndof, nequil, second_call)

    use numbers
    use input_module,   only: leqi
    use io_module,      only: read_velocity
    use md_model,       only: type_md_model
    use md_control,     only: lattice_vec, type_thermostat, type_barostat, read_md_checkpoint
    use GenComms,       only: inode, ionode, gcopy, cq_warn
    use memory_module,  only: reg_alloc_mem, type_dbl
    use global_module,  only: rcellx, rcelly, rcellz, temp_ion, ni_in_cell, &
         flag_MDcontinue, flag_read_velocity, &
         flag_MDdebug, iprint_MD, flag_atomic_stress, &
         atomic_stress, area_moveatoms, &
         id_glob, atom_vels, &
         flag_full_stress, flag_heat_flux, min_layer
    use move_atoms,     only: init_velocity
    use io_module,      only: return_prefix

    ! passed variables
    type(type_barostat), intent(inout)    :: baro
    type(type_thermostat), intent(inout)  :: thermo
    type(type_md_model), intent(inout)    :: mdl
    integer, intent(in)                   :: md_ndof
    integer, intent(in)                   :: nequil
    logical, optional                     :: second_call

    ! local variables
    character(50)  :: file_velocity='velocity.dat'
    integer       :: stat, i
    character(len=12) :: subname = "init_md: "
    character(len=120) :: prefix
    real(double) :: ion_KE

    prefix = return_prefix(subname, min_layer)
    if (inode==ionode .and. iprint_MD + min_layer > 1) &
         write(io_lun,'(4x,a)') trim(prefix)//" MD ensemble is "//md_ensemble

    if (.not. allocated(ion_velocity)) then
       allocate(ion_velocity(3,ni_in_cell), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating velocity in init_md: ", &
            ni_in_cell, stat)
       call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)
    end if

    !  2020/10/7 Tsuyoshi Miyazaki
    !   Order of atoms is different between ion_velocity and atom_vels:
    !   we are planning to remove ion_velocity in the future.
    !       global_labelling :  atom_coord & atom_vels
    !     partition_labelling:  x(or y or z)_atom_cell & ion_velocity
    !
    if (.not. allocated(atom_vels)) then
       allocate(atom_vels(3,ni_in_cell), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating atom_vels in init_md: ", &
            ni_in_cell, stat)
       call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)
       atom_vels = zero
    end if

    if (.not. flag_MDcontinue) then
       ! I've moved the velocity initialisation here to make reading the new
       ! unified md checkpoint file easier - zamaan
       ion_velocity = zero
       if (flag_read_velocity) then
          call read_velocity(ion_velocity, file_velocity)
       else
          if(temp_ion > RD_ERR) then
             call init_velocity(ni_in_cell, temp_ion, ion_velocity, ion_ke)
          end if
       end if

       ! atom_vels : 2020/Jul/30
       do i=1,ni_in_cell
          atom_vels(1,id_glob(i)) = ion_velocity(1,i)
          atom_vels(2,id_glob(i)) = ion_velocity(2,i)
          atom_vels(3,id_glob(i)) = ion_velocity(3,i)
       end do

    end if

    if(leqi(md_thermo_type,'svr').AND.md_tau_T<zero) then
       call cq_warn("init_md","No value set for SVR time constant; defaulting to 50fs")
       md_tau_T = 50.0_double
    end if

    select case(md_ensemble)
    case('nve')
       ! Just for computing temperature
       call thermo%init_thermo('none', 'none', MDtimestep, md_ndof, md_tau_T, &
            ion_ke)
       call baro%init_baro('none', MDtimestep, md_ndof, ion_velocity, &
            md_tau_P, ion_ke)
    case('nvt')
       if (nequil > 0) then ! Equilibrate ?
          if (inode==ionode .and. iprint_MD + min_layer > 1) then
             write (io_lun, '(4x,a,i8," steps")') &
                  trim(prefix)//"Equilibrating using SVR thermostat for ",nequil
          end if
          call thermo%init_thermo('svr', 'none', MDtimestep, md_ndof, &
               md_tau_T, ion_ke)
       else
          call thermo%init_thermo(md_thermo_type, 'none', MDtimestep, md_ndof, &
               md_tau_T, ion_ke)
       end if
       call baro%init_baro('none', MDtimestep, md_ndof, ion_velocity, &
            md_tau_P, ion_ke)
    case('nph')
       if (nequil > 0) then ! Equilibrate ?
          if (inode==ionode .and. iprint_MD + min_layer > 1) then
             write (io_lun, '(4x,a,i8," steps")') &
                  trim(prefix)//"Equilibrating using PR barostat for ",nequil
          end if
          call thermo%init_thermo('none', 'pr', MDtimestep, &
               md_ndof, md_tau_T_equil, &
               ion_ke)
          call baro%init_baro(md_baro_type, MDtimestep, md_ndof, &
               ion_velocity, md_tau_P_equil, &
               ion_ke)
       else
          call thermo%init_thermo('none', md_baro_type, MDtimestep, &
               md_ndof, md_tau_T, ion_ke)
          call baro%init_baro(md_baro_type, MDtimestep, md_ndof, &
               ion_velocity, md_tau_P, ion_ke)
       end if
    case('npt')
       if (nequil > 0) then ! Equilibrate ?
          if (inode==ionode .and. iprint_MD + min_layer > 1) then
             write (io_lun, '(4x,a,i8," steps")') &
                  trim(prefix)//"Equilibrating using SVR baro/thermostat for ",nequil
          end if
          call thermo%init_thermo('svr', 'pr', MDtimestep, &
               md_ndof, md_tau_T_equil, &
               ion_ke)
          call baro%init_baro('pr', MDtimestep, md_ndof, &
               ion_velocity, md_tau_P_equil, &
               ion_ke)
       else
          call thermo%init_thermo(md_thermo_type, md_baro_type, MDtimestep, &
               md_ndof, md_tau_T, ion_ke)
          call baro%init_baro(md_baro_type, MDtimestep, md_ndof, &
               ion_velocity, md_tau_P, ion_ke)
       end if
    end select

    ! NB Moved from previous position before thermo/baro initialisation
    if (.not. present(second_call)) then
       ! Initialise the model only once per run
       lattice_vec = zero
       lattice_vec(1,1) = rcellx
       lattice_vec(2,2) = rcelly
       lattice_vec(3,3) = rcellz
       call mdl%init_model(md_ensemble, MDtimestep, thermo, baro)
    end if
    mdl%ion_kinetic_energy = ion_ke

    ! N.B. atomic stress is allocated in initialisation_module/initialise! - zamaan

    if (flag_MDcontinue) call read_md_checkpoint(thermo, baro)

  end subroutine init_md
  !!***

  !!****m* md_md_misc/end_md *
  !!  NAME
  !!   end_md
  !!  PURPOSE
  !!   Deallocate MD arrays that are not deallocated elsewhere
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/05/08
  !!  MODIFICATION HISTORY
  !!   2020/10/07 tsuyoshi
  !!    added deallocation of atom_vels
  !!  SOURCE
  !!  
  subroutine end_md(th, baro)

    use GenComms,         only: inode, ionode
    use global_module,    only: area_moveatoms, atomic_stress, &
         flag_atomic_stress, flag_MDdebug, &
         iprint_MD, ni_in_cell, area_moveatoms, atom_vels
    use memory_module,    only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use md_control,     only: type_thermostat, type_barostat

    ! passed variables
    class(type_thermostat), intent(inout)   :: th
    class(type_barostat), intent(inout)     :: baro

    ! local variables
    integer :: stat

    if (inode==ionode .and. flag_MDdebug .and. iprint_MD > 2) &
         write(io_lun,'(2x,a)') "end_md"

    if (flag_nhc) then
       deallocate(md_nhc_mass, md_nhc_cell_mass) ! allocated in initial_read
       deallocate(th%eta, th%v_eta, th%G_nhc, th%m_nhc)
       deallocate(th%eta_cell, th%v_eta_cell, th%G_nhc_cell, th%m_nhc_cell)
       call reg_dealloc_mem(area_moveatoms, 8*th%n_nhc, type_dbl)
    end if

    deallocate(atom_vels, STAT=stat)
    if (stat /= 0) call cq_abort("Error deallocating atom_vels in end_md: ", &
         ni_in_cell, stat)
    call reg_dealloc_mem(area_moveatoms, 3 * ni_in_cell, type_dbl)
    deallocate(ion_velocity, STAT=stat)
    if (stat /= 0) call cq_abort("Error deallocating velocity in end_md: ", &
         ni_in_cell, stat)
    call reg_dealloc_mem(area_moveatoms, 3 * ni_in_cell, type_dbl)
    if (flag_atomic_stress) then
       deallocate(atomic_stress, STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error deallocating atomic_stress in end_md: ", &
            ni_in_cell, stat)
       call reg_dealloc_mem(area_moveatoms, 3*3*ni_in_cell, type_dbl)
    end if

  end subroutine end_md

  !!****m* md_misc/integrate_pt *
  !!  NAME
  !!   integrate_pt
  !!  PURPOSE
  !!   Integrate the thermostat and barostat depending on the ensemble
  !!   and thermostat/barostat type
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/04/23 11:49
  !!  SOURCE
  !!  
  subroutine integrate_pt(baro, thermo, mdl, velocity, second_call)

    use numbers
    use global_module,    only: iprint_MD, min_layer
    use GenComms,         only: inode, ionode
    use md_model,         only: type_md_model
    use md_control,     only: type_thermostat, type_barostat
    use io_module,      only: return_prefix

    ! passed variables
    type(type_barostat), intent(inout)          :: baro
    type(type_thermostat), intent(inout)        :: thermo
    type(type_md_model), intent(inout)          :: mdl
    real(double), dimension(:,:), intent(inout) :: velocity
    logical, optional                           :: second_call

    ! local variables
    character(len=16) :: subname = "integrate_pt: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    if (inode==ionode .and. iprint_MD + min_layer > 1) then
       if (present(second_call)) then
          write(io_lun,'(4x,a)') trim(prefix)//" second call"
       else
          write(io_lun,'(4x,a)') trim(prefix)//" first call"
       end if
    end if

    select case(md_ensemble)
    case('nvt')
       select case(thermo%thermo_type)
       case('nhc')
          call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
          if (present(second_call)) call thermo%get_thermostat_energy
       case('svr')
          if (present(second_call)) then
             call thermo%v_rescale(velocity)
             call thermo%get_thermostat_energy
          else
             call thermo%get_svr_thermo_sf(MDtimestep)
          end if
       end select
    case('nph')
       select case(baro%baro_type)
       case('pr')
          if (present(second_call)) then
             call baro%couple_box_particle_velocity(thermo, velocity)
             call thermo%get_temperature_and_ke(baro, velocity, &
                  mdl%ion_kinetic_energy)
             call baro%get_pressure_and_stress
             call baro%integrate_box(thermo)
          else
             call thermo%get_temperature_and_ke(baro, velocity, &
                  mdl%ion_kinetic_energy)
             call baro%get_pressure_and_stress
             call baro%integrate_box(thermo)
             call baro%couple_box_particle_velocity(thermo, velocity)
          end if
       end select
    case('npt')
       select case(baro%baro_type)
       case('mttk')
          call baro%propagate_npt_mttk(thermo, mdl%ion_kinetic_energy, &
               velocity)
       case('pr')
          if (present(second_call)) then
             call baro%couple_box_particle_velocity(thermo, velocity)
             call thermo%get_temperature_and_ke(baro, velocity, &
                  mdl%ion_kinetic_energy)
             call baro%get_pressure_and_stress
             call baro%integrate_box(thermo)
             select case(thermo%thermo_type)
             case('nhc')
                call baro%couple_box_particle_velocity(thermo, velocity)
                call thermo%get_temperature_and_ke(baro, velocity, &
                     mdl%ion_kinetic_energy)
                call baro%get_pressure_and_stress
                call baro%integrate_box(thermo)
                call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
                call thermo%get_thermostat_energy
             case('svr')
                call thermo%get_temperature_and_ke(baro, velocity, &
                     mdl%ion_kinetic_energy)
                call baro%get_pressure_and_stress
                call baro%couple_box_particle_velocity(thermo, velocity)
                call baro%integrate_box(thermo)
                call thermo%get_svr_thermo_sf(MDtimestep/two, baro)
                call thermo%v_rescale(velocity)
                call thermo%get_thermostat_energy
                call baro%scale_box_velocity(thermo)
             end select
          else
             select case(thermo%thermo_type)
             case('nhc')
                call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
                call thermo%get_temperature_and_ke(baro, velocity, &
                     mdl%ion_kinetic_energy)
                call baro%get_pressure_and_stress
                call baro%integrate_box(thermo)
                call baro%couple_box_particle_velocity(thermo, velocity)
             case('svr')
                call thermo%get_svr_thermo_sf(MDtimestep/two, baro)
                call thermo%v_rescale(velocity)
                call baro%scale_box_velocity(thermo)
                call thermo%get_temperature_and_ke(baro, velocity, &
                     mdl%ion_kinetic_energy)
                call baro%get_pressure_and_stress
                call baro%integrate_box(thermo)
                call baro%couple_box_particle_velocity(thermo, velocity)
             end select
             call thermo%get_temperature_and_ke(baro, velocity, &
                  mdl%ion_kinetic_energy)
             call baro%get_pressure_and_stress
             call baro%integrate_box(thermo)
             call baro%couple_box_particle_velocity(thermo, velocity)
          end if
       end select
    end select
  end subroutine integrate_pt
  !!*****

  !!****m* md_misc/update_pos_and_box *
  !!  NAME
  !!   update_pos_and_box
  !!  PURPOSE
  !!   Perform the (modified) velocity Verlet position and box updates, &
  !!   depending on the ensemble and integrator
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/08/11 10:27
  !!  SOURCE
  !!  
  subroutine update_pos_and_box(baro, nequil, flag_movable)

    use Integrators,   only: vVerlet_r_dt
    use io_module,     only: leqi
    use GenComms,      only: inode, ionode
    use global_module, only: iprint_MD, min_layer
    use io_module,      only: return_prefix
    use md_control,     only: type_barostat

    ! passed variables
    type(type_barostat), intent(inout)  :: baro
    integer, intent(in)                 :: nequil
    logical, dimension(:), intent(in)   :: flag_movable

    ! local variables
    character(len=20) :: subname = "update_pos_and_box: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    if (inode==ionode .and. iprint_MD + min_layer > 1) &
         write(io_lun,'(4x,a)') trim(prefix)//" starting"

    if (leqi(md_ensemble, 'npt')) then
       if (nequil > 0) then
          call vVerlet_r_dt(MDtimestep,ion_velocity,flag_movable)
       else
          ! Modified velocity Verlet position step for NPT ensemble
          select case(md_baro_type)
          case('pr')
             call baro%propagate_box_ssm
             call vVerlet_r_dt(MDtimestep,ion_velocity,flag_movable)
             call baro%propagate_box_ssm
          case('mttk')
             call baro%propagate_r_mttk(MDtimestep, ion_velocity, flag_movable)
             call baro%propagate_box_mttk(MDtimestep)
          end select
       end if
    else
       ! For NVE or NVT
       call vVerlet_r_dt(MDtimestep,ion_velocity,flag_movable)
    end if

  end subroutine update_pos_and_box
  !!*****

  !!****m* md_misc/get_heat_flux *
  !!  NAME
  !!   get_heat_flux
  !!  PURPOSE
  !!   Compute the heat flux according to Green-Kubo formalism as in Carbogno 
  !!   et al. PRL 118, 175901 (2017). Note that in this implementation, we
  !!   ignore the *convective* contribution to the heat flux, so it is only
  !!   applicable to solids! 
  !!      J = sum_i(sigma_i . v_i) 
  !!        sigma_i = stress contribution from atom i
  !!        v_i = velocity of atom i).
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/05/08
  !!  SOURCE
  !!  
  subroutine get_heat_flux(atomic_stress, velocity, heat_flux)

    use numbers
    use GenComms,      only: inode, ionode
    use global_module, only: iprint_MD, ni_in_cell, flag_heat_flux, min_layer
    use io_module,     only: return_prefix

    ! Passed variables
    real(double), dimension(3,3,ni_in_cell), intent(in)   :: atomic_stress
    real(double), dimension(3,ni_in_cell), intent(in)     :: velocity
    real(double), dimension(3), intent(out)               :: heat_flux

    ! local variables
    integer :: i
    character(len=16) :: subname = "get_heat_flux: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    if (inode==ionode .and. iprint_MD + min_layer > 2) &
         write(io_lun,'(6x,a)') "Welcome to get_heat_flux"

    heat_flux = zero
    do i=1,ni_in_cell
       heat_flux = heat_flux + matmul(atomic_stress(:,:,i), velocity(:,i))
    end do

  end subroutine get_heat_flux
  !!*****

  !!****m* md_misc/write_md_data *
  !!  NAME
  !!   write_md_data
  !!  PURPOSE
  !!   Write MD data to various files at the end of an ionic step
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/08/11 10:27
  !!  MODIFIED
  !!   2021/10/19 Jianbo Lin
  !!    Added call for extended XYZ output (includes forces)
  !!  SOURCE
  !!  
  subroutine write_md_data(iter, thermo, baro, mdl, nequil, MDfreq)

    use GenComms,      only: inode, ionode
    use io_module,     only: write_xsf, write_extxyz, return_prefix
    use global_module, only: iprint_MD, flag_baroDebug, flag_thermoDebug, &
         flag_heat_flux, min_layer
    use md_model,      only: type_md_model, md_tdep
    use md_control,     only: type_thermostat, type_barostat, write_md_checkpoint

    ! Passed variables
    type(type_barostat), intent(inout)    :: baro
    type(type_thermostat), intent(inout)  :: thermo
    type(type_md_model), intent(inout)    :: mdl
    integer, intent(in)                   :: iter, nequil, MDfreq

    ! local variables
    character(len=16) :: subname = "write_md_data: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    if (inode==ionode .and. iprint_MD + min_layer > 1) &
         write(io_lun,'(4x,a)') trim(prefix)//" starting"

    call write_md_checkpoint(thermo, baro)
    call mdl%dump_stats(md_stats_file, nequil)
    if (flag_write_xsf) call write_xsf(md_trajectory_file, iter)
    if (flag_write_extxyz) &
         call write_extxyz('trajectory.xyz', mdl%dft_total_energy, mdl%atom_force)
    if (flag_heat_flux) call mdl%dump_heat_flux(md_heat_flux_file)
    if (mod(iter, MDfreq) == 0) then
       call mdl%dump_frame(md_frames_file)
       if (md_tdep) call mdl%dump_tdep
    end if
    if (flag_thermoDebug) &
         call thermo%dump_thermo_state(iter, md_thermo_file)
    if (flag_baroDebug) &
         call baro%dump_baro_state(iter, md_baro_file)
    mdl%append = .true.

  end subroutine write_md_data
  !!*****

end module md_misc
