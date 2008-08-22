! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module minimise
! ------------------------------------------------------------------------------
! Code area 6: Energy minimisation
! ------------------------------------------------------------------------------

!!****h* Conquest/minimise *
!!  NAME
!!   minimise
!!  PURPOSE
!!   Controls energy minimisation and force finding
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07:54, 2003/02/04 dave
!!  MODIFICATION HISTORY
!!   Added call to get_H_matrix to build new H and charge in get_E_and_F
!!   15:02, 12/03/2003 drb 
!!    Added call to get_energy after FindMinDM
!!   15:42, 2003/06/09 dave
!!    Changed to use new blip minimisation scheme
!!   2006/09/25 17:15 dave
!!    Added option to use pulay for SF minimisation
!!   2008/05/25 ast
!!    Added timers
!!  SOURCE
!!
module minimise

  use datatypes
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_eminimisation

  implicit none

  save

  logical :: UsePulay   

  integer :: n_L_iterations, n_support_iterations
  real(double) :: L_tolerance, sc_tolerance, energy_tolerance, expected_reduction

  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id$"
!!***

contains

!!****f* minimise/get_E_and_F *
!!
!!  NAME 
!!   get_E_and_F
!!  USAGE
!! 
!!  PURPOSE
!!   Finds the ground-state energy and force
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   1998 sometime I think
!!  MODIFICATION HISTORY
!!   24/05/2001 dave
!!    ROBODoc header, indentation and stripping calls
!!   13/06/2001 dave
!!    Adapted to use force_module and added RCS Id and Log tags
!!   07:43, 2003/02/04 dave
!!    Changed call to minimise to pick and choose between full run (varying DM, self-consistency and blips), 
!!    self-consistent run (varying DM, self-consistency) and simple TB (varying DM only)
!!   14:39, 26/02/2003 drb 
!!    Added call to get_H_matrix to build new H and charge
!!   08:35, 2003/03/12 dave
!!    Added call to get_energy after FindMinDM
!!   08:30, 2003/10/01 dave
!!    Changed flag_vary_blips to basis and added basis switch for outer loop minimiser
!!   2006/11/13 18:20 dave
!!    Added new flags to test for finding and/or writing forces
!!   2007/04/17 09:34 dave
!!    Passed no. of L iterations to vary_support
!!   2008/05/25 ast
!!    Added timer
!!  SOURCE
!!
  subroutine get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, total_energy, find_forces, write_forces)

    use datatypes
    use force_module, ONLY : force
    use DMMin, ONLY: FindMinDM
    use SelfCon, ONLY: new_SC_potl, atomch_output, get_atomic_charge
    use global_module, ONLY: flag_vary_basis, flag_self_consistent, flag_basis_set, blips, PAOs, &
                             IPRINT_TIME_THRES1
    use energy, ONLY: get_energy
    use GenComms, ONLY: cq_abort, inode, ionode
    use blip_minimisation, ONLY: vary_support
    use pao_minimisation, ONLY: vary_pao, pulay_min_pao
    use timer_module

    implicit none

    ! Shared variables
    logical :: vary_mu, fixed_potential, find_forces, write_forces

    integer :: n_save_freq, n_run

    character(len=40) :: output_file

    real(double) :: number_of_bands, mu
    real(double) :: total_energy

    ! Local variables
    logical :: reset_L
    real(double) :: electrons
    type(cq_timer) :: tmr_l_energy, tmr_l_force

    call start_timer(tmr_std_eminimisation)
    ! reset_L = .true.  ! changed by TM, Aug 2008 
    reset_L = .false.
    call start_timer(tmr_l_energy,WITH_LEVEL)      ! Start timing the energy calculation
    ! Now choose what we vary
    if(flag_vary_basis) then ! Vary everything: DM, charge density, basis set
       if(flag_basis_set==blips) then
          call vary_support( n_support_iterations, fixed_potential, vary_mu, n_L_iterations, &
               number_of_bands, L_tolerance, sc_tolerance, energy_tolerance, mu, total_energy, &
               expected_reduction)
       else if(flag_basis_set==PAOs) then
          if(UsePulay) then
             call pulay_min_pao( n_support_iterations, fixed_potential, vary_mu, n_L_iterations, &
                  number_of_bands, L_tolerance, sc_tolerance, energy_tolerance, mu, total_energy, &
                  expected_reduction)
          else
             call vary_pao( n_support_iterations, fixed_potential, vary_mu, n_L_iterations, &
                  number_of_bands, L_tolerance, sc_tolerance, energy_tolerance, mu, total_energy, &
                  expected_reduction)
          end if
       else 
          call cq_abort("get_E_and_F: basis set undefined: ",flag_basis_set)
       end if
    else if(flag_self_consistent) then ! Vary only DM and charge density
       call new_SC_potl( .false., sc_tolerance, reset_L, fixed_potential, vary_mu, n_L_iterations, &
            number_of_bands, L_tolerance, mu, total_energy)
    else ! Ab initio TB: vary only DM
       call FindMinDM(n_L_iterations, number_of_bands, vary_mu, L_tolerance, mu, inode, ionode, reset_L, .false.)
       call get_energy(total_energy)
    end if
    call stop_print_timer(tmr_l_energy, "calculating ENERGY", IPRINT_TIME_THRES1)
    if(atomch_output) call get_atomic_charge()
    if(find_forces) then 
      call start_timer(tmr_l_force,WITH_LEVEL)    ! Start timing the force calculation
      call force(fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tolerance, &
                 sc_tolerance, mu, total_energy, expected_reduction, write_forces)
      call stop_print_timer(tmr_l_force, "calculating FORCE", IPRINT_TIME_THRES1)   ! Stop timing the force calculation
    endif
    !  Print results of local timers
    call stop_timer(tmr_std_eminimisation)
    return
  end subroutine get_E_and_F
!!***

end module minimise
