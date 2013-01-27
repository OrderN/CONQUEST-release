! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module test_force_module
! ------------------------------------------------------------------------------
! Code area 7: moving atoms
! ------------------------------------------------------------------------------

!!****h* Conquest/test_force_module *
!!  NAME
!!   test_force_module
!!  PURPOSE
!!   Puts together routines to test forces in Conquest
!!
!!   While at any point we can test the total force rather easily
!!   (calculate energy and force at point x and x+delta and then test
!!   if 0.5*(f(x)+f(x+delta)) = -(E(x+delta)-E(x))/delta to some
!!   number of decimal places) testing the individual components of
!!   the force requires more care.  This is (will be) discussed on the
!!   Conquest web server in detail and in notes in the Conquest CQDocs
!!   module.
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   11:40, 02/09/2003 
!!  MODIFICATION HISTORY
!!   12:03, 30/09/2003 drb 
!!    Major re-working: fixed S-pulay; debugged Hellmann-Feynman
!!    forces for non-SC case; amalgamated local and non-local HF
!!    forces; created a full Pulay force testing routine; made
!!    possible to test individual forces easily
!!   14:08, 2003/12/19 dave
!!    Improved to allow PAO testing
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote
!!    in terms of new
!!    matrix routines
!!   2008/02/06 11:03 dave
!!    Changed for output to file not stdout
!!   2008/05/25 ast
!!    Added timers
!!   2011/10/06 13:55 dave
!!    Changes for cDFT and forces
!!  SOURCE
!!
module test_force_module

  use datatypes
  use global_module,          only: io_lun, area_moveatoms
  use timer_stdclocks_module, only: start_timer, stop_timer, &
                                    tmr_std_moveatoms,       &
                                    tmr_std_allocation

  implicit none

  logical      :: flag_test_all_forces
  integer      :: TF_direction, TF_atom_moved, flag_which_force
  real(double) :: TF_delta

!!***

contains

  !!****f* test_force_module/test_forces *
  !!
  !!  NAME 
  !!   test_forces
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the forces - overall controlling routine
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   07:49, 2003/09/03
  !!  MODIFICATION HISTORY
  !!   2008/05/25 ast
  !!    Added timers
  !!   2011/09/19 L.Tong
  !!    Added Spin polarisation
  !!   2011/10/06 13:55 dave
  !!    Added call to test cDFT forces
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_forces(fixed_potential, vary_mu, n_L_iterations, &
                         L_tolerance, tolerance, total_energy,     &
                         expected_reduction)

    use datatypes
    use move_atoms,             only: primary_update, cover_update,    &
                                      update_atom_coord, update_H
    use group_module,           only: parts
    use cover_module,           only: BCS_parts, DCS_parts
    use primary_module,         only: bundle
    use global_module,          only: iprint_MD, x_atom_cell,          &
                                      y_atom_cell, z_atom_cell,        &
                                      id_glob_inv, flag_perform_cDFT,  &
                                      flag_self_consistent,            &
                                      ni_in_cell,                      &
                                      nspin, spin_factor
    use GenComms,               only: myid, inode, ionode, cq_abort
    use io_module,              only: dump_charge, grab_charge,        &
                                      dump_matrix, grab_matrix,        &
                                      dump_blips, grab_blips,          &
                                      grab_locps, dump_locps,          &
                                      grab_projs, dump_projs
    use mult_module,            only: matH, matK, matNL, matKE
    use energy,                 only: hartree_energy, local_ps_energy, &
                                      xc_energy, nl_energy,            &
                                      kinetic_energy, band_energy,     &
                                      delta_E_hartree, delta_E_xc,     &
                                      get_energy
    use ewald_module,           only: ewald_energy
    use pseudopotential_common, only: core_correction, pseudopotential
    use H_matrix_module,        only: get_H_matrix
    use density_module,         only: set_density, density
    use maxima_module,          only: maxngrid
    use functions_on_grid,      only: supportfns
    use SelfCon,                only: new_SC_potl
    use DMMin,                  only: FindMinDM
    use force_module,           only: force, tot_force
    use memory_module,          only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L

    integer :: n_L_iterations

    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer      :: stat, spin
    logical      :: reset_L = .false.
    real(double) :: E0, F0, E1, F1, analytic_force, numerical_force, &
                    electrons_tot
    real(double) :: H0, L0, XC0, NL0, K0, B0, dH0, dXC0, EW0, T0
    real(double), dimension(nspin) :: electrons
    real(double), dimension(:,:), allocatable :: HF_force

    call start_timer(tmr_std_moveatoms)

    allocate(HF_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) call cq_abort("test_forces: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    if (inode == ionode) &
         write (io_lun,fmt='(2x,"******************"/,&
                            &2x,"* Testing Forces *"/,&
                            &2x,"******************"/)')
    if (TF_atom_moved > ni_in_cell) then
       if (inode == ionode) &
            write (io_lun,fmt='(2x,"Error: specified atom &
                               &out-of-range: ",2i8)') &
                  TF_atom_moved, ni_in_cell
       TF_atom_moved = 1
    end if
    ! Here we want to back up the state of the system by saving H, K,
    ! density, blips, chis, local ps Energies
    H0 = hartree_energy
    L0 = local_ps_energy
    XC0 = xc_energy
    NL0 = nl_energy
    K0 = kinetic_energy
    B0 = band_energy
    dH0 = delta_E_hartree 
    dXC0 = delta_E_xc
    EW0 = ewald_energy
    T0 = band_energy + delta_E_hartree + delta_E_xc + ewald_energy + &
         core_correction
    !call dump_blips("O",supportfns, inode)
    !call dump_locps(pseudopotential,n_my_grid_points,inode)
    !call dump_projs(pseudofunctions,NCF*CORE_SIZE,inode)
    !call dump_charge(density,n_my_grid_points,inode)
    !call dump_matrix("H", matH,inode)
    !call dump_matrix("KE",matKE,inode)
    !call dump_matrix("NL",matNL,inode)
    !call dump_matrix("K", matK,inode)
    ! Full forces
    if (flag_test_all_forces .or. flag_which_force == 1) then
       if (inode == ionode) write (io_lun,*) '*** Full ***'

       call test_full(fixed_potential, vary_mu, n_L_iterations, &
                      L_tolerance, tolerance, total_energy,     &
                      expected_reduction)

       !call grab_blips("O",supportfns,inode)
       !call grab_locps(pseudopotential,n_my_grid_points,inode)
       !call grab_projs(pseudofunctions,NCF*CORE_SIZE,inode)
       !call grab_charge(density,n_my_grid_points,inode)
       !call grab_matrix("H", matH,inode)
       !call grab_matrix("KE",matKE,inode)
       !call grab_matrix("NL",matNL,inode)
       !call grab_matrix("K", matK,inode)
       ! call get_H_matrix(.false., fixed_potential, electrons, density,&
       !      n_my_grid_points)

       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if (flag_self_consistent) then ! Vary only DM and charge density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          ! call force( fixed_potential, vary_mu, n_L_iterations, &
          !      number_of_bands, L_tolerance, tolerance, mu, &
          !      total_energy, expected_reduction,.true.)
       end if
       ! call get_H_matrix(.true., fixed_potential, electrons, density, &
       !      n_my_grid_points)
       !delta_E_hartree = dH0  
       !delta_E_xc      = dXC0 
       !hartree_energy  = H0  
       !local_ps_energy = L0  
       !xc_energy       = XC0 
       !nl_energy       = NL0 
       !kinetic_energy  = K0  
       !band_energy     = B0
       !ewald_energy    = EW0
       !total_energy    = T0
    end if
    ! *** Hellmann-Feynman ***
    if (flag_test_all_forces .or. flag_which_force == 2) then
       if (inode == ionode) write (io_lun,*) '*** Full HF ***'
       call test_HF(fixed_potential, vary_mu, n_L_iterations, &
                    L_tolerance, tolerance, total_energy,     &
                    expected_reduction)
       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if (flag_self_consistent) then ! Vary only DM and charge density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          ! call force( fixed_potential, vary_mu, n_L_iterations, &
          !      number_of_bands, L_tolerance, tolerance, mu, &
          !      total_energy, expected_reduction,.true.)
       end if
       !call grab_matrix("H", matH,inode)
       !call grab_matrix("KE",matKE,inode)
       !call grab_matrix("NL",matNL,inode)
       !call grab_locps(pseudopotential,n_my_grid_points,inode)
       !!call grab_projs(pseudofunctions,NCF*CORE_SIZE,inode)
       !hartree_energy  = H0  
       !local_ps_energy = L0  
       !xc_energy       = XC0 
       !nl_energy       = NL0 
       !kinetic_energy  = K0  
       !band_energy     = B0
    end if
    ! *** Total Pulay ***
    if (flag_test_all_forces .or. flag_which_force == 3) then
       if (inode == ionode) write (io_lun, *) '*** Full Pulay ***'
       call test_FullPulay(fixed_potential, vary_mu, n_L_iterations, &
                           L_tolerance, tolerance, total_energy,     &
                           expected_reduction)
       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if (flag_self_consistent) then ! Vary only DM and charge density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          ! call force( fixed_potential, vary_mu, n_L_iterations, &
          !      number_of_bands, L_tolerance, tolerance, mu, &
          !      total_energy, expected_reduction,.true.)
       end if
       !call grab_matrix("H", matH,inode)
       !call grab_matrix("KE",matKE,inode)
       !call grab_matrix("NL",matNL,inode)
       !call grab_matrix("K", matK,inode)
       !call grab_locps(pseudopotential,n_my_grid_points,inode)
       !!call grab_projs(pseudofunctions,NCF*CORE_SIZE,inode)
       !hartree_energy  = H0  
       !local_ps_energy = L0  
       !xc_energy       = XC0 
       !nl_energy       = NL0 
       !kinetic_energy  = K0  
       !band_energy     = B0
    end if
    ! *** NSC ***
    if (flag_test_all_forces .or. flag_which_force == 4) then
       if (.not. flag_self_consistent) then
          if (inode == ionode) &
               write(io_lun,*) '*** Non Self-Consistent ***'
          delta_E_hartree = dH0  
          delta_E_xc      = dXC0 
          call set_density()
          call test_nonSC(fixed_potential, vary_mu, n_L_iterations, &
                          L_tolerance, tolerance, total_energy,     &
                          expected_reduction)
          if (flag_test_all_forces) then
             call update_H(fixed_potential)
             if (flag_self_consistent) then
                ! Vary only DM and charge density
                reset_L = .true.
                call new_SC_potl(.true., tolerance, reset_L,  &
                                 fixed_potential, vary_mu,    &
                                 n_L_iterations, L_tolerance, &
                                 total_energy)
             else ! Ab initio TB: vary only DM
                call get_H_matrix(.true., fixed_potential, electrons, &
                                  density, maxngrid)
                call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                               inode, ionode, reset_L, .false.)
                call get_energy(total_energy)
             end if
             !call force( fixed_potential, vary_mu, n_L_iterations,
             !number_of_bands, L_tolerance, tolerance, mu,
             !total_energy, expected_reduction,.true.)
          end if
       end if
    end if
    ! *** phi Pulay Non-local ***
    if (flag_test_all_forces .or. flag_which_force == 5) then
       if (inode == ionode) &
            write (io_lun, *) '*** Non-local phi Pulay ***'
       call test_PhiPulay_nonlocal(fixed_potential, vary_mu,    &
                                   n_L_iterations, L_tolerance, &
                                   tolerance, total_energy,     &
                                   expected_reduction)
       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if (flag_self_consistent) then ! Vary only DM and charge
             ! density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          !call force( fixed_potential, vary_mu, n_L_iterations,
          !number_of_bands, L_tolerance, tolerance, mu, total_energy,
          !expected_reduction,.true.)
       end if
       !call grab_matrix("H", matH,inode)
       !call grab_matrix("KE",matKE,inode)
       !call grab_matrix("NL",matNL,inode)
       !call grab_blips("O",supportfns,inode)
       !hartree_energy  = H0  
       !local_ps_energy = L0  
       !xc_energy       = XC0 
       !nl_energy       = NL0 
       !kinetic_energy  = K0  
       !band_energy     = B0
    end if
    ! *** phi Pulay KE ***
    if (flag_test_all_forces .or. flag_which_force == 6) then
       if (inode == ionode) &
            write(io_lun,*) '*** KE phi Pulay ***'
       call test_PhiPulay_KE(fixed_potential, vary_mu, n_L_iterations, &
                             L_tolerance, tolerance, total_energy,     &
                             expected_reduction)
       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if (flag_self_consistent) then ! Vary only DM and charge
             ! density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          !call force( fixed_potential, vary_mu, n_L_iterations,
          !number_of_bands, L_tolerance, tolerance, mu, total_energy,
          !expected_reduction,.true.)
       end if
       !call grab_matrix("H", matH,inode)
       !call grab_matrix("KE",matKE,inode)
       !call grab_matrix("NL",matNL,inode)
       !call grab_matrix("K", matK,inode)
       !call grab_blips("O",supportfns,inode)
       !call grab_locps(pseudopotential,n_my_grid_points,inode)
       !call grab_charge(density,n_my_grid_points,inode)
       !hartree_energy  = H0  
       !local_ps_energy = L0  
       !xc_energy       = XC0 
       !nl_energy       = NL0 
       !kinetic_energy  = K0  
       !band_energy     = B0
    end if
    ! *** phi Pulay local ***
    if (flag_test_all_forces .or. flag_which_force == 7) then
       if (inode == ionode) &
            write(io_lun,*) '*** Local (Ha, XC, loc ps) phi Pulay ***'
       call test_PhiPulay_Local(fixed_potential, vary_mu,    &
                                n_L_iterations, L_tolerance, &
                                tolerance, total_energy,     &
                                expected_reduction)
       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if (flag_self_consistent) then ! Vary only DM and charge
             ! density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          !call force( fixed_potential, vary_mu, n_L_iterations,
          !number_of_bands, L_tolerance, tolerance, mu, total_energy,
          !expected_reduction,.true.)
       end if
       !call grab_matrix("H", matH,inode)
       !call grab_matrix("KE",matKE,inode)
       !call grab_matrix("NL",matNL,inode)
       !call grab_blips("O",supportfns,inode)
       !hartree_energy  = H0  
       !local_ps_energy = L0  
       !xc_energy       = XC0 
       !nl_energy       = NL0 
       !kinetic_energy  = K0  
       !band_energy     = B0
    end if
    ! *** S Pulay ***
    if (flag_test_all_forces .or. flag_which_force == 8) then
       if (inode == ionode) write(io_lun,*) '*** S Pulay ***'
       call test_SPulay(fixed_potential, vary_mu, n_L_iterations, &
                        L_tolerance, tolerance, total_energy,     &
                        expected_reduction)
       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if (flag_self_consistent) then ! Vary only DM and charge
             ! density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          !call force( fixed_potential, vary_mu, n_L_iterations,
          !number_of_bands, L_tolerance, tolerance, mu, total_energy,
          !expected_reduction,.true.)
       end if
       !call grab_matrix("H", matH,inode)
       !call grab_matrix("KE",matKE,inode)
       !call grab_matrix("NL",matNL,inode)
       !call grab_matrix("K", matK,inode)
       !call grab_blips("O",supportfns, inode)
       !hartree_energy  = H0  
       !local_ps_energy = L0  
       !xc_energy       = XC0 
       !nl_energy       = NL0 
       !kinetic_energy  = K0  
       !band_energy     = B0
    end if
    ! *** cDFT *** 
    if (flag_test_all_forces .or. flag_which_force == 9) then
       if (inode == ionode) write (io_lun, *) '*** cDFT ***'
       call test_cdft(fixed_potential, vary_mu, n_L_iterations, &
                      L_tolerance, tolerance, total_energy,     &
                      expected_reduction)
       if (flag_test_all_forces) then
          call update_H(fixed_potential)
          if(flag_self_consistent) then ! Vary only DM and charge
             ! density
             reset_L = .true.
             call new_SC_potl(.true., tolerance, reset_L,  &
                              fixed_potential, vary_mu,    &
                              n_L_iterations, L_tolerance, &
                              total_energy)
          else ! Ab initio TB: vary only DM
             call get_H_matrix(.true., fixed_potential, electrons, &
                               density, maxngrid)
             call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
                            inode, ionode, reset_L, .false.)
             call get_energy(total_energy)
          end if
          !call force( fixed_potential, vary_mu, n_L_iterations,
          !number_of_bands, L_tolerance, tolerance, mu, total_energy,
          !expected_reduction,.true.)
       end if
    end if

    deallocate(HF_force, STAT=stat)
    if (stat /= 0) call cq_abort("test_forces: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    call stop_timer(tmr_std_moveatoms)

    return
  end subroutine test_forces
!!***


  !!****f* test_force_module/test_FullPulay *
  !!
  !!  NAME 
  !!   test_FullPulay
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the full Pulay forces
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:15, 02/09/2003
  !!  MODIFICATION HISTORY
  !!   2011/11/29 L.Tong
  !!     Added spin polarisation
  !!   2011/12/09 L.Tong
  !!     Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!  2012/04/03 10:13 dave
  !!   Update for analytic blips
  !!  SOURCE
  !!
  subroutine test_FullPulay(fixed_potential, vary_mu, n_L_iterations, &
                            L_tolerance, tolerance, total_energy,     &
                            expected_reduction)

    use datatypes
    use numbers
    use move_atoms,      only: primary_update, cover_update,         &
                               update_atom_coord
    use group_module,    only: parts
    use cover_module,    only: BCS_parts, DCS_parts
    use primary_module,  only: bundle
    use global_module,   only: iprint_MD, x_atom_cell, y_atom_cell,  &
                               z_atom_cell, id_glob_inv, WhichPulay, &
                               BothPulay, flag_self_consistent,      &
                               flag_basis_set, PAOs, blips,          &
                               ni_in_cell, nspin, flag_analytic_blip_int
    use energy,          only: local_ps_energy, hartree_energy,      &
                               xc_energy, get_energy, band_energy
    use force_module,    only: pulay_force, Pulay, HF_and_Pulay,     &
                               get_HF_non_local_force, get_KE_force
    use GenComms,        only: myid, inode, ionode, cq_abort
    use H_matrix_module, only: get_H_matrix
    use S_matrix_module, only: get_S_matrix
    use mult_module,     only: matK, allocate_temp_matrix,           &
                               free_temp_matrix, matrix_sum
    use matrix_data,     only: Hrange
    use SelfCon,         only: new_SC_potl
    use DMMin,           only: FindMinDM
    use density_module,  only: density
    use maxima_module,   only: maxngrid
    use memory_module,   only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L
    integer :: n_L_iterations

    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    logical      :: reset_L
    integer      :: spin, stat
    real(double) :: E0, F0, E1, F1, analytic_force, numerical_force
    integer,      dimension(nspin)        :: matKold
    real(double), dimension(nspin)        :: electrons
    real(double), dimension(:,:), allocatable :: p_force, HF_NL_force, KE_force

    allocate(p_force(3,ni_in_cell), HF_NL_force(3,ni_in_cell), &
             KE_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_FullPulay: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 9*ni_in_cell, type_dbl)

    ! Warn user that we're NOT getting just HF with PAOs !
    if (flag_basis_set==PAOs) then
       if (myid==0) write(io_lun,&
            fmt='(2x,"********************************************************")')
       if(myid==0) write(io_lun,&
            fmt='(2x,"*                                                      *")')
       if(myid==0) write(io_lun,&
            fmt='(2x,"*    WARNING * WARNING * WARNING * WARNING * WARNING   *")')
       if(myid==0) write(io_lun,&
            fmt='(2x,"*                                                      *")')
       if(myid==0) write(io_lun,&
            fmt='(2x,"* With a PAO basis we calculate NL HF AND Pulay forces *")')
       if(myid==0) write(io_lun,&
            fmt='(2x,"*                                                      *")')
       if(myid==0) write(io_lun,&
            fmt='(2x,"********************************************************")')
    end if
    ! We're coming in from initial_H: assume that initial E found
    E0 = total_energy
    ! Find force
    WhichPulay = BothPulay
    call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
                     n_L_iterations, L_tolerance, tolerance, &
                     total_energy, expected_reduction, ni_in_cell)
    ! This routine deals with the movement of the nonlocal
    ! pseudopotential.
    if (flag_basis_set == PAOs) then
       call get_HF_non_local_force(HF_NL_force, HF_and_Pulay, &
                                   ni_in_cell)
    else if (flag_basis_set==blips) then
       call get_HF_non_local_force(HF_NL_force, Pulay, ni_in_cell)
    end if
    ! Get the kinetic energy force component
    do spin = 1, nspin
       matKold(spin) = allocate_temp_matrix(Hrange, 0)
       call matrix_sum(zero, matKold(spin), one, matK(spin))
    end do
    ! get_KE_force changes matK, that is why we need to use matKold to
    ! get back the correct matK
    if(.NOT.flag_analytic_blip_int) call get_KE_force(KE_force, ni_in_cell)
    do spin = 1, nspin
       call matrix_sum(zero, matK(spin), one, matKold(spin))
    end do
    ! end LT 2011/11/29
    ! Store local energy
    ! E0 = hartree_energy + xc_energy + local_ps_energy
    ! Find out direction and atom for displacement
    if (myid == 0) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," in direction ",&
                             &i2," by ",f10.6," bohr")') &
               TF_atom_moved, TF_direction, TF_delta
    F0 = p_force(TF_direction,TF_atom_moved) +      &
         HF_NL_force(TF_direction, TF_atom_moved) + &
         KE_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Initial energy: ",f20.12,/,2x,&
                             &"Initial pulay force: ",f20.12)') E0, F0
    ! Move the specified atom 
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    ! Regenerate S
    call get_S_matrix(inode, ionode)
   if (flag_self_consistent) then ! Vary only DM and charge density
      reset_L = .true.
      call new_SC_potl(.true., tolerance, reset_L, fixed_potential, &
                       vary_mu, n_L_iterations, L_tolerance,        &
                       total_energy)
    else ! Ab initio TB: vary only DM
       call get_H_matrix(.true., fixed_potential, electrons, density, &
                         maxngrid)
       call FindMinDM(n_L_iterations, vary_mu, L_tolerance, inode, &
                      ionode, reset_L, .false.)
       call get_energy(total_energy)
    end if
    ! Note that we've held K fixed but allow potential to vary ? Yes:
    ! this way we get h_on_support in workspace_support
    ! call get_H_matrix(.false., fixed_potential, electrons, potential,
    !                   density, pseudopotential, N_GRID_MAX)
    ! call get_energy(total_energy)
    ! Find force
    call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
                     n_L_iterations, L_tolerance, tolerance, &
                     total_energy, expected_reduction, ni_in_cell)
    if (flag_basis_set == PAOs) then
       call get_HF_non_local_force(HF_NL_force, HF_and_Pulay, &
                                   ni_in_cell)
    else if (flag_basis_set == blips) then
       call get_HF_non_local_force(HF_NL_force, Pulay, ni_in_cell)
    end if
    ! Get the kinetic energy force component
    do spin = 1, nspin
       call matrix_sum(zero, matKold(spin), one, matK(spin))
    end do
    if(.NOT.flag_analytic_blip_int) call get_KE_force(KE_force, ni_in_cell)
    do spin = 1, nspin
       call matrix_sum(zero, matK(spin), one, matKold(spin))
    end do
    do spin = nspin, 1, -1
       call free_temp_matrix(matKold(spin))
    end do
    ! E1 = hartree_energy + xc_energy + local_ps_energy
    E1 = total_energy
    F1 = p_force(TF_direction,TF_atom_moved) +     &
         HF_NL_force(TF_direction,TF_atom_moved) + &
         KE_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Final energy: ",f20.12,/,2x,&
                             &"Final pulay force: ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = half * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical Force: ",f20.12,/,2x,&
                             &"Analytic Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    ! Regenerate S
    call get_S_matrix(inode, ionode)
    ! Note that we've held K fixed but allow potential to vary ? Yes:
    ! this way we get h_on_support in workspace_support
    call get_H_matrix(.true., fixed_potential, electrons, density, &
                      maxngrid)
    call get_energy(total_energy)

    deallocate(p_force, HF_NL_force, KE_force, STAT=stat)
    if (stat /= 0) call cq_abort("test_FullPulay: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 9*ni_in_cell, type_dbl)

    return
  end subroutine test_FullPulay
  !!***


  !!****f* test_force_module/test_HF *
  !!
  !!  NAME 
  !!   test_HF
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the Hellmann-Feynman non-local pseudopotential force
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:40, 02/09/2003
  !!  MODIFICATION HISTORY
  !!   2011/11/29 L.Tong
  !!     - Added spin polarisation
  !!     - fixed memory leak due to not freeing memories
  !!   2011/12/09 L.Tong
  !!     Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_HF(fixed_potential, vary_mu, n_L_iterations, &
                     L_tolerance, tolerance, total_energy,     &
                     expected_reduction)

    use datatypes
    use numbers
    use move_atoms,             only: primary_update, cover_update, &
                                      update_atom_coord
    use group_module,           only: parts
    use cover_module,           only: BCS_parts, DCS_parts
    use primary_module,         only: bundle
    use global_module,          only: iprint_MD, x_atom_cell,       &
                                      y_atom_cell, z_atom_cell,     &
                                      id_glob_inv,                  &
                                      flag_self_consistent,         &
                                      flag_basis_set, PAOs, blips,  &
                                      ni_in_cell,                   &
                                      nspin, spin_factor
    use pseudopotential_data,   only: init_pseudo
    use pseudo_tm_module,       only: set_tm_pseudo,                &
                                      loc_pp_derivative_tm
    use pseudopotential_common, only: pseudo_type, OLDPS, SIESTA,   &
                                      STATE, ABINIT, core_correction
    use energy,                 only: nl_energy, get_energy,        &
                                      local_ps_energy, band_energy, &
                                      kinetic_energy
    use force_module,           only: get_HF_force,                 &
                                      get_HF_non_local_force,       &
                                      HF_and_Pulay, HF
    use GenComms,               only: myid, inode, ionode, cq_abort
    use H_matrix_module,        only: get_H_matrix
    use density_module,         only: set_density,                  &
                                      get_electronic_density, density
    use functions_on_grid,      only: supportfns, H_on_supportfns
    use maxima_module,          only: maxngrid
    use memory_module,          only: reg_alloc_mem, reg_dealloc_mem, &
                                      type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L
    integer :: n_L_iterations
    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer      :: stat, spin
    real(double) :: E0, F0, E1, F1, analytic_force, numerical_force, &
                    Enl0, Enl1, Fnl0, Fnl1, KE0
    real(double), dimension(nspin)          :: electrons
    real(double), dimension(:,:), allocatable :: HF_force
    real(double), dimension(:,:), allocatable :: HF_NL_force
    real(double), dimension(:),   allocatable :: density_out_tot
    real(double), dimension(:,:), allocatable :: density_out

    allocate(HF_force(3,ni_in_cell), HF_NL_force(3,ni_in_cell), &
             density_out_tot(maxngrid), density_out(maxngrid,nspin), &
             STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_HF: Error alloc mem: ", ni_in_cell, maxngrid)
    call reg_alloc_mem(area_moveatoms, 6*ni_in_cell+(1+nspin)*maxngrid, &
                       type_dbl)
    
    ! We're coming in from initial_H: assume that initial E found
    
    ! Warn user that we're NOT getting just HF with PAOs !
    if (flag_basis_set == PAOs) then
       if (myid == 0) write (io_lun, &
            fmt='(2x,"********************************************************")')
       if (myid == 0) write (io_lun, &
            fmt='(2x,"*                                                      *")')
       if (myid == 0) write (io_lun, &
            fmt='(2x,"*    WARNING * WARNING * WARNING * WARNING * WARNING   *")')
       if (myid == 0) write (io_lun, &
            fmt='(2x,"*                                                      *")')
       if(myid == 0) write (io_lun, &
            fmt='(2x,"* With a PAO basis we calculate NL HF AND Pulay forces *")')
       if(myid == 0) write (io_lun, &
            fmt='(2x,"*                                                      *")')
       if(myid == 0) write (io_lun, &
            fmt='(2x,"********************************************************")')
    end if
    ! We're coming in from initial_H: assume that initial E found
    ! Non-local
    Enl0 = nl_energy
    ! Full band energy
    E0 = band_energy
    ! Store KE for later correction
    KE0 = kinetic_energy
    ! Find force: local
    ! If necessary, find output density (for HF)
    if (.not. flag_self_consistent) then ! Harris-Foulkes requires
       ! output density
       call get_electronic_density(density_out, electrons, supportfns,&
                                   H_on_supportfns(1), inode, ionode, &
                                   maxngrid)
    else
       do spin = 1, nspin 
          density_out(:,spin) = density(:,spin)  
       end do
    end if
    density_out_tot = spin_factor * sum(density_out, 2)

    select case (pseudo_type)
     case (OLDPS)
      call get_HF_force(hf_force, density_out_tot, ni_in_cell, maxngrid)
     case (SIESTA)
      call loc_pp_derivative_tm(hf_force, density_out_tot, maxngrid)
     case (ABINIT)
      call loc_pp_derivative_tm(hf_force, density_out_tot, maxngrid)
    end select
    ! Find force
    if (flag_basis_set == PAOs) then
       call get_HF_non_local_force(HF_NL_force, HF_and_Pulay, &
                                   ni_in_cell)
    else if (flag_basis_set == blips) then
       call get_HF_non_local_force(HF_NL_force, HF, ni_in_cell)
    end if
    ! Find out direction and atom for displacement
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," in direction ",&
                             &i2," by ",f10.6," bohr")') &
               TF_atom_moved, TF_direction, TF_delta
    Fnl0 = HF_NL_force(TF_direction,TF_atom_moved)
    if(inode == ionode) &
         write (io_lun,fmt='(2x,"Initial NL energy: ", f20.12,/,2x,&
                             &"Initial NL force : ",f20.12)') Enl0, Fnl0
    ! LT 2011/11/29: this is redundant???
    ! F0 = HF_NL_force(TF_direction,TF_atom_moved) + &
    !      HF_force(TF_direction,TF_atom_moved)
    ! end LT 2011/11/29
    F0 = HF_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Initial band energy: ",f20.12,/,2x,&
                             &"Initial HF force : ",f20.12)') E0, F0
    ! Move the specified atom 
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    ! Recalculate pseudopotentials
    select case (pseudo_type) 
    case (OLDPS)
       call init_pseudo(core_correction)
    case (SIESTA)
       call set_tm_pseudo
    case (ABINIT)
       call set_tm_pseudo
    end select
    ! Note that we've held K and |phi> fixed
    ! Calculate new energy
    call get_H_matrix(.true., fixed_potential, electrons, density, &
                      maxngrid)
    call get_energy(total_energy)
    Enl1 = nl_energy
    E1 = band_energy - kinetic_energy + KE0 ! Fix change in KE
    E1 = band_energy - kinetic_energy - nl_energy + Enl0 + KE0 ! Fix change
    ! in KE AND NL
    ! Find force
    ! Now that the atoms have moved, calculate the terms again
    select case (pseudo_type)
     case(OLDPS)
        call get_HF_force(hf_force, density_out_tot, ni_in_cell, maxngrid)
     case(SIESTA)
        call loc_pp_derivative_tm(hf_force, density_out_tot, maxngrid)
     case(ABINIT)
        call loc_pp_derivative_tm(hf_force, density_out_tot, maxngrid)
    end select
    ! Find force
    if (flag_basis_set == PAOs) then
       call get_HF_non_local_force(HF_NL_force, HF_and_Pulay, &
                                   ni_in_cell)
    else if(flag_basis_set == blips) then
       call get_HF_non_local_force(HF_NL_force, HF, ni_in_cell)
    end if
    Fnl1 = HF_NL_force(TF_direction,TF_atom_moved)
    ! LT 2011/11/29: is this redundant?
    ! F1 = HF_NL_force(TF_direction,TF_atom_moved) + &
    !      HF_force(TF_direction,TF_atom_moved)
    ! end LT 2011/11/29
    F1 = HF_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write(io_lun,fmt='(2x,"Final NL energy: ",f20.12,/,2x,"Final &
                            &NL force : ",f20.12)') Enl1, Fnl1
    numerical_force = -(Enl1 - Enl0) / TF_delta
    analytic_force = half * (Fnl1 + Fnl0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical NL Force: ",f20.12,/,2x,&
                             &"Analytic NL Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Final band energy: ",f20.12,/,2x,&
                             &"Final HF force : ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = half * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical local Force: ",f20.12,/,2x,&
                             &"Analytic local Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode==ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    ! Now regenerate the pseudos and h_on_support
    select case (pseudo_type) 
    case (OLDPS)
       call init_pseudo(core_correction)
    case (SIESTA)
       call set_tm_pseudo
    case (ABINIT)
       call set_tm_pseudo
    end select
    ! Restore h_on_support
    ! Calculate new energy
    call get_H_matrix(.true., fixed_potential, electrons, density, &
                      maxngrid)
    call get_energy(total_energy)

    deallocate(HF_force, HF_NL_force, density_out_tot, density_out, STAT=stat)
    if (stat /= 0) call cq_abort("test_HF: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 6*ni_in_cell+(1+nspin)*maxngrid, &
                         type_dbl)

    return
  end subroutine test_HF
  !!***


  !!****f* test_force_module/test_PhiPulay_nonlocal *
  !!
  !!  NAME 
  !!   test_PhiPulay_nonlocal
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the phi Pulay non-local pseudopotential force
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:15, 02/09/2003
  !!  MODIFICATION HISTORY
  !!   2011/11/29 L.Tong
  !!     Added spin polarisation
  !!   2011/12/09 L.Tong
  !!     Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_PhiPulay_nonlocal(fixed_potential, vary_mu,    &
                                    n_L_iterations, L_tolerance, &
                                    tolerance, total_energy,     &
                                    expected_reduction)

    use numbers
    use datatypes
    use move_atoms,                 only: primary_update,              &
                                          cover_update,                &
                                          update_atom_coord
    use group_module,               only: parts
    use cover_module,               only: BCS_parts, DCS_parts
    use primary_module,             only: bundle
    use global_module,              only: iprint_MD, x_atom_cell,      &
                                          y_atom_cell, z_atom_cell,    &
                                          id_glob_inv, flag_basis_set, &
                                          blips, PAOs, ni_in_cell,     &
                                          nspin
    use energy,                     only: nl_energy, get_energy
    use force_module,               only: get_HF_non_local_force,      &
                                          Pulay
    use GenComms,                   only: myid, inode, ionode, cq_abort
    use H_matrix_module,            only: get_H_matrix
    use blip_grid_transform_module, only: blip_to_support_new
    use functions_on_grid,          only: supportfns
    use density_module,             only: density
    use maxima_module,              only: maxngrid
    use memory_module,              only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L
    integer :: n_L_iterations
    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer :: stat
    real(double) ::  E0, F0, E1, F1, analytic_force, numerical_force
    real(double), dimension(nspin) :: electrons
    real(double), dimension(:,:), allocatable :: HF_NL_force

    allocate(HF_NL_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_PhiPulay_nonlocal: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    ! Warn user that we're don't get NL phi Pulay with PAOs !
    if (flag_basis_set == PAOs) then
       if (myid == 0) &
            write (io_lun,fmt='(2x,"*****************************&
                                &********************************")')
       if (myid == 0) &
            write (io_lun,fmt='(2x,"*                           &
                                &                                *")')
       if (myid == 0) &
            write (io_lun,fmt='(2x,"* WARNING * WARNING * WARNING &
                                &* WARNING * WARNING * WARNING *")')
       if (myid == 0) &
            write (io_lun,fmt='(2x,"*                           &
                                &                                *")')
       if (myid == 0) &
            write (io_lun,fmt='(2x,"* With a PAO basis we DO &
                                &NOT calculate NL phi Pulay forces  *")')
       if (myid == 0) &
            write (io_lun,fmt='(2x,"*                           &
                                &                                *")')
       if (myid == 0) &
            write (io_lun,fmt='(2x,"*****************************&
                                &********************************")')
       return
    end if
    ! We're coming in from initial_H: assume that initial E found
    ! Find force
    call get_HF_non_local_force(HF_NL_force, Pulay, ni_in_cell)
    ! Store local energy
    E0 = nl_energy
    ! Find out direction and atom for displacement
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," in direction ",&
                             &i2," by ",f10.6," bohr")') &
               TF_atom_moved, TF_direction,TF_delta
    F0 = HF_NL_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write(io_lun,fmt='(2x,"Initial NL energy: ",f20.12,/,2x,&
                            &"Initial NL force : ",f20.12)') E0, F0
    ! Move the specified atom 
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    if (flag_basis_set == blips) then
       ! Reproject blips
       call blip_to_support_new(inode-1, supportfns)    
    end if
    ! Note that we've held K and |chi> fixed
    ! Calculate new energy
    call get_H_matrix(.true., fixed_potential, electrons, &
                      density, maxngrid)
    call get_energy(total_energy)
    ! Find force
    call get_HF_non_local_force(HF_NL_force, Pulay, ni_in_cell)
    E1 = nl_energy
    F1 = HF_NL_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write(io_lun,fmt='(2x,"Final NL energy: ",f20.12,/,2x,"Final &
                            &NL force : ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = half * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical Force: ",f20.12,/,2x,&
                             &"Analytic Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)

    deallocate(HF_NL_force, STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_PhiPulay_nonlocal: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    return
  end subroutine test_PhiPulay_nonlocal
  !!***


  !!****f* test_force_module/test_PhiPulay_KE *
  !!
  !!  NAME 
  !!   test_PhiPulay_KE
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the phi Pulay non-local pseudopotential force
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:15, 02/09/2003
  !!  MODIFICATION HISTORY
  !!  L.Tong
  !!   Removed redundant parameter number_of_bands
  !!  2012/03/27 L.Tong
  !!  - Changed spin implementation
  !!  - Removed redundant input parameter real(double) mu
  !!  2012/04/03 10:13 dave
  !!   Update for analytic blips
  !!  SOURCE
  !!
  subroutine test_PhiPulay_KE(fixed_potential, vary_mu,               &
                              n_L_iterations, L_tolerance, tolerance, &
                              total_energy, expected_reduction)

    use datatypes
    use numbers
    use move_atoms,                 only: primary_update, &
                                          cover_update, &
                                          update_atom_coord
    use group_module,               only: parts
    use cover_module,               only : BCS_parts, DCS_parts
    use primary_module,             only : bundle
    use global_module,              only: iprint_MD, x_atom_cell, &
                                          y_atom_cell, z_atom_cell, &
                                          id_glob_inv, flag_basis_set,&
                                          blips, ni_in_cell, &
                                          nspin, flag_analytic_blip_int
    use energy,                     only: kinetic_energy, get_energy
    use force_module,               only: get_KE_force, pulay_force
    use GenComms,                   only: myid, inode, ionode, cq_abort
    use H_matrix_module,            only: get_H_matrix
    use blip_grid_transform_module, only: blip_to_support_new
    use mult_module,                only: matK, allocate_temp_matrix, &
                                          free_temp_matrix, matrix_sum
    use matrix_data,                only: Hrange
    use functions_on_grid,          only: supportfns
    use density_module,             only: density
    use maxima_module,              only: maxngrid
    use memory_module,              only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L
    integer :: n_L_iterations

    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer      :: spin, stat
    real(double) ::  E0, F0, E1, F1, analytic_force, numerical_force
    integer,      dimension(nspin) :: matKold
    real(double), dimension(nspin) :: electrons
    real(double), dimension(:,:), allocatable :: KE_force, p_force

    allocate(KE_force(3,ni_in_cell), p_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_PhiPulay_KE: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 6*ni_in_cell, type_dbl)

    ! We're coming in from initial_H: assume that initial E found
    ! Find force
    do spin = 1, nspin
       matKold(spin) = allocate_temp_matrix(Hrange,0)
       call matrix_sum(zero, matKold(spin), one, matK(spin))
    end do
    if(flag_analytic_blip_int) then
       call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
            n_L_iterations, L_tolerance, tolerance, &
            total_energy, expected_reduction, ni_in_cell)       
    else
       call get_KE_force(KE_force, ni_in_cell)
    end if
    do spin = 1, nspin
       call matrix_sum(zero, matK(spin), one, matKold(spin))
    end do
    do spin = nspin, 1, -1
       call free_temp_matrix(matKold(spin))
    end do
    ! Store local energy
    E0 = kinetic_energy
    ! Find out direction and atom for displacement
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," in direction ",i2,&
                             &" by ",f10.6," bohr")') &
               TF_atom_moved, TF_direction,TF_delta
    F0 = KE_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Initial KE energy: ",f20.12,/,2x,&
                             &"Initial KE force : ",f20.12)') E0, F0
    ! Move the specified atom 
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    if (flag_basis_set == blips) then
       ! Reproject blips
       call blip_to_support_new(inode-1, supportfns)    
    end if
    ! Note that we've held K and |chi> fixed
    ! Calculate new energy
    call get_H_matrix(.true., fixed_potential, electrons, density, &
                      maxngrid)
    call get_energy (total_energy)
    ! Find force
    if(flag_analytic_blip_int) then
       call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
            n_L_iterations, L_tolerance, tolerance, &
            total_energy, expected_reduction, ni_in_cell)       
    else
       call get_KE_force(KE_force, ni_in_cell)
    end if
    E1 = kinetic_energy
    F1 = KE_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Final KE energy: ",f20.12,/,&
                             &2x,"Final KE force : ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = half * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical Force: ",f20.12,/,2x,&
                             &"Analytic Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)

    deallocate(KE_force, p_force, STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_PhiPulay_KE: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 6*ni_in_cell, type_dbl)

    return
  end subroutine test_PhiPulay_KE
  !!***


  !!****f* test_force_module/test_PhiPulay_local *
  !!
  !!  NAME 
  !!   test_PhiPulay_local
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the phi Pulay non-local pseudopotential force
  !!
  !!   This is a knotty one ! The most general way to get the energy
  !!   change is to consider the Harris-Foulkes functional (which 
  !!   at self-consistency is the same as the "alternative" energy),
  !!   whose value is 2Tr[KH].  Then we fix K_{ij} but allow H_{ij} 
  !!   to vary (i.e. rebuild from \hat{H} and |\phi_i>).  BUT because
  !!   we're being careful and looking at all the different bits 
  !!   separately we'll FIX the NL and KE parts of H (they don't
  !!   matter).
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:15, 02/09/2003
  !!  MODIFICATION HISTORY
  !!   08:20, 2003/09/03 dave
  !!    Changed to use the change in 2Tr[K.H].
  !!   2011/11/29 L.Tong
  !!    Added spin polarisation
  !!   2011/12/09 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - removed input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_PhiPulay_local(fixed_potential, vary_mu,    &
                                 n_L_iterations, L_tolerance, &
                                 tolerance, total_energy,     &
                                 expected_reduction)

    use datatypes
    use move_atoms,                 only: primary_update, &
                                          cover_update, &
                                          update_atom_coord
    use group_module,               only: parts
    use cover_module,               only: BCS_parts, DCS_parts
    use primary_module,             only: bundle
    use global_module,              only: iprint_MD, x_atom_cell, &
                                          y_atom_cell, z_atom_cell, &
                                          id_glob_inv, WhichPulay, &
                                          PhiPulay, flag_basis_set, &
                                          blips, PAOs, ni_in_cell, &
                                          nspin
    use energy,                     only: local_ps_energy, &
                                          hartree_energy, xc_energy, &
                                          get_energy, band_energy
    use force_module,               only: pulay_force
    use GenComms,                   only: myid, inode, ionode, cq_abort
    use H_matrix_module,            only: get_H_matrix
    use blip_grid_transform_module, only: blip_to_support_new
    use PAO_grid_transform_module,  only: PAO_to_grid
    use functions_on_grid,          only: supportfns
    use density_module,             only: density
    use maxima_module,              only: maxngrid
    use memory_module,              only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L
    integer :: n_L_iterations

    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer :: stat
    real(double) ::  E0, F0, E1, F1, analytic_force, numerical_force
    real(double), dimension(nspin) :: electrons
    real(double), dimension(:,:), allocatable :: p_force, KE_force

    allocate(KE_force(3,ni_in_cell), p_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_PhiPulay_local: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 6*ni_in_cell, type_dbl)

    ! We're coming in from initial_H: assume that initial E found
    ! Find force
    call get_H_matrix(.false., fixed_potential, electrons, density, &
                      maxngrid)

    WhichPulay = PhiPulay
    call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
                     n_L_iterations, L_tolerance, tolerance, &
                     total_energy, expected_reduction, ni_in_cell)
    ! Store local energy
!    E0 = hartree_energy + xc_energy + local_ps_energy
    E0 = band_energy
    ! Find out direction and atom for displacement
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," in direction ",&
                             &i2," by ",f10.6," bohr")')&
               TF_atom_moved, TF_direction,TF_delta
    F0 = p_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Initial energy: ",f20.12,/,2x,&
                            &"Initial phi Pulay force : ",f20.12)') E0, F0
    ! Move the specified atom 
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    if (flag_basis_set == blips) then
       ! Reproject blips
       call blip_to_support_new(inode-1, supportfns)    
    else if (flag_basis_set == PAOs) then
       ! Regenerate support with a call to PAO_to_grid
       call PAO_to_grid(inode-1, supportfns)
    end if
    ! Note that we've held K fixed but allow potential to vary ? Yes:
    ! this way we get h_on_support in workspace_support
    call get_H_matrix(.false., fixed_potential, electrons, density, &
                      maxngrid)
    call get_energy(total_energy)
    ! Find force
    call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
                     n_L_iterations, L_tolerance, tolerance, &
                     total_energy, expected_reduction, ni_in_cell)
!    E1 = hartree_energy + xc_energy + local_ps_energy
    E1 = band_energy
    F1 = p_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write(io_lun,fmt='(2x,"Final energy: ",f20.12,/,2x,&
                            &"Final phi Pulay force: ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = 0.5_double * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical Force: ",f20.12,/,2x,&
                             &"Analytic Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)

    deallocate(KE_force, p_force, STAT=stat)
    if (stat /= 0) call cq_abort("test_PhiPulay_local: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 6*ni_in_cell, type_dbl)

    return
  end subroutine test_PhiPulay_local
  !!***


  !!****f* test_force_module/test_SPulay *
  !!
  !!  NAME 
  !!   test_SPulay
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the phi Pulay non-local pseudopotential force
  !!
  !!   This is a knotty one ! The most general way to get the energy
  !!   change is to consider the Harris-Foulkes functional (which 
  !!   at self-consistency is the same as the "alternative" energy),
  !!   whose value is 2Tr[KH].  Then we fix K_{ij} but allow H_{ij} 
  !!   to vary (i.e. rebuild from \hat{H} and |\phi_i>).  BUT because
  !!   we're being careful and looking at all the different bits 
  !!   separately we'll FIX the NL and KE parts of H (they don't
  !!   matter).
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:15, 02/09/2003
  !!  MODIFICATION HISTORY
  !!   08:20, 2003/09/03 dave
  !!    Changed to use the change in 2Tr[K.H].
  !!   2011/11/29 L.Tong
  !!    - Removed redundant dependence on H_matrix_module
  !!      updated calls to pulay_force, new_SC_potl and FindMinDM
  !!    - Removed redundant parameter on number_of_bands
  !!    - Removed redundant local variable electrons
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - removed input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_SPulay(fixed_potential, vary_mu, n_L_iterations, &
                         L_tolerance, tolerance, total_energy,     &
                         expected_reduction)

    use datatypes
    use move_atoms,      only: primary_update, cover_update,        &
                               update_atom_coord
    use group_module,    only: parts
    use cover_module,    only: BCS_parts, DCS_parts
    use primary_module,  only: bundle
    use global_module,   only: iprint_MD, x_atom_cell, y_atom_cell, &
                               z_atom_cell, id_glob_inv,            &
                               flag_self_consistent, ni_in_cell
    use energy,          only: local_ps_energy, hartree_energy,     &
                               xc_energy, get_energy, band_energy
    use force_module,    only: pulay_force
    use GenComms,        only: myid, inode, ionode, cq_abort
    use S_matrix_module, only: get_S_matrix
    use global_module,   only: WhichPulay, SPulay, flag_basis_set,  &
                               blips, PAOs
    use SelfCon,         only: new_SC_potl
    use DMMin,           only: FindMinDM
    use memory_module,   only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L
    integer :: n_L_iterations

    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer :: stat
    real(double) ::  E0, F0, E1, F1, analytic_force, numerical_force
    real(double), dimension(:,:), allocatable :: p_force, KE_force
    logical :: reset_L

    allocate(KE_force(3,ni_in_cell), p_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_SPulay: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 6*ni_in_cell, type_dbl)

    ! We're coming in from initial_H: assume that initial E found
    ! Find force
    WhichPulay = SPulay
    call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
                     n_L_iterations, L_tolerance, tolerance, &
                     total_energy, expected_reduction, ni_in_cell)
    ! Store local energy
    ! E0 = hartree_energy + xc_energy + local_ps_energy
    E0 = band_energy
    ! Find out direction and atom for displacement
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," in direction ",i2,"&
                             & by ",f10.6," bohr")') &
               TF_atom_moved, TF_direction, TF_delta
    F0 = p_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Initial energy: ",f20.12,/,2x,&
                             &"Initial S-pulay force: ",f20.12)') E0, F0
    ! Move the specified atom 
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    ! Regenerate S
    call get_S_matrix(inode, ionode)
    ! Now we diagonalise
   if (flag_self_consistent) then ! Vary only DM and charge density
      reset_L = .true.
      call new_SC_potl(.true., tolerance, reset_L, fixed_potential, &
                       vary_mu, n_L_iterations, L_tolerance,        &
                       total_energy)
    else ! Ab initio TB: vary only DM
!       call get_H_matrix(.true., fixed_potential, electrons,
       !       potential, density, pseudopotential, N_GRID_MAX)
       call FindMinDM(n_L_iterations, vary_mu, L_tolerance, inode, &
                      ionode, reset_L, .false.)
       call get_energy(total_energy)
    end if
    ! Project H onto support
    ! call get_H_matrix(.false., fixed_potential, electrons, potential,
    !                   density, pseudopotential, N_GRID_MAX)
    ! call get_energy(total_energy)
    ! Find force
    if (flag_basis_set == PAOs) then
       ! Move the specified atom back
       if (TF_direction == 1) then
          x_atom_cell(id_glob_inv(TF_atom_moved)) = &
               x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
       else if (TF_direction == 2) then
          y_atom_cell(id_glob_inv(TF_atom_moved)) = &
               y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
       else if (TF_direction == 3) then
          z_atom_cell(id_glob_inv(TF_atom_moved)) = &
               z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
       end if
       call update_atom_coord
       ! Update positions and indices
       call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                           bundle, parts, myid)
       call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                         BCS_parts, parts)
       call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                         DCS_parts, parts)
    end if
    call pulay_force(p_force, KE_force, fixed_potential, vary_mu,      &
                     n_L_iterations, L_tolerance, tolerance, &
                     total_energy, expected_reduction, ni_in_cell)
!    E1 = hartree_energy + xc_energy + local_ps_energy
    E1 = band_energy
    F1 = p_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Final energy: ",f20.12,/,2x,"Final S-&
                            &pulay force: ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = 0.5_double * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical Force: ",f20.12,/,2x,&
                             &"Analytic Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (flag_basis_set == blips) then
       if (TF_direction == 1) then
          x_atom_cell(id_glob_inv(TF_atom_moved)) = &
               x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
       else if (TF_direction == 2) then
          y_atom_cell(id_glob_inv(TF_atom_moved)) = &
               y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
       else if (TF_direction == 3) then
          z_atom_cell(id_glob_inv(TF_atom_moved)) = &
               z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
       end if
       call update_atom_coord
       ! Update positions and indices
       call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                           bundle, parts, myid)
       call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                         BCS_parts, parts)
       call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                         DCS_parts, parts)
    end if

    deallocate(KE_force, p_force, STAT=stat)
    if (stat /= 0) call cq_abort("test_SPulay: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 6*ni_in_cell, type_dbl)

    return
  end subroutine test_SPulay
  !!***


  !!****f* test_force_module/test_nonSC *
  !!
  !!  NAME 
  !!   test_nonSC
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the non-SC force.
  !!
  !!   We're interested in how the band energy and the Hartree and XC
  !!   correction energies vary as the pseudo-atomic densities
  !!   move over the grid.
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   08:32, 2003/09/03 dave
  !!  MODIFICATION HISTORY
  !!   2011/11/29 L.Tong
  !!     Added spin polarisation
  !!   2011/12/09 L.Tong
  !!     Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_nonSC(fixed_potential, vary_mu, n_L_iterations, &
                        L_tolerance, tolerance, total_energy,     &
                        expected_reduction)

    use datatypes
    use numbers
    use move_atoms,        only: primary_update, cover_update,         &
                                 update_atom_coord
    use group_module,      only: parts
    use cover_module,      only: BCS_parts, DCS_parts
    use primary_module,    only: bundle
    use global_module,     only: iprint_MD, x_atom_cell, y_atom_cell,  &
                                 z_atom_cell, id_glob_inv, ni_in_cell, &
                                 nspin
    use energy,            only: local_ps_energy, hartree_energy,      &
                                 xc_energy, get_energy, band_energy,   &
                                 delta_E_hartree, delta_E_xc
    use force_module,      only: get_nonSC_correction_force
    use GenComms,          only: myid, inode, ionode, cq_abort
    use H_matrix_module,   only: get_H_matrix
    use density_module,    only: set_density, get_electronic_density,  &
                                 density
    use functions_on_grid, only: supportfns, H_on_supportfns
    use maxima_module,     only: maxngrid
    use memory_module,     only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical      :: vary_mu, find_chdens, fixed_potential
    logical      :: start, start_L
    integer      :: n_L_iterations
    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer      :: stat
    real(double) :: E0, F0, E1, F1, analytic_force, numerical_force, &
                    electrons_tot
    real(double), dimension(nspin) :: electrons
    real(double), dimension(:,:), allocatable :: nonSC_force
    real(double), dimension(:,:), allocatable :: density_out
    
    allocate(nonSC_force(3,ni_in_cell), density_out(maxngrid,nspin), &
             STAT=stat)
    if (stat /= 0) &
         call cq_abort("test_nonSC: Error alloc mem: ", ni_in_cell, maxngrid)
    call reg_alloc_mem(area_moveatoms, 3*ni_in_cell+nspin*maxngrid, type_dbl)

    ! We're coming in from initial_H: assume that initial E found
    call get_electronic_density(density_out, electrons, supportfns, &
                                H_on_supportfns(1), inode, ionode,  &
                                maxngrid)
    ! Find force
    call get_nonSC_correction_force(nonSC_force, density_out, inode, &
                                    ionode, ni_in_cell, maxngrid)
    ! Store local energy
    E0 = band_energy + delta_E_hartree + delta_E_xc
    ! Find out direction and atom for displacement
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," &
                             &in direction ",i2," by ",f10.6," bohr")') &
               TF_atom_moved, TF_direction, TF_delta
    F0 = nonSC_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Initial energy: ",f20.12,/,2x,&
                             &"Initial NSC force: ",f20.12)') E0, F0
    ! Move the specified atom
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    ! Recalculate atomic densities
    call set_density
    ! Note that we've held K fixed but allow potential to vary ? Yes:
    ! this way we get h_on_support in workspace_support
    call get_H_matrix(.false., fixed_potential, electrons, density, &
                      maxngrid)
    call get_energy(total_energy)
    ! Find force
    call get_nonSC_correction_force(nonSC_force, density_out, inode, &
                                    ionode, ni_in_cell, maxngrid)
    E1 = band_energy + delta_E_hartree + delta_E_xc
    F1 = nonSC_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Final energy: ",f20.12,/,2x,"Final &
                             &NSC force: ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = half * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical Force: ",f20.12,/,2x,&
                             &"Analytic Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)

    deallocate(nonSC_force, density_out, STAT=stat)
    if (stat /= 0) call cq_abort("test_nonSC: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 3*ni_in_cell+nspin*maxngrid, type_dbl)

    return
  end subroutine test_nonSC
  !!***

  !!****f* test_force_module/test_full *
  !!
  !!  NAME 
  !!   test_full
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the full force
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   10:02, 05/09/2003 
  !!  MODIFICATION HISTORY
  !!   2007/02/21 17:13 dave
  !!    Minor bug: ewald_CS update onl needed for new ewald
  !!   2011/09/19 L.Tong
  !!    Added Spin polarisation
  !!   2011/11/29 L.Tong
  !!    now uses numbers module
  !!   2011/12/09 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_full(fixed_potential, vary_mu, n_L_iterations, &
                       L_tolerance, tolerance, total_energy,     &
                       expected_reduction)

    use datatypes
    use numbers
    use move_atoms,      only: primary_update, cover_update, &
                               update_atom_coord, update_H
    use group_module,    only: parts
    use cover_module,    only: BCS_parts, DCS_parts, ewald_CS
    use primary_module,  only: bundle
    use global_module,   only: iprint_MD, x_atom_cell, y_atom_cell, &
                               z_atom_cell, id_glob_inv, ni_in_cell
    use energy,          only: local_ps_energy, hartree_energy, &
                               xc_energy, get_energy, band_energy
    use force_module,    only: force, tot_force
    use GenComms,        only: myid, inode, ionode, cq_abort
    use H_matrix_module, only: get_H_matrix
    use global_module,   only: flag_self_consistent, nspin
    use SelfCon,         only: new_SC_potl
    use DMMin,           only: FindMinDM
    use density_module,  only: density
    use ewald_module,    only: flag_old_ewald
    use maxima_module,   only: maxngrid
    use memory_module,   only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical      :: vary_mu, find_chdens, fixed_potential
    logical      :: start, start_L
    integer      :: n_L_iterations
    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer      :: stat
    logical      :: reset_L
    real(double) :: E0, F0, E1, F1, analytic_force, numerical_force
    real(double) :: electrons_tot
    real(double), dimension(nspin) :: electrons
    real(double), dimension(:,:), allocatable :: p_force

    allocate(p_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) call cq_abort("test_full: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    ! We're coming in from initial_H: assume that initial E found
    ! Ensure that h_on_support is in workspace_support
    call get_H_matrix(.false., fixed_potential, electrons, density, &
                      maxngrid)
    ! Store local energy
    call get_energy(total_energy)
    E0 = total_energy
    ! Find force
    call force(fixed_potential, vary_mu, n_L_iterations, L_tolerance, &
               tolerance, total_energy, expected_reduction, .true.)
    ! Find out direction and atom for displacement
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Moving atom ",i5," in &
                             &direction ",i2," by ",f10.6," bohr")') &
               TF_atom_moved, TF_direction, TF_delta
    F0 = tot_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun, fmt='(2x,"Initial energy: ",&
                              &f20.12,/,2x,"Initial force : ",f20.12)') &
               E0, F0
    ! Move the specified atom 
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                        bundle, parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    if (.not. flag_old_ewald) &
         call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                           ewald_CS, parts)
    ! Reproject blips
    call update_H(fixed_potential)
    if (flag_self_consistent) then ! Vary only DM and charge density
       reset_L = .true.
       call new_SC_potl(.true., tolerance, reset_L, fixed_potential, &
                        vary_mu, n_L_iterations, L_tolerance, total_energy)
    else ! Ab initio TB: vary only DM
       call get_H_matrix(.true., fixed_potential, electrons, density, &
                         maxngrid)
       call FindMinDM(n_L_iterations, vary_mu, L_tolerance, inode, &
                      ionode, reset_L, .false.)
       call get_energy(total_energy)
    end if
    E1 = total_energy
    call force(fixed_potential, vary_mu, n_L_iterations, L_tolerance, &
               tolerance, total_energy, expected_reduction, .true.)
    F1 = tot_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Final energy: ",f20.12,/,2x,"Final &
                             &force : ",f20.12)') E1, F1
    numerical_force = -(E1 - E0) / TF_delta
    analytic_force = half * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical Force: ",f20.12,/,2x,&
                             &"Analytic Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)

    deallocate(p_force, STAT=stat)
    if (stat /= 0) call cq_abort("test_full: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    return
  end subroutine test_full
  !!***


  !!****f* test_force_module/test_cdft *
  !!
  !!  NAME 
  !!   test_cdft
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Tests the cdft force
  !!
  !!   
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14:15, 02/09/2003
  !!  MODIFICATION HISTORY
  !!   08:20, 2003/09/03 dave
  !!    Changed to use the change in 2Tr[K.H].
  !!   2011/11/29 L.Tong
  !!    - updated calls to build_Becke_weight_forces
  !!    - removed redundant local variable electrons
  !!    - removed redundant dependence on H_matrix_module
  !!    - removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine test_cdft(fixed_potential, vary_mu, n_L_iterations, &
                       L_tolerance, tolerance, total_energy,     &
                       expected_reduction)

    use datatypes
    use move_atoms,      only: primary_update, cover_update,        &
                               update_atom_coord
    use group_module,    only: parts
    use cover_module,    only: BCS_parts, DCS_parts
    use primary_module,  only:  bundle
    use global_module,   only: iprint_MD, x_atom_cell, y_atom_cell, &
                               z_atom_cell, id_glob_inv,            &
                               flag_self_consistent, ni_in_cell
    use energy,          only: local_ps_energy, hartree_energy,     &
                               xc_energy, get_energy, band_energy,  &
                               cdft_energy
    use GenComms,        only: myid, inode, ionode, cq_abort
    use S_matrix_module, only: get_S_matrix
    use global_module,   only: WhichPulay, SPulay, flag_basis_set,  &
                               blips, PAOs
    use SelfCon,         only: new_SC_potl
    use DMMin,           only: FindMinDM
    use density_module,  only: build_Becke_weight_forces,           &
                               build_Becke_weights
    use maxima_module,   only: maxngrid
    use cdft_module,     only: make_weights, get_cdft_constraint
    use memory_module,   only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: vary_mu, find_chdens, fixed_potential
    logical :: start, start_L
    integer :: n_L_iterations

    real(double) :: tolerance, L_tolerance
    real(double) :: expected_reduction
    real(double) :: total_energy

    ! Local variables
    integer :: stat
    logical :: reset_L
    real(double) ::  E0, F0, E1, F1, analytic_force, numerical_force
    real(double), dimension(:,:), allocatable :: c_force

    allocate(c_force(3,ni_in_cell), STAT=stat)
    if (stat /= 0) call cq_abort("test_cdft: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    ! We're coming in from initial_H: assume that initial E found
    ! Find force
    call get_S_matrix(inode, ionode)
    call get_energy(total_energy)

    call build_Becke_weight_forces(c_force)

    E0 = cdft_energy
    ! Find out direction and atom for displacement
    write (io_lun, fmt='(2x,"Moving atom ",i5," in direction ",i2," by &
                         &",f10.6," bohr")') &
          TF_atom_moved, TF_direction,TF_delta
    F0 = c_force(TF_direction,TF_atom_moved)
    if (inode == ionode) &
         write (io_lun, fmt='(2x,"Initial cDFT energy: ",f20.12,/,2x,&
                              &"Initial cDFT force: ",f20.12)') E0, F0
    ! Move the specified atom 
    if(TF_direction==1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if(TF_direction==2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    else if(TF_direction==3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) + TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    ! Regenerate W
    call build_Becke_weights
    call make_weights
    call get_cdft_constraint
    call get_energy(total_energy) ! Gets the new constraint energy
    ! before Vc changes
    call build_Becke_weight_forces(c_force)
    E1 = cdft_energy
    F1 = c_force(TF_direction,TF_atom_moved)

    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Final cDFT energy: ",f20.12,/,2x,&
                             &"Final cDFT force: ",f20.12)') E1, F1
    numerical_force = -(E1-E0)/TF_delta
    analytic_force = 0.5_double * (F1 + F0)
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"Numerical cDFT Force: ",f20.12,/,2x,&
                             &"Analytic cDFT Force : ",f20.12)') &
               numerical_force, analytic_force
    if (inode == ionode) &
         write (io_lun,fmt='(2x,"cDFT Force error: ",e20.12)') &
               numerical_force - analytic_force
    ! Move the specified atom back
    if (TF_direction == 1) then
       x_atom_cell(id_glob_inv(TF_atom_moved)) = &
            x_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 2) then
       y_atom_cell(id_glob_inv(TF_atom_moved)) = &
            y_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    else if (TF_direction == 3) then
       z_atom_cell(id_glob_inv(TF_atom_moved)) = &
            z_atom_cell(id_glob_inv(TF_atom_moved)) - TF_delta
    end if
    call update_atom_coord
    ! Update positions and indices
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle,&
                        parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, &
                      DCS_parts, parts)
    call get_S_matrix(inode, ionode) !Builds new weight matrix
    !automatically.

    deallocate(c_force, STAT=stat)
    if (stat /= 0) call cq_abort("test_cdft: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

    return
  end subroutine test_cdft
  !!***

end module test_force_module
