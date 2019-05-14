! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module force_module
! ------------------------------------------------------------------------------
! Code area 7: moving atoms
! ------------------------------------------------------------------------------

!!****h* Conquest/force_module *
!!  NAME
!!   force_module
!!  PURPOSE
!!   Gathers together the force routines and variables
!!  USES
!!   atoms, blip_grid_transform_module, blip,
!!   calc_matrix_elements_module, common, datatypes, DiagModule,
!!   dimens, ion_electrostatic, generic_blas, generic_comms, global_module,
!!   grid_index, H_matrix_module, logicals, matrix_data,
!!   maxima_module, multiply_module, mult_module, numbers,
!!   primary_module, pseudopotential_data, set_bucket_module,
!!   S_matrix_module, species_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   12/06/2001 dave
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Included matrix_diagonal
!!    Removed use matrix_diag in get_HF_non_local_force
!!    Removed matrix scaling from same
!!   13/05/2002 dave
!!    Added RCS id variable and TODO to get_HF_force - this needs fixing !
!!   29/05/2002 dave
!!    Added flag to check for ordern before getting K, also RCS static
!!    object and tweaked headers
!!   31/07/2002 dave
!!    Changed pulay_force to use data_M12 from matrix_data and not to
!!    pass it to get_support_gradient
!!   15:42, 08/04/2003 drb
!!    Bug fixes (factor of 1/3 in XC, sign in nonSC), added nonSC
!!    forces and differential of GTH XC LDA parameterisation
!!   14:20, 02/09/2003 drb
!!    Added parameters to allow the non-local routine to find HF, phi
!!    Pulay or both.
!!   11:51, 30/09/2003 drb
!!    Mainly changes to correct non-SC forces
!!   10:09, 13/02/2006 drb
!!    Removed all explicit references to data_ variables and rewrote
!!    in terms of new
!!    matrix routines
!!   2008/02/06 08:14 dave
!!    Changed for output to file not stdout
!!   2011/03/31 11:52 M.Arita
!!    Added sbrt get_pcc_force and modified sbrt get_nonSC_correction
!!    for P.C.C.
!!   2011/10/03 08:17 dave
!!    Adding cDFT forces
!!   2012/04/16 L.Tong
!!   - Moved all the dxc_potential subroutines to XC_module
!!   2012/04/03 09:38 dave
!!    Changes for analytic blips
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2015/05/01 13:55 dave and sym
!!    Implementing stress
!!   2015/08/10 08:04 dave
!!    Adding non-SCF and PCC stress components
!!   2015/11/26 15:24 dave
!!    Changing ewald_force and ewald_stress to ion_interaction_force and _stress
!!   2017/11/8 10:44 zamaan
!!    added target attribute to stress
!!   2018/01/24 11:45 JST dave
!!    Added NA integral & projector approach to forces (flag_neutral_atom_projector)
!!   2019/03/28 zamaan
!!    Changed stress 3-vectors to 3x3 matrices. Modified calculations to
!!    include off-diagonal elements.
!!   2019/05/08 zamaan
!!    Added atomic contributions to stress for computing heat flux
!!  SOURCE
!!
module force_module

  use datatypes
  use global_module,          only: io_lun, atomic_stress
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_allocation,      &
                                    tmr_std_matrices

  implicit none

  save

  real(double), dimension(:,:), allocatable, target :: tot_force
  real(double), dimension(:,:), allocatable :: s_pulay_for, phi_pulay_for

  ! On-site part of stress tensor as Conquest uses orthorhombic cells (easily extended)
  real(double), dimension(3,3), target :: stress
  real(double), dimension(3,3) :: SP_stress, KE_stress, NL_stress, &
                                  PP_stress, GPV_stress, XC_stress, &
                                  nonSCF_stress, pcc_stress, NA_stress

  ! Useful parameters for selecting force calculations in NL part
  integer, parameter :: HF = 1
  integer, parameter :: Pulay = 2
  integer, parameter :: HF_and_Pulay = 3

  ! Area identification
  integer, parameter, private :: area = 7

  ! RCS tag for object file identification
  character(len=80), private :: &
       RCSid = "$Id$"

!!***

contains

  ! ------------------------------------------------------------------------------
  ! Subroutine force
  ! ------------------------------------------------------------------------------

  !!****f* force_module/force *
  !!
  !!  NAME
  !!   force
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Collects all force calculations in one place
  !!
  !! ********************
  !! ** IMPORTANT NOTE **
  !! ********************
  !!
  !!   I know it shouldn't be like this, but it is: the order of the
  !!   calls is VITAL ! pulay_force assumes that workspace_support
  !!   contains h_on_atomfns (as generated by get_H_matrix) so should be
  !!   first (or should regenerate h_on_atomfns).  get_KE_force
  !!   overwrites the on-site blocks of data_K with zeroes, and so
  !!   should LAST !
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   Sometime in 1998, I think
  !!  MODIFICATION HISTORY
  !!   15/05/2001 dave
  !!    F90, indenting and ROBODoc
  !!   21/05/2001 dave
  !!    Reduced calls to get_HF_force and get_KE_force
  !!    Reduced call to pulay_force
  !!   22/05/2001 dave
  !!    Bug fixes
  !!    Shortened overall subroutine call and comments
  !!    before subroutine calls
  !!   08/06/2001 dave
  !!    Added GenComms for my_barrier and RCS Id and Log tags
  !!   12/06/2001 dave
  !!    Included in force_module
  !!   13/05/2002 dave
  !!    Tidied format statements
  !!   11:52, 30/09/2003 drb
  !!    Added density_out variable to store output density for non-SC
  !!    calculations and changed call to get_nonSC_correction_force.
  !!    Also added call to get_electronic_density.  Also placed
  !!    pulay_force call
  !!    FIRST (as it needs h_on_support in workspace_support).
  !!   12:00, 31/03/2011 M.Arita
  !!    Added the contributions from P.C.C.
  !!   2011/09/19 L.Tong
  !!    Added Spin Polarisation
  !!    - removed reference to potential of potential_module, this is
  !!      now referenced in get_nonSC_correction_force subroutine
  !!   2011/10/03 08:17 dave
  !!    Calls for cDFT
  !!   2011/11/28 L.Tong
  !!    Removed redundant dependency on number_of_bands
  !!   2012/03/26 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!   2012/04/03 09:39 dave
  !!    Added KE force to pulay force call (for analytic blips)
  !!   2013/07/10 11:29 dave
  !!    Bug fix for sum over two components of rho even without spin
  !!   2015/05/01 13:56 dave and sym
  !!    Adding stress sum and output
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2015/09/03 17:10 dave
  !!    Correcting Jacobian terms for non-SCF stress and displaying
  !!    PCC and non-SCF stresses
  !!   2015/11/24 08:38 dave
  !!    Adjusted name of hartree_energy to hartree_energy_total_rho for neutral atom implementation
  !!   2015/11/26 15:24 dave
  !!    Name change for ewald_force and stress to ion_interaction_force and stress
  !!   2015/12/09 17:31 dave
  !!    Rationalised force sum and output (esp. PCC) and added changes for neutral atom
  !!   2016/01/28 16:43 dave
  !!    Changed to use ion_electrostatic (instead of ewald_module)
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2018/03/06 15:46 dave
  !!    Added output in GPa for pressure (along with conventional sign change)
  !!   2019/05/08 zamaan
  !!    Initialise atomic_stress, the array for storing atomic stress contributions
  !!  SOURCE
  !!
  subroutine force(fixed_potential, vary_mu, n_cg_L_iterations, &
                   tolerance, con_tolerance, total_energy,      &
                   expected_reduction, write_forces, level )

    use datatypes
    use numbers
    use units
    use timer_module
    use ion_electrostatic,      only: ion_interaction_force, ion_interaction_stress, &
                                      screened_ion_force, screened_ion_stress
    use pseudopotential_data,   only: non_local
    use GenComms,               only: my_barrier, inode, ionode,       &
                                      cq_abort, gsum
    ! TM new pseudo
    use pseudopotential_common, only: pseudo_type, OLDPS, SIESTA,      &
                                      STATE, ABINIT, core_correction, flag_neutral_atom_projector
    use pseudo_tm_module,       only: loc_pp_derivative_tm, loc_HF_stress, loc_G_stress
    use global_module,          only: flag_self_consistent,            &
                                      flag_move_atom, id_glob,         &
                                      WhichPulay, BothPulay, PhiPulay, &
                                      SPulay, flag_basis_set, PAOs,    &
                                      atomf, sf,                       &
                                      blips, ni_in_cell, iprint_MD,    &
                                      IPRINT_TIME_THRES2,              &
                                      area_moveatoms, flag_pcc_global, &
                                      flag_perform_cdft, flag_dft_d2,  &
                                      nspin, spin_factor, &
                                      flag_analytic_blip_int, &
                                      flag_neutral_atom, flag_stress, &
                                      rcellx, rcelly, rcellz, &
                                      flag_atomic_stress, non_atomic_stress, &
                                      flag_heat_flux
    use density_module,         only: get_electronic_density, density, &
                                      build_Becke_weight_forces
    use functions_on_grid,      only: atomfns, H_on_atomfns
    use dimens,                 only: n_my_grid_points
    use matrix_data,            only: Hrange
    use mult_module,            only: matK, matKatomf, SF_to_AtomF_transform
    use maxima_module,          only: maxngrid
    use memory_module,          only: reg_alloc_mem, reg_dealloc_mem,  &
                                      type_dbl
    use DFT_D2,                 only: disp_force
    use energy,                  only: hartree_energy_total_rho, local_ps_energy, &
                                       delta_E_xc, xc_energy, hartree_energy_drho
    use hartree_module, only: Hartree_stress
    use XC, ONLY: XC_GGA_stress

    implicit none

    ! Passed variables
    logical      :: vary_mu, fixed_potential, write_forces
    integer      :: n_cg_L_iterations
    real(double) :: tolerance, con_tolerance, total_energy
    real(double) :: expected_reduction
    integer, optional :: level 

    ! Local variables
    integer        :: i, j, ii, stat, max_atom, max_compt, ispin, &
                      direction, dir1, dir2
    real(double)   :: max_force, volume, scale
    type(cq_timer) :: tmr_l_tmp1
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

    real(double), dimension(nspin)            :: electrons
    real(double), dimension(:),   allocatable :: density_total
    real(double), dimension(:,:), allocatable :: p_force,         &
                                                 cdft_force,      &
                                                 HF_force,        &
                                                 HF_NL_force,     &
                                                 KE_force,        &
                                                 nonSC_force,     &
                                                 pcc_force, NA_force
    real(double), dimension(:),   allocatable :: density_out_tot
    real(double), dimension(:,:), allocatable :: density_out

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='force',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    call start_timer(tmr_std_allocation)
    allocate(density_total(maxngrid), STAT=stat)
    if (stat /= 0) call cq_abort("force: Error alloc mem: ", maxngrid)
    call reg_alloc_mem(area_moveatoms, maxngrid, type_dbl)
    if (flag_pcc_global) then
       allocate(p_force(3,ni_in_cell), KE_force(3,ni_in_cell),      &
                HF_force(3,ni_in_cell), HF_NL_force(3,ni_in_cell),  &
                nonSC_force(3,ni_in_cell), pcc_force(3,ni_in_cell), &
                STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating forces: ", ni_in_cell)
       call reg_alloc_mem(area_moveatoms, 6 * 3 * ni_in_cell, type_dbl)
    else
       allocate(p_force(3,ni_in_cell), KE_force(3,ni_in_cell),     &
                HF_force(3,ni_in_cell), HF_NL_force(3,ni_in_cell), &
                nonSC_force(3,ni_in_cell), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating forces: ", ni_in_cell)
       call reg_alloc_mem(area_moveatoms, 5 * 3 * ni_in_cell, type_dbl)
    end if
    if (flag_perform_cdft) then
       allocate(cdft_force(3,ni_in_cell), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating cdft_force: ", ni_in_cell)
       call reg_alloc_mem(area_moveatoms, 3 * ni_in_cell, type_dbl)
       cdft_force = zero
    end if
    if(iprint_MD>3) then
       allocate(s_pulay_for(3,ni_in_cell),phi_pulay_for(3,ni_in_cell))
       s_pulay_for = zero
       phi_pulay_for = zero
    end if
    allocate(NA_force(3,ni_in_cell))
    NA_force = zero
    call stop_timer (tmr_std_allocation)
    ! get total density
    density_total = zero
    if(nspin==1) then
       density_total(:) = spin_factor * density(:,1)
    else
       density_total    = spin_factor * sum(density, nspin)
    end if
    ! The 'pulay force' is the force due to the change in energy caused
    ! by the change in the basis functions as the atoms move.
    p_force     = zero
    KE_force    = zero
    HF_force    = zero
    HF_NL_force = zero
    nonSC_force = zero
    if (flag_pcc_global) pcc_force = zero ! for P.C.C.
    tot_force   = zero
    ! Zero stresses
    stress  = zero
    KE_stress = zero
    SP_stress = zero
    PP_stress = zero
    NL_stress = zero
    GPV_stress = zero
    XC_stress = zero
    pcc_stress = zero
    nonSCF_stress = zero
    if(flag_neutral_atom_projector) NA_stress = zero
    ! Probably wrong to call this GPV; it's really a jacobian term, from the change in the
    ! integration volume element for grid-based integrals
    ! Different definitions for non-SCF and SCF
    if (flag_stress) then
      do dir1 = 1,3
         if(flag_neutral_atom) then
            if(flag_neutral_atom_projector) then
               GPV_stress(dir1,dir1) = hartree_energy_drho ! NA is not done on grid
            else
               GPV_stress(dir1,dir1) = (hartree_energy_drho + local_ps_energy)
            end if
         else
            GPV_stress(dir1,dir1) = (hartree_energy_total_rho + &
              local_ps_energy - core_correction) ! core contains 1/V term
         end if
         if(flag_self_consistent) then
            XC_stress(dir1,dir1) = xc_energy + &
              spin_factor*XC_GGA_stress(dir1,dir1)
         else ! nonSCF XC found later, along with corrections to Hartree
            XC_stress(dir1,dir1) = delta_E_xc !xc_energy + spin_factor*XC_GGA_stress(direction)
         end if
      end do    
    end if
    WhichPulay  = BothPulay

    ! matK->matKatomf backtransformation for contracted SFs
    if (atomf.ne.sf) then
       do ispin = 1, nspin
          call SF_to_AtomF_transform(matK(ispin), matKatomf(ispin), ispin, Hrange)
       enddo
    endif

    ! This ASSUMES that H_on_atomfns contains the values of H
    ! acting on atomic functions (SF or PAO)
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    call pulay_force(p_force, KE_force, fixed_potential, vary_mu, &
                     n_cg_L_iterations, tolerance, con_tolerance, &
                     total_energy, expected_reduction, ni_in_cell)
    call stop_print_timer(tmr_l_tmp1, "Pulay force", IPRINT_TIME_THRES2)
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    if (flag_perform_cdft) then
       call build_Becke_weight_forces(cdft_force)
    end if
    call stop_print_timer(tmr_l_tmp1, "cDFT force", IPRINT_TIME_THRES2)

    ! Different forces depending on whether we're doing Harris-Foulkes
    ! or self-consistent
    if (.not. flag_self_consistent) then

       call start_timer(tmr_std_allocation)
       allocate(density_out(maxngrid,nspin), &
                density_out_tot(maxngrid), STAT=stat)
       if (stat /= 0) &
            call cq_abort ("Error allocating output density: ", stat)
       call reg_alloc_mem(area_moveatoms, maxngrid * (nspin + 1), type_dbl)

       density_out     = zero
       density_out_tot = zero
       call stop_timer(tmr_std_allocation)

       call start_timer(tmr_l_tmp1,WITH_LEVEL)
       call get_electronic_density(density_out, electrons,   &
                                   atomfns, H_on_atomfns(1), &
                                   inode, ionode, maxngrid)
       ! get the total density_out
       if(nspin==1) then
          density_out_tot(:) = spin_factor * density_out(:,1)
       else
          density_out_tot    = spin_factor * sum(density_out, nspin)
       end if
       call stop_print_timer(tmr_l_tmp1, "get_electronic_density", &
                             IPRINT_TIME_THRES2)
       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       ! for P.C.C.
       if (flag_pcc_global) then
          call get_pcc_force(pcc_force, inode, ionode, ni_in_cell, &
                             maxngrid, density_out) ! Pass output density for non-SCF stress
       end if
       
       call get_nonSC_correction_force(nonSC_force, density_out,  &
                                       inode, ionode, ni_in_cell, &
                                       maxngrid)
       call stop_print_timer(tmr_l_tmp1, "NSC force", IPRINT_TIME_THRES2)

       ! Local HF force
       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       select case (pseudo_type)
       case (OLDPS)
          call get_HF_force(hf_force, density_out_tot, ni_in_cell, maxngrid)
       case (SIESTA)
          call loc_pp_derivative_tm(hf_force, density_out_tot, maxngrid)
       case (ABINIT)
          call loc_pp_derivative_tm(hf_force, density_out_tot, maxngrid)
       end select
       call stop_print_timer(tmr_l_tmp1, "local pseudopotential force", &
                             IPRINT_TIME_THRES2)

       ! deallocate the temporary arrays
       call start_timer(tmr_std_allocation)
       deallocate(density_out_tot, density_out, STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error deallocating output density: ", &
                          maxngrid)
       call reg_dealloc_mem(area_moveatoms, maxngrid * (nspin + 1), type_dbl)
       call stop_timer(tmr_std_allocation)

    else ! for SCF

       ! Local HF force
       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       ! for P.C.C.
       if (flag_pcc_global) then
          call get_pcc_force(pcc_force, inode, ionode, ni_in_cell, maxngrid)
       end if
       select case (pseudo_type)
       case (OLDPS)
          call get_HF_force(hf_force, density_total, ni_in_cell, maxngrid)
       case (SIESTA)
          call loc_pp_derivative_tm(hf_force, density_total, maxngrid)
       case (ABINIT)
          call loc_pp_derivative_tm(hf_force, density_total, maxngrid)
       end select
       call stop_print_timer(tmr_l_tmp1, "local pseudopotential force", &
                             IPRINT_TIME_THRES2)
    end if ! (.not. flag_self_consistent)

    ! This routine deals with the movement of the nonlocal pseudopotential.
    if (non_local) then
       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       call get_HF_non_local_force(HF_NL_force, HF_and_Pulay, ni_in_cell)
       call stop_print_timer(tmr_l_tmp1, "nonlocal pseudopotential force", &
                             IPRINT_TIME_THRES2)
    else
       HF_NL_force = zero
    end if
    ! Get the kinetic energy force component
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if(flag_basis_set==PAOs.OR.(flag_basis_set==blips.AND.(.NOT.flag_analytic_blip_int))) call get_KE_force(KE_force, ni_in_cell)
    call stop_print_timer(tmr_l_tmp1, "kinetic energy force", &
                          IPRINT_TIME_THRES2)
    if(flag_neutral_atom_projector) call get_HNA_force(NA_force)
    max_force = zero
    max_atom  = 0
    max_compt = 0
    if (inode == ionode .and. write_forces) then
       write (io_lun, fmt='(/,20x,"Forces on atoms (",a2,"/",a2,")"/)') &
             en_units(energy_units), d_units(dist_units)
       write (io_lun, fmt='(18x,"    Atom   X              Y              Z")')
    end if
    ! Calculate forces and write out
    do i = 1, ni_in_cell
       do j = 1, 3
          ! Force components that are always needed
          tot_force(j,i) = HF_force(j,i) + HF_NL_force(j,i) + &
               p_force(j,i) +    KE_force(j,i)
          ! Non-self-consistent
          if (.not.flag_self_consistent) &
               tot_force(j,i) = tot_force(j,i) + nonSC_force(j,i)
          ! Neutral atom or conventional pseudopotential
          if(flag_neutral_atom) then
             tot_force(j,i) = tot_force(j,i) + screened_ion_force(j,i)
          else
             tot_force(j,i) = tot_force(j,i) + ion_interaction_force(j,i)
          end if
          ! PCC
          if (flag_pcc_global) &
               tot_force(j,i) = tot_force(j,i) + pcc_force(j,i)
          ! CDFT
          if (flag_perform_cdft) &
               tot_force(j,i) = tot_force(j,i) + cdft_force(j,i)
          ! DFT-D2
          if (flag_dft_d2) &
               tot_force(j,i) = tot_force(j,i) + disp_force(j,i)
          if(flag_neutral_atom_projector) &
               tot_force(j,i) = tot_force(j,i) + NA_force(j,i)
          ! Zero force on fixed atoms
          if (.not. flag_move_atom(j,i)) then
             tot_force(j,i) = zero
          end if
          if (abs (tot_force(j,i)) > max_force) then
             max_force = abs (tot_force(j,i))
             max_atom  = i
             max_compt = j
          end if
       end do ! j
       if (inode == ionode) then
          if(iprint_MD > 2) then
             write(io_lun, 101) i
             write(io_lun, 102) (for_conv *   HF_force(j,i),  j = 1, 3)
             if(flag_neutral_atom_projector) write (io_lun, fmt='("Force NA     : ",3f15.10)') (for_conv*NA_force(j,i),j=1,3)
             write(io_lun, 112) (for_conv * HF_NL_force(j,i), j = 1, 3)
             write(io_lun, 103) (for_conv *     p_force(j,i), j = 1, 3)
             if(iprint_MD>3) then
                write (io_lun, fmt='("  Phi pulay  : ",3f15.10)') (for_conv*phi_pulay_for(j,i),j=1,3)
                write (io_lun, fmt='("  S pulay    : ",3f15.10)') (for_conv*s_pulay_for(j,i),j=1,3)
             end if
             write(io_lun, 104) (for_conv *    KE_force(j,i), j = 1, 3)
             if(flag_neutral_atom) then
                write(io_lun, 106) (for_conv * screened_ion_force(j,i), j = 1, 3)
             else
                write(io_lun, 106) (for_conv * ion_interaction_force(j,i), j = 1, 3)
             end if
             if (flag_pcc_global) write(io_lun, 108) (for_conv *   pcc_force(j,i), j = 1, 3)
             if (flag_dft_d2) write (io_lun, 109) (for_conv * disp_force(j,i), j = 1, 3)
             if (flag_perform_cdft) write (io_lun, fmt='("Force cDFT : ",3f15.10)') &
                  (for_conv*cdft_force(j,i),j=1,3)
             if (flag_self_consistent) then
                write (io_lun, 105) (for_conv * tot_force(j,i),   j = 1, 3)
             else
                write (io_lun, 107) (for_conv * nonSC_force(j,i), j = 1, 3)
                write (io_lun, 105) (for_conv * tot_force(j,i),   j = 1, 3)
             end if
          else if (write_forces) then
             write (io_lun,fmt='(20x,i6,3f15.10)') &
                  i, (for_conv * tot_force(j,i), j = 1, 3)
          end if ! (iprint_MD > 2)
       end if ! (inode == ionode)
    end do ! i
    if (inode == ionode) &
         write (io_lun,                                      &
                fmt='(4x,"Maximum force : ",f15.8,"(",a2,"/",&
                      &a2,") on atom, component ",2i9)')     &
               for_conv * max_force, en_units(energy_units), &
               d_units(dist_units), max_atom, max_compt
    ! We will add PCC and nonSCF stresses even if the flags are not set, as they are
    ! zeroed at the start
    if (flag_stress) then
      if(flag_neutral_atom) then
         do dir1 = 1, 3
            do dir2 = 1,3
               stress(dir1,dir2) = KE_stress(dir1,dir2) + &
                 SP_stress(dir1,dir2) + PP_stress(dir1,dir2) + &
                 NL_stress(dir1,dir2) + GPV_stress(dir1,dir2) + &
                 XC_stress(dir1,dir2) + screened_ion_stress(dir1,dir2) + &
                 Hartree_stress(dir1,dir2) + loc_HF_stress(dir1,dir2) +  &
                 pcc_stress(dir1,dir2) + nonSCF_stress(dir1,dir2)
               if (flag_atomic_stress) then
                 non_atomic_stress(dir1,dir2) = &
                   non_atomic_stress(dir1,dir2) + GPV_stress(dir1,dir2) + &
                   XC_stress(dir1,dir2) + Hartree_stress(dir1,dir2)
               end if
               if (flag_neutral_atom_projector) then
                 stress(dir1,dir2) = stress(dir1,dir2) + NA_stress(dir1,dir2)
               end if

            end do
         end do
      else
         do dir1 = 1, 3
            do dir2 = 1, 3
               stress(dir1,dir2) = KE_stress(dir1,dir2) + &
                 SP_stress(dir1,dir2) + PP_stress(dir1,dir2) + &
                 NL_stress(dir1,dir2) + GPV_stress(dir1,dir2) + &
                 XC_stress(dir1,dir2) + Ion_Interaction_stress(dir1,dir2) + &
                 Hartree_stress(dir1,dir2) + loc_HF_stress(dir1,dir2) + &
                 loc_G_stress(dir1,dir2) + pcc_stress(dir1,dir2) + &
                 nonSCF_stress(dir1,dir2)
               if (flag_atomic_stress) then
                 non_atomic_stress(dir1,dir2) = &
                   non_atomic_stress(dir1,dir2) + GPV_stress(dir1,dir2) + &
                   XC_stress(dir1,dir2) + Hartree_stress(dir1,dir2) + &
                   loc_G_stress(dir1,dir2)
               end if
            end do
         end do
      end if
      if (flag_atomic_stress) call gsum(atomic_stress,3,3,ni_in_cell)

      if (inode == ionode) then       
         write (io_lun,fmt='(/4x,"                  ",3a15)') "X","Y","Z"
         write(io_lun,fmt='(4x,"Stress contributions:")')
      end if
      call print_stress("K.E. stress:      ", KE_stress, 3)
      call print_stress("S-Pulay stress:   ", SP_stress, 3)
      call print_stress("Phi-Pulay stress: ", PP_stress, 3)
      call print_stress("Local stress:     ", loc_HF_stress, 3)
      call print_stress("Local G stress:   ", loc_G_stress, 3)
      call print_stress("Non-local stress: ", NL_stress, 3)
      call print_stress("Jacobian stress:  ", GPV_stress, 3)
      call print_stress("XC stress:        ", XC_stress, 3)
      if (flag_neutral_atom) then
        call print_stress("Ion-ion stress:   ", screened_ion_stress, 3)
        call print_stress("N.A. stress:      ", NA_stress, 3)
      else
        call print_stress("Ion-ion stress:   ", ion_interaction_stress, 3)
      end if
      call print_stress("Hartree stress:   ", Hartree_stress, 3)
      call print_stress("PCC stress:       ", pcc_stress, 3)
      call print_stress("non-SCF stress:   ", nonSCF_stress, 3)
      call print_stress("Total stress:     ", stress, 0)
      volume = rcellx*rcelly*rcellz
      ! Include Ha/cubic bohr to GPa conversion and 1/volume factor
      ! Factor of 1e21 comes from Ang to m (1e30) and Pa to GPa (1e-9) 
      scale = -(HaToeV*eVToJ*1e21_double)/(volume*BohrToAng*BohrToAng*BohrToAng)
      call print_stress("Total pressure:   ", stress*scale, 0)
      if (flag_atomic_stress .and. iprint_MD > 2) call check_atomic_stress
    end if

    call my_barrier()
    if(iprint_MD>3) deallocate(s_pulay_for,phi_pulay_for)
    if (inode == ionode .and. iprint_MD > 1 .and. write_forces) &
         write (io_lun, fmt='(4x,"Finished force")')

    call start_timer(tmr_std_allocation)
    if (flag_pcc_global) then
       deallocate(p_force, KE_force, HF_force, HF_NL_force, &
                  nonSC_force, pcc_force, STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error deallocating forces: ", ni_in_cell)
       call reg_dealloc_mem(area_moveatoms, 6 * 3 * ni_in_cell, type_dbl)
    else
       deallocate(p_force, KE_force, HF_force, HF_NL_force, &
                  nonSC_force, STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error deallocating forces: ", ni_in_cell)
       call reg_dealloc_mem(area_moveatoms, 5 * 3 * ni_in_cell, type_dbl)
    end if
    deallocate(density_total, STAT=stat)
    if (stat /= 0) call cq_abort("force: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, maxngrid, type_dbl)
    call stop_timer(tmr_std_allocation)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='force',echo=.true.)
!****lat>$

    return

101 format('Force on atom ',i9)
102 format('Force H-F    : ',3f15.10)
112 format('Force H-Fnl  : ',3f15.10)
103 format('Force pulay  : ',3f15.10)
104 format('Force KE     : ',3f15.10)
106 format('Force Ion-Ion: ',3f15.10)
107 format('Force nonSC  : ',3f15.10)
108 format('Force PCC    : ',3f15.10)
105 format('Force Total  : ',3f15.10)
109 format('Force disp   : ',3f15.10)

  end subroutine force
  !!***


  ! -----------------------------------------------------------
  ! Subroutine pulay_force
  ! -----------------------------------------------------------

  !!****f* force_module/pulay_force *
  !!
  !!  NAME
  !!   pulay_force
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Evaluates the Pulay contribution to the forces on
  !!   the atoms. This is the change in energy due to
  !!   the movement of the blip functions causing a change in support
  !!   function. The change in energy wrt blip functions is zero (variational)
  !!   but the change in energy wrt the variation of a given support function
  !!   at a given point is not in general zero. If the atom moves, all the
  !!   blip functions associated with it move, and there is a corresponding
  !!   change in the values of the support functions. The gradient of the
  !!   support functions wrt movement of the atoms with which they are
  !!   associated is held in gradsupport. This need to be integrated with
  !!   the gradient of energy wrt change in support functions, which is
  !!   given by the subroutine get_support_gradient.
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   C.M.Goringe
  !!  CREATION DATE
  !!   13/09/95
  !!  MODIFICATION HISTORY
  !!   24/4/96 Chris Goringe
  !!    We have changed the data format, and these changes are required to
  !!    deal with it. blip_to_gradsupport is now replaced by the general
  !!    blip_transform routine
  !!
  !!   24/2/97 Chris Goringe
  !!    Quite a lot of other changes have occured in the code, in particular
  !!    we no longer ever get the grad of the support functions; we just
  !!    get it one component at a time.
  !!
  !!    workspace_support  is used to hold the gradient of energy wrt
  !!                       support functions
  !!    workspace2_support is used to hold gradient of support functions
  !!                       with respect to atomi position
  !!   17/5/97 DRB to implement HeadGordon
  !!   02/02/2001 TM
  !!    blip_to_support_new, get_matrix_elements_new, blip_to_grad_new
  !!   16/05/2001 dave
  !!    F90 format, ROBODoc, new blip_to_grad call
  !!   18/05/2001 dave
  !!    Fixed the get_H_matrix call
  !!   22/05/2001 dave
  !!    Stripped subroutine call
  !!   23/05/2001 dave
  !!    Various small fixes
  !!   12/06/2001 dave
  !!    Included in force_module
  !!   29/05/2002 dave
  !!    Added flag from DiagModule so that K isn't made if we're diagonalising
  !!   31/07/2002 dave
  !!    Changed to use data_M12 from matrix_data and changed call to
  !!    get_support_gradient to remove data_M12
  !!   14:49, 2003/06/09 dave
  !!    Added use statement for blip_gradient (get_support_gradient moved)
  !!   08:06, 2003/10/01 dave
  !!    Added basic structure for PAO basis set
  !!   2011/06/15 13:04 dave
  !!    Added definition of WhichPulay before get_support_gradient call as check
  !!   2011/09/19 L.Tong
  !!    Added spin polarisation
  !!   2011/11/28 L.Tong
  !!    Removed redundant dependency on number_of_bands
  !!   2012/03/25 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!   2012/04/03 09:41 dave
  !!    Added code for analytic blip integrals
  !!   2013/04/24 15:06 dave
  !!    Updates to correct analytic blip forces
  !!   2015/05/08 08:26 dave and sym
  !!    Adding local phi Pulay and S-Pulay stress calculations
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_sf_rem -> atomf_atomf_rem
  !!   2016/07/29 18:30 nakata
  !!    Renamed supports_on_atom -> blips_on_atom
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2016/11/09 21:30 nakata
  !!    Introduce atomf-based calculations
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2017/03/28 dave
  !!    Small bug fix on order of deallocation of matM12atomf
  !!   2019/04/04 zamaan
  !!    Added calculations for off-diagonal elements of PP_stress and SP_stress
  !!   2019/05/08 zamaan
  !!    Added atomic stress contributions
  !!  SOURCE
  !!
  subroutine pulay_force(p_force, KE_force, fixed_potential, vary_mu,  &
                         n_cg_L_iterations, L_tol, self_tol, &
                         total_energy, expected_reduction, n_atoms, level)

    use datatypes
    use logicals
    use numbers
    use primary_module,              only: bundle
    use matrix_module,               only: matrix, matrix_halo
    use matrix_data,                 only: mat, halo, blip_trans, Srange, aSa_range
    use mult_module,                 only: LNV_matrix_multiply,      &
                                           matM12,                   &
                                           allocate_temp_matrix,     &
                                           free_temp_matrix,         &
                                           return_matrix_value,      &
                                           matrix_pos, ltrans,       &
                                           scale_matrix_value,       &
                                           matKatomf,                &
                                           return_matrix_block_pos,  &
                                           matrix_scale,             &
                                           SF_to_AtomF_transform
    use global_module,               only: iprint_MD, WhichPulay,    &
                                           BothPulay, PhiPulay,      &
                                           SPulay, flag_basis_set,   &
                                           blips, PAOs, sf, atomf,   &
                                           flag_onsite_blip_ana,     &
                                           nspin, spin_factor,       &
                                           flag_analytic_blip_int,   &
                                           id_glob, species_glob,    &
                                           flag_diagonalisation,     &
                                           flag_full_stress, flag_stress, &
                                           flag_atomic_stress
    use set_bucket_module,           only: rem_bucket, atomf_atomf_rem
    use blip_grid_transform_module,  only: blip_to_support_new,      &
                                           blip_to_grad_new
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use S_matrix_module,             only: get_S_matrix, get_dS_analytic_oneL, get_r_on_atomfns
    use H_matrix_module,             only: get_H_matrix
    use GenComms,                    only: gsum, cq_abort, mtime,    &
                                           inode, ionode
    !use DiagModule,                  only: diagon
!    use blip_gradient,               only: get_support_gradient
    use PAO_grid_transform_module,   only: single_PAO_to_grad
    use build_PAO_matrices,          only: assemble_deriv_2
    use cover_module,                only: BCS_parts
    use functions_on_grid,           only: atomfns,                  &
                                           H_on_atomfns,             &
                                           allocate_temp_fn_on_grid, &
                                           free_temp_fn_on_grid, gridfunctions
    use species_module,              only: nsf_species, natomf_species
    ! Temp
    use dimens,                      only: grid_point_volume,        &
                                           n_my_grid_points
    use density_module,              only: get_electronic_density,   &
                                           density
    use maxima_module,               only: maxngrid
    use comms_module,                ONLY: start_blip_transfer, fetch_blips
    use group_module,                ONLY: parts
    use blip,                        ONLY: blip_info
    use GenComms,                    ONLY: myid, cq_abort, my_barrier, mtime
    use mpi
    use support_spec_format,         ONLY: support_function, coefficient_array, &
                                           blips_on_atom
    use GenBlas,                     only: axpy, scal
    use calc_matrix_elements_module, only: act_on_vectors_new

    implicit none

    ! Passed variables
    logical      :: vary_mu, fixed_potential
    integer      :: n_atoms
    integer      :: n_cg_L_iterations
    real(double) :: L_tol, self_tol, total_energy, expected_reduction
    real(double), dimension(3,n_atoms) :: p_force, KE_force
    integer, optional :: level

    ! Local variables
    logical      :: test
    integer      :: dir1, dir2, count, nb, na, nsf1, point, place, i, &
                    neigh, ist, wheremat, jsf, gcspart, tmp_fn, n1,  &
                    n2, this_nsf, spin, tmp_fn2
    integer :: spec, icall, jpart, ind_part, j_in_halo, pb_len, nab
    integer :: j,jseq,specj,i_in_prim,speci, pb_st, nblipsj, nod, nsfj,sends,ierr
    integer :: neigh_global_num, neigh_global_part, neigh_species, neigh_prim, this_nsfi, this_nsfj
    integer, allocatable, dimension(:) :: nreqs

    real(double) :: energy_in, t0, t1
    real(double) :: dx, dy, dz, time0, time1
    real(double), dimension(nspin) :: electrons, energy_tmp
    real(double), allocatable, dimension(:), target :: part_blips
    type(support_function) :: supp_on_j
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    real(double), allocatable, dimension(:,:,:) :: this_data_K
    real(double), allocatable, dimension(:,:,:) :: this_data_M12
    real(double) :: forS, forKE
    real(double) :: thisG_dS_dR, thisK_dH_dR
    real(double), dimension(3) :: dr, r_str

    ! Workarray for New Version
    !    I am not sure the present way for calculating p_force is the right way.
    !   Now, we calculate matrix elements data_tmp (= Sum_l phi_ialpha(rl)
    !   * psi_jbeta(rl)), whose off-diagonal elements are not needed.
    !   Then, calculate p_force from summing up diagonal elements of data_tmp.
    !    In the old version, every processor has to communicate with all other
    !   processors. The communication in the new version is local, but
    !   it calculates unneeded off-diagonal elements. Similar situation
    !   occurs in the normalisation of support functions at (initial_phis.f90)
    !    I am planning to implement the subroutine which only calculates
    !   the diagonal elements of matrix (=Sum_l phi_ialpha(rl) * psi_ialpha(rl))
    !   with local communication. (I am not sure it is important or not)
    !    --  02/Feb/2001  Tsuyoshi Miyazaki

    integer        :: iprim, np, ni, isf, mat_tmp, mat_tmp2, n_stress, iatom
    integer, dimension(nspin) :: matM12atomf
    real(double)   :: matM12_value, matK_value
    type(cq_timer) :: tmr_std_loc
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='pulay_force', &
         where=area,level=backtrace_level,echo=.true.)
!****lat>$    

    ! the force due to the change in T matrix elements is done differently...
    if (inode == ionode .and. iprint_MD > 2) &
         write (io_lun, fmt='(4x,"Starting pulay_force()")')
    !p_force = zero

    ! first, lets make sure we have everything up to date.
    ! Update H will give us support, S, K, H
    ! and h(local)_on_support [in workspace_support]

    ! This should be checked and changed: we don't want to do all this
    ! again if we can avoid it !  Probably we will (a) have the data
    ! already, or can (b) try get_SC_potl/get_new_rho or even just
    ! get_H_matrix...
    test = .false.
    if (test) then
       ! (1) get S matrix
       call get_S_matrix(inode, ionode, build_AtomF_matrix=.false.)
       ! (2) get K matrix
       if (.not. flag_diagonalisation) then
          call LNV_matrix_multiply(electrons, energy_tmp, doK, dontM1,&
                                   dontM2, dontM3, dontM4, dontphi, dontE)
       end if
       call get_electronic_density(density, electrons, atomfns,    &
                                   H_on_atomfns(1), inode, ionode, &
                                   maxngrid)
       ! H_on_atomfns is given the correct values after calls of
       ! get_H_matrix
       call get_H_matrix(.false., fixed_potential, electrons, density,&
                         maxngrid)
       if (.not. flag_diagonalisation) then
          call LNV_matrix_multiply(electrons, energy_tmp, dontK, doM1,&
                                   doM2, dontM3, dontM4, dontphi,     &
                                   dontE, mat_M12=matM12)
          ! note both Ne and E are not calculated here (dontphi
          ! and dontE)
       end if
    end if ! if (test)

    t0 = mtime()
    ! for the energy wrt support function we need  M1 and M2...
    ! If we're diagonalising, we've already build data_M12
    if (.not. flag_diagonalisation) then
       call LNV_matrix_multiply(electrons, energy_tmp, dontK, doM1,   &
                                doM2, dontM3, dontM4, dontphi, dontE, &
                                mat_M12=matM12)
    end if
    t1 = mtime()
    t0 = t1

    ! N.B. IT IS VITAL TO HAVE H_on_atomfns to storing the H (and
    ! corresponding spin down component) acting on support functions!

    ! Analytic blip KE and S-pulay force and stress; otherwise zero on-site M12 for on-site analytic blips
    if(flag_basis_set == blips .and. flag_analytic_blip_int) then 
       allocate(nreqs(blip_trans%npart_send))
       ! For speed, we should have a blip_trans%max_len and max_nsf and allocate once
       time0 = mtime()
       call start_blip_transfer(nreqs,sends,parts%mx_ngonn)
       time1 = mtime()
       do jpart=1,halo(aSa_range)%np_in_halo
          pb_len = blip_trans%len_recv(jpart)
          ind_part = halo(aSa_range)%lab_hcell(jpart)
          nod = parts%i_cc2node(ind_part)
          ! Fetch remote blip coefficients for partition
          ! or copy local blip coefficients
          time0 = mtime()
          if(jpart>1) then
             if(ind_part/=halo(aSa_range)%lab_hcell(jpart-1)) then
                allocate(part_blips(pb_len))
                part_blips = zero
                if(nod==myid+1) then
                   pb_st = blip_trans%partst(parts%i_cc2seq(ind_part))
                   part_blips(1:pb_len) = coefficient_array(pb_st:pb_st+pb_len-1)
                else
                   call fetch_blips(part_blips,pb_len,nod-1,(myid)*parts%mx_ngonn + parts%i_cc2seq(ind_part))
                end if
             end if
          else
             allocate(part_blips(pb_len))
             part_blips = zero
             if(nod==myid+1) then
                pb_st = blip_trans%partst(parts%i_cc2seq(ind_part))
                part_blips(1:pb_len) = coefficient_array(pb_st:pb_st+pb_len-1)
             else
                call fetch_blips(part_blips,pb_len,nod-1,(myid)*parts%mx_ngonn + parts%i_cc2seq(ind_part))
             end if
          endif
          time1 = mtime()
          gcspart = halo(aSa_range)%i_hbeg(halo(aSa_range)%lab_hcover(jpart))
          pb_st = 1
          time0 = mtime()
          do j=1,halo(aSa_range)%nh_part(jpart) ! Loop over atoms j in partition
             j_in_halo = halo(aSa_range)%j_beg(jpart)+j-1
             jseq = halo(aSa_range)%j_seq(j_in_halo)
             specj = species_glob( id_glob( parts%icell_beg(halo(aSa_range)%lab_hcell(jpart))+jseq-1) )
             nblipsj = blip_info(specj)%NBlipsRegion
             this_nsfj = nsf_species(specj)
             allocate(supp_on_j%supp_func(this_nsfj))
             do nsfj=1,this_nsfj
                supp_on_j%supp_func(nsfj)%ncoeffs = nblipsj
                supp_on_j%supp_func(nsfj)%coefficients => part_blips(pb_st:pb_st+nblipsj-1)
                pb_st = pb_st+nblipsj
             end do
             do i=1,ltrans(aSa_range)%n_hnab(j_in_halo) ! Loop over atoms i: primary set neighbours of j
                i_in_prim=ltrans(aSa_range)%i_prim(ltrans(aSa_range)%i_beg(j_in_halo)+i-1)
                speci = bundle%species(i_in_prim)
                iatom = bundle%ig_prim(i_in_prim)
                this_nsfi = nsf_species(speci)
                dr(1)=BCS_parts%xcover(gcspart+jseq-1)-bundle%xprim(i_in_prim)
                dr(2)=BCS_parts%ycover(gcspart+jseq-1)-bundle%yprim(i_in_prim)
                dr(3)=BCS_parts%zcover(gcspart+jseq-1)-bundle%zprim(i_in_prim)
                if((dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))>RD_ERR) then
                   allocate(this_data_K(this_nsfi,this_nsfj,nspin),this_data_M12(this_nsfi,this_nsfj,nspin))
                   this_data_K = zero
                   this_data_M12 = zero
                   do spin=1,nspin
                      wheremat = matrix_pos(matKatomf(spin),i_in_prim,j_in_halo,1,1)
                      call return_matrix_block_pos(matKatomf(spin),wheremat,this_data_K(:,:,spin),this_nsfi*this_nsfj)
                      wheremat = matrix_pos(matM12(spin),i_in_prim,j_in_halo,1,1)
                      call return_matrix_block_pos(matM12(spin),wheremat,this_data_M12(:,:,spin),this_nsfi*this_nsfj)
                   end do
                   do dir1=1,3
                      call get_dS_analytic_oneL(blips_on_atom(i_in_prim),supp_on_j, &
                           forS,forKE,this_data_M12,this_data_K, i_in_prim, &
                           j_in_halo,dr(1),dr(2),dr(3),speci,specj, &
                           this_nsfi,this_nsfj,dir1)
                      p_force(dir1,iatom) = p_force(dir1,iatom) - forS!(dir1)
                      KE_force(dir1,iatom) = KE_force(dir1,iatom) - forKE!(dir1)
                      if (flag_stress) then
                        if (flag_full_stress) then
                          do dir2=1,3
                            SP_stress(dir1,dir2) = SP_stress(dir1,dir2) - &
                              forS*dr(dir2)
                            KE_stress(dir1,dir2) = KE_stress(dir1,dir2) - &
                              forKE*dr(dir2)
                            if (flag_atomic_stress) then
                              atomic_stress(dir1,dir2,iatom) = &
                                atomic_stress(dir1,dir2,iatom) - &
                                (forS*dr(dir2) + forKE*dr(dir2))*half
                            end if
                          end do
                        else
                          SP_stress(dir1,dir1) = SP_stress(dir1,dir1) - &
                            forS*dr(dir1)
                          KE_stress(dir1,dir1) = KE_stress(dir1,dir1) - &
                            forKE*dr(dir1)
                        end if ! flag_full_stress
                      end if ! flag_stress
                   end do
                   deallocate(this_data_M12,this_data_K)
                end if
             end do
             do nsfj=1,this_nsfj
                nullify(supp_on_j%supp_func(nsfj)%coefficients)
             end do
             deallocate(supp_on_j%supp_func)
          end do
          if(jpart<halo(aSa_range)%np_in_halo) then
             if(ind_part/=halo(aSa_range)%lab_hcell(jpart+1)) then
                deallocate(part_blips)
             end if
          else
             deallocate(part_blips)
          end if
          time1 = mtime()
       end do
       if(sends>0) then
          do i=1,sends
             call MPI_Wait(nreqs(i),mpi_stat,ierr)
             if(ierr/=0) call cq_abort("Error waiting for blip send to finish",i)
          end do
       end if
       call my_barrier
       deallocate(nreqs)
       call gsum (KE_force, 3, n_atoms)
       if (flag_stress) call gsum(KE_stress, 3, 3)
       if(WhichPulay==PhiPulay) p_force = zero
       if(WhichPulay==BothPulay) WhichPulay = PhiPulay ! We've DONE S-pulay above
       if(WhichPulay==SPulay) then
          !  In principle, the summation below is not needed.
          !  p_force should be calculated only for my primary set of atoms.
          call gsum(p_force, 3, n_atoms)
          return
       end if
    else ! Zero on-site terms if necessary
       if (flag_basis_set == blips .and. flag_onsite_blip_ana) then
          iprim = 0
          do np = 1, bundle%groups_on_node
             do ni = 1, bundle%nm_nodgroup(np)
                iprim = iprim+1
                this_nsf = nsf_species(bundle%species(iprim))
                do n1 = 1, this_nsf
                   do n2 = 1, this_nsf
                      do spin = 1, nspin
                         call scale_matrix_value(matM12(spin), np, ni, &
                              iprim, 0, n1, n2, zero, 1)
                      end do ! spin
                   end do ! n2
                end do ! n1
             end do ! ni
          end do ! np
       end if
       !WhichPulay = BothPulay
    end if
    ! Originally, we evaluated support gradient, but this is only
    ! applicable for blips; moreover it is not appropriate for stresses
    ! so I have commented it out (DRB 2015/05/13 08:13)  We can restore
    ! it if we introduce a no-stress-calculation flag; if so, we may need
    ! to add a clause (basis_set==blips.AND.WhichPulay==SPulay) to the
    ! initial if statement
    !
    ! I would tend to prefer NOT to restore the call to get_support_gradient
    ! but instead to calculate the full gradient here; it seems poor practice
    ! to mix a call to a BLIP-specific gradient into a PAO/blip force routine
    !%%!  This calculates K|H phi> (PAOs) or G|phi> + K|H phi> (blips)
    !%%! call get_support_gradient(H_on_atomfns(1), inode, ionode)
    !%%! t1 = mtime()
    !%%! if (inode == ionode .and. iprint_MD > 3) then
    !%%!    write (io_lun, fmt='(4x,"get_support_gradient time: ",f12.5)') &
    !%%!          t1 - t0
    !%%! end if
    !%%! t0 = t1
    ! NB we could combine the loops below with the S-pulay calculation very easily
    ! ---------------------------------------
    ! Calculate local phi Pulay forces (<grad \phi_i|K_ij|H\phi_j>)
    if (WhichPulay == BothPulay .or. WhichPulay == PhiPulay) then
       mat_tmp = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
       tmp_fn = allocate_temp_fn_on_grid(atomf)
       gridfunctions(tmp_fn)%griddata = zero
       ! Act on H\phi_j with K and store result in H_on_atomfns(1) to save memory
       ! act_on_vectors accumulates
       do spin = 1, nspin
          call act_on_vectors_new(inode-1, rem_bucket(3), matKatomf(spin), &
               tmp_fn, H_on_atomfns(spin))
       end do
       call scal(gridfunctions(tmp_fn)%size, minus_two, &
            gridfunctions(tmp_fn)%griddata, 1)
       gridfunctions(H_on_atomfns(1))%griddata = zero
       call axpy(gridfunctions(H_on_atomfns(1))%size, spin_factor, &
            gridfunctions(tmp_fn)%griddata, 1,           &
            gridfunctions(H_on_atomfns(1))%griddata, 1)
       gridfunctions(tmp_fn)%griddata = zero
       t1 = mtime()
       if (inode == ionode .and. iprint_MD > 3) then
          write (io_lun, fmt='(4x,"get_support_gradient time: ",f12.5)') &
               t1 - t0
       end if
       t0 = t1
       ! Allocate more temporary variables
       tmp_fn2 = allocate_temp_fn_on_grid(atomf)
       gridfunctions(tmp_fn2)%griddata = zero
       mat_tmp2 = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
       !temp_local = allocate_temp_fn_on_grid(atomf) ! Edited SYM 2014/09/01 17:55 
       ! now for each direction in turn
       do dir1 = 1, 3
          ! we get the grad of the support functions
          if (flag_basis_set == blips) then
             call blip_to_grad_new(inode-1, dir1, tmp_fn)
          else if (flag_basis_set == PAOs) then
             call single_PAO_to_grad(dir1, tmp_fn)
          else
             call cq_abort("pulay_force: basis set undefined ", flag_basis_set)
          end if
          t1 = mtime()
          if (inode == ionode .and. iprint_MD > 3)&
               write (io_lun, fmt='(10x,"Phi Pulay grad ",i4," time: ",f12.5)')&
                     dir1, t1 - t0
          t0 = t1
          ! Now scale the gradient by r
          t1 = mtime()
          if (inode == ionode .and. iprint_MD > 3)&
               write (io_lun, fmt='(10x,"Phi Pulay r_on_supp ",i4," time: ",f12.5)')&
                     dir1, t1 - t0
          t0 = t1
          ! Make matrix elements
          call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
                                       mat_tmp, H_on_atomfns(1), tmp_fn)

          n_stress = 1 ! 1 direction for diagonal elements only
          ! Loop over 3 directions for full stress tensor
          if (flag_stress .and. flag_full_stress) n_stress = 3 
          do dir2=1,n_stress
            if (flag_stress) then
              if (flag_full_stress) then
                call get_r_on_atomfns(dir2,tmp_fn,tmp_fn2)
              else 
                call get_r_on_atomfns(dir1,tmp_fn,tmp_fn2)
              end if
            end if
            call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
                                         mat_tmp2, H_on_atomfns(1), tmp_fn2)
            t1 = mtime()
            if (inode == ionode .and. iprint_MD > 3) &
                 write (io_lun, fmt='(10x,"Phi Pulay int ",i4," time: ",f12.5)')&
                       dir1, t1 - t0
            t0 = t1

            ! now store calculated values to p_force
            iprim = 0
            call start_timer(tmr_std_matrices)
            do np = 1, bundle%groups_on_node
               if (bundle%nm_nodgroup(np) > 0) then
                  do ni = 1, bundle%nm_nodgroup(np)
                     iprim = iprim + 1
                     !i=id_glob(index_my_atoms(iprim))
                     !i=index_my_atoms(iprim)
                     i = bundle%ig_prim(iprim)
                     do isf = 1, natomf_species(bundle%species(iprim))
                        p_force(dir1, i) = p_force(dir1, i) - &
                             return_matrix_value(mat_tmp, np, ni, 0, 0, isf, isf, 1)
                        if (flag_stress) then
                          if (flag_full_stress) then
                            PP_stress(dir1,dir2) = PP_stress(dir1,dir2) - &  
                               return_matrix_value(mat_tmp2, np, ni, 0, 0, isf, isf, 1)
                            if (flag_atomic_stress) then
                              atomic_stress(dir1,dir2,i) = &
                                atomic_stress(dir1,dir2,i) - &
                                return_matrix_value(mat_tmp2, np, ni, 0, 0, isf, isf, 1)
                            end if
                          else
                            PP_stress(dir1,dir1) = PP_stress(dir1,dir1) - &  
                              return_matrix_value(mat_tmp2, np, ni, 0, 0, isf, isf, 1)
                          end if
                        end if

                        if(iprint_MD>3) phi_pulay_for(dir1, i) = phi_pulay_for(dir1, i) - &
                             return_matrix_value(mat_tmp, np, ni, 0, 0, isf, isf, 1)
                     end do ! isf
                  end do ! ni
               end if ! if the partition has atoms
            end do ! np
            call stop_timer(tmr_std_matrices)
            t1 = mtime()
            if (inode == ionode .and. iprint_MD > 3) &
                 write (io_lun, fmt='(10x,"Phi Pulay sum ",i4," time: ",f12.5)')&
                       dir1, t1-t0
            t0 = t1
          end do ! dir2
       end do ! dir1
       call free_temp_fn_on_grid(tmp_fn2)
       call free_temp_fn_on_grid(tmp_fn)
       call free_temp_matrix(mat_tmp2)
       call free_temp_matrix(mat_tmp)
    end if ! if (WhichPulay == BothPulay .OR. WhichPulay == PhiPulay .OR.&
           !  & (WhichPulay == SPulay .AND. flag_basis_set == blips)) then

    ! NB this used to be PAOs only (flag_basis_set == PAOs)
    if (WhichPulay == BothPulay .or. WhichPulay == SPulay) then
       mat_tmp = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
       tmp_fn = allocate_temp_fn_on_grid(atomf)
       if (atomf==sf) then
          do spin = 1, nspin
             matM12atomf(spin) = matM12(spin)
          enddo
       else
          do spin = 1, nspin
             matM12atomf(spin) = allocate_temp_matrix(aSa_range, 0, atomf, atomf)
             call SF_to_AtomF_transform(matM12(spin), matM12atomf(spin), spin, Srange)
          enddo
       endif
       do dir1 = 1, 3
          ! Call assemble to generate dS_ij/dR_kl
          if (flag_basis_set == blips) then
             call blip_to_grad_new(inode-1, dir1, tmp_fn)
             call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem),&
                  mat_tmp, atomfns, tmp_fn)
             call matrix_scale(minus_one, mat_tmp)
          else if (flag_basis_set == PAOs) then
             call assemble_deriv_2(dir1, aSa_range, mat_tmp, 1)
          end if
          ! For each primary set atom, we want \sum_j dS_ij G_ij (I think)
          iprim = 0
          call start_timer(tmr_std_matrices)
          do np = 1, bundle%groups_on_node
             if (bundle%nm_nodgroup(np) > 0) then
                do ni = 1, bundle%nm_nodgroup(np)
                   iprim = iprim + 1
                   !i=id_glob(index_my_atoms(iprim))
                   !i=index_my_atoms(iprim)
                   i = bundle%ig_prim(iprim)
                   do neigh = 1, mat(np,aSa_range)%n_nab(ni)
                      ist = mat(np,aSa_range)%i_acc(ni) + neigh - 1
                      gcspart = &
                           BCS_parts%icover_ibeg(mat(np,aSa_range)%i_part(ist)) + &
                           mat(np,aSa_range)%i_seq(ist) - 1
                      r_str(1)=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                      r_str(2)=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                      r_str(3)=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                      ! matM12(1) here is just used to work out the
                      ! position in matrix, this position will be
                      ! identical for matM12(nspin), so only matM12(1)
                      ! needed here
                      wheremat = matrix_pos(matM12atomf(1), iprim,           &
                                            halo(aSa_range)%i_halo(gcspart), &
                                            1, 1)
                      if (wheremat /= mat(np,aSa_range)%onsite(ni)) then
                         do isf = 1, mat(np,aSa_range)%ndimj(ist)
                            do jsf = 1, mat(np,aSa_range)%ndimi(ni)
                               matM12_value = zero
                               do spin = 1, nspin
                                  matM12_value = &
                                       matM12_value + spin_factor * &
                                       return_matrix_value(matM12atomf(spin), &
                                       np, ni, iprim, neigh, jsf, isf)
                               end do ! spin
                               ! the factor of two here comes from
                               ! chain rule in derivatives of K
                               ! respect to phi
                               thisG_dS_dR = two*matM12_value * return_matrix_value(mat_tmp, np, ni, iprim, neigh, jsf, isf)
                               p_force(dir1,i) = p_force(dir1,i) + thisG_dS_dR
                               if(iprint_MD>3) s_pulay_for(dir1,i) = s_pulay_for(dir1,i) + thisG_dS_dR
                               if (flag_stress) then
                                 if (flag_full_stress) then
                                   do dir2=1,3
                                     SP_stress(dir1,dir2) = &
                                       SP_stress(dir1,dir2) + &
                                       thisG_dS_dR * r_str(dir2)
                                       if (flag_atomic_stress) then
                                         atomic_stress(dir1,dir2,i) = &
                                           atomic_stress(dir1,dir2,i) + &
                                           thisG_dS_dR * r_str(dir2) * half
                                       end if
                                   end do
                                 else
                                   SP_stress(dir1,dir1) = &
                                     SP_stress(dir1,dir1) + &
                                     thisG_dS_dR * r_str(dir1)
                                 end if ! flag_full_stress
                               end if ! flag_stress
                            end do ! jsf
                         end do ! isf
                      end if ! (wheremat /= mat(np,aSa_range)%onsite(ni))
                   end do ! neigh
                end do ! ni
             end if ! (bundle%nm_nodgroup(np) > 0)
          end do ! np
          call stop_timer(tmr_std_matrices)
          t1 = mtime ()
          if (inode == ionode .AND. iprint_MD > 3) &
               write (io_lun,*) 'S Pulay ', dir1, ' time: ', t1 - t0
          t0 = t1
       end do ! dir1

       if (atomf.ne.sf) then
          do spin = nspin,1,-1 ! Inverse order required
             call free_temp_matrix(matM12atomf(spin))
          enddo
       endif
       call free_temp_fn_on_grid(tmp_fn)
       call free_temp_matrix(mat_tmp)
    end if

    !  In principle, the summation below is not needed.
    !  p_force should be calculated only for my primary set of atoms.
    call gsum(p_force, 3, n_atoms)
    if(iprint_MD>3) then
       call gsum(s_pulay_for,3,n_atoms)
       call gsum(phi_pulay_for,3,n_atoms)
    end if
    if (flag_stress) then
      call gsum(PP_stress,3,3)
      call gsum(SP_stress,3,3)
      SP_stress = half*SP_stress
      if (flag_basis_set == blips .and. flag_analytic_blip_int) &
        KE_stress = half*KE_stress
    end if
    ! NB we do NOT need to halve PP_stress because it is from an integral on the grid

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='pulay_force',echo=.true.)
!****lat>$

    return
  end subroutine pulay_force
  !!***

! -----------------------------------------------------------
! Subroutine get_HF_force
! -----------------------------------------------------------

!!****f* force_module/get_HF_force *
!!
!!  NAME
!!   get_HF_force
!!  USAGE
!!
!!  PURPOSE
!!   Gets the Hellman-Feynman contribution to the
!!   atomic forces. For this both the electron density
!!   and Hartree potential are needed.
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   Summer 1995
!!  MODIFICATION HISTORY
!!   12/9/95 by Chris Goringe to use new data distribution
!!   22/1/96 by EH to correct.
!!   24/4/96 by EH to adapt to new data structure.
!!   27/11/97 by DRB to correct expressions
!!   21/05/2001 dave
!!    Added ROBODoc header, converted to F90 and stripped
!!    subroutine calls
!!   08/06/2001 dave
!!    Changed to use gsum from GenComms and added RCS Id
!!    and Log tags
!!   12/06/2001 dave
!!    Included in force_module
!!   2004/10/05 drb
!!    Added hartree_module use
!!   2008/03/03 18:46 dave
!!    Changed float to real
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!   2017/08/29 jack baker & dave
!!    Removed r_super_x references (redundant)
!!  TODO
!!    Fix this so that it doesn't loop over all processors ! Follow
!!    set_pseudo 13/05/2002 dave
!!  SOURCE
!!
  subroutine get_HF_force(hf_force, density, n_atoms, size)

    use datatypes
    use numbers
    use dimens, only: n_my_grid_points, x_grid, y_grid, z_grid, &
         grid_point_volume
    use grid_index, only: grid_point_block, grid_point_x, &
         grid_point_y, grid_point_z
    use species_module, only: species, charge
    use pseudopotential_data, only: ps_exponent, core_radius_2, &
         radius_max, n_points_max, local_pseudopotential, &
         d2_local_pseudopotential
    use global_module, only: rcellx,rcelly,rcellz,id_glob, &
         species_glob, nlpf, ni_in_cell
    use block_module, only : nx_in_block,ny_in_block,nz_in_block, &
         n_pts_in_block
    use group_module, only : blocks, parts
    use primary_module, only: domain
    use cover_module, only: DCS_parts
    use set_blipgrid_module, only : naba_atoms_of_blocks
    use GenComms, only: my_barrier, cq_abort, inode, ionode, gsum
    use hartree_module, only: hartree
    use maxima_module, only: maxngrid

    implicit none

    ! Passed variables
    integer :: n_atoms, size

    real(double), dimension(:,:) :: hf_force
    real(double), dimension(:) :: density

    ! Local variables
    integer :: i, j, my_block, n, the_species, iatom

    real(double) :: a, b, da, db, dc, dd, r, alpha, beta, derivative, &
         fx_1, fx_2, fy_1, fy_2, fz_1, fz_2, gauss, h_energy, rx, ry, &
         rz, x_point, y_point, z_point, r2, r_from_i, x, y, z, step, &
         q, elec_here, local_potential, c, d

    real(double):: dcellx_block, dcelly_block, dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    integer :: ipart, jpart, ind_part, ia, ii, icover, ig_atom
    real(double):: xatom, yatom, zatom
    real(double):: xblock, yblock, zblock
    real(double) :: dx, dy, dz
    integer :: ix, iy, iz, iblock, ipoint, igrid

    ! Automatic
    real(double) :: h_potential(size)
    type(cq_timer) :: backtrace_timer

!    the_species = 1
!    step = radius_max(1)/real(n_points_max(1)-1,double)
!    do i=1,n_points_max(1)-1
!       r = i*step+0.5*step
!       a = ( real(i+1,double)*step - r ) / step
!       b = one - a
!       c = a * ( a * a - one ) * step * step / six
!       d = b * ( b * b - one ) * step * step / six
!       da = -one / step
!       db =  one / step
!       dc = -step * ( three * a * a - one ) / six
!       dd =  step * ( three * b * b - one ) / six
!       local_potential =   &
!            a * local_pseudopotential(i,the_species)+  &
!            b * local_pseudopotential(i+1,the_species)+ &
!            c * d2_local_pseudopotential(i,the_species)+  &
!            d * d2_local_pseudopotential(i+1,the_species)
!       derivative =  &
!            da * local_pseudopotential(i,the_species) +  &
!            db * local_pseudopotential(i+1,the_species) + &
!            dc * d2_local_pseudopotential(i,the_species) +  &
!            dd * d2_local_pseudopotential(i+1,the_species)
!    end do

    call start_backtrace(t=backtrace_timer,who='get_HF_force',where=7,level=3,echo=.true.)
    ! get Hartree potential
    HF_force = zero

    call hartree (density, h_potential, maxngrid, h_energy)

    dcellx_block = rcellx / blocks%ngcellx
    dcelly_block = rcelly / blocks%ngcelly
    dcellz_block = rcellz / blocks%ngcellz

    dcellx_grid = dcellx_block / nx_in_block
    dcelly_grid = dcelly_block / ny_in_block
    dcellz_grid = dcellz_block / nz_in_block

    call my_barrier()

    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock = (domain%idisp_primx(iblock) + domain%nx_origin - 1) *&
            & dcellx_block
       yblock = (domain%idisp_primy(iblock) + domain%ny_origin - 1) *&
            & dcelly_block
       zblock = (domain%idisp_primz(iblock) + domain%nz_origin - 1) *&
            & dcellz_block

       if (naba_atoms_of_blocks(nlpf)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom = 0
          do ipart = 1, naba_atoms_of_blocks(nlpf)%no_of_part(iblock)
             jpart = naba_atoms_of_blocks(nlpf)%list_part(ipart, iblock)
             if(jpart > DCS_parts%mx_gcover) then
                call cq_abort ('set_ps: JPART ERROR ', ipart, jpart)
             end if
             ind_part = DCS_parts%lab_cell(jpart)
             do ia = 1, naba_atoms_of_blocks(nlpf)%no_atom_on_part(ipart, iblock)
                iatom = iatom + 1
                ii = naba_atoms_of_blocks(nlpf)%list_atom(iatom, iblock)
                icover = DCS_parts%icover_ibeg(jpart) + ii - 1
                ig_atom = id_glob(parts%icell_beg(ind_part) + ii - 1)
                if (parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                   call cq_abort ('set_ps: globID ERROR ',&
                        & ii, parts %icell_beg(ind_part))
                end if
                if (icover > DCS_parts%mx_mcover) then
                   call cq_abort ('set_ps: icover ERROR ',&
                        & icover, DCS_parts%mx_mcover)
                end if
                xatom = DCS_parts%xcover(icover)
                yatom = DCS_parts%ycover(icover)
                zatom = DCS_parts%zcover(icover)
                ! the_species = species(ig_atom)
                the_species = species_glob(ig_atom)
                q = charge(the_species)
                alpha = ps_exponent(the_species)
                beta = (alpha / pi)**1.5_double * grid_point_volume
                step = radius_max(the_species) /&
                     & real (n_points_max(the_species)-1, double)
                ipoint = 0
                do iz = 1, nz_in_block
                   do iy = 1, ny_in_block
                      do ix = 1, nx_in_block
                         ipoint = ipoint + 1
                         igrid = n_pts_in_block * (iblock-1) + ipoint
                         if (igrid > n_my_grid_points) then
                            call cq_abort ('set_ps: igrid error ',&
                                 & igrid, n_my_grid_points)
                         end if
                         elec_here = density(igrid) * grid_point_volume
                         dx = dcellx_grid * (ix-1)
                         dy = dcelly_grid * (iy-1)
                         dz = dcellz_grid * (iz-1)
                         rx = xblock + dx - xatom
                         ry = yblock + dy - yatom
                         rz = zblock + dz - zatom
                         r2 = rx * rx + ry * ry + rz * rz
                         gauss = exp (-alpha * r2)
                         if (r2 < core_radius_2(the_species)) then
                            r_from_i = sqrt( r2 )
                            !if ( r_from_i > zero ) then
                            if (r_from_i > RD_ERR) then
                               x = rx / r_from_i
                               y = ry / r_from_i
                               z = rz / r_from_i
                            else
                               x = zero
                               y = zero
                               z = zero
                            end if
                            j = aint (r_from_i / step) + 1
                            ! check j
                            if (j > n_points_max(the_species) - 1) then
                               call cq_abort ('set_ps: overrun problem', j)
                            end if
                            r = real (j, double) * step
                            a = (r - r_from_i) / step
                            b = one - a
                            da = -one / step
                            db =  one / step
                            dc = -step * (three * a * a - one) / six
                            dd =  step * (three * b * b - one) / six
                            derivative = da * &
                                 local_pseudopotential(j, &
                                 the_species) + db * &
                                 local_pseudopotential(j + 1, &
                                 the_species) + dc * &
                                 d2_local_pseudopotential(j, &
                                 the_species) + dd * &
                                 d2_local_pseudopotential(j + 1, &
                                 the_species)
                            ! This is a little obscure, but the
                            ! multiples below have had the minus sign
                            ! which should be there removed to correct
                            ! an earlier on dropped in the da, db, dc,
                            ! dd expressions.DRB 27.11.97
                            fx_1 = x * derivative
                            fy_1 = y * derivative
                            fz_1 = z * derivative
                            fx_2 = two * alpha * beta * rx * gauss * &
                                   h_potential(igrid) * q
                            fy_2 = two * alpha * beta * ry * gauss * &
                                   h_potential(igrid) * q
                            fz_2 = two * alpha * beta * rz * gauss * &
                                   h_potential(igrid) * q
                         else
                            fx_1 = zero
                            fy_1 = zero
                            fz_1 = zero
                            fx_2 = zero
                            fy_2 = zero
                            fz_2 = zero
                         end if ! if (r2 < core_radius_2(the_species)) then
                         HF_force(1,ig_atom) = HF_force(1,ig_atom) + &
                                fx_1 * elec_here
                         HF_force(2,ig_atom) = HF_force(2,ig_atom) + &
                                fy_1 * elec_here
                         HF_force(3,ig_atom) = HF_force(3,ig_atom) + &
                                fz_1 * elec_here
                         HF_force(1,ig_atom) = HF_force(1,ig_atom) + &
                                fx_2
                         HF_force(2,ig_atom) = HF_force(2,ig_atom) + &
                                fy_2
                         HF_force(3,ig_atom) = HF_force(3,ig_atom) + &
                                fz_2
                      end do !ix
                   end do  !iy
                end do   !iz
             end do ! naba_atoms
          end do ! naba_part
       end if !(naba_atoms_of_blocks(nlpf)%no_of_part(iblock) > 0) !naba atoms?
    end do ! iblock : primary set of blocks

!    ! now loop over grid points and accumulate HF part of the force
!    do n=1, n_my_grid_points
!       my_block = grid_point_block(n)
!       elec_here = density(n) * grid_point_volume
!       x_point = (grid_point_x(n) - 1) * x_grid
!       y_point = (grid_point_y(n) - 1) * y_grid
!       z_point = (grid_point_z(n) - 1) * z_grid
!       do i=1, n_atoms
!          q = charge(species(i))
!          alpha = exponent(species(i))
!          beta = (alpha / pi)**1.5_double * grid_point_volume
!          rx = ( x_point - rionx(i) )
!          ry = ( y_point - riony(i) )
!          rz = ( z_point - rionz(i) )
!          rx = rx - anint(rx)
!          ry = ry - anint(ry)
!          rz = rz - anint(rz)
!          rx = rx * r_super_x
!          ry = ry * r_super_y
!          rz = rz * r_super_z
!          r2 = rx * rx + ry * ry + rz * rz
!          gauss = dexp( -alpha * r2 )
!          if ( r2 .lt. core_radius_2(species(i)) ) then
!             r_from_i = sqrt( r2 )
!             if ( r_from_i .gt. zero ) then
!                x = rx / r_from_i
!                y = ry / r_from_i
!                z = rz / r_from_i
!             else
!                x = zero
!                y = zero
!                z = zero
!             end if
!             ! now we construct the derivative of the local part of the
!             ! pseudopotential by spline interpolation from a table
!             step = radius_max(species(i)) / &
!                  real( n_points_max(species(i)) - 1, double)
!             j = aint( r_from_i / step ) + 1
!
!             ! with this we can now use the spline interpolation tables to
!             ! construct the derivative of the local part of the
!             ! pseudopotential on the grid
!             r = real(j,double) * step
!             a = ( r - r_from_i ) / step
!             b = one - a
!             da = -one / step
!             db =  one / step
!             dc = -step * ( three * a * a - one ) / six
!             dd =  step * ( three * b * b - one ) / six
!             derivative =  &
!                  da * local_pseudopotential(j,species(i)) +  &
!                  db * local_pseudopotential(j+1,species(i)) + &
!                  dc * d2_local_pseudopotential(j,species(i)) +  &
!                  dd * d2_local_pseudopotential(j+1,species(i))
!
!             ! This is a little obscure, but the multiples below have had
!             ! the minus sign which should be there removed to correct an
!             ! earlier on dropped in the da, db, dc, dd expressions.DRB 27.11.97
!             fx_1 = x * derivative
!             fy_1 = y * derivative
!             fz_1 = z * derivative
!             fx_2 = two * alpha * beta * rx * gauss * h_potential( n ) * q
!             fy_2 = two * alpha * beta * ry * gauss * h_potential( n ) * q
!             fz_2 = two * alpha * beta * rz * gauss * h_potential( n ) * q
!          else
!             fx_1 = zero
!             fy_1 = zero
!             fz_1 = zero
!             fx_2 = zero
!             fy_2 = zero
!             fz_2 = zero
!          end if
!          HF_force(1,i) = HF_force(1,i) + fx_1 * elec_here
!          HF_force(2,i) = HF_force(2,i) + fy_1 * elec_here
!          HF_force(3,i) = HF_force(3,i) + fz_1 * elec_here
!!            fx_2 = two * alpha * beta * rx * gauss * h_potential( n ) * q
!!            fy_2 = two * alpha * beta * ry * gauss * h_potential( n ) * q
!!            fz_2 = two * alpha * beta * rz * gauss * h_potential( n ) * q
!          HF_force(1,i) = HF_force(1,i) + fx_2
!          HF_force(2,i) = HF_force(2,i) + fy_2
!          HF_force(3,i) = HF_force(3,i) + fz_2
!       end do
!    end do

    ! and add contributions from all nodes
!    do i = 1, n_atoms
!       write(io_lun,*) inode,i,(HF_force(j,i),j=1,3)
!    end do
    call gsum (HF_force, 3, n_atoms)
!    do i = 1, n_atoms
!       write(io_lun,*) inode,i,(HF_force(j,i),j=1,3)
!    end do
    call stop_backtrace(t=backtrace_timer,who='get_HF_force',echo=.true.)

    return
  end subroutine get_HF_force
!!***

  ! -----------------------------------------------------------
  ! Subroutine get_HF_non_local_force
  ! -----------------------------------------------------------

  !!****f* force_module/get_HF_non_local_force *
  !!
  !!  NAME
  !!   get_HF_non_local_force
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Gets the non-local part of the HF force
  !!
  !!   There is a little wrangling about this - it's not entirely clear
  !!   whether it's a Hellmann-Feynman force or a Pulay force. Whatever
  !!   it's called, it is due to the change in the support/projector
  !!   overlap matrix with no dependence on the charge density.
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   29/4/96
  !!  MODIFICATION HISTORY
  !!   24/2/97 CMG
  !!    We no longer carry grad support around, it needs to be generated on the fly
  !!   04/09/2000 TM
  !!    new get_matrix_elements
  !!   15/01/2001 TM
  !!    new blip_to_grad
  !!   11/05/01 and 14/05/01 DRB
  !!    Removed many of the arguments passed and redundant local variables
  !!   17/05/2001 dave
  !!    Changed subroutine call to blip_to_grad
  !!   21/05/2001 dave
  !!    Reduced subroutine call and use statements
  !!   22/05/2001 dave
  !!    Added more discussion to the PURPOSE comment field
  !!   11/06/2001 dave
  !!    Added RCS Id and Log tags and GenComms
  !!   11/06/2001 dave
  !!    Corrected to add data_U finding lines (required for correct forces)
  !!   12/06/2001 dave
  !!    Included in force_module
  !!   21/06/2001 dave
  !!    Removed use matrix_diag (now in this module)
  !!    Removed matrix_scaling
  !!   14:20, 02/09/2003 drb
  !!    Added flag to select HF or Pulay or both (for easier force testing)
  !!   08:12, 2003/10/01 dave
  !!    Added structure to prepare for PAO basis set
  !!   2011/11/22 L.Tong
  !!    Added spin polarisation
  !!    - changed the definition of HF_NL_force to be more aligned to
  !!      F90 syntax style and as an assumed shape array
  !!   2011/12/12 17:26 dave
  !!    Added force calculation for analytic blips
  !!   2012/03/26 L.Tong
  !!   - Changed spin implementation
  !!   2015/05/08 14:33 dave and sym
  !!    Adding phi-Pulay NL stress
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_nlpf_rem -> atomf_nlpf_rem
  !!   2016/07/29 18:30 nakata
  !!    Renamed supports_on_atom -> blips_on_atom
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2017/02/21 16:30 nakata
  !!    Removed unused PAO_to_grad
  !!   2019/05/08 zamaan
  !!    Added atomic stress contributions
  !!  SOURCE
  !!
  subroutine get_HF_non_local_force(HF_NL_force, what_force, n_atoms)

    use datatypes
    use numbers
    use dimens
    use GenBlas
    use matrix_data,                 only: APrange, PArange, mat, halo
    use mult_module,                 only: AP_trans, mult, aHa_AP_AP,  &
                                           matrix_product,             &
                                           matrix_scale,               &
                                           matrix_transpose,           &
                                           allocate_temp_matrix,       &
                                           free_temp_matrix,           &
                                           matU, matUT, matPA, matAP,  &
                                           matKatomf,                  &
                                           return_matrix_value,        &
                                           store_matrix_value
    use species_module,              only: species
    use pseudopotential_data,        only: pseudopotential_derivatives
    use set_bucket_module,           only: rem_bucket, atomf_nlpf_rem
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use blip_grid_transform_module,  only: blip_to_grad_new
    use GenComms,                    only: gsum, cq_abort, inode,      &
                                           ionode
    ! TM new pseudo
    use pseudopotential_common,      only: pseudo_type, OLDPS, SIESTA, &
                                           STATE, ABINIT
    use pseudo_tm_module,            only: nonloc_pp_derivative_tm
    use global_module,               only: iprint_MD, flag_basis_set,  &
                                           blips, PAOs, atomf, nlpf,   &
                                           nspin, spin_factor,         &
                                           id_glob, species_glob,      &
                                           flag_analytic_blip_int,     &
                                           ni_in_cell, flag_full_stress, &
                                           flag_stress, flag_atomic_stress, &
                                           area_moveatoms
    ! TEMP
    use build_PAO_matrices,          only: assemble_deriv_2
    use group_module,                only: parts
    use primary_module,              only: bundle
    use cover_module,                only: BCS_parts
    use GenComms,                    only: gsum, myid
    use functions_on_grid,           only: allocate_temp_fn_on_grid,   &
                                           free_temp_fn_on_grid,       &
                                           atomfns,                    &
                                           H_on_atomfns, pseudofns
    use nlpf2blip,                   only: get_dSP, nlpf_on_atom
    use primary_module ,             only: bundle
    use support_spec_format,         only: blips_on_atom
    use group_module,                only: parts
    use primary_module,              only: bundle
    use cover_module,                only: BCS_parts
    use species_module,              only: nsf_species, nlpf_species
    use memory_module,               only: reg_alloc_mem, reg_dealloc_mem,  &
                                           type_dbl

    implicit none

    ! Passed variables
    integer :: n_atoms, what_force
    real(double), dimension(3,n_atoms) :: HF_NL_force

    ! Local variables

    ! data_dPA will be obtained by a global trans of data_dPA_t
    !  this trans is same as the one used in U => UT
    !     Tsuyoshi Miyazaki 28/12/2000
    integer, dimension(3) :: matdAP, matdPA
    integer, dimension(3,3) :: matdAPr, matdPAr
    integer      :: dir1, dir2, k, l, stat, dpseudofns, np, nn, i, i1, i2, &
                    spec, this_nsf, this_nlpf, ni, spin
    integer      :: iprim, gcspart, ist, nab, neigh_global_num,        &
                    neigh_global_part, neigh_species, wheremat, isf, jsf
    real(double) :: dx, dy, dz, thisdAP
    real(double), dimension(3,3) :: NL_P_stress, NL_HF_stress
    real(double), dimension(:,:,:), allocatable :: NL_P_atomic, NL_HF_atomic
    real(double), dimension(3) :: r_str
    type(cq_timer) :: backtrace_timer

    call start_backtrace(t=backtrace_timer,who='get_HF_non_local_force',where=7,level=3,echo=.true.)

    do k = 1, 3
       matdAP(k)  = allocate_temp_matrix (APrange, AP_trans, atomf, nlpf)
       matdPA(k)  = allocate_temp_matrix (PArange, AP_trans, nlpf, atomf)
       if (flag_stress) then
         if (flag_full_stress) then
           do l = 1,3
              matdAPr(k,l) = &
                allocate_temp_matrix(APrange, AP_trans, atomf, nlpf)
              matdPAr(k,l) = &
                allocate_temp_matrix(PArange, AP_trans, nlpf, atomf)
           end do
         else
           matdAPr(k,k) = allocate_temp_matrix(APrange, AP_trans, atomf, nlpf)
           matdPAr(k,k) = allocate_temp_matrix(PArange, AP_trans, nlpf, atomf)
         end if
       end if
    end do
    if (flag_basis_set == blips .and. (.not. flag_analytic_blip_int)) then
       dpseudofns = allocate_temp_fn_on_grid(nlpf)
    end if
    if (flag_atomic_stress) then
      allocate(NL_HF_atomic(3,3,ni_in_cell), STAT=stat)
      if (stat /= 0) &
        call cq_abort("Error allocating NL atomic stress: ", ni_in_cell)
      allocate(NL_P_atomic(3,3,ni_in_cell), STAT=stat)
      if (stat /= 0) &
        call cq_abort("Error allocating NL atomic stress: ", ni_in_cell)
      call reg_alloc_mem(area_moveatoms, 2*3*3*ni_in_cell, type_dbl)
    end if
    HF_NL_force = zero
    NL_P_stress = zero
    NL_HF_stress = zero
    NL_stress = zero

    ! to save memory we do each direction in turn...
    do dir1 = 1, 3
       call matrix_scale(zero, matdAP(dir1))
       call matrix_scale(zero, matdPA(dir1))
       if (flag_stress) then
         if (flag_full_stress) then
           do dir2 = 1, 3
             call matrix_scale(zero, matdAPr(dir1,dir2))
             call matrix_scale(zero, matdPAr(dir1,dir2))
           end do
         else
           call matrix_scale(zero, matdAPr(dir1,dir1))
           call matrix_scale(zero, matdPAr(dir1,dir1))
         end if
       end if
       if ((flag_basis_set == blips) .and. flag_analytic_blip_int) then
          iprim=0
          do np = 1, bundle%groups_on_node
             if (bundle%nm_nodgroup(np) > 0) then
                do ni = 1, bundle%nm_nodgroup(np)
                   iprim = iprim + 1
                   spec = bundle%species(iprim)
                   ! Loop over neighbours of atom
                   do nab = 1, mat(np,APrange)%n_nab(ni)
                      ist = mat(np,APrange)%i_acc(ni) + nab - 1
                      ! Build the distances between atoms - needed for phases
                      gcspart = &
                           BCS_parts%icover_ibeg(mat(np,APrange)%i_part(ist)) + &
                           mat(np,APrange)%i_seq(ist) - 1
                      ! Displacement vector
                      dx = BCS_parts%xcover(gcspart) - bundle%xprim(iprim)
                      dy = BCS_parts%ycover(gcspart) - bundle%yprim(iprim)
                      dz = BCS_parts%zcover(gcspart) - bundle%zprim(iprim)
                      ! We need to know the species of neighbour
                      neigh_global_part = &
                           BCS_parts%lab_cell(mat(np,APrange)%i_part(ist))
                      neigh_global_num  = &
                           id_glob(parts%icell_beg(neigh_global_part) + &
                                   mat(np,APrange)%i_seq(ist) - 1)
                      neigh_species = species_glob(neigh_global_num)
                      if ((dx*dx + dy*dy + dz*dz) > RD_ERR) then
                         call get_dSP(blips_on_atom(iprim),            &
                                      nlpf_on_atom(neigh_species),     &
                                      matdAP(dir1), iprim,        &
                                      halo(APrange)% i_halo(gcspart),  &
                                      dx, dy, dz, spec, neigh_species, &
                                      dir1)
                      end if
                   end do ! nab
                end do ! ni
             end if ! (bundle%nm_nodgroup(np) > 0)
          end do ! np
          call matrix_transpose(matdAP(dir1), matdPA(dir1))
          call matrix_scale(-one, matdAP(dir1))
       else if (flag_basis_set == blips .and. &
                (.not. flag_analytic_blip_int)) then
          ! first of all construct the dir1al derivatives of the core
          ! functions on the grid
          select case (pseudo_type)
          case (OLDPS)
             call pseudopotential_derivatives(dir1, dpseudofns)
          case (SIESTA)
             call nonloc_pp_derivative_tm(dir1, dpseudofns)
          case (ABINIT)
             call nonloc_pp_derivative_tm(dir1, dpseudofns)
          end select
          !if(dir1.eq.3) then
          !   stop
          !else
          !   continue
          !end if
          ! now calculate overlap between these derivatives and support functions
          ! New get_matrix_elements  TSUYOSHI MIYAZAKI 28/Dec/2000
          call get_matrix_elements_new(inode-1,                    &
                                       rem_bucket(atomf_nlpf_rem), &
                                       matdAP(dir1), atomfns, &
                                       dpseudofns)
          call matrix_transpose(matdAP(dir1), matdPA(dir1))
          call matrix_scale(-one, matdPA(dir1))
          ! and the overlap between the derivatives of the support and
          ! the core functions, H_on_atomfns(1) is used as a work array
          call blip_to_grad_new(inode-1, dir1, H_on_atomfns(1))
          call get_matrix_elements_new(inode-1,                    &
                                       rem_bucket(atomf_nlpf_rem), &
                                       matdAP(dir1),          &
                                       H_on_atomfns(1),            &
                                       pseudofns)
          call matrix_scale(-one, matdAP(dir1))
          ! actually, because the derivatives we want are with respect to the
          ! position of the atom supporting each support function, we must
          ! scale the above integrals by -1
       else if (flag_basis_set == PAOs) then
          ! Get matrix elements between derivative of projectors and
          ! atom functions (support function or PAO)
          call assemble_deriv_2(dir1, APrange, matdAP(dir1), 3)
          call matrix_scale(-one, matdAP(dir1))

          ! Get matrix elements between projectors and derivative of
          ! support functions BY TRANSPOSE. This is fine (I think)
          ! because the transpose ought to be exact
          call matrix_transpose(matdAP(dir1), matdPA(dir1))
          call matrix_scale(-one, matdPA(dir1))
       else
          call cq_abort("get_HF_NL_force: basis set undefined ", &
                        flag_basis_set)
       end if ! ((flag_basis_set == blips) .and. flag_analytic_blip_int)
       ! Now scale dAP and dPA by R_{ji} for stress
       ! 2014/08/06 11:45 Shereif
       iprim = 0
       do np = 1, bundle%groups_on_node
          if (bundle%nm_nodgroup(np) > 0) then
             do ni = 1, bundle%nm_nodgroup(np)
                iprim = iprim + 1
                do nab = 1, mat(np,APrange)%n_nab(ni)
                   ist = mat(np,APrange)%i_acc(ni) + nab - 1
                   gcspart = BCS_parts%icover_ibeg(mat(np,APrange)%i_part(ist)) + mat(np,APrange)%i_seq(ist) - 1
                   r_str(1)=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                   r_str(2)=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                   r_str(3)=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                   do isf = 1, mat(np,APrange)%ndimj(ist)
                      do jsf = 1, mat(np,APrange)%ndimi(ni)
                        thisdAP = return_matrix_value(matdAP(dir1), np, ni, iprim, nab, jsf, isf)
                        if (flag_stress) then
                          if (flag_full_stress) then
                            do dir2=1,3
                               call store_matrix_value(matdAPr(dir1,dir2), np, ni, iprim, nab, jsf, isf, r_str(dir2)*thisdAP)
                            end do
                          else
                            call store_matrix_value(matdAPr(dir1,dir1), np, ni, iprim, nab, jsf, isf, r_str(dir1)*thisdAP)
                          end if
                        end if
                           ! Reuse variable
                           !thisdAP = return_matrix_value(matdPA(dir1), np, ni, iprim, nab, jsf, isf)
                           !call store_matrix_value(matdPAr(dir1), np, ni, iprim, nab, jsf, isf, r_str*thisdAP)
                     end do
                   end do
                end do
                !do nab = 1, mat(np,PArange)%n_nab(ni)
                !   ist = mat(np,PArange)%i_acc(ni) + nab - 1
                !   gcspart = BCS_parts%icover_ibeg(mat(np,PArange)%i_part(ist)) + mat(np,PArange)%i_seq(ist) - 1
                !   if(dir1==1) then
                !      r_str=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                !   else if(dir1==2) then
                !      r_str=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                !   else if(dir1==3) then
                !      r_str=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                !   end if
                !   do isf = 1, mat(np,PArange)%ndimj(ist)
                !      do jsf = 1, mat(np,PArange)%ndimi(ni)
                !         thisdAP = return_matrix_value(matdPA(dir1), np, ni, iprim, nab, jsf, isf)
                !         call store_matrix_value(matdPAr(dir1), np, ni, iprim, nab, jsf, isf, r_str*thisdAP)
                !      end do
                !   end do
                !end do
             end do
          end if
       end do
       if (flag_stress) then
         if (flag_full_stress) then
           do dir2=1,3
             call matrix_transpose(matdAPr(dir1,dir2), matdPAr(dir1,dir2))
           end do
         else
           call matrix_transpose(matdAPr(dir1,dir1), matdPAr(dir1,dir1))
         end if
       end if
       !call matrix_scale(-one, matdPAr(direction))
    end do ! dir1
    if (flag_basis_set == blips .and. (.not.flag_analytic_blip_int)) &
         call free_temp_fn_on_grid(dpseudofns)
    ! First of all, find U (= K.AP)
    call matrix_transpose(matAP, matPA)
    do spin = 1, nspin
       call matrix_product(matKatomf(spin), matPA, matU(spin), mult(aHa_AP_AP))
       ! we factor matU(1) with spin_factor for spin non-polarised
       ! calculation. This may cause slight confusion, as the factor
       ! are not applied to other matrices of the spin channels, and
       ! matU(1) is supposed to be U for spin up, but doing this way
       ! we make the calculations more efficient.
       call matrix_scale(minus_two * spin_factor, matU(spin))
       ! Find the transpose of U (= K.AP)
       call matrix_transpose(matU(spin), matUT(spin))
    end do

    ! For PAOs, we need to zero the on-site terms

    ! Evaluate the pulay term - due to the phis changing
    ! NB We want sum_k dAP_ik.UT_ki, so let's do sum_k dAP_ik.U_ik
    if (what_force == Pulay .OR. what_force == HF_and_Pulay) then
       do k = 1, 3
          do spin = 1, nspin
             ! Note that matrix_diagonal accumulates HF_NL_force(k,:)
             call matrix_diagonal(matdAP(k), matU(spin), &
                                  HF_NL_force(k,:), APrange, inode)
             if (flag_stress) then
               if (flag_full_stress) then
                 do l = 1, 3
                   if (flag_atomic_stress) then
                     call matrix_diagonal_stress(matdAPr(k,l), matU(spin), &
                          NL_P_stress(k,l), APrange, inode,atomic_stress(k,l,:))
                   else
                     call matrix_diagonal_stress(matdAPr(k,l), matU(spin), &
                          NL_P_stress(k,l), APrange, inode)
                   end if
                 end do
               else
                 call matrix_diagonal_stress(matdAPr(k,k), matU(spin), &
                                      NL_P_stress(k,k), APrange, inode)
               end if
             end if
          end do ! spin
       end do ! k
    end if
    ! Evaluate the Hellmann-Feynman term - due to the chis changing
    if (what_force == HF .or. what_force == HF_and_Pulay) then
       do k = 1, 3
          do spin = 1, nspin
             ! Note that matrix_diagonal accumulates HF_NL_force(k,:)
             call matrix_diagonal(matdPA(k), matUT(spin), &
                                  HF_NL_force(k,:), PArange,inode)
             if (flag_stress) then
               if (flag_full_stress) then
                 do l = 1, 3
                   if (flag_atomic_stress) then
                     call matrix_diagonal_stress(matdPAr(k,l), matUT(spin), &
                          NL_HF_stress(k,l), PArange,inode,atomic_stress(k,l,:))
                   else
                     call matrix_diagonal_stress(matdPAr(k,l), matUT(spin), &
                          NL_HF_stress(k,l), PArange, inode)
                   end if
                 end do
               else
                 call matrix_diagonal_stress(matdPAr(k,k), matUT(spin), &
                      NL_HF_stress(k,k), PArange, inode)
               end if
             end if
          end do ! spin
       end do ! k
    end if

    call gsum(HF_NL_force, 3, n_atoms)
    if (flag_atomic_stress) then
      deallocate(NL_HF_atomic, NL_P_atomic, STAT=stat)
      if (stat /= 0) &
        call cq_abort("Error deallocating NL atomic stress: ", ni_in_cell)
      call reg_dealloc_mem(area_moveatoms, 2*3*3*ni_in_cell, type_dbl)
    end if
    if (flag_stress) then
       do k=1,3
          if (flag_full_stress) then
            do l=1,3
              NL_stress(k,l) = NL_stress(k,l) + &
                               half*(NL_P_stress(k,l) + NL_HF_stress(k,l))
            end do
          else
            NL_stress(k,k) = NL_stress(k,k) + &
                             half*(NL_P_stress(k,k) + NL_HF_stress(k,k))
          end if
       end do
      call gsum(NL_stress,3,3)
      ! Note that the atomic contributions are computed in flag_diagonal_stress
      ! for memory efficiency
    end if
    do k = 3, 1, -1
       if (flag_stress) then
         if (flag_full_stress) then
           do l = 3, 1, -1
             call free_temp_matrix(matdPAr(k,l))
             call free_temp_matrix(matdAPr(k,l))
           end do
         else
           call free_temp_matrix(matdPAr(k,k))
           call free_temp_matrix(matdAPr(k,k))
         end if
       end if
       call free_temp_matrix(matdPA(k))
       call free_temp_matrix(matdAP(k))
    end do

    call stop_backtrace(t=backtrace_timer,who='get_HF_non_local_force',echo=.true.)

    return
  end subroutine get_HF_non_local_force
  !!***

  ! -----------------------------------------------------------
  ! Subroutine get_KE_force
  ! -----------------------------------------------------------

  !!****f* force_module/get_KE_force *
  !!
  !!  NAME
  !!   get_KE_force
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Evaluates the force due to the change in the kinetic
  !!   energy matrix when an atom moves. This is therefore
  !!   the kinetic pulay type I force. As the onsite KE is evaluated
  !!   on the blip grid, there is _no_ contribution; the contribution
  !!   is due to the change in the intersite terms.
  !!
  !!   If we just have a single support function per atom...
  !! (the extension to NSF!=1 is standard...)
  !!   The gradient of the energy when atom i moves in direction \alpha
  !!   is given by a sum over atoms j, directions \beta and the grid
  !!
  !! (d/dR_{i\alpha}) KE = 2 \sum_j \sum_{\beta} \sum_{grid}
  !!                         Grad_{\beta} psi_j  *
  !!                       (d/dR_{i\alpha}) Grad_{\beta} psi_i  *
  !!                         K_{ij}   * grid_point_volume             [i!=j]
  !!
  !!   This can be written alternatively (more usefully) as
  !!
  !! (d/dR_{i\alpha}) KE = (d/dR_{i\alpha}) T_{ij} K_{ij}           [i!=j]
  !!
  !!   where
  !!
  !! (d/dR_{i\alpha}) T_{ij} = 2 \sum_{\beta} \sum_{grid}
  !!                             grid_point_volume * Grad_{\beta} psi_j  *
  !!                           (d/dR_{i\alpha}) Grad_{\beta} psi_i
  !!
  !!   So within the loop over direction (for the force) we need a loop
  !!   over the components of the grad (\beta). For each of these we
  !!   evaluate Grad_{\beta} psi_j [blip_to_grad] and
  !! (d/dR_{i\alpha}) Grad_{\beta} psi_i [blip_to_gradgrad], and thus
  !! (d/dR_{i\alpha}) T_{ij} [get_matrix_elements].
  !!   We then loop through the list of overlap range interactions, and
  !!   accumulate the force due to that interaction onto the atom i.
  !!
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   C.M.Goringe
  !!  CREATION DATE
  !!   03/03/97
  !!  MODIFICATION HISTORY
  !!   15/08/2000 TM
  !!    New get_matrix_elements
  !!   15/01/2001 TM
  !!    New blip_to_grad
  !!   15/05/2001 dave
  !!    ROBODoc headers, new get_matrix_elements
  !!   17/05/2001 dave
  !!    Shortened calls to blip_to_grad and blip_to_gradgrad
  !!   21/05/2001 dave
  !!    Shortened overall subroutine call
  !!   11/06/2001 dave
  !!    Added RCS Id and Log tags and GenComms
  !!   12/06/2001 dave
  !!    Included in force_module
  !!   08:16, 2003/10/01 dave
  !!    Added structure for preparing for PAO basis set
  !!   2011/11/24 L.Tong
  !!    Added spin polarisation
  !!   2012/03/26 L.Tong
  !!   - Changed spin implementation
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_H_sf_rem -> atomf_H_atomf_rem
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!   2019/04/04 zamaan
  !!    Added off-diagonal elements of KE_stress
  !!   2019/05/08 zamaan
  !!    Added atomic stress contributions
  !!  SOURCE
  !!
  subroutine get_KE_force(KE_force, n_atoms)

    use datatypes
    use numbers
    use primary_module,              only: bundle
    use matrix_module,               only: matrix, matrix_halo
    use matrix_data,                 only: mat, aHa_range, halo
    use cover_module,                only: BCS_parts
    use mult_module,                 only: allocate_temp_matrix,      &
                                           free_temp_matrix,          &
                                           return_matrix_value, matKatomf, &
                                           matrix_pos
    use GenBlas
    use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use blip_grid_transform_module,  only: blip_to_grad_new,          &
                                           blip_to_gradgrad_new
    use GenComms,                    only: my_barrier, gsum, inode,   &
                                           ionode
    use global_module,               only: iprint_MD, flag_basis_set, &
                                           blips, PAOs, atomf,        &
                                           flag_onsite_blip_ana,      &
                                           nspin, spin_factor,        &
                                           flag_full_stress, flag_stress, &
                                           flag_atomic_stress
    use build_PAO_matrices,          only: assemble_deriv_2
    use functions_on_grid,           only: H_on_atomfns,              &
                                           allocate_temp_fn_on_grid,  &
                                           free_temp_fn_on_grid

    implicit none

    ! Passed variables
    integer :: n_atoms
    real(double), dimension(3, n_atoms) :: KE_force
    ! local variables
    integer :: i, j, grad_direction, dir1, dir2, element, np, nn,&
               atom, n1, n2, mat_grad_T, ist, gcspart, iprim, tmp_fn, &
               spin
    real(double) :: thisK_gradT
    real(double), dimension(3) :: r_str
    type(cq_timer) :: backtrace_timer

    call start_backtrace(t=backtrace_timer,who='get_KE_force',where=7,level=3,echo=.true.)

!    ! First, clear the diagonal blocks of data K; this is the easiest way
!    ! to avoid doing the onsite terms
    KE_force = zero
    mat_grad_T = allocate_temp_matrix(aHa_range,0)
    ! Now, for the offsite part, done on the integration grid.
    if (flag_basis_set == blips) then

       tmp_fn = allocate_temp_fn_on_grid(atomf)
       call start_timer(tmr_std_matrices)
       do grad_direction = 1, 3

          ! get the Grad_{\beta} psi_j term
          call blip_to_grad_new(inode-1, grad_direction, H_on_atomfns(1))
          do dir1 = 1, 3
             ! get the (d/dR_{i\alpha}) Grad_{\beta} psi_i term
             call blip_to_gradgrad_new(inode-1, grad_direction, &
                                       dir1, tmp_fn)
             ! now get the (d/dR_{i\alpha}) T_{ij} term from the above
             !new matrix_elements
             call get_matrix_elements_new(inode-1,                       &
                                          rem_bucket(atomf_H_atomf_rem), &
                                          mat_grad_T, tmp_fn,            &
                                          H_on_atomfns(1))
             iprim = 0
             do np = 1,bundle%groups_on_node
                if (bundle%nm_nodgroup(np) > 0) then
                   do i = 1, bundle%nm_nodgroup(np)
                      ! The numbering of the ig_prim index is the same
                      ! as a sequential index running over atoms in
                      ! the primary set (seen in pulay force for
                      ! instance) DRB and TM 10:38, 29/08/2003
                      iprim = iprim+1
                      atom = bundle%ig_prim(bundle%nm_nodbeg(np)+i-1)
                      do j = 1,mat(np,aHa_range)%n_nab(i)
                         ist = mat(np,aHa_range)%i_acc(i)+j-1
                         gcspart =                                                &
                              BCS_parts%icover_ibeg(mat(np,aHa_range)%i_part(ist)) + &
                              mat(np,aHa_range)%i_seq(ist) - 1
                         r_str(1)=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                         r_str(2)=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                         r_str(3)=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                         ! matKatomf(1) here is just to work out the matrix
                         ! position, it is going to be the same as
                         ! matKatomf(nspin)
                         element =                       &
                              matrix_pos(matKatomf(1), iprim, &
                                         halo(aHa_range)%i_halo(gcspart), 1, 1)
                         if ((.not. flag_onsite_blip_ana) .or. &
                             (element /= mat(np,aHa_range)%onsite(i))) then
                            do n1 = 1, mat(np,aHa_range)%ndimj(ist)
                               do n2 = 1, mat(np,aHa_range)%ndimi(i)
                                  do spin = 1, nspin
                                     thisK_gradT = spin_factor * & 
                                        return_matrix_value(matKatomf(spin), & 
                                        np, i, iprim, j, n2, n1) *  &
                                        return_matrix_value(mat_grad_T,   &
                                        np, i, iprim, &
                                        j, n2, n1)
                                     KE_force(dir1,atom) = &
                                          KE_force(dir1,atom) + thisK_gradT
                                     if (flag_stress) then
                                       if (flag_full_stress) then
                                         do dir2=1,3
                                           KE_stress(dir1,dir2) = &
                                                KE_stress(dir1,dir2) + &
                                                thisK_gradT * r_str(dir2)
                                           if (flag_atomic_stress) then
                                             atomic_stress(dir1,dir2,atom) = &
                                               atomic_stress(dir1,dir2,atom) +&
                                               thisK_gradT * r_str(dir2) * half

                                           end if
                                         end do
                                       else
                                         KE_stress(dir1,dir1) = &
                                              KE_stress(dir1,dir1) + &
                                              thisK_gradT * r_str(dir1)
                                       end if ! flag_full_stress
                                     end if ! flag_stress
                                  end do ! spin
                               end do ! n2
                            end do ! n1
                         end if
                      end do ! j = mat%n_nab
                   end do ! i = bundle%nm_nodgroup
                end if  ! (bundle%nm_nodgroup(np) > 0)
             end do ! np = bundle%groups_on_node
          end do ! force directions
       end do ! grad directions
       call stop_timer(tmr_std_matrices)
       call free_temp_fn_on_grid(tmp_fn)

    else if (flag_basis_set==PAOs) then

       do dir1 = 1, 3
          ! Build derivatives
          call assemble_deriv_2(dir1,aHa_range, mat_grad_T, 2)
          call start_timer(tmr_std_matrices)
          iprim = 0
          do np = 1, bundle%groups_on_node
             if (bundle%nm_nodgroup(np) > 0) then
                do i = 1, bundle%nm_nodgroup(np)
                   ! The numbering of the ig_prim index is the same as
                   ! a sequential index running over atoms in the
                   ! primary set (seen in pulay force for instance)
                   ! DRB and TM 10:38, 29/08/2003
                   atom = bundle%ig_prim(bundle%nm_nodbeg(np)+i-1)
                   iprim = iprim + 1
                   do j = 1, mat(np,aHa_range)%n_nab(i)
                      ist = mat(np,aHa_range)%i_acc(i) + j - 1
                      gcspart = &
                        BCS_parts%icover_ibeg(mat(np,aHa_range)%i_part(ist)) +&
                        mat(np,aHa_range)%i_seq(ist) - 1
                      r_str(1)=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                      r_str(2)=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                      r_str(3)=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                      ! matKatomf(1) is used to just to get the position
                      element = matrix_pos(matKatomf(1), iprim, &
                                halo(aHa_range)%i_halo(gcspart), 1, 1)
                      if (element /= mat(np,aHa_range)%onsite(i)) then
                         do n1 = 1, mat(np,aHa_range)%ndimj(ist)
                            do n2 = 1, mat(np,aHa_range)%ndimi(i)
                               do spin = 1, nspin
                                  thisK_gradT = spin_factor * &
                                       return_matrix_value(matKatomf(spin), &
                                       np, i, iprim, j, n2, n1) *  &
                                       return_matrix_value(mat_grad_T,   &
                                       np, i, iprim, j, n2, n1)
                                  KE_force(dir1,atom) =       &
                                       KE_force(dir1,atom) +  &
                                       thisK_gradT
                                  if (flag_stress) then
                                    if (flag_full_stress) then
                                      do dir2=1,3
                                        KE_stress(dir1,dir2) = &
                                             KE_stress(dir1,dir2) + &
                                             thisK_gradT * r_str(dir2)
                                        if (flag_atomic_stress) then
                                          atomic_stress(dir1,dir2,atom) = &
                                            atomic_stress(dir1,dir2,atom) + &
                                            thisK_gradT * r_str(dir2) * half
                                        end if
                                      end do ! dir2
                                    else
                                      KE_stress(dir1,dir1) = &
                                           KE_stress(dir1,dir1) + &
                                           thisK_gradT * r_str(dir1)
                                    end if ! flag_full_stress
                                  end if ! flag_stress
                               end do ! spin
                            end do ! n2
                         end do ! n1
                      end if ! (element /= mat(np,aHa_range)%onsite(i))
                   end do ! j
                end do ! i
             end if  ! (bundle%nm_nodgroup(np) > 0)
          end do ! np
          call stop_timer(tmr_std_matrices)
       end do ! force directions
    end if ! (flag_basis_set == blips)

    call gsum(KE_force, 3, n_atoms)
    if (flag_stress) then
      KE_stress = half*KE_stress
      call gsum(KE_stress, 3, 3)
    end if

      call free_temp_matrix(mat_grad_T)
      call stop_backtrace(t=backtrace_timer,who='get_KE_force',echo=.true.)

      return
    end subroutine get_KE_force
    !!***

    ! -----------------------------------------------------------
    ! Subroutine get_HNA_force
    ! -----------------------------------------------------------

    !!****f* force_module/get_HNA_force *
    !!
    !!  NAME
    !!   get_HNA_force
    !!  USAGE
    !!
    !!  PURPOSE
    !!   Gets the neutral atom part of the HF force if using projectors
    !!   This mixes Hellman-Feynman and Pulay forces (as with NL part above)
    !!  INPUTS
    !!
    !!
    !!  USES
    !!
    !!  AUTHOR
    !!   D. R. Bowler
    !!  CREATION DATE
    !!   2018/01/10
    !!  MODIFICATION HISTORY
  !!   2018/01/25 12:52 JST dave
  !!    Changed transpose type for mat_dNA to aNAa_trans
  !!   2018/01/30 10:06 dave
  !!    Bug fix (found by Jack Baker): removed premature gsum call on NA_stress
  !!    which led to erroneous stress values on multiple processors
  !!   2019/05/08 zamaan
  !!    Added atomic stress contributions
  !!  SOURCE
  !!
  subroutine get_HNA_force(NA_force)

    use datatypes
    use numbers
    use primary_module,              only: bundle
    use matrix_module,               only: matrix, matrix_halo
    use matrix_data,                 only: mat, aHa_range, halo, aNArange, NAarange, Hrange
    use cover_module,                only: BCS_parts
    use mult_module,                 only: allocate_temp_matrix,      &
                                           free_temp_matrix,          &
                                           return_matrix_value, matKatomf, &
                                           matrix_pos, matrix_transpose, matrix_sum, &
                                           matrix_product, matrix_scale,               &
                                           S_trans, scale_matrix_value, aHa_aNA_aNA, aNA_trans, &
                                           matU, matUT, matNAa, mataNA, aNAa_trans, &
                                           return_matrix_value, store_matrix_value, mult, matK, matrix_product_trace
    use GenBlas
    use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use blip_grid_transform_module,  only: blip_to_grad_new,          &
                                           blip_to_gradgrad_new
    use GenComms,                    only: my_barrier, gsum, inode,   &
                                           ionode
    use global_module,               only: iprint_MD, flag_basis_set, &
                                           blips, PAOs, atomf,        &
                                           flag_onsite_blip_ana,      &
                                           nspin, spin_factor, ni_in_cell, &
                                           napf, id_glob, flag_full_stress, &
                                           flag_stress, flag_atomic_stress
    use build_PAO_matrices,          only: assemble_deriv_2
    use functions_on_grid,           only: H_on_atomfns,              &
                                           allocate_temp_fn_on_grid,  &
                                           free_temp_fn_on_grid
    use group_module, ONLY: parts

    implicit none

    ! Passed variables
    real(double), dimension(3, ni_in_cell) :: NA_force
    ! local variables
    integer :: i, j, grad_direction, element, np, nn,&
               atom, n1, n2, mat_grad_T, ist, gcspart, iprim, tmp_fn, &
               spin, mat_dNA, mat_dNAT, j_atom, dir1, dir2
    integer :: ip, nsf1, nsf2
    real(double) :: thisK_gradT

    real(double), dimension(3) :: r_str
    integer, dimension(3) :: matdaNA, matdNAa
    integer, dimension(3,3) :: matdaNAr, matDNAar
    integer      :: k, l, stat, dpseudofns, i1, i2, &
         spec, this_nsf, this_nlpf, ni, isf, jsf
    integer, dimension(2) :: matKzero
    integer      :: nab, neigh_global_num,        &
                    neigh_global_part, neigh_species, wheremat, matU_NA, matUT_NA
    real(double) :: dx, dy, dz, thisdAP, locforce
    real(double), dimension(3,3) :: NA_P_stress, NA_HF_stress
    real(double), dimension(:), allocatable :: force_contrib, f_c2
    type(cq_timer) :: backtrace_timer
    
    call start_backtrace(t=backtrace_timer,who='get_HNA_force',where=7,level=3,echo=.true.)

!    ! First, clear the diagonal blocks of data K; this is the easiest way
!    ! to avoid doing the onsite terms
    NA_force = zero
    allocate(force_contrib(ni_in_cell))
    force_contrib = zero
    ! 3-centre terms
    do k = 1, 3
       matdaNA(k) = allocate_temp_matrix (aNArange, aNA_trans, atomf, napf)
       matdNAa(k) = allocate_temp_matrix (NAarange, aNA_trans, napf, atomf)
       if (flag_stress) then
         if (flag_full_stress) then
           do l = 1, 3
             matdaNAr(k,l) = &
               allocate_temp_matrix (aNArange, aNA_trans, atomf, napf)
             matdNAar(k,l) = &
               allocate_temp_matrix (NAarange, aNA_trans, napf, atomf)
           end do
         else
             matdaNAr(k,k) = &
               allocate_temp_matrix (aNArange, aNA_trans, atomf, napf)
             matdNAar(k,k) = &
               allocate_temp_matrix (NAarange, aNA_trans, napf, atomf)
         end if
       end if
    end do
    matU_NA = allocate_temp_matrix(aNArange,aNA_trans,atomf,napf)
    matUT_NA = allocate_temp_matrix(NAarange,aNA_trans,napf,atomf)
    do spin = 1,nspin
       matKzero(spin) = allocate_temp_matrix(aHa_range,0,atomf,atomf)
    end do
    NA_P_stress = zero
    NA_HF_stress = zero
    NA_stress = zero

    do spin = 1,nspin
       call matrix_sum(zero,matKzero(spin),one,matKatomf(spin))
    end do
    ip = 0
    do np = 1,bundle%groups_on_node
       if (bundle%nm_nodgroup(np) > 0) then
          do i = 1,bundle%nm_nodgroup(np)
             ip = ip+1
             do nsf1=1,mat(np,aHa_range)%ndimi(i)
                do nsf2=1,mat(np,aHa_range)%ndimi(i)
                   call scale_matrix_value(matKzero(1), np, i, ip, 0, nsf1, nsf2, zero,1)
                   if(nspin==2) call scale_matrix_value(matKzero(2), np, i, ip, 0, nsf1, nsf2, zero,1)
                end do
             end do
          end do
       end if
    end do
    ! to save memory we do each direction in turn...
    do dir1 = 1, 3
       call matrix_scale(zero, matdaNA(dir1))
       call matrix_scale(zero, matdNAa(dir1))
       if (flag_stress) then
         if (flag_full_stress) then
           do dir2 = 1, 3
              call matrix_scale(zero, matdaNAr(dir1,dir2))
              call matrix_scale(zero, matdNAar(dir1,dir2))
           end do
         else
           call matrix_scale(zero, matdaNAr(dir1,dir1))
           call matrix_scale(zero, matdNAar(dir1,dir1))
         end if
       end if
          ! Get matrix elements between derivative of projectors and
          ! support functions
       call assemble_deriv_2(dir1, aNArange, matdaNA(dir1), 5)
       call matrix_scale(-one, matdaNA(dir1))

       ! Get matrix elements between projectors and derivative of
       ! support functions BY TRANSPOSE. This is fine (I think)
       ! because the transpose ought to be exact
       call matrix_transpose(matdaNA(dir1), matdNAa(dir1))
       call matrix_scale(-one, matdNAa(dir1))
       ! Now scale dAP and dPA by R_{ji} for stress
       ! 2014/08/06 11:45 Shereif
       iprim = 0
       do np = 1, bundle%groups_on_node
          if (bundle%nm_nodgroup(np) > 0) then
             do ni = 1, bundle%nm_nodgroup(np)
                iprim = iprim + 1
                do nab = 1, mat(np,aNArange)%n_nab(ni)
                   ist = mat(np,aNArange)%i_acc(ni) + nab - 1
                   gcspart = BCS_parts%icover_ibeg(mat(np,aNArange)%i_part(ist)) + mat(np,aNArange)%i_seq(ist) - 1
                   r_str(1)=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                   r_str(2)=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                   r_str(3)=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                   do isf = 1, mat(np,aNArange)%ndimj(ist)
                      do jsf = 1, mat(np,aNArange)%ndimi(ni)
                         thisdAP = return_matrix_value(matdaNA(dir1), np, ni, iprim, nab, jsf, isf)
                         if (flag_stress) then
                           if (flag_full_stress) then
                             do dir2=1,3
                               call store_matrix_value(matdaNAr(dir1,dir2), &
                                 np, ni, iprim, nab, jsf, isf, &
                                 r_str(dir2)*thisdAP)
                             end do
                           else
                             call store_matrix_value(matdaNAr(dir1,dir1), &
                               np, ni, iprim, nab, jsf, isf, &
                               r_str(dir1)*thisdAP)
                           end if
                         end if
                      end do
                   end do
                end do
             end do
          end if
       end do
       if (flag_stress) then
         if (flag_full_stress) then
           do dir2=1,3
             call matrix_transpose(matdaNAr(dir1,dir2), matdNAar(dir1,dir2))
           end do
         else
           call matrix_transpose(matdaNAr(dir1,dir1), matdNAar(dir1,dir1))
         end if
       end if
       
       
    end do ! Now end the direction loop

    
    
    ! First of all, find U (= K.SC)
    call matrix_transpose(mataNA, matNAa)
    do spin = 1, nspin
       call matrix_scale(zero, matU_NA)
       call matrix_product(matKzero(spin), matNAa, matU_NA, mult(aHa_aNA_aNA))
       !call matrix_product(matKatomf(spin), matNAa, matU_NA, mult(aHa_aNA_aNA))

      ! we factor matU(1) with spin_factor for spin non-polarised
      ! calculation. This may cause slight confusion, as the factor
      ! are not applied to other matrices of the spin channels, and
      ! matU(1) is supposed to be U for spin up, but doing this way
      ! we make the calculations more efficient.
      call matrix_scale(minus_two * spin_factor, matU_NA)
      ! Find the transpose of U (= K.SC)
      call matrix_transpose(matU_NA, matUT_NA)

      ! For PAOs, we need to zero the on-site terms
      ! Evaluate the pulay term - due to the phis changing
      ! NB We want sum_k dSC_ik.UT_ki, so let's do sum_k dSC_ik.U_ik
      !if (what_force == Pulay .OR. what_force == HF_and_Pulay) then
        do k = 1, 3
          ! Note that matrix_diagonal accumulates HF_NA_force(k,:)
          call matrix_diagonal(matdaNA(k), matU_NA, &
            NA_force(k,:), aNArange, inode)
          if (flag_stress) then
            if (flag_full_stress) then
              do l = 1, 3
                if (flag_atomic_stress) then
                  call matrix_diagonal_stress(matdaNAr(k,l), matU_NA, &
                    NA_P_stress(k,l), aNArange, inode, atomic_stress(k,l,:))
                else
                  call matrix_diagonal_stress(matdaNAr(k,l), matU_NA, &
                    NA_P_stress(k,l), aNArange, inode)
                end if
              end do
            else
              call matrix_diagonal_stress(matdaNAr(k,k), matU_NA, &
                NA_P_stress(k,k), aNArange, inode)
            end if
          end if
        end do ! k
        !end if
      ! Evaluate the Hellmann-Feynman term - due to the chis changing

      !if (what_force == HF .or. what_force == HF_and_Pulay) then
        do k = 1, 3
          ! Note that matrix_diagonal accumulates HF_NA_force(k,:)
          call matrix_diagonal(matdNAa(k), matUT_NA, &
            NA_force(k,:), NAarange,inode)
          if (flag_stress) then
            if (flag_full_stress) then
              do l = 1, 3
                if (flag_atomic_stress) then
                  call matrix_diagonal_stress(matdNAar(k,l), matUT_NA, &
                    NA_HF_stress(k,l), NAarange,inode, atomic_stress(k,l,:))
                else
                  call matrix_diagonal_stress(matdNAar(k,l), matUT_NA, &
                    NA_HF_stress(k,l), NAarange,inode)
                end if
              end do
            else
              call matrix_diagonal_stress(matdNAar(k,k), matUT_NA, &
                NA_HF_stress(k,k), NAarange,inode)
            end if
          end if
        end do ! k
      !end if
    end do ! spin

    if (flag_stress) then
      do k=1,3
        if (flag_full_stress) then
          do l=1,3
            NA_stress(k,l) = NA_stress(k,l) + &
                             half*(NA_P_stress(k,l) + NA_HF_stress(k,l))
          end do
        else
          NA_stress(k,k) = NA_stress(k,k) + &
                           half*(NA_P_stress(k,k) + NA_HF_stress(k,k))
        end if
      end do
    end if

    do spin = nspin,1,-1
       call free_temp_matrix(matKzero(spin))
    end do
    call free_temp_matrix(matUT_NA)
    call free_temp_matrix(matU_NA)
    do k = 3, 1, -1
       if (flag_stress) then
         if (flag_full_stress) then
           do l = 3, 1, -1
             call free_temp_matrix(matdNAar(k,l))
             call free_temp_matrix(matdaNAr(k,l))
           end do
         else
           call free_temp_matrix(matdNAar(k,k))
           call free_temp_matrix(matdaNAr(k,k))
         end if
       end if
       call free_temp_matrix(matdNAa(k))
       call free_temp_matrix(matdaNA(k))
    end do

    ! 1- and 2-centre terms
    mat_dNA = allocate_temp_matrix(aHa_range,aNAa_trans,atomf,atomf)
    mat_dNAT = allocate_temp_matrix(aHa_range,aNAa_trans,atomf,atomf)
    ! Now, for the offsite part, done on the integration grid.
    do dir1 = 1, 3
       force_contrib = zero
       ! Build derivatives
       call assemble_deriv_2(dir1,aHa_range, mat_dNA, 4)
       call matrix_transpose(mat_dNA, mat_dNAT)
       call matrix_sum(one,mat_dNA,-one,mat_dNAT)
       ! Add < i | VNA_j | i >
       call start_timer(tmr_std_matrices)
       iprim = 0
       do np = 1, bundle%groups_on_node
          if (bundle%nm_nodgroup(np) > 0) then
             do i = 1, bundle%nm_nodgroup(np)
                ! The numbering of the ig_prim index is the same as
                ! a sequential index running over atoms in the
                ! primary set (seen in pulay force for instance)
                ! DRB and TM 10:38, 29/08/2003
                atom = bundle%ig_prim(bundle%nm_nodbeg(np)+i-1)
                iprim = iprim + 1
                do j = 1, mat(np,aHa_range)%n_nab(i)
                   ist = mat(np,aHa_range)%i_acc(i) + j - 1
                   gcspart =                                                &
                        BCS_parts%icover_ibeg(mat(np,aHa_range)%i_part(ist)) + &
                        mat(np,aHa_range)%i_seq(ist) - 1
                   r_str(1)=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                   r_str(2)=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                   r_str(3)=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                   ! matKatomf(1) is used to just to get the position
                   element =                    &
                        matrix_pos(matKatomf(1), iprim, &
                        halo(aHa_range)%i_halo(gcspart), 1, 1)
                   if (element /= mat(np,aHa_range)%onsite(i)) then
                      do n1 = 1, mat(np,aHa_range)%ndimj(ist)
                         do n2 = 1, mat(np,aHa_range)%ndimi(i)
                            do spin = 1, nspin
                               thisK_gradT = spin_factor * &
                                    return_matrix_value(matKatomf(spin),   &
                                    np, i, iprim, j, n2, n1) *  &
                                    return_matrix_value(mat_dNA,   &
                                    np, i, iprim, j, n2, n1)
                               NA_force(dir1,atom) =       &
                                    NA_force(dir1,atom) + two*thisK_gradT
                               force_contrib(atom) = force_contrib(atom) + two*thisK_gradT
                               if (flag_stress) then
                                 if (flag_full_stress) then
                                   do dir2=1,3
                                     NA_stress(dir1,dir2) = &
                                       NA_stress(dir1,dir2) + &
                                       thisK_gradT * r_str(dir2)
                                     if (flag_atomic_stress) then
                                       atomic_stress(dir1,dir2,atom) = &
                                         atomic_stress(dir1,dir2,atom) + &
                                         thisK_gradT * r_str(dir2)
                                     end if
                                   end do ! dir2
                                 else
                                   NA_stress(dir1,dir1) = &
                                     NA_stress(dir1,dir1) + &
                                     thisK_gradT * r_str(dir1)
                                 end if ! flag_full_stress
                               end if ! flag_stress
                            end do ! spin
                         end do ! n2
                      end do ! n1
                   end if
                end do ! j
             end do ! i
          end if  ! (bundle%nm_nodgroup(np) > 0)
       end do ! np
       call matrix_scale(zero,mat_dNAT)
       call assemble_deriv_2(dir1,aHa_range, mat_dNAT, 6)
       force_contrib = zero
       !f_c2 = zero
       iprim = 0
       do np = 1, bundle%groups_on_node
          if (bundle%nm_nodgroup(np) > 0) then
             do i = 1, bundle%nm_nodgroup(np)
                ! The numbering of the ig_prim index is the same as
                ! a sequential index running over atoms in the
                ! primary set (seen in pulay force for instance)
                ! DRB and TM 10:38, 29/08/2003
                atom = bundle%ig_prim(bundle%nm_nodbeg(np)+i-1)
                iprim = iprim + 1
                do j = 1, mat(np,aHa_range)%n_nab(i)
                   ist = mat(np,aHa_range)%i_acc(i) + j - 1
                   gcspart =                                                &
                        BCS_parts%icover_ibeg(mat(np,aHa_range)%i_part(ist)) + &
                        mat(np,aHa_range)%i_seq(ist) - 1
                   r_str(1)=BCS_parts%xcover(gcspart)-bundle%xprim(iprim)
                   r_str(2)=BCS_parts%ycover(gcspart)-bundle%yprim(iprim)
                   r_str(3)=BCS_parts%zcover(gcspart)-bundle%zprim(iprim)
                   ! Global number of neighbour: id_glob( parts%icell_beg(gcs%lab_cell(np)) +ni-1 )
                   j_atom = id_glob( parts%icell_beg( BCS_parts%lab_cell(mat(np,aHa_range)%i_part(ist))) &
                        + mat(np,aHa_range)%i_seq(ist)-1 )
                   NA_force(dir1,j_atom) =       &
                      NA_force(dir1,j_atom) - &
                      spin_factor * return_matrix_value(mat_dNAT,   &
                      np, i, iprim, j, 1, 1)
                   if (flag_stress) then
                     if (flag_full_stress) then
                       do dir2=1,3
                         NA_stress(dir1,dir2) = NA_stress(dir1,dir2) + &
                            spin_factor * return_matrix_value(mat_dNAT,   &
                            np, i, iprim, j, 1, 1) * r_str(dir2)
                         if (flag_atomic_stress) then
                           atomic_stress(dir1,dir2,j_atom) = &
                             atomic_stress(dir1,dir2,j_atom) + &
                             spin_factor * return_matrix_value(mat_dNAT, &
                             np, i, iprim, j, 1, 1) * r_str(dir2)
                         end if
                       end do
                     else
                       NA_stress(dir1,dir1) = NA_stress(dir1,dir1) + &
                          spin_factor * return_matrix_value(mat_dNAT,   &
                          np, i, iprim, j, 1, 1) * r_str(dir1)
                     end if !flag_full_stress
                   end if ! flag_stress
                   force_contrib(j_atom) = force_contrib(j_atom)-spin_factor * return_matrix_value(mat_dNAT,   &
                        np, i, iprim, j, 1, 1)
                   force_contrib(atom) = force_contrib(atom)+spin_factor * return_matrix_value(mat_dNAT,   &
                         np, i, iprim, j, 1, 1)
                   end do ! j
                end do ! i
             end if  ! (bundle%nm_nodgroup(np) > 0)
          end do ! np
       call stop_timer(tmr_std_matrices)
    end do ! dir1

    call gsum(NA_force, 3, ni_in_cell)
    !KE_stress = half*KE_stress
    if (flag_stress) call gsum(NA_stress,3,3)
    ! Note that the atomic contributions are computed in flag_diagonal_stress
    ! for memory efficiency

    call free_temp_matrix(mat_dNAT)
    call free_temp_matrix(mat_dNA)
    deallocate(force_contrib)
    call stop_backtrace(t=backtrace_timer,who='get_HNA_force',echo=.true.)
    return
  end subroutine get_HNA_force
  !!***


! -----------------------------------------------------------
! Subroutine matrix_diagonal
! -----------------------------------------------------------

!!****f* force_module/matrix_diagonal *
!!
!!  NAME
!!   matrix_diagonal
!!  USAGE
!!
!!  PURPOSE
!!   Basically, this is designed for the non-local contribution to the
!!   forces, and does sum_k one_ik.two^T_ki as sum_k one_ik two_ik
!!   We're being _very_ bad, and assuming that i_elements is sorted on
!!   atom on node, so that the first n entries are for atom one on this
!!   node, etc.
!!
!!   Now (a) complies with new matrix mults and (b) is used by lots of
!!   the force routines
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   7/7/98 (Happy 5 months Rob !)
!!  MODIFICATION HISTORY
!!   6/6/00 dave
!!    Changed to use new matrix mults
!!   21/06/2001 dave
!!    Added ROBODoc header and included in force_module
!!  SOURCE
!!
  subroutine matrix_diagonal(matA, matB, diagonal, range, inode)

    use datatypes
    use primary_module, only : bundle
    use matrix_module, only: matrix, matrix_halo
    use matrix_data, only : mat, halo
    use cover_module, only: BCS_parts
    use mult_module, only: return_matrix_value_pos, matrix_pos
    use GenBlas

    implicit none

    ! Shared variables
    real(double) :: diagonal( :)

    integer :: range, inode,matA,matB

    ! Local variables
    integer :: iprim, np, i, j, atom,element, n1,n2, ist, gcspart
    
    iprim = 0
    call start_timer(tmr_std_matrices)
    do np = 1,bundle%groups_on_node
       if(bundle%nm_nodgroup(np) > 0) then
          do i=1,bundle%nm_nodgroup(np)
             iprim = iprim + 1
             atom = bundle%ig_prim(iprim)
             do j = 1,mat(np,range)%n_nab(i)
                ist = mat(np,range)%i_acc(i)+j-1
                gcspart = BCS_parts%icover_ibeg(mat(np,range)%i_part(ist))+mat(np,range)%i_seq(ist)-1
                do n2 = 1,mat(np,range)%ndimj(ist)
                   do n1 = 1,mat(np,range)%ndimi(i)
                      element = matrix_pos(matA,iprim,halo(range)%i_halo(gcspart),n1,n2)
                      diagonal(atom) = diagonal(atom) + &
                           return_matrix_value_pos(matA,element)*return_matrix_value_pos(matB,element)
                   end do
                end do
             end do
          end do
       end if  !  (bundle%nm_nodgroup(np) > 0)
    end do
    call stop_timer(tmr_std_matrices)
    return
  end subroutine matrix_diagonal
!!***

!!****f* force_module/matrix_diagonal_stress *
!!
!!  NAME
!!   matrix_diagonal
!!  USAGE
!!
!!  PURPOSE
!!   Does the same thing as matrix_diagonal, except it sums the 
!!   atomic contributions. Adapted from matrix_diagonal, for 
!!   memory efficiency
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2019/05/07
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine matrix_diagonal_stress(matA, matB, stress, range, inode, diagonal)

    use datatypes
    use global_module, only: atomic_stress
    use primary_module, only : bundle
    use matrix_module, only: matrix, matrix_halo
    use matrix_data, only : mat, halo
    use cover_module, only: BCS_parts
    use mult_module, only: return_matrix_value_pos, matrix_pos
    use GenBlas

    implicit none

    ! Shared variables
    real(double) :: stress
    real(double), dimension(:), optional :: diagonal
    integer :: range, inode, matA, matB

    ! Local variables
    integer :: iprim, np, i, j, element, n1, n2, ist, gcspart, atom
    
    iprim = 0
    call start_timer(tmr_std_matrices)
    do np = 1,bundle%groups_on_node
       if(bundle%nm_nodgroup(np) > 0) then
          do i=1,bundle%nm_nodgroup(np)
             iprim = iprim + 1
             atom = bundle%ig_prim(iprim)
             do j = 1,mat(np,range)%n_nab(i)
                ist = mat(np,range)%i_acc(i)+j-1
                gcspart = BCS_parts%icover_ibeg(mat(np,range)%i_part(ist))+mat(np,range)%i_seq(ist)-1
                do n2 = 1,mat(np,range)%ndimj(ist)
                   do n1 = 1,mat(np,range)%ndimi(i)
                      element = matrix_pos(matA,iprim,halo(range)%i_halo(gcspart),n1,n2)
                      stress = stress + &
                           return_matrix_value_pos(matA,element) * &
                           return_matrix_value_pos(matB,element)
                      if (present(diagonal)) then
                        diagonal(atom) = diagonal(atom) + &
                          return_matrix_value_pos(matA,element) * &
                          return_matrix_value_pos(matB,element) * half
                        ! the factor of half arises because the total 
                        ! non-local/neutral atom is 
                        ! half*(N*_P_stress + N*_HF_stress), but for memory
                        ! efficiency, I'm not storing the atomic contributions
                        ! in an array - zamaan
                      end if
                   end do
                end do
             end do
          end do
       end if  !  (bundle%nm_nodgroup(np) > 0)
    end do
    call stop_timer(tmr_std_matrices)
    return
  end subroutine matrix_diagonal_stress
!!***

  ! -----------------------------------------------------------
  ! Subroutine get_nonSC_correction_force
  ! -----------------------------------------------------------

  !!****f* force_module/get_nonSC_correction_force_nospin *
  !!
  !!  NAME
  !!   get_nonSC_correction_force_nospin
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Gets the correction to the forces for when we're doing
  !!   Harris-Foulkes (i.e. non self-consistent ab initio tight binding)
  !!   calculations.  Operates in exactly the same manner as the
  !!   Hellmann-Feynman forces for the local part of the pseudopotential
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:12, 20/03/2003 drb
  !!  MODIFICATION HISTORY
  !!   11:53, 30/09/2003 drb
  !!    Removed get_electronic_density (done above in force) and
  !!    simplified subroutine call
  !!   2004/10/05 drb
  !!    Added hartree_module use
  !!   14:08, 2006/07/17
  !!    Added functional type selector and PBE force
  !!   2008/11/13
  !!    Added new PBE functional types (revPBE, RPBE)
  !!   2011/03/31 M.Arita
  !!    Added the contribution from P.C.C.
  !!   2011/10/17 L.Tong
  !!    - Added spin polarisation
  !!    - Added optional parameter density_out_dn. The density_out array will
  !!      be used to store the total out density, where as a new local
  !!      arrays density_out_up is used to store the spin up component
  !!      of the out put density. BUT, IMPORTANT, upon input, the
  !!      density_out is used for the spin up component of the output
  !!      density, we copy this into density_out_up immediately inside
  !!      the routine and density_out is then used as total output. And
  !!      at the end of the array density_out is reverted back to the
  !!      values stored in density_out_up and density_out_up is
  !!      destroyed.
  !!    - Removed dependence on maxngrid, and moved all dependence on
  !!      the grid size to the input parameter size
  !!    - made potential to be taken from potential_module directly
  !!      instead of as an input parameter
  !!   2012/02/12 L.Tong
  !!    - Changed procedure name to get_nonSC_correction_force_nospin.
  !!      This subroutine is to be called from the interface
  !!      get_nonSC_correction_force_nospin.
  !!    - Removed all spin polarised calculations, this subroutine is to
  !!      be used for spin non-polarised calculations only.
  !!   2012/02/13 L.Tong
  !!    - Made h_potential automatic array.
  !!   2012/03/26 L.Tong
  !!    - Changed spin implementation.
  !!    - merged subroutines get_nonSC_correction_force_nospin and
  !!      get_nonSC_correction_force_spin
  !!    - renamed subroutine name back to get_nonSC_correction_force
  !!    - deleted the interface get_nonSC_correction_force
  !!   2012/05/29 L.Tong
  !!    - Cleaned up xc-functonal selector. Now spin and non-spin
  !!      calculations share the same calls. However so far only case
  !!      that works for spin polarised calculations is for LSDA-PW92.
  !!   2013/02/28 L.Tong
  !!    - Added warning for non-implemented PBE non-SC forces (for all
  !!      PBE variants.) The code still runs, but will set non SC
  !!      correction forces to 0
  !!   2013/07/10 11:35 dave
  !!    - Bug fix for sum over two components of rho even without spin
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2015/08/07 16:40 dave
  !!    - Fixed bug for GGA PBE force (used work_potential not potential)
  !!   2015/08/10 08:07 dave
  !!    Adding stress calculation
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/11/02 10:24 dave
  !!    - Subtle GGA error: return is inside if(inode==ionode) loop causing hang !
  !!   2018/02/13 11:52 dave
  !!    Changes for new, universal XC interface
  !!   2018/02/14 13:26 dave
  !!    More subtle errors ! The PCC, non-SCF XC stress did not have spin
  !!    factor applied (as above in non-SCF routine)
  !!   2019/05/08 zamaan
  !!    Added atomic stress contributions
  !!  SOURCE
  !!
  subroutine get_nonSC_correction_force(HF_force, density_out, inode, &
                                        ionode, n_atoms, nsize)

    use datatypes
    use numbers
    use species_module,      only: species
    use GenComms,            only: gsum
    use global_module,       only: rcellx, rcelly, rcellz, id_glob,    &
                                   ni_in_cell, species_glob, dens,     &
                                   area_moveatoms, IPRINT_TIME_THRES3, &
                                   flag_pcc_global, nspin, spin_factor, &
                                   flag_full_stress, flag_stress,      &
                                   flag_atomic_stress
    use XC,                  only: get_xc_potential,                   &
                                   get_dxc_potential,                  &
                                   flag_is_GGA
    use block_module,        only: nx_in_block, ny_in_block,           &
                                   nz_in_block, n_pts_in_block
    use group_module,        only: blocks, parts
    use primary_module,      only: domain
    use cover_module,        only: DCS_parts
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use GenComms,            only: my_barrier, cq_abort
    use atomic_density,      only: atomic_density_table
    use pseudo_tm_info,      only: pseudo
    use spline_module,       only: dsplint
    use dimens,              only: grid_point_volume, n_my_grid_points
    use GenBlas,             only: axpy
    use density_module,      only: density, density_scale, density_pcc
    use hartree_module,      only: hartree, Hartree_stress
    use potential_module,    only: potential
    use memory_module,       only: reg_alloc_mem, reg_dealloc_mem,     &
                                   type_dbl
    use timer_module,        only: cq_timer, start_timer, stop_timer,  &
                                   print_timer, stop_print_timer,      &
                                   WITH_LEVEL, TIME_ACCUMULATE_NO,     &
                                   TIME_ACCUMULATE_YES
    use pseudopotential_common,      only: pseudopotential

    implicit none

    ! Passed variables
    integer :: n_atoms, nsize
    integer :: inode, ionode
    real(double), dimension(:,:) :: density_out
    real(double), dimension(:,:) :: HF_force

    ! Local variables
    integer        :: i, j, my_block, n, the_species, iatom, spin, spin_2
    integer        :: ix, iy, iz, iblock, ipoint, igrid, stat, dir1, dir2
    integer        :: ipart, jpart, ind_part, ia, ii, icover, ig_atom
    real(double)   :: derivative, h_energy, rx, ry, rz, r2, r_from_i,  &
                      x, y, z, step
    real(double)   :: dcellx_block, dcelly_block, dcellz_block
    real(double)   :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double)   :: xatom, yatom, zatom
    real(double)   :: xblock, yblock, zblock
    real(double)   :: dx, dy, dz, loc_cutoff, loc_cutoff2, v
    real(double)   :: pcc_cutoff, pcc_cutoff2, step_pcc, x_pcc, y_pcc, &
                      z_pcc, derivative_pcc, v_pcc, jacobian
    real(double), dimension(3)       :: r, r_1, r_pcc
    real(double), dimension(3,nspin) :: fr_1, fr_pcc
    logical        :: range_flag
    type(cq_timer) :: tmr_l_tmp1, tmr_l_tmp2
    type(cq_timer) :: backtrace_timer

    real(double), dimension(3,3) :: loc_stress
    real(double), dimension(:),     allocatable :: h_potential,   &
                                                   density_total, &
                                                   density_out_total
    real(double), dimension(:,:,:), allocatable :: dVxc_drho
    real(double), dimension(nspin) :: pot_here, pot_here_pcc
    ! only for GGA with P.C.C.
    real(double), allocatable, dimension(:)   :: h_potential_in,       &
                                                 wk_grid_total,        &
                                                 density_out_GGA_total
    real(double), allocatable, dimension(:,:) :: wk_grid,              &
                                                 density_out_GGA    


!****lat<$
    call start_backtrace(t=backtrace_timer,who='get_nonSC_correction_force',where=7,level=3,echo=.true.)
!****lat>$ 
    ! Spin-polarised PBE non-SCF forces not implemented, so exit if necessary
    if ((nspin == 2) .and. flag_is_GGA) then ! Only true for CQ not LibXC
       if (inode == ionode) then
          write (io_lun, fmt='(10x,a)') &
               "*****************************************************"
          write (io_lun, fmt='(10x,a)') &
               "**** WARNING!!! WARNING!!! WARNING!!! WARNING!!! ****"
          write (io_lun, fmt='(10x,a)') &
               "**** non-SC correction forces are not implemented ***"
          write (io_lun, fmt='(10x,a)') &
               "**** for spin polarised version of PBE functionals **"
          write (io_lun, fmt='(10x,a)') &
               "**** correction forces will be set to ZERO !!!   ****"
          write (io_lun, fmt='(10x,a)') &
               "**** The forces are NOT to be TRUSTED !!!        ****"
          write (io_lun, fmt='(10x,a)') &
               "*****************************************************"
       end if
       HF_force = zero
       return
    end if

    stat = 0
    call start_timer(tmr_std_allocation)
    allocate(h_potential(nsize),           STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_nonSC_correction_force: Error alloc mem: ", nsize)
    allocate(density_total(nsize),         STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_nonSC_correction_force: Error alloc mem: ", nsize)
    allocate(density_out_total(nsize),     STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_nonSC_correction_force: Error alloc mem: ", nsize)
    allocate(dVxc_drho(nsize,nspin,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_nonSC_correction_force: Error alloc mem: ", nsize)
    call reg_alloc_mem(area_moveatoms, (3+nspin*nspin)*nsize, type_dbl)
    call stop_timer(tmr_std_allocation)

!****lat<$
    !print*, size(density_total,    dim=1)
    !print*, size(dVxc_drho,        dim=1)
    !print*, size(density_out_total,dim=1)
    !print*, size(potential,        dim=1)
    !print*, size(density_out_total,dim=1)
!****lat>$

    HF_force = zero

    dcellx_block = rcellx / blocks%ngcellx
    dcelly_block = rcelly / blocks%ngcelly
    dcellz_block = rcellz / blocks%ngcellz

    dcellx_grid = dcellx_block / nx_in_block
    dcelly_grid = dcelly_block / ny_in_block
    dcellz_grid = dcellz_block / nz_in_block

    ! We need a potential-like array containing the appropriate
    ! differences (between VHa for input and output charge densities
    ! and between Vxc' for input and output densities) Something like
    ! this:

    ! First the PAD n(r), density
    call start_timer(tmr_l_tmp2)

    h_potential = zero
    if(nspin==1) then
       density_total(:)     = spin_factor * density(:,1)
       density_out_total(:) = spin_factor * density_out(:,1)
    else
       density_total     = spin_factor * sum(density,    nspin)
       density_out_total = spin_factor * sum(density_out, nspin)
    end if

    ! Find the Hartree stress (loc_stress) of the output charge in the input potential for Jacobian
    call hartree(density_total, h_potential, nsize, h_energy,density_out_total,loc_stress)
    jacobian = zero
    ! use density_total (input charge) WITHOUT a factor of half so that this term corrects the GPV
    ! found in the main force routine
    do ipoint = 1,nsize
       jacobian = jacobian + (density_out_total(ipoint) - density_total(ipoint))* &
            (h_potential(ipoint) + pseudopotential(ipoint)) ! Also calculate the correction to GPV for local potential
    end do
    ! Don't apply gsum to jacobian because nonSCF_stress will be summed (end of this routine)
    jacobian = jacobian*grid_point_volume
    ! loc_stress is \int n^{out} V^{PAD}_{Har}
    ! NB The Hartree routine does NOT apply a factor of 0.5 to the second stress !
    ! Hartree_stress is 0.5 \int n^{PAD} V^{PAD}_{Har}
    ! Hartree_stress and loc_stress are the terms coming from the change of reciprocal space
    ! lattice vectors, 2GaGb/G^4
    ! Writing it this way gives out (n^{out} - 0.5 n^{PAD})*V^{PAD} as required
    if (flag_stress) then
      Hartree_stress(:,:) = loc_stress(:,:) - Hartree_stress(:,:)
      nonSCF_stress(1,1) = jacobian 
      nonSCF_stress(2,2) = jacobian 
      nonSCF_stress(3,3) = jacobian
    end if
    ! Find the PAD XC potential to calculate the correct Jacobian
    ! The correct Jacobian comes from \int n^{out} V_{XC}(n^{PAD}) + DeltaXC
    ! DeltaXC is added in the main force routine
    ! For PCC we will do this in the PCC force routine (easier)
    if (.NOT.flag_pcc_global) then
       call get_xc_potential(density, dVxc_drho(:,:,1),    &
               potential(:,1), y_pcc, nsize)
       jacobian = zero
       do spin = 1, nspin
          do ipoint = 1,nsize
             jacobian = jacobian + spin_factor*density_out(ipoint,spin)*dVxc_drho(ipoint,spin,1)
          end do
       end do
       jacobian = jacobian*grid_point_volume
       call gsum(jacobian) ! gsum as XC_stress isn't summed elsewhere
       ! Correct XC stress 
       if (flag_stress) then
         XC_stress(1,1) = XC_stress(1,1) + jacobian
         XC_stress(2,2) = XC_stress(2,2) + jacobian
         XC_stress(3,3) = XC_stress(3,3) + jacobian
       end if
    end if
    ! Accumulate PAD Hartree potential
    potential = zero
    do spin = 1, nspin
       call axpy(nsize, one, h_potential, 1, potential(:,spin), 1)
    end do
    dVxc_drho = zero
    call stop_timer(tmr_l_tmp2, TIME_ACCUMULATE_NO)

    ! for P.C.C.
    if (flag_pcc_global) then
       allocate(wk_grid_total(nsize), wk_grid(nsize,nspin), STAT=stat)
       wk_grid_total = zero
       wk_grid       = zero
       if (stat /= 0) &
            call cq_abort('Error allocating wk_grids in &
                           &get_nonSC_correction ', stat)
       call reg_alloc_mem(area_moveatoms, (nspin + 1) * nsize, type_dbl)
       do spin = 1, nspin
          wk_grid(:,spin)  = density(:,spin)  + half * density_pcc(:)
          wk_grid_total(:) = wk_grid_total(:) + spin_factor * wk_grid(:,spin)
       end do
       ! only for GGA
       if (flag_is_GGA) then
          allocate(density_out_GGA_total(nsize), density_out_GGA(nsize,nspin), STAT=stat)
          if (stat /= 0) call cq_abort ('Error allocating &
                               &density_out_GGAs in get_nonSC_force ', stat)
          call reg_alloc_mem(area_moveatoms, (nspin + 1) * nsize, type_dbl)
          density_out_GGA_total = zero
          density_out_GGA       = zero
          do spin = 1, nspin
             density_out_GGA(:,spin)  = density_out(:,spin)      + half * density_pcc(:)
             density_out_GGA_total(:) = density_out_GGA_total(:) + spin_factor * density_out_GGA(:,spin)
          end do
          ! copy hartree potential
          allocate(h_potential_in(nsize), STAT=stat)
          if (stat /= 0) &
               call cq_abort('Error allocating h_potential_in in &
                              &get_nonSC_force ', stat)
          call reg_alloc_mem(area_moveatoms, nsize, type_dbl)
          h_potential_in = h_potential
       end if
    end if

    call start_timer (tmr_l_tmp1, WITH_LEVEL)
    if (flag_pcc_global) then
       if(flag_is_GGA) then
          call get_dxc_potential(wk_grid, dVxc_drho, nsize, density_out_GGA)
          ! GGA with spin not implemented ! 
          potential(:,1) = potential(:,1) + dVxc_drho(:,1,1) 
       else
          call get_dxc_potential(wk_grid, dVxc_drho, nsize)
       end if
    else
       if(flag_is_GGA) then
          call get_dxc_potential(density, dVxc_drho, nsize, density_out)
          ! GGA with spin not implemented ! 
          potential(:,1) = potential(:,1) + dVxc_drho(:,1,1) 
       else
          call get_dxc_potential(density, dVxc_drho, nsize)
       end if
    end if
    ! deallocating density_out_GGA: only for P.C.C.
    if (flag_pcc_global)  then
       if (flag_is_GGA) then
          deallocate(density_out_GGA_total, density_out_GGA, STAT=stat)
          if (stat /= 0) call cq_abort('Error deallocating density_out_GGAs in &
                              &get_nonSC_force ', stat)
          call reg_dealloc_mem(area_moveatoms, (nspin + 1) * nsize, type_dbl)
          ! make a copy of potential at this point
          ! use wk_grid as a temporary storage
          do spin = 1, nspin
             wk_grid(:,spin) = potential(:,spin)
          end do
       end if
    end if !flag_pcc_global

    ! for LDA
    if (.NOT.(flag_is_GGA)) then
       do spin = 1, nspin
          do spin_2 = 1, nspin
             do i = 1, n_my_grid_points
                ! contribution from (n^PAD - n^out) * mu'(n^PAD)
                potential(i,spin) =                                &
                     potential(i,spin) + spin_factor *             &
                     (density(i,spin_2) - density_out(i,spin_2)) * &
                     dVxc_drho(i,spin_2,spin)
             end do ! i
          end do ! spin_2
       end do ! spin
    end if ! for LDA

    call stop_print_timer(tmr_l_tmp1, "NSC force - XC", IPRINT_TIME_THRES3)

    ! Restart of the timer; level assigned here
    call start_timer(tmr_l_tmp2, WITH_LEVEL)
    h_potential = zero
    ! Preserve the Hartree stress we've calculated
    loc_stress = Hartree_stress
    call hartree(density_out_total, h_potential, nsize, h_energy)
    ! And restore
    Hartree_stress = loc_stress
    do spin = 1, nspin
       call axpy(nsize, -one, h_potential, 1, potential(:,spin), 1)
    end do
    ! now, potential = - delta V_H - ( n_out - n_PAD ) * dxc_potential
    call stop_print_timer(tmr_l_tmp2, "NSC force - Hartree", &
                          IPRINT_TIME_THRES3, TIME_ACCUMULATE_YES)
    call my_barrier()
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock = (domain%idisp_primx(iblock) + domain%nx_origin - 1) * dcellx_block
       yblock = (domain%idisp_primy(iblock) + domain%ny_origin - 1) * dcelly_block
       zblock = (domain%idisp_primz(iblock) + domain%nz_origin - 1) * dcellz_block
       if (naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom = 0
          do ipart = 1, naba_atoms_of_blocks(dens)%no_of_part(iblock)
             jpart = naba_atoms_of_blocks(dens)%list_part(ipart,iblock)
             if (jpart > DCS_parts%mx_gcover)&
                  call cq_abort('set_ps: JPART ERROR ', ipart, jpart)
             ind_part = DCS_parts%lab_cell(jpart)
             do ia = 1, naba_atoms_of_blocks(dens)%no_atom_on_part(ipart, iblock)
                iatom = iatom + 1
                ii = naba_atoms_of_blocks(dens)%list_atom(iatom, iblock)
                icover = DCS_parts%icover_ibeg(jpart) + ii - 1
                ig_atom = id_glob(parts%icell_beg(ind_part) + ii - 1)
                if (parts%icell_beg(ind_part) + ii - 1 > ni_in_cell) &
                     call cq_abort('set_ps: globID ERROR ', &
                                   ii, parts%icell_beg(ind_part))
                if (icover > DCS_parts%mx_mcover) &
                     call cq_abort ('set_ps: icover ERROR ', &
                                    icover, DCS_parts%mx_mcover)
                xatom = DCS_parts%xcover(icover)
                yatom = DCS_parts%ycover(icover)
                zatom = DCS_parts%zcover(icover)
                the_species = species_glob(ig_atom)
                !the_species=species(ig_atom)
                loc_cutoff = atomic_density_table(the_species)%cutoff
                loc_cutoff2 = loc_cutoff*loc_cutoff
                step = loc_cutoff / &
                       real(atomic_density_table(the_species)%length - 1, double)
                ipoint = 0
                do iz = 1, nz_in_block
                   do iy = 1, ny_in_block
                      do ix = 1, nx_in_block
                         ipoint = ipoint + 1
                         igrid = n_pts_in_block * (iblock - 1) + ipoint
                         if (igrid > n_my_grid_points) &
                              call cq_abort('get_nonSC_force: igrid error ',&
                                            igrid,n_my_grid_points)
                         ! scale with integration measure
                         do spin = 1, nspin
                            pot_here(spin) = &
                                 potential(igrid,spin) * grid_point_volume
                         end do
                         dx = dcellx_grid * (ix - 1)
                         dy = dcelly_grid * (iy - 1)
                         dz = dcellz_grid * (iz - 1)
                         r(1) = xblock + dx - xatom
                         r(2) = yblock + dy - yatom
                         r(3) = zblock + dz - zatom
                         r2 = rx * rx + ry * ry + rz * rz
                         if (r2 < loc_cutoff2) then
                            r_from_i = sqrt(r2)
                            if (r_from_i > RD_ERR) then
                               r_1(1) = r(1) / r_from_i
                               r_1(2) = r(2) / r_from_i
                               r_1(3) = r(3) / r_from_i
                            else
                               r_1 = zero
                            end if
                            call dsplint(                                       &
                                 step,                                          &
                                 atomic_density_table(the_species)%table(:),    &
                                 atomic_density_table(the_species)%d2_table(:), &
                                 atomic_density_table(the_species)%length,      &
                                 r_from_i, v, derivative, range_flag)
                            if (range_flag) &
                                 call cq_abort('get_nonSC_force: overrun problem')
                            ! We assumed the atomic densities were evenly devided in spin channels at
                            ! start, (in set_density of density module). So we assume the same to be
                            ! consistent, and then apply density_scale calculated from set_density
                            ! NB This means that spin_factor cancels out the half for non-spin polarised
                            do dir1=1,3
                              do spin = 1, nspin
                                 fr_1(dir1,spin) = -r_1(dir1) * half * &
                                   derivative * density_scale(spin)
                              end do
                            end do
                         else
                            fr_1 = zero
                         end if
                         ! could be written in a simpler form, but written this way gives more clear idea
                         ! on what we are doing here.
                         do spin = 1, nspin
                            do dir1 = 1, 3
                                HF_force(dir1,ig_atom) = &
                                  HF_force(dir1,ig_atom) + spin_factor * &
                                  fr_1(dir1,spin) * pot_here(spin)
                              if (flag_stress) then
                                if (flag_full_stress) then
                                  do dir2 = 1, 3
                                    nonSCF_stress(dir1,dir2) = &
                                      nonSCF_stress(dir1,dir2) + r(dir1) * &
                                      spin_factor * fr_1(dir2,spin) * &
                                      pot_here(spin)
                                    if (flag_atomic_stress) then
                                      atomic_stress(dir1,dir2,ig_atom) = &
                                        atomic_stress(dir1,dir2,ig_atom) + &
                                        r(dir1) * spin_factor * &
                                        fr_1(dir2,spin) * pot_here(spin)
                                    end if ! flag_atomic_stress
                                  end do ! dir2
                                else
                                  nonSCF_stress(dir1,dir1) = &
                                    nonSCF_stress(dir1,dir1) + r(dir1) * &
                                    spin_factor * fr_1(dir1,spin) * &
                                    pot_here(spin)
                                end if ! flag_full_stress
                              end if ! flag_stress
                            end do ! dir1
                         end do ! spin
                      end do !ix
                   end do  !iy
                end do   !iz
             end do ! naba_atoms
          end do ! naba_part
       end if !(naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) !naba atoms?
    end do ! iblock : primary set of blocks

    call stop_print_timer(tmr_l_tmp1, "NSC force - Orbital part", &
                          IPRINT_TIME_THRES3)

    ! only called for P.C.C.
    ! compute - int d^3r ( delta n_{v} * dxc(n_{c} + n_{v} ) * dn_{c} )
    if (flag_pcc_global) then
       if (flag_is_GGA) then
          ! for GGA
          potential = zero
          do spin = 1, nspin
             do i = 1, n_my_grid_points
                ! -delta n * dxc_potential
                potential(i,spin) = wk_grid(i,spin) - h_potential_in(i)
             end do
          end do
       else
          ! For LDA
          potential = zero
          do spin = 1, nspin
             do spin_2 = 1, nspin
                do i = 1, n_my_grid_points
                   potential(i,spin) =                                &
                        potential(i,spin) +                           &
                        spin_factor *                                 &
                        (density(i,spin_2) - density_out(i,spin_2)) * &
                        dVxc_drho(i,spin_2,spin)
                end do
             end do
          end do
       end if

       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       do iblock = 1, domain%groups_on_node ! primary set of blocks
          xblock = (domain%idisp_primx(iblock) + domain%nx_origin - 1) * &
                   dcellx_block
          yblock = (domain%idisp_primy(iblock) + domain%ny_origin - 1) * &
                   dcelly_block
          zblock = (domain%idisp_primz(iblock) + domain%nz_origin - 1) * &
                   dcellz_block
          if (naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) then ! if there are naba atoms
             iatom = 0
             do ipart = 1, naba_atoms_of_blocks(dens)%no_of_part(iblock)
                jpart = naba_atoms_of_blocks(dens)%list_part(ipart, iblock)
                if (jpart > DCS_parts%mx_gcover) &
                     call cq_abort('set_ps: JPART ERROR ', ipart, jpart)
                ind_part = DCS_parts%lab_cell(jpart)
                do ia = 1, naba_atoms_of_blocks(dens)%no_atom_on_part(ipart, iblock)
                   iatom = iatom + 1
                   ii = naba_atoms_of_blocks(dens)%list_atom(iatom, iblock)
                   icover = DCS_parts%icover_ibeg(jpart) + ii - 1
                   ig_atom = id_glob(parts%icell_beg(ind_part) + ii - 1)
                   the_species = species_glob(ig_atom)
                   ! for P.C.C. treatment
                   if (.not. pseudo(the_species)%flag_pcc) cycle
                   if (parts%icell_beg(ind_part) + ii - 1 > ni_in_cell) &
                        call cq_abort('set_ps: globID ERROR ', &
                                      ii, parts%icell_beg(ind_part))
                   if (icover > DCS_parts%mx_mcover) &
                        call cq_abort('set_ps: icover ERROR ', &
                                      icover, DCS_parts%mx_mcover)
                   xatom = DCS_parts%xcover(icover)
                   yatom = DCS_parts%ycover(icover)
                   zatom = DCS_parts%zcover(icover)
                   pcc_cutoff = pseudo(the_species)%chpcc%cutoff
                   pcc_cutoff2 = pcc_cutoff * pcc_cutoff
                   step_pcc = pseudo(the_species)%chpcc%delta
                   ipoint = 0
                   do iz = 1, nz_in_block
                      do iy = 1, ny_in_block
                         do ix = 1, nx_in_block
                            ipoint = ipoint+1
                            igrid = n_pts_in_block * (iblock - 1) + ipoint
                            if (igrid > n_my_grid_points) &
                                 call cq_abort('get_nonSC_force: igrid error ',&
                                               igrid, n_my_grid_points)
                            do spin = 1, nspin
                               pot_here_pcc(spin) = &
                                    potential(igrid,spin) * grid_point_volume
                            end do
                            dx = dcellx_grid * (ix - 1)
                            dy = dcelly_grid * (iy - 1)
                            dz = dcellz_grid * (iz - 1)
                            r(1) = xblock + dx - xatom
                            r(2) = yblock + dy - yatom
                            r(3) = zblock + dz - zatom
                            r2 = rx * rx + ry * ry + rz * rz
                            if (r2 < pcc_cutoff2) then
                               r_from_i = sqrt( r2 )
                               if ( r_from_i > RD_ERR ) then
                                  r_pcc(1) = r(1) / r_from_i
                                  r_pcc(2) = r(2) / r_from_i
                                  r_pcc(3) = r(3) / r_from_i
                               else
                                  r_pcc = zero
                               end if
                               call dsplint(step_pcc, &
                                            pseudo(the_species)%chpcc%f(:), &
                                            pseudo(the_species)%chpcc%d2(:), &
                                            pseudo(the_species)%chpcc%n, &
                                            r_from_i, v_pcc, derivative_pcc, &
                                            range_flag)
                               if (range_flag) &
                                    call cq_abort('get_nonSC_force: &
                                                  &overrun problem')
                               ! We assumed the atomic densities were
                               ! evenly devided in spin channels at
                               ! start, (in set_density of density
                               ! module). So we assume the same to be
                               ! consistent, and then apply density_scale
                               ! calculated from set_density
                               do spin = 1, nspin
                                  do dir1 = 1, 3
                                    fr_pcc(dir1,spin) = r_pcc(dir1) * half * derivative_pcc * density_scale(spin)
                                  end do
                               end do
                            else
                               fr_pcc = zero
                            end if
                            ! assuming derivative of ppc charge
                            ! (assumed same for different spin
                            ! components)
                            do spin = 1, nspin
                               do dir1 = 1, 3
                                 HF_force(dir1,ig_atom) = &
                                   HF_force(dir1,ig_atom) + spin_factor * &
                                   fr_pcc(dir1,spin) * pot_here_pcc(spin)
                                 if (flag_stress) then
                                   if (flag_full_stress) then
                                     do dir2 = 1, 3
                                       nonSCF_stress(dir1,dir2) = &
                                         nonSCF_stress(dir1,dir2) + r(dir1) * &
                                         spin_factor * fr_pcc(dir2,spin) * &
                                         pot_here_pcc(spin)
                                       if (flag_atomic_stress) then
                                         atomic_stress(dir1,dir2,ig_atom) = &
                                           atomic_stress(dir1,dir2,ig_atom) + &
                                           r(dir1) * spin_factor * &
                                           fr_pcc(dir2,spin) * &
                                           pot_here_pcc(spin)
                                       end if ! flag_atomic_stress
                                     end do ! dir1
                                   else
                                     nonSCF_stress(dir1,dir1) = &
                                       nonSCF_stress(dir1,dir1) + r(dir1) * &
                                       spin_factor * fr_pcc(dir1,spin) * &
                                       pot_here_pcc(spin)
                                   end if ! flag_full_stress
                                 end if ! flag_stress
                               end do ! dir2
                            end do ! spin
                         end do !ix
                      end do  !iy
                   end do   !iz
                end do ! naba_atoms
             end do ! naba_part
          end if !(naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) !naba atoms?
       end do ! iblock : primary set of blocks
    end if ! (flag_pcc_global)

    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    call gsum(HF_force, 3, n_atoms)
    if (flag_stress) call gsum(nonSCF_stress,3,3)
    call stop_print_timer(tmr_l_tmp1, "NSC force - Compilation", &
                          IPRINT_TIME_THRES3)

    ! deallocating temporary arrays
    call start_timer(tmr_std_allocation)
    if (flag_pcc_global) then
       deallocate(wk_grid_total, wk_grid, STAT=stat)
       if (stat /= 0) &
            call cq_abort('Error deallocating wk_grid in &
                           &get_nonSC_correction_force ', stat)
       call reg_dealloc_mem(area_moveatoms, (nspin + 1) * nsize, type_dbl)
       if (flag_is_GGA) then
          deallocate(h_potential_in, STAT=stat)
          if (stat /= 0) &
               call cq_abort('Error deallocating h_potential_in in &
                             &get_nonSC_correction_force ', stat)
          call reg_dealloc_mem(area_moveatoms, nsize, type_dbl)
       end if ! for GGA
    end if ! flag_pcc_global
    deallocate(h_potential, density_total, density_out_total, dVxc_drho, &
               STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_nonSC_correction_force: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, (3+nspin*nspin)*nsize, type_dbl)
    call stop_timer(tmr_std_allocation)


!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_nonSC_correction_force',echo=.true.)
!****lat>$  

    return
  end subroutine get_nonSC_correction_force
  !!***


  ! -----------------------------------------------------------
  ! Subroutine get_pcc_force
  ! -----------------------------------------------------------

  !!****f* force_module/get_pcc_force *
  !!
  !!  NAME
  !!   get_pcc_force
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Gets the P.C.C. force in the case where the partial core
  !!   correction is taken into account.  This contribution works on
  !!   both SCF/NSC calculations.
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   M.Arita
  !!  CREATION DATE
  !!   12:03, 2011/03/31 M.Arita
  !!  MODIFICATION HISTORY
  !!   2011/10/12 L.Tong
  !!     Corrected memory allocation and deallocation registers
  !!     Added spin polarisation
  !!     - Assumed the pcc correction to density is the same for both
  !!       spin components, and hence the pcc correction will be half *
  !!       density_pcc for all spin
  !!     - xc_potential for spin polarised case is used to store the
  !!       total exchange-correlation potential. First the spin up
  !!       component is stored in xc_potential, then xc_potential =
  !!       xc_potential + xc_potential_dn
  !!   2012/03/25 L.Tong
  !!   - Changed spin implementation
  !!   2012/05/29 L.Tong
  !!   - Cleaned up xc-functonal selector. Now spin and non-spin
  !!     calculations share the same calls more or less.
  !!   2015/08/10 08:12 dave
  !!    Adding stress for PCC and non-SCF PCC XC term
  !!    N.B. this relies on the output density, which is passed as an optional argument
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2017/10/20 12:08 dave
  !!    Added extra optional argument to allow return of XC energy (for force testing)
  !!   2018/02/09 14:41 dave
  !!    Adding call for LibXC integration
  !!   2018/02/13 11:52 dave
  !!    Changes for new, universal XC interface
  !!   2018/02/14 13:26 dave
  !!    More subtle errors ! The PCC, non-SCF XC stress did not have spin
  !!    factor applied (as above in non-SCF routine)
  !!   2019/05/08 zamaan
  !!    Added atomic stress contributions
  !!  SOURCE
  !!
  subroutine get_pcc_force(pcc_force, inode, ionode, n_atoms, size, density_out,xc_energy_ret)

    use datatypes
    use numbers
    use species_module,      only: species
    use GenComms,            only: gsum
    use global_module,       only: rcellx, rcelly, rcellz, id_glob,    &
                                   ni_in_cell, species_glob, dens,     &
                                   area_moveatoms, IPRINT_TIME_THRES3, &
                                   nspin, spin_factor, flag_self_consistent, &
                                   flag_full_stress, flag_stress,      &
                                   flag_atomic_stress
    use block_module,        only: nx_in_block,ny_in_block,            &
                                   nz_in_block, n_pts_in_block
    use group_module,        only: blocks, parts
    use primary_module,      only: domain
    use cover_module,        only: DCS_parts
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use GenComms,            only: my_barrier, cq_abort
    use pseudo_tm_info,      only: pseudo
    use spline_module,       only: dsplint
    use dimens,              only: grid_point_volume, n_my_grid_points
    use GenBlas,             only: axpy
    use density_module,      only: density, density_scale, density_pcc
    use XC,                  only: get_xc_potential
    use maxima_module,       only: maxngrid
    use memory_module,       only: reg_alloc_mem, reg_dealloc_mem,     &
                                   type_dbl
    use timer_module,        only: cq_timer, start_timer, stop_timer,  &
                                   print_timer, stop_print_timer,      &
                                   WITH_LEVEL, TIME_ACCUMULATE_NO,     &
                                   TIME_ACCUMULATE_YES

    implicit none

    ! Passed variables
    integer :: n_atoms, size
    integer :: inode, ionode
    real(double), dimension(3,n_atoms) :: pcc_force
    real(double), dimension(:,:), OPTIONAL :: density_out
    real(double), OPTIONAL :: xc_energy_ret
    
    ! Local variables
    integer        :: i, j, my_block, n, the_species, iatom, spin
    integer        :: ix, iy, iz, iblock, ipoint, igrid, stat, dir1, dir2
    integer        :: ipart, jpart, ind_part, ia, ii, icover, ig_atom
    real(double)   :: derivative_pcc, xc_energy, r2,      &
                      r_from_i, x_pcc, y_pcc, z_pcc, step_pcc
    real(double)   :: dcellx_block, dcelly_block, dcellz_block
    real(double)   :: dcellx_grid, dcelly_grid, dcellz_grid
    real(double)   :: xatom, yatom, zatom
    real(double)   :: xblock, yblock, zblock
    real(double)   :: dx, dy, dz, pcc_cutoff, pcc_cutoff2, electrons, &
                      v_pcc, jacobian
    logical        :: range_flag
    type(cq_timer) :: tmr_l_tmp1, tmr_l_tmp2
    ! automatic arrays
    real(double), dimension(nspin) :: pot_here_pcc
    real(double), dimension(3)   :: r_pcc, fr_pcc, r
    ! allocatable arrays
    real(double), dimension(:),   allocatable :: xc_epsilon, density_wk_tot
    real(double), dimension(:,:), allocatable :: xc_potential, density_wk
    type(cq_timer) :: backtrace_timer

    call start_backtrace(t=backtrace_timer,who='get_PCC_force',where=7,level=3,echo=.true.)
    allocate(xc_epsilon(size), density_wk_tot(size), &
             xc_potential(size,nspin), density_wk(size,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("get_pcc_force: Error alloc mem: ", size)
    call reg_alloc_mem(area_moveatoms, (2+2*nspin)*size, type_dbl)

    ! initialise arrays
    pcc_force = zero
    xc_potential = zero
    xc_epsilon = zero

    dcellx_block = rcellx / blocks%ngcellx
    dcelly_block = rcelly / blocks%ngcelly
    dcellz_block = rcellz / blocks%ngcellz

    dcellx_grid = dcellx_block / nx_in_block
    dcelly_grid = dcelly_block / ny_in_block
    dcellz_grid = dcellz_block / nz_in_block

    call start_timer (tmr_l_tmp2)
    density_wk = zero
    density_wk_tot = zero
    do spin = 1, nspin
       density_wk(:,spin) = density(:,spin) + half * density_pcc(:)
       density_wk_tot(:) = density_wk_tot(:) + spin_factor * density_wk(:,spin)
    end do

    call get_xc_potential(density_wk, xc_potential,     &
         xc_epsilon, xc_energy, size)
    if(PRESENT(xc_energy_ret)) xc_energy_ret = xc_energy
    ! We do this here to re-use xc_potential - for non-PCC we do it in get_nonSC_correction_force
    if(.NOT.flag_self_consistent) then
       if(.NOT.present(density_out)) call cq_abort("Output density not passed to PCC force for nonSCF calculation")
       if (flag_stress) then
         jacobian = zero
         do spin=1,nspin
            do ipoint = 1,size
               jacobian = jacobian + spin_factor*density_out(ipoint,spin)*xc_potential(ipoint,spin)
            end do
         end do
         jacobian = jacobian*grid_point_volume
         call gsum(jacobian) ! gsum as XC_stress isn't summed elsewhere
         ! Correct XC stress 
         XC_stress(1,1) = XC_stress(1,1) + jacobian
         XC_stress(2,2) = XC_stress(2,2) + jacobian
         XC_stress(3,3) = XC_stress(3,3) + jacobian
       end if
    end if
    ! This restarts the count for this timer
    call stop_timer(tmr_l_tmp2, TIME_ACCUMULATE_NO)

    call my_barrier()

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    do iblock = 1, domain%groups_on_node ! primary set of blocks
       xblock = (domain%idisp_primx(iblock) + domain%nx_origin - 1) * dcellx_block
       yblock = (domain%idisp_primy(iblock) + domain%ny_origin - 1) * dcelly_block
       zblock = (domain%idisp_primz(iblock) + domain%nz_origin - 1) * dcellz_block
       if (naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) then ! if there are naba atoms
          iatom = 0
          do ipart = 1, naba_atoms_of_blocks(dens)%no_of_part(iblock)
             jpart = naba_atoms_of_blocks(dens)%list_part(ipart,iblock)
             if (jpart > DCS_parts%mx_gcover) &
                  call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
             ind_part = DCS_parts%lab_cell(jpart)
             do ia = 1, naba_atoms_of_blocks(dens)%no_atom_on_part(ipart,iblock)
                iatom = iatom + 1
                ii = naba_atoms_of_blocks(dens)%list_atom(iatom,iblock)
                icover = DCS_parts%icover_ibeg(jpart) + ii - 1
                ig_atom = id_glob(parts%icell_beg(ind_part) + ii - 1)
                the_species = species_glob(ig_atom)
                !the_species=species(ig_atom)
                ! for P.C.C. treatment
                if (.not. pseudo(the_species)%flag_pcc) cycle
                if (parts%icell_beg(ind_part) + ii-1 > ni_in_cell) &
                     call cq_abort('set_ps: globID ERROR ', &
                                   ii, parts%icell_beg(ind_part))
                if (icover > DCS_parts%mx_mcover) &
                     call cq_abort('set_ps: icover ERROR ', &
                                   icover, DCS_parts%mx_mcover)
                xatom = DCS_parts%xcover(icover)
                yatom = DCS_parts%ycover(icover)
                zatom = DCS_parts%zcover(icover)
                pcc_cutoff = pseudo(the_species)%chpcc%cutoff
                pcc_cutoff2 = pcc_cutoff * pcc_cutoff
                step_pcc = pseudo(the_species)%chpcc%delta
                ipoint = 0
                do iz = 1, nz_in_block
                   do iy = 1, ny_in_block
                      do ix = 1, nx_in_block
                         ipoint = ipoint + 1
                         igrid = n_pts_in_block*(iblock-1)+ipoint
                         if (igrid > n_my_grid_points) &
                              call cq_abort('get_nonSC_force: igrid error ', &
                                            igrid, n_my_grid_points)
                         do spin = 1, nspin
                            pot_here_pcc(spin) = &
                                 xc_potential(igrid,spin) * grid_point_volume
                         end do
                         dx = dcellx_grid * (ix - 1)
                         dy = dcelly_grid * (iy - 1)
                         dz = dcellz_grid * (iz - 1)
                         r(1) = xblock + dx - xatom
                         r(2) = yblock + dy - yatom
                         r(3) = zblock + dz - zatom
                         r2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
                         if (r2 < pcc_cutoff2) then
                            r_from_i = sqrt(r2)
                            if (r_from_i > RD_ERR) then
                               do dir1=1,3
                                  r_pcc(dir1) = r(dir1)/r_from_i
                               end do
                            else
                               r_pcc = zero
                            end if
                            call dsplint(step_pcc,                        &
                                         pseudo(the_species)%chpcc%f(:),  &
                                         pseudo(the_species)%chpcc%d2(:), &
                                         pseudo(the_species)%chpcc%n,     &
                                         r_from_i, v_pcc, derivative_pcc, &
                                         range_flag)
                            if (range_flag) &
                                 call cq_abort('get_pcc_force: overrun problem')
                            ! the factor of half here is because for
                            ! spin polarised calculations, I have
                            ! assumed contribution from pcc_density is
                            ! exactly half of the total in each spin
                            ! channel.
                            derivative_pcc = half * derivative_pcc
                            do dir1=1,3
                               fr_pcc(dir1) = r_pcc(dir1)*derivative_pcc
                            end do
                         else
                            fr_pcc = zero
                         end if
                         do spin = 1, nspin
                            do dir1=1,3
                               pcc_force(dir1,ig_atom) = &
                                 pcc_force(dir1,ig_atom) + spin_factor * &
                                 fr_pcc(dir1) * pot_here_pcc(spin)
                               if (flag_stress) then
                                 if (flag_full_stress) then
                                   do dir2=1,3
                                     pcc_stress(dir1,dir2) = &
                                       pcc_stress(dir1,dir2) + &
                                       r(dir1) * spin_factor * fr_pcc(dir2) * &
                                       pot_here_pcc(spin)
                                     if (flag_atomic_stress) then
                                       atomic_stress(dir1,dir2,ig_atom) = &
                                         atomic_stress(dir1,dir2,ig_atom) + &
                                         r(dir1) * spin_factor * &
                                         fr_pcc(dir2)  * pot_here_pcc(spin)
                                     end if
                                   end do ! dir2
                                 else
                                    pcc_stress(dir1,dir1) = &
                                      pcc_stress(dir1,dir1) + &
                                      r(dir1) * spin_factor * fr_pcc(dir1) * &
                                      pot_here_pcc(spin)
                                end if
                              end if
                           end do ! dir1
                         end do ! spin
                      end do !ix
                   end do  !iy
                end do   !iz
             end do ! naba_atoms
          end do ! naba_part
       end if !(naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) !naba atoms?
    end do ! iblock : primary set of blocks
    call stop_print_timer(tmr_l_tmp1, "PCC force - Orbital part", &
                          IPRINT_TIME_THRES3)
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    call gsum(pcc_force, 3, n_atoms)
    if (flag_stress) call gsum(pcc_stress, 3, 3)
    call stop_print_timer(tmr_l_tmp1, "PCC force - Compilation", &
                          IPRINT_TIME_THRES3)

    deallocate(xc_epsilon, density_wk_tot, xc_potential, density_wk, STAT=stat)
    if (stat /= 0) call cq_abort("get_pcc_force: Error dealloc mem")
    call reg_dealloc_mem(area_moveatoms, (2+2*nspin)*size, type_dbl)
    call stop_backtrace(t=backtrace_timer,who='get_PCC_force',echo=.true.)

    return
  end subroutine get_pcc_force
  !*****

!!****f* force_module/print_stress *
!!
!!  NAME 
!!   print_stress
!!  PURPOSE
!!   Print the stress or stress contribution, format depending on print 
!!   level and whether the full tensor was computed
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   28 March 2019
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
subroutine print_stress(label, str_mat, print_level)

  use datatypes
  use numbers
  use units
  use GenComms,       only: inode, ionode
  use global_module,  only: iprint_MD, flag_full_stress

  ! Passed variables
  character(*), intent(in)                  :: label
  real(double), dimension(3,3), intent(in)  :: str_mat
  integer, intent(in)                       :: print_level

  ! local variables
  character(20) :: fmt = '(4x,a18,3f15.8,a3)'
  character(20) :: blank = ''

  if (inode==ionode) then
    if (iprint_MD >= print_level) then
      if (flag_full_stress) then
        write(io_lun,fmt=fmt) label, str_mat(1,:), en_units(energy_units)
        write(io_lun,fmt=fmt) blank, str_mat(2,:), blank
        write(io_lun,fmt=fmt) blank, str_mat(3,:), blank
      else
        write(io_lun,fmt=fmt) label, str_mat(1,1), str_mat(2,2), &
                              str_mat(3,3), en_units(energy_units)
      end if
    end if
  end if
end subroutine print_stress
!!*****

!!****f* force_module/check_atomic_stress *
!!
!!  NAME 
!!   check_atomic_stress
!!  PURPOSE
!!   Sanity check the atomic stress contributions
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   28 March 2019
!!  MODIFICATION HISTORY
!!  
!!  SOURCE
!!  
subroutine check_atomic_stress

  use GenComms,       only: inode, ionode
  use global_module,  only: iprint_MD, flag_full_stress, atomic_stress, &
                            non_atomic_stress

  ! Passed variables

  ! local variables
  real(double), dimension(3,3) :: total_atomic, total
  integer :: dir1, dir2

  do dir1=1,3
    do dir2=1,3
      total_atomic(dir1,dir2) = sum(atomic_stress(dir1,dir2,:))
      total(dir1,dir2)  = total_atomic(dir1,dir2) + &
                          non_atomic_stress(dir1,dir2)
    end do
  end do

  if (inode==ionode) then
    write(io_lun,'(2x,a)') &
      "Checking sum of atomic and non-atomic contributions to stress:"
  end if
  call print_stress("Atomic total:     ", total_atomic, 3)
  call print_stress("Non-atomic total: ", non_atomic_stress, 3)
  call print_stress("Total:            ", total, 3)

end subroutine check_atomic_stress
!!*****

end module force_module
