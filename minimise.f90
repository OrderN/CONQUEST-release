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
!!   2015/06/10 15:48 cor & dave
!!    Wavefunction output
!!  SOURCE
!!
module minimise

  use datatypes
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_stdclocks_module, only: tmr_std_eminimisation

  implicit none

  save

  logical      :: UsePulay
  integer      :: n_L_iterations, n_support_iterations
  real(double) :: L_tolerance, sc_tolerance, energy_tolerance, &
                  expected_reduction

  ! Area identification
  integer, parameter, private :: area = 6

  ! RCS tag for object file identification
  character(len=80), private :: &
       RCSid = "$Id$"
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
  !!    Changed call to minimise to pick and choose between full run
  !!    (varying DM, self-consistency and blips), self-consistent run
  !!    (varying DM, self-consistency) and simple TB (varying DM only)
  !!   14:39, 26/02/2003 drb
  !!    Added call to get_H_matrix to build new H and charge
  !!   08:35, 2003/03/12 dave
  !!    Added call to get_energy after FindMinDM
  !!   08:30, 2003/10/01 dave
  !!    Changed flag_vary_blips to basis and added basis switch for
  !!    outer loop minimiser
  !!   2006/11/13 18:20 dave
  !!    Added new flags to test for finding and/or writing forces
  !!   2007/04/17 09:34 dave
  !!    Passed no. of L iterations to vary_support
  !!   2008/05/25 ast
  !!    Added timer
  !!   2011/12/07 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2012/03/26 L.Tong
  !!   - Removed redundant input parameter real(double) mu
  !!   2012/04/29 L.Tong
  !!   - Added calculations for vdWDF xc-energy correction
  !!   2012/05/29 L.Tong
  !!   - Added timer for vdWDFT energy correction
  !!   2013/01/30 10:32 dave
  !!   - Added call for deltaSCF (with U. Terranova)
  !!   2013/08/20 M.Arita
  !!   -  Added 'iter' as a dummy variable
  !!   -  Added grequency to go back to McW
  !!   -  Added call for writing out L-matrix (this will be deleted later)
  !!   2013/12/03 M.Arita
  !!   - Removed calls for writing out L-matrix
  !!   2015/06/08 lat
  !!   - Added experimental backtrace
  !!   2015/06/10 15:48 cor & dave
  !!    Added call for wavefunction output at self-consistency
  !!   2015/06/26 13:46 dave
  !!    Turned off force calculation when writing out bands
  !!   2017/02/15 (or earlier) nakata
  !!    Added LFD minimisation for multisite support functions
  !!   2017/11/10 dave
  !!    Removed calls to dump K matrix (now done in DMMinModule)
  !!  SOURCE
  !!
  !subroutine get_E_and_F(fixed_potential, vary_mu, total_energy, &
  !                       find_forces, write_forces)
  subroutine get_E_and_F(fixed_potential, vary_mu, total_energy, &
                         find_forces, write_forces, iter, level)

    use datatypes
    use force_module,      only: force
    use DMMin,             only: FindMinDM
    use SelfCon,           only: new_SC_potl, atomch_output,           &
                                 get_atomic_charge
    use global_module,     only: flag_vary_basis,                      &
                                 flag_self_consistent, flag_basis_set, &
                                 blips, PAOs, IPRINT_TIME_THRES1,      &
                                 runtype, flag_vdWDFT, io_lun,         &
                                 flag_DeltaSCF, flag_excite, runtype,  &
                                 flag_LmatrixReuse,McWFreq,            &
                                 flag_multisite,                       &
                                 io_lun, flag_out_wf, wf_self_con, flag_write_DOS, &
                                 flag_diagonalisation, nspin
    use energy,            only: get_energy, xc_energy, final_energy
    use GenComms,          only: cq_abort, inode, ionode
    use blip_minimisation, only: vary_support
    use pao_minimisation,  only: vary_pao, pulay_min_pao
    use timer_module
    use input_module,      only: leqi
    use vdWDFT_module,     only: vdWXC_energy, vdWXC_energy_slow
    use density_module,    only: density
    use multisiteSF_module,only: flag_LFD_minimise, LFD_minimise
    use units

    implicit none

    ! Shared variables
    integer, intent(in), optional :: iter
    integer, intent(in), optional :: level
    logical           :: vary_mu, fixed_potential
    logical           :: find_forces, write_forces
    integer           :: n_save_freq, n_run
    character(len=40) :: output_file
    real(double)      :: total_energy

    ! Local variables
    logical           :: reset_L
    real(double)      :: vdW_energy_correction, vdW_xc_energy
    type(cq_timer)    :: tmr_l_energy, tmr_l_force, tmr_vdW
    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level
    
!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_E_and_F',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    call start_timer(tmr_std_eminimisation)
    ! reset_L = .true.  ! changed by TM, Aug 2008
    !   if (leqi(runtype,'static')) then
    !    reset_L = .false.
    !   else
    !    reset_L = .true.   ! temporary for atom movements
    !   end if

    ! Terrible coding.. Should be modified later. [2013/08/20 michi]
    reset_L = .false.
    if (.NOT. leqi(runtype,'static')) then
      !old if (.NOT. flag_MDold .AND. flag_LmatrixReuse) then
      if (flag_LmatrixReuse) then
        reset_L = .false.
        if (McWFreq.NE.0) then
          if (present(iter)) then
            if (mod(iter,McWFreq).EQ.0) then
              if (inode.EQ.ionode) write (io_lun,*) &
                   "Go back to McWeeny! Iteration:", iter
              reset_L = .true.
            endif
          endif
        endif
      ! Using an old-fashioned updates
      else
        reset_L = .true.
      endif
    endif

    ! Start timing the energy calculation
    call start_timer(tmr_l_energy, WITH_LEVEL)
    ! Now choose what we vary
    if (flag_LFD_minimise) then ! Vary everything, this flag is only for PAO-based multi-site SFs
       ! minimise by repeating LFD with updated SCF density
       call LFD_minimise(fixed_potential, vary_mu, n_L_iterations, L_tolerance, &
                         sc_tolerance, expected_reduction, total_energy, density)
       ! Numerical optimisation subsequently 
       if (flag_vary_basis) then
          if (UsePulay) then
             call pulay_min_pao(n_support_iterations, fixed_potential,&
                                vary_mu, n_L_iterations, L_tolerance, &
                                sc_tolerance, energy_tolerance,       &
                                total_energy, expected_reduction)
          else
             call vary_pao(n_support_iterations, fixed_potential, &
                           vary_mu, n_L_iterations, L_tolerance,  &
                           sc_tolerance, energy_tolerance,        &
                           total_energy, expected_reduction)
          end if
       endif
    else if (flag_vary_basis) then ! Vary everything: DM, charge density, basis set
       if (flag_basis_set == blips) then
          call vary_support(n_support_iterations, fixed_potential, &
                            vary_mu, n_L_iterations, L_tolerance,  &
                            sc_tolerance, energy_tolerance,        &
                            total_energy, expected_reduction)
       else if (flag_basis_set == PAOs) then
          if (flag_multisite .and. .not.flag_LFD_minimise) then
             if (inode==ionode) write(io_lun,'(/3x,A/)') &
                'WARNING: Numerical PAO minimisation will be performed without doing LFD_minimisation!'   
             !2017.Dec.28 TM: We need Selfconsistent Hamiltonian if this routine is called from control
             call new_SC_potl(.false., sc_tolerance, reset_L,             &
                              fixed_potential, vary_mu, n_L_iterations, &
                              L_tolerance, total_energy, backtrace_level)
          endif

          if (UsePulay) then
             call pulay_min_pao(n_support_iterations, fixed_potential,&
                                vary_mu, n_L_iterations, L_tolerance, &
                                sc_tolerance, energy_tolerance,       &
                                total_energy, expected_reduction)
          else
             call vary_pao(n_support_iterations, fixed_potential, &
                           vary_mu, n_L_iterations, L_tolerance,  &
                           sc_tolerance, energy_tolerance,        &
                           total_energy, expected_reduction)
          end if
       else
          call cq_abort("get_E_and_F: basis set undefined: ", &
                        flag_basis_set)
       end if
    else if (flag_self_consistent) then ! Vary only DM and charge density
       call new_SC_potl(.false., sc_tolerance, reset_L,           &
                        fixed_potential, vary_mu, n_L_iterations, &
                        L_tolerance, total_energy, backtrace_level)
    else ! Ab initio TB: vary only DM
       call FindMinDM(n_L_iterations, vary_mu, L_tolerance, inode, &
                      ionode, reset_L, .false., backtrace_level)
       call get_energy(total_energy, level=backtrace_level)
    end if
    ! Once ground state is reached, if we are doing deltaSCF, perform excitation
    ! and solve for the new ground state (on excited Born-Oppenheimer surface)
    if(flag_DeltaSCF.and.(.not.flag_excite)) then
       flag_excite = .true.
       if(inode==ionode.AND.iprint>2) write(io_lun,fmt='(2x,"Starting excitation loop")')
       if (flag_vary_basis) then ! Vary everything: DM, charge density, basis set
          if (flag_basis_set == blips) then
             call vary_support(n_support_iterations, fixed_potential, &
                  vary_mu, n_L_iterations, L_tolerance,  &
                  sc_tolerance, energy_tolerance,        &
                  total_energy, expected_reduction)
          else if (flag_basis_set == PAOs) then
             if (UsePulay) then
                call pulay_min_pao(n_support_iterations, fixed_potential,&
                     vary_mu, n_L_iterations, L_tolerance, &
                     sc_tolerance, energy_tolerance,       &
                     total_energy, expected_reduction)
             else
                call vary_pao(n_support_iterations, fixed_potential, &
                     vary_mu, n_L_iterations, L_tolerance,  &
                     sc_tolerance, energy_tolerance,        &
                     total_energy, expected_reduction)
             end if
          else
             call cq_abort("get_E_and_F: basis set undefined: ", &
                  flag_basis_set)
          end if
       else if (flag_self_consistent) then ! Vary only DM and charge density
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, total_energy, backtrace_level)

       else ! Ab initio TB: vary only DM
          call FindMinDM(n_L_iterations, vary_mu, L_tolerance, inode, &
               ionode, reset_L, .false., backtrace_level)
       end if
    end if
    ! calculate vdW energy correction to xc energy
    if (flag_vdWDFT) then
       call get_energy(total_energy)
       if (inode == ionode) &
            write (io_lun, '(/,10x,a,/)') &
                   'Calculating van der Waals correction...'
       call start_timer(tmr_vdW, WITH_LEVEL)
       call vdWXC_energy(density, vdW_xc_energy)
       vdW_energy_correction = vdW_xc_energy - xc_energy
       call stop_print_timer(tmr_vdW, "calculating vdW energy correction", &
                             IPRINT_TIME_THRES1)
       if (inode == ionode) then
          write (io_lun, '(10x,a,f25.15," ",a2)') &
                'van der Waals correction to XC-energy : ', &
                vdW_energy_correction * en_conv, en_units(energy_units)
          write (io_lun, '(10x,a,f25.15," ",a2)') &
                'Harris-Foulkes Energy after vdW correction : ', &
                (total_energy + vdW_energy_correction) * en_conv, &
                en_units(energy_units)
       end if
! LT_debug 2012/04/30 begin
!        flag_vdWDFT_slow = .false.
!        if (flag_vdWDFT_slow) then
!           if (inode == ionode) &
!                write (io_lun, '(/,8x,a,/)') &
!                      'Calculating van der Waals correction the slow way...'
!           call vdWXC_energy_slow(density, vdW_xc_energy)
!           vdW_energy_correction = vdW_xc_energy - xc_energy
!           if (inode == ionode) then
!              write (io_lun, '(10x,a,f25.15," ",a2)') &
!                    'van der Waals correction to XC-energy : ', &
!                    vdW_energy_correction * en_conv, en_units(energy_units)
!              write (io_lun, '(10x,a,f25.15," ",a2)') &
!                    'Harris-Foulkes Energy after vdW correction : ', &
!                    (total_energy + vdW_energy_correction) * en_conv, &
!                    en_units(energy_units)
!           end if
!        end if
! LT_debug 2012/04/30 end
    end if

!****lat<$
    call final_energy(backtrace_level)
!****lat>$

    ! output WFs or DOS
    if (flag_self_consistent.AND.(flag_out_wf.OR.flag_write_DOS)) then
       wf_self_con=.true.
       call FindMinDM(n_L_iterations, vary_mu, L_tolerance, inode,&
            ionode, reset_L, .false.)
       wf_self_con=.false.
    end if

    call stop_print_timer(tmr_l_energy, "calculating ENERGY", &
                          IPRINT_TIME_THRES1)
    if (atomch_output) call get_atomic_charge()

    if (find_forces) then
       ! Start timing the force calculation
      call start_timer(tmr_l_force, WITH_LEVEL)

      call force(fixed_potential, vary_mu, n_L_iterations, &
                 L_tolerance, sc_tolerance, total_energy,  &
                 expected_reduction, write_forces, backtrace_level)
 
      ! Stop timing the force calculation
      call stop_print_timer(tmr_l_force, "calculating FORCE", &
                            IPRINT_TIME_THRES1)
    end if

    !% NOTE: This call should be outside this subroutine [2013/08/20 michi]
    ! Writes out L-matrix at the PREVIOUS step
    !   --> Removed and moved to FindMinDM [2013/12/03 michi]

    !  Print results of local timers
    call stop_timer(tmr_std_eminimisation)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_E_and_F',echo=.true.)
!****lat>$

    return
  end subroutine get_E_and_F
  !!***

end module minimise
