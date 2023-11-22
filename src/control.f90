! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module control
! ------------------------------------------------------------------------------
! Code area 9: General
! ------------------------------------------------------------------------------

!!****h* Conquest/control *
!!  NAME
!!   control
!!  PURPOSE
!!   controls the run
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16:40, 2003/02/03 dave
!!  MODIFICATION HISTORY
!!   17:31, 2003/02/05 dave
!!    Changed flow of control so that the switch occurs on the local variable runtype
!!   14:26, 26/02/2003 drb 
!!    Small bug fixes
!!   10:06, 12/03/2003 drb 
!!    Added simple MD
!!   13:19, 22/09/2003 drb 
!!    Bug fixes in cg
!!   10:52, 13/02/2006 drb 
!!    General tidying related to new matrices
!!   2008/02/04 17:10 dave
!!    Changes for output to file not stdout
!!   2015/06/08 lat
!!    - Added experimental backtrace
!!   2017/10/24 zamaan
!!    Added md_ensemble variable for md_run
!!   2018/05/30 zamaan
!!    Implemented Berendsen thermostat and BerendsenEquil in NVT ensemble
!!   2018/07/12 zamaan
!!    Added SSM integrator, adapted from W. Shinoda et al., PRB 69:134103
!!    (2004)
!!   2019/05/09 zamaan
!!    Renamed init_ensemble to init_md, created new end_md subroutine 
!!    to deallocate any stray allocated arrays. Added heat flux calculation
!!    to md_run.
!!   2019/02/28 zamaan
!!    2 new subroutines for cell optimisation (double loop and single vector)
!!   2020/01/02 12:14 dave
!!    Added new module variable for L-BFGS history length
!!   2020/01/06 11:44 dave
!!    Removing Berendsen thermostat
!!   2022/09/29 16:47 dave
!!    Moved various subroutines to md_misc_module
!!   2022/12/12 11:41 dave
!!    Added SQNM maximum step size (sqnm_trust_step) as user-adjustable parameter
!!   2023/09/13 lu
!!    Added XSF and XSF output frequency as user-adjustable parameter
!!  SOURCE
!!
module control

  use datatypes
  use global_module, only: io_lun
  use GenComms,      only: cq_abort
  use timer_module,  only: start_timer, stop_timer, cq_timer 
  use timer_module,  only: start_backtrace, stop_backtrace
  use md_control,    only: md_ensemble, MDtimestep

  implicit none

  integer      :: MDn_steps 
  integer      :: MDfreq
  integer      :: XSFfreq
  integer      :: XYZfreq
  real(double) :: MDcgtol, sqnm_trust_step
  logical      :: CGreset
  integer      :: LBFGS_history

  ! Area identification
  integer, parameter, private :: area = 9

!!***

contains

!!****f* control/control_run *
!!
!!  NAME 
!!   control
!!  USAGE
!! 
!!  PURPOSE
!!   Very simple routine to control execution of Conquest
!!  INPUTS
!! 
!! 
!!  USES
!!   atoms, common, datatypes, dimens, ewald_module, force_module, GenComms, 
!!   matrix_data, pseudopotential_data
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   November 1998
!!  MODIFICATION HISTORY
!!   24/05/2001 dave
!!   - Indented, ROBODoc, stripped get_E_and_F calls
!!   11/06/2001 dave
!!   - Added RCS Id and Log tags and GenComms
!!   13/06/2001 dave
!!   - Adapted to use force_module
!!   17/06/2002 dave
!!   - Improved headers slightly
!!   16:42, 2003/02/03 dave
!!   - Changed to be control_run as part of control module
!!   17:31, 2003/02/05 dave
!!   - Removed old, silly MD and replaced with call to either static
!!    (i.e. just ground state) or CG
!!   10:06, 12/03/2003 drb 
!!   - Added MD
!!   2004/11/10, drb
!!   - Removed common use
!!   2008/02/13 08:23 dave
!!   - Added Pulay optimiser for atomic positions
!!   2011/12/11 L.Tong
!!   - Removed redundant parameter number_of_bands
!!   2012/03/27 L.Tong
!!   - Removed redundant input parameter real(double) mu
!!   2012/03/27 L.Tong
!!   - Removed redundant input parameter real(double) mu
!!   2014/10/05 L.Truflandier
!!   - Removed return for if(runtype,'static')
!!   2015/11/24 10:25 dave
!!    Removed redundant ewald use
!!   2017/08/29 jack baker & dave
!!    Added cell optimisation controls
!!   2019/06/10 zamaan
!!    Added two new methods for relaxation of both cell and ionic
!!    coordinates.
!!   2022/09/16 16:51 dave
!!    Added SQNM options for ion/cell optimisation (only partly complete)
!!  SOURCE
!!
  subroutine control_run(fixed_potential, vary_mu, total_energy)

    use datatypes
    use dimens,               only: r_core_squared, r_h
    use GenComms,             only: my_barrier, cq_abort, cq_warn
    use pseudopotential_data, only: set_pseudopotential
    use force_module,         only: tot_force
    use minimise,             only: get_E_and_F
    use global_module,        only: runtype, flag_self_consistent, &
                                    flag_out_wf, flag_write_DOS, wf_self_con, &
                                    flag_opt_cell, optcell_method, min_layer, flag_DM_converged
    use input_module,         only: leqi
    use store_matrix,         only: dump_pos_and_matrices

    implicit none

    ! Shared variables
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_energy

    ! Local variables
    integer, parameter:: backtrace_level = 0
    type(cq_timer)    :: backtrace_timer
    logical           :: NoMD, flag_ff, flag_wf
    integer           :: i, j

    real(double) :: spr, e_r_OLD

!****lat<$
    call start_backtrace(t=backtrace_timer,who='control_run',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$
    flag_DM_converged = .false.
    if ( leqi(runtype,'static') ) then
       !if(.NOT.flag_self_consistent.AND.(flag_out_wf.OR.flag_write_DOS)) return
       flag_ff = .true.
       flag_wf = .true.
       if (flag_out_wf.OR.flag_write_DOS) then
          ! This is done within get_E_and_F
          !wf_self_con=.true.
          flag_ff = .false.
          flag_wf = .false.
       endif
       call get_E_and_F(fixed_potential, vary_mu, total_energy,&
                        flag_ff, flag_wf, level=backtrace_level)
       !
    else if ( leqi(runtype, 'cg')    ) then
       if (flag_opt_cell) then
          select case(optcell_method)
          case(1)
             call cell_cg_run(fixed_potential, vary_mu, total_energy)
           case(2)
             call full_double_loop(fixed_potential, vary_mu, &
                                          total_energy)
           case(3)
             call full_cg_run_single_vector(fixed_potential, vary_mu, &
                  total_energy)
          end select
       else
          call cg_run(fixed_potential, vary_mu, total_energy)
       end if
       !
    else if ( leqi(runtype, 'md')    ) then
       call md_run(fixed_potential,     vary_mu, total_energy)
       !
    else if ( leqi(runtype, 'pulay') ) then
       call pulay_relax(fixed_potential,vary_mu, total_energy)
       !
    else if ( leqi(runtype, 'lbfgs') ) then
       call lbfgs(fixed_potential,vary_mu, total_energy)
       !
    else if ( leqi(runtype, 'sqnm') ) then
       if (flag_opt_cell) then
          select case(optcell_method)
          case(1)
             call cq_warn("control_run","SQNM cell optimisation not implemented; using CG.")
             call cell_cg_run(fixed_potential, vary_mu, total_energy)
             !call cell_sqnm(fixed_potential, vary_mu, total_energy)
          case(2)
             call cq_warn("control_run","SQNM cell optimisation not implemented; using CG.")
             call full_double_loop(fixed_potential, vary_mu, &
                  total_energy)
          case(3)
             call cq_abort("SQNM cell+ion optimisation not available yet.")
             !call full_sqnm(fixed_potential, vary_mu, total_energy)
          end select
       else
          call sqnm(fixed_potential, vary_mu, total_energy)
       end if
       !
    else if ( leqi(runtype, 'dummy') ) then
       call dummy_run(fixed_potential,  vary_mu, total_energy)
       !
    else
       call cq_abort('control: Runtype not specified !')
       !
    end if

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='control_run',echo=.true.)
!****lat>$

    ! Added if, otherwise step numbering is broken on MD restart - zamaan
    if (.not. leqi(runtype, 'md')) call dump_pos_and_matrices
    return
  end subroutine control_run
  !!***


  !!****f* control/cg_run *
  !!
  !!  NAME 
  !!   cg_run - Does CG minimisation
  !!  USAGE
  !!   cg_run(velocity,atmforce)
  !!  PURPOSE
  !!   Performs CG minimisation by repeated calling of minimiser
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:43, 2003/02/03 dave
  !!  MODIFICATION HISTORY
  !!   08:06, 2003/02/04 dave
  !!    Added vast numbers of passed arguments for get_E_and_F
  !!   17:26, 2003/02/05 dave
  !!    Corrected atom position arguments
  !!   14:27, 26/02/2003 drb 
  !!    Other small bug fixes
  !!   13:20, 22/09/2003 drb 
  !!    Added id_glob referencing for forces
  !!   2006/09/04 07:53 dave
  !!    Dynamic allocation for positions
  !!   2011/12/11 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Removed redundant input parameter real(double) mu
  !!   2013/08/21 M.Arita
  !!   - Added call for safemin2 necessary in reusing L-matrix
  !!    2017/08/29 jack baker & dave
  !!     Removed rcellx references (redundant)
  !!   2017/11/10 14:06 dave
  !!    Removing dump_InfoGlobal calls
  !!   2019/12/04 11:47 dave
  !!    Tweak to write convergence only on output process
  !!   2019/12/20 15:59 dave
  !!    Added choice of line minimiser: standard safemin2 or backtracking
  !!   2021/10/15 17:36 dave
  !!    total_energy now returns final energy
  !!   2022/08/03 15:23 dave
  !!    Add option for backtrack or adaptive backtracking
  !!   2022/09/16 17:10 dave
  !!    Added test for forces below threshold at the start of the run
  !!   2022/09/22 16:49 dave
  !!    Output tidying and implementation of FR for backtrack (PR gives gamma < 0)
  !!  SOURCE
  !!
  subroutine cg_run(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, only: iprint_gen, ni_in_cell, x_atom_cell,  &
                             y_atom_cell, z_atom_cell, id_glob,    &
                             atom_coord, &
                             area_general, iprint_MD,              &
                             IPRINT_TIME_THRES1, min_layer
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: adapt_backtrack_linemin, backtrack_linemin, &
                             safemin2, cg_line_min, safe, adapt_backtrack, backtrack
    use GenComms,      only: gsum, myid, inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: tot_force
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop, write_xsf, write_extxyz, &
                             print_atomic_positions, return_prefix
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  ONLY: dump_InfoMatGlobal, dump_pos_and_matrices
    use md_control,    only: flag_write_xsf, flag_write_extxyz

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy
    
    ! Local variables
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma,gg1, test_dot
    integer        :: i,j,k,iter,length, jj, lun, stat
    logical        :: done
    type(cq_timer) :: tmr_l_iter
    real(double), allocatable, dimension(:,:) :: cg, old_force
    real(double), allocatable, dimension(:)   :: x_new_pos, y_new_pos,&
         z_new_pos

    character(len=12) :: subname = "cg_run: "
    character(len=120) :: prefix, prefixGO
    character(len=10) :: prefixF = "          "

    prefix = return_prefix(subname, min_layer)
    prefixGO = return_prefix("GeomOpt", min_layer)
    allocate(cg(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ",&
                       ni_in_cell, stat)
    allocate(old_force(3,ni_in_cell), STAT=stat)
    allocate(x_new_pos(ni_in_cell), y_new_pos(ni_in_cell), &
             z_new_pos(ni_in_cell), STAT=stat)
    if(stat/=0) &
         call cq_abort("Error allocating _new_pos in control: ", &
                       ni_in_cell,stat)
    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    if (inode == ionode) then
       if(cg_line_min==safe) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with safemin line minimisation"
       else if (cg_line_min==backtrack) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with backtracking line minimisation"
       else if (cg_line_min==adapt_backtrack) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with adaptive backtracking line minimisation"
       end if
       if(iprint_gen + min_layer>0) write(io_lun, fmt='(4x,a,f8.4,x,a2,"/",a2,a,i4)') &
            trim(prefix)//" Tolerance: ",for_conv*MDcgtol,en_units(energy_units), d_units(dist_units),&
            " Maximum steps: ",MDn_steps
    end if
    cg = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3 * ni_in_cell
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    min_layer = min_layer - 1
    if (iprint_MD + min_layer > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.,0)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    end if
    min_layer = min_layer + 1
    call dump_pos_and_matrices
    call get_maxf(max)
    if (inode==ionode) then
      write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," E: ",f16.8,x,a2," dE: ",f14.8,x,a2/)') & 
           trim(prefixGO)//" - Iter: ",0, for_conv*max, en_units(energy_units), d_units(dist_units), &
           en_conv*energy0, en_units(energy_units), en_conv*dE, en_units(energy_units)
    end if
    ! Check for trivial case where forces are converged
    if (abs(max) < MDcgtol) then
       done = .true.
       if (myid == 0) then
          write(io_lun,'(4x,a,i4,a)') trim(prefixGO)//" converged in ", iter, " iterations"
          write(io_lun, fmt='(4x,a,f12.5,x,a2,"/",a2)') &
               trim(prefix)//" maximum force below threshold: ", &
               for_conv*max, en_units(energy_units), d_units(dist_units)
       end if
       return
    end if

    iter = 1
    ggold = zero
    old_force = zero
    energy1 = energy0
    do while (.not. done)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       ! Construct ratio for conjugacy
       gg = zero
       gg1 = zero! PR
       do j = 1, ni_in_cell
          gg = gg +                              &
               tot_force(1,j) * tot_force(1,j) + &
               tot_force(2,j) * tot_force(2,j) + &
               tot_force(3,j) * tot_force(3,j)
          gg1 = gg1 +                              &
               tot_force(1,j) * old_force(1,j) + &
               tot_force(2,j) * old_force(2,j) + &
               tot_force(3,j) * old_force(3,j)
       end do
       if (abs(ggold) < 1.0e-6_double) then
          gamma = zero
       else
          gamma = (gg-gg1)/ggold ! PR - change to gg/ggold for FR
       end if
       if(gamma<zero) then
          if(inode==ionode .and. iprint_MD + min_layer > 2) &
               write(io_lun,fmt='(4x,a)') trim(prefix)//" gamma < 0.  Setting to gg/ggold."
          gamma = gg/ggold ! Default to FR if PR gives negative
       end if
       if (CGreset) then
          if (gamma > one) then
             if (inode == ionode .and. iprint_MD + min_layer > 1) &
                  write(io_lun,fmt='(4x,a)') &
                  trim(prefix)//" gamma>1, so resetting direction"
             gamma = zero
          end if
       end if
       if (inode == ionode .and. iprint_MD + min_layer > 1) then
          write(io_lun,fmt='(4x,a,f12.8)') trim(prefix)//' gamma = ', gamma
       end if
       ggold = gg
       ! Build search direction
       test_dot = zero
       do j = 1, ni_in_cell
          jj = id_glob(j)
          cg(1,j) = gamma*cg(1,j) + tot_force(1,jj)
          cg(2,j) = gamma*cg(2,j) + tot_force(2,jj)
          cg(3,j) = gamma*cg(3,j) + tot_force(3,jj)
          test_dot = test_dot + cg(1,j)*tot_force(1,jj) + cg(2,j)*tot_force(2,jj) &
               + cg(3,j)*tot_force(3,jj)
          x_new_pos(j) = x_atom_cell(j)
          y_new_pos(j) = y_atom_cell(j)
          z_new_pos(j) = z_atom_cell(j)
       end do
       ! Test
       if(test_dot<zero) then
          if(inode==ionode.AND.iprint_MD + min_layer > 0) &
               write(io_lun,fmt='(4x,a,e14.6)') &
               trim(prefix)//" search dot gradient < zero. Reset direction ",test_dot
          gamma = zero
          do j = 1, ni_in_cell
             jj = id_glob(j)
             cg(1,j) = tot_force(1,jj)
             cg(2,j) = tot_force(2,jj)
             cg(3,j) = tot_force(3,jj)
          end do
       end if
       old_force = tot_force
       ! Minimise in this direction
       min_layer = min_layer - 1
       if(cg_line_min==safe) then
          call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0,&
               energy1, fixed_potential, vary_mu)
       else if(cg_line_min==backtrack) then
          call backtrack_linemin(cg, energy0, energy1, fixed_potential, vary_mu)
       else if(cg_line_min==adapt_backtrack) then
          call adapt_backtrack_linemin(cg, energy0, energy1, fixed_potential, vary_mu)
       end if
       min_layer = min_layer + 1
       ! Output positions
       if (myid == 0 .and. iprint_gen + min_layer > 1) then
          call print_atomic_positions
       end if
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)
       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       call get_maxf(max)
       ! Output and energy changes
       dE = energy1 - energy0

       if (inode==ionode) then
          write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," E: ",f16.8,x,a2," dE: ",f14.8,x,a2/)') & 
               trim(prefixGO)//" - Iter: ", iter, for_conv*max, &
               en_units(energy_units), d_units(dist_units), &
               en_conv*energy1, en_units(energy_units), en_conv*dE, en_units(energy_units)
          if (iprint_MD + min_layer > 1) then
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force Residual:     ", &
                  for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Maximum force:      ", &
                  for_conv*max, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force tolerance:    ", &
                  for_conv*MDcgtol, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2)') prefixF(1:-2*min_layer)//"Energy change:      ", &
                  en_conv*dE, en_units(energy_units)
          end if
       end if

       energy0 = energy1
       total_energy = energy0
       !energy1 = abs(dE)
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write(io_lun, fmt='(4x,a,i6)') &
                     trim(prefix)//" Exceeded number of MD steps: ",iter
       end if
       if (abs(max) < MDcgtol) then
          done = .true.
          if (myid == 0) then
             write(io_lun,'(/4x,a,i4,a)') trim(prefixGO)//" converged in ", iter, " iterations"
             write(io_lun, fmt='(4x,a,f12.5,x,a2,"/",a2)') &
                  trim(prefix)//" maximum force below threshold: ",&
                  for_conv*max, en_units(energy_units), d_units(dist_units)
          end if
       end if
       iter = iter + 1

       !call dump_pos_and_matrices
       
       call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
       if (.not. done) call check_stop(done, iter)
    end do
    total_energy = energy0
    ! Output final positions
!    if(myid==0) call write_positions(parts)
    deallocate(z_new_pos, y_new_pos, x_new_pos, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating _new_pos in control: ", &
                       ni_in_cell,stat)
    deallocate(old_force, STAT=stat)
    deallocate(cg, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating cg in control: ", ni_in_cell,stat)
    call reg_dealloc_mem(area_general, 6*ni_in_cell, type_dbl)

  end subroutine cg_run
  !!***

!!****f* control/md_run *
!!
!!  NAME 
!!   md_run
!!  USAGE
!! 
!!  PURPOSE
!!   Does a QUENCHED MD run
!!   Now, this can also employ NVE-MD
!!
!!   Note that we require fire_N_below_thresh consecutive steps with the force below threshold
!!   before convergence is accepted
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   10:07, 12/03/2003 drb 
!!  MODIFICATION HISTORY
!!   2011/12/07 L.Tong
!!    - changed 0.0_double to zero in numbers module
!!    - updated calls to get_E_and_F
!!   2011/12/11 L.Tong
!!    Removed redundant parameter number_of_bands
!!   2012/03/27 L.Tong
!!   - Removed redundant input parameter real(double) mu
!!   2013/08/20 M.Arita
!!   - Added 'iter' in velocityVerlet
!!   - Added calls for reading/dumping global information
!!   - Modified initial and final MD steps
!!   2013/12/02 M.Arita
!!   - Added calls for L-matrix reconstruction & update_H
!!   2013/12/03 M.Arita
!!   - Added calls for Ready_XLBOMD and Do_XLBOMD
!!   2014/02/03 M.Arita
!!   - Adopted new integration scheme instead of velocityVerlet
!!     with routines in Integrators_module
!!   - Added calls for constraint-MD
!!   2015/06/19 13:52 dave
!!   - Included FIRE routines but moved reading of parameters to io_module
!!   2016/07/21 10:17 dave
!!   - Changed format for MD Step output line to add more space
!!   2017/02/23 dave
!!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
!!   2017/05/09 dave
!!    Added calls to read L-matrix for both spin channels
!!   2017/10/20 09:32 dave
!!    Removed fire_step_max from call to fire_qMD (now module variable set by user)
!!   2017/10/25 zamaan
!!    Added refactored NVT
!!   2017/11/02 zamaan
!!    Changed velocity local variable to ion_velocity from md_control
!!   2017/11/06 zamaan
!!    added lattice_vec for md_model to track cell vectors
!!   2017/11/10 dave
!!    Removed calls to dump K matrix (now done in DMMinModule)
!!   2017/11/10 15:15 dave
!!    Created variable fire_N_below_thresh to set number of consecutive steps below threshold
!!   2018/01/22 tsuyoshi (with dave)
!!    Changes to use new atom update routines
!!   2018/04/20 nakata
!!    Changed to initialise n_stop_qMD
!!   2018/04/24 zamaan
!!    Initial NPT implementation (isotropic fluctuations only)
!!   2018/05/29 zamaan
!!    Added time step argument to init_thermo calls
!!   2018/05/30 zamaan
!!    Switch on NHC after using BerendsenEquil in NVT ensemble
!!   2018/7/23 zamaan
!!    Replaced calculate_kinetic_energy calls with get_temperature_and_ke,
!!    which also computes the kinetic stress
!!   2018/8/11 zamaan
!!    Moved ionic position and box update to its own subroutine; moved 
!!    velocity initialisation to init_md
!!   2019/05/09 zamaan
!!    Moved velocity array allocation/deallocation to init_md/end_md
!!   2020/01/06 15:40 dave
!!    Add pressure-based termination for equilibration and remove Berendsen thermostat
!!   2023/07/22 J.Lin
!!    Added machine learning statements
!!   2023/09/13 lu
!!    Add parameters for xsf and xyz output frequency
!!  SOURCE
!!
  subroutine md_run (fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module,  only: iprint_gen, ni_in_cell, x_atom_cell,    &
                              y_atom_cell, z_atom_cell, area_general, &
                              flag_read_velocity, flag_quench_MD,     &
                              temp_ion, flag_MDcontinue, MDinit_step, &
                              flag_LmatrixReuse,flag_XLBOMD,          &
                              flag_dissipation,flag_FixCOM,           &
                              flag_fire_qMD, flag_diagonalisation,    &
                              nspin, flag_thermoDebug, flag_baroDebug,&
                              flag_move_atom,rcellx, rcelly, rcellz,  &
                              flag_Multisite,flag_SFcoeffReuse, &
                              atom_coord, flag_quench_MD, atomic_stress, &
                              non_atomic_stress, flag_heat_flux, flag_MLFF, min_layer
    use group_module,   only: parts
    use minimise,       only: get_E_and_F
    use move_atoms,     only: velocityVerlet, updateIndices,           &
                              init_velocity,update_H,check_move_atoms, &
                              update_atom_coord,wrap_xyz_atom_cell,    &
                              zero_COM_velocity,                       &
                              fac_Kelvin2Hartree
    use GenComms,       only: gsum, my_barrier, inode, ionode, gcopy
    use GenBlas,        only: dot
    use force_module,   only: tot_force, stress
    use io_module,      only: read_fire, check_stop, print_atomic_positions, return_prefix
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use move_atoms,     only: fac_Kelvin2Hartree, update_pos_and_matrices, &
                              updateL, updateLorK, updateSFcoeff
    use store_matrix,   ONLY: dump_InfoMatGlobal,grab_InfoMatGlobal, &
                    matrix_store_global, InfoMatrixFile, dump_pos_and_matrices
    use mult_module,    ONLY: matL,L_trans,matK,S_trans, matrix_scale
    use matrix_data,    ONLY: Lrange,Hrange, rcut, max_range
    use XLBOMD_module,  ONLY: Ready_XLBOMD, Do_XLBOMD
    use Integrators,    ONLY: vVerlet_v_dthalf,vVerlet_r_dt, fire_qMD, &
                              fire_N_below_thresh
    use constraint_module, ONLY: correct_atomic_position,correct_atomic_velocity, &
         ready_constraint,flag_RigidBonds
    use input_module,   ONLY: leqi
    use cover_module,   only: BCS_parts, make_cs, deallocate_cs, make_iprim, &
                              send_ncover
    use md_model,       only: type_md_model, heat_flux
    use md_control,     only: type_thermostat, type_barostat, md_n_nhc, &
                              md_n_ys, md_n_mts, ion_velocity, lattice_vec, &
                              md_baro_type, target_pressure, md_ndof_ions, &
                              md_equil_steps, md_equil_press, md_tau_T, md_tau_P, &
                              md_thermo_type
    use md_misc,        only: write_md_data, get_heat_flux, &
                              update_pos_and_box, integrate_pt, init_md, end_md
    use atoms,          only: distribute_atoms,deallocate_distribute_atom
    use global_module,  only: atom_coord_diff, iprint_MD, area_moveatoms
    use mlff,           only: get_MLFF
    use mlff_type,      only: flag_debug_mlff,flag_time_mlff
    use mpi

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical      :: vary_mu, fixed_potential
    real(double) :: total_energy
    
    ! Local variables
    integer       ::  iter, i, j, k, length, stat, i_first, i_last, md_ndof
    integer       :: nequil ! number of equilibration steps - zamaan
    real(double)  :: energy1, energy0, dE, max, g0
    real(double)  :: t0,t1,t2
    logical       :: done,second_call
    logical,allocatable,dimension(:) :: flag_movable

    type(matrix_store_global) :: InfoGlob
    type(InfoMatrixFile),pointer:: InfoL(:)

    !! quenched MD optimisation is stopped
    !! if the maximum force component is bellow threshold
    !! for 3 consecutive iterations

    ! FIRE parameters
    integer :: step_qMD, n_stop_qMD, fire_N, fire_N2
    real(double) :: fire_step_max, fire_P0, fire_alpha
    integer :: iter_MD, final_call

    ! Storage for MD data
    type(type_md_model)           :: mdl

    ! thermostat, barostat
    type(type_thermostat), target :: thermo
    type(type_barostat), target   :: baro

    character(len=12) :: subname = "md_run: "
    character(len=120) :: prefix
    character(len=10) :: prefixF = "          "

    prefix = return_prefix(subname, min_layer)
    n_stop_qMD = 0
    final_call = 1
    energy0 = zero
    energy1 = zero
    dE = zero
    length = 3*ni_in_cell

    ! initialisation/contintuation when we're using equilibration
    if (flag_MDcontinue) then
      if (MDinit_step >= md_equil_steps) then
        nequil = 0
      else
        nequil = md_equil_steps - MDinit_step
      end if
    else
      nequil = md_equil_steps 
    end if

    ! Initialise number of degrees of freedom
    md_ndof = 3*ni_in_cell
    md_ndof_ions = md_ndof
    do i=1,ni_in_cell
      do j=1,3
        if (flag_move_atom(j,i) .eqv. .false.) md_ndof = md_ndof-1
      end do
    end do
    !2020/Jul/30 TM
     if(flag_FixCOM) md_ndof = md_ndof-3

    if (inode==ionode) &
      write(io_lun,'(4x,a,i6," steps")') trim(prefix)//" starting MD run with ",MDn_steps

    ! Thermostat/barostat initialisation
    call init_md(baro, thermo, mdl, md_ndof, nequil)
    if (inode==ionode .and. flag_debug_mlff) then
      write(*,'(4x,a,i6," steps")') trim(prefix)//" starting MD run with ",MDn_steps
      write(*,*) 'check stress after init_md ini:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
    end if
    call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                      mdl%ion_kinetic_energy)
    if (inode==ionode .and. flag_debug_mlff) &
       write(*,*) 'check stress after thermo%get_temperature_and_ke ini:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
    call baro%get_pressure_and_stress
    !! have pressure from here
    if (inode==ionode .and. flag_debug_mlff) &
       write(*,*) 'check stress after baro%get_pressure_and_stress ini:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
    call mdl%get_cons_qty

    ! Find energy and forces
    min_layer = min_layer - 1
    if (flag_MLFF) then
      call get_MLFF
    else
      if (flag_fire_qMD) then
        call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.)
      else
        call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.)
      end if
    end if
    min_layer = min_layer + 1
    mdl%dft_total_energy = energy0

    ! XL-BOMD
    if (flag_XLBOMD .AND. flag_dissipation) &
      call Ready_XLBOMD()

    ! Get converted 1-D array for flag_atom_move
    allocate (flag_movable(3*ni_in_cell), STAT=stat)
    if (stat.NE.0) call cq_abort('Error allocating flag_movable: ',3*ni_in_cell)
    call check_move_atoms(flag_movable)

    ! constraint-MD
    if (flag_RigidBonds) call ready_constraint()

    ! Get an initial MD step
    if (.NOT. flag_MDcontinue) then
      i_first = 1
      i_last  = MDn_steps
    else
      i_first = MDinit_step + 1
      i_last = i_first + MDn_steps - 1
    endif
    if (flag_MLFF) then
      call dump_InfoMatGlobal(index=0,velocity=ion_velocity,MDstep=i_first)
    else
      call dump_pos_and_matrices(index=0,MDstep=i_first,velocity=ion_velocity)
    end if

    if (inode==ionode) then
       write(io_lun,'(/4x,"MD step: ",i6," KE: ",f18.8,x,a2," &
            &IntEnergy: ",f20.8,x,a2," TotalEnergy: ",f20.8,x,a2)') &
            i_first-1, en_conv*mdl%ion_kinetic_energy, en_units(energy_units), &
            en_conv*mdl%dft_total_energy, en_units(energy_units), &
            en_conv*(mdl%ion_kinetic_energy+mdl%dft_total_energy), en_units(energy_units)
       if(iprint_MD + min_layer > 0) write(io_lun, fmt='(4x,"Ionic temperature       : ",f15.8," K")') &
            mdl%ion_kinetic_energy / (three_halves * real(ni_in_cell,double) * fac_Kelvin2Hartree) 
    end if

    if (flag_fire_qMD) then
       step_qMD = i_first ! SA 20150201
       done = .false.     ! SA 20150201
       ! reading FIRE parameters of the last run
       call read_fire(fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha)
    endif

    if (flag_heat_flux) &
      call get_heat_flux(atomic_stress, ion_velocity, heat_flux)

    if (inode==ionode .and. flag_debug_mlff) &
       write(*,*) 'check stress before mdl%get_cons_qty:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
    if (.not. flag_MDcontinue) then
       ! Check this: it fixes H' for NVE but needs NVT confirmation
       ! DRB & TM 2020/01/24 12:03
       call mdl%get_cons_qty
       call write_md_data(i_first-1, thermo, baro, mdl, nequil, MDfreq, XSFfreq, XYZfreq)
    end if

    do iter = i_first, i_last ! Main MD loop
       t0=MPI_wtime()
       mdl%step = iter
       if (inode==ionode .and. iprint_MD + min_layer > 0) &
            write(io_lun,fmt='(/4x,a,i5)') trim(prefix)//" iteration ",iter

       if (inode==ionode .and. (flag_debug_mlff .or. flag_time_mlff)) &
            write(*,fmt='(/4x,a,i5)') trim(prefix)//" iteration ",iter

       if (flag_heat_flux) then
         atomic_stress = zero
         non_atomic_stress = zero
       end if
       call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                          mdl%ion_kinetic_energy)
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress after thermo%get_temperature_and_ke 1:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
       call baro%get_pressure_and_stress
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress after baro%get_pressure_and_stress 1:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa

       ! thermostat/barostat (MTTK splitting of Liouvillian)
       call integrate_pt(baro, thermo, mdl, ion_velocity)
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress after integrate_pt 1:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa

       !! For Debuggging !!
       !     call dump_pos_and_matrices(index=1,MDstep=i_first)
       !! For Debuggging !!

       !%%! Evolve atoms - either FIRE (quenched MD) or velocity Verlet
       if (flag_fire_qMD) then
          call fire_qMD(MDtimestep,ion_velocity,tot_force,flag_movable,iter, &
                        fire_N,fire_N2,fire_P0,fire_alpha) ! SA 20150204
       else
          call vVerlet_v_dthalf(MDtimestep,ion_velocity,tot_force,flag_movable)
          ! The velocity Verlet dt position update plus box update(s)
          call update_pos_and_box(baro, nequil, flag_movable)
       end if
       ! Constrain position
       if (flag_RigidBonds) &
         call correct_atomic_position(ion_velocity,MDtimestep)
       ! Update box if equilibration has finished
       if (md_ensemble(2:2) == 'p') then
          if (nequil < 1) then
             call baro%update_cell
          end if
       end if
       t1=MPI_wtime()
       if(flag_SFcoeffReuse) then
         call update_pos_and_matrices(updateSFcoeff,ion_velocity)
       else
         call update_pos_and_matrices(updateLorK,ion_velocity)
       endif

       t2=MPI_wtime()
       if (inode==ionode .and. flag_time_mlff) &
            write(*,2023) 'Time at update_pos_and_matrices in MD:', t2-t1
       2023 format(a,e16.6)
       if (flag_XLBOMD) call Do_XLBOMD(iter,MDtimestep)
       if (.not. flag_MLFF) call update_H(fixed_potential)

       !2018.Jan.4 TM 
       !   We need to update flag_movable, since the order of x_atom_cell (or id_glob) 
       !   may change after the atomic positions are updated.
       call check_move_atoms(flag_movable)
       min_layer = min_layer - 1
       if (flag_fire_qMD) then
          if (flag_MLFF) then
            call get_MLFF
            if (inode==ionode .and. flag_debug_mlff) &
               write(*,*) 'check stress after gret_MLFF:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
            call check_stop(done, iter)
            call dump_InfoMatGlobal(index=0,velocity=ion_velocity,MDstep=iter)
          else
            call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .true.,iter)
            call check_stop(done, iter)   !2019/Nov/14
            call dump_pos_and_matrices(index=0,MDstep=iter,velocity=ion_velocity)
          end if ! flag_MLFF
       else
          if (flag_MLFF) then
            t1=MPI_wtime()
            call get_MLFF
            t2=MPI_wtime()
            if (inode==ionode .and. flag_time_mlff) &
               write(*,2023) 'Time at get_E_and_F_ML in MD:', t2-t1
            if (inode==ionode .and. flag_debug_mlff) &
               write(*,*) 'check stress after get_MLFF:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
            call check_stop(done, iter)
            ! Here so that the kinetic stress is reported just after the
            ! static stress - zamaan
            if (inode == ionode .and. iprint_MD > 2) then
              write(io_lun,fmt='(/4x,a, 3f15.8,a3)') trim(prefix)//" Kinetic stress    ",&
                baro%ke_stress(1,1), baro%ke_stress(2,2), baro%ke_stress(3,3), &
                en_units(energy_units)
            end if
            call dump_InfoMatGlobal(index=0,velocity=ion_velocity,MDstep=iter)
          else
            call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .false.,iter)
            call check_stop(done, iter)   !2019/Nov/14
            ! Here so that the kinetic stress is reported just after the
            ! static stress - zamaan
            if (inode == ionode .and. iprint_MD > 2) then
              write(io_lun,fmt='(/4x,a, 3f15.8,a3)') trim(prefix)//" Kinetic stress    ",&
                baro%ke_stress(1,1), baro%ke_stress(2,2), baro%ke_stress(3,3), &
                en_units(energy_units)
            end if
            call dump_pos_and_matrices(index=0,MDstep=iter,velocity=ion_velocity)
          end if ! flag_MLFF
          call vVerlet_v_dthalf(MDtimestep,ion_velocity,tot_force,flag_movable,second_call)
       end if
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress after if flag_fire_qMD:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext*HaBohr3ToGPa
       t1=MPI_wtime()
       min_layer = min_layer + 1
       ! Update DFT energy
       mdl%dft_total_energy = energy1
       !******
       !2019/Nov/14 tsuyoshi
       ! we would dump the files(InfoGlobal and *matrix2.i**.p******)  once at the end of each job.
       ! We should call "dump_pos_and_matrices" after calling "check_stop", since 
       ! this subroutine would change the format (binary or ascii) of the files,
       ! if user sets different "IO.MatrixFile.BinaryFormat.Grab" and "IO.MatrixFile.BinaryFormat.Dump".
       !  (reading binary files at the beginning and writing ascii files at the end is possible.)
       ! 
       ! In the near future,  (depending on "IO.DumpEveryStep")
       !  we would use "store_pos_and_matrices" at each MD step, 
       !  while call "dump_pos_and_matrices" if "done" is true.
       !******

       thermo%ke_ions = mdl%ion_kinetic_energy
       call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                          mdl%ion_kinetic_energy)
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress after thermo%get_temperature_and_ke 2:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext/HaBohr3ToGPa
       call baro%get_pressure_and_stress
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress after baro%get_pressure_and_stress 2:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext/HaBohr3ToGPa

       ! thermostat/barostat (MTTK splitting of Liouvillian)
       call mdl%get_cons_qty
       call integrate_pt(baro, thermo, mdl, ion_velocity, second_call)
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress after integrate_pt 2:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext/HaBohr3ToGPa

       ! Constrain velocity
       if (flag_RigidBonds) call correct_atomic_velocity(ion_velocity)
       if (flag_FixCOM) call zero_COM_velocity(ion_velocity)
       !%%! END of Evolve atoms

       ! Print out energy
       if (inode==ionode) then
         write(io_lun,'(/4x,"MD step: ",i6," KE: ",f18.8,x,a2," &
                       &IntEnergy: ",f20.8,x,a2," TotalEnergy: ",f20.8,x,a2)') &
                       iter, en_conv*mdl%ion_kinetic_energy, en_units(energy_units), &
                       en_conv*mdl%dft_total_energy, en_units(energy_units), &
                       en_conv*(mdl%ion_kinetic_energy+mdl%dft_total_energy), en_units(energy_units)
         if(iprint_MD + min_layer > 0) write(io_lun, fmt='(4x,"Ionic temperature       : ",f15.8," K")') &
              mdl%ion_kinetic_energy / (three_halves * real(ni_in_cell,double) * fac_Kelvin2Hartree) 
       end if

       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
          end do
       end do

       ! Output and energy changes
       dE = energy1 - energy0
       energy0 = energy1
       energy1 = abs(dE)
       if(inode==ionode .and. iprint_MD>1) then
         write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') &
              prefixF(1:-2*min_layer)//"Maximum force:      ", &
              for_conv*max, en_units(energy_units), d_units(dist_units)
         write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') &
              prefixF(1:-2*min_layer)//"Force Residual:     ", &
              for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), d_units(dist_units)
         write(io_lun,'(4x,a,f19.8," ",a2)') &
              prefixF(1:-2*min_layer)//"Energy change:      ", &
              en_conv*dE, en_units(energy_units)
       end if

       ! Compute and print the conserved quantity and its components
       if (leqi(thermo%thermo_type, 'nhc')) call thermo%get_thermostat_energy
       if (leqi(baro%baro_type, 'pr')) call baro%get_barostat_energy(final_call)
       call mdl%print_md_energy()

       ! Output positions
       if (inode==ionode .and. iprint_gen > 1) then
          call print_atomic_positions
       endif
       call my_barrier

       ! The kinetic component of stress changes after the second velocity
       ! update 
       call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                          mdl%ion_kinetic_energy, &
                                          final_call)
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress thermo%get_temperature_and_ke second velocity:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext/HaBohr3ToGPa
       call baro%get_pressure_and_stress(final_call)
       if (inode==ionode .and. flag_debug_mlff) &
          write(*,*) 'check stress baro%get_pressure_and_stress second velocity:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext/HaBohr3ToGPa
 
       if (nequil > 0) then
          nequil = nequil - 1
          if (abs(baro%P_int - baro%P_ext) < md_equil_press) nequil = 0
          if (nequil == 0) then
             if (inode==ionode) &
                  write(io_lun, '(4x,a)') trim(prefix)//" equilibration finished, &
                  &starting extended system dynamics."
             call init_md(baro, thermo, mdl, md_ndof, nequil, second_call)
          end if
       end if

       ! Compute heat flux
       if (flag_heat_flux) &
         call get_heat_flux(atomic_stress, ion_velocity, heat_flux)

       t2=MPI_wtime()
       if (inode==ionode .and. flag_time_mlff) &
           write(*,2023) 'Time after get_E_and_F_ML in MD:', t2-t1
       t1=MPI_wtime()
       ! Write all MD data and checkpoints to disk
       call write_md_data(iter, thermo, baro, mdl, nequil, MDfreq, XSFfreq, XYZfreq)
       if (inode==ionode .and. flag_debug_mlff) &
           write(*,*) 'check stress write_md_data second velocity:', stress,baro%P_int*HaBohr3ToGPa,&
               baro%P_ext/HaBohr3ToGPa
       t2=MPI_wtime()
       if (inode==ionode .and. flag_time_mlff) &
            write(*,2023) 'Time at save_MD_data in MD:', t2-t1
       !call check_stop(done, iter) ! moved above. 2019/Nov/14 tsuyoshi
       if (flag_fire_qMD.OR.flag_quench_MD) then
          if (abs(max) < MDcgtol) then
             if ((iter - step_qMD) > 1) then
                n_stop_qMD = 0
             end if
             n_stop_qMD = n_stop_qMD + 1
             step_qMD = iter
             if (inode==ionode) then
                write(io_lun, fmt='(4x,a,f12.5,x,a2,"/",a2)') &
                     trim(prefix)//" maximum force below threshold: ", &
                     for_conv*max, en_units(energy_units), d_units(dist_units)
             end if
             if (n_stop_qMD > fire_N_below_thresh) then
                done = .true.
             end if
          end if
       end if

       if (done) exit
    end do ! Main MD loop

    call end_md(thermo, baro)

    return
  end subroutine md_run
  !!***

  !!****f* control/dummy_run *
  !! PURPOSE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!
  !! MODIFICATION HISTORY
  !!   2011/12/07 L.Tong
  !!     - Updated calls to get_E_and_F
  !!     - changed 0.0_double to zero in numbers module
  !!   2011/12/11 L.Tong
  !!     Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Removed redundant input parameter real(double) mu
  !!   2012/06/18 L.Tong
  !!   - Removed unused variable velocity
  !! SOURCE
  !!
  subroutine dummy_run(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use global_module,  only: iprint_gen, ni_in_cell, x_atom_cell, &
                              y_atom_cell, z_atom_cell
    use group_module,   only: parts
    use primary_module, only: bundle
    use minimise,       only: get_E_and_F
    use move_atoms,     only: velocityVerlet, updateIndices2
    use GenComms,       only: gsum, myid, my_barrier
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_positions

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical      :: vary_mu, fixed_potential
    real(double) :: total_energy
    
    ! Local variables
    integer      :: iter, i, k, length, stat
    real(double) :: temp, KE, energy1, energy0, dE, max, g0

    MDfreq = 100
    if (myid == 0 .and. iprint_gen > 0) write (io_lun, 2) MDn_steps
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., &
                     .true.)
    do iter = 1, MDn_steps
       if (myid == 0) &
            write (io_lun, fmt='(4x,"Dummy run, iteration ",i5)') iter
       ! Output positions
       if (myid == 0 .and. iprint_gen > 1) then
          do i = 1, ni_in_cell
             write (io_lun, 1) i, x_atom_cell(i), y_atom_cell(i), &
                               z_atom_cell(i)
          end do
       end if
       call updateIndices2(.true., fixed_potential)
       call get_E_and_F(fixed_potential, vary_mu, energy1, .true., &
                        .true.)
       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = tot_force(k,i)
          end do
       end do
       ! Output and energy changes
       dE = energy1 - energy0
       if (myid == 0) write (io_lun, 6) max
       if (myid == 0) write (io_lun, 4) dE
       if (myid == 0) write (io_lun, 5) sqrt(g0/ni_in_cell)
       energy0 = energy1
       energy1 = abs(dE)
       if (myid == 0 .and. mod(iter, MDfreq) == 0) &
            call write_positions(iter, parts)
       call my_barrier
    end do
    
    return
1   format(4x,'Atom ',i4,' Position ',3f15.8)
2   format(4x,'Welcome to dummy_run. Doing ',i4,' steps')
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Energy change           : ',f15.8)
5   format(4x,'Force Residual          : ',f15.8)
6   format(4x,'Maximum force component : ',f15.8)
  end subroutine dummy_run
  !!*****


  !!****f* control/pulay_relax *
  !!
  !!  NAME 
  !!   pulay_relax
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Performs pulay relaxation
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2001 in ParaDens
  !!  MODIFICATION HISTORY
  !!   2007/05/24 16:12 dave
  !!    Incorporated into Conquest
  !!   2011/12/07 L.Tong
  !!    - updated calls to get_E_and_F
  !!    - changed 0.0_double to zero
  !!   2011/12/11 L.Tong
  !!    Removed redundant paramter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Removed redundant input parameter real(double) mu
  !!   2013/03/06 17:08 dave
  !!   - Minor bug fix (length of force array)
  !!   2014/02/04 M.Arita
  !!   - Added call for update_H since it is not called at
  !!     updateIndices any longer
  !!   2017/08/29 jack baker & dave
  !!    Removed rcellx references (redundant)
  !!   2017/11/10 dave
  !!    Removed calls to dump K matrix (now done in DMMinModule)
  !!   2018/01/22 tsuyoshi (with dave)
  !!    Changes to use new atom update routines
  !!  SOURCE
  !!
  subroutine pulay_relax(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module,  only: iprint_MD, ni_in_cell, x_atom_cell,   &
                              y_atom_cell, z_atom_cell, id_glob,    &
                              atom_coord, area_general, flag_pulay_simpleStep, &
                              flag_diagonalisation, nspin, flag_LmatrixReuse, &
                              flag_SFcoeffReuse
    use group_module,   only: parts
    use minimise,       only: get_E_and_F
    use move_atoms,     only: pulayStep, velocityVerlet,            &
                              updateIndices, updateIndices3, update_atom_coord,     &
                              safemin2, update_H, update_pos_and_matrices
    use move_atoms,     only: updateL, updateLorK, updateSFcoeff
    use GenComms,       only: gsum, myid, inode, ionode, gcopy, my_barrier
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_atomic_positions, pdb_template, &
                              check_stop
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use primary_module, only: bundle
    use store_matrix,   only: dump_pos_and_matrices
    use mult_module, ONLY: matK, S_trans, matrix_scale, matL, L_trans
    use matrix_data, ONLY: Hrange, Lrange

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double), allocatable, dimension(:,:)   :: velocity
    real(double), allocatable, dimension(:,:)   :: cg
    real(double), allocatable, dimension(:)     :: x_new_pos, &
                                                   y_new_pos, &
                                                   z_new_pos
    real(double), allocatable, dimension(:,:,:) :: posnStore
    real(double), allocatable, dimension(:,:,:) :: forceStore
    real(double) :: energy0, energy1, max, g0, dE, gg, ggold, gamma, &
                    temp, KE, guess_step, step
    integer      :: i,j,k,iter,length, jj, lun, stat, npmod, pul_mx, &
                    mx_pulay, i_first, i_last, &
                    nfile, symm
    logical      :: done!, simpleStep

    step = MDtimestep
    !simpleStep = .false.
    allocate(velocity(3,ni_in_cell),STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating velocity in md_run: ", &
                       ni_in_cell, stat)
    call reg_alloc_mem(area_general, 3 * ni_in_cell, type_dbl)
    mx_pulay = 5
    allocate(posnStore(3,ni_in_cell,mx_pulay), &
             forceStore(3,ni_in_cell,mx_pulay), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(cg(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(x_new_pos(ni_in_cell), y_new_pos(ni_in_cell), &
             z_new_pos(ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating _new_pos in control: ", &
                       ni_in_cell, stat)
    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    if (myid == 0) &
         write (io_lun, fmt='(/4x,"Starting Pulay atomic relaxation"/)')
    if (myid == 0 .and. iprint_MD > 1) then
       do i = 1, ni_in_cell
          write (io_lun, 1) i, x_atom_cell(i), y_atom_cell(i), &
               z_atom_cell(i)
       end do
    end if
    posnStore = zero
    forceStore = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3 * ni_in_cell!bundle%n_prim
    if (myid == 0 .and. iprint_MD > 0) &
         write (io_lun, 2) MDn_steps, MDcgtol
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    call dump_pos_and_matrices
    iter = 1
    ggold = zero
    energy1 = energy0
    ! Store positions and forces on entry (primary set only)
    do i=1,ni_in_cell
       jj = id_glob(i)
       posnStore (1,jj,1) = x_atom_cell(i)!i)
       posnStore (2,jj,1) = y_atom_cell(i)!i)
       posnStore (3,jj,1) = z_atom_cell(i)!i)
       forceStore(:,jj,1) = tot_force(:,jj)
    end do
    do while (.not. done)
       ! Book-keeping
       npmod = mod(iter, mx_pulay)+1
       pul_mx = min(iter+1, mx_pulay)
       if (myid == 0 .and. iprint_MD > 2) &
            write(io_lun,fmt='(2x,"Pulay relaxation iteration ",i4)') iter
       x_atom_cell = zero
       y_atom_cell = zero
       z_atom_cell = zero
       tot_force = zero
       do i = 1, ni_in_cell
          jj = id_glob(i)
          if (npmod > 1) then
             x_atom_cell(i) = posnStore(1,jj,npmod-1)
             y_atom_cell(i) = posnStore(2,jj,npmod-1)
             z_atom_cell(i) = posnStore(3,jj,npmod-1)
             jj = id_glob(i)
             tot_force(:,jj) = forceStore(:,jj,npmod-1)
             !tot_force(:,jj) = forceStore(:,i,npmod-1)
          else
             x_atom_cell(i) = posnStore(1,jj,pul_mx)
             y_atom_cell(i) = posnStore(2,jj,pul_mx)
             z_atom_cell(i) = posnStore(3,jj,pul_mx)
             jj = id_glob(i)
             tot_force(:,jj) = forceStore(:,jj,pul_mx)
             !tot_force(:,jj) = forceStore(:,i,pul_mx)
          end if
       end do
       ! Take a trial step
       velocity = zero
       ! Move the atoms using velocity Verlet (why not ?)
       !call velocityVerlet(bundle,MDtimestep,temp,KE, &
       !     .false.,velocity,tot_force)
       ! Moved so we can pass cg to updateIndices3
       do j=1,ni_in_cell
          jj=id_glob(j)
          cg(1,j) = tot_force(1,jj)
          cg(2,j) = tot_force(2,jj)
          cg(3,j) = tot_force(3,jj)
          x_new_pos(j) = x_atom_cell(j)
          y_new_pos(j) = y_atom_cell(j)
          z_new_pos(j) = z_atom_cell(j)
       end do
       if(myid==0) write(io_lun,*) 'Taking trial step'
       if (flag_pulay_simpleStep) then
          if (myid == 0) write (io_lun, *) 'Step is ', step
          do i = 1, ni_in_cell
             x_atom_cell(i) = x_atom_cell(i) + step*cg(1,i)
             y_atom_cell(i) = y_atom_cell(i) + step*cg(2,i)
             z_atom_cell(i) = z_atom_cell(i) + step*cg(3,i)
          end do
          ! Output positions
          if (myid == 0 .and. iprint_MD > 1) then
             do i = 1, ni_in_cell
                write (io_lun, 1) i, x_atom_cell(i), y_atom_cell(i), &
                                  z_atom_cell(i)
             end do
          end if
          if(flag_SFcoeffReuse) then
             call update_pos_and_matrices(updateSFcoeff,cg)
          else
             call update_pos_and_matrices(updateLorK,cg)
          endif
          call update_H(fixed_potential)
          call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .false.)
          call dump_pos_and_matrices
       else
          call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
               energy1, fixed_potential, vary_mu)
       end if
       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
          end do
       end do
       if (myid == 0) &
            write (io_lun, 5) for_conv * sqrt(g0/real(ni_in_cell,double)), & !) / ni_in_cell, &
                              en_units(energy_units), d_units(dist_units)
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Exceeded number of MD steps: ",i4)') iter
       end if
       if (abs(max) < MDcgtol) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Maximum force below threshold: ",f12.5)') max
       end if
       if (.not. done) call check_stop(done, iter)
       if(done) exit

       ! Pass forces and positions down to pulayStep routine
       if(myid==0) write(io_lun,*) 'Building optimal structure'
       cg = zero
       do i = 1, ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i)!i)
          posnStore (2,jj,npmod) = y_atom_cell(i)!i)
          posnStore (3,jj,npmod) = z_atom_cell(i)!i)
          forceStore(:,jj,npmod) = tot_force(:,jj) !jj)
       end do
       call pulayStep(npmod, posnStore, forceStore, x_atom_cell, &
                      y_atom_cell, z_atom_cell, mx_pulay, pul_mx)
       ! Output positions
       if (myid == 0 .and. iprint_MD > 1) then
          do i = 1, ni_in_cell
             write (io_lun, 1) i, x_atom_cell(i), y_atom_cell(i), &
                               z_atom_cell(i)
          end do
       end if
       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff,cg)
       else
          call update_pos_and_matrices(updateLorK,cg)
       endif
       call update_H(fixed_potential)
       call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .false.)
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       call dump_pos_and_matrices

       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
          end do
       end do
       if (myid == 0) &
            write (io_lun, 5) for_conv * sqrt(g0/real(ni_in_cell,double)), & !) / ni_in_cell, &
                              en_units(energy_units), d_units(dist_units)
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Exceeded number of MD steps: ",i4)') &
                     iter
       end if
       if (abs(max) < MDcgtol) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Maximum force below threshold: ",f12.5)') &
                     max
       end if
       do i = 1, ni_in_cell!bundle%n_prim !!! ERROR ?? !!!
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i)
          posnStore (2,jj,npmod) = y_atom_cell(i)
          posnStore (3,jj,npmod) = z_atom_cell(i)
          forceStore(:,jj,npmod) = tot_force(:,jj)
       end do       
       if (.not. done) call check_stop(done, iter)
       iter = iter + 1
       dE = energy1 - energy0
       guess_step = dE*sqrt(real(ni_in_cell,double)/g0)!real(ni_in_cell,double)/sqrt(g0)
       if(myid==0) write(*,*) 'Guessed step is ',guess_step
       ! We really need a proper trust radius or something
       if(guess_step>2.0_double) guess_step = 2.0_double 
       if(guess_step>MDtimestep) then
          step=guess_step
       else
          step=MDtimestep
       end if
       energy0 = energy1
    end do
    ! Output final positions
    !    if (myid == 0) call write_positions(parts)
    deallocate(z_new_pos, y_new_pos, x_new_pos, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating _new_pos in control: ", &
                       ni_in_cell, stat)
    deallocate(cg, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating cg in control: ", &
                       ni_in_cell, stat)
    call reg_dealloc_mem(area_general, 6 * ni_in_cell, type_dbl)

1   format(4x,'Atom ',i8,' Position ',3f15.8)
2   format(4x,'Welcome to pulay_run. Doing ',i4,&
           ' steps with tolerance of ',f8.4,' ev/A')
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Energy change: ',f15.8,' ',a2)
5   format(4x,'Force Residual: ',f15.10,' ',a2,'/',a2)
6   format(4x,'Maximum force component: ',f15.8,' ',a2,'/',a2)
7   format(4x,3f15.8)
  end subroutine pulay_relax
  !!***

  !!****f* control/lbfgs *
  !!
  !!  NAME 
  !!   lbfgs
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Performs L-BFGS relaxation
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2019/12/10
  !!  MODIFICATION HISTORY
  !!   2020/01/17 smujahed
  !!    - Added call to subroutine write_xsf to write trajectory file
  !!   2022/05/27 10:05 dave
  !!    Bug fix: define g0 on all processes, not just ionode
  !!   2022/09/16 17:10 dave
  !!    Added test for forces below threshold at the start of the run
  !!  SOURCE
  !!
  subroutine lbfgs(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module,  only: iprint_MD, ni_in_cell, x_atom_cell,   &
         y_atom_cell, z_atom_cell, id_glob,    &
         atom_coord, area_general, flag_pulay_simpleStep, &
         flag_diagonalisation, nspin, flag_LmatrixReuse, &
         flag_SFcoeffReuse
    use group_module,   only: parts
    use minimise,       only: get_E_and_F
    use move_atoms,     only: pulayStep, velocityVerlet,            &
         updateIndices, updateIndices3, update_atom_coord,     &
         update_H, update_pos_and_matrices, backtrack_linemin
    use move_atoms,     only: updateL, updateLorK, updateSFcoeff
    use GenComms,       only: gsum, myid, inode, ionode, gcopy, my_barrier
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_atomic_positions, pdb_template, &
         check_stop, write_xsf
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use primary_module, only: bundle
    use store_matrix,   only: dump_pos_and_matrices
    use mult_module, ONLY: matK, S_trans, matrix_scale, matL, L_trans
    use matrix_data, ONLY: Hrange, Lrange
    use dimens,        only: r_super_x, r_super_y, r_super_z
    use md_control,    only: flag_write_xsf

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double), allocatable, dimension(:,:)   :: cg, cg_new
    real(double), allocatable, dimension(:)     :: x_new_pos, y_new_pos, z_new_pos
    real(double), allocatable, dimension(:)     :: alpha, beta, rho
    real(double), allocatable, dimension(:,:,:) :: posnStore
    real(double), allocatable, dimension(:,:,:) :: forceStore
    real(double) :: energy0, energy1, max, g0, dE, gg, ggold, gamma, &
         temp, KE, guess_step, step, test_dot
    integer      :: i,j,k,iter,length, jj, lun, stat, npmod, pul_mx, &
         i_first, i_last, &
         nfile, symm, iter_low, iter_high, this_iter
    logical      :: done

    step = MDtimestep
    allocate(posnStore(3,ni_in_cell,LBFGS_history), &
         forceStore(3,ni_in_cell,LBFGS_history), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(cg(3,ni_in_cell), cg_new(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(x_new_pos(ni_in_cell), y_new_pos(ni_in_cell), &
         z_new_pos(ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating _new_pos in control: ", &
         ni_in_cell, stat)
    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    allocate(alpha(LBFGS_history), beta(LBFGS_history), rho(LBFGS_history))
    if (myid == 0) &
         write (io_lun, fmt='(/4x,"Starting L-BFGS atomic relaxation"/)')
    if (myid == 0 .and. iprint_MD > 1) then
       do i = 1, ni_in_cell
          write (io_lun, 1) i, x_atom_cell(i), y_atom_cell(i), &
               z_atom_cell(i)
       end do
    end if
    posnStore = zero
    forceStore = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3 * ni_in_cell
    if (myid == 0 .and. iprint_MD > 0) &
         write (io_lun, 2) MDn_steps, MDcgtol, en_units(energy_units), & 
         d_units(dist_units)
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    !min_layer = min_layer - 1
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    call dump_pos_and_matrices
    call get_maxf(max)
    !min_layer = min_layer + 1
    ! Check for trivial case where forces are converged
    if (abs(max) < MDcgtol) then
       done = .true.
       if (inode==ionode) then
          write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," E: ",e16.8," dE: ",f12.8)') & 
               iter, max, energy1, en_conv*dE
          if (iprint_MD > 1) then
             write(io_lun,'(4x,"Force Residual:     ",f19.8," ",a2,"/",a2)') &
                  for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
                  d_units(dist_units)
             write(io_lun,'(4x,"Maximum force:      ",f19.8)') max
             write(io_lun,'(4x,"Force tolerance:    ",f19.8)') MDcgtol
             write(io_lun,'(4x,"Energy change:      ",f19.8," ",a2)') &
                  en_conv*dE, en_units(energy_units)
          end if
          write(io_lun, fmt='(4x,"Maximum force below threshold: ",f12.5)') max
          write(io_lun,'(2x,a,i4,a)') "GeomOpt converged in ", iter, " iterations"
       end if
       return
    end if
    iter = 0
    ggold = zero
    energy1 = energy0
    cg_new = -tot_force ! The L-BFGS is in terms of grad E
    do while (.not. done)
       ! Define g0 consistently on all processes
       g0 = dot(length, cg_new, 1, cg_new, 1)
       if (inode==ionode) then
          write(io_lun,'(/4x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," E: ",e16.8," dE: ",f12.8/)') & 
               iter, max, energy1, en_conv*dE
          if (iprint_MD > 1) then
             temp = dot(length, tot_force, 1, tot_force, 1)
             write(io_lun,'(4x,"Force Residual:     ",f19.8," ",a2,"/",a2)') &
                  for_conv*sqrt(temp/ni_in_cell), en_units(energy_units), & 
                  d_units(dist_units)
             write(io_lun,'(4x,"Maximum force:      ",f19.8)') max
             write(io_lun,'(4x,"Force tolerance:    ",f19.8)') MDcgtol
             write(io_lun,'(4x,"Energy change:      ",f19.8," ",a2)') &
                  en_conv*dE, en_units(energy_units)
             write(io_lun,'(4x,"Search direction has magnitude ",f19.8)') sqrt(g0/ni_in_cell)
          end if
       end if
       ! Book-keeping
       iter = iter + 1
       npmod = mod(iter, LBFGS_history)
       if(npmod==0) npmod = LBFGS_history
       do i=1,ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i)
          posnStore (2,jj,npmod) = y_atom_cell(i)
          posnStore (3,jj,npmod) = z_atom_cell(i)
          forceStore(:,jj,npmod) = -tot_force(:,jj)
          x_new_pos(i) = x_atom_cell(i)
          y_new_pos(i) = y_atom_cell(i)
          z_new_pos(i) = z_atom_cell(i)
          cg(:,i) = -cg_new(:,jj) ! Search downhill
       end do
       ! Set up limits for sums
       if(iter>LBFGS_history) then
          iter_low = iter-LBFGS_history+1
       else
          iter_low = 1
       end if
       if (myid == 0 .and. iprint_MD > 2) &
            write(io_lun,fmt='(2x,"L-BFGS iteration ",i4)') iter
       ! Line search
       call backtrack_linemin(cg, energy0, energy1, fixed_potential, vary_mu)
       ! Update stored position difference and force difference
       do i=1,ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i) - posnStore (1,jj,npmod)
          if(abs(posnStore(1,jj,npmod)/r_super_x)>0.7_double) posnStore(1,jj,npmod) &
               = posnStore(1,jj,npmod) &
               - nint(posnStore(1,jj,npmod)/r_super_x)*r_super_x
          posnStore (2,jj,npmod) = y_atom_cell(i) - posnStore (2,jj,npmod)
          if(abs(posnStore(2,jj,npmod)/r_super_y)>0.7_double) posnStore(2,jj,npmod) &
               = posnStore(2,jj,npmod) &
               - nint(posnStore(2,jj,npmod)/r_super_y)*r_super_y
          posnStore (3,jj,npmod) = z_atom_cell(i) - posnStore (3,jj,npmod)
          if(abs(posnStore(3,jj,npmod)/r_super_z)>0.7_double) posnStore(3,jj,npmod) &
               = posnStore(3,jj,npmod) &
               - nint(posnStore(3,jj,npmod)/r_super_z)*r_super_z
          forceStore(:,jj,npmod) = -tot_force(:,jj) - forceStore(:,jj,npmod)
          ! New search direction
       end do
       cg_new = -tot_force ! The L-BFGS is in terms of grad E
       ! Add call to write_atomic_positions and write_xsf (2020/01/17: smujahed)
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)
       ! Build q
       do i=iter, iter_low, -1
          ! Indexing
          this_iter = mod(i,LBFGS_history)
          if(this_iter==0) this_iter = LBFGS_history
          ! parameters
          rho(this_iter) = one/dot(length,forceStore(:,:,this_iter),1,posnStore(:,:,this_iter),1)
          alpha(this_iter) = rho(this_iter)*dot(length,posnStore(:,:,this_iter),1,cg_new,1)
          cg_new = cg_new - alpha(this_iter)*forceStore(:,:,this_iter)
       end do
       ! Apply preconditioning in future
       ! Build z
       do i=iter_low, iter
          ! Indexing
          this_iter = mod(i,LBFGS_history)
          if(this_iter==0) this_iter = LBFGS_history
          ! Build
          beta(this_iter) = rho(this_iter)*dot(length,forceStore(:,:,this_iter),1,cg_new,1)
          if(inode==ionode.AND.iprint_MD>2) &
               write(io_lun,fmt='(4x,"L-BFGS iter ",i2," rho, alpha, beta: ",3e15.7)') &
               this_iter, rho(this_iter), alpha(this_iter), beta(this_iter)
          cg_new = cg_new + (alpha(this_iter) - beta(this_iter))*posnStore(:,:,this_iter)
       end do
       gg = dot(length, tot_force, 1, cg_new, 1)
       if(gg>zero) then
          if(inode==ionode.AND.iprint_MD>1) then
             write(io_lun,fmt='(4x,"L-BFGS Search direction uphill; resetting to force")')
             write(io_lun,fmt='(4x,"Force residual and force.dir: ",2e15.7)') sqrt(g0/ni_in_cell), &
                  sqrt(gg/ni_in_cell)
          end if
          cg_new = -tot_force
       else if(abs(gg/g0)>ten) then
          if(inode==ionode.AND.iprint_MD>1) then
             write(io_lun,fmt='(4x,"L-BFGS Search direction anomalous; resetting to force")')
             write(io_lun,fmt='(4x,"Force residual and force.dir: ",2e15.7)') sqrt(g0/ni_in_cell), &
                  sqrt(gg/ni_in_cell)
          end if
          cg_new = -tot_force
       end if
       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
          end do
       end do
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Exceeded number of MD steps: ",i4)') iter
       end if
       if (abs(max) < MDcgtol) then
          done = .true.
          if (inode==ionode) then
             write(io_lun,'(/4x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," E: ",f16.8," dE: ",f12.8/)') & 
                  iter, max, energy1, en_conv*dE
             if (iprint_MD > 1) then
                write(io_lun,'(4x,"Force Residual:     ",f19.8," ",a2,"/",a2)') &
                     for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
                     d_units(dist_units)
                write(io_lun,'(4x,"Maximum force:      ",f19.8)') max
                write(io_lun,'(4x,"Force tolerance:    ",f19.8)') MDcgtol
                write(io_lun,'(4x,"Energy change:      ",f19.8," ",a2)') &
                     en_conv*dE, en_units(energy_units)
             end if
             write(io_lun, fmt='(4x,"Maximum force below threshold: ",f12.5)') max
             write(io_lun,'(/4x,a,i4,a)') "GeomOpt converged in ", iter, " iterations"
          end if
       end if
       if (.not. done) call check_stop(done, iter)
       if(done) exit
       if (.not. done) call check_stop(done, iter)
       dE = energy1 - energy0
       energy0 = energy1
    end do
    deallocate(cg, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating cg in control: ", &
         ni_in_cell, stat)
    call reg_dealloc_mem(area_general, 6 * ni_in_cell, type_dbl)

1   format(4x,'Atom ',i8,' Position ',3f15.8)
2   format(4x,'L-BFGS structural relaxation. Maximum of ',i4,&
         ' steps with tolerance of ',f8.4,a2,"/",a2)
  end subroutine lbfgs
  !!***

  !!****f* control/sqnm *
  !!
  !!  NAME 
  !!   sqnm
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Performs stabilised Quasi-Newton minimisation
  !!   Based on J. Chem. Phys. 142, 034112 (2015)
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2021/05/25
  !!  MODIFICATION HISTORY
  !!   2021/09/15 14:36 dave
  !!    Bug fix for force restoration when energy rises
  !!   2022/09/16 17:10 dave
  !!    Added test for forces below threshold at the start of the run
  !!   2022/12/12 11:42 dave
  !!    Added check for maximum distance moved by atoms in increase of alpha
  !!  SOURCE
  !!
  subroutine sqnm(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module,  only: iprint_MD, ni_in_cell, x_atom_cell,   &
                              y_atom_cell, z_atom_cell, id_glob,    &
                              atom_coord, area_general, flag_pulay_simpleStep, &
                              flag_diagonalisation, nspin, flag_LmatrixReuse, &
                              flag_SFcoeffReuse, flag_move_atom, min_layer
    use group_module,   only: parts
    use minimise,       only: get_E_and_F
    use move_atoms,     only: single_step, backtrack_linemin
    use GenComms,       only: myid, inode, ionode
    use GenBlas,        only: dot, syev
    use force_module,   only: tot_force
    use io_module,      only: write_atomic_positions, pdb_template, &
                              check_stop, write_xsf, print_atomic_positions, return_prefix
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use store_matrix,   only: dump_pos_and_matrices
    use dimens,        only: r_super_x, r_super_y, r_super_z
    use md_control,    only: flag_write_xsf

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double), allocatable, dimension(:,:)   :: cg, cg_new, Sij, Hij, omega, vi, ri_vec
    real(double), allocatable, dimension(:)     :: x_new_pos, y_new_pos, z_new_pos
    real(double), allocatable, dimension(:)     :: mod_dr, lambda, kappa, kappa_prime
    real(double), allocatable, dimension(:,:,:) :: posnStore, dr_tilde, dg_tilde, vi_tilde
    real(double), allocatable, dimension(:,:,:) :: forceStore
    real(double) :: energy0, energy1, max, g0, dE, gg, ggold, gamma, &
         temp, KE, guess_step, step, test_dot, lambda_max, alpha, f_dot_sd
    integer      :: i,j,k,iter,length, jj, lun, stat, npmod, pul_mx, &
         i_first, i_last, n_store, n_dim, info, n_hist, &
         nfile, symm, iter_loc, iter_high, this_iter
    logical      :: done

    character(len=12) :: subname = "sqnm: "
    character(len=120) :: prefix, prefixGO
    character(len=10) :: prefixF = "          "

    prefix = return_prefix(subname, min_layer)
    prefixGO = return_prefix("GeomOpt", min_layer)
    step = MDtimestep
    allocate(posnStore(3,ni_in_cell,LBFGS_history), &
         forceStore(3,ni_in_cell,LBFGS_history), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(cg(3,ni_in_cell), cg_new(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(x_new_pos(ni_in_cell), y_new_pos(ni_in_cell), &
         z_new_pos(ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating _new_pos in control: ", &
         ni_in_cell, stat)
    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    if (myid == 0) then
       write(io_lun, fmt='(/4x,a)') trim(prefix)//" starting SQNM atomic relaxation"
       if(iprint_MD + min_layer>0) write(io_lun, fmt='(4x,a,f8.4,1x,a2,"/",a2,a,i4)') &
            trim(prefix)//" Tolerance: ",MDcgtol,en_units(energy_units), d_units(dist_units),&
            " Maximum steps: ",MDn_steps
    end if
    posnStore = zero
    forceStore = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3 * ni_in_cell
    alpha = one
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    min_layer = min_layer - 1
    if (iprint_MD + min_layer > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.,0)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    end if
    call dump_pos_and_matrices
    call get_maxf(max)
    min_layer = min_layer + 1
    iter = 0
    if (inode==ionode) then
       write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," E: ",f16.8,x,a2," dE: ",f14.8,x,a2/)') & 
            trim(prefixGO)//" - Iter: ", iter, for_conv*max, &
            en_units(energy_units), d_units(dist_units), &
            en_conv*energy0, en_units(energy_units), en_conv*dE, en_units(energy_units)
    end if
    ! Check for trivial case where forces are converged
    if (abs(max) < MDcgtol) then
       done = .true.
       if (inode==ionode) then
          write(io_lun, fmt='(4x,a,f12.5,x,a2,"/",a2)') &
               prefixF(1:-2*min_layer)//"maximum force below threshold: ", &
               for_conv*max, en_units(energy_units), d_units(dist_units)
       end if
       return
    end if
    iter_loc = 0
    ggold = zero
    energy1 = energy0
    cg_new = -tot_force ! The L-BFGS is in terms of grad E
    g0 = dot(length,cg_new,1,cg_new,1)
    if (inode==ionode .and. iprint_MD + min_layer > 2) then
       write(io_lun,'(4x,a,f19.8)') &
            trim(prefix)//" search direction has magnitude ", sqrt(g0/ni_in_cell)
    end if
    do while (.not. done)
       ! Book-keeping
       iter = iter + 1
       iter_loc = iter_loc + 1
       npmod = mod(iter_loc, LBFGS_history)
       if(npmod==0) npmod = LBFGS_history
       do i=1,ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i)
          posnStore (2,jj,npmod) = y_atom_cell(i)
          posnStore (3,jj,npmod) = z_atom_cell(i)
          forceStore(:,jj,npmod) = -tot_force(:,jj)
          x_new_pos(i) = x_atom_cell(i)
          y_new_pos(i) = y_atom_cell(i)
          z_new_pos(i) = z_atom_cell(i)
          cg(:,i) = -cg_new(:,jj) ! Search downhill
       end do
       ! Set up limits for sums
       min_layer = min_layer - 1
       ! Line search
       if(iter==1) then
          call backtrack_linemin(cg, energy0, energy1, fixed_potential, vary_mu)
       else
          call single_step(cg, energy0, energy1, fixed_potential, vary_mu)
       end if
       if(energy1>energy0) then
          if(inode==ionode.AND.iprint_MD>2) &
               write(io_lun,fmt='(4x,a)') trim(prefix)//" energy rise: resetting history"
          cg_new = -tot_force
          do i=1,ni_in_cell
             jj = id_glob(i)
             cg(:,i) = -cg_new(:,jj) ! Search downhill
          end do
          call backtrack_linemin(cg, energy0, energy1, fixed_potential, vary_mu)
          !call single_step(cg, energy0, energy1, fixed_potential, vary_mu)
          if(energy1>energy0) call cq_abort("Energy rise twice in succession: check SCF and other tolerances")
          npmod = 1
          iter_loc = 1
          ! In the original paper, this is alpha = alpha/2 but heuristically this seems better
          ! DRB 2021/09/15
          alpha = one
       endif
       min_layer = min_layer + 1
       ! Update stored position difference and force difference
       do i=1,ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i) - posnStore (1,jj,npmod)
          if(abs(posnStore(1,jj,npmod)/r_super_x)>0.7_double) posnStore(1,jj,npmod) &
               = posnStore(1,jj,npmod) &
               - nint(posnStore(1,jj,npmod)/r_super_x)*r_super_x
          posnStore (2,jj,npmod) = y_atom_cell(i) - posnStore (2,jj,npmod)
          if(abs(posnStore(2,jj,npmod)/r_super_y)>0.7_double) posnStore(2,jj,npmod) &
               = posnStore(2,jj,npmod) &
               - nint(posnStore(2,jj,npmod)/r_super_y)*r_super_y
          posnStore (3,jj,npmod) = z_atom_cell(i) - posnStore (3,jj,npmod)
          if(abs(posnStore(3,jj,npmod)/r_super_z)>0.7_double) posnStore(3,jj,npmod) &
               = posnStore(3,jj,npmod) &
               - nint(posnStore(3,jj,npmod)/r_super_z)*r_super_z
          forceStore(:,jj,npmod) = -tot_force(:,jj) - forceStore(:,jj,npmod)
          ! New search direction
       end do
       n_store = min(iter_loc,LBFGS_history) ! Number of stored states
       if(inode==ionode.AND.iprint_MD + min_layer > 2) &
            write(io_lun,fmt='(4x,a,i3)') trim(prefix)//" number of stored histories: ",n_store
       allocate(mod_dr(n_store),Sij(n_store,n_store),lambda(n_store),&
            omega(n_store,n_store))
       mod_dr = zero
       ! Normalise dR and dg
       do i=1,n_store
          mod_dr(i) = sqrt(dot(length,posnStore(:,:,i),1,posnStore(:,:,i),1))
       end do
       posnStore(:,:,npmod) = posnStore(:,:,npmod)/mod_dr(npmod)
       forceStore(:,:,npmod) = forceStore(:,:,npmod)/mod_dr(npmod)
       cg_new = -tot_force ! The L-BFGS is in terms of grad E
       ! Output positions
       if (myid == 0 .and. iprint_md > 1) then
          call print_atomic_positions
       end if
       ! Add call to write_atomic_positions and write_xsf (2020/01/17: smujahed)
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)
       ! Build significant subspace
       Sij = zero
       omega = zero
       do i = 1, n_store
          do j=i,n_store
             Sij(j,i) = dot(length,posnStore(:,:,j),1,posnStore(:,:,i),1)
             if(j>i) Sij(i,j) = Sij(j,i)
          end do
       end do
       ! Solve for eigenvectors of Sij
       omega = Sij
       if(n_store>1) then
          call syev('U',n_store,omega,n_store,lambda,info)
          if(info<0) call cq_abort("Error in SQNM calling dsyev: ",info)
          if(info>0.and.inode==ionode) &
               write(io_lun,fmt='(4x,a,i4)') trim(prefix)//" issue in SQNM; dsyev returned ",info
       else
          lambda = one
          omega = one
       end if
       if(inode==ionode.AND.iprint_MD + min_layer > 3) then
          write(io_lun,fmt='(4x,a)') trim(prefix)//" eigenvalues of Sij: "
          do i = 1,n_store
             write(io_lun,fmt='(6x,f7.4)') lambda(i)
          end do
       end if
       lambda_max = maxval(lambda)
       n_dim = n_store
       do i=1, n_store
          if(lambda(i)/lambda_max<1.0e-4_double) then
             n_dim = n_dim - 1
          end if
       end do
       if(inode==ionode.AND.iprint_MD + min_layer > 3) &
            write(io_lun,fmt='(4x,a,i3)') trim(prefix)//" number of eigenstates kept: ",n_dim
       ! Build dr_tilde and dg_tilde
       allocate(dr_tilde(3,ni_in_cell,n_dim), dg_tilde(3,ni_in_cell,n_dim))
       dr_tilde = zero
       dg_tilde = zero
       do i=1,n_dim
          do j=1,n_store
             dr_tilde(:,:,i) = dr_tilde(:,:,i) + omega(j,i+(n_store-n_dim))*posnStore(:,:,j)
             dg_tilde(:,:,i) = dg_tilde(:,:,i) + omega(j,i+(n_store-n_dim))*forceStore(:,:,j)
          end do
          dr_tilde(:,:,i) = dr_tilde(:,:,i)/sqrt(lambda(i+(n_store-n_dim)))
          dg_tilde(:,:,i) = dg_tilde(:,:,i)/sqrt(lambda(i+(n_store-n_dim)))
       end do
       ! Construct approximate Hessian projected onto subspace
       allocate(Hij(n_dim,n_dim))
       Hij = zero
       do i=1,n_dim
          do j=i,n_dim
             Hij(j,i) = half*(dot(length,dr_tilde(:,:,i),1,dg_tilde(:,:,j),1) + &
                  dot(length,dr_tilde(:,:,j),1,dg_tilde(:,:,i),1))
             if(j>i) Hij(i,j) = Hij(j,i)
          end do
       end do
       ! Find eigenvectors
       allocate(kappa(n_dim),vi(n_dim,n_dim),kappa_prime(n_dim))
       kappa = zero
       kappa_prime = zero
       vi = Hij
       if(n_dim>1) then
          call syev('U',n_dim,vi,n_dim,kappa,info)
          if(info<0) call cq_abort("Error in SQNM calling dsyev: ",info)
          if(info>0.and.inode==ionode) &
               write(io_lun,fmt='(4x,a,i4)') trim(prefix)//" issue in SQNM; dsyev returned ",info
       else
          kappa = one !Hij(1,1)
          vi = one
       end if
       if(inode==ionode.AND.iprint_MD + min_layer > 3) then
          write(io_lun,fmt='(4x,a)') trim(prefix)//" kappa: "
          do i=1,n_dim
             write(io_lun,fmt='(6x,f7.4)') kappa(i)
          end do
       end if
       ! Build v tilde
       allocate(vi_tilde(3,ni_in_cell,n_dim),ri_vec(3,ni_in_cell))
       vi_tilde = zero
       do i=1,n_dim
          ri_vec = zero
          do j=1,n_dim
             vi_tilde(:,:,i) = vi_tilde(:,:,i) + vi(j,i)*dr_tilde(:,:,j)
             ri_vec(:,:) = ri_vec(:,:) + vi(j,i)*dg_tilde(:,:,j)
          end do
          ri_vec(:,:) = ri_vec(:,:) - kappa(i)*vi_tilde(:,:,i)
          !kappa_prime(i) = sqrt(0.0025_double + kappa(i)*kappa(i))
          kappa_prime(i) = sqrt(dot(length,ri_vec,1,ri_vec,1) + kappa(i)*kappa(i))
       end do
       if(inode==ionode.AND.iprint_MD>3) then
          write(io_lun,fmt='(4x,a)') trim(prefix)//" kappa prime: "
          do i=1,n_dim
             write(io_lun,fmt='(6x,f7.4)') kappa_prime(i)
          end do
       end if
       ! Build preconditioned search
       cg_new = -alpha*tot_force
       do i=1,n_dim
          cg_new = cg_new - (one/kappa_prime(i) - alpha)*dot(length,tot_force,1,vi_tilde(:,:,i),1)*vi_tilde(:,:,i)
       end do
       ! Zero search direction for fixed atoms and find maximum force
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
             if (.not. flag_move_atom(k,i)) then
                cg_new(k,i) = zero
             end if
          end do
       end do
       gg = dot(length, cg_new, 1, cg_new, 1)
       ! Analyse forces
       if(inode==ionode.AND.iprint_MD + min_layer > 1) write(io_lun,'(4x,a,f19.8)') &
            trim(prefix)//" search direction has magnitude ", sqrt(gg/ni_in_cell)
       g0 = dot(length, tot_force, 1, tot_force, 1)
       f_dot_sd = dot(length,cg_new,1,-tot_force,1)/sqrt(g0*gg)
       if(inode==ionode.AND.iprint_MD + min_layer > 2) &
            write(io_lun,fmt='(4x,a,f12.7)') trim(prefix)//" force dot search direction: ",f_dot_sd
       ! Now adjust alpha
       if(f_dot_sd>0.2_double) then
          alpha = alpha*1.1_double
          max = maxval(abs(cg_new))
          ! Hard limit on alpha so that atoms are not moved by more than sqnm_trust_step Bohr
          if(max*alpha>sqnm_trust_step) alpha = sqnm_trust_step/max
       else
          alpha = alpha*0.85_double
       end if
       if(inode==ionode.AND.iprint_MD + min_layer > 2) &
            write(io_lun,fmt='(4x,a,f9.5)') trim(prefix)//" alpha set to: ",alpha
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
          end do
       end do
       ! Output and energy changes
       dE = energy1 - energy0
       energy0 = energy1
       if (inode==ionode) then
          write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," E: ",f16.8,x,a2," dE: ",f14.8,x,a2/)') & 
               trim(prefixGO)//" - Iter: ", iter, for_conv*max, &
               en_units(energy_units), d_units(dist_units), &
               en_conv*energy1, en_units(energy_units), en_conv*dE, en_units(energy_units)
          if (iprint_MD + min_layer > 1) then
             g0 = dot(length, tot_force, 1, tot_force, 1)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force Residual:     ", &
                  for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Maximum force:      ", &
                  for_conv*max, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force tolerance:    ", &
                  for_conv*MDcgtol, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2)') prefixF(1:-2*min_layer)//"Energy change:      ", &
                  en_conv*dE, en_units(energy_units)
          end if
       end if
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write(io_lun, fmt='(4x,a,i6)') &
                     trim(prefix)//"Exceeded number of MD steps: ",iter
       end if
       if (abs(max) < MDcgtol) then
          done = .true.
          if (inode==ionode) then
             write(io_lun,'(/4x,a,i4,a)') trim(prefixGO)//" converged in ", iter, " iterations"
             write(io_lun, fmt='(4x,a,f12.5,x,a2,"/",a2)') &
                  trim(prefix)//" maximum force below threshold: ",&
                  for_conv*max, en_units(energy_units), d_units(dist_units)
          end if
       end if
       deallocate(dr_tilde,dg_tilde,Hij,kappa,vi,kappa_prime,vi_tilde,ri_vec)
       deallocate(mod_dr,Sij,lambda,omega)
       if (.not. done) call check_stop(done, iter)
       if(done) exit
    end do ! .not. done i.e. until max iterations or force tolerance reached
    deallocate(cg, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating cg in control: ", &
         ni_in_cell, stat)
    call reg_dealloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    deallocate(posnStore, forceStore, cg_new, x_new_pos, y_new_pos, z_new_pos)
  end subroutine sqnm
  !!***

  !!****f* control/cell_sqnm *
  !!
  !!  NAME
  !!   cell_sqnm
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs stabilised Quasi-Newton minimisation to optimise
  !!   simulation cell size.  Based on algorithm for ionic
  !!   optimisation in J. Chem. Phys. 142, 034112 (2015) and updated
  !!   paper arXiv 2206.07339
  !!
  !!   Note that the arXiv paper works in terms of dE/dA_{lat} which
  !!   is converted from stress as dE/dA = stress/{a|b|c}
  !!
  !!  ** NB not yet functional: do not use 2022/09/16 16:52 dave **
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2022/08/11
  !!  MODIFICATION HISTORY
  !!   2022/08/17 15:18 dave
  !!    Introduced scaling to improve conditioning in arxiv/2206.07339
  !!  SOURCE
  !!
  subroutine cell_sqnm(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, only: iprint_gen, ni_in_cell, x_atom_cell,  &
                             y_atom_cell, z_atom_cell, id_glob,    &
                             atom_coord, rcellx, rcelly, rcellz,   &
                             area_general, iprint_MD,              &
                             IPRINT_TIME_THRES1, cell_en_tol,      &
                             cell_constraint_flag, cell_stress_tol, min_layer
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: single_step_cell, backtrack_linemin_cell, enthalpy, enthalpy_tolerance
    use GenComms,      only: gsum, myid, inode, ionode
    use GenBlas,       only: dot, syev
    use force_module,  only: stress, tot_force
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop, write_xsf, leqi
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use dimens, ONLY: r_super_x, r_super_y, r_super_z
    use store_matrix,  only: dump_pos_and_matrices
    use md_control,    only: target_pressure, flag_write_xsf

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double), allocatable, dimension(:,:)   :: omega, vi, Sij, Hij
    real(double), dimension(3)     :: cg, cg_new, orcell
    real(double), allocatable, dimension(:)     :: mod_dr, lambda, kappa, kappa_prime, ri_vec
    real(double), allocatable, dimension(:,:) :: posnStore, dr_tilde, dg_tilde, vi_tilde
    real(double), allocatable, dimension(:,:) :: forceStore
    real(double) :: energy0, energy1, max, g0, dE, gg, ggold, gamma, &
         temp, KE, guess_step, step, test_dot, lambda_max, alpha, f_dot_sd, &
         enthalpy0, enthalpy1, dH, press
    integer      :: i,j,k,iter,length, jj, lun, stat, npmod, pul_mx, &
         i_first, i_last, n_store, n_dim, info, n_hist, &
         nfile, symm, iter_loc, iter_high, this_iter
    logical      :: done
    real(double) :: search_dir_x, search_dir_y,&
         search_dir_z, stressx, stressy, stressz, RMSstress, newRMSstress,&
         dRMSstress, search_dir_mean, mean_stress, max_stress, &
         stress_diff, volume, stress_target, orcellx, orcelly, orcellz, wscal

    ! Store original cell size
    orcellx = rcellx
    orcelly = rcelly
    orcellz = rcellz
    orcell(1) = rcellx
    orcell(2) = rcelly
    orcell(3) = rcellz
    ! Scaling: w = 2 Bohr x sqrt(Natoms)
    wscal = two*sqrt(real(ni_in_cell,double))
    step = MDtimestep
    allocate(posnStore(3,LBFGS_history), &
         forceStore(3,LBFGS_history), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating in cell_sqnm: ", ni_in_cell, stat)
    call reg_alloc_mem(area_general, 6 * LBFGS_history, type_dbl)
    if (myid == 0) &
         write(io_lun, fmt='(/4x,"Starting SQNM cell relaxation"/)')
    posnStore = zero
    forceStore = zero
    search_dir_x = zero
    search_dir_y = zero
    search_dir_z = zero
    search_dir_mean = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3
    alpha = one
    if (myid == 0 .and. iprint_MD > 0) &
         write(io_lun, fmt='(4x,"SQNM cell optimisation. Maximum of ",i4, " steps with tolerance of ",&
         &f8.4,a2,"/",a2)') MDn_steps, MDcgtol, en_units(energy_units), d_units(dist_units)
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    min_layer = min_layer - 1
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.,0)
    min_layer = min_layer + 1
    call dump_pos_and_matrices
    call get_maxf(max)
    iter = 0
    press = target_pressure/HaBohr3ToGPa
    ! Stress tolerance in Ha/Bohr3
    stress_target = cell_stress_tol/HaBohr3ToGPa
    enthalpy0 = enthalpy(energy0, press)
    enthalpy1 = enthalpy0
    dH = zero
    max_stress = zero
    volume = rcellx*rcelly*rcellz
    do i=1,3
       stress_diff = abs(press*volume + stress(i,i))/volume
       if (stress_diff > max_stress) max_stress = stress_diff
    end do
    if (inode==ionode) then
      write(io_lun,'(/4x,"GeomOpt - Iter: ",i4," MaxStr: ",f12.8," H: ",f16.8," dH: ",f12.8/)') &
           0, max_stress, enthalpy0, zero
    end if
    iter_loc = 0
    ggold = zero
    energy1 = energy0
    if (leqi(cell_constraint_flag, 'volume')) then
       cg_new(1) = third*(stress(1,1)+stress(2,2)+stress(3,3)) - press*volume
    else
       !do i = 1,3
       cg_new(1) = (stress(1,1) - press*volume)*orcellx/(rcellx*wscal)
       cg_new(2) = (stress(2,2) - press*volume)*orcelly/(rcelly*wscal)
       cg_new(3) = (stress(3,3) - press*volume)*orcellz/(rcellz*wscal)
       !end do
    end if
    if (inode==ionode .and. iprint_MD > 1) then
       g0 = dot(length,cg_new,1,cg_new,1)
       write(io_lun,'(4x,"Search direction has magnitude ",f19.8)') sqrt(g0/three)
    end if
    do while (.not. done)
       ! Book-keeping
       iter = iter + 1
       iter_loc = iter_loc + 1
       npmod = mod(iter_loc, LBFGS_history)
       if(npmod==0) npmod = LBFGS_history
       volume = rcellx*rcelly*rcellz
       posnStore (1,npmod) = wscal*rcellx/orcellx
       posnStore (2,npmod) = wscal*rcelly/orcelly
       posnStore (3,npmod) = wscal*rcellz/orcellz
       if (leqi(cell_constraint_flag, 'volume')) then
          forceStore(1,npmod) = -third*(stress(1,1)+stress(2,2)+stress(3,3)) + press*volume
          cg(1) = -cg_new(1)
       else
          forceStore(1,npmod) = (-stress(1,1) + press*volume)*orcellx/(wscal*rcellx)
          forceStore(2,npmod) = (-stress(2,2) + press*volume)*orcelly/(wscal*rcelly)
          forceStore(3,npmod) = (-stress(3,3) + press*volume)*orcellz/(wscal*rcellz)
          cg(1) = -cg_new(1)*wscal/orcellx
          cg(2) = -cg_new(2)*wscal/orcelly
          cg(3) = -cg_new(3)*wscal/orcellz
       end if
       write(*,*) 'Search direction: ',cg
       ! Set up limits for sums
       if (myid == 0 .and. iprint_MD > 2) &
            write(io_lun,fmt='(2x,"SQNM iteration ",i4)') iter
       ! Take a step downhill
       if(iter==1) then
          call backtrack_linemin_cell(cg, press, enthalpy0, enthalpy1, fixed_potential, vary_mu)
       else
          call single_step_cell(cg, press, enthalpy0, enthalpy1, fixed_potential, vary_mu)
       end if
       if(enthalpy1>enthalpy0) then
          if(inode==ionode.AND.iprint_MD>1) write(io_lun,fmt='(4x,"Energy rise: resetting history")')
          cg_new(1) = (stress(1,1) - press*volume)*orcellx/(wscal*rcellx)
          cg_new(2) = (stress(2,2) - press*volume)*orcelly/(wscal*rcelly)
          cg_new(3) = (stress(3,3) - press*volume)*orcellz/(wscal*rcellz)
          cg(1) = -cg_new(1)*wscal/orcellx
          cg(2) = -cg_new(2)*wscal/orcelly
          cg(3) = -cg_new(3)*wscal/orcellz
          !do i=1,3
          !cg_new(i) = (stress(i,i) - press*volume)*orcell(i)/wscal
          !cg(i) = -cg_new(i)*wscal/orcell(i) ! Search downhill
          !end do
          call backtrack_linemin_cell(cg, press, enthalpy0, enthalpy1, fixed_potential, vary_mu)
          if(enthalpy1>enthalpy0) call cq_abort("Energy rise twice in succession: check SCF and other tolerances")
          npmod = 1
          iter_loc = 1
          ! In the original paper, this is alpha = alpha/2 but heuristically this seems better
          ! DRB 2021/09/15
          alpha = one
       endif
       ! Update stored position difference and force difference
       posnStore (1,npmod) = rcellx*wscal/orcellx - posnStore (1,npmod)
       posnStore (2,npmod) = rcelly*wscal/orcelly - posnStore (2,npmod)
       posnStore (3,npmod) = rcellz*wscal/orcellz - posnStore (3,npmod)
       forceStore(1,npmod) = (-stress(1,1) + press*volume)*orcellx/(wscal*rcellx) - forceStore(1,npmod)
       forceStore(2,npmod) = (-stress(2,2) + press*volume)*orcelly/(wscal*rcelly) - forceStore(2,npmod)
       forceStore(3,npmod) = (-stress(3,3) + press*volume)*orcellz/(wscal*rcellz) - forceStore(3,npmod)
       n_store = min(iter_loc,LBFGS_history) ! Number of stored states
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,fmt='(4x,"Number of stored histories ",i3)') n_store
       allocate(mod_dr(n_store),Sij(n_store,n_store),lambda(n_store),&
            omega(n_store,n_store))
       mod_dr = zero
       ! Normalise dR and dg
       do i=1,n_store
          mod_dr(i) = sqrt(dot(length,posnStore(:,i),1,posnStore(:,i),1))
       end do
       posnStore(:,npmod) = posnStore(:,npmod)/mod_dr(npmod)
       forceStore(:,npmod) = forceStore(:,npmod)/mod_dr(npmod)
       ! Add call to write_atomic_positions and write_xsf (2020/01/17: smujahed)
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)
       ! Build significant subspace
       Sij = zero
       omega = zero
       do i = 1, n_store
          do j=i,n_store
             Sij(j,i) = dot(length,posnStore(:,j),1,posnStore(:,i),1)
             if(j>i) Sij(i,j) = Sij(j,i)
          end do
       end do
       ! Solve for eigenvectors of Sij
       omega = Sij
       if(n_store>1) then
          call syev('U',n_store,omega,n_store,lambda,info)
          if(info<0) call cq_abort("Error in SQNM calling dsyev: ",info)
          if(info>0.and.inode==ionode) write(io_lun,fmt='(4x,"Possible error in SQNM; dsyev returned ",i4)') info
       else
          lambda = one
          omega = one
       end if
       if(inode==ionode.AND.iprint_MD>2) then
          write(io_lun,fmt='(4x,"Eigenvalues of Sij: ",(f7.4))') lambda
       end if
       lambda_max = maxval(lambda)
       n_dim = n_store
       do i=1, n_store
          if(lambda(i)/lambda_max<1.0e-4_double) then
             n_dim = n_dim - 1
          end if
       end do
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,fmt='(4x,"Number of eigenstates kept: ",i3)') n_dim
       ! Build dr_tilde and dg_tilde
       allocate(dr_tilde(3,n_dim), dg_tilde(3,n_dim))
       dr_tilde = zero
       dg_tilde = zero
       do i=1,n_dim
          do j=1,n_store
             dr_tilde(:,i) = dr_tilde(:,i) + omega(j,i+(n_store-n_dim))*posnStore(:,j)
             dg_tilde(:,i) = dg_tilde(:,i) + omega(j,i+(n_store-n_dim))*forceStore(:,j)
          end do
          dr_tilde(:,i) = dr_tilde(:,i)/sqrt(lambda(i+(n_store-n_dim)))
          dg_tilde(:,i) = dg_tilde(:,i)/sqrt(lambda(i+(n_store-n_dim)))
       end do
       ! Construct approximate Hessian projected onto subspace
       allocate(Hij(n_dim,n_dim))
       Hij = zero
       do i=1,n_dim
          do j=i,n_dim
             Hij(j,i) = half*(dot(length,dr_tilde(:,i),1,dg_tilde(:,j),1) + &
                  dot(length,dr_tilde(:,j),1,dg_tilde(:,i),1))
             if(j>i) Hij(i,j) = Hij(j,i)
          end do
       end do
       ! Find eigenvectors
       allocate(kappa(n_dim),vi(n_dim,n_dim),kappa_prime(n_dim))
       kappa = zero
       kappa_prime = zero
       vi = Hij
       if(n_dim>1) then
          call syev('U',n_dim,vi,n_dim,kappa,info)
          if(info<0) call cq_abort("Error in SQNM calling dsyev: ",info)
          if(info>0.and.inode==ionode) write(io_lun,fmt='(4x,"Possible error in SQNM; dsyev returned ",i4)') info
       else
          kappa = one !Hij(1,1)
          vi = one
       end if
       if(inode==ionode.AND.iprint_MD>3) write(io_lun,fmt='(4x,"Kappa: ",(f7.4))') kappa
       ! Build v tilde
       allocate(vi_tilde(3,n_dim),ri_vec(3))
       vi_tilde = zero
       do i=1,n_dim
          ri_vec = zero
          do j=1,n_dim
             vi_tilde(:,i) = vi_tilde(:,i) + vi(j,i)*dr_tilde(:,j)
             ri_vec(:) = ri_vec(:) + vi(j,i)*dg_tilde(:,j)
          end do
          ri_vec(:) = ri_vec(:) - kappa(i)*vi_tilde(:,i)
          !kappa_prime(i) = sqrt(0.0025_double + kappa(i)*kappa(i))
          kappa_prime(i) = sqrt(dot(length,ri_vec,1,ri_vec,1) + kappa(i)*kappa(i))
       end do
       if(inode==ionode.AND.iprint_MD>3) write(io_lun,fmt='(4x,"Kappa prime: ",(f7.4))') kappa_prime
       ! Build preconditioned search
       if (leqi(cell_constraint_flag, 'volume')) then
          cg_new(1) = third*(stress(1,1)+stress(2,2)+stress(3,3)) - press*volume
       else
          cg_new(1) = (stress(1,1) - press*volume)*orcellx/(wscal*rcellx)
          cg_new(2) = (stress(2,2) - press*volume)*orcelly/(wscal*rcelly)
          cg_new(3) = (stress(3,3) - press*volume)*orcellz/(wscal*rcellz)
          !do i = 1,3
          !   cg_new(i) = (stress(i,i) - press*volume)*orcell(i)/wscal
          !end do
       end if
       cg_new = alpha*cg_new
       do i=1,n_dim
          temp = zero
          !do j=1,3
          temp = temp + (stress(1,1)-press*volume)*orcell(1)*vi_tilde(1,i)/(wscal*rcellx)
          temp = temp + (stress(2,2)-press*volume)*orcell(2)*vi_tilde(2,i)/(wscal*rcelly)
          temp = temp + (stress(3,3)-press*volume)*orcell(3)*vi_tilde(3,i)/(wscal*rcellz)
          !end do
          do j=1,3
             cg_new(j) = cg_new(j) - (one/kappa_prime(i) - alpha)*temp*vi_tilde(j,i)
          end do
       end do
       ! Zero search direction for fixed atoms and find maximum force
       gg = dot(length, cg_new, 1, cg_new, 1)
       ! Analyse forces
       g0 = stress(1,1)*stress(1,1) + stress(2,2)*stress(2,2) + stress(3,3)*stress(3,3)
       temp = -cg_new(1)*stress(1,1) - cg_new(2)*stress(2,2) - cg_new(3)*stress(3,3)
       f_dot_sd = temp/sqrt(g0*gg)
       if(inode==ionode.AND.iprint_MD>2) &
            write(io_lun,fmt='(4x,"Dot product of search direction and force: ",f12.7)') f_dot_sd
       ! Now adjust alpha
       if(f_dot_sd>0.2_double) then
          alpha = alpha*1.1_double
       else
          alpha = alpha*0.85_double
       end if
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,fmt='(4x,"Alpha set to: ",f9.5)') alpha
       max_stress = zero
       volume = rcellx*rcelly*rcellz
       do i=1,3
          stress_diff = abs(press*volume + stress(i,i))/volume
          if (stress_diff > max_stress) max_stress = stress_diff
       end do
       dH = enthalpy1 - enthalpy0
       newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
            (stress(2,2)*stress(2,2)) + &
            (stress(3,3)*stress(3,3)))/3)
       dRMSstress = (RMSstress - newRMSstress)/volume
       enthalpy0 = enthalpy1
       ! Check exit criteria
       volume = rcellx*rcelly*rcellz
       max_stress = zero
       do i=1,3
          stress_diff = abs(press*volume + stress(i,i))/volume
          if (stress_diff > max_stress) max_stress = stress_diff
       end do
       if (inode==ionode) then
          write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxStr: ",f12.8," H: ",f16.8," dH: ",f12.8)') &
               iter, max_stress*volume, enthalpy1, en_conv*dH
          if (iprint_MD > 1) then
             write(io_lun,'(4x,"Maximum stress         ",e14.6," Ha/Bohr**3")') max_stress
             write(io_lun,'(4x,"Simulation cell volume ",e14.6," Bohr**3")') volume
             write(io_lun,'(4x,"Maximum stress         ",f14.6," GPa")') &
                  max_stress*HaBohr3ToGPa
             write(io_lun,'(4x,"Stress tolerance:      ",f14.6," GPa")') &
                  cell_stress_tol
             !write(io_lun,'(4x,"Change in stress:      ",e14.6," ",a2)') max_stress*volume, &
             !     en_units(energy_units)
             write(io_lun,'(4x,"Enthalpy change:       ",e14.6," ",a2)') &
                  en_conv*dH, en_units(energy_units)
             write(io_lun,'(4x,"Enthalpy tolerance:    ",e14.6," ",a2)') &
                  en_conv*enthalpy_tolerance, en_units(energy_units)
          else if (iprint_MD > 0) then
             write(io_lun, fmt='(4x,"Enthalpy change: ",f15.8," ",a2)') en_conv*dH, en_units(energy_units)
             write(io_lun, fmt='(4x,"RMS Stress change: ",f15.8," ",a2)') dRMSstress, "Ha"
             write(io_lun,'(4x,"Maximum stress: ",f15.8," ",a2)') max_stress*volume, &
                  en_units(energy_units)
          end if
       end if
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write(io_lun, fmt='(4x,"Exceeded number of MD steps: ",i4)') iter
       end if
       if (abs(dH)<enthalpy_tolerance .and. max_stress < stress_target) then
          if (inode==ionode) &
               write(io_lun,'(2x,a,i4,a)') "GeomOpt converged in ", &
               iter, " iterations"
          done = .true.
          if (myid == 0 .and. iprint_md > 0) &
               write(io_lun, fmt='(4x,"Enthalpy change below threshold: ",f19.8," ",a2)') &
               dH*en_conv, en_units(energy_units)
          if(myid==0) write(io_lun, fmt='(4x,"Maximum stress below threshold:   ",f19.8," GPa")') &
               max_stress*HaBohr3ToGPa
       end if
       deallocate(dr_tilde,dg_tilde,Hij,kappa,vi,kappa_prime,vi_tilde,ri_vec)
       deallocate(mod_dr,Sij,lambda,omega)
       if (.not. done) call check_stop(done, iter)
       if(done) exit
    end do ! .not. done i.e. until max iterations or force tolerance reached
    deallocate(posnStore, forceStore)
  end subroutine cell_sqnm
  !!***

  !!****f* control/full_sqnm *
  !!
  !!  NAME 
  !!   sqnm
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Performs stabilised Quasi-Newton minimisation
  !!   optimisation for ionic positions and simulation cell
  !!   Based on J. Chem. Phys. 142, 034112 (2015) and
  !!   arXiv/2206.07339
  !!
  !!  ** NB not yet functional: do not use 2022/09/16 16:52 dave **
  !!
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2022/08/19
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine full_sqnm(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module,  only: iprint_MD, ni_in_cell, x_atom_cell,   &
         y_atom_cell, z_atom_cell, id_glob,    &
         atom_coord, area_general, flag_pulay_simpleStep, &
         flag_diagonalisation, nspin, flag_LmatrixReuse, &
         flag_SFcoeffReuse, flag_move_atom, rcellx, rcelly, rcellz, &
         cell_constraint_flag, cell_stress_tol, min_layer
    use group_module,   only: parts
    use minimise,       only: get_E_and_F
    use move_atoms,     only: pulayStep, velocityVerlet,            &
         updateIndices, updateIndices3, update_atom_coord,     &
         update_H, update_pos_and_matrices, single_step_full, enthalpy
    use move_atoms,     only: updateL, updateLorK, updateSFcoeff, backtrack_linemin_full
    use GenComms,       only: gsum, myid, inode, ionode, gcopy, my_barrier
    use GenBlas,        only: dot, syev
    use force_module,   only: tot_force, stress
    use io_module,      only: write_atomic_positions, pdb_template, &
         check_stop, write_xsf
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use primary_module, only: bundle
    use store_matrix,   only: dump_pos_and_matrices
    use mult_module, ONLY: matK, S_trans, matrix_scale, matL, L_trans
    use matrix_data, ONLY: Hrange, Lrange
    use dimens,        only: r_super_x, r_super_y, r_super_z
    use md_control,    only: flag_write_xsf, target_pressure

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double), allocatable, dimension(:,:)   :: cg, cg_new, Sij, Hij, omega, vi, ri_vec
    real(double), allocatable, dimension(:)     :: x_new_pos, y_new_pos, z_new_pos
    real(double), allocatable, dimension(:)     :: mod_dr, lambda, kappa, kappa_prime
    real(double), allocatable, dimension(:,:,:) :: posnStore, dr_tilde, dg_tilde, vi_tilde
    real(double), allocatable, dimension(:,:,:) :: forceStore
    real(double), dimension(3)     :: orcell, rcell
    real(double) :: energy0, energy1, max, g0, dE, gg, ggold, gamma, &
         temp, KE, guess_step, step, test_dot, lambda_max, alpha, f_dot_sd, &
         orcellx, orcelly, orcellz, wscal, press, volume, enthalpy0, enthalpy1, dH, &
         max_stress, stress_diff, stress_target
    integer      :: i,j,k,iter,length, jj, lun, stat, npmod, pul_mx, &
         i_first, i_last, n_store, n_dim, info, n_hist, &
         nfile, symm, iter_loc, iter_high, this_iter
    logical      :: done

    ! Store original cell size
    orcellx = rcellx
    orcelly = rcelly
    orcellz = rcellz
    orcell(1) = rcellx
    orcell(2) = rcelly
    orcell(3) = rcellz
    ! Scaling: w = 2 Bohr x sqrt(Natoms)
    wscal = two*sqrt(real(ni_in_cell,double))
    step = MDtimestep
    allocate(posnStore(3,ni_in_cell+1,LBFGS_history), &
         forceStore(3,ni_in_cell+1,LBFGS_history), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(cg(3,ni_in_cell+1), cg_new(3,ni_in_cell+1), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ", ni_in_cell, stat)
    allocate(x_new_pos(ni_in_cell), y_new_pos(ni_in_cell), &
         z_new_pos(ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating _new_pos in control: ", &
         ni_in_cell, stat)
    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    if (myid == 0) &
         write(io_lun, fmt='(/4x,"Starting SQNM cell/atomic relaxation"/)')
    if (myid == 0 .and. iprint_MD > 1) then
       do i = 1, ni_in_cell
          write(io_lun, fmt='(4x,"Atom ",i8," Position ",3f15.8)') i, x_atom_cell(i), y_atom_cell(i), &
               z_atom_cell(i)
       end do
    end if
    posnStore = zero
    forceStore = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3 * ni_in_cell + 3
    alpha = one
    if (myid == 0 .and. iprint_MD > 0) &
         write(io_lun, fmt='(4x,"SQNM structural relaxation. Maximum of ",i4, " steps with tolerance of ",&
         &f8.4,a2,"/",a2)') MDn_steps, MDcgtol, en_units(energy_units), d_units(dist_units)
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    min_layer = min_layer - 1
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    call dump_pos_and_matrices
    call get_maxf(max)
    min_layer = min_layer + 1
    iter = 0
    press = target_pressure/HaBohr3ToGPa
    ! Stress tolerance in Ha/Bohr3
    stress_target = cell_stress_tol/HaBohr3ToGPa
    enthalpy0 = enthalpy(energy0, press)
    enthalpy1 = enthalpy0
    dH = zero
    max_stress = zero
    volume = rcellx*rcelly*rcellz
    do i=1,3
       stress_diff = abs(press*volume + stress(i,i))/volume
       if (stress_diff > max_stress) max_stress = stress_diff
    end do
    if (inode==ionode) then
       write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," E: ",e18.10)') & 
            iter, for_conv*max, en_conv*energy0
       write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxStr: ",f12.8," MaxF: ",f12.8," H: ",f16.8," dH: ",f12.8)') &
            iter, max_stress*volume, for_conv*max, en_conv*enthalpy1, en_conv*dH
    end if
    iter_loc = 0
    ggold = zero
    energy1 = energy0
    cg_new(1,1:ni_in_cell) = -tot_force(1,:)*rcellx/orcellx
    cg_new(2,1:ni_in_cell) = -tot_force(2,:)*rcellx/orcellx
    cg_new(3,1:ni_in_cell) = -tot_force(3,:)*rcellx/orcellx
    cg_new(1,ni_in_cell+1) = (stress(1,1) - press*volume)*orcell(1)/wscal
    cg_new(2,ni_in_cell+1) = (stress(2,2) - press*volume)*orcell(2)/wscal
    cg_new(3,ni_in_cell+1) = (stress(3,3) - press*volume)*orcell(3)/wscal
    if (inode==ionode .and. iprint_MD > 1) then
       ! This needs scaling
       g0 = dot(length,cg_new,1,cg_new,1)
       write(io_lun,'(4x,"Search direction has magnitude ",f19.8)') sqrt(g0/ni_in_cell)
    end if
    do while (.not. done)
       ! Book-keeping
       iter = iter + 1
       iter_loc = iter_loc + 1
       npmod = mod(iter_loc, LBFGS_history)
       if(npmod==0) npmod = LBFGS_history
       ! Ions
       do i=1,ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i)*orcellx/rcellx
          posnStore (2,jj,npmod) = y_atom_cell(i)*orcelly/rcelly
          posnStore (3,jj,npmod) = z_atom_cell(i)*orcellz/rcellz
          forceStore(1,jj,npmod) = -tot_force(1,jj)*rcellx/orcellx
          forceStore(2,jj,npmod) = -tot_force(2,jj)*rcelly/orcelly
          forceStore(3,jj,npmod) = -tot_force(3,jj)*rcellz/orcellz
          x_new_pos(i) = x_atom_cell(i)
          y_new_pos(i) = y_atom_cell(i)
          z_new_pos(i) = z_atom_cell(i)
          cg(1,i) = -cg_new(1,jj)*orcellx/rcellx ! Search downhill
          cg(2,i) = -cg_new(2,jj)*orcelly/rcelly ! Search downhill
          cg(3,i) = -cg_new(3,jj)*orcellz/rcellz ! Search downhill
       end do
       ! Cell
       posnStore(1,ni_in_cell+1,npmod) = wscal*rcellx/orcellx
       posnStore(2,ni_in_cell+1,npmod) = wscal*rcelly/orcelly
       posnStore(3,ni_in_cell+1,npmod) = wscal*rcellz/orcellz
       forceStore(1,ni_in_cell+1,npmod) = (-stress(1,1) + press*volume)*orcellx/wscal
       forceStore(2,ni_in_cell+1,npmod) = (-stress(2,2) + press*volume)*orcelly/wscal
       forceStore(3,ni_in_cell+1,npmod) = (-stress(3,3) + press*volume)*orcellz/wscal
       cg(1,ni_in_cell+1) = -cg_new(1,ni_in_cell+1)*wscal/orcellx
       cg(2,ni_in_cell+1) = -cg_new(2,ni_in_cell+1)*wscal/orcelly
       cg(3,ni_in_cell+1) = -cg_new(3,ni_in_cell+1)*wscal/orcellz
       ! Set up limits for sums
       if (myid == 0 .and. iprint_MD > 2) &
            write(io_lun,fmt='(2x,"SQNM iteration ",i4)') iter
       ! Line search
       if(iter==1) then
          !call backtrack_linemin_full(cg, energy0, energy1, fixed_potential, vary_mu)
       else
          call single_step_full(cg, energy0, energy1, fixed_potential, vary_mu)
       end if
       if(energy1>energy0) then
          if(inode==ionode.AND.iprint_MD>1) write(io_lun,fmt='(4x,"Energy rise: resetting history")')
          cg_new(:,1:ni_in_cell) = -tot_force
          do i=1,ni_in_cell
             jj = id_glob(i)
             cg(1,i) = -cg_new(1,jj)*orcellx/rcellx ! Search downhill
             cg(2,i) = -cg_new(2,jj)*orcelly/rcelly ! Search downhill
             cg(3,i) = -cg_new(3,jj)*orcellz/rcellz ! Search downhill
          end do
          do i=1,3
             cg_new(i,ni_in_cell+1) = (stress(i,i) - press*volume)*orcell(i)/wscal
             cg(i,ni_in_cell+1) = -cg_new(i,ni_in_cell+1)*wscal/orcell(i)
          end do
          !call backtrack_linemin_full(cg, energy0, energy1, fixed_potential, vary_mu)
          !call single_step(cg, energy0, energy1, fixed_potential, vary_mu)
          if(energy1>energy0) call cq_abort("Energy rise twice in succession: check SCF and other tolerances")
          npmod = 1
          iter_loc = 1
          ! In the original paper, this is alpha = alpha/2 but heuristically this seems better
          ! DRB 2021/09/15
          alpha = one
       endif
       ! Update stored position difference and force difference
       do i=1,ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i)*orcellx/rcellx - posnStore (1,jj,npmod)
          if(abs(posnStore(1,jj,npmod)/r_super_x)>0.7_double) posnStore(1,jj,npmod) &
               = posnStore(1,jj,npmod) &
               - nint(posnStore(1,jj,npmod)/r_super_x)*r_super_x
          posnStore (2,jj,npmod) = y_atom_cell(i)*orcelly/rcelly - posnStore (2,jj,npmod)
          if(abs(posnStore(2,jj,npmod)/r_super_y)>0.7_double) posnStore(2,jj,npmod) &
               = posnStore(2,jj,npmod) &
               - nint(posnStore(2,jj,npmod)/r_super_y)*r_super_y
          posnStore (3,jj,npmod) = z_atom_cell(i)*orcellz/rcellz - posnStore (3,jj,npmod)
          if(abs(posnStore(3,jj,npmod)/r_super_z)>0.7_double) posnStore(3,jj,npmod) &
               = posnStore(3,jj,npmod) &
               - nint(posnStore(3,jj,npmod)/r_super_z)*r_super_z
          forceStore(1,jj,npmod) = -tot_force(1,jj)*rcellx/orcellx - forceStore(1,jj,npmod)
          forceStore(2,jj,npmod) = -tot_force(2,jj)*rcelly/orcelly - forceStore(2,jj,npmod)
          forceStore(3,jj,npmod) = -tot_force(3,jj)*rcellz/orcellz - forceStore(3,jj,npmod)
          ! New search direction
       end do
       posnStore (1,ni_in_cell+1,npmod) = rcellx*wscal/orcellx - posnStore (1,ni_in_cell+1,npmod)
       posnStore (2,ni_in_cell+1,npmod) = rcelly*wscal/orcelly - posnStore (2,ni_in_cell+1,npmod)
       posnStore (3,ni_in_cell+1,npmod) = rcellz*wscal/orcellz - posnStore (3,ni_in_cell+1,npmod)
       forceStore(1,ni_in_cell+1,npmod) = (-stress(1,1) + press*volume)*orcellx/wscal - forceStore(1,ni_in_cell+1,npmod)
       forceStore(2,ni_in_cell+1,npmod) = (-stress(2,2) + press*volume)*orcelly/wscal - forceStore(2,ni_in_cell+1,npmod)
       forceStore(3,ni_in_cell+1,npmod) = (-stress(3,3) + press*volume)*orcellz/wscal - forceStore(3,ni_in_cell+1,npmod)
       n_store = min(iter_loc,LBFGS_history) ! Number of stored states
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,fmt='(4x,"Number of stored histories ",i3)') n_store
       allocate(mod_dr(n_store),Sij(n_store,n_store),lambda(n_store),&
            omega(n_store,n_store))
       mod_dr = zero
       ! Normalise dR and dg
       do i=1,n_store
          mod_dr(i) = sqrt(dot(length,posnStore(:,:,i),1,posnStore(:,:,i),1))
       end do
       posnStore(:,:,npmod) = posnStore(:,:,npmod)/mod_dr(npmod)
       forceStore(:,:,npmod) = forceStore(:,:,npmod)/mod_dr(npmod)
       !cg_new = -tot_force ! The L-BFGS is in terms of grad E
       ! Add call to write_atomic_positions and write_xsf (2020/01/17: smujahed)
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)
       ! Build significant subspace
       Sij = zero
       omega = zero
       do i = 1, n_store
          do j=i,n_store
             Sij(j,i) = dot(length,posnStore(:,:,j),1,posnStore(:,:,i),1)
             if(j>i) Sij(i,j) = Sij(j,i)
          end do
       end do
       ! Solve for eigenvectors of Sij
       omega = Sij
       if(n_store>1) then
          call syev('U',n_store,omega,n_store,lambda,info)
          if(info<0) call cq_abort("Error in SQNM calling dsyev: ",info)
          if(info>0.and.inode==ionode) write(io_lun,fmt='(4x,"Possible error in SQNM; dsyev returned ",i4)') info
       else
          lambda = one
          omega = one
       end if
       if(inode==ionode.AND.iprint_MD>2) then
          write(io_lun,fmt='(4x,"Eigenvalues of Sij: ",(f7.4))') lambda
       end if
       lambda_max = maxval(lambda)
       n_dim = n_store
       do i=1, n_store
          if(lambda(i)/lambda_max<1.0e-4_double) then
             n_dim = n_dim - 1
          end if
       end do
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,fmt='(4x,"Number of eigenstates kept: ",i3)') n_dim
       ! Build dr_tilde and dg_tilde
       allocate(dr_tilde(3,ni_in_cell+1,n_dim), dg_tilde(3,ni_in_cell+1,n_dim))
       dr_tilde = zero
       dg_tilde = zero
       do i=1,n_dim
          do j=1,n_store
             dr_tilde(:,:,i) = dr_tilde(:,:,i) + omega(j,i+(n_store-n_dim))*posnStore(:,:,j)
             dg_tilde(:,:,i) = dg_tilde(:,:,i) + omega(j,i+(n_store-n_dim))*forceStore(:,:,j)
          end do
          dr_tilde(:,:,i) = dr_tilde(:,:,i)/sqrt(lambda(i+(n_store-n_dim)))
          dg_tilde(:,:,i) = dg_tilde(:,:,i)/sqrt(lambda(i+(n_store-n_dim)))
       end do
       ! Construct approximate Hessian projected onto subspace
       allocate(Hij(n_dim,n_dim))
       Hij = zero
       do i=1,n_dim
          do j=i,n_dim
             Hij(j,i) = half*(dot(length,dr_tilde(:,:,i),1,dg_tilde(:,:,j),1) + &
                  dot(length,dr_tilde(:,:,j),1,dg_tilde(:,:,i),1))
             if(j>i) Hij(i,j) = Hij(j,i)
          end do
       end do
       ! Find eigenvectors
       allocate(kappa(n_dim),vi(n_dim,n_dim),kappa_prime(n_dim))
       kappa = zero
       kappa_prime = zero
       vi = Hij
       if(n_dim>1) then
          call syev('U',n_dim,vi,n_dim,kappa,info)
          if(info<0) call cq_abort("Error in SQNM calling dsyev: ",info)
          if(info>0.and.inode==ionode) write(io_lun,fmt='(4x,"Possible error in SQNM; dsyev returned ",i4)') info
       else
          kappa = one !Hij(1,1)
          vi = one
       end if
       if(inode==ionode.AND.iprint_MD>3) write(io_lun,fmt='(4x,"Kappa: ",(f7.4))') kappa
       ! Build v tilde
       allocate(vi_tilde(3,ni_in_cell+1,n_dim),ri_vec(3,ni_in_cell+1))
       vi_tilde = zero
       do i=1,n_dim
          ri_vec = zero
          do j=1,n_dim
             vi_tilde(:,:,i) = vi_tilde(:,:,i) + vi(j,i)*dr_tilde(:,:,j)
             ri_vec(:,:) = ri_vec(:,:) + vi(j,i)*dg_tilde(:,:,j)
          end do
          ri_vec(:,:) = ri_vec(:,:) - kappa(i)*vi_tilde(:,:,i)
          !kappa_prime(i) = sqrt(0.0025_double + kappa(i)*kappa(i))
          kappa_prime(i) = sqrt(dot(length,ri_vec,1,ri_vec,1) + kappa(i)*kappa(i))
       end do
       if(inode==ionode.AND.iprint_MD>3) write(io_lun,fmt='(4x,"Kappa prime: ",(f7.4))') kappa_prime
       ! Build preconditioned search
       cg_new(1,1:ni_in_cell) = -alpha*tot_force(1,:)*rcellx/orcellx
       cg_new(2,1:ni_in_cell) = -alpha*tot_force(2,:)*rcellx/orcellx
       cg_new(3,1:ni_in_cell) = -alpha*tot_force(3,:)*rcellx/orcellx
       cg_new(1,ni_in_cell+1) =  alpha*(stress(1,1) - press*volume)*orcell(1)/wscal
       cg_new(2,ni_in_cell+1) =  alpha*(stress(2,2) - press*volume)*orcell(2)/wscal
       cg_new(3,ni_in_cell+1) =  alpha*(stress(3,3) - press*volume)*orcell(3)/wscal
       rcell(1) = rcellx
       rcell(2) = rcelly
       rcell(3) = rcellz
       do i=1,n_dim
          temp = zero
          do j=1,3
             temp = temp + dot(ni_in_cell,tot_force(j,:),1,vi_tilde(j,:,i),1)*rcell(j)/orcell(j)
             temp = temp + (stress(j,j)-press*volume)*orcell(j)*vi_tilde(j,ni_in_cell+1,i)/wscal
          end do
          do j=1,3
             cg_new(j,1:ni_in_cell) = cg_new(j,1:ni_in_cell) - (one/kappa_prime(i) - alpha)* &
                  temp*vi_tilde(j,1:ni_in_cell,i)
             cg_new(j,ni_in_cell+1) = cg_new(j,ni_in_cell+1) - (one/kappa_prime(i) - alpha)* &
                  temp*vi_tilde(j,ni_in_cell+1,i)
          end do
       end do
       ! Zero search direction for fixed atoms and find maximum force
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
             if (.not. flag_move_atom(k,i)) then
                cg_new(k,i) = zero
             end if
          end do
       end do
       gg = dot(length, cg_new, 1, cg_new, 1)
       ! Analyse forces
       g0 = zero
       g0 = g0 + dot(ni_in_cell,tot_force(1,:),1,tot_force(1,:),1)*rcellx*rcellx/(orcellx*orcellx)
       g0 = g0 + dot(ni_in_cell,tot_force(2,:),1,tot_force(2,:),1)*rcelly*rcelly/(orcelly*orcelly)
       g0 = g0 + dot(ni_in_cell,tot_force(3,:),1,tot_force(3,:),1)*rcellz*rcellz/(orcellz*orcellz)
       g0 = g0 + (stress(1,1)-press*volume)*(stress(1,1)-press*volume)*orcellx*orcellx/(wscal*wscal)
       g0 = g0 + (stress(2,2)-press*volume)*(stress(2,2)-press*volume)*orcelly*orcelly/(wscal*wscal)
       g0 = g0 + (stress(3,3)-press*volume)*(stress(3,3)-press*volume)*orcellz*orcellz/(wscal*wscal)
       f_dot_sd = zero
       f_dot_sd = f_dot_sd - dot(ni_in_cell,cg_new(1,1:ni_in_cell),1,tot_force(1,:),1)*rcellx/orcellx
       f_dot_sd = f_dot_sd - dot(ni_in_cell,cg_new(2,1:ni_in_cell),1,tot_force(2,:),1)*rcelly/orcelly
       f_dot_sd = f_dot_sd - dot(ni_in_cell,cg_new(3,1:ni_in_cell),1,tot_force(3,:),1)*rcellz/orcellz
       f_dot_sd = f_dot_sd + cg_new(1,ni_in_cell+1)*(stress(1,1)-press*volume)*orcellx/wscal
       f_dot_sd = f_dot_sd + cg_new(2,ni_in_cell+1)*(stress(2,2)-press*volume)*orcelly/wscal
       f_dot_sd = f_dot_sd + cg_new(3,ni_in_cell+1)*(stress(3,3)-press*volume)*orcellz/wscal
       f_dot_sd = f_dot_sd/sqrt(g0*gg)
       !f_dot_sd = dot(length,cg_new,1,-tot_force,1)/sqrt(g0*gg)
       if(inode==ionode.AND.iprint_MD>2) &
            write(io_lun,fmt='(4x,"Dot product of search direction and force: ",f12.7)') f_dot_sd
       ! Now adjust alpha
       if(f_dot_sd>0.2_double) then
          alpha = alpha*1.1_double
       else
          alpha = alpha*0.85_double
       end if
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,fmt='(4x,"Alpha set to: ",f9.5)') alpha
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
          end do
       end do
       dE = energy1 - energy0
       energy0 = energy1
       if (inode==ionode) then
          write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," E: ",e18.10," dE: ",f12.8)') & 
               iter, for_conv*max, en_conv*energy1, en_conv*dE
          if (iprint_MD > 1) then
             g0 = dot(length, tot_force, 1, tot_force, 1)
             write(io_lun,'(4x,"Force Residual:     ",f19.8," ",a2,"/",a2)') &
                  for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
                  d_units(dist_units)
             write(io_lun,'(4x,"Maximum force:      ",f19.8)') for_conv*max
             write(io_lun,'(4x,"Force tolerance:    ",f19.8)') for_conv*MDcgtol
             write(io_lun,'(4x,"Energy change:      ",f19.8," ",a2)') &
                  en_conv*dE, en_units(energy_units)
             g0 = dot(length,cg_new,1,cg_new,1)
             write(io_lun,'(4x,"Search direction has magnitude ",f19.8)') sqrt(g0/ni_in_cell)
          end if
       end if
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write(io_lun, fmt='(4x,"Exceeded number of MD steps: ",i4)') iter
       end if
       if (abs(max) < MDcgtol) then
          done = .true.
          if (inode==ionode) then
             write(io_lun, fmt='(4x,"Maximum force below threshold: ",f12.5)') max
             write(io_lun,'(2x,a,i4,a)') "GeomOpt converged in ", iter, " iterations"
          end if
       end if
       deallocate(dr_tilde,dg_tilde,Hij,kappa,vi,kappa_prime,vi_tilde,ri_vec)
       deallocate(mod_dr,Sij,lambda,omega)
       if (.not. done) call check_stop(done, iter)
       if(done) exit
    end do ! .not. done i.e. until max iterations or force tolerance reached
    deallocate(cg, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating cg in control: ", &
         ni_in_cell, stat)
    call reg_dealloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    deallocate(posnStore, forceStore, cg_new, x_new_pos, y_new_pos, z_new_pos)
  end subroutine full_sqnm
  !!***

  !!****f* control/cell_cg_run *
  !!
  !!  NAME
  !!   cg_run - Does CG minimisation of simulation box dimensions
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs CG minimisation of the total energy (WRT a, b and c) by
  !!   repeated calls to the minimizer
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   J.S.Baker
  !!   S.Mujahed
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   01/06/17 J.S.Baker
  !!  MODIFICATION HISTORY
  !!   2021/10/13 08:44 dave
  !!    Changes to account for stress tolerance from user being in GPa
  !!   2021/10/15 17:37 dave
  !!    Tweak output and iteration update
  !!   2022/09/16 16:57 dave
  !!    Added backtrack line minimiser
  !!   2022/10/05 08:50 dave
  !!    Improving output
  !!  SOURCE
  !!
  subroutine cell_cg_run(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, only: iprint_gen, ni_in_cell, x_atom_cell,  &
         y_atom_cell, z_atom_cell, id_glob,    &
         atom_coord, rcellx, rcelly, rcellz,   &
         area_general, iprint_MD,              &
         IPRINT_TIME_THRES1, cell_en_tol,      &
         cell_constraint_flag, cell_stress_tol, min_layer
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin_cell, enthalpy, enthalpy_tolerance, &
         backtrack_linemin_cell, adapt_backtrack, backtrack, safe, cg_line_min
    use GenComms,      only: gsum, myid, inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: stress, tot_force
    use io_module,     only: write_atomic_positions, pdb_template, &
         check_stop, print_atomic_positions, return_prefix
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use io_module,      only: leqi
    use dimens, ONLY: r_super_x, r_super_y, r_super_z
    use store_matrix,  only: dump_pos_and_matrices
    use md_control,    only: target_pressure

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma, &
         enthalpy0, enthalpy1, dH, press
    integer        :: i,j,k,iter,length, jj, lun, stat, reset_iter
    logical        :: done
    type(cq_timer) :: tmr_l_iter
    real(double) :: new_rcellx, new_rcelly, new_rcellz, search_dir_x, search_dir_y,&
         search_dir_z, stressx, stressy, stressz, RMSstress, newRMSstress,&
         dRMSstress, search_dir_mean, mean_stress, max_stress, &
         stress_diff, volume, stress_target
    real(double), dimension(3) :: cg

    character(len=20) :: subname = "cell_cg_run: "
    character(len=120) :: prefix, prefixGO
    character(len=10) :: prefixF = "          "

    prefix = return_prefix(subname, min_layer)
    prefixGO = return_prefix("GeomOpt", min_layer)
    if (myid == 0) then
       if(cg_line_min==safe) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with safemin line minimisation"
       else if (cg_line_min==backtrack) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with backtracking line minimisation"
       end if
       if(iprint_gen + min_layer>0) write(io_lun, fmt='(4x,a,f8.4,a,i4)') &
            trim(prefix)//" tolerance: ",cell_stress_tol," GPa Maximum steps: ",MDn_steps
    end if
    search_dir_x = zero
    search_dir_y = zero
    search_dir_z = zero
    search_dir_mean = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    min_layer = min_layer - 1
    if (iprint_MD + min_layer > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.,0)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    end if
    min_layer = min_layer + 1
    iter = 1
    reset_iter = 1
    ggold = zero
    energy1 = energy0

    press = target_pressure/HaBohr3ToGPa
    ! Stress tolerance in Ha/Bohr3
    stress_target = cell_stress_tol/HaBohr3ToGPa
    enthalpy0 = enthalpy(energy0, press)
    enthalpy1 = enthalpy0
    dH = zero
    max_stress = zero
    volume = rcellx*rcelly*rcellz
    do i=1,3
       stress_diff = abs(press*volume + stress(i,i))/volume
       if (stress_diff > max_stress) max_stress = stress_diff
    end do
    if (inode==ionode) then
       write(io_lun,'(/4x,a,i4," MaxStr: ",f12.8," GPa H: ",f16.8,x,a2," dH: ",f12.8,a2/)') &
            trim(prefixGO)//" - Iter: ",0, max_stress*HaBohr3ToGPa, en_conv*enthalpy0, &
            en_units(energy_units), zero, en_units(energy_units)
    end if
    ! Check for trivial case where pressure is converged
    if(max_stress < stress_target) then
       if (inode==ionode) &
            write(io_lun,'(4x,a,i4,a)') trim(prefixGO)//" converged in ", &
            iter, " iterations"
       done = .true.
       if (myid == 0 .and. iprint_gen > 0) &
            write(io_lun, fmt='(4x,a,f19.8," GPa")') &
            trim(prefix)//" maximum stress below threshold:   ",max_stress*HaBohr3ToGPa
       return
    end if
    call dump_pos_and_matrices(index=0,MDstep=iter)
    do while (.not. done)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       volume = rcellx*rcelly*rcellz
       ! Keep these as stresses
       stressx = -stress(1,1)!/volume
       stressy = -stress(2,2)!/volume
       stressz = -stress(3,3)!/volume
       mean_stress = (stressx + stressy + stressz)/3
       RMSstress = sqrt(((stressx*stressx) + (stressy*stressy) + (stressz*stressz))/3)

       ! Construct ratio for conjugacy. Constraints are initially applied within
       ! get_gamma_cell_cg.
       call get_gamma_cell_cg(ggold, gg, gamma, stressx, stressy, stressz)
       if (CGreset) then
          if (gamma > one .or. reset_iter > length) then
             if (inode == ionode .and. iprint_MD + min_layer > 1) &
                  write(io_lun,fmt='(4x,a)') &
                  trim(prefix)//" gamma>1, so resetting direction"
             gamma = zero
             reset_iter = 0
          end if
       end if
       if (inode == ionode .and. iprint_MD + min_layer > 2) then
          write(io_lun,fmt='(4x,a,f12.8)') &
               trim(prefix)//' gamma = ', gamma
       end if
       ggold = gg

       !Build search direction.
       ! If the volume constraint is set, there is only one search direction!
       ! This is the direction which minimises the mean stress.
       if (leqi(cell_constraint_flag, 'volume')) then
          search_dir_mean = gamma*search_dir_mean + mean_stress - press*volume
       else
          search_dir_x = gamma*search_dir_x + stressx - press*volume
          search_dir_y = gamma*search_dir_y + stressy - press*volume
          search_dir_z = gamma*search_dir_z + stressz - press*volume
       end if

       new_rcellx = rcellx
       new_rcelly = rcelly
       new_rcellz = rcellz

       ! Minimise in this direction. Constraint information is also used within
       ! safemin_cell. Look in move_atoms.module.f90 for further information.
       min_layer = min_layer - 1
       if(cg_line_min==safe) then
          call safemin_cell(new_rcellx, new_rcelly, new_rcellz, search_dir_x, &
               search_dir_y, search_dir_z, search_dir_mean, press, &
               enthalpy0, enthalpy1, fixed_potential, vary_mu)
       else if(cg_line_min==backtrack.OR.cg_line_min==adapt_backtrack) then
          cg(1) = search_dir_x
          cg(2) = search_dir_y
          cg(3) = search_dir_z
          call backtrack_linemin_cell(cg, press, enthalpy0, enthalpy1, fixed_potential, vary_mu)
       end if
       min_layer = min_layer + 1
       ! Output positions to UpdatedAtoms.dat
       if (myid == 0 .and. iprint_gen + min_layer > 1) then
          write(io_lun, fmt='(/4x,a)') trim(prefix)//" simulation cell dimensions: "
          write(io_lun, fmt='(6x,f12.5,1x,a2," x ",f12.5,1x,a2," x ",f12.5,1x,a2)') &
            rcellx, d_units(dist_units), rcelly, d_units(dist_units), rcellz, d_units(dist_units)
       end if
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))

       ! Analyse Stresses and energies
       dH = enthalpy1 - enthalpy0
       volume = rcellx*rcelly*rcellz
       newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
            (stress(2,2)*stress(2,2)) + &
            (stress(3,3)*stress(3,3)))/3)
       dRMSstress = (RMSstress - newRMSstress)/volume

       enthalpy0 = enthalpy1
       total_energy = enthalpy1
       ! Check exit criteria
       max_stress = zero
       do i=1,3
          stress_diff = abs(press*volume + stress(i,i))/volume
          if (stress_diff > max_stress) max_stress = stress_diff
       end do

       reset_iter = reset_iter +1

       if (inode==ionode) then
          write(io_lun,'(/4x,a,i4," MaxStr: ",f12.8," GPa H: ",f16.8,x,a2," dH: ",f12.8,a2/)') &
               trim(prefixGO)//" - Iter: ",iter, max_stress*HaBohr3ToGPa, en_conv*enthalpy1, &
               en_units(energy_units), en_conv*dH, en_units(energy_units)
          if (iprint_MD + min_layer > 1) then
             write(io_lun,'(4x,a,f21.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Maximum stress:     ", max_stress*HaBohr3ToGPa
             write(io_lun,'(4x,a,f21.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Stress tolerance:   ", cell_stress_tol
             write(io_lun,'(4x,a,f21.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Change in stress:   ", dRMSstress*HaBohr3ToGPa
             write(io_lun,'(4x,a,f21.8," ",a2)') &
                  prefixF(1:-2*min_layer)//"Enthalpy change:    ", en_conv*dH, en_units(energy_units)
             write(io_lun,'(4x,a,f21.8," ",a2)') &
                  prefixF(1:-2*min_layer)//"Enthalpy tolerance: ", en_conv*enthalpy_tolerance, &
                  en_units(energy_units)
          else if (iprint_MD + min_layer > -1) then
             write(io_lun,'(4x,a,f21.8," ",a2)') &
                  prefixF(1:-2*min_layer)//"Enthalpy change:    ", en_conv*dH, en_units(energy_units)
             write(io_lun,'(4x,a,f21.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Maximum stress:     ", max_stress*HaBohr3ToGPa
          end if
       end if

       ! First exit is if too many steps have been taken. Default is 50.
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write(io_lun, fmt='(4x,a,i4)') trim(prefix)//" exceeded number of MD steps: ",iter
       end if

       ! Second exit is if the desired enthalpy and stress tolerancex have ben reached
       if (abs(dH)<enthalpy_tolerance .and. max_stress < stress_target) then
          done = .true.
          if (inode==ionode) then
             write(io_lun,'(/4x,a,i4,a)') trim(prefixGO)//" converged in ", &
                  iter, " iterations"
             if (iprint_MD + min_layer > 0) &
                  write(io_lun, fmt='(4x,a,f19.8," ",a2)') &
                  trim(prefix)//"Enthalpy change below threshold: ", dH*en_conv, en_units(energy_units)
             write(io_lun, fmt='(4x,a,f19.8," GPa")') &
                  trim(prefix)//"Maximum stress below threshold:  ", max_stress*HaBohr3ToGPa
          end if
       end if

       call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
       if (.not. done) call check_stop(done, iter)
       iter = iter + 1
    end do

    if (myid == 0 .and. iprint_gen + min_layer > 0) then
       write(io_lun, fmt='(/4x,a)') trim(prefix)//" final simulation box dimensions are: "
       write(io_lun, fmt='(8x,f12.5,1x,a2," x ",f12.5,1x,a2," x ",f12.5,1x,a2)') &
            rcellx, d_units(dist_units), rcelly, d_units(dist_units), rcellz, d_units(dist_units)
    end if

    call reg_dealloc_mem(area_general, 6*ni_in_cell, type_dbl)

  end subroutine cell_cg_run

  !!***

  !!****f* control/get_gamma_cell_cg *
  !!
  !!  NAME
  !!   get_gamma_cell_cg
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Gets the Gamma value for the CG algorithm in cell_cg_run based on the
  !!   user set constraints.
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   J.S.Baker
  !!
  !!  CREATION DATE
  !!  01/06/17 J.S.Baker
  !!  MODIFICATION HISTORY
  !!   2017/08/29 dave
  !!    Indentation fixed (emacs standard) and added check for inode to write statement
  !!  SOURCE
  !!
  subroutine get_gamma_cell_cg(ggold, gg, gamma, stressx, stressy, stressz)

    !Module usage
    use io_module,      only: leqi
    use global_module, only: cell_constraint_flag, iprint_gen
    use numbers
    use GenComms,      only: gsum, myid, inode, ionode

    implicit none

    !Passed arguments
    real(double) :: ggold, gg, stressx, stressy, stressz, gamma, mean_stress

    ! gamma is zero if the first iteration
    if (abs(ggold) < 1.0e-6_double) then
       gamma = zero
    else
       ! use mean stress to scale by volume
       if (leqi(cell_constraint_flag, 'volume')) then
          mean_stress = (stressx + stressy + stressz)/3
          gg = mean_stress*mean_stress

          ! no contraints or single ratio fixed = all three stress components used
       else if (leqi(cell_constraint_flag, 'none') .or. leqi(cell_constraint_flag, 'c/a') &
            .or. leqi(cell_constraint_flag, 'a/c') .or. leqi(cell_constraint_flag, 'a/b') & 
            .or. leqi(cell_constraint_flag, 'b/a') .or. leqi(cell_constraint_flag, 'b/c') &
            .or. leqi(cell_constraint_flag, 'c/b')) then

          gg = stressx*stressx + stressy*stressy + stressz*stressz

          ! Fix a single lattice parameter (no longer need a stress component)
       else if (leqi(cell_constraint_flag, 'a')) then
          gg = stressy*stressy + stressz*stressz
       else if (leqi(cell_constraint_flag, 'b')) then
          gg = stressx*stressx + stressz*stressz
       else if (leqi(cell_constraint_flag, 'c')) then
          gg = stressy*stressy + stressx*stressx

          ! Fix more than one lattice parameter (only one stress component required)
       else if (leqi(cell_constraint_flag, 'c a') .or. leqi(cell_constraint_flag, 'a c')) then
          gg = stressy*stressy
       else if (leqi(cell_constraint_flag, 'a b') .or. leqi(cell_constraint_flag, 'b a')) then
          gg = stressz*stressz
       else if (leqi(cell_constraint_flag, 'b c') .or. leqi(cell_constraint_flag, 'c b')) then
          gg = stressx*stressx

       end if

       gamma = gg/ggold

    end if
    return
  end subroutine get_gamma_cell_cg
!!***

  !!****f* control/full_double_loop *
  !!
  !!  NAME
  !!   full_double_loop
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Simple implementation of cell + ionic geometry optimisation using a full ionic relaxation followed
  !!   by a full cell optimisation, and repeat until converged.  It's not
  !!   clear if this is more efficient than performing full ionic relaxation
  !!   after each line minimisation of the cell relaxation but does seem so
  !!   in tests.
  !!
  !!   Now allows choice of CG vs SQNM
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Z Raza
  !!  CREATION DATE
  !!   2018/02/06
  !!  MODIFICATION HISTORY
  !!   2021/10/15 10:45 dave
  !!    Rewrite to call cg_run for ions and cell_cg_run for cell sequentially
  !!   2022/08/15 10:31 dave
  !!    Add possibility of SQNM instead of CG
  !!  SOURCE
  !!
  subroutine full_double_loop(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, only: iprint_gen, ni_in_cell, x_atom_cell,  &
         y_atom_cell, z_atom_cell, id_glob,    &
         atom_coord, area_general, iprint_MD,  &
         IPRINT_TIME_THRES1,                   &
         cell_en_tol, cell_stress_tol,         &
         rcellx, rcelly, rcellz, runtype, min_layer
    use input_module,         only: leqi
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: enthalpy, enthalpy_tolerance
    use GenComms,      only: inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: tot_force, stress
    use io_module,     only: write_atomic_positions, pdb_template, &
         check_stop, write_xsf, return_prefix
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  only: dump_InfoMatGlobal, dump_pos_and_matrices
    use md_control,    only: flag_write_xsf, target_pressure

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables for ionic relaxation
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma,gg1, &
         enthalpy0, enthalpy1, dH, press
    integer        :: i,j,k,iter,length, jj, lun, stat
    logical        :: done_ions, done_cell
    type(cq_timer) :: tmr_l_iter
    real(double), allocatable, dimension(:,:) :: cg, old_force
    real(double), allocatable, dimension(:)   :: x_new_pos, y_new_pos,&
         z_new_pos

    ! Local variables for cell relaxation
    real(double) :: new_rcellx, new_rcelly, new_rcellz, search_dir_x, &
         search_dir_y, search_dir_z, stressx, stressy, stressz, &
         RMSstress, newRMSstress, dRMSstress, search_dir_mean, &
         mean_stress, energy2, volume, stress_diff, max_stress, stress_target
    integer      :: iter_cell
    character(len=20) :: subname = "full_double_loop: "
    character(len=120) :: prefix, prefixGO
    character(len=10) :: prefixF = "          "

    prefix = return_prefix(subname, min_layer)
    prefixGO = return_prefix("GeomOpt", min_layer)
    allocate(cg(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ",&
         ni_in_cell, stat)
    allocate(old_force(3,ni_in_cell), STAT=stat)
    allocate(x_new_pos(ni_in_cell), y_new_pos(ni_in_cell), &
         z_new_pos(ni_in_cell), STAT=stat)
    if (stat/=0) &
         call cq_abort("Error allocating _new_pos in control: ", &
         ni_in_cell,stat)
    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    search_dir_z = zero
    search_dir_x = zero
    search_dir_y = zero
    search_dir_mean = zero
    cg = zero
    done_cell = .false.
    length = 3 * ni_in_cell
    if (inode == ionode) then
       if ( leqi(runtype, 'cg')    ) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with CG for ions and cell"
       else if ( leqi(runtype, 'sqnm')    ) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with SQNM for ions and CG for cell"
       end if
       if(iprint_gen + min_layer>0) then
          write(io_lun, fmt='(4x,a,f8.4," ",a2,"/",a2,a,i4)') &
               trim(prefix)//" force tolerance:  ",MDcgtol,en_units(energy_units), d_units(dist_units),&
               " Maximum steps: ",MDn_steps
          write(io_lun,'(4x,a,f8.4," GPa")') trim(prefix)//" stress tolerance: ",cell_stress_tol
       end if
    end if
    energy0 = total_energy
    energy1 = zero
    energy2 = zero
    dE = zero

    ! Find energy and forces
    min_layer = min_layer - 1
    if (iprint_MD + min_layer > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.,0)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    end if
    call dump_pos_and_matrices
    call get_maxf(max)
    min_layer = min_layer + 1
    press = target_pressure/HaBohr3ToGPa
    ! Stress tolerance in Ha/Bohr3
    stress_target = cell_stress_tol/HaBohr3ToGPa
    enthalpy0 = enthalpy(energy0, press)
    dH = zero
    volume = rcellx*rcelly*rcellz
    max_stress = zero
    do i=1,3
       stress_diff = abs(press*volume + stress(i,i))/volume
       if (stress_diff > max_stress) max_stress = stress_diff
    end do
    if (inode==ionode) then
       write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," H: " ,f16.8,x,a2," MaxS: ",f12.8," GPa")') &
            trim(prefixGO)//" + Iter: ",0, for_conv*max, en_units(energy_units), d_units(dist_units), &
            enthalpy0, en_units(energy_units), max_stress*HaBohr3ToGPa
    end if
    ! Check for trivial case where pressure is converged
    if(max_stress < stress_target .and. abs(max) < MDcgtol) then
       if (inode==ionode) &
            write(io_lun,'(4x,a,i4,a)') trim(prefixGO)//" converged in ", &
            iter, " iterations"
       if (inode==ionode .and. iprint_gen > 0) then
          write(io_lun, fmt='(4x,a,f19.8," GPa")') &
               trim(prefix)//" maximum stress below threshold:   ",max_stress*HaBohr3ToGPa
          write(io_lun, fmt='(4x,a,f19.8,x,a2,"/",a2)') &
               trim(prefix)//" maximum force below threshold:    ",for_conv*max_stress, &
               en_units(energy_units), d_units(dist_units)
       end if
       return
    end if

    iter = 1
    iter_cell = 1
    ! Cell loop
    do while (.not. done_cell)
       ! Relax ions
       if(abs(max)>MDcgtol) then
          if(inode==ionode) write(io_lun,fmt='(4x,a)') trim(prefix)//" ionic relaxation"
          !min_layer = min_layer - 1
          if ( leqi(runtype, 'cg')    ) then
             call cg_run(fixed_potential, vary_mu, energy1)
          else if ( leqi(runtype, 'sqnm')    ) then
             call sqnm(fixed_potential, vary_mu, energy1)
          end if
          !min_layer = min_layer + 1
          ! Analyse forces, stresses and energies
          call get_maxf(max)
       else
          energy1 = energy0
          if(inode==ionode) write(io_lun,fmt='(4x,a)') trim(prefix)//" no ionic relaxation (force converged)"
       end if
       enthalpy1 = enthalpy(energy1, press)
       dH = enthalpy1 - enthalpy0
       volume = rcellx*rcelly*rcellz
       newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
            (stress(2,2)*stress(2,2)) + &
            (stress(3,3)*stress(3,3)))/three)
       dRMSstress = (RMSstress - newRMSstress)/volume

       max_stress = zero
       do i=1,3
          stress_diff = abs(press*volume + stress(i,i))/volume
          if (stress_diff > max_stress) max_stress = stress_diff
       end do
       ! Test for finish
       if (abs(max) < MDcgtol .and. max_stress < stress_target) then
          done_cell = .true.
          if (inode==ionode) then
             write(io_lun,'(/4x,a,i4,a)') trim(prefixGO)//" converged in ", &
                  iter," combined ionic and cell steps"
             write(io_lun, fmt='(4x,a," Maximum stress below threshold:   ",f19.8," GPa")') &
                  trim(prefix),max_stress*HaBohr3ToGPa
             write(io_lun, fmt='(4x,a," Maximum force below threshold:    ",f19.8," ",a2,"/",a2)') &
                  trim(prefix),for_conv*max, en_units(energy_units), d_units(dist_units)
          end if
          exit
       end if
       ! Relax cell
       if(max_stress > stress_target) then
          if(inode==ionode) write(io_lun,fmt='(4x,a)') trim(prefix)//" cell relaxation"
          !min_layer = min_layer - 1
          !if ( leqi(runtype, 'cg')    ) then
          call cell_cg_run(fixed_potential, vary_mu, energy1)
          !else if ( leqi(runtype, 'sqnm')    ) then
          !   call cell_sqnm(fixed_potential, vary_mu, energy1)
          !end if
          !min_layer = min_layer + 1
          ! Analyse forces, stresses and energies
          call get_maxf(max)
       else
          if(inode==ionode) write(io_lun,fmt='(4x,a)') trim(prefix)//" no cell relaxation (stress converged)"
       end if
       enthalpy1 = enthalpy(energy1, press)
       dH = enthalpy1 - enthalpy0
       newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
            (stress(2,2)*stress(2,2)) + &
            (stress(3,3)*stress(3,3)))/three)
       dRMSstress = (RMSstress - newRMSstress)/volume

       volume = rcellx*rcelly*rcellz
       max_stress = zero
       do i=1,3
          stress_diff = abs(press*volume + stress(i,i))/volume
          if (stress_diff > max_stress) max_stress = stress_diff
       end do
       ! Test for finish
       if (inode==ionode) then
          write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," H: " ,f16.8,x,a2," MaxS: ",f12.8," GPa/")') &
               trim(prefixGO)//" + Iter: ",iter, for_conv*max, en_units(energy_units), d_units(dist_units), &
               enthalpy1, en_units(energy_units), max_stress*HaBohr3ToGPa
          if (iprint_MD + min_layer > 1) then
             g0 = dot(length, tot_force, 1, tot_force, 1)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force Residual:     ", &
                  for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Maximum force:      ", &
                  for_conv*max, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force tolerance:    ", &
                  for_conv*MDcgtol, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Maximum stress:     ", max_stress*HaBohr3ToGPa
             write(io_lun,'(4x,a,f19.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Stress tolerance:   ", cell_stress_tol
             write(io_lun,'(4x,a,f19.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Change in stress:   ", dRMSstress*HaBohr3ToGPa
             ! I don't think that these are helpful: they refer to the change
             ! over the complete step, and aren't used in the stopping criteria
             ! 2022/10/11 15:17 dave
             !write(io_lun,'(6x,"Enthalpy change:    ",f19.8," ",a2)') &
             !     en_conv*dH, en_units(energy_units)
             !write(io_lun,'(6x,"Enthalpy tolerance: ",f19.8," ",a2)') &
             !     en_conv*enthalpy_tolerance, en_units(energy_units)
          end if
       end if
       if (abs(max) < MDcgtol .and. max_stress < stress_target) then
          done_cell = .true.
          if (inode==ionode) then
             write(io_lun,'(/4x,a,i4,a)') trim(prefixGO)//" converged in ", iter, " iterations"
             write(io_lun, fmt='(4x,a," Maximum stress below threshold:   ",f19.8," GPa")') &
                  trim(prefix),max_stress*HaBohr3ToGPa
             write(io_lun, fmt='(4x,a," Maximum force below threshold:    ",f19.8," ",a2,"/",a2)') &
                  trim(prefix),for_conv*max, en_units(energy_units), d_units(dist_units)
          end if
       end if

       iter = iter + 1

       energy0 = energy1
       enthalpy0 = enthalpy1

       ! Check exit criteria
       if (iter > MDn_steps) then
          done_cell = .true.
          if (inode==ionode) &
               write(io_lun, fmt='(4x,a,i4)') trim(prefix)//" exceeded number of CG steps: ",iter
       end if
       
       call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
       if (.not. done_cell) call check_stop(done_cell, iter_cell)
    end do

    if (inode == ionode .and. iprint_gen + min_layer > 0) then
       write(io_lun, fmt='(/4x,a)') trim(prefix)//" final simulation box dimensions are: "
       write(io_lun, fmt='(4x,a,f12.5,1x,a2," x ",f12.5,1x,a2," x ",f12.5,1x,a2)') &
            prefixF(1:-2*min_layer), rcellx, d_units(dist_units), rcelly, d_units(dist_units), &
            rcellz, d_units(dist_units)
    end if

    deallocate(z_new_pos, y_new_pos, x_new_pos, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating _new_pos in control: ", &
         ni_in_cell,stat)
    deallocate(old_force, STAT=stat)
    deallocate(cg, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating cg in control: ", ni_in_cell,stat)
    call reg_dealloc_mem(area_general, 6*ni_in_cell, type_dbl)

  end subroutine full_double_loop
!!***

  ! Keeping this temporarily so that we can test the relative efficiency of
  ! alternating full ionic and full cell optimisation (full_cg_run_double_loop)
  ! and full ionic with single line minimisation cell optimisation (this routine)
  ! Use cell optimisation method 4 for this
  subroutine full_cg_run_double_loop_alt(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, only: iprint_gen, ni_in_cell, x_atom_cell,  &
                             y_atom_cell, z_atom_cell, id_glob,    &
                             atom_coord, area_general, iprint_MD,  &
                             IPRINT_TIME_THRES1,                   &
                             cell_en_tol, cell_stress_tol,         &
                             rcellx, rcelly, rcellz, min_layer
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin_cell, enthalpy, enthalpy_tolerance, safemin2
    use GenComms,      only: inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: tot_force, stress
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop, write_xsf
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  only: dump_InfoMatGlobal, dump_pos_and_matrices
    use md_control,    only: flag_write_xsf, target_pressure

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables for ionic relaxation
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma,gg1, &
                      enthalpy0, enthalpy1, dH, press
    integer        :: i,j,k,iter,length, jj, lun, stat
    logical        :: done_ions, done_cell
    type(cq_timer) :: tmr_l_iter
    real(double), allocatable, dimension(:,:) :: cg, old_force
    real(double), allocatable, dimension(:)   :: x_new_pos, y_new_pos,&
                                                 z_new_pos

    ! Local variables for cell relaxation
    real(double) :: new_rcellx, new_rcelly, new_rcellz, search_dir_x, &
                    search_dir_y, search_dir_z, stressx, stressy, stressz, &
                    RMSstress, newRMSstress, dRMSstress, search_dir_mean, &
                    mean_stress, energy2, volume, stress_diff, max_stress, stress_target
    integer      :: iter_cell


    if (inode==ionode .and. iprint_MD > 2) &
      write(io_lun,'(2x,a)') "control/full_cg_run_double_loop"

    allocate(cg(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
      call cq_abort("Error allocating cg in control: ",&
                    ni_in_cell, stat)
    allocate(old_force(3,ni_in_cell), STAT=stat)
    allocate(x_new_pos(ni_in_cell), y_new_pos(ni_in_cell), &
             z_new_pos(ni_in_cell), STAT=stat)
    if (stat/=0) &
      call cq_abort("Error allocating _new_pos in control: ", &
                    ni_in_cell,stat)
    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)
    search_dir_z = zero
    search_dir_x = zero
    search_dir_y = zero
    search_dir_mean = zero
    cg = zero
    done_cell = .false.
    length = 3 * ni_in_cell
    if (inode==ionode .and. iprint_gen > 0) then
      write(io_lun,'(4x,"Welcome to cg_run. Doing ",i4," steps")') MDn_steps
      write(io_lun,'(4x,"Force tolerance:    ",f19.8)') MDcgtol
      write(io_lun,'(4x,"Enthalpy tolerance: ",f19.8)') enthalpy_tolerance
    end if
    energy0 = total_energy
    energy1 = zero
    energy2 = zero
    dE = zero

    ! Find energy and forces
    min_layer = min_layer - 1
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.,0)
    call dump_pos_and_matrices
    call get_maxf(max)
    min_layer = min_layer + 1
    press = target_pressure/HaBohr3ToGPa
    ! Stress tolerance in Ha/Bohr3
    stress_target = cell_stress_tol/HaBohr3ToGPa
    enthalpy0 = enthalpy(energy0, press)
    dH = zero
    if (inode==ionode) then
      write(io_lun,'(/4x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," H: ", f16.8," dH: ",f12.8/)') & 
           0, max, enthalpy0, dH
    end if

    if (inode == ionode .and. iprint_gen > 0) then
       write(io_lun, fmt='(/4x,"Starting full cell optimisation"/)')
       write(io_lun,*)  "Initial cell dims", rcellx, rcelly, rcellz
    end if

    iter = 1
    iter_cell = 1
    ! Cell loop
    do while (.not. done_cell)
      ggold = zero
      old_force = zero
      energy1 = energy0
      enthalpy1 = enthalpy0
      done_ions = .false.

      do while (.not. done_ions) ! ionic loop
        call start_timer(tmr_l_iter, WITH_LEVEL)
        ! Construct ratio for conjugacy
        gg = zero
        gg1 = zero! PR
        do j = 1, ni_in_cell
          gg = gg +                              &
               tot_force(1,j) * tot_force(1,j) + &
               tot_force(2,j) * tot_force(2,j) + &
               tot_force(3,j) * tot_force(3,j)
          gg1 = gg1 +                              &
                tot_force(1,j) * old_force(1,j) + &
                tot_force(2,j) * old_force(2,j) + &
                tot_force(3,j) * old_force(3,j)
        end do
        if (abs(ggold) < 1.0e-6_double) then
          gamma = zero
        else
          gamma = (gg-gg1)/ggold ! PR - change to gg/ggold for FR
        end if
        if(gamma<zero) gamma = zero
        if (inode == ionode .and. iprint_MD > 2) &
          write (io_lun,*) ' CHECK :: Force Residual = ', &
                               for_conv * sqrt(gg)/ni_in_cell
        if (inode == ionode .and. iprint_MD > 2) &
          write (io_lun,*) ' CHECK :: gamma = ', gamma
        if (CGreset) then
          if (gamma > one) then
            if (inode == ionode) &
              write(io_lun,*) ' CG direction is reset! '
            gamma = zero
          end if
        end if
        if (inode == ionode) &
          write (io_lun, fmt='(/4x,"Atomic relaxation CG iteration: ",i5)') iter
        ggold = gg
        ! Build search direction
        do j = 1, ni_in_cell
          jj = id_glob(j)
          cg(1,j) = gamma*cg(1,j) + tot_force(1,jj)
          cg(2,j) = gamma*cg(2,j) + tot_force(2,jj)
          cg(3,j) = gamma*cg(3,j) + tot_force(3,jj)
          x_new_pos(j) = x_atom_cell(j)
          y_new_pos(j) = y_atom_cell(j)
          z_new_pos(j) = z_atom_cell(j)
        end do
        old_force = tot_force
        min_layer = min_layer - 1
        ! Minimise in this direction
        call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
                        energy1, fixed_potential, vary_mu)
        min_layer = min_layer + 1
        ! Output positions
        if (inode==ionode .and. iprint_gen > 1) then
          write(io_lun,'(4x,a4,a15)') "Atom", "Position"
          do i = 1, ni_in_cell
            write(io_lun,'(4x,i8,3f15.8)') i,atom_coord(:,i)
          end do
        end if
        call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
        if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)
        ! Analyse forces
        g0 = dot(length, tot_force, 1, tot_force, 1)
        call get_maxf(max)
        ! Output and energy changes
        enthalpy0 = enthalpy(energy0, press)
        enthalpy1 = enthalpy(energy1, press)
        dE = energy1 - energy0
        dH = enthalpy1 - enthalpy0
        newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
             (stress(2,2)*stress(2,2)) + &
             (stress(3,3)*stress(3,3)))/3)
        dRMSstress = (RMSstress - newRMSstress)
        enthalpy0 = enthalpy1

        if (inode==ionode) then
          write(io_lun,'(/4x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," H: ",f16.8," dH: ",f12.8/)') &
               iter, max, enthalpy1, en_conv*dH
          if (iprint_MD > 1) then
            write(io_lun,'(4x,"Force Residual:     ",f19.8," ",a2,"/",a2)') &
              for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
              d_units(dist_units)
            write(io_lun,'(4x,"Maximum force:      ",f19.8)') max
            write(io_lun,'(4x,"Force tolerance:    ",f19.8)') MDcgtol
            write(io_lun,'(4x,"Maximum stress         ",3f10.6)') &
              max_stress
            write(io_lun,'(4x,"Stress tolerance:   ",f19.8)') &
              cell_stress_tol
            write(io_lun,'(4x,"Enthalpy change:    ",f19.8," ",a2)') &
              en_conv*dH, en_units(energy_units)
            write(io_lun,'(4x,"Enthalpy tolerance: ",f19.8)') &
              enthalpy_tolerance
          end if
        end if

        iter = iter + 1
        if (iter > MDn_steps) then
          done_ions = .true.
          if (inode==ionode) &
            write (io_lun, fmt='(4x,"Exceeded number of CG steps: ",i6)') &
                   iter
        end if
        if (abs(max) < MDcgtol) then
          done_ions = .true.
          if (inode==ionode) &
            write (io_lun, &
              fmt='(4x,"Maximum force below threshold: ",f12.5)') max
        end if

        call dump_pos_and_matrices
        call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
        if (.not. done_ions) call check_stop(done_ions, iter)
      end do ! ionic loop

      enthalpy0 = enthalpy(energy0, press)
      volume = rcellx*rcelly*rcellz
      stressx = -stress(1,1)!/volume
      stressy = -stress(2,2)!/volume
      stressz = -stress(3,3)!/volume
      mean_stress = (stressx + stressy + stressz)/3
      RMSstress = sqrt(((stressx*stressx) + (stressy*stressy) + (stressz*stressz))/3)

      gamma = zero ! steepest descent only for cell iteration
      if (inode == ionode .and. iprint_gen > 0) &
        write (io_lun, fmt='(/4x,"Lattice vector relaxation iteration: ",i5)') iter_cell
      ggold = gg

      !Build search direction.
      ! If the volume constraint is set, there is only one search direction!
      ! This is the direction which minimises the mean stress.
      ! if (leqi(cell_constraint_flag, 'volume')) then
      !   search_dir_mean = gamma*search_dir_mean + mean_stress
      ! else
      search_dir_x = gamma*search_dir_x + stressx - press*volume
      search_dir_y = gamma*search_dir_y + stressy - press*volume
      search_dir_z = gamma*search_dir_z + stressz - press*volume
      ! end if

      new_rcellx = rcellx
      new_rcelly = rcelly
      new_rcellz = rcellz

      ! Minimise in this direction. Constraint information is also used within
      ! safemin_cell. Look in move_atoms.module.f90 for further information.
      min_layer = min_layer - 1
      call safemin_cell(new_rcellx, new_rcelly, new_rcellz, search_dir_x, &
                        search_dir_y, search_dir_z, search_dir_mean, press, &
                        enthalpy0, enthalpy1, fixed_potential, vary_mu)
      min_layer = min_layer + 1
      ! Output positions to UpdatedAtoms.dat
      if (inode==ionode .and. iprint_MD > 1) then
        write(io_lun,'(4x,a4,a15)') "Atom", "Position"
        do i = 1, ni_in_cell
          write(io_lun,'(4x,i8,3f15.8)') i,atom_coord(:,i)
        end do
      end if
      call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))

      ! Analyse Stresses and energies
      dH = enthalpy1 - enthalpy0
      call get_maxf(max)
      newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
                           (stress(2,2)*stress(2,2)) + &
                           (stress(3,3)*stress(3,3)))/three)
      dRMSstress = RMSstress - newRMSstress

      volume = rcellx*rcelly*rcellz
      max_stress = zero
      do i=1,3
        stress_diff = abs(press*volume + stress(i,i))/volume
        if (stress_diff > max_stress) max_stress = stress_diff
      end do

      if (inode==ionode) then
        write(io_lun,'(/4x,"GeomOpt + Iter: ",i4," MaxF: ",f12.8," H: " ,f16.8," dH: ",f12.8/)') &
             iter, max, enthalpy1, en_conv*dH
        if (iprint_MD > 1) then
          write(io_lun,'(4x,"Force Residual:     ",f19.8," ",a2,"/",a2)') &
            for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
            d_units(dist_units)
          write(io_lun,'(4x,"Maximum force:         ",f19.8)') max
          write(io_lun,'(4x,"Force tolerance:       ",f19.8)') MDcgtol
          write(io_lun,'(4x,"Maximum stress         ",e14.6," Ha/Bohr**3")') &
            max_stress
          write(io_lun,'(4x,"Simulation cell volume ",e14.6," Bohr**3")') volume
          write(io_lun,'(4x,"Maximum stress         ",f14.6," GPa")') &
               max_stress*HaBohr3ToGPa
          write(io_lun,'(4x,"Stress tolerance: ",f14.6," GPa")') &
            cell_stress_tol
          write(io_lun,'(4x,"Enthalpy change:       ",f19.8," ",a2)') &
            en_conv*dH, en_units(energy_units)
          write(io_lun,'(4x,"Enthalpy tolerance:    ",f19.8)') &
            enthalpy_tolerance
        end if
      end if

      iter = iter + 1
      iter_cell = iter_cell + 1

      energy0 = energy1
      enthalpy0 = enthalpy1

      ! Check exit criteria
      if (iter > MDn_steps) then
        done_cell = .true.
        if (inode==ionode) &
          write (io_lun, fmt='(4x,"Exceeded number of SD steps: ",i4)') iter
      end if
      ! If the force convergence criterion has been met after this step,
      ! we can assume the ionic positions are correct, THEN check cell
      ! convergence
      if (abs(dH) < enthalpy_tolerance) then
        if (abs(max) < MDcgtol .and. max_stress < stress_target) then
          done_cell = .true.
          if (inode==ionode) then
            write(io_lun,'(/4x,a,i4,a,i4,a)') "GeomOpt converged in ", &
              iter-iter_cell," ionic steps and ", iter_cell, " cell steps"
            write(io_lun, fmt='(4x,"Maximum force:      ",f19.8)') max
            write(io_lun, fmt='(4x,"Force tolerance:    ",f19.8)') MDcgtol
            write(io_lun, fmt='(4x,"Enthalpy change:    ",f19.8)') dH
            write(io_lun, fmt='(4x,"Enthalpy tolerance: ",f19.8)') &
              enthalpy_tolerance
            write(io_lun, fmt='(4x,"Maximum stress:     ",f19.8," GPa")') &
              max_stress*HaBohr3ToGPa
            write(io_lun, fmt='(4x,"Stress tolerance:   ",f19.8," GPa")') &
              cell_stress_tol
          end if
        end if
      else
        continue
      end if

      call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
      if (.not. done_cell) call check_stop(done_cell, iter_cell)
    end do

   if (inode == ionode .and. iprint_gen > 0) then
     write(io_lun, fmt='(/4x,"Finished full cell optimisation"/)')
     write(io_lun,*)  "Final cell dims", rcellx, rcelly, rcellz
   end if

    deallocate(z_new_pos, y_new_pos, x_new_pos, STAT=stat)
    if (stat /= 0) &
      call cq_abort("Error deallocating _new_pos in control: ", &
                       ni_in_cell,stat)
    deallocate(old_force, STAT=stat)
    deallocate(cg, STAT=stat)
    if (stat /= 0) &
      call cq_abort("Error deallocating cg in control: ", ni_in_cell,stat)
    call reg_dealloc_mem(area_general, 6*ni_in_cell, type_dbl)

  3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
  4   format(4x,'Enthalpy change: ',f15.8,' ',a2)
  5   format(4x,'Force Residual: ',f15.10,' ',a2,'/',a2)
  6   format(4x,'Maximum force component: ',f15.8,' ',a2,'/',a2)
  7   format(4x,3f15.8)

  end subroutine full_cg_run_double_loop_alt

  !!****f* control/full_cg_run_single_vector *
  !!
  !!  NAME
  !!   full_cg_run_single_vector
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Conjugate gradients opnimisation of cell vectors and ionic positions 
  !!   by minimisation of force vector defined in Pfrommer et al.
  !!   J. Comput. Phys. 131, 233 (1997)
  !!
  !!   Beware! This routine uses force ordering (or atom_coord ordering) for
  !!   the search direction.  The updates to the coordinates use this ordering
  !!   as well, so it is consistent, but this is different to other routines
  !!   in this module and requires re-ordering in safemin_full
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Z Raza
  !!  CREATION DATE
  !!   2018/02/06
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  subroutine full_cg_run_single_vector(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, only: iprint_gen, ni_in_cell, x_atom_cell,  &
                             y_atom_cell, z_atom_cell, id_glob,    &
                             atom_coord, rcellx, rcelly, rcellz,   &
                             area_general, iprint_MD,              &
                             IPRINT_TIME_THRES1,                   &
                             cell_stress_tol, min_layer
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin_full, cq_to_vector, enthalpy, &
                             enthalpy_tolerance, backtrack_linemin_full, &
                             cg_line_min, safe, backtrack
    use GenComms,      only: inode, ionode, cq_warn
    use GenBlas,       only: dot
    use force_module,  only: tot_force, stress
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop, write_xsf, return_prefix, print_atomic_positions
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  ONLY: dump_InfoMatGlobal, dump_pos_and_matrices
    use md_control,    only: flag_write_xsf, target_pressure

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma, gg1, &
                      press, rcellx_ref, rcelly_ref, rcellz_ref, &
                      enthalpy0, enthalpy1, dH, volume, max_stress, &
                      stress_diff, dRMSstress, grad_f_dot_p, RMSstress, newRMSstress
    integer        :: i,j,k,iter,length, jj, stat
    logical        :: done
    type(cq_timer) :: tmr_l_iter
    real(double), allocatable, dimension(:,:) :: cg, force, force_old, &
                                                 config
    real(double), dimension(3) :: one_plus_strain, strain, cell_ref

    character(len=20) :: subname = "full_cg_run: "
    character(len=120) :: prefix, prefixGO
    character(len=10) :: prefixF = "          "

    prefix = return_prefix(subname, min_layer)
    prefixGO = return_prefix("GeomOpt", min_layer)
    allocate(cg(3,ni_in_cell+1), STAT=stat)
    allocate(force(3,ni_in_cell+1), STAT=stat)
    allocate(force_old(3,ni_in_cell+1), STAT=stat)
    allocate(config(3,ni_in_cell+1), STAT=stat)
    length = 3*(ni_in_cell+1)
    call reg_alloc_mem(area_general, 5*length, type_dbl)
    if (inode == ionode) then
       if(cg_line_min==safe) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with safemin line minimisation"
       else if (cg_line_min==backtrack) then
          write(io_lun, fmt='(/4x,a)') &
               trim(prefix)//" starting relaxation with safemin line minimisation"
          write(io_lun, fmt='(4x,a)') &
               trim(prefix)//" as backtrack line minimisation not yet implemented"
       end if
       if(iprint_gen + min_layer>0) write(io_lun, fmt='(4x,a,f8.4,a2,"/",a2,a,i4)') &
            trim(prefix)//" tolerance: ",MDcgtol,en_units(energy_units), d_units(dist_units),&
            " maximum steps: ",MDn_steps
    end if
    cg = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3*(ni_in_cell+1)
    press = target_pressure/HaBohr3ToGPa
    energy0 = total_energy
    enthalpy0 = enthalpy(energy0, press)
    dH = zero

    ! reference cell to compute strain
    cell_ref(1) = rcellx
    cell_ref(2) = rcelly
    cell_ref(3) = rcellz

    ! Find energy and forces
    min_layer = min_layer - 1
    if (iprint_MD + min_layer > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.,0)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.,0)
    end if
    call dump_pos_and_matrices
    call get_maxf(max)
    min_layer = min_layer + 1
    volume = rcellx*rcelly*rcellz
    max_stress = zero
    do i=1,3
       stress_diff = abs(press*volume + stress(i,i))/volume
       if (stress_diff > max_stress) max_stress = stress_diff
    end do
    enthalpy0 = enthalpy(energy0, press)
    if (inode==ionode) then
       write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," H: " ,f16.8,x,a2," MaxS: ",f12.8," GPa")') &
            trim(prefixGO)//" - Iter: ",0, for_conv*max, en_units(energy_units), d_units(dist_units), &
            enthalpy0, en_units(energy_units), max_stress*HaBohr3ToGPa
    end if
    ! Check for trivial case where forces are converged
    if (abs(max) < MDcgtol .and. max_stress*HaBohr3ToGPa < cell_stress_tol) then
       done = .true.
       if (inode == ionode) then
          write(io_lun,'(4x,a,i4,a)') trim(prefixGO)//" converged in ", iter, " iterations"
          write(io_lun, fmt='(4x,a," Maximum stress below threshold:   ",f12.5," GPa")') &
               trim(prefix),max_stress*HaBohr3ToGPa
          write(io_lun, fmt='(4x,a," Maximum force below threshold:    ",f12.5," ",a2,"/",a2)') &
               trim(prefix),for_conv*max, en_units(energy_units), d_units(dist_units)
       end if
       return
    else if(abs(max) < MDcgtol) then
       if(inode==ionode) &
            call cq_warn(trim(prefix)," Maximum force below threshold; consider only optimising cell")
    end if
    iter = 1
    ggold = zero
    force_old = zero
    energy1 = energy0
    enthalpy1 = enthalpy0
    do while (.not. done)
      call start_timer(tmr_l_iter, WITH_LEVEL) ! Construct ratio for conjugacy
      RMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
           (stress(2,2)*stress(2,2)) + &
           (stress(3,3)*stress(3,3)))/3)
      ! Construct the vector to optimise (config) and its force vector (force)
      call cq_to_vector(force, config, cell_ref, press)
      gg = zero
      gg1 = zero
      do j = 1, ni_in_cell+1 ! extra row for cell "force"
        gg = gg + force(1,j)*force(1,j) + &
                  force(2,j)*force(2,j) + &
                  force(3,j)*force(3,j)
        gg1 = gg1 + force(1,j)*force_old(1,j) + &
                    force(2,j)*force_old(2,j) + &
                    force(3,j)*force_old(3,j)
      end do
      if (abs(ggold) < 1.0e-6_double) then
        gamma = zero
      else
         !gamma = (gg-gg1)/ggold ! PR - change to gg/ggold for FR
         gamma = gg/ggold ! PR - change to gg/ggold for FR
      end if
       if(gamma<zero) then
          if(inode==ionode .and. iprint_MD + min_layer > 2) &
               write(io_lun,fmt='(4x,a)') trim(prefix)//" gamma < 0.  Setting to gg/ggold."
          gamma = gg/ggold ! Default to FR if PR gives negative
       end if
      if (CGreset) then
        if (gamma > one) then
           if (inode == ionode .and. iprint_MD + min_layer > 1) &
                write(io_lun,fmt='(4x,a)') &
                trim(prefix)//" gamma>1, so resetting direction"
           gamma = zero
        end if
      end if
      if (inode == ionode .and. iprint_MD + min_layer > 2) then
         write(io_lun,fmt='(4x,a,f12.8)') &
              trim(prefix)//' gamma = ', gamma
      end if
      ggold = gg
      ! Build search direction - note that indexing is correct (cf safemin)
      ! as we update atom_coord in safemin_full
      do j=1,ni_in_cell+1
         do i=1,3
            cg(i,j) = gamma*cg(i,j) + force(i,j)
         end do
      end do
      force_old = force
      ! Minimise in this direction
      min_layer = min_layer - 1
      if(cg_line_min==safe) then
         call safemin_full(config, cg, cell_ref, enthalpy0, enthalpy1, &
              press, fixed_potential, vary_mu)
      else if(cg_line_min==backtrack) then
         call safemin_full(config, cg, cell_ref, enthalpy0, enthalpy1, &
              press, fixed_potential, vary_mu)
         ! Backtrack does not (yet) work - don't use ! 
         !grad_f_dot_p = -dot(3*ni_in_cell+3,force,1,cg,1)
         !call backtrack_linemin_full(config, cg, cell_ref, enthalpy0, enthalpy1, &
         !     press, grad_f_dot_p, fixed_potential, vary_mu)         
      end if
      min_layer = min_layer + 1
      ! Output positions
       if (inode == ionode .and. iprint_gen + min_layer > 1) then
          call print_atomic_positions
       end if
      call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
      if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)

      ! Analyse forces and stress
      g0 = dot(length-3, tot_force, 1, tot_force, 1)
      call get_maxf(max)
      volume = rcellx*rcelly*rcellz
      max_stress = zero
      do i=1,3
        stress_diff = abs(press*volume + stress(i,i))/volume
        if (stress_diff > max_stress) max_stress = stress_diff
      end do

      ! Output and energy changes
      dH = enthalpy1 - enthalpy0
      newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
           (stress(2,2)*stress(2,2)) + &
           (stress(3,3)*stress(3,3)))/3)
      dRMSstress = (RMSstress - newRMSstress)/volume
      if (inode==ionode) then
         write(io_lun,'(/4x,a,i4," MaxF: ",f12.8,x,a2,"/",a2," H: " ,f16.8,x,a2," MaxS: ",f12.8," GPa/")') &
              trim(prefixGO)//" - Iter: ",iter, for_conv*max, en_units(energy_units), d_units(dist_units), &
              enthalpy1, en_units(energy_units), max_stress*HaBohr3ToGPa
         if (iprint_MD + min_layer > 1) then
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force Residual:     ", &
                  for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Maximum force:      ", &
                  for_conv*max, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," ",a2,"/",a2)') prefixF(1:-2*min_layer)//"Force tolerance:    ", &
                  for_conv*MDcgtol, en_units(energy_units), d_units(dist_units)
             write(io_lun,'(4x,a,f19.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Maximum stress:     ", max_stress*HaBohr3ToGPa
             write(io_lun,'(4x,a,f19.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Stress tolerance:   ", cell_stress_tol
             write(io_lun,'(4x,a,f19.8," GPa")') &
                  prefixF(1:-2*min_layer)//"Change in stress:   ", dRMSstress*HaBohr3ToGPa
             write(io_lun,'(4x,a,f19.8," ",a2)') &
                  prefixF(1:-2*min_layer)//"Enthalpy change:    ", en_conv*dH, en_units(energy_units)
             write(io_lun,'(4x,a,f19.8," ",a2)') &
                  prefixF(1:-2*min_layer)//"Enthalpy tolerance: ", en_conv*enthalpy_tolerance, en_units(energy_units)
         end if
      end if

      enthalpy0 = enthalpy1
      if (iter > MDn_steps) then
         done = .true.
         if (inode == ionode) &
              write(io_lun, fmt='(4x,a,i6)') &
              trim(prefix)//" exceeded number of MD steps: ",iter
      end if
      if (abs(dH) < enthalpy_tolerance) then
        if (abs(max) < MDcgtol .and. max_stress*HaBohr3ToGPa < cell_stress_tol) then
          done = .true.
          if (inode==ionode) then
             write(io_lun,'(/4x,a,i4,a)') trim(prefixGO)//" converged in ", iter, " iterations"
             write(io_lun, fmt='(4x,a," Maximum stress below threshold:   ",f19.8," GPa")') &
                  trim(prefix),max_stress*HaBohr3ToGPa
             write(io_lun, fmt='(4x,a," Maximum force below threshold:    ",f19.8," ",a2,"/",a2)') &
                  trim(prefix),for_conv*max, en_units(energy_units), d_units(dist_units)
          end if
        end if
      end if

      call dump_pos_and_matrices
      iter = iter + 1
      call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
      if (.not. done) call check_stop(done, iter)
    end do
    ! Output final positions
    !    if(myid==0) call write_positions(parts)
    deallocate(config, force, force_old, cg, STAT=stat)
    if (stat /= 0) &
      call cq_abort("Error deallocating vectors in control: ", &
                   ni_in_cell,stat)
    call reg_dealloc_mem(area_general, 5*length, type_dbl)

  end subroutine full_cg_run_single_vector
  !!***

  !!****f* control/get_maxf *
  !!
  !!  NAME
  !!   get_maxf
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Find the maximum component of force in tot_force
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Z Raza
  !!  CREATION DATE
  !!   2019/06/07
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  subroutine get_maxf(maxf)

    use numbers
    use force_module,  only: tot_force
    use global_module, only: ni_in_cell

    ! Passed varaibles
    real(double), intent(out) :: maxf

    ! local variables
    integer      :: i,k

    maxf = zero
    do i=1,ni_in_cell
      do k=1,3
        if (abs(tot_force(k,i)) > maxf) maxf = abs(tot_force(k,i))
      end do
    end do

  end subroutine get_maxf
  !!***

end module control
