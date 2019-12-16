! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!  SOURCE
!!
module control

  use datatypes
  use global_module, only: io_lun
  use GenComms,      only: cq_abort
  use timer_module,  only: start_timer, stop_timer, cq_timer 
  use timer_module,  only: start_backtrace, stop_backtrace

  implicit none

  integer      :: MDn_steps 
  integer      :: MDfreq 
  real(double) :: MDtimestep 
  real(double) :: MDcgtol 
  logical      :: CGreset
  character(3) :: md_ensemble

  ! Area identification
  integer, parameter, private :: area = 9

  ! RCS tag for object file identification
  character(len=80), save, private :: &
       RCSid = "$Id$"
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
!!  SOURCE
!!
  subroutine control_run(fixed_potential, vary_mu, total_energy)

    use datatypes
    use dimens,               only: r_core_squared, r_h
    use GenComms,             only: my_barrier, cq_abort
    use pseudopotential_data, only: set_pseudopotential
    use force_module,         only: tot_force
    use minimise,             only: get_E_and_F
    use global_module,        only: runtype, flag_self_consistent, &
                                    flag_out_wf, flag_write_DOS, &
                                    flag_opt_cell, optcell_method
    use input_module,         only: leqi
    use store_matrix,         only: dump_pos_and_matrices

    implicit none

    ! Shared variables
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_energy

    ! Local variables
    integer, parameter:: backtrace_level = 0
    type(cq_timer)    :: backtrace_timer
    logical           :: NoMD
    integer           :: i, j

    real(double) :: spr, e_r_OLD

!****lat<$
    call start_backtrace(t=backtrace_timer,who='control_run',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    if ( leqi(runtype,'static') ) then
       if(.NOT.flag_self_consistent.AND.(flag_out_wf.OR.flag_write_DOS)) return
       call get_E_and_F(fixed_potential,vary_mu, total_energy,&
                        .true.,.true.,level=backtrace_level)
       !
    else if ( leqi(runtype, 'cg')    ) then
        if (flag_opt_cell) then
           select case(optcell_method)
           case(1)
             call cell_cg_run(fixed_potential, vary_mu, total_energy)
           case(2)
             call full_cg_run_double_loop(fixed_potential, vary_mu, &
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
  !!   - Added call for safemin2 necessary in reusibg L-matrix
  !!    2017/08/29 jack baker & dave
  !!     Removed rcellx references (redundant)
  !!   2017/11/10 14:06 dave
  !!    Removing dump_InfoGlobal calls
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
                             IPRINT_TIME_THRES1
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin2, backtrack_linemin
    use GenComms,      only: gsum, myid, inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: tot_force
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop, write_xsf
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  ONLY: dump_InfoMatGlobal, dump_pos_and_matrices
    use md_control,    only: flag_write_xsf

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
    if (myid == 0) &
         write (io_lun, fmt='(/4x,"Starting CG atomic relaxation"/)')
    cg = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3 * ni_in_cell
    if (myid == 0 .and. iprint_gen > 0) &
         write (io_lun, 2) MDn_steps, MDcgtol
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.)
    call dump_pos_and_matrices
    call get_maxf(max)
    if (inode==ionode) then
      write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," E: "   e16.8," dE: ",f12.8)') & 
           0, max, energy0, dE
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
       !test_dot = dot(3*ni_in_cell, cg, 1, tot_force, 1)
       !call gsum(test_dot)
       if(test_dot<zero) then
          if(inode==ionode) write(io_lun,fmt='(2x,"Reset direction ",f20.12)') test_dot
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
       !call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
       call backtrack_linemin(cg, energy0, energy1, fixed_potential, vary_mu, energy1)
       ! Output positions
       if (myid == 0 .and. iprint_gen > 1) then
          do i = 1, ni_in_cell
             write (io_lun, 1) i, atom_coord(1,i), atom_coord(2,i), &
                               atom_coord(3,i)
          end do
       end if
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)
       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       call get_maxf(max)
       ! Output and energy changes
       dE = energy0 - energy1

       if (inode==ionode) then
         write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," E: "e16.8," dE: ",f12.8)') & 
              iter, max, energy1, en_conv*dE
         if (iprint_MD > 1) then
           write(io_lun,'(4x,"Force Residual:     ",f20.10," ",a2,"/",a2)') &
             for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
             d_units(dist_units)
           write(io_lun,'(4x,"Maximum force:      ",f20.10)') max
           write(io_lun,'(4x,"Force tolerance:    ",f20.10)') MDcgtol
           write(io_lun,'(4x,"Energy change:      ",f20.10," ",a2)') &
             en_conv*dE, en_units(energy_units)
           ! write(io_lun,'(4x,"Energy tolerance:   ",f20.10)') &
           !   enthalpy_tolerance
           ! Why is there no energy threshold? - zamaan
         end if
       end if

       iter = iter + 1
       energy0 = energy1
       !energy1 = abs(dE)
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Exceeded number of MD steps: ",i6)') &
                     iter
       end if
       if (abs(max) < MDcgtol) then
          write(io_lun,'(2x,a,i4,a)') "GeomOpt converged in ", iter, " iterations"
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Maximum force below threshold: ",f12.5)') &
                     max
       end if

       call dump_pos_and_matrices
       
       call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
       if (.not. done) call check_stop(done, iter)
    end do
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

1   format(4x,'Atom ',i8,' Position ',3f15.8)
2   format(4x,'Welcome to cg_run. Doing ',i4,&
           ' steps with tolerance of ',f8.4,' ev/A')
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Energy change: ',f15.8,' ',a2)
5   format(4x,'Force Residual: ',f15.10,' ',a2,'/',a2)
6   format(4x,'Maximum force component: ',f15.8,' ',a2,'/',a2)
7   format(4x,3f15.8)
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
                              non_atomic_stress, flag_heat_flux
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
    use io_module,      only: read_fire, check_stop
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
    use md_model,       only: type_md_model
    use md_control,     only: type_thermostat, type_barostat, md_n_nhc, &
                              md_n_ys, md_n_mts, ion_velocity, lattice_vec, &
                              md_baro_type, target_pressure, md_ndof_ions, &
                              md_berendsen_equil, md_tau_T, md_tau_P, &
                              md_thermo_type, heat_flux

    use atoms,          only: distribute_atoms,deallocate_distribute_atom
    use global_module,  only: atom_coord_diff, iprint_MD, area_moveatoms

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical      :: vary_mu, fixed_potential
    real(double) :: total_energy
    
    ! Local variables
    integer       ::  iter, i, j, k, length, stat, i_first, i_last, md_ndof
    integer       :: nequil ! number of Berendsen equilibration steps - zamaan
    real(double)  :: energy1, energy0, dE, max, g0
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

    n_stop_qMD = 0
    final_call = 1
    energy0 = zero
    energy1 = zero
    dE = zero
    length = 3*ni_in_cell

    ! initialisation/contintuation when we're using Berendsen equilibration
    if (flag_MDcontinue) then
      if (MDinit_step >= md_berendsen_equil) then
        nequil = 0
      else
        nequil = md_berendsen_equil - MDinit_step
      end if
    else
      nequil = md_berendsen_equil 
    end if

    ! Initialise number of degrees of freedom
    md_ndof = 3*ni_in_cell
    md_ndof_ions = md_ndof
    do i=1,ni_in_cell
      do j=1,3
        if (flag_move_atom(j,i) .eqv. .false.) md_ndof = md_ndof-1
      end do
    end do

    if (inode==ionode .and. iprint_gen > 0) &
      write(io_lun,'(4x,"Welcome to md_run. Doing ",i6," steps")') &
            MDn_steps

    ! Thermostat/barostat initialisation
    call init_md(baro, thermo, mdl, md_ndof, nequil)
    call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                      mdl%ion_kinetic_energy)
    call baro%get_pressure_and_stress
    call mdl%get_cons_qty

    ! Find energy and forces
    if (flag_fire_qMD) then
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.)
    end if
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
    call dump_pos_and_matrices(index=0,MDstep=i_first,velocity=ion_velocity)

    if (flag_fire_qMD) then
       step_qMD = i_first ! SA 20150201
       done = .false.     ! SA 20150201
       ! reading FIRE parameters of the last run
       call read_fire(fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha)
    endif

    if (flag_heat_flux) &
      call get_heat_flux(atomic_stress, ion_velocity, heat_flux)

    if (.not. flag_MDcontinue) &
      call write_md_data(i_first-1, thermo, baro, mdl, nequil)

    do iter = i_first, i_last ! Main MD loop
       mdl%step = iter
       if (inode==ionode) &
            write(io_lun,fmt='(4x,"MD run, iteration ",i5)') iter

       if (flag_heat_flux) then
         atomic_stress = zero
         non_atomic_stress = zero
       end if

       call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                          mdl%ion_kinetic_energy)
       call baro%get_pressure_and_stress

       ! thermostat/barostat (MTTK splitting of Liouvillian)
       call integrate_pt(baro, thermo, mdl, ion_velocity)

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
       ! Reset-up
          if (md_ensemble(2:2) == 'p') then
            if (.not. leqi(md_baro_type, "berendsen") .or. nequil < 1) then
              call baro%update_cell
            end if
          end if
          if(flag_SFcoeffReuse) then
             call update_pos_and_matrices(updateSFcoeff,ion_velocity)
          else
             call update_pos_and_matrices(updateLorK,ion_velocity)
          endif

       if (flag_XLBOMD) call Do_XLBOMD(iter,MDtimestep)
       call update_H(fixed_potential)

       !2018.Jan.4 TM 
       !   We need to update flag_movable, since the order of x_atom_cell (or id_glob) 
       !   may change after the atomic positions are updated.
       call check_move_atoms(flag_movable)
       
       !! For Debuggging !!
       !     call dump_pos_and_matrices(index=2,MDstep=i_first)
       !! For Debuggging !!
       
       if (flag_fire_qMD) then
          call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .true.,iter)
          call check_stop(done, iter)   !2019/Nov/14
          call dump_pos_and_matrices(index=0,MDstep=iter,velocity=ion_velocity)
       else
          call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .false.,iter)
          call check_stop(done, iter)   !2019/Nov/14
          ! Here so that the kinetic stress is reported just after the 
          ! static stress - zamaan
          if (inode == ionode .and. iprint_MD > 2) then
            write(io_lun,fmt='(/4x,"Kinetic stress    ", 3f15.8,a3)') &
              baro%ke_stress(1,1), baro%ke_stress(2,2), baro%ke_stress(3,3), &
              en_units(energy_units)
          end if
          call dump_pos_and_matrices(index=0,MDstep=iter,velocity=ion_velocity)
          call vVerlet_v_dthalf(MDtimestep,ion_velocity,tot_force,flag_movable,second_call)
       end if
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

       ! Rescale the ionic positions for the berendsen barostat AFTER the 
       ! velocity updates to avoid rescaling velocities
       if (leqi(md_ensemble, 'npt')) then
         if (leqi(md_baro_type, 'berendsen') .or. nequil > 0) &
           call baro%propagate_berendsen(flag_movable)
       end if
       thermo%ke_ions = mdl%ion_kinetic_energy
       call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                          mdl%ion_kinetic_energy)
       call baro%get_pressure_and_stress

       !! For Debuggging !!
       !!call cq_abort(" STOP FOR DEBUGGING")
       !! For Debuggging !!

       ! thermostat/barostat (MTTK splitting of Liouvillian)
       call integrate_pt(baro, thermo, mdl, ion_velocity, second_call)

       ! Constrain velocity
       if (flag_RigidBonds) call correct_atomic_velocity(ion_velocity)
       if (flag_FixCOM) call zero_COM_velocity(ion_velocity)
       !%%! END of Evolve atoms

       ! Print out energy
       mdl%dft_total_energy = energy1
       if (inode==ionode) then
         write(io_lun,'(4x,"***MD step ",i6," KE: ",f18.8," &
                       &IntEnergy: ",f20.8," TotalEnergy: ",f20.8)') &
               iter, mdl%ion_kinetic_energy, mdl%dft_total_energy, &
               mdl%ion_kinetic_energy+mdl%dft_total_energy
         write(io_lun, fmt='(4x,"Kinetic Energy in K     : ",f15.8)') &
               mdl%ion_kinetic_energy / (three / two * ni_in_cell) / &
               fac_Kelvin2Hartree
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
       dE = energy0 - energy1
       energy0 = energy1
       energy1 = abs(dE)
       if(inode==ionode) then
         write (io_lun,'(4x,"Maximum force component : ",f15.8)') max
         write (io_lun,'(4x,"Force Residual          : ",f15.8)') &
                sqrt(g0/ni_in_cell)
         write (io_lun,'(4x,"Energy change           : ",f15.8)') dE
       end if

       ! Compute and print the conserved quantity and its components
       call thermo%get_nhc_energy
       call baro%get_box_energy(final_call)
       call mdl%get_cons_qty
       call mdl%print_md_energy()

       ! Output positions
       if (inode==ionode .and. iprint_gen > 1) then
          do i = 1, ni_in_cell
            write(io_lun,'(4x,"Atom ",i4," Position ",3f15.8)') &
                  i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
          enddo
       endif
       call my_barrier

       ! The kinetic component of stress changes after the second velocity
       ! update 
       call thermo%get_temperature_and_ke(baro, ion_velocity, &
                                          mdl%ion_kinetic_energy, &
                                          final_call)
       call baro%get_pressure_and_stress(final_call)
 
       if (nequil > 0) then
         nequil = nequil - 1
         if (nequil == 0) then
           if (inode==ionode) &
             write (io_lun, '(4x,a)') "Berendsen equilibration finished, &
                                      &starting extended system dynamics."
           call init_md(baro, thermo, mdl, md_ndof, nequil, second_call)
         end if
       end if

       ! Compute heat flux
       if (flag_heat_flux) &
         call get_heat_flux(atomic_stress, ion_velocity, heat_flux)

       ! Write all MD data and checkpoints to disk
       call write_md_data(iter, thermo, baro, mdl, nequil)

       !call check_stop(done, iter) ! moved above. 2019/Nov/14 tsuyoshi
       if (flag_fire_qMD.OR.flag_quench_MD) then
          if (abs(max) < MDcgtol) then
             if ((iter - step_qMD) > 1) then
                n_stop_qMD = 0
             end if
             n_stop_qMD = n_stop_qMD + 1
             step_qMD = iter
             if (inode==ionode) then
                write(io_lun,fmt='(4x,i4,4x,"Maximum force below &
                     &threshold: ",f12.6)') n_stop_qMD, max
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

  !!****m* control/init_md *
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
  !!  SOURCE
  !!  
  subroutine init_md(baro, thermo, mdl, md_ndof, nequil, second_call)

    use numbers
    use input_module,   only: leqi
    use io_module,      only: read_velocity
    use md_model,       only: type_md_model
    use GenComms,       only: inode, ionode, gcopy
    use memory_module,  only: reg_alloc_mem, type_dbl
    use global_module,  only: rcellx, rcelly, rcellz, temp_ion, ni_in_cell, &
                              flag_MDcontinue, flag_read_velocity, &
                              flag_MDdebug, iprint_MD, flag_atomic_stress, &
                              atomic_stress, area_moveatoms, &
                              flag_full_stress, flag_heat_flux
    use move_atoms,     only: init_velocity
    use md_control,     only: type_thermostat, type_barostat, md_tau_T, &
                              md_tau_P, md_n_nhc, md_n_ys, md_n_mts, &
                              md_nhc_mass, ion_velocity, lattice_vec, &
                              md_thermo_type, md_baro_type, target_pressure, &
                              md_ndof_ions, md_tau_P_equil, md_tau_T_equil, &
                              write_md_checkpoint, read_md_checkpoint, &
                              flag_extended_system, heat_flux

    ! passed variableariables
    type(type_barostat), intent(inout)    :: baro
    type(type_thermostat), intent(inout)  :: thermo
    type(type_md_model), intent(inout)    :: mdl
    integer, intent(in)                   :: md_ndof
    integer, intent(in)                   :: nequil
    logical, optional                     :: second_call

    ! local variables
    character(50) :: file_velocity='velocity.dat'
    integer       :: stat

    if (inode==ionode .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "Welcome to init_md"

    if (.not. allocated(ion_velocity)) then
      allocate(ion_velocity(3,ni_in_cell), STAT=stat)
      if (stat /= 0) &
           call cq_abort("Error allocating velocity in init_md: ", &
                         ni_in_cell, stat)
      call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)
    end if

    if (.not. flag_MDcontinue) then
      ! I've moved the velocity initialisation here to make reading the new
      ! unified md checkpoint file easier - zamaan
      ion_velocity = zero
      if (flag_read_velocity) then
        call read_velocity(ion_velocity, file_velocity)
      else
        if(temp_ion > RD_ERR) then
           call init_velocity(ni_in_cell, temp_ion, ion_velocity)
        end if
      end if
    end if

    if (.not. present(second_call)) then
    ! Initialise the model only once per run
      lattice_vec = zero
      lattice_vec(1,1) = rcellx
      lattice_vec(2,2) = rcelly
      lattice_vec(3,3) = rcellz
      call mdl%init_model(md_ensemble, MDtimestep, thermo, baro)
    end if

    select case(md_ensemble)
    case('nve')
      ! Just for computing temperature
      call thermo%init_thermo('none', 'none', MDtimestep, md_ndof, md_tau_T, &
                              mdl%ion_kinetic_energy)
      call baro%init_baro('none', MDtimestep, md_ndof, ion_velocity, &
                          md_tau_P, mdl%ion_kinetic_energy) !to get the pressure
    case('nvt')
      if (nequil > 0) then ! Equilibrate using Berendsen?
        if (inode==ionode) then
          write (io_lun, '(4x,"Equilibrating using Berendsen &
                          &baro/thermostat for ",i8," steps")') nequil
        end if
        call thermo%init_thermo('berendsen', 'none', MDtimestep, md_ndof, &
                                md_tau_T, mdl%ion_kinetic_energy)
      else
        call thermo%init_thermo(md_thermo_type, 'none', MDtimestep, md_ndof, &
                                md_tau_T, mdl%ion_kinetic_energy)
      end if
      call baro%init_baro('none', MDtimestep, md_ndof, ion_velocity, &
                          md_tau_P, mdl%ion_kinetic_energy) !to get the pressure
    case('npt')
      if (leqi(md_baro_type, 'berendsen')) then
        md_thermo_type = 'berendsen'
        call thermo%init_thermo('berendsen', 'berendsen', MDtimestep, &
                                md_ndof, md_tau_T, mdl%ion_kinetic_energy)
        call baro%init_baro('berendsen', MDtimestep, md_ndof, &
                            ion_velocity, md_tau_P, mdl%ion_kinetic_energy)
      else
        md_thermo_type = 'nhc'
        if (nequil > 0) then ! Equilibrate using Berendsen?
          if (inode == ionode) then
            write (io_lun, '(4x,"Equilibrating using Berendsen &
                            &baro/thermostat for ",i8," steps")') nequil
          end if
          call thermo%init_thermo('berendsen', 'berendsen', MDtimestep, &
                                  md_ndof, md_tau_T_equil, &
                                  mdl%ion_kinetic_energy)
          call baro%init_baro('berendsen', MDtimestep, md_ndof, &
                              ion_velocity, md_tau_P_equil, &
                              mdl%ion_kinetic_energy)
        else
          call thermo%init_thermo(md_thermo_type, md_baro_type, MDtimestep, &
                                  md_ndof, md_tau_T, mdl%ion_kinetic_energy)
          call baro%init_baro(md_baro_type, MDtimestep, md_ndof, &
                              ion_velocity, md_tau_P, mdl%ion_kinetic_energy)
        end if
      end if
    end select

    ! N.B. atomic stress is allocated in initialisation_module/initialise! - zamaan

    if (flag_MDcontinue) call read_md_checkpoint(thermo, baro)

  end subroutine init_md

  !!****m* md_control/end_md *
  !!  NAME
  !!   end_md
  !!  PURPOSE
  !!   Deallocate MD arrays that are not deallocated elsewhere
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/05/08
  !!  SOURCE
  !!  
  subroutine end_md(th, baro)

    use GenComms,         only: inode, ionode
    use global_module,    only: area_moveatoms, atomic_stress, &
                                flag_atomic_stress, flag_MDdebug, &
                                iprint_MD, ni_in_cell, area_moveatoms
    use memory_module,    only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use md_control,       only: type_thermostat, type_barostat, &
                                md_nhc_mass, md_nhc_cell_mass, &
                                flag_nhc, ion_velocity

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

    deallocate(ion_velocity, STAT=stat)
    if (stat /= 0) call cq_abort("Error deallocating velocity in md_run: ", &
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

  !!****m* control/integrate_pt *
  !!  NAME
  !!   integrate_pt
  !!  PURPOSE
  !!   Integrate the thermostat and barostat using the
  !!    Martyna/Tobias/Tuckerman/Klein splitting of the Liouvillian
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/04/23 11:49
  !!  SOURCE
  !!  
  subroutine integrate_pt(baro, thermo, mdl, velocity, second_call)

    use global_module,    only: iprint_MD
    use GenComms,         only: inode, ionode
    use md_model,         only: type_md_model
    use md_control,       only: type_thermostat, type_barostat, &
                                md_thermo_type, md_baro_type

    ! passed variables
    type(type_barostat), intent(inout)          :: baro
    type(type_thermostat), intent(inout)        :: thermo
    type(type_md_model), intent(inout)          :: mdl
    real(double), dimension(:,:), intent(inout) :: velocity
    logical, optional                           :: second_call

    ! local variables
  
    if (inode==ionode .and. iprint_MD > 1) then
      if (present(second_call)) then
        write(io_lun,'(2x,a)') "Welcome to integrate_pt, second call"
      else
        write(io_lun,'(2x,a)') "Welcome to integrate_pt"
      end if
    end if

    select case(md_ensemble)
    case('nvt')
      select case(thermo%thermo_type)
      case('nhc')
        call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
        if (present(second_call)) call thermo%get_nhc_energy
      case('berendsen')
        if (present(second_call)) then
          call thermo%berendsen_v_rescale(velocity)
        else
          call thermo%get_berendsen_thermo_sf(MDtimestep)
        end if
     end select
   case('npt')
      select case(baro%baro_type)
      case('berendsen')
        if (present(second_call)) then
          call thermo%berendsen_v_rescale(velocity)
          ! Berendsen barostat propagation occurs after second
          ! vVerlet_v_dthalf step
        else
          call thermo%get_berendsen_thermo_sf(MDtimestep)
          call baro%get_berendsen_baro_sf(MDtimestep)
        end if
!            call baro%propagate_berendsen(flag_movable)
      case('iso-mttk')
        if (mdl%nequil  > 0) then
          if (present(second_call)) then
            call thermo%berendsen_v_rescale(velocity)
          else
            call thermo%get_berendsen_thermo_sf(MDtimestep)
            call baro%get_berendsen_baro_sf(MDtimestep)
          end if
        else
          call baro%propagate_npt_mttk(thermo, mdl%ion_kinetic_energy, &
                                       velocity)
        end if
      case('ortho-mttk')
      case('mttk')
      case('iso-ssm')
        if (present(second_call)) then
          call baro%couple_box_particle_velocity(thermo, velocity)
          call thermo%get_temperature_and_ke(baro, velocity, &
                                             mdl%ion_kinetic_energy)
          call baro%get_pressure_and_stress
          call baro%integrate_box(thermo)
          call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
        else
          call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
          call thermo%get_temperature_and_ke(baro, velocity, &
                                             mdl%ion_kinetic_energy)
          call baro%get_pressure_and_stress
          call baro%integrate_box(thermo)
          call baro%couple_box_particle_velocity(thermo, velocity)
        end if
      case('ortho-ssm')
        if (present(second_call)) then
          call baro%couple_box_particle_velocity(thermo, velocity)
          call thermo%get_temperature_and_ke(baro, velocity, &
                                             mdl%ion_kinetic_energy)
          call baro%get_pressure_and_stress
          call baro%integrate_box(thermo)
          call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
        else
          call thermo%integrate_nhc(baro, velocity, mdl%ion_kinetic_energy)
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

!!****m* control/update_pos_and_box *
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
  use global_module, only: iprint_MD
  use md_control,    only: type_thermostat, type_barostat, ion_velocity, &
                           md_baro_type

  ! passed variables
  type(type_barostat), intent(inout)  :: baro
  integer, intent(in)                 :: nequil
  logical, dimension(:), intent(in)   :: flag_movable

  ! local variables

    if (inode==ionode .and. iprint_MD > 1) &
        write(io_lun,'(2x,a)') "Welcome to update_pos_and_box"

  if (leqi(md_ensemble, 'npt')) then
    if (leqi(md_baro_type, "berendsen") .or. nequil > 0) then
      call vVerlet_r_dt(MDtimestep,ion_velocity,flag_movable)
    else
      ! Modified velocity Verlet position step for NPT ensemble
      select case(md_baro_type)
      case('iso-ssm')
        call baro%propagate_box_ssm
        call vVerlet_r_dt(MDtimestep,ion_velocity,flag_movable)
        call baro%propagate_box_ssm
      case('ortho-ssm')
        call baro%propagate_box_ssm
        call vVerlet_r_dt(MDtimestep,ion_velocity,flag_movable)
        call baro%propagate_box_ssm
      case('iso-mttk')
        call baro%propagate_r_mttk(MDtimestep, ion_velocity, flag_movable)
        call baro%propagate_box_mttk(MDtimestep)
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

!!****m* control/get_heat_flux *
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
  use global_module, only: iprint_MD, ni_in_cell, flag_heat_flux

  ! Passed variables
  real(double), dimension(3,3,ni_in_cell), intent(in)   :: atomic_stress
  real(double), dimension(3,ni_in_cell), intent(in)     :: velocity
  real(double), dimension(3), intent(out)               :: heat_flux

  ! local variables
  integer :: i

  if (inode==ionode .and. iprint_MD > 1) &
    write(io_lun,'(2x,a)') "Welcome to get_heat_flux"

  heat_flux = zero
  do i=1,ni_in_cell
    heat_flux = heat_flux + matmul(atomic_stress(:,:,i), velocity(:,i))
  end do

end subroutine get_heat_flux
!!*****

!!****m* control/write_md_data *
!!  NAME
!!   write_md_data
!!  PURPOSE
!!   Write MD data to various files at the end of an ionic step
!!  AUTHOR
!!   Zamaan Raza
!!  CREATION DATE
!!   2018/08/11 10:27
!!  SOURCE
!!  
subroutine write_md_data(iter, thermo, baro, mdl, nequil)

  use GenComms,      only: inode, ionode
  use io_module,     only: write_xsf
  use global_module, only: iprint_MD, flag_baroDebug, flag_thermoDebug, &
                           flag_heat_flux
  use md_model,      only: type_md_model, md_tdep
  use md_control,    only: type_barostat, type_thermostat, &
                           write_md_checkpoint, flag_write_xsf, &
                           md_thermo_file, md_baro_file, &
                           md_trajectory_file, md_frames_file, &
                           md_stats_file, md_heat_flux_file

  ! Passed variables
  type(type_barostat), intent(inout)    :: baro
  type(type_thermostat), intent(inout)  :: thermo
  type(type_md_model), intent(inout)    :: mdl
  integer, intent(in)                   :: iter, nequil

    if (inode==ionode .and. iprint_MD > 1) &
        write(io_lun,'(2x,a)') "Welcome to write_md_data"

  call write_md_checkpoint(thermo, baro)
  call mdl%dump_stats(md_stats_file, nequil)
  if (flag_write_xsf) call write_xsf(md_trajectory_file, iter)
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
       dE = energy0 - energy1
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
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., &
                     .false.)
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
          call get_E_and_F(fixed_potential, vary_mu, energy1, .true., &
               .false.)
          call dump_pos_and_matrices
       else
          call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
               energy1, fixed_potential, vary_mu, energy1)
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
       call get_E_and_F(fixed_potential, vary_mu, energy1, .true., &
                        .false.)
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
       dE = energy0 - energy1
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
                              safemin2, update_H, update_pos_and_matrices, backtrack_linemin
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
    use dimens,        only: r_super_x, r_super_y, r_super_z

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
                    mx_pulay, i_first, i_last, &
                    nfile, symm, iter_low, iter_high, this_iter
    logical      :: done!, simpleStep

    step = MDtimestep
    !simpleStep = .false.
    mx_pulay = 5
    allocate(posnStore(3,ni_in_cell,mx_pulay), &
             forceStore(3,ni_in_cell,mx_pulay), STAT=stat)
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
    allocate(alpha(mx_pulay), beta(mx_pulay), rho(mx_pulay))
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
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., &
                     .false.)
    call dump_pos_and_matrices
    iter = 0
    ggold = zero
    energy1 = energy0
    ! Store positions and forces on entry (primary set only)
    !do i=1,ni_in_cell
       !jj = id_glob(i)
       !cg_new(:,i) = tot_force(:,jj)
    !end do
    cg_new = -tot_force ! The L-BFGS is in terms of grad E
    do while (.not. done)
       test_dot = dot(length, cg_new, 1, tot_force, 1)
       if(inode==ionode.AND.iprint_MD>1) write(io_lun,fmt='(2x,"grad f dot p ",f15.7)') test_dot
       if(test_dot>zero) then
          iter = 0
          if(inode==ionode.AND.iprint_MD>1) write(io_lun,fmt='(2x,"Resetting history")')
          cg_new = -tot_force
       end if
       ! Book-keeping
       iter = iter + 1
       npmod = mod(iter, mx_pulay)
       if(npmod==0) npmod = mx_pulay
       do i=1,ni_in_cell
          jj = id_glob(i)
          posnStore (1,jj,npmod) = x_atom_cell(i)!i)
          posnStore (2,jj,npmod) = y_atom_cell(i)!i)
          posnStore (3,jj,npmod) = z_atom_cell(i)!i)
          forceStore(:,jj,npmod) = -tot_force(:,jj)
          x_new_pos(i) = x_atom_cell(i)
          y_new_pos(i) = y_atom_cell(i)
          z_new_pos(i) = z_atom_cell(i)
          cg(:,i) = -cg_new(:,jj) ! Search downhill
       end do
       ! Set up limits for sums
       if(iter>mx_pulay) then
          iter_low = iter-mx_pulay
       else
          iter_low = 1
       end if
       if (myid == 0 .and. iprint_MD > 2) &
            write(io_lun,fmt='(2x,"L-BFGS iteration ",i4)') iter
       ! Line search
       call backtrack_linemin(cg, energy0, energy1, fixed_potential, vary_mu, energy1)
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
          cg_new = -tot_force ! The L-BFGS is in terms of grad E
       end do
       ! Build q
       do i=iter, iter_low, -1
          ! Indexing
          this_iter = mod(i,mx_pulay)
          if(this_iter==0) this_iter = mx_pulay
          ! parameters
          rho(this_iter) = one/dot(length,forceStore(:,:,this_iter),1,posnStore(:,:,this_iter),1)
          alpha(this_iter) = rho(this_iter)*dot(length,posnStore(:,:,this_iter),1,cg_new,1)
          cg_new = cg_new - alpha(this_iter)*forceStore(:,:,this_iter)
       end do
       ! Apply preconditioning in future
       ! Build z
       do i=iter_low, iter
          ! Indexing
          this_iter = mod(i,mx_pulay)
          if(this_iter==0) this_iter = mx_pulay
          ! Build
          beta(this_iter) = rho(this_iter)*dot(length,forceStore(:,:,this_iter),1,cg_new,1)
          cg_new = cg_new + (alpha(this_iter) - beta(this_iter))*posnStore(:,:,this_iter)
       end do
       ! Analyse forces
       g0 = dot(length, cg_new, 1, tot_force, 1)
       write(io_lun,*) 'z dot -grad E is ',g0
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
       if (.not. done) call check_stop(done, iter)
       dE = energy0 - energy1
       energy0 = energy1
    end do
    ! Output final positions
    !    if (myid == 0) call write_positions(parts)
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
  end subroutine lbfgs
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
  !!  01/06/17 J.S.Baker
  !!  MODIFICATION HISTORY
  !!
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
                             cell_constraint_flag, cell_stress_tol
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin_cell, enthalpy, enthalpy_tolerance
    use GenComms,      only: gsum, myid, inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: stress, tot_force
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use io_module,      only: leqi
    use dimens, ONLY: r_super_x, r_super_y, r_super_z
    use store_matrix,  only: dump_pos_and_matrices
    use md_control,    only: target_pressure, fac_HaBohr32GPa

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
                    stress_diff, volume

    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)

    if (myid == 0 .and. iprint_gen > 0) &
         write (io_lun, fmt='(/4x,"Starting CG lattice vector relaxation"/)')
    search_dir_x = zero
    search_dir_y = zero
    search_dir_z = zero
    search_dir_mean = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3
    if (myid == 0 .and. iprint_gen > 0) &
         write (io_lun, 2) MDn_steps, cell_en_tol,en_units(energy_units)
    energy0 = total_energy
    energy1 = zero
    dE = zero
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.)
    iter = 1
    reset_iter = 1
    ggold = zero
    energy1 = energy0

    press = target_pressure/fac_HaBohr32GPa
    enthalpy0 = enthalpy(energy0, press)
    enthalpy1 = enthalpy0
    dH = zero
    max_stress = zero
    volume = rcellx*rcelly*rcellz
    do i=1,3
      stress_diff = abs(press + stress(i,i))/volume
      if (stress_diff > max_stress) max_stress = stress_diff
    end do
    if (inode==ionode) then
      write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxStr: ",f12.8," H: "e16.8," dH: ",f12.8)') &
           0, max_stress, enthalpy0, zero
    end if

    call dump_pos_and_matrices(index=0,MDstep=iter)
    do while (.not. done)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       volume = rcellx*rcelly*rcellz
       stressx = -stress(1,1)/volume
       stressy = -stress(2,2)/volume
       stressz = -stress(3,3)/volume
       mean_stress = (stressx + stressy + stressz)/3
       RMSstress = sqrt(((stressx*stressx) + (stressy*stressy) + (stressz*stressz))/3)

       ! Construct ratio for conjugacy. Constraints are initially applied within
       ! get_gamma_cell_cg.
       call get_gamma_cell_cg(ggold, gg, gamma, stressx, stressy, stressz)
       if (myid == 0 .and. iprint_gen > 0) &
           write(io_lun, 3) iter, gamma

       if (inode == ionode .and. iprint_MD > 2) &
            write (io_lun,*) ' CHECK :: energy residual = ', &
                               dE
       if (inode == ionode .and. iprint_MD > 2) &
            write (io_lun,*) ' CHECK :: gamma = ', gamma

       if (CGreset) then
          if (gamma > one .or. reset_iter > length) then
             if (inode == ionode .and. iprint_gen > 0) &
                  write(io_lun,*) ' CG direction is reset to steepest descents! '
             gamma = zero
             reset_iter = 0
          end if
       end if
       if (inode == ionode .and. iprint_gen > 0) &
            write (io_lun, fmt='(/4x,"Lattice vector relaxation CG iteration: ",i5)') iter
       ggold = gg

       !Build search direction.
       ! If the volume constraint is set, there is only one search direction!
       ! This is the direction which minimises the mean stress.
       if (leqi(cell_constraint_flag, 'volume')) then
         search_dir_mean = gamma*search_dir_mean + mean_stress
       else
         search_dir_x = gamma*search_dir_x + stressx - press
         search_dir_y = gamma*search_dir_y + stressy - press
         search_dir_z = gamma*search_dir_z + stressz - press
       end if

       if (inode == ionode .and. iprint_gen > 0) &
           write(io_lun,*)  "Initial cell dims", rcellx, rcelly, rcellz

       new_rcellx = rcellx
       new_rcelly = rcelly
       new_rcellz = rcellz

       ! Minimise in this direction. Constraint information is also used within
       ! safemin_cell. Look in move_atoms.module.f90 for further information.
       call safemin_cell(new_rcellx, new_rcelly, new_rcellz, search_dir_x, &
                         search_dir_y, search_dir_z, search_dir_mean, press, &
                         enthalpy0, enthalpy1, fixed_potential, vary_mu)
       ! Output positions to UpdatedAtoms.dat
       if (myid == 0 .and. iprint_gen > 1) then
          do i = 1, ni_in_cell
             write (io_lun, 1) i, atom_coord(1,i), atom_coord(2,i), &
                               atom_coord(3,i)
          end do
       end if
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))

       ! Analyse Stresses and energies
       dH = enthalpy0 - enthalpy1
       newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
                            (stress(2,2)*stress(2,2)) + &
                            (stress(3,3)*stress(3,3)))/3)
       dRMSstress = RMSstress - newRMSstress

       enthalpy0 = enthalpy1
       ! Check exit criteria
       volume = rcellx*rcelly*rcellx
       max_stress = zero
       do i=1,3
         stress_diff = abs(press + stress(i,i))/volume
         if (stress_diff > max_stress) max_stress = stress_diff
       end do

       reset_iter = reset_iter +1

       if (inode==ionode) then
         write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxStr: ",f12.8," H: "e16.8," dH: ",f12.8)') &
              iter, max_stress, enthalpy1, en_conv*dH
         if (iprint_MD > 1) then
           write(io_lun,'(4x,"Maximum stress         ",3f10.6)') &
             max_stress
           write(io_lun,'(4x,"Stress tolerance:   ",f20.10)') &
             cell_stress_tol
           write(io_lun,'(4x,"Change in stress: ",f15.8," ",a2)') max_stress, &
             en_units(energy_units)
           write(io_lun,'(4x,"Enthalpy change:    ",f20.10," ",a2)') &
             en_conv*dH, en_units(energy_units)
           write(io_lun,'(4x,"Enthalpy tolerance: ",f20.10)') &
             enthalpy_tolerance
         end if
       end if

       iter = iter + 1
       if (myid == 0 .and. iprint_MD > 0) then
           write (io_lun, 4) en_conv*dH, en_units(energy_units)
           write (io_lun, 5) dRMSstress, "Ha"
           write(io_lun,'(4x,"Change in stress: ",f15.8," ",a2)') max_stress, &
                 en_units(energy_units)
       end if

       ! First exit is if too many steps have been taken. Default is 50.
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Exceeded number of MD steps: ",i4)') &
                     iter
       end if

       ! Second exit is if the desired energy tolerence has ben reached
       ! Will replace with stress tolerance when more reliable
       if (abs(dH)<enthalpy_tolerance .and. max_stress < cell_stress_tol) then
          if (inode==ionode) &
               write(io_lun,'(2x,a,i4,a)') "GeomOpt converged in ", &
               iter, " iterations"
          done = .true.
          if (myid == 0 .and. iprint_gen > 0) &
               write (io_lun, fmt='(4x,"Enthalpy change below threshold: ",f20.10)') &
                     dH
               write (io_lun, fmt='(4x,"Stress change below threshold:   ",f20.10)') &
                     max_stress
       end if

       call stop_print_timer(tmr_l_iter, "a CG iteration", IPRINT_TIME_THRES1)
       if (.not. done) call check_stop(done, iter)
    end do

    if (myid == 0 .and. iprint_gen > 0) then
        write(io_lun, fmt='("Final simulation box dimensions are: ")')
        write(io_lun, fmt='(2x,"a = ",f12.5,1x,a2)') rcellx, d_units(dist_units)
        write(io_lun, fmt='(2x,"b = ",f12.5,1x,a2)') rcelly, d_units(dist_units)
        write(io_lun, fmt='(2x,"c = ",f12.5,1x,a2)') rcellz, d_units(dist_units)
    end if

    call reg_dealloc_mem(area_general, 6*ni_in_cell, type_dbl)

1   format(4x,'Atom ',i8,' Position ',3f15.8)
2   format(4x,'Welcome to cell_cg_run. Doing ',i4,&
           ' steps with tolerance of ',f12.5,a2)
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Enthalpy change: ',f15.8,' ',a2)
5   format(4x,'RMS Stress change: ',f15.8,' ',a2)

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
    if(inode==ionode.AND.iprint_gen>1) write(io_lun,fmt='(2x,"Stress gg is ",f12.5)') gg

  end subroutine get_gamma_cell_cg
!!***

  !!****f* control/full_cg_run_double_loop *
  !!
  !!  NAME
  !!   full_cg_run_double_loop
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Simple implementation of cell + ionic geometry optimisation using a
  !!   nested loop: the outer loop does one cell steepest descent step, the
  !!   inner does a full ionic conjugate gradients optimisation
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
  subroutine full_cg_run_double_loop(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, only: iprint_gen, ni_in_cell, x_atom_cell,  &
                             y_atom_cell, z_atom_cell, id_glob,    &
                             atom_coord, area_general, iprint_MD,  &
                             IPRINT_TIME_THRES1,                   &
                             cell_en_tol, cell_stress_tol,         &
                             rcellx, rcelly, rcellz
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin2, safemin_cell, enthalpy, &
                             enthalpy_tolerance
    use GenComms,      only: inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: tot_force, stress
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop, write_xsf
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  only: dump_InfoMatGlobal, dump_pos_and_matrices
    use md_control,    only: flag_write_xsf, target_pressure, &
                             fac_HaBohr32GPa

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
                    mean_stress, energy2, volume, stress_diff, max_stress
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
      write(io_lun,'(4x,"Force tolerance:    ",f20.10)') MDcgtol
      write(io_lun,'(4x,"Enthalpy tolerance: ",f20.10)') enthalpy_tolerance
    end if
    energy0 = total_energy
    energy1 = zero
    energy2 = zero
    dE = zero

    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.)
    call dump_pos_and_matrices
    call get_maxf(max)
    press = target_pressure/fac_HaBohr32GPa
    enthalpy0 = enthalpy(energy0, press)
    dH = zero
    if (inode==ionode) then
      write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," H: " e16.8," dH: ",f12.8)') & 
           0, max, enthalpy0, dH
    end if

   if (inode == ionode .and. iprint_gen > 0) then
     write (io_lun, fmt='(/4x,"Starting full cell optimisation"/)')
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
        ! Minimise in this direction
          call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
                        energy1, fixed_potential, vary_mu, energy1)
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
        dE = energy0 - energy1
        dH = enthalpy0 - enthalpy1
        dH = enthalpy0 - enthalpy1
        energy0 = energy1

        if (inode==ionode) then
          write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," H: "e16.8," dH: ",f12.8)') &
               iter, max, enthalpy1, en_conv*dH
          if (iprint_MD > 1) then
            write(io_lun,'(4x,"Force Residual:     ",f20.10," ",a2,"/",a2)') &
              for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
              d_units(dist_units)
            write(io_lun,'(4x,"Maximum force:      ",f20.10)') max
            write(io_lun,'(4x,"Force tolerance:    ",f20.10)') MDcgtol
            write(io_lun,'(4x,"Maximum stress         ",3f10.6)') &
              max_stress
            write(io_lun,'(4x,"Stress tolerance:   ",f20.10)') &
              cell_stress_tol
            write(io_lun,'(4x,"Enthalpy change:    ",f20.10," ",a2)') &
              en_conv*dH, en_units(energy_units)
            write(io_lun,'(4x,"Enthalpy tolerance: ",f20.10)') &
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
      stressx = -stress(1,1)/volume
      stressy = -stress(2,2)/volume
      stressz = -stress(3,3)/volume
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
      search_dir_x = gamma*search_dir_x + stressx - press
      search_dir_y = gamma*search_dir_y + stressy - press
      search_dir_z = gamma*search_dir_z + stressz - press
      ! end if

      new_rcellx = rcellx
      new_rcelly = rcelly
      new_rcellz = rcellz

      ! Minimise in this direction. Constraint information is also used within
      ! safemin_cell. Look in move_atoms.module.f90 for further information.
      call safemin_cell(new_rcellx, new_rcelly, new_rcellz, search_dir_x, &
                        search_dir_y, search_dir_z, search_dir_mean, press, &
                        enthalpy0, enthalpy1, fixed_potential, vary_mu)
      ! Output positions to UpdatedAtoms.dat
      if (inode==ionode .and. iprint_gen > 1) then
        write(io_lun,'(4x,a4,a15)') "Atom", "Position"
        do i = 1, ni_in_cell
          write(io_lun,'(4x,i8,3f15.8)') i,atom_coord(:,i)
        end do
      end if
      call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))

      ! Analyse Stresses and energies
      dH = enthalpy0 - enthalpy1
      newRMSstress = sqrt(((stress(1,1)*stress(1,1)) + &
                           (stress(2,2)*stress(2,2)) + &
                           (stress(3,3)*stress(3,3)))/three)
      dRMSstress = RMSstress - newRMSstress

      volume = rcellx*rcelly*rcellx
      max_stress = zero
      do i=1,3
        stress_diff = abs(press + stress(i,i))/volume
        if (stress_diff > max_stress) max_stress = stress_diff
      end do

      if (inode==ionode) then
        write(io_lun,'(2x,"GeomOpt + Iter: ",i4," MaxF: ",f12.8," H: " e16.8," dH: ",f12.8)') &
             iter, max, enthalpy1, en_conv*dH
        if (iprint_MD > 1) then
          write(io_lun,'(4x,"Force Residual:     ",f20.10," ",a2,"/",a2)') &
            for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
            d_units(dist_units)
          write(io_lun,'(4x,"Maximum force:         ",f20.10)') max
          write(io_lun,'(4x,"Force tolerance:       ",f20.10)') MDcgtol
          write(io_lun,'(4x,"Maximum stress         ",3f10.6)') &
            max_stress
          write(io_lun,'(4x,"Stress tolerance: ",f20.10)') &
            cell_stress_tol
          write(io_lun,'(4x,"Enthalpy change:       ",f20.10," ",a2)') &
            en_conv*dH, en_units(energy_units)
          write(io_lun,'(4x,"Enthalpy tolerance:    ",f20.10)') &
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
        if (abs(max) < MDcgtol .and. max_stress < cell_stress_tol) then
          done_cell = .true.
          if (inode==ionode) then
            write(io_lun,'(2x,a,i4,a,i4,a)') "GeomOpt converged in ", &
              iter-iter_cell," ionic steps and ", iter_cell, " cell steps"
            write(io_lun, fmt='(4x,"Maximum force:      ",f20.10)') max
            write(io_lun, fmt='(4x,"Force tolerance:    ",f20.10)') MDcgtol
            write(io_lun, fmt='(4x,"Enthalpy change:    ",f20.10)') dH
            write(io_lun, fmt='(4x,"Enthalpy tolerance: ",f20.10)') &
              enthalpy_tolerance
            write(io_lun, fmt='(4x,"Maximum stress:     ",f20.10)') &
              max_stress
            write(io_lun, fmt='(4x,"Stress tolerance:   ",f20.10)') &
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

  end subroutine full_cg_run_double_loop

  !!****f* control/full_cg_run_single_vector *
  !!
  !!  NAME
  !!   full_cg_run_single_vector
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Conjugate gradients opnimisation of cell vectors and ionic positions 
  !!  by minimisation of force vector defined in Pfrommer et al.
  !!  J. Comput. Phys. 131, 233 (1997)
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
                             cell_stress_tol
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin_full, cq_to_vector, enthalpy, &
                             enthalpy_tolerance
    use GenComms,      only: inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: tot_force, stress
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop, write_xsf
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  ONLY: dump_InfoMatGlobal, dump_pos_and_matrices
    use md_control,    only: flag_write_xsf, target_pressure, fac_HaBohr32GPa

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma, gg1, &
                      press, rcellx_ref, rcelly_ref, rcellz_ref, &
                      enthalpy0, enthalpy1, dH, volume, max_stress, &
                      stress_diff
    integer        :: i,j,k,iter,length, jj, stat
    logical        :: done
    type(cq_timer) :: tmr_l_iter
    real(double), allocatable, dimension(:,:) :: cg, force, force_old, &
                                                 config
    real(double), dimension(3) :: one_plus_strain, strain, cell_ref

    if (inode==ionode .and. iprint_MD > 2) &
      write(io_lun,'(2x,a)') "control/full_cg_run_single_vector"

    allocate(cg(3,ni_in_cell+1), STAT=stat)
    allocate(force(3,ni_in_cell+1), STAT=stat)
    allocate(force_old(3,ni_in_cell+1), STAT=stat)
    allocate(config(3,ni_in_cell+1), STAT=stat)
    length = 3*(ni_in_cell+1)
    call reg_alloc_mem(area_general, 5*length, type_dbl)

    if (inode==ionode) &
      write (io_lun, fmt='(/4x,"Starting CG cell and atomic relaxation"/)')
    cg = zero
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3*(ni_in_cell+1)
    if (inode==ionode .and. iprint_gen > 0) then
      write(io_lun,'(4x,"Welcome to cg_run. Doing ",i4," steps")') MDn_steps
      write(io_lun,'(4x,"Force tolerance:    ",f20.10)') MDcgtol
      write(io_lun,'(4x,"Enthalpy tolerance: ",f20.10)') enthalpy_tolerance
    end if
    press = target_pressure/fac_HaBohr32GPa
    energy0 = total_energy
    enthalpy0 = enthalpy(energy0, press)
    dH = zero

    ! reference cell to compute strain
    cell_ref(1) = rcellx
    cell_ref(2) = rcelly
    cell_ref(3) = rcellz

    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.)
    call dump_pos_and_matrices
    call get_maxf(max)
    enthalpy0 = enthalpy(energy0, press)
    if (inode==ionode) then
      write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," H: "e16.8," dH: ",f12.8)') &
           0, max, enthalpy0, zero
    end if

    iter = 1
    ggold = zero
    force_old = zero
    energy1 = energy0
    enthalpy1 = enthalpy0
    do while (.not. done)
      call start_timer(tmr_l_iter, WITH_LEVEL) ! Construct ratio for conjugacy
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
        gamma = (gg-gg1)/ggold ! PR - change to gg/ggold for FR
      end if
      if(gamma<zero) gamma = zero
      if (inode == ionode .and. iprint_MD > 2) &
        write(io_lun,*) ' CHECK :: Force Residual = ', &
                        for_conv * sqrt(gg)/ni_in_cell
      if (inode == ionode .and. iprint_MD > 2) &
        write(io_lun,*) ' CHECK :: gamma = ', gamma
      if (CGreset) then
        if (gamma > one) then
          if (inode==ionode) &
            write(io_lun,*) ' CG direction is reset! '
          gamma = zero
      end if
      end if
      if (inode == ionode) &
        write (io_lun,'(4x,"CG iteration: ",i5," Gamma: ",f12.6)') iter, gamma
      ggold = gg
      ! Build search direction
      do j=1,ni_in_cell+1
        do i=1,3
          cg(i,j) = gamma*cg(i,j) + force(i,j)
        end do
      end do
      force_old = force
      ! Minimise in this direction
      call safemin_full(config, cg, cell_ref, enthalpy0, enthalpy1, &
                        press, fixed_potential, vary_mu, enthalpy1)

      ! Output positions
      if (inode==ionode .and. iprint_gen > 1) then
        write(io_lun,'(4x,a4,a15)') "Atom", "Position"
        do i = 1, ni_in_cell
          write(io_lun,'(4x,i8,3f15.8)') i,atom_coord(:,i)
        end do
      end if
      call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
      if (flag_write_xsf) call write_xsf('trajectory.xsf', iter)

      ! Analyse forces and stress
      g0 = dot(length, force, 1, force, 1)
      call get_maxf(max)
      volume = rcellx*rcelly*rcellx
      max_stress = zero
      do i=1,3
        stress_diff = abs(press + stress(i,i))/volume
        if (stress_diff > max_stress) max_stress = stress_diff
      end do

      ! Output and energy changes
      dH = enthalpy0 - enthalpy1
      if (inode==ionode) then
        write(io_lun,'(2x,"GeomOpt - Iter: ",i4," MaxF: ",f12.8," H: "e16.8," dH: ",f12.8)') &
             iter, max, enthalpy1, en_conv*dH
        if (iprint_MD > 1) then
          write(io_lun,'(4x,"Force Residual:     ",f20.10," ",a2,"/",a2)') &
            for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), & 
            d_units(dist_units)
          write(io_lun,'(4x,"Maximum force:      ",f20.10)') max
          write(io_lun,'(4x,"Force tolerance:    ",f20.10)') MDcgtol
          write(io_lun,'(4x,"Maximum stress         ",3f10.6)') &
            max_stress
          write(io_lun, fmt='(4x,"Stress tolerance:   ",f20.10)') &
            cell_stress_tol
          write(io_lun,'(4x,"Enthalpy change:    ",f20.10," ",a2)') &
            en_conv*dH, en_units(energy_units)
          write(io_lun,'(4x,"Enthalpy tolerance: ",f20.10)') &
            enthalpy_tolerance
        end if
      end if

      enthalpy0 = enthalpy1
      if (iter > MDn_steps) then
        done = .true.
        if (inode==ionode) &
          write(io_lun, fmt='(4x,"Exceeded number of MD steps: ",i6)') iter
      end if
      if (abs(dH) < enthalpy_tolerance) then
        if (abs(max) < MDcgtol .and. max_stress < cell_stress_tol) then
          done = .true.
          if (inode==ionode) then
            write(io_lun,'(2x,a,i4,a)') "GeomOpt converged in ", iter, " iterations"
            write(io_lun, fmt='(4x,"Maximum force:      ",f20.10)') max
            write(io_lun, fmt='(4x,"Force tolerance:    ",f20.10)') MDcgtol
            write(io_lun, fmt='(4x,"Enthalpy change:    ",f20.10)') dH
            write(io_lun, fmt='(4x,"Enthalpy tolerance: ",f20.10)') &
              enthalpy_tolerance
            write(io_lun, fmt='(4x,"Maximum stress:     ",f20.10)') &
              max_stress
            write(io_lun, fmt='(4x,"Stress tolerance:   ",f20.10)') &
              cell_stress_tol
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
