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
!!  SOURCE
!!
  subroutine control_run(fixed_potential, vary_mu, total_energy)

    use datatypes
    use dimens,               only: r_core_squared, r_h
    use GenComms,             only: my_barrier, cq_abort
    use pseudopotential_data, only: set_pseudopotential
    use force_module,         only: tot_force
    use minimise,             only: get_E_and_F
    use global_module,        only: runtype, flag_self_consistent, flag_out_wf, &
                                    flag_write_DOS, flag_opt_cell
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
            call cell_cg_run(fixed_potential, vary_mu, total_energy)
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

    call dump_pos_and_matrices
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
                             IPRINT_TIME_THRES1, flag_MDold
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin, safemin2
    use GenComms,      only: gsum, myid, inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: tot_force
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use store_matrix,  ONLY: dump_InfoMatGlobal, dump_pos_and_matrices

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy
    
    ! Local variables
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma
    integer        :: i,j,k,iter,length, jj, lun, stat
    logical        :: done
    type(cq_timer) :: tmr_l_iter
    real(double), allocatable, dimension(:,:) :: cg
    real(double), allocatable, dimension(:)   :: x_new_pos, y_new_pos,&
                                                 z_new_pos

    allocate(cg(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating cg in control: ",&
                       ni_in_cell, stat)
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
    if (.NOT. flag_MDold) then
      call dump_InfoMatGlobal()
    endif
    iter = 1
    ggold = zero
    energy1 = energy0
    do while (.not. done)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       ! Construct ratio for conjugacy
       gg = zero
       do j = 1, ni_in_cell
          gg = gg +                              &
               tot_force(1,j) * tot_force(1,j) + &
               tot_force(2,j) * tot_force(2,j) + &
               tot_force(3,j) * tot_force(3,j)
       end do
       if (abs(ggold) < 1.0e-6_double) then
          gamma = zero
       else
          gamma = gg/ggold
       end if
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
       ! Minimise in this direction
       !ORI call safemin(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
       !ORI              energy1, fixed_potential, vary_mu, energy1)
       if (.NOT. flag_MDold) then
         call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
                      energy1, fixed_potential, vary_mu, energy1)
       else
         call safemin(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
                      energy1, fixed_potential, vary_mu, energy1)
       endif
       ! Output positions
       if (myid == 0 .and. iprint_gen > 1) then
          do i = 1, ni_in_cell
             write (io_lun, 1) i, atom_coord(1,i), atom_coord(2,i), &
                               atom_coord(3,i)
          end do
       end if
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       ! Analyse forces
       g0 = dot(length, tot_force, 1, tot_force, 1)
       max = zero
       do i = 1, ni_in_cell
          do k = 1, 3
             if (abs(tot_force(k,i)) > max) max = abs(tot_force(k,i))
          end do
       end do
       ! Output and energy changes
       iter = iter + 1
       dE = energy0 - energy1
       !if(myid==0) write(io_lun,6) for_conv*max, en_units(energy_units), d_units(dist_units)
       if (myid == 0) write (io_lun, 4) en_conv*dE, en_units(energy_units)
       if (myid == 0) write (io_lun, 5) for_conv*sqrt(g0/ni_in_cell), &
                                        en_units(energy_units),       &
                                        d_units(dist_units)
       energy0 = energy1
       !energy1 = abs(dE)
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
!!  SOURCE
!!
  subroutine md_run (fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use global_module,  only: iprint_gen, ni_in_cell, x_atom_cell,    &
                              y_atom_cell, z_atom_cell, area_general, &
                              flag_read_velocity, flag_quench_MD,     &
                              temp_ion, flag_MDcontinue, MDinit_step, &
                              flag_MDold,n_proc_old,glob2node_old,    &
                              flag_LmatrixReuse,flag_XLBOMD,          &
                              flag_dissipation,flag_FixCOM,           &
                              flag_fire_qMD, flag_diagonalisation,    &
                              nspin, flag_Multisite, flag_SFcoeffReuse
    use group_module,   only: parts
    use primary_module, only: bundle
    use minimise,       only: get_E_and_F
    use move_atoms,     only: velocityVerlet, updateIndices,           &
                              init_velocity,update_H,check_move_atoms, &
                              update_atom_coord,wrap_xyz_atom_cell,    &
                              zero_COM_velocity, calculate_kinetic_energy
    use GenComms,       only: gsum, myid, my_barrier, inode, ionode,  &
                              gcopy
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_positions, read_velocity,         &
                              write_velocity, read_fire
    use io_module,      only: write_atomic_positions, pdb_template,   &
                              check_stop
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use move_atoms,     only: fac_Kelvin2Hartree, update_pos_and_matrices, updateL, updateLorK, updateSFcoeff
    !use io_module2,     ONLY: InfoL
    use store_matrix,   ONLY: dump_InfoMatGlobal,grab_InfoMatGlobal, &
                    matrix_store_global, grab_matrix2, InfoMatrixFile, dump_pos_and_matrices
    !use DiagModule,     ONLY: diagon
    use mult_module,    ONLY: matL,L_trans
    use matrix_data,    ONLY: Lrange
    use UpdateInfo_module, ONLY: Matrix_CommRebuild
    use XLBOMD_module,  ONLY: Ready_XLBOMD, Do_XLBOMD
    use Integrators,    ONLY: vVerlet_v_dthalf,vVerlet_r_dt, fire_qMD
    use constraint_module, ONLY: correct_atomic_position,correct_atomic_velocity, &
         ready_constraint,flag_RigidBonds
    use input_module,   ONLY: io_assign, io_close

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical      :: vary_mu, fixed_potential
    real(double) :: total_energy
    
    ! Local variables
    real(double), allocatable, dimension(:,:) :: velocity
    integer       ::  iter, i, k, length, stat, i_first, i_last, &
         nfile, symm
    integer       :: lun, ios, md_steps ! SA 150204; 150213 md_steps: counter for MD steps
    real(double)  :: temp, KE, energy1, energy0, dE, max, g0
    real(double)  :: energy_md
    character(50) :: file_velocity='velocity.dat'
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
    
    allocate(velocity(3,ni_in_cell), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating velocity in md_run: ", &
                       ni_in_cell, stat)
    call reg_alloc_mem(area_general, 3*ni_in_cell, type_dbl)
    velocity = zero
    if (flag_read_velocity) then
       call read_velocity(velocity, file_velocity)
    else
       if(temp_ion > RD_ERR) then
          if(inode == ionode) &
               call init_velocity(ni_in_cell, temp_ion, velocity)
          call gcopy(velocity, 3, ni_in_cell)
       else
          velocity = zero
       end if
    end if
    energy0 = zero
    energy1 = zero
    dE = zero
    length = 3*ni_in_cell
    if (myid == 0 .and. iprint_gen > 0) write(io_lun, 2) MDn_steps
    ! Find energy and forces
    if (flag_fire_qMD) then
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .true.)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy0, .true., .false.)
    end if

    ! XL-BOMD
    if (flag_XLBOMD .AND. flag_dissipation .AND. .NOT.flag_MDold) &
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
    ! Dump global data
    if (.NOT. flag_MDold) then
      !call dump_InfoMatGlobal(index=0,MDstep=i_first)
      call dump_pos_and_matrices(index=0,MDstep=i_first,velocity=velocity)
    endif

    energy_md = energy0

    if (flag_fire_qMD) then
       step_qMD = i_first ! SA 20150201
       done = .false.     ! SA 20150201
       fire_step_max = 10.0_double * MDtimestep
       ! reading FIRE parameters of the last run
       call read_fire(fire_N, fire_N2, fire_P0, MDtimestep, fire_alpha)
    endif

    !ORI do iter = 1, MDn_steps
    do iter = i_first, i_last
       if (myid == 0) &
            write (io_lun, fmt='(4x,"MD run, iteration ",i5)') iter
       ! Fetch old relations
       if (.NOT. flag_MDold) then
         if (.NOT. allocated(glob2node_old)) then
           allocate (glob2node_old(ni_in_cell), STAT=stat)
           if (stat.NE.0) call cq_abort('Error allocating glob2node_old: ', ni_in_cell)
         endif

        !commented out Nov. 12, 2017, TM
        ! call grab_InfoMatGlobal(InfoGlob,0)
        ! n_proc_old = InfoGlob%numprocs
        ! glob2node_old(:) = InfoGlob%glob_to_node(:)
         
       endif

  !! For Debuggging !!
  !     call dump_pos_and_matrices(index=1,MDstep=i_first)
  !! For Debuggging !!

       !%%! Evolve atoms - either FIRE (quenched MD) or velocity Verlet
       if (flag_fire_qMD) then
          call fire_qMD(fire_step_max,MDtimestep,velocity,tot_force,flag_movable,iter,&
                        fire_N,fire_N2,fire_P0,fire_alpha) ! SA 20150204
       else
          call vVerlet_v_dthalf(MDtimestep,velocity,tot_force,flag_movable)
          call vVerlet_r_dt(MDtimestep,velocity,flag_movable)
       end if
       ! Constrain position
       if (flag_RigidBonds) call correct_atomic_position(velocity,MDtimestep)
       ! Reset-up
       if(flag_SFcoeffReuse) then
        call update_pos_and_matrices(updateSFcoeff,velocity)
       else
        call update_pos_and_matrices(updateLorK,velocity)
       endif

       if (flag_XLBOMD) call Do_XLBOMD(iter,MDtimestep)
       call update_H(fixed_potential)

  !2018.Jan.4 TM 
  !   We need to update flag_movable, since the order of x_atom_cell (or id_glob) 
  !   may change after the atomic positions are updated.
  !
       call check_move_atoms(flag_movable)

  !! For Debuggging !!
  !     call dump_pos_and_matrices(index=2,MDstep=i_first)
  !! For Debuggging !!

       if (flag_fire_qMD) then
          call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .true.,iter)
          call dump_pos_and_matrices(index=0,MDstep=iter,velocity=velocity)
       else
          call get_E_and_F(fixed_potential, vary_mu, energy1, .true., .false.,iter)
          call dump_pos_and_matrices(index=0,MDstep=iter,velocity=velocity)
          call vVerlet_v_dthalf(MDtimestep,velocity,tot_force,flag_movable,second_call)
       end if

  !! For Debuggging !!
       !!call cq_abort(" STOP FOR DEBUGGING")
  !! For Debuggging !!

       ! Constrain velocity
       if (flag_RigidBonds) call correct_atomic_velocity(velocity)
       !%%! END of Evolve atoms

       ! Let's analyse
       if (flag_FixCOM) call zero_COM_velocity(velocity)
       call calculate_kinetic_energy(velocity,KE)

       ! Print out energy
       energy_md = energy1
       if (myid == 0) &
         write (io_lun, fmt='(4x,"Kinetic Energy in K     : ",f15.8)') &
                KE / (three / two * ni_in_cell) / fac_Kelvin2Hartree
       if (myid == 0) write (io_lun, 8) iter, KE, energy_md, KE+energy_md
       ! Output positions
       if (myid == 0 .and. iprint_gen > 1) then
         do i = 1, ni_in_cell
           write (io_lun, 1) i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
         enddo
       endif

       ! Dump global data
       !if (.NOT. flag_MDold) then
       !  call dump_InfoMatGlobal(index=0,MDstep=iter)
       !endif
       
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
       if(myid==0) then
          write (io_lun, 6) max
          write (io_lun, 4) dE
          write (io_lun, 5) sqrt(g0/ni_in_cell)
       end if
       energy0 = energy1
       energy1 = abs(dE)
       if (myid == 0 .and. mod(iter, MDfreq) == 0) &
            call write_positions(iter, parts)
       call my_barrier
       !to check IO of velocity files
       call write_velocity(velocity, file_velocity)
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       !to check IO of velocity files

       call check_stop(done, iter)
       if (flag_fire_qMD) then
          if (abs(max) < MDcgtol) then
             if ((iter - step_qMD) > 1) then
                n_stop_qMD = 0
             end if
             n_stop_qMD = n_stop_qMD + 1
             step_qMD = iter
             if (myid == 0) then
                write (io_lun, fmt='(4x,i4,4x,"Maximum force below threshold: ",f12.6)') n_stop_qMD, max
             end if
             if (n_stop_qMD > 2) then
                done = .true.
             end if
          end if
       end if

       if (done) exit

!%%!   if (myid == 0) &
!%%!        write (io_lun, fmt='(4x,"MD run, iteration ",i5)') iter
!%%!   ! Fetch old relations
!%%!   if (.NOT. flag_MDold) then
!%%!     if (.NOT. allocated(glob2node_old)) then
!%%!       allocate (glob2node_old(ni_in_cell), STAT=stat)
!%%!       if (stat.NE.0) call cq_abort('Error allocating glob2node_old: ', ni_in_cell)
!%%!     endif
!%%!     if (inode.EQ.ionode) call grab_InfoGlobal(n_proc_old,glob2node_old)
!%%!     call gcopy(n_proc_old)
!%%!     call gcopy(glob2node_old,ni_in_cell)
!%%!   endif

!%%!   call velocityVerlet(fixed_potential, bundle, MDtimestep, temp, &
!%%!                       KE, flag_quench_MD, velocity, tot_force, iter)
!%%!   if (myid == 0) &
!%%!        write (io_lun, fmt='(4x,"Kinetic Energy in K     : ",f15.8)') &
!%%!              KE / (three / two * ni_in_cell) / fac_Kelvin2Hartree
!%%!   if (myid == 0) write (io_lun, 8) iter, KE, energy_md, KE+energy_md
!%%!   ! Output positions
!%%!   if (myid == 0 .and. iprint_gen > 1) then
!%%!      do i = 1, ni_in_cell
!%%!         write (io_lun, 1) i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
!%%!      end do
!%%!   end if
!%%!   !Now, updateIndices and update_atom_coord are done in velocityVerlet 
!%%!   !call updateIndices(.false.,fixed_potential, number_of_bands) 
!%%!   ! L-matrix reconstruction (used to be called at updateIndices3)
!%%!   if (.NOT.flag_MDold .AND. &
!%%!       .NOT.diagon     .AND. &
!%%!       flag_LmatrixReuse       ) then
!%%!     call grab_matrix2('L',inode,nfile,InfoL)
!%%!     call my_barrier()
!%%!     call Matrix_CommRebuild(InfoL,Lrange,L_trans,matL(1),nfile,symm)
!%%!   endif
!%%!   ! For XL-BOMD
!%%!   if (flag_XLBOMD .AND. .NOT.diagon) call Do_XLBOMD(iter,MDtimestep)
!%%!   ! Updates hamiltonian (used to be called at updateIndices3)
!%%!   if (.NOT.flag_MDold) call update_H(fixed_potential)
!%%!   !ORI call get_E_and_F(fixed_potential, vary_mu, energy1, .true., &
!%%!   !ORI                  .false.)
!%%!   call get_E_and_F(fixed_potential, vary_mu, energy1, .true., &
!%%!                    .false., iter)
!%%!    energy_md = energy1
!%%!   ! Dump global data
!%%!   if (.NOT. flag_MDold) then
!%%!     if (inode.EQ.ionode) call dump_InfoGlobal(iter)
!%%!   endif

!%%!   ! Analyse forces
!%%!   g0 = dot(length, tot_force, 1, tot_force, 1)
!%%!   max = zero
!%%!   do i = 1, ni_in_cell
!%%!      do k = 1, 3
!%%!         if (abs(tot_force(k,i)) > max) max = tot_force(k,i)
!%%!      end do
!%%!   end do
!%%!   ! Output and energy changes
!%%!   dE = energy0 - energy1
!%%!   if (myid == 0) write (io_lun, 6) max
!%%!   if (myid == 0) write (io_lun, 4) dE
!%%!   if (myid == 0) write (io_lun, 5) sqrt(g0/ni_in_cell)
!%%!   energy0 = energy1
!%%!   energy1 = abs(dE)
!%%!   if (myid == 0 .and. mod(iter, MDfreq) == 0) &
!%%!        call write_positions(iter, parts)
!%%!   call my_barrier
!%%!   !to check IO of velocity files
!%%!   call write_velocity(velocity, file_velocity)
!%%!   call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
!%%!   !to check IO of velocity files
!%%!   call check_stop(done, iter)
!%%!   if (done) exit
    end do
    deallocate(velocity, STAT=stat)
    if (stat /= 0) call cq_abort("Error deallocating velocity in md_run: ", &
                                 ni_in_cell, stat)
    call reg_dealloc_mem(area_general, 3 * ni_in_cell, type_dbl)
    return
1   format(4x,'Atom ',i4,' Position ',3f15.8)
2   format(4x,'Welcome to md_run. Doing ',i4,' steps')
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Energy change           : ',f15.8)
5   format(4x,'Force Residual          : ',f15.8)
6   format(4x,'Maximum force component : ',f15.8)
8   format(4x,'*** MD step ',i4,' KE: ',f18.8,&
           ' IntEnergy',f20.8,' TotalEnergy',f20.8)
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
  !!  SOURCE
  !!
  subroutine pulay_relax(fixed_potential, vary_mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module,  only: iprint_MD, ni_in_cell, x_atom_cell,   &
                              y_atom_cell, z_atom_cell, id_glob,    &
                              atom_coord, area_general, flag_pulay_simpleStep, &
                              glob2node_old, n_proc_old, flag_MDold, &
                              flag_diagonalisation, nspin, flag_LmatrixReuse, &
                              flag_SFcoeffReuse
    use group_module,   only: parts
    use minimise,       only: get_E_and_F
    use move_atoms,     only: pulayStep, velocityVerlet,            &
                              updateIndices, updateIndices3, update_atom_coord,     &
                              safemin, safemin2, update_H, update_pos_and_matrices
    use move_atoms,     only: updateL, updateLorK, updateSFcoeff
    use GenComms,       only: gsum, myid, inode, ionode, gcopy, my_barrier
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_atomic_positions, pdb_template, &
                              check_stop
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use primary_module, only: bundle
    use mult_module, ONLY: matK, S_trans, matrix_scale, matL, L_trans
    use matrix_data, ONLY: Hrange, Lrange
    use UpdateInfo_module, ONLY: Matrix_CommRebuild
    use store_matrix,   only: dump_pos_and_matrices

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
    if(.NOT.flag_MDold) call dump_pos_and_matrices
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
          if(.NOT.flag_MDold) then
           if(flag_SFcoeffReuse) then
            call update_pos_and_matrices(updateSFcoeff,cg)
           else
            call update_pos_and_matrices(updateLorK,cg)
           endif
          else
             call update_atom_coord
             call updateIndices(.true., fixed_potential)
          end if
          call update_H(fixed_potential)
          call get_E_and_F(fixed_potential, vary_mu, energy1, .true., &
                           .false.)
          if(.NOT.flag_MDold) call dump_pos_and_matrices
       else
          if (.NOT. flag_MDold) then
             call safemin2(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
                  energy1, fixed_potential, vary_mu, energy1)
          else
             call safemin(x_new_pos, y_new_pos, z_new_pos, cg, energy0, &
                  energy1, fixed_potential, vary_mu, energy1)
          end if
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
       if(.NOT.flag_MDold) then
          if(flag_SFcoeffReuse) then
           call update_pos_and_matrices(updateSFcoeff,cg)
          else
           call update_pos_and_matrices(updateLorK,cg)
          endif
       else
          call update_atom_coord
          call updateIndices(.true., fixed_potential)
       end if
       call update_H(fixed_potential)
       call get_E_and_F(fixed_potential, vary_mu, energy1, .true., &
                        .false.)
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))
       if(.NOT.flag_MDold) call dump_pos_and_matrices

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
                             IPRINT_TIME_THRES1, flag_MDold, cell_en_tol, &
                             cell_constraint_flag
    use group_module,  only: parts
    use minimise,      only: get_E_and_F
    use move_atoms,    only: safemin_cell
    use GenComms,      only: gsum, myid, inode, ionode
    use GenBlas,       only: dot
    use force_module,  only: stress, tot_force
    use io_module,     only: write_atomic_positions, pdb_template, &
                             check_stop
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use timer_module
    use io_module,      only: leqi
    use dimens, ONLY: r_super_x, r_super_y, r_super_z

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential
    real(double) :: total_energy

    ! Local variables
    real(double)   :: energy0, energy1, max, g0, dE, gg, ggold, gamma
    integer        :: i,j,k,iter,length, jj, lun, stat, reset_iter
    logical        :: done
    type(cq_timer) :: tmr_l_iter
    real(double) :: new_rcellx, new_rcelly, new_rcellz, search_dir_x, search_dir_y,&
                    search_dir_z, stressx, stressy, stressz, RMSstress, newRMSstress,&
                    dRMSstress, search_dir_mean, mean_stress

    call reg_alloc_mem(area_general, 6 * ni_in_cell, type_dbl)

    if (myid == 0 .and. iprint_gen > 0) &
         write (io_lun, fmt='(/4x,"Starting CG lattice vector relaxation"/)')
    search_dir_z = zero
    search_dir_x = zero
    search_dir_y = zero
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
    do while (.not. done)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       stressx = -stress(1)
       stressy = -stress(2)
       stressz = -stress(3)
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
         search_dir_x = gamma*search_dir_x + stressx
         search_dir_y = gamma*search_dir_y + stressy
         search_dir_z = gamma*search_dir_z + stressz
       end if

       if (inode == ionode .and. iprint_gen > 0) &
           write(io_lun,*)  "Initial cell dims", rcellx, rcelly, rcellz

       new_rcellx = rcellx
       new_rcelly = rcelly
       new_rcellz = rcellz

       ! Minimise in this direction. Constraint information is also used within
       ! safemin_cell. Look in move_atoms.module.f90 for further information.
       call safemin_cell(new_rcellx, new_rcelly, new_rcellz, search_dir_x, &
                         search_dir_y, search_dir_z, energy0, energy1, &
                         fixed_potential, vary_mu, energy1, search_dir_mean)

       ! Output positions to UpdatedAtoms.dat
       if (myid == 0 .and. iprint_gen > 1) then
          do i = 1, ni_in_cell
             write (io_lun, 1) i, atom_coord(1,i), atom_coord(2,i), &
                               atom_coord(3,i)
          end do
       end if
       call write_atomic_positions("UpdatedAtoms.dat", trim(pdb_template))

       ! Analyse Stresses and energies
       dE = energy0 - energy1
       newRMSstress = sqrt(((stress(1)*stress(1)) + (stress(2)*stress(2)) + (stress(3)*stress(3)))/3)
       dRMSstress = RMSstress - newRMSstress

       iter = iter + 1
       reset_iter = reset_iter +1

       if (myid == 0 .and. iprint_gen > 0) then
           write (io_lun, 4) en_conv*dE, en_units(energy_units)
           write (io_lun, 5) dRMSstress, "Ha"
       end if

       energy0 = energy1

       ! Check exit criteria

       ! First exit is if too many steps have been taken. Default is 50.
       if (iter > MDn_steps) then
          done = .true.
          if (myid == 0) &
               write (io_lun, fmt='(4x,"Exceeded number of MD steps: ",i4)') &
                     iter
       end if

       ! Second exit is if the desired energy tolerence has ben reached
       ! Will replace with stress tolerance when more reliable
       if (abs(dE)<cell_en_tol) then
          done = .true.
          if (myid == 0 .and. iprint_gen > 0) &
               write (io_lun, fmt='(4x,"Energy change below threshold: ",f20.10)') &
                     max
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
4   format(4x,'Energy change: ',f15.8,' ',a2)
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

! SUB FLASH_POS_AND_MATRICES has been moved to "store_matrix"
!  subroutine dump_pos_and_matrices
!    use global_module, ONLY: nspin, nspin_SF, flag_diagonalisation, flag_Multisite
!    use matrix_data, ONLY: Lrange, Hrange, SFcoeff_range, SFcoeffTr_range, HTr_range
!    use mult_module, ONLY: matL,L_trans, matK, matSFcoeff
!    use store_matrix, ONLY:dump_matrix_update
!    use io_module, ONLY: append_coords, write_atomic_positions, pdb_template
!
!    implicit none
!    integer:: both=0 , mat=1
!    logical :: append_coords_bkup
!
!    !!! Check whether we should write out the files or not.  !!!
!     !   1. check elapsed time 
!     !   2. check some specific file
!
!    !!! Write out Files to restart..   
!     !   InfoGlob.dat
!     !     PAO coefficients of supports  or Blips
!     !     L matrix of K matrix
!     !     (for XL-BOMD, files in previous steps should be also printed out)
!     !     coodinates file ?
!
!       if(flag_Multisite) then
!        call dump_matrix_update('SFcoeff',matSFcoeff(1),SFcoeff_range,index_in=0,iprint_mode=mat)
!        if(nspin_SF .eq. 2) call dump_matrix_update('SFcoeff2',matSFcoeff(2),SFcoeff_range,index_in=0,iprint_mode=mat)
!       endif
!
!       if(flag_diagonalisation) then
!        call dump_matrix_update('K',matK(1),Hrange,index_in=0,iprint_mode=both)
!        if(nspin .eq. 2) call dump_matrix_update('K2',matK(2),Hrange,index_in=0,iprint_mode=both)
!       else
!        call dump_matrix_update('L',matL(1),Lrange,index_in=0,iprint_mode=both)
!        if(nspin .eq. 2) call dump_matrix_update('L2',matL(2),Lrange,index_in=0,iprint_mode=both)
!       endif
!
!       append_coords_bkup = append_coords; append_coords = .false.
!        call write_atomic_positions('coord_next.dat',trim(pdb_template))
!       append_coords = append_coords_bkup
!
!   return
!  end subroutine dump_pos_and_matrices

end module control
