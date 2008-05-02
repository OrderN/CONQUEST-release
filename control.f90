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
!!  SOURCE
!!
module control

  use datatypes
  use global_module, ONLY: io_lun
  use GenComms, ONLY: cq_abort
  
  implicit none

  integer :: MDn_steps 
  integer :: MDfreq 
  real(double) :: MDtimestep 
  real(double) :: MDcgtol 
  logical :: CGreset

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
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
!!    Indented, ROBODoc, stripped get_E_and_F calls
!!   11/06/2001 dave
!!    Added RCS Id and Log tags and GenComms
!!   13/06/2001 dave
!!    Adapted to use force_module
!!   17/06/2002 dave
!!    Improved headers slightly
!!   16:42, 2003/02/03 dave
!!    Changed to be control_run as part of control module
!!   17:31, 2003/02/05 dave
!!    Removed old, silly MD and replaced with call to either static (i.e. just ground state) or CG
!!   10:06, 12/03/2003 drb 
!!    Added MD
!!   2004/11/10, drb
!!    Removed common use
!!   2008/02/13 08:23 dave
!!    Added Pulay optimiser for atomic positions
!!  SOURCE
!!
  subroutine control_run(fixed_potential, vary_mu, number_of_bands, mu, total_energy)

    use datatypes
    use dimens, ONLY: r_core_squared,r_h
    use GenComms, ONLY: my_barrier, cq_abort
    use ewald_module, ONLY: ewald
    use pseudopotential_data, ONLY: set_pseudopotential
    use force_module, ONLY: tot_force
    use minimise, ONLY: get_E_and_F
    use global_module, ONLY: runtype

    implicit none

    ! Shared variables
    logical :: vary_mu, fixed_potential

    real(double) :: number_of_bands, mu
    real(double) :: total_energy

    ! Local variables
    logical :: NoMD
    logical, external :: leqi

    integer :: i,j

    real(double) :: spr, e_r_OLD

    if(leqi(runtype,'static')) then
       call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, total_energy, .true., .true.)
       return
    else if(leqi(runtype,'cg')) then
       call cg_run(fixed_potential, vary_mu, number_of_bands, mu, total_energy)
    else if(leqi(runtype,'md')) then
       call md_run(fixed_potential, vary_mu, number_of_bands, mu, total_energy)
    else if(leqi(runtype,'pulay')) then
       call pulay_relax(fixed_potential, vary_mu, number_of_bands, mu, total_energy)
    else if(leqi(runtype,'dummy')) then
       call dummy_run(fixed_potential, vary_mu, number_of_bands, mu, total_energy)
    else
       call cq_abort('control: Runtype not specified !')
    end if
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
!!  SOURCE
!!
  subroutine cg_run(fixed_potential, vary_mu, number_of_bands, mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, ONLY: iprint_gen, ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell, id_glob, &
         atom_coord, rcellx, rcelly, rcellz, area_general, iprint_MD
    use group_module, ONLY: parts
    use minimise, ONLY: get_E_and_F
    use move_atoms, ONLY: safemin
    use GenComms, ONLY: gsum, myid, inode, ionode
    use GenBlas, ONLY: dot
    use force_module, ONLY: tot_force
    use io_module, ONLY: write_atomic_positions, pdb_template
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential

    real(double) :: number_of_bands, mu
    real(double) :: total_energy
    
    ! Local variables
    real(double), allocatable, dimension(:,:) :: cg
    real(double), allocatable, dimension(:) :: x_new_pos, y_new_pos, z_new_pos
    real(double) :: energy0, energy1, max, g0, dE, gg, ggold, gamma
    integer :: i,j,k,iter,length, jj, lun, stat
    logical :: done

    allocate(cg(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating cg in control: ",ni_in_cell,stat)
    allocate(x_new_pos(ni_in_cell),y_new_pos(ni_in_cell),z_new_pos(ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating _new_pos in control: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,6*ni_in_cell,type_dbl)
    if(myid==0) write(io_lun,fmt='(/4x,"Starting CG atomic relaxation"/)')
    cg = 0.0_double
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3*ni_in_cell
    if(myid==0.AND.iprint_gen>0) write(io_lun,2) MDn_steps, MDcgtol
    energy0 = total_energy
    energy1 = 0.0_double
    dE = 0.0_double
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy0, .true., .true.)
    iter = 1
    ggold = 0.0_double
    energy1 = energy0
    do while(.NOT.done)
       ! Construct ratio for conjugacy
       gg = 0.0_double
       do j=1,ni_in_cell
          gg = gg+tot_force(1,j)*tot_force(1,j)+ &
               tot_force(2,j)*tot_force(2,j)+ &
               tot_force(3,j)*tot_force(3,j)
       end do
       if(abs(ggold).lt.1.0e-6_double) then
          gamma = 0.0_double
       else
          gamma = gg/ggold
       end if
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,*) ' CHECK :: Force Residual = ', for_conv*sqrt(gg)/ni_in_cell
       if(inode==ionode.AND.iprint_MD>2) write(io_lun,*) ' CHECK :: gamma = ', gamma
       if(CGreset) then
          if(gamma>one) then
             if(inode==ionode) write(io_lun,*) ' CG direction is reset! '
             gamma=0.0_double
          end if
       end if
       if(myid==0) write(io_lun,fmt='(/4x,"Atomic relaxation CG iteration: ",i5)') iter
       ggold = gg
       ! Build search direction
       do j=1,ni_in_cell
          jj=id_glob(j)
          cg(1,j) = gamma*cg(1,j)+tot_force(1,jj)
          cg(2,j) = gamma*cg(2,j)+tot_force(2,jj)
          cg(3,j) = gamma*cg(3,j)+tot_force(3,jj)
          x_new_pos(j) = x_atom_cell(j)
          y_new_pos(j) = y_atom_cell(j)
          z_new_pos(j) = z_atom_cell(j)
       end do
       ! Minimise in this direction
       call safemin(x_new_pos,y_new_pos,z_new_pos,cg,energy0,energy1,&
            fixed_potential, vary_mu, number_of_bands, mu, energy1)
       ! Output positions
       if(myid==0.AND.iprint_gen>1) then
          do i=1,ni_in_cell
             write(io_lun,1) i,atom_coord(1,i),atom_coord(2,i),atom_coord(3,i)
          end do
       end if
       call write_atomic_positions("UpdatedAtoms.dat",trim(pdb_template))
       ! Analyse forces
       g0 = dot(length,tot_force,1,tot_force,1)
       max = 0.0_double
       do i=1,ni_in_cell
          do k=1,3
             if(abs(tot_force(k,i))>max) max = abs(tot_force(k,i))
          end do
       end do
       ! Output and energy changes
       iter = iter + 1
       dE = energy0 - energy1
       !if(myid==0) write(io_lun,6) for_conv*max, en_units(energy_units), d_units(dist_units)
       if(myid==0) write(io_lun,4) en_conv*dE, en_units(energy_units)
       if(myid==0) write(io_lun,5) for_conv*sqrt(g0/ni_in_cell), en_units(energy_units), d_units(dist_units)
       energy0 = energy1
       !energy1 = abs(dE)
       if(iter>MDn_steps) then
          done = .true.
          if(myid==0) write(io_lun,fmt='(4x,"Exceeded number of MD steps: ",i4)') iter
       endif
       if(abs(max)<MDcgtol) then
          done = .true.
          if(myid==0) write(io_lun,fmt='(4x,"Maximum force below threshold: ",f12.5)') max
       end if
    enddo
    ! Output final positions
!    if(myid==0) call write_positions(parts)
    deallocate(z_new_pos,y_new_pos,x_new_pos,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating _new_pos in control: ",ni_in_cell,stat)
    deallocate(cg,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating cg in control: ",ni_in_cell,stat)
    call reg_dealloc_mem(area_general,6*ni_in_cell,type_dbl)
1   format(4x,'Atom ',i8,' Position ',3f15.8)
2   format(4x,'Welcome to cg_run. Doing ',i4,' steps with tolerance of ',f8.4,' ev/A')
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
!!
!!  SOURCE
!!
  subroutine md_run(fixed_potential, vary_mu, number_of_bands, mu, total_energy)

    ! Module usage
    use numbers
    use global_module, ONLY: iprint_gen, ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell, area_general
    use group_module, ONLY: parts
    use primary_module, ONLY : bundle
    use minimise, ONLY: get_E_and_F
    use move_atoms, ONLY: velocityVerlet, updateIndices
    use GenComms, ONLY: gsum, myid, my_barrier
    use GenBlas, ONLY: dot
    use force_module, ONLY: tot_force
    use io_module, ONLY: write_positions
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential

    real(double) :: number_of_bands, mu
    real(double) :: total_energy
    
    ! Local variables
    real(double), allocatable, dimension(:,:) :: velocity
    integer ::  iter, i, k, length, stat
    real(double) :: temp, KE, energy1, energy0, dE, max, g0

    allocate(velocity(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating velocity in md_run: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,3*ni_in_cell,type_dbl)
    velocity = 0.0_double
    energy0 = 0.0_double
    energy1 = 0.0_double
    dE = 0.0_double
    length = 3*ni_in_cell
    if(myid==0.AND.iprint_gen>0) write(io_lun,2) MDn_steps
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy0, .true., .false.)
    do iter = 1,MDn_steps
       if(myid==0) write(io_lun,fmt='(4x,"MD run, iteration ",i5)') iter
       call velocityVerlet(bundle,MDtimestep,temp,KE,.true.,velocity,tot_force)
       if(myid==0) write(io_lun,fmt='(4x,"Kinetic Energy          : ",f15.8)') KE
       ! Output positions
       if(myid==0.AND.iprint_gen>1) then
          do i=1,ni_in_cell
             write(io_lun,1) i,x_atom_cell(i),y_atom_cell(i),z_atom_cell(i)
          end do
       end if
       call updateIndices(.false.,fixed_potential, number_of_bands)
       call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy1, .true., .false.)
       ! Analyse forces
       g0 = dot(length,tot_force,1,tot_force,1)
       max = 0.0_double
       do i=1,ni_in_cell
          do k=1,3
             if(abs(tot_force(k,i))>max) max = tot_force(k,i)
          end do
       end do
       ! Output and energy changes
       dE = energy0 - energy1
       if(myid==0) write(io_lun,6) max
       if(myid==0) write(io_lun,4) dE
       if(myid==0) write(io_lun,5) sqrt(g0/ni_in_cell)
       energy0 = energy1
       energy1 = abs(dE)
       if(myid==0.AND.mod(iter,MDfreq)==0) call write_positions(iter,parts)
       call my_barrier
    end do
    deallocate(velocity,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating velocity in md_run: ",ni_in_cell,stat)
    call reg_dealloc_mem(area_general,3*ni_in_cell,type_dbl)
    return
1   format(4x,'Atom ',i4,' Position ',3f15.8)
2   format(4x,'Welcome to md_run. Doing ',i4,' steps')
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Energy change           : ',f15.8)
5   format(4x,'Force Residual          : ',f15.8)
6   format(4x,'Maximum force component : ',f15.8)
  end subroutine md_run
!!***

  subroutine dummy_run(fixed_potential, vary_mu, number_of_bands, mu, total_energy)

    ! Module usage
    use numbers
    use global_module, ONLY: iprint_gen, ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell
    use group_module, ONLY: parts
    use primary_module, ONLY : bundle
    use minimise, ONLY: get_E_and_F
    use move_atoms, ONLY: velocityVerlet, updateIndices2
    use GenComms, ONLY: gsum, myid, my_barrier
    use GenBlas, ONLY: dot
    use force_module, ONLY: tot_force
    use io_module, ONLY: write_positions

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential

    real(double) :: number_of_bands, mu
    real(double) :: total_energy
    
    ! Local variables
    real(double), allocatable, dimension(:,:) :: velocity
    integer :: iter, i, k, length, stat
    real(double) :: temp, KE, energy1, energy0, dE, max, g0

    MDfreq = 100
    if(myid==0.AND.iprint_gen>0) write(io_lun,2) MDn_steps
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy0, .true., .true.)
    do iter = 1,MDn_steps
       if(myid==0) write(io_lun,fmt='(4x,"Dummy run, iteration ",i5)') iter
       ! Output positions
       if(myid==0.AND.iprint_gen>1) then
          do i=1,ni_in_cell
             write(io_lun,1) i,x_atom_cell(i),y_atom_cell(i),z_atom_cell(i)
          end do
       end if
       call updateIndices2(.true.,fixed_potential, number_of_bands)
       call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy1, .true., .true.)
       ! Analyse forces
       g0 = dot(length,tot_force,1,tot_force,1)
       max = 0.0_double
       do i=1,ni_in_cell
          do k=1,3
             if(abs(tot_force(k,i))>max) max = tot_force(k,i)
          end do
       end do
       ! Output and energy changes
       dE = energy0 - energy1
       if(myid==0) write(io_lun,6) max
       if(myid==0) write(io_lun,4) dE
       if(myid==0) write(io_lun,5) sqrt(g0/ni_in_cell)
       energy0 = energy1
       energy1 = abs(dE)
       if(myid==0.AND.mod(iter,MDfreq)==0) call write_positions(iter,parts)
       call my_barrier
    end do
    if(stat/=0) call cq_abort("Error deallocating velocity in md_run: ",ni_in_cell,stat)
    return
1   format(4x,'Atom ',i4,' Position ',3f15.8)
2   format(4x,'Welcome to dummy_run. Doing ',i4,' steps')
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Energy change           : ',f15.8)
5   format(4x,'Force Residual          : ',f15.8)
6   format(4x,'Maximum force component : ',f15.8)
  end subroutine dummy_run

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
!!  SOURCE
!!
  subroutine pulay_relax(fixed_potential, vary_mu, number_of_bands, mu, total_energy)

    ! Module usage
    use numbers
    use units
    use global_module, ONLY: iprint_MD, ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell, id_glob, &
         atom_coord, rcellx, rcelly, rcellz, area_general, flag_pulay_simpleStep
    use group_module, ONLY: parts
    use minimise, ONLY: get_E_and_F
    use move_atoms, ONLY: pulayStep, velocityVerlet, updateIndices, update_atom_coord, safemin
    use GenComms, ONLY: gsum, myid
    use GenBlas, ONLY: dot
    use force_module, ONLY: tot_force
    use io_module, ONLY: write_atomic_positions, pdb_template
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use primary_module, ONLY: bundle

    implicit none

    ! Passed variables
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential

    real(double) :: number_of_bands, mu
    real(double) :: total_energy
    
    ! Local variables
    real(double), allocatable, dimension(:,:) :: velocity
    real(double), allocatable, dimension(:,:) :: cg
    real(double), allocatable, dimension(:) :: x_new_pos, y_new_pos, z_new_pos
    real(double), allocatable, dimension(:,:,:) :: posnStore
    real(double), allocatable, dimension(:,:,:) :: forceStore
    real(double) :: energy0, energy1, max, g0, dE, gg, ggold, gamma, temp, KE
    integer :: i,j,k,iter,length, jj, lun, stat, npmod, pul_mx, mx_pulay
    logical :: done!, simpleStep

    !simpleStep = .false.
    allocate(velocity(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating velocity in md_run: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,3*ni_in_cell,type_dbl)
    mx_pulay = 5
    allocate(posnStore(3,ni_in_cell,mx_pulay),forceStore(3,ni_in_cell,mx_pulay),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating cg in control: ",ni_in_cell,stat)
    allocate(cg(3,ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating cg in control: ",ni_in_cell,stat)
    allocate(x_new_pos(ni_in_cell),y_new_pos(ni_in_cell),z_new_pos(ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating _new_pos in control: ",ni_in_cell,stat)
    call reg_alloc_mem(area_general,6*ni_in_cell,type_dbl)
    if(myid==0) write(io_lun,fmt='(/4x,"Starting Pulay atomic relaxation"/)')
    posnStore = 0.0_double
    forceStore = 0.0_double
    ! Do we need to add MD.MaxCGDispl ?
    done = .false.
    length = 3*bundle%n_prim
    if(myid==0.AND.iprint_MD>0) write(io_lun,2) MDn_steps, MDcgtol
    energy0 = total_energy
    energy1 = 0.0_double
    dE = 0.0_double
    ! Find energy and forces
    call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy0, .true., .false.)
    iter = 1
    ggold = 0.0_double
    energy1 = energy0
    ! Store positions and forces on entry (primary set only)
    do i=1,ni_in_cell
       posnStore(1,i,1) = x_atom_cell(i)
       posnStore(2,i,1) = y_atom_cell(i)
       posnStore(3,i,1) = z_atom_cell(i)
       jj = id_glob(i)
       forceStore(:,i,1) = tot_force(:,jj)
    end do
    do while(.NOT.done)
       ! Book-keeping
       npmod = mod(iter, mx_pulay)+1
       pul_mx = min(iter+1, mx_pulay)
       if(iprint_MD>2) write(io_lun,*) 'npmod: ',npmod,pul_mx
       x_atom_cell = zero
       y_atom_cell = zero
       z_atom_cell = zero
       tot_force = zero
       do i=1,ni_in_cell
          if(npmod>1) then
             x_atom_cell(i) = posnStore(1,i,npmod-1)
             y_atom_cell(i) = posnStore(2,i,npmod-1)
             z_atom_cell(i) = posnStore(3,i,npmod-1)
             jj = id_glob(i)
             tot_force(:,jj) = forceStore(:,i,npmod-1)
          else
             x_atom_cell(i) = posnStore(1,i,pul_mx)
             y_atom_cell(i) = posnStore(2,i,pul_mx)
             z_atom_cell(i) = posnStore(3,i,pul_mx)
             jj = id_glob(i)
             tot_force(:,jj) = forceStore(:,i,pul_mx)
          end if
       end do
       !call gsum(x_atom_cell,ni_in_cell)
       !call gsum(y_atom_cell,ni_in_cell)
       !call gsum(z_atom_cell,ni_in_cell)
       !call gsum(tot_force,3,ni_in_cell)
       ! Take a trial step
       velocity = zero
       ! Move the atoms using velocity Verlet (why not ?)
       !call velocityVerlet(bundle,MDtimestep,temp,KE, &
       !     .false.,velocity,tot_force)
       if(flag_pulay_simpleStep) then
          write(io_lun,*) 'Step is ',MDtimestep
          do i=1,ni_in_cell
             jj = id_glob(i)
             x_atom_cell(i) = x_atom_cell(i) + MDtimestep*tot_force(1,jj)
             y_atom_cell(i) = y_atom_cell(i) + MDtimestep*tot_force(2,jj)
             z_atom_cell(i) = z_atom_cell(i) + MDtimestep*tot_force(3,jj)
          end do
          ! Do we want gsum or a scatter/gather ?
          !call gsum(x_atom_cell,ni_in_cell)
          !call gsum(y_atom_cell,ni_in_cell)
          !call gsum(z_atom_cell,ni_in_cell)
          ! Output positions
          if(myid==0.AND.iprint_MD>1) then
             do i=1,ni_in_cell
                write(io_lun,1) i,x_atom_cell(i),y_atom_cell(i),z_atom_cell(i)
             end do
          end if
          call update_atom_coord
          call updateIndices(.true.,fixed_potential, number_of_bands)
          call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy1, .true., .false.)
       else
          do j=1,ni_in_cell
             jj=id_glob(j)
             cg(1,j) = tot_force(1,jj)
             cg(2,j) = tot_force(2,jj)
             cg(3,j) = tot_force(3,jj)
             x_new_pos(j) = x_atom_cell(j)
             y_new_pos(j) = y_atom_cell(j)
             z_new_pos(j) = z_atom_cell(j)
          end do
          call safemin(x_new_pos,y_new_pos,z_new_pos,cg,energy0,energy1,&
               fixed_potential, vary_mu, number_of_bands, mu, energy1)
       end if
       ! Pass forces and positions down to pulayStep routine
       do i=1,ni_in_cell
          posnStore(1,i,npmod) = x_atom_cell(i)
          posnStore(2,i,npmod) = y_atom_cell(i)
          posnStore(3,i,npmod) = z_atom_cell(i)
          jj = id_glob(i)
         forceStore(:,i,npmod) = tot_force(:,jj)
       end do
       call pulayStep(posnStore,forceStore,x_atom_cell,y_atom_cell,z_atom_cell,mx_pulay,pul_mx)
       !call gsum(x_atom_cell,ni_in_cell)
       !call gsum(y_atom_cell,ni_in_cell)
       !call gsum(z_atom_cell,ni_in_cell)
       ! Output positions
       if(myid==0.AND.iprint_MD>1) then
          do i=1,ni_in_cell
             write(io_lun,1) i,x_atom_cell(i),y_atom_cell(i),z_atom_cell(i)
          end do
       end if
       call update_atom_coord
       call updateIndices(.true.,fixed_potential, number_of_bands)
       call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy1, .true., .false.)
       call write_atomic_positions("UpdatedAtoms.dat",trim(pdb_template))
       ! Analyse forces
       g0 = dot(length,tot_force,1,tot_force,1)
       max = 0.0_double
       do i=1,ni_in_cell
          do k=1,3
             if(abs(tot_force(k,i))>max) max = abs(tot_force(k,i))
          end do
       end do
       if(myid==0) write(io_lun,5) for_conv*sqrt(g0)/ni_in_cell, en_units(energy_units), d_units(dist_units)
       if(iter>MDn_steps) then
          done = .true.
          if(myid==0) write(io_lun,fmt='(4x,"Exceeded number of MD steps: ",i4)') iter
       endif
       if(abs(max)<MDcgtol) then
          done = .true.
          if(myid==0) write(io_lun,fmt='(4x,"Maximum force below threshold: ",f12.5)') max
       end if
       do i=1,ni_in_cell!bundle%n_prim !!! ERROR ?? !!!
          posnStore(1,i,npmod) = x_atom_cell(i)
          posnStore(2,i,npmod) = y_atom_cell(i)
          posnStore(3,i,npmod) = z_atom_cell(i)
          jj = id_glob(i)
         forceStore(:,i,npmod) = tot_force(:,jj)
       end do       
       iter = iter + 1
       energy0 = energy1
    enddo
    ! Output final positions
!    if(myid==0) call write_positions(parts)
    deallocate(z_new_pos,y_new_pos,x_new_pos,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating _new_pos in control: ",ni_in_cell,stat)
    deallocate(cg,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating cg in control: ",ni_in_cell,stat)
    call reg_dealloc_mem(area_general,6*ni_in_cell,type_dbl)
1   format(4x,'Atom ',i8,' Position ',3f15.8)
2   format(4x,'Welcome to pulay_run. Doing ',i4,' steps with tolerance of ',f8.4,' ev/A')
3   format(4x,'*** CG step ',i4,' Gamma: ',f14.8)
4   format(4x,'Energy change: ',f15.8,' ',a2)
5   format(4x,'Force Residual: ',f15.10,' ',a2,'/',a2)
6   format(4x,'Maximum force component: ',f15.8,' ',a2,'/',a2)
7   format(4x,3f15.8)
  end subroutine pulay_relax
!!***

end module control

