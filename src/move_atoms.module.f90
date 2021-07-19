! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module move_atoms
! ------------------------------------------------------------------------------
! Code area 7: Moving atoms
! ------------------------------------------------------------------------------

!!****h* Conquest/move_atoms *
!!  NAME
!!   move_atoms
!!  PURPOSE
!!   Move atoms, and update various lists
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08/01/2003
!!  MODIFICATION HISTORY
!!   07:39, 2003/01/29 dave
!!    Added constants for movement, velocity Verlet, primary and cover update
!!   08:55, 2003/02/05 dave
!!    Created update_H to reproject blips, build S, reproject pseudos,
!!    build n(r) and H
!!   14:41, 26/02/2003 drb 
!!    Added n_atoms to safemin, gsum on check to updateIndices
!!   15:57, 27/02/2003 drb & tm 
!!    Sorted out charge densities in update_H
!!   11:05, 2003/02/28 dave
!!    Added deallocation call for tm pseudopotentials
!!   10:25, 06/03/2003 drb 
!!    Corrected updating of tm pseudos
!!   15:02, 12/03/2003 drb 
!!    Tidied use statements in updateIndices
!!   13:27, 22/09/2003 drb 
!!    Added TM's changes to update positions in different arrays
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote
!!    in terms of new
!!    matrix routines
!!   2008/02/06 08:24 dave
!!    Changed for output to file not stdout
!!   2008/05/25
!!    Added timers
!!   2013/07/01 M.Arita
!!    Added sbrt: wrap_xyz_atom_cell
!!   2013/08/21 M.Arita
!!    Added sbrt: safemin2 & update_start_xyz
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2018/07/16 16:32 dave
!!    Added user-control flag_stop_on_empty_bundle
!!   2018/09/07 tsuyoshi
!!    introduced flag_debug_move_atoms for debugging
!!   2019/02/28 zamaan
!!    New subroutinee safemin_full plus dependencies for cell optimisation by 
!!    minimising a single vector
!!   2019/05/21 zamaan
!!    Removed old RNG, replaced calls with new one from rng module
!!   2019/11/18 tsuyoshi
!!    Removed the places related to flag_MDold 
!!  SOURCE
!!
module move_atoms

  use datatypes
  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_moveatoms, &
                                    tmr_std_indexing, &
                                    tmr_std_allocation

  logical :: flag_debug_move_atoms = .false.

  ! Useful physical constants
  real(double), parameter:: amu = 1.660566e-27_double
!  real(double), parameter:: ang = 1.0e-10_double
  real(double), parameter:: ang = 0.529177e-10_double
  real(double), parameter:: tscale = 1.0e-15_double
!  real(double), parameter:: ev = 1.602189e-19_double
  real(double), parameter:: ev = 2.0_double * 13.6058_double * 1.602189e-19_double
  real(double), parameter:: fac = amu*ang*ang/(tscale*tscale*ev)
  real(double), parameter:: kB = 1.3806503e-23_double
  real(double), parameter:: fac_Kelvin2Hartree = kB/ev
  !real(double), parameter:: fac_Kelvin2Hartree = 2.92126269e-6_double

  real(double) :: threshold_resetCD

  real(double) :: enthalpy_tolerance

  logical :: flag_stop_on_empty_bundle

  ! Choose line minimiser for CG
  integer :: cg_line_min
  integer, parameter :: safe = 0
  integer, parameter :: backtrack = 1
  ! Table to show the methods to update  (for update_pos_and_matrix)
   integer, parameter :: updatePos  = 0
   integer, parameter :: updateL    = 1
   integer, parameter :: updateK    = 2
   integer, parameter :: updateLorK = 3
   integer, parameter :: updateSFcoeff = 4
   integer, parameter :: extrplL    = 5
   integer, parameter :: updateX    = 6
   integer, parameter :: updateXdiss= 7
   integer, parameter :: updateS    = 8

!!***

contains

  ! --------------------------------------------------------------------
  ! Subroutine finish_blipgrid
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/finish_blipgrid *
  !!
  !!  NAME 
  !!   finish_blipgrid
  !!  USAGE
  !! 
  !!  PURPOSE
  !! 
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   08:04, 08/01/2003 dave
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine finish_blipgrid

    use set_blipgrid_module, only: free_blipgrid
    use set_bucket_module,   only: free_bucket
    use functions_on_grid,   only: dissociate_fn_on_grid

    call dissociate_fn_on_grid
    call free_bucket
    call free_blipgrid
  end subroutine finish_blipgrid
  !!***


  ! --------------------------------------------------------------------
  ! Subroutine velocityVerlet
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/velocityVerlet *
  !!
  !!  NAME 
  !!   velocityVerlet
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Moves atoms according to forces using the velocity
  !!   Verlet algorithm.  If the quenchflag is set, then quench
  !!   the motion - when v.F<0, set v=0.
  !! 
  !!   The velocity Verlet algorithm is an adaption of the Verlet
  !!   algorithm which allows calculation of the velocities in a 
  !! "better" way - it's described in "Understanding Molecular
  !!   Simulation" by Frenkel and Smit (though in a rather confusing
  !!   way - see below) or "Computer Simulation of Liquids" by
  !!   Allen and Tildesley.  The formal algorithm (as given by both
  !!   A&T and F&S) is as follows (remembering that we start with 
  !!   r(t) and v(t) and enter the routine with f(t) - a(t) = f(t)/m):
  !!
  !!   r(t+dt) = r(t) + dt.v(t) + half.a(t).dt.dt 
  !!   v(t+dt) = v(t) + half.dt.(f(t+dt)+f(t))
  !!
  !!   This is not how it's implemented - instead (as described certainly
  !!   in A&T) we do:
  !! 
  !!   v(t) = v(t-dt/2) + half.dt.a(t)
  !! [Perform analysis and output requiring v(t)]
  !!   r(t+dt) = r(t) + dt.v(t) + half.a(t).dt.dt 
  !!   v(t+dt/2) = v(t) + half.dt.a(t)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   17:09, 2003/02/04 dave (imported to Conquest from ParaDens)
  !!  MODIFICATION HISTORY
  !!   17:08, 2003/02/04 dave
  !!    Changed position to x_atom_cell
  !!   2007/08/16 15:40 dave
  !!    Bug fix for indexing of force
  !!   2008/05/25
  !!    Added timers
  !!   2011/12/09 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2013/07/01 M.Arita
  !!    The new process of wrapping atoms was introduced along with member updates
  !!   2013/08/21 M.Arita
  !!    - Added iter as a dummy variable
  !!    - Bug fix on call for update_atom_coord
  !!   2016/01/13 08:31 dave
  !!    Removed call to set_density (now included in update_H)
  !!   2020/07/28 tsuyoshi
  !!    Velocities for fixed atoms are forced to be zero.
  !!    (though this subroutine is not used now.)
  !!  TODO
  !!   Proper buffer zones for matrix mults so initialisation doesn't have
  !!   to be done at every step 03/07/2001 dave
  !!  SOURCE
  !!
  subroutine velocityVerlet(fixed_potential, prim, step, T, KE, &
                            quenchflag, velocity, force, iter)
  !ORI subroutine velocityVerlet(fixed_potential, prim, step, T, KE, &
  !ORI                           quenchflag, velocity, force)

    use datatypes
    use numbers
    use basic_types
    use global_module,  only: iprint_MD, x_atom_cell, y_atom_cell, &
                              z_atom_cell, ni_in_cell, id_glob,    &
                              flag_reset_dens_on_atom_move,        &
                              flag_move_atom
    use species_module, only: species, mass
    use GenComms,       only: myid

    implicit none

    ! Passed variables
    logical, intent(in) :: fixed_potential
    logical             :: quenchflag
    real(double)        :: step, T, KE
    type(primary_set)   :: prim
    real(double), dimension(3,ni_in_cell) :: velocity
    real(double), dimension(3,ni_in_cell) :: force
    integer             :: iter

    ! Local variables
    logical      :: flagx, flagy, flagz
    integer      :: part, memb, atom, speca, k, gatom
    real(double) :: massa, acc
    real(double) :: dx, dy, dz

    call start_timer(tmr_std_moveatoms)
    if (myid == 0 .and. iprint_MD > 0) write (io_lun,1) step, quenchflag
1   format(4x,'In velocityVerlet, timestep is ',f10.5/, &
           'Quench is ',l3)
    do atom = 1, ni_in_cell
       speca = species(atom) 
       massa = mass(speca)*fac
       gatom = id_glob(atom)
       flagx = flag_move_atom(1,gatom)
       flagy = flag_move_atom(2,gatom)
       flagz = flag_move_atom(3,gatom)
       if(quenchflag) then
          do k=1,3
             if(velocity(k,atom)*force(k,gatom)<zero) &
                  velocity(k,atom) = zero
             velocity(k,atom) = velocity(k,atom)+ &
                  step * force(k,gatom) / (two * massa)
          end do
       else
          do k=1,3
             velocity(k,atom) = velocity(k,atom)+ &
                  step*force(k,gatom) / (two * massa)
          end do
       end if
       !Now, we assume forces are forced to be zero, when
       ! flagx, y or z is false. But, I(TM) think we should
       ! have the followings, in the future. 
       !  2020/Jul/28 TM activated the following three lines, 
       !   though this subroutine is not used now.
       if(.not.flagx) velocity(1,atom) = zero
       if(.not.flagy) velocity(2,atom) = zero
       if(.not.flagz) velocity(3,atom) = zero
    end do
    ! Maybe fiddle with KE
    KE = zero
    do atom = 1, ni_in_cell
       speca = species(atom) 
       massa = mass(speca)*fac
      do k = 1, 3
       KE = KE + half * massa * velocity(k,atom) * velocity(k,atom)
      end do
    end do
    ! Update positions and velocities
    do atom = 1, ni_in_cell
       gatom = id_glob(atom)
       speca = species(atom) 
       massa = mass(speca) * fac
       flagx = flag_move_atom(1,gatom)
       flagy = flag_move_atom(2,gatom)
       flagz = flag_move_atom(3,gatom)
       ! X
       if (flagx) then
        acc = force(1,gatom) / massa
        !atom_coord_diff(1,gatom)=step*velocity(1,atom)+half*step*step*acc
        !x_atom_cell(atom) = x_atom_cell(atom) + atom_coord_diff(1,gatom)
        dx=step*velocity(1,atom)+half*step*step*acc
        x_atom_cell(atom) = x_atom_cell(atom) + dx
        velocity(1,atom)  = velocity(1,atom) + half * step * acc
       end if
       ! Y
       if (flagy) then
        acc = force(2,gatom) / massa
        !atom_coord_diff(2,gatom)=step*velocity(2,atom)+half*step*step*acc
        !y_atom_cell(atom) = y_atom_cell(atom) + atom_coord_diff(2,gatom)
        dy=step*velocity(2,atom)+half*step*step*acc
        y_atom_cell(atom) = y_atom_cell(atom) + dy
        velocity(2,atom) = velocity(2,atom) + half * step * acc
       end if
       ! Z
       if (flagz) then
        acc = force(3,gatom) / massa
        !atom_coord_diff(3,gatom)=step*velocity(3,atom)+half*step*step*acc
        !z_atom_cell(atom) = z_atom_cell(atom) + atom_coord_diff(3,gatom)
        dz=step*velocity(3,atom)+half*step*step*acc
        z_atom_cell(atom) = z_atom_cell(atom) + dz
        velocity(3,atom) = velocity(3,atom) + half * step * acc
       end if
    end do

    ! IMPORTANT: You MUST wrap atoms BEFORE updating members if they get out of the cell.
    !            Otherwise, you will get an error message at BtoG-transformation.
      call wrap_xyz_atom_cell
      call update_atom_coord
      call updateIndices3(fixed_potential,velocity)

    ! DRB 2016/01/13
    ! This line removed because this call is done in update_H
    ! NB this routine seems to be no longer called
    !%%! ! 25/Jun/2010 TM : calling set_density for SCF-MD
    !%%! if (flag_reset_dens_on_atom_move) call set_density ()
    !%%! ! 25/Jun/2010 TM : calling set_density for SCF-MD
    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine velocityVerlet
  !!***

! --------------------------------------------------------------------
! Subroutine safemin
! --------------------------------------------------------------------

!!****f* move_atoms/safemin *
!!
!!  NAME 
!!   safemin
!!  USAGE
!! 
!!  PURPOSE
!!   Finds a minimum in energy given a search direction
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   07:50, 2003/01/29 dave
!!  MODIFICATION HISTORY
!!   08:07, 2003/02/04 dave
!!    Added calls to get_E_and_F and added passed variables for
!!    get_E_and_F
!!   08:57, 2003/02/05 dave
!!    Sorted out arguments to pass to updateIndices
!!   14:41, 26/02/2003 drb 
!!    Added n_atoms from atoms
!!   08:50, 11/05/2005 dave 
!!    Added code to write out atomic positions during minimisation;
!!    also added code to subtract off atomic densities for old atomic
!!    positions and add it back on for new ones after atoms moved;
!!    this is commented out as it's not been tested or checked
!!    rigorously
!!   09:51, 25/10/2005 drb 
!!    Added correction so that the present energy is passed to
!!    get_E_and_F for blip minimisation loop
!!   15:13, 27/04/2007 drb 
!!    Reworked minimiser to be more robust; added (but not
!!    implemented) corrections to permit VERY SIMPLE charge density
!!    prediction
!!   07/11/2007 vb
!!    Added cq_abort when the trial step in safemin gets too small
!!    Changed output format for energies and brackets so that the
!!    numbers are not out of range
!!   2008/05/25
!!    Added timers
!!   2011/03/31 M.Arita
!!    Added the statements to recall set_density_pcc as atoms move.
!!   2011/12/07 L.Tong
!!    - Removed redundant dependency of density from density_module, as
!!      all usage of density are commented out. (So no spin polarisation
!!      modifications are needed.)
!!    - Changed 0.0_double to zero from numbers module
!!    - Updated calls to get_E_and_F and new_SC_potl
!!   2011/12/06 17:03 dave
!!    Bug fix for format statement
!!   2011/12/09 L.Tong
!!    Removed redundant parameter number_of_bands
!!   2012/03/27 L.Tong
!!   - Removed redundant input parameter real(double) mu
!!   2014/02/03 M.Arita
!!   - Added call for update_H because this is no longer called at updateIndices
!!   2016/01/13 08:31 dave
!!    Removed call to set_density (now included in update_H)
!!  SOURCE
!!
  subroutine safemin(start_x, start_y, start_z, direction, energy_in, &
                     energy_out, fixed_potential, vary_mu, total_energy)

    ! Module usage
    use datatypes
    use numbers
    use units
    use global_module,      only: iprint_MD, x_atom_cell, y_atom_cell,    &
                                  z_atom_cell, flag_vary_basis,           &
                                  atom_coord, ni_in_cell, rcellx, rcelly, &
                                  rcellz, flag_self_consistent,           &
                                  flag_reset_dens_on_atom_move,           &
                                  IPRINT_TIME_THRES1, flag_pcc_global
    use minimise,           only: get_E_and_F, sc_tolerance, L_tolerance, &
                                  n_L_iterations
    use GenComms,           only: my_barrier, myid, inode, ionode,        &
                                  cq_abort
    use SelfCon,            only: new_SC_potl
    use GenBlas,            only: dot
    use force_module,       only: tot_force
    use io_module,          only: write_atomic_positions, pdb_template
    use density_module,     only: density, set_density_pcc
    use maxima_module,      only: maxngrid
    use multisiteSF_module, only: flag_LFD_NonSCF
    use timer_module

    implicit none

    ! Passed variables
    real(double) :: energy_in, energy_out
    real(double), dimension(3,ni_in_cell) :: direction
    real(double), dimension(ni_in_cell)   :: start_x, start_y, start_z
    ! Shared variables needed by get_E_and_F for now (!)
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_energy
    character(len=40) :: output_file
        

    ! Local variables
    integer        :: i, j, iter, lun
    logical        :: reset_L = .false.
    logical        :: done
    type(cq_timer) :: tmr_l_iter, tmr_l_tmp1
    real(double)   :: k0, k1, k2, k3, lambda, k3old
    real(double)   :: e0, e1, e2, e3, tmp, bottom
    real(double), save :: kmin = zero, dE = zero
    real(double), dimension(:), allocatable :: store_density

    call start_timer(tmr_std_moveatoms)
    !allocate(store_density(maxngrid))
    e0 = total_energy
    if (inode == ionode .and. iprint_MD > 0) &
         write (io_lun, &
                fmt='(4x,"In safemin, initial energy is ",f20.10," ",a2)') &
               en_conv * energy_in, en_units(energy_units)
    if (inode == ionode) &
         write (io_lun, fmt='(/4x,"Seeking bracketing triplet of points"/)')
    ! Unnecessary and over cautious !
    k0 = zero
    !do i=1,ni_in_cell
    !   x_atom_cell(i) = start_x(i) + k0*direction(1,i)
    !   y_atom_cell(i) = start_y(i) + k0*direction(2,i)
    !   z_atom_cell(i) = start_z(i) + k0*direction(3,i)
    !   !write(io_lun,*) 'Position: ',i,x_atom_cell(i),y_atom_cell(i),z_atom_cell(i)
    !end do
    !!  Update atom_coord : TM 27Aug2003
    !call update_atom_coord
    !!  Update atom_coord : TM 27Aug2003
    !
    !!   Get energy and forces
    !call my_barrier
    !call updateIndices(.false.,fixed_potential, number_of_bands, &
    !     potential, density, pseudopotential, &
    !     N_GRID_MAX)
    !call get_E_and_F(output_file, n_save_freq, n_run, n_minimisation_iterations, n_support_iterations,&
    !     fixed_potential, vary_mu, n_CG_L_iterations, number_of_bands, L_tolerance, sc_tolerance, energy_tolerance, mu, &
    !     e0, potential, pseudopotential, density, expected_reduction, N_GRID_MAX)

    iter = 1
    k1 = zero
    e1 = energy_in
    k2 = k0
    e2 = e0
    e3 = e2
    !k3 = zero
    !k3old = k3
    if (kmin < 1.0e-3) then
       kmin = 0.7_double
    else
       kmin = 0.75_double * kmin
    end if
    k3 = kmin
    lambda = two
    done = .false.
    ! Loop to find a bracketing triplet
    do while (.not. done) !e3<=e2)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       !if (k2==k0) then
       !   !k3 = 0.001_double
       !   !if(abs(kmin) < RD_ERR) then
       !   if(abs(dE) < RD_ERR) then
       !      if(k3<RD_ERR) then ! First guess
       !         k3 = 0.70_double!k3old/lambda
       !      end if
       !   else
       !      tmp = dot(3*ni_in_cell,direction,1,tot_force,1)
       !      k3 = abs(0.5_double*dE/tmp)!kmin/lambda
       !      if(abs(k3)>abs(kmin)) k3 = kmin
       !   endif
       !   !   k3 = 0.1_double
       !   !else
       !   !   k3 = kmin/lambda
       !   !endif
       !elseif (k2==0.01_double) then
       !   k3 = 0.01_double
       !else
       !   k3 = lambda*k2          
       !endif
!       k3 = 0.032_double
       ! These lines calculate the difference between atomic densities and total density
       !%%!if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       !%%!   ! Subtract off atomic densities
       !%%!   store_density = density
       !%%!   call set_density()
       !%%!   density = store_density - density
       !%%!end if
       ! Move atoms
       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       do i = 1, ni_in_cell
          x_atom_cell(i) = start_x(i) + k3 * direction(1,i)
          y_atom_cell(i) = start_y(i) + k3 * direction(2,i)
          z_atom_cell(i) = start_z(i) + k3 * direction(3,i)
          if (inode == ionode .and. iprint_MD > 2) &
               write (io_lun,*) 'Position: ', i, x_atom_cell(i), &
                                y_atom_cell(i), z_atom_cell(i)
       end do
       !Update atom_coord : TM 27Aug2003
       call update_atom_coord
       !Update atom_coord : TM 27Aug2003
       ! Update indices and find energy and forces
       !call updateIndices(.false.,fixed_potential, number_of_bands)
       call updateIndices(.true., fixed_potential)
       call update_H(fixed_potential)
       ! These lines add back on the atomic densities for NEW atomic positions
       !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
          ! Add on atomic densities
          !store_density = density
          !call set_density()
          !density = store_density + density
       !end if
       ! Write out atomic positions
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
                                      trim(pdb_template))
       end if
       ! Now in update_H DRB 2016/01/13
       !if (flag_reset_dens_on_atom_move) call set_density()
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, "atom updates", IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
                           fixed_potential, vary_mu, n_L_iterations, &
                           L_tolerance, e3)
       end if
       call get_E_and_F(fixed_potential, vary_mu, e3, .false., &
                        .false.)
       if (inode == ionode .and. iprint_MD > 1) &
            write (io_lun, &
                   fmt='(4x,"In safemin, iter ",i3," step and energy &
                         &are ",2f20.10," ",a2)') &
                  iter, k3, en_conv * e3, en_units(energy_units)
       if (e3 < e2) then ! We're still going down hill
          k1 = k2
          e1 = e2
          k2 = k3
          e2 = e3
          ! New DRB 2007/04/18
          k3 = lambda * k3
          iter = iter + 1
       else if (abs(k2) < RD_ERR) then ! We've gone too far
          !k3old = k3
          !if(abs(dE)<RD_ERR) then 
          !   k3 = k3old/2.0_double
          !   dE = 1.0_double
          !else
          !   dE = 0.5_double*dE
          !end if
          !e3 = e2
          k3 = k3/lambda
       else
          done = .true.
       endif
       if (k3 <= very_small) call cq_abort("Step too small: safemin failed!")
       call stop_print_timer(tmr_l_iter, "a safemin iteration", &
                             IPRINT_TIME_THRES1)
    end do ! while (.not. done)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)  ! Final interpolation and updates
    if (inode == ionode) write(io_lun, fmt='(/4x,"Interpolating minimum"/)')
    ! Interpolate to find minimum.
    if (inode == ionode .and. iprint_MD > 1) &
            write (io_lun, fmt='(4x,"In safemin, brackets are: ",6f18.10)') &
                  k1, e1, k2, e2, k3, e3
    bottom = ((k1-k3)*(e1-e2)-(k1-k2)*(e1-e3))
    if (abs(bottom) > RD_ERR) then
       kmin = 0.5_double * (((k1*k1 - k3*k3)*(e1 - e2) -    &
                             (k1*k1 - k2*k2) * (e1 - e3)) / &
                            ((k1-k3)*(e1-e2) - (k1-k2)*(e1-e3)))
    else
       if (inode == ionode) then
          write (io_lun, fmt='(4x,"Error in safemin !")')
          write (io_lun, fmt='(4x,"Interpolation failed: ",6f15.10)') &
                k1, e1, k2, e2, k3, e3
       end if
       kmin = k2
    end if
    !%%!if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
    !%%!   ! Subtract off atomic densities
    !%%!   store_density = density
    !%%!   call set_density()
    !%%!   density = store_density - density
    !%%!end if
    do i=1,ni_in_cell
       x_atom_cell(i) = start_x(i) + kmin*direction(1,i)
       y_atom_cell(i) = start_y(i) + kmin*direction(2,i)
       z_atom_cell(i) = start_z(i) + kmin*direction(3,i)
    end do
    !Update atom_coord : TM 27Aug2003
    call update_atom_coord
    !Update atom_coord : TM 27Aug2003
    ! Check minimum: update indices and find energy and forces
    !call updateIndices(.false.,fixed_potential, number_of_bands)
    call updateIndices(.true., fixed_potential)
    call update_H(fixed_potential)
    !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       ! Add on atomic densities
       !store_density = density
       !call set_density()
       !density = store_density + density
    !end if
    if (iprint_MD > 2) then
       call write_atomic_positions("UpdatedAtoms_tmp.dat", trim(pdb_template))
    end if
    ! Now in update_H
    ! if(flag_reset_dens_on_atom_move) call set_density()
    if (flag_pcc_global) call set_density_pcc()
    call stop_print_timer(tmr_l_tmp1, &
                          "safemin - Final interpolation and updates", &
                          IPRINT_TIME_THRES1)
    ! We've just moved the atoms - we need a self-consistent ground state before we can minimise blips !
    if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
       call new_SC_potl(.false., sc_tolerance, reset_L,           &
                        fixed_potential, vary_mu, n_L_iterations, &
                        L_tolerance, e3)
    end if
    energy_out = e3
    if (iprint_MD > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .true.)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .false.)
    end if
    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, &
                fmt='(4x,"In safemin, Interpolation step and energy &
                      &are ",f15.10,f20.10," ",a2)') &
               kmin, en_conv*energy_out, en_units(energy_units)
    if (energy_out > e2 .and. abs(bottom) > RD_ERR) then
       ! The interpolation failed - go back
       call start_timer(tmr_l_tmp1,WITH_LEVEL) 
       if (inode == ionode) &
            write (io_lun,fmt='(/4x,"Interpolation failed; reverting"/)')
       kmin = k2
       !%%!if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       !%%!   ! Subtract off atomic densities
       !%%!   store_density = density
       !%%!   call set_density()
       !%%!   density = store_density - density
       !%%!end if
       do i=1,ni_in_cell
          x_atom_cell(i) = start_x(i) + kmin*direction(1,i)
          y_atom_cell(i) = start_y(i) + kmin*direction(2,i)
          z_atom_cell(i) = start_z(i) + kmin*direction(3,i)
       end do
!Update atom_coord : TM 27Aug2003
       call update_atom_coord
!Update atom_coord : TM 27Aug2003
       call updateIndices(.true., fixed_potential)
       call update_H(fixed_potential)
       !call updateIndices(.false.,fixed_potential, number_of_bands)
       !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
          ! Add on atomic densities
          !store_density = density
          !call set_density()
          !density = store_density + density
       !end if
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
                                      trim(pdb_template))
       end if
       ! Now in update_H
       !if(flag_reset_dens_on_atom_move) call set_density()
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, &
                             "safemin - Failed interpolation + Retry", &
                             IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if(flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
                           fixed_potential, vary_mu, n_L_iterations, &
                           L_tolerance, e3)
       end if
       energy_out = e3
       if (iprint_MD > 0) then
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
                           .true., .true.)
       else
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
                           .true., .false.)
       end if
    end if
    dE = e0 - energy_out
7   format(4x,3f15.8)
    if (inode == ionode .and. iprint_MD > 0) then
       write (io_lun, &
              fmt='(4x,"In safemin, exit after ",i4," &
                    &iterations with energy ",f20.10," ",a2)') &
            iter, en_conv * energy_out, en_units(energy_units)
    else if (inode == ionode) then
       write (io_lun, fmt='(/4x,"Final energy: ",f20.10," ",a2)') &
             en_conv * energy_out, en_units(energy_units)
    end if
    !deallocate(store_density)
    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine safemin
  !!***

  !!****f* move_atoms/safemin2 *
  !! PURPOSE
  !!  Carry out line minimisation in conjunction with
  !!  reusing L-matrix
  !! INPUTS
  !!
  !! AUTHOR
  !!   Michiaki Arita
  !! CREATION DATE 
  !!   2013/08/21
  !! MODIFICATION HISTORY
  !!   2013/12/02 M.Arita
  !!    - Added calls for L-matrix reconstruction & update_H
  !!   2014/02/03 M.Arita
  !!    - update_H moved outside if statement
  !!   2015/06/08 lat
  !!    - Corrected bug by adding barriers: grab_InfoGlobal
  !!   2016/01/13 08:31 dave
  !!    Removed call to set_density (now included in update_H)
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2017/05/09 dave
  !!    Adding code to load both L matrix for both spin channels
  !!   2017/11/10 dave
  !!    Removed calls to dump K matrix (now done in DMMinModule)
  !!   2018/07/11 12:07 dave
  !!    Tidying: only output on ionode, remove redundant call to wrap_xyz
  !! SOURCE
  !!
  subroutine safemin2(start_x, start_y, start_z, direction, energy_in, &
                      energy_out, fixed_potential, vary_mu, total_energy)

    ! Module usage
    use datatypes
    use numbers
    use units
    use global_module,  only: iprint_MD, x_atom_cell, y_atom_cell,    &
         z_atom_cell, flag_vary_basis,           &
         atom_coord, ni_in_cell, rcellx, rcelly, &
         rcellz, flag_self_consistent,           &
         flag_reset_dens_on_atom_move,           &
         IPRINT_TIME_THRES1, flag_pcc_global,    &
         id_glob,                                &
         flag_LmatrixReuse, flag_diagonalisation, nspin, &
         flag_SFcoeffReuse 
    use minimise,       only: get_E_and_F, sc_tolerance, L_tolerance, &
         n_L_iterations
    use GenComms,       only: my_barrier, myid, inode, ionode,        &
         cq_abort, gcopy
    use SelfCon,        only: new_SC_potl
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_atomic_positions, pdb_template
    use density_module, only: density, set_density_pcc
    use maxima_module,  only: maxngrid
    use matrix_data, ONLY: Lrange, Hrange, SFcoeff_range, SFcoeffTr_range, HTr_range
    use mult_module, ONLY: matL,L_trans, matK, matSFcoeff
    use timer_module
    use dimens, ONLY: r_super_x, r_super_y, r_super_z
    use store_matrix, ONLY: dump_pos_and_matrices
    use multisiteSF_module, only: flag_LFD_NonSCF
    !for Debugging
    use mult_module, ONLY: allocate_temp_matrix, free_temp_matrix, matrix_sum
    use global_module, ONLY: atomf, sf
    use io_module, ONLY: dump_matrix

    implicit none

    ! Passed variables
    real(double) :: energy_in, energy_out
    real(double), dimension(3,ni_in_cell) :: direction
    real(double), dimension(ni_in_cell)   :: start_x, start_y, start_z
    ! Shared variables needed by get_E_and_F for now (!)
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_energy
    character(len=40) :: output_file


    ! Local variables
    integer        :: i, j, iter, lun, gatom, stat, nfile, symm
    logical        :: reset_L = .false.
    logical        :: done
    type(cq_timer) :: tmr_l_iter, tmr_l_tmp1
    real(double)   :: k0, k1, k2, k3, lambda, k3old
    real(double)   :: e0, e1, e2, e3, tmp, bottom
    real(double), save :: kmin = zero, dE = zero
    real(double), dimension(:), allocatable :: store_density
    real(double) :: k3_old, k3_local, kmin_old

    integer :: ig, both, mat

    ! for debugging
    integer :: mat_SFcoeff_old, mat_K_old
    ! for debugging

    call start_timer(tmr_std_moveatoms)

    ! for debugging
    !    mat_SFcoeff_old = allocate_temp_matrix(SFcoeff_range,SFcoeffTr_range,atomf,sf)
    !    mat_K_old = allocate_temp_matrix(Hrange,Htr_range,sf,sf)
    ! for debugging

    !allocate(store_density(maxngrid))
    e0 = total_energy
    if (inode == ionode .and. iprint_MD > 0) &
         write (io_lun, &
         fmt='(4x,"In safemin2, initial energy is ",f20.10," ",a2)') &
         en_conv * energy_in, en_units(energy_units)
    if (inode == ionode) &
         write (io_lun, fmt='(/4x,"Seeking bracketing triplet of points"/)')
    ! Unnecessary and over cautious !

    k0 = zero

    iter = 1
    k1 = zero
    e1 = energy_in
    k2 = k0
    e2 = e0
    e3 = e2
    !k3 = zero
    !k3old = k3
    if (kmin < 1.0e-3) then
       kmin = 0.7_double
    else
       kmin = 0.75_double * kmin
    end if
    k3 = kmin
    k3_local = k3
    lambda = two
    done = .false.
    ! Loop to find a bracketing triplet
    do while (.not. done) !e3<=e2)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       call start_timer(tmr_l_tmp1, WITH_LEVEL)

       !CHECK READING K_MATRIX AND SFCOEFF   2017/12/04
       !  k3 = zero
       !     if (inode == ionode) write(io_lun,*) ' K3 = ZERO TO CHECK READING K_MATIRX AND SFCOEFF '
       ! if(flag_SFCoeffReuse) then
       !   call matrix_sum(zero, mat_SFcoeff_old, one, matSFcoeff(1))
       !   call matrix_sum(zero, mat_K_old, one, matK(1))
       !    call dump_matrix("SFcoeff_before_",  matSFcoeff(1), inode)
       ! endif
       !CHECK READING K_MATRIX AND SFCOEFF   2017/12/04

       do i = 1, ni_in_cell
          x_atom_cell(i) = start_x(i) + k3 * direction(1,i)
          y_atom_cell(i) = start_y(i) + k3 * direction(2,i)
          z_atom_cell(i) = start_z(i) + k3 * direction(3,i)
       end do

       if (ionode.EQ.inode) write (io_lun,*) "CG: 1st stage, call updateIndices3"
       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff,direction)
          !CHECK READING K_MATRIX AND SFCOEFF   2017/12/04
          !    call dump_matrix("SFcoeff_after_",  matSFcoeff(1), inode)
          !    call matrix_sum(-one, matSFcoeff(1), one, mat_SFcoeff_old)
          !    call dump_matrix("SFcoeff_diff",    matSFcoeff(1), inode)
          !    call matrix_sum(-one, mat_K_old, one, matK(1))
          !    call dump_matrix("K_diff",  mat_K_old, inode)
          !    call matrix_sum(zero, matSFcoeff(1), one, mat_SFcoeff_old)
          !CHECK READING K_MATRIX AND SFCOEFF   2017/12/04
       else
          call update_pos_and_matrices(updateLorK,direction)
       endif
       if (inode == ionode .and. iprint_MD > 2) then
          do i=1,ni_in_cell
             write (io_lun,*) 'Position: ', i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
          end do
       end if
       call update_H(fixed_potential)
       !Update start_x,start_y & start_z
       call update_start_xyz(start_x,start_y,start_z)
       ! These lines add back on the atomic densities for NEW atomic positions
       !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       ! Add on atomic densities
       !store_density = density
       !call set_density()
       !density = store_density + density
       !end if
       ! Write out atomic positions
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! Now in update_H
       !if (flag_reset_dens_on_atom_move) call set_density()
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, "atom updates", IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       call get_E_and_F(fixed_potential, vary_mu, e3, .false., &
            .false.)
       ! Now, we call dump_pos_and_matrices here. : 2018.Jan19 TM
       !  but if we want to use the information of the matrices in the beginning of this line minimisation
       !  you can comment the following line, in the future. 
       call dump_pos_and_matrices

       if (inode == ionode .and. iprint_MD > 1) &
            write (io_lun, &
            fmt='(4x,"In safemin2, iter ",i3," step and energy &
            &are ",2f20.10," ",a2)') &
            iter, k3, en_conv * e3, en_units(energy_units)
       k3_old = k3
       if (e3 < e2) then ! We're still going down hill
          k1 = k2
          e1 = e2
          k2 = k3
          e2 = e3
          ! New DRB 2007/04/18
          k3 = lambda * k3
          iter = iter + 1
          !db if (inode.EQ.ionode) write (io_lun,*) "INCREASE k3!"
       else if (abs(k2) < RD_ERR) then ! We've gone too far
          !k3old = k3
          !if(abs(dE)<very_small) then
          !   k3 = k3old/2.0_double
          !   dE = 1.0_double
          !else
          !   dE = 0.5_double*dE
          !end if
          !e3 = e2
          k3 = k3/lambda
          !db if (inode.EQ.ionode) write (io_lun,*) "DECREASE k3!"
       else
          done = .true.
          !db if (inode.EQ.ionode) write (io_lun,*) "TRUE!"
       endif
       k3_local = k3 - k3_old
       if (inode.EQ.ionode) write (io_lun,'(a,1x,3f15.10)') "k3,k3_old,k3_local:", &
            k3,k3_old,k3_local
       if (k3 <= very_small) call cq_abort("Step too small: safemin2 failed!")
       call stop_print_timer(tmr_l_iter, "a safemin2 iteration", &
            IPRINT_TIME_THRES1)
       if (inode.EQ.ionode) write (io_lun,*) "Cycle the loop! -- CG"
       if (inode.EQ.ionode) write (io_lun,*) "iter & k3:", iter, k3
    end do ! while (.not. done)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)  ! Final interpolation and updates
    if (inode == ionode) write(io_lun, fmt='(/4x,"Interpolating minimum"/)')
    ! Interpolate to find minimum.
    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, fmt='(4x,"In safemin2, brackets are: ",6f18.10)') &
         k1, e1, k2, e2, k3, e3
    bottom = ((k1-k3)*(e1-e2)-(k1-k2)*(e1-e3))
    if (abs(bottom) > very_small) then
       kmin = 0.5_double * (((k1*k1 - k3*k3)*(e1 - e2) -    &
            (k1*k1 - k2*k2) * (e1 - e3)) / &
            ((k1-k3)*(e1-e2) - (k1-k2)*(e1-e3)))
    else
       if (inode == ionode) then
          write (io_lun, fmt='(4x,"Error in safemin2 !")')
          write (io_lun, fmt='(4x,"Interpolation failed: ",6f15.10)') &
               k1, e1, k2, e2, k3, e3
       end if
       kmin = k2
    end if
    !%%!if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
    !%%!   ! Subtract off atomic densities
    !%%!   store_density = density
    !%%!   call set_density()
    !%%!   density = store_density - density
    !%%!end if
    do i=1,ni_in_cell
       x_atom_cell(i) = start_x(i) + kmin*direction(1,i)
       y_atom_cell(i) = start_y(i) + kmin*direction(2,i)
       z_atom_cell(i) = start_z(i) + kmin*direction(3,i)
    end do
    ! Get atomic displacements: atom_coord_diff(1:3, ni_in_cell)
    k3_local = kmin - k3

    if(inode==ionode.AND.iprint_MD>0) write (io_lun,*) "CG: 2nd stage"
    if(flag_SFcoeffReuse) then
       call update_pos_and_matrices(updateSFcoeff,direction)
    else
       call update_pos_and_matrices(updateLorK,direction)
    endif
    if (inode == ionode .and. iprint_MD > 2) then
       do i=1,ni_in_cell
          write (io_lun,*) 'Position: ', i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
       end do
    end if
    call update_H(fixed_potential)

    !Update start_x,start_y & start_z
    call update_start_xyz(start_x,start_y,start_z)
    !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
    ! Add on atomic densities
    !store_density = density
    !call set_density()
    !density = store_density + density
    !end if
    if (iprint_MD > 2) then
       call write_atomic_positions("UpdatedAtoms_tmp.dat", trim(pdb_template))
    end if
    ! Now in update_H
    !if(flag_reset_dens_on_atom_move) call set_density()
    if (flag_pcc_global) call set_density_pcc()
    call stop_print_timer(tmr_l_tmp1, &
         "safemin2 - Final interpolation and updates", &
         IPRINT_TIME_THRES1)
    ! We've just moved the atoms - we need a self-consistent ground state before
    ! we can minimise blips !
    if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
       call new_SC_potl(.false., sc_tolerance, reset_L,           &
            fixed_potential, vary_mu, n_L_iterations, &
            L_tolerance, e3)
    end if
    energy_out = e3
    if (iprint_MD > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .true.)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .false.)
    end if

    ! 2018.Jan19  TM
    call dump_pos_and_matrices

    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, &
         fmt='(4x,"In safemin2, Interpolation step and energy &
         &are ",f15.10,f20.10," ",a2)') &
         kmin, en_conv*energy_out, en_units(energy_units)
    ! If interpolation step failed, do interpolation AGAIN
    if (energy_out > e2 .and. abs(bottom) > RD_ERR) then
       if(e1<e3) then ! Keep k1
          if(k2<kmin) then
             k3 = kmin
             e3 = energy_out
          else
             k3 = k2
             e3 = e2
             k2 = kmin
             e2 = energy_out
          end if
       else ! Keep k3
          if(k2<kmin) then
             k1 = k2
             e1 = e2
             k2 = kmin
             e2 = energy_out
          end if
       end if
       kmin_old = kmin
       if (inode == ionode .and. iprint_MD > 1) &
            write (io_lun, fmt='(4x,"In safemin2, brackets are: ",6f18.10)') &
            k1, e1, k2, e2, k3, e3
       bottom = ((k1-k3)*(e1-e2)-(k1-k2)*(e1-e3))
       if (abs(bottom) > very_small) then
          kmin = 0.5_double * (((k1*k1 - k3*k3)*(e1 - e2) -    &
               (k1*k1 - k2*k2) * (e1 - e3)) / &
               ((k1-k3)*(e1-e2) - (k1-k2)*(e1-e3)))
          if (inode == ionode .and. iprint_MD > 1) &
               write (io_lun, &
               fmt='(4x,"In safemin2, second interpolation step is ", f15.10)') kmin
          if(kmin<k1.OR.kmin>k3) then
             if(inode == ionode .and. iprint_MD > 0) &
                  write(io_lun,*) 'Second interpolation outside limits: ',k1,k3,kmin
             dE = e0 - energy_out
             kmin = kmin_old
             if (inode == ionode .and. iprint_MD > 0) then
                write (io_lun, &
                     fmt='(4x,"In safemin2, exit after ",i4," &
                     &iterations with energy ",f20.10," ",a2)') &
                     iter, en_conv * energy_out, en_units(energy_units)
             else if (inode == ionode) then
                write (io_lun, fmt='(/4x,"Final energy: ",f20.10," ",a2)') &
                     en_conv * energy_out, en_units(energy_units)
             end if
             if (inode.EQ.ionode) write (io_lun,*) "Get out of safemin2 !" !db
             return
          end if
       else
          dE = e0 - energy_out
          if (inode == ionode .and. iprint_MD > 0) then
             write (io_lun, &
                  fmt='(4x,"In safemin2, exit after ",i4," &
                  &iterations with energy ",f20.10," ",a2)') &
                  iter, en_conv * energy_out, en_units(energy_units)
          else if (inode == ionode) then
             write (io_lun, fmt='(/4x,"Final energy: ",f20.10," ",a2)') &
                  en_conv * energy_out, en_units(energy_units)
          end if
          if (inode.EQ.ionode) write (io_lun,*) "Get out of safemin2 !" !db
          return
       end if
       do i=1,ni_in_cell
          x_atom_cell(i) = start_x(i) + kmin*direction(1,i)
          y_atom_cell(i) = start_y(i) + kmin*direction(2,i)
          z_atom_cell(i) = start_z(i) + kmin*direction(3,i)
       end do
       ! Get atomic displacements: atom_coord_diff(1:3, ni_in_cell)
       k3_local = kmin-kmin_old!03/07/2013
       !if (inode.EQ.ionode) write (io_lun,'(a,1x,3f15.10)') "k3, kmin,k3_local:", k3,kmin,k3_local

       write (io_lun,*) "CG: 3rd stage"
       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff,direction)
       else
          call update_pos_and_matrices(updateLorK,direction)
       endif
       if (inode == ionode .and. iprint_MD > 2) then
          do i=1,ni_in_cell
             write (io_lun,*) 'Position: ', i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
          end do
       end if
       call update_H(fixed_potential)

       ! Update start_x,start_y & start_z
       call update_start_xyz(start_x,start_y,start_z)!25/01/2013
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! Now in update_H
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, &
            "safemin2 - Failed interpolation + Retry", &
            IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if(flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       energy_out = e3
       if (iprint_MD > 0) then
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
               .true., .true.)
       else
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
               .true., .false.)
       end if

       ! 2018.Jan19  TM : probably we don't need to call dump_pos_and_matrices here, since
       !                  we will call it after calling safemin2
       call dump_pos_and_matrices

    end if
    dE = e0 - energy_out
7   format(4x,3f15.8)
    if (inode == ionode .and. iprint_MD > 0) then
       write (io_lun, &
            fmt='(4x,"In safemin2, exit after ",i4," &
            &iterations with energy ",f20.10," ",a2)') &
            iter, en_conv * energy_out, en_units(energy_units)
    else if (inode == ionode) then
       write (io_lun, fmt='(/4x,"Final energy: ",f20.10," ",a2)') &
            en_conv * energy_out, en_units(energy_units)
    end if
    !deallocate(store_density)

    !Debugging
    !    call free_temp_matrix(mat_K_old)
    !    call free_temp_matrix(mat_SFcoeff_old)
    !Debugging

    if (inode.EQ.ionode) write (io_lun,*) "Get out of safemin2 !" !db
    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine safemin2
 !!***

  !!****f* move_atoms/backtrack_linemin *
  !! PURPOSE
  !!  Carry out back-tracking line minimisation
  !! INPUTS
  !!
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   2019/12/09
  !! MODIFICATION HISTORY
  !!  2020/01/08 12:52 dave
  !!   Bug fix: reset alpha to one on entry
  !! SOURCE
  !!
  subroutine backtrack_linemin(direction, energy_in, &
                      energy_out, fixed_potential, vary_mu, total_energy)

    ! Module usage
    use datatypes
    use numbers
    use units
    use global_module,  only: iprint_MD, x_atom_cell, y_atom_cell,    &
         z_atom_cell, flag_vary_basis,           &
         atom_coord, ni_in_cell, rcellx, rcelly, &
         rcellz, flag_self_consistent,           &
         flag_reset_dens_on_atom_move,           &
         IPRINT_TIME_THRES1, flag_pcc_global,    &
         id_glob,                                &
         flag_LmatrixReuse, flag_diagonalisation, nspin, &
         flag_SFcoeffReuse 
    use minimise,       only: get_E_and_F, sc_tolerance, L_tolerance, &
         n_L_iterations
    use GenComms,       only: my_barrier, myid, inode, ionode,        &
         cq_abort, gcopy
    use SelfCon,        only: new_SC_potl
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_atomic_positions, pdb_template
    use density_module, only: density, set_density_pcc
    use maxima_module,  only: maxngrid
    use matrix_data, ONLY: Lrange, Hrange, SFcoeff_range, SFcoeffTr_range, HTr_range
    use mult_module, ONLY: matL,L_trans, matK, matSFcoeff
    use timer_module
    use dimens, ONLY: r_super_x, r_super_y, r_super_z
    use store_matrix, ONLY: dump_pos_and_matrices
    use multisiteSF_module, only: flag_LFD_NonSCF
    use mult_module, ONLY: allocate_temp_matrix, free_temp_matrix, matrix_sum
    use global_module, ONLY: atomf, sf
    use io_module, ONLY: dump_matrix
    use force_module,      only: force

    implicit none

    ! Passed variables
    real(double) :: energy_in, energy_out
    real(double), dimension(3,ni_in_cell) :: direction
    ! Shared variables needed by get_E_and_F for now (!)
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_energy

    ! Local variables
    integer        :: i, j, iter, lun, gatom, stat, nfile, symm
    logical        :: reset_L = .false.
    logical        :: done
    type(cq_timer) :: tmr_l_iter, tmr_l_tmp1
    real(double)   :: alpha_new, armijo, grad_f_dot_p, grad_fp_dot_p, old_alpha
    real(double)   :: e0, e1, e2, e3, tmp, bottom
    real(double), save :: kmin = zero, dE = zero
    real(double), dimension(:), allocatable :: store_density
    real(double) :: k3_old, k3_local, kmin_old
    real(double) :: alpha = one
    real(double) :: c1, c2

    integer :: ig, both, mat

    call start_timer(tmr_std_moveatoms)

    iter = 0
    old_alpha = zero
    alpha = one
    e0 = total_energy
    if (inode == ionode .and. iprint_MD > 0) &
         write (io_lun, &
         fmt='(4x,"In backtrack_linemin, initial energy is ",f16.6," ",a2)') &
         en_conv * energy_in, en_units(energy_units)

    c1 = 0.1_double
    c2 = 0.9_double
    ! grad f dot p  Note that the ordering of direction and tot_force is different
    grad_f_dot_p = zero
    do i=1, ni_in_cell
       j = id_glob(i)
       grad_f_dot_p = grad_f_dot_p - direction(1,i)*tot_force(1,j)
       grad_f_dot_p = grad_f_dot_p - direction(2,i)*tot_force(2,j)
       grad_f_dot_p = grad_f_dot_p - direction(3,i)*tot_force(3,j)
    end do
    if(inode==ionode.AND.iprint_MD>1) &
         write(io_lun, fmt='(2x,"Starting backtrack_linemin, magnitude of grad_f.p is ",e16.6)') &
         sqrt(-grad_f_dot_p/ni_in_cell)
    done = .false.
    do while (.not. done)
       iter = iter+1
       ! Take a step along search direction
       do i = 1, ni_in_cell
          x_atom_cell(i) = x_atom_cell(i) + (alpha - old_alpha) * direction(1,i)
          y_atom_cell(i) = y_atom_cell(i) + (alpha - old_alpha) * direction(2,i)
          z_atom_cell(i) = z_atom_cell(i) + (alpha - old_alpha) * direction(3,i)
       end do

       ! Update and find new energy
       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff,direction)
       else
          call update_pos_and_matrices(updateLorK,direction)
       endif
       if (inode == ionode .and. iprint_MD > 2) then
          do i=1,ni_in_cell
             write (io_lun,fmt='(2x,"Position: ",i3,3f13.8)') i, &
                  x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
          end do
       end if
       call update_H(fixed_potential)
       ! Write out atomic positions
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       call get_E_and_F(fixed_potential, vary_mu, e3, .false., &
            .false.)
       !call dump_pos_and_matrices
       ! e3 is f(x + alpha p)
       armijo = e0 + c1 * alpha * grad_f_dot_p

       if (inode == ionode .and. iprint_MD > 1) then
          write (io_lun, &
               fmt='(4x,"In backtrack_linemin, iter ",i3," step and energy &
               &are ",2f16.6," ",a2)') &
               iter, alpha, en_conv * e3, en_units(energy_units)
          write(io_lun, fmt='(6x,"Armijo threshold is ",f16.6," ",a2)') armijo, en_units(energy_units)
       end if
       if(e3<armijo) then ! success
          done = .true.
       else
          old_alpha = alpha
          alpha_new = (-half * alpha * grad_f_dot_p) / ((e3 - e0)/alpha - grad_f_dot_p)
          alpha = max(alpha_new, 0.1_double*alpha)
       end if
    end do ! while (.not. done)
    energy_out = e3
    call dump_pos_and_matrices
    ! Now find forces
    call force(fixed_potential, vary_mu, n_L_iterations, &
         L_tolerance, sc_tolerance, energy_out, .true.)
    ! Evaluate new grad f dot p
    grad_fp_dot_p = zero
    do i=1, ni_in_cell
       j = id_glob(i)
       grad_fp_dot_p = grad_f_dot_p - direction(1,i)*tot_force(1,j)
       grad_fp_dot_p = grad_f_dot_p - direction(2,i)*tot_force(2,j)
       grad_fp_dot_p = grad_f_dot_p - direction(3,i)*tot_force(3,j)
    end do
    if(inode==ionode.AND.iprint_MD>3) &
         write(io_lun,fmt='(6x,"In backtrack_linemin, second Wolfe: ",e11.4," < ",e11.4)') &
         abs(grad_fp_dot_p), c2*abs(grad_f_dot_p)

    dE = e0 - energy_out
    if (inode == ionode .and. iprint_MD > 2) then
       write (io_lun, &
            fmt='(4x,"In backtrack_linemin, exit after ",i4," &
            &iterations with energy ",f16.6," ",a2)') &
            iter, en_conv * energy_out, en_units(energy_units)
    else if (inode == ionode .and. iprint_MD > 0) then
       write (io_lun, fmt='(/4x,"In backtrack_linemin, final energy is   ",f16.6," ",a2)') &
            en_conv * energy_out, en_units(energy_units)
    end if

    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine backtrack_linemin
!!***
  
  !!****f* move_atoms/adapt_backtrack_linemin *
  !! PURPOSE
  !!  Carry out back-tracking line minimisation
  !!  Tweak to search for large step size initially
  !! INPUTS
  !!
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !!   2019/12/09
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine adapt_backtrack_linemin(direction, energy_in, &
                      energy_out, fixed_potential, vary_mu, total_energy)

    ! Module usage
    use datatypes
    use numbers
    use units
    use global_module,  only: iprint_MD, x_atom_cell, y_atom_cell,    &
         z_atom_cell, flag_vary_basis,           &
         atom_coord, ni_in_cell, rcellx, rcelly, &
         rcellz, flag_self_consistent,           &
         flag_reset_dens_on_atom_move,           &
         IPRINT_TIME_THRES1, flag_pcc_global,    &
         id_glob,                                &
         flag_LmatrixReuse, flag_diagonalisation, nspin, &
         flag_SFcoeffReuse 
    use minimise,       only: get_E_and_F, sc_tolerance, L_tolerance, &
         n_L_iterations
    use GenComms,       only: my_barrier, myid, inode, ionode,        &
         cq_abort, gcopy
    use SelfCon,        only: new_SC_potl
    use GenBlas,        only: dot
    use force_module,   only: tot_force
    use io_module,      only: write_atomic_positions, pdb_template
    use density_module, only: density, set_density_pcc
    use maxima_module,  only: maxngrid
    use matrix_data, ONLY: Lrange, Hrange, SFcoeff_range, SFcoeffTr_range, HTr_range
    use mult_module, ONLY: matL,L_trans, matK, matSFcoeff
    use timer_module
    use dimens, ONLY: r_super_x, r_super_y, r_super_z
    use store_matrix, ONLY: dump_pos_and_matrices
    use multisiteSF_module, only: flag_LFD_NonSCF
    use mult_module, ONLY: allocate_temp_matrix, free_temp_matrix, matrix_sum
    use global_module, ONLY: atomf, sf
    use io_module, ONLY: dump_matrix
    use force_module,      only: force

    implicit none

    ! Passed variables
    real(double) :: energy_in, energy_out
    real(double), dimension(3,ni_in_cell) :: direction
    ! Shared variables needed by get_E_and_F for now (!)
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_energy

    ! Local variables
    integer        :: i, j, iter, lun, gatom, stat, nfile, symm
    logical        :: reset_L = .false.
    logical        :: done
    type(cq_timer) :: tmr_l_iter, tmr_l_tmp1
    real(double)   :: alpha_new, armijo, grad_f_dot_p, grad_fp_dot_p, old_alpha
    real(double)   :: e0, e1, e2, e3, tmp, bottom
    real(double), save :: kmin = zero, dE = zero
    real(double), dimension(:), allocatable :: store_density
    real(double) :: k3_old, k3_local, kmin_old
    real(double), save :: alpha = one
    real(double), save :: scale = 0.9_double
    real(double) :: c1, c2

    integer :: ig, both, mat

    call start_timer(tmr_std_moveatoms)

    iter = 0
    old_alpha = zero
    e0 = total_energy
    if (inode == ionode .and. iprint_MD > 0) &
         write (io_lun, &
         fmt='(4x,"In backtrack_linemin, initial energy is ",f16.6," ",a2)') &
         en_conv * energy_in, en_units(energy_units)

    c1 = 0.1_double
    c2 = 0.9_double
    ! grad f dot p
    grad_f_dot_p = zero
    do i=1, ni_in_cell
       j = id_glob(i)
       grad_f_dot_p = grad_f_dot_p - direction(1,i)*tot_force(1,j)
       grad_f_dot_p = grad_f_dot_p - direction(2,i)*tot_force(2,j)
       grad_f_dot_p = grad_f_dot_p - direction(3,i)*tot_force(3,j)
    end do
    if(inode==ionode.AND.iprint_MD>1) &
         write(io_lun, fmt='(2x,"Starting backtrack_linemin, grad_f.p is ",f16.6)') grad_f_dot_p
    done = .false.
    iter = 0
    do while (.not. done)
       iter = iter+1
       ! Take a step along sesarch direction
       do i = 1, ni_in_cell
          x_atom_cell(i) = x_atom_cell(i) + (alpha - old_alpha) * direction(1,i)
          y_atom_cell(i) = y_atom_cell(i) + (alpha - old_alpha) * direction(2,i)
          z_atom_cell(i) = z_atom_cell(i) + (alpha - old_alpha) * direction(3,i)
       end do

       ! Update and find new energy
       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff,direction)
       else
          call update_pos_and_matrices(updateLorK,direction)
       endif
       if (inode == ionode .and. iprint_MD > 2) then
          do i=1,ni_in_cell
             write (io_lun,fmt='(2x,"Pos: ",i3,3f13.8)') i, &
                  x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
          end do
       end if
       call update_H(fixed_potential)
       ! Write out atomic positions
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       call get_E_and_F(fixed_potential, vary_mu, e3, .false., &
            .false.)
       !call dump_pos_and_matrices
       ! e3 is f(x + alpha p)
       armijo = e0 + c1 * alpha * grad_f_dot_p

       if (inode == ionode .and. iprint_MD > 1) then
          write (io_lun, &
               fmt='(4x,"In backtrack_linemin, iter ",i3," step and energy &
               &are ",2f20.10," ",a2)') &
               iter, alpha, en_conv * e3, en_units(energy_units)
          write(io_lun, fmt='(6x,"Armijo threshold is ",f16.6," ",a2)') armijo, en_units(energy_units)
       end if
       if(e3<armijo) then ! success
          done = .true.
       else
          old_alpha = alpha
          alpha_new = (-half * alpha * grad_f_dot_p) / ((e3 - e0)/alpha - grad_f_dot_p)
          alpha = max(alpha_new, 0.1_double*alpha)
       end if
    end do ! while (.not. done)
    energy_out = e3
    call dump_pos_and_matrices
    ! Test increase of alpha
    if(iter==1) then
       scale = scale*0.9_double
       alpha = (one+scale)*alpha
    end if
    ! Now find forces
    call force(fixed_potential, vary_mu, n_L_iterations, &
         L_tolerance, sc_tolerance, energy_out, .true.)
    grad_fp_dot_p = zero
    do i=1, ni_in_cell
       j = id_glob(i)
       grad_fp_dot_p = grad_f_dot_p - direction(1,i)*tot_force(1,j)
       grad_fp_dot_p = grad_f_dot_p - direction(2,i)*tot_force(2,j)
       grad_fp_dot_p = grad_f_dot_p - direction(3,i)*tot_force(3,j)
    end do
    if(inode==ionode.AND.iprint_MD>3) &
         write(io_lun,fmt='(6x,"In backtrack_linemin, second Wolfe: ",e11.4," < ",e11.4)') &
         abs(grad_fp_dot_p), c2*abs(grad_f_dot_p)
    dE = e0 - energy_out
    if (inode == ionode .and. iprint_MD > 0) then
       write (io_lun, &
            fmt='(4x,"In backtrack_linemin, exit after ",i4," &
            &iterations with energy ",f20.10," ",a2)') &
            iter, en_conv * energy_out, en_units(energy_units)
    else if (inode == ionode .and. iprint_MD > 0) then
       write (io_lun, fmt='(/4x,"Final energy: ",f20.10," ",a2)') &
            en_conv * energy_out, en_units(energy_units)
    end if

    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine adapt_backtrack_linemin
!!***
  
  !!****f* move_atoms/safemin_cell *
  !! PURPOSE
  !! Optimize the simulation cell dimensions a b and c
  !! Heavily borrowed from previous safemin subroutines
  !! INPUTS
  !!
  !! AUTHOR
  !!   Jack Baker
  !!   Shereif Mujahed
  !! CREATION DATE
  !!   2017/05/12
  !! MODIFICATION HISTORY
  !!   2017/05/25 dave
  !!    Added more variables that need updating following cell vector changes,
  !!    notably k point locations and reciprocal lattice vectors for FFTs
  !!   2017/06/20 J.S.B
  !!    Moved Dave's changes to a separate subroutine "update_cell_dims" for
  !!    clarity.
  !!   2017/11/10 dave
  !!    Removed calls to dump K matrix (now done in DMMinModule)
  !!    Sometime since September implemented fractional coordinates for coord diffs
  !!    so that uniform volume scaling of cell doesn't generate atom position changes
  !!   2018/01/22 tsuyoshi (with dave)
  !!    Updated to use new matrix rebuilding following atom movement
  !!   2019/02/28 zamaan
  !!    Modified to minimise enthalpy instead of energy, relax to target
  !!    pressure
  !! SOURCE

  subroutine safemin_cell(start_rcellx, start_rcelly, start_rcellz, &
                          search_dir_x, search_dir_y, search_dir_z, &
                          search_dir_mean, target_press, enthalpy_in, enthalpy_out, &
                          fixed_potential, vary_mu)

    ! Module usage


    use datatypes
    use numbers
    use units
    use global_module,      only: iprint_MD, x_atom_cell, y_atom_cell,    &
                                  z_atom_cell, flag_vary_basis,           &
                                  atom_coord, ni_in_cell, rcellx, rcelly, &
                                  rcellz, flag_self_consistent,           &
                                  flag_reset_dens_on_atom_move,           &
                                  IPRINT_TIME_THRES1, flag_pcc_global, &
                                  flag_diagonalisation, cell_constraint_flag, &
                                  flag_SFcoeffReuse
    use minimise,           only: get_E_and_F, sc_tolerance, L_tolerance, &
                                  n_L_iterations
    use GenComms,           only: my_barrier, myid, inode, ionode, cq_abort, &
                                  gcopy
    use SelfCon,            only: new_SC_potl
    use GenBlas,            only: dot
    use force_module,       only: tot_force
    use io_module,          only: write_atomic_positions, pdb_template
    use density_module,     only: density, set_density_pcc
    use maxima_module,      only: maxngrid
    use multisiteSF_module, only: flag_LFD_NonSCF
    use timer_module
    use dimens,             only: r_super_x, r_super_y, r_super_z, &
                                  r_super_x_squared, r_super_y_squared, &
                                  r_super_z_squared, volume, &
                                  grid_point_volume, &
                                  one_over_grid_point_volume, n_grid_x, &
                                  n_grid_y, n_grid_z
    use fft_module,         only: recip_vector, hartree_factor, i0
    use DiagModule,         only: kk, nkp
    use store_matrix,       only: dump_pos_and_matrices

    implicit none

    ! Passed variables
    real(double), intent(in)  :: enthalpy_in, target_press, start_rcellx, &
                                 start_rcelly, start_rcellz, search_dir_x, &
                                 search_dir_y, search_dir_z, search_dir_mean
    real(double), intent(out) :: enthalpy_out
    ! Shared variables needed by get_E_and_F for now (!)
    logical, intent(in)       :: vary_mu, fixed_potential

    ! Local variables
    integer        :: i, j, iter, lun, stat, nfile, symm
    logical        :: reset_L = .false.
    logical        :: done
    type(cq_timer) :: tmr_l_iter, tmr_l_tmp1
    real(double)   :: k0, k1, k2, k3, lambda, k3old, orcellx, orcelly, orcellz, scale, ratio
    real(double)   :: e0, e1, e2, e3, tmp, bottom, xvec, yvec, zvec, r2, &
                      h0, h1, h2, h3, dH, energy_out
    real(double), save :: kmin = zero, dE = zero
    real(double), dimension(:), allocatable :: store_density
    real(double), dimension(3,ni_in_cell) :: direction

    call start_timer(tmr_std_moveatoms)
    direction = zero
    !allocate(store_density(maxngrid))
    h0 = enthalpy_in
    if (inode == ionode .and. iprint_MD > 0) &
         write (io_lun, &
         fmt='(4x,"In safemin_cell, initial enthalpy is ",f20.10," ",a2)') &
         en_conv * enthalpy_in, en_units(energy_units)
    if (inode == ionode .and. iprint_MD > 0) &
         write (io_lun, fmt='(/4x,"Seeking bracketing triplet of points"/)')
    ! Unnecessary and over cautious !
    k0 = zero
    iter = 1
    k1 = zero
    h1 = enthalpy_in
    k2 = k0
    h2 = h0
    h3 = h2
    !k3 = zero
    !k3old = k3
    if (kmin < 1.0e-3) then
       kmin = 0.7_double
    else
       kmin = 0.75_double * kmin
    end if
    k3 = kmin
    lambda = 1.5_double ! Reduce severity
    done = .false.
    ! Loop to find a bracketing triplet
    do while (.not. done) !e3<=e2)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       ! get new lattice vectors
       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       ! DRB added 2017/05/24 17:13
       ! Keep previous cell to allow scaling
       ! update_cell_dims updates the cell according to the user set
       ! constraints.
       call update_cell_dims(start_rcellx, start_rcelly, &
            start_rcellz, search_dir_x, search_dir_y, search_dir_z,&
            k3, iter, search_dir_mean)

       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff,direction)
       else
          call update_pos_and_matrices(updateLorK,direction)
       endif
       call update_H(fixed_potential)

       ! These lines add back on the atomic densities for NEW atomic positions
       ! Write out atomic positions
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! Now in update_H DRB 2016/01/13
       !if (flag_reset_dens_on_atom_move) call set_density()
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, "atom updates", IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       call get_E_and_F(fixed_potential, vary_mu, e3, .false., &
            .false.)
       call dump_pos_and_matrices

       h3 = enthalpy(e3, target_press)
       if (inode == ionode .and. iprint_MD > 1) &
            write (io_lun, &
            fmt='(4x,"In safemin_cell, iter ",i3," step and enthalpy &
            &are ",2f20.10," ",a2)') &
            iter, k3, en_conv * h3, en_units(energy_units)
       write(io_lun,*) "h3 is", h3, "h2 is", h2
       write(io_lun,*) "k1 is", k1, "k2 is", k2, "k3 is", k3
       if (h3 < h2) then ! We're still going down hill
          if (inode == ionode .and. iprint_MD > 0) write(io_lun,*) "h3 larger than e2. Going downhill"
          k1 = k2
          h1 = h2
          k2 = k3
          h2 = h3
          ! New DRB 2007/04/18
          k3 = lambda * k3
          iter = iter + 1
       else if (abs(k2) < RD_ERR) then ! We've gone too far
          k3 = k3/lambda
          if (inode == ionode .and. iprint_MD > 0) write(io_lun,*) "Gone too far"
       else
          done = .true.
       endif
       if (k3 <= very_small) call cq_abort("Step too small: safemin_cell failed!")
       call stop_print_timer(tmr_l_iter, "a safemin_cell iteration", &
            IPRINT_TIME_THRES1)
    end do !while (.not. done)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)  ! Final interpolation and updates
    if (inode == ionode .and. iprint_MD > 0) write(io_lun, fmt='(/4x,"Interpolating minimum"/)')
    ! Interpolate to find minimum.
    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, fmt='(4x,"In safemin_cell, brackets are: ",6f18.10)') &
         k1, h1, k2, h2, k3, h3
    bottom = ((k1-k3)*(h1-h2)-(k1-k2)*(h1-h3))
    if (abs(bottom) > RD_ERR) then
       kmin = 0.5_double * (((k1*k1 - k3*k3)*(h1 - h2) -    &
            (k1*k1 - k2*k2) * (h1 - h3)) / &
            ((k1-k3)*(h1-h2) - (k1-k2)*(h1-h3)))
    else
       if (inode == ionode .and. iprint_MD > 0) then
          write (io_lun, fmt='(4x,"Error in safemin_cell !")')
          write (io_lun, fmt='(4x,"Interpolation failed: ",6f15.10)') &
               k1, h1, k2, h2, k3, h3
       end if
       kmin = k2
    end if

    if (inode == ionode .and. iprint_MD > 0) write(io_lun,*) 'kmin is ',kmin
    call update_cell_dims(start_rcellx, start_rcelly, &
         start_rcellz, search_dir_x, search_dir_y, search_dir_z,&
         kmin, iter, search_dir_mean)
    if(flag_SFcoeffReuse) then
       call update_pos_and_matrices(updateSFcoeff,direction)
    else
       call update_pos_and_matrices(updateLorK,direction)
    endif
    call update_H(fixed_potential)

    !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
    ! Add on atomic densities
    !store_density = density
    !call set_density()
    !density = store_density + density
    !end if
    if (iprint_MD > 2) then
       call write_atomic_positions("UpdatedAtoms_tmp.dat", trim(pdb_template))
    end if
    ! Now in update_H
    ! if(flag_reset_dens_on_atom_move) call set_density()
    if (flag_pcc_global) call set_density_pcc()
    call stop_print_timer(tmr_l_tmp1, &
         "safemin_cell - Final interpolation and updates", &
         IPRINT_TIME_THRES1)
    ! We've just moved the atoms - we need a self-consistent ground state before we can minimise blips !
    if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
       call new_SC_potl(.false., sc_tolerance, reset_L,           &
            fixed_potential, vary_mu, n_L_iterations, &
            L_tolerance, e3)
    end if
    h3 = enthalpy(e3, target_press)
    enthalpy_out = h3
    if (iprint_MD > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .true.)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .false.)
    end if
    call dump_pos_and_matrices
    enthalpy_out = enthalpy(energy_out, target_press)

    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, &
         fmt='(4x,"In safemin_cell, Interpolation step and enthalpy &
         &are ",f15.10,f20.10," ",a2)') &
         kmin, en_conv*enthalpy_out, en_units(energy_units)
    if (enthalpy_out > h2 .and. abs(bottom) > RD_ERR) then
       ! The interpolation failed - go back
       call start_timer(tmr_l_tmp1,WITH_LEVEL)
       if (inode == ionode .and. iprint_MD > 0) &
            write (io_lun,fmt='(/4x,"Interpolation failed; reverting"/)')
       kmin = k2
       ! DRB added 2017/05/24 17:13
       ! Keep previous cell to allow scaling
       call update_cell_dims(start_rcellx, start_rcelly, &
            start_rcellz, search_dir_x, search_dir_y, search_dir_z,&
            kmin, iter, search_dir_mean)
       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff,direction)
       else
          call update_pos_and_matrices(updateLorK,direction)
       endif

       call update_H(fixed_potential)
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! Now in update_H
       !if(flag_reset_dens_on_atom_move) call set_density()
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, &
            "safemin_cell - Failed interpolation + Retry", &
            IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if(flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       h3 = enthalpy(e3, target_press)
       enthalpy_out = e3
       if (iprint_MD > 0) then
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
               .true., .true.)
       else
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
               .true., .false.)
       end if
       ! we may not need to call dump_pos_and_matrices here. (if it would be called in the part after calling safemin_cell)
       call dump_pos_and_matrices  
    end if
    dH = h0 - enthalpy_out
7   format(4x,3f15.8)
    if (inode == ionode .and. iprint_MD > 0) then
       write (io_lun, &
            fmt='(4x,"In safemin_cell, exit after ",i4," &
            &iterations with enthalpy ",f20.10," ",a2)') &
            iter, en_conv * enthalpy_out, en_units(energy_units)
    else if (inode == ionode) then
       write (io_lun, fmt='(/4x,"Final enthalpy: ",f20.10," ",a2)') &
            en_conv * enthalpy_out, en_units(energy_units)
    end if
    !deallocate(store_density)
    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine safemin_cell
  !!***

  !!****f* move_atoms/safemin_full *
  !! PURPOSE
  !!  Carry out line minimisation of cell + ionic degrees of freedom 
  !! in conjunction with reusing L-matrix (adapted from safemin2)
  !! INPUTS
  !!
  !! AUTHOR
  !!   Zamaan Raza
  !! CREATION DATE 
  !!   2019/02/06
  !! MODIFICATION HISTORY
  !!
  !! SOURCE
  !!
  subroutine safemin_full(config, force, cell_ref, enthalpy_in, enthalpy_out, &
                          target_press, fixed_potential, vary_mu, &
                          total_enthalpy)

    ! Module usage
    use datatypes
    use numbers
    use units
    use global_module,  only: iprint_MD, x_atom_cell, y_atom_cell,    &
                              z_atom_cell, flag_vary_basis,           &
                              ni_in_cell, rcellx, rcelly, rcellz,     &
                              flag_self_consistent,                   &
                              IPRINT_TIME_THRES1, flag_pcc_global,    &
                              flag_LmatrixReuse, flag_SFcoeffReuse,   &
                              atom_coord
    use minimise,       only: get_E_and_F, sc_tolerance, L_tolerance, &
                              n_L_iterations
    use GenComms,       only: inode, ionode, cq_abort
    use SelfCon,        only: new_SC_potl
    use io_module,      only: write_atomic_positions, pdb_template
    use density_module, only: set_density_pcc
    use timer_module
    use store_matrix, ONLY: dump_pos_and_matrices
    use multisiteSF_module, only: flag_LFD_NonSCF

    implicit none

    ! Passed variables
    real(double) :: enthalpy_in, enthalpy_out, target_press
    real(double), dimension(:,:) :: config, force
    real(double), dimension(:)   :: cell_ref
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_enthalpy


    ! Local variables
    integer        :: i, j, iter, lun, gatom, stat, nfile, symm
    logical        :: reset_L = .false.
    logical        :: done
    type(cq_timer) :: tmr_l_iter, tmr_l_tmp1
    real(double)   :: k0, k1, k2, k3, lambda, k3old, energy_out
    real(double)   :: e0, e1, e2, e3, h0, h1, h2, h3, tmp, bottom
    real(double), save :: kmin = zero, dH = zero
    real(double), dimension(:), allocatable :: store_density
    real(double) :: k3_old, k3_local, kmin_old
    real(double), dimension(:,:), allocatable :: config_start
    real(double), dimension(3,ni_in_cell) :: dummy
    real(double)  :: orcellx, orcelly, orcellz

    real(double) :: dx, dy, dz, d

    integer :: ig, both, mat

    ! for debugging
    integer :: mat_SFcoeff_old, mat_K_old
    ! for debugging

    allocate(config_start(3,ni_in_cell+1))

    if (inode==ionode .and. iprint_MD > 2) &
      write(io_lun,'(2x,a)') "move_atoms/safemin_full"

    call start_timer(tmr_std_moveatoms)

    h0 = total_enthalpy
    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, &
         fmt='(4x,"In safemin_full, initial enthalpy is ",f20.10," ",a2)') &
         en_conv * enthalpy_in, en_units(energy_units)
    if (inode == ionode) &
         write (io_lun, fmt='(/4x,"Seeking bracketing triplet of points"/)')
    ! Unnecessary and over cautious !

    config_start = config

    k0 = zero

    iter = 1
    k1 = zero
    h1 = enthalpy_in
    k2 = k0
    h2 = h0
    h3 = h2
    if (kmin < 1.0e-3) then
       kmin = 0.3_double
    else
       kmin = 0.75_double * kmin
    end if
    k3 = kmin
    k3_local = k3
    lambda = two
    done = .false.

    ! Loop to find a bracketing triplet
    do while (.not. done) !e3<=e2)
       call start_timer(tmr_l_iter, WITH_LEVEL)
       call start_timer(tmr_l_tmp1, WITH_LEVEL)

       ! Change the configuration
       call propagate_vector(force, config_start, config, cell_ref, k3)
       call vector_to_cq(config, cell_ref, orcellx, orcelly, orcellz)

       if (ionode==inode) write(io_lun,*) "CG: 1st stage"

       if (flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff, force(:,1:ni_in_cell))
       else
          call update_pos_and_matrices(updateLorK, force(:,1:ni_in_cell))
       end if

       if (inode==ionode .and. iprint_MD>2) then
          write(io_lun,'(4x,a)') "atom_coord:"
          do i=1,ni_in_cell
             write (io_lun,'(6x,i8,3f16.8)') i, atom_coord(:,i)
             ! write (io_lun,*) 'Position: ', i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
          end do
       end if
       call update_H(fixed_potential)

       ! Write out atomic positions
       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! Now in update_H
       !if (flag_reset_dens_on_atom_move) call set_density()
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, "atom updates", IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       call get_E_and_F(fixed_potential, vary_mu, e3, .false., &
            .false.)
       ! Now, we call dump_pos_and_matrices here. : 2018.Jan19 TM
       !  but if we want to use the information of the matrices in the beginning of this line minimisation
       !  you can comment the following line, in the future. 
       call dump_pos_and_matrices

       h3 = enthalpy(e3, target_press)
       if (inode == ionode .and. iprint_MD > 1) &
            write (io_lun, &
            fmt='(4x,"In safemin_full, iter ",i3," step and enthalpy &
            &are ",2f20.10," ",a2)') &
            iter, k3, en_conv * h3, en_units(energy_units)
       k3_old = k3
       if (h3 < h2) then ! We're still going down hill
          k1 = k2
          h1 = h2
          k2 = k3
          h2 = h3
          k3 = lambda * k3
          iter = iter + 1
          ! if (inode.EQ.ionode) write (io_lun,*) "INCREASE k3!"
       ! zamaan - should be k2 < small surely?
       else if (abs(k2) < RD_ERR) then ! We've gone too far
          k3 = k3/lambda
          ! if (inode.EQ.ionode) write (io_lun,*) "DECREASE k3!"
       else
          done = .true.
          ! if (inode.EQ.ionode) write (io_lun,*) "TRUE!"
       end if
       k3_local = k3 - k3_old
       if (inode.EQ.ionode) write (io_lun,'(a,1x,3f15.10)') "k3,k3_old,k3_local:", &
            k3,k3_old,k3_local
       if (k3 <= very_small) call cq_abort("Step too small: safemin_full failed!")
       call stop_print_timer(tmr_l_iter, "a safemin_full iteration", &
            IPRINT_TIME_THRES1)
       if (inode.EQ.ionode) write (io_lun,*) "Cycle the loop! -- CG"
       if (inode.EQ.ionode) write (io_lun,*) "iter & k3:", iter, k3
    end do ! while (.not. done)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)  ! Final interpolation and updates
    if (inode == ionode) write(io_lun, fmt='(/4x,"Interpolating minimum"/)')
    ! Interpolate to find minimum.
    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, fmt='(4x,"In safemin_full, brackets are: ",6f18.10)') &
         k1, h1, k2, h2, k3, h3
    bottom = ((k1-k3)*(h1-h2)-(k1-k2)*(h1-h3))
    if (abs(bottom) > very_small) then
       kmin = 0.5_double * (((k1*k1 - k3*k3)*(h1 - h2) -    &
            (k1*k1 - k2*k2) * (h1 - h3)) / &
            ((k1-k3)*(h1-h2) - (k1-k2)*(h1-h3)))
    else
       if (inode == ionode) then
          write (io_lun, fmt='(4x,"Error in safemin_full !")')
          write (io_lun, fmt='(4x,"Interpolation failed: ",6f15.10)') &
               k1, h1, k2, h2, k3, h3
       end if
       kmin = k2
    end if

    call propagate_vector(force, config_start, config, cell_ref, kmin)
    call vector_to_cq(config, cell_ref, orcellx, orcelly, orcellz)
    k3_local = kmin - k3

    if(inode==ionode.AND.iprint_MD>0) write (io_lun,*) "CG: 2nd stage"
    if(flag_SFcoeffReuse) then
       call update_pos_and_matrices(updateSFcoeff,force(:,1:ni_in_cell))
    else
       call update_pos_and_matrices(updateLorK,force(:,1:ni_in_cell))
    endif

    if (inode == ionode .and. iprint_MD > 1) then
       write(io_lun,'(4x,a)') "atom_coord:"
       do i=1,ni_in_cell
          write (io_lun,'(6x,i8,3f16.8)') i, atom_coord(:,i)
          ! write (io_lun,*) 'Position: ', i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
       end do
    end if
    call update_H(fixed_potential)

    if (iprint_MD > 2) then
       call write_atomic_positions("UpdatedAtoms_tmp.dat", trim(pdb_template))
    end if
    ! Now in update_H
    !if(flag_reset_dens_on_atom_move) call set_density()
    if (flag_pcc_global) call set_density_pcc()
    call stop_print_timer(tmr_l_tmp1, &
         "safemin_full - Final interpolation and updates", &
         IPRINT_TIME_THRES1)
    ! We've just moved the atoms - we need a self-consistent ground state before
    ! we can minimise blips !
    if (flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
       call new_SC_potl(.false., sc_tolerance, reset_L,           &
            fixed_potential, vary_mu, n_L_iterations, &
            L_tolerance, e3)
    end if
    energy_out = e3
    if (iprint_MD > 0) then
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .true.)
    else
       call get_E_and_F(fixed_potential, vary_mu, energy_out, .true., .false.)
    end if

    ! 2018.Jan19  TM
    call dump_pos_and_matrices

    enthalpy_out = enthalpy(energy_out, target_press)
    if (inode == ionode .and. iprint_MD > 1) &
         write (io_lun, &
         fmt='(4x,"In safemin_full, Interpolation step and enthalpy &
         &are ",f15.10,f20.10," ",a2)') &
         kmin, en_conv*enthalpy_out, en_units(energy_units)
    ! If interpolation step failed, do interpolation AGAIN
    if (enthalpy_out > h2 .and. abs(bottom) > RD_ERR) then
       if(h1<h3) then ! Keep k1
          if(k2<kmin) then
             k3 = kmin
             h3 = enthalpy_out
          else
             k3 = k2
             h3 = h2
             k2 = kmin
             h2 = enthalpy_out
          end if
       else ! Keep k3
          if(k2<kmin) then
             k1 = k2
             h1 = h2
             k2 = kmin
             h2 = enthalpy_out
          end if
       end if
       kmin_old = kmin
       if (inode == ionode .and. iprint_MD > 1) &
            write (io_lun, fmt='(4x,"In safemin_full, brackets are: ",6f18.10)') &
            k1, h1, k2, h2, k3, h3
       bottom = ((k1-k3)*(h1-h2)-(k1-k2)*(h1-h3))
       if (abs(bottom) > very_small) then
          kmin = 0.5_double * (((k1*k1 - k3*k3)*(h1 - h2) -    &
               (k1*k1 - k2*k2) * (h1 - h3)) / &
               ((k1-k3)*(h1-h2) - (k1-k2)*(h1-h3)))
          if (inode == ionode .and. iprint_MD > 1) &
               write (io_lun, &
               fmt='(4x,"In safemin_full, second interpolation step is ", f15.10)') kmin
          if(kmin<k1.OR.kmin>k3) then
             if(inode == ionode .and. iprint_MD > 1) &
                  write(io_lun,*) 'Second interpolation outside limits: ',k1,k3,kmin
             dH = h0 - enthalpy_out
             kmin = kmin_old
             if (inode == ionode .and. iprint_MD > 1) then
                write (io_lun, &
                     fmt='(4x,"In safemin_full, exit after ",i4," &
                     &iterations with enthlapy ",f20.10," ",a2)') &
                     iter, en_conv * enthalpy_out, en_units(energy_units)
             else if (inode == ionode) then
                write (io_lun, fmt='(/4x,"Final enthalpy: ",f20.10," ",a2)') &
                     en_conv * enthalpy_out, en_units(energy_units)
             end if
             if (inode.EQ.ionode) write (io_lun,*) "Get out of safemin_full !" !db
             return
          end if
       else
          dH = h0 - enthalpy_out
          if (inode == ionode .and. iprint_MD > 1) then
             write (io_lun, &
                  fmt='(4x,"In safemin_full, exit after ",i4," &
                  &iterations with enthalpy ",f20.10," ",a2)') &
                  iter, en_conv * enthalpy_out, en_units(energy_units)
          else if (inode == ionode) then
             write (io_lun, fmt='(/4x,"Final enthalpy: ",f20.10," ",a2)') &
                  en_conv * enthalpy_out, en_units(energy_units)
          end if
          if (inode==ionode .and. iprint_MD > 1) &
            write (io_lun,'(2x,a)') "Get out of safemin_full!" 
          return
       end if

       call propagate_vector(force, config_start, config, cell_ref, kmin)
       call vector_to_cq(config, cell_ref, orcellx, orcelly, orcellz)
       k3_local = kmin-kmin_old!03/07/2013

       write (io_lun,*) "CG: 3rd stage"
       if(flag_SFcoeffReuse) then
          call update_pos_and_matrices(updateSFcoeff, force(:,1:ni_in_cell))
       else
          call update_pos_and_matrices(updateLorK, force(:,1:ni_in_cell))
       endif

       if (inode == ionode .and. iprint_MD > 2) then
          write(io_lun,'(4x,a)') "atom_coord:"
          do i=1,ni_in_cell
             write (io_lun,'(6x,i8,3f16.8)') i, atom_coord(:,i)
             ! write (io_lun,*) 'Position: ', i, x_atom_cell(i), y_atom_cell(i), z_atom_cell(i)
          end do
       end if
       call update_H(fixed_potential)

       if (iprint_MD > 2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat", &
               trim(pdb_template))
       end if
       ! Now in update_H
       if (flag_pcc_global) call set_density_pcc()
       call stop_print_timer(tmr_l_tmp1, &
            "safemin_full - Failed interpolation + Retry", &
            IPRINT_TIME_THRES1)
       ! We've just moved the atoms - we need a self-consistent ground
       ! state before we can minimise blips !
       if(flag_vary_basis .or. .NOT.flag_LFD_NonSCF) then
          call new_SC_potl(.false., sc_tolerance, reset_L,           &
               fixed_potential, vary_mu, n_L_iterations, &
               L_tolerance, e3)
       end if
       energy_out = e3
       if (iprint_MD > 0) then
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
               .true., .true.)
       else
          call get_E_and_F(fixed_potential, vary_mu, energy_out, &
               .true., .false.)
       end if

       ! 2018.Jan19  TM : probably we don't need to call dump_pos_and_matrices here, since
       !                  we will call it after calling safemin2
       call dump_pos_and_matrices

    end if
    dH = h0 - enthalpy_out
7   format(4x,3f15.8)
    if (inode==ionode .and. iprint_MD > 1) then
       write (io_lun, &
            fmt='(4x,"In safemin_full, exit after ",i4," &
            &iterations with enthalpy ",f20.10," ",a2)') &
            iter, en_conv * enthalpy_out, en_units(energy_units)
    else if (inode == ionode) then
       write (io_lun, fmt='(/4x,"Final enthalpy: ",f20.10," ",a2)') &
            en_conv * enthalpy_out, en_units(energy_units)
    end if

    deallocate(config_start)
    if (inode==ionode) write (io_lun,*) "Get out of safemin_full !" !db
    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine safemin_full
  !!***

  !!****f* move_atoms/update_start_xyz *
  !!
  !!NAME 
  !! update_start_xyz
  !!USAGE
  !! 
  !!PURPOSE
  !! Updates start_x, start_y and start_z after updating member info.
  !!INPUTS
  !! 
  !!USES
  !! 
  !!AUTHOR
  !! Michiaki Arita
  !!CREATION DATE
  !! 2013/08/21
  !!MODIFICATION HISTORY
  !!
  !!SOURCE
  !!
  subroutine update_start_xyz(x,y,z)

    ! Module usage
    use datatypes
    use global_module, ONLY: ni_in_cell,id_glob,id_glob_inv_old, &
                             flag_MDdebug,Iprint_MDdebug, iprint_MD
    use GenComms, ONLY: cq_abort
    ! DB
    use input_module, ONLY: io_assign, io_close
    use GenComms, ONLY: inode,ionode

    implicit none

    ! passed variable
    real(double) :: x(ni_in_cell),y(ni_in_cell),z(ni_in_cell)

    ! local variables
    integer :: ni,stat_alloc
    integer :: id_global,ni_old
    real(double), allocatable :: x_tmp(:),y_tmp(:),z_tmp(:)
    ! DB
    integer :: lun_db
    character(7) :: file_name

    if (inode==ionode .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "move_atoms/update_start_xyz"

    allocate (x_tmp(ni_in_cell),y_tmp(ni_in_cell),z_tmp(ni_in_cell), &
              STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error allocating x_tmp:', ni_in_cell)

    x_tmp=x ; y_tmp=y ; z_tmp=z
    ! Update x,y & z
    do ni = 1, ni_in_cell
      id_global = id_glob(ni)
      ni_old = id_glob_inv_old(id_global)
      x(ni) = x_tmp(ni_old)
      y(ni) = y_tmp(ni_old)
      z(ni) = z_tmp(ni_old)
    enddo

    deallocate (x_tmp,y_tmp,z_tmp, STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error deallocating x_tmp,y_tmp and z_tmp:', ni_in_cell)

    !! ---- DEBUG: 25/01/2013 ---- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.2) then
      if (inode.EQ.ionode) then
        call io_assign(lun_db)
        open (lun_db,file='xyz.dat',position='append')
        !write (lun_db,*) "safemin2 iter:", iter
        write (lun_db,*) "safemin2:"
        do ni = 1, ni_in_cell
          write (lun_db,'(a,1x,3f15.10)') "ni, start_x,y,z:", x(ni), y(ni), z(ni)
        enddo
        write (lun_db,*) ""
        call io_close(lun_db)
      endif
    endif
    !! ---- DEBUG: 25/01/2013 ---- !!

    return
  end subroutine update_start_xyz
  !!***

  !!****f* move_atoms/pulayStep *
  !!
  !!NAME 
  !! pulayStep
  !!USAGE
  !! 
  !!PURPOSE
  !! Relaxes the atoms to their minimum energy positions using the
  !! guaranteed reduction Pulay algorithm (see Chem. Phys. Lett. 325, 796
  !! (2000) for more details - also minimise).  Take a step with the atoms
  !! based on the timestep and then minimise the norm of the force vector.
  !!INPUTS
  !! 
  !! 
  !!USES
  !! 
  !!AUTHOR
  !! D.R.Bowler
  !!CREATION DATE
  !! 17/07/2001
  !!MODIFICATION HISTORY
  !! 20/07/2001 dave
  !!  Changed so that loops only go over primary set atoms
  !! 2008/05/25
  !!  Added timers
  !! 2012/05/26 L.Tong
  !! - Added input npmod, this is used by the new version of DoPulay
  !!   2019/10/24 11:52 dave
  !!    Changed function calls to DoPulay
  !!SOURCE
  !!
  subroutine pulayStep(npmod, posnStore, forceStore, x_atom_cell, &
                       y_atom_cell, z_atom_cell, mx_pulay, pul_mx)

    use datatypes
    use global_module,  only: iprint_MD, ni_in_cell, id_glob_inv
    use numbers
    use GenBlas,        only: dot, axpy
    use GenComms,       only: gsum, myid, inode, ionode
    use Pulay,          only: DoPulay
    use primary_module, only: bundle

    implicit none

    ! Passed variables
    integer :: npmod, mx_pulay, pul_mx
    real(double), dimension(3,ni_in_cell,mx_pulay) :: forceStore
    real(double), dimension(3,ni_in_cell,mx_pulay) :: posnStore
    real(double), dimension(ni_in_cell)            :: x_atom_cell
    real(double), dimension(ni_in_cell)            :: y_atom_cell
    real(double), dimension(ni_in_cell)            :: z_atom_cell

    ! Local variables
    integer      :: i,j, length, jj
    real(double) :: gg
    real(double), dimension(mx_pulay,mx_pulay) :: Aij
    real(double), dimension(mx_pulay)          :: alph

    call start_timer(tmr_std_moveatoms)
    length = 3*ni_in_cell
    Aij = zero
    do i=1,pul_mx
       do j=1,pul_mx
          gg = dot(length, forceStore(1:,1:,j),1, &
               forceStore(1:,1:,i),1)
          Aij(j,i) = gg
          !write(io_lun,fmt='(4x,"A is : ",2i3,f22.17)') i,j,Aij(j,i)
       enddo
    enddo
    !call gsum(Aij,mx_pulay,mx_pulay)
    call DoPulay(npmod,Aij,alph,pul_mx,mx_pulay)
    if(myid==0.AND.iprint_MD>2) write(io_lun,*) 'Alpha: ', alph
    x_atom_cell(:) = 0.0_double
    y_atom_cell(:) = 0.0_double
    z_atom_cell(:) = 0.0_double
    do i=1,pul_mx
       do j=1,ni_in_cell
          jj = id_glob_inv(j)
          x_atom_cell(jj) = x_atom_cell(jj) + alph(i)*posnStore(1,j,i)
          y_atom_cell(jj) = y_atom_cell(jj) + alph(i)*posnStore(2,j,i)
          z_atom_cell(jj) = z_atom_cell(jj) + alph(i)*posnStore(3,j,i)
       enddo
    enddo
    call stop_timer(tmr_std_moveatoms)
    return
  end subroutine pulayStep
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine updateIndices
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/updateIndices *
  !!
  !!  NAME 
  !!   updateIndices
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Updates the indices for matrices, saves relevant information
  !!   and stores (if necessary) old L matrix
  !!
  !!   At the simplest, all this does is update the positions of the
  !!   atoms in the primary and covering sets, and rebuild the
  !!   Hamiltonian
  !!  INPUTS
  !!   logical :: matrix_update Flags whether the user wants ALL
  !!   matrix information updated
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   08:55, 2003/01/29 dave
  !!  MODIFICATION HISTORY
  !!   08:13, 2003/02/04 dave
  !!    Added blipgrid initialisation and reinitialisation calls
  !!   14:42, 26/02/2003 drb 
  !!    Added gsum on check
  !!   08:36, 2003/03/12 dave
  !!    Removed unnecessary use of H_matrix_module
  !!   09:01, 2003/11/10 dave
  !!    D'oh ! Put in a call to cover_update for ewald_CS so that the
  !!    new ewald routines work
  !!   08:49, 11/05/2005 dave 
  !!    Added lines which check for change in band energy, and reset
  !!    DM if too large; these are commented out as these ideas are
  !!    not rigorously tested
  !!   2006/09/08 07:59 dave
  !!    Various changes for dynamic allocation
  !!   2008/05/25
  !!    Added timers
  !!   2011/09/29 16:48 M. Arita
  !!    CS is updated for DFT-D2
  !!   2011/11/17 10:18 dave
  !!    Updated call to set_blipgrid
  !!   2011/12/09 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2014/02/03 M.Arita
  !!    Removed call for update_H
  !!   2015/11/24 08:32 dave
  !!    Removed old ewald calls
  !!   2016/09/16 17:00 nakata
  !!    Used RadiusAtomf instead of RadiusSupport
  !!  TODO
  !!   Think about updating radius component of matrix derived type,
  !!   or eliminating it !
  !!  SOURCE
  !!
  subroutine updateIndices(matrix_update, fixed_potential)

    ! Module usage
    use datatypes
    use mult_module,            only: fmmi, immi
    use matrix_module,          only: allocate_matrix,                &
                                      deallocate_matrix,              &
                                      set_matrix_pointers2, matrix
    use group_module,           only: parts
    use cover_module,           only: BCS_parts, DCS_parts, ion_ion_CS, &
                                      D2_CS
    use primary_module,         only: bundle
    use global_module,          only: iprint_MD, x_atom_cell,         &
                                      y_atom_cell, z_atom_cell,       &
                                      IPRINT_TIME_THRES2,             &
                                      flag_Becke_weights, flag_dft_d2, flag_diagonalisation
    use matrix_data,            only: Hrange, mat, rcut
    use maxima_module,          only: maxpartsproc
    use set_blipgrid_module,    only: set_blipgrid
    use set_bucket_module,      only: set_bucket
    use dimens,                 only: r_core_squared,r_h,             &
                                      RadiusAtomf
    use pseudopotential_common, only: core_radius
    use GenComms,               only: myid, cq_abort, gsum
    use functions_on_grid,      only: associate_fn_on_grid
    use numbers
    use timer_module
    use density_module,         only: build_Becke_weights
    use DiagModule, only: end_scalapack_format, init_scalapack_format
    
    implicit none

    ! Passed variables
    logical, intent(in) :: matrix_update

    ! Shared variables needed by get_H_matrix for now (!)
    logical :: fixed_potential

    ! Local variables
    logical        :: check
    integer        :: i,k,stat
    type(cq_timer) :: tmr_l_tmp1,tmr_l_tmp2

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    ! Update positions in primary and covering sets
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle, parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, DCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, ion_ion_CS, parts)
    if (flag_dft_d2) &
         call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, D2_CS, parts)
    ! If there's a new interaction of Hamiltonian range, then we REALLY need to rebuild the matrices etc
    call checkBonds(check,bundle,BCS_parts,mat(1,Hrange),maxpartsproc,rcut(Hrange))
    ! If one processor gets a new bond, they ALL need to redo the indices
    call gsum(check)
    ! There's also an option for the user to force it via matrix_update (which could be set to every n iterations ?)
    if(check.OR.matrix_update) then
       call start_timer(tmr_l_tmp2,WITH_LEVEL)
       if(flag_diagonalisation) call end_scalapack_format
       ! Deallocate all matrix storage
       ! finish blip-grid indexing
       call finish_blipgrid
       ! finish matrix multiplication indexing
       call fmmi(bundle)
       ! Reallocate and find new indices
       call immi(parts,bundle,BCS_parts,myid+1)
       ! Reallocate for blip grid
       call set_blipgrid(myid, RadiusAtomf, core_radius)
       !call set_blipgrid(myid,r_h,sqrt(r_core_squared))
       call set_bucket(myid)
       call associate_fn_on_grid
       if(flag_diagonalisation) call init_scalapack_format
       call stop_print_timer(tmr_l_tmp2,"matrix reindexing",IPRINT_TIME_THRES2)
    end if
    if (flag_Becke_weights) call build_Becke_weights
    ! Rebuild S, n(r) and hamiltonian based on new positions
    !call update_H(fixed_potential)
    call stop_print_timer(tmr_l_tmp1,"indices update",IPRINT_TIME_THRES2)
    return
  end subroutine updateIndices
  !!***


  !!****f* move_atoms/updateIndices2 *
  !! PURPOSE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   David Bowler
  !! CREATION DATE 
  !! MODIFICATION HISTORY
  !!   2011/12/09 L.Tong
  !!     Removed redundant parameter number_of_bands
  !!   2014/02/03 M.Arita
  !!     Removed call for update_H
  !!   2015/11/24 08:32 dave
  !!    Removed old ewald calls
  !!   2016/09/16 17:00 nakata
  !!    Used RadiusAtomf instead of RadiusSupport
  !! SOURCE
  !!
  subroutine updateIndices2(matrix_update, fixed_potential)

    ! Module usage
    use datatypes
    use mult_module,            only: fmmi, immi
    use matrix_module,          only: allocate_matrix,                &
                                      deallocate_matrix,              &
                                      set_matrix_pointers2, matrix
    use group_module,           only: parts
    use cover_module,           only: BCS_parts, DCS_parts, ion_ion_CS, &
                                      D2_CS
    use primary_module,         only: bundle
    use global_module,          only: iprint_MD, x_atom_cell,         &
                                      y_atom_cell, z_atom_cell,       &
                                      flag_Becke_weights, flag_dft_d2, flag_diagonalisation
    use matrix_data,            only: Hrange, mat, rcut
    use maxima_module,          only: maxpartsproc
    use set_blipgrid_module,    only: set_blipgrid
    use set_bucket_module,      only: set_bucket
    use dimens,                 only: r_core_squared,r_h,             &
                                      RadiusAtomf
    use pseudopotential_common, only: core_radius
    use GenComms,               only: myid, cq_abort, gsum
    use functions_on_grid,      only: associate_fn_on_grid
    use density_module,         only: build_Becke_weights
    use numbers
    use DiagModule, only: end_scalapack_format, init_scalapack_format

    implicit none

    ! Passed variables
    logical, intent(in) :: matrix_update

    ! Shared variables needed by get_H_matrix for now (!)
    logical :: fixed_potential

    ! Local variables
    logical :: check
    integer :: i, k, stat

    ! Update positions in primary and covering sets
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle, parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, DCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, &
         z_atom_cell, ion_ion_CS, parts)
    if (flag_dft_d2) call cover_update(x_atom_cell, y_atom_cell, &
         z_atom_cell, D2_CS, parts)
    check = .false.
    ! If there's a new interaction of Hamiltonian range, then we
    ! REALLY need to rebuild the matrices etc
    call checkBonds(check,bundle,BCS_parts,mat(1,Hrange),maxpartsproc,rcut(Hrange))
    ! If one processor gets a new bond, they ALL need to redo the indices
    call gsum(check)
    ! There's also an option for the user to force it via
    ! matrix_update (which could be set to every n iterations ?)
    if(check.OR.matrix_update) then
       if(flag_diagonalisation) call end_scalapack_format
       ! Deallocate all matrix storage
       ! finish blip-grid indexing
       call finish_blipgrid
       ! finish matrix multiplication indexing
       call fmmi(bundle)
       ! Reallocate and find new indices
       call immi(parts,bundle,BCS_parts,myid+1,1)
       ! Reallocate for blip grid
       call set_blipgrid(myid, RadiusAtomf, core_radius)
       !call set_blipgrid(myid,r_h,sqrt(r_core_squared))
       call set_bucket(myid)
       call associate_fn_on_grid
       if(flag_diagonalisation) call init_scalapack_format
    end if
    if(flag_Becke_weights) call build_Becke_weights
    ! Rebuild S, n(r) and hamiltonian based on new positions
    !call update_H (fixed_potential)
    return
  end subroutine updateIndices2
  !!*****

  !!****f* move_atoms/updateIndices3 *
  !! PURPOSE
  !!  Updates the member information in each partition
  !! INPUTS
  !!  fixed_potential,velocity,step,iteration
  !!   - iteration is optional
  !!   - iteration will be deleted in the next update
  !! AUTHOR
  !!   Michiaki Arita
  !! CREATION DATE 
  !!   2013/07/02
  !! MODIFICATION HISTORY
  !!   2013/08/20 M.Arita
  !!    -  Implemented L-matrix reconstruction
  !!    -  Deleted step and iteration
  !!   2013/12/02 M.Arita
  !!    -  Deleted calls for L-matrix reconstruction and update_H. They are
  !!       called at md_run and safemin2
  !!    -  Added calls for initialising and finalising matrix indexing for XL-BOMD
  !!   2015/11/24 08:32 dave
  !!    Removed old ewald calls
  !!   2016/09/16 17:00 nakata
  !!    Used RadiusAtomf instead of RadiusSupport
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2018/07/11 12:08 dave
  !!    Added routines to redistribute atoms to partitions and partitions to
  !!    processes if an empty bundle is found
  !!   2019/11/18 14:37 dave
  !!    Updates to rebuild covering sets if cell varies during run
  !! SOURCE
  !!
  subroutine updateIndices3(fixed_potential,velocity)

    ! Module usage
    use datatypes
    use global_module, ONLY: flag_Becke_weights,flag_dft_d2, flag_variable_cell, id_glob, &
                             ni_in_cell,x_atom_cell,y_atom_cell,z_atom_cell,      &
                             IPRINT_TIME_THRES2,glob2node, io_lun,     &
                             flag_XLBOMD, flag_diagonalisation, flag_neutral_atom, &
                             numprocs, atom_coord, species_glob, iprint_MD
    use GenComms, ONLY: inode,ionode,my_barrier,myid,gcopy, cq_abort
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts, DCS_parts, ion_ion_CS, D2_CS, BCS_blocks, &
         make_cs,make_iprim,send_ncover, deallocate_cs
    use mult_module,            ONLY: fmmi,immi
    use set_blipgrid_module, ONLY: set_blipgrid
    use set_bucket_module, ONLY: set_bucket
    use dimens, ONLY: RadiusAtomf
    use pseudopotential_common, ONLY: core_radius
    use functions_on_grid, ONLy: associate_fn_on_grid
    use density_module, ONLY: build_Becke_weights
    use UpdateMember_module, ONLY: updateMembers_group, updateMembers_cs
    use atoms, ONLY: distribute_atoms,deallocate_distribute_atom
    use timer_module
    use numbers
    use io_module, ONLY: append_coords,write_atomic_positions,pdb_template
    use UpdateInfo, ONLY: make_glob2node
    use XLBOMD_module, ONLY: immi_XL,fmmi_XL
    use group_module,   ONLY: blocks, deallocate_group_set, make_cc2
    use primary_module, ONLY: deallocate_primary_set, bundle, make_prim, domain
    use construct_module, ONLY: init_primary
    use sfc_partitions_module, ONLY: sfc_partitions_to_processors
    use ion_electrostatic, ONLY: ewald_real_cutoff, ion_ion_cutoff
    use species_module, ONLY: species
    use matrix_data,    ONLY: rcut,max_range
    use dimens,         ONLY: r_core_squared,r_h, r_dft_d2
    use DiagModule, only: end_scalapack_format, init_scalapack_format
    use maxima_module, ONLY: maxpartsproc, maxatomsproc, maxatomspart

    implicit none

     ! Passed variables
    real(double) :: velocity(3,ni_in_cell)
    logical :: fixed_potential

    ! Local variables
    integer :: nfile,symm,np, ni, id_global, ni_old
    logical :: append_coords_bkup, flag_empty_bundle
    type(cq_timer) :: tmr_l_tmp1,tmr_l_tmp2
    real(double) :: rcut_max
    real(double) :: velocity_tmp(3,ni_in_cell)


    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    ! Update members in bundle and check for empty bundle
    call updateMembers_group(velocity, flag_empty_bundle)
    if(flag_empty_bundle.and.flag_stop_on_empty_bundle) &
       call cq_abort("Empty bundle detected: user set stop_on_empty_bundle, so stopping...")
    ! Update CS member locations
    !if( (.NOT.(flag_variable_cell)) .AND. (.NOT.flag_empty_bundle)) &
    !     call updateMembers_cs(velocity)
    ! Start updates
    call start_timer(tmr_l_tmp2,WITH_LEVEL)
    if(flag_diagonalisation) call end_scalapack_format
    ! finish blip-grid indexing
    call finish_blipgrid
    ! finish matrix multiplication indexing
    if (flag_XLBOMD) call fmmi_XL()
    call fmmi(bundle)
    ! Now we need to redistribute; if the cell is changing or one process
    ! has no atoms then we must rebuild the covering sets
    if(flag_empty_bundle.OR.flag_variable_cell) then
       ! Deallocate parts and covering sets
       call deallocate_cs(BCS_parts,.true.)
       call deallocate_cs(DCS_parts,.true.)
       call deallocate_cs(BCS_blocks,.false.)
       call deallocate_cs(ion_ion_CS,.true.)
       if(flag_dft_d2) call deallocate_cs(D2_CS,.true.)
       call deallocate_distribute_atom
       ! If one process has no atoms then we have to redistribute the
       ! overall workload; in the longer term, we could trigger this
       ! if the load balancing becomes poor
       if(flag_empty_bundle) then
          if(inode==ionode) &
               write(io_lun,fmt='(2x,"Empty bundle detected: redistributing atoms between processes")')
          call deallocate_primary_set(bundle)
          call deallocate_group_set(parts)
          ! Call Hilbert curve
          call sfc_partitions_to_processors(parts)
          ! inverse table to npnode
          do np=1,parts%ngcellx*parts%ngcelly*parts%ngcellz
             parts%inv_ngnode(parts%ngnode(np))=np
          end do
          call make_cc2(parts,numprocs)
          ! NB  velocity update is done in update_pos_and_matrices
          do ni = 1, ni_in_cell
             id_global= id_glob(ni)
             x_atom_cell(ni) = atom_coord(1,id_global)
             y_atom_cell(ni) = atom_coord(2,id_global)
             z_atom_cell(ni) = atom_coord(3,id_global)
             species(ni)     = species_glob(id_global)
          end do
          ! Covering sets are made in setgrid
          ! Create primary set for atoms: bundle of partitions
          call init_primary(bundle, maxatomsproc, maxpartsproc, .true.)
          call make_prim(parts, bundle, inode-1, id_glob, x_atom_cell, &
               y_atom_cell, z_atom_cell, species)
       end if
       ! Sorts out which processor owns which atoms
       call distribute_atoms(inode, ionode)
       call make_cs(inode-1, rcut(max_range), BCS_parts, parts, bundle, &
            ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell)
       call make_iprim(BCS_parts, bundle)
       call send_ncover(BCS_parts, inode)
       call my_barrier
       ! Reallocate and find new indices
       call immi(parts,bundle,BCS_parts,myid+1)
       if (flag_XLBOMD) call immi_XL(parts,bundle,BCS_parts,myid+1)
       rcut_max = max(sqrt(r_core_squared),r_h) + very_small
       call make_cs(myid,rcut_max, DCS_parts , parts , domain, &
            ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell)
       call make_cs(myid,rcut_max, BCS_blocks, blocks, bundle)
       call send_ncover(DCS_parts, myid + 1)
       call send_ncover(BCS_blocks, myid + 1)
       ! Initialise the routines to calculate ion-ion interactions
       if(flag_neutral_atom) then
          call make_cs(inode-1,ion_ion_cutoff,ion_ion_CS,parts,bundle,&
               ni_in_cell, x_atom_cell,y_atom_cell,z_atom_cell)
       else
          call make_cs(inode-1,ewald_real_cutoff,ion_ion_CS,parts,bundle,&
               ni_in_cell, x_atom_cell,y_atom_cell,z_atom_cell)
       end if
       if (flag_dft_d2) call make_cs(inode-1, r_dft_d2, D2_CS, parts, bundle, ni_in_cell, &
               x_atom_cell, y_atom_cell, z_atom_cell)
    else
       call updateMembers_cs
       call deallocate_distribute_atom
       ! Reallocate and find new indices
       call distribute_atoms(inode,ionode)
       call immi(parts,bundle,BCS_parts,myid+1)
       if (flag_XLBOMD) call immi_XL(parts,bundle,BCS_parts,myid+1)
    end if
    ! Write out new coordinates
    append_coords_bkup = append_coords
    append_coords = .false.
    call write_atomic_positions('coord_next.dat',trim(pdb_template))
    append_coords = append_coords_bkup

    ! Update glob2node
    if (inode.EQ.ionode) call make_glob2node
    call gcopy(glob2node,ni_in_cell)
    ! Reallocate for blip grid
    call set_blipgrid(myid, RadiusAtomf, core_radius)
    call set_bucket(inode-1)
    call associate_fn_on_grid
    if(flag_diagonalisation) call init_scalapack_format
    call stop_print_timer(tmr_l_tmp2,"matrix reindexing",IPRINT_TIME_THRES2)
    if (flag_Becke_weights) call build_Becke_weights
    call stop_print_timer(tmr_l_tmp1,"indices update",IPRINT_TIME_THRES2)

    return
  end subroutine updateIndices3
  !!*****

  ! --------------------------------------------------------------------
  ! Subroutine update_H
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/update_H *
  !!
  !!  NAME 
  !!   update_H
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Updates various quantities when atoms move: blips, S, n(r), H
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   08:24, 2003/02/05 dave
  !!  MODIFICATION HISTORY
  !!   15:04, 27/02/2003 drb & tm 
  !!    Added call to set_density for Harris-Foulkes type calculations;
  !!    completely sorted out charge density questions
  !!   18:28, 2003/02/27 dave
  !!    Added call to deallocate Tm pseudopotential
  !!   10:25, 06/03/2003 drb 
  !!    Corrected Tm pseudo updating (alloc/dealloc not needed)
  !!   13:12, 22/10/2003 mjg & drb 
  !!    Added old/new ewald calls
  !!   12:13  31/03/2011 M.Arita
  !!    Added the statement to recall sbrt: set_density_pcc for NSC cg calculations
  !!   2011/09/29 16:50 M. Arita
  !!    Dispersions are calculated with a new set of atoms
  !!   2011/11/28 L.Tong
  !!    Added spin polarisation
  !!   2011/12/09 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2012/03/26 L.Tong
  !!   - Changed spin implementation
  !!   2013/08/26 M.Arita
  !!   - Added call for get_electronic_density to calculate charge density
  !!     from L-matrix
  !!   2013/12/02 M.Arita
  !!   - Corrected calls to generate charge density
  !!   - Added call for get_initiaL_XL to calculate an initial guess for
  !!     L-matrix when XL-BOMD applies
  !!   2015/11/24 08:32 dave
  !!    Removed old ewald calls
  !!   2016/01/13 08:20 dave
  !!    Added call to set_atomic_density for reset, non-SCF and NA
  !!   2016/01/29 15:01 dave
  !!    Removed prefix for ewald call
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2017/01/18 10:00 nakata
  !!    Added call to initail_SFcoeff_SSSF/MSSF
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2017/09/06 19:00 nakata
  !!    Changed to call set_atomic_density with ".true." before calling initial_SFcoeff
  !!   2017/11/07 10:02 dave
  !!    Added scaling of electron density after atom move
  !!   2017/11/17 14:51 dave
  !!    Bug fix: removed erroneous spin scaling on electron density
  !!  SOURCE
  !!
  subroutine update_H(fixed_potential)

    use numbers
    use logicals
    use timer_module    
    use S_matrix_module,        only: get_S_matrix
    use H_matrix_module,        only: get_H_matrix
    use mult_module,            only: LNV_matrix_multiply, matrix_scale, matrix_transpose, &
                                      matSFcoeff,matSFcoeff_tran 
    use ion_electrostatic,      only: ewald, screened_ion_interaction
    use pseudopotential_data,   only: init_pseudo
    use pseudo_tm_module,       only: set_tm_pseudo
    use pseudopotential_common, only: pseudo_type, OLDPS, SIESTA,      &
                                      STATE, ABINIT, core_correction
    use global_module,          only: iprint_MD, flag_self_consistent, &
                                      IPRINT_TIME_THRES2,              &
                                      flag_pcc_global, flag_dft_d2,    &
                                      nspin, io_lun,                   &
                                      flag_mix_L_SC_min, flag_XLBOMD,  &
                                      flag_reset_dens_on_atom_move,    &
                                      flag_LmatrixReuse,               &
                                      flag_neutral_atom,               &
                                      atomf, sf, nspin_SF, flag_LFD,   &
                                      flag_SFcoeffReuse, flag_diagonalisation, &
                                      ne_spin_in_cell,                 &
                                      ne_in_cell, spin_factor,         &
                                      ni_in_cell, area_moveatoms
    use density_module,         only: set_atomic_density,              &
                                      density, set_density_pcc,        &
                                      get_electronic_density
    use GenComms,               only: cq_abort, inode, ionode
    use maxima_module,          only: maxngrid
    use DFT_D2,                 only: dispersion_D2
    use functions_on_grid,      ONLY: atomfns, H_on_atomfns
    use XLBOMD_module,          ONLY: get_initialL_XL
    use multisiteSF_module,     ONLY: initial_SFcoeff, flag_LFD_MD_UseAtomicDensity
    use memory_module,          only: reg_alloc_mem, type_dbl

    implicit none

    ! Shared variables needed by get_H_matrix for now (!)
    logical :: fixed_potential

    ! Local variables
    type(cq_timer) :: tmr_l_tmp1
    real(double), dimension(nspin) :: electrons, energy_tmp
    real(double) :: scale
    integer :: spin_SF, spin, stat

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    ! (0) Pseudopotentials: choose correct form
    select case (pseudo_type)
    case (OLDPS)
       call init_pseudo(core_correction)
    case (SIESTA)
       call set_tm_pseudo
    case (ABINIT)
       call set_tm_pseudo
    end select
    ! (0) Prepare SF-PAO coefficients for contracted SFs
    if (atomf.ne.sf) then
       if (flag_SFcoeffReuse) then
       ! Use the coefficients in the previous step   
       ! SF coeffs are already updated before calling update_H, but we need its transpose
         do spin_SF = 1,nspin_SF
          call matrix_transpose(matSFcoeff(spin_SF), matSFcoeff_tran(spin_SF))
         enddo
       else
        do spin_SF = 1, nspin_SF
          call matrix_scale(zero,matSFcoeff(spin_SF))
        enddo
       ! Make SF coefficients newly
          ! Use the atomic density if flag_LFD_MD_UseAtomicDensity=T,
          ! otherwise, use the density in the previous step
          if (flag_LFD_MD_UseAtomicDensity) call set_atomic_density(.true.)
          call initial_SFcoeff(.true., .true., fixed_potential, .true.)
       endif
    endif
    ! (1) Get S matrix (includes blip-to-grid transform)
    if (flag_LFD .and. .not.flag_SFcoeffReuse) then
       ! Spao was already made in sub:initial_SFcoeff
       call get_S_matrix(inode, ionode, build_AtomF_matrix=.false.)
    else
       call get_S_matrix(inode, ionode)
    endif
    ! (2) Make L
    if (flag_XLBOMD) call get_initialL_XL()
    ! (3) get K matrix if O(N)
    if (.not. flag_diagonalisation) then
       call LNV_matrix_multiply(electrons, energy_tmp, doK, dontM1, &
                                dontM2, dontM3, dontM4, dontphi,    &
                                dontE)
    end if
    ! (4) core correction?
    ! (5) Find the Ewald energy for the initial set of atoms
    if(flag_neutral_atom) then
       call screened_ion_interaction
    else
       call ewald
    end if
    ! (6) Find the dispersion for the initial set of atoms
    if (flag_dft_d2) call dispersion_D2
    ! Now we call set_density if (a) we have non-SCF and atomic densities or
    ! (b) the flag_reset_dens_on_atom_move is set
        !2019Dec27 tsuyoshi flag_no_atomic_densities was removed
        ! For (non-SCF or reset density by atomic densities)
    if(((.NOT. flag_self_consistent) .AND. (.NOT. flag_mix_L_SC_min))&
         .OR.flag_reset_dens_on_atom_move) then
        call set_atomic_density(.true.)
    ! For SCF-O(N) calculations
    elseif (.NOT.flag_diagonalisation) then
       if (flag_self_consistent .OR. flag_mix_L_SC_min) then
          if(flag_neutral_atom .and. .not.flag_LFD_MD_UseAtomicDensity) call set_atomic_density(.false.)
          if(flag_LmatrixReuse) then
             if (inode.EQ.ionode.AND.iprint_MD>2) write (io_lun,*) "update_H: Get charge density from L-matrix"
             call get_electronic_density(density,electrons,atomfns,H_on_atomfns(1), &
                  inode,ionode,maxngrid)
            do spin=1,nspin
               scale = ne_spin_in_cell(spin)/electrons(spin)
               if(abs(scale-one)<threshold_resetCD .AND. abs(scale-one)>RD_ERR) then
                  density(:,spin) = density(:,spin)*scale
               else if (abs(scale-one)>=threshold_resetCD) then
                  if(inode == ionode) write(io_lun,*) &
                       ' WARNING! in update_H: constructed charge density is strange... ; scale = ', &
                       scale, threshold_resetCD
                  call set_atomic_density(.true.)
                  exit
               end if
            end do
             ! if flag_LFD=T, update SF-PAO coefficients with the obtained density
             ! and update S with the coefficients
             !ORI if (flag_LFD) then
             if (flag_LFD .and. .not.flag_SFcoeffReuse) then
                call initial_SFcoeff(.false., .true., fixed_potential, .false.)
                call get_S_matrix(inode, ionode, build_AtomF_matrix=.false.)
             endif
          end if
       endif
    else if(flag_diagonalisation.AND.flag_neutral_atom) then
       if (.not.flag_LFD_MD_UseAtomicDensity) call set_atomic_density(.false.)
    end if
    ! If we have read K and are predicting density from it, then rebuild
    if(flag_diagonalisation.AND.flag_LmatrixReuse) then
       call get_electronic_density(density,electrons,atomfns,H_on_atomfns(1), &
            inode,ionode,maxngrid)
       do spin=1,nspin
          scale = ne_spin_in_cell(spin)/electrons(spin)
          if(abs(scale-one)<threshold_resetCD .AND.abs(scale-one)>RD_ERR) then
             density(:,spin) = density(:,spin)*scale
          else if (abs(scale-one)>=threshold_resetCD) then
             if(inode == ionode) write(io_lun,*)  ' WARNING2! in update_H: &
                  & constructed charge density is strange... ; scale,threshold = ', &
                  scale, threshold_resetCD
             call set_atomic_density(.true.)
             exit
          end if
       end do
    end if
    if (flag_pcc_global) call set_density_pcc()
    ! (7) Now generate a new H matrix, including a new charge density
    if (flag_LFD .and. .not.flag_SFcoeffReuse) then
       ! Hpao was already made in sub:initial_SFcoeff
       call get_H_matrix(.false., fixed_potential, electrons, density, &
                         maxngrid, build_AtomF_matrix=.false.)
    else
       call get_H_matrix(.true., fixed_potential, electrons, density, &
                         maxngrid)
    endif
    call stop_print_timer(tmr_l_tmp1, "update_H", IPRINT_TIME_THRES2)
    return
  end subroutine update_H
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine checkBonds
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/checkBonds *
  !!
  !!  NAME 
  !!   checkBonds
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Checks to see if there are new H range interactions
  !!
  !!   Relies heavily on the methodology of get_naba: loops over GCS
  !!   partitions and atoms, checking for atoms within range and
  !!   compares to known atoms.  Note various things:
  !!
  !!    i)   We can rely on the order of the GCS not changing
  !!    ii)  We only want to compare the partition and sequence
  !!         numbers of the atoms, not their separation
  !!    iii) As an extra check, we compare the number of neighbours
  !!         of primary set atoms
  !!
  !!  INPUTS
  !!   logiccal :: newAtom Flags if a new atom is found
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   09:19, 2003/01/29 dave
  !!  MODIFICATION HISTORY
  !!   2008/07/18 ast
  !!     Added timers
  !!  SOURCE
  !!
  subroutine checkBonds(newAtom, prim, gcs, amat, partsproc, rcut)

    use datatypes
    use basic_types,   only: primary_set, cover_set
    use matrix_module, only: matrix
    use global_module, only: IPRINT_TIME_THRES2
    use timer_module

    implicit none

    ! Passed variables
    logical, intent(out) :: newAtom
    integer           :: partsproc
    type(primary_set) :: prim
    type(cover_set)   :: gcs
    real(double)      :: rcut
    type(matrix)      :: amat(partsproc)

    ! Local variables
    real(double)   :: rcutsq, dx, dy, dz
    real(double)   :: tol = 1.0e-8_double
    integer        :: inp, nn, j, np, ni, ist, n_nab
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    rcutsq = rcut*rcut
    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    newAtom = .false.
    do nn=1,prim%groups_on_node ! Partitions in primary set
       if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
          do j=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
             n_nab = 0
             do np=1,gcs%ng_cover  ! Loop over partitions in GCS
                if(gcs%n_ing_cover(np).gt.0) then  ! Are there atoms ?
                   do ni=1,gcs%n_ing_cover(np)
                      dx=gcs%xcover(gcs%icover_ibeg(np)+ni-1)-prim%xprim(inp)
                      dy=gcs%ycover(gcs%icover_ibeg(np)+ni-1)-prim%yprim(inp)
                      dz=gcs%zcover(gcs%icover_ibeg(np)+ni-1)-prim%zprim(inp)
                      if(dx*dx+dy*dy+dz*dz.lt.rcutsq-tol) then ! Neighbour
                         n_nab = n_nab + 1
                         ist = amat(nn)%i_acc(j)+n_nab-1
                         if(ist>amat(nn)%part_nabs) then
                            newAtom = .true.
                         else
                            ! Check - is this one we've seen before ?
                            if (np /= amat(nn)%i_part(ist) .and. &
                                ni /= amat(nn)%i_seq(ist)) &
                                newAtom = .true.
                         end if
                      end if
                   end do ! End n_inp_cover
                end if
             end do ! End np_cover
             inp = inp + 1  ! Indexes primary-set atoms
             if (n_nab /= amat(nn)%n_nab(j)) newAtom = .true.
          end do ! End prim%nm_nodgroup
       end if ! End if(prim%nm_nodgroup>0)
    end do ! End part_on_node
    call stop_print_timer(tmr_l_tmp1,"checking bonds",IPRINT_TIME_THRES2)
    return
  end subroutine checkBonds
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine primary_update
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/primary_update *
  !!
  !!  NAME 
  !!   primary_update
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Updates the atomic positions in primary set after atom 
  !!   movement
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   07:41, 2003/01/29 dave (from ParaDens)
  !!  MODIFICATION HISTORY
  !!   2008/07/18 ast
  !!     Added timers
  !!  SOURCE
  !!
  subroutine primary_update(x_position, y_position, z_position, prim, &
                            groups, myid)

    use datatypes
    use basic_types
    use global_module, only: ni_in_cell, rcellx, rcelly, rcellz, &
                             IPRINT_TIME_THRES3
    use timer_module

    implicit none

    ! Passed variables
    real(double), dimension(ni_in_cell) :: x_position,y_position,z_position
    type(primary_set) :: prim
    type(group_set)   :: groups
    integer           :: myid

    ! Local variables
    integer        :: ng, ind_group, nx, ny, nz, nx1, ny1, nz1, nnd, &
                      n_prim, ni
    real(double)   :: dcellx, dcelly, dcellz
    real(double)   :: xadd, yadd, zadd
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_std_indexing)    ! NOTE: This will be annotated in area 8
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    n_prim = 0
    nnd = myid+1
    dcellx=rcellx/groups%ngcellx
    dcelly=rcelly/groups%ngcelly
    dcellz=rcellz/groups%ngcellz
    do ng = 1,groups%ng_on_node(nnd)
       ind_group=groups%ngnode(groups%inode_beg(nnd)+ng-1)
       nx=1+(ind_group-1)/(groups%ngcelly*groups%ngcellz)
       ny=1+(ind_group-1-(nx-1)*groups%ngcelly*groups%ngcellz)/groups%ngcellz
       nz=ind_group-(nx-1)*groups%ngcelly*groups%ngcellz-(ny-1)*groups%ngcellz
       nx1=prim%nx_origin+prim%idisp_primx(ng)
       ny1=prim%ny_origin+prim%idisp_primy(ng)
       nz1=prim%nz_origin+prim%idisp_primz(ng)
       xadd=real(nx1-nx,double)*dcellx
       yadd=real(ny1-ny,double)*dcelly
       zadd=real(nz1-nz,double)*dcellz
       if(prim%nm_nodgroup(ng).gt.0) then
          do ni=1,prim%nm_nodgroup(ng)
             n_prim = n_prim + 1
             prim%xprim(n_prim)= &
                  x_position(groups%icell_beg(ind_group)+ni-1)+xadd
             prim%yprim(n_prim)= &
                  y_position(groups%icell_beg(ind_group)+ni-1)+yadd
             prim%zprim(n_prim)= &
                  z_position(groups%icell_beg(ind_group)+ni-1)+zadd
          end do ! Atoms in partition
       end if ! If atoms in partition
    end do ! Groups on node
    call stop_print_timer(tmr_l_tmp1,"primary update",IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_indexing)
  end subroutine primary_update
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine cover_update
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/cover_update *
  !!
  !!  NAME 
  !!   cover_update
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Updates the atomic positions in cover set after atom 
  !!   movement
  !!
  !!   The details of how atoms are offset in covering sets
  !!   are rather horrific (18 certificate certainly) and
  !!   are described in graphic detail in a set of notes by
  !!   Dave Bowler, referred to in cover_module - see there
  !!   for details.
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   07:41, 2003/01/29 dave (from ParaDens)
  !!  MODIFICATION HISTORY
  !!   2008/05/25
  !!    Added timers
  !!   2019/11/04 15:14 dave
  !!    Replace call to indexx with call to heapsort_integer_index
  !!  SOURCE
  !!
  subroutine cover_update(x_position, y_position, z_position, set, groups)

    use datatypes
    use basic_types
    use global_module, only: ni_in_cell, rcellx, rcelly, rcellz, &
                             IPRINT_TIME_THRES3
    use functions,  only: heapsort_integer_index
    use GenComms,      only: cq_abort, myid
    use timer_module

    implicit none

    ! Passed variables
    real(double), dimension(ni_in_cell) :: x_position,y_position,z_position
    type(cover_set) :: set
    type(group_set) :: groups

    ! Local variables
    integer        :: ind_cover, nx, ny, nz, nx_o, ny_o, nz_o, ni
    integer        :: nsx,nsy,nsz,nqx,nqy,nqz,nmodx,nmody,nmodz
    integer        :: cover_part, ind_qart
    real(double)   :: dcellx, dcelly, dcellz
    real(double)   :: xadd, yadd, zadd
    type(cq_timer) :: tmr_l_tmp1

    integer :: nrepx(groups%mx_gedge)
    integer :: nrepy(groups%mx_gedge)
    integer :: nrepz(groups%mx_gedge)
    ! x,y,z numbering of CS groups (for CC labels)
    integer, allocatable, dimension(:)  :: nx_in_cover
    integer, allocatable, dimension(:)  :: ny_in_cover
    integer, allocatable, dimension(:)  :: nz_in_cover
    ! Variables for irreducible CS
    integer :: ind_min(groups%mx_gcell)
    integer :: ngcx_min(groups%mx_gcell)
    integer :: ngcy_min(groups%mx_gcell)
    integer :: ngcz_min(groups%mx_gcell)
    integer :: min_sort(groups%mx_gcell)
    integer :: noccx,nremx,minx,ngcx,noccy,nremy,miny,ngcy,noccz,nremz
    integer :: minz,ngcz,ng_in_min,ind,ino, stat, ng_in_cell
    integer :: nrx,nry,nrz, nnd

    call start_timer(tmr_std_indexing)    ! NOTE: This will be annotated in area 8
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    nnd = myid+1
    call start_timer(tmr_std_allocation)
    allocate(nx_in_cover(set%ng_cover), ny_in_cover(set%ng_cover), &
             nz_in_cover(set%ng_cover), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating nx_in_cover: ", set%ng_cover,stat)
    call stop_timer(tmr_std_allocation)
    ! Conversion factors from unit cell lengths->groups
    dcellx=rcellx/real(groups%ngcellx,double)
    dcelly=rcelly/real(groups%ngcelly,double)
    dcellz=rcellz/real(groups%ngcellz,double)
    ! Used in calculating offsets of groups in CS
    nmodx=((groups%ngcellx+set%nspanlx-1)/groups%ngcellx)*groups%ngcellx
    nmody=((groups%ngcelly+set%nspanly-1)/groups%ngcelly)*groups%ngcelly
    nmodz=((groups%ngcellz+set%nspanlz-1)/groups%ngcellz)*groups%ngcellz
    ! Origin of CS
    nx_o = set%nx_origin
    ny_o = set%ny_origin
    nz_o = set%nz_origin

    noccx=set%ncoverx/groups%ngcellx
    nremx=set%ncoverx-noccx*groups%ngcellx
    minx=min(set%ncoverx,groups%ngcellx)
    if (minx>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in x-edge', minx)
    end if
    do ngcx=1,minx
       if(ngcx<=nremx) then
          nrepx(ngcx)=noccx+1
       else
          nrepx(ngcx)=noccx
       end if
    end do
    ! ... y-direction
    noccy=set%ncovery/groups%ngcelly
    nremy=set%ncovery-noccy*groups%ngcelly
    miny=min(set%ncovery,groups%ngcelly)
    if (miny>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in y-edge',miny)
    end if
    do ngcy=1,miny
       if(ngcy<=nremy) then
          nrepy(ngcy)=noccy+1
       else
          nrepy(ngcy)=noccy
       end if
    end do
    ! ... z-direction
    noccz=set%ncoverz/groups%ngcellz
    nremz=set%ncoverz-noccz*groups%ngcellz
    minz=min(set%ncoverz,groups%ngcellz)
    if(minz>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in z-edge', minz)
    end if
    do ngcz=1,minz
       if(ngcz<=nremz) then
          nrepz(ngcz)=noccz+1
       else
          nrepz(ngcz)=noccz
       end if
    end do
    ! go over groups in GCS periodic-irreducible set, calculating
    ! simulation-cell (node-order, home-start) label of each 
    ng_in_cell = groups%ngcellx*groups%ngcelly*groups%ngcellz
    ng_in_min = minx*miny*minz
    ind=0
    do ngcx=1,minx
       do ngcy=1,miny
          do ngcz=1,minz
             ind=ind+1
             nqx=1+mod(nx_o+ngcx-set%nspanlx-2+nmodx,groups%ngcellx)
             nqy=1+mod(ny_o+ngcy-set%nspanly-2+nmody,groups%ngcelly)
             nqz=1+mod(nz_o+ngcz-set%nspanlz-2+nmodz,groups%ngcellz)
             ind_qart = (nqx-1) * groups%ngcelly * groups%ngcellz + &
                        (nqy-1) * groups%ngcellz+nqz
             ino=groups%inv_ngnode(ind_qart)
             ind_min(ind)=1+mod(ino-groups%inode_beg(nnd)+ ng_in_cell,ng_in_cell)
             ngcx_min(ind)=ngcx
             ngcy_min(ind)=ngcy
             ngcz_min(ind)=ngcz
          enddo
       enddo
    enddo
    ! sort minimum CS by nodes 
    call heapsort_integer_index(ng_in_min,ind_min,min_sort)
    ! go over all GCS groups in node-periodic-grouped order 
    ind_cover=0
    do ind=1,ng_in_min
       ngcx=ngcx_min(min_sort(ind))
       ngcy=ngcy_min(min_sort(ind))
       ngcz=ngcz_min(min_sort(ind))
       do nrx=1,nrepx(ngcx)
          do nry=1,nrepy(ngcy)
             do nrz=1,nrepz(ngcz)
                ind_cover=ind_cover+1
                nx_in_cover(ind_cover)=ngcx-1-set%nspanlx+(nrx-1)*groups%ngcellx
                ny_in_cover(ind_cover)=ngcy-1-set%nspanly+(nry-1)*groups%ngcelly
                nz_in_cover(ind_cover)=ngcz-1-set%nspanlz+(nrz-1)*groups%ngcellz
             enddo
          enddo
       enddo
    enddo


    do ind_cover=1,set%ng_cover
       !cover_part = set%lab_cover(ind_cover)
       !nx=1+(cover_part-1)/(set%ncovery*set%ncoverz)
       !ny=1+(cover_part-1-(nx-1)*set%ncovery*&
       !     set%ncoverz)/set%ncoverz
       !nz=cover_part-(nx-1)*set%ncovery*set%ncoverz-&
       !     (ny-1)*set%ncoverz
       !nsx=nx-1-set%nspanlx
       !nsy=ny-1-set%nspanly
       !nsz=nz-1-set%nspanlz
       nsx=nx_in_cover(ind_cover)
       nsy=ny_in_cover(ind_cover)
       nsz=nz_in_cover(ind_cover)
       nqx=1+mod(nx_o+nsx+nmodx-1,groups%ngcellx)
       nqy=1+mod(ny_o+nsy+nmody-1,groups%ngcelly)
       nqz=1+mod(nz_o+nsz+nmodz-1,groups%ngcellz)
       xadd=(nx_o+nsx-nqx)*dcellx
       yadd=(ny_o+nsy-nqy)*dcelly
       zadd=(nz_o+nsz-nqz)*dcellz
       !ind_qart=(nqx-1)*groups%ngcelly*groups%ngcellz+&
       !     (nqy-1)*groups%ngcellz+nqz
       ind_qart= set%lab_cell(ind_cover)
       do ni=1,groups%nm_group(ind_qart)
          set%xcover(set%icover_ibeg(ind_cover)+ni-1)= &
               x_position(groups%icell_beg(ind_qart)+ni-1)+xadd
          set%ycover(set%icover_ibeg(ind_cover)+ni-1)= &
               y_position(groups%icell_beg(ind_qart)+ni-1)+yadd
          set%zcover(set%icover_ibeg(ind_cover)+ni-1)= &
               z_position(groups%icell_beg(ind_qart)+ni-1)+zadd
       end do
    end do
    call start_timer(tmr_std_allocation)
    deallocate(nx_in_cover,ny_in_cover,nz_in_cover,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating nx_in_cover: ",set%ng_cover,stat)
    call stop_timer(tmr_std_allocation)
    call stop_print_timer(tmr_l_tmp1,"cover update",IPRINT_TIME_THRES3)
    call stop_timer(tmr_std_indexing)
  end subroutine cover_update
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine update_atom_coord
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/update_atom_coord *
  !!  
  !!  NAME 
  !!   update_atom_coord
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Updates the atomic positions (atom_coord) in global_module after atom 
  !!   movement 
  !!  
  !!  INPUTS
  !!  
  !!  
  !!  USES
  !!   global_module
  !!  AUTHOR
  !!   T. Miyazaki
  !!  CREATION DATE
  !!   27 Aug 2003
  !!  MODIFICATION HISTORY
  !!   2008/05/25
  !!    Added timers
  !!   2018/07/11 12:11 dave
  !!    Changed iprint level for output of partition boundary crossing to > 3
  !!  SOURCE
  !!
  subroutine update_atom_coord
    
    use datatypes
    use global_module, only: x_atom_cell, y_atom_cell, z_atom_cell,   &
                             id_glob, atom_coord, ni_in_cell, io_lun, &
                             iprint_MD, IPRINT_TIME_THRES2
    use dimens,        only: r_super_x, r_super_y, r_super_z
    use group_module,  only: parts
    use timer_module
    use GenComms, only: inode, ionode
    
    implicit none

    integer        :: ni, id_global
    real(double)   :: dx, dy, dz
    type(cq_timer) :: tmr_l_tmp1
    
    call start_timer(tmr_std_indexing)    ! NOTE: This will be annotated in area 8
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    dx = r_super_x / parts%ngcellx
    dy = r_super_y / parts%ngcelly
    dz = r_super_z / parts%ngcellz

    do ni = 1, ni_in_cell
       id_global = id_glob(ni)
       if (iprint_MD > 3) then
          if (floor(atom_coord(1,id_global)/dx) /= &
              floor(x_atom_cell(ni)/dx)) then
             write (io_lun, *) inode, id_global, &
                               ' Partition boundary crossed in x ! ', &
                               dx, atom_coord(1,id_global), x_atom_cell(ni)
          end if
          if (floor(atom_coord(2,id_global)/dy) /= &
              floor(y_atom_cell(ni)/dy)) then
             write (io_lun, *) inode, id_global, &
                               'Partition boundary crossed in y ! ', &
                               dy, atom_coord(2,id_global), y_atom_cell(ni)
          end if
          if (floor(atom_coord(3,id_global)/dz) /= &
              floor(z_atom_cell(ni)/dz)) then
             write (io_lun, *) inode, id_global, &
                               'Partition boundary crossed in z ! ', &
                               dz, atom_coord(3,id_global), z_atom_cell(ni)
          end if
       end if
       atom_coord(1,id_global)= x_atom_cell(ni)
       atom_coord(2,id_global)= y_atom_cell(ni)
       atom_coord(3,id_global)= z_atom_cell(ni)
    end do
    call stop_print_timer(tmr_l_tmp1, "coordinates update", &
                          IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_indexing)
    return
  end subroutine update_atom_coord
  !!***

  !!****f* move_atoms/update_r_atom_cell *
  !!  
  !!  NAME 
  !!   update_r_atom_cell
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Updates the x/y/z_atom_cell in global_module after atom movement
  !!   Adapted from update_atom_coord
  !!
  !!  INPUTS
  !!
  !!  USES
  !!   global_module
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/06/04
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine update_r_atom_cell
    
    use datatypes
    use global_module, only: x_atom_cell, y_atom_cell, z_atom_cell,   &
                             id_glob, atom_coord, ni_in_cell, io_lun, &
                             iprint_MD, IPRINT_TIME_THRES2, id_glob_inv
    use dimens,        only: r_super_x, r_super_y, r_super_z
    use group_module,  only: parts
    use timer_module
    use GenComms, only: inode, ionode
    
    implicit none

    integer        :: ni, id_global
    real(double)   :: dx, dy, dz
    type(cq_timer) :: tmr_l_tmp1
    
    call start_timer(tmr_std_indexing)    ! NOTE: This will be annotated in area 8
    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    dx = r_super_x / parts%ngcellx
    dy = r_super_y / parts%ngcelly
    dz = r_super_z / parts%ngcellz

    do id_global = 1, ni_in_cell
       ni = id_glob_inv(id_global)
       if (iprint_MD > 3) then
          if (floor(atom_coord(1,id_global)/dx) /= &
              floor(x_atom_cell(ni)/dx)) then
             write (io_lun, *) inode, id_global, &
                               ' Partition boundary crossed in x ! ', &
                               dx, atom_coord(1,id_global), x_atom_cell(ni)
          end if
          if (floor(atom_coord(2,id_global)/dy) /= &
              floor(y_atom_cell(ni)/dy)) then
             write (io_lun, *) inode, id_global, &
                               'Partition boundary crossed in y ! ', &
                               dy, atom_coord(2,id_global), y_atom_cell(ni)
          end if
          if (floor(atom_coord(3,id_global)/dz) /= &
              floor(z_atom_cell(ni)/dz)) then
             write (io_lun, *) inode, id_global, &
                               'Partition boundary crossed in z ! ', &
                               dz, atom_coord(3,id_global), z_atom_cell(ni)
          end if
       end if
       x_atom_cell(ni) = atom_coord(1,id_global)
       y_atom_cell(ni) = atom_coord(2,id_global)
       z_atom_cell(ni) = atom_coord(3,id_global)
    end do
    call stop_print_timer(tmr_l_tmp1, "coordinates update", &
                          IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_indexing)
    return
  end subroutine update_r_atom_cell
  !!***

  !!****f*  move_atoms/init_velocity *
  !!
  !!  NAME 
  !!   init_velocity
  !!  USAGE
  !!   
  !!  PURPOSE
  !!   Initialise ionic velocities for MD via normal distribution
  !!  INPUTS
  !!   ni_in_cell : no. of atoms in the cell
  !!   velocity(3, ni_in_cell) : velocity in (fs * Ha/bohr)/amu unit
  !!   temp_ion  : temperature for atoms (ions)
  !!  USES
  !!   datatypes, numbers, species_module, global_module
  !!  AUTHOR
  !!   T. Miyazaki
  !!  CREATION DATE
  !!   2010/6/30 
  !!  MODIFICATION HISTORY
  !!   2019/05/21 zamaan
  !!    Replaced old rng calls with new one from rng module
  !!   2019/05/22 14:40 dave & tsuyoshi
  !!    Moved ionode criterion for generation of velocities from init_ensemble
  !!   2019/05/23 zamaan
  !!    Zeroed COM velocity after initialisation
  !!   2020/07/28 tsuyoshi
  !!    Zeroed velocity for the fixed degree of freedom
  !!  SOURCE
  !!
  subroutine init_velocity(ni_in_cell, temp, velocity)

    use datatypes,      only: double
    use numbers,        only: three,two,twopi, zero, one, RD_ERR, half, three_halves
    use species_module, only: species, mass
    use global_module,  only: id_glob_inv, flag_move_atom, species_glob, &
                              iprint_MD, flag_FixCOM
    use GenComms,       only: cq_abort, inode, ionode, gcopy
    use rng,            only: type_rng

    implicit none

    integer,intent(in) :: ni_in_cell
    real(double),intent(in) :: temp
    real(double),intent(out):: velocity(3,ni_in_cell)
    integer :: dir, ia, iglob
    real(double) :: xx, yy, zz, u0, ux, uy, uz, v0
    real(double) :: massa, KE
    integer :: speca

    type(type_rng) :: myrng

    if (inode == ionode) then
       KE = zero
       if(iprint_MD > 2) write(io_lun,'(2x,a,f16.8)') ' Welcome to init_velocity, fac = ',fac
       velocity(:,:) = zero
       call myrng%init_rng
       call myrng%init_normal(one, zero)

       !  We would like to use the order of global labelling in the following.
       !(since we use random numbers, the order of atoms is probably relevant
       ! if we want to have a same distribution of velocities as in other codes.)

       do iglob=1,ni_in_cell
          ia= id_glob_inv(iglob)
          speca= species(ia)
          massa= mass(speca)
          if(ia < 1 .or. ia > ni_in_cell) &
               call cq_abort('ERROR in init_velocity : ia,iglob ',ia,iglob)

          ! -- (Important Notes) ----
          ! it is tricky, but velocity is in the unit, bohr/fs, transforming from
          ! (fs * Har/bohr)/ amu, with the factor (fac) defined in the beginning of 
          ! this module.  This factor comes from that v is calculated as (dt*F/mass), 
          ! and we want to express dt in femtosecond, force in Hartree/bohr, 
          ! and m in atomic mass units. 
          ! (it should be equivalent to express dt and m in atomic units, I think.)
          ! Kinetic Energy is calculated as m/2*v^2 *fac in Hartree unit, and
          ! Positions are calculated as v*dt in bohr unit. (m in amu, dt in fs)
          v0 = sqrt(temp*fac_Kelvin2Hartree/(massa*fac)) 
          do dir=1,3
             if(flag_move_atom(dir,iglob)) then
              u0 = myrng%rng_normal()
             else
              u0 = zero
             endif
             velocity(dir,ia) = v0 * u0
             KE = KE + half * massa * fac * velocity(dir,ia)**2
          end do
       enddo
       if (flag_FixCOM) call zero_COM_velocity(velocity)
       KE = KE*three_halves
       KE = KE/(real(ni_in_cell,double)*fac_Kelvin2Hartree)
       if(iprint_MD > 2) write(io_lun,*) ' init_velocity: Kinetic Energy in K = ',KE
    end if
    call gcopy(velocity, 3, ni_in_cell)
    return
  end subroutine init_velocity
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine wrap_xyz_atom_cell
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/wrap_xyz_atom_cell *
  !!  
  !!  NAME 
  !!   wrap_xyz_atom_cell
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Wrapping atomic positions ("x,y,z_atom_cell": bohr units, partition labelling)
  !!   into the unit cell. This is necessary for 'partition' technology in Conquest.
  !!   In order to have a common distribution of atoms into partitions, we need
  !!   to use the same shift_in_bohr as used in atom2part or allatom2part.
  !!   This is important for the atoms on the bondary of partitions.
  !!  INPUTS
  !!  
  !!  USES
  !!   global_module
  !!  AUTHOR
  !!   M.Arita & T.Miyazaki
  !!  CREATION DATE
  !!   2013/07/01
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine wrap_xyz_atom_cell
    
    use datatypes
    use global_module, only: x_atom_cell, y_atom_cell, z_atom_cell,   &
                             shift_in_bohr, ni_in_cell, io_lun, iprint_MD
    use dimens,        only: r_super_x, r_super_y, r_super_z

    implicit none
    integer        :: atom
    real(double)   :: eps

    eps=shift_in_bohr
    do atom = 1, ni_in_cell
      x_atom_cell(atom) = x_atom_cell(atom) - floor((x_atom_cell(atom)+eps)/r_super_x)*r_super_x
      y_atom_cell(atom) = y_atom_cell(atom) - floor((y_atom_cell(atom)+eps)/r_super_y)*r_super_y
      z_atom_cell(atom) = z_atom_cell(atom) - floor((z_atom_cell(atom)+eps)/r_super_z)*r_super_z
    enddo
      
    return
  end subroutine wrap_xyz_atom_cell
  !!***

  ! --------------------------------------------------------------------
  ! Subroutine calculate_kinetic_energy
  ! --------------------------------------------------------------------
  
  !!****f* move_atoms/calculate_kinetic_energy *
  !!  NAME 
  !!   calculate_kinetic_energy
  !!  USAGE
  !!   call calculate_kinetic_energy(v,KE)
  !!  PURPOSE
  !!   Calculates the ionic kinetic energy
  !!  INPUTS
  !!   real(double), v : particle velocity
  !!   real(double), KE: kinetic energy
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2014/02/03
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine calculate_kinetic_energy(v,KE)
    ! Module usage
    use datatypes
    use numbers, ONLY: zero,half
    use global_module, ONLY: ni_in_cell
    use species_module, ONLY: species,mass

    implicit none
    ! passed variables
    real(double),dimension(3,ni_in_cell),intent(in) :: v
    real(double),intent(out) :: KE
    ! local variables
    integer :: atom,k,speca
    real(double) :: massa

    KE = zero
    do atom = 1, ni_in_cell
      speca = species(atom)
      massa = mass(speca)*fac
      do k = 1, 3
        KE = KE + massa*v(k,atom)*v(k,atom)
      enddo
    enddo
    KE = half*KE

    return
  end subroutine calculate_kinetic_energy
  !!***

  !!****f* move_atoms/zero_COM_velocity *
  !!  NAME 
  !!   zero_COM_velocity
  !!  USAGE
  !!   call zero_COM_velocity(v)
  !!  PURPOSE
  !!   Fixes the centre-of-mass of the system
  !!  INPUT
  !!   real(double), v : particle velocity
  !!  OUTPUT
  !!   real(double), v : particle velocity subtracted by COM velocity
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2014/02/03
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine zero_COM_velocity(v)
    ! Module usage
    use datatypes
    use numbers, ONLY: zero
    use global_module, ONLY: ni_in_cell
    use species_module, ONLY: species,mass

    implicit none
    ! passed variable
    real(double),dimension(3,ni_in_cell),intent(inout) :: v
    ! local variables
    integer :: atom,k,speca
    real(double) :: massa,M
    real(double),dimension(3) :: COMv

    ! Calculates centre-of-mass velocity
    M = zero
    COMv = zero
    do atom = 1, ni_in_cell
      speca = species(atom)
      massa = mass(speca)
      M = M + massa
      do k = 1, 3
        COMv(k) = COMv(k) + massa*v(k,atom)
      enddo
    enddo
    COMv = COMv / M

    ! Subtracts centre-of-mass velocity from particle velocity
    do atom = 1, ni_in_cell
      do k = 1, 3
        v(k,atom) = v(k,atom) - COMv(k)
      enddo
    enddo

    return
  end subroutine zero_COM_velocity
  !!***

  !!****f* move_atoms/check_move_atoms *
  !!  NAME 
  !!   check_move_atoms
  !!  USAGE
  !!   call check_move_atoms(flag_movable)
  !!  PURPOSE
  !!   Converts flag_move_atom to 1-D array
  !!  INPUT
  !!   logical,flag_movable: converted 1-D array to tell if atoms move
  !!  OUTPUT
  !!   logical,flag_movable: converted 1-D array to tell if atoms move
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2014/02/03
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine check_move_atoms(flag_movable)
    ! Module usage
    use global_module, ONLY: ni_in_cell,id_glob,flag_move_atom

    implicit none
    ! passed variable
    logical,dimension(3*ni_in_cell) :: flag_movable
    ! local variables
    integer :: atom,k,gatom,ibeg_atom

    ibeg_atom = 1
    do atom = 1, ni_in_cell
      gatom = id_glob(atom)
      do k = 1, 3
        flag_movable(ibeg_atom+k-1) = flag_move_atom(k,gatom)
      enddo
      ibeg_atom = ibeg_atom + 3
    enddo

    return
  end subroutine check_move_atoms
  !!***

  !!****f* move_atoms/update_cell_dims *
  !!  NAME
  !!   update_cell_dims
  !!  USAGE
  !!   call update_cell_dims()
  !!  PURPOSE
  !!   Updates the simulation cell dimensions subject to constraints on ratios.
  !!   e.g c/a = const. Upon a change in a, b or c, grids are updated and the
  !!   density is scaled.
  !!  INPUT
  !!  OUTPUT
  !!  AUTHOR
  !!  Jack Baker
  !!  David Bowler
  !!  CREATION DATE
  !!   30/05/17
  !!  MODIFICATION HISTORY
  !!   2020/04/24 08:15 dave
  !!    Bug fix for constrained ratios
  !!   2020/05/15 12:26 dave
  !!    Update to remove unnecessary code (a/c and c/a etc are the same)
  !!  SOURCE
  !!
  subroutine update_cell_dims(start_rcellx, start_rcelly, start_rcellz, &
                              search_dir_x, search_dir_y, search_dir_z, &
                              k, iter, search_dir_mean)
    use datatypes
    use numbers
    use units
    use global_module,      only: iprint_MD, x_atom_cell, y_atom_cell, z_atom_cell, &
         atom_coord, ni_in_cell, rcellx, rcelly, &
         rcellz, flag_self_consistent,           &
         flag_reset_dens_on_atom_move,           &
         IPRINT_TIME_THRES1, flag_pcc_global, &
         flag_diagonalisation, cell_constraint_flag
    use minimise,           only: get_E_and_F, sc_tolerance, L_tolerance, &
         n_L_iterations
    use GenComms,           only: my_barrier, myid, inode, ionode,        &
         cq_abort
    use io_module,          only: write_atomic_positions, pdb_template
    use density_module,     only: density, set_density_pcc
    use maxima_module,      only: maxngrid
    use multisiteSF_module, only: flag_LFD_NonSCF
    use timer_module
    use dimens, ONLY: r_super_x, r_super_y, r_super_z, &
         r_super_x_squared, r_super_y_squared, r_super_z_squared, volume, &
         grid_point_volume, one_over_grid_point_volume, n_grid_x, n_grid_y, n_grid_z
    use fft_module, ONLY: recip_vector, hartree_factor, i0
    use DiagModule, ONLY: kk, nkp
    use input_module,         only: leqi

    implicit none

    ! Passed variables
    real(double) :: start_rcellx, start_rcelly, start_rcellz,&
         search_dir_x, search_dir_y, search_dir_z, k, search_dir_mean
    integer iter
    ! Shared variables needed by get_E_and_F for now (!)
    logical           :: vary_mu, fixed_potential
    real(double)      :: total_energy
    character(len=40) :: output_file

    ! local variables
    real(double) :: orcellx, orcelly, orcellz, xvec, yvec, zvec, r2, scale
    integer :: i, j

    orcellx = rcellx
    orcelly = rcelly
    orcellz = rcellz
    ! Update based on constraints.
    ! none => Unconstrained case
    if (leqi(cell_constraint_flag, 'none')) then
        rcellx = start_rcellx + k * search_dir_x
        rcelly = start_rcelly + k * search_dir_y
        rcellz = start_rcellz + k * search_dir_z

    else if (leqi(cell_constraint_flag, 'volume')) then
        rcellx = start_rcellx + k * search_dir_mean
        rcelly = start_rcelly + k * search_dir_mean
        rcellz = start_rcellz + k * search_dir_mean

    ! Fix a single dimension?
    else if (leqi(cell_constraint_flag, 'a')) then
        rcelly = start_rcelly + k * search_dir_y
        rcellz = start_rcellz + k * search_dir_z
    else if (leqi(cell_constraint_flag, 'b')) then
        rcellx = start_rcellx + k * search_dir_x
        rcellz = start_rcellz + k * search_dir_z
    else if (leqi(cell_constraint_flag, 'c')) then
        rcelly = start_rcelly + k * search_dir_y
        rcellx = start_rcellx + k * search_dir_x

    ! Fix two dimensions?
    else if (leqi(cell_constraint_flag, 'c a') .or. leqi(cell_constraint_flag, 'a c')) then
        rcelly = start_rcelly + k * search_dir_y
    else if (leqi(cell_constraint_flag, 'a b') .or. leqi(cell_constraint_flag, 'b a')) then
        rcellz = start_rcellz + k * search_dir_z
    else if (leqi(cell_constraint_flag, 'b c') .or. leqi(cell_constraint_flag, 'c b')) then
        rcellx = start_rcellx + k * search_dir_x

    ! Fix a single ratio?
    else if (leqi(cell_constraint_flag, 'a/c') .OR. leqi(cell_constraint_flag, 'c/a')) then
       rcellx = start_rcellx + k * search_dir_x
       rcelly = start_rcelly + k * search_dir_y
       rcellz = start_rcellz + k * (start_rcellz/start_rcellx)*search_dir_x
    else if (leqi(cell_constraint_flag, 'a/b') .OR. leqi(cell_constraint_flag, 'b/a')) then
       rcellx = start_rcellx + k * search_dir_x
       rcelly = start_rcelly + k * (start_rcelly/start_rcellx)*search_dir_x
       rcellz = start_rcellz + k * search_dir_z
    else if (leqi(cell_constraint_flag, 'b/c') .OR. leqi(cell_constraint_flag, 'c/b')) then
       rcellx = start_rcellx + k * search_dir_x
       rcelly = start_rcelly + k * search_dir_y
       rcellz = start_rcellz + k * (start_rcellz/start_rcelly)*search_dir_y
    end if

    r_super_x = rcellx
    r_super_y = rcelly
    r_super_z = rcellz
    ! DRB added 2017/05/24 17:05
    ! We've changed the simulation cell. Now we must update grids and the density
    r_super_x_squared = r_super_x * r_super_x
    r_super_y_squared = r_super_y * r_super_y
    r_super_z_squared = r_super_z * r_super_z
    volume = r_super_x * r_super_y * r_super_z
    grid_point_volume = volume/(n_grid_x*n_grid_y*n_grid_z)
    one_over_grid_point_volume = one / grid_point_volume
    scale = (orcellx*orcelly*orcellz)/volume
    density = density * scale
    if(flag_diagonalisation) then
       do i = 1, nkp
          kk(1,i) = kk(1,i) * orcellx / rcellx
          kk(2,i) = kk(2,i) * orcelly / rcelly
          kk(3,i) = kk(3,i) * orcellz / rcellz
       end do
    end if
    do j = 1, maxngrid
       recip_vector(j,1) = recip_vector(j,1) * orcellx / rcellx
       recip_vector(j,2) = recip_vector(j,2) * orcelly / rcelly
       recip_vector(j,3) = recip_vector(j,3) * orcellz / rcellz
       xvec = recip_vector(j,1)/(two*pi)
       yvec = recip_vector(j,2)/(two*pi)
       zvec = recip_vector(j,3)/(two*pi)
       r2 = xvec*xvec + yvec*yvec + zvec*zvec
       if(j/=i0) hartree_factor(j) = one/r2 ! i0 notates gamma point
    end do
    do j = 1, ni_in_cell
       x_atom_cell(j) = (rcellx/orcellx)*x_atom_cell(j)
       y_atom_cell(j) = (rcelly/orcelly)*y_atom_cell(j)
       z_atom_cell(j) = (rcellz/orcellz)*z_atom_cell(j)
       if (inode == ionode .and. iprint_MD > 2) &
            write (io_lun,*) 'Position: ', j, x_atom_cell(j), &
            y_atom_cell(j), z_atom_cell(j)
    end do
    write(io_lun,*) "Iteration ", iter
    write(io_lun,*) "rcellx/start_rcellx = ", rcellx/start_rcellx
    write(io_lun,*) "rcelly/start_rcelly = ", rcelly/start_rcelly
    write(io_lun,*) "rcellz/start_rcellz = ", rcellz/start_rcellz
    !write(io_lun,*) 'new sim cell dims', start_rcellx, start_rcellx, start_rcellx
    write(io_lun,*) 'current sim cell dims', rcellx, rcelly, rcellz

  end subroutine update_cell_dims
  !!***

  !!****f* move_atoms/rescale_grids_and_density *
  !!  NAME
  !!   rescale_grids_and_density
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   scale the grid and density when the cell is updated
  !!   (adapted from move_atoms/update_cell_dims)
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2018/02/07
  !!  MODIFICATION HISTORY
  !! 
  !!  SOURCE
  !!
  subroutine rescale_grids_and_density(orcellx, orcelly, orcellz)

    use datatypes
    use numbers
    use global_module,  only: x_atom_cell, y_atom_cell, z_atom_cell, &
                              rcellx, rcelly, rcellz, flag_diagonalisation, &
                              iprint_MD
    use GenComms,       only: inode, ionode
    use maxima_module,  only: maxngrid
    use dimens,         only: r_super_x, r_super_y, r_super_z, &
                              r_super_x_squared, r_super_y_squared, &
                              r_super_z_squared, volume, grid_point_volume, &
                              one_over_grid_point_volume, n_grid_x, n_grid_y, &
                              n_grid_z
    use DiagModule,     only: kk, nkp
    use fft_module,     only: recip_vector, hartree_factor, i0
    use density_module, only: density
    use input_module,   only: leqi
    
    implicit none

    ! passed variables
    real(double), intent(in)  :: orcellx, orcelly, orcellz ! old cell
    ! local variables
    integer                   :: i, j
    real(double)              :: xvec, yvec, zvec, r2, scale

    if (inode==ionode .and. iprint_MD > 1) &
      write(io_lun,'(2x,a)') "move_atoms/rescale_grid_and_density"

    r_super_x = rcellx
    r_super_y = rcelly
    r_super_z = rcellz
    ! DRB added 2017/05/24 17:05
    ! We've changed the simulation cell. Now we must update grids and the density
    r_super_x_squared = r_super_x * r_super_x
    r_super_y_squared = r_super_y * r_super_y
    r_super_z_squared = r_super_z * r_super_z
    volume = r_super_x * r_super_y * r_super_z
    grid_point_volume = volume/(n_grid_x*n_grid_y*n_grid_z)
    one_over_grid_point_volume = one / grid_point_volume
    scale = (orcellx*orcelly*orcellz)/volume
    density = density * scale
    if(flag_diagonalisation) then
       do i = 1, nkp
          kk(1,i) = kk(1,i) * orcellx / rcellx
          kk(2,i) = kk(2,i) * orcelly / rcelly
          kk(3,i) = kk(3,i) * orcellz / rcellz
       end do
    end if
    do j = 1, maxngrid
       recip_vector(j,1) = recip_vector(j,1) * orcellx / rcellx
       recip_vector(j,2) = recip_vector(j,2) * orcelly / rcelly
       recip_vector(j,3) = recip_vector(j,3) * orcellz / rcellz
       xvec = recip_vector(j,1)/(two*pi)
       yvec = recip_vector(j,2)/(two*pi)
       zvec = recip_vector(j,3)/(two*pi)
       r2 = xvec*xvec + yvec*yvec + zvec*zvec
       if(j/=i0) hartree_factor(j) = one/r2 ! i0 notates gamma point
    end do

  end subroutine rescale_grids_and_density
  !!***

  ! -----------------------------------------------------------------------
  ! Subroutine update_pos_and_matrices
  ! -----------------------------------------------------------------------

  !!****f* move_atoms/update_pos_and_matrices *
  !!
  !!  NAME
  !!    update_pos_and_matrices
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Update information of atomic positions, neighbour lists, ...
  !!     AND...     matrices  
  !!
  !!  INPUTS
  !!    update_method :: showing which matrices will be updated.
  !!  USES
  !!
  !!  AUTHOR
  !!   Tsuyoshi Miyazaki
  !!  CREATION DATE
  !!   2017/Nov/12
  !!  MODIFICATION
  !!   2018/Sep/07  tsuyoshi
  !!       added calling ReportUpdateMatrix when flag_debug_move_atoms is true.
  !!   2018/Nov/13 17:30 nakata
  !!       changed matS to be spin_SF dependent
  !!   2019/Nov/14  tsuyoshi
  !!       removed glob2node_old, n_proc_old
  !!   2019/Jul/27  tsuyoshi
  !!       added atom_vels (from global_module), and removed local velocity_global
  !!
  !!  SOURCE
  !!
 subroutine update_pos_and_matrices(update_method, velocity)
  use datatypes
  use numbers,         only: half, zero, one, very_small
  use global_module,   only: flag_diagonalisation, atom_coord, atom_vels, atom_coord_diff, &
                             rcellx, rcelly, rcellz, ni_in_cell, nspin, nspin_SF, id_glob, &
                             area_moveatoms
    ! n_proc_old and glob2node_old have been removed
  use GenComms,        only: my_barrier, inode, ionode, cq_abort, gcopy
  use mult_module,     only: matL, L_trans, matK, matS, S_trans, matSFcoeff, SFcoeff_trans, &
                             matrix_scale, matrix_transpose, matSFcoeff_tran
  use matrix_data,     only: Lrange, Hrange, Srange, SFcoeff_range
  use store_matrix,    only: matrix_store_global, InfoMatrixFile, grab_InfoMatGlobal, grab_matrix2, &
                             set_atom_coord_diff
  use UpdateInfo, only: Matrix_CommRebuild, Report_UpdateMatrix
  use memory_module,   only: reg_alloc_mem, type_dbl, reg_dealloc_mem


  implicit none
  integer, intent(in) :: update_method
  real(double), intent(inout) :: velocity(3, ni_in_cell)
  logical :: flag_L, flag_K, flag_S, flag_SFcoeff, flag_X
  integer :: nfile, symm, ig, spin_SF
  logical :: fixed_potential 
   ! should be removed in the future (calling update_H is outside of this routine)
  real(double) :: scale_x, scale_y, scale_z, rms_change
  real(double) :: small_change = 0.3_double
  
 !H_trans is not prepared. If we need to symmetrise K, we need H_trans
  integer :: H_trans = 1
 
 !InfoGlob and Info can be defined locally.
 ! for extrapolation, we need to prepare multiple InfoGlob and Info.
  type(matrix_store_global) :: InfoGlob
  type(InfoMatrixFile),pointer :: InfoMat(:)

  integer :: i, stat

 !Switch on Debugging
 !  flag_debug_move_atoms = .true.

 !!! Note: for developers  !!!
 !  if you want to update some new matrix, you should
 !   1)define  updateXXX
 !   2)define  flagXXX
 !   3)check mult_module, and add "use mult_module, ONLY: matXXX, XXX_trans (if symm is needed).
 !   4)add if(flagXXX) statement
 !

 flag_L=.false.
 flag_K=.false.
 flag_S=.false.
 flag_X=.false.
 flag_SFcoeff=.false.

 select case(update_method)
  case(updatePos)
   if(inode .eq. ionode) write(io_lun,*) 'Update_Pos_and_Matrices:: only Positions are updated '
  case(updateL)
   flag_L = .true.
  case(updateK)
   flag_K = .true.
  case(updateLorK)
   if(flag_diagonalisation) then
     flag_K=.true.
   else
     flag_L=.true.
   endif
  case(updateS)
   flag_S=.true.
  case(updateX)
   if(inode .eq. ionode) write(io_lun,*) 'Update_Pos_and_Matrices:: updateX is not implemented yet.'
  case(updateSFcoeff)
   flag_SFcoeff=.true.
   if(flag_diagonalisation) then
     flag_K=.true.
   else
     flag_L=.true.
   endif
  case(extrplL)
   if(inode .eq. ionode) write(io_lun,*) 'Update_Pos_and_Matrices:: updateX is not implemented yet.'
   ! TM plans to implement  L(n)=2L(n-1)-L(n-2)
  case default
   if(inode .eq. ionode) write(io_lun,*) 'Update_Pos_and_Matrices:: Invalid Flag !',update_method
 end select ! case (update_method)


 !First updating information of atomic positions, neighbour lists, etc...
 !  2020/10/7 Tsuyoshi Miyazaki
 !   we are planning to use `atom_vels` and remove `velocity` (= direction in CG).
    if (.not. allocated(atom_vels)) then
       allocate(atom_vels(3,ni_in_cell), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating atom_vels in init_md: ", &
            ni_in_cell, stat)
       call reg_alloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)
       atom_vels = zero
    end if

     call wrap_xyz_atom_cell
     call update_atom_coord
     do i=1,ni_in_cell
        atom_vels(1,id_glob(i)) = velocity(1,i)
        atom_vels(2,id_glob(i)) = velocity(2,i)
        atom_vels(3,id_glob(i)) = velocity(3,i)
     end do
  !Before calling this routine, we need 1) call dump_InfoMatGlobal or 2) call set_InfoMatGlobal
  ! Then, we use InfoGlob read from the file or use InfoGlob as it is (in the case of 2))
  !       Now, we just assume 1).
     call grab_InfoMatGlobal(InfoGlob,index=0)
     call set_atom_coord_diff(InfoGlob)

  ! Since the order of the atoms in x,y,z_atom_cell and velocity (or direction in CG) changes 
  ! depending on their partitions, they are rearranged in updateIndices3.
  ! (these arrays should be replaced by atom_coord and atom_veloc in the future.)
  !    updateIndices3 : deallocates member of parts, bundles, covering sets, domain, ...
  !                     and allocates them following the new atomic positions.
  !   Now, old informaiton is stored in InfoGlob, thus some of the information is redundant 
  !    I (TM) should make new routine like "updateIndices4" using InfoGlob, soon.  ! 2017.Nov.13  Tsuyoshi Miyazaki
  !   (NOTE) If CONQUEST stops before calling "finalise", coord_next.dat should be used in the next job.
  !       coord_next.dat is made in "updateIndices3", at present.
     call updateIndices3(fixed_potential, velocity)
     do i=1,ni_in_cell
        velocity(1,i) = atom_vels(1,id_glob(i))
        velocity(2,i) = atom_vels(2,id_glob(i))
        velocity(3,i) = atom_vels(3,id_glob(i))
     end do

! 2020/Oct/12 TM   
!  if we want to reduce the memory size ...  
!      deallocate(atom_vels, STAT=stat)
!       if (stat /= 0) &
!            call cq_abort("Error deallocating atom_vels in init_md: stat=", stat)
!       call reg_dealloc_mem(area_moveatoms, 3*ni_in_cell, type_dbl)

 !
 ! Then, matrices will be read from the corresponding files
 !
  if(flag_L) then
     call grab_matrix2('L',inode,nfile,InfoMat,InfoGlob,index=0,n_matrix=nspin)
     call my_barrier()
     call Matrix_CommRebuild(InfoGlob,InfoMat,Lrange,L_trans,matL,nfile,symm,n_matrix=nspin)
      if(flag_debug_move_atoms) call Report_UpdateMatrix("Lmat")
  endif

  if(flag_K) then
     call grab_matrix2('K',inode,nfile,InfoMat,InfoGlob,index=0,n_matrix=nspin)
     call my_barrier()
     call Matrix_CommRebuild(InfoGlob,InfoMat,Hrange,H_trans,matK,nfile,n_matrix=nspin)
      if(flag_debug_move_atoms) call Report_UpdateMatrix("Kmat")
  endif

  if(flag_S) then
    ! If we introduce spin-dependent support, matS -> matS(nspin_SF)
     call grab_matrix2('S',inode,nfile,InfoMat,InfoGlob,index=0,n_matrix=nspin_SF)
     call my_barrier()
     call Matrix_CommRebuild(InfoGlob,InfoMat,Srange,S_trans,matS,nfile,symm,n_matrix=nspin_SF)
      if(flag_debug_move_atoms) call Report_UpdateMatrix("Smat")
  endif

  if(flag_SFcoeff) then
     do spin_SF = 1,nspin_SF
      call matrix_scale(zero,matSFcoeff(spin_SF))
     enddo !spin_SF = 1,nspin_SF

     call grab_matrix2('SFcoeff',inode,nfile,InfoMat,InfoGlob,index=0,n_matrix=nspin_SF)
     call my_barrier()
     call Matrix_CommRebuild(InfoGlob,InfoMat,SFcoeff_range,SFcoeff_trans,matSFcoeff,nfile,n_matrix=nspin_SF)
      if(flag_debug_move_atoms) call Report_UpdateMatrix("SFc1")

     do spin_SF = 1,nspin_SF
      call matrix_scale(zero,matSFcoeff_tran(spin_SF))
      call matrix_transpose(matSFcoeff(spin_SF), matSFcoeff_tran(spin_SF))
     enddo
  endif

 !Switch off Debugging
 !  flag_debug_move_atoms = .false.

 return
 end subroutine update_pos_and_matrices
 !!***

  !!****f* control/propagate_vector *
  !!
  !!  NAME
  !!   propagate_vector
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Generate new configuration by shifting the position vector config in 
  !!   direction force by distance k
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Z Raza
  !!  CREATION DATE
  !!   2019/02/08
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  subroutine propagate_vector(force, config, config_new, cell_ref, k)

    use GenComms,      only: inode, ionode
    use global_module, only: iprint_MD, ni_in_cell, id_glob

    implicit none

    ! passed variables
    real(double), dimension(:,:), intent(in)  :: force, config
    real(double), dimension(:,:), intent(out) :: config_new
    real(double), dimension(3), intent(in)    :: cell_ref
    real(double), intent(in)                  :: k

    ! local variables
    integer       :: i, j, i_glob
    real(double)  :: d

    if (inode==ionode .and. iprint_MD > 2) &
      write(io_lun,'(2x,a,f12.6)') "move_atoms/propagate_vector: k=",k

    ! config_new = config + k*force
    do i=1,ni_in_cell+1
      do j=1,3
        config_new(j,i) = config(j,i) + k*force(j,i)
      end do
    end do

  end subroutine propagate_vector
  !!***

  !!****f* control/vector_to_cq *
  !!
  !!  NAME
  !!   vector_to_cq
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Convert optimisable vector config (containing strains and fractional
  !!   coordinates) to Conquest variables (Cartesian coordinates) for
  !!   electronic structure calculation. Increment the positions using the
  !!   force vector, containing scaled stresses and ionic forces in lattice
  !!   coordinates, as defined in !!   Pfrommer et al. J. Comput. Phys. 131,
  !!    233 (1997)
  !!
  !!   Stress components: 
  !!      f_sigma,i = -(sigma_i + pV)(1 + epsilon_i)^-1
  !!   f_sigma,i = ith component of force on cell as defined above
  !!   sigma_i = ith component of stress
  !!   epsilon = ith component of strain
  !!   ionic force components:
  !!      F_i = h^T h f_i
  !!   h = matrix of lattice vectors (h^T h = g, metric tensor)
  !!   f_i = force on atom i, as computed from gradient of energy
  !! 
  !!   force contains scaled stresses force(:,1) and ionic forces in lattice
  !!   coordinates force(:,2:)
  !!
  !!   config contains strains config(:,1) and fractional coordinates
  !!   config(:,2:)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Z Raza
  !!  CREATION DATE
  !!   2019/02/06
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  subroutine vector_to_cq(config, cell_ref, orcellx, orcelly, orcellz)

    use numbers
    use GenComms,      only: inode, ionode
    use global_module, only: rcellx, rcelly, rcellz, ni_in_cell, &
                             iprint_MD, id_glob_inv, atom_coord

    implicit none

    ! passed variables
    real(double), dimension(:,:), intent(in)  :: config
    real(double), dimension(3), intent(in)    :: cell_ref
    real(double), intent(out) :: orcellx, orcelly, orcellz

    ! local variables
    integer       :: i, i_global
    real(double)  :: dx, dy, dz

    if (inode==ionode .and. iprint_MD > 2) &
      write(io_lun,'(2x,a)') "move_atoms/vector_to_cq"


    orcellx = rcellx
    orcelly = rcelly
    orcellz = rcellz

    rcellx = (one + config(1,ni_in_cell+1))*cell_ref(1)
    rcelly = (one + config(2,ni_in_cell+1))*cell_ref(2)
    rcellz = (one + config(3,ni_in_cell+1))*cell_ref(3)
    do i=1,ni_in_cell
      atom_coord(1,i) = config(1,i)*rcellx
      atom_coord(2,i) = config(2,i)*rcelly
      atom_coord(3,i) = config(3,i)*rcellz
    end do

    ! Now we've changed atom_coord, we want these changes to be reflected
    ! in x/y/z_atom_cell
    call update_r_atom_cell
    call rescale_grids_and_density(orcellx, orcelly, orcellz)

  end subroutine vector_to_cq
  !!***

  !!****f* control/cq_to_vector *
  !!
  !!  NAME
  !!   cq_to_vector
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Convert Conquest variables x_atom_cell, y_atom_cell, z_atom_cell,
  !!   rcellx, rcelly, rcellz and tot_force to the vectors required for full
  !!   cell optimisation, as defined in Pfrommer et al. 
  !!   J. Comput. Phys. 131, 233 (1997) (see cq_to_vector)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Z Raza
  !!  CREATION DATE
  !!   2019/02/06
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine cq_to_vector(force, config, cell_ref, target_press)
 
    use numbers
    use global_module, only: rcellx, rcelly, rcellz, ni_in_cell, &
                             iprint_MD, id_glob, atom_coord
    use GenComms,      only: inode, ionode
    use force_module,  only: stress, tot_force

    implicit none
 
    ! passed variables
    real(double), dimension(:,:), intent(out) :: force, config
    real(double), dimension(3), intent(in)    :: cell_ref
    real(double), intent(in)                  :: target_press

    ! local variables
    integer                     :: i, i_global
    real(double)                :: vol
    real(double), dimension(3)  :: one_plus_strain

    if (inode==ionode .and. iprint_MD>2) &
      write(io_lun,'(2x,a)') "move_atoms/cq_to_vector"
 
    vol = rcellx*rcelly*rcellz
    one_plus_strain(1) = rcellx/cell_ref(1)
    one_plus_strain(2) = rcelly/cell_ref(2)
    one_plus_strain(3) = rcellz/cell_ref(3)
    do i=1,3
      config(i,ni_in_cell+1) = one_plus_strain(i) - one
      force(i,ni_in_cell+1) = &
        -(stress(i,i) + target_press*vol)/one_plus_strain(i)
    end do
    do i=1,ni_in_cell
      config(1,i) = atom_coord(1,i)/rcellx
      config(2,i) = atom_coord(2,i)/rcelly
      config(3,i) = atom_coord(3,i)/rcellz
      force(1,i) = tot_force(1,i)*rcellx
      force(2,i) = tot_force(2,i)*rcelly
      force(3,i) = tot_force(3,i)*rcellz
    end do

    if (inode==ionode .and. iprint_MD>1) then
      write(io_lun,'(4x,a)') "Position vector:"
      do i=1,ni_in_cell
        write(io_lun,'(4x,i8,3f20.12)') i, config(:,i)
      end do
      write(io_lun,'(4x,a8,3f20.12)') "strain: ", config(:,ni_in_cell+1)
      write(io_lun,'(4x,a)') "Force vector:"
      do i=1,ni_in_cell
        write(io_lun,'(4x,i8,3f20.12)') i, force(:,i)
      end do
      write(io_lun,'(4x,a8,3f20.12)') "stress: ", force(:,ni_in_cell+1)
      write(io_lun,*)
    end if

  end subroutine cq_to_vector
  !!***

  !!****f* move_atoms/enthalpy *
  !!
  !!  NAME 
  !!   enthalpy
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Compute the enthalpy for cell optimisation
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   Zamaan Raza
  !!  CREATION DATE
  !!   2019/02/07
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  function enthalpy(e, p) result(h)

    use datatypes
    use GenComms,       only: inode, ionode
    use global_module,  only: rcellx, rcelly, rcellz, iprint_MD, flag_MDdebug

    implicit none

    ! passed variables
    real(double), intent(in)  :: e, p    
    real(double)              :: h

    ! local variables
    real(double)              :: pv

    pv = p*(rcellx*rcelly*rcellz)
    h = e + pv

    if (inode==ionode .and. iprint_MD > 2) then
      write(io_lun,'(2x,a)') "move_atoms/enthalpy"
      if (flag_MDdebug) then
        write(io_lun,'(4x,"energy   = ",f16.8)') e
        write(io_lun,'(4x,"P        = ",f16.8)') p
        write(io_lun,'(4x,"V        = ",f16.8)') rcellx*rcelly*rcellz
        write(io_lun,'(4x,"PV       = ",f16.8)') pv
        write(io_lun,'(4x,"enthalpy = ",f16.8)') h
      end if 
    end if

  end function enthalpy
  !!***

end module move_atoms
