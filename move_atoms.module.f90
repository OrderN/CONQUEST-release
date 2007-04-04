! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!    Created update_H to reproject blips, build S, reproject pseudos, build n(r) and H
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
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!  SOURCE
!!
module move_atoms

  use datatypes

  ! Useful physical constants
  real(double), parameter:: amu = 1.660566e-27_double
!  real(double), parameter:: ang = 1.0e-10_double
  real(double), parameter:: ang = 0.529177e-10_double
  real(double), parameter:: tscale = 1.0e-15_double
!  real(double), parameter:: ev = 1.602189e-19_double
  real(double), parameter:: ev = 2.0_double*13.6058_double*1.602189e-19_double
  real(double), parameter:: fac = amu*ang*ang/(tscale*tscale*ev)

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), save, private :: RCSid = "$Id$"
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

    use set_blipgrid_module, ONLY: free_blipgrid
    use set_bucket_module, ONLY: free_bucket
    use functions_on_grid, ONLY: dissociate_fn_on_grid

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
!!  TODO
!!   Proper buffer zones for matrix mults so initialisation doesn't have
!!   to be done at every step 03/07/2001 dave
!!  SOURCE
!!
  subroutine velocityVerlet(prim, step,T,KE,quenchflag,velocity,force)

    use datatypes
    use basic_types
    use global_module, ONLY: iprint_MD, x_atom_cell, y_atom_cell, z_atom_cell, ni_in_cell
    use species_module, ONLY: species, mass
    use numbers
    use GenComms, ONLY: myid

    implicit none

    ! Passed variables
    real(double) :: step, T, KE
    type(primary_set) :: prim
    real(double), dimension(3,ni_in_cell) :: velocity
    real(double), dimension(3,ni_in_cell) :: force
    logical :: quenchflag

    ! Local variables
    integer :: part, memb, atom, speca, k
    real(double) :: massa, acc

    if(myid==0.AND.iprint_MD>0) write(*,1) step,quenchflag
1   format(4x,'In velocityVerlet, timestep is ',f10.5/&
         'Quench is ',l3)
    do atom = 1, ni_in_cell
       speca = species(atom) 
       massa = mass(speca)*fac
       if(quenchflag) then
          do k=1,3
             if(velocity(k,atom)*force(k,atom)<zero) &
                  velocity(k,atom) = zero
             velocity(k,atom) = velocity(k,atom)+ &
                  step*force(k,atom)/(two*massa)
          end do
       else
          do k=1,3
             velocity(k,atom) = velocity(k,atom)+ &
                  step*force(k,atom)/(two*massa)
          enddo
       endif
    end do
    ! Maybe fiddle with KE
    KE = 0.0_double
    do atom = 1, ni_in_cell
       speca = species(atom) 
       massa = mass(speca)*fac
       KE = KE + half*massa*velocity(k,atom)*velocity(k,atom)
    end do
    ! Update positions and velocities
    do atom = 1, ni_in_cell
       speca = species(atom) 
       massa = mass(speca)*fac
       ! X
       acc = force(1,atom)/massa
       x_atom_cell(atom) = x_atom_cell(atom) + &
            step*velocity(1,atom) + half*step*step*acc
       velocity(1,atom) = velocity(1,atom)+half*step*acc
       ! Y
       acc = force(2,atom)/massa
       y_atom_cell(atom) = y_atom_cell(atom) + &
            step*velocity(2,atom) + half*step*step*acc
       velocity(2,atom) = velocity(2,atom)+half*step*acc
       ! Z
       acc = force(3,atom)/massa
       z_atom_cell(atom) = z_atom_cell(atom) + &
            step*velocity(3,atom) + half*step*step*acc
       velocity(3,atom) = velocity(3,atom)+half*step*acc
    end do
!Update atom_coord : TM 27Aug2003
    call update_atom_coord
!Update atom_coord : TM 27Aug2003
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
!!    Added calls to get_E_and_F and added passed variables for get_E_and_F
!!   08:57, 2003/02/05 dave
!!    Sorted out arguments to pass to updateIndices
!!   14:41, 26/02/2003 drb 
!!    Added n_atoms from atoms
!!   08:50, 11/05/2005 dave 
!!    Added code to write out atomic positions during minimisation; also added code to subtract off atomic
!!    densities for old atomic positions and add it back on for new ones after atoms moved; this is commented
!!    out as it's not been tested or checked rigorously
!!   09:51, 25/10/2005 drb 
!!    Added correction so that the present energy is passed to get_E_and_F for blip minimisation loop
!!  SOURCE
!!
  subroutine safemin(start_x,start_y,start_z,direction,energy_in,energy_out,&
       fixed_potential, vary_mu, number_of_bands, mu, total_energy)

    ! Module usage
    use datatypes
    use numbers
    use units
    use global_module, ONLY: iprint_MD, x_atom_cell, y_atom_cell, z_atom_cell, &
         flag_vary_basis, atom_coord, ni_in_cell, rcellx, rcelly, rcellz, flag_self_consistent
    use minimise, ONLY: get_E_and_F, sc_tolerance, L_tolerance, n_L_iterations
    use GenComms, ONLY: my_barrier, myid, inode, ionode
    use SelfCon, ONLY: new_SC_potl
    use GenBlas, ONLY: dot
    use force_module, ONLY: tot_force
    use io_module, ONLY: write_atomic_positions, pdb_template

    implicit none

    ! Passed variables
    real(double), dimension(3,ni_in_cell) :: direction
    real(double), dimension(ni_in_cell) :: start_x, start_y, start_z
    real(double) :: energy_in, energy_out
    ! Shared variables needed by get_E_and_F for now (!)
    logical :: vary_mu, fixed_potential

    character(len=40) :: output_file

    real(double) :: number_of_bands, mu
    real(double) :: total_energy

    ! Local variables
    real(double) :: k0, k1, k2, k3, lambda, k3old
    real(double), save :: kmin = 0.0_double, dE = 0.0_double
    real(double) :: e0, e1, e2, e3, tmp
    integer :: i,j, iter, lun
    logical :: reset_L = .false.

    e0 = total_energy
    if(inode==ionode.AND.iprint_MD>0) &
         write(*,fmt='(4x,"In safemin, initial energy is ",f15.10," ",a2)') en_conv*energy_in,en_units(energy_units)
    if(inode==ionode) write(*,fmt='(/4x,"Seeking bracketing triplet of points"/)')
    ! Unnecessary and over cautious !
    k0 = 0.000001_double
    !do i=1,ni_in_cell
    !   x_atom_cell(i) = start_x(i) + k0*direction(1,i)
    !   y_atom_cell(i) = start_y(i) + k0*direction(2,i)
    !   z_atom_cell(i) = start_z(i) + k0*direction(3,i)
    !   !write(*,*) 'Position: ',i,x_atom_cell(i),y_atom_cell(i),z_atom_cell(i)
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
    k1 = 0.0_double
    e1 = energy_in
    k2 = k0
    e2 = e0
    e3 = e2
    k3 = 0.0_double
    k3old = k3
    lambda = 2.0_double
    ! Loop to find a bracketing triplet
    do while(e3<=e2)
       if (k2==k0) then
          !k3 = 0.001_double
          !if(abs(kmin) < very_small) then
          if(abs(dE) < very_small) then
             if(k3<very_small) then ! First guess
                k3 = 0.70_double!k3old/lambda
             end if
          else
             tmp = dot(3*ni_in_cell,direction,1,tot_force,1)
             k3 = 0.5_double*dE/tmp!kmin/lambda
             if(abs(k3)>abs(kmin)) k3 = kmin
          endif
          !   k3 = 0.1_double
          !else
          !   k3 = kmin/lambda
          !endif
       elseif (k2==0.01_double) then
          k3 = 0.01_double
       else
          k3 = lambda*k2          
       endif
!       k3 = 0.032_double
       ! These lines calculate the difference between atomic densities and total density
       !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       !   ! Subtract off atomic densities
       !   store_density = density
       !   call set_density()
       !   density = store_density - density
       !end if
       ! Move atoms
       do i=1,ni_in_cell
          x_atom_cell(i) = start_x(i) + k3*direction(1,i)
          y_atom_cell(i) = start_y(i) + k3*direction(2,i)
          z_atom_cell(i) = start_z(i) + k3*direction(3,i)
          !write(*,*) 'Position: ',i,x_atom_cell(i),y_atom_cell(i),z_atom_cell(i)
       end do
!Update atom_coord : TM 27Aug2003
       call update_atom_coord
!Update atom_coord : TM 27Aug2003
       ! Update indices and find energy and forces
       !call updateIndices(.false.,fixed_potential, number_of_bands)
       call updateIndices(.true.,fixed_potential, number_of_bands)
       ! These lines add back on the atomic densities for NEW atomic positions
       !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       !   ! Add on atomic densities
       !   store_density = density
       !   call set_density()
       !   density = store_density + density
       !end if
       ! We've just moved the atoms - we need a self-consistent ground state before we can minimise blips !
       ! Write out atomic positions
       if(iprint_MD>2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat",trim(pdb_template))
       end if
       if(flag_vary_basis) then
          call new_SC_potl( .false., sc_tolerance, reset_L, fixed_potential, vary_mu, n_L_iterations, &
               number_of_bands, L_tolerance, mu, e3)
       end if
       call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, e3, .false., .false.)
       if(inode==ionode.AND.iprint_MD>1) &
            write(*,fmt='(4x,"In safemin, iter ",i3," step and energy are ",2f15.10)') iter,k3,e3
       if (e3<e2) then ! We're still going down hill
          k1 = k2
          e1 = e2
          k2 = k3
          e2 = e3
          iter=iter+1
       else if(k2==0.0_double) then
          k3old = k3
          if(abs(dE)<very_small) then 
             k3 = k3old/2.0_double
             dE = 1.0_double
          else
             dE = 0.5_double*dE
          end if
          e3 = e2
       endif
    end do
    if(inode==ionode) write(*,fmt='(/4x,"Interpolating minimum"/)')
    ! Interpolate to find minimum.
    kmin = 0.5_double*(((k1*k1-k3*k3)*(e1-e2)-(k1*k1-k2*k2)*(e1-e3))/((k1-k3)*(e1-e2)-(k1-k2)*(e1-e3)))
    !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
    !   ! Subtract off atomic densities
    !   store_density = density
    !   call set_density()
    !   density = store_density - density
    !end if
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
    call updateIndices(.true.,fixed_potential, number_of_bands)
    !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
    !   ! Add on atomic densities
    !   store_density = density
    !   call set_density()
    !   density = store_density + density
    !end if
    if(iprint_MD>2) then
       call write_atomic_positions("UpdatedAtoms_tmp.dat",trim(pdb_template))
    end if
    ! We've just moved the atoms - we need a self-consistent ground state before we can minimise blips !
    if(flag_vary_basis) then
       call new_SC_potl( .false., sc_tolerance, reset_L, fixed_potential, vary_mu, n_L_iterations, &
            number_of_bands, L_tolerance, mu, e3)
    end if
    energy_out = e3
    if(iprint_MD>0) then
       call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy_out, .true., .true.)
    else
       call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy_out, .true., .false.)
    end if
    if (energy_out>e2) then ! The interpolation failed - go back
       if(inode==ionode) write(*,fmt='(/4x,"Interpolation failed; reverting"/)')
       kmin = k2
       !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       !   ! Subtract off atomic densities
       !   store_density = density
       !   call set_density()
       !   density = store_density - density
       !end if
       do i=1,ni_in_cell
          x_atom_cell(i) = start_x(i) + kmin*direction(1,i)
          y_atom_cell(i) = start_y(i) + kmin*direction(2,i)
          z_atom_cell(i) = start_z(i) + kmin*direction(3,i)
       end do
!Update atom_coord : TM 27Aug2003
       call update_atom_coord
!Update atom_coord : TM 27Aug2003
       call updateIndices(.true.,fixed_potential, number_of_bands)
       !call updateIndices(.false.,fixed_potential, number_of_bands)
       !if(flag_self_consistent.AND.(.NOT.flag_no_atomic_densities)) then
       !   ! Add on atomic densities
       !   store_density = density
       !   call set_density()
       !   density = store_density + density
       !end if
       if(iprint_MD>2) then
          call write_atomic_positions("UpdatedAtoms_tmp.dat",trim(pdb_template))
       end if
       ! We've just moved the atoms - we need a self-consistent ground state before we can minimise blips !
       if(flag_vary_basis) then
          call new_SC_potl( .false., sc_tolerance, reset_L, fixed_potential, vary_mu, n_L_iterations, &
               number_of_bands, L_tolerance, mu, e3)
       end if
       energy_out = e3
       if(iprint_MD>0) then
          call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy_out, .true., .true.)
       else
          call get_E_and_F(fixed_potential, vary_mu, number_of_bands, mu, energy_out, .true., .false.)
       end if
    end if
    dE = e0 - energy_out
7   format(4x,3f15.8)
    if(inode==ionode.AND.iprint_MD>0) then
       write(*,fmt='(4x,"In safemin, exit after ",i4," iterations with energy ",f20.10)') iter,energy_out
    else
       write(*,fmt='(/4x,"Final energy: ",f15.10)/') energy_out
    end if
    return
  end subroutine safemin
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
!!   Updates the indices for matrices, saves relevant information and stores (if necessary) old L matrix
!!
!!   At the simplest, all this does is update the positions of the atoms in the primary and covering sets, 
!!   and rebuild the Hamiltonian
!!  INPUTS
!!   logical :: matrix_update Flags whether the user wants ALL matrix information updated
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
!!    D'oh ! Put in a call to cover_update for ewald_CS so that the new ewald routines work
!!   08:49, 11/05/2005 dave 
!!    Added lines which check for change in band energy, and reset DM if too large; these are commented out
!!    as these ideas are not rigorously tested
!!   2006/09/08 07:59 dave
!!    Various changes for dynamic allocation
!!  TODO
!!   Think about updating radius component of matrix derived type, or eliminating it !
!!  SOURCE
!!
  subroutine updateIndices(matrix_update,fixed_potential, number_of_bands)

    ! Module usage
    use datatypes
    use mult_module, ONLY: fmmi, immi
    use matrix_module, ONLY: allocate_matrix, deallocate_matrix, set_matrix_pointers2, matrix
    use group_module, ONLY: parts
    use cover_module, ONLY : BCS_parts, DCS_parts, ewald_CS
    use primary_module, ONLY : bundle
    use global_module, ONLY: iprint_MD, x_atom_cell, y_atom_cell, z_atom_cell
    use matrix_data, ONLY: Hrange, mat, rcut
    use maxima_module, ONLY: maxpartsproc
    use set_blipgrid_module, ONLY: set_blipgrid
    use set_bucket_module,   ONLY: set_bucket
    use dimens, ONLY: r_core_squared,r_h, RadiusSupport
    use pseudopotential_common, ONLY: core_radius
    use GenComms, ONLY: myid, cq_abort, gsum
    use functions_on_grid, ONLY: associate_fn_on_grid
    use ewald_module, ONLY: flag_old_ewald
    use blip, ONLY: Extent
    use numbers

    implicit none

    ! Passed variables
    logical, intent(in) :: matrix_update

    ! Shared variables needed by get_H_matrix for now (!)
    logical :: fixed_potential

    real(double) :: number_of_bands

    ! Local variables
    logical :: check
    integer :: i,k,stat

    ! Update positions in primary and covering sets
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle, parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, DCS_parts, parts)
    if(.NOT.flag_old_ewald) call cover_update(x_atom_cell, y_atom_cell, &
         z_atom_cell, ewald_CS, parts)
    ! If there's a new interaction of Hamiltonian range, then we REALLY need to rebuild the matrices etc
    call checkBonds(check,bundle,BCS_parts,mat(1,Hrange),maxpartsproc,rcut(Hrange))
    ! If one processor gets a new bond, they ALL need to redo the indices
    call gsum(check)
    ! There's also an option for the user to force it via matrix_update (which could be set to every n iterations ?)
    if(check.OR.matrix_update) then
       ! Deallocate all matrix storage
       ! finish blip-grid indexing
       call finish_blipgrid
       ! finish matrix multiplication indexing
       call fmmi(bundle)
       ! Reallocate and find new indices
       call immi(parts,bundle,BCS_parts,myid+1)
       ! Reallocate for blip grid
       call set_blipgrid(myid,RadiusSupport,core_radius,Extent)
       !call set_blipgrid(myid,r_h,sqrt(r_core_squared))
       call set_bucket(myid)
       call associate_fn_on_grid
    end if
    ! Rebuild S, n(r) and hamiltonian based on new positions
    call update_H(fixed_potential, number_of_bands)
    return
  end subroutine updateIndices
!!***

  subroutine updateIndices2(matrix_update,fixed_potential, number_of_bands)

    ! Module usage
    use datatypes
    use mult_module, ONLY: fmmi, immi
    use matrix_module, ONLY: allocate_matrix, deallocate_matrix, set_matrix_pointers2, matrix
    use group_module, ONLY: parts
    use cover_module, ONLY : BCS_parts, DCS_parts, ewald_CS
    use primary_module, ONLY : bundle
    use global_module, ONLY: iprint_MD, x_atom_cell, y_atom_cell, z_atom_cell
    use matrix_data, ONLY: Hrange, mat, rcut
    use maxima_module, ONLY: maxpartsproc
    use set_blipgrid_module, ONLY: set_blipgrid
    use set_bucket_module,   ONLY: set_bucket
    use dimens, ONLY: r_core_squared,r_h, RadiusSupport
    use pseudopotential_common, ONLY: core_radius
    use GenComms, ONLY: myid, cq_abort, gsum
    use functions_on_grid, ONLY: associate_fn_on_grid
    use ewald_module, ONLY: flag_old_ewald
    use blip, ONLY: Extent
    use numbers

    implicit none

    ! Passed variables
    logical, intent(in) :: matrix_update

    ! Shared variables needed by get_H_matrix for now (!)
    logical :: fixed_potential

    real(double) :: number_of_bands

    ! Local variables
    logical :: check
    integer :: i,k,stat

    ! Update positions in primary and covering sets
    call primary_update(x_atom_cell, y_atom_cell, z_atom_cell, bundle, parts, myid)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, BCS_parts, parts)
    call cover_update(x_atom_cell, y_atom_cell, z_atom_cell, DCS_parts, parts)
    if(.NOT.flag_old_ewald) call cover_update(x_atom_cell, y_atom_cell, &
         z_atom_cell, ewald_CS, parts)
    check = .false.
    ! If there's a new interaction of Hamiltonian range, then we REALLY need to rebuild the matrices etc
    call checkBonds(check,bundle,BCS_parts,mat(1,Hrange),maxpartsproc,rcut(Hrange))
    ! If one processor gets a new bond, they ALL need to redo the indices
    call gsum(check)
    ! There's also an option for the user to force it via matrix_update (which could be set to every n iterations ?)
    if(check.OR.matrix_update) then
       ! Deallocate all matrix storage
       ! finish blip-grid indexing
       call finish_blipgrid
       ! finish matrix multiplication indexing
       call fmmi(bundle)
       ! Reallocate and find new indices
       call immi(parts,bundle,BCS_parts,myid+1,1)
       ! Reallocate for blip grid
       call set_blipgrid(myid,RadiusSupport,core_radius,Extent)
       !call set_blipgrid(myid,r_h,sqrt(r_core_squared))
       call set_bucket(myid)
       call associate_fn_on_grid
    end if
    ! Rebuild S, n(r) and hamiltonian based on new positions
    call update_H(fixed_potential, number_of_bands)
    return
  end subroutine updateIndices2


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
!!  SOURCE
!!
  subroutine update_H(fixed_potential, number_of_bands)

    use S_matrix_module, ONLY: get_S_matrix
    use H_matrix_module, ONLY: get_H_matrix
    use DiagModule, ONLY: diagon
    use mult_module, ONLY: LNV_matrix_multiply
    use ewald_module, ONLY: ewald, mikes_ewald, flag_old_ewald
    use pseudopotential_data, ONLY: init_pseudo
    use pseudo_tm_module, ONLY: set_tm_pseudo
    use pseudopotential_common, ONLY: pseudo_type, OLDPS, SIESTA, STATE, ABINIT, core_correction
    use logicals
    use global_module, ONLY: iprint_MD, flag_self_consistent
    use density_module, ONLY: set_density, flag_no_atomic_densities, density
    use GenComms, ONLY: cq_abort, inode, ionode
    use maxima_module, ONLY: maxngrid
    
    implicit none

    ! Shared variables needed by get_H_matrix for now (!)
    logical :: fixed_potential

    real(double) :: number_of_bands

    ! Local variables
    real(double) :: tmp

    ! (1) Get S matrix (includes blip-to-grid transform)
    call get_S_matrix(inode, ionode)

    ! (2) get K matrix if O(N)
    if(.NOT.diagon) then 
       call LNV_matrix_multiply(tmp, tmp, doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
    end if
    ! (3) Find the Ewald energy for the initial set of atoms
    if(flag_old_ewald) then
       call ewald
    else
       call mikes_ewald
    end if
    ! (4) Pseudopotentials: choose correct form
    select case(pseudo_type) 
    case(OLDPS)
       call init_pseudo(number_of_bands, core_correction)
    case(SIESTA)
       call set_tm_pseudo
    case(ABINIT)
       call set_tm_pseudo
    end select
    ! Now we call set_density if we're using atomic densities
    if((.NOT.flag_self_consistent).AND.(.NOT.flag_no_atomic_densities)) then
       call set_density
    else if((.NOT.flag_self_consistent).AND.(flag_no_atomic_densities)) then
       call cq_abort("update_H: Can't run non-self-consistent without PAOs !")
    end if
    ! (5) Now generate a new H matrix, including a new charge density
    call get_H_matrix(.true., fixed_potential, tmp, density, maxngrid)
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
!!   Relies heavily on the methodology of get_naba: loops over GCS partitions and atoms, checking for 
!!   atoms within range and compares to known atoms.  Note various things:
!!
!!    i) We can rely on the order of the GCS not changing
!!    ii) We only want to compare the partition and sequence numbers of the atoms, not their separation
!!    iii) As an extra check, we compare the number of neighbours of primary set atoms
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
!!
!!  SOURCE
!!
  subroutine checkBonds(newAtom,prim,gcs,amat,partsproc,rcut)

    use datatypes
    use basic_types, ONLY: primary_set, cover_set
    use matrix_module, ONLY: matrix

    implicit none

    ! Passed variables
    logical, intent(out) :: newAtom
    integer :: partsproc
    type(primary_set) :: prim
    type(cover_set) :: gcs
    real(double) :: rcut
    type(matrix) :: amat(partsproc)

    ! Local variables
    real(double) :: rcutsq, dx, dy, dz
    real(double) :: tol = 1.0e-8_double
    integer :: inp, nn, j, np, ni, ist, n_nab

    rcutsq=rcut*rcut
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
                            if(np/=amat(nn)%i_part(ist).AND.ni/=amat(nn)%i_seq(ist)) newAtom = .true.
                         end if
                      endif
                   enddo ! End n_inp_cover
                endif
             enddo ! End np_cover
             inp=inp+1  ! Indexes primary-set atoms
             if(n_nab/=amat(nn)%n_nab(j)) newAtom = .true.
          enddo ! End prim%nm_nodgroup
       endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node
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
!!
!!  SOURCE
!!
  subroutine primary_update(x_position,y_position,z_position,prim,groups,myid)

    use datatypes
    use basic_types
    use global_module, ONLY: ni_in_cell, rcellx, rcelly, rcellz

    implicit none

    ! Passed variables
    real(double), dimension(ni_in_cell) :: x_position,y_position,z_position
    type(primary_set) :: prim
    type(group_set) :: groups
    integer :: myid

    ! Local variables
    integer :: ng, ind_group, nx, ny, nz, nx1, ny1, nz1, nnd, n_prim, ni
    real(double) :: dcellx, dcelly, dcellz
    real(double) :: xadd, yadd, zadd

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
!!
!!  SOURCE
!!
  subroutine cover_update(x_position,y_position,z_position,set,groups)

    use datatypes
    use basic_types
    use global_module, ONLY: ni_in_cell, rcellx, rcelly, rcellz
    use cover_module, ONLY: indexx
    use GenComms, ONLY: cq_abort, myid

    implicit none

    ! Passed variables
    real(double), dimension(ni_in_cell) :: x_position,y_position,z_position
    type(cover_set) :: set
    type(group_set) :: groups

    ! Local variables
    integer :: ind_cover, nx, ny, nz, nx_o, ny_o, nz_o, ni
    integer :: nsx,nsy,nsz,nqx,nqy,nqz,nmodx,nmody,nmodz
    integer :: cover_part, ind_qart
    real(double) :: dcellx, dcelly, dcellz
    real(double) :: xadd, yadd, zadd

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

    nnd = myid+1
    allocate(nx_in_cover(set%ng_cover),ny_in_cover(set%ng_cover),nz_in_cover(set%ng_cover),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating nx_in_cover: ",set%ng_cover,stat)
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
    if(minx>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in x-edge',minx)
    endif
    do ngcx=1,minx
       if(ngcx<=nremx) then
          nrepx(ngcx)=noccx+1
       else
          nrepx(ngcx)=noccx
       endif
    enddo
    ! ... y-direction
    noccy=set%ncovery/groups%ngcelly
    nremy=set%ncovery-noccy*groups%ngcelly
    miny=min(set%ncovery,groups%ngcelly)
    if(miny>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in y-edge',miny)
    endif
    do ngcy=1,miny
       if(ngcy<=nremy) then
          nrepy(ngcy)=noccy+1
       else
          nrepy(ngcy)=noccy
       endif
    enddo
    ! ... z-direction
    noccz=set%ncoverz/groups%ngcellz
    nremz=set%ncoverz-noccz*groups%ngcellz
    minz=min(set%ncoverz,groups%ngcellz)
    if(minz>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in z-edge',minz)
    endif
    do ngcz=1,minz
       if(ngcz<=nremz) then
          nrepz(ngcz)=noccz+1
       else
          nrepz(ngcz)=noccz
       endif
    enddo
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
             ind_qart=(nqx-1)*groups%ngcelly*groups%ngcellz+(nqy-1)*groups%ngcellz+nqz
             ino=groups%inv_ngnode(ind_qart)
             ind_min(ind)=1+mod(ino-groups%inode_beg(nnd)+ ng_in_cell,ng_in_cell)
             ngcx_min(ind)=ngcx
             ngcy_min(ind)=ngcy
             ngcz_min(ind)=ngcz
          enddo
       enddo
    enddo
    ! sort minimum CS by nodes 
    call indexx(groups%mx_gcell,ng_in_min,ind_min,min_sort)
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
    deallocate(nx_in_cover,ny_in_cover,nz_in_cover,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating nx_in_cover: ",set%ng_cover,stat)
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
!!
!!  SOURCE
!!
  subroutine update_atom_coord
    
    use global_module, ONLY: x_atom_cell, y_atom_cell, z_atom_cell, &
         id_glob, atom_coord, ni_in_cell

    implicit none

    integer :: ni, id_global

    do ni=1, ni_in_cell
       id_global=id_glob(ni)
       atom_coord(1, id_global)= x_atom_cell(ni)
       atom_coord(2, id_global)= y_atom_cell(ni)
       atom_coord(3, id_global)= z_atom_cell(ni)
    enddo
    return
  end subroutine update_atom_coord
!!***
end module move_atoms
