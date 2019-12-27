! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module initialisation
! ------------------------------------------------------------------------------
! Code area 1: initialisation
! ------------------------------------------------------------------------------

!!****h* Conquest/initialisation *
!!  NAME
!!   initialisation
!!  PURPOSE
!!   Hold initialisation routines
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006/10/16 (bringing together existing routines)
!!  MODIFICATION HISTORY
!!   2008/02/06 08:06 dave
!!    Changed for output to file not stdout
!!   2008/05/15 ast
!!    Added some timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module initialisation

  use datatypes
  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer, cq_timer 
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_initialisation,  &
                                    tmr_std_densitymat,      &
                                    tmr_std_matrices
  use timer_module,           only: init_timing_system

  implicit none

  ! Area identification
  integer, parameter, private :: area = 1

  ! RCS tag for object file identification
  character(len=80), save, private :: &
       RCSid = "$Id$"

!!***

contains

  !!****f* initialisation/initialise *
  !!
  !!  NAME
  !!   initialise
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Controls initialisation process for a run - reads in the
  !!   parameters and writes out information for the user, performs
  !!   various tedious setting up operations on pseudopotentials,
  !!   grids and matrices, sorts out the initial support functions
  !!   and finally gets a self-consistent Hamiltonian and potential.
  !!   At that point, the job is ready to go !
  !!  INPUTS
  !!
  !!
  !!  USES
  !!   common, datatypes, dimens, GenComms, initial_read, matrix_data
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   November 1998
  !!  MODIFICATION HISTORY
  !!   25/05/2001 dave
  !!    ROBODoc header, indenting and stripping subroutine calls
  !!   29/05/2001 dave
  !!    Stripped subroutine call
  !!   08/06/2001 dave
  !!    Added RCS Id and Log tags and GenComms for my_barrier
  !!   13/06/2001 dave
  !!    Changed call to set_up for init_pseudo
  !!   10/05/2002 dave
  !!    Added use statement for initial_read to get read_and_write
  !!   05/09/2002 mjg & drb
  !!    Added ionic data calls and uses
  !!   15:58, 25/09/2002 mjg & drb
  !!    Added careful flags for checking to see if atomic densities have been
  !!    specified, and whether the user wants to initialise the initial charge
  !!    density from the initial K, or from atomic densities, or hasn't told us
  !!    at all !
  !!   15:24, 27/02/2003 drb & tm
  !!    Moved flag_no_atomic_densities into density_module
  !!   11:03, 24/03/2003 drb
  !!    Simplified call to read_and_write, removed initial_phi
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new
  !!    matrix routines
  !!   2008/05/15 ast
  !!    Added some timers
  !!   2011/09/29 M. Arita
  !!    Statements were added for calculating only disperions
  !!   2011/12/11 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2016/03/15 13:54 dave
  !!    Removed restart_file from call to read_and_write and initial_phis (redundant)
  !!   2018/01/22 12:39 JST dave
  !!    Changes to find maximum angular momentum for PAOs and pseudopotentials (for factorial function)
  !!   2018/02/13 12:18 dave
  !!    Changes to new XC interface
  !!   2019/12/26 tsuyoshi
  !!    Removed flag_no_atomic_densities
  !!  SOURCE
  !!
  subroutine initialise(vary_mu, fixed_potential, mu, total_energy)

    use datatypes
    use numbers
    use global_module,     only: x_atom_cell, y_atom_cell, &
                                 z_atom_cell, ni_in_cell,  &
                                 flag_only_dispersion, flag_neutral_atom, &
                                 flag_atomic_stress, flag_heat_flux, &
                                 flag_full_stress, area_moveatoms, &
                                 atomic_stress, non_atomic_stress
    use GenComms,          only: inode, ionode, my_barrier, end_comms, &
                                 cq_abort
    use initial_read,      only: read_and_write
    use ionic_data,        only: get_ionic_data
    use memory_module,     only: init_reg_mem, reg_alloc_mem, type_dbl
    use group_module,      only: parts
    use primary_module,    only: bundle
    use cover_module,      only: make_cs, D2_CS
    use dimens,            only: r_dft_d2
    use DFT_D2
    use pseudo_tm_module,   only: make_neutral_atom
    use angular_coeff_routines, only: set_fact
    use maxima_module,          only: lmax_ps, lmax_pao
    use XC, only: init_xc
    
    implicit none

    ! Passed variables
    logical           :: vary_mu, find_chdens, fixed_potential
    character(len=40) :: output_file
    real(double)      :: mu
    real(double)      :: total_energy

    ! Local variables
    integer, parameter:: std_level_loc = 0
    type(cq_timer)    :: backtrace_timer
    logical           :: start, start_L
    logical           :: read_phi
    integer :: lmax_tot, stat

    call init_timing_system(inode)

    call init_reg_mem

!****lat<$
    call start_backtrace(t=backtrace_timer,who='initialise',&
         where=area,level=std_level_loc,echo=.true.)
!****lat>$
    lmax_pao = 0
    lmax_ps = 0
    ! Read input
    call read_and_write(start, start_L, inode, ionode, &
                        vary_mu, mu, find_chdens, read_phi)
    call init_xc
    ! IMPORTANT!!!!! No timers allowed before this point
    !                We need to know if the user wants them or not
    !                  before actually starting one
    call start_timer(tmr_std_initialisation)

    ! Only calculate the dispersion
    if (flag_only_dispersion) then
      call make_cs(inode-1, r_dft_d2, D2_CS, parts, bundle, ni_in_cell, &
                   x_atom_cell, y_atom_cell, z_atom_cell)
      call read_para_D2
      call dispersion_D2
      call end_comms()
      stop
    end if

    ! Call routines to read or make data for isolated ions
    call get_ionic_data(inode, ionode)
   
    ! 2019/Dec/26 TM
    ! Since flag_no_atomic_densities is always .false. now, 
    ! we don't need the following lines
    !if (flag_no_atomic_densities .and. (.not. find_chdens)) then
    !   if (inode == ionode) &
    !        write (io_lun, *) 'No initial charge density specified - &
    !                           &building from initial K'
    !   find_chdens = .true.
    !end if

    lmax_tot = lmax_pao+lmax_ps
    if(2*lmax_pao>lmax_tot) lmax_tot = 2*lmax_pao
    if(lmax_tot<8) lmax_tot = 8
    call set_fact(lmax_tot)
    if(flag_neutral_atom) call make_neutral_atom
    call set_up(find_chdens,std_level_loc+1)
    
    call my_barrier()

    call initial_phis(read_phi, start, fixed_potential, std_level_loc+1)

    ! ewald/screened_ion force and stress is computed in initial_H, so we
    ! need to allocate the atomic_stress array here - zamaan
    if (flag_heat_flux) then
      if (.not. flag_full_stress) then
        flag_full_stress = .true.
        if (inode==ionode) write(io_lun,'(2x,a)') &
          "WARNING: setting AtomMove.FullStress T for heat flux calculation"
      end if
      if (.not. flag_atomic_stress) then
        flag_atomic_stress = .true.
        if (inode==ionode) write(io_lun,'(2x,a)') &
          "WARNING: setting AtomMove.AtomicStress T for heat flux calculation"
      end if
    end if
    if (flag_atomic_stress) then
      if (.not. flag_neutral_atom) then
        call cq_abort("Atomic stress contributions not implemented for Ewald electrostatics (yet). Set General.NeutralAtom T")
      end if
      allocate(atomic_stress(3,3,ni_in_cell), STAT=stat)
      if (stat /= 0) &
        call cq_abort("Error allocating atomic_stress: ", ni_in_cell)
      call reg_alloc_mem(area_moveatoms, 3*3*ni_in_cell, type_dbl)
      atomic_stress = zero
      non_atomic_stress = zero
    end if

    call initial_H(start, start_L, find_chdens, fixed_potential, &
                   vary_mu, total_energy,std_level_loc+1)

    call stop_timer(tmr_std_initialisation)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='initialise',echo=.true.)
!****lat>$

    return
  end subroutine initialise
  !!***


  !!****f* initialisation/set_up *
  !!
  !!  NAME
  !!   set_up
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs various calls needed to set up the calculation
  !!   like pseudopotentials and grids
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   November 1998 ?
  !!  MODIFICATION HISTORY
  !!   25/05/2001 dave
  !!    Stripped various subroutine calls, added ROBODoc header
  !!    Stripped overall subroutine call
  !!   11/06/2001 dave
  !!    Added RCS Id and Log tags and GenComms
  !!   13/06/2001 dave
  !!    Changed to use pseudopotential_data for init_pseudo and
  !!    passed core_correction
  !!   15:58, 25/09/2002 mjg & drb
  !!    Added set_density for initialising charge density
  !!   11:51, 04/02/2003 drb
  !!    Removed cq_exit from GenComms use
  !!   13:11, 22/10/2003 mjg & drb
  !!    Changed set_ewald call to read appropriate flag and call old or new routine
  !!   11:59, 12/11/2004 dave
  !!    Changed to get nodes from GenComms not common
  !!   03/30/2011 19:22 M.Arita
  !!    Added the contribution from P.C.C.
  !!   2011/04/01 L.Tong
  !!    Added implementation for Spin polarisation: allocation for
  !!    density_dn and potential_dn Added numbers module dependence,
  !!    replacing for example 0.0_double to zero
  !!   2011/07/21 17:41 dave
  !!    Initialisation for cDFT
  !!   29/09/2011 16:24 M. Arita
  !!    Generate the covering set for DT-D2
  !!   2011/11/16 15:52 dave
  !!    Removing call to set up Extent (part of blip data changes)
  !!   2011/12/09 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   2013/08/20 M.Arita
  !!    Add call for make_glob2node (necessary for matrix reconstruction)
  !!   2013/12/03 M.Arita
  !!    Added call for immi_XL for XL-BOMD
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2015/11/09 17:12 dave with TM and NW (Mizuho)
  !!    - Added allocation for atomic density array, density_atom
  !!   2015/11/24 08:31 dave
  !!    - Removed old ewald calls
  !!   2015/11/30 17:09 dave
  !!    - Adding branches for neutral atom/ewald
  !!   2016/01/28 16:45 dave
  !!    Updated module name to ion_electrostatic
  !!   2016/01/29 15:00 dave
  !!    Removed prefix for ewald call
  !!   2016/09/16 17:00 nakata
  !!    Removed unused RadiusSupport
  !!   2017/06/22 11:04 dave
  !!    Adding diagonalisation initialisation calls
  !!   2017/08/29 jack baker & dave
  !!    Removed r_super_x references (redundant)
  !!   2018/01/22 12:41 JST dave
  !!    Adding check for maximum angular momentum for Bessel functions
  !!  SOURCE
  !!
  subroutine set_up(find_chdens,level)

    use datatypes
    use global_module,          only: iprint_init, flag_read_blocks,   &
                                      x_atom_cell, y_atom_cell,        &
                                      z_atom_cell, ni_in_cell,         &
                                      area_init, area_index,           &
                                      flag_Becke_weights,              &
                                      flag_pcc_global, flag_dft_d2,    &
                                      iprint_gen, flag_perform_cDFT,   &
                                      nspin, runtype,                  &
                                      glob2node, flag_XLBOMD,          &
                                      flag_neutral_atom, flag_diagonalisation
    use memory_module,          only: reg_alloc_mem, reg_dealloc_mem,  &
                                      type_dbl, type_int
    use group_module,           only: parts
    use primary_module,         only: bundle
    use cover_module,           only: BCS_parts, make_cs, make_iprim,  &
                                      send_ncover, D2_CS
    use mult_module,            only: immi
    use construct_module
    use matrix_data,            only: rcut, Lrange, Srange,            &
                                      mx_matrices, max_range
    use ion_electrostatic,      only: set_ewald, setup_screened_ion_interaction
    use atoms,                  only: distribute_atoms
    use dimens,                 only: n_grid_x, n_grid_y, n_grid_z,    &
                                      r_core_squared, r_h, &
                                      n_my_grid_points, r_dft_d2
    use fft_module,             only: set_fft_map, fft3
    use GenComms,               only: cq_abort, my_barrier, inode,     &
                                      ionode, gcopy
    use pseudopotential_data,   only: init_pseudo
    ! Troullier-Martin pseudos    15/11/2002 TM
    use pseudo_tm_module,       only: init_pseudo_tm
    use pseudopotential_common, only: pseudo_type, OLDPS, SIESTA,      &
                                      STATE, ABINIT, core_correction,  &
                                      pseudopotential
    ! Troullier-Martin pseudos    15/11/2002 TM
    use density_module,         only: set_atomic_density, density,            &
                                      density_scale, atomcharge,       &
                                      build_Becke_weights,             &
                                      build_Becke_charges,             &
                                      set_density_pcc, density_pcc,    &
                                      density_atom
    use block_module,           only: nx_in_block,ny_in_block,         &
                                      nz_in_block, n_pts_in_block,     &
                                      set_blocks_from_new,             &
                                      set_blocks_from_old,             &
                                      set_domains, n_blocks
    use grid_index,             only: grid_point_x, grid_point_y,      &
                                      grid_point_z, grid_point_block,  &
                                      grid_point_position
    use primary_module,         only: domain
    use group_module,           only: blocks
    use io_module,              only: read_blocks
    use functions_on_grid,      only: associate_fn_on_grid
    use potential_module,       only: potential
    use maxima_module,          only: maxngrid, lmax_ps, lmax_pao
    use species_module,         only: n_species
    use angular_coeff_routines, only: set_fact, set_prefac, set_prefac_real
    use numbers,                only: zero
    use cDFT_module,            only: init_cdft
    use DFT_D2,                 only: read_para_D2
    use input_module,           ONLY: leqi
    use UpdateInfo,             ONLY: make_glob2node
    use XLBOMD_module,          ONLY: immi_XL
    use DiagModule,             only: init_blacs_pg, init_scalapack_format
    
    implicit none

    ! Passed variables
    logical          :: find_chdens
    integer,optional :: level

    ! Local variables
    complex(double_cplx), allocatable, &
         dimension(:) :: chdenr
    integer           :: i, stat, spec, lmax_tot
    real(double)      :: rcut_BCS  !TM 26/Jun/2003
    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='set_up',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    ! Set organisation of blocks of grid-points.
    ! set_blocks determines the number of blocks on this node,
    ! and makes a list of these blocks.
    if (flag_read_blocks) then
       call read_blocks(blocks)
    else
       call set_blocks_from_old(n_grid_x, n_grid_y, n_grid_z)
    endif
    call set_blocks_from_new()
    ! Allocate ?
    !call set_blocks(inode, ionode)
    if (inode == ionode .and. iprint_init > 1) &
         write(io_lun,*) 'Completed set_blocks()'
    n_my_grid_points = n_blocks*n_pts_in_block
    ! allocate(grid_point_x(n_my_grid_points),&
    !      grid_point_y(n_my_grid_points),&
    !      grid_point_z(n_my_grid_points), &
    !      grid_point_block(n_my_grid_points),&
    !      grid_point_position(n_my_grid_points),STAT=stat)
    allocate(grid_point_x(maxngrid), grid_point_y(maxngrid),     &
             grid_point_z(maxngrid), grid_point_block(maxngrid), &
             grid_point_position(maxngrid), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating grid_point variables: ", &
                       maxngrid, stat)
    call reg_alloc_mem(area_index, 5 * maxngrid, type_int)
    ! Construct list of grid-points in each domain (i.e. grid-points belonging
    ! to each node). In present version, grid-points are organised into
    ! blocks, with each node responsible for a cluster of blocks.
    call set_domains(inode)
    allocate(density(maxngrid,nspin), potential(maxngrid,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating density and potential: ", &
                       maxngrid, stat)
    call reg_alloc_mem(area_index, 2 * nspin * maxngrid, type_dbl)
    if(flag_neutral_atom) then ! Allocate atomic density array
       allocate(density_atom(maxngrid), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating density_atom: ", &
            maxngrid, stat)
       call reg_alloc_mem(area_index, maxngrid, type_dbl)
    end if
    allocate(density_scale(nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating density_scale: ", nspin, stat)
    call reg_alloc_mem(area_index, nspin, type_dbl)
    allocate(pseudopotential(maxngrid), STAT=stat)
    if (stat/=0) &
         call cq_abort("Error allocating grids: ", maxngrid, stat)
    call reg_alloc_mem(area_index, maxngrid, type_dbl)
    ! extra local potential for second spin channel for spin polarised calculation
    call my_barrier()
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Completed set_domains()'

    ! Sorts out which processor owns which atoms
    call distribute_atoms(inode, ionode)
    call my_barrier
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Completed distribute_atoms()'
    ! Create a covering set
    call my_barrier
    !Define rcut_BCS  !TM 26/Jun/2003
    !rcut_BCS= 2.0_double*rcut(Lrange)+rcut(Srange)
    !do i=1, mx_matrices
    !   if(rcut_BCS < rcut(i)) rcut_BCS= rcut(i)
    !enddo !  i=1, mx_matrices
    rcut_BCS = rcut(max_range)
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) ' rcut for BCS_parts =', rcut_BCS

    call make_cs(inode-1, rcut_BCS, BCS_parts, parts, bundle, &
                 ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell)
    call my_barrier
    call make_iprim(BCS_parts, bundle, inode-1)
    call send_ncover(BCS_parts, inode)
    call my_barrier
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Made covering set for matrix multiplications'

    ! Create all of the indexing required to perform matrix multiplications
    ! at a later point. This routine also identifies all the density
    ! matrix range interactions and hamiltonian range interactions
    ! associated with any atom being handled by this processor.
    call immi(parts, bundle, BCS_parts, inode)
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Completed immi()'
    if (flag_XLBOMD) call immi_XL(parts,bundle,BCS_parts,inode)

    ! set up all the data block by block for atoms overlapping any
    ! point on block and similar
    !call set_blocks_from_old(                                       &
    !     in_block_x, in_block_y, in_block_z, n_grid_x, n_grid_y, n_grid_z, &
    !     inode, ionode, NODES)
    call setgrid(inode-1, r_core_squared, r_h)

    call my_barrier()
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Completed set_grid()'

    call associate_fn_on_grid
    call my_barrier()
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Completed associate_fn_on_grid()'

    ! The FFT requires the data to be reorganised into columns parallel to
    ! each axis in turn. The data for this organisation is help in map.inc,
    ! and is initialised by set_fft_map.
    !
    ! The FFT calls themselves require value tables, which are held in
    ! ffttable.inc, and are initialised by calling fft3 with isign=0.
    !
    ! this is just for initialising fft, no need to do it again for
    ! other spin component.
    call set_fft_map()
    density = zero
    ! chdenr is the temporary storage for Fourier transform of density
    allocate(chdenr(maxngrid), STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating chdenr: ", maxngrid, stat)
    call reg_alloc_mem(area_init, maxngrid, type_dbl)
    call fft3(density(:,1), chdenr, maxngrid, 0)
    deallocate(chdenr, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating chdenr: ", maxngrid, stat)
    call reg_dealloc_mem(area_init, maxngrid, type_dbl)
    call my_barrier()
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Completed fft init'

    ! Initialise the routines to calculate ion-ion interactions
    if(flag_neutral_atom) then
       call setup_screened_ion_interaction
       call my_barrier
       if (inode == ionode .and. iprint_init > 1) &
            write (io_lun, *) 'Completed setup_ion_interaction()'
    else
       ! set up the Ewald sumation: find out how many superlatices
       ! in the real space sum and how many reciprocal latice vectors in the
       ! reciprocal space sum are needed for a given energy tolerance. 
       call set_ewald(inode,ionode)
       call my_barrier
       if (inode == ionode .and. iprint_init > 1) &
            write (io_lun, *) 'Completed set_ewald()'
    end if
    ! +++

    ! Generate D2CS
    if (flag_dft_d2) then
      if ((inode == ionode) .and. (iprint_gen > 1) ) &
           write (io_lun, '(/1x,"The dispersion is considered in the &
                           &DFT-D2 level.")')
      call make_cs(inode-1, r_dft_d2, D2_CS, parts, bundle, ni_in_cell, &
                   x_atom_cell, y_atom_cell, z_atom_cell)
      if ( (inode == ionode) .and. (iprint_gen > 1) ) then
        write (io_lun, '(/8x,"+++ D2_CS%ng_cover:",i10)')       &
              D2_CS%ng_cover
        write (io_lun, '(8x,"+++ D2_CS%ncoverx, y, z:",3i8)')   &
              D2_CS%ncoverx, D2_CS%ncovery, D2_CS%ncoverz
        write (io_lun, '(8x,"+++ D2_CS%nspanlx, y, z:",3i8)')   &
              D2_CS%nspanlx, D2_CS%nspanly, D2_CS%nspanlz
        write (io_lun, '(8x,"+++ D2_CS%nx_origin, y, z:",3i8)') &
              D2_CS%nx_origin, D2_CS%ny_origin, D2_CS%nz_origin
      end if
      call read_para_D2
      if (inode == ionode) then                               !! DEBUG !!
         write (io_lun, '(a, f10.5)') &                       !! DEBUG !!
               "Sbrt: make_cs for DFT-D2, the cutoff is ", &  !! DEBUG !!
               r_dft_d2                                       !! DEBUG !!
      end if                                                   !! DEBUG !!
   end if

   ! external potential - first set up angular momentum bits
   !call set_fact
   lmax_tot = lmax_pao+lmax_ps
   if(lmax_tot<2*lmax_pao) lmax_tot = 2*lmax_pao
   if(lmax_tot<8) lmax_tot = 8
   call set_prefac(lmax_tot+1)
   call set_prefac_real(lmax_tot+1)
    !  TM's pseudo or not : 15/11/2002 TM
   select case (pseudo_type)
   case(OLDPS)
      call init_pseudo(core_correction)
   case(SIESTA)
      call init_pseudo_tm(core_correction)
   case(ABINIT)
      call init_pseudo_tm(core_correction)
   end select

   ! For P.C.C.
   if (flag_pcc_global) then
      allocate(density_pcc(maxngrid), STAT=stat)
      call set_density_pcc
      if (stat/=0) then
         call cq_abort("ERROR allocating grids: ", maxngrid, stat)
         call reg_alloc_mem(area_init, maxngrid, type_dbl)
      end if
   end if

   ! Get the table showing the relation between atoms & processors
     allocate (glob2node(ni_in_cell), STAT=stat)
     if (stat.NE.0) call cq_abort('Error allocating glob2node: ', ni_in_cell)
     if (inode.EQ.ionode) call make_glob2node
     call gcopy(glob2node,ni_in_cell)

   ! Create initial density from superposition of atomic densities
   if (.not. find_chdens) then
      call set_atomic_density(.true.)  ! Set density and NA atomic density
   else if( flag_neutral_atom ) then
      call set_atomic_density(.false.) ! Need atomic density for neutral atom potential
   end if

   if (flag_perform_cDFT) then
      call init_cdft
   end if
   if (flag_Becke_weights) then
      allocate(atomcharge(ni_in_cell,nspin), STAT=stat)
      if (stat /= 0) &
           call cq_abort("Error allocating atomcharge: ", &
                          ni_in_cell, stat)
      call reg_alloc_mem(area_init, nspin * ni_in_cell, type_dbl)
      call build_Becke_weights
      call build_Becke_charges(atomcharge, density, maxngrid)
   end if
   if (inode == ionode .and. iprint_init > 1) &
        write (io_lun, *) 'Done init_pseudo '

   if(flag_diagonalisation) then
      call init_blacs_pg
      call init_scalapack_format
   end if
!****lat<$
   call stop_backtrace(t=backtrace_timer,who='set_up',echo=.true.)
!****lat>$

   return
 end subroutine set_up
 !!***


 !!****f* initialisation/initial_phis *
 !!
 !!  NAME
 !!   initial_phis
 !!  USAGE
 !!
 !!  PURPOSE
 !!   Provides initial values for the blip coefficients
 !!   representing the support functions. Two ways of
 !!   doing this are provided, these ways being sp/ecified
 !!   by the character-valued variable init_blip_flag, as follows:
 !!
 !!    init_blip_flag = 'gauss': blips coeffs set as values of a Gaussian
 !!      function at the distance of each blip from the atom position.
 !!      Can only be used if there are four supports on every atom,
 !!      in which case the supports have the form of an s-function
 !!      and three p-functions. The Gaussian exponents are alph and beta
 !!    (atomic units). This way of intitiating blip coefficients is
 !!      a relic of the very early history of Conquest, and its
 !!      use is strongly discouraged. Don't put init_blip_flag = 0
 !!      unless you have thought very carefully about why you
 !!      need to do this.
 !!
 !!    init_blip_flag = 'pao': blip coeffs set to give best fit to
 !!      support functions represented as a given linear combination
 !!      of pseudo-atomic orbitals (PAO's). The PAO data itself
 !!      is held in the module pao_format, and the PAO representations
 !!      of support functions in the module support_spec_format.
 !!
 !!   These two ways of initiating the blip coeffs are implemented
 !!   by calls to the subroutines:
 !!
 !!    init_blip_flag = 'gauss': gauss2blip
 !!    init_blip_flag = 'pao': make_blips_from_paos
 !!
 !!  INPUTS
 !!
 !!
 !!  USES
 !!
 !!  AUTHOR
 !!   D.R.Bowler
 !!  CREATION DATE
 !!   Probably late 1998
 !!  MODIFICATION HISTORY
 !!   28/05/2001 dave
 !!    ROBODoc header, stripped normalise_support call
 !!
 !!   11/06/2001 dave
 !!    Added RCS Id and Log tags and GenComms
 !!
 !!   13/06/02 mjg
 !!    Rewritten to allow 2 ways of initiation (see above PURPOSE)
 !!   15:17, 04/02/2003 drb
 !!    Added call to do on-site S matrix elements analytically
 !!   11:04, 24/03/2003 drb
 !!    Removed call to interp_phi
 !!   15:44, 08/04/2003 drb
 !!    Added another blip_to_support and get_matrix_elements after
 !!    normalisation to ensure correctly
 !!    normalised blips used
 !!   14:41, 29/07/2003 drb
 !!    Small, safe changes: no call to get_onsite_S until fixed, and
 !!    normalisation based on NUMERICAL integration
 !!   08:32, 2003/09/22 dave
 !!    Changed to not do blip stuff for PAO basis set
 !!   14:05, 22/09/2003 drb
 !!    Added test for division by zero
 !!   11:02, 23/09/2003 drb
 !!    Bug fix to write statement
 !!   12:58, 31/10/2003 drb
 !!    Added NLPF call for generating SP matrix from PAOs
 !!   11:40, 12/11/2004 dave
 !!    Removed inappropriate common variables
 !!   10:09, 13/02/2006 drb
 !!    Removed all explicit references to data_ variables and rewrote
 !!    in terms of new
 !!    matrix routines
 !!   2006/03/06 05:17 dave
 !!    Tidied call and passed variables
 !!   2006/09/13 07:57 dave
 !!    Changed to get number of coefficients for blips from
 !!    support_function structure
 !!   2007/05/01 08:31 dave
 !!    Changed start_blip into read_option for unified naming
 !!   2011/11/28 07:55 dave
 !!    Added call for representing NLPFs with blips
 !!   2015/06/08 lat
 !!    Added experimental backtrace
 !!   2016/03/15 13:55 dave
 !!    Removed restart_file (redundant) and mu
 !!   2016/07/29 18:30 nakata
 !!    Renamed supports_on_atom -> blips_on_atom
 !!   2016/08/08 15:30 nakata
 !!    Renamed supportfns -> atomfns
 !!   2017/02/21 16:00 nakata
 !!    Removed unused get_support_pao_rep
 !!   2019/12/26 tsuyoshi
 !!    Removed unused find_chdens
 !!  SOURCE
 !!
 subroutine initial_phis(read_phi, start, fixed_potential, level)

    use datatypes
    use blip,                        only: init_blip_flag, make_pre,   &
                                           set_blip_index, gauss2blip
    use blip_grid_transform_module,  only: blip_to_support_new
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use dimens,                      only: grid_point_volume, r_h
    !use fdf,                        only : fdf_boolean
    use GenComms,                    only: cq_abort, my_barrier,       &
                                           gcopy, inode, ionode
    use global_module,               only: iprint_init,                &
                                           flag_basis_set, blips,      &
                                           PAOs, flag_onsite_blip_ana, &
                                           flag_analytic_blip_int, flag_neutral_atom
    use matrix_data,                 only: Srange, mat
    use numbers,                     only: zero, RD_ERR, one
    use pao2blip,                    only: make_blips_from_paos
    use primary_module ,             only: bundle
    use set_bucket_module,           only: rem_bucket
    use species_module,              only: n_species, nsf_species
    use io_module,                   only: grab_blip_coeffs,           &
                                           dump_matrix
    use mult_module,                 only: return_matrix_value, matS
    ! Temp
    use S_matrix_module,             only: get_onsite_S, get_S_matrix
    use make_rad_tables,             only: gen_rad_tables,             &
         gen_nlpf_supp_tbls, gen_paoNApao_tbls, &
         gen_napf_supp_tbls
    use angular_coeff_routines,      only: make_ang_coeffs, set_fact,  &
                                           set_prefac, set_prefac_real
    use read_support_spec,           only: read_support
    use functions_on_grid,           only: atomfns
    use support_spec_format,         only: blips_on_atom,              &
                                           coefficient_array,          &
                                           coeff_array_size,           &
                                           read_option
    use input_module,                only: leqi
    use nlpf2blip,                   only: make_blips_from_nlpfs
    use pseudopotential_common,      only: flag_neutral_atom_projector

    implicit none

    ! Passed variables
    logical           :: read_phi, start, fixed_potential
    integer, optional :: level

    ! Local variables
    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level
    real(double)      :: factor
    integer           :: isf, np, ni, iprim, n_blip
    integer           :: n_run, spec, this_nsf
    integer           :: spin_SF

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='initial_phis',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    ! Used by pseudopotentials as well as PAOs
    if (flag_basis_set==blips) then

       spin_SF = 1
!       if (flag_SpinDependentSF) spin_SF = spin   ! SpinDependentSF is not available with blips at present

       if (inode==ionode) &
            write(io_lun,fmt='(10x,"Using blips as basis set for &
                               &support functions")')
       call set_blip_index(inode, ionode)

       call my_barrier
       if((inode == ionode) .and. (iprint_init > 1)) &
            write(io_lun,*) 'initial_phis: completed set_blip_index()'

       !if((inode == ionode).and.(iprint_init >= 0)) then
       !   write(unit=io_lun,fmt='(/" initial_phis: n_species:",i3)') n_species
       !   write(unit=io_lun,fmt='(/" initial_phis: r_h:",f12.6)') r_h
       !   write(unit=io_lun,fmt='(/" initial_phis: support_grid_spacing:"&
       !        &,f12.6)') support_grid_spacing
       !end if
       if(.NOT.read_option) then
          if(leqi(init_blip_flag,'gauss')) then
             call gauss2blip
          else if(leqi(init_blip_flag,'pao')) then
             call make_blips_from_paos(inode,ionode,n_species)
          else
             call cq_abort('initial_phis: init_blip_flag no good')
          end if
       end if
       call my_barrier
       if((inode == ionode).AND.(iprint_init > 0)) then
          write(unit=io_lun,&
                fmt='(10x,"initial_phis: initial blip coeffs created")')
       end if
       ! Create blip representation of non-local projector functions
       if(flag_analytic_blip_int) call make_blips_from_nlpfs

       ! Length scale Preconditioning.
       call make_pre(inode, ionode)
       ! Restart files
       ! Tweak DRB 2007/03/34 Remove need for start
       !if (start.and.(.not.start_blips)) then
       if (read_option) then
          if(inode==ionode.AND.iprint_init>0) &
               write(io_lun,fmt='(10x,"Loading blips")')
          call grab_blip_coeffs(coefficient_array,coeff_array_size, inode)
       end if
       ! Normalisation
       call blip_to_support_new(inode-1, atomfns)
       if((inode == ionode).AND.(iprint_init > 1)) then
          write(unit=io_lun,&
                fmt='(10x,"initial_phis: completed blip_to_support()")')
       end if
       if (start .or. (.NOT.read_option)) then
          n_run = 0
          !     call normalise_support(support, inode, ionode,&
          !          NSF, SUPPORT_SIZE)
          !     call blip_to_support_new(inode-1, support, data_blip, &
          !          NSF, SUPPORT_SIZE, MAX_N_BLIPS)
          !     write(io_lun,*) 'S matrix for normalisation on Node= ',inode
          call get_matrix_elements_new(inode-1,rem_bucket(1),matS(spin_SF),&
                                       atomfns,atomfns)
          ! Do the onsite elements analytically
          if(flag_onsite_blip_ana) then
             iprim=0
             do np=1,bundle%groups_on_node
                if(bundle%nm_nodgroup(np) > 0) then
                   do ni=1,bundle%nm_nodgroup(np)
                      iprim=iprim+1
                      spec = bundle%species(iprim)
                      this_nsf = nsf_species(spec)
                      call get_onsite_S(blips_on_atom(iprim), matS(spin_SF),&
                                        np, ni, iprim, this_nsf, spec)
                   end do
                end if
             end do
          end if
          iprim=0
          call start_timer(tmr_std_matrices)
          do np=1,bundle%groups_on_node
             if(bundle%nm_nodgroup(np) > 0) then
                do ni=1,bundle%nm_nodgroup(np)
                   iprim=iprim+1
                   do isf=1,mat(np,Srange)%ndimi(ni)
                      factor = return_matrix_value(matS(spin_SF),np,ni,iprim,0,isf,isf,1)
                      if(factor>RD_ERR) then
                         factor=one/sqrt(factor)
                      else
                         factor = zero
                      end if
                      do n_blip=1,blips_on_atom(iprim)%supp_func(isf)%ncoeffs
                         blips_on_atom(iprim)%&
                              supp_func(isf)%coefficients(n_blip) = &
                              factor*blips_on_atom(iprim)%&
                              supp_func(isf)%coefficients(n_blip)
                      enddo ! n_blip
                   enddo ! isf
                enddo ! ni
             endif ! if the partition has atoms
          enddo ! np
          call stop_timer(tmr_std_matrices)
          call my_barrier()
          if(inode==ionode.AND.iprint_init>1) &
               write(io_lun,fmt='(10x,"Completed normalise_support")')
          call blip_to_support_new(inode-1, atomfns)
          call get_matrix_elements_new(inode-1,rem_bucket(1),matS(spin_SF),&
                                       atomfns,atomfns)
          ! Do the onsite elements analytically
          if(flag_onsite_blip_ana) then
             iprim=0
             do np=1,bundle%groups_on_node
                if(bundle%nm_nodgroup(np) > 0) then
                   do ni=1,bundle%nm_nodgroup(np)
                      iprim=iprim+1
                      spec = bundle%species(iprim)
                      this_nsf = nsf_species(spec)
                      call get_onsite_S(blips_on_atom(iprim), matS(spin_SF),&
                                        np, ni, iprim, this_nsf, spec)
                   end do
                end if
             end do
          end if
          !call dump_matrix("NS",matS(1),inode)
       else
          if(inode==ionode.AND.iprint_init>1) &
               write(io_lun,fmt='(10x,"Skipped normalise_support")')
       end if
    else if(flag_basis_set==PAOs) then
       call gen_rad_tables(inode,ionode)
       call gen_nlpf_supp_tbls(inode,ionode)
       ! Neutral atom Projector functions
       if( flag_neutral_atom .and. flag_neutral_atom_projector ) then
          call gen_paoNApao_tbls(inode,ionode)
          call gen_napf_supp_tbls(inode,ionode)
       end if
       call make_ang_coeffs
       if(inode==ionode) &
            write(io_lun,&
                  fmt='(10x,"Using PAOs as basis set for support functions")')
       !call make_pre_paos
       if((inode == ionode).and.(iprint_init >0)) then
          write(unit=io_lun,fmt='(10x,"initial_phis: n_species:",i3)') n_species
          write(unit=io_lun,fmt='(10x,"initial_phis: r_h:",f12.6)') r_h
       end if
       ! We don't need a PAO equivalent of blip_to_support here: this
       ! is done by get_S_matrix_PAO
    end if

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='initial_phis',echo=.true.)
!****Lat>$

    return
  end subroutine initial_phis
  !!***


  !!****f* initialisation/initial_H *
  !!
  !!  NAME
  !!   initial_H
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Makes an initial, self-consistent Hamiltonian (and potential)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!   datatypes, DMMin, ion_electrostatic, GenComms, global_module, logicals,
  !!   matrix_data, mult_module, numbers, SelfCon, S_matrix_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   14/05/99
  !!  MODIFICATION HISTORY
  !!   18/05/2001 dave
  !!    Stripped down call to new_SC_potl
  !!   25/05/2001 dave
  !!    Used S_matrix_module for get_S_matrix
  !!    Shortened overall subroutine call
  !!   08/06/2001 dave
  !!    Used GenComms and added RCS Id and Log tags
  !!   13/06/2001 dave
  !!    Removed get_core_correction (a routine which does
  !!    not need to be passed a pig)
  !!   17/06/2002 dave
  !!    Improved headers slightly and changed call to
  !!    main_matrix_multiply to remove unnecessary work
  !!   15:38, 04/02/2003 drb
  !!    Added all sorts of calls for testing forces: left in for now,
  !!    but commented out, so that we can see the kinds of things that
  !!    need to be done
  !!   14:33, 26/02/2003 drb
  !!    Added appropriate calls based on whether or not we're doing
  !!    self-consistency
  !!   15:55, 27/02/2003 drb & tm
  !!    Changed call to get_H_matrix to turn off charge generation
  !!   08:34, 2003/03/12 dave
  !!    Added call to get_energy after FindMinDM
  !!   12:01, 30/09/2003 drb
  !!    Added force-testing call, tidied
  !!   13:12, 22/10/2003 mjg & drb
  !!    Added old/new ewald call
  !!   09:13, 11/05/2005 dave
  !!    Stopped FindMinDM if restart_L flag set
  !!   2005/07/11 10:24 dave
  !!    Removed redeclaration of restart flags
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2011/07/18 L.Tong
  !!    Added Spin polarisation
  !!    Added grab matrix for L_dn if spin polarised
  !!    Implemented Spin polarisation for calcuating K and Phi
  !!   2011/09/19 L.Tong
  !!    Removed some dependence on number_of_bands as no longer required
  !!    by some updated subroutines
  !!   29/09/2011 16:26 M. Arita
  !!    Calculate the dispersion in the DFT-D2 level
  !!   2011/12/09 L.Tong
  !!    Removed redundant parameter number_of_bands
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   - Removed redundant input parameter real(double) mu
  !!   2013/07/05 dave
  !!    Changed L matrix reload to use up/dn names for spin
  !!   2013/08/20 M.Arita
  !!    Added calls for reading global data and reconstructing L-matrix
  !!    and T-matrix
  !!   2013/12/03 M.Arita
  !!    Added calls for reading X-matrix files for XL-BOMD
  !!   2014/02/03 michi
  !!    Bug fix for changing the number of processors at the sequential job
  !!   2015/06/08 lat
  !!    Added experimental backtrace
  !!   2015/11/24 08:32 dave
  !!    Removed old ewald calls
  !!   2015/11/30 17:10 dave
  !!    Added branches for neutral atom/ewald
  !!   2016/01/28 16:44 dave
  !!    Updated module name to ion_electrostatic
  !!   2016/01/29 15:00 dave
  !!    Removed prefix for ewald call
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2016/12/19 17:15 nakata
  !!    Removed unused flag_vary_basis
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2017/04/10 dave
  !!    Small bug fix: transpose SF coefficients before transforming (fixes issue 26)
  !!   2017/05/09 dave
  !!    Removed restriction spin and on L-matrix update 
  !!   2017/10/11 dave
  !!    Bug fix: changed call to get_electronic_density to use inode, ionode (not ionode, ionode)
  !!    Removed restriction spin and on L-matrix update
  !!   2017/11/10 15:24 dave
  !!    Added allocation of glob2node_old
  !!   2018/01/22 tsuyoshi (with dave)
  !!    Updates for new atom movement routines
  !!   2018/02/08 tsuyoshi (with dave)
  !!    Bug fix for reading K after atom movement
  !!   2018/11/13 17:30 nakata
  !!    Changed matT to be spin_SF dependent
  !!   2019/10/24 11:52 dave
  !!    Changed function calls to FindMinDM
  !!   2019/11/14 tsuyoshi 
  !!    The part setting atom_coord_diff is replaced by subroutine call to "set_atom_coord_diff".
  !!    Removed n_proc_old and glob2node_old
  !!   2019/11/18 tsuyoshi 
  !!    Removed flag_MDold
  !!  SOURCE
  !!
  subroutine initial_H(start, start_L, find_chdens, fixed_potential, &
       vary_mu, total_energy,level)

    use datatypes
    use numbers
    use logicals
    use mult_module,         only: LNV_matrix_multiply, matL, matphi, &
         matT, T_trans, L_trans, LS_trans,  &
         SFcoeff_trans, matK,               &
         matrix_scale, matSFcoeff, matSFcoeff_tran, matrix_transpose
    use SelfCon,             only: new_SC_potl
    use global_module,       only: iprint_init, flag_self_consistent, &
         flag_basis_set, blips, PAOs,       &
         restart_LorK,                      &
         restart_rho, flag_test_forces,     &
         flag_dft_d2, nspin, spin_factor,   &
         flag_MDcontinue,                   &
         restart_T,glob2node,               &
         MDinit_step,ni_in_cell,            &
         flag_XLBOMD, flag_dissipation,     &
         flag_propagateX, flag_propagateL, restart_X, &
         flag_exx, exx_scf, flag_out_wf, wf_self_con, &
         flag_write_DOS, flag_neutral_atom, &
         atomf, sf, flag_LFD, nspin_SF, flag_diagonalisation, &
         atom_coord, atom_coord_diff, rcellx, rcelly, rcellz, &
         ne_in_cell
    use ion_electrostatic,   only: ewald, screened_ion_interaction
    use S_matrix_module,     only: get_S_matrix
    use GenComms,            only: my_barrier,end_comms,inode,ionode, &
         cq_abort,gcopy
    use DMMin,               only: correct_electron_number, FindMinDM
    use H_matrix_module,     only: get_H_matrix
    use energy,              only: get_energy
    use test_force_module,   only: test_forces
    use io_module,           only: grab_matrix, grab_charge
    !use DiagModule,          only: diagon
    use density_module,      only: get_electronic_density, density
    use functions_on_grid,   only: atomfns, H_on_atomfns
    use dimens,              only: n_my_grid_points
    use maxima_module,       only: maxngrid
    use minimise,            only: sc_tolerance, L_tolerance,         &
         n_L_iterations, expected_reduction
    use DFT_D2,              only: dispersion_D2
    use matrix_data,         ONLY: Lrange,Trange,LSrange,SFcoeff_range,Hrange
    use store_matrix,        ONLY: matrix_store_global, grab_InfoMatGlobal, grab_matrix2, &
         InfoMatrixFile, set_atom_coord_diff
    use UpdateInfo,          ONLY: make_glob2node,Matrix_CommRebuild, Report_UpdateMatrix
    use XLBOMD_module,       ONLY: grab_XXvelS,grab_Xhistories
    use support_spec_format, only: read_option
    use multisiteSF_module,  only: initial_SFcoeff

    implicit none

    ! Passed variables
    logical           :: vary_mu, find_chdens, fixed_potential
    logical           :: start, start_L
    real(double)      :: total_energy
    integer, optional :: level

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level
    logical        :: reset_L, record, rebuild_KE_NL, build_X ! charge, store
    integer        :: force_to_test, stat, nfile, symm
    real(double)   :: electrons_tot, bandE
    real(double), dimension(nspin) :: electrons, energy_tmp
    integer        :: spin_SF
    !H_trans is not prepared. If we need to symmetrise K, we need H_trans
    integer        :: H_trans = 1

    type(matrix_store_global) :: InfoGlob
    type(InfoMatrixFile),pointer :: Info(:)

    ! Dummy vars for MMM

    !****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='initial_H',&
         where=area,level=backtrace_level,echo=.true.)
    !****lat>$

    ! (0) Get the global information
    !      --> Fetch and distribute date on old job

    if (flag_MDcontinue.or. &
         restart_LorK.or. &
         restart_T   .or. &
         read_option  ) then
       if (inode.eq.ionode) write (io_lun,*) "Get global info to load matrices"
       if (inode.eq.ionode) call make_glob2node
       call gcopy(glob2node, ni_in_cell)
       call grab_InfoMatGlobal(InfoGlob,MDinit_step)  
       call set_atom_coord_diff(InfoGlob)
       MDinit_step = InfoGlob%MDstep
     
      ! 2019Dec26 TM
      ! Now, since the default value of find_chdens is 'true' if we set restart_LorK 
      ! as .true., the following line should be commented out. 
      !  Even when we set restart_LorK, we should be able to choose the option 
      ! where find_chdens = .false. (initial charge = atomic charge)
      ! But.. since this change will affect the result, we will issue this change later.
      !
       if(restart_LorK) find_chdens=.true.  ! 2018JFeb12 TM 

       call my_barrier()
    endif

    ! (0) If we use PAOs and contract them, prepare SF-PAO coefficients here
    if (atomf.ne.sf) then
       if (restart_rho) then
          ! Read density from input files here to make SF-PAO coefficients
          if (nspin == 2) then
             call grab_charge(density(:,1), n_my_grid_points, inode, spin=1)
             call grab_charge(density(:,2), n_my_grid_points, inode, spin=2)
          else
             call grab_charge(density(:,1), n_my_grid_points, inode)
             density = density / spin_factor
          end if
       endif
       do spin_SF = 1, nspin_SF
          call matrix_scale(zero,matSFcoeff(spin_SF))
       enddo
       if (read_option) then
          if (inode == ionode) write (io_lun,*) 'Read supp_pao coefficients from SFcoeff files'
          call grab_matrix2('SFcoeff',inode,nfile,Info,InfoGlob,index=0,n_matrix=nspin_SF)
          call my_barrier()
          call Matrix_CommRebuild(InfoGlob,Info,SFcoeff_range,SFcoeff_trans,matSFcoeff,nfile,n_matrix=nspin_SF)

          ! Added DRB 2017/04/10 to fix issue 26: transpose required before transformation can occur
          ! Transpose
          do spin_SF = 1,nspin_SF
             call matrix_scale(zero,matSFcoeff_tran(spin_SF))
             call matrix_transpose(matSFcoeff(spin_SF), matSFcoeff_tran(spin_SF))
          enddo
       else
          ! make SF-PAO coefficients
          call initial_SFcoeff(.true., .true., fixed_potential, .true.)
       endif
       if (inode == ionode .and. iprint_init > 1) write (io_lun, *) 'Got SFcoeff'
       call my_barrier
    endif
!!$
!!$
!!$
!!$
    total_energy = zero
    ! If we're vary PAOs, allocate memory
    ! (1) Get S matrix
    if (restart_T) then
       call grab_matrix2('T',inode,nfile,Info,InfoGlob,index=0,n_matrix=nspin_SF)
       call my_barrier()
       call Matrix_CommRebuild(InfoGlob,Info,Trange,T_trans,matT,nfile,symm,n_matrix=nspin_SF)
    endif
    if (flag_LFD .and. .not.read_option) then
       ! Spao was already made in sub:initial_SFcoeff
       call get_S_matrix(inode, ionode, build_AtomF_matrix=.false.)
    else
       call get_S_matrix(inode, ionode)
    endif
    if (inode == ionode .and. iprint_init > 1) write (io_lun, *) 'Got S'
    call my_barrier
!!$
!!$
!!$
!!$
    ! (2) Make an inital estimate for the density matrix, L, which is an
    !     approximation to L = S^-1. Then use correct_electron_number()
    !     to modify L so that the electron number is correct.
    if (.not. flag_diagonalisation .and. find_chdens .and. (start .or. start_L)) then
       call initial_L()
       call my_barrier()
       if (inode == ionode .and. iprint_init > 1) &
            write (io_lun, *) 'Got L  matrix'
       if (vary_mu) then
          ! This cannot be timed within the routine
          call start_timer(tmr_std_densitymat)
          call correct_electron_number
          call stop_timer(tmr_std_densitymat)
       end if
    end if
    if (restart_LorK) then
       if(.not.flag_diagonalisation) then
          call grab_matrix2('L',inode,nfile,Info,InfoGlob,index=0,n_matrix=nspin)
          call my_barrier()
          call Matrix_CommRebuild(InfoGlob,Info,Lrange,L_trans,matL,nfile,symm,n_matrix=nspin)
          if (inode == ionode .and. iprint_init > 1) write (io_lun, *) 'Grabbed L  matrix'
       else
          call grab_matrix2('K',inode,nfile,Info,InfoGlob,index=0,n_matrix=nspin)
          call my_barrier()
          call Matrix_CommRebuild(InfoGlob,Info,Hrange,H_trans,matK,nfile,n_matrix=nspin)
          if (inode == ionode .and. iprint_init > 1) write (io_lun, *) 'Grabbed K  matrix'
          !DEBUG call Report_UpdateMatrix("Kmat")  
       end if
    end if
    ! XL-BOMD
    if (restart_X) then
       if (flag_propagateX) then
          call grab_XXvelS(LSrange,LS_trans,InfoGlob)
          if (flag_dissipation) call grab_Xhistories(LSrange,LS_trans,InfoGlob)
       elseif (flag_propagateL) then
          call grab_XXvelS(Lrange,L_trans,InfoGlob)
          if (flag_dissipation) call grab_Xhistories(Lrange,L_trans,InfoGlob)
       endif
    endif
!!$
!!$
!!$
!!$
    ! (3) get K matrix (and also get phi matrix)
    if (.not. flag_diagonalisation .and. (find_chdens .or. restart_LorK)) then
       call LNV_matrix_multiply(electrons, energy_tmp, doK, dontM1,   &
            dontM2, dontM3, dontM4, dophi, dontE, &
            mat_phi=matphi)
       electrons_tot = spin_factor * sum(electrons)
       if (inode == ionode .and. iprint_init > 1)              &
            write (io_lun,*) 'Got elect: (Nup, Ndn, Ntotal) ', &
            electrons(1), electrons(nspin),   &
            electrons_tot
    end if
!!$
!!$
!!$
!!$
    ! (4) get the core correction to the pseudopotential energy
    !     (note that this routine does not need to be passed a pig)
!!$
!!$
!!$
!!$
    ! (5) Find the Ewald energy for the initial set of atoms
    if (inode == ionode .and. iprint_init > 1) &
         write (io_lun, *) 'Ionic electrostatics'
    if(flag_neutral_atom) then
       if (inode == ionode .and. iprint_init > 1) &
            write (io_lun, *) 'Calling screened_ion_interaction'
       call screened_ion_interaction
    else
       if (inode == ionode .and. iprint_init > 1) &
            write (io_lun, *) 'Calling ewald'
       call ewald
    end if
!!$
!!$
!!$
!!$
    ! +++
    ! (6) Find the dispersion energy for the initial set of atoms
    if (flag_dft_d2) then
       if ((inode == ionode) .and. (iprint_init > 1) ) &
            write (io_lun, *) 'Calling DFT-D2'
       call dispersion_D2
    end if
    call my_barrier
!!$
!!$
!!$
!!$
    if (inode == ionode .and. iprint_init > 2) &
         write (io_lun, *) 'Find_chdens is ', find_chdens
!!$
!!$
!!$
!!$
    ! (7) Make a self-consistent H matrix and potential
    if (find_chdens) then
       call get_electronic_density(density, electrons, atomfns,     &
            H_on_atomfns(1), inode, ionode,  &
            maxngrid)
       electrons_tot = spin_factor * sum(electrons)
       density = density * ne_in_cell/electrons_tot
       if (inode == ionode .and. iprint_init > 1) &
            write (io_lun, *) 'In initial_H, electrons: ', electrons_tot
       ! if flag_LFD=T, update SF-PAO coefficients with the obtained density unless they have been read
       ! and update S with the coefficients
       if ((.NOT.read_option).AND.flag_LFD) then
          call initial_SFcoeff(.false., .true., fixed_potential, .false.)
          call get_S_matrix(inode, ionode, build_AtomF_matrix=.false.)
       endif
    else if (restart_rho .and. atomf==sf) then
       ! when atomf/=sf, density was already grabbed in (0).
       if (nspin == 2) then
          call grab_charge(density(:,1), n_my_grid_points, inode, spin=1)
          call grab_charge(density(:,2), n_my_grid_points, inode, spin=2)
       else
          ! grab_charge without spin optional variable should fetch
          ! chden.* which should contain the total charge
          ! density. However density(:,1) should only contain the spin
          ! up component of the charge density. Therefore one need to
          ! divide by the factor of spin_factor
          call grab_charge(density(:,1), n_my_grid_points, inode)
          density = density / spin_factor
       end if
    end if
!!$
!!$
!!$  S C F
!!$
!!$
    if ( flag_self_consistent ) then ! Vary only DM and charge density
       !
       if ( restart_LorK ) then
          record  = .true.
          reset_L = .false.
          call new_SC_potl(record, sc_tolerance, reset_L, &
               fixed_potential, vary_mu, n_L_iterations,  &
               L_tolerance, total_energy, backtrace_level)
          !
       else
          if (flag_LFD .and. .not.read_option) then
             ! Hpao was already made in sub:initial_SFcoeff
             rebuild_KE_NL = .false. 
             call get_H_matrix(rebuild_KE_NL, fixed_potential, electrons, &
                  density, maxngrid, level=backtrace_level, build_AtomF_matrix=.false.)
          else
             rebuild_KE_NL = .true. 
             call get_H_matrix(rebuild_KE_NL, fixed_potential, electrons, &
                  density, maxngrid, level=backtrace_level)
          endif
          !
          electrons_tot = spin_factor * sum(electrons)
          !
          record  = .false.
          reset_L = .true.                
          call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
               reset_L, record, backtrace_level)
          !
          record  = .true.             
          reset_L = .false.
          call new_SC_potl(record, sc_tolerance, reset_L, &
               fixed_potential, vary_mu, n_L_iterations,  &
               L_tolerance, total_energy, backtrace_level)
          !
       end if
       !
    else ! Ab initio TB: vary only DM

       rebuild_KE_NL = .true.
       !build_X = .false
       if (flag_LFD .and. .not.read_option) then
          ! Hpao was already made in sub:initial_SFcoeff
          rebuild_KE_NL = .false.
          call get_H_matrix(rebuild_KE_NL, fixed_potential, electrons, &
               density, maxngrid, level=backtrace_level, build_AtomF_matrix=.false.)
       else
          rebuild_KE_NL = .true.
          call get_H_matrix(rebuild_KE_NL, fixed_potential, electrons, &
               density, maxngrid, level=backtrace_level)
       endif
       electrons_tot = spin_factor * sum(electrons)
       if (flag_out_wf.OR.flag_write_DOS) then
          wf_self_con=.true.
       endif

       if ( .not. restart_LorK ) then
          record  = .false.   
          reset_L = .true.
          call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
               reset_L, record, backtrace_level)
       else
          record  = .false.
          reset_L = .false.
          call FindMinDM(n_L_iterations, vary_mu, L_tolerance, &
               reset_L, record, backtrace_level)
       end if
       if (flag_out_wf.OR.flag_write_DOS) then
          wf_self_con=.false.
       endif
       call get_energy(total_energy=total_energy,level=backtrace_level)
    end if
!!$
!!$
!!$
!!$
    ! Do we want to just test the forces ?
    if (flag_test_forces) then
       call test_forces(fixed_potential, vary_mu, n_L_iterations, &
            L_tolerance, sc_tolerance, total_energy,  &
            expected_reduction)
       call end_comms
       stop
    end if

    !****lat<$
    call stop_backtrace(t=backtrace_timer,who='initial_H',echo=.true.)
    !****lat>$

    return

5   format(2x,'Energy as 2Tr[K.H]: ',f18.11,' eV')
6   format(2x,'2Tr[S.G]: ',f18.11,' eV')

  end subroutine initial_H
  !!***


  ! ------------------------------------------------------------------------------
  ! setgrid
  ! ------------------------------------------------------------------------------

  !!****f* initialisation/setgrid *
  !!
  !!  NAME
  !!   setgrid_new
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Overall control of everything grid-related
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   T. Miyazaki
  !!  CREATION DATE
  !!   Sometime 2000-2001
  !!  MODIFICATION HISTORY
  !!   09:33, 11/05/2005 dave
  !!    Added ROBODoc header, indented code, added 1.1 factor to rcut_max
  !!   2011/11/16 15:53 dave
  !!    Removed Extent from set_blipgrid call
  !!   2016/09/16 17:00 nakata
  !!    Added atomf and RadiusAtomf
  !!  SOURCE
  !!
  subroutine setgrid(myid, r_core_squared, r_h)

    !Modules and Dummy Arguments
    use datatypes
    use numbers,                only: very_small
    use global_module,          only: x_atom_cell,y_atom_cell,        &
                                      z_atom_cell, ni_in_cell,        &
                                      iprint_index,                   &
                                      atomf, sf
    use block_module,           only: n_pts_in_block
    use maxima_module,          only: maxblocks
    use group_module,           only: blocks, parts
    use construct_module,       only: init_primary, init_cover
    use primary_module,         only: domain, bundle, make_prim
    use cover_module,           only: DCS_parts, BCS_blocks, make_cs, &
                                      make_iprim, send_ncover,        &
                                      BCS_parts
    use set_blipgrid_module,    only: set_blipgrid
    use set_bucket_module,      only: set_bucket
    use GenComms,               only: my_barrier
    use dimens,                 only: RadiusAtomf
    use pseudopotential_common, only: core_radius
    use input_module,           only: leqi

    implicit none
    integer,intent(in)      :: myid
    real(double),intent(in) :: r_core_squared, r_h

    ! Local variables
    ! Temporary
    real(double) :: rcut_max, r_core

    !-- Start of the subroutine (set_grid_new)
    if(myid == 0 .and. iprint_index > 1) &
         write (io_lun, *) 'setgrid_new starts'
    !if(iprint_index > 4) write(io_lun,*) ' setgrid_new starts for myid= ',myid

    !Sets up domain
    call init_primary(domain, n_pts_in_block*maxblocks, maxblocks, .false.)
    !if(iprint_index > 4) write(io_lun,*) 'init_primary end for myid = ',myid,' par = ',n_pts_in_block, maxblocks
    ! call my_barrier()  !TMP
    call make_prim(blocks, domain, myid)
    !if(iprint_index > 4) write(io_lun,*) 'make_prim end for myid = ',myid
    !Sets up DCS_parts& BCS_blocks
    r_core = sqrt(r_core_squared)
    ! No longer necessary as this is done in dimens_module
    !if(leqi(runtype,'static')) then
    rcut_max = max(r_core,r_h) + very_small
    !else
    !   rcut_max=1.1_double*(max(r_core,r_h)+very_small)
    !endif
    !if(iprint_index > 4) write(io_lun,*) ' rcut_max for DCS_parts and BCS_blocks = ',rcut_max, r_core, r_h

    call make_cs(myid,rcut_max, DCS_parts , parts , domain, &
                 ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell)
    !if(iprint_index > 4) write(io_lun,*) 'Node ',myid+1,' Done make_DCSparts'
    call make_cs(myid,rcut_max, BCS_blocks, blocks, bundle)
    !if(iprint_index > 4) write(io_lun,*) 'Node ',myid+1,' Done make_BCSblocks'

    !call make_iprim(DCS_parts,bundle,myid) !primary number for members
    !  write(io_lun,*) 'Node ',myid+1,' Done make_iprim'
    call my_barrier

    call send_ncover(DCS_parts, myid + 1)
    !if(iprint_index > 4) write(io_lun,*) 'Node ',myid+1,' Done send_ncover for DCS_parts'
    call my_barrier
    call send_ncover(BCS_blocks, myid + 1)
    !if(iprint_index > 4) write(io_lun,*) 'Node ',myid+1,' Done send_ncover for BCS_blocks'

    ! Commented out DRB 2008/09/10: if we really need this, then uncomment and recompile
    !if(iprint_index > 4) call check_setgrid

    !if(iprint_index > 4) write(io_lun,*) ' DCS & BCS has been prepared for myid = ',myid
    !Makes variables used in Blip-Grid transforms
    ! See (naba_blk_module.f90), (set_blipgrid_module.f90), (make_table.f90)
    call set_blipgrid(myid, RadiusAtomf, core_radius)
    !if(iprint_index > 4) write(io_lun,*) 'Node ',myid+1,' Done set_blipgrid'

    !Makes variables used in calculation (integration) of matrix elements
    ! See (bucket_module.f90) and (set_bucket_module.f90)
    call set_bucket(myid)
    !if(iprint_index > 4) write(io_lun,*) 'Node ',myid+1,' Done set_bucket'

    return
  end subroutine setgrid
  !!***


  subroutine check_setgrid

    use global_module,  only: numprocs,rcellx,rcelly,rcellz,  &
                              iprint_index
    use GenComms,       only: myid, my_barrier
    use group_module,   only: blocks, parts
    use primary_module, only: domain, bundle, make_prim
    use cover_module,   only: DCS_parts, BCS_blocks, make_cs, &
                              make_iprim, send_ncover, BCS_parts

    implicit none

    integer      :: nnd, nnd2
    real(double) :: xx, yy, zz, dcellx, dcelly, dcellz
    integer      :: nnx, nny, nnz
    integer      :: no_of_naba_atom, ip
    integer      :: ierror = 0

    !FOR DEBUGGING
    if(iprint_index > 6) then
       call my_barrier()
       !-- CHECK -- bundle
       dcellx=rcellx/parts%ngcellx
       dcelly=rcelly/parts%ngcelly
       dcellz=rcellz/parts%ngcellz
       do nnd=1,numprocs
          if(myid == nnd-1) then
             write(io_lun,*)
             write(io_lun,*) ' Node ',nnd,' CHECK bundle '
             write(io_lun,*) ' n_prim, groups_on_node = ',&
                             bundle%n_prim,bundle%groups_on_node
             write(io_lun,101) bundle%nx_origin,bundle%ny_origin,bundle%nz_origin
             write(io_lun,102) bundle%nw_primx ,bundle%nw_primy ,bundle%nw_primz
             write(io_lun,103) bundle%nleftx   ,bundle%nlefty   ,bundle%nleftz
             nnx=bundle%nx_origin-bundle%nleftx
             nny=bundle%ny_origin-bundle%nlefty
             nnz=bundle%nz_origin-bundle%nleftz
             xx=(nnx-1)*dcellx
             yy=(nny-1)*dcelly
             zz=(nnz-1)*dcellz
             write(io_lun,104) xx,yy,zz
             nnx=nnx+bundle%nw_primx-1
             nny=nny+bundle%nw_primy-1
             nnz=nnz+bundle%nw_primz-1
             xx=nnx*dcellx
             yy=nny*dcelly
             zz=nnz*dcellz
             write(io_lun,105) xx,yy,zz
101          format(3x,' origin    ',3i5)
102          format(3x,' width     ',3i5)
103          format(3x,' left_span ',3i5)
104          format(3x,' Left  Down Bottom',3f15.6)
105          format(3x,' Right  Up   TOP  ',3f15.6)
          end if
          call my_barrier()
       end do

       !-- CHECK -- domain
       dcellx = rcellx / blocks%ngcellx
       dcelly = rcelly / blocks%ngcelly
       dcellz = rcellz / blocks%ngcellz
       do nnd = 1, numprocs
          if(myid == nnd-1) then
             write(io_lun,*)
             write(io_lun,*) ' Node ',nnd,' CHECK domain '
             write(io_lun,*) ' n_prim, groups_on_node = ', &
                             domain%n_prim,domain%groups_on_node
             write(io_lun,101) domain%nx_origin,domain%ny_origin,domain%nz_origin
             write(io_lun,102) domain%nw_primx ,domain%nw_primy ,domain%nw_primz
             write(io_lun,103) domain%nleftx   ,domain%nlefty   ,domain%nleftz
             nnx=domain%nx_origin-domain%nleftx
             nny=domain%ny_origin-domain%nlefty
             nnz=domain%nz_origin-domain%nleftz
             xx=(nnx-1)*dcellx
             yy=(nny-1)*dcelly
             zz=(nnz-1)*dcellz
             write(io_lun,104) xx,yy,zz
             nnx=nnx+domain%nw_primx-1
             nny=nny+domain%nw_primy-1
             nnz=nnz+domain%nw_primz-1
             xx=nnx*dcellx
             yy=nny*dcelly
             zz=nnz*dcellz
             write(io_lun,105) xx,yy,zz
             !101 format(3x,' origin    ',3i5)
             !102 format(3x,' width     ',3i5)
             !103 format(3x,' left_span ',3i5)
             !104 format(3x,' Left  Down Bottom',3f15.6)
             !105 format(3x,' Right  Up   TOP  ',3f15.6)
             !write(io_lun,*) ' @@ mx_nbonn for node ',nnd,' >= ',domain%groups_on_node
             !if(domain%groups_on_node > mx_nbonn) ierror=ierror+1
          end if
          call my_barrier()
       end do

       !-- CHECK -- BCS_parts
       dcellx = rcellx / parts%ngcellx
       dcelly = rcelly / parts%ngcelly
       dcellz = rcellz / parts%ngcellz
       do nnd2 = 1, numprocs
          if(myid == nnd2-1 ) then
             write(io_lun,*)
             write(io_lun,*) ' Node ', nnd2, ' CHECK BCS_parts'
             write(io_lun,101) BCS_parts%nx_origin, BCS_parts%ny_origin,&
                               BCS_parts%nz_origin
             write(io_lun,102) BCS_parts%ncoverx  , BCS_parts%ncovery  ,&
                               BCS_parts%ncoverz
             write(io_lun,103) BCS_parts%nspanlx  , BCS_parts%nspanly  ,&
                               BCS_parts%nspanlz
             nnx = BCS_parts%nx_origin - BCS_parts%nspanlx
             nny = BCS_parts%ny_origin - BCS_parts%nspanly
             nnz = BCS_parts%nz_origin - BCS_parts%nspanlz
             xx = (nnx-1) * dcellx
             yy = (nny-1) * dcelly
             zz = (nnz-1) * dcellz
             write(io_lun,104) xx, yy, zz
             nnx = nnx + BCS_parts%ncoverx - 1
             nny = nny + BCS_parts%ncovery - 1
             nnz = nnz + BCS_parts%ncoverz - 1
             xx = nnx * dcellx
             yy = nny * dcelly
             zz = nnz * dcellz
             write(io_lun,105) xx,yy,zz
             do nnd = 1, numprocs
                write(io_lun,*) ' for Node = ',nnd,' BCSparts_ncover_remote ',&
                                BCS_parts%ncover_rem(1+3*(nnd-1)),&
                                BCS_parts%ncover_rem(2+3*(nnd-1)),&
                                BCS_parts%ncover_rem(3+3*(nnd-1))
             end do
          end if
          call my_barrier()
       end do
       !-- CHECK -- DCS_parts
       do nnd2 = 1, numprocs
          if(myid == nnd2-1) then
             write(io_lun,*)
             write(io_lun,*) ' Node ', nnd2, ' CHECK DCS_parts'
             write(io_lun,101) DCS_parts%nx_origin, DCS_parts%ny_origin, &
                               DCS_parts%nz_origin
             write(io_lun,102) DCS_parts%ncoverx  , DCS_parts%ncovery  , &
                               DCS_parts%ncoverz
             write(io_lun,103) DCS_parts%nspanlx  , DCS_parts%nspanly  , &
                               DCS_parts%nspanlz
             nnx = DCS_parts%nx_origin - DCS_parts%nspanlx
             nny = DCS_parts%ny_origin - DCS_parts%nspanly
             nnz = DCS_parts%nz_origin - DCS_parts%nspanlz
             xx = (nnx-1) * dcellx
             yy = (nny-1) * dcelly
             zz = (nnz-1) * dcellz
             write(io_lun,104) xx, yy, zz
             nnx = nnx + DCS_parts%ncoverx - 1
             nny = nny + DCS_parts%ncovery - 1
             nnz = nnz + DCS_parts%ncoverz - 1
             xx = nnx * dcellx
             yy = nny * dcelly
             zz = nnz * dcellz
             write(io_lun,105) xx, yy, zz
             do nnd = 1, numprocs
                write(io_lun,*) ' for Node = ',nnd,' BCSparts_ncover_remote ',&
                                DCS_parts%ncover_rem(1+3*(nnd-1)),&
                                DCS_parts%ncover_rem(2+3*(nnd-1)),&
                                DCS_parts%ncover_rem(3+3*(nnd-1))
             end do
             !write(io_lun,*) ' @@ mx_pcover_DCS for node ',nnd2,' >= ',DCS_parts%ng_cover
             !if(DCS_parts%ng_cover > mx_pcover_DCS) ierror=ierror+2
             no_of_naba_atom=0
             do ip=1,DCS_parts%ng_cover
                no_of_naba_atom=no_of_naba_atom+DCS_parts%n_ing_cover(ip)
             end do
             !write(io_lun,*) ' @@ mx_icover_DCS for node ',nnd2,' >= ',&
             !     DCS_parts%icover_ibeg(DCS_parts%ng_cover)+DCS_parts%n_ing_cover(DCS_parts%ng_cover)-1, &
             !     no_of_naba_atom
             !if(mx_icover_DCS < no_of_naba_atom) ierror=ierror+4
          end if
          call my_barrier()
       end do
       !-- CHECK -- BCS_blocks
       dcellx=rcellx/blocks%ngcellx
       dcelly=rcelly/blocks%ngcelly
       dcellz=rcellz/blocks%ngcellz
       do nnd2=1,numprocs
          if(myid ==nnd2-1 ) then
             write(io_lun,*)
             write(io_lun,*) ' Node ', nnd2, ' CHECK BCS_blocks'
             write(io_lun,101) BCS_blocks%nx_origin, BCS_blocks%ny_origin, &
                               BCS_blocks%nz_origin
             write(io_lun,102) BCS_blocks%ncoverx  , BCS_blocks%ncovery  , &
                               BCS_blocks%ncoverz
             write(io_lun,103) BCS_blocks%nspanlx  , BCS_blocks%nspanly  , &
                               BCS_blocks%nspanlz
             nnx=BCS_blocks%nx_origin-BCS_blocks%nspanlx
             nny=BCS_blocks%ny_origin-BCS_blocks%nspanly
             nnz=BCS_blocks%nz_origin-BCS_blocks%nspanlz
             xx=(nnx-1)*dcellx
             yy=(nny-1)*dcelly
             zz=(nnz-1)*dcellz
             write(io_lun,104) xx,yy,zz
             nnx=nnx+BCS_blocks%ncoverx-1
             nny=nny+BCS_blocks%ncovery-1
             nnz=nnz+BCS_blocks%ncoverz-1
             xx=nnx*dcellx
             yy=nny*dcelly
             zz=nnz*dcellz
             write(io_lun,105) xx,yy,zz
             do nnd=1,numprocs
                write(io_lun,*) ' for Node = ',nnd,' BCSparts_ncover_remote ',&
                     BCS_blocks%ncover_rem(1+3*(nnd-1)),&
                     BCS_blocks%ncover_rem(2+3*(nnd-1)),&
                     BCS_blocks%ncover_rem(3+3*(nnd-1))
             enddo
             write(io_lun,*) ' @@ mx_bcover for node ',nnd2,&
                             ' >= ',BCS_blocks%ng_cover
             !if(BCS_blocks%ng_cover > mx_bcover) ierror=ierror+8
          endif
          call my_barrier()
       enddo
    end if
    !END OF DEBUGGING

  end subroutine check_setgrid


  !!****f* initialisation/initial_L *
  !!
  !!  NAME
  !!   initial_L
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Finds initial L (set equal to 1/2 S^-1)
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler/C.M.Goringe
  !!  CREATION DATE
  !!   07/03/95
  !!  MODIFICATION HISTORY
  !!   04/05/01 dave
  !!    Takes S^-1 from Hotelling's method
  !!   21/06/2001 dave
  !!    Added ROBODoc header and indented
  !!   12:20, 2004/06/09 dave
  !!    Fixed bug: Srange not Trange in final option
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new
  !!    matrix routines
  !!   2006/11/14 07:58 dave
  !!    Included in initialisation
  !!   2011/07/01 L.Tong
  !!    Added initialisation for matL_dn, for spin polarisation
  !!   2012/03/27 L.Tong
  !!   - Changed spin implementation
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2018/11/13 17:30 nakata
  !!    Changed matT to be spin_SF dependent
  !!   2018/11/15 15:45 nakata
  !!    Bug fix: matL(1) should be matL(spin)
  !!  SOURCE
  !!
  subroutine initial_L(level)

    use datatypes
    use numbers,       only: half, zero
    use mult_module,   only: matL, matT, matrix_sum
    use global_module, only: nspin, flag_SpinDependentSF

    implicit none

    integer, optional :: level
    integer           :: spin, spin_SF
    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level 

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='initial_L',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    spin_SF = 1
    do spin = 1, nspin
       if (flag_SpinDependentSF) spin_SF = spin
       ! set L for the second spin component also equal to 1/2 S^-1
       call matrix_sum(zero, matL(spin), half, matT(spin_SF))
    end do

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='initial_L',echo=.true.)
!****lat>$

    return
  end subroutine initial_L
  !!***


end module initialisation
