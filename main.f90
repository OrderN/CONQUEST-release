! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Conquest: O(N) DFT
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/Conquest *
!!
!!  NAME 
!!   Conquest
!!  PURPOSE
!!   Main routine for the Conquest code
!!  USES
!!   common, dimens, logicals, numbers, matrix_data,
!!   pseudopotential_data, global_module, mpi, io_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   Way, way back in the deep history of Conquest - 1995 or so
!!  MODIFICATION HISTORY
!!   07/05/01 by DRB
!!    Added ROBODoc header, changed indentation
!!   28/05/2001 dave
!!    Changed ROBODoc header to a module so that it shows up at the top :)
!!    Stripped control call
!!   29/05/2001 dave
!!    Stripped initialise call
!!    Stripped and inserted main.inc
!!   31/05/2001 dave
!!    Added RCS Id tag to header
!!   07/06/2001 dave
!!    Added use fft_module, datatypes and removed includes and mpi
!!    Also changed comms to cq_init and cq_exit
!!   13/06/2001 dave
!!    Changed call to cq_init to use numprocs
!!   15/03/2002 dave
!!    Added RCS id tag for object file identification and tidied 
!!    comment structure a little (added header)
!!   15:20, 25/09/2002 mjg & drb 
!!    Changed so that the variable density is used from density_module
!!   12:43, 04/02/2003 drb 
!!    Changed to use init_comms and end_comms from GenComms
!!   13:50, 10/02/2003 drb 
!!    Changed to use control_run
!!   15:01, 12/03/2003 drb 
!!    Removed reference to workspace
!!   10:18, 2004/01/09 dave
!!    Added use pseudopotential_common
!!   2008/02/06 08:21 dave
!!    Changed for output to file not stdout
!!   2011/12/11 L.Tong
!!    Removed redundant variable number_of_bands, use ne_in_cell
!!    instead
!!   2012/03/27 L.Tong
!!   - Changed spin implementation
!!   - updated some calls that no longer passes mu
!!   2016/01/28 16:45 dave
!!    Updated module name to ion_electrostatic
!!  SOURCE
!!
program Conquest

  use datatypes
  use blip
  use dimens
  use logicals
  use numbers
  use matrix_data
  use group_module
  use primary_module
  use cover_module
  use mult_module
  use pseudopotential_data
  use pseudopotential_common
  use pseudo_tm_info, only: pseudo
  use pao_format,     only: pao
  use initialisation
  use potential_module
  use global_module
  use force_module,   only: tot_force
  use grid_index
  use atoms
  use species_module
  use fft_module
  use GenComms
  use density_module, only: density
  use control,        only: control_run
  use support_spec_format
  use ol_int_datatypes
  use functions_on_grid
  use ion_electrostatic
  use DiagModule
  use ScalapackFormat
  use atomic_density, only: atomic_density_table, rcut_dens, &
                            flag_atomic_density_from_pao
  use block_module,   only: nx_in_block,ny_in_block,nz_in_block, &
                            n_pts_in_block
  use memory_module
  use minimise
  use timer_stdclocks_module
  !use output_module
  !use complete_module


  implicit none

  ! RCS tag for object file identification 
  character(len=80), save :: &
       CQid = "$Id$"

  ! Global variables (for now !)
  logical      :: fixed_potential

  ! Energies and electrons
  real(double) :: total_energy
  
  ! Chemical potential
  real(double) :: mu
  logical      :: vary_mu

  fixed_potential = .false.
  ! identify what node we are on
  call init_comms(myid, numprocs)
  !call init_timing_system(inode)
 
  !call init_reg_mem
 
  !**<lat>** can be used here to open std output
  !          for printing at the very beginning
  !call init_user_output()

  ! Initialise reads in data and finds an initial, self-consistent
  ! Hamiltonian
  call initialise(vary_mu, fixed_potential, mu, total_energy)
  if (inode == ionode .and. iprint_gen > 0) &
       write (io_lun, '(4x,"Finished initialise")')

  ! control does (at the moment) blip minimisation and simple MD
  call control_run(fixed_potential, vary_mu, total_energy)
  if (inode == ionode .and. iprint_gen > 0) &
       write (io_lun, '(4x,"Finished control")')

  call write_mem_use

  call print_time_report
  
  !**<lat>** as initialise routine we may need a complete routine
  !          cheking if everything is fine before ending ; write_mem_use
  !          and print_time_report can be within this routine 
  !call complete()
  !call close_user_output()

  ! Wrap up
  call my_barrier ()
  !call save_state( inode, ionode, mu, expected_reduction, n_run, output_file)
  !call save_blip(inode, ionode)
  call end_comms ()

end program Conquest
!!*****
