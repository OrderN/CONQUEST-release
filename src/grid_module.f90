! $Id$
! -----------------------------------------------------------
! Module grid_module
! -----------------------------------------------------------

!!****h* Conquest/block_module
!!  NAME
!!   grid_module
!!  PURPOSE
!!   This module is preared for non-orthorhombic calculations.
!!   It should include 
!!   (we may merge this module to some other modules in the future.)
!!
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   13/Nov/2025  
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module grid_module

  use datatypes
  use global_module, only: cell_vector, inv_cell_vector,io_lun, shift_in_frac, shift_in_bohr
  use GenComms, only: myid
  !use dimens, only: grid_point_volume
  implicit none
  real(double) :: dcell_block(3)
  real(double) :: dcell_grid(3)
  real(double) :: grid_point_volume

contains
 subroutine set_grid_parameters
  use group_module, only:blocks
  use block_module, only:nx_in_block, ny_in_block, nz_in_block
  use lattice_module, only: cell_length, volume
  implicit none
  integer :: ii, ngridx, ngridy, ngridz
 
  !dcell_block
   dcell_block (1) = cell_length(1)/blocks%ngcellx
   dcell_block (2) = cell_length(2)/blocks%ngcelly
   dcell_block (3) = cell_length(3)/blocks%ngcellz
   ! xblock, yblock, zblock should be defined with frac_block(1:3)

  !dcell_grid 
   dcell_grid (1) = dcell_block(1)/nx_in_block
   dcell_grid (2) = dcell_block(2)/ny_in_block
   dcell_grid (3) = dcell_block(3)/ny_in_block

  !grid_point_volume
    ngridx = blocks%ngcellx * nx_in_block
    ngridy = blocks%ngcelly * ny_in_block
    ngridz = blocks%ngcellz * nz_in_block
    !to 
   grid_point_volume = volume/ngridx/ngridy/ngridz

  return
 end subroutine set_grid_parameters

end module grid_module
