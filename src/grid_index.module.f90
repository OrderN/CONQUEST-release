! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module grid_index
! ------------------------------------------------------------------------------
! Code area 8: indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/grid_index
!!  NAME
!!   grid_index
!!  PURPOSE
!!   Collects together variables used in indexing the grid
!!  USE
!!   common
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   06/04/01
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    ROBODoc header
!!   18/03/2002 dave
!!    Added RCS Id and Log tags and static tag for object file id
!!   11:38, 12/11/2004 dave 
!!    Pared down to eight arrays; these will go when ffts are rewritten
!!   2006/11/02 10:25 dave
!!    Added setgrid as subroutine
!!   2022/05/18 14:04 dave
!!    Made grid_point_x/y/z targets for pointers in surface dipole correction
!!  SOURCE
!!
module grid_index

  implicit none
  save

  integer, allocatable, dimension(:), target :: grid_point_x, grid_point_y, grid_point_z
  integer, allocatable, dimension(:) :: grid_point_block, grid_point_position

  integer, allocatable, dimension(:) :: ind_block_x, ind_block_y, ind_block_z !(mx_nbonn)

end module grid_index
!!***
