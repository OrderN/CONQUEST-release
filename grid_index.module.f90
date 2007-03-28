! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: grid_index.module.f90,v 1.3 2004/11/12 02:38:31 drb Exp $
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
!!  SOURCE
!!
module grid_index

  implicit none
  save

  ! RCS tag for object file identification 
  character(len=80), private :: RCSid = "$Id: grid_index.module.f90,v 1.3 2004/11/12 02:38:31 drb Exp $"

  integer, allocatable, dimension(:) :: grid_point_x, grid_point_y, grid_point_z, grid_point_block, &
          grid_point_position

  integer, allocatable, dimension(:) :: ind_block_x, ind_block_y, ind_block_z !(mx_nbonn)

end module grid_index
!!***
