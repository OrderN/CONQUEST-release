! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module potential_module
! ------------------------------------------------------------------------------
! Code area 3: operators
! ------------------------------------------------------------------------------

!!****h* Conquest/potential_module *
!!  NAME
!!   potential_module
!!  PURPOSE
!!   Hold routines and data associated with local potential (this means the total potential
!!   on the grid at the moment, but may grow to include neutral atom potential in the 
!!   future).  Will hold local potential expansion ideas of Ozaki in the near future.
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006/08/02 17:28 dave
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!
module potential_module

  use datatypes

  implicit none
  save

  ! RCS tag for object file identification
  character(len=80), private :: RCSid = "$Id$"

  real(double), allocatable, dimension(:) :: potential

!!***

end module potential_module
