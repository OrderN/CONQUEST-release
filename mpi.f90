! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: mpi.f90,v 1.2 2002/05/07 13:46:54 drb Exp $
! ------------------------------------------------------------------------------
! Module mpi
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/mpi *
!!  NAME
!!   mpi
!!  PURPOSE
!!   Contains variable definitions for mpi calls
!!
!!   The include statement assumes that there is a version of MPI available that
!!   has a freeform fortran include file.  If necessary, this will need to be 
!!   replaced
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   24/04/2002 dave
!!    Added headers and ROBODoc fields
!!  SOURCE
!!
module mpi

  implicit none

  include 'mpif.h'

end module mpi
!!***
