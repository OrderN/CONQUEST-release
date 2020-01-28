! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module comm_array_module
! ------------------------------------------------------------------------------
! Code area 8: indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/comm_array_module *
!!  NAME
!!   comm_array_module
!!  PURPOSE
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   2000/1
!!  MODIFICATION HISTORY
!!   10:01, 2003/06/11 dave
!!    Added headers
!!   2008/02/04 17:02 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module comm_array_module

  use datatypes
  use global_module, ONLY: io_lun
  use numbers, ONLY: zero
  use GenComms, ONLY: cq_abort

  implicit none

  save
  real(double), allocatable, dimension(:), target :: send_array
  real(double), allocatable, dimension(:), target :: recv_array
  integer, allocatable, dimension(:), target      :: isend_array
  integer, allocatable, dimension(:), target      :: irecv_array

!!***

end module comm_array_module
