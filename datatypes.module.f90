! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module datatypes
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/datatypes
!!  NAME
!!   datatypes -- defines various useful datatypes
!!  PURPOSE
!!   Define standard, portable datatypes of same precision
!!  AUTHOR
!!   I.J.Bush
!!  CREATION DATE
!!   Sometime in 1996 or 1997
!!  MODIFICATION HISTORY
!!   By D.R.Bowler from time to time
!!   30/05/2001 dave
!!    Added ROBODoc header
!!   18/03/2002 dave
!!    Added static id for object files and tidied header
!!   2014/07/24 11:40 dave
!!    Added type wide for 64-bit integer (i.e. that can include beyond 10^9)
!!  SOURCE
!!
module datatypes

  ! RCS tag for object file identification 
  character(len=80), save, private :: RCSid = "$Id$"

  ! This picks a the most memory efficient integer kind
  ! that can hold number from -10**6 to 10**6. This
  ! will typically be a 32 bit integer.
  integer, parameter :: short = selected_int_kind( 6 )

  ! The default types.
  integer, parameter :: double = selected_real_kind( 6, 70 )
  integer, parameter :: double_cplx = selected_real_kind( 6, 70 )
  integer, parameter :: integ  = selected_int_kind( 9 )
  integer, parameter :: wide = selected_int_kind(15) 

end module datatypes
!!***

      

