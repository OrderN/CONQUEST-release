! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: datatypes.module.f90,v 1.2 2002/04/19 13:56:24 drb Exp $
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
!!  SOURCE
!!
module datatypes

  ! RCS tag for object file identification 
  character(len=80), save, private :: RCSid = "$Id: datatypes.module.f90,v 1.2 2002/04/19 13:56:24 drb Exp $"

  ! This picks a the most memory efficient integer kind
  ! that can hold number from -10**6 to 10**6. This
  ! will typically be a 32 bit integer.
  integer, parameter :: short = selected_int_kind( 6 )

  ! The default types.
  integer, parameter :: double = selected_real_kind( 6, 70 )
  integer, parameter :: double_cplx = selected_real_kind( 6, 70 )
  integer, parameter :: integ  = selected_int_kind( 9 )

end module datatypes
!!***

      

