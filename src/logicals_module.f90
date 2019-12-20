! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module logicals
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/logicals *
!!  NAME
!!   logicals
!!  PURPOSE
!!   Collects together useful aliases for true and false
!!  AUTHOR
!!   E.H.Hernandez/C.M.Goringe/I.J.Bush
!!  CREATION DATE
!!   Somewhere from 1995 to 1997
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    Added ROBODoc header
!!   18/03/2002 dave
!!    Increased header, added RCS tag for object file id
!!  SOURCE
!!
module logicals

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"

  logical, parameter :: doK = .true. 
  logical, parameter :: dontK = .false. 
  logical, parameter :: doM1 = .true. 
  logical, parameter :: dontM1 = .false. 
  logical, parameter :: doM2 = .true. 
  logical, parameter :: dontM2 = .false. 
  logical, parameter :: doM3 = .true. 
  logical, parameter :: dontM3 = .false. 
  logical, parameter :: doM4 = .true. 
  logical, parameter :: dontM4 = .false. 
  logical, parameter :: dophi = .true. 
  logical, parameter :: dontphi = .false. 
  logical, parameter :: doE = .true. 
  logical, parameter :: dontE = .false. 
  logical, parameter :: doA = .true. 
  logical, parameter :: dontA = .false. 
  logical, parameter :: doT1S = .true. 
  logical, parameter :: dontT1S = .false. 
  logical, parameter :: doGrad = .true. 
  logical, parameter :: dontGrad = .false. 
  logical, parameter :: dox = .true. 
  logical, parameter :: dontx = .false. 
  logical, parameter :: doomega = .true. 
  logical, parameter :: dontomega = .false. 

end module logicals
!!***


