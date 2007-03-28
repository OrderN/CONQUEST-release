! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module 
! ------------------------------------------------------------------------------
! Code area 8: General
! ------------------------------------------------------------------------------

!!****h* Conquest/units
!!  NAME
!!   units
!!  PURPOSE
!!   Holds units for converting energy, force etc
!!  AUTHOR
!!   D. R. Bowler
!!  CREATION DATE
!!   2006/11/03
!!  MODIFICATION HISTORY
!! 
!!  SOURCE
!!
module units

  use datatypes
  use numbers

  implicit none
  save

  real(double), parameter :: HaToRy = two
  real(double), parameter :: HaToeV = 27.211383_double
  real(double), parameter :: RyToeV = 13.605691_double
  real(double), parameter :: BohrToAng = 0.52917721_double
  real(double), parameter :: AngToBohr = 1.88972613_double
  real(double), parameter :: kB = 0.0009765625_double
  real(double), parameter :: MB = 9.5367431640625e-07_double
  real(double), parameter :: GB = 9.3132257461547852e-10_double
  

  integer, parameter :: har = 1
  integer, parameter :: ryd = 2
  integer, parameter :: ev = 3
  integer, parameter :: bohr = 1
  integer, parameter :: ang = 2
  integer, parameter :: kbytes = 1
  integer, parameter :: mbytes = 2
  integer, parameter :: gbytes = 3

  character(len=2), dimension(3) :: en_units = (/"Ha","Ry","eV"/)
  character(len=2), dimension(2) :: d_units = (/"a0","A "/)
  character(len=2), dimension(3) :: mem_units = (/"kB","MB","GB"/)

  integer :: dist_units, energy_units, m_units
  real(double) :: for_conv, en_conv, dist_conv, mem_conv

end module units
!!***
