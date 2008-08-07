! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module numbers
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/numbers
!!  NAME
!!   numbers
!!  PURPOSE
!!   Defines various numerical constants
!!  USES
!!   datatypes
!!  AUTHOR
!!   C.M.Goringe/E.H.Hernandez/I.J.Bush/D.R.Bowler
!!  CREATION DATE
!!   Sometime between 1995 and 1999
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    Added ROBODoc header
!!   31/05/2001 dave
!!    Added minus_i and RCS Id tag
!!   18/03/2002 dave
!!    Added static RCS id and more to header
!!   04/07/2002 mike
!!    Added fifteen, sixteen, three_quarters and three_halves
!!  SOURCE
!!
module numbers

  use datatypes

  implicit none      

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"



  real(double), parameter :: zero = 0.0_double 
  real(double), parameter :: one = 1.0_double 
  real(double), parameter :: two = 2.0_double 
  real(double), parameter :: three = 3.0_double 
  real(double), parameter :: four = 4.0_double 
  real(double), parameter :: five = 5.0_double 
  real(double), parameter :: six = 6.0_double 
  real(double), parameter :: fifteen = 15.0_double
  real(double), parameter :: sixteen = 16.0_double
  real(double), parameter :: minus_four = -4.0_double 
  real(double), parameter :: quarter = 0.25_double
  real(double), parameter :: half = 0.5_double 
  real(double), parameter :: three_quarters = 0.75_double
  real(double), parameter :: three_halves = 1.5_double
  real(double), parameter :: one_third = 1.0_double/3.0_double
  real(double), parameter :: two_thirds = 2.0_double/3.0_double
  real(double), parameter :: four_thirds = 4.0_double/3.0_double
  real(double), parameter :: one_sixth = 1.0_double/6.0_double
  real(double), parameter :: five_sixths = 5.0_double/6.0_double
  real(double), parameter :: seven_sixths = 7.0_double/6.0_double
  real(double), parameter :: one_ninth = 1.0_double/9.0_double
  real(double), parameter :: seven_thirtysixths = 7.0_double/36.0_double
  real(double), parameter :: pi = 3.14159265358979323846_double 
  real(double), parameter :: twopi = 6.28318530717958647692_double 
  real(double), parameter :: sqrt_pi = 1.77245385090551588192_double
  real(double), parameter :: very_small = 1.0e-8_double
  real(double), parameter :: BIG = 1.0e12_double 

  complex(double_cplx), parameter :: minus_i = (0.0_double_cplx, -1.0_double_cplx)

end module numbers
!!***
