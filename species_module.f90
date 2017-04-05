! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module species_module
! ------------------------------------------------------------------------------
! Code area 9: General
! ------------------------------------------------------------------------------

!!****h* Conquest/species_module
!!  NAME
!!   species_module
!!  PURPOSE
!!   Holds variables relating to different species
!!  USES
!!   common, datatypes
!!  AUTHOR
!!   E.H.Hernandez/C.M.Goringe/I.J.Bush
!!  CREATION DATE
!!   1997 or 1998 sometime
!!  MODIFICATION HISTORY
!!   31/05/2001 dave
!!    Added ROBODoc header and RCS Id tag to top line
!!   18/03/2002 dave
!!    Added RCS log and static tag for object file id
!!   10:05, 04/02/2003 drb 
!!    Changed length of ps_file to 50
!!   08:20, 2003/06/11 dave
!!    Added nlpf_species
!!   12:00, 12/11/2004 dave 
!!    Removed common use and replaced with global for mx_icell
!!   2006/03/07 23:52 dave
!!    Added nsf_species
!!   2006/07/14 19:09 dave
!!    Made all arrays (apart from species) allocatable
!!   2011/09/16 10:58 dave
!!    Changed species_label to length 30
!!   2014/10/12 10:58 lat
!!    Added species_file
!!   2016/11/09 21:00 nakata
!!    Added natomf_species
!!  SOURCE
!!
module species_module

  use datatypes

  implicit none

  save

  ! RCS tag for object file identification 
  character(len=80), private :: RCSid = "$Id$"

  integer :: n_species
  logical :: species_from_files
  !
  integer,           allocatable, dimension(:) :: species
  real(double),      allocatable, dimension(:) :: charge, mass
  character(len=50), allocatable, dimension(:) :: ps_file
  character(len=40), allocatable, dimension(:) :: ch_file
  character(len=40), allocatable, dimension(:) :: phi_file
  !
  character(len=30), allocatable, dimension(:) :: species_label
  character(len=80), allocatable, dimension(:) :: species_file
  !
  logical, allocatable, dimension(:)  :: non_local_species
  integer, allocatable, dimension(:)  :: nlpf_species
  integer, allocatable, dimension(:)  :: nsf_species
  integer, allocatable, dimension(:)  :: npao_species
  integer, allocatable, dimension(:)  :: natomf_species
  integer, allocatable, dimension(:)  :: type_species

!!***

end module species_module
