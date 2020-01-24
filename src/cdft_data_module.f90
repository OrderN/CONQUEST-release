! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cdft_data_module
! ------------------------------------------------------------------------------
! Code area 6: energy minimisation
! ------------------------------------------------------------------------------

!!***h* Conquest/cdft_data_module *
!!  NAME
!!   cdft_data_module
!!  CREATION DATE
!!   2009/02/01 A. M. P. Sena and D. R. Bowler
!!  MODIFICATION HISTORY
!!   2011/07/21 16:41 dave
!!    Preparing for inclusion in main trunk
!!   2011/12/10 L.Tong
!!    Added matHzero_dn for spin polarised calculation
!!   2012/03/13 L.Tong
!!    Changed matHzero and matHzero_dn pair to matHzero(:) array,
!!    with index denoting spin channel
!!  SOURCE
!!
module cdft_data

  use datatypes

  implicit none

  integer :: cDFT_Type
  integer, parameter :: cDFT_Fix_Charge = 1
  integer, parameter :: cDFT_Fix_ChargeDifference = 2
  integer :: cDFT_MaxIterations
  integer :: cDFT_NumberAtomGroups
  integer, dimension(:), allocatable :: matHzero
  integer, dimension(:), allocatable :: cDFT_NAtoms
  integer, dimension(:), allocatable :: matWc
  ! Either charge on atoms or charge difference
  real(double), dimension(:), allocatable :: cDFT_Target 
  real(double) :: cDFT_Tolerance
  real(double), dimension(:), allocatable :: cDFT_Vc, cDFT_W
  character(len=10), allocatable, dimension(:) :: cDFT_BlockLabel

  type atom_list
     integer, dimension(:), pointer :: Numbers
  end type atom_list
  type(atom_list), dimension(:), allocatable :: cDFT_AtomList

!!***


end module cdft_data
