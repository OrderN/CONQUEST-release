! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module maxima_module
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/maxima_module
!!  NAME
!!   maxima_module
!!  PURPOSE
!!   Contains maxima for matrix arrays in particular and most other arrays
!!   along with routines to assign appropriate maxima to derived types
!!  USES
!!   matrix_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    ROBODoc header, removed group_, prim_ and cover_maxima
!!   18/03/2002 dave
!!    Added RCS Id and Log tags, and static tag for object file id
!!   12:40, 25/09/2002 mjg & drb 
!!    Added maxima for neighbour tables of atomic charge densities
!!   2018/01/19 14:11 JST dave
!!    Added maxima for PAO and PS l values
!!  SOURCE
!!
module maxima_module

  implicit none

  ! RCS tag for object file identification 
  character(len=80), save, private :: RCSid = "$Id$"

  !include 'maxima.inc'

  integer :: maxnsf, maxngrid, maxblocks
  integer :: maxpartsproc, maxatomspart, maxatomsproc, maxpartscell
  integer :: maxnabaprocs
  integer :: lmax_pao, lmax_ps
  
  ! Temporary !
  integer,parameter :: max_blip_nu_int       =        5

!!***

contains

! -----------------------------------------------------------
! Subroutine matrix_maxima
! -----------------------------------------------------------

!!****f* maxima_module/matrix_maxima *
!!
!!  NAME 
!!   matrix_maxima
!!  USAGE
!! 
!!  PURPOSE
!!   Assigns maxima to matrix derived type
!!  INPUTS
!! 
!! 
!!  USES
!!   matrix_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    ROBODoc header
!!  SOURCE
!!
  subroutine matrix_maxima(mat,mx_nab,mx_abs)

    ! Module usage
    use matrix_module

    implicit none

    ! Passed variables
    type(matrix), dimension(:) :: mat
    integer :: mx_nab,mx_abs

    mat%mx_nab = mx_nab  ! Average neighbours of atom
    mat%mx_abs = mx_abs  ! Maximum neighbours of atom
    return
  end subroutine matrix_maxima
!!***

! -----------------------------------------------------------
! Subroutine ltrans_maxima
! -----------------------------------------------------------

!!****f* maxima_module/ltrans_maxima *
!!
!!  NAME 
!!   ltrans_maxima
!!  USAGE
!! 
!!  PURPOSE
!!   Assigns maxima to local transpose derived type
!!  INPUTS
!! 
!! 
!!  USES
!!   matrix_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    ROBODoc header
!!  SOURCE
!!
  subroutine ltrans_maxima(tr,mx_halo,mx_nab)

    ! Module usage
    use matrix_module

    implicit none

    ! Passed variables
    type(matrix_trans) :: tr
    integer :: mx_halo,mx_nab

    tr%mx_halo = mx_halo  ! Maximum atoms in halo
    tr%mx_nab = mx_nab    ! Average neighbours of atom
    return
  end subroutine ltrans_maxima
!!***

! -----------------------------------------------------------
! Subroutine halo_maxima
! -----------------------------------------------------------

!!****f* maxima_module/halo_maxima *
!!
!!  NAME 
!!   halo_maxima
!!  USAGE
!! 
!!  PURPOSE
!!   Assigns maxima to halo derived type
!!  INPUTS
!! 
!! 
!!  USES
!!   matrix_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    ROBODoc header
!!  SOURCE
!!
  subroutine halo_maxima(halo,mx_part,mx_halo)

    ! Module usage
    use matrix_module

    implicit none

    ! Passed variables
    type(matrix_halo) :: halo
    integer :: mx_part,mx_halo

    halo%mx_part = mx_part  ! Maximum partitions in halo
    halo%mx_halo = mx_halo  ! Maximum atoms in halo
  end subroutine halo_maxima
!!***
end module maxima_module
