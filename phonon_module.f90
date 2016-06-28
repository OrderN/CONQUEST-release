! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module phonon_module
! ------------------------------------------------------------------------------
! Code area XX: spectroscopy
! ------------------------------------------------------------------------------

!!***h* Conquest/phonon_module *
!!  NAME
!!   phonon_module
!!  CREATION DATE
!!   2016/28/06 Z. Raza and L. A. Truflandier
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module phonon_module

  use datatypes
  use GenComms,  only: inode, ionode
  use phonon_io, only: phonon_global_write
  
  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"

!!***

contains

  !!***
  
  !!****f* cdft_module/dummy *
  !!
  !!  NAME 
  !!
  !!  USAGE
  !!
  !!  PURPOSE
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   
  !!  CREATION DATE
  !!
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!  
  subroutine dummy_phonon()

    
    implicit none

    call phonon_global_write()
    
    return
  end subroutine dummy_phonon
  !!***

end module phonon_module


