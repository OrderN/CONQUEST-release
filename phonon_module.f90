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

    use mult_module,    only: mat_p, matD, matH
    
    implicit none

    !type matrix_pointer
    !  integer :: length
    !  integer :: sf1_type, sf2_type
    !  real(double), pointer, dimension(:) :: matrix
    !end type matrix_pointer

    print*, 'Dmat'
    print*, mat_p(matD)%length
    print*, mat_p(matD)%sf1_type
    print*, mat_p(matD)%sf2_type
    
    print*, 'Hmat'
    print*, mat_p(matH)%length
    print*, mat_p(matH)%sf1_type
    print*, mat_p(matH)%sf2_type
    
    !call phonon_global_write()
    
    return
  end subroutine dummy_phonon
  !!***

end module phonon_module


