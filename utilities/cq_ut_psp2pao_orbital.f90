! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_orbital
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_orbital *
!!
!!NAME
!! cq_ut_psp2pao_orbital
!!PURPOSE
!! Treatment of orbitals
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_orbital

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine orthogonalise
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_orbital/orthogonalise *
!!
!!NAME
!! orthogonalise
!!USAGE
!!
!!PURPOSE
!! Orthogonalise a basis
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine orthogonalise

     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_numeric

     implicit none

     integer :: i, j, k, stat
     real(double), dimension(:), allocatable :: tmp
     real(double) :: product

     ! Parameters
     real(double), parameter :: four_pi = 4.0_double * 3.1415926535897932_double

     allocate(tmp(gl_basis%points_mesh), stat=STAT)

     ! Normalisation
     do i=1, gl_basis%no_orbitals
        do j=1, gl_basis%points_mesh
           tmp(j) = four_pi * gl_basis%r(j) * gl_basis%r(j) &
                  * gl_basis%orb_ul(i,j) * gl_basis%orb_ul(i,j)
        end do
        call simpson_integral_prod(gl_basis%r, tmp, gl_basis%points_mesh, &
                                   gl_basis%dnu, product)
        product=1.0_double / sqrt(product)
        do j=1, gl_basis%points_mesh
           gl_basis%orb_ul(i,j)=gl_basis%orb_ul(i,j) * product
        end do
     end do

     ! Orthogonalisation
     do i=2, gl_basis%no_orbitals
        do j=1, i-1
           if (gl_basis%orb_l(i) == gl_basis%orb_l(j)) then
              do k=1, gl_basis%points_mesh
                 tmp(k) = four_pi * gl_basis%r(k) * gl_basis%r(k) &
                        * gl_basis%orb_ul(i,k) * gl_basis%orb_ul(j,k)
              end do
              call simpson_integral_prod(gl_basis%r, tmp, &
                                         gl_basis%points_mesh, &
                                         gl_basis%dnu, product)
              do k=1, gl_basis%points_mesh
                 gl_basis%orb_ul(i,k) = gl_basis%orb_ul(i,k) &
                                      - gl_basis%orb_ul(j,k) * product
              end do

              ! Normalisation
              do k=1, gl_basis%points_mesh
                 tmp(k) = four_pi * gl_basis%r(k) * gl_basis%r(k) &
                        * gl_basis%orb_ul(i,k) * gl_basis%orb_ul(i,k)
              end do
              call simpson_integral_prod(gl_basis%r, tmp, &
                                         gl_basis%points_mesh, &
                                         gl_basis%dnu, product)
              product = 1.0_double / sqrt(product)
              do k=1, gl_basis%points_mesh
                 gl_basis%orb_ul(i,k) = gl_basis%orb_ul(i,k) * product
              end do
           end if
        end do
     end do


     if(allocated(tmp)) deallocate(tmp)   
 
  end subroutine orthogonalise
!!***

! -----------------------------------------------------------------------------
! Subroutine smooth_basis
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_orbital/smooth_basis *
!!
!!NAME
!! smooth_basis
!!USAGE
!!
!!PURPOSE
!! Smooth the tail of the basis orbitals
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 19/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine smooth_basis(range)

     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_numeric

     implicit none

     ! Input parameters

     real(double) :: range

     ! Local variables

     integer :: i, j, stat
     real(double), dimension(:), allocatable :: tmp
     real(double) :: factor, product

     ! Parameters
     real(double), parameter :: four_pi = 4.0_double * 3.1415926535897932_double


     allocate(tmp(gl_basis%points_mesh), STAT=stat)

     do i=1, gl_basis%no_orbitals
        do j=1, gl_basis%points_mesh
           factor = (gl_basis%r(j) - gl_basis%orb_cutoff_radius(i)) / range
           gl_basis%orb_ul(i,j) = gl_basis%orb_ul(i, j) &
                   * (1.0_double - exp(-0.5_double * factor * factor))
        end do

        ! Normalisation
        do j=1, gl_basis%points_mesh
           tmp(j) = four_pi * gl_basis%r(i) * gl_basis%r(i) &
                  * gl_basis%orb_ul(i,j) * gl_basis%orb_ul(i,j)
        end do
        call simpson_integral_prod(gl_basis%r, tmp, gl_basis%points_mesh, &
                                   gl_basis%dnu, product)
        product = 1 / sqrt(product)

        do j=1, gl_basis%points_mesh
           gl_basis%orb_ul(i,j) = gl_basis%orb_ul(i,j) * product
        end do
     end do

     if(allocated(tmp)) deallocate(tmp)
 
  end subroutine smooth_basis
!!***

end module cq_ut_psp2pao_orbital
