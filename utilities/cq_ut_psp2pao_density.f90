! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_density
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_density *
!!
!!NAME
!! cq_ut_psp2pao_density
!!PURPOSE
!! Build the density
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 15/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_density

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine get_rho
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_density/get_rho *
!!
!!NAME
!! get_rho
!!USAGE
!! 
!!PURPOSE
!! Build the density
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 15/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine get_rho

     use cq_ut_psp2pao_global

     implicit none

     integer :: i, j
     real(double) :: psi

     do i=1, 2*gl_points_mesh
        gl_rho(i) = 0.0_double
     end do

     do i=0, gl_no_orbitals-1
        do j=1, gl_points_mesh
           psi = gl_ul(i * gl_points_mesh + j) / gl_r(j)
           gl_rho(j) = gl_rho(j) + gl_occ(i+1) * psi *psi

           psi = gl_ul(gl_no_orbitals * gl_points_mesh &
                     + i * gl_points_mesh + j) / gl_r(j)
           gl_rho(gl_points_mesh + j) = gl_rho(gl_points_mesh + j) &
                                      + gl_occ(gl_no_orbitals + i + 1) * psi *psi
        end do
     end do

  end subroutine get_rho
!!***

! -----------------------------------------------------------------------------
! Subroutine get_rho_basis
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_density/get_rho_basis *
!!
!!NAME
!! get_rho_basis
!!USAGE
!!
!!PURPOSE
!! Build the density for the basis
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 17/05/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine get_rho_basis

     use cq_ut_psp2pao_global

     implicit none

     integer :: i, j
     real(double) :: psi

     do i=1, gl_basis%points_mesh
        gl_basis%rho(i) = 0.0_double
     end do

     do i=1, gl_basis%no_orbitals
        do j=1, gl_basis%points_mesh
           psi = gl_basis%orb_ul(i,j)
           gl_basis%rho(j) = gl_basis%rho(j) + gl_basis%orb_occ(i) * psi *psi
           gl_basis%rho(j+gl_basis%points_mesh) = gl_basis%rho(j)
        end do
     end do

  end subroutine get_rho_basis
!!***

end module cq_ut_psp2pao_density
