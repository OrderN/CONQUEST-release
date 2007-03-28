! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_hartree
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_hartree *
!!
!!NAME
!! cq_ut_psp2pao_hartree
!!PURPOSE
!! Hartree potential and energy related tasks
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_hartree

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine get_v_hartree
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_hartree/get_v_hartree *
!!
!!NAME
!! get_v_hartree
!!USAGE
!!
!!PURPOSE
!! Build the hartree potential
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_numeric
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 03/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine get_v_hartree

     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_numeric

     implicit none

     ! Local variables
     integer :: stat, i
     real(double) :: rho, result

     ! Parameters
     real(double), parameter :: four_pi = 4.0_double * 3.1415926535897932_double

     do i=1,gl_points_mesh
       rho = gl_rho(i) + gl_rho(i + gl_points_mesh)

       gl_tmp1(i) = four_pi * gl_r(i) * rho
       gl_tmp2(i) = gl_tmp1(i) * gl_r(i)
     end do

     ! Calculate derivatives
     call finite_differentiation (gl_tmp1, gl_dtmp1)
     call finite_differentiation (gl_tmp2, gl_dtmp2)

     gl_tmp3(1) = 0.0_double
     gl_dtmp3(1) = 0.0_double

     do i=2,gl_points_mesh
        call integral_cubic(gl_r, gl_tmp2, gl_dtmp2, i - 1, i, result)
        gl_tmp3(i) = gl_tmp3(i-1) + result

        call integral_cubic(gl_r, gl_tmp1, gl_dtmp1, i - 1, i, result)
        gl_dtmp3(i) = gl_dtmp3(i-1) + result
     end do

     do i=1,gl_points_mesh
        gl_v_hartree(i) = gl_dtmp3(gl_points_mesh) - gl_dtmp3(i) &
                        + gl_tmp3(i) / gl_r(i)
     end do

  end subroutine get_v_hartree
!!***

end module cq_ut_psp2pao_hartree

