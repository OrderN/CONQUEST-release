! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_ev
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_ev *
!!
!!NAME
!! cq_ut_psp2pao_ev
!!PURPOSE
!! Recalculation of energy and potential
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_xc
!!  cq_ut_psp2pao_hartree
!!  cq_ut_psp2pao_density
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 14/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_ev

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine recalculate_energy
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_ev/recalculate_energy *
!!
!!NAME
!! recalculate_energy
!!USAGE
!!
!!PURPOSE
!! Recalculation of the energy
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_numeric
!!  cq_ut_psp2pao_xc
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 14/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine recalculate_energy()

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_numeric
     use cq_ut_psp2pao_xc

     implicit none

     !! Local variables

     integer :: i
     real(double) :: tmp_rho
     real(double) :: result

     ! Band energy
     real(double) :: e_band

     ! Double counting terms
     real(double) :: dcE_hartree, dcE_xc

     !! Parameters
     real(double), parameter :: two_pi = 2.0_double * 3.1415926535897932_double

     e_band = 0.0_double
     do i=1, 2*gl_no_orbitals
        e_band = e_band + gl_occ(i) * gl_eigenvalues(i)
     end do

     ! Double counting: Hartree
     do i=1, gl_points_mesh
        tmp_rho = gl_rho(i) + gl_rho(i + gl_points_mesh)
        gl_tmp1(i) = two_pi * gl_r(i) * gl_r(i) *tmp_rho * gl_v_hartree(i)
     end do

     call simpson_integral_prod(gl_r, gl_tmp1, gl_points_mesh, gl_dnu, result)
     dcE_hartree = -result

     ! Double counting: Exchange-correlation
     ! (Result in dcE_xc)
     call get_dcE_xc(dcE_xc)


     ! Total energy
     gl_e_total = e_band + dcE_hartree + dcE_xc

  end subroutine recalculate_energy
!!***


! -----------------------------------------------------------------------------
! Subroutine recalculate_potential
! -----------------------------------------------------------------------------

!!****f* cq_ut_psp2pao_ev/recalculate_potential *
!!
!!NAME
!! recalculate_potential
!!USAGE
!!
!!PURPOSE
!! Recalculation of the potential
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_hartree
!!  cq_ut_psp2pao_xc
!!  cq_ut_psp2pao_density
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 14/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine recalculate_potential

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_hartree
     use cq_ut_psp2pao_xc
     use cq_ut_psp2pao_density

     implicit none

     integer :: i

     call get_rho
     call get_v_hartree
     call get_v_xc

     do i=1,gl_points_mesh
        gl_v(i) = gl_v(i)  &
                + ( gl_v_nuclear(i) &
                  + gl_v_hartree(i) &
                  + gl_v_xc(i) &
                  + gl_v_ext(i) &
                  - gl_v(i) ) * gl_potential_mix
        gl_v(gl_points_mesh + i) = gl_v(gl_points_mesh + i)  &
                                 + ( gl_v_nuclear(i) &
                                   + gl_v_hartree(i) &
                                   + gl_v_xc(gl_points_mesh + i) &
                                   + gl_v_ext(i) &
                                   - gl_v(gl_points_mesh + i) ) * gl_potential_mix
     end do

  end subroutine recalculate_potential
!!***

end module cq_ut_psp2pao_ev


