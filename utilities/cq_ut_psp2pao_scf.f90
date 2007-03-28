! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_scf
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_scf *
!!
!!NAME
!! cq_ut_psp2pao_scf
!!PURPOSE
!! Self-consistent field related tasks
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_solve
!!  cq_ut_psp2pao_ev
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_scf

  implicit none

contains

! -----------------------------------------------------------------------------
! Subroutine scf_loop
! -----------------------------------------------------------------------------


!!****f* cq_ut_psp2pao_scf/scf_loop *
!!
!!NAME
!! scf_loop
!!USAGE
!!
!!PURPOSE
!! Iterations until self-consistence is achieved or
!!   the maximum number of iterations is reached
!!INPUTS
!!
!!USES
!!  cq_ut_psp2pao_global
!!  cq_ut_psp2pao_types
!!  cq_ut_psp2pao_solve
!!  cq_ut_psp2pao_ev
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 08/03/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
  subroutine scf_loop

     use cq_ut_psp2pao_global
     use cq_ut_psp2pao_types
     use cq_ut_psp2pao_solve
     use cq_ut_psp2pao_ev

     implicit none

     integer :: scf_cycle, i, tmp
     real(double) :: e_total_old

     ! Constants
     real(double), parameter :: tolerance = 1.0e-10

     ! Make sure that the first cycle happens
     gl_e_total = 0.0_double
     e_total_old = tolerance + 1.0_double
     scf_cycle = 0

     do while ((scf_cycle < gl_max_scf_cycles) &
         .and. (abs(gl_e_total - e_total_old) > tolerance))

        e_total_old = gl_e_total

        ! Only NON-POLARISED, for the moment
        call solve_eigenstates
        tmp = gl_no_orbitals * gl_points_mesh
        do i=1, tmp
           gl_ul(i + tmp) = gl_ul(i)
        end do
        do i=1, gl_no_orbitals
           gl_eigenvalues(i + gl_no_orbitals) = gl_eigenvalues(i)
        end do

        call recalculate_potential

        call recalculate_energy()

        write(*,*) "scf-cycle ",scf_cycle, " - Energy: ",2.0*gl_e_total

        scf_cycle = scf_cycle + 1
     end do

  end subroutine scf_loop
!!***

end module cq_ut_psp2pao_scf


