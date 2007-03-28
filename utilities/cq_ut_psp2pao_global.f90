! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_global
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_global *
!!
!!NAME
!! cq_ut_psp2pao_global
!!PURPOSE
!! Global variables and data structures
!!USES
!! cq_ut_psp2pao_types
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 27/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_global

   use cq_ut_psp2pao_types

   implicit none

   !! Global variables

   ! Setup variables
   integer :: gl_xcflag
   integer :: gl_max_scf_cycles
   real(double) :: gl_potential_mix
   real(double) :: gl_orbital_mix

   ! Grid data
   real(double) :: gl_grid_rmin, gl_grid_rmax, gl_dnu
   integer :: gl_points_psp, gl_points_mesh

   ! Energy data
   real(double) :: gl_e_total


   ! Pseudopotential data
   type(in_pseudopotential) :: gl_psp_in

   ! Orbitals
   integer :: gl_no_orbitals
   type(orbital_info), dimension(:), pointer :: gl_orbitals
   real(double), dimension(:), pointer :: gl_ul
   real(double), dimension(:), pointer :: gl_eigenvalues

   ! Density, potentials, integrals...
   real(double), dimension(:), pointer :: gl_r
   real(double), dimension(:), pointer :: gl_rho
   real(double), dimension(:), pointer :: gl_v
   real(double), dimension(:), pointer :: gl_v_nuclear
   real(double), dimension(:), pointer :: gl_v_nonlocal
   real(double), dimension(:), pointer :: gl_v_nonlocal_int_in
   real(double), dimension(:), pointer :: gl_v_nonlocal_int_out
   real(double), dimension(:), pointer :: gl_v_ext
   real(double), dimension(:), pointer :: gl_v_hartree
   real(double), dimension(:), pointer :: gl_v_xc
   real(double), dimension(:), pointer :: gl_rhopc
   real(double), dimension(:), pointer :: gl_norm
   real(double), dimension(:), pointer :: gl_occ

   ! Cutoff
   real(double) :: gl_cutoff_radius                       ! Maximum cutoff, used for PAOs
   integer :: gl_local_cutoff_tolerance                   ! Defines the truncation of the local part tail

   ! Temporary vectors
   real(double), dimension(:), pointer :: gl_tmp1
   real(double), dimension(:), pointer :: gl_tmp2
   real(double), dimension(:), pointer :: gl_tmp3
   real(double), dimension(:), pointer :: gl_dtmp1
   real(double), dimension(:), pointer :: gl_dtmp2
   real(double), dimension(:), pointer :: gl_dtmp3

   ! Basis data
   type(basis_data), dimension(:), pointer :: gl_wf_set   ! For the separate wavefunction sets
   type(basis_data) :: gl_basis                           ! For the basis

end module cq_ut_psp2pao_global

