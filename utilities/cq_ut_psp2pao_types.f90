! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module cq_ut_psp2pao_types
! ------------------------------------------------------------------------------
! $Log$
! ------------------------------------------------------------------------------

!!****h* Conquest/utilities/cq_ut_psp2pao_types *
!!
!!NAME
!! cq_ut_psp2pao_types
!!PURPOSE
!! Derived types for the utility
!!USES
!!
!!AUTHOR
!! Antonio S. Torralba
!!CREATION DATE
!! 23/02/2006
!!MODIFICATION HISTORY
!!
!!SOURCE
!!
module cq_ut_psp2pao_types

   implicit none

   integer, parameter :: double = selected_real_kind( 6, 70 )
!   integer, parameter :: double = selected_real_kind( 15, 307 )

!!****s* cq_ut_psp2pao_types/orbital_info *
!!NAME
!! orbital_info
!!PURPOSE
!! Defines a derived type that contains some information about an orbital
!!AUTHOR
!! Antonio S. Torralba
!!SOURCE
!!
   type orbital_info
      integer :: n
      integer :: l
      real(double) :: occ
      logical :: keep 
   end type orbital_info
!!***

!!****s* cq_ut_psp2pao_types/psp_in_table_fhi *
!!NAME
!! psp_in_table_fhi
!!PURPOSE
!! Table for storing the pseudopotential data as read from a file
!! The file must have Fritz-Haber-Institute Abinit type 6 format
!!AUTHOR
!! Antonio S. Torralba
!!SOURCE
!!
   type psp_in_table_fhi
     ! Mesh data
     integer :: nmesh
     real(double) :: amesh
     real(double) :: lmesh

     ! Kleiman-Bylander integral
     real(double) :: kb_integral

     ! Radii and potentials
     real(double) :: cutoff
     real(double), dimension(:), pointer :: r
     real(double), dimension(:), pointer :: u
     real(double), dimension(:), pointer :: u2
     real(double), dimension(:), pointer :: v
     real(double), dimension(:), pointer :: v2

     ! Partial core
     real(double), dimension(:), pointer :: rpc
     real(double), dimension(:), pointer :: rhopc
     real(double), dimension(:), pointer :: rhopc2
   end type psp_in_table_fhi
!!***

!!****s* cq_ut_psp2pao_types/in_pseudopotential *
!!NAME
!! in_pseudopotential
!!PURPOSE
!! Pseudopotential data
!!AUTHOR
!! Antonio S. Torralba
!!SOURCE
!!
   type in_pseudopotential
     integer :: gri_points

     ! Radii
     real(double) :: r_max
     real(double), dimension(:), pointer :: r
     real(double), dimension(:), pointer :: v_nonlocal_cutoff
     real(double) :: v_local_cutoff

     ! Atom info
     integer :: valence_elec   ! Number of valence electrons
     integer :: psp_comp       ! Number of pseudopotential components
     logical :: partial_core   ! Is there a partial core?

     ! Non-local info
     integer, dimension(:), pointer :: l_nonlocal
     integer, dimension(:), pointer :: sign_nonlocal

     ! Potentials
     real(double), dimension(:), pointer :: v_nuclear
     real(double), dimension(:), pointer :: v_nuclear2
     real(double), dimension(:), pointer :: v_nonlocal
     real(double), dimension(:), pointer :: v_nonlocal2

     ! Partial core
     real(double), dimension(:), pointer :: rhopc
     real(double), dimension(:), pointer :: rhopc2
   end type in_pseudopotential
!!***

!!****s* cq_ut_psp2pao_types/basis_data *
!!NAME
!! basis_data
!!PURPOSE
!! Defines a derived type that contains the description of the basis
!!AUTHOR
!! Antonio S. Torralba
!!SOURCE
!!
   type basis_data
      integer :: no_orbitals
      integer :: points_mesh
      logical :: partial_core
      real(double) :: e_total
      real(double) :: core_charge
      real(double) :: dnu
      real(double), dimension(:), pointer :: r
      real(double), dimension(:), pointer :: rho
      real(double), dimension(:), pointer :: rhopc
      real(double), dimension(:), pointer :: rho2
      real(double), dimension(:), pointer :: rhopc2
      
      ! Potentials
      integer :: psp_comp  ! Number of nonlocal components
      integer, dimension(:), pointer :: l_nonlocal
      integer, dimension(:), pointer :: sign_nonlocal
      real(double), dimension(:), pointer :: v_local
      real(double), dimension(:), pointer :: v_nonlocal
      real(double), dimension(:), pointer :: v_local2
      real(double), dimension(:), pointer :: v_nonlocal2

      ! Cutoffs
      real(double) :: cutoff_radius  ! Maximum cutoff
      real(double) :: v_local_cutoff
      real(double), dimension(:), pointer :: v_nonlocal_cutoff

      ! Orbitals
      integer, dimension(:), pointer :: orb_n
      integer, dimension(:), pointer :: orb_l
      logical, dimension(:), pointer :: orb_keep
      real(double), dimension(:), pointer :: orb_cutoff_radius
      real(double), dimension(:), pointer :: orb_occ
      real(double), dimension(:), pointer :: orb_eigenvalues
      real(double), dimension(:,:), pointer :: orb_ul
      real(double), dimension(:,:), pointer :: orb_ul2
   end type basis_data
!!***


end module cq_ut_psp2pao_types

