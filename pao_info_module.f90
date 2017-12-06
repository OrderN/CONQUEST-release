module pao_info

  use datatypes
  
  implicit none

  type potential
     integer :: ngrid
     real(double), pointer, dimension(:) :: rr, charge, pcc
     real(double), pointer, dimension(:,:) :: potential  
  end type potential

  type potential_vkb
     integer :: ngrid, ngrid_vkb
     integer :: nloc ! Index for local potential
     real(double) :: r_vkb
     real(double), pointer, dimension(:) :: rr, charge, pcc, local
     real(double), pointer, dimension(:,:,:) :: projector ! i, proj, ell
     real(double), pointer, dimension(:,:) :: semilocal_potential  
  end type potential_vkb

  type radial_table
     real(double), pointer, dimension(:) :: x ! Regular mesh for interpolation
     real(double), pointer, dimension(:) :: f
     real(double) :: delta
     integer :: n
  end type radial_table
  
  type pao
     integer :: n_shells, lmax
     integer :: total_paos
     integer :: polarised_l, polarised_n, polarised_shell, inner_shell
     integer :: flag_cutoff ! 1 for radii, 2 for energies, 3 for default (mix of all)
     logical :: flag_perturb_polarise
     integer :: flag_zetas ! 1 for split-norm, 2 for compressed energies
     integer, pointer, dimension(:) :: nzeta, l, n, npao, l_no, inner
     real(double), pointer, dimension(:,:) :: cutoff, energy
     type(radial_table), pointer, dimension(:,:) :: psi, psi_reg
     type(radial_table), pointer, dimension(:) :: pol, pol_reg
     logical, pointer, dimension(:) :: has_semicore
     real(double), pointer, dimension(:) :: width, prefac ! Parameters for exponential confinement
  end type pao

  type(pao), allocatable, dimension(:) :: paos

  type valence_info
     integer, pointer, dimension(:) :: n, l, semicore, nkb, npao, inner
     integer :: n_shells, n_occ, functional
     real(double), pointer, dimension(:) :: occ, en_ps, en_pao
     logical, pointer, dimension(:) :: has_semicore
  end type valence_info

  logical :: flag_default_cutoffs

end module pao_info
