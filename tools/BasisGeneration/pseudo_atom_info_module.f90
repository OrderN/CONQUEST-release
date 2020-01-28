! Store atomic data - semi-local potentials, KB projectors, local potential,
! PAOs, data on valence electrons needed for PAO creation
module pseudo_atom_info

  use datatypes
  
  implicit none

  type potential_vkb
     integer :: ngrid, ngrid_vkb
     integer :: nloc ! Index for local potential
     integer, dimension(0:6) :: n_proj
     integer :: n_nl_proj
     real(double) :: r_vkb
     real(double), pointer, dimension(:) :: rr, charge, pcc, local
     real(double), pointer, dimension(:,:,:) :: projector ! i, proj, ell
     real(double), pointer, dimension(:,:) :: semilocal_potential
     real(double), dimension(0:6) :: core_radius ! We won't go to l=6 but want to store easily
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
     integer :: polarised_l, polarised_n, polarised_shell
     integer :: flag_cutoff ! 1 for radii, 2 for energies, 3 for default (mix of all)
     logical :: flag_perturb_polarise
     integer :: flag_zetas ! 1 for split-norm, 2 for compressed energies
     integer, pointer, dimension(:) :: nzeta, l, n, npao, l_no, inner
     real(double) :: pol_pf
     real(double), pointer, dimension(:,:) :: cutoff, energy
     type(radial_table), pointer, dimension(:,:) :: psi, psi_reg
     type(radial_table), pointer, dimension(:) :: pol, pol_reg
     real(double), pointer, dimension(:) :: width, prefac ! Parameters for exponential confinement
  end type pao

  type valence_info
     integer, pointer, dimension(:) :: n, l, semicore, nkb, npao, inner
     integer :: n_shells, n_occ, functional
     real(double), pointer, dimension(:) :: occ, en_ps, en_pao
  end type valence_info

  type(pao) :: paos
  type(potential_vkb) :: local_and_vkb
  type(valence_info) :: val
  logical :: flag_default_cutoffs

  integer :: pseudo_type

  logical :: flag_plot_output, flag_use_Vl

  real(double) :: deltaE_large_radius! = 0.00073498_double
  real(double) :: deltaE_small_radius! = 0.073498_double
  
  ! Useful parameters to improve code readability
  integer, parameter :: pao_cutoff_energies = 1
  integer, parameter :: pao_cutoff_radii = 2
  integer, parameter :: pao_cutoff_default = 3 ! This will be set properly at some point
  
  character(len=80), dimension(80) :: input_file ! Store input file and append at end of ion file
  integer :: input_file_length
  integer :: hamann_version

contains

  subroutine allocate_pao(n_shells)

    use numbers, only: zero, one
    
    implicit none

    integer :: n_shells

    paos%n_shells = n_shells
    allocate(paos%nzeta(n_shells),paos%l(n_shells), &
            paos%n(n_shells),paos%npao(n_shells), &
            paos%inner(n_shells),&
            paos%width(n_shells),paos%prefac(n_shells),paos%l_no(n_shells))
    paos%n = 0
    paos%npao = 0
    paos%nzeta = 0
    paos%l = 0
    paos%l_no = 0
    paos%lmax = 0
    paos%width = one
    paos%prefac = zero
    paos%inner = 0
    paos%pol_pf = one
    return
  end subroutine allocate_pao

  subroutine allocate_pao_z(maxz)

    use numbers, only: zero
    
    implicit none

    integer :: maxz

    allocate(paos%cutoff(maxz,paos%n_shells))
    paos%cutoff = zero
    allocate(paos%energy(maxz,paos%n_shells))
    paos%energy = zero
    allocate(paos%psi(maxz,paos%n_shells))
    allocate(paos%psi_reg(maxz,paos%n_shells))
    return
  end subroutine allocate_pao_z
  
  subroutine deallocate_pao

    implicit none

    integer :: i_shell, zeta
    
    ! Deallocate internal tables
    do i_shell = 1,paos%n_shells
       do zeta = 1,paos%nzeta(i_shell)
          deallocate(paos%psi(zeta,i_shell)%f,paos%psi_reg(zeta,i_shell)%f,paos%psi_reg(zeta,i_shell)%x)
       end do
    end do
    ! Deallocate maxz terms
    deallocate(paos%cutoff,paos%energy,paos%psi,paos%psi_reg)
    ! Deallocate single terms
    deallocate(paos%nzeta,paos%l, paos%n,paos%npao, &
            paos%inner, paos%width,paos%prefac)
    paos%n_shells = 0
    return
  end subroutine deallocate_pao

  subroutine allocate_val(n_shells)

    use numbers, only: zero

    implicit none

    integer :: n_shells

        val%n_shells = n_shells
    allocate(val%n(n_shells),val%npao(n_shells),val%l(n_shells), &
         val%occ(n_shells),val%semicore(n_shells),val%en_ps(n_shells), &
         val%en_pao(n_shells),val%inner(n_shells))
    val%n = 0
    val%npao = 0
    val%l = 0
    val%occ = zero
    val%semicore = 0
    val%inner = 0
    val%en_ps = zero
    val%en_pao = zero
    return
  end subroutine allocate_val

  subroutine deallocate_val

    implicit none

    deallocate(val%n,val%npao,val%l, val%occ,val%semicore,val%en_ps, &
         val%en_pao,val%inner)
    return
  end subroutine deallocate_val

  subroutine allocate_vkb(ngrid,i_species)

    use numbers, only: zero
    use pseudo_tm_info, ONLY: pseudo

    implicit none

    ! Passed variables
    integer :: ngrid, i_species

    ! Local variables
    integer :: i

    local_and_vkb%ngrid = ngrid
    allocate(local_and_vkb%rr(ngrid))
    local_and_vkb%rr = zero
    allocate(local_and_vkb%charge(ngrid))
    local_and_vkb%charge = zero
    allocate(local_and_vkb%local(ngrid))
    local_and_vkb%local = zero
    i = maxval(local_and_vkb%n_proj)
    allocate(local_and_vkb%projector(ngrid,i,0:pseudo(i_species)%lmax))
    local_and_vkb%projector = zero
    allocate(local_and_vkb%semilocal_potential(ngrid,0:pseudo(i_species)%lmax))
    local_and_vkb%semilocal_potential = zero
    if(pseudo(i_species)%flag_pcc) then
       allocate(local_and_vkb%pcc(ngrid))
       local_and_vkb%pcc = zero
    end if
    return
  end subroutine allocate_vkb

  subroutine deallocate_vkb(i_species)

    use pseudo_tm_info, ONLY: pseudo

    implicit none

    ! Passed variables
    integer :: i_species

    deallocate(local_and_vkb%rr)
    deallocate(local_and_vkb%charge)
    deallocate(local_and_vkb%local)
    deallocate(local_and_vkb%projector)
    deallocate(local_and_vkb%semilocal_potential)
    if(pseudo(i_species)%flag_pcc) then
       deallocate(local_and_vkb%pcc)
    end if
    return
  end subroutine deallocate_vkb
  
end module pseudo_atom_info
