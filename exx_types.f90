! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: $
! -----------------------------------------------------------
! Module exx_types.f90
! -----------------------------------------------------------
! Code area 13: EXX
! -----------------------------------------------------------

!!****h* Conquest/exx_types *
!!  NAME
!!   exx_types
!!
!!  PURPOSE
!!   Contain all the global variables and 
!!   derived-types for EXX
!!
!!  USES
!!   datatypes and cq_timer derived-type 
!!   from timer_module
!! 
!!  AUTHOR
!!   L.A. Truflandier
!!  CREATION DATE
!!   2011/02/11
!!
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
module exx_types 

  use datatypes
  use timer_module,           only: cq_timer
!  use timer_stdclocks_module, only: tmr_std_exx

  type fftw3d
     !     complex(double_cplx), allocatable, dimension(:,:,:) :: auxin
     !     complex(double_cplx), allocatable, dimension(:,:,:) :: auxout
     complex(double_cplx), allocatable, dimension(:,:,:) :: arrayin
     complex(double_cplx), allocatable, dimension(:,:,:) :: arrayout
     integer*8 :: planF, planR
  end type fftw3d

  !  type fftw2d
  !     complex(double_cplx), allocatable, dimension(:,:) :: arrayin
  !     complex(double_cplx), allocatable, dimension(:,:) :: arrayout
  !     integer*8 :: planF, planR
  !  end type fftw2d
  
  !  type fftw1d
  !     complex(double_cplx), allocatable, dimension(:) :: arrayin
  !     complex(double_cplx), allocatable, dimension(:) :: arrayout
  !     integer*8 :: planF, planR
  !  end type fftw1d

  ! PAOs on grid
  real(double), dimension(:,:,:,:),     allocatable :: phi_i
  real(double), dimension(:,:,:,:),     allocatable :: phi_j
  real(double), dimension(:,:,:,:),     allocatable :: phi_k
  real(double), dimension(:,:,:,:),     allocatable :: phi_l  

  ! Auxiliary densities and potentials
  real(double), dimension(:,:,:,:,:),   allocatable :: rho_kj
  real(double), dimension(:,:,:,:,:),   allocatable :: Ome_kj
  real(double), dimension(:,:,:,:,:),   allocatable :: vhf_kj
  real(double), dimension(:,:,:,:),     allocatable :: Phy_k

  ! Work matrix
  real(double), dimension(:,:,:),       allocatable :: work_in_3d
  real(double), dimension(:,:,:),       allocatable :: work_out_3d

  ! For the Poisson equation/FFTW in reciprocal space
  type(fftw3d)          :: fftwrho3d

  complex(double), dimension(:,:,:),    allocatable :: reckernel_3d
  real(double),    dimension(:,:,:),    allocatable :: ewald_rho
  real(double),    dimension(:,:,:),    allocatable :: ewald_pot
  real(double) :: ewald_charge


  ! For the Poisson equation/ISF in real space
  real(double), dimension(:,:,:),       allocatable :: isf_rho
  real(double), dimension(:,:,:),       allocatable :: isf_pot_ion
  real(double), pointer   :: kernel(:)   

  ! Poisson solver settings
  character(100), parameter :: p_scheme_default = 'v(G=0)=0'
  character(100)            :: p_scheme
  character(100)            :: exx_psolver
  real(double)              :: p_cutoff
  real(double)              :: p_factor
  real(double)              :: p_omega

  real(double)              :: pulay_factor
  real(double)              :: pulay_radius

  integer,        parameter :: p_ngauss = 89     ! do not touch !
  real(double)              :: p_gauss(p_ngauss)
  real(double)              :: w_gauss(p_ngauss)

  real(double),   parameter :: magic_number = 0.1864d0 ! 32/43 = nombre magic !
  real(double)              :: ewald_alpha
  integer                   :: isf_order

  real(double)   :: edge, volume, screen

  ! Grid settings
  integer                 :: extent
  integer                 :: ngrid
  real(double)            :: r_int
  real(double)            :: grid_spacing

  ! Timers
  !type(cq_timer), save :: tmr_std_exx
  type(cq_timer), save :: tmr_std_exx_setup
  type(cq_timer), save :: tmr_std_exx_write
  type(cq_timer), save :: tmr_std_exx_kernel
  type(cq_timer), save :: tmr_std_exx_fetch
  type(cq_timer), save :: tmr_std_exx_fetch_K
  type(cq_timer), save :: tmr_std_exx_accumul

  type(cq_timer), save :: tmr_std_exx_matmult
  type(cq_timer), save :: tmr_std_exx_quadrat
  type(cq_timer), save :: tmr_std_exx_poisson
  type(cq_timer), save :: tmr_std_exx_allocat
  type(cq_timer), save :: tmr_std_exx_dealloc

  type(cq_timer), save :: tmr_std_exx_evalpao
  type(cq_timer), save :: tmr_std_exx_evalgto
  type(cq_timer), save :: tmr_std_exx_splitpao
  real(double)         :: exx_total_time
  
  ! User settings
  integer :: exx_scheme      ! 4center ERIs or 3center reduction integrals
  integer :: exx_mem         ! reduced memory allocation 
  logical :: exx_phil        ! on-the-fly or storage of phi_l
  logical :: exx_phik        ! on-the-fly or storage of phi_k
  logical :: exx_overlap     ! compute overlap local boxes
  logical :: exx_alloc       ! on-the-fly or global memory allocation
  logical :: exx_cartesian   ! cartesian or spherial representation for PAO
  logical :: exx_screen      ! screening
  logical :: exx_screen_pao  ! method for screening
  logical :: exx_gto         ! testing
  real(double) :: exx_cutoff ! cutoff for screening (experimental)
  real(double) :: exx_radius ! radius for integration
  real(double) :: exx_hgrid  ! radius for integration

  ! For debuging/testing
  logical :: exx_debug 
  logical :: exx_Kij
  logical :: exx_Kkl

  ! I/O
  integer :: unit_global_write 
  integer :: unit_matrix_write 
  integer :: unit_output_write 
  integer :: unit_timers_write
  integer :: unit_memory_write 
  integer :: unit_screen_write 
  !=================================================================<<

  type prim_atomic_data
     integer  :: pr

     integer  :: ip
     integer  :: num
     integer  :: spec
     integer  :: nsup
     integer  :: labcell
     real(double)               :: radi
     character(len=2)           :: name
     real(double)               :: radi_ip
     real(double), dimension(3) :: xyz_ip
     !
     real(double), dimension(3) :: xyz
     real(double) :: r
     real(double) :: d
  end type prim_atomic_data

  type neigh_atomic_data
     integer  :: nb

     integer  :: ist
     integer  :: npc
     integer  :: nic
     integer  :: global_part
     integer  :: global_num
     integer  :: spec
     integer  :: nsup
     real(double)               :: radi
     character(len=2)           :: name
     real(double), dimension(3) :: xyz_nb
     real(double), dimension(3) :: xyz_hl
     real(double), dimension(3) :: xyz_cv
     real(double), dimension(3) :: xyz_ipnb
     real(double) :: r_ipnb
     real(double) :: d_ipnb
     !
     real(double), dimension(3) :: xyz
     real(double) :: r
     real(double) :: d
  end type neigh_atomic_data

end module exx_types
