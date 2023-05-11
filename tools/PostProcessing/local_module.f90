module local

  use datatypes

  ! These give the number of blocks in x, y and z
  integer :: nblockx, nblocky, nblockz 
  real(double) :: block_size_x, block_size_y, block_size_z, grid_x, grid_y, grid_z
  
  ! Processes used
  integer :: nprocs

  integer :: nkp, n_eval_window, n_bands_active, n_bands_total, n_bands_process
  integer, dimension(:), allocatable :: band_no, active_bands, band_proc_no, band_full_to_active
  real(double), dimension(2) :: efermi ! Allow two Fermi levels with spin
  real(double), dimension(:), allocatable :: kx,ky,kz, wtk
  real(double), dimension(:,:,:), allocatable :: eigenvalues
  integer, dimension(:,:,:), allocatable :: band_active_kp
  integer, dimension(:,:), allocatable   :: band_active_all
  real(double), dimension(:,:,:), allocatable :: current

  ! Store eigenvector coefficients
  complex(double_cplx), allocatable, dimension(:,:,:,:,:), save :: evec_coeff ! PAOs, atoms, bands, kpoints, spin

  
  character(len=50) :: root_file
  
  real(double) :: stm_bias, fermi_offset, stm_z_min, stm_z_max, stm_x_min, stm_x_max, &
       stm_y_min, stm_y_max, stm_broad, gpv, E_wf_min, E_wf_max, E_procwf_min, E_procwf_max
  integer :: nptsx, nptsy, nptsz, nxmin, nymin, nzmin
  integer :: nrptx, nrpty, nrptz, nsampx, nsampy, nsampz

  integer :: flag_output
  integer, parameter :: dx = 1
  integer, parameter :: cube = 2
  
  logical :: flag_only_charge, flag_by_kpoint, flag_wf_range, flag_proc_range, flag_procwf_range_Ef
  logical :: flag_total_iDOS, flag_write_forces, flag_write_spin_moments
  character(len=80) :: charge_stub

  integer :: i_job ! Job type
  integer :: coord_format ! Output format: xyz (1) or cell (2)

  type block_set
     integer :: process, num_blocks
     integer, pointer, dimension(:) :: num, nx, ny, nz, active
  end type block_set

  type(block_set), allocatable, dimension(:) :: block_store

  ! From DiagModule
  real(double) :: kT
  ! Flags controlling Methfessel-Paxton approximation to step-function
  integer :: flag_smear_type, iMethfessel_Paxton
  
end module local
