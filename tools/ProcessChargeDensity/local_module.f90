module local

  use datatypes

  ! These give the number of blocks in x, y and z
  integer :: nblockx, nblocky, nblockz 
  real(double) :: block_size_x, block_size_y, block_size_z, grid_x, grid_y, grid_z
  
  ! Processes used
  integer :: nprocs

  integer :: nkp, n_eval_window, n_bands_active
  integer, dimension(:), allocatable :: band_no
  real(double) :: efermi
  real(double), dimension(:), allocatable :: kx,ky,kz, wtk
  real(double), dimension(:,:), allocatable :: eigenvalues
  real(double), dimension(:,:,:), allocatable :: current

  character(len=50) :: root_file
  
  real(double) :: stm_bias, fermi_offset, stm_z_min, stm_z_max, stm_x_min, stm_x_max, &
       stm_y_min, stm_y_max, stm_broad, gpv, E_wf_min, E_wf_max
  integer :: nptsx, nptsy, nptsz, nxmin, nymin, nzmin
  integer :: nrptx, nrpty, nrptz, nsampx, nsampy, nsampz

  integer :: flag_output
  integer, parameter :: dx = 1
  integer, parameter :: cube = 2
  
  logical :: flag_only_charge, flag_by_kpoint, flag_range
  character(len=80) :: charge_stub
  
  type block_set
     integer :: process, num_blocks
     integer, pointer, dimension(:) :: num, nx, ny, nz, active
  end type block_set

  type(block_set), allocatable, dimension(:) :: block_store
end module local
