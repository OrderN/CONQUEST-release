! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module global_module
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/global_module *
!!  NAME
!!   global_module
!!  PURPOSE
!!   Holds various global variables
!!  USES
!!   datatypes
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    ROBODoc header, removed unnecessary variables
!!   18/03/2002 dave
!!    Added RCS Id and Log tags and static tag for object file id
!!   13:49, 10/02/2003 drb 
!!    Added flags for minimisation control
!!   10:58, 2003/06/10 dave
!!    Added new flags for different options and new atomic coordinate variables
!!   12:19, 29/08/2003 drb 
!!    Added flag_move_atom
!!   13:33, 22/09/2003 drb 
!!    Added flags to allow separate testing of S-Pulay and phi-Pulay forces
!!   08:31, 2003/10/01 dave
!!    Changed flag_vary_blips to flag_vary_basis
!!   14:56, 02/05/2005 dave 
!!    Added global ne_in_cell variable for total electron number in cell
!!   09:11, 11/05/2005 dave 
!!    Added max L iterations parameter
!!   2006/09/04 08:06 dave
!!    Dynamic allocation implemented for final variables
!!   2006/10/10 10:05 ast
!!    Flag for selecting functional and string for its description
!!   2007/03/23 17:16 dave
!!    Added new variables for automatic partitioning
!!   2007/04/18 17:26 dave
!!    Added flag for block assignment
!!   2008/02/01 03:43 dave
!!    Added output unit (output to file rather than stdout)
!!   12:19, 14/02/2008 drb 
!!    Added flag for Pulay relaxation algorithm
!!   2008/07/16 ast
!!    New iprint levels for timing 
!!   2009/07/24 16:41 dave
!!    Added new flag for global or per atom tolerances
!!   2011/03/30 18:59 M.Arita
!!    Added new flag for for P.C.C.
!!   2011/04/01 L.Tong
!!    Added flag_spin_polarisation as a switch for spin polarised calculation
!!    Added flag functional_lsda_pw92 for LSDA
!!   2011/07/26 L.Tong
!!    Added flag_fix_spin_population as switch for fixing spin population
!!   Friday, 2011/08/05 L.Tong
!!    Added initial (fixed) electron numbers for spin up and down
!!    channels, used when flag_fix_spin_population is true
!!    Moved all new variables for spin polarisation calculations together
!!   2011/09/29 14:51 M. Arita
!!    Added new flags for DFT-D2
!!   2011/07/21 16:35 dave
!!    Flags for cDFT
!!   2011/12/12 17:26 dave
!!    Flag for analytic blip integrals
!!   2012/03/07 L.Tong
!!    Added some more flags for spin polarisation, and uses numbers module
!!   2012/03/27 L.Tong
!!   - Added variable nspin
!!   - Added variable ne_spin_in_cell(nspin). This replaces
!!     ne_up_in_cell and ne_dn_in_dell
!!   - Added variable spin_factor
!!   - Default values are:
!!     nspin = 1
!!     spin_factor = two
!!   - removed now obsolete flag: flag_spin_polarisation
!!   2012/05/29 L.Tong
!!   - removed functional_lsda_pw92, now redundant. Just use
!!     functional_lda_pw92 for PW92 LDA.
!!   2012/06/24 L.Tong
!!   - Added flag flag_dump_L for controlling if L is to be dumped
!!   2013/01/30 10:30 dave
!!   - Adding deltaSCF variables (with U. Terranova)
!!   2013/07/01 M.Arita
!!   - Added flags and parameters for the efficient MD scheme
!!   2013/08/20 M.Arita
!!   - Added flags and variables for matrix reconstruction
!!   2013/12/02 M.Arita
!!   - Added flags and variables for XL-BOMD
!!   2014/01/17 lat 
!!    Added new area and flag for EXX 
!!   2014/09/20 lat
!!    Added flags for PBE0, Xalpha and Hartree-Fock functional
!!   2014/10/03 lat
!!    Added parameters for SCF control of EXX and iprint_exx
!!   2015/05/11 L.Truflandier
!!   - Added optional total spin magnetization ne_magn_in_cell
!!   2015/05/29
!!    Wavefunction output flags (COR and dave)
!!   2015/06/19
!!    FIRE implmementation (SA, COR, dave)
!!   2015/07/08 08:03 dave
!!    DOS and k-point by k-point wavefunction output (for STM)
!!   2015/11/09 08:23 dave (with TM, NW of Mizuho)
!!    Added neutral atom flag
!!   2016/02/16 JKS
!!    Added Wu-Cohen XC functional  (PRB 73, 235116  (2006) )
!!   2016/08/01 17:30 nakata
!!    Introduced atomf
!!   2016/08/09 21:30 nakata
!!    Added parameters for Contracted SFs and multi-site SFs
!!   2017/02/23 dave
!!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
!!   2017/04/05 18:00 nakata
!!    Added flag_readAtomicSpin to initialise spin
!!   2017/08/29 jack baker & dave
!!    Adding variables for cell optimisation
!!   2017/10/20 09:19 dave
!!    Moved fire variables to Integrators_module
!!   2017/11/13 18:15 nakata
!!    Added a flag to normalise pDOS
!!   2017/12/05 09:59 dave with TM & NW (MIZUHO)
!!    Added new function type - NA projector function (napf)
!!   2018/04/25 10:00 zamaan
!!    Added target attribute to rcellx, x_atom_cell etc.
!!   2018/05/17 12:51 dave with Ayako Nakata
!!    Changed flag_readAtomicSpin to flag_InitialAtomicSpin (more descriptive) and moved to density_module
!!   2018/09/19 18:30 nakata
!!    Added a flag for orbital angular momentum resolved PDOS
!!   2018/10/22 14:25 dave & jsb
!!    Adding (l,m)-projection for PDOS
!!   2019/02/28 zamaan
!!    Added enthalpy and stress tolerances for cell optimisation
!!   2019/03/28 zamaan
!!    Added flag_stress and flag_full_stress
!!   2019/05/08 zamaan
!!    Added flag_atomic_stress and atomic_stress for atomic contributions to 
!!    stress and heat flux
!!   2019/05/21 zamaan
!!    Added RNG seed
!!   2019/11/14 tsuyoshi
!!    Removed n_proc_old and glob2node_old
!!   2019/11/18 tsuyoshi
!!    Removed flag_MDold
!!   2019/11/18 14:37 dave
!!    Added flag_variable_cell
!!  SOURCE
!!
module global_module

  ! Module usage
  use datatypes
  use numbers

  implicit none

  integer :: iprint                 ! Level of output
  integer :: io_lun                 ! Output unit
  integer, allocatable, dimension(:) :: id_glob      ! global label of atom in sim cell (CC)
  integer, allocatable, dimension(:) :: id_glob_inv  ! gives global number for a CC atom
  integer, dimension(:), allocatable :: species_glob ! gives species 
  integer :: numprocs               ! number of processors
  real(double), target :: rcellx,rcelly,rcellz  ! cell side lengths
  real(double), allocatable, dimension(:), target :: x_atom_cell ! position of atom in sim cell (CC)
  real(double), allocatable, dimension(:), target :: y_atom_cell
  real(double), allocatable, dimension(:), target :: z_atom_cell
  integer,      allocatable, dimension(:), target :: coord_to_glob
  integer      :: ni_in_cell ! Atoms in cell
  real(double) :: ne_in_cell ! Electrons in cell
  ! atom_coord : Use global labelling, in the future this array should
  ! be used instead of x, y, z_atom_cell. by T. Miyazaki
  real(double), dimension(:,:), allocatable, target :: atom_coord ! Atomic coordinates
  integer,      dimension(:,:), allocatable :: sorted_coord ! Atom IDs of atoms sorted according to x, y, z coords
  logical,      dimension(:,:), allocatable :: flag_move_atom  ! Move atoms ?
  integer,      dimension(:),   allocatable :: flag_cdft_atom
  logical :: restart_DM, restart_rho, restart_T, restart_X

  integer :: global_maxatomspart ! Maximum atoms per partition, if exceeded, triggers partitioning refinement

  integer :: rng_seed

  integer :: load_balance
  logical :: many_processors ! Selects appropriate algorithm for partitioning

  character(len=20), save :: runtype ! What type of run is it ?

  logical :: flag_stress   ! Compute the stress tensor?
  logical :: flag_full_stress ! Compute the off-diagonal elements?
  logical :: flag_atomic_stress ! Compute atomic contributions to stress?
  logical :: flag_heat_flux ! Compute heat flux during MD?

  ! Atomic contributions to total stress
  ! I would rather not put this in global module, but it is required by enough
  ! different force computing modules that it's impossible to put it in a more
  ! sensible file without circular dependencies - zamaan
  real(double), dimension(:,:,:), allocatable :: atomic_stress
  real(double), dimension(3,3)                :: non_atomic_stress

  logical :: flag_opt_cell ! optimize the simulation cell?
  logical :: flag_variable_cell ! Global indicator of whether cell will change
  integer :: optcell_method ! method for cell optimiisation 1 = cell, fixed fractional coords, 2 = nested loop cell + geometry optimisation, 3 = single vector full optimisation
  ! specify sim cell dims/ratios of dims to be held constant.
  character(len=20), save :: cell_constraint_flag
  ! Termination condition to exit cell_cg_run
  real(double) :: cell_en_tol, cell_stress_tol


  ! Logical flags controlling run
  logical :: flag_vary_basis, flag_self_consistent, flag_residual_done
  ! Logical flags controlling preconditioning
  logical :: flag_precondition_blips
  logical :: flag_analytic_blip_int
  ! Logical flag controlling atomic coordinate format
  logical :: flag_fractional_atomic_coords
  logical :: flag_old_partitions
  logical :: flag_read_blocks ! Do we read make_blk.prt or raster ?
  logical :: flag_test_forces
  logical :: flag_reset_dens_on_atom_move
  logical :: flag_continue_on_SC_fail
  logical :: flag_SCconverged
  logical :: UseGemm
  logical :: flag_pulay_simpleStep
  logical :: flag_global_tolerance
  logical :: flag_Becke_weights
  logical :: flag_Becke_atomic_radii
  logical :: flag_perform_cDFT
  logical :: flag_mix_L_SC_min
  logical :: flag_onsite_blip_ana
  logical :: flag_read_velocity   ! 16/06/2010 TM
  logical :: flag_quench_MD       ! 25/06/2010 TM
  logical :: flag_fire_qMD        ! 2014/07/31 SA
  real(double) :: temp_ion        ! 25/06/2010 TM

  ! How should blocks be assigned ? See block_module.f90
  integer :: flag_assign_blocks

  ! Number of L iterations
  integer :: max_L_iterations

  ! Numerical flag choosing basis sets
  integer :: flag_basis_set
  integer, parameter :: blips = 1
  integer, parameter :: PAOS  = 2

  ! Switch for variation of blips in get_support_gradient
  integer :: WhichPulay
  integer, parameter :: PhiPulay  = 1
  integer, parameter :: SPulay    = 2
  integer, parameter :: BothPulay = 3 

  ! What are the local functions ? 
  integer, parameter :: sf   = 1 ! Support functions
  integer, parameter :: nlpf = 2 ! Projector functions
  integer, parameter :: paof = 3 ! Pseudo-atomic orbitals
  integer, parameter :: dens = 4 ! Atomic charge density
  integer, parameter :: napf = 5 ! Neutral atom projector functions
  integer            :: atomf    ! 1(=sf) for blips and primitive paos, 3(=paof) for contracted paos

  ! Define areas of the code
  integer, parameter :: n_areas        = 13
  integer, parameter :: area_init      = 1
  integer, parameter :: area_matrices  = 2
  integer, parameter :: area_ops       = 3
  integer, parameter :: area_DM        = 4 
  integer, parameter :: area_SC        = 5
  integer, parameter :: area_minE      = 6
  integer, parameter :: area_moveatoms = 7
  integer, parameter :: area_index     = 8
  integer, parameter :: area_general   = 9
  integer, parameter :: area_pseudo    = 10
  integer, parameter :: area_basis     = 11
  integer, parameter :: area_integn    = 12
  integer, parameter :: area_exx       = 13
  integer :: iprint_init, iprint_mat,     iprint_ops,   iprint_DM,    &
             iprint_SC,   iprint_minE,    iprint_MD,    iprint_index, &
             iprint_gen,  iprint_pseudo,  iprint_basis, iprint_intgn, &
             iprint_time, iprint_MDdebug, iprint_exx

  integer, parameter :: IPRINT_TIME_THRES0 = 0  ! Always print
  integer, parameter :: IPRINT_TIME_THRES1 = 2  ! Important local timers
  integer, parameter :: IPRINT_TIME_THRES2 = 4  ! Not that important
  integer, parameter :: IPRINT_TIME_THRES3 = 6  ! For special purposes

  integer :: min_layer ! Layer of minimisation algorithm (from 0 to -n)
  ! For P.C.C.
  logical :: flag_pcc_global = .false.

  !! For Spin polarised calculations (L.Tong)
  ! default to spin nonpolarised calculation
  ! number of spin channels
  integer      :: nspin = 1
  real(double) :: spin_factor = two
  ! Logical flag determine if spin populations are fixed (fixed magnetic moment) 
  logical      :: flag_fix_spin_population = .false.
  ! fixed electron numbers for different spin channels. This is used
  ! even for spin non-polarised case
  real(double), dimension(2) :: ne_spin_in_cell ! 1 = up, 2 = down
  real(double)               :: ne_magn_in_cell

  ! For DFT-D2
  logical :: flag_dft_d2
  logical :: flag_SCconverged_D2 = .false.
  logical :: flag_only_dispersion

  ! For vdwDFT
  logical :: flag_vdWDFT          ! selector for turning on vdW energy correction
  integer :: vdW_LDA_functional   ! selector for LDA functional

  ! DeltaSCF
  logical :: flag_DeltaSCF
  logical :: flag_excite = .false.
  logical :: flag_local_excitation
  integer :: dscf_source_level, dscf_target_level, dscf_source_spin, &
       dscf_target_spin, dscf_source_nfold, dscf_target_nfold, &
       dscf_homo_limit, dscf_lumo_limit
  real(double) :: dscf_HOMO_thresh, dscf_LUMO_thresh

  ! For EXX
  logical      :: flag_exx      = .false. ! switch on/off EXX
  integer      :: exx_scf       = 0       ! method used during the SCF using hybrid functional or Hartree-Fock
  real(double) :: exx_alpha     = zero    ! mixing factor for hybrid Exc
  
  integer      :: exx_niter     = 1       ! for EXX control during SCF
  integer      :: exx_siter     = 1       ! for EXX control during SCF
  real(double) :: exx_pulay_r0  = zero    ! get the R0 pulay residual for control
  real(double) :: exx_scf_ratio = zero    ! for EXX control during SCF
  real(double) :: exx_scf_tol   = zero    ! for EXX control during SCF

  ! pre-defined grid spacing in Bohr
  real(double), parameter :: exx_hgrid_coarse = 0.60_double
  real(double), parameter :: exx_hgrid_medium = 0.50_double
  real(double), parameter :: exx_hgrid_fine   = 0.20_double

  ! Flag to control if matrix L is dumped to files
  logical :: flag_dump_L
  logical :: flag_DumpMatrices

  ! Hold an old relation between global & partition labels
  integer,allocatable :: id_glob_old(:),id_glob_inv_old(:)

  ! For MD
  logical :: flag_LmatrixReuse
  logical :: flag_TmatrixReuse
  logical :: flag_SkipEarlyDM
  logical :: flag_MDcontinue
  logical :: flag_MDdebug
  logical :: flag_thermoDebug
  logical :: flag_baroDebug
  logical :: flag_FixCOM 
  integer :: McWFreq
  integer :: MDinit_step  
  !ORI real(double),parameter   :: shift_in_bohr = 1.0E-03_double
  real(double),parameter   :: shift_in_bohr = 1.0E-06_double
  ! Table showing atoms (global) in nodes
  integer,allocatable :: glob2node(:)        ! size: ni_in_cell
  ! Displacement of atoms from a previous step
  real(double),allocatable :: atom_coord_diff(:,:)
  ! XL-BOMD
  logical :: flag_XLBOMD
  logical :: flag_propagateX,flag_propagateL
  logical :: flag_dissipation
  character(20) :: integratorXL
  
  ! Wavefunction output
  logical :: flag_out_wf                        !output WFs?
  logical :: flag_out_wf_by_kp                  !output WFs k-point by k-point
  integer,allocatable,dimension(:)::out_wf      !which bands to output  
  integer::max_wf                               !total no of bands
  logical :: wf_self_con                        !flag to select output at the end of SCF cycle
  real(double) :: E_wf_min, E_wf_max            ! Limits for energy range
  logical :: flag_wf_range_Ef                   ! Is the energy range relative to Ef (T) or absolute (F)

  ! This is in the WF output section as I introduced it for WF output, but
  ! it more properly applies to matrices (specifically how many temporary matrices we can store)
  integer :: mx_temp_matrices                   ! Defaults to 100; used in mult_module (immi)
  
  ! DOS output (NB Maybe move these into DiagModule and revisit names)
  logical :: flag_write_DOS, flag_write_projected_DOS, flag_normalise_pDOS, flag_pDOS_angmom, flag_pDOS_lm
  real(double) :: E_DOS_min, E_DOS_max, sigma_DOS
  integer :: n_DOS

  ! Neutral atom potential
  logical :: flag_neutral_atom

  ! Contracted SF
  logical :: flag_SpinDependentSF
  integer :: nspin_SF
  logical :: flag_SFcoeffReuse

  ! Multisite
  logical :: flag_Multisite
  logical :: flag_LFD

  ! diagonalise or linear scaling
  logical :: flag_diagonalisation

end module global_module
!!***
