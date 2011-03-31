! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!  SOURCE
!!
module global_module

  ! Module usage
  use datatypes

  implicit none

  ! RCS tag for object file identification 
  character(len=80), save, private :: RCSid = "$Id$"

  integer :: iprint                 ! Level of output
  integer :: io_lun                 ! Output unit
  integer, allocatable, dimension(:) :: id_glob      ! global label of atom in sim cell (CC)
  integer, allocatable, dimension(:) :: id_glob_inv  ! gives global number for a CC atom
  integer, dimension(:), allocatable :: species_glob ! gives species 
  integer :: numprocs               ! number of processors
  real(double) :: rcellx,rcelly,rcellz  ! cell side lengths
  real(double), allocatable, dimension(:) :: x_atom_cell ! position of atom in sim cell (CC)
  real(double), allocatable, dimension(:) :: y_atom_cell
  real(double), allocatable, dimension(:) :: z_atom_cell
  integer, allocatable, dimension(:) :: coord_to_glob
  integer :: ni_in_cell ! Atoms in cell
  integer :: ne_in_cell ! Electrons in cell
  ! atom_coord : Use global labelling, in the future this array should
  ! be used instead of x, y, z_atom_cell. by T. Miyazaki
  real(double), dimension(:,:), allocatable :: atom_coord ! Atomic coordinates
  integer, dimension(:,:), allocatable :: sorted_coord ! Atom IDs of atoms sorted according to x, y, z coords
  logical, dimension(:,:), allocatable :: flag_move_atom  ! Move atoms ?
  logical :: restart_L, restart_rho

  integer :: global_maxatomspart ! Maximum atoms per partition, if exceeded, triggers partitioning refinement

  integer :: load_balance
  logical :: many_processors ! Selects appropriate algorithm for partitioning

  character(len=20), save :: runtype ! What type of run is it ?

  ! Logical flags controlling run
  logical :: flag_vary_basis, flag_self_consistent, flag_residual_done
  ! Logical flags controlling preconditioning
  logical :: flag_precondition_blips
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
  logical :: flag_mix_L_SC_min
  logical :: flag_onsite_blip_ana
  logical :: flag_read_velocity   ! 16/06/2010 TM
  logical :: flag_quench_MD       ! 25/06/2010 TM
  real(double) :: temp_ion        ! 25/06/2010 TM

  ! How should blocks be assigned ? See block_module.f90
  integer :: flag_assign_blocks

  ! Number of L iterations
  integer :: max_L_iterations

  ! Numerical flag choosing basis sets
  integer :: flag_basis_set
  integer, parameter :: blips = 1
  integer, parameter :: PAOS = 2

  ! Numerical flag choosing functional type
  integer :: flag_functional_type
  character(len=15) :: functional_description 
  integer, parameter :: functional_lda_pz81        = 1
  integer, parameter :: functional_lda_gth96       = 2
  integer, parameter :: functional_lda_pw92        = 3
  integer, parameter :: functional_gga_pbe96       = 101  ! Standard PBE
  integer, parameter :: functional_gga_pbe96_rev98 = 102  ! revPBE (PBE + Zhang-Yang 1998)
  integer, parameter :: functional_gga_pbe96_r99   = 103  ! RPBE (PBE + Hammer-Hansen-Norskov 1999)

  ! Switch for variation of blips in get_support_gradient
  integer :: WhichPulay
  integer, parameter :: PhiPulay = 1
  integer, parameter :: SPulay = 2
  integer, parameter :: BothPulay = 3 

  ! What are the local functions ? 
  integer, parameter :: sf = 1  ! Support functions
  integer, parameter :: nlpf = 2  ! Projector functions
  integer, parameter :: paof = 3 ! Pseudo-atomic orbitals
  integer, parameter :: dens = 4 ! Atomic charge density

  ! Define areas of the code
  integer, parameter :: n_areas = 12
  integer, parameter :: area_init = 1
  integer, parameter :: area_matrices = 2
  integer, parameter :: area_ops = 3
  integer, parameter :: area_DM = 4 
  integer, parameter :: area_SC = 5
  integer, parameter :: area_minE = 6
  integer, parameter :: area_moveatoms = 7
  integer, parameter :: area_index = 8
  integer, parameter :: area_general = 9
  integer, parameter :: area_pseudo = 10
  integer, parameter :: area_basis = 11
  integer, parameter :: area_integn = 12
  integer :: iprint_init, iprint_mat, iprint_ops, iprint_DM, iprint_SC, iprint_minE, &
             iprint_MD, iprint_index, iprint_gen, iprint_pseudo, iprint_basis, &
             iprint_intgn, iprint_time

  integer, parameter :: IPRINT_TIME_THRES0 = 0  ! Always print
  integer, parameter :: IPRINT_TIME_THRES1 = 2  ! Important local timers
  integer, parameter :: IPRINT_TIME_THRES2 = 4  ! Not that important
  integer, parameter :: IPRINT_TIME_THRES3 = 6  ! For special purposes

  ! For P.C.C.
  logical :: flag_pcc_global = .false.
  
end module global_module
!!***
