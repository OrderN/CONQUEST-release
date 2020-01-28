!-*- mode: F90; mode: font-lock -*-
! -----------------------------------------------------------------------------
! $Id$
! -----------------------------------------------------------------------------
! DiagModule
! -----------------------------------------------------------------------------
! Code area 4: density matrix
! -----------------------------------------------------------------------------

!!****h* Conquest/DiagModule *
!!  NAME
!!   DiagModule - contains all routines needed to diagonalise the
!!     Hamiltonian
!!  PURPOSE
!!   Forms the Hamiltonian (by summing over images of local atoms),
!!   distributes it appropriately (i.e. according to Scalapack format)
!!   calls Scalapack and redistributes eigenvectors to build K matrix
!!
!!   This whole module is discussed in exhaustive detail in the
!!   Conquest notes "Implementation of Diagonalisation in Conquest"
!!   (Diagonalisation.tex)
!!  USES
!!   common, cover_module, datatypes, fdf, GenBlas, GenComms,
!!   global_module, group_module, matrix_data, matrix_module,
!!   maxima_module, mpi, numbers, primary_module, ScalapackFormat
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   15/02/2002
!!  MODIFICATION HISTORY
!!   19/02/2002 dave
!!    Happy 1st Birthday Christopher !
!!    Added more detail to rough outlining
!!   04/03/2002 dave
!!    Debugging
!!   07/03/2002 dave
!!    Rewrote various bits - previous implementation was badly thought
!!    out
!!   08/03/2002 dave
!!    Added Scalapack calls to FindEvals
!!   05/04/2002 dave
!!    Many changes over the last month - the code now works for 8 and
!!    64 atoms on 1,2,3 and 4 processors in every grid combination I
!!    can think of.  Main changes were to PrepareSend and
!!    DistributeCQtoSC.  There will be more changes as k-point
!!    strategy is implemented and as we build the K matrix.  Also
!!    generalised to work row-by-row (i.e. not making blocks and atoms
!!    commensurate).
!!   16/04/2002 dave
!!    Added distance between atoms to loc, and added k-point reading
!!    and loops to FindEvals
!!   16/04/2002 dave
!!    Shifted reading out of FindEvals into readDiagInfo
!!   23/04/2002 dave
!!    Started creation of occupancy and K matrix building routines
!!   23/04/2002 dave
!!    Defined a new derived type to hold arrays relating to
!!    distribution - this will allow easy incorporation into Conquest
!!   01/05/2002 drb
!!    Imported to Conquest
!!   17/06/2002 dave
!!    Fixed bugs, tidied (moved initialisation into separate routine
!!    etc), added logical diagon for solution method choice (maybe
!!    should be elsewhere !)
!!   31/07/2002 dave
!!    Added calculation for the matrix M12 which makes the Pulay force
!!    (contribution to dE/dphi_i for the variation of S)
!!   15:37, 03/02/2003 drb
!!    Various bits of tidying, increase of precision and general
!!    improvement
!!   08:20, 2003/07/28 dave
!!    Fairly major reworking to reflect need for multiple CQ elements
!!    to fold into one SC element.  New derived type (element),
!!    renamed and reworked parts of DistributeData type, reworked
!!    PrepareSend and DistributeCQ_to_SC to reflect all this.
!!   11:48, 30/09/2003 drb
!!    Changed iprint levels throughout
!!   2007/08/13 17:27 dave
!!    Changed kT to be user-set parameter in input file (Diag.kT keyword)
!!   2007/10/15 Veronika
!!    Added keyword maxefermi: Max number of iteration when searching
!!    for E_Fermi
!!   2008/02/01 17:46 dave
!!    Changes for output to file not stdout
!!   2008/07/31 ast
!!    Added timers
!!   2010/06/14 21:53 lt
!!    Added flags for Methfessel-Paxton approximation for step-function
!!    and a further flag to choose smearing type
!!   2010/06/15 17:03 lt
!!    Added the max_brkt_iterations flag so that in the bracket search
!!    algorithm in FindFermi() if the search for uppe/lower bound is
!!    unsucessful above the max_brkt_iterations, then the search
!!    restarts with a smaller incEf
!!   2010/07/26 lt
!!    Added erfc function.  This is a modified verion of the one in
!!    ewald_module, and works for all x
!!   2011/11/29 L.Tong
!!    Added spin polarisation
!!    - added SCH_dnmat, occ_up, occ_dn, w_up, w_dn, local_w_up,
!!      local_w_dn, Efermi_up, Efermi_dn
!!   2012/02/29 L.Tong
!!   - Added interface for findFermi
!!   2012/03/08 L.Tong
!!   - Major rewrite of implementation of spin polarisation
!!   - Making SCHmat, SCSmat, z, occ, w, and local_w to have one extra
!!     index for spin
!!   - Efermi(nspin) now become arrays storing the corresponding
!!     values in each spin channel.
!!   - Make Scalapack dimensions row_size and col_size module
!!     global. It makes more sense this way, otherwise the matrix
!!     dimensions have to be calculated in every subroutine that needs
!!     them.
!!   - Added new subroutine endDiag to handle deallocations of module
!!     global arrays.
!!   2013/02/08 08:31 dave (with UT)
!!   - Adding band_ef variable to track HOMO location (defined as 0K HOMO)
!!     along with several other DeltaSCF changes
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2015/06/05 16:43 dave and cor
!!    Writing out band-by-band charge density
!!   2016/02/09 08:16 dave
!!    Moved erfc to functions module
!!   2017/02/23 dave
!!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
!!   2017/06/22 13:41 dave
!!    Rearranged initialisation so that BLACS is only started once, and allocations are now split
!!    into once only (init_blacs_pg), once for every atom move (init/end_scalapack_format) now
!!    called from updateIndices (all routines) and once for each diagonalisation session (initDiag).
!!    To perform a set of diagonalisations call initDiag, distrib_and_diag, endDiag
!!   2017/10/18 11:06 dave
!!    Change to allow arbitrary k-points to be passed to distrib_and_diag
!!   2017/11/13 18:15 nakata
!!    Added subroutine get_weight_pDOS
!!   2018/09/19 18:30 nakata
!!    Added orbital angular momentum resolved DOS (pDOS_angmom)
!!   2018/10/22 14:24 dave & jsb
!!    Implementing (l,m)-projected DOS and moving flag into this routine
!!   2018/10/30 11:50 dave
!!    Added flag_pDOS_include_semicore and lines to exclude semi-core states from pDOS if required
!!***
module DiagModule

  use datatypes
  use global_module,          only: io_lun, area_DM, iprint_DM, flag_diagonalisation
  use GenComms,               only: cq_abort, inode, ionode, myid
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_matrices

  implicit none

  save

  !!****s* DiagModule/location *
  !!  NAME
  !!   location
  !!  PURPOSE
  !!   Holds the location of a matrix element in data_H.  When
  !!   converting between Conquest storage and Scalapack storage, we
  !!   need to know where in data_H a given row and column of the matrix
  !!   are, along with support function positions in that block.  This
  !!   is stored in loc.
  !!
  !!   Now holds the distance between the atoms (true distance, not FSC)
  !!   for general use with k-points (16/04/2002 dave)
  !!  AUTHOR
  !!   D.R.Bowler
  !!  SOURCE
  !!
  type location
     integer      :: supfn_row, supfn_col
     integer      :: loci
     integer      :: locj
     real(double) :: dx,dy,dz
  end type location
  !!***

  !!****s* DiagModule/elements *
  !!  NAME
  !!   elements
  !!  PURPOSE
  !!   Holds an array of locations of matrix elements.  When
  !!   redistributing Conquest matrix elements to Scalapack matrix
  !!   elements, this will be required when the support function radius
  !!   is more than quarter of the unit cell side (i.e. when an atom and
  !!   its periodic image are both neighbours of the same primary set
  !!   atom).  Allows more than one CQ element to contribute to a given
  !!   SC element
  !!  AUTHOR
  !!   D.R.Bowler
  !!  SOURCE
  !!
  type element
     integer                               :: n_elements
     type(location), dimension(:), pointer :: where
  end type element
  !!***

  !!****s* DiagModule/DistributeData *
  !!  NAME
  !!   DistributeData
  !!  PURPOSE
  !!   Holds arrays related to distribution of a matrix from Conquest
  !!   compressed row storage to ScaLAPACK 2D block cyclic storage
  !!  AUTHOR
  !!   D.R.Bowler
  !!  SOURCE
  !!
  type DistributeData
     integer,       dimension(:),     pointer :: num_rows, start_row
     type(element), dimension(:,:,:), pointer :: images
     integer,       dimension(:),     pointer :: send_rows, firstrow
  end type DistributeData
  !!***

  !!****s* DiagModule/Krecv_data *
  !!  NAME
  !!   Krecv_data
  !!  PURPOSE
  !!   Holds various data about matrix elements we're going to build
  !!   from eigenvectors we'll receive
  !!  ATTRIBUTES
  !!   integer orbs                  : Number of orbitals
  !!   integer, pointer ints(:)      : Number of interactions
  !!   integer, pointer ndimj(:)     : Dimension of j
  !!   integer, pointer prim_atom(:) : Primary set atom for a given interaction
  !!   integer, pointer locj(:)      : Where to put matrix element in data_Matrix
  !!   real(double), pointer dx, dy, dz : Vector between atoms
  !!  AUTHOR
  !!   D.R.Bowler
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  type Krecv_data
     integer                               :: orbs
     integer,      dimension(:),   pointer :: ints
     !integer,     dimension(:),   pointer :: ndimi
     integer,      dimension(:),   pointer :: ndimj
     integer,      dimension(:,:), pointer :: prim_atom
     integer,      dimension(:,:), pointer :: locj
     real(double), dimension(:,:), pointer :: dx, dy, dz
  end type Krecv_data
  !!***


  ! The matrix that holds the SC data - the HAMILTONIAN
  complex(double_cplx), dimension(:,:,:), allocatable :: SCHmat, SCSmat, z
  ! Buffer for receiving data
  complex(double_cplx), dimension(:,:),   allocatable :: RecvBuffer
  ! Buffer for sending data
  complex(double_cplx), dimension(:,:),   allocatable :: SendBuffer
  ! Data type that stores details of where to put and where to get elements
  type(DistributeData) :: DistribH, DistribS

  ! Fermi Energy
  real(double), dimension(2) :: Efermi
  integer, dimension(2) :: band_ef

  ! K-point data - here so that reading of k-points can take place in
  ! different routine to FindEvals
  integer :: nkp
  real(double), dimension(:),     allocatable :: wtk
  real(double), dimension(:,:),   allocatable :: kk
  real(double), dimension(:,:,:), allocatable :: occ
  ! 2007/08/13 dave changed this to be set by user
  real(double) :: kT
  logical :: first = .true.

  !logical :: diagon ! Do we diagonalise or use O(N) ?

  ! Local scratch data
  real(double), dimension(:,:,:), allocatable :: w ! matrix_size, nkp, nspin
  real(double), dimension(:,:),   allocatable :: local_w ! matrix_size, nspin
  !complex(double_cplx), dimension(:),allocatable :: work, rwork, gap
  complex(double_cplx), dimension(:), allocatable :: work
  real(double),         dimension(:), allocatable :: rwork, gap
  integer,              dimension(:), allocatable :: iwork, ifail, iclustr
  integer :: lwork, lrwork, liwork

  ! dimensions of local chunk of scalapack matrix
  integer :: row_size, col_size

  ! BLACS variables
  integer :: me
  integer :: context
  integer, dimension(50) :: desca, descz, descb
  real(double) :: abstol = 1e-30_double ! pdlamch (context,'U')

  ! Max number of iterations when searching for E_Fermi
  integer :: maxefermi

  ! Flags controlling Methfessel-Paxton approximation to step-function
  integer :: flag_smear_type, iMethfessel_Paxton

  ! Flags controlling the algorithms for finding Fermi Energy when
  ! using Methfessel-Paxton smearing
  real(double) :: gaussian_height, finess, NElec_less

  ! Maximum number of steps in the bracking search allowed before
  ! halfing incEf (introduced guarantee success in the very rare case
  ! that Methfessel-Paxton approximation may casue the bracket search
  ! algorithm to fail.)
  integer :: max_brkt_iterations

  ! DOS-related variables
  real(double), allocatable, dimension(:,:) :: total_DOS
  real(double), allocatable, dimension(:,:,:) :: pDOS
  real(double), allocatable, dimension(:,:,:,:,:) :: pDOS_angmom ! Bin, atom, l, m, spin
  real(double), allocatable, dimension(:,:) :: w_pDOS
  real(double) :: dE_DOS, pf_DOS
  integer :: n_DOS_max, n_DOS_wid
  logical :: flag_pDOS_include_semicore
  
contains

  ! -------------------------------------------------------------------
  ! Subroutine FindEvals
  ! -------------------------------------------------------------------

  !!****f* DiagModule/FindEvals *
  !!
  !!  NAME
  !!   FindEvals - finds the eigenvalues
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Call the ScalapackFormat routines and DiagModule routines to find
  !!   the eigenvalues by exact diagonalisation.  See the Conquest notes
  !!   "Implementation of Diagonalisation within Conquest" for a
  !!   detailed discussion of this routine.
  !!
  !!   Note added 31/07/2002 on Pulay force (drb)
  !!
  !!   The Pulay contribution to the force (which is also the rate of
  !!   change of the total energy with respect to a given support
  !!   function) requires the matrix M12 (which has range S).  When
  !!   using LNV order N this is given by 3LHL - 2(LHLSL + LSLHL), but
  !!   during diagonalisation it's given by:
  !!
  !!   M12_ij = -\sum_k w_k \sum_n f_n \epsilon_n c^n_i c^n_j
  !!
  !!   This can be trivially implemented using buildK if we scale the
  !!   occupancies by the eigenvalues after building K and before
  !!   building M12, and this is what is done.  Within buildK, the c_j
  !!   coefficients are scaled by f_n before a dot product is taken with
  !!   c_i, so scaling the occupancies by the eigenvalues
  !!   effectively builds M12.
  !!  INPUTS
  !!   real(double) :: electrons - number of electrons in system
  !!  USES
  !!   common, datatypes, GenComms, global_module, matrix_data,
  !!   maxima_module, primary_module, ScalapackFormat
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   04/03/2002
  !!  MODIFICATION HISTORY
  !!   11/03/2002 dave
  !!    Fixed calls to ScaLAPACK and increased scratch space used
  !!    Added fdf calls to get block sizes and processor grid
  !!   05/04/2002 dave (and before)
  !!    Main change was to encapsulate the Scalapack call only on processors
  !!    involved in diagonalisation.  Also generalised so that we work by
  !!    rows of the matrix (rather than by atoms)
  !!   16/04/2002 dave
  !!    Added reading of k-points (using fdf) and loop over k-points for
  !!    the diagonalisation calls.  Also now passes k-point to
  !!    DistributeCQtoSC and calls a new ScaLAPACK routine (pzheev) to
  !!    do the diagonalisation (as the Hamiltonian is now complex)
  !!   16/04/2002 dave
  !!    Moved all I/O and allocation for k-points down to readDiagInfo
  !!   18/04/2002 dave
  !!    Changed DistributeCQ_to_SC so that the matrix we're distributing
  !!    into is passed - this will allow distribution of S as well as H
  !!   22/04/2002 dave
  !!    Tidied and moved code around in preparation for building K from
  !!    eigenvectors
  !!   23/04/2002 dave
  !!    Added calls to find Fermi level and occupancies of k points, and
  !!    moved output of eigenvalues
  !!    so that occupancies can be written out as well
  !!   01/05/2002 dave
  !!    Moved allocation of SCHmat and z out of PrepareRecv
  !!   29/05/2002 dave
  !!    Moved initialisation (BLACS start up, descinit, DistribH,
  !!    allocation of memory) to initDiag
  !!   31/07/2002 dave
  !!    Added Pulay force calculation (see comments for discussion)
  !!   13:49, 24/01/2003 drb
  !!    Moved location of -M12 shift and tidied
  !!   2004/10/29 drb
  !!    Added check on size of Distrib before deallocating
  !!   2004/11/10 drb
  !!    Changed nsf to come from maxima, not common
  !!   09:10, 11/05/2005 dave
  !!    Added check on block sizes and matrix size
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new
  !!    matrix routines
  !!   2007/08/14 17:31 dave
  !!    Added entropy calculation
  !!   2010/02/13 L.Tong
  !!    Added k-point parallelisation
  !!   2011/06/14 16:41 dave
  !!    Small tweak to remove unnecessary gcopy calls and introduce
  !!    local group kpt scaling as variable
  !!   2011/12/01 L.Tong
  !!    Added spin polarisation
  !!    - added optional parameters electrons_up and electrons_dn
  !!      for correct electrons in the different spin channels
  !!    - added additional arrays expH_up and expH_dn (eignvector
  !!      matrices in Conquest format) for spin, in case of spin
  !!      polarised calculations the array expH is not allocated
  !!    - added SCS_dnmat for spin polarised calculations. In theory
  !!      from construction, matS and SCSmat should not be spin
  !!      dependent, but we need as copy of SCSmat for spin down
  !!      component because the original SCSmat will be changed as a
  !!      result of calling pzhegvx once. So another is needed for the
  !!      other spin component. This extra storage avoids the
  !!      communication intensive format conversion from matS to SCSmat
  !!      again.
  !!   2012/02/29 L.Tong
  !!   - uses findFermi interface for the different findFermi calls
  !!   2012/03/04 L.Tong
  !!   - Changed spin implementation
  !!   - electrons now an array of dimension nspin, which stores the
  !!     required constraint in each spin channel.
  !!   - Moved most memory deallocations to subroutine endDiag.
  !!   2012/06/17 L.Tong
  !!   - Found print_info was undefined when used in condition
  !!     (print_info == 0), so initialise it to 0 at beginning
  !!   2013/01/28 17:21 dave
  !!    Including changes for DeltaSCF from Umberto Terranova
  !!   2015/06/29 17:14 dave
  !!    Added DOS output
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf and paof
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2016/08/09 14:00 nakata
  !!    Renamed support_K -> atom_fns_K
  !!   2016/09/15 22:00 nakata
  !!    Introduce matBand -> matBand_atomf transformations
  !!   2017/06/21 17:19 dave
  !!    Tidying: remove block size check, calculation of matrix size
  !!    Moving calls to DistribH/S and pzhegvx into a new routine
  !!    Establishing desca/b/z as globals (calculated once)
  !!   2017/06/22 08:35 dave
  !!    Continued tidying: removing print_info
  !!   2017/06/29 14:00 nakata
  !!    Changed the position to deallocate matBand_atomf
  !!   2017/11/09 16:00 nakata
  !!    Introduced PDOS with MSSFs (tentative)
  !!   2017/11/13 18:15 nakata
  !!    Added optional normalization of each eigenstate of PDOS
  !!   2018/10/22 14:28 dave & jsb
  !!    Adding (l,m)-projection for pDOS
  !!   2018/11/13 17:30 nakata
  !!    Changed matS to be spin_SF dependent
  !!   2019/03/18 17:00 nakata
  !!    Added wf_self_con for accumulate_DOS
  !!   2019/05/09 dave
  !!    Bug fix for band density output (and added write_eigenvalues call)
  !!   2019/10/23 10:50 dave
  !!    Bug fix for pDOS output
  !!   2019/11/21 10:40 nakata
  !!    Bug fix for pDOS normalisation for spin-polarised calculations 
  !!    (added the dimension of spin to w_pDOS)
  !!  SOURCE
  !!
  subroutine FindEvals(electrons)

    use datatypes
    use numbers
    use units
    use global_module,   only: iprint_DM, ni_in_cell, numprocs,       &
         area_DM, flag_fix_spin_population,     &
         nspin, spin_factor, flag_DeltaSCF, flag_excite, &
         flag_local_excitation, dscf_HOMO_thresh, &
         dscf_LUMO_thresh, dscf_source_level, dscf_target_level, &
         dscf_target_spin, dscf_source_spin, flag_cdft_atom,  &
         dscf_HOMO_limit, dscf_LUMO_limit, &
         flag_out_wf,wf_self_con, max_wf, paof, sf, atomf, flag_out_wf_by_kp, &
         out_wf, n_DOS, E_DOS_max, E_DOS_min, flag_write_DOS, sigma_DOS, &
         flag_write_projected_DOS, flag_normalise_pDOS, flag_pDOS_angmom, flag_pDOS_lm, &
         E_wf_min, E_wf_max, flag_wf_range_Ef, &
         flag_SpinDependentSF
    use GenComms,        only: my_barrier, cq_abort, mtime, gsum, myid
    use ScalapackFormat, only: matrix_size, proc_rows, proc_cols,     &
         block_size_r,       &
         block_size_c, pg_kpoints, proc_groups, &
         nkpoints_max, pgid, N_procs_in_pg,     &
         N_kpoints_in_pg
    use mult_module,     only: matH, matS, matK, matM12, SF_to_AtomF_transform, &
         matrix_scale, matrix_product_trace, allocate_temp_matrix, free_temp_matrix
    use matrix_data,     only: Hrange, Srange, aHa_range
    use primary_module,  only: bundle
    use species_module,  only: species, nsf_species, species_label, n_species
    use memory_module,   only: type_dbl, type_int, type_cplx,         &
         reg_alloc_mem, reg_dealloc_mem
    use energy,          only: entropy
    use cdft_data, only: cDFT_NumberAtomGroups
    use maxima_module,   ONLY: maxngrid
    use functions_on_grid,           only: atomfns, &
         allocate_temp_fn_on_grid,    &
         free_temp_fn_on_grid
    use density_module, ONLY: get_band_density
    use io_module, ONLY: dump_DOS, dump_projected_DOS, write_eigenvalues
    use pao_format, ONLY: pao

    implicit none

    ! Passed variables
    real(double), dimension(:), intent(in) :: electrons

    ! Local variables
    real(double)                   :: a, time0, time1, vl, vu, &
         orfac, scale, entropy_total,     &
         bandE_total, coeff, setA, setB, Eband
    real(double), dimension(nspin) :: locc, bandE, entropy_local
    real(double), external         :: dlamch
    complex(double_cplx), dimension(:,:,:), allocatable :: expH
    complex(double_cplx) :: c_n_alpha2, c_n_setA2, c_n_setB2
    integer :: info, stat, il, iu, i, j, m, mz, prim_size, ng, wf_no, &
         kp, spin, spin_SF, iacc, iprim, l, band, cdft_group, atom_fns_K, &
         n_band_min, n_band_max
    integer, allocatable, dimension(:) :: matBand
    integer, allocatable, dimension(:,:) :: matBand_kp
    integer :: matBand_atomf
    integer :: iatom_spec, Nangmom ! number of orbital angular momentum to be dumped (ex. (s,p,d)=3)

    logical :: flag_keepexcite, flag_full_DOS

    real(double), dimension(:),allocatable :: abs_wf

    if (iprint_DM >= 2 .AND. myid == 0) &
         write (io_lun, fmt='(10x,"Entering FindEvals")')
    prim_size = 0
    do i = 1, bundle%n_prim
       prim_size = prim_size + nsf_species(bundle%species(i))
    end do

    ! Initialise - start BLACS, sort out matrices, allocate memory
    call initDiag

    scale = one / real(N_procs_in_pg(pgid), double)

    ! ------------------------------------------------------------------------
    ! Start diagonalisation
    ! ------------------------------------------------------------------------
    ! First diagonalisation - get eigenvalues only (so that we can find Efermi)
    time0 = mtime ()

    if (iprint_DM >= 2 .and. (inode == ionode)) &
         write (io_lun, fmt='(10x,"In FindEvals, tolerance is ", g20.12)') &
         abstol

    ! zero the global and local eigenvalues
    w = zero
    local_w = zero

    do spin = 1, nspin
       do i = 1, nkpoints_max ! Loop over the kpoints within each process group
          call distrib_and_diag(spin,i,'N',.true.) ! Eigenvalues
          !call my_barrier()
       end do ! End do i = 1, nkpoints_max
       ! sum the w on each node together to give the whole w on each
       ! node, note that the repeating of the same eigenvalues in each
       ! proc_group is taken care of by the additional factor
       ! 1 / N_procs_in_pg
       call gsum(w(:,:,spin), matrix_size, nkp)
    end do ! spin
    ! Allocate matrices to store band K matrices
    time1 = mtime()
    if (iprint_DM >= 2 .AND. myid == 0) &
         write (io_lun, 2) myid, time1 - time0
    ! If we are trying to localise a level, we do NOT want to excite just yet
    if(flag_DeltaSCF.AND.flag_excite.AND.flag_local_excitation) then
       flag_keepexcite = .true.
       flag_excite = .false.
    end if
    ! Find Fermi level, given the eigenvalues at all k-points (in w)
    ! if (me < proc_rows*proc_cols) then
    call findFermi(electrons, w, matrix_size, nkp, Efermi, occ)
    if (wf_self_con .and. flag_out_wf) then
       ! Has the user specified an energy range ? If so, work out the band limits
       if(max_wf==0) then
          n_band_min = 1e9
          n_band_max = 0
          do spin=1,nspin
             do i=1,nkpoints_max
                if (i <= N_kpoints_in_pg(pgid)) then
                   kp = pg_kpoints(pgid, i)
                   do j = 1, matrix_size
                      if(flag_wf_range_Ef) then
                         Eband = w(j,kp,spin) - Efermi(spin)
                      else
                         Eband = w(j,kp,spin)
                      end if
                      if((Eband>E_wf_min).AND.(j<n_band_min)) n_band_min = j
                      if((Eband<E_wf_max).AND.(j>n_band_max)) n_band_max = j
                   end do
                end if
             end do
          end do
          if(myid==0.AND.iprint_DM>=2) write(io_lun,fmt='(2x,"WF band limits set to: ",2i6)') &
               n_band_min, n_band_max
          max_wf = n_band_max - n_band_min + 1
          allocate(out_wf(max_wf))
          do i=1,max_wf
             out_wf(i) = n_band_min + i-1
          end do
       end if
       if(flag_out_wf_by_kp) then
          allocate(matBand_kp(max_wf,nkp))
          do j=1,nkp
             do i=1,max_wf
                matBand_kp(i,j) = allocate_temp_matrix(Hrange,0,sf,sf)
             end do
          end do
       else
          allocate(matBand(max_wf))
          do i=1,max_wf
             matBand(i) = allocate_temp_matrix(Hrange,0,sf,sf)
          end do
       end if
       if (atomf.ne.sf) matBand_atomf = allocate_temp_matrix(aHa_range,0,atomf,atomf)
    end if
    ! Preparatory work for DOS
    if(wf_self_con.AND.flag_write_DOS) then
       allocate(total_DOS(n_DOS,nspin))
       total_DOS = zero
       ! Only if projecting DOS onto atoms
       if(flag_write_projected_DOS) then
          if (atomf==sf) allocate(pDOS(n_DOS,bundle%n_prim,nspin))
          if (atomf/=sf) allocate(pDOS(n_DOS,ni_in_cell,   nspin))
          pDOS = zero
          if (flag_pDOS_angmom) then
             Nangmom = 0
             do iatom_spec = 1, n_species
                Nangmom = max(Nangmom, pao(iatom_spec)%greatest_angmom)
             enddo
             Nangmom = Nangmom + 1 ! Nangmom should be 1 larger than greatest_angmom
             if(flag_pDOS_lm) then
                ! NB we want (2*l+1) for m component but Nangmom = l+1 so 2*l+1 = 2*Nangmom-1
                if (atomf==sf) allocate(pDOS_angmom(n_DOS,bundle%n_prim,Nangmom,2*Nangmom-1,nspin)) 
                if (atomf/=sf) allocate(pDOS_angmom(n_DOS,ni_in_cell,   Nangmom,2*Nangmom-1,nspin))
             else
                if (atomf==sf) allocate(pDOS_angmom(n_DOS,bundle%n_prim,Nangmom,1,nspin)) ! 1 for m component (not used)
                if (atomf/=sf) allocate(pDOS_angmom(n_DOS,ni_in_cell,   Nangmom,1,nspin))
             end if
             pDOS_angmom = zero
          endif
          if (flag_normalise_pDOS) then
             allocate(w_pDOS(matrix_size,nspin))
             w_pDOS = zero
          end if
       end if
       ! If the user hasn't specified limits
       if(E_DOS_min==zero.AND.E_DOS_max==zero) then
          E_DOS_min = 1e30_double
          E_DOS_max = -1e30_double
          do i=1,nkp
             if(E_DOS_min>w(1,i,1)) E_DOS_min = w(1,i,1)
             if(E_DOS_max<w(matrix_size,i,1)) E_DOS_max = w(matrix_size,i,1)
          end do
          if(myid==0.AND.iprint_DM>=2) write(io_lun,fmt='(2x,"DOS limits set automatically: ",2f12.5)') &
               E_DOS_min, E_DOS_max
       end if
       dE_DOS = (E_DOS_max - E_DOS_min)/real(n_DOS-1,double)
       n_DOS_wid = floor(6.0_double*sigma_DOS/dE_DOS) ! How many bins either side of state we consider
       pf_DOS = one/(sigma_DOS*sqrt(twopi))
    end if

    ! Allocate space to expand eigenvectors into (i.e. when reversing
    ! ScaLAPACK distribution)
    allocate(expH(matrix_size,prim_size,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort('FindEvals: failed to alloc expH', stat)
    call reg_alloc_mem(area_DM, matrix_size * prim_size * nspin, type_cplx)

    !call gcopy(Efermi)
    !call gcopy(occ,matrix_size,nkp)
    ! else
    !   call gcopy(Efermi)
    !   call gcopy(occ,matrix_size,nkp)
    ! end if
    if(flag_DeltaSCF.AND.flag_local_excitation.AND.flag_keepexcite) then
       ! Set band_ef: this works because with no spin a factor of two is normally applied
       if(nspin>1) then
          do spin=1,nspin
             band_ef(spin) = int(electrons(spin))
          end do
       else
          band_ef(:) = int(electrons(1))
       end if
       flag_excite = .true.
       if(dscf_homo_limit/=0) then
          ! Diagonalise
          spin = dscf_source_spin
          i = 1 ! Just use first k-point - local states should not show dispersion
          call distrib_and_diag(spin,i,'V',.false.)
          call my_barrier()
          ! Reverse the CQ to SC distribution so that eigenvector
          ! coefficients for atoms are on the appropriate processor.
          do ng = 1, proc_groups
             if (i <= N_kpoints_in_pg(ng)) then
                kp = pg_kpoints(ng, i)
                call DistributeSC_to_ref(DistribH, ng, z(:,:,spin), &
                     expH(:,:,spin))       ! Find local excitation
             end if
          end do
          ! Locate state
          dscf_source_level = band_ef(spin)
          do band=band_ef(spin),band_ef(spin)-dscf_homo_limit,-1
             setA = zero
             setB = zero
             iacc = 1
             do iprim=1,bundle%n_prim
                coeff = zero
                if(flag_cdft_atom(bundle%ig_prim(iprim))==1) then ! HOMO group
                   do l=1,nsf_species(bundle%species(iprim))
                      coeff = coeff+expH(band,iacc+l-1,spin)*conjg(expH(band,iacc+l-1,spin))
                   end do
                   setA = setA+coeff
                else
                   do l=1,nsf_species(bundle%species(iprim))
                      coeff = coeff+expH(band,iacc+l-1,spin)*conjg(expH(band,iacc+l-1,spin))
                   end do
                   setB = setB+coeff
                end if
                iacc = iacc + nsf_species(bundle%species(iprim))
             end do
             call gsum(setA)
             call gsum(setB)
             if(inode==ionode.AND.iprint_DM>2) &
                  write(io_lun,fmt='(4x,"DeltaSCF HOMO search: Band, coefficients A/B: ",i5,2f12.5)') band,setA,setB
             if(setA>dscf_HOMO_thresh) then
                dscf_source_level = band
                exit
             end if
          end do
       end if
       if(dscf_source_spin/=dscf_target_spin.OR.dscf_homo_limit==0) then ! Need to re-diagonalise
          ! Diagonalise
          spin = dscf_target_spin
          i = 1 ! Just use first k-point - local states should not show dispersion
          call distrib_and_diag(spin,i,'V',.false.)
          ! Reverse the CQ to SC distribution so that eigenvector
          ! coefficients for atoms are on the appropriate processor.
          do ng = 1, proc_groups
             if (i <= N_kpoints_in_pg(ng)) then
                kp = pg_kpoints(ng, i)
                call DistributeSC_to_ref(DistribH, ng, z(:,:,spin), &
                     expH(:,:,spin))       ! Find local excitation
             end if
          end do
       end if
       ! Locate LUMO
       dscf_target_level = band_ef(spin)+1
       if(cDFT_NumberAtomGroups==1) then
          cdft_group = 1
       else if(cDFT_NumberAtomGroups==2) then
          cdft_group = 2
       else
          call cq_abort("You must specify a group of atoms using cDFT.NumberAtomsGroup")
       end if
       do band=band_ef(spin)+1,band_ef(spin)+dscf_lumo_limit
          setA = zero
          setB = zero
          iacc = 1
          do iprim=1,bundle%n_prim
             coeff = zero
             if(flag_cdft_atom(bundle%ig_prim(iprim))==cdft_group) then ! HOMO group
                do l=1,nsf_species(bundle%species(iprim))
                   coeff = coeff+expH(band,iacc+l-1,spin)*conjg(expH(band,iacc+l-1,spin))
                end do
                setA = setA+coeff
             else
                do l=1,nsf_species(bundle%species(iprim))
                   coeff = coeff+expH(band,iacc+l-1,spin)*conjg(expH(band,iacc+l-1,spin))
                end do
                setB = setB+coeff
             end if
             iacc = iacc + nsf_species(bundle%species(iprim))
          end do
          call gsum(setA)
          call gsum(setB)
          if(inode==ionode.AND.iprint_DM>2) &
               write(io_lun,fmt='(4x,"DeltaSCF LUMO search: Band, coefficients A/B: ",i5,2f12.5)') band,setA,setB
          if(setA>dscf_LUMO_thresh) then
             dscf_target_level = band
             exit
          end if
       end do
       ! Find Fermi level, given the eigenvalues at all k-points (in w)
       call findFermi(electrons, w, matrix_size, nkp, Efermi, occ)
    end if ! DeltaSCF localised excitation
    ! Now write out eigenvalues and occupancies
    if (iprint_DM == 2 .AND. myid == 0) then
       bandE = zero
       do i = 1, nkp
          write (io_lun, 7) i, kk(1,i), kk(2,i), kk(3,i)
          do spin = 1, nspin
             if (nspin == 2) &
                  write (io_lun, '(10x,"For spin = ",i1)') spin
             do j = 1, matrix_size, 3
                if (j == matrix_size) then
                   write (io_lun, 8) w(j,i,spin), occ(j,i,spin)
                   bandE(spin) = bandE(spin) + w(j,i,spin) * occ(j,i,spin)
                else if (j == matrix_size - 1) then
                   write (io_lun, 9) w(j,i,spin), occ(j,i,spin), &
                        w(j+1,i,spin), occ(j+1,i,spin)
                   bandE(spin) = bandE(spin) + w(j,i,spin) * occ(j,i,spin) + &
                        w(j+1,i,spin) * occ(j+1,i,spin)
                else
                   write (io_lun, 10) w(j,i,spin), occ(j,i,spin), &
                        w(j+1,i,spin), occ(j+1,i,spin), &
                        w(j+2,i,spin), occ(j+2,i,spin)
                   bandE(spin) = bandE(spin) + w(j,i,spin) * occ(j,i,spin) + &
                        w(j+1,i,spin) * occ(j+1,i,spin) + &
                        w(j+2,i,spin) * occ(j+2,i,spin)
                endif
             end do ! j=matrix_size
             write (io_lun, &
                  fmt='("Sum of eigenvalues for spin = ", &
                  &i1, ": ", f18.11," ", a2)') &
                  spin, en_conv * bandE(spin), en_units(energy_units)
          end do ! spin
          if (nspin == 2) then
             write (io_lun, &
                  fmt='("Total sum of eigenvalues: ", f18.11, " ",a2)') &
                  en_conv * (bandE(1) + bandE(2)), en_units(energy_units)
          else
             write(io_lun, 4) en_conv * two * bandE(1), en_units(energy_units)
          end if
       end do ! do i = 1, nkp
    else if (iprint_DM >= 3 .AND. myid == 0) then
       bandE = zero
       do i = 1, nkp
          write (io_lun, 7) i, kk(1,i), kk(2,i), kk(3,i)
          do spin = 1, nspin
             if (nspin == 2) &
                  write (io_lun, '(10x,"For spin = ",i1)') spin
             do j = 1, matrix_size
                write (io_lun, fmt='(10x,i5,f12.5,f6.3)') j, w(j,i,spin), occ(j,i,spin)
                bandE(spin) = bandE(spin) + w(j,i,spin) * occ(j,i,spin)
             end do ! j=matrix_size
             write (io_lun, &
                  fmt='("Sum of eigenvalues for spin = ", &
                  &i1, ": ", f18.11," ", a2)') &
                  spin, en_conv * bandE(spin), en_units(energy_units)
          end do ! spin
          if (nspin == 2) then
             write (io_lun, &
                  fmt='("Total sum of eigenvalues: ", f18.11, " ",a2)') &
                  en_conv * (bandE(1) + bandE(2)), en_units(energy_units)
          else
             write(io_lun, 4) en_conv * two * bandE(1), en_units(energy_units)
          end if
       end do ! do i = 1, nkp
    end if ! if(iprint_DM>=1.AND.myid==0)

    time0 = mtime()
    do spin = 1, nspin
       call matrix_scale(zero, matK(spin))
       call matrix_scale(zero, matM12(spin))
    end do
    ! Second diagonalisation - get eigenvectors and build K
    entropy = zero
    do spin = 1, nspin
       do i = 1, nkpoints_max
          call distrib_and_diag(spin,i,'V',.false.)
          if (iprint_DM >= 5 .and. inode == ionode) &
               write (io_lun, *) myid, ' Calling barrier'
          call my_barrier()
          if(iprint_DM >= 5 .AND. inode == ionode) &
               write (io_lun, *) myid, ' Calling DistributeSC_to_Ref'
          ! Reverse the CQ to SC distribution so that eigenvector
          ! coefficients for atoms are on the appropriate processor.
          ! Loop over the process-node groups, we build K one k-point at
          ! a time
          do ng = 1, proc_groups
             if (i <= N_kpoints_in_pg(ng)) then
                kp = pg_kpoints(ng, i)
                call DistributeSC_to_ref(DistribH, ng, z(:,:,spin), &
                     expH(:,:,spin))
                if(wf_self_con.AND.flag_write_DOS) then
                   if(flag_write_projected_DOS) then
                      if (flag_normalise_pDOS) then
                         call get_weight_pDOS(expH(:,:,spin),w_pDOS(:,spin))
                         call gsum(w_PDOS(:,spin), matrix_size)
                      endif
                      if (flag_pDOS_angmom) then
                         if (flag_normalise_pDOS) then
                            call accumulate_DOS(wtk(kp),w(:,kp,spin),expH(:,:,spin),total_DOS(:,spin),spin, &
                                                projDOS=pDOS(:,:,spin),projDOS_angmom=pDOS_angmom(:,:,:,:,spin), &
                                                weight_pDOS=w_pDOS(:,spin))
                         else
                            call accumulate_DOS(wtk(kp),w(:,kp,spin),expH(:,:,spin),total_DOS(:,spin),spin, &
                                                projDOS=pDOS(:,:,spin),projDOS_angmom=pDOS_angmom(:,:,:,:,spin))
                         endif
                      else
                         if (flag_normalise_pDOS) then
                            call accumulate_DOS(wtk(kp),w(:,kp,spin),expH(:,:,spin),total_DOS(:,spin),spin, &
                                                projDOS=pDOS(:,:,spin),weight_pDOS=w_pDOS(:,spin))
                         else
                            call accumulate_DOS(wtk(kp),w(:,kp,spin),expH(:,:,spin),total_DOS(:,spin),spin, &
                                                projDOS=pDOS(:,:,spin))
                         endif
                      endif
                   else
                      call accumulate_DOS(wtk(kp),w(:,kp,spin),expH(:,:,spin),total_DOS(:,spin),spin)
                   end if
                end if
                ! Build K and K_dn from the eigenvectors
                if (iprint_DM >= 4 .and. inode == ionode) &
                     write (io_lun, *) myid, ' Calling buildK ', &
                     Hrange, matK(spin)
                ! Pass band-by-band K matrices if we are outputting densities
                if(wf_self_con .and. flag_out_wf) then
                   if(flag_out_wf_by_kp) then
                      call buildK(Hrange, matK(spin), occ(:,kp,spin), &
                           kk(:,kp), wtk(kp), expH(:,:,spin),matBand_kp(:,kp))
                   else
                      call buildK(Hrange, matK(spin), occ(:,kp,spin), &
                           kk(:,kp), wtk(kp), expH(:,:,spin),matBand)
                   end if
                else
                   call buildK(Hrange, matK(spin), occ(:,kp,spin), &
                        kk(:,kp), wtk(kp), expH(:,:,spin))
                end if
                ! Build matrix needed for Pulay force
                ! We scale the occupation number for this k-point by the
                ! eigenvalues in order to build the matrix M12
                ! We can do this simply because we won't use them again
                ! (though we could use a dummy variable if we wanted to
                ! use them again)
                do j = 1, matrix_size
                   ! Calculate entropic contribution to electronic energy
                   select case (flag_smear_type)
                   case (0) ! Fermi smearing
                      if (occ(j,kp,spin) > RD_ERR .and. &
                           (wtk(kp) - occ(j,kp,spin)) > RD_ERR) then
                         locc(spin) = occ(j,kp,spin) / wtk(kp)
                         entropy_local(spin) = &
                              locc(spin) * log(locc(spin)) + &
                              (one - locc(spin)) * log (one - locc(spin))
                         if (iprint_DM > 3 .and. inode == ionode) &
                              write (io_lun, &
                              fmt='(2x,"Spin, Occ, wt: ", i1, 2f12.8, &
                              &" ent: ", f20.12)') &
                              spin, locc(spin), wtk(kp), entropy_local(spin)
                         entropy = entropy - &
                              spin_factor * wtk(kp) * entropy_local(spin)
                      end if
                   case (1) ! Methfessel-Paxton smearing
                      entropy = entropy + spin_factor * wtk(kp) * &
                           MP_entropy((w(j,kp,spin) - Efermi(spin)) / kT, &
                           iMethfessel_Paxton)
                   case default
                      call cq_abort ("FindEvals: Smearing flag not recognised",&
                           flag_smear_type)

                   end select
                   ! occ is now used to construct matM12, factor by eps^n
                   ! to allow reuse of buildK
                   occ(j,kp,spin) = - occ(j,kp,spin) * w(j,kp,spin)
                end do ! j = 1, matrix_size
                ! Now build data_M12_ij (=-\sum_n eps^n c^n_i c^n_j -
                ! hence scaling occs by eps allows reuse of buildK)
                call buildK(Srange, matM12(spin), occ(:,kp,spin), &
                     kk(:,kp), wtk(kp), expH(:,:,spin))
             end if ! End if (i <= N_kpoints_in_pg(ng)) then
          end do ! End do ng = 1, proc_groups
       end do ! End do i = 1, nkpoints_max
    end do ! spin

    !------ output WFs  --------
    if (wf_self_con .and. flag_out_wf) then
       allocate(abs_wf(maxngrid),STAT=stat)
       if (stat /= 0) call cq_abort('wf_out: Failed to allocate wfs', stat)
       call reg_alloc_mem(area_DM, maxngrid, type_dbl)
       atom_fns_K = allocate_temp_fn_on_grid(atomf)
       if(flag_out_wf_by_kp) then
          if(inode==ionode) call write_eigenvalues(w,matrix_size,nkp,nspin,kk,wtk,Efermi)
          do i=1,nkp
             do wf_no=1,max_wf
                write(io_lun,fmt='(2x,"Band : ",i4," k-point: ",i3)') out_wf(wf_no),i
                if(nspin>1) then
                   do spin = 1, nspin
                      abs_wf(:)=zero
                      if (atomf.ne.sf) then
                         call SF_to_AtomF_transform(matBand_kp(wf_no,i), matBand_atomf, spin, Hrange)
                         call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand_atomf,maxngrid)
                      else                     
                         call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand_kp(wf_no,i),maxngrid)
                      endif
                      if(i==1) then
                         call wf_output(spin,abs_wf,wf_no,kk(:,i),w(out_wf(wf_no),i,spin),i)
                      else
                         call wf_output(spin,abs_wf,wf_no,kk(:,i),w(out_wf(wf_no),i,spin),i)
                      end if
                      call my_barrier()
                   end do
                else
                   abs_wf(:)=zero
                   spin = 1
                   if (atomf.ne.sf) then
                      call SF_to_AtomF_transform(matBand_kp(wf_no,i), matBand_atomf, spin, Hrange)
                      call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand_atomf,maxngrid)
                   else
                      call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand_kp(wf_no,i),maxngrid)
                   endif
                   if(i==1) then
                      call wf_output(0,abs_wf,wf_no,kk(:,i),w(out_wf(wf_no),i,spin),i)
                   else
                      call wf_output(0,abs_wf,wf_no,kk(:,i),w(out_wf(wf_no),i,spin))
                   end if
                   call my_barrier()
                end if
             end do
          end do
       else
          if(inode==ionode) call write_eigenvalues(w,matrix_size,nkp,nspin,kk,wtk,Efermi)
          if(nspin>1) then
             do wf_no=1,max_wf
                do spin = 1, nspin
                   abs_wf(:)=zero
                   if (atomf.ne.sf) then
                      call SF_to_AtomF_transform(matBand(wf_no), matBand_atomf, spin, Hrange)
                      call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand_atomf,maxngrid)
                   else
                      call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand(wf_no),maxngrid)
                   endif
                   call wf_output(spin,abs_wf,wf_no)
                   call my_barrier()
                end do
             end do
          else
             spin = 1
             do wf_no=1,max_wf
                abs_wf(:)=zero
                if (atomf.ne.sf) then
                   call SF_to_AtomF_transform(matBand(wf_no), matBand_atomf, spin, Hrange)
                   call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand_atomf,maxngrid)
                else
                   call get_band_density(abs_wf,spin,atomfns,atom_fns_K,matBand(wf_no),maxngrid)
                endif
                call wf_output(0,abs_wf,wf_no)
                call my_barrier()
             end do
          end if
       end if
       deallocate(abs_wf,STAT=stat)
       if (stat /= 0) call cq_abort('Find Evals: Failed to deallocate wfs',stat)
       call reg_dealloc_mem(area_DM, maxngrid, type_dbl)
       call free_temp_fn_on_grid(atom_fns_K)
       if (atomf.ne.sf) call free_temp_matrix(matBand_atomf)
       if(flag_out_wf_by_kp) then
          do j=nkp,1,-1
             do i=max_wf,1,-1
                call free_temp_matrix(matBand_kp(i,j))
             end do
          end do
       else
          do i=max_wf,1,-1
             call free_temp_matrix(matBand(i))
          end do
       end if
    end if
    if(wf_self_con.AND.flag_write_DOS) then
       ! output DOS
       if(inode==ionode) call dump_DOS(total_DOS,Efermi)
       call my_barrier()
       if(flag_write_projected_DOS) then
          if (atomf==sf) then
             if (.not.flag_pDOS_angmom) call dump_projected_DOS(pDOS,Efermi)
             if (     flag_pDOS_angmom) call dump_projected_DOS(pDOS,Efermi,pDOS_angmom=pDOS_angmom,Nangmom=Nangmom)
          else
             call gsum(pDOS(:,:,:),n_DOS,ni_in_cell,nspin)
             if(flag_pDOS_angmom) then
                if(flag_pDOS_lm) then
                   do spin = 1, nspin
                      call gsum(pDOS_angmom(:,:,:,:,spin),n_DOS,ni_in_cell,Nangmom,2*Nangmom-1)
                   enddo
                else
                   do spin = 1, nspin
                      call gsum(pDOS_angmom(:,:,:,1,spin),n_DOS,ni_in_cell,Nangmom)
                   enddo
                end if
                if (inode==ionode) then
                   call dump_projected_DOS(pDOS,Efermi,pDOS_angmom=pDOS_angmom,Nangmom=Nangmom)
                end if
             else
                if (inode==ionode) call dump_projected_DOS(pDOS,Efermi)
             endif
          endif
          if (flag_normalise_pDOS) deallocate(w_pDOS)
          if (flag_pDOS_angmom) deallocate(pDOS_angmom)
          deallocate(pDOS)
       end if
       deallocate(total_DOS)
    end if

    if (iprint_DM > 3 .and. inode == ionode) &
         write (io_lun, *) "Entropy, TS: ", entropy, kT * entropy
    ! store entropy as TS instead of S
    entropy = entropy * kT
    time1 = mtime()
    if (iprint_DM >= 2 .and. inode == ionode) &
         write (io_lun, 3) myid, time1 - time0

    ! -------------------------------------------------------------
    ! End diagonalisation
    ! -------------------------------------------------------------
    ! Write out the Fermi Energy
    if (iprint_DM >= 1 .and. inode == ionode) then
       do spin = 1, nspin
          write (io_lun, 13) spin, en_conv * Efermi(spin), &
               en_units(energy_units)
       end do
    end if
    ! Write out the band energy and trace of K
    if (iprint_DM >= 1) then
       ! for tr(K.H)
       bandE_total = zero
       do spin = 1, nspin
          bandE(spin) = matrix_product_trace(matK(spin), matH(spin))
          bandE_total = bandE_total + spin_factor * bandE(spin)
       end do
       if (inode == ionode) then
          if (nspin == 1) then
             write (io_lun, 5)  en_conv * bandE_total, en_units(energy_units)
          else
             write (io_lun, 14) en_conv * bandE(1), en_units(energy_units), &
                  en_conv * bandE(2), en_units(energy_units), &
                  en_conv * bandE_total, en_units(energy_units)
          end if
       end if
       ! for tr(S.G)
       bandE_total = zero
       spin_SF = 1
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          bandE(spin) = matrix_product_trace(matS(spin_SF), matM12(spin))
          bandE_total = bandE_total + spin_factor * bandE(spin)
       end do
       if (inode == ionode) then
          if (nspin == 1) then
             write (io_lun, 6) en_conv * bandE_total, en_units(energy_units)
          else
             write (io_lun, 15) en_conv * bandE(1), en_units(energy_units), &
                  en_conv * bandE(2), en_units(energy_units), &
                  en_conv * bandE_total, en_units(energy_units)
          end if
       end if
    end if ! iprint_DM >= 1

    call my_barrier()

    ! deallocate
    ! local
    deallocate(expH, STAT=stat)
    if (stat /= 0) call cq_abort('FindEvals: failed to deallocate expH', stat)
    call reg_dealloc_mem(area_DM, matrix_size * prim_size * nspin, type_cplx)
    ! global
    call endDiag

    return

2   format(10x,'Proc: ',i5, ' Time taken for eval diag: ',f20.8,' ms')
3   format(10x,'Proc: ', i5,' Time taken for evec diag: ',f20.8,' ms')
4   format(10x,'Sum of eigenvalues: ',f18.11,' ',a2)
5   format(10x,'Energy as 2Tr[K.H]: ',f18.11, ' ',a2)
6   format(10x,'2Tr[S.G]: ',f18.11,' ',a2)
7   format(10x,'Eigenvalues and occupancies for k-point ',i3,' : ',3f12.5)
8   format(10x,f12.5,f6.3,2x)
9   format(10x,f12.5,f6.3,2x,f12.5,f6.3,2x)
10  format(10x,f12.5,f6.3,2x,f12.5,f6.3,2x,f12.5,f6.3,2x)
11  format(10x,'Proc: ',i5,' Diagonalising for eigenvectors')
12  format(10x,'Proc: ',i5,' row, col size: ',2i5)
13  format(10x,'Fermi energy for spin = ',i1,' is ',f18.11,' ',a2)
14  format(10x,'Energies Tr[K.H_up], Tr[K.H_down], and their sum: ',&
         f18.11,' ',a2,/,60x,f18.11,' ',a2,/,60x,f18.11,' ',a2)
15  format(19x,'Tr[K.G_up], Tr[K.G_down], and their sum: ',&
         f18.11,' ',a2,/,60x,f18.11,' ',a2,/,60x,f18.11,' ',a2)
  end subroutine FindEvals
  !!***


  ! -----------------------------------------------------------------------------
  ! Subroutine initDiag
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/initDiag *
  !!
  !!  NAME
  !!   initDiag
  !!  USAGE
  !!   call initDiag (desca, descb, descz)
  !!  PURPOSE
  !!   Contains various routines and assignments that initialise the
  !!   diagonalisation
  !!  INPUTS
  !!   integer, dimension(50) :: desca, descb, descz - descriptors for
  !!                             H, S and eigenvectors
  !!  OUTPUT
  !!   All variables below are ScaLAPACK scratch space and scratch sizes
  !!   (found and allocated here)
  !!
  !!   complex(double_cplx), dimension(:),allocatable :: work, rwork, gap
  !!   integer, dimension(:), allocatable :: iwork, ifail, iclustr
  !!  USES
  !!   ScalapackFormat, matrix_data, GenComms
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   29/05/2002
  !!  MODIFICATION HISTORY
  !!   2011/12/01 L.Tong
  !!   - Updated call to PrepareSend, removed unnecessay passing of matH
  !!     and matS, not used in PrepareSend
  !!   - Added spin polarisation
  !!   - Added the usage of numbers
  !!   - Added memory registers
  !!   2012/03/08 L.Tong
  !!   - Made lwork, lrwork and liwork module globals
  !!   - Removed redundant dependency on matH and matS
  !!   - Removed redundant dependency on mat from matrix_data
  !!   2015/10.06 dave
  !!   - Changed size of rwo to three (as required by pzhegvx workspace query)
  !!   2017/06/22 dave
  !!   Made descriptors module variables
  !!   Moved many operations out to new routines so this contains only work needed each time
  !!  SOURCE
  !!
  subroutine initDiag

    use numbers
    use ScalapackFormat, only: proc_rows, proc_cols, matrix_size
    use global_module,   only: numprocs, nspin
    use GenComms,        only: my_barrier, cq_abort, myid
    use memory_module,   only: type_dbl, type_int, type_cplx,         &
         reg_alloc_mem, reg_dealloc_mem

    implicit none

    ! Local variables
    integer         :: nump, merow, mecol, numrows, numcols, info, &
         stat, m, mz, ng
    real(double)    :: rwo(3)
    complex(double) :: wo(1)
    integer         :: iwo(1)

    ! First, work out how much data we're going to receive How many
    ! rows and columns do we have ? only works if the rows are exact
    ! integer multiples of block rows

    ! Allocate space for the distributed Scalapack matrices
    stat = 0
    allocate(SCHmat(row_size,col_size,nspin), SCSmat(row_size,col_size,nspin), &
         z(row_size,col_size,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("initDiag: failed to allocate SCHmat, SCSmat and z", stat)
    call reg_alloc_mem(area_DM, 3 * row_size * col_size * nspin, type_cplx)
    SCHmat = zero
    SCSmat = zero
    z = zero

    allocate(w(matrix_size,nkp,nspin), occ(matrix_size,nkp,nspin), STAT=stat)
    if (stat /= 0) call cq_abort('initDiag: failed to allocate w and occ', stat)
    call reg_alloc_mem(area_DM, 2 * matrix_size * nkp * nspin, type_dbl)

    allocate(local_w(matrix_size, nspin), STAT=stat)
    if (stat /= 0) call cq_abort('initDiag: failed to allocate local_w', stat)
    call reg_alloc_mem(area_DM, matrix_size * nspin, type_dbl)

    allocate(ifail(matrix_size), iclustr(2 * proc_rows * proc_cols), STAT=stat)
    if (stat /= 0) call cq_abort("initDiag: failed to allocate ifail and iclustr", stat)
    call reg_alloc_mem(area_DM, matrix_size + 2 * proc_rows * proc_cols, type_int)

    allocate(gap(proc_rows * proc_cols), STAT=stat)
    if (stat /= 0) call cq_abort("initDiag: failed to allocate gap", stat)
    call reg_alloc_mem(area_DM, proc_rows * proc_cols, type_dbl)

    ! the pzhegvx is only called here to get the optimal work array
    call pzhegvx(1, 'V', 'A', 'U', matrix_size, SCHmat(:,:,1), 1, 1,  &
         desca, SCSmat(:,:,1), 1, 1, descb, zero, zero, 0, 0, &
         1.0e-307_double, m, mz, w(1,1,1), -one, z(:,:,1), 1, &
         1, descz, wo, -1, rwo, -1, iwo, -1, ifail, iclustr,  &
         gap, info)

    ! Allocate scratch space for ScaLAPACK
    lwork  = 2 * wo(1)
    lrwork = 2 * rwo(1)
    liwork = iwo(1)
    allocate(work(lwork), rwork(lrwork), STAT=stat)
    if (stat /= 0) call cq_abort ('initDiag: failed to allocate work and rwork', stat)
    call reg_alloc_mem(area_DM, lwork + lrwork, type_cplx)

    allocate(iwork(liwork), STAT=stat)
    if (stat /= 0) call cq_abort ('initDiag: failed to allocate iwork', stat)
    call reg_alloc_mem(area_DM, liwork, type_int)

    return

  end subroutine initDiag
  !!***


  !!****f* DiagModule/endDiag
  !! PURPOSE
  !!   Dealloate arrays used for direct diagonalisation calculations
  !! USAGE
  !!   call endDiag
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/08
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine endDiag

    use memory_module,   only: type_dbl, type_int, type_cplx, &
         reg_dealloc_mem
    use global_module,   only: numprocs, nspin
    use ScalapackFormat, only: deallocate_arrays, proc_rows,  &
         proc_cols, matrix_size

    implicit none

    integer :: i, j, k, stat

    ! Deallocate memory

    deallocate(SCHmat, SCSmat, z, STAT=stat)
    if (stat /= 0) &
         call cq_abort("endDiag: failed to deallocate SCHmat, SCSmat and z", stat)
    call reg_dealloc_mem(area_DM, 3 * row_size * col_size * nspin, type_cplx)

    deallocate(w, occ, STAT=stat)
    if (stat /= 0) call cq_abort('endDiag: failed to deallocate w and occ', stat)
    call reg_dealloc_mem(area_DM, 2 * matrix_size * nkp * nspin, type_dbl)

    deallocate(local_w, STAT=stat)
    if (stat /= 0) call cq_abort('endDiag: failed to allocate local_w', stat)
    call reg_dealloc_mem(area_DM, matrix_size * nspin, type_dbl)

    ! Shut down BLACS

    !    if(me<proc_rows*proc_cols) then
    !call blacs_gridexit(context) ! this is a BLACS library subroutine

    deallocate(ifail, iclustr, STAT=stat)
    if (stat /= 0) &
         call cq_abort("endDiag: failed to deallocate ifail and iclustr", stat)
    call reg_dealloc_mem(area_DM, matrix_size + 2 * proc_rows * &
         proc_cols, type_int)

    deallocate(gap, STAT=stat)
    if (stat /= 0) call cq_abort("endDiag: failed to deallocate gap", stat)
    call reg_dealloc_mem(area_DM, proc_rows * proc_cols, type_dbl)

    deallocate(work, rwork, STAT=stat)
    if (stat /= 0) &
         call cq_abort ('endDiag: failed to deallocate work and rwork', stat)
    call reg_dealloc_mem(area_DM, lwork + lrwork, type_cplx)

    deallocate(iwork, STAT=stat)
    if (stat /= 0) call cq_abort('endDiag: failed to deallocate iwork', stat)
    call reg_dealloc_mem(area_DM, liwork, type_int)
    ! end if

    ! deallocate ScalapackFormat arrays
    !call deallocate_arrays

    return
  end subroutine endDiag
  !!*****


  ! -----------------------------------------------------------------------------
  ! Subroutine PrepareRecv
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/PrepareRecv *
  !!
  !!  NAME
  !!   PrepareRecv
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Prepares to receive data for Scalapack diagonalisation.
  !!
  !!   Scalapack has its own format (discussed in more detail in
  !!   ScalapackFormat module) and we need to distribute the Conquest
  !!   data (stored by row, compressed) to the appropriate processors
  !!   in order to perform Scalapack work.  This routine works out
  !!   how much data we're getting from each processor, and allocates
  !!   memory to store that data.
  !!  INPUTS
  !!   type(DistributeData) :: Distrib - holds arrays created here
  !!  USES
  !!   ScalapackFormat
  !!   group_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   01/03/2002
  !!  MODIFICATION HISTORY
  !!   06/03/2002 dave
  !!    Rewrote in large part so that it makes sense !
  !!   16/04/2002 dave
  !!    Added initial welcome
  !!   23/04/2002 dave
  !!    Changed to use the DistributeData derived type which is passed
  !!    down. Also removed max_rows.
  !!   2004/11/10 drb
  !!    Changed nsf to come from maxima, not common
  !!   2006/08/30 16:49 dave
  !!    Added allocate for arrays in Distrib
  !!  SOURCE
  !!
  subroutine PrepareRecv(Distrib)

    use ScalapackFormat, only: proc_start, SC_row_block_atom,        &
         block_size_r, block_size_c, blocks_r, &
         my_row, matrix_size
    use global_module,   only: iprint_DM, ni_in_cell, numprocs
    use group_module,    only: parts
    use GenComms,        only: myid, cq_abort

    implicit none

    ! Passed variables
    type(DistributeData), intent(out) :: Distrib

    ! Local variables
    integer :: i,j,stat
    integer :: count
    integer :: row, rowblock,proc

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,"Entering PrepareRecv")')
    allocate(Distrib%num_rows(numprocs),Distrib%start_row(numprocs),Distrib%send_rows(numprocs),&
         Distrib%firstrow(numprocs),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating Distrib arrays: ",numprocs,stat)
    ! Prepare counters
    Distrib%num_rows = 0
    count = 1
    i = 1
    ! looping over the local SC rows responsible by local node
    do rowblock=1,blocks_r
       if(my_row(rowblock)>0) then ! If this row block is part of my chunk
          do row = 1,block_size_r
             if(iprint_DM>=5.AND.myid==0) write(io_lun,4) myid,i,rowblock,row
             ! Find processor and increment processors and rows from proc
             proc = parts%i_cc2node(SC_row_block_atom(row,rowblock)%part) ! find proc on which the partition containing the row is stored
             ! remember at the moment data is still stored in CQ format
             if(Distrib%num_rows(proc)==0) Distrib%start_row(proc)=count ! Where data from proc goes
             Distrib%num_rows(proc) = Distrib%num_rows(proc)+1
             if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,i,proc,Distrib%start_row(proc),Distrib%num_rows(proc)
             count = count + 1  ! Position within my chunk of matrix
             i = i+1
             if(i>matrix_size+1) call cq_abort('Matrix too large !')
          end do
       end if
    end do
    return
2   format(10x,'PR Proc: ',i5,' start atom, block: ',2i5)
3   format(10x,'PR Proc: ',i5,' row: ',i5,' proc, start, rows: ',3i5)
4   format(10x,'PR Proc: ',i5,' Row, block, blockrow: ',3i5)
  end subroutine PrepareRecv
  !!***

  ! -----------------------------------------------------------------------------
  ! Subroutine PrepareSend
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/PrepareSend *
  !!
  !!  NAME
  !!   PrepareSend
  !!  USAGE
  !!   PrepareSend(matrix structure of matrix to send)
  !!  PURPOSE
  !!   Sets up arrays which are vital to the sending of data when
  !!   going from CQ format (compressed arrays) to SC format.  Works
  !!   in three stages:
  !!
  !!   i)   Loops over SC matrix rows by going over atoms and support
  !!        functions.  Notes which processors each row is sent to, and
  !!        where the data for a processor starts.
  !!   ii)  Works out how many rows and columns are sent, and allocates
  !!        data storage
  !!   iii) For each non-zero matrix element, works out where in SC
  !!        format matrix it goes
  !!  INPUTS
  !!   type(matrix) :: mat
  !!   type(DistributeData) :: Distrib
  !!  USES
  !!   global_module, common, maxima_module, matrix_module,
  !!   group_module, primary_module, cover_module, ScalapackFormat
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   04/03/2002
  !!  MODIFICATION HISTORY
  !!   08/03/2002 dave
  !!    Changed and reworked so that the number of rows to be sent to
  !!    a processor is found, and loc is stored offset from the start
  !!    of the send buffer (where it should be !) rather than the start
  !!    of the SC block, or even the matrix, as it was.
  !!   15/03/2002 dave
  !!    Tidied a little - mainly improved the allocation of loc and its
  !!    assignment
  !!   05/04/2002 dave
  !!    Many changes, mainly aimed at making it work properly (!).  Now works
  !!    row-by-row (rather than atom-by-atom) so that blocks not equal to nsf
  !!    and variable nsf can be used.  Not convinced that the first loop is
  !!    quite necessary in its current form.
  !!   15/04/2002 dave
  !!    Updated the early stages where we're trying to work out which
  !!    processors we send each row to, and where the data for each
  !!    processor starts.  Simplified the algorithm
  !!    significantly. Also added further comments and changed matrix
  !!    from Hmat to mat.
  !!   16/04/2002 dave
  !!    Added distance vector between atoms to loc so that the phases
  !!    to be applied to the Hamiltonian can be calculated later
  !!   16/04/2002 dave
  !!    Added welcome statement
  !!   23/04/2002 dave
  !!    Added intent to mat
  !!   23/04/2002 dave
  !!    Added a derived datatype (Distrib) which is passed down to
  !!    store all the arrays related to communications - it will allow
  !!    easier generalisation to Conquest
  !!   01/05/2002 dave
  !!    Added level 4 iprint_DM to a write statement
  !!   08:32, 2003/07/28 dave
  !!    Major reworking to allow multiple CQ elements to contribute to
  !!    a single SC elements (should have been there from the start
  !!    but overlooked !)
  !!   2004/10/29 drb
  !!    Added check on size of n_elements before allocating Distrib
  !!   2004/11/10 drb
  !!    Changed nsf to come from maxima, not common
  !!   09:09, 11/05/2005 dave
  !!    Added check on j,k dimensions
  !!   2011/02/13 L.Tong
  !!    Added k-point parallelisation modifications
  !!   2011/12/01 L.Tong
  !!    Removed the redundant parameter matA, it is not used in the subroutine
  !!   2012/06/17 L.Tong
  !!   - Since Distrib contains members that has been allocated
  !!     outside this subrouint, it cannot be set to intent(out) in
  !!     this subroutine, because otherwise Distrib becomes undefined
  !!     upon entering this subroutine, and the array members all
  !!     becomes undefined and unallocated. So set this to intent(inout)
  !!  SOURCE
  !!
  subroutine PrepareSend(range,Distrib)

    use global_module,   only: numprocs, iprint_DM, x_atom_cell,     &
         y_atom_cell, z_atom_cell
    use GenComms,        only: myid, my_barrier
    use maxima_module,   only: maxnsf
    use matrix_module,   only: matrix, matrix_halo
    use group_module,    only: parts
    use primary_module,  only: bundle
    use cover_module,    only: BCS_parts
    use ScalapackFormat, only: CC_to_SC, maxrow, maxcol, proc_block, &
         SC_to_refx, SC_to_refy, block_size_r, &
         block_size_c, blocks_c, proc_start,   &
         proc_groups
    use matrix_data,     only: mat, halo
    use species_module,  only: nsf_species

    implicit none

    ! Passed variables
    integer :: range
    type(DistributeData), intent(inout) :: Distrib

    ! Local variables
    integer, allocatable, dimension(:,:,:) :: sendlist
    integer, allocatable, dimension(:,:,:) :: ele_list
    integer, allocatable, dimension(:)     :: firstcol, sendto
    integer :: part, memb, neigh, ist, atom_num
    integer :: Row_FSC_part, Row_FSC_seq, Col_FSC_part, Col_FSC_seq
    integer :: SCblockr, SCblockc, SCrowc, i, j, k, l, proc, stat, supfn_r, supfn_c
    integer :: maxr, maxc, CC, row, currow, start, gcspart
    integer :: FSCpart, i_acc_prim_nsf, prim_nsf_so_far
    integer :: ng

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,"Entering PrepareSend")')
    ! Initialise and allocate memory
    allocate(sendlist(maxnsf,bundle%n_prim,numprocs),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Failed to alloc sendlist",stat)
    allocate(firstcol(numprocs),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Failed to alloc firstcol",stat)
    allocate(sendto(numprocs),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Failed to alloc sendto",stat)
    ! Zero
    sendto = 0
    Distrib%firstrow = 0
    firstcol = 0
    sendlist = 0
    ! Part (i)
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'Part i ',myid
    call my_barrier()
    ! Work out which processors I send to
    i_acc_prim_nsf = 0   ! CQ format row index accumulator
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions (we loop over the CQ matrix stored on local node)
       if(bundle%nm_nodgroup(part)>0) then
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms in each partition
             atom_num = bundle%nm_nodbeg(part)+memb-1
             Row_FSC_part = CC
             Row_FSC_seq  = memb
             do supfn_r = 1,nsf_species(bundle%species(atom_num)) ! loop over primary set atom support functions
                i_acc_prim_nsf = i_acc_prim_nsf + 1
                ! Here we work out which row of the unpacked SC matrix we're in
                SCblockr = CC_to_SC(Row_FSC_part,Row_FSC_seq,supfn_r)%block_r
                ! Loop over columns of blocks in SC format
                do i=1,blocks_c
                   do ng = 1, proc_groups ! loop over the process groups
                      ! get proc id of the proc responsible for the SC element
                      proc = proc_block(ng, SC_to_refx(SCblockr,i), SC_to_refy(SCblockr,i))
                      sendlist(supfn_r,atom_num,proc) = 1
                      ! If this is the first row and column sent to this processor, then record where data starts
                      if(sendto(proc)==0) then
                         sendto(proc) = 1
                         Distrib%firstrow(proc) = i_acc_prim_nsf ! the first row of the CQ format matrix on local node to send to remote proc
                         firstcol(proc)=(i-1)*block_size_c+1 ! the global SC format col index of the first element to be from local to remote
                      end if ! End if sendto == 0
                   end do ! End do ng = 1, proc_groups
                end do ! End do i=1,blocks_c
             end do ! End do supfn_r
          end do ! End do memb = 1,bundle%nm_nodgroup
       end if ! End if nm_nodgroup>0
    end do ! End do part = 1,groups_on_node
    ! Part (ii)
    ! Work out how many rows to send to each remote processor
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'Part ii ',myid
    call my_barrier()
    Distrib%send_rows = 0
    do proc=1,numprocs ! Loop over processors
       start = 0
       if(sendto(proc)/=0) then ! If we send to this processor
          do i=1,bundle%n_prim  ! loop over primary atoms
             do supfn_r = 1,nsf_species(bundle%species(i))
                if(sendlist(supfn_r,i,proc)==1) then  ! if this row has to be sent to remote proc
                   Distrib%send_rows(proc)=Distrib%send_rows(proc)+1 ! increment the send_rows counter for remote proc
                   if(iprint_DM>=5.AND.myid==0) write(io_lun,11) myid,proc,i,supfn_r
                end if
             end do ! End do supfn_r
          end do ! End do i=bundle%n_prim
       end if ! End if sendto/=0
    end do ! End do proc
    ! Now find out maximum rows and columns sent (for allocating indexing space)
    maxr = 1
    maxc = 1
    do proc=1,numprocs
       if(Distrib%send_rows(proc)>maxr) maxr = Distrib%send_rows(proc)
       if(proc_start(proc)%cols*block_size_c>maxc) maxc = proc_start(proc)%cols*block_size_c
    end do
    if(iprint_DM>=5) then
       if(myid==0) write(io_lun,*) myid,' Allocating ele_list: ',numprocs,maxr,maxc
       call my_barrier()
    end if
    stat = 0
    allocate(ele_list(numprocs,maxr,maxc),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Failed to alloc ele_list ",stat)
    ele_list = 0
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'Part iia ',myid
    call my_barrier()
    ! Added to count elements: how many CQ elements fall onto this one SC one ?
    i_acc_prim_nsf = 0
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(bundle%nm_nodgroup(part)>0) then
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms in each partition
             atom_num = bundle%nm_nodbeg(part)+memb-1
             Row_FSC_part = CC
             Row_FSC_seq  = memb
             prim_nsf_so_far = i_acc_prim_nsf
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom (cols in CQ format on local)
                if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1  ! accumulative atomic index in CQ format (viewed as 1D mem array) on local
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist))
                Col_FSC_seq  = mat(part,range)%i_seq(ist)
                i_acc_prim_nsf = prim_nsf_so_far
                do supfn_r = 1,nsf_species(bundle%species(atom_num))
                   i_acc_prim_nsf = i_acc_prim_nsf + 1
                   row = i_acc_prim_nsf
                   ! Here we work out which row of the unpacked SC matrix we're in
                   SCblockr = CC_to_SC(Row_FSC_part,Row_FSC_seq,supfn_r)%block_r
                   do supfn_c=1,mat(part,range)%ndimj(ist)
                      ! Use CC_to_SC to get SC column and row
                      SCblockc = CC_to_SC(Col_FSC_part,Col_FSC_seq,supfn_c)%block_c
                      SCrowc   = CC_to_SC(Col_FSC_part,Col_FSC_seq,supfn_c)%row_c
                      do ng = 1, proc_groups
                         ! Find processor that this block belongs to
                         proc = proc_block(ng, SC_to_refx(SCblockr,SCblockc), &
                              SC_to_refy(SCblockr,SCblockc))
                         ! Create atom numbers and store posn in array at point
                         j = row-Distrib%firstrow(proc)+1
                         k = (SCblockc-1)*block_size_c + SCrowc-firstcol(proc)+1
                         if(j>maxr.OR.k>maxc.OR.proc>numprocs) then
                            write(io_lun,*) 'Error ! Maxr/c exceeded: ',j,maxr,k,maxc,proc,numprocs
                            call cq_abort('Error in counting elements !')
                         end if
                         if(iprint_DM>=5.AND.myid==0) then
                            write(io_lun,*) myid,' Loop: ',part,memb,supfn_r,supfn_c
                            write(io_lun,*) myid,' j,k,proc: ',j,k,proc
                         end if
                         ele_list(proc,j,k) = ele_list(proc,j,k) + 1
                      end do ! End do ng = 1, proc_groups
                   end do ! End do supfn_c = 1, mat(part,range)%ndimj(ist)
                end do ! End do supfn_r
             end do ! End do neigh
          end do ! End do memb = 1,bundle%nm_nodgroup
       end if ! End if nm_nodgroup>0
    end do ! End do part = 1,groups_on_node
    call my_barrier()
    if(iprint_DM>=5) then
       if(myid==0) write(io_lun,10) myid,maxr,maxc,maxrow,maxcol
       call my_barrier()
    end if
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'Part iib ',myid
    call my_barrier()
    allocate(Distrib%images(numprocs,maxr,maxc),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Failed to alloc images",stat,numprocs*maxr*maxc)
    do i=1,maxc
       do j=1,maxr
          do k = 1,numprocs
             ! I think that I need this as a counter for where to put the data in part (iii)
             Distrib%images(k,j,i)%n_elements = 0
             stat = 0
             if(ele_list(k,j,i)>0) then
                allocate(Distrib%images(k,j,i)%where(ele_list(k,j,i)),STAT=stat)
                if(stat/=0) call cq_abort("DiagModule: Failed to alloc images%where",ele_list(k,j,i),stat)
                do l=1,ele_list(k,j,i)
                   Distrib%images(k,j,i)%where(l)%loci = 0
                   Distrib%images(k,j,i)%where(l)%locj = 0
                   Distrib%images(k,j,i)%where(l)%supfn_row = 0
                   Distrib%images(k,j,i)%where(l)%supfn_col = 0
                   Distrib%images(k,j,i)%where(l)%dx = 0.0_double
                   Distrib%images(k,j,i)%where(l)%dy = 0.0_double
                   Distrib%images(k,j,i)%where(l)%dz = 0.0_double
                end do
             end if ! (ele_list(k,j,i)>0)
          end do
       end do
    end do
    if(iprint_DM>=5) then
       if(myid==0) write(io_lun,*) myid," calling part iii"
       call my_barrier()
    end if
    ! Part (iii)
    ! Build the loc index - for each non-zero element in data_H, record where it goes in SC matrix
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'Part iii ',myid
    call my_barrier()
    i_acc_prim_nsf = 0
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(iprint_DM>=5.AND.myid==0) write(io_lun,1) myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             if(iprint_DM>=5.AND.myid==0) write(io_lun,2) myid,memb
             atom_num = bundle%nm_nodbeg(part)+memb-1
             Row_FSC_part = CC
             Row_FSC_seq  = memb
             prim_nsf_so_far = i_acc_prim_nsf
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist))
                Col_FSC_seq  = mat(part,range)%i_seq(ist)
                gcspart = BCS_parts%icover_ibeg(mat(part,range)%i_part(ist))+Col_FSC_seq-1
                ! Debugging information
                if(iprint_DM>=5.AND.myid==0) write(io_lun,5) myid,Col_FSC_part,Col_FSC_seq
                if(iprint_DM>=5.AND.myid==0) write(io_lun,6) myid,BCS_parts%xcover(gcspart),&
                     BCS_parts%ycover(gcspart), BCS_parts%zcover(gcspart)
                i_acc_prim_nsf = prim_nsf_so_far
                do supfn_r=1,nsf_species(bundle%species(atom_num))
                   i_acc_prim_nsf = i_acc_prim_nsf + 1
                   row = i_acc_prim_nsf
                   do supfn_c=1,mat(part,range)%ndimj(ist)
                      ! Use CC_to_SC to get SC column and row
                      SCblockc = CC_to_SC(Col_FSC_part,Col_FSC_seq,supfn_c)%block_c
                      SCrowc   = CC_to_SC(Col_FSC_part,Col_FSC_seq,supfn_c)%row_c
                      SCblockr = CC_to_SC(Row_FSC_part,Row_FSC_seq,supfn_r)%block_r
                      do ng = 1, proc_groups
                         ! Find processor that this block belongs to
                         proc = proc_block(ng, SC_to_refx(SCblockr,SCblockc), SC_to_refy(SCblockr,SCblockc))
                         if(iprint_DM>=5.AND.myid==0) write(io_lun,4) myid,proc,SCblockr,SCblockc
                         ! Create atom numbers and store posn in array at point
                         j = row-Distrib%firstrow(proc)+1
                         k = (SCblockc-1)*block_size_c + SCrowc-firstcol(proc)+1
                         if(iprint_DM>=5.AND.myid==0) write(io_lun,9) myid,proc,SCblockc,SCrowc,block_size_c,firstcol(proc)
                         ! Create location
                         ! First increment (and check) number of elements
                         Distrib%images(proc,j,k)%n_elements = Distrib%images(proc,j,k)%n_elements + 1
                         if(Distrib%images(proc,j,k)%n_elements > ele_list(proc,j,k)) &
                              call cq_abort('Overrun in Distrib%images !',ele_list(proc,j,k))
                         ! Store location
                         Distrib%images(proc,j,k)%where(Distrib%images(proc,j,k)%n_elements)%loci = atom_num
                         Distrib%images(proc,j,k)%where(Distrib%images(proc,j,k)%n_elements)%locj = halo(range)%i_halo(gcspart)
                         Distrib%images(proc,j,k)%where(Distrib%images(proc,j,k)%n_elements)%supfn_row = supfn_r
                         Distrib%images(proc,j,k)%where(Distrib%images(proc,j,k)%n_elements)%supfn_col = supfn_c
                         ! Build the distances between atoms - needed for phases
                         Distrib%images(proc,j,k)%where(Distrib%images(proc,j,k)%n_elements)%dx = &
                              BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
                         Distrib%images(proc,j,k)%where(Distrib%images(proc,j,k)%n_elements)%dy = &
                              BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
                         Distrib%images(proc,j,k)%where(Distrib%images(proc,j,k)%n_elements)%dz = &
                              BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
                         ! Change distances here: we now use displacement between supercells
                         !                      FSCpart = BCS_parts%lab_cell(mat(part,range)%i_part(ist))!gcspart)
                         ! Here we assume that j_0 is in the FSC, as is i_0
                         !                      write(io_lun,*) myid,' FSCpart, atom and xyz: ', &
                         !                           FSCpart,parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1,&
                         !                           x_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1),&
                         !                           y_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1),&
                         !                           z_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
                         !                      Distrib%loc(proc,j,k)%dx = BCS_parts%xcover(gcspart)- &
                         !                           x_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
                         !                      Distrib%loc(proc,j,k)%dy = BCS_parts%ycover(gcspart)- &
                         !                           y_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
                         !                      Distrib%loc(proc,j,k)%dz = BCS_parts%zcover(gcspart)- &
                         !                           z_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
                         !                      if(iprint_DM>=5) write(io_lun,7) myid,proc,j,k,Distrib%loc(proc,j,k)
                         if((j>maxr.OR.k>maxc).AND.myid==0) write(io_lun,*) 'Problem! j,k,max row,col: ',j,k,maxr,maxc
                      end do ! End do ng = 1, proc_groups
                   end do ! End do supfn_c = 1,nsf
                end do ! End do supfn_r = 1,nsf
             end do ! End do neigh=1,mat%n_nab
          end do ! End do memb =1,nm_nodgroup
       end if ! End if nm_nodgroup > 0
    end do ! End do part=1,groups_on_node
    call my_barrier()
    deallocate(ele_list,STAT=stat)
    deallocate(sendto,STAT=stat)
    deallocate(firstcol,STAT=stat)
    deallocate(sendlist,STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Failed to alloc sendto",stat)
    return
1   format(10x,'Processor: ',i7,' Partition: ',i5)
2   format(10x,'Processor: ',i7,' Atom: ',i9)
3   format(10x,'Processor: ',i7,' Neighbour: ',i5)
4   format(10x,'Processor: ',i7,' Recv: ',i5,' SC blocks: ',2i5)
5   format(10x,'Processor: ',i7,' FSC part, seq: ',2i5)
6   format(10x,'Processor: ',i7,' Neigh xyz: ',3f15.10)
7   format(10x,'Processor: ',i7,' Recv Proc: ',i7,' PS Atom j,k,loc: ',5i5,' dr: ',3f10.5)
8   format(10x,'Processor: ',i7,' Recv Proc: ',i7,' First row, col: ',2i5)
9   format(10x,'Processor: ',i7,' Recv Proc: ',i7,' SC block, atom, size: ',3i5,' firstcol: ',i5)
10  format(10x,'Processor: ',i7,' Maxr,c: ',2i5,' Maxrow,col: ',2i5)
11  format(10x,'Processor: ',i7,' PrepSend to : ',i7,' Atom, Supfn: ',2i9)
  end subroutine PrepareSend
  !!***

  ! -----------------------------------------------------------------------------
  ! Subroutine DistributeCQ_to_SC
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/DistributeCQ_to_SC *
  !!
  !!  NAME
  !!   DistributeCQ_to_SC
  !!  USAGE
  !!   DistributeCQ_to_SC(matrix maximum, matrix)
  !!  PURPOSE
  !!   Distributes data stored in Conquest compressed row format to
  !!   processors ready for a Scalapack operation (e.g. diagonalisation)
  !!
  !!   Operates cyclically - each processor starts by redistributing local
  !!   data, then increments the processor they send to and decrements the
  !!   processor they receive from at each iteration.  This is an N^2
  !!   process, but we're going to diagonalise, which is N^3, so it's not
  !!   really very important.
  !!  INPUTS
  !!   integer :: matrix maximum
  !!   integer :: matrix - matrix to be sent
  !!  USES
  !!   datatypes, common, global_module, mpi, numbers,
  !!   ScalapackFormat, maxima_module, GenComms
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   01/03/2002
  !!  MODIFICATION HISTORY
  !!   08/03/2002 dave
  !!    Changed for debugging - rewrote where the data comes from
  !!    and where it goes to !
  !!   11/03/2002 dave
  !!    Added many further write statements and changed how start_row
  !!    was used - it refers to an atom row, while ParaDens works in
  !!    terms of orbitals (!) - so I had to scale by nsf.  This needs
  !!    generalisation
  !!   15/03/2002 dave
  !!    Tidied output and comments
  !!   05/04/2002 dave
  !!    Major rewriting to ensure that it works - proceeds row-by-row
  !!    now and sends chunks which are the correct size on both sending
  !!    and receiving processors.
  !!   15/04/2002 dave
  !!    Changed so that the matrix to be distributed is passed as an
  !!    argument (allows more than one type of matrix to be sent) and
  !!    added more comments
  !!   16/04/2002 dave
  !!    Added k-point to arguments passed
  !!   16/04/2002 dave
  !!    Added phase calculation to construction of Hamiltonian and
  !!    made all appropriate variables complex
  !!   18/04/2002 dave
  !!    Added matrix we distribute into to list of arguments
  !!    (i.e. allow distribution of more than one matrix) Also removed
  !!    TODO entry, as I'd done it
  !!   23/04/2002 dave
  !!    Added intent to passed variables
  !!   23/04/2002 dave
  !!    Added DistributeData variable to carry details of
  !!    matrix-specific arrays for distribution
  !!   24/04/2002 dave
  !!    Fixed a problem with accumulation of H for on-processor transfers
  !!   08:18, 2003/07/28 dave
  !!    Completely rewrote the accumulation loop to account for
  !!    multiple CQ elements adding onto a single SC element.  Changed
  !!    Distrib structure, introduced loop over CQ elements, reworked
  !!    various checks
  !!   2004/11/10 drb
  !!    Changed nsf to come from maxima, not common
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2017/10/18 14:46 dave
  !!    Added optional argument kpassed to allow for calling of diagonalisation for
  !!    arbitrary k-points (NB at present ONLY works for ONE process group).  Changed
  !!    to define k-point location at start of routine
  !!   2018/10/08 16:38 dave
  !!    Changing MPI tags to conform to MPI standard
  !!  SOURCE
  !!
  subroutine DistributeCQ_to_SC(Distrib,matA,pgk,SCmat,kpassed)

    use datatypes
    use global_module,   only: numprocs, iprint_DM
    use mpi
    use numbers,         only: zero, minus_i
    use ScalapackFormat, only: proc_start, block_size_r, block_size_c,&
         pgid, proc_groups, pg_kpoints, pgroup, &
         N_kpoints_in_pg
    use GenComms,        only: my_barrier, myid
    use GenBlas,         only: copy
    use mult_module,     only: return_matrix_value_pos, matrix_pos

    implicit none

    ! Passed variables
    type(DistributeData), intent(in)  :: Distrib
    integer,              intent(in)  :: matA  ! the CQ format matrix on local node
    integer,              intent(in)  :: pgk  ! k-point index within process group
    complex(double_cplx), intent(out) :: SCmat(:,:)  ! the SC format sub matrix on local node
    real(double), dimension(3), OPTIONAL :: kpassed

    ! Local variables
    integer :: send_proc,recv_proc,send_size,recv_size,send_pgid
    integer :: sendtag, recvtag, stat
    integer :: srow_size,scol_size,rrow_size,rcol_size
    integer :: req, i, j, k, l, ierr, wheremat
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    real(double) :: phase,rfac,ifac,hreal,himag, kx, ky, kz, send_kx, send_ky, send_kz

    ! Work out k-point - has it been passed in explicitly ?
    if(PRESENT(kpassed)) then
       kx = kpassed(1)
       ky = kpassed(2)
       kz = kpassed(3)
    else if (pgk <= N_kpoints_in_pg(pgid)) then
       kx = kk(1,pg_kpoints(pgid,pgk))
       ky = kk(2,pg_kpoints(pgid,pgk))
       kz = kk(3,pg_kpoints(pgid,pgk))
    else
       kx = zero
       ky = zero
       kz = zero
    end if
    send_proc = myid
    recv_proc = myid
    ! Distribute data: loop over processors and issue sends and
    ! receives as appropriate
    SCmat = zero
    send_size = 0
    recv_size = 0
    ! We only send one message to each process
    sendtag = 1
    recvtag = 1
    call start_timer(tmr_std_matrices)
    do i=1,numprocs
       if(iprint_DM>=4.AND.myid==0) write(io_lun,1) myid,i,send_proc,recv_proc
       ! Sizes
       ! work out send_proc group id
       send_pgid = pgroup(send_proc+1)
       if (pgk <= N_kpoints_in_pg(send_pgid)) then
          srow_size = Distrib%send_rows(send_proc+1)            ! number of CQ rows to send to remote
          scol_size = proc_start(send_proc+1)%cols*block_size_c ! number of SC cols responsible by remote, and hence the no. of cols to send
          send_size = srow_size*scol_size     ! size of send buffer
          allocate(SendBuffer(srow_size,scol_size),STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Can't alloc SendBuffer",stat)
          ! Zero SendBuffer
          SendBuffer = cmplx(zero,zero,double_cplx)
          send_kx = kk(1,pg_kpoints(send_pgid,pgk))
          send_ky = kk(2,pg_kpoints(send_pgid,pgk))
          send_kz = kk(3,pg_kpoints(send_pgid,pgk))
       end if
       if (pgk <= N_kpoints_in_pg(pgid)) then
          rrow_size = Distrib%num_rows(recv_proc+1)             ! numbef of rows to be received by local from remote
          rcol_size = proc_start(myid+1)%cols*block_size_c      ! number of SC cols responsible by local and hence the no. of cols to receive
          recv_size = rrow_size*rcol_size     ! size of remote buffer
          allocate(RecvBuffer(rrow_size,rcol_size),STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Can't alloc RecvBuffer",stat)
          ! Zero RecvBuffer
          RecvBuffer = cmplx(zero,zero,double_cplx)
       end if

       ! On-site
       if(send_proc==myid.AND.recv_proc==myid) then
          ! no need to send and receive data to and from onsite processor
          if (pgk <= N_kpoints_in_pg(pgid)) then
             if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'num_rows, send_rows: ',Distrib%num_rows(myid+1),Distrib%send_rows(myid+1)
             if(iprint_DM>=5.AND.myid==0) write(io_lun,11) myid,srow_size,scol_size,send_size,rrow_size,rcol_size,recv_size
             ! Fill local copy of SCmat
             do j=1,srow_size     ! loop over the send buffer rows and cols (i.e. elements),
                do k=1,scol_size  ! the content destined for send buffer goes directly to SCmat
                   if(Distrib%images(send_proc+1,j,k)%n_elements>0) then  ! if the element is a part of atom that is on neighbourlist
                      do l=1,Distrib%images(send_proc+1,j,k)%n_elements   ! loop over the periodic images of the element
                         wheremat = matrix_pos(matA,Distrib%images(send_proc+1,j,k)%where(l)%loci, &
                              Distrib%images(send_proc+1,j,k)%where(l)%locj,&
                              Distrib%images(send_proc+1,j,k)%where(l)%supfn_row,&
                              Distrib%images(send_proc+1,j,k)%where(l)%supfn_col) ! work out where is the (image) element located on local node
                         if(iprint_DM>=5.AND.myid==0) write(io_lun,7) myid,send_proc,j,k,wheremat
                         ! for onsite terms, we need to work with k-points in the proc_group pgid
                         ! more precisely we work with pgk-th kpoint responsible by proc_group pgid
                         phase = kx*Distrib%images(send_proc+1,j,k)%where(l)%dx + &
                              ky*Distrib%images(send_proc+1,j,k)%where(l)%dy + &
                              kz*Distrib%images(send_proc+1,j,k)%where(l)%dz
                         rfac = cos(phase)* return_matrix_value_pos(matA,wheremat)
                         ifac = sin(phase)* return_matrix_value_pos(matA,wheremat)
                         ! Care here - we need to accumulate
                         SCmat(Distrib%start_row(myid+1)+j-1,k) = SCmat(Distrib%start_row(myid+1)+j-1,k) + &
                              cmplx(rfac,ifac,double_cplx)
                      end do ! Distrib%images%n_elements
                   end if ! n_elements>0
                end do ! k=1,scol_size
             end do ! j=1,srow_size
             if(iprint_DM>=5.AND.myid==0) write(io_lun,8) myid,Distrib%start_row(recv_proc+1)
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) '  Done on-proc'
             deallocate(RecvBuffer,STAT=stat)  ! the send and recv buffer sizes are different for different remote processors
             deallocate(SendBuffer,STAT=stat)  ! so we need to allocate and deallocate for every proc in the i loop.
             ! Send and receive data to/from remote processors
          end if ! End if (pgk <= N_kpoints_in_pg(pgid)) then
       else ! if(send_proc==myid.AND.recv_proc==myid)
          ! ---------------
          ! Fill SendBuffer
          ! ---------------
          if (pgk <= N_kpoints_in_pg(send_pgid)) then
             if(send_size>0) then
                do j=1,srow_size
                   do k=1,scol_size
                      if(Distrib%images(send_proc+1,j,k)%n_elements>0) then
                         do l=1,Distrib%images(send_proc+1,j,k)%n_elements
                            wheremat = matrix_pos(matA,Distrib%images(send_proc+1,j,k)%where(l)%loci, &
                                 Distrib%images(send_proc+1,j,k)%where(l)%locj,&
                                 Distrib%images(send_proc+1,j,k)%where(l)%supfn_row,&
                                 Distrib%images(send_proc+1,j,k)%where(l)%supfn_col)
                            if(iprint_DM>=5.AND.myid==0) write(io_lun,7) myid,send_proc,j,k,wheremat
                            ! the kpoints used should be that responsible by the remote node, i.e.
                            ! corresponding to the pgk-th kpoint in proc_group that contains send_proc
                            phase = send_kx*Distrib%images(send_proc+1,j,k)%where(l)%dx + &
                                    send_ky*Distrib%images(send_proc+1,j,k)%where(l)%dy + &
                                    send_kz*Distrib%images(send_proc+1,j,k)%where(l)%dz
                            rfac = cos(phase)*return_matrix_value_pos(matA,wheremat)
                            ifac = sin(phase)*return_matrix_value_pos(matA,wheremat)
                            ! Accumulate the data
                            SendBuffer(j,k)= SendBuffer(j,k)+cmplx(rfac,ifac,double_cplx)
                         end do ! l=Distrib%images%n_elements
                      end if ! n_elements>0
                   end do ! j=srow_size
                end do ! k=scol_size
             end if ! if send_size>0
          end if ! End if (pgk <= N_kpoints_in_pg(send_pgid)) then
          ! ---------------
          ! Do the transfer
          ! ---------------
          if (pgk <= N_kpoints_in_pg(send_pgid)) then
             ! Debugging output
             if(iprint_DM>=5.AND.myid==0) then
                write(io_lun,*) 'Proc: ',myid,' Sizes: ',send_size,recv_size
                write(io_lun,10) i,myid,send_proc,recv_proc,srow_size, scol_size, rrow_size, rcol_size
                write(io_lun,*) 'About to send...'
             end if
             ! Issue non-blocking send and then receive
             if(send_size>0) then
                call MPI_isend(SendBuffer,send_size,MPI_DOUBLE_COMPLEX,&
                     send_proc,sendtag,MPI_COMM_WORLD,req,ierr)
             end if
          end if
          if (pgk <= N_kpoints_in_pg(pgid)) then
             if(recv_size>0) then
                call MPI_recv(RecvBuffer,recv_size,MPI_DOUBLE_COMPLEX,&
                     recv_proc,recvtag,MPI_COMM_WORLD,mpi_stat,ierr)
             end if
             ! More debugging output
             if(iprint_DM>=5.AND.myid==0) then
                write(io_lun,8) myid,Distrib%start_row(recv_proc+1)
             end if
             if(iprint_DM>=4.AND.myid==0) write(io_lun,2) myid,i
             ! Put the data in place from the receive buffer
             if(iprint_DM>=4.AND.myid==0) &
                  write(io_lun,*) myid,'rowsize, start, col: ',rrow_size, Distrib%start_row(recv_proc+1), rcol_size
             if(rrow_size > 0) then
                do j=1,rrow_size
                   SCmat(Distrib%start_row(recv_proc+1)+j-1,1:rcol_size) = RecvBuffer(j,1:rcol_size)
                end do
             end if ! (rrow_size > 0)
          end if
          ! Now wait for the non-blocking send to finish before deallocating !
          if(iprint_DM>=4.AND.myid==0) write(io_lun,13) myid,i
          ! only call MPI_Wait if the isend is called before
          if (pgk <= N_kpoints_in_pg(send_pgid)) then
             if(send_size>0) call MPI_Wait(req,mpi_stat,ierr)
          end if
          if (pgk <= N_kpoints_in_pg(pgid)) then
             deallocate(RecvBuffer,STAT=stat)
             if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
          end if
          if (pgk <= N_kpoints_in_pg(send_pgid)) then
             deallocate(SendBuffer,STAT=stat)
             if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
          end if
          if(iprint_DM>=4.AND.myid==0) write(io_lun,12) myid,i
       end if ! else part of if(send_proc==myid.AND.recv_proc==myid)
       ! Increment/decrement recv and send, and wrap
       ! Remember that we go from 0->numprocs-1
       send_proc = send_proc +1
       if(send_proc.GT.numprocs-1) send_proc = 0
       recv_proc = recv_proc -1
       if(recv_proc.LT.0) recv_proc = numprocs-1
       call my_barrier()
    end do ! End loop over processors
    call stop_timer(tmr_std_matrices)
    return
1   format(10x,'Proc: ',i5,' Iter: ',i5,' Send/Recv: ',2i5)
2   format(10x,'Proc: ',i5,' done send/recv',i5)
3   format(10x,'Proc: ',i5,' i,j: ',2i5,' Data: ',4f15.10)
4   format(10x,'CQ2SC Proc: ',i5,' i,j: ',2i5,' Hmat: ',4f15.10)
7   format(10x,'Processor: ',i5,' Recv Proc: ',i5,' D Atom j,k,loc: ',3i5)
8   format(10x,'Proc: ',i5,' Starting row for data: ',i5)
9   format(10x,'CQ2SC Proc: ',i5,' i,j: ',2i5,' RecvBuff: ',4f15.10)
10  format(10x,i4,'CQ2SCProc: ',i5,' To/From ',2i5,' Rows, Cols: ',4i5)
11  format(10x,'On-site Proc: ',i5,' Send row,col,size: ',3i8,' Recv row,col,size: ',3i8)
12  format(10x,'Proc: ',i5,' done dealloc ',i5)
13  format(10x,'Proc: ',i5,' calling MPI_Wait ',i5)
  end subroutine DistributeCQ_to_SC
  !!***

  ! -----------------------------------------------------------------------------
  ! Subroutine DistributeSC_to_ref
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/DistributeSC_to_ref *
  !!
  !!  NAME
  !!   DistributeSC_to_ref - send data back to processors for reference format
  !!  USAGE
  !!   DistributeSC_to_ref(Scalapack formatted eigenvector chunk,
  !!                       Conquest formatted eigenvectore chunk)
  !!  PURPOSE
  !!   Once a diagonalisation has been performed, we want eigenvector
  !!   coefficients for every band for a given support function to be
  !!   returned, in order of energy, to the processor responsible for
  !!   the support function.  In other words, we have \psi_n =
  !!   \sum_{i\alpha} c^n_{i\alpha} \phi_{i\alpha} (roughly speaking)
  !!   and we want all n values of c^n_{i\alpha} on the processor
  !!   responsible for i.
  !!
  !!   This is reasonably simple to accomplish - the data transfer is
  !!   done just by reversing the calls of DistributeCQ_to_SC.
  !!   However, the data that arrives back is then ordered correctly
  !!   in terms of rows but has its columns in SC format - we want
  !!   them in reference format ! So we use mapy (I think) to do this.
  !!   For a chunk which has come back, we know that we have all
  !!   columns, so we work out which row block we're in (for a given
  !!   row) and loop over columns, mapping to reference block and
  !!   hence position in chunk.
  !!  INPUTS
  !!   real(double), dimension(:,:) :: SCeig - the piece of eigenvector
  !!                                           matrix created by Scalapack
  !!   real(double), dimension(:,:) :: localEig - local storage for
  !!                                              eigenvector coefficients
  !!                                              for atoms
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   17/04/2002
  !!  MODIFICATION HISTORY
  !!   18/04/2002 dave
  !!    More documentation and finished off decoding of RecvBuffer (I
  !!    think)
  !!   23/04/2002 dave
  !!    Added intents
  !!   23/04/2002 dave
  !!    Added DistributeData derived type to carry information about
  !!    matrix-specific data - though actually here I think that the
  !!    only bits used (firstrow) are matrix-independent.
  !!   2004/11/10 drb
  !!    Changed nsf to come from maxima, not common
  !!   2006/11/02 16:59 dave
  !!    Bug fix: only call MPI_Wait if MPI_isend called
  !!    (i.e. send_size greater than zero)
  !!   2011/02/13 L.Tong
  !!    Added k-point parallelisation
  !!   2018/10/08 16:38 dave
  !!    Changing MPI tags to conform to MPI standard
  !!  SOURCE
  !!
  subroutine DistributeSC_to_ref(Distrib,ng,SCeig,localEig)

    use datatypes
    use global_module,   only: numprocs, iprint_DM
    use mpi
    use numbers,         only: zero, minus_i
    use ScalapackFormat, only: proc_start, block_size_r, block_size_c,&
         mapy, pgroup, pgid
    use GenComms,        only: my_barrier, myid

    implicit none

    ! Passed variables
    type(DistributeData), intent(in) :: Distrib
    integer, intent(in) :: ng  ! process group index
    complex(double_cplx), dimension(:,:), intent(in)  :: SCeig    ! the local submatrix in SC storage format for eigenvector matrix
    complex(double_cplx), dimension(:,:), intent(out) :: localEig ! the local submatrix in CQ storage format for eigenvector matrix

    ! Local variables
    integer :: send_proc, recv_proc, send_size, recv_size
    integer :: sendtag, recvtag, stat, rblock, cblock, refblock, roff,&
         coff, req1, req2, ierr
    integer :: srow_size, scol_size, rrow_size, rcol_size
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    integer :: i, j, k

    send_proc = myid
    recv_proc = myid
    ! Distribute data: loop over processors and issue sends and
    ! receives as appropriate. We are only going to receive data from
    ! process group ng, but will send to all processors
    localEig = zero
    send_size = 0
    recv_size = 0
    ! We send two messages; below, we use tag and tag+1
    sendtag = 1
    recvtag = 1
    do i=1,numprocs
       if(iprint_DM>=4.AND.myid==0) write(io_lun,1) myid,i,send_proc,recv_proc
       ! Sizes
       if (pgroup(recv_proc+1) == ng) then ! we only need to receive from group ng
          ! note that the number of rows to receive from remote is the same as number
          ! of rows sent to remote
          rrow_size = Distrib%send_rows(recv_proc+1)                    ! number of rows to receive from remote
          rcol_size = proc_start(recv_proc+1)%cols*block_size_c         ! number of cols to receive from remote
          recv_size = rrow_size*rcol_size                               ! RecvBuffer size
          allocate(RecvBuffer(rrow_size,rcol_size),STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Can't alloc RecvBuffer",stat)
          ! Zero RecvBuffer
          RecvBuffer = cmplx(zero,zero,double_cplx)
       end if
       if (pgid == ng) then ! only if the local processor is in group ng do we have to send
          srow_size = Distrib%num_rows(send_proc+1)                     ! number of rows to send to remote
          scol_size = proc_start(myid+1)%cols*block_size_c              ! number of cols to send to remote
          send_size = srow_size*scol_size                               ! SendBuffer size
          ! Allocate memory
          allocate(SendBuffer(srow_size,scol_size),STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Can't alloc SendBuffer",stat)
          ! Zero SendBuffer
          SendBuffer = cmplx(zero,zero,double_cplx)
       end if

       if(iprint_DM>=5.AND.myid==0) write(io_lun,8) myid,Distrib%start_row(recv_proc+1)

       ! On-site
       if(send_proc==myid.AND.recv_proc==myid) then
          if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'num_rows, send_rows: ',Distrib%num_rows(myid+1),Distrib%send_rows(myid+1)
          if(iprint_DM>=5.AND.myid==0) write(io_lun,11) myid,srow_size,scol_size,send_size,rrow_size,rcol_size,recv_size
          ! Fill send buffer
          ! Normally we'd then send this, and copy our receive buffer out
          if (pgroup(recv_proc + 1) == ng) then ! we only receive from processes in group ng
             do j=1,srow_size
                ! SendBuffer(j,1:scol_size) = SCeig(start_row(send_proc+1)+j-1,1:scol_size)
                RecvBuffer(j,1:scol_size) = SCeig(Distrib%start_row(send_proc+1)+j-1,1:scol_size)
                ! note that Distrib%start_row is not to be confused with Distribute%first_row
             end do
             roff = Distrib%start_row(send_proc+1) ! note what used to be recv_proc in DistribCQ_to_SC is now send_proc
             ! Decode receive buffer into local eigenvector store
             do j=1,rrow_size
                rblock = floor(real((roff-1+j-1)/block_size_r))+1  ! get block row ind corresponding to reveiving row
                do k=1,rcol_size,block_size_c
                   cblock = floor(real((k-1)/block_size_c))+1
                   refblock = mapy(recv_proc+1,rblock,cblock)
                   coff = (refblock-1)*block_size_c + 1
                   if(iprint_DM>=5.AND.myid==0) &
                        write(io_lun,3) myid,j,k,rblock,cblock,refblock,coff,Distrib%firstrow(recv_proc+1),RecvBuffer(j,k)
                   ! localEig(Distrib%firstrow(recv_proc+1)+j-1,coff:coff+block_size_c-1) = RecvBuffer(j,k:k+block_size_c-1)
                   localEig(coff:coff+block_size_c-1,Distrib%firstrow(recv_proc+1)+j-1) = RecvBuffer(j,k:k+block_size_c-1)
                end do
             end do
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) '  Done on-proc'
             deallocate(RecvBuffer,STAT=stat)
          end if
          if (pgid == ng) deallocate(SendBuffer,STAT=stat)
          ! Send and receive data to/from remote processors
       else ! if(send_proc==myid.AND.recv_proc==myid)
          ! ---------------
          ! Fill SendBuffer
          ! ---------------
          if (pgid == ng)  then ! we only send if local is one of processes in proc_group ng
             if(send_size>0) then
                do j=1,srow_size
                   SendBuffer(j,1:scol_size) = SCeig(Distrib%start_row(send_proc+1)+j-1,1:scol_size)
                end do
             end if ! if send_size>0
             ! ---------------
             ! Do the transfer
             ! ---------------
             ! Debugging output
             if(iprint_DM>=5.AND.myid==0) then
                write(io_lun,*) 'Proc: ',myid,' Sizes: ',send_size,recv_size
                write(io_lun,10) i,myid,send_proc,recv_proc,srow_size, scol_size, rrow_size, rcol_size
                write(io_lun,*) 'About to send...'
             end if
             ! Issue non-blocking send and then receive
             if(send_size>0) then
                call MPI_isend(Distrib%start_row(send_proc+1),1,MPI_INTEGER,send_proc,sendtag,MPI_COMM_WORLD,req1,ierr)
                ! we need to send Distrib%start_row to send_proc because the remote need to know this info to get roff
                call MPI_isend(SendBuffer,send_size,MPI_DOUBLE_COMPLEX,send_proc,sendtag+1,MPI_COMM_WORLD,req2,ierr)
             end if
          end if ! End if (pgid == ng)

          if (pgroup(recv_proc + 1) == ng) then ! we only receive from procs in proc_group ng
             if(recv_size>0) then
                call MPI_recv(roff,1,MPI_INTEGER,recv_proc,recvtag,MPI_COMM_WORLD,mpi_stat,ierr)
                call MPI_recv(RecvBuffer,recv_size,MPI_DOUBLE_COMPLEX,recv_proc,recvtag+1,MPI_COMM_WORLD,mpi_stat,ierr)
             end if
             if(iprint_DM>=4.AND.myid==0) write(io_lun,2) myid,roff
             ! Put the data in place from the receive buffer
             if(rrow_size > 0) then
                do j=1,rrow_size
                   rblock = floor(real((roff-1+j-1)/block_size_r))+1
                   do k=1,rcol_size,block_size_c
                      cblock = floor(real((k-1)/block_size_c))+1
                      refblock = mapy(recv_proc+1,rblock,cblock)
                      coff = (refblock-1)*block_size_c + 1
                      if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,j,k,rblock,cblock,refblock,coff,&
                           Distrib%firstrow(recv_proc+1),RecvBuffer(j,k)
                      !localEig(Distrib%firstrow(recv_proc+1)+j-1,coff:coff+block_size_c-1) = RecvBuffer(j,k:k+block_size_c-1)
                      localEig(coff:coff+block_size_c-1,Distrib%firstrow(recv_proc+1)+j-1) = &
                           RecvBuffer(j,k:k+block_size_c-1)
                   end do
                end do
             end if ! (rrow_size > 0)
          end if ! End if (pgroup(recv_proc + 1) == ng)
          ! Now wait for the non-blocking send to finish before deallocating !
          if (pgid == ng)  then
             if(send_size>0) then
                call MPI_Wait(req1,mpi_stat,ierr)
                call MPI_Wait(req2,mpi_stat,ierr)
             end if
          end if
          if (pgroup(recv_proc+1) == ng) then
             deallocate(RecvBuffer,STAT=stat)
             if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
          end if
          if (pgid == ng) then
             deallocate(SendBuffer,STAT=stat)
             if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
          end if
       end if ! else part of if(send_proc==myid.AND.recv_proc==myid)
       ! Increment/decrement recv and send, and wrap
       ! Remember that we go from 0->numprocs-1
       send_proc = send_proc +1
       if(send_proc.GT.numprocs-1) send_proc = 0
       recv_proc = recv_proc -1
       if(recv_proc.LT.0) recv_proc = numprocs-1
       call my_barrier()
    end do ! End loop over processors
    return
1   format(10x,'SC2Ref Proc: ',i5,' Iter: ',i5,' Send/Recv: ',2i5)
2   format(10x,'SC2Ref Proc: ',i5,' done send/recv Offset: ',i5)
3   format(10x,'SC2Ref Proc: ',i5,' Coords: ',2i5,' Block: ',2i5,' Refbloc,coff,firstro: ',3i15,2f15.8)
8   format(10x,'Proc: ',i5,' Starting row for data: ',i5)
10  format(10x,i4,'SC2Ref Proc: ',i5,' To/From ',2i5,' Rows, Cols: ',4i5)
11  format(10x,'SC2Ref On-site Proc: ',i5,' Send row,col,size: ',3i8,' Recv row,col,size: ',3i8)
  end subroutine DistributeSC_to_ref
  !!***

  ! -----------------------------------------------------------------------------
  ! Subroutine findFermi
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/findFermi
  !! PURPOSE
  !!   Finds the fermi level given a set of eigenvalues at a number of k
  !!   points.
  !!
  !!   This is replaces the original findFermi subroutine by D.Bowler
  !! USAGE
  !!   call findFermi(electrons, eig, nbands, nkp, Ef, occ)
  !! INPUTS
  !!   integer nbands                     : number of bands
  !!   integer nkp                        : number of kpoints
  !!   real(double) eig(nbands,nkp,nspin) : eigenvalues
  !!   real(double) electrons(nspin)      : number of electrons in each spin
  !!                                        channel which the fermi energy is
  !!                                        required to get.
  !! OUTPUT
  !!   real(double) Ef(nspin) : Fermi energy for each spin channel
  !!   real(double) occ(nbands,nkp,nspin) : occupancies
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/06
  !! MODIFICATION HISTORY
  !!  2013/02/08 08:09 dave and UT
  !!  - Adding simple DeltaSCF excitations
  !! SOURCE
  !!
  subroutine findFermi(electrons, eig, nbands, nkp, Ef, occ)

    use datatypes
    use numbers
    use global_module, only: iprint_DM, nspin, spin_factor, &
         flag_fix_spin_population, flag_DeltaSCF, flag_excite, &
         dscf_source_level, dscf_target_level, dscf_source_spin, &
         dscf_target_spin, dscf_source_nfold, dscf_target_nfold
    use GenComms,      only: myid

    implicit none

    ! passed variables
    integer,                        intent(in)  :: nbands, nkp
    real(double), dimension(:,:,:), intent(in)  :: eig
    real(double), dimension(:),     intent(in)  :: electrons
    real(double), dimension(:),     intent(out) :: Ef
    real(double), dimension(:,:,:), intent(out) :: occ

    ! local variables
    real(double)            :: electrons_total
    real(double), parameter :: tolElec = 1.0e-6_double
    real(double) :: locals_occ, localt_occ
    integer :: ikp, i, ispin

    if (nspin == 2) then
       electrons_total = electrons(1) + electrons(2)
    else
       electrons_total = two * electrons(1)
    end if
    if (flag_fix_spin_population .or. nspin == 1) then
       call findFermi_fixspin(electrons, eig, nbands, nkp, Ef, occ)
    else
       call findFermi_varspin(electrons_total, eig, nbands, nkp, Ef, occ)
    end if
    if(flag_DeltaSCF.AND.flag_excite) then
       if(nspin==1) then
          if(dscf_target_spin==2) dscf_target_spin=1
          if(dscf_source_spin==2) dscf_source_spin=1
       end if
       localt_occ = one/real(spin_factor*dscf_target_nfold,double)
       locals_occ = one/real(spin_factor*dscf_source_nfold,double)
       do ikp = 1,nkp
          do i=1,dscf_target_nfold
             occ(dscf_target_level+i-1,ikp,dscf_target_spin) = wtk(ikp)*localt_occ
          end do
          do i=1,dscf_source_nfold
             occ(dscf_source_level+i-1,ikp,dscf_source_spin) = wtk(ikp)*(one - locals_occ)
          end do
       end do
    end if
    return
  end subroutine findFermi
  !!*****


  !!****f* DiagModule/findFermi_fixspin
  !! PURPOSE
  !!   Finds Fermi energy for the case when the populations in the
  !!   different spin channels are fixed.
  !!
  !!   Note that spin non-polarised case is equivalent to the case
  !!   with spin population fixed at 50% each. So it also uses this
  !!   subroutine.
  !! USAGE
  !!   call findFermi_fixspin(electrons, eig, nbands, nkp, Ef, occ)
  !! INPUTS
  !!   integer nbands                     : number of bands
  !!   integer nkp                        : number of kpoints
  !!   real(double) eig(nbands,nkp,nspin) : eigenvalues
  !!   real(double) electrons(nspin)      : number of electrons in each spin
  !!                                        channel which the fermi energy is
  !!                                        required to get.
  !! OUTPUT
  !!   real(double) Ef(nspin)             : Fermi energy for each spin channel
  !!   real(double) occ(nbands,nkp,nspin) : occupancies
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/06
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine findFermi_fixspin(electrons, eig, nbands, nkp, Ef, occ)

    use datatypes
    use numbers
    use global_module, only: iprint_DM, nspin
    use GenComms,      only: myid

    implicit none

    ! passed variables
    integer,                        intent(in)  :: nbands, nkp
    real(double), dimension(:),     intent(in)  :: electrons
    real(double), dimension(:,:,:), intent(in)  :: eig
    real(double), dimension(:),     intent(out) :: Ef
    real(double), dimension(:,:,:), intent(out) :: occ

    ! local variables
    real(double), dimension(nspin) :: thisElec, lowElec, highElec, testElec
    real(double), dimension(nspin) :: lowEf, highEf, testEf, incEf
    integer,      dimension(nspin) :: ne
    real(double) :: electrons_total, gaussian_width
    integer      :: counter, ibrkt, lband, lkp, iband, ikp, spin
    real(double), parameter :: tolElec = 1.0e-6_double

    if (nspin == 2) then
       electrons_total = electrons(1) + electrons(2)
    else
       electrons_total = two * electrons(1)
    end if

    if (iprint_DM >= 2 .and. (inode == ionode)) then
       if (nspin == 1) then
          write (io_lun, 1) myid, electrons_total
       else
          write (io_lun, 2) myid, electrons(1), electrons(2), electrons_total
       end if
    end if

    ! find the correct bracket trapping Ef
    labspin: do spin = 1, nspin

       if (iprint_DM >= 2 .and. nspin == 2 .and. inode == ionode) &
            write (io_lun, 3) myid, spin
       select case (flag_smear_type)

       case (0) ! Fermi smearing

          ! Take first guess as double filling each band at first k
          ! point. Note that electrons(spin) stores number of electrons
          ! in each spin channel
          ne(spin) = int(electrons(spin))
          if (ne(spin) < 1) ne(spin) = 1
          Ef(spin) = eig(ne(spin),1,spin)

          ! occupy for a single spin channel
          call occupy(occ, eig, Ef, thisElec, nbands, nkp, spin=spin)
          ! Find two values than bracket true Ef
          incEf(spin) = one

          if (thisElec(spin) < electrons(spin)) then ! found a lower bound

             if (iprint_DM >= 4 .and. (inode == ionode)) &
                  write (io_lun, 4) myid, Ef(spin)
             lowEf(spin) = Ef(spin)
             lowElec(spin) = thisElec(spin)
             highEf(spin) = lowEf(spin) + incEf(spin)
             ibrkt = 1
             call occupy(occ, eig, highEf, highElec, nbands, nkp, spin=spin)
             do while (highElec(spin) < electrons(spin))
                ! reset if necessary
                if (ibrkt == max_brkt_iterations) then
                   ibrkt = 0
                   highEf(spin) = Ef(spin)  ! start from begining
                   highElec(spin) = thisElec(spin)
                   ! half the searching step-size
                   incEf(spin) = half * incEf(spin)
                end if
                ! increase upper bound
                lowEf(spin) = highEf(spin)
                lowElec(spin) = highElec(spin)
                highEf(spin) = highEf(spin) + incEf(spin)
                call occupy(occ, eig, highEf, highElec, nbands, nkp, spin=spin)
                if (iprint_DM >= 4 .and. inode == ionode) &
                     write (io_lun, 5) myid, highEf(spin), highElec(spin)
                ibrkt = ibrkt + 1
             end do ! while (highElec(spin) < electrons(spin))

          else ! found an upper bound

             if (iprint_DM >= 4 .AND. inode == ionode) &
                  write (io_lun, 6) myid, Ef(spin)
             highEf(spin) = Ef(spin)
             highElec(spin) = thisElec(spin)
             lowEf(spin) = highEf(spin) - incEf(spin)
             ibrkt = 1
             call occupy(occ, eig, lowEf, lowElec, nbands, nkp, spin=spin)
             do while (lowElec(spin) > electrons(spin))
                ! reset if necessary
                if (ibrkt == max_brkt_iterations) then
                   ibrkt = 0
                   lowEf(spin) = Ef(spin)  ! start from begining
                   lowElec(spin) = thisElec(spin)
                   ! half the searching step-size
                   incEf(spin) = half * incEf(spin)
                end if
                ! decrease lower bound
                highEf(spin) = lowEf(spin)
                highElec(spin) = lowElec(spin)
                lowEf(spin) = lowEf(spin) - incEf(spin)
                call occupy(occ, eig, lowEf, lowElec, nbands, nkp, spin=spin)
                if (iprint_DM >= 4 .AND. inode == ionode) &
                     write (io_lun, 5) myid, lowEf(spin), lowElec(spin)
                ibrkt = ibrkt + 1
             end do

          end if ! if (thisElec(spin) < electrons(spin))

          if (iprint_DM > 3 .and. inode == ionode) &
               write (io_lun, 12) myid, lowEf(spin), highEf(spin)

       case (1) ! Methfessel-Paxton smearing

          ! Fill the bands for the first (electrons-NElec_less) electrons
          if (NElec_less >= electrons(spin)) then
             if (inode == ionode) write (io_lun, 7) myid, spin
             NELec_less = electrons(spin)
          end if
          thisElec(spin) = zero
          band1: do iband = 1, nbands
             kpoint1: do ikp = 1, nkp
                thisElec(spin) = thisElec(spin) + wtk(ikp)
                if (thisElec(spin) >= (electrons(spin) - NElec_less)) then
                   lband = iband
                   lkp = ikp
                   exit band1
                end if
             end do kpoint1
          end do band1
          lowEf(spin) = eig(lband,lkp,spin)
          call occupy(occ, eig, lowEf, lowElec, nbands, nkp, spin=spin)
          ! check if we indeed have a good lower bound
          if ((electrons(spin) - lowElec(spin)) < one) then
             if (inode == ionode) &
                  write (io_lun, 8) myid, lowElec(spin), spin, NElec_less
             ! find the lowest energy and start from there
             lband = 1
             lkp = 1
             band : do iband = 1, nbands
                kpoint : do ikp = 1, nkp
                   if (eig(lband,lkp,spin) > eig(iband,ikp,spin)) then
                      lband = iband
                      lkp = ikp
                   end if
                end do kpoint
             end do band
             lowEf = eig(lband,lkp,spin)
             call occupy(occ, eig, lowEf, lowElec, nbands, nkp, spin=spin)
          end if
          ! now that we have a lower-bound, find upper bound
          ! get gaussian width
          if (gaussian_height >= one) then
             if (myid == 0) write (io_lun, 9) myid
             gaussian_height = 0.1_double
          end if
          gaussian_width = two * sqrt(-log(gaussian_height)) * kT
          if (iMethfessel_Paxton < 2) then   !! prevents division by zero using default iMethfessel value
             incEf(spin) = gaussian_width / &
                  (two * finess)
          else 
             incEf(spin) = gaussian_width / &
                  (two * real(iMethfessel_Paxton, double) * finess)
          end if
          highEf(spin) = lowEf(spin) + incEf(spin)
          call occupy(occ, eig, highEf, highElec, nbands, nkp, spin=spin)
          do while (highElec(spin) < electrons(spin))
             ! find upper bound
             lowEf(spin) = highEf(spin)
             lowElec(spin) = highElec(spin)
             highEf(spin) = lowEf(spin) + incEf(spin)
             call occupy(occ, eig, highEf, highElec, nbands, nkp, spin=spin)
             if (iprint_DM >= 4 .and. inode == ionode) &
                  write (io_lun, 5) myid, highEf(spin), highElec(spin)
          end do

          if (iprint_DM > 3 .and. inode == ionode) &
               write (io_lun, 12) myid, lowEf(spin), highEf(spin)

       case default
          call cq_abort ("FindEvals: Smearing flag not recognised",&
               flag_smear_type)
       end select

       ! Starting Bisection
       Ef(spin) = half * (lowEf(spin) + highEf(spin))
       call occupy(occ, eig, Ef, thisElec, nbands, nkp, spin=spin)
       counter = 0
       do while ((abs(thisElec(spin) - electrons(spin))) > tolElec .and. &
            (counter <= maxefermi))
          counter = counter + 1
          if (thisElec(spin) > electrons(spin)) then
             highElec(spin) = thisElec(spin)
             highEf(spin) = Ef(spin)
          else
             lowElec(spin) = thisElec(spin)
             lowEf(spin) = Ef(spin)
          end if
          Ef(spin) = half * (lowEf(spin) + highEf(spin))
          call occupy(occ, eig, Ef, thisElec, nbands, nkp, spin=spin)
       end do

       if (iprint_DM >= 2 .AND. inode == ionode) then
          if (nspin == 1) then
             write (io_lun, 10) Ef(spin)
          else
             write (io_lun, 11) spin, Ef(spin)
          end if
       end if

    end do labspin

    return

1   format(10x, 'Proc: ', i5, ' findFermi_fixspin: searching for Ne: ', f12.5)
2   format(10x, 'Proc: ', i5, &
         ' findFermi_fixspin: searching for Ne (up, dn, total): ', 3f12.5)
3   format(10x, 'Proc: ', i5, ' findFermi_fixspin: Finding Ef for spin = ', i2)
4   format(10x, 'Proc: ', i5, ' findFermi_fixspin: found lower bound', f12.5)
5   format(10x, 'Proc: ', i5, ' findFermi_fixspin: level, Ne: ', 2f12.5)
6   format(10x, 'Proc: ', i5, ' findFermi_fixspin: found upper bound', f12.5)
7   format(10x, 'Proc: ', i5, ' findFermi_fixspin: Warning! Diag.NElecLess >= &
         &total number of electrons for spin channel ', i2, &
         ' setting it equal to number of electrons, but this is slow &
         &and you may want to change it to something smaller.')
8   format(10x, 'Proc: ', i5, ' findFermi_fixspin: Warning! the &
         &calculated number of electrons (',f12.5, &
         ') > electron_number (for spin ', i2, ' ) - 1.0. May be you &
         &should increase the value of Diag.NElecLess (at the moment =&
         & ',f12.5,')')
9   format(10x, 'Proc: ', i5, ' findFermi_fixspin: Warning! Diag.gaussianHeight &
         &must be less than one, reset to 0.1 as default')
10  format(10x, 'Fermi level is ', f12.5)
11  format(10x, 'Fermi level for spin ', i2, ' is ', f12.5)
12  format(10x, 'Proc: ', i5, ' bracketed Ef: ', 2f12.5)

  end subroutine findFermi_fixspin
  !!*****


  !!****f* DiagModule/findFermi_varspin
  !! PURPOSE
  !!   Finds Fermi energy for the case when the populations in the
  !!   different spin channels are allowed to vary. In this case Fermi
  !!   levels in both spin channels are equal.
  !! USAGE
  !!   call findFermi_varspin(electrons_total, eig, nbands, nkp, Ef, occ)
  !! INPUTS
  !!   integer nbands                     : number of bands
  !!   integer nkp                        : number of kpoints
  !!   real(double) eig(nbands,nkp,nspin) : eigenvalues
  !!   real(double) electrons_total       : total number of electrons
  !!                                        which the fermi energy is
  !!                                        required to get.
  !! OUTPUT
  !!   real(double) Ef(nspin)             : Fermi energy for each spin channel
  !!   real(double) occ(nbands,nkp,nspin) : occupancies
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/06
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine findFermi_varspin(electrons_total, eig, nbands, nkp, Ef, occ)

    use datatypes
    use numbers
    use global_module, only: iprint_DM, nspin, spin_factor
    use GenComms,      only: myid

    implicit none

    ! passed variables
    integer,                        intent(in)  :: nbands, nkp
    real(double),                   intent(in)  :: electrons_total
    real(double), dimension(:,:,:), intent(in)  :: eig
    real(double), dimension(:),     intent(out) :: Ef
    real(double), dimension(:,:,:), intent(out) :: occ

    ! local variables
    real(double), dimension(nspin) :: electrons
    real(double), dimension(nspin) :: lowEf, highEf, incEf
    real(double) :: gaussian_width, thisElec, lowElec, highElec
    integer      :: counter, ne, ibrkt, lband, lkp, iband, ikp, spin, lspin
    real(double), parameter :: tolElec = 1.0e-6_double

    ! Finding the correct bracket trapping Ef
    select case (flag_smear_type)

    case (0) ! Fermi smearing

       if (iprint_DM >= 2 .AND. inode == ionode) &
            write (io_lun, 1) myid, electrons_total
       ! Take first guess as double filling each band at first k point
       ne = int(electrons_total / two)
       if (ne < 1) ne = 1
       ! choose Ef to be the minimum of both spin channels
       ! Ef will be the same for all spin channels for variable spin
       do spin = 1, nspin
          Ef(spin) = minval(eig(ne,1,:))
       end do
       call occupy(occ, eig, Ef, electrons, nbands, nkp)
       thisElec = spin_factor * sum(electrons(:))
       ! Find two values than bracket true Ef
       incEf(:) = one
       if (thisElec < electrons_total) then ! found lower bound
          if (iprint_DM >= 4 .and. inode == ionode) &
               write (io_lun, 2) myid, Ef(1) ! Ef(1) = Ef(2) always
          lowEf(:) = Ef(:)
          lowElec = thisElec
          highEf(:) = lowEf(:) + incEf(:)
          ibrkt = 1
          call occupy(occ, eig, highEf, electrons, nbands, nkp)
          highElec = spin_factor * sum(electrons(:))
          do while (highElec < electrons_total) ! Increase upper bound
             if (ibrkt == max_brkt_iterations) then
                ibrkt = 0
                highEf(:) = Ef(:)  ! start from begining
                highElec = thisElec
                ! half the searching step-size
                incEf(:) = half * incEf(:)
             end if
             lowEf(:) = highEf(:)
             lowElec = highElec
             highEf(:) = highEf(:) + incEf(:)
             call occupy(occ, eig, highEf, electrons, nbands, nkp)
             highElec = spin_factor * sum(electrons(:))
             if (iprint_DM >= 4 .and. inode == ionode) &
                  write (io_lun, 3) myid, highEf(1), highElec
             ibrkt = ibrkt + 1
          end do
       else ! found upper bound
          if (iprint_DM >= 4 .and. inode == ionode) &
               write (io_lun, 6) myid, Ef(1)
          highEf(:) = Ef(:)
          highElec = thisElec
          lowEf(:) = highEf(:) - incEf(:)
          ibrkt = 1
          call occupy(occ, eig, lowEf, electrons, nbands, nkp)
          lowElec = spin_factor * sum(electrons(:))
          do while (lowElec > electrons_total) ! Decrease lower bound
             if (ibrkt == max_brkt_iterations) then
                ibrkt = 0
                lowEf(:) = Ef(:) ! start from begining
                lowElec = thisElec
                ! half the searching step-size
                incEf(:) = half * incEf(:)
             end if
             highEf(:) = lowEf(:)
             highElec = lowElec
             lowEf(:) = lowEf(:) - incEf(:)
             call occupy(occ, eig, lowEf, electrons, nbands, nkp)
             lowElec = spin_factor * sum(electrons(:))
             if (iprint_DM >= 4 .and. inode == ionode) &
                  write (io_lun, 3) myid, lowEf(1), lowElec
             ibrkt = ibrkt + 1
          end do
       end if
       if (iprint_DM > 3 .AND. inode == ionode) &
            write (io_lun, 5) myid, lowEf(1), highEf(1)
       ! search method to use in the case of Methmessel Paxton smearing

    case (1) ! Methfessel-Paxton smearing

       ! Fill the bands for the first (electrons_toal - NElec_less) electrons
       if (NElec_less >= electrons_total) then
          if (inode == ionode) write (io_lun, 6)
          NELec_less = electrons_total
       end if
       thisElec = zero
       band1 : do iband = 1, nbands
          kpoint1 : do ikp = 1, nkp
             thisElec = thisElec + two * wtk(ikp)
             if (thisElec >= (electrons_total - NElec_less)) then
                lband = iband
                lkp = ikp
                exit band1
             end if
          end do kpoint1
       end do band1
       lowEf(:) = min(eig(lband,lkp,1), eig(lband,lkp,2))
       call occupy(occ, eig, lowEf, electrons, nbands, nkp)
       lowElec = spin_factor * sum(electrons(:))
       ! check if we indeed have a good lower bound
       if ((electrons_total - lowElec) < two) then
          if (inode == ionode) write (io_lun, 8) lowElec, NElec_less
          ! find the lowest energy and start from there
          lband = 1
          lkp = 1
          lspin = 1
          do spin = 1, nspin
             do iband = 1, nbands
                do ikp = 1, nkp
                   if (eig(lband,lkp,lspin) > eig(iband,ikp,spin)) then
                      lband = iband
                      lkp = ikp
                      lspin = spin
                   end if
                end do
             end do
          end do
          ! get the overall minimum
          lowEf(:) = eig(lband,lkp,lspin)
          call occupy(occ, eig, lowEf, electrons, nbands, nkp)
          lowElec = spin_factor * sum(electrons)
       end if
       ! now that we have a lower-bound, find upper bound
       ! get gaussian width
       if (gaussian_height >= one) then
          if (inode == ionode) write (io_lun, 7)
          gaussian_height = 0.1_double
       end if
       gaussian_width = two * sqrt(-log(gaussian_height)) * kT
       incEf(:) = gaussian_width / &
            (two * real(iMethfessel_Paxton, double) * finess)
       highEf(:) = lowEf(:) + incEf(:)
       call occupy(occ, eig, highEf, electrons, nbands, nkp)
       highElec = spin_factor * sum(electrons(:))
       do while (highElec < electrons_total) ! find upperbound
          lowEf(:) = highEf(:)
          lowElec = highElec
          highEf(:) = lowEf(:) + incEf(:)
          call occupy(occ, eig, highEf, electrons, nbands, nkp)
          highElec = spin_factor * sum(electrons(:))
          if (iprint_DM >= 4 .and. inode == ionode) &
               write (io_lun, 3) myid, highEf(1), highElec
       end do

    case default
       call cq_abort ("FindEvals: Smearing flag not recognised",&
            flag_smear_type)
    end select

    ! Starting Bisection
    Ef(:) = half * (lowEf(:) + highEf(:))
    call occupy(occ, eig, Ef, electrons, nbands, nkp)
    thisElec = spin_factor * sum(electrons(:))
    counter = 0
    do while ((abs(thisElec - electrons_total)) > tolElec .and. &
         (counter <= maxefermi))
       counter = counter + 1
       if (thisElec > electrons_total) then
          highElec = thisElec
          highEf(:) = Ef(:)
       else
          lowElec = thisElec
          lowEf(:) = Ef(:)
       end if
       Ef(:) = half * (lowEf(:) + highEf(:))
       call occupy(occ, eig, Ef, electrons, nbands, nkp)
       thisElec = spin_factor * sum(electrons(:))
    end do

    if (iprint_DM >= 2 .and. inode == ionode) write (io_lun, 8) Ef(1)

    return

1   format(10x, 'Proc: ', i5, ' findFermi_varspin: searching for Ne: ', f12.5)
2   format(10x, 'Proc: ', i5, ' findFermi_varspin: found lower bound ', f12.5)
3   format(10x, 'Proc: ', i5, ' findFermi_varspin: level, Ne: ', 2f12.5)
4   format(10x, 'Proc: ', i5, ' findFermi_varspin: found upper bound ', f12.5)
5   format(10x, 'Proc: ', i5, ' bracketed Ef: ', 2f12.5)
6   format(10x, 'In findFermi, Warning! Diag.NElecLess >= total number &
         &of electrons, setting it equal to number of electrons, but &
         &this is slow and you may want to change it to something &
         &smaller.')
7   format(10x, 'In findFermi, Warning! Diag.gaussianHeight must be &
         &less than one, reset to 0.1 as default')
8   format(10x, 'Fermi level is ', f12.5)

  end subroutine findFermi_varspin
  !!*****


  ! -----------------------------------------------------------------------------
  ! Subroutine occupy
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/occupy *
  !!
  !!  NAME
  !!   occupy - occupy eigenstates
  !!  USAGE
  !!   occupy()
  !!  PURPOSE
  !!   Populates the eigenstates up to a given fermi level at each
  !!   k-point.  Called by findFermi to find the fermi level.
  !!  INPUTS
  !!   real(double), dimension(nbands,nkp) :: occ - occupancies
  !!   real(double), dimension(nbands,nkp) :: ebands - eigenvalues
  !!   real(double) :: Ef - fermi energy
  !!   integer :: nbands, nkp - numbers of bands and k-points
  !!  USES
  !!   datatypes, numbers, global_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   23/04/2002
  !!  MODIFICATION HISTORY
  !!   2010/06/14 23:25 L.Tong
  !!    Added option for using Methfessel-Paxton approximation for
  !!    step-function
  !!   2010/07/26 L.Tong
  !!    Realised occ may be overlapping with the module variable occ,
  !!    hence change the name to occu
  !!   2011/11/30 L.Tong
  !!    Added spin polarisation
  !!     - Added optional parameter spin, and spin = 1 for up and 2 for
  !!       down
  !!     - If spin is present then only the occupancy for the (single)
  !!       given spin channel is calculated.
  !!   2012/03/04 L.Tong
  !!    Major rewrite of the spin implementation.
  !!    - occu(nbands,nkp,nspin), ebands(nbands,nkp,nspin), Ef(nspin),
  !!      electrons(nspin) should only store the corresponding values in
  !!      a single spin channel. This is also the case for spin
  !!      non-polarised cases.
  !!    - If optional parameter spin is present (spin = 1|2), then only
  !!      a single spin component is calculated
  !!  SOURCE
  !!
  subroutine occupy(occu, ebands, Ef, electrons, nbands, nkp, spin)
    use datatypes
    use numbers
    use global_module, only: iprint_DM, nspin
    use GenComms,      only: myid

    implicit none

    ! passed variables
    integer,                        intent(in)  :: nbands, nkp
    real(double), dimension(:),     intent(in)  :: Ef
    real(double), dimension(:,:,:), intent(in)  :: ebands
    real(double), dimension(:),     intent(out) :: electrons
    real(double), dimension(:,:,:), intent(out) :: occu
    integer,      optional,         intent(in)  :: spin

    ! local variables
    integer :: ikp, iband, ss
    integer :: ss_start, ss_end

    if (nspin == 2 .and. present(spin)) then
       ss_start = spin
       ss_end = spin
    else
       ss_start = 1
       ss_end = nspin
    end if

    electrons = zero
    labspin: do ss = ss_start, ss_end
       kp: do ikp = 1, nkp
          band: do iband = 1, nbands
             select case (flag_smear_type)
             case (0) ! Fermi smearing
                occu(iband,ikp,ss) = &
                     wtk(ikp) * fermi(ebands(iband,ikp,ss) - Ef(ss), kT)
             case (1) ! Methfessel Paxton smearing
                occu(iband,ikp,ss) = &
                     wtk(ikp) * MP_step(ebands(iband,ikp,ss) - Ef(ss), &
                     iMethfessel_Paxton, kT)
             case default
                call cq_abort ("FindEvals: Smearing flag not recognised",&
                     flag_smear_type)
             end select
             electrons(ss) = electrons(ss) + occu(iband,ikp,ss)
          end do band
       end do kp
    end do labspin

    if (iprint_DM >=5 .and. (inode == ionode)) then
       if (nspin == 1) then
          write (io_lun, 1) myid, Ef(1), two * electrons(1)
       else if (present (spin)) then
          write (io_lun, 2) myid, spin, Ef(spin), electrons(spin)
       else
          write (io_lun, 3) myid, Ef(1), Ef(2), &
               electrons(1), electrons(2), &
               electrons(1) + electrons(2)
       end if
    end if

    return

1   format(10x, 'In occupy on proc: ', i5, &
         ' For Ef of ', f8.5,            &
         ' we get ', f12.5, ' electrons')
2   format(10x, 'In occupy on proc: ', i5,           &
         ' For Spin = ', i1, ' with Ef of ', f8.5, &
         ' we get ', f12.5, ' electrons')
3   format(10x, 'In occupy on proc: ', i5,                   &
         ' For Ef (up, down) of (', f8.5, ', ', f8.5, ')', &
         ' we get electrons (up, down, total)', 3f12.5)

  end subroutine occupy
  !!***


  ! -----------------------------------------------------------------------------
  ! Function fermi
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/fermi *
  !!
  !!  NAME
  !!   fermi - evaluate fermi function
  !!  USAGE
  !!   fermi(E,kT)
  !!  PURPOSE
  !!   Evaluates the fermi occupation of an energy
  !!
  !!   I'm assuming (for the sake of argument) that if both the energy and
  !!   the smearing (kT) are zero then we get an occupation of 0.5 - this is
  !!   certainly the limit if E and kT are equal and heading to zero, or if
  !!   E is smaller than kT and both head for zero.
  !!  INPUTS
  !!   real(double), intent(in) :: E - energy
  !!   real(double), intent(in) :: kT - smearing energy
  !!  USES
  !!   datatypes, numbers
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   23/04/2002
  !!  MODIFICATION HISTORY
  !!   2006/10/02 17:54 dave
  !!    Small fix to prevent maths overflows by only calculating
  !!    exponential if x well bounded
  !!   2012/01/22 L.Tong
  !!    Small change to use FORTRAN 90 function declaration notation
  !!    This works better with etags
  !!  SOURCE
  !!
  function fermi(E,kT)

    use datatypes
    use numbers, only: zero, one, half

    implicit none

    ! result
    real(double) :: fermi

    ! Passed variables
    real(double), intent(in) :: E
    real(double), intent(in) :: kT

    ! Local variables
    real(double) :: x
    real(double), parameter :: cutoff = 10.0_double

    if(kT==zero) then
       if(E>zero) then
          fermi = zero
       else if(E<zero) then
          fermi = one
       else if(E==zero) then
          fermi = half
       end if
    else
       x = E/kT
       if(x > cutoff) then
          fermi = zero
       elseif(x < -cutoff) then
          fermi = one
       else
          fermi = one/(one + exp(x))
       endif
    end if
  end function fermi
  !!***


  ! -----------------------------------------------------------------------------
  ! Function MP_step
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/MP_step *
  !!
  !!  NAME
  !!   MP_step - evaluate Methfessel-Paxton step function
  !!  USAGE
  !!   MP_step(E,order,smear)
  !!  PURPOSE
  !!   Evaluates the order (order) Methfessel-Paxton approximation to
  !!   step function
  !!  INPUTS
  !!   real(double), intent(in) :: E - energy
  !!   integer, intent(in) :: order - order of Methfessel expansion
  !!   real(double), intent(in) :: smear - smearing energy, nothing to
  !!                               do with physical temperature
  !!  USES
  !!   datatypes, numbers
  !!  AUTHOR
  !!   L.Tong (lt)
  !!  CREATION DATE
  !!   2010/06/15 00:17
  !!  MODIFICATION HISTORY
  !!   2012/01/22 L.Tong
  !!     - Small change to use FORTRAN 90 function declaration notation
  !!       This works better with etags
  !!  SOURCE
  !!
  function MP_step(E,order,smear)

    use datatypes
    use numbers, only: zero, one, half, two, four, pi
    use functions, ONLY: erfc_cq

    implicit none

    ! Result
    real(double) :: MP_step

    ! Passed variables
    real(double), intent(in) :: E
    integer,      intent(in) :: order
    real(double), intent(in) :: smear

    ! Internal variables
    real(double) :: x, A, H0, H1, H2, nd, x2
    integer      :: n

    ! in case of smear==0, we have the exact step function
    if(smear==zero) then
       if(E>zero) then
          MP_step = zero
       else if(E<zero) then
          MP_step = one
       else if(E==zero) then
          MP_step = half
       end if
    else if(smear>zero) then
       x = E/smear
       if(order==0) then
          MP_step = half*erfc_cq(x)
       else
          x2 = x*x
          A = one/sqrt(pi)
          H0 = one
          H1 = two*x
          MP_step = half*erfc_cq(x)
          nd = one
          do n=1,order
             A = A/((-four)*real(n,double))
             MP_step = MP_Step + A*H1*exp(-x2)
             H2 = two*x*H1 - two*nd*H0
             H0 = H1
             H1 = H2
             nd = nd + one
             H2 = two*x*H1 - two*nd*H0
             H0 = H1
             H1 = H2
             nd = nd + one
          end do
       end if
    end if
  end function MP_step
  !!***

  ! -----------------------------------------------------------------------------
  ! Function MP_entropy
  ! -----------------------------------------------------------------------------

  !!****f* DiagModule/MP_entropy *
  !!
  !!  NAME
  !!   MP_entropy - evaluate the function SN(x) = 0.5*A_N*H_2N(x)*exp(-x^2)
  !!                where A_N = (-1)^N / (N!*4^N*sqrt(PI))
  !!  USAGE
  !!   MP_entropy(x,order)
  !!  PURPOSE
  !!   Evaluate the function SN(x) = 0.5*A_N*H_2N(x)*exp(-x^2)
  !!                where A_N = (-1)^N / (N!*4^N*sqrt(PI))
  !!  INPUTS
  !!   real(double), intent(in) :: x
  !!   integer, intent(in) :: order - order of Methfessel expansion
  !!  USES
  !!   datatypes, numbers
  !!  AUTHOR
  !!   L.Tong (lt)
  !!  CREATION DATE
  !!   2010/07/21 13:50
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  function MP_entropy (x, order)

    use datatypes
    use numbers, only: one, two, half, four, pi

    implicit none

    ! Passed variables
    real(double), intent(in) :: x
    integer,      intent(in) :: order
    ! Result
    real(double) :: MP_entropy
    ! Internal variables
    real(double) :: A, H1, H2, H3, nd
    integer      :: n

    if(order==0) then
       MP_entropy = exp(-(x*x))/sqrt(pi)
    else
       ! Evaluate A_n and Hermite Polynomial order 2n
       A = one/(sqrt(pi)*(-four))
       H1 = two*x
       H2 = two*x*H1 - two
       nd = two
       do n=2,order
          A = A/((-four)*real(n,double))
          H3 = two*x*H2 - two*nd*H1
          H1 = H2
          H2 = H3
          nd = nd + one
          H3 = two*x*H2 - two*nd*H1
          H1 = H2
          H2 = H3
          nd = nd + one
       end do
       MP_entropy = half*A*H2*exp(-(x*x))
    end if
  end function MP_entropy
  !!***


  !!****f* DiagModule/buildK *
  !!
  !!  NAME
  !!   buildK - makes K
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Builds K from eigenvectors - this involves working out which
  !!   processors we're going to need which eigenvectors from, fetching
  !!   the data and building the matrix
  !!
  !!   N.B. The conjugation of one set of eigenvector coefficients takes
  !!   place when the dot product is performed - dot maps onto zdotc
  !!   which conjugates the FIRST vector.  Really we should conjugate
  !!   the second, but as K is real, it shouldn't matter.
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   24/04/2002
  !!  MODIFICATION HISTORY
  !!   01/05/2002 dave
  !!    This routine is now working, so tidied - added deallocate
  !!    statements, iprint_DM levels to write statements and deleted
  !!    unnecessary rubbish
  !!   2004/11/10 drb
  !!    Changed nsf to come from maxima, not common
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2006/10/02 17:52 dave
  !!    Added deallocate for norb_send (thanks TM, TO)
  !!   2012/01/23 L.Tong
  !!    - Added modification for spin polarisation. This mainly to fix the
  !!      fact that occ for spin polarised case will be always from 0 to 1
  !!      and hence the factor of 0.5 is not required.
  !!    - used numbers module
  !!   2012/03/07 L.Tong
  !!    - Changed implementation in spin polarisation (and spin
  !!      non-polarised case). Now occ(spin) should only contain
  !!      occupancies of single spin channels, no matter the calculation
  !!      is spin polarised or not. Hence occ_correction should be one
  !!      in all cases.
  !!   2013/02/13 17:10 dave
  !!    - Changed length of sent array to include all bands up to last occupied
  !!   2015/06/05 16:43 dave
  !!    Added code to calculate K matrix for individual bands (to allow output of charge density from bands)
  !!   2018/10/08 16:38 dave
  !!    Changing MPI tags to conform to MPI standard
  !!  SOURCE
  !!
  subroutine buildK(range, matA, occs, kps, weight, localEig,matBand,locw)

    !use maxima_module, only: mx_nponn, mx_at_prim
    use numbers
    use matrix_module,   only: matrix, matrix_halo
    use group_module,    only: parts
    use primary_module,  only: bundle
    use cover_module,    only: BCS_parts
    use ScalapackFormat, only: CC_to_SC,maxrow,maxcol,proc_block,    &
         SC_to_refx,SC_to_refy, block_size_r,  &
         block_size_c, blocks_c, proc_start,   &
         matrix_size
    use global_module,   only: numprocs, iprint_DM, id_glob,         &
         ni_in_cell, x_atom_cell, y_atom_cell, &
         z_atom_cell, max_wf, out_wf
    use mpi
    use GenBlas,         only: dot
    use GenComms,        only: myid
    use mult_module,     only: store_matrix_value_pos, matrix_pos
    use matrix_data,     only: mat, halo
    use species_module,  only: nsf_species

    implicit none

    ! Passed variables
    real(double), dimension(:), intent(in) :: occs
    real(double), dimension(3), intent(in) :: kps
    real(double) :: weight
    integer :: matA, range
    complex(double_cplx), dimension(:,:), intent(in) :: localEig
    integer, OPTIONAL, dimension(max_wf) :: matBand
    real(double), OPTIONAL, dimension(matrix_size) :: locw
    ! Local variables
    type(Krecv_data), dimension(:), allocatable :: recv_info
    integer :: part, memb, neigh, ist, prim_atom, owning_proc, locatom
    integer :: Row_FSC_part, Row_FSC_seq, Col_FSC_part, Col_FSC_seq, &
         FSC_atom
    integer :: SCblockr, SCblockc, SCrowc, i, j, k, proc, stat, &
         supfn_r, supfn_c
    integer :: maxloc, maxint, maxsend, curr, gcspart, CC, orb_count
    integer :: len, send_size, recv_size, send_proc, recv_proc, nsf1, &
         sendtag, recvtag
    integer :: req1, req2, ierr, atom, inter, prim, wheremat, row_sup,&
         col_sup
    integer, dimension(:,:), allocatable :: ints, atom_list, &
         send_prim, send_info, send_orbs, send_off
    integer, dimension(:), allocatable :: current_loc_atoms, &
         LocalAtom, num_send, norb_send, send_FSC, recv_to_FSC, &
         mapchunk, prim_orbs
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    real(double) :: phase, rfac, ifac, rcc, icc, rsum
    complex(double_cplx) :: zsum
    complex(double_cplx), dimension(:,:), allocatable :: RecvBuffer, &
         SendBuffer
    logical :: flag, flag_write_out
    integer :: FSCpart, ipart, iband, iwf, len_occ
    ! for spin polarisation
    real(double) :: occ_correction

    call start_timer(tmr_std_matrices)
    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,"Entering &
         &buildK ",i4)') matA

    if(PRESENT(matBand)) then
       flag_write_out = .true.
    else
       flag_write_out = .false.
    end if
    ! get occ_correction
    occ_correction = one

    ! Allocate data and zero arrays
    allocate(ints(numprocs,bundle%mx_iprim),&
         current_loc_atoms(numprocs),atom_list(numprocs,bundle%&
         mx_iprim), LocalAtom(ni_in_cell),send_prim(numprocs,bundle%&
         mx_iprim), num_send(numprocs),norb_send(numprocs),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating ints, &
         &current_loc_atoms and atom_list !',stat)
    ints = 0
    current_loc_atoms = 0
    atom_list = 0
    LocalAtom = 0
    send_prim = 0
    num_send = 0
    norb_send = 0
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'buildK: Stage one'
    ! Step one - work out which processors we need to exchange data with
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(iprint_DM>=5.AND.myid==0) write(io_lun,1) myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             if(iprint_DM>=5.AND.myid==0) write(io_lun,2) myid,memb
             prim_atom = bundle%nm_nodbeg(part)+memb-1
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist))
                Col_FSC_seq  = mat(part,range)%i_seq(ist)
                ! Now use i_cc2node(Col_FSC_part) to establish which processor owns this atom
                owning_proc = parts%i_cc2node(Col_FSC_part)
                ! Find Fundamental Simulation Cell atom
                FSC_atom = id_glob(parts%icell_beg(Col_FSC_part)+Col_FSC_seq-1)
                ! Find if we have seen this before
                flag = .false.
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'prim, neigh, FSC: ',prim_atom, neigh, FSC_atom
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'curr_loc_atoms: ',current_loc_atoms(owning_proc)
                if(current_loc_atoms(owning_proc)>0) then
                   do i=1,current_loc_atoms(owning_proc)
                      if(atom_list(owning_proc,i)==FSC_atom) then
                         if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'Loc atom: ',i, LocalAtom(FSC_atom)
                         ints(owning_proc,LocalAtom(FSC_atom)) = ints(owning_proc,LocalAtom(FSC_atom)) + 1
                         send_prim(owning_proc,prim_atom) = nsf_species(bundle%species(prim_atom))
                         flag = .true.
                         exit
                      end if
                   end do
                end if ! current_loc_atoms(owning_proc)>0
                if(flag) then
                   cycle
                end if
                ! Record
                current_loc_atoms(owning_proc) = current_loc_atoms(owning_proc) + 1
                atom_list(owning_proc,current_loc_atoms(owning_proc)) = FSC_atom
                LocalAtom(FSC_atom) = current_loc_atoms(owning_proc)
                ints(owning_proc,LocalAtom(FSC_atom)) = ints(owning_proc,LocalAtom(FSC_atom)) + 1
                send_prim(owning_proc,prim_atom) = nsf_species(bundle%species(prim_atom))
             end do ! End do neigh=1,mat%n_nab
          end do ! End do memb =1,nm_nodgroup
       end if ! End if nm_nodgroup > 0
    end do ! End do part=1,groups_on_node
    ! Find max value of current_loc_atoms and interactions
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'buildK: Stage two'
    maxloc = 0
    maxint = 0
    maxsend = 0
    do i=1,numprocs
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) myid,' Curr loc atoms: ',i,current_loc_atoms(i)
       if(current_loc_atoms(i)>maxloc) maxloc = current_loc_atoms(i)
       do j=1,bundle%mx_iprim ! Needs to be mx_iprim because goes over primary atoms on REMOTE processors
          if(ints(i,j)>maxint) maxint = ints(i,j)
          if(send_prim(i,j)>0) num_send(i) = num_send(i) + 1
          norb_send(i) = norb_send(i) + send_prim(i,j)
          if(iprint_DM>=5.AND.myid==0) write(io_lun,4) myid,j,send_prim(i,j),num_send(i)
       end do
       if(num_send(i)>maxsend) maxsend = num_send(i)
    end do
    if(iprint_DM>=4.AND.myid==0) write(io_lun,*) myid,' Maxima: ',maxloc, maxint, maxsend
    ! Allocate recv_info
    allocate(send_info(numprocs,maxsend),send_orbs(numprocs,maxsend),send_off(numprocs,maxsend), &
         prim_orbs(bundle%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating send_info !',stat)
    send_info = 0
    send_orbs = 0
    send_off = 0
    prim_orbs = 0
    orb_count = 0
    do j=1,bundle%n_prim
       prim_orbs(j) = orb_count
       orb_count = orb_count + nsf_species(bundle%species(j))
    end do
    allocate(recv_info(numprocs),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating recv_info !',stat)
    do i=1,numprocs
       ! Build the list of which primary set atoms we send to which processor
       curr = 0
       orb_count = 0
       do j=1,bundle%n_prim
          if(send_prim(i,j)>0) then
             curr = curr+1
             send_info(i,curr)=j
             send_off(i,curr)=orb_count
             send_orbs(i,curr)=nsf_species(bundle%species(j))
          end if
          orb_count = orb_count + nsf_species(bundle%species(j))
       end do
       allocate(recv_info(i)%ints(maxloc),recv_info(i)%ndimj(maxloc), &
            recv_info(i)%prim_atom(maxint,maxloc), recv_info(i)%locj(maxint,maxloc), &
            recv_info(i)%dx(maxint,maxloc),recv_info(i)%dy(maxint,maxloc), recv_info(i)%dz(maxint,maxloc),STAT=stat)
       recv_info(i)%orbs = 0
       recv_info(i)%ints = 0
       !recv_info(i)%ndimi = 0
       recv_info(i)%ndimj = 0
       recv_info(i)%prim_atom = 0
       recv_info(i)%locj = 0
       recv_info(i)%dx = 0.0_double
       recv_info(i)%dy = 0.0_double
       recv_info(i)%dz = 0.0_double
       if(stat/=0) call cq_abort('buildK: Error allocating recv_info !',stat)
    end do
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(iprint_DM>=5.AND.myid==0) write(io_lun,1) myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             if(iprint_DM>=5.AND.myid==0) write(io_lun,2) myid,memb
             prim_atom = bundle%nm_nodbeg(part)+memb-1
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist))
                Col_FSC_seq  = mat(part,range)%i_seq(ist)
                ! Now use i_cc2node(Col_FSC_part) to establish which processor owns this atom
                owning_proc = parts%i_cc2node(Col_FSC_part)
                ! Find Fundamental Simulation Cell atom
                FSC_atom = id_glob(parts%icell_beg(Col_FSC_part)+Col_FSC_seq-1)
                ! Work out a map from primary atom + FSC + identifier to distance and position in data_Matrix
                locatom = LocalAtom(FSC_atom) ! Which atom in the list on the remote proc is this ?
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) myid,' own, FSC, loc: ',owning_proc, FSC_atom, locatom, &
                     recv_info(owning_proc)%ints(locatom)
                recv_info(owning_proc)%ints(locatom) = recv_info(owning_proc)%ints(locatom) + 1
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) myid,' ints: ',recv_info(owning_proc)%ints(locatom)
                gcspart = BCS_parts%icover_ibeg(mat(part,range)%i_part(ist))+mat(part,range)%i_seq(ist)-1
                !recv_info(owning_proc)%ndimi(locatom) = mat(part,range)%ndimi(memb)
                recv_info(owning_proc)%ndimj(locatom) = mat(part,range)%ndimj(ist)
                if(recv_info(owning_proc)%ints(locatom)==1) &
                     recv_info(owning_proc)%orbs = recv_info(owning_proc)%orbs + mat(part,range)%ndimj(ist)
                recv_info(owning_proc)%prim_atom(recv_info(owning_proc)%ints(locatom),locatom) = prim_atom
                recv_info(owning_proc)%locj(recv_info(owning_proc)%ints(locatom),locatom) = halo(range)%i_halo(gcspart)
                ! Build the distances between atoms - needed for phases
                FSCpart = BCS_parts%lab_cell(mat(part,range)%i_part(ist))!gcspart)
                recv_info(owning_proc)%dx(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%xcover(gcspart)-bundle%xprim(prim_atom)
                recv_info(owning_proc)%dy(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%ycover(gcspart)-bundle%yprim(prim_atom)
                recv_info(owning_proc)%dz(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%zcover(gcspart)-bundle%zprim(prim_atom)
             end do ! End do neigh=1,mat%n_nab
          end do ! End do memb =1,nm_nodgroup
       end if ! End if nm_nodgroup > 0
    end do ! End do part=1,groups_on_node
    ! Work out length
    ! NB we send all bands up to last occupied one (for deltaSCF this will include some empty)
    len = 0
    do i=1,matrix_size ! Effectively all bands
       if(abs(occs(i))>RD_ERR) then
          len = i !len+1
          if(myid==0.AND.iprint_DM>=4) write(io_lun,*) 'Occ is ',occs(i)
       end if
    end do
    ! DRB for WF output - means that we transfer the coefficients into the conduction band
    ! but leave the occupancies along
    if(flag_write_out) then
       len_occ = len
       do iwf = 1,max_wf
          iband = out_wf(iwf)
          if(iband>len) len = iband
       end do
       write(io_lun,*) 'Number of bands: ',len_occ,len
    else
       len_occ = len
    end if
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'buildK: Stage three len:',len, matA
    ! Step three - loop over processors, send and recv data and build K
    allocate(send_fsc(bundle%mx_iprim),recv_to_FSC(bundle%mx_iprim),mapchunk(bundle%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating send_fsc, recv_to_FSC and mapchunk',stat)
    send_fsc = 0
    recv_to_FSC = 0
    mapchunk = 0
    send_proc = myid
    recv_proc = myid
    ! Two messages per process; use tag and tag + 1 below
    sendtag = 1
    recvtag = 1
    do i=1,numprocs
       send_size = len*norb_send(send_proc+1)!num_send(send_proc+1)*nsf
       recv_size = len*recv_info(recv_proc+1)%orbs!current_loc_atoms(recv_proc+1)*nsf
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Send and recv sizes: ',send_size, recv_size
       ! Fill SendBuffer
       allocate(SendBuffer(len,norb_send(send_proc+1)),STAT=stat)
       if(stat/=0) call cq_abort('buildK: Unable to allocate SendBuffer !',stat)
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Filling SendBuffer'
       orb_count = 0
       do j=1,num_send(send_proc+1)
          do nsf1=1,send_orbs(send_proc+1,j)
             orb_count = orb_count+1
             SendBuffer(1:len,orb_count) = localEig(1:len,send_off(send_proc+1,j)+nsf1)
          end do
          ! We also need to send a list of what FSC each primary atom sent corresponds to - use bundle%ig_prim
          send_FSC(j) = bundle%ig_prim(send_info(send_proc+1,j))
          if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Building send_FSC: ',send_info(send_proc+1,j), &
               bundle%ig_prim(send_info(send_proc+1,j)),send_FSC(j)
       end do
       if(orb_count/=norb_send(send_proc+1)) call cq_abort("Orbital mismatch in buildK: ",orb_count,norb_send(send_proc+1))
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Sending'
       ! Now send
       if(send_size>0) then
          if(send_proc/=myid) then
             call MPI_issend(send_FSC,num_send(send_proc+1),MPI_INTEGER,send_proc,sendtag,MPI_COMM_WORLD,req1,ierr)
             call MPI_issend(SendBuffer,send_size,MPI_DOUBLE_COMPLEX,send_proc,sendtag+1,MPI_COMM_WORLD,req2,ierr)
          end if
       end if
       ! Now receive data
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Alloc RecvBuffer ',len,recv_info(recv_proc+1)%orbs
       !allocate(RecvBuffer(len,current_loc_atoms(recv_proc+1)*nsf),STAT=stat)
       allocate(RecvBuffer(len,recv_info(recv_proc+1)%orbs),STAT=stat)
       if(stat/=0) call cq_abort('buildK: Unable to allocate RecvBuffer !',stat)
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Recving'
       if(recv_size>0) then
          if(recv_proc/=myid) then
             call MPI_recv(recv_to_FSC,current_loc_atoms(recv_proc+1),MPI_INTEGER,recv_proc,recvtag,MPI_COMM_WORLD,mpi_stat,ierr)
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Got recv_to_FSC'
             call MPI_recv(RecvBuffer,recv_size,MPI_DOUBLE_COMPLEX,&
                  recv_proc,recvtag+1,MPI_COMM_WORLD,mpi_stat,ierr)
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Got RecvBuffer'
          else
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'On-proc: getting recv_to_FSC'
             recv_to_FSC(1:current_loc_atoms(recv_proc+1)) = send_FSC(1:current_loc_atoms(recv_proc+1))
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'On-proc: getting RecvBuffer'
             RecvBuffer(1:len,1:recv_info(recv_proc+1)%orbs) = SendBuffer(1:len,1:recv_info(recv_proc+1)%orbs)
          end if
          if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Doing the mapchunk', recv_to_FSC
          do j=1,current_loc_atoms(recv_proc+1)
             mapchunk(j) = LocalAtom(recv_to_FSC(j))
          end do
          if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'filling buffer'
          do j=1,len_occ ! This is a loop over eigenstates
             RecvBuffer(j,1:recv_info(recv_proc+1)%orbs) = RecvBuffer(j,1:recv_info(recv_proc+1)%orbs)*occ_correction*occs(j)
          end do
          orb_count = 0
          do atom = 1,current_loc_atoms(recv_proc+1)
             locatom = mapchunk(atom)
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Atom, loc: ',atom,locatom,recv_info(recv_proc+1)%ints(locatom)
             ! Scale the eigenvector coefficients we've received
             ! The factor of 0.5 is because the occupation numbers are from 0->2 (we expect 0->1 in K)
             ! The occupation numbers contain the k-point weight
             ! When we're doing M12 not K then occ also contains the eigenvalue
             !do col_sup = 1,nsf
             !   do j=1,len
             !      RecvBuffer(j,(atom-1)*nsf+col_sup) = RecvBuffer(j,(atom-1)*nsf+col_sup)*0.5_double*occs(j)
             !   end do
             !end do
             ! N.B. the routine used for dot is zdotc which takes the complex conjugate of the first vector
             do inter = 1,recv_info(recv_proc+1)%ints(locatom)
                prim = recv_info(recv_proc+1)%prim_atom(inter,locatom)
                if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Inter: ',inter,prim
                phase = kps(1)*recv_info(recv_proc+1)%dx(inter,locatom) + kps(2)*recv_info(recv_proc+1)%dy(inter,locatom) + &
                     kps(3)*recv_info(recv_proc+1)%dz(inter,locatom)
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'Prim, where, phase: ',prim, whereMat, phase
                rfac = cos(phase)
                ifac = sin(phase)
                do row_sup = 1,recv_info(recv_proc+1)%ndimj(locatom)
                   do col_sup = 1,nsf_species(bundle%species(prim))
                      whereMat = matrix_pos(matA,recv_info(recv_proc+1)%prim_atom(inter,locatom), &
                           recv_info(recv_proc+1)%locj(inter,locatom),col_sup,row_sup)
                      zsum = dot(len_occ,localEig(1:len_occ,prim_orbs(prim)+col_sup),1,RecvBuffer(1:len_occ,orb_count+row_sup),1)
                      call store_matrix_value_pos(matA,whereMat,real(zsum*cmplx(rfac,ifac,double_cplx),double))
                      if(flag_write_out) then
                         do iwf = 1,max_wf
                            iband = out_wf(iwf)
                            zsum = conjg(localEig(iband,prim_orbs(prim)+col_sup))*RecvBuffer(iband,orb_count+row_sup)
                            call store_matrix_value_pos(matBand(iwf),whereMat,real(zsum*cmplx(rfac,ifac,double_cplx),double))
                         end do ! iwf = max_wf
                      end if
                   end do ! col_sup=nsf
                end do ! row_sup=nsf
             end do ! inter=recv_info%ints
             ! Careful - we only want to increment after ALL interactions done
             orb_count = orb_count + recv_info(recv_proc+1)%ndimj(locatom)
          end do ! atom=current_loc_atoms
       end if ! recv_size>0
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Calling MPI_Wait'
       if(send_size>0.AND.myid/=send_proc) then
          call MPI_Wait(req1,mpi_stat,ierr)
          call MPI_Wait(req2,mpi_stat,ierr)
       end if
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Calling dealloc'
       deallocate(RecvBuffer,STAT=stat)
       if(stat/=0) call cq_abort("buildK: Failed to dealloc buffer",stat)
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Calling dealloc'
       deallocate(SendBuffer,STAT=stat)
       if(stat/=0) call cq_abort("buildK: Failed to dealloc buffer",stat)
       ! Increment/decrement recv and send, and wrap
       ! Remember that we go from 0->numprocs-1
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Doing proc thang'
       send_proc = send_proc +1
       if(send_proc.GT.numprocs-1) send_proc = 0
       recv_proc = recv_proc -1
       if(recv_proc.LT.0) recv_proc = numprocs-1
    end do ! do i=numprocs
    ! Now deallocate all arrays
    deallocate(send_fsc,recv_to_FSC,mapchunk,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error deallocating send_fsc, recv_to_FSC and mapchunk !',stat)
    do i=numprocs,1,-1
       deallocate(recv_info(i)%ints,recv_info(i)%prim_atom,recv_info(i)%locj,&
            recv_info(i)%dx,recv_info(i)%dy,recv_info(i)%dz,STAT=stat)
       if(stat/=0) call cq_abort('buildK: Error deallocating recvinfo !',i,stat)
    end do
    deallocate(prim_orbs,send_off,send_orbs,send_info,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating send_info !',stat)
    deallocate(recv_info,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating recv_info !',stat)
    deallocate(ints,current_loc_atoms,atom_list,LocalAtom,send_prim,num_send,norb_send,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating ints etc !',stat)
    call stop_timer(tmr_std_matrices)
    return
1   format(10x,'Processor: ',i5,' Partition: ',i5)
2   format(10x,'Processor: ',i5,' Atom: ',i5)
3   format(10x,'Processor: ',i5,' Neighbour: ',i5)
4   format(10x,'Proc: ',i5,' Prim, send_prim, num_send: ',3i5)
  end subroutine buildK
  !!***

  !!****f*  DiagModule/wf_output 
  !!
  !!  NAME 
  !!   wf_output
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Outputs the KS wavefuntion charge
  !!
  !!      
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   C. O'Rourke
  !!  CREATION DATE
  !!   2015/05/29 
  !!  MODIFICATION HISTORY
  !!   2015/06/05 16:42 dave
  !!    Tidied and removed sum over grid points (done in get_band_density)
  !!   2015/07/02 08:21 dave
  !!    Changing to write by k-point (same file but introduce header)
  !!  SOURCE
  !!
  subroutine wf_output(spin,abs_wf,wf_no,kp,energy,first)

    use datatypes
    use GenComms,       ONLY: gsum, my_barrier
    use global_module,  ONLY: out_wf, nspin
    use io_module,      ONLY: dump_band_charge
    use maxima_module,  ONLY: maxngrid
    use numbers,        ONLY: zero

    implicit NONE

    real(double), dimension(:) :: abs_wf
    integer(integ)    :: spin,wf_no
    real(double), optional :: energy
    real(double), optional, dimension(3) :: kp
    integer, optional :: first

    character(len=50) :: ci,cspin

    !ci=adjustl(ci)
    if(PRESENT(kp)) then 
       write(ci,'("Band",I0.6,"K")') out_wf(wf_no)
    else
       write(ci,'("Band",I0.6)') out_wf(wf_no)
    end if
    !ci=trim(ci)
    if(present(kp).AND.present(energy)) then
       if(present(first)) then
          call dump_band_charge(trim(ci)//"wf",abs_wf(:),maxngrid,inode,spin,kp,energy,nkp)
       else
          call dump_band_charge(trim(ci)//"wf",abs_wf(:),maxngrid,inode,spin,kp,energy)
       end if
    else
       call dump_band_charge(trim(ci)//"wf",abs_wf(:),maxngrid,inode,spin)
    end if
  end subroutine wf_output
  !!***

  !!****f*  DiagModule/accumulate_DOS
  !!
  !!  NAME 
  !!   wf_output
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Accumulates DOS
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2016 ?
  !!  MODIFICATION HISTORY
  !!   2017/11/01 18:00 nakata
  !!    Introduced PDOS with MSSFs, projecting on neighbor atoms with global ID.
  !!    Added spinSF to specify the spin of the MSSF coefficients.
  !!    SpinSF is not used for primitive PAOs.
  !!   2017/11/13 18:15 nakata
  !!    Introduced the normalization of each eigenstate of PDOS.
  !!    Added optional argument weight_pDOS, the normalisation weight.
  !!   2018/09/19 18:30 nakata
  !!    Introduced orbital angular momentum resolved DOS.
  !!    Added optional projDOS_angmom, l-projected PDOS.
  !!   2018/10/22 14:18 dave & jsb
  !!    Adding (l,m) projection for pDOS
  !!   2018/10/30 11:43 dave
  !!    Implementing semi-core exclusion given right flags
  !!    (NB at present, semi-core states are not flagged but will be !)
  !!   2018/11/02 16:30 nakata
  !!    Bug fix: changed atom_spec to neigh_species for semicore of neighbour atoms
  !!  SOURCE
  !!
  subroutine accumulate_DOS(weight,eval,evec,DOS,spinSF,projDOS,projDOS_angmom,weight_pDOS)

    use datatypes
    use numbers,         only: half, zero
    use global_module,   only: n_DOS, E_DOS_max, E_DOS_min, flag_write_DOS, sigma_DOS, flag_write_projected_DOS, &
                               sf, atomf, id_glob, species_glob, flag_normalise_pDOS, flag_pDOS_angmom, flag_pDOS_lm, &
                               flag_SpinDependentSF
    use ScalapackFormat, only: matrix_size
    use species_module,  only: nsf_species, natomf_species
    use group_module,    only: parts
    use primary_module,  only: bundle
    use cover_module,    only: BCS_parts
    use matrix_data,     only: mat, halo, SFcoeff_range
    use mult_module,     only: matSFcoeff, matrix_pos, mat_p
    use GenComms,        only: cq_abort
    use pao_format
    
    implicit none

    ! Passed variables
    real(double) :: weight
    complex(double_cplx), dimension(:,:), intent(in) :: evec
    real(double), dimension(:) :: eval
    real(double), dimension(n_DOS) :: DOS
    integer, intent(in) :: spinSF ! used only if atomf/=sf
    real(double), OPTIONAL, dimension(:,:) :: projDOS
    real(double), OPTIONAL, dimension(:,:,:,:) :: projDOS_angmom ! Dimensions are bin, atom, l, m
    real(double), OPTIONAL, dimension(:) :: weight_pDOS

    ! Local variables
    integer :: iwf, n_band, n_min, n_max, i, acc, spin_SF, &
               atom, isf1, nsf1, atom_spec, l1, nacz1, m1, &
               neigh_global_num, iatomf2, natomf2, neigh_species, l2, nacz2, m2, &
               atom_num, gcspart, neigh_global_part, j_in_halo, wheremat
    integer :: iprim, part, memb, neigh, ist
    real(double) :: Ebin, a, fac, fac1, fac2, val
    real(double), dimension(6,13) :: fac_angmom, fac2_angmom ! up to h-orbital and 2l+1
    real(double), dimension(n_DOS) :: tmp

    if(present(projDOS).AND.(.NOT.flag_write_projected_DOS)) call cq_abort("Called pDOS without flag")
    if(present(projDOS_angmom).AND.(.NOT.flag_pDOS_angmom))  call cq_abort("Called pDOS_angmom without flag")
    if(present(weight_pDOS) .AND.(.NOT.flag_normalise_pDOS)) call cq_abort("Normalised pDOS without flag")
    ! ---------------
    ! DOS calculation
    ! ---------------
    ! Now accumulate DOS for this band
    do iwf=1,matrix_size ! Effectively all bands
       tmp = zero
       n_band = floor((eval(iwf) - E_DOS_min)/dE_DOS) + 1
       n_min = n_band - n_DOS_wid
       if(n_min<1) n_min = 1
       n_max = n_band + n_DOS_wid
       if(n_max>n_DOS) n_max = n_DOS
       do i = n_min, n_max
          Ebin = real(i-1,double)*dE_DOS + E_DOS_min
          a = (Ebin-eval(iwf))/sigma_DOS
          tmp(i) = weight*pf_DOS*exp(-half*a*a)
          DOS(i) = DOS(i) + tmp(i)
       end do

       ! Having found DOS, we now project onto atoms
       if(flag_write_projected_DOS) then
          if (atomf == sf) then
             acc = 0
             do atom=1,bundle%n_prim
                atom_spec = bundle%species(atom)
                fac = zero
                if (.not.flag_pDOS_angmom) then
                   do isf1 = 1,nsf_species(atom_spec)
                      fac = fac + real(evec(iwf,acc+isf1)*conjg(evec(iwf,acc+isf1)),double)
                   end do
                else if(flag_pDOS_lm) then
                   fac_angmom(:,:) = zero
                   isf1 = 0
                   do l1 = 0, pao(atom_spec)%greatest_angmom
                      do nacz1 = 1, pao(atom_spec)%angmom(l1)%n_zeta_in_angmom
                         if((pao(atom_spec)%angmom(l1)%semicore(nacz1)==0) .OR. &
                              (flag_pDOS_include_semicore)) then
                            do m1 = -l1,l1
                               isf1 = isf1 + 1
                               fac = fac + real(evec(iwf,acc+isf1)*conjg(evec(iwf,acc+isf1)),double)
                               ! l, m so shift m1 by l1+1 so it runs from 1 to 2*l1+1
                               fac_angmom(l1+1,m1+l1+1) = fac_angmom(l1+1,m1+l1+1) + &
                                    real(evec(iwf,acc+isf1)*conjg(evec(iwf,acc+isf1)),double)
                            enddo
                         end if
                      enddo
                   enddo
                else
                   fac_angmom(:,:) = zero
                   isf1 = 0
                   do l1 = 0, pao(atom_spec)%greatest_angmom
                      do nacz1 = 1, pao(atom_spec)%angmom(l1)%n_zeta_in_angmom
                         if((pao(atom_spec)%angmom(l1)%semicore(nacz1)==0) .OR. &
                              (flag_pDOS_include_semicore)) then
                            do m1 = -l1,l1
                               isf1 = isf1 + 1
                               fac = fac + real(evec(iwf,acc+isf1)*conjg(evec(iwf,acc+isf1)),double)
                               fac_angmom(l1+1,1) = fac_angmom(l1+1,1) &
                                    + real(evec(iwf,acc+isf1)*conjg(evec(iwf,acc+isf1)),double)
                            enddo
                         end if
                      enddo
                   enddo
                   if (isf1.ne.nsf_species(atom_spec)) call cq_abort("Error in NSF in the PDOS calculation.")
                endif
                if (flag_normalise_pDOS) then
                   fac = fac / weight_pDOS(iwf)
                   if (flag_pDOS_angmom) fac_angmom(:,:) = fac_angmom(:,:) / weight_pDOS(iwf)
                endif
                do i=n_min,n_max
                   projDOS(i,atom) = projDOS(i,atom) + tmp(i)*fac
                end do
                if (flag_pDOS_angmom) then
                   if(flag_pDOS_lm) then
                      do l1 = 0, pao(atom_spec)%greatest_angmom
                         do m1=-l1,l1
                            do i=n_min,n_max
                               projDOS_angmom(i,atom,l1+1,m1+l1+1) = projDOS_angmom(i,atom,l1+1,m1+l1+1) + &
                                    tmp(i)*fac_angmom(l1+1,m1+l1+1)
                            end do
                         end do
                      end do
                   else
                      do l1 = 0, pao(atom_spec)%greatest_angmom
                         do i=n_min,n_max
                            projDOS_angmom(i,atom,l1+1,1) = projDOS_angmom(i,atom,l1+1,1) + tmp(i)*fac_angmom(l1+1,1)
                         end do
                      end do
                   end if
                endif
                acc = acc + nsf_species(atom_spec)
             end do ! atom
          else
             spin_SF = 1
             if (flag_SpinDependentSF) spin_SF = spinSF
             acc = 0 
             iprim = 0
             do part = 1,bundle%groups_on_node ! Loop over primary set partitions
                if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
                   do memb = 1,bundle%nm_nodgroup(part) ! Loop over primary atoms
                      atom_num = bundle%nm_nodbeg(part)+memb-1
                      iprim=iprim+1
                      nsf1 = nsf_species(bundle%species(atom_num)) ! = mat(part,SFcoeff_range)%ndimi(memb)
                      do neigh = 1, mat(part,SFcoeff_range)%n_nab(memb) ! Loop over neighbours of atom
                         fac = zero
                         if (flag_pDOS_angmom) fac_angmom(:,:) = zero
                         ist = mat(part,SFcoeff_range)%i_acc(memb)+neigh-1
                         gcspart = BCS_parts%icover_ibeg(mat(part,SFcoeff_range)%i_part(ist))+ &
                                   mat(part,SFcoeff_range)%i_seq(ist)-1
                         neigh_global_part = BCS_parts%lab_cell(mat(part,SFcoeff_range)%i_part(ist))
                         neigh_global_num  = id_glob(parts%icell_beg(neigh_global_part)+ &
                                             mat(part,SFcoeff_range)%i_seq(ist)-1)
                         neigh_species = species_glob(neigh_global_num)
                         j_in_halo = halo(SFcoeff_range)%i_halo(gcspart)
                         natomf2 =  natomf_species(neigh_species)
                         ! Now loop over support functions and atomf (basically PAOs)
                         do isf1 = 1, nsf1
                            fac1 = real(evec(iwf,acc+isf1)*conjg(evec(iwf,acc+isf1)),double)
                            fac2 = zero
                            if (.not.flag_pDOS_angmom) then
                               do iatomf2 = 1, natomf2
                                  wheremat = matrix_pos(matSFcoeff(spin_SF),iprim,j_in_halo,isf1,iatomf2)
                                  val = mat_p(matSFcoeff(spin_SF))%matrix(wheremat)
                                  fac2 = fac2 + val*val
                               enddo
                               fac = fac + fac1 * fac2
                            else if(flag_pDOS_lm) then ! m and l resolved
                               fac2_angmom(:,:) = zero
                               iatomf2 = 0
                               do l2 = 0, pao(neigh_species)%greatest_angmom
                                  do nacz2 = 1, pao(neigh_species)%angmom(l2)%n_zeta_in_angmom
                                     if((pao(neigh_species)%angmom(l2)%semicore(nacz2)==0) .OR. &
                                          (flag_pDOS_include_semicore)) then
                                        do m2 = -l2,l2
                                           iatomf2 = iatomf2 + 1
                                           wheremat = matrix_pos(matSFcoeff(spin_SF),iprim,j_in_halo,isf1,iatomf2)
                                           val = mat_p(matSFcoeff(spin_SF))%matrix(wheremat)
                                           fac2 = fac2 + val*val
                                           fac2_angmom(l2+1,m2+l2+1) = fac2_angmom(l2+1,m2+l2+1) + val*val
                                        enddo ! m2
                                     end if
                                  enddo ! nacz2
                               enddo ! l2
                               fac = fac + fac1 * fac2
                               fac_angmom(:,:) = fac_angmom(:,:) + fac1 * fac2_angmom(:,:)
                            else ! NOT m resolved
                               fac2_angmom(:,:) = zero
                               iatomf2 = 0
                               do l2 = 0, pao(neigh_species)%greatest_angmom
                                  do nacz2 = 1, pao(neigh_species)%angmom(l2)%n_zeta_in_angmom
                                     if((pao(neigh_species)%angmom(l2)%semicore(nacz2)==0) .OR. &
                                          (flag_pDOS_include_semicore)) then
                                        do m2 = -l2,l2
                                           iatomf2 = iatomf2 + 1
                                           wheremat = matrix_pos(matSFcoeff(spin_SF),iprim,j_in_halo,isf1,iatomf2)
                                           val = mat_p(matSFcoeff(spin_SF))%matrix(wheremat)
                                           fac2 = fac2 + val*val
                                           fac2_angmom(l2+1,1) = fac2_angmom(l2+1,1) + val*val
                                        enddo ! m2
                                     end if
                                  enddo ! nacz2
                               enddo ! l2
                               fac = fac + fac1 * fac2
                               fac_angmom(:,1) = fac_angmom(:,1) + fac1 * fac2_angmom(:,1)
                               if (iatomf2.ne.natomf2) call cq_abort("Error in NATOMF in the PDOS calculation.")
                            endif ! flag_pDOS_angmom
                         enddo ! isf1
                         if (flag_normalise_pDOS) then
                            fac = fac / weight_pDOS(iwf)
                            if (flag_pDOS_angmom) fac_angmom(:,:) = fac_angmom(:,:) / weight_pDOS(iwf)
                         endif
                         ! project on the neighbour atom 
                         do i=n_min,n_max
                            projDOS(i,neigh_global_num) = projDOS(i,neigh_global_num) + tmp(i)*fac
                         end do
                         if (flag_pDOS_angmom) then
                            if(flag_pDOS_lm) then ! (l,m) resolved
                               do l2 = 0, pao(neigh_species)%greatest_angmom
                                  do m2 = -l2,l2
                                     do i=n_min,n_max
                                        projDOS_angmom(i,neigh_global_num,l2+1,m2+l2+1) = &
                                             projDOS_angmom(i,neigh_global_num,l2+1,m2+l2+1) + tmp(i)*fac_angmom(l2+1,l2+m2+1)
                                     end do ! i
                                  end do ! m2
                               end do ! l2
                            else ! l resolved only
                               do l2 = 0, pao(neigh_species)%greatest_angmom
                                  do i=n_min,n_max
                                     projDOS_angmom(i,neigh_global_num,l2+1,1) = projDOS_angmom(i,neigh_global_num,l2+1,1) &
                                          + tmp(i)*fac_angmom(l2+1,1)
                                  end do ! i
                               end do ! l2
                            end if
                         endif ! flag_pDOS_angmom
                      end do ! neigh
                      acc = acc + nsf1
                   end do ! memb
                end if ! nm_nodgroup
             end do ! part
          endif ! atomf
       end if ! flag_write_projected_DOS
    end do ! iwf
  end subroutine accumulate_DOS
  !!***
 
  !!****f*  DiagModule/weight_pDOS
  !!
  !!  NAME
  !!   weight_pDOS
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Normalise projected DOS
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   A. Nakata
  !!  CREATION DATE
  !!   07/09/2017
  !!  MODIFICATION HISTORY
  !!   21/11/2019
  !!    Renamed w_pDOS to weight_pDOS 
  !!  SOURCE
  !!
  subroutine get_weight_pDOS(evec,weight_pDOS)

    use datatypes
    use numbers, ONLY: zero
    use ScalapackFormat, only: matrix_size
    use species_module,  only: nsf_species
    use primary_module,  only: bundle

    implicit none

    ! Passed variables
    complex(double_cplx), dimension(:,:), intent(in) :: evec
    real(double), dimension(:) :: weight_pDOS

    ! Local variables
    integer :: iwf, acc, atom, nsf
    real(double) :: fac

    ! ---------------
    ! Calculate weight for each band to normalise pDOS
    ! ---------------
    weight_pDOS = zero
    do iwf=1,matrix_size ! Effectively all bands
       acc = 0
       fac = zero
       do atom=1,bundle%n_prim
          do nsf = 1,nsf_species(bundle%species(atom))
             fac = fac + real(evec(iwf,acc+nsf)*conjg(evec(iwf,acc+nsf)),double)
          end do
          acc = acc + nsf_species(bundle%species(atom))
       end do
       weight_pDOS(iwf) = fac
    end do
  end subroutine get_weight_pDOS
  !!***
 
  !!****f*  DiagModule/distrib_and_diag
  !!
  !!  NAME 
  !!   distrib_and_diag
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   For k-point i, distributes H and S and diagonalises
  !!
  !!      
  !!  INPUTS
  !!   spin: spin component (1 or 2)
  !!   i: k-point index
  !!   mode: 'N' (eigenvalues) or 'V' (eigenvectors)
  !!   flag_store_w: do we accumulate w ?
  !!  USES
  !! 
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2017/06/21
  !!  MODIFICATION HISTORY
  !!   2017/10/18 11:11 dave
  !!    Changed input variable i to index_kpoint
  !!   2017/10/18 14:50 dave
  !!    Added optional kpassed argument to allow calling of routine for arbitrary
  !!    (passed) kpoint.  In this case, set index_kpoint to one.
  !!   2017/10/19 09:05 dave
  !!    Added abort if we pass k-point with multiple process groups (not possible at present)
  !!   2018/09/05 14:25 dave
  !!    Updating the behaviour when info/=0
  !!   2018/11/13 17:30 nakata
  !!    Changed matS to be spin_SF dependent
  !!  SOURCE
  !!
  subroutine distrib_and_diag(spin,index_kpoint,mode,flag_store_w,kpassed)

    use datatypes
    use numbers
    use global_module,   only: iprint_DM, flag_SpinDependentSF
    use mult_module,     only: matH, matS
    use ScalapackFormat, only: matrix_size, proc_rows, proc_cols,     &
         nkpoints_max, pgid, N_kpoints_in_pg, pg_kpoints, N_procs_in_pg, proc_groups

    implicit none

    ! Passed
    integer :: spin ! Spin component
    integer :: index_kpoint ! K-point index in the process group list
    character(len=1) :: mode
    logical :: flag_store_w
    real(double), dimension(3), OPTIONAL :: kpassed

    ! Local
    real(double) :: vl, vu, orfac, scale
    integer :: il, iu, m, mz, info, spin_SF

    spin_SF = 1
    if (flag_SpinDependentSF) spin_SF = spin

    scale = one / real(N_procs_in_pg(pgid), double)
    vl = zero
    vu = zero
    orfac = -one
    il = 0
    iu = 0
    if (iprint_DM > 3 .and. (inode == ionode)) &
         write (io_lun, *) myid, ' Calling DistributeCQ_to_SC for H'
    ! Form the Hamiltonian and overlap for this k-point and send them to appropriate processors
    if(PRESENT(kpassed)) then
       if(proc_groups>1) call cq_abort("Coding error: can't have more than one PG and pass k-point to distrib_and_diag")
       call DistributeCQ_to_SC(DistribH, matH(spin), index_kpoint, SCHmat(:,:,spin),kpassed)
       call DistributeCQ_to_SC(DistribS, matS(spin_SF), index_kpoint, SCSmat(:,:,spin),kpassed)
    else
       call DistributeCQ_to_SC(DistribH, matH(spin), index_kpoint, SCHmat(:,:,spin))
       call DistributeCQ_to_SC(DistribS, matS(spin_SF), index_kpoint, SCSmat(:,:,spin))
    end if
    ! Now, if this processor is involved, do the diagonalisation
    if (iprint_DM > 3 .and. inode == ionode) &
         write (io_lun, *) myid, 'Proc row, cols, me: ', &
         proc_rows, proc_cols, me, index_kpoint, nkpoints_max
    if (index_kpoint <= N_kpoints_in_pg(pgid)) then
       ! Call the diagonalisation routine for generalised problem
       ! H.psi = E.S.psi
       call pzhegvx(1, mode, 'A', 'U', matrix_size, SCHmat(:,:,spin), &
            1, 1, desca, SCSmat(:,:,spin), 1, 1, descb,      &
            vl, vu, il, iu, abstol, m, mz, local_w(:,spin),  &
            orfac, z(:,:,spin), 1, 1, descz, work, lwork,    &
            rwork, lrwork, iwork, liwork, ifail, iclustr,    &
            gap, info)
       if (info /= 0) then
          if(info==2.OR.info==4) then ! These are safe to continue
             if(inode==ionode) then
                write(io_lun,fmt='(2x,"************************************")')
                write(io_lun,fmt='(2x,"** ScaLAPACK pzhegvx evec warning **")')
                write(io_lun,fmt='(2x,"** INFO=",i2," but continuing  **")') info
                write(io_lun,fmt='(2x,"************************************")')
             end if
          else
             call cq_abort ("FindEvals: pzhegvx failed for mode "//mode//" with INFO=", info)
          end if
       end if
       ! Copy local_w into appropriate place in w
       if(flag_store_w) w(1:matrix_size, pg_kpoints(pgid, index_kpoint), spin) = &
            scale * local_w(1:matrix_size, spin)
    end if ! End if (i<=N_kpoints_in_pg(pgid))
  end subroutine distrib_and_diag
  !!***

  !!****f*  DiagModule/init_blacs_pg
  !!
  !!  NAME 
  !!   init_blacs_pg
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   New subroutine which is run once per calculation to initialise BLACS and set up
  !!   progress groups for Conquest
  !!      
  !!  INPUTS
  !!
  !!  USES
  !! 
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2017/06/22
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine init_blacs_pg

    use global_module,   only: numprocs
    use GenComms,        only: my_barrier, cq_abort, myid
    use memory_module,   only: type_dbl, type_int, type_cplx,         &
         reg_alloc_mem, reg_dealloc_mem
    use ScalapackFormat, only: allocate_arrays, pg_initialise, ref_to_SC_blocks, make_maps, &
         block_size_r, block_size_c, proc_rows, &
         proc_cols, matrix_size, pgid, procid, proc_start
    
    implicit none

    ! Passed

    ! Local
    integer, dimension(:,:), allocatable :: imap
    integer :: nump, merow, mecol, numrows, numcols, info, stat
    
    call allocate_arrays(nkp)
    call pg_initialise(nkp) ! defined the process group parameters
    call ref_to_SC_blocks
    call make_maps
    ! Start up BLACS
    ! First check if there are enough process nodes
    call blacs_pinfo(me, nump)  ! get the total number of nodes
    ! avaliable for BLACS
    if (nump < numprocs) call cq_abort('initDiag: There are not enough nodes for BLACS', nump, numprocs)
    if (me /= myid) call cq_abort('initDiag: me and myid is not the same', me, myid)

    ! allocate imap for defining BLACS grid on the group
    allocate(imap(proc_rows, proc_cols), STAT=stat)
    if (stat /= 0) call cq_abort('initDiag: Failed to allocate imap', stat)
    call reg_alloc_mem(area_DM, proc_rows * proc_cols, type_int)

    ! assign the process grid map from ScalapackFormat procid, note
    ! that each context is local to the node, which associates to the
    ! map of the group the node belongs to
    imap(1:proc_rows, 1:proc_cols) = procid(pgid, 1:proc_rows, 1:proc_cols) - 1

    ! get the default system context
    call blacs_get(0, 0, context)
    ! replace the default context with the context defined for imap
    call blacs_gridmap(context, imap, proc_rows, proc_rows, proc_cols)

    ! imap is no longer required
    deallocate(imap, STAT=stat)
    if (stat /= 0) call cq_abort('initDiag: Failed to deallocate imap', stat)
    call reg_dealloc_mem(area_DM, proc_rows * proc_cols, type_int)

    ! check if we get the correct map
    call blacs_gridinfo(context, numrows, numcols, merow, mecol)
    if (iprint_DM > 3 .AND. myid == 0) then
       write (io_lun, fmt="(10x, 'process_grid info: ', i5, i5)") numrows, numcols
       write (io_lun, 1) myid, me, merow, mecol
    end if
    ! Sizes of local "chunk", used to initialise submatrix info for ScaLAPACK
    row_size = proc_start(myid+1)%rows * block_size_r
    col_size = proc_start(myid+1)%cols * block_size_c
    if (iprint_DM > 3 .AND. myid == 0) write (io_lun, 12) myid, row_size, col_size

    ! Register the description of the distribution of H
    call descinit(desca, matrix_size, matrix_size, block_size_r, block_size_c, 0, 0, context, row_size, info)
    if (info /= 0) call cq_abort("initDiag: descinit(a) failed !", info)
    ! Register the description of the distribution of S
    call descinit(descb, matrix_size, matrix_size, block_size_r, block_size_c, 0, 0, context, row_size, info)
    if (info /= 0) call cq_abort("initDiag: descinit(a) failed !", info)
    ! And register eigenvector distribution
    call descinit(descz, matrix_size, matrix_size, block_size_r, block_size_c, 0, 0, context, row_size, info)
    ! Find scratch space requirements for ScaLAPACk
    if (info /= 0) call cq_abort("initDiag: descinit(z) failed !", info)
1   format(10x, 'Proc: ', i5, ' BLACS proc, row, col: ', 3i5)
12  format(10x, 'Proc: ', i5, ' row, col size: ', 2i5)
    
  end subroutine init_blacs_pg
  !!***

  !!****f*  DiagModule/init_scalapack_format
  !!
  !!  NAME 
  !!   init_scalapack_format
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   New subroutine which is run once per atomic configuration to define the mapping
  !!   between Conquest and ScaLAPACK
  !!  INPUTS
  !!
  !!  USES
  !! 
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2017/06/22
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine init_scalapack_format

    use matrix_data,     only: Hrange, Srange
    use ScalapackFormat, only: find_SC_row_atoms, find_ref_row_atoms, find_SC_col_atoms, CC_to_SC
    use maxima_module, ONLY: maxnsf, maxpartscell, maxatomspart
         
    implicit none

    integer :: stat
    
    ! ScalapackFormat - do various translation tasks between Conquest
    ! compressed row storage and distributed Scalapack form. The
    ! arrays created will be used to form the Hamiltonian and K
    ! matrices.
    allocate(CC_to_SC(maxpartscell,maxatomspart,maxnsf), STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc CC2SC",stat)    
    call find_SC_row_atoms
    call find_ref_row_atoms
    call find_SC_col_atoms
    ! Now that's done, we can prepare to distribute things
    ! For Hamiltonian
    call PrepareRecv(DistribH)
    call PrepareSend(Hrange, DistribH)
    ! For Overlap
    call PrepareRecv(DistribS)
    call PrepareSend(Srange, DistribS)

  end subroutine init_scalapack_format
  !!***

  !!****f*  DiagModule/end_scalapack_format
  !!
  !!  NAME 
  !!   end_scalapack_format
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   New subroutine which is run once per atomic configuration to deallocate variables
  !!   from init_scalapack_format
  !!  INPUTS
  !!
  !!  USES
  !! 
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2017/06/22
  !!  MODIFICATION HISTORY
  !!  SOURCE
  !!
  subroutine end_scalapack_format

    use memory_module,   only: type_dbl, type_int, type_cplx, &
         reg_dealloc_mem
    use global_module,   only: numprocs
    use ScalapackFormat, only: CC_to_SC

    implicit none

    integer :: i, j, k, stat

    deallocate(CC_to_SC, STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc CC2SC",stat)    
    ! For S
    do i = 1, size(DistribS%images, 3)
       do j = 1, size(DistribS%images, 2)
          do k = 1, size(DistribS%images, 1)
             if (DistribS%images(k,j,i)%n_elements > 0) then
                deallocate(DistribS%images(k,j,i)%where, STAT=stat)
                if (stat /= 0) &
                     call cq_abort('endDiag: failed to deallocate (2)', stat)
             end if
          end do
       end do
    end do
    deallocate (DistribS%images, STAT=stat)
    if (stat /= 0) call cq_abort('endDiag: failed to deallocate (2)', stat)
    deallocate (DistribS%num_rows, DistribS%start_row, DistribS%send_rows, &
         DistribS%firstrow, STAT=stat)
    if (stat /= 0) call cq_abort('endDiag: failed to deallocate (2a)', stat)
    call reg_dealloc_mem(area_DM, 4 * numprocs, type_int)
    ! For H
    do i = 1, size(DistribH%images, 3)
       do j = 1, size(DistribH%images, 2)
          do k = 1, size(DistribH%images, 1)
             if (DistribH%images(k,j,i)%n_elements > 0) then
                deallocate (DistribH%images(k,j,i)%where, STAT=stat)
                if (stat /= 0) &
                     call cq_abort ('endDiag: failed to deallocate (2)', stat)
             end if
          end do
       end do
    end do
    deallocate (DistribH%images, STAT=stat)
    if (stat /= 0) call cq_abort('endDiag: failed to deallocate (2)', stat)
    deallocate (DistribH%num_rows, DistribH%start_row, DistribH%send_rows, &
         DistribH%firstrow, STAT=stat)
    if (stat /= 0) call cq_abort('endDiag: failed to deallocate (2a)', stat)
    call reg_dealloc_mem(area_DM, 4 * numprocs, type_int)

  end subroutine end_scalapack_format
  !!***
  
end module DiagModule
