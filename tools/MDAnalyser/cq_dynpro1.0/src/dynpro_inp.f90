MODULE dynpro_inp
  
  USE kinds, ONLY : DP
  
  IMPLICIT NONE
  
  SAVE
  
  !==---- Dynpro MPI ------------------------------------------------------------------------==
  
  INTEGER                       :: nproc, mpime, world, root, ios, ndr

  
  !==---- Dynpro XML ------------------------------------------------------------------------==
  !
  !  PRINCIPAL VARIABLE DEFINITION (from data-file.xml)
  !  
  !    nat_       == total number of atoms within the unit cell
  !    nsp_       == total number of species
  !    iteration_ == total number of frames
  !    time_      == total duration of the aiMD simulation
  !    atm_       == atom type
  !    label_     == atom label (name)
  !    ityp_      == atom index
  !    stau_      == ion postions in fractionel coordinates
  !    svel_      == ion velocities 
  !    force_     == ion forces
  !    ht_        == cell parameters       
  !==----------------------------------------------------------------------------------------==

  
  !REAL(DP)                      :: time_xml    
  !INTEGER                       :: iteration_xml, xmlunit

  INTEGER                       :: n_xml 
  INTEGER                       :: nat_xml, nsp_xml
  INTEGER                       :: nr1, nr2, nr3, tmp
  LOGICAL                       :: found

  CHARACTER(len=3),ALLOCATABLE  :: label_xml(:), atm_xml(:)
  INTEGER,         ALLOCATABLE  :: ityp_xml(:)
  !REAL(DP),        ALLOCATABLE  :: stau_xml(:,:), svel_xml(:,:), xyz_xml(:,:)
  !REAL(DP),        ALLOCATABLE  :: force_xml(:,:)
  !REAL(DP)                      :: ht_xml(3, 3)

  INTEGER :: stdin = 5, stderr = 6, ierr = 0
  
  !==---- Conquest input --------------------------------------------------------------------==
  INTEGER                        :: nat_cqinp, nsp_cqinp
  INTEGER,          ALLOCATABLE  :: ityp_cqcoord(:)
  CHARACTER(len=3), ALLOCATABLE  :: atm_cqinp(:)
  CHARACTER(len=3), ALLOCATABLE  :: label_cqinp(:)
  REAL(DP)                       :: cel_cqinp(3, 3)
  REAL(DP)                       :: evp_cqinp(10)

  LOGICAL                        :: lcq_coordforce
  LOGICAL                        :: lcq_mdframes
  
  !==---- Dynpro input ----------------------------------------------------------------------==
  LOGICAL              :: lau_unit  = .false. ! atomic units
  LOGICAL              :: lgrace    = .false. ! for Xmgrace format
  LOGICAL              :: lcenter   = .false. ! Centering the cell around an atom
  LOGICAL              :: lrdf      = .false. ! Radial distribution function
  LOGICAL              :: lsample   = .false. ! Sampling
  LOGICAL              :: lsphere   = .false. ! Solvation sphere reconstruction
  LOGICAL              :: lsub      = .false. ! Atom substitution
  LOGICAL              :: lvacf     = .false. ! Atomic velocity auto-correlation function
  LOGICAL              :: lnmr      = .false.
  LOGICAL              :: ldebug    = .false.

  REAL(DP)             :: time_step = 5.0d0   ! time step used in the FPMD in ua
  INTEGER              :: iprint    = 10      ! same as for CP parameters (conf. INPUT_CP.txt)
  INTEGER              :: nat, nsp, iteration
  LOGICAL              :: lspe      = .false.
  LOGICAL              :: lpsd      = .false.


  CHARACTER(len=3),ALLOCATABLE  :: label(:), atm(:), tmp_atm(:)
  INTEGER,         ALLOCATABLE  :: ityp(:),  na(:)
  INTEGER,         ALLOCATABLE  :: atomic_number(:)
  REAL(DP),        ALLOCATABLE  :: xyz(:,:), abc(:,:)
  REAL(DP),        ALLOCATABLE  :: vel(:,:), for(:,:)
  REAL(DP),        ALLOCATABLE  :: cel(:,:), evp(:,:)

  INTEGER,         ALLOCATABLE  :: ityp_read(:)
  INTEGER,         ALLOCATABLE  :: atm_read(:)
 
  
  !==---- Dynpro read/write --==
  CHARACTER(len=256)            :: fileout   = 'out'

  CHARACTER(len=256)            :: filepp, filemol, filerdf, filener, filedbg
  CHARACTER(len=256)            :: filecel, filepos, filefor, filevel, fileevp, filesort
  CHARACTER(len=256)            :: filesph, filevac, filepsd, filesph_tmp
  CHARACTER(len=256)            :: fileclas, filesub_in, filesub_out
  CHARACTER(len=256)            :: filedip, filespe, filesph_pc, filedip_tcf
  CHARACTER(len=256)            :: fileeig, filevec, filetens, fileran
  CHARACTER(len=256)            :: fileR, fileRiso, fileRant, fileRsym
  CHARACTER(len=256)            :: fileRsph

  CHARACTER(len=256)            :: prefix_cqdata
  CHARACTER(len=256)            :: filecqd

  
  CHARACTER(len=256)            :: filetoto

  !==---- CENTER cell --------==
  INTEGER              :: ct_atom   = 1       ! first atom of the matrix of coordinates    

  CHARACTER(len=2),ALLOCATABLE  :: syst_ID_lab(:)
 
  !==---- RDF calculation ----==
  REAL(DP)             :: rmin      = 1.0d0   ! Bohr unit
  REAL(DP)             :: rmax      = 0.0d0   ! rmax default value is set to the cell dimension see below
  REAL(DP)             :: rint(100) = 0.0d0   ! integraion radius for RDF; must below <= rmax  

  INTEGER              :: nr        = 151     ! nb. of points for the grid used for the RDFs calculations
  INTEGER              :: rd_atom   = 1       ! 

  !==----- Sampling ----------==
  INTEGER              :: sp_start  = 0
  INTEGER              :: sp_end    = 0
  INTEGER              :: sp_n      = 0
  INTEGER              :: start, end, interval

  REAL(DP)             :: time_dt
  
  !------ Sphere -------------==
  INTEGER              :: sph_atn(100) = 0
  INTEGER              :: sph_atom  = 0
  INTEGER              :: sph_n     = 0
  INTEGER              :: sph_m     = 0
  INTEGER              :: sph_q     = 0
  INTEGER              :: sph_p     = 0
  REAL(DP)             :: sph_rcut  = 1.2d0
  INTEGER              :: sph_natoms=0

  LOGICAL              :: langle      = .false.
  LOGICAL              :: langle_spef = .false.
  LOGICAL              :: ldihed_spef = .false.  
  LOGICAL              :: ldist       = .false.
  LOGICAL              :: ldist_spef  = .false.
  INTEGER              :: dist_atom   = 2
  INTEGER              :: angle_atom1 = 1
  INTEGER              :: angle_atom2 = 2
  INTEGER              :: angle_atom3 = 3

  INTEGER              :: ang_grid  = 512
  REAL(DP)             :: ang_start = 180.0d0
  REAL(DP)             :: ang_stop  =   0.0d0
  REAL(DP)             :: ang_test  = 120.0d0
 

  INTEGER,         ALLOCATABLE  :: sph_spe(:)
  REAL(DP),        ALLOCATABLE  :: sph_xyz(:,:)
  REAL(DP),        ALLOCATABLE  :: sph_vel(:,:)
  CHARACTER(len=3),ALLOCATABLE  :: sph_label(:)
  CHARACTER(len=6)              :: solvent_name


  INTEGER                       :: perm_xyz_tmp(1000)
  INTEGER                       :: perm_n = 1
  INTEGER,         ALLOCATABLE  :: perm_xyz(:)
  
  !==----- Substitute --------==
  LOGICAL              :: lstruct = .true.
  INTEGER              :: sub_num_atom(100) = 0
  INTEGER              :: sub_nsp           = 1
  INTEGER              :: sub_nsp_n(100)    = 1
  INTEGER              :: nsp_old
  CHARACTER(len=2)     :: sub_lab_atom(100) = 'XX'
  
  INTEGER,         ALLOCATABLE  :: ityp_old(:)
  REAL(DP),        ALLOCATABLE  :: xyz_sub(:,:)
  CHARACTER(len=3),ALLOCATABLE  :: label_old(:)

  INTEGER,         ALLOCATABLE  :: syst_ID_spe(:,:)
  REAL(DP),        ALLOCATABLE  :: atom_ID_xyz(:,:)

  REAL(DP),        ALLOCATABLE  :: sort_xyz(:,:), sort_vel(:,:)
 

  !------ VACF ---------------==
  !LOGICAL              :: lpsd      = .false.
  !LOGICAL              :: lspe      = .false.
  INTEGER              :: unit_vac

  !------ PSD ----------------==

  INTEGER              :: m_cof     = 0
  INTEGER              :: win       = 1
  INTEGER              :: n_seg     = 1

  LOGICAL              :: loverlap  = .false.
  REAL(DP)             :: temp      = 300.0d0
  INTEGER              :: qcfactor  = 0
  INTEGER              :: unit_psd
 

  REAL(DP),        ALLOCATABLE  :: C_nsp(:,:), C_tot(:), tmp_vac(:), tmp_1(:)

  !------ Dipole -----------==

  LOGICAL              :: ldip      = .false.
    
  !------ NMR Relaxation -----------==
  LOGICAL                       :: lrandom = .false.
  LOGICAL                       :: lstat   = .false.
  INTEGER                       :: verbose_nmr = 2
  REAL(DP)                      :: spin
  REAL(DP)                      :: quad
  CHARACTER(len=256)            :: prefix_nmr
  CHARACTER(len=256)            :: cor_function = 'default'
  INTEGER                       :: cor_order    = 2

END MODULE dynpro_inp
